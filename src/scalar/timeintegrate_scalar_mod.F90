MODULE timeintegrate_scalar_mod
    USE core_mod
    USE fieldmapper_mod
    USE ib_mod, ONLY: parent, ftoc, ib
    USE scacore_mod
    USE flow_mod, ONLY: ilesmodel, gmol, rho
    USE bound_scalar_mod
    USE itinfo_scalar_mod
    USE gc_scastencils_mod
    USE realfield_mod, ONLY: get_grid1_new, get_grid3_new

    IMPLICIT NONE(type, external)
    PRIVATE

    PUBLIC :: timeintegrate_scalar, itinfo_scalar

CONTAINS
    SUBROUTINE timeintegrate_scalar(itstep, ittot, timeph, dt, irk, rkscheme)
        ! Subroutine arguments
        INTEGER(intk), INTENT(in) :: itstep
        INTEGER(intk), INTENT(in) :: ittot
        REAL(realk), INTENT(in) :: timeph
        REAL(realk), INTENT(in) :: dt
        INTEGER(intk), INTENT(in) :: irk
        TYPE(rk_2n_t), INTENT(in) :: rkscheme

        ! Local variables
        INTEGER(intk) :: ilevel, l, i
        REAL(realk) :: frhs, fu, dtrk, dtrki
        TYPE(field_t), POINTER :: t, told, dt_f, qtt, qtu, qtv, qtw

        IF (.NOT. solve_scalar) RETURN
        CALL start_timer(400)

        CALL start_timer(401)
        ! (Cray) Frequently reallocating memory on an USM system is very
        !        costly and should be avoided. Thus, we keep the scalar
        !        fluxes allocated throughout the full simulation, which
        !        increases the total memory footprint
        CALL get_field(qtt, "QTT")
        CALL get_field(qtu, "QTU")
        CALL get_field(qtv, "QTV")
        CALL get_field(qtw, "QTW")
        CALL stop_timer(401)

        CALL start_timer(402)
        DO l = 1, nsca
            CALL start_timer(429)
            ! Fetch fields
            CALL get_field(t, scalar(l)%name)
            CALL get_field(dt_f, "D"//TRIM(scalar(l)%name))
            CALL get_field(told, TRIM(scalar(l)%name)//"_OLD")

            !!$omp target update to(mapper(maparr): t)
            CALL start_timer(430)
            ! Copy to "T_OLD" (not needed here but used for Boussinesq)
            !$omp target teams loop
            DO i = 1, SIZE(told%arr)
                told%arr(i) = t%arr(i)
            END DO
            !$omp end target teams loop
            CALL stop_timer(430)
            !!$omp target update from(mapper(maparr): told)

            !!$omp target update to(mapper(maparr): t)
            CALL start_timer(435)
            ! TSTSCA4 zeroize qtu, qtv, qtw before use internally
            CALL tstsca4(qtu, qtv, qtw, t, scalar(l))
            CALL stop_timer(435)
            !!$omp target update from(mapper(maparr): qtu, qtv, qtw)

            ! This operation apply boundary conditions to qtu, qtv, qtw ONLY!
            ! Does not modify t-field at all!
            DO ilevel = minlevel, maxlevel
                CALL start_timer(440)
                CALL parent(ilevel, qtu, qtv, qtw)
                CALL stop_timer(440)
                CALL start_timer(445)
                CALL bound_scaflux%bound(ilevel, qtu, qtv, qtw, t)
                CALL stop_timer(445)
            END DO
            !!$omp target update to(mapper(maparr): t, qtu, qtv, qtw)

            CALL start_timer(450)
            ! fluxbalance zeroize qtt before use internally
            CALL fluxbalance(qtt, qtu, qtv, qtw)
            CALL stop_timer(450)
            !!$omp target update from(mapper(maparr): qtt)

            CALL start_timer(455)
            ! Additional source terms
            CALL sourceterm(qtt, scalar(l))
            CALL stop_timer(455)

            CALL start_timer(460)
            ! Ghost cell "flux" boundary condition applied to qtt field
            IF (ib%type == "GHOSTCELL") THEN
                CALL set_scastencils("P", scalar(l), qtt=qtt)
            END IF
            CALL stop_timer(460)

            ! In IRK 1, FRHS is zero, therefore we do not need to zeroize
            ! the dt field before each step
            CALL rkscheme%get_coeffs(frhs, fu, dtrk, dtrki, irk)

            !!$omp target update to(mapper(maparr): qtt)
            CALL start_timer(465)
            ! dT_j = A_j*dT_(j-1) + QTT
            ! T_j = T_(j-1) + B_j*dT_j
            ! -------------------------
            ! t: read and write, dt_f: read and write, qtt: read.
            ! dt_f is not used anywhere else. Don't need to update
            CALL rkstep(t%arr, dt_f%arr, qtt%arr, frhs, dt*fu)
            ! qtt does not need to be copied back, it will be overridden next
            ! iteration anyways.
            CALL stop_timer(465)

            CALL start_timer(470)
            ! Mask blocked cells
            CALL maskbt(t)
            CALL stop_timer(470)
            !!$omp target update from(mapper(maparr): t)

            ! Ghost cell "value" boundary condition applied to t field
            IF (ib%type == "GHOSTCELL") THEN
                CALL start_timer(475)
                CALL connect(layers=2, s1=t, corners=.TRUE.)
                CALL stop_timer(475)
                CALL start_timer(480)
                CALL set_scastencils("P", scalar(l), t=t)
                CALL stop_timer(480)
            END IF

            DO ilevel = maxlevel, minlevel+1, -1
                CALL start_timer(485)
                CALL ftoc(ilevel, t, t, 'T')
                CALL stop_timer(485)
            END DO

            CALL start_timer(490)
            CALL connect(layers=2, s1=t, corners=.TRUE.)
            CALL stop_timer(490)

            ! TODO: Fill ghost layers of T (maybe only at last IRK?)
            CALL stop_timer(429)
        END DO
        CALL stop_timer(402)

        CALL stop_timer(400)
    END SUBROUTINE timeintegrate_scalar


    SUBROUTINE itinfo_scalar(itstep, ittot, timeph, dt, exploded)
        ! Subroutine arguments
        INTEGER(intk), INTENT(in) :: itstep
        INTEGER(intk), INTENT(in) :: ittot
        REAL(realk), INTENT(in) :: timeph
        REAL(realk), INTENT(in) :: dt
        INTEGER(intk), INTENT(inout) :: exploded

        ! Local variables
        INTEGER(intk) :: i, l, igrid
        INTEGER(intk) :: kk, jj, ii
        TYPE(field_t), POINTER :: ddx_f, ddy_f, ddz_f

        REAL(realk), POINTER, CONTIGUOUS :: t(:, :, :)
        REAL(realk), POINTER, CONTIGUOUS :: ddx(:), ddy(:), ddz(:)

        REAL(realk) :: tmean(nsca), tmeansqr(nsca)

        IF (.NOT. solve_scalar) RETURN
        CALL start_timer(400)
        CALL start_timer(420)

        ! Get fields
        CALL get_field(ddx_f, "DDX")
        CALL get_field(ddy_f, "DDY")
        CALL get_field(ddz_f, "DDZ")

        ! Compute CFL, divergence and kinetic energy
        DO i = 1, nmygrids
            igrid = mygrids(i)

            CALL get_mgdims(kk, jj, ii, igrid)

            CALL ddx_f%get_ptr(ddx, igrid)
            CALL ddy_f%get_ptr(ddy, igrid)
            CALL ddz_f%get_ptr(ddz, igrid)

            DO l = 1, nsca
                CALL get_fieldptr(t, scalar(l)%name, igrid)
                CALL comp_tmean(tmean(l), tmeansqr(l), kk, jj, ii, t, ddx, &
                    ddy, ddz)
            END DO

            CALL itinfo_scalar_sample(igrid, tmean, tmeansqr)
        END DO

        CALL itinfo_scalar_print(itstep, ittot, timeph, exploded)

        CALL stop_timer(420)
        CALL stop_timer(400)
    END SUBROUTINE itinfo_scalar


    SUBROUTINE tstsca4(qtu_f, qtv_f, qtw_f, t_f, sca)
        ! Subroutine arguments
        TYPE(field_t), INTENT(inout) :: qtu_f, qtv_f, qtw_f
        TYPE(field_t), INTENT(in) :: t_f
        TYPE(scalar_t), INTENT(in) :: sca

        ! Local variables
        INTEGER(intk) :: i, igrid
        INTEGER(intk) :: kk, jj, ii, ip3, ipx, ipy, ipz
        INTEGER(intk) :: nfro, nbac, nrgt, nlft, nbot, ntop
        TYPE(field_t), POINTER :: u_f, v_f, w_f, g_f
        TYPE(field_t), POINTER :: bt_f, ddx_f, ddy_f, ddz_f, rdx_f, rdy_f, rdz_f
        REAL, CONTIGUOUS, POINTER, DIMENSION(:, :, :) :: qtu, qtv, qtw, t, u, v, w, g, bt
        REAL, CONTIGUOUS, POINTER, DIMENSION(:) :: ddx, ddy, ddz, rdx, rdy, rdz

        CALL start_timer(410)

        CALL get_field(u_f, "U")
        CALL get_field(v_f, "V")
        CALL get_field(w_f, "W")
        CALL get_field(g_f, "G")

        CALL get_field(bt_f, "BT")
        CALL get_field(ddx_f, "DDX")
        CALL get_field(ddy_f, "DDY")
        CALL get_field(ddz_f, "DDZ")
        CALL get_field(rdx_f, "RDX")
        CALL get_field(rdy_f, "RDY")
        CALL get_field(rdz_f, "RDZ")

        !$omp target teams loop
        DO i = 1, SIZE(qtu_f%arr)
            qtu_f%arr(i) = 0.0_realk
        END DO
        !$omp end target teams loop

        !$omp target teams loop
        DO i = 1, SIZE(qtv_f%arr)
            qtv_f%arr(i) = 0.0_realk
        END DO
        !$omp end target teams loop

        !$omp target teams loop
        DO i = 1, SIZE(qtw_f%arr)
            qtw_f%arr(i) = 0.0_realk
        END DO
        !$omp end target teams loop

        !$omp target teams distribute private(igrid, kk, jj, ii, nfro, nbac, nrgt, nlft, nbot, ntop, &
        !$omp& qtu, qtv, qtw, t, u, v, w, g, bt, ddx, ddy, ddz, rdx, rdy, rdz)
        DO i = 1, nmygrids
            igrid = mygrids(i)

            CALL get_mgdims(kk, jj, ii, igrid)
            CALL get_mgbasb(nfro, nbac, nrgt, nlft, nbot, ntop, igrid)

            CALL get_grid3_new(qtu_f, qtu, igrid)
            CALL get_grid3_new(qtv_f, qtv, igrid)
            CALL get_grid3_new(qtw_f, qtw, igrid)

            CALL get_grid3_new(t_f, t, igrid)
            CALL get_grid3_new(u_f, u, igrid)
            CALL get_grid3_new(v_f, v, igrid)
            CALL get_grid3_new(w_f, w, igrid)
            CALL get_grid3_new(g_f, g, igrid)
            CALL get_grid3_new(bt_f, bt, igrid)

            CALL get_grid1_new(ddx_f, ddx, igrid)
            CALL get_grid1_new(ddy_f, ddy, igrid)
            CALL get_grid1_new(ddz_f, ddz, igrid)
            CALL get_grid1_new(rdx_f, rdx, igrid)
            CALL get_grid1_new(rdy_f, rdy, igrid)
            CALL get_grid1_new(rdz_f, rdz, igrid)

            CALL tstsca4_grid(kk, jj, ii, qtu, qtv, qtw, t, u, v, w, g, bt, &
                ddx, ddy, ddz, rdx, rdy, rdz, sca, nfro, nbac, nrgt, nlft, &
                nbot, ntop)
        END DO
        !$omp end target teams distribute

        CALL stop_timer(410)
    END SUBROUTINE tstsca4


    SUBROUTINE tstsca4_grid(kk, jj, ii, qtu, qtv, qtw, t, u, v, w, g, bt, &
            ddx, ddy, ddz, rdx, rdy, rdz, sca, nfro, nbac, nrgt, nlft, &
            nbot, ntop)
        !$omp declare target
        ! Subroutine arguments
        INTEGER(intk), INTENT(IN) :: kk, jj, ii
        REAL(realk), INTENT(OUT), DIMENSION(kk, jj, ii) :: qtu, qtv, qtw
        REAL(realk), INTENT(IN), DIMENSION(kk, jj, ii) :: t, u, v, w, g, bt
        REAL(realk), INTENT(IN) :: ddx(ii), ddy(jj), ddz(kk), &
            rdx(ii), rdy(jj), rdz(kk)
        TYPE(scalar_t), INTENT(in) :: sca
        INTEGER(intk), INTENT(IN) :: nfro, nbac, nrgt, nlft, nbot, ntop

        ! Local variables
        INTEGER(intk) :: i, j, k
        INTEGER(intk) :: nfu, nbu, nrv, nlv, nbw, ntw
        INTEGER(intk) :: iles
        REAL(realk) :: gsca, adv, diff, area
        ! REAL(realk) :: gscamol, gtgmolp, gtgmoln

        ! Usually, the fluxes across the grid boundaries are already set
        nfu = 0
        nbu = 0
        nrv = 0
        nlv = 0
        nbw = 0
        ntw = 0

        ! Only for CON boundaries, fluxes are computed for one more layer
        ! (this avoids a connect on qtu, qtv, qtw)
        IF (nfro == 7) nfu = 1
        IF (nbac == 7) nbu = 1
        IF (nrgt == 7) nrv = 1
        IF (nlft == 7) nlv = 1
        IF (nbot == 7) nbw = 1
        IF (ntop == 7) ntw = 1

        iles = 1
        IF (ilesmodel == 0) iles = 0

        ! X direction
        !$omp parallel do collapse(3) private(gsca, adv, area, diff)
        DO i = 3-nfu, ii-3+nbu
            DO j = 3, jj-2
                ! Scalar diffusivity LES/DNS computation
                DO k = 3, kk - 2
                    ! IF (iles == 1) THEN
                    !     gscamol = gmol/rho/sca%prmol
                    !     gtgmolp = (g(k, j, i) - gmol)/gmol
                    !     gtgmoln = (g(k, j, i+1) - gmol)/gmol

                    !     ! 1/Re * 1/Pr + 1/Re_t * 1/Pr_t:
                    !     gsca = gscamol &
                    !         + (g(k, j, i+1) + g(k, j, i) - 2.0*gmol) / rho &
                    !         / (sca%prt(gtgmoln) + sca%prt(gtgmolp))

                    !     ! Limit gsca here MAX(..., 0): no negative diffusion!
                    !     gsca = MAX(gscamol, gsca)
                    ! ELSE
                        gsca = gmol/rho/sca%prmol
                    ! END IF

                    ! Convective fluxes
                    ! It is assumed that the velocity field is already masked
                    ! with BU, BV, BW = no new masking necessary (!)
                    adv = (ddy(j)*ddz(k)) * u(k, j, i) &
                        * 0.5 * (t(k, j, i) + t(k, j, i+1))

                    ! Depending on the knowledge about the the cell and its
                    ! neighbours it is determined if faces are blocked (=0)
                    ! or open (=1)
                    area = bt(k, j, i)*bt(k, j, i+1)*(ddy(j)*ddz(k))
                    diff = -gsca*rdx(i)*(t(k, j, i+1) - t(k, j, i))*area

                    ! Final result
                    qtu(k, j, i) = adv + diff
                END DO
            END DO
        END DO
        !$omp end parallel do

        ! Y direction
        !$omp parallel do collapse(3) private(gsca, adv, area, diff)
        DO i = 3, ii-2
            DO j = 3-nrv, jj-3+nlv
                ! Scalar diffusivity LES/DNS computation
                DO k = 3, kk-2
                    ! IF (iles == 1) THEN
                    !     gscamol = gmol/rho/sca%prmol
                    !     gtgmolp = (g(k, j, i) - gmol)/gmol
                    !     gtgmoln = (g(k, j+1, i) - gmol)/gmol

                    !     ! 1/Re * 1/Pr + 1/Re_t * 1/Pr_t:
                    !     gsca = gscamol &
                    !         + (g(k, j+1, i) + g(k, j, i) - 2.0*gmol) / rho &
                    !         / (sca%prt(gtgmoln) + sca%prt(gtgmolp))

                    !     ! Limit gsca here MAX(..., 0): no negative diffusion!
                    !     gsca = MAX(gscamol, gsca)
                    ! ELSE
                        gsca = gmol/rho/sca%prmol
                    ! END IF

                    ! Convective fluxes
                    ! It is assumed that the velocity field is already masked
                    ! with BU, BV, BW = no new masking necessary (!)
                    adv = (ddx(i)*ddz(k)) * v(k, j, i) &
                        * 0.5 * (t(k, j, i) + t(k, j+1, i))

                    ! Depending on the knowledge about the the cell and its
                    ! neighbours it is determined if faces are blocked (=0)
                    ! or open (=1)
                    area = bt(k, j, i)*bt(k, j+1, i)*(ddx(i)*ddz(k))
                    diff = -gsca*rdy(j)*(t(k, j+1, i) - t(k, j, i))*area

                    ! Final result
                    qtv(k, j, i) = adv + diff
                END DO
            END DO
        END DO
        !$omp end parallel do

        ! Z direction
        !$omp parallel do collapse(3) private(gsca, adv, area, diff)
        DO i = 3, ii-2
            DO j = 3, jj-2
                ! Scalar diffusivity LES/DNS computation
                DO k = 3-nbw, kk-3+ntw
                    ! IF (iles == 1) THEN
                    !     gscamol = gmol/rho/sca%prmol
                    !     gtgmolp = (g(k, j, i) - gmol)/gmol
                    !     gtgmoln = (g(k+1, j, i) - gmol)/gmol

                    !     ! 1/Re * 1/Pr + 1/Re_t * 1/Pr_t:
                    !     gsca = gscamol &
                    !         + (g(k+1, j, i) + g(k, j, i) - 2.0*gmol) / rho &
                    !         / (sca%prt(gtgmoln) + sca%prt(gtgmolp))

                    !     ! Limit gsca here MAX(..., 0): no negative diffusion!
                    !     gsca = MAX(gscamol, gsca)
                    ! ELSE
                        gsca = gmol/rho/sca%prmol
                    ! END IF

                    ! Convective fluxes
                    ! It is assumed that the velocity field is already masked
                    ! with BU, BV, BW = no new masking necessary (!)
                    adv = (ddx(i)*ddy(j)) * w(k, j, i) &
                        * 0.5 * (t(k, j, i) + t(k+1, j, i))

                    ! Depending on the knowledge about the the cell and its
                    ! neighbours it is determined if faces are blocked (=0)
                    ! or open (=1)
                    area = bt(k, j, i)*bt(k+1, j, i)*(ddx(i)*ddy(j))
                    diff = -gsca*rdz(k)*(t(k+1, j, i) - t(k, j, i))*area

                    ! Final result
                    qtw(k, j, i) = adv + diff
                END DO
            END DO
        END DO
        !$omp end parallel do


        ! Special treatment at par boundaries
        ! Substraction of downwind and addition of upwind T-value
        ! to finally get an upwind scheme in case of flow towards coarse grid
        IF (nfro == 8) THEN
            i =  3
            !$omp parallel do collapse(2) private(adv)
            DO j = 3, jj-2
                DO k = 3, kk-2
                    adv = (ddy(j)*ddz(k)) * (u(k, j, i) - ABS(u(k, j, i))) &
                        * 0.5 * 0.5 * (-t(k, j, i) + t(k, j, i+1))
                    qtu(k, j, i) = qtu(k, j, i) + adv
                END DO
            END DO
            !$omp end parallel do
        END IF

        IF (nbac == 8) THEN
            i = ii-3
            !$omp parallel do collapse(2) private(adv)
            DO j = 3, jj-2
                DO k = 3, kk-2
                    adv = (ddy(j)*ddz(k)) * (u(k, j, i) + ABS(u(k, j, i))) &
                        * 0.5 * 0.5 * (t(k, j, i) - t(k, j, i+1))
                    qtu(k, j, i) = qtu(k, j, i) + adv
                END DO
            END DO
            !$omp end parallel do
        END IF

        IF (nrgt == 8) THEN
            j = 3
            !$omp parallel do collapse(2) private(adv)
            DO i = 3, ii-2
                DO k = 3, kk-2
                    adv = (ddx(i)*ddz(k)) * (v(k, j, i) - ABS(v(k, j, i))) &
                        * 0.5 * 0.5 * (-t(k, j, i) + t(k, j+1, i))
                    qtv(k, j, i) = qtv(k, j, i) + adv
                END DO
            END DO
            !$omp end parallel do
        END IF

        IF (nlft == 8) THEN
            j = jj-3
            !$omp parallel do collapse(2) private(adv)
            DO i = 3, ii-2
                DO k = 3, kk-2
                    adv = (ddx(i)*ddz(k)) * (v(k, j, i) + ABS(v(k, j, i))) &
                        * 0.5 * 0.5 * (t(k, j, i) - t(k, j+1, i))
                    qtv(k, j, i) = qtv(k, j, i) + adv
                END DO
            END DO
            !$omp end parallel do
        END IF

        IF (nbot == 8) THEN
            k = 3
            !$omp parallel do collapse(2) private(adv)
            DO i = 3, ii-2
                DO j = 3, jj-2
                    adv = (ddx(i)*ddy(j)) * (w(k, j, i) - ABS(w(k, j, i))) &
                        * 0.5 * 0.5 * (-t(k, j, i) + t(k+1, j, i))
                    qtw(k, j, i) = qtw(k, j, i) + adv
                END DO
            END DO
            !$omp end parallel do
        END IF

        IF (ntop == 8) THEN
            k = kk-3
            !$omp parallel do collapse(2) private(adv)
            DO i = 3, ii-2
                DO j = 3, jj-2
                    adv = (ddx(i)*ddy(j)) * (w(k, j, i) + ABS(w(k, j, i))) &
                        * 0.5 * 0.5 * (t(k, j, i) - t(k+1, j, i))
                    qtw(k, j, i) = qtw(k, j, i) + adv
                END DO
            END DO
            !$omp end parallel do
        END IF
    END SUBROUTINE tstsca4_grid


    SUBROUTINE fluxbalance(qtt_f, qtu_f, qtv_f, qtw_f)
        ! Subroutine arguments
        TYPE(field_t), INTENT(inout) :: qtt_f
        TYPE(field_t), INTENT(in) :: qtu_f, qtv_f, qtw_f

        ! Local variables
        INTEGER(intk) :: i, igrid
        INTEGER(intk) :: kk, jj, ii, ip3, ipx, ipy, ipz
        TYPE(field_t), POINTER :: rddx_f, rddy_f, rddz_f
        REAL, CONTIGUOUS, POINTER, DIMENSION(:, :, :) :: qtt, qtu, qtv, qtw
        REAL, CONTIGUOUS, POINTER, DIMENSION(:) :: rddx, rddy, rddz

        CALL start_timer(411)

        CALL get_field(rddx_f, "RDDX")
        CALL get_field(rddy_f, "RDDY")
        CALL get_field(rddz_f, "RDDZ")

        !$omp target teams loop
        DO i = 1, SIZE(qtt_f%arr)
            qtt_f%arr(i) = 0.0_realk
        END DO
        !$omp end target teams loop

        !$omp target teams distribute private(igrid, kk, jj, ii, qtt, qtu, qtv, qtw, rddx, rddy, rddz)
        DO i = 1, nmygrids
            igrid = mygrids(i)

            CALL get_mgdims(kk, jj, ii, igrid)
            CALL get_grid3_new(qtt_f, qtt, igrid)
            CALL get_grid3_new(qtu_f, qtu, igrid)
            CALL get_grid3_new(qtv_f, qtv, igrid)
            CALL get_grid3_new(qtw_f, qtw, igrid)
            CALL get_grid1_new(rddx_f, rddx, igrid)
            CALL get_grid1_new(rddy_f, rddy, igrid)
            CALL get_grid1_new(rddz_f, rddz, igrid)

            CALL fluxbalance_grid(kk, jj, ii, qtt, qtu, qtv, qtw, &
                rddx, rddy, rddz)
        END DO
        !$omp end target teams distribute
        CALL stop_timer(411)
    END SUBROUTINE fluxbalance


    SUBROUTINE fluxbalance_grid(kk, jj, ii, qtt, qtu, qtv, qtw, &
            rddx, rddy, rddz)
        !$omp declare target
        ! Subroutine arguments
        INTEGER(intk), INTENT(IN) :: kk, jj, ii
        REAL(realk), INTENT(OUT), DIMENSION(kk, jj, ii) :: qtt
        REAL(realk), INTENT(IN), DIMENSION(kk, jj, ii) :: qtu, qtv, qtw
        REAL(realk), INTENT(IN) :: rddx(ii), rddy(jj), rddz(kk)

        ! Local variables
        INTEGER(intk) :: i, j, k
        REAL(realk) :: netflux

        !$omp parallel do collapse(3) private(netflux)
        DO i = 3, ii-2
            DO j = 3, jj-2
                DO k = 3, kk-2
                    ! Computing netflux resulting from exchange with neighbors
                    netflux = qtu(k, j, i-1) - qtu(k, j, i) + qtv(k, j-1, i) &
                        - qtv(k, j, i) + qtw(k-1, j, i) - qtw(k, j, i)

                    qtt(k, j, i) = rddz(k)*rddy(j)*rddx(i)*netflux
                END DO
            END DO
        END DO
        !$omp end parallel do
    END SUBROUTINE fluxbalance_grid


    SUBROUTINE sourceterm(qtt_f, sca)
        ! Subroutine arguments
        TYPE(field_t), INTENT(inout) :: qtt_f
        TYPE(scalar_t), INTENT(in), TARGET :: sca

        ! Local variables
        INTEGER(intk) :: i, igrid, isource, nsource
        INTEGER(intk) :: kk, jj, ii
        TYPE(scalar_source_t), POINTER :: source
        TYPE(field_t), POINTER :: field_f
        REAL(realk), POINTER, CONTIGUOUS :: qtt(:, :, :), field(:, :, :)

        CALL start_timer(413)

        nsource = SIZE(sca%sources)

        DO isource = 1, nsource
            source => sca%sources(isource)

            IF (LEN_TRIM(source%field) == 0) THEN
                DO i = 1, nmygrids
                    igrid = mygrids(i)
                    CALL get_mgdims(kk, jj, ii, igrid)
                    CALL qtt_f%get_ptr(qtt, igrid)

                    CALL sourceterm_const(kk, jj, ii, qtt, source%value)
                END DO
            ELSE
                CALL get_field(field_f, source%field)
                IF (field_f%istag /= 0 .OR. field_f%jstag /= 0 .OR. &
                        field_f%kstag /= 0) THEN
                    CALL errr(__FILE__, __LINE__)
                END IF

                DO i = 1, nmygrids
                    igrid = mygrids(i)
                    CALL get_mgdims(kk, jj, ii, igrid)
                    CALL qtt_f%get_ptr(qtt, igrid)
                    CALL field_f%get_ptr(field, igrid)

                    CALL sourceterm_field(kk, jj, ii, qtt, field, source%value)
                END DO
            END IF
        END DO

        CALL stop_timer(413)
    END SUBROUTINE sourceterm


    PURE SUBROUTINE sourceterm_const(kk, jj, ii, qtt, sourceval)
        ! Subroutine arguments
        INTEGER(intk), INTENT(IN) :: kk, jj, ii
        REAL(realk), INTENT(INOUT), DIMENSION(kk, jj, ii) :: qtt
        REAL(realk), INTENT(IN) :: sourceval

        ! Local variables
        INTEGER(intk) :: i, j, k

        DO i = 3, ii-2
            DO j = 3, jj-2
                DO k = 3, kk-2
                    qtt(k, j, i) = qtt(k, j, i) + sourceval
                END DO
            END DO
        END DO
    END SUBROUTINE sourceterm_const


    PURE SUBROUTINE sourceterm_field(kk, jj, ii, qtt, field, sourceval)
        ! Subroutine arguments
        INTEGER(intk), INTENT(IN) :: kk, jj, ii
        REAL(realk), INTENT(INOUT) :: qtt(kk, jj, ii)
        REAL(realk), INTENT(IN) :: field(kk, jj, ii)
        REAL(realk), INTENT(IN) :: sourceval

        ! Local variables
        INTEGER(intk) :: i, j, k

        DO i = 3, ii-2
            DO j = 3, jj-2
                DO k = 3, kk-2
                    qtt(k, j, i) = qtt(k, j, i) + sourceval*field(k, j, i)
                END DO
            END DO
        END DO
    END SUBROUTINE sourceterm_field


    SUBROUTINE comp_tmean(tmean, tmeansqr, kk, jj, ii, t, ddx, ddy, ddz)
        ! Subroutine arguments
        REAL(realk), INTENT(out) :: tmean, tmeansqr
        INTEGER(intk), INTENT(in) :: kk, jj, ii
        REAL(realk), INTENT(in) :: t(kk, jj, ii)
        REAL(realk), INTENT(in) :: ddx(ii), ddy(jj), ddz(kk)

        ! Local variables
        REAL(realk) :: vsum, vol
        INTEGER(intk) :: k, j, i

        tmean = 0.0
        tmeansqr = 0.0
        vsum = 0.0

        DO i = 3, ii-2
            DO j = 3, jj-2
                DO k = 3, kk-2
                    vol = ddx(i)*ddy(j)*ddz(k)
                    tmean = tmean + vol*t(k, j, i)
                    tmeansqr = tmeansqr + vol*t(k, j, i)**2
                    vsum = vsum + vol
                END DO
            END DO
        END DO

        tmean = tmean/vsum
        tmeansqr = tmeansqr/vsum
    END SUBROUTINE comp_tmean
END MODULE timeintegrate_scalar_mod
