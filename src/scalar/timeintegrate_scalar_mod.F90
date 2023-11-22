MODULE timeintegrate_scalar_mod
    USE core_mod
    USE ib_mod, ONLY: parent, ftoc, ib
    USE scacore_mod
    USE flow_mod, ONLY: ilesmodel, gmol, rho
    USE bound_scalar_mod
    USE itinfo_scalar_mod
    USE gc_scastencils_mod

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
        INTEGER(intk) :: ilevel, l
        REAL(realk) :: frhs, fu, dtrk, dtrki
        TYPE(field_t), POINTER :: t, told, dt_f
        TYPE(field_t) :: qtt, qtu, qtv, qtw

        IF (.NOT. solve_scalar) RETURN
        CALL start_timer(400)
        CALL start_timer(401)

        ! Local temporary storage ("scrap")
        CALL qtt%init("QTT")
        CALL qtu%init("QTU", istag=1)
        CALL qtv%init("QTV", jstag=1)
        CALL qtw%init("QTW", kstag=1)

        CALL qtu%init_buffers()
        CALL qtv%init_buffers()
        CALL qtw%init_buffers()
        CALL stop_timer(401)

        CALL start_timer(402)
        DO l = 1, nsca
            ! Fetch fields
            CALL get_field(t, scalar(l)%name)
            CALL get_field(dt_f, "D"//TRIM(scalar(l)%name))
            CALL get_field(told, TRIM(scalar(l)%name)//"_OLD")

            ! Copy to "T_OLD"
            told%arr = t%arr

            ! TSTSCA4 zeroize qtu, qtv, qtw before use internally
            CALL tstsca4(qtu, qtv, qtw, t, scalar(l))

            ! This operation apply boundary conditions to qtu, qtv, qtw ONLY!
            ! Does not modify t-field at all!
            DO ilevel = minlevel, maxlevel
                CALL parent(ilevel, qtu, qtv, qtw)
                CALL bound_scaflux%bound(ilevel, qtu, qtv, qtw, t)
            END DO

            ! fluxbalance zeroize qtt before use internally
            CALL fluxbalance(qtt, qtu, qtv, qtw)

            ! Ghost cell "flux" boundary condition applied to qtt field
            IF (ib%type == "GHOSTCELL") THEN
                CALL set_scastencils("P", scalar(l), qtt=qtt)
            END IF

            ! In IRK 1, FRHS is zero, therefore we do not need to zeroize
            ! the dt field before each step
            CALL rkscheme%get_coeffs(frhs, fu, dtrk, dtrki, irk)

            ! dT_j = A_j*dT_(j-1) + QTT
            dt_f%arr = frhs*dt_f%arr + qtt%arr

            ! T_j = T_(j-1) + B_j*dT_j
            t%arr = t%arr + dt*fu*dt_f%arr

            ! Mask blocked cells
            CALL maskbt(t)

            ! Ghost cell "value" boundary condition applied to t field
            IF (ib%type == "GHOSTCELL") THEN
                CALL connect(layers=2, s1=t, corners=.TRUE.)
                CALL set_scastencils("P", scalar(l), t=t)
            END IF

            DO ilevel = maxlevel, minlevel+1, -1
                CALL ftoc(ilevel, t%arr, t%arr, 'T')
            END DO

            CALL connect(layers=2, s1=t, corners=.TRUE.)

            ! TODO: Fill ghost layers of T (maybe only at last IRK?)
        END DO
        CALL stop_timer(402)

        CALL qtt%finish()
        CALL qtu%finish()
        CALL qtv%finish()
        CALL qtw%finish()

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
        INTEGER(intk) :: kk, jj, ii
        INTEGER(intk) :: nfro, nbac, nrgt, nlft, nbot, ntop
        TYPE(field_t), POINTER :: u_f, v_f, w_f, g_f
        TYPE(field_t), POINTER :: bt_f, ddx_f, ddy_f, ddz_f, rdx_f, rdy_f, rdz_f
        REAL(realk), POINTER, CONTIGUOUS, DIMENSION(:, :, :) :: qtu, qtv, qtw, &
            t, u, v, w, g, bt
        REAL(realk), POINTER, CONTIGUOUS, DIMENSION(:) :: ddx, ddy, ddz, &
            rdx, rdy, rdz

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

        DO i = 1, nmygrids
            igrid = mygrids(i)

            CALL get_mgdims(kk, jj, ii, igrid)
            CALL get_mgbasb(nfro, nbac, nrgt, nlft, nbot, ntop, igrid)

            CALL qtu_f%get_ptr(qtu, igrid)
            CALL qtv_f%get_ptr(qtv, igrid)
            CALL qtw_f%get_ptr(qtw, igrid)

            CALL t_f%get_ptr(t, igrid)
            CALL u_f%get_ptr(u, igrid)
            CALL v_f%get_ptr(v, igrid)
            CALL w_f%get_ptr(w, igrid)
            CALL g_f%get_ptr(g, igrid)
            CALL bt_f%get_ptr(bt, igrid)

            CALL ddx_f%get_ptr(ddx, igrid)
            CALL ddy_f%get_ptr(ddy, igrid)
            CALL ddz_f%get_ptr(ddz, igrid)
            CALL rdx_f%get_ptr(rdx, igrid)
            CALL rdy_f%get_ptr(rdy, igrid)
            CALL rdz_f%get_ptr(rdz, igrid)

            CALL tstsca4_grid(kk, jj, ii, qtu, qtv, qtw, t, u, v, w, g, bt, &
                ddx, ddy, ddz, rdx, rdy, rdz, sca, nfro, nbac, nrgt, nlft, &
                nbot, ntop)
        END DO

        CALL stop_timer(410)
    END SUBROUTINE tstsca4


    SUBROUTINE tstsca4_grid(kk, jj, ii, qtu, qtv, qtw, t, u, v, w, g, bt, &
            ddx, ddy, ddz, rdx, rdy, rdz, sca, nfro, nbac, nrgt, nlft, &
            nbot, ntop)

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
        REAL(realk) :: gsca(kk)
        REAL(realk) :: adv, diff, area
        REAL(realk) :: gscamol, gtgmolp, gtgmoln

        ! Set INTENT(out) to zero
        qtu = 0.0
        qtv = 0.0
        qtw = 0.0

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
        DO i = 3-nfu, ii-3+nbu
            DO j = 3, jj-2
                ! Scalar diffusivity LES/DNS computation
                IF (iles == 1) THEN
                    DO k = 3, kk-2
                        gscamol = gmol/rho/sca%prmol
                        gtgmolp = (g(k, j, i) - gmol)/gmol
                        gtgmoln = (g(k, j, i+1) - gmol)/gmol

                        ! 1/Re * 1/Pr + 1/Re_t * 1/Pr_t:
                        gsca(k) = gscamol &
                            + (g(k, j, i+1) + g(k, j, i) - 2.0*gmol) / rho &
                            / (sca%prt(gtgmoln) + sca%prt(gtgmolp))

                        ! Limit gsca here MAX(...,0): no negative diffusion!
                        gsca(k) = MAX(gscamol, gsca(k))
                    END DO
                ELSE
                    DO k = 3, kk-2
                        gsca(k) = gmol/rho/sca%prmol
                    END DO
                END IF

                ! Final asembly
                DO k = 3, kk-2
                    ! Convective fluxes
                    ! It is assumed that the velocity field is already masked
                    ! with BU, BV, BW = no new masking necessary (!)
                    adv = (ddy(j)*ddz(k)) * u(k, j, i) &
                        * 0.5 * (t(k, j, i) + t(k, j, i+1))

                    ! Depending on the knowledge about the the cell and its
                    ! neighbours it is determined if faces are blocked (=0)
                    ! or open (=1)
                    area = bt(k, j, i)*bt(k, j, i+1)*(ddy(j)*ddz(k))
                    diff = -gsca(k)*rdx(i)*(t(k, j, i+1) - t(k, j, i))*area

                    ! Final result
                    qtu(k, j, i) = adv + diff
                END DO
            END DO
        END DO

        ! Y direction
        DO i = 3, ii-2
            DO j = 3-nrv, jj-3+nlv
                ! Scalar diffusivity LES/DNS computation
                IF (iles == 1) THEN
                    DO k = 3, kk-2
                        gscamol = gmol/rho/sca%prmol
                        gtgmolp = (g(k, j, i) - gmol)/gmol
                        gtgmoln = (g(k, j+1, i) - gmol)/gmol

                        ! 1/Re * 1/Pr + 1/Re_t * 1/Pr_t:
                        gsca(k) = gscamol &
                            + (g(k, j+1, i) + g(k, j, i) - 2.0*gmol) / rho &
                            / (sca%prt(gtgmoln) + sca%prt(gtgmolp))

                        ! Limit gsca here MAX(...,0): no negative diffusion!
                        gsca(k) = MAX(gscamol, gsca(k))
                    END DO
                ELSE
                    DO k = 3, kk-2
                        gsca(k) = gmol/rho/sca%prmol
                    END DO
                END IF

                ! Final asembly
                DO k = 3, kk-2
                    ! Convective fluxes
                    ! It is assumed that the velocity field is already masked
                    ! with BU, BV, BW = no new masking necessary (!)
                    adv = (ddx(i)*ddz(k)) * v(k, j, i) &
                        * 0.5 * (t(k, j, i) + t(k, j+1, i))

                    ! Depending on the knowledge about the the cell and its
                    ! neighbours it is determined if faces are blocked (=0)
                    ! or open (=1)
                    area = bt(k, j, i)*bt(k, j+1, i)*(ddx(i)*ddz(k))
                    diff = -gsca(k)*rdy(j)*(t(k, j+1, i) - t(k, j, i))*area

                    ! Final result
                    qtv(k, j, i) = adv + diff
                END DO
            END DO
        END DO

        ! Z direction
        DO i = 3, ii-2
            DO j = 3, jj-2
                ! Scalar diffusivity LES/DNS computation
                IF (iles == 1) THEN
                    DO k = 3-nbw, kk-3+ntw
                        gscamol = gmol/rho/sca%prmol
                        gtgmolp = (g(k, j, i) - gmol)/gmol
                        gtgmoln = (g(k+1, j, i) - gmol)/gmol

                        ! 1/Re * 1/Pr + 1/Re_t * 1/Pr_t:
                        gsca(k) = gscamol &
                            + (g(k+1, j, i) + g(k, j, i) - 2.0*gmol) / rho &
                            / (sca%prt(gtgmoln) + sca%prt(gtgmolp))

                        ! Limit gsca here MAX(...,0): no negative diffusion!
                        gsca(k) = MAX(gscamol, gsca(k))
                    END DO
                ELSE
                    DO k = 3-nbw, kk-3+ntw
                        gsca(k) = gmol/rho/sca%prmol
                    END DO
                END IF

                ! Final asembly
                DO k = 3-nbw, kk-3+ntw
                    ! Convective fluxes
                    ! It is assumed that the velocity field is already masked
                    ! with BU, BV, BW = no new masking necessary (!)
                    adv = (ddx(i)*ddy(j)) * w(k, j, i) &
                        * 0.5 * (t(k, j, i) + t(k+1, j, i))

                    ! Depending on the knowledge about the the cell and its
                    ! neighbours it is determined if faces are blocked (=0)
                    ! or open (=1)
                    area = bt(k, j, i)*bt(k+1, j, i)*(ddx(i)*ddy(j))
                    diff = -gsca(k)*rdy(j)*(t(k+1, j, i) - t(k, j, i))*area

                    ! Final result
                    qtw(k, j, i) = adv + diff
                END DO
            END DO
        END DO


        ! Special treatment at par boundaries
        ! Substraction of downwind and addition of upwind T-value
        ! to finally get an upwind scheme in case of flow towards coarse grid
        IF (nfro == 8) THEN
            i =  3
            DO j = 3, jj-2
                DO k = 3, kk-2
                    adv = (ddy(j)*ddz(k)) * (u(k, j, i) - ABS(u(k, j, i))) &
                        * 0.5 * 0.5 * (-t(k, j, i) + t(k, j, i+1))
                    qtu(k, j, i) = qtu(k, j, i) + adv
                END DO
            END DO
        END IF

        IF (nbac == 8) THEN
            i = ii-3
            DO j = 3, jj-2
                DO k = 3, kk-2
                    adv = (ddy(j)*ddz(k)) * (u(k, j, i) + ABS(u(k, j, i))) &
                        * 0.5 * 0.5 * (t(k, j, i) - t(k, j, i+1))
                    qtu(k, j, i) = qtu(k, j, i) + adv
                END DO
            END DO
        END IF

        IF (nrgt == 8) THEN
            DO i = 3, ii-2
                j = 3
                DO k = 3, kk-2
                    adv = (ddx(i)*ddz(k)) * (v(k, j, i) - ABS(v(k, j, i))) &
                        * 0.5 * 0.5 * (-t(k, j, i) + t(k, j+1, i) )
                    qtv(k, j, i) = qtv(k, j, i) + adv
                END DO
            END DO
        END IF

        IF (nlft == 8) THEN
            DO i = 3, ii-2
                j = jj-3
                DO k = 3, kk-2
                    adv = (ddx(i)*ddz(k)) * (v(k, j, i) + ABS(v(k, j, i))) &
                        * 0.5 * 0.5 * (t(k, j, i) - t(k, j+1, i))
                    qtv(k, j, i) = qtv(k, j, i) + adv
                END DO
            END DO
        END IF

        IF (nbot == 8) THEN
            DO i = 3, ii-2
                DO j = 3, jj-2
                    k = 3
                    adv = (ddx(i)*ddy(j)) * (w(k, j, i) - ABS(w(k, j, i))) &
                        * 0.5 * 0.5 * (-t(k, j, i) + t(k+1, j, i))
                    qtw(k, j, i) = qtw(k, j, i) + adv
                END DO
            END DO
        END IF

        IF (ntop == 8) THEN
            DO i = 3, ii-2
                DO j = 3, jj-2
                    k = kk-3
                    adv = (ddx(i)*ddy(j)) * (w(k, j, i) + ABS(w(k, j, i))) &
                        * 0.5 * 0.5 * (t(k, j, i) - t(k+1, j, i))
                    qtw(k, j, i) = qtw(k, j, i) + adv
                END DO
            END DO
        END IF
    END SUBROUTINE tstsca4_grid


    SUBROUTINE fluxbalance(qtt_f, qtu_f, qtv_f, qtw_f)
        ! Subroutine arguments
        TYPE(field_t), INTENT(inout) :: qtt_f
        TYPE(field_t), INTENT(in) :: qtu_f, qtv_f, qtw_f

        ! Local variables
        INTEGER(intk) :: i, igrid
        INTEGER(intk) :: kk, jj, ii
        TYPE(field_t), POINTER :: rddx_f, rddy_f, rddz_f
        REAL(realk), POINTER, CONTIGUOUS, DIMENSION(:, :, :) :: qtu, qtv, qtw, &
            qtt
        REAL(realk), POINTER, CONTIGUOUS:: rddx(:), rddy(:), rddz(:)

        CALL start_timer(411)

        CALL get_field(rddx_f, "RDDX")
        CALL get_field(rddy_f, "RDDY")
        CALL get_field(rddz_f, "RDDZ")

        DO i = 1, nmygrids
            igrid = mygrids(i)

            CALL get_mgdims(kk, jj, ii, igrid)

            CALL qtt_f%get_ptr(qtt, igrid)
            CALL qtu_f%get_ptr(qtu, igrid)
            CALL qtv_f%get_ptr(qtv, igrid)
            CALL qtw_f%get_ptr(qtw, igrid)

            CALL rddx_f%get_ptr(rddx, igrid)
            CALL rddy_f%get_ptr(rddy, igrid)
            CALL rddz_f%get_ptr(rddz, igrid)

            CALL fluxbalance_grid(kk, jj, ii, qtt, qtu, qtv, qtw, &
                rddx, rddy, rddz)
        END DO

        CALL stop_timer(411)
    END SUBROUTINE fluxbalance


    SUBROUTINE fluxbalance_grid(kk, jj, ii, qtt, qtu, qtv, qtw, &
            rddx, rddy, rddz)
        ! Subroutine arguments
        INTEGER(intk), INTENT(IN) :: kk, jj, ii
        REAL(realk), INTENT(OUT), DIMENSION(kk, jj, ii) :: qtt
        REAL(realk), INTENT(IN), DIMENSION(kk, jj, ii) :: qtu, qtv, qtw
        REAL(realk), INTENT(IN) :: rddx(ii), rddy(jj), rddz(kk)

        ! Local variables
        INTEGER(intk) :: i, j, k
        REAL(realk) :: netflux

        ! Set INTENT(out) to zero
        qtt = 0.0

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
    END SUBROUTINE fluxbalance_grid


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
