MODULE timeintegrate_scalar_mod
    USE core_mod
    USE ib_mod, ONLY: parent, ftoc, ib
    USE scacore_mod
    USE flow_mod, ONLY: ilesmodel, gmol, rho
    USE bound_scalar_mod
    USE itinfo_scalar_mod
    USE gc_scastencils_mod
    USE conn_mod, ONLY: conn

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
        TYPE(field_t), POINTER :: qtt, qtu, qtv, qtw

        IF (.NOT. solve_scalar) RETURN
        CALL start_timer(400)
        CALL start_timer(401)

        ! Local temporary storage ("scrap")
        CALL push_field(qtt, "QTT")
        CALL push_field(qtu, "QTU", istag=1)
        CALL push_field(qtv, "QTV", jstag=1)
        CALL push_field(qtw, "QTW", kstag=1)

        CALL stop_timer(401)

        CALL start_timer(402)
        ! In IRK 1, FRHS is zero, therefore we do not need to zeroize
        ! the dt field before each step
        CALL rkscheme%get_coeffs(frhs, fu, dtrk, dtrki, irk)

        DO l = 1, nsca
            ! Fetch fields
            CALL get_field(t, scalar(l)%name)
            CALL get_field(dt_f, "D"//TRIM(scalar(l)%name))
            CALL get_field(told, TRIM(scalar(l)%name)//"_OLD")

            ! Copy to "T_OLD" (not needed here but used for Boussinesq)
            CALL copy_arr(told%arr, t%arr)

            ! TSTSCA4 zeroize qtu, qtv, qtw before use internally
            CALL tstsca4(qtu, qtv, qtw, t, scalar(l))

            ! This operation apply boundary conditions to qtu, qtv, qtw ONLY!
            ! Does not modify t-field at all!
            DO ilevel = minlevel, maxlevel
                CALL parent(ilevel, qtu, qtv, qtw, device=.TRUE.)
                CALL apply_bound_scalar(ilevel, l, qtu, qtv, qtw, t)
            END DO

            ! fluxbalance zeroize qtt before use internally
            CALL fluxbalance(qtt, qtu, qtv, qtw)

            ! Additional source terms
            CALL sourceterm(qtt, scalar(l))

            ! Ghost cell "flux" boundary condition applied to qtt field
            IF (ib%type == "GHOSTCELL") THEN
                CALL set_scastencils("P", scalar(l), qtt=qtt)
            END IF

            ! dT_j = A_j*dT_(j-1) + QTT
            ! T_j = T_(j-1) + B_j*dT_j
            CALL rkstep(t%arr, dt_f%arr, qtt%arr, frhs, dt*fu)

            ! Mask blocked cells
            CALL maskbt(t)

            ! Ghost cell "value" boundary condition applied to t field
            IF (ib%type == "GHOSTCELL") THEN
                CALL conn(layers=2, s1=t, corners=.TRUE.)
                CALL set_scastencils("P", scalar(l), t=t)
            END IF

            DO ilevel = maxlevel, minlevel+1, -1
                CALL ftoc(ilevel, t, t, 'T', device=.TRUE.)
            END DO

            CALL conn(layers=2, s1=t, corners=.TRUE.)

            ! TODO: Fill ghost layers of T (maybe only at last IRK?)
        END DO
        CALL stop_timer(402)

        CALL pop_field(qtw)
        CALL pop_field(qtv)
        CALL pop_field(qtu)
        CALL pop_field(qtt)

        IF (irk == rkscheme%nrk) THEN
            CALL sync_scalars_to_host()
        END IF

        CALL stop_timer(400)
    END SUBROUTINE timeintegrate_scalar


    SUBROUTINE sync_scalars_to_host()
        INTEGER(intk) :: l
        TYPE(field_t), POINTER :: t

        DO l = 1, nsca
            CALL get_field(t, scalar(l)%name)
            CALL map_arr_from_device(t, message="from:scalar")
        END DO
    END SUBROUTINE sync_scalars_to_host


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
        TYPE(field_t), POINTER :: u_f, v_f, w_f, g_f
        TYPE(field_t), POINTER :: bt_f, ddx_f, ddy_f, ddz_f, rdx_f, rdy_f, rdz_f

        CALL start_timer(410)

        CALL zero_field_arr(qtu_f, device=.TRUE.)
        CALL zero_field_arr(qtv_f, device=.TRUE.)
        CALL zero_field_arr(qtw_f, device=.TRUE.)

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

        CALL tstsca4_impl(qtu_f%arr, qtv_f%arr, qtw_f%arr, t_f%arr, &
            u_f%arr, v_f%arr, w_f%arr, g_f%arr, bt_f%arr, ddx_f%arr, &
            ddy_f%arr, ddz_f%arr, rdx_f%arr, rdy_f%arr, rdz_f%arr, &
            sca%prmol, sca%kayscrawford, prturb)

        CALL stop_timer(410)
    END SUBROUTINE tstsca4


    SUBROUTINE tstsca4_impl(qtu, qtv, qtw, t, u, v, w, g, bt, ddx, ddy, &
            ddz, rdx, rdy, rdz, prmol, kayscrawford, prturb2)
        REAL(realk), INTENT(inout) :: qtu(*), qtv(*), qtw(*)
        REAL(realk), INTENT(in) :: t(*), u(*), v(*), w(*), g(*), bt(*)
        REAL(realk), INTENT(in) :: ddx(*), ddy(*), ddz(*)
        REAL(realk), INTENT(in) :: rdx(*), rdy(*), rdz(*)
        REAL(realk), INTENT(in) :: prmol, prturb2
        INTEGER(intk), INTENT(in) :: kayscrawford

        INTEGER(intk) :: i, igrid, ip3, ipx, ipy, ipz
        INTEGER(intk) :: kk, jj, ii, nfro, nbac, nrgt, nlft, nbot, ntop

        !$omp target teams distribute &
        !$omp& private(i, igrid, ip3, ipx, ipy, ipz, kk, jj, ii, nfro, &
        !$omp& nbac, nrgt, nlft, nbot, ntop)
        DO i = 1, nmygrids
            igrid = mygrids(i)
            CALL get_mgdims(kk, jj, ii, igrid)
            CALL get_mgbasb(nfro, nbac, nrgt, nlft, nbot, ntop, igrid)
            CALL get_ip3(ip3, igrid)
            CALL get_ip1x(ipx, igrid)
            CALL get_ip1y(ipy, igrid)
            CALL get_ip1z(ipz, igrid)

            !$omp parallel
            CALL tstsca4_grid(kk, jj, ii, qtu(ip3), qtv(ip3), qtw(ip3), &
                t(ip3), u(ip3), v(ip3), w(ip3), g(ip3), bt(ip3), ddx(ipx), &
                ddy(ipy), ddz(ipz), rdx(ipx), rdy(ipy), rdz(ipz), prmol, &
                kayscrawford, prturb2, nfro, nbac, nrgt, nlft, nbot, ntop)
            !$omp end parallel
        END DO
        !$omp end target teams distribute
    END SUBROUTINE tstsca4_impl


    SUBROUTINE tstsca4_grid(kk, jj, ii, qtu, qtv, qtw, t, u, v, w, g, bt, &
            ddx, ddy, ddz, rdx, rdy, rdz, prmol, kayscrawford, prturb2, &
            nfro, nbac, nrgt, nlft, nbot, ntop)
        !$omp declare target

        ! Subroutine arguments
        INTEGER(intk), INTENT(IN) :: kk, jj, ii
        REAL(realk), INTENT(INOUT), DIMENSION(kk, jj, ii) :: qtu, qtv, qtw
        REAL(realk), INTENT(IN), DIMENSION(kk, jj, ii) :: t, u, v, w, g, bt
        REAL(realk), INTENT(IN) :: ddx(ii), ddy(jj), ddz(kk), &
            rdx(ii), rdy(jj), rdz(kk)
        REAL(realk), INTENT(in) :: prmol, prturb2
        INTEGER(intk), INTENT(in) :: kayscrawford
        INTEGER(intk), INTENT(IN) :: nfro, nbac, nrgt, nlft, nbot, ntop

        ! Local variables
        INTEGER(intk) :: i, j, k
        INTEGER(intk) :: nfu, nbu, nrv, nlv, nbw, ntw
        INTEGER(intk) :: iles
        REAL(realk) :: gscak
        REAL(realk) :: adv, diff, area
        REAL(realk) :: gscamol, gtgmolp, gtgmoln

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
        !$omp do collapse(3) private(i, j, k, gscamol, gtgmolp, gtgmoln, &
        !$omp& adv, diff, area, gscak)
        DO i = 3-nfu, ii-3+nbu
            DO j = 3, jj-2
                DO k = 3, kk-2
                    ! Scalar diffusivity LES/DNS computation
                    IF (iles == 1) THEN
                        gscamol = gmol/rho/prmol
                        gtgmolp = (g(k, j, i) - gmol)/gmol
                        gtgmoln = (g(k, j, i+1) - gmol)/gmol

                        ! 1/Re * 1/Pr + 1/Re_t * 1/Pr_t:
                        gscak = gscamol &
                            + (g(k, j, i+1) + g(k, j, i) - 2.0*gmol) / rho &
                            / (scalar_prt(gtgmoln, prmol, kayscrawford, &
                            prturb2) + scalar_prt(gtgmolp, prmol, &
                            kayscrawford, prturb2))

                        ! Limit gsca here MAX(..., 0): no negative diffusion!
                        gscak = MAX(gscamol, gscak)
                    ELSE
                        gscak = gmol/rho/prmol
                    END IF

                    ! Convective fluxes
                    ! It is assumed that the velocity field is already masked
                    ! with BU, BV, BW = no new masking necessary (!)
                    adv = (ddy(j)*ddz(k)) * u(k, j, i) &
                        * 0.5 * (t(k, j, i) + t(k, j, i+1))

                    ! Depending on the knowledge about the the cell and its
                    ! neighbours it is determined if faces are blocked (=0)
                    ! or open (=1)
                    area = bt(k, j, i)*bt(k, j, i+1)*(ddy(j)*ddz(k))
                    diff = -gscak*rdx(i)*(t(k, j, i+1) - t(k, j, i))*area

                    ! Final result
                    qtu(k, j, i) = adv + diff
                END DO
            END DO
        END DO
        !$omp end do

        ! Y direction
        !$omp do collapse(3) private(i, j, k, gscamol, gtgmolp, gtgmoln, &
        !$omp& adv, diff, area, gscak)
        DO i = 3, ii-2
            DO j = 3-nrv, jj-3+nlv
                DO k = 3, kk-2
                    ! Scalar diffusivity LES/DNS computation
                    IF (iles == 1) THEN
                        gscamol = gmol/rho/prmol
                        gtgmolp = (g(k, j, i) - gmol)/gmol
                        gtgmoln = (g(k, j+1, i) - gmol)/gmol

                        ! 1/Re * 1/Pr + 1/Re_t * 1/Pr_t:
                        gscak = gscamol &
                            + (g(k, j+1, i) + g(k, j, i) - 2.0*gmol) / rho &
                            / (scalar_prt(gtgmoln, prmol, kayscrawford, &
                            prturb2) + scalar_prt(gtgmolp, prmol, &
                            kayscrawford, prturb2))

                        ! Limit gsca here MAX(..., 0): no negative diffusion!
                        gscak = MAX(gscamol, gscak)
                    ELSE
                        gscak = gmol/rho/prmol
                    END IF

                    ! Convective fluxes
                    ! It is assumed that the velocity field is already masked
                    ! with BU, BV, BW = no new masking necessary (!)
                    adv = (ddx(i)*ddz(k)) * v(k, j, i) &
                        * 0.5 * (t(k, j, i) + t(k, j+1, i))

                    ! Depending on the knowledge about the the cell and its
                    ! neighbours it is determined if faces are blocked (=0)
                    ! or open (=1)
                    area = bt(k, j, i)*bt(k, j+1, i)*(ddx(i)*ddz(k))
                    diff = -gscak*rdy(j)*(t(k, j+1, i) - t(k, j, i))*area

                    ! Final result
                    qtv(k, j, i) = adv + diff
                END DO
            END DO
        END DO
        !$omp end do

        ! Z direction
        !$omp do collapse(3) private(i, j, k, gscamol, gtgmolp, gtgmoln, &
        !$omp& adv, diff, area, gscak)
        DO i = 3, ii-2
            DO j = 3, jj-2
                DO k = 3-nbw, kk-3+ntw
                    ! Scalar diffusivity LES/DNS computation
                    IF (iles == 1) THEN
                        gscamol = gmol/rho/prmol
                        gtgmolp = (g(k, j, i) - gmol)/gmol
                        gtgmoln = (g(k+1, j, i) - gmol)/gmol

                        ! 1/Re * 1/Pr + 1/Re_t * 1/Pr_t:
                        gscak = gscamol &
                            + (g(k+1, j, i) + g(k, j, i) - 2.0*gmol) / rho &
                            / (scalar_prt(gtgmoln, prmol, kayscrawford, &
                            prturb2) + scalar_prt(gtgmolp, prmol, &
                            kayscrawford, prturb2))

                        ! Limit gsca here MAX(..., 0): no negative diffusion!
                        gscak = MAX(gscamol, gscak)
                    ELSE
                        gscak = gmol/rho/prmol
                    END IF

                    ! Convective fluxes
                    ! It is assumed that the velocity field is already masked
                    ! with BU, BV, BW = no new masking necessary (!)
                    adv = (ddx(i)*ddy(j)) * w(k, j, i) &
                        * 0.5 * (t(k, j, i) + t(k+1, j, i))

                    ! Depending on the knowledge about the the cell and its
                    ! neighbours it is determined if faces are blocked (=0)
                    ! or open (=1)
                    area = bt(k, j, i)*bt(k+1, j, i)*(ddx(i)*ddy(j))
                    diff = -gscak*rdz(k)*(t(k+1, j, i) - t(k, j, i))*area

                    ! Final result
                    qtw(k, j, i) = adv + diff
                END DO
            END DO
        END DO
        !$omp end do


        ! Special treatment at par boundaries
        ! Substraction of downwind and addition of upwind T-value
        ! to finally get an upwind scheme in case of flow towards coarse grid
        IF (nfro == 8) THEN
            i =  3
            !$omp do collapse(2) private(j, k, adv)
            DO j = 3, jj-2
                DO k = 3, kk-2
                    adv = (ddy(j)*ddz(k)) * (u(k, j, i) - ABS(u(k, j, i))) &
                        * 0.5 * 0.5 * (-t(k, j, i) + t(k, j, i+1))
                    qtu(k, j, i) = qtu(k, j, i) + adv
                END DO
            END DO
            !$omp end do
        END IF

        IF (nbac == 8) THEN
            i = ii-3
            !$omp do collapse(2) private(j, k, adv)
            DO j = 3, jj-2
                DO k = 3, kk-2
                    adv = (ddy(j)*ddz(k)) * (u(k, j, i) + ABS(u(k, j, i))) &
                        * 0.5 * 0.5 * (t(k, j, i) - t(k, j, i+1))
                    qtu(k, j, i) = qtu(k, j, i) + adv
                END DO
            END DO
            !$omp end do
        END IF

        IF (nrgt == 8) THEN
            !$omp do private(i, j, k, adv)
            DO i = 3, ii-2
                j = 3
                DO k = 3, kk-2
                    adv = (ddx(i)*ddz(k)) * (v(k, j, i) - ABS(v(k, j, i))) &
                        * 0.5 * 0.5 * (-t(k, j, i) + t(k, j+1, i))
                    qtv(k, j, i) = qtv(k, j, i) + adv
                END DO
            END DO
            !$omp end do
        END IF

        IF (nlft == 8) THEN
            !$omp do private(i, j, k, adv)
            DO i = 3, ii-2
                j = jj-3
                DO k = 3, kk-2
                    adv = (ddx(i)*ddz(k)) * (v(k, j, i) + ABS(v(k, j, i))) &
                        * 0.5 * 0.5 * (t(k, j, i) - t(k, j+1, i))
                    qtv(k, j, i) = qtv(k, j, i) + adv
                END DO
            END DO
            !$omp end do
        END IF

        IF (nbot == 8) THEN
            !$omp do collapse(2) private(i, j, k, adv)
            DO i = 3, ii-2
                DO j = 3, jj-2
                    k = 3
                    adv = (ddx(i)*ddy(j)) * (w(k, j, i) - ABS(w(k, j, i))) &
                        * 0.5 * 0.5 * (-t(k, j, i) + t(k+1, j, i))
                    qtw(k, j, i) = qtw(k, j, i) + adv
                END DO
            END DO
            !$omp end do
        END IF

        IF (ntop == 8) THEN
            !$omp do collapse(2) private(i, j, k, adv)
            DO i = 3, ii-2
                DO j = 3, jj-2
                    k = kk-3
                    adv = (ddx(i)*ddy(j)) * (w(k, j, i) + ABS(w(k, j, i))) &
                        * 0.5 * 0.5 * (t(k, j, i) - t(k+1, j, i))
                    qtw(k, j, i) = qtw(k, j, i) + adv
                END DO
            END DO
            !$omp end do
        END IF
    END SUBROUTINE tstsca4_grid


    PURE ELEMENTAL REAL(realk) FUNCTION scalar_prt(gtgmol, prmol, &
            kayscrawford_flag, prturb2)
        !$omp declare target
        REAL(realk), INTENT(in) :: gtgmol, prmol, prturb2
        INTEGER(intk), INTENT(in) :: kayscrawford_flag

        REAL(realk) :: kayscrawford

        IF (kayscrawford_flag == 0) THEN
            scalar_prt = prturb2
        ELSE IF (gtgmol > 0.0) THEN
            kayscrawford = 0.5882 + 0.228*gtgmol &
                - 0.0441*gtgmol**2*(1.0 - EXP(-5.165/gtgmol))
            scalar_prt = kayscrawford
        ELSE
            scalar_prt = prmol
        END IF
    END FUNCTION scalar_prt


    SUBROUTINE fluxbalance(qtt_f, qtu_f, qtv_f, qtw_f)
        ! Subroutine arguments
        TYPE(field_t), INTENT(inout) :: qtt_f
        TYPE(field_t), INTENT(in) :: qtu_f, qtv_f, qtw_f

        ! Local variables
        TYPE(field_t), POINTER :: rddx_f, rddy_f, rddz_f

        CALL start_timer(411)

        CALL zero_field_arr(qtt_f, device=.TRUE.)

        CALL get_field(rddx_f, "RDDX")
        CALL get_field(rddy_f, "RDDY")
        CALL get_field(rddz_f, "RDDZ")

        CALL fluxbalance_impl(qtt_f%arr, qtu_f%arr, qtv_f%arr, qtw_f%arr, &
            rddx_f%arr, rddy_f%arr, rddz_f%arr)

        CALL stop_timer(411)
    END SUBROUTINE fluxbalance


    SUBROUTINE fluxbalance_impl(qtt, qtu, qtv, qtw, rddx, rddy, rddz)
        REAL(realk), INTENT(inout) :: qtt(*)
        REAL(realk), INTENT(in) :: qtu(*), qtv(*), qtw(*)
        REAL(realk), INTENT(in) :: rddx(*), rddy(*), rddz(*)

        INTEGER(intk) :: i, igrid, kk, jj, ii, ip3, ipx, ipy, ipz

        !$omp target teams distribute &
        !$omp& private(i, igrid, kk, jj, ii, ip3, ipx, ipy, ipz)
        DO i = 1, nmygrids
            igrid = mygrids(i)
            CALL get_mgdims(kk, jj, ii, igrid)
            CALL get_ip3(ip3, igrid)
            CALL get_ip1x(ipx, igrid)
            CALL get_ip1y(ipy, igrid)
            CALL get_ip1z(ipz, igrid)

            !$omp parallel
            CALL fluxbalance_grid(kk, jj, ii, qtt(ip3), qtu(ip3), &
                qtv(ip3), qtw(ip3), rddx(ipx), rddy(ipy), rddz(ipz))
            !$omp end parallel
        END DO
        !$omp end target teams distribute
    END SUBROUTINE fluxbalance_impl


    SUBROUTINE fluxbalance_grid(kk, jj, ii, qtt, qtu, qtv, qtw, &
            rddx, rddy, rddz)
        !$omp declare target
        ! Subroutine arguments
        INTEGER(intk), INTENT(IN) :: kk, jj, ii
        REAL(realk), INTENT(INOUT), DIMENSION(kk, jj, ii) :: qtt
        REAL(realk), INTENT(IN), DIMENSION(kk, jj, ii) :: qtu, qtv, qtw
        REAL(realk), INTENT(IN) :: rddx(ii), rddy(jj), rddz(kk)

        ! Local variables
        INTEGER(intk) :: i, j, k
        REAL(realk) :: netflux

        !$omp do collapse(3) private(i, j, k, netflux)
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
        !$omp end do
    END SUBROUTINE fluxbalance_grid


    SUBROUTINE sourceterm(qtt_f, sca)
        ! Subroutine arguments
        TYPE(field_t), INTENT(inout) :: qtt_f
        TYPE(scalar_t), INTENT(in), TARGET :: sca

        ! Local variables
        INTEGER(intk) :: isource, nsource
        TYPE(scalar_source_t), POINTER :: source
        TYPE(field_t), POINTER :: field_f

        CALL start_timer(413)

        nsource = SIZE(sca%sources)

        DO isource = 1, nsource
            source => sca%sources(isource)

            IF (LEN_TRIM(source%field) == 0) THEN
                CALL sourceterm_const_impl(qtt_f%arr, source%value)
            ELSE
                CALL get_field(field_f, source%field)
                IF (field_f%istag /= 0 .OR. field_f%jstag /= 0 .OR. &
                        field_f%kstag /= 0) THEN
                    CALL errr(__FILE__, __LINE__)
                END IF
                CALL sourceterm_field_impl(qtt_f%arr, field_f%arr, &
                    source%value)
            END IF
        END DO

        CALL stop_timer(413)
    END SUBROUTINE sourceterm


    SUBROUTINE sourceterm_const_impl(qtt, sourceval)
        REAL(realk), INTENT(inout) :: qtt(*)
        REAL(realk), INTENT(in) :: sourceval

        INTEGER(intk) :: i, igrid, kk, jj, ii, ip3

        !$omp target teams distribute private(i, igrid, kk, jj, ii, ip3)
        DO i = 1, nmygrids
            igrid = mygrids(i)
            CALL get_mgdims(kk, jj, ii, igrid)
            CALL get_ip3(ip3, igrid)
            !$omp parallel
            CALL sourceterm_const(kk, jj, ii, qtt(ip3), sourceval)
            !$omp end parallel
        END DO
        !$omp end target teams distribute
    END SUBROUTINE sourceterm_const_impl


    SUBROUTINE sourceterm_field_impl(qtt, field, sourceval)
        REAL(realk), INTENT(inout) :: qtt(*)
        REAL(realk), INTENT(in) :: field(*)
        REAL(realk), INTENT(in) :: sourceval

        INTEGER(intk) :: i, igrid, kk, jj, ii, ip3

        !$omp target teams distribute private(i, igrid, kk, jj, ii, ip3)
        DO i = 1, nmygrids
            igrid = mygrids(i)
            CALL get_mgdims(kk, jj, ii, igrid)
            CALL get_ip3(ip3, igrid)
            !$omp parallel
            CALL sourceterm_field(kk, jj, ii, qtt(ip3), field(ip3), &
                sourceval)
            !$omp end parallel
        END DO
        !$omp end target teams distribute
    END SUBROUTINE sourceterm_field_impl


    SUBROUTINE sourceterm_const(kk, jj, ii, qtt, sourceval)
        !$omp declare target
        ! Subroutine arguments
        INTEGER(intk), INTENT(IN) :: kk, jj, ii
        REAL(realk), INTENT(INOUT), DIMENSION(kk, jj, ii) :: qtt
        REAL(realk), INTENT(IN) :: sourceval

        ! Local variables
        INTEGER(intk) :: i, j, k

        !$omp do collapse(3) private(i, j, k)
        DO i = 3, ii-2
            DO j = 3, jj-2
                DO k = 3, kk-2
                    qtt(k, j, i) = qtt(k, j, i) + sourceval
                END DO
            END DO
        END DO
        !$omp end do
    END SUBROUTINE sourceterm_const


    SUBROUTINE sourceterm_field(kk, jj, ii, qtt, field, sourceval)
        !$omp declare target
        ! Subroutine arguments
        INTEGER(intk), INTENT(IN) :: kk, jj, ii
        REAL(realk), INTENT(INOUT) :: qtt(kk, jj, ii)
        REAL(realk), INTENT(IN) :: field(kk, jj, ii)
        REAL(realk), INTENT(IN) :: sourceval

        ! Local variables
        INTEGER(intk) :: i, j, k

        !$omp do collapse(3) private(i, j, k)
        DO i = 3, ii-2
            DO j = 3, jj-2
                DO k = 3, kk-2
                    qtt(k, j, i) = qtt(k, j, i) + sourceval*field(k, j, i)
                END DO
            END DO
        END DO
        !$omp end do
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
