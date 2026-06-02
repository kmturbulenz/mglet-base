MODULE pressuresolver_mod
    USE bound_flow_mod
    USE core_mod
    USE flowcore_mod
    USE ib_mod
    USE itinfo_mod, ONLY: itinfo_sample
    USE plog_mod
    USE laplacephi_mod
    USE sip_hyperplane_mod, sipiter1_hp => sipiter1, sipiter2_hp => sipiter2
    USE sip_classic_mod, sipiter1_cl => sipiter1, sipiter2_cl => sipiter2
    USE sor_mod

    IMPLICIT NONE (type, external)
    PRIVATE

    ! Type of pressure solver
    !   0 : Hyperplane SIP
    !   1 : SIP on coarsest level, then SOR on subsequent levels
    !   2 : Classic SIP
    INTEGER(intk) :: ityp

    ! Number of inner iterations/sweeps in pressure solver
    INTEGER(intk) :: ninner

    ! Maximum number of outer iterations allowed
    INTEGER(intk) :: nouter

    ! Minimum number of outer iterations to always run independent of
    ! divergence
    INTEGER(intk) :: nouter_min

    ! Convergence threshold
    REAL(realk) :: epcorr = 0.0

    ! Loglevel:
    !   0: no special logging
    !   1: print number of iterations and final divergence after each mgpoisl
    !      (i.e. once per RK flow substep)
    !   2: Additionally prints the max residual per level after each outer
    !      pressure iteration
    !   3: Also compute and print the max residual before any iteration is done.
    INTEGER(intk), PROTECTED :: loglevel = 0

    ! A-versions: simple versions not considering the BP field
    ! B-versions: IB versions using the BP field

    INTERFACE pressureftocone
        MODULE PROCEDURE :: pressureftocone_A, pressureftocone_B
    END INTERFACE pressureftocone

    ! Bound operation
    TYPE, EXTENDS(bound_t) :: bound_pressure_t
    CONTAINS
        PROCEDURE, NOPASS :: front => bfront
        PROCEDURE, NOPASS :: back => bfront
        PROCEDURE, NOPASS :: right => bright
        PROCEDURE, NOPASS :: left => bright
        PROCEDURE, NOPASS :: bottom => bbottom
        PROCEDURE, NOPASS :: top => bbottom
    END TYPE bound_pressure_t
    TYPE(bound_pressure_t) :: bound_pressure

    PUBLIC :: init_pressuresolver, finish_pressuresolver, mgpoisl

CONTAINS

    SUBROUTINE init_pressuresolver()
        ! Subroutine arguments
        ! none...

        ! Local variables
        TYPE(config_t) :: psolveconf
        CHARACTER(len=16) :: type
        REAL(realk) :: omg

        ! Required values
        CALL fort7%get(psolveconf, "/flow/pressuresolver")

        CALL psolveconf%get_value("/ninner", ninner, default_value=5)
        ! Allowing zero pressure iterations are useful for debugging purposes
        IF (ninner < 0) THEN
            WRITE(*, '("Invalid number of ninner: ", I15)') ninner
            CALL errr(__FILE__, __LINE__)
        END IF

        CALL psolveconf%get_value("/nouter", nouter, default_value=9)
        ! Allowing zero pressure iterations are useful for debugging purposes
        IF (nouter < 0) THEN
            WRITE(*, '("Invalid number of nouter: ", I15)') nouter
            CALL errr(__FILE__, __LINE__)
        END IF

        CALL psolveconf%get_value("/nouter_min", nouter_min, default_value=0)
        IF (nouter_min < 0) THEN
            WRITE(*, '("Invalid number of nouter_min: ", I15)') nouter_min
            CALL errr(__FILE__, __LINE__)
        END IF

        CALL psolveconf%get_value("/loglevel", loglevel, default_value=0)
        IF (loglevel < 0) THEN
            WRITE(*, '("Invalid number of loglevel: ", I15)') loglevel
            CALL errr(__FILE__, __LINE__)
        END IF
        IF (loglevel > 0) THEN
            CALL init_plog(dcont)
        END IF

        ! There are no way we can set a default value for this - so this
        ! must always be provided
        CALL psolveconf%get_value("/epcorr", epcorr)
        IF (epcorr < 0.0) THEN
            WRITE(*, '("Invalid epcorr: ", F15.7)') epcorr
            CALL errr(__FILE__, __LINE__)
        END IF

        CALL psolveconf%get_value("/type", type, default_value="sip")
        SELECT CASE (LOWER(type))
        CASE("sip")
            ityp = 0
        CASE ("sor")
            ityp = 1
        CASE("classic")
            ityp = 2
        CASE DEFAULT
            WRITE(*, '("Invalid pressure solver type: ", A)') type
            CALL errr(__FILE__, __LINE__)
        END SELECT

        ! Always initialize SIP - it's always used at the coarsest level
        CALL set_field("SIPLW")
        CALL set_field("SIPLS")
        CALL set_field("SIPLB")
        CALL set_field("SIPUE")
        CALL set_field("SIPUN")
        CALL set_field("SIPUT")
        CALL set_field("SIPLPR")


        IF (ityp == 2) THEN
            ! Using the classic SIP solver
            CALL sip_classic_init()

        ELSE
            ! Using the hyperplane solver
            CALL sip_hyperplane_init()
        END IF

        ! Initialize SOR if enabled
        IF (ityp == 1) THEN
            CALL psolveconf%get_value("/omega", omg, default_value=1.1)
            IF (omg <= 0.0 .OR. omg >= 2.0) THEN
                WRITE(*, '("Invalid omega: ", F15.7)') omg
                CALL errr(__FILE__, __LINE__)
            END IF

            CALL sor_init(omg)
        END IF

    END SUBROUTINE init_pressuresolver


    SUBROUTINE finish_pressuresolver()
        CALL finish_plog()

        IF (ityp == 2) THEN
            CALL sip_classic_finish()
        ELSE
            CALL sip_hyperplane_finish()
        END IF

    END SUBROUTINE finish_pressuresolver


    SUBROUTINE mgpoisl(u, v, w, p, dt, ittot, irk)
        USE MPI_f08

        ! Subroutine arguments
        TYPE(field_t), INTENT(inout) :: u
        TYPE(field_t), INTENT(inout) :: v
        TYPE(field_t), INTENT(inout) :: w
        TYPE(field_t), INTENT(inout) :: p
        REAL(realk), INTENT(in) :: dt
        INTEGER(intk), INTENT(in) :: ittot
        INTEGER(intk), INTENT(in) :: irk

        ! Local variables
        TYPE(field_t), POINTER :: dp, hilf, rhs, res
        TYPE(field_t), POINTER :: bp
        INTEGER(intk) :: ilevel, ipcount, ipc, i
        REAL(realk) :: prefak, maxrhs, maxrhsall
        REAL(realk), ALLOCATABLE :: maxrhslvl(:)

        CALL start_timer(320)

        ALLOCATE(maxrhslvl(minlevel:maxlevel))

        NULLIFY(bp)
        IF (ib%type /= "NONE") THEN
            CALL get_field(bp, "BP")
        END IF

        ! All of these fields are initialized to zero automatically
        CALL push_field(dp, "DP")
        CALL push_field(hilf, "HILF")
        CALL push_field(rhs, "RHS")
        CALL push_field(res, "RES")

        ! laplace(dp) = prefak * div(u) is the underlying equation
        prefak = rho/dt
        CALL ib%divcal(rhs, u, v, w, prefak)

        DO ilevel = maxlevel, minlevel, -1
            CALL ftoc(ilevel, rhs, rhs, 'R')
        END DO

        ! For debug logging
        IF (loglevel >= 3) THEN
            CALL maxabscal(maxrhs, maxrhslvl, rhs)
            CALL MPI_Allreduce(MPI_IN_PLACE, maxrhslvl, &
                maxlevel-minlevel+1, mglet_mpi_real, MPI_MAX, &
                MPI_COMM_WORLD)
            DO ilevel = minlevel, maxlevel
                CALL sample_plog(ilevel, 0, maxrhslvl(ilevel))
            END DO
            CALL print_plog(ittot, irk, 0)
        END IF

        ipc = 0
        outer: DO ipcount = 1, nouter
            ! Inner pressure iterations
            ! HINT: 'res' is passed into mgpoisit as a temporary storage!!
            CALL start_timer(322)
            DO ilevel = minlevel, maxlevel
                CALL ctof(ilevel, hilf%arr, hilf%arr)
                CALL parent(ilevel, s1=hilf)
                CALL mgpoisit(ilevel, hilf, rhs, res, bp)
            END DO
            CALL stop_timer(322)

            ! --- intermediate state ---
            ! every grid level has an inner solution
            ! stored in hilf

            ! vom feinsten zum groebsten (fine to coarse)
            ! hilf <- hilf^f
            ! dp = dp + hilf
            ! wenn ltst: rhs = rhs - laplace(dp)
            ! rhs  <-  rhs^f
            ! vom feinsten zum groebsten (fine to coarse)
            ! hilf = 0.0
            DO ilevel = maxlevel, minlevel, -1
                CALL ftoc(ilevel, hilf, hilf, 'P')
            END DO

            ! --- intermediate state ---
            ! every grid level has the best solution
            ! for hilf retrieved from the locally
            ! finest grid.

            ! Connect needed due to prior ftoc, since this does not do
            ! anything on the finest level, no need to connect finest level
            ! either.
            CALL connect(layers=1, s1=hilf)

            ! res <- laplace(hilf)
            CALL laplacephi(res, hilf, bp)
            ! rhs <- rhs + res
            CALL rescal(rhs, res)

            DO ilevel = maxlevel, minlevel+1, -1
                CALL ftoc(ilevel, rhs, rhs, 'R')
            END DO

            ! Max of RHS scaled according to levels
            CALL maxabscal(maxrhs, maxrhslvl, rhs)

            ! dp = dp + hilf
            dp%arr = dp%arr + hilf%arr
            hilf%arr = 0.0
            ipc = ipc + ninner

            ! Pressure solver debug logging
            IF (loglevel >= 2) THEN
                CALL MPI_Allreduce(MPI_IN_PLACE, maxrhslvl, &
                    maxlevel-minlevel+1, mglet_mpi_real, MPI_MAX, &
                    MPI_COMM_WORLD)
                DO ilevel = minlevel, maxlevel
                    CALL sample_plog(ilevel, ninner, maxrhslvl(ilevel))
                END DO
                CALL print_plog(ittot, irk, ipcount)
            END IF

            ! Check outer iteration stopping criterion
            IF (ipcount >= nouter_min) THEN
                CALL MPI_Allreduce(maxrhs, maxrhsall, 1, mglet_mpi_real, &
                    MPI_MAX, MPI_COMM_WORLD)

                IF (maxrhsall/prefak < epcorr) THEN
                    EXIT outer
                END IF
            END IF
        END DO outer

        ! For general.log and screen output. Sample number of pressure
        ! iterations. Currently the same for all grids...
        DO i = 1, nmygrids
            CALL itinfo_sample(mygrids(i), ipc=ipc)
        END DO

        ! Pressure solver debug logging
        IF (loglevel == 1) THEN
            CALL MPI_Allreduce(MPI_IN_PLACE, maxrhslvl, maxlevel-minlevel+1, &
                mglet_mpi_real, MPI_MAX, MPI_COMM_WORLD)
            DO ilevel = minlevel, maxlevel
                CALL sample_plog(ilevel, ipc, maxrhslvl(ilevel))
            END DO
            CALL print_plog(ittot, irk, ipcount-1)
        END IF

        ! --- intermediate state ---
        ! The outer iteration has been left after
        ! a value of dp was found that leads to a acceptably small residual
        DO ilevel = minlevel, maxlevel
            CALL parent(ilevel, s1=dp)
            CALL bound_pressure%bound(ilevel, dp, bp)
        END DO

        ! Pressure correction: P = P + dtrk/rho*DP
        ! Velocity fields are modified and become solenoidal based on DP
        CALL mgpcorr(u, v, w, p, dp, dt/rho, bp)
        DO ilevel = maxlevel, minlevel, -1
            CALL ftoc(ilevel, u, v, w, p)
        END DO

        ! All levels (coarse to fine)
        ! Propagation of the solution to neighbours and childs
        ! The order of the calls is crucial:
        ! - First parent and boundmg, to fill the ghost layers at PAR-boundaries
        ! - Second to connect inside the level. Now also the correct information
        !   in the the ghost layers at PAR-boundaries is used in connect.
        DO ilevel = minlevel, maxlevel
            CALL parent(ilevel, u, v, w, p)
            CALL bound_flow%bound(ilevel, u, v, w, p)
            CALL connect(ilevel, 2, v1=u, v2=v, v3=w, s1=p, corners=.TRUE.)
        END DO

        CALL pop_field(res)
        CALL pop_field(rhs)
        CALL pop_field(hilf)
        CALL pop_field(dp)

        DEALLOCATE(maxrhslvl)
        CALL stop_timer(320)
    END SUBROUTINE mgpoisl


    ! 'res' is a temporary storage for the SIP algorithm,
    SUBROUTINE mgpoisit(ilevel, dp, rhs, res, bp)
        ! Subroutine arguments
        INTEGER(intk), INTENT(in) :: ilevel
        TYPE(field_t), INTENT(inout) :: dp
        TYPE(field_t), INTENT(in) :: rhs
        TYPE(field_t), INTENT(inout) :: res
        TYPE(field_t), INTENT(inout), OPTIONAL :: bp

        ! Local variables
        INTEGER(intk) :: iloop
        TYPE(field_t), POINTER :: gsaw, gsae, gsas, gsan, gsab, gsat, &
            gsap, gsrap
        TYPE(field_t), POINTER :: siplw, sipls, siplb, sipue, sipun, siput, &
            siplpr

        CALL start_timer(321)

        CALL get_field(gsaw, "GSAW")
        CALL get_field(gsae, "GSAE")
        CALL get_field(gsas, "GSAS")
        CALL get_field(gsan, "GSAN")
        CALL get_field(gsab, "GSAB")
        CALL get_field(gsat, "GSAT")
        CALL get_field(gsap, "GSAP")
        IF (ityp == 1) THEN
            ! Relax solver needs this as well
            CALL get_field(gsrap, "SOR_RAP")
        END IF

        ! Getting the adapted coefficients for SIP
        ! (same fields used for both classic and hyperplane SIP)
        CALL get_field(siplw, "SIPLW")
        CALL get_field(sipls, "SIPLS")
        CALL get_field(siplb, "SIPLB")
        CALL get_field(sipue, "SIPUE")
        CALL get_field(sipun, "SIPUN")
        CALL get_field(siput, "SIPUT")
        CALL get_field(siplpr, "SIPLPR")

        DO iloop = 1, ninner
            CALL bound_pressure%bound(ilevel, dp, bp)

            IF (ityp == 1 .AND. ilevel > minlevel) THEN
                ! The SOR relaxation is usually not efficient at the
                ! coarsest level, hence only apply at the finer levels
                CALL sor(ilevel, dp, rhs, gsaw, gsae, gsas, gsan, gsab, gsat, &
                    gsrap, bp)
            ELSE
                ! Use the SIP solver
                CALL sip(ilevel, iloop, dp, res, rhs, siplw, sipls, siplb, &
                    sipue, sipun, siput, siplpr, bp)
            END IF

            CALL connect(ilevel, 1, s1=dp)
        END DO

        CALL bound_pressure%bound(ilevel, dp, bp)

        CALL stop_timer(321)
    END SUBROUTINE mgpoisit


    SUBROUTINE sip(ilevel, iloop, dp, res, rhs, siplw, sipls, siplb, &
            sipue, sipun, siput, siplpr, bp)
        ! Subroutine arguments
        INTEGER(intk), INTENT(in) :: ilevel
        INTEGER(intk), INTENT(in) :: iloop
        TYPE(field_t), INTENT(inout) :: dp
        TYPE(field_t), INTENT(inout) :: res
        TYPE(field_t), INTENT(in) :: rhs
        TYPE(field_t), INTENT(in) :: siplw
        TYPE(field_t), INTENT(in) :: sipls
        TYPE(field_t), INTENT(in) :: siplb
        TYPE(field_t), INTENT(in) :: sipue
        TYPE(field_t), INTENT(in) :: sipun
        TYPE(field_t), INTENT(in) :: siput
        TYPE(field_t), INTENT(in) :: siplpr
        TYPE(field_t), INTENT(in), OPTIONAL :: bp

        ! Local variables
        INTEGER(intk) :: i, igrid
        INTEGER(intk) :: kk, jj, ii
        REAL(realk), POINTER, CONTIGUOUS :: lw(:, :, :), ls(:, :, :), &
            lb(:, :, :), ue(:, :, :), un(:, :, :), ut(:, :, :), &
            lpr(:, :, :)
        REAL(realk), POINTER, CONTIGUOUS :: dp_p(:, :, :), res_p(:, :, :), &
            rhs_p(:, :, :)
        INTEGER(ifk), CONTIGUOUS, POINTER :: mip_ptr(:), idx_ptr(:)

        CALL laplacephi_level(ilevel, res, dp, bp)

        DO i = 1, nmygridslvl(ilevel)
            igrid = mygridslvl(i, ilevel)
            CALL get_mgdims(kk, jj, ii, igrid)

            CALL res%get_ptr(res_p, igrid)
            CALL rhs%get_ptr(rhs_p, igrid)

            CALL siplw%get_ptr(lw, igrid)
            CALL sipls%get_ptr(ls, igrid)
            CALL siplb%get_ptr(lb, igrid)
            CALL siplpr%get_ptr(lpr, igrid)

            IF (ityp == 2) THEN
                CALL sipiter1_cl(kk, jj, ii, rhs_p, res_p, lw, ls, lb, lpr)
            ELSE
                CALL get_grid3_ifk_linear(mip_ptr, mip_hp_f, igrid)
                CALL get_grid3_ifk_linear(idx_ptr, idx_hp_f, igrid)
                CALL sipiter1_hp(kk, jj, ii, rhs_p, res_p, lw, ls, lb, lpr, &
                    mip_ptr, idx_ptr)
            END IF
        END DO

        IF (iloop < ninner) THEN
            CALL connect(ilevel, 1, s1=res)
        ELSE
            CALL connect(ilevel, 1, s1=res, forward=-1)
        END IF

        DO i = 1, nmygridslvl(ilevel)
            igrid = mygridslvl(i, ilevel)
            CALL get_mgdims(kk, jj, ii, igrid)

            CALL dp%get_ptr(dp_p, igrid)
            CALL res%get_ptr(res_p, igrid)

            CALL sipue%get_ptr(ue, igrid)
            CALL sipun%get_ptr(un, igrid)
            CALL siput%get_ptr(ut, igrid)

            IF (ityp == 2) THEN
                CALL sipiter2_cl(kk, jj, ii, dp_p, res_p, ue, un, ut)
            ELSE
                CALL get_grid3_ifk_linear(mip_ptr, mip_hp_f, igrid)
                CALL get_grid3_ifk_linear(idx_ptr, idx_hp_f, igrid)
                CALL sipiter2_hp(kk, jj, ii, dp_p, res_p, ue, un, ut, &
                    mip_ptr, idx_ptr)
            END IF
        END DO
    END SUBROUTINE sip


    SUBROUTINE bfront(igrid, iface, ibocd, ctyp, f1, f2, f3, f4, timeph)
        ! Subroutine arguments
        INTEGER(intk), INTENT(in) :: igrid, iface, ibocd
        CHARACTER(len=*), INTENT(in) :: ctyp
        TYPE(field_t), INTENT(inout) :: f1
        TYPE(field_t), INTENT(inout), OPTIONAL :: f2, f3, f4
        REAL(realk), INTENT(in), OPTIONAL :: timeph

        ! Local variables
        INTEGER(intk) :: kk, jj, ii
        INTEGER(intk) :: k, j, i, i2, i3, i4, istag2
        INTEGER(intk) :: m, n
        REAL(realk) :: pcnew, bpc, fak
        REAL(realk) :: sb11, sb12, sb13, sb14
        REAL(realk), POINTER, CONTIGUOUS :: p(:, :, :)
        REAL(realk), POINTER, CONTIGUOUS :: pbuffer(:, :)
        REAL(realk), POINTER, CONTIGUOUS :: bp(:, :, :)
        REAL(realk), POINTER, CONTIGUOUS :: dx(:), dy(:), dz(:)
        REAL(realk), POINTER, CONTIGUOUS :: ddx(:), ddy(:), ddz(:)

        ! Only works on PAR boundaries, should do nothing otherwise
        IF (ctyp /= 'PAR') RETURN

        CALL f1%get_ptr(p, igrid)
        CALL f1%get_buffer(pbuffer, igrid, iface)
        CALL get_mgdims(kk, jj, ii, igrid)

        CALL get_fieldptr(dx, "DX", igrid)
        CALL get_fieldptr(dy, "DY", igrid)
        CALL get_fieldptr(dz, "DZ", igrid)
        CALL get_fieldptr(ddx, "DDX", igrid)
        CALL get_fieldptr(ddy, "DDY", igrid)
        CALL get_fieldptr(ddz, "DDZ", igrid)

        SELECT CASE (iface)
        CASE (1)
            ! Front
            i2 = 2
            i3 = 3
            i4 = 4
            istag2 = 2
        CASE (2)
            ! Back
            i2 = ii - 1
            i3 = ii - 2
            i4 = ii - 3
            istag2 = ii - 2
        CASE DEFAULT
            CALL errr(__FILE__, __LINE__)
        END SELECT

        i = MIN(i3, i4)
        IF (PRESENT(f2)) THEN
            CALL f2%get_ptr(bp, igrid)

            DO j = 3, jj-2, 2
                DO k = 3, kk-2, 2
                    CALL pressureftocone(k, j, i, kk, jj, ii, p, bp, &
                        ddx, ddy, ddz, pcnew, bpc)
                    IF (bpc < 0.5) pcnew = pbuffer(k, j)

                    sb11 = bp(k, j, i2)*bp(k, j, i3)
                    sb12 = bp(k, j+1, i2)*bp(k, j+1, i3)
                    sb13 = bp(k+1, j, i2)*bp(k+1, j, i3)
                    sb14 = bp(k+1, j+1, i2)*bp(k+1, j+1, i3)

                    fak = (sb11*ddy(j)*ddz(k) + sb12*ddy(j+1)*ddz(k) &
                        + sb13*ddy(j)*ddz(k+1) + sb14*ddy(j+1)*ddz(k+1)) &
                        /((ddy(j)+ddy(j+1))*(ddz(k)+ddz(k+1)))
                    IF (fak < 0.1) fak = 1.0
                    fak = 1.0/fak

                    DO m = 0, 1
                        DO n = 0, 1
                            p(k+n, j+m, i2) = p(k+n, j+m, i3) &
                                + fak*dx(istag2)/(ddx(i3)+ddx(i2)) &
                                *(pbuffer(k, j) - pcnew)
                        END DO
                    END DO
                END DO
            END DO
        ELSE
            DO j = 3, jj-2, 2
                DO k = 3, kk-2, 2
                    CALL pressureftocone(k, j, i, kk, jj, ii, p, &
                        ddx, ddy, ddz, pcnew, bpc)
                    DO m = 0, 1
                        DO n = 0, 1
                            p(k+n, j+m, i2) = p(k+n, j+m, i3) &
                                + dx(istag2)/(ddx(i3)+ddx(i2)) &
                                *(pbuffer(k, j) - pcnew)
                        END DO
                    END DO
                END DO
            END DO
        END IF
    END SUBROUTINE bfront


    SUBROUTINE bright(igrid, iface, ibocd, ctyp, f1, f2, f3, f4, timeph)
        ! Subroutine arguments
        INTEGER(intk), INTENT(in) :: igrid, iface, ibocd
        CHARACTER(len=*), INTENT(in) :: ctyp
        TYPE(field_t), INTENT(inout) :: f1
        TYPE(field_t), INTENT(inout), OPTIONAL :: f2, f3, f4
        REAL(realk), INTENT(in), OPTIONAL :: timeph

        ! Local variables
        INTEGER(intk) :: kk, jj, ii
        INTEGER(intk) :: k, j, i, j2, j3, j4, jstag2
        INTEGER(intk) :: n, l
        REAL(realk) :: pcnew, bpc, fak
        REAL(realk) :: sb11, sb12, sb13, sb14
        REAL(realk), POINTER, CONTIGUOUS :: p(:, :, :)
        REAL(realk), POINTER, CONTIGUOUS :: pbuffer(:, :)
        REAL(realk), POINTER, CONTIGUOUS :: bp(:, :, :)
        REAL(realk), POINTER, CONTIGUOUS :: dx(:), dy(:), dz(:)
        REAL(realk), POINTER, CONTIGUOUS :: ddx(:), ddy(:), ddz(:)

        ! Only works on PAR boundaries, should do nothing otherwise
        IF (ctyp /= 'PAR') RETURN

        CALL f1%get_ptr(p, igrid)
        CALL f1%get_buffer(pbuffer, igrid, iface)
        CALL get_mgdims(kk, jj, ii, igrid)

        CALL get_fieldptr(dx, "DX", igrid)
        CALL get_fieldptr(dy, "DY", igrid)
        CALL get_fieldptr(dz, "DZ", igrid)
        CALL get_fieldptr(ddx, "DDX", igrid)
        CALL get_fieldptr(ddy, "DDY", igrid)
        CALL get_fieldptr(ddz, "DDZ", igrid)

        SELECT CASE (iface)
        CASE (3)
            ! Right
            j2 = 2
            j3 = 3
            j4 = 4
            jstag2 = 2
        CASE (4)
            ! Left
            j2 = jj - 1
            j3 = jj - 2
            j4 = jj - 3
            jstag2 = jj - 2
        CASE DEFAULT
            CALL errr(__FILE__, __LINE__)
        END SELECT

        j = MIN(j3, j4)
        IF (PRESENT(f2)) THEN
            CALL f2%get_ptr(bp, igrid)

            DO i = 3, ii-2, 2
                DO k = 3, kk-2, 2
                    CALL pressureftocone(k, j, i, kk, jj, ii, p, bp, &
                        ddx, ddy, ddz, pcnew, bpc)
                    IF (bpc < 0.5) pcnew = pbuffer(k, i)

                    sb11 = bp(k, j2, i)*bp(k, j3, i)
                    sb12 = bp(k, j2, i+1)*bp(k, j3, i+1)
                    sb13 = bp(k+1, j2, i)*bp(k+1, j3, i)
                    sb14 = bp(k+1, j2, i+1)*bp(k+1, j3, i+1)

                    fak = (sb11*ddx(i)*ddz(k) + sb12*ddx(i+1)*ddz(k) &
                        + sb13*ddx(i)*ddz(k+1) + sb14*ddx(i+1)*ddz(k+1)) &
                        /((ddx(i)+ddx(i+1))*(ddz(k)+ddz(k+1)))
                    IF (fak < 0.1) fak = 1.0
                    fak = 1.0/fak

                    DO l = 0, 1
                        DO n = 0, 1
                            p(k+n, j2, i+l) = p(k+n, j3, i+l) &
                                + fak*dy(jstag2)/(ddy(j3)+ddy(j2)) &
                                *(pbuffer(k, i) - pcnew)
                        END DO
                    END DO
                END DO
            END DO
        ELSE
            DO i = 3, ii-2, 2
                DO k = 3, kk-2, 2
                    CALL pressureftocone(k, j, i, kk, jj, ii, p, &
                        ddx, ddy, ddz, pcnew, bpc)
                    DO l = 0, 1
                        DO n = 0, 1
                            p(k+n, j2, i+l) = p(k+n, j3, i+l) &
                                + dy(jstag2)/(ddy(j3)+ddy(j2)) &
                                *(pbuffer(k, i) - pcnew)
                        END DO
                    END DO
                END DO
            END DO
        END IF
    END SUBROUTINE bright


    SUBROUTINE bbottom(igrid, iface, ibocd, ctyp, f1, f2, f3, f4, timeph)
        ! Subroutine arguments
        INTEGER(intk), INTENT(in) :: igrid, iface, ibocd
        CHARACTER(len=*), INTENT(in) :: ctyp
        TYPE(field_t), INTENT(inout) :: f1
        TYPE(field_t), INTENT(inout), OPTIONAL :: f2, f3, f4
        REAL(realk), INTENT(in), OPTIONAL :: timeph

        ! Local variables
        INTEGER(intk) :: kk, jj, ii
        INTEGER(intk) :: k, j, i, k2, k3, k4, kstag2
        INTEGER(intk) :: l, m
        REAL(realk) :: pcnew, bpc, fak
        REAL(realk) :: sb11, sb12, sb13, sb14
        REAL(realk), POINTER, CONTIGUOUS :: p(:, :, :)
        REAL(realk), POINTER, CONTIGUOUS :: pbuffer(:, :)
        REAL(realk), POINTER, CONTIGUOUS :: bp(:, :, :)
        REAL(realk), POINTER, CONTIGUOUS :: dx(:), dy(:), dz(:)
        REAL(realk), POINTER, CONTIGUOUS :: ddx(:), ddy(:), ddz(:)

        ! Only works on PAR boundaries, should do nothing otherwise
        IF (ctyp /= 'PAR') RETURN

        CALL f1%get_ptr(p, igrid)
        CALL f1%get_buffer(pbuffer, igrid, iface)
        CALL get_mgdims(kk, jj, ii, igrid)

        CALL get_fieldptr(dx, "DX", igrid)
        CALL get_fieldptr(dy, "DY", igrid)
        CALL get_fieldptr(dz, "DZ", igrid)
        CALL get_fieldptr(ddx, "DDX", igrid)
        CALL get_fieldptr(ddy, "DDY", igrid)
        CALL get_fieldptr(ddz, "DDZ", igrid)

        SELECT CASE (iface)
        CASE (5)
            ! Bottom
            k2 = 2
            k3 = 3
            k4 = 4
            kstag2 = 2
        CASE (6)
            ! Top
            k2 = kk - 1
            k3 = kk - 2
            k4 = kk - 3
            kstag2 = kk - 2
        CASE DEFAULT
            CALL errr(__FILE__, __LINE__)
        END SELECT

        k = MIN(k3, k4)
        IF (PRESENT(f2)) THEN
            CALL f2%get_ptr(bp, igrid)

            DO i = 3, ii-2, 2
                DO j = 3, jj-2, 2
                    CALL pressureftocone(k, j, i, kk, jj, ii, p, bp, &
                        ddx, ddy, ddz, pcnew, bpc)
                    IF (bpc < 0.5) pcnew = pbuffer(j, i)

                    sb11 = bp(k2, j, i)*bp(k3, j, i)
                    sb12 = bp(k2, j, i+1)*bp(k3, j, i+1)
                    sb13 = bp(k2, j+1, i)*bp(k3, j+1, i)
                    sb14 = bp(k2, j+1, i+1)*bp(k3, j+1, i+1)

                    fak = (sb11*ddx(i)*ddy(j) + sb12*ddx(i+1)*ddy(j) &
                        + sb13*ddx(i)*ddy(j+1) + sb14*ddx(i+1)*ddy(j+1)) &
                        /((ddx(i)+ddx(i+1))*(ddy(j)+ddy(j+1)))
                    IF (fak < 0.1) fak = 1.0
                    fak = 1.0/fak

                    DO l = 0, 1
                        DO m = 0, 1
                            p(k2, j+m, i+l) = p(k3, j+m, i+l) &
                                + fak*dz(kstag2)/(ddz(k3)+ddz(k2)) &
                                *(pbuffer(j, i) - pcnew)
                        END DO
                    END DO
                END DO
            END DO
        ELSE
            DO i = 3, ii-2, 2
                DO j = 3, jj-2, 2
                    CALL pressureftocone(k, j, i, kk, jj, ii, p, &
                        ddx, ddy, ddz, pcnew, bpc)
                    DO l = 0, 1
                        DO m = 0, 1
                            p(k2, j+m, i+l) = p(k3, j+m, i+l) &
                                + dz(kstag2)/(ddz(k3)+ddz(k2)) &
                                *(pbuffer(j, i) - pcnew)
                        END DO
                    END DO
                END DO
            END DO
        END IF
    END SUBROUTINE bbottom


    PURE SUBROUTINE pressureftocone_A(k, j, i, kk, jj, ii, p, ddx, ddy, &
            ddz, pc, bpc)
        ! Subroutine arguments
        INTEGER(intk), INTENT(in) :: k, j, i
        INTEGER(intk), INTENT(in) :: kk, jj, ii
        REAL(realk), INTENT(in) :: p(kk, jj, ii)
        REAL(realk), INTENT(in) :: ddx(ii), ddy(jj), ddz(kk)
        REAL(realk), INTENT(out) :: pc, bpc

        ! Local variables
        INTEGER(intk) :: l, m, n
        REAL(realk) :: vol, sump, sumvol

        bpc = 1.0

        sump = 0.0
        sumvol = 0.0
        DO l = 0, 1
            DO m = 0, 1
                DO n = 0, 1
                    vol = ddz(k+n)*ddy(j+m)*ddx(i+l)
                    sump = sump + p(k+n, j+m, i+l)*vol
                    sumvol = sumvol + vol
                END DO
            END DO
        END DO
        pc = sump/sumvol
    END SUBROUTINE pressureftocone_A


    PURE SUBROUTINE pressureftocone_B(k, j, i, kk, jj, ii, p, bp, ddx, ddy, &
            ddz, pc, bpc)
        ! Subroutine arguments
        INTEGER(intk), INTENT(in) :: k, j, i
        INTEGER(intk), INTENT(in) :: kk, jj, ii
        REAL(realk), INTENT(in) :: p(kk, jj, ii)
        REAL(realk), INTENT(in) :: bp(kk, jj, ii)
        REAL(realk), INTENT(in) :: ddx(ii), ddy(jj), ddz(kk)
        REAL(realk), INTENT(out) :: pc, bpc

        ! Local variables
        INTEGER(intk) :: l, m, n
        REAL(realk) :: vol, sumbp, sump, sumvol

        sumbp = 0.0
        DO l = 0, 1
            DO m = 0, 1
                DO n = 0, 1
                    sumbp = sumbp + bp(k+n, j+m, i+l)
                END DO
            END DO
        END DO
        bpc = MIN(sumbp, 1.0_realk)

        IF (bpc < 0.5) THEN
            pc = 0.0
        ELSE
            sump = 0.0
            sumvol = 0.0
            DO l = 0, 1
                DO m = 0, 1
                    DO n = 0, 1
                        vol = bp(k+n, j+m, i+l)*ddz(k+n)*ddy(j+m)*ddx(i+l)
                        sump = sump + p(k+n, j+m, i+l)*vol
                        sumvol = sumvol + vol
                    END DO
                END DO
            END DO
            pc = sump/sumvol
        END IF
    END SUBROUTINE pressureftocone_B


    SUBROUTINE maxabscal(maxabs, maxabslevel, phi)
        ! Subroutine arguments
        REAL(realk), INTENT(out) :: maxabs
        REAL(realk), INTENT(inout), ALLOCATABLE :: maxabslevel(:)
        TYPE(field_t), INTENT(in) :: phi

        ! Local variables
        REAL(realk) :: maxabsgrid
        INTEGER(intk) :: i, igrid, ilevel, ip3
        INTEGER(intk) :: kk, jj, ii

        maxabs = 0.0
        maxabslevel = 0.0

        DO i = 1, nmygrids
            igrid = mygrids(i)
            ilevel = level(igrid)
            CALL get_mgdims(kk, jj, ii, igrid)
            CALL get_ip3(ip3, igrid)
            CALL maxabscal_grid(kk, jj, ii, maxabsgrid, phi%arr(ip3))

            maxabs = MAX(ABS(maxabsgrid*(2.0**(maxlevel-ilevel))), maxabs)
            maxabslevel(ilevel) = MAX(maxabslevel(ilevel), maxabsgrid)
        END DO
    END SUBROUTINE maxabscal


    PURE SUBROUTINE maxabscal_grid(kk, jj, ii, maxabs, phi)
        ! Subroutine arguments
        INTEGER(intk), INTENT(in) :: kk, jj, ii
        REAL(realk), INTENT(out) :: maxabs
        REAL(realk), INTENT(in) :: phi(kk, jj, ii)

        ! Local variables
        INTEGER(intk) :: k, j, i

        maxabs = 0.0
        DO i = 3, ii-2
            DO j = 3, jj-2
                DO k = 3, kk-2
                    maxabs = MAX(ABS(phi(k, j, i)), maxabs)
                END DO
            END DO
        END DO
    END SUBROUTINE maxabscal_grid


    SUBROUTINE rescal(rhs, res)
        ! Subroutine arguments
        TYPE(field_t), INTENT(inout) :: rhs
        TYPE(field_t), INTENT(in) :: res

        ! Local variables
        INTEGER(intk) :: i, igrid, ip3
        INTEGER(intk) :: kk, jj, ii

        DO i = 1, nmygrids
            igrid = mygrids(i)
            CALL get_mgdims(kk, jj, ii, igrid)
            CALL get_ip3(ip3, igrid)
            CALL rescal_grid(kk, jj, ii, rhs%arr(ip3), res%arr(ip3))
        END DO
    END SUBROUTINE rescal


    PURE SUBROUTINE rescal_grid(kk, jj, ii, rhs, res)
        ! Subroutine arguments
        INTEGER(intk), INTENT(in) :: kk, jj, ii
        REAL(realk), INTENT(inout) :: rhs(kk, jj, ii)
        REAL(realk), INTENT(in) :: res(kk, jj, ii)

        ! Local variables
        INTEGER(intk) :: k, j, i

        ! TODO: Check if indices can be extended
        DO i = 3, ii-2
            DO j = 3, jj-2
                DO k = 3, kk-2
                    rhs(k, j, i) = rhs(k, j, i) + res(k, j, i)
                END DO
            END DO
        END DO
    END SUBROUTINE rescal_grid


    SUBROUTINE mgpcorr(u, v, w, p, dp, rfak, bp_f)
        ! Subroutine arguments
        TYPE(field_t), INTENT(inout) :: u, v, w, p
        TYPE(field_t), INTENT(in) :: dp
        REAL(realk), INTENT(in) :: rfak
        TYPE(field_t), INTENT(in), OPTIONAL :: bp_f

        ! Local variables
        INTEGER(intk) :: i, igrid, ip3
        INTEGER(intk) :: kk, jj, ii

        TYPE(field_t), POINTER :: rdx_f
        TYPE(field_t), POINTER :: rdy_f
        TYPE(field_t), POINTER :: rdz_f

        REAL(realk), POINTER, CONTIGUOUS :: rdx(:), rdy(:), rdz(:), bp(:, :, :)

        NULLIFY(bp)

        CALL get_field(rdx_f, "RDX")
        CALL get_field(rdy_f, "RDY")
        CALL get_field(rdz_f, "RDZ")

        DO i = 1, nmygrids
            igrid = mygrids(i)
            CALL get_mgdims(kk, jj, ii, igrid)
            CALL get_ip3(ip3, igrid)

            CALL rdx_f%get_ptr(rdx, igrid)
            CALL rdy_f%get_ptr(rdy, igrid)
            CALL rdz_f%get_ptr(rdz, igrid)
            IF (PRESENT(bp_f)) THEN
                CALL bp_f%get_ptr(bp, igrid)
            END IF

            CALL mgpcorr_grid(kk, jj, ii, u%arr(ip3), v%arr(ip3), w%arr(ip3), &
                p%arr(ip3), dp%arr(ip3), rdx, rdy, rdz, rfak, bp)
        END DO
    END SUBROUTINE mgpcorr


    PURE SUBROUTINE mgpcorr_grid(kk, jj, ii, u, v, w, p, dp, rdx, rdy, rdz, &
            rfak, bp)
        ! Subroutine arguments
        INTEGER(intk), INTENT(in) :: kk, jj, ii
        REAL(realk), INTENT(inout) :: u(kk, jj, ii)
        REAL(realk), INTENT(inout) :: v(kk, jj, ii)
        REAL(realk), INTENT(inout) :: w(kk, jj, ii)
        REAL(realk), INTENT(inout) :: p(kk, jj, ii)
        REAL(realk), INTENT(in) :: dp(kk, jj, ii)
        REAL(realk), INTENT(in) :: rdx(ii)
        REAL(realk), INTENT(in) :: rdy(jj)
        REAL(realk), INTENT(in) :: rdz(kk)
        REAL(realk), INTENT(in) :: rfak
        REAL(realk), INTENT(in), OPTIONAL :: bp(kk, jj, ii)

        ! Local variables
        INTEGER(intk) :: k, j, i

        IF (PRESENT(bp)) THEN
            DO i = 2, ii-1
                DO j = 2, jj-1
                    DO k = 2, kk-1
                        p(k, j, i) = p(k, j, i) + dp(k, j, i)*bp(k, j, i)
                    END DO
                END DO
            END DO

            DO i = 2, ii-2
                DO j = 3, jj-2
                    DO k = 3, kk-2
                        u(k, j, i) = u(k, j, i) &
                            + (dp(k, j, i) - dp(k, j, i+1)) &
                            *bp(k, j, i)*bp(k, j, i+1)*rdx(i)*rfak
                    END DO
                END DO
            END DO

            DO i = 3, ii-2
                DO j = 2, jj - 2
                    DO k = 3, kk-2
                        v(k, j, i) = v(k, j, i) &
                            + (dp(k, j, i) - dp(k, j+1, i)) &
                            *bp(k, j, i)*bp(k, j+1, i)*rdy(j)*rfak
                    END DO
                END DO
            END DO

            DO i = 3, ii-2
                DO j = 3, jj-2
                    DO k = 2, kk-2
                        w(k, j, i) = w(k, j, i) &
                            + (dp(k, j, i) - dp(k+1, j, i)) &
                            *bp(k, j, i)*bp(k+1, j, i)*rdz(k)*rfak
                    END DO
                END DO
            END DO
        ELSE
            DO i = 2, ii-1
                DO j = 2, jj-1
                    DO k = 2, kk-1
                        p(k, j, i) = p(k, j, i) + dp(k, j, i)
                    END DO
                END DO
            END DO

            DO i = 2, ii-2
                DO j = 3, jj-2
                    DO k = 3, kk-2
                        u(k, j, i) = u(k, j, i) &
                            + (dp(k, j, i) - dp(k, j, i+1))*rdx(i)*rfak
                    END DO
                END DO
            END DO

            DO i = 3, ii-2
                DO j = 2, jj - 2
                    DO k = 3, kk-2
                        v(k, j, i) = v(k, j, i) &
                            + (dp(k, j, i) - dp(k, j+1, i))*rdy(j)*rfak
                    END DO
                END DO
            END DO

            DO i = 3, ii-2
                DO j = 3, jj-2
                    DO k = 2, kk-2
                        w(k, j, i) = w(k, j, i) &
                            + (dp(k, j, i) - dp(k+1, j, i))*rdz(k)*rfak
                    END DO
                END DO
            END DO
        END IF
    END SUBROUTINE mgpcorr_grid
END MODULE pressuresolver_mod
