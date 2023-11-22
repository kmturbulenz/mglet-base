MODULE pressuresolver_mod
    USE bound_flow_mod
    USE core_mod
    USE flowcore_mod
    USE ib_mod
    USE itinfo_mod, ONLY: itinfo_sample
    USE plog_mod

    IMPLICIT NONE (type, external)
    PRIVATE

    ! Type of pressure solver
    !   0 : Standard SIP
    !   1 : SIP on coarsest level, then SOR on subsequent levels
    INTEGER(intk) :: ityp

    ! Number of inner iterations/sweeps in pressure solver
    INTEGER(intk) :: ninner

    ! Maximum number of outer iterations allowed
    INTEGER(intk) :: nouter

    ! Minimum number of outer iterations to always run independent of
    ! divergence
    INTEGER(intk) :: nouter_min

    ! Relaxation factor for SOR
    REAL(realk) :: omg

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
    INTERFACE sipiter1
        MODULE PROCEDURE :: sipiter1_A, sipiter1_B
    END INTERFACE sipiter1

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
        CASE DEFAULT
            WRITE(*, '("Invalid pressure solver type: ", A)') type
            CALL errr(__FILE__, __LINE__)
        END SELECT

        ! Always intialize SIP - it's always used at the coarsest level
        CALL init_sip()

        ! Initialize SOR if enabled
        IF (ityp == 1) THEN
            CALL psolveconf%get_value("/omega", omg, default_value=1.1)
            IF (omg <= 0.0 .OR. omg >= 2.0) THEN
                WRITE(*, '("Invalid omega: ", F15.7)') omg
                CALL errr(__FILE__, __LINE__)
            END IF

            CALL init_sor()
        END IF
    END SUBROUTINE init_pressuresolver


    SUBROUTINE finish_pressuresolver()
        CALL finish_plog()
    END SUBROUTINE finish_pressuresolver


    SUBROUTINE init_sip()
        ! Subroutine arguments
        ! none...

        ! Local variables
        INTEGER(intk) :: igr, igrid
        INTEGER(intk) :: k, j, i
        INTEGER(intk) :: kk, jj, ii
        REAL(realk), PARAMETER :: alfa = 0.92
        REAL(realk) :: p1, p2, p3
        REAL(realk), POINTER, CONTIGUOUS :: bp(:, :, :)
        REAL(realk), POINTER, CONTIGUOUS :: lb(:, :, :), lw(:, :, :), &
            ls(:, :, :), ue(:, :, :), un(:, :, :), ut(:, :, :), lpr(:, :, :)
        REAL(realk), POINTER, CONTIGUOUS :: ae(:),  aw(:),  an(:), as(:), &
            at(:), ab(:), ap(:, :, :)

        CALL set_field("SIPLW")
        CALL set_field("SIPLS")
        CALL set_field("SIPLB")
        CALL set_field("SIPUE")
        CALL set_field("SIPUN")
        CALL set_field("SIPUT")
        CALL set_field("SIPLPR")

        DO igr = 1, nmygrids
            igrid = mygrids(igr)

            CALL get_mgdims(kk, jj, ii, igrid)

            CALL get_fieldptr(lw, "SIPLW", igrid)
            CALL get_fieldptr(ls, "SIPLS", igrid)
            CALL get_fieldptr(lb, "SIPLB", igrid)
            CALL get_fieldptr(ue, "SIPUE", igrid)
            CALL get_fieldptr(un, "SIPUN", igrid)
            CALL get_fieldptr(ut, "SIPUT", igrid)
            CALL get_fieldptr(lpr, "SIPLPR", igrid)

            CALL get_fieldptr(aw, "GSAW", igrid)
            CALL get_fieldptr(ae, "GSAE", igrid)
            CALL get_fieldptr(as, "GSAS", igrid)
            CALL get_fieldptr(an, "GSAN", igrid)
            CALL get_fieldptr(ab, "GSAB", igrid)
            CALL get_fieldptr(at, "GSAT", igrid)

            CALL get_fieldptr(ap, "GSAP", igrid)

            CALL get_fieldptr(bp, "BP", igrid)

            DO i = 3, ii-2
                DO j = 3, jj-2
                    DO k = 3, kk-2
                        lw(k, j, i) = aw(i)*bu(k, j, i-1) &
                            /(1.0 + alfa*(un(k, j, i-1) + ut(k, j, i-1)))

                        ls(k, j, i) = as(j)*bv(k, j-1, i) &
                            /(1.0 + alfa*(ue(k, j-1, i) + ut(k, j-1, i)))

                        lb(k,j,i) = ab(k)*bw(k-1, j, i) &
                            /(1.0 + alfa*(un(k-1, j, i) + ue(k-1, j, i)))

                        p1 = alfa*(lb(k, j, i)*ue(k-1, j, i) &
                            + ls(k, j, i)*ue(k, j-1, i))
                        p2 = alfa*(lb(k, j, i)*un(k-1, j, i) &
                            + lw(k, j, i)*un(k, j, i-1))
                        p3 = alfa*(lw(k, j, i)*ut(k, j, i-1) &
                            + ls(k, j, i)*ut(k, j-1, i))

                        lpr(k, j, i) = 1.0/(ap(k, j, i) + p1 + p2 + p3 &
                            - lb(k, j, i)*ut(k-1, j, i) &
                            - lw(k, j, i)*ue(k, j, i-1) &
                            - ls(k, j, i)*un(k, j-1, i) &
                            + 1.0e-20)

                        ue(k, j, i) = (ae(i)*bu(k, j, i) - p1)*lpr(k, j, i)
                        un(k, j, i) = (an(j)*bv(k, j, i) - p2)*lpr(k, j, i)
                        ut(k, j, i) = (at(k)*bw(k, j, i) - p3)*lpr(k, j, i)
                    END DO
                END DO
             END DO

            DO i = 3, ii-2
                DO j = 3, jj-2
                    DO k = 3, kk-2
                        lw(k, j, i) = lw(k, j, i)*lpr(k, j, i)
                    END DO
                    DO k = 3, kk-2
                        ls(k, j, i) = ls(k, j, i)*lpr(k, j, i)
                    END DO
                    DO k = 3, kk-2
                        lb(k, j, i) = lb(k, j, i)*lpr(k, j, i)
                    END DO
                END DO
            END DO
        END DO

    CONTAINS
        ! These routines are not protected against out-of-bounds and not valid
        ! for the edges ii, jj, kk - but in the scope above that is OK
        PURE REAL(realk) FUNCTION bu(k, j, i)
            INTEGER(intk), INTENT(in) :: k, j, i
            bu = bp(k, j, i)*bp(k, j, i+1)
        END FUNCTION bu

        PURE REAL(realk) FUNCTION bv(k, j, i)
            INTEGER(intk), INTENT(in) :: k, j, i
            bv = bp(k, j, i)*bp(k, j+1, i)
        END FUNCTION bv

        PURE REAL(realk) FUNCTION bw(k, j, i)
            INTEGER(intk), INTENT(in) :: k, j, i
            bw = bp(k, j, i)*bp(k+1, j, i)
        END FUNCTION bw
    END SUBROUTINE init_sip


    SUBROUTINE init_sor()
        ! Subroutine arguments
        ! none...

        ! Local variables
        INTEGER(intk) :: igr, igrid
        INTEGER(intk) :: k, j, i
        INTEGER(intk) :: kk, jj, ii
        REAL(realk), POINTER, CONTIGUOUS :: ap(:, :, :), rap(:, :, :)

        CALL set_field("SOR_RAP")

        DO igr = 1, nmygrids
            igrid = mygrids(igr)

            CALL get_mgdims(kk, jj, ii, igrid)
            CALL get_fieldptr(rap, "SOR_RAP", igrid)
            CALL get_fieldptr(ap, "GSAP", igrid)
            DO i = 3, ii-2
                DO j = 3, jj-2
                    DO k = 3, kk-2
                        rap(k, j, i) = divide0(1.0_realk, ap(k, j, i))
                    END DO
                END DO
             END DO
        END DO
    END SUBROUTINE init_sor


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
        TYPE(field_t) :: dp, hilf, rhs, res
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
        CALL dp%init("DP")
        CALL hilf%init("HILF")
        CALL rhs%init("RHS")
        CALL res%init("RES")

        CALL dp%init_buffers()
        CALL hilf%init_buffers()

        ! laplace(dp) = prefak * div(u) is the underlying equation
        prefak = rho/dt
        CALL ib%divcal(rhs, u, v, w, prefak)

        DO ilevel = maxlevel, minlevel, -1
            CALL ftoc(ilevel, rhs%arr, rhs%arr, 'R')
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
            DO ilevel = minlevel, maxlevel
                CALL ctof(ilevel, hilf%arr, hilf%arr)
                CALL parent(ilevel, s1=hilf)
                CALL mgpoisit(ilevel, hilf, rhs, res, bp)
            END DO

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
                CALL ftoc(ilevel, hilf%arr, hilf%arr, 'P')
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

            DO ilevel = maxlevel, minlevel, -1
                CALL ftoc(ilevel, rhs%arr, rhs%arr, 'R')
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
            CALL ftoc(ilevel, u%arr, u%arr, 'U')
            CALL ftoc(ilevel, v%arr, v%arr, 'V')
            CALL ftoc(ilevel, w%arr, w%arr, 'W')
            CALL ftoc(ilevel, p%arr, p%arr, 'P')
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

        CALL res%finish()
        CALL rhs%finish()
        CALL hilf%finish()
        CALL dp%finish()

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

            CALL sipiter1(kk, jj, ii, rhs_p, res_p, lw, ls, lb, lpr)
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

            CALL sipiter2(kk, jj, ii, dp_p, res_p, ue, un, ut)
        END DO
    END SUBROUTINE sip


    SUBROUTINE sor(ilevel, dp, rhs, gsaw, gsae, gsas, gsan, gsab, gsat, &
            gsrap, bp)
        ! Subroutine arguments
        INTEGER(intk), INTENT(in) :: ilevel
        TYPE(field_t), INTENT(inout) :: dp
        TYPE(field_t), INTENT(in) :: rhs
        TYPE(field_t), INTENT(in) :: gsaw
        TYPE(field_t), INTENT(in) :: gsae
        TYPE(field_t), INTENT(in) :: gsas
        TYPE(field_t), INTENT(in) :: gsan
        TYPE(field_t), INTENT(in) :: gsab
        TYPE(field_t), INTENT(in) :: gsat
        TYPE(field_t), INTENT(in) :: gsrap
        TYPE(field_t), INTENT(in), OPTIONAL :: bp

        ! Local variables
        INTEGER(intk) :: i, igrid
        INTEGER(intk) :: kk, jj, ii
        REAL(realk), POINTER, CONTIGUOUS :: aw(:), ae(:), as(:), an(:), &
            ab(:), at(:), rap(:, :, :)
        REAL(realk), POINTER, CONTIGUOUS :: dp_p(:, :, :), &
            rhs_p(:, :, :), bp_p(:, :, :)

        ! Ensure this does not point to anything
        NULLIFY(bp_p)

        DO i = 1, nmygridslvl(ilevel)
            igrid = mygridslvl(i, ilevel)

            CALL dp%get_ptr(dp_p, igrid)
            CALL rhs%get_ptr(rhs_p, igrid)

            CALL gsaw%get_ptr(aw, igrid)
            CALL gsae%get_ptr(ae, igrid)
            CALL gsas%get_ptr(as, igrid)
            CALL gsan%get_ptr(an, igrid)
            CALL gsab%get_ptr(ab, igrid)
            CALL gsat%get_ptr(at, igrid)
            CALL gsrap%get_ptr(rap, igrid)

            IF (PRESENT(bp)) CALL bp%get_ptr(bp_p, igrid)

            CALL get_mgdims(kk, jj, ii, igrid)
            CALL relax(kk, jj, ii, dp_p, rhs_p, aw, ae, as, an, &
                ab, at, rap, bp_p)
        END DO
    END SUBROUTINE sor


    SUBROUTINE laplacephi(res, phi, bp)
        ! Subroutine arguments
        TYPE(field_t), INTENT(inout) :: res
        TYPE(field_t), INTENT(in) :: phi
        TYPE(field_t), INTENT(in), OPTIONAL :: bp

        ! Local variables
        INTEGER(intk) :: ilevel

        DO ilevel = minlevel, maxlevel
            CALL laplacephi_level(ilevel, res, phi, bp)
        END DO
    END SUBROUTINE laplacephi


    SUBROUTINE laplacephi_level(ilevel, res_f, phi_f, bp_f)
        ! Subroutine arguments
        INTEGER(intk), INTENT(in):: ilevel
        TYPE(field_t), INTENT(inout) :: res_f
        TYPE(field_t), INTENT(in) :: phi_f
        TYPE(field_t), INTENT(in), OPTIONAL :: bp_f

        ! Local variables
        INTEGER(intk) :: i, igrid
        INTEGER(intk) :: kk, jj, ii

        TYPE(field_t), POINTER :: gsaw, gsae, gsas, gsan, gsab, gsat, gsap
        REAL(realk), POINTER, CONTIGUOUS :: phi(:, :, :), res(:, :, :)
        REAL(realk), POINTER, CONTIGUOUS :: aw(:), ae(:), as(:), an(:), &
            ab(:), at(:), ap(:, :, :)
        REAL(realk), POINTER, CONTIGUOUS :: bp(:, :, :)

        NULLIFY(bp)

        CALL get_field(gsaw, "GSAW")
        CALL get_field(gsae, "GSAE")
        CALL get_field(gsas, "GSAS")
        CALL get_field(gsan, "GSAN")
        CALL get_field(gsab, "GSAB")
        CALL get_field(gsat, "GSAT")
        CALL get_field(gsap, "GSAP")

        DO i = 1, nmygridslvl(ilevel)
            igrid = mygridslvl(i, ilevel)
            CALL get_mgdims(kk, jj, ii, igrid)

            CALL phi_f%get_ptr(phi, igrid)
            CALL res_f%get_ptr(res, igrid)

            CALL gsaw%get_ptr(aw, igrid)
            CALL gsae%get_ptr(ae, igrid)
            CALL gsas%get_ptr(as, igrid)
            CALL gsan%get_ptr(an, igrid)
            CALL gsab%get_ptr(ab, igrid)
            CALL gsat%get_ptr(at, igrid)
            CALL gsap%get_ptr(ap, igrid)

            IF (PRESENT(bp_f)) THEN
                CALL bp_f%get_ptr(bp, igrid)
            END IF

            CALL laplacephi_grid(kk, jj, ii, res, phi, aw, ae, an, as, &
                at, ab, ap, bp)
        END DO
    END SUBROUTINE laplacephi_level


    PURE SUBROUTINE laplacephi_grid(kk, jj, ii, res, phi, aw, ae, an, as, &
            at, ab, ap, bp)
        ! Subroutine arguments
        INTEGER(intk), INTENT(in) :: kk, jj, ii
        REAL(realk), INTENT(out) :: res(kk, jj, ii)
        REAL(realk), INTENT(in) :: phi(kk, jj, ii)
        REAL(realk), INTENT(in) :: aw(ii), ae(ii), an(jj), as(jj), &
            at(kk), ab(kk)
        REAL(realk), INTENT(in) :: ap(kk, jj, ii)
        REAL(realk), INTENT(in), OPTIONAL :: bp(kk, jj, ii)

        ! Local variables
        INTEGER :: k, j, i

        ! To properly define an INTENT(out) since the loop does not iterate
        ! over all indices of res
        res = 0.0

        IF (PRESENT(bp)) THEN
            DO i = 3, ii-2
                DO j = 3, jj-2
                    !$omp simd
                    DO k = 3, kk-2
                        res(k, j, i) = &
                            - aw(i)*phi(k, j, i-1)*bp(k, j, i-1)*bp(k, j, i) &
                            - ae(i)*phi(k, j, i+1)*bp(k, j, i)*bp(k, j, i+1) &
                            - as(j)*phi(k, j-1, i)*bp(k, j-1, i)*bp(k, j, i) &
                            - an(j)*phi(k, j+1, i)*bp(k, j, i)*bp(k, j+1, i) &
                            - ab(k)*phi(k-1, j, i)*bp(k-1, j, i)*bp(k, j, i) &
                            - at(k)*phi(k+1, j, i)*bp(k, j, i)*bp(k+1, j, i) &
                            - ap(k, j, i)*phi(k, j, i)
                    END DO
                END DO
            END DO
        ELSE
            DO i = 3, ii-2
                DO j = 3, jj-2
                    !$omp simd
                    DO k = 3, kk-2
                        res(k, j, i) = &
                            - aw(i) * phi(k, j, i-1) &
                            - ae(i) * phi(k, j, i+1) &
                            - as(j) * phi(k, j-1, i) &
                            - an(j) * phi(k, j+1, i) &
                            - ab(k) * phi(k-1, j, i) &
                            - at(k) * phi(k+1, j, i) &
                            - ap(k, j, i) * phi(k, j, i)
                    END DO
                END DO
            END DO
        END IF
    END SUBROUTINE laplacephi_grid


    PURE SUBROUTINE sipiter0(kk, jj, ii, rhs, res, lpr)
        ! Subroutine arguments
        INTEGER(intk), INTENT(in) :: kk, jj, ii
        REAL(realk), INTENT(in) :: rhs(kk, jj, ii)
        REAL(realk), INTENT(inout) :: res(kk, jj, ii)
        REAL(realk), INTENT(in) :: lpr(kk, jj, ii)

        ! Local variables
        INTEGER(intk) :: k, j, i

        DO i = 3, ii-2
            DO j = 3, jj-2
                DO k = 3, kk-2
                    res(k, j, i) = (rhs(k, j, i) + res(k, j, i))*lpr(k, j, i)
                END DO
            END DO
        END DO
    END SUBROUTINE sipiter0


    PURE SUBROUTINE sipiter1_A(kk, jj, ii, rhs, res, lw, ls, lb)
        ! Subroutine arguments
        INTEGER(intk), INTENT(in) :: kk, jj, ii
        REAL(realk), INTENT(inout) :: res(kk, jj, ii)
        REAL(realk), INTENT(in) :: rhs(kk, jj, ii)
        REAL(realk), INTENT(in) :: lw(kk, jj, ii), ls(kk, jj, ii), &
            lb(kk, jj, ii)

        ! Local variables
        INTEGER(intk) :: k, j, i

        DO i = 3, ii-2
            DO j = 3, jj-2
                DO k = 3, kk-2
                    res(k, j, i) = res(k, j, i) - lw(k, j, i)*res(k, j, i-1) &
                        - ls(k, j, i)*res(k, j-1, i)
                END DO
                DO k = 3, kk-2
                    res(k, j, i) = res(k, j, i) - lb(k, j, i)*res(k-1, j, i)
                END DO
            END DO
        END DO
    END SUBROUTINE sipiter1_A


    PURE SUBROUTINE sipiter1_B(kk, jj, ii, rhs, res, lw, ls, lb, lpr)
        ! Subroutine arguments
        INTEGER(intk), INTENT(in) :: kk, jj, ii
        REAL(realk), INTENT(inout) :: res(kk, jj, ii)
        REAL(realk), INTENT(in) :: rhs(kk, jj, ii)
        REAL(realk), INTENT(in) :: lw(kk, jj, ii), ls(kk, jj, ii), &
            lb(kk, jj, ii), lpr(kk, jj, ii)

        ! Local variables
        INTEGER(intk) :: k, j, i

        DO i = 3, ii-2
            DO j = 3, jj-2
                DO k = 3, kk-2
                    res(k, j, i) = (rhs(k, j, i) + res(k, j, i))*lpr(k, j, i)
                    res(k, j, i) = res(k, j, i) - lw(k, j, i)*res(k, j, i-1) &
                        - ls(k, j, i)*res(k, j-1, i)
                END DO
                DO k = 3, kk-2
                    res(k, j, i) = res(k, j, i) - lb(k, j, i)*res(k-1, j, i)
                END DO
            END DO
        END DO
    END SUBROUTINE sipiter1_B


    PURE SUBROUTINE sipiter2(kk, jj, ii, phi, res, ue, un, ut)
        ! Subroutine arguments
        INTEGER(intk), INTENT(in) :: kk, jj, ii
        REAL(realk), INTENT(inout) :: phi(kk, jj, ii), res(kk, jj, ii)
        REAL(realk), INTENT(in) :: ue(kk, jj, ii), un(kk, jj, ii), &
            ut(kk, jj, ii)

        ! Local variables
        INTEGER(intk) :: k, j, i

        DO i = ii-2, 3, -1
            DO j = jj-2, 3, -1
                DO k = 3, kk-2
                    res(k, j, i) = res(k, j, i) - un(k, j, i)*res(k, j+1, i) &
                        - ue(k, j, i)*res(k, j, i+1)
                END DO
                DO k = kk-2, 3, -1
                    res(k, j, i) = res(k, j, i) - ut(k, j, i)*res(k+1, j, i)
                END DO
            END DO
        END DO

        DO i = 3, ii-2
            DO j = 3, jj-2
                DO k = 3, kk-2
                    phi(k, j, i) = phi(k, j, i) + res(k, j, i)
                END DO
            END DO
        END DO
    END SUBROUTINE sipiter2


    PURE SUBROUTINE relax(kk, jj, ii, dp, rhs, gsaw, gsae, gsas, gsan, &
            gsab, gsat, gsrap, bp)

        ! Subroutine arguments
        INTEGER(intk), INTENT(in) :: kk, jj, ii
        REAL(realk), INTENT(inout) :: dp(kk, jj, ii)
        REAL(realk), INTENT(in) :: rhs(kk, jj, ii)
        REAL(realk), INTENT(in) :: gsaw(ii), gsae(ii), gsas(jj), gsan(jj), &
            gsab(kk), gsat(kk)
        REAL(realk), INTENT(in) :: gsrap(kk, jj, ii)
        REAL(realk), INTENT(in), OPTIONAL :: bp(kk, jj, ii)

        ! Local variables
        INTEGER(intk) :: i, j, k, rb
        INTEGER(intk) :: kstart, kstop
        REAL(realk) :: res
        REAL(realk) :: aw, ae, as, an, ab, at, rap

        IF (PRESENT(bp)) THEN
            DO rb = 0, 1
                DO i = 3, ii-2
                    DO j = 3, jj-2
                        IF (MOD(i + j, 2) == rb) THEN
                            kstart = 4
                            kstop = kk - 2
                        ELSE IF (MOD(i + j, 2) /= rb) THEN
                            kstart = 3
                            kstop = kk - 3
                        END IF

                        !$omp simd private(aw, ae, as, an, ab, at, rap, res)
                        DO k = kstart, kstop, 2
                            ! Variations in numerical formulation, please
                            ! keep for future reference. Should be the same
                            ! numerics, but performance vary.
                            !
                            ! aw = bu(k, j, i-1)/(dx(i-1)*dx(i))
                            ! ae = bu(k, j, i  )/(dx(i  )*dx(i))
                            ! as = bv(k, j-1, i)/(dy(j-1)*dy(j))
                            ! an = bv(k, j  , i)/(dy(j  )*dy(j))
                            ! ab = bw(k-1, j, i)/(dz(k-1)*dz(k))
                            ! at = bw(k  , j, i)/(dz(k  )*dz(k))
                            !
                            ! aw = bp(k, j, i-1)*bp(k, j, i  )/(dx(i-1)*dx(i))
                            ! ae = bp(k, j, i  )*bp(k, j, i+1)/(dx(i  )*dx(i))
                            ! as = bp(k, j-1, i)*bp(k, j  , i)/(dy(j-1)*dy(j))
                            ! an = bp(k, j  , i)*bp(k, j+1, i)/(dy(j  )*dy(j))
                            ! ab = bp(k-1, j, i)*bp(k  , j, i)/(dz(k-1)*dz(k))
                            ! at = bp(k  , j, i)*bp(k+1, j, i)/(dz(k  )*dz(k))

                            aw = bp(k, j, i-1)*bp(k, j, i  )*gsaw(i)
                            ae = bp(k, j, i  )*bp(k, j, i+1)*gsae(i)
                            as = bp(k, j-1, i)*bp(k, j  , i)*gsas(j)
                            an = bp(k, j  , i)*bp(k, j+1, i)*gsan(j)
                            ab = bp(k-1, j, i)*bp(k  , j, i)*gsab(k)
                            at = bp(k  , j, i)*bp(k+1, j, i)*gsat(k)
                            rap = gsrap(k, j, i)

                            res = (aw * dp(k, j, i-1) &
                                + ae * dp(k, j, i+1) &
                                + as * dp(k, j-1, i) &
                                + an * dp(k, j+1, i) &
                                + ab * dp(k-1, j, i) &
                                + at * dp(k+1, j, i) &
                                - rhs(k, j, i)) * rap

                            dp(k, j, i) = (1.0 - omg)*dp(k, j, i) - omg*res
                        END DO
                    END DO
                END DO
            END DO
        ELSE
            DO rb = 0, 1
                DO i = 3, ii-2
                    DO j = 3, jj-2
                        IF (MOD(i + j, 2) == rb) THEN
                            kstart = 4
                            kstop = kk - 2
                        ELSE IF (MOD(i + j, 2) /= rb) THEN
                            kstart = 3
                            kstop = kk - 3
                        END IF

                        !$omp simd private(res)
                        DO k = kstart, kstop, 2
                            res = (gsaw(i) * dp(k, j, i-1) &
                                  + gsae(i) * dp(k, j, i+1) &
                                  + gsas(j) * dp(k, j-1, i) &
                                  + gsan(j) * dp(k, j+1, i) &
                                  + gsab(k) * dp(k-1, j, i) &
                                  + gsat(k) * dp(k+1, j, i) &
                                  - rhs(k, j, i)) * gsrap(k, j, i)

                            dp(k, j, i) = (1.0 - omg)*dp(k, j, i) - omg*res
                        END DO
                    END DO
                END DO
            END DO
        END IF
    END SUBROUTINE relax


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
        REAL(realk), POINTER, CONTIGUOUS :: pbuffer(:, :, :)
        REAL(realk), POINTER, CONTIGUOUS :: bp(:, :, :)
        REAL(realk), POINTER, CONTIGUOUS :: dx(:), dy(:), dz(:)
        REAL(realk), POINTER, CONTIGUOUS :: ddx(:), ddy(:), ddz(:)

        ! Only works on PAR boundaries, should do nothing otherwise
        IF (ctyp /= 'PAR') RETURN

        CALL f1%get_ptr(p, igrid)
        CALL f1%buffers%get_buffer(pbuffer, igrid, iface)
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
                    IF (bpc < 0.5) pcnew = pbuffer(k, j, 1)

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
                                *(pbuffer(k, j, 1) - pcnew)
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
                                *(pbuffer(k, j, 1) - pcnew)
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
        REAL(realk), POINTER, CONTIGUOUS :: pbuffer(:, :, :)
        REAL(realk), POINTER, CONTIGUOUS :: bp(:, :, :)
        REAL(realk), POINTER, CONTIGUOUS :: dx(:), dy(:), dz(:)
        REAL(realk), POINTER, CONTIGUOUS :: ddx(:), ddy(:), ddz(:)

        ! Only works on PAR boundaries, should do nothing otherwise
        IF (ctyp /= 'PAR') RETURN

        CALL f1%get_ptr(p, igrid)
        CALL f1%buffers%get_buffer(pbuffer, igrid, iface)
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
                    IF (bpc < 0.5) pcnew = pbuffer(k, i, 1)

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
                                *(pbuffer(k, i, 1) - pcnew)
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
                                *(pbuffer(k, i, 1) - pcnew)
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
        REAL(realk), POINTER, CONTIGUOUS :: pbuffer(:, :, :)
        REAL(realk), POINTER, CONTIGUOUS :: bp(:, :, :)
        REAL(realk), POINTER, CONTIGUOUS :: dx(:), dy(:), dz(:)
        REAL(realk), POINTER, CONTIGUOUS :: ddx(:), ddy(:), ddz(:)

        ! Only works on PAR boundaries, should do nothing otherwise
        IF (ctyp /= 'PAR') RETURN

        CALL f1%get_ptr(p, igrid)
        CALL f1%buffers%get_buffer(pbuffer, igrid, iface)
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
                    IF (bpc < 0.5) pcnew = pbuffer(j, i, 1)

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
                                *(pbuffer(j, i, 1) - pcnew)
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
                                *(pbuffer(j, i, 1) - pcnew)
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
