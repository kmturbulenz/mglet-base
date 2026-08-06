MODULE pressuresolver_mod
    USE bound_flow_mod
    USE bound_pressure_mod
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

#ifdef _MGLET_OFFLOAD_
        IF (ityp /= 0) THEN
            WRITE(*, *) "MGLET_OFFLOAD only supports hyperplane sip."
            CALL errr(__FILE__, __LINE__)
        END IF
#endif

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

        CALL init_bound_pressure()
    END SUBROUTINE init_pressuresolver


    SUBROUTINE finish_pressuresolver()
        CALL finish_bound_pressure()
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

        CALL push_field(dp, "DP")
        CALL push_field(hilf, "HILF")
        CALL push_field(rhs, "RHS")
        CALL push_field(res, "RES")
        CALL set_field_arr(dp, 0.0_realk, device=.TRUE.)
        CALL set_field_arr(hilf, 0.0_realk, device=.TRUE.)
        CALL set_field_arr(rhs, 0.0_realk)
        CALL set_field_arr(res, 0.0_realk, device=.TRUE.)


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

        ! TODO(offload): Remove once surrounding subroutines are offloaded
        CALL map_arr_to_device(rhs, message="to:rhs%arr")

        ipc = 0
        outer: DO ipcount = 1, nouter
            ! Inner pressure iterations
            ! HINT: 'res' is passed into mgpoisit as a temporary storage!!
            CALL start_timer(322)
            DO ilevel = minlevel, maxlevel
                CALL ctof(ilevel, hilf, hilf, device=.TRUE.)
                CALL parent(ilevel, s1=hilf, device=.TRUE.)
                CALL mgpoisit(ilevel, hilf, rhs, res, bp)
            END DO
            CALL stop_timer(322)

            ! TODO(offload): Remove once surrounding subroutines offloaded
            CALL map_arr_from_device(hilf, message="from:hilf%arr")

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

            ! TODO(offload): Remove once surrounding subroutines are offloaded
            CALL map_arr_to_device(hilf, message="to:hilf%arr")

            ! --- intermediate state ---
            ! every grid level has the best solution
            ! for hilf retrieved from the locally
            ! finest grid.

            ! Connect needed due to prior ftoc, since this does not do
            ! anything on the finest level, no need to connect finest level
            ! either.
            CALL conn(layers=1, s1=hilf)

            ! res <- laplace(hilf)
            CALL laplacephi(res, hilf)
            ! rhs <- rhs + res
            CALL rescal(rhs, res)
            ! TODO(offload): Remove once surrounding subroutines are offloaded
            CALL map_arr_from_device(rhs, message="from:rhs%arr")

            DO ilevel = maxlevel, minlevel+1, -1
                CALL ftoc(ilevel, rhs, rhs, 'R')
            END DO

            ! TODO(offload): Remove once surrounding subroutines are offloaded
            CALL map_arr_to_device(rhs, message="to:rhs%arr")

            ! Max of RHS scaled according to levels
            CALL maxabscal(maxrhs, maxrhslvl, rhs)

            ! dp = dp + hilf
            CALL accumulate_pcorr(dp, hilf)
            CALL set_field_arr(hilf, 0.0_realk, device=.TRUE.)
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
            CALL parent(ilevel, s1=dp, device=.TRUE.)
            CALL bound_pressure(ilevel, dp, bp)
        END DO

        ! TODO(offload): Remove once surrounding subroutines are offloaded
        CALL map_arr_to_device(u, v, w, p, message="to:u|v|w|p%arr")

        ! Pressure correction: P = P + dtrk/rho*DP
        ! Velocity fields are modified and become solenoidal based on DP
        CALL mgpcorr(u, v, w, p, dp, dt/rho)

        ! TODO(offload): Remove once surrounding subroutines are offloaded
        CALL map_arr_from_device(u, v, w, p, message="from:u|v|w|p%arr")

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
            CALL bound_pressure(ilevel, dp, bp)

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

            CALL conn(ilevel, 1, s1=dp)
        END DO

        CALL bound_pressure(ilevel, dp, bp)

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
        ! none...

        CALL laplacephi_level(ilevel, res, dp)

        IF (ityp == 2) THEN
            CALL sipiter1_classic_level(ilevel, res, rhs, siplw, sipls, siplb, &
                siplpr)
        ELSE
            CALL sipiter1_hyperplane_level(ilevel, res, rhs, siplw, sipls, &
                siplb, siplpr)
        END IF

        IF (iloop < ninner) THEN
            CALL conn(ilevel, 1, s1=res)
        ELSE
            CALL conn(ilevel, 1, s1=res, forward=-1)
        END IF

        IF (ityp == 2) THEN
            CALL sipiter2_classic_level(ilevel, dp, res, sipue, sipun, siput)
        ELSE
            CALL sipiter2_hyperplane_level(ilevel, dp, res, sipue, sipun, &
                siput)
        END IF
    END SUBROUTINE sip


    SUBROUTINE sipiter1_classic_level(ilevel, res, rhs, siplw, sipls, siplb, &
            siplpr)
        ! Subroutine arguments
        INTEGER(intk), INTENT(in) :: ilevel
        TYPE(field_t), INTENT(inout) :: res
        TYPE(field_t), INTENT(in) :: rhs
        TYPE(field_t), INTENT(in) :: siplw
        TYPE(field_t), INTENT(in) :: sipls
        TYPE(field_t), INTENT(in) :: siplb
        TYPE(field_t), INTENT(in) :: siplpr

        ! Local variables
        INTEGER(intk) :: i, igrid
        INTEGER(intk) :: kk, jj, ii
        REAL(realk), POINTER, CONTIGUOUS :: lw(:, :, :), ls(:, :, :), &
            lb(:, :, :), lpr(:, :, :)
        REAL(realk), POINTER, CONTIGUOUS :: res_p(:, :, :), rhs_p(:, :, :)

        DO i = 1, nmygridslvl(ilevel)
            igrid = mygridslvl(i, ilevel)
            CALL get_mgdims(kk, jj, ii, igrid)

            CALL res%get_ptr(res_p, igrid)
            CALL rhs%get_ptr(rhs_p, igrid)

            CALL siplw%get_ptr(lw, igrid)
            CALL sipls%get_ptr(ls, igrid)
            CALL siplb%get_ptr(lb, igrid)
            CALL siplpr%get_ptr(lpr, igrid)

            CALL sipiter1_cl(kk, jj, ii, rhs_p, res_p, lw, ls, lb, lpr)
        END DO
    END SUBROUTINE sipiter1_classic_level


    SUBROUTINE sipiter1_hyperplane_level(ilevel, res_f, rhs_f, siplw, sipls, &
            siplb, siplpr)
        ! Subroutine arguments
        INTEGER(intk), INTENT(in) :: ilevel
        TYPE(field_t), INTENT(inout) :: res_f
        TYPE(field_t), INTENT(in) :: rhs_f
        TYPE(field_t), INTENT(in) :: siplw
        TYPE(field_t), INTENT(in) :: sipls
        TYPE(field_t), INTENT(in) :: siplb
        TYPE(field_t), INTENT(in) :: siplpr

        ! Local variables
        INTEGER(intk) :: i, igrid
        INTEGER(intk) :: kk, jj, ii, ip3

        ASSOCIATE( &
            lw => siplw%arr, &
            ls => sipls%arr, &
            lb => siplb%arr, &
            lpr => siplpr%arr, &
            res => res_f%arr, &
            rhs => rhs_f%arr, &
            mip => mip_hp_f%arr, &
            idx => idx_hp_f%arr)

#ifdef _MGLET_PROFILE_ANNOTATIONS_
        CALL profile_range_push("sipiter1_hp")
#endif
        !$omp target teams distribute private(i, igrid, kk, jj, ii, ip3)
        DO i = 1, nmygridslvl(ilevel)
            igrid = mygridslvl(i, ilevel)
            CALL get_mgdims(kk, jj, ii, igrid)
            CALL get_ip3(ip3, igrid)

            CALL sipiter1_hp(kk, jj, ii, rhs(ip3), res(ip3), lw(ip3), ls(ip3), &
                lb(ip3), lpr(ip3), mip(ip3), idx(ip3))
        END DO
        !$omp end target teams distribute
#ifdef _MGLET_PROFILE_ANNOTATIONS_
        CALL profile_range_pop()
#endif
        END ASSOCIATE
    END SUBROUTINE sipiter1_hyperplane_level


    SUBROUTINE sipiter2_classic_level(ilevel, dp, res, sipue, sipun, siput)
        ! Subroutine arguments
        INTEGER(intk), INTENT(in) :: ilevel
        TYPE(field_t), INTENT(inout) :: dp
        TYPE(field_t), INTENT(inout) :: res
        TYPE(field_t), INTENT(in) :: sipue
        TYPE(field_t), INTENT(in) :: sipun
        TYPE(field_t), INTENT(in) :: siput

        ! Local variables
        INTEGER(intk) :: i, igrid
        INTEGER(intk) :: kk, jj, ii
        REAL(realk), POINTER, CONTIGUOUS :: ue(:, :, :), un(:, :, :), &
            ut(:, :, :)
        REAL(realk), POINTER, CONTIGUOUS :: dp_p(:, :, :), res_p(:, :, :)

        DO i = 1, nmygridslvl(ilevel)
            igrid = mygridslvl(i, ilevel)
            CALL get_mgdims(kk, jj, ii, igrid)

            CALL dp%get_ptr(dp_p, igrid)
            CALL res%get_ptr(res_p, igrid)

            CALL sipue%get_ptr(ue, igrid)
            CALL sipun%get_ptr(un, igrid)
            CALL siput%get_ptr(ut, igrid)

            CALL sipiter2_cl(kk, jj, ii, dp_p, res_p, ue, un, ut)
        END DO
    END SUBROUTINE sipiter2_classic_level


    SUBROUTINE sipiter2_hyperplane_level(ilevel, dp_f, res_f, sipue_f, &
        sipun_f, siput_f)
        ! Subroutine arguments
        INTEGER(intk), INTENT(in) :: ilevel
        TYPE(field_t), INTENT(inout) :: dp_f
        TYPE(field_t), INTENT(inout) :: res_f
        TYPE(field_t), INTENT(in) :: sipue_f
        TYPE(field_t), INTENT(in) :: sipun_f
        TYPE(field_t), INTENT(in) :: siput_f

        ! Local variables
        INTEGER(intk) :: i, igrid
        INTEGER(intk) :: kk, jj, ii, ip3

        ASSOCIATE( &
            dp => dp_f%arr, &
            res => res_f%arr, &
            sipue => sipue_f%arr, &
            sipun => sipun_f%arr, &
            siput => siput_f%arr, &
            mip => mip_hp_f%arr, &
            idx => idx_hp_f%arr)

#ifdef _MGLET_PROFILE_ANNOTATIONS_
        CALL profile_range_push("sipiter2_hp")
#endif
        !$omp target teams distribute private(i, igrid, kk, jj, ii, ip3)
        DO i = 1, nmygridslvl(ilevel)
            igrid = mygridslvl(i, ilevel)
            CALL get_mgdims(kk, jj, ii, igrid)
            CALL get_ip3(ip3, igrid)

            CALL sipiter2_hp(kk, jj, ii, dp(ip3), res(ip3), sipue(ip3), &
                sipun(ip3), siput(ip3), mip(ip3), idx(ip3))
        END DO
        !$omp end target teams distribute
#ifdef _MGLET_PROFILE_ANNOTATIONS_
        CALL profile_range_pop()
#endif
        END ASSOCIATE
    END SUBROUTINE sipiter2_hyperplane_level


    SUBROUTINE maxabscal(maxabs, maxabslevel, phi_f)
        ! Subroutine arguments
        REAL(realk), INTENT(out) :: maxabs
        REAL(realk), INTENT(inout), ALLOCATABLE :: maxabslevel(:)
        TYPE(field_t), INTENT(in) :: phi_f

        ! Local variables
        REAL(realk) :: maxabsgrid(nmygrids)
        INTEGER(intk) :: imygrid, igrid, ilevel, ip3
        INTEGER(intk) :: kk, jj, ii

        maxabs = 0.0
        maxabslevel = 0.0

        ASSOCIATE(phi => phi_f%arr)

#ifdef _MGLET_PROFILE_ANNOTATIONS_
        CALL profile_range_push("maxabscal")
#endif

        !$omp target teams distribute private(imygrid, igrid, kk, jj, ii, ip3) &
        !$omp& map(from: maxabsgrid)
        DO imygrid = 1, nmygrids
            igrid = mygrids(imygrid)
            CALL get_mgdims(kk, jj, ii, igrid)
            CALL get_ip3(ip3, igrid)
            CALL maxabscal_grid(kk, jj, ii, maxabsgrid(imygrid), phi(ip3))
        END DO
        !$omp end target teams distribute

#ifdef _MGLET_PROFILE_ANNOTATIONS_
        CALL profile_range_pop()
#endif

        END ASSOCIATE

        DO imygrid = 1, nmygrids
            igrid = mygrids(imygrid)
            ilevel = level(igrid)
            maxabs = MAX(ABS(maxabsgrid(imygrid)*(2.0**(maxlevel-ilevel))), &
                maxabs)
            maxabslevel(ilevel) = MAX(maxabslevel(ilevel), maxabsgrid(imygrid))
        END DO
    END SUBROUTINE maxabscal


    PURE SUBROUTINE maxabscal_grid(kk, jj, ii, maxabs, phi)
        ! Subroutine arguments
        !$omp declare target
        INTEGER(intk), INTENT(in) :: kk, jj, ii
        REAL(realk), INTENT(out) :: maxabs
        REAL(realk), INTENT(in) :: phi(kk, jj, ii)

        ! Local variables
        INTEGER(intk) :: k, j, i

        maxabs = 0.0
        !$omp parallel do collapse(3) private(i, j, k) reduction(max:maxabs)
        DO i = 3, ii-2
            DO j = 3, jj-2
                DO k = 3, kk-2
                    maxabs = MAX(ABS(phi(k, j, i)), maxabs)
                END DO
            END DO
        END DO
        !$omp end parallel do
    END SUBROUTINE maxabscal_grid


    SUBROUTINE rescal(rhs_f, res_f)
        ! Subroutine arguments
        TYPE(field_t), INTENT(inout) :: rhs_f
        TYPE(field_t), INTENT(in) :: res_f

        ! Local variables
        INTEGER(intk) :: i, igrid
        INTEGER(intk) :: kk, jj, ii, ip3

        ASSOCIATE(rhs => rhs_f%arr, res => res_f%arr)

#ifdef _MGLET_PROFILE_ANNOTATIONS_
        CALL profile_range_push("rescal")
#endif
        !$omp target teams distribute private(i, igrid, kk, jj, ii, ip3)
        DO i = 1, nmygrids
            igrid = mygrids(i)
            CALL get_mgdims(kk, jj, ii, igrid)
            CALL get_ip3(ip3, igrid)

            CALL rescal_grid(kk, jj, ii, rhs(ip3), res(ip3))
        END DO
        !$omp end target teams distribute
#ifdef _MGLET_PROFILE_ANNOTATIONS_
        CALL profile_range_pop()
#endif
        END ASSOCIATE
    END SUBROUTINE rescal


    PURE SUBROUTINE rescal_grid(kk, jj, ii, rhs, res)
        !$omp declare target
        ! Subroutine arguments
        INTEGER(intk), INTENT(in) :: kk, jj, ii
        REAL(realk), INTENT(inout) :: rhs(kk, jj, ii)
        REAL(realk), INTENT(in) :: res(kk, jj, ii)

        ! Local variables
        INTEGER(intk) :: k, j, i

        ! TODO: Check if indices can be extended
        !$omp parallel do collapse(3) private(i, j, k)
        DO i = 3, ii-2
            DO j = 3, jj-2
                DO k = 3, kk-2
                    rhs(k, j, i) = rhs(k, j, i) + res(k, j, i)
                END DO
            END DO
        END DO
        !$omp end parallel do
    END SUBROUTINE rescal_grid


    SUBROUTINE mgpcorr(u_f, v_f, w_f, p_f, dp_f, rfak)
        ! Subroutine arguments
        TYPE(field_t), INTENT(inout) :: u_f, v_f, w_f, p_f
        TYPE(field_t), INTENT(in) :: dp_f
        REAL(realk), INTENT(in) :: rfak

        ! Local variables
        INTEGER(intk) :: i, igrid, ip3, ip1x, ip1y, ip1z
        INTEGER(intk) :: kk, jj, ii

        TYPE(field_t), POINTER :: rdx_f
        TYPE(field_t), POINTER :: rdy_f
        TYPE(field_t), POINTER :: rdz_f
        TYPE(field_t), POINTER :: bp_f

        CALL get_field(rdx_f, "RDX")
        CALL get_field(rdy_f, "RDY")
        CALL get_field(rdz_f, "RDZ")
        CALL get_field(bp_f, "BP")

        ASSOCIATE(u => u_f%arr, v => v_f%arr, w => w_f%arr, p => p_f%arr, &
            dp => dp_f%arr, bp => bp_f%arr, &
            rdx => rdx_f%arr, rdy => rdy_f%arr, rdz => rdz_f%arr)

#ifdef _MGLET_PROFILE_ANNOTATIONS_
        CALL profile_range_push("mgpcorr")
#endif

        !$omp target teams distribute private(i, igrid, kk, jj, ii, ip3, ip1x, &
        !$omp& ip1y, ip1z)
        DO i = 1, nmygrids
            igrid = mygrids(i)
            CALL get_mgdims(kk, jj, ii, igrid)
            CALL get_ip3(ip3, igrid)
            CALL get_ip1x(ip1x, igrid)
            CALL get_ip1y(ip1y, igrid)
            CALL get_ip1z(ip1z, igrid)

            !$omp parallel
            CALL mgpcorr_grid(kk, jj, ii, u(ip3), v(ip3), w(ip3), p(ip3), &
                dp(ip3), bp(ip3), rdx(ip1x), rdy(ip1y), rdz(ip1z), rfak)
            !$omp end parallel
        END DO
        !$omp end target teams distribute

#ifdef _MGLET_PROFILE_ANNOTATIONS_
        CALL profile_range_pop()
#endif

        END ASSOCIATE
    END SUBROUTINE mgpcorr


    PURE SUBROUTINE mgpcorr_grid(kk, jj, ii, u, v, w, p, dp, bp, rdx, rdy, &
            rdz, rfak)
        !$omp declare target
        ! Subroutine arguments
        INTEGER(intk), INTENT(in) :: kk, jj, ii
        REAL(realk), INTENT(inout) :: u(kk, jj, ii)
        REAL(realk), INTENT(inout) :: v(kk, jj, ii)
        REAL(realk), INTENT(inout) :: w(kk, jj, ii)
        REAL(realk), INTENT(inout) :: p(kk, jj, ii)
        REAL(realk), INTENT(in) :: dp(kk, jj, ii)
        REAL(realk), INTENT(in) :: bp(kk, jj, ii)
        REAL(realk), INTENT(in) :: rdx(ii)
        REAL(realk), INTENT(in) :: rdy(jj)
        REAL(realk), INTENT(in) :: rdz(kk)
        REAL(realk), INTENT(in) :: rfak

        ! Local variables
        INTEGER(intk) :: k, j, i

        !$omp do collapse(3)
        DO i = 2, ii-1
            DO j = 2, jj-1
                DO k = 2, kk-1
                    p(k, j, i) = p(k, j, i) + dp(k, j, i)*bp(k, j, i)
                END DO
            END DO
        END DO
        !$omp end do

        !$omp do collapse(3)
        DO i = 2, ii-2
            DO j = 3, jj-2
                DO k = 3, kk-2
                    u(k, j, i) = u(k, j, i) &
                        + (dp(k, j, i) - dp(k, j, i+1)) &
                        *bp(k, j, i)*bp(k, j, i+1)*rdx(i)*rfak
                END DO
            END DO
        END DO
        !$omp end do

        !$omp do collapse(3)
        DO i = 3, ii-2
            DO j = 2, jj - 2
                DO k = 3, kk-2
                    v(k, j, i) = v(k, j, i) &
                        + (dp(k, j, i) - dp(k, j+1, i)) &
                        *bp(k, j, i)*bp(k, j+1, i)*rdy(j)*rfak
                END DO
            END DO
        END DO
        !$omp end do

        !$omp do collapse(3)
        DO i = 3, ii-2
            DO j = 3, jj-2
                DO k = 2, kk-2
                    w(k, j, i) = w(k, j, i) &
                        + (dp(k, j, i) - dp(k+1, j, i)) &
                        *bp(k, j, i)*bp(k+1, j, i)*rdz(k)*rfak
                END DO
            END DO
        END DO
        !$omp end do
    END SUBROUTINE mgpcorr_grid


    SUBROUTINE accumulate_pcorr(dp, hilf)
        ! Subroutine arguments
        TYPE(field_t), INTENT(inout) :: dp
        TYPE(field_t), INTENT(in) :: hilf

        ! Local variables
        INTEGER(intk) :: i, n

        IF (SIZE(dp%arr) /= SIZE(hilf%arr)) CALL errr(__FILE__, __LINE__)

        n = SIZE(dp%arr)

        ASSOCIATE(dp => dp%arr, hilf => hilf%arr)

#ifdef _MGLET_PROFILE_ANNOTATIONS_
        CALL profile_range_push("accumulate_pcorr")
#endif
        !$omp target teams loop
        DO i = 1, n
            dp(i) = dp(i) + hilf(i)
        END DO
        !$omp end target teams loop
#ifdef _MGLET_PROFILE_ANNOTATIONS_
        CALL profile_range_pop()
#endif
        END ASSOCIATE
    END SUBROUTINE accumulate_pcorr
END MODULE pressuresolver_mod
