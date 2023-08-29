MODULE timeloop_mod
    USE MPI_f08
    USE core_mod
    USE flow_mod, ONLY: timeintegrate_flow, itinfo_flow
    USE scalar_mod, ONLY: timeintegrate_scalar
    USE runinfo_mod

    IMPLICIT NONE(type, external)
    PRIVATE

    ! Flag used to skip timeintegration entirely
    LOGICAL :: skip_timeloop = .FALSE.

    ! Filename of file to check for runtime contraint
    CHARACTER(len=*), PARAMETER :: fort57name = 'fort.57'

    ! Signal number to ctach:
    ! Options:
    !   2: SIGINT (Ctrl+C)
    !  10: SIGUSR1
    !  24: SIGXCPU
    INTEGER(int32), PARAMETER :: signum(2) = [10, 24]

    ! Integer flag to indicate if signal is received
    INTEGER(int32), VOLATILE :: received_signal = 0

    ! For walltime book-keeping
    REAL(real64), PROTECTED :: wtime0 = 0.0
    REAL(real64), PROTECTED :: wtime1 = 0.0
    INTEGER(intk), PROTECTED :: it1 = 0
    INTEGER(int64) :: ncellstot

    ! Stopping criterions
    INTEGER(intk) :: mtstep, itint
    REAL(realk) :: tend
    LOGICAL :: has_tend, has_itint

    ! For timekeeping
    REAL(realk), PROTECTED :: timeph = 0.0
    REAL(realk), PROTECTED :: tstat = 0.0
    REAL(realk), PROTECTED :: dt = 0.0
    REAL(realk), PROTECTED :: targetcflmax = 0.0
    INTEGER(intk), PROTECTED :: ittot = 0
    INTEGER(intk), PROTECTED :: itstep = 0
    INTEGER(intk), PROTECTED :: itinfo = 0
    INTEGER(intk), PROTECTED :: itsamp = 0
    INTEGER(intk), PROTECTED :: itstop = 0
    TYPE(timekeeper_t) :: timekeeper

    ! For checkpointing
    INTEGER(intk), PROTECTED :: itcheck = 0   ! Min. checkpointing frequency
    INTEGER(intk), PROTECTED :: istepchk = 0  ! Steps since last checkpoint

    ! Time log file
    CHARACTER(len=*), PARAMETER :: timefile = logdir//"/time.log"

    ! Main time integration scheme - flow and scalar follow this
    TYPE(rk_2n_t), PROTECTED :: rkscheme

    PUBLIC :: init_timeloop, finish_timeloop, timeloop, timeph, dt, &
        ittot, itstep, itinfo, itsamp, itstop, wtime0
CONTAINS
    SUBROUTINE init_timeloop()
#if defined __INTEL_COMPILER
        USE ifport

        INTEGER(int64) :: ierr
#endif
#if defined __GFORTRAN__
        INTEGER(intk) :: ierr
#endif

        TYPE(config_t) :: timeconf
        INTEGER(intk) :: i, igrid, logunit
        INTEGER(intk) :: kk, jj, ii
        CHARACTER(len=16) :: rkmethod
        LOGICAL :: exists
        REAL(realk) :: dt2

        ! Initialize signals
        received_signal = 0
        DO i = 1, SIZE(signum)
#if defined __INTEL_COMPILER
            ierr = SIGNAL(signum(i), sighandler, -1)
#endif
#if defined __GFORTRAN__
            CALL SIGNAL(signum(i), sighandler, ierr)
#endif
        END DO

        ! Read configuration values - if not exists no timeintegration is
        ! performed
        skip_timeloop = .FALSE.
        IF (.NOT. fort7%exists("/time")) THEN
            skip_timeloop = .TRUE.
            IF (myid == 0) THEN
                WRITE(*, '("SKIPPING TIMELOOP")')
                WRITE(*, '()')
            END IF
            RETURN
        END IF
        CALL fort7%get(timeconf, "/time")

        ! Required values
        CALL timeconf%get_value("/mtstep", mtstep)
        CALL timeconf%get_value("/dt", dt)
        IF (dt <= 0) THEN
            WRITE(*, '("Timestep must be positive: ", F15.7)') dt
            CALL errr(__FILE__, __LINE__)
        END IF

        ! Optional variables
        CALL timeconf%get_value("/itsamp", itsamp, 1)
        CALL timeconf%get_value("/itinfo", itinfo, 10)
        CALL timeconf%get_value("/itint", itint, 0, has_itint)
        CALL timeconf%get_value("/tend", tend, 0.0, has_tend)
        CALL timeconf%get_value("/tstat", tstat, 0.0)
        CALL timeconf%get_value("/itcheck", itcheck, 0)
        CALL timeconf%get_value("/targetcflmax", targetcflmax, 0.0)
        CALL timeconf%get_value("/rkmethod", rkmethod, "williamson")
        CALL rkscheme%init(rkmethod)

        IF (dcont) THEN
            CALL read_runinfo(ittot, timeph, dt2)
            ! The timestep is only taken from the restart file when the
            ! targetcfl-feature is in use!
            IF (targetcflmax > 0.0) THEN
                dt = dt2
            END IF
        ELSE
            timeph = 0.0
            ittot = 0
        END IF

        ! Timekeeper for accurate time-keeping
        CALL timekeeper%init(timeph)

        ! itstop is the interval in which the stopping criterions (besides
        ! MTSTEP which is a hard criterion) is checked
        ! TODO: itstop = lcm(itstop, itsamp_probes) and do on neccesary???
        itstop = itsamp

        ! Compute overall number of (active) grid cells for performance
        ! number calculation
        ncellstot = 0
        DO i = 1, nmygrids
            igrid = mygrids(i)
            CALL get_mgdims(kk, jj, ii, igrid)
            ncellstot = ncellstot + INT((kk-4)*(jj-4)*(ii-4), int64)
        END DO
        CALL MPI_Allreduce(MPI_IN_PLACE, ncellstot, 1, MPI_INTEGER8, MPI_SUM, &
            MPI_COMM_WORLD)

        ! Sometimes one is not careful when restarting, and forget to
        ! increase tend in case of a restart. Then we have a situation where
        ! timeph > tend beofre we have integrated anyting. Set mtstep to 0
        ! in this case and write a warning:
        IF (timeph >= tend .AND. has_tend) THEN
            IF (myid == 0) THEN
                WRITE(*, '("WARNING: tend > timeph at start of simulation, ")')
                WRITE(*, '("skipping time integration by setting mtstep = 0")')
                WRITE(*, '("tend = ", F15.6)') tend
                WRITE(*, '("timeph = ", F15.6)') timeph
                WRITE(*, '()')
            END IF
            mtstep = 0
            RETURN
        END IF

        ! Time log file
        IF (myid == 0) THEN
            ! Create log file if log file does not exist
            INQUIRE(FILE=timefile, EXIST=exists)
            IF (.NOT. exists .OR. .NOT. dcont) THEN
                OPEN(NEWUNIT=logunit, FILE=timefile)
                WRITE(logunit, '(A9, A18, 4A14)') &
                    "#  ITTOT", "TIME", "DT", "WTIME", "DELTAWT", "TPC"
                CLOSE(logunit)
            END IF
        END IF

        ! Late plugin initialization
        CALL init_late_plugins(ittot, mtstep, itint, timeph, dt, tend)

        ! Initialize statistics
        CALL init_statistics()
    END SUBROUTINE init_timeloop


    SUBROUTINE finish_timeloop()
        ! Local variables
        REAL(real64) :: wtime, deltawt

        IF (skip_timeloop) RETURN

        wtime = MPI_Wtime()
        deltawt = wtime - wtime0
        CALL write_runinfo(itstep, ittot, timeph, dt, targetcflmax, deltawt, &
            ncellstot)

        CALL finish_statistics()
    END SUBROUTINE finish_timeloop


    SUBROUTINE timeloop()
        ! Local variables
        LOGICAL :: stop_now, allow_checkpoint
        INTEGER(intk) :: irk, exploded
        REAL(realk) :: cflmax

        IF (skip_timeloop) RETURN

        ! Initialize and set timer
        CALL set_timer(2, 'TIMELOOP')
        CALL start_timer(2)

        ! Start walltime
        wtime0 = MPI_Wtime()
        wtime1 = wtime0
        it1 = 0

        ! Initial printout
        IF (ittot == 0) THEN
            CALL itinfo_time(itstep, ittot, timeph, dt)
            CALL itinfo_flow(itstep, ittot, timeph, dt, cflmax, exploded)
            IF (myid == 0) WRITE(*, '()')
        END IF

        timeintegration: DO itstep = 1, mtstep
            ! Global RK loop for tightly coupled quantities like flow and
            ! scalar transport
            rkloop: DO irk = 1, rkscheme%nrk
                CALL timeintegrate_scalar(itstep, ittot, timeph, dt, irk, &
                    rkscheme)
                CALL timeintegrate_flow(itstep, ittot, timeph, dt, irk, &
                    rkscheme)
            END DO rkloop

            ! Timeintegrate, _before_ time is globally incremented
            CALL timeintegrate_plugins(itstep, ittot, timeph, dt)

            ! Increment time values - timestepping finished!
            ittot = ittot + 1
            CALL timekeeper%add_to_time(dt)
            CALL timekeeper%get_time(timeph)

            ! Print to terminal (itinfo frequency)
            !
            ! Comes _before_ general plugin postprocessing, to have the
            ! terminal output _after_ the general itinfo information
            IF (MOD(ittot, itinfo) == 0) THEN
                CALL itinfo_time(itstep, ittot, timeph, dt)
                ! CALL itinfo_scalar(itstep, ittot, timeph, dt, exploded)
                CALL itinfo_flow(itstep, ittot, timeph, dt, cflmax, exploded)
                IF (myid == 0) WRITE(*, '()')

                ! Call plugins
                CALL itinfo_plugins(itstep, ittot, timeph, dt)
            END IF

            ! Sample statistics
            IF (timeph > tstat .AND. MOD(ittot, itsamp) == 0) THEN
                CALL sample_statistics(dt)
            END IF

            ! Postprocess, _after_ time is globally incremented
            CALL postprocess_plugins(itstep, ittot, timeph, dt)

            ! Explosions trigger immediate stop
            IF (exploded > 0) THEN
                IF (myid == 0) THEN
                    WRITE(*, '("STOP SIMULATION, SOLUTION EXPLODED")')
                    WRITE(*, '()')
                END IF
                EXIT timeintegration
            END IF

            ! Check stop criterions
            IF (MOD(ittot, itstop) == 0) THEN
                stop_now = check_stop(itstep, ittot, timeph)
                IF (stop_now) EXIT timeintegration
            END IF

            ! Checkpointing after all plugins have processed, and after all stop
            ! criterions have been checked (to avoid writing a checkpoint just
            ! before stopping the simulation), but before the new DT for the
            ! next timestep is set.
            istepchk = istepchk + 1
            IF (istepchk >= itcheck .AND. itcheck > 0) THEN
                CALL can_checkpoint_plugins(allow_checkpoint)
                IF (allow_checkpoint) THEN
                    IF (myid == 0) WRITE(*, '("Writing checkpoint...")')
                    CALL checkpoint_fields()
                    CALL checkpoint_plugins()
                    istepchk = 0
                END IF
            END IF

            ! Adjust timestep at end of all postprocessing. Some tools (probes,
            ! statistics (computing PtPt_AVG)) need the DT, and then we must
            ! wait until here to adjust DT
            IF (MOD(ittot, itinfo) == 0) THEN
                IF (targetcflmax > 0.0 .AND. timeph < tstat) THEN
                    CALL adjust_timestep(cflmax)
                END IF
            END IF
        END DO timeintegration

        CALL stop_timer(2)
    END SUBROUTINE timeloop


    SUBROUTINE itinfo_time(itstep, ittot, timeph, dt)
        ! Subroutine arguments
        INTEGER(intk), INTENT(in) :: itstep
        INTEGER(intk), INTENT(in) :: ittot
        REAL(realk), INTENT(in) :: timeph
        REAL(realk), INTENT(in) :: dt

        ! Local variables
        REAL(real64) :: wtime, delta, tpc
        INTEGER :: logunit

        wtime = MPI_Wtime()
        IF (itstep > 0) THEN
            delta = (wtime - wtime1)/REAL(itstep-it1, real64)
        ELSE
            delta = 0.0
        END IF

        ! Compute performance number
        tpc = delta*REAL(numprocs, real64)/REAL(ncellstot, real64)

        IF (myid == 0) THEN
            ! CLI print
            WRITE(*,'("ITSTEP, WTIME, DELTA: ", I8, F16.3, F16.3)') &
                itstep, wtime - wtime0, delta
            WRITE(*, '("ITTOT, TIME, DT: ", I13, F16.8, ES16.4)') &
                ittot, timeph, dt

            ! Logfile print
            OPEN(NEWUNIT=logunit, FILE=timefile, POSITION="APPEND")
            WRITE(logunit, '(I9, ES18.10, 4ES14.5)') &
                ittot, timeph, dt, wtime - wtime0, delta, tpc
            CLOSE(logunit)
        END IF
        wtime1 = wtime
        it1 = itstep
    END SUBROUTINE itinfo_time


    LOGICAL FUNCTION check_stop(itstep, ittot, timeph)
        ! Subroutine arguments
        INTEGER(intk), INTENT(in) :: itstep
        INTEGER(intk), INTENT(in) :: ittot
        REAL(realk), INTENT(in) :: timeph

        ! Local variables
        LOGICAL :: fileexist
        INTEGER(intk) :: ios, istop
        REAL(real64) :: wtime, cpdiff, tstop

        check_stop = .FALSE.

        IF (itstep >= mtstep) THEN
            check_stop = .TRUE.
            IF (myid == 0) THEN
                WRITE(*, '("STOP SIMULATION, MTSTEP REACHED", I13, I13)') &
                    itstep, mtstep
                WRITE(*, '()')
            END IF
        END IF

        IF (has_itint .AND. ittot >= itint) THEN
            check_stop = .TRUE.
            IF (myid == 0) THEN
                WRITE(*, '("STOP SIMULATION, ITINT REACHED", I13, I13)') &
                    itint, ittot
                WRITE(*, '()')
            END IF
        END IF

        IF (has_tend .AND. timeph >= tend) THEN
            check_stop = .TRUE.
            IF (myid == 0) THEN
                WRITE(*, '("STOP SIMULATION, TEND REACHED", F16.8, F16.8)') &
                    tend, timeph
                WRITE(*, '()')
            END IF
        END IF

        ! Tee following stopping criterions are checked only by rank 0
        istop = 0

        ! fort.57 stopping criterion, only rank 0 check
        IF (myid == 0) THEN
            fileexist = .FALSE.
            INQUIRE(FILE=fort57name, EXIST=fileexist)
            IF (fileexist) THEN
                OPEN(57, FILE=fort57name, ACTION="READ")
                READ(57, *, IOSTAT=ios) tstop
                CLOSE(57)

                wtime = MPI_Wtime()
                cpdiff = wtime - wtime0

                IF ((tstop <= cpdiff) .AND. (ios == 0)) THEN
                    istop = 1
                END IF
            END IF
        END IF

        ! Any rank that receive a signal will trigger a stop
        IF (received_signal > 0) THEN
            istop = 2
        END IF

        ! Broadcast stop message to all processes
        CALL MPI_Allreduce(MPI_IN_PLACE, istop, 1, mglet_mpi_int, MPI_MAX, &
            MPI_COMM_WORLD)

        IF (istop == 1) THEN
            check_stop = .TRUE.
            IF (myid == 0) THEN
                WRITE(*, '("STOP SIMULATION, WTIME REACHED")')
                WRITE(*, '()')
            END IF
        END IF

        IF (istop == 2) THEN
            check_stop = .TRUE.
            IF (myid == 0) THEN
                WRITE(*, '("STOP SIMULATION, SIGNAL RECEIVED")')
                WRITE(*, '()')
            END IF
        END IF
    END FUNCTION check_stop


    SUBROUTINE adjust_timestep(cflmax)
        ! Subroutine arguments
        REAL(realk), INTENT(in) :: cflmax

        ! Local variables
        INTEGER(intk), PARAMETER :: ndtrelax = 500
        REAL(realk) :: ks, ra, ki

        ks = cflmax/dt
        ki = 1.0/ks/ndtrelax*itinfo
        ra = targetcflmax - cflmax
        dt = dt + ki*ra
    END SUBROUTINE adjust_timestep


    ! Routine to hanle signals
    INTEGER(int32) FUNCTION sighandler(sig)
        INTEGER(kind=int32) :: sig
        received_signal = sig
        sighandler = 0
    END FUNCTION sighandler


    SUBROUTINE checkpoint_fields()
        ! Subroutine arguments
        ! none...

        ! Local variables
        REAL(real64) :: wtime, deltawt

        ! Open output file for writing
        CALL fields_begin_write()

        ! Write runinfo-table
        ! TODO: discuss the values written here. Maybe better to use the
        ! walltime since last checkpoint? But then, what to write in the end?
        wtime = MPI_Wtime()
        deltawt = wtime - wtime0
        CALL write_runinfo(itstep, ittot, timeph, dt, targetcflmax, deltawt, &
            ncellstot)

        ! Write field data
        CALL fields_write()

        ! Close file
        CALL fields_end_rw()
    END SUBROUTINE checkpoint_fields
END MODULE timeloop_mod
