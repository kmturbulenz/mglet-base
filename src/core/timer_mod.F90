MODULE timer_mod
    USE comms_mod
    USE err_mod
    USE precision_mod
    USE buildinfo_mod, ONLY: mglet_dbg_envvar
    USE, INTRINSIC :: ISO_C_BINDING, ONLY: c_double, c_int, c_long_long, &
        c_f_pointer
    USE MPI_f08

    IMPLICIT NONE (type, external)
    PRIVATE

    INTEGER(intk), PARAMETER :: maxtimers = 1000
    INTEGER(intk), PARAMETER :: stackdepth = 32
    INTEGER(intk), PARAMETER :: desclen = 32

    CHARACTER(len=*), PARAMETER :: filename = 'mglet-perf-report.txt'

    TYPE timer_t
        ! For measuring "inclusive" time
        REAL(real64) :: tic = -1.0
        REAL(real64) :: meanTime = 0.0

        ! For measuring "exclusive" time
        REAL(real64) :: exTic = -1.0
        REAL(real64) :: exSum = 0.0
        REAL(real64) :: exMeanTime = 0.0

        ! Number of instances
        INTEGER(int64) :: n = 0
        CHARACTER(len=desclen) :: desc
    END TYPE timer_t

    TYPE, BIND(C) :: timer_stat_t
        REAL(c_double) :: sumTime
        REAL(c_double) :: minTime
        REAL(c_double) :: maxTime
        REAL(c_double) :: avgTime
        REAL(c_double) :: avgN
        INTEGER(c_long_long) :: n
    END TYPE timer_stat_t

    TYPE(timer_t) :: timers(maxtimers)

    INTEGER :: stack(stackdepth)
    INTEGER :: stackpos = 0

    TYPE(MPI_Datatype) :: mpitype
    TYPE(MPI_Op) :: mpiop

    PUBLIC :: init_timer, set_timer, start_timer, stop_timer, finish_timer, &
        stop_all

CONTAINS
    SUBROUTINE init_timer()
        INTEGER(intk) :: idx

        DO idx = 1, maxtimers
            timers(idx)%tic = -1.0
            timers(idx)%meanTime = 0.0
            timers(idx)%exTic = -1.0
            timers(idx)%exSum = 0.0
            timers(idx)%exMeanTime = 0.0
            timers(idx)%n = 0
            timers(idx)%desc = REPEAT(' ', desclen)
        END DO
    END SUBROUTINE init_timer


    ! Initialize and restart timer
    SUBROUTINE set_timer(idx, desc)
        ! Subroutine arguments
        INTEGER, INTENT(in) :: idx
        CHARACTER(len=*), INTENT(in) :: desc

        ! Local variables
        ! none...

        IF (idx > maxtimers) CALL errr(__FILE__, __LINE__)

        timers(idx)%tic = -1.0
        timers(idx)%meanTime = 0.0
        timers(idx)%exTic = -1.0
        timers(idx)%exSum = 0.0
        timers(idx)%exMeanTime = 0.0
        timers(idx)%n = 0
        timers(idx)%desc = desc(1:MIN(LEN(desc), desclen))
    END SUBROUTINE set_timer


    ! Start timer
    SUBROUTINE start_timer(idx)
        ! Subroutine arguments
        INTEGER, INTENT(in) :: idx

        ! Local variables
        INTEGER :: stackIdx
        REAL(real64) :: tic

        IF (idx > maxtimers) CALL errr(__FILE__, __LINE__)

        IF (timers(idx)%tic < 0.0) THEN
            tic = MPI_Wtime()

            timers(idx)%tic = tic
            timers(idx)%exTic = tic

            ! Stop previous stack
            IF (stackpos > 0) THEN
                stackIdx = stack(stackpos)
                timers(stackIdx)%exSum = timers(stackIdx)%exSum + &
                    (tic - timers(stackIdx)%exTic)
            END IF

            stackpos = stackpos + 1
            stack(stackpos) = idx
        ELSE
            CALL err_abort(255, "Timer"//CHAR(idx)//"not stopped", __FILE__, &
                __LINE__)
        END IF
    END SUBROUTINE start_timer


    ! Stop timer and update running average
    SUBROUTINE stop_timer(idx)
        ! Subroutine arguments
        INTEGER, INTENT(in) :: idx

        ! Local variables
        INTEGER :: stackIdx
        REAL(real64) :: toc, elapsed

        IF (idx > maxtimers) CALL errr(__FILE__, __LINE__)
        IF (stack(stackpos) /= idx) THEN
           write(*, *) "Timer", idx, "stack",stack(stackpos)
           CALL errr(__FILE__, __LINE__)
        END IF
        IF (timers(idx)%tic < 0.0) THEN
            write(*, *) "Timer", idx, "not started."
            CALL errr(__FILE__, __LINE__)
        END IF

        toc = MPI_Wtime()
        timers(idx)%n = timers(idx)%n + 1

        ! Update current "inclusivce" time
        elapsed = toc - timers(idx)%tic
        timers(idx)%meanTime = timers(idx)%meanTime + &
            (elapsed - timers(idx)%meanTime)/REAL(timers(idx)%n, real64)
        timers(idx)%tic = -1.0

        ! Update current "exclusive" time
        elapsed = toc - timers(idx)%exTic
        timers(idx)%exSum = timers(idx)%exSum + elapsed
        timers(idx)%exMeanTime = timers(idx)%exMeanTime + &
            (timers(idx)%exSum - timers(idx)%exMeanTime)/REAL(timers(idx)%n, real64)
        timers(idx)%exSum = 0.0

        ! Start previous stack timer again
        stackpos = stackpos - 1
        IF (stackpos > 0) THEN
            stackIdx = stack(stackpos)
            timers(stackIdx)%exTic = toc
        END IF
    END SUBROUTINE stop_timer


    ! Stop all timers (for use from clean_exit)
    SUBROUTINE stop_all()
        ! Local variables
        INTEGER(intk) :: i

        DO i = 1, stackpos
            CALL stop_timer(stack(stackpos))
        END DO
    END SUBROUTINE stop_all


    ! Print result
    SUBROUTINE finish_timer()
        ! Local variables
        INTEGER(intk) :: i
        TYPE(timer_stat_t) :: inStat, stat
        INTEGER :: wu

        IF (stackpos > 0) THEN
            WRITE(*, *) "Timers not stopped. Stackpos:", stackpos
            WRITE(*, *) "Open timers", stack(1:stackpos)
            CALL errr(__FILE__, __LINE__)
        END IF

        CALL create_reduction()

        IF (myid == 0) THEN
            OPEN(newunit=wu, file=filename)
            WRITE(wu, *) "TIMING RESULTS:"
            WRITE(wu, *) "  Average time per instance."
            WRITE(wu, *) "  Minimum and maximum average instance time per process."
            WRITE(wu, *) "  Average number of instances per process."
            WRITE(wu, *) "  Average total time."
            WRITE(wu, *) ""
            WRITE(wu, '(A13, ES12.3)') "   MPI_Wtick:", MPI_Wtick()
            WRITE(wu, *) ""
            WRITE(wu, *) "INCLUSIVE TIME:"
            WRITE(wu, '(A6, 1X, A32, 1X, A12, A12, A12, A12, A12)') "idx", &
                "region", "average", "minimum", "maximum", "instances", "total"
        END IF

        DO i = 1, maxtimers
            inStat%sumTime = timers(i)%meanTime*timers(i)%n
            inStat%minTime = timers(i)%meanTime
            inStat%maxTime = timers(i)%meanTime
            inStat%n = timers(i)%n

            IF (inStat%n > 0) THEN
                inStat%avgTime = inStat%sumTime/inStat%n
            ELSE
                inStat%avgTime = 0
            END IF
            inStat%avgN = inStat%n/numprocs

            ! Explicitly initialize stat%n to avoid valgrind errors in the
            ! IF below (otherwise it is uninitialized for ranks /= 0)
            stat%n = 0
            CALL MPI_Reduce(inStat, stat, 1, mpitype, mpiop, 0, &
                MPI_COMM_WORLD)

            IF (stat%n > 0 .AND. myid == 0) THEN
                WRITE(wu, '(I6, 1X, A32, 1X, ES12.3, ES12.3, ES12.3, ES12.3, ES12.3)') &
                    i, TRIM(timers(i)%desc), stat%avgTime, stat%minTime, &
                    stat%maxTime, stat%avgN, stat%avgTime*stat%avgN
            END IF
        END DO

        IF (myid == 0) THEN
            WRITE(wu, *) ""
            WRITE(wu, *) "EXCLUSIVE TIME:"
            WRITE(wu, '(A6, 1X, A32, 1X, A12, A12, A12, A12, A12)') "idx", &
                "region", "average", "minimum", "maximum", "instances", "total"
        END IF

        DO i = 1, maxtimers
            inStat%sumTime = timers(i)%exMeanTime*timers(i)%n
            inStat%minTime = timers(i)%exMeanTime
            inStat%maxTime = timers(i)%exMeanTime
            inStat%n = timers(i)%n

            IF (inStat%n > 0) THEN
                inStat%avgTime = inStat%sumTime/inStat%n
            ELSE
                inStat%avgTime = 0
            END IF
            inStat%avgN = inStat%n/numprocs

            ! Explicitly initialize stat%n to avoid valgrind errors in the
            ! IF below (otherwise it is uninitialized for ranks /= 0)
            stat%n = 0
            CALL MPI_Reduce(inStat, stat, 1, mpitype, mpiop, 0, &
                MPI_COMM_WORLD)

            IF (stat%n > 0 .AND. myid == 0) THEN
                WRITE(wu, '(I6, 1X, A32, 1X, ES12.3, ES12.3, ES12.3, ES12.3, ES12.3)') &
                    i, TRIM(timers(i)%desc), stat%avgTime, stat%minTime, &
                    stat%maxTime, stat%avgN, stat%avgTime*stat%avgN
            END IF
        END DO

        IF (myid == 0) THEN
            close(wu)
        END IF

        CALL MPI_Op_free(mpiop)
        CALL MPI_Type_free(mpitype)

        IF (INDEX(mglet_dbg_envvar, "individual_timers") > 0) THEN
            CALL print_individual_timers()
        END IF
    END SUBROUTINE finish_timer


    SUBROUTINE print_individual_timers()
        ! Subroutine arguments
        ! none...

        ! Local variables
        INTEGER(intk) :: i
        CHARACTER(len=256) :: filename
        INTEGER :: unit

        WRITE(filename, '("mglet-perf-rank-", I0, ".txt")') myid
        OPEN(newunit=unit, file=filename)

        WRITE(unit, '("TIMING RESULT FOR RANK: ", I0)') myid
        WRITE(unit, '()')
        WRITE(unit, '("INCLUSIVE TIME:")')
        WRITE(unit, '(A6, 1X, A32, 1X, A12, A12, A12, A12, A12)') "idx", &
            "region", "average", "instances", "total"
        DO i = 1, maxtimers
            IF (timers(i)%n > 0) THEN
                WRITE(unit, '(I6, 1X, A32, 1X, ES12.3, I12, ES12.3)') &
                    i, TRIM(timers(i)%desc), timers(i)%meanTime, &
                    timers(i)%n, timers(i)%meanTime*timers(i)%n
            END IF
        END DO
        WRITE(unit, '()')
        WRITE(unit, '("EXCLUSIVE TIME:")')
        WRITE(unit, '(A6, 1X, A32, 1X, A12, A12, A12, A12, A12)') "idx", &
            "region", "average", "instances", "total"
        DO i = 1, maxtimers
            IF (timers(i)%n > 0) THEN
                WRITE(unit, '(I6, 1X, A32, 1X, ES12.3, I12, ES12.3)') &
                    i, TRIM(timers(i)%desc), timers(i)%exMeanTime, &
                    timers(i)%n, timers(i)%exMeanTime*timers(i)%n
            END IF
        END DO

        CLOSE(unit)
    END SUBROUTINE print_individual_timers


    SUBROUTINE create_reduction()
        ! Local variables
        INTEGER(int32) :: blocklen(6)
        TYPE(MPI_Datatype) :: types(6)
        INTEGER(mpi_address_kind) :: base, disp(6)
        TYPE(timer_stat_t) :: foo

        blocklen = 1

        CALL MPI_Get_address(foo%sumTime, disp(1))
        CALL MPI_Get_address(foo%minTime, disp(2))
        CALL MPI_Get_address(foo%maxTime, disp(3))
        CALL MPI_Get_address(foo%avgTime, disp(4))
        CALL MPI_Get_address(foo%avgN, disp(5))
        CALL MPI_Get_address(foo%n, disp(6))

        base = disp(1)
        disp(1) = disp(1) - base
        disp(2) = disp(2) - base
        disp(3) = disp(3) - base
        disp(4) = disp(4) - base
        disp(5) = disp(5) - base
        disp(6) = disp(6) - base

        types(1) = MPI_REAL8
        types(2) = MPI_REAL8
        types(3) = MPI_REAL8
        types(4) = MPI_REAL8
        types(5) = MPI_REAL8
        types(6) = MPI_INTEGER8

        CALL MPI_Type_create_struct(6, blocklen, disp, types, mpitype)
        CALL MPI_Type_commit(mpitype)

        CALL MPI_Op_create(reduce_time, .TRUE., mpiop)
    END SUBROUTINE create_reduction


    SUBROUTINE reduce_time(invec, inoutvec, length, datatype)
        ! Subroutine arguments
        TYPE(C_PTR), VALUE :: invec
        TYPE(C_PTR), VALUE :: inoutvec
        INTEGER(int32) :: length        ! Declaring INTENT makes it incompatible
        TYPE(MPI_Datatype) :: datatype  ! with MPI_f08 MPI_User_function

        ! Local variables
        INTEGER(c_int) :: i
        TYPE(timer_stat_t), POINTER :: idata(:)
        TYPE(timer_stat_t), POINTER :: iodata(:)

        CALL C_F_POINTER(invec, idata, [length])
        CALL C_F_POINTER(inoutvec, iodata, [length])

        IF (datatype /= mpitype) CALL errr(__FILE__, __LINE__)

        DO i = 1, length
            iodata(i)%sumTime = iodata(i)%sumTime + idata(i)%sumTime
            iodata(i)%minTime = MIN(iodata(i)%minTime, idata(i)%minTime)
            iodata(i)%maxTime = MAX(iodata(i)%maxTime, idata(i)%maxTime)
            iodata(i)%n = iodata(i)%n + idata(i)%n

            IF (iodata(i)%n > 0) THEN
                iodata(i)%avgTime = &
                iodata(i)%sumTime/REAL(iodata(i)%n, real64)
            ELSE
                iodata(i)%avgTime = 0.0
            END IF

            iodata(i)%avgN = &
                REAL(iodata(i)%n, real64)/REAL(numprocs, real64)
        END DO
    END SUBROUTINE reduce_time

END MODULE timer_mod
