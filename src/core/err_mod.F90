MODULE err_mod
    USE, INTRINSIC :: ISO_FORTRAN_ENV, ONLY: error_unit

    USE MPI_f08

    USE precision_mod

    IMPLICIT NONE (type, external)
    PRIVATE

    INTEGER(intk), PARAMETER :: maxcode = 255

    INTEGER(int32) :: aux_flag = 0_int32
    !$omp declare target(aux_flag)

    INTEGER(int32), PARAMETER :: maxfilelength = 512
    INTEGER(int32), PROTECTED :: err_kernel_line = 0_int32
    CHARACTER(len=maxfilelength), PROTECTED :: err_kernel_file = ""
    !$omp declare target(err_kernel_line, err_kernel_file)

    PUBLIC :: errr, err_print, err_abort, err_mpi, &
        errr_kernel, err_kernel_check, init_err_kernel, finish_err_kernel
CONTAINS
    SUBROUTINE errr_kernel(fname, line)
        !$omp declare target
        ! Kernel error to record error information on the device side

        ! Subroutine arguments
        CHARACTER(len=*), INTENT(in) :: fname
        INTEGER(int32), INTENT(in) :: line

        ! Local variables
        INTEGER(int32) :: access_flag

        ! Set kernel error information (note that "fname" may be truncated)

        !$omp atomic capture
        access_flag = aux_flag
        aux_flag = 1_int32
        !$omp end atomic

        IF (access_flag == 0) THEN
            err_kernel_line = line
            err_kernel_file = fname
        END IF
    END SUBROUTINE errr_kernel


    SUBROUTINE err_kernel_check()
        ! Kernel check error to evaluate the status from the host side

        ! Local variables
        INTEGER(intk), PARAMETER :: default_code = 254
        CHARACTER(len=*), PARAMETER :: default_message = "MGLET Kernel Error"

        !$omp target update from(err_kernel_line)
        IF (err_kernel_line > 0_int32) THEN
            !$omp target update from(err_kernel_file)
            CALL err_abort(default_code, default_message, &
                TRIM(err_kernel_file), INT(err_kernel_line, intk))
        END IF
    END SUBROUTINE err_kernel_check


    SUBROUTINE init_err_kernel()
        ! Initialize kernel error handling on the host side

        err_kernel_line = -1_int32
        err_kernel_file = REPEAT(" ", maxfilelength)
        aux_flag = 0_int32  ! = setting active

        !$omp target update to(err_kernel_line, err_kernel_file, aux_flag)
    END SUBROUTINE init_err_kernel


    SUBROUTINE finish_err_kernel()
        ! Finalize kernel error handling on the host side

        CALL err_kernel_check()

        err_kernel_line = -1_int32
        err_kernel_file = REPEAT(" ", maxfilelength)
        aux_flag = 1_int32  ! = setting inactive

        !$omp target update to(err_kernel_line, err_kernel_file, aux_flag)
    END SUBROUTINE finish_err_kernel


    SUBROUTINE err_print(errorcode, message, fname, line)
        ! Prints error message to standard error. Does not abort execution.
        USE comms_mod, ONLY: myid
#if defined __INTEL_COMPILER
        USE IFCORE
#endif

        ! Subroutine arguments
        INTEGER(intk), INTENT(in) :: errorcode
        CHARACTER(len=*), INTENT(in) :: message
        CHARACTER(len=*), INTENT(in), OPTIONAL :: fname
        INTEGER(intk), INTENT(in), OPTIONAL :: line

        write(error_unit, '()')
        write(error_unit, '("ERRR: ", I0, " MGLET FATAL ERROR")') myid
        write(error_unit, '("ERRR: ", I0, " CODE:    ", I0)') myid, errorcode
        write(error_unit, '("ERRR: ", I0, " MESSAGE: ", A)') myid, message
        IF (PRESENT(fname)) THEN
            write(error_unit, '("ERRR: ", I0, " FILE:    ", A)') myid, fname
        END IF
        IF (PRESENT(line)) THEN
            write(error_unit, '("ERRR: ", I0, " LINE:    ", I0)') myid, line
        END IF
        write(error_unit, '()')

        ! Produce backtrace to indicate from where the error comming from
#ifdef __INTEL_COMPILER
        CALL TRACEBACKQQ(user_exit_code=-1)
#endif
#ifdef __GFORTRAN__
        CALL BACKTRACE()
#endif
    END SUBROUTINE err_print


    SUBROUTINE err_abort(errorcode, message, fname, line)
        ! Generic error handler if something happens on an arbitrary single
        ! rank. Aborts execution of MGLET.

        ! Subroutine arguments
        INTEGER(intk), INTENT(in) :: errorcode
        CHARACTER(len=*), INTENT(in) :: message
        CHARACTER(len=*), INTENT(in), OPTIONAL :: fname
        INTEGER(intk), INTENT(in), OPTIONAL :: line

        ! Local variables
        INTEGER(kind=int32) :: mpicode

        CALL err_print(errorcode, message, fname, line)

        mpicode = INT(errorcode, int32)
        IF (mpicode > maxcode) THEN
            mpicode = maxcode
        END IF
#ifdef _MGLET_OFFLOAD_
        CALL EXIT(mpicode)
#else
        CALL MPI_Abort(MPI_COMM_WORLD, mpicode)
#endif
    END SUBROUTINE err_abort


    SUBROUTINE errr(fname, line)
        ! Generic error handler if something happens on an arbitrary single
        ! rank. Aborts execution of MGLET.

        ! Subroutine arguments
        CHARACTER(len=*), INTENT(in) :: fname
        INTEGER(int32), INTENT(in) :: line

        ! Local variables
        INTEGER(intk), PARAMETER :: default_code = 254
        CHARACTER(len=*), PARAMETER :: default_message = "MGLET Error"

        CALL err_abort(default_code, default_message, fname, INT(line, intk))
    END SUBROUTINE errr


    SUBROUTINE err_mpi(fname, line)
        ! MPI error

        ! Subroutine arguments
        CHARACTER(len=*), INTENT(in), OPTIONAL :: fname
        INTEGER(intk), INTENT(in), OPTIONAL :: line

        CALL err_abort(254_intk, "MPI communication error", fname, line)
    END SUBROUTINE err_mpi
END MODULE err_mod
