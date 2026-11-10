MODULE err_mod
    USE, INTRINSIC :: ISO_FORTRAN_ENV, ONLY: error_unit

    USE MPI_f08

    USE precision_mod

    IMPLICIT NONE (type, external)
    PRIVATE

    INTEGER(intk), PARAMETER :: maxcode = 255

    ! Public data items
    PUBLIC :: err_print, err_abort, errr, err_mpi

CONTAINS
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
        CALL MPI_Abort(MPI_COMM_WORLD, mpicode)
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
