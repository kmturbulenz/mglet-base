MODULE plog_mod
    USE core_mod

    IMPLICIT NONE (type, external)
    PRIVATE

    TYPE, BIND(C) :: plog_t
        INTEGER(c_intk) :: ipc
        REAL(c_realk) :: maxrhs
    END TYPE plog_t

    TYPE(plog_t), ALLOCATABLE :: plog_level(:)

    CHARACTER(len=*), PARAMETER :: logfile = logdir//"/pressuresolver.log"

    PUBLIC :: init_plog, finish_plog, sample_plog, print_plog

CONTAINS
    ! Tag the grids, create MPI datatype and reduction operator
    SUBROUTINE init_plog(dcont)
        ! Subroutine arguments
        LOGICAL, INTENT(in) :: dcont

        ! Local variables
        LOGICAL :: exists
        INTEGER :: logunit

        ALLOCATE(plog_level(minlevel:maxlevel))

        IF (myid == 0) THEN
            INQUIRE(FILE=logfile, EXIST=exists)
            IF (.NOT. exists .OR. .NOT. dcont) THEN
                OPEN(NEWUNIT=logunit, FILE=logfile)
                WRITE(logunit, '(5A9, 1A14)') &
                    "#  ITTOT", "IRK", "IPCOUNT", "LEVEL", "IPC", "MAXRES"
                CLOSE(logunit)
            END IF
        END IF
    END SUBROUTINE init_plog


    SUBROUTINE finish_plog()
        ! Subroutine arguments
        ! None...

        ! Local variables
        ! None...

        IF (ALLOCATED(plog_level)) DEALLOCATE(plog_level)
    END SUBROUTINE finish_plog


    SUBROUTINE sample_plog(ilevel, ipc, maxrhs)
        ! Subroutine arguments
        INTEGER(intk), INTENT(in) :: ilevel
        INTEGER(intk), INTENT(in) :: ipc
        REAL(realk), INTENT(in) :: maxrhs

        plog_level(ilevel)%ipc = ipc
        plog_level(ilevel)%maxrhs = maxrhs
    END SUBROUTINE sample_plog


    SUBROUTINE print_plog(ittot, irk, ipcount)
        ! Subroutine arguments
        INTEGER(intk), INTENT(in) :: ittot
        INTEGER(intk), INTENT(in) :: irk
        INTEGER(intk), INTENT(in) :: ipcount

        ! Local variables
        INTEGER(intk) :: ilevel
        INTEGER :: logunit

        IF (myid == 0) THEN
            OPEN(NEWUNIT=logunit, FILE=logfile, POSITION="APPEND")
            DO ilevel = minlevel, maxlevel
                WRITE(logunit, '(5I9, ES14.5)') &
                    ittot, irk, ipcount, ilevel, plog_level(ilevel)%ipc, &
                    plog_level(ilevel)%maxrhs
            END DO
            CLOSE(logunit)
        END IF
    END SUBROUTINE print_plog
END MODULE plog_mod
