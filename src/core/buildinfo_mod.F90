MODULE buildinfo_mod
    USE comms_mod, ONLY: myid
    USE precision_mod, ONLY: real_bytes, int_bytes, ifk_bytes, intk
    USE envvars_mod, ONLY: getenv_char_coll

    IMPLICIT NONE(type, external)
    PRIVATE

    INTEGER(intk) :: year, month, day, hour, minute, second
    CHARACTER(len=1024) :: mglet_dbg_envvar  ! Exported and globally available

    PUBLIC :: init_buildinfo, finish_buildinfo, get_starttime, mglet_dbg_envvar

CONTAINS
    SUBROUTINE init_buildinfo()
        USE ISO_FORTRAN_ENV, ONLY: COMPILER_VERSION, COMPILER_OPTIONS
        ! Local variables
        INTEGER :: values(8)
        CHARACTER(len=64) :: jobid
        INTEGER :: length, status

        CALL DATE_AND_TIME(VALUES=values)
        year = values(1)
        month = values(2)
        day = values(3)
        hour = values(5)
        minute = values(6)
        second = values(7)

        ! Needs to be called by all processes
        CALL getenv_char_coll(mglet_dbg_envvar, "MGLET_DBG", "")

        ! Non-root proceeses can return here
        IF (myid /= 0) RETURN

        WRITE(*, '("BUILD INFORMATION:")')
        WRITE(*, '("    Date:          ", A)') __DATE__
        WRITE(*, '("    Time:          ", A)') __TIME__
        WRITE(*, '("    Compiler:      ", A)') COMPILER_VERSION()
        WRITE(*, '("    Options:       ", A)') COMPILER_OPTIONS()
        WRITE(*, '()')

        WRITE(*, '("COMPILED PRECISION:")')
        WRITE(*, '("    REAL:          ", I0, " bytes")') real_bytes
        WRITE(*, '("    INTEGER:       ", I0, " bytes")') int_bytes
        WRITE(*, '("    INTEGER(ifk):  ", I0, " bytes")') ifk_bytes
        WRITE(*, '()')

        WRITE(*, '("JOB INFORMATION:")')
        WRITE(*, '("    Date:          ", I0.4, "-", I0.2, "-", I0.2)') &
            year, month, day
        WRITE(*, '("    Time:          ", I0.2, ":", I0.2, ":", I0.2)') &
            hour, minute, second

        ! Try to fetch batch system job ID
        CALL get_environment_variable("PBS_JOBID", jobid, length, status)
        IF (status < 1) THEN
            WRITE(*, '("    PBS_JOBID:     ", A)') TRIM(jobid)
        END IF

        CALL get_environment_variable("LOADL_STEP_ID", jobid, length, status)
        IF (status < 1) THEN
            WRITE(*, '("    LOADL_STEP_ID: ", A)') TRIM(jobid)
        END IF

        CALL get_environment_variable("SLURM_JOBID", jobid, length, status)
        IF (status < 1) THEN
            WRITE(*, '("    SLURM_JOBID:   ", A)') TRIM(jobid)
        END IF

        IF (LEN_TRIM(mglet_dbg_envvar) > 0) THEN
            WRITE(*, '("    MGLET_DBG:     ", A)') TRIM(mglet_dbg_envvar)
        END IF

        WRITE(*, '()')
    END SUBROUTINE init_buildinfo


    SUBROUTINE finish_buildinfo()
        CONTINUE
    END SUBROUTINE finish_buildinfo


    SUBROUTINE get_starttime(values)
        INTEGER(intk), INTENT(OUT) :: values(8)

        values(1) = year
        values(2) = month
        values(3) = day
        values(4) = 0
        values(5) = hour
        values(6) = minute
        values(7) = second
        values(7) = 0
    END SUBROUTINE get_starttime
END MODULE buildinfo_mod
