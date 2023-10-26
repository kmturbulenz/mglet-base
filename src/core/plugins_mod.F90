MODULE plugins_mod
    USE precision_mod, ONLY: intk, realk

    IMPLICIT NONE(type, external)
    PRIVATE

    ABSTRACT INTERFACE
        SUBROUTINE init_plugin_i()
            ! No arguments
        END SUBROUTINE init_plugin_i

        SUBROUTINE finish_plugin_i()
            ! No arguments
        END SUBROUTINE finish_plugin_i

        SUBROUTINE init_late_plugin_i(ittot, mtstep, itint, timeph, dt, tend)
            IMPORT :: intk, realk
            INTEGER(intk), INTENT(in) :: ittot
            INTEGER(intk), INTENT(in) :: mtstep
            INTEGER(intk), INTENT(in) :: itint
            REAL(realk), INTENT(in) :: timeph
            REAL(realk), INTENT(in) :: dt
            REAL(realk), INTENT(in) :: tend
        END SUBROUTINE init_late_plugin_i

        SUBROUTINE timeintegrate_plugin_i(itstep, ittot, timeph, dt)
            IMPORT :: intk, realk
            INTEGER(intk), INTENT(in) :: itstep
            INTEGER(intk), INTENT(in) :: ittot
            REAL(realk), INTENT(in) :: timeph
            REAL(realk), INTENT(in) :: dt
        END SUBROUTINE timeintegrate_plugin_i

        SUBROUTINE can_checkpoint_plugin_i(can_checkpoint)
            LOGICAL, INTENT(out) :: can_checkpoint
        END SUBROUTINE can_checkpoint_plugin_i

        SUBROUTINE checkpoint_plugin_i()
            ! No arguments
        END SUBROUTINE checkpoint_plugin_i
    END INTERFACE

    ! For storing general run-time plugins
    TYPE :: plugin_t
        CHARACTER(len=16) :: name = REPEAT(" ", 16)
        PROCEDURE(init_plugin_i), POINTER, NOPASS :: init => NULL()
        PROCEDURE(init_late_plugin_i), POINTER, NOPASS :: init_late => NULL()
        PROCEDURE(timeintegrate_plugin_i), POINTER, NOPASS :: &
            timeintegrate => NULL()
        PROCEDURE(timeintegrate_plugin_i), POINTER, NOPASS :: &
            itinfo => NULL()
        PROCEDURE(timeintegrate_plugin_i), POINTER, NOPASS :: &
            postprocess => NULL()
        PROCEDURE(can_checkpoint_plugin_i), POINTER, NOPASS :: &
            can_checkpoint => NULL()
        PROCEDURE(checkpoint_plugin_i), POINTER, NOPASS :: checkpoint => NULL()
        PROCEDURE(finish_plugin_i), POINTER, NOPASS :: finish => NULL()
    END TYPE plugin_t
    TYPE(plugin_t) :: plugins(32)
    INTEGER(intk) :: n_plugins = 0

    PUBLIC :: init_plugins, init_late_plugins, timeintegrate_plugins, &
        itinfo_plugins, postprocess_plugins, can_checkpoint_plugins, &
        checkpoint_plugins, finish_plugins, register_plugin

CONTAINS
    SUBROUTINE init_plugins()
        ! Local variables
        INTEGER(intk) :: i

        DO i = 1, n_plugins
            IF (ASSOCIATED(plugins(i)%init)) THEN
                CALL plugins(i)%init()
            END IF
        END DO
    END SUBROUTINE init_plugins


    SUBROUTINE init_late_plugins(ittot, mtstep, itint, timeph, dt, tend)
        ! Subroutine arguments
        INTEGER(intk), INTENT(in) :: ittot
        INTEGER(intk), INTENT(in) :: mtstep
        INTEGER(intk), INTENT(in) :: itint
        REAL(realk), INTENT(in) :: timeph
        REAL(realk), INTENT(in) :: dt
        REAL(realk), INTENT(in) :: tend

        ! Local variables
        INTEGER(intk) :: i

        DO i = 1, n_plugins
            IF (ASSOCIATED(plugins(i)%init_late)) THEN
                CALL plugins(i)%init_late(ittot, mtstep, itint, timeph, &
                    dt, tend)
            END IF
        END DO
    END SUBROUTINE init_late_plugins


    SUBROUTINE timeintegrate_plugins(itstep, ittot, timeph, dt)
        ! Subroutine arguments
        INTEGER(intk), INTENT(in) :: itstep
        INTEGER(intk), INTENT(in) :: ittot
        REAL(realk), INTENT(in) :: timeph
        REAL(realk), INTENT(in) :: dt

        ! Local variables
        INTEGER(intk) :: i

        DO i = 1, n_plugins
            IF (ASSOCIATED(plugins(i)%timeintegrate)) THEN
                CALL plugins(i)%timeintegrate(itstep, ittot, timeph, dt)
            END IF
        END DO
    END SUBROUTINE timeintegrate_plugins


    SUBROUTINE itinfo_plugins(itstep, ittot, timeph, dt)
        ! Subroutine arguments
        INTEGER(intk), INTENT(in) :: itstep
        INTEGER(intk), INTENT(in) :: ittot
        REAL(realk), INTENT(in) :: timeph
        REAL(realk), INTENT(in) :: dt

        ! Local variables
        INTEGER(intk) :: i

        DO i = 1, n_plugins
            IF (ASSOCIATED(plugins(i)%itinfo)) THEN
                CALL plugins(i)%itinfo(itstep, ittot, timeph, dt)
            END IF
        END DO
    END SUBROUTINE itinfo_plugins


    SUBROUTINE postprocess_plugins(itstep, ittot, timeph, dt)
        ! Subroutine arguments
        INTEGER(intk), INTENT(in) :: itstep
        INTEGER(intk), INTENT(in) :: ittot
        REAL(realk), INTENT(in) :: timeph
        REAL(realk), INTENT(in) :: dt

        ! Local variables
        INTEGER(intk) :: i

        DO i = 1, n_plugins
            IF (ASSOCIATED(plugins(i)%postprocess)) THEN
                CALL plugins(i)%postprocess(itstep, ittot, timeph, dt)
            END IF
        END DO
    END SUBROUTINE postprocess_plugins


    SUBROUTINE can_checkpoint_plugins(can_checkpoint)
        ! Subroutine arguments
        LOGICAL, INTENT(out) :: can_checkpoint

        ! Local variables
        INTEGER(intk) :: i
        LOGICAL :: can_checkpoint_this

        can_checkpoint = .TRUE.
        DO i = 1, n_plugins
            IF (ASSOCIATED(plugins(i)%checkpoint)) THEN
                CALL plugins(i)%can_checkpoint(can_checkpoint_this)
                IF (.NOT. can_checkpoint_this) can_checkpoint = .FALSE.
            END IF
        END DO
    END SUBROUTINE can_checkpoint_plugins


    SUBROUTINE checkpoint_plugins()
        ! Local variables
        INTEGER(intk) :: i

        DO i = 1, n_plugins
            IF (ASSOCIATED(plugins(i)%checkpoint)) THEN
                CALL plugins(i)%checkpoint()
            END IF
        END DO
    END SUBROUTINE checkpoint_plugins


    SUBROUTINE finish_plugins()
        ! Local variables
        INTEGER(intk) :: i

        DO i = 1, n_plugins
            IF (ASSOCIATED(plugins(i)%finish)) THEN
                CALL plugins(i)%finish()
            END IF
        END DO
    END SUBROUTINE finish_plugins


    SUBROUTINE register_plugin(name, init, init_late, timeintegrate, &
            itinfo, postprocess, can_checkpoint, checkpoint, finish)
        ! Subroutine arguments
        CHARACTER(len=*), INTENT(in) :: name
        PROCEDURE(init_plugin_i), OPTIONAL :: init
        PROCEDURE(init_late_plugin_i), OPTIONAL :: init_late
        PROCEDURE(timeintegrate_plugin_i), OPTIONAL :: timeintegrate
        PROCEDURE(timeintegrate_plugin_i), OPTIONAL :: itinfo
        PROCEDURE(timeintegrate_plugin_i), OPTIONAL :: postprocess
        PROCEDURE(can_checkpoint_plugin_i), OPTIONAL :: can_checkpoint
        PROCEDURE(checkpoint_plugin_i), OPTIONAL :: checkpoint
        PROCEDURE(finish_plugin_i), OPTIONAL :: finish

        IF (n_plugins >= SIZE(plugins)) ERROR STOP
        IF (LEN_TRIM(name) > 16) ERROR STOP
        n_plugins = n_plugins + 1
        plugins(n_plugins)%name = TRIM(name)

        IF (PRESENT(init)) THEN
            plugins(n_plugins)%init => init
        END IF

        IF (PRESENT(init_late)) THEN
            plugins(n_plugins)%init_late => init_late
        END IF

        IF (PRESENT(finish)) THEN
            plugins(n_plugins)%finish => finish
        END IF

        IF (PRESENT(timeintegrate)) THEN
            plugins(n_plugins)%timeintegrate => timeintegrate
        END IF

        IF (PRESENT(itinfo)) THEN
            plugins(n_plugins)%itinfo => itinfo
        END IF

        IF (PRESENT(postprocess)) THEN
            plugins(n_plugins)%postprocess => postprocess
        END IF

        IF (PRESENT(can_checkpoint)) THEN
            plugins(n_plugins)%can_checkpoint => can_checkpoint
        END IF

        IF (PRESENT(checkpoint)) THEN
            plugins(n_plugins)%checkpoint => checkpoint
        END IF
    END SUBROUTINE register_plugin
END MODULE plugins_mod
