MODULE fort7_mod
    USE config_mod, ONLY: config_t
    USE precision_mod, ONLY: mglet_filename_max
    USE comms_mod, ONLY: myid
    USE err_mod

    IMPLICIT NONE(type, external)
    PRIVATE

    TYPE(config_t) :: fort7

    ! Some "global" variables that are very useful to have in a common
    ! place for all plugins and modules to use
    LOGICAL, PROTECTED :: dread
    LOGICAL, PROTECTED :: dwrite
    LOGICAL, PROTECTED :: dcont

    CHARACTER(len=mglet_filename_max) :: filename = 'parameters.json'

    PUBLIC :: fort7, dread, dwrite, dcont, init_fort7, finish_fort7

CONTAINS
    SUBROUTINE init_fort7()
        ! Local variables
        LOGICAL :: exists
        INTEGER :: nargs, status

        ! Name of parameter file can optionally be passed on the CLI as the
        ! first argument.
        nargs = COMMAND_ARGUMENT_COUNT()
        SELECT CASE(nargs)
        CASE (0)
            CONTINUE
        CASE (1)
            CALL GET_COMMAND_ARGUMENT(1, value=filename, status=status)
            IF (status /= 0) THEN
                WRITE(*,*) "filename, status: ", filename, status
                CALL errr(__FILE__, __LINE__)
            END IF
            IF (myid == 0) THEN
                WRITE(*, '("Reading parameter file: ", A)') TRIM(filename)
                WRITE(*, '()')
            END IF
        CASE DEFAULT
            CALL errr(__FILE__, __LINE__)
        END SELECT

        INQUIRE(FILE=filename, EXIST=exists)
        IF (.NOT. exists) THEN
            WRITE(*,*) "File not found: ", filename
            CALL errr(__FILE__, __LINE__)
        END IF
        CALL fort7%read(filename)

        CALL fort7%get_value("/time/read", dread, .FALSE.)
        CALL fort7%get_value("/time/write", dwrite, .TRUE.)
        CALL fort7%get_value("/time/continue", dcont, .TRUE.)
        IF (.NOT. dread) THEN
            dcont = .FALSE.
        END IF

        ! is_type = fort7%is_int("/bananas")
        ! WRITE(*,*) "bananas is integer: ", is_type
        !
        ! inumber = 0
        ! CALL fort7%get_value("/bananas", inumber, 10)
        ! WRITE(*,*) "bananans: ", inumber
        !
        ! rnumber = 0.0
        ! CALL fort7%get_value("/apples", rnumber)
        ! WRITE(*,*) "apples: ", rnumber
        !
        ! shortstring = REPEAT(" ", LEN(shortstring))
        ! CALL fort7%get_value("/fruit", shortstring)
        ! WRITE(*,*) "fruit: ", shortstring
        !
        ! is_type = fort7%is_logical("/ufrcon")
        ! WRITE(*,*) "ufrcon is logical: ", is_type
        !
        ! is_type = fort7%is_array("/ufrcon")
        ! WRITE(*,*) "ufrcon is array: ", is_type
        !
        ! CALL fort7%get_size("/ufrcon", arrsize)
        ! WRITE(*,*) "size of ufrcon: ", arrsize
        !
        ! CALL fort7%get_array("/ufrcon", ufrcon)
        ! WRITE(*,*) "ufrcon: ", ufrcon
        !
        ! is_present = fort7%exists("/gradpx")
        ! WRITE(*,*) "gradpx is present: ", is_present
        !
        ! acoustics = fort7%get("/acoustics")
        ! CALL acoustics%get_value("/c", rnumber)
        ! WRITE(*,*) "speed of sound: ", rnumber
        !
        ! WRITE(*,*) "Finishing init_fort.7"
        ! WRITE(*,*) ""
    END SUBROUTINE init_fort7


    SUBROUTINE finish_fort7()
        CALL fort7%finish()
    END SUBROUTINE finish_fort7
END MODULE fort7_mod
