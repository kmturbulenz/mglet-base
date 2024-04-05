MODULE envvars_mod
    USE MPI_f08

    USE precision_mod
    USE err_mod
    USE charfunc_mod
    USE comms_mod, ONLY: myid

    IMPLICIT NONE (type, external)
    PRIVATE

#ifdef _MGLET_ENVPREFIX_
    CHARACTER(len=*), PARAMETER :: envprefix = _MGLET_ENVPREFIX_
#else
    CHARACTER(len=*), PARAMETER :: envprefix = ""
#endif

    INTEGER(intk), PARAMETER :: maxlength = 64

    ! Public data items
    PUBLIC :: getenv_char, getenv_char_coll, getenv_int, getenv_int_coll, &
        getenv_bool, getenv_bool_coll

CONTAINS
    SUBROUTINE getenv_char(envvalue, varname, defaultval, found)
        ! Subroutine arguments
        CHARACTER(len=*), INTENT(OUT) :: envvalue
        CHARACTER(len=*), INTENT(IN) :: varname
        CHARACTER(len=*), INTENT(IN), OPTIONAL :: defaultval
        LOGICAL, INTENT(OUT), OPTIONAL :: found

        ! Local variables
        CHARACTER(len=maxlength) :: prefixname
        INTEGER :: length, status
        LOGICAL :: found2

        found2 = .FALSE.

        IF (LEN(envprefix) + LEN(varname) > maxlength) THEN
            WRITE(*,*) "Too long variable name: ", varname
            CALL errr(__FILE__, __LINE__)
        END IF
        prefixname = envprefix//varname

        CALL get_environment_variable(prefixname, envvalue, length, status)
        IF (status == 0) THEN
            found2 = .TRUE.
        ELSE IF (status == 1 .AND. PRESENT(defaultval)) THEN
            ! Variable not present - but we use default
            found2 = .FALSE.
            envvalue = defaultval
        ELSE IF (status == 1 .AND. PRESENT(found)) THEN
            ! Variable not present - but found is passed so we know
            found2 = .FALSE.
            envvalue = REPEAT(" ", LEN(envvalue))
        ELSE IF (status == -1) THEN
            ! Variable contents did not fit into "envvalue"
            found2 = .TRUE.
            WRITE(*,*) "Error reading: ", varname
            WRITE(*,*) "Has length: ", length
            CALL errr(__FILE__, __LINE__)
        ELSE
            WRITE(*,*) "Unknown error, got status: ", status
            CALL errr(__FILE__, __LINE__)
        END IF

        ! Return the found status if required
        IF (PRESENT(found)) THEN
            found = found2
        END IF
    END SUBROUTINE getenv_char


    SUBROUTINE getenv_char_coll(envvalue, varname, defaultval, found)
        ! Gets a character environment variable and broadcast it from rank 0
        ! to all other ranks

        ! Subroutine arguments
        CHARACTER(len=*), INTENT(OUT) :: envvalue
        CHARACTER(len=*), INTENT(IN) :: varname
        CHARACTER(len=*), INTENT(IN), OPTIONAL :: defaultval
        LOGICAL, INTENT(OUT), OPTIONAL :: found

        ! Local variables
        INTEGER(int32) :: nchars

        IF (myid == 0) THEN
            CALL getenv_char(envvalue, varname, defaultval, found)
        END IF

        nchars = LEN(envvalue)
        CALL MPI_Bcast(envvalue, nchars, MPI_CHARACTER, 0, MPI_COMM_WORLD)

        IF (PRESENT(found)) THEN
            CALL MPI_Bcast(found, 1, MPI_LOGICAL, 0, MPI_COMM_WORLD)
        END IF
    END SUBROUTINE getenv_char_coll


    SUBROUTINE getenv_int(envvalue, varname, defaultval)
        ! Subroutine arguments
        INTEGER(intk), INTENT(OUT) :: envvalue
        CHARACTER(len=*), INTENT(IN) :: varname
        INTEGER(intk), INTENT(IN), OPTIONAL :: defaultval

        ! Local variables
        CHARACTER(len=maxlength) :: charval
        CHARACTER(len=maxlength) :: defaultchar

        IF (PRESENT(defaultval)) THEN
            WRITE(defaultchar, '(I0)') defaultval
            CALL getenv_char(charval, varname, defaultchar)
        ELSE
            CALL getenv_char(charval, varname)
        END IF

        READ(charval, '(I11)') envvalue
    END SUBROUTINE getenv_int


    SUBROUTINE getenv_int_coll(envvalue, varname, defaultval)
        ! Subroutine arguments
        INTEGER(intk), INTENT(OUT) :: envvalue
        CHARACTER(len=*), INTENT(IN) :: varname
        INTEGER(intk), INTENT(IN), OPTIONAL :: defaultval

        ! Local variables
        ! none...

        IF (myid == 0) THEN
            CALL getenv_int(envvalue, varname, defaultval)
        END IF
        CALL MPI_Bcast(envvalue, 1, mglet_mpi_int, 0, MPI_COMM_WORLD)
    END SUBROUTINE getenv_int_coll


    SUBROUTINE getenv_bool(envvalue, varname, defaultval)
        ! Subroutine arguments
        LOGICAL, INTENT(OUT) :: envvalue
        CHARACTER(len=*), INTENT(IN) :: varname
        LOGICAL, INTENT(IN), OPTIONAL :: defaultval

        ! Local variables
        CHARACTER(len=16) :: charval
        CHARACTER(len=16) :: upperval
        CHARACTER(len=16) :: defaultchar

        IF (PRESENT(defaultval)) THEN
            IF (defaultval .EQV. .TRUE.) THEN
                defaultchar = "TRUE"
            ELSE
                defaultchar = "FALSE"
            END IF
            CALL getenv_char(charval, varname, defaultchar)
        ELSE
            CALL getenv_char(charval, varname)
        END IF

        upperval = UPPER(charval)
        IF (TRIM(upperval) == "TRUE" .OR. TRIM(upperval) == "YES" .OR. &
                TRIM(upperval) == "1" .OR. TRIM(upperval) == "T") THEN
            envvalue = .TRUE.
        ELSE IF (TRIM(upperval) == "FALSE" .OR. TRIM(upperval) == "NO" .OR. &
                TRIM(upperval) == "0" .OR. TRIM(upperval) == "F") THEN
            envvalue = .FALSE.
        ELSE
            WRITE(*,*) "Unknown bool, got: ", charval
            CALL errr(__FILE__, __LINE__)
        END IF
    END SUBROUTINE getenv_bool


    SUBROUTINE getenv_bool_coll(envvalue, varname, defaultval)
        ! Subroutine arguments
        LOGICAL, INTENT(OUT) :: envvalue
        CHARACTER(len=*), INTENT(IN) :: varname
        LOGICAL, INTENT(IN), OPTIONAL :: defaultval

        ! Local variables
        ! none...

        CALL getenv_bool(envvalue, varname, defaultval)
        CALL MPI_Bcast(envvalue, 1, MPI_LOGICAL, 0, MPI_COMM_WORLD)
    END SUBROUTINE getenv_bool_coll
END MODULE envvars_mod
