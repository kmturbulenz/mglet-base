MODULE envvars_mod
    USE MPI_f08

    USE precision_mod
    USE err_mod
    USE charfunc_mod

    IMPLICIT NONE (type, external)
    PRIVATE

#ifdef _MGLET_ENVPREFIX_
    CHARACTER(len=*), PARAMETER :: envprefix = _MGLET_ENVPREFIX_
#else
    CHARACTER(len=*), PARAMETER :: envprefix = "MGLET"
#endif

    INTEGER(intk), PARAMETER :: maxlength = 64

    ! Public data items
    PUBLIC :: getenv_char, getenv_char_coll, getenv_int, getenv_int_coll, &
        getenv_bool, getenv_bool_coll

CONTAINS
    SUBROUTINE getenv_char(varname, envvalue, defaultval)
        ! Subroutine arguments
        CHARACTER(len=*), INTENT(IN) :: varname
        CHARACTER(len=*), INTENT(OUT) :: envvalue
        CHARACTER(len=*), INTENT(IN), OPTIONAL :: defaultval

        ! Local variables
        CHARACTER(len=maxlength) :: prefixname
        INTEGER :: length, status

        IF (LEN(envprefix) + LEN("_") + LEN(varname) > maxlength) THEN
            CALL err_abort(255, "Too long variable name: "//varname, &
                __FILE__, __LINE__)
        END IF
        prefixname = envprefix//"_"//varname

        CALL get_environment_variable(prefixname, envvalue, length, status)
        IF (status == 0) THEN
            ! No error - continue
            CONTINUE
        ELSE IF (status == 1 .AND. PRESENT(defaultval)) THEN
            ! Variable not present - but we use default
            envvalue = defaultval
        ELSE
            CALL err_abort(255, "Unknown error, got status: "//CHAR(status), &
                __FILE__, __LINE__)
        END IF
    END SUBROUTINE getenv_char


    SUBROUTINE getenv_char_coll(varname, value, defaultval)
        ! Gets a character environment variable and broadcast it from rank 0
        ! to all other ranks

        ! Subroutine arguments
        CHARACTER(len=*), INTENT(IN) :: varname
        CHARACTER(len=*), INTENT(OUT) :: value
        CHARACTER(len=*), INTENT(IN), OPTIONAL :: defaultval

        ! Local variables
        INTEGER(int32) :: nchars

        CALL getenv_char(varname, value, defaultval)
        nchars = LEN(value)
        CALL MPI_Bcast(value, nchars, MPI_CHARACTER, 0, MPI_COMM_WORLD)
    END SUBROUTINE getenv_char_coll


    SUBROUTINE getenv_int(varname, value, defaultval)
        ! Subroutine arguments
        CHARACTER(len=*), INTENT(IN) :: varname
        INTEGER(intk), INTENT(OUT) :: value
        INTEGER(intk), INTENT(IN), OPTIONAL :: defaultval

        ! Local variables
        CHARACTER(len=maxlength) :: charval
        CHARACTER(len=maxlength) :: defaultchar

        IF (PRESENT(defaultval)) THEN
            WRITE(defaultchar, '(I0)') defaultval
            CALL getenv_char(varname, charval, defaultchar)
        ELSE
            CALL getenv_char(varname, charval)
        END IF

        READ(charval, '(I11)') value
    END SUBROUTINE getenv_int


    SUBROUTINE getenv_int_coll(varname, value, defaultval)
        ! Subroutine arguments
        CHARACTER(len=*), INTENT(IN) :: varname
        INTEGER(intk), INTENT(OUT) :: value
        INTEGER(intk), INTENT(IN), OPTIONAL :: defaultval

        ! Local variables
        ! none...

        CALL getenv_int(varname, value, defaultval)
        CALL MPI_Bcast(value, 1, mglet_mpi_int, 0, MPI_COMM_WORLD)
    END SUBROUTINE getenv_int_coll


    SUBROUTINE getenv_bool(varname, value, defaultval)
        ! Subroutine arguments
        CHARACTER(len=*), INTENT(IN) :: varname
        LOGICAL, INTENT(OUT) :: value
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
            CALL getenv_char(varname, charval, defaultchar)
        ELSE
            CALL getenv_char(varname, charval)
        END IF

        upperval = UPPER(charval)
        IF (upperval(1:5) == "TRUE " .OR. upperval(1:5) == "YES  " .OR. &
                upperval(1:5) == "1    ") THEN
            value = .TRUE.
        ELSE IF (upperval(1:5) == "FALSE" .OR. upperval(1:5) == "NO   " .OR. &
                upperval(1:5) == "0    ") THEN
            value = .FALSE.
        ELSE
            CALL err_abort(255, "Unknown bool, got: "//charval, &
                __FILE__, __LINE__)
        END IF
    END SUBROUTINE getenv_bool


    SUBROUTINE getenv_bool_coll(varname, value, defaultval)
        ! Subroutine arguments
        CHARACTER(len=*), INTENT(IN) :: varname
        LOGICAL, INTENT(OUT) :: value
        LOGICAL, INTENT(IN), OPTIONAL :: defaultval

        ! Local variables
        ! none...

        CALL getenv_bool(varname, value, defaultval)
        CALL MPI_Bcast(value, 1, MPI_LOGICAL, 0, MPI_COMM_WORLD)
    END SUBROUTINE getenv_bool_coll
END MODULE envvars_mod
