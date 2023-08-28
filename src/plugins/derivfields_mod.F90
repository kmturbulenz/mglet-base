MODULE derivfields_mod
    USE core_mod
    USE MPI_f08
    USE ib_mod, ONLY: ib

    IMPLICIT NONE(type, external)
    PRIVATE

    INTEGER(intk) :: nfields
    CHARACTER(len=nchar_name), ALLOCATABLE :: derivfields(:)

    PUBLIC :: init_derivfields, calc_derivfields, finish_derivfields

CONTAINS
    SUBROUTINE init_derivfields(ittot, mtstep, itint, timeph, dt, tend)
        ! Subroutine arguments
        INTEGER(intk), INTENT(in) :: ittot
        INTEGER(intk), INTENT(in) :: mtstep
        INTEGER(intk), INTENT(in) :: itint
        REAL(realk), INTENT(in) :: timeph
        REAL(realk), INTENT(in) :: dt
        REAL(realk), INTENT(in) :: tend

        ! Local variables
        CHARACTER(len=64) :: jsonptr
        INTEGER(intk) :: i

        nfields = 0
        IF (.NOT. fort7%exists("/derivfields")) THEN
            RETURN
        END IF

        CALL fort7%get_size("/derivfields", nfields)
        IF (nfields <= 0) RETURN

        ALLOCATE(derivfields(nfields))

        DO i = 1, nfields
            WRITE(jsonptr, '("/derivfields/", I0)') i-1
            CALL fort7%get_value(jsonptr, derivfields(i))
            CALL set_field(derivfields(i), dwrite=.TRUE.)
            CALL register_statfield(TRIM(derivfields(i))//"_AVG", comp_avg)
            CALL register_statfield( &
                TRIM(derivfields(i))//TRIM(derivfields(i))//"_AVG", &
                comp_sqr_avg)
        END DO
    END SUBROUTINE init_derivfields


    SUBROUTINE finish_derivfields()
        IF (nfields <= 0) RETURN

        DEALLOCATE(derivfields)
    END SUBROUTINE finish_derivfields


    SUBROUTINE calc_derivfields(itstep, ittot, timeph, dt)
        ! Subroutine arguments
        INTEGER(intk), INTENT(in) :: itstep
        INTEGER(intk), INTENT(in) :: ittot
        REAL(realk), INTENT(in) :: timeph
        REAL(realk), INTENT(in) :: dt

        ! Local variables
        INTEGER(intk) :: i
        TYPE(field_t), POINTER :: field

        IF (nfields <= 0) RETURN

        DO i = 1, nfields
            CALL get_field(field, derivfields(i))

            SELECT CASE (TRIM(derivfields(i)))
            CASE ("DIV")
                BLOCK
                    TYPE(field_t), POINTER :: u, v, w
                    CALL get_field(u, "U")
                    CALL get_field(v, "V")
                    CALL get_field(w, "W")
                    CALL ib%divcal(field, u, v, w, 1.0_realk)
                END BLOCK
            CASE DEFAULT
                CALL errr(__FILE__, __LINE__)
            END SELECT
        END DO
    END SUBROUTINE calc_derivfields
END MODULE derivfields_mod
