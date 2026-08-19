MODULE fieldpool_mod
    USE precision_mod
    USE err_mod, ONLY: errr
    USE field_mod, ONLY: basefield_t, field_t, intfield_t, nchar_name

    IMPLICIT NONE(type, external)
    PRIVATE

    INTEGER(intk), PARAMETER :: max_realfields = 16
    INTEGER(intk), PARAMETER :: max_intfields = 16

    TYPE(field_t), TARGET :: realfields(max_realfields)
    TYPE(intfield_t), TARGET :: intfields(max_intfields)

    ! stackptr always points to the next free slot
    INTEGER(intk) :: stackptr_real = 1
    INTEGER(intk) :: stackptr_int = 1

    INTERFACE push_field
        MODULE PROCEDURE :: push_realfield
        MODULE PROCEDURE :: push_intfield
    END INTERFACE push_field

    INTERFACE pop_field
        MODULE PROCEDURE :: pop_realfield
        MODULE PROCEDURE :: pop_intfield
    END INTERFACE pop_field

    PUBLIC :: push_field, pop_field, init_fieldpool, finish_fieldpool
CONTAINS
    SUBROUTINE init_fieldpool()
        CONTINUE
    END SUBROUTINE init_fieldpool


    SUBROUTINE finish_fieldpool()
        ! Subroutine arguments
        ! none...

        ! Local variables
        INTEGER(intk) :: i

        IF (stackptr_real > 1) THEN
            WRITE(*, *) "Stack for realfields is not empty when finishing"
            WRITE(*, *) "Top of stack: ", TRIM(realfields(stackptr_real-1)%name)
            CALL errr(__FILE__, __LINE__)
        END IF

        IF (stackptr_int > 1) THEN
            WRITE(*, *) "Stack for intfields is not empty when finishing"
            WRITE(*, *) "Top of stack: ", TRIM(intfields(stackptr_int-1)%name)
            CALL errr(__FILE__, __LINE__)
        END IF

        DO i = 1, max_realfields
            IF (realfields(i)%is_init) THEN
                !$omp target exit data map(delete: realfields(i)%arr)
                !$omp target exit data map(delete: realfields(i)%buffers)
                CALL realfields(i)%finish()
            END IF
        END DO
        stackptr_real = 1

        DO i = 1, max_intfields
            IF (intfields(i)%is_init) THEN
                !$omp target exit data map(delete: intfields(i)%arr)
                CALL intfields(i)%finish()
            END IF
        END DO
        stackptr_int = 1
    END SUBROUTINE finish_fieldpool


    ! Returns a pointer to a field with allocated arr and buffers that have
    ! unspecified values. The arr member is NOT zeroed by default and the
    ! programmer need to take care to explicitly set the values of the
    ! array before using it.
    SUBROUTINE push_realfield(field, name, istag, jstag, kstag, units)
        ! Subroutine arguments
        TYPE(field_t), INTENT(out), POINTER :: field
        CHARACTER(len=*), INTENT(in) :: name
        INTEGER(intk), INTENT(in), OPTIONAL :: istag, jstag, kstag
        INTEGER(intk), INTENT(in), OPTIONAL :: units(7)

        ! Local variables
        ! none...

        IF (stackptr_real > max_realfields) THEN
            WRITE(*, *) "Exceeded maximum number of realfields in fieldpool."
            WRITE(*, *) "Maximum limit is: ", max_realfields
            WRITE(*, *) "Stack pointer: ", stackptr_real
            CALL errr(__FILE__, __LINE__)
        END IF

        field => realfields(stackptr_real)

        IF (.NOT. field%is_init) THEN
            CALL field%init(name=name, istag=istag, jstag=jstag, kstag=kstag, &
                units=units, zero=.FALSE.)
            CALL field%init_buffers()
            !$omp target enter data map(alloc: realfields(stackptr_real)%arr, &
            !$omp& realfields(stackptr_real)%buffers)
        ELSE
            CALL set_field_properties(field, name, istag, jstag, kstag, units)
        END IF

#ifdef _MGLET_DEBUG_
        ! Fieldpool never set memory to zero. In debug builds, we set the
        ! array to NaN to catch uninitialized memory usage. In this way,
        ! if some code tries to make use of some memory that have not been
        ! initialized, it will be easy to catch the bug.
        BLOCK
            USE, INTRINSIC :: IEEE_ARITHMETIC
            USE, INTRINSIC :: IEEE_EXCEPTIONS

            LOGICAL :: saved_fpe_mode(SIZE(ieee_all))
            REAL(realk) :: nan

            ! Make sure we do not trigger floating point exceptions when
            ! setting the array to NaN
            CALL IEEE_GET_HALTING_MODE(IEEE_ALL, saved_fpe_mode)
            CALL IEEE_SET_HALTING_MODE(IEEE_ALL, .FALSE.)

            ! Define NaN and set that value in the array
            nan = IEEE_VALUE(0.0_realk, IEEE_SIGNALING_NAN)
            field%arr = nan
            !$omp target update to(realfields(stackptr_real)%arr)

            IF (ALLOCATED(field%buffers)) THEN
                field%buffers = nan
                !$omp target update to(realfields(stackptr_real)%buffers)
            END IF

            ! Restore the previous floating point exception mode
            CALL IEEE_SET_FLAG(IEEE_ALL, .FALSE.)
            CALL IEEE_SET_HALTING_MODE(IEEE_ALL, saved_fpe_mode)
        END BLOCK
#endif

        stackptr_real = stackptr_real+1
    END SUBROUTINE push_realfield


    ! Returns a pointer to a field with allocated arr that has unspecified
    ! values.
    ! The arr member is NOT zeroed by default
    SUBROUTINE push_intfield(field, name, istag, jstag, kstag, units)
        ! Subroutine arguments
        TYPE(intfield_t), INTENT(out), POINTER :: field
        CHARACTER(len=*), INTENT(in) :: name
        INTEGER(intk), INTENT(in), OPTIONAL :: istag, jstag, kstag
        INTEGER(intk), INTENT(in), OPTIONAL :: units(7)

        ! Local variables
        ! none...

        IF (stackptr_int > max_intfields) THEN
            WRITE(*, *) "Exceeded maximum number of intfields in fieldpool."
            WRITE(*, *) "Maximum limit is: ", max_intfields
            WRITE(*, *) "Stack pointer: ", stackptr_int
            CALL errr(__FILE__, __LINE__)
        END IF

        field => intfields(stackptr_int)

        IF (.NOT. field%is_init) THEN
            CALL field%init(name=name, istag=istag, jstag=jstag, kstag=kstag, &
                units=units, zero=.FALSE.)
                !$omp target enter data map(alloc: intfields(stackptr_int)%arr)
        ELSE
            CALL set_field_properties(field, name, istag, jstag, kstag, units)
        END IF

        stackptr_int = stackptr_int+1
    END SUBROUTINE push_intfield


    SUBROUTINE pop_realfield(field)
        ! Subroutine arguments
        TYPE(field_t), INTENT(inout), POINTER :: field

        ! Local variables
        ! none...

        IF (stackptr_real <= 1) THEN
            WRITE(*, *) "Attempting to pop from empty realfield stack."
            CALL errr(__FILE__, __LINE__)
        END IF

        IF (.NOT. ASSOCIATED(field, realfields(stackptr_real-1))) THEN
            WRITE(*, *) "Attempting to pop field that is not on top of stack"
            WRITE(*, *) "or unassociated."
            CALL errr(__FILE__, __LINE__)
        END IF

        stackptr_real = stackptr_real - 1
        CALL set_field_properties(field, "UNASSOCIATED")
        NULLIFY(field)
    END SUBROUTINE pop_realfield


    SUBROUTINE pop_intfield(field)
        ! Subroutine arguments
        TYPE(intfield_t), INTENT(inout), POINTER :: field

        ! Local variables
        ! none...

        IF (stackptr_int <= 1) THEN
            WRITE(*, *) "Attempting to pop from empty intfield stack."
            CALL errr(__FILE__, __LINE__)
        END IF

        IF (.NOT. ASSOCIATED(field, intfields(stackptr_int-1))) THEN
            WRITE(*, *) "Attempting to pop field that is not on top of stack"
            WRITE(*, *) "or unassociated."
            CALL errr(__FILE__, __LINE__)
        END IF

        stackptr_int = stackptr_int - 1
        CALL set_field_properties(field, "UNASSOCIATED")
        NULLIFY(field)
    END SUBROUTINE pop_intfield


    SUBROUTINE set_field_properties(field, name, istag, jstag, kstag, units)
        ! Subroutine arguments
        CLASS(basefield_t), INTENT(inout) :: field
        CHARACTER(len=*), INTENT(in) :: name
        INTEGER(intk), INTENT(in), OPTIONAL :: istag, jstag, kstag
        INTEGER(intk), INTENT(in), OPTIONAL :: units(7)

        ! Local variables
        ! none...

        IF (LEN_TRIM(name) > nchar_name) CALL errr(__FILE__, __LINE__)
        field%name = name

        IF (PRESENT(istag)) THEN
            field%istag = istag
        ELSE
            field%istag = 0
        END IF

        IF (PRESENT(jstag)) THEN
            field%jstag = jstag
        ELSE
            field%jstag = 0
        END IF

        IF (PRESENT(kstag)) THEN
            field%kstag = kstag
        ELSE
            field%kstag = 0
        END IF

        IF (PRESENT(units)) THEN
            field%units = units
        ELSE
            field%units = 0
        END IF
    END SUBROUTINE set_field_properties
END MODULE fieldpool_mod
