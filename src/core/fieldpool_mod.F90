MODULE fieldpool_mod
    USE precision_mod
    USE err_mod, ONLY: errr
    USE field_mod, ONLY: basefield_t, field_t, intfield_t, nchar_name
    USE fieldhelper_mod, ONLY: set_field_arr

    IMPLICIT NONE(type, external)
    PRIVATE

    INTEGER(intk), PARAMETER :: max_realfields = 16
    INTEGER(intk), PARAMETER :: max_intfields = 16

    TYPE(field_t), TARGET :: realfields(max_realfields)
    TYPE(intfield_t), TARGET :: intfields(max_intfields)

    ! stackptr always points to the next free slot
    INTEGER(intk) :: stackptr_real
    INTEGER(intk) :: stackptr_int

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
        ! Subroutine arguments
        ! none...

        ! Local variables
        ! none...

        ! 1-based indexing into field arrays
        stackptr_real = 1
        stackptr_int = 1
    END SUBROUTINE init_fieldpool


    SUBROUTINE finish_fieldpool()
        USE fieldmapper_mod
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
                !$omp target exit data map(mapper(default), &
                !$omp& delete: realfields(i))
                CALL realfields(i)%finish()
            END IF
        END DO
        stackptr_real = 1

        DO i = 1, max_intfields
            IF (intfields(i)%is_init) THEN
                !$omp target exit data map(mapper(default), &
                !$omp& delete: intfields(i))
                CALL intfields(i)%finish()
            END IF
        END DO
        stackptr_int = 1
    END SUBROUTINE finish_fieldpool


    SUBROUTINE push_realfield(field, name, istag, jstag, kstag, units, zero)
        USE fieldmapper_mod
        ! Subroutine arguments
        TYPE(field_t), INTENT(inout), POINTER :: field
        CHARACTER(len=*), INTENT(in) :: name
        INTEGER(intk), INTENT(in), OPTIONAL :: istag, jstag, kstag
        INTEGER(intk), INTENT(in), OPTIONAL :: units(7)
        LOGICAL, INTENT(in), OPTIONAL :: zero

        ! Local variables
        LOGICAL :: zero2, did_init

        IF (PRESENT(zero)) THEN
            zero2 = zero
        ELSE
            zero2 = .TRUE.
        END IF

        IF (stackptr_real > max_realfields) THEN
            WRITE(*, *) "Exceeded maximum number of realfields in fieldpool."
            WRITE(*, *) "Maximum limit is: ", max_realfields
            WRITE(*, *) "Stack pointer: ", stackptr_real
            CALL errr(__FILE__, __LINE__)
        END IF

        field => realfields(stackptr_real)

        did_init = .FALSE.
        IF (.NOT. field%is_init) THEN
            CALL field%init(name=name, istag=istag, jstag=jstag, kstag=kstag, &
                units=units, zero=zero2)
            CALL field%init_buffers()
            !$omp target enter data map(mapper(default), &
            !$omp& to: realfields(stackptr_real))
            did_init = .TRUE.
        ELSE
            CALL set_field_properties(field, name, istag, jstag, kstag, units)
        END IF

        IF (zero2 .AND. .NOT. did_init) THEN
            ! TODO(offload): Once larger portions of the code are ported,
            ! remove the device=.FALSE. setter and preprocessor ifdef.
            ! Currently, we can not generalize whether all fields in the
            ! fieldpool are only required on the device and are thus required
            ! to set both host and device buffers to zero.
#ifdef _MGLET_OFFLOAD_
            CALL set_field_arr(field, 0.0_realk, device=.TRUE.)
#endif
            CALL set_field_arr(field, 0.0_realk, device=.FALSE.)
        END IF

        stackptr_real = stackptr_real+1
    END SUBROUTINE push_realfield


    SUBROUTINE push_intfield(field, name, istag, jstag, kstag, units, zero)
        USE fieldmapper_mod
        ! Subroutine arguments
        TYPE(intfield_t), INTENT(inout), POINTER :: field
        CHARACTER(len=*), INTENT(in) :: name
        INTEGER(intk), INTENT(in), OPTIONAL :: istag, jstag, kstag
        INTEGER(intk), INTENT(in), OPTIONAL :: units(7)
        LOGICAL, INTENT(in), OPTIONAL :: zero

        ! Local variables
        LOGICAL :: zero2, did_init

        IF (PRESENT(zero)) THEN
            zero2 = zero
        ELSE
            zero2 = .TRUE.
        END IF

        IF (stackptr_int > max_intfields) THEN
            WRITE(*, *) "Exceeded maximum number of intfields in fieldpool."
            WRITE(*, *) "Maximum limit is: ", max_intfields
            WRITE(*, *) "Stack pointer: ", stackptr_int
            CALL errr(__FILE__, __LINE__)
        END IF

        field => intfields(stackptr_int)

        did_init = .FALSE.
        IF (.NOT. field%is_init) THEN
            CALL field%init(name=name, istag=istag, jstag=jstag, kstag=kstag, &
                units=units, zero=zero2)
            !$omp target enter data map(mapper(default), &
            !$omp& to: intfields(stackptr_int))
            did_init = .TRUE.
        ELSE
            CALL set_field_properties(field, name, istag, jstag, kstag, units)
        END IF

        IF (zero2 .AND. .NOT. did_init) THEN
            ! TODO(offload): Once larger portions of the code are ported,
            ! remove the device=.FALSE. setter and preprocessor ifdef.
            ! Currently, we can not generalize whether all fields in the
            ! fieldpool are only required on the device and are thus required
            ! to set both host and device buffers to zero.
#ifdef _MGLET_OFFLOAD_
            CALL set_field_arr(field, 0_ifk, device=.TRUE.)
#endif
            CALL set_field_arr(field, 0_ifk, device=.FALSE.)
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
