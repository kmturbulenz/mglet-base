MODULE statistics_mod
    USE err_mod
    USE field_mod
    USE fields_mod
    USE fort7_mod
    USE precision_mod

    IMPLICIT NONE(type, external)
    PRIVATE

    ABSTRACT INTERFACE
        SUBROUTINE comp_stat_i(field, name, dt)
            IMPORT :: field_t, realk
            TYPE(field_t), INTENT(inout) :: field
            CHARACTER(len=*), INTENT(in) :: name
            REAL(realk), INTENT(in) :: dt
        END SUBROUTINE comp_stat_i
    END INTERFACE

    ! For storing available statfields
    TYPE :: statfield_t
        CHARACTER(len=nchar_name) :: name = REPEAT(" ", nchar_name)
        PROCEDURE(comp_stat_i), POINTER, NOPASS :: func => NULL()
    END TYPE statfield_t

    ! This list store all available statistical fields that can be produced
    INTEGER(intk) :: n_statfields = 0
    TYPE(statfield_t) :: statfields(1000)

    ! Flag to see if init_statistics have been called - after this no fields
    ! can be registered
    LOGICAL :: is_init = .FALSE.

    ! These lists store the actual statistical fields that are being computed
    INTEGER(intk), ALLOCATABLE :: active_fields(:)

    PUBLIC :: init_statistics, finish_statistics, register_statfield, &
        sample_statistics, comp_avg, comp_sqr_avg

CONTAINS
    SUBROUTINE init_statistics()
        ! Subroutine arguments
        ! none...

        ! Local variables
        INTEGER(intk) :: nfields, i, idx
        CHARACTER(len=nchar_name) :: statname
        CHARACTER(len=64) :: jsonptr
        TYPE(field_t) :: field               ! Field to be averaged
        TYPE(field_t), POINTER :: statfield  ! Field to store average in
        REAL(realk) :: tsamp

        is_init = .TRUE.
        IF (.NOT. fort7%is_array("/statistics")) RETURN
        CALL fort7%get_size("/statistics", nfields)
        ALLOCATE(active_fields(nfields))

        ! Create list of fields that are active
        DO i = 1, nfields
            WRITE(jsonptr, '("/statistics/", I0)') i-1
            CALL fort7%get_value(jsonptr, statname)
            CALL get_statidx(idx, statname)
            active_fields(i) = idx
        END DO

        ! Read in fields when restarting, allocate memory
        DO i = 1, nfields
            idx = active_fields(i)

            ! Initialize statfield by creating a field and copy the parameters
            CALL statfields(idx)%func(field, statfields(idx)%name, 1.0_realk)
            CALL set_field(statfields(idx)%name, istag=field%istag, &
                jstag=field%jstag, kstag=field%kstag, units=field%units, &
                dread=dcont, dwrite=dwrite, active_level=field%active_level)
            CALL field%finish()

            ! Initialize TSAMP when not reading in
            IF (.NOT. dcont) THEN
                CALL get_field(statfield, statfields(idx)%name)
                tsamp = 0.0
                CALL statfield%set_attr(tsamp, "TSAMP")
            END IF
        END DO
    END SUBROUTINE init_statistics


    SUBROUTINE finish_statistics()
        IF (ALLOCATED(active_fields)) DEALLOCATE(active_fields)
    END SUBROUTINE finish_statistics


    SUBROUTINE sample_statistics(dt)
        ! Subroutine arguments
        REAL(realk), INTENT(in) :: dt

        ! Local variables
        INTEGER(intk) :: i, idx
        TYPE(field_t) :: field               ! Field to be averaged
        TYPE(field_t), POINTER :: statfield  ! Field to store average in
        REAL(real64) :: fac1, fac2, tsum
        REAL(realk) :: tsamp

        IF (.NOT. ALLOCATED(active_fields)) RETURN

        DO i = 1, SIZE(active_fields)
            idx = active_fields(i)

            ! This possibly compute and return the quantity to be averaged
            CALL statfields(idx)%func(field, statfields(idx)%name, dt)

            ! This retrieve the current average field (the one that will
            ! be updated)
            CALL get_field(statfield, statfields(idx)%name)

            ! Get sampling time and compute weighting factors
            CALL statfield%get_attr(tsamp, "TSAMP")
            tsum = REAL(tsamp, real64) + REAL(dt, real64)
            fac1 = REAL(tsamp, real64)/tsum
            fac2 = REAL(dt, real64)/tsum

            ! Compute updated statistical field
            statfield%arr = REAL(fac1, realk)*statfield%arr &
                + REAL(fac2, realk)*field%arr

            ! Update sampling time
            tsamp = REAL(tsum, realk)
            CALL statfield%set_attr(tsamp, "TSAMP")

            ! Prepare for next field
            CALL field%finish()
        END DO
    END SUBROUTINE sample_statistics


    SUBROUTINE get_statidx(idx, statname)
        ! Subroutine arguments
        INTEGER(intk), INTENT(out) :: idx
        CHARACTER(len=*), INTENT(in) :: statname

        ! Local variables
        INTEGER(intk) :: i

        idx = 0
        DO i = 1, n_statfields
            IF (TRIM(statfields(i)%name) == TRIM(statname)) THEN
                idx = i
                RETURN
            END IF
        END DO

        WRITE(*,*) "Could not find statfield: ", TRIM(statname)
        CALL errr(__FILE__, __LINE__)
    END SUBROUTINE get_statidx


    SUBROUTINE register_statfield(name, funcptr)
        ! Subroutine arguments
        CHARACTER(len=*), INTENT(in) :: name
        PROCEDURE(comp_stat_i), OPTIONAL :: funcptr

        IF (is_init) ERROR STOP
        IF (n_statfields >= SIZE(statfields)) ERROR STOP
        IF (LEN_TRIM(name) > nchar_name) ERROR STOP
        n_statfields = n_statfields + 1
        statfields(n_statfields)%name = TRIM(name)
        statfields(n_statfields)%func => funcptr
    END SUBROUTINE register_statfield


    ! General routine to compute "ordinary averages" such as U_AVG, V_AVG etc.
    ! Does not do any interpolation of staggered quantities.
    !
    ! Use the naming convention U_AVG -> average of U etc.
    SUBROUTINE comp_avg(field, name, dt)
        ! Subroutine arguments
        TYPE(field_t), INTENT(inout) :: field
        CHARACTER(len=*), INTENT(in) :: name
        REAL(realk), INTENT(in) :: dt

        ! Local variables
        TYPE(field_t), POINTER :: infield
        CHARACTER(len=nchar_name) :: base_name
        INTEGER(intk) :: nchar

        ! Strip off "_AVG" at end of name to get field to compute average from
        nchar = LEN_TRIM(name)

        ! Sanity checks
        IF (nchar < 5) CALL errr(__FILE__, __LINE__)
        IF (nchar > nchar_name) CALL errr(__FILE__, __LINE__)
        IF (name(nchar-3:nchar) /= "_AVG") CALL errr(__FILE__, __LINE__)

        base_name = name(1:nchar-4)
        CALL get_field(infield, base_name)

        CALL field%copy_from(infield)
    END SUBROUTINE comp_avg


    ! General routine to compute squares such as UU_AVG, VV_AVG etc.
    ! Does not do any interpolation of staggered quantities.
    !
    ! Use the naming convention UU_AVG -> average of U*U etc.
    SUBROUTINE comp_sqr_avg(field, name, dt)
        ! Subroutine arguments
        TYPE(field_t), INTENT(inout) :: field
        CHARACTER(len=*), INTENT(in) :: name
        REAL(realk), INTENT(in) :: dt

        ! Local variables
        TYPE(field_t), POINTER :: infield
        CHARACTER(len=nchar_name) :: base_name
        INTEGER(intk) :: nchar

        ! Strip off "_AVG" at end of name to get field to compute average from
        nchar = LEN_TRIM(name)

        ! Sanity checks
        IF (nchar < 6) CALL errr(__FILE__, __LINE__)
        IF (nchar > nchar_name) CALL errr(__FILE__, __LINE__)
        IF (name(nchar-3:nchar) /= "_AVG") CALL errr(__FILE__, __LINE__)
        IF (MOD(nchar-4, 2) /= 0) CALL errr(__FILE__, __LINE__)

        base_name = name(1:(nchar-4)/2)
        CALL get_field(infield, base_name)

        CALL field%copy_from(infield)
        field%units = infield%units*2
        field%arr = field%arr(:)**2
    END SUBROUTINE comp_sqr_avg

END MODULE statistics_mod
