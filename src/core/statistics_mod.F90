MODULE statistics_mod
    USE err_mod
    USE field_mod
    USE fields_mod
    USE fieldpool_mod, ONLY: pop_field, push_field
    USE fort7_mod
    USE grids_mod
    USE pointers_mod, ONLY: get_ip1x, get_ip1y, get_ip1z, get_ip3
    USE precision_mod
    USE timer_mod

    IMPLICIT NONE(type, external)
    PRIVATE

    ABSTRACT INTERFACE
        SUBROUTINE comp_stat_i(field, name, dt)
            IMPORT :: field_t, realk
            TYPE(field_t), POINTER, INTENT(out) :: field
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
        sample_statistics, sync_statistics_to_host, comp_avg, comp_sqr_avg, &
        comp_cube_avg, differentiate, power_stat_field, product_stat_fields, &
        add_stat_fields

CONTAINS
    SUBROUTINE init_statistics()
        ! Subroutine arguments
        ! none...

        ! Local variables
        INTEGER(intk) :: nfields, i, idx
        CHARACTER(len=nchar_name) :: statname
        CHARACTER(len=64) :: jsonptr
        TYPE(field_t), POINTER :: field      ! Field to be averaged
        TYPE(field_t), POINTER :: statfield  ! Field to store average in
        REAL(realk) :: tsamp

        is_init = .TRUE.
        CALL set_timer(170, "STATISTICS")

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
            CALL pop_field(field)

            ! Initialize TSAMP when not reading in
            IF (.NOT. dcont) THEN
                CALL get_field(statfield, statfields(idx)%name)
                tsamp = 0.0
                CALL statfield%set_attr(tsamp, "TSAMP")
            END IF
        END DO
    END SUBROUTINE init_statistics


    SUBROUTINE finish_statistics()
        CALL sync_statistics_to_host()
        IF (ALLOCATED(active_fields)) DEALLOCATE(active_fields)
    END SUBROUTINE finish_statistics


    SUBROUTINE sample_statistics(dt)
        ! Subroutine arguments
        REAL(realk), INTENT(in) :: dt

        ! Local variables
        INTEGER(intk) :: i, idx
        INTEGER(intk) :: iarr
        TYPE(field_t), POINTER :: field      ! Field to be averaged
        TYPE(field_t), POINTER :: statfield  ! Field to store average in
        REAL(real64) :: fac1, fac2, tsum
        REAL(realk) :: fac1_realk, fac2_realk
        REAL(realk) :: tsamp

        IF (.NOT. ALLOCATED(active_fields)) RETURN
        CALL start_timer(170)

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
            fac1_realk = REAL(fac1, realk)
            fac2_realk = REAL(fac2, realk)

            ! Compute updated statistical field
            ASSOCIATE(stat_arr => statfield%arr, sample_arr => field%arr)
                !$omp target teams loop firstprivate(fac1_realk, fac2_realk)
                DO iarr = 1, SIZE(stat_arr)
                    stat_arr(iarr) = fac1_realk*stat_arr(iarr) &
                        + fac2_realk*sample_arr(iarr)
                END DO
                !$omp end target teams loop
            END ASSOCIATE

            ! Update sampling time
            tsamp = REAL(tsum, realk)
            CALL statfield%set_attr(tsamp, "TSAMP")

            ! Prepare for next field
            CALL pop_field(field)
        END DO

        CALL stop_timer(170)
    END SUBROUTINE sample_statistics


    SUBROUTINE sync_statistics_to_host()
        INTEGER(intk) :: i, idx
        TYPE(field_t), POINTER :: statfield

        IF (.NOT. ALLOCATED(active_fields)) RETURN

        DO i = 1, SIZE(active_fields)
            idx = active_fields(i)
            CALL get_field(statfield, statfields(idx)%name)
            !$omp target update from(statfield%arr)
        END DO
    END SUBROUTINE sync_statistics_to_host


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

        WRITE(*, *) "Could not find statfield: ", TRIM(statname)
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
        TYPE(field_t), POINTER, INTENT(out) :: field
        CHARACTER(len=*), INTENT(in) :: name
        REAL(realk), INTENT(in) :: dt

        ! Local variables
        TYPE(field_t), POINTER :: infield
        CHARACTER(len=nchar_name) :: base_name
        INTEGER(intk) :: iarr, nchar

        ! Strip off "_AVG" at end of name to get field to compute average from
        nchar = LEN_TRIM(name)

        ! Sanity checks
        IF (nchar < 5) CALL errr(__FILE__, __LINE__)
        IF (nchar > nchar_name) CALL errr(__FILE__, __LINE__)
        IF (name(nchar-3:nchar) /= "_AVG") CALL errr(__FILE__, __LINE__)

        base_name = name(1:nchar-4)
        CALL get_field(infield, base_name)
        IF (infield%ndim /= 3 .OR. .NOT. ALL(infield%active_level)) THEN
            CALL errr(__FILE__, __LINE__)
        END IF

        CALL push_field(field, name, istag=infield%istag, &
            jstag=infield%jstag, kstag=infield%kstag, units=infield%units)
        IF (SIZE(field%arr) /= SIZE(infield%arr)) CALL errr(__FILE__, __LINE__)

        ASSOCIATE(out => field%arr, input => infield%arr)
            !$omp target teams loop
            DO iarr = 1, SIZE(out)
                out(iarr) = input(iarr)
            END DO
            !$omp end target teams loop
        END ASSOCIATE
    END SUBROUTINE comp_avg


    ! General routine to compute squares such as UU_AVG, VV_AVG etc.
    ! Does not do any interpolation of staggered quantities.
    !
    ! Use the naming convention UU_AVG -> average of U*U etc.
    SUBROUTINE comp_sqr_avg(field, name, dt)
        ! Subroutine arguments
        TYPE(field_t), POINTER, INTENT(out) :: field
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
        IF (infield%ndim /= 3 .OR. .NOT. ALL(infield%active_level)) THEN
            CALL errr(__FILE__, __LINE__)
        END IF

        CALL push_field(field, name, istag=infield%istag, &
            jstag=infield%jstag, kstag=infield%kstag, units=2*infield%units)
        IF (SIZE(field%arr) /= SIZE(infield%arr)) CALL errr(__FILE__, __LINE__)

        CALL product_stat_fields(field, infield, infield)
    END SUBROUTINE comp_sqr_avg


    ! General routine to compute cubes such as UUU_AVG, VVV_AVG etc.
    ! Does not do any interpolation of staggered quantities.
    !
    ! Use the naming convention UUU_AVG -> average of U*U*U etc.
    SUBROUTINE comp_cube_avg(field, name, dt)
        ! Subroutine arguments
        TYPE(field_t), POINTER, INTENT(out) :: field
        CHARACTER(len=*), INTENT(in) :: name
        REAL(realk), INTENT(in) :: dt

        ! Local variables
        TYPE(field_t), POINTER :: infield
        CHARACTER(len=nchar_name) :: base_name
        INTEGER(intk) :: nchar

        ! Strip off "_AVG" at end of name to get field to compute average from
        nchar = LEN_TRIM(name)

        ! Sanity checks
        IF (nchar < 7) CALL errr(__FILE__, __LINE__)
        IF (nchar > nchar_name) CALL errr(__FILE__, __LINE__)
        IF (name(nchar-3:nchar) /= "_AVG") CALL errr(__FILE__, __LINE__)
        IF (MOD(nchar-4, 3) /= 0) CALL errr(__FILE__, __LINE__)

        base_name = name(1:(nchar-4)/3)
        CALL get_field(infield, base_name)
        IF (infield%ndim /= 3 .OR. .NOT. ALL(infield%active_level)) THEN
            CALL errr(__FILE__, __LINE__)
        END IF

        CALL push_field(field, name, istag=infield%istag, &
            jstag=infield%jstag, kstag=infield%kstag, units=3*infield%units)
        IF (SIZE(field%arr) /= SIZE(infield%arr)) CALL errr(__FILE__, __LINE__)

        CALL product_stat_fields(field, infield, infield, infield)
    END SUBROUTINE comp_cube_avg


    SUBROUTINE power_stat_field(field, power)
        TYPE(field_t), INTENT(inout) :: field
        INTEGER(intk), INTENT(in) :: power

        INTEGER(intk) :: iarr

        ASSOCIATE(arr => field%arr)
            !$omp target teams loop firstprivate(power)
            DO iarr = 1, SIZE(arr)
                arr(iarr) = arr(iarr)**power
            END DO
            !$omp end target teams loop
        END ASSOCIATE
    END SUBROUTINE power_stat_field


    SUBROUTINE product_stat_fields(field, a, b, c)
        TYPE(field_t), INTENT(inout) :: field
        TYPE(field_t), INTENT(in) :: a, b
        TYPE(field_t), INTENT(in), OPTIONAL :: c

        INTEGER(intk) :: iarr

        IF (SIZE(field%arr) /= SIZE(a%arr) .OR. &
                SIZE(field%arr) /= SIZE(b%arr)) CALL errr(__FILE__, __LINE__)

        IF (PRESENT(c)) THEN
            IF (SIZE(field%arr) /= SIZE(c%arr)) CALL errr(__FILE__, __LINE__)
            ASSOCIATE(outf => field%arr, in1 => a%arr, in2 => b%arr, &
                    in3 => c%arr)
                !$omp target teams loop
                DO iarr = 1, SIZE(outf)
                    outf(iarr) = in1(iarr)*in2(iarr)*in3(iarr)
                END DO
                !$omp end target teams loop
            END ASSOCIATE
        ELSE
            ASSOCIATE(outf => field%arr, in1 => a%arr, in2 => b%arr)
                !$omp target teams loop
                DO iarr = 1, SIZE(outf)
                    outf(iarr) = in1(iarr)*in2(iarr)
                END DO
                !$omp end target teams loop
            END ASSOCIATE
        END IF
    END SUBROUTINE product_stat_fields


    SUBROUTINE add_stat_fields(field, a, b)
        TYPE(field_t), INTENT(inout) :: field
        TYPE(field_t), INTENT(in) :: a, b

        INTEGER(intk) :: iarr

        IF (SIZE(field%arr) /= SIZE(a%arr) .OR. &
                SIZE(field%arr) /= SIZE(b%arr)) CALL errr(__FILE__, __LINE__)

        ASSOCIATE(outf => field%arr, in1 => a%arr, in2 => b%arr)
            !$omp target teams loop
            DO iarr = 1, SIZE(outf)
                outf(iarr) = in1(iarr) + in2(iarr)
            END DO
            !$omp end target teams loop
        END ASSOCIATE
    END SUBROUTINE add_stat_fields


    ! General routine for spatial differentiation of fields.
    ! Staggeredness of field is not checked internally (!)
    !
    ! naming convention of "ivar":
    ! DD* = diff. in direction of variable's staggering
    ! D*S = diff. in other direction than variable's staggering
    ! D*T = central difference (just for scalar)
    !
    SUBROUTINE differentiate(field, a, ivar)
        ! Subroutine arguments
        CLASS(field_t), INTENT(inout) :: field
        CLASS(field_t), INTENT(in) :: a
        CHARACTER(len=3), INTENT(in) :: ivar

        ! Local Variables
        INTEGER(intk) :: idiff, n
        TYPE(field_t), POINTER :: rdxyz

        ! checking for dimensional consistency
        DO n = 1, 7
            IF (n == 2) THEN
                IF (field%units(n) /= a%units(n)-1) THEN
                    CALL errr(__FILE__, __LINE__)
                END IF
            ELSE
                IF (field%units(n) /= a%units(n)) THEN
                    CALL errr(__FILE__, __LINE__)
                END IF
            END IF
        END DO

        SELECT CASE (ivar)
        CASE ("DXS")
            ! non-x-staggered variable in x
            IF (a%istag == 0 .AND. field%istag == 1 .AND. &
                    a%jstag == field%jstag .AND. &
                    a%kstag == field%kstag) THEN
                CALL get_field(rdxyz, "RDX")
                idiff = 1
            ELSE
                CALL errr(__FILE__, __LINE__)
            END IF
        CASE ("DYS")
            ! non-y-staggered variable in y
            IF (a%istag == field%istag .AND. &
                    a%jstag == 0 .AND. field%jstag == 1 .AND. &
                    a%kstag == field%kstag) THEN
                CALL get_field(rdxyz, "RDY")
                idiff = 2
            ELSE
                CALL errr(__FILE__, __LINE__)
            END IF
        CASE ("DZS")
            ! non-z-staggered variable in z
            IF (a%istag == field%istag .AND. &
                    a%jstag == field%jstag .AND. &
                    a%kstag == 0 .AND. field%kstag == 1) THEN
                CALL get_field(rdxyz, "RDZ")
                idiff = 3
            ELSE
                CALL errr(__FILE__, __LINE__)
            END IF
        CASE ("DXT")
            ! scalar (non-staggered) in x
            IF (a%istag == 0 .AND. field%istag == 0 .AND. &
                    a%jstag == 0 .AND. field%jstag == 0 .AND. &
                    a%kstag == 0 .AND. field%kstag == 0) THEN
                CALL get_field(rdxyz, "DX")
                idiff = 4
            ELSE
                CALL errr(__FILE__, __LINE__)
            END IF
        CASE ("DYT")
            ! scalar (non-staggered) in y
            IF (a%istag == 0 .AND. field%istag == 0 .AND. &
                    a%jstag == 0 .AND. field%jstag == 0 .AND. &
                    a%kstag == 0 .AND. field%kstag == 0) THEN
                CALL get_field(rdxyz, "DY")
                idiff = 5
            ELSE
                CALL errr(__FILE__, __LINE__)
            END IF
        CASE ("DZT")
            ! scalar (non-staggered) in z
            IF (a%istag == 0 .AND. field%istag == 0 .AND. &
                    a%jstag == 0 .AND. field%jstag == 0 .AND. &
                    a%kstag == 0 .AND. field%kstag == 0) THEN
                CALL get_field(rdxyz, "DZ")
                idiff = 6
            ELSE
                CALL errr(__FILE__, __LINE__)
            END IF
        CASE ("DDX")
            ! x-staggered variable in x
            IF (a%istag == 1 .AND. field%istag == 0 .AND. &
                    a%jstag == 0 .AND. field%jstag == 0 .AND. &
                    a%kstag == 0 .AND. field%kstag == 0) THEN
                CALL get_field(rdxyz, "RDDX")
                idiff = 7
            ELSE
                CALL errr(__FILE__, __LINE__)
            END IF
        CASE ("DDY")
            ! y-staggered variable in y
            IF (a%istag == 0 .AND. field%istag == 0 .AND. &
                    a%jstag == 1 .AND. field%jstag == 0 .AND. &
                    a%kstag == 0 .AND. field%kstag == 0) THEN
                CALL get_field(rdxyz, "RDDY")
                idiff = 8
            ELSE
                CALL errr(__FILE__, __LINE__)
            END IF
        CASE ("DDZ")
            ! z-staggered variable in z
            IF (a%istag == 0 .AND. field%istag == 0 .AND. &
                    a%jstag == 0 .AND. field%jstag == 0 .AND. &
                    a%kstag == 1 .AND. field%kstag == 0) THEN
                CALL get_field(rdxyz, "RDDZ")
                idiff = 9
            ELSE
                CALL errr(__FILE__, __LINE__)
            END IF
        CASE DEFAULT
            ! undefined derivative type
            WRITE(*, *) "ivar: ", ivar
            CALL errr(__FILE__, __LINE__)
        END SELECT

        CALL differentiate_impl(field%arr, a%arr, rdxyz%arr, idiff)
    END SUBROUTINE differentiate


    SUBROUTINE differentiate_impl(outf, phi, metric, idiff)
        REAL(realk), INTENT(inout) :: outf(*)
        REAL(realk), INTENT(in) :: phi(*), metric(*)
        INTEGER(intk), INTENT(in) :: idiff

        INTEGER(intk) :: igr, igrid, kk, jj, ii, ip3, ipmetric

        !$omp target teams distribute firstprivate(idiff) &
        !$omp& private(igr, igrid, kk, jj, ii, ip3, ipmetric)
        DO igr = 1, nmygrids
            igrid = mygrids(igr)
            CALL get_mgdims(kk, jj, ii, igrid)
            CALL get_ip3(ip3, igrid)

            SELECT CASE (idiff)
            CASE (1, 4, 7)
                CALL get_ip1x(ipmetric, igrid)
            CASE (2, 5, 8)
                CALL get_ip1y(ipmetric, igrid)
            CASE (3, 6, 9)
                CALL get_ip1z(ipmetric, igrid)
            END SELECT

            !$omp parallel
            CALL differentiate_grid(kk, jj, ii, outf(ip3), phi(ip3), &
                metric(ipmetric), idiff)
            !$omp end parallel
        END DO
        !$omp end target teams distribute
    END SUBROUTINE differentiate_impl


    SUBROUTINE differentiate_grid(kk, jj, ii, outf, phi, metric, idiff)
        !$omp declare target

        INTEGER(intk), INTENT(in) :: kk, jj, ii, idiff
        REAL(realk), INTENT(inout) :: outf(kk, jj, ii)
        REAL(realk), INTENT(in) :: phi(kk, jj, ii), metric(*)

        INTEGER(intk) :: k, j, i

        SELECT CASE (idiff)
        CASE (1)
            !$omp do collapse(3) private(k, j, i)
            DO i = 2, ii-1
                DO j = 2, jj-1
                    DO k = 2, kk-1
                        outf(k, j, i) = (phi(k, j, i+1) - phi(k, j, i)) &
                            *metric(i)
                    END DO
                END DO
            END DO
            !$omp end do
        CASE (2)
            !$omp do collapse(3) private(k, j, i)
            DO i = 2, ii-1
                DO j = 2, jj-1
                    DO k = 2, kk-1
                        outf(k, j, i) = (phi(k, j+1, i) - phi(k, j, i)) &
                            *metric(j)
                    END DO
                END DO
            END DO
            !$omp end do
        CASE (3)
            !$omp do collapse(3) private(k, j, i)
            DO i = 2, ii-1
                DO j = 2, jj-1
                    DO k = 2, kk-1
                        outf(k, j, i) = (phi(k+1, j, i) - phi(k, j, i)) &
                            *metric(k)
                    END DO
                END DO
            END DO
            !$omp end do
        CASE (4)
            !$omp do collapse(3) private(k, j, i)
            DO i = 3, ii-2
                DO j = 3, jj-2
                    DO k = 3, kk-2
                        outf(k, j, i) = (phi(k, j, i+1) - phi(k, j, i-1)) &
                            /(metric(i) + metric(i-1))
                    END DO
                END DO
            END DO
            !$omp end do
        CASE (5)
            !$omp do collapse(3) private(k, j, i)
            DO i = 3, ii-2
                DO j = 3, jj-2
                    DO k = 3, kk-2
                        outf(k, j, i) = (phi(k, j+1, i) - phi(k, j-1, i)) &
                            /(metric(j) + metric(j-1))
                    END DO
                END DO
            END DO
            !$omp end do
        CASE (6)
            !$omp do collapse(3) private(k, j, i)
            DO i = 3, ii-2
                DO j = 3, jj-2
                    DO k = 3, kk-2
                        outf(k, j, i) = (phi(k+1, j, i) - phi(k-1, j, i)) &
                            /(metric(k) + metric(k-1))
                    END DO
                END DO
            END DO
            !$omp end do
        CASE (7)
            !$omp do collapse(3) private(k, j, i)
            DO i = 2, ii-1
                DO j = 2, jj-1
                    DO k = 2, kk-1
                        outf(k, j, i) = (phi(k, j, i) - phi(k, j, i-1)) &
                            *metric(i)
                    END DO
                END DO
            END DO
            !$omp end do
        CASE (8)
            !$omp do collapse(3) private(k, j, i)
            DO i = 2, ii-1
                DO j = 2, jj-1
                    DO k = 2, kk-1
                        outf(k, j, i) = (phi(k, j, i) - phi(k, j-1, i)) &
                            *metric(j)
                    END DO
                END DO
            END DO
            !$omp end do
        CASE (9)
            !$omp do collapse(3) private(k, j, i)
            DO i = 2, ii-1
                DO j = 2, jj-1
                    DO k = 2, kk-1
                        outf(k, j, i) = (phi(k, j, i) - phi(k-1, j, i)) &
                            *metric(k)
                    END DO
                END DO
            END DO
            !$omp end do
        END SELECT
    END SUBROUTINE differentiate_grid

END MODULE statistics_mod
