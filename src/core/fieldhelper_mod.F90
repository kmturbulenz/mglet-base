
MODULE fieldhelper_mod

    ! The fieldhelper_mod module provides a set of utility functions
    ! to reduce the complexity inherent to accessing the data arrays
    ! for individiul grids. As derived-type objects are badly handled
    ! by compilers on GPU, a helper infrastructure is created to reduced
    ! access to several attributes etc. at runtime.

    USE err_mod, ONLY: errr
    USE field_mod, ONLY: field_t, intfield_t
    USE grids_mod, ONLY: get_mgdims, mygrids, nmygrids, level, get_imygrid
    USE pointers_mod, ONLY: get_ibb
    USE precision_mod, ONLY: realk, intk, ifk

    IMPLICIT NONE(type, external)
    PRIVATE

    INTEGER(intk), ALLOCATABLE, DIMENSION(:) :: ip1x, ip1y, ip1z, ip3
    INTEGER(intk), ALLOCATABLE, DIMENSION(:) :: len1x, len1y, len1z, len3
    !$omp declare target(ip1x, ip1y, ip1z, ip3, len1x, len1y, len1z, len3)

    INTERFACE set_field_arr
        PROCEDURE set_field_arr_realk
        PROCEDURE set_field_arr_ifk
    END INTERFACE set_field_arr

    PUBLIC :: init_fieldhelper, finish_fieldhelper, &
        get_grid1x_real, get_grid1y_real, get_grid1z_real, &
        get_grid3_real, get_grid3_real_linear, &
        get_grid3_ifk, get_grid3_ifk_linear, get_grid3_buffer, set_field_arr, &
        map_arr_to_device, map_arr_from_device, map_buf_to_device

CONTAINS

    SUBROUTINE init_fieldhelper()

        ! The initialization obtains data from dummy field_t objects
        ! to fill local look-up arrays. The objective is to avoid contact with
        ! the field_t objects within offloaded kernels.

        ! Local variables
        TYPE(field_t) :: dummy_3d
        INTEGER(intk) :: imygrid, ii, jj, kk, igrid

        ! Using pressure and coordiante as dummy
        CALL dummy_3d%init("DUMMY")

        ALLOCATE(ip1x(nmygrids))
        ALLOCATE(ip1y(nmygrids))
        ALLOCATE(ip1z(nmygrids))
        ALLOCATE(ip3(nmygrids))
        ALLOCATE(len1x(nmygrids))
        ALLOCATE(len1y(nmygrids))
        ALLOCATE(len1z(nmygrids))
        ALLOCATE(len3(nmygrids))

        DO imygrid = 1, nmygrids

            ! Getting the grid dimensions for the current grid
            igrid = mygrids(imygrid)
            CALL get_mgdims(kk, jj, ii, igrid)

            len1x(imygrid) = ii
            len1y(imygrid) = jj
            len1z(imygrid) = kk
            len3(imygrid) = ii * jj * kk

            IF (imygrid == 1) THEN
                ip1x(1) = 1
                ip1y(1) = 1
                ip1z(1) = 1
                ip3(1) = 1
            ELSE
                ip1x(imygrid) = ip1x(imygrid-1) + len1x(imygrid-1)
                ip1y(imygrid) = ip1y(imygrid-1) + len1y(imygrid-1)
                ip1z(imygrid) = ip1z(imygrid-1) + len1z(imygrid-1)
                ip3(imygrid) = ip3(imygrid-1) + len3(imygrid-1)
            END IF
        END DO

        !$omp target enter data map(always, to: ip1x, ip1y, ip1z, &
        !$omp&  ip3, len1x, len1y, len1z, len3)

        CALL dummy_3d%finish()

    END SUBROUTINE init_fieldhelper


    SUBROUTINE finish_fieldhelper()

        !$omp target exit data map(delete: ip1x, ip1y, ip1z, &
        !$omp&  ip3, len1x, len1y, len1z, len3)
        DEALLOCATE(ip1x)
        DEALLOCATE(ip1y)
        DEALLOCATE(ip1z)
        DEALLOCATE(ip3)
        DEALLOCATE(len1x)
        DEALLOCATE(len1y)
        DEALLOCATE(len1z)
        DEALLOCATE(len3)

    END SUBROUTINE finish_fieldhelper


    SUBROUTINE get_grid1x_real(ptr, field, igrid)
        !$omp declare target
        ! Subroutine arguments
        REAL(realk), POINTER, CONTIGUOUS, INTENT(out) :: ptr(:)
        TYPE(field_t), INTENT(in), TARGET :: field
        INTEGER(intk), INTENT(in) :: igrid

        ! Local variables
        INTEGER(intk) :: ip, len, imygrid

        CALL get_imygrid(imygrid, igrid)
        ip = ip1x(imygrid)
        len = len1x(imygrid)

#ifndef _MGLET_OFFLOAD_PERFORMANCE_
        IF (len <= 0) CALL errr(__FILE__, __LINE__)
#endif

        ptr(1:len) => field%arr(ip:ip+len-1)
    END SUBROUTINE get_grid1x_real

    SUBROUTINE get_grid1y_real(ptr, field, igrid)
        !$omp declare target
        ! Subroutine arguments
        REAL(realk), POINTER, CONTIGUOUS, INTENT(out) :: ptr(:)
        TYPE(field_t), INTENT(in), TARGET :: field
        INTEGER(intk), INTENT(in) :: igrid

        ! Local variables
        INTEGER(intk) :: ip, len, imygrid

        CALL get_imygrid(imygrid, igrid)
        ip = ip1y(imygrid)
        len = len1y(imygrid)

#ifndef _MGLET_OFFLOAD_PERFORMANCE_
        IF (len <= 0) CALL errr(__FILE__, __LINE__)
#endif

        ptr(1:len) => field%arr(ip:ip+len-1)
    END SUBROUTINE get_grid1y_real

    SUBROUTINE get_grid1z_real(ptr, field, igrid)
        !$omp declare target
        ! Subroutine arguments
        REAL(realk), POINTER, CONTIGUOUS, INTENT(out) :: ptr(:)
        TYPE(field_t), INTENT(in), TARGET :: field
        INTEGER(intk), INTENT(in) :: igrid

        ! Local variables
        INTEGER(intk) :: ip, len, imygrid

        CALL get_imygrid(imygrid, igrid)
        ip = ip1z(imygrid)
        len = len1z(imygrid)

#ifndef _MGLET_OFFLOAD_PERFORMANCE_
        IF (len <= 0) CALL errr(__FILE__, __LINE__)
#endif

        ptr(1:len) => field%arr(ip:ip+len-1)
    END SUBROUTINE get_grid1z_real


    SUBROUTINE get_grid3_real(ptr, field, igrid)
        !$omp declare target
        ! Subroutine arguments
        REAL(realk), POINTER, CONTIGUOUS, INTENT(out) :: ptr(:, :, :)
        TYPE(field_t), INTENT(in), TARGET :: field
        INTEGER(intk), INTENT(in) :: igrid

        ! Local variables
        INTEGER(intk) :: kk, jj, ii, ip, imygrid, len

        CALL get_imygrid(imygrid, igrid)
        CALL get_mgdims(kk, jj, ii, igrid)
        len = len3(imygrid)
        ip = ip3(imygrid)

#ifndef _MGLET_OFFLOAD_PERFORMANCE_
        IF (len <= 0) CALL errr(__FILE__, __LINE__)
        IF (len /= kk*jj*ii) CALL errr(__FILE__, __LINE__)
#endif

        ptr(1:kk, 1:jj, 1:ii) => field%arr(ip:ip+kk*jj*ii-1)
    END SUBROUTINE get_grid3_real


    SUBROUTINE get_grid3_real_linear(ptr, field, igrid)
        !$omp declare target
        ! Subroutine arguments
        REAL(realk), POINTER, CONTIGUOUS, INTENT(out) :: ptr(:)
        TYPE(field_t), INTENT(in), TARGET :: field
        INTEGER(intk), INTENT(in) :: igrid

        ! Local variables
        INTEGER(intk) :: kk, jj, ii, ip, imygrid, len

        CALL get_imygrid(imygrid, igrid)
        CALL get_mgdims(kk, jj, ii, igrid)
        len = len3(imygrid)
        ip = ip3(imygrid)

#ifndef _MGLET_OFFLOAD_PERFORMANCE_
        IF (len <= 0) CALL errr(__FILE__, __LINE__)
        IF (len /= kk*jj*ii) CALL errr(__FILE__, __LINE__)
#endif

        ptr(1:len) => field%arr(ip:ip+len-1)
    END SUBROUTINE get_grid3_real_linear


    SUBROUTINE get_grid3_ifk(ptr, field, igrid)
        !$omp declare target
        ! Subroutine arguments
        INTEGER(ifk), POINTER, CONTIGUOUS, INTENT(out) :: ptr(:, :, :)
        TYPE(intfield_t), INTENT(in), TARGET :: field
        INTEGER(intk), INTENT(in) :: igrid

        ! Local variables
        INTEGER(intk) :: kk, jj, ii, ip, imygrid, len

        CALL get_imygrid(imygrid, igrid)
        CALL get_mgdims(kk, jj, ii, igrid)
        len = len3(imygrid)
        ip = ip3(imygrid)

#ifndef _MGLET_OFFLOAD_PERFORMANCE_
        IF (len <= 0) CALL errr(__FILE__, __LINE__)
        IF (len /= kk*jj*ii) CALL errr(__FILE__, __LINE__)
#endif

        ptr(1:kk, 1:jj, 1:ii) => field%arr(ip:ip+kk*jj*ii-1)
    END SUBROUTINE get_grid3_ifk


    SUBROUTINE get_grid3_ifk_linear(ptr, field, igrid)
        !$omp declare target
        ! Subroutine arguments
        INTEGER(ifk), POINTER, CONTIGUOUS, INTENT(out) :: ptr(:)
        TYPE(intfield_t), INTENT(in), TARGET :: field
        INTEGER(intk), INTENT(in) :: igrid

        ! Local variables
        INTEGER(intk) :: kk, jj, ii, ip, imygrid, len

        CALL get_imygrid(imygrid, igrid)
        CALL get_mgdims(kk, jj, ii, igrid)
        len = len3(imygrid)
        ip = ip3(imygrid)

#ifndef _MGLET_OFFLOAD_PERFORMANCE_
        IF (len <= 0) CALL errr(__FILE__, __LINE__)
        IF (len /= kk*jj*ii) CALL errr(__FILE__, __LINE__)
#endif

        ptr(1:len) => field%arr(ip:ip+len-1)
    END SUBROUTINE get_grid3_ifk_linear


    SUBROUTINE get_grid3_buffer(ptr, field, igrid, iface)
        !$omp declare target
        ! Subroutine arguments
        REAL(realk), POINTER, CONTIGUOUS, INTENT(out) :: ptr(:, :)
        TYPE(field_t), INTENT(in), TARGET :: field
        INTEGER(intk), INTENT(in) :: igrid, iface

        ! Local variables
        INTEGER(intk) :: kk, jj, ii, ibb

#ifndef _MGLET_OFFLOAD_PERFORMANCE_
        IF (.NOT. ALLOCATED(field%buffers)) THEN
            WRITE(*, *) "Buffers not initialized"
            CALL errr(__FILE__, __LINE__)
        END IF
#endif

        CALL get_mgdims(kk, jj, ii, igrid)
        CALL get_ibb(ibb, iface, igrid)

#ifndef _MGLET_OFFLOAD_PERFORMANCE_
        ! Buffers are only allocated on FIX, OP1 and PAR boundaries. If the
        ! returned ibb is zero, this means that get_buffer was called on
        ! another boundary condition which does not have a buffer
        IF (ibb == 0) THEN
            WRITE(*, *) "Buffer not allocated for this boundary condition"
            WRITE(*, *) "  iface: ", iface, " igrid: ", igrid
            CALL errr(__FILE__, __LINE__)
        END IF
#endif

        SELECT CASE (iface)
        CASE (1, 2)
            ptr(1:kk, 1:jj) => field%buffers(ibb:ibb+kk*jj-1)
        CASE (3, 4)
            ptr(1:kk, 1:ii) => field%buffers(ibb:ibb+kk*ii-1)
        CASE (5, 6)
            ptr(1:jj, 1:ii) => field%buffers(ibb:ibb+jj*ii-1)
#ifndef _MGLET_OFFLOAD_PERFORMANCE_
        CASE DEFAULT
            WRITE(*, '("Invalid face: ", I0)') iface
            CALL errr(__FILE__, __LINE__)
#endif
        END SELECT
    END SUBROUTINE get_grid3_buffer


    SUBROUTINE set_field_arr_realk(field, val, device)
        ! Subroutine arguments
        TYPE(field_t), INTENT(inout) :: field
        REAL(realk), INTENT(in) :: val
        LOGICAL, OPTIONAL, INTENT(in) :: device

        ! Local variables
        INTEGER(intk) :: i, n
        LOGICAL :: device2

        IF (PRESENT(device)) THEN
            device2 = device
        ELSE
            device2 = .FALSE.
        END IF

        IF (device2) THEN
            n = SIZE(field%arr)
            !$omp target teams loop
            DO i = 1, n
                field%arr(i) = val
            END DO
            !$omp end target teams loop
        ELSE
            ! Faster than loop on CPU
            field%arr = val
        END IF
    END SUBROUTINE set_field_arr_realk


    SUBROUTINE set_field_arr_ifk(field, val, device)
        ! Subroutine arguments
        TYPE(intfield_t), INTENT(inout) :: field
        INTEGER(ifk), INTENT(in) :: val
        LOGICAL, OPTIONAL, INTENT(in) :: device

        ! Local variables
        INTEGER(intk) :: i, n
        LOGICAL :: device2

        IF (PRESENT(device)) THEN
            device2 = device
        ELSE
            device2 = .FALSE.
        END IF

        IF (device2) THEN
            n = SIZE(field%arr)
            !$omp target teams loop
            DO i = 1, n
                field%arr(i) = val
            END DO
            !$omp end target teams loop
        ELSE
            ! Faster than loop on CPU
            field%arr = val
        END IF
    END SUBROUTINE set_field_arr_ifk


    SUBROUTINE map_arr_to_device(f1, f2, f3, f4, f5, f6, message)
        ! Ugly wrapper while code is only partly offloaded
        USE fieldmapper_mod
        USE profile_tools_mod

        ! Subroutine arguments
        TYPE(field_t), INTENT(in) :: f1
        TYPE(field_t), INTENT(in), OPTIONAL :: f2, f3, f4, f5, f6
        CHARACTER(*), INTENT(in), OPTIONAL :: message

        ! Local variables
        LOGICAL :: has_message

        has_message = PRESENT(message)

        IF (has_message) THEN
        ! Can not wrap the IF statement or the compiler will complain about
        ! an unused variable...
#if defined(_MGLET_PROFILE_ANNOTATIONS_) && defined(_MGLET_OFFLOAD_)
            CALL profile_range_push(message)
#endif
        END IF

        !$omp target update to(mapper(field_t__map_arr): f1)

        IF (PRESENT(f2)) THEN
            !$omp target update to(mapper(field_t__map_arr): f2)
        END IF
        IF (PRESENT(f3)) THEN
            !$omp target update to(mapper(field_t__map_arr): f3)
        END IF
        IF (PRESENT(f4)) THEN
            !$omp target update to(mapper(field_t__map_arr): f4)
        END IF
        IF (PRESENT(f5)) THEN
            !$omp target update to(mapper(field_t__map_arr): f5)
        END IF
        IF (PRESENT(f6)) THEN
            !$omp target update to(mapper(field_t__map_arr): f6)
        END IF

#if defined(_MGLET_PROFILE_ANNOTATIONS_) && defined(_MGLET_OFFLOAD_)
        IF (has_message) THEN
            CALL profile_range_pop()
        END IF
#endif
    END SUBROUTINE map_arr_to_device


    SUBROUTINE map_buf_to_device(f1, message)
        ! Ugly wrapper while code is only partly offloaded
        USE fieldmapper_mod
        USE profile_tools_mod

        ! Subroutine arguments
        TYPE(field_t), INTENT(in) :: f1
        CHARACTER(*), INTENT(in), OPTIONAL :: message

        ! Local variables
        LOGICAL :: has_message

        has_message = PRESENT(message)

        IF (has_message) THEN
        ! Can not wrap the IF statement or the compiler will complain about
        ! an unused variable...
#if defined(_MGLET_PROFILE_ANNOTATIONS_) && defined(_MGLET_OFFLOAD_)
            CALL profile_range_push(message)
#endif
        END IF

        !$omp target update to(mapper(field_t__map_buffers): f1)

#if defined(_MGLET_PROFILE_ANNOTATIONS_) && defined(_MGLET_OFFLOAD_)
        IF (has_message) THEN
            CALL profile_range_pop()
        END IF
#endif
    END SUBROUTINE map_buf_to_device


    SUBROUTINE map_arr_from_device(f1, f2, f3, f4, f5, f6, message)
        ! Ugly wrapper while code is only partly offloaded
        USE fieldmapper_mod
        USE profile_tools_mod

        ! Subroutine arguments
        TYPE(field_t), INTENT(in) :: f1
        TYPE(field_t), INTENT(in), OPTIONAL :: f2, f3, f4, f5, f6
        CHARACTER(*), INTENT(in), OPTIONAL :: message

        ! Local variables
        LOGICAL :: has_message

        has_message = PRESENT(message)

        IF (has_message) THEN
        ! Can not wrap the IF statement or the compiler will complain about
        ! an unused variable...
#if defined(_MGLET_PROFILE_ANNOTATIONS_) && defined(_MGLET_OFFLOAD_)
            CALL profile_range_push(message)
#endif
        END IF

        !$omp target update from(mapper(field_t__map_arr): f1)

        IF (PRESENT(f2)) THEN
            !$omp target update from(mapper(field_t__map_arr): f2)
        END IF
        IF (PRESENT(f3)) THEN
            !$omp target update from(mapper(field_t__map_arr): f3)
        END IF
        IF (PRESENT(f4)) THEN
            !$omp target update from(mapper(field_t__map_arr): f4)
        END IF
        IF (PRESENT(f5)) THEN
            !$omp target update from(mapper(field_t__map_arr): f5)
        END IF
        IF (PRESENT(f6)) THEN
            !$omp target update from(mapper(field_t__map_arr): f6)
        END IF

#if defined(_MGLET_PROFILE_ANNOTATIONS_) && defined(_MGLET_OFFLOAD_)
        IF (has_message) THEN
            CALL profile_range_pop()
        END IF
#endif
    END SUBROUTINE map_arr_from_device

END MODULE fieldhelper_mod
