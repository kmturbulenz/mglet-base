
MODULE fieldhelper_mod
    USE err_mod, ONLY: errr
    USE field_mod, ONLY: field_t, intfield_t
    USE grids_mod, ONLY: get_mgdims, mygrids, nmygrids, level, get_imygrid
    USE pointers_mod, ONLY: get_ibb
    USE precision_mod, ONLY: realk, intk, ifk

    IMPLICIT NONE(type, external)
    PRIVATE

    INTERFACE set_field_arr
        PROCEDURE set_field_arr_realk
        PROCEDURE set_field_arr_ifk
    END INTERFACE set_field_arr

    PUBLIC :: get_grid1_real, get_grid3_real, get_grid3_real_linear, &
        get_grid3_ifk, get_grid3_ifk_linear, get_grid3_buffer, set_field_arr, &
        map_arr_to_device, map_arr_from_device, map_buf_to_device
CONTAINS
    SUBROUTINE get_grid1_real(ptr, field, igrid)
        !$omp declare target
        ! Subroutine arguments
        REAL(realk), POINTER, CONTIGUOUS, INTENT(out) :: ptr(:)
        TYPE(field_t), INTENT(in), TARGET :: field
        INTEGER(intk), INTENT(in) :: igrid

        ! Local variables
        INTEGER(intk) :: ip, len, imygrid

        CALL get_imygrid(imygrid, igrid)
        ip = field%ptr(imygrid)
        len = field%length(imygrid)
#ifndef _MGLET_OFFLOAD_PERFORMANCE_
        IF (len <= 0) CALL errr(__FILE__, __LINE__)
#endif

        ptr(1:len) => field%arr(ip:ip+len-1)
    END SUBROUTINE get_grid1_real


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
        len = field%length(imygrid)
#ifndef _MGLET_OFFLOAD_PERFORMANCE_
        IF (len <= 0) CALL errr(__FILE__, __LINE__)
        IF (len /= kk*jj*ii) CALL errr(__FILE__, __LINE__)
#endif

        ip = field%ptr(imygrid)
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
        len = field%length(imygrid)
#ifndef _MGLET_OFFLOAD_PERFORMANCE_
        IF (len <= 0) CALL errr(__FILE__, __LINE__)
        IF (len /= kk*jj*ii) CALL errr(__FILE__, __LINE__)
#endif

        ip = field%ptr(imygrid)
        ptr(1:kk*jj*ii) => field%arr(ip:ip+kk*jj*ii-1)
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
        len = field%length(imygrid)
#ifndef _MGLET_OFFLOAD_PERFORMANCE_
        IF (len <= 0) CALL errr(__FILE__, __LINE__)
        IF (len /= kk*jj*ii) CALL errr(__FILE__, __LINE__)
#endif

        ip = field%ptr(imygrid)
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
        len = field%length(imygrid)
#ifndef _MGLET_OFFLOAD_PERFORMANCE_
        IF (len <= 0) CALL errr(__FILE__, __LINE__)
        IF (len /= kk*jj*ii) CALL errr(__FILE__, __LINE__)
#endif

        ip = field%ptr(imygrid)
        ptr(1:kk*jj*ii) => field%arr(ip:ip+kk*jj*ii-1)
    END SUBROUTINE get_grid3_ifk_linear


    SUBROUTINE get_grid3_buffer(ptr, field, igrid, iface)
        !$omp declare target
        ! Subroutine arguments
        REAL(realk), POINTER, CONTIGUOUS, INTENT(out) :: ptr(:, :)
        TYPE(field_t), INTENT(in), TARGET :: field
        INTEGER(intk), INTENT(in) :: igrid, iface

        ! Local variables
        INTEGER(intk) :: kk, jj, ii, ibb

        IF (.NOT. ALLOCATED(field%buffers)) THEN
            WRITE(*, *) "Buffers not initialized"
            CALL errr(__FILE__, __LINE__)
        END IF

        CALL get_mgdims(kk, jj, ii, igrid)
        CALL get_ibb(ibb, iface, igrid)

        ! Buffers are only allocated on FIX, OP1 and PAR boundaries. If the
        ! returned ibb is zero, this means that get_buffer was called on
        ! another boundary condition which does not have a buffer
        IF (ibb == 0) THEN
            WRITE(*, *) "Buffer not allocated for this boundary condition"
            WRITE(*, *) "  iface: ", iface, " igrid: ", igrid
            CALL errr(__FILE__, __LINE__)
        END IF

        SELECT CASE (iface)
        CASE (1, 2)
            ptr(1:kk, 1:jj) => field%buffers(ibb:ibb+kk*jj-1)
        CASE (3, 4)
            ptr(1:kk, 1:ii) => field%buffers(ibb:ibb+kk*ii-1)
        CASE (5, 6)
            ptr(1:jj, 1:ii) => field%buffers(ibb:ibb+jj*ii-1)
        CASE DEFAULT
            WRITE(*, '("Invalid face: ", I0)') iface
            CALL errr(__FILE__, __LINE__)
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
        USE roctxprofile_mod

        ! Subroutine arguments
        TYPE(field_t), INTENT(in) :: f1
        TYPE(field_t), INTENT(in), OPTIONAL :: f2, f3, f4, f5, f6
        CHARACTER(*), INTENT(in), OPTIONAL :: message

        ! Local variables
        LOGICAL :: has_message

        has_message = PRESENT(message)

        IF (has_message) THEN
            CALL roctxrangepush(message)
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

        IF (has_message) THEN
            CALL roctxrangepop()
        END IF
    END SUBROUTINE map_arr_to_device


    SUBROUTINE map_buf_to_device(f1, message)
        ! Ugly wrapper while code is only partly offloaded
        USE fieldmapper_mod
        USE roctxprofile_mod

        ! Subroutine arguments
        TYPE(field_t), INTENT(in) :: f1
        CHARACTER(*), INTENT(in), OPTIONAL :: message

        ! Local variables
        LOGICAL :: has_message

        has_message = PRESENT(message)

        IF (has_message) THEN
            CALL roctxrangepush(message)
        END IF

        !$omp target update to(mapper(field_t__map_buffers): f1)

        IF (has_message) THEN
            CALL roctxrangepop()
        END IF
    END SUBROUTINE map_buf_to_device


    SUBROUTINE map_arr_from_device(f1, f2, f3, f4, f5, f6, message)
        ! Ugly wrapper while code is only partly offloaded
        USE fieldmapper_mod
        USE roctxprofile_mod

        ! Subroutine arguments
        TYPE(field_t), INTENT(in) :: f1
        TYPE(field_t), INTENT(in), OPTIONAL :: f2, f3, f4, f5, f6
        CHARACTER(*), INTENT(in), OPTIONAL :: message

        ! Local variables
        LOGICAL :: has_message

        has_message = PRESENT(message)

        IF (has_message) THEN
            CALL roctxrangepush(message)
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

        IF (has_message) THEN
            CALL roctxrangepop()
        END IF
    END SUBROUTINE map_arr_from_device
END MODULE fieldhelper_mod
