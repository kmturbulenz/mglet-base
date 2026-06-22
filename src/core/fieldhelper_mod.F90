
MODULE fieldhelper_mod
    USE err_mod, ONLY: errr
    USE field_mod, ONLY: field_t, intfield_t
    USE grids_mod, ONLY: get_mgdims, mygrids, nmygrids, level, get_imygrid
    USE precision_mod, ONLY: realk, intk, ifk

    IMPLICIT NONE(type, external)
    PRIVATE

    INTERFACE set_field_arr
        PROCEDURE set_field_arr_realk
        PROCEDURE set_field_arr_ifk
    END INTERFACE set_field_arr

    PUBLIC :: get_grid1_real, get_grid3_real, get_grid3_real_linear, &
        get_grid3_ifk, get_grid3_ifk_linear, set_field_arr
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

        ! Getting and checking properties
        CALL get_imygrid(imygrid, igrid)
        CALL get_mgdims(kk, jj, ii, igrid)
        len = field%length(imygrid)
        IF (len <= 0) CALL errr(__FILE__, __LINE__)
        IF (len /= kk*jj*ii) CALL errr(__FILE__, __LINE__)

        ! Setting the pointer
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

        ! Getting and checking properties
        CALL get_imygrid(imygrid, igrid)
        CALL get_mgdims(kk, jj, ii, igrid)
        len = field%length(imygrid)
        IF (len <= 0) CALL errr(__FILE__, __LINE__)
        IF (len /= kk*jj*ii) CALL errr(__FILE__, __LINE__)

        ! Setting the pointer
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

        ! Getting and checking properties
        CALL get_imygrid(imygrid, igrid)
        CALL get_mgdims(kk, jj, ii, igrid)
        len = field%length(imygrid)
        IF (len <= 0) CALL errr(__FILE__, __LINE__)
        IF (len /= kk*jj*ii) CALL errr(__FILE__, __LINE__)

        ! Setting the pointer
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

        ! Getting and checking properties
        CALL get_imygrid(imygrid, igrid)
        CALL get_mgdims(kk, jj, ii, igrid)
        len = field%length(imygrid)
        IF (len <= 0) CALL errr(__FILE__, __LINE__)
        IF (len /= kk*jj*ii) CALL errr(__FILE__, __LINE__)

        ! Setting the pointer
        ip = field%ptr(imygrid)
        ptr(1:kk*jj*ii) => field%arr(ip:ip+kk*jj*ii-1)
    END SUBROUTINE get_grid3_ifk_linear


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
END MODULE fieldhelper_mod
