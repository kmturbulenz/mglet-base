
MODULE get_grid_mod

    ! Module provides simplier interfaces to get grid pointers
    ! for intfield_t and realfield_t. This measure is supposed to
    ! reduce complexity for offloading compilers

    USE err_mod, ONLY: errr
    USE grids_mod, ONLY: get_mgdims, mygrids, nmygrids, level, get_imygrid
    USE precision_mod
    USE basefield_mod
    USE intfield_mod
    USE realfield_mod

    IMPLICIT NONE(type, external)

    PUBLIC :: get_grid1_real, get_grid3_real, get_grid3_real_linear, &
        get_grid3_ifk, get_grid3_ifk_linear

CONTAINS

    SUBROUTINE get_grid1_real(ptr, field, igrid)
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


    ! >>> Integer fields (ifk, not intk)

    SUBROUTINE get_grid3_ifk(ptr, field, igrid)
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

END MODULE get_grid_mod