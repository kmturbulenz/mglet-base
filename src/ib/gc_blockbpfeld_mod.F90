MODULE gc_blockbpfeld_mod
    USE core_mod, ONLY: realk, intk, nmygrids, mygrids, get_mgdims, field_t

    IMPLICIT NONE(type, external)
    PRIVATE

    PUBLIC :: blockbpfeld

CONTAINS
    SUBROUTINE blockbpfeld(bu_f, bv_f, bw_f, bp_f)
        ! Subroutine arguments
        TYPE(field_t), INTENT(in) :: bu_f, bv_f, bw_f
        TYPE(field_t), INTENT(inout) :: bp_f

        ! Local variables
        INTEGER(intk) :: i, igrid, kk, jj, ii
        REAL(realk), POINTER, CONTIGUOUS, DIMENSION(:, :, :) :: bu, bv, bw, bp

        DO i = 1, nmygrids
            igrid = mygrids(i)
            CALL get_mgdims(kk, jj, ii, igrid)

            CALL bu_f%get_ptr(bu, igrid)
            CALL bv_f%get_ptr(bv, igrid)
            CALL bw_f%get_ptr(bw, igrid)
            CALL bp_f%get_ptr(bp, igrid)

            CALL blockbpfeld_grid(kk, jj, ii, bu, bv, bw, bp)
        END DO
    END SUBROUTINE blockbpfeld


    SUBROUTINE blockbpfeld_grid(kk, jj, ii, bu, bv, bw, bp)
        ! Subroutine arguments
        INTEGER(intk), INTENT(in) :: kk, jj, ii
        REAL(realk), INTENT(in) :: bu(kk, jj, ii), bv(kk, jj, ii), &
            bw(kk, jj, ii)
        REAL(realk), INTENT(out) :: bp(kk, jj, ii)

        ! Local variables
        INTEGER(intk) :: i, j, k

        bp = 0.0

        DO i = 2, ii-1
            DO j = 2, jj-1
                DO k = 2, kk-1
                    bp(k, j, i) = MIN( &
                        bu(k, j, i-1) + bu(k, j, i), &
                        bv(k, j-1, i) + bv(k, j, i), &
                        bw(k-1, j, i) + bw(k, j, i), &
                        1.0_realk)
                END DO
            END DO
        END DO
    END SUBROUTINE blockbpfeld_grid
END MODULE gc_blockbpfeld_mod
