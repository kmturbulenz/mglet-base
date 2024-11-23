MODULE gc_blockface_mod
    USE core_mod, ONLY: intk, realk, field_t, nmygrids, mygrids, get_mgdims
    IMPLICIT NONE(type, external)
    PRIVATE

    PUBLIC :: blockface
CONTAINS
    SUBROUTINE blockface(knoten_f, bu_f, bv_f, bw_f)
        ! Subroutine arguments
        TYPE(field_t), INTENT(in) :: knoten_f
        TYPE(field_t), INTENT(inout) :: bu_f, bv_f, bw_f

        ! Local variables
        INTEGER(intk) :: i, igrid, kk, jj, ii
        REAL(realk), POINTER, CONTIGUOUS :: knoten(:, :, :), bu(:, :, :), &
            bv(:, :, :), bw(:, :, :)

        DO i = 1, nmygrids
            igrid = mygrids(i)
            CALL get_mgdims(kk, jj, ii, igrid)

            CALL knoten_f%get_ptr(knoten, igrid)
            CALL bu_f%get_ptr(bu, igrid)
            CALL bv_f%get_ptr(bv, igrid)
            CALL bw_f%get_ptr(bw, igrid)

            CALL blockface_grid(kk, jj, ii, knoten, bu, bv, bw)
        END DO
    END SUBROUTINE blockface


    SUBROUTINE blockface_grid(kk, jj, ii, knoten, bu, bv, bw)
        ! Subroutine arguments
        INTEGER(intk), INTENT(in) :: kk, jj, ii
        REAL(realk), INTENT(in) :: knoten(kk, jj, ii)
        REAL(realk), INTENT(out) :: bu(kk, jj, ii), bv(kk, jj, ii), &
            bw(kk, jj, ii)

        ! Local variables
        INTEGER(intk) :: k, j, i

        bu = 0
        bv = 0
        bw = 0

        ! Sobald ein Knoten Koerper ist, ist die Flaeche geblockt
        DO i = 1, ii
            DO j = 2, jj
                DO k = 2, kk
                    bu(k, j, i) = MIN( &
                        knoten(k-1, j-1, i), &
                        knoten(k, j-1, i), &
                        knoten(k-1, j, i), &
                        knoten(k, j, i))
                END DO
            END DO
        END DO

        DO i = 2, ii
            DO j = 1, jj
                DO k = 2, kk
                    bv(k, j, i) = MIN( &
                        knoten(k-1, j, i-1), &
                        knoten(k, j, i-1), &
                        knoten(k-1, j, i), &
                        knoten(k, j, i))
                END DO
            END DO
        END DO

        DO i = 2, ii
            DO j = 2, jj
                DO k = 1, kk
                    bw(k, j, i) = MIN( &
                        knoten(k, j-1, i-1), &
                        knoten(k, j-1, i), &
                        knoten(k, j, i-1), &
                        knoten(k, j, i))
                END DO
            END DO
        END DO
    END SUBROUTINE blockface_grid
END MODULE gc_blockface_mod
