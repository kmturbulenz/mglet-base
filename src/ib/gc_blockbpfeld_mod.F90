MODULE gc_blockbpfeld_mod
    USE core_mod, ONLY: realk, intk, errr, minlevel, maxlevel, nmygridslvl, &
        mygridslvl, get_ip3, get_mgdims, field_t

    IMPLICIT NONE(type, external)
    PRIVATE

    PUBLIC :: blockbpfeld

CONTAINS
    SUBROUTINE blockbpfeld(bu, bv, bw, bp)
        ! Subroutine arguments
        TYPE(field_t), INTENT(in) :: bu, bv, bw
        TYPE(field_t), INTENT(inout) :: bp

        ! Local variables
        INTEGER(intk) :: ilevel

        DO ilevel = minlevel, maxlevel
            CALL blockbpfeld_level(ilevel, bu, bv, bw, bp)
        END DO
    END SUBROUTINE blockbpfeld


    SUBROUTINE blockbpfeld_level(ilevel, bu, bv, bw, bp)
        ! Subroutine arguments
        INTEGER(intk), INTENT(in) :: ilevel
        TYPE(field_t), INTENT(in) :: bu, bv, bw
        TYPE(field_t), INTENT(inout) :: bp

        ! Local variables
        INTEGER(intk) :: i, igrid, kk, jj, ii, ip3

        DO i = 1, nmygridslvl(ilevel)
            igrid = mygridslvl(i, ilevel)

            CALL get_mgdims(kk, jj, ii, igrid)
            CALL get_ip3(ip3, igrid)
            CALL blockbpfeld_grid(kk, jj, ii, bu%arr(ip3), bv%arr(ip3), &
                bw%arr(ip3), bp%arr(ip3))
        END DO
    END SUBROUTINE blockbpfeld_level


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
                        bu(k, j, i-1) + bu(k,j,i), &
                        bv(k, j-1, i) + bv(k, j, i), &
                        bw(k-1, j, i) + bw(k, j, i), &
                        1.0_realk)
                END DO
            END DO
        END DO
    END SUBROUTINE blockbpfeld_grid
END MODULE gc_blockbpfeld_mod
