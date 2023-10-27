MODULE bubvbw_mod
    USE core_mod

    IMPLICIT NONE(type, external)
    PRIVATE

    PUBLIC :: bubvbw

CONTAINS
    SUBROUTINE bubvbw(bp, bu, bv, bw)
        ! Subroutine arguments
        TYPE(field_t), INTENT(in) :: bp
        TYPE(field_t), INTENT(inout) :: bu, bv, bw

        ! Local variables
        INTEGER(intk) :: ilevel

        DO ilevel = minlevel, maxlevel
            CALL bubvbw_level(ilevel, bp, bu, bv, bw)
        END DO
    END SUBROUTINE bubvbw


    SUBROUTINE bubvbw_level(ilevel, bp, bu, bv, bw)
        ! Subroutine arguments
        INTEGER(intk), INTENT(in) :: ilevel
        TYPE(field_t), INTENT(in) :: bp
        TYPE(field_t), INTENT(inout) :: bu, bv, bw

        ! Local variables
        INTEGER(intk) :: i, igrid, kk, jj, ii, ip3

        DO i = 1, nmygridslvl(ilevel)
            igrid = mygridslvl(i, ilevel)

            CALL get_mgdims(kk, jj, ii, igrid)
            CALL get_ip3(ip3, igrid)
            CALL bubvbw_grid(kk, jj, ii, bp%arr(ip3), bu%arr(ip3), &
                bv%arr(ip3), bw%arr(ip3))
        END DO

        ! Previously directly in blockbp - this is last place where these
        ! fields are touched
        CALL connect(ilevel, 2, v1=bu, v2=bv, v3=bw, corners=.TRUE.)
    END SUBROUTINE bubvbw_level


    PURE SUBROUTINE bubvbw_grid(kk, jj, ii, bp, bu, bv, bw)
        ! Subroutine arguments
        INTEGER(intk), INTENT(in) :: kk, jj, ii
        REAL(realk), INTENT(in) :: bp(kk, jj, ii)
        REAL(realk), INTENT(out) :: bu(kk, jj, ii), bv(kk, jj, ii), &
            bw(kk, jj, ii)

        ! Local variables
        INTEGER(intk) :: i, j, k

        DO i = 1, ii-1
            DO j = 1, jj
                DO k = 1, kk
                    bu(k, j, i) = bp(k, j, i)*bp(k, j, i+1)
                END DO
            END DO
        END DO
        DO j = 1, jj
            DO k = 1, kk
                bu(k, j, ii) = bp(k, j, ii)
            END DO
        END DO

        DO i = 1, ii
            DO j = 1, jj-1
                DO k = 1, kk
                    bv(k, j, i) = bp(k, j, i)*bp(k, j+1, i)
                END DO
            END DO
        END DO
        DO i = 1, ii
            DO k = 1, kk
                bv(k, jj, i) = bp(k, jj, i)
            END DO
        END DO

        DO i = 1, ii
            DO j = 1, jj
                DO k = 1, kk-1
                    bw(k, j, i) = bp(k, j, i)*bp(k+1, j, i)
                END DO
            END DO
        END DO
        DO i = 1, ii
            DO j = 1, jj
                bw(kk, j, i) = bp(kk, j, i)
            END DO
        END DO
    END SUBROUTINE bubvbw_grid
END MODULE bubvbw_mod
