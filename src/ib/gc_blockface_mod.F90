MODULE gc_blockface_mod
    USE core_mod, ONLY: intk, realk, errr, minlevel, field_t, &
        maxlevel, nmygridslvl, mygridslvl, get_mgdims, get_ip3
    IMPLICIT NONE(type, external)
    PRIVATE

    PUBLIC :: blockface
CONTAINS
    SUBROUTINE blockface(knoten, bu, bv, bw)
        ! Subroutine arguments
        TYPE(field_t), INTENT(inout) :: knoten
        TYPE(field_t), INTENT(inout) :: bu, bv, bw

        ! Local variables
        INTEGER(intk) :: ilevel

        DO ilevel = minlevel, maxlevel
            CALL blockface_level(ilevel, knoten, bu, bv, bw)
        END DO
    END SUBROUTINE blockface


    SUBROUTINE blockface_level(ilevel, knoten, bu, bv, bw)
        ! Subroutine arguments
        INTEGER(intk), INTENT(in) :: ilevel
        TYPE(field_t), INTENT(inout) :: knoten
        TYPE(field_t), INTENT(inout) :: bu, bv, bw

        ! Local variables
        INTEGER(intk) :: i, igrid, kk, jj, ii, ip3

        DO i = 1, nmygridslvl(ilevel)
            igrid = mygridslvl(i, ilevel)

            CALL get_mgdims(kk, jj, ii, igrid)
            CALL get_ip3(ip3, igrid)
            CALL blockface_grid(kk, jj, ii, knoten%arr(ip3), &
                bu%arr(ip3), bv%arr(ip3), bw%arr(ip3))
        END DO
    END SUBROUTINE blockface_level


    SUBROUTINE blockface_grid(kk, jj, ii, knoten, bu, bv, bw)
        ! Subroutine arguments
        INTEGER(intk), INTENT(IN) :: kk, jj, ii
        REAL(realk), INTENT(inout) :: knoten(kk, jj, ii)
        REAL(realk), INTENT(OUT) :: bu(kk, jj, ii), bv(kk, jj, ii), &
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
