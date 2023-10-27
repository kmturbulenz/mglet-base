MODULE blocknodes_mod
    USE core_mod, ONLY: realk, intk, idim3d, mygridslvl, nmygridslvl, &
        minlevel, maxlevel, get_fieldptr, get_mgdims, get_ip3, &
        field_t, connect
    USE topol_mod, ONLY: topol_t

    IMPLICIT NONE (type, external)
    PRIVATE

    PUBLIC :: blocknodes

CONTAINS
    SUBROUTINE blocknodes(kanteu, kantev, kantew, knoten)
        ! Subroutine arguments
        TYPE(field_t), INTENT(in) :: kanteu, kantev, kantew
        TYPE(field_t), INTENT(inout) :: knoten

        ! Local variables
        INTEGER(intk) :: ilevel

        DO ilevel = minlevel, maxlevel
            CALL blocknodes_level(ilevel, kanteu, kantev, kantew, knoten)
            CALL connect(ilevel, 2, s1=knoten, corners=.TRUE.)
        END DO
    END SUBROUTINE blocknodes


    SUBROUTINE blocknodes_level(ilevel, kanteu, kantev, kantew, knoten)
        ! Subroutine arguments
        INTEGER(intk), INTENT(in) :: ilevel
        TYPE(field_t), INTENT(in) :: kanteu, kantev, kantew
        TYPE(field_t), INTENT(inout) :: knoten

        ! Local variables
        INTEGER(intk) :: i, igrid, kk, jj, ii, ip3

        DO i = 1, nmygridslvl(ilevel)
            igrid = mygridslvl(i, ilevel)

            CALL get_mgdims(kk, jj, ii, igrid)
            CALL get_ip3(ip3, igrid)
            CALL blocknodes_grid(kk, jj, ii, kanteu%arr(ip3), kantev%arr(ip3), &
                kantew%arr(ip3), knoten%arr(ip3))
        END DO
    END SUBROUTINE blocknodes_level


    PURE SUBROUTINE blocknodes_grid(kk, jj, ii, kanteu, kantev, kantew, knoten)
        ! Subroutine arguments
        INTEGER(intk), INTENT(in) :: kk, jj, ii
        REAL(realk), INTENT(in) :: kanteu(kk, jj, ii), kantev(kk, jj, ii), &
            kantew(kk, jj, ii)
        REAL(realk), INTENT(out) :: knoten(kk, jj, ii)

        ! Local variables
        INTEGER(intk) :: i, j, k, i2, j2, k2

        ! Initialize intent(out) variables
        knoten = 0.0

        ! Alle Knoten einer Druckzelle, von der mind. eine Kante geschnitten
        ! wird, werden geschlossen
        DO i = 1, ii
            DO j = 1, jj
                DO k = 1, kk
                    IF (kanteu(k, j, i) < -0.5) THEN
                        DO i2 = MAX(i-1, 1), i
                            DO j2 = MAX(j-1, 1), MIN(j+1, jj)
                                DO k2 = MAX(k-1, 1), MIN(k+1, kk)
                                    knoten(k2, j2, i2) = -1.0
                                END DO
                            END DO
                        END DO
                    END IF
                    if (kantev(k, j, i) < -0.5) THEN
                        DO i2 = MAX(i-1, 1), MIN(i+1, ii)
                            DO j2 = MAX(j-1, 1), j
                                DO k2 = MAX(k-1, 1), MIN(k+1, kk)
                                    knoten(k2, j2, i2) = -1.0
                                END DO
                            END DO
                        END DO
                    END IF
                    IF (kantew(k, j, i) < -0.5) THEN
                        DO i2 = MAX(i-1, 1), MIN(i+1, ii)
                            DO j2 = MAX(j-1, 1), MIN(j+1, jj)
                                DO k2 = MAX(k-1, 1), k
                                    knoten(k2, j2, i2) = -1.0
                                END DO
                            END DO
                        END DO
                    END IF
                END DO
            END DO
        END DO
    END SUBROUTINE blocknodes_grid
END MODULE blocknodes_mod
