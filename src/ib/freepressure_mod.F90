MODULE freepressure_mod
    USE core_mod, ONLY: realk, intk, int32, mygrids, nmygrids, &
        mygridslvl, nmygridslvl, minlevel, maxlevel, errr, &
        field_t, get_mgdims, get_ip3, connect

    IMPLICIT NONE(type, external)
    PRIVATE

    PUBLIC :: freepressure

CONTAINS
    SUBROUTINE freepressure(kanteu, kantev, kantew, knoten)
        ! Subroutine arguments
        TYPE(field_t), INTENT(in) :: kanteu, kantev, kantew
        TYPE(field_t), INTENT(inout) :: knoten

        ! Local variables
        INTEGER(intk) :: ilevel, iloop

        DO iloop = 1, 4
            DO ilevel = minlevel, maxlevel
                CALL freepressure_level(ilevel, kanteu, kantev, &
                    kantew, knoten)
            END DO
        END DO

        DO ilevel = minlevel, maxlevel
            CALL connect(ilevel, 2, s1=knoten, corners=.TRUE.)
        END DO

        DO ilevel = minlevel, maxlevel
            CALL correct_knoten_level(ilevel, knoten)
        END DO
    END SUBROUTINE freepressure


    SUBROUTINE freepressure_level(ilevel, kanteu, kantev, kantew, knoten)
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
            CALL freepressure_grid(kk, jj, ii, kanteu%arr(ip3), &
                kantev%arr(ip3), kantew%arr(ip3), knoten%arr(ip3))
        END DO
    END SUBROUTINE freepressure_level


    SUBROUTINE freepressure_grid(kk, jj, ii, kanteu, kantev, kantew, knoten)
        ! Subroutine arguments
        INTEGER(intk), INTENT(in) :: kk, jj, ii
        REAL(realk), INTENT(in) :: kanteu(kk, jj, ii), kantev(kk, jj, ii), &
            kantew(kk, jj, ii)
        REAL(realk), INTENT(inout) :: knoten(kk, jj, ii)

        ! Local variables
        INTEGER(intk) :: i, j, k
        INTEGER(intk) :: i0, i1, j0, j1, k0, k1

        DO i = 1, ii
            DO j = 1, jj
                DO k = 1, kk
                    knoten(k, j, i) = MAX(0.0_realk, knoten(k, j, i))
                END DO
            END DO
        END DO

        DO i = 1, ii
            DO j = 1, jj
                DO k = 1, kk
                    i1 = MIN(i+1, ii)
                    IF (NINT(knoten(k, j, i1)*kanteu(k, j, i1)) == 1) THEN
                        IF (NINT(knoten(k, j, i)) /= 1) knoten(k, j, i) = 2.0
                    END IF

                    i0 = MAX(i-1, 1)
                    IF (NINT(knoten(k, j, i0)*kanteu(k, j, i)) == 1) THEN
                        IF (NINT(knoten(k, j, i)) /= 1) knoten(k, j, i) = 2.0
                    END IF

                    j1 = MIN(j+1, jj)
                    IF (NINT(knoten(k, j1, i)*kantev(k, j1, i)) == 1) THEN
                        IF (NINT(knoten(k, j, i)) /= 1) knoten(k, j, i) = 2.0
                    END IF

                    j0 = MAX(j-1, 1)
                    IF (NINT(knoten(k, j0, i)*kantev(k, j, i)) == 1) THEN
                        IF (NINT(knoten(k, j, i)) /= 1) knoten(k, j, i) = 2.0
                    END IF

                    k1 = MIN(k+1, kk)
                    IF (NINT(knoten(k1, j, i)*kantew(k1, j, i)) == 1) THEN
                        IF (NINT(knoten(k, j, i)) /= 1) knoten(k, j, i) = 2.0
                    END IF

                    k0 = MAX(k-1, 1)
                    IF (NINT(knoten(k0, j, i)*kantew(k, j, i)) == 1) THEN
                        IF (NINT(knoten(k, j, i)) /= 1) knoten(k, j, i) = 2.0
                    END IF
                END DO
            END DO
        END DO

        DO i = 1, ii
            DO j = 1, jj
                DO k = 1, kk
                    knoten(k, j, i) = MIN(1.0_realk, knoten(k, j, i))
                END DO
            END DO
        END DO
    END SUBROUTINE freepressure_grid


    SUBROUTINE correct_knoten_level(ilevel, knoten)
        ! Subroutine arguments
        INTEGER(intk), INTENT(in) :: ilevel
        TYPE(field_t), INTENT(inout) :: knoten

        ! Local variables
        INTEGER(intk) :: i, igrid, kk, jj, ii, ip3

        DO i = 1, nmygridslvl(ilevel)
            igrid = mygridslvl(i, ilevel)

            CALL get_mgdims(kk, jj, ii, igrid)
            CALL get_ip3(ip3, igrid)
            CALL correct_knoten_grid(kk, jj, ii, knoten%arr(ip3))
        END DO
    END SUBROUTINE correct_knoten_level


    SUBROUTINE correct_knoten_grid(kk, jj, ii, knoten)
        ! Subroutine arguments
        INTEGER(intk), INTENT(in) :: kk, jj, ii
        REAL(realk), INTENT(inout) :: knoten(kk, jj, ii)

        ! Local variables
        INTEGER(intk) :: i, j, k

        ! Code previously found in blockbp.F and blockbpcc_subroutines.F
        ! TODO: Ask Johannes about this
        DO i = 1, ii
            DO j = 1, jj
                DO k = 1, kk
                    knoten(k, j, i) = MIN(knoten(k, j, i), 1.0_realk)
                END DO
            END DO
        END DO
    END SUBROUTINE correct_knoten_grid
END MODULE freepressure_mod
