MODULE blockcheck_mod
    USE core_mod, ONLY: realk, intk, int32, mygrids, nmygrids, &
        mygridslvl, nmygridslvl, minlevel, maxlevel, errr, connect, &
        field_t, get_mgdims, get_ip3, get_ip3n, get_fieldptr, idim3d

    IMPLICIT NONE(type, external)
    PRIVATE

    PUBLIC :: blockcheck, blockcheck_grid

CONTAINS
    SUBROUTINE blockcheck(kanteu, kantev, kantew, knoten, flag)
        ! Subroutine arguments
        TYPE(field_t), INTENT(inout) :: kanteu, kantev, kantew
        TYPE(field_t), INTENT(in) :: knoten
        INTEGER(intk), INTENT(in) :: flag

        ! Local variables
        INTEGER(intk) :: ilevel

        DO ilevel = minlevel, maxlevel
            CALL blockcheck_level(ilevel, kanteu, kantev, kantew, knoten, flag)
        END DO
    END SUBROUTINE blockcheck


    SUBROUTINE blockcheck_level(ilevel, kanteu, kantev, kantew, knoten, flag)
        ! Subroutine arguments
        INTEGER(intk), INTENT(in) :: ilevel
        TYPE(field_t), INTENT(inout) :: kanteu, kantev, kantew
        TYPE(field_t), INTENT(in) :: knoten
        INTEGER(intk), INTENT(in) :: flag

        ! Local variables
        INTEGER(intk) :: i, igrid, kk, jj, ii, ip3

        DO i = 1, nmygridslvl(ilevel)
            igrid = mygridslvl(i, ilevel)

            CALL get_mgdims(kk, jj, ii, igrid)
            CALL get_ip3(ip3, igrid)
            CALL blockcheck_grid(kk, jj, ii, kanteu%arr(ip3), &
                kantev%arr(ip3), kantew%arr(ip3), knoten%arr(ip3), flag)
        END DO
    END SUBROUTINE blockcheck_level


    SUBROUTINE blockcheck_grid(kk, jj, ii, kanteu, kantev, kantew, &
            knoten, flag)
        ! Subroutine arguments
        INTEGER(intk), INTENT(in) :: kk, jj, ii
        REAL(realk), INTENT(inout) :: kanteu(kk, jj, ii), &
            kantev(kk, jj, ii), kantew(kk, jj, ii)
        REAL(realk), INTENT(in) :: knoten(kk, jj, ii)
        INTEGER(intk), INTENT(in) :: flag

        ! Local variables
        INTEGER(intk) :: i, j, k

        DO i = 2, ii
            DO j = 1, jj
                DO k = 1, kk
                    IF (flag == 0 .AND. NINT(kanteu(k, j, i)) == 0) THEN
                        WRITE(*, *) 'blockcheck, kanteu=0', k, j, i
                        CALL errr(__FILE__, __LINE__)
                    END IF
                    IF (kanteu(k, j, i) >= 1.5) THEN
                        WRITE(*, *) 'blockcheck, kanteu>1', k, j, i
                        CALL errr(__FILE__, __LINE__)
                    END IF
                    ! kanteu_i = 1, aber ungleiche knoten
                    ! massnahme: kante = 0 (s. punktekoordinaten)
                    IF (NINT(kanteu(k, j, i)) == 1 .AND. &
                            NINT(knoten(k, j, i)) &
                            /= NINT(knoten(k, j, i-1))) THEN
                        kanteu(k, j, i) = 0.0
                    END IF
                    ! kanteu_i =< 0, aber gleiche knoten
                    ! massnahme: kante = 1
                    IF (kanteu(k, j, i) <= 0.5 .AND. &
                            NINT(knoten(k, j, i)) &
                            == NINT(knoten(k, j, i-1))) THEN
                        kanteu(k, j, i) = 1.0
                    END IF
                END DO
            END DO
        END DO

        DO i = 1, ii
            DO j = 2, jj
                DO k = 1, kk
                    IF (flag == 0 .AND. NINT(kantev(k, j, i)) == 0) THEN
                        WRITE(*, *) 'blockcheck, kantev=0', k, j, i
                        CALL errr(__FILE__, __LINE__)
                    END IF
                    IF (kantev(k, j, i) >= 1.5) THEN
                        WRITE(*, *) 'blockcheck, kantev>1', k, j, i
                        CALL errr(__FILE__, __LINE__)
                    END IF
                    ! kanteu_i = 1, aber ungleiche knoten
                    ! massnahme: kante = 0 (s. punktekoordinaten)
                    IF (NINT(kantev(k, j, i)) == 1 .AND. &
                            NINT(knoten(k, j, i)) &
                            /= NINT(knoten(k, j-1, i))) THEN
                        kantev(k, j, i) = 0.0
                    END IF
                    ! kanteu_i =< 0, aber gleiche knoten
                    IF (kantev(k, j, i) <= 0.5 .AND. &
                            NINT(knoten(k, j, i)) &
                            == NINT(knoten(k, j-1, i))) THEN
                        kantev(k, j, i) = 1.0
                    END IF
                END DO
            END DO
        END DO

        DO i = 1, ii
            DO j = 1, jj
                DO k = 2, kk
                    IF (flag == 0 .AND. NINT(kantew(k, j, i)) == 0) THEN
                        WRITE(*, *) 'blockcheck, kantew=0', k, j, i
                        CALL errr(__FILE__, __LINE__)
                    END IF
                    ! kanteu_i > 1
                    IF (kantew(k, j, i) >= 1.5) THEN
                        WRITE(*, *) 'blockcheck, kantew>1', k, j, i
                        CALL errr(__FILE__, __LINE__)
                    END IF
                    ! kanteu_i = 1, aber ungleiche knoten
                    ! massnahme: kante = 0 (s. punktekoordinaten)
                    IF (NINT(kantew(k, j, i)) == 1 .AND. &
                            NINT(knoten(k, j, i)) &
                            /= NINT(knoten(k-1, j, i))) THEN
                        kantew(k, j, i) = 0.0
                    END IF
                    ! kanteu_i =< 0, aber gleiche knoten
                    IF (kantew(k, j, i) <= 0.5 .AND. &
                            NINT(knoten(k, j, i)) &
                            == NINT(knoten(k-1, j, i))) THEN
                        kantew(k, j, i) = 1.0
                    END IF
                END DO
            END DO
        END DO
    END SUBROUTINE blockcheck_grid
END MODULE blockcheck_mod
