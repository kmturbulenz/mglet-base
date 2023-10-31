MODULE gc_zelltyp_mod
    USE core_mod, ONLY: intk, realk, errr, get_fieldptr, ngrid, minlevel, &
        maxlevel, nmygridslvl, mygridslvl, get_mgdims, get_ip3, field_t
    IMPLICIT NONE(type, external)
    PRIVATE

    PUBLIC :: zelltyp, zelltyp_grid
CONTAINS
    SUBROUTINE zelltyp(knoten, bzelltyp, icells)
        ! Subroutine arguments
        TYPE(field_t), INTENT(in) :: knoten
        INTEGER(intk), INTENT(out) :: bzelltyp(*)
        INTEGER(intk), INTENT(out), OPTIONAL :: icells(:)

        ! Local variables
        INTEGER(intk) :: ilevel

        ! Sanity check
        IF (PRESENT(icells)) THEN
            IF (SIZE(icells) /= ngrid) CALL errr(__FILE__, __LINE__)
            icells = 0
        END IF

        DO ilevel = minlevel, maxlevel
            CALL zelltyp_level(ilevel, knoten, bzelltyp, icells)
        END DO
    END SUBROUTINE zelltyp


    SUBROUTINE zelltyp_level(ilevel, knoten, bzelltyp, icells)
        ! Subroutine arguments
        INTEGER(intk), INTENT(in) :: ilevel
        TYPE(field_t), INTENT(in) :: knoten
        INTEGER(intk), INTENT(out) :: bzelltyp(*)
        INTEGER(intk), INTENT(out), OPTIONAL :: icells(:)

        ! Local variables
        INTEGER(intk) :: i, igrid, kk, jj, ii, ip3, icells2

        DO i = 1, nmygridslvl(ilevel)
            igrid = mygridslvl(i, ilevel)

            CALL get_mgdims(kk, jj, ii, igrid)
            CALL get_ip3(ip3, igrid)
            CALL zelltyp_grid(kk, jj, ii, knoten%arr(ip3), &
                bzelltyp(ip3), icells2)

            IF (PRESENT(icells)) icells(igrid) = icells2
        END DO
    END SUBROUTINE zelltyp_level


    SUBROUTINE zelltyp_grid(kk, jj, ii, knoten, bzelltyp, icells)
        INTEGER(intk), INTENT(IN) :: kk, jj, ii
        REAL(realk), INTENT(IN) :: knoten(kk, jj, ii)
        INTEGER(intk), INTENT(OUT) :: bzelltyp(kk, jj, ii)
        INTEGER(intk), INTENT(OUT) :: icells

        INTEGER(intk) :: k, j, i
        REAL(realk) :: sumknoten

        !  Alle acht Knoten anschauen
        DO i = 1, ii
            DO j = 1, jj
                DO k = 1, kk
                    sumknoten = &
                        knoten(k, j, i) &
                        + knoten(k, j, MAX(i-1, 1)) &
                        + knoten(k, MAX(j-1, 1), i) &
                        + knoten(k, MAX(j-1, 1), MAX(i-1, 1)) &
                        + knoten(MAX(k-1, 1), j, i) &
                        + knoten(MAX(k-1, 1), j, MAX(i-1, 1)) &
                        + knoten(MAX(k-1, 1), MAX(j-1, 1), i) &
                        + knoten(MAX(k-1, 1), MAX(j-1, 1), MAX(i-1, 1))
                    bzelltyp(k, j, i) = NINT(sumknoten)
                END DO
            END DO
        END DO

        ! 1: Fluid
        ! 0: Koerper
        ! -1: Interface
        DO i = 1, ii
            DO j = 1, jj
                DO k = 1, kk
                    IF (bzelltyp(k, j, i) > 0 .AND. bzelltyp(k, j, i) < 8) THEN
                        bzelltyp(k, j, i) = -1
                    ELSE IF (bzelltyp(k, j, i) == 8) THEN
                        bzelltyp(k, j, i) = 1
                    ELSE IF (bzelltyp(k, j, i) == 0) THEN
                        CONTINUE
                    ELSE
                        CALL errr(__FILE__, __LINE__)
                    END IF
                END DO
            END DO
        END DO

        icells = 0
        DO i = 2, ii
            DO j = 2, jj
                DO k = 2, kk
                    IF (bzelltyp(k, j, i) == -1) THEN
                        icells = icells + 1
                        bzelltyp(k, j, i) = -icells
                    END IF
                END DO
            END DO
        END DO
    END SUBROUTINE zelltyp_grid
END MODULE gc_zelltyp_mod
