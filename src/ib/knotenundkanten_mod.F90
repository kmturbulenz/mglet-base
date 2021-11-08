MODULE knotenundkanten_mod
    USE core_mod, ONLY: realk, intk, errr

    IMPLICIT NONE(type, external)
    PRIVATE

    PUBLIC :: knotenundkanten

CONTAINS

    SUBROUTINE knotenundkanten(k, j, i, kk, jj, ii, icheck, knoten, &
            kanteu, kantev, kantew, knotenzelle, kanten)
        ! Berechnet fuer einen Zelle:
        ! knotenzelle =1: Fluid-Knoten, =0: Koerper-Knoten
        ! kanten =1: geschnittene Kante, =0: ungeschnittene Kante

        ! Subroutine arguments
        INTEGER(intk), INTENT(in) :: k, j, i, kk, jj, ii, icheck
        REAL(realk), INTENT(in) :: knoten(kk, jj, ii)
        REAL(realk), INTENT(in) :: kanteu(kk, jj, ii), kantev(kk, jj, ii), &
            kantew(kk, jj, ii)
        INTEGER(intk), intent(out) :: knotenzelle(8), kanten(12)

        ! icheck = 1: first ghost layer is relevant
        ! icheck = 0: first ghost layer is not relevant

        ! Local variables
        INTEGER(intk) :: nkantecutted, ikante

        ! Initialize INTENT(out)
        kanten = 0

        knotenzelle(1) = NINT(knoten(k-1, j-1, i-1))
        knotenzelle(2) = NINT(knoten(k-1, j-1, i))
        knotenzelle(3) = NINT(knoten(k, j-1, i))
        knotenzelle(4) = NINT(knoten(k, j-1, i-1))
        knotenzelle(5) = NINT(knoten(k-1, j, i-1))
        knotenzelle(6) = NINT(knoten(k-1, j, i))
        knotenzelle(7) = NINT(knoten(k, j, i))
        knotenzelle(8) = NINT(knoten(k, j, i-1))


        IF (kanteu(k-1, j-1, i) <= 0.5) kanten(1) = 1
        IF (kantew(k, j-1, i) <= 0.5) kanten(2) = 1
        IF (kanteu(k, j-1, i) <= 0.5) kanten(3) = 1
        IF (kantew(k, j-1, i-1) <= 0.5) kanten(4) = 1

        IF (kanteu(k-1, j, i) <= 0.5) kanten(5) = 1
        IF (kantew(k, j, i) <= 0.5) kanten(6) = 1
        IF (kanteu(k, j, i) <= 0.5) kanten(7) = 1
        IF (kantew(k, j, i-1) <= 0.5) kanten(8) = 1

        IF (kantev(k-1, j, i-1) <= 0.5) kanten(9) = 1
        IF (kantev(k-1, j, i) <= 0.5) kanten(10) = 1
        IF (kantev(k, j, i) <= 0.5) kanten(11) = 1
        IF (kantev(k, j, i-1) <= 0.5) kanten(12) = 1

        ! Sicherheitsabfragen
        nkantecutted = 0
        DO ikante = 1, 12
            IF (kanten(ikante) == 1) THEN
                nkantecutted = nkantecutted + 1
            END IF
        END DO
        IF (nkantecutted < 3 .OR. nkantecutted > 6) THEN
            if (k == kk .OR. j == jj .OR. i == ii) then
                ! in second ghost layer not relevant
                kanten=0
                knotenzelle=0
            ELSE IF (MIN(1, MAX(0, MIN(k-2, kk-1-k))) &
                    + MIN(1, MAX(0, MIN(j-2, jj-1-j))) &
                    + MIN(1, MAX(0, MIN(i-2, ii-1-i))) <= 1) THEN
                ! in corner gost layer not relevant
                kanten=0
                knotenzelle=0
            ELSE IF (icheck == 0 .AND. &
                    (k == 2 .OR. j == 2 .OR. i == 2 .OR. &
                    k == kk-1 .OR. j == jj-1 .OR. i == ii-1)) THEN
                ! in first ghost layer not relevant
                kanten=0
                knotenzelle=0
            ELSE IF (MIN(1, MAX(0, MIN(k-2, kk-1-k))) &
                    + MIN(1, MAX(0, MIN(j-2, jj-1-j))) &
                    + MIN(1, MAX(0, MIN(i-2, ii-1-i))) <= 1) THEN
                ! in corner gost layer not relevant
                WRITE(*, *) "Warning 301 in knotenundkanten"
                kanten = 0
                knotenzelle = 0
            ELSE
                WRITE(*, *) "k, j, i", k, j, i
                WRITE(*, *) "nkantecutted", nkantecutted
                DO ikante = 1, 12
                    WRITE(*, *) "kante(i)", ikante, kanten(ikante)
                END DO
                CALL errr(__FILE__, __LINE__)
            END IF
        END IF
    END SUBROUTINE knotenundkanten

END MODULE knotenundkanten_mod
