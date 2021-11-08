MODULE flzelle_mod
    USE core_mod, ONLY: realk, intk, errr, divide0

    IMPLICIT NONE(type, external)
    PRIVATE

    PUBLIC :: flzelle

CONTAINS

    SUBROUTINE flzelle(k, j, i, kk, jj, ii, icheck, &
            knotenzelle, kanten, nka, connectka, nkn, connectkn, &
            nwa, connectwa)

        ! Bekommt Knoten und Kanten,
        ! (knoten(i) = 0: Koerper, = 1: Fluid
        !  kanten(i) = 0: kein Schnittpunkt, = 1: Schnittpunkt)
        ! liefert die Punkte der Begrenzungsflaechen:
        !    connectka: Kanten-Schnittpunkte fuer jede angeschnittene Ebene
        !    connectkn: Knoten fuer jede angeschnittene Ebene
        !    connectwa: Wandflaeche (nur Kanten-Schnittpunkte)
        !
        !
        ! Nummerierung der Knoten und Kanten:
        !
        !           (8)---7---(7)
        !           /|        /|
        !            8         6
        !            |         |
        !           (5)---5---(6)
        !           /         /
        !
        !      12        11
        !      /         /
        !     /         /
        !   (4)---3---(3)               Z
        !    |  9      | 10             ^   Y
        !    4 /       2 /              | /
        !    |/        |/               |/
        !   (1)---1---(2)               |----> X
        !
        !
        ! Nummerierung der Flaechen:
        !   1 FRO (x-)
        !   2 BAC (x+)
        !   3 RGT (y-)
        !   4 LFT (y+)
        !   5 BOT (z-)
        !   6 TOP (z+)

        INTEGER(intk), INTENT(in) :: k, j, i
        INTEGER(intk), INTENT(in) :: kk, jj, ii
        INTEGER(intk), INTENT(in) :: icheck
        INTEGER(intk), INTENT(in) :: knotenzelle(8), kanten(12)
        INTEGER(intk), INTENT(out) :: nka(6), connectka(2, 6)
        INTEGER(intk), INTENT(out) :: nkn(6), connectkn(4, 6)
        INTEGER(intk), INTENT(out) :: nwa, connectwa(7)

        ! icheck = 1: first ghost layer is relevant
        ! icheck = 0: first ghost layer is not relevant

        ! Local variables
        INTEGER(intk) :: kantenenden(2, 12)
        INTEGER(intk) :: nknoten, nkanten0, nkanten1
        INTEGER(intk) :: indkanten(4), indknoten(4), kantenfl(4), knotenfl(4)
        INTEGER(intk) :: connectkafl(2), connectknfl(4), nkafl, nknfl
        INTEGER(intk) :: ifl, idx
        INTEGER(intk) :: ierrfloffencon, ierrflwandcon

        ! Zeroize all inten(out)
        nka = 0
        connectka = 0
        nkn = 0
        connectkn = 0
        nwa = 0
        connectwa = 0

        kantenenden(1, 1) = 1
        kantenenden(2, 1) = 2
        kantenenden(1, 2) = 2
        kantenenden(2, 2) = 3
        kantenenden(1, 3) = 3
        kantenenden(2, 3) = 4
        kantenenden(1, 4) = 4
        kantenenden(2, 4) = 1
        kantenenden(1, 5) = 5
        kantenenden(2, 5) = 6
        kantenenden(1, 6) = 6
        kantenenden(2, 6) = 7
        kantenenden(1, 7) = 7
        kantenenden(2, 7) = 8
        kantenenden(1, 8) = 8
        kantenenden(2, 8) = 5
        kantenenden(1, 9) = 1
        kantenenden(2, 9) = 5
        kantenenden(1, 10) = 2
        kantenenden(2, 10) = 6
        kantenenden(1, 11) = 3
        kantenenden(2, 11) = 7
        kantenenden(1, 12) = 4
        kantenenden(2, 12) = 8


        nknoten = 0
        DO idx = 1, 8
            nknoten = nknoten + knotenzelle(idx)
        END DO

        nkanten0 = 0
        nkanten1 = 0
        DO idx = 1, 12
            IF (knotenzelle(kantenenden(1, idx)) == &
                    knotenzelle(kantenenden(2, idx))) THEN
                IF (knotenzelle(kantenenden(1, idx)) == 1) THEN
                    nkanten1 = nkanten1 + 1
                ELSE IF (knotenzelle(kantenenden(1, idx)) == 0) THEN
                    nkanten0 = nkanten0 + 1
                END IF
            END IF
        END DO

        IF (nknoten > nkanten1 + 1) THEN
            IF (k == kk .OR. j == jj .OR. i == ii) then
                ! in second ghost layer not relevant
                nka = 0
                nkn = 0
                nwa = 0
                RETURN
            ELSE
                WRITE(*, *) 'flzelle: offene knoten nicht zusammenhaengend'
                WRITE(*, *) 'flzelle: nknoten, nkanten1', nknoten, nkanten1
                WRITE(*, *) 'flzelle: i, j, k', i, j, k
                CALL errr(__FILE__, __LINE__)
            END IF
        END IF

        IF (8-nknoten > nkanten0+1) THEN
            if (k == kk .OR. j == jj .OR. i == ii) THEN
                ! in second ghost layer not relevant
                nka = 0
                nkn = 0
                nwa = 0
                RETURN
            ELSE
                WRITE(*, *) 'flzelle: geschl. knoten nicht zusammenhaengend'
                WRITE(*, *) 'flzelle: nknoten, nkanten0', 8-nknoten, nkanten0
                WRITE(*, *) 'flzelle: i, j, k', i, j, k
                CALL errr(__FILE__, __LINE__)
            END IF
        END IF


        ! fro (x-)
        ifl = 1

        indkanten(1) = 4
        indkanten(2) = 12
        indkanten(3) = 8
        indkanten(4) = 9

        indknoten(1) = 1
        indknoten(2) = 4
        indknoten(3) = 8
        indknoten(4) = 5

        DO idx = 1, 4
            kantenfl(idx) = kanten(indkanten(idx))
            knotenfl(idx) = knotenzelle(indknoten(idx))
        END DO

        CALL floffencon(knotenfl, kantenfl, nkafl, connectkafl, nknfl, &
            connectknfl, ierrfloffencon)

        IF (ierrfloffencon > 0) THEN
            IF (k == kk .OR. j == jj .OR. i == ii) THEN
                ! in second ghost layer not relevant
                nkafl = 0
                nknfl = 0
            ELSE IF (MIN(1, MAX(0, MIN(k-2, kk-1-k))) &
                    + MIN(1, MAX(0, MIN(j-2, jj-1-j))) &
                    + MIN(1, MAX(0, MIN(i-2, ii-1-i))) <= 1) THEN
                ! in corner gost layer not relevant
                nkafl = 0
                nknfl = 0
            ELSE IF (icheck == 0 .AND. &
                    (k == 2 .OR. j == 2 .OR. i == 2 .OR. &
                    k == kk-1 .OR. j == jj-1 .OR. i == ii-1)) THEN
                ! in first ghost layer not relevant
                nkafl = 0
                nknfl = 0
            ELSE
                WRITE(*, *) "floffencon", k, j, i, ifl
                CALL errr(__FILE__, __LINE__)
            END IF
        END IF

        nka(ifl) = nkafl
        nkn(ifl) = nknfl
        DO idx = 1, nkafl
            connectka(idx, ifl) =  indkanten(connectkafl(idx))
        END DO
        DO idx = 1, nknfl
            connectkn(idx, ifl) =  indknoten(connectknfl(idx))
        END DO

        !!! bac (x+)
        ifl = 2

        indkanten(1) = 10
        indkanten(2) = 6
        indkanten(3) = 11
        indkanten(4) = 2

        indknoten(1) = 2
        indknoten(2) = 6
        indknoten(3) = 7
        indknoten(4) = 3

        DO idx = 1, 4
            kantenfl(idx) = kanten(indkanten(idx))
            knotenfl(idx) = knotenzelle(indknoten(idx))
        END DO

        CALL floffencon(knotenfl, kantenfl, nkafl, connectkafl, nknfl, &
            connectknfl, ierrfloffencon)
        IF (ierrfloffencon > 0) then
            IF (k == kk .OR. j == jj .OR. i == ii) THEN
                ! in second ghost layer not relevant
                nkafl = 0
                nknfl = 0
            ELSE IF (MIN(1, MAX(0, MIN(k-2, kk-1-k))) &
                    +MIN(1, MAX(0, MIN(j-2, jj-1-j))) &
                    +MIN(1, MAX(0, MIN(i-2, ii-1-i))) <= 1) THEN
                ! in corner gost layer not relevant
                nkafl = 0
                nknfl = 0
            ELSE IF (icheck == 0 .AND. &
                    (k == 2 .OR. j == 2 .OR. i == 2 .OR. &
                    k == kk-1 .OR. j == jj-1 .OR. i == ii-1)) THEN
                ! in first ghost layer not relevant
                nkafl = 0
                nknfl = 0
            ELSE
                WRITE(*, *) "floffencon", k, j, i, ifl
                CALL errr(__FILE__, __LINE__)
            END IF
        END IF

        nka(ifl) = nkafl
        nkn(ifl) = nknfl
        DO idx = 1, nkafl
            connectka(idx, ifl) =  indkanten(connectkafl(idx))
        END DO
        DO idx = 1, nknfl
            connectkn(idx, ifl) =  indknoten(connectknfl(idx))
        END DO

        !!! rgt (y-)
        ifl = 3

        indkanten(1) = 1
        indkanten(2) = 2
        indkanten(3) = 3
        indkanten(4) = 4

        indknoten(1) = 1
        indknoten(2) = 2
        indknoten(3) = 3
        indknoten(4) = 4

        DO idx = 1,4
            kantenfl(idx) = kanten(indkanten(idx))
            knotenfl(idx) = knotenzelle(indknoten(idx))
        END DO

        CALL floffencon(knotenfl, kantenfl, nkafl, connectkafl, nknfl, &
            connectknfl, ierrfloffencon)
        IF (ierrfloffencon > 0) then
            IF (k == kk .OR. j == jj .OR. i == ii) THEN
                ! in second ghost layer not relevant
                nkafl = 0
                nknfl = 0
            ELSE IF (MIN(1, MAX(0, MIN(k-2, kk-1-k))) &
                    + MIN(1, MAX(0, MIN(j-2, jj-1-j))) &
                    + MIN(1, MAX(0, MIN(i-2, ii-1-i))) <= 1) THEN
                ! in corner gost layer not relevant
                nkafl = 0
                nknfl = 0
            ELSE IF (icheck == 0 .AND. &
                    (k == 2 .OR. j == 2 .OR. i == 2 .OR. &
                    k == kk-1 .OR. j == jj-1 .OR. i == ii-1)) THEN
                ! in first ghost layer not relevant
                nkafl = 0
                nknfl = 0
            ELSE
                WRITE(*, *) "floffencon", k, j, i, ifl
                CALL errr(__FILE__, __LINE__)
            END IF
        END IF

        nka(ifl) = nkafl
        nkn(ifl) = nknfl
        DO idx = 1, nkafl
            connectka(idx, ifl) =  indkanten(connectkafl(idx))
        END DO
        DO idx = 1, nknfl
            connectkn(idx, ifl) =  indknoten(connectknfl(idx))
        END DO

        !!! lft (y+)
        ifl = 4

        indkanten(1) = 8
        indkanten(2) = 7
        indkanten(3) = 6
        indkanten(4) = 5

        indknoten(1) = 5
        indknoten(2) = 8
        indknoten(3) = 7
        indknoten(4) = 6

        DO idx = 1, 4
            kantenfl(idx) = kanten(indkanten(idx))
            knotenfl(idx) = knotenzelle(indknoten(idx))
        END DO

        CALL floffencon(knotenfl, kantenfl, nkafl, connectkafl, nknfl, &
            connectknfl, ierrfloffencon)
        IF (ierrfloffencon > 0) THEN
            IF (k == kk .OR. j == jj .OR. i == ii) THEN
                ! in second ghost layer not relevant
                nkafl = 0
                nknfl = 0
            ELSE IF (MIN(1, MAX(0, MIN(k-2, kk-1-k))) &
                    + MIN(1, MAX(0, MIN(j-2, jj-1-j))) &
                    + MIN(1, MAX(0, MIN(i-2, ii-1-i))) <= 1) THEN
                ! in corner gost layer not relevant
                nkafl = 0
                nknfl = 0
            ELSE IF (icheck == 0 .AND. &
                    (k == 2 .OR. j == 2 .OR. i == 2 .OR. &
                    k == kk-1 .OR. j == jj-1 .OR. i == ii-1)) THEN
                ! in first ghost layer not relevant
                nkafl = 0
                nknfl = 0
            ELSE
                WRITE(*, *) "floffencon", k, j, i, ifl
                CALL errr(__FILE__, __LINE__)
            END IF
        END IF

        nka(ifl) = nkafl
        nkn(ifl) = nknfl
        DO idx = 1, nkafl
            connectka(idx, ifl) =  indkanten(connectkafl(idx))
        END DO
        DO idx = 1, nknfl
            connectkn(idx, ifl) =  indknoten(connectknfl(idx))
        END DO

        !!! bot (z-)
        ifl = 5

        indkanten(1) = 9
        indkanten(2) = 5
        indkanten(3) = 10
        indkanten(4) = 1

        indknoten(1) = 1
        indknoten(2) = 5
        indknoten(3) = 6
        indknoten(4) = 2

        DO idx = 1, 4
            kantenfl(idx) = kanten(indkanten(idx))
            knotenfl(idx) = knotenzelle(indknoten(idx))
        END DO

        CALL floffencon(knotenfl, kantenfl, nkafl, connectkafl, nknfl, &
            connectknfl, ierrfloffencon)
        IF (ierrfloffencon > 0) THEN
            IF (k == kk .OR. j == jj .OR. i == ii) THEN
                ! in second ghost layer not relevant
                nkafl = 0
                nknfl = 0
            ELSE IF (MIN(1, MAX(0, MIN(k-2, kk-1-k))) &
                    + MIN(1, MAX(0, MIN(j-2, jj-1-j))) &
                    + MIN(1, MAX(0, MIN(i-2, ii-1-i))) <= 1) THEN
                ! in corner gost layer not relevant
                nkafl = 0
                nknfl = 0
            ELSE IF (icheck == 0 .AND. &
                    (k == 2 .OR. j == 2 .OR. i == 2 .OR. &
                    k == kk-1 .OR. j == jj-1 .OR. i == ii-1)) THEN
                ! in first ghost layer not relevant
                nkafl = 0
                nknfl = 0
            ELSE
                WRITE(*, *) "floffencon", k, j, i, ifl
                CALL errr(__FILE__, __LINE__)
            END IF
        END IF

        nka(ifl) = nkafl
        nkn(ifl) = nknfl
        DO idx = 1, nkafl
            connectka(idx, ifl) =  indkanten(connectkafl(idx))
        END DO
        DO idx = 1,nknfl
            connectkn(idx, ifl) =  indknoten(connectknfl(idx))
        END DO

        !!! top (z+)
        ifl = 6

        indkanten(1) = 11
        indkanten(2) = 7
        indkanten(3) = 12
        indkanten(4) = 3

        indknoten(1) = 3
        indknoten(2) = 7
        indknoten(3) = 8
        indknoten(4) = 4

        DO idx = 1,4
            kantenfl(idx) = kanten(indkanten(idx))
            knotenfl(idx) = knotenzelle(indknoten(idx))
        END DO

        CALL floffencon(knotenfl, kantenfl, nkafl, connectkafl, nknfl, &
            connectknfl, ierrfloffencon)
        IF (ierrfloffencon > 0) THEN
            IF (k == kk .OR. j == jj .OR. i == ii) THEN
                ! in second ghost layer not relevant
                nkafl = 0
                nknfl = 0
            ELSE IF (MIN(1, MAX(0, MIN(k-2, kk-1-k))) &
                    + MIN(1, MAX(0, MIN(j-2, jj-1-j))) &
                    + MIN(1, MAX(0, MIN(i-2, ii-1-i))) <= 1) THEN
                ! in corner gost layer not relevant
                nkafl = 0
                nknfl = 0
            ELSE IF (icheck == 0 .AND. &
                    (k == 2 .OR. j == 2 .OR. i == 2 .OR. &
                    k == kk-1 .OR. j == jj-1 .OR. i == ii-1)) THEN
                ! in first ghost layer not relevant
                nkafl = 0
                nknfl = 0
            ELSE
                WRITE(*, *) "floffencon", k, j, i, ifl
                CALL errr(__FILE__, __LINE__)
            END IF
        END IF

        nka(ifl) = nkafl
        nkn(ifl) = nknfl
        DO idx = 1, nkafl
            connectka(idx, ifl) =  indkanten(connectkafl(idx))
        END DO
        DO idx = 1, nknfl
            connectkn(idx, ifl) =  indknoten(connectknfl(idx))
        END DO


        CALL flwandcon(nka, connectka, nwa, connectwa, ierrflwandcon)
        IF (ierrflwandcon > 0) THEN
            IF (k == kk .OR. j == jj .OR. i == ii) THEN
                ! in second ghost layer not relevant
                nwa = 0
            ELSE IF (MIN(1, MAX(0, MIN(k-2, kk-1-k))) &
                    + MIN(1, MAX(0, MIN(j-2, jj-1-j))) &
                    + MIN(1, MAX(0, MIN(i-2, ii-1-i))) <=1 ) THEN
                ! in corner gost layer not relevant
                nwa = 0
            ELSE IF (icheck == 0 .AND. &
                    (k == 2 .OR. j == 2 .OR. i == 2 .OR. &
                    k == kk-1 .OR. j == jj-1 .OR. i == ii-1)) THEN
                ! in first ghost layer not relevant
                nwa = 0
            ELSE
                WRITE(*, *) "flwandcon", k, j, i, ifl
                CALL errr(__FILE__, __LINE__)
            END IF
        END IF
    END SUBROUTINE flzelle


    SUBROUTINE floffencon(knotenfl, kanten, nka, connectka, nkn, &
            connectkn, ierr)
        ! Verbindet die punkte einer flaeche
        !
        !   (4)---3---(3)
        !    |         |
        !    4         2
        !    |         |
        !   (1)---1---(2)
        !
        ! rechte hand nach aussen

        ! Subroutine arguments
        INTEGER(intk), INTENT(in) :: knotenfl(4), kanten(4)
        INTEGER(intk), INTENT(out) :: nka, connectka(2), nkn, connectkn(4)
        INTEGER(intk), INTENT(out) :: ierr

        ! Local variables
        INTEGER(intk) :: idx, nschnitt, nknoten

        ! Zeroize INTENT(out)
        nka = 0
        connectka = 0
        nkn = 0
        connectkn = 0
        ierr = 0

        DO idx = 1, 4
            IF (knotenfl(idx) /= knotenfl(MODULO(idx, 4)+1)) THEN
                IF (kanten(idx) /= 1) THEN
                    ! Aufeinanderfolgende knoten verschieden, aber kante
                    ! ungeschnitten
                    WRITE(*, *) 'floffencon, knotenfl', knotenfl
                    WRITE(*, *) 'floffencon, kanten', kanten
                    ierr = 502
                    RETURN
                END IF
            ELSE
                IF (kanten(idx) /= 0) THEN
                    ! Aufeinanderfolgende knoten gleich, aber kante geschnitten
                    WRITE(*, *) 'floffencon, knotenfl', knotenfl
                    WRITE(*, *) 'floffencon, kanten', kanten
                    ierr = 503
                    RETURN
                END IF
            END IF
        END DO

        ! Welchen fall
        nknoten = 0
        DO idx = 1, 4
            nknoten = nknoten + knotenfl(idx)
        END DO

        IF (nknoten == 0) THEN
            RETURN
        END IF

        IF (nknoten == 4) THEN
            nka = 0
            connectka = 0
            nkn = 4
            DO idx = 1, 4
                connectkn(idx) = idx
            END DO
            RETURN
        END IF

        nschnitt = 0
        do idx = 1, 4
            nschnitt = nschnitt + kanten(idx)
        END DO

        IF (nschnitt == 2) THEN
            ! Anzahl geschnittene kanten darf nur 2 sein, wenn die flaeche
            ! angeschniteten ist
            nka = 2
        ELSE
            WRITE(*, *) 'floffencon, kanten', kanten
            ierr = 501
            RETURN
        END IF

        ! Anfangspunkt: offener knotenfl nach kante mit schnittpunkt
        DO idx = 1, 4
            IF (kanten(idx) == 1 .AND. knotenfl(idx) == 1) THEN
                connectka(1) = idx
                EXIT
            END IF
        END DO

        DO idx = 1, 4
            IF (kanten(idx) == 1 .AND. idx /= connectka(1)) THEN
                connectka(2) = idx
            END IF
        END DO

        ! Connect belegen
        nkn = 0
        DO idx = connectka(2), 4-1+connectka(2)
            IF (knotenfl(MODULO(idx, 4)+1) == 1) THEN
                nkn = nkn + 1
                IF (nkn > 3) THEN
                    ! Auf geschnittener Flaeche koennen max. 3 knoten offen
                    ! sein.
                    ierr = 504
                    RETURN
                END IF
                connectkn(nkn) = MODULO(idx, 4) + 1
            ELSE
                EXIT
            END IF
        END DO
    END SUBROUTINE floffencon


    SUBROUTINE flwandcon(nka, connectka, nwa, connectwa, ierr)
        ! Sucht die Kantenschnittpunkte, die die Wandflaeche bilden
        ! (connectwa) braucht die Kantenschnittpunkte der Seitenflaechen
        ! (connectka)

        ! Subroutine arguments
        INTEGER(intk), INTENT(in) :: nka(6)
        INTEGER(intk), INTENT(in) :: connectka(2,6)
        INTEGER(intk), INTENT(out) :: nwa
        INTEGER(intk), INTENT(out) :: connectwa(7)
        INTEGER(intk), INTENT(out) :: ierr

        ! Local variables
        INTEGER(intk) :: ii, istart, found, iii

        ! Zeroize INTENT(out)
        nwa = 0
        connectwa = 0
        ierr = 0

        istart = 0
        DO ii = 1,6
            IF (nka(ii) == 2) THEN
                istart = ii
                EXIT
            END IF
        END DO

        IF (istart == 0) THEN
            nwa = 0
            RETURN
        END IF

        ! Beginne bei den zwei Kanten-Schnittpunkten einer Seitenflaeche
        connectwa(1) = connectka(2, istart)
        connectwa(2) = connectka(1, istart)
        nwa = 2
        IF (connectwa(1) == 0) CALL errr(__FILE__, __LINE__)
        IF (connectwa(2) == 0) CALL errr(__FILE__, __LINE__)

        ! Gehe weiter zum Kantenschnittpunkt der benachbarten Seitenflaeche,
        ! bis der Ausgangspunkt erreicht ist
        DO
            found = 0
            DO ii = 1, 6
                IF (connectka(2, ii) == connectwa(nwa)) THEN
                    nwa = nwa + 1

                    IF (nwa > 7) THEN
                        DO iii = 1, 7
                            WRITE(*, *) 'flwandcon, i, connectwa(i)', &
                                iii, connectwa(iii)
                        END DO
                        DO iii = 1, 6
                            WRITE(*, *) 'flwandcon, i, connectka(1:2,i)', &
                                iii, connectka(1:2, iii)
                        END DO
                        CALL errr(__FILE__, __LINE__)
                    END IF
                    connectwa(nwa) = connectka(1, ii)
                    found = 1
                END IF
            END DO

            IF (found == 0) THEN
                ierr = 502
                nwa = 0
                RETURN
            END IF

            IF (connectwa(nwa) == connectwa(1)) THEN
                nwa = nwa - 1
                EXIT
            END IF
        END DO
    END SUBROUTINE flwandcon
END MODULE flzelle_mod
