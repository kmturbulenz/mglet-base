MODULE freekante_mod
    USE core_mod, ONLY: realk, intk, int32, mygrids, nmygrids, &
        mygridslvl, nmygridslvl, minlevel, maxlevel, errr, &
        field_t, get_mgdims, get_ip3, get_ip3n, get_fieldptr
    USE topol_mod, ONLY: topol_t
    USE punktekoordinaten_mod, ONLY: punkteeinekante2

    IMPLICIT NONE(type, external)
    PRIVATE

    PUBLIC :: freekante

CONTAINS
    SUBROUTINE freekante(topol, ntrimax, knoten, kanteu, kantev, &
        kantew, triau, triav, triaw)

        ! Subroutine arguments
        TYPE(topol_t), INTENT(in) :: topol
        INTEGER(intk), INTENT(in) :: ntrimax
        TYPE(field_t), INTENT(inout) :: knoten
        TYPE(field_t), INTENT(inout) :: kanteu
        TYPE(field_t), INTENT(inout) :: kantev
        TYPE(field_t), INTENT(inout) :: kantew
        INTEGER(intk), INTENT(in) :: triau(*)
        INTEGER(intk), INTENT(in) :: triav(*)
        INTEGER(intk), INTENT(in) :: triaw(*)

        ! Local variables
        INTEGER(intk) :: ilevel

        DO ilevel = minlevel, maxlevel
            CALL freekante_level(ilevel, topol, ntrimax, knoten, kanteu, &
                    kantev, kantew, triau, triav, triaw)
        END DO
    END SUBROUTINE freekante


    SUBROUTINE freekante_level(ilevel, topol, ntrimax, knoten, kanteu, &
            kantev, kantew, triau, triav, triaw)
        ! Subroutine arguments
        INTEGER(intk), INTENT(in) :: ilevel
        TYPE(topol_t), INTENT(in) :: topol
        INTEGER(intk), INTENT(in) :: ntrimax
        TYPE(field_t), INTENT(inout) :: knoten
        TYPE(field_t), INTENT(inout) :: kanteu
        TYPE(field_t), INTENT(inout) :: kantev
        TYPE(field_t), INTENT(inout) :: kantew
        INTEGER(intk), INTENT(in) :: triau(*)
        INTEGER(intk), INTENT(in) :: triav(*)
        INTEGER(intk), INTENT(in) :: triaw(*)

        ! Local variables
        INTEGER(intk) :: i, igrid, kk, jj, ii, ip3, ip3n
        REAL(realk), POINTER, CONTIGUOUS :: xstag(:), ystag(:), zstag(:), &
            ddx(:), ddy(:), ddz(:)

        DO i = 1, nmygridslvl(ilevel)
            igrid = mygridslvl(i, ilevel)

            CALL get_fieldptr(xstag, "XSTAG", igrid)
            CALL get_fieldptr(ystag, "YSTAG", igrid)
            CALL get_fieldptr(zstag, "ZSTAG", igrid)

            CALL get_fieldptr(ddx, "DDX", igrid)
            CALL get_fieldptr(ddy, "DDY", igrid)
            CALL get_fieldptr(ddz, "DDZ", igrid)

            CALL get_mgdims(kk, jj, ii, igrid)
            CALL get_ip3(ip3, igrid)
            CALL get_ip3n(ip3n, ntrimax, igrid)
            CALL freekante_grid(kk, jj, ii, xstag, ystag, zstag, &
                ddx, ddy, ddz, topol%n, topol%topol, &
                ntrimax, triau(ip3n), triav(ip3n), triaw(ip3n), &
                kanteu%arr(ip3), kantev%arr(ip3), kantew%arr(ip3), &
                knoten%arr(ip3))
        END DO
    END SUBROUTINE freekante_level


    SUBROUTINE freekante_grid(kk, jj, ii, xstag, ystag, zstag, &
            ddx, ddy, ddz, ntopol, topol, &
            ntrimax, triau, triav, triaw, &
            kanteu, kantev, kantew, knoten)
        ! Subroutine arguments
        INTEGER(intk), INTENT(in) :: kk, jj, ii
        REAL(realk), INTENT(in) :: xstag(ii), ystag(jj), zstag(kk)
        REAL(realk), INTENT(in) :: ddx(ii), ddy(jj), ddz(kk)
        INTEGER(intk), INTENT(in) :: ntopol
        REAL(realk), INTENT(in) :: topol(3, 3, ntopol)
        INTEGER(intk), INTENT(in) :: ntrimax
        INTEGER(intk), INTENT(in) :: triau(ntrimax, kk, jj, ii), &
            triav(ntrimax, kk, jj, ii), triaw(ntrimax, kk, jj, ii)
        REAL(realk), INTENT(inout) :: kanteu(kk, jj, ii), kantev(kk, jj, ii), &
            kantew(kk, jj, ii)
        REAL(realk), INTENT(inout) :: knoten(kk, jj, ii)

        ! Local variables
        INTEGER(intk) :: i, j, k, direction
        REAL(realk) :: x1, x2
        LOGICAL :: closed

        ! Schnittpunkte zwischen zwei geschlossenen Knoten werden geloescht
        !            (Kante ungeschnitten)
        ! Wenn eine geschnittene Kante zwischen zwei offenen Knoten liegt,
        !            wird jener mit den meisten geschnittenen
        !            oder geschlossenen Kanten in
        !            der Nachbarschaft geschlossen, falls beide dieselbe
        !            Anzahl haben, werden diese geschlossen, falls sie nicht
        !            aussenrum ueber offene Kanten und Knoten verbinden sind
        DO i = 2, ii-1
            DO j = 2, jj-1
                DO k = 2, kk-1
                    IF (NINT(knoten(k, j, i-1)) == NINT(knoten(k, j, i))) THEN
                        IF (NINT(knoten(k, j, i)) == 0) THEN
                            kanteu(k, j, i) = 1.0
                        ELSE
                            IF (kanteu(k, j, i) < 0.5) THEN
                                ! zwei offen Knoten mit geschnittener Kante
                                ! dazwischen
                                CALL checksurroundingedges(k, j, i, &
                                    kk, jj, ii, 'X', &
                                    kanteu, kantev, kantew, knoten, closed)
                                IF (closed .EQV. .FALSE.) THEN
                                    direction = 1
                                    CALL punkteeinekante2(k, j, i, &
                                        kk, jj, ii, xstag, ystag, zstag, &
                                        ddx, ddy, ddz, ntrimax, &
                                        ntopol, topol, &
                                        kanteu, kantev, kantew, &
                                        triau, triav, triaw, &
                                        direction, x1, x2)
                                    IF (x1-xstag(i-1) < xstag(i)-x2) THEN
                                        ! naechster Schnittpunkt bei knoten 1
                                        knoten(k, j, i-1) = 0.0
                                    ELSE IF (x1-xstag(i-1) > xstag(i)-x2) THEN
                                        ! naechster Schnittpunkt bei knoten 2
                                        knoten(k, j, i) = 0.0
                                    ELSE
                                        knoten(k, j, i) = 0.0
                                        knoten(k, j, i-1) = 0.0
                                    END IF
                                END IF
                            END IF
                        END IF
                    END IF
                END DO
            END DO
        END DO

        DO i = 2, ii-1
            DO j = 2, jj-1
                DO k = 2, kk-1
                    IF (NINT(knoten(k, j-1, i)) == NINT(knoten(k, j, i))) THEN
                        IF (NINT(knoten(k, j, i)) == 0) THEN
                            kantev(k, j, i) = 1.0
                        ELSE
                            IF (kantev(k, j, i) < 0.5) THEN
                                ! zwei offen Knoten mit geschnittener Kante
                                ! dazwischen
                                CALL checksurroundingedges(k, j, i, &
                                    kk, jj, ii, 'Y', &
                                    kanteu, kantev, kantew, knoten, closed)
                                IF (closed .EQV. .FALSE.) THEN
                                    direction = 2
                                    CALL punkteeinekante2(k, j, i, &
                                        kk, jj, ii, xstag, ystag, zstag, &
                                        ddx, ddy, ddz, ntrimax, &
                                        ntopol, topol, &
                                        kanteu, kantev, kantew, &
                                        triau, triav, triaw, &
                                        direction, x1, x2)
                                    IF (x1-ystag(j-1) < ystag(j)-x2) THEN
                                        ! naechster Schnittpunkt bei knoten 1
                                        knoten(k, j-1, i) = 0.0
                                    ELSE IF (x1-ystag(j-1) > ystag(j)-x2) THEN
                                        ! naechster Schnittpunkt bei knoten 2
                                        knoten(k, j, i) = 0.0
                                    ELSE
                                        knoten(k, j, i) = 0.0
                                        knoten(k, j-1, i) = 0.0
                                    END IF
                                END IF
                            END IF
                        END IF
                    END IF
                END DO
            END DO
        END DO

        DO i = 2, ii-1
            DO j = 2, jj-1
                DO k = 2, kk-1
                    IF (NINT(knoten(k-1, j, i)) == NINT(knoten(k, j, i))) THEN
                        IF (NINT(knoten(k, j, i)) == 0) THEN
                            kantew(k, j, i) = 1.0
                        ELSE
                            IF (kantew(k, j, i) < 0.5) THEN
                                ! zwei offen Knoten mit geschnittener Kante
                                ! dazwischen
                                CALL checksurroundingedges(k, j, i, &
                                    kk, jj, ii, 'Z', &
                                    kanteu, kantev, kantew, knoten, closed)
                                IF (closed .EQV. .FALSE.) THEN
                                    direction = 3
                                    CALL punkteeinekante2(k, j, i, &
                                        kk, jj, ii, xstag, ystag, zstag, &
                                        ddx, ddy, ddz, ntrimax, &
                                        ntopol, topol, &
                                        kanteu, kantev, kantew, &
                                        triau, triav, triaw, &
                                        direction, x1, x2)
                                    IF (x1-zstag(k-1) < zstag(k)-x2) THEN
                                        ! naechster Schnittpunkt bei knoten 1
                                        knoten(k-1, j, i) = 0.0
                                    ELSE IF (x1-zstag(k-1) > zstag(k)-x2) THEN
                                        ! naechster Schnittpunkt bei knoten 2
                                        knoten(k, j, i) = 0.0
                                    ELSE
                                        knoten(k, j, i) = 0.0
                                        knoten(k-1, j, i) = 0.0
                                    END IF
                                END IF
                            END IF
                        END IF
                    END IF
                END DO
            END DO
        END DO
    END SUBROUTINE freekante_grid


    SUBROUTINE checksurroundingedges(k, j, i, kk, jj, ii, cdir, &
            kanteu, kantev, kantew, knoten, closed)

        ! Subroutine arguments
        INTEGER(intk), INTENT(in) :: k, j, i
        INTEGER(intk), INTENT(in) :: kk, jj, ii
        CHARACTER, INTENT(in) :: cdir
        REAL(realk), INTENT(in) :: kanteu(kk, jj, ii), kantev(kk, jj, ii), &
            kantew(kk, jj, ii)
        REAL(realk), INTENT(inout) :: knoten(kk, jj, ii)
        LOGICAL, INTENT(out) :: closed

        ! Local variables
        REAL(realk) :: node1, node2

        closed = .TRUE.
        node1 = 0.0
        node2 = 0.0

        IF (cdir == 'X') THEN
            node1 = node1 + MAX(0.0_realk, kanteu(k, j, i+1)) &
                + MAX(0.0_realk, kantev(k, j, i)) &
                + MAX(0.0_realk, kantev(k, j+1, i)) &
                + MAX(0.0_realk, kantew(k, j, i)) &
                + MAX(0.0_realk, kantew(k+1, j, i))
            node2 = node2 + MAX(0.0_realk, kanteu(k, j, i-1)) &
                + MAX(0.0_realk, kantev(k, j, i-1)) &
                + MAX(0.0_realk, kantev(k, j+1, i-1)) &
                + MAX(0.0_realk, kantew(k, j, i-1)) &
                + MAX(0.0_realk, kantew(k+1, j, i-1))

            IF (NINT(node1) > NINT(node2) .AND. &
                    NINT(MIN(node1, node2)) <= 2) THEN
                knoten(k, j, i-1) = 0.0
            ELSE IF (NINT(node2) > NINT(node1) .AND. &
                    NINT(MIN(node1, node2)) <= 2) THEN
                knoten(k, j, i) = 0.0
            ELSE
                closed = .FALSE.
            END IF
        ELSE IF (cdir == 'Y') THEN
            node1 = node1 + MAX(0.0_realk, kantev(k, j+1, i)) &
                + MAX(0.0_realk, kanteu(k, j, i)) &
                + MAX(0.0_realk, kanteu(k, j, i+1)) &
                + MAX(0.0_realk, kantew(k, j, i)) &
                + MAX(0.0_realk, kantew(k+1, j, i))
            node2 = node2 + MAX(0.0_realk,kantev(k, j-1, i)) &
                + MAX(0.0_realk, kanteu(k, j-1, i)) &
                + MAX(0.0_realk, kanteu(k, j-1, i+1)) &
                + MAX(0.0_realk, kantew(k, j-1, i)) &
                + MAX(0.0_realk, kantew(k+1, j-1, i))

            IF (NINT(node1) > NINT(node2) .AND. &
                    NINT(MIN(node1, node2)) <= 2) THEN
                knoten(k, j-1, i) = 0.0
            ELSE IF (NINT(node2) > NINT(node1) .AND. &
                    NINT(MIN(node1, node2)) <= 2) THEN
                knoten(k, j, i) = 0.0
            ELSE
                closed = .FALSE.
            END IF
        ELSE IF (cdir == 'Z') THEN
            node1 = node1 + MAX(0.0_realk, kantew(k+1, j, i)) &
                + MAX(0.0_realk, kanteu(k, j, i)) &
                + MAX(0.0_realk, kanteu(k, j, i+1)) &
                + MAX(0.0_realk, kantev(k, j, i)) &
                + MAX(0.0_realk, kantev(k, j+1, i))

            node2 = node2 + MAX(0.0_realk, kantew(k-1, j, i)) &
                + MAX(0.0_realk, kanteu(k-1, j, i)) &
                + MAX(0.0_realk, kanteu(k-1, j, i+1)) &
                + MAX(0.0_realk, kantev(k-1, j, i)) &
                + MAX(0.0_realk, kantev(k-1, j+1, i))

            IF (NINT(node1) > NINT(node2) .AND. &
                    NINT(MIN(node1, node2)) <= 2) THEN
                knoten(k-1, j, i) = 0.0
            ELSE IF (NINT(node2) > NINT(node1) .AND. &
                    NINT(MIN(node1, node2)) <= 2) THEN
                knoten(k, j, i) = 0.0
            ELSE
                closed = .FALSE.
            END IF
        ELSE
            CALL  errr(__FILE__, __LINE__)
        END IF

    END SUBROUTINE checksurroundingedges
END MODULE freekante_mod
