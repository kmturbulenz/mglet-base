MODULE punktekoordinaten_mod
    USE core_mod, ONLY: realk, intk, errr
    USE cutcorner_mod, ONLY: calcint

    IMPLICIT NONE(type, external)
    PRIVATE

    PUBLIC :: punktekoordinaten, punkteeinekante, punkteeinekante2

CONTAINS
    SUBROUTINE punktekoordinaten(k, j, i, kk, jj, ii, xstag, ystag, zstag, &
            ntopol, topol, topolbodyid, ntrimax, triau, triav, triaw, &
            knoten, kanteu, kantev, kantew, npoints, npmax, xx, bodyidkanten)

        ! Liefert fuer eine Zelle
        !         xx(1:3,npoints+1:npoints+20)
        !              Koordinaten der Eckpunkte (8)
        !              Koordinaten der Kantenschnittpunkte (12)
        !                     kanteu neg. : Schnittpunkt, der offenem Knoten am
        !                                   nächsten ist
        !                     kanteu = 0  : Schnittpunkt wird kuenstlich gelegt
        !                                   je nach _NEWCUTMITTE_
        !                                 : bei #undef _NEWCUTMITTE_ wird durch
        !                                   den defualt algo
        !                                : der Schnittpunkt auf den
        !                                  geschlossenen Knoten gelegt
        !                     kanteu = 1  : Koordinate des Kantenendes
        !         npoints = npoints + 20

        ! Subroutine arguments
        INTEGER(intk), INTENT(in) :: k, j, i
        INTEGER(intk), INTENT(in) :: kk, jj, ii
        REAL(realk), INTENT(in) :: xstag(ii), ystag(jj), zstag(kk)
        INTEGER(intk), INTENT(in) :: ntopol
        REAL(realk), INTENT(in) :: topol(3, 3, ntopol)
        INTEGER(intk), INTENT(in) :: topolbodyid(ntopol)
        INTEGER(intk), INTENT(in) :: ntrimax
        INTEGER(intk), INTENT(in) :: triau(ntrimax, kk, jj, ii), &
            triav(ntrimax, kk, jj, ii), triaw(ntrimax, kk, jj, ii)
        REAL(realk), INTENT(in) :: knoten(kk, jj, ii)
        REAL(realk), INTENT(in) :: kanteu(kk, jj, ii), kantev(kk, jj, ii), &
            kantew(kk, jj, ii)
        INTEGER(intk), INTENT(inout) :: npoints
        INTEGER(intk), INTENT(in) :: npmax
        REAL(realk), INTENT(inout) :: xx(3, npmax)
        INTEGER(intk), intent(out) :: bodyidkanten(12)

        ! Local variables
        INTEGER(intk) :: direction, ikante, ip, jp, kp

        ! Zeroize INTENT(out) variables
        ! ukanten = 0.0
        bodyidkanten = 0
        ! #ifdef _SCABODY_
        ! scakanten = 0.0
        ! #endif

        ! Schnittpunkte auf den Kanten
        ikante = 1
        direction = 1
        ! Startpunkt der Kante in der Zelle
        ip = 0-1
        jp = 0-1
        kp = 0-1
        CALL punkteeinekante(k, j, i, kp, jp, ip, kk, jj, ii, &
            xstag, ystag, zstag, direction, ikante, ntopol, topol, &
            topolbodyid, ntrimax, triau, triav, triaw, knoten, &
            kanteu, kantev, kantew, npoints, npmax, xx, bodyidkanten)

        ikante = 2
        direction = 3
        ip = 1-1
        jp = 0-1
        kp = 0-1
        CALL punkteeinekante(k, j, i, kp, jp, ip, kk, jj, ii, &
            xstag, ystag, zstag, direction, ikante, ntopol, topol, &
            topolbodyid, ntrimax, triau, triav, triaw, knoten, &
            kanteu, kantev, kantew, npoints, npmax, xx, bodyidkanten)

        ikante = 3
        direction = 1
        ip = 0-1
        jp = 0-1
        kp = 1-1
        CALL punkteeinekante(k, j, i, kp, jp, ip, kk, jj, ii, &
            xstag, ystag, zstag, direction, ikante, ntopol, topol, &
            topolbodyid, ntrimax, triau, triav, triaw, knoten, &
            kanteu, kantev, kantew, npoints, npmax, xx, bodyidkanten)

        ikante = 4
        direction = 3
        ip = 0-1
        jp = 0-1
        kp = 0-1
        CALL punkteeinekante(k, j, i, kp, jp, ip, kk, jj, ii, &
            xstag, ystag, zstag, direction, ikante, ntopol, topol, &
            topolbodyid, ntrimax, triau, triav, triaw, knoten, &
            kanteu, kantev, kantew, npoints, npmax, xx, bodyidkanten)

        ikante = 5
        direction = 1
        ip = 0-1
        jp = 1-1
        kp = 0-1
        CALL punkteeinekante(k, j, i, kp, jp, ip, kk, jj, ii, &
            xstag, ystag, zstag, direction, ikante, ntopol, topol, &
            topolbodyid, ntrimax, triau, triav, triaw, knoten, &
            kanteu, kantev, kantew, npoints, npmax, xx, bodyidkanten)

        ikante = 6
        direction = 3
        ip = 1-1
        jp = 1-1
        kp = 0-1
        CALL punkteeinekante(k, j, i, kp, jp, ip, kk, jj, ii, &
            xstag, ystag, zstag, direction, ikante, ntopol, topol, &
            topolbodyid, ntrimax, triau, triav, triaw, knoten, &
            kanteu, kantev, kantew, npoints, npmax, xx, bodyidkanten)

        ikante = 7
        direction = 1
        ip = 0-1
        jp = 1-1
        kp = 1-1
        CALL punkteeinekante(k, j, i, kp, jp, ip, kk, jj, ii, &
            xstag, ystag, zstag, direction, ikante, ntopol, topol, &
            topolbodyid, ntrimax, triau, triav, triaw, knoten, &
            kanteu, kantev, kantew, npoints, npmax, xx, bodyidkanten)

        ikante = 8
        direction = 3
        ip = 0-1
        jp = 1-1
        kp = 0-1
        CALL punkteeinekante(k, j, i, kp, jp, ip, kk, jj, ii, &
            xstag, ystag, zstag, direction, ikante, ntopol, topol, &
            topolbodyid, ntrimax, triau, triav, triaw, knoten, &
            kanteu, kantev, kantew, npoints, npmax, xx, bodyidkanten)

        ikante = 9
        direction = 2
        ip = 0-1
        jp = 0-1
        kp = 0-1
        CALL punkteeinekante(k, j, i, kp, jp, ip, kk, jj, ii, &
            xstag, ystag, zstag, direction, ikante, ntopol, topol, &
            topolbodyid, ntrimax, triau, triav, triaw, knoten, &
            kanteu, kantev, kantew, npoints, npmax, xx, bodyidkanten)

        ikante = 10
        direction = 2
        ip = 1-1
        jp = 0-1
        kp = 0-1
        CALL punkteeinekante(k, j, i, kp, jp, ip, kk, jj, ii, &
            xstag, ystag, zstag, direction, ikante, ntopol, topol, &
            topolbodyid, ntrimax, triau, triav, triaw, knoten, &
            kanteu, kantev, kantew, npoints, npmax, xx, bodyidkanten)

        ikante = 11
        direction = 2
        ip = 1-1
        jp = 0-1
        kp = 1-1
        CALL punkteeinekante(k, j, i, kp, jp, ip, kk, jj, ii, &
            xstag, ystag, zstag, direction, ikante, ntopol, topol, &
            topolbodyid, ntrimax, triau, triav, triaw, knoten, &
            kanteu, kantev, kantew, npoints, npmax, xx, bodyidkanten)

        ikante = 12
        direction = 2
        ip = 0-1
        jp = 0-1
        kp = 1-1
        CALL punkteeinekante(k, j, i, kp, jp, ip, kk, jj, ii, &
            xstag, ystag, zstag, direction, ikante, ntopol, topol, &
            topolbodyid, ntrimax, triau, triav, triaw, knoten, &
            kanteu, kantev, kantew, npoints, npmax, xx, bodyidkanten)

        ! Knoten an der Quaderecken
        xx(1, npoints+1) = xstag(i-1)
        xx(1, npoints+2) = xstag(i)
        xx(1, npoints+3) = xstag(i)
        xx(1, npoints+4) = xstag(i-1)
        xx(1, npoints+5) = xstag(i-1)
        xx(1, npoints+6) = xstag(i)
        xx(1, npoints+7) = xstag(i)
        xx(1, npoints+8) = xstag(i-1)

        xx(2, npoints+1) = ystag(j-1)
        xx(2, npoints+2) = ystag(j-1)
        xx(2, npoints+3) = ystag(j-1)
        xx(2, npoints+4) = ystag(j-1)
        xx(2, npoints+5) = ystag(j)
        xx(2, npoints+6) = ystag(j)
        xx(2, npoints+7) = ystag(j)
        xx(2, npoints+8) = ystag(j)

        xx(3, npoints+1) = zstag(k-1)
        xx(3, npoints+2) = zstag(k-1)
        xx(3, npoints+3) = zstag(k)
        xx(3, npoints+4) = zstag(k)
        xx(3, npoints+5) = zstag(k-1)
        xx(3, npoints+6) = zstag(k-1)
        xx(3, npoints+7) = zstag(k)
        xx(3, npoints+8) = zstag(k)

        npoints = npoints + 20
    END SUBROUTINE punktekoordinaten


    SUBROUTINE punkteeinekante(k, j, i, kp, jp, ip, kk, jj, ii, &
            xstag, ystag, zstag, direction, ikante, &
            ntopol, topol, topolbodyid, &
            ntrimax, triau, triav, triaw, knoten, kanteu, kantev, kantew, &
            npoints, npmax, xx, bodyidkanten)

        ! Subroutine arguments
        INTEGER(intk), INTENT(in) :: k, j, i
        INTEGER(intk), INTENT(in) :: kp, jp, ip
        INTEGER(intk), INTENT(in) :: kk, jj, ii
        REAL(realk), INTENT(in) :: xstag(ii), ystag(jj), zstag(kk)
        INTEGER(intk), INTENT(in) :: direction, ikante
        INTEGER(intk), INTENT(in) :: ntopol
        REAL(realk), INTENT(in) :: topol(3, 3, ntopol)
        INTEGER(intk), INTENT(in) :: topolbodyid(ntopol)
        INTEGER(intk), INTENT(in) :: ntrimax
        INTEGER(intk), INTENT(in) :: triau(ntrimax, kk, jj, ii), &
            triav(ntrimax, kk, jj, ii), triaw(ntrimax, kk, jj, ii)
        REAL(realk), INTENT(in) :: knoten(kk, jj, ii)
        REAL(realk), INTENT(in) :: kanteu(kk, jj, ii), kantev(kk, jj, ii), &
            kantew(kk, jj, ii)
        INTEGER(intk), INTENT(in) :: npoints
        INTEGER(intk), INTENT(in) :: npmax
        REAL(realk), INTENT(inout) :: xx(3, npmax)
        INTEGER(intk), intent(inout) :: bodyidkanten(12)

        ! Local variables
        INTEGER(intk) :: idx, itri, ntri
        REAL(realk) :: xs, ys, zs
        REAL(realk) :: xp, yp, zp
        REAL(realk) :: xswahl, xstmp, xo

        ! Sanity check
        IF (npoints + ikante + 8 > npmax) THEN
            WRITE(*, *) "npoints, ikante, npmax: ", npoints, ikante, npmax
            CALL errr(__FILE__, __LINE__)
        END IF

        xx(1, npoints+ikante+8) = xstag(i+ip)
        xx(2, npoints+ikante+8) = ystag(j+jp)
        xx(3, npoints+ikante+8) = zstag(k+kp)

        ! xo : Koordinate des offenen Knotens
        ! xswahl: initialisiert mit Koordinate des geschlossenen Knotens
        SELECT CASE(direction)
        CASE (1)
            xo = xstag(i+ip+1-NINT(knoten(k+kp, j+jp, i+ip)))
            xswahl = xstag(i+ip+NINT(knoten(k+kp, j+jp, i+ip)))
            ntri = -NINT(kanteu(k+kp, j+jp, i))
        CASE (2)
            xo = ystag(j+jp+1-NINT(knoten(k+kp, j+jp, i+ip)))
            xswahl = ystag(j+jp  +NINT(knoten(k+kp, j+jp, i+ip)))
            ntri = -NINT(kantev(k+kp, j, i+ip))
        CASE (3)
            xo = zstag(k+kp+1-NINT(knoten(k+kp, j+jp, i+ip)))
            xswahl = zstag(k+kp  +NINT(knoten(k+kp, j+jp, i+ip)))
            ntri = -NINT(kantew(k, j+jp, i+ip))
        CASE DEFAULT
            CALL errr(__FILE__, __LINE__)
        END SELECT

        DO idx = 1, ntri
            xp = xstag(i+ip)
            yp = ystag(j+jp)
            zp = zstag(k+kp)

            SELECT CASE(direction)
            CASE (1)
                itri = triau(idx, k+kp, j+jp, i)
            CASE (2)
                itri = triav(idx, k+kp, j, i+ip)
            CASE (3)
                itri = triaw(idx, k, j+jp, i+ip)
            END SELECT

            IF (itri > 0) THEN
                CALL calcint(xp, yp, zp, xs, ys, zs, &
                    itri, ntopol, direction, topol)

                SELECT CASE(direction)
                CASE (1)
                    xstmp = xs
                CASE (2)
                    xstmp = ys
                CASE (3)
                    xstmp = zs
                END SELECT

                ! First iteration of loop
                IF (idx == 1) THEN
                    ! DO idir = 1, 3
                    !     ukanten(idir, ikante) = topol(4, idir, itri)
                    !     #ifdef _SCABODY_
                    !     scakanten(idir, ikante) = topol(5, idir, itri)
                    !     #endif
                    ! END DO
                    bodyidkanten(ikante) = topolbodyid(itri)
                END IF

                ! Subsequent iterations/corrections
                IF (ABS(xstmp-xo) < ABS(xswahl-xo)) THEN
                    xswahl = xstmp
                    ! DO idir = 1, 3
                    !     ukanten(idir, ikante) = topol(4, idir, itri)
                    !     #ifdef _SCABODY_
                    !     scakanten(idir, ikante) = topol(5, idir, itri)
                    !     #endif
                    ! END DO
                    bodyidkanten(ikante) = topolbodyid(itri)
                END IF
            ELSE IF (itri == -1) THEN
                xswahl = xo
            ELSE
                IF (idx > 1) THEN
                    WRITE(*, *) 'PUNKTEEINEKANTE Info: ein Dreieck gleich Null'
                    WRITE(*, *) 'PUNKTEEINEKANTE Info: ', i, j, k, &
                        ip, jp, kp, ntri, itri, idx
                ELSE
                    WRITE(*, *)'PUNKTEEINEKANTE err: ', i, j, k, &
                        ip, jp, kp, ntri, itri, idx
                    CALL  errr(__FILE__, __LINE__)
                END IF
            END IF
        END DO

        xx(direction, npoints+ikante+8) = xswahl
    END SUBROUTINE punkteeinekante


    SUBROUTINE punkteeinekante2(k, j, i, kk, jj, ii, xstag, ystag, zstag, &
            ddx, ddy, ddz, ntrimax, ntopol, topol, kanteu, kantev, kantew, &
            triau, triav, triaw, direction, x1, x2)

        ! Subroutine arguments
        INTEGER(intk), INTENT(in) :: k, j, i, kk, jj, ii
        REAL(realk), INTENT(in) :: xstag(ii), ystag(jj), zstag(kk)
        REAL(realk), INTENT(in) :: ddx(ii), ddy(jj), ddz(kk)
        INTEGER(intk), INTENT(in) :: ntrimax, ntopol
        REAL(realk), INTENT(in) :: topol(3, 3, ntopol)
        REAL(realk), INTENT(in) :: kanteu(kk, jj, ii), kantev(kk, jj, ii), &
            kantew(kk, jj, ii)
        INTEGER(intk), INTENT(in) :: triau(ntrimax, kk, jj, ii), &
            triav(ntrimax, kk, jj, ii), triaw(ntrimax, kk, jj, ii)
        INTEGER(intk), INTENT(in) :: direction
        REAL(realk), INTENT(out) :: x1, x2

        ! Local variables
        INTEGER(intk) :: idx, ntri, itri
        REAL(realk) :: xs, ys, zs, xstmp, xp, yp, zp

        ! x1,x2: erster / letzter schnittpunkt
        ! xp, yp, zp: startecke der kante
        xp = xstag(i)
        yp = ystag(j)
        zp = zstag(k)

        SELECT CASE(direction)
        CASE (1)
            x1 = xstag(i)
            x2 = xstag(i) - ddx(i)
            ntri = -NINT(kanteu(k, j, i))
            xp = x2
        CASE (2)
            x1 = ystag(j)
            x2 = ystag(j) - ddy(j)
            ntri = -NINT(kantev(k, j, i))
            yp = x2
        CASE (3)
            x1 = zstag(k)
            x2 = zstag(k) - ddz(k)
            ntri = -NINT(kantew(k, j, i))
            zp = x2
        CASE DEFAULT
            CALL errr(__FILE__, __LINE__)
        END SELECT

        ! xs : schnittpunkt
        DO idx = 1, ntri
            SELECT CASE(direction)
            CASE (1)
                itri = triau(idx, k, j, i)
            CASE (2)
                itri = triav(idx, k, j, i)
            CASE (3)
                itri = triaw(idx, k, j, i)
            END SELECT

            IF (itri > 0) THEN
                CALL calcint(xp, yp, zp, xs, ys, zs, &
                    itri, ntopol, direction, topol)
                SELECT CASE(direction)
                CASE (1)
                    xstmp = xs
                CASE (2)
                    xstmp = ys
                CASE (3)
                    xstmp = zs
                END SELECT
                x1 = MIN(xstmp, x1)
                x2 = MAX(xstmp, x2)
            ELSE IF (itri == -1) then
                CALL  errr(__FILE__, __LINE__)
            ELSE
                WRITE(*, *) 'punkteeinekante2 info: ein dreieck gleich null'
                WRITE(*, *) 'punkteeinekante2 err:  ', i, j, k, ntri, itri, &
                    idx
                IF (.NOT. (i == 1 .AND. direction == 1) .AND. .NOT. &
                        (j  == 1 .AND. direction == 2) .AND. .NOT. &
                        (k == 1.and. direction == 3)) THEN
                    CALL  errr(__FILE__, __LINE__)
                END IF
            END IF
        END DO
    END SUBROUTINE punkteeinekante2
END MODULE punktekoordinaten_mod
