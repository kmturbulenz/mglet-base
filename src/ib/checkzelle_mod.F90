MODULE checkzelle_mod
    USE core_mod, ONLY: realk, intk, int32, mygrids, nmygrids, &
        mygridslvl, nmygridslvl, minlevel, maxlevel, errr, connect, &
        field_t, get_mgdims, get_ip3, get_ip3n, get_fieldptr, &
        ind2sub
    USE cutcorner_mod, ONLY: calcint
    USE blockcheck_mod, ONLY: blockcheck_grid
    USE punktekoordinaten_mod, ONLY: punkteeinekante2
    USE topol_mod, ONLY: topol_t

    IMPLICIT NONE(type, external)
    PRIVATE

    PUBLIC :: checkzelle, checkzelle_grid

CONTAINS
    SUBROUTINE checkzelle(topol, ntrimax, knoten, kanteu, kantev, &
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
            CALL checkzelle_level(ilevel, topol, ntrimax, knoten, kanteu, &
                kantev, kantew, triau, triav, triaw)
        END DO
    END SUBROUTINE checkzelle


    SUBROUTINE checkzelle_level(ilevel, topol, ntrimax, knoten, kanteu, &
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
        INTEGER(intk) :: i, igrid, kk, jj, ii, ip3, ip3n, iloop, suminfo
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

            DO iloop = 1, 20
                CALL checkzelle_grid(kk, jj, ii, kanteu%arr(ip3), &
                    kantev%arr(ip3), kantew%arr(ip3), knoten%arr(ip3), &
                    ntrimax, triau(ip3n), triav(ip3n), triaw(ip3n), &
                    xstag, ystag, zstag, ddx, ddy, ddz, topol%n, &
                    topol%topol, suminfo)
                IF (suminfo == 0) EXIT
            END DO

            IF (suminfo /= 0) THEN
                WRITE(*, *) 'checkzelle suminfo', sumInfo, igrid
                CALL errr(__FILE__, __LINE__)
            END IF
        END DO
    END SUBROUTINE checkzelle_level


    SUBROUTINE checkzelle_grid(kk, jj, ii, kanteu, kantev, kantew, knoten, &
            ntrimax, triau, triav, triaw, xstag, ystag, zstag, ddx, ddy, ddz, &
            ntopol, topol, suminfo)
        ! Subroutine arguments
        INTEGER(intk), INTENT(in) :: kk, jj, ii
        REAL(realk), INTENT(inout) :: kanteu(kk, jj, ii), kantev(kk, jj, ii), &
            kantew(kk, jj, ii)
        REAL(realk), INTENT(inout) :: knoten(kk, jj, ii)
        INTEGER(intk), INTENT(in) :: ntrimax
        INTEGER(intk), INTENT(in) :: triau(ntrimax, kk, jj, ii), &
            triav(ntrimax, kk, jj, ii), triaw(ntrimax, kk, jj, ii)
        REAL(realk), INTENT(in) :: xstag(ii), ystag(jj), zstag(kk)
        REAL(realk), INTENT(in) :: ddx(ii), ddy(jj), ddz(kk)
        INTEGER(intk), INTENT(in) :: ntopol
        REAL(realk), INTENT(in) :: topol(3, 3, ntopol)
        INTEGER(intk), INTENT(out) :: suminfo

        ! Local variables
        INTEGER(intk) :: i, j, k
        INTEGER(intk) :: in, jn, kn
        INTEGER(intk) :: direction, info, n, nmin
        INTEGER(intk) :: knotenzelle(8)

        REAL(realk) :: x1, x2
        REAL(realk) :: volk(8), lx, ly, lz
        REAL(realk) :: vmin

        suminfo = 0

        DO i = 2, ii
            DO j = 2, jj
                DO k = 2, kk

                    innerloop: DO
                        ! Nummeriegung der Knoten s. CHECKVERB
                        knotenzelle(1) = NINT(knoten(k-1, j-1, i-1))
                        knotenzelle(2) = NINT(knoten(k-1, j-1, i))
                        knotenzelle(3) = NINT(knoten(k, j-1, i))
                        knotenzelle(4) = NINT(knoten(k, j-1, i-1))
                        knotenzelle(5) = NINT(knoten(k-1, j, i-1))
                        knotenzelle(6) = NINT(knoten(k-1, j, i))
                        knotenzelle(7) = NINT(knoten(k, j, i))
                        knotenzelle(8) = NINT(knoten(k, j, i-1))

                        ! Unklar: Was bedeutet knoten=2 ?
                        DO n = 1, 8
                            IF (knotenzelle(n) == 2) THEN
                                knotenzelle(n) = 0
                            END IF
                        END DO

                        CALL checkverb(knotenzelle, info)
                        suminfo = suminfo + info
                        ! info = 0 : offene/geschlossene Knoten
                        !            zusammenhaengend
                        ! info = 1 : offene Knoten nicht zusammenhaengend
                        ! info = 2 : geschlossene Knoten nicht zusammenhaengend
                        ! info = 3 : offene und geschlossene Knoten nicht
                        !            zusammenhaengend
                        IF (info == 0) EXIT

                        ! offene Knoten nicht zusammenhaengend
                        ! geschlossene Knoten zusammenhaengend oder auch nicht
                        ! => Knoten neben geschlossenen schliessen
                        !    danach nochmals pruefen

                        n = 0
                        DO in = i-1, i
                            DO jn = j-1, j
                                DO kn = k-1, k
                                    ! n: Knotenindex, Nummerierung nicht wie in
                                    ! CHECKVERB
                                    n = n + 1
                                    knotenzelle(n) =  NINT(knoten(kn,jn,in))

                                    ! volk: jedem Knoten zugeordnetes Volumen
                                    IF (knotenzelle(n) == 0) THEN
                                        volk(n) = 0.0
                                        CYCLE
                                    END IF

                                    ! freie Laenge in x-Richtung
                                    direction = 1
                                    IF (kanteu(kn, jn, in) < 0.5) THEN
                                        CALL punkteeinekante2(kn, jn, in, &
                                            kk, jj, ii, &
                                            xstag, ystag, zstag, &
                                            ddx, ddy, ddz, &
                                            ntrimax, ntopol, topol, &
                                            kanteu, kantev, kantew,&
                                            triau, triav ,triaw, &
                                            direction, x1, x2)
                                        lx = xstag(in) - x2
                                    ELSE
                                        lx = ddx(in)
                                    END IF
                                    IF (in < ii .AND. kanteu(kn, jn, &
                                            MIN(in+1, ii)) < 0.5) THEN
                                        CALL punkteeinekante2(kn, jn, in+1, &
                                            kk, jj, ii, &
                                            xstag, ystag, zstag, &
                                            ddx, ddy, ddz, &
                                            ntrimax, ntopol, topol, &
                                            kanteu, kantev, kantew,&
                                            triau, triav ,triaw, &
                                            direction, x1, x2)
                                        lx = lx + x1 - xstag(in)
                                    ELSE
                                        lx = lx + ddx(MIN(in+1, ii))
                                    END IF

                                    ! freie Laenge in y-Richtung
                                    direction = 2
                                    IF (kantev(kn, jn, in) < 0.5) THEN
                                        CALL punkteeinekante2(kn, jn, in, &
                                            kk, jj, ii, &
                                            xstag, ystag, zstag, &
                                            ddx, ddy, ddz, &
                                            ntrimax, ntopol, topol, &
                                            kanteu, kantev, kantew,&
                                            triau, triav ,triaw, &
                                            direction, x1, x2)
                                        ly = ystag(jn) - x2
                                    ELSE
                                        ly = ddy(jn)
                                    END IF
                                    IF (jn < jj .AND. kantew(kn, &
                                            MIN(jn+1, jj), in) < 0.5) THEN
                                        CALL punkteeinekante2(kn, jn+1, in, &
                                            kk, jj, ii, &
                                            xstag, ystag, zstag, &
                                            ddx, ddy, ddz, &
                                            ntrimax, ntopol, topol, &
                                            kanteu, kantev, kantew,&
                                            triau, triav ,triaw, &
                                            direction, x1, x2)
                                        ly = ly + x1 - ystag(jn)
                                    ELSE
                                        ly = ly + ddy(MIN(jn+1, jj))
                                    END IF


                                    ! freie Laenge in z-Richtung
                                    direction = 3
                                    IF (kantew(kn, jn, in) < 0.5) THEN
                                        CALL punkteeinekante2(kn, jn, in, &
                                            kk, jj, ii, &
                                            xstag, ystag, zstag, &
                                            ddx, ddy, ddz, &
                                            ntrimax, ntopol, topol, &
                                            kanteu, kantev, kantew,&
                                            triau, triav ,triaw, &
                                            direction, x1, x2)
                                        lz = zstag(kn) - x2
                                    ELSE
                                        lz = ddz(kn)
                                    END IF
                                    IF (kn < kk .AND. kantew(MIN(kn+1, kk), &
                                            jn, in) < 0.5) THEN
                                        CALL punkteeinekante2(kn+1, jn, in, &
                                            kk, jj, ii, &
                                            xstag, ystag, zstag, &
                                            ddx, ddy, ddz, &
                                            ntrimax, ntopol, topol, &
                                            kanteu, kantev, kantew,&
                                            triau, triav ,triaw, &
                                            direction, x1, x2)
                                        lz = lz + x1 - zstag(kn)
                                    ELSE
                                        lz = lz + ddz(MIN(kn+1, kk))
                                    END IF

                                    volk(n) = lx*ly*lz
                                END DO
                            END DO
                        END DO

                        nmin = 0
                        vmin = 0.0
                        DO n = 1, 8
                            IF (volk(n) /= 0.0) THEN
                                vmin = volk(n)
                                EXIT
                            END IF
                        END DO

                        DO n = 1, 8
                            IF (volk(n) /= 0.0 .AND. volk(n) <= vmin) THEN
                                vmin = volk(n)
                                nmin = n
                            END IF
                        END DO

                        IF (nmin == 0) THEN
                            IF (vmin == 0.0) THEN
                                knoten(k-1:k, j-1:j, i-1:i) = 0.0
                                WRITE(*, *) 'CHECKZELLE A.5 all', knotenzelle
                                WRITE(*, *) 'CHECKZELLE A.6 all', k, j, i
                            ELSE
                                CALL errr(__FILE__, __LINE__)
                            END IF
                        ELSE
                            CALL ind2sub(nmin, kn, jn, in, 2, 2, 2)
                            knoten(k-1+kn-1, j-1+jn-1, i-1+in-1) = 0.0
                        END IF
                    END DO innerloop

                END DO
            END DO
        END DO

        DO i = 1, ii
            DO j = 1, jj
                DO k = 1, kk
                    knoten(k, j, i) = MIN(knoten(k, j, i), 1.0_realk)
                END DO
            END DO
        END DO

        CALL blockcheck_grid(kk, jj, ii, kanteu, kantev, kantew, knoten, 1)
    END SUBROUTINE checkzelle_grid


    PURE SUBROUTINE checkverb(knotenzelle, info)
        ! bekommt Knoten,
        !  (knoten(i) = 0: Koerper, = 1: Fluid)
        ! liefert Information, ob die offenen/geschlossenen Knoten ueber
        !   Kanten verbunden sind:
        !   info = 0 : offene/geschlossene Knoten zusammenhaengend
        !   info = 1 : offene Knoten nicht zusammenhaengend
        !   info = 2 : geschlossene Knoten nicht zusammenhaengend
        !   info = 3 : offene und geschlossene Knoten nicht zusammenhaengend
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
        !

        INTEGER(intk), INTENT(in) :: knotenzelle(8)
        INTEGER(intk), INTENT(out) :: info

        INTEGER(intk) :: kantenenden(2, 12)
        INTEGER(intk) :: nknoten, nkanten0, nkanten1
        INTEGER(intk) :: ii

        info = 0

        ! Knotennummern der Kantenenden
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
        DO ii = 1, 8
            nknoten = nknoten + knotenzelle(ii)
        END DO

        nkanten0 = 0
        nkanten1 = 0
        DO ii = 1,12
            IF (knotenzelle(kantenenden(1, ii)) == &
                    knotenzelle(kantenenden(2, ii))) THEN
                IF (knotenzelle(kantenenden(1, ii)) == 1) THEN
                    nkanten1 = nkanten1 + 1
                ELSE IF (knotenzelle(kantenenden(1, ii)) == 0) then
                    nkanten0 = nkanten0 + 1
                END IF
            END IF
        END DO

        IF (nknoten > nkanten1+1) THEN
            info = info + 1
        END IF

        IF (8-nknoten > nkanten0+1) THEN
            info = info + 2
        END IF
    END SUBROUTINE checkverb
END MODULE checkzelle_mod
