MODULE calcauavaw_mod
    USE core_mod, ONLY: realk, intk, errr, minlevel, maxlevel, nmygrids, &
        mygrids, nmygridslvl, mygridslvl, get_mgdims, get_ip3, get_ip3n, &
        get_fieldptr, field_t, connect
    USE blockcheck_mod, ONLY: blockcheck_grid
    USE calcfacearea_mod, ONLY: calcfacedata, calcwallfacecenter, &
        calcwallfacecenterrescue
    USE flzelle_mod, ONLY: flzelle
    USE ibconst_mod, ONLY: maccur
    USE knotenundkanten_mod, ONLY: knotenundkanten
    USE punktekoordinaten_mod, ONLY: punktekoordinaten
    USE topol_mod, ONLY: topol_t

    IMPLICIT NONE(type, external)
    PRIVATE

    PUBLIC :: calcauavaw

CONTAINS
    SUBROUTINE calcauavaw(topol, ntrimax, triau, triav, triaw, &
            knoten, kanteu, kantev, kantew, bzelltyp, au, av, aw, &
            icells, icellspointer, xpsw, yus, zus, xvs, zvs, xws, yws)
        ! Subroutine arguments
        TYPE(topol_t), INTENT(in) :: topol
        INTEGER(intk), INTENT(in) :: ntrimax
        INTEGER(intk), INTENT(in) :: triau(*), triav(*), triaw(*)
        TYPE(field_t), INTENT(in) :: knoten
        TYPE(field_t), INTENT(inout) :: kanteu, kantev, kantew
        INTEGER(intk), INTENT(in) :: bzelltyp(*)
        TYPE(field_t), INTENT(inout) :: au, av, aw
        INTEGER(intk), INTENT(in), OPTIONAL :: icells(:)
        INTEGER(intk), INTENT(in), OPTIONAL :: icellspointer(:)
        REAL(realk), INTENT(out), CONTIGUOUS, OPTIONAL :: xpsw(:, :)
        REAL(realk), INTENT(out), OPTIONAL :: yus(*), zus(*), xvs(*), zvs(*), &
            xws(*), yws(*)

        ! Local variables
        INTEGER(intk) :: ilevel

        DO ilevel = minlevel, maxlevel
            CALL calcauavaw_level(ilevel, topol, ntrimax, triau, triav, triaw, &
                knoten, kanteu, kantev, kantew, bzelltyp, au, av, aw, &
                icells, icellspointer, xpsw, yus, zus, xvs, zvs, xws, yws)
        END DO
    END SUBROUTINE calcauavaw


    SUBROUTINE calcauavaw_level(ilevel, topol, ntrimax, triau, triav, triaw, &
            knoten, kanteu, kantev, kantew, bzelltyp, au, av, aw, &
            icells, icellspointer, xpsw, &
            yus, zus, xvs, zvs, xws, yws)
        ! Subroutine arguments
        INTEGER(intk), INTENT(in) :: ilevel
        TYPE(topol_t), INTENT(in) :: topol
        INTEGER(intk), INTENT(in) :: ntrimax
        INTEGER(intk), INTENT(in) :: triau(*), triav(*), triaw(*)
        TYPE(field_t), INTENT(in) :: knoten
        TYPE(field_t), INTENT(inout) :: kanteu, kantev, kantew
        INTEGER(intk), INTENT(in) :: bzelltyp(*)
        TYPE(field_t), INTENT(inout) :: au, av, aw
        INTEGER(intk), INTENT(in), OPTIONAL :: icells(:)
        INTEGER(intk), INTENT(in), OPTIONAL :: icellspointer(:)
        REAL(realk), INTENT(out), OPTIONAL, CONTIGUOUS, TARGET :: xpsw(:, :)
        REAL(realk), INTENT(out), OPTIONAL :: yus(*), zus(*), xvs(*), zvs(*), &
            xws(*), yws(*)

        ! Local variables
        INTEGER(intk) :: i, igrid, kk, jj, ii, ip3, ip3n, ipp, ncells
        REAL(realk), POINTER, CONTIGUOUS :: xstag(:), ystag(:), zstag(:), &
            ddx(:), ddy(:), ddz(:)
        REAL(realk), POINTER, CONTIGUOUS :: xpsw_p(:, :)
        LOGICAL :: has_yus

        IF (PRESENT(yus) .AND. PRESENT(zus) .AND. PRESENT(xvs) .AND. &
                PRESENT(zvs) .AND. PRESENT(xws) .AND. PRESENT(yws)) THEN
            ! All optional arguments are present - fine!
            has_yus = .TRUE.
        ELSE IF (PRESENT(yus) .OR. PRESENT(zus) .OR. PRESENT(xvs) .OR. &
                PRESENT(zvs) .OR. PRESENT(xws) .OR. PRESENT(yws)) THEN
            ! One or more optional arguments are present but not all - wrong!
            CALL errr(__FILE__, __LINE__)
        ELSE
            ! No optional arguments
            has_yus = .FALSE.
        END IF

        ! ix xpsw is present, we also need icellspointer and icells
        IF (PRESENT(xpsw)) THEN
            IF (.NOT. PRESENT(icellspointer)) CALL errr(__FILE__, __LINE__)
            IF (.NOT. PRESENT(icells)) CALL errr(__FILE__, __LINE__)
        END IF

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

            IF (PRESENT(xpsw)) THEN
                ipp = icellspointer(igrid)
                ncells = icells(igrid)
                xpsw_p => xpsw(:, ipp:ipp+ncells-1)
            ELSE
                NULLIFY(xpsw_p)
            END IF

            ! Ugly code, because I'm not allowed to do yus(ip3) if yus
            ! is not present, therefore I must do this
            IF (has_yus) THEN
                CALL calcauavaw_grid(kk, jj, ii, xstag, ystag, zstag, &
                    ddx, ddy, ddz, topol%n, topol%topol, topol%bodyid, &
                    ntrimax, triau(ip3n), triav(ip3n), triaw(ip3n), &
                    knoten%arr(ip3), kanteu%arr(ip3), kantev%arr(ip3), &
                    kantew%arr(ip3), bzelltyp(ip3), au%arr(ip3), &
                    av%arr(ip3), aw%arr(ip3), xpsw_p, &
                    yus(ip3), zus(ip3), xvs(ip3), zvs(ip3), xws(ip3), yws(ip3))
            ELSE
                CALL calcauavaw_grid(kk, jj, ii, xstag, ystag, zstag, &
                    ddx, ddy, ddz, topol%n, topol%topol, topol%bodyid, &
                    ntrimax, triau(ip3n), triav(ip3n), triaw(ip3n), &
                    knoten%arr(ip3), kanteu%arr(ip3), kantev%arr(ip3), &
                    kantew%arr(ip3), bzelltyp(ip3), au%arr(ip3), av%arr(ip3), &
                    aw%arr(ip3), xpsw_p)
            END IF
        END DO

        ! Originally found on blockbp directly - now moved here becuase this is
        ! the last spot where AU, AV, AW are touched
        CALL connect(ilevel, 2, v1=au, v2=av, v3=aw, corners=.TRUE.)
    END SUBROUTINE calcauavaw_level


    SUBROUTINE calcauavaw_grid(kk, jj, ii, xstag, ystag, zstag, &
            ddx, ddy, ddz, ntopol, topol, topolbodyid, &
            ntrimax, triau, triav, triaw, knoten, kanteu, kantev, kantew, &
            bzelltyp, au, av, aw, xpsw, yus, zus, xvs, zvs, xws, yws)

        ! Subroutine arguments
        INTEGER(intk), INTENT(in) :: kk, jj, ii
        REAL(realk), INTENT(in) :: xstag(ii), ystag(jj), zstag(kk)
        REAL(realk), INTENT(in) :: ddx(ii), ddy(jj), ddz(kk)
        INTEGER(intk), INTENT(in) :: ntopol
        REAL(realk), INTENT(in) :: topol(3, 3, ntopol)
        INTEGER(intk), INTENT(in) :: topolbodyid(ntopol)
        INTEGER(intk), INTENT(in) :: ntrimax
        INTEGER(intk), INTENT(in) :: triau(ntrimax, kk, jj, ii), &
            triav(ntrimax, kk, jj, ii), triaw(ntrimax, kk, jj, ii)
        REAL(realk), INTENT(in) :: knoten(kk, jj, ii)
        REAL(realk), INTENT(inout) :: kanteu(kk, jj, ii), kantev(kk, jj, ii), &
            kantew(kk, jj, ii)
        INTEGER(intk), INTENT(in) :: bzelltyp(kk, jj, ii)
        REAL(realk), INTENT(out) :: au(kk, jj, ii), av(kk, jj, ii), &
            aw(kk, jj, ii)
        REAL(realk), INTENT(out), OPTIONAL :: xpsw(:, :)
        REAL(realk), INTENT(out), OPTIONAL :: yus(kk, jj, ii), &
            zus(kk, jj, ii), xvs(kk, jj, ii), zvs(kk, jj, ii), &
            xws(kk, jj, ii), yws(kk, jj, ii)

        ! Local variables
        ! TODO: Warning: this parameter is not fully honored. Please
        ! do not change from 20 without carefully reviewing all code.
        ! There are lots of places in this file where 20 is still hardcoded.
        INTEGER(intk), PARAMETER :: npmax = 20
        INTEGER(intk) :: knotenzelle(8), kanten(12)
        INTEGER(intk) :: nka(6), connectka(2, 6), nkn(6), connectkn(4, 6), &
            nwa, connectwa(7)
        INTEGER(intk) :: bodyidkanten(12)
        INTEGER(intk) :: connect(6), ncon, fltyp
        INTEGER(intk) :: ika, ikn, npoints, i, j, k, icell, idx
        REAL(realk) :: xx(3, npmax)
        REAL(realk) :: area, cartarea, s1, s2
        LOGICAL :: has_opt

        ! AU,AV,AW: Anteil der offenen Flaeche
        ! yus, zus, xvs, zvs, xws, yws: Flaechneschwerpunkte der offenen Flaechen

        au = 0.0
        av = 0.0
        aw = 0.0

        connect = 0
        ncon = 0

        IF (PRESENT(yus) .AND. PRESENT(zus) .AND. PRESENT(xvs) .AND. &
                PRESENT(zvs) .AND. PRESENT(xws) .AND. PRESENT(yws)) THEN
            ! All optional arguments are present - fine!
            has_opt = .TRUE.
        ELSE IF (PRESENT(yus) .OR. PRESENT(zus) .OR. PRESENT(xvs) .OR. &
                PRESENT(zvs) .OR. PRESENT(xws) .OR. PRESENT(yws)) THEN
            ! One or more optional arguments are present but not all - wrong!
            CALL errr(__FILE__, __LINE__)
        ELSE
            ! No optional arguments
            has_opt = .FALSE.
        END IF

        CALL blockcheck_grid(kk, jj, ii, kanteu, kantev, kantew, &
                knoten, 1)

        ! Indices must be the same as in calc_nvecs_grid in gc_mod.F90!!!
        DO i = 2, ii
            DO j = 2, jj
                DO k = 2, kk
                    npoints = 0

                    IF (bzelltyp(k, j, i) == 1) THEN
                        au(k, j, i) = 1.0
                        av(k, j, i) = 1.0
                        aw(k, j, i) = 1.0
                        IF (i == 2) THEN
                            ! FRO Flaeche
                            au(k, j, i-1) = 1.0
                        END IF
                        IF (j == 2) THEN
                            ! RGT Flaeche
                            av(k, j-1, i) = 1.0
                        END IF
                        IF (k == 2) THEN
                            ! BOT Flaeche
                            aw(k-1, j, i) = 1.0
                        END IF

                        IF (has_opt) THEN
                            yus(k, j, i) = ystag(j) - 0.5*ddy(j)
                            zus(k, j, i) = zstag(k) - 0.5*ddz(k)
                            xvs(k, j, i) = xstag(i) - 0.5*ddx(i)
                            zvs(k, j, i) = zstag(k) - 0.5*ddz(k)
                            xws(k, j, i) = xstag(i) - 0.5*ddx(i)
                            yws(k, j, i) = ystag(j) - 0.5*ddy(j)

                            IF (i == 2) THEN
                                ! FRO Flaeche
                                yus(k, j, i-1) = ystag(j) - 0.5*ddy(j)
                                zus(k, j, i-1) = zstag(k) - 0.5*ddz(k)
                            END IF
                            IF (j == 2) THEN
                                ! RGT Flaeche
                                xvs(k, j-1, i) = xstag(i) - 0.5*ddx(i)
                                zvs(k, j-1, i) = zstag(k) - 0.5*ddz(k)
                            END IF
                            IF (k == 2) THEN
                                ! BOT Flaeche
                                xws(k-1, j, i) = xstag(i) - 0.5*ddx(i)
                                yws(k-1, j, i) = ystag(j) - 0.5*ddy(j)
                            END IF
                        END IF

                    ELSE IF (bzelltyp(k, j, i) == 0) THEN
                        au(k, j, i) = 0.0
                        av(k, j, i) = 0.0
                        aw(k, j, i) = 0.0
                        IF (i == 2) THEN
                            ! FRO Flaeche
                            AU(k, j, i-1) = 0.0
                        END IF
                        IF (j == 2) THEN
                            ! RGT Flaeche
                            av(k, j-1, i) = 0.0
                        END IF
                        IF (k == 2) THEN
                            ! BOT Flaeche
                            aw(k-1, j, i) = 0.0
                        END IF

                        IF (has_opt) THEN
                            yus(k, j, i) = ystag(j) - 0.5*ddy(j)
                            zus(k, j, i) = zstag(k) - 0.5*ddz(k)
                            xvs(k, j, i) = xstag(i) - 0.5*ddx(i)
                            zvs(k, j, i) = zstag(k) - 0.5*ddz(k)
                            xws(k, j, i) = xstag(i) - 0.5*ddx(i)
                            yws(k, j, i) = ystag(j) - 0.5*ddy(j)

                            IF (i == 2) THEN
                                ! FRO Flaeche
                                yus(k, j, i-1) = ystag(j) - 0.5*ddy(j)
                                zus(k, j, i-1) = zstag(k) - 0.5*ddz(k)
                            END IF
                            IF (j == 2) THEN
                                ! RGT Flaeche
                                xvs(k, j-1, i) = xstag(i) - 0.5*ddx(i)
                                zvs(k, j-1, i) = zstag(k) - 0.5*ddz(k)
                            END IF
                            IF (k == 2) THEN
                                ! BOT Flaeche
                                xws(k-1, j, i) = xstag(i) - 0.5*ddx(i)
                                yws(k-1, j, i) = ystag(j) - 0.5*ddy(j)
                            END IF
                        END IF

                    ELSE IF (bzelltyp(k, j, i) < 0) THEN
                        icell = -bzelltyp(k, j, i)

                        CALL punktekoordinaten(k, j, i, kk, jj, ii, &
                            xstag, ystag, zstag, &
                            ntopol, topol, topolbodyid, &
                            ntrimax, triau, triav, triaw, &
                            knoten, kanteu, kantev, kantew, &
                            npoints, npmax, xx, bodyidkanten)

                        CALL knotenundkanten(k, j, i, kk, jj, ii, 1, &
                            knoten, kanteu, kantev, kantew, &
                            knotenzelle, kanten)

                        CALL flzelle(k, j, i, kk, jj, ii, 1, &
                            knotenzelle, kanten, nka, connectka, &
                            nkn, connectkn, nwa, connectwa)

                        ! ueber alle Quaderseiten
                        DO idx = 1, 6
                            IF (nka(idx) + nkn(idx) > 0) THEN
                                ncon = nka(idx) + nkn(idx)
                                fltyp = 1
                                DO ika = 1, nka(idx)
                                    fltyp = 0
                                    connect(ika) = connectka(ika, idx) + 8 &
                                        + npoints - 20
                                END DO
                                DO ikn = 1, nkn(idx)
                                    connect(ikn + nka(idx)) = &
                                        connectkn(ikn, idx) + npoints - 20
                                END DO

                                IF (i == 2 .AND. idx == 1) THEN
                                    ! FRO Flaeche
                                    IF ((knotenzelle(1) + knotenzelle(5) &
                                            + knotenzelle(8) &
                                            + knotenzelle(4)) == 4) THEN
                                        au(k, j, i-1) = 1.0
                                        s1 = ystag(j) - 0.5*ddy(j)
                                        s2 = zstag(k) - 0.5*ddz(k)
                                    ELSE
                                        CALL calcfacedata(k, j, i, 1, &
                                            ncon, connect, npoints, xx, &
                                            au(k, j, i-1), s1, s2)
                                        ! MIN to avoid open face fractions
                                        ! bigger than 1 caused by roundoff
                                        ! errors
                                        au(k, j, i-1) = MIN(au(k, j, i-1) &
                                            /(ddy(j)*ddz(k)), 1.0_realk)
                                    END IF
                                    IF (has_opt) THEN
                                        yus(k, j, i-1) = s1
                                        zus(k, j, i-1) = s2
                                    END IF
                                END IF

                                IF (j == 2 .AND. idx == 3) THEN
                                    ! RGT Flaeche
                                    IF ((knotenzelle(1) + knotenzelle(2) &
                                            + knotenzelle(3) &
                                            + knotenzelle(4)) == 4) THEN
                                        av(k, j-1, i) = 1.0
                                        s1 = zstag(k) - 0.5*ddz(k)
                                        s2 = xstag(i) - 0.5*ddx(i)
                                    ELSE
                                        CALL calcfacedata(k, j, i, 2, &
                                            ncon, connect, npoints, xx, &
                                            av(k, j-1, i), s1, s2)
                                        av(k, j-1, i) = MIN(av(k, j-1, i) &
                                            /(ddx(i)*ddz(k)), 1.0_realk)
                                    END IF
                                    IF (has_opt) THEN
                                        zvs(k, j-1, i) = s1
                                        xvs(k, j-1, i) = s2
                                    END IF
                                END IF
                                IF (k == 2 .AND. idx == 5) THEN
                                    ! BOT Flaeche
                                    IF ((knotenzelle(1) + knotenzelle(2) &
                                            + knotenzelle(6) &
                                            + knotenzelle(5)) == 4) THEN
                                        aw(k-1, j, i) = 1.0
                                        s1 = xstag(i) - 0.5*ddx(i)
                                        s2 = ystag(j) - 0.5*ddy(j)
                                    ELSE
                                        CALL calcfacedata(k, j, i, 3, &
                                            ncon, connect, npoints, xx, &
                                            aw(k-1, j, i), s1, s2)
                                        aw(k-1, j, i) = MIN(aw(k-1, j, i) &
                                            /(ddy(j)*ddx(i)), 1.0_realk)
                                    END IF
                                    IF (has_opt) THEN
                                        xws(k-1, j, i) = s1
                                        yws(k-1, j, i) = s2
                                    END IF
                                END IF
                                IF (idx == 2) THEN
                                    ! BAC Flaeche
                                    IF ((knotenzelle(2) + knotenzelle(6) + &
                                            knotenzelle(7) + &
                                            knotenzelle(3)) == 4) THEN
                                        au(k, j, i) = 1.0
                                        s1 = ystag(j) - 0.5*ddy(j)
                                        s2 = zstag(k) - 0.5*ddz(k)
                                    ELSE
                                        CALL calcfacedata(k, j, i, 1, &
                                            ncon, connect, npoints, xx, &
                                            au(k, j, i), s1, s2)
                                        au(k, j, i) = MIN(au(k, j, i) &
                                            /(ddy(j)*ddz(k)), 1.0_realk)
                                    END IF
                                    IF (has_opt) THEN
                                        yus(k, j, i) = s1
                                        zus(k, j, i) = s2
                                    END IF
                                END IF
                                IF (idx == 4) THEN
                                    ! LFT Flaeche
                                    IF ((knotenzelle(5) + knotenzelle(6) &
                                            + knotenzelle(7) &
                                            + knotenzelle(8)) == 4) THEN
                                        av(k, j, i) = 1.
                                        s1 = zstag(k) - 0.5*ddz(k)
                                        s2 = xstag(i) - 0.5*ddx(i)
                                    ELSE
                                        CALL calcfacedata(k, j, i, 2, &
                                            ncon, connect, npoints, xx, &
                                            av(k, j, i), s1, s2)
                                        av(k, j, i) = MIN(av(k, j, i) &
                                            /(ddx(i)*ddz(k)), 1.0_realk)
                                    END IF
                                    IF (has_opt) THEN
                                        zvs(k, j, i) = s1
                                        xvs(k, j, i) = s2
                                    END IF
                                END IF
                                IF (idx == 6) THEN
                                    ! TOP Flaeche
                                    IF ((knotenzelle(3) + knotenzelle(4) &
                                            + knotenzelle(8) &
                                            + knotenzelle(7)) == 4) THEN
                                        aw(k, j, i) = 1.
                                        s1 = xstag(i) - 0.5*ddx(i)
                                        s2 = ystag(j) - 0.5*ddy(j)
                                    ELSE
                                        CALL calcfacedata(k, j, i, 3, &
                                            ncon, connect, npoints, xx, &
                                            aw(k, j, i), s1, s2)
                                        aw(k, j, i) = MIN(aw(k, j, i) &
                                            /(ddy(j)*ddx(i)), 1.0_realk)
                                    END IF
                                    IF (has_opt) THEN
                                        xws(k, j, i) = s1
                                        yws(k, j, i) = s2
                                    END IF
                                END IF
                            ELSE
                                IF (has_opt) THEN
                                    IF (idx == 1 .AND. i == 2) THEN
                                        yus(k, j, i-1) = ystag(j) - 0.5*ddy(j)
                                        zus(k, j, i-1) = zstag(k) - 0.5*ddz(k)
                                    ELSE IF (idx == 3 .AND. j == 2) THEN
                                        xvs(k, j-1, i) = xstag(i) - 0.5*ddx(i)
                                        zvs(k, j-1, i) = zstag(k) - 0.5*ddz(k)
                                    ELSE IF (idx == 5 .AND. k == 2) THEN
                                        xws(k-1, j, i) = xstag(i) - 0.5*ddx(i)
                                        yws(k-1, j, i) = ystag(j) - 0.5*ddy(j)
                                    ELSE IF (idx == 2) THEN
                                        yus(k, j, i) = ystag(j) - 0.5*ddy(j)
                                        zus(k, j, i) = zstag(k) - 0.5*ddz(k)
                                    ELSE IF (idx == 4) THEN
                                        xvs(k, j, i) = xstag(i) - 0.5*ddx(i)
                                        zvs(k, j, i) = zstag(k) - 0.5*ddz(k)
                                    ELSE IF (idx == 6) THEN
                                        xws(k, j, i) = xstag(i) - 0.5*ddx(i)
                                        yws(k, j, i) = ystag(j) - 0.5*ddy(j)
                                    END IF
                                END IF
                            END IF
                        END DO

                        IF (PRESENT(xpsw)) THEN
                            xpsw = 0.0

                            area = SQRT( &
                                (ddz(k)*ddy(j)*(au(k, j, i) &
                                - au(k, j, i-1)))**2 &
                                + (ddx(i)*ddz(k)*(av(k, j, i)&
                                - av(k, j-1, i)))**2 &
                                + (ddx(i)*ddy(j)*(aw(k, j, i)&
                                - aw(k-1, j, i)))**2)

                            ! Referenzflaeche: Mittel der kartesischen
                            ! Flaechen aus allen drei Richtungen.
                            cartarea = (1.0/3.0)*(ddx(i)*ddy(j) + &
                                ddx(i)*ddz(k) + ddz(k)*ddy(j))

                            ! Kann vorkommen, dass Flaeche zu klein,
                            ! dann muss Notfallnormale berechnet werden.
                            IF (area/cartarea < maccur**2) THEN
                                CALL calcwallfacecenterrescue(nwa, npoints, &
                                    connectwa, xx, xpsw(1, icell), &
                                    xpsw(2, icell), xpsw(3, icell))
                            ELSE
                                CALL calcwallfacecenter(k, j, i, nwa, npoints, &
                                    connectwa, xx, cartarea, area, &
                                    xpsw(1, icell), xpsw(2, icell), &
                                    xpsw(3, icell))
                            END IF
                        END IF
                    ELSE
                        CALL errr(__FILE__, __LINE__)
                    END IF
                END DO
            END DO
        END DO
    END SUBROUTINE calcauavaw_grid

END MODULE calcauavaw_mod
