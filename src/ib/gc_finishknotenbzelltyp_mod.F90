MODULE gc_finishknotenbezelltyp_mod
    USE core_mod, ONLY: realk, intk, errr, get_fieldptr, ngrid, minlevel, &
        maxlevel, nmygridslvl, mygridslvl, idim3d, get_mgdims, get_ip3, &
        get_ip3n, field_t
    USE checkzelle_mod, ONLY: checkzelle_grid
    USE gc_zelltyp_mod, ONLY: zelltyp_grid
    USE blockcheck_mod, ONLY: blockcheck_grid
    USE findinterface_mod, ONLY: findinterface2
    USE topol_mod, ONLY: topol_t

    IMPLICIT NONE(type, external)
    PRIVATE

    PUBLIC :: finishknotenbezelltyp

CONTAINS
    SUBROUTINE finishknotenbezelltyp(topol, ntrimax, triau, triav, triaw, &
            bzelltyp, bp, knoten, kanteu, kantev, kantew, icells)

        ! Subroutine arguments
        TYPE(topol_t), INTENT(in) :: topol
        INTEGER(intk), INTENT(in) :: ntrimax
        INTEGER(intk), INTENT(in) :: triau(*)
        INTEGER(intk), INTENT(in) :: triav(*)
        INTEGER(intk), INTENT(in) :: triaw(*)
        TYPE(field_t), INTENT(in) :: bp
        INTEGER(intk), INTENT(inout) :: bzelltyp(*)
        TYPE(field_t), INTENT(inout) :: knoten
        TYPE(field_t), INTENT(inout) :: kanteu
        TYPE(field_t), INTENT(inout) :: kantev
        TYPE(field_t), INTENT(inout) :: kantew
        INTEGER(intk), INTENT(out) :: icells(:)

        ! Local variables
        INTEGER(intk) :: ilevel

        ! Sanity check
        IF (SIZE(icells) /= ngrid) CALL errr(__FILE__, __LINE__)
        icells = 0

        DO ilevel = minlevel, maxlevel
            CALL finishknotenbezelltyp_level(ilevel, topol, ntrimax, &
                triau, triav, triaw, bp, bzelltyp, knoten, &
                kanteu, kantev, kantew, icells)
        END DO
    END SUBROUTINE finishknotenbezelltyp


    SUBROUTINE finishknotenbezelltyp_level(ilevel, topol, ntrimax, triau, &
            triav, triaw, bp, bzelltyp, knoten, kanteu, kantev, kantew, icells)
        ! Subroutine arguments
        INTEGER(intk), INTENT(in) :: ilevel
        TYPE(topol_t), INTENT(in) :: topol
        INTEGER(intk), INTENT(in) :: ntrimax
        INTEGER(intk), INTENT(in) :: triau(*)
        INTEGER(intk), INTENT(in) :: triav(*)
        INTEGER(intk), INTENT(in) :: triaw(*)
        TYPE(field_t), INTENT(in) :: bp
        INTEGER(intk), INTENT(inout) :: bzelltyp(*)
        TYPE(field_t), INTENT(inout) :: knoten
        TYPE(field_t), INTENT(inout) :: kanteu
        TYPE(field_t), INTENT(inout) :: kantev
        TYPE(field_t), INTENT(inout) :: kantew
        INTEGER(intk), INTENT(out) :: icells(:)

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

            CALL finishknotenbezelltyp_grid(kk, jj, ii, xstag, ystag, zstag, &
                ddx, ddy, ddz, topol%n, topol%topol, ntrimax, triau(ip3n), &
                triav(ip3n), triaw(ip3n), bp%arr(ip3), bzelltyp(ip3), &
                knoten%arr(ip3), kanteu%arr(ip3), kantev%arr(ip3), &
                kantew%arr(ip3), icells(igrid))
        END DO
    END SUBROUTINE finishknotenbezelltyp_level


    SUBROUTINE finishknotenbezelltyp_grid(kk, jj, ii, xstag, ystag, zstag, &
            ddx, ddy, ddz, ntopol, topol, ntrimax, triau, triav, triaw, bp, &
            bzelltyp, knoten, kanteu, kantev, kantew, icells)

        ! Subroutine arguments
        INTEGER(intk), INTENT(in) :: kk, jj, ii
        REAL(realk), INTENT(in) :: xstag(ii), ystag(jj), zstag(kk)
        REAL(realk), INTENT(in) :: ddx(ii), ddy(jj), ddz(kk)
        INTEGER(intk), INTENT(in) :: ntopol
        REAL(realk), INTENT(in) :: topol(3, 3, ntopol)
        INTEGER(intk), INTENT(in) :: ntrimax
        INTEGER(intk), INTENT(in) :: triau(ntrimax, kk, jj, ii), &
            triav(ntrimax, kk, jj, ii), triaw(ntrimax, kk, jj, ii)
        REAL(realk), INTENT(in) :: bp(kk, jj, ii)
        INTEGER(intk), INTENT(inout) :: bzelltyp(kk, jj, ii)
        REAL(realk), INTENT(inout) :: knoten(kk, jj, ii)
        REAL(realk), INTENT(inout) :: kanteu(kk, jj, ii), kantev(kk, jj, ii), &
            kantew(kk, jj, ii)
        INTEGER(intk), INTENT(out) :: icells

        ! Local variables
        INTEGER(intk) :: k, j, i
        INTEGER(intk) :: suminfo, flag, iloop
        INTEGER(intk) :: found, foundx1, foundx2, foundy1, foundy2, &
            foundz1, foundz2, foundnr

        ! Knoten, die nur geblockte Druckzellen haben,
        ! werden geblockt
        DO i = 2, ii-1
            DO j = 2, jj-1
                DO k = 2, kk-1
                    knoten(k, j, i) = MIN(knoten(k,j,i), &
                        bp(k+1, j+1, i+1) &
                        + bp(k+1, j+1, i) &
                        + bp(k+1, j, i+1) &
                        + bp(k+1, j, i) &
                        + bp(k, j+1, i+1) &
                        + bp(k, j+1, i) &
                        + bp(k, j, i+1) &
                        + bp(k, j, i))
                END DO
            END DO
        END DO

        ! Knoten, der Druckzellen die Durchlaesse verstopfen
        DO i = 2, ii-1
            DO j = 2, jj-1
                DO k = 2, kk-1
                    CALL findinterface2(k, j, i, kk, jj, ii, bp, found, &
                        foundx1, foundx2, foundy1, foundy2, foundz1, &
                        foundz2, foundnr)
                    IF (bp(k, j, i) < 0.5) THEN
                        IF (foundx1*foundx2 + foundy1*foundy2 &
                                + foundz1*foundz2 == 1 .AND. foundnr == 2) THEN
                            knoten(k, j, i) = 0.0
                            knoten(k, j, i-1) = 0.0
                            knoten(k-1, j, i-1) = 0.0
                            knoten(k-1, j, i) = 0.0
                            knoten(k, j-1, i) = 0.0
                            knoten(k, j-1, i-1) = 0.0
                            knoten(k-1, j-1, i-1) = 0.0
                            knoten(k-1, j-1, i) = 0.0
                        END IF
                    END IF
                END DO
            END DO
        END DO

        ! Konsistentes Knoten und Kantenfeld, nur Kanten werden modifiziert
        flag = 1
        CALL blockcheck_grid(kk, jj, ii, kanteu, kantev, kantew, knoten, flag)

        ! Behandlung von Zellen mit nichtzusammenhaengenden
        ! offenen / geschlossenen Bereichen
        DO iloop = 1, 20
            CALL checkzelle_grid(kk, jj, ii, kanteu, kantev, kantew, knoten, &
                ntrimax, triau, triav, triaw, xstag, ystag, zstag, &
                ddx, ddy, ddz, ntopol, topol, suminfo)
            IF (suminfo == 0) EXIT
        END DO
        IF (suminfo /= 0) THEN
            WRITE(*, *) 'checkzelle suminfo', suminfo
            CALL errr(__FILE__, __LINE__)
        END IF

        CALL zelltyp_grid(kk, jj, ii, knoten, bzelltyp, icells)
    end SUBROUTINE finishknotenbezelltyp_grid

END MODULE gc_finishknotenbezelltyp_mod
