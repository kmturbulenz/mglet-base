MODULE calcnormals_mod
    USE core_mod, ONLY: realk, intk, errr, get_fieldptr, ngrid, minlevel, &
        maxlevel, nmygridslvl, mygridslvl, idim3d, get_mgdims, get_ip3, &
        get_ip3n, sub2ind, most_frequent_nonzero, field_t, &
        connect
    USE blockcheck_mod, ONLY: blockcheck_grid
    USE punktekoordinaten_mod, ONLY: punktekoordinaten
    USE knotenundkanten_mod, ONLY: knotenundkanten
    USE flzelle_mod, ONLY: flzelle
    USE topol_mod, ONLY: topol_t
    USE stencils_mod, ONLY: stencils_t

    IMPLICIT NONE(type, external)
    PRIVATE

    LOGICAL, PARAMETER :: write_geom_with_gl = .TRUE.
    LOGICAL, PARAMETER :: write_geom_onlyactivecells = .FALSE.

    PUBLIC :: calcnormals

CONTAINS
    SUBROUTINE calcnormals(velocity, topol, ntrimax, triau, triav, triaw, &
            knoten, kanteu, kantev, kantew, bzelltyp, au, av, aw, &
            icells, icellspointer, ncellstot, bodyid, ucell, &
            stencils, writegeom)

        ! Subroutine arguments
        REAL(realk), INTENT(in) :: velocity(:, :)
        TYPE(topol_t), INTENT(in) :: topol
        INTEGER(intk), INTENT(in) :: ntrimax
        INTEGER(intk), INTENT(in) :: triau(*), triav(*), triaw(*)
        TYPE(field_t), INTENT(in) :: knoten
        TYPE(field_t), INTENT(inout) :: kanteu, kantev, kantew
        INTEGER(intk), INTENT(in) :: bzelltyp(idim3d)
        TYPE(field_t), INTENT(in) :: au, av, aw
        INTEGER(intk), INTENT(in) :: icells(:)
        INTEGER(intk), INTENT(in) :: icellspointer(:)
        INTEGER(intk), INTENT(in) :: ncellstot
        INTEGER(intk), INTENT(out) :: bodyid(:)
        REAL(realk), INTENT(out) :: ucell(:, :)
        CLASS(stencils_t), INTENT(inout) :: stencils
        LOGICAL, INTENT(in), OPTIONAL :: writegeom

        ! Local variables
        INTEGER(intk) :: ilevel

        DO ilevel = minlevel, maxlevel
            CALL calcnormals_level(ilevel, velocity, topol, ntrimax, &
                triau, triav, triaw, knoten, kanteu, kantev, kantew, &
                bzelltyp, au, av, aw, icells, icellspointer, ncellstot, &
                bodyid, ucell, stencils, writegeom)
        END DO
    END SUBROUTINE calcnormals


    SUBROUTINE calcnormals_level(ilevel, velocity, topol, ntrimax, &
            triau, triav, triaw, knoten, kanteu, kantev, kantew, &
            bzelltyp, au, av, aw, icells, icellspointer, ncellstot, &
            bodyid, ucell, stencils, writegeom)

        ! Subroutine arguments
        INTEGER(intk), INTENT(in) :: ilevel
        REAL(realk), INTENT(in) :: velocity(:, :)
        TYPE(topol_t), INTENT(in) :: topol
        INTEGER(intk), INTENT(in) :: ntrimax
        INTEGER(intk), INTENT(in) :: triau(*), triav(*), triaw(*)
        TYPE(field_t), INTENT(in) :: knoten
        TYPE(field_t), INTENT(inout) :: kanteu, kantev, kantew
        INTEGER(intk), INTENT(in) :: bzelltyp(idim3d)
        TYPE(field_t), INTENT(in) :: au, av, aw
        INTEGER(intk), INTENT(in) :: icells(:)
        INTEGER(intk), INTENT(in) :: icellspointer(:)
        INTEGER(intk), INTENT(in) :: ncellstot
        INTEGER(intk), INTENT(out) :: bodyid(:)
        REAL(realk), INTENT(out) :: ucell(:, :)
        CLASS(stencils_t), INTENT(inout) :: stencils
        LOGICAL, INTENT(in), OPTIONAL :: writegeom

        ! Local variables
        INTEGER(intk) :: i, igrid, kk, jj, ii, ip3, ip3n, ipp, ncells
        REAL(realk), POINTER, CONTIGUOUS :: xstag(:), ystag(:), zstag(:), &
            ddx(:), ddy(:), ddz(:), finecell(:, :, :)


        DO i = 1, nmygridslvl(ilevel)
            igrid = mygridslvl(i, ilevel)

            CALL get_fieldptr(xstag, "XSTAG", igrid)
            CALL get_fieldptr(ystag, "YSTAG", igrid)
            CALL get_fieldptr(zstag, "ZSTAG", igrid)

            CALL get_fieldptr(ddx, "DDX", igrid)
            CALL get_fieldptr(ddy, "DDY", igrid)
            CALL get_fieldptr(ddz, "DDZ", igrid)

            CALL get_fieldptr(finecell, "FINECELL", igrid)

            CALL get_mgdims(kk, jj, ii, igrid)
            CALL get_ip3(ip3, igrid)
            CALL get_ip3n(ip3n, ntrimax, igrid)

            ipp = icellspointer(igrid)
            ncells = icells(igrid)

            ! Ugly code, because I'm not allowed to do yus(ip3) if yus
            ! is not present, therefore I must do this
            CALL calcnormals_grid(igrid, kk, jj, ii, xstag, ystag, zstag, &
                ddx, ddy, ddz, velocity, topol%n, topol%topol, topol%bodyid, &
                ntrimax, triau(ip3n), triav(ip3n), triaw(ip3n), finecell, &
                knoten%arr(ip3), kanteu%arr(ip3), kantev%arr(ip3), &
                kantew%arr(ip3), bzelltyp(ip3), au%arr(ip3), av%arr(ip3), &
                aw%arr(ip3), icells(igrid), bodyid(ipp:ipp+ncells-1), &
                ucell(:, ipp:ipp+ncells-1), stencils, writegeom)
        END DO
    END SUBROUTINE calcnormals_level


    SUBROUTINE calcnormals_grid(igrid, kk, jj, ii, xstag, ystag, zstag, &
            ddx, ddy, ddz, velocity, ntopol, topol, topolbodyid, &
            ntrimax, triau, triav, triaw, finecell, &
            knoten, kanteu, kantev, kantew, bzelltyp, au, av, aw, &
            ncells, bodyid, ucell, stencils, writegeom)

        ! Subroutine arguments
        INTEGER(intk), INTENT(in) :: igrid
        INTEGER(intk), INTENT(in) :: kk, jj, ii
        REAL(realk), INTENT(in) :: xstag(ii), ystag(jj), zstag(kk)
        REAL(realk), INTENT(in) :: ddx(ii), ddy(jj), ddz(kk)
        REAL(realk), INTENT(in) :: velocity(:, :)
        INTEGER(intk), INTENT(in) :: ntopol
        REAL(realk), INTENT(in) :: topol(3, 3, ntopol)
        INTEGER(intk), INTENT(in) :: topolbodyid(ntopol)
        INTEGER(intk), INTENT(in) :: ntrimax
        INTEGER(intk), INTENT(in) :: triau(ntrimax, kk, jj, ii), &
            triav(ntrimax, kk, jj, ii), triaw(ntrimax, kk, jj, ii)
        REAL(realk), INTENT(in) :: finecell(kk, jj, ii)
        REAL(realk), INTENT(in) :: knoten(kk, jj, ii)
        REAL(realk), INTENT(inout) :: kanteu(kk, jj, ii), kantev(kk, jj, ii), &
            kantew(kk, jj, ii)
        INTEGER(intk), INTENT(in) :: bzelltyp(kk, jj, ii)
        REAL(realk), INTENT(in) :: au(kk, jj, ii), av(kk, jj, ii), &
            aw(kk, jj, ii)
        INTEGER(intk), INTENT(in) :: ncells
        INTEGER(intk), INTENT(out) :: bodyid(:)
        REAL(realk), INTENT(out) :: ucell(:, :)
        CLASS(stencils_t), INTENT(inout) :: stencils
        LOGICAL, INTENT(in), OPTIONAL :: writegeom

        ! Local variables
        INTEGER(intk) :: i, j, k
        INTEGER(intk) :: nka(6), connectka(2, 6), nkn(6), connectkn(4, 6), &
            nwa, connectwa(7)

        INTEGER(intk) :: icell, icellsgeom
        INTEGER(intk) :: ika, npoints
        INTEGER(intk) :: knotenzelle(8), kanten(12)
        INTEGER(intk) :: bodyidkanten(12)
        INTEGER(intk) :: nfixed_cells

        REAL(realk) :: umean(3), sum, this_vel(3)
        INTEGER(intk) :: ikante, idir, idmin, irepeat

        INTEGER(intk), ALLOCATABLE :: icelllist(:), cellind(:)
        INTEGER(intk), ALLOCATABLE :: connect(:, :), ncon(:)
        REAL(realk), ALLOCATABLE :: xx(:, :), arealist(:)

        INTEGER(intk) :: indi, indj, indk
        REAL(realk) :: indnn

        ALLOCATE(icelllist(ncells))
        ALLOCATE(arealist(ncells))
        ALLOCATE(connect(6, ncells*7))
        ALLOCATE(ncon(ncells*7))
        ALLOCATE(xx(3, ncells*20))
        ALLOCATE(cellind(ncells))

        icelllist = 0
        arealist = 0.0
        connect = 0
        ncon = 0
        xx = 0.0
        cellind = 0
        npoints = 0

        ! Korrektur von kanteu, kantev, kantew
        CALL blockcheck_grid(kk, jj, ii, kanteu, kantev, kantew, knoten, 1)

        ! Loop indices in this loop must be identical to indices in
        ! stencils_mod.F90
        icellsgeom = 0
        DO i = 2, ii
            DO j = 2, jj
                DO k = 2, kk
                    IF (bzelltyp(k, j, i) >= 0) CYCLE
                    icell = -bzelltyp(k, j, i)

                    ! Berechnung xx, ukanten, bodyidkanten
                    CALL punktekoordinaten(k, j, i, kk, jj, ii, &
                        xstag, ystag, zstag, &
                        ntopol, topol, topolbodyid, &
                        ntrimax, triau, triav, triaw, &
                        knoten, kanteu, kantev, kantew, &
                        npoints, 20*ncells, xx, bodyidkanten)

                    ! Belegen von knotenzelle, kanten
                    ! Konsistenzcheck ob Anzahl der Punkte > 3
                    CALL knotenundkanten(k, j, i, kk, jj, ii, 0, knoten, &
                        kanteu, kantev, kantew, knotenzelle, kanten)

                    ! belegen von nka,connectka,nkn,connectkn,nwa,connectwa
                    CALL flzelle(k, j, i, kk, jj, ii, 0, &
                        knotenzelle, kanten, nka, connectka, nkn, connectkn, &
                        nwa, connectwa)

                    ! Connect for real geometry plane parts
                    ncon(icell) = nwa
                    DO ika = 1, nwa
                        connect(ika, icell) = connectwa(ika)+8-1 + npoints-20
                    END DO

                    arealist(icell) = SQRT( &
                        (ddz(k)*ddy(j)*(au(k, j, i) - au(k, j, i-1)))**2 &
                        + (ddx(i)*ddz(k)*(av(k, j, i) - av(k, j-1, i)))**2 &
                        + (ddx(i)*ddy(j)*(aw(k, j, i) - aw(k-1, j, i)))**2)

                    umean = 0.0
                    idmin = 0

                    sum = 0.0
                    DO ikante = 1, 12
                        IF (kanten(ikante) == 1) THEN

                            IF (bodyidkanten(ikante) > 0) THEN
                                this_vel = velocity(:, bodyidkanten(ikante))
                            ELSE
                                this_vel = 0.0
                            END IF

                            DO idir = 1, 3
                                umean(idir) = umean(idir) + this_vel(idir)
                            END DO

                            ! Mittelung nur ueber Kanten mit Geschwdinigkeit
                            ! /= Null
                            IF (SQRT(this_vel(1)**2 + this_vel(2)**2 &
                                    + this_vel(3)**2) > 0.0) THEN
                                sum = sum + 1.0
                            END IF
                        END IF
                        IF (bodyidkanten(ikante) > 0) THEN
                            idmin = MAX(idmin, bodyidkanten(ikante))
                        END IF
                    END DO

                    ! TODO: implement better alternative
                    ! see util_mod.F90 most_frequent_nonzero
                    ! bodyid(icell) = most_frequent_nonzero(bodyidkanten)
                    bodyid(icell) = idmin

                    IF (sum < 0.5) sum = 1.0
                    DO idir = 1, 3
                        ucell(idir, icell) = umean(idir)/sum
                    END DO

                    IF (write_geom_with_gl) THEN
                        indnn = 1.0
                    ELSE
                        indi = MAX(0, MIN(1, i-3+1))*MAX(0, MIN(1, ii-2-i+1))
                        indj = MAX(0, MIN(1, j-3+1))*MAX(0, MIN(1, jj-2-j+1))
                        indk = MAX(0, MIN(1, k-3+1))*MAX(0, MIN(1, kk-2-k+1))
                        indnn = REAL(indi*indj*indk)
                    END IF
                    IF (write_geom_onlyactivecells) THEN
                        indnn = indnn*finecell(k, j, i)
                    END IF

                    IF (indnn > 0.5) THEN
                        icellsgeom = icellsgeom + 1
                        icelllist(icellsgeom) = icell
                    END IF

                    CALL sub2ind(cellind(icell), k, j, i, kk, jj, ii)
                END DO
            END DO
        END DO

        ! Try to fix missing bodyid's
        DO irepeat = 1, 3
            CALL fixup_bodyid(kk, jj, ii, bzelltyp, bodyid, nfixed_cells)
            IF (nfixed_cells == 0) EXIT
        END DO

        IF (PRESENT(writegeom)) THEN
            IF (writegeom) THEN
                CALL stencils%set_geometry(icellsgeom, ncells, icelllist, &
                    ncon, connect, xx, cellind, arealist, igrid)
                ! CALL stencils%set_bodyid(igrid, ncells, bodyid, &
                !     icellsgeom, icelllist)
            END IF
        END IF

        DEALLOCATE(cellind)
        DEALLOCATE(xx)
        DEALLOCATE(ncon)
        DEALLOCATE(connect)
        DEALLOCATE(arealist)
        DEALLOCATE(icelllist)
    END SUBROUTINE calcnormals_grid


    SUBROUTINE fixup_bodyid(kk, jj, ii, bzelltyp, bodyid, nfixed_cells)
        ! Subroutine arguments
        INTEGER(intk), INTENT(in) :: kk, jj, ii
        INTEGER(intk), INTENT(in) :: bzelltyp(kk, jj, ii)
        INTEGER(intk), INTENT(inout) :: bodyid(:)
        INTEGER(intk), INTENT(out) :: nfixed_cells

        ! Local variables
        INTEGER(intk) :: i, j, k, icell, ip, jp, kp, icell_lauf
        ! INTEGER(intk) :: nval, values(27)

        nfixed_cells = 0
        DO i = 2, ii
            DO j = 2, jj
                DO k = 2, kk
                    IF (bzelltyp(k, j, i) >= 0) CYCLE
                    icell = -bzelltyp(k, j, i)

                    IF (bodyid(icell) == 0) THEN
                        ! TODO: implement this as better alternative
                        ! see util_mod.F90 most_frequent_nonzero
                        ! nval = 0
                        ! values = 0
                        ! DO ip = -1, MIN(ii-i, 1)
                        !     DO jp = -1, MIN(jj-j, 1)
                        !         DO kp = -1, MIN(kk-k, 1)
                        !             IF (bzelltyp(k+kp, j+jp, i+ip) >= 0) CYCLE
                        !             icell_lauf = -bzelltyp(k+kp, j+jp, i+ip)
                        !             nval = nval + 1
                        !             values(nval) = bodyid(icell_lauf)
                        !         END DO
                        !     END DO
                        ! END DO
                        ! bodyid(icell) = most_frequent_nonzero(values(1:nval))

                        DO ip = MAX(2-i, -1), MIN(ii-i, 1)
                            DO jp = MAX(2-j, -1), MIN(jj-j, 1)
                                DO kp = MAX(2-k, -1), MIN(kk-k, 1)
                                    IF (bzelltyp(k+kp, j+jp, i+ip) >= 0) CYCLE
                                    icell_lauf = -bzelltyp(k+kp, j+jp, i+ip)
                                    bodyid(icell) = MAX(bodyid(icell), &
                                        bodyid(icell_lauf))
                                END DO
                            END DO
                        END DO

                        ! If bodyid of cell is not zero anymore, it was fixed.
                        ! Count it.
                        IF (bodyid(icell) /= 0) nfixed_cells = nfixed_cells + 1
                    END IF
                END DO
            END DO
        END DO
    END SUBROUTINE fixup_bodyid
END MODULE calcnormals_mod
