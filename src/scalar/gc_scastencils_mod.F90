MODULE gc_scastencils_mod
    USE core_mod
    USE ib_mod
    USE scacore_mod

    USE gc_createstencils_mod, ONLY: wmcheckneighbor, choosestencil, &
        wmindexlistn

    IMPLICIT NONE (type, external)
    PRIVATE

    TYPE(int_stencils_t), ALLOCATABLE, TARGET :: scaxpoli(:)
    TYPE(real_stencils_t), ALLOCATABLE, TARGET :: scaxpolr(:), scaxpolrvel(:)
    INTEGER(intk), ALLOCATABLE :: scanblg(:)

    PUBLIC :: create_scastencils, finish_scastencils, set_scastencils

CONTAINS
    SUBROUTINE create_scastencils(gc)

        ! Subroutine arguments
        TYPE(gc_t), INTENT(in) :: gc

        ! Local variables
        INTEGER(intk) :: ilevel
        REAL(realk), POINTER, CONTIGUOUS :: bp(:)

        CALL get_fieldptr(bp, "BP")

        ! Always allocate - createstencils should be called only once.
        ALLOCATE(scaxpoli(nmygrids))
        ALLOCATE(scaxpolr(nmygrids))
        ALLOCATE(scaxpolrvel(nmygrids))
        ALLOCATE(scanblg(nmygrids))

        DO ilevel = minlevel, maxlevel
            CALL createstencils_level(ilevel, bp, gc%bzelltyp, gc%icells, &
                gc%icellspointer, gc%nvecs, gc%arealist, gc%bodyid)
        END DO

        ! mglet_dbg_envvar is in buildinfo_mod and initialized at startup
        IF (INDEX(mglet_dbg_envvar, "stencilvtk") > 0) THEN
            CALL writestencils()
        END IF
    END SUBROUTINE create_scastencils


    SUBROUTINE finish_scastencils()
        DEALLOCATE(scaxpoli)
        DEALLOCATE(scaxpolr)
        DEALLOCATE(scaxpolrvel)
        DEALLOCATE(scanblg)
    END SUBROUTINE finish_scastencils


    SUBROUTINE createstencils_level(ilevel, bp, bzelltyp, icells, &
            icellspointer, nvecs, arealist, bodyid)

        ! Subroutine arguments
        INTEGER(intk), INTENT(in) :: ilevel
        REAL(realk), INTENT(in) :: bp(*)
        INTEGER(intk), INTENT(in) :: bzelltyp(*)
        INTEGER(intk), CONTIGUOUS, INTENT(in) :: icells(:)
        INTEGER(intk), CONTIGUOUS, INTENT(in) :: icellspointer(:)
        REAL(realk), CONTIGUOUS, INTENT(in) :: nvecs(:, :)
        REAL(realk), CONTIGUOUS, INTENT(in) :: arealist(:)
        INTEGER(intk), CONTIGUOUS, INTENT(in) :: bodyid(:)

        ! Local variables
        INTEGER(intk) :: i, igrid, kk, jj, ii, ip3, ipp, ncells
        INTEGER(intk) :: bconds(6)
        REAL(realk), POINTER, CONTIGUOUS :: x(:), y(:), z(:)
        REAL(realk), POINTER, CONTIGUOUS :: xstag(:), ystag(:), zstag(:)
        REAL(realk), POINTER, CONTIGUOUS :: ddx(:), ddy(:), ddz(:)

        DO i = 1, nmygridslvl(ilevel)
            igrid = mygridslvl(i, ilevel)
            ipp = icellspointer(igrid)
            ncells = icells(igrid)

            CALL get_fieldptr(x, "X", igrid)
            CALL get_fieldptr(y, "Y", igrid)
            CALL get_fieldptr(z, "Z", igrid)

            CALL get_fieldptr(xstag, "XSTAG", igrid)
            CALL get_fieldptr(ystag, "YSTAG", igrid)
            CALL get_fieldptr(zstag, "ZSTAG", igrid)

            CALL get_fieldptr(ddx, "DDX", igrid)
            CALL get_fieldptr(ddy, "DDY", igrid)
            CALL get_fieldptr(ddz, "DDZ", igrid)

            CALL get_mgdims(kk, jj, ii, igrid)
            CALL get_ip3(ip3, igrid)

            CALL get_mgbasb(bconds, igrid)

            CALL tscastencil(igrid, kk, jj, ii, x, y, z, xstag, ystag, zstag, &
                ddx, ddy, ddz, bp(ip3), bzelltyp(ip3), bconds, ncells, &
                nvecs(:, ipp:ipp+ncells-1), arealist(ipp:ipp+ncells-1), &
                bodyid(ipp:ipp+ncells-1))
        END DO
    END SUBROUTINE createstencils_level


    SUBROUTINE tscastencil(igrid, kk, jj, ii, x, y, z, xstag, ystag, zstag, &
            ddx, ddy, ddz, bp, bzelltyp, bconds, icells, nvecs, arealist, &
            bodyid)

        USE findinterface_mod, ONLY: findinterface2

        ! Subroutine arguments
        INTEGER(intk), INTENT(in) :: igrid
        INTEGER(intk), INTENT(in) :: kk, jj, ii
        REAL(realk), INTENT(in) :: x(ii), y(jj), z(kk)
        REAL(realk), INTENT(in) :: xstag(ii), ystag(jj), zstag(kk)
        REAL(realk), INTENT(in) :: ddx(ii), ddy(jj), ddz(kk)
        REAL(realk), INTENT(in) :: bp(kk, jj, ii)
        INTEGER(intk), INTENT(in) :: bzelltyp(kk, jj, ii)
        INTEGER(intk), INTENT(in) :: bconds(6)
        INTEGER(intk), INTENT(in) :: icells
        REAL(realk), CONTIGUOUS, INTENT(in) :: nvecs(:, :)
        REAL(realk), CONTIGUOUS, INTENT(in) :: arealist(:)
        INTEGER(intk), CONTIGUOUS, INTENT(in) :: bodyid(:)

        ! Local variables
        INTEGER(intk), PARAMETER :: lsize = 122
        INTEGER(intk), PARAMETER :: velpts = 6
        INTEGER(intk), PARAMETER :: ldofa = 20

        INTEGER(intk) :: i, j, k, it, jt, kt
        INTEGER(intk) :: pntxpoli, pntxpolr, counter, ngeomnbr
        INTEGER(intk) :: found, found2, found3, foundnr
        INTEGER(intk) :: foundx1, foundx2, foundy1, foundy2, foundz1, foundz2
        INTEGER(intk) :: xpolisize, xpolrsize, imygrid
        INTEGER(intk) :: inlst(lsize), jnlst(lsize), knlst(lsize), nnlst
        INTEGER(intk) :: body
        REAL(realk) :: area, areabyvol
        INTEGER(intk), ALLOCATABLE :: xpoli(:)
        REAL(realk), ALLOCATABLE :: xpolr(:), xpolrvel(:)

        ! Allocating the arrays with possible overhead
        xpolisize = NINT(icells*1.1 + MAX(kk, jj, ii)**2)*(2+velpts)
        xpolrsize = NINT(icells*1.1 + MAX(kk, jj, ii)**2)*velpts
        ALLOCATE(xpoli(xpolisize))
        ALLOCATE(xpolr(xpolrsize))
        ALLOCATE(xpolrvel(xpolrsize))

        ! Initialization
        counter = 0
        pntxpoli = 0
        xpoli = 0
        pntxpolr = 0.0
        xpolr = 0.0

        ! Index stencil for looking for neighbour cells.
        CALL wmindexlistn(inlst, jnlst, knlst, nnlst)

        DO i = 3, ii-2
            DO j = 3, jj-2
                DO k = 3, kk-2

                    CALL findinterface2(k, j, i, kk, jj, ii, bp, found, &
                        foundx1, foundx2, foundy1, foundy2, foundz1, &
                        foundz2, foundnr)

                    ! If interface cell not found, skip rest of loop
                    IF (found == 0) THEN
                        ! Cell itself is open but was marked as interface cell
                        IF (bzelltyp(k, j, i) < 0 .AND. bp(k, j, i) > 0.5) THEN
                            counter = counter + 1

                            ! Retrieving the values for surface-intersected cell
                            area = arealist(-bzelltyp(k, j, i))
                            areabyvol = area/(ddx(i)*ddy(j)*ddz(k))
                            body = bodyid(-bzelltyp(k, j, i))

                            ! Creating a stencil that only stores area for flux
                            CALL tscastencilcoeffarea(k, j, i, kk, jj, ii, &
                                pntxpoli, xpoli, pntxpolr, xpolr, xpolrvel, &
                                areabyvol, body)
                        END IF
                    ELSE
                        counter = counter + 1

                        found2 = 0
                        IF (bzelltyp(k, j, i) < 0) THEN
                            ! bzelltyp < 0: cell has information about the
                            ! boundary condition
                            it = i
                            jt = j
                            kt = k
                            found2 = 1
                        ELSE
                            ! looking of neighboring cells have information
                            ngeomnbr = 0
                            if (foundx1 == 1 .AND. bzelltyp(k, j, i+1) < 0) &
                                ngeomnbr = ngeomnbr + 1
                            if (foundx2 == 1 .AND. bzelltyp(k, j, i-1) < 0) &
                                ngeomnbr = ngeomnbr + 1
                            if (foundy1 == 1 .AND. bzelltyp(k, j+1, i) < 0) &
                                ngeomnbr = ngeomnbr + 1
                            if (foundy2 == 1 .AND. bzelltyp(k, j-1, i) < 0) &
                                ngeomnbr = ngeomnbr + 1
                            if (foundz1 == 1 .AND. bzelltyp(k+1, j, i) < 0) &
                                ngeomnbr = ngeomnbr + 1
                            if (foundz2 == 1 .AND. bzelltyp(k-1, j, i) < 0) &
                                ngeomnbr = ngeomnbr + 1

                            IF (ngeomnbr > 0) THEN
                                ! there is information in (at least) one
                                ! neighboring cells. In case of multiple
                                ! eligible cells, the biggest bodyid is
                                ! selected.
                                body = 0
                                IF (foundx1 == 1 .AND. bzelltyp(k, j, i+1) < 0) THEN
                                    IF (bodyid(-bzelltyp(k, j, i+1)) > body) THEN
                                        body = bodyid(-bzelltyp(k, j, i+1))
                                        it = i+1
                                        jt = j
                                        kt = k
                                    END IF
                                END IF
                                IF (foundx2 == 1 .AND. bzelltyp(k, j, i-1) < 0) THEN
                                    IF (bodyid(-bzelltyp(k, j, i-1)) > body) THEN
                                        body = bodyid(-bzelltyp(k, j, i-1))
                                        it = i-1
                                        jt = j
                                        kt = k
                                    END IF
                                END IF
                                IF (foundy1 == 1 .AND. bzelltyp(k, j+1, i) < 0) THEN
                                    IF (bodyid(-bzelltyp(k, j+1, i)) > body) THEN
                                        body = bodyid(-bzelltyp(k, j+1, i))
                                        it = i
                                        jt = j+1
                                        kt = k
                                    END IF
                                END IF
                                IF (foundy2 == 1 .AND. bzelltyp(k, j-1, i) < 0) THEN
                                    IF (bodyid(-bzelltyp(k, j-1, i  )) > body) THEN
                                        body = bodyid(-bzelltyp(k, j-1, i))
                                        it = i
                                        jt = j-1
                                        kt = k
                                    END IF
                                END IF
                                IF (foundz1 == 1 .AND. bzelltyp(k+1, j, i) < 0) THEN
                                    IF (bodyid(-bzelltyp(k+1, j, i)) > body) THEN
                                        body = bodyid(-bzelltyp(k+1, j, i))
                                        it = i
                                        jt = j
                                        kt = k+1
                                    END IF
                                END IF
                                IF (foundz2 == 1 .AND. bzelltyp(k-1, j, i) < 0) THEN
                                    IF (bodyid(-bzelltyp(k-1, j, i)) > body) THEN
                                        body = bodyid(-bzelltyp(k-1, j, i))
                                        it = i
                                        jt = j
                                        kt = k-1
                                    END IF
                                END IF
                            END IF

                            ! A valid bodyid was found
                            IF (body > 0) found2 = 1
                        END IF

                        IF (found2 == 1) THEN
                            ! Check if area comes from this cell
                            IF (k == kt .AND. j == jt .AND. i == it) THEN
                                area = arealist(-bzelltyp(kt, jt, it))
                                areabyvol = area/(ddx(i)*ddy(j)*ddz(k))
                            ELSE
                                areabyvol = 0.0
                            END IF
                            body = bodyid(-bzelltyp(kt, jt, it))

                            ! Computation of the "proper" stencil for scalar
                            CALL tscastencilcoeff(k, j, i, kt, jt, it, kk, jj, &
                                ii, x, y, z, xstag, ystag, zstag, ddx, ddy, &
                                ddz, bp, bzelltyp, bconds, nnlst, inlst, &
                                jnlst, knlst, velpts, nvecs, ldofa, areabyvol, &
                                body, pntxpoli, xpoli, pntxpolr, xpolr, &
                                xpolrvel, found3)

                            ! Check successful generation of a stencil by
                            ! tscastencilcoeff
                            IF (found3 <= 0) THEN
                                CALL tscastencilcoeffrescue(k, j, i, kk, jj, &
                                    ii, pntxpoli, xpoli, pntxpolr, xpolr, &
                                    xpolrvel, areabyvol, body)
                            END IF
                        ELSE
                           ! rescuestencil, because no triangle was found
                            ! (emergency solution)
                            areabyvol = 0.0
                            body = 0
                            CALL tscastencilcoeffrescue(k, j, i, kk, jj, &
                                ii, pntxpoli, xpoli, pntxpolr, xpolr, &
                                xpolrvel, areabyvol, body)
                        END IF
                    END IF

                END DO
            END DO
        END DO

        CALL get_imygrid(imygrid, igrid)
        ALLOCATE(scaxpoli(imygrid)%arr(pntxpoli))
        ALLOCATE(scaxpolr(imygrid)%arr(pntxpolr))
        ALLOCATE(scaxpolrvel(imygrid)%arr(pntxpolr))

        scaxpoli(imygrid)%arr = xpoli(1:pntxpoli)
        scaxpolr(imygrid)%arr = xpolr(1:pntxpolr)
        scaxpolrvel(imygrid)%arr = xpolrvel(1:pntxpolr)
        scanblg(imygrid) = counter

        DEALLOCATE(xpoli)
        DEALLOCATE(xpolr)
        DEALLOCATE(xpolrvel)
    END SUBROUTINE tscastencil


    SUBROUTINE tscastencilcoeffarea(k, j, i, kk, jj, ii, pntxpoli, xpoli, &
            pntxpolr, xpolr, xpolrvel, areabyvol, body)

        ! Creation of special stencil, that only considers fluxes.
        ! The intersected area divided by the open cell volume
        ! ddx(i)*ddy(j)*ddz(k) is stored.
        ! All coefficients of Dirichlet parameters are set such that value
        ! is not changed.

        ! Subroutine arguments
        INTEGER(intk), INTENT(in) :: k, j, i, kk, jj, ii
        INTEGER(intk), INTENT(inout) :: xpoli(:)
        REAL(realk), INTENT(inout) :: xpolr(:)
        REAL(realk), INTENT(inout) :: xpolrvel(:)
        INTEGER(intk), INTENT(inout) :: pntxpoli, pntxpolr
        REAL(realk), intent(in) :: areabyvol
        INTEGER(intk), INTENT(in) :: body

        ! Local variables
        INTEGER(intk), PARAMETER :: velpts = 1
        INTEGER(intk) :: xpolisize, xpolrsize
        INTEGER(intk) :: cellind

        xpolisize = SIZE(xpoli)
        xpolrsize = SIZE(xpolr)

        IF (pntxpoli+velpts+2 > xpolisize) THEN
            CALL errr(__FILE__, __LINE__)
        END IF
        IF (pntxpolr+velpts+1 > xpolrsize) THEN
            CALL errr(__FILE__, __LINE__)
        END IF

        ! TODO: sub2ind
        cellind = 1+(k-1)+(j-1)*kk+(i-1)*jj*kk

        ! 1st integer: linearized index of cell to be interpolated
        pntxpoli = pntxpoli + 1
        xpoli(pntxpoli) = cellind

        ! 2nd integer: body id of this stencils
        pntxpoli = pntxpoli + 1
        xpoli(pntxpoli) = body

        ! 3rd integer: number of the peripheral stencil points (N)
        pntxpoli = pntxpoli + 1
        xpoli(pntxpoli) = velpts  ! = 1

        ! 1st real: setting area/volume of the interface within the cell
        ! (additional in comparison to the velocity stencil)
        pntxpolr = pntxpolr + 1
        xpolr(pntxpolr) = areabyvol
        xpolrvel(pntxpolr) = areabyvol

        ! -- replaces iteration over velpoints --
        pntxpoli = pntxpoli + 1
        xpoli(pntxpoli) = cellind

        pntxpolr = pntxpolr + 1
        xpolr(pntxpolr) = 1.0
        xpolrvel(pntxpolr) = 1.0

        ! 3rd real: setting the additive constant
        pntxpolr = pntxpolr + 1
        xpolr(pntxpolr) = 0.0
        xpolrvel(pntxpolr) =  0.0
    END SUBROUTINE tscastencilcoeffarea


    SUBROUTINE tscastencilcoeff(k, j, i, kt, jt, it, kk, jj, ii, x, y, z, &
        xstag, ystag, zstag, ddx, ddy, ddz, bp, bzelltyp, bconds, nnlst, &
        inlst, jnlst, knlst, velpts, nvecs, ldofa, areabyvol, body, pntxpoli, &
        xpoli, pntxpolr, xpolr, xpolrvel, found)

        ! NOTE: This subroutine is very close to calcflux from
        ! gc_flowstencils_mod.F90

        ! Subroutine arguments
        INTEGER(intk), INTENT(in) :: k, j, i
        INTEGER(intk), INTENT(in) :: kt, jt, it
        INTEGER(intk), INTENT(in) :: kk, jj, ii
        REAL(realk), INTENT(in) :: x(ii), y(jj), z(kk)
        REAL(realk), INTENT(in) :: xstag(ii), ystag(jj), zstag(kk)
        REAL(realk), INTENT(in) :: ddx(ii), ddy(jj), ddz(kk)
        REAL(realk), INTENT(in) :: bp(kk, jj, ii)
        INTEGER(intk), INTENT(in) :: bzelltyp(kk, jj, ii)
        INTEGER(intk), INTENT(in) :: bconds(6)
        INTEGER(intk), INTENT(in) :: nnlst
        INTEGER(intk), CONTIGUOUS, INTENT(in) :: inlst(:), jnlst(:), knlst(:)
        INTEGER(intk), INTENT(in) :: velpts
        REAL(realk), CONTIGUOUS, INTENT(in) :: nvecs(:, :)
        INTEGER(intk), INTENT(in) :: ldofa
        REAL(realk), INTENT(in) :: areabyvol
        INTEGER(intk), INTENT(in) :: body
        INTEGER(intk), INTENT(inout) :: pntxpoli
        INTEGER(intk), CONTIGUOUS, INTENT(inout) :: xpoli(:)
        INTEGER(intk), INTENT(inout) :: pntxpolr
        REAL(realk), CONTIGUOUS, INTENT(inout) :: xpolr(:)
        REAL(realk), CONTIGUOUS, INTENT(inout) :: xpolrvel(:)
        INTEGER(intk), INTENT(inout) :: found

        ! Local variables
        INTEGER(intk) :: istc(SIZE(inlst)), jstc(SIZE(jnlst)), kstc(SIZE(knlst))
        INTEGER(intk) :: nstc, chosenvelpts
        REAL(realk) :: coeffstc(ldofa), acoeffstc
        REAL(realk) :: coeffstcvel(ldofa), acoeffstcvel
        INTEGER(intk) :: mm

        ! Find non-blocked neighbors.
        CALL wmcheckneighbor(k, j, i, kk, jj, ii, bp, bconds, nnlst, &
            inlst, jnlst, knlst, istc, jstc, kstc, nstc)

        ! Auswaehlen des Stencils
        CALL choosestencil(k, j, i, nstc, istc, jstc, kstc, chosenvelpts)
        IF (chosenvelpts > velpts) THEN
            found = -2
            RETURN
        END IF

        ! keine offenen Nachbarn: rescue-Stencil
        IF (chosenvelpts == 0) THEN
            found = -3
            RETURN
        END IF

        CALL calcfluxtsca(k, j, i, kt, jt, it, kk, jj, ii, &
            x, y, z, xstag, ystag, zstag, ddx, ddy, ddz, bzelltyp, &
            istc, jstc, kstc, chosenvelpts, nvecs, coeffstc, acoeffstc, &
            coeffstcvel, acoeffstcvel, found)

        ! catching bad values, rescue-Stencil (similar to velocity stencil)
        DO mm = 1, chosenvelpts
            IF (ABS(coeffstcvel(mm)) > 12.0) THEN
                found = -4
                RETURN
            END IF
        END DO

        ! --- adding to the coefficient list ---
        IF (found == 1) THEN
            CALL tscastencilcoefffull(k, j, i, kk, jj, ii, pntxpoli, xpoli, &
                pntxpolr, xpolr, xpolrvel, chosenvelpts, istc, jstc, kstc, &
                coeffstc, acoeffstc, coeffstcvel, acoeffstcvel, areabyvol, body)
        END IF
    END SUBROUTINE tscastencilcoeff


    SUBROUTINE calcfluxtsca(k, j, i, kt, jt, it, kk, jj, ii, x, y, z, &
            xstag, ystag, zstag, ddx, ddy, ddz, bzelltyp, istc, jstc, kstc, &
            velpts, nvecs, coeffstc, acoeffstc, coeffstctsca, &
            acoeffstctsca, found)

        !  Function computes the interpolation coefficients for the stencils:
        !  - coeffstc = required for scalar cell value (these coefficients are
        !    multiplied with the field values at peripheral points of the
        !    stencil)
        !  - acoeffstc = required for scalar cell value (these coefficients are
        !    added to the resulting expression, a = additive)
        !  - coeffstctsca = required for scalar point value (these coefficients
        !    are multiplied with the field values at peripheral points of the
        !    stencil)
        !  - acoeffstctsca = required for scalar point value (these
        !    coefficients are added to the resulting expression, a = additive)

        ! NOTE: This subroutine is very close to wmmultimatrix from
        ! gc_createstencils_mod.F90

        ! Subroutine arguments
        INTEGER(intk), INTENT(in) :: k, j, i
        INTEGER(intk), INTENT(in) :: kt, jt, it
        INTEGER(intk), INTENT(in) :: kk, jj, ii
        REAL(realk), INTENT(in) :: x(ii), y(jj), z(kk)
        REAL(realk), INTENT(in) :: xstag(ii), ystag(jj), zstag(kk)
        REAL(realk), INTENT(in) :: ddx(ii), ddy(jj), ddz(kk)
        INTEGER(intk), INTENT(in) :: bzelltyp(kk, jj, ii)
        INTEGER(intk), CONTIGUOUS, INTENT(in) :: istc(:), jstc(:), kstc(:)
        INTEGER(intk), INTENT(in) :: velpts
        REAL(realk), CONTIGUOUS, INTENT(in) :: nvecs(:, :)
        REAL(realk), CONTIGUOUS, INTENT(out) :: coeffstc(:)
        REAL(realk), INTENT(out) :: acoeffstc
        REAL(realk), CONTIGUOUS, INTENT(out) :: coeffstctsca(:)
        REAL(realk), INTENT(out) :: acoeffstctsca
        INTEGER(intk), INTENT(out) :: found

        ! Local variables
        INTEGER(intk), PARAMETER :: parts = 20
        INTEGER(intk) :: kt2, jt2, it2
        INTEGER(intk) :: n, nanstenc
        INTEGER(intk) :: idx
        INTEGER(intk) :: cvol1, cvol2     ! counters (1=fluid, 2=solid)
        REAL(realk) :: a, b, c, dpl
        REAL(realk) :: xs, ys, zs
        REAL(realk) :: ddi, ddj, ddk, xm, ym, zm, dvol, dir, sum
        REAL(realk) :: dist(velpts), ce(velpts), cesum, dall, dpoint

        ! Same as in the beginning of wmmultimatrix
        idx = -bzelltyp(kt, jt, it)
        a = nvecs(1, idx)
        b = nvecs(2, idx)
        c = nvecs(3, idx)
        dpl = nvecs(4, idx)

        coeffstc = 0.0
        acoeffstc = 0.0

        coeffstctsca = 0.0
        acoeffstctsca = 0.0

        found = -11

        ! Distance of cell center to plane (= interface of solid body)
        dpoint = a*x(i) + b*y(j) + c*z(k) - dpl

        dall = 0.0
        dist = 0.0

        DO n = 1, velpts
            ! Computing distance for each peripheral stencil point
            it2 = istc(n)
            jt2 = jstc(n)
            kt2 = kstc(n)
            dist(n) = a*x(it2) + b*y(jt2) + c*z(kt2) - dpl

            ! setting to zero if significantly smaller than central point
            ! distance
            ! JK: Fraglich, ob das fuer den Skalar sinnvoll ist.
            IF (dist(n) < ABS(dpoint)*0.2) THEN
                dist(n) = 0.0
            END IF

            ! summing up the (corrected) distances
            dall = dall + dist(n)**2.0
        END DO

        ! if (a) normal vector has absolute of 0 or
        !    (b) all points for the stencil are on the wrong side

        IF (ABS(dall) < TINY(0.0_realk)) THEN
            ! do not generate any stencil
            ce = 0.0
            DO n = 1, velpts
                ! computing distance for each peripheral stencil point
                it2 = istc(n)
                jt2 = jstc(n)
                kt2 = kstc(n)
                dist(n) = a*x(it2) + b*y(jt2) + c*z(kt2) - dpl
            END DO
            RETURN
        END IF

        cesum = 0.0
        DO n = 1, velpts
            ce(n) = dist(n)/dall
            cesum = cesum + ce(n)
        END DO

        ! TODO: replace the volume integration here with something better,
        ! preferably from the blocking/stencils.
        ddi = ddx(i)/parts
        ddj = ddy(j)/parts
        ddk = ddz(k)/parts

        xs = xstag(i-1) - ddi/2.0
        ys = ystag(j-1) - ddj/2.0
        zs = zstag(k-1) - ddk/2.0

        cvol1 = 0
        cvol2 = 0

        ! Volume of one single small part of the cell
        dvol = ddi*ddj*ddk

        ! iteration over all small dV = dx*dy*dz (volume integration)
        DO it2 = 1, parts
            xm = xs + ddi*REAL(it2, realk)
            DO jt2 = 1, parts
                ym = ys + ddj*REAL(jt2, realk)
                DO kt2 = 1, parts
                    zm = zs + ddk*REAL(kt2, realk)
                    ! (xm,ym,zm) = center of a small volume part dx * dy * dz

                    dir = (a*xm + b*ym + c*zm) - dpl
                    ! distance of this part from the plane
                    ! (serves also as a direction indicator - inside or outside)

                    IF (dir > 0.0) THEN
                        ! part dV lies in fluid
                        dpoint = dir
                        DO n = 1, velpts
                            coeffstc(n) = coeffstc(n) + dpoint * ce(n)
                        END DO
                        ! multiplied with tval later
                        acoeffstc = acoeffstc + (1.0 - dpoint*cesum)
                        cvol1 = cvol1 + 1   ! fluid counter++
                    ELSE
                        ! part dV lies in solid
                        cvol2 = cvol2 + 1   ! solid counter++
                    END IF
                END DO
            END DO
        END DO

        !!! --- Cell value stencil --- (comparable to flux stencil)
        ! stencil is applied to set the mean value of the fluid-filled part
        ! of the cell as the cell value for proper statistics and
        ! postprocessing (masked with VOLP)
        IF (cvol1 > 0) THEN
            ! multiplicative constants (=> fluid volume mean)
            DO n = 1, velpts
                coeffstc(n) = coeffstc(n) * dvol / (cvol1*dvol)
            END DO

            ! additive constant (=> fluid volume mean)
            ! IMPORTANT: this additive constant is NOT yet multiplied with
            ! the Dirichlet value
            acoeffstc = ( acoeffstc * dvol / (cvol1 * dvol) )
        ELSE
            ! multiplicative constants
            DO n = 1, velpts
                coeffstc(n) = 0.0
            END DO
            ! additive constant
            acoeffstc = 0.0
        END IF

        ! --- Point value stencil ---
        ! (required for the scalar point interpolation)
        ! Stencil is applied to set a second-order accurate point value in
        ! the Ghost Cell such that the demanded value is met exactly on
        ! at the interface position assuming linear interpolation
        dpoint = a*x(i) + b*y(j) + c*z(k) - dpl
        DO n = 1, velpts
            coeffstctsca(n) = dpoint * ce(n)
        END DO

        ! IMPORTANT: this additive constant is NOT yet multiplied with the
        ! Dirichlet value
        acoeffstctsca = 1.0 - dpoint*cesum

        ! Stencil generated
        found = 1

        ! Removing the invalid coefficients
        sum = 0.0
        DO n = 1, velpts
            sum = sum + coeffstc(n)
        END DO

        ! jk 7.2.2008
        ! Grenze auf 3/5: Wert fuer Wuerfelgitter, wenn Ebene durch Ecke geht

        IF (sum > 0.61) THEN
            found = -12
        END IF

        ! test coeffient for NAN
        nanstenc = 0
        DO n = 1,velpts
            IF (ABS(coeffstctsca(n)) < 0.0) THEN
                nanstenc = 1
            END IF
            IF (ABS(coeffstc(n)) < 0.0) THEN
                nanstenc = 1
            END IF
        END DO
        IF (ABS(acoeffstctsca) < 0.0) THEN
            nanstenc = 1
        END IF
        IF (ABS(acoeffstc) < 0.0) THEN
            nanstenc = 1
        END IF

        IF (nanstenc == 1) THEN
            WRITE(*,*) "-----------------------------------------------"
            WRITE(*,*) "i, j, k", i, j, k
            DO n = 1, velpts
                WRITE(*,*) "dist(n)", n, dist(n)
            END DO
            WRITE(*,*) "dall",dall
            DO n = 1, velpts
                write(*,*) "ce(n)", n, ce(n)
            END DO
            WRITE(*,*) "cesum", cesum
            WRITE(*,*) "velpts", velpts
            WRITE(*,*) "a, b, c, dpl", a, b, c, dpl
            WRITE(*,*) "dpoint", dpoint
            WRITE(*,*)
            DO n = 1, velpts
                write(*,*) "coeffstc(n)", n, coeffstc(n)
            END DO
            WRITE(*,*) "acoeffstc", acoeffstc
            WRITE(*,*)
            DO n = 1, velpts
                write(*,*) "coeffstctsca(n)", n, coeffstctsca(n)
            END DO
            WRITE(*,*) "acoeffstctsca", acoeffstctsca
            WRITE(*,*) "sr calcfluxtsca: acoeffstc", acoeffstc
            WRITE(*,*) "cvol1, cvol2", cvol1, cvol2
            WRITE(*,*) "-----------------------------------------------"

            CALL errr(__FILE__, __LINE__)
        END IF
    END SUBROUTINE calcfluxtsca


    SUBROUTINE tscastencilcoefffull(k, j, i, kk, jj, ii, pntxpoli, xpoli, &
            pntxpolr, xpolr, xpolrvel, velpts, istc, jstc, kstc, coeffstc, &
            acoeffstc, coeffstcvel, acoeffstcvel, areabyvol, body)

        ! Subroutine arguments
        INTEGER(intk), INTENT(in) :: k, j, i
        INTEGER(intk), INTENT(in) :: kk, jj, ii
        INTEGER(intk), INTENT(inout) :: pntxpoli
        INTEGER(intk), INTENT(inout) :: xpoli(:)
        INTEGER(intk), INTENT(inout) :: pntxpolr
        REAL(realk), INTENT(inout) :: xpolr(:), xpolrvel(:)
        INTEGER(intk), INTENT(in) :: velpts
        INTEGER(intk), INTENT(in) :: istc(:), jstc(:), kstc(:)
        REAL(realk), INTENT(in) :: coeffstc(:), acoeffstc
        REAL(realk), INTENT(in) :: coeffstcvel(:), acoeffstcvel
        REAL(realk), INTENT(in) :: areabyvol
        INTEGER(intk), INTENT(in) :: body

        ! Local variables
        INTEGER(intk) :: n
        INTEGER(intk) :: xpolisize, xpolrsize

        xpolisize = SIZE(xpoli)
        xpolrsize = SIZE(xpolr)

        IF (pntxpoli+velpts+3 > xpolisize) THEN
            CALL errr(__FILE__, __LINE__)
        END IF
        IF (pntxpolr+velpts+2.gt.xpolrsize) THEN
            CALL errr(__FILE__, __LINE__)
        END IF

        ! TODO: sub2ind
        pntxpoli = pntxpoli + 1
        xpoli(pntxpoli) = 1+(k-1)+(j-1)*kk+(i-1)*jj*kk

        ! 2nd integer: body id of this stencils
        pntxpoli = pntxpoli + 1
        xpoli(pntxpoli) = body

        ! 3rd integer: number of the peripheral stencil points (N)
        pntxpoli = pntxpoli + 1
        xpoli(pntxpoli) = velpts

        ! 1st real: setting area/volume of the interface within the cell
        ! (additional in comparison to the velocity stencil)
        pntxpolr = pntxpolr + 1
        xpolr(pntxpolr) = areabyvol
        xpolrvel(pntxpolr) = areabyvol

        ! here with iteration over all N peripheral stencil points
        ! (velpts here carries the value of "chosenvelpts")
        DO n = 1, velpts
            pntxpoli = pntxpoli + 1
            xpoli(pntxpoli) = 1+(kstc(n)-1)+(jstc(n)-1)*kk+(istc(n)-1)*jj*kk

            pntxpolr = pntxpolr + 1
            xpolr(pntxpolr) = coeffstc(n)
            xpolrvel(pntxpolr) = coeffstcvel(n)
        END DO

        ! last real: setting the additive constant
        pntxpolr = pntxpolr + 1
        xpolr(pntxpolr) = acoeffstc
        xpolrvel(pntxpolr) = acoeffstcvel
    END SUBROUTINE tscastencilcoefffull


    SUBROUTINE tscastencilcoeffrescue(k, j, i, kk, jj, ii, pntxpoli, xpoli, &
            pntxpolr, xpolr, xpolrvel, areabyvol, body)

        ! Subroutine arguments
        INTEGER(intk), INTENT(in) :: k, j, i, kk, jj, ii
        INTEGER(intk), INTENT(inout) :: xpoli(:)
        REAL(realk), INTENT(inout) :: xpolr(:)
        REAL(realk), INTENT(inout) :: xpolrvel(:)
        INTEGER(intk), INTENT(inout) :: pntxpoli, pntxpolr
        REAL(realk), intent(in) :: areabyvol
        INTEGER(intk), INTENT(in) :: body

        ! Local variables
        INTEGER(intk), PARAMETER :: velpts = 0
        INTEGER(intk) :: xpolisize, xpolrsize

        xpolisize = SIZE(xpoli)
        xpolrsize = SIZE(xpolr)

        IF (pntxpoli+velpts+2 > xpolisize) THEN
            CALL errr(__FILE__, __LINE__)
        END IF
        IF (pntxpolr+velpts+1 > xpolrsize) THEN
            CALL errr(__FILE__, __LINE__)
        END IF

        ! increasing the index counter only by one
        ! (no iteration over points as in function above)

        ! TODO: sub2ind
        pntxpoli = pntxpoli + 1
        xpoli(pntxpoli) = 1+(k-1)+(j-1)*kk+(i-1)*jj*kk

        ! 2nd integer: body id of this stencils
        pntxpoli = pntxpoli + 1
        xpoli(pntxpoli) = body

        ! 3rd integer: number of the peripheral stencil points (N)
        pntxpoli = pntxpoli + 1
        xpoli(pntxpoli) = velpts

        ! 1st real: setting area/volume of the interface within the cell
        ! (additional in comparison to the velocity stencil)
        pntxpolr = pntxpolr + 1
        xpolr(pntxpolr) = areabyvol
        xpolrvel(pntxpolr) = areabyvol

        ! -- no iteration over velpoints --

        ! 3rd real: setting the additive constant (here: equal to the final cell value)
        pntxpolr = pntxpolr + 1
        xpolr(pntxpolr) = 1.0
        xpolrvel(pntxpolr) =  1.0
    END SUBROUTINE tscastencilcoeffrescue


    SUBROUTINE set_scastencils(ctyp, sca, t, qtt)
        ! Subroutine arguments
        CHARACTER(len=1), INTENT(in) :: ctyp
        TYPE(scalar_t), INTENT(in) :: sca
        TYPE(field_t), INTENT(inout), OPTIONAL :: t, qtt

        ! Local variables
        INTEGER(intk) :: i, igrid
        INTEGER(intk) :: kk, jj, ii
        REAL(realk), POINTER, CONTIGUOUS :: t_p(:, :, :), qtt_p(:, :, :)
        REAL(realk), POINTER, CONTIGUOUS :: xpolr(:)

        CALL start_timer(412)

        NULLIFY(t_p)
        NULLIFY(qtt_p)

        DO i = 1, nmygrids
            igrid = mygrids(i)

            CALL get_mgdims(kk, jj, ii, igrid)

            IF (PRESENT(t)) CALL t%get_ptr(t_p, igrid)
            IF (PRESENT(qtt)) CALL qtt%get_ptr(qtt_p, igrid)

            SELECT CASE(ctyp)
            CASE("C")
                xpolr => scaxpolr(i)%arr
            CASE("P")
                xpolr => scaxpolrvel(i)%arr
            CASE DEFAULT
                CALL errr(__FILE__, __LINE__)
            END SELECT

            CALL wmxpolsoltsca(kk, jj, ii, scanblg(i), sca, scaxpoli(i)%arr, &
                xpolr, t_p, qtt_p)
        END DO

        CALL stop_timer(412)
    END SUBROUTINE set_scastencils


    SUBROUTINE wmxpolsoltsca(kk, jj, ii, ncells, sca, xpoli, xpolr, t, qtt)
        ! Subroutine arguments
        INTEGER(intk), INTENT(in) :: kk, jj, ii, ncells
        TYPE(scalar_t), INTENT(in) :: sca
        INTEGER(intk), INTENT(in), CONTIGUOUS :: xpoli(:)
        REAL(realk), INTENT(in), CONTIGUOUS :: xpolr(:)
        REAL(realk), INTENT(inout), OPTIONAL :: t(kk*jj*ii), qtt(kk*jj*ii)

        ! Local variables
        INTEGER(intk) :: cellcount, pntxpoli, pntxpolr

        IF (PRESENT(t)) THEN
            pntxpoli = 1
            pntxpolr = 1
            DO cellcount = 1, ncells
                CALL wmxpoltsca_t(kk, jj, ii, sca, pntxpoli, pntxpolr, &
                    xpoli, xpolr, t)
            END DO
        END IF

        IF (PRESENT(qtt)) THEN
            pntxpoli = 1
            pntxpolr = 1
            DO cellcount = 1, ncells
                CALL wmxpoltsca_qtt(kk, jj, ii, sca, pntxpoli, pntxpolr, &
                    xpoli, xpolr, qtt)
            END DO
        END IF
    END SUBROUTINE wmxpolsoltsca


    SUBROUTINE wmxpoltsca_t(kk, jj, ii, sca, pntxpoli, pntxpolr, xpoli, &
        xpolr, var)

        ! Subroutine arguments
        INTEGER(intk), INTENT(in) :: kk, jj, ii
        TYPE(scalar_t), INTENT(in) :: sca
        INTEGER(intk), INTENT(inout) :: pntxpoli
        INTEGER(intk), INTENT(inout) :: pntxpolr
        INTEGER(intk), INTENT(in), CONTIGUOUS :: xpoli(:)
        REAL(realk), INTENT(in), CONTIGUOUS :: xpolr(:)
        REAL(realk), INTENT(inout) :: var(kk*jj*ii)

        ! Local variables
        INTEGER(intk) :: intcell, body, bctype, stencils, n, stcell
        REAL(realk) :: bcval, coeff, val

        ! 1st integer: identifier of the interface cell
        intcell = xpoli(pntxpoli)
        pntxpoli = pntxpoli + 1

        ! 2nd integer: body ID
        body = xpoli(pntxpoli)
        pntxpoli = pntxpoli + 1

        ! Look up boundary condition at bodyid
        IF (body > 0) THEN
            bctype = sca%geometries(body)%flag
            bcval = sca%geometries(body)%value
        ELSE
            ! default: adiabatic
            ! TODO: Add body 0 to above list
            bctype = 1
            bcval = 0.0
        END IF

        ! 3rd integer: cells used for the stencil (may also be 0)
        stencils = xpoli(pntxpoli)
        pntxpoli = pntxpoli + 1

        ! 1st real: areabyvol required for the flux
        ! areabyvol = xpolr(pntxpolr)
        pntxpolr = pntxpolr + 1

        ! Computing the fixed value
        val = 0.0
        DO n = 1, stencils
            stcell = xpoli(pntxpoli)
            pntxpoli = pntxpoli + 1

            coeff = xpolr(pntxpolr)
            pntxpolr = pntxpolr + 1

            val = val + var(stcell)*coeff
        END DO

        coeff = xpolr(pntxpolr)   ! additive constant "acoeffvel"
        pntxpolr = pntxpolr + 1

        ! Set value at cell if BC is fixed value
        IF (bctype == 0) THEN
            val = val + coeff*bcval
            var(intcell) = val
        END IF
    END SUBROUTINE wmxpoltsca_t


    SUBROUTINE wmxpoltsca_qtt(kk, jj, ii, sca, pntxpoli, pntxpolr, xpoli, &
            xpolr, var)

        ! Subroutine arguments
        INTEGER(intk), INTENT(in) :: kk, jj, ii
        TYPE(scalar_t), INTENT(in) :: sca
        INTEGER(intk), INTENT(inout) :: pntxpoli
        INTEGER(intk), INTENT(inout) :: pntxpolr
        INTEGER(intk), INTENT(in), CONTIGUOUS :: xpoli(:)
        REAL(realk), INTENT(in), CONTIGUOUS :: xpolr(:)
        REAL(realk), INTENT(inout) :: var(kk*jj*ii)

        ! Local variables
        INTEGER(intk) :: intcell, body, bctype, stencils
        REAL(realk) :: areabyvol, bcval

        ! 1st integer: identifier of the interface cell
        intcell = xpoli(pntxpoli)
        pntxpoli = pntxpoli + 1

        ! 2nd integer: body ID
        body = xpoli(pntxpoli)
        pntxpoli = pntxpoli + 1

        ! Look up boundary condition at bodyid
        IF (body > 0) THEN
            bctype = sca%geometries(body)%flag
            bcval = sca%geometries(body)%value
        ELSE
            ! default: adiabatic
            ! TODO: Add body 0 to above list
            bctype = 1
            bcval = 0.0
        END IF

        ! 3rd integer: number of cells used for the stencil (may also be 0)
        stencils = xpoli(pntxpoli)
        pntxpoli = pntxpoli + 1

        ! 1st real: areabyvol required for the flux
        areabyvol = xpolr(pntxpolr)
        pntxpolr = pntxpolr + 1

        ! Increment counter for unused variables
        pntxpoli = pntxpoli + stencils     ! 'stcell'
        pntxpolr = pntxpolr + stencils     ! 'coeff'
        pntxpolr = pntxpolr + 1            ! Additive constant "acoeffvel"

        ! Add flux to list of fluxes
        IF (bctype == 1) THEN
            var(intcell) = var(intcell) + areabyvol*bcval
        END IF
    END SUBROUTINE wmxpoltsca_qtt


    SUBROUTINE writestencils()
        USE MPI_f08, ONLY: MPI_Barrier, MPI_COMM_WORLD

        ! Subroutine arguments
        ! none...

        ! Local variables
        INTEGER(intk) :: i, igrid

        IF (myid == 0) THEN
            CALL create_directory("STENCILS")
        END IF
        CALL MPI_Barrier(MPI_COMM_WORLD)

        DO i = 1, nmygrids
            igrid = mygrids(i)
            CALL writestencils_grid(igrid)
        END DO
        CALL MPI_Barrier(MPI_COMM_WORLD)

        IF (myid == 0) THEN
            CALL pvtk_directory("STENCILS", "stecil-sca")
        END IF
    END SUBROUTINE writestencils


    ! TODO: Merge with flow stencil writer. If we add the bodyid to the
    ! flow stencils they should be compatible.
    SUBROUTINE writestencils_grid(igrid)
        ! Subroutine arguments
        INTEGER(intk), INTENT(in) :: igrid

        ! Local variables
        INTEGER(intk) :: kk, jj, ii, imygrid, ncells
        INTEGER(intk), POINTER, CONTIGUOUS :: xpoli(:)
        REAL(realk), POINTER, CONTIGUOUS :: xpolr(:)
        REAL(realk), POINTER, CONTIGUOUS :: x(:), y(:), z(:)
        REAL(realk), POINTER, CONTIGUOUS :: dx(:), dy(:), dz(:)

        CALL get_imygrid(imygrid, igrid)
        xpoli => scaxpoli(imygrid)%arr
        xpolr => scaxpolrvel(imygrid)%arr
        ncells = scanblg(imygrid)

        ! If no cells are present we return here...
        IF (ncells == 0) RETURN

        CALL get_fieldptr(x, "X", igrid)
        CALL get_fieldptr(y, "Y", igrid)
        CALL get_fieldptr(z, "Z", igrid)
        CALL get_fieldptr(dx, "DX", igrid)
        CALL get_fieldptr(dy, "DY", igrid)
        CALL get_fieldptr(dz, "DZ", igrid)
        CALL get_mgdims(kk, jj, ii, igrid)
        CALL writestencilsvtk(igrid, kk, jj, ii, ncells, xpoli, xpolr, &
            x, y, z, dz, dy, dz)
    END SUBROUTINE writestencils_grid


    SUBROUTINE writestencilsvtk(igrid, kk, jj, ii, nblgcells, xpoli, &
            xpolr, x, y, z, dx, dy, dz)

        ! Subroutine arguments
        INTEGER(intk), INTENT(in) :: igrid
        INTEGER(intk), INTENT(in) :: kk, jj, ii
        INTEGER(intk), INTENT(in) :: nblgcells
        INTEGER(intk), INTENT(in), CONTIGUOUS :: xpoli(:)
        REAL(realk), INTENT(in), CONTIGUOUS :: xpolr(:)
        REAL(realk), INTENT(in), CONTIGUOUS :: x(:), y(:), z(:)
        REAL(realk), INTENT(in), CONTIGUOUS :: dx(:), dy(:), dz(:)

        ! Local variables
        INTEGER(intk) :: cellcount, intcell, stencils, nstencils, n, stcell, &
            pntxpoli, pntxpolr, icell, istag(3)
        INTEGER(intk) :: unit
        INTEGER(intk) :: k, j, i
        CHARACTER(len=mglet_filename_max) :: filename
        CHARACTER(len=10) :: prefix

        istag = 0
        prefix = "stecil-sca"

        WRITE(filename, '("STENCILS/", A, "-igrid-", I0, ".vtk")') &
                TRIM(prefix), igrid

        OPEN(newunit=unit, file=TRIM(filename))
        WRITE(unit, '("# vtk DataFile Version 3.0")')
        WRITE(unit, '("MGLET SUBROUTINE writestencilsvtk")')
        WRITE(unit, '("ASCII")')
        WRITE(unit, '("DATASET UNSTRUCTURED_GRID")')

        ! Count nstencils
        nstencils = 0
        pntxpoli = 1
        DO cellcount = 1, nblgcells
            ! not needing this one
            ! intcell = xpoli(pntxpoli)
            pntxpoli = pntxpoli + 1

            ! not needing this one either
            ! body = xpoli(pntxpoli)
            pntxpoli = pntxpoli + 1

            stencils = xpoli(pntxpoli)
            pntxpoli = pntxpoli + 1

            nstencils = nstencils + stencils
            pntxpoli = pntxpoli + stencils
        END DO

        WRITE(unit, '("POINTS ", I0, " float")') nblgcells + nstencils

        ! Write points
        pntxpoli = 1
        DO cellcount = 1, nblgcells
            ! intcell the same for all components
            intcell = xpoli(pntxpoli)
            pntxpoli = pntxpoli + 1

            ! body = xpoli(pntxpoli)
            pntxpoli = pntxpoli + 1

            CALL ind2sub(intcell, k, j, i, kk, jj, ii)
            WRITE(unit, '(G0, 1X, G0, 1X, G0)') &
                x(i) + istag(1)*dx(i)/2.0, &
                y(j) + istag(2)*dy(j)/2.0, &
                z(k) + istag(3)*dz(k)/2.0

            stencils = xpoli(pntxpoli)
            pntxpoli = pntxpoli + 1

            DO n = 1, stencils
                stcell = xpoli(pntxpoli)
                pntxpoli = pntxpoli + 1

                CALL ind2sub(stcell, k, j, i, kk, jj, ii)
                WRITE(unit, '(G0, 1X, G0, 1X, G0)') &
                    x(i) + istag(1)*dx(i)/2.0, &
                    y(j) + istag(2)*dy(j)/2.0, &
                    z(k) + istag(3)*dz(k)/2.0
            END DO
        END DO

        ! Write cells
        WRITE(unit,'("CELLS ", I0, 1X, I0)') nblgcells + nstencils*2, &
            (nblgcells + nstencils)*2 + nstencils*3

        ! Points
        DO n = 1, nblgcells + nstencils
            WRITE(unit, '(I0, 1X, I0)') 1, n-1
        END DO

        ! Lines
        icell = 0
        pntxpoli = 1
        DO cellcount = 1, nblgcells
            ! not needing this one
            !intcell = xpoli(pntxpoli)
            pntxpoli = pntxpoli + 1

            ! body = xpoli(pntxpoli)
            pntxpoli = pntxpoli + 1

            stencils = xpoli(pntxpoli)
            pntxpoli = pntxpoli + 1

            DO n = 1, stencils
                WRITE(unit, '(I0, 1X, I0, 1X, I0)') 2, icell, icell+n
                pntxpoli = pntxpoli + 1
            END DO
            icell = icell + 1 + stencils
        END DO

        ! Cell types
        WRITE(unit, '("CELL_TYPES ", I0)') nblgcells + nstencils*2
        DO n = 1, nblgcells+nstencils
            WRITE(unit, '("1")')
        END DO
        DO n = 1, nstencils
            WRITE(unit, '("3")')
        END DO

        ! Write cell data stencils
        WRITE(unit,'("CELL_DATA ", I0)') nblgcells+nstencils*2
        WRITE(unit,'("SCALARS stencils int 1")')
        WRITE(unit,'("LOOKUP_TABLE default")')

        pntxpoli = 1
        DO cellcount = 1, nblgcells
            ! not needing this one
            !intcell = xpoli(pntxpoli)
            pntxpoli = pntxpoli + 1

            ! body = xpoli(pntxpoli)
            pntxpoli = pntxpoli + 1

            stencils = xpoli(pntxpoli)
            WRITE(unit, '(I0)') stencils
            pntxpoli = pntxpoli + 1

            DO n = 1, stencils
                WRITE(unit, '("-1")')
                pntxpoli = pntxpoli + 1
            END DO
        END DO

        pntxpoli = 1
        DO cellcount = 1,nblgcells
            ! not needing this one
            !intcell = xpoli(pntxpoli)
            pntxpoli = pntxpoli + 1

            ! body = xpoli(pntxpoli)
            pntxpoli = pntxpoli + 1

            stencils = xpoli(pntxpoli)
            pntxpoli = pntxpoli + 1

            DO n = 1, stencils
                WRITE(unit, '(I0)') stencils
                pntxpoli = pntxpoli + 1
            END DO
        END DO

        ! Write cell data coefficient
        WRITE(unit,'("SCALARS coefficient float 1")')
        WRITE(unit,'("LOOKUP_TABLE default")')

        pntxpoli = 1
        pntxpolr = 1
        DO cellcount = 1, nblgcells
            ! not needing this one
            !intcell = xpoli(pntxpoli)
            pntxpoli = pntxpoli + 1

            ! body = xpoli(pntxpoli)
            pntxpoli = pntxpoli + 1

            stencils = xpoli(pntxpoli)
            pntxpoli = pntxpoli + 1

            ! >>> ghost cell
            WRITE(unit, '(G0)') xpolr(pntxpolr+stencils)

            DO n = 1, stencils
                !stcell = xpoli(pntxpoli)
                pntxpoli = pntxpoli + 1

                !coeff = xpolr(pntxpolr)
                ! stencil cell
                WRITE(unit, '(G0)') xpolr(pntxpolr)
                pntxpolr = pntxpolr + 1
            END DO

            !coeff = xpolr(pntxpolr)
            pntxpolr = pntxpolr + 1
        END DO

        pntxpoli = 1
        pntxpolr = 1
        DO cellcount = 1, nblgcells
            ! not needing this one
            !intcell = xpoli(pntxpoli)
            pntxpoli = pntxpoli + 1

            ! body = xpoli(pntxpoli)
            pntxpoli = pntxpoli + 1

            stencils = xpoli(pntxpoli)
            pntxpoli = pntxpoli + 1

            DO n = 1, stencils
                !stcell = xpoli(pntxpoli)
                pntxpoli = pntxpoli + 1

                !coeff = xpolr(pntxpolr)
                ! >>> stencil cell
                WRITE(unit, '(G0)') xpolr(pntxpolr)
                pntxpolr = pntxpolr + 1
            END DO

            !coeff = xpolr(pntxpolr)
            pntxpolr = pntxpolr + 1
        END DO

        ! Phew...
        CLOSE(unit)
    END SUBROUTINE writestencilsvtk
END MODULE gc_scastencils_mod
