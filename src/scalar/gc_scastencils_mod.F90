MODULE gc_scastencils_mod
    USE core_mod
    USE ib_mod, ONLY: gc_t, wmindexlistn, findinterface2, &
        wmcheckneighbor, choosestencil
    USE scacore_mod

    IMPLICIT NONE (type, external)
    PRIVATE

    INTEGER(intk), ALLOCATABLE, TARGET :: scaxpoli(:, :)
    REAL(realk), ALLOCATABLE, TARGET :: scaxpolr(:, :), scaxpolrvel(:, :)
    INTEGER(intk), ALLOCATABLE :: scanblg(:)

    INTEGER(intk), PARAMETER :: nperi = 6 ! = up to 6 peripheral stencil points
    INTEGER(intk), PARAMETER :: isize = 3 + nperi
    INTEGER(intk), PARAMETER :: rsize = 2 + nperi

    PUBLIC :: create_scastencils, finish_scastencils, set_scastencils

CONTAINS
    SUBROUTINE create_scastencils(gc)
        ! Subroutine arguments
        TYPE(gc_t), INTENT(in) :: gc

        ! Local variables
        INTEGER(intk) :: nstencils_estimate
        INTEGER(intk) :: i, igrid, kk, jj, ii, ip3, ipp, ncells
        INTEGER(intk) :: bconds(6)
        REAL(realk), POINTER, CONTIGUOUS :: x(:), y(:), z(:)
        REAL(realk), POINTER, CONTIGUOUS :: xstag(:), ystag(:), zstag(:)
        REAL(realk), POINTER, CONTIGUOUS :: ddx(:), ddy(:), ddz(:)
        REAL(realk), POINTER, CONTIGUOUS :: bp(:, :, :)

        ! Allocate space for stencils with possible overhead
        CALL estimate_nstencils(nstencils_estimate, gc%icells)
        ALLOCATE(scaxpoli(nstencils_estimate*isize, nmygrids))
        ALLOCATE(scaxpolr(nstencils_estimate*rsize, nmygrids))
        ALLOCATE(scaxpolrvel(nstencils_estimate*rsize, nmygrids))
        ALLOCATE(scanblg(nmygrids))

        ! Stencil structure (standard-sized stencil):

        ! I1 = cell to be treated
        ! I2 = body id of cell to be treated
        ! I3 = number of peripheral stencil points (N)
        ! I4 = 1D-index of peripheral stencil point 1
        ! I5 = 1D-index of peripheral stencil point 2 (or -1 if N<2)
        ! ...
        ! I(N+3) = 1D-index of peripheral stencil point N
        ! -> therefore: isize = 3 + nperi

        ! R1 = surface area divided by cell volume
        ! R2 = coefficient for peripheral stencil point 1 (or dummy if N<1)
        ! R3 = coefficient for peripheral stencil point 2 (or dummy if N<2)
        ! ...
        ! R(N+1) = coefficient for peripheral stencil point N (or dummy)
        ! R(N+2) = additive constant (e.g. for Dirichlet conditions)
        ! -> therefore: rsize = 2 + nperi

        DO i = 1, nmygrids
            igrid = mygrids(i)
            ipp = gc%icellspointer(igrid)
            ncells = gc%icells(igrid)

            CALL get_fieldptr(x, "X", igrid)
            CALL get_fieldptr(y, "Y", igrid)
            CALL get_fieldptr(z, "Z", igrid)

            CALL get_fieldptr(xstag, "XSTAG", igrid)
            CALL get_fieldptr(ystag, "YSTAG", igrid)
            CALL get_fieldptr(zstag, "ZSTAG", igrid)

            CALL get_fieldptr(ddx, "DDX", igrid)
            CALL get_fieldptr(ddy, "DDY", igrid)
            CALL get_fieldptr(ddz, "DDZ", igrid)

            CALL get_fieldptr(bp, "BP", igrid)

            CALL get_mgdims(kk, jj, ii, igrid)
            CALL get_ip3(ip3, igrid)

            CALL get_mgbasb(bconds, igrid)

            CALL tscastencil(igrid, kk, jj, ii, x, y, z, xstag, ystag, zstag, &
                ddx, ddy, ddz, bp, gc%bzelltyp(ip3), bconds, ncells, &
                gc%nvecs(:, ipp:ipp+ncells-1), gc%arealist(ipp:ipp+ncells-1), &
                gc%bodyid(ipp:ipp+ncells-1))
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


    SUBROUTINE estimate_nstencils(nstencils, icells)
        ! Local variables
        INTEGER(intk), INTENT(out) :: nstencils
        INTEGER(intk), CONTIGUOUS, INTENT(in) :: icells(:)

        ! Local variables
        INTEGER(intk) :: i, igrid, kk, jj, ii
        INTEGER(intk) :: istencils

        nstencils = 0
        DO i = 1, nmygrids
            igrid = mygrids(i)

            CALL get_mgdims(kk, jj, ii, igrid)

            istencils = estimate_istencils(icells(igrid), kk, jj, ii)
            nstencils = MAX(nstencils, istencils)
        END DO
    END SUBROUTINE estimate_nstencils


    PURE FUNCTION estimate_istencils(icells, kk, jj, ii) RESULT(istencils)
        ! Subroutine arguments
        INTEGER(intk), INTENT(in) :: icells
        INTEGER(intk), INTENT(in) :: kk, jj, ii
        INTEGER(intk) :: istencils

        istencils = NINT(icells*1.1 + MAX(kk, jj, ii)**2)
    END FUNCTION estimate_istencils


    SUBROUTINE tscastencil(igrid, kk, jj, ii, x, y, z, xstag, ystag, zstag, &
            ddx, ddy, ddz, bp, bzelltyp, bconds, icells, nvecs, arealist, &
            bodyid)

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
        INTEGER(intk) :: pntxpoli, pntxpolr, counter
        INTEGER(intk) :: found, found2, found3, foundnr
        INTEGER(intk) :: foundx1, foundx2, foundy1, foundy2, foundz1, foundz2
        INTEGER(intk) :: xpolisize, xpolrsize, imygrid
        INTEGER(intk) :: inlst(lsize), jnlst(lsize), knlst(lsize), nnlst
        INTEGER(intk) :: body
        REAL(realk) :: area, areabyvol
        INTEGER(intk), ALLOCATABLE :: xpoli(:)
        REAL(realk), ALLOCATABLE :: xpolr(:), xpolrvel(:)
        INTEGER(intk) :: initial_pntxpoli, initial_pntxpolr

        ! Allocating the arrays with possible overhead
        xpolisize = estimate_istencils(icells, kk, jj, ii)*isize
        xpolrsize = estimate_istencils(icells, kk, jj, ii)*rsize
        ALLOCATE(xpoli(xpolisize))
        ALLOCATE(xpolr(xpolrsize))
        ALLOCATE(xpolrvel(xpolrsize))

        ! Initialization
        counter = 0
        pntxpoli = 0
        xpoli = 0
        pntxpolr = 0
        xpolr = 0.0

        ! Index stencil for looking for neighbour cells.
        CALL wmindexlistn(inlst, jnlst, knlst, nnlst)

        DO i = 3, ii-2
            DO j = 3, jj-2
                DO k = 3, kk-2
                    ! Storing the initial values of the index counters
                    initial_pntxpoli = pntxpoli
                    initial_pntxpolr = pntxpolr

                    ! Starting the stencil generation
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
                            body = 0

                            ! In case of multiple eligible cells, the biggest
                            ! bodyid is selected.
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
                                IF (bodyid(-bzelltyp(k, j-1, i)) > body) THEN
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

                    ! Cycle if no new stencil was created at all
                    IF (pntxpoli - initial_pntxpoli == 0 .AND. &
                        pntxpolr - initial_pntxpolr == 0) THEN
                            CYCLE
                    END IF

                    ! Ensuring that standard-sized stencils were created
                    IF (pntxpoli - initial_pntxpoli /= isize) THEN
                        WRITE(*, *) "Error: No standard-sized stencil ", &
                            pntxpoli, " initial: ", initial_pntxpoli
                        CALL errr(__FILE__, __LINE__)
                    END IF
                    IF (pntxpolr - initial_pntxpolr /= rsize) THEN
                        WRITE(*, *) "Error: No standard-sized stencil ", &
                            pntxpolr, " initial: ", initial_pntxpolr
                        CALL errr(__FILE__, __LINE__)
                    END IF

                END DO
            END DO
        END DO

        CALL get_imygrid(imygrid, igrid)
        scaxpoli(1:pntxpoli, imygrid) = xpoli(1:pntxpoli)
        scaxpolr(1:pntxpolr, imygrid) = xpolr(1:pntxpolr)
        scaxpolrvel(1:pntxpolr, imygrid) = xpolrvel(1:pntxpolr)
        scanblg(imygrid) = counter

        ! Final size check
        IF (counter > 0) THEN
            IF (pntxpoli /= isize * counter) THEN
                WRITE(*, *) "Error: pntxpoli / counter = ", pntxpoli, counter
                CALL errr(__FILE__, __LINE__)
            END IF
            IF (pntxpolr /= rsize * counter) THEN
                WRITE(*, *) "Error: pntxpolr / counter = ", pntxpolr, counter
                CALL errr(__FILE__, __LINE__)
            END IF
        ELSE
            IF (pntxpoli /= 0 .OR. pntxpolr /= 0) THEN
                CALL errr(__FILE__, __LINE__)
            END IF
        END IF

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
        INTEGER(intk) :: xpolisize, xpolrsize, n

        xpolisize = SIZE(xpoli)
        xpolrsize = SIZE(xpolr)

        ! Check if there is enough space in the arrays to store the new stencil
        IF (pntxpoli+isize > xpolisize) THEN
            CALL errr(__FILE__, __LINE__)
        END IF
        IF (pntxpolr+rsize > xpolrsize) THEN
            CALL errr(__FILE__, __LINE__)
        END IF

        ! 1st integer: linearized index of cell to be interpolated
        pntxpoli = pntxpoli + 1
        CALL sub2ind(xpoli(pntxpoli), k, j, i, kk, jj, ii)

        ! 2nd integer: body id of this stencil
        pntxpoli = pntxpoli + 1
        xpoli(pntxpoli) = body

        ! 3rd integer: number of the peripheral stencil points (N)
        pntxpoli = pntxpoli + 1
        xpoli(pntxpoli) = velpts  ! = 1 (formally)

        ! 1st real: setting area/volume of the interface within the cell
        ! (additional in comparison to the velocity stencil)
        pntxpolr = pntxpolr + 1
        xpolr(pntxpolr) = areabyvol
        xpolrvel(pntxpolr) = areabyvol

        ! here with iteration over all N peripheral stencil points
        ! (velpts here carries the value of "chosenvelpts")
        DO n = 1, nperi
            ! Keep incrementing the index counter anyway
            pntxpoli = pntxpoli + 1
            pntxpolr = pntxpolr + 1
            IF (n <= velpts) THEN
                IF (n /= 1) CALL errr(__FILE__, __LINE__)
                CALL sub2ind(xpoli(pntxpoli), k, j, i, kk, jj, ii)
                xpolr(pntxpolr) = 1.0
                xpolrvel(pntxpolr) = 1.0
            ELSE
                ! Setting unused dummy values (yields standard-sized stencils)
                xpoli(pntxpoli) = -1
                xpolr(pntxpolr) = -HUGE(0.0_realk)
                xpolrvel(pntxpolr) = -HUGE(0.0_realk)
            END IF
        END DO

        ! 3rd real: setting the additive constant
        pntxpolr = pntxpolr + 1
        xpolr(pntxpolr) = 0.0
        xpolrvel(pntxpolr) =  0.0

        ! >>> Integer counter += isize
        ! >>> Real counter += rsize

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

        USE, INTRINSIC :: IEEE_ARITHMETIC, ONLY: IEEE_IS_FINITE

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

        ! --- Cell value stencil --- (comparable to flux stencil)
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
            acoeffstc = (acoeffstc * dvol / (cvol1 * dvol))
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

        ! test coefficients for invalid floating-point values
        nanstenc = 0
        DO n = 1, velpts
            IF (.NOT. IEEE_IS_FINITE(coeffstctsca(n))) THEN
                nanstenc = 1
            END IF
            IF (.NOT. IEEE_IS_FINITE(coeffstc(n))) THEN
                nanstenc = 1
            END IF
        END DO
        IF (.NOT. IEEE_IS_FINITE(acoeffstctsca)) THEN
            nanstenc = 1
        END IF
        IF (.NOT. IEEE_IS_FINITE(acoeffstc)) THEN
            nanstenc = 1
        END IF

        IF (nanstenc == 1) THEN
            WRITE(*, *) "-----------------------------------------------"
            WRITE(*, *) "i, j, k", i, j, k
            DO n = 1, velpts
                WRITE(*, *) "dist(n)", n, dist(n)
            END DO
            WRITE(*, *) "dall", dall
            DO n = 1, velpts
                write(*, *) "ce(n)", n, ce(n)
            END DO
            WRITE(*, *) "cesum", cesum
            WRITE(*, *) "velpts", velpts
            WRITE(*, *) "a, b, c, dpl", a, b, c, dpl
            WRITE(*, *) "dpoint", dpoint
            WRITE(*, *)
            DO n = 1, velpts
                write(*, *) "coeffstc(n)", n, coeffstc(n)
            END DO
            WRITE(*, *) "acoeffstc", acoeffstc
            WRITE(*, *)
            DO n = 1, velpts
                write(*, *) "coeffstctsca(n)", n, coeffstctsca(n)
            END DO
            WRITE(*, *) "acoeffstctsca", acoeffstctsca
            WRITE(*, *) "sr calcfluxtsca: acoeffstc", acoeffstc
            WRITE(*, *) "cvol1, cvol2", cvol1, cvol2
            WRITE(*, *) "-----------------------------------------------"

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

        ! Check if the arrays are large enough to hold the new stencil data
        IF (pntxpoli+isize > xpolisize) THEN
            CALL errr(__FILE__, __LINE__)
        END IF
        IF (pntxpolr+rsize > xpolrsize) THEN
            CALL errr(__FILE__, __LINE__)
        END IF

        ! 1st integer: linearized index of cell to be treated
        pntxpoli = pntxpoli + 1
        CALL sub2ind(xpoli(pntxpoli), k, j, i, kk, jj, ii)

        ! 2nd integer: body id of this stencil
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

        ! Check that the number of peripheral stencil points <= nperi
        IF (velpts > nperi) THEN
            WRITE(*, *) "Number of peripheral stencil points > ", nperi
            CALL errr(__FILE__, __LINE__)
        END IF

        ! here with iteration over all N peripheral stencil points
        ! (velpts here carries the value of "chosenvelpts")
        DO n = 1, nperi
            ! Keep incrementing the index counter anyway
            pntxpoli = pntxpoli + 1
            pntxpolr = pntxpolr + 1
            IF (n <= velpts) THEN
                ! Setting the used values
                xpoli(pntxpoli) = kstc(n) + (jstc(n)-1)*kk + (istc(n)-1)*jj*kk
                xpolr(pntxpolr) = coeffstc(n)
                xpolrvel(pntxpolr) = coeffstcvel(n)
            ELSE
                ! Setting unused dummy values (yields standard-sized stencils)
                xpoli(pntxpoli) = -1
                xpolr(pntxpolr) = -HUGE(0.0_realk)
                xpolrvel(pntxpolr) = -HUGE(0.0_realk)
            END IF
        END DO

        ! last real: setting the additive constant
        pntxpolr = pntxpolr + 1
        xpolr(pntxpolr) = acoeffstc
        xpolrvel(pntxpolr) = acoeffstcvel

        ! >>> Integer counter += isize
        ! >>> Real counter += rsize

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
        INTEGER(intk) :: xpolisize, xpolrsize, n

        xpolisize = SIZE(xpoli)
        xpolrsize = SIZE(xpolr)

        ! Checking if the arrays are large enough to hold the new stencil data
        IF (pntxpoli+isize > xpolisize) THEN
            CALL errr(__FILE__, __LINE__)
        END IF
        IF (pntxpolr+rsize > xpolrsize) THEN
            CALL errr(__FILE__, __LINE__)
        END IF

        ! increasing the index counter only by one
        ! (no iteration over points as in function above)

        ! 1st integer: linearized index of cell to be treated
        pntxpoli = pntxpoli + 1
        CALL sub2ind(xpoli(pntxpoli), k, j, i, kk, jj, ii)

        ! 2nd integer: body id of this stencil
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
        DO n = 1, nperi
            ! Keep incrementing the index counter anyway
            pntxpoli = pntxpoli + 1
            pntxpolr = pntxpolr + 1
            IF (n <= velpts) THEN
                ! Setting the used values
                WRITE(*, *) "This case never applies as VELPTS = 0"
                CALL errr(__FILE__, __LINE__)
            ELSE
                ! Setting unused dummy values (yields standard-sized stencils)
                xpoli(pntxpoli) = -1
                xpolr(pntxpolr) = -HUGE(0.0_realk)
                xpolrvel(pntxpolr) = -HUGE(0.0_realk)
            END IF
        END DO

        ! 3rd real: setting the additive constant (here: equal to the final cell value)
        pntxpolr = pntxpolr + 1
        xpolr(pntxpolr) = 1.0
        xpolrvel(pntxpolr) =  1.0

        ! >>> Integer counter += isize
        ! >>> Real counter += rsize

    END SUBROUTINE tscastencilcoeffrescue


    SUBROUTINE set_scastencils(ctyp, sca, t, qtt)
        ! Subroutine arguments
        CHARACTER(len=1), INTENT(in) :: ctyp
        TYPE(scalar_t), INTENT(in) :: sca
        TYPE(field_t), INTENT(inout), OPTIONAL :: t, qtt

        ! Local variables
        REAL(realk), CONTIGUOUS, POINTER :: xpolr(:, :)

        CALL start_timer(412)

        IF (.NOT. PRESENT(t) .AND. .NOT. PRESENT(qtt)) THEN
            CALL errr(__FILE__, __LINE__)
        END IF

        SELECT CASE(ctyp)
        CASE("C")
            ! Warning: Setting stencils for scalar cell values is no longer
            ! supported and tested.
            xpolr => scaxpolr
        CASE("P")
            xpolr => scaxpolrvel
        CASE DEFAULT
            CALL errr(__FILE__, __LINE__)
        END SELECT

        IF (PRESENT(t)) THEN
            CALL set_scastencils_t(xpolr, sca, t)
        END IF

        IF (PRESENT(qtt)) THEN
            CALL set_scastencils_qtt(xpolr, sca, qtt)
        END IF

        CALL stop_timer(412)
    END SUBROUTINE set_scastencils


    SUBROUTINE set_scastencils_t(xpolr, sca, t_f)
        ! Subroutine arguments
        REAL(realk), CONTIGUOUS, INTENT(in) :: xpolr(:, :)
        TYPE(scalar_t), INTENT(in) :: sca
        TYPE(field_t), INTENT(inout) :: t_f

        ! Local variables
        INTEGER(intk) :: i, igrid, kk, jj, ii
        REAL(realk), POINTER, CONTIGUOUS :: t(:, :, :)

        DO i = 1, nmygrids
            igrid = mygrids(i)
            CALL get_mgdims(kk, jj, ii, igrid)

            CALL get_grid3_real(t, t_f, igrid)

            ! Pass stencil i as index to avoid slicing (:, i), which the
            ! compiler does not accept with a CONTIGUOUS shape
            CALL wmxpoltsca_t(kk, jj, ii, scanblg(i), sca, scaxpoli, xpolr, i, &
                t)
        END DO
    END SUBROUTINE set_scastencils_t


    SUBROUTINE wmxpoltsca_t(kk, jj, ii, ncells, sca, xpoli, xpolr, istencil, &
        var)
        ! Subroutine arguments
        INTEGER(intk), INTENT(in) :: kk, jj, ii
        INTEGER(intk), INTENT(in) :: ncells
        TYPE(scalar_t), INTENT(in) :: sca
        INTEGER(intk), CONTIGUOUS, INTENT(in) :: xpoli(:, :)
        REAL(realk), CONTIGUOUS, INTENT(in) :: xpolr(:, :)
        INTEGER(intk), INTENT(in) :: istencil
        REAL(realk), INTENT(inout) :: var(kk*jj*ii)

        ! Local variables
        INTEGER(intk) :: cellcount, start_pntxpoli, start_pntxpolr
        INTEGER(intk) :: intcell, body, bctype, stencils, n, stcell
        REAL(realk) :: bcval, coeff, val

        DO cellcount = 1, ncells
            start_pntxpoli = (cellcount-1)*isize
            start_pntxpolr = (cellcount-1)*rsize

            ! 1st integer: identifier of the interface cell
            intcell = xpoli(start_pntxpoli + 1, istencil)

            ! 2nd integer: body ID
            body = xpoli(start_pntxpoli + 2, istencil)

            ! Look up boundary condition at bodyid
            bctype = 1
            bcval = 0.0
            IF (body > 0) THEN
                bctype = sca%geometries(body)%flag
                bcval = sca%geometries(body)%value
            END IF

            ! 3rd integer: cells used for the stencil (may also be 0)
            stencils = xpoli(start_pntxpoli + 3, istencil)

            ! Computing the fixed value
            val = 0.0
            DO n = 1, stencils
                stcell = xpoli(start_pntxpoli + 3 + n, istencil)
                coeff = xpolr(start_pntxpolr + 1 + n, istencil)
                val = val + var(stcell)*coeff
            END DO

            ! additive constant "acoeffvel"
            coeff = xpolr(start_pntxpolr + (2 + nperi), istencil)

            ! Set value at cell if BC is fixed value
            IF (bctype == 0) THEN
                val = val + coeff*bcval
                var(intcell) = val
            END IF
        END DO
    END SUBROUTINE wmxpoltsca_t


    SUBROUTINE set_scastencils_qtt(xpolr, sca, qtt_f)
        ! Subroutine arguments
        REAL(realk), CONTIGUOUS, INTENT(in) :: xpolr(:, :)
        TYPE(scalar_t), INTENT(in) :: sca
        TYPE(field_t), INTENT(inout) :: qtt_f

        ! Local variables
        INTEGER(intk) :: i, igrid, kk, jj, ii
        REAL(realk), POINTER, CONTIGUOUS :: qtt(:, :, :)

        DO i = 1, nmygrids
            igrid = mygrids(i)
            CALL get_mgdims(kk, jj, ii, igrid)

            CALL get_grid3_real(qtt, qtt_f, igrid)

            ! Pass stencil i as index to avoid slicing (:, i), which the
            ! compiler does not accept with a CONTIGUOUS shape
            CALL wmxpoltsca_qtt(kk, jj, ii, scanblg(i), sca, scaxpoli, xpolr, &
                i, qtt)
        END DO
    END SUBROUTINE set_scastencils_qtt


    SUBROUTINE wmxpoltsca_qtt(kk, jj, ii, ncells, sca, xpoli, xpolr, istencil, &
        var)
        ! Subroutine arguments
        INTEGER(intk), INTENT(in) :: kk, jj, ii
        INTEGER(intk), INTENT(in) :: ncells
        TYPE(scalar_t), INTENT(in) :: sca
        INTEGER(intk), CONTIGUOUS, INTENT(in) :: xpoli(:, :)
        REAL(realk), CONTIGUOUS, INTENT(in) :: xpolr(:, :)
        INTEGER(intk), INTENT(in) :: istencil
        REAL(realk), INTENT(inout) :: var(kk*jj*ii)

        ! Local variables
        INTEGER(intk) :: cellcount, start_pntxpoli, start_pntxpolr
        INTEGER(intk) :: intcell, body, bctype
        REAL(realk) :: areabyvol, bcval

        DO cellcount = 1, ncells
            start_pntxpoli = (cellcount-1)*isize
            start_pntxpolr = (cellcount-1)*rsize

            ! 1st integer: identifier of the interface cell
            intcell = xpoli(start_pntxpoli + 1, istencil)

            ! 2nd integer: body ID
            body = xpoli(start_pntxpoli + 2, istencil)

            ! Look up boundary condition at bodyid
            bctype = 1
            bcval = 0.0
            IF (body > 0) THEN
                bctype = sca%geometries(body)%flag
                bcval = sca%geometries(body)%value
            END IF

            ! 1st real: areabyvol required for the flux
            areabyvol = xpolr(start_pntxpolr + 1, istencil)

            ! Add flux to list of fluxes
            IF (bctype == 1) THEN
                var(intcell) = var(intcell) + areabyvol*bcval
            END IF
        END DO
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
        xpoli => scaxpoli(:, imygrid)
        xpolr => scaxpolrvel(:, imygrid)
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
            x, y, z, dx, dy, dz)
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
            icell, istag(3)
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

        ! I1 = cell to be treated
        ! I2 = body id of cell to be treated
        ! I3 = number of peripheral stencil points (N = nperi)
        ! I4 = 1D-index of peripheral stencil point 1 (or -1 if N<1)
        ! I5 = 1D-index of peripheral stencil point 2 (or -1 if N<2)
        ! ...
        ! I(N+3) = 1D-index of peripheral stencil point N (or -1 if N<N)
        ! -> therefore: isize = 3 + nperi

        ! R1 = surface area divided by cell volume
        ! R2 = coefficient for peripheral stencil point 1 (or dummy if N<1)
        ! R3 = coefficient for peripheral stencil point 2 (or dummy if N<2)
        ! ...
        ! R(N+1) = coefficient for peripheral stencil point N (or dummy if N<N)
        ! R(N+2) = additive constant (e.g. for Dirichlet conditions)
        ! -> therefore: rsize = 2 + nperi


        ! Count nstencils
        nstencils = 0
        DO cellcount = 0, nblgcells-1
            stencils = xpoli(isize*cellcount+3)
            nstencils = nstencils + stencils
        END DO

        WRITE(unit, '("POINTS ", I0, " float")') nblgcells + nstencils

        ! Write points
        DO cellcount = 0, nblgcells-1
            intcell = xpoli(isize*cellcount+1)
            CALL ind2sub(intcell, k, j, i, kk, jj, ii)
            WRITE(unit, '(G0, 1X, G0, 1X, G0)') &
                x(i) + istag(1)*dx(i)/2.0, &
                y(j) + istag(2)*dy(j)/2.0, &
                z(k) + istag(3)*dz(k)/2.0

            stencils = xpoli(isize*cellcount+3)
            DO n = 1, stencils
                stcell = xpoli(isize*cellcount+3+n)
                CALL ind2sub(stcell, k, j, i, kk, jj, ii)
                WRITE(unit, '(G0, 1X, G0, 1X, G0)') &
                    x(i) + istag(1)*dx(i)/2.0, &
                    y(j) + istag(2)*dy(j)/2.0, &
                    z(k) + istag(3)*dz(k)/2.0
            END DO
        END DO

        ! Write cells
        WRITE(unit, '("CELLS ", I0, 1X, I0)') nblgcells + nstencils*2, &
            (nblgcells + nstencils)*2 + nstencils*3

        ! Points
        DO n = 1, nblgcells + nstencils
            WRITE(unit, '(I0, 1X, I0)') 1, n-1
        END DO

        ! Lines
        icell = 0
        DO cellcount = 0, nblgcells-1
            stencils = xpoli(isize*cellcount+3)
            DO n = 1, stencils
                WRITE(unit, '(I0, 1X, I0, 1X, I0)') 2, icell, icell+n
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
        WRITE(unit, '("CELL_DATA ", I0)') nblgcells+nstencils*2
        WRITE(unit, '("SCALARS stencils int 1")')
        WRITE(unit, '("LOOKUP_TABLE default")')

        DO cellcount = 0, nblgcells-1
            stencils = xpoli(isize*cellcount+3)
            WRITE(unit, '(I0)') stencils
            DO n = 1, stencils
                WRITE(unit, '("-1")')
            END DO
        END DO

        DO cellcount = 0, nblgcells-1
            stencils = xpoli(isize*cellcount+3)
            DO n = 1, stencils
                WRITE(unit, '(I0)') stencils
            END DO
        END DO

        ! Write cell data coefficient
        WRITE(unit, '("SCALARS coefficient float 1")')
        WRITE(unit, '("LOOKUP_TABLE default")')


        DO cellcount = 0, nblgcells-1
            stencils = xpoli(isize*cellcount+3)
            WRITE(unit, '(G0)') xpolr(rsize*cellcount+2+nperi)
            DO n = 1, stencils
                WRITE(unit, '(G0)') xpolr(rsize*cellcount+1+n)
            END DO
        END DO

        DO cellcount = 0, nblgcells-1
            stencils = xpoli(isize*cellcount+3)
            DO n = 1, stencils
                WRITE(unit, '(G0)') xpolr(rsize*cellcount+1+n)
            END DO
        END DO

        ! Phew...
        CLOSE(unit)
    END SUBROUTINE writestencilsvtk
END MODULE gc_scastencils_mod
