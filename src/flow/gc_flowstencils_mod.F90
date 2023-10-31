MODULE gc_flowstencils_mod
    USE bound_flow_mod
    USE core_mod
    USE ib_mod

    USE gc_createstencils_mod, ONLY: wmcheckneighbor, choosestencil, &
        wmmultimatrix, wmaddcoefflist, wmindexlistn

    IMPLICIT NONE (type, external)
    PRIVATE

    TYPE(int_stencils_t), ALLOCATABLE, TARGET :: fxpoli(:)
    TYPE(real_stencils_t), ALLOCATABLE, TARGET :: fxpolr(:)
    INTEGER(intk), ALLOCATABLE :: fnblg(:)

    TYPE(int_stencils_t), ALLOCATABLE, TARGET :: uxpoli(:)
    TYPE(real_stencils_t), ALLOCATABLE, TARGET :: uxpolr(:), uxpolrvel(:)
    INTEGER(intk), ALLOCATABLE :: unblg(:)
    TYPE(int_stencils_t), ALLOCATABLE, TARGET :: upoldsol(:), upoldsolvel(:)
    TYPE(real_stencils_t), ALLOCATABLE, TARGET :: uoldsol(:), uoldsolvel(:)

    TYPE(int_stencils_t), ALLOCATABLE, TARGET :: vxpoli(:)
    TYPE(real_stencils_t), ALLOCATABLE, TARGET ::vxpolr(:), vxpolrvel(:)
    INTEGER(intk), ALLOCATABLE :: vnblg(:)
    TYPE(int_stencils_t), ALLOCATABLE, TARGET :: vpoldsol(:), vpoldsolvel(:)
    TYPE(real_stencils_t), ALLOCATABLE, TARGET :: voldsol(:), voldsolvel(:)

    TYPE(int_stencils_t), ALLOCATABLE, TARGET :: wxpoli(:)
    TYPE(real_stencils_t), ALLOCATABLE, TARGET :: wxpolr(:), wxpolrvel(:)
    INTEGER(intk), ALLOCATABLE :: wnblg(:)
    TYPE(int_stencils_t), ALLOCATABLE, TARGET :: wpoldsol(:), wpoldsolvel(:)
    TYPE(real_stencils_t), ALLOCATABLE, TARGET :: woldsol(:), woldsolvel(:)

    PUBLIC :: create_flowstencils, setpointvalues, setibvalues, getibvalues, &
        setsdivfield, finish_flowstencils

CONTAINS
    SUBROUTINE create_flowstencils(gc)

        ! Subroutine arguments
        TYPE(gc_t), INTENT(in) :: gc

        ! Local variables
        INTEGER(intk) :: ilevel
        REAL(realk), POINTER, CONTIGUOUS :: bp(:), bu(:), bv(:), bw(:)
        REAL(realk), POINTER, CONTIGUOUS :: au(:), av(:), aw(:)

        CALL get_fieldptr(bp, "BP")
        CALL get_fieldptr(bu, "BU")
        CALL get_fieldptr(bv, "BV")
        CALL get_fieldptr(bw, "BW")

        CALL get_fieldptr(au, "AU")
        CALL get_fieldptr(av, "AV")
        CALL get_fieldptr(aw, "AW")

        ! Always allocate - createstencils should be called only once.
        ALLOCATE(fxpoli(nmygrids))
        ALLOCATE(fxpolr(nmygrids))
        ALLOCATE(fnblg(nmygrids))

        ALLOCATE(uxpoli(nmygrids))
        ALLOCATE(uxpolr(nmygrids))
        ALLOCATE(uxpolrvel(nmygrids))
        ALLOCATE(unblg(nmygrids))
        ALLOCATE(upoldsol(nmygrids))
        ALLOCATE(uoldsol(nmygrids))
        ALLOCATE(upoldsolvel(nmygrids))
        ALLOCATE(uoldsolvel(nmygrids))

        ALLOCATE(vxpoli(nmygrids))
        ALLOCATE(vxpolr(nmygrids))
        ALLOCATE(vxpolrvel(nmygrids))
        ALLOCATE(vnblg(nmygrids))
        ALLOCATE(vpoldsol(nmygrids))
        ALLOCATE(voldsol(nmygrids))
        ALLOCATE(vpoldsolvel(nmygrids))
        ALLOCATE(voldsolvel(nmygrids))

        ALLOCATE(wxpoli(nmygrids))
        ALLOCATE(wxpolr(nmygrids))
        ALLOCATE(wxpolrvel(nmygrids))
        ALLOCATE(wnblg(nmygrids))
        ALLOCATE(wpoldsol(nmygrids))
        ALLOCATE(woldsol(nmygrids))
        ALLOCATE(wpoldsolvel(nmygrids))
        ALLOCATE(woldsolvel(nmygrids))

        DO ilevel = minlevel, maxlevel
            CALL createstencils_level(ilevel, bp, bu, bv, bw, au, av, aw, &
                gc%bzelltyp, gc%icells, gc%icellspointer, gc%bodyid, &
                gc%nvecs, gc%ucell)
        END DO

        ! mglet_dbg_envvar is in buildinfo_mod and initialized at startup
        IF (INDEX(mglet_dbg_envvar, "stencilvtk") > 0) THEN
            CALL writestencils()
        END IF
    END SUBROUTINE create_flowstencils


    SUBROUTINE finish_flowstencils()
        DEALLOCATE(fxpoli)
        DEALLOCATE(fxpolr)
        DEALLOCATE(fnblg)

        DEALLOCATE(uxpoli)
        DEALLOCATE(uxpolr)
        DEALLOCATE(uxpolrvel)
        DEALLOCATE(unblg)
        DEALLOCATE(upoldsol)
        DEALLOCATE(uoldsol)
        DEALLOCATE(upoldsolvel)
        DEALLOCATE(uoldsolvel)

        DEALLOCATE(vxpoli)
        DEALLOCATE(vxpolr)
        DEALLOCATE(vxpolrvel)
        DEALLOCATE(vnblg)
        DEALLOCATE(vpoldsol)
        DEALLOCATE(voldsol)
        DEALLOCATE(vpoldsolvel)
        DEALLOCATE(voldsolvel)

        DEALLOCATE(wxpoli)
        DEALLOCATE(wxpolr)
        DEALLOCATE(wxpolrvel)
        DEALLOCATE(wnblg)
        DEALLOCATE(wpoldsol)
        DEALLOCATE(woldsol)
        DEALLOCATE(wpoldsolvel)
        DEALLOCATE(woldsolvel)
    END SUBROUTINE finish_flowstencils


    SUBROUTINE createstencils_level(ilevel, bp, bu, bv, bw, au, av, aw, &
            bzelltyp, icells, icellspointer, bodyid, nvecs, ucell)

        ! Subroutine arguments
        INTEGER(intk), INTENT(in) :: ilevel
        REAL(realk), INTENT(in) :: bp(*), bu(*), bv(*), bw(*), &
            au(*), av(*), aw(*)
        INTEGER(intk), INTENT(in) :: bzelltyp(*)
        INTEGER(intk), CONTIGUOUS, INTENT(in) :: icells(:)
        INTEGER(intk), CONTIGUOUS, INTENT(in) :: icellspointer(:)
        INTEGER(intk), CONTIGUOUS, INTENT(in) :: bodyid(:)
        REAL(realk), CONTIGUOUS, INTENT(in) :: nvecs(:, :)
        REAL(realk), CONTIGUOUS, INTENT(in) :: ucell(:, :)

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

            CALL fluxcorrection(igrid, kk, jj, ii, ddx, ddy, ddz, &
                bp(ip3), bu(ip3), bv(ip3), bw(ip3), au(ip3), av(ip3), &
                aw(ip3), bzelltyp(ip3), icells(igrid), &
                ucell(:, ipp:ipp+ncells-1))

            CALL fluxstencil(igrid, kk, jj, ii, x, y, z, xstag, ystag, &
                zstag, bu(ip3), bzelltyp(ip3), bconds, icells(igrid), &
                nvecs(:, ipp:ipp+ncells-1), ucell(:, ipp:ipp+ncells-1), 1)
            CALL fluxstencil(igrid, kk, jj, ii, x, y, z, xstag, ystag, &
                zstag, bv(ip3), bzelltyp(ip3), bconds, icells(igrid), &
                nvecs(:, ipp:ipp+ncells-1), ucell(:, ipp:ipp+ncells-1), 2)
            CALL fluxstencil(igrid, kk, jj, ii, x, y, z, xstag, ystag, &
                zstag, bw(ip3), bzelltyp(ip3), bconds, icells(igrid), &
                nvecs(:, ipp:ipp+ncells-1), ucell(:, ipp:ipp+ncells-1), 3)
        END DO
    END SUBROUTINE createstencils_level


    SUBROUTINE fluxcorrection(igrid, kk, jj, ii, ddx, ddy, ddz, &
            bp, bu, bv, bw, au, av, aw, bzelltyp, icells, ucell)
        !----------------------------------------------------------------------
        ! SUBROUTINE FLUXCORRECTION
        !
        !    description:
        !       creates correction factors for fluxes.
        !
        !       c-Koeffzient bestimmt ob die Flaeche zur DIV berechnung
        !       benutzt wird:
        !       C = 0 -> nein C = 1 -> ja
        !       cx1  BAC Zellflaeche
        !       cx2  FRO Zellflaeche
        !       cy1  LFT
        !       cy2  RGT
        !       cz1  TOP
        !       cz2  BOT
        !
        !
        !       a = Zellfaeche
        !         = Anteil der Flaeche im Fluid, wenn Interface zu offener
        !           Zelle (veraenderbare Geschwindigkeit)
        !         = 0, wenn Interface zu geschlossener Zelle
        !
        !       ax1  BAC Zellflaeche
        !       ax2  FRO Zellflaeche
        !       ay1  LFT
        !       ay2  RGT
        !       az1  TOP
        !       az2  BOT
        !
        !       acoeffstc = Fluss durch die Flaechen komplett im Koerper
        !
        !----------------------------------------------------------------------

        USE findinterface_mod, ONLY: findinterface2
        USE ibconst_mod, ONLY: maccur

        ! Subroutine arguments
        INTEGER(intk), INTENT(in) :: igrid
        INTEGER(intk), INTENT(in) :: kk, jj, ii
        REAL(realk), INTENT(in) :: ddx(ii), ddy(jj), ddz(kk)
        REAL(realk), INTENT(in) :: bp(kk, jj, ii)
        REAL(realk), INTENT(in) :: bu(kk, jj, ii), bv(kk, jj, ii), &
            bw(kk, jj, ii)
        REAL(realk), INTENT(in) :: au(kk, jj, ii), av(kk, jj, ii), &
            aw(kk, jj, ii)
        INTEGER(intk), INTENT(in) :: bzelltyp(kk, jj, ii)
        INTEGER(intk), INTENT(in) :: icells
        REAL(realk), CONTIGUOUS, INTENT(in) :: ucell(:, :)

        ! Local variables
        INTEGER(intk) :: i, j, k
        INTEGER(intk) :: pntxpoli, pntxpolr
        INTEGER(intk) :: foundnr
        INTEGER(intk) :: foundx1, foundx2, foundy1, found
        INTEGER(intk) :: foundy2, foundz1, foundz2, counter
        INTEGER(intk) :: imax, jmax, kmax
        INTEGER(intk) :: found_dum
        INTEGER(intk) :: ngeomnbr
        INTEGER(intk) :: nfro, nbac, nrgt, nlft, nbot, ntop

        REAL(realk) :: xvel, yvel, zvel, ax1, ax2, ay1, ay2, az1, az2, refarea
        REAL(realk) :: cx1, cx2, cy1, cy2, cz1, cz2, dx, dy, dz, sarea
        REAL(realk) :: acoeffstc

        INTEGER(intk) :: xpolisize, xpolrsize, imygrid
        INTEGER(intk), ALLOCATABLE :: xpoli(:)
        REAL(realk), ALLOCATABLE :: xpolr(:)

        INTEGER(intk), PARAMETER :: additionalcells  = 10000

        CALL get_mgbasb(nfro, nbac, nrgt, nlft, nbot, ntop, igrid)

        xpolisize = NINT((icells*1.3 + additionalcells)*2)
        xpolrsize = NINT((icells*1.3 + additionalcells)*13)
        ALLOCATE(xpoli(xpolisize), xpolr(xpolrsize))

        pntxpoli = 0
        pntxpolr = 0

        counter = 0

        kmax = kk - 2
        jmax = jj - 2
        imax = ii - 2

        ! bei CONNECT-RB
        IF (nbac == 7) imax = ii - 1
        IF (nlft == 7) jmax = jj - 1
        IF (ntop == 7) kmax = kk - 1

        DO i = 3, imax
            DO j = 3, jmax
                DO k = 3, kmax

                    CALL findinterface2(k, j, i, kk, jj, ii, bp, found, &
                        foundx1, foundx2, foundy1, foundy2, foundz1, &
                        foundz2, foundnr)

                    ! If interface cell not found, skip rest of loop
                    IF (found /= 1) CYCLE

                    dx = ddx(i)
                    dy = ddy(j)
                    dz = ddz(k)

                    ! found = 1 wenn Geschw. geblockt, aber Nachbargeschw. offen
                    CALL findinterface2(k, j, i-1, kk, jj, ii , bu, found_dum)
                    cx2 = REAL(found_dum)

                    CALL findinterface2(k, j, i, kk, jj, ii, bu, found_dum)
                    cx1 = REAL(found_dum)

                    CALL findinterface2(k, j-1, i, kk, jj, ii, bv, found_dum)
                    cy2 = REAL(found_dum)

                    CALL findinterface2(k, j, i, kk, jj, ii, bv, found_dum)
                    cy1 = REAL(found_dum)

                    CALL findinterface2(k-1, j, i, kk, jj, ii, bw, found_dum)
                    cz2 = REAL(found_dum)

                    CALL findinterface2(k, j, i, kk, jj, ii, bw, found_dum)
                    cz1 = REAL(found_dum)

                    ! Koerpergeschwindigkeit
                    IF (bzelltyp(k, j, i) < 0) THEN
                        ! in der Zelle selbst ist Geometrie
                        IF (bzelltyp(k, j, i) < -icells) THEN
                            WRITE(*, *) 'fluxcorrection: -bzelltyp(k, j, i)', &
                                bzelltyp(k, j, i), k, j, i
                            CALL errr(__FILE__, __LINE__)
                        END IF
                        xvel = ucell(1, -bzelltyp(k, j, i))
                        yvel = ucell(2, -bzelltyp(k, j, i))
                        zvel = ucell(3, -bzelltyp(k, j, i))
                    ELSE
                        ngeomnbr = 0

                        IF (foundx1 == 1 .AND. bzelltyp(k, j, i+1) < 0) &
                            ngeomnbr = ngeomnbr+1
                        IF (foundx2 == 1 .AND. bzelltyp(k, j, i-1) < 0) &
                            ngeomnbr = ngeomnbr+1
                        IF (foundy1 == 1 .AND. bzelltyp(k, j+1, i) < 0) &
                            ngeomnbr = ngeomnbr+1
                        IF (foundy2 == 1 .AND. bzelltyp(k, j-1, i) < 0) &
                            ngeomnbr = ngeomnbr+1
                        IF (foundz1 == 1 .AND. bzelltyp(k+1, j, i) < 0) &
                            ngeomnbr = ngeomnbr+1
                        IF (foundz2 == 1 .AND. bzelltyp(k-1, j, i) < 0) &
                            ngeomnbr = ngeomnbr+1

                        IF (ngeomnbr == 0) THEN
                            IF (foundx1 == 0 .AND. bzelltyp(k, j, i+1) < 0) &
                                ngeomnbr = ngeomnbr-1
                            IF (foundx2 == 0 .AND. bzelltyp(k, j, i-1) < 0) &
                                ngeomnbr = ngeomnbr-1
                            IF (foundy1 == 0 .AND. bzelltyp(k, j+1, i) < 0) &
                                ngeomnbr = ngeomnbr-1
                            IF (foundy2 == 0 .AND. bzelltyp(k, j-1, i) < 0) &
                                ngeomnbr = ngeomnbr-1
                            IF (foundz1 == 0 .AND. bzelltyp(k+1, j, i) < 0) &
                                ngeomnbr = ngeomnbr-1
                            IF (foundz2 == 0 .AND. bzelltyp(k-1, j, i) < 0) &
                                ngeomnbr = ngeomnbr-1
                        END IF

                        xvel = 0.0
                        yvel = 0.0
                        zvel = 0.0
                        IF (ngeomnbr > 0) THEN
                            ! in einer offenen Nachbarzelle ist Geometrie
                            IF (foundx1 == 1 .AND. bzelltyp(k, j, i+1) < 0) THEN
                                xvel = xvel + ucell(1, -bzelltyp(k, j, i+1))
                                yvel = yvel + ucell(2, -bzelltyp(k, j, i+1))
                                zvel = zvel + ucell(3, -bzelltyp(k, j, i+1))
                            END IF
                            IF (foundx2 == 1 .AND. bzelltyp(k, j, i-1) < 0) THEN
                                xvel = xvel + ucell(1, -bzelltyp(k, j, i-1))
                                yvel = yvel + ucell(2, -bzelltyp(k, j, i-1))
                                zvel = zvel + ucell(3, -bzelltyp(k, j, i-1))
                            END IF
                            IF (foundy1 == 1 .AND. bzelltyp(k, j+1, i) < 0) THEN
                                xvel = xvel + ucell(1, -bzelltyp(k, j+1, i))
                                yvel = yvel + ucell(2, -bzelltyp(k, j+1, i))
                                zvel = zvel + ucell(3, -bzelltyp(k, j+1, i))
                            END IF
                            IF (foundy2 == 1 .AND. bzelltyp(k, j-1, i) < 0) THEN
                                xvel = xvel + ucell(1, -bzelltyp(k, j-1, i))
                                yvel = yvel + ucell(2, -bzelltyp(k, j-1, i))
                                zvel = zvel + ucell(3, -bzelltyp(k, j-1, i))
                            END IF
                            IF (foundz1 == 1 .AND. bzelltyp(k+1, j, i) < 0) THEN
                                xvel = xvel + ucell(1, -bzelltyp(k+1, j, i))
                                yvel = yvel + ucell(2, -bzelltyp(k+1, j, i))
                                zvel = zvel + ucell(3, -bzelltyp(k+1, j, i))
                            END IF
                            IF (foundz2 == 1 .AND. bzelltyp(k-1, j, i) < 0) THEN
                                xvel = xvel + ucell(1, -bzelltyp(k-1, j, i))
                                yvel = yvel + ucell(2, -bzelltyp(k-1, j, i))
                                zvel = zvel + ucell(3, -bzelltyp(k-1, j, i))
                            END IF
                            xvel = xvel/ngeomnbr
                        ELSE IF (ngeomnbr < 0) THEN
                            ! in einer geschlossenen Nachbarzelle ist Geometrie
                            IF (foundx1 == 0 .AND. bzelltyp(k, j, i+1) < 0) THEN
                                xvel = xvel + ucell(1, -bzelltyp(k, j, i+1))
                                yvel = yvel + ucell(2, -bzelltyp(k, j, i+1))
                                zvel = zvel + ucell(3, -bzelltyp(k, j, i+1))
                            END IF
                            IF (foundx2 == 0 .AND. bzelltyp(k, j, i-1) < 0) THEN
                                xvel = xvel + ucell(1, -bzelltyp(k, j, i-1))
                                yvel = yvel + ucell(2, -bzelltyp(k, j, i-1))
                                zvel = zvel + ucell(3, -bzelltyp(k, j, i-1))
                            END IF
                            IF (foundy1 == 0 .AND. bzelltyp(k, j+1, i) < 0) THEN
                                xvel = xvel + ucell(1, -bzelltyp(k, j+1, i))
                                yvel = yvel + ucell(2, -bzelltyp(k, j+1, i))
                                zvel = zvel + ucell(3, -bzelltyp(k, j+1, i))
                            END IF
                            IF (foundy2 == 0 .AND. bzelltyp(k, j-1, i) < 0) THEN
                                xvel = xvel + ucell(1, -bzelltyp(k, j-1, i))
                                yvel = yvel + ucell(2, -bzelltyp(k, j-1, i))
                                zvel = zvel + ucell(3, -bzelltyp(k, j-1, i))
                            END IF
                            IF (foundz1 == 0 .AND. bzelltyp(k+1, j, i) < 0) THEN
                                xvel = xvel + ucell(1, -bzelltyp(k+1, j, i))
                                yvel = yvel + ucell(2, -bzelltyp(k+1, j, i))
                                zvel = zvel + ucell(3, -bzelltyp(k+1, j, i))
                            END IF
                            IF (foundz2 == 0 .AND. bzelltyp(k-1, j, i) < 0) THEN
                                xvel = xvel + ucell(1, -bzelltyp(k-1, j, i))
                                yvel = yvel + ucell(2, -bzelltyp(k-1, j, i))
                                zvel = zvel + ucell(3, -bzelltyp(k-1, j, i))
                            END IF
                            xvel = -xvel/ngeomnbr
                        END IF
                    END IF

                    acoeffstc = (1.0-cx1)*xvel*dy*dz - (1.0-cx2)*xvel*dy*dz + &
                        (1.0-cy1)*yvel*dx*dz - (1.0-cy2)*yvel*dx*dz + &
                        (1.0-cz1)*zvel*dx*dy - (1.0-cz2)*zvel*dx*dy

                    ! Korrigierbare Fluesse/Flaechen
                    ax1 = REAL(foundx1)
                    ax2 = REAL(foundx2)
                    ay1 = REAL(foundy1)
                    ay2 = REAL(foundy2)
                    az1 = REAL(foundz1)
                    az2 = REAL(foundz2)

                    ! Randbedingungen
                    ! "FIX", 2 : FIXED-CONDITION
                    ! "NOS", 5 : NOSLIP
                    ! "SLI", 6 : SLIP
                    IF (i == 3) THEN
                        IF (nfro == 2 .OR. nfro == 5 .OR. nfro == 6) THEN
                            ax2 = 0.0
                        END IF
                    END IF
                    IF (i == ii - 2) THEN
                        IF (nbac == 2 .OR. nbac == 5 .OR. nbac == 6) THEN
                            ax1 = 0.0
                        END IF
                    END IF
                    IF (j == 3) THEN
                        IF (nrgt == 2 .OR. nrgt == 5 .OR. nrgt == 6) THEN
                            ay2 = 0.0
                        END IF
                    END IF
                    IF (j == jj - 2) THEN
                        IF (nlft == 2 .OR. nlft == 5 .OR. nlft == 6) THEN
                            ay1 = 0.0
                        END IF
                    END IF
                    IF (k == 3) THEN
                        IF (nbot == 2 .OR. nbot == 5 .OR. nbot == 6) THEN
                            az2 = 0.0
                        END IF
                    END IF
                    IF (k == kk - 2) THEN
                        IF (ntop == 2 .OR. ntop == 5 .OR. ntop == 6) THEN
                            az1 = 0.0
                        END IF
                    END IF

                    sarea = ax1*au(k, j, i)*ddy(j)*ddz(k) &
                        + ax2*au(k, j  , i-1)*ddy(j)*ddz(k) &
                        + ay1*av(k, j, i)*ddx(i)*ddz(k) &
                        + ay2*av(k, j-1, i)*ddx(i)*ddz(k) &
                        + az1*aw(k, j, i)*ddx(i)*ddy(j) &
                        + az2*aw(k-1, j, i)*ddx(i)*ddy(j)

                    refarea = ((dx + dy + dz)/3.0)**2
                    ! Nur wenn seara > 0 werden die ax... mit den
                    ! tatsaechlichen Flächen belegt ansonsten sind sie
                    ! weiterhin 0 oder 1
                    IF (sarea > refarea*maccur) THEN
                        ax1 = ax1*au(k, j, i)*ddy(j)*ddz(k)
                        ax2 = ax2*au(k, j, i-1)*ddy(j)*ddz(k)
                        ay1 = ay1*av(k, j, i)*ddx(i)*ddz(k)
                        ay2 = ay2*av(k, j-1, i)*ddx(i)*ddz(k)
                        az1 = az1*aw(k, j, i)*ddx(i)*ddy(j)
                        az2 = az2*aw(k-1, j, i)*ddx(i)*ddy(j)
                    END IF

                    IF (pntxpoli + 2 > xpolisize) THEN
                        CALL errr(__FILE__, __LINE__)
                    END IF
                    IF (pntxpolr + 13 > xpolrsize) THEN
                        CALL errr(__FILE__, __LINE__)
                    END IF

                    counter = counter + 1

                    pntxpoli = pntxpoli + 1
                    xpoli(pntxpoli) = foundnr

                    ! TODO: sub2ind
                    pntxpoli = pntxpoli + 1
                    xpoli(pntxpoli) = 1 + (k-1) + (j-1)*kk + (i-1)*jj*kk

                    pntxpolr = pntxpolr + 1
                    xpolr(pntxpolr) = ax1

                    pntxpolr = pntxpolr + 1
                    xpolr(pntxpolr) = ax2

                    pntxpolr = pntxpolr + 1
                    xpolr(pntxpolr) = ay1

                    pntxpolr = pntxpolr + 1
                    xpolr(pntxpolr) = ay2

                    pntxpolr = pntxpolr + 1
                    xpolr(pntxpolr) = az1

                    pntxpolr = pntxpolr + 1
                    xpolr(pntxpolr) = az2

                    pntxpolr = pntxpolr + 1
                    xpolr(pntxpolr) = acoeffstc
                END DO
            END DO
        END DO

        CALL get_imygrid(imygrid, igrid)
        ALLOCATE(fxpoli(imygrid)%arr(pntxpoli))
        ALLOCATE(fxpolr(imygrid)%arr(pntxpolr))
        fxpoli(imygrid)%arr = xpoli(1:pntxpoli)
        fxpolr(imygrid)%arr = xpolr(1:pntxpolr)
        fnblg(imygrid) = counter

        DEALLOCATE(xpoli)
        DEALLOCATE(xpolr)
    END SUBROUTINE fluxcorrection


    SUBROUTINE fluxstencil(igrid, kk, jj, ii, x, y, z, xstag, ystag, &
            zstag, bp, bzelltyp, bconds, icells, nvecs, ucell, &
            compon)

        USE findinterface_mod, ONLY: findinterface2

        ! Subroutine arguments
        INTEGER(intk), INTENT(IN) :: igrid
        INTEGER(intk), INTENT(IN) :: kk, jj, ii
        REAL(realk), INTENT(IN) :: x(ii), y(jj), z(kk)
        REAL(realk), INTENT(IN) :: xstag(ii), ystag(jj), zstag(kk)
        REAL(realk), INTENT(IN) :: bp(kk, jj, ii)
        INTEGER(intk), INTENT(IN) :: bzelltyp(kk, jj, ii)
        INTEGER(intk), INTENT(IN) :: bconds(6)
        INTEGER(intk), INTENT(IN) :: icells
        REAL(realk), CONTIGUOUS, INTENT(IN) :: nvecs(:, :)
        REAL(realk), CONTIGUOUS, INTENT(IN) :: ucell(:, :)
        INTEGER(intk), INTENT(IN) :: compon

        ! Local variables
        INTEGER(intk), PARAMETER :: lsize = 122
        INTEGER(intk), PARAMETER :: velpts = 6
        INTEGER(intk), PARAMETER :: ldofa = 20
        INTEGER(intk) :: found, foundone, imygrid
        INTEGER(intk) :: inlst(lsize), jnlst(lsize), knlst(lsize), nnlst
        INTEGER(intk) :: i, j, k, it, jt, kt, ib, jb, kb
        INTEGER(intk) :: calc, add
        INTEGER(intk) :: pntxpoli, pntxpolr
        INTEGER(intk) :: xpolisize, xpolrsize
        INTEGER(intk), ALLOCATABLE :: xpoli(:)
        REAL(realk), ALLOCATABLE :: xpolr(:), xpolrvel(:)
        REAL(realk) :: avel
        INTEGER(intk) :: counter, itemp(2), jtemp(2), ktemp(2)

        ! To zero for easier debugging
        inlst = 0
        jnlst = 0
        knlst = 0

        ! jk, 17.10.2007: Faktor 1.1, da sonst Dimensionierung zu knapp
        ! jk, 23.03.2008: +max(ni, nj, nk)**2, falls Gitter komplett geblockt,
        !                 aber am Rand
        xpolisize = NINT(icells*1.1 + MAX(kk, jj, ii)**2)*(2+velpts)
        xpolrsize = NINT(icells*1.1 + MAX(kk, jj, ii)**2)*velpts

        ALLOCATE(xpoli(xpolisize))
        ALLOCATE(xpolr(xpolrsize))
        ALLOCATE(xpolrvel(xpolrsize))

        counter = 0

        ! Index stencil for looking for neighbour cells.
        CALL wmindexlistn(inlst, jnlst, knlst, nnlst)

        pntxpoli = 0
        pntxpolr = 0

        kb = 3
        jb = 3
        ib = 3
        IF (compon == 3) kb = 2
        IF (compon == 2) jb = 2
        IF (compon == 1) ib = 2
        DO i = ib, ii-2
            DO j = jb, jj-2
                DO k = kb, kk-2

                    CALL findinterface2(k, j, i, kk, jj, ii, bp, found)

                    ! If interface cell not found, skip rest of loop
                    IF (found /= 1) CYCLE

                    ! If attached p-cell has normal vector
                    itemp = i
                    jtemp = j
                    ktemp = k
                    calc = 0
                    IF (compon == 1) THEN
                        IF (bzelltyp(k, j, i) < 0) THEN
                            calc = calc + 1
                        END IF
                        IF (bzelltyp(k, j, i+1) < 0) THEN
                            calc = calc + 1
                            itemp(calc) = i + 1
                        END IF
                    ELSE IF (compon == 2) THEN
                        IF (bzelltyp(k, j, i) < 0) THEN
                            calc = calc + 1
                        END IF
                        IF (bzelltyp(k, j+1, i) < 0) THEN
                            calc = calc + 1
                            jtemp(calc) = j + 1
                        END IF
                    ELSE IF (compon == 3) then
                        IF (bzelltyp(k, j, i) < 0) THEN
                            calc = calc + 1
                        END IF
                        IF (bzelltyp(k+1, j, i) < 0) THEN
                            calc = calc + 1
                            ktemp(calc) = k + 1
                        END IF
                    END IF

                    ! Check if plane found...
                    IF (calc == 0) THEN
                        ! rescuestencil, da kein Dreieck gefunden wurde
                        CALL wmdocoefflist0(k, j, i, kk, jj, ii, 0.0_realk, &
                            0.0_realk, pntxpoli, xpoli, pntxpolr, xpolr, &
                            xpolrvel)
                        counter = counter + 1
                    ELSE
                        foundone = 0
                        DO add = 1, calc
                            it = itemp(add)
                            jt = jtemp(add)
                            kt = ktemp(add)
                            CALL calcflux(compon, k, j, i, kt, jt, it, &
                                kk, jj, ii, x, y, z, xstag, ystag, zstag, &
                                bp, bzelltyp, bconds, nnlst, inlst, jnlst, &
                                knlst, velpts, nvecs, ucell, ldofa, pntxpoli, &
                                xpoli, pntxpolr, xpolr, xpolrvel, found, &
                                foundone)
                            IF (found == 1) foundone = 1
                        END DO

                        IF (foundone == 1) THEN
                            counter = counter + 1
                        ELSE
                            ! rescuestencil, da zuwenig Punkte gefunden wurden
                            avel = ucell(compon, &
                                -bzelltyp(ktemp(1), jtemp(1), itemp(1)))
                            CALL wmdocoefflist0(k, j, i, kk, jj, ii, avel, &
                                avel, pntxpoli, xpoli, pntxpolr, xpolr, &
                                xpolrvel)
                            counter = counter + 1
                        END IF
                    END IF
                END DO
            END DO
        END DO

        CALL get_imygrid(imygrid, igrid)
        SELECT CASE(compon)
        CASE(1)
            ALLOCATE(uxpoli(imygrid)%arr(pntxpoli))
            ALLOCATE(uxpolr(imygrid)%arr(pntxpolr))
            ALLOCATE(uxpolrvel(imygrid)%arr(pntxpolr))
            ALLOCATE(upoldsol(imygrid)%arr(counter))
            ALLOCATE(uoldsol(imygrid)%arr(counter))
            ALLOCATE(upoldsolvel(imygrid)%arr(counter))
            ALLOCATE(uoldsolvel(imygrid)%arr(counter))
            uxpoli(imygrid)%arr = xpoli(1:pntxpoli)
            uxpolr(imygrid)%arr = xpolr(1:pntxpolr)
            uxpolrvel(imygrid)%arr = xpolrvel(1:pntxpolr)
            upoldsol(imygrid)%arr = 0
            uoldsol(imygrid)%arr = 0.0
            upoldsolvel(imygrid)%arr = 0
            uoldsolvel(imygrid)%arr = 0.0
            unblg(imygrid) = counter
        CASE(2)
            ALLOCATE(vxpoli(imygrid)%arr(pntxpoli))
            ALLOCATE(vxpolr(imygrid)%arr(pntxpolr))
            ALLOCATE(vxpolrvel(imygrid)%arr(pntxpolr))
            ALLOCATE(vpoldsol(imygrid)%arr(counter))
            ALLOCATE(voldsol(imygrid)%arr(counter))
            ALLOCATE(vpoldsolvel(imygrid)%arr(counter))
            ALLOCATE(voldsolvel(imygrid)%arr(counter))
            vxpoli(imygrid)%arr = xpoli(1:pntxpoli)
            vxpolr(imygrid)%arr = xpolr(1:pntxpolr)
            vxpolrvel(imygrid)%arr = xpolrvel(1:pntxpolr)
            vpoldsol(imygrid)%arr = 0
            voldsol(imygrid)%arr = 0.0
            vpoldsolvel(imygrid)%arr = 0
            voldsolvel(imygrid)%arr = 0.0
            vnblg(imygrid) = counter
        CASE(3)
            ALLOCATE(wxpoli(imygrid)%arr(pntxpoli))
            ALLOCATE(wxpolr(imygrid)%arr(pntxpolr))
            ALLOCATE(wxpolrvel(imygrid)%arr(pntxpolr))
            ALLOCATE(wpoldsol(imygrid)%arr(counter))
            ALLOCATE(woldsol(imygrid)%arr(counter))
            ALLOCATE(wpoldsolvel(imygrid)%arr(counter))
            ALLOCATE(woldsolvel(imygrid)%arr(counter))
            wxpoli(imygrid)%arr = xpoli(1:pntxpoli)
            wxpolr(imygrid)%arr = xpolr(1:pntxpolr)
            wxpolrvel(imygrid)%arr = xpolrvel(1:pntxpolr)
            wpoldsol(imygrid)%arr = 0
            woldsol(imygrid)%arr = 0.0
            wpoldsolvel(imygrid)%arr = 0
            woldsolvel(imygrid)%arr = 0.0
            wnblg(imygrid) = counter
        CASE DEFAULT
            CALL errr(__FILE__, __LINE__)
        END SELECT

        DEALLOCATE(xpolrvel)
        DEALLOCATE(xpoli)
        DEALLOCATE(xpolr)
    END SUBROUTINE fluxstencil


    SUBROUTINE wmdocoefflist(k, j, i, kk, jj, ii, acoeffstc, acoeffstcvel, &
            pntxpoli, xpoli, pntxpolr, xpolr, xpolrvel, velpts, istc, jstc, &
            kstc, coeffstc, coeffstcvel)

        ! Subroutine arguments
        ! Same header as wmdocoefflist0 + additional arguments in the end
        INTEGER(intk), INTENT(in) :: k, j, i
        INTEGER(intk), INTENT(in) :: kk, jj, ii
        REAL(realk), INTENT(in) :: acoeffstc, acoeffstcvel

        INTEGER(intk), INTENT(inout) :: pntxpoli
        INTEGER(intk), CONTIGUOUS, INTENT(inout) :: xpoli(:)

        INTEGER(intk), INTENT(inout) :: pntxpolr
        REAL(realk), CONTIGUOUS, INTENT(inout) :: xpolr(:)
        REAL(realk), CONTIGUOUS, INTENT(inout) :: xpolrvel(:)

        INTEGER(intk), INTENT(in) :: velpts
        INTEGER(intk), CONTIGUOUS, INTENT(in) :: istc(:), jstc(:), kstc(:)
        REAL(realk), CONTIGUOUS, INTENT(in) :: coeffstc(:)
        REAL(realk), CONTIGUOUS, INTENT(in) :: coeffstcvel(:)

        ! Local variables
        INTEGER(intk) :: n, xpolisize, xpolrsize

        xpolisize = SIZE(xpoli)
        xpolrsize = SIZE(xpolr)
        IF (SIZE(xpolr) /= SIZE(xpolrvel)) THEN
            CALL errr(__FILE__, __LINE__)
        END IF

        IF (pntxpoli+1+1+velpts > xpolisize) THEN
            WRITE(*, *) 'WMDOCOEFFLIST pntxpoli > xpolisize', &
                pntxpoli+1+1+velpts, xpolisize
            CALL errr(__FILE__, __LINE__)
        END IF
        IF (pntxpolr+velpts+1 > xpolrsize) THEN
            WRITE(*, *) 'WMDOCOEFFLIST pntxpolr > xpolrsize', &
                pntxpolr+velpts+1, xpolrsize
            CALL errr(__FILE__, __LINE__)
        END IF

        pntxpoli = pntxpoli + 1
        ! TODO: sub2ind
        xpoli(pntxpoli) = 1+(k-1)+(j-1)*kk+(i-1)*jj*kk

        pntxpoli = pntxpoli + 1
        xpoli(pntxpoli) = velpts

        DO n = 1, velpts
            pntxpoli = pntxpoli + 1
            ! TODO: sub2ind
            xpoli(pntxpoli) = 1+(kstc(n)-1)+(jstc(n)-1)*kk+(istc(n)-1)*jj*kk

            pntxpolr = pntxpolr + 1
            xpolr(pntxpolr)    = coeffstc(n)
            xpolrvel(pntxpolr) = coeffstcvel(n)
        END DO

        pntxpolr = pntxpolr + 1
        xpolr(pntxpolr) = acoeffstc
        xpolrvel(pntxpolr) = acoeffstcvel
    END SUBROUTINE wmdocoefflist


    SUBROUTINE wmdocoefflist0(k, j, i, kk, jj, ii, acoeffstc, acoeffstcvel, &
            pntxpoli, xpoli, pntxpolr, xpolr, xpolrvel)

        ! Subroutine arguments
        INTEGER(intk), INTENT(in) :: k, j, i
        INTEGER(intk), INTENT(in) :: kk, jj, ii
        REAL(realk), INTENT(in) :: acoeffstc, acoeffstcvel

        INTEGER(intk), INTENT(inout) :: pntxpoli
        INTEGER(intk), CONTIGUOUS, INTENT(inout) :: xpoli(:)

        INTEGER(intk), INTENT(inout) :: pntxpolr
        REAL(realk), CONTIGUOUS, INTENT(inout) :: xpolr(:)
        REAL(realk), CONTIGUOUS, INTENT(inout) :: xpolrvel(:)

        CALL wmdocoefflist(k, j, i, kk, jj, ii, acoeffstc, acoeffstcvel, &
            pntxpoli, xpoli, pntxpolr, xpolr, xpolrvel, 0, [0], [0], [0], &
            [0.0_realk], [0.0_realk])
    END SUBROUTINE wmdocoefflist0


    SUBROUTINE calcflux(compon, k, j, i, kt, jt, it, kk, jj, ii, x, y, z, &
            xstag, ystag, zstag, bp, bzelltyp, bconds, nnlst, inlst, jnlst, &
            knlst, velpts, nvecs, ucell, ldofa, pntxpoli, xpoli, pntxpolr, &
            xpolr, xpolrvel, found, foundone)

        ! NOTE: This subroutine is very close to tscastencilcoeff from
        ! gc_scastencils_mod.F90

        ! Subroutine arguments
        INTEGER(intk), INTENT(in) :: compon
        INTEGER(intk), INTENT(in) :: k, j, i
        INTEGER(intk), INTENT(in) :: kt, jt, it
        INTEGER(intk), INTENT(in) :: kk, jj, ii
        REAL(realk), INTENT(in) :: x(ii), y(jj), z(kk)
        REAL(realk), INTENT(in) :: xstag(ii), ystag(jj), zstag(kk)
        REAL(realk), INTENT(in) :: bp(kk, jj, ii)
        INTEGER(intk), INTENT(in) :: bzelltyp(kk, jj, ii)
        INTEGER(intk), INTENT(in) :: bconds(6)
        INTEGER(intk), INTENT(in) :: nnlst
        INTEGER(intk), CONTIGUOUS, INTENT(in) :: inlst(:), jnlst(:), knlst(:)
        INTEGER(intk), INTENT(in) :: velpts
        REAL(realk), CONTIGUOUS, INTENT(in) :: nvecs(:, :), ucell(:, :)
        INTEGER(intk), INTENT(in) :: ldofa
        INTEGER(intk), INTENT(inout) :: pntxpoli
        INTEGER(intk), CONTIGUOUS, INTENT(inout) :: xpoli(:)
        INTEGER(intk), INTENT(inout) :: pntxpolr
        REAL(realk), CONTIGUOUS, INTENT(inout) :: xpolr(:)
        REAL(realk), CONTIGUOUS, INTENT(inout) :: xpolrvel(:)
        INTEGER(intk), INTENT(inout) :: found, foundone

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
            WRITE(*, *) 'calcflux: ijk, compon, ', i, j, k, compon
            WRITE(*, *) 'calcflux:  chosenvelpts > velpts', chosenvelpts, velpts
            found = -2
            RETURN
        END IF

        ! keine offenen Nachbarn: rescue-Stencil
        IF (chosenvelpts == 0) THEN
            found = -3
            RETURN
        END IF

        ! Multiply matrix and calculate coefficients --> found!
        CALL wmmultimatrix(compon, k, j, i, kt, jt, it, kk, jj, ii, &
            x, y, z, xstag, ystag, zstag, bzelltyp, istc, jstc, kstc, &
            chosenvelpts, nvecs, ucell, coeffstc, coeffstcvel, &
            acoeffstc, acoeffstcvel, found)

        ! Abfangen schlechter Werte, rescue-Stencil
        DO mm = 1, chosenvelpts
            IF (ABS(coeffstcvel(mm)) > 12.0) found = -4
        END DO

        ! Extend coefficient list
        IF (found == 1) THEN
            IF (foundone == 0) THEN
                CALL wmdocoefflist(k, j, i, kk, jj, ii, acoeffstc, &
                    acoeffstcvel, pntxpoli, xpoli, pntxpolr, xpolr, xpolrvel, &
                    chosenvelpts, istc, jstc, kstc, coeffstc, coeffstcvel)
            ELSE IF (foundone == 1) THEN
                CALL wmaddcoefflist(chosenvelpts, xpolr, xpolrvel, pntxpolr, &
                    acoeffstc, acoeffstcvel, coeffstc, coeffstcvel)
            END IF
        END IF
    END SUBROUTINE calcflux


    PURE SUBROUTINE wmxpol(kk, jj, ii, pntxpoli, pntxpolr, xpoli, xpolr, &
            var, intcell, flux)

        ! Subroutine arguments
        INTEGER(intk), INTENT(in) :: kk, jj, ii
        INTEGER(intk), INTENT(inout) :: pntxpoli, pntxpolr
        INTEGER(intk), INTENT(in), CONTIGUOUS :: xpoli(:)
        REAL(realk), INTENT(in), CONTIGUOUS :: xpolr(:)
        REAL(realk), INTENT(in) :: var(kk*jj*ii)
        INTEGER(intk), INTENT(out) :: intcell
        REAL(realk), INTENT(out) :: flux

        ! Local variables
        INTEGER(intk) :: stencils, n, stcell
        REAL(realk) :: coeff

        intcell = xpoli(pntxpoli)
        pntxpoli = pntxpoli + 1

        stencils = xpoli(pntxpoli)
        pntxpoli = pntxpoli + 1

        flux = 0.0

        DO n = 1, stencils
            stcell = xpoli(pntxpoli)
            pntxpoli = pntxpoli + 1

            coeff = xpolr(pntxpolr)
            pntxpolr = pntxpolr + 1

            flux = flux + var(stcell)*coeff
        END DO

        coeff = xpolr(pntxpolr)
        pntxpolr = pntxpolr + 1

        flux = flux + coeff
    END SUBROUTINE wmxpol


    SUBROUTINE wmxpolsol(cmp, kk, jj, ii, u, v, w, xpoli, xpolr, ncells)

        ! Subroutine arguments
        INTEGER(intk), INTENT(in) :: cmp
        INTEGER(intk), INTENT(in) :: kk, jj, ii
        REAL(realk), INTENT(inout) :: u(kk*jj*ii)
        REAL(realk), INTENT(inout) :: v(kk*jj*ii)
        REAL(realk), INTENT(inout) :: w(kk*jj*ii)
        INTEGER(intk), INTENT(in), CONTIGUOUS :: xpoli(:)
        REAL(realk), INTENT(in), CONTIGUOUS :: xpolr(:)
        INTEGER(intk), INTENT(in) :: ncells

        ! Local variables
        INTEGER(intk) :: i, intcell
        INTEGER(intk) :: pntxpoli, pntxpolr
        REAL(realk) :: flux

        pntxpoli = 1
        pntxpolr = 1

        SELECT CASE(cmp)
        CASE (1)
            DO i = 1, ncells
                CALL wmxpol(kk, jj, ii, pntxpoli, pntxpolr, &
                    xpoli, xpolr, u, intcell, flux)
                u(intcell) = flux
            END DO
        CASE (2)
            DO i = 1, ncells
                CALL wmxpol(kk, jj, ii, pntxpoli, pntxpolr, &
                    xpoli, xpolr, v, intcell, flux)
                v(intcell) = flux
            END DO
        CASE (3)
            DO i = 1, ncells
                CALL wmxpol(kk, jj, ii, pntxpoli, pntxpolr, &
                    xpoli, xpolr, w, intcell, flux)
                w(intcell) = flux
            END DO
        CASE DEFAULT
            CALL errr(__FILE__, __LINE__)
        END SELECT
    END SUBROUTINE wmxpolsol


    SUBROUTINE wmxpolsolrlx(cmp, kk, jj, ii, u, v, w, ncells, &
            oldsol, poldsol)

        ! Subroutine arguments
        INTEGER(intk), INTENT(in) :: cmp
        INTEGER(intk), INTENT(in) :: kk, jj, ii
        REAL(realk), INTENT(inout) :: u(kk*jj*ii)
        REAL(realk), INTENT(inout) :: v(kk*jj*ii)
        REAL(realk), INTENT(inout) :: w(kk*jj*ii)
        INTEGER(intk), INTENT(in) :: ncells
        REAL(realk), INTENT(in) :: oldsol(ncells)
        INTEGER(intk), INTENT(in) :: poldsol(ncells)

        ! Local variables
        INTEGER(intk) :: i, idx

        SELECT CASE(cmp)
        CASE (1)
            DO i = 1, ncells
                idx = poldsol(i)
                u(idx) = oldsol(i)
            END DO
        CASE (2)
            DO i = 1, ncells
                idx = poldsol(i)
                v(idx) = oldsol(i)
            END DO
        CASE (3)
            DO i = 1, ncells
                idx = poldsol(i)
                w(idx) = oldsol(i)
            END DO
        CASE DEFAULT
            CALL errr(__FILE__, __LINE__)
        END SELECT
    END SUBROUTINE wmxpolsolrlx


    SUBROUTINE wmxpolsolsav(cmp, kk, jj, ii, u, v, w, xpoli, ncells, &
            oldsol, poldsol)

        ! Subroutine arguments
        INTEGER(intk), INTENT(in) :: cmp
        INTEGER(intk), INTENT(in) :: kk, jj, ii
        REAL(realk), INTENT(in) :: u(kk*jj*ii)
        REAL(realk), INTENT(in) :: v(kk*jj*ii)
        REAL(realk), INTENT(in) :: w(kk*jj*ii)
        INTEGER(intk), INTENT(in), CONTIGUOUS :: xpoli(:)
        INTEGER(intk), INTENT(in) :: ncells
        REAL(realk), INTENT(inout) :: oldsol(ncells)
        INTEGER(intk), INTENT(inout) :: poldsol(ncells)

        ! Local variables
        INTEGER(intk) :: i, idx
        INTEGER(intk) :: pntxpoli
        INTEGER(intk) :: stencils

        pntxpoli = 1

        SELECT CASE(cmp)
        CASE (1)
            DO i = 1, ncells
                idx = xpoli(pntxpoli)
                pntxpoli = pntxpoli + 1

                stencils = xpoli(pntxpoli)
                pntxpoli = pntxpoli + 1

                pntxpoli = pntxpoli + stencils

                poldsol(i) = idx
                oldsol(i) = u(idx)
            END DO
        CASE (2)
            DO i = 1, ncells
                idx = xpoli(pntxpoli)
                pntxpoli = pntxpoli + 1

                stencils = xpoli(pntxpoli)
                pntxpoli = pntxpoli + 1

                pntxpoli = pntxpoli + stencils

                poldsol(i) = idx
                oldsol(i) = v(idx)
            END DO
        CASE (3)
            DO i = 1, ncells
                idx = xpoli(pntxpoli)
                pntxpoli = pntxpoli + 1

                stencils = xpoli(pntxpoli)
                pntxpoli = pntxpoli + 1

                pntxpoli = pntxpoli + stencils

                poldsol(i) = idx
                oldsol(i) = w(idx)
            END DO
        CASE DEFAULT
            CALL errr(__FILE__, __LINE__)
        END SELECT
    END SUBROUTINE wmxpolsolsav


    SUBROUTINE wmxpolsolcorr(kk, jj, ii, u, v, w, ddx, ddy, ddz, &
            xpoli, xpolr, ncells)

        ! Subroutine arguments
        INTEGER(intk), INTENT(in) :: kk, jj, ii
        REAL(realk), INTENT(inout) :: u(kk*jj*ii)
        REAL(realk), INTENT(inout) :: v(kk*jj*ii)
        REAL(realk), INTENT(inout) :: w(kk*jj*ii)
        REAL(realk), INTENT(in) :: ddx(ii), ddy(jj), ddz(kk)
        INTEGER(intk), INTENT(in), CONTIGUOUS :: xpoli(:)
        REAL(realk), INTENT(in), CONTIGUOUS :: xpolr(:)
        INTEGER(intk), INTENT(in) :: ncells

        ! Local variables
        INTEGER(intk) :: k, j, i
        INTEGER(intk) :: cellcount, intcell
        INTEGER(intk) :: ishift, jshift, kshift, foundnr
        INTEGER(intk) :: pntxpoli, pntxpolr
        REAL(realk) :: div, acoeffstc, sarea
        REAL(realk) :: ax1, ax2, ay1, ay2, az1, az2

        ishift = jj*kk
        jshift = kk
        kshift = 1

        pntxpoli = 1
        pntxpolr = 1
        DO cellcount = 1, ncells
            foundnr = xpoli(pntxpoli)
            pntxpoli = pntxpoli + 1

            intcell = xpoli(pntxpoli)
            pntxpoli = pntxpoli + 1

            ax1 = xpolr(pntxpolr)
            pntxpolr = pntxpolr + 1
            ax2 = xpolr(pntxpolr)
            pntxpolr = pntxpolr + 1
            ay1 = xpolr(pntxpolr)
            pntxpolr = pntxpolr + 1
            ay2 = xpolr(pntxpolr)
            pntxpolr = pntxpolr + 1
            az1 = xpolr(pntxpolr)
            pntxpolr = pntxpolr + 1
            az2 = xpolr(pntxpolr)
            pntxpolr = pntxpolr + 1

            acoeffstc = xpolr(pntxpolr)
            pntxpolr = pntxpolr + 1

            CALL ind2sub(intcell, k, j, i, kk, jj, ii)

            div = ddy(j)*ddz(k)*(u(intcell) - u(intcell - ishift)) &
                + ddx(i)*ddz(k)*(v(intcell) - v(intcell - jshift)) &
                + ddx(i)*ddy(j)*(w(intcell) - w(intcell - kshift)) &
                + acoeffstc

            sarea = ax1 + ax2 + ay1 + ay2 + az1 + az2
            IF (sarea < TINY(1.0_realk)) THEN
                CALL errr(__FILE__, __LINE__)
            END IF
            div = div/sarea

            u(intcell) = u(intcell) - ax1*div/(ddy(j)*ddz(k))
            u(intcell-ishift) = u(intcell-ishift) + ax2*div/(ddy(j)*ddz(k))
            v(intcell) = v(intcell) - ay1*div/(ddx(i)*ddz(k))
            v(intcell-jshift) = v(intcell-jshift) + ay2*div/(ddx(i)*ddz(k))
            w(intcell) = w(intcell) - az1*div/(ddx(i)*ddy(j))
            w(intcell-kshift) = w(intcell-kshift) + az2*div/(ddx(i)*ddy(j))
        END DO
    END SUBROUTINE wmxpolsolcorr


    SUBROUTINE wmxpolquadvel(igrid, cmp, kk, jj, ii, ityp, u, v, w)

        ! Subroutine arguments
        INTEGER(intk), INTENT(in) :: igrid, cmp, kk, jj, ii
        CHARACTER(len=1), INTENT(in) :: ityp
        REAL(realk), INTENT(inout) :: u(kk*jj*ii)
        REAL(realk), INTENT(inout) :: v(kk*jj*ii)
        REAL(realk), INTENT(inout) :: w(kk*jj*ii)

        ! Local variables
        INTEGER(intk) :: imygrid
        INTEGER(intk) :: ncells
        INTEGER(intk), POINTER, CONTIGUOUS :: xpoli(:), poldsol(:)
        REAL(realk), POINTER, CONTIGUOUS :: xpolr(:), oldsol(:)

        CALL get_imygrid(imygrid, igrid)
        SELECT CASE (cmp)
        CASE(1)
            xpoli => uxpoli(imygrid)%arr
            xpolr => uxpolrvel(imygrid)%arr
            poldsol => upoldsolvel(imygrid)%arr
            oldsol => uoldsolvel(imygrid)%arr
            ncells = unblg(imygrid)
        CASE(2)
            xpoli => vxpoli(imygrid)%arr
            xpolr => vxpolrvel(imygrid)%arr
            poldsol => vpoldsolvel(imygrid)%arr
            oldsol => voldsolvel(imygrid)%arr
            ncells = vnblg(imygrid)
        CASE(3)
            xpoli => wxpoli(imygrid)%arr
            xpolr => wxpolrvel(imygrid)%arr
            poldsol => wpoldsolvel(imygrid)%arr
            oldsol => woldsolvel(imygrid)%arr
            ncells = wnblg(imygrid)
        CASE DEFAULT
            CALL errr(__FILE__, __LINE__)
        END SELECT

        SELECT CASE (ityp)
        CASE("X")
            CALL wmxpolsol(cmp, kk, jj, ii, u, v, w, xpoli, xpolr, ncells)
        CASE("Y")
            CALL wmxpolsolrlx(cmp, kk, jj, ii, u, v, w, ncells, &
                oldsol, poldsol)
        CASE("Z")
            CALL wmxpolsolsav(cmp, kk, jj, ii, u, v, w, xpoli, &
                ncells, oldsol, poldsol)
        CASE DEFAULT
            CALL errr(__FILE__, __LINE__)
        END SELECT
    END SUBROUTINE wmxpolquadvel


    SUBROUTINE wmxpolquad(igrid, cmp, kk, jj, ii, ityp, u, v, w)

        ! Subroutine arguments
        INTEGER(intk), INTENT(in) :: igrid, cmp, kk, jj, ii
        CHARACTER(len=1), INTENT(in) :: ityp
        REAL(realk), INTENT(inout) :: u(kk*jj*ii)
        REAL(realk), INTENT(inout) :: v(kk*jj*ii)
        REAL(realk), INTENT(inout) :: w(kk*jj*ii)

        ! Local variables
        INTEGER(intk) :: imygrid
        INTEGER(intk) :: ncells
        INTEGER(intk), POINTER, CONTIGUOUS :: xpoli(:), poldsol(:)
        REAL(realk), POINTER, CONTIGUOUS :: xpolr(:), oldsol(:)

        CALL get_imygrid(imygrid, igrid)
        SELECT CASE (cmp)
        CASE(1)
            xpoli => uxpoli(imygrid)%arr
            xpolr => uxpolr(imygrid)%arr
            poldsol => upoldsol(imygrid)%arr
            oldsol => uoldsol(imygrid)%arr
            ncells = unblg(imygrid)
        CASE(2)
            xpoli => vxpoli(imygrid)%arr
            xpolr => vxpolr(imygrid)%arr
            poldsol => vpoldsol(imygrid)%arr
            oldsol => voldsol(imygrid)%arr
            ncells = vnblg(imygrid)
        CASE(3)
            xpoli => wxpoli(imygrid)%arr
            xpolr => wxpolr(imygrid)%arr
            poldsol => wpoldsol(imygrid)%arr
            oldsol => woldsol(imygrid)%arr
            ncells = wnblg(imygrid)
        CASE DEFAULT
            CALL errr(__FILE__, __LINE__)
        END SELECT

        SELECT CASE (ityp)
        CASE("F")
            CALL wmxpolsol(cmp, kk, jj, ii, u, v, w, xpoli, xpolr, ncells)
        CASE("R")
            CALL wmxpolsolrlx(cmp, kk, jj, ii, u, v, w, ncells, &
                oldsol, poldsol)
        CASE("S")
            CALL wmxpolsolsav(cmp, kk, jj, ii, u, v, w, xpoli, &
                ncells, oldsol, poldsol)
        CASE DEFAULT
            CALL errr(__FILE__, __LINE__)
        END SELECT
    END SUBROUTINE wmxpolquad


    SUBROUTINE wmxpolquadfcorr(igrid, kk, jj, ii, ityp, u, v, w)
        ! Subroutine arguments
        INTEGER(intk), INTENT(in) :: igrid, kk, jj, ii
        CHARACTER(len=1), INTENT(in) :: ityp
        REAL(realk), INTENT(inout) :: u(kk*jj*ii)
        REAL(realk), INTENT(inout) :: v(kk*jj*ii)
        REAL(realk), INTENT(inout) :: w(kk*jj*ii)

        ! Local variables
        INTEGER(intk) :: imygrid
        INTEGER(intk) :: ncells
        INTEGER(intk), POINTER, CONTIGUOUS :: xpoli(:)
        REAL(realk), POINTER, CONTIGUOUS :: xpolr(:)
        REAL(realk), POINTER, CONTIGUOUS :: ddx(:), ddy(:), ddz(:)

        CALL get_imygrid(imygrid, igrid)
        xpoli => fxpoli(imygrid)%arr
        xpolr => fxpolr(imygrid)%arr
        ncells = fnblg(imygrid)

        CALL get_fieldptr(ddx, "DDX", igrid)
        CALL get_fieldptr(ddy, "DDY", igrid)
        CALL get_fieldptr(ddz, "DDZ", igrid)

        CALL wmxpolsolcorr(kk, jj, ii, u, v, w, ddx, ddy, ddz, &
            xpoli, xpolr, ncells)
    END SUBROUTINE wmxpolquadfcorr


    SUBROUTINE setpointvalues(pwu, pwv, pww, u, v, w, comp_new)
        ! Subroutine arguments
        TYPE(field_t), INTENT(inout) :: pwu, pwv, pww
        TYPE(field_t), INTENT(in) :: u, v, w
        LOGICAL, INTENT(in) :: comp_new

        ! Local variables
        INTEGER(intk) :: ilevel

        CALL start_timer(342)

        ! Copy the velocity field u, v, w into the point value
        ! velocity field pwu, pwv, pww
        pwu%arr = u%arr
        pwv%arr = v%arr
        pww%arr = w%arr

        ! Copy buffers from u, v, w to pwu, pwv, pww
        pwu%buffers = u%buffers
        pwv%buffers = v%buffers
        pww%buffers = w%buffers

        ! Set the immersed boundary ghost cell values of the
        ! point value velocity field
        IF (comp_new) THEN
            ! Extrapolate the immersed boundary ghost cell values from
            ! the velocity field.
            CALL setpointvalues_all('X', pwu, pwv, pww)
        ELSE
            ! Restore the ghost cell values from the point value buffer
            CALL setpointvalues_all('Y', pwu, pwv, pww)
        END IF

        DO ilevel = minlevel, maxlevel
            CALL connect(ilevel, 1, v1=pwu, v2=pwv, v3=pww, geom=.TRUE.)
            CALL bound_flow%bound(ilevel, pwu, pwv, pww)
        END DO

        IF (comp_new) THEN
            ! Save the immersed boundary ghost cell values into a buffer
            CALL setpointvalues_all('Z', pwu, pwv, pww)
        END IF

        CALL stop_timer(342)
    END SUBROUTINE setpointvalues


    SUBROUTINE setpointvalues_all(ityp, pwu, pwv, pww)
        ! Subroutine arguments
        CHARACTER(len=1), INTENT(in) :: ityp
        TYPE(field_t), INTENT(inout) :: pwu, pwv, pww

        ! Local variables
        INTEGER(intk) :: i, igrid, ip3
        INTEGER(intk) :: kk, jj, ii

        DO i = 1, nmygrids
            igrid = mygrids(i)
            CALL get_mgdims(kk, jj, ii, igrid)
            CALL get_ip3(ip3, igrid)
            CALL wmxpolquadvel(igrid, 1, kk, jj, ii, ityp, pwu%arr(ip3), &
                pwv%arr(ip3), pww%arr(ip3))
            CALL wmxpolquadvel(igrid, 2, kk, jj, ii, ityp, pwu%arr(ip3), &
                pwv%arr(ip3), pww%arr(ip3))
            CALL wmxpolquadvel(igrid, 3, kk, jj, ii, ityp, pwu%arr(ip3), &
                pwv%arr(ip3), pww%arr(ip3))
        END DO
    END SUBROUTINE setpointvalues_all


    SUBROUTINE setibvalues(u, v, w)
        ! Subroutine arguments
        TYPE(field_t), INTENT(inout) :: u, v, w

        ! Local variables
        INTEGER(intk) :: ilevel, irepeat

        CALL start_timer(340)

        DO ilevel = minlevel, maxlevel
            CALL parent(ilevel, u, v, w)
            CALL bound_flow%bound(ilevel, u, v, w)
            CALL connect(ilevel, 2, v1=u, v2=v, v3=w)

            CALL setibvalues_level(ilevel, 'F', u, v, w)
            CALL bound_flow%bound(ilevel, u, v, w)

            CALL connect(ilevel, 1, v1=u, v2=v, v3=w, corners=.TRUE.)

            CALL setibvalues_level(ilevel, 'C', u, v, w)
            CALL bound_flow%bound(ilevel, u, v, w)

            CALL connect(ilevel, 1, v1=u, v2=v, v3=w, corners=.TRUE.)
        END DO

        DO ilevel = maxlevel, minlevel, -1
            DO irepeat = 1, 2
                CALL ftoc(ilevel, u%arr, u%arr, 'U')
                CALL ftoc(ilevel, v%arr, v%arr, 'V')
                CALL ftoc(ilevel, w%arr, w%arr, 'W')

                IF (ilevel > minlevel) THEN
                    IF (irepeat == 1) THEN
                        CALL connect(ilevel-1, 1, v1=u, v2=v, &
                            v3=w, normal=.TRUE., forward=-1, ityp='Y')
                    ELSE IF (irepeat == 2) THEN
                        CALL connect(ilevel-1, 1, v1=u, v2=v, &
                            v3=w, corners=.TRUE.)
                    ELSE
                        CALL errr(__FILE__, __LINE__)
                    END IF
                END IF
            END DO
        END DO

        DO ilevel = minlevel, maxlevel
            CALL parent(ilevel, u, v, w)
            CALL bound_flow%bound(ilevel, u, v, w)

            CALL connect(ilevel, 1, v1=u, v2=v, v3=w, corners=.TRUE.)

            CALL setibvalues_level(ilevel, 'C', u, v, w)
            CALL bound_flow%bound(ilevel, u, v, w)
            CALL connect(ilevel, 1, v1=u, v2=v, v3=w, corners=.TRUE.)

            ! The computed fluxes are saved
            CALL setibvalues_level(ilevel, 'S', u, v, w)
        END DO

        CALL stop_timer(340)
    END SUBROUTINE setibvalues


    SUBROUTINE getibvalues(u, v, w)
        ! Subroutine arguments
        TYPE(field_t), INTENT(inout) :: u, v, w

        ! Local variables
        INTEGER(intk) :: ilevel

        CALL start_timer(341)
        DO ilevel = minlevel, maxlevel
            CALL setibvalues_level(ilevel, 'R', u, v, w)
        END DO
        CALL stop_timer(341)
    END SUBROUTINE getibvalues


    SUBROUTINE setibvalues_level(ilevel, ityp, u, v, w)
        ! Subroutine arguments
        INTEGER(intk), INTENT(in) :: ilevel
        CHARACTER(len=1), INTENT(in) :: ityp
        TYPE(field_t), INTENT(inout) :: u, v, w

        ! Local variables
        INTEGER(intk) :: i, igrid, ip3
        INTEGER(intk) :: kk, jj, ii

        DO i = 1, nmygridslvl(ilevel)
            igrid = mygridslvl(i, ilevel)
            CALL get_mgdims(kk, jj, ii, igrid)
            CALL get_ip3(ip3, igrid)
            SELECT CASE(ityp)
            CASE("F", "R", "S")
                CALL wmxpolquad(igrid, 1, kk, jj, ii, ityp, u%arr(ip3), &
                    v%arr(ip3), w%arr(ip3))
                CALL wmxpolquad(igrid, 2, kk, jj, ii, ityp, u%arr(ip3), &
                    v%arr(ip3), w%arr(ip3))
                CALL wmxpolquad(igrid, 3, kk, jj, ii, ityp, u%arr(ip3), &
                    v%arr(ip3), w%arr(ip3))
            CASE ("C")
                CALL wmxpolquadfcorr(igrid, kk, jj, ii, ityp, u%arr(ip3), &
                    v%arr(ip3), w%arr(ip3))
            CASE DEFAULT
                CALL errr(__FILE__, __LINE__)
            END SELECT
        END DO
    END SUBROUTINE setibvalues_level


    SUBROUTINE setsdivfield(sdiv_f)
        ! Subroutine arguments
        TYPE(field_t), INTENT(inout) :: sdiv_f

        ! Local variables
        TYPE(field_t), POINTER :: ddx_f, ddy_f, ddz_f
        REAL(realk), POINTER, CONTIGUOUS :: sdiv(:, :, :), ddx(:), &
            ddy(:), ddz(:)
        INTEGER(intk) :: i, igrid, ilevel
        INTEGER(intk) :: kk, jj, ii

        CALL get_field(ddx_f, "DDX")
        CALL get_field(ddy_f, "DDY")
        CALL get_field(ddz_f, "DDZ")

        DO i = 1, nmygrids
            igrid = mygrids(i)
            CALL get_mgdims(kk, jj, ii, igrid)
            CALL sdiv_f%get_ptr(sdiv, igrid)
            CALL ddx_f%get_ptr(ddx, igrid)
            CALL ddy_f%get_ptr(ddy, igrid)
            CALL ddz_f%get_ptr(ddz, igrid)
            CALL setsdivfield_grid(igrid, kk, jj, ii, sdiv, ddx, ddy, ddz)
        END DO

        DO ilevel = maxlevel, minlevel+1, -1
            CALL ftoc(ilevel, sdiv_f%arr, sdiv_f%arr, 'S')
        END DO
    END SUBROUTINE setsdivfield


    SUBROUTINE setsdivfield_grid(igrid, kk, jj, ii, sdiv, ddx, ddy, ddz)
        ! Subroutine arguments
        INTEGER(intk), INTENT(in) :: igrid, kk, jj, ii
        REAL(realk), INTENT(out) :: sdiv(kk, jj, ii)
        REAL(realk), INTENT(in) :: ddx(ii)
        REAL(realk), INTENT(in) :: ddy(jj)
        REAL(realk), INTENT(in) :: ddz(kk)

        ! Local variables
        INTEGER(intk) :: pntxpoli, pntxpolr, cellcount, intcell
        INTEGER(intk) :: imygrid
        INTEGER(intk) :: k, j, i
        REAL(realk) :: acoeffstc

        ! Initialize INTENT(out)
        sdiv = 0.0

        pntxpoli = 1
        pntxpolr = 1
        CALL get_imygrid(imygrid, igrid)
        DO cellcount = 1, fnblg(imygrid)
            pntxpoli = pntxpoli + 1

            intcell = fxpoli(imygrid)%arr(pntxpoli)
            pntxpoli = pntxpoli + 1

            ! skip ax1, ...
            pntxpolr = pntxpolr + 6
            acoeffstc = fxpolr(imygrid)%arr(pntxpolr)
            pntxpolr = pntxpolr + 1

            CALL ind2sub(intcell, k, j, i, kk, jj, ii)
            sdiv(k, j, i) = acoeffstc/(ddx(i)*ddy(j)*ddz(k))
        END DO
    END SUBROUTINE setsdivfield_grid


    SUBROUTINE writestencils()
        USE MPI_f08, ONLY: MPI_Barrier, MPI_COMM_WORLD

        ! Subroutine arguments
        ! none...

        ! Local variables
        INTEGER(intk) :: i, igrid, comp

        IF (myid == 0) THEN
            CALL create_directory("STENCILS")
        END IF
        CALL MPI_Barrier(MPI_COMM_WORLD)

        DO i = 1, nmygrids
            igrid = mygrids(i)
            DO comp = 1, 4
                CALL writestencils_grid(comp, igrid)
            END DO
        END DO
        CALL MPI_Barrier(MPI_COMM_WORLD)

        IF (myid == 0) THEN
            CALL pvtk_directory("STENCILS", "stecil-u")
            CALL pvtk_directory("STENCILS", "stecil-v")
            CALL pvtk_directory("STENCILS", "stecil-w")
            CALL pvtk_directory("STENCILS", "stecil-f")
        END IF
    END SUBROUTINE writestencils


    SUBROUTINE writestencils_grid(comp, igrid)
        ! Subroutine arguments
        INTEGER(intk), INTENT(in) :: comp, igrid

        ! Local variables
        INTEGER(intk) :: kk, jj, ii, imygrid, ncells
        INTEGER(intk), POINTER, CONTIGUOUS :: xpoli(:)
        REAL(realk), POINTER, CONTIGUOUS :: xpolr(:)
        REAL(realk), POINTER, CONTIGUOUS :: x(:), y(:), z(:)
        REAL(realk), POINTER, CONTIGUOUS :: dx(:), dy(:), dz(:)

        CALL get_imygrid(imygrid, igrid)
        SELECT CASE (comp)
        CASE(1)
            xpoli => uxpoli(imygrid)%arr
            xpolr => uxpolrvel(imygrid)%arr
            ncells = unblg(imygrid)
        CASE(2)
            xpoli => vxpoli(imygrid)%arr
            xpolr => vxpolrvel(imygrid)%arr
            ncells = vnblg(imygrid)
        CASE(3)
            xpoli => wxpoli(imygrid)%arr
            xpolr => wxpolrvel(imygrid)%arr
            ncells = wnblg(imygrid)
        CASE(4)
            xpoli => fxpoli(imygrid)%arr
            xpolr => fxpolr(imygrid)%arr
            ncells = fnblg(imygrid)
        CASE DEFAULT
            CALL errr(__FILE__, __LINE__)
        END SELECT

        ! If no cells are present we return here...
        IF (ncells == 0) RETURN

        CALL get_fieldptr(x, "X", igrid)
        CALL get_fieldptr(y, "Y", igrid)
        CALL get_fieldptr(z, "Z", igrid)
        CALL get_fieldptr(dx, "DX", igrid)
        CALL get_fieldptr(dy, "DY", igrid)
        CALL get_fieldptr(dz, "DZ", igrid)
        CALL get_mgdims(kk, jj, ii, igrid)
        CALL writestencilsvtk(comp, igrid, kk, jj, ii, ncells, xpoli, xpolr, &
            x, y, z, dz, dy, dz)
    END SUBROUTINE writestencils_grid


    SUBROUTINE writestencilsvtk(comp, igrid, kk, jj, ii, nblgcells, xpoli, &
            xpolr, x, y, z, dx, dy, dz)

        ! Subroutine arguments
        INTEGER(intk), INTENT(in) :: comp, igrid
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
        CHARACTER(len=8) :: prefix

        istag = 0
        SELECT CASE (comp)
        CASE(1)
            prefix = "stecil-u"
            istag(comp) = 1
        CASE(2)
            prefix = "stecil-v"
            istag(comp) = 1
        CASE(3)
            prefix = "stecil-w"
            istag(comp) = 1
        CASE(4)
            prefix = "stecil-f"
        CASE DEFAULT
            CALL errr(__FILE__, __LINE__)
        END SELECT

        WRITE(filename, '("STENCILS/", A, "-igrid-", I0, ".vtk")') &
                TRIM(prefix), igrid

        OPEN(newunit=unit, file=TRIM(filename))
        WRITE(unit, '("# vtk DataFile Version 3.0")')
        WRITE(unit, '("MGLET SUBROUTINE writestencilsvtk")')
        WRITE(unit, '("ASCII")')
        WRITE(unit, '("DATASET UNSTRUCTURED_GRID")')

        IF (comp <= 3) THEN
            ! Count nstencils
            nstencils = 0
            pntxpoli = 1
            DO cellcount = 1, nblgcells
                ! not needing this one
                !intcell = xpoli(pntxpoli)
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

        ELSE IF (comp == 4) THEN
            WRITE(unit, '("POINTS ", I0, " float")') nblgcells

            ! Write points
            pntxpoli = 1
            DO cellcount = 1, nblgcells
                pntxpoli = pntxpoli + 1

                intcell = xpoli(pntxpoli)
                CALL ind2sub(intcell, k, j, i, kk, jj, ii)
                WRITE(unit, '(G0, 1X, G0, 1X, G0)') x(i), y(j), z(k)
                pntxpoli = pntxpoli + 1
            END DO

            ! Write cells
            WRITE(unit,'("CELLS ", I0, 1X, I0)') nblgcells, nblgcells*2

            ! Points
            DO n = 1, nblgcells
                WRITE(unit, '(I0, 1X, I0)') 1, n-1
            END DO

            WRITE(unit, '("CELL_TYPES ", I0)') nblgcells
            DO n = 1, nblgcells
                WRITE(unit, '("1")')
            END DO

            WRITE(unit,'("CELL_DATA ", I0)') nblgcells
            WRITE(unit,'(A)') 'SCALARS fluxsource float 1'
            WRITE(unit,'(A)') 'LOOKUP_TABLE default'
            pntxpolr = 1
            DO n = 1, nblgcells
                !ax1 = xpolr(pntxpolr)
                pntxpolr = pntxpolr + 1
                !ax2 = xpolr(pntxpolr)
                pntxpolr = pntxpolr + 1
                !ay1 = xpolr(pntxpolr)
                pntxpolr = pntxpolr + 1
                !ay2 = xpolr(pntxpolr)
                pntxpolr = pntxpolr + 1
                !az1 = xpolr(pntxpolr)
                pntxpolr = pntxpolr + 1
                !az2 = xpolr(pntxpolr)
                pntxpolr = pntxpolr + 1

                !acoeffstc = xpolr(pntxpolr)
                WRITE(unit, '(G0)') xpolr(pntxpolr)
                pntxpolr = pntxpolr + 1
            END DO
        END IF

        ! Phew...
        CLOSE(unit)
    END SUBROUTINE writestencilsvtk
END MODULE gc_flowstencils_mod
