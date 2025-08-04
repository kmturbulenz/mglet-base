MODULE gc_totwasser_mod
    USE MPI_f08
    USE core_mod, ONLY: realk, intk, int64, mygridslvl, nmygridslvl, &
        nmygrids, mygrids, minlevel, maxlevel, errr, connect, field_t, &
        get_mgdims, get_field, get_mgbasb, get_ip3, myid
    USE ftoc_mod, ONLY: ftoc
    USE ibconst_mod, ONLY: nloopmax
    USE parent_mod, ONLY: parent
    USE filling_mod, ONLY: blockparentboundary_p, blockquad_search

    IMPLICIT NONE (type, external)
    PRIVATE

    PUBLIC :: totwasser

CONTAINS
    SUBROUTINE totwasser(fluidpoints, bp)
        ! Subroutine arguments
        REAL(realk), INTENT(in) :: fluidpoints(:, :)
        TYPE(field_t), INTENT(inout) :: bp

        ! Local variables
        INTEGER(intk) :: iloop, ilevel, igrid, i, kk, jj, ii, ip3
        INTEGER(intk) :: nfluidpoints, nfound
        INTEGER(int64) :: nfilled_tot, nfilled_old
        INTEGER(int64) :: cellcounter(3)
        INTEGER(intk), ALLOCATABLE :: ifluidpoints(:, :)
        LOGICAL :: converged
        INTEGER(intk) :: nfro, nbac, nrgt, nlft, nbot, ntop

        TYPE(field_t), POINTER :: x_f, y_f, z_f, finecell_f
        REAL(realk), POINTER, CONTIGUOUS :: x(:), y(:), z(:)
        REAL(realk), POINTER, CONTIGUOUS :: finecell(:)

        ! Get gridspacing
        CALL get_field(x_f, "X")
        CALL get_field(y_f, "Y")
        CALL get_field(z_f, "Z")
        CALL get_field(finecell_f, "FINECELL")

        nfluidpoints = SIZE(fluidpoints, 2)
        ALLOCATE(ifluidpoints(3, nfluidpoints))
        ifluidpoints = 0

        converged = .FALSE.
        nfilled_old = 0
        nfilled_tot = 0
        DO iloop = 1, nloopmax
            DO ilevel = maxlevel, minlevel+1, -1
                CALL ftoc(ilevel, bp%arr, bp%arr, 'F')
            END DO

            DO ilevel = minlevel, maxlevel
                CALL parent(ilevel, s1=bp)

                DO i = 1, nmygridslvl(ilevel)
                    igrid = mygridslvl(i, ilevel)

                    CALL blockparentboundary_p(bp, igrid)

                    CALL get_mgdims(kk, jj, ii, igrid)
                    CALL get_ip3(ip3, igrid)
                    CALL get_mgbasb(nfro, nbac, nrgt, nlft, nbot, ntop, igrid)
                    CALL x_f%get_ptr(x, igrid)
                    CALL y_f%get_ptr(y, igrid)
                    CALL z_f%get_ptr(z, igrid)
                    CALL finecell_f%get_ptr(finecell, igrid)

                    CALL blockquad_search(nfluidpoints, nfound, ifluidpoints, &
                        fluidpoints, kk, jj, ii, x, y, z)

                    CALL totwasser_grid(kk, jj, ii, nfro, nbac, nrgt, nlft, &
                        nbot, ntop, nfound, ifluidpoints, finecell, &
                        bp%arr(ip3))
                END DO

                CALL connect(ilevel, 2, s1=bp, corners=.TRUE.)
            END DO

            ! This counts the overall number of filled cells in the domain
            ! after each iterattion. If the number of filled cells does not
            ! change anymore, the totwasser is converged.
            nfilled_old = nfilled_tot
            CALL count_filled(nfilled_tot, bp)

            IF (nfilled_old == nfilled_tot) THEN
                converged = .TRUE.
                EXIT
            END IF
        END DO

        IF (.NOT. converged) THEN
            WRITE(*, *) "Totwasser algorithm did not converge: ", iloop
            CALL errr(__FILE__, __LINE__)
        END IF

        IF (myid == 0) THEN
            WRITE(*, '(" Totwasser finished in ", I0, " iterations")') iloop
        END IF

        DEALLOCATE(ifluidpoints)

        ! Finish BP field and count cells
        cellcounter = 0
        DO ilevel = minlevel, maxlevel
            DO i = 1, nmygridslvl(ilevel)
                igrid = mygridslvl(i, ilevel)

                CALL get_mgdims(kk, jj, ii, igrid)
                CALL get_ip3(ip3, igrid)
                CALL totwasser_finish(kk, jj, ii, cellcounter(1), &
                    cellcounter(2), cellcounter(3), bp%arr(ip3))
            END DO
        END DO
        CALL MPI_Allreduce(MPI_IN_PLACE, cellcounter, 3, MPI_INTEGER8, &
            MPI_SUM, MPI_COMM_WORLD)
        IF (myid == 0) THEN
            WRITE(*, '(" Found ", I0, " blocked, ", I0, " totwasser and ", ' &
                //'I0, " open cells")') cellcounter(1), cellcounter(2), &
                cellcounter(3)
        END IF

        ! Previously directly in blockbp - this is last place where these
        ! fields are touched
        DO ilevel = minlevel, maxlevel
            CALL connect(ilevel, 2, s1=bp, corners=.TRUE.)
        END DO
    END SUBROUTINE totwasser


    SUBROUTINE totwasser_grid(kk, jj, ii, nfro, nbac, nrgt, nlft, nbot, ntop, &
            npts, i0, finecell, bp)
        ! Das BP-Feld hat zu diesem Zeitpunkt entweder 1 (Fluid) oder 0 (Koerper)
        ! Zum Finden der Totwassergebiete wird das BP-Feld ausgehend von der
        ! Painting-Area um 1 erhöht. Abschliessend werden alle BP=1 Werte auf 0
        ! gesetzt und alle BP>1 auf 1

        ! Subroutine arguments
        INTEGER(intk), INTENT(in) :: kk, jj, ii
        INTEGER(intk), INTENT(in) :: nfro, nbac, nrgt, nlft, nbot, ntop
        INTEGER(intk), INTENT(in) :: npts
        INTEGER(intk), INTENT(in) :: i0(3, npts)
        REAL(realk), INTENT(in) :: finecell(kk, jj, ii)
        REAL(realk), INTENT(inout) :: bp(kk, jj, ii)

        ! Local variables
        INTEGER(intk) :: idx, k, j, i
        INTEGER(intk) :: ista, isto, jsta, jsto, ksta, ksto
        INTEGER(intk) :: ncount, iloop
        INTEGER(intk), PARAMETER :: itermax = 10000
        LOGICAL :: converged

        ! Initial starting point (seed point)
        DO idx = 1, npts
            i = i0(1, idx)
            j = i0(2, idx)
            k = i0(3, idx)

            IF (bp(k, j, i) > 0.0) THEN
                bp(k, j, i) = bp(k, j, i) + 1.0
            END IF
        END DO

        ista = 0
        isto = 0
        jsta = 0
        jsto = 0
        ksta = 0
        ksto = 0
        if (nfro == 7 .OR. nfro == 8) ista = 2
        if (nbac == 7 .OR. nbac == 8) isto = 2
        if (nrgt == 7 .OR. nrgt == 8) jsta = 2
        if (nlft == 7 .OR. nlft == 8) jsto = 2
        if (nbot == 7 .OR. nbot == 8) ksta = 2
        if (ntop == 7 .OR. ntop == 8) ksto = 2

        converged = .FALSE.
        DO iloop = 1, itermax
            ncount = 0

            DO i = 1+ista, ii-isto
                DO j = 1+jsta, jj-jsto
                    DO k = MAX(2, 1+ksta), kk-ksto
                        ! Using NINT is not performant here.
                        ! IF (NINT(bp(k-1, j, i)) >= 2) THEN
                        IF (bp(k-1, j, i) > 1.5_realk) THEN
                            ! Various ways to check if the bp value is 1.0
                            ! without performing a flaoting point equality
                            ! check. NINT is the slowest, the two others are
                            ! more or less equivalent.
                            ! IF (NINT(bp(k, j, i)) == 1) THEN
                            ! IF (ABS(bp(k, j, i) - 1.0) < 0.001) THEN
                            IF (bp(k, j, i) > 0.9 .AND. bp(k, j, i) < 1.1) THEN
                                bp(k, j, i) = 1.0 &
                                    + bp(k-1, j, i)*finecell(k, j, i)
                                ncount = ncount + NINT(finecell(k, j, i))
                            END IF
                        END IF
                    END DO
                END DO
            END DO

            DO i = 1+ista, ii-isto
                DO j = 1+jsta, jj-jsto
                    DO k = MIN(kk-1, kk-ksto), 1+ksta, -1
                        IF (bp(k+1, j, i) > 1.5_realk) THEN
                            IF (bp(k, j, i) > 0.9 .AND. bp(k, j, i) < 1.1) THEN
                                bp(k, j, i) = 1.0 &
                                    + bp(k+1, j, i)*finecell(k, j, i)
                                ncount = ncount + NINT(finecell(k, j, i))
                            END IF
                        END IF
                    END DO
                END DO
            END DO

            DO i = 1+ista, ii-isto
                DO j = MAX(2, 1+jsta), jj-jsto
                    DO k = 1+ksta, kk-ksto
                        IF (bp(k, j-1, i) > 1.5_realk) THEN
                            IF (bp(k, j, i) > 0.9 .AND. bp(k, j, i) < 1.1) THEN
                                bp(k, j, i) = 1.0 &
                                    + bp(k, j-1, i)*finecell(k, j, i)
                                ncount = ncount + NINT(finecell(k, j, i))
                            END IF
                        END IF
                    END DO
                END DO
            END DO

            DO i = 1+ista, ii-isto
                DO j = MIN(jj-1, jj-jsto), 1+jsta, -1
                    DO k = 1+ksta, kk-ksto
                        IF (bp(k, j+1, i) > 1.5_realk) THEN
                            IF (bp(k, j, i) > 0.9 .AND. bp(k, j, i) < 1.1) THEN
                                bp(k, j, i) = 1.0 &
                                    + bp(k, j+1, i)*finecell(k, j, i)
                                ncount = ncount + NINT(finecell(k, j, i))
                            END IF
                        END IF
                    END DO
                END DO
            END DO

            DO i = MAX(1+ista, 2), ii-isto
                DO j = 1+jsta, jj-jsto
                    DO k = 1+ksta, kk-ksto
                        IF (bp(k, j, i-1) > 1.5_realk) THEN
                            IF (bp(k, j, i) > 0.9 .AND. bp(k, j, i) < 1.1) THEN
                                bp(k, j, i) = 1.0 &
                                    + bp(k, j, i-1)*finecell(k, j, i)
                                ncount = ncount + NINT(finecell(k, j, i))
                            END IF
                        END IF
                    END DO
                END DO
            END DO

            DO i = MIN(ii-1, ii-isto), 1+ista, -1
                DO j = 1+jsta, jj-jsto
                    DO k = 1+ksta, kk-ksto
                        IF (bp(k, j, i+1) > 1.5_realk) THEN
                            IF (bp(k, j, i) > 0.9 .AND. bp(k, j, i) < 1.1) THEN
                                bp(k, j, i) = 1.0 &
                                    + bp(k, j, i+1)*finecell(k, j, i)
                                ncount = ncount + NINT(finecell(k, j, i))
                            END IF
                        END IF
                    END DO
                END DO
            END DO

            ! If no cells were filled during one iteration, leave
            IF (ncount == 0) THEN
                converged = .TRUE.
                EXIT
            END IF
        END DO

        IF (.NOT. converged) THEN
            WRITE(*, *) "Totwasser algorithm did not converge: ", iloop
            CALL errr(__FILE__, __LINE__)
        END IF
    END SUBROUTINE totwasser_grid


    SUBROUTINE totwasser_finish(kk, jj, ii, count_block, count_tot, &
            count_fil, bp)
        ! Subroutine arguments
        INTEGER(intk), INTENT(in) :: kk, jj, ii
        INTEGER(int64), INTENT(inout) :: count_block, count_tot, count_fil
        REAL(realk), INTENT(inout) :: bp(kk, jj, ii)

        ! Local variables
        INTEGER(intk) :: k, j, i

        ! Finish BP field and count blocked, totwasser and open cells
        DO i = 1, ii
            DO j = 1, jj
                DO k = 1, kk
                    IF (bp(k, j, i) <= 0.5) THEN
                        count_block = count_block + 1
                    ELSE IF (bp(k, j, i) <= 1.5) THEN
                        bp(k, j, i) = 0.0
                        count_tot = count_tot + 1
                    ELSE
                        bp(k, j, i) = 1.0
                        count_fil = count_fil + 1
                    END IF
                END DO
            END DO
        END DO
    END SUBROUTINE totwasser_finish


    SUBROUTINE count_filled(nall, bp_f)
        ! To make it easier to make 64-bit integer literals
        USE precision_mod, ONLY: i8 => int64

        ! Subroutine arguments
        INTEGER(i8), INTENT(out) :: nall
        TYPE(field_t), INTENT(inout) :: bp_f

        ! Local variables
        INTEGER(i8) :: nsingle
        INTEGER(intk) :: igr, igrid, bpi
        INTEGER(intk) :: k, j, i
        INTEGER(intk) :: kk, jj, ii
        REAL(realk), POINTER, CONTIGUOUS :: bp(:, :, :)

        nsingle = 0

        DO igr = 1, nmygrids
            igrid = mygrids(igr)
            CALL get_mgdims(kk, jj, ii, igrid)
            CALL bp_f%get_ptr(bp, igrid)

            DO i = 3, ii-2
                DO j = 3, jj-2
                    DO k = 3, kk-2
                        bpi = NINT(bp(k, j, i))
                        nsingle = nsingle + MIN(MAX(bpi - 1_i8, 0_i8), 1_i8)
                    END DO
                END DO
            END DO
        END DO

        CALL MPI_Allreduce(nsingle, nall, 1, MPI_INTEGER8, MPI_SUM, &
            MPI_COMM_WORLD)
    END SUBROUTINE count_filled

END MODULE gc_totwasser_mod
