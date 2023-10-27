MODULE filling_mod
    USE MPI_f08
    USE core_mod, ONLY: realk, intk, int32, int64, mygridslvl, nmygridslvl, &
        minlevel, maxlevel, errr, connect, mglet_mpi_real, field_t, &
        get_mgdims, get_field, get_mgbasb, get_ip3, myid
    USE ftoc_mod, ONLY: ftoc
    USE ibconst_mod, ONLY: nloopmax, itermax
    USE parent_mod, ONLY: parent

    IMPLICIT NONE (type, external)
    PRIVATE

    PUBLIC :: fillfluid, blockquad_search, &
        blockparentboundary_n, blockparentboundary_p

CONTAINS
    SUBROUTINE fillfluid(fluidpoints, knoten)
        ! Subroutine arguments
        REAL(realk), INTENT(in) :: fluidpoints(:, :)
        TYPE(field_t), INTENT(inout) :: knoten

        ! Local variables
        INTEGER(intk) :: iloop, ilevel, igrid, i, kk, jj, ii, ip3
        INTEGER(intk) :: nfluidpoints, nfound, nfilled
        INTEGER(int64) :: nfilled_tot
        INTEGER(intk), ALLOCATABLE :: ifluidpoints(:, :)
        LOGICAL :: converged
        INTEGER(intk) :: nfro, nbac, nrgt, nlft, nbot, ntop

        TYPE(field_t), POINTER :: x_f, y_f, z_f
        REAL(realk), POINTER, CONTIGUOUS :: x(:), y(:), z(:)

        ! Get gridspacing
        CALL get_field(x_f, "X")
        CALL get_field(y_f, "Y")
        CALL get_field(z_f, "Z")

        ! Initialize parent buffers
        CALL knoten%init_buffers()

        nfluidpoints = SIZE(fluidpoints, 2)
        ALLOCATE(ifluidpoints(3, nfluidpoints))
        ifluidpoints = 0

        converged = .FALSE.
        DO iloop = 1, nloopmax
            DO ilevel = maxlevel, minlevel+1, -1
                CALL ftoc(ilevel, knoten%arr, knoten%arr, 'N')
            END DO

            nfilled_tot = 0
            DO ilevel = minlevel, maxlevel
                CALL parent(ilevel, s1=knoten)

                DO i = 1, nmygridslvl(ilevel)
                    igrid = mygridslvl(i, ilevel)

                    CALL blockparentboundary_n(knoten, igrid)

                    CALL get_mgdims(kk, jj, ii, igrid)
                    CALL get_ip3(ip3, igrid)
                    CALL get_mgbasb(nfro, nbac, nrgt, nlft, nbot, ntop, igrid)
                    CALL x_f%get_ptr(x, igrid)
                    CALL y_f%get_ptr(y, igrid)
                    CALL z_f%get_ptr(z, igrid)

                    CALL blockquad_search(nfluidpoints, nfound, ifluidpoints, &
                        fluidpoints, kk, jj, ii, x, y, z)

                    CALL blockingadd(kk, jj, ii, nfro, nbac, nrgt, nlft, &
                        nbot, ntop, nfound, ifluidpoints, knoten%arr(ip3), &
                        nfilled)

                    nfilled_tot = nfilled_tot + INT(nfilled, int64)
                END DO

                CALL connect(ilevel, 2, s1=knoten, corners=.TRUE.)
            END DO

            CALL MPI_Allreduce(MPI_IN_PLACE, nfilled_tot, 1, MPI_INTEGER8, &
                MPI_SUM, MPI_COMM_WORLD)

            IF (nfilled_tot < 1) THEN
                converged = .TRUE.
                EXIT
            END IF
        END DO

        IF (.NOT. converged) THEN
            WRITE(*,*) "Filling algorithm did not converge: ", iloop
            CALL errr(__FILE__, __LINE__)
        END IF

        IF (myid == 0) THEN
           WRITE(*, '(" Filling finished in ", I0, " iterations")') &
               iloop
        END IF

        DEALLOCATE(ifluidpoints)
    END SUBROUTINE fillfluid


    PURE SUBROUTINE blockquad_search(npts, nfound, i0, x0, kk, jj, ii, x, y, z)
        ! Subroutine arguments
        INTEGER(intk), INTENT(in) :: npts
        INTEGER(intk), INTENT(out) :: nfound
        INTEGER(intk), INTENT(out) :: i0(3, npts)
        REAL(realk), INTENT(in) :: x0(3, npts)
        INTEGER(intk), INTENT(in) :: kk, jj, ii
        REAL(realk), INTENT(in) :: x(ii), y(jj), z(kk)

        ! Local variables
        INTEGER(intk) :: idx, k, j, i, ipos, jpos, kpos
        REAL(realk) :: xpos, ypos, zpos

        nfound = 0
        i0 = 0

        DO idx = 1, npts
            xpos = x0(1, idx)
            ypos = x0(2, idx)
            zpos = x0(3, idx)

            ! Skip points that are completely outside the grid
            IF (xpos < x(1) .OR. xpos > x(ii) &
                    .OR. ypos < y(1) .OR. ypos > y(jj) &
                    .OR. zpos < z(1) .OR. zpos > z(kk)) THEN
                CYCLE
            END IF

            ! Point is inside grid, determine index
            ipos = 1
            DO i = 1, ii
                IF (xpos <= x(i)) THEN
                    ipos = i
                    EXIT
                END IF
            END DO

            jpos = 1
            DO j = 1, jj
                IF (ypos <= y(j)) THEN
                    jpos = j
                    EXIT
                END IF
            END DO

            kpos = 1
            DO k = 1, kk
                IF (zpos <= z(k)) THEN
                    kpos = k
                    EXIT
                END IF
            END DO

            nfound = nfound + 1
            i0(1, nfound) = ipos
            i0(2, nfound) = jpos
            i0(3, nfound) = kpos
        END DO
    END SUBROUTINE blockquad_search


    SUBROUTINE blockingadd(kk, jj, ii, nfro, nbac, &
            nrgt, nlft, nbot, ntop, npts, i0, knoten, nfilled)

        ! Subroutine arguments
        INTEGER(intk), INTENT(in) :: kk, jj, ii
        INTEGER(intk), INTENT(in) :: nfro, nbac, nrgt, nlft, nbot, ntop
        INTEGER(intk), INTENT(in) :: npts
        INTEGER(intk), INTENT(in) :: i0(3, npts)
        REAL(realk), INTENT(inout) :: knoten(kk, jj, ii)
        INTEGER(intk), INTENT(out) :: nfilled

        ! Local variables
        INTEGER(intk) :: idx, k, j, i
        INTEGER(intk) :: ista, isto, jsta, jsto, ksta, ksto
        INTEGER(intk) :: ncount, iloop
        INTEGER(intk), PARAMETER :: itermax = 10000
        LOGICAL :: converged

        nfilled = 0

        ! Initial starting point (seed point)
        DO idx = 1, npts
            i = i0(1, idx)
            j = i0(2, idx)
            k = i0(3, idx)

            IF (knoten(k, j, i) >= -0.5_realk) THEN
                knoten(k, j, i) = 1.0_realk
            END IF
        END DO

        ! Ab hier werden ausgehend vom zuvor bestimmten Knoten
        ! Feld mit Geometrieschnittpunkten die Regionen
        ! innerhalb und ausserhalb des Koerpers
        ! bestimmt. --> Fuellalgorithmus

        ! Schau in Sternform um die nicht-geblockte Zelle
        ! Wenn die betrachtete Zelle nicht geblockt
        ! ist, dann setzte auch die Nachbarzelle
        ! auf nicht-geblockt, wenn diese keine
        ! Geometrie enthaelt

        ! Also fill connect buffers
        ista = 0
        isto = 0
        jsta = 0
        jsto = 0
        ksta = 0
        ksto = 0
        if (nfro == 7) ista = 2
        if (nbac == 7) isto = 2
        if (nrgt == 7) jsta = 2
        if (nlft == 7) jsto = 2
        if (nbot == 7) ksta = 2
        if (ntop == 7) ksto = 2

        converged = .FALSE.
        DO iloop = 1, itermax
            ncount = 0

            DO i = 1+ista, ii-isto
                DO j = 1+jsta, jj-jsto
                    DO k = MAX(2, 1+ksta), kk-ksto
                        IF (knoten(k-1, j, i) >= 0.5_realk) THEN
                            IF (NINT(knoten(k, j, i)) == 0) THEN
                                knoten(k, j, i) = 1.0_realk + knoten(k-1, j, i)
                                ncount = ncount + 1
                            END IF
                        END IF
                    END DO
                END DO
            END DO

            DO i = 1+ista, ii-isto
                DO j = 1+jsta, jj-jsto
                    DO k = MIN(kk-1, kk-ksto), 1+ksta, -1
                        IF (knoten(k+1, j, i) >= 0.5_realk) THEN
                            IF (NINT(knoten(k, j, i)) == 0) THEN
                                knoten(k, j, i) = 1.0_realk + knoten(k+1, j, i)
                                ncount = ncount + 1
                            END IF
                        END IF
                    END DO
                END DO
            END DO

            DO i = 1+ista, ii-isto
                DO j = MAX(2, 1+jsta), jj-jsto
                    DO k = 1+ksta, kk-ksto
                        IF (knoten(k, j-1, i) >= 0.5_realk) THEN
                            IF (NINT(knoten(k, j, i)) == 0) THEN
                                knoten(k, j, i) = 1.0_realk + knoten(k, j-1, i)
                                ncount = ncount + 1
                            END IF
                        END IF
                    END DO
                END DO
            END DO

            DO i = 1+ista, ii-isto
                DO j = MIN(jj-1, jj-jsto), 1+jsta, -1
                    DO k = 1+ksta, kk-ksto
                        IF (knoten(k, j+1, i) >= 0.5_realk) THEN
                            IF (NINT(knoten(k, j, i)) == 0) THEN
                                knoten(k, j, i) = 1.0_realk + knoten(k, j+1, i)
                                ncount = ncount + 1
                            END IF
                        END IF
                    END DO
                END DO
            END DO

            DO i = MAX(1+ista, 2), ii-isto
                DO j = 1+jsta, jj-jsto
                    DO k = 1+ksta, kk-ksto
                        IF (knoten(k, j, i-1) >= 0.5_realk) THEN
                            IF (NINT(knoten(k, j, i)) == 0) THEN
                                knoten(k, j, i) = 1.0_realk + knoten(k, j, i-1)
                                ncount = ncount + 1
                            END IF
                        END IF
                    END DO
                END DO
            END DO

            DO i = MIN(ii-1, ii-isto), 1+ista, -1
                DO j = 1+jsta, jj-jsto
                    DO k = 1+ksta, kk-ksto
                        IF (knoten(k, j, i+1) >= 0.5_realk) THEN
                            IF (NINT(knoten(k, j, i)) == 0) THEN
                                knoten(k, j, i) = 1.0_realk + knoten(k, j, i+1)
                                ncount = ncount + 1
                            END IF
                        END IF
                    END DO
                END DO
            END DO

            nfilled = nfilled + ncount

            ! If no cells were filled during one iteration, leave
            IF (ncount == 0) THEN
                converged = .TRUE.
                EXIT
            END IF
        END DO

        IF (.NOT. converged) THEN
            WRITE(*,*) "Filling algorithm did not converge: ", iloop
            CALL errr(__FILE__, __LINE__)
        END IF
    END SUBROUTINE blockingadd


    SUBROUTINE blockparentboundary_n(knoten_p, igrid)
        ! Subroutine arguments
        TYPE(field_t), INTENT(inout) :: knoten_p
        INTEGER(intk), INTENT(in) :: igrid

        ! Local variables
        INTEGER(intk) :: kk, jj, ii
        INTEGER(intk) :: k, j, i
        INTEGER(intk) :: nfro, nbac, nrgt, nlft, nbot, ntop

        REAL(realk), POINTER, CONTIGUOUS :: knoten(:, :, :)
        REAL(realk), POINTER, CONTIGUOUS :: parbuf(:, :, :)

        CALL knoten_p%get_ptr(knoten, igrid)
        CALL get_mgdims(kk, jj, ii, igrid)
        CALL get_mgbasb(nfro, nbac, nrgt, nlft, nbot, ntop, igrid)

        IF (nfro == 8) THEN
            CALL knoten_p%buffers%get_buffer(parbuf, igrid, 1)
            DO k = 2, kk-2, 2
                DO j = 2, jj-2, 2
                    IF (ABS(parbuf(k, j, 1)) > 0.5_realk) THEN
                        knoten(k, j, 2) = parbuf(k, j, 1)
                    END IF
                END DO
            END DO
        END IF

        IF (nbac == 8) THEN
            CALL knoten_p%buffers%get_buffer(parbuf, igrid, 2)
            DO k = 2, kk-2, 2
                DO j = 2, jj-2, 2
                    IF (ABS(parbuf(k, j, 1)) > 0.5_realk) THEN
                        knoten(k, j, ii-2) = parbuf(k, j, 1)
                    END IF
                END DO
            END DO
        END IF

        IF (nrgt == 8) THEN
            CALL knoten_p%buffers%get_buffer(parbuf, igrid, 3)
            DO k = 2, kk-2, 2
                DO i = 2, ii-2, 2
                    IF (ABS(parbuf(k, i, 1)) > 0.5_realk) THEN
                        knoten(k, 2, i) = parbuf(k, i, 1)
                    END IF
                END DO
            END DO
        END IF

        IF (nlft == 8) THEN
            CALL knoten_p%buffers%get_buffer(parbuf, igrid, 4)
            DO k = 2, kk-2, 2
                DO i = 2, ii-2, 2
                    IF (ABS(parbuf(k, i, 1)) > 0.5_realk) THEN
                        knoten(k, jj-2, i) = parbuf(k, i, 1)
                    END IF
                END DO
            END DO
        END IF

        IF (nbot == 8) THEN
            CALL knoten_p%buffers%get_buffer(parbuf, igrid, 5)
            DO j = 2, jj-2, 2
                DO i = 2, ii-2, 2
                    IF (ABS(parbuf(j, i, 1)) > 0.5_realk) THEN
                        knoten(2, j, i) = parbuf(j, i, 1)
                    END IF
                END DO
            END DO
        END IF

        IF (ntop == 8) THEN
            CALL knoten_p%buffers%get_buffer(parbuf, igrid, 6)
            DO j = 2, jj-2, 2
                DO i = 2, ii-2, 2
                    IF (ABS(parbuf(j, i, 1)) > 0.5_realk) THEN
                        knoten(kk-2, j, i) = parbuf(j, i, 1)
                    END IF
                END DO
            END DO
        END IF
    END SUBROUTINE blockparentboundary_n


    SUBROUTINE blockparentboundary_p(bp_p, igrid)
        ! Subroutine arguments
        TYPE(field_t), INTENT(inout) :: bp_p
        INTEGER(intk), INTENT(in) :: igrid

        ! Local variables
        INTEGER(intk) :: kk, jj, ii
        INTEGER(intk) :: ibbuf1, ibbuf2, ibbuf3, ibbuf4
        INTEGER(intk) :: k, j, i
        INTEGER(intk) :: nfro, nbac, nrgt, nlft, nbot, ntop

        REAL(realk), POINTER, CONTIGUOUS :: bp(:, :, :)
        REAL(realk), POINTER, CONTIGUOUS :: parbuf(:, :, :)

        CALL bp_p%get_ptr(bp, igrid)
        CALL get_mgdims(kk, jj, ii, igrid)
        CALL get_mgbasb(nfro, nbac, nrgt, nlft, nbot, ntop, igrid)

        ibbuf1 = 1
        ibbuf2 = 2
        ibbuf3 = 1
        ibbuf4 = 0

        IF (nfro == 8) THEN
            CALL bp_p%buffers%get_buffer(parbuf, igrid, 1)
            DO k = 1, kk
                DO j = 1, jj
                    DO i = ibbuf1, ibbuf2
                        bp(k, j, i) = parbuf(k, j, 1)
                    END DO
                END DO
            END DO
        END IF

        IF (nbac == 8) THEN
            CALL bp_p%buffers%get_buffer(parbuf, igrid, 2)
            DO k = 1, kk
                DO j = 1, jj
                    DO i = ii-ibbuf3, ii-ibbuf4
                        bp(k, j, i) = parbuf(k, j, 1)
                    END DO
                END DO
            END DO
        END IF

        IF (nrgt == 8) THEN
            CALL bp_p%buffers%get_buffer(parbuf, igrid, 3)
            DO k = 1, kk
                DO j = ibbuf1, ibbuf2
                    DO i = 1, ii
                        bp(k, j, i) = parbuf(k, i, 1)
                    END DO
                END DO
            END DO
        END IF

        IF (nlft == 8) THEN
            CALL bp_p%buffers%get_buffer(parbuf, igrid, 4)
            DO k = 1, kk
                DO j = jj-ibbuf3, jj-ibbuf4
                    DO i = 1, ii
                        bp(k, j, i) = parbuf(k, i, 1)
                    END DO
                END DO
            END DO
        END IF

        IF (nbot == 8) THEN
            CALL bp_p%buffers%get_buffer(parbuf, igrid, 5)
            DO k = ibbuf1, ibbuf2
                DO j = 1, jj
                    DO i = 1, ii
                        bp(k, j, i) = parbuf(j, i, 1)
                    END DO
                END DO
            END DO
        END IF

        IF (ntop == 8) THEN
            CALL bp_p%buffers%get_buffer(parbuf, igrid, 6)
            DO k = kk-ibbuf3, kk-ibbuf4
                DO j = 1, jj
                    DO i = 1, ii
                        bp(k, j, i) = parbuf(j, i, 1)
                    END DO
                END DO
            END DO
        END IF
    END SUBROUTINE blockparentboundary_p

END MODULE filling_mod
