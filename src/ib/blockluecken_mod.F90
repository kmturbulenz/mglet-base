MODULE blockluecken_mod
    USE core_mod, ONLY: realk, intk, mygrids, nmygrids, &
        get_mgdims, get_mgbasb, errr, field_t

    IMPLICIT NONE (type, external)
    PRIVATE

    PUBLIC :: blockluecken, blockluecken_grid

CONTAINS
    SUBROUTINE blockluecken(bp_f, use_bc)
        ! Subroutine arguments
        TYPE(field_t), INTENT(inout) :: bp_f
        LOGICAL, INTENT(in) :: use_bc

        ! Local variables
        INTEGER(intk) :: i, igrid, kk, jj, ii
        INTEGER(intk) :: nfro, nbac, nrgt, nlft, nbot, ntop
        REAL(realk), POINTER, CONTIGUOUS :: bp(:, :, :)

        DO i = 1, nmygrids
            igrid = mygrids(i)
            CALL get_mgdims(kk, jj, ii, igrid)

            IF (use_bc) THEN
                CALL get_mgbasb(nfro, nbac, nrgt, nlft, nbot, ntop, igrid)
            ELSE
                nfro = 0
                nbac = 0
                nrgt = 0
                nlft = 0
                nbot = 0
                ntop = 0
            END IF
            CALL bp_f%get_ptr(bp, igrid)
            CALL blockluecken_grid(kk, jj, ii, nfro, nbac, nrgt, nlft, &
                nbot, ntop, bp)
        END DO
    END SUBROUTINE blockluecken


    SUBROUTINE blockluecken_grid(kk, jj, ii, nfro, nbac, nrgt, nlft, &
            nbot, ntop, bp)
        ! Subroutine arguments
        INTEGER(intk), INTENT(in) :: kk, jj, ii
        INTEGER(intk), INTENT(in) :: nfro, nbac, nrgt, nlft, nbot, ntop
        REAL(realk), INTENT(inout) :: bp(kk, jj, ii)

        ! Local variables
        INTEGER(intk) :: i, j, k
        INTEGER(intk) :: iloop, nmodb
        INTEGER(intk), PARAMETER :: nloopmax = 320
        INTEGER(intk) :: nxm, nxp, nym, nyp, nzm, nzp

        nxm = -1
        nxp = -1
        nym = -1
        nyp = -1
        nzm = -1
        nzp = -1

        ! Spezial bei CON und NIX
        IF (nfro == 7 .OR. nfro == 99) nxm = 0
        IF (nbac == 7 .OR. nbac == 99) nxp = 0
        IF (nrgt == 7 .OR. nrgt == 99) nym = 0
        IF (nlft == 7 .OR. nlft == 99) nyp = 0
        IF (nbot == 7 .OR. nbot == 99) nzm = 0
        IF (ntop == 7 .OR. ntop == 99) nzp = 0

        ! NIX boundary: no open cells allowed
        IF (nfro == 99) THEN
            DO i = 1, 2
                DO j = 3, jj-2
                    DO k = 3, kk-2
                        bp(k, j, i) = 0.0
                    END DO
                END DO
            END DO
        END IF

        IF (nbac == 99) THEN
            DO i = ii-1, ii
                DO j = 3, jj-2
                    DO k = 3, kk-2
                        bp(k, j, i) = 0.0
                    END DO
                END DO
            END DO
        END IF

        IF (nrgt == 99) THEN
            DO i = 3, ii-2
                DO j = 1, 2
                    DO k = 3, kk-2
                        bp(k, j, i)= 0.0
                    END DO
                END DO
            END DO
        END IF

        IF (nlft == 99) THEN
            DO i = 3, ii-2
                DO j = jj-1, jj
                    DO k = 3, kk-2
                        bp(k, j, i)= 0.0
                    END DO
                END DO
            END DO
        END IF

        IF (nbot == 99) THEN
            DO i = 3, ii-2
                DO j = 3, jj-2
                    DO k = 1, 2
                        bp(k, j, i)= 0.0
                    END DO
                END DO
            END DO
        END IF

        IF (ntop == 99) THEN
            DO i = 3, ii-2
                DO j = 3, jj-2
                    DO k = kk-1, kk
                        bp(k, j, i) = 0.0
                    END DO
                END DO
            END DO
        END IF

        DO iloop = 1, nloopmax
            nmodb = 0
            ! Z
            DO i = 3, ii-2
                DO j = 3, jj-2
                    DO k = 3, kk-2+nzp
                        IF (NINT(bp(k, j, i)) == 1) THEN
                            IF (bp(k+2, j, i) <= 0.5 .AND. &
                                    bp(k-1, j, i) <= 0.5) THEN
                                bp(k, j, i) = 0.0
                                nmodb = nmodb + 1
                            END IF
                        END IF
                    END DO
                END DO
            END DO

            DO i = 3, ii-2
                DO j = 3, jj-2
                    DO k = 3-nzm, kk-2
                        IF (NINT(bp(k, j, i)) == 1) THEN
                            IF (bp(k+1, j, i) <= 0.5 .AND. &
                                    bp(k-2, j, i) <= 0.5) THEN
                                bp(k, j, i) = 0.0
                                nmodb = nmodb + 1
                            END IF
                        END IF
                    END DO
                END DO
            END DO

            ! Y
            DO i = 3, ii-2
                DO j = 3, jj-2+nyp
                    DO k = 3, kk-2
                        IF (NINT(bp(k, j, i)) == 1) THEN
                            IF (bp(k, j+2, i) <= 0.5 .AND. &
                                    bp(k, j-1, i) <= 0.5) THEN
                                bp(k, j, i) = 0.0
                                nmodb = nmodb + 1
                            END IF
                        END IF
                    END DO
                END DO
            END DO

            DO i = 3, ii-2
                DO j = 3-nym, jj-2
                    DO k = 3, kk-2
                        IF (NINT(bp(k, j, i)) == 1) THEN
                            IF (bp(k, j+1, i) <= 0.5 .AND. &
                                    bp(k, j-2, i) <= 0.5) THEN
                                bp(k, j, i) = 0.0
                                nmodb = nmodb + 1
                            END IF
                        END IF
                    END DO
                END DO
            END DO

            ! X
            DO i = 3, ii-2+nxp
                DO j = 3, jj-2
                    DO k = 3, kk-2
                        IF (NINT(bp(k, j, i)) == 1) THEN
                            IF (bp(k, j, i+2) <= 0.5 .AND. &
                                    bp(k, j, i-1) <= 0.5) THEN
                                bp(k, j, i) = 0.0
                                nmodb = nmodb + 1
                            END IF
                        END IF
                    END DO
                END DO
            END DO

            DO i = 3-nxm, ii-2
                DO j = 3, jj-2
                    DO k = 3, kk-2
                        IF (NINT(bp(k, j, i)) == 1) THEN
                            IF (bp(k, j, i+1) <= 0.5 .AND. &
                                    bp(k, j, i-2) <= 0.5) THEN
                                bp(k, j, i) = 0.0
                                nmodb = nmodb + 1
                            END IF
                        END IF
                    END DO
                END DO
            END DO

            ! einser Luecken
            ! Z
            DO i = 3, ii-2
                DO j = 3, jj-2
                    DO k = 3, kk-2
                        IF (NINT(bp(k, j, i)) == 1) THEN
                            IF (bp(k+1, j, i) <= 0.5 .AND. &
                                    bp(k-1, j, i) <= 0.5) THEN
                                bp(k, j, i) = 0.0
                                nmodb = nmodb + 1
                            END IF
                        END IF
                    END DO
                END DO
            END DO

            ! Y
            DO i = 3, ii-2
                DO j = 3, jj-2
                    DO k = 3, kk-2
                        IF (NINT(bp(k, j, i)) == 1) THEN
                            IF (bp(k, j+1, i) <= 0.5 .AND. &
                                    bp(k, j-1, i) <= 0.5) THEN
                                bp(k, j, i) = 0.0
                                nmodb = nmodb + 1
                            END IF
                        END IF
                    END DO
                END DO
            END DO

            ! X
            DO i = 3, ii-2
                DO j = 3, jj-2
                    DO k = 3, kk-2
                        IF (NINT(bp(k, j, i)) == 1) THEN
                            IF (bp(k, j, i+1) <= 0.5 .AND. &
                                    bp(k, j, i-1) <= 0.5) THEN
                                bp(k, j, i) = 0.0
                                nmodb = nmodb + 1
                            END IF
                        END IF
                    END DO
                END DO
            END DO

            ! special treatment of domain boundaries with fixed normal velocity
            ! FIX, NOS, SLI
            IF (nfro == 2 .OR. nfro == 5 .OR. nfro == 6) THEN
                DO i = MIN(4, ii-2), 2, -1
                    DO j = 3, jj-2
                        DO k = 3, kk-2
                            IF (NINT(bp(k, j, i)) == 1) THEN
                                IF (bp(k, j, i+1) <= 0.5) THEN
                                    bp(k, j, i) = 0.0
                                    nmodb = nmodb + 1
                                END IF
                            END IF
                        END DO
                    END DO
                END DO
            END IF

            IF (nbac == 2 .OR. nbac == 5 .OR. nbac == 6) THEN
                DO i = MAX(ii-3, 3), ii-1
                    DO j = 3, jj-2
                        DO k = 3, kk-2
                            IF (NINT(bp(k, j, i)) == 1) THEN
                                IF (bp(k, j, i-1) <= 0.5) THEN
                                    bp(k, j, i) = 0.0
                                    nmodb = nmodb + 1
                                END IF
                            END IF
                        END DO
                    END DO
                END DO
            END IF

            IF (nrgt == 2 .OR. nrgt == 5 .OR. nrgt == 6) THEN
                DO i = 3, ii-2
                    DO j = MIN(4, jj-2), 2, -1
                        DO k = 3, kk-2
                            IF (NINT(bp(k, j, i)) == 1) THEN
                                IF (bp(k, j+1, i) <= 0.5) THEN
                                    bp(k, j, i) = 0.0
                                    nmodb = nmodb + 1
                                END IF
                            END IF
                        END DO
                    END DO
                END DO
            END IF

            IF (nlft == 2 .OR. nlft == 5 .OR. nlft == 6) THEN
                DO i = 3, ii-2
                    DO j = MAX(jj-3, 3), jj-1
                        DO k = 3, kk-2
                            IF (NINT(bp(k, j, i)) == 1) THEN
                                IF (bp(k, j-1, i) <= 0.5) THEN
                                    bp(k, j, i) = 0.0
                                    nmodb = nmodb + 1
                                END IF
                            END IF
                        END DO
                    END DO
                END DO
            END IF

            IF (nbot == 2 .OR. nbot == 5 .OR. nbot == 6) THEN
                DO i = 3, ii-2
                    DO j = 3, jj-2
                        DO k = MIN(4, kk-2), 2, -1
                            IF (NINT(bp(k, j, i)) == 1) THEN
                                IF (bp(k+1, j, i) <= 0.5) THEN
                                    bp(k, j, i) = 0.0
                                    nmodb = nmodb + 1
                                END IF
                            END IF
                        END DO
                    END DO
                END DO
            END IF

            IF (ntop == 2 .OR. ntop == 5 .OR. ntop == 6) THEN
                DO i = 3, ii-2
                    DO j = 3, jj-2
                        DO k = MAX(kk-3, 3), kk-1
                            IF (NINT(bp(k, j, i)) == 1) THEN
                                IF (bp(k-1, j, i) <= 0.5) THEN
                                    bp(k, j, i) = 0.0
                                    nmodb = nmodb + 1
                                END IF
                            END IF
                        END DO
                    END DO
                END DO
            END IF

            IF (nmodb == 0) EXIT

            IF (iloop == nloopmax) THEN
                CALL errr(__FILE__, __LINE__)
            END IF
        END DO
    END SUBROUTINE blockluecken_grid
END MODULE blockluecken_mod
