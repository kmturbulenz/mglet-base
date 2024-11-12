MODULE blockclosetoboundary_mod
    USE core_mod, ONLY: realk, intk

    IMPLICIT NONE (type, external)
    PRIVATE

    PUBLIC :: blockclosetoboundary

CONTAINS
    PURE SUBROUTINE blockclosetoboundary(kk, jj, ii, nfro, nbac, nrgt, nlft, &
            nbot, ntop, skip2, bp)
        ! Subroutine arguments
        INTEGER(intk), INTENT(in) :: kk, jj, ii
        INTEGER(intk), INTENT(in) :: nfro, nbac, nrgt, nlft, nbot, ntop
        LOGICAL, INTENT(in) :: skip2
        REAL(realk), INTENT(inout) :: bp(kk, jj, ii)

        ! Local variables
        INTEGER(intk) :: k, j, i
        INTEGER(intk) :: i2, i3, i4
        INTEGER(intk) :: j2, j3, j4
        INTEGER(intk) :: k2, k3, k4
        REAL(realk) :: flag, bpca, bui

        IF (nfro == 8) THEN
            DO k = 3, kk-2, 2
                DO j = 3, jj-2, 2
                    IF ((bp(k, j, 4) < 0.1 .AND. bp(k, j, 3) > 0.9) &
                            .OR. &
                            (bp(k+1, j, 4) < 0.1 .AND. bp(k+1, j, 3) > 0.9) &
                            .OR. &
                            (bp(k, j+1, 4) < 0.1 .AND. bp(k, j+1, 3) > 0.9) &
                            .OR. &
                            (bp(k+1, j+1, 4) < 0.1 .AND. &
                            bp(k+1, j+1, 3) > 0.9)) THEN
                        bp(k, j, 4) = 0.0
                        bp(k+1, j, 4) = 0.0
                        bp(k, j+1, 4) = 0.0
                        bp(k+1, j+1, 4) = 0.0
                        bp(k, j, 3) = 0.0
                        bp(k+1, j, 3) = 0.0
                        bp(k, j+1, 3) = 0.0
                        bp(k+1, j+1, 3) = 0.0
                    END IF
                END DO
            END DO
        END IF

        IF (nbac == 8) THEN
            DO k = 3, kk-2, 2
                DO j = 3, jj-2, 2
                    IF ((bp(k, j, ii-3) < 0.1 .AND. bp(k, j, ii-2) > 0.9) &
                            .OR. (bp(k+1, j, ii-3) < 0.1 .AND. &
                            bp(k+1, j, ii-2) > 0.9) .OR. &
                            (bp(k, j+1, ii-3) < 0.1 .AND. &
                            bp(k, j+1, ii-2) > 0.9) .OR. &
                            (bp(k+1, j+1, ii-3) < 0.1 .AND. &
                            bp(k+1, j+1, ii-2) > 0.9)) THEN
                        bp(k, j, ii-2) = 0.0
                        bp(k+1, j, ii-2) = 0.0
                        bp(k, j+1, ii-2) = 0.0
                        bp(k+1, j+1, ii-2) = 0.0
                        bp(k, j, ii-3) = 0.0
                        bp(k+1, j, ii-3) = 0.0
                        bp(k, j+1, ii-3) = 0.0
                        bp(k+1, j+1, ii-3) = 0.0
                    END IF
                END DO
            END DO
        END IF

        IF (nrgt == 8) THEN
            DO i = 3, ii-2, 2
                DO k = 3, kk-2, 2
                    IF ((bp(k, 4, i) < 0.1 .AND. bp(k, 3, i) > 0.9) &
                            .OR. &
                            (bp(k+1, 4, i) < 0.1 .AND. bp(k+1, 3, i) > 0.9) &
                            .OR. &
                            (bp(k, 4, i+1) < 0.1 .AND. bp(k, 3, i+1) > 0.9) &
                            .OR. &
                            (bp(k+1, 4, i+1) < 0.1 .AND. &
                            bp(k+1, 3, i+1) > 0.9)) THEN
                        bp(k, 3, i) = 0.0
                        bp(k+1, 3, i) = 0.0
                        bp(k, 3, i+1) = 0.0
                        bp(k+1, 3, i+1) = 0.0
                        bp(k, 4, i) = 0.0
                        bp(k+1, 4, i) = 0.0
                        bp(k, 4, i+1) = 0.0
                        bp(k+1, 4, i+1) = 0.0
                    END IF
                END DO
            END DO
        END IF

        IF (nlft == 8) THEN
            DO i = 3, ii-2, 2
                DO k = 3, kk-2, 2
                    IF ((bp(k, jj-3, i) < 0.1 .AND. bp(k, jj-2, i) > 0.9) &
                            .OR. (bp(k+1, jj-3, i) < 0.1 &
                            .AND. bp(k+1, jj-2, i) > 0.9) .OR. &
                            (bp(k, jj-3, i+1) < 0.1 .AND. &
                            bp(k, jj-2, i+1) > 0.9) .OR. &
                            (bp(k+1, jj-3, i+1) < 0.1 .AND. &
                            bp(k+1, jj-2, i+1) > 0.9)) THEN
                        bp(k, jj-2, i) = 0.0
                        bp(k+1, jj-2, i) = 0.0
                        bp(k, jj-2, i+1) = 0.0
                        bp(k+1, jj-2, i+1) = 0.0
                        bp(k, jj-3, i) = 0.0
                        bp(k+1, jj-3, i) = 0.0
                        bp(k, jj-3, i+1) = 0.0
                        bp(k+1, jj-3, i+1) = 0.0
                    END IF
                END DO
            END DO
        END IF

        IF (nbot == 8) THEN
            DO i = 3, ii-2, 2
                DO j = 3, jj-2, 2
                    IF ((bp(4, j, i) < 0.1 .AND. bp(3, j, i) > 0.9) &
                            .OR. &
                            (bp(4, j+1, i) < 0.1 .AND. bp(3, j+1, i) > 0.9) &
                            .OR. &
                            (bp(4, j, i+1) < 0.1 .AND. bp(3, j, i+1) > 0.9) &
                            .OR. &
                            (bp(4, j+1, i+1) < 0.1 .AND. &
                            bp(3, j+1, i+1) > 0.9)) THEN
                        bp(3, j, i) = 0.0
                        bp(3, j+1, i) = 0.0
                        bp(3, j, i+1) = 0.0
                        bp(3, j+1, i+1) = 0.0
                        bp(4, j, i) = 0.0
                        bp(4, j+1, i) = 0.0
                        bp(4, j, i+1) = 0.0
                        bp(4, j+1, i+1) = 0.0
                    END IF
                END DO
            END DO
        END IF

        IF (ntop == 8) THEN
            DO i = 3, ii-2, 2
                DO j = 3, jj-2, 2
                    IF ((bp(kk-3, j, i) < 0.1 .AND. bp(kk-2, j, i) > 0.9) &
                            .OR. (bp(kk-3, j+1, i) < 0.1 .AND. &
                            bp(kk-2, j+1, i) > 0.9) .OR.&
                            (bp(kk-3, j, i+1) < 0.1 .AND. &
                            bp(kk-2, j, i+1) > 0.9) .OR. &
                            (bp(kk-3, j+1, i+1) < 0.1 .AND. &
                            bp(kk-2, j+1, i+1) > 0.9)) THEN
                        bp(kk-2, j, i) = 0.0
                        bp(kk-2, j+1, i) = 0.0
                        bp(kk-2, j, i+1) = 0.0
                        bp(kk-2, j+1, i+1) = 0.0
                        bp(kk-3, j, i) = 0.0
                        bp(kk-3, j+1, i) = 0.0
                        bp(kk-3, j, i+1) = 0.0
                        bp(kk-3, j+1, i+1) = 0.0
                    END IF
                END DO
            END DO
        END IF

        IF (skip2) RETURN

        ! Fall 2: Interface zu, aber Grobgitter immer BP=1
        IF (nfro == 8) THEN
            i2 = 2
            i3 = 3
            i4 = 4
            DO k = 3, kk-2, 2
                DO j = 3, jj-2, 2
                    bpca = bp(k, j, i2)
                    bui = MIN(1.0_realk, bp(k, j, i2)*bp(k, j, i3) &
                        + bp(k+1, j, i2)*bp(k+1, j, i3) &
                        + bp(k, j+1, i2)*bp(k, j+1, i3) &
                        + bp(k+1, j+1, i2)*bp(k+1, j+1, i3))

                    ! flag = 0 wenn bpca=1 und bui=0
                    flag = 1.0 - bpca*(1.0-bui)

                    bp(k, j, i4) = bp(k, j, i4)*flag
                    bp(k+1, j, i4) = bp(k+1, j, i4)*flag
                    bp(k, j+1, i4) = bp(k, j+1, i4)*flag
                    bp(k+1, j+1, i4) = bp(k+1, j+1, i4)*flag
                END DO
            END DO
        END IF

        IF (nbac == 8) THEN
            i2 = ii-1
            i3 = ii-2
            i4 = ii-3
            DO k = 3, kk-2, 2
                DO j = 3, jj-2, 2
                    bpca = bp(k, j, i2)
                    bui = MIN(1.0_realk, bp(k, j, i2)*bp(k, j, i3) &
                        + bp(k+1, j, i2)*bp(k+1, j, i3) &
                        + bp(k, j+1, i2)*bp(k, j+1, i3) &
                        + bp(k+1, j+1, i2)*bp(k+1, j+1, i3))

                    ! flag = 0 wenn bpca=1 und bui=0
                    flag = 1.0 - bpca*(1.0-bui)

                    bp(k, j, i4) = bp(k, j, i4)*flag
                    bp(k+1, j, i4) = bp(k+1, j, i4)*flag
                    bp(k, j+1, i4) = bp(k, j+1, i4)*flag
                    bp(k+1, j+1, i4) = bp(k+1, j+1, i4)*flag
                END DO
            END DO
        END IF

        IF (nrgt == 8) THEN
            j2 = 2
            j3 = 3
            j4 = 4
            DO k = 3, kk-2, 2
                DO i = 3, ii-2, 2
                    bpca = bp(k, j2, i)
                    bui = MIN(1.0_realk, bp(k, j2, i)*bp(k, j3, i) &
                        + bp(k+1, j2, i)*bp(k+1, j3, i) &
                        + bp(k, j2, i+1)*bp(k, j3, i+1) &
                        + bp(k+1, j2, i+1)*bp(k+1, j3, i+1))

                    ! flag = 0 wenn bpca=1 und bui=0
                    flag = 1.0 - bpca*(1.0-bui)

                    bp(k, j4, i) = bp(k, j4, i)*flag
                    bp(k+1, j4, i) = bp(k+1, j4, i)*flag
                    bp(k, j4, i+1) = bp(k, j4, i+1)*flag
                    bp(k+1, j4, i+1) = bp(k+1, j4, i+1)*flag
                END DO
            END DO
        END IF

        IF (nlft == 8) THEN
            j2 = jj-1
            j3 = jj-2
            j4 = jj-3
            DO k = 3, kk-2, 2
                DO i = 3, ii-2, 2
                    bpca = bp(k, j2, i)
                    bui = MIN(1.0_realk, bp(k, j2, i)*bp(k, j3, i) &
                        + bp(k+1, j2, i)*bp(k+1, j3, i) &
                        + bp(k, j2, i+1)*bp(k, j3, i+1) &
                        + bp(k+1, j2, i+1)*bp(k+1, j3, i+1))

                    ! flag = 0 wenn bpca=1 und bui=0
                    flag = 1.0 - bpca*(1.0-bui)

                    bp(k, j4, i) = bp(k, j4, i)*flag
                    bp(k+1, j4, i) = bp(k+1, j4, i)*flag
                    bp(k, j4, i+1) = bp(k, j4, i+1)*flag
                    bp(k+1, j4, i+1) = bp(k+1, j4, i+1)*flag
                END DO
            END DO
        END IF

        IF (nbot == 8) THEN
            k2 = 2
            k3 = 3
            k4 = 4
            DO j = 3, jj-2, 2
                DO i = 3, ii-2, 2
                    bpca = bp(k2, j, i)
                    bui = MIN(1.0_realk, bp(k2, j, i)*bp(k3, j, i) &
                        + bp(k2, j+1, i)*bp(k3, j+1, i) &
                        + bp(k2, j, i+1)*bp(k3, j, i+1) &
                        + bp(k2, j+1, i+1)*bp(k3, j+1, i+1))

                    ! flag = 0 wenn bpca=1 und bui=0
                    flag = 1.0 - bpca*(1.0-bui)

                    bp(k4, j, i) = bp(k4, j, i)*flag
                    bp(k4, j+1, i) = bp(k4, j+1, i)*flag
                    bp(k4, j, i+1) = bp(k4, j, i+1)*flag
                    bp(k4, j+1, i+1) = bp(k4, j+1, i+1)*flag
                END DO
            END DO
        END IF

        IF (ntop == 8) THEN
            k2 = kk-1
            k3 = kk-2
            k4 = kk-3
            DO j = 3, jj-2, 2
                DO i = 3, ii-2, 2
                    bpca = bp(k2, j, i)
                    bui = MIN(1.0_realk, bp(k2, j, i)*bp(k3, j, i) &
                        + bp(k2, j+1, i)*bp(k3, j+1, i) &
                        + bp(k2, j, i+1)*bp(k3, j, i+1) &
                        + bp(k2, j+1, i+1)*bp(k3, j+1, i+1))

                    ! flag = 0 wenn bpca=1 und bui=0
                    flag = 1.0 - bpca*(1.0-bui)

                    bp(k4, j, i) = bp(k4, j, i)*flag
                    bp(k4, j+1, i) = bp(k4, j+1, i)*flag
                    bp(k4, j, i+1) = bp(k4, j, i+1)*flag
                    bp(k4, j+1, i+1) = bp(k4, j+1, i+1)*flag
                END DO
            END DO
        END IF
    END SUBROUTINE blockclosetoboundary
END MODULE blockclosetoboundary_mod
