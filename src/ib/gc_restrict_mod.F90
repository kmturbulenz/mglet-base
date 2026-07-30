MODULE gc_restrict_mod
    USE core_mod

    IMPLICIT NONE(type, external)
    PRIVATE

    PUBLIC :: restrict_gc_p_t, restrict_gc_r, restrict_gc_n, restrict_gc_e, &
        restrict_gc_f, restrict_gc_i
CONTAINS
    SUBROUTINE restrict_gc_p_t(kk, jj, ii, ff, ddx, ddy, ddz, b, &
            messagelength, sbuf, istart, istop, jstart, jstop, kstart, kstop)
        !$omp declare target
        ! Subroutine arguments
        INTEGER(intk), INTENT(in) :: kk, jj, ii
        REAL(realk), INTENT(in) :: ff(kk, jj, ii)
        REAL(realk), INTENT(in) :: ddx(ii), ddy(jj), ddz(kk)
        REAL(realk), INTENT(in) :: b(kk, jj, ii)
        INTEGER(int32), INTENT(in) :: messagelength
        REAL(realk), INTENT(inout) :: sbuf(messagelength)
        INTEGER(intk), INTENT(in) :: istart, istop, jstart, jstop, kstart, kstop

        ! Local variables
        INTEGER(intk) :: i, j, k, icount
        INTEGER(intk) :: nk, ink

        REAL(realk) :: vol1(kk), vol2(kk), vol3(kk), vol4(kk)
        REAL(realk) :: pv1(kk), pv2(kk), pv3(kk), pv4(kk)
        REAL(realk) :: sum_pv, sum_v

        ! Number of k iterations
        ! (not to be confused with kk - these are not the same!!!)
        nk = (kstop-kstart)/2+1

        icount = 0
        DO i = istart, istop, 2
            DO j = jstart, jstop, 2
                ! Vomume
                DO k = kstart, kstop+1
                    vol1(k) = b(k, j, i)*ddz(k)*ddy(j)*ddx(i)
                END DO
                DO k = kstart, kstop+1
                    vol2(k) = b(k, j, i+1)*ddz(k)*ddy(j)*ddx(i+1)
                END DO
                DO k = kstart, kstop+1
                    vol3(k) = b(k, j+1, i)*ddz(k)*ddy(j+1)*ddx(i)
                END DO
                DO k = kstart, kstop+1
                    vol4(k) = b(k, j+1, i+1)*ddz(k)*ddy(j+1)*ddx(i+1)
                END DO

                ! Pressure times volume
                DO k = kstart, kstop+1
                    pv1(k) = ff(k, j, i)*vol1(k)
                END DO
                DO k = kstart, kstop+1
                    pv2(k) = ff(k, j, i+1)*vol2(k)
                END DO
                DO k = kstart, kstop+1
                    pv3(k) = ff(k, j+1, i)*vol3(k)
                END DO
                DO k = kstart, kstop+1
                    pv4(k) = ff(k, j+1, i+1)*vol4(k)
                END DO

                ! Sum up and divide
                DO ink = 1, nk
                    k = kstart + 2*(ink-1)

                    sum_pv = pv1(k) + pv2(k) + pv3(k) + pv4(k) &
                        + pv1(k+1) + pv2(k+1) + pv3(k+1) + pv4(k+1)

                    sum_v = vol1(k) + vol2(k) + vol3(k) + vol4(k) &
                        + vol1(k+1) + vol2(k+1) + vol3(k+1) + vol4(k+1)

                    IF (sum_v > 0.0) THEN
                        sbuf(icount + ink) = sum_pv/sum_v
                    ELSE
                        sbuf(icount + ink) = 0.0
                    END IF
                END DO

                ! Increment counter
                icount = icount + nk

                ! Legacy code - keep for reference
                ! DO k = kstart, kstop, 2
                !     sum_pv = ff(k, j, i)*b(k, j, i)*ddz(k)*ddy(j)*ddx(i) &
                !         + ff(k, j, i+1)*b(k, j, i+1)*ddz(k)*ddy(j)*ddx(i+1) &
                !         + ff(k, j+1, i)*b(k, j+1, i)*ddz(k)*ddy(j+1)*ddx(i) &
                !         + ff(k, j+1, i+1)*b(k, j+1, i+1)*ddz(k)*ddy(j+1)*ddx(i+1) &
                !         + ff(k+1, j, i)*b(k+1, j, i)*ddz(k+1)*ddy(j)*ddx(i) &
                !         + ff(k+1, j, i+1)*b(k+1, j, i+1)*ddz(k+1)*ddy(j)*ddx(i+1) &
                !         + ff(k+1, j+1, i)*b(k+1, j+1, i)*ddz(k+1)*ddy(j+1)*ddx(i) &
                !         + ff(k+1, j+1, i+1)*b(k+1, j+1, i+1)*ddz(k+1)*ddy(j+1)*ddx(i+1)

                !     sum_v = b(k, j, i)*ddz(k)*ddy(j)*ddx(i) &
                !         + b(k, j, i+1)*ddz(k)*ddy(j)*ddx(i+1) &
                !         + b(k, j+1, i)*ddz(k)*ddy(j+1)*ddx(i) &
                !         + b(k, j+1, i+1)*ddz(k)*ddy(j+1)*ddx(i+1) &
                !         + b(k+1, j, i)*ddz(k+1)*ddy(j)*ddx(i) &
                !         + b(k+1, j, i+1)*ddz(k+1)*ddy(j)*ddx(i+1) &
                !         + b(k+1, j+1, i)*ddz(k+1)*ddy(j+1)*ddx(i) &
                !         + b(k+1, j+1, i+1)*ddz(k+1)*ddy(j+1)*ddx(i+1)

                !     icount = icount + 1
                !     sendbuf(icount) = divide0(sum_pv, sum_v)
                ! END DO
            END DO
        END DO
    END SUBROUTINE restrict_gc_p_t


    SUBROUTINE restrict_gc_r(kk, jj, ii, ff, ddx, ddy, ddz, bp, &
            messagelength, sbuf, istart, istop, jstart, jstop, kstart, kstop)
        !$omp declare target
        ! Subroutine arguments
        INTEGER(intk), INTENT(in) :: kk, jj, ii
        REAL(realk), INTENT(in) :: ff(kk, jj, ii)
        REAL(realk), INTENT(in) :: ddx(ii), ddy(jj), ddz(kk)
        REAL(realk), INTENT(in) :: bp(kk, jj, ii)
        INTEGER(int32), INTENT(in) :: messagelength
        REAL(realk), INTENT(inout) :: sbuf(messagelength)
        INTEGER(intk), INTENT(in) :: istart, istop, jstart, jstop, kstart, kstop

        ! Local variables
        INTEGER(intk) :: i, j, k, icount
        INTEGER(intk) :: nk, ink

        REAL(realk) :: pv1(kk), pv2(kk), pv3(kk), pv4(kk)
        REAL(realk) :: sum_pv, sum_v

        ! Number of k iterations
        ! (not to be confused with kk - these are not the same!!!)
        nk = (kstop-kstart)/2+1

        icount = 0
        DO i = istart, istop, 2
            DO j = jstart, jstop, 2

                ! Vomume
                DO k = kstart, kstop+1
                    pv1(k) = ff(k, j, i)*bp(k, j, i)*ddz(k)*ddy(j)*ddx(i)
                END DO
                DO k = kstart, kstop+1
                    pv2(k) = ff(k, j, i+1)*bp(k, j, i+1)*ddz(k)*ddy(j)*ddx(i+1)
                END DO
                DO k = kstart, kstop+1
                    pv3(k) = ff(k, j+1, i)*bp(k, j+1, i)*ddz(k)*ddy(j+1)*ddx(i)
                END DO
                DO k = kstart, kstop+1
                    pv4(k) = ff(k, j+1, i+1)*bp(k, j+1, i+1) &
                        *ddz(k)*ddy(j+1)*ddx(i+1)
                END DO

                ! Sum up and divide
                DO ink = 1, nk
                    k = kstart + 2*(ink-1)

                    sum_pv = pv1(k) + pv2(k) + pv3(k) + pv4(k) &
                        + pv1(k+1) + pv2(k+1) + pv3(k+1) + pv4(k+1)

                    sum_v = (ddx(i) + ddx(i+1))*(ddy(j) + ddy(j+1)) &
                        *(ddz(k) + ddz(k+1))

                    sbuf(icount + ink) = sum_pv/sum_v
                END DO

                ! Increment counter
                icount = icount + nk

                ! Legacy code - keep for reference
                ! DO k = kstart, kstop, 2
                !     sum_pv = ff(k, j, i)*bp(k, j, i)*ddz(k)*ddy(j)*ddx(i) &
                !         + ff(k, j, i+1)*bp(k, j, i+1)*ddz(k)*ddy(j)*ddx(i+1) &
                !         + ff(k, j+1, i)*bp(k, j+1, i)*ddz(k)*ddy(j+1)*ddx(i) &
                !         + ff(k, j+1, i+1)*bp(k, j+1, i+1)*ddz(k)*ddy(j+1)*ddx(i+1) &
                !         + ff(k+1, j, i)*bp(k+1, j, i)*ddz(k+1)*ddy(j)*ddx(i) &
                !         + ff(k+1, j, i+1)*bp(k+1, j, i+1)*ddz(k+1)*ddy(j)*ddx(i+1) &
                !         + ff(k+1, j+1, i)*bp(k+1, j+1, i)*ddz(k+1)*ddy(j+1)*ddx(i) &
                !         + ff(k+1, j+1, i+1)*bp(k+1, j+1, i+1)*ddz(k+1)*ddy(j+1)*ddx(i+1)

                !     sum_v = (ddx(i) + ddx(i+1))*(ddy(j) + ddy(j+1)) &
                !         *(ddz(k) + ddz(k+1))

                !     icount = icount + 1
                !     sendbuf(icount) = sum_pv/sum_v
                ! END DO
            END DO
        END DO
    END SUBROUTINE restrict_gc_r


    SUBROUTINE restrict_gc_n(kk, jj, ii, ff, messagelength, sbuf, &
            istart, istop, jstart, jstop, kstart, kstop)
        !$omp declare target
        ! Subroutine arguments
        INTEGER(intk), INTENT(in) :: kk, jj, ii
        REAL(realk), INTENT(in) :: ff(kk, jj, ii)
        INTEGER(int32), INTENT(in) :: messagelength
        REAL(realk), INTENT(inout) :: sbuf(messagelength)
        INTEGER(intk), INTENT(in) :: istart, istop, jstart, jstop, kstart, kstop

        ! Local variables
        INTEGER(intk) :: i, j, k, icount

        icount = 0
        DO i = istart, istop, 2
            DO j = jstart, jstop, 2
                DO k = kstart, kstop, 2
                    icount = icount + 1
                    sbuf(icount) = ff(k, j, i)
                END DO
            END DO
        END DO
    END SUBROUTINE restrict_gc_n


    SUBROUTINE restrict_gc_e(kk, jj, ii, ff, messagelength, sbuf, &
            istart, istop, jstart, jstop, kstart, kstop)
        !$omp declare target
        ! Subroutine arguments
        INTEGER(intk), INTENT(in) :: kk, jj, ii
        REAL(realk), INTENT(in) :: ff(kk, jj, ii)
        INTEGER(int32), INTENT(in) :: messagelength
        REAL(realk), INTENT(inout) :: sbuf(messagelength)
        INTEGER(intk), INTENT(in) :: istart, istop, jstart, jstop, kstart, kstop

        ! Local variables
        INTEGER(intk) :: i, j, k, icount

        icount = 0
        DO i = istart, istop, 2
            DO j = jstart, jstop, 2
                DO k = kstart, kstop, 2
                    icount = icount + 1
                    ! TODO: Check purpose of MIN(1.0, ...) MAX(...) is never
                    ! greater than 1 ???
                    sbuf(icount) = MIN(1.0_realk, MAX(ff(k, j, i), &
                        ff(k, j, i+1), &
                        ff(k, j+1, i), &
                        ff(k, j+1, i+1), &
                        ff(k+1, j, i), &
                        ff(k+1, j, i+1), &
                        ff(k+1, j+1, i), &
                        ff(k+1, j+1, i+1)))
                END DO
            END DO
        END DO
    END SUBROUTINE restrict_gc_e


    SUBROUTINE restrict_gc_f(kk, jj, ii, ff, messagelength, sbuf, &
            istart, istop, jstart, jstop, kstart, kstop)
        !$omp declare target
        ! Subroutine arguments
        INTEGER(intk), INTENT(in) :: kk, jj, ii
        REAL(realk), INTENT(in) :: ff(kk, jj, ii)
        INTEGER(int32), INTENT(in) :: messagelength
        REAL(realk), INTENT(inout) :: sbuf(messagelength)
        INTEGER(intk), INTENT(in) :: istart, istop, jstart, jstop, kstart, kstop

        ! Local variables
        INTEGER(intk) :: i, j, k, icount

        icount = 0
        DO i = istart, istop, 2
            DO j = jstart, jstop, 2
                DO k = kstart, kstop, 2
                    icount = icount + 1
                    sbuf(icount) = MIN(2.0_realk, MAX(ff(k, j, i), &
                        ff(k, j, i+1), &
                        ff(k, j+1, i), &
                        ff(k, j+1, i+1), &
                        ff(k+1, j, i), &
                        ff(k+1, j, i+1), &
                        ff(k+1, j+1, i), &
                        ff(k+1, j+1, i+1)))
                END DO
            END DO
        END DO
    END SUBROUTINE restrict_gc_f


    SUBROUTINE restrict_gc_i(kk, jj, ii, ff, messagelength, sbuf, &
            istart, istop, jstart, jstop, kstart, kstop)
        !$omp declare target
        ! Subroutine arguments
        INTEGER(intk), INTENT(in) :: kk, jj, ii
        REAL(realk), INTENT(in) :: ff(kk, jj, ii)
        INTEGER(int32), INTENT(in) :: messagelength
        REAL(realk), INTENT(inout) :: sbuf(messagelength)
        INTEGER(intk), INTENT(in) :: istart, istop, jstart, jstop, kstart, kstop

        ! Local variables
        INTEGER(intk) :: i, j, k, icount

        icount = 0
        DO i = istart, istop, 2
            DO j = jstart, jstop, 2
                DO k = kstart, kstop, 2
                    icount = icount + 1
                    sbuf(icount) = MAX(ff(k, j, i), &
                        ff(k, j, i+1), &
                        ff(k, j+1, i), &
                        ff(k, j+1, i+1), &
                        ff(k+1, j, i), &
                        ff(k+1, j, i+1), &
                        ff(k+1, j+1, i), &
                        ff(k+1, j+1, i+1))

                    ! TODO: implement this as better alternative
                    ! see util_mod.F90 most_frequent_nonzero
                    ! list(1) =  NINT(ff(k, j, i))
                    ! list(2) =  NINT(ff(k, j, i+1))
                    ! list(3) =  NINT(ff(k, j+1, i))
                    ! list(4) =  NINT(ff(k, j+1, i+1))
                    ! list(5) =  NINT(ff(k+1, j, i))
                    ! list(6) =  NINT(ff(k+1, j, i+1))
                    ! list(7) =  NINT(ff(k+1, j+1, i))
                    ! list(8) =  NINT(ff(k+1, j+1, i+1))
                    ! sendbuf(icount) = REAL(most_frequent_nonzero(list), realk)
                END DO
            END DO
        END DO
    END SUBROUTINE restrict_gc_i
END MODULE gc_restrict_mod
