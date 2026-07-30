MODULE noib_restrict_mod
    USE core_mod

    IMPLICIT NONE(type, external)
    PRIVATE

    PUBLIC :: restrict_noib_u, restrict_noib_v, restrict_noib_w, restrict_noib_s
CONTAINS
    SUBROUTINE restrict_noib_u(kk, jj, ii, ff, ddx, ddy, ddz, messagelength, &
            sbuf, istart, istop, jstart, jstop, kstart, kstop)
        !$omp declare target
        ! Subroutine arguments
        INTEGER(intk), INTENT(in) :: kk, jj, ii
        REAL(realk), INTENT(in) :: ff(kk, jj, ii)
        REAL(realk), INTENT(in) :: ddx(ii), ddy(jj), ddz(kk)
        INTEGER(int32), INTENT(in) :: messagelength
        REAL(realk), DIMENSION(messagelength), INTENT(inout) :: sbuf
        INTEGER(intk), INTENT(in) :: istart, istop, jstart, jstop, kstart, kstop

        ! Local variables
        INTEGER(intk) :: i, j, k, icount
        REAL(realk) :: sum_ua, sum_a

        icount = 0
        DO i = istart, istop, 2
            DO j = jstart, jstop, 2
                DO k = kstart, kstop, 2
                    sum_ua = ff(k, j, i)*ddy(j)*ddz(k) &
                        + ff(k, j+1, i)*ddy(j+1)*ddz(k) &
                        + ff(k+1, j, i)*ddy(j)*ddz(k+1) &
                        + ff(k+1, j+1, i)*ddy(j+1)*ddz(k+1)

                    sum_a = (ddy(j) + ddy(j+1))*(ddz(k) + ddz(k+1))

                    icount = icount + 1
                    sbuf(icount) = sum_ua/sum_a
                END DO
            END DO
        END DO
    END SUBROUTINE restrict_noib_u


    SUBROUTINE restrict_noib_v(kk, jj, ii, ff, ddx, ddy, ddz, messagelength, &
            sbuf, istart, istop, jstart, jstop, kstart, kstop)
        !$omp declare target
        ! Subroutine arguments
        INTEGER(intk), INTENT(in) :: kk, jj, ii
        REAL(realk), INTENT(in) :: ff(kk, jj, ii)
        REAL(realk), INTENT(in) :: ddx(ii), ddy(jj), ddz(kk)
        INTEGER(int32), INTENT(in) :: messagelength
        REAL(realk), INTENT(inout) :: sbuf(messagelength)
        INTEGER(intk), INTENT(in) :: istart, istop, jstart, jstop, kstart, kstop

        ! Local variables
        INTEGER(intk) :: i, j, k, icount
        REAL(realk) :: sum_va, sum_a

        icount = 0
        DO i = istart, istop, 2
            DO j = jstart, jstop, 2
                DO k = kstart, kstop, 2
                    sum_va = ff(k, j, i)*ddx(i)*ddz(k) &
                        + ff(k, j, i+1)*ddx(i+1)*ddz(k) &
                        + ff(k+1, j, i)*ddx(i)*ddz(k+1) &
                        + ff(k+1, j, i+1)*ddx(i+1)*ddz(k+1)

                    sum_a = (ddx(i) + ddx(i+1))*(ddz(k) + ddz(k+1))

                    icount = icount + 1
                    sbuf(icount) = sum_va/sum_a
                END DO
            END DO
        END DO
    END SUBROUTINE restrict_noib_v


    SUBROUTINE restrict_noib_w(kk, jj, ii, ff, ddx, ddy, ddz, messagelength, &
            sbuf, istart, istop, jstart, jstop, kstart, kstop)
        !$omp declare target
        ! Subroutine arguments
        INTEGER(intk), INTENT(in) :: kk, jj, ii
        REAL(realk), INTENT(in) :: ff(kk, jj, ii)
        REAL(realk), INTENT(in) :: ddx(ii), ddy(jj), ddz(kk)
        INTEGER(int32), INTENT(in) :: messagelength
        REAL(realk), INTENT(inout) :: sbuf(messagelength)
        INTEGER(intk), INTENT(in) :: istart, istop, jstart, jstop, kstart, kstop

        ! Local variables
        INTEGER(intk) :: i, j, k, icount
        REAL(realk) :: sum_wa, sum_a

        icount = 0
        DO i = istart, istop, 2
            DO j = jstart, jstop, 2
                DO k = kstart, kstop, 2
                    sum_wa = ff(k, j, i)*ddx(i)*ddy(j) &
                        + ff(k, j, i+1)*ddx(i+1)*ddy(j) &
                        + ff(k, j+1, i)*ddx(i)*ddy(j+1) &
                        + ff(k, j+1, i+1)*ddx(i+1)*ddy(j+1)

                    sum_a = (ddx(i) + ddx(i+1))*(ddy(j) + ddy(j+1))

                    icount = icount + 1
                    sbuf(icount) = sum_wa/sum_a
                END DO
            END DO
        END DO
    END SUBROUTINE restrict_noib_w


    SUBROUTINE restrict_noib_s(kk, jj, ii, ff, ddx, ddy, ddz, messagelength, &
            sbuf, istart, istop, jstart, jstop, kstart, kstop)
        !$omp declare target
        ! Subroutine arguments
        INTEGER(intk), INTENT(in) :: kk, jj, ii
        REAL(realk), INTENT(in) :: ff(kk, jj, ii)
        REAL(realk), INTENT(in) :: ddx(ii), ddy(jj), ddz(kk)
        INTEGER(int32), INTENT(in) :: messagelength
        REAL(realk), INTENT(inout) :: sbuf(messagelength)
        INTEGER(intk), INTENT(in) :: istart, istop, jstart, jstop, kstart, kstop

        ! Local variables
        INTEGER(intk) :: i, j, k, icount
        REAL(realk) :: sum_pv, sum_v

        icount = 0
        DO i = istart, istop, 2
            DO j = jstart, jstop, 2
                DO k = kstart, kstop, 2
                    sum_pv = ff(k, j, i)*ddz(k)*ddy(j)*ddx(i) &
                        + ff(k, j, i+1)*ddz(k)*ddy(j)*ddx(i+1) &
                        + ff(k, j+1, i)*ddz(k)*ddy(j+1)*ddx(i) &
                        + ff(k, j+1, i+1)*ddz(k)*ddy(j+1)*ddx(i+1) &
                        + ff(k+1, j, i)*ddz(k+1)*ddy(j)*ddx(i) &
                        + ff(k+1, j, i+1)*ddz(k+1)*ddy(j)*ddx(i+1) &
                        + ff(k+1, j+1, i)*ddz(k+1)*ddy(j+1)*ddx(i) &
                        + ff(k+1, j+1, i+1)*ddz(k+1)*ddy(j+1)*ddx(i+1)

                    sum_v = (ddx(i) + ddx(i+1))*(ddy(j) + ddy(j+1)) &
                        *(ddz(k) + ddz(k+1))

                    icount = icount + 1
                    sbuf(icount) = sum_pv/sum_v
                END DO
            END DO
        END DO
    END SUBROUTINE restrict_noib_s
END MODULE noib_restrict_mod
