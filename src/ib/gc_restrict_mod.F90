MODULE gc_restrict_mod
    USE core_mod, ONLY: realk, intk, get_fieldptr, get_mgdims, errr
    USE noib_restrict_mod, ONLY: noib_restrict_t

    IMPLICIT NONE(type, external)
    PRIVATE

    TYPE, EXTENDS(noib_restrict_t) :: gc_restrict_t
    CONTAINS
        PROCEDURE :: restrict
        PROCEDURE :: restrict_p
        PROCEDURE :: restrict_r
        PROCEDURE :: restrict_n
        PROCEDURE :: restrict_e
        PROCEDURE :: restrict_f
        PROCEDURE :: restrict_i
    END TYPE gc_restrict_t

    PUBLIC :: gc_restrict_t

CONTAINS
    SUBROUTINE restrict(this, kk, jj, ii, ff, sendbuf, ctyp, igrid)
        ! Subroutine arguments
        CLASS(gc_restrict_t), INTENT(inout) :: this
        INTEGER(intk), INTENT(IN) :: kk, jj, ii
        REAL(realk), INTENT(IN) :: ff(kk, jj, ii)
        REAL(realk), CONTIGUOUS, INTENT(INOUT) :: sendbuf(:)
        CHARACTER(len=1), INTENT(in) :: ctyp
        INTEGER(intk), INTENT(in) :: igrid

        SELECT CASE (ctyp)
        CASE ("U")
            CALL this%restrict_u(kk, jj, ii, ff, sendbuf, ctyp, igrid)
        CASE ("V")
            CALL this%restrict_v(kk, jj, ii, ff, sendbuf, ctyp, igrid)
        CASE ("W")
            CALL this%restrict_w(kk, jj, ii, ff, sendbuf, ctyp, igrid)
        CASE ("P", "T")
            CALL this%restrict_p(kk, jj, ii, ff, sendbuf, ctyp, igrid)
        CASE ("R")
            CALL this%restrict_r(kk, jj, ii, ff, sendbuf, ctyp, igrid)
        CASE ("S")
            CALL this%restrict_s(kk, jj, ii, ff, sendbuf, ctyp, igrid)
        CASE ("N")
            CALL this%restrict_n(kk, jj, ii, ff, sendbuf, ctyp, igrid)
        CASE ("E")
            CALL this%restrict_e(kk, jj, ii, ff, sendbuf, ctyp, igrid)
        CASE ("F")
            CALL this%restrict_f(kk, jj, ii, ff, sendbuf, ctyp, igrid)
        CASE ("I")
            CALL this%restrict_i(kk, jj, ii, ff, sendbuf, ctyp, igrid)
        CASE DEFAULT
            CALL errr(__FILE__, __LINE__)
        END SELECT
    END SUBROUTINE restrict


    SUBROUTINE restrict_p(this, kk, jj, ii, ff, sendbuf, ctyp, igrid)
        USE core_mod, ONLY: divide0

        ! Subroutine arguments
        CLASS(gc_restrict_t), INTENT(inout) :: this
        INTEGER(intk), INTENT(IN) :: kk, jj, ii
        REAL(realk), INTENT(IN) :: ff(kk, jj, ii)
        REAL(realk), CONTIGUOUS, INTENT(INOUT) :: sendbuf(:)
        CHARACTER(len=1), INTENT(in) :: ctyp
        INTEGER(intk), INTENT(in) :: igrid

        ! Local variables
        INTEGER(intk) :: i, j, k, icount
        INTEGER(intk) :: istart, istop, jstart, jstop, kstart, kstop
        REAL(realk) :: sum_pv, sum_v
        REAL(realk), POINTER, CONTIGUOUS :: ddx(:), ddy(:), ddz(:), bp(:, :, :)

        CALL this%start_and_stop(istart, istop, jstart, jstop, &
            kstart, kstop, ctyp, igrid)

        CALL get_fieldptr(ddx, "DDX", igrid)
        CALL get_fieldptr(ddy, "DDY", igrid)
        CALL get_fieldptr(ddz, "DDZ", igrid)

        IF (ctyp == 'T') THEN
            CALL get_fieldptr(bp, "BT", igrid)
        ELSE
            CALL get_fieldptr(bp, "BP", igrid)
        END IF

        icount = 0
        DO i = istart, istop, 2
            DO j = jstart, jstop, 2
                DO k = kstart, kstop, 2
                    sum_pv = ff(k, j, i)*bp(k, j, i)*ddz(k)*ddy(j)*ddx(i) &
                        + ff(k, j, i+1)*bp(k, j, i+1)*ddz(k)*ddy(j)*ddx(i+1) &
                        + ff(k, j+1, i)*bp(k, j+1, i)*ddz(k)*ddy(j+1)*ddx(i) &
                        + ff(k, j+1, i+1)*bp(k, j+1, i+1)*ddz(k)*ddy(j+1)*ddx(i+1) &
                        + ff(k+1, j, i)*bp(k+1, j, i)*ddz(k+1)*ddy(j)*ddx(i) &
                        + ff(k+1, j, i+1)*bp(k+1, j, i+1)*ddz(k+1)*ddy(j)*ddx(i+1) &
                        + ff(k+1, j+1, i)*bp(k+1, j+1, i)*ddz(k+1)*ddy(j+1)*ddx(i) &
                        + ff(k+1, j+1, i+1)*bp(k+1, j+1, i+1)*ddz(k+1)*ddy(j+1)*ddx(i+1)

                    sum_v = bp(k, j, i)*ddz(k)*ddy(j)*ddx(i) &
                        + bp(k, j, i+1)*ddz(k)*ddy(j)*ddx(i+1) &
                        + bp(k, j+1, i)*ddz(k)*ddy(j+1)*ddx(i) &
                        + bp(k, j+1, i+1)*ddz(k)*ddy(j+1)*ddx(i+1) &
                        + bp(k+1, j, i)*ddz(k+1)*ddy(j)*ddx(i) &
                        + bp(k+1, j, i+1)*ddz(k+1)*ddy(j)*ddx(i+1) &
                        + bp(k+1, j+1, i)*ddz(k+1)*ddy(j+1)*ddx(i) &
                        + bp(k+1, j+1, i+1)*ddz(k+1)*ddy(j+1)*ddx(i+1)

                    icount = icount + 1
                    sendbuf(icount) = divide0(sum_pv, sum_v)
                END DO
            END DO
        END DO
    END SUBROUTINE restrict_p


    SUBROUTINE restrict_r(this, kk, jj, ii, ff, sendbuf, ctyp, igrid)
        ! Subroutine arguments
        CLASS(gc_restrict_t), INTENT(inout) :: this
        INTEGER(intk), INTENT(IN) :: kk, jj, ii
        REAL(realk), INTENT(IN) :: ff(kk, jj, ii)
        REAL(realk), CONTIGUOUS, INTENT(INOUT) :: sendbuf(:)
        CHARACTER(len=1), INTENT(in) :: ctyp
        INTEGER(intk), INTENT(in) :: igrid

        ! Local variables
        INTEGER(intk) :: i, j, k, icount
        INTEGER(intk) :: istart, istop, jstart, jstop, kstart, kstop
        REAL(realk) :: sum_pv, sum_v
        REAL(realk), POINTER, CONTIGUOUS :: ddx(:), ddy(:), ddz(:), bp(:, :, :)

        CALL this%start_and_stop(istart, istop, jstart, jstop, &
            kstart, kstop, ctyp, igrid)

        CALL get_fieldptr(ddx, "DDX", igrid)
        CALL get_fieldptr(ddy, "DDY", igrid)
        CALL get_fieldptr(ddz, "DDZ", igrid)
        CALL get_fieldptr(bp, "BP", igrid)

        icount = 0
        DO i = istart, istop, 2
            DO j = jstart, jstop, 2
                DO k = kstart, kstop, 2
                    sum_pv = ff(k, j, i)*bp(k, j, i)*ddz(k)*ddy(j)*ddx(i) &
                        + ff(k, j, i+1)*bp(k, j, i+1)*ddz(k)*ddy(j)*ddx(i+1) &
                        + ff(k, j+1, i)*bp(k, j+1, i)*ddz(k)*ddy(j+1)*ddx(i) &
                        + ff(k, j+1, i+1)*bp(k, j+1, i+1)*ddz(k)*ddy(j+1)*ddx(i+1) &
                        + ff(k+1, j, i)*bp(k+1, j, i)*ddz(k+1)*ddy(j)*ddx(i) &
                        + ff(k+1, j, i+1)*bp(k+1, j, i+1)*ddz(k+1)*ddy(j)*ddx(i+1) &
                        + ff(k+1, j+1, i)*bp(k+1, j+1, i)*ddz(k+1)*ddy(j+1)*ddx(i) &
                        + ff(k+1, j+1, i+1)*bp(k+1, j+1, i+1)*ddz(k+1)*ddy(j+1)*ddx(i+1)

                    sum_v = (ddx(i) + ddx(i+1))*(ddy(j) + ddy(j+1)) &
                        *(ddz(k) + ddz(k+1))

                    icount = icount + 1
                    sendbuf(icount) = sum_pv/sum_v
                END DO
            END DO
        END DO
    END SUBROUTINE restrict_r


    SUBROUTINE restrict_n(this, kk, jj, ii, ff, sendbuf, ctyp, igrid)
        ! Subroutine arguments
        CLASS(gc_restrict_t), INTENT(inout) :: this
        INTEGER(intk), INTENT(IN) :: kk, jj, ii
        REAL(realk), INTENT(IN) :: ff(kk, jj, ii)
        REAL(realk), CONTIGUOUS, INTENT(INOUT) :: sendbuf(:)
        CHARACTER(len=1), INTENT(in) :: ctyp
        INTEGER(intk), INTENT(in) :: igrid

        ! Local variables
        INTEGER(intk) :: i, j, k, icount
        INTEGER(intk) :: istart, istop, jstart, jstop, kstart, kstop

        CALL this%start_and_stop(istart, istop, jstart, jstop, &
            kstart, kstop, ctyp, igrid)

        icount = 0
        DO i = istart, istop, 2
            DO j = jstart, jstop, 2
                DO k = kstart, kstop, 2
                    icount = icount + 1
                    sendbuf(icount) = ff(k, j, i)
                END DO
            END DO
        END DO
    END SUBROUTINE restrict_n


    SUBROUTINE restrict_e(this, kk, jj, ii, ff, sendbuf, ctyp, igrid)
        ! Subroutine arguments
        CLASS(gc_restrict_t), INTENT(inout) :: this
        INTEGER(intk), INTENT(IN) :: kk, jj, ii
        REAL(realk), INTENT(IN) :: ff(kk, jj, ii)
        REAL(realk), CONTIGUOUS, INTENT(INOUT) :: sendbuf(:)
        CHARACTER(len=1), INTENT(in) :: ctyp
        INTEGER(intk), INTENT(in) :: igrid

        ! Local variables
        INTEGER(intk) :: i, j, k, icount
        INTEGER(intk) :: istart, istop, jstart, jstop, kstart, kstop

        CALL this%start_and_stop(istart, istop, jstart, jstop, &
            kstart, kstop, ctyp, igrid)

        icount = 0
        DO i = istart, istop, 2
            DO j = jstart, jstop, 2
                DO k = kstart, kstop, 2
                    icount = icount + 1
                    ! TODO: Check purpose of MIN(1.0, ...) MAX(...) is never
                    ! greater than 1 ???
                    sendbuf(icount) = MIN(1.0_realk, MAX(ff(k, j, i), &
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
    END SUBROUTINE restrict_e


    SUBROUTINE restrict_f(this, kk, jj, ii, ff, sendbuf, ctyp, igrid)
        ! Subroutine arguments
        CLASS(gc_restrict_t), INTENT(inout) :: this
        INTEGER(intk), INTENT(IN) :: kk, jj, ii
        REAL(realk), INTENT(IN) :: ff(kk, jj, ii)
        REAL(realk), CONTIGUOUS, INTENT(INOUT) :: sendbuf(:)
        CHARACTER(len=1), INTENT(in) :: ctyp
        INTEGER(intk), INTENT(in) :: igrid

        ! Local variables
        INTEGER(intk) :: i, j, k, icount
        INTEGER(intk) :: istart, istop, jstart, jstop, kstart, kstop

        CALL this%start_and_stop(istart, istop, jstart, jstop, &
            kstart, kstop, ctyp, igrid)

        icount = 0
        DO i = istart, istop, 2
            DO j = jstart, jstop, 2
                DO k = kstart, kstop, 2
                    icount = icount + 1
                    sendbuf(icount) = MIN(2.0_realk, MAX(ff(k, j, i), &
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
    END SUBROUTINE restrict_f


    SUBROUTINE restrict_i(this, kk, jj, ii, ff, sendbuf, ctyp, igrid)
        ! Subroutine arguments
        CLASS(gc_restrict_t), INTENT(inout) :: this
        INTEGER(intk), INTENT(IN) :: kk, jj, ii
        REAL(realk), INTENT(IN) :: ff(kk, jj, ii)
        REAL(realk), CONTIGUOUS, INTENT(INOUT) :: sendbuf(:)
        CHARACTER(len=1), INTENT(in) :: ctyp
        INTEGER(intk), INTENT(in) :: igrid

        ! Local variables
        INTEGER(intk) :: i, j, k, icount
        INTEGER(intk) :: istart, istop, jstart, jstop, kstart, kstop

        CALL this%start_and_stop(istart, istop, jstart, jstop, &
            kstart, kstop, ctyp, igrid)

        icount = 0
        DO i = istart, istop, 2
            DO j = jstart, jstop, 2
                DO k = kstart, kstop, 2
                    icount = icount + 1
                    sendbuf(icount) = MAX(ff(k, j, i), &
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
    END SUBROUTINE restrict_i
END MODULE gc_restrict_mod
