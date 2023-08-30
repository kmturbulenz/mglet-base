MODULE noib_restrict_mod
    USE core_mod, ONLY: realk, intk, get_mgdims, errr, get_fieldptr
    USE ibmodel_mod, ONLY: restrict_t

    IMPLICIT NONE(type, external)
    PRIVATE

    TYPE, EXTENDS(restrict_t) :: noib_restrict_t
    CONTAINS
        PROCEDURE, NOPASS :: start_and_stop
        PROCEDURE :: restrict
        PROCEDURE :: restrict_u
        PROCEDURE :: restrict_v
        PROCEDURE :: restrict_w
        PROCEDURE :: restrict_s
    END TYPE noib_restrict_t

    PUBLIC :: noib_restrict_t

CONTAINS
    SUBROUTINE restrict(this, kk, jj, ii, ff, sendbuf, ctyp, igrid)
        ! Subroutine arguments
        CLASS(noib_restrict_t), INTENT(inout) :: this
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
        CASE ("P", "R", "S", "T")
            CALL this%restrict_s(kk, jj, ii, ff, sendbuf, ctyp, igrid)
        CASE DEFAULT
            CALL errr(__FILE__, __LINE__)
        END SELECT
    END SUBROUTINE restrict


    SUBROUTINE start_and_stop(ista, isto, jsta, jsto, &
            ksta, ksto, ctyp, igrid)
        USE core_mod, ONLY: get_mgdims

        ! Subroutine arguments
        INTEGER(intk), INTENT(out) :: ista, isto, jsta, jsto, ksta, ksto
        CHARACTER(len=1), INTENT(in) :: ctyp
        INTEGER(intk), INTENT(in) :: igrid

        ! Local variables
        INTEGER(intk) :: ii, jj, kk
        INTEGER(intk) :: istart, iend, jstart, jend, kstart, kend

        ! Loop ranges are DO i = istart, ii-iend
        ! returned as DO i = ista, isto, 2
        SELECT CASE (ctyp)
        CASE ("U")
            istart = 2
            iend = 2
            jstart = 3
            jend = 3
            kstart = 3
            kend = 3
        CASE ("V")
            istart = 3
            iend = 3
            jstart = 2
            jend = 2
            kstart = 3
            kend = 3
        CASE ("W")
            istart = 3
            iend = 3
            jstart = 3
            jend = 3
            kstart = 2
            kend = 2
        CASE ("E", "F", "P", "R", "S", 'I', 'T')
            istart = 3
            iend = 3
            jstart = 3
            jend = 3
            kstart = 3
            kend = 3
        CASE ("N")
            istart = 2
            iend = 2
            jstart = 2
            jend = 2
            kstart = 2
            kend = 2
        CASE DEFAULT
            CALL errr(__FILE__, __LINE__)
        END SELECT

        CALL get_mgdims(kk, jj, ii, igrid)

        ista = istart
        isto = ii - iend
        jsta = jstart
        jsto = jj - jend
        ksta = kstart
        ksto = kk - kend

        ! Since the loop has strides of 2, do a sanity check that
        ! the upper index is correct. The stop-index MUST be the actual
        ! index of the last iteration for the message size to be computed
        ! correctly
        !
        ! E.g. The loop:
        !   DO i = 3, kk-3, 2
        ! with kk = 28 will iterate like:
        !   i = 3, 5, 7, 9, 11, 13, 15, 17, 19, 21, 23, 25
        ! with 12 iterations.
        !
        ! If the loop had been:
        !   DO i = 3, kk-2, 2
        ! the iterations would have been the same (the last one being i = 25)
        ! but the message size would be wrong.
        !
        ! Therefore require that MOD(isto-ista, 2) == 0
        IF (MOD(isto-ista, 2) > 0) CALL errr(__FILE__, __LINE__)
        IF (MOD(jsto-jsta, 2) > 0) CALL errr(__FILE__, __LINE__)
        IF (MOD(ksto-ksta, 2) > 0) CALL errr(__FILE__, __LINE__)
    END SUBROUTINE start_and_stop


    SUBROUTINE restrict_u(this, kk, jj, ii, ff, sendbuf, ctyp, igrid)
        ! Subroutine arguments
        CLASS(noib_restrict_t), INTENT(inout) :: this
        INTEGER(intk), INTENT(IN) :: kk, jj, ii
        REAL(realk), INTENT(IN) :: ff(kk, jj, ii)
        REAL(realk), CONTIGUOUS, INTENT(INOUT) :: sendbuf(:)
        CHARACTER(len=1), INTENT(in) :: ctyp
        INTEGER(intk), INTENT(in) :: igrid

        ! Local variables
        INTEGER(intk) :: i, j, k, icount
        INTEGER(intk) :: istart, istop, jstart, jstop, kstart, kstop
        REAL(realk) :: sum_ua, sum_a
        REAL(realk), POINTER, CONTIGUOUS :: ddx(:), ddy(:), ddz(:)

        CALL this%start_and_stop(istart, istop, jstart, jstop, &
            kstart, kstop, ctyp, igrid)

        CALL get_fieldptr(ddx, "DDX", igrid)
        CALL get_fieldptr(ddy, "DDY", igrid)
        CALL get_fieldptr(ddz, "DDZ", igrid)

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
                    sendbuf(icount) = sum_ua/sum_a
                END DO
            END DO
        END DO
    END SUBROUTINE restrict_u


    SUBROUTINE restrict_v(this, kk, jj, ii, ff, sendbuf, ctyp, igrid)
        ! Subroutine arguments
        CLASS(noib_restrict_t), INTENT(inout) :: this
        INTEGER(intk), INTENT(IN) :: kk, jj, ii
        REAL(realk), INTENT(IN) :: ff(kk, jj, ii)
        REAL(realk), CONTIGUOUS, INTENT(INOUT) :: sendbuf(:)
        CHARACTER(len=1), INTENT(in) :: ctyp
        INTEGER(intk), INTENT(in) :: igrid

        ! Local variables
        INTEGER(intk) :: i, j, k, icount
        INTEGER(intk) :: istart, istop, jstart, jstop, kstart, kstop
        REAL(realk) :: sum_va, sum_a
        REAL(realk), POINTER, CONTIGUOUS :: ddx(:), ddy(:), ddz(:)

        CALL this%start_and_stop(istart, istop, jstart, jstop, &
            kstart, kstop, ctyp, igrid)

        CALL get_fieldptr(ddx, "DDX", igrid)
        CALL get_fieldptr(ddy, "DDY", igrid)
        CALL get_fieldptr(ddz, "DDZ", igrid)

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
                    sendbuf(icount) = sum_va/sum_a
                END DO
            END DO
        END DO
    END SUBROUTINE restrict_v


    SUBROUTINE restrict_w(this, kk, jj, ii, ff, sendbuf, ctyp, igrid)
        ! Subroutine arguments
        CLASS(noib_restrict_t), INTENT(inout) :: this
        INTEGER(intk), INTENT(IN) :: kk, jj, ii
        REAL(realk), INTENT(IN) :: ff(kk, jj, ii)
        REAL(realk), CONTIGUOUS, INTENT(INOUT) :: sendbuf(:)
        CHARACTER(len=1), INTENT(in) :: ctyp
        INTEGER(intk), INTENT(in) :: igrid

        ! Local variables
        INTEGER(intk) :: i, j, k, icount
        INTEGER(intk) :: istart, istop, jstart, jstop, kstart, kstop
        REAL(realk) :: sum_wa, sum_a
        REAL(realk), POINTER, CONTIGUOUS :: ddx(:), ddy(:), ddz(:)

        CALL this%start_and_stop(istart, istop, jstart, jstop, &
            kstart, kstop, ctyp, igrid)

        CALL get_fieldptr(ddx, "DDX", igrid)
        CALL get_fieldptr(ddy, "DDY", igrid)
        CALL get_fieldptr(ddz, "DDZ", igrid)

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
                    sendbuf(icount) = sum_wa/sum_a
                END DO
            END DO
        END DO
    END SUBROUTINE restrict_w


    SUBROUTINE restrict_s(this, kk, jj, ii, ff, sendbuf, ctyp, igrid)
        ! Subroutine arguments
        CLASS(noib_restrict_t), INTENT(inout) :: this
        INTEGER(intk), INTENT(IN) :: kk, jj, ii
        REAL(realk), INTENT(IN) :: ff(kk, jj, ii)
        REAL(realk), CONTIGUOUS, INTENT(INOUT) :: sendbuf(:)
        CHARACTER(len=1), INTENT(in) :: ctyp
        INTEGER(intk), INTENT(in) :: igrid

        ! Local variables
        INTEGER(intk) :: i, j, k, icount
        INTEGER(intk) :: istart, istop, jstart, jstop, kstart, kstop
        REAL(realk) :: sum_pv, sum_v
        REAL(realk), POINTER, CONTIGUOUS :: ddx(:), ddy(:), ddz(:)

        CALL this%start_and_stop(istart, istop, jstart, jstop, &
            kstart, kstop, ctyp, igrid)

        CALL get_fieldptr(ddx, "DDX", igrid)
        CALL get_fieldptr(ddy, "DDY", igrid)
        CALL get_fieldptr(ddz, "DDZ", igrid)

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
                    sendbuf(icount) = sum_pv/sum_v
                END DO
            END DO
        END DO
    END SUBROUTINE restrict_s
END MODULE noib_restrict_mod
