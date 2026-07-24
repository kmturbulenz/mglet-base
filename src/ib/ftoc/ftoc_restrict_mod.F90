

MODULE ftoc_restrict_mod

    USE core_mod

CONTAINS

    ! --- Routines used for noib_t --- (from noib_restrict_mod.F90)

    SUBROUTINE restrict_noib_u(kk, jj, ii, ff, sendbuf, ctyp, igrid)
        ! Subroutine arguments
        INTEGER(intk), INTENT(IN) :: kk, jj, ii
        REAL(realk), INTENT(IN) :: ff(kk, jj, ii)
        REAL(realk), CONTIGUOUS, INTENT(INOUT) :: sendbuf(:)
        CHARACTER(len=1), INTENT(in) :: ctyp
        INTEGER(intk), INTENT(in) :: igrid

        ! Local variables
        INTEGER(intk) :: i, j, k, icount, is, js, ks, isteps, jsteps, ksteps
        INTEGER(intk) :: istart, istop, jstart, jstop, kstart, kstop
        REAL(realk) :: sum_ua, sum_a
        REAL(realk), POINTER, CONTIGUOUS :: ddx(:), ddy(:), ddz(:)

        CALL start_and_stop_restrict(istart, istop, jstart, jstop, &
            kstart, kstop, ctyp, igrid)

        CALL get_fieldptr(ddx, "DDX", igrid)
        CALL get_fieldptr(ddy, "DDY", igrid)
        CALL get_fieldptr(ddz, "DDZ", igrid)

        ! Computing number of steps in strided loop
        isteps = (istop - istart)/2 + 1
        jsteps = (jstop - jstart)/2 + 1
        ksteps = (kstop - kstart)/2 + 1

        DO i = istart, istop, 2
            DO j = jstart, jstop, 2
                DO k = kstart, kstop, 2
                    ! Weighting fine velocities with fine areas
                    sum_ua = ff(k, j, i)*ddy(j)*ddz(k) &
                        + ff(k, j+1, i)*ddy(j+1)*ddz(k) &
                        + ff(k+1, j, i)*ddy(j)*ddz(k+1) &
                        + ff(k+1, j+1, i)*ddy(j+1)*ddz(k+1)
                    ! Computing total area
                    sum_a = (ddy(j) + ddy(j+1))*(ddz(k) + ddz(k+1))
                    ! Compute the index of the coarse cell in the send buffer
                    ks = (k - kstart)/2
                    js = (j - jstart)/2
                    is = (i - istart)/2
                    icount = 1 + ks + js*ksteps + is*jsteps*ksteps
                    ! Assign area-averages velocity to coarse cell buffer entry
                    sendbuf(icount) = sum_ua / sum_a
                END DO
            END DO
        END DO

        ! Sanity check that the number of iterations is correct
        IF (icount /= isteps*jsteps*ksteps) CALL errr(__FILE__, __LINE__)

    END SUBROUTINE restrict_noib_u


    SUBROUTINE restrict_noib_v(kk, jj, ii, ff, sendbuf, ctyp, igrid)
        ! Subroutine arguments
        INTEGER(intk), INTENT(IN) :: kk, jj, ii
        REAL(realk), INTENT(IN) :: ff(kk, jj, ii)
        REAL(realk), CONTIGUOUS, INTENT(INOUT) :: sendbuf(:)
        CHARACTER(len=1), INTENT(in) :: ctyp
        INTEGER(intk), INTENT(in) :: igrid

        ! Local variables
        INTEGER(intk) :: i, j, k, icount, is, js, ks, isteps, jsteps, ksteps
        INTEGER(intk) :: istart, istop, jstart, jstop, kstart, kstop
        REAL(realk) :: sum_va, sum_a
        REAL(realk), POINTER, CONTIGUOUS :: ddx(:), ddy(:), ddz(:)

        CALL start_and_stop_restrict(istart, istop, jstart, jstop, &
            kstart, kstop, ctyp, igrid)

        CALL get_fieldptr(ddx, "DDX", igrid)
        CALL get_fieldptr(ddy, "DDY", igrid)
        CALL get_fieldptr(ddz, "DDZ", igrid)

        ! Computing number of steps in strided loop
        isteps = (istop - istart)/2 + 1
        jsteps = (jstop - jstart)/2 + 1
        ksteps = (kstop - kstart)/2 + 1

        DO i = istart, istop, 2
            DO j = jstart, jstop, 2
                DO k = kstart, kstop, 2
                    ! Weighting fine velocities with fine areas
                    sum_va = ff(k, j, i)*ddx(i)*ddz(k) &
                        + ff(k, j, i+1)*ddx(i+1)*ddz(k) &
                        + ff(k+1, j, i)*ddx(i)*ddz(k+1) &
                        + ff(k+1, j, i+1)*ddx(i+1)*ddz(k+1)
                    ! Computing total area
                    sum_a = (ddx(i) + ddx(i+1))*(ddz(k) + ddz(k+1))
                    ! Compute the index of the coarse cell in the send buffer
                    ks = (k - kstart)/2
                    js = (j - jstart)/2
                    is = (i - istart)/2
                    icount = 1 + ks + js*ksteps + is*jsteps*ksteps
                    ! Assign area-averages velocity to coarse cell buffer entry
                    sendbuf(icount) = sum_va / sum_a
                END DO
            END DO
        END DO

        ! Sanity check that the number of iterations is correct
        IF (icount /= isteps*jsteps*ksteps) CALL errr(__FILE__, __LINE__)

    END SUBROUTINE restrict_noib_v


    SUBROUTINE restrict_noib_w(kk, jj, ii, ff, sendbuf, ctyp, igrid)
        ! Subroutine arguments
        INTEGER(intk), INTENT(IN) :: kk, jj, ii
        REAL(realk), INTENT(IN) :: ff(kk, jj, ii)
        REAL(realk), CONTIGUOUS, INTENT(INOUT) :: sendbuf(:)
        CHARACTER(len=1), INTENT(in) :: ctyp
        INTEGER(intk), INTENT(in) :: igrid

        ! Local variables
        INTEGER(intk) :: i, j, k, icount, is, js, ks, isteps, jsteps, ksteps
        INTEGER(intk) :: istart, istop, jstart, jstop, kstart, kstop
        REAL(realk) :: sum_wa, sum_a
        REAL(realk), POINTER, CONTIGUOUS :: ddx(:), ddy(:), ddz(:)

        CALL start_and_stop_restrict(istart, istop, jstart, jstop, &
            kstart, kstop, ctyp, igrid)

        CALL get_fieldptr(ddx, "DDX", igrid)
        CALL get_fieldptr(ddy, "DDY", igrid)
        CALL get_fieldptr(ddz, "DDZ", igrid)

        ! Computing number of steps in strided loop
        isteps = (istop - istart)/2 + 1
        jsteps = (jstop - jstart)/2 + 1
        ksteps = (kstop - kstart)/2 + 1

        DO i = istart, istop, 2
            DO j = jstart, jstop, 2
                DO k = kstart, kstop, 2
                    ! Weighting fine velocities with fine areas
                    sum_wa = ff(k, j, i)*ddx(i)*ddy(j) &
                        + ff(k, j, i+1)*ddx(i+1)*ddy(j) &
                        + ff(k, j+1, i)*ddx(i)*ddy(j+1) &
                        + ff(k, j+1, i+1)*ddx(i+1)*ddy(j+1)
                    ! Computing total area
                    sum_a = (ddx(i) + ddx(i+1))*(ddy(j) + ddy(j+1))
                    ! Compute the index of the coarse cell in the send buffer
                    ks = (k - kstart)/2
                    js = (j - jstart)/2
                    is = (i - istart)/2
                    icount = 1 + ks + js*ksteps + is*jsteps*ksteps
                    ! Assign area-averages velocity to coarse cell buffer entry
                    sendbuf(icount) = sum_wa/sum_a
                END DO
            END DO
        END DO

        ! Sanity check that the number of iterations is correct
        IF (icount /= isteps*jsteps*ksteps) CALL errr(__FILE__, __LINE__)

    END SUBROUTINE restrict_noib_w


    SUBROUTINE restrict_noib_s(kk, jj, ii, ff, sendbuf, ctyp, igrid)
        ! Subroutine arguments
        INTEGER(intk), INTENT(IN) :: kk, jj, ii
        REAL(realk), INTENT(IN) :: ff(kk, jj, ii)
        REAL(realk), CONTIGUOUS, INTENT(INOUT) :: sendbuf(:)
        CHARACTER(len=1), INTENT(in) :: ctyp
        INTEGER(intk), INTENT(in) :: igrid

        ! Local variables
        INTEGER(intk) :: i, j, k, icount, is, js, ks, isteps, jsteps, ksteps
        INTEGER(intk) :: istart, istop, jstart, jstop, kstart, kstop
        REAL(realk) :: sum_pv, sum_v
        REAL(realk), POINTER, CONTIGUOUS :: ddx(:), ddy(:), ddz(:)

        CALL start_and_stop_restrict(istart, istop, jstart, jstop, &
            kstart, kstop, ctyp, igrid)

        CALL get_fieldptr(ddx, "DDX", igrid)
        CALL get_fieldptr(ddy, "DDY", igrid)
        CALL get_fieldptr(ddz, "DDZ", igrid)

        ! Computing number of steps in strided loop
        isteps = (istop - istart)/2 + 1
        jsteps = (jstop - jstart)/2 + 1
        ksteps = (kstop - kstart)/2 + 1

        DO i = istart, istop, 2
            DO j = jstart, jstop, 2
                DO k = kstart, kstop, 2
                    ! Weighting fine velocities with fine volumes
                    sum_pv = ff(k, j, i)*ddz(k)*ddy(j)*ddx(i) &
                        + ff(k, j, i+1)*ddz(k)*ddy(j)*ddx(i+1) &
                        + ff(k, j+1, i)*ddz(k)*ddy(j+1)*ddx(i) &
                        + ff(k, j+1, i+1)*ddz(k)*ddy(j+1)*ddx(i+1) &
                        + ff(k+1, j, i)*ddz(k+1)*ddy(j)*ddx(i) &
                        + ff(k+1, j, i+1)*ddz(k+1)*ddy(j)*ddx(i+1) &
                        + ff(k+1, j+1, i)*ddz(k+1)*ddy(j+1)*ddx(i) &
                        + ff(k+1, j+1, i+1)*ddz(k+1)*ddy(j+1)*ddx(i+1)
                    ! Computing total volume
                    sum_v = (ddx(i) + ddx(i+1))*(ddy(j) + ddy(j+1)) &
                        *(ddz(k) + ddz(k+1))
                    ! Compute the index of the coarse cell in the send buffer
                    ks = (k - kstart)/2
                    js = (j - jstart)/2
                    is = (i - istart)/2
                    icount = 1 + ks + js*ksteps + is*jsteps*ksteps
                    ! Assign volume-averaged scalar to coarse cell buffer entry
                    sendbuf(icount) = sum_pv / sum_v
                END DO
            END DO
        END DO

        ! Sanity check that the number of iterations is correct
        IF (icount /= isteps*jsteps*ksteps) CALL errr(__FILE__, __LINE__)

    END SUBROUTINE restrict_noib_s







    ! --- Routines used for gc_t ---  (from gc_restrict_mod.F90)





    SUBROUTINE restrict_p(kk, jj, ii, ff, sendbuf, ctyp, igrid)
        ! Subroutine arguments
        INTEGER(intk), INTENT(IN) :: kk, jj, ii
        REAL(realk), INTENT(IN) :: ff(kk, jj, ii)
        REAL(realk), CONTIGUOUS, INTENT(INOUT) :: sendbuf(:)
        CHARACTER(len=1), INTENT(in) :: ctyp
        INTEGER(intk), INTENT(in) :: igrid

        ! Local variables
        INTEGER(intk) :: i, j, k, icount, is, js, ks, isteps, jsteps, ksteps
        INTEGER(intk) :: istart, istop, jstart, jstop, kstart, kstop
        REAL(realk) :: sum_pv, sum_v
        REAL(realk), POINTER, CONTIGUOUS :: ddx(:), ddy(:), ddz(:)

        CALL start_and_stop_restrict(istart, istop, jstart, jstop, &
            kstart, kstop, ctyp, igrid)

        CALL get_fieldptr(ddx, "DDX", igrid)
        CALL get_fieldptr(ddy, "DDY", igrid)
        CALL get_fieldptr(ddz, "DDZ", igrid)

        ! Computing number of steps in strided loop
        isteps = (istop - istart)/2 + 1
        jsteps = (jstop - jstart)/2 + 1
        ksteps = (kstop - kstart)/2 + 1

        DO i = istart, istop, 2
            DO j = jstart, jstop, 2
                DO k = kstart, kstop, 2
                    ! Weighting fine pressures with the open fine volumes
                    sum_pv = ff(k, j, i)*bp(k, j, i)* &
                            ddz(k)*ddy(j)*ddx(i) &
                        + ff(k, j, i+1)*bp(k, j, i+1)* &
                            ddz(k)*ddy(j)*ddx(i+1) &
                        + ff(k, j+1, i)*bp(k, j+1, i)* &
                            ddz(k)*ddy(j+1)*ddx(i) &
                        + ff(k, j+1, i+1)*bp(k, j+1, i+1)* &
                            ddz(k)*ddy(j+1)*ddx(i+1) &
                        + ff(k+1, j, i)*bp(k+1, j, i)* &
                            ddz(k+1)*ddy(j)*ddx(i) &
                        + ff(k+1, j, i+1)*bp(k+1, j, i+1)* &
                            ddz(k+1)*ddy(j)*ddx(i+1) &
                        + ff(k+1, j+1, i)*bp(k+1, j+1, i)* &
                            ddz(k+1)*ddy(j+1)*ddx(i) &
                        + ff(k+1, j+1, i+1)*bp(k+1, j+1, i+1)* &
                            ddz(k+1)*ddy(j+1)*ddx(i+1)
                    ! Computing total volume of open fine cells (can be zero)
                    sum_v = bp(k, j, i)*ddz(k)*ddy(j)*ddx(i) &
                        + bp(k, j, i+1)*ddz(k)*ddy(j)*ddx(i+1) &
                        + bp(k, j+1, i)*ddz(k)*ddy(j+1)*ddx(i) &
                        + bp(k, j+1, i+1)*ddz(k)*ddy(j+1)*ddx(i+1) &
                        + bp(k+1, j, i)*ddz(k+1)*ddy(j)*ddx(i) &
                        + bp(k+1, j, i+1)*ddz(k+1)*ddy(j)*ddx(i+1) &
                        + bp(k+1, j+1, i)*ddz(k+1)*ddy(j+1)*ddx(i) &
                        + bp(k+1, j+1, i+1)*ddz(k+1)*ddy(j+1)*ddx(i+1)
                    ! Compute the index of the coarse cell in the send buffer
                    ks = (k - kstart)/2
                    js = (j - jstart)/2
                    is = (i - istart)/2
                    icount = 1 + ks + js*ksteps + is*jsteps*ksteps
                    ! Assign volume-averaged pressure to coarse buffer entry
                    sendbuf(icount) = 0.0
                    IF (sum_v > 0.0) THEN
                        sendbuf(icount) = sum_pv / sum_v
                    END IF
                END DO
            END DO
        END DO

        ! Sanity check that the number of iterations is correct
        IF (icount /= isteps*jsteps*ksteps) CALL errr(__FILE__, __LINE__)

    END SUBROUTINE restrict_p



    SUBROUTINE restrict_r(kk, jj, ii, ff, sendbuf, ctyp, igrid)
        ! Subroutine arguments
        INTEGER(intk), INTENT(IN) :: kk, jj, ii
        REAL(realk), INTENT(IN) :: ff(kk, jj, ii)
        REAL(realk), CONTIGUOUS, INTENT(INOUT) :: sendbuf(:)
        CHARACTER(len=1), INTENT(in) :: ctyp
        INTEGER(intk), INTENT(in) :: igrid

        ! Local variables
        INTEGER(intk) :: i, j, k, icount, is, js, ks, isteps, jsteps, ksteps
        INTEGER(intk) :: istart, istop, jstart, jstop, kstart, kstop
        REAL(realk) :: sum_pv, sum_v

        REAL(realk), POINTER, CONTIGUOUS :: ddx(:), ddy(:), ddz(:), bp(:, :, :)

        CALL start_and_stop_restrict(istart, istop, jstart, jstop, &
            kstart, kstop, ctyp, igrid)

        CALL get_fieldptr(ddx, "DDX", igrid)
        CALL get_fieldptr(ddy, "DDY", igrid)
        CALL get_fieldptr(ddz, "DDZ", igrid)
        CALL get_fieldptr(bp, "BP", igrid)

        ! Computing number of steps in strided loop
        isteps = (istop - istart)/2 + 1
        jsteps = (jstop - jstart)/2 + 1
        ksteps = (kstop - kstart)/2 + 1

        DO i = istart, istop, 2
            DO j = jstart, jstop, 2
                DO k = kstart, kstop, 2
                    ! Weighting fine pressures with the open fine volumes
                    sum_pv = ff(k, j, i)*bp(k, j, i)* &
                            ddz(k)*ddy(j)*ddx(i) &
                        + ff(k, j, i+1)*bp(k, j, i+1)* &
                            ddz(k)*ddy(j)*ddx(i+1) &
                        + ff(k, j+1, i)*bp(k, j+1, i)* &
                            ddz(k)*ddy(j+1)*ddx(i) &
                        + ff(k, j+1, i+1)*bp(k, j+1, i+1)* &
                            ddz(k)*ddy(j+1)*ddx(i+1) &
                        + ff(k+1, j, i)*bp(k+1, j, i)* &
                            ddz(k+1)*ddy(j)*ddx(i) &
                        + ff(k+1, j, i+1)*bp(k+1, j, i+1)* &
                            ddz(k+1)*ddy(j)*ddx(i+1) &
                        + ff(k+1, j+1, i)*bp(k+1, j+1, i)* &
                            ddz(k+1)*ddy(j+1)*ddx(i) &
                        + ff(k+1, j+1, i+1)*bp(k+1, j+1, i+1)* &
                            ddz(k+1)*ddy(j+1)*ddx(i+1)
                    ! Computing total volume of open + closed fine cells (> 0)
                    sum_v = (ddx(i) + ddx(i+1))*(ddy(j) + ddy(j+1)) &
                        *(ddz(k) + ddz(k+1))
                    ! Compute the index of the coarse cell in the send buffer
                    ks = (k - kstart)/2
                    js = (j - jstart)/2
                    is = (i - istart)/2
                    icount = 1 + ks + js*ksteps + is*jsteps*ksteps
                    ! Assign volume-averaged pressure to coarse buffer entry
                    sendbuf(icount) = sum_pv / sum_v
                END DO
            END DO
        END DO

        ! Sanity check that the number of iterations is correct
        IF (icount /= isteps*jsteps*ksteps) CALL errr(__FILE__, __LINE__)

    END SUBROUTINE restrict_r


    SUBROUTINE restrict_n(kk, jj, ii, ff, sendbuf, ctyp, igrid)
        ! Subroutine arguments
        INTEGER(intk), INTENT(IN) :: kk, jj, ii
        REAL(realk), INTENT(IN) :: ff(kk, jj, ii)
        REAL(realk), CONTIGUOUS, INTENT(INOUT) :: sendbuf(:)
        CHARACTER(len=1), INTENT(in) :: ctyp
        INTEGER(intk), INTENT(in) :: igrid

        ! Local variables
        INTEGER(intk) :: i, j, k, icount
        INTEGER(intk) :: istart, istop, jstart, jstop, kstart, kstop
        INTEGER(intk) :: isteps, jsteps, ksteps, is, js, ks

        CALL start_and_stop_restrict(istart, istop, jstart, jstop, &
            kstart, kstop, ctyp, igrid)

        ! Computing number of steps in strided loop
        isteps = (istop - istart)/2 + 1
        jsteps = (jstop - jstart)/2 + 1
        ksteps = (kstop - kstart)/2 + 1

        DO i = istart, istop, 2
            DO j = jstart, jstop, 2
                DO k = kstart, kstop, 2
                    ! Compute the index of the coarse cell in the send buffer
                    ks = (k - kstart)/2
                    js = (j - jstart)/2
                    is = (i - istart)/2
                    icount = 1 + ks + js*ksteps + is*jsteps*ksteps
                    ! Assign one fine value to coarse cell buffer entry
                    sendbuf(icount) = ff(k, j, i)
                END DO
            END DO
        END DO

        ! Sanity check that the number of iterations is correct
        IF (icount /= isteps*jsteps*ksteps) CALL errr(__FILE__, __LINE__)

    END SUBROUTINE restrict_n


    SUBROUTINE restrict_e(kk, jj, ii, ff, sendbuf, ctyp, igrid)
        ! Subroutine arguments
        INTEGER(intk), INTENT(IN) :: kk, jj, ii
        REAL(realk), INTENT(IN) :: ff(kk, jj, ii)
        REAL(realk), CONTIGUOUS, INTENT(INOUT) :: sendbuf(:)
        CHARACTER(len=1), INTENT(in) :: ctyp
        INTEGER(intk), INTENT(in) :: igrid

        ! Local variables
        INTEGER(intk) :: i, j, k, icount, is, js, ks, isteps, jsteps, ksteps
        INTEGER(intk) :: istart, istop, jstart, jstop, kstart, kstop

        CALL start_and_stop_restrict(istart, istop, jstart, jstop, &
            kstart, kstop, ctyp, igrid)

        ! Computing number of steps in strided loop
        isteps = (istop - istart)/2 + 1
        jsteps = (jstop - jstart)/2 + 1
        ksteps = (kstop - kstart)/2 + 1

        DO i = istart, istop, 2
            DO j = jstart, jstop, 2
                DO k = kstart, kstop, 2
                    ! Compute the index of the coarse cell in the send buffer
                    ks = (k - kstart)/2
                    js = (j - jstart)/2
                    is = (i - istart)/2
                    icount = 1 + ks + js*ksteps + is*jsteps*ksteps
                    ! Very specific here...
                    sendbuf(icount) = MIN(1.0_realk, MAX(ff(k, j, i), &
                        ff(k, j, i+1), ff(k, j+1, i), ff(k, j+1, i+1), &
                        ff(k+1, j, i), ff(k+1, j, i+1), ff(k+1, j+1, i), &
                        ff(k+1, j+1, i+1)))
                END DO
            END DO
        END DO

        ! Sanity check that the number of iterations is correct
        IF (icount /= isteps*jsteps*ksteps) CALL errr(__FILE__, __LINE__)

    END SUBROUTINE restrict_e


    SUBROUTINE restrict_f(kk, jj, ii, ff, sendbuf, ctyp, igrid)
        ! Subroutine arguments
        INTEGER(intk), INTENT(IN) :: kk, jj, ii
        REAL(realk), INTENT(IN) :: ff(kk, jj, ii)
        REAL(realk), CONTIGUOUS, INTENT(INOUT) :: sendbuf(:)
        CHARACTER(len=1), INTENT(in) :: ctyp
        INTEGER(intk), INTENT(in) :: igrid

        ! Local variables
        INTEGER(intk) :: i, j, k, icount, is, js, ks, isteps, jsteps, ksteps
        INTEGER(intk) :: istart, istop, jstart, jstop, kstart, kstop

        CALL start_and_stop_restrict(istart, istop, jstart, jstop, &
            kstart, kstop, ctyp, igrid)

        ! Computing number of steps in strided loop
        ksteps = (kstop - kstart)/2 + 1
        jsteps = (jstop - jstart)/2 + 1
        isteps = (istop - istart)/2 + 1

        DO i = istart, istop, 2
            DO j = jstart, jstop, 2
                DO k = kstart, kstop, 2
                    ! Compute the index of the coarse cell in the send buffer
                    ks = (k - kstart)/2
                    js = (j - jstart)/2
                    is = (i - istart)/2
                    icount = 1 + ks + js*ksteps + is*jsteps*ksteps
                    ! Very specific here...
                    sendbuf(icount) = MIN(2.0_realk, MAX(ff(k, j, i), &
                        ff(k, j, i+1), ff(k, j+1, i), ff(k, j+1, i+1), &
                        ff(k+1, j, i), ff(k+1, j, i+1), ff(k+1, j+1, i), &
                        ff(k+1, j+1, i+1)))
                END DO
            END DO
        END DO

        ! Sanity check that the number of iterations is correct
        IF (icount /= isteps*jsteps*ksteps) CALL errr(__FILE__, __LINE__)

    END SUBROUTINE restrict_f


    SUBROUTINE restrict_i(kk, jj, ii, ff, sendbuf, ctyp, igrid)
        ! Subroutine arguments
        INTEGER(intk), INTENT(IN) :: kk, jj, ii
        REAL(realk), INTENT(IN) :: ff(kk, jj, ii)
        REAL(realk), CONTIGUOUS, INTENT(INOUT) :: sendbuf(:)
        CHARACTER(len=1), INTENT(in) :: ctyp
        INTEGER(intk), INTENT(in) :: igrid

        ! Local variables
        INTEGER(intk) :: i, j, k, icount, is, js, ks
        INTEGER(intk) :: istart, istop, jstart, jstop, kstart, kstop
        INTEGER(intk) :: isteps, jsteps, ksteps

        CALL start_and_stop_restrict(istart, istop, jstart, jstop, &
            kstart, kstop, ctyp, igrid)

        ! Computing number of steps in strided loop
        ksteps = (kstop - kstart)/2 + 1
        jsteps = (jstop - jstart)/2 + 1
        isteps = (istop - istart)/2 + 1

        DO i = istart, istop, 2
            DO j = jstart, jstop, 2
                DO k = kstart, kstop, 2
                    ! Compute the index of the coarse cell in the send buffer
                    ks = (k - kstart)/2
                    js = (j - jstart)/2
                    is = (i - istart)/2
                    icount = 1 + ks + js*ksteps + is*jsteps*ksteps
                    ! Assign maximum fine value to coarse cell buffer entry
                    sendbuf(icount) = MAX(ff(k, j, i), &
                        ff(k, j, i+1), &
                        ff(k, j+1, i), &
                        ff(k, j+1, i+1), &
                        ff(k+1, j, i), &
                        ff(k+1, j, i+1), &
                        ff(k+1, j+1, i), &
                        ff(k+1, j+1, i+1))
                END DO
            END DO
        END DO

        ! Sanity check that the number of iterations is correct
        IF (icount /= isteps*jsteps*ksteps) CALL errr(__FILE__, __LINE__)

    END SUBROUTINE restrict_i






    ! --- Helper subroutine for computing start and stop indices ---

    SUBROUTINE start_and_stop_restrict(ista, isto, jsta, jsto, ksta, ksto, &
        ctyp, igrid)

        ! Subroutine arguments
        INTEGER(intk), INTENT(out) :: ista, isto, jsta, jsto, ksta, ksto
        CHARACTER(len=1), INTENT(in) :: ctyp
        INTEGER(intk), INTENT(in) :: igrid

        ! Local variables
        INTEGER(intk) :: ii, jj, kk
        INTEGER(intk) :: istart, iend, jstart, jend, kstart, kend

        ! Initialize start and stop indices to default values
        istart = 3
        iend = 3
        jstart = 3
        jend = 3
        kstart = 3
        kend = 3

        ! Loop ranges are DO i = istart, ii-iend
        ! returned as DO i = ista, isto, 2
        SELECT CASE (ctyp)
        CASE ("U")
            istart = 2
            iend = 2
        CASE ("V")
            jstart = 2
            jend = 2
        CASE ("W")
            kstart = 2
            kend = 2
        CASE ("N")
            istart = 2
            iend = 2
            jstart = 2
            jend = 2
            kstart = 2
            kend = 2
        CASE ("E", "F", "P", "R", "S", 'I', 'T')
            ! No changes to start and stop indices
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

        IF (MOD(isto-ista, 2) > 0) CALL errr(__FILE__, __LINE__)
        IF (MOD(jsto-jsta, 2) > 0) CALL errr(__FILE__, __LINE__)
        IF (MOD(ksto-ksta, 2) > 0) CALL errr(__FILE__, __LINE__)

    END SUBROUTINE start_and_stop_restrict

END MODULE ftoc_restrict_mod