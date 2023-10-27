MODULE blockbt_mod
    USE core_mod
    USE ib_mod, ONLY: ctof, parent, ib

    IMPLICIT NONE(type, external)
    PRIVATE

    PUBLIC :: blockbt

CONTAINS
    SUBROUTINE blockbt(bt_f)
        ! This routine derives a blocking field BT for the scalar
        ! from the the blocking field BP for the pressure.
        !
        ! This implementation has a strong focus on robustness and is
        ! supposed to ensure mass conservation under all circumstances.
        ! There is surely space for refinement (of course not now...)

        ! Subroutine arguments
        TYPE(field_t), INTENT(inout) :: bt_f

        ! Local variables
        TYPE(field_t) :: hilf_f
        TYPE(field_t), POINTER :: bp_f
        REAL(realk), POINTER, CONTIGUOUS :: bt(:, :, :), bp(:, :, :), &
            hilf(:, :, :)
        INTEGER(intk) :: i, igrid, ilevel
        INTEGER(intk) :: kk, jj, ii
        INTEGER(intk) :: nfro, nbac, nrgt, nlft, nbot, ntop

        CALL get_field(bp_f, "BP")
        CALL hilf_f%init("HILF")

        ! Checking if any neighboring cell (shared face) has bp=1
        ! [only those cells can have a flux stencil associated with them]
        DO i = 1, nmygrids
            igrid = mygrids(i)
            CALL get_mgdims(kk, jj, ii, igrid)
            CALL bt_f%get_ptr(bt, igrid)
            CALL bp_f%get_ptr(bp, igrid)
            CALL blockbt_grid(kk, jj, ii, bt, bp)
        END DO

        DO ilevel = minlevel, maxlevel
            CALL connect(ilevel, 2, s1=bt_f, corners=.TRUE.)
        END DO

        ! Opening the cells near the parent boundary where the coarse grid
        ! allows a scalar flux
        DO ilevel = minlevel, maxlevel
            hilf_f%arr = -1.0
            CALL ctof(ilevel, hilf_f%arr, bt_f%arr)

            DO i = 1, nmygridslvl(ilevel)
                igrid = mygridslvl(i,ilevel)

                CALL get_mgbasb(nfro, nbac, nrgt, nlft, nbot, ntop, igrid)
                CALL get_mgdims(kk, jj, ii, igrid)
                CALL bt_f%get_ptr(bt, igrid)
                CALL hilf_f%get_ptr(hilf, igrid)

                CALL open_scalarflux_boundary(kk, jj, ii, bt, hilf, &
                    nfro, nbac, nrgt, nlft, nbot, ntop)
            END DO

            CALL connect(ilevel, 2, s1=bt_f, corners=.TRUE.)
        END DO

        CALL hilf_f%finish()
    END SUBROUTINE blockbt


    PURE SUBROUTINE blockbt_grid(kk, jj, ii, bt, bp)
        ! Subroutine arguments
        INTEGER(intk), INTENT(in) :: kk, jj, ii
        REAL(realk), INTENT(out) :: bt(kk, jj, ii)
        REAL(realk), INTENT(in) :: bp(kk, jj, ii)

        ! Local variables
        INTEGER(intk) :: i, j, k

        DO i = 2, ii-1
            DO j = 2, jj-1
                DO k = 2, kk-1
                    bt(k, j, i) = MAX(bp(k,j,i), bp(k,j,i+1), bp(k,j,i-1), &
                        bp(k,j+1,i), bp(k,j-1,i), bp(k+1,j,i), bp(k-1,j,i))
                END DO
            END DO
        END DO
    END SUBROUTINE blockbt_grid


    SUBROUTINE open_scalarflux_boundary(kk, jj, ii, bt, hilf, nfro, nbac, &
            nrgt, nlft, nbot, ntop)

        ! Subroutine arguments
        INTEGER(intk), INTENT(in) :: ii, jj, kk
        REAL(realk), INTENT(inout) :: bt(kk, jj, ii)
        REAL(realk), INTENT(in) :: hilf(kk, jj, ii)
        INTEGER(intk), INTENT(in) :: nfro, nbac, nrgt, nlft, nbot, ntop

        ! Local variables
        INTEGER(intk) :: i, j, k, i2, i3, j2, j3, k2, k3
        REAL(realk) :: c2, c3

        ! front = parent
        IF (nfro == 8) THEN
            i2 = 2
            i3 = 3
            DO j = 2, jj-1
                DO k = 2, kk-1
                    c2 = hilf(k, j, i2)
                    c3 = hilf(k, j, i3)
                    ! sanity check if correct ctof
                    IF (c2 < -0.5) CALL errr(__FILE__, __LINE__)
                    IF (c3 < -0.5) CALL errr(__FILE__, __LINE__)
                    ! if scalar flux on coarse, open receiver cells in fine
                    bt(k, j, i3) = MAX(c2*c3, bt(k, j, i3))
                END DO
            END DO
        END IF

        ! back = parent
        IF (nbac == 8) THEN
            i2 = ii-1
            i3 = ii-2
            DO j = 2, jj-1
                DO k = 2, kk-1
                    c2 = hilf(k, j, i2)
                    c3 = hilf(k, j, i3)
                    ! sanity check if correct ctof
                    IF (c2 < -0.5) CALL errr(__FILE__, __LINE__)
                    IF (c3 < -0.5) CALL errr(__FILE__, __LINE__)
                    ! if scalar flux on coarse, open receiver cells in fine
                    bt(k, j, i3) = MAX(c2*c3, bt(k, j, i3))
                END DO
            END DO
        END IF

        ! right = parent
        IF (nrgt == 8) THEN
            j2 = 2
            j3 = 3
            DO i = 2, ii-1
                DO k = 2, kk-1
                    c2 = hilf(k, j2, i)
                    c3 = hilf(k, j3, i)
                    ! sanity check if correct ctof
                    IF (c2 < -0.5) CALL errr(__FILE__, __LINE__)
                    IF (c3 < -0.5) CALL errr(__FILE__, __LINE__)
                    ! if scalar flux on coarse, open receiver cells in fine
                    bt(k, j3, i) = MAX(c2*c3, bt(k, j3, i))
                END DO
            END DO
        END IF

        ! left = parent
        IF (nlft == 8) THEN
            j2 = jj-1
            j3 = jj-2
            DO i = 2, ii-1
                DO k = 2, kk-1
                    c2 = hilf(k, j2, i)
                    c3 = hilf(k, j3, i)
                    ! sanity check if correct ctof
                    IF (c2 < -0.5) CALL errr(__FILE__, __LINE__)
                    IF (c3 < -0.5) CALL errr(__FILE__, __LINE__)
                    ! if scalar flux on coarse, open receiver cells in fine
                    bt(k, j3, i) = MAX(c2*c3, bt(k, j3, i))
                END DO
            END DO
        END IF

        ! bottom = parent
        IF (nbot == 8) THEN
            k2 = 2
            k3 = 3
            DO i = 2, ii-1
                DO j = 2, jj-1
                    c2 = hilf(k2, j, i)
                    c3 = hilf(k3, j, i)
                    ! sanity check i f correct ctof
                    IF (c2 < -0.5) CALL errr(__FILE__, __LINE__)
                    IF (c3 < -0.5) CALL errr(__FILE__, __LINE__)
                    ! if scalar flux on coarse, open receiver cells in fine
                    bt(k3, j, i) = MAX(c2*c3, bt(k3, j, i))
                END DO
            END DO
        END IF

        ! top = parent
        IF (ntop == 8) THEN
            k2 = kk-1
            k3 = kk-2
            DO i = 2, ii-1
                DO j = 2, jj-1
                    c2 = hilf(k2, j, i)
                    c3 = hilf(k3, j, i)
                    ! sanity check if correct ctof
                    IF (c2 < -0.5) CALL errr(__FILE__, __LINE__)
                    IF (c3 < -0.5) CALL errr(__FILE__, __LINE__)
                    ! if scalar flux on coarse, open receiver cells in fine
                    bt(k3, j, i) = MAX(c2*c3, bt(k3, j, i))
                END DO
            END DO
        END IF

    END SUBROUTINE open_scalarflux_boundary
END MODULE blockbt_mod
