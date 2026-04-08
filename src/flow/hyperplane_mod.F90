
MODULE hyperplane_mod

    USE core_mod

    IMPLICIT NONE

    ! Infrastructure for hyperplane traversal
    TYPE(int_stencils_t), ALLOCATABLE, PROTECTED, TARGET :: mip_sip_list(:)
    TYPE(int_stencils_t), ALLOCATABLE, PROTECTED, TARGET :: idx_sip_list(:)

    ! Coefficients in [L] of ILU (including diagonal)
    TYPE(real_stencils_t), ALLOCATABLE, PROTECTED, TARGET :: lw_sip_list(:)
    TYPE(real_stencils_t), ALLOCATABLE, PROTECTED, TARGET :: ls_sip_list(:)
    TYPE(real_stencils_t), ALLOCATABLE, PROTECTED, TARGET :: lb_sip_list(:)
    TYPE(real_stencils_t), ALLOCATABLE, PROTECTED, TARGET :: lpr_sip_list(:)

    ! Coefficients in [L] of ILU
    TYPE(real_stencils_t), ALLOCATABLE, PROTECTED, TARGET :: ue_sip_list(:)
    TYPE(real_stencils_t), ALLOCATABLE, PROTECTED, TARGET :: un_sip_list(:)
    TYPE(real_stencils_t), ALLOCATABLE, PROTECTED, TARGET :: ut_sip_list(:)

    PUBLIC :: hyperplane_init, hyperplane_finish

CONTAINS

    SUBROUTINE hyperplane_init()

        ! Local variables
        INTEGER(intk) :: i, igrid, kk, jj, ii
        REAL(realk), POINTER, CONTIGUOUS :: lw(:), ls(:), lb(:), lpr(:)
        REAL(realk), POINTER, CONTIGUOUS :: ue(:), un(:), ut(:)

        ! Allocation of the lists of stencils for all grids
        ALLOCATE(mip_sip_list(nmygrids))
        ALLOCATE(idx_sip_list(nmygrids))

        ! Iteration over all local grids
        DO i = 1, nmygrids

            igrid = mygrids(i)
            CALL get_mgdims(kk, jj, ii, igrid)

            ! Allocation of the arrays for one grid
            ALLOCATE(mip_sip_list(i)%arr(ii+jj+kk))
            ALLOCATE(idx_sip_list(i)%arr(ii*jj*kk))

            ! Populating the arrays for one grid
            CALL hyperplane_init_grid(kk, jj, ii, mip_sip_list(i)%arr, &
                idx_sip_list(i)%arr)

            ! Sorting the indices within each hyperplane for one grid
            CALL hyperplane_sort_grid(kk, jj, ii, mip_sip_list(i)%arr, &
                idx_sip_list(i)%arr)

            ! Checking the correctness of the arrays for one grid
            CALL hyperplane_check_grid(kk, jj, ii, mip_sip_list(i)%arr, &
                idx_sip_list(i)%arr)
        END DO

        ! Allocation of the lists of stencils for all grids
        ALLOCATE(lw_sip_list(nmygrids))
        ALLOCATE(ls_sip_list(nmygrids))
        ALLOCATE(lb_sip_list(nmygrids))
        ALLOCATE(lpr_sip_list(nmygrids))

        ALLOCATE(ue_sip_list(nmygrids))
        ALLOCATE(un_sip_list(nmygrids))
        ALLOCATE(ut_sip_list(nmygrids))

        ! Iteration over all local grids
        DO i = 1, nmygrids

            igrid = mygrids(i)
            CALL get_mgdims(kk, jj, ii, igrid)

            ! Getting the already computed SIP coefficients for one grid
            CALL get_fieldptr(lw, "SIPLW", igrid)
            CALL get_fieldptr(ls, "SIPLS", igrid)
            CALL get_fieldptr(lb, "SIPLB", igrid)
            CALL get_fieldptr(lpr, "SIPLPR", igrid)
            CALL get_fieldptr(ue, "SIPUE", igrid)
            CALL get_fieldptr(un, "SIPUN", igrid)
            CALL get_fieldptr(ut, "SIPUT", igrid)

            ALLOCATE(lw_sip_list(i)%arr(ii*jj*kk))
            ALLOCATE(ls_sip_list(i)%arr(ii*jj*kk))
            ALLOCATE(lb_sip_list(i)%arr(ii*jj*kk))
            ALLOCATE(lpr_sip_list(i)%arr(ii*jj*kk))
            ALLOCATE(ue_sip_list(i)%arr(ii*jj*kk))
            ALLOCATE(un_sip_list(i)%arr(ii*jj*kk))
            ALLOCATE(ut_sip_list(i)%arr(ii*jj*kk))

            ! Store the SIP coefficients for contiguous memory access
            CALL reorder_for_access(kk, jj, ii, lw_sip_list(i)%arr, lw, &
                mip_sip_list(i)%arr, idx_sip_list(i)%arr)
            CALL reorder_for_access(kk, jj, ii, ls_sip_list(i)%arr, ls, &
                mip_sip_list(i)%arr, idx_sip_list(i)%arr)
            CALL reorder_for_access(kk, jj, ii, lb_sip_list(i)%arr, lb, &
                mip_sip_list(i)%arr, idx_sip_list(i)%arr)
            CALL reorder_for_access(kk, jj, ii, lpr_sip_list(i)%arr, lpr, &
                mip_sip_list(i)%arr, idx_sip_list(i)%arr)
            CALL reorder_for_access(kk, jj, ii, ue_sip_list(i)%arr, ue, &
                mip_sip_list(i)%arr, idx_sip_list(i)%arr)
            CALL reorder_for_access(kk, jj, ii, un_sip_list(i)%arr, un, &
                mip_sip_list(i)%arr, idx_sip_list(i)%arr)
            CALL reorder_for_access(kk, jj, ii, ut_sip_list(i)%arr, ut, &
                mip_sip_list(i)%arr, idx_sip_list(i)%arr)

        END DO

    END SUBROUTINE hyperplane_init


    PURE SUBROUTINE reorder_for_access(kk, jj, ii, acc, arr, mip, idxsip)

        ! Subroutine arguments
        INTEGER(intk), INTENT(in) :: kk, jj, ii
        REAL(realk), INTENT(inout) :: acc(kk*jj*ii)
        REAL(realk), INTENT(in) :: arr(kk*jj*ii)

        ! For the hyperplane traversal
        INTEGER(intk), INTENT(in) :: mip(ii+jj+kk)
        INTEGER(intk), INTENT(in) :: idxsip(ii*jj*kk)

        ! Local variables
        INTEGER(intk) :: n3dmin, n3dmax, m, lm, len, ip, idx, iacc

        ! Subroutine body
        n3dmin = 3 + 3 + 3
        n3dmax = (ii-2) + (jj-2) + (kk-2)
        acc(kk*jj*ii) = -HUGE(1.0_realk)

        ! Iterating over the hyperplanes H(k, j, i) = m
        DO m = n3dmin, n3dmax

            lm = mip(m)
            len = mip(m+1) - lm

            ! > Parallel operations on the hyperplane (k, j, i) = m
            DO ip = 1, len

                iacc = lm + ip
                idx = idxsip(iacc)
                acc(iacc) = arr(idx)

            END DO
        END DO

    END SUBROUTINE reorder_for_access



    SUBROUTINE hyperplane_finish()

        ! Local variables
        INTEGER(intk) :: i

        ! Iteration over all local grids
        DO i = 1, nmygrids
            ! Deallocation of the arrays for one grid
            DEALLOCATE(mip_sip_list(i)%arr)
            DEALLOCATE(idx_sip_list(i)%arr)
            DEALLOCATE(lw_sip_list(i)%arr)
            DEALLOCATE(ls_sip_list(i)%arr)
            DEALLOCATE(lb_sip_list(i)%arr)
            DEALLOCATE(lpr_sip_list(i)%arr)
            DEALLOCATE(ue_sip_list(i)%arr)
            DEALLOCATE(un_sip_list(i)%arr)
            DEALLOCATE(ut_sip_list(i)%arr)
        END DO

        ! Deallocation of the lists of stencils for all grids
        DEALLOCATE(mip_sip_list)
        DEALLOCATE(idx_sip_list)
        DEALLOCATE(lw_sip_list)
        DEALLOCATE(ls_sip_list)
        DEALLOCATE(lb_sip_list)
        DEALLOCATE(lpr_sip_list)
        DEALLOCATE(ue_sip_list)
        DEALLOCATE(un_sip_list)
        DEALLOCATE(ut_sip_list)

    END SUBROUTINE hyperplane_finish


    ! ---------------------------------------------------------------------
    ! > Initialize index vectors for hyperplane traversal
    SUBROUTINE hyperplane_init_grid(kk, jj, ii, mip, idxsip)

        ! Subroutine arguments
        INTEGER(intk), INTENT(in) :: kk, jj, ii
        INTEGER(intk), INTENT(inout) :: mip(ii+jj+kk), idxsip(ii*jj*kk)

        ! Local variables
        INTEGER(intk) :: n3dmin, n3dmax
        INTEGER(intk) :: lm, nkmin, nkmax, ic1, nimin, nimax
        INTEGER(intk) :: m, k, j, i, idx

        ! Only cover [i=3,ii-2], [j=3,jj-2], [k=3,kk-2]
        n3dmin = 3 + 3 + 3
        n3dmax = (ii-2) + (jj-2) + (kk-2)

        ! Initialization (to facilitate error detection)
        mip = 0
        idxsip = -1

        ! Iterate over hyperplanes H(k,j,i)=m
        DO m = n3dmin, n3dmax

            ic1 = 0
            lm  = mip(m)
            nkmin = MAX(m - (ii-2) - (jj-2), 3)
            nkmax = MIN((m - 2*3), (kk-2))

            DO k = nkmin, nkmax
                nimin = MAX(m - (jj-2) - k, 3)
                nimax = MIN(m - k - 1*3, (ii-2))

                DO i = nimin, nimax

                    ! Calculate remaining index
                    j = m - k - i

                    ! Check bounds
                    IF (i < 3 .OR. i > (ii-2)) THEN
                        WRITE(*, *) 'hyperplane_init_grid: i out of bounds', i
                        CALL errr(__FILE__, __LINE__)
                    END IF
                    IF (j < 3 .OR. j > (jj-2)) THEN
                        WRITE(*, *) 'hyperplane_init_grid: j out of bounds', j
                        CALL errr(__FILE__, __LINE__)
                    END IF
                    IF (k < 3 .OR. k > (kk-2)) THEN
                        WRITE(*, *) 'hyperplane_init_grid: k out of bounds', k
                        CALL errr(__FILE__, __LINE__)
                    END IF

                    ! Increment counter
                    ic1 = ic1 + 1

                    ! Store linear index
                    CALL sub2ind(idx, k, j, i, kk, jj, ii)
                    idxsip(lm + ic1) = idx

                END DO
            END DO

            mip(m+1) = mip(m) + ic1

        END DO

    END SUBROUTINE hyperplane_init_grid


    ! ---------------------------------------------------------------------
    ! > Sort indices within each hyperplane
    SUBROUTINE hyperplane_sort_grid(kk, jj, ii, mip, idxsip)

        ! Subroutine arguments
        INTEGER(intk), INTENT(in) :: kk, jj, ii
        INTEGER(intk), INTENT(in) :: mip(ii+jj+kk)
        INTEGER(intk), INTENT(inout) :: idxsip(ii*jj*kk)

        ! Local variables
        INTEGER(intk) :: n3dmin, n3dmax, len, m, s, e

        n3dmin = 3 + 3 + 3
        n3dmax = (ii-2) + (jj-2) + (kk-2)

        ! Iterate over hyperplanes H(k,j,i)=m and sorting within each plane
        DO m = n3dmin, n3dmax
            len = mip(m+1) - mip(m)
            s   = mip(m) + 1
            e   = mip(m+1)
            CALL sort(len, idxsip(s:e))
        END DO

    END SUBROUTINE hyperplane_sort_grid


    ! ---------------------------------------------------------------------
    ! > Simple selection sort (for testing -- not performant)
    SUBROUTINE sort(n, vec)

        ! Subroutine arguments
        INTEGER(intk), INTENT(in) :: n
        INTEGER(intk), INTENT(inout) :: vec(n)

        ! Local variables
        INTEGER(intk) :: ix, v_rand(n), v_asc(n)
        LOGICAL :: mask(n)

        v_rand(:) = vec(:)
        mask = .TRUE.

        DO ix = 1, n
            v_asc(ix) = MINVAL(v_rand, mask)
            mask(MINLOC(v_rand, mask)) = .FALSE.
        END DO

        DO ix = 1, n-1
            IF (v_asc(ix+1) <= v_asc(ix)) THEN
                WRITE(*, *) 'sort: sorting failed'
                CALL errr(__FILE__, __LINE__)
            END IF
        END DO

        vec(:) = v_asc(:)

    END SUBROUTINE sort


    ! ---------------------------------------------------------------------
    ! > Check index vector correctness
    SUBROUTINE hyperplane_check_grid(kk, jj, ii, mip, idxsip)

        ! Subroutine arguments
        INTEGER(intk), INTENT(in) :: kk, jj, ii
        INTEGER(intk), INTENT(in) :: mip(ii+jj+kk), idxsip(ii*jj*kk)

        ! Local variables
        INTEGER(intk) :: k, j, i, len, m, idx, n3dmin, n3dmax, lm, ip
        INTEGER(intk), ALLOCATABLE :: test(:, :, :)

        ! Creating test grid
        ALLOCATE(test(kk, jj, ii))
        test = 0

        n3dmin = 3 + 3 + 3
        n3dmax = (ii-2) + (jj-2) + (kk-2)

        DO m = n3dmin, n3dmax
            len = mip(m+1) - mip(m)
            lm  = mip(m)
            DO ip = 1, len
                idx = idxsip(lm + ip)
                CALL ind2sub(idx, k, j, i, kk, jj, ii)
                test(k, j, i) = test(k, j, i) + 1
                IF (k+j+i /= m) THEN
                    WRITE(*, *) 'hyperplane_check_grid: index error'
                    CALL errr(__FILE__, __LINE__)
                END IF
            END DO
        END DO

        ! Performing tests and deallocating test grid
        IF (ANY(test > 1)) THEN
            WRITE(*, *) 'hyperplane_check_grid: some indices are duplicated'
            CALL errr(__FILE__, __LINE__)
        END IF
        IF (ANY(test(3:kk-2, 3:jj-2, 3:ii-2) < 1)) THEN
            WRITE(*, *) 'hyperplane_check_grid: some indices are missing'
            CALL errr(__FILE__, __LINE__)
        END IF
        IF (SUM(test) /= (kk-4)*(jj-4)*(ii-4)) THEN
            WRITE(*, *) 'hyperplane_check_grid: summation failed'
            CALL errr(__FILE__, __LINE__)
        END IF
        DEALLOCATE(test)

    END SUBROUTINE hyperplane_check_grid

END MODULE hyperplane_mod