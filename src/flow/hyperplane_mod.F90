
MODULE hyperplane_mod

    USE core_mod
    USE realfield_mod, ONLY: get_grid3_linear
    USE intfield_mod, ONLY: get_grid3_ifk_linear
    USE fieldmapper_mod

    IMPLICIT NONE

    ! Infrastructure for hyperplane traversal
    TYPE(int_stencils_t), ALLOCATABLE, PROTECTED, TARGET :: mip_sip_list(:)
    !$omp declare target(mip_sip_list)
    TYPE(int_stencils_t), ALLOCATABLE, PROTECTED, TARGET :: idx_sip_list(:)
    !$omp declare target(idx_sip_list)

    TYPE(intfield_t) :: mip_hp_f
    !$omp declare target(mip_hp_f)
    TYPE(intfield_t) :: idx_hp_f
    !$omp declare target(idx_hp_f)

    PUBLIC :: hyperplane_init, hyperplane_finish, mip_hp_f, idx_hp_f

CONTAINS

    SUBROUTINE hyperplane_init()

        ! Local variables
        INTEGER(intk) :: i, igrid, kk, jj, ii
        REAL(realk), POINTER, CONTIGUOUS :: lw(:), ls(:), lb(:), lpr(:), &
            ue(:), un(:), ut(:)

        REAL(realk), POINTER, CONTIGUOUS :: lw_hp(:), ls_hp(:), lb_hp(:), &
            lpr_hp(:), ue_hp(:), un_hp(:), ut_hp(:)

        INTEGER(ifk), POINTER, CONTIGUOUS :: mip_hp(:), idx_hp(:)

        TYPE(field_t), POINTER :: lw_f, ls_f, lb_f, lpr_f, ue_f, un_f, ut_f
        TYPE(field_t), POINTER :: lw_hp_f, ls_hp_f, lb_hp_f, lpr_hp_f, &
            ue_hp_f, un_hp_f, ut_hp_f


        ! Allocation of the lists of stencils for all grids
        ALLOCATE(mip_sip_list(nmygrids))
        ALLOCATE(idx_sip_list(nmygrids))

        ! Declarign the Hyperlane fields
        ! Coefficients in [L] of ILU (including diagonal)
        CALL set_field("SIPLW_HP")
        CALL set_field("SIPLS_HP")
        CALL set_field("SIPLB_HP")
        CALL set_field("SIPLPR_HP")
        ! Coefficients in [U] of ILU
        CALL set_field("SIPUE_HP")
        CALL set_field("SIPUN_HP")
        CALL set_field("SIPUT_HP")

        ! Hyperplane traversal indices (preliminary allocation)
        CALL mip_hp_f%init("SIP_MIP_HP")
        CALL idx_hp_f%init("SIP_IDX_HP")

        ! Getting the original
        CALL get_field(lw_f, "SIPLW")
        CALL get_field(ls_f, "SIPLS")
        CALL get_field(lb_f, "SIPLB")
        CALL get_field(lpr_f, "SIPLPR")
        CALL get_field(ue_f, "SIPUE")
        CALL get_field(un_f, "SIPUN")
        CALL get_field(ut_f, "SIPUT")

        ! Getting the adapted field
        CALL get_field(lw_hp_f, "SIPLW_HP")
        CALL get_field(ls_hp_f, "SIPLS_HP")
        CALL get_field(lb_hp_f, "SIPLB_HP")
        CALL get_field(lpr_hp_f, "SIPLPR_HP")
        CALL get_field(ue_hp_f, "SIPUE_HP")
        CALL get_field(un_hp_f, "SIPUN_HP")
        CALL get_field(ut_hp_f, "SIPUT_HP")

        ! Iteration over all local grids
        DO i = 1, nmygrids

            igrid = mygrids(i)
            CALL get_mgdims(kk, jj, ii, igrid)

            ! Allocation of the arrays for one grid
            ALLOCATE(mip_sip_list(i)%arr(ii+jj+kk))
            ALLOCATE(idx_sip_list(i)%arr(ii*jj*kk))

            CALL get_grid3_ifk_linear(mip_hp, mip_hp_f, igrid)
            CALL get_grid3_ifk_linear(idx_hp, idx_hp_f, igrid)
            mip_hp = 0
            idx_hp = 1

            ! Populating the arrays for one grid
            CALL hyperplane_init_grid(kk, jj, ii, mip_hp, idx_hp)
            ! Sorting the indices within each hyperplane for one grid
            CALL hyperplane_sort_grid(kk, jj, ii, mip_hp, idx_hp)
            ! Checking the correctness of the arrays for one grid
            CALL hyperplane_check_grid(kk, jj, ii, mip_hp, idx_hp)

            ! ! Populating the arrays for one grid
            ! CALL hyperplane_init_grid(kk, jj, ii, mip_sip_list(i)%arr, &
            !     idx_sip_list(i)%arr)

            ! ! Sorting the indices within each hyperplane for one grid
            ! CALL hyperplane_sort_grid(kk, jj, ii, mip_sip_list(i)%arr, &
            !     idx_sip_list(i)%arr)

            ! ! Checking the correctness of the arrays for one grid
            ! CALL hyperplane_check_grid(kk, jj, ii, mip_sip_list(i)%arr, &
            !     idx_sip_list(i)%arr)
        END DO


        ! Iteration over all local grids
        DO i = 1, nmygrids

            igrid = mygrids(i)
            CALL get_mgdims(kk, jj, ii, igrid)

            ! Getting the already computed SIP coefficients for one grid
            CALL get_grid3_linear(lw, lw_f, igrid)
            CALL get_grid3_linear(ls, ls_f, igrid)
            CALL get_grid3_linear(lb, lb_f, igrid)
            CALL get_grid3_linear(lpr, lpr_f, igrid)
            CALL get_grid3_linear(ue, ue_f, igrid)
            CALL get_grid3_linear(un, un_f, igrid)
            CALL get_grid3_linear(ut, ut_f, igrid)

            CALL get_grid3_linear(lw_hp, lw_hp_f, igrid)
            CALL get_grid3_linear(ls_hp, ls_hp_f, igrid)
            CALL get_grid3_linear(lb_hp, lb_hp_f, igrid)
            CALL get_grid3_linear(lpr_hp, lpr_hp_f, igrid)
            CALL get_grid3_linear(ue_hp, ue_hp_f, igrid)
            CALL get_grid3_linear(un_hp, un_hp_f, igrid)
            CALL get_grid3_linear(ut_hp, ut_hp_f, igrid)

            CALL get_grid3_ifk_linear(mip_hp, mip_hp_f, igrid)
            CALL get_grid3_ifk_linear(idx_hp, idx_hp_f, igrid)

            ! Store the SIP coefficients for contiguous memory access
            CALL reorder_for_access(kk, jj, ii, lw_hp, lw, mip_hp, idx_hp)
            CALL reorder_for_access(kk, jj, ii, ls_hp, ls, mip_hp, idx_hp)
            CALL reorder_for_access(kk, jj, ii, lb_hp, lb, mip_hp, idx_hp)
            CALL reorder_for_access(kk, jj, ii, lpr_hp, lpr, mip_hp, idx_hp)
            CALL reorder_for_access(kk, jj, ii, ue_hp, ue, mip_hp, idx_hp)
            CALL reorder_for_access(kk, jj, ii, un_hp, un, mip_hp, idx_hp)
            CALL reorder_for_access(kk, jj, ii, ut_hp, ut, mip_hp, idx_hp)

        END DO

        !$omp target enter data map(mapper(map_intfield), always,to: mip_hp_f, idx_hp_f)

    END SUBROUTINE hyperplane_init


    PURE SUBROUTINE reorder_for_access(kk, jj, ii, acc, arr, mip, idxsip)

        ! Subroutine arguments
        INTEGER(intk), INTENT(in) :: kk, jj, ii
        REAL(realk), INTENT(inout) :: acc(kk*jj*ii)
        REAL(realk), INTENT(in) :: arr(kk*jj*ii)

        ! For the hyperplane traversal
        INTEGER(ifk), INTENT(in) :: mip(ii+jj+kk)
        INTEGER(ifk), INTENT(in) :: idxsip(ii*jj*kk)

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

        CALL mip_hp_f%finish()
        CALL idx_hp_f%finish()

    END SUBROUTINE hyperplane_finish


    ! ---------------------------------------------------------------------
    ! > Initialize index vectors for hyperplane traversal
    SUBROUTINE hyperplane_init_grid(kk, jj, ii, mip, idxsip)

        ! Subroutine arguments
        INTEGER(intk), INTENT(in) :: kk, jj, ii
        INTEGER(ifk), INTENT(inout) :: mip(ii+jj+kk), idxsip(ii*jj*kk)

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
        INTEGER(ifk), INTENT(in) :: mip(ii+jj+kk)
        INTEGER(ifk), INTENT(inout) :: idxsip(ii*jj*kk)

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
        INTEGER(ifk), INTENT(inout) :: vec(n)

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
        INTEGER(ifk), INTENT(in) :: mip(ii+jj+kk), idxsip(ii*jj*kk)

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