
MODULE sip_hyperplane_mod

    USE, INTRINSIC :: iso_fortran_env
    USE, INTRINSIC :: ieee_arithmetic

    USE core_mod
    USE laplacephi_mod

    IMPLICIT NONE(type, external)

    PRIVATE

    TYPE(intfield_t), PROTECTED, TARGET :: mip_hp_f
    TYPE(intfield_t), PROTECTED, TARGET :: idx_hp_f

    PUBLIC :: sip_hyperplane_init, sip_hyperplane_finish, &
        sipiter1, sipiter2, mip_hp_f, idx_hp_f

CONTAINS

    SUBROUTINE sip_hyperplane_init()

        ! Local variables
        INTEGER(intk) :: i, igrid, kk, jj, ii
        REAL(realk), POINTER, CONTIGUOUS :: lw(:), ls(:), lb(:), lpr(:), &
            ue(:), un(:), ut(:), aw(:), ae(:), as(:), an(:), ab(:), at(:), &
            ap(:, :, :), bp(:, :, :)
        INTEGER(ifk), POINTER, CONTIGUOUS :: mip_hp(:), idx_hp(:)

        TYPE(field_t), POINTER :: lw_f, ls_f, lb_f, lpr_f, ue_f, un_f, ut_f
        TYPE(field_t), POINTER :: aw_f, ae_f, as_f, an_f, ab_f, at_f, ap_f, bp_f

        ! Setting up the hyperplane traversal infrastructure
        CALL mip_hp_f%init("SIP_MIP_HP")
        CALL idx_hp_f%init("SIP_IDX_HP")

        ! Iteration over all local grids
        DO i = 1, nmygrids

            igrid = mygrids(i)
            CALL get_mgdims(kk, jj, ii, igrid)
            CALL get_grid3_ifk_linear(mip_hp, mip_hp_f, igrid)
            CALL get_grid3_ifk_linear(idx_hp, idx_hp_f, igrid)
            mip_hp = 0
            idx_hp = -1

            ! Populating the integer arrays for one grid
            CALL sip_hyperplane_init_grid(kk, jj, ii, mip_hp, idx_hp)
            ! Sorting the indices within each hyperplane for one grid
            CALL sip_hyperplane_sort_grid(kk, jj, ii, mip_hp, idx_hp)
            ! Checking the correctness of the arrays for one grid
            CALL sip_hyperplane_check_grid(kk, jj, ii, mip_hp, idx_hp)
        END DO


        ! Getting the coefficient fields
        ! Coefficients in [L] of ILU (including diagonal)
        CALL get_field(lw_f, "SIPLW")
        CALL get_field(ls_f, "SIPLS")
        CALL get_field(lb_f, "SIPLB")
        CALL get_field(lpr_f, "SIPLPR")
        CALL get_field(ue_f, "SIPUE")
        CALL get_field(un_f, "SIPUN")
        CALL get_field(ut_f, "SIPUT")

        CALL get_field(aw_f, "GSAW")
        CALL get_field(ae_f, "GSAE")
        CALL get_field(as_f, "GSAS")
        CALL get_field(an_f, "GSAN")
        CALL get_field(ab_f, "GSAB")
        CALL get_field(at_f, "GSAT")
        CALL get_field(ap_f, "GSAP")
        CALL get_field(bp_f, "BP")


        ! Performing a hyperplane-parallel LU decomposition
        DO i = 1, nmygrids

            igrid = mygrids(i)
            CALL get_mgdims(kk, jj, ii, igrid)

            ! Getting the already computed SIP coefficients for one grid
            CALL get_grid3_real_linear(lw, lw_f, igrid)
            CALL get_grid3_real_linear(ls, ls_f, igrid)
            CALL get_grid3_real_linear(lb, lb_f, igrid)
            CALL get_grid3_real_linear(lpr, lpr_f, igrid)
            CALL get_grid3_real_linear(ue, ue_f, igrid)
            CALL get_grid3_real_linear(un, un_f, igrid)
            CALL get_grid3_real_linear(ut, ut_f, igrid)

            CALL get_grid1_real(aw, aw_f, igrid)
            CALL get_grid1_real(ae, ae_f, igrid)
            CALL get_grid1_real(as, as_f, igrid)
            CALL get_grid1_real(an, an_f, igrid)
            CALL get_grid1_real(ab, ab_f, igrid)
            CALL get_grid1_real(at, at_f, igrid)

            CALL get_grid3_real(ap, ap_f, igrid)
            CALL get_grid3_real(bp, bp_f, igrid)

            CALL get_grid3_ifk_linear(mip_hp, mip_hp_f, igrid)
            CALL get_grid3_ifk_linear(idx_hp, idx_hp_f, igrid)

            ! Hyperplane variants of the SIP LU decomposition for one grid
            CALL sip_hyperplane_lu(kk, jj, ii, lw, ls, lb, lpr, &
                ue, un, ut, aw, ae, as, an, ab, at, ap, bp, mip_hp, idx_hp)

        END DO

    END SUBROUTINE sip_hyperplane_init


    SUBROUTINE sip_hyperplane_lu(kk, jj, ii, lw, ls, lb, lpr, &
        ue, un, ut, aw, ae, as, an, ab, at, ap, bp, mip, idxsip)

        ! Subroutine arguments
        INTEGER(intk), INTENT(in) :: kk, jj, ii

        REAL(realk), INTENT(inout) :: lw(kk*jj*ii), ls(kk*jj*ii), &
            lb(kk*jj*ii), lpr(kk*jj*ii), ue(kk*jj*ii), un(kk*jj*ii), &
            ut(kk*jj*ii)

        REAL(realk), INTENT(in) :: ae(ii), an(jj), at(kk), aw(ii), as(jj), &
            ab(kk), ap(kk, jj, ii), bp(kk, jj, ii)

        INTEGER(ifk), INTENT(in) :: mip(ii+jj+kk)
        INTEGER(ifk), INTENT(in) :: idxsip(ii*jj*kk)

        ! Local variables
        REAL(realk), PARAMETER :: alfa = 0.92
        REAL(realk) :: p1, p2, p3
        INTEGER(intk) :: n3dmin, n3dmax, m, lm, lp, ip, iacc, &
            idx, idx_km, idx_jm, idx_im, i, j, k

        REAL(realk), ALLOCATABLE :: lw_tmp(:), ls_tmp(:), &
            lb_tmp(:), lpr_tmp(:), ue_tmp(:), un_tmp(:), ut_tmp(:)


        ! Allocating temporary arrays
        ALLOCATE(lw_tmp(kk*jj*ii), ls_tmp(kk*jj*ii), lb_tmp(kk*jj*ii), &
            lpr_tmp(kk*jj*ii), ue_tmp(kk*jj*ii), un_tmp(kk*jj*ii), &
            ut_tmp(kk*jj*ii))

        lw_tmp = 0.0_realk
        ls_tmp = 0.0_realk
        lb_tmp = 0.0_realk
        lpr_tmp = 0.0_realk
        ue_tmp = 0.0_realk
        un_tmp = 0.0_realk
        ut_tmp = 0.0_realk

        n3dmin = 3 + 3 + 3
        n3dmax = (ii-2) + (jj-2) + (kk-2)
        iacc = -1

        ! Iterating over the hyperplanes h(k, j, i) = m
        DO m = n3dmin, n3dmax

            lm = INT(mip(m), intk)
            lp = INT(mip(m+1), intk) - lm

            ! > Parallel operations within the hyperplane h(k, j, i) = m
            DO ip = 1, lp

                ! Computing the contiguous access index
                iacc = lm + ip

                ! Computing the required indices
                idx = INT(idxsip(iacc), intk)
                idx_km = idx - 1
                idx_jm = idx - kk
                idx_im = idx - (kk*jj)
                CALL ind2sub(idx, k, j, i, kk, jj, ii)

                lw_tmp(idx) = aw(i) * bu(k, j, i-1) &
                    / (1.0 + alfa * (un_tmp(idx_im) + ut_tmp(idx_im)))
                ! lw(k, j, i) = aw(i)*bu(k, j, i-1) &
                !     /(1.0 + alfa*(un(k, j, i-1) + ut(k, j, i-1)))

                ls_tmp(idx) = as(j) * bv(k, j-1, i) &
                    / (1.0 + alfa * (ue_tmp(idx_jm) + ut_tmp(idx_jm)))
                ! ls(k, j, i) = as(j)*bv(k, j-1, i) &
                !     /(1.0 + alfa*(ue(k, j-1, i) + ut(k, j-1, i)))

                lb_tmp(idx) = ab(k) * bw(k-1, j, i) &
                    / (1.0 + alfa * (un_tmp(idx_km) + ue_tmp(idx_km)))
                ! lb(k, j, i) = ab(k)*bw(k-1, j, i) &
                !     /(1.0 + alfa*(un(k-1, j, i) + ue(k-1, j, i)))

                p1 = alfa * (lb_tmp(idx) * ue_tmp(idx_km) &
                    + ls_tmp(idx) * ue_tmp(idx_jm))
                ! p1 = alfa*(lb(k, j, i)*ue(k-1, j, i) &
                !   + ls(k, j, i)*ue(k, j-1, i))

                p2 = alfa * (lb_tmp(idx) * un_tmp(idx_km) &
                    + lw_tmp(idx) * un_tmp(idx_im))
                ! p2 = alfa*(lb(k, j, i)*un(k-1, j, i) &
                !     + lw(k, j, i)*un(k, j, i-1))

                p3 = alfa * (lw_tmp(idx) * ut_tmp(idx_im) &
                    + ls_tmp(idx) * ut_tmp(idx_jm))
                ! p3 = alfa*(lw(k, j, i)*ut(k, j, i-1) &
                !     + ls(k, j, i)*ut(k, j-1, i))

                lpr_tmp(idx) = 1.0 / (ap(k, j, i) + p1 + p2 + p3 &
                    - lb_tmp(idx) * ut_tmp(idx_km) &
                    - lw_tmp(idx) * ue_tmp(idx_im) &
                    - ls_tmp(idx) * un_tmp(idx_jm) &
                    + 1.0e-20)
                ! lpr(k, j, i) = 1.0/(ap(k, j, i) + p1 + p2 + p3 &
                !     - lb(k, j, i)*ut(k-1, j, i) &
                !     - lw(k, j, i)*ue(k, j, i-1) &
                !     - ls(k, j, i)*un(k, j-1, i) &
                !     + 1.0e-20)

                ue_tmp(idx) = (ae(i) * bu(k, j, i) - p1) * lpr_tmp(idx)
                un_tmp(idx) = (an(j) * bv(k, j, i) - p2) * lpr_tmp(idx)
                ut_tmp(idx) = (at(k) * bw(k, j, i) - p3) * lpr_tmp(idx)
                ! ue(k, j, i) = (ae(i)*bu(k, j, i) - p1)*lpr(k, j, i)
                ! un(k, j, i) = (an(j)*bv(k, j, i) - p2)*lpr(k, j, i)
                ! ut(k, j, i) = (at(k)*bw(k, j, i) - p3)*lpr(k, j, i)

                lw_tmp(idx) = lw_tmp(idx) * lpr_tmp(idx)
                ls_tmp(idx) = ls_tmp(idx) * lpr_tmp(idx)
                lb_tmp(idx) = lb_tmp(idx) * lpr_tmp(idx)
                ! lw(k, j, i) = lw(k, j, i)*lpr(k, j, i)
                ! ls(k, j, i) = ls(k, j, i)*lpr(k, j, i)
                ! lb(k, j, i) = lb(k, j, i)*lpr(k, j, i)

                ! Feeding computed values into access-optimized arrays
                lw(iacc) = lw_tmp(idx)
                ls(iacc) = ls_tmp(idx)
                lb(iacc) = lb_tmp(idx)
                lpr(iacc) = lpr_tmp(idx)
                ue(iacc) = ue_tmp(idx)
                un(iacc) = un_tmp(idx)
                ut(iacc) = ut_tmp(idx)

            END DO
        END DO

        ! Deallocating temporary arrays
        DEALLOCATE(lw_tmp, ls_tmp, lb_tmp, lpr_tmp, ue_tmp, un_tmp, ut_tmp)

    CONTAINS

        ! These routines are not protected against out-of-bounds and not valid
        ! for the edges ii, jj, kk - but in the scope above that is OK
        PURE REAL(realk) FUNCTION bu(k, j, i)
            INTEGER(intk), INTENT(in) :: k, j, i
            bu = bp(k, j, i)*bp(k, j, i+1)
        END FUNCTION bu

        PURE REAL(realk) FUNCTION bv(k, j, i)
            INTEGER(intk), INTENT(in) :: k, j, i
            bv = bp(k, j, i)*bp(k, j+1, i)
        END FUNCTION bv

        PURE REAL(realk) FUNCTION bw(k, j, i)
            INTEGER(intk), INTENT(in) :: k, j, i
            bw = bp(k, j, i)*bp(k+1, j, i)
        END FUNCTION bw

    END SUBROUTINE sip_hyperplane_lu


    SUBROUTINE sip_hyperplane_finish()

        CALL mip_hp_f%finish()
        CALL idx_hp_f%finish()

    END SUBROUTINE sip_hyperplane_finish


    ! ---------------------------------------------------------------------
    ! > Initialize index vectors for hyperplane traversal
    SUBROUTINE sip_hyperplane_init_grid(kk, jj, ii, mip, idxsip)

        ! Subroutine arguments
        INTEGER(intk), INTENT(in) :: kk, jj, ii
        INTEGER(ifk), INTENT(inout) :: mip(ii+jj+kk), idxsip(ii*jj*kk)

        ! Local variables
        INTEGER(intk) :: n3dmin, n3dmax
        INTEGER(intk) :: nkmin, nkmax, ic1, nimin, nimax
        INTEGER(intk) :: m, k, j, i, idx
        INTEGER(ifk) :: lm

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

    END SUBROUTINE sip_hyperplane_init_grid


    ! ---------------------------------------------------------------------
    ! > Sort indices within each hyperplane
    SUBROUTINE sip_hyperplane_sort_grid(kk, jj, ii, mip, idxsip)

        ! Subroutine arguments
        INTEGER(intk), INTENT(in) :: kk, jj, ii
        INTEGER(ifk), INTENT(in) :: mip(ii+jj+kk)
        INTEGER(ifk), INTENT(inout) :: idxsip(ii*jj*kk)

        ! Local variables
        INTEGER(intk) :: n3dmin, n3dmax
        INTEGER(ifk) :: len, m, s, e

        n3dmin = 3 + 3 + 3
        n3dmax = (ii-2) + (jj-2) + (kk-2)

        ! Iterate over hyperplanes H(k,j,i)=m and sorting within each plane
        DO m = n3dmin, n3dmax
            len = mip(m+1) - mip(m)
            s = mip(m) + 1
            e = mip(m+1)
            CALL sip_hyperplane_sort(len, idxsip(s:e))
        END DO

    END SUBROUTINE sip_hyperplane_sort_grid


    ! ---------------------------------------------------------------------
    ! > Simple selection sort (not performant, not for frequent execution)
    SUBROUTINE sip_hyperplane_sort(n, vec)

        ! Subroutine arguments
        INTEGER(ifk), INTENT(in) :: n
        INTEGER(ifk), INTENT(inout) :: vec(n)

        ! Local variables
        INTEGER(ifk) :: ix, v_rand(n), v_asc(n)
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

    END SUBROUTINE sip_hyperplane_sort


    ! ---------------------------------------------------------------------
    ! > Check index vector correctness
    SUBROUTINE sip_hyperplane_check_grid(kk, jj, ii, mip, idxsip)

        ! Subroutine arguments
        INTEGER(intk), INTENT(in) :: kk, jj, ii
        INTEGER(ifk), INTENT(in) :: mip(ii+jj+kk), idxsip(ii*jj*kk)

        ! Local variables
        INTEGER(intk) :: k, j, i, m, n3dmin, n3dmax, idx
        INTEGER(intk), ALLOCATABLE :: test(:, :, :)
        INTEGER(ifk) :: len, lm, ip

        ! Creating test grid
        ALLOCATE(test(kk, jj, ii), source=0_intk)

        ! Visiting  indicies within (i, j, k) in [3,ii-2] x [3,jj-2] x [3,kk-2]
        ! and incrementing the corresponding test grid point by 1

        n3dmin = 3 + 3 + 3
        n3dmax = (ii-2) + (jj-2) + (kk-2)

        DO m = n3dmin, n3dmax
            len = mip(m+1) - mip(m)
            lm  = mip(m)
            DO ip = 1, len
                idx = INT(idxsip(lm + ip), kind=intk)
                CALL ind2sub(idx, k, j, i, kk, jj, ii)
                test(k, j, i) = test(k, j, i) + 1
                IF (k+j+i /= m) THEN
                    WRITE(*, *) 'hyperplane_check_grid: index error'
                    CALL errr(__FILE__, __LINE__)
                END IF
            END DO
        END DO

        ! Performing tests for correct coverage and deallocating test grid
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

    END SUBROUTINE sip_hyperplane_check_grid


    SUBROUTINE sipiter1(kk, jj, ii, rhs, res, lw, ls, lb, lpr, mip, idxsip)

        ! Subroutine arguments
        INTEGER(intk), INTENT(in) :: kk, jj, ii
        REAL(realk), INTENT(inout) :: res(kk*jj*ii)
        REAL(realk), INTENT(in) :: rhs(kk*jj*ii), lw(kk*jj*ii), ls(kk*jj*ii), &
            lb(kk*jj*ii), lpr(kk*jj*ii)

        ! For the hyperplane traversal
        INTEGER(ifk), INTENT(in) :: mip(ii+jj+kk)
        INTEGER(ifk), INTENT(in) :: idxsip(ii*jj*kk)

        ! Local variables
        INTEGER(intk) :: n3dmin, n3dmax, m, lm, lp, ip, iacc, &
            idx, idx_km, idx_jm, idx_im

        ! Subroutine body
        n3dmin = 3 + 3 + 3
        n3dmax = (ii-2) + (jj-2) + (kk-2)
        iacc = -1

        ! Iterating over the hyperplanes H(k, j, i) = m
        DO m = n3dmin, n3dmax

            lm = INT(mip(m), intk)
            lp = INT(mip(m+1), intk) - lm

            ! > Parallel operations on the hyperplane (k, j, i) = m
            DO ip = 1, lp

                ! Computing the contiguous access index
                iacc = lm + ip

                ! Computing the required indices
                idx = INT(idxsip(iacc), intk)
                idx_km = idx - 1
                idx_jm = idx - kk
                idx_im = idx - (kk*jj)

                ! Accounting for RHS
                res(idx) = (rhs(idx) + res(idx)) * lpr(iacc)

                ! Performing the forward substitution
                res(idx) = res(idx) - lb(iacc)*res(idx_km) - &
                    ls(iacc)*res(idx_jm) - lw(iacc)*res(idx_im)

                ! Original code:
                ! res(k, j, i) = res(k, j, i) - lw(k, j, i)*res(k, j, i-1) &
                !     - ls(k, j, i)*res(k, j-1, i) - lb(k, j, i)*res(k-1, j, i)

            END DO
        END DO

    END SUBROUTINE sipiter1


    SUBROUTINE sipiter2(kk, jj, ii, phi, res, ue, un, ut, mip, idxsip)

        ! Subroutine arguments
        INTEGER(intk), INTENT(in) :: kk, jj, ii
        REAL(realk), INTENT(inout) :: phi(kk*jj*ii), res(kk*jj*ii)
        REAL(realk), INTENT(in) :: ue(kk*jj*ii), un(kk*jj*ii), ut(kk*jj*ii)

        ! For the hyperplane traversal
        INTEGER(ifk), INTENT(in) :: mip(ii+jj+kk)
        INTEGER(ifk), INTENT(in) :: idxsip(ii*jj*kk)

        ! Local variables
        INTEGER(intk) :: n3dmin, n3dmax, m, lm, lp, ip, iacc, &
            idx, idx_kp, idx_jp, idx_ip, i, j, k

        ! Subroutine body
        n3dmin = 3 + 3 + 3
        n3dmax = (ii-2) + (jj-2) + (kk-2)
        iacc = -1

        ! Iterating (REVERSE) over the hyperplanes H(k, j, i) = m
        DO m = n3dmax, n3dmin, -1

            lm = INT(mip(m), intk)
            lp = INT(mip(m+1), intk) - lm

            ! > Parallel operations on the hyperplane (k, j, i) = m
            DO ip = 1, lp

                ! Computing the contiguous access index
                iacc = lm + ip

                ! Computing the required indices
                idx = INT(idxsip(iacc), intk)
                idx_kp = idx + 1
                idx_jp = idx + kk
                idx_ip = idx + (kk*jj)

                ! Performing the backward substitution
                res(idx) = res(idx) - ut(iacc)*res(idx_kp) - &
                    un(iacc)*res(idx_jp) - ue(iacc)*res(idx_ip)

                ! Original code:
                ! res(k, j, i) = res(k, j, i) - un(k, j, i)*res(k, j+1, i) &
                !     - ue(k, j, i)*res(k, j, i+1) - ut(k, j, i)*res(k+1, j, i)

            END DO
        END DO

        DO i = 3, ii-2
            DO j = 3, jj-2
                DO k = 3, kk-2
                    idx = k + (j-1)*kk + (i-1)*kk*jj
                    phi(idx) = phi(idx) + res(idx)
                END DO
            END DO
        END DO

    END SUBROUTINE sipiter2


END MODULE sip_hyperplane_mod
