MODULE sor_mod

    USE core_mod

    IMPLICIT NONE(type, external)

    PRIVATE

    LOGICAL :: is_init = .FALSE.

    ! Relaxation factor for SOR
    REAL(realk) :: omg

    PUBLIC :: sor_init, sor, sor_finish

CONTAINS

    SUBROUTINE sor_init(omg_in)
        ! Subroutine arguments
        REAL(realk), INTENT(in) :: omg_in

        ! Local variable
        INTEGER(intk) :: igr, igrid
        INTEGER(intk) :: k, j, i
        INTEGER(intk) :: kk, jj, ii
        REAL(realk), POINTER, CONTIGUOUS :: ap(:, :, :), rap(:, :, :)

        omg = omg_in

        CALL set_field("SOR_RAP")

        DO igr = 1, nmygrids
            igrid = mygrids(igr)

            CALL get_mgdims(kk, jj, ii, igrid)
            CALL get_fieldptr(rap, "SOR_RAP", igrid)
            CALL get_fieldptr(ap, "GSAP", igrid)
            DO i = 3, ii-2
                DO j = 3, jj-2
                    DO k = 3, kk-2
                        rap(k, j, i) = divide0(1.0_realk, ap(k, j, i))
                    END DO
                END DO
            END DO
        END DO

        is_init = .TRUE.

    END SUBROUTINE sor_init


    SUBROUTINE sor(ilevel, dp, rhs, gsaw, gsae, gsas, gsan, gsab, gsat, &
            gsrap, bp)
        ! Subroutine arguments
        INTEGER(intk), INTENT(in) :: ilevel
        TYPE(field_t), INTENT(inout) :: dp
        TYPE(field_t), INTENT(in) :: rhs
        TYPE(field_t), INTENT(in) :: gsaw
        TYPE(field_t), INTENT(in) :: gsae
        TYPE(field_t), INTENT(in) :: gsas
        TYPE(field_t), INTENT(in) :: gsan
        TYPE(field_t), INTENT(in) :: gsab
        TYPE(field_t), INTENT(in) :: gsat
        TYPE(field_t), INTENT(in) :: gsrap
        TYPE(field_t), INTENT(in), OPTIONAL :: bp

        ! Local variables
        INTEGER(intk) :: i, igrid
        INTEGER(intk) :: kk, jj, ii
        REAL(realk), POINTER, CONTIGUOUS :: aw(:), ae(:), as(:), an(:), &
            ab(:), at(:), rap(:, :, :)
        REAL(realk), POINTER, CONTIGUOUS :: dp_p(:, :, :), &
            rhs_p(:, :, :), bp_p(:, :, :)

        IF (.NOT. is_init) CALL errr(__FILE__, __LINE__)
        ! Ensure this does not point to anything
        NULLIFY(bp_p)

        DO i = 1, nmygridslvl(ilevel)
            igrid = mygridslvl(i, ilevel)

            CALL dp%get_ptr(dp_p, igrid)
            CALL rhs%get_ptr(rhs_p, igrid)

            CALL gsaw%get_ptr(aw, igrid)
            CALL gsae%get_ptr(ae, igrid)
            CALL gsas%get_ptr(as, igrid)
            CALL gsan%get_ptr(an, igrid)
            CALL gsab%get_ptr(ab, igrid)
            CALL gsat%get_ptr(at, igrid)
            CALL gsrap%get_ptr(rap, igrid)

            IF (PRESENT(bp)) CALL bp%get_ptr(bp_p, igrid)

            CALL get_mgdims(kk, jj, ii, igrid)
            CALL relax(kk, jj, ii, dp_p, rhs_p, aw, ae, as, an, &
                ab, at, rap, bp_p)
        END DO
    END SUBROUTINE sor


    PURE SUBROUTINE relax(kk, jj, ii, dp, rhs, gsaw, gsae, gsas, gsan, &
            gsab, gsat, gsrap, bp)

        ! Subroutine arguments
        INTEGER(intk), INTENT(in) :: kk, jj, ii
        REAL(realk), INTENT(inout) :: dp(kk, jj, ii)
        REAL(realk), INTENT(in) :: rhs(kk, jj, ii)
        REAL(realk), INTENT(in) :: gsaw(ii), gsae(ii), gsas(jj), gsan(jj), &
            gsab(kk), gsat(kk)
        REAL(realk), INTENT(in) :: gsrap(kk, jj, ii)
        REAL(realk), INTENT(in), OPTIONAL :: bp(kk, jj, ii)

        ! Local variables
        INTEGER(intk) :: i, j, k, rb
        INTEGER(intk) :: kstart, kstop
        REAL(realk) :: res
        REAL(realk) :: aw, ae, as, an, ab, at, rap

        IF (PRESENT(bp)) THEN
            DO rb = 0, 1
                DO i = 3, ii-2
                    DO j = 3, jj-2
                        IF (MOD(i + j, 2) == rb) THEN
                            kstart = 4
                            kstop = kk - 2
                        ELSE IF (MOD(i + j, 2) /= rb) THEN
                            kstart = 3
                            kstop = kk - 3
                        END IF
#ifndef _MGLET_OFFLOAD_
                        !$omp simd private(aw, ae, as, an, ab, at, rap, res)
#endif
                        DO k = kstart, kstop, 2
                            ! Variations in numerical formulation, please
                            ! keep for future reference. Should be the same
                            ! numerics, but performance vary.
                            !
                            ! aw = bu(k, j, i-1)/(dx(i-1)*dx(i))
                            ! ae = bu(k, j, i  )/(dx(i  )*dx(i))
                            ! as = bv(k, j-1, i)/(dy(j-1)*dy(j))
                            ! an = bv(k, j  , i)/(dy(j  )*dy(j))
                            ! ab = bw(k-1, j, i)/(dz(k-1)*dz(k))
                            ! at = bw(k  , j, i)/(dz(k  )*dz(k))
                            !
                            ! aw = bp(k, j, i-1)*bp(k, j, i  )/(dx(i-1)*dx(i))
                            ! ae = bp(k, j, i  )*bp(k, j, i+1)/(dx(i  )*dx(i))
                            ! as = bp(k, j-1, i)*bp(k, j  , i)/(dy(j-1)*dy(j))
                            ! an = bp(k, j  , i)*bp(k, j+1, i)/(dy(j  )*dy(j))
                            ! ab = bp(k-1, j, i)*bp(k  , j, i)/(dz(k-1)*dz(k))
                            ! at = bp(k  , j, i)*bp(k+1, j, i)/(dz(k  )*dz(k))

                            aw = bp(k, j, i-1)*bp(k, j, i)*gsaw(i)
                            ae = bp(k, j, i)*bp(k, j, i+1)*gsae(i)
                            as = bp(k, j-1, i)*bp(k, j, i)*gsas(j)
                            an = bp(k, j, i)*bp(k, j+1, i)*gsan(j)
                            ab = bp(k-1, j, i)*bp(k, j, i)*gsab(k)
                            at = bp(k, j, i)*bp(k+1, j, i)*gsat(k)
                            rap = gsrap(k, j, i)

                            res = (aw * dp(k, j, i-1) &
                                + ae * dp(k, j, i+1) &
                                + as * dp(k, j-1, i) &
                                + an * dp(k, j+1, i) &
                                + ab * dp(k-1, j, i) &
                                + at * dp(k+1, j, i) &
                                - rhs(k, j, i)) * rap

                            dp(k, j, i) = (1.0 - omg)*dp(k, j, i) - omg*res
                        END DO
                    END DO
                END DO
            END DO
        ELSE
            DO rb = 0, 1
                DO i = 3, ii-2
                    DO j = 3, jj-2
                        IF (MOD(i + j, 2) == rb) THEN
                            kstart = 4
                            kstop = kk - 2
                        ELSE IF (MOD(i + j, 2) /= rb) THEN
                            kstart = 3
                            kstop = kk - 3
                        END IF
#ifndef _MGLET_OFFLOAD_
                        !$omp simd private(res)
#endif
                        DO k = kstart, kstop, 2
                            res = (gsaw(i) * dp(k, j, i-1) &
                                  + gsae(i) * dp(k, j, i+1) &
                                  + gsas(j) * dp(k, j-1, i) &
                                  + gsan(j) * dp(k, j+1, i) &
                                  + gsab(k) * dp(k-1, j, i) &
                                  + gsat(k) * dp(k+1, j, i) &
                                  - rhs(k, j, i)) * gsrap(k, j, i)

                            dp(k, j, i) = (1.0 - omg)*dp(k, j, i) - omg*res
                        END DO
                    END DO
                END DO
            END DO
        END IF
    END SUBROUTINE relax


    SUBROUTINE sor_finish()
        ! Subroutine arguments
        ! none...

        ! Local variables
        ! none...

        is_init = .FALSE.
    END SUBROUTINE sor_finish

END MODULE sor_mod