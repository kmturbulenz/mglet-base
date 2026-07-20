
MODULE laplacephi_mod

    USE core_mod

    IMPLICIT NONE(type, external)

    PRIVATE

    PUBLIC :: laplacephi, laplacephi_level

CONTAINS

    SUBROUTINE laplacephi(res_f, phi_f)
        ! Subroutine arguments
        TYPE(field_t), INTENT(inout) :: res_f
        TYPE(field_t), INTENT(in) :: phi_f

        ! Local variables
        INTEGER(intk) :: i, igrid, kk, jj, ii
        TYPE(field_t), POINTER :: gsaw, gsae, gsas, gsan, gsab, gsat, gsap, bp_f
        REAL(realk), POINTER, CONTIGUOUS :: phi(:, :, :), res(:, :, :)
        REAL(realk), POINTER, CONTIGUOUS :: aw(:), ae(:), as(:), an(:), &
            ab(:), at(:), ap(:, :, :)
        REAL(realk), POINTER, CONTIGUOUS :: bp(:, :, :)

        CALL get_field(gsaw, "GSAW")
        CALL get_field(gsae, "GSAE")
        CALL get_field(gsas, "GSAS")
        CALL get_field(gsan, "GSAN")
        CALL get_field(gsab, "GSAB")
        CALL get_field(gsat, "GSAT")
        CALL get_field(gsap, "GSAP")

        ! BP for noib is 1.0 anyways. Take the few multiplications instead of
        ! branching in the kernel.
        CALL get_field(bp_f, "BP")

        !$omp target teams distribute private(i, igrid, kk, jj, ii, phi, res, &
        !$omp& aw, ae, an, as, at, ab, ap, bp)
        DO i = 1, nmygrids
            igrid = mygrids(i)
            CALL get_mgdims(kk, jj, ii, igrid)

            CALL get_grid3_real(phi, phi_f, igrid)
            CALL get_grid3_real(res, res_f, igrid)
            CALL get_grid3_real(ap, gsap, igrid)
            CALL get_grid3_real(bp, bp_f, igrid)

            CALL get_grid1_real(aw, gsaw, igrid)
            CALL get_grid1_real(ae, gsae, igrid)
            CALL get_grid1_real(as, gsas, igrid)
            CALL get_grid1_real(an, gsan, igrid)
            CALL get_grid1_real(ab, gsab, igrid)
            CALL get_grid1_real(at, gsat, igrid)

            CALL laplacephi_grid(kk, jj, ii, res, phi, aw, ae, an, as, &
                at, ab, ap, bp)
        END DO
        !$omp end target teams distribute
    END SUBROUTINE laplacephi


    SUBROUTINE laplacephi_level(ilevel, res_f, phi_f)
        ! Subroutine arguments
        INTEGER(intk), INTENT(in):: ilevel
        TYPE(field_t), INTENT(inout) :: res_f
        TYPE(field_t), INTENT(in) :: phi_f

        ! Local variables
        INTEGER(intk) :: i, igrid
        INTEGER(intk) :: kk, jj, ii

        TYPE(field_t), POINTER :: gsaw, gsae, gsas, gsan, gsab, gsat, gsap, bp_f
        REAL(realk), POINTER, CONTIGUOUS :: phi(:, :, :), res(:, :, :)
        REAL(realk), POINTER, CONTIGUOUS :: aw(:), ae(:), as(:), an(:), &
            ab(:), at(:), ap(:, :, :)
        REAL(realk), POINTER, CONTIGUOUS :: bp(:, :, :)

        CALL get_field(gsaw, "GSAW")
        CALL get_field(gsae, "GSAE")
        CALL get_field(gsas, "GSAS")
        CALL get_field(gsan, "GSAN")
        CALL get_field(gsab, "GSAB")
        CALL get_field(gsat, "GSAT")
        CALL get_field(gsap, "GSAP")

        ! BP for noib is 1.0 anyways. Take the few multiplications instead of
        ! branching in the kernel.
        CALL get_field(bp_f, "BP")

        !$omp target teams distribute private(i, igrid, kk, jj, ii, phi, res, &
        !$omp& aw, ae, an, as, at, ab, ap, bp)
        DO i = 1, nmygridslvl(ilevel)
            igrid = mygridslvl(i, ilevel)
            CALL get_mgdims(kk, jj, ii, igrid)

            CALL get_grid3_real(phi, phi_f, igrid)
            CALL get_grid3_real(res, res_f, igrid)
            CALL get_grid3_real(ap, gsap, igrid)
            CALL get_grid3_real(bp, bp_f, igrid)

            CALL get_grid1_real(aw, gsaw, igrid)
            CALL get_grid1_real(ae, gsae, igrid)
            CALL get_grid1_real(as, gsas, igrid)
            CALL get_grid1_real(an, gsan, igrid)
            CALL get_grid1_real(ab, gsab, igrid)
            CALL get_grid1_real(at, gsat, igrid)

            CALL laplacephi_grid(kk, jj, ii, res, phi, aw, ae, an, as, &
                at, ab, ap, bp)
        END DO
        !$omp end target teams distribute
    END SUBROUTINE laplacephi_level


    PURE SUBROUTINE laplacephi_grid(kk, jj, ii, res, phi, aw, ae, an, as, &
            at, ab, ap, bp)
        ! Subroutine arguments
        INTEGER(intk), INTENT(in) :: kk, jj, ii
        REAL(realk), INTENT(inout) :: res(kk, jj, ii)
        REAL(realk), INTENT(in) :: phi(kk, jj, ii)
        REAL(realk), INTENT(in) :: aw(ii), ae(ii), an(jj), as(jj), &
            at(kk), ab(kk)
        REAL(realk), INTENT(in) :: ap(kk, jj, ii)
        REAL(realk), INTENT(in), OPTIONAL :: bp(kk, jj, ii)

        ! Local variables
        INTEGER :: k, j, i

        !$omp parallel do collapse(3) private(i, j, k)
        DO i = 3, ii-2
            DO j = 3, jj-2
                DO k = 3, kk-2
                    res(k, j, i) = &
                        - aw(i)*phi(k, j, i-1)*bp(k, j, i-1)*bp(k, j, i) &
                        - ae(i)*phi(k, j, i+1)*bp(k, j, i)*bp(k, j, i+1) &
                        - as(j)*phi(k, j-1, i)*bp(k, j-1, i)*bp(k, j, i) &
                        - an(j)*phi(k, j+1, i)*bp(k, j, i)*bp(k, j+1, i) &
                        - ab(k)*phi(k-1, j, i)*bp(k-1, j, i)*bp(k, j, i) &
                        - at(k)*phi(k+1, j, i)*bp(k, j, i)*bp(k+1, j, i) &
                        - ap(k, j, i)*phi(k, j, i)
                END DO
            END DO
        END DO
        !$omp end parallel do
    END SUBROUTINE laplacephi_grid
END MODULE laplacephi_mod
