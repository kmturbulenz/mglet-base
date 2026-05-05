
MODULE laplacephi_mod

    USE core_mod

    IMPLICIT NONE(type, external)

    PRIVATE

    PUBLIC :: laplacephi, laplacephi_level

CONTAINS

    SUBROUTINE laplacephi(res, phi, bp)
        ! Subroutine arguments
        TYPE(field_t), INTENT(inout) :: res
        TYPE(field_t), INTENT(in) :: phi
        TYPE(field_t), INTENT(in), OPTIONAL :: bp

        ! Local variables
        INTEGER(intk) :: ilevel

        DO ilevel = minlevel, maxlevel
            CALL laplacephi_level(ilevel, res, phi, bp)
        END DO
    END SUBROUTINE laplacephi


    SUBROUTINE laplacephi_level(ilevel, res_f, phi_f, bp_f)
        ! Subroutine arguments
        INTEGER(intk), INTENT(in):: ilevel
        TYPE(field_t), INTENT(inout) :: res_f
        TYPE(field_t), INTENT(in) :: phi_f
        TYPE(field_t), INTENT(in), OPTIONAL :: bp_f

        ! Local variables
        INTEGER(intk) :: i, igrid
        INTEGER(intk) :: kk, jj, ii

        TYPE(field_t), POINTER :: gsaw, gsae, gsas, gsan, gsab, gsat, gsap
        REAL(realk), POINTER, CONTIGUOUS :: phi(:, :, :), res(:, :, :)
        REAL(realk), POINTER, CONTIGUOUS :: aw(:), ae(:), as(:), an(:), &
            ab(:), at(:), ap(:, :, :)
        REAL(realk), POINTER, CONTIGUOUS :: bp(:, :, :)

        NULLIFY(bp)

        CALL get_field(gsaw, "GSAW")
        CALL get_field(gsae, "GSAE")
        CALL get_field(gsas, "GSAS")
        CALL get_field(gsan, "GSAN")
        CALL get_field(gsab, "GSAB")
        CALL get_field(gsat, "GSAT")
        CALL get_field(gsap, "GSAP")

        DO i = 1, nmygridslvl(ilevel)
            igrid = mygridslvl(i, ilevel)
            CALL get_mgdims(kk, jj, ii, igrid)

            CALL phi_f%get_ptr(phi, igrid)
            CALL res_f%get_ptr(res, igrid)

            CALL gsaw%get_ptr(aw, igrid)
            CALL gsae%get_ptr(ae, igrid)
            CALL gsas%get_ptr(as, igrid)
            CALL gsan%get_ptr(an, igrid)
            CALL gsab%get_ptr(ab, igrid)
            CALL gsat%get_ptr(at, igrid)
            CALL gsap%get_ptr(ap, igrid)

            IF (PRESENT(bp_f)) THEN
                CALL bp_f%get_ptr(bp, igrid)
            END IF

            CALL laplacephi_grid(kk, jj, ii, res, phi, aw, ae, an, as, &
                at, ab, ap, bp)
        END DO
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

        IF (PRESENT(bp)) THEN
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
        ELSE
            DO i = 3, ii-2
                DO j = 3, jj-2
                    DO k = 3, kk-2
                        res(k, j, i) = &
                            - aw(i) * phi(k, j, i-1) &
                            - ae(i) * phi(k, j, i+1) &
                            - as(j) * phi(k, j-1, i) &
                            - an(j) * phi(k, j+1, i) &
                            - ab(k) * phi(k-1, j, i) &
                            - at(k) * phi(k+1, j, i) &
                            - ap(k, j, i) * phi(k, j, i)
                    END DO
                END DO
            END DO
        END IF

    END SUBROUTINE laplacephi_grid

END MODULE laplacephi_mod