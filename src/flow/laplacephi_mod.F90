
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
        INTEGER(intk) :: i, igrid, kk, jj, ii, ip3, ipx, ipy, ipz
        TYPE(field_t), POINTER :: gsaw, gsae, gsas, gsan, gsab, gsat, gsap, bp_f
        REAL(realk), ALLOCATABLE, DIMENSION(:) :: phi, res, ap, bp, &
            aw, ae, an, as, at, ab

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

        ASSOCIATE ( &
            phi => phi_f%arr, &
            res => res_f%arr, &
            ap  => gsap%arr, &
            aw  => gsaw%arr, &
            ae  => gsae%arr, &
            an  => gsan%arr, &
            as  => gsas%arr, &
            at  => gsat%arr, &
            ab  => gsab%arr, &
            bp  => bp_f%arr)

        !$omp target teams distribute private(i, igrid, kk, jj, ii, ip3, ipx, &
        !$omp& ipy, ipz)
        DO i = 1, nmygrids
            igrid = mygrids(i)
            CALL get_mgdims(kk, jj, ii, igrid)

            CALL get_ip3(ip3, igrid)
            CALL get_ipx(ipx, igrid)
            CALL get_ipy(ipy, igrid)
            CALL get_ipz(ipz, igrid)

            CALL laplacephi_grid(kk, jj, ii, res(ip3), phi(ip3), &
                aw(ipx), ae(ipx), an(ipy), as(ipy), at(ipz), ab(ipz), &
                ap(ip3), bp(ip3))
        END DO
        !$omp end target teams distribute

        END ASSOCIATE
    END SUBROUTINE laplacephi


    SUBROUTINE laplacephi_level(ilevel, res_f, phi_f)
        ! Subroutine arguments
        INTEGER(intk), INTENT(in):: ilevel
        TYPE(field_t), TARGET, INTENT(inout) :: res_f
        TYPE(field_t), TARGET, INTENT(in) :: phi_f

        ! Local variables
        INTEGER(intk) :: i, igrid
        INTEGER(intk) :: kk, jj, ii, ip3, ipx, ipy, ipz

        TYPE(field_t), POINTER :: gsaw, gsae, gsas, gsan, gsab, gsat, gsap, bp_f
        REAL(realk), ALLOCATABLE, DIMENSION(:) :: phi, res, ap, bp, &
            aw, ae, an, as, at, ab

        CALL get_field(gsaw, "GSAW")
        CALL get_field(gsae, "GSAE")
        CALL get_field(gsas, "GSAS")
        CALL get_field(gsan, "GSAN")
        CALL get_field(gsab, "GSAB")
        CALL get_field(gsat, "GSAT")
        CALL get_field(gsap, "GSAP")
        CALL get_field(bp_f, "BP")

        ASSOCIATE ( &
            phi => phi_f%arr, &
            res => res_f%arr, &
            ap  => gsap%arr, &
            aw  => gsaw%arr, &
            ae  => gsae%arr, &
            an  => gsan%arr, &
            as  => gsas%arr, &
            at  => gsat%arr, &
            ab  => gsab%arr, &
            bp  => bp_f%arr)

        ! use phi, res, ap, aw, ae, an, as, at, ab here
        ! BP for noib is 1.0 anyways. Take the few multiplications instead of
        ! branching in the kernel.

        CALL profile_range_push("laplacephi_level")
        !$omp target teams distribute private(i, igrid, kk, jj, ii, ip3, ipx, &
        !$omp& ipy, ipz)
        DO i = 1, nmygridslvl(ilevel)
            igrid = mygridslvl(i, ilevel)
            CALL get_mgdims(kk, jj, ii, igrid)

            CALL get_ip3(ip3, igrid)
            CALL get_ipx(ipx, igrid)
            CALL get_ipy(ipy, igrid)
            CALL get_ipz(ipz, igrid)

            CALL laplacephi_grid(kk, jj, ii, res(ip3), phi(ip3), &
                aw(ipx), ae(ipx), an(ipy), as(ipy), at(ipz), ab(ipz), &
                ap(ip3), bp(ip3))
        END DO
        !$omp end target teams distribute
        CALL profile_range_pop()
        END ASSOCIATE
    END SUBROUTINE laplacephi_level


    PURE SUBROUTINE laplacephi_grid(kk, jj, ii, res, phi, aw, ae, an, as, &
            at, ab, ap, bp)
        !$omp declare target
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
