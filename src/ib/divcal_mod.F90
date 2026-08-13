MODULE divcal_mod
    USE core_mod

    IMPLICIT NONE(type, external)
    PRIVATE

    PUBLIC :: divcal
CONTAINS
    SUBROUTINE divcal(div_f, u_f, v_f, w_f, fak, device)
        ! Subroutine arguments
        TYPE(field_t), INTENT(inout) :: div_f
        TYPE(field_t), INTENT(in) :: u_f, v_f, w_f
        REAL(realk), INTENT(in) :: fak
        LOGICAL, OPTIONAL, INTENT(in) :: device

        ! Local variables
        TYPE(field_t), POINTER :: rddx_f, rddy_f, rddz_f, bp_f
        LOGICAL, ALLOCATABLE :: active_levels(:)
        INTEGER(intk) :: i, igrid, ilevel, ip3, ip1x, ip1y, ip1z, kk, jj, ii
        LOGICAL :: device2

        CALL start_timer(240)

        ALLOCATE(active_levels(1:maxlevel-minlevel+1))

        IF (PRESENT(device)) THEN
            device2 = device
        ELSE
            device2 = .FALSE.
        END IF

        CALL get_field(rddx_f, "RDDX")
        CALL get_field(rddy_f, "RDDY")
        CALL get_field(rddz_f, "RDDZ")
        CALL get_field(bp_f, "BP")

        active_levels = u_f%active_level(minlevel:maxlevel)

        ASSOCIATE(div => div_f%arr(:), &
            u => u_f%arr(:), v => v_f%arr(:), w => w_f%arr(:), &
            rddx => rddx_f%arr(:), &
            rddy => rddy_f%arr(:), &
            rddz => rddz_f%arr(:), &
            bp => bp_f%arr(:))

        ! active_level must be handled specially because it might start at
        ! minlvl=0 or minlvl=1. If we want to avoid copying descriptors, then
        ! it should be ensured that active_level always starts at minlvl=1.
        ! This explains the allocatable array active_levels, which is a copy of
        ! the active_level array of u_f.

        !$omp target teams distribute private(kk, jj, ii, ip3, ip1x, ip1y, &
        !$omp& ip1z, igrid, ilevel) &
        !$omp& map(to: active_levels) &
        !$omp& if(device2)
        DO i = 1, nmygrids
            igrid = mygrids(i)
            ilevel = level(igrid)

            ! Assume that U, V, W and DIV are defined on the same levels!!!
            IF (.NOT. active_levels(ilevel-minlevel+1)) CYCLE

            CALL get_mgdims(kk, jj, ii, igrid)
            CALL get_ip3(ip3, igrid)
            CALL get_ip1x(ip1x, igrid)
            CALL get_ip1y(ip1y, igrid)
            CALL get_ip1z(ip1z, igrid)

            !$omp parallel
            CALL divcal_grid(kk, jj, ii, div(ip3), u(ip3), v(ip3), w(ip3), &
                bp(ip3), rddx(ip1x), rddy(ip1y), rddz(ip1z), fak)
            !$omp end parallel
        END DO
        !$omp end target teams distribute

        END ASSOCIATE

        CALL stop_timer(240)
    END SUBROUTINE divcal


    SUBROUTINE divcal_grid(kk, jj, ii, div, u, v, w, bp, rddx, rddy, rddz, fak)
        !$omp declare target
        ! Subroutine arguments
        INTEGER(intk), INTENT(in) :: kk, jj, ii
        REAL(realk), INTENT(inout) :: div(kk, jj, ii)
        REAL(realk), INTENT(in) :: u(kk, jj, ii)
        REAL(realk), INTENT(in) :: v(kk, jj, ii)
        REAL(realk), INTENT(in) :: w(kk, jj, ii)
        REAL(realk), INTENT(in) :: bp(kk, jj, ii)
        REAL(realk), INTENT(in) :: rddx(ii), rddy(jj), rddz(kk)
        REAL(realk), INTENT(in) :: fak

        ! Local variables
        INTEGER(intk) :: k, j, i

        !$omp do collapse(3) private(k, j, i)
        DO i = 3, ii-2
            DO j = 3, jj-2
                DO k = 3, kk-2
                    div(k, j, i) = fak*((u(k, j, i) - u(k, j, i-1))*rddx(i) &
                        + (v(k, j, i) - v(k, j-1, i))*rddy(j) &
                        + (w(k, j, i) - w(k-1, j, i))*rddz(k)) &
                        * bp(k, j, i)
                END DO
            END DO
        END DO
        !$omp end do
    END SUBROUTINE divcal_grid
END MODULE divcal_mod
