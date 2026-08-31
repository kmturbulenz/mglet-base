
MODULE coriolisterm_mod
    USE core_mod

    IMPLICIT NONE(type, external)
    PRIVATE

    ! activity flag and rotation vector "omega"
    LOGICAL, PROTECTED :: has_coriolis = .FALSE.
    REAL(realk), PROTECTED :: omega(3)

    ! module functions
    PUBLIC :: init_coriolisterm, finish_coriolisterm, coriolisterm

CONTAINS

    SUBROUTINE init_coriolisterm()

        ! leaving inactive if no parameters specified
        has_coriolis = .FALSE.
        IF (.NOT. fort7%exists("/flow/coriolis")) RETURN

        ! retrieving rotation rate vector from parameters.json
        CALL fort7%get_array("/flow/coriolis/omega", omega)

        ! display obtained parameters
        IF (myid == 0) THEN
            WRITE(*, '("CORIOLIS TERM:")')
            WRITE(*, '(2X, "omega (rotation vector):   ", 3(G0, 1X))') omega
            WRITE(*, '()')
        END IF

        ! set active
        has_coriolis = .TRUE.

    END SUBROUTINE init_coriolisterm


    SUBROUTINE finish_coriolisterm

        ! revoking activity
        has_coriolis = .FALSE.

        RETURN

    END SUBROUTINE finish_coriolisterm


    SUBROUTINE coriolisterm(uo_f, vo_f, wo_f)

        ! subroutine arguments (fields storing momentum balance)
        TYPE(field_t), INTENT(inout) :: uo_f, vo_f, wo_f

        ! local variables
        TYPE(field_t), POINTER :: u_f, v_f, w_f

        ! checking activity
        IF (.NOT. has_coriolis) RETURN

        CALL start_timer(370)

        ! getting the velocity fields
        CALL get_field(u_f, "U")
        CALL get_field(v_f, "V")
        CALL get_field(w_f, "W")

        CALL coriolisterm_impl(uo_f%arr, vo_f%arr, wo_f%arr, u_f%arr, &
            v_f%arr, w_f%arr, omega(1), omega(2), omega(3))

        CALL stop_timer(370)
    END SUBROUTINE coriolisterm


    SUBROUTINE coriolisterm_impl(uo, vo, wo, u, v, w, omega1, omega2, omega3)
        REAL(realk), INTENT(inout) :: uo(*), vo(*), wo(*)
        REAL(realk), INTENT(in) :: u(*), v(*), w(*)
        REAL(realk), INTENT(in) :: omega1, omega2, omega3

        INTEGER(intk) :: i, igrid, kk, jj, ii, ip3
        INTEGER(intk) :: nfro, nbac, nrgt, nlft, nbot, ntop

        !$omp target teams distribute private(i, igrid, kk, jj, ii, ip3, &
        !$omp& nfro, nbac, nrgt, nlft, nbot, ntop)
        DO i = 1, nmygrids
            igrid = mygrids(i)
            CALL get_mgdims(kk, jj, ii, igrid)
            CALL get_mgbasb(nfro, nbac, nrgt, nlft, nbot, ntop, igrid)
            CALL get_ip3(ip3, igrid)
            !$omp parallel
            CALL coriolis_grid(kk, jj, ii, uo(ip3), vo(ip3), wo(ip3), &
                u(ip3), v(ip3), w(ip3), omega1, omega2, omega3, nfro, &
                nbac, nrgt, nlft, nbot, ntop)
            !$omp end parallel
        END DO
        !$omp end target teams distribute
    END SUBROUTINE coriolisterm_impl


    SUBROUTINE coriolis_grid(kk, jj, ii, uo, vo, wo, u, v, w, omega1, omega2, &
            omega3, nfro, nbac, nrgt, nlft, nbot, ntop)
        !$omp declare target

        ! subroutine arguments
        INTEGER(intk), INTENT(in) :: kk, jj, ii
        REAL(realk), INTENT(inout) :: uo(kk, jj, ii), &
            vo(kk, jj, ii), wo(kk, jj, ii)
        REAL(realk), INTENT(in) :: u(kk, jj, ii), &
            v(kk, jj, ii), w(kk, jj, ii)
        REAL(realk), INTENT(in) :: omega1, omega2, omega3
        INTEGER(intk), INTENT(in) :: nfro, nbac, nrgt, nlft, nbot, ntop

        ! local variables
        INTEGER(intk) :: k, j, i
        INTEGER(intk) :: nbu, nfu, nrv, nbw, ntw, nlv
        REAL(realk) :: vlocal1, vlocal2, vlocal3, cterm

        nfu = 0; nbu = 0; nrv = 0
        nlv = 0; nbw = 0; ntw = 0

        ! CON = 7
        IF (nbac == 7) nbu = 1
        IF (nlft == 7) nlv = 1
        IF (ntop == 7) ntw = 1

        ! OP1 = 3
        IF (nfro == 3) nfu = 1
        IF (nbac == 3) nbu = 1
        IF (nrgt == 3) nrv = 1
        IF (nlft == 3) nlv = 1
        IF (nbot == 3) nbw = 1
        IF (ntop == 3) ntw = 1

        ! comment: Not yet implemented with a linear
        ! interploation to the exact position

        IF (ABS(omega2) + ABS(omega3) > 0.0) THEN
            !$omp do collapse(3) private(i, j, k, vlocal2, vlocal3, cterm)
            DO i = 3-nfu, ii-3+nbu
                DO j = 3, jj-2
                    DO k = 3, kk-2
                        ! averaging to U-velocity point (stag=1,0,0)
                        vlocal2 = 0.25 * (v(k, j-1, i) + &
                            v(k, j-1, i+1) + v(k, j, i) + v(k, j, i+1))
                        vlocal3 = 0.25 * (w(k-1, j, i) + &
                            w(k-1, j, i+1) + w(k, j, i) + w(k, j, i+1))
                        ! computing the cross product
                        cterm = -2.0 * (omega2*vlocal3 - omega3*vlocal2)
                        ! adding to the momentum balance
                        uo(k, j, i) = uo(k, j, i) + cterm
                    END DO
                END DO
            END DO
            !$omp end do
        END IF

        IF (ABS(omega1) + ABS(omega3) > 0.0) THEN
            !$omp do collapse(3) private(i, j, k, vlocal1, vlocal3, cterm)
            DO i = 3, ii-2
                DO j = 3-nrv, jj-3+nlv
                    DO k = 3, kk-2
                        ! averaging to V-velocity point (stag=0,1,0)
                        vlocal1 = 0.25 * (u(k, j, i-1) + &
                            u(k, j+1, i-1) + u(k, j, i) + u(k, j+1, i))
                        vlocal3 = 0.25 * (w(k-1, j, i) + &
                            w(k-1, j+1, i) + w(k, j, i) + w(k, j+1, i))
                        ! computing the cross product
                        cterm = -2.0 * (omega3*vlocal1 - omega1*vlocal3)
                        ! adding to the momentum balance
                        vo(k, j, i) = vo(k, j, i) + cterm
                    END DO
                END DO
            END DO
            !$omp end do
        END IF

        IF (ABS(omega1) + ABS(omega2) > 0.0) THEN
            !$omp do collapse(3) private(i, j, k, vlocal1, vlocal2, cterm)
            DO i = 3, ii-2
                DO j = 3, jj-2
                    DO k = 3-nbw, kk-3+ntw
                        ! averaging to W-velocity point (stag=0,0,1)
                        vlocal1 = 0.25 * (u(k, j, i-1) + &
                            u(k+1, j, i-1) + u(k, j, i) + u(k+1, j, i))
                        vlocal2 = 0.25 * (v(k, j-1, i) + &
                            v(k+1, j-1, i) + v(k, j, i) + v(k+1, j, i))
                        ! computing the cross product
                        cterm = -2.0 * (omega1*vlocal2 - omega2*vlocal1)
                        ! adding to the momentum balance
                        wo(k, j, i) = wo(k, j, i) + cterm
                    END DO
                END DO
            END DO
            !$omp end do
        END IF

    END SUBROUTINE coriolis_grid

END MODULE coriolisterm_mod
