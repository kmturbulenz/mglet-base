
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
        INTEGER(intk) :: i, igrid, kk, jj, ii
        INTEGER(intk) :: nfro, nbac, nrgt, nlft, nbot, ntop
        TYPE(field_t), POINTER :: u_f, v_f, w_f
        REAL(realk), POINTER, CONTIGUOUS, DIMENSION(:, :, :) :: uo, vo, wo
        REAL(realk), POINTER, CONTIGUOUS, DIMENSION(:, :, :) :: u, v, w

        ! checking activity
        IF (.NOT. has_coriolis) RETURN

        CALL start_timer(370)

        ! getting the velocity fields
        CALL get_field(u_f, "U")
        CALL get_field(v_f, "V")
        CALL get_field(w_f, "W")

        ! iterating over all grids on processor (from grids_mod.F90)
        DO i = 1, nmygrids
            igrid = mygrids(i)

            ! getting grid parameters (dimensions and boundary conditions)
            CALL get_mgdims(kk, jj, ii, igrid)
            CALL get_mgbasb(nfro, nbac, nrgt, nlft, nbot, ntop, igrid)

            CALL uo_f%get_ptr(uo, igrid)
            CALL vo_f%get_ptr(vo, igrid)
            CALL wo_f%get_ptr(wo, igrid)

            CALL u_f%get_ptr(u, igrid)
            CALL v_f%get_ptr(v, igrid)
            CALL w_f%get_ptr(w, igrid)

            ! calling the function kernel for each grid
            CALL coriolis_grid(kk, jj, ii, &
                uo, vo, wo, u, v, w, &
                nfro, nbac, nrgt, nlft, nbot, ntop)
        END DO

        CALL stop_timer(370)

        RETURN

    END SUBROUTINE coriolisterm


    SUBROUTINE coriolis_grid(kk, jj, ii, &
        uo, vo, wo, u, v, w, &
        nfro, nbac, nrgt, nlft, nbot, ntop)

        ! subroutine arguments
        INTEGER(intk), INTENT(in) :: kk, jj, ii
        REAL(realk), INTENT(inout) :: uo(kk, jj, ii), &
            vo(kk, jj, ii), wo(kk, jj, ii)
        REAL(realk), INTENT(in) :: u(kk, jj, ii), &
            v(kk, jj, ii), w(kk, jj, ii)
        INTEGER(intk), INTENT(in) :: nfro, nbac, nrgt, nlft, nbot, ntop

        ! local variables
        INTEGER(intk) :: k, j, i
        INTEGER(intk) :: nbu, nfu, nrv, nbw, ntw, nlv
        REAL(realk) :: vlocal(3), cterm

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

        IF (ABS(omega(2)) + ABS(omega(3)) > 0.0) THEN
            DO i = 3-nfu, ii-3+nbu
                DO j = 3, jj-2
                    DO k = 3, kk-2
                        ! averaging to U-velocity point (stag=1,0,0)
                        vlocal(1) = u(k, j, i)
                        vlocal(2) = 0.25 * (v(k, j-1, i) + &
                            v(k, j-1, i+1) + v(k, j, i) + v(k, j, i+1))
                        vlocal(3) = 0.25 * (w(k-1, j, i) + &
                            w(k-1, j, i+1) + w(k, j, i) + w(k, j, i+1))
                        ! computing the cross product
                        cterm = -2.0 * (omega(2) * vlocal(3) - &
                            omega(3) * vlocal(2))
                        ! adding to the momentum balance
                        uo(k, j, i) = uo(k, j, i) + cterm
                    END DO
                END DO
            END DO
        END IF

        IF (ABS(omega(1)) + ABS(omega(3)) > 0.0) THEN
            DO i = 3, ii-2
                DO j = 3-nrv, jj-3+nlv
                    DO k = 3, kk-2
                        ! averaging to V-velocity point (stag=0,1,0)
                        vlocal(1) = 0.25 * (u(k, j, i-1) + &
                            u(k, j+1, i-1) + u(k, j, i) + u(k, j+1, i))
                        vlocal(2) = v(k, j, i)
                        vlocal(3) = 0.25 * (w(k-1, j, i) + &
                            w(k-1, j+1, i) + w(k, j, i) + w(k, j+1, i))
                        ! computing the cross product
                        cterm = -2.0 * (omega(3) * vlocal(1) - &
                            omega(1) * vlocal(3))
                        ! adding to the momentum balance
                        vo(k, j, i) = vo(k, j, i) + cterm
                    END DO
                END DO
            END DO
        END IF

        IF (ABS(omega(1)) + ABS(omega(2)) > 0.0) THEN
            DO i = 3, ii-2
                DO j = 3, jj-2
                    DO k = 3-nbw, kk-3+ntw
                        ! averaging to W-velocity point (stag=0,0,1)
                        vlocal(1) = 0.25 * (u(k, j, i-1) + &
                            u(k+1, j, i-1) + u(k, j, i) + u(k+1, j, i))
                        vlocal(2) = 0.25 * (v(k, j-1, i) + &
                            v(k+1, j-1, i) + v(k, j, i) + v(k+1, j, i))
                        vlocal(3) = w(k, j, i)
                        ! computing the cross product
                        cterm = -2.0 * (omega(1) * vlocal(2) - &
                            omega(2) * vlocal(1))
                        ! adding to the momentum balance
                        wo(k, j, i) = wo(k, j, i) + cterm
                    END DO
                END DO
            END DO
        END IF

    END SUBROUTINE coriolis_grid

END MODULE coriolisterm_mod
