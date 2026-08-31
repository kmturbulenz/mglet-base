MODULE boussinesqterm_mod
    USE core_mod

    IMPLICIT NONE(type, external)
    PRIVATE

    LOGICAL, PROTECTED :: has_buoyancy = .FALSE.
    REAL(realk), PROTECTED :: gravity(3)
    REAL(realk), PROTECTED :: expansioncoefficient
    CHARACTER(len=nchar_name), PROTECTED :: scaname

    PUBLIC :: init_boussinesqterm, finish_boussinesqterm, boussinesqterm

CONTAINS

    SUBROUTINE init_boussinesqterm()
        ! Subroutine arguments
        ! none...

        ! Local variables
        ! none...

        has_buoyancy = .FALSE.
        IF (.NOT. fort7%exists("/flow/buoyancy")) THEN
            RETURN
        END IF
        has_buoyancy = .TRUE.

        ! TODO: check IB types, does not work with CC

        CALL fort7%get_array("/flow/buoyancy/gravity", gravity)
        CALL fort7%get_value("/flow/buoyancy/expansioncoefficient", &
            expansioncoefficient)
        CALL fort7%get_value("/flow/buoyancy/scalar", scaname)

        ! TODO: relationship between gravity, rho and expansion for rho /= 1.0

        IF (myid == 0) THEN
            WRITE(*, '("BUOYANCY TERM:")')
            WRITE(*, '(2X, "Gravity:              ", 3(G0, 1X))') gravity
            WRITE(*, '(2X, "Expansioncoefficient: ", G0)') &
                expansioncoefficient
            WRITE(*, '(2X, "Scalar name:          ", A)') TRIM(scaname)
            WRITE(*, '()')
        END IF
    END SUBROUTINE init_boussinesqterm


    SUBROUTINE finish_boussinesqterm
        ! Nothing to do right now...
        CONTINUE
    END SUBROUTINE finish_boussinesqterm


    SUBROUTINE boussinesqterm(uo_f, vo_f, wo_f)
        ! Subroutine arguments
        TYPE(field_t), INTENT(inout) :: uo_f, vo_f, wo_f

        ! Local variables
        TYPE(field_t), POINTER :: scafield

        IF (.NOT. has_buoyancy) RETURN
        CALL start_timer(360)

        ! Use "old" value for scalar field since scalar field already have
        ! been timestepped to next time level!
        CALL get_field(scafield, TRIM(scaname)//"_OLD")

        CALL boussinesqterm_impl(uo_f%arr, vo_f%arr, wo_f%arr, &
            scafield%arr, gravity, expansioncoefficient)

        CALL stop_timer(360)
    END SUBROUTINE boussinesqterm


    SUBROUTINE boussinesqterm_impl(uo, vo, wo, t, gravity2, expansion)
        REAL(realk), INTENT(inout) :: uo(*), vo(*), wo(*)
        REAL(realk), INTENT(in) :: t(*), gravity2(3), expansion

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
            CALL boussinesqterm_grid(kk, jj, ii, uo(ip3), vo(ip3), &
                wo(ip3), t(ip3), gravity2, expansion, nfro, nbac, nrgt, &
                nlft, nbot, ntop)
            !$omp end parallel
        END DO
        !$omp end target teams distribute
    END SUBROUTINE boussinesqterm_impl


        SUBROUTINE boussinesqterm_grid(kk, jj, ii, uo, vo, wo, t, gravity2, &
            expansion, nfro, nbac, nrgt, nlft, nbot, ntop)
        !$omp declare target

        ! Subroutine arguments
        INTEGER(intk), INTENT(in) :: kk, jj, ii
        REAL(realk), INTENT(inout) :: uo(kk, jj, ii), vo(kk, jj, ii), &
            wo(kk, jj, ii)
        REAL(realk), INTENT(in) :: t(kk, jj, ii)
        REAL(realk), INTENT(in) :: gravity2(3), expansion
        INTEGER(intk), INTENT(in) :: nfro, nbac, nrgt, nlft, nbot, ntop

        ! Local variables
        INTEGER(intk) :: k, j, i
        INTEGER(intk) :: nbu, nfu, nrv, nbw, ntw, nlv
        REAL(realk) :: term, tboundary

        nbw = 0
        ntw = 0

        nfu = 0
        nbu = 0
        nrv = 0
        nlv = 0
        nbw = 0
        ntw = 0

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

        ! U-velocioty
        IF (ABS(gravity2(1)) > 0.0) THEN
            !$omp do collapse(3) private(i, j, k, tboundary, term)
            DO i = 3-nfu, ii-3+nbu
                DO j = 3, jj-2
                    DO k = 3, kk-2
                        ! Interpolation of the scalar field to U-velocity point
                        tboundary = 0.5*(t(k, j, i) + t(k, j, i+1))

                        ! Computation of the additional force
                        term = -gravity2(1)*expansion*tboundary

                        uo(k, j, i) = uo(k, j, i) + term
                    END DO
                END DO
            END DO
            !$omp end do
        END IF

        ! V-velocioty
        IF (ABS(gravity2(2)) > 0.0) THEN
            !$omp do collapse(3) private(i, j, k, tboundary, term)
            DO i = 3, ii-2
                DO j = 3-nrv, jj-3+nlv
                    DO k = 3, kk-2
                        ! Interpolation of the scalar field to V-velocity point
                        tboundary = 0.5*(t(k, j, i) + t(k, j+1, i))

                        ! Computation of the additional force
                        term = -gravity2(2)*expansion*tboundary

                        vo(k, j, i) = vo(k, j, i) + term
                    END DO
                END DO
            END DO
            !$omp end do
        END IF


        ! W-velocioty
        IF (ABS(gravity2(3)) > 0.0) THEN
            !$omp do collapse(3) private(i, j, k, tboundary, term)
            DO i = 3, ii-2
                DO j = 3, jj-2
                    DO k = 3-nbw, kk-3+ntw
                        ! Interpolation of the scalar field to W-velocity point
                        tboundary = 0.5*(t(k, j, i) + t(k+1, j, i))

                        ! Computation of the additional force
                        term = -gravity2(3)*expansion*tboundary

                        wo(k, j, i) = wo(k, j, i) + term
                    END DO
                END DO
            END DO
            !$omp end do
        END IF
    END SUBROUTINE boussinesqterm_grid
END MODULE boussinesqterm_mod
