MODULE rungekutta_mod
    USE charfunc_mod, ONLY: lower
    USE err_mod, ONLY: errr
    USE precision_mod, ONLY: intk, realk

    IMPLICIT NONE(type, external)
    PRIVATE

    INTEGER(intk), PARAMETER :: maxnrk = 6
    REAL(realk), PARAMETER :: oneeps = (1.0_realk - EPSILON(1.0_realk))

    TYPE, ABSTRACT :: rk_t
        INTEGER(intk) :: nrk = 0
        REAL(realk) :: cflmax = 0.0
        REAL(realk) :: c(maxnrk) = 0.0
    CONTAINS
        PROCEDURE(init_i), DEFERRED :: init
    END TYPE rk_t

    ABSTRACT INTERFACE
        SUBROUTINE init_i(this, ctyp)
            IMPORT :: rk_t
            CLASS(rk_t), INTENT(out) :: this
            CHARACTER(len=*), INTENT(in) :: ctyp
        END SUBROUTINE init_i
    END INTERFACE

    ! 2N storage versions
    ! A more logical name would be 2n_rk_t but it is forbidden to have a
    ! variable name that starts with a digit
    TYPE, EXTENDS(rk_t) :: rk_2n_t
        REAL(realk) :: a(maxnrk) = 0.0
        REAL(realk) :: b(maxnrk) = 0.0
    CONTAINS
        PROCEDURE :: init => init_2n
        GENERIC :: get_coeffs => get_coeffs_a, get_coeffs_b

        PROCEDURE, PRIVATE :: get_coeffs_a, get_coeffs_b
        PROCEDURE, PRIVATE :: init_williamson
        PROCEDURE, PRIVATE :: init_berland
        PROCEDURE, PRIVATE :: init_carpenter
        PROCEDURE, PRIVATE :: init_bernardini
        PROCEDURE, PRIVATE :: comp_c
    END TYPE rk_2n_t

    PUBLIC :: rk_2n_t, rkstep
CONTAINS

    SUBROUTINE init_2n(this, ctyp)
        ! Function arguments
        CLASS(rk_2n_t), INTENT(out) :: this
        CHARACTER(len=*), INTENT(in) :: ctyp

        SELECT CASE(lower(TRIM(ctyp)))
        CASE ("williamson")
            CALL this%init_williamson()
        CASE ("berland")
            CALL this%init_berland()
        CASE ("carpenter")
            CALL this%init_carpenter()
        CASE ("bernardini")
            CALL this%init_bernardini()
        CASE DEFAULT
            WRITE(*, *) "Invalid: ", lower(TRIM(ctyp))
            CALL errr(__FILE__, __LINE__)
        END SELECT
    END SUBROUTINE init_2n


    SUBROUTINE init_williamson(rk)
        ! Classic Williamson RK3 scheme
        ! From: Williamson, John H, Low-storage Runge-Kutta schemes. Journal
        ! of computational physics 35.1 (1980)
        ! DOI: https://doi.org/10.1016/0021-9991(80)90033-9

        ! Subroutine arguments
        CLASS(rk_2n_t), INTENT(out) :: rk

        ! Set to zero
        rk%c = 0.0
        rk%a = 0.0
        rk%b = 0.0

        ! Used coefficients
        rk%nrk = 3
        rk%cflmax = oneeps*SQRT(3.0)
        rk%c(1:rk%nrk) = [0.0, 1.0/3.0, 3.0/4.0]
        rk%a(1:rk%nrk) = [0.0, -5.0/9.0, -153.0/128.0]
        rk%b(1:rk%nrk) = [1.0/3.0, 15.0/16.0, 8.0/15.0]

        ! Compute C from A and B
        ! CALL rk%comp_c()
    END SUBROUTINE init_williamson


    SUBROUTINE init_berland(rk)
        ! 4th order 6-stage low disipation scheme
        ! From: Berland, Bogey and Bailly, Low-dissipation and low-dispersion
        ! fourth-order Runge–Kutta algorithm, Computers & Fluids 35 (2006),
        ! DOI: 10.1016/j.compfluid.2005.04.003

        ! Subroutine arguments
        CLASS(rk_2n_t), INTENT(out) :: rk

        ! Set to zero
        rk%a = 0.0
        rk%b = 0.0

        ! Used coefficients
        rk%nrk = 6
        rk%cflmax = 3.80 - 0.01  ! 3.80 from paper, 0.01 round off margin
        rk%a(1:rk%nrk) = [0.0, -0.737101392796, -1.634740794341, &
            -0.744739003780, -1.469897351522, -2.813971388035]
        rk%b(1:rk%nrk) = [0.032918605146, 0.823256998200, 0.381530948900, &
            0.200092213184, 1.718581042715, 0.27]

        ! Compute C from A and B
        CALL rk%comp_c()
    END SUBROUTINE init_berland


    SUBROUTINE init_carpenter(rk)
        ! 4th order 5-stage scheme
        ! Ref: Carpenter & Kennedy, Fourth-Order 2N-Storage Runge-Kutta Schemes,
        ! NASA Technical Memorandum 109112, June 1994.
        !
        ! Online version:
        ! http://www.ece.uvic.ca/~bctill/papers/numacoust/Carpenter_Kennedy_1994.pdf

        ! Subroutine arguments
        CLASS(rk_2n_t), INTENT(out) :: rk

        ! Set to zero
        rk%a = 0.0
        rk%b = 0.0

        ! Used coefficients
        rk%nrk = 5
        rk%cflmax = 3.34 - 0.01  ! 3.34 given in paper, 0.01 round off margin
        rk%a(1:rk%nrk) = [0.0, -0.4178904745, -1.192151694643, &
            -1.697784692471, -1.514183444257]
        rk%b(1:rk%nrk) = [0.1496590219993, 0.3792103129999, 0.8229550293869, &
            0.6994504559488, 0.1530572479681]

        ! Compute C from A and B
        CALL rk%comp_c()
    END SUBROUTINE init_carpenter


    SUBROUTINE init_bernardini(rk)
        ! 2th order 5-stage scheme
        ! From: Matteo Bernardini and Sergio Pirozzoli, A general strategy for
        ! the optimization of Runge–Kutta schemes for wave propagation
        ! phenomena, Journal of Computational Physics 228 (2009)
        ! DOI: 10.1016/j.jcp.2009.02.032
        !
        ! Coefficients from table 4

        ! Subroutine arguments
        CLASS(rk_2n_t), INTENT(out) :: rk

        ! Set to zero
        rk%a = 0.0
        rk%b = 0.0

        ! Used coefficients
        rk%nrk = 5
        rk%cflmax = 3.47 - 0.01  ! 3.47 given in paper, 0.01 round off margin
        rk%a(1:rk%nrk) = [0.0, -1.0, -1.55798, -1.0, -0.45031]
        rk%b(1:rk%nrk) = [0.2, 0.83204, 0.6, 0.35394, 0.2]

        ! Compute C from A and B
        CALL rk%comp_c()
    END SUBROUTINE init_bernardini


    PURE SUBROUTINE get_coeffs_a(this, frhs, fu, dtrk, dtrki, irk)
        CLASS(rk_2n_t), INTENT(in) :: this
        REAL(realk), INTENT(out) :: frhs, fu, dtrk, dtrki
        INTEGER(intk), INTENT(in) :: irk

        IF (irk > this%nrk .OR. irk < 0) ERROR STOP

        CALL this%get_coeffs_b(frhs, fu, irk)

        IF (irk < this%nrk) THEN
            dtrk = this%c(irk+1)
        ELSE
            dtrk = 1.0
        END IF

        dtrki = dtrk - this%c(irk)
    END SUBROUTINE get_coeffs_a


    PURE SUBROUTINE get_coeffs_b(this, frhs, fu, irk)
        CLASS(rk_2n_t), INTENT(in) :: this
        REAL(realk), INTENT(out) :: frhs, fu
        INTEGER(intk), INTENT(in) :: irk

        IF (irk > this%nrk .OR. irk < 0) ERROR STOP

        frhs = this%a(irk)
        fu = this%b(irk)
    END SUBROUTINE get_coeffs_b


    ! Compute coefficient C
    !
    ! Example computation for a six-stage scheme:
    ! C(1) = 0.0
    ! C(2) = B(1)
    ! C(3) = B(1) + B(2)*(A(2) + 1)
    ! C(4) = B(1) + B(2)*(A(2) + 1) + B(3)*(A(3)*(A(2) + 1) + 1)
    ! C(5) = B(1) + B(2)*(A(2) + 1) + B(3)*(A(3)*(A(2) + 1) + 1) + B(4)*(A(4)*(A(3)*(A(2) + 1) + 1) + 1)
    ! C(6) = B(1) + B(2)*(A(2) + 1) + B(3)*(A(3)*(A(2) + 1) + 1) + B(4)*(A(4)*(A(3)*(A(2) + 1) + 1) + 1) + B(5)*(A(5)*(A(4)*(A(3)*(A(2) + 1) + 1) + 1) + 1)
    !
    ! This implementation is generally valid for any number of stages.
    !
    ! Ref: Carpenter & Kennedy, Fourth-Order 2N-Storage Runge-Kutta Schemes,
    ! NASA Technical Memorandum 109112, June 1994.
    !
    ! Online version:
    ! http://www.ece.uvic.ca/~bctill/papers/numacoust/Carpenter_Kennedy_1994.pdf
    SUBROUTINE comp_c(this)
        CLASS(rk_2n_t), INTENT(inout) :: this

        INTEGER(intk) i, j, k
        REAL(realk) :: fak

        ! Make sure C is zero
        this%c = 0.0

        ! Compute C(i)
        DO i = 1, this%nrk
            ! Compute individual terms, i.e. like B(3)*(A(3)*(A(2) + 1) + 1)
            DO j = 1, i-1
                ! Compute factors of individual terms, i.e. like (A(3)*(A(2) + 1) + 1)
                fak = 1.0
                DO k = 2, j
                    fak = this%a(k)*fak + 1.0
                END DO
                this%c(i) = this%c(i) + this%b(j)*fak
            END DO
        END DO
    END SUBROUTINE comp_c


    ! Perform an update of fields in the RK time integration scheme
    ! dU_j = A_j*dU_(j-1) + dt*uo
    ! U_j = U_(j-1) + B_j*dU_j
    PURE SUBROUTINE rkstep(p, dp, rhsp, frhs, dtfu)
        ! Subroutine arguments
        REAL(realk), CONTIGUOUS, INTENT(inout) :: p(:)
        REAL(realk), CONTIGUOUS, INTENT(inout) :: dp(:)
        REAL(realk), CONTIGUOUS, INTENT(in) :: rhsp(:)
        REAL(realk), INTENT(in) :: frhs
        REAL(realk), INTENT(in) :: dtfu

        ! Local variables
        INTEGER(intk) :: i

        ! Perform the update in a manually crafted loop is faster than using an
        ! implicit loop, because of cache effects (dp(i) is already in cache
        ! when p(i) is updated)
        DO i = 1, SIZE(p)
            dp(i) = frhs*dp(i) + rhsp(i)
            p(i) = p(i) + dtfu*dp(i)
        END DO
    END SUBROUTINE rkstep
END MODULE rungekutta_mod
