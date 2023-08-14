MODULE flowcore_mod
    USE core_mod
    USE ib_mod, ONLY: ftoc, parent, ib

    IMPLICIT NONE(type, external)
    PRIVATE

    ! Control parameters
    LOGICAL, PROTECTED :: has_flow = .FALSE.
    LOGICAL, PROTECTED :: solve_flow = .FALSE.

    ! Fluid/physical paramters
    REAL(realk), PROTECTED :: gmol
    REAL(realk), PROTECTED :: rho
    REAL(realk), PROTECTED :: uinf(3) = 0.0
    REAL(realk), PROTECTED :: tu_level
    REAL(realk), PROTECTED :: targetcflmax
    REAL(realk), PROTECTED :: gradp(3)

    ! TODO: Allocatable length - some expressions can be LONG!
    CHARACTER(len=10240), PROTECTED :: uinf_expr(3) = ""
    LOGICAL, PROTECTED :: uinf_is_expr = .FALSE.
    LOGICAL, PROTECTED :: uinf_is_time = .FALSE.

    PUBLIC :: init_flowcore, finish_flowcore, has_flow, solve_flow, gmol, &
        rho, uinf, uinf_expr, uinf_is_expr, uinf_is_time, tu_level, &
        targetcflmax, gradp

CONTAINS
    SUBROUTINE init_flowcore()
        ! Subroutine arguments
        ! None...

        ! Local variables
        TYPE(config_t) :: flowconf
        INTEGER(intk), PARAMETER :: units_v(7) = [0, 1, -1, 0, 0, 0, 0]
        INTEGER(intk), PARAMETER :: units_p(7) = [1, -1, -2, 0, 0, 0, 0]
        INTEGER(intk), PARAMETER :: units_g(7) = [1, -1, -1, 0, 0, 0, 0]

        ! Read configuration values - if not exists no timeintegration is
        ! performed
        has_flow = .FALSE.
        IF (.NOT. fort7%exists("/flow")) THEN
            IF (myid == 0) THEN
                WRITE(*, '("NO FLOW")')
                WRITE(*, '()')
            END IF
            RETURN
        END IF
        has_flow = .TRUE.

        ! Required values
        flowconf = fort7%get("/flow")
        CALL flowconf%get_value("/gmol", gmol)

        ! Either uinf is real or expression
        IF (flowconf%is_real("/uinf/0")) THEN
            CALL flowconf%get_array("/uinf", uinf)
            uinf_is_expr = .FALSE.
        ELSE IF (flowconf%is_char("/uinf/0")) THEN
            CALL flowconf%get_value("/uinf/0", uinf_expr(1))
            CALL flowconf%get_value("/uinf/1", uinf_expr(2))
            CALL flowconf%get_value("/uinf/2", uinf_expr(3))
            uinf_is_expr = .TRUE.

            ! If *any* expression contains 'timeph' or 'rand' the expression is
            ! considered time-dependent
            uinf_is_time = .FALSE.
            IF (INDEX(lower(uinf_expr(1)), "timeph") > 0) uinf_is_time = .TRUE.
            IF (INDEX(lower(uinf_expr(2)), "timeph") > 0) uinf_is_time = .TRUE.
            IF (INDEX(lower(uinf_expr(3)), "timeph") > 0) uinf_is_time = .TRUE.
            IF (INDEX(lower(uinf_expr(1)), "rand") > 0) uinf_is_time = .TRUE.
            IF (INDEX(lower(uinf_expr(2)), "rand") > 0) uinf_is_time = .TRUE.
            IF (INDEX(lower(uinf_expr(3)), "rand") > 0) uinf_is_time = .TRUE.

            IF (myid == 0) THEN
                WRITE(*, '("INITIAL CONDITION:")')
                WRITE(*, '(2X, "Using expression for U:")')
                WRITE(*, '(2X, A)') TRIM(uinf_expr(1))
                WRITE(*, '("  Using expression for V:")')
                WRITE(*, '(2X, A)') TRIM(uinf_expr(2))
                WRITE(*, '("  Using expression for W:")')
                WRITE(*, '(2X, A)') TRIM(uinf_expr(3))
                WRITE(*, '("  Expression is time-dependent: ", L1)') uinf_is_time
                WRITE(*, '()')
            END IF
        ELSE
            WRITE(*, *) "uinf must be real or character"
            CALL errr(__FILE__, __LINE__)
        END IF

        ! Optional values
        CALL flowconf%get_value("/rho", rho, 1.0)
        CALL flowconf%get_value("/tu_level", tu_level, 0.1)
        CALL flowconf%get_value("/solve", solve_flow, .TRUE.)

        ! TODO: optional gradp not working!!! Fix!
        gradp = 0.0
        CALL flowconf%get_array("/gradp", gradp, required=.FALSE.)

        CALL set_field("U", istag=1, units=units_v, dread=dread, &
            required=dread, dwrite=dwrite, buffers=.TRUE.)
        CALL set_field("V", jstag=1, units=units_v, dread=dread, &
            required=dread, dwrite=dwrite, buffers=.TRUE.)
        CALL set_field("W", kstag=1, units=units_v, dread=dread, &
            required=dread, dwrite=dwrite, buffers=.TRUE.)
        CALL set_field("P", units=units_p, dread=dread, &
            required=dread, dwrite=dwrite, buffers=.TRUE.)
        CALL set_field("G", units=units_g, dread=dread, &
            required=dread, dwrite=dwrite, buffers=.TRUE.)

        ! For RK time integration
        CALL set_field("DU", istag=1)
        CALL set_field("DV", jstag=1)
        CALL set_field("DW", kstag=1)
    END SUBROUTINE init_flowcore


    SUBROUTINE finish_flowcore
        CONTINUE
    END SUBROUTINE finish_flowcore

END MODULE flowcore_mod
