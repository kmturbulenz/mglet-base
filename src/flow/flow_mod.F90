MODULE flow_mod
    ! Do not 'USE' core_mod or ib_mod here, because they will get exported
    ! (this module had no 'PRIVATE' default)

    USE flowcore_mod
    USE flowstat_mod
    USE timeintegration_mod
    USE wernerwengle_mod

    IMPLICIT NONE(type, external)

    PRIVATE :: init_uvwp

CONTAINS
    SUBROUTINE init_flow()
        ! These symbold do not need to be exported from this module
        USE core_mod
        USE gc_flowstencils_mod
        USE ib_mod
        USE lesmodel_mod
        USE pressuresolver_mod
        USE itinfo_mod, ONLY: init_itinfo

        ! Local variables
        TYPE(field_t), POINTER :: u, v, w
        TYPE(field_t), POINTER :: pwu, pwv, pww, sdiv

        CALL init_flowcore()
        IF (.NOT. has_flow) RETURN

        IF (.NOT. dread) THEN
            CALL init_uvwp()
        END IF

        CALL init_flowstat()

        ! Wall- and LES models are needed for pure scalar simulation
        CALL init_wernerwengle()
        CALL init_lesmodel()

        ! These sould only be needed when flow is actually solved
        IF (solve_flow) THEN
            CALL init_pressuresolver()
            CALL init_itinfo(dcont)

            ! Need to call this here - cannot be in flowcore because that
            ! create a circular dependency
            SELECT TYPE(ib)
            TYPE IS (gc_t)
                CALL create_flowstencils(ib)

                CALL set_field("PWU", istag=1, buffers=.TRUE.)
                CALL set_field("PWV", jstag=1, buffers=.TRUE.)
                CALL set_field("PWW", kstag=1, buffers=.TRUE.)

                CALL get_field(pwu, "PWU")
                CALL get_field(pwv, "PWV")
                CALL get_field(pww, "PWW")
                CALL get_field(u, "U")
                CALL get_field(v, "V")
                CALL get_field(w, "W")
                CALL setpointvalues(pwu, pwv, pww, u, v, w, .TRUE.)
                CALL setibvalues(u, v, w)

                CALL get_field(sdiv, "SDIV")
                CALL setsdivfield(sdiv)
            END SELECT
        END IF
    END SUBROUTINE init_flow


    SUBROUTINE finish_flow
        ! These symbold do not need to be exported from this module
        USE lesmodel_mod
        USE pressuresolver_mod
        USE itinfo_mod, ONLY: finish_itinfo
        USE gc_flowstencils_mod
        USE ib_mod

        IF (.NOT. has_flow) RETURN

        IF (solve_flow) THEN
            SELECT TYPE(ib)
            TYPE IS (gc_t)
                IF (solve_flow) THEN
                    CALL finish_flowstencils()
                END IF
            END SELECT

            CALL finish_itinfo
            CALL finish_pressuresolver()
        END IF

        CALL finish_lesmodel()
        CALL finish_wernerwengle()
        CALL finish_flowcore()
    END SUBROUTINE finish_flow


    SUBROUTINE init_uvwp()
        USE bound_flow_mod
        USE core_mod
        USE ib_mod
        USE setboundarybuffers_mod, ONLY: setboundarybuffers

        TYPE(field_t), POINTER :: u_f, v_f, w_f, p_f
        REAL(realk), POINTER, CONTIGUOUS :: u(:), v(:), w(:), p(:)
        INTEGER(intk) :: ilevel

        CALL get_field(u_f, "U")
        CALL get_field(v_f, "V")
        CALL get_field(w_f, "W")
        CALL get_field(p_f, "P")

        ! Associate pinters
        u => u_f%arr
        v => v_f%arr
        w => w_f%arr
        p => p_f%arr

        ! Set inflow buffers for FIX bc's
        DO ilevel = minlevel, maxlevel
            CALL setboundarybuffers%bound(ilevel, u_f, v_f, w_f, &
                timeph=0.0_realk)
        END DO

        ! Set initial condition
        IF (uinf_is_expr) THEN
            CALL init_uvw_expr(u_f, v_f, w_f)
        ELSE
            CALL init_uvw_uinf(u_f, v_f, w_f)
        END IF
        p = 0.0

        CALL zero_ghostlayers(u_f)
        CALL zero_ghostlayers(v_f)
        CALL zero_ghostlayers(w_f)

        DO ilevel = minlevel, maxlevel
            CALL connect(ilevel, 2, u, v, w, p, corners=.TRUE.)
        END DO

        DO ilevel = minlevel+1, maxlevel
            CALL parent(ilevel, u_f, v_f, w_f, p_f)
            CALL bound_flow%bound(ilevel, u_f, v_f, w_f, p_f)
        END DO

        DO ilevel = maxlevel, minlevel+1, -1
            CALL ftoc(ilevel, u, u, 'U')
            CALL ftoc(ilevel, v, v, 'V')
            CALL ftoc(ilevel, w, w, 'W')
            CALL ftoc(ilevel, p, p, 'P')
        END DO
    END SUBROUTINE init_uvwp


    SUBROUTINE init_uvw_uinf(u_f, v_f, w_f)
        USE core_mod

        ! Subroutine arguments
        TYPE(field_t), INTENT(inout) :: u_f
        TYPE(field_t), INTENT(inout) :: v_f
        TYPE(field_t), INTENT(inout) :: w_f

        ! Local variables
        REAL(realk) :: uinf_mag

        ! Initialize fields - random noise between 0.0 .. 1.0
        CALL RANDOM_NUMBER(u_f%arr)
        CALL RANDOM_NUMBER(v_f%arr)
        CALL RANDOM_NUMBER(w_f%arr)

        uinf_mag = SQRT(uinf(1)*uinf(1) + uinf(2)*uinf(2) + uinf(3)*uinf(3))

        ! Subtract 0.5 from random noise to have a noise field with zero mean
        ! value
        u_f%arr = (u_f%arr-0.5)*tu_level*uinf_mag + uinf(1)
        v_f%arr = (v_f%arr-0.5)*tu_level*uinf_mag + uinf(2)
        w_f%arr = (w_f%arr-0.5)*tu_level*uinf_mag + uinf(3)
    END SUBROUTINE init_uvw_uinf


    SUBROUTINE init_uvw_expr(u_f, v_f, w_f)
        USE core_mod

        ! Subroutine arguments
        TYPE(field_t), INTENT(inout) :: u_f
        TYPE(field_t), INTENT(inout) :: v_f
        TYPE(field_t), INTENT(inout) :: w_f

        ! Local variables
        INTEGER(intk) :: i, igrid
        REAL(realk), POINTER, CONTIGUOUS :: u(:, :, :), v(:, :, :), w(:, :, :)
        REAL(realk), POINTER, CONTIGUOUS :: x(:), y(:), z(:)
        REAL(realk), POINTER, CONTIGUOUS :: dx(:), dy(:), dz(:)
        REAL(realk), POINTER, CONTIGUOUS :: ddx(:), ddy(:), ddz(:)
        REAL(realk), POINTER, CONTIGUOUS :: xstag(:), ystag(:), zstag(:)

        DO i = 1, nmygrids
            igrid = mygrids(i)

            CALL u_f%get_ptr(u, igrid)
            CALL v_f%get_ptr(v, igrid)
            CALL w_f%get_ptr(w, igrid)

            CALL get_fieldptr(x, "X", igrid)
            CALL get_fieldptr(y, "Y", igrid)
            CALL get_fieldptr(z, "Z", igrid)

            CALL get_fieldptr(dx, "DX", igrid)
            CALL get_fieldptr(dy, "DY", igrid)
            CALL get_fieldptr(dz, "DZ", igrid)

            CALL get_fieldptr(ddx, "DDX", igrid)
            CALL get_fieldptr(ddy, "DDY", igrid)
            CALL get_fieldptr(ddz, "DDZ", igrid)

            CALL get_fieldptr(xstag, "XSTAG", igrid)
            CALL get_fieldptr(ystag, "YSTAG", igrid)
            CALL get_fieldptr(zstag, "ZSTAG", igrid)

            CALL initial_condition(u, "u", uinf_expr(1), rho, gmol, tu_level, &
                0.0_realk, xstag, y, z, dx, dy, dz,  ddx, ddy, ddz)
            CALL initial_condition(v, "v", uinf_expr(2), rho, gmol, tu_level, &
                0.0_realk, x, ystag, z, dx, dy, dz,  ddx, ddy, ddz)
            CALL initial_condition(w, "w", uinf_expr(3), rho, gmol, tu_level, &
                0.0_realk, x, y, zstag, dx, dy, dz,  ddx, ddy, ddz)
        END DO

    END SUBROUTINE init_uvw_expr
END MODULE flow_mod
