MODULE flow_mod
    ! Do not 'USE' core_mod or ib_mod here, because they will get exported
    ! (this module had no 'PRIVATE' default)

    USE flowcore_mod
    USE flowstat_mod
    USE timeintegration_mod
    USE wernerwengle_mod
    USE lesmodel_mod, ONLY: ilesmodel

    IMPLICIT NONE(type, external)

    PRIVATE :: init_uvwp

CONTAINS
    SUBROUTINE init_flow()
        ! These symbold do not need to be exported from this module
        USE core_mod
        USE gc_compbodyforce_mod, ONLY: init_compbodyforce
        USE gc_flowstencils_mod
        USE ib_mod
        USE lesmodel_mod
        USE pressuresolver_mod
        USE itinfo_mod, ONLY: init_itinfo
        USE boussinesqterm_mod, ONLY: init_boussinesqterm
        USE coriolisterm_mod, ONLY: init_coriolisterm

        ! Local variables
        TYPE(field_t), POINTER :: u, v, w
        TYPE(field_t), POINTER :: pwu, pwv, pww
        TYPE(field_t), POINTER :: pwub, pwvb, pwwb

        CALL init_flowcore()
        IF (.NOT. has_flow) RETURN

        CALL init_uvwp()
        CALL init_flowstat()

        ! Wall- and LES models are needed for pure scalar simulation
        CALL init_wernerwengle()
        CALL init_lesmodel()

        ! Werner wengle and LES models are often used by other models (scalar)
        ! even when no flow is solved.
        IF (.NOT. solve_flow) RETURN

        ! These sould only be needed when flow is actually solved
        CALL set_timer(300, "FLOW")
        CALL set_timer(310, "FLOW_TSTLE4")
        CALL set_timer(320, "FLOW_MGPOISL")
        CALL set_timer(321, "FLOW_MGPOISIT")
        CALL set_timer(322, "FLOW_MGPOISL_INNER")
        CALL set_timer(330, "FLOW_LESMODEL")
        CALL set_timer(340, "FLOW_SETIBVALUES")
        CALL set_timer(341, "FLOW_GETIBVALUES")
        CALL set_timer(342, "FLOW_SETPOINTVALUES")
        CALL set_timer(350, "FLOW_ITINFO")
        CALL set_timer(351, "FLOW_COMPBODYFORCE")
        CALL set_timer(360, "FLOW_BOUSSINESQTERM")
        CALL set_timer(370, "FLOW_CORIOLISTERM")

        CALL init_pressuresolver()
        CALL init_boussinesqterm()
        CALL init_coriolisterm()
        CALL init_itinfo(dcont)

        ! Need to call this here - cannot be in flowcore because that
        ! create a circular dependency
        SELECT TYPE(ib)
        TYPE IS (gc_t)
            CALL create_flowstencils(ib)
            IF (compbodyforce) CALL init_compbodyforce(dcont)

            CALL get_field(u, "U")
            CALL get_field(v, "V")
            CALL get_field(w, "W")

            ! Fields for the point values (german: PunktWert = PW)
            CALL set_field("PWU", istag=1, buffers=.TRUE., &
                units=u%units, dwrite=writepntvalues)
            CALL set_field("PWV", jstag=1, buffers=.TRUE., &
                units=v%units, dwrite=writepntvalues)
            CALL set_field("PWW", kstag=1, buffers=.TRUE., &
                units=w%units, dwrite=writepntvalues)

            IF (writepntvalues) THEN
                ! Inidicator field if point values are present
                CALL set_field("PWUB", istag=1, dwrite=.TRUE.)
                CALL set_field("PWVB", jstag=1, dwrite=.TRUE.)
                CALL set_field("PWWB", kstag=1, dwrite=.TRUE.)
                CALL get_field(pwub, "PWUB")
                CALL get_field(pwvb, "PWVB")
                CALL get_field(pwwb, "PWWB")
                CALL setpointvaluemarkers(pwub, pwvb, pwwb)
            END IF

            CALL get_field(pwu, "PWU")
            CALL get_field(pwv, "PWV")
            CALL get_field(pww, "PWW")
            CALL setpointvalues(pwu, pwv, pww, u, v, w, .TRUE.)
            CALL setibvalues(u, v, w)
        END SELECT
    END SUBROUTINE init_flow


    SUBROUTINE finish_flow
        ! These symbold do not need to be exported from this module
        USE lesmodel_mod
        USE pressuresolver_mod
        USE itinfo_mod, ONLY: finish_itinfo
        USE gc_compbodyforce_mod, ONLY: finish_compbodyforce
        USE gc_flowstencils_mod
        USE ib_mod
        USE boussinesqterm_mod, ONLY: finish_boussinesqterm
        USE coriolisterm_mod, ONLY: finish_coriolisterm

        IF (.NOT. has_flow) RETURN

        IF (solve_flow) THEN
            SELECT TYPE(ib)
            TYPE IS (gc_t)
                IF (solve_flow) THEN
                    IF (compbodyforce) CALL finish_compbodyforce()
                    CALL finish_flowstencils()
                END IF
            END SELECT

            CALL finish_itinfo
            CALL finish_boussinesqterm()
            CALL finish_coriolisterm()
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

        TYPE(field_t), POINTER :: u, v, w, p
        INTEGER(intk) :: ilevel

        IF (ib%type /= "NONE") THEN
            CALL init_fix_openfraction()
        END IF

        CALL get_field(u, "U")
        CALL get_field(v, "V")
        CALL get_field(w, "W")
        CALL get_field(p, "P")

        ! Set inflow buffers for FIX bc's
        DO ilevel = minlevel, maxlevel
            CALL setboundarybuffers%bound(ilevel, u, v, w, timeph=0.0_realk)
        END DO

        ! The rest of this routine is initializing the flow field - we do not
        ! want to overwrite results read from a restart file
        IF (dread) RETURN

        ! Set initial condition
        IF (uinf_is_expr) THEN
            CALL init_uvw_expr(u, v, w)
        ELSE
            CALL init_uvw_uinf(u, v, w)
        END IF
        p = 0.0_realk

        CALL zero_ghostlayers(u)
        CALL zero_ghostlayers(v)
        CALL zero_ghostlayers(w)

        DO ilevel = minlevel, maxlevel
            CALL connect(ilevel, 2, u, v, w, p, corners=.TRUE.)
        END DO

        DO ilevel = minlevel+1, maxlevel
            CALL parent(ilevel, u, v, w, p)
            CALL bound_flow%bound(ilevel, u, v, w, p)
        END DO

        DO ilevel = maxlevel, minlevel+1, -1
            CALL ftoc(ilevel, u%arr, u%arr, 'U')
            CALL ftoc(ilevel, v%arr, v%arr, 'V')
            CALL ftoc(ilevel, w%arr, w%arr, 'W')
            CALL ftoc(ilevel, p%arr, p%arr, 'P')
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
                0.0_realk, xstag, y, z, dx, dy, dz, ddx, ddy, ddz)
            CALL initial_condition(v, "v", uinf_expr(2), rho, gmol, tu_level, &
                0.0_realk, x, ystag, z, dx, dy, dz, ddx, ddy, ddz)
            CALL initial_condition(w, "w", uinf_expr(3), rho, gmol, tu_level, &
                0.0_realk, x, y, zstag, dx, dy, dz, ddx, ddy, ddz)
        END DO

    END SUBROUTINE init_uvw_expr


    ! This subroutine compute the open fraction for the FIX boundary condition
    ! used to create a divergence-free inflow condition also when an immersed
    ! boundary intersect the inflow
    SUBROUTINE init_fix_openfraction()
        USE core_mod
        USE ib_mod

        TYPE(field_t) :: bucoarse, bvcoarse, bwcoarse
        TYPE(field_t), POINTER :: bu, bv, bw
        INTEGER(intk) :: i, igrid, iface, ipic, ilevel
        INTEGER(intk) :: kk, jj, ii
        CHARACTER(len=8) :: ctyp
        REAL(realk), POINTER, CONTIGUOUS :: buf(:, :, :), field(:, :, :)

        CALL get_field(bu, "BU")
        CALL get_field(bv, "BV")
        CALL get_field(bw, "BW")

        CALL bu%init_buffers()
        CALL bv%init_buffers()
        CALL bw%init_buffers()

        CALL bucoarse%copy_from(bu)
        CALL bvcoarse%copy_from(bv)
        CALL bwcoarse%copy_from(bw)

        DO ilevel = maxlevel, minlevel+1, -1
            CALL ftoc(ilevel, bucoarse%arr, bucoarse%arr, 'U')
            CALL ftoc(ilevel, bvcoarse%arr, bvcoarse%arr, 'V')
            CALL ftoc(ilevel, bwcoarse%arr, bwcoarse%arr, 'W')
        END DO

        DO i = 1, nmygrids
            igrid = mygrids(i)
            CALL get_mgdims(kk, jj, ii, igrid)

            DO iface = 1, 6
                ! It is assumed that FIX is in the first position ibocd = 1
                CALL get_bc_ctyp(ctyp, 1, iface, igrid)
                IF (ctyp /= "FIX") CYCLE

                SELECT CASE(iface)
                CASE(1, 3, 5)
                    ipic = 2
                CASE(2)
                    ipic = ii-2
                CASE(4)
                    ipic = jj-2
                CASE(6)
                    ipic = kk-2
                END SELECT

                ! Both slices in buf are set to the same value, makes it easier
                ! to multiply correctly
                SELECT CASE(iface)
                CASE(1, 2)
                    CALL bu%buffers%get_buffer(buf, igrid, iface)
                    CALL bucoarse%get_ptr(field, igrid)
                    buf(:, :, 1) = field(:, :, ipic)
                    buf(:, :, 2) = field(:, :, ipic)
                CASE(3, 4)
                    CALL bv%buffers%get_buffer(buf, igrid, iface)
                    CALL bvcoarse%get_ptr(field, igrid)
                    buf(:, :, 1) = field(:, ipic, :)
                    buf(:, :, 2) = field(:, ipic, :)
                CASE(5, 6)
                    CALL bw%buffers%get_buffer(buf, igrid, iface)
                    CALL bwcoarse%get_ptr(field, igrid)
                    buf(:, :, 1) = field(ipic, :, :)
                    buf(:, :, 2) = field(:, ipic, :)
                END SELECT
            END DO
        END DO

        CALL bucoarse%finish()
        CALL bvcoarse%finish()
        CALL bwcoarse%finish()
    END SUBROUTINE init_fix_openfraction
END MODULE flow_mod
