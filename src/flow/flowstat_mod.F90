MODULE flowstat_mod
    USE core_mod
    USE flowcore_mod

    IMPLICIT NONE(type, external)
    PRIVATE

    PUBLIC :: init_flowstat, finish_flowstat

CONTAINS
    SUBROUTINE init_flowstat()
        ! Subroutine arguments
        ! none...

        ! Local variables
        ! none....

        CALL register_statfield("U_AVG", comp_avg)
        CALL register_statfield("V_AVG", comp_avg)
        CALL register_statfield("W_AVG", comp_avg)
        CALL register_statfield("P_AVG", comp_avg)
        CALL register_statfield("G_AVG", comp_avg)

        CALL register_statfield("UU_AVG", comp_sqr_avg)
        CALL register_statfield("VV_AVG", comp_sqr_avg)
        CALL register_statfield("WW_AVG", comp_sqr_avg)
        CALL register_statfield("PP_AVG", comp_sqr_avg)
        CALL register_statfield("GG_AVG", comp_sqr_avg)

        CALL register_statfield("UV_AVG", comp_uv_avg)
        CALL register_statfield("UW_AVG", comp_uv_avg)
        CALL register_statfield("VW_AVG", comp_uv_avg)

        CALL register_statfield("UUU_AVG", comp_cube_avg)
        CALL register_statfield("VVV_AVG", comp_cube_avg)
        CALL register_statfield("WWW_AVG", comp_cube_avg)

        CALL register_statfield("UUV_AVG", comp_uvw_avg)
        CALL register_statfield("UUW_AVG", comp_uvw_avg)
        CALL register_statfield("UVV_AVG", comp_uvw_avg)
        CALL register_statfield("UVW_AVG", comp_uvw_avg)
        CALL register_statfield("UWW_AVG", comp_uvw_avg)
        CALL register_statfield("VVW_AVG", comp_uvw_avg)
        CALL register_statfield("VWW_AVG", comp_uvw_avg)

        CALL register_statfield("laplaceP_AVG", comp_laplacep_avg)
        CALL register_statfield("laplaceP_SQR_AVG", comp_laplacep_sqr_avg)

        CALL register_statfield("DISSIP_AVG", comp_dissip_avg)

        CALL register_statfield("UP_AVG", comp_up_avg)
        CALL register_statfield("VP_AVG", comp_up_avg)
        CALL register_statfield("WP_AVG", comp_up_avg)

        CALL register_statfield("UxP_AVG", comp_uxp_avg)
        CALL register_statfield("VyP_AVG", comp_uxp_avg)
        CALL register_statfield("WzP_AVG", comp_uxp_avg)

        CALL register_statfield("UxUx_AVG", comp_uxux_avg)
        CALL register_statfield("UyUy_AVG", comp_uxux_avg)
        CALL register_statfield("UzUz_AVG", comp_uxux_avg)
        CALL register_statfield("VxVx_AVG", comp_uxux_avg)
        CALL register_statfield("VyVy_AVG", comp_uxux_avg)
        CALL register_statfield("VzVz_AVG", comp_uxux_avg)
        CALL register_statfield("WxWx_AVG", comp_uxux_avg)
        CALL register_statfield("WyWy_AVG", comp_uxux_avg)
        CALL register_statfield("WzWz_AVG", comp_uxux_avg)

        CALL register_statfield("UxVx_AVG", comp_uxvx_avg)
        CALL register_statfield("UxWx_AVG", comp_uxvx_avg)
        CALL register_statfield("UyVy_AVG", comp_uxvx_avg)
        CALL register_statfield("UyWy_AVG", comp_uxvx_avg)
        CALL register_statfield("UzVz_AVG", comp_uxvx_avg)
        CALL register_statfield("UzWz_AVG", comp_uxvx_avg)
        CALL register_statfield("VxWx_AVG", comp_uxvx_avg)
        CALL register_statfield("VyWy_AVG", comp_uxvx_avg)
        CALL register_statfield("VzWz_AVG", comp_uxvx_avg)

        CALL register_statfield("UyP+VxP_AVG", comp_uyvxp_avg)
        CALL register_statfield("UzP+WxP_AVG", comp_uyvxp_avg)
        CALL register_statfield("VzP+WyP_AVG", comp_uyvxp_avg)
    END SUBROUTINE init_flowstat


    SUBROUTINE finish_flowstat
        CONTINUE
    END SUBROUTINE finish_flowstat


    ! Routine to compute the UV_AVG, UW_AVG and VW_AVG fields
    SUBROUTINE comp_uv_avg(field, name, dt)
        ! Subroutine arguments
        TYPE(field_t), INTENT(inout) :: field
        CHARACTER(len=*), INTENT(in) :: name
        REAL(realk), INTENT(in) :: dt

        ! Local variables
        INTEGER(intk), PARAMETER :: units(*) = [0, 2, -2, 0, 0, 0, 0]
        TYPE(field_t), POINTER :: in1, in2
        INTEGER(intk) :: istag, jstag, kstag

        SELECT CASE (TRIM(name))
        CASE ("UV_AVG")
            CALL get_field(in1, "U")
            CALL get_field(in2, "V")
            istag = 1
            jstag = 1
            kstag = 0
        CASE ("UW_AVG")
            CALL get_field(in1, "U")
            CALL get_field(in2, "W")
            istag = 1
            jstag = 0
            kstag = 1
        CASE ("VW_AVG")
            CALL get_field(in1, "V")
            CALL get_field(in2, "W")
            istag = 0
            jstag = 1
            kstag = 1
        CASE DEFAULT
            CALL errr(__FILE__, __LINE__)
        END SELECT

        CALL field%init(name, istag=istag, jstag=jstag, kstag=kstag, &
            units=units)
        CALL field%multiply(in1, in2)
    END SUBROUTINE comp_uv_avg


    ! Routine to compute the UUV_AVG, UUW_AVG, UVV_AVG,
    ! UVW_AVG, UWW_AVG, VVW_AVG and VWW_AVG fields
    SUBROUTINE comp_uvw_avg(field, name, dt)
        ! Subroutine arguments
        TYPE(field_t), INTENT(inout) :: field
        CHARACTER(len=*), INTENT(in) :: name
        REAL(realk), INTENT(in) :: dt

        ! Local variables
        INTEGER(intk), PARAMETER :: units(*) = [0, 3, -3, 0, 0, 0, 0]
        TYPE(field_t), POINTER :: in1, in2, in3
        INTEGER(intk) :: istag, jstag, kstag

        SELECT CASE (TRIM(name))
        CASE ("UUV_AVG")
            CALL get_field(in1, "U")
            CALL get_field(in2, "U")
            CALL get_field(in3, "V")
            istag = 1
            jstag = 1
            kstag = 0
        CASE ("UUW_AVG")
            CALL get_field(in1, "U")
            CALL get_field(in2, "U")
            CALL get_field(in3, "W")
            istag = 1
            jstag = 0
            kstag = 1
        CASE ("UVV_AVG")
            CALL get_field(in1, "U")
            CALL get_field(in2, "V")
            CALL get_field(in3, "V")
            istag = 1
            jstag = 1
            kstag = 0
        CASE ("UVW_AVG")
            CALL get_field(in1, "U")
            CALL get_field(in2, "V")
            CALL get_field(in3, "W")
            istag = 1
            jstag = 1
            kstag = 1
        CASE ("UWW_AVG")
            CALL get_field(in1, "U")
            CALL get_field(in2, "W")
            CALL get_field(in3, "W")
            istag = 1
            jstag = 0
            kstag = 1
        CASE ("VVW_AVG")
            CALL get_field(in1, "V")
            CALL get_field(in2, "V")
            CALL get_field(in3, "W")
            istag = 0
            jstag = 1
            kstag = 1
        CASE ("VWW_AVG")
            CALL get_field(in1, "V")
            CALL get_field(in2, "W")
            CALL get_field(in3, "W")
            istag = 0
            jstag = 1
            kstag = 1
        CASE DEFAULT
            CALL errr(__FILE__, __LINE__)
        END SELECT

        CALL field%init(name, istag=istag, jstag=jstag, kstag=kstag, &
            units=units)
        CALL field%multiply(in1, in2, in3)
    END SUBROUTINE comp_uvw_avg


    SUBROUTINE comp_laplacep_avg(field, name, dt)
        ! Subroutine arguments
        TYPE(field_t), INTENT(inout) :: field
        CHARACTER(len=*), INTENT(in) :: name
        REAL(realk), INTENT(in) :: dt

        ! Local variables
        TYPE(field_t), POINTER :: p_f, bp_f, dx_f, dy_f, dz_f
        REAL(realk), CONTIGUOUS, POINTER :: p(:, :, :), bp(:, :, :), &
            lpp(:, :, :), dx(:), dy(:), dz(:)
        INTEGER(intk), PARAMETER :: units(*) = [1, -3, -2, 0, 0, 0, 0]
        INTEGER(intk) :: i, igrid
        INTEGER(intk) :: nfro, nbac, nlft, nrgt, ntop, nbot
        INTEGER(intk) :: kk, jj, ii

        IF (name /= "laplaceP_AVG") CALL errr(__FILE__, __LINE__)

        CALL get_field(p_f, "P")
        CALL get_field(bp_f, "BP")
        CALL get_field(dx_f, "DX")
        CALL get_field(dy_f, "DY")
        CALL get_field(dz_f, "DZ")

        ! Create a field to store the result, this creates an empty field
        CALL field%init(name, units=units)

        ! Compute laplaceP
        DO i = 1, nmygrids
            igrid = mygrids(i)

            CALL field%get_ptr(lpp, igrid)
            CALL p_f%get_ptr(p, igrid)
            CALL bp_f%get_ptr(bp, igrid)
            CALL dx_f%get_ptr(dx, igrid)
            CALL dy_f%get_ptr(dy, igrid)
            CALL dz_f%get_ptr(dz, igrid)
            CALL get_mgbasb(nfro, nbac, nrgt, nlft, nbot, ntop, igrid)
            CALL get_mgdims(kk, jj, ii, igrid)

            CALL calclpp_grid(kk, jj, ii, lpp, p, bp, dx, dy, dz, &
                nfro, nbac, nrgt, nlft, nbot, ntop)
        END DO
    END SUBROUTINE comp_laplacep_avg


    SUBROUTINE comp_laplacep_sqr_avg(field, name, dt)
        ! Subroutine arguments
        TYPE(field_t), INTENT(inout) :: field
        CHARACTER(len=*), INTENT(in) :: name
        REAL(realk), INTENT(in) :: dt

        ! Local variables
        ! nobe...

        IF (name /= "laplaceP_SQR_AVG") CALL errr(__FILE__, __LINE__)

        CALL comp_laplacep_avg(field, "laplaceP_AVG", dt)
        field%arr = field%arr(:)**2
        field%name = "laplaceP_SQR_AVG"
        field%units = field%units*2
    END SUBROUTINE comp_laplacep_sqr_avg


    SUBROUTINE calclpp_grid(kk, jj, ii, lpp, p, bp, dx, dy, dz, &
            nfro, nbac, nrgt, nlft, nbot, ntop)
        ! Subroutine arguments
        INTEGER(intk), INTENT(in) :: kk, jj, ii
        REAL(realk), INTENT(inout) :: lpp(kk, jj, ii)
        REAL(realk), INTENT(in) :: p(kk, jj, ii)
        REAL(realk), INTENT(in) :: bp(kk, jj, ii)
        REAL(realk), INTENT(in) :: dx(ii), dy(jj), dz(kk)
        INTEGER(intk), INTENT(in) :: nfro, nbac, nrgt, nlft, nbot, ntop

        ! Local variables
        INTEGER(intk) :: k, j, i
        INTEGER(intk) :: ifr, iba, jri, jle, kbo, kto
        REAL(realk) :: flfr, flba, flri, flle, flto, flbo
        REAL(realk) :: d2pdx2, d2pdy2, d2pdz2

        ! 2 = FIX, 5 = NOS, 6 = SLI, 19 = CO1
        IF (nfro == 2 .OR. nfro == 5 .OR. nfro == 6 .OR. nfro == 19) ifr = 1
        IF (nrgt == 2 .OR. nrgt == 5 .OR. nrgt == 6 .OR. nrgt == 19) jri = 1
        IF (nbot == 2 .OR. nbot == 5 .OR. nbot == 6 .OR. nbot == 19) kbo = 1
        IF (nbac == 2 .OR. nbac == 5 .OR. nbac == 6 .OR. nbac == 19) iba = 1
        IF (nlft == 2 .OR. nlft == 5 .OR. nlft == 6 .OR. nlft == 19) jle = 1
        IF (ntop == 2 .OR. ntop == 5 .OR. ntop == 6 .OR. ntop == 19) kto = 1

        ! lpp is declared INTENT(inout), therefore there are no need to set
        ! buffers (where no lpp is computed) to zero. If it were INTENT(out)
        ! the entire field should have been completely defined inside this
        ! routine

        DO i = 3, ii-2
            IF (ifr == 1 .AND. i == 3) THEN
                flfr = 1.0
            ELSE
                flfr = 0.0
            END IF
            IF (iba == 1 .AND. i == ii-2) THEN
                flba = 1.0
            ELSE
                flba = 0.0
            END IF

            DO j = 3, jj-2
                IF (jri == 1 .AND. j == 3) THEN
                    flri = 1.0
                ELSE
                    flri = 0.0
                END IF
                IF (jle == 1 .AND. j == jj-2) THEN
                    flle = 1.0
                ELSE
                    flle = 0.0
                END IF

                DO k = 3, kk-2
                    IF (kbo == 1 .AND. k == 3) THEN
                        flbo = 1.0
                    ELSE
                        flbo = 0.0
                    END IF
                    IF(kto == 1 .AND. k == kk-2) THEN
                        flto = 1.0
                    ELSE
                        flto = 0.0
                    END IF

                    d2pdx2 = ((1.0-flba)*bp(k, j, i+1)*p(k, j, i+1) &
                        + (-bp(k, j, i+1)-bp(k, j, i-1)+flfr+flba)*p(k, j, i) &
                        + (1.0-flfr)*bp(k, j, i-1)*p(k, j, i-1)) &
                        /(dx(i)*dx(i-1))

                    d2pdy2 = ((1.0-flle)*bp(k, j+1, i)*p(k, j+1, i) &
                        + (-bp(k, j+1, i)-bp(k, j-1, i)+flri+flle)*p(k, j, i) &
                        + (1.0-flri)*bp(k, j-1, i)*p(k, j-1, i)) &
                        /(dy(j)*dy(j-1))

                    d2pdz2 = ((1.0-flto)*bp(k+1, j, i)*p(k+1, j, i) &
                        + (-bp(k+1, j, i)-bp(k-1, j, i)+flto+flbo)*p(k, j, i) &
                        + (1.0-flbo)*bp(k-1, j, i)*p(k-1, j, i)) &
                        /(dz(k)*dz(k-1))

                    lpp(k, j, i) = bp(k, j, i)*(d2pdx2 + d2pdy2 + d2pdz2)
                END DO
            END DO
        END DO
    END SUBROUTINE calclpp_grid


    SUBROUTINE comp_dissip_avg(field, name, dt)
        ! Subroutine arguments
        TYPE(field_t), INTENT(inout) :: field
        CHARACTER(len=*), INTENT(in) :: name
        REAL(realk), INTENT(in) :: dt

        ! Local Variables
        TYPE(field_t), POINTER :: u_f, v_f, w_f, bp_f, g_f
        TYPE(field_t), POINTER :: dx_f, dy_f, dz_f
        TYPE(field_t), POINTER :: ddx_f, ddy_f, ddz_f
        TYPE(field_t), POINTER :: rddx_f, rddy_f, rddz_f

        REAL(realk), POINTER, CONTIGUOUS :: u(:, :, :), v(:, :, :), &
            w(:, :, :), dfg(:, :, :), bp(:, :, :), g(:, :, :)
        REAL(realk), POINTER, CONTIGUOUS :: dx(:), dy(:), dz(:)
        REAL(realk), POINTER, CONTIGUOUS :: ddx(:), ddy(:), ddz(:)
        REAL(realk), POINTER, CONTIGUOUS :: rddx(:), rddy(:), rddz(:)

        INTEGER(intk), PARAMETER :: units(*) = [0, 2, -3, 0, 0, 0, 0]
        INTEGER(intk) :: i, igrid
        INTEGER(intk) :: kk, jj, ii

        IF (name /= "DISSIP_AVG") CALL errr(__FILE__, __LINE__)

        CALL get_field(u_f, "U")
        CALL get_field(v_f, "V")
        CALL get_field(w_f, "W")
        CALL get_field(bp_f, "BP")

        ! "G" is always defined and holds "gmol" for DNS
        CALL get_field(g_f, "G")

        CALL get_field(dx_f, "DX")
        CALL get_field(dy_f, "DY")
        CALL get_field(dz_f, "DZ")

        CALL get_field(ddx_f, "DDX")
        CALL get_field(ddy_f, "DDY")
        CALL get_field(ddz_f, "DDZ")

        CALL get_field(rddx_f, "RDDX")
        CALL get_field(rddy_f, "RDDY")
        CALL get_field(rddz_f, "RDDZ")

        ! Create a field to store the result, this creates an empty field
        CALL field%init(name, units=units)

        ! Compute Dissipation (It considers only molecular viscosity so far)
        DO i = 1, nmygrids
            igrid = mygrids(i)

            CALL field%get_ptr(dfg, igrid)

            CALL u_f%get_ptr(u, igrid)
            CALL v_f%get_ptr(v, igrid)
            CALL w_f%get_ptr(w, igrid)
            CALL bp_f%get_ptr(bp, igrid)
            CALL g_f%get_ptr(g, igrid)

            CALL dx_f%get_ptr(dx, igrid)
            CALL dy_f%get_ptr(dy, igrid)
            CALL dz_f%get_ptr(dz, igrid)

            CALL ddx_f%get_ptr(ddx, igrid)
            CALL ddy_f%get_ptr(ddy, igrid)
            CALL ddz_f%get_ptr(ddz, igrid)

            CALL rddx_f%get_ptr(rddx, igrid)
            CALL rddy_f%get_ptr(rddy, igrid)
            CALL rddz_f%get_ptr(rddz, igrid)

            CALL get_mgdims(kk, jj, ii, igrid)
            CALL calc_dissip(kk, jj, ii, dfg, u, v, w, bp, g, &
                dx, dy, dz, ddx, ddy, ddz, rddx, rddy, rddz)
        END DO
    END SUBROUTINE comp_dissip_avg


    SUBROUTINE calc_dissip(kk, jj, ii, dfg, u, v, w, bp, g, &
        dx, dy, dz, ddx, ddy, ddz, rddx, rddy, rddz)

        ! Subroutine arguments
        INTEGER(intk), INTENT(in) :: kk, jj, ii
        REAL(realk), INTENT(inout) :: dfg (kk, jj, ii)
        REAL(realk), INTENT(in) :: u(kk, jj, ii), v(kk, jj, ii), w(kk, jj, ii)
        REAL(realk), INTENT(in) :: bp(kk, jj, ii)
        REAL(realk), INTENT(in) :: g(kk, jj, ii)
        REAL(realk), INTENT(in) :: dx(ii), dy(jj), dz(kk)
        REAL(realk), INTENT(in) :: ddx(ii), ddy(jj), ddz(kk)
        REAL(realk), INTENT(in) :: rddx(ii), rddy(jj), rddz(kk)

        ! Local Variables
        INTEGER(intk) :: k, j, i
        REAL(realk) :: rddxpl, rddypl, rddzpl
        REAL(realk) :: dxf, dyf, dzf
        REAL(realk) :: dudx, dudy, dudz, dvdx, dvdy, dvdz, dwdx, dwdy, dwdz

        DO i = 3, ii-2
            dxf = 0.5*dx(i-1)*rddx(i)

            DO j = 3, jj-2
                dyf = 0.5*dy(j-1)*rddy(j)

                DO k = 3, kk-2
                    dzf = 0.5*dz(k-1)*rddz(k)

                    rddxpl = 1.0/(ddx(i) + 0.5*dx(i+1) + 0.5*dx(i-1))
                    rddypl = 1.0/(ddy(j) + 0.5*dy(j+1) + 0.5*dy(j-1))
                    rddzpl = 1.0/(ddz(k) + 0.5*dz(k+1) + 0.5*dz(k-1))

                    dudx = rddx(i)*(u(k, j, i) - u(k, j, i-1))

                    dudy = rddypl * &
                        ((u(k, j+1, i) - u(k, j-1, i)) * dxf &
                        + (u(k, j+1, i-1) - u(k, j-1, i-1)) * (1.0-dxf))

                    dudz = rddzpl * &
                        ((u(k+1, j, i) - u(k-1, j, i)) * dxf &
                        + (u(k+1, j, i-1) - u(k-1, j, i-1)) * (1.0-dxf))

                    dvdx = rddxpl * &
                        ((v(k, j, i+1) - v(k, j, i-1)) * dyf &
                        + (v(k, j-1, i+1) - v(k, j-1, i-1)) * (1.0-dyf))

                    dvdy = rddy(j)*(v(k, j, i) - v(k, j-1, i))

                    dvdz = rddzpl * &
                        ((v(k+1, j, i) - v(k-1, j, i)) * dyf &
                        + (v(k+1, j-1, i) - v(k-1, j-1, i)) * (1.0-dyf))

                    dwdx = rddxpl * &
                        ((w(k, j, i+1) - w(k, j, i-1)) * dzf &
                        + (w(k-1, j, i+1) - w(k-1, j, i-1)) * (1.0-dzf))

                    dwdy = rddypl * &
                        ((w(k, j+1, i) - w(k, j-1, i)) * dzf &
                        + (w(k-1, j+1, i) - w(k-1, j-1, i)) * (1.0-dzf))

                    dwdz = rddz(k)*(w(k, j, i) - w(k-1, j, i))

                    ! masking with "bp" removes IB-intersected cells
                    dfg(k, j, i) = g(k, j, i) * bp(k, j, i) * &
                        (2.0*(dudx*dudx + dvdy*dvdy + dwdz*dwdz) &
                        + (dudy+dvdx) * (dudy+dvdx) &
                        + (dudz+dwdx) * (dudz+dwdx) &
                        + (dvdz+dwdy) * (dvdz+dwdy))
                END DO
            END DO
        END DO
    END SUBROUTINE calc_dissip


    SUBROUTINE comp_up_avg(field, name, dt)
        ! Subroutine arguments
        TYPE(field_t), INTENT(inout) :: field
        CHARACTER(len=*), INTENT(in) :: name
        REAL(realk), INTENT(in) :: dt

        ! Local variables
        INTEGER(intk), PARAMETER :: units(*) = [1, 0, -3, 0, 0, 0, 0]
        TYPE(field_t), POINTER :: in1, in2
        INTEGER(intk) :: istag, jstag, kstag

        SELECT CASE (TRIM(name))
        CASE ("UP_AVG")
            CALL get_field(in1, "U")
            CALL get_field(in2, "P")
            istag = 0
            jstag = 0
            kstag = 0
        CASE ("VP_AVG")
            CALL get_field(in1, "V")
            CALL get_field(in2, "P")
            istag = 0
            jstag = 0
            kstag = 0
        CASE ("WP_AVG")
            CALL get_field(in1, "W")
            CALL get_field(in2, "P")
            istag = 0
            jstag = 0
            kstag = 0
        CASE DEFAULT
            CALL errr(__FILE__, __LINE__)
        END SELECT

        CALL field%init(name, istag=istag, jstag=jstag, kstag=kstag, &
            units=units)
        CALL field%multiply(in1, in2)
    END SUBROUTINE comp_up_avg


    ! This subroutine computes the term du_i/dx_i * P
    SUBROUTINE comp_uxp_avg(field, name, dt)
        ! Subroutine arguments
        TYPE(field_t), INTENT(inout) :: field
        CHARACTER(len=*), INTENT(in) :: name
        REAL(realk), INTENT(in) :: dt

        ! Local variables
        TYPE(field_t) :: ux_f
        INTEGER(intk), PARAMETER :: units(*) = [1, -1, -3, 0, 0, 0, 0]
        INTEGER(intk), PARAMETER :: units_ux(*) = [0, 0, -1, 0, 0, 0, 0]
        TYPE(field_t), POINTER :: u_f, p_f
        INTEGER(intk) :: istag, jstag, kstag
        CHARACTER(len=3) :: ivar
        CHARACTER(len=2) :: name_ux

        SELECT CASE (TRIM(name))
        CASE ("UxP_AVG")
            CALL get_field(u_f, "U")
            CALL get_field(p_f, "P")
            istag = 0
            jstag = 0
            kstag = 0
            name_ux = 'UX'
            ivar = 'DDX'
        CASE ("VyP_AVG")
            CALL get_field(u_f, "V")
            CALL get_field(p_f, "P")
            istag = 0
            jstag = 0
            kstag = 0
            name_ux = 'VY'
            ivar = 'DDY'
        CASE ("WzP_AVG")
            CALL get_field(u_f, "W")
            CALL get_field(p_f, "P")
            istag = 0
            jstag = 0
            kstag = 0
            name_ux= 'WZ'
            ivar = 'DDZ'
        CASE DEFAULT
            CALL errr(__FILE__, __LINE__)
        END SELECT

        CALL field%init(name, istag=istag, jstag=jstag, kstag=kstag, &
            units=units)
        CALL ux_f%init(name_ux, istag=istag, jstag=jstag, kstag=kstag, &
            units=units_ux)

        CALL differentiate(ux_f, u_f, ivar)

        CALL field%multiply(ux_f, p_f)
        CALL ux_f%finish()
    END SUBROUTINE comp_uxp_avg


    SUBROUTINE comp_uxux_avg(field, name, dt)
        ! Subroutine arguments
        TYPE(field_t), INTENT(inout) :: field
        CHARACTER(len=*), INTENT(in) :: name
        REAL(realk), INTENT(in) :: dt

        ! Local variables
        TYPE(field_t) :: ux_f  ! It can represent any velocity component
        INTEGER(intk), PARAMETER :: units(*) = [0, 0, -2, 0, 0, 0, 0]
        INTEGER(intk), PARAMETER :: units_ux(*) = [0, 0, -1, 0, 0, 0, 0]
        TYPE(field_t), POINTER :: u_f
        INTEGER(intk) :: istag, jstag, kstag
        CHARACTER(len=3) :: ivar
        CHARACTER(len=2) :: name_ux

        SELECT CASE (TRIM(name))
        CASE ("UxUx_AVG")
            CALL get_field(u_f, "U")
            istag = 0
            jstag = 0
            kstag = 0
            name_ux = 'UX'
            ivar = 'DDX'
        CASE ("UyUy_AVG")
            CALL get_field(u_f, "U")
            istag = 1
            jstag = 1
            kstag = 0
            name_ux = 'UY'
            ivar = 'DYS'
        CASE ("UzUz_AVG")
            CALL get_field(u_f, "U")
            istag = 1
            jstag = 0
            kstag = 1
            name_ux = 'UZ'
            ivar = 'DZS'
        CASE ("VxVx_AVG")
            CALL get_field(u_f, "V")
            istag = 1
            jstag = 1
            kstag = 0
            name_ux = 'VX'
            ivar = 'DXS'
        CASE ("VyVy_AVG")
            CALL get_field(u_f, "V")
            istag = 0
            jstag = 0
            kstag = 0
            name_ux = 'VY'
            ivar = 'DDY'
        CASE ("VzVz_AVG")
            CALL get_field(u_f, "V")
            istag = 0
            jstag = 1
            kstag = 1
            name_ux = 'VZ'
            ivar = 'DZS'
        CASE ("WxWx_AVG")
            CALL get_field(u_f, "W")
            istag = 1
            jstag = 0
            kstag = 1
            name_ux = 'WX'
            ivar = 'DXS'
        CASE ("WyWy_AVG")
            CALL get_field(u_f, "W")
            istag = 0
            jstag = 1
            kstag = 1
            name_ux = 'WY'
            ivar = 'DYS'
        CASE ("WzWz_AVG")
            CALL get_field(u_f, "W")
            istag = 0
            jstag = 0
            kstag = 0
            name_ux = 'WZ'
            ivar = 'DDZ'
        CASE DEFAULT
            CALL errr(__FILE__, __LINE__)
        END SELECT

        CALL field%init(name, istag=istag, jstag=jstag, kstag=kstag, &
            units=units)
        CALL ux_f%init(name_ux, istag=istag, jstag=jstag, kstag=kstag, &
            units=units_ux)
        CALL differentiate(ux_f, u_f, ivar)

        field%arr = ux_f%arr**2
        CALL ux_f%finish()
    END SUBROUTINE comp_uxux_avg


    SUBROUTINE comp_uxvx_avg(field, name, dt)
        ! Subroutine arguments
        TYPE(field_t), INTENT(inout) :: field
        CHARACTER(len=*), INTENT(in) :: name
        REAL(realk), INTENT(in) :: dt

        ! Local Variables
        TYPE(field_t) :: ux_f, vx_f
        INTEGER(intk), PARAMETER :: units(*) = [0, 0, -2, 0, 0, 0, 0]
        INTEGER(intk), PARAMETER :: units_ux(*) = [0, 0, -1, 0, 0, 0, 0]
        ! v_f can be any velocity component different than u_f
        TYPE(field_t), POINTER :: u_f, v_f
        INTEGER(intk) :: istag, jstag, kstag
        CHARACTER(len=3) :: ivar1, ivar2
        CHARACTER(len=2) :: name_ux, name_vx

        ! Temporarily staggeredness similar to Legacy version (0, 0, 0).
        ! But it should be taken into account when creating
        ! both du_i/dk_j (Ux), and du_i/dk_j*du_k/dx_j (UxVx).

        SELECT CASE (TRIM(name))
        CASE ("UxVx_AVG")
            CALL get_field(u_f, "U")
            CALL get_field(v_f, "V")
            istag = 0
            jstag = 0
            kstag = 0
            name_ux = 'UX'
            name_vx = 'VX'
            ivar1 = 'DDX'
            ivar2 = 'DXS'
        CASE ("UxWx_AVG")
            CALL get_field(u_f, "U")
            CALL get_field(v_f, "W")
            istag = 0
            jstag = 0
            kstag = 0
            name_ux = 'UX'
            name_vx = 'WX'
            ivar1 = 'DDX'
            ivar2 = 'DXS'
        CASE ("UyVy_AVG")
            CALL get_field(u_f, "U")
            CALL get_field(v_f, "V")
            istag = 0
            jstag = 0
            kstag = 0
            name_ux = 'UY'
            name_vx = 'VY'
            ivar1 = 'DYS'
            ivar2 = 'DDY'
        CASE ("UyWy_AVG")
            CALL get_field(u_f, "U")
            CALL get_field(v_f, "W")
            istag = 0
            jstag = 0
            kstag = 0
            name_ux = 'UY'
            name_vx = 'WY'
            ivar1 = 'DYS'
            ivar2 = 'DYS'
        CASE ("UzVz_AVG")
            CALL get_field(u_f, "U")
            CALL get_field(v_f, "V")
            istag = 0
            jstag = 0
            kstag = 0
            name_ux = 'UZ'
            name_vx = 'VZ'
            ivar1 = 'DZS'
            ivar2 = 'DZS'
        CASE ("UzWz_AVG")
            CALL get_field(u_f, "U")
            CALL get_field(v_f, "W")
            istag = 0
            jstag = 0
            kstag = 0
            name_ux = 'UZ'
            name_vx = 'WZ'
            ivar1 = 'DZS'
            ivar2 = 'DDZ'
        CASE ("VxWx_AVG")
            CALL get_field(u_f, "V")
            CALL get_field(v_f, "W")
            istag = 0
            jstag = 0
            kstag = 0
            name_ux = 'VX'
            name_vx = 'WX'
            ivar1 = 'DXS'
            ivar2 = 'DXS'
        CASE ("VyWy_AVG")
            CALL get_field(u_f, "V")
            CALL get_field(v_f, "W")
            istag = 0
            jstag = 0
            kstag = 0
            name_ux = 'VY'
            name_vx = 'WY'
            ivar1 = 'DDY'
            ivar2 = 'DYS'
        CASE ("VzWz_AVG")
            CALL get_field(u_f, "V")
            CALL get_field(v_f, "W")
            istag = 0
            jstag = 0
            kstag = 0
            name_ux = 'VZ'
            name_vx = 'WZ'
            ivar1 = 'DZS'
            ivar2 = 'DDZ'
        CASE DEFAULT
            CALL errr(__FILE__, __LINE__)
        END SELECT

        CALL field%init(name, istag=istag, jstag=jstag, kstag=kstag, &
            units=units)
        CALL ux_f%init(name_ux, istag=istag, jstag=jstag, kstag=kstag, &
            units=units_ux)
        CALL vx_f%init(name_vx, istag=istag, jstag=jstag, kstag=kstag, &
            units=units_ux)

        CALL differentiate(ux_f, u_f, ivar1)
        CALL differentiate(vx_f, v_f, ivar2)

        ! Since temporarily ux and vx are located at the same point
        ! multiply method won't be used.
        ! TO DO: if location is not the same multiply method should be used.

        field%arr = ux_f%arr * vx_f%arr
        CALL ux_f%finish()
        CALL vx_f%finish()
    END SUBROUTINE comp_uxvx_avg


    SUBROUTINE comp_uyvxp_avg(field, name, dt)
        ! Subroutine arguments
        TYPE(field_t), INTENT(inout) :: field
        CHARACTER(len=*), INTENT(in) :: name
        REAL(realk), INTENT(in) :: dt

        ! Local Variables
        TYPE(field_t) :: uy_f, vx_f, temp_f
        INTEGER(intk), PARAMETER :: units(*) = [0, 0, -2, 0, 0, 0, 0]
        INTEGER(intk), PARAMETER :: units_ux(*) = [0, 0, -1, 0, 0, 0, 0]
        TYPE(field_t), POINTER :: u_f, v_f, p_f
        INTEGER(intk) :: istag, jstag, kstag
        CHARACTER(len=3) :: ivar1, ivar2
        CHARACTER(len=2) :: name_uy, name_vx

        SELECT CASE (TRIM(name))
        CASE ("UyP+VxP_AVG")
            CALL get_field(u_f, "U")
            CALL get_field(v_f, "V")
            istag = 1
            jstag = 1
            kstag = 0
            name_uy = 'UY'
            name_vx = 'VX'
            ivar1 = 'DYS'
            ivar2 = 'DXS'
        CASE("UzP+WxP_AVG")
            CALL get_field(u_f, "U")
            CALL get_field(v_f, "W")
            istag = 1
            jstag = 0
            kstag = 1
            name_uy = 'UZ'
            name_vx = 'WX'
            ivar1 = 'DZS'
            ivar2 = 'DXS'
        CASE ("VzP+WyP_AVG")
            CALL get_field(u_f, "V")
            CALL get_field(v_f, "W")
            istag = 0
            jstag = 1
            kstag = 1
            name_uy = 'VZ'
            name_vx = 'WY'
            ivar1 = 'DZS'
            ivar2 = 'DYS'
        CASE DEFAULT
            CALL errr(__FILE__, __LINE__)
        END SELECT

        CALL field%init(name, istag=istag, jstag=jstag, kstag=kstag, &
            units=units)
        CALL uy_f%init(name_uy, istag=istag, jstag=jstag, kstag=kstag, &
            units=units_ux)
        CALL vx_f%init(name_vx, istag=istag, jstag=jstag, kstag=kstag, &
            units=units_ux)
        CALL temp_f%init('tmp', istag=istag, jstag=jstag, kstag=kstag, &
            units=units_ux)

        CALL differentiate(uy_f, u_f, ivar1)
        CALL differentiate(vx_f, v_f, ivar2)
        temp_f%arr = uy_f%arr + vx_f%arr

        CALL get_field(p_f, "P")
        CALL field%multiply(temp_f, p_f)

        CALL uy_f%finish()
        CALL vx_f%finish()
        CALL temp_f%finish()
    END SUBROUTINE comp_uyvxp_avg

END MODULE flowstat_mod
