MODULE timeintegration_mod
    USE bound_flow_mod
    USE core_mod
    USE flowcore_mod
    USE gc_flowstencils_mod, ONLY: setpointvalues, setibvalues, getibvalues
    USE ib_mod
    USE itinfo_mod, ONLY: itinfo_sample, itinfo_print
    USE lesmodel_mod
    USE pressuresolver_mod
    USE tstle4_mod
    USE setboundarybuffers_mod
    USE boussinesqterm_mod, ONLY: boussinesqterm

    IMPLICIT NONE(type, external)
    PRIVATE

    PUBLIC :: timeintegrate_flow, itinfo_flow

CONTAINS
    SUBROUTINE timeintegrate_flow(itstep, ittot, timeph, dt, irk, rkscheme)
        ! Subroutine arguments
        INTEGER(intk), INTENT(in) :: itstep
        INTEGER(intk), INTENT(in) :: ittot
        REAL(realk), INTENT(in) :: timeph
        REAL(realk), INTENT(in) :: dt
        INTEGER(intk), INTENT(in) :: irk
        TYPE(rk_2n_t), INTENT(in) :: rkscheme

        ! Local variables
        LOGICAL :: lastrk
        INTEGER(intk) :: ilevel
        REAL(realk) :: frhs, fu, dtrk, dtrki, timerk
        TYPE(field_t), POINTER :: u, v, w, ut, vt, wt, pwu, pwv, pww, p, g
        TYPE(field_t), POINTER :: du, dv, dw
        TYPE(field_t) :: uo, vo, wo

        ! Just return if no flow is to be solved
        IF (.NOT. solve_flow) RETURN
        CALL start_timer(300)

        CALL get_field(u, "U")
        CALL get_field(v, "V")
        CALL get_field(w, "W")
        CALL get_field(p, "P")
        CALL get_field(g, "G")

        ! In all implemented RK schemes FRHS is 0.0 for IRK 1, this means
        ! that the method itself takes care of "initializing" these fields
        ! to zero, and in case we implement schemes with non-zero FRHS for
        ! the first step in the future, we should not zeroize them either,
        ! becuase then the method is not self-starting... So absolutely do not
        ! set these to zero here!
        CALL get_field(du, "DU")
        CALL get_field(dv, "DV")
        CALL get_field(dw, "DW")

        CALL uo%init("UO")
        CALL vo%init("VO")
        CALL wo%init("WO")

        ! Transporting velocities for the convective terms
        ! Only CC use a different transporting velocity
        CALL get_field(ut, "U")
        CALL get_field(vt, "V")
        CALL get_field(wt, "W")

        ! The transported velocities are different between GC and NONE ib types
        SELECT CASE (ib%type)
        CASE ("GHOSTCELL")
            ! "PWU" (PunktWerte U) etc. for GC
            CALL get_field(pwu, "PWU")
            CALL get_field(pwv, "PWV")
            CALL get_field(pww, "PWW")
        CASE ("NONE")
            ! No IB use u, w, w for the transported quantities
            CALL get_field(pwu, "U")
            CALL get_field(pwv, "V")
            CALL get_field(pww, "W")
        CASE DEFAULT
            CALL errr(__FILE__, __LINE__)
        END SELECT

        ! Compute the GC fluxes
        IF (irk == 1 .AND. ib%type == "GHOSTCELL") THEN
            CALL setibvalues(u, v, w)
        END IF

        ! TSTLE4 zeroize uo, vo, wo before use internally
        CALL tstle4(uo, vo, wo, pwu, pwv, pww, ut, vt, wt, p, g)
        CALL boussinesqterm(uo, vo, wo)

        ! In IRK 1, FRHS is zero, therefore we do not need to zeroize
        ! the du, dv, dw fields before each step
        CALL rkscheme%get_coeffs(frhs, fu, dtrk, dtrki, irk)

        ! dU_j = A_j*dU_(j-1) + dt*uo
        du%arr = frhs*du%arr + uo%arr
        dv%arr = frhs*dv%arr + vo%arr
        dw%arr = frhs*dw%arr + wo%arr

        ! U_j = U_(j-1) + B_j*dU_j
        u%arr = u%arr + dt*fu*du%arr
        v%arr = v%arr + dt*fu*dv%arr
        w%arr = w%arr + dt*fu*dw%arr

        IF (ib%type == "GHOSTCELL") THEN
            ! Equivalent to old "cop3dzero"
            CALL maskbp(u, v, w, p)

            ! Equivalent to old boundmg with ityp 'R'
            CALL getibvalues(u, v, w)
        END IF

        IF (uinf_is_time) THEN
            DO ilevel = minlevel, maxlevel
                timerk = timeph + dt*dtrk
                CALL setboundarybuffers%bound(ilevel, u, v, w, timeph=timerk)
            END DO
        END IF

        ! For divergence computation
        DO ilevel = minlevel, maxlevel
            CALL connect(ilevel, 1, v1=u, v2=v, v3=w, &
                normal=.true., forward=1)
            CALL parent(ilevel, u, v, w, p)
            CALL bound_flow%bound(ilevel, u, v, w, p)
        END DO

        ! TODO: check dtrk
        CALL mgpoisl(u, v, w, p, dtrk*dt, ittot, irk)
        CALL lesmodel(g)

        IF (ib%type == "GHOSTCELL") THEN
            lastrk = (irk == rkscheme%nrk)
            CALL setpointvalues(pwu, pwv, pww, u, v, w, lastrk)
        END IF

        ! TODO: mgplevel

        CALL stop_timer(300)
    END SUBROUTINE timeintegrate_flow


    SUBROUTINE itinfo_flow(itstep, ittot, timeph, dt, globalcflmax, exploded)
        ! Subroutine arguments
        INTEGER(intk), INTENT(in) :: itstep
        INTEGER(intk), INTENT(in) :: ittot
        REAL(realk), INTENT(in) :: timeph
        REAL(realk), INTENT(in) :: dt
        REAL(realk), INTENT(out) :: globalcflmax
        INTEGER(intk), INTENT(inout) :: exploded

        ! Local variables
        INTEGER(intk) :: i, igrid
        INTEGER(intk) :: kk, jj, ii
        TYPE(field_t), POINTER :: u_f, v_f, w_f, g_f, bp_f, sdiv_f
        TYPE(field_t), POINTER :: x_f, y_f, z_f
        TYPE(field_t), POINTER :: dx_f, dy_f, dz_f
        TYPE(field_t), POINTER :: ddx_f, ddy_f, ddz_f
        TYPE(field_t), POINTER :: rddx_f, rddy_f, rddz_f

        REAL(realk), POINTER, CONTIGUOUS :: u(:, :, :), v(:, :, :), &
            w(:, :, :), g(:, :, :), bp(:, :, :), sdiv(:, :, :)
        REAL(realk), POINTER, CONTIGUOUS :: x(:), y(:), z(:)
        REAL(realk), POINTER, CONTIGUOUS :: dx(:), dy(:), dz(:)
        REAL(realk), POINTER, CONTIGUOUS :: ddx(:), ddy(:), ddz(:)
        REAL(realk), POINTER, CONTIGUOUS :: rddx(:), rddy(:), rddz(:)

        REAL(realk) :: cflmax, cflmax_pos(3)
        REAL(realk) :: divmax, divmax_pos(3)
        REAL(realk) :: esumg, esums

        IF (.NOT. solve_flow) RETURN
        CALL start_timer(300)
        CALL start_timer(350)

        ! Get fields
        CALL get_field(u_f, "U")
        CALL get_field(v_f, "V")
        CALL get_field(w_f, "W")
        CALL get_field(g_f, "G")
        CALL get_field(bp_f, "BP")
        CALL get_field(sdiv_f, "SDIV")

        CALL get_field(x_f, "X")
        CALL get_field(y_f, "Y")
        CALL get_field(z_f, "Z")

        CALL get_field(dx_f, "DX")
        CALL get_field(dy_f, "DY")
        CALL get_field(dz_f, "DZ")

        CALL get_field(ddx_f, "DDX")
        CALL get_field(ddy_f, "DDY")
        CALL get_field(ddz_f, "DDZ")

        CALL get_field(rddx_f, "RDDX")
        CALL get_field(rddy_f, "RDDY")
        CALL get_field(rddz_f, "RDDZ")

        ! Compute CFL, divergence and kinetic energy
        DO i = 1, nmygrids
            igrid = mygrids(i)

            CALL get_mgdims(kk, jj, ii, igrid)

            CALL u_f%get_ptr(u, igrid)
            CALL v_f%get_ptr(v, igrid)
            CALL w_f%get_ptr(w, igrid)
            CALL g_f%get_ptr(g, igrid)
            CALL bp_f%get_ptr(bp, igrid)
            CALL sdiv_f%get_ptr(sdiv, igrid)

            CALL x_f%get_ptr(x, igrid)
            CALL y_f%get_ptr(y, igrid)
            CALL z_f%get_ptr(z, igrid)

            CALL dx_f%get_ptr(dx, igrid)
            CALL dy_f%get_ptr(dy, igrid)
            CALL dz_f%get_ptr(dz, igrid)

            CALL ddx_f%get_ptr(ddx, igrid)
            CALL ddy_f%get_ptr(ddy, igrid)
            CALL ddz_f%get_ptr(ddz, igrid)

            CALL rddx_f%get_ptr(rddx, igrid)
            CALL rddy_f%get_ptr(rddy, igrid)
            CALL rddz_f%get_ptr(rddz, igrid)

            CALL compcflmax(cflmax, cflmax_pos, dt, kk, jj, ii, &
                u, v, w, bp, x, y, z, dx, dy, dz, ddx, ddy, ddz)
            CALL compdivmax(divmax, divmax_pos, kk, jj, ii, &
                u, v, w, x, y, z, rddx, rddy, rddz, bp, sdiv)
            CALL enerfg(esumg, kk, jj, ii, u, v, w, dx, dy, dz, ddx, ddy, ddz)
            CALL enerfs(esums, kk, jj, ii, g, ddx, ddy, ddz)

            CALL itinfo_sample(igrid, divmax, divmax_pos, cflmax, cflmax_pos, &
                esumg, esums)
        END DO

        CALL itinfo_print(itstep, ittot, timeph, globalcflmax, exploded)

        CALL stop_timer(350)
        CALL stop_timer(300)
    END SUBROUTINE itinfo_flow


    SUBROUTINE maskbp(u, v, w, p)
        ! Subroutine arguments
        TYPE(field_t), INTENT(inout) :: u, v, w, p

        ! Local variables
        INTEGER(intk) :: i, ip3, igrid
        INTEGER(intk) :: kk, jj, ii
        TYPE(field_t), POINTER :: bp

        CALL get_field(bp, "BP")

        DO i = 1, nmygrids
            igrid = mygrids(i)
            CALL get_mgdims(kk, jj, ii, igrid)
            CALL get_ip3(ip3, igrid)

            CALL maskbp_grid(kk, jj, ii, u%arr(ip3), v%arr(ip3), w%arr(ip3), &
                p%arr(ip3), bp%arr(ip3))
        END DO
    END SUBROUTINE maskbp


    PURE SUBROUTINE maskbp_grid(kk, jj, ii, u, v, w, p, bp)
        ! Subroutine arguments
        INTEGER(intk), INTENT(in) :: kk, jj, ii
        REAL(realk), INTENT(inout) :: u(kk, jj, ii), v(kk, jj, ii), &
            w(kk, jj, ii), p(kk, jj, ii)
        REAL(realk), INTENT(in) :: bp(kk, jj, ii)

        ! Local variables
        INTEGER :: k, j, i

        DO i = 1, ii-1
            DO j = 1, jj-1
                DO k = 1, kk-1
                    u(k, j, i) = u(k, j, i)*bp(k, j, i)*bp(k, j, i+1)
                END DO
                DO k = 1, kk-1
                    v(k, j, i) = v(k, j, i)*bp(k, j, i)*bp(k, j+1, i)
                END DO
                DO k = 1, kk-1
                    w(k, j, i) = w(k, j, i)*bp(k, j, i)*bp(k+1, j, i)
                END DO
                DO k = 1, kk-1
                    p(k, j, i) = p(k, j, i)*bp(k, j, i)
                END DO
            END DO
        END DO
    END SUBROUTINE maskbp_grid


    PURE SUBROUTINE compcflmax(cflmax, cflmax_pos, dt, kk, jj, ii, &
            u, v, w, bp, x, y, z, dx, dy, dz, ddx, ddy, ddz)
        ! Subroutine arguments
        REAL(realk), INTENT(out) :: cflmax, cflmax_pos(3)
        REAL(realk), INTENT(in) :: dt
        INTEGER(intk), INTENT(in) :: kk, jj, ii
        REAL(realk), INTENT(in) :: u(kk, jj, ii), v(kk, jj, ii), w(kk, jj, ii)
        REAL(realk), INTENT(in) :: bp(kk, jj, ii)
        REAL(realk), INTENT(in) :: x(ii), y(jj), z(kk)
        REAL(realk), INTENT(in) :: dx(ii), dy(jj), dz(kk)
        REAL(realk), INTENT(in) :: ddx(ii), ddy(jj), ddz(kk)

        ! Local variables
        REAL(realk) :: cflmaxtemp, cflu, cflv, cflw
        INTEGER(intk) :: k, j, i

        cflmaxtemp = 0.0
        cflmax = 0.0
        cflmax_pos = 0.0

        DO i = 3, ii-2
            DO j = 3, jj-2
                DO k = 3, kk-2
                    cflu = ABS(2.0*u(k, j, i-1) + (u(k, j, i) &
                            - u(k, j, i-1))*dx(i-1)/ddx(i))/ddx(i) &
                        + ABS(2.0*u(k, j, i) + (u(k, j, i+1) &
                            - u(k, j, i))*dx(i)/ddx(i+1))/ddx(i) &
                        + ABS(v(k, j, i) + v(k, j, i+1))/ddy(j) &
                        + ABS(v(k, j-1, i) + v(k, j-1, i+1))/ddy(j) &
                        + ABS(w(k, j, i) + w(k, j, i+1))/ddz(k) &
                        + ABS(w(k-1, j, i) + w(k-1, j, i+1))/ddz(k)
                    cflu = cflu*bp(k, j, i)*bp(k, j, i+1)

                    cflv = ABS(u(k, j, i ) + u(k, j+1, i))/ddx(i) &
                        + ABS(u(k, j, i-1) + u(k, j+1, i-1))/ddx(i) &
                        + ABS(2.0*v(k, j-1, i) + (v(k, j, i) &
                            - v(k, j-1, i))*dy(j-1)/ddy(j))/ddy(j) &
                        + ABS(2.0*v(k, j, i) + (v(k, j+1, i) &
                            - v(k, j, i))*dy(j)/ddy(j+1))/ddy(j) &
                        + ABS(w(k, j, i) + w(k, j+1, i))/ddz(k) &
                        + ABS(w(k-1, j, i) + w(k-1, j+1, i))/ddz(k)
                    cflv = cflv*bp(k, j, i)*bp(k, j+1, i)

                    cflw = ABS(u(k, j, i) + u(k+1, j, i))/ddx(i) &
                        + ABS(u(k, j, i-1) + u(k+1, j, i-1))/ddx(i) &
                        + ABS(v(k, j, i) + v(k+1, j, i  ))/ddy(j) &
                        + ABS(v(k, j-1, i) + v(k+1, j-1, i))/ddy(j) &
                        + ABS(2.0*w(k-1, j, i) + (w(k, j, i) &
                            - w(k-1, j, i))*dz(k-1)/ddz(k))/ddz(k) &
                        + ABS(2.0*w(k, j, i) + (w(k+1, j, i) &
                            - w(k, j, i))*dz(k)/ddz(k+1))/ddz(k)
                    cflw = cflw*bp(k, j, i)*bp(k+1, j, i)

                    IF (MAX(cflu, cflv, cflw) > cflmaxtemp) THEN
                        cflmaxtemp = MAX(cflu, cflv, cflw)
                        cflmax_pos(1) = x(i)
                        cflmax_pos(2) = y(j)
                        cflmax_pos(3) = z(k)
                    END IF
                END DO
            END DO
        END DO

        cflmax = cflmaxtemp * 0.25 * dt
    END SUBROUTINE compcflmax


    ! Simpler routine than "divcal", because this one does only compute the
    ! /max/ divergence and its location, and does not explicitly fill in
    ! a full 3-D divergcence field.
    PURE SUBROUTINE compdivmax(divmax, divmax_pos, kk, jj, ii, u, v, w, &
            x, y, z, rddx, rddy, rddz, bp, sdiv)
        ! Subroutine arguments
        REAL(realk), INTENT(out) :: divmax, divmax_pos(3)
        INTEGER(intk), INTENT(in) :: kk, jj, ii
        REAL(realk), INTENT(in) :: u(kk, jj, ii)
        REAL(realk), INTENT(in) :: v(kk, jj, ii)
        REAL(realk), INTENT(in) :: w(kk, jj, ii)
        REAL(realk), INTENT(in) :: x(ii), y(jj), z(kk)
        REAL(realk), INTENT(in) :: rddx(ii), rddy(jj), rddz(kk)
        REAL(realk), INTENT(in) :: bp(kk, jj, ii)
        REAL(realk), INTENT(in) :: sdiv(kk, jj, ii)

        ! Local variables
        INTEGER(intk) :: k, j, i
        REAL(realk) :: div

        divmax = 0.0
        divmax_pos = 0.0

        DO i = 3, ii-2
            DO j = 3, jj-2
                DO k = 3, kk-2
                    div = sdiv(k, j, i) &
                        + (u(k, j, i) - u(k, j, i-1))*rddx(i) &
                        + (v(k, j, i) - v(k, j-1, i))*rddy(j) &
                        + (w(k, j, i) - w(k-1, j, i))*rddz(k)
                    ! div = bp(k, j, i)*div

                    IF (ABS(div) > divmax) THEN
                        divmax = ABS(div)
                        divmax_pos(1) = x(i)
                        divmax_pos(2) = y(j)
                        divmax_pos(3) = z(k)
                    END IF
                END DO
            END DO
        END DO
    END SUBROUTINE compdivmax


    SUBROUTINE enerfg(esum, kk, jj, ii, u, v, w, dx, dy, dz, ddx, ddy, ddz)
        ! Subroutine arguments
        REAL(realk), INTENT(out) :: esum
        INTEGER(intk), INTENT(in) :: kk, jj, ii
        REAL(realk), INTENT(in) :: u(kk, jj, ii), v(kk, jj, ii), w(kk, jj, ii)
        REAL(realk), INTENT(in) :: dx(ii), dy(jj), dz(kk)
        REAL(realk), INTENT(in) :: ddx(ii), ddy(jj), ddz(kk)

        ! Local variables
        REAL(realk) :: up, vp, wp, vsum, vol
        INTEGER(intk) :: k, j, i

        esum = 0.0
        vsum = 0.0

        DO i = 3, ii-2
            DO j = 3, jj-2
                DO k = 3, kk-2
                    up = u(k, j, i-1) &
                        + (u(k, j, i) - u(k, j, i-1))*0.5*dx(i-1)/ddx(i)

                    vp = v(k, j-1, i) &
                        + (v(k, j, i) - v(k, j-1, i))*0.5*dy(j-1)/ddy(j)

                    wp = (w(k-1, j, i) &
                        + (w(k, j, i) - w(k-1, j, i))*0.5*dz(k-1)/ddz(k))

                    vol = ddx(i)*ddy(j)*ddz(k)
                    esum = esum + vol*0.5*(up**2 + vp**2 + wp**2)
                    vsum = vsum + vol
                END DO
            END DO
        END DO

        esum = esum/vol
    END SUBROUTINE enerfg


    SUBROUTINE enerfs(esum, kk, jj, ii, g, ddx, ddy, ddz)
        ! Subroutine arguments
        REAL(realk), INTENT(out) :: esum
        INTEGER(intk), INTENT(in) :: kk, jj, ii
        REAL(realk), INTENT(in) :: g(kk, jj, ii)
        REAL(realk), INTENT(in) :: ddx(ii), ddy(jj), ddz(kk)

        ! Local variables
        REAL(realk), PARAMETER :: conv2s = 0.094
        REAL(realk) :: delta, ener, vsum, vol
        INTEGER(intk) :: k, j, i

        esum = 0.0
        vsum = 0.0

        DO i = 3, ii-2
            DO j = 3, jj-2
                DO k = 3, kk-2
                    delta = cube_root(ddx(i)*ddy(j)*ddz(k))
                    ener = ((g(k, j, i) - gmol)/(rho*conv2s*delta))**2

                    vol = ddx(i)*ddy(j)*ddz(k)
                    esum = esum + vol*ener
                    vsum = vsum + vol
                END DO
            END DO
        END DO

        esum = esum/vol
    END SUBROUTINE enerfs
END MODULE timeintegration_mod
