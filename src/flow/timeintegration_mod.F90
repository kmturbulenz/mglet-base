MODULE timeintegration_mod
    USE bound_flow_mod
    USE core_mod
    USE flowcore_mod
    USE gc_compbodyforce_mod, ONLY: sample_compbodyforce
    USE gc_flowstencils_mod, ONLY: setpointvalues, setibvalues, getibvalues
    USE ib_mod
    USE itinfo_mod, ONLY: itinfo_sample, itinfo_print
    USE lesmodel_mod
    USE pressuresolver_mod
    USE tstle4_mod
    USE setboundarybuffers_mod
    USE boussinesqterm_mod, ONLY: boussinesqterm
    USE coriolisterm_mod, ONLY: coriolisterm

    IMPLICIT NONE(type, external)
    PRIVATE

    ! Custom reductions to bundle a value with its cell indices on the grid.
    ! Surprisingly, using derived types here does not come with runtime cost.
    TYPE :: valpos_t
        REAL(realk) :: val
        REAL(realk) :: x, y, z
    END TYPE valpos_t
    !$omp declare reduction(valpos: valpos_t : &
    !$omp& omp_out = get_maxvalpos(omp_out, omp_in)) &
    !$omp& initializer(valpos_init(omp_priv))

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
        TYPE(field_t), POINTER :: du, dv, dw, uo, vo, wo

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

        CALL push_field(uo, "UO")
        CALL push_field(vo, "VO")
        CALL push_field(wo, "WO")
        CALL set_field_arr(uo, 0.0_realk)
        CALL set_field_arr(vo, 0.0_realk)
        CALL set_field_arr(wo, 0.0_realk)

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

        ! In IRK 1, FRHS is zero, therefore we do not need to zeroize
        ! the du, dv, dw fields before each step
        CALL rkscheme%get_coeffs(frhs, fu, dtrk, dtrki, irk)

        ! Compute the GC fluxes
        IF (irk == 1 .AND. ib%type == "GHOSTCELL") THEN
            CALL setibvalues(u, v, w)
        END IF

        ! TSTLE4 zeroize uo, vo, wo before use internally
        CALL tstle4(uo, vo, wo, pwu, pwv, pww, ut, vt, wt, p, g)
        CALL boussinesqterm(uo, vo, wo)
        CALL coriolisterm(uo, vo, wo)

        ! dU_j = A_j*dU_(j-1) + dt*uo
        ! U_j = U_(j-1) + B_j*dU_j
        CALL rkstep(u%arr, du%arr, uo%arr, frhs, dt*fu)
        CALL rkstep(v%arr, dv%arr, vo%arr, frhs, dt*fu)
        CALL rkstep(w%arr, dw%arr, wo%arr, frhs, dt*fu)

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
            CALL connect(ilevel, 1, v1=u, v2=v, v3=w, normal=.true., forward=1)
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

        CALL pop_field(wo)
        CALL pop_field(vo)
        CALL pop_field(uo)

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
        TYPE(field_t), POINTER :: u, v, w, p, g, bp, sdiv
        TYPE(field_t), POINTER :: x, y, z
        TYPE(field_t), POINTER :: dx, dy, dz
        TYPE(field_t), POINTER :: ddx, ddy, ddz
        TYPE(field_t), POINTER :: rddx, rddy, rddz

        IF (.NOT. solve_flow) RETURN
        CALL start_timer(300)
        CALL start_timer(350)

        ! Get fields
        CALL get_field(u, "U")
        CALL get_field(v, "V")
        CALL get_field(w, "W")
        CALL get_field(g, "G")
        CALL get_field(bp, "BP")
        CALL get_field(sdiv, "SDIV")

        CALL get_field(x, "X")
        CALL get_field(y, "Y")
        CALL get_field(z, "Z")

        CALL get_field(dx, "DX")
        CALL get_field(dy, "DY")
        CALL get_field(dz, "DZ")

        CALL get_field(ddx, "DDX")
        CALL get_field(ddy, "DDY")
        CALL get_field(ddz, "DDZ")

        CALL get_field(rddx, "RDDX")
        CALL get_field(rddy, "RDDY")
        CALL get_field(rddz, "RDDZ")

        ! TODO(offload): Remove once surrounding subroutines are offloaded
        CALL map_arr_to_device(u, v, w, g)

        CALL compcflmax(dt, u, v, w, bp, x, y, z, dx, dy, dz, ddx, ddy, ddz)
        CALL compdivmax(u, v, w, bp, x, y, z, rddx, rddy, rddz, sdiv)
        CALL compenergy(u, v, w, g, dx, dy, dz, ddx, ddy, ddz)

        CALL itinfo_print(itstep, ittot, timeph, globalcflmax, exploded)

        CALL stop_timer(350)

        IF (compbodyforce .AND. ib%type == "GHOSTCELL") THEN
            CALL get_field(p, "P")
            CALL sample_compbodyforce(u, v, w, p, g, ittot, timeph)
        END IF

        CALL stop_timer(300)
    END SUBROUTINE itinfo_flow


    SUBROUTINE compcflmax(dt, u_f, v_f, w_f, bp_f, x_f, y_f, z_f, &
            dx_f, dy_f, dz_f, ddx_f, ddy_f, ddz_f)
        ! Subroutine arguments
        REAL(realk), INTENT(in) :: dt
        TYPE(field_t), INTENT(in) :: u_f, v_f, w_f
        TYPE(field_t), INTENT(in) :: bp_f
        TYPE(field_t), INTENT(in) :: x_f, y_f, z_f
        TYPE(field_t), INTENT(in) :: dx_f, dy_f, dz_f
        TYPE(field_t), INTENT(in) :: ddx_f, ddy_f, ddz_f

        ! Local variables
        REAL(realk), ALLOCATABLE :: cflmax_grid(:)
        REAL(realk), ALLOCATABLE :: cflmax_pos_grid(:, :)
        INTEGER(intk) :: i, igrid

        ALLOCATE(cflmax_grid(nmygrids))
        ALLOCATE(cflmax_pos_grid(3, nmygrids))

#ifdef _MGLET_PROFILE_ANNOTATIONS_
        CALL profile_range_push("compcflmax")
#endif

        CALL compcflmax_impl(u_f%arr, v_f%arr, w_f%arr, bp_f%arr, &
            x_f%arr, y_f%arr, z_f%arr, dx_f%arr, dy_f%arr, dz_f%arr, &
            ddx_f%arr, ddy_f%arr, ddz_f%arr, cflmax_grid, &
            cflmax_pos_grid, dt)

#ifdef _MGLET_PROFILE_ANNOTATIONS_
        CALL profile_range_pop()
#endif

        DO i = 1, nmygrids
            igrid = mygrids(i)
            CALL itinfo_sample(igrid, cflmax=cflmax_grid(i), &
                cflmax_pos=cflmax_pos_grid(:, i))
        END DO

        DEALLOCATE(cflmax_grid)
        DEALLOCATE(cflmax_pos_grid)
    END SUBROUTINE compcflmax


    SUBROUTINE compcflmax_impl(u, v, w, bp, x, y, z, dx, dy, dz, ddx, ddy, &
            ddz, cflmaxgrid, cflmaxposgrid, dt)
        ! Subroutine arguments
        REAL(realk), INTENT(in) :: u(*), v(*), w(*), bp(*)
        REAL(realk), INTENT(in) :: x(*), y(*), z(*)
        REAL(realk), INTENT(in) :: dx(*), dy(*), dz(*)
        REAL(realk), INTENT(in) :: ddx(*), ddy(*), ddz(*)
        REAL(realk), INTENT(inout) :: cflmaxgrid(*)
        REAL(realk), INTENT(inout) :: cflmaxposgrid(3, *)
        REAL(realk), INTENT(in) :: dt

        ! Local variables
        INTEGER(intk) :: i, igrid, kk, jj, ii, ip3, ip1x, ip1y, ip1z
        TYPE(valpos_t) :: maxvalpos

        !$omp target teams distribute private(i, igrid, kk, jj, ii, ip3, ip1x, &
        !$omp& ip1y, ip1z, maxvalpos) &
        !$omp& map(from: cflmaxgrid(1:nmygrids), &
        !$omp& cflmaxposgrid(1:3, 1:nmygrids))
        DO i = 1, nmygrids
            igrid = mygrids(i)
            cflmaxgrid(i) = -HUGE(1.0_realk)
            cflmaxposgrid(1, i) = -HUGE(1.0_realk)
            cflmaxposgrid(2, i) = -HUGE(1.0_realk)
            cflmaxposgrid(3, i) = -HUGE(1.0_realk)

            CALL get_mgdims(kk, jj, ii, igrid)
            CALL get_ip3(ip3, igrid)
            CALL get_ip1x(ip1x, igrid)
            CALL get_ip1y(ip1y, igrid)
            CALL get_ip1z(ip1z, igrid)

            CALL valpos_init(maxvalpos)
            !$omp parallel
            CALL compcflmax_grid(maxvalpos, dt, &
                kk, jj, ii, u(ip3), v(ip3), w(ip3), bp(ip3), &
                x(ip1x), y(ip1y), z(ip1z), dx(ip1x), dy(ip1y), dz(ip1z), &
                ddx(ip1x), ddy(ip1y), ddz(ip1z))
            !$omp end parallel

            IF (maxvalpos%val >= 0.0_realk) THEN
                cflmaxgrid(i) = maxvalpos%val * 0.25 * dt
                cflmaxposgrid(1, i) = maxvalpos%x
                cflmaxposgrid(2, i) = maxvalpos%y
                cflmaxposgrid(3, i) = maxvalpos%z
            END IF
        END DO
        !$omp end target teams distribute
    END SUBROUTINE compcflmax_impl


    SUBROUTINE compcflmax_grid(maxvalpos, dt, kk, jj, ii, &
            u, v, w, bp, x, y, z, dx, dy, dz, ddx, ddy, ddz)
        !$omp declare target
        ! Subroutine arguments
        TYPE(valpos_t), INTENT(inout) :: maxvalpos
        REAL(realk), INTENT(in) :: dt
        INTEGER(intk), INTENT(in) :: kk, jj, ii
        REAL(realk), INTENT(in) :: u(kk, jj, ii), v(kk, jj, ii), w(kk, jj, ii)
        REAL(realk), INTENT(in) :: bp(kk, jj, ii)
        REAL(realk), INTENT(in) :: x(ii), y(jj), z(kk)
        REAL(realk), INTENT(in) :: dx(ii), dy(jj), dz(kk)
        REAL(realk), INTENT(in) :: ddx(ii), ddy(jj), ddz(kk)

        ! Local variables
        TYPE(valpos_t) :: localvalpos
        REAL(realk) :: cflmaxtemp, cflu, cflv, cflw
        INTEGER(intk) :: k, j, i

        cflmaxtemp = -HUGE(1.0_realk)

        !$omp do collapse(3) private(i, j, k, cflu, cflv, cflw, &
        !$omp& cflmaxtemp, localvalpos) &
        !$omp& reduction(valpos:maxvalpos)
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

                    cflv = ABS(u(k, j, i) + u(k, j+1, i))/ddx(i) &
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
                        + ABS(v(k, j, i) + v(k+1, j, i))/ddy(j) &
                        + ABS(v(k, j-1, i) + v(k+1, j-1, i))/ddy(j) &
                        + ABS(2.0*w(k-1, j, i) + (w(k, j, i) &
                            - w(k-1, j, i))*dz(k-1)/ddz(k))/ddz(k) &
                        + ABS(2.0*w(k, j, i) + (w(k+1, j, i) &
                            - w(k, j, i))*dz(k)/ddz(k+1))/ddz(k)
                    cflw = cflw*bp(k, j, i)*bp(k+1, j, i)

                    cflmaxtemp = MAX(cflu, cflv, cflw)

                    localvalpos%val = cflmaxtemp
                    localvalpos%x = x(i)
                    localvalpos%y = y(j)
                    localvalpos%z = z(k)
                    maxvalpos = get_maxvalpos(maxvalpos, localvalpos)
                END DO
            END DO
        END DO
        !$omp end do
    END SUBROUTINE compcflmax_grid


    SUBROUTINE compdivmax(u_f, v_f, w_f, bp_f, x_f, y_f, z_f, &
            rddx_f, rddy_f, rddz_f, sdiv_f)
        ! Subroutine arguments
        TYPE(field_t), INTENT(in) :: u_f, v_f, w_f
        TYPE(field_t), INTENT(in) :: bp_f
        TYPE(field_t), INTENT(in) :: x_f, y_f, z_f
        TYPE(field_t), INTENT(in) :: rddx_f, rddy_f, rddz_f
        TYPE(field_t), INTENT(in) :: sdiv_f

        ! Local variables
        REAL(realk), ALLOCATABLE :: divmax_grid(:)
        REAL(realk), ALLOCATABLE :: divmax_pos_grid(:, :)
        INTEGER(intk) :: i, igrid

        ALLOCATE(divmax_grid(nmygrids))
        ALLOCATE(divmax_pos_grid(3, nmygrids))

#ifdef _MGLET_PROFILE_ANNOTATIONS_
        CALL profile_range_push("compdivmax")
#endif

        CALL compdivmax_impl(u_f%arr, v_f%arr, w_f%arr, bp_f%arr, &
            x_f%arr, y_f%arr, z_f%arr, rddx_f%arr, rddy_f%arr, &
            rddz_f%arr, sdiv_f%arr, divmax_grid, divmax_pos_grid)

#ifdef _MGLET_PROFILE_ANNOTATIONS_
        CALL profile_range_pop()
#endif

        DO i = 1, nmygrids
            igrid = mygrids(i)
            CALL itinfo_sample(igrid, divmax=divmax_grid(i), &
                divmax_pos=divmax_pos_grid(:, i))
        END DO

        DEALLOCATE(divmax_grid)
        DEALLOCATE(divmax_pos_grid)
    END SUBROUTINE compdivmax


    SUBROUTINE compdivmax_impl(u, v, w, bp, x, y, z, rddx, rddy, rddz, &
            sdiv, divmaxgrid, divmaxposgrid)
        ! Subroutine arguments
        REAL(realk), INTENT(in) :: u(*), v(*), w(*), bp(*)
        REAL(realk), INTENT(in) :: x(*), y(*), z(*)
        REAL(realk), INTENT(in) :: rddx(*), rddy(*), rddz(*), sdiv(*)
        REAL(realk), INTENT(inout) :: divmaxgrid(*)
        REAL(realk), INTENT(inout) :: divmaxposgrid(3, *)

        ! Local variables
        INTEGER(intk) :: i, igrid, kk, jj, ii, ip3, ip1x, ip1y, ip1z
        TYPE(valpos_t) :: maxvalpos

        !$omp target teams distribute private(i, igrid, kk, jj, ii, ip3, ip1x, &
        !$omp& ip1y, ip1z, maxvalpos) &
        !$omp& map(from: divmaxgrid(1:nmygrids), &
        !$omp& divmaxposgrid(1:3, 1:nmygrids))
        DO i = 1, nmygrids
            igrid = mygrids(i)
            divmaxgrid(i) = -HUGE(1.0_realk)
            divmaxposgrid(1, i) = -HUGE(1.0_realk)
            divmaxposgrid(2, i) = -HUGE(1.0_realk)
            divmaxposgrid(3, i) = -HUGE(1.0_realk)

            CALL get_mgdims(kk, jj, ii, igrid)
            CALL get_ip3(ip3, igrid)
            CALL get_ip1x(ip1x, igrid)
            CALL get_ip1y(ip1y, igrid)
            CALL get_ip1z(ip1z, igrid)

            CALL valpos_init(maxvalpos)
            !$omp parallel
            CALL compdivmax_grid(maxvalpos, &
                kk, jj, ii, u(ip3), v(ip3), w(ip3), x(ip1x), y(ip1y), z(ip1z), &
                rddx(ip1x), rddy(ip1y), rddz(ip1z), bp(ip3), sdiv(ip3))
            !$omp end parallel

            IF (maxvalpos%val >= 0.0_realk) THEN
                divmaxgrid(i) = maxvalpos%val
                divmaxposgrid(1, i) = maxvalpos%x
                divmaxposgrid(2, i) = maxvalpos%y
                divmaxposgrid(3, i) = maxvalpos%z
            END IF
        END DO
        !$omp end target teams distribute
    END SUBROUTINE compdivmax_impl


    ! Simpler routine than "divcal", because this one does only compute the
    ! /max/ divergence and its location, and does not explicitly fill in
    ! a full 3-D divergcence field.
    SUBROUTINE compdivmax_grid(maxvalpos, kk, jj, ii, u, v, w, &
            x, y, z, rddx, rddy, rddz, bp, sdiv)
        !$omp declare target
        ! Subroutine arguments
        TYPE(valpos_t), INTENT(inout) :: maxvalpos
        INTEGER(intk), INTENT(in) :: kk, jj, ii
        REAL(realk), INTENT(in) :: u(kk, jj, ii)
        REAL(realk), INTENT(in) :: v(kk, jj, ii)
        REAL(realk), INTENT(in) :: w(kk, jj, ii)
        REAL(realk), INTENT(in) :: x(ii), y(jj), z(kk)
        REAL(realk), INTENT(in) :: rddx(ii), rddy(jj), rddz(kk)
        REAL(realk), INTENT(in) :: bp(kk, jj, ii)
        REAL(realk), INTENT(in) :: sdiv(kk, jj, ii)

        ! Local variables
        TYPE(valpos_t) :: localvalpos
        INTEGER(intk) :: k, j, i
        REAL(realk) :: div

        !$omp do collapse(3) private(i, j, k, div, localvalpos) &
        !$omp& reduction(valpos:maxvalpos)
        DO i = 3, ii-2
            DO j = 3, jj-2
                DO k = 3, kk-2
                    div = sdiv(k, j, i) &
                        + (u(k, j, i) - u(k, j, i-1))*rddx(i) &
                        + (v(k, j, i) - v(k, j-1, i))*rddy(j) &
                        + (w(k, j, i) - w(k-1, j, i))*rddz(k)
                    div = bp(k, j, i)*div

                    localvalpos%val = ABS(div)
                    localvalpos%x = x(i)
                    localvalpos%y = y(j)
                    localvalpos%z = z(k)
                    maxvalpos = get_maxvalpos(maxvalpos, localvalpos)
                END DO
            END DO
        END DO
        !$omp end do
    END SUBROUTINE compdivmax_grid


    PURE SUBROUTINE valpos_init(priv)
        !$omp declare target
        ! Subroutine arguments
        TYPE(valpos_t), INTENT(out) :: priv

        ! Local variables
        ! none...

        priv%val = -HUGE(1.0_realk)
        priv%x = HUGE(1.0_realk)
        priv%y = HUGE(1.0_realk)
        priv%z = HUGE(1.0_realk)
    END SUBROUTINE valpos_init


    PURE TYPE(valpos_t) FUNCTION get_maxvalpos(lhs, rhs) RESULT(best)
        !$omp declare target
        ! Function arguments
        TYPE(valpos_t), INTENT(in) :: lhs, rhs

        ! Local variables
        ! none...

        best = lhs

        IF (rhs%val > lhs%val) THEN
            best = rhs
        END IF
    END FUNCTION get_maxvalpos


    SUBROUTINE compenergy(u_f, v_f, w_f, g_f, dx_f, dy_f, dz_f, &
            ddx_f, ddy_f, ddz_f)
        ! Subroutine arguments
        TYPE(field_t), INTENT(in) :: u_f, v_f, w_f
        TYPE(field_t), INTENT(in) :: g_f
        TYPE(field_t), INTENT(in) :: dx_f, dy_f, dz_f
        TYPE(field_t), INTENT(in) :: ddx_f, ddy_f, ddz_f

        ! Local variables
        REAL(realk), ALLOCATABLE :: esumg_grid(:), vsumg_grid(:)
        REAL(realk), ALLOCATABLE :: esums_grid(:), vsums_grid(:)
        INTEGER(intk) :: i, igrid

        ALLOCATE(esumg_grid(nmygrids))
        ALLOCATE(esums_grid(nmygrids))
        ALLOCATE(vsumg_grid(nmygrids))
        ALLOCATE(vsums_grid(nmygrids))

        CALL compenergy_impl(u_f%arr, v_f%arr, w_f%arr, g_f%arr, dx_f%arr, &
            dy_f%arr, dz_f%arr, ddx_f%arr, ddy_f%arr, ddz_f%arr, &
            esumg_grid, esums_grid, vsumg_grid, vsums_grid)

        DO i = 1, nmygrids
            igrid = mygrids(i)
            CALL itinfo_sample(igrid, esumg=esumg_grid(i), &
                esums=esums_grid(i))
        END DO

        DEALLOCATE(esumg_grid)
        DEALLOCATE(vsumg_grid)
        DEALLOCATE(esums_grid)
        DEALLOCATE(vsums_grid)
    END SUBROUTINE compenergy


    SUBROUTINE compenergy_impl(u, v, w, g, dx, dy, dz, ddx, ddy, ddz, &
            esumggrid, esumsgrid, vsumggrid, vsumsgrid)
        ! Subroutine arguments
        REAL(realk), INTENT(in) :: u(*), v(*), w(*), g(*)
        REAL(realk), INTENT(in) :: dx(*), dy(*), dz(*)
        REAL(realk), INTENT(in) :: ddx(*), ddy(*), ddz(*)
        REAL(realk), INTENT(inout) :: esumggrid(*), esumsgrid(*)
        REAL(realk), INTENT(inout) :: vsumggrid(*), vsumsgrid(*)

        ! Local variables
        INTEGER(intk) :: i, igrid, kk, jj, ii, ip3, ip1x, ip1y, ip1z

#ifdef _MGLET_PROFILE_ANNOTATIONS_
        CALL profile_range_push("compenergy")
#endif

        !$omp target teams distribute private(i, igrid, kk, jj, ii, ip3, ip1x, &
        !$omp& ip1y, ip1z) &
        !$omp& map(from: esumggrid(1:nmygrids), esumsgrid(1:nmygrids)) &
        !$omp& map(alloc: vsumggrid(1:nmygrids), vsumsgrid(1:nmygrids))
        DO i = 1, nmygrids
            igrid = mygrids(i)
            esumggrid(i) = 0.0_realk
            esumsgrid(i) = 0.0_realk
            vsumggrid(i) = 0.0_realk
            vsumsgrid(i) = 0.0_realk

            CALL get_mgdims(kk, jj, ii, igrid)
            CALL get_ip3(ip3, igrid)
            CALL get_ip1x(ip1x, igrid)
            CALL get_ip1y(ip1y, igrid)
            CALL get_ip1z(ip1z, igrid)

            !$omp parallel
            CALL enerfg_grid(esumggrid(i), vsumggrid(i), kk, jj, ii, &
                u(ip3), v(ip3), w(ip3), dx(ip1x), dy(ip1y), dz(ip1z), &
                ddx(ip1x), ddy(ip1y), ddz(ip1z))
            CALL enerfs_grid(esumsgrid(i), vsumsgrid(i), kk, jj, ii, g(ip3), &
                ddx(ip1x), ddy(ip1y), ddz(ip1z))
            !$omp end parallel

            IF (vsumggrid(i) > 0.0_realk) THEN
                esumggrid(i) = esumggrid(i)/vsumggrid(i)
            END IF
            IF (vsumsgrid(i) > 0.0_realk) THEN
                esumsgrid(i) = esumsgrid(i)/vsumsgrid(i)
            END IF
        END DO
        !$omp end target teams distribute

#ifdef _MGLET_PROFILE_ANNOTATIONS_
        CALL profile_range_pop()
#endif

    END SUBROUTINE compenergy_impl


    SUBROUTINE enerfg_grid(esum, vsum, kk, jj, ii, u, v, w, dx, dy, dz, &
            ddx, ddy, ddz)
        !$omp declare target
        ! Subroutine arguments
        REAL(realk), INTENT(inout) :: esum, vsum
        INTEGER(intk), INTENT(in) :: kk, jj, ii
        REAL(realk), INTENT(in) :: u(kk, jj, ii), v(kk, jj, ii), w(kk, jj, ii)
        REAL(realk), INTENT(in) :: dx(ii), dy(jj), dz(kk)
        REAL(realk), INTENT(in) :: ddx(ii), ddy(jj), ddz(kk)

        ! Local variables
        REAL(realk) :: up, vp, wp, vol
        INTEGER(intk) :: k, j, i

        !$omp do collapse(3) private(i, j, k, up, vp, wp, vol) &
        !$omp& reduction(+:esum, vsum)
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
        !$omp end do
    END SUBROUTINE enerfg_grid


    SUBROUTINE enerfs_grid(esum, vsum, kk, jj, ii, g, ddx, ddy, ddz)
        !$omp declare target
        ! Subroutine arguments
        REAL(realk), INTENT(inout) :: esum, vsum
        INTEGER(intk), INTENT(in) :: kk, jj, ii
        REAL(realk), INTENT(in) :: g(kk, jj, ii)
        REAL(realk), INTENT(in) :: ddx(ii), ddy(jj), ddz(kk)

        ! Local variables
        REAL(realk), PARAMETER :: conv2s = 0.094
        REAL(realk) :: delta, ener, vol
        INTEGER(intk) :: k, j, i

        !$omp do collapse(3) private(i, j, k, delta, ener, vol) &
        !$omp& reduction(+:esum, vsum)
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
        !$omp end do
    END SUBROUTINE enerfs_grid


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
END MODULE timeintegration_mod
