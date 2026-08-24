MODULE tstle4_mod
    USE core_mod
    USE flowcore_mod
    USE lesmodel_mod, ONLY: ilesmodel
    USE wernerwengle_mod, ONLY: tauwin

    IMPLICIT NONE(type, external)
    PRIVATE

    PUBLIC :: tstle4

    INTERFACE
        SUBROUTINE tstle4_kon_u_c(kk, jj, ii, uo, u, v, w, ut, vt, wt, dx, &
            dy, dz, ddx, ddy, ddz, rdx, rdy, rdz, rddx, rddy, rddz, nfro, &
            nbac) BIND(C, name="tstle4_kon_u_c")
            IMPORT :: c_intk, c_realk
            INTEGER(c_intk), VALUE, INTENT(in) :: kk, jj, ii
            REAL(c_realk), INTENT(inout) :: uo(*)
            REAL(c_realk), INTENT(in) :: u(*), v(*), w(*), ut(*), vt(*), wt(*)
            REAL(c_realk), INTENT(in) :: dx(*), dy(*), dz(*)
            REAL(c_realk), INTENT(in) :: ddx(*), ddy(*), ddz(*)
            REAL(c_realk), INTENT(in) :: rdx(*), rdy(*), rdz(*)
            REAL(c_realk), INTENT(in) :: rddx(*), rddy(*), rddz(*)
            INTEGER(c_intk), VALUE, INTENT(in) :: nfro, nbac
        END SUBROUTINE tstle4_kon_u_c

        SUBROUTINE tstle4_kon_v_c(kk, jj, ii, vo, u, v, w, ut, vt, wt, dx, &
            dy, dz, ddx, ddy, ddz, rdx, rdy, rdz, rddx, rddy, rddz, nrgt, &
            nlft) BIND(C, name="tstle4_kon_v_c")
            IMPORT :: c_intk, c_realk
            INTEGER(c_intk), VALUE, INTENT(in) :: kk, jj, ii
            REAL(c_realk), INTENT(inout) :: vo(*)
            REAL(c_realk), INTENT(in) :: u(*), v(*), w(*), ut(*), vt(*), wt(*)
            REAL(c_realk), INTENT(in) :: dx(*), dy(*), dz(*)
            REAL(c_realk), INTENT(in) :: ddx(*), ddy(*), ddz(*)
            REAL(c_realk), INTENT(in) :: rdx(*), rdy(*), rdz(*)
            REAL(c_realk), INTENT(in) :: rddx(*), rddy(*), rddz(*)
            INTEGER(c_intk), VALUE, INTENT(in) :: nrgt, nlft
        END SUBROUTINE tstle4_kon_v_c

        SUBROUTINE tstle4_kon_w_c(kk, jj, ii, wo, u, v, w, ut, vt, wt, dx, &
            dy, dz, ddx, ddy, ddz, rdx, rdy, rdz, rddx, rddy, rddz, nbot, &
            ntop) BIND(C, name="tstle4_kon_w_c")
            IMPORT :: c_intk, c_realk
            INTEGER(c_intk), VALUE, INTENT(in) :: kk, jj, ii
            REAL(c_realk), INTENT(inout) :: wo(*)
            REAL(c_realk), INTENT(in) :: u(*), v(*), w(*), ut(*), vt(*), wt(*)
            REAL(c_realk), INTENT(in) :: dx(*), dy(*), dz(*)
            REAL(c_realk), INTENT(in) :: ddx(*), ddy(*), ddz(*)
            REAL(c_realk), INTENT(in) :: rdx(*), rdy(*), rdz(*)
            REAL(c_realk), INTENT(in) :: rddx(*), rddy(*), rddz(*)
            INTEGER(c_intk), VALUE, INTENT(in) :: nbot, ntop
        END SUBROUTINE tstle4_kon_w_c

        SUBROUTINE tstle4_diff_u_c(kk, jj, ii, uo, u, v, w, g, dx, dy, dz, &
            ddx, ddy, ddz, rdx, rdy, rdz, rddx, rddy, rddz, nfro, nbac, &
            ilesmodel_in, gmol_in, rho_in) BIND(C, name="tstle4_diff_u_c")
            IMPORT :: c_intk, c_realk
            INTEGER(c_intk), VALUE, INTENT(in) :: kk, jj, ii
            REAL(c_realk), INTENT(inout) :: uo(*)
            REAL(c_realk), INTENT(in) :: u(*), v(*), w(*), g(*)
            REAL(c_realk), INTENT(in) :: dx(*), dy(*), dz(*)
            REAL(c_realk), INTENT(in) :: ddx(*), ddy(*), ddz(*)
            REAL(c_realk), INTENT(in) :: rdx(*), rdy(*), rdz(*)
            REAL(c_realk), INTENT(in) :: rddx(*), rddy(*), rddz(*)
            INTEGER(c_intk), VALUE, INTENT(in) :: nfro, nbac, ilesmodel_in
            REAL(c_realk), VALUE, INTENT(in) :: gmol_in, rho_in
        END SUBROUTINE tstle4_diff_u_c

        SUBROUTINE tstle4_diff_v_c(kk, jj, ii, vo, u, v, w, g, dx, dy, dz, &
            ddx, ddy, ddz, rdx, rdy, rdz, rddx, rddy, rddz, nrgt, nlft, &
            ilesmodel_in, gmol_in, rho_in) BIND(C, name="tstle4_diff_v_c")
            IMPORT :: c_intk, c_realk
            INTEGER(c_intk), VALUE, INTENT(in) :: kk, jj, ii
            REAL(c_realk), INTENT(inout) :: vo(*)
            REAL(c_realk), INTENT(in) :: u(*), v(*), w(*), g(*)
            REAL(c_realk), INTENT(in) :: dx(*), dy(*), dz(*)
            REAL(c_realk), INTENT(in) :: ddx(*), ddy(*), ddz(*)
            REAL(c_realk), INTENT(in) :: rdx(*), rdy(*), rdz(*)
            REAL(c_realk), INTENT(in) :: rddx(*), rddy(*), rddz(*)
            INTEGER(c_intk), VALUE, INTENT(in) :: nrgt, nlft, ilesmodel_in
            REAL(c_realk), VALUE, INTENT(in) :: gmol_in, rho_in
        END SUBROUTINE tstle4_diff_v_c

        SUBROUTINE tstle4_diff_w_c(kk, jj, ii, wo, u, v, w, g, dx, dy, dz, &
            ddx, ddy, ddz, rdx, rdy, rdz, rddx, rddy, rddz, nbot, ntop, &
            ilesmodel_in, gmol_in, rho_in) BIND(C, name="tstle4_diff_w_c")
            IMPORT :: c_intk, c_realk
            INTEGER(c_intk), VALUE, INTENT(in) :: kk, jj, ii
            REAL(c_realk), INTENT(inout) :: wo(*)
            REAL(c_realk), INTENT(in) :: u(*), v(*), w(*), g(*)
            REAL(c_realk), INTENT(in) :: dx(*), dy(*), dz(*)
            REAL(c_realk), INTENT(in) :: ddx(*), ddy(*), ddz(*)
            REAL(c_realk), INTENT(in) :: rdx(*), rdy(*), rdz(*)
            REAL(c_realk), INTENT(in) :: rddx(*), rddy(*), rddz(*)
            INTEGER(c_intk), VALUE, INTENT(in) :: nbot, ntop, ilesmodel_in
            REAL(c_realk), VALUE, INTENT(in) :: gmol_in, rho_in
        END SUBROUTINE tstle4_diff_w_c

        SUBROUTINE tstle4_diff_swc_c(kk, jj, ii, uo, vo, wo, u, v, w, ddx, &
            ddy, ddz, nfro, nbac, nrgt, nlft, nbot, ntop, gmol_in, rho_in) &
            BIND(C, name="tstle4_diff_swc_c")
            IMPORT :: c_intk, c_realk
            INTEGER(c_intk), VALUE, INTENT(in) :: kk, jj, ii
            REAL(c_realk), INTENT(inout) :: uo(*), vo(*), wo(*)
            REAL(c_realk), INTENT(in) :: u(*), v(*), w(*)
            REAL(c_realk), INTENT(in) :: ddx(*), ddy(*), ddz(*)
            INTEGER(c_intk), VALUE, INTENT(in) :: nfro, nbac, nrgt, nlft
            INTEGER(c_intk), VALUE, INTENT(in) :: nbot, ntop
            REAL(c_realk), VALUE, INTENT(in) :: gmol_in, rho_in
        END SUBROUTINE tstle4_diff_swc_c
    END INTERFACE

    !$omp declare target(tstle4_kon_u_c, tstle4_kon_v_c, tstle4_kon_w_c, &
    !$omp&   tstle4_diff_u_c, tstle4_diff_v_c, tstle4_diff_w_c, &
    !$omp&   tstle4_diff_swc_c)

CONTAINS
    SUBROUTINE tstle4(uo_f, vo_f, wo_f, u_f, v_f, w_f, ut_f, vt_f, wt_f, &
            p_f, g_f)
        ! Subroutine arguments
        TYPE(field_t), INTENT(inout) :: uo_f
        TYPE(field_t), INTENT(inout) :: vo_f
        TYPE(field_t), INTENT(inout) :: wo_f
        TYPE(field_t), INTENT(in) :: u_f
        TYPE(field_t), INTENT(in) :: v_f
        TYPE(field_t), INTENT(in) :: w_f
        TYPE(field_t), INTENT(in) :: ut_f
        TYPE(field_t), INTENT(in) :: vt_f
        TYPE(field_t), INTENT(in) :: wt_f
        TYPE(field_t), INTENT(in) :: p_f
        TYPE(field_t), INTENT(in) :: g_f

        ! Local variables
        TYPE(field_t), POINTER :: dx_f, dy_f, dz_f, ddx_f, ddy_f, ddz_f
        TYPE(field_t), POINTER :: rdx_f, rdy_f, rdz_f, rddx_f, rddy_f, rddz_f
        TYPE(field_t), POINTER :: wcu_f, wcv_f, wcw_f
        CALL start_timer(310)

        CALL set_field_arr(uo_f, 0.0_realk, device=.TRUE.)
        CALL set_field_arr(vo_f, 0.0_realk, device=.TRUE.)
        CALL set_field_arr(wo_f, 0.0_realk, device=.TRUE.)

        CALL get_field(dx_f, "DX")
        CALL get_field(dy_f, "DY")
        CALL get_field(dz_f, "DZ")

        CALL get_field(ddx_f, "DDX")
        CALL get_field(ddy_f, "DDY")
        CALL get_field(ddz_f, "DDZ")

        CALL get_field(rdx_f, "RDX")
        CALL get_field(rdy_f, "RDY")
        CALL get_field(rdz_f, "RDZ")

        CALL get_field(rddx_f, "RDDX")
        CALL get_field(rddy_f, "RDDY")
        CALL get_field(rddz_f, "RDDZ")

        CALL push_field(wcu_f, "TSTLE4_WCU")
        CALL push_field(wcv_f, "TSTLE4_WCV")
        CALL push_field(wcw_f, "TSTLE4_WCW")

        CALL start_timer(311)
        CALL tstle4_kon_impl(uo_f%arr, vo_f%arr, wo_f%arr, u_f%arr, &
            v_f%arr, w_f%arr, ut_f%arr, vt_f%arr, wt_f%arr, dx_f%arr, &
            dy_f%arr, dz_f%arr, ddx_f%arr, ddy_f%arr, ddz_f%arr, &
            rdx_f%arr, rdy_f%arr, rdz_f%arr, rddx_f%arr, rddy_f%arr, &
            rddz_f%arr)
        CALL stop_timer(311)

        CALL start_timer(312)
        CALL tstle4_diff_impl(uo_f%arr, vo_f%arr, wo_f%arr, u_f%arr, &
            v_f%arr, w_f%arr, g_f%arr, dx_f%arr, dy_f%arr, dz_f%arr, &
            ddx_f%arr, ddy_f%arr, ddz_f%arr, rdx_f%arr, rdy_f%arr, &
            rdz_f%arr, rddx_f%arr, rddy_f%arr, rddz_f%arr)
        CALL stop_timer(312)

        CALL start_timer(313)
        CALL tstle4_gradp_impl(uo_f%arr, vo_f%arr, wo_f%arr, p_f%arr, &
            dx_f%arr, dy_f%arr, dz_f%arr)
        CALL stop_timer(313)

        CALL start_timer(314)
        CALL tstle4_par_impl(uo_f%arr, vo_f%arr, wo_f%arr, u_f%arr, &
            v_f%arr, w_f%arr, ut_f%arr, vt_f%arr, wt_f%arr, dx_f%arr, &
            dy_f%arr, dz_f%arr, ddx_f%arr, ddy_f%arr, ddz_f%arr, &
            rdx_f%arr, rdy_f%arr, rdz_f%arr, rddx_f%arr, rddy_f%arr, &
            rddz_f%arr, wcu_f%arr, wcv_f%arr, wcw_f%arr)
        CALL stop_timer(314)

        CALL pop_field(wcw_f)
        CALL pop_field(wcv_f)
        CALL pop_field(wcu_f)

        CALL stop_timer(310)
    END SUBROUTINE tstle4


    SUBROUTINE tstle4_kon_impl(uo, vo, wo, u, v, w, ut, vt, wt, dx, dy, &
            dz, ddx, ddy, ddz, rdx, rdy, rdz, rddx, rddy, rddz)
        REAL(realk), INTENT(inout) :: uo(*), vo(*), wo(*)
        REAL(realk), INTENT(in) :: u(*), v(*), w(*), ut(*), vt(*), wt(*)
        REAL(realk), INTENT(in) :: dx(*), dy(*), dz(*), ddx(*), ddy(*), ddz(*)
        REAL(realk), INTENT(in) :: rdx(*), rdy(*), rdz(*)
        REAL(realk), INTENT(in) :: rddx(*), rddy(*), rddz(*)

        INTEGER(intk) :: i, igrid, ip3, ipx, ipy, ipz
        INTEGER(intk) :: kk, jj, ii, nfro, nbac, nrgt, nlft, nbot, ntop

        !$omp target teams distribute &
        !$omp& private(i, igrid, ip3, ipx, ipy, ipz, &
        !$omp&  kk, jj, ii, nfro, nbac, nrgt, nlft, nbot, ntop)
        DO i = 1, nmygrids
            igrid = mygrids(i)

            CALL get_mgdims(kk, jj, ii, igrid)
            CALL get_mgbasb(nfro, nbac, nrgt, nlft, nbot, ntop, igrid)
            CALL get_ip3(ip3, igrid)
            CALL get_ip1x(ipx, igrid)
            CALL get_ip1y(ipy, igrid)
            CALL get_ip1z(ipz, igrid)

            !$omp parallel
            CALL tstle4_kon_u_c(INT(kk, c_intk), INT(jj, c_intk), &
                INT(ii, c_intk), uo(ip3), u(ip3), v(ip3), w(ip3), ut(ip3), &
                vt(ip3), wt(ip3), dx(ipx), dy(ipy), dz(ipz), ddx(ipx), &
                ddy(ipy), ddz(ipz), rdx(ipx), rdy(ipy), rdz(ipz), rddx(ipx), &
                rddy(ipy), rddz(ipz), INT(nfro, c_intk), INT(nbac, c_intk))

            CALL tstle4_kon_v_c(INT(kk, c_intk), INT(jj, c_intk), &
                INT(ii, c_intk), vo(ip3), u(ip3), v(ip3), w(ip3), ut(ip3), &
                vt(ip3), wt(ip3), dx(ipx), dy(ipy), dz(ipz), ddx(ipx), &
                ddy(ipy), ddz(ipz), rdx(ipx), rdy(ipy), rdz(ipz), rddx(ipx), &
                rddy(ipy), rddz(ipz), INT(nrgt, c_intk), INT(nlft, c_intk))

            CALL tstle4_kon_w_c(INT(kk, c_intk), INT(jj, c_intk), &
                INT(ii, c_intk), wo(ip3), u(ip3), v(ip3), w(ip3), ut(ip3), &
                vt(ip3), wt(ip3), dx(ipx), dy(ipy), dz(ipz), ddx(ipx), &
                ddy(ipy), ddz(ipz), rdx(ipx), rdy(ipy), rdz(ipz), rddx(ipx), &
                rddy(ipy), rddz(ipz), INT(nbot, c_intk), INT(ntop, c_intk))
            !$omp end parallel
        END DO
        !$omp end target teams distribute

    END SUBROUTINE tstle4_kon_impl


    SUBROUTINE tstle4_diff_impl(uo, vo, wo, u, v, w, g, dx, dy, dz, ddx, &
            ddy, ddz, rdx, rdy, rdz, rddx, rddy, rddz)
        REAL(realk), INTENT(inout) :: uo(*), vo(*), wo(*)
        REAL(realk), INTENT(in) :: u(*), v(*), w(*), g(*)
        REAL(realk), INTENT(in) :: dx(*), dy(*), dz(*), ddx(*), ddy(*), ddz(*)
        REAL(realk), INTENT(in) :: rdx(*), rdy(*), rdz(*)
        REAL(realk), INTENT(in) :: rddx(*), rddy(*), rddz(*)

        INTEGER(intk) :: i, igrid, ip3, ipx, ipy, ipz
        INTEGER(intk) :: kk, jj, ii, nfro, nbac, nrgt, nlft, nbot, ntop

        !$omp target teams distribute &
        !$omp& private(i, igrid, ip3, ipx, ipy, ipz, &
        !$omp&  kk, jj, ii, nfro, nbac, nrgt, nlft, nbot, ntop)
        DO i = 1, nmygrids
            igrid = mygrids(i)

            CALL get_mgdims(kk, jj, ii, igrid)
            CALL get_mgbasb(nfro, nbac, nrgt, nlft, nbot, ntop, igrid)
            CALL get_ip3(ip3, igrid)
            CALL get_ip1x(ipx, igrid)
            CALL get_ip1y(ipy, igrid)
            CALL get_ip1z(ipz, igrid)

            !$omp parallel
            CALL tstle4_diff_swc_c(INT(kk, c_intk), INT(jj, c_intk), &
                INT(ii, c_intk), uo(ip3), vo(ip3), wo(ip3), u(ip3), v(ip3), &
                w(ip3), ddx(ipx), ddy(ipy), ddz(ipz), INT(nfro, c_intk), &
                INT(nbac, c_intk), INT(nrgt, c_intk), INT(nlft, c_intk), &
                INT(nbot, c_intk), INT(ntop, c_intk), REAL(gmol, c_realk), &
                REAL(rho, c_realk))

            CALL tstle4_diff_u_c(INT(kk, c_intk), INT(jj, c_intk), &
                INT(ii, c_intk), uo(ip3), u(ip3), v(ip3), w(ip3), g(ip3), &
                dx(ipx), dy(ipy), dz(ipz), ddx(ipx), ddy(ipy), ddz(ipz), &
                rdx(ipx), rdy(ipy), rdz(ipz), rddx(ipx), rddy(ipy), rddz(ipz), &
                INT(nfro, c_intk), INT(nbac, c_intk), INT(ilesmodel, c_intk), &
                REAL(gmol, c_realk), REAL(rho, c_realk))

            CALL tstle4_diff_v_c(INT(kk, c_intk), INT(jj, c_intk), &
                INT(ii, c_intk), vo(ip3), u(ip3), v(ip3), w(ip3), g(ip3), &
                dx(ipx), dy(ipy), dz(ipz), ddx(ipx), ddy(ipy), ddz(ipz), &
                rdx(ipx), rdy(ipy), rdz(ipz), rddx(ipx), rddy(ipy), rddz(ipz), &
                INT(nrgt, c_intk), INT(nlft, c_intk), INT(ilesmodel, c_intk), &
                REAL(gmol, c_realk), REAL(rho, c_realk))

            CALL tstle4_diff_w_c(INT(kk, c_intk), INT(jj, c_intk), &
                INT(ii, c_intk), wo(ip3), u(ip3), v(ip3), w(ip3), g(ip3), &
                dx(ipx), dy(ipy), dz(ipz), ddx(ipx), ddy(ipy), ddz(ipz), &
                rdx(ipx), rdy(ipy), rdz(ipz), rddx(ipx), rddy(ipy), rddz(ipz), &
                INT(nbot, c_intk), INT(ntop, c_intk), INT(ilesmodel, c_intk), &
                REAL(gmol, c_realk), REAL(rho, c_realk))
            !$omp end parallel

        END DO
        !$omp end target teams distribute

    END SUBROUTINE tstle4_diff_impl


    SUBROUTINE tstle4_gradp_impl(uo, vo, wo, p, dx, dy, dz)
        REAL(realk), INTENT(inout) :: uo(*), vo(*), wo(*)
        REAL(realk), INTENT(in) :: p(*), dx(*), dy(*), dz(*)

        INTEGER(intk) :: i, igrid, ip3, ipx, ipy, ipz
        INTEGER(intk) :: kk, jj, ii, nfro, nbac, nrgt, nlft, nbot, ntop

        !$omp target teams distribute &
        !$omp& private(i, igrid, ip3, ipx, ipy, ipz, &
        !$omp&  kk, jj, ii, nfro, nbac, nrgt, nlft, nbot, ntop)
        DO i = 1, nmygrids
            igrid = mygrids(i)

            CALL get_mgdims(kk, jj, ii, igrid)
            CALL get_mgbasb(nfro, nbac, nrgt, nlft, nbot, ntop, igrid)
            CALL get_ip3(ip3, igrid)
            CALL get_ip1x(ipx, igrid)
            CALL get_ip1y(ipy, igrid)
            CALL get_ip1z(ipz, igrid)

            !$omp parallel
            CALL tstle4_gradp(kk, jj, ii, uo(ip3), vo(ip3), wo(ip3), p(ip3), &
                dx(ipx), dy(ipy), dz(ipz), &
                nfro, nbac, nrgt, nlft, nbot, ntop, igrid)
            !$omp end parallel
        END DO
        !$omp end target teams distribute
    END SUBROUTINE tstle4_gradp_impl


    SUBROUTINE tstle4_par_impl(uo, vo, wo, u, v, w, ut, vt, wt, dx, dy, dz, &
            ddx, ddy, ddz, rdx, rdy, rdz, rddx, rddy, rddz, wcu, wcv, wcw)
        REAL(realk), INTENT(inout) :: uo(*), vo(*), wo(*)
        REAL(realk), INTENT(inout) :: wcu(*), wcv(*), wcw(*)
        REAL(realk), INTENT(in) :: u(*), v(*), w(*), ut(*), vt(*), wt(*)
        REAL(realk), INTENT(in) :: dx(*), dy(*), dz(*), ddx(*), ddy(*), ddz(*)
        REAL(realk), INTENT(in) :: rdx(*), rdy(*), rdz(*)
        REAL(realk), INTENT(in) :: rddx(*), rddy(*), rddz(*)

        INTEGER(intk) :: i, igrid, ip3, ipx, ipy, ipz
        INTEGER(intk) :: kk, jj, ii, nfro, nbac, nrgt, nlft, nbot, ntop

        !$omp target teams distribute &
        !$omp& private(i, igrid, ip3, ipx, ipy, ipz, &
        !$omp&  kk, jj, ii, nfro, nbac, nrgt, nlft, nbot, ntop)
        DO i = 1, nmygrids
            igrid = mygrids(i)

            CALL get_mgdims(kk, jj, ii, igrid)
            CALL get_mgbasb(nfro, nbac, nrgt, nlft, nbot, ntop, igrid)
            CALL get_ip3(ip3, igrid)
            CALL get_ip1x(ipx, igrid)
            CALL get_ip1y(ipy, igrid)
            CALL get_ip1z(ipz, igrid)

            !$omp parallel
            CALL tstle4_par(kk, jj, ii, uo(ip3), vo(ip3), wo(ip3), u(ip3), &
                v(ip3), w(ip3), ut(ip3), vt(ip3), wt(ip3), dx(ipx), dy(ipy), &
                dz(ipz), ddx(ipx), ddy(ipy), ddz(ipz), rdx(ipx), rdy(ipy), &
                rdz(ipz), rddx(ipx), rddy(ipy), rddz(ipz), &
                wcu(ip3), wcv(ip3), wcw(ip3), &
                nfro, nbac, nrgt, nlft, nbot, ntop)
            !$omp end parallel
        END DO
        !$omp end target teams distribute
    END SUBROUTINE tstle4_par_impl


    ! The convective terms are computed in two steps:
    ! First, the mass fluxes (transporting velocities) are interpolated to
    ! the faces of the momentum cell. This interpolation is performed in a
    ! way which ensures mass conservation at the momentum cell if the
    ! velocity field is divergence-free on the adjacent pressure cells.
    ! Second, the transported velocities are interpolated in a symmetry-
    ! preserving manner (the convective term has to be skew-symmetric in
    ! order to conserve energy).
    !
    ! Details can be found in:
    ! [1] Heinz Werner, Grobstruktursimulation der turbulenten Strömung
    !     über eine querliegende Rippe in einem Plattenkanal bei hoher
    !     Reynolds-Zahl, PhD Thesis, Technical University of Munich, 1991
    ! [2] Verstappen et al., SYMMETRY-PRESERVING DISCRETIZATIONS OF THE
    !     INCOMPRESSIBLE NAVIER-STOKES EQUATIONS, European Conference on
    !     Computational Fluid Dynamics, ECCOMAS CFD 2006
    SUBROUTINE tstle4_kon_u(kk, jj, ii, uo, u, v, w, ut, vt, wt, dx, dy, dz, &
            ddx, ddy, ddz, rdx, rdy, rdz, rddx, rddy, rddz, nfro, nbac, ntop)
        !$omp declare target
        ! Subroutine arguments
        INTEGER(intk), INTENT(in) :: kk, jj, ii
        REAL(realk), INTENT(inout) :: uo(kk, jj, ii)
        REAL(realk), INTENT(in) :: u(kk, jj, ii), v(kk, jj, ii), w(kk, jj, ii)
        REAL(realk), INTENT(in) :: ut(kk, jj, ii), vt(kk, jj, ii), &
            wt(kk, jj, ii)
        REAL(realk), INTENT(in) :: dx(ii), dy(jj), dz(kk)
        REAL(realk), INTENT(in) :: ddx(ii), ddy(jj), ddz(kk)
        REAL(realk), INTENT(in) :: rdx(ii), rdy(jj), rdz(kk)
        REAL(realk), INTENT(in) :: rddx(ii), rddy(jj), rddz(kk)
        INTEGER, INTENT(in) :: nfro, nbac, ntop

        ! Local variables
        INTEGER(intk) :: k, j, i
        INTEGER(intk) :: nbu, nfu
        REAL(realk) :: ax, ay, az
        REAL(realk) :: fw, fe, ft, fb, fn, fs
        REAL(realk) :: qw, qe, qt, qb, qn, qs

        nfu = 0
        nbu = 0

        ! CON = 7
        IF (nbac == 7) nbu = 1

        ! OP1 = 3
        IF (nfro == 3) nfu = 1
        IF (nbac == 3) nbu = 1

        !$omp do collapse(3) private(i, j, k, ax, ay, az, fe, fw, &
        !$omp& fn, fs, ft, fb, qe, qw, qn, qs, qt, qb)
        DO i = 3-nfu, ii-3+nbu
            DO j = 3, jj-2
                DO k = 3, kk-2
                    ax = ddy(j)*ddz(k)
                    ay = dx(i)*ddz(k)
                    az = dx(i)*ddy(j)

                    fe = ax*(ut(k, j, i) + (ut(k, j, i+1) - ut(k, j, i)) &
                        * 0.5*dx(i)/ddx(i+1))
                    fw = ax*(ut(k, j, i-1) + (ut(k, j, i) - ut(k, j, i-1)) &
                        * 0.5*dx(i-1)/ddx(i))
                    fn = ay*(vt(k, j, i) + vt(k, j, i+1))*0.5
                    fs = ay*(vt(k, j-1, i) + vt(k, j-1, i+1))*0.5
                    ft = az*(wt(k, j, i) + wt(k, j, i+1))*0.5
                    fb = az*(wt(k-1, j, i) + wt(k-1, j, i+1))*0.5

                    qe = 0.5*fe*(u(k, j, i) + u(k, j, i+1))
                    qw = 0.5*fw*(u(k, j, i-1) + u(k, j, i))
                    qn = 0.5*fn*(u(k, j, i) + u(k, j+1, i))
                    qs = 0.5*fs*(u(k, j-1, i) + u(k, j, i))
                    qt = 0.5*ft*(u(k, j, i) + u(k+1, j, i))
                    qb = 0.5*fb*(u(k-1, j, i) + u(k, j, i))

                    uo(k, j, i) = -(qe-qw+qn-qs+qt-qb)*rdx(i)*rddy(j)*rddz(k)
                END DO
            END DO
        END DO
        !$omp end do

    END SUBROUTINE tstle4_kon_u


    SUBROUTINE tstle4_kon_v(kk, jj, ii, vo, u, v, w, ut, vt, wt, dx, dy, dz, &
            ddx, ddy, ddz, rdx, rdy, rdz, rddx, rddy, rddz, nfro, nbac, nrgt, &
            nlft, nbot, ntop)
        !$omp declare target
        ! Subroutine arguments
        INTEGER(intk), INTENT(in) :: kk, jj, ii
        REAL(realk), INTENT(inout) :: vo(kk, jj, ii)
        REAL(realk), INTENT(in) :: u(kk, jj, ii), v(kk, jj, ii), w(kk, jj, ii)
        REAL(realk), INTENT(in) :: ut(kk, jj, ii), vt(kk, jj, ii), &
            wt(kk, jj, ii)
        REAL(realk), INTENT(in) :: dx(ii), dy(jj), dz(kk)
        REAL(realk), INTENT(in) :: ddx(ii), ddy(jj), ddz(kk)
        REAL(realk), INTENT(in) :: rdx(ii), rdy(jj), rdz(kk)
        REAL(realk), INTENT(in) :: rddx(ii), rddy(jj), rddz(kk)
        INTEGER, INTENT(in) :: nfro, nbac, nrgt, nlft, nbot, ntop

        ! Local variables
        INTEGER(intk) :: k, j, i
        INTEGER(intk) :: nrv, nlv
        REAL(realk) :: ax, ay, az
        REAL(realk) :: fw, fe, ft, fb, fn, fs
        REAL(realk) :: qw, qe, qt, qb, qn, qs

        nrv = 0
        nlv = 0

        ! CON = 7
        IF (nlft == 7) nlv = 1

        ! OP1 = 3
        IF (nrgt == 3) nrv = 1
        IF (nlft == 3) nlv = 1

        !$omp do collapse(3) private(i, j, k, ax, ay, az, fe, fw, &
        !$omp& fn, fs, ft, fb, qe, qw, qn, qs, qt, qb)
        DO i = 3, ii-2
            DO j = 3-nrv, jj-3+nlv
                DO k = 3, kk-2
                    ax = dy(j)*ddz(k)
                    ay = ddx(i)*ddz(k)
                    az = ddx(i)*dy(j)

                    fe = ax*(ut(k, j, i) + ut(k, j+1, i))*0.5
                    fw = ax*(ut(k, j, i-1) + ut(k, j+1, i-1))*0.5
                    fn = ay*(vt(k, j, i) + (vt(k, j+1, i) - vt(k, j, i)) &
                        * 0.5*dy(j)/ddy(j+1))
                    fs = ay*(vt(k, j-1, i) + (vt(k, j, i) -vt(k, j-1, i)) &
                        * 0.5*dy(j-1)/ddy(j))
                    ft = az*(wt(k, j, i) + wt(k, j+1, i))*0.5
                    fb = az*(wt(k-1, j, i) + wt(k-1, j+1, i))*0.5

                    qe = 0.5*fe*(v(k, j, i) + v(k, j, i+1))
                    qw = 0.5*fw*(v(k, j, i-1) + v(k, j, i))
                    qn = 0.5*fn*(v(k, j, i) + v(k, j+1, i))
                    qs = 0.5*fs*(v(k, j-1, i) + v(k, j, i))
                    qt = 0.5*ft*(v(k, j, i) + v(k+1, j, i))
                    qb = 0.5*fb*(v(k-1, j, i) + v(k, j, i))

                    vo(k, j, i) = -(qe-qw+qn-qs+qt-qb)*rddx(i)*rdy(j)*rddz(k)
                END DO
            END DO
        END DO
        !$omp end do

    END SUBROUTINE tstle4_kon_v


    SUBROUTINE tstle4_kon_w(kk, jj, ii, wo, u, v, w, ut, vt, wt, dx, dy, dz, &
            ddx, ddy, ddz, rdx, rdy, rdz, rddx, rddy, rddz, nfro, nbac, nrgt, &
            nlft, nbot, ntop)
        !$omp declare target
        ! Subroutine arguments
        INTEGER(intk), INTENT(in) :: kk, jj, ii
        REAL(realk), INTENT(inout) :: wo(kk, jj, ii)
        REAL(realk), INTENT(in) :: u(kk, jj, ii), v(kk, jj, ii), w(kk, jj, ii)
        REAL(realk), INTENT(in) :: ut(kk, jj, ii), vt(kk, jj, ii), &
            wt(kk, jj, ii)
        REAL(realk), INTENT(in) :: dx(ii), dy(jj), dz(kk)
        REAL(realk), INTENT(in) :: ddx(ii), ddy(jj), ddz(kk)
        REAL(realk), INTENT(in) :: rdx(ii), rdy(jj), rdz(kk)
        REAL(realk), INTENT(in) :: rddx(ii), rddy(jj), rddz(kk)
        INTEGER, INTENT(in) :: nfro, nbac, nrgt, nlft, nbot, ntop

        ! Local variables
        INTEGER(intk) :: k, j, i
        INTEGER(intk) :: nbw, ntw
        REAL(realk) :: ax, ay, az
        REAL(realk) :: fw, fe, ft, fb, fn, fs
        REAL(realk) :: qw, qe, qt, qb, qn, qs

        nbw = 0
        ntw = 0

        ! CON = 7
        IF (ntop == 7) ntw = 1

        ! OP1 = 3
        IF (nbot == 3) nbw = 1
        IF (ntop == 3) ntw = 1

        !$omp do collapse(3) private(i, j, k, ax, ay, az, fe, fw, &
        !$omp& fn, fs, ft, fb, qe, qw, qn, qs, qt, qb)
        DO i = 3, ii-2
            DO j = 3, jj-2
                DO k = 3-nbw, kk-3+ntw
                    ax = ddy(j)*dz(k)
                    ay = ddx(i)*dz(k)
                    az = ddx(i)*ddy(j)

                    fe = ax*(ut(k, j, i) + ut(k+1, j, i))*0.5
                    fw = ax*(ut(k, j, i-1)+ ut(k+1, j, i-1))*0.5
                    fn = ay*(vt(k, j, i) + vt(k+1, j, i))*0.5
                    fs = ay*(vt(k, j-1, i)+ vt(k+1, j-1, i))*0.5
                    ft = az*(wt(k, j, i) + (wt(k+1, j, i) - wt(k, j, i)) &
                        * 0.5*dz(k)/ddz(k+1))
                    fb = az*(wt(k-1, j, i) + (wt(k, j, i) - wt(k-1, j, i)) &
                        * 0.5*dz(k-1)/ddz(k))

                    qe = 0.5*fe*(w(k, j, i) + w(k, j, i+1))
                    qw = 0.5*fw*(w(k, j, i-1) + w(k, j, i))
                    qn = 0.5*fn*(w(k, j, i) + w(k, j+1, i))
                    qs = 0.5*fs*(w(k, j-1, i) + w(k, j, i))
                    qt = 0.5*ft*(w(k, j, i) + w(k+1, j, i))
                    qb = 0.5*fb*(w(k-1, j, i) + w(k, j, i))

                    wo(k, j, i) = -(qe-qw+qn-qs+qt-qb)*rddx(i)*rddy(j)*rdz(k)
                END DO
            END DO
        END DO
        !$omp end do

    END SUBROUTINE tstle4_kon_w


    SUBROUTINE tstle4_diff_swc(kk, jj, ii, uo, vo, wo, u, v, w, ddx, ddy, &
            ddz, nfro, nbac, nrgt, nlft, nbot, ntop)
        !$omp declare target
        ! Subroutine arguments
        INTEGER(intk), INTENT(in) :: kk, jj, ii
        REAL(realk), INTENT(inout) :: uo(kk, jj, ii), vo(kk, jj, ii), &
            wo(kk, jj, ii)
        REAL(realk), INTENT(in) :: u(kk, jj, ii), v(kk, jj, ii), w(kk, jj, ii)
        REAL(realk), INTENT(in) :: ddx(ii), ddy(jj), ddz(kk)
        INTEGER, INTENT(in) :: nfro, nbac, nrgt, nlft, nbot, ntop

        CALL swcle3d(kk, jj, ii, uo, vo, wo, u, v, w, &
            ddx, ddy, ddz, nfro, nbac, nrgt, nlft, nbot, ntop)

    END SUBROUTINE tstle4_diff_swc


    ! SUBROUTINE tstle4_diff_u(kk, jj, ii, uo, u, v, w, g, dx, dy, dz, ddx, &
    !         ddy, ddz, rdx, rdy, rdz, rddx, rddy, rddz, nfro, nbac, nrgt, &
    !         nlft, nbot, ntop)
    !     !$omp declare target

    !     ! Subroutine arguments
    !     INTEGER(intk), INTENT(in) :: kk, jj, ii
    !     REAL(realk), INTENT(inout) :: uo(kk, jj, ii)
    !     REAL(realk), INTENT(in) :: u(kk, jj, ii), v(kk, jj, ii), w(kk, jj, ii)
    !     REAL(realk), INTENT(in) :: g(kk, jj, ii)
    !     REAL(realk), INTENT(in) :: dx(ii), dy(jj), dz(kk)
    !     REAL(realk), INTENT(in) :: ddx(ii), ddy(jj), ddz(kk)
    !     REAL(realk), INTENT(in) :: rdx(ii), rdy(jj), rdz(kk)
    !     REAL(realk), INTENT(in) :: rddx(ii), rddy(jj), rddz(kk)
    !     INTEGER, INTENT(in) :: nfro, nbac, nrgt, nlft, nbot, ntop

    !     CALL tstle4_diff_u_c(INT(kk, c_intk), INT(jj, c_intk), &
    !         INT(ii, c_intk), uo, u, v, w, g, dx, dy, dz, ddx, ddy, ddz, &
    !         rdx, rdy, rdz, rddx, rddy, rddz, INT(nfro, c_intk), &
    !         INT(nbac, c_intk), INT(ilesmodel, c_intk), REAL(gmol, c_realk), &
    !         REAL(rho, c_realk))

    ! END SUBROUTINE tstle4_diff_u


    SUBROUTINE tstle4_diff_v(kk, jj, ii, vo, u, v, w, g, dx, dy, dz, ddx, &
            ddy, ddz, rdx, rdy, rdz, rddx, rddy, rddz, nfro, nbac, nrgt, &
            nlft, nbot, ntop)
        !$omp declare target
        ! Subroutine arguments
        INTEGER(intk), INTENT(in) :: kk, jj, ii
        REAL(realk), INTENT(inout) :: vo(kk, jj, ii)
        REAL(realk), INTENT(in) :: u(kk, jj, ii), v(kk, jj, ii), w(kk, jj, ii)
        REAL(realk), INTENT(in) :: g(kk, jj, ii)
        REAL(realk), INTENT(in) :: dx(ii), dy(jj), dz(kk)
        REAL(realk), INTENT(in) :: ddx(ii), ddy(jj), ddz(kk)
        REAL(realk), INTENT(in) :: rdx(ii), rdy(jj), rdz(kk)
        REAL(realk), INTENT(in) :: rddx(ii), rddy(jj), rddz(kk)
        INTEGER, INTENT(in) :: nfro, nbac, nrgt, nlft, nbot, ntop

        ! Local variables
        INTEGER(intk) :: k, j, i
        INTEGER(intk) :: nrv, nlv
        INTEGER(intk) :: iles
        REAL(realk) :: ax, ay, az
        REAL(realk) :: ge, gw, gn, gs, gt, gb
        REAL(realk) :: qw, qe, qt, qb, qn, qs
        REAL(realk) :: st, qc, fak

        nrv = 0
        nlv = 0

        ! CON = 7
        IF (nlft == 7) nlv = 1

        ! OP1 = 3
        IF (nrgt == 3) nrv = 1
        IF (nlft == 3) nlv = 1

        iles = 1
        IF (ilesmodel == 0) iles = 0

        !$omp do collapse(3) private(i, j, k, ax, ay, az, ge, gw, &
        !$omp& gn, gs, gt, gb, qe, qw, qn, qs, qt, qb, st, qc, fak)
        DO i = 3, ii-2
            DO j = 3-nrv, jj-3+nlv
                DO k = 3, kk-2
                    ax = dy(j)*ddz(k)
                    ay = ddx(i)*ddz(k)
                    az = ddx(i)*dy(j)

                    ge = g(k, j, i)*g(k, j, i+1) &
                        /MAX(g(k, j, i) + g(k, j, i+1), gmol) &
                        + g(k, j+1, i)*g(k, j+1, i+1) &
                        /MAX(g(k, j+1, i) + g(k, j+1, i+1), gmol)
                    gw = g(k, j, i-1)*g(k, j, i) &
                        /MAX(g(k, j, i-1) + g(k, j, i), gmol) &
                        + g(k, j+1, i-1)*g(k, j+1, i) &
                        /MAX(g(k, j+1, i-1) + g(k, j+1, i), gmol)
                    gn = g(k, j+1, i)
                    gs = g(k, j, i)
                    gt = g(k, j, i)*g(k+1, j, i) &
                        /MAX(g(k, j, i) + g(k+1, j, i), gmol) &
                        + g(k, j+1, i)*g(k+1, j+1, i) &
                        /MAX(g(k, j+1, i) + g(k+1, j+1, i), gmol)
                    gb = g(k-1, j, i)*g(k, j, i) &
                        /MAX(g(k-1, j, i) + g(k, j, i), gmol) &
                        + g(k-1, j+1, i)*g(k, j+1, i) &
                        /MAX(g(k-1, j+1, i) + g(k, j+1, i), gmol)

                    qe = -ge*ax*rdx(i) * (v(k, j, i+1) - v(k, j, i))
                    qw = -gw*ax*rdx(i-1) * (v(k, j, i) - v(k, j, i-1))
                    qn = -gn*ay*rddy(j+1) * (v(k, j+1, i) - v(k, j, i))
                    qs = -gs*ay*rddy(j) * (v(k, j, i) - v(k, j-1, i))
                    qt = -gt*az*rdz(k) * (v(k+1, j, i) - v(k, j, i))
                    qb = -gb*az*rdz(k-1) * (v(k, j, i) - v(k-1, j, i))

                    st = ((ge*(u(k, j+1, i) - u(k, j, i))) &
                        - (gw*(u(k, j+1, i-1) - u(k, j, i-1))))*ddz(k) &
                        + ((gn*(v(k, j+1, i) - v(k, j, i))*rddy(j+1)) &
                        - (gs*(v(k, j, i) - v(k, j-1, i))*rddy(j)))*ay &
                        + ((gt*(w(k, j+1, i) - w(k, j, i))) &
                        - (gb*(w(k-1, j+1, i)- w(k-1, j, i))))*ddx(i)
                    qc = st * iles

                    fak = 1.0/rho*rddx(i)*rdy(j)*rddz(k)
                    vo(k, j, i) = vo(k, j, i) - fak*(qe-qw+qn-qs+qt-qb-qc)
                END DO
            END DO
        END DO
        !$omp end do
    END SUBROUTINE tstle4_diff_v


    SUBROUTINE tstle4_diff_w(kk, jj, ii, wo, u, v, w, g, dx, dy, dz, ddx, &
            ddy, ddz, rdx, rdy, rdz, rddx, rddy, rddz, nfro, nbac, nrgt, &
            nlft, nbot, ntop)
        !$omp declare target
        ! Subroutine arguments
        INTEGER(intk), INTENT(in) :: kk, jj, ii
        REAL(realk), INTENT(inout) :: wo(kk, jj, ii)
        REAL(realk), INTENT(in) :: u(kk, jj, ii), v(kk, jj, ii), w(kk, jj, ii)
        REAL(realk), INTENT(in) :: g(kk, jj, ii)
        REAL(realk), INTENT(in) :: dx(ii), dy(jj), dz(kk)
        REAL(realk), INTENT(in) :: ddx(ii), ddy(jj), ddz(kk)
        REAL(realk), INTENT(in) :: rdx(ii), rdy(jj), rdz(kk)
        REAL(realk), INTENT(in) :: rddx(ii), rddy(jj), rddz(kk)
        INTEGER, INTENT(in) :: nfro, nbac, nrgt, nlft, nbot, ntop

        ! Local variables
        INTEGER(intk) :: k, j, i
        INTEGER(intk) :: nbw, ntw
        INTEGER(intk) :: iles
        REAL(realk) :: ax, ay, az
        REAL(realk) :: ge, gw, gn, gs, gt, gb
        REAL(realk) :: qw, qe, qt, qb, qn, qs
        REAL(realk) :: st, qc, fak

        nbw = 0
        ntw = 0

        ! CON = 7
        IF (ntop == 7) ntw = 1

        ! OP1 = 3
        IF (nbot == 3) nbw = 1
        IF (ntop == 3) ntw = 1

        iles = 1
        IF (ilesmodel == 0) iles = 0

        !$omp do collapse(3) private(i, j, k, ax, ay, az, ge, gw, &
        !$omp& gn, gs, gt, gb, qe, qw, qn, qs, qt, qb, st, qc, fak)
        DO i = 3, ii-2
            DO j = 3, jj-2
                DO k = 3-nbw, kk-3+ntw
                    ax = ddy(j)*dz(k)
                    ay = ddx(i)*dz(k)
                    az = ddx(i)*ddy(j)

                    ge = g(k, j, i)*g(k, j, i+1) &
                        /MAX(g(k, j, i) + g(k, j, i+1), gmol) &
                        + g(k+1, j, i)*g(k+1, j, i+1) &
                        /MAX(g(k+1, j, i) + g(k+1, j, i+1), gmol)
                    gw = g(k, j, i-1)*g(k, j, i) &
                        /MAX(g(k, j, i-1) + g(k, j, i), gmol) &
                        + g(k+1, j, i-1)*g(k+1, j, i) &
                        /MAX(g(k+1, j, i-1) + g(k+1, j, i), gmol)
                    gn = g(k, j, i)*g(k, j+1, i) &
                        /MAX(g(k, j, i) + g(k, j+1, i), gmol) &
                        + g(k+1, j, i)*g(k+1, j+1, i) &
                        /MAX(g(k+1, j, i) + g(k+1, j+1, i), gmol)
                    gs = g(k, j-1, i)*g(k, j, i) &
                        /MAX(g(k, j-1, i) + g(k, j, i), gmol) &
                        + g(k+1, j-1, i)*g(k+1, j, i) &
                        /MAX(g(k+1, j-1, i) + g(k+1, j, i), gmol)
                    gt = g(k+1, j, i)
                    gb = g(k, j, i)

                    qe = -ge*ax*rdx(i) * (w(k, j, i+1) - w(k, j, i))
                    qw = -gw*ax*rdx(i-1) * (w(k, j, i) - w(k, j, i-1))
                    qn = -gn*ay*rdy(j) * (w(k, j+1, i) - w(k, j, i))
                    qs = -gs*ay*rdy(j-1) * (w(k, j, i) - w(k, j-1, i))
                    qt = -gt*az*rddz(k+1)* (w(k+1, j, i) - w(k, j, i))
                    qb = -gb*az*rddz(k) * (w(k, j, i) - w(k-1, j, i))

                    st = ((ge*(u(k+1, j, i) - u(k, j, i))) &
                        - (gw*(u(k+1, j, i-1) - u(k, j, i-1))))*ddy(j) &
                        + ((gn*(v(k+1, j, i) - v(k, j, i))) &
                        - (gs*(v(k+1, j-1, i) - v(k, j-1, i))))*ddx(i) &
                        + ((gt*(w(k+1, j, i) - w(k, j, i))*rddz(k+1)) &
                        - (gb*(w(k, j, i) - w(k-1, j, i))*rddz(k)))*az
                    qc = st * iles

                    fak = 1.0/rho*rddx(i)*rddy(j)*rdz(k)
                    wo(k, j, i) = wo(k, j, i) - fak*(qe-qw+qn-qs+qt-qb-qc)
                END DO
            END DO
        END DO
        !$omp end do
    END SUBROUTINE tstle4_diff_w
    SUBROUTINE tstle4_gradp(kk, jj, ii, uo, vo, wo, p, dx, dy, dz, &
            nfro, nbac, nrgt, nlft, nbot, ntop, igrid)
        !$omp declare target
        ! Subroutine arguments
        INTEGER(intk), INTENT(in) :: kk, jj, ii
        REAL(realk), INTENT(inout) :: uo(kk, jj, ii), vo(kk, jj, ii), &
            wo(kk, jj, ii)
        REAL(realk), INTENT(in) :: p(kk, jj, ii)
        REAL(realk), INTENT(in) :: dx(ii), dy(jj), dz(kk)
        INTEGER, INTENT(in) :: nfro, nbac, nrgt, nlft, nbot, ntop
        INTEGER, INTENT(in) :: igrid

        ! Local variables
        INTEGER(intk) :: k, j, i
        INTEGER(intk) :: nbu, nfu, nrv, nbw, ntw, nlv
        INTEGER(intk) :: gradpflag
        REAL(realk) :: gpx, gpy, gpz

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

        CALL get_gradpxflag(gradpflag, igrid)
        gpx = gradp(1)*gradpflag
        gpy = gradp(2)*gradpflag
        gpz = gradp(3)*gradpflag

        !$omp do collapse(3) private(i, j, k)
        DO i = 3-nfu, ii-3+nbu
            DO j = 3, jj-2
                DO k = 3, kk-2
                    uo(k, j, i) = uo(k, j, i) - 1.0/(rho*dx(i)) &
                        *(p(k, j, i+1) - p(k, j, i) + gpx*dx(i))
                END DO
            END DO
        END DO
        !$omp end do

        !$omp do collapse(3) private(i, j, k)
        DO i = 3, ii-2
            DO j = 3-nrv, jj-3+nlv
                DO k = 3, kk-2
                    vo(k, j, i) = vo(k, j, i) - 1.0/(rho*dy(j)) &
                        *(p(k, j+1, i) - p(k, j, i) + gpy*dy(j))
                END DO
            END DO
        END DO
        !$omp end do

        !$omp do collapse(3) private(i, j, k)
        DO i = 3, ii-2
            DO j = 3, jj-2
                DO k = 3-nbw, kk-3+ntw
                    wo(k, j, i) = wo(k, j, i) - 1.0/(rho*dz(k)) &
                        *(p(k+1, j, i) - p(k, j, i) + gpz*dz(k))
                END DO
            END DO
        END DO
        !$omp end do
    END SUBROUTINE tstle4_gradp


    SUBROUTINE tstle4_par(kk, jj, ii, uo, vo, wo, u, v, w, ut, vt, wt, &
            dx, dy, dz, ddx, ddy, ddz, rdx, rdy, rdz, rddx, rddy, rddz, &
            wcu, wcv, wcw, nfro, nbac, nrgt, nlft, nbot, ntop)
        !$omp declare target
        ! Subroutine arguments
        INTEGER(intk), INTENT(in) :: kk, jj, ii
        REAL(realk), INTENT(inout) :: uo(kk, jj, ii), vo(kk, jj, ii), &
            wo(kk, jj, ii)
        REAL(realk), INTENT(in) :: u(kk, jj, ii), v(kk, jj, ii), w(kk, jj, ii)
        REAL(realk), INTENT(in) :: ut(kk, jj, ii), vt(kk, jj, ii), &
            wt(kk, jj, ii)
        REAL(realk), INTENT(in) :: dx(ii), dy(jj), dz(kk)
        REAL(realk), INTENT(in) :: ddx(ii), ddy(jj), ddz(kk)
        REAL(realk), INTENT(in) :: rdx(ii), rdy(jj), rdz(kk)
        REAL(realk), INTENT(in) :: rddx(ii), rddy(jj), rddz(kk)
        REAL(realk), INTENT(inout) :: wcu(kk, jj, ii), wcv(kk, jj, ii), &
            wcw(kk, jj, ii)
        INTEGER, INTENT(in) :: nfro, nbac, nrgt, nlft, nbot, ntop

        ! Local variables
        INTEGER(intk) :: k, j, i
        REAL(realk) :: fkdtu, fkdtv, fkdtw
        REAL(realk) :: qkubadd, qkusadd, qkvbadd
        REAL(realk) :: qkvwadd, qkwsadd, qkwwadd

        REAL(realk) :: qkut, qkub, qkun, qkus, qkvw, qkve, qkvt, qkvb, &
            qkww, qkwe, qkwn, qkws, &
            fut, fub, fun, fus, auy, auz, &
            fvw, fve, fvt, fvb, avx, avz, &
            fww, fwe, fwn, fws, awx, awy
        REAL(realk) :: dxi, ddxi, dyj, ddyj, dzk, ddzk, rdzk, rddzk
        REAL(realk), PARAMETER :: wkon = 1.0

        ! Upwind in vorletzter schicht bei PAR-randbedingung
        IF (nfro == 8) THEN
            i = 4
            !$omp do private(j, k, dyj, ddyj, fkdtu, fkdtv, fkdtw, &
            !$omp& dzk, rdzk, ddzk, rddzk, avx, awx, fvw, qkvwadd, fww, &
            !$omp& qkwwadd)
            DO j = 3, jj-2
                dyj = dy(j)
                ddyj = ddy(j)
                fkdtu = -1.0*rddy(j)*rdx(i)*wkon
                fkdtv = -1.0*rddx(i)*rdy(j)*wkon
                fkdtw = -1.0*rddx(i)*rddy(j)*wkon

                DO k = 3, kk-2
                    dzk = dz(k)
                    rdzk = rdz(k)
                    ddzk = ddz(k)
                    rddzk = rddz(k)
                    avx = dyj*ddzk
                    awx = ddyj*dzk
                    fvw = avx*(ut(k, j, i-1) + ut(k, j+1, i-1))*0.5
                    qkvwadd = 0.5*(fvw-ABS(fvw)) &
                        *(0.5*v(k, j, i)-1.5*v(k, j, i-1)+v(k, j, i-2))*0.5*0.5

                    vo(k, j, i) = vo(k, j, i) + fkdtv * rddzk * (-qkvwadd)
                    vo(k, j, i-1) = vo(k, j, i-1) + fkdtv * rddzk * (+qkvwadd)

                    fww = awx*(ut(k, j, i-1)+ ut(k+1, j, i-1))*0.5
                    qkwwadd = 0.5 *(fww-ABS(fww)) &
                        *(0.5*w(k, j, i)-1.5*w(k, j, i-1)+w(k, j, i-2))*0.5*0.5

                    wo(k, j, i) = wo(k, j, i) + fkdtw * rdzk * (-qkwwadd)
                    wo(k, j, i-1) = wo(k, j, i-1) + fkdtw * rdzk * (+qkwwadd)
                END DO
            END DO
            !$omp end do
            !$omp barrier
        END IF

        IF (nbac == 8) THEN
            i = ii-2
            !$omp do private(j, k, dyj, ddyj, fkdtu, fkdtv, fkdtw, &
            !$omp& dzk, rdzk, ddzk, rddzk, avx, awx, fvw, qkvwadd, fww, &
            !$omp& qkwwadd)
            DO j = 3, jj-2
                dyj = dy(j)
                ddyj = ddy(j)
                fkdtu = -1.0*rddy(j)*rdx(i)*wkon
                fkdtv = -1.0*rddx(i)*rdy(j)*wkon
                fkdtw = -1.0*rddx(i)*rddy(j)*wkon

                DO k = 3, kk-2
                    dzk = dz(k)
                    rdzk = rdz(k)
                    ddzk = ddz(k)
                    rddzk = rddz(k)
                    avx = dyj*ddzk
                    awx = ddyj*dzk
                    fvw = avx*(ut(k, j, i-1) + ut(k, j+1, i-1))*0.5

                    qkvwadd = 0.5*(fvw+ABS(fvw)) &
                        *(0.5*v(k, j, i-1)-1.5*v(k, j, i)+v(k, j, i+1))*0.5*0.5

                    vo(k, j, i) = vo(k, j, i) + fkdtv * rddzk * (-qkvwadd)
                    vo(k, j, i-1) = vo(k, j, i-1) + fkdtv * rddzk * (+qkvwadd)

                    fww = awx*(ut(k, j, i-1)+ ut(k+1, j, i-1))*0.5
                    qkwwadd = 0.5*(fww+ABS(fww)) &
                        *(0.5*w(k, j, i-1)-1.5*w(k, j, i)+w(k, j, i+1))*0.5*0.5

                    wo(k, j, i) = wo(k, j, i) + fkdtw * rdzk * (-qkwwadd)
                    wo(k, j, i-1) = wo(k, j, i-1) + fkdtw * rdzk * (+qkwwadd)
                END DO
            END DO
            !$omp end do
            !$omp barrier
        END IF

        IF (nrgt == 8) THEN
            !$omp do private(i, j, k, dxi, ddxi, fkdtu, fkdtv, &
            !$omp& fkdtw, dzk, rdzk, ddzk, rddzk, auy, awy, fus, qkusadd, &
            !$omp& fws, qkwsadd)
            DO i = 3, ii-2
                j = 4
                dxi = dx(i)
                ddxi = ddx(i)
                fkdtu = -1.0*rddy(j)*rdx(i)*wkon
                fkdtv = -1.0*rddx(i)*rdy(j)*wkon
                fkdtw = -1.0*rddx(i)*rddy(j)*wkon

                DO k = 3, kk-2
                    dzk = dz(k)
                    rdzk = rdz(k)
                    ddzk = ddz(k)
                    rddzk = rddz(k)
                    auy = dxi*ddzk
                    awy = ddxi*dzk

                    fus = auy*(vt(k, j-1, i) + vt(k, j-1, i+1))*0.5
                    qkusadd = 0.5 *(fus-ABS(fus)) &
                        *(0.5*u(k, j, i)-1.5*u(k, j-1, i)+u(k, j-2, i))*0.5*0.5

                    uo(k, j, i) = uo(k, j, i) + fkdtu * rddzk * (-qkusadd)
                    uo(k, j-1, i) = uo(k, j-1, i) + fkdtu * rddzk * (+qkusadd)

                    fws = awy*(vt(k, j-1, i)+ vt(k+1, j-1, i))*0.5
                    qkwsadd = 0.5 *(fws-ABS(fws)) &
                        *(0.5*w(k, j, i)-1.5*w(k, j-1, i)+w(k, j-2, i))*0.5*0.5

                    wo(k, j, i) = wo(k, j, i) + fkdtw * rdzk * (-qkwsadd)
                    wo(k, j-1, i) = wo(k, j-1, i) + fkdtw * rdzk * (+qkwsadd)
                END DO
            END DO
            !$omp end do
            !$omp barrier
        END IF

        IF (nlft == 8) THEN
            !$omp do private(i, j, k, dxi, ddxi, fkdtu, fkdtv, &
            !$omp& fkdtw, dzk, rdzk, ddzk, rddzk, auy, awy, fus, qkusadd, &
            !$omp& fws, qkwsadd)
            DO i = 3, ii-2
                j = jj-2
                dxi = dx(i)
                ddxi = ddx(i)
                fkdtu = -1.0*rddy(j)*rdx(i)*wkon
                fkdtv = -1.0*rddx(i)*rdy(j)*wkon
                fkdtw = -1.0*rddx(i)*rddy(j)*wkon

                DO k = 3, kk-2
                    dzk = dz(k)
                    rdzk = rdz(k)
                    ddzk = ddz(k)
                    rddzk = rddz(k)
                    auy = dxi*ddzk
                    awy = ddxi*dzk

                    fus = auy*(vt(k, j-1, i) + vt(k, j-1, i+1))*0.5
                    qkusadd = 0.5*(fus+ABS(fus)) &
                        *(0.5*u(k, j-1, i)-1.5*u(k, j, i)+u(k, j+1, i))*0.5*0.5

                    uo(k, j, i) = uo(k, j, i) + fkdtu * rddzk * (-qkusadd)
                    uo(k, j-1, i) = uo(k, j-1, i) + fkdtu * rddzk * (+qkusadd)
                    fws = awy*(vt(k, j-1, i) + vt(k+1, j-1, i))*0.5
                    qkwsadd = 0.5*(fws+ABS(fws)) &
                        *(0.5*w(k, j-1, i)-1.5*w(k, j, i)+w(k, j+1, i))*0.5*0.5

                    wo(k, j, i) = wo(k, j, i) + fkdtw * rdzk * (-qkwsadd)
                    wo(k, j-1, i) = wo(k, j-1, i) + fkdtw * rdzk * (+qkwsadd)
                END DO
            END DO
            !$omp end do
            !$omp barrier
        END IF

        IF (nbot == 8) THEN
            !$omp do private(i, j, k, dxi, ddxi, fkdtu, fkdtv, &
            !$omp& fkdtw, dyj, ddyj, auz, avz, rdzk, rddzk, fub, qkubadd, &
            !$omp& fvb, qkvbadd)
            DO i = 3, ii-2
                dxi = dx(i)
                ddxi = ddx(i)

                DO j = 3, jj-2
                    k = 4

                    fkdtu = -1.0*rddy(j)*rdx(i)*wkon
                    fkdtv = -1.0*rddx(i)*rdy(j)*wkon
                    fkdtw = -1.0*rddx(i)*rddy(j)*wkon

                    dyj = dy(j)
                    ddyj = ddy(j)
                    auz = dxi*ddyj
                    avz = ddxi*dyj
                    rdzk = rdz(k)
                    rddzk = rddz(k)

                    fub = auz*(wt(k-1, j, i)+ wt(k-1, j, i+1))*0.5
                    qkubadd = 0.5*(fub-ABS(fub)) &
                        *(0.5*u(k, j, i)-1.5*u(k-1, j, i)+u(k-2, j, i))*0.5*0.5
                    uo(k, j, i) = uo(k, j, i) + fkdtu * rddzk * (-qkubadd)
                    uo(k-1, j, i) = uo(k-1, j, i) + fkdtu * rddzk * (+qkubadd)

                    fvb = avz*(wt(k-1, j, i) + wt(k-1, j+1, i))*0.5
                    qkvbadd = 0.5*(fvb-ABS(fvb)) &
                        *(0.5*v(k, j, i)-1.5*v(k-1, j, i)+v(k-2, j, i))*0.5*0.5

                    vo(k, j, i) = vo(k, j, i) + fkdtv * rddzk * (-qkvbadd)
                    vo(k-1, j, i) = vo(k-1, j, i) + fkdtv * rddzk * (+qkvbadd)
                END DO
            END DO
            !$omp end do
            !$omp barrier
        END IF

        IF (ntop == 8) THEN
            !$omp do private(i, j, k, dxi, ddxi, fkdtu, fkdtv, &
            !$omp& fkdtw, dyj, ddyj, auz, avz, rdzk, rddzk, fub, qkubadd, &
            !$omp& fvb, qkvbadd)
            DO i = 3, ii-2
                dxi = dx(i)
                ddxi = ddx(i)

                DO j = 3, jj-2
                    k = kk-2

                    fkdtu = -1.0*rddy(j)*rdx(i)*wkon
                    fkdtv = -1.0*rddx(i)*rdy(j)*wkon
                    fkdtw = -1.0*rddx(i)*rddy(j)*wkon

                    dyj = dy(j)
                    ddyj = ddy(j)
                    auz = dxi*ddyj
                    avz = ddxi* dyj

                    rdzk = rdz(k)
                    rddzk = rddz(k)

                    fub = auz*(wt(k-1, j, i) + wt(k-1, j, i+1))*0.5
                    qkubadd = 0.5 *(fub+ABS(fub)) &
                        *(0.5*u(k-1, j, i)-1.5*u(k, j, i)+u(k+1, j, i))*0.5*0.5
                    uo(k, j, i) = uo(k, j, i) + fkdtu * rddzk * (-qkubadd)
                    uo(k-1, j, i) = uo(k-1, j, i) + fkdtu * rddzk * (+qkubadd)

                    fvb = avz*(wt(k-1, j, i)+ wt(k-1, j+1, i))*0.5
                    qkvbadd = 0.5 *(fvb+ABS(fvb)) &
                        *(0.5*v(k-1, j, i)-1.5*v(k, j, i)+v(k+1, j, i))*0.5*0.5

                    vo(k, j, i) = vo(k, j, i) + fkdtv * rddzk * (-qkvbadd)
                    vo(k-1, j, i) = vo(k-1, j, i) + fkdtv * rddzk * (+qkvbadd)
                END DO
            END DO
            !$omp end do
            !$omp barrier
        END IF

        ! PAR-RB Impulserhaltend BACK
        IF (nbac == 8) THEN
            ! W-Impulszelle
            ! STANDART-BERECHNUNG DES KONVEKTIVEN FLUSSES
            i = ii-2
            !$omp do private(j, k, ddyj, dzk, awx, fwe, qkwe)
            DO j = 3, jj-2
                ddyj = ddy(j)
                DO k = 2, kk-2
                    dzk = dz(k)
                    awx = ddyj*dzk
                    fwe = awx*(ut(k, j, i) + ut(k+1, j, i))*0.5
                    qkwe = 0.5*fwe*(w(k, j, i) + w(k, j, i+1))
                    wcw(k, j, i+1) = qkwe
                END DO
            END DO
            !$omp end do
            !$omp barrier

            ! NEUE BERECHNUNG DES KONVEKTIVEN FLUSSES
            ! FUER JEDE GROBGITTERGESCHWINDIGKEITT
            i = ii-2
            !$omp do private(j, k, ddyj, dzk, awx, fwe, qkwe)
            DO j = 3, jj-2, 2
                ddyj = ddy(j)
                k = 2
                dzk = dz(k)
                awx = ddyj*dzk
                fwe = awx*(2*ut(k, j, i) + 2*ut(k, j+1, i) &
                    + ut(k+1, j, i) + ut(k+1, j+1, i) &
                    + ut(k+2, j, i) + ut(k+2, j+1, i))*0.125
                qkwe = 0.5*fwe*(w(k, j, i) + w(k, j, i+1))
                wcw(k, j, i) = qkwe
                qkwe = 0.5*fwe*(w(k, j+1, i) + w(k, j+1, i+1))
                wcw(k, j+1, i) = qkwe
                DO k = 4, kk-4, 2
                    dzk = dz(k)
                    awx = ddyj*dzk
                    fwe = awx*(ut(k-1, j, i) + ut(k-1, j+1, i) &
                        + ut(k, j, i) + ut(k, j+1, i) &
                        + ut(k+1, j, i) + ut(k+1, j+1, i) &
                        + ut(k+2, j, i)+ ut(k+2, j+1, i))*0.125
                    qkwe = 0.5*fwe*(w(k, j, i) + w(k, j, i+1))
                    wcw(k, j, i) = qkwe
                    qkwe = 0.5*fwe*(w(k, j+1, i) + w(k, j+1, i+1))
                    wcw(k, j+1, i) = qkwe
                END DO
                k = kk-2
                dzk = dz(k)
                awx = ddyj*dzk
                fwe = awx*(ut(k-1, j, i) + ut(k-1, j+1, i) &
                    + ut(k, j, i) + ut(k, j+1, i)&
                    + 2*ut(k+1, j, i) + 2*ut(k+1, j+1, i)) *0.125
                qkwe = 0.5*fwe*(w(k, j, i) + w(k, j, i+1))
                wcw(k, j, i) = qkwe
                qkwe = 0.5*fwe*(w(k, j+1, i) + w(k, j+1, i+1))
                wcw(k, j+1, i) = qkwe
            END DO
            !$omp end do
            !$omp barrier

            ! VERTEILUNG
            i = ii-2
            !$omp do collapse(2) private(j, k)
            DO j = 3, jj-2
                DO k = 2, kk-4, 2
                    wcw(k+1, j, i) = 0.5*(wcw(k, j, i) + wcw(k+2, j, i))
                END DO
            END DO
            !$omp end do
            !$omp barrier

            ! AUF WO SCHREIBEN
            !$omp do private(j, k, fkdtw, rdzk)
            DO j = 3, jj-2
                fkdtw = -1.0*rddx(i)*rddy(j)*wkon
                DO k = 2, kk-2
                    rdzk = rdz(k)
                    wo(k, j, i) = wo(k, j, i) &
                        + fkdtw*rdzk*(-wcw(k, j, i+1) + wcw(k, j, i))
                END DO
            END DO
            !$omp end do
            !$omp barrier

            ! V-Impulszelle
            ! STANDARD-BERECHNUNG DES KONVEKTIVEN FLUSSES
            i = ii-2
            !$omp do private(j, k, dyj, ddzk, avx, fve, qkve)
            DO j = 2, jj-2
                dyj = dy(j)
                DO k = 3, kk-2
                    ddzk = ddz(k)
                    avx = ddzk*dyj
                    fve = avx *(ut(k, j, i) + ut(k, j+1, i))*0.5
                    qkve = 0.5*fve*(v(k, j, i) + v(k, j, i+1))
                    wcv(k, j, i+1) = qkve
                END DO
            END DO
            !$omp end do
            !$omp barrier

            ! NEUE BERECHNUNG DES KONVEKTIVEN FLUSSES
            ! FUER JEDE GROBGITTERGESCHWINDIGKEIT
            i = ii-2

            ! YM-RAND
            j = 2
            dyj = dy(j)
            !$omp do private(k, ddzk, avx, fve, qkve)
            DO k = 3, kk-2, 2
                ddzk = ddz(k)
                avx = ddzk*dyj
                fve = avx*(2*ut(k, j, i) + 2*ut(k+1, j, i) &
                    + ut(k, j+1, i) + ut(k+1, j+1, i) &
                    + ut(k, j+2, i) + ut(k+1, j+2, i))*0.125
                qkve = 0.5*fve*(v(k, j, i) + v(k, j, i+1))
                wcv(k, j, i) = qkve
                qkve = 0.5*fve*(v(k+1, j, i) + v(k+1, j, i+1))
                wcv(k+1, j, i) = qkve
            END DO
            !$omp end do
            !$omp barrier

            ! IM GEBIET
            !$omp do private(j, k, dyj, ddzk, avx, fve, qkve)
            DO j = 4, jj-4, 2
                dyj = dy(j)
                DO k = 3, kk-2, 2
                    ddzk = ddz(k)
                    avx = ddzk*dyj
                    fve = avx*(ut(k, j-1, i) + ut(k+1, j-1, i) &
                        + ut(k, j, i) + ut(k+1, j, i) &
                        + ut(k, j+1, i) + ut(k+1, j+1, i) &
                        + ut(k, j+2, i) + ut(k+1, j+2, i))*0.125
                    qkve = 0.5*fve*(v(k, j, i) + v(k, j, i+1))
                    wcv(k, j, i) = qkve
                    qkve = 0.5*fve*(v(k+1, j, i) + v(k+1, j, i+1))
                    wcv(k+1, j, i) = qkve
                END DO
            END DO
            !$omp end do
            !$omp barrier

            ! YP-RAND
            j = jj-2
            dyj = dy(j)
            !$omp do private(k, ddzk, avx, fve, qkve)
            DO k = 3, kk-2, 2
                ddzk = ddz(k)
                avx = ddzk*dyj
                fve = avx*(ut(k, j-1, i) + ut(k+1, j-1, i) &
                    + ut(k, j, i) + ut(k+1, j, i) &
                    + 2*ut(k, j+1, i) + 2*ut(k+1, j+1, i))*0.125
                qkve = 0.5*fve*(v(k, j, i) + v(k, j, i+1))
                wcv(k, j, i) = qkve
                qkve = 0.5*fve*(v(k+1, j, i) + v(k+1, j, i+1))
                wcv(k+1, j, i) = qkve
            END DO
            !$omp end do
            !$omp barrier

            ! VERTEILUNG
            i = ii-2
            !$omp do collapse(2) private(j, k)
            DO j = 2, jj-4, 2
                DO k = 3, kk-2
                    wcv(k, j+1, i) = 0.5*(wcv(k, j, i) + wcv(k, j+2, i))
                END DO
            END DO
            !$omp end do
            !$omp barrier

            ! AUF VO SCHREIBEN
            !$omp do private(j, k, fkdtv, rddzk)
            DO j = 2, jj-2
                fkdtv = -1.0*rddx(i)* rdy(j)*wkon
                DO k = 3, kk-2
                    rddzk = rddz(k)
                    vo(k, j, i) = vo(k, j, i) &
                        + fkdtv*rddzk*(wcv(k, j, i) - wcv(k, j, i+1))
                END DO
            END DO
            !$omp end do
            !$omp barrier
        END IF

        ! PAR-RB Impulserhaltend FRONT
        IF (nfro == 8) THEN
            ! W-Impulszelle
            ! STANDART-BERECHNUNG DES KONVEKTIVEN FLUSSES
            i = 3
            !$omp do private(j, k, ddyj, dzk, awx, fww, qkww)
            DO j = 3, jj-2
                ddyj = ddy(j)
                DO k = 2, kk-2
                    dzk = dz(k)
                    awx = ddyj*dzk
                    fww = awx*(ut(k, j, i-1) + ut(k+1, j, i-1))*0.5
                    qkww = 0.5*fww*(w(k, j, i) + w(k, j, i-1))
                    wcw(k, j, i-1) = qkww
                END DO
            END DO
            !$omp end do
            !$omp barrier

            ! NEUE BERECHNUNG DES KONVEKTIVEN FLUSSES
            ! FUER JEDE GROBGITTERGESCHWINDIGKEIT
            i = 3
            !$omp do private(j, k, ddyj, dzk, awx, fww, qkww)
            DO j = 3, jj-2, 2
                ddyj = ddy(j)
                k = 2
                dzk = dz(k)
                awx = ddyj*dzk
                fww = awx*(2*ut(k, j, i-1) + 2*ut(k, j+1, i-1) &
                    + ut(k+1, j, i-1) + ut(k+1, j+1, i-1) &
                    + ut(k+2, j, i-1) + ut(k+2, j+1, i-1))*0.125
                qkww = 0.5*fww*(w(k, j, i) + w(k, j, i-1))
                wcw(k, j, i) = qkww
                qkww = 0.5*fww*(w(k, j+1, i) + w(k, j+1, i-1))
                wcw(k, j+1, i) = qkww
                DO k = 4, kk-4, 2
                    dzk = dz(k)
                    awx = ddyj*dzk
                    fww = awx*(ut(k-1, j, i-1) + ut(k-1, j+1, i-1) &
                        + ut(k, j, i-1) + ut(k, j+1, i-1) &
                        + ut(k+1, j, i-1) + ut(k+1, j+1, i-1) &
                        + ut(k+2, j, i-1) + ut(k+2, j+1, i-1))*0.125
                    qkww = 0.5*fww*(w(k, j, i) + w(k, j, i-1))
                    wcw(k, j, i) = qkww
                    qkww = 0.5*fww*(w(k, j+1, i) + w(k, j+1, i-1))
                    wcw(k, j+1, i) = qkww
                END DO
                k = kk-2
                dzk = dz(k)
                awx = ddyj*dzk
                fww = awx*(ut(k-1, j, i-1) + ut(k-1, j+1, i-1) &
                    + ut(k, j, i-1) + ut(k, j+1, i-1) &
                    + 2*ut(k+1, j, i-1) + 2*ut(k+1, j+1, i-1))*0.125
                qkww = 0.5*fww*(w(k, j, i) + w(k, j, i-1))
                wcw(k, j, i) = qkww
                qkww = 0.5*fww*(w(k, j+1, i) + w(k, j+1, i-1))
                wcw(k, j+1, i) = qkww
            END DO
            !$omp end do
            !$omp barrier

            ! VERTEILUNG
            i = 3
            !$omp do collapse(2) private(j, k)
            DO j = 3, jj-2
                DO k = 2, kk-4, 2
                    wcw(k+1, j, i) = 0.5*(wcw(k, j, i) + wcw(k+2, j, i))
                END DO
            END DO
            !$omp end do
            !$omp barrier

            ! AUF WO SCHREIBEN
            !$omp do private(j, k, fkdtw, rdzk)
            DO j = 3, jj-2
                fkdtw = -1.0*rddx(i)*rddy(j)*wkon
                DO k = 2, kk-2
                    rdzk = rdz(k)
                    wo(k, j, i) = wo(k, j, i) &
                        + fkdtw*rdzk*(wcw(k, j, i-1) - wcw(k, j, i))
                END DO
            END DO
            !$omp end do
            !$omp barrier

            ! V-Impulszelle
            ! STANDARD-BERECHNUNG DES KONVEKTIVEN FLUSSES
            i = 3
            !$omp do private(j, k, dyj, ddzk, avx, fvw, qkvw)
            DO j = 2, jj-2
                dyj = dy(j)
                DO k = 3, kk-2
                    ddzk = ddz(k)
                    avx = ddzk*dyj
                    fvw = avx *(ut(k, j, i-1) + ut(k, j+1, i-1))*0.5
                    qkvw = 0.5*fvw*(v(k, j, i) + v(k, j, i-1))
                    wcv(k, j, i-1) = qkvw
                END DO
            END DO
            !$omp end do
            !$omp barrier

            ! NEUE BERECHNUNG DES KONVEKTIVEN FLUSSES
            ! FUER JEDE GROBGITTERGESCHWINDIGKEIT
            i = 3

            ! YM-RAND
            j = 2
            dyj = dy(j)
            !$omp do private(k, ddzk, avx, fvw, qkvw)
            DO k = 3, kk-2, 2
                ddzk = ddz(k)
                avx = ddzk*dyj
                fvw = avx*(2*ut(k, j, i-1) + 2*ut(k+1, j, i-1) &
                    + ut(k, j+1, i-1) + ut(k+1, j+1, i-1) &
                    + ut(k, j+2, i-1) + ut(k+1, j+2, i-1))*0.125
                qkvw = 0.5*fvw*(v(k, j, i) + v(k, j, i-1))
                wcv(k, j, i) = qkvw
                qkvw = 0.5*fvw*(v(k+1, j, i) + v(k+1, j, i-1))
                wcv(k+1, j, i) = qkvw
            END DO
            !$omp end do
            !$omp barrier

            ! IM GEBIET
            !$omp do private(j, k, dyj, ddzk, avx, fvw, qkvw)
            DO j = 4, jj-4, 2
                dyj = dy(j)
                DO k = 3, kk-2, 2
                    ddzk = ddz(k)
                    avx = ddzk*dyj
                    fvw = avx*(ut(k, j-1, i-1) + ut(k+1, j-1, i-1) &
                        + ut(k, j, i-1) + ut(k+1, j, i-1) &
                        + ut(k, j+1, i-1) + ut(k+1, j+1, i-1) &
                        + ut(k, j+2, i-1) + ut(k+1, j+2, i-1))*0.125
                    qkvw = 0.5*fvw*(v(k, j, i) + v(k, j, i-1))
                    wcv(k, j, i) = qkvw
                    qkvw = 0.5*fvw*(v(k+1, j, i) + v(k+1, j, i-1))
                    wcv(k+1, j, i) = qkvw
                END DO
            END DO
            !$omp end do
            !$omp barrier

            ! YP-RAND
            j = jj-2
            dyj = dy(j)
            !$omp do private(k, ddzk, avx, fvw, qkvw)
            DO k = 3, kk-2, 2
                ddzk = ddz(k)
                avx = ddzk*dyj
                fvw = avx*(ut(k, j-1, i-1) + ut(k+1, j-1, i-1) &
                    + ut(k, j, i-1) + ut(k+1, j, i-1) &
                    + 2*ut(k, j+1, i-1) + 2*ut(k+1, j+1, i-1))*0.125
                qkvw = 0.5*fvw*(v(k, j, i) + v(k, j, i-1))
                wcv(k, j, i) = qkvw
                qkvw = 0.5*fvw*(v(k+1, j, i) + v(k+1, j, i-1))
                wcv(k+1, j, i) = qkvw
            END DO
            !$omp end do
            !$omp barrier

            ! VERTEILUNG
            i = 3
            !$omp do collapse(2) private(j, k)
            DO j = 2, jj-4, 2
                DO k = 3, kk-2
                    wcv(k, j+1, i) = 0.5*(wcv(k, j, i) + wcv(k, j+2, i))
                END DO
            END DO
            !$omp end do
            !$omp barrier

            ! AUF VO SCHREIBEN
            !$omp do private(j, k, fkdtv, rddzk)
            DO j = 2, jj-2
                fkdtv = -1.0*rddx(i)*rdy(j)*wkon
                DO k = 3, kk-2
                    rddzk = rddz(k)
                    vo(k, j, i) = vo(k, j, i) &
                        + fkdtv*rddzk*(-wcv(k, j, i) + wcv(k, j, i-1))
                END DO
            END DO
            !$omp end do
            !$omp barrier
        END IF

        ! PAR-RB Impulserhaltend TOP
        IF (ntop == 8) THEN

            ! U-Impulszelle
            ! STANDARD-BERECHNUNG DES KONVEKTIVEN FLUSSES
            k = kk-2
            !$omp do private(i, j, dxi, ddyj, auz, fut, qkut)
            DO i = 2, ii-2
                dxi = dx(i)
                DO j = 3, jj-2
                    ddyj = ddy(j)
                    auz = dxi*ddyj
                    fut = auz*(wt(k, j, i) + wt(k, j, i+1))*0.5
                    qkut = 0.5*fut*(u(k, j, i) + u(k+1, j, i))
                    wcu(k+1, j, i) = qkut
                END DO
            END DO
            !$omp end do
            !$omp barrier

            ! NEUE BERECHNUNG DES KONVEKTIVEN FLUSSES
            ! FUER JEDE GROBGITTERGESCHWINDIGKEIT
            k = kk-2

            ! XM-RAND
            i = 2
            dxi = dx(i)
            !$omp do private(j, ddyj, auz, fut, qkut)
            DO j = 3, jj-2, 2
                ddyj = ddy(j)
                auz = dxi*ddyj
                fut = auz*(2*wt(k, j, i) + 2*wt(k, j+1, i) &
                    + wt(k, j, i+1) + wt(k, j+1, i+1) &
                    + wt(k, j, i+2) + wt(k, j+1, i+2))*0.125
                qkut = 0.5*fut*(u(k, j, i) + u(k+1, j, i))
                wcu(k, j, i) = qkut
                qkut = 0.5*fut*(u(k, j+1, i) + u(k+1, j+1, i))
                wcu(k, j+1, i) = qkut
            END DO
            !$omp end do
            !$omp barrier

            ! IM GEBIET
            !$omp do private(i, j, dxi, ddyj, auz, fut, qkut)
            DO i = 4, ii-4, 2
                dxi = dx(i)
                DO j = 3, jj-2, 2
                    ddyj = ddy(j)
                    auz = dxi*ddyj
                    fut = auz*(wt(k, j, i-1) + wt(k, j+1, i-1) &
                        + wt(k, j, i) + wt(k, j+1, i) &
                        + wt(k, j, i+1) + wt(k, j+1, i+1) &
                        + wt(k, j, i+2) + wt(k, j+1, i+2))*0.125
                    qkut = 0.5*fut*(u(k, j, i) + u(k+1, j, i))
                    wcu(k, j, i) = qkut
                    qkut = 0.5*fut*(u(k, j+1, i) + u(k+1, j+1, i))
                    wcu(k, j+1, i) = qkut
                END DO
            END DO
            !$omp end do
            !$omp barrier

            ! XP-RAND
            i = ii-2
            dxi = dx(i)
            !$omp do private(j, ddyj, auz, fut, qkut)
            DO j = 3, jj-2, 2
                ddyj = ddy(j)
                auz = dxi*ddyj
                fut = auz*(wt(k, j, i-1) + wt(k, j+1, i-1) &
                    + wt(k, j, i) + wt(k, j+1, i) &
                    + 2*wt(k, j, i+1) + 2*wt(k, j+1, i+1))*0.125
                qkut = 0.5*fut*(u(k, j, i) + u(k+1, j, i))
                wcu(k, j, i) = qkut
                qkut = 0.5*fut*(u(k, j+1, i) + u(k+1, j+1, i))
                wcu(k, j+1, i) = qkut
            END DO
            !$omp end do
            !$omp barrier

            ! VERTEILUNG
            k = kk-2
            !$omp do collapse(1) private(i, j)
            DO i = 2, ii-4, 2
                DO j = 3, jj-2
                    wcu(k, j, i+1) = 0.5*(wcu(k, j, i) + wcu(k, j, i+2))
                END DO
            END DO
            !$omp end do
            !$omp barrier

            ! AUF U0 SCHREIBEN
            rddzk = rddz(k)
            !$omp do collapse(2) private(i, j, fkdtu)
            DO i = 2, ii-2
                DO j = 3, jj-2
                    fkdtu = -1.0*rddy(j)*rdx(i)*wkon
                    uo(k, j, i) = uo(k, j, i) &
                        + fkdtu*rddzk*(wcu(k, j, i) - wcu(k+1, j, i))
                END DO
            END DO
            !$omp end do
            !$omp barrier

            ! V-Impulszelle
            ! STANDARD-BERECHNUNG DES KONVEKTIVEN FLUSSES
            k = kk-2
            !$omp do private(i, j, ddxi, dyj, avz, fvt, qkvt)
            DO i = 3, ii-2
                ddxi = ddx(i)
                DO j = 2, jj-2
                    dyj = dy(j)
                    avz = ddxi*dyj
                    fvt = avz*(wt(k, j, i) + wt(k, j+1, i))*0.5
                    qkvt = 0.5*fvt*(v(k, j, i) + v(k+1, j, i))
                    wcv(k+1, j, i) = qkvt
                END DO
            END DO
            !$omp end do
            !$omp barrier

            ! NEUE BERECHNUNG DES KONVEKTIVEN FLUSSES
            ! FUER JEDE GROBGITTERGESCHWINDIGKEIT
            k = kk-2
            !$omp do private(i, j, ddxi, dyj, avz, fvt, qkvt)
            DO i = 3, ii-2, 2
                ddxi = ddx(i)

                ! ym-rand
                j = 2
                dyj = dy(j)
                avz = ddxi*dyj
                fvt = avz*(2*wt(k, j, i) + 2*wt(k, j, i+1) &
                    + wt(k, j+1, i) + wt(k, j+1, i+1) &
                    + wt(k, j+2, i) + wt(k, j+2, i+1))*0.125
                qkvt = 0.5*fvt*(v(k, j, i) + v(k+1, j, i))
                wcv(k, j, i) = qkvt
                qkvt = 0.5*fvt*(v(k, j, i+1) + v(k+1, j, i+1))
                wcv(k, j, i+1) = qkvt

                ! IM GEBIET
                DO j = 4, jj-4, 2
                    dyj = dy(j)
                    avz = ddxi*dyj
                    fvt = avz*(wt(k, j-1, i) + wt(k, j-1, i+1) &
                        + wt(k, j, i) + wt(k, j, i+1) &
                        + wt(k, j+1, i) + wt(k, j+1, i+1) &
                        + wt(k, j+2, i) + wt(k, j+2, i+1))*0.125
                    qkvt = 0.5*fvt*(v(k, j, i) + v(k+1, j, i))
                    wcv(k, j, i) = qkvt
                    qkvt = 0.5*fvt*(v(k, j, i+1) + v(k+1, j, i+1))
                    wcv(k, j, i+1) = qkvt
                END DO

                ! yp-rand
                j = jj-2
                dyj = dy(j)
                avz = ddxi*dyj
                fvt = avz*(wt(k, j-1, i) + wt(k, j-1, i+1) &
                    + wt(k, j, i) + wt(k, j, i+1) &
                    + 2*wt(k, j+1, i) + 2*wt(k, j+1, i+1))*0.125
                qkvt = 0.5*fvt*(v(k, j, i) + v(k+1, j, i))
                wcv(k, j, i) = qkvt
                qkvt = 0.5*fvt*(v(k, j, i+1) + v(k+1, j, i+1))
                wcv(k, j, i+1) = qkvt
            END DO
            !$omp end do
            !$omp barrier

            ! VERTEILUNG
            k = kk-2
            !$omp do collapse(2) private(i, j)
            DO i = 3, ii-2
                DO j = 2, jj-4, 2
                    wcv(k, j+1, i) = 0.5*(wcv(k, j, i) + wcv(k, j+2, i))
                END DO
            END DO
            !$omp end do
            !$omp barrier

            ! AUF VO SCHREIBEN
            rddzk = rddz(k)
            !$omp do collapse(2) private(i, j, fkdtv)
            DO i = 3, ii-2
                DO j = 2, jj-2
                    fkdtv = -1.0*rddx(i)*rdy(j)*wkon
                    vo(k, j, i) = vo(k, j, i) &
                        + fkdtv*rddzk*(wcv(k, j, i) - wcv(k+1, j, i))
                END DO
            END DO
            !$omp end do
            !$omp barrier
        END IF

        ! PAR-RB Impulserhaltend BOTTOM
        IF (nbot == 8) THEN
            ! U-Impulszelle
            ! STANDARD-BERECHNUNG DES KONVEKTIVEN FLUSSES
            k = 3
            !$omp do private(i, j, dxi, ddyj, auz, fub, qkub)
            DO i = 2, ii-2
                dxi = dx(i)
                DO j = 3, jj-2
                    ddyj = ddy(j)
                    auz = dxi*ddyj
                    fub = auz*(wt(k-1, j, i) + wt(k-1, j, i+1))*0.5
                    qkub = 0.5*fub*(u(k, j, i) + u(k-1, j, i))
                    wcu(k-1, j, i) = qkub
                END DO
            END DO
            !$omp end do
            !$omp barrier

            ! NEUE BERECHNUNG DES KONVEKTIVEN FLUSSES
            ! FUER JEDE GROBGITTERGESCHWINDIGKEIT
            k = 3

            ! XM-RAND
            i = 2
            dxi = dx(i)
            !$omp do private(j, ddyj, auz, fub, qkub)
            DO j = 3, jj-2, 2
                ddyj = ddy(j)
                auz = dxi*ddyj
                fub = auz*(2*wt(k-1, j, i) + 2*wt(k-1, j+1, i) &
                    + wt(k-1, j, i+1) + wt(k-1, j+1, i+1) &
                    + wt(k-1, j, i+2) + wt(k-1, j+1, i+2))*0.125
                qkub = 0.5*fub*(u(k, j, i) + u(k-1, j, i))
                wcu(k, j, i) = qkub
                qkub = 0.5*fub*(u(k, j+1, i) + u(k-1, j+1, i))
                wcu(k, j+1, i) = qkub
            END DO
            !$omp end do
            !$omp barrier

            ! IM GEBIET
            !$omp do private(i, j, dxi, ddyj, auz, fub, qkub)
            DO i = 4, ii-4, 2
                dxi = dx(i)
                DO j = 3, jj-2, 2
                    ddyj = ddy(j)
                    auz = dxi*ddyj
                    fub = auz*(wt(k-1, j, i-1) + wt(k-1, j+1, i-1) &
                        + wt(k-1, j, i) + wt(k-1, j+1, i) &
                        + wt(k-1, j, i+1) + wt(k-1, j+1, i+1) &
                        + wt(k-1, j, i+2) + wt(k-1, j+1, i+2))*0.125
                    qkub = 0.5*fub*(u(k, j, i) + u(k-1, j, i))
                    wcu(k, j, i) = qkub
                    qkub = 0.5*fub*(u(k, j+1, i) + u(k-1, j+1, i))
                    wcu(k, j+1, i) = qkub
                END DO
            END DO
            !$omp end do
            !$omp barrier

            ! XP-RAND
            i = ii-2
            dxi = dx(i)
            !$omp do private(j, ddyj, auz, fub, qkub)
            DO j = 3, jj-2, 2
                ddyj = ddy(j)
                auz = dxi*ddyj
                fub = auz*(wt(k-1, j, i-1) + wt(k-1, j+1, i-1) &
                    + wt(k-1, j, i) + wt(k-1, j+1, i) &
                    + 2*wt(k-1, j, i+1) + 2*wt(k-1, j+1, i+1))*0.125
                qkub = 0.5*fub*(u(k, j, i) + u(k-1, j, i))
                wcu(k, j, i) = qkub
                qkub = 0.5*fub*(u(k, j+1, i) + u(k-1, j+1, i))
                wcu(k, j+1, i) = qkub
            END DO
            !$omp end do
            !$omp barrier

            ! VERTEILUNG
            k = 3
            !$omp do collapse(2) private(i, j)
            DO i = 2, ii-4, 2
                DO j = 3, jj-2
                    wcu(k, j, i+1) = 0.5*(wcu(k, j, i) + wcu(k, j, i+2))
                END DO
            END DO
            !$omp end do
            !$omp barrier

            ! AUF U0 SCHREIBEN
            rddzk = rddz(k)
            !$omp do collapse(2) private(i, j, fkdtu)
            DO i = 2, ii-2
                DO j = 3, jj-2
                    fkdtu = -1.0*rddy(j)*rdx(i)*wkon
                    uo(k, j, i) = uo(k, j, i) &
                        + fkdtu*rddzk*(-wcu(k, j, i) + wcu(k-1, j, i))
                END DO
            END DO
            !$omp end do
            !$omp barrier

            ! V-Impulszelle
            ! STANDARD-BERECHNUNG DES KONVEKTIVEN FLUSSES
            k = 3
            !$omp do private(i, j, ddxi, dyj, avz, fvb, qkvb)
            DO i = 3, ii-2
                ddxi = ddx(i)
                DO j = 2, jj-2
                    dyj = dy(j)
                    avz = ddxi*dyj
                    fvb = avz*(wt(k-1, j, i) + wt(k-1, j+1, i))*0.5
                    qkvb = 0.5*fvb*(v(k, j, i) + v(k-1, j, i))
                    wcv(k-1, j, i) = qkvb
                END DO
            END DO
            !$omp end do
            !$omp barrier

            ! NEUE BERECHNUNG DES KONVEKTIVEN FLUSSES
            ! FUER JEDE GROBGITTERGESCHWINDIGKEIT
            k = 3
            !$omp do private(i, j, ddxi, dyj, avz, fvb, qkvb)
            DO i = 3, ii-2, 2
                ddxi = ddx(i)

                ! ym-rand
                j = 2
                dyj = dy(j)
                avz = ddxi* dyj
                fvb = avz*(2*wt(k-1, j, i) + 2*wt(k-1, j, i+1) &
                    + wt(k-1, j+1, i) + wt(k-1, j+1, i+1) &
                    + wt(k-1, j+2, i) + wt(k-1, j+2, i+1))*0.125
                qkvb = 0.5*fvb*(v(k, j, i) + v(k-1, j, i))
                wcv(k, j, i) = qkvb
                qkvb = 0.5*fvb*(v(k, j, i+1) + v(k-1, j, i+1))
                wcv(k, j, i+1) = qkvb

                ! IM GEBIET
                DO j = 4, jj-4, 2
                    dyj = dy(j)
                    avz = ddxi* dyj
                    fvb = avz*(wt(k-1, j-1, i) + wt(k-1, j-1, i+1) &
                        + wt(k-1, j, i) + wt(k-1, j, i+1) &
                        + wt(k-1, j+1, i) + wt(k-1, j+1, i+1) &
                        + wt(k-1, j+2, i) + wt(k-1, j+2, i+1))*0.125
                    qkvb = 0.5*fvb*(v(k, j, i) + v(k-1, j, i))
                    wcv(k, j, i) = qkvb
                    qkvb = 0.5*fvb*(v(k, j, i+1) + v(k-1, j, i+1))
                    wcv(k, j, i+1) = qkvb

                END DO

                ! yp-rand
                j = jj-2
                dyj = dy(j)
                avz = ddxi* dyj
                fvb = avz*(wt(k-1, j-1, i) + wt(k-1, j-1, i+1) &
                    + wt(k-1, j, i) + wt(k-1, j, i+1) &
                    + 2*wt(k-1, j+1, i) + 2*wt(k-1, j+1, i+1))*0.125
                qkvb = 0.5*fvb*(v(k, j, i) + v(k-1, j, i))
                wcv(k, j, i) = qkvb
                qkvb = 0.5*fvb*(v(k, j, i+1) + v(k-1, j, i+1))
                wcv(k, j, i+1) = qkvb
            END DO
            !$omp end do
            !$omp barrier

            ! VERTEILUNG
            k = 3
            !$omp do collapse(2) private(i, j)
            DO i = 3, ii-2
                DO j = 2, jj-4, 2
                    wcv(k, j+1, i) = 0.5*(wcv(k, j, i) + wcv(k, j+2, i))
                END DO
            END DO
            !$omp end do
            !$omp barrier

            ! AUF VO SCHREIBEN
            rddzk = rddz(k)
            !$omp do collapse(2) private(i, j, fkdtv)
            DO i = 3, ii-2
                DO j = 2, jj-2
                    fkdtv = -1.0*rddx(i)*rdy(j)*wkon
                    vo(k, j, i) = vo(k, j, i) &
                        + fkdtv*rddzk*(-wcv(k, j, i) + wcv(k-1, j, i))
                END DO
            END DO
            !$omp end do
            !$omp barrier
        END IF

        ! PAR-RB Impulserhaltend LEFT
        IF (nlft == 8) THEN
            ! U-Impulszelle
            ! STANDARD-BERECHNUNG DES KONVEKTIVEN FLUSSES
            j = jj-2
            !$omp do private(i, k, dxi, ddzk, auy, fun, qkun)
            DO i = 2, ii-2
                dxi = dx(i)
                DO k = 3, kk-2
                    ddzk = ddz(k)
                    auy = dxi*ddzk
                    fun = auy*(vt(k, j, i) + vt(k, j, i+1))*0.5
                    qkun = 0.5*fun*(u(k, j, i) + u(k, j+1, i))
                    wcu(k, j+1, i) = qkun
                END DO
            END DO
            !$omp end do
            !$omp barrier

            ! NEUE BERECHNUNG DES KONVEKTIVEN FLUSSES
            ! FUER JEDE GROBGITTERGESCHWINDIGKEIT
            j = jj-2

            ! XM-RAND
            i = 2
            dxi = dx(i)
            !$omp do private(k, ddzk, auy, fun, qkun)
            DO k = 3, kk-2, 2
                ddzk = ddz(k)
                auy = dxi*ddzk
                fun = auy*(2*vt(k, j, i) + 2*vt(k+1, j, i) &
                    + vt(k, j, i+1) + vt(k+1, j, i+1) &
                    + vt(k, j, i+2) + vt(k+1, j, i+2))*0.125
                qkun = 0.5*fun*(u(k, j, i) + u(k, j+1, i))
                wcu(k, j, i) = qkun
                qkun = 0.5*fun*(u(k+1, j, i) + u(k+1, j+1, i))
                wcu(k+1, j, i) = qkun
            END DO
            !$omp end do
            !$omp barrier

            ! IM GEBIET
            !$omp do private(i, k, dxi, ddzk, auy, fun, qkun)
            DO i = 4, ii-4, 2
                dxi = dx(i)
                DO k = 3, kk-2, 2
                    ddzk = ddz(k)
                    auy = dxi*ddzk
                    fun = auy*(vt(k, j, i-1) + vt(k+1, j, i-1) &
                        + vt(k, j, i) + vt(k+1, j, i) &
                        + vt(k, j, i+1) + vt(k+1, j, i+1) &
                        + vt(k, j, i+2) + vt(k+1, j, i+2))*0.125
                    qkun = 0.5*fun*(u(k, j, i) + u(k, j+1, i))
                    wcu(k, j, i) = qkun
                    qkun = 0.5*fun*(u(k+1, j, i) + u(k+1, j+1, i))
                    wcu(k+1, j, i) = qkun
                END DO
            END DO
            !$omp end do
            !$omp barrier

            ! XP-RAND
            i = ii-2
            dxi = dx(i)
            !$omp do private(k, ddzk, auy, fun, qkun)
            DO k = 3, kk-2, 2
                ddzk = ddz(k)
                auy = dxi*ddzk
                fun = auy*(vt(k, j, i-1) + vt(k+1, j, i-1) &
                    + vt(k, j, i) + vt(k+1, j, i) &
                    + 2*vt(k, j, i+1) + 2*vt(k+1, j, i+1))*0.125
                qkun = 0.5*fun*(u(k, j, i) + u(k, j+1, i))
                wcu(k, j, i) = qkun
                qkun = 0.5*fun*(u(k+1, j, i) + u(k+1, j+1, i))
                wcu(k+1, j, i) = qkun
            END DO
            !$omp end do
            !$omp barrier

            ! VERTEILUNG
            j = jj-2
            !$omp do collapse(2) private(i, k)
            DO i = 2, ii-4, 2
                DO k = 3, kk-2
                    wcu(k, j, i+1) = 0.5*(wcu(k, j, i) + wcu(k, j, i+2))
                END DO
            END DO
            !$omp end do
            !$omp barrier

            ! AUF U0 SCHREIBEN
            j = jj-2
            !$omp do private(i, k, fkdtu, rddzk)
            DO i = 2, ii-2
                fkdtu = -1.0*rddy(j)*rdx(i)*wkon
                DO k = 3, kk-2
                    rddzk = rddz(k)
                    uo(k, j, i) = uo(k, j, i) &
                        + fkdtu*rddzk*(wcu(k, j, i) - wcu(k, j+1, i))
                END DO
            END DO
            !$omp end do
            !$omp barrier

            ! W-IMPULSZELLE
            ! STANDARD-BERECHNUNG DES KONVEKTIVEN FLUSSES
            j = jj-2
            !$omp do private(i, k, ddxi, dzk, awy, fwn, qkwn)
            DO i = 3, ii-2
                ddxi = ddx(i)
                DO k = 2, kk-2
                    dzk = dz(k)
                    awy = ddxi*dzk
                    fwn = awy*(vt(k, j, i) + vt(k+1, j, i))*0.5
                    qkwn = 0.5*fwn*(w(k, j, i) + w(k, j+1, i))
                    wcw(k, j+1, i) = qkwn
                END DO
            END DO
            !$omp end do
            !$omp barrier

            ! NEUE BERECHNUNG DES KONVEKTIVEN FLUSSES
            ! FUER JEDE GROBGITTERGESCHWINDIGKEIT
            j = jj-2
            !$omp do private(i, k, ddxi, dzk, awy, fwn, qkwn)
            DO i = 3, ii-2, 2
                ddxi = ddx(i)

                ! ZM-RAND
                k = 2
                dzk = dz(k)
                awy = ddxi*dzk
                fwn = awy*(2*vt(k, j, i) + 2*vt(k, j, i+1) &
                    + vt(k+1, j, i) + vt(k+1, j, i+1) &
                    + vt(k+2, j, i) + vt(k+2, j, i+1))*0.125
                qkwn = 0.5*fwn*(w(k, j, i) + w(k, j+1, i))
                wcw(k, j, i) = qkwn
                qkwn = 0.5*fwn*(w(k, j, i+1) + w(k, j+1, i+1))
                wcw(k, j, i+1) = qkwn

                ! IM-GEBIET
                DO k = 4, kk-4
                    dzk = dz(k)
                    awy = ddxi*dzk
                    fwn = awy*(vt(k-1, j, i) + vt(k-1, j, i+1) &
                        + vt(k, j, i) + vt(k, j, i+1) &
                        + vt(k+1, j, i) + vt(k+1, j, i+1) &
                        + vt(k+2, j, i) + vt(k+2, j, i+1))*0.125
                    qkwn = 0.5*fwn*(w(k, j, i) + w(k, j+1, i))
                    wcw(k, j, i) = qkwn
                    qkwn = 0.5*fwn*(w(k, j, i+1) + w(k, j+1, i+1))
                    wcw(k, j, i+1) = qkwn
                END DO

                ! ZP-RAND
                k = kk-2
                dzk = dz(k)
                awy = ddxi*dzk
                fwn = awy*(vt(k-1, j, i) + vt(k-1, j, i+1) &
                    + vt(k, j, i) + vt(k, j, i+1) &
                    + 2*vt(k+1, j, i) + 2*vt(k+1, j, i+1))*0.125
                qkwn = 0.5*fwn*(w(k, j, i) + w(k, j+1, i))
                wcw(k, j, i) = qkwn
                qkwn = 0.5*fwn*(w(k, j, i+1) + w(k, j+1, i+1))
                wcw(k, j, i+1) = qkwn
            END DO
            !$omp end do
            !$omp barrier

            ! VERTEILUNG
            j = jj-2
            !$omp do collapse(2) private(i, k)
            DO i = 3, ii-2
                DO k = 2, kk-4, 2
                    wcw(k+1, j, i) = 0.5*(wcw(k, j, i) + wcw(k+2, j, i))
                END DO
            END DO
            !$omp end do
            !$omp barrier

            ! AUF W0 SCHREIBEN
            j = jj-2
            !$omp do private(i, k, fkdtw, rdzk)
            DO i = 3, ii-2
                fkdtw = -1.0*rddx(i)*rddy(j)*wkon
                DO k = 2, kk-2
                    rdzk = rdz(k)
                    wo(k, j, i) = wo(k, j, i) &
                        + fkdtw*rdzk*(-wcw(k, j+1, i) + wcw(k, j, i))
                END DO
            END DO
            !$omp end do
            !$omp barrier
        END IF

        ! PAR-RB Impulserhaltend RIGHT
        IF (nrgt == 8) THEN
            ! U-Impulszelle
            ! STANDARD-BERECHNUNG DES KONVEKTIVEN FLUSSES
            j = 3
            !$omp do private(i, k, dxi, ddzk, auy, fus, qkus)
            DO i = 2, ii-2
                dxi = dx(i)
                DO k = 3, kk-2
                    ddzk = ddz(k)
                    auy = dxi*ddzk
                    fus = auy*(vt(k, j-1, i) + vt(k, j-1, i+1))*0.5
                    qkus = 0.5*fus*(u(k, j, i) + u(k, j-1, i))
                    wcu(k, j-1, i) = qkus
                END DO
            END DO
            !$omp end do
            !$omp barrier

            ! NEUE BERECHNUNG DES KONVEKTIVEN FLUSSES
            ! FUER JEDE GROBGITTERGESCHWINDIGKEIT
            j = 3

            ! XM-RAND
            i = 2
            dxi = dx(i)
            !$omp do private(k, ddzk, auy, fus, qkus)
            DO k = 3, kk-2, 2
                ddzk = ddz(k)
                auy = dxi*ddzk
                fus = auy*(2*vt(k, j-1, i) + 2*vt(k+1, j-1, i) &
                    + vt(k, j-1, i+1) + vt(k+1, j-1, i+1) &
                    + vt(k, j-1, i+2) + vt(k+1, j-1, i+2))*0.125
                qkus = 0.5*fus*(u(k, j, i) + u(k, j-1, i))
                wcu(k, j, i) = qkus
                qkus = 0.5*fus*(u(k+1, j, i) + u(k+1, j-1, i))
                wcu(k+1, j, i) = qkus
            END DO
            !$omp end do
            !$omp barrier

            ! IM GEBIET
            !$omp do private(i, k, dxi, ddzk, auy, fus, qkus)
            DO i = 4, ii-4, 2
                dxi = dx(i)
                DO k = 3, kk-2, 2
                    ddzk = ddz(k)
                    auy = dxi*ddzk
                    fus = auy*(vt(k, j-1, i-1) + vt(k+1, j-1, i-1) &
                        + vt(k, j-1, i) + vt(k+1, j-1, i) &
                        + vt(k, j-1, i+1) + vt(k+1, j-1, i+1) &
                        + vt(k, j-1, i+2) + vt(k+1, j-1, i+2))*0.125
                    qkus = 0.5*fus*(u(k, j, i) + u(k, j-1, i))
                    wcu(k, j, i) = qkus
                    qkus = 0.5*fus*(u(k+1, j, i) + u(k+1, j-1, i))
                    wcu(k+1, j, i) = qkus
                END DO
            END DO
            !$omp end do
            !$omp barrier

            ! XP-RAND
            i = ii-2
            dxi = dx(i)
            !$omp do private(k, ddzk, auy, fus, qkus)
            DO k = 3, kk-2, 2
                ddzk = ddz(k)
                auy = dxi*ddzk
                fus = auy*(vt(k, j-1, i-1) + vt(k+1, j-1, i-1) &
                    + vt(k, j-1, i) + vt(k+1, j-1, i) &
                    +2 *vt(k, j-1, i+1) + 2*vt(k+1, j-1, i+1))*0.125
                qkus = 0.5*fus*(u(k, j, i) + u(k, j-1, i))
                wcu(k, j, i) = qkus
                qkus = 0.5*fus*(u(k+1, j, i) + u(k+1, j-1, i))
                wcu(k+1, j, i) = qkus
            END DO
            !$omp end do
            !$omp barrier

            ! VERTEILUNG
            j = 3
            !$omp do collapse(2) private(i, k)
            DO i = 2, ii-4, 2
                DO k = 3, kk-2
                    wcu(k, j, i+1) = 0.5*(wcu(k, j, i) + wcu(k, j, i+2))
                END DO
            END DO
            !$omp end do
            !$omp barrier

            ! AUF U0 SCHREIBEN
            j = 3
            !$omp do private(i, k, fkdtu, rddzk)
            DO i = 2, ii-2
                fkdtu = -1.0*rddy(j)*rdx(i)*wkon
                DO k = 3, kk-2
                    rddzk = rddz(k)
                    uo(k, j, i) = uo(k, j, i) &
                        + fkdtu*rddzk*(-wcu(k, j, i) + wcu(k, j-1, i))
                END DO
            END DO
            !$omp end do
            !$omp barrier

            ! W-IMPULSZELLE
            ! STANDARD-BERECHNUNG DES KONVEKTIVEN FLUSSES
            j = 3
            !$omp do private(i, k, ddxi, dzk, awy, fws, qkws)
            DO i = 3, ii-2
                ddxi = ddx(i)
                DO k = 2, kk-2
                    dzk = dz(k)
                    awy = ddxi*dzk
                    fws = awy*(vt(k, j-1, i) + vt(k+1, j-1, i))*0.5
                    qkws = 0.5*fws*(w(k, j, i) + w(k, j-1, i))
                    wcw(k, j-1, i) = qkws
                END DO
            END DO
            !$omp end do
            !$omp barrier

            ! NEUE BERECHNUNG DES KONVEKTIVEN FLUSSES
            ! FUER JEDE GROBGITTERGESCHWINDIGKEIT
            j = 3
            !$omp do private(i, k, ddxi, dzk, awy, fws, qkws)
            DO i = 3, ii-2, 2
                ddxi = ddx(i)

                ! zm-rand
                k = 2
                dzk = dz(k)
                awy = ddxi*dzk
                fws = awy*(2*vt(k, j-1, i) + 2*vt(k, j-1, i+1) &
                    + vt(k+1, j-1, i) + vt(k+1, j-1, i+1) &
                    + vt(k+2, j-1, i) + vt(k+2, j-1, i+1))*0.125
                qkws = 0.5*fws*(w(k, j, i) + w(k, j-1, i))
                wcw(k, j, i) = qkws
                qkws = 0.5*fws*(w(k, j, i+1) + w(k, j-1, i+1))
                wcw(k, j, i+1) = qkws

                ! IM-GEBIET
                DO k = 4, kk-4
                    dzk = dz(k)
                    awy = ddxi*dzk
                    fws = awy*(vt(k-1, j-1, i) + vt(k-1, j-1, i+1) &
                        + vt(k, j-1, i) + vt(k, j-1, i+1) &
                        + vt(k+1, j-1, i) + vt(k+1, j-1, i+1) &
                        + vt(k+2, j-1, i) + vt(k+2, j-1, i+1))*0.125
                    qkws = 0.5*fws*(w(k, j, i) + w(k, j-1, i))
                    wcw(k, j, i) = qkws
                    qkws = 0.5*fws*(w(k, j, i+1) + w(k, j-1, i+1))
                    wcw(k, j, i+1) = qkws
                END DO

                ! ZP-RAND
                k = kk-2
                dzk = dz(k)
                awy = ddxi*dzk
                fws = awy*(vt(k-1, j-1, i) + vt(k-1, j-1, i+1) &
                    + vt(k, j-1, i) + vt(k, j-1, i+1) &
                    + 2*vt(k+1, j-1, i) + 2*vt(k+1, j-1, i+1))*0.125
                qkws = 0.5*fws*(w(k, j, i) + w(k, j-1, i))
                wcw(k, j, i) = qkws
                qkws = 0.5*fws*(w(k, j, i+1) + w(k, j-1, i+1))
                wcw(k, j, i+1) = qkws
            END DO
            !$omp end do
            !$omp barrier

            ! VERTEILUNG
            j = 3
            !$omp do collapse(2) private(i, k)
            DO i = 3, ii-2
                DO k = 2, kk-4, 2
                    wcw(k+1, j, i) = 0.5*(wcw(k, j, i) + wcw(k+2, j, i))
                END DO
            END DO
            !$omp end do
            !$omp barrier

            ! AUF W0 SCHREIBEN
            j = 3
            !$omp do private(i, k, fkdtw, rdzk)
            DO i = 3, ii-2
                fkdtw = -1.0*rddx(i)*rddy(j)*wkon
                DO k = 2, kk-2
                    rdzk = rdz(k)
                    wo(k, j, i) = wo(k, j, i) &
                        + fkdtw*rdzk*(wcw(k, j-1, i) - wcw(k, j, i))
                END DO
            END DO
            !$omp end do
            ! no need for a barrier here since this is the last operation in
            ! the subroutine
        END IF
    END SUBROUTINE tstle4_par


    SUBROUTINE swcle3d(kk, jj, ii, uo, vo, wo, u, v, w, ddx, ddy, ddz, &
            nfro, nbac, nrgt, nlft, nbot, ntop)
        !$omp declare target

        ! Subroutine arguments
        INTEGER(intk), INTENT(in) :: kk, jj, ii
        REAL(realk), INTENT(inout) :: uo(kk, jj, ii), vo(kk, jj, ii), &
            wo(kk, jj, ii)
        REAL(realk), INTENT(in) :: u(kk, jj, ii), v(kk, jj, ii), w(kk, jj, ii)
        REAL(realk), INTENT(in) :: ddx(ii), ddy(jj), ddz(kk)
        INTEGER(intk), INTENT(in) :: nfro, nbac, nrgt, nlft, nbot, ntop

        ! Local variables
        INTEGER(intk) :: k, j, i

        IF (nfro == 5) THEN
            i = 3
            !$omp do collapse(2) private(j, k)
            DO j = 2, jj-1
                DO k = 2, kk-1
                    vo(k, j, i) = vo(k, j, i) - swcle3d_one(ddx(i), v(k, j, i))
                    wo(k, j, i) = wo(k, j, i) - swcle3d_one(ddx(i), w(k, j, i))
                END DO
            END DO
            !$omp end do
            !$omp barrier
        END IF

        IF (nbac == 5) THEN
            i = ii-2
            !$omp do collapse(2) private(j, k)
            DO j = 2, jj-1
                DO k = 2, kk-1
                    vo(k, j, i) = vo(k, j, i) - swcle3d_one(ddx(i), v(k, j, i))
                    wo(k, j, i) = wo(k, j, i) - swcle3d_one(ddx(i), w(k, j, i))
                END DO
            END DO
            !$omp end do
            !$omp barrier
        END IF

        IF (nrgt == 5) THEN
            j = 3
            !$omp do collapse(2) private(i, k)
            DO i = 2, ii-1
                DO k = 2, kk-1
                    uo(k, j, i) = uo(k, j, i) - swcle3d_one(ddy(j), u(k, j, i))
                    wo(k, j, i) = wo(k, j, i) - swcle3d_one(ddy(j), w(k, j, i))
                END DO
            END DO
            !$omp end do
            !$omp barrier
        END IF

        IF (nlft == 5) THEN
            j = jj-2
            !$omp do collapse(2) private(i, k)
            DO i = 2, ii-1
                DO k = 2, kk-1
                    uo(k, j, i) = uo(k, j, i) - swcle3d_one(ddy(j), u(k, j, i))
                    wo(k, j, i) = wo(k, j, i) - swcle3d_one(ddy(j), w(k, j, i))
                END DO
            END DO
            !$omp end do
            !$omp barrier
        END IF

        IF (nbot == 5) THEN
            k = 3
            !$omp do collapse(2) private(i, j)
            DO i = 2, ii-1
                DO j = 2, jj-1
                    uo(k, j, i) = uo(k, j, i) - swcle3d_one(ddz(k), u(k, j, i))
                    vo(k, j, i) = vo(k, j, i) - swcle3d_one(ddz(k), v(k, j, i))
                END DO
            END DO
            !$omp end do
            !$omp barrier
        END IF

        IF (ntop == 5) THEN
            k = kk-2
            !$omp do collapse(2) private(i, j)
            DO i = 2, ii-1
                DO j = 2, jj-1
                    uo(k, j, i) = uo(k, j, i) - swcle3d_one(ddz(k), u(k, j, i))
                    vo(k, j, i) = vo(k, j, i) - swcle3d_one(ddz(k), v(k, j, i))
                END DO
            END DO
            !$omp end do
            ! no need for a barrier here since this is the last operation in
            ! the subroutine
        END IF
    END SUBROUTINE swcle3d


    PURE ELEMENTAL REAL(realk) FUNCTION swcle3d_one(ddz, u) RESULT(uo)
        !$omp declare target
#ifndef _MGLET_OFFLOAD_
        !$omp declare simd(swcle3d_one)
#endif

        ! Function arguments
        REAL(realk), INTENT(in) :: ddz  ! wall normal
        REAL(realk), INTENT(in) :: u    ! velocity

        ! Local variables
        ! none...

        uo = tauwin(u, ddz)/rho/ddz
    END FUNCTION swcle3d_one
END MODULE tstle4_mod
