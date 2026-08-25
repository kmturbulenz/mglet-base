MODULE tstle4_mod
    USE core_mod
    USE flowcore_mod
    USE lesmodel_mod, ONLY: ilesmodel

    IMPLICIT NONE(type, external)
    PRIVATE

    ! Bundling grid information for C processing

    TYPE, BIND(C) :: tstle4_grid_t
        ! Grid dimensions
        INTEGER(c_intk) :: ii
        INTEGER(c_intk) :: jj
        INTEGER(c_intk) :: kk
        ! Array indices for the 3D / 1D arrays
        INTEGER(c_intk) :: ip3
        INTEGER(c_intk) :: ipx
        INTEGER(c_intk) :: ipy
        INTEGER(c_intk) :: ipz
        ! Boundary information for the grid
        INTEGER(c_intk) :: nfro
        INTEGER(c_intk) :: nbac
        INTEGER(c_intk) :: nrgt
        INTEGER(c_intk) :: nlft
        INTEGER(c_intk) :: nbot
        INTEGER(c_intk) :: ntop
        ! Information on pressure gradient
        REAL(c_realk) :: gpx
        REAL(c_realk) :: gpy
        REAL(c_realk) :: gpz
    END TYPE tstle4_grid_t

    TYPE(tstle4_grid_t), ALLOCATABLE :: tstle4_grids(:)
    !$omp declare target(tstle4_grids)

    LOGICAL, PROTECTED :: is_initialized = .FALSE.


    PUBLIC :: tstle4, initialize_tstle4, finalize_tstle4

    INTERFACE

        SUBROUTINE tstle4_kon_u_c(kk, jj, ii, uo, u, v, w, ut, vt, wt, dx, &
            dy, dz, ddx, ddy, ddz, nfro, &
            nbac) BIND(C, name="tstle4_kon_u_c")
            IMPORT :: c_intk, c_realk
            INTEGER(c_intk), VALUE, INTENT(in) :: kk, jj, ii
            REAL(c_realk), INTENT(inout) :: uo(*)
            REAL(c_realk), INTENT(in) :: u(*), v(*), w(*), ut(*), vt(*), wt(*)
            REAL(c_realk), INTENT(in) :: dx(*), dy(*), dz(*)
            REAL(c_realk), INTENT(in) :: ddx(*), ddy(*), ddz(*)
            INTEGER(c_intk), VALUE, INTENT(in) :: nfro, nbac
        END SUBROUTINE tstle4_kon_u_c

        SUBROUTINE tstle4_kon_v_c(kk, jj, ii, vo, u, v, w, ut, vt, wt, dx, &
            dy, dz, ddx, ddy, ddz, nrgt, &
            nlft) BIND(C, name="tstle4_kon_v_c")
            IMPORT :: c_intk, c_realk
            INTEGER(c_intk), VALUE, INTENT(in) :: kk, jj, ii
            REAL(c_realk), INTENT(inout) :: vo(*)
            REAL(c_realk), INTENT(in) :: u(*), v(*), w(*), ut(*), vt(*), wt(*)
            REAL(c_realk), INTENT(in) :: dx(*), dy(*), dz(*)
            REAL(c_realk), INTENT(in) :: ddx(*), ddy(*), ddz(*)
            INTEGER(c_intk), VALUE, INTENT(in) :: nrgt, nlft
        END SUBROUTINE tstle4_kon_v_c

        SUBROUTINE tstle4_kon_w_c(kk, jj, ii, wo, u, v, w, ut, vt, wt, dx, &
            dy, dz, ddx, ddy, ddz, nbot, &
            ntop) BIND(C, name="tstle4_kon_w_c")
            IMPORT :: c_intk, c_realk
            INTEGER(c_intk), VALUE, INTENT(in) :: kk, jj, ii
            REAL(c_realk), INTENT(inout) :: wo(*)
            REAL(c_realk), INTENT(in) :: u(*), v(*), w(*), ut(*), vt(*), wt(*)
            REAL(c_realk), INTENT(in) :: dx(*), dy(*), dz(*)
            REAL(c_realk), INTENT(in) :: ddx(*), ddy(*), ddz(*)
            INTEGER(c_intk), VALUE, INTENT(in) :: nbot, ntop
        END SUBROUTINE tstle4_kon_w_c

        SUBROUTINE tstle4_par_c(kk, jj, ii, uo, vo, wo, u, v, w, ut, vt, wt, &
            dx, dy, dz, ddx, ddy, ddz, rdx, rdy, rdz, rddx, rddy, rddz, &
            wcu, wcv, wcw, nfro, nbac, nrgt, nlft, nbot, ntop) &
            BIND(C, name="tstle4_par_c")
            IMPORT :: c_intk, c_realk
            INTEGER(c_intk), VALUE, INTENT(in) :: kk, jj, ii
            REAL(c_realk), INTENT(inout) :: uo(*), vo(*), wo(*)
            REAL(c_realk), INTENT(in) :: u(*), v(*), w(*), ut(*), vt(*), wt(*)
            REAL(c_realk), INTENT(in) :: dx(*), dy(*), dz(*)
            REAL(c_realk), INTENT(in) :: ddx(*), ddy(*), ddz(*)
            REAL(c_realk), INTENT(in) :: rdx(*), rdy(*), rdz(*)
            REAL(c_realk), INTENT(in) :: rddx(*), rddy(*), rddz(*)
            REAL(c_realk), INTENT(inout) :: wcu(*), wcv(*), wcw(*)
            INTEGER(c_intk), VALUE, INTENT(in) :: nfro, nbac, nrgt, nlft
            INTEGER(c_intk), VALUE, INTENT(in) :: nbot, ntop
        END SUBROUTINE tstle4_par_c

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

        SUBROUTINE tstle4_kon_impl_c(uo, vo, wo, u, v, w, ut, vt, wt, dx, dy, &
            dz, ddx, ddy, ddz, nmygrids_in, &
            grids) BIND(C, name="tstle4_kon_impl_c")
            IMPORT :: c_intk, c_realk, tstle4_grid_t
            REAL(c_realk), INTENT(inout) :: uo(*), vo(*), wo(*)
            REAL(c_realk), INTENT(in) :: u(*), v(*), w(*), ut(*), vt(*), wt(*)
            REAL(c_realk), INTENT(in) :: dx(*), dy(*), dz(*)
            REAL(c_realk), INTENT(in) :: ddx(*), ddy(*), ddz(*)
            INTEGER(c_intk), VALUE, INTENT(in) :: nmygrids_in
            TYPE(tstle4_grid_t), INTENT(in) :: grids(*)
        END SUBROUTINE tstle4_kon_impl_c

        SUBROUTINE tstle4_diff_impl_c(uo, vo, wo, u, v, w, g, dx, dy, dz, &
            ddx, ddy, ddz, ilesmodel_in, gmol_in, rho_in, nmygrids_in, grids) &
            BIND(C, name="tstle4_diff_impl_c")
            IMPORT :: c_intk, c_realk, tstle4_grid_t
            REAL(c_realk), INTENT(inout) :: uo(*), vo(*), wo(*)
            REAL(c_realk), INTENT(in) :: u(*), v(*), w(*), g(*)
            REAL(c_realk), INTENT(in) :: dx(*), dy(*), dz(*)
            REAL(c_realk), INTENT(in) :: ddx(*), ddy(*), ddz(*)
            INTEGER(c_intk), VALUE, INTENT(in) :: ilesmodel_in
            REAL(c_realk), VALUE, INTENT(in) :: gmol_in, rho_in
            INTEGER(c_intk), VALUE, INTENT(in) :: nmygrids_in
            TYPE(tstle4_grid_t), INTENT(in) :: grids(*)
        END SUBROUTINE tstle4_diff_impl_c


        SUBROUTINE tstle4_gradp_impl_c(uo, vo, wo, p, dx, dy, dz, rho_in, &
            nmygrids_in, grids) &
            BIND(C, name="tstle4_gradp_impl_c")
            IMPORT :: c_intk, c_realk, tstle4_grid_t
            REAL(c_realk), INTENT(inout) :: uo(*), vo(*), wo(*)
            REAL(c_realk), INTENT(in) :: p(*), dx(*), dy(*), dz(*)
            REAL(c_realk), VALUE, INTENT(in) :: rho_in
            INTEGER(c_intk), VALUE, INTENT(in) :: nmygrids_in
            TYPE(tstle4_grid_t), INTENT(in) :: grids(*)
        END SUBROUTINE tstle4_gradp_impl_c
    END INTERFACE

CONTAINS

    SUBROUTINE initialize_tstle4()

        ! Local variables
        INTEGER(intk) :: i, gradpflag
        REAL(c_realk) :: rflag

        ! Initialize the TSTLE4 module and create array for C processing
        ALLOCATE(tstle4_grids(nmygrids))

        DO i = 1, nmygrids
            ! Feeding dimensions
            CALL get_mgdims(tstle4_grids(i)%kk, tstle4_grids(i)%jj, &
                tstle4_grids(i)%ii, mygrids(i))
            ! Feeding boundary information
            CALL get_mgbasb(tstle4_grids(i)%nfro, tstle4_grids(i)%nbac, &
                tstle4_grids(i)%nrgt, tstle4_grids(i)%nlft, &
                tstle4_grids(i)%nbot, tstle4_grids(i)%ntop, mygrids(i))
            ! Feeding array indices
            CALL get_ip3(tstle4_grids(i)%ip3, mygrids(i))
            CALL get_ip1x(tstle4_grids(i)%ipx, mygrids(i))
            CALL get_ip1y(tstle4_grids(i)%ipy, mygrids(i))
            CALL get_ip1z(tstle4_grids(i)%ipz, mygrids(i))
            ! Feeding pressure gradient
            CALL get_gradpxflag(gradpflag, mygrids(i))
            rflag = REAL(gradpflag, c_realk)
            tstle4_grids(i)%gpx = REAL(gradp(1)*rflag, c_realk)
            tstle4_grids(i)%gpy = REAL(gradp(2)*rflag, c_realk)
            tstle4_grids(i)%gpz = REAL(gradp(3)*rflag, c_realk)
        END DO
        !$omp target enter data map(to: tstle4_grids)
        is_initialized = .TRUE.

    END SUBROUTINE initialize_tstle4

    SUBROUTINE finalize_tstle4

        !$omp target exit data map(delete: tstle4_grids)
        DEALLOCATE(tstle4_grids)

        is_initialized = .FALSE.

    END SUBROUTINE finalize_tstle4




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
        TYPE(field_t), POINTER :: wcu_f, wcv_f, wcw_f
        CALL start_timer(310)

        IF (.NOT. is_initialized) THEN
            CALL initialize_tstle4()
        END IF

        CALL set_field_arr(uo_f, 0.0_realk, device=.TRUE.)
        CALL set_field_arr(vo_f, 0.0_realk, device=.TRUE.)
        CALL set_field_arr(wo_f, 0.0_realk, device=.TRUE.)

        CALL get_field(dx_f, "DX")
        CALL get_field(dy_f, "DY")
        CALL get_field(dz_f, "DZ")

        CALL get_field(ddx_f, "DDX")
        CALL get_field(ddy_f, "DDY")
        CALL get_field(ddz_f, "DDZ")

        ! SW: The kernels do not use the reciprocal grid spacing rdx, rdy, &
        ! rdz, rddx, rddy, rddz. As we are mainly memory bound, we can afford
        ! divisions within the kernel at the benefit of reduced memory traffic.

        CALL push_field(wcu_f, "TSTLE4_WCU")
        CALL push_field(wcv_f, "TSTLE4_WCV")
        CALL push_field(wcw_f, "TSTLE4_WCW")

        ! Convective terms
        CALL start_timer(311)
        CALL tstle4_kon_impl_c(uo_f%arr, vo_f%arr, wo_f%arr, u_f%arr, &
            v_f%arr, w_f%arr, ut_f%arr, vt_f%arr, wt_f%arr, dx_f%arr, &
            dy_f%arr, dz_f%arr, ddx_f%arr, ddy_f%arr, ddz_f%arr, &
            INT(nmygrids, c_intk), tstle4_grids)
        CALL stop_timer(311)

        ! Diffusion terms
        CALL start_timer(312)
        CALL tstle4_diff_impl_c(uo_f%arr, vo_f%arr, wo_f%arr, u_f%arr, &
            v_f%arr, w_f%arr, g_f%arr, dx_f%arr, dy_f%arr, dz_f%arr, &
            ddx_f%arr, ddy_f%arr, ddz_f%arr, INT(ilesmodel, c_intk), &
            REAL(gmol, c_realk), REAL(rho, c_realk), &
            INT(nmygrids, c_intk), tstle4_grids)
        CALL stop_timer(312)

        CALL start_timer(313)
        CALL tstle4_gradp_impl_c(uo_f%arr, vo_f%arr, wo_f%arr, p_f%arr, &
            dx_f%arr, dy_f%arr, dz_f%arr, REAL(rho, c_realk), &
            INT(nmygrids, c_intk), tstle4_grids)
        CALL stop_timer(313)

        ! CALL start_timer(314)
        ! CALL tstle4_par_impl(uo_f%arr, vo_f%arr, wo_f%arr, u_f%arr, &
        !     v_f%arr, w_f%arr, ut_f%arr, vt_f%arr, wt_f%arr, dx_f%arr, &
        !     dy_f%arr, dz_f%arr, ddx_f%arr, ddy_f%arr, ddz_f%arr, &
        !     rdx_f%arr, rdy_f%arr, rdz_f%arr, rddx_f%arr, rddy_f%arr, &
        !     rddz_f%arr, wcu_f%arr, wcv_f%arr, wcw_f%arr)
        ! CALL stop_timer(314)

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
                ddy(ipy), ddz(ipz), INT(nfro, c_intk), INT(nbac, c_intk))

            CALL tstle4_kon_v_c(INT(kk, c_intk), INT(jj, c_intk), &
                INT(ii, c_intk), vo(ip3), u(ip3), v(ip3), w(ip3), ut(ip3), &
                vt(ip3), wt(ip3), dx(ipx), dy(ipy), dz(ipz), ddx(ipx), &
                ddy(ipy), ddz(ipz), INT(nrgt, c_intk), INT(nlft, c_intk))

            CALL tstle4_kon_w_c(INT(kk, c_intk), INT(jj, c_intk), &
                INT(ii, c_intk), wo(ip3), u(ip3), v(ip3), w(ip3), ut(ip3), &
                vt(ip3), wt(ip3), dx(ipx), dy(ipy), dz(ipz), ddx(ipx), &
                ddy(ipy), ddz(ipz), INT(nbot, c_intk), INT(ntop, c_intk))
            !$omp end parallel
        END DO
        !$omp end target teams distribute

    END SUBROUTINE tstle4_kon_impl


    ! SUBROUTINE tstle4_gradp_impl(uo, vo, wo, p, dx, dy, dz)
    !     REAL(realk), INTENT(inout) :: uo(*), vo(*), wo(*)
    !     REAL(realk), INTENT(in) :: p(*), dx(*), dy(*), dz(*)

    !     INTEGER(intk) :: i
    !     INTEGER(intk) :: gradpflag

    !     DO i = 1, nmygrids
    !         CALL get_gradpxflag(gradpflag, mygrids(i))
    !         tstle4_grids(i)%gpx = REAL(gradp(1)*gradpflag, c_realk)
    !         tstle4_grids(i)%gpy = REAL(gradp(2)*gradpflag, c_realk)
    !         tstle4_grids(i)%gpz = REAL(gradp(3)*gradpflag, c_realk)
    !     END DO

    !     CALL tstle4_gradp_impl_c(uo, vo, wo, p, dx, dy, dz, &
    !         REAL(rho, c_realk), INT(nmygrids, c_intk), tstle4_grids)

    ! END SUBROUTINE tstle4_gradp_impl


END MODULE tstle4_mod
