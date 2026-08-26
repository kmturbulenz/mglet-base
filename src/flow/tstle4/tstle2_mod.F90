MODULE tstle2_mod

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


    PUBLIC :: tstle2, tstle2_init, tstle2_finish

    INTERFACE

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

        SUBROUTINE tstle4_par_impl_c(uo, vo, wo, u, v, w, ut, vt, wt, &
            dx, dy, dz, ddx, ddy, ddz, wcu, wcv, wcw, nmygrids_in, grids) &
            BIND(C, name="tstle4_par_impl_c")
            IMPORT :: c_intk, c_realk, tstle4_grid_t
            REAL(c_realk), INTENT(inout) :: uo(*), vo(*), wo(*)
            REAL(c_realk), INTENT(in) :: u(*), v(*), w(*), ut(*), vt(*), &
                wt(*)
            REAL(c_realk), INTENT(in) :: dx(*), dy(*), dz(*)
            REAL(c_realk), INTENT(in) :: ddx(*), ddy(*), ddz(*)
            REAL(c_realk), INTENT(inout) :: wcu(*), wcv(*), wcw(*)
            INTEGER(c_intk), VALUE, INTENT(in) :: nmygrids_in
            TYPE(tstle4_grid_t), INTENT(in) :: grids(*)
        END SUBROUTINE tstle4_par_impl_c

    END INTERFACE

CONTAINS

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


    SUBROUTINE tstle2(uo_f, vo_f, wo_f, u_f, v_f, w_f, ut_f, vt_f, wt_f, &
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
            CALL tstle2_init()
        END IF

        CALL zero_field_arr(uo_f, device=.TRUE.)
        CALL zero_field_arr(vo_f, device=.TRUE.)
        CALL zero_field_arr(wo_f, device=.TRUE.)

        CALL get_field(dx_f, "DX")
        CALL get_field(dy_f, "DY")
        CALL get_field(dz_f, "DZ")

        CALL get_field(ddx_f, "DDX")
        CALL get_field(ddy_f, "DDY")
        CALL get_field(ddz_f, "DDZ")

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

        ! Pressure gradient terms
        CALL start_timer(313)
        CALL tstle4_gradp_impl_c(uo_f%arr, vo_f%arr, wo_f%arr, p_f%arr, &
            dx_f%arr, dy_f%arr, dz_f%arr, REAL(rho, c_realk), &
            INT(nmygrids, c_intk), tstle4_grids)
        CALL stop_timer(313)

        ! Special treatment near the parent grid boundaries
        CALL start_timer(314)
        CALL tstle4_par_impl_c(uo_f%arr, vo_f%arr, wo_f%arr, u_f%arr, &
            v_f%arr, w_f%arr, ut_f%arr, vt_f%arr, wt_f%arr, dx_f%arr, &
            dy_f%arr, dz_f%arr, ddx_f%arr, ddy_f%arr, ddz_f%arr, &
            wcu_f%arr, wcv_f%arr, wcw_f%arr, &
            INT(nmygrids, c_intk), tstle4_grids)
        CALL stop_timer(314)

        CALL pop_field(wcw_f)
        CALL pop_field(wcv_f)
        CALL pop_field(wcu_f)

        CALL stop_timer(310)
    END SUBROUTINE tstle2



    SUBROUTINE tstle2_init()

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

    END SUBROUTINE tstle2_init



    SUBROUTINE tstle2_finish()

        !$omp target exit data map(delete: tstle4_grids)
        DEALLOCATE(tstle4_grids)
        is_initialized = .FALSE.

    END SUBROUTINE tstle2_finish


END MODULE tstle2_mod
