MODULE lesmodel_mod
    USE core_mod
    USE flowcore_mod
    USE ib_mod
    USE wernerwengle_mod

    IMPLICIT NONE(type, external)
    PRIVATE

    INTEGER(intk), PARAMETER :: nchar = 16
    CHARACTER(len=nchar) :: clesmodel
    INTEGER(intk), PROTECTED :: ilesmodel
    !$omp declare target(ilesmodel)

    TYPE, EXTENDS(bound_t) :: boundg_t
    CONTAINS
        PROCEDURE, NOPASS :: front => bfront
        PROCEDURE, NOPASS :: back => bfront
        PROCEDURE, NOPASS :: right => bright
        PROCEDURE, NOPASS :: left => bright
        PROCEDURE, NOPASS :: bottom => bbottom
        PROCEDURE, NOPASS :: top => bbottom
    END TYPE boundg_t

    TYPE, BIND(C) :: lesmodel_grid_t
        INTEGER(c_intk) :: ii
        INTEGER(c_intk) :: jj
        INTEGER(c_intk) :: kk
        INTEGER(c_intk) :: ip3
        INTEGER(c_intk) :: ipx
        INTEGER(c_intk) :: ipy
        INTEGER(c_intk) :: ipz
        INTEGER(c_intk) :: nfro
        INTEGER(c_intk) :: nbac
        INTEGER(c_intk) :: nrgt
        INTEGER(c_intk) :: nlft
        INTEGER(c_intk) :: nbot
        INTEGER(c_intk) :: ntop
    END TYPE lesmodel_grid_t

    ! LES model constant
    REAL(realk) :: Cm
    !$omp declare target(Cm)

    ! Bound operation
    TYPE(boundg_t) :: bound

    INTEGER(intk), PARAMETER :: boundgtasksize = 3
    INTEGER(intk), ALLOCATABLE :: boundgtasks(:, :, :)
    INTEGER(intk), ALLOCATABLE :: nboundgtaskslvl(:)
    TYPE(lesmodel_grid_t), ALLOCATABLE :: lesmodel_grids(:)
    !$omp declare target(boundgtasks, lesmodel_grids)

    INTERFACE
        SUBROUTINE lesmodel_gc_impl_c(g, u, v, w, bp, dx, dy, dz, ddx, ddy, &
            ddz, rddx, rddy, rddz, ilesmodel_in, cm_in, gmol_in, rho_in, &
            nmygrids_in, grids) BIND(C, name="lesmodel_gc_impl_c")
            IMPORT :: c_intk, c_realk, lesmodel_grid_t
            REAL(c_realk), INTENT(inout) :: g(*)
            REAL(c_realk), INTENT(in) :: u(*), v(*), w(*), bp(*)
            REAL(c_realk), INTENT(in) :: dx(*), dy(*), dz(*)
            REAL(c_realk), INTENT(in) :: ddx(*), ddy(*), ddz(*)
            REAL(c_realk), INTENT(in) :: rddx(*), rddy(*), rddz(*)
            INTEGER(c_intk), VALUE, INTENT(in) :: ilesmodel_in
            REAL(c_realk), VALUE, INTENT(in) :: cm_in, gmol_in, rho_in
            INTEGER(c_intk), VALUE, INTENT(in) :: nmygrids_in
            TYPE(lesmodel_grid_t), INTENT(in) :: grids(*)
        END SUBROUTINE lesmodel_gc_impl_c
    END INTERFACE

    !$omp declare target(lesmodel_gc_impl_c)

    PUBLIC :: init_lesmodel, finish_lesmodel, lesmodel, ilesmodel

CONTAINS

    SUBROUTINE init_lesmodel()
        ! Subroutine arguments
        ! none...

        ! Local variables
        INTEGER(intk) :: ilevel
        TYPE(config_t) :: lesconf
        TYPE(field_t), POINTER :: g

        CALL fort7%get(lesconf, "/flow/lesmodel")
        CALL lesconf%get_value("/model", clesmodel, "smagorinsky")

        SELECT CASE (lower(clesmodel))
        CASE ("none")
            ilesmodel = 0
            Cm = 0.0
        CASE ("smagorinsky")
            ilesmodel = 1
            Cm = 0.1
        CASE ("wale")
            ilesmodel = 2
            Cm = 0.5
        CASE("sigma")
#ifdef _MGLET_OFFLOAD_
            WRITE(*, *) "LES model 'sigma' is not supported with offloading"
            CALL errr(__FILE__, __LINE__)
#endif
            ilesmodel = 5
            Cm = 1.35
        CASE DEFAULT
            WRITE(*, *) "Invalid LES model:", clesmodel
            CALL errr(__FILE__, __LINE__)
        END SELECT

        ! Override default model parameter
        CALL lesconf%get_value("/Cm", Cm, Cm)
        !$omp target update to(ilesmodel, Cm)

        CALL init_boundg()
        CALL init_lesmodel_gc()

        ! Compute viscosity corresponding to initial condition
        CALL get_field(g, "G")
        IF (ilesmodel == 0) THEN
            g%arr = gmol
            CALL zero_ghostlayers(g)

            DO ilevel = minlevel, maxlevel
                CALL connect(ilevel, 2, s1=g, corners=.TRUE.)
            END DO

            DO ilevel = minlevel+1, maxlevel
                CALL parent(ilevel, s1=g)
                CALL bound%bound(ilevel, g)
            END DO

            CALL map_arr_to_device(g, message="to:init_lesmodel")
        ELSE
            CALL map_arr_to_device(g, message="to:init_lesmodel")
            CALL lesmodel(g)
        END IF
    END SUBROUTINE init_lesmodel


    SUBROUTINE finish_lesmodel()
        CALL finish_lesmodel_gc()
        !$omp target exit data map(always, delete: boundgtasks)
        DEALLOCATE(boundgtasks)
        DEALLOCATE(nboundgtaskslvl)
    END SUBROUTINE finish_lesmodel


    SUBROUTINE lesmodel(g)
        ! Subroutine arguments
        TYPE(field_t), INTENT(inout) :: g

        ! Local variables
        ! none...

        ! In case of no LESMODEL
        ! IF (ilesmodel == 0) RETURN
        CALL start_timer(330)
        CALL lesmodel_gc(g)
        CALL stop_timer(330)
    END SUBROUTINE lesmodel


    SUBROUTINE lesmodel_gc(g_f)

        ! Subroutine arguments
        TYPE(field_t), INTENT(inout) :: g_f

        ! Local variables
        TYPE(field_t), POINTER :: u_f, v_f, w_f, bp_f
        TYPE(field_t), POINTER :: dx_f, dy_f, dz_f
        TYPE(field_t), POINTER :: ddx_f, ddy_f, ddz_f
        TYPE(field_t), POINTER :: rddx_f, rddy_f, rddz_f

        INTEGER(intk) :: ilevel

        CALL get_field(u_f, "U")
        CALL get_field(v_f, "V")
        CALL get_field(w_f, "W")
        CALL get_field(bp_f, "BP")

        CALL get_field(dx_f, "DX")
        CALL get_field(dy_f, "DY")
        CALL get_field(dz_f, "DZ")

        CALL get_field(ddx_f, "DDX")
        CALL get_field(ddy_f, "DDY")
        CALL get_field(ddz_f, "DDZ")

        CALL get_field(rddx_f, "RDDX")
        CALL get_field(rddy_f, "RDDY")
        CALL get_field(rddz_f, "RDDZ")

        CALL lesmodel_gc_impl_c(g_f%arr, u_f%arr, v_f%arr, w_f%arr, bp_f%arr, &
            dx_f%arr, dy_f%arr, dz_f%arr, ddx_f%arr, ddy_f%arr, ddz_f%arr, &
            rddx_f%arr, rddy_f%arr, rddz_f%arr, INT(ilesmodel, c_intk), &
            REAL(Cm, c_realk), REAL(gmol, c_realk), REAL(rho, c_realk), &
            INT(nmygrids, c_intk), lesmodel_grids)

        DO ilevel = minlevel, maxlevel
            CALL parent(ilevel, s1=g_f, device=.TRUE.)
            ! CALL apply_boundg(ilevel, g_f, bp_f)   !! <<< Fucking SHiat
            CALL conn(ilevel, 1, s1=g_f)
        END DO

        CALL setginbody_impl(g_f%arr, bp_f%arr)

        ! TSTLE4 access corner values of viscosity such as (k, j+1, i+1),
        ! therefore connect with corners
        CALL conn(layers=1, s1=g_f, corners=.TRUE.)

        ! !$omp target update from(g_f%arr, bp_f%arr, u_f%arr, v_f%arr, w_f%arr)
        ! WRITE(*, *) minval(g_f%arr), maxval(g_f%arr)
        ! WRITE(*, *) minval(bp_f%arr), maxval(bp_f%arr)
        ! WRITE(*, *) minval(u_f%arr), maxval(u_f%arr)
        ! WRITE(*, *) minval(v_f%arr), maxval(v_f%arr)
        ! WRITE(*, *) minval(w_f%arr), maxval(w_f%arr)

    END SUBROUTINE lesmodel_gc


    SUBROUTINE setginbody_impl(g, bp)
        REAL(realk), INTENT(inout) :: g(*)
        REAL(realk), INTENT(in) :: bp(*)

        INTEGER(intk) :: i, igrid, ip3, kk, jj, ii

        !$omp target teams distribute private(i, igrid, ip3, kk, jj, ii)
        DO i = 1, nmygrids
            igrid = mygrids(i)
            CALL get_mgdims(kk, jj, ii, igrid)
            CALL get_ip3(ip3, igrid)
            !$omp parallel
            CALL setginbody(kk, jj, ii, bp(ip3), g(ip3))
            !$omp end parallel
        END DO
        !$omp end target teams distribute
    END SUBROUTINE setginbody_impl


    SUBROUTINE init_lesmodel_gc()
        INTEGER(intk) :: i

        ALLOCATE(lesmodel_grids(nmygrids))
        DO i = 1, nmygrids
            CALL get_mgdims(lesmodel_grids(i)%kk, lesmodel_grids(i)%jj, &
                lesmodel_grids(i)%ii, mygrids(i))
            CALL get_mgbasb(lesmodel_grids(i)%nfro, lesmodel_grids(i)%nbac, &
                lesmodel_grids(i)%nrgt, lesmodel_grids(i)%nlft, &
                lesmodel_grids(i)%nbot, lesmodel_grids(i)%ntop, mygrids(i))
            CALL get_ip3(lesmodel_grids(i)%ip3, mygrids(i))
            CALL get_ip1x(lesmodel_grids(i)%ipx, mygrids(i))
            CALL get_ip1y(lesmodel_grids(i)%ipy, mygrids(i))
            CALL get_ip1z(lesmodel_grids(i)%ipz, mygrids(i))
        END DO

        !$omp target enter data map(to: lesmodel_grids)
    END SUBROUTINE init_lesmodel_gc


    SUBROUTINE finish_lesmodel_gc()
        !$omp target exit data map(delete: lesmodel_grids)
        DEALLOCATE(lesmodel_grids)
    END SUBROUTINE finish_lesmodel_gc


    SUBROUTINE init_boundg()
        INTEGER(intk) :: nlevels, ilevel, ilevel_index
        INTEGER(intk) :: i, igrid, iface, ibocd, nbocd, itask
        CHARACTER(len=8) :: ctyp

        nlevels = maxlevel - minlevel + 1
        ALLOCATE(nboundgtaskslvl(nlevels), source=0_intk)

        DO ilevel = minlevel, maxlevel
            ilevel_index = ilevel - minlevel + 1
            DO i = 1, nmygridslvl(ilevel)
                igrid = mygridslvl(i, ilevel)
                DO iface = 1, 6
                    nbocd = nboconds(iface, igrid)
                    DO ibocd = 1, nbocd
                        CALL get_bc_ctyp(ctyp, ibocd, iface, igrid)
                        IF (boundg_ctyp(ctyp) == 0) CYCLE
                        nboundgtaskslvl(ilevel_index) = &
                            nboundgtaskslvl(ilevel_index) + 1
                    END DO
                END DO
            END DO
        END DO

        ALLOCATE(boundgtasks(boundgtasksize, MAXVAL(nboundgtaskslvl), &
            nlevels), source=-1_intk)

        DO ilevel = minlevel, maxlevel
            ilevel_index = ilevel - minlevel + 1
            itask = 0
            DO i = 1, nmygridslvl(ilevel)
                igrid = mygridslvl(i, ilevel)
                DO iface = 1, 6
                    nbocd = nboconds(iface, igrid)
                    DO ibocd = 1, nbocd
                        CALL get_bc_ctyp(ctyp, ibocd, iface, igrid)
                        IF (boundg_ctyp(ctyp) == 0) CYCLE
                        itask = itask + 1
                        boundgtasks(:, itask, ilevel_index) = &
                            [igrid, iface, boundg_ctyp(ctyp)]
                    END DO
                END DO
            END DO
        END DO

        !$omp target enter data map(always, to: boundgtasks)
    END SUBROUTINE init_boundg


    INTEGER(intk) FUNCTION boundg_ctyp(ctyp)
        CHARACTER(len=*), INTENT(in) :: ctyp

        SELECT CASE (ctyp)
        CASE ("NOS", "SLI")
            boundg_ctyp = 1
        CASE ("FIX", "OP1")
            boundg_ctyp = 2
        CASE ("PAR")
            boundg_ctyp = 3
        CASE DEFAULT
            boundg_ctyp = 0
        END SELECT
    END FUNCTION boundg_ctyp


    SUBROUTINE apply_boundg(ilevel, g_f, bp_f)
        INTEGER(intk), INTENT(in) :: ilevel
        TYPE(field_t), INTENT(inout) :: g_f
        TYPE(field_t), INTENT(in) :: bp_f

        INTEGER(intk) :: ilevel_index, ntasks

        ilevel_index = ilevel - minlevel + 1
        ntasks = nboundgtaskslvl(ilevel_index)

        CALL apply_boundg_impl(ilevel_index, ntasks, g_f%arr, g_f%buffers, &
            bp_f%arr)
    END SUBROUTINE apply_boundg


    SUBROUTINE apply_boundg_impl(ilevel_index, ntasks, g, gbuffer, bp)
        INTEGER(intk), INTENT(in) :: ilevel_index, ntasks
        REAL(realk), INTENT(inout) :: g(*), gbuffer(*)
        REAL(realk), INTENT(in) :: bp(*)

        INTEGER(intk) :: itask, igrid, iface, ityp, kk, jj, ii, ip3, ipbb

        !$omp target teams distribute private(itask, igrid, iface, ityp, &
        !$omp& kk, jj, ii, ip3, ipbb)
        DO itask = 1, ntasks
            igrid = boundgtasks(1, itask, ilevel_index)
            iface = boundgtasks(2, itask, ilevel_index)
            ityp = boundgtasks(3, itask, ilevel_index)

            CALL get_mgdims(kk, jj, ii, igrid)
            CALL get_ip3(ip3, igrid)
            IF (ityp /= 3) THEN
                !$omp parallel
                CALL boundg_nobuffer_device(kk, jj, ii, iface, ityp, g(ip3))
                !$omp end parallel
                CYCLE
            END IF
            CALL get_ipbb(ipbb, iface, igrid)

            !$omp parallel
            SELECT CASE (iface)
            CASE (1)
                CALL bfront_device(kk, jj, ii, 2, 3, ityp, &
                    gbuffer(ipbb), g(ip3), bp(ip3))
            CASE (2)
                CALL bfront_device(kk, jj, ii, ii-1, ii-2, ityp, &
                    gbuffer(ipbb), g(ip3), bp(ip3))
            CASE (3)
                CALL bright_device(kk, jj, ii, 2, 3, ityp, &
                    gbuffer(ipbb), g(ip3), bp(ip3))
            CASE (4)
                CALL bright_device(kk, jj, ii, jj-1, jj-2, ityp, &
                    gbuffer(ipbb), g(ip3), bp(ip3))
            CASE (5)
                CALL bbottom_device(kk, jj, ii, 2, 3, ityp, &
                    gbuffer(ipbb), g(ip3), bp(ip3))
            CASE (6)
                CALL bbottom_device(kk, jj, ii, kk-1, kk-2, ityp, &
                    gbuffer(ipbb), g(ip3), bp(ip3))
            END SELECT
            !$omp end parallel
        END DO
        !$omp end target teams distribute
    END SUBROUTINE apply_boundg_impl


    SUBROUTINE setginbody(kk, jj, ii, bp, g)
        !$omp declare target
        ! Subroutine arguments
        INTEGER(intk), INTENT(in) :: kk, jj, ii
        REAL(realk), INTENT(in) :: bp(kk, jj, ii)
        REAL(realk), INTENT(inout) :: g(kk, jj, ii)

        ! Local variables
        INTEGER(intk) :: k, j, i
        REAL(realk) :: gbpn, nn

        !$omp do collapse(3) private(i, j, k, gbpn, nn)
        DO i = 3, ii-2
            DO j = 3, jj-2
                DO k = 3, kk-2
                    IF (bp(k, j, i) > 0.5_realk) CYCLE
                    nn = 0.0_realk
                    gbpn = 0.0_realk
                    IF (bp(k, j, i+1) > 0.5_realk) THEN
                        nn = nn + 1.0_realk
                        gbpn = gbpn + g(k, j, i+1)
                    END IF
                    IF (bp(k, j, i-1) > 0.5_realk) THEN
                        nn = nn + 1.0_realk
                        gbpn = gbpn + g(k, j, i-1)
                    END IF
                    IF (bp(k, j+1, i) > 0.5_realk) THEN
                        nn = nn + 1.0_realk
                        gbpn = gbpn + g(k, j+1, i)
                    END IF
                    IF (bp(k, j-1, i) > 0.5_realk) THEN
                        nn = nn + 1.0_realk
                        gbpn = gbpn + g(k, j-1, i)
                    END IF
                    IF (bp(k+1, j, i) > 0.5_realk) THEN
                        nn = nn + 1.0_realk
                        gbpn = gbpn + g(k+1, j, i)
                    END IF
                    IF (bp(k-1, j, i) > 0.5_realk) THEN
                        nn = nn + 1.0_realk
                        gbpn = gbpn + g(k-1, j, i)
                    END IF
                    g(k, j, i) = gbpn/MAX(nn, 1.0_realk)
                END DO
            END DO
        END DO
        !$omp end do
    END SUBROUTINE setginbody


    SUBROUTINE boundg_nobuffer_device(kk, jj, ii, iface, ityp, g)
        !$omp declare target
        INTEGER(intk), INTENT(in) :: kk, jj, ii, iface, ityp
        REAL(realk), INTENT(inout) :: g(kk, jj, ii)

        INTEGER(intk) :: k, j, i, i2, i3, j2, j3, k2, k3

        SELECT CASE (iface)
        CASE (1, 2)
            IF (iface == 1) THEN
                i2 = 2
                i3 = 3
            ELSE
                i2 = ii - 1
                i3 = ii - 2
            END IF
            !$omp do collapse(2) private(j, k)
            DO j = 1, jj
                DO k = 1, kk
                    IF (ityp == 1) THEN
                        g(k, j, i2) = -EPSILON(1.0_realk)*gmol
                    ELSE
                        g(k, j, i2) = g(k, j, i3)
                    END IF
                END DO
            END DO
            !$omp end do
        CASE (3, 4)
            IF (iface == 3) THEN
                j2 = 2
                j3 = 3
            ELSE
                j2 = jj - 1
                j3 = jj - 2
            END IF
            !$omp do collapse(2) private(i, k)
            DO i = 1, ii
                DO k = 1, kk
                    IF (ityp == 1) THEN
                        g(k, j2, i) = -EPSILON(1.0_realk)*gmol
                    ELSE
                        g(k, j2, i) = g(k, j3, i)
                    END IF
                END DO
            END DO
            !$omp end do
        CASE (5, 6)
            IF (iface == 5) THEN
                k2 = 2
                k3 = 3
            ELSE
                k2 = kk - 1
                k3 = kk - 2
            END IF
            !$omp do collapse(2) private(i, j)
            DO i = 1, ii
                DO j = 1, jj
                    IF (ityp == 1) THEN
                        g(k2, j, i) = -EPSILON(1.0_realk)*gmol
                    ELSE
                        g(k2, j, i) = g(k3, j, i)
                    END IF
                END DO
            END DO
            !$omp end do
        END SELECT
    END SUBROUTINE boundg_nobuffer_device


    SUBROUTINE bfront_device(kk, jj, ii, i2, i3, ityp, buffer, g, bp)
        !$omp declare target
        INTEGER(intk), INTENT(in) :: kk, jj, ii, i2, i3, ityp
        REAL(realk), INTENT(in) :: buffer(kk, jj)
        REAL(realk), INTENT(inout) :: g(kk, jj, ii)
        REAL(realk), INTENT(in) :: bp(kk, jj, ii)

        INTEGER(intk) :: k, j
        REAL(realk) :: sbp

        SELECT CASE (ityp)
        CASE (1)
            !$omp do collapse(2) private(j, k)
            DO j = 1, jj
                DO k = 1, kk
                    g(k, j, i2) = -EPSILON(1.0_realk)*gmol
                END DO
            END DO
            !$omp end do
        CASE (2)
            !$omp do collapse(2) private(j, k)
            DO j = 1, jj
                DO k = 1, kk
                    g(k, j, i2) = g(k, j, i3)
                END DO
            END DO
            !$omp end do
        CASE (3)
            !$omp do collapse(2) private(j, k, sbp)
            DO j = 2, jj-1
                DO k = 2, kk-1
                    sbp = bp(k, j, i2)
                    g(k, j, i2) = buffer(k, j)*sbp &
                        + (1.0-sbp)*g(k, j, i2)
                END DO
            END DO
            !$omp end do
        END SELECT
    END SUBROUTINE bfront_device


    SUBROUTINE bright_device(kk, jj, ii, j2, j3, ityp, buffer, g, bp)
        !$omp declare target
        INTEGER(intk), INTENT(in) :: kk, jj, ii, j2, j3, ityp
        REAL(realk), INTENT(in) :: buffer(kk, ii)
        REAL(realk), INTENT(inout) :: g(kk, jj, ii)
        REAL(realk), INTENT(in) :: bp(kk, jj, ii)

        INTEGER(intk) :: k, i
        REAL(realk) :: sbp

        SELECT CASE (ityp)
        CASE (1)
            !$omp do collapse(2) private(i, k)
            DO i = 1, ii
                DO k = 1, kk
                    g(k, j2, i) = -EPSILON(1.0_realk)*gmol
                END DO
            END DO
            !$omp end do
        CASE (2)
            !$omp do collapse(2) private(i, k)
            DO i = 1, ii
                DO k = 1, kk
                    g(k, j2, i) = g(k, j3, i)
                END DO
            END DO
            !$omp end do
        CASE (3)
            !$omp do collapse(2) private(i, k, sbp)
            DO i = 2, ii-1
                DO k = 2, kk-1
                    sbp = bp(k, j2, i)
                    g(k, j2, i) = buffer(k, i)*sbp &
                        + (1.0-sbp)*g(k, j2, i)
                END DO
            END DO
            !$omp end do
        END SELECT
    END SUBROUTINE bright_device


    SUBROUTINE bbottom_device(kk, jj, ii, k2, k3, ityp, buffer, g, bp)
        !$omp declare target
        INTEGER(intk), INTENT(in) :: kk, jj, ii, k2, k3, ityp
        REAL(realk), INTENT(in) :: buffer(jj, ii)
        REAL(realk), INTENT(inout) :: g(kk, jj, ii)
        REAL(realk), INTENT(in) :: bp(kk, jj, ii)

        INTEGER(intk) :: j, i
        REAL(realk) :: sbp

        SELECT CASE (ityp)
        CASE (1)
            !$omp do collapse(2) private(i, j)
            DO i = 1, ii
                DO j = 1, jj
                    g(k2, j, i) = -EPSILON(1.0_realk)*gmol
                END DO
            END DO
            !$omp end do
        CASE (2)
            !$omp do collapse(2) private(i, j)
            DO i = 1, ii
                DO j = 1, jj
                    g(k2, j, i) = g(k3, j, i)
                END DO
            END DO
            !$omp end do
        CASE (3)
            !$omp do collapse(2) private(i, j, sbp)
            DO i = 2, ii-1
                DO j = 2, jj-1
                    sbp = bp(k2, j, i)
                    g(k2, j, i) = buffer(j, i)*sbp &
                        + (1.0-sbp)*g(k2, j, i)
                END DO
            END DO
            !$omp end do
        END SELECT
    END SUBROUTINE bbottom_device


    SUBROUTINE bfront(igrid, iface, ibocd, ctyp, f1, f2, f3, f4, timeph)
        ! Subroutine arguments
        INTEGER(intk), INTENT(in) :: igrid, iface, ibocd
        CHARACTER(len=*), INTENT(in) :: ctyp
        TYPE(field_t), INTENT(inout) :: f1
        TYPE(field_t), INTENT(inout), OPTIONAL :: f2, f3, f4
        REAL(realk), INTENT(in), OPTIONAL :: timeph

        ! Local variables
        INTEGER(intk) :: k, j, i2, i3
        REAL(realk) :: sbp
        REAL(realk), POINTER, CONTIGUOUS :: g(:, :, :), buffer(:, :), &
            bp(:, :, :)

        ! Return early when no action is to be taken
        SELECT CASE (ctyp)
        CASE ("FIX", "OP1", "NOS", "SLI", "PAR")
            CONTINUE
        CASE DEFAULT
            RETURN
        END SELECT

        CALL f1%get_ptr(g, igrid)

        SELECT CASE (iface)
        CASE (1)
            ! Front
            i2 = 2
            i3 = 3
        CASE (2)
            ! Back
            i2 = SIZE(g, 3) - 1
            i3 = SIZE(g, 3) - 2
        CASE DEFAULT
            CALL errr(__FILE__, __LINE__)
        END SELECT

        SELECT CASE (ctyp)
        CASE ("NOS", "SLI")
            DO j = 1, SIZE(g, 2)
                DO k = 1, SIZE(g, 1)
                    g(k, j, i2) = -EPSILON(1.0_realk)*gmol
                END DO
            END DO
        CASE ("FIX", "OP1")
            DO j = 1, SIZE(g, 2)
                DO k = 1, SIZE(g, 1)
                    g(k, j, i2) = g(k, j, i3)
                END DO
            END DO
        CASE ("PAR")
            CALL get_fieldptr(bp, "BP", igrid)
            CALL f1%get_buffer(buffer, igrid, iface)
            DO j = 2, SIZE(g, 2)-1
                DO k = 2, SIZE(g, 1)-1
                    sbp = bp(k, j, i2)
                    g(k, j, i2) = buffer(k, j)*sbp &
                        + (1.0-sbp)*g(k, j, i2)
                END DO
            END DO
        END SELECT
    END SUBROUTINE bfront


    SUBROUTINE bright(igrid, iface, ibocd, ctyp, f1, f2, f3, f4, timeph)
        ! Subroutine arguments
        INTEGER(intk), INTENT(in) :: igrid, iface, ibocd
        CHARACTER(len=*), INTENT(in) :: ctyp
        TYPE(field_t), INTENT(inout) :: f1
        TYPE(field_t), INTENT(inout), OPTIONAL :: f2, f3, f4
        REAL(realk), INTENT(in), OPTIONAL :: timeph

        ! Local variables
        INTEGER(intk) :: k, i, j2, j3
        REAL(realk) :: sbp
        REAL(realk), POINTER, CONTIGUOUS :: g(:, :, :), buffer(:, :), &
            bp(:, :, :)

        ! Return early when no action is to be taken
        SELECT CASE (ctyp)
        CASE ("FIX", "OP1", "NOS", "SLI", "PAR")
            CONTINUE
        CASE DEFAULT
            RETURN
        END SELECT

        CALL f1%get_ptr(g, igrid)

        SELECT CASE (iface)
        CASE (3)
            ! Right
            j2 = 2
            j3 = 3
        CASE (4)
            ! Left
            j2 = SIZE(g, 2) - 1
            j3 = SIZE(g, 2) - 2
        CASE DEFAULT
            CALL errr(__FILE__, __LINE__)
        END SELECT

        SELECT CASE (ctyp)
        CASE ("NOS", "SLI")
            DO i = 1, SIZE(g, 3)
                DO k = 1, SIZE(g, 1)
                    g(k, j2, i) = -EPSILON(1.0_realk)*gmol
                END DO
            END DO
        CASE ("FIX", "OP1")
            DO i = 1, SIZE(g, 3)
                DO k = 1, SIZE(g, 1)
                    g(k, j2, i) = g(k, j3, i)
                END DO
            END DO
        CASE ("PAR")
            CALL get_fieldptr(bp, "BP", igrid)
            CALL f1%get_buffer(buffer, igrid, iface)
            DO i = 2, SIZE(g, 3)-1
                DO k = 2, SIZE(g, 1)-1
                    sbp = bp(k, j2, i)
                    g(k, j2, i) = buffer(k, i)*sbp &
                        + (1.0-sbp)*g(k, j2, i)
                END DO
            END DO
        END SELECT
    END SUBROUTINE bright


    SUBROUTINE bbottom(igrid, iface, ibocd, ctyp, f1, f2, f3, f4, timeph)
        ! Subroutine arguments
        INTEGER(intk), INTENT(in) :: igrid, iface, ibocd
        CHARACTER(len=*), INTENT(in) :: ctyp
        TYPE(field_t), INTENT(inout) :: f1
        TYPE(field_t), INTENT(inout), OPTIONAL :: f2, f3, f4
        REAL(realk), INTENT(in), OPTIONAL :: timeph

        ! Local variables
        INTEGER(intk) :: j, i, k2, k3
        REAL(realk) :: sbp
        REAL(realk), POINTER, CONTIGUOUS :: g(:, :, :), buffer(:, :), &
            bp(:, :, :)

        ! Return early when no action is to be taken
        SELECT CASE (ctyp)
        CASE ("FIX", "OP1", "NOS", "SLI", "PAR")
            CONTINUE
        CASE DEFAULT
            RETURN
        END SELECT

        CALL f1%get_ptr(g, igrid)

        SELECT CASE (iface)
        CASE (5)
            ! Bottom
            k2 = 2
            k3 = 3
        CASE (6)
            ! Top
            k2 = SIZE(g, 1) - 1
            k3 = SIZE(g, 1) - 2
        CASE DEFAULT
            CALL errr(__FILE__, __LINE__)
        END SELECT

        SELECT CASE (ctyp)
        CASE ("NOS", "SLI")
            DO i = 1, SIZE(g, 3)
                DO j = 1, SIZE(g, 2)
                    g(k2, j, i) = -EPSILON(1.0_realk)*gmol
                END DO
            END DO
        CASE ("FIX", "OP1")
            DO i = 1, SIZE(g, 3)
                DO j = 1, SIZE(g, 2)
                    g(k2, j, i) = g(k3, j, i)
                END DO
            END DO
        CASE ("PAR")
            CALL get_fieldptr(bp, "BP", igrid)
            CALL f1%get_buffer(buffer, igrid, iface)
            DO i = 2, SIZE(g, 3)-1
                DO j = 2, SIZE(g, 2)-1
                    sbp = bp(k2, j, i)
                    g(k2, j, i) = buffer(j, i)*sbp &
                        + (1.0-sbp)*g(k2, j, i)
                END DO
            END DO
        END SELECT
    END SUBROUTINE bbottom

END MODULE lesmodel_mod
