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

    TYPE, EXTENDS(bound_t) :: boundg_t
    CONTAINS
        PROCEDURE, NOPASS :: front => bfront
        PROCEDURE, NOPASS :: back => bfront
        PROCEDURE, NOPASS :: right => bright
        PROCEDURE, NOPASS :: left => bright
        PROCEDURE, NOPASS :: bottom => bbottom
        PROCEDURE, NOPASS :: top => bbottom
    END TYPE boundg_t

    ! LES model constant
    REAL(realk) :: Cm

    ! Bound operation
    TYPE(boundg_t) :: bound

    INTEGER(intk), PARAMETER :: boundgtasksize = 3
    INTEGER(intk), ALLOCATABLE :: boundgtasks(:, :, :)
    INTEGER(intk), ALLOCATABLE :: nboundgtaskslvl(:)
    !$omp declare target(boundgtasks)

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
            ilesmodel = 5
            Cm = 1.35
        CASE DEFAULT
            WRITE(*, *) "Invalid LES model:", clesmodel
            CALL errr(__FILE__, __LINE__)
        END SELECT

        ! Override default model parameter
        CALL lesconf%get_value("/Cm", Cm, Cm)

        CALL init_boundg()

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
        IF (ilesmodel == 0) RETURN
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

        CALL lesmodel_gc_impl(g_f%arr, u_f%arr, v_f%arr, w_f%arr, bp_f%arr, &
            dx_f%arr, dy_f%arr, dz_f%arr, ddx_f%arr, ddy_f%arr, ddz_f%arr, &
            rddx_f%arr, rddy_f%arr, rddz_f%arr)

        DO ilevel = minlevel, maxlevel
            CALL parent(ilevel, s1=g_f, device=.TRUE.)
            CALL apply_boundg(ilevel, g_f, bp_f)
            CALL conn(ilevel, 1, s1=g_f)
        END DO

        CALL setginbody_impl(g_f%arr, bp_f%arr)

        ! TSTLE4 access corner values of viscosity such as (k, j+1, i+1),
        ! therefore connect with corners
        CALL conn(layers=1, s1=g_f, corners=.TRUE.)
    END SUBROUTINE lesmodel_gc


    SUBROUTINE lesmodel_gc_impl(g, u, v, w, bp, dx, dy, dz, ddx, ddy, ddz, &
            rddx, rddy, rddz)
        REAL(realk), INTENT(inout) :: g(*)
        REAL(realk), INTENT(in) :: u(*), v(*), w(*), bp(*)
        REAL(realk), INTENT(in) :: dx(*), dy(*), dz(*), ddx(*), ddy(*), ddz(*)
        REAL(realk), INTENT(in) :: rddx(*), rddy(*), rddz(*)

        SELECT CASE (ilesmodel)
        CASE (1)
            CALL lesmodel_gc_smagorinsky_impl(g, u, v, w, bp, dx, dy, dz, &
                ddx, ddy, ddz, rddx, rddy, rddz, Cm)
        CASE (2)
            CALL lesmodel_gc_wale_impl(g, u, v, w, bp, dx, dy, dz, ddx, &
                ddy, ddz, rddx, rddy, rddz, Cm)
        CASE (5)
            CALL lesmodel_gc_sigma_impl(g, u, v, w, bp, dx, dy, dz, ddx, &
                ddy, ddz, rddx, rddy, rddz, Cm)
        CASE DEFAULT
            CALL errr(__FILE__, __LINE__)
        END SELECT
    END SUBROUTINE lesmodel_gc_impl


#define LESMODEL_IMPL lesmodel_gc_smagorinsky_impl
#define LESMODEL_GRID efvisc_gc_smagorinsky
#define LESMODEL_EVALUATE smagorinsky
#include "lesmodel_kernel.inc"
#undef LESMODEL_EVALUATE
#undef LESMODEL_GRID
#undef LESMODEL_IMPL

#define LESMODEL_IMPL lesmodel_gc_wale_impl
#define LESMODEL_GRID efvisc_gc_wale
#define LESMODEL_EVALUATE wale
#include "lesmodel_kernel.inc"
#undef LESMODEL_EVALUATE
#undef LESMODEL_GRID
#undef LESMODEL_IMPL

#define LESMODEL_IMPL lesmodel_gc_sigma_impl
#define LESMODEL_GRID efvisc_gc_sigma
#define LESMODEL_EVALUATE sigma
#include "lesmodel_kernel.inc"
#undef LESMODEL_EVALUATE
#undef LESMODEL_GRID
#undef LESMODEL_IMPL


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


    PURE ELEMENTAL REAL(realk) FUNCTION smagorinsky(dudx, dudy, dudz, dvdx, &
            dvdy, dvdz, dwdx, dwdy, dwdz, model_constant)
        !$omp declare target
#ifndef _MGLET_OFFLOAD_
        !$omp declare simd(smagorinsky)
#endif

        ! Function arguments
        REAL(realk), INTENT(in) :: dudx, dudy, dudz, dvdx, &
            dvdy, dvdz, dwdx, dwdy, dwdz
        REAL(realk), INTENT(in) :: model_constant

        ! Local variables
        REAL(realk) :: S_abs
        REAL(realk), PARAMETER :: root2 = SQRT(2.0_realk)

        ! S_ij = 0.5*(du_i/dx_j + du_j/dx_i)
        ! sabs = |S_ij| = sqrt(Sij*Sij)
        S_abs = sabs(dudx, dudy, dudz, &
            dvdx, dvdy, dvdz, dwdx, dwdy, dwdz)
        smagorinsky = model_constant**2*root2*S_abs
    END FUNCTION smagorinsky


    PURE ELEMENTAL REAL(realk) FUNCTION wale(dudx, dudy, dudz, dvdx, &
            dvdy, dvdz, dwdx, dwdy, dwdz, model_constant)
        !$omp declare target

        ! Function arguments
        REAL(realk), INTENT(in) :: dudx, dudy, dudz, dvdx, &
            dvdy, dvdz, dwdx, dwdy, dwdz
        REAL(realk), INTENT(in) :: model_constant

        ! Local variables
        REAL(realk) :: tr
        REAL(realk) :: Sd00, Sd01, Sd02, Sd11, Sd12, Sd22
        REAL(realk) :: S_abs, Sd_abs, Sd_abs_root

        ! For the expression (SABS/absGaRoot)**5 we must ensure it does not
        ! overflow an REAL(realk)
        !
        ! Ideally the tolerance is exactly 1.0/(HUGE**(1./5.))
        ! (HUGE = 3.4e38 in single prec.), this lead to exactly
        ! (HUGE**(1./5.))**5
        !
        ! However, the overall expression it is used in is anyways 0 in this
        ! case... Therefore give some numerical "slack" and set the
        ! exponent to (1./6.)
        !
        ! With this tolerance you get HUGE**(5./6.) as the resulting cut-off.
        REAL(realk), PARAMETER :: sabs_tol = 1.0_realk/ &
            (huge(0.0_realk)**(1.0_realk/6.0_realk))

        ! Sd_ij = (1/2)*(g_ij**2 + g_ji**2) - (1/3)*delta_ij*g_kk**2
        ! The final result is a symmetric tensor, so only need to compute
        ! 6 out of 9 components
        tr = dudx*dudx + dudy*dvdx + dudz*dwdx + dvdx*dudy + dvdy*dvdy &
            + dvdz*dwdy + dwdx*dudz + dwdy*dvdz + dwdz*dwdz
        Sd00 = 0.5*(dudx*dudx + dudy*dvdx + dudz*dwdx + dudx*dudx &
            + dudy*dvdx + dudz*dwdx) - (1.0/3.0)*tr
        Sd01 = 0.5*(dudx*dudy + dudy*dvdy + dudz*dwdy + dvdx*dudx &
            + dvdy*dvdx + dvdz*dwdx)
        Sd02 = 0.5*(dudx*dudz + dudy*dvdz + dudz*dwdz + dwdx*dudx &
            + dwdy*dvdx + dwdz*dwdx)
        ! Sd10 = 0.5*(dvdx*dudx + dvdy*dvdx + dvdz*dwdx + dudx*dudy &
        !     + dudy*dvdy + dudz*dwdy)
        Sd11 = 0.5*(dvdx*dudy + dvdy*dvdy + dvdz*dwdy + dvdx*dudy &
            + dvdy*dvdy + dvdz*dwdy) - (1.0/3.0)*tr
        Sd12 = 0.5*(dvdx*dudz + dvdy*dvdz + dvdz*dwdz + dwdx*dudy &
            + dwdy*dvdy + dwdz*dwdy)
        ! Sd20 = 0.5*(dwdx*dudx + dwdy*dvdx + dwdz*dwdx &
        !     + dudx*dudz + dudy*dvdz + dudz*dwdz)
        ! Sd21 = 0.5*(dwdx*dudy + dwdy*dvdy + dwdz*dwdy &
        !     + dvdx*dudz + dvdy*dvdz + dvdz*dwdz)
        Sd22 = 0.5*(dwdx*dudz + dwdy*dvdz + dwdz*dwdz + dwdx*dudz &
            + dwdy*dvdz + dwdz*dwdz) - (1.0/3.0)*tr

        ! |Sd| = SQRT(Sd_ij*Sd_ij)
        ! We need to take the off-diagonal terms twice due to symmetry,
        ! i.e. Sd12**2 + Sd21**2 = 2*Sd12**2
        Sd_abs = SQRT(Sd00**2 + 2*Sd01**2 + 2*Sd02**2 &
            + Sd11**2 + 2*Sd12**2 + Sd22**2)

#if 0
        ! Alternative and slightly less performant implementation
        gij%g = RESHAPE([dudx, dudy, dudz, dvdx, dvdy, dvdz, &
            dwdx, dwdy, dwdz], [3, 3])

        ! gg = g_ik*g_kj
        gg = gij%sqr()

        ! Sd_ij = (1/2)*(g_ij**2 + g_ji**2) - (1/3)*delta_ij*g_kk**2
        Sd = 0.5*(gg + gg%t()) - (1.0/3.0)*gg%trace()*ident

        ! |Sd| = SQRT(Sd_ij*Sd_ij)
        Sd_abs = Sd%abst()
#endif

        ! As in Smagorinsky model
        ! |S| = SQRT(S_ij*S_ij)
        S_abs = sabs(dudx, dudy, dudz, dvdx, dvdy, dvdz, dwdx, dwdy, dwdz)

        ! Now here comes a dirty trick:
        !
        ! The naive formulation right from the Nicoud et al. paper is:
        !     D_w = Sd_abs**3.0/(S_abs**5.0 + Sd_abs**(5.0/2.0))
        ! This is not fesible in single prec. When the griadients become large
        ! the power of 5 become numerically troublesome since it overflows
        !
        ! A better formulation is to divide both sides of the fraction by
        ! Sd_abs**(5.0/2.0), then we get:
        !     D_w = SQRT(Sd_abs)/((S_abs/SQRT(Sd_abs))**5.0+1.0)
        ! Sd_abs_root = SQRT(Sd_abs) = (Sd_ij*Sd_ij)**(1.0/4.0) saves us one
        ! square root computation (if the compiler did not optimize it away
        ! already):
        Sd_abs_root = SQRT(Sd_abs)

        ! Final viscosity computation
        IF (Sd_abs_root < sabs_tol*S_abs) THEN
            wale = 0.0
        ELSE
            wale = model_constant**2*Sd_abs_root/ &
                (divide0(S_abs, Sd_abs_root)**5 + 1.0)
        END IF

    END FUNCTION wale


    PURE ELEMENTAL REAL(realk) FUNCTION sigma(dudx, dudy, dudz, &
            dvdx, dvdy, dvdz, dwdx, dwdy, dwdz, model_constant)
        !$omp declare target
        ! Function arguments
        REAL(realk), INTENT(in) :: dudx, dudy, dudz, dvdx, &
            dvdy, dvdz, dwdx, dwdy, dwdz
        REAL(realk), INTENT(in) :: model_constant

        ! Local variables
        TYPE(tensor_t) :: gij, G
        REAL(realk) :: eig1, eig2, eig3
        REAL(realk) :: sigma1, sigma2, sigma3

        gij%g = RESHAPE([dudx, dudy, dudz, dvdx, dvdy, dvdz, &
            dwdx, dwdy, dwdz], [3, 3])

        G = gij%t()*gij
        CALL G%eig_b(eig1, eig2, eig3)

        ! Sometimes I see eigenvalues such as -4.19095159e-09
        eig1 = MAX(eig1, 0.0_realk)
        eig2 = MAX(eig2, 0.0_realk)
        eig3 = MAX(eig3, 0.0_realk)

        sigma1 = SQRT(eig1)
        sigma2 = SQRT(eig2)
        sigma3 = SQRT(eig3)

        sigma = model_constant**2*divide0( &
            sigma3*(sigma1-sigma2)*(sigma2-sigma3), sigma1**2)
    END FUNCTION sigma


    PURE ELEMENTAL REAL(realk) FUNCTION sabs(dudx, dudy, dudz, dvdx, &
            dvdy, dvdz, dwdx, dwdy, dwdz)
        !$omp declare target
#ifndef _MGLET_OFFLOAD_
        !$omp declare simd(sabs)
#endif

        ! Function arguments
        REAL(realk), INTENT(in) :: dudx, dudy, dudz, dvdx, &
            dvdy, dvdz, dwdx, dwdy, dwdz

        ! S_ij = 0.5*(du_i/dx_j + du_j/dx_i)
        ! sabs = |S_ij| = sqrt(Sij*Sij)
        sabs = SQRT(dudx**2 + dvdy**2 + dwdz**2 &
            + 0.5*(dvdx + dudy)**2 &
            + 0.5*(dwdx + dudz)**2 &
            + 0.5*(dwdy + dvdz)**2)
    END FUNCTION sabs


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
