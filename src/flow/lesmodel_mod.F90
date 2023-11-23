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
        CALL lesconf%get_value("/model", clesmodel, "Smagorinsky")

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
        ELSE
            CALL lesmodel(g)
        END IF
    END SUBROUTINE init_lesmodel


    SUBROUTINE finish_lesmodel()
        ! Does nothing right now...
        CONTINUE
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

        REAL(realk), POINTER, CONTIGUOUS :: u(:, :, :), v(:, :, :), &
            w(:, :, :), bp(:, :, :), g(:, :, :)
        REAL(realk), POINTER, CONTIGUOUS :: dx(:), dy(:), dz(:)
        REAL(realk), POINTER, CONTIGUOUS :: ddx(:), ddy(:), ddz(:)
        REAL(realk), POINTER, CONTIGUOUS :: rddx(:), rddy(:), rddz(:)

        INTEGER(intk) :: ilevel, i, igrid
        INTEGER(intk) :: kk, jj, ii
        INTEGER(intk) :: nfro, nbac, nrgt, nlft, nbot, ntop

        CALL get_field(u_f, "U")
        CALL get_field(v_f, "V")
        CALL get_field(w_f, "W")
        CALL get_field(bp_f, "BP")

        CALL get_field(dx_f, "DX")
        CALL get_field(dy_f, "DY")
        CALL get_field(dz_f, "DZ")

        CALL get_field(ddx_f, "DX")
        CALL get_field(ddy_f, "DY")
        CALL get_field(ddz_f, "DZ")

        CALL get_field(rddx_f, "RDDX")
        CALL get_field(rddy_f, "RDDY")
        CALL get_field(rddz_f, "RDDZ")

        DO ilevel = maxlevel, minlevel, -1
            DO i = 1, nmygridslvl(ilevel)
                igrid = mygridslvl(i, ilevel)

                CALL get_mgdims(kk, jj, ii, igrid)
                CALL get_mgbasb(nfro, nbac, nrgt, nlft, nbot, ntop, igrid)

                CALL g_f%get_ptr(g, igrid)

                CALL u_f%get_ptr(u, igrid)
                CALL v_f%get_ptr(v, igrid)
                CALL w_f%get_ptr(w, igrid)
                CALL bp_f%get_ptr(bp, igrid)

                CALL dx_f%get_ptr(dx, igrid)
                CALL dy_f%get_ptr(dy, igrid)
                CALL dz_f%get_ptr(dz, igrid)

                CALL ddx_f%get_ptr(ddx, igrid)
                CALL ddy_f%get_ptr(ddy, igrid)
                CALL ddz_f%get_ptr(ddz, igrid)

                CALL rddx_f%get_ptr(rddx, igrid)
                CALL rddy_f%get_ptr(rddy, igrid)
                CALL rddz_f%get_ptr(rddz, igrid)

                CALL efvisc_gc(kk, jj, ii, nfro, nbac, nrgt, nlft, nbot, ntop, &
                    dx, dy, dz, ddx, ddy, ddz, rddx, rddy, rddz, &
                    u, v, w, gmol, rho, bp, g)
            END DO
        END DO

        DO ilevel = minlevel, maxlevel
            CALL parent(ilevel, s1=g_f)
            CALL bound%bound(ilevel, g_f)
            CALL connect(ilevel, 1, s1=g_f)
        END DO

        DO ilevel = minlevel, maxlevel
            DO i = 1, nmygridslvl(ilevel)
                igrid = mygridslvl(i, ilevel)
                CALL get_mgdims(kk, jj, ii, igrid)
                CALL g_f%get_ptr(g, igrid)
                CALL bp_f%get_ptr(bp, igrid)
                CALL setginbody(kk, jj, ii, bp, g)
            END DO
        END DO

        ! TSTLE4 access corner values of viscosity such as (k, j+1, i+1),
        ! therefore connect with corners
        CALL connect(layers=1, s1=g_f, corners=.TRUE.)
    END SUBROUTINE lesmodel_gc


    SUBROUTINE efvisc_gc(kk, jj, ii, nfro, nbac, nrgt, nlft, nbot, ntop, &
            dx, dy, dz, ddx, ddy, ddz, rddx, rddy, rddz, u, v, w, gmol, rho, &
            bp, g)

        ! Subroutine arguments
        INTEGER(intk), INTENT(in) :: kk, jj, ii
        INTEGER(intk), INTENT(in) :: nfro, nbac, nrgt, nlft, nbot, ntop
        REAL(realk), INTENT(in) :: dx(ii), dy(jj), dz(kk)
        REAL(realk), INTENT(in) :: ddx(ii), ddy(jj), ddz(kk)
        REAL(realk), INTENT(in) :: rddx(ii), rddy(jj), rddz(kk)
        REAL(realk), INTENT(in) :: u(kk, jj, ii), v(kk, jj, ii), w(kk, jj, ii)
        REAL(realk), INTENT(in) :: gmol, rho
        REAL(realk), INTENT(in) :: bp(kk, jj, ii)
        REAL(realk), INTENT(inout) :: g(kk, jj, ii)

        ! Local variables
        INTEGER(intk) :: k, j, i
        REAL(realk) :: dudx(kk-4), dudy(kk-4), dudz(kk-4), &
            dvdx(kk-4), dvdy(kk-4), dvdz(kk-4), &
            dwdx(kk-4), dwdy(kk-4), dwdz(kk-4)
        REAL(realk) :: dxf, dyf, dzf
        REAL(realk) :: dxf2, dyf2, dzf2
        REAL(realk) :: delta(kk-4), dm

        DO i = 3, ii-2
            DO j = 3, jj-2
                !$omp simd private(dxf, dyf, dzf)
                DO k = 3, kk-2
                    dxf = 0.25*dx(i-1)*rddx(i)
                    ! dU/dX
                    dudx(k-2) = (u(k, j, i) - u(k, j, i-1))*rddx(i)
                    ! dU/dY
                    dudy(k-2) = ((u(k, j+1, i) - u(k, j-1, i))*dxf &
                        + (u(k, j+1, i-1) - u(k, j-1, i-1))*(0.5-dxf))*rddy(j)
                    ! dU/dZ
                    dudz(k-2) = ((u(k+1, j, i) - u(k-1, j, i))*dxf &
                        + (u(k+1, j, i-1) - u(k-1, j, i-1))*(0.5-dxf))*rddz(k)

                    dyf = 0.25*dy(j-1)*rddy(j)
                    ! dV/dX
                    dvdx(k-2) = ((v(k, j, i+1) - v(k, j, i-1))*dyf &
                        + (v(k, j-1, i+1) - v(k, j-1, i-1))*(0.5-dyf))*rddx(i)
                    ! dV/dY
                    dvdy(k-2) = (v(k, j, i) - v(k, j-1, i))*rddy(j)
                    ! dV/dZ
                    dvdz(k-2) = ((v(k+1, j, i) - v(k-1, j, i))*dyf &
                        + (v(k+1, j-1, i) - v(k-1, j-1, i))*(0.5-dyf))*rddz(k)

                    dzf = 0.25*dz(k-1)*rddz(k)
                    ! dW/dX
                    dwdx(k-2) = ((w(k, j, i+1) - w(k, j, i-1))*dzf &
                        + (w(k-1, j, i+1) - w(k-1, j, i-1))*(0.5-dzf))*rddx(i)
                    ! dW/dY
                    dwdy(k-2) = ((w(k, j+1, i) - w(k, j-1, i))*dzf &
                        + (w(k-1, j+1, i) - w(k-1, j-1, i))*(0.5-dzf))*rddy(j)
                    ! dW/dZ
                    dwdz(k-2) = (w(k, j, i) - w(k-1, j, i))*rddz(k)
                END DO

                !$omp simd
                DO k = 3, kk-2
                    delta(k-2) = cube_root(ddx(i)*ddy(j)*ddz(k))
                    delta(k-2) = delta(k-2)*bp(k, j, i)
                END DO

                ! NOS corrections to gradients
                ! GRADP2 uses thw WW wall model to compute the gradient
                ! ddx/2 away from the wall
                IF (nfro == 5 .AND. i == 3) THEN
                    DO k = 3, kk-2
                        dyf2 = 0.5*dy(j-1)*rddy(j)
                        dvdx(k-2) = gradp2(dyf2*v(k, j, i) &
                            + (1.0-dyf2)*v(k, j-1, i), ddx(i))

                        dzf2 = 0.5*dz(k-1)*rddz(k)
                        dwdx(k-2) = gradp2(dzf2*w(k, j, i) &
                            + (1.0-dzf2)*w(k-1, j, i), ddx(i))
                    END DO
                END IF
                IF (nbac == 5 .AND. i == ii-2) THEN
                    DO k = 3, kk-2
                        dyf2 = 0.5*dy(j-1)*rddy(j)
                        dvdx(k-2) = -gradp2(dyf2*v(k, j, i) &
                            + (1.0-dyf2)*v(k, j-1, i), ddx(i))

                        dzf2 = 0.5*dz(k-1)*rddz(k)
                        dwdx(k-2) = -gradp2(dzf2*w(k, j, i) &
                            + (1.0-dzf2)*w(k-1, j, i), ddx(i))
                    END DO
                END IF
                IF (nrgt == 5 .AND. j == 3) THEN
                    DO k = 3, kk-2
                        dxf2 = 0.5*dx(i-1)*rddx(i)
                        dudy(k-2) = gradp2(dxf2*u(k, j, i) &
                            + (1.0-dxf2)*u(k, j, i-1), ddy(j))

                        dzf2 = 0.5*dz(k-1)*rddz(k)
                        dwdy(k-2) = gradp2(dzf2*w(k, j, i) &
                            + (1.0-dzf2)*w(k-1, j, i), ddy(j))
                    END DO
                END IF
                IF (nlft == 5 .AND. j == jj-2) THEN
                    DO k = 3, kk-2
                        dxf2 = 0.5*dx(i-1)*rddx(i)
                        dudy(k-2) = -gradp2(dxf2*u(k, j, i) &
                            + (1.0-dxf2)*u(k, j, i-1), ddy(j))

                        dzf2 = 0.5*dz(k-1)*rddz(k)
                        dwdy(k-2) = -gradp2(dzf2*w(k, j, i) &
                            + (1.0-dzf2)*w(k-1, j, i), ddy(j))
                    END DO
                END IF
                IF (nbot == 5) THEN
                    ! The programmer is lazy and the compiler intelligent...
                    DO k = 3, 3
                        dxf2 = 0.5*dx(i-1)*rddx(i)
                        dudz(k-2) = gradp2(dxf2*u(k, j, i) &
                            + (1.0-dxf2)*u(k, j, i-1), ddz(k))

                        dyf2 = 0.5*dy(j-1)*rddy(j)
                        dvdz(k-2) = gradp2(dyf2*v(k, j, i) &
                            + (1.0-dyf2)*v(k, j-1, i), ddz(k))
                    END DO
                END IF
                IF (ntop == 5) THEN
                    DO k = kk-2, kk-2
                        dxf2 = 0.5*dx(i-1)*rddx(i)
                        dudz(k-2) = -gradp2(dxf2*u(k, j, i) &
                            + (1.0-dxf2)*u(k, j, i-1), ddz(k))

                        dyf2 = 0.5*dy(j-1)*rddy(j)
                        dvdz(k-2) = -gradp2(dyf2*v(k, j, i) &
                            + (1.0-dyf2)*v(k, j-1, i), ddz(k))
                    END DO
                END IF

                SELECT CASE (ilesmodel)
                CASE (1)
                    !$omp simd private(dm)
                    DO k = 3, kk-2
                        dm = smagorinsky(dudx(k-2), dudy(k-2), dudz(k-2), &
                            dvdx(k-2), dvdy(k-2), dvdz(k-2), &
                            dwdx(k-2), dwdy(k-2), dwdz(k-2))
                        g(k, j, i) = rho*delta(k-2)**2*dm + gmol
                    END DO
                CASE (2)
                    DO k = 3, kk-2
                        dm = wale(dudx(k-2), dudy(k-2), dudz(k-2), &
                            dvdx(k-2), dvdy(k-2), dvdz(k-2), &
                            dwdx(k-2), dwdy(k-2), dwdz(k-2))
                        g(k, j, i) = rho*delta(k-2)**2*dm + gmol
                    END DO
                CASE (5)
                    DO k = 3, kk-2
                        dm = sigma(dudx(k-2), dudy(k-2), dudz(k-2), &
                            dvdx(k-2), dvdy(k-2), dvdz(k-2), &
                            dwdx(k-2), dwdy(k-2), dwdz(k-2))
                        g(k, j, i) = rho*delta(k-2)**2*dm + gmol
                    END DO
                END SELECT
            END DO
        END DO
    END SUBROUTINE efvisc_gc


    PURE ELEMENTAL REAL(realk) FUNCTION smagorinsky(dudx, dudy, dudz, dvdx, &
    dvdy, dvdz, dwdx, dwdy, dwdz)
        !$omp declare simd(smagorinsky)

        ! Function arguments
        REAL(realk), INTENT(in) :: dudx, dudy, dudz, dvdx, &
            dvdy, dvdz, dwdx, dwdy, dwdz

        ! Local variables
        REAL(realk) :: S_abs
        REAL(realk), PARAMETER :: root2 = SQRT(2.0_realk)

        ! S_ij = 0.5*(du_i/dx_j + du_j/dx_i)
        ! sabs = |S_ij| = sqrt(Sij*Sij)
        S_abs = sabs(dudx, dudy, dudz, &
            dvdx, dvdy, dvdz, dwdx, dwdy, dwdz)
        smagorinsky = Cm**2*root2*S_abs
    END FUNCTION smagorinsky


    PURE ELEMENTAL REAL(realk) FUNCTION wale(dudx, dudy, dudz, dvdx, &
            dvdy, dvdz, dwdx, dwdy, dwdz)

        ! Function arguments
        REAL(realk), INTENT(in) :: dudx, dudy, dudz, dvdx, &
            dvdy, dvdz, dwdx, dwdy, dwdz

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
        Sd_abs = SQRT(Sd00**2 + 2*Sd01**2 + 2*Sd02**2 \
            + Sd11**2 + 2*Sd12**2 + Sd22**2)

#if 0
        ! Alternative and slightly less performant implementation
        gij%g = RESHAPE((/dudx, dudy, dudz, dvdx, dvdy, dvdz, &
            dwdx, dwdy, dwdz/), (/3, 3/))

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
            wale = Cm**2*Sd_abs_root/ &
                (divide0(S_abs, Sd_abs_root)**5 + 1.0)
        END IF

    END FUNCTION wale


    PURE ELEMENTAL REAL(realk) FUNCTION sigma(dudx, dudy, dudz, &
            dvdx, dvdy, dvdz, dwdx, dwdy, dwdz)

        ! Function arguments
        REAL(realk), INTENT(in) :: dudx, dudy, dudz, dvdx, &
            dvdy, dvdz, dwdx, dwdy, dwdz

        ! Local variables
        TYPE(tensor_t) :: gij, G
        REAL(realk) :: eig1, eig2, eig3
        REAL(realk) :: sigma1, sigma2, sigma3

        gij%g = RESHAPE((/dudx, dudy, dudz, dvdx, dvdy, dvdz, &
            dwdx, dwdy, dwdz/), (/3, 3/))

        G = gij%t()*gij
        CALL G%eig_b(eig1, eig2, eig3)

        ! Sometimes I see eigenvalues such as -4.19095159e-09
        eig1 = MAX(eig1, 0.0_realk)
        eig2 = MAX(eig2, 0.0_realk)
        eig3 = MAX(eig3, 0.0_realk)

        sigma1 = SQRT(eig1)
        sigma2 = SQRT(eig2)
        sigma3 = SQRT(eig3)

        sigma = Cm**2*divide0( &
            sigma3*(sigma1-sigma2)*(sigma2-sigma3), sigma1**2)
    END FUNCTION sigma


    PURE ELEMENTAL REAL(realk) FUNCTION sabs(dudx, dudy, dudz, dvdx, &
            dvdy, dvdz, dwdx, dwdy, dwdz)
        !$omp declare simd(sabs)

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


    PURE SUBROUTINE setginbody(kk, jj, ii, bp, g)
        ! Subroutine arguments
        INTEGER(intk), INTENT(in) :: kk, jj, ii
        REAL(realk), INTENT(in) :: bp(kk, jj, ii)
        REAL(realk), INTENT(inout) :: g(kk, jj, ii)

        ! Local variables
        INTEGER(intk) :: k, j, i
        REAL(realk) :: gbpn, nn

        DO i = 3, ii-2
            DO j = 3, jj-2
                DO k = 3, kk-2
                    nn = bp(k, j, i+1) + bp(k, j, i-1) + bp(k, j+1, i) + &
                        bp(k, j-1, i) + bp(k+1, j, i) + bp(k-1, j, i)
                    nn = MAX(nn, 1.0_realk)
                    gbpn = (g(k, j, i+1)*bp(k, j, i+1) &
                        + g(k, j, i-1)*bp(k, j, i-1) &
                        + g(k, j+1, i)*bp(k, j+1, i) &
                        + g(k, j-1, i)*bp(k, j-1, i) &
                        + g(k+1, j, i)*bp(k+1, j, i) &
                        + g(k-1, j, i)*bp(k-1, j, i))/nn
                    g(k, j, i) = bp(k, j, i)*g(k, j, i) + &
                        (1.0 - bp(k, j, i))*gbpn
                END DO
            END DO
        END DO
    END SUBROUTINE setginbody


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
        REAL(realk), POINTER, CONTIGUOUS :: g(:, :, :), buffer(:, :, :), &
            bp(:, :, :)

        ! Return early when no action is to be taken
        SELECT CASE (ctyp)
        CASE ("FIX", "OP1", "NOS", "SLI", "PAR")
            CONTINUE
        CASE DEFAULT
            RETURN
        END SELECT

        CALL f1%get_ptr(g, igrid)
        CALL f1%buffers%get_buffer(buffer, igrid, iface)
        CALL get_fieldptr(bp, "BP", igrid)

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
            DO j = 2, SIZE(g, 2)-1
                DO k = 2, SIZE(g, 1)-1
                    sbp = bp(k, j, i2)
                    g(k, j, i2) = buffer(k, j, 1)*sbp &
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
        REAL(realk), POINTER, CONTIGUOUS :: g(:, :, :), buffer(:, :, :), &
            bp(:, :, :)

        ! Return early when no action is to be taken
        SELECT CASE (ctyp)
        CASE ("FIX", "OP1", "NOS", "SLI", "PAR")
            CONTINUE
        CASE DEFAULT
            RETURN
        END SELECT

        CALL f1%get_ptr(g, igrid)
        CALL f1%buffers%get_buffer(buffer, igrid, iface)
        CALL get_fieldptr(bp, "BP", igrid)

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
            DO i = 2, SIZE(g, 3)-1
                DO k = 2, SIZE(g, 1)-1
                    sbp = bp(k, j2, i)
                    g(k, j2, i) = buffer(k, i, 1)*sbp &
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
        REAL(realk), POINTER, CONTIGUOUS :: g(:, :, :), buffer(:, :, :), &
            bp(:, :, :)

        ! Return early when no action is to be taken
        SELECT CASE (ctyp)
        CASE ("FIX", "OP1", "NOS", "SLI", "PAR")
            CONTINUE
        CASE DEFAULT
            RETURN
        END SELECT

        CALL f1%get_ptr(g, igrid)
        CALL f1%buffers%get_buffer(buffer, igrid, iface)
        CALL get_fieldptr(bp, "BP", igrid)

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
            DO i = 2, SIZE(g, 3)-1
                DO j = 2, SIZE(g, 2)-1
                    sbp = bp(k2, j, i)
                    g(k2, j, i) = buffer(j, i, 1)*sbp &
                        + (1.0-sbp)*g(k2, j, i)
                END DO
            END DO
        END SELECT
    END SUBROUTINE bbottom

END MODULE lesmodel_mod
