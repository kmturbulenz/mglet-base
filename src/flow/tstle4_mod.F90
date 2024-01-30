MODULE tstle4_mod
    USE core_mod
    USE flowcore_mod
    USE lesmodel_mod, ONLY: ilesmodel
    USE wernerwengle_mod, ONLY: tauwin

    IMPLICIT NONE(type, external)
    PRIVATE

    PUBLIC :: tstle4

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
        REAL(realk), POINTER, CONTIGUOUS, DIMENSION(:, :, :) :: uo, vo, wo
        REAL(realk), POINTER, CONTIGUOUS, DIMENSION(:, :, :) :: u, v, w
        REAL(realk), POINTER, CONTIGUOUS, DIMENSION(:, :, :) :: ut, vt, wt
        REAL(realk), POINTER, CONTIGUOUS, DIMENSION(:, :, :) :: p, g
        REAL(realk), POINTER, CONTIGUOUS :: dx(:), dy(:), dz(:)
        REAL(realk), POINTER, CONTIGUOUS :: ddx(:), ddy(:), ddz(:)
        REAL(realk), POINTER, CONTIGUOUS :: rdx(:), rdy(:), rdz(:)
        REAL(realk), POINTER, CONTIGUOUS :: rddx(:), rddy(:), rddz(:)
        INTEGER(intk) :: i, igrid
        INTEGER(intk) :: kk, jj, ii
        INTEGER(intk) :: nfro, nbac, nrgt, nlft, nbot, ntop

        CALL start_timer(310)

        ! Set all the output to zero everywhere before we start!
        uo_f = 0.0_realk
        vo_f = 0.0_realk
        wo_f = 0.0_realk

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

        DO i = 1, nmygrids
            igrid = mygrids(i)

            CALL get_mgdims(kk, jj, ii, igrid)
            CALL get_mgbasb(nfro, nbac, nrgt, nlft, nbot, ntop, igrid)

            CALL uo_f%get_ptr(uo, igrid)
            CALL vo_f%get_ptr(vo, igrid)
            CALL wo_f%get_ptr(wo, igrid)

            CALL u_f%get_ptr(u, igrid)
            CALL v_f%get_ptr(v, igrid)
            CALL w_f%get_ptr(w, igrid)

            CALL ut_f%get_ptr(ut, igrid)
            CALL vt_f%get_ptr(vt, igrid)
            CALL wt_f%get_ptr(wt, igrid)

            CALL p_f%get_ptr(p, igrid)
            CALL g_f%get_ptr(g, igrid)

            CALL dx_f%get_ptr(dx, igrid)
            CALL dy_f%get_ptr(dy, igrid)
            CALL dz_f%get_ptr(dz, igrid)

            CALL ddx_f%get_ptr(ddx, igrid)
            CALL ddy_f%get_ptr(ddy, igrid)
            CALL ddz_f%get_ptr(ddz, igrid)

            CALL rdx_f%get_ptr(rdx, igrid)
            CALL rdy_f%get_ptr(rdy, igrid)
            CALL rdz_f%get_ptr(rdz, igrid)

            CALL rddx_f%get_ptr(rddx, igrid)
            CALL rddy_f%get_ptr(rddy, igrid)
            CALL rddz_f%get_ptr(rddz, igrid)

            CALL tstle4_kon(kk, jj, ii, uo, vo, wo, u, v, w, ut, vt, wt, &
                dx, dy, dz, ddx, ddy, ddz, rdx, rdy, rdz, rddx, rddy, rddz, &
                nfro, nbac, nrgt, nlft, nbot, ntop)

            CALL tstle4_diff(kk, jj, ii, uo, vo, wo, u, v, w, g, &
                dx, dy, dz, ddx, ddy, ddz, rdx, rdy, rdz, rddx, rddy, rddz, &
                nfro, nbac, nrgt, nlft, nbot, ntop)

            CALL tstle4_gradp(kk, jj, ii, uo, vo, wo, p, dx, dy, dz, &
                nfro, nbac, nrgt, nlft, nbot, ntop, igrid)

            CALL tstle4_par(kk, jj, ii, uo, vo, wo, u, v, w, ut, vt, wt, &
                dx, dy, dz, ddx, ddy, ddz, rdx, rdy, rdz, rddx, rddy, rddz, &
                nfro, nbac, nrgt, nlft, nbot, ntop)
        END DO

        CALL stop_timer(310)
    END SUBROUTINE tstle4


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
    SUBROUTINE tstle4_kon(kk, jj, ii, uo, vo, wo, u, v, w, ut, vt, wt, &
            dx, dy, dz, ddx, ddy, ddz, rdx, rdy, rdz, rddx, rddy, rddz, &
            nfro, nbac, nrgt, nlft, nbot, ntop, bu, bv, bw)
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
        INTEGER, INTENT(in) :: nfro, nbac, nrgt, nlft, nbot, ntop
        REAL(realk), INTENT(in), OPTIONAL :: bu(kk, jj, ii), bv(kk, jj, ii), &
            bw(kk, jj, ii)

        ! Local variables
        INTEGER(intk) :: k, j, i
        INTEGER(intk) :: nbu, nfu, nrv, nbw, ntw, nlv
        REAL(realk) :: ax, ay, az
        REAL(realk) :: fw, fe, ft, fb, fn, fs
        REAL(realk) :: qw, qe, qt, qb, qn, qs

        ! Sanity check
        IF (PRESENT(bu) .NEQV. PRESENT(bv) .OR. &
                PRESENT(bu) .NEQV. PRESENT(bw)) THEN
            CALL errr(__FILE__, __LINE__)
        END IF

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

                    uo(k, j, i) = -(qe-qw+qn-qs+qt-qb)
                END DO

                IF (PRESENT(bu)) THEN
                    DO k = 3, kk-2
                        uo(k, j, i) = bu(k, j, i)*uo(k, j, i)
                    END DO
                ELSE
                    DO k = 3, kk-2
                        uo(k, j, i) = rdx(i)*rddy(j)*rddz(k)*uo(k, j, i)
                    END DO
                END IF
            END DO
        END DO

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

                    vo(k, j, i) = -(qe-qw+qn-qs+qt-qb)
                END DO

                IF (PRESENT(bv)) THEN
                    DO k = 3, kk-2
                        vo(k, j, i) = bv(k, j, i)*vo(k, j, i)
                    END DO
                ELSE
                    DO k = 3, kk-2
                        vo(k, j, i) = rddx(i)*rdy(j)*rddz(k)*vo(k, j, i)
                    END DO
                END IF
            END DO
        END DO

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

                    wo(k, j, i) = -(qe-qw+qn-qs+qt-qb)
                END DO

                IF (PRESENT(bw)) THEN
                    DO k = 3-nbw, kk-3+ntw
                        wo(k, j, i) = bw(k, j, i)*wo(k, j, i)
                    END DO
                ELSE
                    DO k = 3-nbw, kk-3+ntw
                        wo(k, j, i) = rddx(i)*rddy(j)*rdz(k)*wo(k, j, i)
                    END DO
                END IF
            END DO
        END DO
    END SUBROUTINE tstle4_kon


    SUBROUTINE tstle4_diff(kk, jj, ii, uo, vo, wo, u, v, w, g, &
            dx, dy, dz, ddx, ddy, ddz, rdx, rdy, rdz, rddx, rddy, rddz, &
            nfro, nbac, nrgt, nlft, nbot, ntop)
        ! Subroutine arguments
        INTEGER(intk), INTENT(in) :: kk, jj, ii
        REAL(realk), INTENT(inout) :: uo(kk, jj, ii), vo(kk, jj, ii), &
            wo(kk, jj, ii)
        REAL(realk), INTENT(in) :: u(kk, jj, ii), v(kk, jj, ii), w(kk, jj, ii)
        REAL(realk), INTENT(in) :: g(kk, jj, ii)
        REAL(realk), INTENT(in) :: dx(ii), dy(jj), dz(kk)
        REAL(realk), INTENT(in) :: ddx(ii), ddy(jj), ddz(kk)
        REAL(realk), INTENT(in) :: rdx(ii), rdy(jj), rdz(kk)
        REAL(realk), INTENT(in) :: rddx(ii), rddy(jj), rddz(kk)
        INTEGER, INTENT(in) :: nfro, nbac, nrgt, nlft, nbot, ntop

        ! Local variables
        INTEGER(intk) :: k, j, i
        INTEGER(intk) :: nbu, nfu, nrv, nbw, ntw, nlv
        INTEGER(intk) :: iles
        REAL(realk) :: ax, ay, az
        REAL(realk) :: ge, gw, gn, gs, gt, gb
        REAL(realk) :: qw, qe, qt, qb, qn, qs
        REAL(realk) :: st, qc, fak

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

        iles = 1
        IF (ilesmodel == 0) iles = 0

        CALL swcle3d(kk, jj, ii, uo, vo, wo, u, v, w, dx, dy, dz, &
            ddx, ddy, ddz, nfro, nbac, nrgt, nlft, nbot, ntop)

        DO i = 3-nfu, ii-3+nbu
            DO j = 3, jj-2
                DO k = 3, kk-2
                    ax = ddy(j)*ddz(k)
                    ay = dx(i)*ddz(k)
                    az = dx(i)*ddy(j)

                    ge = g(k, j, i+1)
                    gw = g(k, j, i)
                    gn = g(k, j, i)*g(k, j+1, i) &
                        /MAX(g(k, j, i) + g(k, j+1, i), gmol) &
                        + g(k, j, i+1)*g(k, j+1, i+1) &
                        /MAX(g(k, j, i+1) + g(k, j+1, i+1), gmol)
                    gs = g(k, j-1, i)*g(k, j, i) &
                        /MAX(g(k, j-1, i) + g(k, j, i), gmol) &
                        + g(k, j-1, i+1)*g(k, j, i+1) &
                        /MAX(g(k, j-1, i+1) + g(k, j, i+1), gmol)
                    gt = g(k, j, i)*g(k+1, j, i) &
                        /MAX(g(k, j, i) + g(k+1, j, i), gmol) &
                        + g(k, j, i+1)*g(k+1, j, i+1) &
                        /MAX(g(k, j, i+1) + g(k+1, j, i+1), gmol)
                    gb = g(k-1, j, i)*g(k, j, i) &
                        /MAX(g(k-1, j, i) + g(k, j, i), gmol) &
                        + g(k-1, j, i+1)*g(k, j, i+1) &
                        /MAX(g(k-1, j, i+1) + g(k, j, i+1), gmol)

                    qe = -ge*ax*rddx(i+1)*(u(k, j, i+1) - u(k, j, i))
                    qw = -gw*ax*rddx(i)*(u(k, j, i) - u(k, j, i-1))
                    qn = -gn*ay*rdy(j)*(u(k, j+1, i) - u(k, j, i))
                    qs = -gs*ay*rdy(j-1)*(u(k, j, i) - u(k, j-1, i))
                    qt = -gt*az*rdz(k)*(u(k+1, j, i) - u(k, j, i))
                    qb = -gb*az*rdz(k-1)*(u(k, j, i) - u(k-1, j, i))

                    st = ((ge*(u(k, j, i+1) - u(k, j, i))*rddx(i+1)) &
                        - (gw*(u(k, j, i) - u(k, j, i-1))*rddx(i)))*ax &
                        + ((gn*(v(k, j, i+1) - v(k, j, i))) &
                        - (gs*(v(k, j-1, i+1) - v(k, j-1, i))))*ddz(k) &
                        + ((gt*(w(k, j, i+1) - w(k, j, i))) &
                        - (gb*(w(k-1, j, i+1) - w(k-1, j, i))))*ddy(j)
                    qc = st*iles

                    fak = 1.0/rho*rddy(j)*rdx(i)*rddz(k)
                    uo(k, j, i) = uo(k, j, i) - fak*(qe-qw+qn-qs+qt-qb-qc)
                END DO
            END DO
        END DO

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
    END SUBROUTINE tstle4_diff


    SUBROUTINE tstle4_gradp(kk, jj, ii, uo, vo, wo, p, dx, dy, dz, &
            nfro, nbac, nrgt, nlft, nbot, ntop, igrid)
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

        DO i = 3-nfu, ii-3+nbu
            DO j = 3, jj-2
                DO k = 3, kk-2
                    uo(k, j, i) = uo(k, j, i) - 1.0/(rho*dx(i)) &
                        *(p(k, j, i+1) - p(k, j, i) + gpx*dx(i))
                END DO
            END DO
        END DO

        DO i = 3, ii-2
            DO j = 3-nrv, jj-3+nlv
                DO k = 3, kk-2
                    vo(k, j, i) = vo(k, j, i) - 1.0/(rho*dy(j)) &
                        *(p(k, j+1, i) - p(k, j, i) + gpy*dy(j))
                END DO
            END DO
        END DO

        DO i = 3, ii-2
            DO j = 3, jj-2
                DO k = 3-nbw, kk-3+ntw
                    wo(k, j, i) = wo(k, j, i) - 1.0/(rho*dz(k)) &
                        *(p(k+1, j, i) - p(k, j, i) + gpz*dz(k))
                END DO
            END DO
        END DO
    END SUBROUTINE tstle4_gradp


    SUBROUTINE tstle4_par(kk, jj, ii, uo, vo, wo, u, v, w, ut, vt, wt, &
            dx, dy, dz, ddx, ddy, ddz, rdx, rdy, rdz, rddx, rddy, rddz, &
            nfro, nbac, nrgt, nlft, nbot, ntop)
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

        ! Temporary storage
        ! TODO: Can this be replaced with a single tmp(:, :, :) array??
        ! I see no places where three arrays are needed at the same time...
        REAL(realk), ALLOCATABLE :: wcu(:, :, :), wcv(:, :, :), wcw(:, :, :)

        ALLOCATE(wcu(kk, jj, ii))
        ALLOCATE(wcv(kk, jj, ii))
        ALLOCATE(wcw(kk, jj, ii))

        ! Upwind in vorletzter schicht bei PAR-randbedingung
        IF (nfro == 8) THEN
            i = 4
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
        END IF

        IF (nbac == 8) THEN
            i = ii-2
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
        END IF

        IF (nrgt == 8) THEN
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
        END IF

        IF (nlft == 8) THEN
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
        END IF

        IF (nbot == 8) THEN
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
        END IF

        IF (ntop == 8) THEN
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

        END IF

        ! PAR-RB Impulserhaltend BACK
        IF (nbac == 8) THEN
            ! W-Impulszelle
            ! STANDART-BERECHNUNG DES KONVEKTIVEN FLUSSES
            i = ii-2
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

            ! NEUE BERECHNUNG DES KONVEKTIVEN FLUSSES
            ! FUER JEDE GROBGITTERGESCHWINDIGKEITT
            i = ii-2
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

            ! VERTEILUNG
            i = ii-2
            DO j = 3, jj-2
                DO k = 2, kk-4, 2
                    wcw(k+1, j, i) = 0.5*(wcw(k, j, i) + wcw(k+2, j, i))
                END DO
            END DO

            ! AUF WO SCHREIBEN
            DO j = 3, jj-2
                fkdtw = -1.0*rddx(i)*rddy(j)*wkon
                DO k = 2, kk-2
                    rdzk = rdz(k)
                    wo(k, j, i) = wo(k, j, i) &
                        + fkdtw*rdzk*(-wcw(k, j, i+1) + wcw(k, j, i))
                END DO
            END DO

            ! V-Impulszelle
            ! STANDARD-BERECHNUNG DES KONVEKTIVEN FLUSSES
            i = ii-2
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

            ! NEUE BERECHNUNG DES KONVEKTIVEN FLUSSES
            ! FUER JEDE GROBGITTERGESCHWINDIGKEIT
            i = ii-2

            ! YM-RAND
            j = 2
            dyj = dy(j)
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

            ! IM GEBIET
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

            ! YP-RAND
            j = jj-2
            dyj = dy(j)
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

            ! VERTEILUNG
            i = ii-2
            DO j = 2, jj-4, 2
                DO k = 3, kk-2
                    wcv(k, j+1, i) = 0.5*(wcv(k, j, i) + wcv(k, j+2, i))
                END DO
            END DO

            ! AUF VO SCHREIBEN
            DO j = 2, jj-2
                fkdtv = -1.0*rddx(i)* rdy(j)*wkon
                DO k = 3, kk-2
                    rddzk = rddz(k)
                    vo(k, j, i) = vo(k, j, i) &
                        + fkdtv*rddzk*(wcv(k, j, i) - wcv(k, j, i+1))
                END DO
            END DO
        END IF

        ! PAR-RB Impulserhaltend FRONT
        IF (nfro == 8) THEN
            ! W-Impulszelle
            ! STANDART-BERECHNUNG DES KONVEKTIVEN FLUSSES
            i = 3
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

            ! NEUE BERECHNUNG DES KONVEKTIVEN FLUSSES
            ! FUER JEDE GROBGITTERGESCHWINDIGKEIT
            i = 3
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

            ! VERTEILUNG
            i = 3
            DO j = 3, jj-2
                DO k = 2, kk-4, 2
                    wcw(k+1, j, i) = 0.5*(wcw(k, j, i) + wcw(k+2, j, i))
                END DO
            END DO

            ! AUF WO SCHREIBEN
            DO j = 3, jj-2
                fkdtw = -1.0*rddx(i)*rddy(j)*wkon
                DO k = 2, kk-2
                    rdzk = rdz(k)
                    wo(k, j, i) = wo(k, j, i) &
                        + fkdtw*rdzk*(wcw(k, j, i-1) - wcw(k, j, i))
                END DO
            END DO

            ! V-Impulszelle
            ! STANDARD-BERECHNUNG DES KONVEKTIVEN FLUSSES
            i = 3
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

            ! NEUE BERECHNUNG DES KONVEKTIVEN FLUSSES
            ! FUER JEDE GROBGITTERGESCHWINDIGKEIT
            i = 3

            ! YM-RAND
            j = 2
            dyj = dy(j)
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

            ! IM GEBIET
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

            ! YP-RAND
            j = jj-2
            dyj = dy(j)
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

            ! VERTEILUNG
            i = 3
            DO j = 2, jj-4, 2
                DO k = 3, kk-2
                    wcv(k, j+1, i) = 0.5*(wcv(k, j, i) + wcv(k, j+2, i))
                END DO
            END DO

            ! AUF VO SCHREIBEN
            DO j = 2, jj-2
                fkdtv = -1.0*rddx(i)*rdy(j)*wkon
                DO k = 3, kk-2
                    rddzk = rddz(k)
                    vo(k, j, i) = vo(k, j, i) &
                        + fkdtv*rddzk*(-wcv(k, j, i) + wcv(k, j, i-1))
                END DO
            END DO
        END IF

        ! PAR-RB Impulserhaltend TOP
        IF (ntop == 8) THEN

            ! U-Impulszelle
            ! STANDARD-BERECHNUNG DES KONVEKTIVEN FLUSSES
            k = kk-2
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

            ! NEUE BERECHNUNG DES KONVEKTIVEN FLUSSES
            ! FUER JEDE GROBGITTERGESCHWINDIGKEIT
            k = kk-2

            ! XM-RAND
            i = 2
            dxi = dx(i)
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

            ! IM GEBIET
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

            ! XP-RAND
            i = ii-2
            dxi = dx(i)
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

            ! VERTEILUNG
            k = kk-2
            DO i = 2, ii-4, 2
                DO j = 3, jj-2
                    wcu(k, j, i+1) = 0.5*(wcu(k, j, i) + wcu(k, j, i+2))
                END DO
            END DO

            ! AUF U0 SCHREIBEN
            rddzk = rddz(k)
            DO i = 2, ii-2
                DO j = 3, jj-2
                    fkdtu = -1.0*rddy(j)*rdx(i)*wkon
                    uo(k, j, i) = uo(k, j, i) &
                        + fkdtu*rddzk*(wcu(k, j, i) - wcu(k+1, j, i))
                END DO
            END DO

            ! V-Impulszelle
            ! STANDARD-BERECHNUNG DES KONVEKTIVEN FLUSSES
            k = kk-2
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

            ! NEUE BERECHNUNG DES KONVEKTIVEN FLUSSES
            ! FUER JEDE GROBGITTERGESCHWINDIGKEIT
            k = kk-2
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

            ! VERTEILUNG
            k = kk-2
            DO i = 3, ii-2
                DO j = 2, jj-4, 2
                    wcv(k, j+1, i) = 0.5*(wcv(k, j, i) + wcv(k, j+2, i))
                END DO
            END DO

            ! AUF VO SCHREIBEN
            rddzk = rddz(k)
            DO i = 3, ii-2
                DO j = 2, jj-2
                    fkdtv = -1.0*rddx(i)*rdy(j)*wkon
                    vo(k, j, i) = vo(k, j, i) &
                        + fkdtv*rddzk*(wcv(k, j, i) - wcv(k+1, j, i))
                END DO
            END DO
        END IF

        ! PAR-RB Impulserhaltend BOTTOM
        IF (nbot == 8) THEN
            ! U-Impulszelle
            ! STANDARD-BERECHNUNG DES KONVEKTIVEN FLUSSES
            k = 3
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

            ! NEUE BERECHNUNG DES KONVEKTIVEN FLUSSES
            ! FUER JEDE GROBGITTERGESCHWINDIGKEIT
            k = 3

            ! XM-RAND
            i = 2
            dxi = dx(i)
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

            ! IM GEBIET
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

            ! XP-RAND
            i = ii-2
            dxi = dx(i)
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

            ! VERTEILUNG
            k = 3
            DO i = 2, ii-4, 2
                DO j = 3, jj-2
                    wcu(k, j, i+1) = 0.5*(wcu(k, j, i) + wcu(k, j, i+2))
                END DO
            END DO

            ! AUF U0 SCHREIBEN
            rddzk = rddz(k)
            DO i = 2, ii-2
                DO j = 3, jj-2
                    fkdtu = -1.0*rddy(j)*rdx(i)*wkon
                    uo(k, j, i) = uo(k, j, i) &
                        + fkdtu*rddzk*(-wcu(k, j, i) + wcu(k-1, j, i))
                END DO
            END DO

            ! V-Impulszelle
            ! STANDARD-BERECHNUNG DES KONVEKTIVEN FLUSSES
            k = 3
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

            ! NEUE BERECHNUNG DES KONVEKTIVEN FLUSSES
            ! FUER JEDE GROBGITTERGESCHWINDIGKEIT
            k = 3
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

            ! VERTEILUNG
            k = 3
            DO i = 3, ii-2
                DO j = 2, jj-4, 2
                    wcv(k, j+1, i) = 0.5*(wcv(k, j, i) + wcv(k, j+2, i))
                END DO
            END DO

            ! AUF VO SCHREIBEN
            rddzk = rddz(k)
            DO i = 3, ii-2
                DO j = 2, jj-2
                    fkdtv = -1.0*rddx(i)*rdy(j)*wkon
                    vo(k, j, i) = vo(k, j, i) &
                        + fkdtv*rddzk*(-wcv(k, j, i) + wcv(k-1, j, i))
                END DO
            END DO
        END IF

        ! PAR-RB Impulserhaltend LEFT
        IF (nlft == 8) THEN
            ! U-Impulszelle
            ! STANDARD-BERECHNUNG DES KONVEKTIVEN FLUSSES
            j = jj-2
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

            ! NEUE BERECHNUNG DES KONVEKTIVEN FLUSSES
            ! FUER JEDE GROBGITTERGESCHWINDIGKEIT
            j = jj-2

            ! XM-RAND
            i = 2
            dxi = dx(i)
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

            ! IM GEBIET
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

            ! XP-RAND
            i = ii-2
            dxi = dx(i)
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

            ! VERTEILUNG
            j = jj-2
            DO i = 2, ii-4, 2
                DO k = 3, kk-2
                    wcu(k, j, i+1) = 0.5*(wcu(k, j, i) + wcu(k, j, i+2))
                END DO
            END DO

            ! AUF U0 SCHREIBEN
            j = jj-2
            DO i = 2, ii-2
                fkdtu = -1.0*rddy(j)*rdx(i)*wkon
                DO k = 3, kk-2
                    rddzk = rddz(k)
                    uo(k, j, i) = uo(k, j, i) &
                        + fkdtu*rddzk*(wcu(k, j, i) - wcu(k, j+1, i))
                END DO
            END DO

            ! W-IMPULSZELLE
            ! STANDARD-BERECHNUNG DES KONVEKTIVEN FLUSSES
            j = jj-2
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

            ! NEUE BERECHNUNG DES KONVEKTIVEN FLUSSES
            ! FUER JEDE GROBGITTERGESCHWINDIGKEIT
            j = jj-2
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

            ! VERTEILUNG
            j = jj-2
            DO i = 3, ii-2
                DO k = 2, kk-4, 2
                    wcw(k+1, j, i) = 0.5*(wcw(k, j, i) + wcw(k+2, j, i))
                END DO
            END DO

            ! AUF W0 SCHREIBEN
            j = jj-2
            DO i = 3, ii-2
                fkdtw = -1.0*rddx(i)*rddy(j)*wkon
                DO k = 2, kk-2
                    rdzk = rdz(k)
                    wo(k, j, i) = wo(k, j, i) &
                        + fkdtw*rdzk*(-wcw(k, j+1, i) + wcw(k, j, i))
                END DO
            END DO
        END IF

        ! PAR-RB Impulserhaltend RIGHT
        IF (nrgt == 8) THEN
            ! U-Impulszelle
            ! STANDARD-BERECHNUNG DES KONVEKTIVEN FLUSSES
            j = 3
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

            ! NEUE BERECHNUNG DES KONVEKTIVEN FLUSSES
            ! FUER JEDE GROBGITTERGESCHWINDIGKEIT
            j = 3

            ! XM-RAND
            i = 2
            dxi = dx(i)
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

            ! IM GEBIET
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

            ! XP-RAND
            i = ii-2
            dxi = dx(i)
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

            ! VERTEILUNG
            j = 3
            DO i = 2, ii-4, 2
                DO k = 3, kk-2
                    wcu(k, j, i+1) = 0.5*(wcu(k, j, i) + wcu(k, j, i+2))
                END DO
            END DO

            ! AUF U0 SCHREIBEN
            j = 3
            DO i = 2, ii-2
                fkdtu = -1.0*rddy(j)*rdx(i)*wkon
                DO k = 3, kk-2
                    rddzk = rddz(k)
                    uo(k, j, i) = uo(k, j, i) &
                        + fkdtu*rddzk*(-wcu(k, j, i) + wcu(k, j-1, i))
                END DO
            END DO

            ! W-IMPULSZELLE
            ! STANDARD-BERECHNUNG DES KONVEKTIVEN FLUSSES
            j = 3
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

            ! NEUE BERECHNUNG DES KONVEKTIVEN FLUSSES
            ! FUER JEDE GROBGITTERGESCHWINDIGKEIT
            j = 3
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

            ! VERTEILUNG
            j = 3
            DO i = 3, ii-2
                DO k = 2, kk-4, 2
                    wcw(k+1, j, i) = 0.5*(wcw(k, j, i) + wcw(k+2, j, i))
                END DO
            END DO

            ! AUF W0 SCHREIBEN
            j = 3
            DO i = 3, ii-2
                fkdtw = -1.0*rddx(i)*rddy(j)*wkon
                DO k = 2, kk-2
                    rdzk = rdz(k)
                    wo(k, j, i) = wo(k, j, i) &
                        + fkdtw*rdzk*(wcw(k, j-1, i) - wcw(k, j, i))
                END DO
            END DO
        END IF

        DEALLOCATE(wcu)
        DEALLOCATE(wcv)
        DEALLOCATE(wcw)
    END SUBROUTINE tstle4_par


    SUBROUTINE swcle3d(kk, jj, ii, uo, vo, wo, u, v, w, dx, dy, dz, &
            ddx, ddy, ddz, nfro, nbac, nrgt, nlft, nbot, ntop)

        ! Subroutine arguments
        INTEGER(intk), INTENT(in) :: kk, jj, ii
        REAL(realk), INTENT(inout) :: uo(kk, jj, ii), vo(kk, jj, ii), &
            wo(kk, jj, ii)
        REAL(realk), INTENT(in) :: u(kk, jj, ii), v(kk, jj, ii), w(kk, jj, ii)
        REAL(realk), INTENT(in) :: dx(ii), dy(jj), dz(kk), &
            ddx(ii), ddy(jj), ddz(kk)
        INTEGER(intk), INTENT(in) :: nfro, nbac, nrgt, nlft, nbot, ntop

        ! Local variables
        INTEGER(intk) :: k, j, i

        IF (nfro == 5) THEN
            i = 3
            DO j = 2, jj-1
                DO k = 2, kk-1
                    vo(k, j, i) = vo(k, j, i) &
                        - swcle3d_one(dy(j), ddz(k), ddx(i), v(k, j, i))
                    wo(k, j, i) = wo(k, j, i) &
                        - swcle3d_one(dz(k), ddy(j), ddx(i), w(k, j, i))
                END DO
            END DO
        END IF

        IF (nbac == 5) THEN
            i = ii-2
            DO j = 2, jj-1
                DO k = 2, kk-1
                    vo(k, j, i) = vo(k, j, i) &
                        - swcle3d_one(dy(j), ddz(k), ddx(i), v(k, j, i))
                    wo(k, j, i) = wo(k, j, i) &
                        - swcle3d_one(dz(k), ddy(j), ddx(i), w(k, j, i))
                END DO
            END DO
        END IF

        IF (nrgt == 5) THEN
            j = 3
            DO i = 2, ii-1
                DO k = 2, kk-1
                    uo(k, j, i) = uo(k, j, i) &
                        - swcle3d_one(dx(i), ddz(k), ddy(j), u(k, j, i))
                    wo(k, j, i) = wo(k, j, i) &
                        - swcle3d_one(dz(k), ddx(i), ddy(j), w(k, j, i))
                END DO
            END DO
        END IF

        IF (nlft == 5) THEN
            j = jj-2
            DO i = 2, ii-1
                DO k = 2, kk-1
                    uo(k, j, i) = uo(k, j, i) &
                        - swcle3d_one(dx(i), ddz(k), ddy(j), u(k, j, i))
                    wo(k, j, i) = wo(k, j, i) &
                        - swcle3d_one(dz(k), ddx(i), ddy(j), w(k, j, i))
                END DO
            END DO
        END IF

        IF (nbot == 5) THEN
            k = 3
            DO i = 2, ii-1
                DO j = 2, jj-1
                    uo(k, j, i) = uo(k, j, i) &
                        - swcle3d_one(dx(i), ddy(j), ddz(k), u(k, j, i))
                    vo(k, j, i) = vo(k, j, i) &
                        - swcle3d_one(dy(j), ddx(i), ddz(k), v(k, j, i))
                END DO
            END DO
        END IF

        IF (ntop == 5) THEN
            k = kk-2
            DO i = 2, ii-1
                DO j = 2, jj-1
                    uo(k, j, i) = uo(k, j, i) &
                        - swcle3d_one(dx(i), ddy(j), ddz(k), u(k, j, i))
                    vo(k, j, i) = vo(k, j, i) &
                        - swcle3d_one(dy(j), ddx(i), ddz(k), v(k, j, i))
                END DO
            END DO
        END IF
    END SUBROUTINE swcle3d


    PURE ELEMENTAL REAL(realk) FUNCTION swcle3d_one(dx, ddy, ddz, u) RESULT(uo)
        !$omp declare simd(swcle3d_one)

        ! Function arguments
        REAL(realk), INTENT(in) :: dx   ! in direction of u
        REAL(realk), INTENT(in) :: ddy  ! spanwise
        REAL(realk), INTENT(in) :: ddz  ! wall normal
        REAL(realk), INTENT(in) :: u    ! velocity

        ! Local variables
        REAL(realk) :: WCUZ

        wcuz = dx*ddy*tauwin(u, ddz)/rho
        uo = wcuz/(dx*ddy*ddz)
    END FUNCTION swcle3d_one
END MODULE tstle4_mod
