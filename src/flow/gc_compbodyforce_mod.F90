MODULE gc_compbodyforce_mod
    USE, INTRINSIC :: ISO_C_BINDING, ONLY: C_F_POINTER, C_PTR
    USE MPI_f08

    USE core_mod
    USE ib_mod
    USE flowcore_mod

    IMPLICIT NONE(type, external)
    PRIVATE

    TYPE, BIND(C) :: bodyforce_t
        REAL(c_realk) :: force(3, 3)
        REAL(c_realk) :: power(3, 3)
        REAL(c_realk) :: torque(3, 3)
    END TYPE bodyforce_t

    TYPE(MPI_Datatype) :: mpitype
    TYPE(MPI_Op) :: mpiop

    CHARACTER(len=*), PARAMETER :: bodyforcefile = "bodyforce.log"
    CHARACTER(len=*), PARAMETER :: bodypowerfile = "bodypower.log"
    CHARACTER(len=*), PARAMETER :: bodytorquefile = "bodytorque.log"

    PUBLIC :: init_compbodyforce, finish_compbodyforce, sample_compbodyforce

CONTAINS

    SUBROUTINE init_compbodyforce(dcont)
        ! Subroutine arguments
        LOGICAL, INTENT(in) :: dcont

        ! Local variables
        ! none...

        CALL create_mpi_dtype()
        CALL MPI_Op_create(bodyforce_reduce_mpi, .TRUE., mpiop)
        CALL bodyforce_init_logs(dcont)
    END SUBROUTINE init_compbodyforce


    SUBROUTINE finish_compbodyforce()
        CALL MPI_Op_free(mpiop)
        CALL MPI_Type_free(mpitype)
    END SUBROUTINE finish_compbodyforce


    SUBROUTINE sample_compbodyforce(u_f, v_f, w_f, p_f, g_f, ittot, time)
        ! Subroutine arguments
        TYPE(field_t), INTENT(in) :: u_f
        TYPE(field_t), INTENT(in) :: v_f
        TYPE(field_t), INTENT(in) :: w_f
        TYPE(field_t), INTENT(in) :: p_f
        TYPE(field_t), INTENT(in) :: g_f
        INTEGER(intk), INTENT(in) :: ittot
        REAL(realk), INTENT(in) :: time

        ! Local variables
        TYPE(field_t), POINTER :: x_f, y_f, z_f
        TYPE(field_t), POINTER :: bu_f, bv_f, bw_f
        TYPE(field_t), POINTER :: dx_f, dy_f, dz_f, ddx_f, ddy_f, ddz_f
        TYPE(field_t), POINTER :: rdx_f, rdy_f, rdz_f, rddx_f, rddy_f, rddz_f
        REAL(realk), POINTER, CONTIGUOUS, DIMENSION(:, :, :) :: u, v, w
        REAL(realk), POINTER, CONTIGUOUS, DIMENSION(:, :, :) :: bu, bv, bw
        REAL(realk), POINTER, CONTIGUOUS, DIMENSION(:, :, :) :: p, g
        REAL(realk), POINTER, CONTIGUOUS :: x(:), y(:), z(:)
        REAL(realk), POINTER, CONTIGUOUS :: dx(:), dy(:), dz(:)
        REAL(realk), POINTER, CONTIGUOUS :: ddx(:), ddy(:), ddz(:)
        REAL(realk), POINTER, CONTIGUOUS :: rdx(:), rdy(:), rdz(:)
        REAL(realk), POINTER, CONTIGUOUS :: rddx(:), rddy(:), rddz(:)
        INTEGER(intk) :: i, igrid, ilevel, nlevels
        INTEGER(intk) :: kk, jj, ii
        INTEGER(intk) :: nfro, nbac, nrgt, nlft, nbot, ntop

        TYPE(bodyforce_t) :: force(minlevel:maxlevel)

        CALL start_timer(351)
        CALL zero_bodyforce(force)

        CALL get_field(bu_f, "BU")
        CALL get_field(bv_f, "BV")
        CALL get_field(bw_f, "BW")

        CALL get_field(x_f, "X")
        CALL get_field(y_f, "Y")
        CALL get_field(z_f, "Z")

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
            ilevel = level(igrid)

            CALL get_mgdims(kk, jj, ii, igrid)
            CALL get_mgbasb(nfro, nbac, nrgt, nlft, nbot, ntop, igrid)

            CALL u_f%get_ptr(u, igrid)
            CALL v_f%get_ptr(v, igrid)
            CALL w_f%get_ptr(w, igrid)

            CALL p_f%get_ptr(p, igrid)
            CALL g_f%get_ptr(g, igrid)

            CALL bu_f%get_ptr(bu, igrid)
            CALL bv_f%get_ptr(bv, igrid)
            CALL bw_f%get_ptr(bw, igrid)

            CALL x_f%get_ptr(x, igrid)
            CALL y_f%get_ptr(y, igrid)
            CALL z_f%get_ptr(z, igrid)

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

            CALL compbodyforce_kon(force(ilevel)%force(:, 1), &
                force(ilevel)%power(:, 1), force(ilevel)%torque(:, 1), &
                kk, jj, ii, u, v, w, bu, bv, bw, x, y, z, &
                dx, dy, dz, ddx, ddy, ddz, rdx, rdy, rdz, rddx, rddy, rddz, &
                nfro, nbac, nrgt, nlft, nbot, ntop)

            CALL compbodyforce_diff(force(ilevel)%force(:, 2), &
                force(ilevel)%power(:, 2), force(ilevel)%torque(:, 2), &
                kk, jj, ii, u, v, w, g, bu, bv, bw, x, y, z, &
                dx, dy, dz, ddx, ddy, ddz, rdx, rdy, rdz, rddx, rddy, rddz, &
                nfro, nbac, nrgt, nlft, nbot, ntop)

            CALL compbodyforce_gradp(force(ilevel)%force(:, 3), &
                force(ilevel)%power(:, 3), force(ilevel)%torque(:, 3), &
                kk, jj, ii, u, v, w, p, bu, bv, bw, x, y, z, ddx, ddy, ddz, &
                nfro, nbac, nrgt, nlft, nbot, ntop, igrid)
        END DO

        ! Reduce and print
        nlevels = maxlevel - minlevel + 1
        ! MPI_Reduce fail with MPI_IN_PLACE... Allreduce works.
        CALL MPI_Allreduce(MPI_IN_PLACE, force, nlevels, &
            mpitype, mpiop, MPI_COMM_WORLD)

        CALL bodyforce_print_logs(force, ittot, time)

        CALL stop_timer(351)
    END SUBROUTINE sample_compbodyforce


    SUBROUTINE compbodyforce_kon(force, power, torque, kk, jj, ii, u, v, w, &
            bu, bv, bw, x, y, z, dx, dy, dz, ddx, ddy, ddz, rdx, rdy, rdz, &
            rddx, rddy, rddz, nfro, nbac, nrgt, nlft, nbot, ntop)
        ! Subroutine arguments
        REAL(realk), INTENT(inout) :: force(:), power(:), torque(:)
        INTEGER(intk), INTENT(in) :: kk, jj, ii
        REAL(realk), INTENT(in) :: u(kk, jj, ii), v(kk, jj, ii), w(kk, jj, ii)
        REAL(realk), INTENT(in) :: bu(kk, jj, ii), bv(kk, jj, ii), &
            bw(kk, jj, ii)
        REAL(realk), INTENT(in) :: x(ii), y(jj), z(kk)
        REAL(realk), INTENT(in) :: dx(ii), dy(jj), dz(kk)
        REAL(realk), INTENT(in) :: ddx(ii), ddy(jj), ddz(kk)
        REAL(realk), INTENT(in) :: rdx(ii), rdy(jj), rdz(kk)
        REAL(realk), INTENT(in) :: rddx(ii), rddy(jj), rddz(kk)
        INTEGER, INTENT(in) :: nfro, nbac, nrgt, nlft, nbot, ntop

        ! Local variables
        INTEGER(intk) :: k, j, i
        INTEGER(intk) :: nbu, nfu, nrv, nbw, ntw, nlv
        REAL(realk) :: ax, ay, az
        REAL(realk) :: fw, fe, ft, fb, fn, fs
        REAL(realk) :: qw, qe, qt, qb, qn, qs

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

                    fe = ax*(u(k, j, i) + (u(k, j, i+1) - u(k, j, i)) &
                        * 0.5*dx(i)/ddx(i+1))
                    fw = ax*(u(k, j, i-1) + (u(k, j, i) - u(k, j, i-1)) &
                        * 0.5*dx(i-1)/ddx(i))
                    fn = ay*(v(k, j, i) + v(k, j, i+1))*0.5
                    fs = ay*(v(k, j-1, i) + v(k, j-1, i+1))*0.5
                    ft = az*(w(k, j, i) + w(k, j, i+1))*0.5
                    fb = az*(w(k-1, j, i) + w(k-1, j, i+1))*0.5

                    qe = 0.5*fe*(u(k, j, i) + u(k, j, i+1))
                    qw = 0.5*fw*(u(k, j, i-1) + u(k, j, i))
                    qn = 0.5*fn*(u(k, j, i) + u(k, j+1, i))
                    qs = 0.5*fs*(u(k, j-1, i) + u(k, j, i))
                    qt = 0.5*ft*(u(k, j, i) + u(k+1, j, i))
                    qb = 0.5*fb*(u(k-1, j, i) + u(k, j, i))

                    force(1) = force(1) &
                        + (1.0 - bu(k, j, i+1))*qe &
                        - (1.0 - bu(k, j, i-1))*qw &
                        + (1.0 - bu(k, j+1, i))*qn &
                        - (1.0 - bu(k, j-1, i))*qs &
                        + (1.0 - bu(k+1, j, i))*qt &
                        - (1.0 - bu(k-1, j, i))*qb

                    power(1) = power(1) &
                        + (1.0 - bu(k, j, i+1))*qe*0.5*(u(k, j, i) + u(k, j, i+1)) &
                        - (1.0 - bu(k, j, i-1))*qw*0.5*(u(k, j, i-1) + u(k, j, i)) &
                        + (1.0 - bu(k, j+1, i))*qn*0.5*(u(k, j, i) + u(k, j+1, i)) &
                        - (1.0 - bu(k, j-1, i))*qs*0.5*(u(k, j-1, i) + u(k, j, i)) &
                        + (1.0 - bu(k+1, j, i))*qt*0.5*(u(k, j, i) + u(k+1, j, i)) &
                        - (1.0 - bu(k-1, j, i))*qb*0.5*(u(k-1, j, i) + u(k, j, i))

                    torque(2) = torque(2) &
                        + (1.0 - bu(k, j, i+1))*qe*z(k) &
                        - (1.0 - bu(k, j, i-1))*qw*z(k) &
                        + (1.0 - bu(k, j+1, i))*qn*z(k) &
                        - (1.0 - bu(k, j-1, i))*qs*z(k) &
                        + (1.0 - bu(k+1, j, i))*qt*(z(k) + dz(k)/2.0) &
                        - (1.0 - bu(k-1, j, i))*qb*(z(k) - dz(k-1)/2.0)

                    torque(3) = torque(3) &
                        + (1.0 - bu(k, j, i+1))*qe*(-y(j)) &
                        - (1.0 - bu(k, j, i-1))*qw*(-y(j)) &
                        + (1.0 - bu(k, j+1, i))*qn*(-(y(j) + dy(j)/2.0)) &
                        - (1.0 - bu(k, j-1, i))*qs*(-(y(j) - dy(j-1)/2.0)) &
                        + (1.0 - bu(k+1, j, i))*qt*(-y(j)) &
                        - (1.0 - bu(k-1, j, i))*qb*(-y(j))
                END DO
            END DO
        END DO

        DO i = 3, ii-2
            DO j = 3-nrv, jj-3+nlv
                DO k = 3, kk-2
                    ax = dy(j)*ddz(k)
                    ay = ddx(i)*ddz(k)
                    az = ddx(i)*dy(j)

                    fe = ax*(u(k, j, i) + u(k, j+1, i))*0.5
                    fw = ax*(u(k, j, i-1) + u(k, j+1, i-1))*0.5
                    fn = ay*(v(k, j, i) + (v(k, j+1, i) - v(k, j, i)) &
                        * 0.5*dy(j)/ddy(j+1))
                    fs = ay*(v(k, j-1, i) + (v(k, j, i) -v(k, j-1, i)) &
                        * 0.5*dy(j-1)/ddy(j))
                    ft = az*(w(k, j, i) + w(k, j+1, i))*0.5
                    fb = az*(w(k-1, j, i) + w(k-1, j+1, i))*0.5

                    qe = 0.5*fe*(v(k, j, i) + v(k, j, i+1))
                    qw = 0.5*fw*(v(k, j, i-1) + v(k, j, i))
                    qn = 0.5*fn*(v(k, j, i) + v(k, j+1, i))
                    qs = 0.5*fs*(v(k, j-1, i) + v(k, j, i))
                    qt = 0.5*ft*(v(k, j, i) + v(k+1, j, i))
                    qb = 0.5*fb*(v(k-1, j, i) + v(k, j, i))

                    force(2) = force(2) &
                        + (1.0 - bv(k, j, i+1))*qe &
                        - (1.0 - bv(k, j, i-1))*qw &
                        + (1.0 - bv(k, j+1, i))*qn &
                        - (1.0 - bv(k, j-1, i))*qs &
                        + (1.0 - bv(k+1, j, i))*qt &
                        - (1.0 - bv(k-1, j, i))*qb

                    power(2) = power(2) &
                        + (1.0 - bv(k, j, i+1))*qe*0.5*(v(k, j, i) + v(k, j, i+1)) &
                        - (1.0 - bv(k, j, i-1))*qw*0.5*(v(k, j, i-1) + v(k, j, i)) &
                        + (1.0 - bv(k, j+1, i))*qn*0.5*(v(k, j, i) + v(k, j+1, i)) &
                        - (1.0 - bv(k, j-1, i))*qs*0.5*(v(k, j-1, i) + v(k, j, i)) &
                        + (1.0 - bv(k+1, j, i))*qt*0.5*(v(k, j, i) + v(k+1, j, i)) &
                        - (1.0 - bv(k-1, j, i))*qb*0.5*(v(k-1, j, i) + v(k, j, i))

                    torque(1) = torque(1) &
                        + (1.0 - bv(k, j, i+1))*qe*(-z(k)) &
                        - (1.0 - bv(k, j, i-1))*qw*(-z(k)) &
                        + (1.0 - bv(k, j+1, i))*qn*(-z(k)) &
                        - (1.0 - bv(k, j-1, i))*qs*(-z(k)) &
                        + (1.0 - bv(k+1, j, i))*qt*(-(z(k) + dz(k)/2.0)) &
                        - (1.0 - bv(k-1, j, i))*qb*(-(z(k) - dz(k-1)/2.0))

                    torque(3) = torque(3) &
                        + (1.0 - bv(k, j, i+1))*qe*((x(i) + dx(i)/2.0)) &
                        - (1.0 - bv(k, j, i-1))*qw*((x(i) - dx(i-1)/2.0)) &
                        + (1.0 - bv(k, j+1, i))*qn*(x(i)) &
                        - (1.0 - bv(k, j-1, i))*qs*(x(i)) &
                        + (1.0 - bv(k+1, j, i))*qt*(x(i)) &
                        - (1.0 - bv(k-1, j, i))*qb*(x(i))
                END DO
            END DO
        END DO

        DO i = 3, ii-2
            DO j = 3, jj-2
                DO k = 3-nbw, kk-3+ntw
                    ax = ddy(j)*dz(k)
                    ay = ddx(i)*dz(k)
                    az = ddx(i)*ddy(j)

                    fe = ax*(u(k, j, i) + u(k+1, j, i))*0.5
                    fw = ax*(u(k, j, i-1)+ u(k+1, j, i-1))*0.5
                    fn = ay*(v(k, j, i) + v(k+1, j, i))*0.5
                    fs = ay*(v(k, j-1, i)+ v(k+1, j-1, i))*0.5
                    ft = az*(w(k, j, i) + (w(k+1, j, i) - w(k, j, i)) &
                        * 0.5*dz(k)/ddz(k+1))
                    fb = az*(w(k-1, j, i) + (w(k, j, i) - w(k-1, j, i)) &
                        * 0.5*dz(k-1)/ddz(k))

                    qe = 0.5*fe*(w(k, j, i) + w(k, j, i+1))
                    qw = 0.5*fw*(w(k, j, i-1) + w(k, j, i))
                    qn = 0.5*fn*(w(k, j, i) + w(k, j+1, i))
                    qs = 0.5*fs*(w(k, j-1, i) + w(k, j, i))
                    qt = 0.5*ft*(w(k, j, i) + w(k+1, j, i))
                    qb = 0.5*fb*(w(k-1, j, i) + w(k, j, i))

                    force(3) = force(3) &
                        + (1.0 - bw(k, j, i+1))*qe &
                        - (1.0 - bw(k, j, i-1))*qw &
                        + (1.0 - bw(k, j+1, i))*qn &
                        - (1.0 - bw(k, j-1, i))*qs &
                        + (1.0 - bw(k+1, j, i))*qt &
                        - (1.0 - bw(k-1, j, i))*qb

                    power(3) = power(3) &
                        + (1.0 - bw(k, j, i+1))*qe*0.5*(w(k, j, i) + w(k, j, i+1)) &
                        - (1.0 - bw(k, j, i-1))*qw*0.5*(w(k, j, i-1) + w(k, j, i)) &
                        + (1.0 - bw(k, j+1, i))*qn*0.5*(w(k, j, i) + w(k, j+1, i)) &
                        - (1.0 - bw(k, j-1, i))*qs*0.5*(w(k, j-1, i) + w(k, j, i)) &
                        + (1.0 - bw(k+1, j, i))*qt*0.5*(w(k, j, i) + w(k+1, j, i)) &
                        - (1.0 - bw(k-1, j, i))*qb*0.5*(w(k-1, j, i) + w(k, j, i))

                    torque(1) = torque(1) &
                        + (1.0 - bw(k, j, i+1))*qe*(y(j)) &
                        - (1.0 - bw(k, j, i-1))*qw*(y(j)) &
                        + (1.0 - bw(k, j+1, i))*qn*((y(j) + dy(j)/2.0)) &
                        - (1.0 - bw(k, j-1, i))*qs*((y(j) - dy(j-1)/2.0)) &
                        + (1.0 - bw(k+1, j, i))*qt*(y(j)) &
                        - (1.0 - bw(k-1, j, i))*qb*(y(j))

                    torque(2) = torque(2) &
                        + (1.0 - bw(k, j, i+1))*qe*(-(x(i) + dx(i)/2.0)) &
                        - (1.0 - bw(k, j, i-1))*qw*(-(x(i) - dx(i-1)/2.0)) &
                        + (1.0 - bw(k, j+1, i))*qn*(-x(i)) &
                        - (1.0 - bw(k, j-1, i))*qs*(-x(i)) &
                        + (1.0 - bw(k+1, j, i))*qt*(-x(i)) &
                        - (1.0 - bw(k-1, j, i))*qb*(-x(i))
                END DO
            END DO
        END DO
    END SUBROUTINE compbodyforce_kon


    SUBROUTINE compbodyforce_diff(force, power, torque, kk, jj, ii, &
            u, v, w, g, bu, bv, bw, x, y, z, &
            dx, dy, dz, ddx, ddy, ddz, rdx, rdy, rdz, rddx, rddy, rddz, &
            nfro, nbac, nrgt, nlft, nbot, ntop)
        ! Subroutine arguments
        REAL(realk), INTENT(inout) :: force(:), power(:), torque(:)
        INTEGER(intk), INTENT(in) :: kk, jj, ii
        REAL(realk), INTENT(in) :: u(kk, jj, ii), v(kk, jj, ii), w(kk, jj, ii)
        REAL(realk), INTENT(in) :: g(kk, jj, ii)
        REAL(realk), INTENT(in) :: bu(kk, jj, ii), bv(kk, jj, ii), &
            bw(kk, jj, ii)
        REAL(realk), INTENT(in) :: x(ii), y(jj), z(kk)
        REAL(realk), INTENT(in) :: dx(ii), dy(jj), dz(kk)
        REAL(realk), INTENT(in) :: ddx(ii), ddy(jj), ddz(kk)
        REAL(realk), INTENT(in) :: rdx(ii), rdy(jj), rdz(kk)
        REAL(realk), INTENT(in) :: rddx(ii), rddy(jj), rddz(kk)
        INTEGER, INTENT(in) :: nfro, nbac, nrgt, nlft, nbot, ntop

        ! Local variables
        INTEGER(intk) :: k, j, i
        INTEGER(intk) :: nbu, nfu, nrv, nbw, ntw, nlv
        REAL(realk) :: ax, ay, az
        REAL(realk) :: ge, gw, gn, gs, gt, gb
        REAL(realk) :: qw, qe, qt, qb, qn, qs

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

                    qe = qe - ((ge*(u(k, j, i+1) - u(k, j, i))*rddx(i+1)))*ax
                    qw = qw - ((gw*(u(k, j, i) - u(k, j, i-1))*rddx(i)))*ax
                    qn = qn - ((gn*(v(k, j, i+1) - v(k, j, i))))*ddz(k)
                    qs = qs - ((gs*(v(k, j-1, i+1) - v(k, j-1, i))))*ddz(k)
                    qt = qt - ((gt*(w(k, j, i+1) - w(k, j, i))))*ddy(j)
                    qb = qb - ((gb*(w(k-1, j, i+1) - w(k-1, j, i))))*ddy(j)

                    force(1) = force(1) &
                        + (1.0 - bu(k, j, i+1))*qe &
                        - (1.0 - bu(k, j, i-1))*qw &
                        + (1.0 - bu(k, j+1, i))*qn &
                        - (1.0 - bu(k, j-1, i))*qs &
                        + (1.0 - bu(k+1, j, i))*qt &
                        - (1.0 - bu(k-1, j, i))*qb

                    power(1) = power(1) &
                        + (1.0 - bu(k, j, i+1))*qe*0.5*(u(k, j, i) + u(k, j, i+1)) &
                        - (1.0 - bu(k, j, i-1))*qw*0.5*(u(k, j, i-1) + u(k, j, i)) &
                        + (1.0 - bu(k, j+1, i))*qn*0.5*(u(k, j, i) + u(k, j+1, i)) &
                        - (1.0 - bu(k, j-1, i))*qs*0.5*(u(k, j-1, i) + u(k, j, i)) &
                        + (1.0 - bu(k+1, j, i))*qt*0.5*(u(k, j, i) + u(k+1, j, i)) &
                        - (1.0 - bu(k-1, j, i))*qb*0.5*(u(k-1, j, i) + u(k, j, i))

                    torque(2) = torque(2) &
                        + (1.0 - bu(k, j, i+1))*qe*z(k) &
                        - (1.0 - bu(k, j, i-1))*qw*z(k) &
                        + (1.0 - bu(k, j+1, i))*qn*z(k) &
                        - (1.0 - bu(k, j-1, i))*qs*z(k) &
                        + (1.0 - bu(k+1, j, i))*qt*(z(k) + dz(k)/2.0) &
                        - (1.0 - bu(k-1, j, i))*qb*(z(k) - dz(k-1)/2.0)

                    torque(3) = torque(3) &
                        + (1.0 - bu(k, j, i+1))*qe*(-y(j)) &
                        - (1.0 - bu(k, j, i-1))*qw*(-y(j)) &
                        + (1.0 - bu(k, j+1, i))*qn*(-(y(j) + dy(j)/2.0)) &
                        - (1.0 - bu(k, j-1, i))*qs*(-(y(j) - dy(j-1)/2.0)) &
                        + (1.0 - bu(k+1, j, i))*qt*(-y(j)) &
                        - (1.0 - bu(k-1, j, i))*qb*(-y(j))
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

                    qe = qe - ((ge*(u(k, j+1, i) - u(k, j, i))))*ddz(k)
                    qw = qw - ((gw*(u(k, j+1, i-1) - u(k, j, i-1))))*ddz(k)
                    qn = qn - ((gn*(v(k, j+1, i) - v(k, j, i))*rddy(j+1)))*ay
                    qs = qs - ((gs*(v(k, j, i) - v(k, j-1, i))*rddy(j)))*ay
                    qt = qt - ((gt*(w(k, j+1, i) - w(k, j, i))))*ddx(i)
                    qb = qb - ((gb*(w(k-1, j+1, i)- w(k-1, j, i))))*ddx(i)

                    force(2) = force(2) &
                        + (1.0 - bv(k, j, i+1))*qe &
                        - (1.0 - bv(k, j, i-1))*qw &
                        + (1.0 - bv(k, j+1, i))*qn &
                        - (1.0 - bv(k, j-1, i))*qs &
                        + (1.0 - bv(k+1, j, i))*qt &
                        - (1.0 - bv(k-1, j, i))*qb

                    power(2) = power(2) &
                        + (1.0 - bv(k, j, i+1))*qe*0.5*(v(k, j, i) + v(k, j, i+1)) &
                        - (1.0 - bv(k, j, i-1))*qw*0.5*(v(k, j, i-1) + v(k, j, i)) &
                        + (1.0 - bv(k, j+1, i))*qn*0.5*(v(k, j, i) + v(k, j+1, i)) &
                        - (1.0 - bv(k, j-1, i))*qs*0.5*(v(k, j-1, i) + v(k, j, i)) &
                        + (1.0 - bv(k+1, j, i))*qt*0.5*(v(k, j, i) + v(k+1, j, i)) &
                        - (1.0 - bv(k-1, j, i))*qb*0.5*(v(k-1, j, i) + v(k, j, i))

                    torque(1) = torque(1) &
                        + (1.0 - bv(k, j, i+1))*qe*(-z(k)) &
                        - (1.0 - bv(k, j, i-1))*qw*(-z(k)) &
                        + (1.0 - bv(k, j+1, i))*qn*(-z(k)) &
                        - (1.0 - bv(k, j-1, i))*qs*(-z(k)) &
                        + (1.0 - bv(k+1, j, i))*qt*(-(z(k) + dz(k)/2.0)) &
                        - (1.0 - bv(k-1, j, i))*qb*(-(z(k) - dz(k-1)/2.0))

                    torque(3) = torque(3) &
                        + (1.0 - bv(k, j, i+1))*qe*((x(i) + dx(i)/2.0)) &
                        - (1.0 - bv(k, j, i-1))*qw*((x(i) - dx(i-1)/2.0)) &
                        + (1.0 - bv(k, j+1, i))*qn*(x(i)) &
                        - (1.0 - bv(k, j-1, i))*qs*(x(i)) &
                        + (1.0 - bv(k+1, j, i))*qt*(x(i)) &
                        - (1.0 - bv(k-1, j, i))*qb*(x(i))
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

                    qe = qe - ((ge*(u(k+1, j, i) - u(k, j, i))))*ddy(j)
                    qw = qw - ((gw*(u(k+1, j, i-1) - u(k, j, i-1))))*ddy(j)
                    qn = qn - ((gn*(v(k+1, j, i) - v(k, j, i))))*ddx(i)
                    qs = qs - ((gs*(v(k+1, j-1, i) - v(k, j-1, i))))*ddx(i)
                    qt = qt - ((gt*(w(k+1, j, i) - w(k, j, i)))*rddz(k+1))*az
                    qb = qb - ((gb*(w(k, j, i) - w(k-1, j, i)))*rddz(k))*az

                    force(3) = force(3) &
                        + (1.0 - bw(k, j, i+1))*qe &
                        - (1.0 - bw(k, j, i-1))*qw &
                        + (1.0 - bw(k, j+1, i))*qn &
                        - (1.0 - bw(k, j-1, i))*qs &
                        + (1.0 - bw(k+1, j, i))*qt &
                        - (1.0 - bw(k-1, j, i))*qb

                    power(3) = power(3) &
                        + (1.0 - bw(k, j, i+1))*qe*0.5*(w(k, j, i) + w(k, j, i+1)) &
                        - (1.0 - bw(k, j, i-1))*qw*0.5*(w(k, j, i-1) + w(k, j, i)) &
                        + (1.0 - bw(k, j+1, i))*qn*0.5*(w(k, j, i) + w(k, j+1, i)) &
                        - (1.0 - bw(k, j-1, i))*qs*0.5*(w(k, j-1, i) + w(k, j, i)) &
                        + (1.0 - bw(k+1, j, i))*qt*0.5*(w(k, j, i) + w(k+1, j, i)) &
                        - (1.0 - bw(k-1, j, i))*qb*0.5*(w(k-1, j, i) + w(k, j, i))

                    torque(1) = torque(1) &
                        + (1.0 - bw(k, j, i+1))*qe*(y(j)) &
                        - (1.0 - bw(k, j, i-1))*qw*(y(j)) &
                        + (1.0 - bw(k, j+1, i))*qn*((y(j) + dy(j)/2.0)) &
                        - (1.0 - bw(k, j-1, i))*qs*((y(j) - dy(j-1)/2.0)) &
                        + (1.0 - bw(k+1, j, i))*qt*(y(j)) &
                        - (1.0 - bw(k-1, j, i))*qb*(y(j))

                    torque(2) = torque(2) &
                        + (1.0 - bw(k, j, i+1))*qe*(-(x(i) + dx(i)/2.0)) &
                        - (1.0 - bw(k, j, i-1))*qw*(-(x(i) - dx(i-1)/2.0)) &
                        + (1.0 - bw(k, j+1, i))*qn*(-x(i)) &
                        - (1.0 - bw(k, j-1, i))*qs*(-x(i)) &
                        + (1.0 - bw(k+1, j, i))*qt*(-x(i)) &
                        - (1.0 - bw(k-1, j, i))*qb*(-x(i))
                END DO
            END DO
        END DO
    END SUBROUTINE compbodyforce_diff


    SUBROUTINE compbodyforce_gradp(force, power, torque, kk, jj, ii, u, v, w, &
            p, bu, bv, bw, x, y, z, ddx, ddy, ddz, nfro, nbac, nrgt, nlft, &
            nbot, ntop, igrid)
        ! Subroutine arguments
        REAL(realk), INTENT(inout) :: force(:), power(:), torque(:)
        INTEGER(intk), INTENT(in) :: kk, jj, ii
        REAL(realk), INTENT(in) :: u(kk, jj, ii), v(kk, jj, ii), w(kk, jj, ii)
        REAL(realk), INTENT(in) :: p(kk, jj, ii)
        REAL(realk), INTENT(in) :: bu(kk, jj, ii), bv(kk, jj, ii), &
            bw(kk, jj, ii)
        REAL(realk), INTENT(in) :: x(ii), y(jj), z(kk)
        REAL(realk), INTENT(in) :: ddx(ii), ddy(jj), ddz(kk)
        INTEGER, INTENT(in) :: nfro, nbac, nrgt, nlft, nbot, ntop
        INTEGER, INTENT(in) :: igrid

        ! Local variables
        INTEGER(intk) :: k, j, i
        INTEGER(intk) :: nbu, nfu, nrv, nbw, ntw, nlv

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
                    force(1) = force(1) &
                        + (1.0 - bu(k, j, i+1))*p(k, j, i+1)*ddy(j)*ddz(k) &
                        - (1.0 - bu(k, j, i-1))*p(k, j, i)*ddy(j)*ddz(k)

                    power(1) = power(1) &
                        + (1.0 - bu(k, j, i+1))*p(k, j, i+1)*ddy(j)*ddz(k) &
                        *0.5*(u(k, j, i) + u(k, j, i+1)) &
                        - (1.0 - bu(k, j, i-1))*p(k, j, i)*ddy(j)*ddz(k) &
                        *0.5*(u(k, j, i-1) + u(k, j, i))

                    torque(2) = torque(2) &
                        + (1.0 - bu(k, j, i+1))*p(k, j, i+1)*ddy(j)*ddz(k)*z(k) &
                        - (1.0 - bu(k, j, i-1))*p(k, j, i)*ddy(j)*ddz(k)*z(k)

                    torque(3) = torque(3) &
                        + (1.0 - bu(k, j, i+1))*p(k, j, i+1)*ddy(j)*ddz(k)*(-y(j)) &
                        - (1.0 - bu(k, j, i-1))*p(k, j, i)*ddy(j)*ddz(k)*(-y(j))
                END DO
            END DO
        END DO

        DO i = 3, ii-2
            DO j = 3-nrv, jj-3+nlv
                DO k = 3, kk-2
                    force(2) = force(2) &
                        + (1.0 - bv(k, j+1, i))*p(k, j+1, i)*ddx(i)*ddz(k) &
                        - (1.0 - bv(k, j-1, i))*p(k, j, i)*ddx(i)*ddz(k)

                    power(2) = power(2) &
                        + (1.0 - bv(k, j+1, i))*p(k, j+1, i)*ddx(i)*ddz(k) &
                        *0.5*(v(k, j, i) + v(k, j+1, i)) &
                        - (1.0 - bv(k, j-1, i))*p(k, j, i)*ddx(i)*ddz(k) &
                        *0.5*(v(k, j-1, i) + v(k, j, i))

                    torque(1) = torque(1) &
                        + (1.0 - bv(k, j+1, i))*p(k, j+1, i)*ddx(i)*ddz(k)*(-z(k)) &
                        - (1.0 - bv(k, j-1, i))*p(k, j, i)*ddx(i)*ddz(k)*(-z(k))

                    torque(3) = torque(3) &
                        + (1.0 - bv(k, j+1, i))*p(k, j+1, i)*ddx(i)*ddz(k)*(x(i)) &
                        - (1.0 - bv(k, j-1, i))*p(k, j, i)*ddx(i)*ddz(k)*(x(i))
                END DO
            END DO
        END DO

        DO i = 3, ii-2
            DO j = 3, jj-2
                DO k = 3-nbw, kk-3+ntw
                    force(3) = force(3) &
                        + (1.0 - bw(k+1, j, i))*p(k+1, j, i)*ddx(i)*ddy(j) &
                        - (1.0 - bw(k-1, j, i))*p(k, j, i)*ddx(i)*ddy(j)

                    power(3) = power(3) &
                        + (1.0 - bw(k+1, j, i))*p(k+1, j, i)*ddx(i)*ddy(j) &
                        *0.5*(w(k, j, i) + w(k+1, j, i)) &
                        - (1.0 - bw(k-1, j, i))*p(k, j, i)*ddx(i)*ddy(j) &
                        *0.5*(w(k-1, j, i) + w(k, j, i))

                    torque(1) = torque(1) &
                        + (1.0 - bw(k+1, j, i))*p(k+1, j, i)*ddx(i)*ddy(j)*(y(j)) &
                        - (1.0 - bw(k-1, j, i))*p(k, j, i)*ddx(i)*ddy(j)*(y(j))

                    torque(2) = torque(2) &
                        + (1.0 - bw(k+1, j, i))*p(k+1, j, i)*ddx(i)*ddy(j)*(-x(i)) &
                        - (1.0 - bw(k-1, j, i))*p(k, j, i)*ddx(i)*ddy(j)*(-x(i))
                END DO
            END DO
        END DO
    END SUBROUTINE compbodyforce_gradp


    SUBROUTINE create_mpi_dtype()
        ! Local variables
        INTEGER(intk) :: i
        INTEGER(int32) :: blocklen(3)
        TYPE(MPI_Datatype) :: types(3)
        INTEGER(MPI_ADDRESS_KIND) :: base, disp(3)
        TYPE(bodyforce_t) :: foo

        blocklen = 9

        CALL MPI_Get_address(foo%force, disp(1))
        CALL MPI_Get_address(foo%power, disp(2))
        CALL MPI_Get_address(foo%torque, disp(3))

        base = disp(1)
        DO i = 1, 3
            disp(i) = disp(i) - base
        END DO

        types(1) = mglet_mpi_real
        types(2) = mglet_mpi_real
        types(3) = mglet_mpi_real

        CALL MPI_Type_create_struct(3, blocklen, disp, types, mpitype)
        CALL MPI_Type_commit(mpitype)
    END SUBROUTINE create_mpi_dtype


    SUBROUTINE bodyforce_reduce(ivec, iovec, length, datatype)
        ! Subroutine arguments
        INTEGER(int32), INTENT(in) :: length
        TYPE(bodyforce_t), INTENT(in) :: ivec(length)
        TYPE(bodyforce_t), INTENT(inout) :: iovec(length)
        TYPE(MPI_Datatype), INTENT(in) :: datatype

        ! Local variables
        INTEGER(intk) :: i

        DO i = 1, length
            iovec(i)%force = iovec(i)%force + ivec(i)%force
            iovec(i)%power = iovec(i)%power + ivec(i)%power
            iovec(i)%torque = iovec(i)%torque + ivec(i)%torque
        END DO
    END SUBROUTINE bodyforce_reduce


    SUBROUTINE bodyforce_reduce_mpi(invec, inoutvec, length, datatype)
        ! Subroutine arguments
        TYPE(C_PTR), VALUE :: invec
        TYPE(C_PTR), VALUE :: inoutvec
        INTEGER(int32) :: length        ! Declaring INTENT makes it incompatible
        TYPE(MPI_Datatype) :: datatype  ! with MPI_f08 MPI_User_function

        ! Local variables
        TYPE(bodyforce_t), POINTER :: idata(:)
        TYPE(bodyforce_t), POINTER :: iodata(:)

        CALL C_F_POINTER(invec, idata, [length])
        CALL C_F_POINTER(inoutvec, iodata, [length])

        CALL bodyforce_reduce(idata, iodata, length, datatype)
    END SUBROUTINE bodyforce_reduce_mpi


    PURE SUBROUTINE zero_bodyforce(force)
        ! Subroutine arguments
        TYPE(bodyforce_t), INTENT(out) :: force(:)

        ! Local variables
        INTEGER(intk) :: i

        DO i = LBOUND(force, 1), UBOUND(force, 1)
            force(i)%force = 0.0
            force(i)%power = 0.0
            force(i)%torque = 0.0
        END DO
    END SUBROUTINE zero_bodyforce


    SUBROUTINE bodyforce_init_logs(dcont)
        ! Subroutine arguments
        LOGICAL, INTENT(in) :: dcont

        ! Local variables
        LOGICAL :: exists
        INTEGER(intk) :: logunit
        CHARACTER(len=mglet_filename_max) :: fname

        IF (myid /= 0) RETURN

        fname = logdir // "/" // bodyforcefile
        INQUIRE(FILE=fname, EXIST=exists)
        IF (.NOT. exists .OR. .NOT. dcont) THEN
            OPEN(NEWUNIT=logunit, FILE=fname)
            WRITE(logunit, '(A9, A18, A6, 12A15)') &
                "#  ITTOT", "TIME", "LEVEL", &
                "F_BODY_TOT_X", "F_BODY_TOT_Y", "F_BODY_TOT_Z", &
                "F_BODY_CONV_X", "F_BODY_CONV_Y", "F_BODY_CONV_Z", &
                "F_BODY_VISC_X", "F_BODY_VISC_Y", "F_BODY_VISC_Z", &
                "F_BODY_PRES_X", "F_BODY_PRES_Y", "F_BODY_PRES_Z"
            CLOSE(logunit)
        END IF

        fname = logdir // "/" // bodytorquefile
        INQUIRE(FILE=fname, EXIST=exists)
        IF (.NOT. exists .OR. .NOT. dcont) THEN
            OPEN(NEWUNIT=logunit, FILE=fname)
            WRITE(logunit, '(A9, A18, A6, 12A15)') &
                "#  ITTOT", "TIME", "LEVEL", &
                "T_BODY_TOT_X", "T_BODY_TOT_Y", "T_BODY_TOT_Z", &
                "T_BODY_CONV_X", "T_BODY_CONV_Y", "T_BODY_CONV_Z", &
                "T_BODY_VISC_X", "T_BODY_VISC_Y", "T_BODY_VISC_Z", &
                "T_BODY_PRES_X", "T_BODY_PRES_Y", "T_BODY_PRES_Z"
            CLOSE(logunit)
        END IF

        fname = logdir // "/" // bodypowerfile
        INQUIRE(FILE=fname, EXIST=exists)
        IF (.NOT. exists .OR. .NOT. dcont) THEN
            OPEN(NEWUNIT=logunit, FILE=fname)
            WRITE(logunit, '(A9, A18, A6, 12A15)') &
                "#  ITTOT", "TIME", "LEVEL", &
                "P_BODY_TOT_X", "P_BODY_TOT_Y", "P_BODY_TOT_Z", &
                "P_BODY_CONV_X", "P_BODY_CONV_Y", "P_BODY_CONV_Z", &
                "P_BODY_VISC_X", "P_BODY_VISC_Y", "P_BODY_VISC_Z", &
                "P_BODY_PRES_X", "P_BODY_PRES_Y", "P_BODY_PRES_Z"
            CLOSE(logunit)
        END IF
    END SUBROUTINE bodyforce_init_logs


    SUBROUTINE bodyforce_print_logs(force, ittot, time)
        ! Subroutine arguments
        TYPE(bodyforce_t), INTENT(in) :: force(minlevel:maxlevel)
        INTEGER(intk), INTENT(in) :: ittot
        REAL(realk), INTENT(in) :: time

        ! Local variables
        INTEGER(intk) :: i, logunit
        CHARACTER(len=mglet_filename_max) :: fname
        REAL(realk) :: tot_x, tot_y, tot_z

        IF (myid /= 0) RETURN

        fname = logdir // "/" // bodyforcefile
        OPEN(NEWUNIT=logunit, FILE=fname, POSITION="APPEND")
        DO i = minlevel, maxlevel
            tot_x = force(i)%force(1, 1) + force(i)%force(1, 2) &
                + force(i)%force(1, 3)
            tot_y = force(i)%force(2, 1) + force(i)%force(2, 2) &
                + force(i)%force(2, 3)
            tot_z = force(i)%force(3, 1) + force(i)%force(3, 2) &
                + force(i)%force(3, 3)

            WRITE(logunit, '(I9, ES18.10, I6, 12ES15.5)') &
                ittot, time, i, &
                tot_x, &
                tot_y, &
                tot_z, &
                force(i)%force(1, 1), &
                force(i)%force(2, 1), &
                force(i)%force(3, 1), &
                force(i)%force(1, 2), &
                force(i)%force(2, 2), &
                force(i)%force(3, 2), &
                force(i)%force(1, 3), &
                force(i)%force(2, 3), &
                force(i)%force(3, 3)
        END DO
        CLOSE(logunit)

        fname = logdir // "/" // bodypowerfile
        OPEN(NEWUNIT=logunit, FILE=fname, POSITION="APPEND")
        DO i = minlevel, maxlevel
            tot_x = force(i)%power(1, 1) + force(i)%power(1, 2) &
                + force(i)%power(1, 3)
            tot_y = force(i)%power(2, 1) + force(i)%power(2, 2) &
                + force(i)%power(2, 3)
            tot_z = force(i)%power(3, 1) + force(i)%power(3, 2) &
                + force(i)%power(3, 3)

            WRITE(logunit, '(I9, ES18.10, I6, 12ES15.5)') &
                ittot, time, i, &
                tot_x, &
                tot_y, &
                tot_z, &
                force(i)%power(1, 1), &
                force(i)%power(2, 1), &
                force(i)%power(3, 1), &
                force(i)%power(1, 2), &
                force(i)%power(2, 2), &
                force(i)%power(3, 2), &
                force(i)%power(1, 3), &
                force(i)%power(2, 3), &
                force(i)%power(3, 3)
        END DO
        CLOSE(logunit)

        fname = logdir // "/" // bodytorquefile
        OPEN(NEWUNIT=logunit, FILE=fname, POSITION="APPEND")
        DO i = minlevel, maxlevel
            tot_x = force(i)%torque(1, 1) + force(i)%torque(1, 2) &
                + force(i)%torque(1, 3)
            tot_y = force(i)%torque(2, 1) + force(i)%torque(2, 2) &
                + force(i)%torque(2, 3)
            tot_z = force(i)%torque(3, 1) + force(i)%torque(3, 2) &
                + force(i)%torque(3, 3)

            WRITE(logunit, '(I9, ES18.10, I6, 12ES15.5)') &
                ittot, time, i, &
                tot_x, &
                tot_y, &
                tot_z, &
                force(i)%torque(1, 1), &
                force(i)%torque(2, 1), &
                force(i)%torque(3, 1), &
                force(i)%torque(1, 2), &
                force(i)%torque(2, 2), &
                force(i)%torque(3, 2), &
                force(i)%torque(1, 3), &
                force(i)%torque(2, 3), &
                force(i)%torque(3, 3)
        END DO
        CLOSE(logunit)
    END SUBROUTINE bodyforce_print_logs
END MODULE gc_compbodyforce_mod
