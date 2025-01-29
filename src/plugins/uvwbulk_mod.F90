MODULE uvwbulk_mod
    USE, INTRINSIC :: ISO_C_BINDING, ONLY: C_PTR, C_F_POINTER
    USE core_mod
    USE MPI_f08
    USE ib_mod, ONLY: ib

    IMPLICIT NONE(type, external)
    PRIVATE

    TYPE, BIND(C) :: uvwbulk_t
        REAL(c_realk) :: ubulk = 0.0
        REAL(c_realk) :: vbulk = 0.0
        REAL(c_realk) :: wbulk = 0.0
        REAL(c_realk) :: uubulk = 0.0
        REAL(c_realk) :: vvbulk = 0.0
        REAL(c_realk) :: wwbulk = 0.0
        REAL(c_realk) :: volu = 0.0
        REAL(c_realk) :: volv = 0.0
        REAL(c_realk) :: volw = 0.0
    END TYPE uvwbulk_t
    INTEGER(int32), PARAMETER :: uvwbulk_elems = 9

    LOGICAL :: enable_uvwbulk

    TYPE(MPI_Datatype) :: mpitype
    TYPE(MPI_Op) :: mpiop

    CHARACTER(len=*), PARAMETER :: logfile = logdir//"/uvwbulk.log"

#ifdef _MGLET_DOUBLE_PRECISION_
    CHARACTER(len=*), PARAMETER :: head_fmt = "(A9, A18, 6A23)"
    CHARACTER(len=*), PARAMETER :: num_fmt = "(I9, ES18.10, 6ES23.14)"
#else
    CHARACTER(len=*), PARAMETER :: head_fmt = "(A9, A18, 6A14)"
    CHARACTER(len=*), PARAMETER :: num_fmt = "(I9, ES18.10, 6ES14.5)"
#endif

    PUBLIC :: init_uvwbulk, itinfo_uvwbulk, finish_uvwbulk

CONTAINS
    SUBROUTINE init_uvwbulk(ittot, mtstep, itint, timeph, dt, tend)
        ! Subroutine arguments
        INTEGER(intk), INTENT(in) :: ittot
        INTEGER(intk), INTENT(in) :: mtstep
        INTEGER(intk), INTENT(in) :: itint
        REAL(realk), INTENT(in) :: timeph
        REAL(realk), INTENT(in) :: dt
        REAL(realk), INTENT(in) :: tend

        ! Local variables
        LOGICAL :: exists
        INTEGER :: logunit

        enable_uvwbulk = .FALSE.
        IF (.NOT. fort7%exists("/uvwbulk")) THEN
            RETURN
        END IF
        CALL fort7%get_value("/uvwbulk", enable_uvwbulk)
        IF (.NOT. enable_uvwbulk) RETURN

        IF (ib%type == "CUTCELL") CALL errr(__FILE__, __LINE__)

        IF (myid == 0) THEN
            INQUIRE(FILE=logfile, EXIST=exists)
            IF (.NOT. exists .OR. .NOT. dcont) THEN
                OPEN(NEWUNIT=logunit, FILE=logfile)
                WRITE(logunit, head_fmt) &
                    "#  ITTOT", "TIME", "UBULK", "VBULK", "WBULK", "UUBULK", &
                    "VVBULK", "WWBULK"
                CLOSE(logunit)
            END IF
        END IF

        ! MPI type and reduction operation
        CALL MPI_Type_contiguous(uvwbulk_elems, mglet_mpi_real, mpitype)
        CALL MPI_Type_commit(mpitype)
        CALL MPI_Op_create(uvwbulk_reduce, .TRUE., mpiop)
    END SUBROUTINE init_uvwbulk


    SUBROUTINE finish_uvwbulk()
        IF (.NOT. enable_uvwbulk) RETURN

        CALL MPI_Type_free(mpitype)
        CALL MPI_Op_free(mpiop)
    END SUBROUTINE finish_uvwbulk


    SUBROUTINE itinfo_uvwbulk(itstep, ittot, timeph, dt)
        ! Subroutine arguments
        INTEGER(intk), INTENT(in) :: itstep
        INTEGER(intk), INTENT(in) :: ittot
        REAL(realk), INTENT(in) :: timeph
        REAL(realk), INTENT(in) :: dt

        ! Local variables
        INTEGER :: logunit
        INTEGER(intk) :: i, igrid
        TYPE(uvwbulk_t), TARGET :: uvwbulk, uvwbulk_this
        TYPE(field_t), POINTER :: u_f, v_f, w_f, finecell_f, &
            dx_f, dy_f, dz_f, ddx_f, ddy_f, ddz_f
        REAL(realk), POINTER, CONTIGUOUS :: u(:, :, :), v(:, :, :), &
            w(:, :, :), finecell(:, :, :)
        REAL(realk), POINTER, CONTIGUOUS :: dx(:), dy(:), dz(:), ddx(:), &
            ddy(:), ddz(:)

        IF (.NOT. enable_uvwbulk) RETURN

        CALL get_field(u_f, "U")
        CALL get_field(v_f, "V")
        CALL get_field(w_f, "W")
        CALL get_field(finecell_f, "FINECELL")

        CALL get_field(dx_f, "DX")
        CALL get_field(dy_f, "DY")
        CALL get_field(dz_f, "DZ")

        CALL get_field(ddx_f, "DDX")
        CALL get_field(ddy_f, "DDY")
        CALL get_field(ddz_f, "DDZ")

        DO i = 1, nmygrids
            igrid = mygrids(i)

            CALL u_f%get_ptr(u, igrid)
            CALL v_f%get_ptr(v, igrid)
            CALL w_f%get_ptr(w, igrid)
            CALL finecell_f%get_ptr(finecell, igrid)

            CALL dx_f%get_ptr(dx, igrid)
            CALL dy_f%get_ptr(dy, igrid)
            CALL dz_f%get_ptr(dz, igrid)

            CALL ddx_f%get_ptr(ddx, igrid)
            CALL ddy_f%get_ptr(ddy, igrid)
            CALL ddz_f%get_ptr(ddz, igrid)

            ! Compute bulk velocities
            CALL uvwbulk_grid(uvwbulk, u, v, w, finecell, dx, dy, dz, &
                ddx, ddy, ddz)

            ! Instantly update local reduction in uvwbulk_this
            CALL uvwbulk_reduce(C_LOC(uvwbulk), C_LOC(uvwbulk_this), &
                1, mpitype)
        END DO

        ! Each process has their local reductions in uvwbulk_this, do a global
        ! reduction into uvwbulk
        CALL MPI_Reduce(uvwbulk_this, uvwbulk, 1, mpitype, mpiop, &
            0, MPI_COMM_WORLD)

        ! Write logfile
        IF (myid == 0) THEN
            OPEN(NEWUNIT=logunit, FILE=logfile, POSITION="APPEND")
            WRITE(logunit, num_fmt) ittot, timeph, &
                uvwbulk%ubulk, uvwbulk%vbulk, uvwbulk%wbulk, &
                uvwbulk%uubulk, uvwbulk%vvbulk, uvwbulk%wwbulk
            CLOSE(logunit)
        END IF
    END SUBROUTINE itinfo_uvwbulk


    SUBROUTINE uvwbulk_grid(uvwbulk, u, v, w, finecell, dx, dy, dz, &
            ddx, ddy, ddz)
        ! Subroutine arguments
        TYPE(uvwbulk_t), INTENT(out) :: uvwbulk
        REAL(realk), INTENT(in), CONTIGUOUS :: u(:, :, :), v(:, :, :), &
            w(:, :, :), finecell(:, :, :)
        REAL(realk), INTENT(in), CONTIGUOUS :: dx(:), dy(:), dz(:), ddx(:), &
            ddy(:), ddz(:)

        ! Local variables
        INTEGER(intk) :: kk, jj, ii
        INTEGER(intk) :: k, j, i

        REAL(real64) :: ucell, vcell, wcell
        REAL(real64) :: uucell, vvcell, wwcell
        REAL(real64) :: volucell, volvcell, volwcell

        REAL(real64) :: ubulk, vbulk, wbulk
        REAL(real64) :: uubulk, vvbulk, wwbulk
        REAL(real64) :: volu, volv, volw

        kk = SIZE(u, 1)
        jj = SIZE(u, 2)
        ii = SIZE(u, 3)

        IF (SIZE(dx) /= ii) ERROR STOP
        IF (SIZE(dy) /= jj) ERROR STOP
        IF (SIZE(dz) /= kk) ERROR STOP
        IF (SIZE(ddx) /= ii) ERROR STOP
        IF (SIZE(ddy) /= jj) ERROR STOP
        IF (SIZE(ddz) /= kk) ERROR STOP

        ! Temporary summation in double precision to avoid rounding errors
        ! when iterating over individual cells
        ubulk = 0.0
        vbulk = 0.0
        wbulk = 0.0

        uubulk = 0.0
        vvbulk = 0.0
        wwbulk = 0.0

        volu = 0.0
        volv = 0.0
        volw = 0.0

        DO i = 3, ii-2
            DO j = 3, jj-2
                DO k = 3, kk-2
                    volucell = finecell(k, j, i)*dx(i)*ddy(j)*ddz(k)
                    volvcell = finecell(k, j, i)*ddx(i)*dy(j)*ddz(k)
                    volwcell = finecell(k, j, i)*ddx(i)*ddy(j)*dz(k)

                    ucell = u(k, j, i)*volucell
                    vcell = v(k, j, i)*volvcell
                    wcell = w(k, j, i)*volwcell

                    uucell = u(k, j, i)*u(k, j, i)*volucell
                    vvcell = v(k, j, i)*v(k, j, i)*volvcell
                    wwcell = w(k, j, i)*w(k, j, i)*volwcell

                    ubulk = ubulk + ucell
                    vbulk = vbulk + vcell
                    wbulk = wbulk + wcell

                    uubulk = uubulk + uucell
                    vvbulk = vvbulk + vvcell
                    wwbulk = wwbulk + wwcell

                    volu = volu + volucell
                    volv = volv + volvcell
                    volw = volw + volwcell
                END DO
            END DO
        END DO

        ! The resulting values for the complete grid is casted back to working
        ! precision
        uvwbulk%ubulk = REAL(divide0(ubulk, volu), realk)
        uvwbulk%vbulk = REAL(divide0(vbulk, volv), realk)
        uvwbulk%wbulk = REAL(divide0(wbulk, volw), realk)

        uvwbulk%uubulk = REAL(divide0(uubulk, volu), realk)
        uvwbulk%vvbulk = REAL(divide0(vvbulk, volv), realk)
        uvwbulk%wwbulk = REAL(divide0(wwbulk, volw), realk)

        uvwbulk%volu = REAL(volu, realk)
        uvwbulk%volv = REAL(volv, realk)
        uvwbulk%volw = REAL(volw, realk)
    END SUBROUTINE uvwbulk_grid


    SUBROUTINE uvwbulk_reduce(invec, inoutvec, length, datatype)
        ! Subroutine arguments
        TYPE(C_PTR), VALUE :: invec
        TYPE(C_PTR), VALUE :: inoutvec
        INTEGER(int32) :: length        ! Declaring INTENT makes it incompatible
        TYPE(MPI_Datatype) :: datatype  ! with MPI_f08 MPI_User_function

        ! Local variables
        INTEGER(int32) :: i
        TYPE(uvwbulk_t), POINTER :: idata(:)
        TYPE(uvwbulk_t), POINTER :: iodata(:)

        CALL C_F_POINTER(invec, idata, [length])
        CALL C_F_POINTER(inoutvec, iodata, [length])

        DO i = 1, length
            iodata(i)%ubulk = divide0((iodata(i)%ubulk*iodata(i)%volu + &
                idata(i)%ubulk*idata(i)%volu), (iodata(i)%volu + idata(i)%volu))
            iodata(i)%vbulk = divide0((iodata(i)%vbulk*iodata(i)%volv + &
                idata(i)%vbulk*idata(i)%volv), (iodata(i)%volv + idata(i)%volv))
            iodata(i)%wbulk = divide0((iodata(i)%wbulk*iodata(i)%volw + &
                idata(i)%wbulk*idata(i)%volw), (iodata(i)%volw + idata(i)%volw))

            iodata(i)%uubulk = divide0((iodata(i)%uubulk*iodata(i)%volu + &
                idata(i)%uubulk*idata(i)%volu), &
                (iodata(i)%volu + idata(i)%volu))
            iodata(i)%vvbulk = divide0((iodata(i)%vvbulk*iodata(i)%volv + &
                idata(i)%vvbulk*idata(i)%volv), &
                (iodata(i)%volv + idata(i)%volv))
            iodata(i)%wwbulk = divide0((iodata(i)%wwbulk*iodata(i)%volw + &
                idata(i)%wwbulk*idata(i)%volw), &
                (iodata(i)%volw + idata(i)%volw))

            iodata(i)%volu = iodata(i)%volu + idata(i)%volu
            iodata(i)%volv = iodata(i)%volv + idata(i)%volv
            iodata(i)%volw = iodata(i)%volw + idata(i)%volw
        END DO
    END SUBROUTINE uvwbulk_reduce

END MODULE uvwbulk_mod
