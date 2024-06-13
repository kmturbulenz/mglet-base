MODULE setboundarybuffers_mod
    USE core_mod
    USE ib_mod, ONLY: ib
    USE flowcore_mod

    IMPLICIT NONE(type, external)
    PRIVATE

    ! Bound operation 'T' operate on U, V, W, P
    TYPE, EXTENDS(bound_t) :: setboundarybuffers_t
    CONTAINS
        PROCEDURE, NOPASS :: front => bfront
        PROCEDURE, NOPASS :: back => bfront
        PROCEDURE, NOPASS :: right => bright
        PROCEDURE, NOPASS :: left => bright
        PROCEDURE, NOPASS :: bottom => bbottom
        PROCEDURE, NOPASS :: top => bbottom
    END TYPE setboundarybuffers_t
    TYPE(setboundarybuffers_t) :: setboundarybuffers

    PUBLIC :: setboundarybuffers

CONTAINS
    SUBROUTINE bfront(igrid, iface, ibocd, ctyp, f1, f2, f3, f4, timeph)
        ! Subroutine arguments
        INTEGER(intk), INTENT(in) :: igrid, iface, ibocd
        CHARACTER(len=*), INTENT(in) :: ctyp
        TYPE(field_t), INTENT(inout) :: f1
        TYPE(field_t), INTENT(inout), OPTIONAL :: f2, f3, f4
        REAL(realk), INTENT(in), OPTIONAL :: timeph

        ! Local variables
        REAL(realk), POINTER, CONTIGUOUS :: xstag(:), y(:), z(:)
        REAL(realk), POINTER, CONTIGUOUS :: dx(:), dy(:), dz(:)
        REAL(realk), POINTER, CONTIGUOUS :: ddy(:), ddz(:)
        REAL(realk), POINTER, CONTIGUOUS :: ubuf(:, :, :), vbuf(:, :, :), &
            wbuf(:, :, :), bubuf(:, :, :)
        REAL(realk) :: xbuf(2), dxbuf(2)
        TYPE(field_t), POINTER :: bu

        ! Only works on FIX and OP1 boundaries, should do nothing otherwise
        ! For OP1 the velocity set is the velocity that is prescribed when
        ! there are backflow/inflow through the boundary.
        SELECT CASE (ctyp)
        CASE ("FIX", "OP1")
            CONTINUE
        CASE DEFAULT
            RETURN
        END SELECT

        ! Assure that required fields are present
        IF (.NOT. PRESENT(f2) .OR. .NOT. PRESENT(f3)) THEN
            CALL errr(__FILE__, __LINE__)
        END IF
        IF (.NOT. PRESENT(timeph)) CALL errr(__FILE__, __LINE__)

        ! Fetch pointers
        CALL f1%buffers%get_buffer(ubuf, igrid, iface)
        CALL f2%buffers%get_buffer(vbuf, igrid, iface)
        CALL f3%buffers%get_buffer(wbuf, igrid, iface)

        IF (uinf_is_expr) THEN
            CALL get_fieldptr(xstag, "XSTAG", igrid)
            CALL get_fieldptr(y, "Y", igrid)
            CALL get_fieldptr(z, "Z", igrid)

            CALL get_fieldptr(dx, "DX", igrid)
            CALL get_fieldptr(dy, "DY", igrid)
            CALL get_fieldptr(dz, "DZ", igrid)

            CALL get_fieldptr(ddy, "DDY", igrid)
            CALL get_fieldptr(ddz, "DDZ", igrid)

            SELECT CASE (iface)
            CASE (1)
                ! Front
                xbuf = xstag(2)
                dxbuf = dx(2)
            CASE (2)
                ! Back
                xbuf = xstag(SIZE(xstag)-2)
                dxbuf = dx(SIZE(dx)-2)
            CASE DEFAULT
                CALL errr(__FILE__, __LINE__)
            END SELECT

            CALL initial_condition(ubuf, "u", uinf_expr(1), rho, gmol, &
                tu_level, timeph, xbuf, y, z, dxbuf, dy, dz, dxbuf, ddy, ddz)
            CALL initial_condition(vbuf, "v", uinf_expr(2), rho, gmol, &
                tu_level, timeph, xbuf, y, z, dxbuf, dy, dz, dxbuf, ddy, ddz)
            CALL initial_condition(wbuf, "w", uinf_expr(3), rho, gmol, &
                tu_level, timeph, xbuf, y, z, dxbuf, dy, dz, dxbuf, ddy, ddz)
        ELSE
            ubuf = uinf(1)
            vbuf = uinf(2)
            wbuf = uinf(3)
        END IF

        ! Without IB or if the boundary condition is not FIX, we are done
        IF (ib%type == "NONE" .OR. ctyp /= "FIX") RETURN

        ! This multiplies with the open fraction of a cell in the FIX boundary
        ! condition. This is done to ensure that the inflow is divergence free
        CALL get_field(bu, "BU")
        CALL bu%buffers%get_buffer(bubuf, igrid, iface)
        ubuf = ubuf*bubuf
    END SUBROUTINE bfront


    SUBROUTINE bright(igrid, iface, ibocd, ctyp, f1, f2, f3, f4, timeph)
        ! Subroutine arguments
        INTEGER(intk), INTENT(in) :: igrid, iface, ibocd
        CHARACTER(len=*), INTENT(in) :: ctyp
        TYPE(field_t), INTENT(inout) :: f1
        TYPE(field_t), INTENT(inout), OPTIONAL :: f2, f3, f4
        REAL(realk), INTENT(in), OPTIONAL :: timeph

        ! Local variables
        REAL(realk), POINTER, CONTIGUOUS :: x(:), ystag(:), z(:)
        REAL(realk), POINTER, CONTIGUOUS :: dx(:), dy(:), dz(:)
        REAL(realk), POINTER, CONTIGUOUS :: ddx(:), ddz(:)
        REAL(realk), POINTER, CONTIGUOUS :: ubuf(:, :, :), vbuf(:, :, :), &
            wbuf(:, :, :), bvbuf(:, :, :)
        REAL(realk) :: ybuf(2), dybuf(2)
        TYPE(field_t), POINTER :: bv

        ! Only works on FIX and OP1 boundaries, should do nothing otherwise
        ! For OP1 the velocity set is the velocity that is prescribed when
        ! there are backflow/inflow through the boundary.
        SELECT CASE (ctyp)
        CASE ("FIX", "OP1")
            CONTINUE
        CASE DEFAULT
            RETURN
        END SELECT

        ! Assure that required fields are present
        IF (.NOT. PRESENT(f2) .OR. .NOT. PRESENT(f3)) THEN
            CALL errr(__FILE__, __LINE__)
        END IF
        IF (.NOT. PRESENT(timeph)) CALL errr(__FILE__, __LINE__)

        ! Fetch pointers
        CALL f1%buffers%get_buffer(ubuf, igrid, iface)
        CALL f2%buffers%get_buffer(vbuf, igrid, iface)
        CALL f3%buffers%get_buffer(wbuf, igrid, iface)

        IF (uinf_is_expr) THEN
            CALL get_fieldptr(x, "X", igrid)
            CALL get_fieldptr(ystag, "YSTAG", igrid)
            CALL get_fieldptr(z, "Z", igrid)

            CALL get_fieldptr(dx, "DX", igrid)
            CALL get_fieldptr(dy, "DY", igrid)
            CALL get_fieldptr(dz, "DZ", igrid)

            CALL get_fieldptr(ddx, "DDX", igrid)
            CALL get_fieldptr(ddz, "DDZ", igrid)

            SELECT CASE (iface)
            CASE (3)
                ! Right
                ybuf = ystag(2)
                dybuf = dy(2)
            CASE (4)
                ! Left
                ybuf = ystag(SIZE(ystag)-2)
                dybuf = dy(SIZE(dy)-2)
            CASE DEFAULT
                CALL errr(__FILE__, __LINE__)
            END SELECT

            CALL initial_condition(ubuf, "u", uinf_expr(1), rho, gmol, &
                tu_level, timeph, ybuf, x, z, dybuf, dx, dz, dybuf, ddx, ddz)
            CALL initial_condition(vbuf, "v", uinf_expr(2), rho, gmol, &
                tu_level, timeph, ybuf, x, z, dybuf, dx, dz, dybuf, ddx, ddz)
            CALL initial_condition(wbuf, "w", uinf_expr(3), rho, gmol, &
                tu_level, timeph, ybuf, x, z, dybuf, dx, dz, dybuf, ddx, ddz)
        ELSE
            ubuf = uinf(1)
            vbuf = uinf(2)
            wbuf = uinf(3)
        END IF

        ! Without IB or if the boundary condition is not FIX, we are done
        IF (ib%type == "NONE" .OR. ctyp /= "FIX") RETURN

        ! This multiplies with the open fraction of a cell in the FIX boundary
        ! condition. This is done to ensure that the inflow is divergence free
        CALL get_field(bv, "BV")
        CALL bv%buffers%get_buffer(bvbuf, igrid, iface)
        vbuf = vbuf*bvbuf
    END SUBROUTINE bright


    SUBROUTINE bbottom(igrid, iface, ibocd, ctyp, f1, f2, f3, f4, timeph)
        ! Subroutine arguments
        INTEGER(intk), INTENT(in) :: igrid, iface, ibocd
        CHARACTER(len=*), INTENT(in) :: ctyp
        TYPE(field_t), INTENT(inout) :: f1
        TYPE(field_t), INTENT(inout), OPTIONAL :: f2, f3, f4
        REAL(realk), INTENT(in), OPTIONAL :: timeph

        ! Local variables
        REAL(realk), POINTER, CONTIGUOUS :: x(:), y(:), zstag(:)
        REAL(realk), POINTER, CONTIGUOUS :: dx(:), dy(:), dz(:)
        REAL(realk), POINTER, CONTIGUOUS :: ddx(:), ddy(:)
        REAL(realk), POINTER, CONTIGUOUS :: ubuf(:, :, :), vbuf(:, :, :), &
            wbuf(:, :, :), bwbuf(:, :, :)
        REAL(realk) :: zbuf(2), dzbuf(2)
        TYPE(field_t), POINTER :: bw

        ! Only works on FIX and OP1 boundaries, should do nothing otherwise
        ! For OP1 the velocity set is the velocity that is prescribed when
        ! there are backflow/inflow through the boundary.
        SELECT CASE (ctyp)
        CASE ("FIX", "OP1")
            CONTINUE
        CASE DEFAULT
            RETURN
        END SELECT

        ! Assure that required fields are present
        IF (.NOT. PRESENT(f2) .OR. .NOT. PRESENT(f3)) THEN
            CALL errr(__FILE__, __LINE__)
        END IF
        IF (.NOT. PRESENT(timeph)) CALL errr(__FILE__, __LINE__)

        ! Fetch pointers
        CALL f1%buffers%get_buffer(ubuf, igrid, iface)
        CALL f2%buffers%get_buffer(vbuf, igrid, iface)
        CALL f3%buffers%get_buffer(wbuf, igrid, iface)

        IF (uinf_is_expr) THEN
            CALL get_fieldptr(x, "X", igrid)
            CALL get_fieldptr(y, "Y", igrid)
            CALL get_fieldptr(zstag, "ZSTAG", igrid)

            CALL get_fieldptr(dx, "DX", igrid)
            CALL get_fieldptr(dy, "DY", igrid)
            CALL get_fieldptr(dz, "DZ", igrid)

            CALL get_fieldptr(ddx, "DDX", igrid)
            CALL get_fieldptr(ddy, "DDY", igrid)

            SELECT CASE (iface)
            CASE (5)
                ! Bottom
                zbuf = zstag(2)
                dzbuf = dz(2)
            CASE (6)
                ! Top
                zbuf = zstag(SIZE(zstag)-2)
                dzbuf = dz(SIZE(dz)-2)
            CASE DEFAULT
                CALL errr(__FILE__, __LINE__)
            END SELECT

            CALL initial_condition(ubuf, "u", uinf_expr(1), rho, gmol, &
                tu_level, timeph, zbuf, x, y, dzbuf, dx, dy, dzbuf, ddx, ddy)
            CALL initial_condition(vbuf, "v", uinf_expr(2), rho, gmol, &
                tu_level, timeph, zbuf, x, y, dzbuf, dx, dy, dzbuf, ddx, ddy)
            CALL initial_condition(wbuf, "w", uinf_expr(3), rho, gmol, &
                tu_level, timeph, zbuf, x, y, dzbuf, dx, dy, dzbuf, ddx, ddy)
        ELSE
            ubuf = uinf(1)
            vbuf = uinf(2)
            wbuf = uinf(3)
        END IF

        ! Without IB or if the boundary condition is not FIX, we are done
        IF (ib%type == "NONE" .OR. ctyp /= "FIX") RETURN

        ! This multiplies with the open fraction of a cell in the FIX boundary
        ! condition. This is done to ensure that the inflow is divergence free
        CALL get_field(bw, "BW")
        CALL bw%buffers%get_buffer(bwbuf, igrid, iface)
        wbuf = wbuf*bwbuf
    END SUBROUTINE bbottom
END MODULE setboundarybuffers_mod
