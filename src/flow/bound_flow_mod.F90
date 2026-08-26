MODULE bound_flow_mod
    USE core_mod
    USE flowcore_mod

    IMPLICIT NONE(type, external)
    PRIVATE

    ! Bound operation 'T' operate on U, V, W, P
    TYPE, EXTENDS(bound_t) :: bound_flow_t
    CONTAINS
        PROCEDURE, NOPASS :: front => bfront
        PROCEDURE, NOPASS :: back => bfront
        PROCEDURE, NOPASS :: right => bright
        PROCEDURE, NOPASS :: left => bright
        PROCEDURE, NOPASS :: bottom => bbottom
        PROCEDURE, NOPASS :: top => bbottom
    END TYPE bound_flow_t
    TYPE(bound_flow_t) :: bound_flow

    INTEGER(intk), PARAMETER :: boundtasksize = 3
    INTEGER(intk), ALLOCATABLE :: boundtasks(:, :, :)
    INTEGER(intk), ALLOCATABLE :: nboundtaskslvl(:)
    REAL(realk), ALLOCATABLE :: boundtaskpinf(:, :)
    !$omp declare target(boundtasks, boundtaskpinf)

    PUBLIC :: bound_flow, init_bound_flow, apply_bound_flow, finish_bound_flow

    INTERFACE
        SUBROUTINE bound_flow_device_c(kk, jj, ii, iface, ityp, pinf_in, u, &
            v, w, p, bp, ubuf, vbuf, wbuf, ddx, ddy, ddz, rho_in) BIND(C, &
            name="bound_flow_device_c")
            IMPORT :: c_intk, c_realk
            INTEGER(c_intk), VALUE, INTENT(in) :: kk, jj, ii, iface, ityp
            REAL(c_realk), VALUE, INTENT(in) :: pinf_in, rho_in
            REAL(c_realk), INTENT(inout) :: u(*), v(*), w(*), p(*)
            REAL(c_realk), INTENT(in) :: bp(*)
            REAL(c_realk), INTENT(inout) :: ubuf(*), vbuf(*), wbuf(*)
            REAL(c_realk), INTENT(in) :: ddx(*), ddy(*), ddz(*)
        END SUBROUTINE bound_flow_device_c
    END INTERFACE

    !$omp declare target(bound_flow_device_c)

CONTAINS
    SUBROUTINE init_bound_flow()
        INTEGER(intk) :: nlevels, ilevel, ilevel_index
        INTEGER(intk) :: imygrid, igrid, iface, ibocd, nbocd, itask, ityp
        INTEGER(intk) :: nop1
        REAL(realk) :: pinf(1)
        CHARACTER(len=8) :: ctyp
        TYPE(field_t), POINTER :: u_f, v_f, w_f

        nlevels = maxlevel - minlevel + 1
        ALLOCATE(nboundtaskslvl(nlevels), source=0_intk)

        DO ilevel = minlevel, maxlevel
            ilevel_index = ilevel - minlevel + 1
            DO imygrid = 1, nmygridslvl(ilevel)
                igrid = mygridslvl(imygrid, ilevel)
                DO iface = 1, 6
                    nbocd = nboconds(iface, igrid)
                    DO ibocd = 1, nbocd
                        CALL get_bc_ctyp(ctyp, ibocd, iface, igrid)
                        IF (bound_type(ctyp) == 0) CYCLE
                        nboundtaskslvl(ilevel_index) = &
                            nboundtaskslvl(ilevel_index) + 1
                    END DO
                END DO
            END DO
        END DO

        ALLOCATE(boundtasks(boundtasksize, MAXVAL(nboundtaskslvl), nlevels), &
            source=-1_intk)
        ALLOCATE(boundtaskpinf(MAXVAL(nboundtaskslvl), nlevels), &
            source=0.0_realk)

        DO ilevel = minlevel, maxlevel
            ilevel_index = ilevel - minlevel + 1
            itask = 0
            DO imygrid = 1, nmygridslvl(ilevel)
                igrid = mygridslvl(imygrid, ilevel)
                DO iface = 1, 6
                    nbocd = nboconds(iface, igrid)
                    DO ibocd = 1, nbocd
                        CALL get_bc_ctyp(ctyp, ibocd, iface, igrid)
                        ityp = bound_type(ctyp)
                        IF (ityp == 0) CYCLE
                        itask = itask + 1
                        boundtasks(:, itask, ilevel_index) = &
                            [igrid, iface, ityp]
                        IF (ityp == 2) THEN
                            CALL get_nbcprms(igrid, iface, ibocd, nreal=nop1)
                            IF (nop1 == 0) THEN
                                pinf(1) = 0.0_realk
                            ELSE IF (nop1 == 1) THEN
                                CALL get_bcprms(pinf, igrid, iface, ibocd)
                            ELSE
                                CALL errr(__FILE__, __LINE__)
                            END IF
                            boundtaskpinf(itask, ilevel_index) = pinf(1)
                        END IF
                    END DO
                END DO
            END DO
        END DO

        !$omp target enter data map(always, to: boundtasks, boundtaskpinf)
        CALL get_field(u_f, "U")
        CALL get_field(v_f, "V")
        CALL get_field(w_f, "W")
    END SUBROUTINE init_bound_flow


    INTEGER(intk) FUNCTION bound_type(ctyp)
        CHARACTER(len=*), INTENT(in) :: ctyp

        SELECT CASE (ctyp)
        CASE ("FIX")
            bound_type = 1
        CASE ("OP1")
            bound_type = 2
        CASE ("NOS")
            bound_type = 3
        CASE ("SLI")
            bound_type = 4
        CASE ("PAR")
            bound_type = 5
        CASE DEFAULT
            bound_type = 0
        END SELECT
    END FUNCTION bound_type


    SUBROUTINE apply_bound_flow(ilevel, u_f, v_f, w_f, p_f)
        INTEGER(intk), INTENT(in) :: ilevel
        TYPE(field_t), INTENT(inout) :: u_f, v_f, w_f
        TYPE(field_t), OPTIONAL, INTENT(inout) :: p_f

        TYPE(field_t), POINTER :: bp_f, ddx_f, ddy_f, ddz_f
        INTEGER(intk) :: ilevel_index, ntasks

        ilevel_index = ilevel - minlevel + 1
        ntasks = nboundtaskslvl(ilevel_index)
        CALL get_field(bp_f, "BP")
        CALL get_field(ddx_f, "DDX")
        CALL get_field(ddy_f, "DDY")
        CALL get_field(ddz_f, "DDZ")

        CALL apply_bound_flow_impl(ilevel_index, ntasks, u_f%arr, v_f%arr, &
            w_f%arr, p_f%arr, u_f%buffers, v_f%buffers, w_f%buffers, &
            bp_f%arr, ddx_f%arr, ddy_f%arr, ddz_f%arr)

    END SUBROUTINE apply_bound_flow


    SUBROUTINE apply_bound_flow_impl(ilevel_index, ntasks, u, v, w, p, &
            ubuffer, vbuffer, wbuffer, bp, ddx, ddy, ddz)
        INTEGER(intk), INTENT(in) :: ilevel_index, ntasks
        REAL(realk), INTENT(inout) :: u(*), v(*), w(*), p(*)
        REAL(realk), INTENT(inout) :: ubuffer(*), vbuffer(*), wbuffer(*)
        REAL(realk), INTENT(in) :: bp(*), ddx(*), ddy(*), ddz(*)

        INTEGER(intk) :: itask, igrid, iface, ityp, kk, jj, ii
        INTEGER(intk) :: ip3, ipx, ipy, ipz, ipbb
        REAL(realk) :: pinf

        IF (ntasks == 0) RETURN

        !$omp target teams distribute private(itask, igrid, iface, ityp, &
        !$omp& kk, jj, ii, ip3, ipx, ipy, ipz, ipbb, pinf)
        DO itask = 1, ntasks

            !$omp parallel
            igrid = boundtasks(1, itask, ilevel_index)
            iface = boundtasks(2, itask, ilevel_index)
            ityp = boundtasks(3, itask, ilevel_index)
            pinf = boundtaskpinf(itask, ilevel_index)

            CALL get_mgdims(kk, jj, ii, igrid)
            CALL get_ip3(ip3, igrid)
            CALL get_ip1x(ipx, igrid)
            CALL get_ip1y(ipy, igrid)
            CALL get_ip1z(ipz, igrid)
            CALL get_ipbb(ipbb, iface, igrid)

            CALL bound_flow_device_c(INT(kk, c_intk), INT(jj, c_intk), &
                INT(ii, c_intk), INT(iface, c_intk), INT(ityp, c_intk), &
                REAL(pinf, c_realk), u(ip3), v(ip3), w(ip3), p(ip3), &
                bp(ip3), ubuffer(ipbb), vbuffer(ipbb), wbuffer(ipbb), &
                ddx(ipx), ddy(ipy), ddz(ipz), REAL(rho, c_realk))
            !$omp end parallel

        END DO
        !$omp end target teams distribute

    END SUBROUTINE apply_bound_flow_impl


    SUBROUTINE finish_bound_flow()
        !$omp target exit data map(always, delete: boundtasks, boundtaskpinf)
        DEALLOCATE(boundtaskpinf)
        DEALLOCATE(boundtasks)
        DEALLOCATE(nboundtaskslvl)
    END SUBROUTINE finish_bound_flow


    SUBROUTINE bfront(igrid, iface, ibocd, ctyp, f1, f2, f3, f4, timeph)
        ! Subroutine arguments
        INTEGER(intk), INTENT(in) :: igrid, iface, ibocd
        CHARACTER(len=*), INTENT(in) :: ctyp
        TYPE(field_t), INTENT(inout) :: f1
        TYPE(field_t), INTENT(inout), OPTIONAL :: f2, f3, f4
        REAL(realk), INTENT(in), OPTIONAL :: timeph

        ! Local variables
        INTEGER(intk) :: kk, jj, ii
        INTEGER(intk) :: k, j, i2, i3, i4, istag1, istag2, dir, nop1
        REAL(realk) :: umagsqr, sbu, sbv, sbw, pinf(1)
        REAL(realk) :: sb11, sb12, sb13, sb14, fak, ubfine
        REAL(realk) :: flag, vi, wi
        REAL(realk), POINTER, CONTIGUOUS :: u(:, :, :), v(:, :, :), &
            w(:, :, :), p(:, :, :), bp(:, :, :)
        REAL(realk), POINTER, CONTIGUOUS :: ubuf(:, :), vbuf(:, :), wbuf(:, :)
        REAL(realk), POINTER, CONTIGUOUS :: ddy(:), ddz(:)

        ! Return early when no action is to be taken
        SELECT CASE (ctyp)
        CASE ("FIX", "OP1", "NOS", "SLI", "PAR")
            CONTINUE
        CASE DEFAULT
            RETURN
        END SELECT

        ! Assure that required fields are present
        IF (.NOT. PRESENT(f2) .OR. .NOT. PRESENT(f3)) THEN
            CALL errr(__FILE__, __LINE__)
        END IF

        ! Fetch pointers
        CALL f1%get_ptr(u, igrid)
        CALL f2%get_ptr(v, igrid)
        CALL f3%get_ptr(w, igrid)

        SELECT CASE (ctyp)
        CASE ("FIX", "OP1", "PAR")
            CALL f1%get_buffer(ubuf, igrid, iface)
            CALL f2%get_buffer(vbuf, igrid, iface)
            CALL f3%get_buffer(wbuf, igrid, iface)
        END SELECT

        CALL get_fieldptr(bp, "BP", igrid)
        CALL get_fieldptr(ddy, "DDY", igrid)
        CALL get_fieldptr(ddz, "DDZ", igrid)
        CALL get_mgdims(kk, jj, ii, igrid)

        SELECT CASE (iface)
        CASE (1)
            ! Front
            i2 = 2
            i3 = 3
            i4 = 4
            istag1 = 1
            istag2 = 2
            dir = -1
        CASE (2)
            ! Back
            i2 = ii - 1
            i3 = ii - 2
            i4 = ii - 3
            istag1 = ii - 1
            istag2 = ii - 2
            dir = 1
        CASE DEFAULT
            CALL errr(__FILE__, __LINE__)
        END SELECT

        SELECT CASE (ctyp)
        CASE ("FIX")
            DO j = 1, jj
                DO k = 1, kk
                    u(k, j, istag2) = ubuf(k, j)
                    v(k, j, i2) = 2.0*vbuf(k, j) - v(k, j, i3)
                    w(k, j, i2) = 2.0*wbuf(k, j) - w(k, j, i3)
                END DO
            END DO
        CASE ("OP1")
            DO j = 1, jj-1
                DO k = 1, kk-1
                    ! Interface-normal velocity is always zero-gradient for OP1
                    u(k, j, istag1) = u(k, j, istag2)

                    ! Tangential velocities are zero-gradient for outflow,
                    ! and prescribed values in case of inflow
                    flag = SIGN(1.0_realk, &
                        dir*(u(k, j, istag2) + u(k, j+1, istag2)))
                    flag = 0.5*(flag + 1.0)
                    v(k, j, i2) = flag*v(k, j, i3) &
                        + (1.0-flag)*(2.0*vbuf(k, j) - v(k, j, i3))

                    flag = SIGN(1.0_realk, &
                        dir*(u(k, j, istag2) + u(k+1, j, istag2)))
                    flag = 0.5*(flag + 1.0)
                    w(k, j, i2) = flag*w(k, j, i3) &
                        + (1.0-flag)*(2.0*wbuf(k, j) - w(k, j, i3))
                END DO
            END DO
        CASE ("NOS")
            DO j = 1, jj
                DO k = 1, kk
                    u(k, j, istag2) = 0.0
                    v(k, j, i2) = 0.0
                    w(k, j, i2) = 0.0
                END DO
            END DO
        CASE ("SLI")
            DO j = 1, jj
                DO k = 1, kk
                    u(k, j, istag2) = 0.0
                    v(k, j, i2) = v(k, j, i3)
                    w(k, j, i2) = w(k, j, i3)
                END DO
            END DO
        CASE ("PAR")
            ! Tangential velocity components
            DO j = 2, jj-1
                DO k = 2, kk-1
                    sbv = bp(k, j, i2)*bp(k, j+1, i2)
                    sbw = bp(k, j, i2)*bp(k+1, j, i2)
                    v(k, j, i2) = (vbuf(k, j) &
                        - 0.5*(v(k, j, i3) - v(k, j, i4)))*sbv &
                        + (1.0-sbv)*v(k, j, i2)
                    w(k, j, i2) = (wbuf(k, j) &
                        - 0.5*(w(k, j, i3)-w(k, j, i4)))*sbw &
                        + (1.0-sbw)*w(k, j, i2)
                END DO
            END DO

            ! Wall-normal velocity component
            DO j = 3, jj-2, 2
                DO k = 3, kk-2, 2
                    sbu = bp(k, j, i2)

                    IF (sbu > 0.5) THEN
                        sb11 = bp(k, j, i2)*bp(k, j, i3)
                        sb12 = bp(k, j+1, i2)*bp(k, j+1, i3)
                        sb13 = bp(k+1, j, i2)*bp(k+1, j, i3)
                        sb14 = bp(k+1, j+1, i2)*bp(k+1, j+1, i3)
                    ELSE
                        sb11 = MIN(bp(k, j, i3) &
                            + bp(k, j, i2) &
                            + bp(k, j, i4) &
                            + bp(k, j-1, i3) &
                            + bp(k, j+1, i3) &
                            + bp(k-1, j, i3) &
                            + bp(k+1, j, i3), 1.0_realk)
                        sb12 = MIN(bp(k, j+1, i3) &
                            + bp(k, j+1, i2) &
                            + bp(k, j+1, i4) &
                            + bp(k, j, i3) &
                            + bp(k, j+2, i3) &
                            + bp(k-1, j+1, i3) &
                            + bp(k+1, j+1, i3), 1.0_realk)
                        sb13 = MIN(bp(k+1, j, i3) &
                            + bp(k+1, j, i2) &
                            + bp(k+1, j, i4) &
                            + bp(k+1, j-1, i3) &
                            + bp(k+1, j+1, i3) &
                            + bp(k, j, i3) &
                            + bp(k+2, j, i3), 1.0_realk)
                        sb14 = MIN(bp(k+1, j+1, i3) &
                            + bp(k+1, j+1, i2) &
                            + bp(k+1, j+1, i4) &
                            + bp(k+1, j, i3) &
                            + bp(k+1, j+2, i3) &
                            + bp(k, j+1, i3) &
                            + bp(k+2, j+1, i3), 1.0_realk)
                    END IF

                    fak = 1.0/((ddy(j)+ddy(j+1))*(ddz(k)+ddz(k+1))) &
                        * (sb11*ddy(j)*ddz(k) &
                        + sb12*ddy(j+1)*ddz(k) &
                        + sb13*ddy(j)*ddz(k+1) &
                        + sb14*ddy(j+1)*ddz(k+1))
                    IF (fak < 0.1) fak = 1.0
                    fak = 1.0/fak

                    ubfine = (u(k, j, istag2)*(1.0-sb11)*ddy(j)*ddz(k) &
                        + u(k, j+1, istag2)*(1.0-sb12)*ddy(j+1)*ddz(k) &
                        + u(k+1, j, istag2)*(1.0-sb13)*ddy(j)*ddz(k+1) &
                        + u(k+1, j+1, istag2)*(1.0-sb14)*ddy(j+1)*ddz(k+1)) &
                        /((ddy(j)+ddy(j+1))*(ddz(k)+ddz(k+1)))

                    u(k, j, istag2) = (ubuf(k, j)-ubfine)*fak*sb11 &
                        + (1.0-sb11)*u(k, j, istag2)*sbu
                    u(k, j+1, istag2) = (ubuf(k, j+1)-ubfine)*fak*sb12 &
                        + (1.0-sb12)*u(k, j+1, istag2)*sbu
                    u(k+1, j, istag2) = (ubuf(k+1, j)-ubfine)*fak*sb13 &
                        + (1.0-sb13)*u(k+1, j, istag2)*sbu
                    u(k+1, j+1, istag2) = (ubuf(k+1, j+1)-ubfine)*fak*sb14 &
                        + (1.0-sb14)*u(k+1, j+1, istag2)*sbu
                END DO
            END DO
        END SELECT

        ! Pressure is optional, when called from setpointvalues for instance,
        ! we do not want to change the pressure
        IF (PRESENT(f4)) THEN
            CALL f4%get_ptr(p, igrid)
            SELECT CASE (ctyp)
            CASE ("OP1")
                CALL get_nbcprms(igrid, iface, ibocd, nreal=nop1)
                IF (nop1 == 0) THEN
                    pinf(1) = 0.0
                ELSE IF (nop1 == 1) THEN
                    CALL get_bcprms(pinf, igrid, iface, ibocd)
                ELSE
                    CALL errr(__FILE__, __LINE__)
                END IF

                DO j = 2, jj
                    DO k = 2, kk
                        vi = 0.25*(v(k, j-1, i2) + v(k, j, i2) &
                            + v(k, j-1, i3) + v(k, j, i3))
                        wi = 0.25*(w(k-1, j, i2) + w(k, j, i2) &
                            + w(k-1, j, i3) + w(k, j, i3))
                        umagsqr = u(k, j, istag2)**2 + vi**2 + wi**2

                        p(k, j, i2) = pinf(1) + MIN(0.0_realk, &
                            0.5*rho*SIGN(umagsqr, dir*u(k, j, istag2)))
                    END DO
                END DO
            CASE ("FIX", "NOS", "SLI")
                DO j = 1, jj
                    DO k = 1, kk
                        p(k, j, i2) = p(k, j, i3)
                    END DO
                END DO
            END SELECT
        END IF
    END SUBROUTINE bfront


    SUBROUTINE bright(igrid, iface, ibocd, ctyp, f1, f2, f3, f4, timeph)
        ! Subroutine arguments
        INTEGER(intk), INTENT(in) :: igrid, iface, ibocd
        CHARACTER(len=*), INTENT(in) :: ctyp
        TYPE(field_t), INTENT(inout) :: f1
        TYPE(field_t), INTENT(inout), OPTIONAL :: f2, f3, f4
        REAL(realk), INTENT(in), OPTIONAL :: timeph

        ! Local variables
        INTEGER(intk) :: kk, jj, ii
        INTEGER(intk) :: k, i, j2, j3, j4, jstag1, jstag2, dir, nop1
        REAL(realk) :: umagsqr, sbu, sbv, sbw, pinf(1)
        REAL(realk) :: sb11, sb12, sb13, sb14, fak, vbfine
        REAL(realk) :: flag, ui, wi
        REAL(realk), POINTER, CONTIGUOUS :: u(:, :, :), v(:, :, :), &
            w(:, :, :), p(:, :, :), bp(:, :, :)
        REAL(realk), POINTER, CONTIGUOUS :: ubuf(:, :), vbuf(:, :), wbuf(:, :)
        REAL(realk), POINTER, CONTIGUOUS :: ddx(:), ddz(:)

        ! Return early when no action is to be taken
        SELECT CASE (ctyp)
        CASE ("FIX", "OP1", "NOS", "SLI", "PAR")
            CONTINUE
        CASE DEFAULT
            RETURN
        END SELECT

        ! Assure that required fields are present
        IF (.NOT. PRESENT(f2) .OR. .NOT. PRESENT(f3)) THEN
            CALL errr(__FILE__, __LINE__)
        END IF

        ! Fetch pointers
        CALL f1%get_ptr(u, igrid)
        CALL f2%get_ptr(v, igrid)
        CALL f3%get_ptr(w, igrid)

        SELECT CASE (ctyp)
        CASE ("FIX", "OP1", "PAR")
            CALL f1%get_buffer(ubuf, igrid, iface)
            CALL f2%get_buffer(vbuf, igrid, iface)
            CALL f3%get_buffer(wbuf, igrid, iface)
        END SELECT

        CALL get_fieldptr(bp, "BP", igrid)
        CALL get_fieldptr(ddx, "DDX", igrid)
        CALL get_fieldptr(ddz, "DDZ", igrid)
        CALL get_mgdims(kk, jj, ii, igrid)

        SELECT CASE (iface)
        CASE (3)
            ! Right
            j2 = 2
            j3 = 3
            j4 = 4
            jstag1 = 1
            jstag2 = 2
            dir = -1
        CASE (4)
            ! Left
            j2 = jj - 1
            j3 = jj - 2
            j4 = jj - 3
            jstag1 = jj - 1
            jstag2 = jj - 2
            dir = 1
        CASE DEFAULT
            CALL errr(__FILE__, __LINE__)
        END SELECT

        SELECT CASE (ctyp)
        CASE ("FIX")
            DO i = 1, ii
                DO k = 1, kk
                    u(k, j2, i) = 2.0*ubuf(k, i) - u(k, j3, i)
                    v(k, jstag2, i) = vbuf(k, i)
                    w(k, j2, i) = 2.0*wbuf(k, i) - w(k, j3, i)
                END DO
            END DO
        CASE ("OP1")
            DO i = 1, ii-1
                DO k = 1, kk-1
                    ! Tangential velocities are zero-gradient for outflow,
                    ! and prescribed values in case of inflow
                    flag = SIGN(1.0_realk, &
                        dir*(v(k, jstag2, i) + v(k, jstag2, i+1)))
                    flag = 0.5*(flag + 1.0)
                    u(k, j2, i) = flag*u(k, j3, i) &
                        + (1.0-flag)*(2.0*ubuf(k, i) - u(k, j3, i))

                    ! Interface-normal velocity is always zero-gradient for OP1
                    v(k, jstag1, i) = v(k, jstag2, i)

                    ! Tangential velocities are zero-gradient for outflow,
                    ! and prescribed values in case of inflow
                    flag = SIGN(1.0_realk, &
                        dir*(v(k, jstag2, i) + v(k+1, jstag2, i)))
                    flag = 0.5*(flag + 1.0)
                    w(k, j2, i) = flag*w(k, j3, i) &
                        + (1.0-flag)*(2.0*wbuf(k, i) - w(k, j3, i))
                END DO
            END DO
        CASE ("NOS")
            DO i = 1, ii
                DO k = 1, kk
                    u(k, j2, i) = 0.0
                    v(k, jstag2, i) = 0.0
                    w(k, j2, i) = 0.0
                END DO
            END DO
        CASE ("SLI")
            DO i = 1, ii
                DO k = 1, kk
                    u(k, j2, i) = u(k, j3, i)
                    v(k, jstag2, i) = 0.0
                    w(k, j2, i) = w(k, j3, i)
                END DO
            END DO
        CASE ("PAR")
            ! Tangential velocity components
            DO i = 2, ii-1
                DO k = 2, kk-1
                    sbu = bp(k, j2, i)*bp(k, j2, i+1)
                    sbw = bp(k, j2, i)*bp(k+1, j2, i)
                    u(k, j2, i) = (ubuf(k, i) &
                        - 0.5*(u(k, j3, i) - u(k, j4, i)))*sbu &
                        + (1.0-sbu)*u(k, j2, i)
                    w(k, j2, i) = (wbuf(k, i) &
                        - 0.5*(w(k, j2, i)-w(k, j4, i)))*sbw &
                        + (1.0-sbw)*w(k, j2, i)
                END DO
            END DO

            ! Wall-normal velocity component
            DO i = 3, ii-2, 2
                DO k = 3, kk-2, 2
                    sbv = bp(k, j2, i)

                    IF (sbv > 0.5) THEN
                        sb11 = bp(k, j2, i)*bp(k, j3, i)
                        sb12 = bp(k, j2, i+1)*bp(k, j3, i+1)
                        sb13 = bp(k+1, j2, i)*bp(k+1, j3, i)
                        sb14 = bp(k+1, j2, i+1)*bp(k+1, j3, i+1)
                    ELSE
                        sb11 = MIN(bp(k, j3, i) &
                            + bp(k, j2, i) &
                            + bp(k, j4, i) &
                            + bp(k, j3, i-1) &
                            + bp(k, j3, i+1) &
                            + bp(k-1, j3, i) &
                            + bp(k+1, j3, i), 1.0_realk)
                        sb12 = MIN(bp(k, j3, i+1) &
                            + bp(k, j2, i+1) &
                            + bp(k, j4, i+1) &
                            + bp(k, j3, i) &
                            + bp(k, j3, i+2) &
                            + bp(k-1, j3, i+1) &
                            + bp(k+1, j3, i+1), 1.0_realk)
                        sb13 = MIN(bp(k+1, j3, i) &
                            + bp(k+1, j2, i) &
                            + bp(k+1, j4, i) &
                            + bp(k+1, j3, i-1) &
                            + bp(k+1, j3, i+1) &
                            + bp(k, j3, i) &
                            + bp(k+2, j3, i), 1.0_realk)
                        sb14 = MIN(bp(k+1, j3, i+1) &
                            + bp(k+1, j2, i+1) &
                            + bp(k+1, j4, i+1) &
                            + bp(k+1, j3, i) &
                            + bp(k+1, j3, i+2) &
                            + bp(k, j3, i+1) &
                            + bp(k+2, j3, i+1), 1.0_realk)
                    END IF

                    fak = 1.0/((ddz(k)+ddz(k+1))*(ddx(i)+ddx(i+1))) &
                        *(sb11*ddz(k)*ddx(i) &
                        + sb12*ddz(k)*ddx(i+1) &
                        + sb13*ddz(k+1)*ddx(i) &
                        + sb14*ddz(k+1)*ddx(i+1))
                    IF (fak < 0.1) fak = 1.0
                    fak = 1.0/fak

                    vbfine = (v(k, jstag2, i)*(1.0-sb11)*ddz(k)*ddx(i) &
                        + v(k, jstag2, i+1)*(1.0-sb12)*ddz(k)*ddx(i+1) &
                        + v(k+1, jstag2, i)*(1.0-sb13)*ddz(k+1)*ddx(i) &
                        + v(k+1, jstag2, i+1)*(1.0-sb14)*ddz(k+1)*ddx(i+1)) &
                        /((ddz(k)+ddz(k+1))*(ddx(i)+ddx(i+1)))

                    v(k, jstag2, i) = (vbuf(k, i)-vbfine)*fak*sb11 &
                        + v(k, jstag2, i)*(1.-sb11)*sbv
                    v(k, jstag2, i+1) = (vbuf(k, i+1)-vbfine)*fak*sb12 &
                        + v(k, jstag2, i+1)*(1.-sb12)*sbv
                    v(k+1, jstag2, i) = (vbuf(k+1, i)-vbfine)*fak*sb13 &
                        + v(k+1, jstag2, i)*(1.-sb13)*sbv
                    v(k+1, jstag2, i+1) = (vbuf(k+1, i+1)-vbfine)*fak*sb14 &
                        + v(k+1, jstag2, i+1)*(1.-sb14)*sbv
                END DO
            END DO
        END SELECT

        ! Pressure is optional, when called from setpointvalues for instance,
        ! we do not want to change the pressure
        IF (PRESENT(f4)) THEN
            CALL f4%get_ptr(p, igrid)
            SELECT CASE (ctyp)
            CASE ("OP1")
                CALL get_nbcprms(igrid, iface, ibocd, nreal=nop1)
                IF (nop1 == 0) THEN
                    pinf(1) = 0.0
                ELSE IF (nop1 == 1) THEN
                    CALL get_bcprms(pinf, igrid, iface, ibocd)
                ELSE
                    CALL errr(__FILE__, __LINE__)
                END IF

                DO i = 2, ii
                    DO k = 2, kk
                        ui = 0.25*(u(k, j2, i-1) + u(k, j2, i) &
                            + u(k, j3, i-1) + u(k, j3, i))
                        wi = 0.25*(w(k-1, j2, i) + w(k, j2, i) &
                            + w(k-1, j3, i) + w(k, j3, i))
                        umagsqr = ui**2 + v(k, jstag2, i)**2 + wi**2

                        p(k, j2, i) = pinf(1) + MIN(0.0_realk, &
                            0.5*rho*SIGN(umagsqr, dir*v(k, jstag2, i)))
                    END DO
                END DO
            CASE ("FIX", "NOS", "SLI")
                DO i = 1, ii
                    DO k = 1, kk
                        p(k, j2, i) = p(k, j3, i)
                    END DO
                END DO
            END SELECT
        END IF
    END SUBROUTINE bright


    SUBROUTINE bbottom(igrid, iface, ibocd, ctyp, f1, f2, f3, f4, timeph)
        ! Subroutine arguments
        INTEGER(intk), INTENT(in) :: igrid, iface, ibocd
        CHARACTER(len=*), INTENT(in) :: ctyp
        TYPE(field_t), INTENT(inout) :: f1
        TYPE(field_t), INTENT(inout), OPTIONAL :: f2, f3, f4
        REAL(realk), INTENT(in), OPTIONAL :: timeph

        ! Local variables
        INTEGER(intk) :: kk, jj, ii
        INTEGER(intk) :: j, i, k2, k3, k4, kstag1, kstag2, dir, nop1
        REAL(realk) :: umagsqr, sbu, sbv, sbw, pinf(1)
        REAL(realk) :: sb11, sb12, sb13, sb14, fak, wbfine
        REAL(realk) :: flag, ui, vi
        REAL(realk), POINTER, CONTIGUOUS :: u(:, :, :), v(:, :, :), &
            w(:, :, :), p(:, :, :), bp(:, :, :)
        REAL(realk), POINTER, CONTIGUOUS :: ubuf(:, :), vbuf(:, :), wbuf(:, :)
        REAL(realk), POINTER, CONTIGUOUS :: ddx(:), ddy(:)

        ! Return early when no action is to be taken
        SELECT CASE (ctyp)
        CASE ("FIX", "OP1", "NOS", "SLI", "PAR")
            CONTINUE
        CASE DEFAULT
            RETURN
        END SELECT

        ! Assure that required fields are present
        IF (.NOT. PRESENT(f2) .OR. .NOT. PRESENT(f3)) THEN
            CALL errr(__FILE__, __LINE__)
        END IF

        ! Fetch pointers
        CALL f1%get_ptr(u, igrid)
        CALL f2%get_ptr(v, igrid)
        CALL f3%get_ptr(w, igrid)

        SELECT CASE (ctyp)
        CASE ("FIX", "OP1", "PAR")
            CALL f1%get_buffer(ubuf, igrid, iface)
            CALL f2%get_buffer(vbuf, igrid, iface)
            CALL f3%get_buffer(wbuf, igrid, iface)
        END SELECT

        CALL get_fieldptr(bp, "BP", igrid)
        CALL get_fieldptr(ddx, "DDX", igrid)
        CALL get_fieldptr(ddy, "DDY", igrid)
        CALL get_mgdims(kk, jj, ii, igrid)

        SELECT CASE (iface)
        CASE (5)
            ! Bottom
            k2 = 2
            k3 = 3
            k4 = 4
            kstag1 = 1
            kstag2 = 2
            dir = -1
        CASE (6)
            ! Top
            k2 = kk - 1
            k3 = kk - 2
            k4 = kk - 3
            kstag1 = kk - 1
            kstag2 = kk - 2
            dir = 1
        CASE DEFAULT
            CALL errr(__FILE__, __LINE__)
        END SELECT

        SELECT CASE (ctyp)
        CASE ("FIX")
            DO i = 1, ii
                DO j = 1, jj
                    u(k2, j, i) = 2.0*ubuf(j, i) - u(k3, j, i)
                    v(k2, j, i) = 2.0*vbuf(j, i) - v(k3, j, i)
                    w(kstag2, j, i) = wbuf(j, i)
                END DO
            END DO
        CASE ("OP1")
            DO i = 1, ii-1
                DO j = 1, jj-1
                    ! Tangential velocities are zero-gradient for outflow,
                    ! and prescribed values in case of inflow
                    flag = SIGN(1.0_realk, &
                        dir*(w(kstag2, j, i) + w(kstag2, j, i+1)))
                    flag = 0.5*(flag + 1.0)
                    u(k2, j, i) = flag*u(k3, j, i) &
                        + (1.0-flag)*(2*ubuf(j, i) - u(k3, j, i))

                    flag = SIGN(1.0_realk, &
                        dir*(w(kstag2, j, i) + w(kstag2, j+1, i)))
                    flag = 0.5*(flag + 1.0)
                    v(k2, j, i) = flag*v(k3, j, i) &
                        + (1.0-flag)*(2.0*vbuf(j, i) - v(k3, j, i))

                    ! Interface-normal velocity is always zero-gradient for OP1
                    w(kstag1, j, i) = w(kstag2, j, i)
                END DO
            END DO
        CASE ("NOS")
            DO i = 1, ii
                DO j = 1, jj
                    u(k2, j, i) = 0.0
                    v(k2, j, i) = 0.0
                    w(kstag2, j, i) = 0.0
                END DO
            END DO
        CASE ("SLI")
            DO i = 1, ii
                DO j = 1, jj
                    u(k2, j, i) = u(k3, j, i)
                    v(k2, j, i) = v(k3, j, i)
                    w(kstag2, j, i) = 0.0
                END DO
            END DO
        CASE ("PAR")
            ! Tangential velocity components
            DO i = 2, ii-1
                DO j = 2, jj-1
                    sbu = bp(k2, j, i)*bp(k2, j, i+1)
                    sbv = bp(k2, j, i)*bp(k2, j+1, i)
                    u(k2, j, i) = (ubuf(j, i) &
                        - 0.5*(u(k3, j, i) - u(k4, j, i)))*sbu &
                        + (1.0-sbu)*u(k2, j, i)
                    v(k2, j, i) = (vbuf(j, i) &
                        - 0.5*(v(k3, j, i) - v(k4, j, i)))*sbv &
                        + (1.0-sbv)*v(k2, j, i)
                END DO
            END DO

            ! Wall-normal velocity component
            DO i = 3, ii-2, 2
                DO j = 3, jj-2, 2
                    sbw = bp(k2, j, i)

                    IF (sbw > 0.5) THEN
                        sb11 = bp(k2, j, i)*bp(k3, j, i)
                        sb12 = bp(k2, j, i+1)*bp(k3, j, i+1)
                        sb13 = bp(k2, j+1, i)*bp(k3, j+1, i)
                        sb14 = bp(k2, j+1, i+1)*bp(k3, j+1, i+1)
                    ELSE
                        sb11 = MIN(bp(k3, j, i) &
                            + bp(k2, j, i) &
                            + bp(k4, j, i) &
                            + bp(k3, j, i-1) &
                            + bp(k3, j, i+1) &
                            + bp(k3, j-1, i) &
                            + bp(k3, j+1, i), 1.0_realk)
                        sb12 = MIN(bp(k3, j, i+1) &
                            + bp(k2, j, i+1) &
                            + bp(k4, j, i+1) &
                            + bp(k3, j, i) &
                            + bp(k3, j, i+2) &
                            + bp(k3, j-1, i+1) &
                            + bp(k3, j+1, i+1), 1.0_realk)
                        sb13 = MIN(bp(k3, j+1, i) &
                            + bp(k2, j+1, i) &
                            + bp(k4, j+1, i) &
                            + bp(k3, j+1, i-1) &
                            + bp(k3, j+1, i+1) &
                            + bp(k3, j, i) &
                            + bp(k3, j+2, i), 1.0_realk)
                        sb14 = MIN(bp(k3, j+1, i+1) &
                            + bp(k2, j+1, i+1) &
                            + bp(k4, j+1, i+1) &
                            + bp(k3, j+1, i) &
                            + bp(k3, j+1, i+2) &
                            + bp(k3, j, i+1) &
                            + bp(k3, j+2, i+1), 1.0_realk)
                    END IF

                    fak = 1./((ddy(j)+ddy(j+1))*(ddx(i)+ddx(i+1))) &
                        * (sb11*ddy(j)*ddx(i) &
                        + sb12*ddy(j)*ddx(i+1) &
                        + sb13*ddy(j+1)*ddx(i) &
                        + sb14*ddy(j+1)*ddx(i+1))
                    IF (fak < 0.1) fak = 1.0
                    fak = 1.0/fak

                    wbfine = (w(kstag2, j, i)*(1.0-sb11)*ddy(j)*ddx(i) &
                        + w(kstag2, j, i+1)*(1.0-sb12)*ddy(j)*ddx(i+1) &
                        + w(kstag2, j+1, i)*(1.0-sb13)*ddy(j+1)*ddx(i) &
                        + w(kstag2, j+1, i+1)*(1.0-sb14)*ddy(j+1)*ddx(i+1)) &
                        /((ddy(j)+ddy(j+1))*(ddx(i)+ddx(i+1)))

                    w(kstag2, j, i) = (wbuf(j, i)-wbfine)*fak*sb11 &
                        + w(kstag2, j, i)*(1.0-sb11)*sbw
                    w(kstag2, j, i+1) = (wbuf(j, i+1)-wbfine)*fak*sb12 &
                        + w(kstag2, j, i+1)*(1.0-sb12)*sbw
                    w(kstag2, j+1, i) = (wbuf(j+1, i)-wbfine)*fak*sb13 &
                        + w(kstag2, j+1, i)*(1.0-sb13)*sbw
                    w(kstag2, j+1, i+1) = (wbuf(j+1, i+1)-wbfine)*fak*sb14 &
                        + w(kstag2, j+1, i+1)*(1.0-sb14)*sbw
                END DO
            END DO
        END SELECT

        ! Pressure is optional, when called from setpointvalues for instance,
        ! we do not want to change the pressure
        IF (PRESENT(f4)) THEN
            CALL f4%get_ptr(p, igrid)
            SELECT CASE (ctyp)
            CASE ("OP1")
                CALL get_nbcprms(igrid, iface, ibocd, nreal=nop1)
                IF (nop1 == 0) THEN
                    pinf(1) = 0.0
                ELSE IF (nop1 == 1) THEN
                    CALL get_bcprms(pinf, igrid, iface, ibocd)
                ELSE
                    CALL errr(__FILE__, __LINE__)
                END IF

                DO i = 2, ii
                    DO j = 2, jj
                        ui = 0.25*(u(k2, j, i-1) + u(k2, j, i) &
                            + u(k3, j, i-1) + u(k3, j, i))
                        vi = 0.25*(v(k2, j-1, i) + v(k2, j, i) &
                            + v(k3, j-1, i) + v(k3, j, i))

                        umagsqr = ui**2 + vi**2 + w(kstag2, j, i)**2
                        p(k2, j, i) = pinf(1) + MIN(0.0_realk, &
                            0.5*rho*SIGN(umagsqr, dir*w(kstag2, j, i)))
                    END DO
                END DO
            CASE ("FIX", "NOS", "SLI")
                DO i = 1, ii
                    DO j = 1, jj
                        p(k2, j, i) = p(k3, j, i)
                    END DO
                END DO
            END SELECT
        END IF
    END SUBROUTINE bbottom
END MODULE bound_flow_mod
