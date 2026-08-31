MODULE bound_scalar_mod
    USE core_mod
    USE scacore_mod
    USE flow_mod, ONLY: ilesmodel, gmol, rho, qwallfix

    IMPLICIT NONE(type, external)
    PRIVATE

    ! Bound operation 'T' operate on U, V, W, P
    TYPE, EXTENDS(bound_t) :: bound_scaflux_t
    CONTAINS
        PROCEDURE, NOPASS :: front => bfront
        PROCEDURE, NOPASS :: back => bfront
        PROCEDURE, NOPASS :: right => bright
        PROCEDURE, NOPASS :: left => bright
        PROCEDURE, NOPASS :: bottom => bbottom
        PROCEDURE, NOPASS :: top => bbottom
    END TYPE bound_scaflux_t
    TYPE(bound_scaflux_t) :: bound_scaflux

    INTEGER(intk), PARAMETER :: boundscatasksize = 4
    INTEGER(intk), ALLOCATABLE :: boundscatasks(:, :, :, :)
    INTEGER(intk), ALLOCATABLE :: nboundscataskslvl(:, :)
    !$omp declare target(boundscatasks)

    PUBLIC :: bound_scaflux, init_bound_scalar, apply_bound_scalar, &
        finish_bound_scalar

CONTAINS
    SUBROUTINE init_bound_scalar()
        INTEGER(intk) :: nlevels, maxntasks
        INTEGER(intk) :: l, ilevel, ilevel_index, imygrid, igrid
        INTEGER(intk) :: iface, ibocd, nbocd, itask, ityp
        INTEGER(intk) :: scbtype(nsca)
        CHARACTER(len=8) :: ctyp

        nlevels = maxlevel - minlevel + 1
        ALLOCATE(nboundscataskslvl(nlevels, nsca), source=0_intk)

        DO l = 1, nsca
            DO ilevel = minlevel, maxlevel
                ilevel_index = ilevel - minlevel + 1
                DO imygrid = 1, nmygridslvl(ilevel)
                    igrid = mygridslvl(imygrid, ilevel)
                    DO iface = 1, 6
                        nbocd = nboconds(iface, igrid)
                        DO ibocd = 1, nbocd
                            CALL get_bc_ctyp(ctyp, ibocd, iface, igrid)
                            IF (boundsca_type(ctyp) == 0) CYCLE
                            nboundscataskslvl(ilevel_index, l) = &
                                nboundscataskslvl(ilevel_index, l) + 1
                        END DO
                    END DO
                END DO
            END DO
        END DO

        maxntasks = MAX(1_intk, MAXVAL(nboundscataskslvl))
        ALLOCATE(boundscatasks(boundscatasksize, maxntasks, nlevels, nsca), &
            source=-1_intk)

        DO l = 1, nsca
            DO ilevel = minlevel, maxlevel
                ilevel_index = ilevel - minlevel + 1
                itask = 0
                DO imygrid = 1, nmygridslvl(ilevel)
                    igrid = mygridslvl(imygrid, ilevel)
                    DO iface = 1, 6
                        nbocd = nboconds(iface, igrid)
                        DO ibocd = 1, nbocd
                            CALL get_bc_ctyp(ctyp, ibocd, iface, igrid)
                            ityp = boundsca_type(ctyp)
                            IF (ityp == 0) CYCLE
                            itask = itask + 1
                            boundscatasks(1:3, itask, ilevel_index, l) = &
                                [igrid, iface, ityp]
                            IF (ityp == 1) THEN
                                boundscatasks(4, itask, ilevel_index, l) = -1
                            ELSE
                                CALL get_bcprms(scbtype, igrid, iface, ibocd)
                                boundscatasks(4, itask, ilevel_index, l) = &
                                    scbtype(l)
                            END IF
                        END DO
                    END DO
                END DO
            END DO
        END DO

        !$omp target enter data map(always, to: boundscatasks)
    END SUBROUTINE init_bound_scalar


    INTEGER(intk) FUNCTION boundsca_type(ctyp)
        CHARACTER(len=*), INTENT(in) :: ctyp

        SELECT CASE (ctyp)
        CASE ("PAR")
            boundsca_type = 1
        CASE ("SWA")
            boundsca_type = 2
        CASE ("SIO")
            boundsca_type = 3
        CASE DEFAULT
            boundsca_type = 0
        END SELECT
    END FUNCTION boundsca_type


    SUBROUTINE apply_bound_scalar(ilevel, isca, qtu_f, qtv_f, qtw_f, t_f)
        INTEGER(intk), INTENT(in) :: ilevel, isca
        TYPE(field_t), INTENT(inout) :: qtu_f, qtv_f, qtw_f, t_f

        INTEGER(intk) :: ilevel_index, ntasks
        TYPE(field_t), POINTER :: u_f, v_f, w_f, bt_f
        TYPE(field_t), POINTER :: dx_f, dy_f, dz_f, ddx_f, ddy_f, ddz_f

        ilevel_index = ilevel - minlevel + 1
        ntasks = nboundscataskslvl(ilevel_index, isca)
        IF (ntasks == 0) RETURN

        CALL get_field(u_f, "U")
        CALL get_field(v_f, "V")
        CALL get_field(w_f, "W")
        CALL get_field(bt_f, "BT")
        CALL get_field(dx_f, "DX")
        CALL get_field(dy_f, "DY")
        CALL get_field(dz_f, "DZ")
        CALL get_field(ddx_f, "DDX")
        CALL get_field(ddy_f, "DDY")
        CALL get_field(ddz_f, "DDZ")

        CALL apply_bound_scalar_impl(ilevel_index, isca, ntasks, &
            qtu_f%arr, qtv_f%arr, qtw_f%arr, t_f%arr, qtu_f%buffers, &
            qtv_f%buffers, qtw_f%buffers, t_f%buffers, u_f%arr, v_f%arr, &
            w_f%arr, bt_f%arr, dx_f%arr, dy_f%arr, dz_f%arr, ddx_f%arr, &
            ddy_f%arr, ddz_f%arr, scalar(isca)%prmol)
    END SUBROUTINE apply_bound_scalar


    SUBROUTINE apply_bound_scalar_impl(ilevel_index, isca, ntasks, qtu, qtv, &
            qtw, t, qtubuf, qtvbuf, qtwbuf, tbuf, u, v, w, bt, dx, dy, dz, &
            ddx, ddy, ddz, prmol)
        INTEGER(intk), INTENT(in) :: ilevel_index, isca, ntasks
        REAL(realk), INTENT(inout) :: qtu(*), qtv(*), qtw(*), t(*)
        REAL(realk), INTENT(inout) :: qtubuf(*), qtvbuf(*), qtwbuf(*), tbuf(*)
        REAL(realk), INTENT(in) :: u(*), v(*), w(*), bt(*)
        REAL(realk), INTENT(in) :: dx(*), dy(*), dz(*)
        REAL(realk), INTENT(in) :: ddx(*), ddy(*), ddz(*), prmol

        INTEGER(intk) :: itask, igrid, iface, ityp, bctype
        INTEGER(intk) :: kk, jj, ii, ip3, ipbb, ipx, ipy, ipz

        !$omp target teams distribute private(itask, igrid, iface, ityp, &
        !$omp& bctype, kk, jj, ii, ip3, ipbb, ipx, ipy, ipz)
        DO itask = 1, ntasks
            igrid = boundscatasks(1, itask, ilevel_index, isca)
            iface = boundscatasks(2, itask, ilevel_index, isca)
            ityp = boundscatasks(3, itask, ilevel_index, isca)
            bctype = boundscatasks(4, itask, ilevel_index, isca)
            CALL get_mgdims(kk, jj, ii, igrid)
            CALL get_ip3(ip3, igrid)
            CALL get_ipbb(ipbb, iface, igrid)
            CALL get_ip1x(ipx, igrid)
            CALL get_ip1y(ipy, igrid)
            CALL get_ip1z(ipz, igrid)

            !$omp parallel
            CALL bound_scalar_device(kk, jj, ii, iface, ityp, bctype, &
                prmol, qtu(ip3), qtv(ip3), qtw(ip3), t(ip3), u(ip3), &
                v(ip3), w(ip3), bt(ip3), qtubuf(ipbb), qtvbuf(ipbb), &
                qtwbuf(ipbb), tbuf(ipbb), dx(ipx), dy(ipy), dz(ipz), &
                ddx(ipx), ddy(ipy), ddz(ipz))
            !$omp end parallel
        END DO
        !$omp end target teams distribute
    END SUBROUTINE apply_bound_scalar_impl


    SUBROUTINE finish_bound_scalar()
        !$omp target exit data map(always, delete: boundscatasks)
        DEALLOCATE(boundscatasks)
        DEALLOCATE(nboundscataskslvl)
    END SUBROUTINE finish_bound_scalar


    SUBROUTINE bound_scalar_device(kk, jj, ii, iface, ityp, bctype, prmol, &
            qtu, qtv, qtw, t, u, v, w, bt, qtubuf, qtvbuf, qtwbuf, tbuf, &
            dx, dy, dz, ddx, ddy, ddz)
        !$omp declare target
        INTEGER(intk), INTENT(in) :: kk, jj, ii, iface, ityp, bctype
        REAL(realk), INTENT(in) :: prmol
        REAL(realk), INTENT(inout) :: qtu(kk, jj, ii)
        REAL(realk), INTENT(inout) :: qtv(kk, jj, ii)
        REAL(realk), INTENT(inout) :: qtw(kk, jj, ii)
        REAL(realk), INTENT(in) :: t(kk, jj, ii)
        REAL(realk), INTENT(in) :: u(kk, jj, ii), v(kk, jj, ii)
        REAL(realk), INTENT(in) :: w(kk, jj, ii), bt(kk, jj, ii)
        REAL(realk), INTENT(in) :: qtubuf(*), qtvbuf(*), qtwbuf(*), tbuf(*)
        REAL(realk), INTENT(in) :: dx(ii), dy(jj), dz(kk)
        REAL(realk), INTENT(in) :: ddx(ii), ddy(jj), ddz(kk)

        INTEGER(intk) :: i, j, k, i3, j3, k3, istag2, jstag2, kstag2
        INTEGER(intk) :: dir, ibuf
        REAL(realk) :: area1, area2, area3, area4, arecvtot, qtot
        REAL(realk) :: area, adv, diff, gamma, gamma2dx, uquer, tout

        SELECT CASE (iface)
        CASE (1, 2)
            IF (iface == 1) THEN
                i3 = 3
                istag2 = 2
                dir = -1
            ELSE
                i3 = ii-2
                istag2 = ii-2
                dir = 1
            END IF

            SELECT CASE (ityp)
            CASE (1)
                !$omp do collapse(2) private(j, k, ibuf, area1, area2, &
                !$omp& area3, area4, arecvtot, qtot)
                DO j = 3, jj-2, 2
                    DO k = 3, kk-2, 2
                        area1 = bt(k, j, i3)
                        area2 = bt(k+1, j, i3)
                        area3 = bt(k, j+1, i3)
                        area4 = bt(k+1, j+1, i3)
                        arecvtot = area1 + area2 + area3 + area4
                        ibuf = k + (j-1)*kk
                        qtot = qtubuf(ibuf)
                        qtu(k, j, istag2) = divide0(area1, arecvtot)*qtot
                        qtu(k+1, j, istag2) = &
                            divide0(area2, arecvtot)*qtot
                        qtu(k, j+1, istag2) = &
                            divide0(area3, arecvtot)*qtot
                        qtu(k+1, j+1, istag2) = &
                            divide0(area4, arecvtot)*qtot
                    END DO
                END DO
                !$omp end do
            CASE (2)
                IF (bctype == 0 .AND. ilesmodel == 0) THEN
                    gamma2dx = 2.0*gmol/rho/prmol/dx(istag2)
                    !$omp do collapse(2) private(j, k, ibuf, diff)
                    DO j = 1, jj
                        DO k = 1, kk
                            ibuf = k + (j-1)*kk
                            diff = gamma2dx*(tbuf(ibuf) - t(k, j, i3))
                            qtu(k, j, istag2) = &
                                -dir*diff*ddy(j)*ddz(k)
                        END DO
                    END DO
                    !$omp end do
                ELSE IF (bctype == 0) THEN
                    !$omp do collapse(2) private(j, k, ibuf, area, uquer)
                    DO j = 2, jj
                        DO k = 2, kk
                            ibuf = k + (j-1)*kk
                            area = ddy(j)*ddz(k)
                            uquer = SQRT( &
                                (w(k-1, j, i3) + (w(k, j, i3) &
                                - w(k-1, j, i3))/ddz(k)*ddz(k-1)*0.5)**2 &
                                + (v(k, j-1, i3) + (v(k, j, i3) &
                                - v(k, j-1, i3))/ddy(j)*ddy(j-1)*0.5)**2)
                            qtu(k, j, istag2) = -dir*qwallfix(tbuf(ibuf), &
                                t(k, j, i3), uquer, ddx(i3), prmol)*area
                        END DO
                    END DO
                    !$omp end do
                ELSE
                    !$omp do collapse(2) private(j, k, ibuf, area)
                    DO j = 1, jj
                        DO k = 1, kk
                            ibuf = k + (j-1)*kk
                            area = ddy(j)*ddz(k)
                            qtu(k, j, istag2) = -dir*tbuf(ibuf)*area
                        END DO
                    END DO
                    !$omp end do
                END IF
            CASE (3)
                gamma = gmol/rho/prmol
                !$omp do collapse(2) private(j, k, ibuf, tout, adv, diff, area)
                DO j = 1, jj
                    DO k = 1, kk
                        ibuf = k + (j-1)*kk
                        IF (bctype == 0) THEN
                            IF (-dir*u(k, j, istag2) >= 0.0) THEN
                                tout = 2.0*tbuf(ibuf) - t(k, j, i3)
                            ELSE
                                tout = t(k, j, i3)
                            END IF
                            adv = -dir*u(k, j, istag2)*0.5 &
                                *(tout + t(k, j, i3))
                            diff = gamma*(tout - t(k, j, i3))/dx(istag2)
                            area = ddy(j)*ddz(k)
                            qtu(k, j, istag2) = -dir*(adv + diff)*area
                        ELSE
                            area = ddy(j)*ddz(k)
                            qtu(k, j, istag2) = -dir*tbuf(ibuf)*area
                        END IF
                    END DO
                END DO
                !$omp end do
            CASE DEFAULT
                CALL errr(__FILE__, __LINE__)
            END SELECT

        CASE (3, 4)
            IF (iface == 3) THEN
                j3 = 3
                jstag2 = 2
                dir = -1
            ELSE
                j3 = jj-2
                jstag2 = jj-2
                dir = 1
            END IF

            SELECT CASE (ityp)
            CASE (1)
                !$omp do collapse(2) private(i, k, ibuf, area1, area2, &
                !$omp& area3, area4, arecvtot, qtot)
                DO i = 3, ii-2, 2
                    DO k = 3, kk-2, 2
                        area1 = bt(k, j3, i)
                        area2 = bt(k+1, j3, i)
                        area3 = bt(k, j3, i+1)
                        area4 = bt(k+1, j3, i+1)
                        arecvtot = area1 + area2 + area3 + area4
                        ibuf = k + (i-1)*kk
                        qtot = qtvbuf(ibuf)
                        qtv(k, jstag2, i) = divide0(area1, arecvtot)*qtot
                        qtv(k+1, jstag2, i) = &
                            divide0(area2, arecvtot)*qtot
                        qtv(k, jstag2, i+1) = &
                            divide0(area3, arecvtot)*qtot
                        qtv(k+1, jstag2, i+1) = &
                            divide0(area4, arecvtot)*qtot
                    END DO
                END DO
                !$omp end do
            CASE (2)
                IF (bctype == 0 .AND. ilesmodel == 0) THEN
                    gamma2dx = 2.0*gmol/rho/prmol/dy(jstag2)
                    !$omp do collapse(2) private(i, k, ibuf, diff)
                    DO i = 1, ii
                        DO k = 1, kk
                            ibuf = k + (i-1)*kk
                            diff = gamma2dx*(tbuf(ibuf) - t(k, j3, i))
                            qtv(k, jstag2, i) = &
                                -dir*diff*ddx(i)*ddz(k)
                        END DO
                    END DO
                    !$omp end do
                ELSE IF (bctype == 0) THEN
                    !$omp do collapse(2) private(i, k, ibuf, area, uquer)
                    DO i = 2, ii
                        DO k = 2, kk
                            ibuf = k + (i-1)*kk
                            area = ddx(i)*ddz(k)
                            uquer = SQRT( &
                                (w(k-1, j3, i) + (w(k, j3, i) &
                                - w(k-1, j3, i))/ddz(k)*ddz(k-1)*0.5)**2 &
                                + (u(k, j3, i-1) + (u(k, j3, i) &
                                - u(k, j3, i-1))/ddx(i)*ddx(i-1)*0.5)**2)
                            qtv(k, jstag2, i) = -dir*qwallfix(tbuf(ibuf), &
                                t(k, j3, i), uquer, ddy(j3), prmol)*area
                        END DO
                    END DO
                    !$omp end do
                ELSE
                    !$omp do collapse(2) private(i, k, ibuf, area)
                    DO i = 1, ii
                        DO k = 1, kk
                            ibuf = k + (i-1)*kk
                            area = ddx(i)*ddz(k)
                            qtv(k, jstag2, i) = -dir*tbuf(ibuf)*area
                        END DO
                    END DO
                    !$omp end do
                END IF
            CASE (3)
                gamma = gmol/rho/prmol
                !$omp do collapse(2) private(i, k, ibuf, tout, adv, diff, area)
                DO i = 1, ii
                    DO k = 1, kk
                        ibuf = k + (i-1)*kk
                        IF (bctype == 0) THEN
                            IF (-dir*v(k, jstag2, i) >= 0.0) THEN
                                tout = 2.0*tbuf(ibuf) - t(k, j3, i)
                            ELSE
                                tout = t(k, j3, i)
                            END IF
                            adv = -dir*v(k, jstag2, i)*0.5 &
                                *(tout + t(k, j3, i))
                            diff = gamma*(tout - t(k, j3, i))/dy(jstag2)
                            area = ddx(i)*ddz(k)
                            qtv(k, jstag2, i) = -dir*(adv + diff)*area
                        ELSE
                            area = ddx(i)*ddz(k)
                            qtv(k, jstag2, i) = -dir*tbuf(ibuf)*area
                        END IF
                    END DO
                END DO
                !$omp end do
            CASE DEFAULT
                CALL errr(__FILE__, __LINE__)
            END SELECT

        CASE (5, 6)
            IF (iface == 5) THEN
                k3 = 3
                kstag2 = 2
                dir = -1
            ELSE
                k3 = kk-2
                kstag2 = kk-2
                dir = 1
            END IF

            SELECT CASE (ityp)
            CASE (1)
                !$omp do collapse(2) private(i, j, ibuf, area1, area2, &
                !$omp& area3, area4, arecvtot, qtot)
                DO i = 3, ii-2, 2
                    DO j = 3, jj-2, 2
                        area1 = bt(k3, j, i)
                        area2 = bt(k3, j+1, i)
                        area3 = bt(k3, j, i+1)
                        area4 = bt(k3, j+1, i+1)
                        arecvtot = area1 + area2 + area3 + area4
                        ibuf = j + (i-1)*jj
                        qtot = qtwbuf(ibuf)
                        qtw(kstag2, j, i) = divide0(area1, arecvtot)*qtot
                        qtw(kstag2, j+1, i) = &
                            divide0(area2, arecvtot)*qtot
                        qtw(kstag2, j, i+1) = &
                            divide0(area3, arecvtot)*qtot
                        qtw(kstag2, j+1, i+1) = &
                            divide0(area4, arecvtot)*qtot
                    END DO
                END DO
                !$omp end do
            CASE (2)
                IF (bctype == 0 .AND. ilesmodel == 0) THEN
                    gamma2dx = 2.0*gmol/rho/prmol/dz(kstag2)
                    !$omp do collapse(2) private(i, j, ibuf, diff)
                    DO i = 1, ii
                        DO j = 1, jj
                            ibuf = j + (i-1)*jj
                            diff = gamma2dx*(tbuf(ibuf) - t(k3, j, i))
                            qtw(kstag2, j, i) = &
                                -dir*diff*ddx(i)*ddy(j)
                        END DO
                    END DO
                    !$omp end do
                ELSE IF (bctype == 0) THEN
                    !$omp do collapse(2) private(i, j, ibuf, area, uquer)
                    DO i = 2, ii
                        DO j = 2, jj
                            ibuf = j + (i-1)*jj
                            area = ddx(i)*ddy(j)
                            uquer = SQRT( &
                                (u(k3, j, i-1) + (u(k3, j, i) &
                                - u(k3, j, i-1))/ddx(i)*ddx(i-1)*0.5)**2 &
                                + (v(k3, j-1, i) + (v(k3, j, i) &
                                - v(k3, j-1, i))/ddy(j)*ddy(j-1)*0.5)**2)
                            qtw(kstag2, j, i) = -dir*qwallfix(tbuf(ibuf), &
                                t(k3, j, i), uquer, ddz(k3), prmol)*area
                        END DO
                    END DO
                    !$omp end do
                ELSE
                    !$omp do collapse(2) private(i, j, ibuf, area)
                    DO i = 1, ii
                        DO j = 1, jj
                            ibuf = j + (i-1)*jj
                            area = ddx(i)*ddy(j)
                            qtw(kstag2, j, i) = -dir*tbuf(ibuf)*area
                        END DO
                    END DO
                    !$omp end do
                END IF
            CASE (3)
                gamma = gmol/rho/prmol
                !$omp do collapse(2) private(i, j, ibuf, tout, adv, diff, area)
                DO i = 1, ii
                    DO j = 1, jj
                        ibuf = j + (i-1)*jj
                        IF (bctype == 0) THEN
                            IF (-dir*w(kstag2, j, i) >= 0.0) THEN
                                tout = 2.0*tbuf(ibuf) - t(k3, j, i)
                            ELSE
                                tout = t(k3, j, i)
                            END IF
                            adv = -dir*w(kstag2, j, i)*0.5 &
                                *(tout + t(k3, j, i))
                            diff = gamma*(tout - t(k3, j, i))/dz(kstag2)
                            area = ddx(i)*ddy(j)
                            qtw(kstag2, j, i) = -dir*(adv + diff)*area
                        ELSE
                            area = ddx(i)*ddy(j)
                            qtw(kstag2, j, i) = -dir*tbuf(ibuf)*area
                        END IF
                    END DO
                END DO
                !$omp end do
            CASE DEFAULT
                CALL errr(__FILE__, __LINE__)
            END SELECT
        CASE DEFAULT
            CALL errr(__FILE__, __LINE__)
        END SELECT
    END SUBROUTINE bound_scalar_device


    SUBROUTINE bfront(igrid, iface, ibocd, ctyp, f1, f2, f3, f4, timeph)
        ! Subroutine arguments
        INTEGER(intk), INTENT(in) :: igrid, iface, ibocd
        CHARACTER(len=*), INTENT(in) :: ctyp
        TYPE(field_t), INTENT(inout) :: f1
        TYPE(field_t), INTENT(inout), OPTIONAL :: f2, f3, f4
        REAL(realk), INTENT(in), OPTIONAL :: timeph

        ! Local variables
        INTEGER(intk) :: kk, jj, ii
        INTEGER(intk) :: k, j, i3, istag2, dir
        REAL(realk) :: area1, area2, area3, area4, arecvtot, qtot
        REAL(realk) :: area, prmol, adv, diff, gamma, gamma2dx, uquer, tout
        INTEGER(intk) :: idx, scbtype(nsca)
        REAL(realk), POINTER, CONTIGUOUS :: qtu(:, :, :), t(:, :, :), &
            bt(:, :, :), u(:, :, :), v(:, :, :), w(:, :, :)
        REAL(realk), POINTER, CONTIGUOUS :: qtubuf(:, :), tbuf(:, :)
        REAL(realk), POINTER, CONTIGUOUS :: dx(:), ddx(:), ddy(:), ddz(:)

        ! Return early when no action is to be taken
        SELECT CASE (ctyp)
        CASE ("PAR", "SIO", "SWA")
            CONTINUE
        CASE DEFAULT
            RETURN
        END SELECT

        ! Assure that required fields are present
        IF (.NOT. PRESENT(f2) .OR. .NOT. PRESENT(f3) &
                .OR. .NOT. PRESENT(f4)) THEN
            CALL errr(__FILE__, __LINE__)
        END IF

        ! Fetch pointers
        CALL f1%get_ptr(qtu, igrid)
        ! CALL f2%get_ptr(qtv, igrid)
        ! CALL f3%get_ptr(qtw, igrid)
        CALL f4%get_ptr(t, igrid)

        CALL f1%get_buffer(qtubuf, igrid, iface)
        ! CALL f2%get_buffer(qtvbuf, igrid, iface)
        ! CALL f3%get_buffer(qtwbuf, igrid, iface)
        CALL f4%get_buffer(tbuf, igrid, iface)

        CALL get_fieldptr(u, "U", igrid)
        CALL get_fieldptr(v, "V", igrid)
        CALL get_fieldptr(w, "W", igrid)

        CALL get_fieldptr(bt, "BT", igrid)

        CALL get_fieldptr(dx, "DX", igrid)
        CALL get_fieldptr(ddx, "DDX", igrid)
        CALL get_fieldptr(ddy, "DDY", igrid)
        CALL get_fieldptr(ddz, "DDZ", igrid)

        CALL get_mgdims(kk, jj, ii, igrid)

        SELECT CASE (iface)
        CASE (1)
            ! Front
            ! i2 = 2
            i3 = 3
            ! i4 = 4
            ! istag1 = 1
            istag2 = 2
            dir = -1
        CASE (2)
            ! Back
            ! i2 = ii - 1
            i3 = ii - 2
            ! i4 = ii - 3
            ! istag1 = ii - 1
            istag2 = ii - 2
            dir = 1
        CASE DEFAULT
            CALL errr(__FILE__, __LINE__)
        END SELECT

        SELECT CASE (ctyp)
        CASE ("PAR")
            ! Wall-normal fluxes
            DO j = 3, jj-2, 2
                DO k = 3, kk-2, 2
                    ! Open areas on sides of the 4 receiver cells
                    area1 = bt(k, j, i3)
                    area2 = bt(k+1, j, i3)
                    area3 = bt(k, j+1, i3)
                    area4 = bt(k+1, j+1, i3)
                    arecvtot = area1 + area2 + area3 + area4

                    ! Coarse grid flux
                    qtot = qtubuf(k, j)

                    ! Distributing flux proportionally to receiver cells
                    qtu(k, j, istag2) = divide0(area1, arecvtot)*qtot
                    qtu(k+1, j, istag2) = divide0(area2, arecvtot)*qtot
                    qtu(k, j+1, istag2) = divide0(area3, arecvtot)*qtot
                    qtu(k+1, j+1, istag2) = divide0(area4, arecvtot)*qtot
                END DO
            END DO
        CASE ("SWA")
            CALL get_bcprms(scbtype, igrid, iface, ibocd)
            CALL f4%get_attr(idx, "SCAIDX")

            SELECT CASE (scbtype(idx))
            CASE (0)  ! Fixed scalar value (wall)
                CALL f4%get_attr(prmol, "PRMOL")

                IF (ilesmodel == 0) THEN
                    gamma2dx = 2.0 * gmol / rho / prmol / dx(istag2)
                    DO j = 1, jj
                        DO k = 1, kk
                            ! Setting the scalar diffusive flux from the wall
                            ! Wall buffer tbuf contains set scalar value
                            diff = gamma2dx*(tbuf(k, j) - t(k, j, i3))
                            qtu(k, j, istag2) = -dir*diff*ddy(j)*ddz(k)
                        END DO
                    END DO
                ELSE
                    DO j = 2, jj
                        DO k = 2, kk
                            ! Setting the scalar flux with a wall model
                            ! Wall buffer tbuf contains set scalar value
                            area = ddy(j)*ddz(k)
                            uquer = SQRT( &
                                (w(k-1, j, i3) + (w(k, j, i3)-w(k-1, j, i3)) &
                                    /ddz(k)*ddz(k-1)*0.5)**2.0 + &
                                (v(k, j-1, i3) + (v(k, j, i3)-v(k, j-1, i3)) &
                                    /ddy(j)*ddy(j-1)*0.5)**2.0)
                            qtu(k, j, istag2) = -dir*qwallfix(tbuf(k, j), &
                                t(k, j, i3), uquer, ddx(i3), prmol)*area
                        END DO
                    END DO
                END IF
            CASE (1)  ! Fixed flux value (wall)
                DO j = 1, jj
                    DO k = 1, kk
                        ! Wall buffer tbuf contains flux at this boundary
                        area = ddy(j)*ddz(k)
                        qtu(k, j, istag2) = -dir*tbuf(k, j)*area
                    END DO
                END DO
            CASE DEFAULT
                CALL errr(__FILE__, __LINE__)
            END SELECT

        CASE ("SIO")
            CALL get_bcprms(scbtype, igrid, iface, ibocd)
            CALL f4%get_attr(idx, "SCAIDX")

            SELECT CASE (scbtype(idx))
            CASE (0)  ! Fixed scalar value (inflow/outflow)
                CALL f4%get_attr(prmol, "PRMOL")
                gamma = gmol / rho / prmol

                DO j = 1, jj
                    DO k = 1, kk
                        IF (-dir*u(k, j, istag2) >= 0.0) THEN
                            ! flow into the domain (requires specified value)
                            tout = 2.0*tbuf(k, j) - t(k, j, i3)
                        ELSE
                            ! flow out of the domain (zero-gradient)
                            tout = t(k, j, i3)
                        END IF

                        ! Adding to the flux (diffusive and convective)
                        ! adv and diff are fluxes /into/ domain!
                        adv = -dir*u(k, j, istag2)*0.5*(tout + t(k, j, i3))
                        diff = gamma*(tout - t(k, j, i3))/dx(istag2)
                        area = ddy(j)*ddz(k)

                        ! flux in /coordinate direction/ (ref. fluxbalance)
                        qtu(k, j, istag2) = -dir*(adv + diff)*area
                    END DO
                END DO
            CASE (1)  ! Fixed flux value
                DO j = 1, jj
                    DO k = 1, kk
                        ! Wall buffer tbuf contains flux at this boundary
                        area = ddy(j)*ddz(k)
                        qtu(k, j, istag2) = -dir*tbuf(k, j)*area
                    END DO
                END DO
            CASE DEFAULT
                CALL errr(__FILE__, __LINE__)
            END SELECT
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
        INTEGER(intk) :: kk, jj, ii
        INTEGER(intk) :: k, i, j3, jstag2, dir
        REAL(realk) :: area1, area2, area3, area4, arecvtot, qtot
        REAL(realk) :: area, prmol, adv, diff, gamma, gamma2dx, uquer, tout
        INTEGER(intk) :: idx, scbtype(nsca)
        REAL(realk), POINTER, CONTIGUOUS :: qtv(:, :, :), t(:, :, :), &
            bt(:, :, :), u(:, :, :), v(:, :, :), w(:, :, :)
        REAL(realk), POINTER, CONTIGUOUS :: qtvbuf(:, :), tbuf(:, :)
        REAL(realk), POINTER, CONTIGUOUS :: dy(:), ddx(:), ddy(:), ddz(:)

        ! Return early when no action is to be taken
        SELECT CASE (ctyp)
        CASE ("PAR", "SIO", "SWA")
            CONTINUE
        CASE DEFAULT
            RETURN
        END SELECT

        ! Assure that required fields are present
        IF (.NOT. PRESENT(f2) .OR. .NOT. PRESENT(f3) &
                .OR. .NOT. PRESENT(f4)) THEN
            CALL errr(__FILE__, __LINE__)
        END IF

        ! Fetch pointers
        ! CALL f1%get_ptr(qtu, igrid)
        CALL f2%get_ptr(qtv, igrid)
        ! CALL f3%get_ptr(qtw, igrid)
        CALL f4%get_ptr(t, igrid)

        ! CALL f1%get_buffer(qtubuf, igrid, iface)
        CALL f2%get_buffer(qtvbuf, igrid, iface)
        ! CALL f3%get_buffer(qtwbuf, igrid, iface)
        CALL f4%get_buffer(tbuf, igrid, iface)

        CALL get_fieldptr(u, "U", igrid)
        CALL get_fieldptr(v, "V", igrid)
        CALL get_fieldptr(w, "W", igrid)

        CALL get_fieldptr(bt, "BT", igrid)

        CALL get_fieldptr(dy, "DY", igrid)
        CALL get_fieldptr(ddx, "DDX", igrid)
        CALL get_fieldptr(ddy, "DDY", igrid)
        CALL get_fieldptr(ddz, "DDZ", igrid)

        CALL get_mgdims(kk, jj, ii, igrid)

        SELECT CASE (iface)
        CASE (3)
            ! Right
            ! j2 = 2
            j3 = 3
            ! j4 = 4
            ! jstag1 = 1
            jstag2 = 2
            dir = -1
        CASE (4)
            ! Left
            ! j2 = jj - 1
            j3 = jj - 2
            ! j4 = jj - 3
            ! jstag1 = jj - 1
            jstag2 = jj - 2
            dir = 1
        CASE DEFAULT
            CALL errr(__FILE__, __LINE__)
        END SELECT

        SELECT CASE (ctyp)
        CASE ("PAR")
            ! Wall-normal fluxes
            DO i = 3, ii-2, 2
                DO k = 3, kk-2, 2
                    ! Open areas on sides of the 4 receiver cells
                    area1 = bt(k, j3, i)
                    area2 = bt(k+1, j3, i)
                    area3 = bt(k, j3, i+1)
                    area4 = bt(k+1, j3, i+1)
                    arecvtot = area1 + area2 + area3 + area4

                    ! Coarse grid flux
                    qtot = qtvbuf(k, i)

                    ! Distributing flux proportionally to receiver cells
                    qtv(k, jstag2, i) = divide0(area1, arecvtot)*qtot
                    qtv(k+1, jstag2, i) = divide0(area2, arecvtot)*qtot
                    qtv(k, jstag2, i+1) = divide0(area3, arecvtot)*qtot
                    qtv(k+1, jstag2, i+1) = divide0(area4, arecvtot)*qtot
                END DO
            END DO
        CASE ("SWA")
            CALL get_bcprms(scbtype, igrid, iface, ibocd)
            CALL f4%get_attr(idx, "SCAIDX")

            SELECT CASE (scbtype(idx))
            CASE (0)  ! Fixed scalar value (wall)
                CALL f4%get_attr(prmol, "PRMOL")

                IF (ilesmodel == 0) THEN
                    gamma2dx = 2.0 * gmol / rho / prmol / dy(jstag2)
                    DO i = 1, ii
                        DO k = 1, kk
                            ! Setting the scalar diffusive flux from the wall
                            ! Wall buffer tbuf contains set scalar value
                            diff = gamma2dx*(tbuf(k, i) - t(k, j3, i))
                            qtv(k, jstag2, i) = -dir*diff*ddx(i)*ddz(k)
                        END DO
                    END DO
                ELSE
                    DO i = 2, ii
                        DO k = 2, kk
                            ! Setting the scalar flux with a wall model
                            ! Wall buffer tbuf contains set scalar value
                            area = ddx(i)*ddz(k)
                            uquer = SQRT( &
                                (w(k-1, j3, i) + (w(k, j3, i)-w(k-1, j3, i)) &
                                    /ddz(k)*ddz(k-1)*0.5)**2.0 + &
                                (u(k, j3, i-1) + (u(k, j3, i)-u(k, j3, i-1)) &
                                    /ddx(i)*ddx(i-1)*0.5)**2.0)
                            qtv(k, jstag2, i) = -dir*qwallfix(tbuf(k, i), &
                                t(k, j3, i), uquer, ddy(j3), prmol)*area
                        END DO
                    END DO
                END IF
            CASE (1)  ! Fixed flux value (wall)
                DO i = 1, ii
                    DO k = 1, kk
                        ! Wall buffer tbuf contains flux at this boundary
                        area = ddx(i)*ddz(k)
                        qtv(k, jstag2, i) = -dir*tbuf(k, i)*area
                    END DO
                END DO
            CASE DEFAULT
                CALL errr(__FILE__, __LINE__)
            END SELECT

        CASE ("SIO")
            CALL get_bcprms(scbtype, igrid, iface, ibocd)
            CALL f4%get_attr(idx, "SCAIDX")

            SELECT CASE (scbtype(idx))
            CASE (0)  ! Fixed scalar value (inflow/outflow)
                CALL f4%get_attr(prmol, "PRMOL")
                gamma = gmol / rho / prmol

                DO i = 1, ii
                    DO k = 1, kk
                        IF (-dir*v(k, jstag2, i) >= 0.0) THEN
                            ! flow into the domain (requires specified value)
                            tout = 2.0*tbuf(k, i) - t(k, j3, i)
                        ELSE
                            ! flow out of the domain (zero-gradient)
                            tout = t(k, j3, i)
                        END IF

                        ! Adding to the flux (diffusive and convective)
                        ! adv and diff are fluxes /into/ domain!
                        adv = -dir*v(k, jstag2, i)*0.5*(tout + t(k, j3, i))
                        diff = gamma*(tout - t(k, j3, i))/dy(jstag2)
                        area = ddx(i)*ddz(k)

                        ! flux in /coordinate direction/ (ref. fluxbalance)
                        qtv(k, jstag2, i) = -dir*(adv + diff)*area
                    END DO
                END DO
            CASE (1)  ! Fixed flux value
                DO i = 1, ii
                    DO k = 1, kk
                        ! Wall buffer tbuf contains flux at this boundary
                        area = ddx(i)*ddz(k)
                        qtv(k, jstag2, i) = -dir*tbuf(k, i)*area
                    END DO
                END DO
            CASE DEFAULT
                CALL errr(__FILE__, __LINE__)
            END SELECT
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
        INTEGER(intk) :: kk, jj, ii
        INTEGER(intk) :: j, i, k3, kstag2, dir
        REAL(realk) :: area1, area2, area3, area4, arecvtot, qtot
        REAL(realk) :: area, prmol, adv, diff, gamma, gamma2dx, uquer, tout
        INTEGER(intk) :: idx, scbtype(nsca)
        REAL(realk), POINTER, CONTIGUOUS :: qtw(:, :, :), t(:, :, :), &
            bt(:, :, :), u(:, :, :), v(:, :, :), w(:, :, :)
        REAL(realk), POINTER, CONTIGUOUS :: qtwbuf(:, :), tbuf(:, :)
        REAL(realk), POINTER, CONTIGUOUS :: dz(:), ddx(:), ddy(:), ddz(:)

        ! Return early when no action is to be taken
        SELECT CASE (ctyp)
        CASE ("PAR", "SIO", "SWA")
            CONTINUE
        CASE DEFAULT
            RETURN
        END SELECT

        ! Assure that required fields are present
        IF (.NOT. PRESENT(f2) .OR. .NOT. PRESENT(f3) &
                .OR. .NOT. PRESENT(f4)) THEN
            CALL errr(__FILE__, __LINE__)
        END IF

        ! Fetch pointers
        ! CALL f1%get_ptr(qtu, igrid)
        ! CALL f2%get_ptr(qtv, igrid)
        CALL f3%get_ptr(qtw, igrid)
        CALL f4%get_ptr(t, igrid)

        ! CALL f1%get_buffer(qtubuf, igrid, iface)
        ! CALL f2%get_buffer(qtvbuf, igrid, iface)
        CALL f3%get_buffer(qtwbuf, igrid, iface)
        CALL f4%get_buffer(tbuf, igrid, iface)

        CALL get_fieldptr(u, "U", igrid)
        CALL get_fieldptr(v, "V", igrid)
        CALL get_fieldptr(w, "W", igrid)

        CALL get_fieldptr(bt, "BT", igrid)

        CALL get_fieldptr(dz, "DZ", igrid)
        CALL get_fieldptr(ddx, "DDX", igrid)
        CALL get_fieldptr(ddy, "DDY", igrid)
        CALL get_fieldptr(ddz, "DDZ", igrid)

        CALL get_mgdims(kk, jj, ii, igrid)

        SELECT CASE (iface)
        CASE (5)
            ! Bottom
            ! k2 = 2
            k3 = 3
            ! k4 = 4
            ! kstag1 = 1
            kstag2 = 2
            dir = -1
        CASE (6)
            ! Top
            ! k2 = kk - 1
            k3 = kk - 2
            ! k4 = kk - 3
            ! kstag1 = kk - 1
            kstag2 = kk - 2
            dir = 1
        CASE DEFAULT
            CALL errr(__FILE__, __LINE__)
        END SELECT

        SELECT CASE (ctyp)
        CASE ("PAR")
            ! Wall-normal fluxes
            DO i = 3, ii-2, 2
                DO j = 3, jj-2, 2
                    ! Open areas on sides of the 4 receiver cells
                    area1 = bt(k3, j, i)
                    area2 = bt(k3, j+1, i)
                    area3 = bt(k3, j, i+1)
                    area4 = bt(k3, j+1, i+1)
                    arecvtot = area1 + area2 + area3 + area4

                    ! Coarse grid flux
                    qtot = qtwbuf(j, i)

                    ! Distributing flux proportionally to receiver cells
                    qtw(kstag2, j, i) = divide0(area1, arecvtot)*qtot
                    qtw(kstag2, j+1, i) = divide0(area2, arecvtot)*qtot
                    qtw(kstag2, j, i+1) = divide0(area3, arecvtot)*qtot
                    qtw(kstag2, j+1, i+1) = divide0(area4, arecvtot)*qtot
                END DO
            END DO
        CASE ("SWA")
            CALL get_bcprms(scbtype, igrid, iface, ibocd)
            CALL f4%get_attr(idx, "SCAIDX")

            SELECT CASE (scbtype(idx))
            CASE (0)  ! Fixed scalar value (wall)
                CALL f4%get_attr(prmol, "PRMOL")

                IF (ilesmodel == 0) THEN
                    gamma2dx = 2.0 * gmol / rho / prmol / dz(kstag2)
                    DO i = 1, ii
                        DO j = 1, jj
                            ! Setting the scalar diffusive flux from the wall
                            ! Wall buffer tbuf contains set scalar value
                            diff = gamma2dx*(tbuf(j, i) - t(k3, j, i))
                            qtw(kstag2, j, i) = -dir*diff*ddx(i)*ddy(j)
                        END DO
                    END DO
                ELSE
                    DO i = 2, ii
                        DO j = 2, jj
                            ! Setting the scalar flux with a wall model
                            ! Wall buffer tbuf contains set scalar value
                            area = ddx(i)*ddy(j)
                            uquer = SQRT( &
                                (u(k3, j, i-1) + (u(k3, j, i)-u(k3, j, i-1)) &
                                    /ddx(i)*ddx(i-1)*0.5)**2.0 + &
                                (v(k3, j-1, i) + (v(k3, j, i)-v(k3, j-1, i)) &
                                    /ddy(j)*ddy(j-1)*0.5)**2.0)
                            qtw(kstag2, j, i) = -dir*qwallfix(tbuf(j, i), &
                                t(k3, j, i), uquer, ddz(k3), prmol)*area
                        END DO
                    END DO
                END IF
            CASE (1)  ! Fixed flux value (wall)
                DO i = 1, ii
                    DO j = 1, jj
                        ! Wall buffer tbuf contains flux at this boundary
                        area = ddx(i)*ddy(j)
                        qtw(kstag2, j, i) = -dir*tbuf(j, i)*area
                    END DO
                END DO
            CASE DEFAULT
                CALL errr(__FILE__, __LINE__)
            END SELECT

        CASE ("SIO")
            CALL get_bcprms(scbtype, igrid, iface, ibocd)
            CALL f4%get_attr(idx, "SCAIDX")

            SELECT CASE (scbtype(idx))
            CASE (0)  ! Fixed scalar value (inflow/outflow)
                CALL f4%get_attr(prmol, "PRMOL")
                gamma = gmol / rho / prmol

                DO i = 1, ii
                    DO j = 1, jj
                        IF (-dir*w(kstag2, j, i) >= 0.0) THEN
                            ! flow into the domain (requires specified value)
                            tout = 2.0*tbuf(j, i) - t(k3, j, i)
                        ELSE
                            ! flow out of the domain (zero-gradient)
                            tout = t(k3, j, i)
                        END IF

                        ! Adding to the flux (diffusive and convective)
                        ! adv and diff are fluxes /into/ domain!
                        adv = -dir*w(kstag2, j, i)*0.5*(tout + t(k3, j, i))
                        diff = gamma*(tout - t(k3, j, i))/dz(kstag2)
                        area = ddx(i)*ddy(j)

                        ! flux in /coordinate direction/ (ref. fluxbalance)
                        qtw(kstag2, j, i) = -dir*(adv + diff)*area
                    END DO
                END DO
            CASE (1)  ! Fixed flux value
                DO i = 1, ii
                    DO j = 1, jj
                        ! Wall buffer tbuf contains flux at this boundary
                        area = ddx(i)*ddy(j)
                        qtw(kstag2, j, i) = -dir*tbuf(j, i)*area
                    END DO
                END DO
            CASE DEFAULT
                CALL errr(__FILE__, __LINE__)
            END SELECT
        END SELECT
    END SUBROUTINE bbottom
END MODULE bound_scalar_mod
