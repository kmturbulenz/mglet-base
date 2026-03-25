MODULE bound_scalar_mod
    USE core_mod
    USE scacore_mod
    USE flow_mod, ONLY: ilesmodel, gmol, rho, qwallfix

    IMPLICIT NONE(type, external)
    PRIVATE

    PUBLIC :: bound_scaflux

CONTAINS
    SUBROUTINE bound_scaflux(ilevel, qtu_f, qtv_f, qtw_f, t_f)
        INTEGER(intk), INTENT(in) :: ilevel
        TYPE(field_t), INTENT(inout) :: qtu_f, qtv_f, qtw_f, t_f

        INTEGER(intk) :: i, igrid, iface, nbocd, ibocd, kk, jj, ii
        CHARACTER(len=8) :: ctyp
        INTEGER(intk) :: scaidx, scb, scbtype(nsca)
        REAL(realk) :: prmol
        TYPE(field_t), POINTER :: u_f, v_f, w_f, bt_f, dx_f, dy_f, dz_f, ddx_f, &
            ddy_f, ddz_f
        REAL(realk), POINTER, CONTIGUOUS :: qtu(:, :, :), qtv(:, :, :), &
            qtw(:, :, :), t(:, :, :)
        REAL(realk), POINTER, CONTIGUOUS :: qtubuf(:, :), qtvbuf(:, :), &
            qtwbuf(:, :), tbuf(:, :)
        REAL(realk), POINTER, CONTIGUOUS :: u(:, :, :), v(:, :, :), &
            w(:, :, :), bt(:, :, :), dx(:), dy(:), dz(:), &
            ddx(:), ddy(:), ddz(:)

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

        CALL t_f%get_attr(scaidx, "SCAIDX")
        CALL t_f%get_attr(prmol, "PRMOL")

        ! TODO: exploit more parallelism here!

        !$omp target teams loop bind(teams) shared(prmol, scaidx) &
        !$omp private(igrid, iface, nbocd, ibocd, ctyp, scb, scbtype, kk, jj, ii, qtu, qtv, &
        !$omp         qtw, t, qtubuf, qtvbuf, qtwbuf, tbuf, u, v, w, bt, dx, &
        !$omp         dy, dz, ddx, ddy, ddz)
        DO i = 1, nmygridslvl(ilevel)
            igrid = mygridslvl(i, ilevel)

            CALL get_mgdims(kk, jj, ii, igrid)

            CALL qtu_f%get_ptr(qtu, igrid)
            CALL qtv_f%get_ptr(qtv, igrid)
            CALL qtw_f%get_ptr(qtw, igrid)
            CALL t_f%get_ptr(t, igrid)

            CALL u_f%get_ptr(u, igrid)
            CALL v_f%get_ptr(v, igrid)
            CALL w_f%get_ptr(w, igrid)
            CALL bt_f%get_ptr(bt, igrid)
            CALL dx_f%get_ptr(dx, igrid)
            CALL dy_f%get_ptr(dy, igrid)
            CALL dz_f%get_ptr(dz, igrid)
            CALL ddx_f%get_ptr(ddx, igrid)
            CALL ddy_f%get_ptr(ddy, igrid)
            CALL ddz_f%get_ptr(ddz, igrid)

            DO iface = 1, 6
                nbocd = nboconds(iface, igrid)
                DO ibocd = 1, nbocd
                    CALL get_bc_ctyp(ctyp, ibocd, iface, igrid)

                    SELECT CASE (ctyp)
                    CASE ("PAR", "SIO", "SWA")
                        CONTINUE
                    CASE DEFAULT
                        CYCLE
                    END SELECT

                    CALL qtu_f%get_buffer(qtubuf, igrid, iface)
                    CALL qtv_f%get_buffer(qtvbuf, igrid, iface)
                    CALL qtw_f%get_buffer(qtwbuf, igrid, iface)
                    CALL t_f%get_buffer(tbuf, igrid, iface)

                    CALL get_bcprms(scbtype, igrid, iface, ibocd)
                    scb = scbtype(scaidx)

                    !$omp parallel
                    SELECT CASE(iface)
                    CASE(1, 2)
                        CALL bfront(igrid, iface, ibocd, ctyp, scb, prmol, &
                            kk, jj, ii, qtu, t, u, v, w, bt, dx, &
                            ddx, ddy, ddz, qtubuf, tbuf)
                    CASE(3, 4)
                        CALL bright(igrid, iface, ibocd, ctyp, scb, prmol, &
                            kk, jj, ii, qtv, t, u, v, w, bt, dy, &
                            ddx, ddy, ddz, qtvbuf, tbuf)
                    CASE(5, 6)
                        CALL bbottom(igrid, iface, ibocd, ctyp, scb, prmol, &
                            kk, jj, ii, qtw, t, u, v, w, bt, dz, &
                            ddx, ddy, ddz, qtwbuf, tbuf)
                    END SELECT
                    !$omp end parallel
                END DO
            END DO
        END DO
    END SUBROUTINE bound_scaflux


    SUBROUTINE bfront(igrid, iface, ibocd, ctyp, scb, prmol, kk, jj, ii, &
            qtu, t, u, v, w, bt, dx, ddx, ddy, ddz, qtubuf, tbuf)
        !$omp declare target

        ! Subroutine arguments
        INTEGER(intk), INTENT(in) :: igrid, iface, ibocd
        CHARACTER(len=*), INTENT(in) :: ctyp
        INTEGER(intk), INTENT(in) :: scb
        REAL(realk), INTENT(in) :: prmol
        INTEGER(intk), INTENT(in) :: kk, jj, ii
        REAL(realk), INTENT(inout) :: qtu(kk, jj, ii)
        REAL(realk), INTENT(in) :: t(kk, jj, ii), u(kk, jj, ii), &
            v(kk, jj, ii), w(kk, jj, ii), bt(kk, jj, ii), dx(ii), &
            ddx(ii), ddy(jj), ddz(kk), qtubuf(kk, jj), tbuf(kk, jj)

        ! Local variables
        INTEGER(intk) :: k, j, i3, istag2, dir
        REAL(realk) :: area1, area2, area3, area4, arecvtot, qtot
        REAL(realk) :: area, adv, diff, gamma, gamma2dx, uquer, tout

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
            ! CALL errr(__FILE__, __LINE__)
        END SELECT

        SELECT CASE (ctyp)
        CASE ("PAR")
            ! Wall-normal fluxes
            !$omp loop collapse(2) bind(parallel) private(k, j, area1, area2, area3, area4, arecvtot, qtot)
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
            SELECT CASE (scb)
            CASE (0)  ! Fixed scalar value (wall)
                IF (ilesmodel == 0) THEN
                    gamma2dx = 2.0 * gmol / rho / prmol / dx(istag2)
                    !$omp loop collapse(2) bind(parallel) private(k, j, diff)
                    DO j = 1, jj
                        DO k = 1, kk
                            ! Setting the scalar diffusive flux from the wall
                            ! Wall buffer tbuf contains set scalar value
                            diff = gamma2dx*(tbuf(k, j) - t(k, j, i3))
                            qtu(k, j, istag2) = -dir*diff*ddy(j)*ddz(k)
                        END DO
                    END DO
                ELSE
                    !$omp loop collapse(2) bind(parallel) private(k, j, area, uquer)
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
                !$omp loop collapse(2) bind(parallel) private(k, j, area)
                DO j = 1, jj
                    DO k = 1, kk
                        ! Wall buffer tbuf contains flux at this boundary
                        area = ddy(j)*ddz(k)
                        qtu(k, j, istag2) = -dir*tbuf(k, j)*area
                    END DO
                END DO
            CASE DEFAULT
                ! CALL errr(__FILE__, __LINE__)
            END SELECT

        CASE ("SIO")
            SELECT CASE (scb)
            CASE (0)  ! Fixed scalar value (inflow/outflow)
                gamma = gmol / rho / prmol

                !$omp loop collapse(2) bind(parallel) private(k, j, tout, adv, diff, area)
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
                !$omp loop collapse(2) bind(parallel) private(k, j, area)
                DO j = 1, jj
                    DO k = 1, kk
                        ! Wall buffer tbuf contains flux at this boundary
                        area = ddy(j)*ddz(k)
                        qtu(k, j, istag2) = -dir*tbuf(k, j)*area
                    END DO
                END DO
            CASE DEFAULT
                ! CALL errr(__FILE__, __LINE__)
            END SELECT
        END SELECT
    END SUBROUTINE bfront


    SUBROUTINE bright(igrid, iface, ibocd, ctyp, scb, prmol, kk, jj, ii, &
            qtv, t, u, v, w, bt, dy, ddx, ddy, ddz, qtvbuf, tbuf)
        !$omp declare target

        ! Subroutine arguments
        INTEGER(intk), INTENT(in) :: igrid, iface, ibocd
        CHARACTER(len=*), INTENT(in) :: ctyp
        INTEGER(intk), INTENT(in) :: scb
        REAL(realk), INTENT(in) :: prmol
        INTEGER(intk), INTENT(in) :: kk, jj, ii
        REAL(realk), INTENT(inout) :: qtv(kk, jj, ii)
        REAL(realk), INTENT(in) :: t(kk, jj, ii), u(kk, jj, ii), &
            v(kk, jj, ii), w(kk, jj, ii), bt(kk, jj, ii), dy(jj), &
            ddx(ii), ddy(jj), ddz(kk), qtvbuf(kk, ii), tbuf(kk, ii)

        ! Local variables
        INTEGER(intk) :: k, i, j3, jstag2, dir
        REAL(realk) :: area1, area2, area3, area4, arecvtot, qtot
        REAL(realk) :: area, adv, diff, gamma, gamma2dx, uquer, tout

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
            ! CALL errr(__FILE__, __LINE__)
        END SELECT

        SELECT CASE (ctyp)
        CASE ("PAR")
            ! Wall-normal fluxes
            !$omp loop collapse(2) bind(parallel) private(k, i, area1, area2, area3, area4, arecvtot, qtot)
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
            SELECT CASE (scb)
            CASE (0)  ! Fixed scalar value (wall)
                IF (ilesmodel == 0) THEN
                    gamma2dx = 2.0 * gmol / rho / prmol / dy(jstag2)
                    !$omp loop collapse(2) bind(parallel) private(k, i, diff)
                    DO i = 1, ii
                        DO k = 1, kk
                            ! Setting the scalar diffusive flux from the wall
                            ! Wall buffer tbuf contains set scalar value
                            diff = gamma2dx*(tbuf(k, i) - t(k, j3, i))
                            qtv(k, jstag2, i) = -dir*diff*ddx(i)*ddz(k)
                        END DO
                    END DO
                ELSE
                    !$omp loop collapse(2) bind(parallel) private(k, i, area, uquer)
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
                !$omp loop collapse(2) bind(parallel) private(k, i, area)
                DO i = 1, ii
                    DO k = 1, kk
                        ! Wall buffer tbuf contains flux at this boundary
                        area = ddx(i)*ddz(k)
                        qtv(k, jstag2, i) = -dir*tbuf(k, i)*area
                    END DO
                END DO
            CASE DEFAULT
                ! CALL errr(__FILE__, __LINE__)
            END SELECT

        CASE ("SIO")
            SELECT CASE (scb)
            CASE (0)  ! Fixed scalar value (inflow/outflow)
                gamma = gmol / rho / prmol

                !$omp loop collapse(2) bind(parallel) private(k, i, tout, adv, diff, area)
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
                !$omp loop collapse(2) bind(parallel) private(k, i, area)
                DO i = 1, ii
                    DO k = 1, kk
                        ! Wall buffer tbuf contains flux at this boundary
                        area = ddx(i)*ddz(k)
                        qtv(k, jstag2, i) = -dir*tbuf(k, i)*area
                    END DO
                END DO
            CASE DEFAULT
                ! CALL errr(__FILE__, __LINE__)
            END SELECT
        END SELECT
    END SUBROUTINE bright


    SUBROUTINE bbottom(igrid, iface, ibocd, ctyp, scb, prmol, kk, jj, ii, &
            qtw, t, u, v, w, bt, dz, ddx, ddy, ddz, qtwbuf, tbuf)
        !$omp declare target

        ! Subroutine arguments
        INTEGER(intk), INTENT(in) :: igrid, iface, ibocd
        CHARACTER(len=*), INTENT(in) :: ctyp
        INTEGER(intk), INTENT(in) :: scb
        REAL(realk), INTENT(in) :: prmol
        INTEGER(intk), INTENT(in) :: kk, jj, ii
        REAL(realk), INTENT(inout) :: qtw(kk, jj, ii)
        REAL(realk), INTENT(in) :: t(kk, jj, ii), u(kk, jj, ii), &
            v(kk, jj, ii), w(kk, jj, ii), bt(kk, jj, ii), dz(kk), &
            ddx(ii), ddy(jj), ddz(kk), qtwbuf(jj, ii), tbuf(jj, ii)

        ! Local variables
        INTEGER(intk) :: j, i, k3, kstag2, dir
        REAL(realk) :: area1, area2, area3, area4, arecvtot, qtot
        REAL(realk) :: area, adv, diff, gamma, gamma2dx, uquer, tout

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
            ! CALL errr(__FILE__, __LINE__)
        END SELECT

        SELECT CASE (ctyp)
        CASE ("PAR")
            ! Wall-normal fluxes
            !$omp loop collapse(2) bind(parallel) private(j, i, area1, area2, area3, area4, arecvtot, qtot)
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
            SELECT CASE (scb)
            CASE (0)  ! Fixed scalar value (wall)
                IF (ilesmodel == 0) THEN
                    gamma2dx = 2.0 * gmol / rho / prmol / dz(kstag2)
                    !$omp loop collapse(2) bind(parallel) private(j, i, diff)
                    DO i = 1, ii
                        DO j = 1, jj
                            ! Setting the scalar diffusive flux from the wall
                            ! Wall buffer tbuf contains set scalar value
                            diff = gamma2dx*(tbuf(j, i) - t(k3, j, i))
                            qtw(kstag2, j, i) = -dir*diff*ddx(i)*ddy(j)
                        END DO
                    END DO
                ELSE
                    !$omp loop collapse(2) bind(parallel) private(j, i, area, uquer)
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
                !$omp loop collapse(2) bind(parallel) private(j, i, area)
                DO i = 1, ii
                    DO j = 1, jj
                        ! Wall buffer tbuf contains flux at this boundary
                        area = ddx(i)*ddy(j)
                        qtw(kstag2, j, i) = -dir*tbuf(j, i)*area
                    END DO
                END DO
            CASE DEFAULT
                ! CALL errr(__FILE__, __LINE__)
            END SELECT

        CASE ("SIO")
            SELECT CASE (scb)
            CASE (0)  ! Fixed scalar value (inflow/outflow)
                gamma = gmol / rho / prmol

                !$omp loop collapse(2) bind(parallel) private(j, i, tout, adv, diff, area)
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
                !$omp loop collapse(2) bind(parallel) private(j, i, area)
                DO i = 1, ii
                    DO j = 1, jj
                        ! Wall buffer tbuf contains flux at this boundary
                        area = ddx(i)*ddy(j)
                        qtw(kstag2, j, i) = -dir*tbuf(j, i)*area
                    END DO
                END DO
            CASE DEFAULT
                ! CALL errr(__FILE__, __LINE__)
            END SELECT
        END SELECT
    END SUBROUTINE bbottom
END MODULE bound_scalar_mod
