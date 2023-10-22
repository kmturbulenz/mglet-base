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

    PUBLIC :: bound_scaflux

CONTAINS
    SUBROUTINE bfront(igrid, iface, ibocd, ctyp, f1, f2, f3, f4, timeph)
        ! Subroutine arguments
        INTEGER(intk), INTENT(in) :: igrid, iface, ibocd
        CHARACTER(len=*), INTENT(in) :: ctyp
        TYPE(field_t), INTENT(inout) :: f1
        TYPE(field_t), INTENT(inout), OPTIONAL :: f2, f3, f4
        REAL(realk), INTENT(in), OPTIONAL :: timeph

        ! Local variables
        INTEGER(intk) :: kk, jj, ii
        INTEGER(intk) :: k, j, i2, i3, i4, istag1, istag2, dir
        REAL(realk) :: area1, area2, area3, area4, arecvtot, qtot
        REAL(realk) :: area, prmol, adv, diff, gamma, gamma2dx, uquer, tout
        INTEGER(intk) :: idx, scbtype(nsca)
        REAL(realk), POINTER, CONTIGUOUS :: qtu(:, :, :), t(:, :, :), &
            bt(:, :, :), u(:, :, :), v(:, :, :), w(:, :, :)
        REAL(realk), POINTER, CONTIGUOUS :: qtubuf(:, :, :), tbuf(:, :, :)
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

        CALL f1%buffers%get_buffer(qtubuf, igrid, iface)
        ! CALL f2%buffers%get_buffer(qtvbuf, igrid, iface)
        ! CALL f3%buffers%get_buffer(qtwbuf, igrid, iface)
        CALL f4%buffers%get_buffer(tbuf, igrid, iface)

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
                    qtot = qtubuf(k, j, 1)

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
                            diff = gamma2dx*(tbuf(k, j, 1) - t(k, j, i3))
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
                                    /ddy(j)*ddy(j-1)*0.5 )**2.0)
                            qtu(k, j, istag2) = -dir*qwallfix(tbuf(k, j, 1), &
                                t(k, j, i3), uquer, ddx(i3), prmol)*area
                        END DO
                    END DO
                END IF
            CASE (1)  ! Fixed flux value (wall)
                DO j = 1, jj
                    DO k = 1, kk
                        ! Wall buffer tbuf contains flux at this boundary
                        area = ddy(j)*ddz(k)
                        qtu(k, j, istag2) = -dir*tbuf(k, j, 1)*area
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
                            tout = 2.0*tbuf(k, j, 1) - t(k, j, i3)
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
                        qtu(k, j, istag2) = -dir*tbuf(k, j, 1)*area
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
        INTEGER(intk) :: k, i, j2, j3, j4, jstag1, jstag2, dir
        REAL(realk) :: area1, area2, area3, area4, arecvtot, qtot
        REAL(realk) :: area, prmol, adv, diff, gamma, gamma2dx, uquer, tout
        INTEGER(intk) :: idx, scbtype(nsca)
        REAL(realk), POINTER, CONTIGUOUS :: qtv(:, :, :), t(:, :, :), &
            bt(:, :, :), u(:, :, :), v(:, :, :), w(:, :, :)
        REAL(realk), POINTER, CONTIGUOUS :: qtvbuf(:, :, :), tbuf(:, :, :)
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

        ! CALL f1%buffers%get_buffer(qtubuf, igrid, iface)
        CALL f2%buffers%get_buffer(qtvbuf, igrid, iface)
        ! CALL f3%buffers%get_buffer(qtwbuf, igrid, iface)
        CALL f4%buffers%get_buffer(tbuf, igrid, iface)

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
                    qtot = qtvbuf(k, i, 1)

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
                            diff = gamma2dx*(tbuf(k, i, 1) - t(k, j3, i))
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
                            qtv(k, jstag2, i) = -dir*qwallfix(tbuf(k, i, 1), &
                                t(k, j3, i), uquer, ddy(j3), prmol)*area
                        END DO
                    END DO
                END IF
            CASE (1)  ! Fixed flux value (wall)
                DO i = 1, ii
                    DO k = 1, kk
                        ! Wall buffer tbuf contains flux at this boundary
                        area = ddx(i)*ddz(k)
                        qtv(k, jstag2, i) = -dir*tbuf(k, i, 1)*area
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
                            tout = 2.0*tbuf(k, i, 1) - t(k, j3, i)
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
                        qtv(k, jstag2, i) = -dir*tbuf(k, i, 1)*area
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
        INTEGER(intk) :: j, i, k2, k3, k4, kstag1, kstag2, dir
        REAL(realk) :: area1, area2, area3, area4, arecvtot, qtot
        REAL(realk) :: area, prmol, adv, diff, gamma, gamma2dx, uquer, tout
        INTEGER(intk) :: idx, scbtype(nsca)
        REAL(realk), POINTER, CONTIGUOUS :: qtw(:, :, :), t(:, :, :), &
            bt(:, :, :), u(:, :, :), v(:, :, :), w(:, :, :)
        REAL(realk), POINTER, CONTIGUOUS :: qtwbuf(:, :, :), tbuf(:, :, :)
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

        ! CALL f1%buffers%get_buffer(qtubuf, igrid, iface)
        ! CALL f2%buffers%get_buffer(qtvbuf, igrid, iface)
        CALL f3%buffers%get_buffer(qtwbuf, igrid, iface)
        CALL f4%buffers%get_buffer(tbuf, igrid, iface)

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
                    qtot = qtwbuf(j, i, 1)

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
                            diff = gamma2dx*(tbuf(j, i, 1) - t(k3, j, i))
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
                                    /ddy(j)*ddy(j-1)*0.5 )**2.0)
                            qtw(kstag2, j, i) = -dir*qwallfix(tbuf(j, i, 1), &
                                t(k3, j, i), uquer, ddz(k3), prmol)*area
                        END DO
                    END DO
                END IF
            CASE (1)  ! Fixed flux value (wall)
                DO i = 1, ii
                    DO j = 1, jj
                        ! Wall buffer tbuf contains flux at this boundary
                        area = ddx(i)*ddy(j)
                        qtw(kstag2, j, i) = -dir*tbuf(j, i, 1)*area
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
                            tout = 2.0*tbuf(j, i, 1) - t(k3, j, i)
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
                        qtw(kstag2, j, i) = -dir*tbuf(j, i, 1)*area
                    END DO
                END DO
            CASE DEFAULT
                CALL errr(__FILE__, __LINE__)
            END SELECT
        END SELECT
    END SUBROUTINE bbottom

END MODULE bound_scalar_mod
