    SUBROUTINE LESMODEL_IMPL(g, u, v, w, bp, dx, dy, dz, ddx, ddy, ddz, &
            rddx, rddy, rddz, model_constant)
        REAL(realk), INTENT(inout) :: g(*)
        REAL(realk), INTENT(in) :: u(*), v(*), w(*), bp(*)
        REAL(realk), INTENT(in) :: dx(*), dy(*), dz(*), ddx(*), ddy(*), ddz(*)
        REAL(realk), INTENT(in) :: rddx(*), rddy(*), rddz(*)
        REAL(realk), INTENT(in) :: model_constant

        INTEGER(intk) :: i, igrid, ip3, ipx, ipy, ipz
        INTEGER(intk) :: kk, jj, ii, nfro, nbac, nrgt, nlft, nbot, ntop

        !$omp target teams distribute firstprivate(model_constant) &
        !$omp& private(i, igrid, ip3, ipx, ipy, ipz, kk, jj, ii, nfro, &
        !$omp& nbac, nrgt, nlft, nbot, ntop)
        DO i = 1, nmygrids
            igrid = mygrids(i)

            CALL get_mgdims(kk, jj, ii, igrid)
            CALL get_mgbasb(nfro, nbac, nrgt, nlft, nbot, ntop, igrid)

            CALL get_ip3(ip3, igrid)
            CALL get_ip1x(ipx, igrid)
            CALL get_ip1y(ipy, igrid)
            CALL get_ip1z(ipz, igrid)

            !$omp parallel
            CALL LESMODEL_GRID(kk, jj, ii, nfro, nbac, nrgt, nlft, nbot, &
                ntop, dx(ipx), dy(ipy), dz(ipz), ddx(ipx), ddy(ipy), &
                ddz(ipz), rddx(ipx), rddy(ipy), rddz(ipz), u(ip3), v(ip3), &
                w(ip3), gmol, rho, bp(ip3), g(ip3), model_constant)
            !$omp end parallel
        END DO
        !$omp end target teams distribute
    END SUBROUTINE LESMODEL_IMPL


    SUBROUTINE LESMODEL_GRID(kk, jj, ii, nfro, nbac, nrgt, nlft, nbot, &
            ntop, dx, dy, dz, ddx, ddy, ddz, rddx, rddy, rddz, u, v, w, &
            gmol, rho, bp, g, model_constant)
        !$omp declare target

        INTEGER(intk), INTENT(in) :: kk, jj, ii
        INTEGER(intk), INTENT(in) :: nfro, nbac, nrgt, nlft, nbot, ntop
        REAL(realk), INTENT(in) :: dx(ii), dy(jj), dz(kk)
        REAL(realk), INTENT(in) :: ddx(ii), ddy(jj), ddz(kk)
        REAL(realk), INTENT(in) :: rddx(ii), rddy(jj), rddz(kk)
        REAL(realk), INTENT(in) :: u(kk, jj, ii), v(kk, jj, ii), w(kk, jj, ii)
        REAL(realk), INTENT(in) :: gmol, rho
        REAL(realk), INTENT(in) :: bp(kk, jj, ii)
        REAL(realk), INTENT(inout) :: g(kk, jj, ii)
        REAL(realk), INTENT(in) :: model_constant

        INTEGER(intk) :: k, j, i
        REAL(realk) :: dudx, dudy, dudz, dvdx, dvdy, dvdz
        REAL(realk) :: dwdx, dwdy, dwdz
        REAL(realk) :: dxf, dyf, dzf
        REAL(realk) :: dxf2, dyf2, dzf2
        REAL(realk) :: delta, dm

        !$omp do collapse(3) private(i, j, k, dudx, dudy, dudz, &
        !$omp& dvdx, dvdy, dvdz, dwdx, dwdy, dwdz, dxf, dyf, dzf, dxf2, &
        !$omp& dyf2, dzf2, delta, dm)
        DO i = 3, ii-2
            DO j = 3, jj-2
                DO k = 3, kk-2
                    dxf = 0.25*dx(i-1)*rddx(i)
                    dudx = (u(k, j, i) - u(k, j, i-1))*rddx(i)
                    dudy = ((u(k, j+1, i) - u(k, j-1, i))*dxf &
                        + (u(k, j+1, i-1) - u(k, j-1, i-1)) &
                        * (0.5-dxf))*rddy(j)
                    dudz = ((u(k+1, j, i) - u(k-1, j, i))*dxf &
                        + (u(k+1, j, i-1) - u(k-1, j, i-1)) &
                        * (0.5-dxf))*rddz(k)

                    dyf = 0.25*dy(j-1)*rddy(j)
                    dvdx = ((v(k, j, i+1) - v(k, j, i-1))*dyf &
                        + (v(k, j-1, i+1) - v(k, j-1, i-1)) &
                        * (0.5-dyf))*rddx(i)
                    dvdy = (v(k, j, i) - v(k, j-1, i))*rddy(j)
                    dvdz = ((v(k+1, j, i) - v(k-1, j, i))*dyf &
                        + (v(k+1, j-1, i) - v(k-1, j-1, i)) &
                        * (0.5-dyf))*rddz(k)

                    dzf = 0.25*dz(k-1)*rddz(k)
                    dwdx = ((w(k, j, i+1) - w(k, j, i-1))*dzf &
                        + (w(k-1, j, i+1) - w(k-1, j, i-1)) &
                        * (0.5-dzf))*rddx(i)
                    dwdy = ((w(k, j+1, i) - w(k, j-1, i))*dzf &
                        + (w(k-1, j+1, i) - w(k-1, j-1, i)) &
                        * (0.5-dzf))*rddy(j)
                    dwdz = (w(k, j, i) - w(k-1, j, i))*rddz(k)

                    delta = cube_root(ddx(i)*ddy(j)*ddz(k))*bp(k, j, i)

                    IF (nfro == 5 .AND. i == 3) THEN
                        dyf2 = 0.5*dy(j-1)*rddy(j)
                        dvdx = gradp2(dyf2*v(k, j, i) &
                            + (1.0-dyf2)*v(k, j-1, i), ddx(i))

                        dzf2 = 0.5*dz(k-1)*rddz(k)
                        dwdx = gradp2(dzf2*w(k, j, i) &
                            + (1.0-dzf2)*w(k-1, j, i), ddx(i))
                    END IF
                    IF (nbac == 5 .AND. i == ii-2) THEN
                        dyf2 = 0.5*dy(j-1)*rddy(j)
                        dvdx = -gradp2(dyf2*v(k, j, i) &
                            + (1.0-dyf2)*v(k, j-1, i), ddx(i))

                        dzf2 = 0.5*dz(k-1)*rddz(k)
                        dwdx = -gradp2(dzf2*w(k, j, i) &
                            + (1.0-dzf2)*w(k-1, j, i), ddx(i))
                    END IF
                    IF (nrgt == 5 .AND. j == 3) THEN
                        dxf2 = 0.5*dx(i-1)*rddx(i)
                        dudy = gradp2(dxf2*u(k, j, i) &
                            + (1.0-dxf2)*u(k, j, i-1), ddy(j))

                        dzf2 = 0.5*dz(k-1)*rddz(k)
                        dwdy = gradp2(dzf2*w(k, j, i) &
                            + (1.0-dzf2)*w(k-1, j, i), ddy(j))
                    END IF
                    IF (nlft == 5 .AND. j == jj-2) THEN
                        dxf2 = 0.5*dx(i-1)*rddx(i)
                        dudy = -gradp2(dxf2*u(k, j, i) &
                            + (1.0-dxf2)*u(k, j, i-1), ddy(j))

                        dzf2 = 0.5*dz(k-1)*rddz(k)
                        dwdy = -gradp2(dzf2*w(k, j, i) &
                            + (1.0-dzf2)*w(k-1, j, i), ddy(j))
                    END IF
                    IF (nbot == 5 .AND. k == 3) THEN
                        dxf2 = 0.5*dx(i-1)*rddx(i)
                        dudz = gradp2(dxf2*u(k, j, i) &
                            + (1.0-dxf2)*u(k, j, i-1), ddz(k))

                        dyf2 = 0.5*dy(j-1)*rddy(j)
                        dvdz = gradp2(dyf2*v(k, j, i) &
                            + (1.0-dyf2)*v(k, j-1, i), ddz(k))
                    END IF
                    IF (ntop == 5 .AND. k == kk-2) THEN
                        dxf2 = 0.5*dx(i-1)*rddx(i)
                        dudz = -gradp2(dxf2*u(k, j, i) &
                            + (1.0-dxf2)*u(k, j, i-1), ddz(k))

                        dyf2 = 0.5*dy(j-1)*rddy(j)
                        dvdz = -gradp2(dyf2*v(k, j, i) &
                            + (1.0-dyf2)*v(k, j-1, i), ddz(k))
                    END IF

                    dm = LESMODEL_EVALUATE(dudx, dudy, dudz, dvdx, dvdy, &
                        dvdz, dwdx, dwdy, dwdz, model_constant)
                    g(k, j, i) = rho*delta**2*dm + gmol
                END DO
            END DO
        END DO
        !$omp end do
    END SUBROUTINE LESMODEL_GRID
