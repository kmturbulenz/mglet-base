MODULE bound_pressure_mod
    USE core_mod

    IMPLICIT NONE (type, external)
    PRIVATE

    ! A-versions: simple versions not considering the BP field
    ! B-versions: IB versions using the BP field
    INTERFACE pressureftocone
        MODULE PROCEDURE :: pressureftocone_A, pressureftocone_B
    END INTERFACE pressureftocone

    ABSTRACT INTERFACE
        SUBROUTINE boundary_callback_t(ilevel_index, igrid, iface, itask)
            IMPORT :: intk
            INTEGER(intk), INTENT(in) :: ilevel_index, igrid, iface, itask
        END SUBROUTINE boundary_callback_t
    END INTERFACE

    INTEGER(intk), PARAMETER :: boundtasksize = 2

    ! Tasks for parallel processing of boundary pressure
    ! Avoids repeated loops over all faces and checks for 'PAR' boundaries
    !   Field 1: grid index
    !   Field 2: face index
    INTEGER(intk), ALLOCATABLE :: boundtasks(:, :, :)
    INTEGER(intk), ALLOCATABLE :: nboundtaskslvl(:)
    !$omp declare target(boundtasks)

    PUBLIC :: init_bound_pressure, bound_pressure, finish_bound_pressure

CONTAINS

    SUBROUTINE init_bound_pressure()
        ! Subroutine arguments
        ! none...

        ! Local variables
        INTEGER(intk) :: nlevels

        nlevels = maxlevel - minlevel + 1
        ALLOCATE(nboundtaskslvl(nlevels), source=0_intk)

        CALL for_all_par_boundaries(boundary_counter=nboundtaskslvl)
        ALLOCATE(boundtasks(boundtasksize, &
            MAXVAL(nboundtaskslvl), nlevels), source=-1_intk)

        CALL for_all_par_boundaries(on_boundary=add_par_boundary_task)

        !$omp target enter data map(always, to: boundtasks)
    END SUBROUTINE init_bound_pressure


    SUBROUTINE for_all_par_boundaries(boundary_counter, on_boundary)
        ! Subroutine arguments
        INTEGER(intk), INTENT(out), OPTIONAL :: boundary_counter(:)
        PROCEDURE(boundary_callback_t), OPTIONAL :: on_boundary

        ! Local variables
        INTEGER(intk) :: ilevel, ilevel_index, imygrid, igrid, iface
        INTEGER(intk) :: ibocd, nbocd
        INTEGER(intk) :: par_boundaries_per_level
        CHARACTER(len=8) :: ctyp

        DO ilevel = minlevel, maxlevel
            CALL level_index(ilevel_index, ilevel)
            par_boundaries_per_level = 0

            DO imygrid = 1, nmygridslvl(ilevel)
                igrid = mygridslvl(imygrid, ilevel)

                DO iface = 1, 6
                    nbocd = nboconds(iface, igrid)
                    DO ibocd = 1, nbocd
                        CALL get_bc_ctyp(ctyp, ibocd, iface, igrid)
                        IF (ctyp /= 'PAR') CYCLE

                        par_boundaries_per_level = par_boundaries_per_level + 1
                        IF (PRESENT(on_boundary)) THEN
                            CALL on_boundary(ilevel_index, igrid, iface, &
                                par_boundaries_per_level)
                        END IF
                    END DO
                END DO
            END DO

            IF (PRESENT(boundary_counter)) THEN
                boundary_counter(ilevel_index) = par_boundaries_per_level
            END IF
        END DO
    END SUBROUTINE for_all_par_boundaries


    SUBROUTINE add_par_boundary_task(ilevel_index, igrid, iface, itask)
        ! Subroutine arguments
        INTEGER(intk), INTENT(in) :: ilevel_index, igrid, iface, itask

        ! Local variables
        ! none...

        boundtasks(1, itask, ilevel_index) = igrid
        boundtasks(2, itask, ilevel_index) = iface
    END SUBROUTINE add_par_boundary_task


    SUBROUTINE bound_pressure(ilevel, dp_f, bp_f)
        ! Subroutine arguments
        INTEGER(intk), INTENT(in) :: ilevel
        TYPE(field_t), INTENT(inout) :: dp_f
        TYPE(field_t), INTENT(in), OPTIONAL :: bp_f

        ! Local variables
        TYPE(field_t), POINTER :: dx_f, dy_f, dz_f, ddx_f, ddy_f, ddz_f
        LOGICAL :: use_bp

        use_bp = PRESENT(bp_f)

        CALL get_field(dx_f, "DX")
        CALL get_field(dy_f, "DY")
        CALL get_field(dz_f, "DZ")
        CALL get_field(ddx_f, "DDX")
        CALL get_field(ddy_f, "DDY")
        CALL get_field(ddz_f, "DDZ")

        IF (use_bp) THEN
            CALL bound_pressure_impl_bp(ilevel, dp_f, bp_f, dx_f, dy_f, dz_f, &
                ddx_f, ddy_f, ddz_f)
        ELSE
            CALL bound_pressure_impl_nobp(ilevel, dp_f, dx_f, dy_f, dz_f, &
                ddx_f, ddy_f, ddz_f)
        END IF
    END SUBROUTINE bound_pressure


    SUBROUTINE bound_pressure_impl_bp(ilevel, dp_f, bp_f, dx_f, dy_f, dz_f, &
        ddx_f, ddy_f, ddz_f)
        ! Subroutine arguments
        INTEGER(intk), INTENT(in) :: ilevel
        TYPE(field_t), INTENT(inout) :: dp_f
        TYPE(field_t), INTENT(in) :: bp_f
        TYPE(field_t), POINTER, INTENT(in) :: dx_f, dy_f, dz_f
        TYPE(field_t), POINTER, INTENT(in) :: ddx_f, ddy_f, ddz_f

        ! Local variables
        INTEGER(intk) :: nboundtasks, ilevel_index, i, igrid, iface
        INTEGER(intk) :: kk, jj, ii, ip3, ipx, ipy, ipz, ipbb
        REAL(realk), ALLOCATABLE, DIMENSION(:) :: p, pbuffer, bp, dx, dy, dz, &
            ddx, ddy, ddz

        CALL level_index(ilevel_index, ilevel)
        nboundtasks = nboundtaskslvl(ilevel_index)

        ASSOCIATE( &
            p => dp_f%arr, &
            pbuffer => dp_f%buffers, &
            bp => bp_f%arr, &
            dx => dx_f%arr, &
            dy => dy_f%arr, &
            dz => dz_f%arr, &
            ddx => ddx_f%arr, &
            ddy => ddy_f%arr, &
            ddz => ddz_f%arr)

        CALL profile_range_push("boundpressurebp")
        !$omp target teams distribute private(i, igrid, iface, kk, jj, ii, &
        !$omp& ip3, ipx, ipy, ipz, ipbb)
        DO i = 1, nboundtasks
            igrid = boundtasks(1, i, ilevel_index)
            iface = boundtasks(2, i, ilevel_index)

            CALL get_mgdims(kk, jj, ii, igrid)
            CALL get_ibb(ipbb, iface, igrid)
            CALL get_ip3(ip3, igrid)
            CALL get_ipx(ipx, igrid)
            CALL get_ipy(ipy, igrid)
            CALL get_ipz(ipz, igrid)

            SELECT CASE (iface)
            CASE (1)
                CALL bfront_bp(kk, jj, ii, 2, 3, 4, 2, &
                    pbuffer(ipbb), p(ip3), bp(ip3), &
                    ddx(ipx), ddy(ipy), ddz(ipz), dx(ipx))
            CASE (2)
                CALL bfront_bp(kk, jj, ii, ii-1, ii-2, ii-3, ii-2, &
                    pbuffer(ipbb), p(ip3), bp(ip3), &
                    ddx(ipx), ddy(ipy), ddz(ipz), dx(ipx))
            CASE (3)
                CALL bright_bp(kk, jj, ii, 2, 3, 4, 2, &
                    pbuffer(ipbb), p(ip3), bp(ip3), &
                    ddx(ipx), ddy(ipy), ddz(ipz), dy(ipy))
            CASE (4)
                CALL bright_bp(kk, jj, ii, jj-1, jj-2, jj-3, jj-2, &
                    pbuffer(ipbb), p(ip3), bp(ip3), &
                    ddx(ipx), ddy(ipy), ddz(ipz), dy(ipy))
            CASE (5)
                CALL bbottom_bp(kk, jj, ii, 2, 3, 4, 2, &
                    pbuffer(ipbb), p(ip3), bp(ip3), &
                    ddx(ipx), ddy(ipy), ddz(ipz), dz(ipz))
            CASE (6)
                CALL bbottom_bp(kk, jj, ii, kk-1, kk-2, kk-3, kk-2, &
                    pbuffer(ipbb), p(ip3), bp(ip3), &
                    ddx(ipx), ddy(ipy), ddz(ipz), dz(ipz))
            END SELECT
        END DO
        !$omp end target teams distribute
        CALL profile_range_pop()

        end associate
    END SUBROUTINE bound_pressure_impl_bp


    SUBROUTINE bound_pressure_impl_nobp(ilevel, dp_f, dx_f, dy_f, dz_f, &
        ddx_f, ddy_f, ddz_f)
        ! Subroutine arguments
        INTEGER(intk), INTENT(in) :: ilevel
        TYPE(field_t), INTENT(inout) :: dp_f
        TYPE(field_t), POINTER, INTENT(in) :: dx_f, dy_f, dz_f
        TYPE(field_t), POINTER, INTENT(in) :: ddx_f, ddy_f, ddz_f

        ! Local variables
        INTEGER(intk) :: nboundtasks, ilevel_index, i, igrid, iface
        INTEGER(intk) :: kk, jj, ii, ip3, ipx, ipy, ipz, ipbb
        REAL(realk), ALLOCATABLE, DIMENSION(:) :: p, pbuffer, bp, dx, dy, dz, &
            ddx, ddy, ddz

        CALL level_index(ilevel_index, ilevel)
        nboundtasks = nboundtaskslvl(ilevel_index)

        ASSOCIATE( &
            p => dp_f%arr, &
            pbuffer => dp_f%buffers, &
            dx => dx_f%arr, &
            dy => dy_f%arr, &
            dz => dz_f%arr, &
            ddx => ddx_f%arr, &
            ddy => ddy_f%arr, &
            ddz => ddz_f%arr)

        !$omp target teams distribute private(i, igrid, iface, kk, jj, ii, &
        !$omp& ip3, ipx, ipy, ipz, ipbb)
        DO i = 1, nboundtasks
            igrid = boundtasks(1, i, ilevel_index)
            iface = boundtasks(2, i, ilevel_index)

            CALL get_mgdims(kk, jj, ii, igrid)

            CALL get_mgdims(kk, jj, ii, igrid)
            CALL get_ibb(ipbb, iface, igrid)
            CALL get_ip3(ip3, igrid)
            CALL get_ipx(ipx, igrid)
            CALL get_ipy(ipy, igrid)
            CALL get_ipz(ipz, igrid)

            SELECT CASE (iface)
            CASE (1)
                CALL bfront(kk, jj, ii, 2, 3, 4, 2, &
                    pbuffer(ipbb), p(ip3), &
                    ddx(ipx), ddy(ipy), ddz(ipz), dx(ipx))
            CASE (2)
                CALL bfront(kk, jj, ii, ii-1, ii-2, ii-3, ii-2, &
                    pbuffer(ipbb), p(ip3), &
                    ddx(ipx), ddy(ipy), ddz(ipz), dx(ipx))
            CASE (3)
                CALL bright(kk, jj, ii, 2, 3, 4, 2, &
                    pbuffer(ipbb), p(ip3), &
                    ddx(ipx), ddy(ipy), ddz(ipz), dy(ipy))
            CASE (4)
                CALL bright(kk, jj, ii, jj-1, jj-2, jj-3, jj-2, &
                    pbuffer(ipbb), p(ip3), &
                    ddx(ipx), ddy(ipy), ddz(ipz), dy(ipy))
            CASE (5)
                CALL bbottom(kk, jj, ii, 2, 3, 4, 2, &
                    pbuffer(ipbb), p(ip3), &
                    ddx(ipx), ddy(ipy), ddz(ipz), dz(ipz))
            CASE (6)
                CALL bbottom(kk, jj, ii, kk-1, kk-2, kk-3, kk-2, &
                    pbuffer(ipbb), p(ip3), &
                    ddx(ipx), ddy(ipy), ddz(ipz), dz(ipz))
            END SELECT
        END DO
        !$omp end target teams distribute

        END ASSOCIATE
    END SUBROUTINE bound_pressure_impl_nobp


    SUBROUTINE bfront(kk, jj, ii, i2, i3, i4, istag2, pbuffer, &
        p, ddx, ddy, ddz, dx)
        !$omp declare target
        ! Subroutine arguments
        INTEGER(intk), INTENT(in) :: kk, jj, ii, i2, i3, i4, istag2
        REAL(realk), INTENT(in) :: pbuffer(kk, jj)
        REAL(realk), INTENT(inout) :: p(kk, jj, ii)
        REAL(realk), INTENT(in) :: ddx(ii), ddy(jj), ddz(kk)
        REAL(realk), INTENT(in) :: dx(ii)

        ! Local variables
        INTEGER(intk) :: k, j, i, m, n
        REAL(realk) :: pcnew, bpc

        i = MIN(i3, i4)

        !$omp parallel do collapse(2) private(j, k, pcnew, bpc, m, n)
        DO j = 3, jj-2, 2
            DO k = 3, kk-2, 2
                CALL pressureftocone(k, j, i, kk, jj, ii, p, &
                    ddx, ddy, ddz, pcnew, bpc)
                DO m = 0, 1
                    DO n = 0, 1
                        p(k+n, j+m, i2) = p(k+n, j+m, i3) &
                            + dx(istag2)/(ddx(i3)+ddx(i2)) &
                            *(pbuffer(k, j) - pcnew)
                    END DO
                END DO
            END DO
        END DO
        !$omp end parallel do
    END SUBROUTINE bfront


    SUBROUTINE bfront_bp(kk, jj, ii, i2, i3, i4, istag2, pbuffer, &
        p, bp, ddx, ddy, ddz, dx)
        !$omp declare target
        ! Subroutine arguments
        INTEGER(intk), INTENT(in) :: kk, jj, ii, i2, i3, i4, istag2
        REAL(realk), INTENT(in) :: pbuffer(kk, jj)
        REAL(realk), INTENT(inout) :: p(kk, jj, ii)
        REAL(realk), INTENT(in) :: bp(kk, jj, ii)
        REAL(realk), INTENT(in) :: ddx(ii), ddy(jj), ddz(kk)
        REAL(realk), INTENT(in) :: dx(ii)

        ! Local variables
        INTEGER(intk) :: k, j, i, m, n
        REAL(realk) :: pcnew, bpc, sb11, sb12, sb13, sb14, fak

        i = MIN(i3, i4)

        !$omp parallel do collapse(2) private(j, k, pcnew, bpc, sb11, sb12, &
        !$omp& sb13, sb14, fak, m, n)
        DO j = 3, jj-2, 2
            DO k = 3, kk-2, 2
                CALL pressureftocone(k, j, i, kk, jj, ii, p, bp, ddx, ddy, &
                    ddz, pcnew, bpc)
                IF (bpc < 0.5) pcnew = pbuffer(k, j)

                sb11 = bp(k, j, i2)*bp(k, j, i3)
                sb12 = bp(k, j+1, i2)*bp(k, j+1, i3)
                sb13 = bp(k+1, j, i2)*bp(k+1, j, i3)
                sb14 = bp(k+1, j+1, i2)*bp(k+1, j+1, i3)

                fak = (sb11*ddy(j)*ddz(k) + sb12*ddy(j+1)*ddz(k) &
                    + sb13*ddy(j)*ddz(k+1) + sb14*ddy(j+1)*ddz(k+1)) &
                    /((ddy(j)+ddy(j+1))*(ddz(k)+ddz(k+1)))
                IF (fak < 0.1) fak = 1.0
                fak = 1.0/fak

                DO m = 0, 1
                    DO n = 0, 1
                        p(k+n, j+m, i2) = p(k+n, j+m, i3) &
                            + fak*dx(istag2)/(ddx(i3)+ddx(i2)) &
                            *(pbuffer(k, j) - pcnew)
                    END DO
                END DO
            END DO
        END DO
        !$omp end parallel do
    END SUBROUTINE bfront_bp


    SUBROUTINE bright(kk, jj, ii, j2, j3, j4, jstag2, pbuffer, &
        p, ddx, ddy, ddz, dy)
        !$omp declare target
        ! Subroutine arguments
        INTEGER(intk), INTENT(in) :: kk, jj, ii, j2, j3, j4, jstag2
        REAL(realk), INTENT(in) :: pbuffer(kk, ii)
        REAL(realk), INTENT(inout) :: p(kk, jj, ii)
        REAL(realk), INTENT(in) :: ddx(ii), ddy(jj), ddz(kk)
        REAL(realk), INTENT(in) :: dy(jj)

        ! Local variables
        INTEGER(intk) :: k, j, i, l, n
        REAL(realk) :: pcnew, bpc

        j = MIN(j3, j4)

        !$omp parallel do collapse(2) private(i, k, pcnew, bpc, l, n)
        DO i = 3, ii-2, 2
            DO k = 3, kk-2, 2
                CALL pressureftocone(k, j, i, kk, jj, ii, p, &
                    ddx, ddy, ddz, pcnew, bpc)
                DO l = 0, 1
                    DO n = 0, 1
                        p(k+n, j2, i+l) = p(k+n, j3, i+l) &
                            + dy(jstag2)/(ddy(j3)+ddy(j2)) &
                            *(pbuffer(k, i) - pcnew)
                    END DO
                END DO
            END DO
        END DO
        !$omp end parallel do
    END SUBROUTINE bright


    SUBROUTINE bright_bp(kk, jj, ii, j2, j3, j4, jstag2, pbuffer, &
        p, bp, ddx, ddy, ddz, dy)
        !$omp declare target
        ! Subroutine arguments
        INTEGER(intk), INTENT(in) :: kk, jj, ii, j2, j3, j4, jstag2
        REAL(realk), INTENT(in) :: pbuffer(kk, ii)
        REAL(realk), INTENT(inout) :: p(kk, jj, ii)
        REAL(realk), INTENT(in) :: bp(kk, jj, ii)
        REAL(realk), INTENT(in) :: ddx(ii), ddy(jj), ddz(kk)
        REAL(realk), INTENT(in) :: dy(jj)

        ! Local variables
        INTEGER(intk) :: k, j, i, l, n
        REAL(realk) :: pcnew, bpc, sb11, sb12, sb13, sb14, fak

        j = MIN(j3, j4)

        !$omp parallel do collapse(2) private(i, k, pcnew, bpc, sb11, sb12, &
        !$omp& sb13, sb14, fak, l, n)
        DO i = 3, ii-2, 2
            DO k = 3, kk-2, 2
                CALL pressureftocone(k, j, i, kk, jj, ii, p, bp, &
                    ddx, ddy, ddz, pcnew, bpc)
                IF (bpc < 0.5) pcnew = pbuffer(k, i)

                sb11 = bp(k, j2, i)*bp(k, j3, i)
                sb12 = bp(k, j2, i+1)*bp(k, j3, i+1)
                sb13 = bp(k+1, j2, i)*bp(k+1, j3, i)
                sb14 = bp(k+1, j2, i+1)*bp(k+1, j3, i+1)

                fak = (sb11*ddx(i)*ddz(k) + sb12*ddx(i+1)*ddz(k) &
                    + sb13*ddx(i)*ddz(k+1) + sb14*ddx(i+1)*ddz(k+1)) &
                    /((ddx(i)+ddx(i+1))*(ddz(k)+ddz(k+1)))
                IF (fak < 0.1) fak = 1.0
                fak = 1.0/fak

                DO l = 0, 1
                    DO n = 0, 1
                        p(k+n, j2, i+l) = p(k+n, j3, i+l) &
                            + fak*dy(jstag2)/(ddy(j3)+ddy(j2)) &
                            *(pbuffer(k, i) - pcnew)
                    END DO
                END DO
            END DO
        END DO
        !$omp end parallel do
    END SUBROUTINE bright_bp


    SUBROUTINE bbottom(kk, jj, ii, k2, k3, k4, kstag2, pbuffer, &
        p, ddx, ddy, ddz, dz)
        !$omp declare target
        ! Subroutine arguments
        INTEGER(intk), INTENT(in) :: kk, jj, ii, k2, k3, k4, kstag2
        REAL(realk), INTENT(in) :: pbuffer(jj, ii)
        REAL(realk), INTENT(inout) :: p(kk, jj, ii)
        REAL(realk), INTENT(in) :: ddx(ii), ddy(jj), ddz(kk)
        REAL(realk), INTENT(in) :: dz(kk)

        ! Local variables
        INTEGER(intk) :: k, j, i, l, m
        REAL(realk) :: pcnew, bpc

        k = MIN(k3, k4)

        !$omp parallel do collapse(2) private(j, i, pcnew, bpc, l, m)
        DO i = 3, ii-2, 2
            DO j = 3, jj-2, 2
                CALL pressureftocone(k, j, i, kk, jj, ii, p, &
                    ddx, ddy, ddz, pcnew, bpc)
                DO l = 0, 1
                    DO m = 0, 1
                        p(k2, j+m, i+l) = p(k3, j+m, i+l) &
                            + dz(kstag2)/(ddz(k3)+ddz(k2)) &
                            *(pbuffer(j, i) - pcnew)
                    END DO
                END DO
            END DO
        END DO
        !$omp end parallel do
    END SUBROUTINE bbottom


    SUBROUTINE bbottom_bp(kk, jj, ii, k2, k3, k4, kstag2, pbuffer, &
        p, bp, ddx, ddy, ddz, dz)
        !$omp declare target
        ! Subroutine arguments
        INTEGER(intk), INTENT(in) :: kk, jj, ii, k2, k3, k4, kstag2
        REAL(realk), INTENT(in) :: pbuffer(jj, ii)
        REAL(realk), INTENT(inout) :: p(kk, jj, ii)
        REAL(realk), INTENT(in) :: bp(kk, jj, ii)
        REAL(realk), INTENT(in) :: ddx(ii), ddy(jj), ddz(kk)
        REAL(realk), INTENT(in) :: dz(kk)

        ! Local variables
        INTEGER(intk) :: k, j, i, l, m
        REAL(realk) :: pcnew, bpc, sb11, sb12, sb13, sb14, fak

        k = MIN(k3, k4)

        !$omp parallel do collapse(2) private(j, i, pcnew, bpc, sb11, sb12, &
        !$omp& sb13, sb14, fak, l, m)
        DO i = 3, ii-2, 2
            DO j = 3, jj-2, 2
                CALL pressureftocone(k, j, i, kk, jj, ii, p, bp, &
                    ddx, ddy, ddz, pcnew, bpc)
                IF (bpc < 0.5) pcnew = pbuffer(j, i)

                sb11 = bp(k2, j, i)*bp(k3, j, i)
                sb12 = bp(k2, j, i+1)*bp(k3, j, i+1)
                sb13 = bp(k2, j+1, i)*bp(k3, j+1, i)
                sb14 = bp(k2, j+1, i+1)*bp(k3, j+1, i+1)

                fak = (sb11*ddx(i)*ddy(j) + sb12*ddx(i+1)*ddy(j) &
                    + sb13*ddx(i)*ddy(j+1) + sb14*ddx(i+1)*ddy(j+1)) &
                    /((ddx(i)+ddx(i+1))*(ddy(j)+ddy(j+1)))
                IF (fak < 0.1) fak = 1.0
                fak = 1.0/fak

                DO l = 0, 1
                    DO m = 0, 1
                        p(k2, j+m, i+l) = p(k3, j+m, i+l) &
                            + fak*dz(kstag2)/(ddz(k3)+ddz(k2)) &
                            *(pbuffer(j, i) - pcnew)
                    END DO
                END DO
            END DO
        END DO
        !$omp end parallel do
    END SUBROUTINE bbottom_bp


    PURE SUBROUTINE pressureftocone_A(k, j, i, kk, jj, ii, p, ddx, ddy, &
            ddz, pc, bpc)
        !$omp declare target
        ! Subroutine arguments
        INTEGER(intk), INTENT(in) :: k, j, i
        INTEGER(intk), INTENT(in) :: kk, jj, ii
        REAL(realk), INTENT(in) :: p(kk, jj, ii)
        REAL(realk), INTENT(in) :: ddx(ii), ddy(jj), ddz(kk)
        REAL(realk), INTENT(out) :: pc, bpc

        ! Local variables
        INTEGER(intk) :: l, m, n
        REAL(realk) :: vol, sump, sumvol

        bpc = 1.0

        sump = 0.0
        sumvol = 0.0
        DO l = 0, 1
            DO m = 0, 1
                DO n = 0, 1
                    vol = ddz(k+n)*ddy(j+m)*ddx(i+l)
                    sump = sump + p(k+n, j+m, i+l)*vol
                    sumvol = sumvol + vol
                END DO
            END DO
        END DO
        pc = sump/sumvol
    END SUBROUTINE pressureftocone_A


    PURE SUBROUTINE pressureftocone_B(k, j, i, kk, jj, ii, p, bp, ddx, ddy, &
            ddz, pc, bpc)
        !$omp declare target
        ! Subroutine arguments
        INTEGER(intk), INTENT(in) :: k, j, i
        INTEGER(intk), INTENT(in) :: kk, jj, ii
        REAL(realk), INTENT(in) :: p(kk, jj, ii)
        REAL(realk), INTENT(in) :: bp(kk, jj, ii)
        REAL(realk), INTENT(in) :: ddx(ii), ddy(jj), ddz(kk)
        REAL(realk), INTENT(out) :: pc, bpc

        ! Local variables
        INTEGER(intk) :: l, m, n
        REAL(realk) :: vol, sumbp, sump, sumvol

        sumbp = 0.0
        DO l = 0, 1
            DO m = 0, 1
                DO n = 0, 1
                    sumbp = sumbp + bp(k+n, j+m, i+l)
                END DO
            END DO
        END DO
        bpc = MIN(sumbp, 1.0_realk)

        IF (bpc < 0.5) THEN
            pc = 0.0
        ELSE
            sump = 0.0
            sumvol = 0.0
            DO l = 0, 1
                DO m = 0, 1
                    DO n = 0, 1
                        vol = bp(k+n, j+m, i+l)*ddz(k+n)*ddy(j+m)*ddx(i+l)
                        sump = sump + p(k+n, j+m, i+l)*vol
                        sumvol = sumvol + vol
                    END DO
                END DO
            END DO
            pc = sump/sumvol
        END IF
    END SUBROUTINE pressureftocone_B


    PURE SUBROUTINE level_index(ilevel_index, ilevel)
        ! Subroutine arguments
        INTEGER(intk), INTENT(out) :: ilevel_index
        INTEGER(intk), INTENT(in) :: ilevel

        ! Local variables
        ! none...

        ilevel_index = ilevel - minlevel + 1
    END SUBROUTINE level_index


    SUBROUTINE finish_bound_pressure()
        ! Subroutine arguments
        ! none...

        ! Local variables
        ! none...

        !$omp target exit data map(always, delete: boundtasks)
        DEALLOCATE(boundtasks)
        DEALLOCATE(nboundtaskslvl)
    END SUBROUTINE finish_bound_pressure
END MODULE bound_pressure_mod
