MODULE par_ftoc_mod
    USE MPI_f08
    USE core_mod

    IMPLICIT NONE (type, external)
    PRIVATE

    ! The information in the first dimension is sorted as follows:
    !   Field 1: Start index in i-direction
    !   Field 2: Stop index in i-direction
    !   Field 3: Start index in j-direction
    !   Field 4: Stop index in j-direction
    !   Field 5: Start index in k-direction
    !   Field 6: Stop index in k-direction
    !   Field 7: ID of grid
    !   Field 8: Level of grid
    !   Field 9: Offset in the send buffer to pack to
    INTEGER(intk), ALLOCATABLE, DIMENSION(:, :) :: packops_v1, packops_v2, &
        packops_v3
    !$omp declare target(packops_v1, packops_v2, packops_v3)

    ! The information in the first dimension is sorted as follows:
    !   Field 1: Start index in i-direction
    !   Field 2: Stop index in i-direction
    !   Field 3: Start index in j-direction
    !   Field 4: Stop index in j-direction
    !   Field 5: Start index in k-direction
    !   Field 6: Stop index in k-direction
    !   Field 7: ID of coarse grid
    !   Field 8: Level of fine grid
    !   Field 9: Offset in the recv buffer to unpack from
    INTEGER(intk), ALLOCATABLE, DIMENSION(:, :) :: unpackops_v1, unpackops_v2, &
        unpackops_v3
    !$omp declare target(unpackops_v1, unpackops_v2, unpackops_v3)

    ! The information in the first dimension is sorted as follows:
    !   Field 1: Rank of sending process
    !   Field 2: Rank of receiving process
    !   Field 3: ID of sending grid
    !   Field 4: ID of receiving grid
    !   Field 5: Which face (1..6) to send
    !   Field 6: Which face (1..6) to receive, -1 for internal
    INTEGER(intk), ALLOCATABLE :: sendconns(:, :), recvconns(:, :)

    ! Lists that hold the send and receive request arrays
    TYPE(MPI_Request), ALLOCATABLE :: sendreqs(:), recvreqs(:)

    ! Number of messages that are ACTUALLY sent and received
    INTEGER(intk) :: nsend, nrecv

    ! Number of send and receive connections
    INTEGER(intk) :: isend = 0, irecv = 0

    ! Variable to indicate if the connection information has
    ! been created.
    LOGICAL :: is_init = .FALSE.

    PUBLIC :: par_ftoc_norm, init_par_ftoc, finish_par_ftoc

CONTAINS

    ! Main parent function
    SUBROUTINE par_ftoc_norm(ilevel, v1, v2, v3, sum)

        ! Subroutine arguments
        INTEGER(intk), INTENT(in) :: ilevel
        TYPE(field_t), INTENT(inout) :: v1, v2, v3
        LOGICAL, OPTIONAL, INTENT(in) :: sum

        ! Local variables
        LOGICAL :: sumflag

        CALL start_timer(211)

        ! Check if the connection information has been created
        IF (.NOT. is_init) THEN
            CALL errr(__FILE__, __LINE__)
        END IF

        IF (PRESENT(sum)) THEN
            sumflag = sum
        ELSE
            sumflag = .FALSE.
        END IF

        CALL recv_all(ilevel)

        CALL pack_v1(ilevel, v1, sumflag)
        CALL pack_v2(ilevel, v2, sumflag)
        CALL pack_v3(ilevel, v3, sumflag)

        CALL send_all(ilevel)
        CALL process_bufs(ilevel, v1, v2, v3)

        nrecv = 0
        nsend = 0

        CALL stop_timer(211)
    END SUBROUTINE par_ftoc_norm


    SUBROUTINE precompute_packops()
        ! Local variables
        INTEGER(intk) :: i, igrid, iface, buf_offset
        INTEGER(intk) :: n1, n2, n3, i1, i2, i3
        INTEGER(intk) :: target_level, maxlev

        ! Count packops per component
        n1 = 0
        n2 = 0
        n3 = 0
        DO i = 1, isend
            iface = sendconns(5, i)
            SELECT CASE (iface)
            CASE (1, 2)
                n1 = n1 + 1
            CASE (3, 4)
                n2 = n2 + 1
            CASE (5, 6)
                n3 = n3 + 1
            END SELECT
        END DO

        ALLOCATE(packops_v1(9, n1))
        ALLOCATE(packops_v2(9, n2))
        ALLOCATE(packops_v3(9, n3))

        ! Find max level
        maxlev = 0
        DO i = 1, isend
            maxlev = MAX(maxlev, level(sendconns(3, i)))
        END DO

        ! Fill packops per component
        i1 = 0
        i2 = 0
        i3 = 0
        DO target_level = 1, maxlev
            buf_offset = 0
            DO i = 1, isend
                igrid = sendconns(3, i)
                IF (level(igrid) /= target_level) CYCLE
                iface = sendconns(5, i)

                SELECT CASE (iface)
                CASE (1, 2)
                    i1 = i1 + 1
                    CALL set_packop(packops_v1(:, i1), igrid, iface, &
                        buf_offset)
                CASE (3, 4)
                    i2 = i2 + 1
                    CALL set_packop(packops_v2(:, i2), igrid, iface, &
                        buf_offset)
                CASE (5, 6)
                    i3 = i3 + 1
                    CALL set_packop(packops_v3(:, i3), igrid, iface, &
                        buf_offset)
                END SELECT

                buf_offset = buf_offset + face_area(igrid, iface)
            END DO
        END DO

        !$omp target enter data map(always, to: packops_v1, packops_v2, &
        !$omp& packops_v3)
    END SUBROUTINE precompute_packops


    SUBROUTINE set_packop(packop, igrid, iface, buf_offset)
        ! Subroutine arguments
        INTEGER(intk), INTENT(out) :: packop(9)
        INTEGER(intk), INTENT(in) :: igrid, iface, buf_offset

        ! Local variables
        INTEGER(intk) :: kstart, kstop, jstart, jstop, istart, istop

        CALL start_and_stop(kstart, kstop, jstart, jstop, istart, istop, &
            igrid, iface)

        packop(1) = istart
        packop(2) = istop
        packop(3) = jstart
        packop(4) = jstop
        packop(5) = kstart
        packop(6) = kstop
        packop(7) = igrid
        packop(8) = level(igrid)
        packop(9) = buf_offset
    END SUBROUTINE set_packop


    SUBROUTINE precompute_unpackops()
        ! Local variables
        INTEGER(intk) :: i, igridf, igridc, iface, ifacerecv, buf_offset
        INTEGER(intk) :: n1, n2, n3, i1, i2, i3
        INTEGER(intk) :: target_level, maxlev

        ! Count unpacks per component
        n1 = 0
        n2 = 0
        n3 = 0
        DO i = 1, irecv
            iface = recvconns(5, i)
            SELECT CASE (iface)
            CASE (1, 2)
                n1 = n1 + 1
            CASE (3, 4)
                n2 = n2 + 1
            CASE (5, 6)
                n3 = n3 + 1
            END SELECT
        END DO

        ALLOCATE(unpackops_v1(9, n1))
        ALLOCATE(unpackops_v2(9, n2))
        ALLOCATE(unpackops_v3(9, n3))

        ! Find max level
        maxlev = 0
        DO i = 1, irecv
            maxlev = MAX(maxlev, level(recvconns(3, i)))
        END DO

        ! Fill unpackops per component, with per-level dense buf_offset
        i1 = 0
        i2 = 0
        i3 = 0
        DO target_level = 1, maxlev
            buf_offset = 0
            DO i = 1, irecv
                igridf = recvconns(3, i)
                IF (level(igridf) /= target_level) CYCLE
                igridc = recvconns(4, i)
                iface = recvconns(5, i)
                ifacerecv = recvconns(6, i)

                SELECT CASE (iface)
                CASE (1, 2)
                    i1 = i1 + 1
                    CALL set_unpackop(unpackops_v1(:, i1), igridf, igridc, &
                        iface, ifacerecv, buf_offset)
                CASE (3, 4)
                    i2 = i2 + 1
                    CALL set_unpackop(unpackops_v2(:, i2), igridf, igridc, &
                        iface, ifacerecv, buf_offset)
                CASE (5, 6)
                    i3 = i3 + 1
                    CALL set_unpackop(unpackops_v3(:, i3), igridf, igridc, &
                        iface, ifacerecv, buf_offset)
                END SELECT

                buf_offset = buf_offset + face_area(igridf, iface)
            END DO
        END DO

        !$omp target enter data map(always, to: unpackops_v1, unpackops_v2, &
        !$omp& unpackops_v3)
    END SUBROUTINE precompute_unpackops


    SUBROUTINE set_unpackop(unpackop, igridf, igridc, iface, ifacerecv, &
            buf_offset)
        ! Subroutine arguments
        INTEGER(intk), INTENT(out) :: unpackop(9)
        INTEGER(intk), INTENT(in) :: igridf, igridc, iface, ifacerecv
        INTEGER(intk), INTENT(in) :: buf_offset

        ! Local variables
        INTEGER(intk) :: ipos, jpos, kpos
        INTEGER(intk) :: kkc, jjc, iic
        INTEGER(intk) :: kkf, jjf, iif
        INTEGER(intk) :: ista, isto, jsta, jsto, ksta, ksto

        CALL get_mgdims(kkf, jjf, iif, igridf)
        ipos = iposition(igridf)  ! Position of fine grid within coarse grid
        jpos = jposition(igridf)
        kpos = kposition(igridf)

        ! Select the entire fine grid in the coarse grid
        ista = ipos
        isto = (ipos + (iif-4)/2 - 1)
        jsta = jpos
        jsto = (jpos + (jjf-4)/2 - 1)
        ksta = kpos
        ksto = (kpos + (kkf-4)/2 - 1)

        IF (ifacerecv < 0) THEN
            ! Internal to the grid. We need to inspect the position of the
            ! fine grid within the coarse grid to determine the start and
            ! stop indices.
            SELECT CASE (iface)  ! The sending face determines position
            CASE (1)
                ista = ipos-1
                isto = ista
            CASE (2)
                ista = isto
                ! isto unchanged
            CASE (3)
                jsta = jpos-1
                jsto = jsta
            CASE (4)
                jsta = jsto
                ! jsto unchanged
            CASE (5)
                ksta = kpos-1
                ksto = ksta
            CASE (6)
                ksta = ksto
                ! ksto unchanged
            END SELECT
        ELSE
            ! The PAR on the fine-grid lies on top of an external face of this
            ! grid (it has to be a CON - everything else would be an error).
            CALL get_mgdims(kkc, jjc, iic, igridc)   ! Dimensions of coarse grid
            SELECT CASE (ifacerecv)  ! The receiving face determines position
            CASE (1)
                ista = 2
                isto = 2
            CASE (2)
                ista = iic-2
                isto = iic-2
            CASE (3)
                jsta = 2
                jsto = 2
            CASE (4)
                jsta = jjc-2
                jsto = jjc-2
            CASE (5)
                ksta = 2
                ksto = 2
            CASE (6)
                ksta = kkc-2
                ksto = kkc-2
            END SELECT
        END IF

        unpackop(1) = ista
        unpackop(2) = isto
        unpackop(3) = jsta
        unpackop(4) = jsto
        unpackop(5) = ksta
        unpackop(6) = ksto
        unpackop(7) = igridc
        unpackop(8) = level(igridf)
        unpackop(9) = buf_offset
    END SUBROUTINE set_unpackop


    SUBROUTINE pack_v1(ilevel, v1, sum)
        ! Subroutine arguments
        INTEGER(intk), INTENT(in) :: ilevel
        TYPE(field_t), INTENT(inout) :: v1
        LOGICAL, INTENT(in) :: sum

        ! Local variables
        INTEGER(intk) :: ipack, npack, i, j, k, idx
        INTEGER(intk) :: jstride, kcount
        INTEGER(intk) :: ista, jsta, ksta, iend, jend, kend
        INTEGER(intk) :: igrid, lvl, buf_offset
        REAL(realk) :: sum_ua, sum_a
        REAL(realk), POINTER, CONTIGUOUS :: ff(:, :, :)
        REAL(realk), POINTER, CONTIGUOUS :: ddy(:), ddz(:)
        TYPE(field_t), POINTER :: ddy_f, ddz_f

        CALL get_field(ddy_f, "DDY")
        CALL get_field(ddz_f, "DDZ")

        npack = SIZE(packops_v1, dim=2)

        !$omp target teams distribute private(ff, ddy, ddz, i, j, k, kcount, &
        !$omp& igrid, lvl, ista, jsta, ksta, iend, jend, kend, buf_offset, &
        !$omp& jstride, sum_ua, sum_a, idx)
        DO ipack = 1, npack
            CALL get_packop(ista, iend, jsta, jend, ksta, kend, &
                igrid, lvl, buf_offset, packops_v1(:, ipack))

            IF (lvl /= ilevel) CYCLE

            CALL get_grid3_real(ff, v1, igrid)
            CALL get_grid1_real(ddy, ddy_f, igrid)
            CALL get_grid1_real(ddz, ddz_f, igrid)

            i = ista
            kcount = (kend - ksta) / 2 + 1

            IF (sum) THEN
                !$omp parallel do collapse(2) private(jstride, idx)
                DO j = jsta, jend, 2
                    DO k = ksta, kend, 2
                        ! Recompute to collapse
                        jstride = ((j - jsta) / 2) * kcount
                        idx = buf_offset + jstride + (k - ksta) / 2 + 1
                        sendbuf(idx) = &
                            ff(k, j, i) + ff(k, j+1, i) &
                            + ff(k+1, j, i) + ff(k+1, j+1, i)
                    END DO
                END DO
                !$omp end parallel do
            ELSE
                !$omp parallel do collapse(2) private(jstride, sum_ua, sum_a, &
                !$omp& idx)
                DO j = jsta, jend, 2
                    DO k = ksta, kend, 2
                        ! Recompute to collapse
                        jstride = ((j - jsta) / 2) * kcount
                        sum_ua = ff(k, j, i) * ddy(j) * ddz(k) &
                            + ff(k, j+1, i) * ddy(j+1) * ddz(k) &
                            + ff(k+1, j, i) * ddy(j) * ddz(k+1) &
                            + ff(k+1, j+1, i) * ddy(j+1) * ddz(k+1)
                        sum_a  = (ddy(j) + ddy(j+1)) * (ddz(k) + ddz(k+1))

                        idx = buf_offset + jstride + (k - ksta)/2 + 1
                        sendbuf(idx) = sum_ua / sum_a
                    END DO
                END DO
                !$omp end parallel do
            END IF
        END DO
        !$omp end target teams distribute
    END SUBROUTINE pack_v1


    SUBROUTINE pack_v2(ilevel, v2, sum)
        ! Subroutine arguments
        INTEGER(intk), INTENT(in) :: ilevel
        TYPE(field_t), INTENT(inout) :: v2
        LOGICAL, INTENT(in) :: sum

        ! Local variables
        INTEGER(intk) :: ipack, npack, i, j, k, idx
        INTEGER(intk) :: istride, kcount
        INTEGER(intk) :: ista, jsta, ksta, iend, jend, kend
        INTEGER(intk) :: igrid, lvl, buf_offset
        REAL(realk) :: sum_ua, sum_a
        REAL(realk), POINTER, CONTIGUOUS :: ff(:, :, :)
        REAL(realk), POINTER, CONTIGUOUS :: ddx(:), ddz(:)
        TYPE(field_t), POINTER :: ddx_f, ddz_f

        CALL get_field(ddx_f, "DDX")
        CALL get_field(ddz_f, "DDZ")

        npack = SIZE(packops_v2, dim=2)

        !$omp target teams distribute private(ff, ddx, ddz, i, j, k, kcount, &
        !$omp& igrid, lvl, ista, jsta, ksta, iend, jend, kend, buf_offset, &
        !$omp& istride, sum_ua, sum_a, idx)
        DO ipack = 1, npack
            CALL get_packop(ista, iend, jsta, jend, ksta, kend, &
                igrid, lvl, buf_offset, packops_v2(:, ipack))

            IF (lvl /= ilevel) CYCLE

            CALL get_grid3_real(ff, v2, igrid)
            CALL get_grid1_real(ddx, ddx_f, igrid)
            CALL get_grid1_real(ddz, ddz_f, igrid)

            j = jsta
            kcount = (kend - ksta) / 2 + 1

            IF (sum) THEN
                !$omp parallel do collapse(2) private(istride, idx)
                DO i = ista, iend, 2
                    DO k = ksta, kend, 2
                        ! Recompute to collapse later
                        istride = ((i - ista) / 2) * kcount
                        idx = buf_offset + istride + (k - ksta) / 2 + 1
                        sendbuf(idx) = &
                            ff(k, j, i) + ff(k, j, i+1) &
                            + ff(k+1, j, i) + ff(k+1, j, i+1)
                    END DO
                END DO
                !$omp end parallel do
            ELSE
                !$omp parallel do collapse(2) private(istride, sum_ua, sum_a, &
                !$omp& idx)
                DO i = ista, iend, 2
                    DO k = ksta, kend, 2
                        ! Recompute to collapse later
                        istride = ((i - ista) / 2) * kcount
                        sum_ua = ff(k, j, i) * ddx(i) * ddz(k) &
                            + ff(k, j, i+1) * ddx(i+1) * ddz(k) &
                            + ff(k+1, j, i) * ddx(i) * ddz(k+1) &
                            + ff(k+1, j, i+1) * ddx(i+1) * ddz(k+1)
                        sum_a = (ddx(i) + ddx(i+1)) * (ddz(k) + ddz(k+1))

                        idx = buf_offset + istride + (k - ksta) / 2 + 1
                        sendbuf(idx) = sum_ua / sum_a
                    END DO
                END DO
                !$omp end parallel do
            END IF
        END DO
        !$omp end target teams distribute
    END SUBROUTINE pack_v2


    SUBROUTINE pack_v3(ilevel, v3, sum)
        ! Subroutine arguments
        INTEGER(intk), INTENT(in) :: ilevel
        TYPE(field_t), INTENT(inout) :: v3
        LOGICAL, INTENT(in) :: sum

        ! Local variables
        INTEGER(intk) :: ipack, npack, i, j, k, idx
        INTEGER(intk) :: istride, jcount
        INTEGER(intk) :: ista, jsta, ksta, iend, jend, kend
        INTEGER(intk) :: igrid, lvl, buf_offset
        REAL(realk) :: sum_ua, sum_a
        REAL(realk), POINTER, CONTIGUOUS :: ff(:, :, :)
        REAL(realk), POINTER, CONTIGUOUS :: ddx(:), ddy(:)
        TYPE(field_t), POINTER :: ddx_f, ddy_f

        CALL get_field(ddx_f, "DDX")
        CALL get_field(ddy_f, "DDY")

        npack = SIZE(packops_v3, dim=2)

        !$omp target teams distribute private(ff, ddx, ddy, i, j, k, jcount, &
        !$omp& igrid, lvl, ista, jsta, ksta, iend, jend, kend, buf_offset, &
        !$omp& istride, sum_ua, sum_a, idx)
        DO ipack = 1, npack
            CALL get_packop(ista, iend, jsta, jend, ksta, kend, &
                igrid, lvl, buf_offset, packops_v3(:, ipack))

            IF (lvl /= ilevel) CYCLE

            CALL get_grid3_real(ff, v3, igrid)
            CALL get_grid1_real(ddx, ddx_f, igrid)
            CALL get_grid1_real(ddy, ddy_f, igrid)

            k = ksta
            jcount = (jend - jsta) / 2 + 1

            IF (sum) THEN
                !$omp parallel do collapse(2) private(istride, idx)
                DO i = ista, iend, 2
                    DO j = jsta, jend, 2
                        ! Recompute to collapse later
                        istride = ((i - ista) / 2) * jcount
                        idx = buf_offset + istride + (j - jsta) / 2 + 1
                        sendbuf(idx) = &
                            ff(k, j, i) + ff(k, j+1, i) &
                            + ff(k, j, i+1) + ff(k, j+1, i+1)
                    END DO
                END DO
                !$omp end parallel do
            ELSE
                !$omp parallel do collapse(2) private(istride, sum_ua, sum_a, &
                !$omp& idx)
                DO i = ista, iend, 2
                    DO j = jsta, jend, 2
                        ! Recompute to collapse later
                        istride = ((i - ista) / 2) * jcount
                        sum_ua = ff(k, j, i) * ddx(i) * ddy(j) &
                            + ff(k, j+1, i) * ddx(i) * ddy(j+1) &
                            + ff(k, j, i+1) * ddx(i+1) * ddy(j) &
                            + ff(k, j+1, i+1) * ddx(i+1) * ddy(j+1)
                        sum_a = (ddx(i) + ddx(i+1)) * (ddy(j) + ddy(j+1))

                        idx = buf_offset + istride + (j - jsta) / 2 + 1
                        sendbuf(idx) = sum_ua / sum_a
                    END DO
                END DO
                !$omp end parallel do
            END IF
        END DO
        !$omp end target teams distribute
    END SUBROUTINE pack_v3


    ! Perform all Recv-calls
    SUBROUTINE recv_all(ilevel)
        ! Subroutine arguments
        INTEGER(intk), INTENT(in) :: ilevel

        ! Local variables
        INTEGER(intk) :: i, iprocnbr, igridf, iface
        INTEGER(int32) :: recvcounter, messagelength

        ! Post all receive calls
        recvcounter = 0
        messagelength = 0
        nrecv = 0

        DO i = 1, irecv
            ! Receiving grid
            igridf = recvconns(3, i)        ! The sender/fine grid - the one with the PAR
            IF (ilevel == level(igridf)) THEN
                iprocnbr = recvconns(1, i)  ! The sender process (fine side)
                iface = recvconns(5, i)     ! The face being sent - used to compute message length

                messagelength = messagelength + face_area(igridf, iface)

                IF (recvcounter + messagelength > idim_mg_bufs) THEN
                    CALL errr(__FILE__, __LINE__)
                END IF
            END IF

            IF (messagelength > 0) THEN
                IF (i == irecv) THEN
                    CALL post_recv(iprocnbr, messagelength, recvcounter)
                ELSE IF (recvconns(1, i + 1) /= iprocnbr) THEN
                    CALL post_recv(iprocnbr, messagelength, recvcounter)
                END IF
            END IF
        END DO
    END SUBROUTINE recv_all


    ! Perform a single Recv
    SUBROUTINE post_recv(iprocnbr, messagelength, recvcounter)
        ! Subroutine arguments
        INTEGER(int32), INTENT(in) :: iprocnbr
        INTEGER(int32), INTENT(inout) :: messagelength
        INTEGER(int32), INTENT(inout) :: recvcounter

        ! Local variables (for convenience)
        ! none...

        nrecv = nrecv + 1
        !$omp target data use_device_addr(recvbuf)
        CALL MPI_Irecv(recvbuf(recvcounter+1), messagelength, &
            mglet_mpi_real, iprocnbr, 1, MPI_COMM_WORLD, recvreqs(nrecv))
        !$omp end target data

        recvcounter = recvcounter + messagelength
        messagelength = 0
    END SUBROUTINE post_recv


    ! Perform all send calls
    SUBROUTINE send_all(ilevel)
        ! Subroutine arguments
        INTEGER(intk), INTENT(in) :: ilevel

        ! Local variables
        INTEGER(intk) :: i, iprocnbr, igridf, iface
        INTEGER(int32) :: sendcounter, messagelength

        sendcounter = 0
        messagelength = 0
        nsend = 0

        DO i = 1, isend
            igridf = sendconns(3, i)
            iprocnbr = sendconns(2, i)

            IF (ilevel == level(igridf)) THEN
                iface = sendconns(5, i)
                messagelength = messagelength + face_area(igridf, iface)

                IF (sendcounter + messagelength > idim_mg_bufs) THEN
                    CALL errr(__FILE__, __LINE__)
                END IF
            END IF

            IF (messagelength > 0) THEN
                IF (i == isend) THEN
                    CALL post_send(iprocnbr, messagelength, sendcounter)
                ELSE IF (sendconns(2, i + 1) /= iprocnbr) THEN
                    CALL post_send(iprocnbr, messagelength, sendcounter)
                END IF
            END IF
        END DO
    END SUBROUTINE send_all


    ! Perform a single send call
    SUBROUTINE post_send(iprocnbr, messagelength, sendcounter)
        ! Subroutine arguments
        INTEGER(int32), INTENT(in) :: iprocnbr
        INTEGER(int32), INTENT(inout) :: messagelength
        INTEGER(int32), INTENT(inout) :: sendcounter

        ! Local variables (for convenience)
        ! none...

        nsend = nsend + 1
        !$omp target data use_device_addr(sendbuf)
        CALL MPI_Isend(sendbuf(sendcounter + 1), messagelength, &
            mglet_mpi_real, iprocnbr, 1, MPI_COMM_WORLD, sendreqs(nsend))
        !$omp end target data

        sendcounter = sendcounter + messagelength
        messagelength = 0
    END SUBROUTINE post_send


    ! Process receive buffers after all receives complete,
    ! then wait for send buffers to be free
    SUBROUTINE process_bufs(ilevel, v1, v2, v3)
        ! Subroutine arguments
        INTEGER(intk), INTENT(in) :: ilevel
        TYPE(field_t), INTENT(inout) :: v1, v2, v3

        IF (nrecv > 0) THEN
            CALL MPI_Waitall(nrecv, recvreqs, MPI_STATUSES_IGNORE)

            CALL unpack_v1(ilevel, v1)
            CALL unpack_v2(ilevel, v2)
            CALL unpack_v3(ilevel, v3)
        END IF

        IF (nsend > 0) THEN
            CALL MPI_Waitall(nsend, sendreqs, MPI_STATUSES_IGNORE)
        END IF
    END SUBROUTINE process_bufs


    SUBROUTINE unpack_v1(ilevel, v1)
        ! Subroutine arguments
        INTEGER(intk), INTENT(in) :: ilevel
        TYPE(field_t), INTENT(inout) :: v1

        ! Local variables
        INTEGER(intk) :: ipack, npack, i, j, k, idx
        INTEGER(intk) :: jstride, kcount
        INTEGER(intk) :: ista, iend, jsta, jend, ksta, kend
        INTEGER(intk) :: igridc, levelf, buf_offset
        REAL(realk), POINTER, CONTIGUOUS :: fc(:, :, :)

        npack = SIZE(unpackops_v1, dim=2)

        !$omp target teams distribute private(fc, kcount, ista, iend, jsta, &
        !$omp& jend, ksta, kend, igridc, levelf, buf_offset, i, j, k, &
        !$omp& jstride, idx)
        DO ipack = 1, npack
            CALL get_unpackop(ista, iend, jsta, jend, ksta, kend, &
                igridc, levelf, buf_offset, unpackops_v1(:, ipack))

            IF (levelf /= ilevel) CYCLE
            CALL get_grid3_real(fc, v1, igridc)
            kcount = (kend - ksta) + 1

            !$omp parallel do collapse(3) private(jstride, idx)
            DO i = ista, iend
                DO j = jsta, jend
                    DO k = ksta, kend
                        ! Recompute to collapse later
                        jstride = (j - jsta) * kcount
                        idx = buf_offset + jstride + (k - ksta) + 1
                        fc(k, j, i) = recvbuf(idx)
                    END DO
                END DO
            END DO
            !$omp end parallel do
        END DO
        !$omp end target teams distribute
    END SUBROUTINE unpack_v1


    SUBROUTINE unpack_v2(ilevel, v2)
        ! Subroutine arguments
        INTEGER(intk), INTENT(in) :: ilevel
        TYPE(field_t), INTENT(inout) :: v2

        ! Local variables
        INTEGER(intk) :: ipack, npack, i, j, k, idx
        INTEGER(intk) :: istride, kcount
        INTEGER(intk) :: ista, iend, jsta, jend, ksta, kend
        INTEGER(intk) :: igridc, levelf, buf_offset
        REAL(realk), POINTER, CONTIGUOUS :: fc(:, :, :)

        npack = SIZE(unpackops_v2, dim=2)

        !$omp target teams distribute private(fc, kcount, ista, iend, jsta, &
        !$omp& jend, ksta, kend, igridc, levelf, buf_offset, i, j, k, &
        !$omp& istride, idx)
        DO ipack = 1, npack
            CALL get_unpackop(ista, iend, jsta, jend, ksta, kend, &
                igridc, levelf, buf_offset, unpackops_v2(:, ipack))

            IF (levelf /= ilevel) CYCLE
            CALL get_grid3_real(fc, v2, igridc)
            kcount = (kend - ksta) + 1

            !$omp parallel do collapse(3) private(istride, idx)
            DO i = ista, iend
                DO j = jsta, jend
                    DO k = ksta, kend
                        istride = (i - ista) * kcount
                        idx = buf_offset + istride + (k - ksta) + 1
                        fc(k, j, i) = recvbuf(idx)
                    END DO
                END DO
            END DO
            !$omp end parallel do
        END DO
        !$omp end target teams distribute
    END SUBROUTINE unpack_v2


    SUBROUTINE unpack_v3(ilevel, v3)
        ! Subroutine arguments
        INTEGER(intk), INTENT(in) :: ilevel
        TYPE(field_t), INTENT(inout) :: v3

        ! Local variables
        INTEGER(intk) :: ipack, npack, i, j, k, idx
        INTEGER(intk) :: istride, jcount
        INTEGER(intk) :: ista, iend, jsta, jend, ksta, kend
        INTEGER(intk) :: igridc, levelf, buf_offset
        REAL(realk), POINTER, CONTIGUOUS :: fc(:, :, :)

        npack = SIZE(unpackops_v3, dim=2)

        !$omp target teams distribute private(fc, jcount, ista, iend, jsta, &
        !$omp& jend, ksta, kend, igridc, levelf, buf_offset, i, j, k, &
        !$omp& istride, idx)
        DO ipack = 1, npack
            CALL get_unpackop(ista, iend, jsta, jend, ksta, kend, &
                igridc, levelf, buf_offset, unpackops_v3(:, ipack))

            IF (levelf /= ilevel) CYCLE
            CALL get_grid3_real(fc, v3, igridc)
            jcount = (jend - jsta) + 1

            !$omp parallel do collapse(3) private(istride, idx)
            DO i = ista, iend
                DO j = jsta, jend
                    DO k = ksta, kend
                        ! Recompute to collapse later
                        istride = (i - ista) * jcount
                        idx = buf_offset + istride + (j - jsta) + 1
                        fc(k, j, i) = recvbuf(idx)
                    END DO
                END DO
            END DO
            !$omp end parallel do
        END DO
        !$omp end target teams distribute
    END SUBROUTINE unpack_v3


    SUBROUTINE init_par_ftoc()
        ! Subroutine arguments
        ! none...

        ! Local variables
        LOGICAL :: one_connect
        INTEGER(intk) :: i, iface, ifacerecv, igrid, inbr, iprocnbr, pos
        INTEGER(intk) :: maxconns
        INTEGER(intk) :: neighbours(26)
        INTEGER(int32) :: ncols

        CHARACTER(len=8) :: ctyp
        INTEGER(int32), ALLOCATABLE :: sendcounts(:), sdispls(:)
        INTEGER(int32), ALLOCATABLE :: recvcounts(:), rdispls(:)

        CALL set_timer(211, "PAR_FTOC")

        ncols = 6
        maxconns = INT((nmygrids+1.0)*6.0*1.2)
        ALLOCATE(sendconns(ncols, maxconns))
        sendconns = 0

        ! The maximum number of concurrent communications are the number
        ! of processes
        ALLOCATE(sendreqs(numprocs))
        ALLOCATE(recvreqs(numprocs))

        ALLOCATE(sendcounts(0:numprocs-1))
        ALLOCATE(sdispls(0:numprocs-1))
        ALLOCATE(recvcounts(0:numprocs-1))
        ALLOCATE(rdispls(0:numprocs-1))
        sendcounts = 0
        sdispls = 0
        recvcounts = 0
        rdispls = 0

        nrecv = 0

        ! On the coarsest level there are not allowed to be any PAR BC's.
        DO i = 1, nmygrids
            igrid = mygrids(i)

            ! Loop over the boundary faces 1..6
            !
            ! Meaning of iterator:
            !   1 : FRONT  ( low X)
            !   2 : BACK   (high X)
            !   3 : RIGHT  ( low Y)
            !   4 : LEFT   (high Y)
            !   5 : BOTTOM ( low Z)
            !   6 : TOP    (high Z)
            DO iface = 1, 6
                ! We assume the PAR is always first/primary BC for that face
                CALL get_bc_ctyp(ctyp, 1, iface, igrid)
                IF (ctyp(1:3) /= "PAR") CYCLE

                ! The PAR BC needs always to transfer data to the parent grid,
                ! independent on the position within this grid.
                !
                ! Additionally, when the PAR is on top of a CON, we transport
                ! data to the neighbour grid as well.

                ! First - always send data to parent grid!
                inbr = iparent(igrid)
                iprocnbr = idprocofgrd(inbr)

                ! Position of fine grid within coarse grid
                SELECT CASE (iface)
                CASE (1, 2)
                    pos = iposition(igrid)
                CASE (3, 4)
                    pos = jposition(igrid)
                CASE (5, 6)
                    pos = kposition(igrid)
                END SELECT

                ! We allow a fine grid to completely refine a coarse grid
                ! in this case both PAR's of the fine grid can sit on top of
                ! a CON.
                !
                ! There are just two valid configurations:
                !   - iif == iic: This means the fine grid refine exactly half
                !     of the coarse grid in that direction. One PAR sits on
                !     top of a CON and another PAR is completely immersed in the
                !     coarse grid.
                !   - iif > iic: Actually, iif = (iic-4)*2 + 4, precisely.
                !     In this case the fine grid cover the entire extent of
                !     the coarse grid in the i-direction. Both PAR BC's sit on
                !     top of a CON in the coarse grid.
                one_connect = .TRUE.
                BLOCK
                    INTEGER :: iif, jjf, kkf, iic, jjc, kkc
                    CALL get_mgdims(kkf, jjf, iif, igrid)
                    CALL get_mgdims(kkc, jjc, iic, inbr)
                    SELECT CASE (iface)
                    CASE (1, 2)
                        IF (iif > iic) one_connect = .FALSE.
                    CASE (3, 4)
                        IF (jjf > jjc) one_connect = .FALSE.
                    CASE (5, 6)
                        IF (kkf > kkc) one_connect = .FALSE.
                    END SELECT
                END BLOCK

                ! In this "first" message, front send to front, back to back
                ! or it goes to somewhere internal to the grid (ifacerecv = -1)
                ifacerecv = iface
                IF (pos == 3 .AND. one_connect) THEN
                    IF (iface == 2 .OR. iface == 4 .OR. iface == 6) &
                        ifacerecv = -1
                ELSE IF (pos > 3) THEN
                    IF (iface == 1 .OR. iface == 3 .OR. iface == 5) &
                        ifacerecv = -1
                END IF

                nsend = nsend + 1
                IF (nsend > maxconns) THEN
                    write(*, *) "Number of PAR's exceeded on process ", myid
                    write(*, *) "maxconns =", maxconns, &
                        "nmygrids =", nmygrids, "nsend = ", nsend
                    CALL errr(__FILE__, __LINE__)
                END IF

                sendconns(1, nsend) = myid       ! Sending process (this process)
                sendconns(2, nsend) = iprocnbr   ! Receiving process (neighbour process)
                sendconns(3, nsend) = igrid      ! Sending grid grid (on current process)
                sendconns(4, nsend) = inbr       ! Receiving grid (on neighbour process)
                sendconns(5, nsend) = iface      ! Which face is sent (1..6)
                sendconns(6, nsend) = ifacerecv  ! Which face is received (1..6), -1 for internal

                sendcounts(iprocnbr) = sendcounts(iprocnbr) + ncols

                ! The same PAR face data can be sent a second time to a
                ! different grid and process.
                IF (pos == 3 .AND. one_connect) THEN
                    IF (iface == 2 .OR. iface == 4 .OR. iface == 6) CYCLE
                ELSE IF (pos > 3) THEN
                    IF (iface == 1 .OR. iface == 3 .OR. iface == 5) CYCLE
                END IF

                ! Second message/transfer follows!

                ! Get neighbours of the parent grid

                CALL get_neighbours(neighbours, inbr)
                inbr = neighbours(iface)
                iprocnbr = idprocofgrd(inbr)

                ! Get receiving face: front receive from back etc.
                SELECT CASE (iface)
                CASE (1)
                    ifacerecv = 2
                CASE (2)
                    ifacerecv = 1
                CASE (3)
                    ifacerecv = 4
                CASE (4)
                    ifacerecv = 3
                CASE (5)
                    ifacerecv = 6
                CASE (6)
                    ifacerecv = 5
                END SELECT

                nsend = nsend + 1
                IF (nsend > maxconns) THEN
                    write(*, *) "Number of PAR's exceeded on process ", myid
                    write(*, *) "maxconns =", maxconns, &
                        "nmygrids =", nmygrids, "nsend = ", nsend
                    CALL errr(__FILE__, __LINE__)
                END IF

                sendconns(1, nsend) = myid       ! Sending process (this process)
                sendconns(2, nsend) = iprocnbr   ! Receiving process (neighbour process)
                sendconns(3, nsend) = igrid      ! Sending grid (on current process)
                sendconns(4, nsend) = inbr       ! Receiving grid (on neighbour process)
                sendconns(5, nsend) = iface      ! Which face is sent (1..6)
                sendconns(6, nsend) = ifacerecv  ! Which face is received (1..6), -1 for internal

                sendcounts(iprocnbr) = sendcounts(iprocnbr) + ncols
            END DO
        END DO

        isend = nsend

        ! Sort sendconns by process ID (col 2)
        CALL sort_conns(sendconns(:, 1:nsend), 2)

        ! Calculate sdispl offset
        DO i = 1, numprocs-1
            sdispls(i) = sdispls(i-1) + sendcounts(i-1)
        END DO

        ! First exchange NUMBER OF ELEMENTS TO RECEIVE, to be able to
        ! calculate rdispls array
        CALL MPI_Alltoall(sendcounts, 1, MPI_INTEGER, recvcounts, 1, &
            MPI_INTEGER, MPI_COMM_WORLD)

        ! Calculate rdispl offset
        DO i = 1, numprocs-1
            rdispls(i) = rdispls(i-1) + recvcounts(i-1)
        END DO

        ! Allocate sendConns array
        irecv = (rdispls(numprocs-1) + recvcounts(numprocs-1))/ncols
        ALLOCATE(recvconns(ncols, irecv))
        recvconns = 0

        ! Exchange connection information
        CALL MPI_Alltoallv(sendconns, sendcounts, sdispls, MPI_INTEGER, &
            recvconns, recvcounts, rdispls, MPI_INTEGER, &
            MPI_COMM_WORLD)

        CALL precompute_packops()
        CALL precompute_unpackops()

        is_init = .TRUE.
        nsend = 0
        nrecv = 0

        ! Clean up
        DEALLOCATE(rdispls)
        DEALLOCATE(recvcounts)
        DEALLOCATE(sdispls)
        DEALLOCATE(sendcounts)
    END SUBROUTINE init_par_ftoc


    ! Calculate the area of a boundary face, i.e. the number of cells
    ! to be exchanged (all planes). Used to calculate message lengths.
    ! Needs to be given *fine grid* igrid and iface and the result is the
    ! number of values in the packed send buffer - i.e. the number of
    ! cells on the coarse side
    FUNCTION face_area(igrid, iface) RESULT(area)

        ! Function arguments
        INTEGER(intk) :: area
        INTEGER(intk), INTENT(in) :: igrid, iface

        ! Local variables
        INTEGER(intk) :: kstart, kstop, jstart, jstop, istart, istop

        CALL start_and_stop(kstart, kstop, jstart, jstop, istart, istop, &
            igrid, iface)

        area = ((istop-istart)/2+1)*((jstop-jstart)/2+1)*((kstop-kstart)/2+1)
    END FUNCTION face_area


    ! Indices from which to send data from the fine grid
    SUBROUTINE start_and_stop(kstart, kstop, jstart, jstop, istart, istop, &
            igrid, iface)

        ! Subroutine arguments
        INTEGER(intk), INTENT(out) :: kstart, kstop, jstart, jstop, &
            istart, istop
        INTEGER(intk), INTENT(in) :: igrid, iface

        ! Local variables
        INTEGER(intk) :: kk, jj, ii

        CALL get_mgdims(kk, jj, ii, igrid)

        ! Select entire grid - later reduce to only contain face
        kstart = 3
        jstart = 3
        istart = 3

        kstop = kk-2
        jstop = jj-2
        istop = ii-2

        ! Reduce selection to only contain face
        ! The velocities here are always solved at the fine grid - so at the
        ! front the fine velocities at the PAR is stored in pos. 2, at the back
        ! in pos. ii-3 and so on for the other faces.
        SELECT CASE (iface)
        CASE (1)
            istart = 2
            istop = 2
        CASE (2)
            istart = ii-2
            istop = ii-2
        CASE (3)
            jstart = 2
            jstop = 2
        CASE (4)
            jstart = jj-2
            jstop = jj-2
        CASE (5)
            kstart = 2
            kstop = 2
        CASE (6)
            kstart = kk-2
            kstop = kk-2
        END SELECT
    END SUBROUTINE start_and_stop


    SUBROUTINE get_packop(ista, iend, jsta, jend, ksta, kend, igrid, &
        lvl, buf_offset, packop)
        ! Subroutine arguments
        INTEGER(intk), INTENT(out) :: ista, iend, jsta, jend, ksta, kend
        INTEGER(intk), INTENT(out) :: igrid, lvl, buf_offset
        INTEGER(intk), INTENT(in) :: packop(9)

        ista = packop(1)
        iend = packop(2)
        jsta = packop(3)
        jend = packop(4)
        ksta = packop(5)
        kend = packop(6)
        igrid = packop(7)
        lvl = packop(8)
        buf_offset = packop(9)
    END SUBROUTINE get_packop


    SUBROUTINE get_unpackop(ista, isto, jsta, jsto, ksta, ksto, igridc, &
        levelf, buf_offset, unpackop)
        ! Subroutine arguments
        INTEGER(intk), INTENT(out) :: ista, isto, jsta, jsto, ksta, ksto
        INTEGER(intk), INTENT(out) :: igridc, levelf, buf_offset
        INTEGER(intk), INTENT(in) :: unpackop(9)

        ista = unpackop(1)
        isto = unpackop(2)
        jsta = unpackop(3)
        jsto = unpackop(4)
        ksta = unpackop(5)
        ksto = unpackop(6)
        igridc = unpackop(7)
        levelf = unpackop(8)
        buf_offset = unpackop(9)
    END SUBROUTINE get_unpackop


    SUBROUTINE finish_par_ftoc()
        DEALLOCATE(sendconns)
        DEALLOCATE(recvconns)
        DEALLOCATE(sendreqs)
        DEALLOCATE(recvreqs)
        !$omp target exit data map(always, delete: packops_v1, packops_v2, &
        !$omp& packops_v3)
        !$omp target exit data map(always, delete: unpackops_v1, unpackops_v2, &
        !$omp& unpackops_v3)
        DEALLOCATE(packops_v1)
        DEALLOCATE(packops_v2)
        DEALLOCATE(packops_v3)
        DEALLOCATE(unpackops_v1)
        DEALLOCATE(unpackops_v2)
        DEALLOCATE(unpackops_v3)
        is_init = .FALSE.
    END SUBROUTINE finish_par_ftoc
END MODULE par_ftoc_mod
