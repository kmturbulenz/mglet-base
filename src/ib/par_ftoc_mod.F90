MODULE par_ftoc_mod
    USE MPI_f08
    USE core_mod

    IMPLICIT NONE (type, external)
    PRIVATE

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

    ! Lists that hold the messages that are ACTUALLY sendt and received
    INTEGER(intk) :: nsend, nrecv
    INTEGER(int32), ALLOCATABLE :: recvlist(:)
    INTEGER(intk), ALLOCATABLE :: recvidxlist(:, :)

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
        CALL send_all(ilevel, v1, v2, v3, sumflag)
        CALL process_bufs(v1, v2, v3)

        nrecv = 0
        nsend = 0

        CALL stop_timer(211)
    END SUBROUTINE par_ftoc_norm


    ! Perform all Recv-calls
    SUBROUTINE recv_all(ilevel)
        ! Subroutine arguments
        INTEGER(intk), INTENT(in) :: ilevel

        ! Local variables
        INTEGER(intk) :: i, iprocnbr, igridf, iface, facearea
        INTEGER(int32) :: recvcounter, messagelength

        ! Post all receive calls
        recvcounter = 0
        messagelength = 0
        nrecv = 0
        recvidxlist = -HUGE(1_intk)
        recvlist = 0

        DO i = 1, irecv
            ! Receiving grid
            igridf = recvconns(3, i)        ! The sender/fine grid - the one with the PAR
            IF (ilevel == level(igridf)) THEN
                iprocnbr = recvconns(1, i)  ! The sender process (fine side)
                iface = recvconns(5, i)     ! The face being sent - used to compute message length

                facearea = face_area(igridf, iface)
                recvidxlist(1, i) = iprocnbr
                recvidxlist(2, i) = facearea
                recvidxlist(3, i) = recvcounter + messagelength
                messagelength = messagelength + facearea

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
        ! Identifier of receive connection
        INTEGER(int32), INTENT(in) :: iprocnbr
        INTEGER(int32), INTENT(inout) :: messagelength
        INTEGER(int32), INTENT(inout) :: recvcounter

        ! Local variables (for convenience)
        ! none...

        nrecv = nrecv + 1
        recvlist(nrecv) = iprocnbr
        CALL MPI_Irecv(recvbuf(recvcounter+1), messagelength, &
            mglet_mpi_real, iprocnbr, 1, MPI_COMM_WORLD, recvreqs(nrecv))

        recvcounter = recvcounter + messagelength
        messagelength = 0
    END SUBROUTINE post_recv


    ! Perform all send calls
    SUBROUTINE send_all(ilevel, v1, v2, v3, sum)
        ! Subroutine arguments
        INTEGER(intk), INTENT(in) :: ilevel
        TYPE(field_t), INTENT(inout) :: v1, v2, v3
        LOGICAL, INTENT(in) :: sum

        ! Local variables
        INTEGER(intk) :: i, iprocnbr, igridf
        INTEGER(int32) :: sendcounter, messagelength

        ! Pack all buffers and send data
        sendcounter = 0
        messagelength = 0
        nSend = 0

        DO i = 1, isend
            igridf = sendconns(3, i)
            iprocnbr = sendconns(2, i)
            IF (ilevel == level(igridf)) THEN
                CALL write_buffer(i, v1, v2, v3, sum, messagelength, &
                    sendcounter)
            END IF

            IF (messagelength > 0) THEN
                IF (i == iSend) THEN
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
        CALL MPI_Isend(sendbuf(sendcounter + 1), messagelength, &
            mglet_mpi_real, iprocnbr, 1, MPI_COMM_WORLD, sendreqs(nsend))

        sendcounter = sendcounter + messagelength
        messagelength = 0
    END SUBROUTINE post_send


    ! Write Send buffers
    !
    ! Write the relevant fields into the send buffers
    SUBROUTINE write_buffer(id, v1, v2, v3, sum, messagelength, sendcounter)
        ! Subroutine arguments
        INTEGER(int32), INTENT(in) :: id
        TYPE(field_t), INTENT(inout), TARGET :: v1, v2, v3
        LOGICAL, INTENT(in) :: sum
        INTEGER(int32), INTENT(inout) :: messagelength
        INTEGER(int32), INTENT(in) :: sendcounter

        ! Local variables
        TYPE(field_t), POINTER :: field
        INTEGER(intk) :: igridf, iface
        INTEGER(int32) :: facearea, offset

        ! Set variables from send table - *fine* grid and face
        igridf = sendconns(3, id)
        iface = sendconns(5, id)

        ! Check that buffer does not overflow
        facearea = face_area(igridf, iface)
        IF (sendcounter + messagelength + facearea > idim_mg_bufs) THEN
            CALL errr(__FILE__, __LINE__)
        END IF

        ! Only send interface-normal vector element
        SELECT CASE (iface)
        CASE (1, 2)
            field => v1
        CASE (3, 4)
            field => v2
        CASE (5, 6)
            field => v3
        END SELECT

        ! Pack
        offset = sendcounter + messagelength + 1
        CALL pack_single(sendbuf(offset:offset+facearea-1), field, &
            igridf, iface, sum)
        messagelength = messagelength + facearea
    END SUBROUTINE write_buffer


    SUBROUTINE pack_single(buf, field, igrid, iface, sum)
        ! Subroutine arguments
        REAL(realk), INTENT(inout), CONTIGUOUS :: buf(:)
        TYPE(field_t), INTENT(inout) :: field
        INTEGER(intk), INTENT(in) :: igrid
        INTEGER(intk), INTENT(in) :: iface
        LOGICAL, INTENT(in) :: sum

        ! Local variables
        INTEGER(intk) :: ista, isto, jsta, jsto, ksta, ksto
        INTEGER(intk) :: k, j, i
        INTEGER(intk) :: icount
        REAL(realk) :: sum_ua, sum_a
        REAL(realk), POINTER, CONTIGUOUS :: ff(:, :, :)
        REAL(realk), POINTER, CONTIGUOUS :: ddx(:), ddy(:), ddz(:)

        icount = 0
        CALL field%get_ptr(ff, igrid)
        CALL start_and_stop(ksta, ksto, jsta, jsto, ista, isto, igrid, iface)

        IF (sum) THEN
            SELECT CASE (iface)
            CASE (1, 2)
                i = ista  ! ista and isto are the same in this case
                DO j = jsta, jsto, 2
                    DO k = ksta, ksto, 2
                        icount = icount + 1
                        buf(icount) = ff(k, j, i) + ff(k, j+1, i) &
                            + ff(k+1, j, i) + ff(k+1, j+1, i)
                    END DO
                END DO
            CASE (3, 4)
                j = jsta  ! jsta and jsto are the same in this case
                DO i = ista, isto, 2
                    DO k = ksta, ksto, 2
                        icount = icount + 1
                        buf(icount) = ff(k, j, i) + ff(k, j, i+1) &
                            + ff(k+1, j, i) + ff(k+1, j, i+1)
                    END DO
                END DO
            CASE (5, 6)
                k = ksta  ! ksta and ksto are the same in this case
                DO i = ista, isto, 2
                    DO j = jsta, jsto, 2
                        icount = icount + 1
                        buf(icount) = ff(k, j, i) + ff(k, j+1, i) &
                            + ff(k, j, i+1) + ff(k, j+1, i+1)
                    END DO
                END DO
            END SELECT
        ELSE
            CALL get_fieldptr(ddx, "DDX", igrid)
            CALL get_fieldptr(ddy, "DDY", igrid)
            CALL get_fieldptr(ddz, "DDZ", igrid)

            SELECT CASE (iface)
            CASE (1, 2)
                i = ista  ! ista and isto are the same in this case
                DO j = jsta, jsto, 2
                    DO k = ksta, ksto, 2
                        sum_ua = ff(k, j, i)*ddy(j)*ddz(k) &
                            + ff(k, j+1, i)*ddy(j+1)*ddz(k) &
                            + ff(k+1, j, i)*ddy(j)*ddz(k+1) &
                            + ff(k+1, j+1, i)*ddy(j+1)*ddz(k+1)

                        sum_a = (ddy(j) + ddy(j+1))*(ddz(k) + ddz(k+1))

                        icount = icount + 1
                        buf(icount) = sum_ua/sum_a
                    END DO
                END DO
            CASE (3, 4)
                j = jsta  ! jsta and jsto are the same in this case
                DO i = ista, isto, 2
                    DO k = ksta, ksto, 2
                        sum_ua = ff(k, j, i)*ddx(i)*ddz(k) &
                            + ff(k, j, i+1)*ddx(i+1)*ddz(k) &
                            + ff(k+1, j, i)*ddx(i)*ddz(k+1) &
                            + ff(k+1, j, i+1)*ddx(i+1)*ddz(k+1)

                        sum_a = (ddx(i) + ddx(i+1))*(ddz(k) + ddz(k+1))

                        icount = icount + 1
                        buf(icount) = sum_ua/sum_a
                    END DO
                END DO
            CASE (5, 6)
                k = ksta  ! ksta and ksto are the same in this case
                DO i = ista, isto, 2
                    DO j = jsta, jsto, 2
                        sum_ua = ff(k, j, i)*ddx(i)*ddy(j) &
                            + ff(k, j+1, i)*ddx(i)*ddy(j+1) &
                            + ff(k, j, i+1)*ddx(i+1)*ddy(j) &
                            + ff(k, j+1, i+1)*ddx(i+1)*ddy(j+1)

                        sum_a = (ddx(i) + ddx(i+1))*(ddy(j) + ddy(j+1))

                        icount = icount + 1
                        buf(icount) = sum_ua/sum_a
                    END DO
                END DO
            END SELECT
        END IF
    END SUBROUTINE pack_single


    ! Read Receive buffers
    !
    ! Write the contents of the receive buffers back in their
    ! matching fields
    SUBROUTINE read_buffer(id, v1, v2, v3)

        ! Subroutine arguments
        INTEGER(intk), INTENT(in) :: id
        TYPE(field_t), INTENT(inout), TARGET :: v1, v2, v3

        ! Local variables
        TYPE(field_t), POINTER :: field
        INTEGER(intk) :: igridf, igridc, iface, ifacerecv
        INTEGER(int32) :: facearea
        INTEGER(int32) :: offset

        ! Set variables from send table - *fine* grid and face
        igridf = recvconns(3, id)

        ! Only send interface-normal vector element
        iface = recvconns(5, id)  ! Face being sent - not received!
        SELECT CASE (iface)
        CASE (1, 2)
            field => v1
        CASE (3, 4)
            field => v2
        CASE (5, 6)
            field => v3
        END SELECT

        ! Receiving face is difference from sending face - can also be
        ! internal - in that case it is -1
        igridc = recvconns(4, id)
        ifacerecv = recvconns(6, id)

        ! Unpack
        offset = recvidxlist(3, id) + 1
        facearea = face_area(igridf, iface)
        CALL unpack_single(recvbuf(offset:offset+facearea-1), field, &
            igridf, igridc, iface, ifacerecv)
    END SUBROUTINE read_buffer


    SUBROUTINE unpack_single(buf, field, igridf, igridc, iface, ifacerecv)
        ! Subroutine arguments
        REAL(realk), INTENT(in), CONTIGUOUS :: buf(:)
        TYPE(field_t), INTENT(inout) :: field
        INTEGER(intk), INTENT(in) :: igridf
        INTEGER(intk), INTENT(in) :: igridc
        INTEGER(intk), INTENT(in) :: iface
        INTEGER(intk), INTENT(in) :: ifacerecv

        ! Local variables
        INTEGER(intk) :: ista, isto, jsta, jsto, ksta, ksto
        INTEGER(intk) :: ipos, jpos, kpos
        INTEGER(intk) :: k, j, i
        INTEGER(intk) :: kkc, jjc, iic
        INTEGER(intk) :: kkf, jjf, iif
        INTEGER(intk) :: icount
        REAL(realk), POINTER, CONTIGUOUS :: fc(:, :, :)

        CALL get_mgdims(kkf, jjf, iif, igridf)   ! Dimensions of fine grid
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
            SELECT CASE (iface)  ! The sending face determine position
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
            SELECT CASE (ifacerecv)  ! The receiving face determine position
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

        ! Unpack
        icount = 0
        CALL field%get_ptr(fc, igridc)
        DO i = ista, isto
            DO j = jsta, jsto
                DO k = ksta, ksto
                    icount = icount + 1
                    fc(k, j, i) = buf(icount)
                END DO
            END DO
        END DO
    END SUBROUTINE unpack_single


    ! Process receive buffers as they arrive, wait for send
    ! buffers to be free
    SUBROUTINE process_bufs(v1, v2, v3)
        ! Subroutine arguments
        TYPE(field_t), INTENT(inout) :: v1, v2, v3

        ! Local variables
        INTEGER(int32) :: idx, i
        TYPE(MPI_Status) :: recvstatus
        INTEGER(int32) :: recvmessagelen
        INTEGER(int32) :: unpacklen

        DO WHILE (.TRUE.)
            IF (nrecv == 0) EXIT
            CALL MPI_Waitany(nrecv, recvreqs, idx, recvstatus)
            IF (idx == MPI_UNDEFINED) EXIT

            CALL MPI_Get_count(recvstatus, mglet_mpi_real, recvmessagelen)

            unpacklen = 0
            DO i = 1, irecv
                IF (recvidxlist(1, i) == recvlist(idx)) THEN
                    CALL read_buffer(i, v1, v2, v3)
                    unpacklen = unpacklen + recvidxlist(2, i)
                END IF
            END DO

            IF (recvmessagelen /= unpacklen) THEN
                CALL errr(__FILE__, __LINE__)
            END IF
        END DO
        CALL MPI_Waitall(nsend, sendreqs, MPI_STATUSES_IGNORE)
    END SUBROUTINE process_bufs


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
        ALLOCATE(recvidxlist(3, maxconns))
        ALLOCATE(recvlist(numprocs))
        ALLOCATE(sendreqs(numprocs))
        ALLOCATE(recvreqs(numprocs))
        recvidxlist = 0
        recvlist = 0

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


    SUBROUTINE finish_par_ftoc()
        DEALLOCATE(sendconns)
        DEALLOCATE(recvconns)
        DEALLOCATE(sendreqs)
        DEALLOCATE(recvreqs)
        DEALLOCATE(recvlist)
        DEALLOCATE(recvidxlist)
        is_init = .FALSE.
    END SUBROUTINE finish_par_ftoc


    SUBROUTINE sort_conns(list, col)
        ! Subroutine arguments
        INTEGER(intk), INTENT(inout) :: list(:, :)
        INTEGER(intk), INTENT(in) :: col

        ! Local variables
        INTEGER(intk) :: i, n
        INTEGER(intk), ALLOCATABLE :: tmplist(:, :)
        INTEGER(intk), ALLOCATABLE :: idx(:)

        ! sortix does not like n = 0 - this happens in one-level testcases
        n = SIZE(list, 2)
        IF (n == 0) RETURN

        ALLOCATE(idx(n))
        CALL sortix(n, list(col, :), idx)

        ! Transfer data to sorted list
        ALLOCATE(tmplist, SOURCE=list)
        DO i = 1, n
            list(:, i) = tmplist(:, idx(i))
        END DO
        DEALLOCATE(tmplist)
    END SUBROUTINE sort_conns

END MODULE par_ftoc_mod
