MODULE par_ftoc1_mod
    USE MPI_f08
    USE core_mod
    USE par_ftoc_core_mod

    IMPLICIT NONE (type, external)
    PRIVATE

    ! Lists that hold the send and receive request arrays
    TYPE(MPI_Request), ALLOCATABLE :: sendreqs(:), recvreqs(:)

    ! Lists that hold the messages that are ACTUALLY sendt and received
    INTEGER(intk) :: nsend, nrecv
    INTEGER(int32), ALLOCATABLE :: recvlist(:)
    INTEGER(intk), ALLOCATABLE :: recvidxlist(:, :)

    ! Variable to indicate if the connection information has
    ! been created.
    LOGICAL :: is_init = .FALSE.

    PUBLIC :: par_ftoc1, init_par_ftoc1, finish_par_ftoc1

CONTAINS

    ! Main parent function
    SUBROUTINE par_ftoc1(ilevel, v1, v2, v3, sum)

        ! Subroutine arguments
        INTEGER(intk), INTENT(in) :: ilevel
        TYPE(field_t), INTENT(inout) :: v1, v2, v3
        LOGICAL, OPTIONAL, INTENT(in) :: sum

        ! Local variables
        LOGICAL :: sumflag

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
    END SUBROUTINE par_ftoc1


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


    SUBROUTINE init_par_ftoc1()
        ! Subroutine arguments
        ! none...

        ! Local variables
        INTEGER(intk) :: maxconns

        IF (is_init) CALL errr(__FILE__, __LINE__)

        maxconns = INT((nmygrids+1.0)*6.0*1.2)

        ! The maximum number of concurrent communications are the number
        ! of processes
        ALLOCATE(recvidxlist(3, maxconns))
        ALLOCATE(recvlist(numprocs))
        ALLOCATE(sendreqs(numprocs))
        ALLOCATE(recvreqs(numprocs))
        recvidxlist = 0
        recvlist = 0

        nsend = 0
        nrecv = 0
        is_init = .TRUE.
    END SUBROUTINE init_par_ftoc1


    SUBROUTINE finish_par_ftoc1()
        ! Subroutine arguments
        ! none...

        ! Local variables
        ! none...

        IF (.NOT. is_init) CALL errr(__FILE__, __LINE__)

        DEALLOCATE(sendreqs)
        DEALLOCATE(recvreqs)
        DEALLOCATE(recvlist)
        DEALLOCATE(recvidxlist)
        is_init = .FALSE.
    END SUBROUTINE finish_par_ftoc1

END MODULE par_ftoc1_mod
