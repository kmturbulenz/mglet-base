MODULE parent1_mod
    USE MPI_f08
    USE core_mod
    USE parent_core_mod

    IMPLICIT NONE (type, external)
    PRIVATE

    ! Lists that hold the send and receive request arrays
    TYPE(MPI_Request), ALLOCATABLE :: sendreqs(:), recvreqs(:)

    ! Lists that hold the messages that are ACTUALLY sendt and received
    INTEGER(intk) :: nsend, nrecv
    INTEGER(int32), ALLOCATABLE :: sendlist(:), recvlist(:)
    INTEGER(intk), ALLOCATABLE :: recvidxlist(:, :)

    ! Variable to indicate if the connection information has been created.
    LOGICAL :: is_init = .FALSE.

    PUBLIC :: parent1, init_parent1, finish_parent1
CONTAINS

    ! Main parent function
    SUBROUTINE parent1(ilevel, v1, v2, v3, s1, s2, s3, normal)

        ! Subroutine arguments
        INTEGER(intk), INTENT(in) :: ilevel
        TYPE(field_t), OPTIONAL, INTENT(inout) :: v1, v2, v3, s1, s2, s3
        LOGICAL, OPTIONAL, INTENT(in) :: normal

        ! Local variables
        LOGICAL :: sn
        INTEGER(intk) :: nvars

        IF (.NOT. is_init) CALL errr(__FILE__, __LINE__)

        ! Check that no other transfers are in progress
        IF (nsend > 0 .OR. nrecv > 0) THEN
            WRITE(*, *) "Other transfer in progress."
            CALL errr(__FILE__, __LINE__)
        END IF

        ! Count number of variables to send/receive
        IF ((PRESENT(v1) .NEQV. PRESENT(v2)) .OR. &
                (PRESENT(v1) .NEQV. PRESENT(v3))) THEN
            CALL errr(__FILE__, __LINE__)
        END IF

        sn = .FALSE.
        IF (PRESENT(normal)) sn = normal

        nvars = 0
        IF (PRESENT(v1)) nvars = nvars + 1
        IF (PRESENT(v2)) nvars = nvars + 1
        IF (PRESENT(v3)) nvars = nvars + 1
        IF (PRESENT(s1)) nvars = nvars + 1
        IF (PRESENT(s2)) nvars = nvars + 1
        IF (PRESENT(s3)) nvars = nvars + 1

        ! If one vector component is present, all must be present
        IF (sn .AND. PRESENT(v1)) nvars = nvars - 2

        CALL recv_all(ilevel, nvars)
        CALL send_all(ilevel, v1, v2, v3, s1, s2, s3, normal)
        CALL process_bufs(v1, v2, v3, s1, s2, s3, normal)

        ! Clear counters and unset pointers. This is important to avoid
        ! problems at next function call
        nrecv = 0
        nsend = 0
    END SUBROUTINE parent1


    ! Perform all Recv-calls
    SUBROUTINE recv_all(ilevel, nvars)
        ! Subroutine arguments
        INTEGER(intk), INTENT(in) :: ilevel
        INTEGER(intk), INTENT(in) :: nvars

        ! Local variables
        INTEGER(intk) :: i, iprocnbr, igrid, iface, ilevelgrid, facearea
        INTEGER(int32) :: recvcounter, messagelength

        ! Post all receive calls
        recvcounter = 0
        messagelength = 0
        nrecv = 0
        recvidxlist = 0

        DO i = 1, irecv
            iprocnbr = recvconns(2, i)
            igrid = recvconns(3, i)
            iface = recvconns(5, i)
            ilevelgrid = level(igrid)

            IF (ilevel == ilevelgrid) THEN
                facearea = face_area(igrid, iface)
                recvidxlist(1, i) = iprocnbr
                recvidxlist(2, i) = nvars*facearea
                recvidxlist(3, i) = recvcounter + messagelength
                messagelength = messagelength + nvars*facearea

                IF (recvcounter + messagelength > idim_mg_bufs) THEN
                    CALL errr(__FILE__, __LINE__)
                END IF
            END IF

            IF (messagelength > 0) THEN
                IF (i == irecv) THEN
                    CALL post_recv(iprocnbr, messagelength, recvcounter)
                ELSE IF (recvconns(2, i + 1) /= iprocnbr) THEN
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
    SUBROUTINE send_all(ilevel, v1, v2, v3, s1, s2, s3, normal)
        ! Subroutine arguments
        INTEGER(intk), INTENT(in) :: ilevel
        TYPE(field_t), OPTIONAL, INTENT(inout) :: &
            v1, v2, v3, s1, s2, s3
        LOGICAL, OPTIONAL, INTENT(in) :: normal

        ! Local variables
        INTEGER(intk) :: i, iprocnbr, igrid, ilevelgrid
        INTEGER(int32) :: sendcounter, messagelength

        ! Pack all buffers and send data
        sendcounter = 0
        messagelength = 0
        nsend = 0

        DO i = 1, isend
            iprocnbr = sendconns(1, i)
            igrid = sendconns(3, i)
            ilevelgrid = level(igrid)

            IF (ilevel == ilevelgrid) THEN
                CALL write_buffer(i, messagelength, sendcounter, &
                    v1, v2, v3, s1, s2, s3, normal)
            END IF

            IF (messagelength > 0) THEN
                IF (i == isend) THEN
                    CALL post_send(iprocnbr, messagelength, sendcounter)
                ELSE IF (sendconns(1, i + 1) /= iprocnbr) THEN
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
        sendlist(nsend) = iprocnbr

        CALL MPI_Isend(sendbuf(sendcounter + 1), messagelength, &
            mglet_mpi_real, iprocnbr, 1, MPI_COMM_WORLD, sendreqs(nsend))

        sendcounter = sendcounter + messagelength
        messagelength = 0
    END SUBROUTINE post_send


    ! Write Send buffers
    !
    ! Write the relevant fields into the send buffers
    SUBROUTINE write_buffer(sendid, messagelength, sendcounter, &
            v1, v2, v3, s1, s2, s3, normal)
        ! Subroutine arguments
        INTEGER(int32), INTENT(in) :: sendid
        INTEGER(int32), INTENT(inout) :: messagelength
        INTEGER(int32), INTENT(in) :: sendcounter
        TYPE(field_t), OPTIONAL, INTENT(inout) :: &
            v1, v2, v3, s1, s2, s3
        LOGICAL, OPTIONAL, INTENT(in) :: normal

        ! Local variables
        INTEGER(intk) :: igrid, igridc, iface, nvars
        INTEGER(int32) :: thismessagelength, facearea, offset
        LOGICAL :: exU, exV, exW, sn

        ! Set variables from send table - *fine* grid and face
        igrid = sendconns(3, sendid)
        igridc = sendconns(4, sendid)
        iface = sendconns(5, sendid)

        sn = .FALSE.
        IF (PRESENT(normal)) sn = normal

        ! Which vectors to exchange
        exU = (sn .AND. iface < 3) .OR. (.NOT. sn)
        exV = (sn .AND. (iface > 2 .AND. iface < 5)) .OR. (.NOT. sn)
        exW = (sn .AND. iface > 4) .OR. (.NOT. sn)

        ! Count number of variables to send
        nvars = 0
        IF (PRESENT(v1) .AND. exU) nvars = nvars + 1
        IF (PRESENT(v2) .AND. exV) nvars = nvars + 1
        IF (PRESENT(v3) .AND. exW) nvars = nvars + 1
        IF (PRESENT(s1)) nvars = nvars + 1
        IF (PRESENT(s2)) nvars = nvars + 1
        IF (PRESENT(s3)) nvars = nvars + 1

        ! Face area
        facearea = face_area(igrid, iface)
        thismessagelength = nvars*facearea

        ! Check that buffer does not overflow
        IF (sendcounter + messagelength + thismessagelength > idim_mg_bufs) THEN
            CALL errr(__FILE__, __LINE__)
        END IF

        ! Reset message size counter
        offset = sendcounter + messagelength + 1

        ! Fill buffers
        IF (PRESENT(v1) .AND. exU) THEN
            CALL pack_single(sendbuf(offset:offset+facearea-1), v1, igrid, &
                igridc, iface)
            offset = offset + facearea
        END IF

        IF (PRESENT(v2) .AND. exV) THEN
            CALL pack_single(sendbuf(offset:offset+facearea-1), v2, igrid, &
                igridc, iface)
            offset = offset + facearea
        END IF

        IF (PRESENT(v3) .AND. exW) THEN
            CALL pack_single(sendbuf(offset:offset+facearea-1), v3, igrid, &
                igridc, iface)
            offset = offset + facearea
        END IF

        IF (PRESENT(s1)) THEN
            CALL pack_single(sendbuf(offset:offset+facearea-1), s1, igrid, &
                igridc, iface)
            offset = offset + facearea
        END IF

        IF (PRESENT(s2)) THEN
            CALL pack_single(sendbuf(offset:offset+facearea-1), s2, igrid, &
                igridc, iface)
            offset = offset + facearea
        END IF

        IF (PRESENT(s3)) THEN
            CALL pack_single(sendbuf(offset:offset+facearea-1), s3, igrid, &
                igridc, iface)
            offset = offset + facearea
        END IF

        IF (offset /= sendcounter + messagelength + thismessagelength + 1) THEN
            WRITE(*, *) "offset:", offset, &
                "expected:", sendcounter + messagelength + thismessagelength + 1
            CALL errr(__FILE__, __LINE__)
        END IF

        messagelength = messagelength + thismessagelength
    END SUBROUTINE write_buffer


    SUBROUTINE pack_single(buf, field, igrid, igridc, iface)
        ! Subroutine arguments
        REAL(realk), INTENT(inout), CONTIGUOUS :: buf(:)
        TYPE(field_t), INTENT(inout) :: field
        INTEGER(intk), INTENT(in) :: igrid, igridc, iface

        ! Local variables
        REAL(realk), POINTER, CONTIGUOUS :: fc(:, :, :)
        INTEGER(intk) :: k, j, i, icomp, icount
        INTEGER(intk) :: ista, isto, jsta, jsto, ksta, ksto

        CALL field%get_ptr(fc, igridc)

        icomp = 0
        SELECT CASE(iface)
        CASE (1, 2)
            IF (field%istag == 1) icomp = 1
        CASE (3, 4)
            IF (field%jstag == 1) icomp = 2
        CASE (5, 6)
            IF (field%kstag == 1) icomp = 3
        END SELECT

        CALL start_and_stop(igrid, iface, icomp, ista, isto, &
            jsta, jsto, ksta, ksto)

        icount = 0
        DO i = ista, isto
            DO j = jsta, jsto
                DO k = ksta, ksto
                    icount = icount + 1
                    buf(icount) = fc(k, j, i)
                END DO
            END DO
        END DO

        IF (icount /= SIZE(buf)) THEN
            WRITE(*, *) "icount:", icount, "SIZE(buf):", SIZE(buf)
            CALL errr(__FILE__, __LINE__)
        END IF
    END SUBROUTINE pack_single


    ! Read Receive buffers
    !
    ! Write the contents of the receive buffers back in their
    ! matching fields
    SUBROUTINE read_buffer(recvid, v1, v2, v3, s1, s2, s3, normal)

        ! Subroutine arguments
        INTEGER(intk), INTENT(in) :: recvid
        TYPE(field_t), OPTIONAL, INTENT(inout) :: &
            v1, v2, v3, s1, s2, s3
        LOGICAL, OPTIONAL, INTENT(in) :: normal

        ! Grid dimensions and pointers
        INTEGER(intk) :: kk, jj, ii
        INTEGER(intk) :: kkc, jjc, iic
        INTEGER(intk) :: idx
        INTEGER(intk) :: jj2d, ii2d, jjc2d, iic2d

        ! Grid to send from
        ! Must be intk because it intreface with MGLET
        INTEGER(intk) :: igrid, iface

        ! Message sizes
        ! Must be int32 because it interface with MPI
        INTEGER(int32) :: offset

        INTEGER(intk) :: facearea, tmp_buf_size
        INTEGER(intk) :: ustag1, ustag2, vstag1, vstag2, wstag1, wstag2
        REAL(realk), POINTER, CONTIGUOUS :: faceptr(:, :)
        REAL(realk), ALLOCATABLE :: tmp_buf(:)
        LOGICAL :: exU, exV, exW, sn

        ! Set variables from send table
        igrid = recvconns(3, recvid)
        iface = recvconns(5, recvid)

        sn = .FALSE.
        IF (PRESENT(normal)) sn = normal

        ! Which vectors to exchange
        exU = (sn .AND. iface < 3) .OR. (.NOT. sn)
        exV = (sn .AND. (iface > 2 .AND. iface < 5)) .OR. (.NOT. sn)
        exW = (sn .AND. iface > 4) .OR. (.NOT. sn)

        ! Get grid dimentsions and pointers
        CALL get_mgdims(kk, jj, ii, igrid)

        ! Get start- and stop indices of grid
        CALL idx2d(kk, jj, ii, iface, kkc, jjc, iic, jj2d, ii2d, jjc2d, iic2d)
        CALL stag(iface, ustag1, ustag2, vstag1, vstag2, wstag1, wstag2)
        facearea = face_area(igrid, iface)

        ! Allocate temporary work buffer
        tmp_buf_size = MAX(ii*jj, ii*kk, jj*kk)
        ALLOCATE(tmp_buf(tmp_buf_size))
        tmp_buf = 0.0

        ! Offset in receive buffer
        offset = recvidxlist(3, recvid) + 1
        idx = 0

        IF (PRESENT(v1) .AND. exU) THEN
            CALL v1%get_buffer(faceptr, igrid, iface)

            CALL prolong1(jjc2d, iic2d, jj2d, ii2d, &
                recvbuf(offset+idx:offset+idx+facearea-1), tmp_buf, ustag1)
            CALL prolong2(jjc2d, iic2d, jj2d, ii2d, tmp_buf, &
                faceptr, ustag2)
            idx = idx + facearea
        END IF

        IF (PRESENT(v2) .AND. exV) THEN
            CALL v2%get_buffer(faceptr, igrid, iface)

            CALL prolong1(jjc2d, iic2d, jj2d, ii2d, &
                recvbuf(offset+idx:offset+idx+facearea-1), tmp_buf, vstag1)
            CALL prolong2(jjc2d, iic2d, jj2d, ii2d, tmp_buf, &
                faceptr, vstag2)
            idx = idx + facearea
        END IF

        IF (PRESENT(v3) .AND. exW) THEN
            CALL v3%get_buffer(faceptr, igrid, iface)

            CALL prolong1(jjc2d, iic2d, jj2d, ii2d, &
                recvbuf(offset+idx:offset+idx+facearea-1), tmp_buf, wstag1)
            CALL prolong2(jjc2d, iic2d, jj2d, ii2d, tmp_buf, &
                faceptr, wstag2)
            idx = idx + facearea
        END IF

        IF (PRESENT(s1)) THEN
            CALL s1%get_buffer(faceptr, igrid, iface)

            CALL prolong1(jjc2d, iic2d, jj2d, ii2d, &
                recvbuf(offset+idx:offset+idx+facearea-1), tmp_buf, 0)
            CALL prolong2(jjc2d, iic2d, jj2d, ii2d, tmp_buf, &
                faceptr, 0)
            idx = idx + facearea
        END IF

        IF (PRESENT(s2)) THEN
            CALL s2%get_buffer(faceptr, igrid, iface)

            CALL prolong1(jjc2d, iic2d, jj2d, ii2d, &
                recvbuf(offset+idx:offset+idx+facearea-1), tmp_buf, 0)
            CALL prolong2(jjc2d, iic2d, jj2d, ii2d, tmp_buf, &
                faceptr, 0)
            idx = idx + facearea
        END IF

        IF (PRESENT(s3)) THEN
            CALL s3%get_buffer(faceptr, igrid, iface)

            CALL prolong1(jjc2d, iic2d, jj2d, ii2d, &
                recvbuf(offset+idx:offset+idx+facearea-1), tmp_buf, 0)
            CALL prolong2(jjc2d, iic2d, jj2d, ii2d, tmp_buf, &
                faceptr, 0)
            idx = idx + facearea
        END IF

        ! Check that message length is calculated correctly
        IF (idx /= recvidxlist(2, recvid)) THEN
            WRITE(*, *) "idx:", idx, &
                "recvidxlist(2, recvid):", recvidxlist(2, recvid)
            CALL errr(__FILE__, __LINE__)
        END IF

        DEALLOCATE(tmp_buf)
    END SUBROUTINE read_buffer


    ! Process receive buffers as they arrive, wait for send
    ! buffers to be free
    SUBROUTINE process_bufs(v1, v2, v3, s1, s2, s3, normal)
        ! Subroutine arguments
        TYPE(field_t), OPTIONAL, INTENT(inout) :: &
            v1, v2, v3, s1, s2, s3
        LOGICAL, OPTIONAL, INTENT(in) :: normal

        INTEGER(int32) :: idx, i
        TYPE(MPI_Status) :: recvstatus
        INTEGER(int32) :: recvmessagelen
        INTEGER(int32) :: unpacklen

        DO WHILE (.TRUE.)
            IF (nrecv == 0) EXIT
            CALL MPI_Waitany(nrecv, recvreqs, idx, recvstatus)

            IF (idx /= MPI_UNDEFINED) THEN
                CALL MPI_Get_count(recvstatus, mglet_mpi_real, &
                    recvmessagelen)

                unpacklen = 0
                DO i = 1, irecv
                    IF (recvidxlist(1, i) == recvlist(idx) &
                            .AND. recvidxlist(2, i) > 0) THEN
                        CALL read_buffer(i, v1, v2, v3, s1, s2, s3, normal)
                        unpacklen = unpacklen + recvidxlist(2, i)
                    END IF
                END DO

                IF (recvmessagelen /= unpacklen) THEN
                    CALL errr(__FILE__, __LINE__)
                END IF
            ELSE
                EXIT
            END IF
        END DO
        CALL MPI_Waitall(nsend, sendreqs, MPI_STATUSES_IGNORE)
    END SUBROUTINE process_bufs


    ! Prolongates the first direction in a 2D field
    SUBROUTINE prolong1(jjc, iic, jj, ii, in, out, istag)
        ! Subroutine arguments
        INTEGER(intk), INTENT(IN) :: jjc, iic, jj, ii
        REAL(realk), INTENT(IN) :: in(jjc, iic)
        REAL(realk), INTENT(OUT) :: out(jj, iic)
        INTEGER(intk), INTENT(IN) :: istag

        ! Local variables
        INTEGER(intk) :: j, ic, jc

        ! Variable non-staggered in first dir.
        IF (istag == 0) THEN
            DO ic = 1, iic
                DO j = 1, jj, 2
                    jc = 2 + (j-1)/2
                    out(j, ic) = in(jc, ic)
                    out(j+1, ic) = in(jc, ic)
                END DO
            END DO

        ! Variable staggered in first dir.
        ELSE IF (istag == 1) THEN
            DO ic = 1, iic
                DO j = 1, jj, 2
                    jc = 2 + (j-1)/2
                    out(j, ic) = 0.5*(in(jc, ic) + in(jc-1, ic))
                    out(j+1, ic) = in(jc, ic)
                END DO
            END DO
        ELSE
            CALL errr(__FILE__, __LINE__)
        END IF
    END SUBROUTINE prolong1


    ! Prolongates the second direction in a 2D field
    SUBROUTINE prolong2(jjc, iic, jj, ii, in, out, istag)
        ! Subroutine arguments
        INTEGER(intk), INTENT(IN) :: jjc, iic, jj, ii
        REAL(realk), INTENT(IN) :: in(jj, iic)
        REAL(realk), INTENT(OUT) :: out(jj, ii)
        INTEGER(intk), INTENT(IN) :: istag

        ! Local variables
        INTEGER(intk) :: i, j, ic

        ! Variable non-staggered in second dir.
        IF (istag == 0) THEN
            DO i = 1, ii, 2
                ic = 2 + (i-1)/2
                DO j = 1, jj
                    out(j, i) = in(j, ic)
                    out(j, i+1) = in(j, ic)
                END DO
            END DO

        ! Variable staggered in second dir.
        ELSE IF (istag == 1) THEN
            DO i = 1, ii, 2
                ic = 2 + (i-1)/2
                DO j = 1, jj
                    out(j, i) = 0.5*(in(j, ic) + in(j, ic-1))
                    out(j, i+1) = in(j, ic)
                END DO
            END DO
        ELSE
            CALL errr(__FILE__, __LINE__)
        END IF
    END SUBROUTINE prolong2


    SUBROUTINE init_parent1()
        ! Subroutine arguments
        ! none...

        ! Local variables
        ! none...

        IF (is_init) CALL errr(__FILE__, __LINE__)

        ! Check if parent_core_mod provides necessary infrastructure
        IF (.NOT. is_parent_core_init) CALL errr(__FILE__, __LINE__)

        ALLOCATE(recvidxlist(3, irecv), source=0_intk)
        ALLOCATE(sendlist(isend), source=0_int32)
        ALLOCATE(recvlist(irecv), source=0_int32)
        ALLOCATE(sendreqs(isend), source=MPI_REQUEST_NULL)
        ALLOCATE(recvreqs(irecv), source=MPI_REQUEST_NULL)
        nrecv = 0
        nsend = 0

        is_init = .TRUE.
    END SUBROUTINE init_parent1


    SUBROUTINE finish_parent1()
        ! Subroutine arguments
        ! none...

        ! Local variables
        ! none...

        IF (.NOT. is_init) CALL errr(__FILE__, __LINE__)

        DEALLOCATE(sendreqs)
        DEALLOCATE(recvreqs)
        DEALLOCATE(sendlist)
        DEALLOCATE(recvlist)
        DEALLOCATE(recvidxlist)
        nrecv = 0
        nsend = 0

        is_init = .FALSE.
    END SUBROUTINE finish_parent1
END MODULE parent1_mod
