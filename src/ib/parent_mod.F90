MODULE parent_mod
    USE MPI_f08
    USE core_mod

    IMPLICIT NONE (type, external)
    PRIVATE

    ! Maximum number of connections on one single process, either
    ! outgoing or incomming, on any single grid level
    INTEGER(intk) :: maxConns

    ! Lists that hold the Send and Recv connections per grid level
    ! This list must be pre-compiled before the first call to 'connect'
    ! is being made. The reason for this is because it is expensive
    ! to compute every single time a connect is being made.
    !
    ! Dimensions contain:
    !   Dim 1: Information about a specific connection
    !   Dim 2: The different connections
    !
    ! The information in the first dimension is sorted as follows:
    !   Field 1: Rank of receiving process
    !   Field 2: Rank of sending process
    !   Field 3: ID of receiving grid
    !   Field 4: ID of sending grid
    !   Field 5: Which face (1..26) to receive
    !   Field 6: Which face (1..26) to send
    !   Field 7: Message tag (for MPI)
    INTEGER(intk), ALLOCATABLE :: sendConns(:, :), recvConns(:, :)

    ! Lists that hold the send and receive request arrays
    TYPE(MPI_Request), ALLOCATABLE :: sendReqs(:), recvReqs(:)

    ! Lists that hold the messages that are ACTUALLY sendt and received
    INTEGER(intk) :: nSend, nRecv, nRecvFaces
    INTEGER(int32), ALLOCATABLE :: sendList(:), recvList(:)
    INTEGER(intk), ALLOCATABLE :: recvIdxList(:,:)

    ! Number of send and receive connections
    INTEGER(intk) :: iSend = 0, iRecv = 0

    ! Counters for send- and receive operations (for locations in
    ! send and receive buffers)
    INTEGER(intk) :: sendCounter, recvCounter

    ! Number of variables per cell to exchange
    INTEGER(intk) :: nVars

    ! Accumulated message length
    INTEGER(int32) :: messageLength

    ! Variable to indicate if the connection information has
    ! been created.
    LOGICAL :: is_init = .FALSE.

    ! Pack also bu, bv, bw with vectors and use for prolongation
    ! TODO: implement!
    LOGICAL :: pack_buvw = .FALSE.

    ! Fields
    TYPE(field_t), POINTER :: u, v, w, p1, p2, p3

    ! If true, exchange only surface normal component of
    ! vector field
    LOGICAL :: sn

    PUBLIC :: parent, init_parent, finish_parent

CONTAINS

    ! Main parent function
    SUBROUTINE parent(ilevel, v1, v2, v3, s1, s2, s3, normal)

        ! Level to set PAR bc
        INTEGER(intk), INTENT(in) :: ilevel

        ! One or more fields
        TYPE(field_t), TARGET, OPTIONAL, INTENT(inout) :: &
            v1, v2, v3, s1, s2, s3

        ! Only exchange face-normal velocity components
        LOGICAL, OPTIONAL, INTENT(in) :: normal

        CALL start_timer(210)

        ! Check if the connection information has been created
        IF (is_init .EQV. .FALSE.) THEN
            WRITE(*,*) "'parent' not initialized"
            CALL errr(__FILE__, __LINE__)
        END IF

        ! Check that no other transfers are in progress
        IF (nSend > 0 .OR. nRecv > 0) THEN
            WRITE(*,*) "Other transfer in progress."
            CALL errr(__FILE__, __LINE__)
        END IF

        ! If one vector argument is given, check that all three is present.
        nVars = 0
        IF (PRESENT(v1)) THEN
            IF (.NOT. (PRESENT(v2) .AND. PRESENT(v3))) THEN
                WRITE(*,*) "If one vector arg is present, all three must be present."
                CALL errr(__FILE__, __LINE__)
            END IF
            u => v1
            v => v2
            w => v3
            nVars = nVars + 3
        ELSE IF (PRESENT(v2) .OR. PRESENT(v3)) THEN
            WRITE(*,*) "If one vector arg is present, all three must be present."
            CALL errr(__FILE__, __LINE__)
        END IF

        ! Set pointers to scalars
        IF (PRESENT(s1)) THEN
            p1 => s1
            nVars = nVars + 1
        END IF
        IF (PRESENT(s2)) THEN
            p2 => s2
            nVars = nVars + 1
        END IF
        IF (PRESENT(s3)) THEN
            p3 => s3
            nVars = nVars + 1
        END IF

        sn = .FALSE.
        IF (PRESENT(normal)) THEN
            IF (normal .eqv. .TRUE.) THEN
                IF (PRESENT(v1)) THEN
                    nVars = nVars - 2
                    sn = .TRUE.
                ELSE
                    WRITE(*,*) "normal=.TRUE. require a vector v1, v2, v3."
                    CALL errr(__FILE__, __LINE__)
                END IF
            END IF
        END IF

        ! Not specifying any fields would be very strange
        IF (nVars == 0) THEN
            WRITE(*,*) "You have not specified any fields to exchange."
            CALL errr(__FILE__, __LINE__)
        END IF

        CALL recv_all(ilevel)
        CALL send_all(ilevel)
        CALL process_bufs()

        ! Clear counters and unset pointers. This is important to avoid
        ! problems at next function call
        nRecv = 0
        nSend = 0
        nVars = 0
        sn = .FALSE.

        NULLIFY(u)
        NULLIFY(v)
        NULLIFY(w)
        NULLIFY(p1)
        NULLIFY(p2)
        NULLIFY(p3)

        CALL stop_timer(210)
    END SUBROUTINE parent


    ! Perform all Recv-calls
    SUBROUTINE recv_all(ilevel)
        ! Subroutine arguments
        INTEGER(intk), INTENT(in) :: ilevel

        ! Local variables
        INTEGER(intk) :: i, iprocnbr, igrid, iface, ilevelgrid, faceArea

        ! Post all receive calls
        recvCounter = 0
        messageLength = 0
        nRecv = 0
        recvIdxList = 0

        DO i = 1, iRecv
            iprocnbr      = recvConns(2, i)
            igrid         = recvConns(3, i)
            iface         = recvConns(5, i)
            ilevelgrid    = level(igrid)

            IF (ilevel == ilevelgrid) THEN
                faceArea      = face_area(igrid, iface)
                nRecvFaces    = nRecvFaces + 1
                recvIdxList(1, i) = iprocnbr
                recvIdxList(2, i) = nVars*faceArea
                recvIdxList(3, i) = recvCounter + messageLength
                messageLength = messageLength + nVars*faceArea

                IF (recvCounter + messageLength > idim_mg_bufs) THEN
                    CALL errr(__FILE__, __LINE__)
                END IF
            END IF

            IF (messageLength > 0) THEN
                IF (i == iRecv) THEN
                    CALL post_recv(iprocnbr)
                ELSE IF (recvConns(2, i + 1) /= iprocnbr) THEN
                    CALL post_recv(iprocnbr)
                END IF
            END IF
        END DO
    END SUBROUTINE recv_all


    ! Perform a single Recv
    SUBROUTINE post_recv(iprocnbr)
        ! Identifier of receive connection
        INTEGER(int32), INTENT(in) :: iprocnbr

        ! Local variables (for convenience)
        ! none...

        nRecv = nRecv + 1
        recvList(nRecv) = iprocnbr

        CALL MPI_Irecv(recvBuf(recvCounter+1) , messageLength, &
            mglet_mpi_real, iprocnbr, 1, MPI_COMM_WORLD, recvReqs(nRecv))

        recvCounter = recvCounter + messageLength
        messageLength = 0
    END SUBROUTINE post_recv


    ! Perform all send calls
    SUBROUTINE send_all(ilevel)
        INTEGER(intk), INTENT(in) :: ilevel
        INTEGER(intk) :: i, iprocnbr, igrid, ilevelgrid, iface, faceArea

        ! Pack all buffers and send data
        sendCounter = 0
        messageLength = 0
        nSend = 0

        DO i = 1, iSend
            iprocnbr      = sendConns(1, i)
            igrid         = sendConns(3, i)
            iface         = sendConns(5, i)
            ilevelgrid    = level(igrid)

            IF (ilevel == ilevelgrid) THEN
                faceArea      = face_area(igrid, iface)
                IF (sendCounter + messageLength + nVars*faceArea > idim_mg_bufs) THEN
                    CALL errr(__FILE__, __LINE__)
                END IF
                CALL write_buffer(i)
                messageLength = messageLength + nVars*faceArea
            END IF

            IF (messageLength > 0) THEN
                IF (i == iSend) THEN
                    CALL post_send(iprocnbr)
                ELSE IF (sendConns(1, i + 1) /= iprocnbr) THEN
                    CALL post_send(iprocnbr)
                END IF
            END IF
        END DO
    END SUBROUTINE send_all


    ! Perform a single send call
    SUBROUTINE post_send(iprocnbr)
        ! Subroutine arguments
        INTEGER(int32), INTENT(in) :: iprocnbr

        ! Local variables (for convenience)
        ! none...

        nSend = nSend + 1
        sendList(nSend) = iprocnbr

        CALL MPI_Isend(sendbuf(sendCounter + 1), messageLength, &
            mglet_mpi_real, iprocnbr, 1, MPI_COMM_WORLD, sendReqs(nSend))

        sendCounter = sendCounter + messageLength
        messageLength = 0
    END SUBROUTINE post_send


    ! Write Send buffers
    !
    ! Write the relevant fields into the send buffers
    SUBROUTINE write_buffer(sendId)
        ! Input parameter
        INTEGER(int32), INTENT(in) :: sendId

        ! Grid dimensions and pointers
        INTEGER(intk) :: kk, jj, ii
        INTEGER(intk) :: ip3

        ! Indices of start- and stop of iteration over boundary face
        INTEGER(intk) :: ista, isto, jsta, jsto, ksta, ksto

        ! Grid to send from
        INTEGER(intk) :: igrid, igridc, iface, icomp

        ! Message sizes
        ! Must be int32 because it iterface with MPI
        INTEGER(int32) :: thisMessageLength, faceArea
        INTEGER(int32) :: iCount, pos, offset

        ! Iterators
        INTEGER(intk) :: i, j, k
        LOGICAL :: exU, exV, exW

        TYPE(field_t), POINTER :: bu

        ! Set variables from send table - *fine* grid and face
        igrid = sendConns(3, sendId)
        igridc = sendConns(4, sendId)
        iface = sendConns(5, sendId)

        ! Get grid dimentsions and pointers for coarse grid
        ! (where data is fetched from)
        CALL get_mgdims(kk, jj, ii, igridc)

        ! Face area
        faceArea = face_area(igrid, iface)
        thisMessageLength = nVars*faceArea

        ! Check that buffer does not overflow
        IF (sendCounter + messageLength + thisMessageLength > idim_mg_bufs) THEN
            CALL errr(__FILE__, __LINE__)
        END IF

        ! Reset message size counter
        offset = sendCounter + messageLength
        iCount = 0

        ! Which vectors to exchange
        exU = (sn .AND. iface < 3) .OR. (.NOT. sn)
        exV = (sn .AND. (iface > 2 .AND. iface < 5)) .OR. (.NOT. sn)
        exW = (sn .AND. iface > 4) .OR. (.NOT. sn)

        ! Pack BU, BV, BW if requested
        IF (exU .AND. pack_buvw) THEN
            CALL get_field(bu, "BU")
            CALL bu%get_ip(ip3, igridc)
            icomp = 0
            SELECT CASE(iface)
                CASE (1, 2)
                    IF (bu%istag == 1) icomp = 1
                CASE (3, 4)
                    IF (bu%jstag == 1) icomp = 2
                CASE (5, 6)
                    IF (bu%kstag == 1) icomp = 3
            END SELECT

            CALL start_and_stop(igrid, iface, icomp, ista, isto, &
                jsta, jsto, ksta, ksto)
            DO i = ista, isto
                DO j = jsta, jsto
                    DO k = ksta, ksto
                        iCount = iCount + 1
                        pos = jj*kk*(i-1) + kk*(j-1) + (k-1)
                        sendbuf(offset + iCount) = bu%arr(ip3 + pos)
                    END DO
                END DO
            END DO
        END IF

        IF (exV .AND. pack_buvw) THEN
            CALL get_field(bu, "BV")
            CALL bu%get_ip(ip3, igridc)
            icomp = 0
            SELECT CASE(iface)
                CASE (1, 2)
                    IF (bu%istag == 1) icomp = 1
                CASE (3, 4)
                    IF (bu%jstag == 1) icomp = 2
                CASE (5, 6)
                    IF (bu%kstag == 1) icomp = 3
            END SELECT

            CALL start_and_stop(igrid, iface, icomp, ista, isto, &
                jsta, jsto, ksta, ksto)
            DO i = ista, isto
                DO j = jsta, jsto
                    DO k = ksta, ksto
                        iCount = iCount + 1
                        pos = jj*kk*(i-1) + kk*(j-1) + (k-1)
                        sendbuf(offset + iCount) = bu%arr(ip3 + pos)
                    END DO
                END DO
            END DO
        END IF

        IF (exW .AND. pack_buvw) THEN
            CALL get_field(bu, "BW")
            CALL bu%get_ip(ip3, igridc)
            icomp = 0
            SELECT CASE(iface)
                CASE (1, 2)
                    IF (bu%istag == 1) icomp = 1
                CASE (3, 4)
                    IF (bu%jstag == 1) icomp = 2
                CASE (5, 6)
                    IF (bu%kstag == 1) icomp = 3
            END SELECT

            CALL start_and_stop(igrid, iface, icomp, ista, isto, &
                jsta, jsto, ksta, ksto)
            DO i = ista, isto
                DO j = jsta, jsto
                    DO k = ksta, ksto
                        iCount = iCount + 1
                        pos = jj*kk*(i-1) + kk*(j-1) + (k-1)
                        sendbuf(offset + iCount) = bu%arr(ip3 + pos)
                    END DO
                END DO
            END DO
        END IF

        ! Fill buffers
        IF (ASSOCIATED(u) .AND. exU) THEN
            CALL u%get_ip(ip3, igridc)
            icomp = 0
            SELECT CASE(iface)
                CASE (1, 2)
                    IF (u%istag == 1) icomp = 1
                CASE (3, 4)
                    IF (u%jstag == 1) icomp = 2
                CASE (5, 6)
                    IF (u%kstag == 1) icomp = 3
            END SELECT

            CALL start_and_stop(igrid, iface, icomp, ista, isto, &
                jsta, jsto, ksta, ksto)
            DO i = ista, isto
                DO j = jsta, jsto
                    DO k = ksta, ksto
                        iCount = iCount + 1
                        pos = jj*kk*(i-1) + kk*(j-1) + (k-1)
                        sendbuf(offset + iCount) = u%arr(ip3 + pos)
                    END DO
                END DO
            END DO
        END IF

        IF (ASSOCIATED(v) .AND. exV) THEN
            CALL v%get_ip(ip3, igridc)
            icomp = 0
            SELECT CASE(iface)
                CASE (1, 2)
                    IF (v%istag == 1) icomp = 1
                CASE (3, 4)
                    IF (v%jstag == 1) icomp = 2
                CASE (5, 6)
                    IF (v%kstag == 1) icomp = 3
            END SELECT

            CALL start_and_stop(igrid, iface, icomp, ista, isto, &
                jsta, jsto, ksta, ksto)
            DO i = ista, isto
                DO j = jsta, jsto
                    DO k = ksta, ksto
                        iCount = iCount + 1
                        pos = jj*kk*(i-1) + kk*(j-1) + (k-1)
                        sendbuf(offset + iCount) = v%arr(ip3 + pos)
                    END DO
                END DO
            END DO
        END IF

        IF (ASSOCIATED(w) .AND. exW) THEN
            CALL w%get_ip(ip3, igridc)
            icomp = 0
            SELECT CASE(iface)
                CASE (1, 2)
                    IF (w%istag == 1) icomp = 1
                CASE (3, 4)
                    IF (w%jstag == 1) icomp = 2
                CASE (5, 6)
                    IF (w%kstag == 1) icomp = 3
            END SELECT

            CALL start_and_stop(igrid, iface, icomp, ista, isto, &
                jsta, jsto, ksta, ksto)
            DO i = ista, isto
                DO j = jsta, jsto
                    DO k = ksta, ksto
                        iCount = iCount + 1
                        pos = jj*kk*(i-1) + kk*(j-1) + (k-1)
                        sendbuf(offset + iCount) = w%arr(ip3 + pos)
                    END DO
                END DO
            END DO
        END IF

        IF (ASSOCIATED(p1)) THEN
            CALL p1%get_ip(ip3, igridc)
            icomp = 0
            SELECT CASE(iface)
                CASE (1, 2)
                    IF (p1%istag == 1) icomp = 1
                CASE (3, 4)
                    IF (p1%jstag == 1) icomp = 2
                CASE (5, 6)
                    IF (p1%kstag == 1) icomp = 3
            END SELECT

            CALL start_and_stop(igrid, iface, icomp, ista, isto, &
                jsta, jsto, ksta, ksto)
            DO i = ista, isto
                DO j = jsta, jsto
                    DO k = ksta, ksto
                        iCount = iCount + 1
                        pos = jj*kk*(i-1) + kk*(j-1) + (k-1)
                        sendbuf(offset + iCount) = p1%arr(ip3 + pos)
                    END DO
                END DO
            END DO
        END IF

        IF (ASSOCIATED(p2)) THEN
            CALL p2%get_ip(ip3, igridc)
            icomp = 0
            SELECT CASE(iface)
                CASE (1, 2)
                    IF (p2%istag == 1) icomp = 1
                CASE (3, 4)
                    IF (p2%jstag == 1) icomp = 2
                CASE (5, 6)
                    IF (p2%kstag == 1) icomp = 3
            END SELECT

            CALL start_and_stop(igrid, iface, icomp, ista, isto, &
                jsta, jsto, ksta, ksto)
            DO i = ista, isto
                DO j = jsta, jsto
                    DO k = ksta, ksto
                        iCount = iCount + 1
                        pos = jj*kk*(i-1) + kk*(j-1) + (k-1)
                        sendbuf(offset + iCount) = p2%arr(ip3 + pos)
                    END DO
                END DO
            END DO
        END IF

        IF (ASSOCIATED(p3)) THEN
            CALL p3%get_ip(ip3, igridc)
            icomp = 0
            SELECT CASE(iface)
                CASE (1, 2)
                    IF (p3%istag == 1) icomp = 1
                CASE (3, 4)
                    IF (p3%jstag == 1) icomp = 2
                CASE (5, 6)
                    IF (p3%kstag == 1) icomp = 3
            END SELECT

            CALL start_and_stop(igrid, iface, icomp, ista, isto, &
                jsta, jsto, ksta, ksto)
            DO i = ista, isto
                DO j = jsta, jsto
                    DO k = ksta, ksto
                        iCount = iCount + 1
                        pos = jj*kk*(i-1) + kk*(j-1) + (k-1)
                        sendbuf(offset + iCount) = p3%arr(ip3 + pos)
                    END DO
                END DO
            END DO
        END IF

        ! Check that message length was calculated correctly
        IF (thisMessageLength /= iCount) THEN
            write(*,*) "nVars:", nVars, " thisMessageLength:", &
                thisMessageLength, " iCount:", iCount
            CALL errr(__FILE__, __LINE__)
        END IF
    END SUBROUTINE write_buffer


    ! Read Receive buffers
    !
    ! Write the contents of the receive buffers back in their
    ! matching fields
    SUBROUTINE read_buffer(recvId)

        ! Input parameter
        INTEGER(intk), INTENT(in) :: recvId

        ! Grid dimensions and pointers
        INTEGER(intk) :: kk, jj, ii
        INTEGER(intk) :: kkc, jjc, iic
        INTEGER(intk) :: idx
        INTEGER(intk) :: jj2d, ii2d, jjc2d, iic2d

        ! Grid to send from
        ! Must be intk because it intreface with MGLET
        INTEGER(intk) :: igrid, igridc, iface

        ! Message sizes
        ! Must be int32 because it iterface with MPI
        INTEGER(int32) :: offset

        INTEGER(intk) :: nelem, tmp_buf_size, bu_len
        INTEGER(intk) :: ustag1, ustag2, vstag1, vstag2, wstag1, wstag2
        REAL(realk), POINTER, CONTIGUOUS :: faceptr(:, :, :)
        REAL(realk), POINTER, CONTIGUOUS :: bu_ptr(:, :), bv_ptr(:, :), &
            bw_ptr(:, :)
        REAL(realk), ALLOCATABLE :: tmp_buf(:)
        LOGICAL :: exU, exV, exW

        ! Set variables from send table
        igrid = recvConns(3, recvId)
        igridc = recvConns(4, recvId)
        iface = recvConns(5, recvId)

        ! Get grid dimentsions and pointers
        CALL get_mgdims(kk, jj, ii, igrid)

        ! Get start- and stop indices of grid
        CALL idx2d(kk, jj, ii, iface, kkc, jjc, iic, jj2d, ii2d, jjc2d, iic2d)
        CALL stag(iface, ustag1, ustag2, vstag1, vstag2, wstag1, wstag2)
        nelem = face_area(igrid, iface)

        ! Allocate temporary work buffer
        tmp_buf_size = MAX(ii*jj, ii*kk, jj*kk)
        ALLOCATE(tmp_buf(tmp_buf_size))
        tmp_buf = 0.0

        ! Offset in receive buffer
        offset = recvIdxList(3, recvId) + 1
        idx = 0

        exU = (sn .AND. iface < 3) .OR. (.NOT. sn)
        exV = (sn .AND. (iface > 2 .AND. iface < 5)) .OR. (.NOT. sn)
        exW = (sn .AND. iface > 4) .OR. (.NOT. sn)

        ! Pack BU, BV, BW if requested
        NULLIFY(bu_ptr)
        IF (exU .AND. pack_buvw) THEN
            bu_len = jjc2d*iic2d
            bu_ptr(1:jjc2d, 1:iic2d) => recvBuf(offset+idx:offset+idx+bu_len)
            idx = idx + nelem
        END IF

        NULLIFY(bv_ptr)
        IF (exV .AND. pack_buvw) THEN
            bu_len = jjc2d*iic2d
            bv_ptr(1:jjc2d, 1:iic2d) => recvBuf(offset+idx:offset+idx+bu_len)
            idx = idx + nelem
        END IF

        NULLIFY(bw_ptr)
        IF (exW .AND. pack_buvw) THEN
            bu_len = jjc2d*iic2d
            bw_ptr(1:jjc2d, 1:iic2d) => recvBuf(offset+idx:offset+idx+bu_len)
            idx = idx + nelem
        END IF


        IF (ASSOCIATED(u) .AND. exU) THEN
            CALL u%buffers%get_buffer(faceptr, igrid, iface)

            CALL prolong1(jjc2d, iic2d, jj2d, ii2d, &
                recvBuf(offset+idx:offset+idx+nelem), tmp_buf, ustag1, &
                bu_ptr)
            CALL prolong2(jjc2d, iic2d, jj2d, ii2d, tmp_buf, &
                faceptr, ustag2, bu_ptr)
            idx = idx + nelem
        END IF


        IF (ASSOCIATED(v) .AND. exV) THEN
            CALL v%buffers%get_buffer(faceptr, igrid, iface)

            CALL prolong1(jjc2d, iic2d, jj2d, ii2d, &
                recvBuf(offset+idx:offset+idx+nelem), tmp_buf, vstag1, &
                bv_ptr)
            CALL prolong2(jjc2d, iic2d, jj2d, ii2d, tmp_buf, &
                faceptr, vstag2, bv_ptr)
            idx = idx + nelem
        END IF


        IF (ASSOCIATED(w) .AND. exW) THEN
            CALL w%buffers%get_buffer(faceptr, igrid, iface)

            CALL prolong1(jjc2d, iic2d, jj2d, ii2d, &
                recvBuf(offset+idx:offset+idx+nelem), tmp_buf, wstag1, &
                bw_ptr)
            CALL prolong2(jjc2d, iic2d, jj2d, ii2d, tmp_buf, &
                faceptr, wstag2, bw_ptr)
            idx = idx + nelem
        END IF

        IF (ASSOCIATED(p1)) THEN
            CALL p1%buffers%get_buffer(faceptr, igrid, iface)

            CALL prolong1(jjc2d, iic2d, jj2d, ii2d, &
                recvBuf(offset+idx:offset+idx+nelem), tmp_buf, 0)
            CALL prolong2(jjc2d, iic2d, jj2d, ii2d, tmp_buf, &
                faceptr, 0)
            idx = idx + nelem
        END IF

        IF (ASSOCIATED(p2)) THEN
            CALL p2%buffers%get_buffer(faceptr, igrid, iface)

            CALL prolong1(jjc2d, iic2d, jj2d, ii2d, &
                recvBuf(offset+idx:offset+idx+nelem), tmp_buf, 0)
            CALL prolong2(jjc2d, iic2d, jj2d, ii2d, tmp_buf, &
                faceptr, 0)
            idx = idx + nelem
        END IF

        IF (ASSOCIATED(p3)) THEN
            CALL p3%buffers%get_buffer(faceptr, igrid, iface)

            CALL prolong1(jjc2d, iic2d, jj2d, ii2d, &
                recvBuf(offset+idx:offset+idx+nelem), tmp_buf, 0)
            CALL prolong2(jjc2d, iic2d, jj2d, ii2d, tmp_buf, &
                faceptr, 0)
            idx = idx + nelem
        END IF

        IF (ASSOCIATED(bu_ptr)) THEN
            NULLIFY(bu_ptr)
        END IF

        IF (ASSOCIATED(bv_ptr)) THEN
            NULLIFY(bv_ptr)
        ENDIF

        IF (ASSOCIATED(bw_ptr)) THEN
            NULLIFY(bw_ptr)
        ENDIF

        ! Check that message length is calculated correctly
        IF (idx /= recvIdxList(2, recvId)) THEN
            WRITE(*,*) "idx:", idx, &
                "recvIdxList(2, recvId):", recvIdxList(2, recvId)
            CALL errr(__FILE__, __LINE__)
        END IF

        DEALLOCATE(tmp_buf)
    END SUBROUTINE read_buffer


    ! Process receive buffers as they arrive, wait for send
    ! buffers to be free
    SUBROUTINE process_bufs()

        INTEGER(int32) :: idx, i
        TYPE(MPI_Status) :: recvstatus
        INTEGER(int32) :: recvMessageLen
        INTEGER(int32) :: unpackLen

        DO WHILE (.TRUE.)
            CALL MPI_Waitany(nRecv, recvReqs, idx, recvstatus)

            IF (idx /= MPI_UNDEFINED) THEN
                CALL MPI_Get_count(recvstatus, mglet_mpi_real, &
                    recvMessageLen)

                unpackLen = 0
                DO i = 1, iRecv
                    IF (recvIdxList(1, i) == recvList(idx) .AND. recvIdxList(2, i) > 0) THEN
                        CALL read_buffer(i)
                        unpackLen = unpackLen + recvIdxList(2, i)
                    END IF
                END DO

                IF (recvMessageLen /= unpackLen) THEN
                    CALL errr(__FILE__, __LINE__)
                END IF
            ELSE
                EXIT
            END IF
        END DO
        CALL MPI_Waitall(nSend, sendReqs, MPI_STATUSES_IGNORE)
    END SUBROUTINE process_bufs


    SUBROUTINE init_parent()
        INTEGER(intk) :: i, iface, igrid, inbr, iprocnbr, itypbc

        INTEGER(int32), ALLOCATABLE :: maxTag(:)
        INTEGER(int32), ALLOCATABLE :: sendcounts(:), sdispls(:)
        INTEGER(int32), ALLOCATABLE :: recvcounts(:), rdispls(:)

        INTEGER(int32) :: ierr

        CALL set_timer(210, "PARENT")

        ! Maximum number of parents for "simple" cases is number
        ! of grids*6. However, due to the possible prescence of
        ! strange grid structures we add a few more
        maxConns = INT((nMyGrids+1.0)*6.0*1.2)
        ALLOCATE(recvConns(7, maxConns))
        recvConns = 0

        ! The maximum number of concurrent communications are the number
        ! of processes
        ALLOCATE(recvIdxList(3, maxConns))
        ALLOCATE(sendList(numprocs))
        ALLOCATE(recvList(numprocs))
        ALLOCATE(sendReqs(numprocs))
        ALLOCATE(recvReqs(numprocs))
        recvIdxList = 0
        sendList = 0
        recvList = 0

        ALLOCATE(maxTag(0:numprocs-1))
        ALLOCATE(sendcounts(0:numprocs-1))
        ALLOCATE(sdispls(0:numprocs-1))
        ALLOCATE(recvcounts(0:numprocs-1))
        ALLOCATE(rdispls(0:numprocs-1))
        maxTag = 0
        sendcounts = 0
        sdispls = 0
        recvcounts = 0
        rdispls = 0

        nRecv = 0

        DO i = 1, nMyGrids
            igrid = myGrids(i)

            ! Loop over the boundary faces 1..6
            !
            ! Meaning of iterator:
            !   1 : FRONT  ( low X)
            !   2 : BACK   (high X)
            !   3 : RIGHT  ( low Y)
            !   4 : LEFT   (high Y)
            !   5 : BOTTOM ( low Z)
            !   6 : TOP    (high Z)
            !
            ! See also setcobone.F
            DO iface = 1, 6
                ! Get type of BC (assuming PAR is ibocond = 1)
                itypbc = itypboconds(1, iface, igrid)

                ! See also setcobone.F
                IF (itypbc == 8) THEN
                    ! Get neighbouring grid and process

                    inbr = iparent(igrid)
                    iprocnbr = idprocofgrd(inbr)

                    nRecv = nRecv + 1

                    IF (nRecv > maxConns) THEN
                        write(*,*) "Number of PAR's exceeded on process ", myid
                        write(*,*) "maxConns =", maxConns, &
                            "nMyGrids =", nMyGrids, "nRecv = ", nRecv
                        CALL errr(__FILE__, __LINE__)
                    END IF

                    maxTag(iprocnbr) = maxTag(iprocnbr) + 1

                    recvConns(1, nRecv) = myid      ! Receiving process (this process)
                    recvConns(2, nRecv) = iprocnbr  ! Sending process (neighbour process)
                    recvConns(3, nRecv) = igrid     ! Receiving grid (on current process)
                    recvConns(4, nRecv) = inbr      ! Sending grid (on neighbour process)
                    recvConns(5, nRecv) = iface     ! Which face receive (1..6)
                    recvConns(6, nRecv) = -1.0      ! unused atm.
                    recvConns(7, nRecv) = maxTag(iprocnbr)  ! Message tag

                    sendcounts(iprocnbr) = sendcounts(iprocnbr) + SIZE(recvConns, 1)
                END IF
            END DO
        END DO

        iRecv = nRecv

        ! Sort recvConns by process ID
        CALL sort_conns(recvConns(:,1:nRecv))

        ! Calculate sdispl offset
        DO i = 1,numprocs-1
            sdispls(i) = sdispls(i-1) + sendcounts(i-1)
        END DO

        ! First exchange NUMBER OF ELEMENTS TO RECEIVE, to be able to
        ! calculate rdispls array
        CALL MPI_Alltoall(sendcounts, 1, MPI_INTEGER, recvcounts, 1, &
            MPI_INTEGER, MPI_COMM_WORLD, ierr)

        ! Calculate rdispl offset
        DO i=1,numprocs-1
            rdispls(i) = rdispls(i-1) + recvcounts(i-1)
        END DO

        ! Allocate sendConns array
        iSend = (rdispls(numprocs-1) + recvcounts(numprocs-1))/SIZE(recvConns, 1)
        ALLOCATE(sendConns(7, iSend))
        sendConns = 0

        ! Exchange connection information
        CALL MPI_Alltoallv(recvConns, sendcounts, sdispls, MPI_INTEGER, &
            sendConns, recvcounts, rdispls, MPI_INTEGER, &
            MPI_COMM_WORLD)

        is_init = .TRUE.

        nRecv = 0

        ! Nullify pointers
        nullify(u)
        nullify(v)
        nullify(w)
        nullify(p1)
        nullify(p2)
        nullify(p3)
    END SUBROUTINE init_parent


    SUBROUTINE finish_parent()
        DEALLOCATE(sendConns)
        DEALLOCATE(recvConns)
        DEALLOCATE(sendReqs)
        DEALLOCATE(recvReqs)
        DEALLOCATE(sendList)
        DEALLOCATE(recvList)
        DEALLOCATE(recvIdxList)
        is_init = .FALSE.
    END SUBROUTINE finish_parent


    SUBROUTINE sort_conns(list)
        ! Input array to be sorted
        INTEGER(int32), INTENT(inout) :: list(:,:)

        INTEGER(intk) :: i,j

        ! Temporary storage
        INTEGER(int32) :: temp(7)

        IF (SIZE(list, 1) /= SIZE(temp)) THEN
            CALL errr(__FILE__, __LINE__)
        END IF

        ! Sort by sending processor number (field 2)
        DO i = 2,SIZE(list, 2)
            j = i - 1
            temp(:) = list(:,i)
            DO WHILE (j >= 1)
                IF (list(2,j) > temp(2)) THEN
                    list(:,j+1) = list(:,j)
                    j = j - 1
                ELSE
                    EXIT
                END IF
            END DO
            list(:,j+1) = temp(:)
        END DO

    END SUBROUTINE sort_conns


    SUBROUTINE idx2d(kk, jj, ii, iface, kkc, jjc, iic, jj2d, ii2d, jjc2d, iic2d)
        ! Subroutine arguments
        INTEGER(intk), INTENT(in) :: kk, jj, ii, iface
        INTEGER(intk), INTENT(out) :: kkc, jjc, iic
        INTEGER(intk), INTENT(out) :: jj2d, ii2d, jjc2d, iic2d

        kkc = kk/2 + 2
        jjc = jj/2 + 2
        iic = ii/2 + 2

        IF (iface == 1 .OR. iface == 2) THEN
            ii2d = jj
            jj2d = kk
            iic2d = jjc
            jjc2d = kkc
        ELSE IF (iface == 3 .OR. iface == 4) THEN
            ii2d = ii
            jj2d = kk
            iic2d = iic
            jjc2d = kkc
        ELSE IF (iface == 5 .OR. iface == 6) THEN
            ii2d = ii
            jj2d = jj
            iic2d = iic
            jjc2d = jjc
        ELSE
            CALL errr(__FILE__, __LINE__)
        END IF
    END SUBROUTINE idx2d


    SUBROUTINE stag(iface, ustag1, ustag2, vstag1, vstag2, wstag1, wstag2)
        ! Subroutine arguments
        INTEGER(intk), INTENT(in) :: iface
        INTEGER(intk), INTENT(out) :: ustag1, ustag2, vstag1, vstag2, &
            wstag1, wstag2

        ustag1 = 0
        ustag2 = 0
        vstag1 = 0
        vstag2 = 0
        wstag1 = 0
        wstag2 = 0

        IF (iface == 1 .OR. iface == 2) THEN
            vstag2 = 1
            wstag1 = 1
        ELSE IF (iface == 3 .OR. iface == 4) THEN
            ustag2 = 1
            wstag1 = 1
        ELSE IF (iface == 5 .OR. iface == 6) THEN
            ustag2 = 1
            vstag1 = 1
        ELSE
            CALL errr(__FILE__, __LINE__)
        END IF
    END SUBROUTINE stag


    SUBROUTINE start_and_stop(igrid, iface, icomp, &
        istart, istop, jstart, jstop, kstart, kstop)

        ! Subroutine arguments
        ! Igrid and iface of *fine* grid is to be given
        INTEGER(intk), INTENT(in) :: igrid, iface, icomp

        ! Returns start- and stop indices of the corresponding *coarse* grid
        ! region to be copied
        INTEGER(intk), INTENT(out) :: istart, istop, jstart, jstop, &
            kstart, kstop

        ! Local variables
        INTEGER(intk) :: kk, jj, ii
        INTEGER(intk) :: kpos, jpos, ipos

        ! Get diensions of (fine)grid
        CALL get_mgdims(kk, jj, ii, igrid)
        kpos = kposition(igrid)
        jpos = jposition(igrid)
        ipos = iposition(igrid)

        ! Select entire slab in coarse grid
        kstart = kpos - 2
        jstart = jpos - 2
        istart = ipos - 2

        kstop = kstart + kk/2 + 1
        jstop = jstart + jj/2 + 1
        istop = istart + ii/2 + 1

        ! Reduce selection to only contain face with 1 or 2 slices
        IF (iface == 1) THEN
            istart = ipos - 1
            istop = ipos - 1
        ELSE IF (iface == 2) THEN
            istart = ipos + (ii-4)/2
            IF (icomp == 1) istart = istart - 1
            istop = istart
        ELSE IF (iface == 3) THEN
            jstart = jpos - 1
            jstop = jpos - 1
        ELSE IF (iface == 4) THEN
            jstart = jpos + (jj-4)/2
            IF (icomp == 2) jstart = jstart - 1
            jstop = jstart
        ELSE IF (iface == 5) THEN
            kstart = kpos - 1
            kstop = kpos - 1
        ELSE IF (iface == 6) THEN
            kstart = kpos + (kk-4)/2
            IF (icomp == 3) kstart = kstart - 1
            kstop = kstart
        END IF
    END SUBROUTINE start_and_stop


    ! Calculate the area of a boundary face, i.e. the number of cells
    ! to be exchanged (all planes). Used to calculate message lengths.
    ! Needs to be given *fine grid* igrid and iface
    FUNCTION face_area(igrid, iface) RESULT(area)

        ! Result = length of message to be passed
        INTEGER(intk) :: area

        ! Input parameters, grid and boundary face information
        ! Must be intk because they are intreface to other MGLET functions
        INTEGER(intk), INTENT(in) :: igrid, iface

        ! Indices of start- and stop of iteration over boundary face
        ! Must be intk because they are intreface to other MGLET functions
        INTEGER(intk) :: istart, istop, jstart, jstop, kstart, kstop

        CALL start_and_stop(igrid, iface, 0, istart, istop, &
            jstart, jstop, kstart, kstop)

        area = (istop-istart+1)*(jstop-jstart+1)*(kstop-kstart+1)

        RETURN
    END FUNCTION face_area


    ! Prolongates the first direction in a 2D field
    SUBROUTINE prolong1(jjc, iic, jj, ii, in, out, istag, hilf)
        ! Subroutine arguments
        INTEGER(intk), INTENT(IN) :: jjc, iic, jj, ii
        REAL(realk), INTENT(IN) :: in(jjc, iic)
        REAL(realk), INTENT(OUT) :: out(jj, iic)
        INTEGER(intk), INTENT(IN) :: istag
        REAL(realk), INTENT(IN), OPTIONAL, POINTER :: hilf(:, :)

        ! Local variables
        INTEGER(intk) :: jf, ic, jc
        REAL(realk) :: bvcs, bvcn

        IF (pack_buvw .AND. istag == 1) THEN
            IF (.NOT. PRESENT(hilf)) THEN
                CALL errr(__FILE__, __LINE__)
            ELSE IF (.NOT. ASSOCIATED(hilf)) THEN
                CALL errr(__FILE__, __LINE__)
            END IF
        END IF

        ! Variable non-staggered in first dir.
        IF (istag == 0) THEN
            DO ic = 1, iic
                DO jf = 1, jj, 2
                    jc = 2 + (jf-1)/2
                    out(jf, ic) = in(jc, ic)
                    out(jf+1, ic) = in(jc, ic)
                END DO
            END DO

        ! Variable staggered in first dir.
        ELSE IF (istag == 1) THEN
            DO ic = 1, iic
                DO jf = 1, jj, 2
                    jc = 2 + (jf-1)/2

                    IF (pack_buvw) THEN
                        bvcn = hilf(jc, ic)*(1.0 - hilf(jc-1, ic))
                        bvcs = hilf(jc-1, ic)*(1.0 - hilf(jc ,ic))

                        out(jf, ic) = 0.5*(in(jc, ic) + in(jc-1, ic)) &
                            *(1.0 - bvcn - bvcs) &
                            + (in(jc, ic) + in(jc-1, ic))*(bvcn + bvcs)

                        out(jf+1, ic) = in(jc, ic)*hilf(jc, ic) &
                            + divide0(in(jc-1, ic) + in(jc+1,ic), &
                                      hilf(jc-1, ic) + hilf(jc+1, ic)) &
                            *(1.0 - hilf(jc, ic))
                    ELSE
                        out(jf, ic) = 0.5*(in(jc, ic) + in(jc-1, ic))
                        out(jf+1, ic) = in(jc, ic)
                    END IF

                END DO
            END DO
        ELSE
            CALL errr(__FILE__, __LINE__)
        END IF
    END SUBROUTINE prolong1


    ! Prolongates the second direction in a 2D field
    SUBROUTINE prolong2(jjc, iic, jj, ii, in, out, istag, hilf)
        ! Subroutine arguments
        INTEGER(intk), INTENT(IN) :: jjc, iic, jj, ii
        REAL(realk), INTENT(IN) :: in(jj, iic)
        REAL(realk), INTENT(OUT) :: out(jj, ii)
        INTEGER(intk), INTENT(IN) :: istag
        REAL(realk), INTENT(IN), OPTIONAL, POINTER :: hilf(:, :)

        ! Local variables
        INTEGER(intk) :: if, jf, ic, jc
        REAL(realk) :: bvcs, bvcn

        IF (pack_buvw .AND. istag == 1) THEN
            IF (.NOT. PRESENT(hilf)) THEN
                CALL errr(__FILE__, __LINE__)
            ELSE IF (.NOT. ASSOCIATED(hilf)) THEN
                CALL errr(__FILE__, __LINE__)
            END IF
        END IF

        ! Variable non-staggered in second dir.
        IF (istag == 0) THEN
            DO if = 1, ii, 2
                ic = 2 + (if-1)/2
                DO jf = 1, jj
                    out(jf, if  ) = in(jf, ic)
                    out(jf, if+1) = in(jf, ic)
                END DO
            END DO

        ! Variable staggered in second dir.
        ELSE IF (istag == 1) THEN
            DO if = 1, ii, 2
                ic = 2 + (if-1)/2
                DO jf = 1, jj
                    jc = 2 + (jf-1)/2

                    IF (pack_buvw) THEN
                        bvcn = hilf(jc, ic)*(1. - hilf(jc, ic-1))
                        bvcs = hilf(jc, ic-1)*(1. - hilf(jc, ic))

                        out(jf, if) = 0.5*(in(jf, ic) + in(jf, ic-1)) &
                            *(1.0 - bvcn - bvcs) &
                            + (in(jf, ic) + in(jf, ic-1))*(bvcn + bvcs)

                        out(jf, if+1) = in(jf, ic)*hilf(jc, ic) &
                            + divide0(in(jf, ic-1) + in(jf, ic+1), &
                                      hilf(jc, ic-1)+hilf(jc, ic+1)) &
                            * (1.0 - hilf(jc, ic))
                    ELSE
                        out(jf, if) = 0.5*(in(jf, ic) + in(jf, ic-1))
                        out(jf, if+1) = in(jf, ic)
                    END IF
                END DO
            END DO
        ELSE
            CALL errr(__FILE__, __LINE__)
        END IF
    END SUBROUTINE prolong2
END MODULE parent_mod
