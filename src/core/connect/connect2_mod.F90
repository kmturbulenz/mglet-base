MODULE connect2_mod
    USE precision_mod
    USE MPI_f08
    USE commbuf_mod, ONLY: sendbuf, recvbuf, isendbuf, irecvbuf
    USE err_mod, ONLY: errr
    USE timer_mod, ONLY: start_timer, set_timer, stop_timer
    ! USE pointers_mod, ONLY: get_ip3
    USE grids_mod, ONLY: mygrids, nmygrids, level, idprocofgrd, itypboconds, &
        maxlevel, minlevel, get_neighbours, get_mgdims
    USE comms_mod, ONLY: myid, numprocs
    USE field_mod
    USE qsort_mod, ONLY: sort_conns
    USE connect_core_mod

    IMPLICIT NONE (type, external)
    PRIVATE

    ! Lists that hold the send and receive request arrays
    TYPE(MPI_Request), ALLOCATABLE :: sendReqs(:), recvReqs(:)

    ! Lists that hold the messages that are ACTUALLY sendt and received
    INTEGER(intk) :: nSend, nRecv, nRecvFaces
    INTEGER(int32), ALLOCATABLE :: sendList(:), recvList(:)
    INTEGER(intk), ALLOCATABLE :: recvIdxList(:, :)

    ! Counters for send- and receive operations (for locations in
    ! send and receive buffers)
    INTEGER(intk) :: sendCounter, recvCounter

    ! Number of variables per cell to exchange
    ! Number of planes to exchange
    INTEGER(intk) :: nVars, nplane

    ! Exchange faces with geom only
    LOGICAL :: geometry

    ! Forward-flag. If 1, only send data forward (in positive
    ! coordinate direction) and if -1 only send data backward
    ! NOt to be used with corners
    INTEGER(intk) :: fwd

    ! If true, exchange only surface normal component of
    ! vector field
    LOGICAL :: sn

    ! Exchange corner lines/points as well
    LOGICAL :: vertices

    ! Flag to indicate very special behaviour
    CHARACTER(len=1) :: flag

    ! Accumulated message length
    INTEGER(int32) :: messageLength

    ! Minimum and maximul level to connect
    INTEGER(intk) :: minConLvl, maxConLvl

    ! Variable to indicate if the connection information has
    ! been created.
    LOGICAL :: isinit = .FALSE.

    ! A connect can be either of either REAL or INTEGER, never both in the
    ! same call. When connect_integer is .TRUE. we use isendbuf, irecvbuf
    ! instead of sendbuf and recvbuf
    LOGICAL :: connect_integer = .FALSE.

    ! Fields
    CLASS(basefield_t), POINTER :: u => NULL(), v => NULL(), w => NULL(), &
        p1 => NULL(), p2 => NULL(), p3 => NULL()

    PUBLIC :: connect, connect_int, init_connect2, finish_connect2

CONTAINS

    SUBROUTINE connect(ilevel, layers, v1, v2, v3, &
            s1, s2, s3, geom, corners, normal, forward, ityp, minlvl, maxlvl)

        ! Subroutine arguments
        INTEGER(intk), INTENT(in), OPTIONAL :: ilevel, layers
        TYPE(field_t), TARGET, OPTIONAL, INTENT(inout) :: &
            v1, v2, v3, s1, s2, s3
        LOGICAL, OPTIONAL, INTENT(in) :: geom, corners, normal
        INTEGER(intk), OPTIONAL, INTENT(in) :: forward
        CHARACTER(len=1), OPTIONAL, INTENT(in) :: ityp
        INTEGER(intk), INTENT(in), OPTIONAL :: minlvl, maxlvl

        ! Local variables
        LOGICAL :: has_v1, has_v2, has_v3, has_s1, has_s2, has_s3

        has_v1 = .FALSE.
        has_v2 = .FALSE.
        has_v3 = .FALSE.
        has_s1 = .FALSE.
        has_s2 = .FALSE.
        has_s3 = .FALSE.

        NULLIFY(u)
        NULLIFY(v)
        NULLIFY(w)
        NULLIFY(p1)
        NULLIFY(p2)
        NULLIFY(p3)

        IF (PRESENT(v1)) THEN
            has_v1 = .TRUE.
            u => v1
        END IF

        IF (PRESENT(v2)) THEN
            has_v2 = .TRUE.
            v => v2
        END IF

        IF (PRESENT(v3)) THEN
            has_v3 = .TRUE.
            w => v3
        END IF

        IF (PRESENT(s1)) THEN
            has_s1 = .TRUE.
            p1 => s1
        END IF

        IF (PRESENT(s2)) THEN
            has_s2 = .TRUE.
            p2 => s2
        END IF

        IF (PRESENT(s3)) THEN
            has_s3 = .TRUE.
            p3 => s3
        END IF

        connect_integer = .FALSE.
        CALL connect_impl(ilevel, layers, has_v1, has_v2, has_v3, &
            has_s1, has_s2, has_s3, geom, corners, normal, forward, ityp, &
            minlvl, maxlvl)
    END SUBROUTINE connect


    SUBROUTINE connect_int(ilevel, layers, v1, v2, v3, &
            s1, s2, s3, geom, corners, normal, forward, ityp, minlvl, maxlvl)

        ! Subroutine arguments
        INTEGER(intk), INTENT(in), OPTIONAL :: ilevel, layers
        TYPE(intfield_t), TARGET, OPTIONAL, INTENT(inout) :: &
            v1, v2, v3, s1, s2, s3
        LOGICAL, OPTIONAL, INTENT(in) :: geom, corners, normal
        INTEGER(intk), OPTIONAL, INTENT(in) :: forward
        CHARACTER(len=1), OPTIONAL, INTENT(in) :: ityp
        INTEGER(intk), INTENT(in), OPTIONAL :: minlvl, maxlvl

        ! Local variables
        LOGICAL :: has_v1, has_v2, has_v3, has_s1, has_s2, has_s3

        has_v1 = .FALSE.
        has_v2 = .FALSE.
        has_v3 = .FALSE.
        has_s1 = .FALSE.
        has_s2 = .FALSE.
        has_s3 = .FALSE.

        NULLIFY(u)
        NULLIFY(v)
        NULLIFY(w)
        NULLIFY(p1)
        NULLIFY(p2)
        NULLIFY(p3)

        IF (PRESENT(v1)) THEN
            has_v1 = .TRUE.
            u => v1
        END IF

        IF (PRESENT(v2)) THEN
            has_v2 = .TRUE.
            v => v2
        END IF

        IF (PRESENT(v3)) THEN
            has_v3 = .TRUE.
            w => v3
        END IF

        IF (PRESENT(s1)) THEN
            has_s1 = .TRUE.
            p1 => s1
        END IF

        IF (PRESENT(s2)) THEN
            has_s2 = .TRUE.
            p2 => s2
        END IF

        IF (PRESENT(s3)) THEN
            has_s3 = .TRUE.
            p3 => s3
        END IF

        connect_integer = .TRUE.
        CALL connect_impl(ilevel, layers, has_v1, has_v2, has_v3, &
            has_s1, has_s2, has_s3, geom, corners, normal, forward, ityp, &
            minlvl, maxlvl)
    END SUBROUTINE connect_int


    ! Main connect function
    SUBROUTINE connect_impl(ilevel, layers, has_v1, has_v2, has_v3, &
            has_s1, has_s2, has_s3, geom, corners, normal, forward, ityp, &
            minlvl, maxlvl)

        ! The ilevel and nplane parameters are the only required
        ! parameters
        INTEGER(intk), INTENT(in), OPTIONAL :: ilevel, layers

        ! I am surprised this is allowed as these are not optional...
        LOGICAL, INTENT(in) :: has_v1, has_v2, has_v3, has_s1, has_s2, has_s3

        ! Optional parameters to control special behaviour
        LOGICAL, OPTIONAL, INTENT(in) :: geom, corners, normal
        INTEGER(intk), OPTIONAL, INTENT(in) :: forward
        CHARACTER(len=1), OPTIONAL, INTENT(in) :: ityp
        INTEGER(intk), INTENT(in), OPTIONAL :: minlvl, maxlvl

        CALL start_timer(150)

        ! Check if the connection information has been created
        IF (isinit .EQV. .FALSE.) THEN
            WRITE(*, *) "'connect' not initialized"
            CALL errr(__FILE__, __LINE__)
        END IF

        ! Check that no other transfers are in progress
        IF (nSend > 0 .OR. nRecv > 0) THEN
            WRITE(*, *) "Other transfer in progress."
            CALL errr(__FILE__, __LINE__)
        END IF

        ! Check that number of connect layers are either 1 or 2
        IF (PRESENT(layers)) THEN
            IF (.NOT. (layers  ==  1 .OR. layers  ==  2)) THEN
                WRITE(*, *) "Invalid layers=", layers
                CALL errr(__FILE__, __LINE__)
            END IF
            nplane = layers
        ELSE
            nplane = 1
        END IF

        ! If one vector argument is given, check that all three is present.
        nVars = 0
        IF (has_v1) THEN
            IF (.NOT. (has_v2 .AND. has_v3)) THEN
                WRITE(*, *) "If one vector arg is present, all three " &
                    // "must be present."
                CALL errr(__FILE__, __LINE__)
            END IF
            nVars = nVars + 3
        ELSE IF (has_v2 .OR. has_v3) THEN
            WRITE(*, *) "If one vector arg is present, all three " &
                // "must be present."
            CALL errr(__FILE__, __LINE__)
        END IF

        ! Set pointers to scalars
        IF (has_s1) THEN
            nVars = nVars + 1
        END IF
        IF (has_s2) THEN
            nVars = nVars + 1
        END IF
        IF (has_s3) THEN
            nVars = nVars + 1
        END IF

        geometry = .FALSE.
        IF (PRESENT(geom)) THEN
            IF (geom .EQV. .TRUE.) THEN
                geometry = .TRUE.
            END IF
        END IF

        vertices = .FALSE.
        IF (PRESENT(corners)) THEN
            IF (corners .EQV. .TRUE.) THEN
                vertices = .TRUE.
            END IF
        END IF

        fwd = 0
        IF (PRESENT(forward)) THEN
            IF (vertices) THEN
                CALL errr(__FILE__, __LINE__)
            END IF
            IF (forward > 0) THEN
                fwd = 1
            ELSE IF (forward < 0) THEN
                fwd = -1
            END IF
        END IF

        ! If exchanging surface normal component, check that a vector is
        ! given.
        sn = .FALSE.
        IF (PRESENT(normal)) THEN
            IF (vertices) THEN
                CALL errr(__FILE__, __LINE__)
            END IF
            IF (normal .eqv. .true.) THEN
                IF (has_v1) THEN
                    nVars = nVars - 2
                    sn = .TRUE.
                ELSE
                    WRITE(*, *) "normal=.TRUE. require a vector v1, v2, v3."
                    CALL errr(__FILE__, __LINE__)
                END IF
            END IF
        END IF

        ! Check special flags
        flag = ' '
        IF (PRESENT(ityp)) THEN
            IF ((ityp == 'W') .OR. (ityp == 'Y')) THEN
                flag = ityp
            ELSE
                write(*, *) "Invalid ITYP=", ityp
                CALL errr(__FILE__, __LINE__)
            END IF
        END IF
        IF (flag == 'W') THEN
            ! Check that required fields are present
            IF (has_v1 .AND. has_v2 .AND. has_v3 .AND. has_s1) THEN
                CONTINUE
            ELSE
                CALL errr(__FILE__, __LINE__)
            END IF

            ! Adjust nvars
            nvars = nvars - 3
        END IF

        ! Not specifying any fields would be very strange
        IF (nVars == 0) THEN
            WRITE(*, *) "You have not specified any fields to exchange."
            CALL errr(__FILE__, __LINE__)
        END IF
        ! TODO: Check why? Maybe because of buffer limitations?
        ! In the meantime, comment out
        ! IF (nVars > maxnVars) THEN
        !     WRITE(*,*) "You specified too many fields to exchange."
        !     CALL errr(__FILE__, __LINE__)
        ! END IF

        ! Either ilevel or both minlvl and maxlvl
        IF (PRESENT(ilevel) .AND. (PRESENT(minlvl) .OR. PRESENT(maxlvl))) THEN
            CALL errr(__FILE__, __LINE__)
        END IF
        IF (PRESENT(minlvl) .NEQV. PRESENT(maxlvl)) THEN
            CALL errr(__FILE__, __LINE__)
        END IF

        minconlvl = minlevel
        maxconlvl = maxlevel
        IF (PRESENT(ilevel)) THEN
            minconlvl = ilevel
            maxconlvl = ilevel
        ELSE IF (PRESENT(minlvl) .AND. PRESENT(maxlvl)) THEN
            minconlvl = minlvl
            maxconlvl = maxlvl
        END IF
        IF ((minconlvl < minlevel) .OR. (maxconlvl > maxlevel)) THEN
            WRITE(*, *) "Invalid levels: ", minconlvl, maxconlvl
            CALL errr(__FILE__, __LINE__)
        END IF

        ! Exchnage data
        CALL recv_all()
        CALL send_all()
        CALL process_bufs()

        ! Clear counters and unset pointers. This is important to avoid
        ! problems at next function call
        nRecv = 0
        nSend = 0
        nplane = 0
        nVars = 0

        sn = .FALSE.
        flag = ' '
        fwd = 0
        vertices = .FALSE.
        geometry = .FALSE.

        NULLIFY(u)
        NULLIFY(v)
        NULLIFY(w)
        NULLIFY(p1)
        NULLIFY(p2)
        NULLIFY(p3)

        CALL stop_timer(150)
    END SUBROUTINE connect_impl


    FUNCTION decide(i, list) RESULT(exchange)
        INTEGER(intk), INTENT(in) :: i
        INTEGER(intk), INTENT(in) :: list(:, :)

        INTEGER(intk) :: ifacerecv, ilevel
        LOGICAL :: exchange

        exchange = .TRUE.
        ifacerecv = list(5, i)

        IF (i > SIZE(list, 2)) THEN
            exchange = .FALSE.
            RETURN
        END IF

        ! Levels
        ilevel = level(list(3, i))
        IF (ilevel > maxConLvl .OR. ilevel < minConLvl) THEN
            exchange = .FALSE.
        END IF

        ! Only exchange geometry
#ifndef _IB_CUTCELL_MOVING_
        IF (geometry .AND. list(8, i) == 0) THEN
            exchange = .FALSE.
        END IF
#endif

        ! Corners
        IF ((.NOT. vertices) .AND. (ifacerecv > 6)) THEN
            exchange = .FALSE.
        END IF

        ! Forward
        IF (fwd == 1 .AND. &
            (ifacerecv == 2 .OR. ifacerecv == 4 .OR. ifacerecv == 6)) THEN
            exchange = .FALSE.
        ELSE IF (fwd == -1 .AND. &
            (ifacerecv == 1 .OR. ifacerecv == 3 .OR. ifacerecv == 5)) THEN
            exchange = .FALSE.
        END IF

        RETURN
    END FUNCTION decide


    ! Perform all Recv-calls
    SUBROUTINE recv_all()

        INTEGER(intk) :: i, iprocnbr, igrid, iface, faceArea
        LOGICAL :: exchange

        ! Post all receive calls
        recvCounter = 0
        messageLength = 0
        nRecv = 0
        recvIdxList = 0

        DO i = 1, iRecv
            exchange = decide(i, recvConns)
            iprocnbr = recvConns(2, i)

            ! Communication with self is handled specially in
            ! send_all - nothing to do here
            IF (iprocnbr == myid) THEN
                CYCLE
            END IF

            IF (exchange) THEN
                igrid         = recvConns(3, i)
                iface         = recvConns(5, i)

                faceArea      = face_area(igrid, iface, nplane, flag)
                nRecvFaces = nRecvFaces + 1
                recvIdxList(1, i) = iprocnbr
                recvIdxList(2, i) = nVars*faceArea
                recvIdxList(3, i) = recvCounter + messageLength
                messageLength = messageLength + nVars*faceArea
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
        ! Subroutine arguments
        INTEGER(int32), INTENT(in) :: iprocnbr

        nRecv = nRecv + 1
        recvList(nRecv) = iprocnbr

        IF (connect_integer) THEN
            IF (recvcounter + messagelength > SIZE(irecvbuf)) THEN
                CALL errr(__FILE__, __LINE__)
            END IF

            CALL MPI_Irecv(irecvbuf(recvcounter+1), messagelength, &
                mglet_mpi_ifk, iprocnbr, 1, MPI_COMM_WORLD, recvreqs(nrecv))
        ELSE
            IF (recvcounter + messagelength > SIZE(recvbuf)) THEN
                CALL errr(__FILE__, __LINE__)
            END IF

            CALL MPI_Irecv(recvbuf(recvcounter+1), messagelength, &
                mglet_mpi_real, iprocnbr, 1, MPI_COMM_WORLD, recvreqs(nrecv))
        END IF

        recvcounter = recvcounter + messagelength
        messagelength = 0
    END SUBROUTINE post_recv


    ! Perform all send calls
    SUBROUTINE send_all()
        INTEGER(intk) :: i, iprocnbr, igrid, ifacerecv, faceArea
        LOGICAL :: exchange

        ! Pack all buffers and send data
        sendCounter = 0
        messageLength = 0
        nSend = 0

        DO i = 1, iSend
            exchange = decide(i, sendConns)
            iprocnbr = sendConns(1, i)

            ! Communication with self copies directly from source to
            ! destination grid - then skip the rest
            IF (iprocnbr == myid .AND. exchange) THEN
                CALL connect_self(i)
                CYCLE
            END IF

            IF (exchange) THEN
                igrid         = sendConns(3, i)
                ifacerecv     = sendConns(5, i)
                faceArea      = face_area(igrid, ifacerecv, nplane, flag)

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

        ! Local variables
        ! none...

        nSend = nSend + 1
        sendList(nSend) = iprocnbr

        IF (connect_integer) THEN
            CALL MPI_Isend(isendbuf(sendCounter + 1), messageLength, &
                mglet_mpi_ifk, iprocnbr, 1, MPI_COMM_WORLD, sendReqs(nSend))
        ELSE
            CALL MPI_Isend(sendbuf(sendCounter + 1), messageLength, &
                mglet_mpi_real, iprocnbr, 1, MPI_COMM_WORLD, sendReqs(nSend))
        END IF

        sendCounter = sendCounter + messageLength
        messageLength = 0
    END SUBROUTINE post_send


    ! Write Send buffers
    !
    ! Write the relevant fields into the send buffers
    SUBROUTINE write_buffer(sendId)

        ! Input parameter
        INTEGER(int32), INTENT(in) :: sendId

        ! Indices of start- and stop of iteration over boundary face
        INTEGER(intk) :: istart, istop, jstart, jstop, kstart, kstop

        ! Grid to send from
        ! Must be intk because it intreface with MGLET
        INTEGER(intk) :: igrid, ifacerecv, ifacesend

        ! Message sizes
        ! Must be int32 because it iterface with MPI
        INTEGER(int32) :: thisMessageLength, faceArea
        INTEGER(int32) :: icount, offset

        ! Flags to indicate exchange of U, V, W
        LOGICAL :: exU, exV, exW, exp1

        ! Set variables from send table
        igrid         = sendConns(4, sendId)
        ifacerecv     = sendConns(5, sendId)
        ifacesend     = sendConns(6, sendId)

        ! Get start- and stop indices of grid
        CALL start_and_stop(igrid, facenbr(ifacerecv), istart, istop, &
            jstart, jstop, kstart, kstop, nplane, flag, nghost=1)
        CALL corr_start_stop(igrid, ifacesend, ifacerecv, &
            istart, istop, jstart, jstop, kstart, kstop, nplane, flag)

        faceArea = (istop-istart+1)*(jstop-jstart+1)*(kstop-kstart+1)
        thisMessageLength = nVars*faceArea

        ! Check that buffer does not overflow
        ! TODO: not correct for integer connects
        IF (sendcounter + messagelength + thismessagelength &
                > SIZE(sendbuf)) THEN
            CALL errr(__FILE__, __LINE__)
        END IF

        ! Reset message size counter
        offset = sendCounter + messageLength
        icount = offset

        ! Fill buffers
        IF (flag == 'W') THEN
            exU = (ifacerecv == 1)
        ELSE
            exU = (sn .AND. ifacerecv < 3) .OR. (.NOT. sn)
        END IF
        IF (ASSOCIATED(u) .AND. exU) THEN
            CALL write_single_buffer(u, icount, igrid, istart, istop, &
                jstart, jstop, kstart, kstop)
        END IF

        IF (flag == 'W') THEN
            exV = (ifacerecv == 3)
        ELSE
            exV = (sn .AND. (ifacerecv > 2 .AND. ifacerecv < 5)) .OR. &
                (.NOT. sn)
        END IF
        IF (ASSOCIATED(v) .AND. exV) THEN
            CALL write_single_buffer(v, icount, igrid, istart, istop, &
                jstart, jstop, kstart, kstop)
        END IF

        IF (flag == 'W') THEN
            exW = (ifacerecv == 5)
        ELSE
            exW = (sn .AND. ifacerecv > 4) .OR. (.NOT. sn)
        END IF
        IF (ASSOCIATED(w) .AND. exW) THEN
            CALL write_single_buffer(w, icount, igrid, istart, istop, &
                jstart, jstop, kstart, kstop)
        END IF

        IF (flag == 'W') THEN
            exp1 = (ifacerecv == 2) .OR. (ifacerecv == 4) .OR. &
                (ifacerecv == 6)
        ELSE
            exp1 = .TRUE.
        END IF
        IF (ASSOCIATED(p1) .AND. exp1) THEN
            CALL write_single_buffer(p1, icount, igrid, istart, istop, &
                jstart, jstop, kstart, kstop)
        END IF

        IF (ASSOCIATED(p2)) THEN
            CALL write_single_buffer(p2, icount, igrid, istart, istop, &
                jstart, jstop, kstart, kstop)
        END IF

        IF (ASSOCIATED(p3)) THEN
            CALL write_single_buffer(p3, icount, igrid, istart, istop, &
                jstart, jstop, kstart, kstop)
        END IF

        ! Check that message length was calculated correctly
        IF (thismessagelength /= (icount - offset)) THEN
            write(*, *) "thismessagelength:", thismessagelength, &
                "icount:", icount
            CALL errr(__FILE__, __LINE__)
        END IF
    END SUBROUTINE write_buffer


    SUBROUTINE write_single_buffer(field, icount, igrid, istart, istop, &
            jstart, jstop, kstart, kstop)

        ! Subroutine arguments
        CLASS(basefield_t), INTENT(in) :: field
        INTEGER(intk), INTENT(inout) :: icount
        INTEGER(intk), INTENT(in) :: igrid
        INTEGER(intk), INTENT(in) :: istart, istop, jstart, jstop, &
            kstart, kstop

        ! Local variables
        INTEGER(intk) :: k, j, i
        REAL(realk), POINTER, CONTIGUOUS :: rarr(:, :, :)
        INTEGER(ifk), POINTER, CONTIGUOUS :: iarr(:, :, :)

        SELECT TYPE(field)
        TYPE IS (field_t)
            CALL field%get_ptr(rarr, igrid)
            DO i = istart, istop
                DO j = jstart, jstop
                    DO k = kstart, kstop
                        icount = icount + 1
                        sendbuf(icount) = rarr(k, j, i)
                    END DO
                END DO
            END DO
        TYPE IS (intfield_t)
            CALL field%get_ptr(iarr, igrid)
            DO i = istart, istop
                DO j = jstart, jstop
                    DO k = kstart, kstop
                        icount = icount + 1
                        isendbuf(icount) = iarr(k, j, i)
                    END DO
                END DO
            END DO
        CLASS DEFAULT
            CALL errr(__FILE__, __LINE__)
        END SELECT
    END SUBROUTINE write_single_buffer


    ! Read Receive buffers
    !
    ! Write the contents of the receive buffers back in their
    ! matching fields
    SUBROUTINE read_buffer(recvId)
        ! Input parameter
        INTEGER(intk), INTENT(in) :: recvId

        ! Indices of start- and stop of iteration over boundary face
        INTEGER(intk) :: istart, istop, jstart, jstop, kstart, kstop

        ! Grid to send from
        ! Must be intk because it intreface with MGLET
        INTEGER(intk) :: igrid, ifacerecv

        ! Message sizes
        ! Must be int32 because it iterface with MPI
        INTEGER(int32) :: offset, icount

        ! Flags to indicate exchange of U, V, W
        LOGICAL :: exU, exV, exW, exp1

        ! Set variables from send table
        igrid         = recvConns(3, recvId)
        ifacerecv     = recvConns(5, recvId)

        ! Get start- and stop indices of grid
        CALL start_and_stop(igrid, ifacerecv, &
            istart, istop, jstart, jstop, kstart, kstop, nplane, flag)

        ! Zeroise message size counter
        offset = recvidxlist(3, recvid)
        icount = offset

        IF (flag == 'W') THEN
            exU = (ifacerecv == 1)
        ELSE
            exU = (sn .AND. ifacerecv < 3) .OR. (.NOT. sn)
        END IF
        IF (ASSOCIATED(u) .AND. exU) THEN
            CALL read_single_buffer(u, icount, igrid, istart, istop, &
                jstart, jstop, kstart, kstop)
        END IF

        IF (flag == 'W') THEN
            exV = (ifacerecv == 3)
        ELSE
            exV = (sn .AND. (ifacerecv > 2 .AND. ifacerecv < 5)) .OR. &
                (.NOT. sn)
        END IF
        IF (ASSOCIATED(v) .AND. exV) THEN
            CALL read_single_buffer(v, icount, igrid, istart, istop, &
                jstart, jstop, kstart, kstop)
        END IF

        IF (flag == 'W') THEN
            exW = (ifacerecv == 5)
        ELSE
            exW = (sn .AND. ifacerecv > 4) .OR. (.NOT. sn)
        END IF
        IF (ASSOCIATED(w) .AND. exW) THEN
            CALL read_single_buffer(w, icount, igrid, istart, istop, &
                jstart, jstop, kstart, kstop)
        END IF

        IF (flag == 'W') THEN
            exp1 = (ifacerecv == 2) .OR. (ifacerecv == 4) .OR. &
                (ifacerecv == 6)
        ELSE
            exp1 = .TRUE.
        END IF
        IF (ASSOCIATED(p1) .AND. exp1) THEN
            CALL read_single_buffer(p1, icount, igrid, istart, istop, &
                jstart, jstop, kstart, kstop)
        END IF

        IF (ASSOCIATED(p2)) THEN
            CALL read_single_buffer(p2, icount, igrid, istart, istop, &
                jstart, jstop, kstart, kstop)
        END IF

        IF (ASSOCIATED(p3)) THEN
            CALL read_single_buffer(p3, icount, igrid, istart, istop, &
                jstart, jstop, kstart, kstop)
        END IF

        ! Check that message length is calculated correctly
        IF ((icount - offset) /= recvidxlist(2, recvid)) THEN
            WRITE(*, *) "icount:", icount, &
                "recvidxlist(2, recvid):", recvidxlist(2, recvid)
            CALL errr(__FILE__, __LINE__)
        END IF
    END SUBROUTINE read_buffer


    SUBROUTINE read_single_buffer(field, icount, igrid, istart, istop, &
            jstart, jstop, kstart, kstop)

        ! Subroutine arguments
        CLASS(basefield_t), INTENT(inout) :: field
        INTEGER(intk), INTENT(inout) :: icount
        INTEGER(intk), INTENT(in) :: igrid
        INTEGER(intk), INTENT(in) :: istart, istop, jstart, jstop, &
            kstart, kstop

        ! Local variables
        INTEGER(intk) :: k, j, i
        REAL(realk), POINTER, CONTIGUOUS :: rarr(:, :, :)
        INTEGER(ifk), POINTER, CONTIGUOUS :: iarr(:, :, :)

        SELECT TYPE(field)
        TYPE IS (field_t)
            CALL field%get_ptr(rarr, igrid)
            DO i = istart, istop
                DO j = jstart, jstop
                    DO k = kstart, kstop
                        icount = icount + 1
                        rarr(k, j, i) = recvbuf(icount)
                    END DO
                END DO
            END DO
        TYPE IS (intfield_t)
            CALL field%get_ptr(iarr, igrid)
            DO i = istart, istop
                DO j = jstart, jstop
                    DO k = kstart, kstop
                        icount = icount + 1
                        iarr(k, j, i) = irecvbuf(icount)
                    END DO
                END DO
            END DO
        CLASS DEFAULT
            CALL errr(__FILE__, __LINE__)
        END SELECT
    END SUBROUTINE read_single_buffer


    ! Connect to self
    !
    ! Directly copes data from source to destination buffer
    SUBROUTINE connect_self(sendId)
        ! Input parameter
        INTEGER(int32), INTENT(in) :: sendId

        ! Indices of start- and stop of iteration over boundary face
        ! Source face
        INTEGER(intk) :: istart, istop, jstart, jstop, kstart, kstop

        ! Indices of start- and stop of iteration over boundary face
        ! Destination face
        INTEGER(intk) :: istart_d, istop_d, jstart_d, jstop_d, kstart_d, kstop_d

        ! Grid to send from
        ! Must be intk because it intreface with MGLET
        INTEGER(intk) :: igrid, igrid_d, ifacerecv, ifacesend

        ! Message sizes
        ! Must be int32 because it iterface with MPI
        INTEGER(int32) :: dest_size, source_size

        ! Flags to indicate exchange of U, V, W
        LOGICAL :: exU, exV, exW, exp1

        ! Set variables from send table
        igrid_d       = sendConns(3, sendId)
        igrid         = sendConns(4, sendId)
        ifacerecv     = sendConns(5, sendId)
        ifacesend     = sendConns(6, sendId)

        ! Get start- and stop indices of source grid
        CALL start_and_stop(igrid, facenbr(ifacerecv), &
            istart, istop, jstart, jstop, kstart, kstop, nplane, flag, nghost=1)
        CALL corr_start_stop(igrid, ifacesend, ifacerecv, &
            istart, istop, jstart, jstop, kstart, kstop, nplane, flag)

        ! Get start- and stop indices of destination grid
        CALL start_and_stop(igrid_d, ifacerecv, &
            istart_d, istop_d, jstart_d, jstop_d, kstart_d, kstop_d, &
            nplane, flag)

        ! Sanity check of message length
        source_size = (istop-istart+1)*(jstop-jstart+1)*(kstop-kstart+1)
        dest_size = (istop_d-istart_d+1) &
            *(jstop_d-jstart_d+1)*(kstop_d-kstart_d+1)
        IF (source_size /= dest_size) THEN
            CALL errr(__FILE__, __LINE__)
        END IF

        IF (flag == 'W') THEN
            exU = (ifacerecv == 1)
        ELSE
            exU = (sn .AND. ifacerecv < 3) .OR. (.NOT. sn)
        END IF
        IF (ASSOCIATED(u) .AND. exU) THEN
            CALL connect_self_single(u, igrid, igrid_d, istart, istop, &
                jstart, jstop, kstart, kstop, istart_d, istop_d, &
                jstart_d, jstop_d, kstart_d, kstop_d)
        END IF

        IF (flag == 'W') THEN
            exV = (ifacerecv == 3)
        ELSE
            exV = (sn .AND. (ifacerecv > 2 .AND. ifacerecv < 5)) .OR. &
                (.NOT. sn)
        END IF
        IF (ASSOCIATED(v) .AND. exV) THEN
            CALL connect_self_single(v, igrid, igrid_d, istart, istop, &
                jstart, jstop, kstart, kstop, istart_d, istop_d, &
                jstart_d, jstop_d, kstart_d, kstop_d)
        END IF

        IF (flag == 'W') THEN
            exW = (ifacerecv == 5)
        ELSE
            exW = (sn .AND. ifacerecv > 4) .OR. (.NOT. sn)
        END IF
        IF (ASSOCIATED(w) .AND. exW) THEN
            CALL connect_self_single(w, igrid, igrid_d, istart, istop, &
                jstart, jstop, kstart, kstop, istart_d, istop_d, &
                jstart_d, jstop_d, kstart_d, kstop_d)
        END IF

        IF (flag == 'W') THEN
            exp1 = (ifacerecv == 2) .OR. (ifacerecv == 4) .OR. &
                (ifacerecv == 6)
        ELSE
            exp1 = .TRUE.
        END IF
        IF (ASSOCIATED(p1) .AND. exp1) THEN
            CALL connect_self_single(p1, igrid, igrid_d, istart, istop, &
                jstart, jstop, kstart, kstop, istart_d, istop_d, &
                jstart_d, jstop_d, kstart_d, kstop_d)
        END IF

        IF (ASSOCIATED(p2)) THEN
            CALL connect_self_single(p2, igrid, igrid_d, istart, istop, &
                jstart, jstop, kstart, kstop, istart_d, istop_d, &
                jstart_d, jstop_d, kstart_d, kstop_d)
        END IF

        IF (ASSOCIATED(p3)) THEN
            CALL connect_self_single(p3, igrid, igrid_d, istart, istop, &
                jstart, jstop, kstart, kstop, istart_d, istop_d, &
                jstart_d, jstop_d, kstart_d, kstop_d)
        END IF
    END SUBROUTINE connect_self


    SUBROUTINE connect_self_single(field, igrid, igrid_d, istart, istop, &
            jstart, jstop, kstart, kstop, istart_d, istop_d, jstart_d, &
            jstop_d, kstart_d, kstop_d)

        ! Subroutine arguments
        CLASS(basefield_t), INTENT(inout) :: field
        INTEGER(intk), INTENT(in) :: igrid, igrid_d
        INTEGER(intk), INTENT(in) :: istart, istop, jstart, jstop, &
            kstart, kstop
        INTEGER(intk), INTENT(in) :: istart_d, istop_d, jstart_d, jstop_d, &
            kstart_d, kstop_d

        ! Local variables
        REAL(realk), POINTER, CONTIGUOUS :: src_rarr(:, :, :), dst_rarr(:, :, :)
        INTEGER(ifk), POINTER, CONTIGUOUS :: src_iarr(:, :, :), &
            dst_iarr(:, :, :)
        INTEGER(intk) :: i, j, k, ioff, joff, koff

        koff = kstart - kstart_d
        joff = jstart - jstart_d
        ioff = istart - istart_d

        SELECT TYPE(field)
        TYPE IS (field_t)
            CALL field%get_ptr(src_rarr, igrid)
            CALL field%get_ptr(dst_rarr, igrid_d)
            DO i = istart_d, istop_d
                DO j = jstart_d, jstop_d
                    DO k = kstart_d, kstop_d
                        dst_rarr(k, j, i) = &
                            src_rarr(k + koff, j + joff, i + ioff)
                    END DO
                END DO
            END DO
        TYPE IS (intfield_t)
            CALL field%get_ptr(src_iarr, igrid)
            CALL field%get_ptr(dst_iarr, igrid_d)
            DO i = istart_d, istop_d
                DO j = jstart_d, jstop_d
                    DO k = kstart_d, kstop_d
                        dst_iarr(k, j, i) = &
                            src_iarr(k + koff, j + joff, i + ioff)
                    END DO
                END DO
            END DO
        CLASS DEFAULT
            CALL errr(__FILE__, __LINE__)
        END SELECT
    END SUBROUTINE connect_self_single


    ! Process receive buffers as they arrive, wait for send
    ! buffers to be free
    SUBROUTINE process_bufs()

        INTEGER(int32) :: idx, i
        TYPE(MPI_Status) :: recvstatus
        INTEGER(int32) :: recvmessagelen
        INTEGER(int32) :: unpacklen

        DO WHILE (.TRUE.)
            IF (nRecv == 0) EXIT
            CALL MPI_Waitany(nRecv, recvReqs, idx, recvstatus)

            IF (idx /= MPI_UNDEFINED) THEN
                IF (connect_integer) THEN
                    CALL MPI_Get_count(recvstatus, mglet_mpi_ifk, &
                        recvmessagelen)
                ELSE
                    CALL MPI_Get_count(recvstatus, mglet_mpi_real, &
                        recvmessagelen)
                END IF

                unpackLen = 0
                DO i = 1, iRecv
                    IF (recvIdxList(1, i) == recvList(idx) &
                            .AND. recvIdxList(2, i) > 0) THEN
                        CALL read_buffer(i)
                        unpacklen = unpacklen + recvIdxList(2, i)
                    END IF
                END DO

                IF (recvmessagelen /= unpacklen) THEN
                    CALL errr(__FILE__, __LINE__)
                END IF
            ELSE
                EXIT
            END IF
        END DO
        CALL MPI_Waitall(nSend, sendReqs, MPI_STATUSES_IGNORE)
    END SUBROUTINE process_bufs


    SUBROUTINE init_connect2()
        CALL init_connect_core()

        ! The maximum number of concurrent communications are the number
        ! of processes
        ALLOCATE(recvidxlist(3, maxconns))
        ALLOCATE(sendlist(numprocs))
        ALLOCATE(recvlist(numprocs))
        ALLOCATE(sendreqs(numprocs))
        ALLOCATE(recvreqs(numprocs))
        recvidxlist = 0
        sendlist = 0
        recvlist = 0
        nrecv = 0
        nsend = 0

        NULLIFY(u)
        NULLIFY(v)
        NULLIFY(w)
        NULLIFY(p1)
        NULLIFY(p2)
        NULLIFY(p3)

        isinit = .TRUE.
    END SUBROUTINE init_connect2


    SUBROUTINE finish_connect2()
        isinit = .FALSE.

        DEALLOCATE(recvIdxList)
        DEALLOCATE(sendList)
        DEALLOCATE(recvList)
        DEALLOCATE(sendReqs)
        DEALLOCATE(recvReqs)

        CALL finish_connect_core()
    END SUBROUTINE finish_connect2
END MODULE connect2_mod
