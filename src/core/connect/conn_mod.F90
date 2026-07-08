MODULE conn_mod
    USE MPI_f08
    USE precision_mod
    USE commbuf_mod, ONLY: sendbuf, recvbuf
    USE err_mod, ONLY: errr
    USE timer_mod, ONLY: start_timer, stop_timer
    USE grids_mod, ONLY: mygrids, nmygrids, level, idprocofgrd, itypboconds, &
        maxlevel, minlevel, get_neighbours, get_mgdims
    USE comms_mod, ONLY: myid, numprocs
    USE field_mod
    USE connect_core_mod

    IMPLICIT NONE (type, external)
    PRIVATE

    ! Lists that hold the send and receive request arrays
    TYPE(MPI_Request), ALLOCATABLE :: sendreqs(:), recvreqs(:)
    INTEGER(int32), ALLOCATABLE :: sendlist(:), recvlist(:)
    INTEGER(intk), ALLOCATABLE :: recvidxlist(:, :)
    INTEGER(intk) :: nsend, nrecv

    PUBLIC :: conn, init_conn, finish_conn

CONTAINS

    SUBROUTINE conn(ilevel, layers, v1, v2, v3, s1, s2, s3, corners, normal, &
            forward, ityp)
        ! conn is a more compact version of connect aiming for CPU offloading.
        ! It will only connect real fields, not integer fields, and does
        ! not have the same capabilities as connect. Over time features
        ! can be transferred, however, initially the goal is to have a
        ! light and easy code for GPU offloading while keeping the
        ! traditional connect in place.

        ! Subroutine arguments
        INTEGER(intk), INTENT(in), OPTIONAL :: ilevel, layers
        TYPE(field_t), OPTIONAL, INTENT(inout) :: &
            v1, v2, v3, s1, s2, s3
        LOGICAL, OPTIONAL, INTENT(in) :: corners, normal
        INTEGER(intk), OPTIONAL, INTENT(in) :: forward
        CHARACTER(len=1), OPTIONAL, INTENT(in) :: ityp

        ! Local variables
        INTEGER(intk) :: minconlvl, maxconlvl, nplane, fwd, nvars
        LOGICAL :: vertices, normal2
        CHARACTER(len=1) :: flag

        IF (PRESENT(ilevel)) THEN
            minconlvl = ilevel
            maxconlvl = ilevel
        ELSE
            minconlvl = minlevel
            maxconlvl = maxlevel
        END IF

        nplane = 1
        IF (PRESENT(layers)) THEN
            nplane = layers
        END IF
        IF (nplane < 1 .OR. nplane > 2) THEN
            CALL errr(__FILE__, __LINE__)
        END IF

        vertices = .FALSE.
        IF (PRESENT(corners)) THEN
            vertices = corners
        END IF

        fwd = 0
        IF (PRESENT(forward)) THEN
            fwd = forward
        END IF

        IF (fwd /= 0 .AND. vertices) THEN
            CALL errr(__FILE__, __LINE__)
        END IF

        flag = ' '
        IF (PRESENT(ityp)) THEN
            IF (ityp == "W") THEN
                flag = 'W'
            ELSE
                CALL errr(__FILE__, __LINE__)
            END IF
        END IF

        nvars = 0
        IF (PRESENT(v1)) THEN
            nvars = nvars + 1
        END IF
        IF (PRESENT(v2)) THEN
            nvars = nvars + 1
        END IF
        IF (PRESENT(v3)) THEN
            nvars = nvars + 1
        END IF
        IF (PRESENT(s1)) THEN
            nvars = nvars + 1
        END IF
        IF (PRESENT(s2)) THEN
            nvars = nvars + 1
        END IF
        IF (PRESENT(s3)) THEN
            nvars = nvars + 1
        END IF
        IF (nvars == 0) THEN
            CALL errr(__FILE__, __LINE__)
        END IF

        normal2 = .FALSE.
        IF (PRESENT(normal)) THEN
            normal2 = normal
        END IF
        IF (normal2) THEN
            IF (.NOT. PRESENT(v1) .OR. .NOT. PRESENT(v2) .OR. &
                    .NOT. PRESENT(v3)) THEN
                CALL errr(__FILE__, __LINE__)
            END IF
            IF (vertices) THEN
                CALL errr(__FILE__, __LINE__)
            END IF
            nvars = nvars - 2
        END IF

        CALL start_timer(150)

        CALL recv_all(minconlvl, maxconlvl, nplane, vertices, normal2, fwd, &
            flag, nvars)
        CALL send_all(minconlvl, maxconlvl, nplane, vertices, normal2, fwd, &
            flag, nvars, v1, v2, v3, s1, s2, s3)
        CALL process_bufs(nplane, normal2, flag, v1, v2, v3, s1, s2, s3)

        CALL stop_timer(150)
    END SUBROUTINE conn


    ! Perform all Recv-calls
    SUBROUTINE recv_all(minconlvl, maxconlvl, nplane, vertices, normal, fwd, &
            flag, nvars)

        ! Subroutine arguments
        INTEGER(intk), INTENT(in) :: minconlvl, maxconlvl, nplane
        LOGICAL, INTENT(in) :: vertices, normal
        INTEGER(intk), INTENT(in) :: fwd
        CHARACTER(len=1), INTENT(in) :: flag
        INTEGER(intk), INTENT(in) :: nvars

        ! Local variables
        INTEGER(intk) :: i, iprocnbr, igrid, iface, facearea
        LOGICAL :: exchange, geometry
        INTEGER(int32) :: recvcounter, messagelength

        geometry = .FALSE.

        ! Post all receive calls
        recvcounter = 0
        messagelength = 0
        nrecv = 0
        recvidxlist = 0

        DO i = 1, SIZE(recvconns, 2)
            exchange = decide(i, recvconns, geometry, vertices, fwd, &
                minconlvl, maxconlvl)
            iprocnbr = recvconns(2, i)

            ! Communication with self is handled specially in
            ! send_all - nothing to do here
            IF (iprocnbr == myid) THEN
                CYCLE
            END IF

            IF (exchange) THEN
                igrid = recvconns(3, i)
                iface = recvconns(5, i)

                facearea = face_area(igrid, iface, nplane, flag)
                recvidxlist(1, i) = iprocnbr
                recvidxlist(2, i) = nvars*facearea
                recvidxlist(3, i) = recvcounter + messagelength
                messagelength = messagelength + nvars*facearea
            END IF

            IF (messagelength > 0) THEN
                IF (i == SIZE(recvconns, 2)) THEN
                    CALL post_recv(iprocnbr, messagelength, recvcounter)
                ELSE IF (recvconns(2, i + 1) /= iprocnbr) THEN
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

        nrecv = nrecv + 1
        recvlist(nrecv) = iprocnbr

        IF (recvcounter + messagelength > SIZE(recvbuf)) THEN
            CALL errr(__FILE__, __LINE__)
        END IF

        CALL MPI_Irecv(recvbuf(recvcounter+1), messagelength, &
            mglet_mpi_real, iprocnbr, 1, MPI_COMM_WORLD, recvreqs(nrecv))

        recvcounter = recvcounter + messagelength
        messagelength = 0
    END SUBROUTINE post_recv


    ! Perform all send calls
    SUBROUTINE send_all(minconlvl, maxconlvl, nplane, vertices, normal, fwd, &
            flag, nvars, v1, v2, v3, s1, s2, s3)
        ! Subroutine arguments
        INTEGER(intk), INTENT(in) :: minconlvl, maxconlvl, nplane
        LOGICAL, INTENT(in) :: vertices, normal
        INTEGER(intk), INTENT(in) :: fwd
        CHARACTER(len=1), INTENT(in) :: flag
        INTEGER(intk), INTENT(in) :: nvars
        TYPE(field_t), OPTIONAL, INTENT(inout) :: v1, v2, v3, s1, s2, s3

        ! Local variables
        INTEGER(intk) :: i, iprocnbr, igrid, ifacerecv, facearea
        LOGICAL :: exchange, geometry
        INTEGER(int32) :: sendcounter, messagelength

        geometry = .FALSE.

        ! Pack all buffers and send data
        sendcounter = 0
        messagelength = 0
        nsend = 0

        DO i = 1, SIZE(sendconns, 2)
            exchange = decide(i, sendconns, geometry, vertices, fwd, &
                minconlvl, maxconlvl)
            iprocnbr = sendconns(1, i)

            ! Communication with self copies directly from source to
            ! destination grid - then skip the rest
            IF (iprocnbr == myid .AND. exchange) THEN
                CALL connect_self(i, nplane, normal, flag, v1, v2, v3, &
                    s1, s2, s3)
                CYCLE
            END IF

            IF (exchange) THEN
                igrid = sendconns(3, i)
                ifacerecv = sendconns(5, i)
                facearea = face_area(igrid, ifacerecv, nplane, flag)

                CALL write_buffer(i, messagelength, sendcounter, nplane, &
                    normal, flag, nvars, v1, v2, v3, s1, s2, s3)

                messagelength = messagelength + nvars*facearea
            END IF

            IF (messagelength > 0) THEN
                IF (i == SIZE(sendconns, 2)) THEN
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

        ! Local variables
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
    SUBROUTINE write_buffer(sendid, messagelength, sendcounter, nplane, &
            normal, flag, nvars, v1, v2, v3, s1, s2, s3)

        ! Subroutine arguments
        INTEGER(int32), INTENT(in) :: sendid
        INTEGER(int32), INTENT(inout) :: messagelength
        INTEGER(int32), INTENT(inout) :: sendcounter
        INTEGER(int32), INTENT(in) :: nplane
        LOGICAL, INTENT(in) :: normal
        CHARACTER(len=1), INTENT(in) :: flag
        INTEGER(intk), INTENT(in) :: nvars
        TYPE(field_t), OPTIONAL, INTENT(inout) :: v1, v2, v3, s1, s2, s3

        ! Indices of start- and stop of iteration over boundary face
        INTEGER(intk) :: istart, istop, jstart, jstop, kstart, kstop

        ! Grid to send from
        ! Must be intk because it intreface with MGLET
        INTEGER(intk) :: igrid, ifacerecv, ifacesend

        ! Message sizes
        ! Must be int32 because it iterface with MPI
        INTEGER(int32) :: thismessagelength, facearea
        INTEGER(int32) :: icount, offset

        ! Flags to indicate exchange of U, V, W
        LOGICAL :: exU, exV, exW, exp1

        ! Set variables from send table
        igrid = sendconns(4, sendid)
        ifacerecv = sendconns(5, sendid)
        ifacesend = sendconns(6, sendid)

        ! Get start- and stop indices of grid
        CALL start_and_stop(igrid, facenbr(ifacerecv), istart, istop, &
            jstart, jstop, kstart, kstop, nplane, flag, nghost=1)
        CALL corr_start_stop(igrid, ifacesend, ifacerecv, &
            istart, istop, jstart, jstop, kstart, kstop, nplane, flag)

        facearea = (istop-istart+1)*(jstop-jstart+1)*(kstop-kstart+1)
        thismessagelength = nvars*facearea

        ! Check that buffer does not overflow
        IF (sendcounter + messagelength + thismessagelength &
                > SIZE(sendbuf)) THEN
            CALL errr(__FILE__, __LINE__)
        END IF

        ! Reset message size counter
        offset = sendcounter + messagelength
        icount = offset

        ! Fill buffers
        IF (flag == 'W') THEN
            exU = (ifacerecv == 1)
        ELSE
            exU = (normal .AND. ifacerecv < 3) .OR. (.NOT. normal)
        END IF
        IF (PRESENT(v1) .AND. exU) THEN
            CALL write_single_buffer(v1, icount, igrid, istart, istop, &
                jstart, jstop, kstart, kstop)
        END IF

        IF (flag == 'W') THEN
            exV = (ifacerecv == 3)
        ELSE
            exV = (normal .AND. (ifacerecv > 2 .AND. ifacerecv < 5)) .OR. &
                (.NOT. normal)
        END IF
        IF (PRESENT(v2) .AND. exV) THEN
            CALL write_single_buffer(v2, icount, igrid, istart, istop, &
                jstart, jstop, kstart, kstop)
        END IF

        IF (flag == 'W') THEN
            exW = (ifacerecv == 5)
        ELSE
            exW = (normal .AND. ifacerecv > 4) .OR. (.NOT. normal)
        END IF
        IF (PRESENT(v3) .AND. exW) THEN
            CALL write_single_buffer(v3, icount, igrid, istart, istop, &
                jstart, jstop, kstart, kstop)
        END IF

        IF (flag == 'W') THEN
            exp1 = (ifacerecv == 2) .OR. (ifacerecv == 4) .OR. &
                (ifacerecv == 6)
        ELSE
            exp1 = .TRUE.
        END IF
        IF (PRESENT(s1) .AND. exp1) THEN
            CALL write_single_buffer(s1, icount, igrid, istart, istop, &
                jstart, jstop, kstart, kstop)
        END IF

        IF (PRESENT(s2)) THEN
            CALL write_single_buffer(s2, icount, igrid, istart, istop, &
                jstart, jstop, kstart, kstop)
        END IF

        IF (PRESENT(s3)) THEN
            CALL write_single_buffer(s3, icount, igrid, istart, istop, &
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
        TYPE(field_t), INTENT(in) :: field
        INTEGER(intk), INTENT(inout) :: icount
        INTEGER(intk), INTENT(in) :: igrid
        INTEGER(intk), INTENT(in) :: istart, istop, jstart, jstop, &
            kstart, kstop

        ! Local variables
        INTEGER(intk) :: k, j, i
        REAL(realk), POINTER, CONTIGUOUS :: rarr(:, :, :)

        CALL field%get_ptr(rarr, igrid)
        DO i = istart, istop
            DO j = jstart, jstop
                DO k = kstart, kstop
                    icount = icount + 1
                    sendbuf(icount) = rarr(k, j, i)
                END DO
            END DO
        END DO
    END SUBROUTINE write_single_buffer


    ! Read Receive buffers
    !
    ! Write the contents of the receive buffers back in their
    ! matching fields
    SUBROUTINE read_buffer(recvid, nplane, normal, flag, v1, v2, v3, &
            s1, s2, s3)
        ! Subroutine arguments
        INTEGER(intk), INTENT(in) :: recvid
        INTEGER(int32), INTENT(in) :: nplane
        LOGICAL, INTENT(in) :: normal
        CHARACTER(len=1), INTENT(in) :: flag
        TYPE(field_t), OPTIONAL, INTENT(inout) :: v1, v2, v3, s1, s2, s3

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
        igrid = recvconns(3, recvid)
        ifacerecv = recvconns(5, recvid)

        ! Get start- and stop indices of grid
        CALL start_and_stop(igrid, ifacerecv, &
            istart, istop, jstart, jstop, kstart, kstop, nplane, flag)

        ! Zeroise message size counter
        offset = recvidxlist(3, recvid)
        icount = offset

        IF (flag == 'W') THEN
            exU = (ifacerecv == 1)
        ELSE
            exU = (normal .AND. ifacerecv < 3) .OR. (.NOT. normal)
        END IF
        IF (PRESENT(v1) .AND. exU) THEN
            CALL read_single_buffer(v1, icount, igrid, istart, istop, &
                jstart, jstop, kstart, kstop)
        END IF

        IF (flag == 'W') THEN
            exV = (ifacerecv == 3)
        ELSE
            exV = (normal .AND. (ifacerecv > 2 .AND. ifacerecv < 5)) .OR. &
                (.NOT. normal)
        END IF
        IF (PRESENT(v2) .AND. exV) THEN
            CALL read_single_buffer(v2, icount, igrid, istart, istop, &
                jstart, jstop, kstart, kstop)
        END IF

        IF (flag == 'W') THEN
            exW = (ifacerecv == 5)
        ELSE
            exW = (normal .AND. ifacerecv > 4) .OR. (.NOT. normal)
        END IF
        IF (PRESENT(v3) .AND. exW) THEN
            CALL read_single_buffer(v3, icount, igrid, istart, istop, &
                jstart, jstop, kstart, kstop)
        END IF

        IF (flag == 'W') THEN
            exp1 = (ifacerecv == 2) .OR. (ifacerecv == 4) .OR. &
                (ifacerecv == 6)
        ELSE
            exp1 = .TRUE.
        END IF
        IF (PRESENT(s1) .AND. exp1) THEN
            CALL read_single_buffer(s1, icount, igrid, istart, istop, &
                jstart, jstop, kstart, kstop)
        END IF

        IF (PRESENT(s2)) THEN
            CALL read_single_buffer(s2, icount, igrid, istart, istop, &
                jstart, jstop, kstart, kstop)
        END IF

        IF (PRESENT(s3)) THEN
            CALL read_single_buffer(s3, icount, igrid, istart, istop, &
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
        TYPE(field_t), INTENT(inout) :: field
        INTEGER(intk), INTENT(inout) :: icount
        INTEGER(intk), INTENT(in) :: igrid
        INTEGER(intk), INTENT(in) :: istart, istop, jstart, jstop, &
            kstart, kstop

        ! Local variables
        INTEGER(intk) :: k, j, i
        REAL(realk), POINTER, CONTIGUOUS :: rarr(:, :, :)

        CALL field%get_ptr(rarr, igrid)
        DO i = istart, istop
            DO j = jstart, jstop
                DO k = kstart, kstop
                    icount = icount + 1
                    rarr(k, j, i) = recvbuf(icount)
                END DO
            END DO
        END DO
    END SUBROUTINE read_single_buffer


    ! Connect to self
    !
    ! Directly copes data from source to destination buffer
    SUBROUTINE connect_self(sendid, nplane, normal, flag, v1, v2, v3, &
            s1, s2, s3)
        ! Subroutine arguments
        INTEGER(int32), INTENT(in) :: sendid
        INTEGER(int32), INTENT(in) :: nplane
        LOGICAL, INTENT(in) :: normal
        CHARACTER(len=1), INTENT(in) :: flag
        TYPE(field_t), INTENT(inout), OPTIONAL :: v1, v2, v3, s1, s2, s3

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
        igrid_d = sendconns(3, sendid)
        igrid = sendconns(4, sendid)
        ifacerecv = sendconns(5, sendid)
        ifacesend = sendconns(6, sendid)

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
            exU = (normal .AND. ifacerecv < 3) .OR. (.NOT. normal)
        END IF
        IF (PRESENT(v1) .AND. exU) THEN
            CALL connect_self_single(v1, igrid, igrid_d, istart, istop, &
                jstart, jstop, kstart, kstop, istart_d, istop_d, &
                jstart_d, jstop_d, kstart_d, kstop_d)
        END IF

        IF (flag == 'W') THEN
            exV = (ifacerecv == 3)
        ELSE
            exV = (normal .AND. (ifacerecv > 2 .AND. ifacerecv < 5)) .OR. &
                (.NOT. normal)
        END IF
        IF (PRESENT(v2) .AND. exV) THEN
            CALL connect_self_single(v2, igrid, igrid_d, istart, istop, &
                jstart, jstop, kstart, kstop, istart_d, istop_d, &
                jstart_d, jstop_d, kstart_d, kstop_d)
        END IF

        IF (flag == 'W') THEN
            exW = (ifacerecv == 5)
        ELSE
            exW = (normal .AND. ifacerecv > 4) .OR. (.NOT. normal)
        END IF
        IF (PRESENT(v3) .AND. exW) THEN
            CALL connect_self_single(v3, igrid, igrid_d, istart, istop, &
                jstart, jstop, kstart, kstop, istart_d, istop_d, &
                jstart_d, jstop_d, kstart_d, kstop_d)
        END IF

        IF (flag == 'W') THEN
            exp1 = (ifacerecv == 2) .OR. (ifacerecv == 4) .OR. &
                (ifacerecv == 6)
        ELSE
            exp1 = .TRUE.
        END IF
        IF (PRESENT(s1) .AND. exp1) THEN
            CALL connect_self_single(s1, igrid, igrid_d, istart, istop, &
                jstart, jstop, kstart, kstop, istart_d, istop_d, &
                jstart_d, jstop_d, kstart_d, kstop_d)
        END IF

        IF (PRESENT(s2)) THEN
            CALL connect_self_single(s2, igrid, igrid_d, istart, istop, &
                jstart, jstop, kstart, kstop, istart_d, istop_d, &
                jstart_d, jstop_d, kstart_d, kstop_d)
        END IF

        IF (PRESENT(s3)) THEN
            CALL connect_self_single(s3, igrid, igrid_d, istart, istop, &
                jstart, jstop, kstart, kstop, istart_d, istop_d, &
                jstart_d, jstop_d, kstart_d, kstop_d)
        END IF
    END SUBROUTINE connect_self


    SUBROUTINE connect_self_single(field, igrid, igrid_d, istart, istop, &
            jstart, jstop, kstart, kstop, istart_d, istop_d, jstart_d, &
            jstop_d, kstart_d, kstop_d)

        ! Subroutine arguments
        TYPE(field_t), INTENT(inout) :: field
        INTEGER(intk), INTENT(in) :: igrid, igrid_d
        INTEGER(intk), INTENT(in) :: istart, istop, jstart, jstop, &
            kstart, kstop
        INTEGER(intk), INTENT(in) :: istart_d, istop_d, jstart_d, jstop_d, &
            kstart_d, kstop_d

        ! Local variables
        REAL(realk), POINTER, CONTIGUOUS :: src_rarr(:, :, :), dst_rarr(:, :, :)
        INTEGER(intk) :: i, j, k, ioff, joff, koff

        koff = kstart - kstart_d
        joff = jstart - jstart_d
        ioff = istart - istart_d

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
    END SUBROUTINE connect_self_single


    ! Process receive buffers as they arrive, wait for send
    ! buffers to be free
    SUBROUTINE process_bufs(nplane, normal, flag, v1, v2, v3, s1, s2, s3)
        ! Subroutine arguments
        INTEGER(int32), INTENT(in) :: nplane
        LOGICAL, INTENT(in) :: normal
        CHARACTER(len=1), INTENT(in) :: flag
        TYPE(field_t), OPTIONAL, INTENT(inout) :: &
            v1, v2, v3, s1, s2, s3

        ! Local variables
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

                unpackLen = 0
                DO i = 1, SIZE(recvidxlist, 2)
                    IF (recvidxlist(1, i) == recvlist(idx) &
                            .AND. recvidxlist(2, i) > 0) THEN
                        CALL read_buffer(i, nplane, normal, flag, v1, v2, v3, &
                            s1, s2, s3)
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


    SUBROUTINE init_conn()
        ! The maximum number of concurrent communications are the number
        ! of processes
        ALLOCATE(recvidxlist(3, SIZE(recvconns, 2)))
        ALLOCATE(sendlist(numprocs))
        ALLOCATE(recvlist(numprocs))
        ALLOCATE(sendreqs(numprocs))
        ALLOCATE(recvreqs(numprocs))
        recvidxlist = 0
        sendlist = 0
        recvlist = 0
        nrecv = 0
        nsend = 0
    END SUBROUTINE init_conn


    SUBROUTINE finish_conn()
        DEALLOCATE(recvidxlist)
        DEALLOCATE(sendlist)
        DEALLOCATE(recvlist)
        DEALLOCATE(sendreqs)
        DEALLOCATE(recvreqs)
    END SUBROUTINE finish_conn
END MODULE conn_mod
