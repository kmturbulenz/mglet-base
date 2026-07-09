MODULE conn_v1_mod

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
    USE fieldhelper_mod
    USE roctxprofile_mod

    IMPLICIT NONE (type, external)
    PRIVATE

    ! Lists that hold the send and receive request arrays
    TYPE(MPI_Request), ALLOCATABLE :: sendreqs(:), recvreqs(:)
    INTEGER(int32), ALLOCATABLE :: sendlist(:), recvlist(:)
    INTEGER(intk), ALLOCATABLE :: recvidxlist(:, :)
    INTEGER(intk) :: nsend, nrecv, maxtasks

    INTEGER(int32), PARAMETER :: mpitasksize = 3
    INTEGER(int32), PARAMETER :: buffertasksize = 9
    INTEGER(int32), PARAMETER :: selftasksize = 15

    ! Workpackages containing individual tasks for packing / unpacking
    INTEGER(int32), ALLOCATABLE :: sendtasks(:, :), recvtasks(:, :), &
        selftasks(:, :), mpisendtasks(:, :)
    !$omp declare target(sendtasks, recvtasks, selftasks)

    PUBLIC :: conn1, init_conn1, finish_conn1

CONTAINS

    SUBROUTINE conn1(ilevel, layers, v1, v2, v3, s1, s2, s3, corners, normal, &
            forward, ityp)

        ! conn is a more compact version of connect aiming for CPU offloading.
        ! It will only connect real fields, not integer fields, and does
        ! not have the same capabilities as connect. Over time features
        ! can be transferred, however, initially the goal is to have a
        ! light and easy code for GPU offloading while keeping the
        ! traditional connect in place.

        ! Subroutine arguments
        INTEGER(intk), INTENT(in), OPTIONAL :: ilevel, layers
        TYPE(field_t), OPTIONAL, TARGET, INTENT(inout) :: &
            v1, v2, v3, s1, s2, s3
        LOGICAL, OPTIONAL, INTENT(in) :: corners, normal
        INTEGER(intk), OPTIONAL, INTENT(in) :: forward
        CHARACTER(len=1), OPTIONAL, INTENT(in) :: ityp

        ! Local variables
        INTEGER(intk) :: minconlvl, maxconlvl, nplane, fwd, nvars
        INTEGER(int32) :: nsendtasks, nselftasks, nmpisendtasks, nrecvtasks
        INTEGER(int8) :: errorcode
        LOGICAL :: vertices, normal2
        CHARACTER(len=1) :: flag

        TYPE(field_t), POINTER :: f1, f2, f3, f4, f5, f6

        ! SIMON: Alternatively, we can point them to some dummy field (!)
        f1 => NULL()
        f2 => NULL()
        f3 => NULL()
        f4 => NULL()
        f5 => NULL()
        f6 => NULL()

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
            f1 => v1
            nvars = nvars + 1
        END IF
        IF (PRESENT(v2)) THEN
            f2 => v2
            nvars = nvars + 1
        END IF
        IF (PRESENT(v3)) THEN
            f3 => v3
            nvars = nvars + 1
        END IF
        IF (PRESENT(s1)) THEN
            f4 => s1
            nvars = nvars + 1
        END IF
        IF (PRESENT(s2)) THEN
            f5 => s2
            nvars = nvars + 1
        END IF
        IF (PRESENT(s3)) THEN
            f6 => s3
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

        ! Prepare the sendtasks, selftasks and mpisendtasks arrays on HOST
        CALL roctxrangepush("prepare_sendtasks_all")
        CALL prepare_sendtasks_all(sendtasks, nsendtasks, &
            selftasks, nselftasks, mpisendtasks, nmpisendtasks, &
            minconlvl, maxconlvl, nplane, vertices, &
            normal2, fwd, flag, nvars, v1, v2, v3, s1, s2, s3)
        CALL roctxrangepop()

        !$omp target update to( &
        !$omp&  sendtasks(1:buffertasksize, 1:nsendtasks+1), &
        !$omp&  selftasks(1:selftasksize, 1:nselftasks+1)) nowait

        ! Posting all non-blocking MPI recv calls with device pointers on HOST
        CALL roctxrangepush("recv_mpi_all")
        CALL recv_mpi_all(minconlvl, maxconlvl, nplane, vertices, normal2, &
            fwd, flag, nvars)
        CALL roctxrangepop()

        ! Wait for the task arrays to be copied to the device
        !$omp taskwait

        ! Packing the send buffer on the DEVICE (kernel must finish)
        CALL roctxrangepush("process_sendtasks")
        CALL process_sendtasks(sendtasks, nsendtasks, f1, f2, f3, f4, f5, f6, &
            errorcode)
        CALL check_error(errorcode, "Error in process_sendtasks")
        CALL roctxrangepop()

        CALL roctxrangepush("send_mpi_all")
        ! Posting all non-blocking MPI send calls with device pointers on HOST
        CALL send_mpi_all(mpisendtasks, nmpisendtasks)
        CALL roctxrangepop()


        ! >>> MPI messages are now in flight, we can do other work...
        ! SIMON: We are here following a conservative variant, which allows
        ! safety checks. For performance reasons, one could already execute
        ! "prepare_recvtasks_all" BEFORE the messages are there. That is
        ! without any checks (!)


        ! Packing the send buffer using on the DEVICE (nowait)
        CALL roctxrangepush("process_selftasks")
        CALL process_selftasks(selftasks, nselftasks, f1, f2, f3, f4, f5, f6, &
            errorcode)
        CALL roctxrangepop()

        ! Prepare the recvtasks arrays on HOST (ensures MPI completion)
        CALL roctxrangepush("prepare_recvtasks_all")
        CALL prepare_recvtasks_all(recvtasks, nrecvtasks, &
            nplane, normal2, flag, v1, v2, v3, s1, s2, s3)
        CALL roctxrangepop()

        !$omp target update to( &
        !$omp&  recvtasks(1:buffertasksize, 1:nrecvtasks+1)) nowait

        !$omp taskwait

        CALL roctxrangepush("process_recvtasks")
        ! Unpacking the recv buffer on the DEVICE (kernel must finish)
        CALL process_recvtasks(recvtasks, nrecvtasks, f1, f2, f3, f4, f5, f6, &
            errorcode)
        CALL roctxrangepop()
        CALL check_error(errorcode, "Error in process_recvtasks")

        CALL stop_timer(150)

    END SUBROUTINE conn1


    SUBROUTINE init_conn1()
        ! The maximum number of concurrent communications are the number
        ! of processes
        ALLOCATE(recvidxlist(3, SIZE(recvconns, 2)))
        ALLOCATE(sendlist(numprocs))
        ALLOCATE(recvlist(numprocs))
        ALLOCATE(sendreqs(numprocs))
        ALLOCATE(recvreqs(numprocs))

        ! Allocating the workpackage arrays at their maximal size
        maxtasks = 6 * SIZE(sendconns, 2) + 1
        ALLOCATE(sendtasks(buffertasksize, maxtasks))
        ALLOCATE(recvtasks(buffertasksize, maxtasks))
        ALLOCATE(selftasks(selftasksize, maxtasks))
        !$omp target enter data map(always, to: sendtasks, recvtasks, selftasks)

        ALLOCATE(mpisendtasks(mpitasksize, SIZE(sendconns, 2)+1))

        recvidxlist = 0
        sendlist = 0
        recvlist = 0
        nrecv = 0
        nsend = 0
    END SUBROUTINE init_conn1


    SUBROUTINE finish_conn1()
        DEALLOCATE(recvidxlist)
        DEALLOCATE(sendlist)
        DEALLOCATE(recvlist)
        DEALLOCATE(sendreqs)
        DEALLOCATE(recvreqs)

        !$omp target exit data map(delete: sendtasks, recvtasks, selftasks)
        DEALLOCATE(sendtasks)
        DEALLOCATE(recvtasks)
        DEALLOCATE(selftasks)

        DEALLOCATE(mpisendtasks)

    END SUBROUTINE finish_conn1


    ! Perform all MPI Send calls
    SUBROUTINE prepare_sendtasks_all(sendtasks, nsendtasks, &
            selftasks, nselftasks, mpisendtasks, nmpisendtasks, &
            minconlvl, maxconlvl, nplane, vertices, &
            normal, fwd, flag, nvars, v1, v2, v3, s1, s2, s3)

        ! Subroutine arguments
        INTEGER(int32), INTENT(inout) :: sendtasks(buffertasksize, maxtasks)
        INTEGER(int32), INTENT(out) :: nsendtasks
        INTEGER(int32), INTENT(inout) :: selftasks(selftasksize, maxtasks)
        INTEGER(int32), INTENT(out) :: nselftasks
        INTEGER(int32), INTENT(inout) :: mpisendtasks(mpitasksize, maxtasks)
        INTEGER(int32), INTENT(out) :: nmpisendtasks

        INTEGER(intk), INTENT(in) :: minconlvl, maxconlvl, nplane
        LOGICAL, INTENT(in) :: vertices, normal
        INTEGER(intk), INTENT(in) :: fwd
        CHARACTER(len=1), INTENT(in) :: flag
        INTEGER(intk), INTENT(in) :: nvars
        TYPE(field_t), OPTIONAL, INTENT(inout) :: v1, v2, v3, s1, s2, s3

        ! Local variables
        INTEGER(intk) :: i, iprocnbr, igrid, ifacerecv, facearea
        INTEGER(intk) :: isendtask, iselftask
        LOGICAL :: exchange, geometry
        INTEGER(int32) :: sendcounter, messagelength, impisendtasks

        geometry = .FALSE.

        ! Pack all buffers and send data
        sendcounter = 0
        messagelength = 0
        nsend = 0

        ! Initializing the task counters to zero
        isendtask = 0
        iselftask = 0
        impisendtasks = 0

        DO i = 1, SIZE(sendconns, 2)

            ! Just-in-time decision whether to exchange data or not
            exchange = decide(i, sendconns, geometry, vertices, fwd, &
                minconlvl, maxconlvl)
            iprocnbr = sendconns(1, i)

            ! Intra-rank communication with direct copy between local grids
            IF (iprocnbr == myid .AND. exchange) THEN

                ! >>> adding entries to selftasks (no execution)
                ! (replaces "connect_self")
                CALL prepare_selftasks(selftasks, iselftask, i, nplane, &
                    normal, flag, v1, v2, v3, s1, s2, s3)
                CYCLE

            END IF

            ! Message is send via MPI and tasks for buffer filling as added
            IF (exchange) THEN

                igrid = sendconns(3, i)
                ifacerecv = sendconns(5, i)
                facearea = face_area(igrid, ifacerecv, nplane, flag)

                ! >>> adding entries to sendtasks (no execution)
                ! (replaces "write_buffer")
                CALL prepare_sendtasks(sendtasks, isendtask, i, messagelength, &
                    sendcounter, nplane, normal, flag, nvars, v1, v2, v3, &
                    s1, s2, s3)

                messagelength = messagelength + nvars*facearea

            END IF

            ! Check if we need to record an MPI task
            ! >>> adding entries to mpisendtasks (no execution)
            ! (replaces "post_send")
            IF (messagelength > 0) THEN
                IF (i == SIZE(sendconns, 2)) THEN
                    CALL add_mpisend_task(mpisendtasks, impisendtasks, &
                        iprocnbr, messagelength, sendcounter)
                ELSE IF (sendconns(1, i + 1) /= iprocnbr) THEN
                    CALL add_mpisend_task(mpisendtasks, impisendtasks, &
                        iprocnbr, messagelength, sendcounter)
                END IF
            END IF

        END DO

        ! Set the output task counters
        nsendtasks = isendtask
        nselftasks = iselftask
        nmpisendtasks = impisendtasks

        ! Add a harmful dummy task at (ntasks+1) to detect execution overshoot
        sendtasks(:, nsendtasks+1) = -1
        selftasks(:, nselftasks+1) = -1
        mpisendtasks(:, nmpisendtasks+1) = -1

    END SUBROUTINE prepare_sendtasks_all


    SUBROUTINE prepare_sendtasks(sendtasks, isendtask, sendid, messagelength, &
            sendcounter, nplane, normal, flag, nvars, v1, v2, v3, s1, s2, s3)

        ! Subroutine arguments
        INTEGER(int32), INTENT(inout) :: sendtasks(buffertasksize, maxtasks)
        INTEGER(int32), INTENT(inout) :: isendtask
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
        INTEGER(int32) :: fieldid

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
            fieldid = 1   ! (for v1)
            isendtask = isendtask + 1
            CALL add_single_task(sendtasks(:, isendtask), fieldid, icount, &
                igrid, istart, istop, jstart, jstop, kstart, kstop)
        END IF

        IF (flag == 'W') THEN
            exV = (ifacerecv == 3)
        ELSE
            exV = (normal .AND. (ifacerecv > 2 .AND. ifacerecv < 5)) .OR. &
                (.NOT. normal)
        END IF
        IF (PRESENT(v2) .AND. exV) THEN
            fieldid = 2   ! (for v2)
            isendtask = isendtask + 1
            CALL add_single_task(sendtasks(:, isendtask), fieldid, icount, &
                igrid, istart, istop, jstart, jstop, kstart, kstop)
        END IF

        IF (flag == 'W') THEN
            exW = (ifacerecv == 5)
        ELSE
            exW = (normal .AND. ifacerecv > 4) .OR. (.NOT. normal)
        END IF
        IF (PRESENT(v3) .AND. exW) THEN
            fieldid = 3   ! (for v3)
            isendtask = isendtask + 1
            CALL add_single_task(sendtasks(:, isendtask), fieldid, icount, &
                igrid, istart, istop, jstart, jstop, kstart, kstop)
        END IF

        IF (flag == 'W') THEN
            exp1 = (ifacerecv == 2) .OR. (ifacerecv == 4) .OR. &
                (ifacerecv == 6)
        ELSE
            exp1 = .TRUE.
        END IF
        IF (PRESENT(s1) .AND. exp1) THEN
            fieldid = 4   ! (for s1)
            isendtask = isendtask + 1
            CALL add_single_task(sendtasks(:, isendtask), fieldid, icount, &
                igrid, istart, istop, jstart, jstop, kstart, kstop)
        END IF

        IF (PRESENT(s2)) THEN
            fieldid = 5   ! (for s2)
            isendtask = isendtask + 1
            CALL add_single_task(sendtasks(:, isendtask), fieldid, icount, &
                igrid, istart, istop, jstart, jstop, kstart, kstop)
        END IF

        IF (PRESENT(s3)) THEN
            fieldid = 6   ! (for s3)
            isendtask = isendtask + 1
            CALL add_single_task(sendtasks(:, isendtask), fieldid, icount, &
                igrid, istart, istop, jstart, jstop, kstart, kstop)
        END IF

        ! Check that message length was calculated correctly
        IF (thismessagelength /= (icount - offset)) THEN
            write(*, *) "thismessagelength:", thismessagelength, &
                "icount:", icount
            CALL errr(__FILE__, __LINE__)
        END IF

    END SUBROUTINE prepare_sendtasks





    ! Perform all MPI Recv calls
    SUBROUTINE recv_mpi_all(minconlvl, maxconlvl, nplane, vertices,  &
            normal, fwd, flag, nvars)

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
    END SUBROUTINE recv_mpi_all


    ! Perform a single MPI Recv
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

        !$omp target data use_device_addr(recvbuf)
        CALL MPI_Irecv(recvbuf(recvcounter+1), messagelength, &
            mglet_mpi_real, iprocnbr, 1, MPI_COMM_WORLD, recvreqs(nrecv))
        !$omp end target data

        recvcounter = recvcounter + messagelength
        messagelength = 0
    END SUBROUTINE post_recv


    ! Perform all MPI Send calls
    SUBROUTINE send_mpi_all(mpisendtasks, nmpisendtasks)

        ! Subroutine arguments
        INTEGER(int32), INTENT(in) :: mpisendtasks(mpitasksize, maxtasks)
        INTEGER(int32), INTENT(in) :: nmpisendtasks

        ! Local variables
        INTEGER(int32) :: itask, iprocnbr, messagelength, sendcounter

        ! Iterate over task and post all non-blocking MPI send calls
        DO itask = 1, nmpisendtasks
            iprocnbr      = mpisendtasks(1, itask)
            messagelength = mpisendtasks(2, itask)
            sendcounter   = mpisendtasks(3, itask)

            ! replaces "post_mpi_send"
            !$omp target data use_device_addr(sendbuf)
            CALL MPI_Isend(sendbuf(sendcounter + 1), messagelength, &
                mglet_mpi_real, iprocnbr, 1, MPI_COMM_WORLD, &
                sendreqs(itask))
            !$omp end target data
        END DO

        ! Safety check based on final dummy entry
        IF (nmpisendtasks < maxtasks) THEN
            IF (.NOT. ALL(mpisendtasks(:, nmpisendtasks+1) == -1)) THEN
                WRITE(*, *) "Did not encounter the expected dummy task."
                CALL errr(__FILE__, __LINE__)
            END IF
        END IF

    END SUBROUTINE send_mpi_all



    ! Perform a single non-blocking MPI send call
    SUBROUTINE add_mpisend_task(mpisendtasks, impisendtask, iprocnbr, &
        messagelength, sendcounter)

        ! Subroutine arguments
        INTEGER(int32), INTENT(inout) :: mpisendtasks(mpitasksize, maxtasks)
        INTEGER(int32), INTENT(inout) :: impisendtask
        INTEGER(int32), INTENT(in) :: iprocnbr
        INTEGER(int32), INTENT(inout) :: messagelength
        INTEGER(int32), INTENT(inout) :: sendcounter

        ! Local variables
        ! none...

        nsend = nsend + 1
        sendlist(nsend) = iprocnbr

        ! Add the MPI send task to the mpisendtasks array
        impisendtask = impisendtask + 1

        mpisendtasks(1, impisendtask) = iprocnbr
        mpisendtasks(2, impisendtask) = messagelength
        mpisendtasks(3, impisendtask) = sendcounter

        sendcounter = sendcounter + messagelength
        messagelength = 0

    END SUBROUTINE add_mpisend_task


    SUBROUTINE process_sendtasks(sendtasks, nsendtasks, f1, f2, f3, &
            f4, f5, f6, errorcode)

        ! Subroutine arguments
        INTEGER(int32), INTENT(in) :: sendtasks(buffertasksize, maxtasks)
        INTEGER(int32), INTENT(in) :: nsendtasks
        TYPE(field_t), POINTER, INTENT(inout) :: f1, f2, f3, f4, f5, f6
        INTEGER(int8), INTENT(out) :: errorcode

        ! Local variables
        INTEGER(int32) :: itask, fieldid, icount, igrid, istart, istop, &
            jstart, jstop, kstart, kstop, jjl, kkl, idx
        INTEGER(intk) :: i, j, k
        TYPE(field_t), POINTER :: field
        REAL(realk), POINTER, CONTIGUOUS :: rarr(:, :, :)

        IF (nsendtasks == 0) THEN
            errorcode = 0
            RETURN
        END IF

        !$omp target teams distribute private(itask, fieldid, icount, &
        !$omp&  igrid, istart, istop, jstart, jstop, kstart, kstop, jjl, kkl, &
        !$omp&  idx, i, j, k, field, rarr) &
        !$omp& reduction(max: errorcode) map(from: errorcode)
        DO itask = 1, nsendtasks
            errorcode = 0
            ! Set variables from sendtasks workpackage
            fieldid = sendtasks(1, itask)
            icount = sendtasks(2, itask)
            igrid = sendtasks(3, itask)
            istart = sendtasks(4, itask)
            istop = sendtasks(5, itask)
            jstart = sendtasks(6, itask)
            jstop = sendtasks(7, itask)
            kstart = sendtasks(8, itask)
            kstop = sendtasks(9, itask)

            ! Assign the correct field pointer based on fieldid
            SELECT CASE (fieldid)
                CASE (1)
                    IF (.NOT. ASSOCIATED(f1)) THEN
                        errorcode = 1
                    ELSE
                        field => f1
                    END IF
                CASE (2)
                    IF (.NOT. ASSOCIATED(f2)) THEN
                        errorcode = 2
                    ELSE
                        field => f2
                    END IF
                CASE (3)
                    IF (.NOT. ASSOCIATED(f3)) THEN
                        errorcode = 3
                    ELSE
                        field => f3
                    END IF
                CASE (4)
                    IF (.NOT. ASSOCIATED(f4)) THEN
                        errorcode = 4
                    ELSE
                        field => f4
                    END IF
                CASE (5)
                    IF (.NOT. ASSOCIATED(f5)) THEN
                        errorcode = 5
                    ELSE
                        field => f5
                    END IF
                CASE (6)
                    IF (.NOT. ASSOCIATED(f6)) THEN
                        errorcode = 6
                    ELSE
                        field => f6
                    END IF
                CASE DEFAULT
                    errorcode = 7
            END SELECT

            IF (errorcode == 0) THEN
                ! The following replaces "write_single_buffer"
                CALL get_grid3_real(rarr, field, igrid)

                ! Dimensions of the subarray to be treated
                kkl = kstop - kstart + 1
                jjl = jstop - jstart + 1

                ! Fully parallelizable copy loop
                !$omp parallel do collapse(3) private(i, j, k, idx)
                DO i = istart, istop
                    DO j = jstart, jstop
                        DO k = kstart, kstop

                            ! Computation of count to avoid incremental
                            idx = 1 + (k - kstart) + (j - jstart)*kkl + &
                                (i - istart)*jjl*kkl + icount

                            sendbuf(idx) = rarr(k, j, i)

                        END DO
                    END DO
                END DO
                !$omp end parallel do
            END IF

            IF (sendtasks(1, nsendtasks+1) /= -1) THEN
                errorcode = 8
            END IF
        END DO
        !$omp end target teams distribute
    END SUBROUTINE process_sendtasks


    SUBROUTINE process_recvtasks(recvtasks, nrecvtasks, f1, f2, f3, &
            f4, f5, f6, errorcode)

        ! Subroutine arguments
        INTEGER(int32), INTENT(in) :: recvtasks(buffertasksize, maxtasks)
        INTEGER(int32), INTENT(in) :: nrecvtasks
        TYPE(field_t), POINTER, INTENT(inout) :: f1, f2, f3, f4, f5, f6
        INTEGER(int8), INTENT(out) :: errorcode

        ! Local variables
        INTEGER(int32) :: itask, fieldid, icount, igrid, istart, istop, &
            jstart, jstop, kstart, kstop, kkl, jjl, idx
        TYPE(field_t), POINTER :: field
        REAL(realk), POINTER, CONTIGUOUS :: rarr(:, :, :)
        INTEGER(intk) :: i, j, k

        IF (nrecvtasks == 0) THEN
            errorcode = 0
            RETURN
        END IF

        !$omp target teams distribute private(itask, fieldid, icount, &
        !$omp&  igrid, istart, istop, jstart, jstop, kstart, kstop, jjl, kkl, &
        !$omp&  idx, i, j, k, field, rarr) &
        !$omp& reduction(max: errorcode) map(from: errorcode)
        DO itask = 1, nrecvtasks

            errorcode = 0
            ! Set variables from recvtasks workpackage
            fieldid = recvtasks(1, itask)
            icount  = recvtasks(2, itask)
            igrid   = recvtasks(3, itask)
            istart  = recvtasks(4, itask)
            istop   = recvtasks(5, itask)
            jstart  = recvtasks(6, itask)
            jstop   = recvtasks(7, itask)
            kstart  = recvtasks(8, itask)
            kstop   = recvtasks(9, itask)

            ! Assign the correct field pointer based on fieldid
            SELECT CASE (fieldid)
                CASE (1)
                    IF (.NOT. ASSOCIATED(f1)) THEN
                        errorcode = 1
                    ELSE
                        field => f1
                    END IF
                CASE (2)
                    IF (.NOT. ASSOCIATED(f2)) THEN
                        errorcode = 2
                    ELSE
                        field => f2
                    END IF
                CASE (3)
                    IF (.NOT. ASSOCIATED(f3)) THEN
                        errorcode = 3
                    ELSE
                        field => f3
                    END IF
                CASE (4)
                    IF (.NOT. ASSOCIATED(f4)) THEN
                        errorcode = 4
                    ELSE
                        field => f4
                    END IF
                CASE (5)
                    IF (.NOT. ASSOCIATED(f5)) THEN
                        errorcode = 5
                    ELSE
                        field => f5
                    END IF
                CASE (6)
                    IF (.NOT. ASSOCIATED(f6)) THEN
                        errorcode = 6
                    ELSE
                        field => f6
                    END IF
                CASE DEFAULT
                    errorcode = 7
            END SELECT

            IF (errorcode == 0) THEN
                ! The following replaces "read_single_buffer"
                CALL get_grid3_real(rarr, field, igrid)

                kkl = kstop - kstart + 1
                jjl = jstop - jstart + 1

                ! Fully parallelizable copy loop
                !$omp parallel do collapse(3) private(i, j, k, idx)
                DO i = istart, istop
                    DO j = jstart, jstop
                        DO k = kstart, kstop

                            idx = 1 + (k - kstart) + (j - jstart)*kkl + &
                                (i - istart)*jjl*kkl + icount

                            rarr(k, j, i) = recvbuf(idx)

                        END DO
                    END DO
                END DO
                !$omp end parallel do
            END IF

            IF (recvtasks(1, nrecvtasks+1) /= -1) THEN
                errorcode = 8
            END IF
        END DO
        !$omp end target teams distribute

    END SUBROUTINE process_recvtasks


    SUBROUTINE process_selftasks(selftasks, nselftasks, f1, f2, f3, &
            f4, f5, f6, errorcode)

        ! Subroutine arguments
        INTEGER(int32), INTENT(in) :: selftasks(selftasksize, maxtasks)
        INTEGER(int32), INTENT(in) :: nselftasks
        TYPE(field_t), POINTER, INTENT(inout) :: f1, f2, f3, f4, f5, f6
        INTEGER(int8), INTENT(out) :: errorcode

        ! Local variables
        INTEGER(int32) :: itask, fieldid
        INTEGER(int32) :: igrid, igrid_d
        INTEGER(int32) :: istart, istop, jstart, jstop, kstart, kstop, &
            istart_d, istop_d, jstart_d, jstop_d, kstart_d, kstop_d
        INTEGER(int32) :: koff, joff, ioff
        INTEGER(int32) :: i, j, k
        TYPE(field_t), POINTER :: field
        REAL(realk), POINTER, CONTIGUOUS :: src_rarr(:, :, :), dst_rarr(:, :, :)

        IF (nselftasks == 0) THEN
            errorcode = 0
            RETURN
        END IF

        !$omp target teams distribute private(itask, fieldid, igrid, igrid_d, &
        !$omp&  istart, istop, jstart, jstop, kstart, kstop, &
        !$omp&  istart_d, istop_d, jstart_d, jstop_d, kstart_d, kstop_d, &
        !$omp&  koff, joff, ioff, i, j, k, field, src_rarr, dst_rarr) &
        !$omp&  reduction(max: errorcode) map(from: errorcode) nowait
        DO itask = 1, nselftasks
            errorcode = 0

            ! Set variables from selftasks workpackage
            fieldid  = selftasks(1, itask)
            igrid    = selftasks(2, itask)
            igrid_d  = selftasks(3, itask)
            istart   = selftasks(4, itask)
            istop    = selftasks(5, itask)
            jstart   = selftasks(6, itask)
            jstop    = selftasks(7, itask)
            kstart   = selftasks(8, itask)
            kstop    = selftasks(9, itask)
            istart_d = selftasks(10, itask)
            istop_d  = selftasks(11, itask)
            jstart_d = selftasks(12, itask)
            jstop_d  = selftasks(13, itask)
            kstart_d = selftasks(14, itask)
            kstop_d  = selftasks(15, itask)

            ! Assign the correct field pointer based on fieldid
            SELECT CASE (fieldid)
                CASE (1)
                    IF (.NOT. ASSOCIATED(f1)) THEN
                        errorcode = 1
                    ELSE
                        field => f1
                    END IF
                CASE (2)
                    IF (.NOT. ASSOCIATED(f2)) THEN
                        errorcode = 2
                    ELSE
                        field => f2
                    END IF
                CASE (3)
                    IF (.NOT. ASSOCIATED(f3)) THEN
                        errorcode = 3
                    ELSE
                        field => f3
                    END IF
                CASE (4)
                    IF (.NOT. ASSOCIATED(f4)) THEN
                        errorcode = 4
                    ELSE
                        field => f4
                    END IF
                CASE (5)
                    IF (.NOT. ASSOCIATED(f5)) THEN
                        errorcode = 5
                    ELSE
                        field => f5
                    END IF
                CASE (6)
                    IF (.NOT. ASSOCIATED(f6)) THEN
                        errorcode = 6
                    ELSE
                        field => f6
                    END IF
                CASE DEFAULT
                    errorcode = 7
            END SELECT

            IF (errorcode == 0) THEN
                ! The following replaces "connect_self_single"
                koff = kstart - kstart_d
                joff = jstart - jstart_d
                ioff = istart - istart_d

                CALL get_grid3_real(src_rarr, field, igrid)
                CALL get_grid3_real(dst_rarr, field, igrid_d)

                !$omp parallel do collapse(3) private(i, j, k)
                DO i = istart_d, istop_d
                    DO j = jstart_d, jstop_d
                        DO k = kstart_d, kstop_d
                            dst_rarr(k, j, i) = &
                                src_rarr(k + koff, j + joff, i + ioff)
                        END DO
                    END DO
                END DO
                !$omp end parallel do
            END IF

            IF (selftasks(1, nselftasks+1) /= -1) THEN
                errorcode = 8
            END IF
        END DO
        !$omp end target teams distribute

    END SUBROUTINE process_selftasks






    ! Process receive buffers as they arrive, wait for send
    ! buffers to be free
    SUBROUTINE prepare_recvtasks_all(recvtasks, nrecvtasks, nplane, normal, &
            flag, v1, v2, v3, s1, s2, s3)

        ! Subroutine arguments
        INTEGER(int32), INTENT(inout) :: recvtasks(buffertasksize, maxtasks)
        INTEGER(int32), INTENT(out) :: nrecvtasks
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
        INTEGER(int32) :: irecvtask

        irecvtask = 0

        DO WHILE (.TRUE.)
            IF (nrecv == 0) EXIT

            ! Get index of an already arrived message
            CALL MPI_Waitany(nrecv, recvreqs, idx, recvstatus)

            IF (idx /= MPI_UNDEFINED) THEN
                CALL MPI_Get_count(recvstatus, mglet_mpi_real, &
                    recvmessagelen)

                unpacklen = 0
                DO i = 1, SIZE(recvidxlist, 2)
                    IF (recvidxlist(1, i) == recvlist(idx) &
                            .AND. recvidxlist(2, i) > 0) THEN

                        ! >>> adding entries to recvtasks (no execution)
                        ! (replaces "read_buffer")
                        CALL prepare_recvtasks(recvtasks, irecvtask, &
                            i, nplane, normal, flag, v1, v2, v3, s1, s2, s3)

                        unpacklen = unpacklen + recvidxlist(2, i)
                    END IF
                END DO

                ! Security check -- not possible before message arrival
                IF (recvmessagelen /= unpacklen) THEN
                    CALL errr(__FILE__, __LINE__)
                END IF
            ELSE
                EXIT
            END IF
        END DO

        ! Enfore the completion of all send requests
        CALL MPI_Waitall(nsend, sendreqs, MPI_STATUSES_IGNORE)

        ! Setting the number of recvtasks to the actual number of tasks
        nrecvtasks = irecvtask

        ! Add a harmful dummy task at (ntasks+1) to detect execution overshoot
        recvtasks(:, nrecvtasks+1) = -1

    END SUBROUTINE prepare_recvtasks_all


    SUBROUTINE prepare_recvtasks(recvtasks, irecvtask, recvid, nplane, &
            normal, flag, v1, v2, v3, s1, s2, s3)

        ! Subroutine arguments
        INTEGER(int32), INTENT(inout) :: recvtasks(buffertasksize, maxtasks)
        INTEGER(int32), INTENT(inout) :: irecvtask
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
        INTEGER(int32) :: offset, icount, fieldid

        ! Flags to indicate exchange of U, V, W
        LOGICAL :: exU, exV, exW, exp1


        ! Array for the workpackage:
        ! -> wp(10, SIZE(recvconns, 2))

        ! Set variables from send table
        igrid = recvconns(3, recvid)
        ifacerecv = recvconns(5, recvid)

        ! ioperation = 0

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
            irecvtask = irecvtask + 1
            fieldid = 1   ! (for v1)
            CALL add_single_task(recvtasks(:, irecvtask), fieldid, icount, &
                igrid, istart, istop, jstart, jstop, kstart, kstop)
        END IF

        IF (flag == 'W') THEN
            exV = (ifacerecv == 3)
        ELSE
            exV = (normal .AND. (ifacerecv > 2 .AND. ifacerecv < 5)) .OR. &
                (.NOT. normal)
        END IF
        IF (PRESENT(v2) .AND. exV) THEN
            irecvtask = irecvtask + 1
            fieldid = 2   ! (for v2)
            CALL add_single_task(recvtasks(:, irecvtask), fieldid, icount, &
                igrid, istart, istop, jstart, jstop, kstart, kstop)
        END IF

        IF (flag == 'W') THEN
            exW = (ifacerecv == 5)
        ELSE
            exW = (normal .AND. ifacerecv > 4) .OR. (.NOT. normal)
        END IF
        IF (PRESENT(v3) .AND. exW) THEN
            irecvtask = irecvtask + 1
            fieldid = 3   ! (for v3)
            CALL add_single_task(recvtasks(:, irecvtask), fieldid, icount, &
                igrid, istart, istop, jstart, jstop, kstart, kstop)
        END IF

        IF (flag == 'W') THEN
            exp1 = (ifacerecv == 2) .OR. (ifacerecv == 4) .OR. &
                (ifacerecv == 6)
        ELSE
            exp1 = .TRUE.
        END IF
        IF (PRESENT(s1) .AND. exp1) THEN
            irecvtask = irecvtask + 1
            fieldid = 4   ! (for s1)
            CALL add_single_task(recvtasks(:, irecvtask), fieldid, icount, &
                igrid, istart, istop, jstart, jstop, kstart, kstop)
        END IF

        IF (PRESENT(s2)) THEN
            irecvtask = irecvtask + 1
            fieldid = 5   ! (for s2)
            CALL add_single_task(recvtasks(:, irecvtask), fieldid, icount, &
                igrid, istart, istop, jstart, jstop, kstart, kstop)
        END IF

        IF (PRESENT(s3)) THEN
            irecvtask = irecvtask + 1
            fieldid = 6   ! (for s3)
            CALL add_single_task(recvtasks(:, irecvtask), fieldid, icount, &
                igrid, istart, istop, jstart, jstop, kstart, kstop)
        END IF

        ! Check that message length is calculated correctly
        IF ((icount - offset) /= recvidxlist(2, recvid)) THEN
            WRITE(*, *) "icount:", icount, &
                "recvidxlist(2, recvid):", recvidxlist(2, recvid)
            CALL errr(__FILE__, __LINE__)
        END IF

    END SUBROUTINE prepare_recvtasks






    ! Routine to prepare the workpackage for intra-rank self connect
    !
    SUBROUTINE prepare_selftasks(selftasks, iselftask, sendid, nplane, &
            normal, flag, v1, v2, v3, s1, s2, s3)

        ! Subroutine arguments
        INTEGER(int32), INTENT(inout) :: selftasks(selftasksize, maxtasks)
        INTEGER(int32), INTENT(inout) :: iselftask
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
        INTEGER(int32) :: dest_size, source_size, fieldid

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
            fieldid = 1   ! (for v1)
            iselftask = iselftask + 1
            CALL add_self_task(selftasks(:, iselftask), fieldid, igrid, &
            igrid_d, istart, istop, jstart, jstop, kstart, kstop, istart_d, &
            istop_d, jstart_d, jstop_d, kstart_d, kstop_d)
        END IF

        IF (flag == 'W') THEN
            exV = (ifacerecv == 3)
        ELSE
            exV = (normal .AND. (ifacerecv > 2 .AND. ifacerecv < 5)) .OR. &
                (.NOT. normal)
        END IF
        IF (PRESENT(v2) .AND. exV) THEN
            fieldid = 2   ! (for v2)
            iselftask = iselftask + 1
            CALL add_self_task(selftasks(:, iselftask), fieldid, igrid, &
            igrid_d, istart, istop, jstart, jstop, kstart, kstop, istart_d, &
            istop_d, jstart_d, jstop_d, kstart_d, kstop_d)
        END IF

        IF (flag == 'W') THEN
            exW = (ifacerecv == 5)
        ELSE
            exW = (normal .AND. ifacerecv > 4) .OR. (.NOT. normal)
        END IF
        IF (PRESENT(v3) .AND. exW) THEN
            fieldid = 3   ! (for v3)
            iselftask = iselftask + 1
            CALL add_self_task(selftasks(:, iselftask), fieldid, igrid, &
            igrid_d, istart, istop, jstart, jstop, kstart, kstop, istart_d, &
            istop_d, jstart_d, jstop_d, kstart_d, kstop_d)
        END IF

        IF (flag == 'W') THEN
            exp1 = (ifacerecv == 2) .OR. (ifacerecv == 4) .OR. &
                (ifacerecv == 6)
        ELSE
            exp1 = .TRUE.
        END IF
        IF (PRESENT(s1) .AND. exp1) THEN
            fieldid = 4   ! (for s1)
            iselftask = iselftask + 1
            CALL add_self_task(selftasks(:, iselftask), fieldid, igrid, &
            igrid_d, istart, istop, jstart, jstop, kstart, kstop, istart_d, &
            istop_d, jstart_d, jstop_d, kstart_d, kstop_d)
        END IF

        IF (PRESENT(s2)) THEN
            fieldid = 5   ! (for s2)
            iselftask = iselftask + 1
            CALL add_self_task(selftasks(:, iselftask), fieldid, igrid, &
            igrid_d, istart, istop, jstart, jstop, kstart, kstop, istart_d, &
            istop_d, jstart_d, jstop_d, kstart_d, kstop_d)
        END IF

        IF (PRESENT(s3)) THEN
            fieldid = 6   ! (for s3)
            iselftask = iselftask + 1
            CALL add_self_task(selftasks(:, iselftask), fieldid, igrid, &
            igrid_d, istart, istop, jstart, jstop, kstart, kstop, istart_d, &
            istop_d, jstart_d, jstop_d, kstart_d, kstop_d)
        END IF
    END SUBROUTINE prepare_selftasks



    ! Core routine to add a single self-copy task to the workpackage
    !
    PURE SUBROUTINE add_self_task(selftask, fieldid, igrid, igrid_d, &
            istart, istop, jstart, jstop, kstart, kstop, istart_d, &
            istop_d, jstart_d, jstop_d, kstart_d, kstop_d)

        ! Subroutine arguments
        INTEGER(intk), INTENT(inout) :: selftask(selftasksize)
        INTEGER(intk), INTENT(in) :: fieldid, igrid, igrid_d
        INTEGER(intk), INTENT(in) :: istart, istop, jstart, jstop, kstart, &
            kstop, istart_d, istop_d, jstart_d, jstop_d, kstart_d, kstop_d

        ! Filling the task
        selftask(1) = fieldid
        selftask(2) = igrid
        selftask(3) = igrid_d
        ! for the source grid (igrid) = for packing
        selftask(4) = istart
        selftask(5) = istop
        selftask(6) = jstart
        selftask(7) = jstop
        selftask(8) = kstart
        selftask(9) = kstop
        ! for the destination grid (igrid_d) = for unpacking
        selftask(10) = istart_d
        selftask(11) = istop_d
        selftask(12) = jstart_d
        selftask(13) = jstop_d
        selftask(14) = kstart_d
        selftask(15) = kstop_d

    END SUBROUTINE add_self_task





    ! Core routine to add a single buffer-related task to the workpackage
    !
    PURE SUBROUTINE add_single_task(task, fieldid, icount, igrid, &
            istart, istop, jstart, jstop, kstart, kstop)

        ! Subroutine arguments
        INTEGER(intk), INTENT(out) :: task(9)
        INTEGER(intk), INTENT(in) :: fieldid
        INTEGER(intk), INTENT(inout) :: icount
        INTEGER(intk), INTENT(in) :: igrid
        INTEGER(intk), INTENT(in) :: istart, istop, jstart, jstop, &
            kstart, kstop

        ! Local variables
        INTEGER(intk) :: tasksize

        ! Filling the task
        task(1) = fieldid
        task(2) = icount
        task(3) = igrid
        task(4) = istart
        task(5) = istop
        task(6) = jstart
        task(7) = jstop
        task(8) = kstart
        task(9) = kstop

        ! Increment the icount by the number of elements in the task
        tasksize = (istop-istart+1)*(jstop-jstart+1)*(kstop-kstart+1)
        icount = icount + tasksize

    END SUBROUTINE add_single_task


    SUBROUTINE check_error(errorcode, msg)
        INTEGER(int8), INTENT(in) :: errorcode
        CHARACTER(*), INTENT(in) :: msg

        IF (errorcode /= 0) THEN
            WRITE(*, *) "ERROR:", TRIM(msg), " errorcode =", errorcode
            IF (errorcode == 1) THEN
                WRITE(*, *) "v1 is not present"
            ELSE IF (errorcode == 2) THEN
                WRITE(*, *) "v2 is not present"
            ELSE IF (errorcode == 3) THEN
                WRITE(*, *) "v3 is not present"
            ELSE IF (errorcode == 4) THEN
                WRITE(*, *) "s1 is not present"
            ELSE IF (errorcode == 5) THEN
                WRITE(*, *) "s2 is not present"
            ELSE IF (errorcode == 6) THEN
                WRITE(*, *) "s3 is not present"
            ELSE IF (errorcode == 7) THEN
                WRITE(*, *) "Invalid fieldid"
            ELSE IF (errorcode == 8) THEN
                WRITE(*, *) "Check by dummy task failed"
            END IF
            CALL errr(__FILE__, __LINE__)
        END IF
    END SUBROUTINE check_error

END MODULE conn_v1_mod
