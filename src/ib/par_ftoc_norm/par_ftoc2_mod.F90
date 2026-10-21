MODULE par_ftoc2_mod
    USE MPI_f08
    USE core_mod
    USE par_ftoc_core_mod

    IMPLICIT NONE (type, external)
    PRIVATE

    TYPE :: sendtask_t
        INTEGER(intk) :: fieldid
        INTEGER(intk) :: sumflag
        INTEGER(intk) :: sendidx
        INTEGER(intk) :: igridf
        INTEGER(intk) :: iface
        INTEGER(intk) :: kstart
        INTEGER(intk) :: kstop
        INTEGER(intk) :: jstart
        INTEGER(intk) :: jstop
        INTEGER(intk) :: istart
        INTEGER(intk) :: istop
    END TYPE sendtask_t

    TYPE :: recvtask_t
        INTEGER(intk) :: fieldid
        INTEGER(intk) :: recvidx
        INTEGER(intk) :: tasksize
        INTEGER(intk) :: igridc
        INTEGER(intk) :: kstart
        INTEGER(intk) :: kstop
        INTEGER(intk) :: jstart
        INTEGER(intk) :: jstop
        INTEGER(intk) :: istart
        INTEGER(intk) :: istop
    END TYPE recvtask_t

    TYPE :: selftask_t
        INTEGER(intk) :: fieldid
        INTEGER(intk) :: sumflag
        INTEGER(intk) :: igridf
        INTEGER(intk) :: igridc
        INTEGER(intk) :: iface
        INTEGER(intk) :: kstart
        INTEGER(intk) :: kstop
        INTEGER(intk) :: jstart
        INTEGER(intk) :: jstop
        INTEGER(intk) :: istart
        INTEGER(intk) :: istop
        INTEGER(intk) :: kpos
        INTEGER(intk) :: jpos
        INTEGER(intk) :: ipos
        INTEGER(intk) :: ifacerecv
    END TYPE selftask_t

    TYPE :: mpitask_t
        INTEGER(intk) :: iprocnbr
        INTEGER(intk) :: messagelength
        INTEGER(intk) :: counter
    END TYPE mpitask_t

    ! Lists that hold the send and receive request arrays
    TYPE(MPI_Request), ALLOCATABLE :: sendreqs(:), recvreqs(:)
    INTEGER(int32), ALLOCATABLE :: recvlist(:)
    INTEGER(intk), ALLOCATABLE :: recvidxlist(:, :)
    INTEGER(intk) :: nsend, nrecv

    ! Task extents
    INTEGER(intk) :: maxsendtasks, maxrecvtasks, maxselftasks
    INTEGER(intk) :: maxmpisendtasks, maxmpirecvtasks

    ! Pointers to fields passed to par_ftoc_norm, point to dummy field if not
    ! present
    TYPE(field_t), POINTER :: f1 => NULL(), f2 => NULL(), f3 => NULL()

    ! Workpackages containing individual tasks for packing / unpacking
    TYPE(sendtask_t), ALLOCATABLE :: sendtasks(:)
    TYPE(recvtask_t), ALLOCATABLE :: recvtasks(:)
    TYPE(selftask_t), ALLOCATABLE :: selftasks(:)
    TYPE(mpitask_t), ALLOCATABLE :: mpisendtasks(:), mpirecvtasks(:)
    !$omp declare target(sendtasks, recvtasks, selftasks)

    TYPE :: work_t
        LOGICAL :: is_init = .FALSE.
        TYPE(sendtask_t), ALLOCATABLE :: sendtasks(:)
        TYPE(recvtask_t), ALLOCATABLE :: recvtasks(:)
        TYPE(selftask_t), ALLOCATABLE :: selftasks(:)
        TYPE(mpitask_t), ALLOCATABLE :: mpisendtasks(:)
        TYPE(mpitask_t), ALLOCATABLE :: mpirecvtasks(:)
    END TYPE work_t

    TYPE(work_t), ALLOCATABLE, TARGET :: workrecords(:, :)
    LOGICAL :: is_recording = .FALSE.

    ! Variable to indicate if the required data structures have been created
    LOGICAL :: is_init = .FALSE.

    PUBLIC :: par_ftoc2, init_par_ftoc2, finish_par_ftoc2
CONTAINS
    SUBROUTINE par_ftoc2(ilevel, v1, v2, v3, sum)
        ! Subroutine arguments
        INTEGER(intk), INTENT(in) :: ilevel
        TYPE(field_t), TARGET, INTENT(inout) :: v1, v2, v3
        LOGICAL, OPTIONAL, INTENT(in) :: sum

        ! Local variables
        INTEGER(intk) :: idx_sum
        LOGICAL :: sum2

        IF (.NOT. is_init) CALL errr(__FILE__, __LINE__)
        IF (ilevel < minlevel .OR. ilevel > maxlevel) THEN
            CALL errr(__FILE__, __LINE__)
        END IF

        sum2 = .FALSE.
        IF (PRESENT(sum)) THEN
            sum2 = sum
        END IF
        idx_sum = MERGE(1, 0, sum2)
        f1 => v1
        f2 => v2
        f3 => v3

        ASSOCIATE(wptr => workrecords(ilevel, idx_sum))
            IF (is_recording) THEN
                CALL recording_pass(wptr, ilevel, sum2)
            ELSE
                IF (wptr%is_init) THEN
                    CALL recorded_par_ftoc2(wptr)
                ELSE
                    CALL jit_par_ftoc2(ilevel, sum2)
                END IF
            END IF
        END ASSOCIATE

        f1 => NULL()
        f2 => NULL()
        f3 => NULL()
        nrecv = 0
        nsend = 0
    END SUBROUTINE par_ftoc2


    SUBROUTINE jit_par_ftoc2(ilevel, sum)
        ! Subroutine arguments
        INTEGER(intk), INTENT(in) :: ilevel
        LOGICAL, INTENT(in) :: sum

        ! Local variables
        INTEGER(intk) :: nsendtasks, nrecvtasks
        INTEGER(intk) :: nselftasks
        INTEGER(intk) :: nmpisendtasks

        ! Manufacturing the workpackage just in-time and execute
        CALL prepare_tasks_all(sendtasks, nsendtasks, selftasks, nselftasks, &
            mpisendtasks, nmpisendtasks, ilevel, sum)

        !$omp target update to( &
        !$omp& sendtasks(1:nsendtasks+1), &
        !$omp& selftasks(1:nselftasks+1)) nowait

        CALL recv_mpi_all(ilevel)

        !$omp taskwait

        CALL process_sendtasks(nsendtasks, sendtasks)
        CALL process_mpisend(nmpisendtasks, mpisendtasks)
        CALL process_selftasks(nselftasks, selftasks)

        CALL prepare_recvtasks_all(recvtasks, nrecvtasks)

        !$omp target update to(recvtasks(1:nrecvtasks+1))

        CALL process_recvtasks(nrecvtasks, recvtasks)
    END SUBROUTINE jit_par_ftoc2


    SUBROUTINE recorded_par_ftoc2(wptr)
        ! Subroutine arguments
        TYPE(work_t), INTENT(in) :: wptr

        ! Local variables
        INTEGER(intk) :: nsendtasks, nrecvtasks
        INTEGER(intk) :: nselftasks
        INTEGER(intk) :: nmpisendtasks, nmpirecvtasks

        ! Using the (offloaded) workpackage from the workrecords array
        nmpirecvtasks = SIZE(wptr%mpirecvtasks) - 1
        nmpisendtasks = SIZE(wptr%mpisendtasks) - 1
        nsendtasks = SIZE(wptr%sendtasks) - 1
        nselftasks = SIZE(wptr%selftasks) - 1
        nrecvtasks = SIZE(wptr%recvtasks) - 1

        CALL process_mpirecv(nmpirecvtasks, wptr%mpirecvtasks)
        CALL process_sendtasks(nsendtasks, wptr%sendtasks)
        CALL process_mpisend(nmpisendtasks, wptr%mpisendtasks)
        CALL process_selftasks(nselftasks, wptr%selftasks)
        CALL MPI_Waitall(nrecv, recvreqs, MPI_STATUSES_IGNORE)
        CALL process_recvtasks(nrecvtasks, wptr%recvtasks)
        CALL MPI_Waitall(nsend, sendreqs, MPI_STATUSES_IGNORE)
    END SUBROUTINE recorded_par_ftoc2


    SUBROUTINE recording_pass(wptr, ilevel, sum)
        ! Subroutine arguments
        TYPE(work_t), INTENT(inout) :: wptr
        INTEGER(intk), INTENT(in) :: ilevel
        LOGICAL, INTENT(in) :: sum

        ! Local variables
        INTEGER(intk) :: nsendtasks, nrecvtasks, nselftasks
        INTEGER(intk) :: nmpisendtasks, nmpirecvtasks

        IF (wptr%is_init) THEN
            WRITE(*, *) "Combination already recorded."
            CALL errr(__FILE__, __LINE__)
        END IF

        ! It is necessary to execute one cycle with communication
        ! as otherwise many valuable checks are not possible

        CALL prepare_mpirecvtasks(mpirecvtasks, nmpirecvtasks, ilevel)
        CALL prepare_tasks_all(sendtasks, nsendtasks, selftasks, &
            nselftasks, mpisendtasks, nmpisendtasks, ilevel, sum)

        !$omp target update to( &
        !$omp& sendtasks(1:nsendtasks+1), &
        !$omp& selftasks(1:nselftasks+1))

        CALL process_mpirecv(nmpirecvtasks, mpirecvtasks)
        CALL process_sendtasks(nsendtasks, sendtasks)
        CALL process_mpisend(nmpisendtasks, mpisendtasks)
        CALL process_selftasks(nselftasks, selftasks)
        CALL prepare_recvtasks_all(recvtasks, nrecvtasks)

        !$omp target update to( &
        !$omp& recvtasks(1:nrecvtasks+1))
        CALL process_recvtasks(nrecvtasks, recvtasks)

        ! Allocate the workpackage arrays in the exact sizes
        ALLOCATE(wptr%sendtasks(nsendtasks+1))
        ALLOCATE(wptr%recvtasks(nrecvtasks+1))
        ALLOCATE(wptr%selftasks(nselftasks+1))
        ALLOCATE(wptr%mpisendtasks(nmpisendtasks+1))
        ALLOCATE(wptr%mpirecvtasks(nmpirecvtasks+1))

        ! Store the workpackage in the workrecords array and map
        wptr%sendtasks = sendtasks(1:nsendtasks+1)
        wptr%recvtasks = recvtasks(1:nrecvtasks+1)
        wptr%selftasks = selftasks(1:nselftasks+1)
        wptr%mpisendtasks = mpisendtasks(1:nmpisendtasks+1)
        wptr%mpirecvtasks = mpirecvtasks(1:nmpirecvtasks+1)

        !$omp target enter data map(to: &
        !$omp& wptr%sendtasks(1:nsendtasks+1), &
        !$omp& wptr%recvtasks(1:nrecvtasks+1), &
        !$omp& wptr%selftasks(1:nselftasks+1))

        ! Mark the workpackage as initialized
        wptr%is_init = .TRUE.
    END SUBROUTINE recording_pass


    SUBROUTINE prepare_mpirecvtasks(mpirtasks, nmpirtasks, ilevel)
        ! Subroutine arguments
        TYPE(mpitask_t), INTENT(inout) :: mpirtasks(maxmpirecvtasks)
        INTEGER(intk), INTENT(out) :: nmpirtasks
        INTEGER(intk), INTENT(in) :: ilevel

        ! Local variables
        INTEGER(intk) :: i, igridf, iface, iprocnbr, facearea
        INTEGER(intk) :: impirtask
        INTEGER(int32) :: recvcounter, messagelength

        recvcounter = 0
        messagelength = 0
        impirtask = 0
        nrecv = 0
        recvidxlist = -HUGE(1_intk)
        recvlist = 0

        DO i = 1, irecv
            iprocnbr = recvconns(1, i)
            igridf = recvconns(3, i)
            iface = recvconns(5, i)
            IF (iprocnbr /= myid .AND. ilevel == level(igridf)) THEN
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
                    CALL add_mpi_task(mpirtasks, impirtask, iprocnbr, &
                        messagelength, recvcounter, 'recv')
                ELSE IF (recvconns(1, i+1) /= iprocnbr) THEN
                    CALL add_mpi_task(mpirtasks, impirtask, iprocnbr, &
                        messagelength, recvcounter, 'recv')
                END IF
            END IF
        END DO

        nmpirtasks = impirtask
        ! Add a harmful dummy task at (nmpirtasks+1) for checking
        mpirtasks(nmpirtasks+1) = mpitask_t(-1, -1, -1)
    END SUBROUTINE prepare_mpirecvtasks


    SUBROUTINE recv_mpi_all(ilevel)
        ! Subroutine arguments
        INTEGER(intk), INTENT(in) :: ilevel

        ! Local variables
        INTEGER(intk) :: i, igridf, iface, iprocnbr, area
        INTEGER(int32) :: counter, messagelength

        counter = 0
        messagelength = 0
        nrecv = 0
        recvidxlist = -HUGE(1_intk)
        recvlist = 0

        DO i = 1, irecv
            iprocnbr = recvconns(1, i)
            igridf = recvconns(3, i)
            iface = recvconns(5, i)

            IF (iprocnbr == myid) CYCLE

            IF (ilevel == level(igridf)) THEN
                area = face_area(igridf, iface)
                recvidxlist(1, i) = iprocnbr
                recvidxlist(2, i) = area
                recvidxlist(3, i) = counter + messagelength
                messagelength = messagelength + area

                IF (counter + messagelength > idim_mg_bufs) THEN
                    CALL errr(__FILE__, __LINE__)
                END IF
            END IF

            IF (messagelength > 0) THEN
                IF (i == irecv) THEN
                    CALL post_recv(iprocnbr, messagelength, counter)
                ELSE IF (recvconns(1, i+1) /= iprocnbr) THEN
                    CALL post_recv(iprocnbr, messagelength, counter)
                END IF
            END IF
        END DO
    END SUBROUTINE recv_mpi_all


    SUBROUTINE post_recv(iprocnbr, messagelength, counter)
        ! Subroutine arguments
        INTEGER(intk), INTENT(in) :: iprocnbr
        INTEGER(int32), INTENT(inout) :: messagelength, counter

        ! Local variables
        ! none...

        nrecv = nrecv + 1
        recvlist(nrecv) = INT(iprocnbr, int32)

        !$omp target data use_device_addr(recvbuf)
        CALL MPI_Irecv(recvbuf(counter+1), messagelength, mglet_mpi_real, &
            iprocnbr, 1, MPI_COMM_WORLD, recvreqs(nrecv))
        !$omp end target data

        counter = counter + messagelength
        messagelength = 0
    END SUBROUTINE post_recv


    SUBROUTINE prepare_tasks_all(stasks, nstasks, etasks, netasks, &
            mpistasks, nmpisendtasks, ilevel, sum)
        ! Subroutine arguments
        TYPE(sendtask_t), INTENT(inout) :: stasks(maxsendtasks)
        INTEGER(intk), INTENT(out) :: nstasks
        TYPE(selftask_t), INTENT(inout) :: etasks(maxselftasks)
        INTEGER(intk), INTENT(out) :: netasks
        TYPE(mpitask_t), INTENT(inout) :: mpistasks(maxmpisendtasks)
        INTEGER(intk), INTENT(out) :: nmpisendtasks
        INTEGER(intk), INTENT(in) :: ilevel
        LOGICAL, INTENT(in) :: sum

        ! Local variables
        INTEGER(intk) :: i, igridf, iprocnbr, istask, ietask, impistask
        INTEGER(int32) :: messagelength, sendcounter

        CALL profile_range_push("prepare_tasks_all")

        sendcounter = 0
        messagelength = 0
        istask = 0
        ietask = 0
        impistask = 0
        nsend = 0
        DO i = 1, isend
            iprocnbr = sendconns(2, i)
            igridf = sendconns(3, i)
            IF (ilevel == level(igridf)) THEN
                IF (iprocnbr == myid) THEN
                    CALL prepare_selftask(etasks, ietask, i, sum)
                ELSE
                    CALL prepare_sendtask(stasks, istask, i, messagelength, &
                        sendcounter, sum)
                END IF
            END IF
            IF (messagelength > 0) THEN
                IF (i == isend) THEN
                    CALL add_mpi_task(mpistasks, impistask, iprocnbr, &
                        messagelength, sendcounter, 'send')
                ELSE IF (sendconns(2, i+1) /= iprocnbr) THEN
                    CALL add_mpi_task(mpistasks, impistask, iprocnbr, &
                        messagelength, sendcounter, 'send')
                END IF
            END IF
        END DO
        nstasks = istask
        netasks = ietask
        nmpisendtasks = impistask
        stasks(nstasks+1) = sendtask_t(-1, -1, -1, -1, -1, -1, -1, -1, &
            -1, -1, -1)
        etasks(netasks+1) = selftask_t(-1, -1, -1, -1, -1, -1, -1, -1, &
            -1, -1, -1, -1, -1, -1, -1)
        mpistasks(nmpisendtasks+1) = mpitask_t(-1, -1, -1)

        CALL profile_range_pop()
    END SUBROUTINE prepare_tasks_all


    SUBROUTINE prepare_sendtask(stasks, istask, sendid, messagelength, &
            sendcounter, sum)
        ! Subroutine arguments
        TYPE(sendtask_t), INTENT(inout) :: stasks(maxsendtasks)
        INTEGER(intk), INTENT(inout) :: istask
        INTEGER(intk), INTENT(in) :: sendid
        INTEGER(int32), INTENT(inout) :: messagelength
        INTEGER(int32), INTENT(in) :: sendcounter
        LOGICAL, INTENT(in) :: sum

        ! Local variables
        INTEGER(intk) :: igridf, iface, facearea
        INTEGER(intk) :: ks, ke, js, je, is, ie

        igridf = sendconns(3, sendid)
        iface = sendconns(5, sendid)
        facearea = face_area(igridf, iface)
        IF (sendcounter + messagelength + facearea > idim_mg_bufs) THEN
            CALL errr(__FILE__, __LINE__)
        END IF
        CALL start_and_stop(ks, ke, js, je, is, ie, igridf, iface)
        istask = istask + 1
        stasks(istask) = sendtask_t(face_field(iface), MERGE(1, 0, sum), &
            sendcounter+messagelength+1, igridf, iface, ks, ke, js, je, &
            is, ie)
        messagelength = messagelength + facearea
    END SUBROUTINE prepare_sendtask


    SUBROUTINE prepare_selftask(etasks, ietask, sendid, sum)
        ! Subroutine arguments
        TYPE(selftask_t), INTENT(inout) :: etasks(maxselftasks)
        INTEGER(intk), INTENT(inout) :: ietask
        INTEGER(intk), INTENT(in) :: sendid
        LOGICAL, INTENT(in) :: sum

        ! Local variables
        INTEGER(intk) :: igridf, igridc, iface, ifacerecv
        INTEGER(intk) :: ks, ke, js, je, is, ie
        INTEGER(intk) :: ksc, kec, jsc, jec, isc, iec

        igridf = sendconns(3, sendid)
        igridc = sendconns(4, sendid)
        iface = sendconns(5, sendid)
        ifacerecv = sendconns(6, sendid)
        CALL start_and_stop(ks, ke, js, je, is, ie, igridf, iface)
        CALL get_unpack_extents(ksc, kec, jsc, jec, isc, iec, igridf, &
            igridc, iface, ifacerecv)
        ietask = ietask + 1
        etasks(ietask) = selftask_t(face_field(iface), MERGE(1, 0, sum), &
            igridf, igridc, iface, ks, ke, js, je, is, ie, ksc, jsc, isc, &
            ifacerecv)
    END SUBROUTINE prepare_selftask


    SUBROUTINE add_mpi_task(mpitasks, impitask, iprocnbr, messagelength, &
            counter, type)
        ! Subroutine arguments
        TYPE(mpitask_t), INTENT(inout) :: mpitasks(*)
        INTEGER(intk), INTENT(inout) :: impitask
        INTEGER(intk), INTENT(in) :: iprocnbr
        INTEGER(int32), INTENT(inout) :: messagelength, counter
        CHARACTER(len=4), INTENT(in) :: type

        ! Local variables
        ! none...

        IF (type == 'send') THEN
            nsend = nsend + 1
        ELSE IF (type == 'recv') THEN
            nrecv = nrecv + 1
            recvlist(nrecv) = INT(iprocnbr, int32)
        ELSE
            CALL errr(__FILE__, __LINE__)
        END IF

        impitask = impitask + 1
        mpitasks(impitask)%iprocnbr = iprocnbr
        mpitasks(impitask)%messagelength = messagelength
        mpitasks(impitask)%counter = counter

        counter = counter + messagelength
        messagelength = 0
    END SUBROUTINE add_mpi_task


    SUBROUTINE process_mpirecv(nmpirtasks, mpirtasks)
        ! Subroutine arguments
        INTEGER(intk), INTENT(in) :: nmpirtasks
        TYPE(mpitask_t), INTENT(in) :: mpirtasks(nmpirtasks+1)

        ! Local variables
        INTEGER(int32) :: itask, iprocnbr, messagelength, recvcounter

        CALL profile_range_push("process_mpirecv")

        DO itask = 1, nmpirtasks
            iprocnbr = INT(mpirtasks(itask)%iprocnbr, kind=int32)
            messagelength = INT(mpirtasks(itask)%messagelength, kind=int32)
            recvcounter = INT(mpirtasks(itask)%counter, kind=int32)
            !$omp target data use_device_addr(recvbuf)
            CALL MPI_Irecv(recvbuf(recvcounter+1), messagelength, &
                mglet_mpi_real, iprocnbr, 1, MPI_COMM_WORLD, &
                recvreqs(itask))
            !$omp end target data
        END DO
        nrecv = nmpirtasks

        ! Safety check based on final dummy entry
        IF (mpirtasks(nmpirtasks+1)%iprocnbr /= -1 .OR. &
                mpirtasks(nmpirtasks+1)%messagelength /= -1 .OR. &
                mpirtasks(nmpirtasks+1)%counter /= -1) THEN
            WRITE(*, *) "Did not encounter the expected dummy task."
            CALL errr(__FILE__, __LINE__)
        END IF

        CALL profile_range_pop()
    END SUBROUTINE process_mpirecv


    SUBROUTINE process_mpisend(nmpistasks, mpistasks)
        ! Subroutine arguments
        INTEGER(intk), INTENT(in) :: nmpistasks
        TYPE(mpitask_t), INTENT(in) :: mpistasks(nmpistasks+1)

        ! Local variables
        INTEGER(intk) :: itask
        INTEGER(int32) :: iprocnbr, messagelength, sendcounter

        CALL profile_range_push("process_mpisend")

        DO itask = 1, nmpistasks
            iprocnbr = INT(mpistasks(itask)%iprocnbr, kind=int32)
            messagelength = INT(mpistasks(itask)%messagelength, kind=int32)
            sendcounter = INT(mpistasks(itask)%counter, kind=int32)
            !$omp target data use_device_addr(sendbuf)
            CALL MPI_Isend(sendbuf(sendcounter+1), messagelength, &
                mglet_mpi_real, iprocnbr, 1, MPI_COMM_WORLD, &
                sendreqs(itask))
            !$omp end target data
        END DO

        nsend = nmpistasks

        ! Safety check based on final dummy entry
        IF (nmpistasks < maxmpisendtasks) THEN
            IF (mpistasks(nmpistasks+1)%iprocnbr /= -1 .OR. &
                    mpistasks(nmpistasks+1)%messagelength /= -1 .OR. &
                    mpistasks(nmpistasks+1)%counter /= -1) THEN
                WRITE(*, *) "Did not encounter the expected dummy task."
                CALL errr(__FILE__, __LINE__)
            END IF
        END IF

        CALL profile_range_pop()
    END SUBROUTINE process_mpisend


    SUBROUTINE prepare_recvtasks_all(rtasks, nrtasks)
        ! Subroutine arguments
        TYPE(recvtask_t), INTENT(inout) :: rtasks(maxrecvtasks)
        INTEGER(intk), INTENT(out) :: nrtasks

        ! Local variables
        INTEGER(intk) :: irtask
        TYPE(MPI_Status) :: recvstatus
        INTEGER(int32) :: idx, i, recvmessagelen, unpacklen

        CALL profile_range_push("prepare_recvtasks_all")

        irtask = 0
        DO WHILE (.TRUE.)
            IF (nrecv == 0) EXIT
            CALL MPI_Waitany(nrecv, recvreqs, idx, recvstatus)
            IF (idx == MPI_UNDEFINED) EXIT
            CALL MPI_Get_count(recvstatus, mglet_mpi_real, recvmessagelen)
            unpacklen = 0
            DO i = 1, irecv
                IF (recvidxlist(1, i) == recvlist(idx) .AND. &
                        recvidxlist(2, i) > 0) THEN
                    CALL prepare_recvtask(rtasks, irtask, i)
                    unpacklen = unpacklen + recvidxlist(2, i)
                END IF
            END DO
            IF (recvmessagelen /= unpacklen) CALL errr(__FILE__, __LINE__)
        END DO
        CALL MPI_Waitall(nsend, sendreqs, MPI_STATUSES_IGNORE)
        nrtasks = irtask
        rtasks(nrtasks+1) = recvtask_t(-1, -1, -1, -1, -1, -1, -1, -1, &
            -1, -1)

        CALL profile_range_pop()
    END SUBROUTINE prepare_recvtasks_all


    SUBROUTINE prepare_recvtask(rtasks, irtask, recvid)
        ! Subroutine arguments
        TYPE(recvtask_t), INTENT(inout) :: rtasks(maxrecvtasks)
        INTEGER(intk), INTENT(inout) :: irtask
        INTEGER(intk), INTENT(in) :: recvid

        ! Local variables
        INTEGER(intk) :: igridf, igridc, iface, ifacerecv
        INTEGER(intk) :: ks, ke, js, je, is, ie

        igridf = recvconns(3, recvid)
        igridc = recvconns(4, recvid)
        iface = recvconns(5, recvid)
        ifacerecv = recvconns(6, recvid)
        CALL get_unpack_extents(ks, ke, js, je, is, ie, igridf, igridc, &
            iface, ifacerecv)
        irtask = irtask + 1
        rtasks(irtask) = recvtask_t(face_field(iface), &
            recvidxlist(3, recvid)+1, recvidxlist(2, recvid), igridc, &
            ks, ke, js, je, is, ie)
    END SUBROUTINE prepare_recvtask


    SUBROUTINE get_unpack_extents(ks, ke, js, je, is, ie, igridf, &
            igridc, iface, ifacerecv)
        ! Subroutine arguments
        INTEGER(intk), INTENT(out) :: ks, ke, js, je, is, ie
        INTEGER(intk), INTENT(in) :: igridf, igridc, iface, ifacerecv

        ! Local variables
        INTEGER(intk) :: kkf, jjf, iif, kkc, jjc, iic

        CALL get_mgdims(kkf, jjf, iif, igridf)
        is = iposition(igridf)
        js = jposition(igridf)
        ks = kposition(igridf)
        ie = is + (iif-4)/2 - 1
        je = js + (jjf-4)/2 - 1
        ke = ks + (kkf-4)/2 - 1
        IF (ifacerecv < 0) THEN
            SELECT CASE (iface)
            CASE (1)
                is = is - 1
                ie = is
            CASE (2)
                is = ie
            CASE (3)
                js = js - 1
                je = js
            CASE (4)
                js = je
            CASE (5)
                ks = ks - 1
                ke = ks
            CASE (6)
                ks = ke
            END SELECT
        ELSE
            CALL get_mgdims(kkc, jjc, iic, igridc)
            SELECT CASE (ifacerecv)
            CASE (1)
                is = 2
                ie = 2
            CASE (2)
                is = iic - 2
                ie = iic - 2
            CASE (3)
                js = 2
                je = 2
            CASE (4)
                js = jjc - 2
                je = jjc - 2
            CASE (5)
                ks = 2
                ke = 2
            CASE (6)
                ks = kkc - 2
                ke = kkc - 2
            END SELECT
        END IF
    END SUBROUTINE get_unpack_extents


    SUBROUTINE process_sendtasks(nstasks, stasks)
        ! Subroutine arguments
        INTEGER(intk), INTENT(in) :: nstasks
        TYPE(sendtask_t), INTENT(in) :: stasks(nstasks+1)

        ! Local variables
        TYPE(field_t), POINTER :: ddx, ddy, ddz

        IF (nstasks == 0) RETURN

        CALL profile_range_push("process_sendtasks")

        CALL get_field(ddx, "DDX")
        CALL get_field(ddy, "DDY")
        CALL get_field(ddz, "DDZ")

        CALL process_sendtasks_impl(nstasks, stasks, f1%arr, f2%arr, &
            f3%arr, ddx%arr, ddy%arr, ddz%arr, sendbuf)

        CALL profile_range_pop()

        ! Safety check based on final dummy entry
        IF (stasks(nstasks+1)%fieldid /= -1 .OR. &
                stasks(nstasks+1)%sumflag /= -1 .OR. &
                stasks(nstasks+1)%sendidx /= -1 .OR. &
                stasks(nstasks+1)%igridf /= -1 .OR. &
                stasks(nstasks+1)%iface /= -1 .OR. &
                stasks(nstasks+1)%kstart /= -1 .OR. &
                stasks(nstasks+1)%kstop /= -1 .OR. &
                stasks(nstasks+1)%jstart /= -1 .OR. &
                stasks(nstasks+1)%jstop /= -1 .OR. &
                stasks(nstasks+1)%istart /= -1 .OR. &
                stasks(nstasks+1)%istop /= -1) THEN
            WRITE(*, *) "Did not encounter the expected dummy task."
            CALL errr(__FILE__, __LINE__)
        END IF
    END SUBROUTINE process_sendtasks


    SUBROUTINE process_sendtasks_impl(nstasks, stasks, a1, a2, a3, &
            ddx, ddy, ddz, sbuf)
        ! Subroutine arguments
        INTEGER(intk), INTENT(in) :: nstasks
        TYPE(sendtask_t), INTENT(in) :: stasks(nstasks+1)
        REAL(realk), INTENT(in) :: a1(*), a2(*), a3(*)
        REAL(realk), INTENT(in) :: ddx(*), ddy(*), ddz(*)
        REAL(realk), INTENT(inout) :: sbuf(*)

        ! Local variables
        INTEGER(intk) :: itask, fieldid, sumflag, sendidx, tasksize
        INTEGER(intk) :: igridf, iface
        INTEGER(intk) :: istart, istop, jstart, jstop, kstart, kstop
        INTEGER(intk) :: kk, jj, ii, ip3, ipx, ipy, ipz

        !$omp target teams distribute private(itask, fieldid, sumflag, &
        !$omp& sendidx, tasksize, igridf, iface, istart, istop, jstart, jstop, &
        !$omp& kstart, kstop, kk, jj, ii, ip3, ipx, ipy, ipz)
        DO itask = 1, nstasks
            fieldid = stasks(itask)%fieldid
            sumflag = stasks(itask)%sumflag
            sendidx = stasks(itask)%sendidx
            igridf = stasks(itask)%igridf
            iface = stasks(itask)%iface
            kstart = stasks(itask)%kstart
            kstop = stasks(itask)%kstop
            jstart = stasks(itask)%jstart
            jstop = stasks(itask)%jstop
            istart = stasks(itask)%istart
            istop = stasks(itask)%istop

            CALL get_mgdims(kk, jj, ii, igridf)
            CALL get_ip3(ip3, igridf)
            CALL get_ip1x(ipx, igridf)
            CALL get_ip1y(ipy, igridf)
            CALL get_ip1z(ipz, igridf)

            tasksize = ((kstop-kstart)/2+1)*((jstop-jstart)/2+1)* &
                ((istop-istart)/2+1)

            !$omp parallel
            SELECT CASE (fieldid)
            CASE (1)
                CALL pack_face(kk, jj, ii, tasksize, a1(ip3), &
                    sbuf(sendidx), ddx(ipx), ddy(ipy), ddz(ipz), iface, &
                    sumflag, kstart, kstop, jstart, jstop, istart, istop)
            CASE (2)
                CALL pack_face(kk, jj, ii, tasksize, a2(ip3), &
                    sbuf(sendidx), ddx(ipx), ddy(ipy), ddz(ipz), iface, &
                    sumflag, kstart, kstop, jstart, jstop, istart, istop)
            CASE (3)
                CALL pack_face(kk, jj, ii, tasksize, a3(ip3), &
                    sbuf(sendidx), ddx(ipx), ddy(ipy), ddz(ipz), iface, &
                    sumflag, kstart, kstop, jstart, jstop, istart, istop)
#ifdef _MGLET_DEBUG_
            CASE DEFAULT
                CALL errr(__FILE__, __LINE__)
#endif
            END SELECT
            !$omp end parallel
        END DO
        !$omp end target teams distribute
    END SUBROUTINE process_sendtasks_impl


    SUBROUTINE pack_face(kk, jj, ii, tasksize, arr, sbuf, ddx, ddy, ddz, &
            iface, sumflag, kstart, kstop, jstart, jstop, istart, istop)
        !$omp declare target
        ! Subroutine arguments
        INTEGER(intk), INTENT(in) :: kk, jj, ii, tasksize
        REAL(realk), INTENT(in) :: arr(kk, jj, ii)
        REAL(realk), INTENT(out) :: sbuf(tasksize)
        REAL(realk), INTENT(in) :: ddx(ii), ddy(jj), ddz(kk)
        INTEGER(intk), INTENT(in) :: iface, sumflag
        INTEGER(intk), INTENT(in) :: kstart, kstop, jstart, jstop
        INTEGER(intk), INTENT(in) :: istart, istop

        ! Local variables
        ! none...

        SELECT CASE (iface)
        CASE (1, 2)
            IF (sumflag == 0) THEN
                CALL pack_face_average_x(kk, jj, ii, tasksize, arr, sbuf, &
                    ddy, ddz, kstart, kstop, jstart, jstop, istart, istop)
            ELSE
                CALL pack_face_sum_x(kk, jj, ii, tasksize, arr, sbuf, &
                    kstart, kstop, jstart, jstop, istart, istop)
            END IF
        CASE (3, 4)
            IF (sumflag == 0) THEN
                CALL pack_face_average_y(kk, jj, ii, tasksize, arr, sbuf, &
                    ddx, ddz, kstart, kstop, jstart, jstop, istart, istop)
            ELSE
                CALL pack_face_sum_y(kk, jj, ii, tasksize, arr, sbuf, &
                    kstart, kstop, jstart, jstop, istart, istop)
            END IF
        CASE (5, 6)
            IF (sumflag == 0) THEN
                CALL pack_face_average_z(kk, jj, ii, tasksize, arr, sbuf, &
                    ddx, ddy, kstart, kstop, jstart, jstop, istart, istop)
            ELSE
                CALL pack_face_sum_z(kk, jj, ii, tasksize, arr, sbuf, &
                    kstart, kstop, jstart, jstop, istart, istop)
            END IF
#ifdef _MGLET_DEBUG_
        CASE DEFAULT
            CALL errr(__FILE__, __LINE__)
#endif
        END SELECT
    END SUBROUTINE pack_face


    SUBROUTINE pack_face_average_x(kk, jj, ii, tasksize, arr, sbuf, ddy, &
            ddz, kstart, kstop, jstart, jstop, istart, istop)
        !$omp declare target
        ! Subroutine arguments
        INTEGER(intk), INTENT(in) :: kk, jj, ii, tasksize
        REAL(realk), INTENT(in) :: arr(kk, jj, ii)
        REAL(realk), INTENT(out) :: sbuf(tasksize)
        REAL(realk), INTENT(in) :: ddy(jj), ddz(kk)
        INTEGER(intk), INTENT(in) :: kstart, kstop
        INTEGER(intk), INTENT(in) :: jstart, jstop, istart, istop

        ! Local variables
        INTEGER(intk) :: k, j, i, idx, nk, nj
        REAL(realk) :: flux, area

        nk = (kstop-kstart)/2 + 1
        nj = (jstop-jstart)/2 + 1

        !$omp do collapse(3) private(i, j, k, idx, flux, area)
        DO i = istart, istop, 2
            DO j = jstart, jstop, 2
                DO k = kstart, kstop, 2
                    idx = 1 + (k-kstart)/2 + (j-jstart)/2*nk + &
                        (i-istart)/2*nj*nk
                    flux = arr(k, j, i)*ddy(j)*ddz(k) + &
                        arr(k, j+1, i)*ddy(j+1)*ddz(k) + &
                        arr(k+1, j, i)*ddy(j)*ddz(k+1) + &
                        arr(k+1, j+1, i)*ddy(j+1)*ddz(k+1)
                    area = (ddy(j)+ddy(j+1))*(ddz(k)+ddz(k+1))
                    sbuf(idx) = flux/area
                END DO
            END DO
        END DO
        !$omp end do
    END SUBROUTINE pack_face_average_x


    SUBROUTINE pack_face_sum_x(kk, jj, ii, tasksize, arr, sbuf, &
            kstart, kstop, jstart, jstop, istart, istop)
        !$omp declare target
        ! Subroutine arguments
        INTEGER(intk), INTENT(in) :: kk, jj, ii, tasksize
        REAL(realk), INTENT(in) :: arr(kk, jj, ii)
        REAL(realk), INTENT(out) :: sbuf(tasksize)
        INTEGER(intk), INTENT(in) :: kstart, kstop, jstart, jstop
        INTEGER(intk), INTENT(in) :: istart, istop

        ! Local variables
        INTEGER(intk) :: k, j, i, idx, nk, nj

        nk = (kstop-kstart)/2 + 1
        nj = (jstop-jstart)/2 + 1

        !$omp do collapse(3) private(i, j, k, idx)
        DO i = istart, istop, 2
            DO j = jstart, jstop, 2
                DO k = kstart, kstop, 2
                    idx = 1 + (k-kstart)/2 + (j-jstart)/2*nk + &
                        (i-istart)/2*nj*nk
                    sbuf(idx) = arr(k, j, i) + arr(k, j+1, i) + &
                        arr(k+1, j, i) + arr(k+1, j+1, i)
                END DO
            END DO
        END DO
        !$omp end do
    END SUBROUTINE pack_face_sum_x


    SUBROUTINE pack_face_average_y(kk, jj, ii, tasksize, arr, sbuf, ddx, &
            ddz, kstart, kstop, jstart, jstop, istart, istop)
        !$omp declare target
        ! Subroutine arguments
        INTEGER(intk), INTENT(in) :: kk, jj, ii, tasksize
        REAL(realk), INTENT(in) :: arr(kk, jj, ii)
        REAL(realk), INTENT(out) :: sbuf(tasksize)
        REAL(realk), INTENT(in) :: ddx(ii), ddz(kk)
        INTEGER(intk), INTENT(in) :: kstart, kstop
        INTEGER(intk), INTENT(in) :: jstart, jstop, istart, istop

        ! Local variables
        INTEGER(intk) :: k, j, i, idx, nk, nj
        REAL(realk) :: flux, area

        nk = (kstop-kstart)/2 + 1
        nj = (jstop-jstart)/2 + 1

        !$omp do collapse(3) private(i, j, k, idx, flux, area)
        DO i = istart, istop, 2
            DO j = jstart, jstop, 2
                DO k = kstart, kstop, 2
                    idx = 1 + (k-kstart)/2 + (j-jstart)/2*nk + &
                        (i-istart)/2*nj*nk
                    flux = arr(k, j, i)*ddx(i)*ddz(k) + &
                        arr(k, j, i+1)*ddx(i+1)*ddz(k) + &
                        arr(k+1, j, i)*ddx(i)*ddz(k+1) + &
                        arr(k+1, j, i+1)*ddx(i+1)*ddz(k+1)
                    area = (ddx(i)+ddx(i+1))*(ddz(k)+ddz(k+1))
                    sbuf(idx) = flux/area
                END DO
            END DO
        END DO
        !$omp end do
    END SUBROUTINE pack_face_average_y


    SUBROUTINE pack_face_sum_y(kk, jj, ii, tasksize, arr, sbuf, &
            kstart, kstop, jstart, jstop, istart, istop)
        !$omp declare target
        ! Subroutine arguments
        INTEGER(intk), INTENT(in) :: kk, jj, ii, tasksize
        REAL(realk), INTENT(in) :: arr(kk, jj, ii)
        REAL(realk), INTENT(out) :: sbuf(tasksize)
        INTEGER(intk), INTENT(in) :: kstart, kstop, jstart, jstop
        INTEGER(intk), INTENT(in) :: istart, istop

        ! Local variables
        INTEGER(intk) :: k, j, i, idx, nk, nj

        nk = (kstop-kstart)/2 + 1
        nj = (jstop-jstart)/2 + 1

        !$omp do collapse(3) private(i, j, k, idx)
        DO i = istart, istop, 2
            DO j = jstart, jstop, 2
                DO k = kstart, kstop, 2
                    idx = 1 + (k-kstart)/2 + (j-jstart)/2*nk + &
                        (i-istart)/2*nj*nk
                    sbuf(idx) = arr(k, j, i) + arr(k, j, i+1) + &
                        arr(k+1, j, i) + arr(k+1, j, i+1)
                END DO
            END DO
        END DO
        !$omp end do
    END SUBROUTINE pack_face_sum_y


    SUBROUTINE pack_face_average_z(kk, jj, ii, tasksize, arr, sbuf, ddx, &
            ddy, kstart, kstop, jstart, jstop, istart, istop)
        !$omp declare target
        ! Subroutine arguments
        INTEGER(intk), INTENT(in) :: kk, jj, ii, tasksize
        REAL(realk), INTENT(in) :: arr(kk, jj, ii)
        REAL(realk), INTENT(out) :: sbuf(tasksize)
        REAL(realk), INTENT(in) :: ddx(ii), ddy(jj)
        INTEGER(intk), INTENT(in) :: kstart, kstop
        INTEGER(intk), INTENT(in) :: jstart, jstop, istart, istop

        ! Local variables
        INTEGER(intk) :: k, j, i, idx, nk, nj
        REAL(realk) :: flux, area

        nk = (kstop-kstart)/2 + 1
        nj = (jstop-jstart)/2 + 1

        !$omp do collapse(3) private(i, j, k, idx, flux, area)
        DO i = istart, istop, 2
            DO j = jstart, jstop, 2
                DO k = kstart, kstop, 2
                    idx = 1 + (k-kstart)/2 + (j-jstart)/2*nk + &
                        (i-istart)/2*nj*nk
                    flux = arr(k, j, i)*ddx(i)*ddy(j) + &
                        arr(k, j+1, i)*ddx(i)*ddy(j+1) + &
                        arr(k, j, i+1)*ddx(i+1)*ddy(j) + &
                        arr(k, j+1, i+1)*ddx(i+1)*ddy(j+1)
                    area = (ddx(i)+ddx(i+1))*(ddy(j)+ddy(j+1))
                    sbuf(idx) = flux/area
                END DO
            END DO
        END DO
        !$omp end do
    END SUBROUTINE pack_face_average_z


    SUBROUTINE pack_face_sum_z(kk, jj, ii, tasksize, arr, sbuf, &
            kstart, kstop, jstart, jstop, istart, istop)
        !$omp declare target
        ! Subroutine arguments
        INTEGER(intk), INTENT(in) :: kk, jj, ii, tasksize
        REAL(realk), INTENT(in) :: arr(kk, jj, ii)
        REAL(realk), INTENT(out) :: sbuf(tasksize)
        INTEGER(intk), INTENT(in) :: kstart, kstop, jstart, jstop
        INTEGER(intk), INTENT(in) :: istart, istop

        ! Local variables
        INTEGER(intk) :: k, j, i, idx, nk, nj

        nk = (kstop-kstart)/2 + 1
        nj = (jstop-jstart)/2 + 1

        !$omp do collapse(3) private(i, j, k, idx)
        DO i = istart, istop, 2
            DO j = jstart, jstop, 2
                DO k = kstart, kstop, 2
                    idx = 1 + (k-kstart)/2 + (j-jstart)/2*nk + &
                        (i-istart)/2*nj*nk
                    sbuf(idx) = arr(k, j, i) + arr(k, j+1, i) + &
                        arr(k, j, i+1) + arr(k, j+1, i+1)
                END DO
            END DO
        END DO
        !$omp end do
    END SUBROUTINE pack_face_sum_z


    SUBROUTINE process_recvtasks(nrtasks, rtasks)
        ! Subroutine arguments
        INTEGER(intk), INTENT(in) :: nrtasks
        TYPE(recvtask_t), INTENT(in) :: rtasks(nrtasks+1)

        ! Local variables
        ! none...

        IF (nrtasks == 0) RETURN

        CALL profile_range_push("process_recvtasks")

        CALL process_recvtasks_impl(nrtasks, rtasks, f1%arr, f2%arr, f3%arr, &
            recvbuf)

        CALL profile_range_pop()

        ! Safety check based on final dummy entry
        IF (rtasks(nrtasks+1)%fieldid /= -1 .OR. &
                rtasks(nrtasks+1)%recvidx /= -1 .OR. &
                rtasks(nrtasks+1)%tasksize /= -1 .OR. &
                rtasks(nrtasks+1)%igridc /= -1 .OR. &
                rtasks(nrtasks+1)%kstart /= -1 .OR. &
                rtasks(nrtasks+1)%kstop /= -1 .OR. &
                rtasks(nrtasks+1)%jstart /= -1 .OR. &
                rtasks(nrtasks+1)%jstop /= -1 .OR. &
                rtasks(nrtasks+1)%istart /= -1 .OR. &
                rtasks(nrtasks+1)%istop /= -1) THEN
            WRITE(*, *) "Did not encounter the expected dummy task."
            CALL errr(__FILE__, __LINE__)
        END IF
    END SUBROUTINE process_recvtasks


    SUBROUTINE process_recvtasks_impl(nrtasks, rtasks, a1, a2, a3, rbuf)
        ! Subroutine arguments
        INTEGER(intk), INTENT(in) :: nrtasks
        TYPE(recvtask_t), INTENT(in) :: rtasks(nrtasks+1)
        REAL(realk), INTENT(inout) :: a1(*), a2(*), a3(*)
        REAL(realk), INTENT(in) :: rbuf(*)

        ! Local variables
        INTEGER(intk) :: itask, fieldid, recvidx, tasksize, igridc
        INTEGER(intk) :: istart, istop, jstart, jstop, kstart, kstop
        INTEGER(intk) :: kk, jj, ii, ip3

        !$omp target teams distribute private(itask, fieldid, recvidx, &
        !$omp& tasksize, igridc, istart, istop, jstart, jstop, kstart, &
        !$omp& kstop, kk, jj, ii, ip3)
        DO itask = 1, nrtasks
            fieldid = rtasks(itask)%fieldid
            recvidx = rtasks(itask)%recvidx
            tasksize = rtasks(itask)%tasksize
            igridc = rtasks(itask)%igridc
            kstart = rtasks(itask)%kstart
            kstop = rtasks(itask)%kstop
            jstart = rtasks(itask)%jstart
            jstop = rtasks(itask)%jstop
            istart = rtasks(itask)%istart
            istop = rtasks(itask)%istop

            CALL get_mgdims(kk, jj, ii, igridc)
            CALL get_ip3(ip3, igridc)

            !$omp parallel
            SELECT CASE (fieldid)
            CASE (1)
                CALL unpack_face(kk, jj, ii, tasksize, a1(ip3), rbuf(recvidx), &
                    kstart, kstop, jstart, jstop, istart, istop)
            CASE (2)
                CALL unpack_face(kk, jj, ii, tasksize, a2(ip3), rbuf(recvidx), &
                    kstart, kstop, jstart, jstop, istart, istop)
            CASE (3)
                CALL unpack_face(kk, jj, ii, tasksize, a3(ip3), rbuf(recvidx), &
                    kstart, kstop, jstart, jstop, istart, istop)
#ifdef _MGLET_DEBUG_
            CASE DEFAULT
                CALL errr(__FILE__, __LINE__)
#endif
            END SELECT
            !$omp end parallel
        END DO
        !$omp end target teams distribute
    END SUBROUTINE process_recvtasks_impl


    SUBROUTINE unpack_face(kk, jj, ii, tasksize, arr, rbuf, kstart, &
            kstop, jstart, jstop, istart, istop)
        !$omp declare target
        ! Subroutine arguments
        INTEGER(intk), INTENT(in) :: kk, jj, ii, tasksize
        REAL(realk), INTENT(inout) :: arr(kk, jj, ii)
        REAL(realk), INTENT(in) :: rbuf(tasksize)
        INTEGER(intk), INTENT(in) :: kstart, kstop, jstart, jstop
        INTEGER(intk), INTENT(in) :: istart, istop

        ! Local variables
        INTEGER(intk) :: k, j, i, idx, nk, nj

        nk = kstop-kstart+1
        nj = jstop-jstart+1

        !$omp do collapse(3) private(i, j, k, idx)
        DO i = istart, istop
            DO j = jstart, jstop
                DO k = kstart, kstop
                    idx = 1 + (k-kstart) + (j-jstart)*nk + (i-istart)*nj*nk
                    arr(k, j, i) = rbuf(idx)
                END DO
            END DO
        END DO
        !$omp end do

#ifdef _MGLET_DEBUG_
        IF (nk*nj*(istop-istart+1) /= tasksize) THEN
            CALL errr(__FILE__, __LINE__)
        END IF
#endif
    END SUBROUTINE unpack_face


    SUBROUTINE process_selftasks(netasks, etasks)
        ! Subroutine arguments
        INTEGER(intk), INTENT(in) :: netasks
        TYPE(selftask_t), INTENT(in) :: etasks(netasks+1)

        ! Local variables
        TYPE(field_t), POINTER :: ddx, ddy, ddz

        IF (netasks == 0) RETURN

        CALL profile_range_push("process_selftasks")

        CALL get_field(ddx, "DDX")
        CALL get_field(ddy, "DDY")
        CALL get_field(ddz, "DDZ")

        CALL process_selftasks_impl(netasks, etasks, f1%arr, f2%arr, &
            f3%arr, ddx%arr, ddy%arr, ddz%arr)

        CALL profile_range_pop()

        ! Safety check based on final dummy entry
        IF (etasks(netasks+1)%fieldid /= -1 .OR. &
                etasks(netasks+1)%sumflag /= -1 .OR. &
                etasks(netasks+1)%igridf /= -1 .OR. &
                etasks(netasks+1)%igridc /= -1 .OR. &
                etasks(netasks+1)%iface /= -1 .OR. &
                etasks(netasks+1)%kstart /= -1 .OR. &
                etasks(netasks+1)%kstop /= -1 .OR. &
                etasks(netasks+1)%jstart /= -1 .OR. &
                etasks(netasks+1)%jstop /= -1 .OR. &
                etasks(netasks+1)%istart /= -1 .OR. &
                etasks(netasks+1)%istop /= -1 .OR. &
                etasks(netasks+1)%kpos /= -1 .OR. &
                etasks(netasks+1)%jpos /= -1 .OR. &
                etasks(netasks+1)%ipos /= -1 .OR. &
                etasks(netasks+1)%ifacerecv /= -1) THEN
            WRITE(*, *) "Did not encounter the expected dummy task."
            CALL errr(__FILE__, __LINE__)
        END IF
    END SUBROUTINE process_selftasks


    SUBROUTINE process_selftasks_impl(netasks, etasks, a1, a2, a3, &
            ddx, ddy, ddz)
        ! Subroutine arguments
        INTEGER(intk), INTENT(in) :: netasks
        TYPE(selftask_t), INTENT(in) :: etasks(netasks+1)
        REAL(realk), INTENT(inout) :: a1(*), a2(*), a3(*)
        REAL(realk), INTENT(in) :: ddx(*), ddy(*), ddz(*)

        ! Local variables
        INTEGER(intk) :: itask, fieldid, sumflag, igridf, igridc, iface
        INTEGER(intk) :: istart, istop, jstart, jstop, kstart, kstop
        INTEGER(intk) :: ipos, jpos, kpos
        INTEGER(intk) :: kkf, jjf, iif, kkc, jjc, iic
        INTEGER(intk) :: ip3f, ip3c, ipx, ipy, ipz

        !$omp target teams distribute private(itask, fieldid, sumflag, &
        !$omp& igridf, igridc, iface, istart, istop, jstart, jstop, kstart, &
        !$omp& kstop, ipos, jpos, kpos, kkf, jjf, iif, kkc, jjc, iic, ip3f, &
        !$omp& ip3c, ipx, ipy, ipz)
        DO itask = 1, netasks
            fieldid = etasks(itask)%fieldid
            sumflag = etasks(itask)%sumflag
            igridf = etasks(itask)%igridf
            igridc = etasks(itask)%igridc
            iface = etasks(itask)%iface
            kstart = etasks(itask)%kstart
            kstop = etasks(itask)%kstop
            jstart = etasks(itask)%jstart
            jstop = etasks(itask)%jstop
            istart = etasks(itask)%istart
            istop = etasks(itask)%istop
            kpos = etasks(itask)%kpos
            jpos = etasks(itask)%jpos
            ipos = etasks(itask)%ipos

            CALL get_mgdims(kkf, jjf, iif, igridf)
            CALL get_mgdims(kkc, jjc, iic, igridc)
            CALL get_ip3(ip3f, igridf)
            CALL get_ip3(ip3c, igridc)
            CALL get_ip1x(ipx, igridf)
            CALL get_ip1y(ipy, igridf)
            CALL get_ip1z(ipz, igridf)

            !$omp parallel
            SELECT CASE (fieldid)
            CASE (1)
                CALL restrict_self(a1, ddx(ipx), ddy(ipy), ddz(ipz), &
                    kkf, jjf, iif, kkc, jjc, iic, ip3f, ip3c, iface, &
                    sumflag, kstart, kstop, jstart, jstop, istart, istop, &
                    kpos, jpos, ipos)
            CASE (2)
                CALL restrict_self(a2, ddx(ipx), ddy(ipy), ddz(ipz), &
                    kkf, jjf, iif, kkc, jjc, iic, ip3f, ip3c, iface, &
                    sumflag, kstart, kstop, jstart, jstop, istart, istop, &
                    kpos, jpos, ipos)
            CASE (3)
                CALL restrict_self(a3, ddx(ipx), ddy(ipy), ddz(ipz), &
                    kkf, jjf, iif, kkc, jjc, iic, ip3f, ip3c, iface, &
                    sumflag, kstart, kstop, jstart, jstop, istart, istop, &
                    kpos, jpos, ipos)
#ifdef _MGLET_DEBUG_
            CASE DEFAULT
                CALL errr(__FILE__, __LINE__)
#endif
            END SELECT
            !$omp end parallel
        END DO
        !$omp end target teams distribute
    END SUBROUTINE process_selftasks_impl


    SUBROUTINE restrict_self(arr, ddx, ddy, ddz, kkf, jjf, iif, &
            kkc, jjc, iic, ip3f, ip3c, iface, sumflag, kstart, kstop, &
            jstart, jstop, istart, istop, kpos, jpos, ipos)
        !$omp declare target
        ! Subroutine arguments
        REAL(realk), INTENT(inout) :: arr(*)
        REAL(realk), INTENT(in) :: ddx(*), ddy(*), ddz(*)
        INTEGER(intk), INTENT(in) :: kkf, jjf, iif, kkc, jjc, iic
        INTEGER(intk), INTENT(in) :: ip3f, ip3c, iface, sumflag
        INTEGER(intk), INTENT(in) :: kstart, kstop, jstart, jstop
        INTEGER(intk), INTENT(in) :: istart, istop, kpos, jpos, ipos

        ! Local variables
        ! none...

        SELECT CASE (iface)
        CASE (1, 2)
            IF (sumflag == 0) THEN
                CALL restrict_self_average_x(kkf, jjf, iif, kkc, jjc, &
                    iic, arr(ip3f), arr(ip3c), ddy, ddz, kstart, kstop, &
                    jstart, jstop, istart, istop, kpos, jpos, ipos)
            ELSE
                CALL restrict_self_sum_x(kkf, jjf, iif, kkc, jjc, iic, &
                    arr(ip3f), arr(ip3c), kstart, kstop, jstart, jstop, &
                    istart, istop, kpos, jpos, ipos)
            END IF
        CASE (3, 4)
            IF (sumflag == 0) THEN
                CALL restrict_self_average_y(kkf, jjf, iif, kkc, jjc, &
                    iic, arr(ip3f), arr(ip3c), ddx, ddz, kstart, kstop, &
                    jstart, jstop, istart, istop, kpos, jpos, ipos)
            ELSE
                CALL restrict_self_sum_y(kkf, jjf, iif, kkc, jjc, iic, &
                    arr(ip3f), arr(ip3c), kstart, kstop, jstart, jstop, &
                    istart, istop, kpos, jpos, ipos)
            END IF
        CASE (5, 6)
            IF (sumflag == 0) THEN
                CALL restrict_self_average_z(kkf, jjf, iif, kkc, jjc, &
                    iic, arr(ip3f), arr(ip3c), ddx, ddy, kstart, kstop, &
                    jstart, jstop, istart, istop, kpos, jpos, ipos)
            ELSE
                CALL restrict_self_sum_z(kkf, jjf, iif, kkc, jjc, iic, &
                    arr(ip3f), arr(ip3c), kstart, kstop, jstart, jstop, &
                    istart, istop, kpos, jpos, ipos)
            END IF
#ifdef _MGLET_DEBUG_
        CASE DEFAULT
            CALL errr(__FILE__, __LINE__)
#endif
        END SELECT
    END SUBROUTINE restrict_self


    SUBROUTINE restrict_self_average_x(kkf, jjf, iif, kkc, jjc, iic, &
            arrf, arrc, ddy, ddz, kstart, kstop, jstart, jstop, istart, &
            istop, kpos, jpos, ipos)
        !$omp declare target
        ! Subroutine arguments
        INTEGER(intk), INTENT(in) :: kkf, jjf, iif, kkc, jjc, iic
        REAL(realk), INTENT(in) :: arrf(kkf, jjf, iif)
        REAL(realk), INTENT(inout) :: arrc(kkc, jjc, iic)
        REAL(realk), INTENT(in) :: ddy(jjf), ddz(kkf)
        INTEGER(intk), INTENT(in) :: kstart, kstop, jstart, jstop
        INTEGER(intk), INTENT(in) :: istart, istop, kpos, jpos, ipos

        ! Local variables
        INTEGER(intk) :: k, j, i, kc, jc, ic
        REAL(realk) :: flux, area

        !$omp do collapse(3) private(i, j, k, ic, jc, kc, flux, area)
        DO i = istart, istop, 2
            DO j = jstart, jstop, 2
                DO k = kstart, kstop, 2
                    ic = ipos + (i-istart)/2
                    jc = jpos + (j-jstart)/2
                    kc = kpos + (k-kstart)/2
                    flux = arrf(k, j, i)*ddy(j)*ddz(k) + &
                        arrf(k, j+1, i)*ddy(j+1)*ddz(k) + &
                        arrf(k+1, j, i)*ddy(j)*ddz(k+1) + &
                        arrf(k+1, j+1, i)*ddy(j+1)*ddz(k+1)
                    area = (ddy(j)+ddy(j+1))*(ddz(k)+ddz(k+1))
                    arrc(kc, jc, ic) = flux/area
                END DO
            END DO
        END DO
        !$omp end do
    END SUBROUTINE restrict_self_average_x


    SUBROUTINE restrict_self_sum_x(kkf, jjf, iif, kkc, jjc, iic, arrf, &
            arrc, kstart, kstop, jstart, jstop, istart, istop, kpos, jpos, &
            ipos)
        !$omp declare target
        ! Subroutine arguments
        INTEGER(intk), INTENT(in) :: kkf, jjf, iif, kkc, jjc, iic
        REAL(realk), INTENT(in) :: arrf(kkf, jjf, iif)
        REAL(realk), INTENT(inout) :: arrc(kkc, jjc, iic)
        INTEGER(intk), INTENT(in) :: kstart, kstop, jstart, jstop
        INTEGER(intk), INTENT(in) :: istart, istop, kpos, jpos, ipos

        ! Local variables
        INTEGER(intk) :: k, j, i, kc, jc, ic

        !$omp do collapse(3) private(i, j, k, ic, jc, kc)
        DO i = istart, istop, 2
            DO j = jstart, jstop, 2
                DO k = kstart, kstop, 2
                    ic = ipos + (i-istart)/2
                    jc = jpos + (j-jstart)/2
                    kc = kpos + (k-kstart)/2
                    arrc(kc, jc, ic) = arrf(k, j, i) + &
                        arrf(k, j+1, i) + arrf(k+1, j, i) + &
                        arrf(k+1, j+1, i)
                END DO
            END DO
        END DO
        !$omp end do
    END SUBROUTINE restrict_self_sum_x


    SUBROUTINE restrict_self_average_y(kkf, jjf, iif, kkc, jjc, iic, &
            arrf, arrc, ddx, ddz, kstart, kstop, jstart, jstop, istart, &
            istop, kpos, jpos, ipos)
        !$omp declare target
        ! Subroutine arguments
        INTEGER(intk), INTENT(in) :: kkf, jjf, iif, kkc, jjc, iic
        REAL(realk), INTENT(in) :: arrf(kkf, jjf, iif)
        REAL(realk), INTENT(inout) :: arrc(kkc, jjc, iic)
        REAL(realk), INTENT(in) :: ddx(iif), ddz(kkf)
        INTEGER(intk), INTENT(in) :: kstart, kstop, jstart, jstop
        INTEGER(intk), INTENT(in) :: istart, istop, kpos, jpos, ipos

        ! Local variables
        INTEGER(intk) :: k, j, i, kc, jc, ic
        REAL(realk) :: flux, area

        !$omp do collapse(3) private(i, j, k, ic, jc, kc, flux, area)
        DO i = istart, istop, 2
            DO j = jstart, jstop, 2
                DO k = kstart, kstop, 2
                    ic = ipos + (i-istart)/2
                    jc = jpos + (j-jstart)/2
                    kc = kpos + (k-kstart)/2
                    flux = arrf(k, j, i)*ddx(i)*ddz(k) + &
                        arrf(k, j, i+1)*ddx(i+1)*ddz(k) + &
                        arrf(k+1, j, i)*ddx(i)*ddz(k+1) + &
                        arrf(k+1, j, i+1)*ddx(i+1)*ddz(k+1)
                    area = (ddx(i)+ddx(i+1))*(ddz(k)+ddz(k+1))
                    arrc(kc, jc, ic) = flux/area
                END DO
            END DO
        END DO
        !$omp end do
    END SUBROUTINE restrict_self_average_y


    SUBROUTINE restrict_self_sum_y(kkf, jjf, iif, kkc, jjc, iic, arrf, &
            arrc, kstart, kstop, jstart, jstop, istart, istop, kpos, jpos, &
            ipos)
        !$omp declare target
        ! Subroutine arguments
        INTEGER(intk), INTENT(in) :: kkf, jjf, iif, kkc, jjc, iic
        REAL(realk), INTENT(in) :: arrf(kkf, jjf, iif)
        REAL(realk), INTENT(inout) :: arrc(kkc, jjc, iic)
        INTEGER(intk), INTENT(in) :: kstart, kstop, jstart, jstop
        INTEGER(intk), INTENT(in) :: istart, istop, kpos, jpos, ipos

        ! Local variables
        INTEGER(intk) :: k, j, i, kc, jc, ic

        !$omp do collapse(3) private(i, j, k, ic, jc, kc)
        DO i = istart, istop, 2
            DO j = jstart, jstop, 2
                DO k = kstart, kstop, 2
                    ic = ipos + (i-istart)/2
                    jc = jpos + (j-jstart)/2
                    kc = kpos + (k-kstart)/2
                    arrc(kc, jc, ic) = arrf(k, j, i) + &
                        arrf(k, j, i+1) + arrf(k+1, j, i) + &
                        arrf(k+1, j, i+1)
                END DO
            END DO
        END DO
        !$omp end do
    END SUBROUTINE restrict_self_sum_y


    SUBROUTINE restrict_self_average_z(kkf, jjf, iif, kkc, jjc, iic, &
            arrf, arrc, ddx, ddy, kstart, kstop, jstart, jstop, istart, &
            istop, kpos, jpos, ipos)
        !$omp declare target
        ! Subroutine arguments
        INTEGER(intk), INTENT(in) :: kkf, jjf, iif, kkc, jjc, iic
        REAL(realk), INTENT(in) :: arrf(kkf, jjf, iif)
        REAL(realk), INTENT(inout) :: arrc(kkc, jjc, iic)
        REAL(realk), INTENT(in) :: ddx(iif), ddy(jjf)
        INTEGER(intk), INTENT(in) :: kstart, kstop, jstart, jstop
        INTEGER(intk), INTENT(in) :: istart, istop, kpos, jpos, ipos

        ! Local variables
        INTEGER(intk) :: k, j, i, kc, jc, ic
        REAL(realk) :: flux, area

        !$omp do collapse(3) private(i, j, k, ic, jc, kc, flux, area)
        DO i = istart, istop, 2
            DO j = jstart, jstop, 2
                DO k = kstart, kstop, 2
                    ic = ipos + (i-istart)/2
                    jc = jpos + (j-jstart)/2
                    kc = kpos + (k-kstart)/2
                    flux = arrf(k, j, i)*ddx(i)*ddy(j) + &
                        arrf(k, j+1, i)*ddx(i)*ddy(j+1) + &
                        arrf(k, j, i+1)*ddx(i+1)*ddy(j) + &
                        arrf(k, j+1, i+1)*ddx(i+1)*ddy(j+1)
                    area = (ddx(i)+ddx(i+1))*(ddy(j)+ddy(j+1))
                    arrc(kc, jc, ic) = flux/area
                END DO
            END DO
        END DO
        !$omp end do
    END SUBROUTINE restrict_self_average_z


    SUBROUTINE restrict_self_sum_z(kkf, jjf, iif, kkc, jjc, iic, arrf, &
            arrc, kstart, kstop, jstart, jstop, istart, istop, kpos, jpos, &
            ipos)
        !$omp declare target
        ! Subroutine arguments
        INTEGER(intk), INTENT(in) :: kkf, jjf, iif, kkc, jjc, iic
        REAL(realk), INTENT(in) :: arrf(kkf, jjf, iif)
        REAL(realk), INTENT(inout) :: arrc(kkc, jjc, iic)
        INTEGER(intk), INTENT(in) :: kstart, kstop, jstart, jstop
        INTEGER(intk), INTENT(in) :: istart, istop, kpos, jpos, ipos

        ! Local variables
        INTEGER(intk) :: k, j, i, kc, jc, ic

        !$omp do collapse(3) private(i, j, k, ic, jc, kc)
        DO i = istart, istop, 2
            DO j = jstart, jstop, 2
                DO k = kstart, kstop, 2
                    ic = ipos + (i-istart)/2
                    jc = jpos + (j-jstart)/2
                    kc = kpos + (k-kstart)/2
                    arrc(kc, jc, ic) = arrf(k, j, i) + &
                        arrf(k, j+1, i) + arrf(k, j, i+1) + &
                        arrf(k, j+1, i+1)
                END DO
            END DO
        END DO
        !$omp end do
    END SUBROUTINE restrict_self_sum_z


    FUNCTION face_field(iface) RESULT(fieldid)
        ! Function arguments
        INTEGER(intk) :: fieldid
        INTEGER(intk), INTENT(in) :: iface

        ! Local variables
        ! none...

        SELECT CASE (iface)
        CASE (1, 2)
            fieldid = 1
        CASE (3, 4)
            fieldid = 2
        CASE (5, 6)
            fieldid = 3
        CASE DEFAULT
            CALL errr(__FILE__, __LINE__)
        END SELECT
    END FUNCTION face_field


    SUBROUTINE init_par_ftoc2()
        ! Subroutine arguments
        ! none...

        ! Local variables
        ! none...

        IF (is_init) CALL errr(__FILE__, __LINE__)

        ! Check if par_ftoc_core_mod provides necessary infrastructure
        IF (.NOT. is_par_ftoc_core_init) CALL errr(__FILE__, __LINE__)

        ALLOCATE(recvidxlist(3, irecv))
        ALLOCATE(recvlist(irecv))
        ALLOCATE(sendreqs(isend))
        ALLOCATE(recvreqs(irecv))
        nsend = 0
        nrecv = 0

        maxsendtasks = isend + 1
        maxrecvtasks = irecv + 1
        maxselftasks = isend + 1
        maxmpisendtasks = isend + 1
        maxmpirecvtasks = irecv + 1
        ALLOCATE(sendtasks(maxsendtasks))
        ALLOCATE(recvtasks(maxrecvtasks))
        ALLOCATE(selftasks(maxselftasks))
        ALLOCATE(mpisendtasks(maxmpisendtasks))
        ALLOCATE(mpirecvtasks(maxmpirecvtasks))
        !$omp target enter data map(always, to: sendtasks, recvtasks, selftasks)

        ALLOCATE(workrecords(minlevel:maxlevel, 0:1))

        is_init = .TRUE.

        ! Record relevant calls for later efficient reuse on device
        CALL run_recording_pass()
    END SUBROUTINE init_par_ftoc2


    SUBROUTINE run_recording_pass()
        ! Subroutine arguments
        ! none...

        ! Local variables
        TYPE(field_t) :: u, v, w
        INTEGER(intk) :: ilevel

        CALL u%init("PAR_FTOC2_U", istag=1)
        CALL v%init("PAR_FTOC2_V", jstag=1)
        CALL w%init("PAR_FTOC2_W", kstag=1)
        !$omp target enter data map(to: u%arr, v%arr, w%arr)
        is_recording = .TRUE.
        DO ilevel = minlevel, maxlevel
            CALL par_ftoc2(ilevel, u, v, w, .FALSE.)
            CALL par_ftoc2(ilevel, u, v, w, .TRUE.)
        END DO
        is_recording = .FALSE.
        !$omp target exit data map(delete: u%arr, v%arr, w%arr)
        CALL u%finish()
        CALL v%finish()
        CALL w%finish()
    END SUBROUTINE run_recording_pass


    SUBROUTINE finish_par_ftoc2()
        ! Subroutine arguments
        ! none...

        ! Local variables
        ! none...

        IF (.NOT. is_init) CALL errr(__FILE__, __LINE__)

        nsend = 0
        nrecv = 0
        DEALLOCATE(recvreqs)
        DEALLOCATE(sendreqs)
        DEALLOCATE(recvlist)
        DEALLOCATE(recvidxlist)

        !$omp target exit data map(always, delete: sendtasks, recvtasks, &
        !$omp& selftasks)
        DEALLOCATE(sendtasks)
        DEALLOCATE(recvtasks)
        DEALLOCATE(selftasks)

        DEALLOCATE(mpisendtasks)
        DEALLOCATE(mpirecvtasks)

        CALL purge_workrecords()
        DEALLOCATE(workrecords)

        is_init = .FALSE.
    END SUBROUTINE finish_par_ftoc2


    SUBROUTINE purge_workrecords()
        ! Subroutine arguments
        ! none...

        ! Local variables
        TYPE(work_t), POINTER :: records(:)
        INTEGER(intk) :: i

        records(1:SIZE(workrecords)) => workrecords
        DO i = 1, SIZE(records)
            IF (.NOT. records(i)%is_init) CYCLE

            !$omp target exit data map(delete: records(i)%sendtasks, &
            !$omp& records(i)%recvtasks, records(i)%selftasks)
            DEALLOCATE(records(i)%sendtasks)
            DEALLOCATE(records(i)%recvtasks)
            DEALLOCATE(records(i)%selftasks)
            DEALLOCATE(records(i)%mpisendtasks)
            DEALLOCATE(records(i)%mpirecvtasks)
            records(i)%is_init = .FALSE.
        END DO
    END SUBROUTINE purge_workrecords
END MODULE par_ftoc2_mod
