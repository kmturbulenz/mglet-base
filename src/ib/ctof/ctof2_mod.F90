MODULE ctof2_mod

    USE core_mod
    USE ctof_core_mod
    USE MPI_f08

    IMPLICIT NONE (type, external)
    PRIVATE

    ! Lists that hold the send and receive request arrays
    TYPE(MPI_Request), ALLOCATABLE :: sendreqs(:), recvreqs(:)

    ! Lists that hold the messages that are ACTUALLY sent and received
    INTEGER(intk) :: nsend, nrecv
    INTEGER(int32), ALLOCATABLE :: recvlist(:)
    INTEGER(intk), ALLOCATABLE :: recvidxlist(:, :)

    ! Variable to indicate if the required data structures have been created
    LOGICAL :: is_init = .FALSE.

    ! Indicator for the recoding mode
    LOGICAL :: is_recording = .FALSE.

    ! Workpackages containing individual tasks for packing / unpacking
    INTEGER(intk), ALLOCATABLE :: sendtasks(:, :), recvtasks(:, :), &
        mpisendtasks(:, :), mpirecvtasks(:, :), selftasks(:, :)
    ! SIMON: I think it is not worthwhile to declare those "target"

    ! Type to hold condensed task arrays to execute a certain type of conn
    TYPE :: work_t
        LOGICAL :: is_init = .FALSE.
        INTEGER(intk), ALLOCATABLE :: sendtasks(:, :)
        INTEGER(intk), ALLOCATABLE :: recvtasks(:, :)
        INTEGER(intk), ALLOCATABLE :: selftasks(:, :)
        INTEGER(intk), ALLOCATABLE :: mpisendtasks(:, :)
        INTEGER(intk), ALLOCATABLE :: mpirecvtasks(:, :)
    END TYPE work_t

    ! Array to store instructions for different values of "ilevel"
    TYPE(work_t), ALLOCATABLE, TARGET :: workrecords(:)
    INTEGER(intk), PARAMETER :: mpitasksize = 3
    INTEGER(intk), PARAMETER :: sendtasksize = 9
    INTEGER(intk), PARAMETER :: recvtasksize = 2
    INTEGER(intk), PARAMETER :: selftasksize = 8

    ! Maximum size of the temporary arrays
    INTEGER(intk) :: maxsize

    ! Internal pointers to the fine and coarse field_t objects
    TYPE(field_t), POINTER :: ff, fc


    PUBLIC :: ctof2, init_ctof2, finish_ctof2

CONTAINS

    SUBROUTINE ctof2(ilevel, ff_f, fc_f)

        ! Subroutine arguments
        INTEGER(intk), INTENT(in) :: ilevel  ! Level of the *fine* side
        TYPE(field_t), TARGET, INTENT(inout) :: ff_f
        TYPE(field_t), TARGET, INTENT(in) :: fc_f

        ! Local variables
        INTEGER(intk) :: nmpirecvtasks, nmpisendtasks
        INTEGER(intk) :: nsendtasks, nrecvtasks, nselftasks

        ! Setting the internal pointers to fine and coarse field_t objects
        ff => ff_f
        fc => fc_f

        ! Looking up the workpackage for this level
        ASSOCIATE(wptr => workrecords(ilevel))

        IF (.NOT. is_recording) THEN

            ! During the execution phase, the workpackage is already initialized
            IF (.NOT. wptr%is_init) THEN
                CALL errr(__FILE__, __LINE__)
            END IF

#ifdef _MGLET_PROFILE_ANNOTATIONS_
            CALL profile_range_push("ctof")
#endif

            ! Obtaining numbers of tasks from the workpackage
            nmpirecvtasks = SIZE(wptr%mpirecvtasks, 2) - 1
            nmpisendtasks = SIZE(wptr%mpisendtasks, 2) - 1
            nsendtasks = SIZE(wptr%sendtasks, 2) - 1
            nrecvtasks = SIZE(wptr%recvtasks, 2) - 1
            nselftasks = SIZE(wptr%selftasks, 2) - 1

            CALL process_mpirecvtasks(wptr%mpirecvtasks, nmpirecvtasks)
            CALL process_sendtasks(wptr%sendtasks, nsendtasks)
            CALL process_mpisendtasks(wptr%mpisendtasks, nmpisendtasks)

            ! While MPI messages travel, local work is done
            CALL process_selftasks(wptr%selftasks, nselftasks)

            CALL MPI_Waitall(nrecv, recvreqs, MPI_STATUSES_IGNORE)
            CALL process_recvtasks(wptr%recvtasks, nrecvtasks)
            CALL MPI_Waitall(nsend, sendreqs, MPI_STATUSES_IGNORE)

#ifdef _MGLET_PROFILE_ANNOTATIONS_
            CALL profile_range_pop()
#endif

        ELSE

            ! During the recording phase, the workpackage is stored
            IF (wptr%is_init) THEN
                WRITE(*, *) "Repeated recording for ilevel = ", ilevel
                CALL errr(__FILE__, __LINE__)
            END IF

            ! Allocating oversized temporary arrays to record the tasks
            maxsize = noflevel(ilevel) + 1
            ALLOCATE(sendtasks(sendtasksize, maxsize))
            ALLOCATE(recvtasks(recvtasksize, maxsize))
            ALLOCATE(selftasks(selftasksize, maxsize))
            ALLOCATE(mpisendtasks(mpitasksize, maxsize))
            ALLOCATE(mpirecvtasks(mpitasksize, maxsize))


            CALL prepare_allsendtasks(&
                sendtasks, nsendtasks, &
                selftasks, nselftasks, &
                mpisendtasks, nmpisendtasks, ilevel)

            CALL prepare_mpirecvtasks(mpirecvtasks, nmpirecvtasks, ilevel)

            !$omp target enter data map(to: &
            !$omp&  sendtasks(1:sendtasksize, 1:nsendtasks+1))
            CALL process_sendtasks(sendtasks, nsendtasks)
            !$omp target exit data map(delete: &
            !$omp&  sendtasks(1:sendtasksize, 1:nsendtasks+1))

            CALL process_mpirecvtasks(mpirecvtasks, nmpirecvtasks)
            CALL process_mpisendtasks(mpisendtasks, nmpisendtasks)

            ! Includes waiting for MPI communication to finish
            CALL prepare_recvtasks(recvtasks, nrecvtasks)

            !$omp target enter data map(to: &
            !$omp&  recvtasks(1:recvtasksize, 1:nrecvtasks+1))
            CALL process_recvtasks(recvtasks, nrecvtasks)
            !$omp target exit data map(delete: &
            !$omp&  recvtasks(1:recvtasksize, 1:nrecvtasks+1))

            !$omp target enter data map(to: &
            !$omp&  selftasks(1:selftasksize, 1:nselftasks+1))
            CALL process_selftasks(selftasks, nselftasks)
            !$omp target exit data map(delete: &
            !$omp&  selftasks(1:selftasksize, 1:nselftasks+1))

            ! At his point, one execution has been performed.

            ! Allocating persistent workpackage with accurate dimensions
            ALLOCATE(wptr%sendtasks(sendtasksize, nsendtasks+1))
            ALLOCATE(wptr%recvtasks(recvtasksize, nrecvtasks+1))
            ALLOCATE(wptr%selftasks(selftasksize, nselftasks+1))

            ALLOCATE(wptr%mpisendtasks(mpitasksize, nmpisendtasks+1))
            ALLOCATE(wptr%mpirecvtasks(mpitasksize, nmpirecvtasks+1))

            ! Tranfering the recorded tasks to the persistent workpackage
            wptr%sendtasks(:, 1:nsendtasks+1) = sendtasks(:, 1:nsendtasks+1)
            wptr%recvtasks(:, 1:nrecvtasks+1) = recvtasks(:, 1:nrecvtasks+1)
            wptr%selftasks(:, 1:nselftasks+1) = selftasks(:, 1:nselftasks+1)

            wptr%mpisendtasks(:, 1:nmpisendtasks+1) = &
                mpisendtasks(:, 1:nmpisendtasks+1)
            wptr%mpirecvtasks(:, 1:nmpirecvtasks+1) = &
                mpirecvtasks(:, 1:nmpirecvtasks+1)

            !$omp target enter data map(to: &
            !$omp&  wptr%sendtasks(1:sendtasksize, 1:nsendtasks+1), &
            !$omp&  wptr%recvtasks(1:recvtasksize, 1:nrecvtasks+1), &
            !$omp&  wptr%selftasks(1:selftasksize, 1:nselftasks+1))

            wptr%is_init = .TRUE.

            ! Deallocating the temporary arrays (already gone on device)
            DEALLOCATE(sendtasks, recvtasks, selftasks, mpisendtasks, &
                mpirecvtasks)

        END IF

        END ASSOCIATE

    END SUBROUTINE ctof2


    SUBROUTINE prepare_mpirecvtasks(mpirtasks, nmpirtasks, ilevel)

        ! Subroutine arguments
        INTEGER(intk), INTENT(inout) :: mpirtasks(mpitasksize, maxsize)
        INTEGER(intk), INTENT(out) :: nmpirtasks
        INTEGER(intk), INTENT(in) :: ilevel   ! Level of the *fine* side

        ! Local variables
        INTEGER(intk) :: i, igridf, iprocnbr
        INTEGER(intk) :: kk, jj, ii, impirtasks, ncells
        INTEGER(int32) :: recvcounter, messagelength

        ! Post all receive calls
        recvcounter = 0
        messagelength = 0
        nrecv = 0
        recvidxlist = -1
        recvlist = 0

        ! Initialize the internal counter
        impirtasks = 0

        ! Iteration over all receive connections
        DO i = 1, irecv
            igridf = recvconns(3, i)
            IF (ilevel == level(igridf)) THEN
                iprocnbr = recvconns(2, i) ! The sender process (coarse side)

                IF (myid == iprocnbr) THEN
                    ! Purely local connection, no MPI communication needed
                    CYCLE
                END IF

                CALL get_mgdims(kk, jj, ii, igridf)
                ncells = kk*jj*ii/8

                ! Entering the information needed for unpacking
                recvidxlist(1, i) = iprocnbr
                recvidxlist(2, i) = ncells
                recvidxlist(3, i) = recvcounter + messagelength
                messagelength = messagelength + ncells

                IF (recvcounter + messagelength > idim_mg_bufs) THEN
                    CALL errr(__FILE__, __LINE__)
                END IF
            END IF

            IF (messagelength > 0) THEN
                IF (i == irecv) THEN
                    ! Adding a new MPI receive task (replaces "post_recv")
                    impirtasks = impirtasks + 1
                    CALL add_mpi_task(mpirtasks(:, impirtasks), iprocnbr, &
                        messagelength, recvcounter)

                ELSE IF (recvconns(2, i + 1) /= iprocnbr) THEN
                    ! Adding a new MPI receive task (replaces "post_recv")
                    impirtasks = impirtasks + 1
                    CALL add_mpi_task(mpirtasks(:, impirtasks), iprocnbr, &
                        messagelength, recvcounter)

                END IF
            END IF
        END DO

        ! Assign the number of tasks to the output variable
        nmpirtasks = impirtasks

        ! Adding a dummy entry at position (end+1)
        mpirtasks(:, nmpirtasks+1) = -1

    END SUBROUTINE prepare_mpirecvtasks



    SUBROUTINE add_mpi_task(mpitasks, iprocnbr, messagelength, recvcounter)
        ! Subroutine arguments
        INTEGER(intk), INTENT(inout) :: mpitasks(mpitasksize)

        INTEGER(int32), INTENT(in) :: iprocnbr
        INTEGER(int32), INTENT(inout) :: messagelength
        INTEGER(int32), INTENT(inout) :: recvcounter

        ! Local variables
        ! none...

        nrecv = nrecv + 1
        recvlist(nrecv) = iprocnbr

        mpitasks(1) = recvcounter+1        ! = index for recvbuf
        mpitasks(2) = messagelength        ! = length of message
        mpitasks(3) = iprocnbr             ! = sender process

        recvcounter = recvcounter + messagelength
        messagelength = 0

    END SUBROUTINE add_mpi_task



    SUBROUTINE process_mpirecvtasks(mpirtasks, nmpirtasks)
        ! Subroutine arguments
        INTEGER(intk), INTENT(in) :: nmpirtasks
        INTEGER(intk), INTENT(in) :: mpirtasks(mpitasksize, nmpirtasks+1)

        ! Local variables
        INTEGER(intk) :: i, offset, messagelength, iprocc

#ifdef _MGLET_PROFILE_ANNOTATIONS_
        CALL profile_range_push("process_mpirecv")
#endif

        !$omp target data use_device_addr(recvbuf)
        DO i = 1, nmpirtasks

            ! Getting connection information
            offset = mpirtasks(1, i)
            messagelength = mpirtasks(2, i)
            iprocc = mpirtasks(3, i)

            IF (iprocc == myid) CALL errr(__FILE__, __LINE__)

            CALL MPI_Irecv(recvbuf(offset), messagelength, mglet_mpi_real, &
                iprocc, 1, MPI_COMM_WORLD, recvreqs(i))
        END DO
        !$omp end target data

#ifdef _MGLET_PROFILE_ANNOTATIONS_
        CALL profile_range_pop()
#endif

        ! Checking for the dummy entry at position (end+1)
        IF (mpirtasks(1, nmpirtasks+1) /= -1) THEN
            CALL errr(__FILE__, __LINE__)
        END IF

        ! Setting the number of posted receives
        nrecv = nmpirtasks

    END SUBROUTINE process_mpirecvtasks


    ! Perform all Send-calls
    SUBROUTINE prepare_allsendtasks(stasks, nstasks, etasks, netasks, &
        mpistasks, nmpistasks, ilevel)

        ! Subroutine arguments
        INTEGER(intk), INTENT(inout) :: stasks(sendtasksize, maxsize)
        INTEGER(intk), INTENT(out) :: nstasks

        INTEGER(intk), INTENT(inout) :: etasks(selftasksize, maxsize)
        INTEGER(intk), INTENT(out) :: netasks

        INTEGER(intk), INTENT(inout) :: mpistasks(mpitasksize, maxsize)
        INTEGER(intk), INTENT(out) :: nmpistasks

        INTEGER(intk), INTENT(in) :: ilevel   ! Level of the *fine* side

        ! Local variables
        INTEGER(intk) :: i, igridf, igridc, iprocnbr
        INTEGER(intk) :: impistasks, istasks, ietasks
        INTEGER(int32) :: sendcounter, messagelength

        ! Post all receive calls
        sendcounter = 0
        messagelength = 0
        nsend = 0

        ! Initialize the internal counter
        istasks = 0
        ietasks = 0
        impistasks = 0

        ! Iteration over all send connections
        DO i = 1, isend

            iprocnbr = sendconns(1, i)
            igridf = sendconns(3, i)
            igridc = sendconns(4, i)

            IF (ilevel == level(igridf)) THEN

                IF (myid == iprocnbr) THEN

                    ! Purely local connection, no MPI communication needed
                    ietasks = ietasks + 1
                    CALL add_selftask(etasks(:, ietasks), igridf, igridc)

                    ! Cycling to avoid adding a new MPI send task
                    CYCLE

                ELSE

                    ! Connection to a different process, MPI communication
                    istasks = istasks + 1
                    CALL add_sendtask(stasks(:, istasks), i, messagelength, &
                        sendcounter)

                END IF
            END IF

            IF (messagelength > 0) THEN
                IF (i == isend) THEN
                    ! Adding a new MPI send task (replaces "post_send")
                    impistasks = impistasks + 1
                    CALL add_mpi_task(mpistasks(:, impistasks), iprocnbr, &
                        messagelength, sendcounter)

                ELSE IF (sendconns(1, i + 1) /= iprocnbr) THEN
                    ! Adding a new MPI send task (replaces "post_send")
                    impistasks = impistasks + 1
                    CALL add_mpi_task(mpistasks(:, impistasks), iprocnbr, &
                        messagelength, sendcounter)

                END IF
            END IF

        END DO

        ! Assign the number of tasks to the output variable
        nmpistasks = impistasks
        nstasks = istasks
        netasks = ietasks

        ! Adding a dummy entry at position (end+1)
        mpistasks(:, nmpistasks+1) = -1
        stasks(:, nstasks+1) = -1
        etasks(:, netasks+1) = -1

    END SUBROUTINE prepare_allsendtasks


    SUBROUTINE add_selftask(etask, igridf, igridc)

        ! Subroutine arguments
        INTEGER(intk), INTENT(inout) :: etask(selftasksize)
        INTEGER(int32), INTENT(in) :: igridf, igridc

        ! Local variables
        INTEGER(intk) :: kkf, jjf, iif
        INTEGER(intk) :: ista, jsta, ksta, isto, jsto, ksto

        ! Compute start- and end-positions in coarse grid
        ista = iposition(igridf) - 1
        jsta = jposition(igridf) - 1
        ksta = kposition(igridf) - 1

        CALL get_mgdims(kkf, jjf, iif, igridf)
        isto = ista + (iif - 4)/2 + 1
        jsto = jsta + (jjf - 4)/2 + 1
        ksto = ksta + (kkf - 4)/2 + 1

        ! Filling the task array
        etask(1) = igridf
        etask(2) = igridc
        etask(3) = ista
        etask(4) = jsta
        etask(5) = ksta
        etask(6) = isto
        etask(7) = jsto
        etask(8) = ksto

    END SUBROUTINE add_selftask


    SUBROUTINE process_selftasks(selftasks, nselftasks)
        ! Subroutine arguments
        INTEGER(intk), INTENT(in) :: nselftasks
        INTEGER(intk), INTENT(in) :: selftasks(selftasksize, nselftasks+1)

        ! Local variables
        INTEGER(intk) :: itask, ista, jsta, ksta, isto, jsto, ksto
        INTEGER(intk) :: iif, jjf, kkf, iic, jjc, kkc
        INTEGER(intk) :: igridc, igridf, ip3c, ip3f

#ifdef _MGLET_PROFILE_ANNOTATIONS_
        CALL profile_range_push("process_selftasks")
#endif

        ASSOCIATE(fc => fc%arr, ff => ff%arr)

        !$omp target teams distribute private(itask, igridf, igridc, &
        !$omp&  ista, jsta, ksta, isto, jsto, ksto, ip3f, ip3c, &
        !$omp&  kkf, jjf, iif, kkc, jjc, iic)
        DO itask = 1, nselftasks

            ! Unpacking the task
            igridf = selftasks(1, itask)
            igridc = selftasks(2, itask)
            ista = selftasks(3, itask)
            jsta = selftasks(4, itask)
            ksta = selftasks(5, itask)
            isto = selftasks(6, itask)
            jsto = selftasks(7, itask)
            ksto = selftasks(8, itask)

            ! Getting the parameters of the grids
            CALL get_ip3(ip3f, igridf)
            CALL get_ip3(ip3c, igridc)
            CALL get_mgdims(kkf, jjf, iif, igridf)
            CALL get_mgdims(kkc, jjc, iic, igridc)

            CALL copy_kernel(kkf, jjf, iif, kkc, jjc, iic, &
                ff(ip3f), fc(ip3c), ista, jsta, ksta, isto, jsto, ksto)
        END DO
        !$omp end target teams distribute

        END ASSOCIATE

#ifdef _MGLET_PROFILE_ANNOTATIONS_
        CALL profile_range_pop()
#endif

        ! Checking for the dummy entry at position (end+1)
        IF (selftasks(1, nselftasks+1) /= -1) THEN
            CALL errr(__FILE__, __LINE__)
        END IF

    END SUBROUTINE process_selftasks



    SUBROUTINE copy_kernel(kkf, jjf, iif, kkc, jjc, iic, &
        ff, fc, ista, jsta, ksta, isto, jsto, ksto)
        !$omp declare target

        ! Subroutine arguments
        INTEGER(intk), INTENT(in) :: kkf, jjf, iif
        INTEGER(intk), INTENT(in) :: kkc, jjc, iic
        REAL(realk), INTENT(inout) :: ff(kkf, jjf, iif)
        REAL(realk), INTENT(in) :: fc(kkc, jjc, iic)
        INTEGER(intk), INTENT(in) :: ista, jsta, ksta
        INTEGER(intk), INTENT(in) :: isto, jsto, ksto

        ! Local variables
        INTEGER(intk) :: i, j, k, ic, jc, kc

        !$omp parallel do collapse(3) private(i, j, k, ic, jc, kc)
        DO i = 1, iif
            DO j = 1, jjf
                DO k = 1, kkf

                    ! Computing the corresponding coarse grid indices
                    ic = (i-1)/2 + ista
                    jc = (j-1)/2 + jsta
                    kc = (k-1)/2 + ksta

                    ! Checking that the coarse grid indices are within bounds
                    IF (ic > isto .OR. jc > jsto .OR. kc > ksto) THEN
                        CALL errr(__FILE__, __LINE__)
                    END IF

                    ! Copying the value from the coarse grid to the fine grid
                    ff(k, j, i) = fc(kc, jc, ic)

                END DO
            END DO
        END DO
        !$omp end parallel do

    END SUBROUTINE copy_kernel


    SUBROUTINE process_mpisendtasks(mpistasks, nmpistasks)
        ! Subroutine arguments
        INTEGER(intk), INTENT(in) :: nmpistasks
        INTEGER(intk), INTENT(in) :: mpistasks(mpitasksize, nmpistasks+1)

        ! Local variables
        INTEGER(intk) :: i, idx_sendbuf, messagelength, iprocf

#ifdef _MGLET_PROFILE_ANNOTATIONS_
        CALL profile_range_push("process_mpisend")
#endif

        !$omp target data use_device_addr(sendbuf)
        DO i = 1, nmpistasks

            ! Getting connection information
            idx_sendbuf = mpistasks(1, i)
            messagelength = mpistasks(2, i)
            iprocf = mpistasks(3, i)

            IF (iprocf == myid) CALL errr(__FILE__, __LINE__)

            ! Non-blocking MPI call with request handle stored in sendreqs
            CALL MPI_Isend(sendbuf(idx_sendbuf), messagelength, &
                mglet_mpi_real, iprocf, 1, MPI_COMM_WORLD, sendreqs(i))
        END DO
        !$omp end target data

        ! Checking for the dummy entry at position (end+1)
        IF (mpistasks(1, nmpistasks+1) /= -1) THEN
            CALL errr(__FILE__, __LINE__)
        END IF

        ! Setting the number of posted receives
        nsend = nmpistasks

#ifdef _MGLET_PROFILE_ANNOTATIONS_
        CALL profile_range_pop()
#endif

    END SUBROUTINE process_mpisendtasks



    SUBROUTINE process_sendtasks(stasks, nstasks)
        ! Subroutine arguments
        INTEGER(intk), INTENT(in) :: nstasks
        INTEGER(intk), INTENT(in) :: stasks(sendtasksize, nstasks+1)

        ! Local variables
        INTEGER(intk) :: igridc, itask, ii, jj, kk, ip3
        INTEGER(intk) :: ista, jsta, ksta, isto, jsto, ksto, idx1, idx2

        IF (nstasks == 0) RETURN

#ifdef _MGLET_PROFILE_ANNOTATIONS_
        CALL profile_range_push("process_sendtasks")
#endif

        ASSOCIATE(coarse => fc%arr)

        !$omp target teams distribute private(itask, ii, jj, kk, ip3, &
        !$omp&  ista, jsta, ksta, isto, jsto, ksto, idx1, idx2, igridc)
        DO itask = 1, nstasks

            ! Unpacking the task
            ista = stasks(1, itask)
            jsta = stasks(2, itask)
            ksta = stasks(3, itask)
            isto = stasks(4, itask)
            jsto = stasks(5, itask)
            ksto = stasks(6, itask)
            igridc = stasks(7, itask)
            idx1 = stasks(8, itask)
            idx2 = stasks(9, itask)

            ! Getting the parameters of the coarse grid
            CALL get_ip3(ip3, igridc)
            CALL get_mgdims(kk, jj, ii, igridc)

            CALL write_buffer(kk, jj, ii, sendbuf(idx1:idx2), &
                coarse(ip3), ista, jsta, ksta, isto, jsto, ksto)

        END DO
        !$omp end target teams distribute

        END ASSOCIATE

#ifdef _MGLET_PROFILE_ANNOTATIONS_
        CALL profile_range_pop()
#endif

    END SUBROUTINE process_sendtasks


    SUBROUTINE write_buffer(kk, jj, ii, buf, fc, ista, jsta, ksta, &
            isto, jsto, ksto)
        !$omp declare target

        ! Subroutine arguments
        INTEGER(intk), INTENT(in) :: kk, jj, ii
        REAL(realk), INTENT(inout) :: buf(:)
        REAL(realk), INTENT(in) :: fc(kk, jj, ii)
        INTEGER(intk), INTENT(in) :: ista, jsta, ksta
        INTEGER(intk), INTENT(in) :: isto, jsto, ksto

        ! Local variables
        INTEGER(intk) :: i, j, k
        INTEGER(intk) :: idx, iic, jjc, kkc

        ! Compute extent of the patch used on the coarse grid
        iic = isto - ista + 1
        jjc = jsto - jsta + 1
        kkc = ksto - ksta + 1

        !$omp parallel do collapse(3) private(i, j, k, idx)
        DO i = ista, isto
            DO j = jsta, jsto
                DO k = ksta, ksto
                    idx = 1 + (k - ksta) + (j - jsta)*kkc + (i - ista)*kkc*jjc
                    buf(idx) = fc(k, j, i)
                END DO
            END DO
        END DO
        !$omp end parallel do

    END SUBROUTINE write_buffer


    SUBROUTINE add_sendtask(stask, sendid, messagelength, sendcounter)
        ! Subroutine arguments
        INTEGER(intk), INTENT(inout) :: stask(sendtasksize)
        INTEGER(int32), INTENT(in) :: sendid
        INTEGER(int32), INTENT(inout) :: messagelength
        INTEGER(int32), INTENT(in) :: sendcounter

        ! Local variables
        INTEGER(intk) :: kkf, jjf, iif
        INTEGER(intk) :: ista, jsta, ksta, isto, jsto, ksto
        INTEGER(intk) :: igridc, igridf
        INTEGER(int32) :: thismessagelength, offset

        igridf = sendconns(3, sendid)
        igridc = sendconns(4, sendid)

        ! Compute start- and end-positions in coarse grid
        ista = iposition(igridf) - 1
        jsta = jposition(igridf) - 1
        ksta = kposition(igridf) - 1

        CALL get_mgdims(kkf, jjf, iif, igridf)
        isto = ista + (iif - 4)/2 + 1
        jsto = jsta + (jjf - 4)/2 + 1
        ksto = ksta + (kkf - 4)/2 + 1

        thismessagelength = (isto-ista+1)*(jsto-jsta+1)*(ksto-ksta+1)

        IF (sendcounter + messagelength + thismessagelength > idim_mg_bufs) THEN
            CALL errr(__FILE__, __LINE__)
        END IF

        offset = sendcounter + messagelength + 1

        ! Writing the task
        stask(1) = ista
        stask(2) = jsta
        stask(3) = ksta
        stask(4) = isto
        stask(5) = jsto
        stask(6) = ksto
        stask(7) = igridc
        stask(8) = offset
        stask(9) = offset + thismessagelength - 1

        messagelength = messagelength + thismessagelength

    END SUBROUTINE add_sendtask


    SUBROUTINE prepare_recvtasks(rtasks, nrtasks)
        ! Subroutine arguments
        INTEGER(intk), INTENT(inout) :: rtasks(:, :)
        INTEGER(intk), INTENT(out) :: nrtasks

        ! Local variables
        TYPE(MPI_Status) :: recvstatus
        INTEGER(intk) :: i, igridf, offset
        INTEGER(int32) :: idx, recvmessagelen, unpacklen, irtasks

        ! Initialize the internal counter
        irtasks = 0

        DO WHILE (.TRUE.)

            IF (nrecv == 0) EXIT
            CALL MPI_Waitany(nrecv, recvreqs, idx, recvstatus)

            IF (idx /= MPI_UNDEFINED) THEN
                ! The index "idx" is returned and can be used to identify the
                ! connection as recvlist(nrecv) = iprocnbr

                unpacklen = 0

                ! Check which non-zero packages belong to this sender and unpack
                DO i = 1, irecv
                    IF (recvidxlist(1, i) == recvlist(idx) &
                            .AND. recvidxlist(2, i) > 0) THEN

                        ! Looking up information for recvid = i
                        igridf = recvconns(3, i)
                        offset = recvidxlist(3, i) + 1

                        ! Adding new recevive task
                        ! (= storing fine grid and offset in recv buffer)
                        irtasks = irtasks + 1
                        rtasks(1, irtasks) = igridf
                        rtasks(2, irtasks) = offset

                        unpacklen = unpacklen + recvidxlist(2, i)
                    END IF
                END DO

                ! Ensure that the complete message has been unpacked
                CALL MPI_Get_count(recvstatus, mglet_mpi_real, recvmessagelen)
                IF (recvmessagelen /= unpacklen) THEN
                    CALL errr(__FILE__, __LINE__)
                END IF

            ELSE
                ! If idx = MPI_UNDEFINED, then all requests are complete
                EXIT
            END IF
        END DO

        ! Make sure that also all non-blocking sends have completed
        CALL MPI_Waitall(nsend, sendreqs, MPI_STATUSES_IGNORE)

        ! Assigning the total count
        nrtasks = irtasks

        ! Adding a dummy entry at position (end+1)
        rtasks(:, nrtasks+1) = -1

    END SUBROUTINE prepare_recvtasks



    SUBROUTINE process_recvtasks(rtasks, nrtasks)
        ! Subroutine arguments
        INTEGER(intk), INTENT(in) :: nrtasks
        INTEGER(intk), INTENT(in) :: rtasks(recvtasksize, nrtasks+1)

        ! Local variables
        INTEGER(intk) :: igridf, idx, len, itask, ii, jj, kk, ip3

        IF (nrtasks == 0) RETURN

#ifdef _MGLET_PROFILE_ANNOTATIONS_
        CALL profile_range_push("process_recvtasks")
#endif

        ASSOCIATE(fine => ff%arr)

        !$omp target teams distribute private(itask, ii, jj, kk, ip3, &
        !$omp&  igridf, idx, len)
        DO itask = 1, nrtasks

            ! Unpacking the task
            igridf = rtasks(1, itask)
            idx = rtasks(2, itask)

            ! Getting parameters of fine grid
            CALL get_mgdims(kk, jj, ii, igridf)
            CALL get_ip3(ip3, igridf)
            len = kk * jj * ii / 8

            ! Copying from buffer to the fine grid
            CALL write_fine(kk, jj, ii, fine(ip3), len, recvbuf(idx:idx+len-1))

        END DO
        !$omp end target teams distribute

        END ASSOCIATE

#ifdef _MGLET_PROFILE_ANNOTATIONS_
        CALL profile_range_pop()
#endif

    END SUBROUTINE process_recvtasks


    SUBROUTINE write_fine(kk, jj, ii, fff, len, fc)
        !$omp declare target

        ! Subroutine arguments
        INTEGER(intk), INTENT(in) :: kk, jj, ii
        REAL(realk), INTENT(inout) :: fff(kk, jj, ii)
        INTEGER(intk), INTENT(in) :: len
        REAL(realk), INTENT(in) :: fc(len)

        ! Local variables
        INTEGER(intk) :: k, j, i, kc, jc, ic, idx
        INTEGER(intk) :: kkc, jjc, iic

        ! The counter for the buffer reflects the coarse frid
        kkc = kk/2
        jjc = jj/2
        iic = ii/2

        !$omp parallel do collapse(3) private(i, j, k, idx, kc, jc, ic)
        DO i = 1, ii
            DO j = 1, jj
                DO k = 1, kk
                    ! Reconstructing the coarse grid indices
                    ic = (i-1)/2 + 1
                    jc = (j-1)/2 + 1
                    kc = (k-1)/2 + 1
                    ! Conversion into buffer index
                    idx = 1 + (kc-1) + (jc-1)*kkc + (ic-1)*kkc*jjc
                    fff(k, j, i) = fc(idx)
                END DO
            END DO
        END DO
        !$omp end parallel do

    END SUBROUTINE write_fine




    ! Initialize arrays and data types
    SUBROUTINE init_ctof2()

        ! Local variables
        INTEGER(intk) :: ilevel
        TYPE(field_t) :: dummy

        ! Check if ctof_core_mod necessary provides infrastructure
        IF (.NOT. has_infrastructure) THEN
            CALL errr(__FILE__, __LINE__)
        END IF

        ! Local variables
        ALLOCATE(sendreqs(isend))
        ALLOCATE(recvreqs(irecv))
        ALLOCATE(recvlist(irecv))
        ALLOCATE(recvidxlist(3, irecv))

        nrecv = 0
        nsend = 0
        is_init = .TRUE.

        ! Allocate the workrecords array for all possible ilevels
        ALLOCATE(workrecords(minlevel:maxlevel))

        CALL dummy%init("DUMMY")

        ! Recording operations for all levels
        is_recording = .TRUE.
        DO ilevel = minlevel, maxlevel
            CALL ctof2(ilevel, dummy, dummy)
        END DO
        is_recording = .FALSE.

        CALL dummy%finish()

    END SUBROUTINE init_ctof2



    SUBROUTINE finish_ctof2()
        INTEGER(intk) :: ilevel

        IF (.NOT. is_init) RETURN
        is_init = .FALSE.

        DEALLOCATE(sendreqs)
        DEALLOCATE(recvreqs)
        DEALLOCATE(recvlist)
        DEALLOCATE(recvidxlist)

        ! Deallocate the workpackage components for each level
        DO ilevel = minlevel+1, maxlevel
            IF (workrecords(ilevel)%is_init) THEN

                !$omp target exit data map(delete: &
                !$omp&  workrecords(ilevel)%sendtasks(:, :), &
                !$omp&  workrecords(ilevel)%recvtasks(:, :))

                DEALLOCATE(workrecords(ilevel)%sendtasks)
                DEALLOCATE(workrecords(ilevel)%recvtasks)
                DEALLOCATE(workrecords(ilevel)%mpisendtasks)
                DEALLOCATE(workrecords(ilevel)%mpirecvtasks)
            END IF
        END DO
        DEALLOCATE(workrecords)

    END SUBROUTINE finish_ctof2

END MODULE ctof2_mod
