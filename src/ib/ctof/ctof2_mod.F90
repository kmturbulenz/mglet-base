MODULE ctof2_mod
    USE core_mod
    USE ctof_core_mod
    USE MPI_f08

    IMPLICIT NONE (type, external)
    PRIVATE

    TYPE :: sendtask_t
        INTEGER(intk) :: ista
        INTEGER(intk) :: jsta
        INTEGER(intk) :: ksta
        INTEGER(intk) :: isto
        INTEGER(intk) :: jsto
        INTEGER(intk) :: ksto
        INTEGER(intk) :: igridc
        INTEGER(intk) :: idx1
        INTEGER(intk) :: idx2
    END TYPE sendtask_t

    TYPE :: recvtask_t
        INTEGER(intk) :: igridf
        INTEGER(intk) :: offset
    END TYPE recvtask_t

    TYPE :: selftask_t
        INTEGER(intk) :: igridf
        INTEGER(intk) :: igridc
        INTEGER(intk) :: ista
        INTEGER(intk) :: jsta
        INTEGER(intk) :: ksta
        INTEGER(intk) :: isto
        INTEGER(intk) :: jsto
        INTEGER(intk) :: ksto
    END TYPE selftask_t

    TYPE :: mpitask_t
        INTEGER(intk) :: offset
        INTEGER(intk) :: messagelength
        INTEGER(intk) :: iprocnbr
    END TYPE mpitask_t

    ! Lists that hold the send and receive request arrays
    TYPE(MPI_Request), ALLOCATABLE :: sendreqs(:), recvreqs(:)
    INTEGER(int32), ALLOCATABLE :: recvlist(:)
    INTEGER(intk), ALLOCATABLE :: recvidxlist(:, :)
    INTEGER(intk) :: nsend, nrecv

    ! Workpackages containing individual tasks for packing / unpacking
    TYPE(sendtask_t), ALLOCATABLE :: sendtasks(:)
    TYPE(recvtask_t), ALLOCATABLE :: recvtasks(:)
    TYPE(selftask_t), ALLOCATABLE :: selftasks(:)
    TYPE(mpitask_t), ALLOCATABLE :: mpisendtasks(:), mpirecvtasks(:)
    ! Not declare target since there is no jit ctof

    ! Type to hold condensed task arrays to execute a certain type of ctof
    TYPE :: work_t
        LOGICAL :: is_init = .FALSE.
        TYPE(sendtask_t), ALLOCATABLE :: sendtasks(:)
        TYPE(recvtask_t), ALLOCATABLE :: recvtasks(:)
        TYPE(selftask_t), ALLOCATABLE :: selftasks(:)
        TYPE(mpitask_t), ALLOCATABLE :: mpisendtasks(:)
        TYPE(mpitask_t), ALLOCATABLE :: mpirecvtasks(:)
    END TYPE work_t

    ! Array to store instructions for different values of "ilevel"
    TYPE(work_t), ALLOCATABLE, TARGET :: workrecords(:)
    LOGICAL :: is_recording = .FALSE.

    ! Maximum size of temporary arrays for recording
    INTEGER(intk) :: maxsize

    ! Variable to indicate if the required data structures have been created
    LOGICAL :: is_init = .FALSE.

    PUBLIC :: ctof2, init_ctof2, finish_ctof2
CONTAINS

    SUBROUTINE ctof2(ilevel, ff, fc)
        ! Subroutine arguments
        INTEGER(intk), INTENT(in) :: ilevel  ! Level of the *fine* side
        TYPE(field_t), TARGET, INTENT(inout) :: ff
        TYPE(field_t), TARGET, INTENT(in) :: fc

        IF (.NOT. is_init) CALL errr(__FILE__, __LINE__)

        ! Looking up the workpackage for this level
        ASSOCIATE(wptr => workrecords(ilevel))

        IF (is_recording) THEN
            CALL recording_pass(wptr, ilevel, ff, fc)
        ELSE
            CALL recorded_ctof(wptr, ff, fc)
        END IF

        END ASSOCIATE
    END SUBROUTINE ctof2


    SUBROUTINE recorded_ctof(wptr, ff, fc)
        ! Subroutine arguments
        TYPE(work_t), INTENT(in) :: wptr
        TYPE(field_t), TARGET, INTENT(inout) :: ff
        TYPE(field_t), TARGET, INTENT(in) :: fc

        ! Local variables
        INTEGER(intk) :: nsendtasks, nrecvtasks
        INTEGER(intk) :: nselftasks
        INTEGER(intk) :: nmpisendtasks, nmpirecvtasks

        ! Obtaining numbers of tasks from the workpackage
        nmpirecvtasks = SIZE(wptr%mpirecvtasks) - 1
        nmpisendtasks = SIZE(wptr%mpisendtasks) - 1
        nsendtasks = SIZE(wptr%sendtasks) - 1
        nrecvtasks = SIZE(wptr%recvtasks) - 1
        nselftasks = SIZE(wptr%selftasks) - 1

        CALL process_mpirecvtasks(wptr%mpirecvtasks, nmpirecvtasks)
        CALL process_sendtasks(fc, nsendtasks, wptr%sendtasks)
        CALL process_mpisendtasks(wptr%mpisendtasks, nmpisendtasks)

        ! While MPI messages travel, local work is done
        CALL process_selftasks(fc, ff, nselftasks, wptr%selftasks)

        CALL MPI_Waitall(nrecv, recvreqs, MPI_STATUSES_IGNORE)
        CALL process_recvtasks(ff, nrecvtasks, wptr%recvtasks)
        CALL MPI_Waitall(nsend, sendreqs, MPI_STATUSES_IGNORE)
    END SUBROUTINE recorded_ctof


    SUBROUTINE recording_pass(wptr, ilevel, ff, fc)
        ! Subroutine arguments
        TYPE(work_t), INTENT(inout) :: wptr
        INTEGER(intk), INTENT(in) :: ilevel
        TYPE(field_t), TARGET, INTENT(inout) :: ff
        TYPE(field_t), TARGET, INTENT(in) :: fc

        ! Local variables
        INTEGER(intk) :: nsendtasks, nrecvtasks
        INTEGER(intk) :: nselftasks
        INTEGER(intk) :: nmpisendtasks, nmpirecvtasks

        ! During the recording phase, the workpackage is stored
        IF (wptr%is_init) THEN
            WRITE(*, *) "Repeated recording for ilevel = ", ilevel
            CALL errr(__FILE__, __LINE__)
        END IF

        ! Allocating oversized temporary arrays to record the tasks
        maxsize = noflevel(ilevel) + 1
        ALLOCATE(sendtasks(maxsize))
        ALLOCATE(recvtasks(maxsize))
        ALLOCATE(selftasks(maxsize))
        ALLOCATE(mpisendtasks(maxsize))
        ALLOCATE(mpirecvtasks(maxsize))

        CALL prepare_allsendtasks(&
            sendtasks, nsendtasks, &
            selftasks, nselftasks, &
            mpisendtasks, nmpisendtasks, ilevel)

        CALL prepare_mpirecvtasks(mpirecvtasks, nmpirecvtasks, ilevel)

        !$omp target enter data map(to: &
        !$omp&  sendtasks(1:nsendtasks+1))
        CALL process_sendtasks(fc, nsendtasks, sendtasks)
        !$omp target exit data map(delete: &
        !$omp&  sendtasks(1:nsendtasks+1))

        CALL process_mpirecvtasks(mpirecvtasks, nmpirecvtasks)
        CALL process_mpisendtasks(mpisendtasks, nmpisendtasks)

        ! Includes waiting for MPI communication to finish
        CALL prepare_recvtasks(nrecvtasks, recvtasks)

        !$omp target enter data map(to: &
        !$omp&  recvtasks(1:nrecvtasks+1))
        CALL process_recvtasks(ff, nrecvtasks, recvtasks)
        !$omp target exit data map(delete: &
        !$omp&  recvtasks(1:nrecvtasks+1))

        !$omp target enter data map(to: &
        !$omp&  selftasks(1:nselftasks+1))
        CALL process_selftasks(fc, ff, nselftasks, selftasks)
        !$omp target exit data map(delete: &
        !$omp&  selftasks(1:nselftasks+1))

        ! At this point, one execution has been performed.

        ! Allocating persistent workpackage with accurate dimensions
        ALLOCATE(wptr%sendtasks(nsendtasks+1))
        ALLOCATE(wptr%recvtasks(nrecvtasks+1))
        ALLOCATE(wptr%selftasks(nselftasks+1))

        ALLOCATE(wptr%mpisendtasks(nmpisendtasks+1))
        ALLOCATE(wptr%mpirecvtasks(nmpirecvtasks+1))

        ! Transferring the recorded tasks to the persistent workpackage
        wptr%sendtasks(1:nsendtasks+1) = sendtasks(1:nsendtasks+1)
        wptr%recvtasks(1:nrecvtasks+1) = recvtasks(1:nrecvtasks+1)
        wptr%selftasks(1:nselftasks+1) = selftasks(1:nselftasks+1)

        wptr%mpisendtasks(1:nmpisendtasks+1) = &
            mpisendtasks(1:nmpisendtasks+1)
        wptr%mpirecvtasks(1:nmpirecvtasks+1) = &
            mpirecvtasks(1:nmpirecvtasks+1)

        !$omp target enter data map(to: &
        !$omp&  wptr%sendtasks(1:nsendtasks+1), &
        !$omp&  wptr%recvtasks(1:nrecvtasks+1), &
        !$omp&  wptr%selftasks(1:nselftasks+1))

        wptr%is_init = .TRUE.

        ! Deallocating the temporary arrays (already gone on device)
        DEALLOCATE(sendtasks, recvtasks, selftasks, mpisendtasks, &
            mpirecvtasks)
    END SUBROUTINE recording_pass


    SUBROUTINE prepare_mpirecvtasks(mpirtasks, nmpirtasks, ilevel)
        ! Subroutine arguments
        TYPE(mpitask_t), INTENT(inout) :: mpirtasks(maxsize)
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
                    CALL add_mpi_task(mpirtasks(impirtasks), iprocnbr, &
                        messagelength, recvcounter, "recv")

                ELSE IF (recvconns(2, i + 1) /= iprocnbr) THEN
                    ! Adding a new MPI receive task (replaces "post_recv")
                    impirtasks = impirtasks + 1
                    CALL add_mpi_task(mpirtasks(impirtasks), iprocnbr, &
                        messagelength, recvcounter, "recv")

                END IF
            END IF
        END DO

        ! Assign the number of tasks to the output variable
        nmpirtasks = impirtasks

        ! Adding a dummy entry at position (end+1)
        mpirtasks(nmpirtasks+1) = mpitask_t(-1, -1, -1)
    END SUBROUTINE prepare_mpirecvtasks


    SUBROUTINE add_mpi_task(mpitasks, iprocnbr, messagelength, counter, type)
        ! Subroutine arguments
        TYPE(mpitask_t), INTENT(inout) :: mpitasks

        INTEGER(int32), INTENT(in) :: iprocnbr
        INTEGER(int32), INTENT(inout) :: messagelength
        INTEGER(int32), INTENT(inout) :: counter
        CHARACTER(len=4), INTENT(in) :: type

        ! Local variables
        ! none...

        IF (type == "recv") THEN
            nrecv = nrecv + 1
            recvlist(nrecv) = iprocnbr
        ELSE IF (type /= "send") THEN
            CALL errr(__FILE__, __LINE__)
        END IF

        mpitasks%offset = counter+1
        mpitasks%messagelength = messagelength
        mpitasks%iprocnbr = iprocnbr

        counter = counter + messagelength
        messagelength = 0
    END SUBROUTINE add_mpi_task


    SUBROUTINE process_mpirecvtasks(mpirtasks, nmpirtasks)
        ! Subroutine arguments
        INTEGER(intk), INTENT(in) :: nmpirtasks
        TYPE(mpitask_t), INTENT(in) :: mpirtasks(nmpirtasks+1)

        ! Local variables
        INTEGER(intk) :: i, offset, messagelength, iprocc

        ! Leaving immediately if there are no tasks to process
        IF (nmpirtasks < 1) RETURN

        CALL profile_range_push("process_mpirecv")

        !$omp target data use_device_addr(recvbuf)
        DO i = 1, nmpirtasks

            ! Getting connection information
            offset = mpirtasks(i)%offset
            messagelength = mpirtasks(i)%messagelength
            iprocc = mpirtasks(i)%iprocnbr

            IF (iprocc == myid) CALL errr(__FILE__, __LINE__)

            CALL MPI_Irecv(recvbuf(offset), messagelength, mglet_mpi_real, &
                iprocc, 1, MPI_COMM_WORLD, recvreqs(i))
        END DO
        !$omp end target data

        CALL profile_range_pop()

        ! Checking for the dummy entry at position (end+1)
        IF (mpirtasks(nmpirtasks+1)%offset /= -1) THEN
            CALL errr(__FILE__, __LINE__)
        END IF

        ! Setting the number of posted receives
        nrecv = nmpirtasks
    END SUBROUTINE process_mpirecvtasks


    SUBROUTINE prepare_allsendtasks(stasks, nstasks, etasks, netasks, &
        mpistasks, nmpistasks, ilevel)
        ! Subroutine arguments
        TYPE(sendtask_t), INTENT(inout) :: stasks(maxsize)
        INTEGER(intk), INTENT(out) :: nstasks
        TYPE(selftask_t), INTENT(inout) :: etasks(maxsize)
        INTEGER(intk), INTENT(out) :: netasks
        TYPE(mpitask_t), INTENT(inout) :: mpistasks(maxsize)
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
                    CALL add_selftask(etasks(ietasks), igridf, igridc)

                    ! Cycling to avoid adding a new MPI send task
                    CYCLE
                ELSE
                    ! Connection to a different process, MPI communication
                    istasks = istasks + 1
                    CALL add_sendtask(stasks(istasks), i, messagelength, &
                        sendcounter)
                END IF
            END IF

            IF (messagelength > 0) THEN
                IF (i == isend) THEN
                    ! Adding a new MPI send task (replaces "post_send")
                    impistasks = impistasks + 1
                    CALL add_mpi_task(mpistasks(impistasks), iprocnbr, &
                        messagelength, sendcounter, "send")

                ELSE IF (sendconns(1, i + 1) /= iprocnbr) THEN
                    ! Adding a new MPI send task (replaces "post_send")
                    impistasks = impistasks + 1
                    CALL add_mpi_task(mpistasks(impistasks), iprocnbr, &
                        messagelength, sendcounter, "send")
                END IF
            END IF
        END DO

        ! Assign the number of tasks to the output variable
        nmpistasks = impistasks
        nstasks = istasks
        netasks = ietasks

        ! Adding a dummy entry at position (end+1)
        mpistasks(nmpistasks+1) = mpitask_t(-1, -1, -1)
        stasks(nstasks+1) = sendtask_t(-1, -1, -1, -1, -1, -1, -1, -1, -1)
        etasks(netasks+1) = selftask_t(-1, -1, -1, -1, -1, -1, -1, -1)
    END SUBROUTINE prepare_allsendtasks


    SUBROUTINE add_selftask(etask, igridf, igridc)
        ! Subroutine arguments
        TYPE(selftask_t), INTENT(inout) :: etask
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
        etask%igridf = igridf
        etask%igridc = igridc
        etask%ista = ista
        etask%jsta = jsta
        etask%ksta = ksta
        etask%isto = isto
        etask%jsto = jsto
        etask%ksto = ksto
    END SUBROUTINE add_selftask


    SUBROUTINE process_selftasks(fc, ff, nselftasks, etasks)
        ! Subroutine arguments
        TYPE(field_t), TARGET, INTENT(in) :: fc
        TYPE(field_t), TARGET, INTENT(inout) :: ff
        INTEGER(intk), INTENT(in) :: nselftasks
        TYPE(selftask_t), INTENT(in) :: etasks(nselftasks+1)

        ! Local variables
        ! none...

        ! Leaving immediately if there are no tasks to process
        IF (nselftasks < 1) RETURN

        CALL profile_range_push("process_selftasks")

        CALL process_selftasks_impl(fc%arr, ff%arr, nselftasks, etasks)

        CALL profile_range_pop()

        ! Checking for the dummy entry at position (end+1)
        IF (etasks(nselftasks+1)%igridf /= -1) THEN
            CALL errr(__FILE__, __LINE__)
        END IF
        END SUBROUTINE process_selftasks


        SUBROUTINE process_selftasks_impl(fc, ff, nselftasks, etasks)
        ! Subroutine arguments
        REAL(realk), INTENT(in) :: fc(*)
        REAL(realk), INTENT(inout) :: ff(*)
        INTEGER(intk), INTENT(in) :: nselftasks
        TYPE(selftask_t), INTENT(in) :: etasks(nselftasks+1)

        ! Local variables
        INTEGER(intk) :: itask, ista, jsta, ksta, isto, jsto, ksto
        INTEGER(intk) :: iif, jjf, kkf, iic, jjc, kkc
        INTEGER(intk) :: igridc, igridf, ip3c, ip3f

        !$omp target teams distribute private(itask, igridf, igridc, &
        !$omp&  ista, jsta, ksta, isto, jsto, ksto, ip3f, ip3c, &
        !$omp&  kkf, jjf, iif, kkc, jjc, iic)
        DO itask = 1, nselftasks

            ! Unpacking the task
            igridf = etasks(itask)%igridf
            igridc = etasks(itask)%igridc
            ista = etasks(itask)%ista
            jsta = etasks(itask)%jsta
            ksta = etasks(itask)%ksta
            isto = etasks(itask)%isto
            jsto = etasks(itask)%jsto
            ksto = etasks(itask)%ksto

            ! Getting the parameters of the grids
            CALL get_ip3(ip3f, igridf)
            CALL get_ip3(ip3c, igridc)
            CALL get_mgdims(kkf, jjf, iif, igridf)
            CALL get_mgdims(kkc, jjc, iic, igridc)

            !$omp parallel
            CALL copy_kernel(kkf, jjf, iif, kkc, jjc, iic, &
                ff(ip3f), fc(ip3c), ista, jsta, ksta, isto, jsto, ksto)
            !$omp end parallel
        END DO
        !$omp end target teams distribute
    END SUBROUTINE process_selftasks_impl


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

        !$omp do collapse(3) private(i, j, k, ic, jc, kc)
        DO i = 1, iif
            DO j = 1, jjf
                DO k = 1, kkf

                    ! Computing the corresponding coarse grid indices
                    ic = (i-1)/2 + ista
                    jc = (j-1)/2 + jsta
                    kc = (k-1)/2 + ksta

#ifdef _MGLET_DEBUG_
                    ! Checking that the coarse grid indices are within bounds
                    IF (ic > isto .OR. jc > jsto .OR. kc > ksto) THEN
                        CALL errr(__FILE__, __LINE__)
                    END IF
#endif

                    ! Copying the value from the coarse grid to the fine grid
                    ff(k, j, i) = fc(kc, jc, ic)

                END DO
            END DO
        END DO
        !$omp end do
    END SUBROUTINE copy_kernel


    SUBROUTINE process_mpisendtasks(mpistasks, nmpistasks)
        ! Subroutine arguments
        INTEGER(intk), INTENT(in) :: nmpistasks
        TYPE(mpitask_t), INTENT(in) :: mpistasks(nmpistasks+1)

        ! Local variables
        INTEGER(intk) :: i, idx_sendbuf, messagelength, iprocf

        ! Leaving immediately if there are no tasks to process
        IF (nmpistasks < 1) RETURN

        CALL profile_range_push("process_mpisend")

        !$omp target data use_device_addr(sendbuf)
        DO i = 1, nmpistasks

            ! Getting connection information
            idx_sendbuf = mpistasks(i)%offset
            messagelength = mpistasks(i)%messagelength
            iprocf = mpistasks(i)%iprocnbr

            IF (iprocf == myid) CALL errr(__FILE__, __LINE__)

            ! Non-blocking MPI call with request handle stored in sendreqs
            CALL MPI_Isend(sendbuf(idx_sendbuf), messagelength, &
                mglet_mpi_real, iprocf, 1, MPI_COMM_WORLD, sendreqs(i))
        END DO
        !$omp end target data

        ! Checking for the dummy entry at position (end+1)
        IF (mpistasks(nmpistasks+1)%offset /= -1) THEN
            CALL errr(__FILE__, __LINE__)
        END IF

        ! Setting the number of posted receives
        nsend = nmpistasks

        CALL profile_range_pop()
    END SUBROUTINE process_mpisendtasks


    SUBROUTINE process_sendtasks(fc, nstasks, stasks)
        ! Subroutine arguments
        TYPE(field_t), TARGET, INTENT(in) :: fc
        INTEGER(intk), INTENT(in) :: nstasks
        TYPE(sendtask_t), INTENT(in) :: stasks(nstasks+1)

        ! Local variables
        ! none...

        ! Leaving immediately if there are no tasks to process
        IF (nstasks < 1) RETURN

        CALL profile_range_push("process_sendtasks")

        CALL process_sendtasks_impl(fc%arr, nstasks, stasks)

        CALL profile_range_pop()
        END SUBROUTINE process_sendtasks


        SUBROUTINE process_sendtasks_impl(coarse, nstasks, stasks)
        ! Subroutine arguments
        REAL(realk), INTENT(in) :: coarse(*)
        INTEGER(intk), INTENT(in) :: nstasks
        TYPE(sendtask_t), INTENT(in) :: stasks(nstasks+1)

        ! Local variables
        INTEGER(intk) :: igridc, itask, ii, jj, kk, ip3
        INTEGER(intk) :: ista, jsta, ksta, isto, jsto, ksto, idx1, idx2

        !$omp target teams distribute private(itask, ii, jj, kk, ip3, &
        !$omp&  ista, jsta, ksta, isto, jsto, ksto, idx1, idx2, igridc)
        DO itask = 1, nstasks

            ! Unpacking the task
            ista = stasks(itask)%ista
            jsta = stasks(itask)%jsta
            ksta = stasks(itask)%ksta
            isto = stasks(itask)%isto
            jsto = stasks(itask)%jsto
            ksto = stasks(itask)%ksto
            igridc = stasks(itask)%igridc
            idx1 = stasks(itask)%idx1
            idx2 = stasks(itask)%idx2

            ! Getting the parameters of the coarse grid
            CALL get_ip3(ip3, igridc)
            CALL get_mgdims(kk, jj, ii, igridc)

            !$omp parallel
            CALL write_buffer(kk, jj, ii, sendbuf(idx1:idx2), &
                coarse(ip3), ista, jsta, ksta, isto, jsto, ksto)
            !$omp end parallel
        END DO
        !$omp end target teams distribute
    END SUBROUTINE process_sendtasks_impl


    SUBROUTINE write_buffer(kk, jj, ii, buf, fc, ista, jsta, ksta, &
            isto, jsto, ksto)
        !$omp declare target

        ! Subroutine arguments
        INTEGER(intk), INTENT(in) :: kk, jj, ii
        REAL(realk), INTENT(inout) :: buf(*)
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

        !$omp do collapse(3) private(i, j, k, idx)
        DO i = ista, isto
            DO j = jsta, jsto
                DO k = ksta, ksto
                    idx = 1 + (k - ksta) + (j - jsta)*kkc + (i - ista)*kkc*jjc
                    buf(idx) = fc(k, j, i)
                END DO
            END DO
        END DO
        !$omp end do
    END SUBROUTINE write_buffer


    SUBROUTINE add_sendtask(stask, sendid, messagelength, sendcounter)
        ! Subroutine arguments
        TYPE(sendtask_t), INTENT(inout) :: stask
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
        stask%ista = ista
        stask%jsta = jsta
        stask%ksta = ksta
        stask%isto = isto
        stask%jsto = jsto
        stask%ksto = ksto
        stask%igridc = igridc
        stask%idx1 = offset
        stask%idx2 = offset + thismessagelength - 1

        messagelength = messagelength + thismessagelength
    END SUBROUTINE add_sendtask


    SUBROUTINE prepare_recvtasks(nrtasks, rtasks)
        ! Subroutine arguments
        INTEGER(intk), INTENT(out) :: nrtasks
        TYPE(recvtask_t), INTENT(inout) :: rtasks(:)

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

                        ! Adding new receive task
                        ! (= storing fine grid and offset in recv buffer)
                        irtasks = irtasks + 1
                        rtasks(irtasks)%igridf = igridf
                        rtasks(irtasks)%offset = offset

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
        rtasks(nrtasks+1) = recvtask_t(-1, -1)
    END SUBROUTINE prepare_recvtasks


    SUBROUTINE process_recvtasks(ff, nrtasks, rtasks)
        ! Subroutine arguments
        TYPE(field_t), TARGET, INTENT(inout) :: ff
        INTEGER(intk), INTENT(in) :: nrtasks
        TYPE(recvtask_t), INTENT(in) :: rtasks(nrtasks+1)

        ! Local variables
        ! none...

        ! Leaving immediately if there are no tasks to process
        IF (nrtasks < 1) RETURN

        CALL profile_range_push("process_recvtasks")

        CALL process_recvtasks_impl(ff%arr, nrtasks, rtasks)

        CALL profile_range_pop()
        END SUBROUTINE process_recvtasks


        SUBROUTINE process_recvtasks_impl(fine, nrtasks, rtasks)
        ! Subroutine arguments
        REAL(realk), INTENT(inout) :: fine(*)
        INTEGER(intk), INTENT(in) :: nrtasks
        TYPE(recvtask_t), INTENT(in) :: rtasks(nrtasks+1)

        ! Local variables
        INTEGER(intk) :: igridf, idx, len, itask, ii, jj, kk, ip3

        !$omp target teams distribute private(itask, ii, jj, kk, ip3, &
        !$omp&  igridf, idx, len)
        DO itask = 1, nrtasks

            ! Unpacking the task
            igridf = rtasks(itask)%igridf
            idx = rtasks(itask)%offset

            ! Getting parameters of fine grid
            CALL get_mgdims(kk, jj, ii, igridf)
            CALL get_ip3(ip3, igridf)
            len = kk * jj * ii / 8


            !$omp parallel
            ! Copying from buffer to the fine grid
            CALL write_fine(kk, jj, ii, fine(ip3), len, &
                recvbuf(idx:idx+len-1))
            !$omp end parallel
        END DO
        !$omp end target teams distribute
    END SUBROUTINE process_recvtasks_impl


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

        ! The counter for the buffer reflects the coarse grid
        kkc = kk/2
        jjc = jj/2
        iic = ii/2

        !$omp do collapse(3) private(i, j, k, idx, kc, jc, ic)
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
        !$omp end do
    END SUBROUTINE write_fine


    SUBROUTINE init_ctof2()
        ! Subroutine arguments
        ! none...

        ! Local variables
        ! none...

        IF (is_init) CALL errr(__FILE__, __LINE__)

        ! Check if ctof_core_mod necessary provides infrastructure
        IF (.NOT. is_ctof_core_init) CALL errr(__FILE__, __LINE__)

        ALLOCATE(sendreqs(isend))
        ALLOCATE(recvreqs(irecv))
        ALLOCATE(recvlist(irecv))
        ALLOCATE(recvidxlist(3, irecv))
        nrecv = 0
        nsend = 0

        ! Allocate the workrecords array for all possible ilevels
        ALLOCATE(workrecords(minlevel:maxlevel))

        is_init = .TRUE.

        ! Initialization is done, now record the workpackages
        CALL run_recording_pass()
    END SUBROUTINE init_ctof2


    SUBROUTINE run_recording_pass()
        ! Subroutine arguments
        ! none...

        ! Local variables
        INTEGER(intk) :: ilevel
        TYPE(field_t) :: dummy

        CALL dummy%init("DUMMY")

        !$omp target enter data map(to: dummy%arr)

        is_recording = .TRUE.

        ! Recording operations for all levels
        DO ilevel = minlevel+1, maxlevel
            CALL ctof2(ilevel, dummy, dummy)
        END DO

        is_recording = .FALSE.

        !$omp target exit data map(delete: dummy%arr)

        CALL dummy%finish()
    END SUBROUTINE run_recording_pass


    SUBROUTINE finish_ctof2()
        ! Subroutine arguments
        ! none...

        ! Local variables
        INTEGER(intk) :: ilevel

        IF (.NOT. is_init) CALL errr(__FILE__, __LINE__)

        DEALLOCATE(sendreqs)
        DEALLOCATE(recvreqs)
        DEALLOCATE(recvlist)
        DEALLOCATE(recvidxlist)

        ! Deallocate the workpackage components for each level
        DO ilevel = minlevel+1, maxlevel
            IF (workrecords(ilevel)%is_init) THEN
                !$omp target exit data map(delete: &
                !$omp&  workrecords(ilevel)%sendtasks(:), &
                !$omp&  workrecords(ilevel)%recvtasks(:), &
                !$omp&  workrecords(ilevel)%selftasks(:))

                DEALLOCATE(workrecords(ilevel)%sendtasks)
                DEALLOCATE(workrecords(ilevel)%recvtasks)
                DEALLOCATE(workrecords(ilevel)%selftasks)
                DEALLOCATE(workrecords(ilevel)%mpisendtasks)
                DEALLOCATE(workrecords(ilevel)%mpirecvtasks)
                workrecords(ilevel)%is_init = .FALSE.
            END IF
        END DO
        DEALLOCATE(workrecords)

        is_init = .FALSE.
    END SUBROUTINE finish_ctof2
END MODULE ctof2_mod
