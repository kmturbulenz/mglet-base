MODULE ctof2_mod

    USE core_mod
    USE MPI_f08

    IMPLICIT NONE (type, external)
    PRIVATE

    ! Variable to indicate if the required data structures and MPI-types
    ! have been created
    LOGICAL :: is_init = .FALSE.

    ! If .TRUE., a prolongation is already in process and you cannot start
    ! another one
    LOGICAL :: in_progress = .FALSE.

    ! Maximum allowed number of childs per parent (i.e. maximum number of
    ! send-conenctions per grid)
    INTEGER(intk), PARAMETER :: maxchilds = 8

    ! Lists that hold the send and receive request arrays
    TYPE(MPI_Request), ALLOCATABLE :: sendreqs(:), recvreqs(:)

    ! Actual number of messages that are to be sendt and received in one
    ! "round" of operations
    INTEGER(intk) :: nsend, nrecv

    ! List of grids to receive data on
    INTEGER(intk), ALLOCATABLE :: recvgrids(:), recvpos(:)

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
    INTEGER(intk), PARAMETER :: mpitasksize = 4
    INTEGER(intk), PARAMETER :: sendtasksize = 9
    INTEGER(intk), PARAMETER :: recvtasksize = 2

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
        INTEGER(intk) :: nmpirecvtasks, nmpisendtasks, nsendtasks, nrecvtasks

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

            ! Obtaining numbers of tasks from the workpackage
            nmpirecvtasks = SIZE(wptr%mpirecvtasks, 2) - 1
            nmpisendtasks = SIZE(wptr%mpisendtasks, 2) - 1
            nsendtasks = SIZE(wptr%sendtasks, 2) - 1
            nrecvtasks = SIZE(wptr%recvtasks, 2) - 1

            CALL process_mpirecvtasks(wptr%mpirecvtasks, nmpirecvtasks)
            CALL process_sendtasks(wptr%sendtasks, nsendtasks)
            CALL process_mpisendtasks(wptr%mpisendtasks, nmpisendtasks)
            CALL MPI_Waitall(nrecv, recvreqs, MPI_STATUSES_IGNORE)
            CALL process_recvtasks(wptr%recvtasks, nrecvtasks)
            CALL MPI_Waitall(nsend, sendreqs, MPI_STATUSES_IGNORE)

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
            ALLOCATE(mpisendtasks(mpitasksize, maxsize))
            ALLOCATE(mpirecvtasks(mpitasksize, maxsize))

            ! Alibi action
            ALLOCATE(selftasks(1, 1))
            DEALLOCATE(selftasks)

            CALL prepare_sendtasks(sendtasks, nsendtasks, ilevel)
            CALL prepare_mpirecvtasks(mpirecvtasks, nmpirecvtasks, ilevel)
            CALL prepare_mpisendtasks(mpisendtasks, nmpisendtasks, ilevel)

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

            ! At his point, one execution has been performed.

            ! Allocating persistent workpackage with accurate dimensions
            ALLOCATE(wptr%sendtasks(sendtasksize, nsendtasks+1))
            ALLOCATE(wptr%recvtasks(recvtasksize, nrecvtasks+1))
            ALLOCATE(wptr%mpisendtasks(mpitasksize, nmpisendtasks+1))
            ALLOCATE(wptr%mpirecvtasks(mpitasksize, nmpirecvtasks+1))

            ! Tranfering the recorded tasks to the persistent workpackage
            wptr%sendtasks(:, 1:nsendtasks+1) = sendtasks(:, 1:nsendtasks+1)
            wptr%recvtasks(:, 1:nrecvtasks+1) = recvtasks(:, 1:nrecvtasks+1)
            wptr%mpisendtasks(:, 1:nmpisendtasks+1) = &
                mpisendtasks(:, 1:nmpisendtasks+1)
            wptr%mpirecvtasks(:, 1:nmpirecvtasks+1) = &
                mpirecvtasks(:, 1:nmpirecvtasks+1)

            !$omp target enter data map(to: &
            !$omp&  wptr%sendtasks(1:sendtasksize, 1:nsendtasks+1), &
            !$omp&  wptr%recvtasks(1:recvtasksize, 1:nrecvtasks+1))

            wptr%is_init = .TRUE.

            ! Deallocating the temporary arrays (already gone on device)
            DEALLOCATE(sendtasks, recvtasks, mpisendtasks, mpirecvtasks)

        END IF

        END ASSOCIATE

    END SUBROUTINE ctof2


    !
    SUBROUTINE prepare_mpirecvtasks(mpirtasks, nmpirtasks, ilevel)

        ! Subroutine arguments
        INTEGER(intk), INTENT(inout) :: mpirtasks(4, maxsize)
        INTEGER(intk), INTENT(out) :: nmpirtasks
        INTEGER(intk), INTENT(in) :: ilevel   ! Level of the *fine* side

        ! Local variables
        INTEGER(intk) :: i, igridf, igridc, iprocc, iprocf
        INTEGER(intk) :: kk, jj, ii, impirtasks, idx_recvbuf
        INTEGER(int32) :: recvcounter, messagelength

        ! Post all receive calls
        recvcounter = 0
        messagelength = 0
        nrecv = 0

        ! Initialize the internal counter
        impirtasks = 0

        DO i = 1, noflevel(ilevel)
            igridf = igrdoflevel(i, ilevel)
            igridc = iparent(igridf)
            IF (igridc == 0) CYCLE

            iprocc = idprocofgrd(igridc)
            iprocf = idprocofgrd(igridf)

            IF (myid == iprocf) THEN

                IF (iprocf == iprocc) THEN
                    ! This is a self-connection, no MPI communication is needed
                    WRITE(*, *) "Self-connection on rank ", myid
                END IF

                nrecv = nrecv + 1

                CALL get_mgdims(kk, jj, ii, igridf)
                messagelength = kk*jj*ii/8
                idx_recvbuf = recvcounter + 1

                IF (recvcounter + messagelength > SIZE(sendbuf)) THEN
                    CALL errr(__FILE__, __LINE__)
                END IF

                ! Adding new recevive task
                impirtasks = impirtasks + 1
                mpirtasks(1, impirtasks) = idx_recvbuf
                mpirtasks(2, impirtasks) = messagelength
                mpirtasks(3, impirtasks) = iprocc
                mpirtasks(4, impirtasks) = igridf

                recvgrids(nrecv) = igridf
                recvpos(nrecv) = recvcounter + 1
                recvcounter = recvcounter + messagelength
            END IF
        END DO

        ! Assign the number of tasks to the output variable
        nmpirtasks = impirtasks

        ! Adding a dummy entry at position (end+1)
        mpirtasks(:, nmpirtasks+1) = -1

    END SUBROUTINE prepare_mpirecvtasks


    SUBROUTINE process_mpirecvtasks(mpirtasks, nmpirtasks)
        ! Subroutine arguments
        INTEGER(intk), INTENT(in) :: nmpirtasks
        INTEGER(intk), INTENT(in) :: mpirtasks(mpitasksize, nmpirtasks+1)


        ! Local variables
        INTEGER(intk) :: i, idx_recvbuf, messagelength, iprocc, igridf

        ! Positing all non-blocking receive calls
        DO i = 1, nmpirtasks

            idx_recvbuf = mpirtasks(1, i)
            messagelength = mpirtasks(2, i)
            iprocc = mpirtasks(3, i)
            igridf = mpirtasks(4, i)

            CALL MPI_Irecv(recvbuf(idx_recvbuf), messagelength, &
                mglet_mpi_real, iprocc, igridf, MPI_COMM_WORLD, &
                recvreqs(i))
        END DO

        ! Checking for the dummy entry at position (end+1)
        IF (mpirtasks(1, nmpirtasks+1) /= -1) THEN
            CALL errr(__FILE__, __LINE__)
        END IF

        ! Setting the number of posted receives
        nrecv = nmpirtasks

    END SUBROUTINE process_mpirecvtasks


    ! Perform all Send-calls
    SUBROUTINE prepare_mpisendtasks(mpistasks, nmpistasks, ilevel)
        ! Subroutine arguments
        INTEGER(intk), INTENT(inout) :: mpistasks(mpitasksize, maxsize)
        INTEGER(intk), INTENT(out) :: nmpistasks
        INTEGER(intk), INTENT(in) :: ilevel   ! Level of the *fine* side

        ! Local variables
        INTEGER(intk) :: i, igridf, igridc, iprocc, iprocf
        INTEGER(intk) :: kk, jj, ii, idx_sendbuf
        INTEGER(intk) :: kkf, jjf, iif, impistasks
        INTEGER(int32) :: sendcounter, messagelength

        ! Post all receive calls
        sendcounter = 0
        messagelength = 0
        nsend = 0

        ! Initialize the internal counter
        impistasks = 0

        DO i = 1, noflevel(ilevel)
            igridf = igrdoflevel(i, ilevel)
            igridc = iparent(igridf)
            IF (igridc == 0) CYCLE

            iprocc = idprocofgrd(igridc)
            iprocf = idprocofgrd(igridf)

            IF (myid == iprocc) THEN

                IF (iprocf == iprocc) THEN
                    ! This is a self-connection, no MPI communication is needed
                    WRITE(*, *) "Self-connection on rank ", myid
                END IF

                nsend = nsend + 1

                CALL get_mgdims(kk, jj, ii, igridc)
                CALL get_mgdims(kkf, jjf, iif, igridf)
                messagelength = kkf*jjf*iif/8
                idx_sendbuf = sendcounter + 1

                IF (sendcounter + messagelength > SIZE(sendbuf)) THEN
                    CALL errr(__FILE__, __LINE__)
                END IF

                ! Adding new recevive task
                impistasks = impistasks + 1
                mpistasks(1, impistasks) = idx_sendbuf
                mpistasks(2, impistasks) = messagelength
                mpistasks(3, impistasks) = iprocf
                mpistasks(4, impistasks) = igridf

                sendcounter = sendcounter + messagelength
            END IF
        END DO

        ! Assign the number of tasks to the output variable
        nmpistasks = impistasks

        ! Adding a dummy entry at position (end+1)
        mpistasks(:, nmpistasks+1) = -1

    END SUBROUTINE prepare_mpisendtasks


    SUBROUTINE process_mpisendtasks(mpistasks, nmpistasks)
        ! Subroutine arguments
        INTEGER(intk), INTENT(in) :: nmpistasks
        INTEGER(intk), INTENT(in) :: mpistasks(mpitasksize, nmpistasks+1)

        ! Local variables
        INTEGER(intk) :: i, idx_sendbuf, messagelength, iprocf, igridf

        ! Positing all non-blocking send calls
        DO i = 1, nmpistasks

            idx_sendbuf = mpistasks(1, i)
            messagelength = mpistasks(2, i)
            iprocf = mpistasks(3, i)
            igridf = mpistasks(4, i)

            CALL MPI_Isend(sendbuf(idx_sendbuf), messagelength, &
                mglet_mpi_real, iprocf, igridf, MPI_COMM_WORLD, &
                sendreqs(i))
        END DO

        ! Checking for the dummy entry at position (end+1)
        IF (mpistasks(1, nmpistasks+1) /= -1) THEN
            CALL errr(__FILE__, __LINE__)
        END IF

        ! Setting the number of posted receives
        nsend = nmpistasks

    END SUBROUTINE process_mpisendtasks




    ! Prepare all tasks needed to prepare the buffers before sending
    !
    SUBROUTINE prepare_sendtasks(stasks, nstasks, ilevel)
        ! Subroutine arguments
        INTEGER(intk), INTENT(inout) :: stasks(:, :)
        INTEGER(intk), INTENT(out) :: nstasks
        INTEGER(intk), INTENT(in) :: ilevel   ! Level of the *fine* grid

        ! Local variables
        INTEGER(intk) :: i, igridf, igridc, iprocc, iprocf, istasks
        INTEGER(intk) :: kk, jj, ii
        INTEGER(intk) :: kkf, jjf, iif
        INTEGER(int32) :: sendcounter, messagelength

        ! Post all receive calls
        sendcounter = 0
        messagelength = 0

        ! Initialize the internal counter
        istasks = 0

        DO i = 1, noflevel(ilevel)
            igridf = igrdoflevel(i, ilevel)
            igridc = iparent(igridf)
            IF (igridc == 0) CYCLE

            iprocc = idprocofgrd(igridc)
            iprocf = idprocofgrd(igridf)

            IF (myid == iprocc) THEN

                CALL get_mgdims(kk, jj, ii, igridc)
                CALL get_mgdims(kkf, jjf, iif, igridf)
                messagelength = kkf*jjf*iif/8

                IF (sendcounter + messagelength > SIZE(sendbuf)) THEN
                    CALL errr(__FILE__, __LINE__)
                END IF

                ! Adding new send task (= writing to the buffer)
                istasks = istasks + 1

                CALL add_sendtask(stasks(:, istasks), sendcounter+1, &
                    sendcounter+messagelength, igridc, igridf)

                sendcounter = sendcounter + messagelength
            END IF
        END DO

        ! Assign the number of tasks to the output variable
        nstasks = istasks

        ! Adding a dummy entry at position (end+1)
        stasks(:, nstasks+1) = -1

    END SUBROUTINE prepare_sendtasks

    SUBROUTINE add_sendtask(stask, idx_sendbuf_start, idx_sendbuf_stop, &
            igridc, igridf)

        ! Subroutine arguments
        INTEGER(intk), INTENT(inout) :: stask(sendtasksize)
        INTEGER(intk), INTENT(in) :: idx_sendbuf_start, idx_sendbuf_stop,&
            igridc, igridf

        ! Local variables
        INTEGER(intk) :: kkf, jjf, iif, counter
        INTEGER(intk) :: ista, jsta, ksta, isto, jsto, ksto

        ! Compute start- and end-positions in coarse grid
        ista = iposition(igridf) - 1
        jsta = jposition(igridf) - 1
        ksta = kposition(igridf) - 1

        CALL get_mgdims(kkf, jjf, iif, igridf)
        isto = ista + (iif - 4)/2 + 1
        jsto = jsta + (jjf - 4)/2 + 1
        ksto = ksta + (kkf - 4)/2 + 1

        ! Sanity checks
        counter = (isto- ista + 1) * (jsto - jsta + 1) * (ksto - ksta + 1)
        IF (counter /= kkf*jjf*iif/8) THEN
            WRITE(*, *) "counter = ", counter
            WRITE(*, *) "kkf = ", kkf
            WRITE(*, *) "jjf = ", jjf
            WRITE(*, *) "iif = ", iif
            WRITE(*, *) "ksta = ", ksta
            WRITE(*, *) "jsta = ", jsta
            WRITE(*, *) "ista = ", ista
            WRITE(*, *) "ksto = ", ksto
            WRITE(*, *) "jsto = ", jsto
            WRITE(*, *) "isto = ", isto
            CALL errr(__FILE__, __LINE__)
        END IF

        IF (counter /= idx_sendbuf_stop - idx_sendbuf_start + 1) THEN
            CALL errr(__FILE__, __LINE__)
        END IF

        ! Writing the task
        stask(1) = ista
        stask(2) = jsta
        stask(3) = ksta
        stask(4) = isto
        stask(5) = jsto
        stask(6) = ksto
        stask(7) = igridc
        stask(8) = idx_sendbuf_start
        stask(9) = idx_sendbuf_stop

    END SUBROUTINE add_sendtask


    SUBROUTINE process_sendtasks(stasks, nstasks)
        ! Subroutine arguments
        INTEGER(intk), INTENT(in) :: nstasks
        INTEGER(intk), INTENT(in) :: stasks(sendtasksize, nstasks+1)

        ! Local variables
        INTEGER(intk) :: igridc, itask, ii, jj, kk, ip3, err
        INTEGER(intk) :: ista, jsta, ksta, isto, jsto, ksto, idx1, idx2

        IF (nstasks == 0) RETURN

        ASSOCIATE(coarse => fc%arr)

        !$omp target teams distribute private(itask, ii, jj, kk, ip3, &
        !$omp&  ista, jsta, ksta, isto, jsto, ksto, idx1, idx2, igridc) &
        !$omp$  map(from: err)
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

            ! Checking for the dummy entry at position (end+1)
            IF (stasks(1, nstasks+1) /= -1) THEN
                err = 1
            ELSE
                err = 0
            END IF

        END DO
        !$omp end target teams distribute

        IF (err /= 0) THEN
            CALL errr(__FILE__, __LINE__)
        END IF

        END ASSOCIATE

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


    SUBROUTINE prepare_recvtasks(rtasks, nrtasks)
        ! Subroutine arguments
        INTEGER(intk), INTENT(inout) :: rtasks(:, :)
        INTEGER(intk), INTENT(out) :: nrtasks

        ! Local variables
        INTEGER(int32) :: idx
        INTEGER(intk) :: irtasks, igridf, idx_recvbuf

        ! Initialize the internal counter
        irtasks = 0

        IF (nrecv > 0) THEN
            DO WHILE (.TRUE.)
                CALL MPI_Waitany(nrecv, recvreqs, idx, MPI_STATUS_IGNORE)

                IF (idx /= MPI_UNDEFINED) THEN

                    ! Adding new recevive task
                    irtasks = irtasks + 1
                    igridf = recvgrids(idx)
                    idx_recvbuf = recvpos(idx)

                    CALL add_recvtask(rtasks(:, irtasks), igridf, idx_recvbuf)
                ELSE
                    EXIT
                END IF
            END DO
        END IF

        CALL MPI_Waitall(nsend, sendreqs, MPI_STATUSES_IGNORE)

        ! Assigning the total count
        nrtasks = irtasks

        ! Adding a dummy entry at position (end+1)
        rtasks(:, nrtasks+1) = -1

    END SUBROUTINE prepare_recvtasks


    SUBROUTINE add_recvtask(rtask, igridf, idx_recvbuf)

        ! Subroutine arguments
        INTEGER(intk), INTENT(inout) :: rtask(2)
        INTEGER(intk), INTENT(in) :: igridf
        INTEGER(intk), INTENT(in) :: idx_recvbuf

        ! Local variables
        INTEGER(intk) :: kk, jj, ii

        CALL get_mgdims(kk, jj, ii, igridf)

        ! Writing the task
        rtask(1) = igridf
        rtask(2) = idx_recvbuf

    END SUBROUTINE add_recvtask



    SUBROUTINE process_recvtasks(rtasks, nrtasks)
        ! Subroutine arguments
        INTEGER(intk), INTENT(in) :: nrtasks
        INTEGER(intk), INTENT(in) :: rtasks(recvtasksize, nrtasks+1)

        ! Local variables
        INTEGER(intk) :: igridf, idx, len, itask, ii, jj, kk, ip3, err

        IF (nrtasks == 0) RETURN

        ASSOCIATE(fine => ff%arr)

        !$omp target teams distribute private(itask, ii, jj, kk, ip3, &
        !$omp&  igridf, idx, len) map(from: err)
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

            ! Checking for the dummy entry at position (end+1)
            IF (rtasks(1, nrtasks+1) /= -1) THEN
                err = 1
            ELSE
                err = 0
            END IF

        END DO
        !$omp end target teams distribute

        IF (err /= 0) THEN
            CALL errr(__FILE__, __LINE__)
        END IF

        END ASSOCIATE

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
        INTEGER(intk) :: igrid, iprocc, ipar, ilevel
        TYPE(field_t) :: dummy

        CALL set_timer(230, "CTOF")
        CALL set_timer(231, "CTOF_BEGIN")
        CALL set_timer(232, "CTOF_END")
        CALL set_timer(235, "CTOF_PROLONG_FINISH")

        IF (.NOT. is_init) THEN
            ALLOCATE(sendreqs(nmygrids*maxchilds))
            ALLOCATE(recvreqs(nmygrids))
            ALLOCATE(recvgrids(nmygrids))
            ALLOCATE(recvpos(nmygrids))
        END IF

        nrecv = 0
        nsend = 0

        DO igrid = 1, ngrid
            ipar = iparent(igrid)
            IF (ipar /= 0) THEN
                iprocc = idprocofgrd(ipar)
                IF (myid == iprocc) THEN
                    nsend = nsend + 1
                    IF (nsend > nmygrids*maxchilds) THEN
                        CALL errr(__FILE__, __LINE__)
                    END IF
                END IF
            END IF
        END DO

        is_init = .TRUE.
        in_progress = .FALSE.

        ! Allocate the workrecords array for all possible ilevels
        ALLOCATE(workrecords(minlevel:maxlevel))
        CALL dummy%init("DUMMY")

        ! Recording operations for all levels
        WRITE(*, *) "Starting to record prolongation operations for all levels"
        is_recording = .TRUE.
        DO ilevel = minlevel, maxlevel
            CALL ctof2(ilevel, dummy, dummy)
        END DO
        is_recording = .FALSE.
        WRITE(*, *) "Finished recording prolongation operations for all levels"

        CALL dummy%finish()

    END SUBROUTINE init_ctof2



    SUBROUTINE finish_ctof2()
        INTEGER(intk) :: ilevel

        IF (is_init .NEQV. .TRUE.) THEN
            RETURN
        END IF

        is_init = .FALSE.

        DEALLOCATE(sendreqs)
        DEALLOCATE(recvreqs)
        DEALLOCATE(recvgrids)
        DEALLOCATE(recvpos)

        ! Deallocate the workpackage components for each level
        DO ilevel = minlevel, maxlevel
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
