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
        selftasks(:, :), mpisendtasks(:, :), mpirecvtasks(:, :)
    !$omp declare target(sendtasks, recvtasks, selftasks)

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


    PUBLIC :: ctof2, init_ctof2, finish_ctof2

CONTAINS

    SUBROUTINE ctof2(ilevel, ff, fc)

        ! Subroutine arguments
        INTEGER(intk), INTENT(in) :: ilevel  ! Level of the *fine* side
        REAL(realk), INTENT(inout) :: ff(:)
        REAL(realk), INTENT(in) :: fc(:)
        ! Local variables
        ! none...

        CALL start_timer(230)

        ! Looking up the workpackage for this level
        ASSOCIATE(wptr => workrecords(ilevel))

        IF (.NOT. is_recording) THEN

            ! During the execution phase, the workpackage is already initialized
            IF (.NOT. wptr%is_init) THEN
                WRITE(*, *) "Routine ctof2 does not provide just-in-time"
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
            maxsize = noflevel(ilevel)
            ALLOCATE(sendtasks(sendtasksize, maxsize))
            ALLOCATE(recvtasks(recvtasksize, maxsize))
            ALLOCATE(mpisendtasks(mpitasksize, maxsize))
            ALLOCATE(mpirecvtasks(mpitasksize, maxsize))

            CALL prepare_sendtasks(sendtasks, nsendtasks, ilevel)
            CALL prepare_mpirecvtasks(mpirecvtasks, nmpirecvtasks, ilevel)
            CALL prepare_mpisendtasks(mpisendtasks, nmpisendtasks, ilevel)

            CALL process_sendtasks(sendtasks, nsendtasks)
            CALL process_mpirecvtasks(mpirecvtasks, nmpirecvtasks)
            CALL process_mpisendtasks(mpisendtasks, nmpisendtasks)

            ! Includes waiting for MPI communication to finish
            CALL prepare_recvtasks(recvtasks, nrecvtasks, ilevel)
            CALL process_recvtasks(recvtasks, nrecvtasks)

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

            !$omp target update to(wptr%sendtasks, wptr%recvtasks)

            wptr%is_init = .TRUE.

            ! Deallocating the temporary arrays
            DEALLOCATE(sendtasks, recvtasks, mpisendtasks, mpirecvtasks)

        END IF

        END ASSOCIATE

        CALL stop_timer(230)

    END SUBROUTINE ctof2


    ! Initiate prolongation
    !
    ! Interpolate the results from a coarse level to a fine level
    ! This function initiate the process, and the ctof_finish
    ! must be called afterwards to clean up.
    SUBROUTINE ctof_begin(ilevel, fc)
        ! Subroutine arguments
        INTEGER(intk), INTENT(in) :: ilevel
        REAL(realk), INTENT(in) :: fc(:)

        ! Local variables
        ! none...

        CALL start_timer(231)

        IF (.NOT. is_init) THEN
            CALL errr(__FILE__, __LINE__)
        END IF
        IF (in_progress) THEN
            CALL errr(__FILE__, __LINE__)
        END IF
        in_progress = .TRUE.

        CALL recv_all(ilevel)
        CALL send_all(ilevel, fc)

        CALL stop_timer(231)
    END SUBROUTINE ctof_begin


    ! Record information for the quick posting of MPI receives
    !
    SUBROUTINE prepare_mpirecvtasks(mpirtasks, nmpirtasks, ilevel)

        ! Subroutine arguments
        INTEGER(intk), INTENT(inout) :: mpirtasks(4, :)
        INTEGER(intk), INTENT(out) :: nmpirtasks
        INTEGER(intk), INTENT(in) :: ilevel   ! Level of the *fine* side

        ! Local variables
        INTEGER(intk) :: i, igridf, igridc, iprocc, iprocf
        INTEGER(intk) :: kk, jj, ii, impirtasks
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
        mpirtasks(nmpirtasks+1, :) = -1

    END SUBROUTINE prepare_mpirecvtasks


    SUBROUTINE process_mpirecvtasks(mpirtasks, nmpirtasks)
        ! Subroutine arguments
        INTEGER(intk), INTENT(in) :: mpirtasks(mpitasksize, :)
        INTEGER(intk), INTENT(in) :: nmpirtasks

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
        IF (mpirtasks(nmpirtasks+1, 1) /= -1) THEN
            CALL errr(__FILE__, __LINE__)
        END IF

        ! Setting the number of posted receives
        nrecv = nmpirtasks

    END SUBROUTINE process_mpirecvtasks


    ! Perform all Send-calls
    SUBROUTINE prepare_mpisendtasks(mpistasks, nmpistasks, ilevel)
        ! Subroutine arguments
        INTEGER(intk), INTENT(inout) :: mpistasks(:, :)
        INTEGER(intk), INTENT(out) :: nmpistasks
        INTEGER(intk), INTENT(in) :: ilevel   ! Level of the *fine* side

        ! Local variables
        INTEGER(intk) :: i, igridf, igridc, iprocc, iprocf, ip3
        INTEGER(intk) :: kk, jj, ii
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

                CALL MPI_Isend(sendbuf(sendcounter+1), messagelength, &
                    mglet_mpi_real, iprocf, igridf, MPI_COMM_WORLD, &
                    sendreqs(nsend))

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
        mpistasks(nmpistasks+1, :) = -1

    END SUBROUTINE prepare_mpisendtasks


    SUBROUTINE process_mpisendtasks(mpistasks, nmpistasks)
        ! Subroutine arguments
        INTEGER(intk), INTENT(in) :: mpistasks(mpitasksize, :)
        INTEGER(intk), INTENT(in) :: nmpistasks

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
        IF (mpistasks(nmpistasks+1, 1) /= -1) THEN
            CALL errr(__FILE__, __LINE__)
        END IF

        ! Setting the number of posted receives
        nsend = nmpistasks

    END SUBROUTINE process_mpisendtasks






    ! Perform all Recv-calls
    SUBROUTINE recv_all(ilevel)
        ! Subroutine arguments
        INTEGER(intk), INTENT(in) :: ilevel   ! Level of the *fine* side

        ! Local variables
        INTEGER(intk) :: i, igridf, igridc, iprocc, iprocf
        INTEGER(intk) :: kk, jj, ii
        INTEGER(int32) :: recvcounter, messagelength

        ! Post all receive calls
        recvcounter = 0
        messagelength = 0
        nrecv = 0

        DO i = 1, noflevel(ilevel)
            igridf = igrdoflevel(i, ilevel)
            igridc = iparent(igridf)
            IF (igridc == 0) CYCLE

            iprocc = idprocofgrd(igridc)
            iprocf = idprocofgrd(igridf)

            IF (myid == iprocf) THEN
                nrecv = nrecv + 1

                CALL get_mgdims(kk, jj, ii, igridf)
                messagelength = kk*jj*ii/8

                IF (recvcounter + messagelength > SIZE(sendbuf)) THEN
                    CALL errr(__FILE__, __LINE__)
                END IF

                CALL MPI_Irecv(recvbuf(recvcounter+1), messagelength, &
                    mglet_mpi_real, iprocc, igridf, MPI_COMM_WORLD, &
                    recvreqs(nrecv))

                recvgrids(nrecv) = igridf
                recvpos(nrecv) = recvcounter + 1
                recvcounter = recvcounter + messagelength
            END IF
        END DO
    END SUBROUTINE recv_all


    ! Perform all Send-calls
    SUBROUTINE send_all(ilevel, fc)
        ! Subroutine arguments
        INTEGER(intk), INTENT(in) :: ilevel   ! Level of the *fine* grid
        REAL(realk), INTENT(in) :: fc(*)      ! Field on the coarse grid

        ! Local variables
        INTEGER(intk) :: i, igridf, igridc, iprocc, iprocf, ip3
        INTEGER(intk) :: kk, jj, ii
        INTEGER(intk) :: kkf, jjf, iif
        INTEGER(int32) :: sendcounter, messagelength

        ! Post all receive calls
        sendcounter = 0
        messagelength = 0
        nsend = 0

        DO i = 1, noflevel(ilevel)
            igridf = igrdoflevel(i, ilevel)
            igridc = iparent(igridf)
            IF (igridc == 0) CYCLE

            iprocc = idprocofgrd(igridc)
            iprocf = idprocofgrd(igridf)

            IF (myid == iprocc) THEN
                nsend = nsend + 1

                CALL get_mgdims(kk, jj, ii, igridc)
                CALL get_mgdims(kkf, jjf, iif, igridf)
                messagelength = kkf*jjf*iif/8

                IF (sendcounter + messagelength > SIZE(sendbuf)) THEN
                    CALL errr(__FILE__, __LINE__)
                END IF

                CALL get_ip3(ip3, igridc)
                CALL pack_send(&
                    sendbuf(sendcounter+1:sendcounter+messagelength), &
                    kk, jj, ii, fc(ip3), igridc, igridf)

                CALL MPI_Isend(sendbuf(sendcounter+1), messagelength, &
                    mglet_mpi_real, iprocf, igridf, MPI_COMM_WORLD, &
                    sendreqs(nsend))

                sendcounter = sendcounter + messagelength
            END IF
        END DO
    END SUBROUTINE send_all


    ! Prepare all tasks needed to prepare the buffers before sending
    !
    SUBROUTINE prepare_sendtasks(stasks, nstasks, ilevel)
        ! Subroutine arguments
        INTEGER(intk), INTENT(inout) :: stasks(:, :)
        INTEGER(intk), INTENT(out) :: nstasks
        INTEGER(intk), INTENT(in) :: ilevel   ! Level of the *fine* grid

        ! Local variables
        INTEGER(intk) :: i, igridf, igridc, iprocc, iprocf
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
        stasks(nstasks+1, :) = -1

    END SUBROUTINE prepare_sendtasks

    SUBROUTINE add_sendtask(stask, idx_sendbuf_start, idx_sendbuf_stop, &
            igridc, igridf)

        ! Subroutine arguments
        INTEGER(intk), INTENT(inout) :: stask(9)
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

        ! Post all receive calls
        sendcounter = 0
        messagelength = 0
        nsend = 0

    END SUBROUTINE add_sendtask


    SUBROUTINE process_sendtasks(stasks, nstasks)
        ! Subroutine arguments
        INTEGER(intk), INTENT(in) :: stasks(9, :)
        INTEGER(intk), INTENT(in) :: nstasks

        ! Local variables
        INTEGER(intk) :: igridf, idx_recvbuf, itask, ii, jj, kk, ip3

        ASSOCIATE(coarse => fc%arr)

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

        END ASSOCIATE

    END SUBROUTINE process_sendtasks


    SUBROUTINE write_buffer(kk, jj, ii, buf, fc, ista, jsta, ksta, &
            isto, jsto, ksto)

        ! Subroutine arguments
        INTEGER(intk), INTENT(in) :: kk, jj, ii
        REAL(realk), INTENT(inout) :: buf(:)
        REAL(realk), INTENT(in) :: fc(kk, jj, ii)
        INTEGER(intk), INTENT(in) :: ista, jsta, ksta
        INTEGER(intk), INTENT(in) :: isto, jsto, ksto

        ! Local variables
        INTEGER(intk) :: i, j, k
        INTEGER(intk) :: counter

        ! Compute extent of the patch used on the coarse grid
        iic = isto - ista + 1
        jjc = jsto - jsta + 1
        kkc = ksto - ksta + 1

        ! Pack buffer
        DO i = ista, isto
            DO j = jsta, jsto
                DO k = ksta, ksto
                    idx = 1 + (k - ksta) + (j - jsta)*kkc + (i - ista)*kkc*jjc
                    buf(idx) = fc(k, j, i)
                END DO
            END DO
        END DO

    END SUBROUTINE write_buffer


    SUBROUTINE pack_send(buf, kk, jj, ii, fc, igridc, igridf)
        ! Subroutine arguments
        REAL(realk), INTENT(inout) :: buf(:)
        INTEGER(intk), INTENT(in) :: kk, jj, ii
        REAL(realk), INTENT(in) :: fc(kk, jj, ii)
        INTEGER(intk), INTENT(in) :: igridc
        INTEGER(intk), INTENT(in) :: igridf

        ! Local variables
        INTEGER(intk) :: i, j, k
        INTEGER(intk) :: kkf, jjf, iif
        INTEGER(intk) :: counter
        INTEGER(intk) :: ista, jsta, ksta, isto, jsto, ksto

        ! Compute start- and end-positions in coarse grid
        ista = iposition(igridf) - 1
        jsta = jposition(igridf) - 1
        ksta = kposition(igridf) - 1

        CALL get_mgdims(kkf, jjf, iif, igridf)
        isto = ista + (iif - 4)/2 + 1
        jsto = jsta + (jjf - 4)/2 + 1
        ksto = ksta + (kkf - 4)/2 + 1

        ! Pack buffer
        counter = 0
        DO i = ista, isto
            DO j = jsta, jsto
                DO k = ksta, ksto
                    counter = counter + 1
                    buf(counter) = fc(k, j, i)
                END DO
            END DO
        END DO

        ! Sanity checks
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
        IF (counter /= SIZE(buf)) THEN
            CALL errr(__FILE__, __LINE__)
        END IF
    END SUBROUTINE pack_send


    ! Finish prolongation
    !
    ! Wait for communication to finish and clean up
    SUBROUTINE ctof_end(ff)
        ! Subroutine arguments
        REAL(realk), INTENT(inout) :: ff(:)

        ! Local variables
        INTEGER(int32) :: idx

        CALL start_timer(232)

        IF (.NOT. in_progress) THEN
            CALL errr(__FILE__, __LINE__)
        END IF

        IF (nrecv > 0) THEN
            DO WHILE (.TRUE.)
                CALL MPI_Waitany(nrecv, recvreqs, idx, MPI_STATUS_IGNORE)

                IF (idx /= MPI_UNDEFINED) THEN
                    CALL start_timer(235)
                    CALL prolong_finish(ff, recvgrids(idx), recvpos(idx))
                    CALL stop_timer(235)
                ELSE
                    EXIT
                END IF
            END DO
        END IF

        CALL MPI_Waitall(nsend, sendreqs, MPI_STATUSES_IGNORE)

        in_progress = .FALSE.

        CALL stop_timer(232)
    END SUBROUTINE ctof_end


    ! ! Finish prolongation, i.e. distribute the data on the grid
    ! SUBROUTINE prolong_finish(ff, igridf, pos)
    !     ! Subroutine arguments
    !     REAL(realk), INTENT(inout), TARGET :: ff(:)
    !     INTEGER(intk), INTENT(in) :: igridf
    !     INTEGER(intk), INTENT(in) :: pos

    !     ! Local variables
    !     INTEGER(intk) :: ip3
    !     INTEGER(intk) :: k, j, i, kc, jc, ic
    !     INTEGER(intk) :: kk, jj, ii, kkc, jjc, iic
    !     REAL(realk), POINTER :: fc(:, :, :)
    !     REAL(realk), POINTER :: fff(:, :, :)

    !     CALL get_ip3(ip3, igridf)
    !     CALL get_mgdims(kk, jj, ii, igridf)
    !     fff(1:kk, 1:jj, 1:ii) => ff(ip3:ip3 + kk*jj*ii - 1)

    !     ! We map the recvbuf to a 3-D field to make lookup easier
    !     ! (remember that only 2..kkc-1 are send - so kkc here is not
    !     ! really the same as kk for the coarse grid)
    !     kkc = kk/2
    !     jjc = jj/2
    !     iic = ii/2
    !     fc(1:kkc, 1:jjc, 1:iic) => recvbuf(pos:pos + kkc*jjc*iic - 1)

    !     DO i = 1, ii
    !         DO j = 1, jj
    !             DO k = 1, kk
    !                 ic = (i-1)/2 + 1
    !                 jc = (j-1)/2 + 1
    !                 kc = (k-1)/2 + 1
    !                 fff(k, j, i) = fc(kc, jc, ic)
    !             END DO
    !         END DO
    !     END DO
    ! END SUBROUTINE prolong_finish




    SUBROUTINE prepare_recvtask(rtasks, nrtasks)
        ! Subroutine arguments
        INTEGER(intk), INTENT(inout) :: rtasks(:, :)
        INTEGER(intk), INTENT(out) :: nrtasks

        ! Local variables
        INTEGER(int32) :: idx
        INTEGER(intk) :: irtasks

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

    END SUBROUTINE prepare_recvtask


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



    SUBROUTINE process_recvtask(rtask, nrtasks)
        ! Subroutine arguments
        INTEGER(intk), INTENT(in) :: rtask(2, :)
        INTEGER(intk), INTENT(in) :: nrtasks

        ! Local variables
        INTEGER(intk) :: igridf, idx_recvbuf, itask, ii, jj, kk, ip3

        ASSOCIATE(fine => ff%arr)

        DO itask = 1, nrtasks

            ! Unpacking the task
            igridf = rtask(1, itask)
            idx_recvbuf = rtask(2, itask)

            ! Getting parameters of fine grid
            CALL get_mgdims(kk, jj, ii, igridf)
            CALL get_ip3(ip3, igridf)

            ! Copying from buffer to the fine grid
            CALL write_fine(kk, jj, ii, fine(ip3), idx_recvbuf)

        END DO

        END ASSOCIATE

    END SUBROUTINE process_recvtask


    SUBROUTINE write_fine(kk, jj, ii, fff, fc)
        ! Subroutine arguments
        INTEGER(intk), INTENT(in) :: kk, jj, ii
        REAL(realk), INTENT(inout) :: fff(kk, jj, ii)
        REAL(realk), INTENT(in) :: fc(:)

        ! Local variables
        INTEGER(intk) :: k, j, i, kc, jc, ic
        INTEGER(intk) :: kkc, jjc, iic

        ! The counter for the buffer reflects the coarse frid
        kkc = kk/2
        jjc = jj/2
        iic = ii/2

        ! Loop for filling the whole fine grid
        DO i = 1, ii
            DO j = 1, jj
                DO k = 1, kk
                    ! Reconstructing the coarse grid indices
                    ic = (i-1)/2 + 1
                    jc = (j-1)/2 + 1
                    kc = (k-1)/2 + 1
                    ! Conversion into buffer index
                    idx = 1 + (kc-1) + (jc-1)*kkc + (ic-1)*kkc*jjc
                    fff(k, j, i) = recvbuf(idx)
                END DO
            END DO
        END DO
    END SUBROUTINE write_fine




    ! ! Finish prolongation, i.e. distribute the data on the grid
    ! SUBROUTINE prolong_finish(ff, igridf, pos)
    !     ! Subroutine arguments
    !     REAL(realk), INTENT(inout), TARGET :: ff(:)
    !     INTEGER(intk), INTENT(in) :: igridf
    !     INTEGER(intk), INTENT(in) :: pos

    !     ! Local variables
    !     INTEGER(intk) :: ip3
    !     INTEGER(intk) :: k, j, i, kc, jc, ic
    !     INTEGER(intk) :: kk, jj, ii, kkc, jjc, iic
    !     REAL(realk), POINTER :: fc(:, :, :)
    !     REAL(realk), POINTER :: fff(:, :, :)

    !     CALL get_ip3(ip3, igridf)
    !     CALL get_mgdims(kk, jj, ii, igridf)
    !     fff(1:kk, 1:jj, 1:ii) => ff(ip3:ip3 + kk*jj*ii - 1)

    !     ! We map the recvbuf to a 3-D field to make lookup easier
    !     ! (remember that only 2..kkc-1 are send - so kkc here is not
    !     ! really the same as kk for the coarse grid)
    !     kkc = kk/2
    !     jjc = jj/2
    !     iic = ii/2
    !     fc(1:kkc, 1:jjc, 1:iic) => recvbuf(pos:pos + kkc*jjc*iic - 1)

    !     DO i = 1, ii
    !         DO j = 1, jj
    !             DO k = 1, kk
    !                 ic = (i-1)/2 + 1
    !                 jc = (j-1)/2 + 1
    !                 kc = (k-1)/2 + 1
    !                 fff(k, j, i) = fc(kc, jc, ic)
    !             END DO
    !         END DO
    !     END DO
    ! END SUBROUTINE prolong_finish





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
        is_recording = .TRUE.
        DO ilevel = minlevel, maxlevel
            CALL ctof2(ilevel, dummy%arr, dummy%arr)
        END DO
        is_recording = .FALSE.

    END SUBROUTINE init_ctof2



    SUBROUTINE finish_ctof2()

        IF (is_init .NEQV. .TRUE.) THEN
            RETURN
        END IF

        is_init = .FALSE.

        DEALLOCATE(sendreqs)
        DEALLOCATE(recvreqs)
        DEALLOCATE(recvgrids)
        DEALLOCATE(recvpos)

    END SUBROUTINE finish_ctof2

END MODULE ctof_mod
