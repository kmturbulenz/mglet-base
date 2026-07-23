MODULE ctof1_mod

    USE core_mod
    USE MPI_f08

    IMPLICIT NONE (type, external)
    PRIVATE

    ! The information in the first dimension is sorted as follows:
    !   Field 1: Rank of sending process
    !   Field 2: Rank of receiving process
    !   Field 3: ID of sending grid (coarse grid)
    !   Field 4: ID of receiving grid (fine grid)
    INTEGER(intk), ALLOCATABLE :: sendconns(:, :), recvconns(:, :)

    ! Number of send and receive connections
    INTEGER(intk) :: isend = 0, irecv = 0

    ! Variable to indicate if the required data structures and MPI-types
    ! have been created
    LOGICAL :: is_init = .FALSE.

    ! Maximum allowed number of childs per parent (i.e. maximum number of
    ! send-conenctions per grid)
    INTEGER(intk), PARAMETER :: maxchilds = 8

    ! Lists that hold the send and receive request arrays
    TYPE(MPI_Request), ALLOCATABLE :: sendreqs(:), recvreqs(:)

    ! Actual number of messages that are to be sendt and received in one
    ! "round" of operations
    INTEGER(intk) :: nsend, nrecv
    INTEGER(int32), ALLOCATABLE :: recvlist(:)
    INTEGER(intk), ALLOCATABLE :: recvidxlist(:, :)

    ! ! List of grids to receive data on
    ! INTEGER(intk), ALLOCATABLE :: recvgrids(:), recvpos(:)

    PUBLIC :: ctof1, init_ctof1, finish_ctof1

CONTAINS

    ! Nota bene: Routine is thought to be called from the fine side, i.e. the
    ! receiving side of the prolongation. The child grid has parent information
    ! but the parent grid does not have information about its children.

    SUBROUTINE ctof1(ilevel, ff, fc)
        ! Subroutine arguments
        INTEGER(intk), INTENT(in) :: ilevel  ! Level of the *fine* side
        TYPE(field_t), INTENT(inout) :: ff
        TYPE(field_t), INTENT(in) :: fc

        CALL start_timer(231)
        IF (.NOT. is_init) THEN
            CALL errr(__FILE__, __LINE__)
        END IF

        ! Posting non-blocking receives and preparing lists for unpacking
        CALL recv_all(ilevel)

        ! Packing data into send buffer and posting non-blocking sends
        CALL send_all(ilevel, fc)

        ! Querying MPI for completed receives and unpacking the recv buffer
        CALL process_bufs(ff)

    END SUBROUTINE ctof1



    ! Perform all Recv-calls
    SUBROUTINE recv_all(ilevel)
        ! Subroutine arguments
        INTEGER(intk), INTENT(in) :: ilevel   ! Level of the *fine* side

        ! Local variables
        INTEGER(intk) :: kk, jj, ii, ncells, i, iprocnbr, igridf
        INTEGER(int32) :: recvcounter, messagelength

        ! Post all receive calls
        recvcounter = 0
        messagelength = 0
        nrecv = 0
        recvidxlist = -1
        recvlist = -1

        ! - 1: Rank of sending process
        ! - 2: Rank of receiving process
        ! - 3: ID of sending grid (coarse grid)
        ! - 4: ID of receiving grid (fine grid)

        ! Iteration over all receive connections
        DO i = 1, irecv

            ! Identifier of the fine grid involved in this connection
            igridf = recvconns(4, i)

            ! Only continue if the fine grid is on the level
            IF (ilevel == level(igridf)) THEN

                ! Getting the sender process and the message size for grid
                iprocnbr = recvconns(1, i)
                CALL get_mgdims(kk, jj, ii, igridf)
                ncells = kk*jj*ii/8

                ! Entering the information into the recvidxlist arrays
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
                    ! Posting at final receive connections
                    CALL post_recv(iprocnbr, messagelength, recvcounter)
                ELSE IF (recvconns(1, i + 1) /= iprocnbr) THEN
                    ! Posting at final receive connections from this sender
                    CALL post_recv(iprocnbr, messagelength, recvcounter)
                END IF
                ! Note that "post_recv" increments nrecv, resets
                ! messagelength, enters recvlist and increments recvcounter
            END IF
        END DO
    END SUBROUTINE recv_all


    SUBROUTINE post_recv(iprocnbr, messagelength, recvcounter)
        ! Subroutine arguments
        INTEGER(int32), INTENT(in) :: iprocnbr
        INTEGER(int32), INTENT(inout) :: messagelength
        INTEGER(int32), INTENT(inout) :: recvcounter

        ! Local variables
        ! none...

        nrecv = nrecv + 1
        recvlist(nrecv) = iprocnbr

        CALL MPI_Irecv(recvbuf(recvcounter+1), messagelength, &
            mglet_mpi_real, iprocnbr, 1, MPI_COMM_WORLD, recvreqs(nrecv))

        recvcounter = recvcounter + messagelength
        messagelength = 0

    END SUBROUTINE post_recv


    SUBROUTINE send_all(ilevel, fc)
        ! Subroutine arguments
        INTEGER(intk), INTENT(in) :: ilevel   ! Level of the *fine* grid
        TYPE(field_t), INTENT(in) :: fc

        ! Local variables
        INTEGER(intk) :: i, igridf, iprocnbr
        INTEGER(int32) :: sendcounter, messagelength

        ! Post all receive calls
        sendcounter = 0
        messagelength = 0
        nsend = 0

        ! - 1: Rank of sending process
        ! - 2: Rank of receiving process
        ! - 3: ID of sending grid (coarse grid)
        ! - 4: ID of receiving grid (fine grid)

        ! Iteration over all send connections
        DO i = 1, isend

            igridf = sendconns(4, i)
            iprocnbr = sendconns(2, i)

            IF (ilevel == level(igridf)) THEN
                ! Extract data from the coarse grid and compact in into buffer
                CALL write_buffer(i, messagelength, sendcounter, fc)
            END IF

            IF (messagelength > 0) THEN
                IF (i == isend) THEN
                    ! Posting at final send connection
                    CALL post_send(iprocnbr, messagelength, sendcounter)
                ELSE IF (sendconns(1, i + 1) /= iprocnbr) THEN
                    ! Posting at final send connection to this receiver
                    CALL post_send(iprocnbr, messagelength, sendcounter)
                END IF
            END IF

        END DO
    END SUBROUTINE send_all


    SUBROUTINE write_buffer(sendid, messagelength, sendcounter, fc)
        ! Subroutine arguments
        INTEGER(int32), INTENT(in) :: sendid
        INTEGER(int32), INTENT(inout) :: messagelength
        INTEGER(int32), INTENT(in) :: sendcounter
        TYPE(field_t), INTENT(in) :: fc

        ! Local variables
        REAL(realk), POINTER, CONTIGUOUS :: fcptr(:, :, :)
        INTEGER(intk) :: kkf, jjf, iif
        INTEGER(intk) :: ista, jsta, ksta, isto, jsto, ksto
        INTEGER(intk) :: igridc, igridf
        INTEGER(int32) :: thismessagelength, offset, counter

        ! - 1: Rank of sending process
        ! - 2: Rank of receiving process
        ! - 3: ID of sending grid (coarse grid)
        ! - 4: ID of receiving grid (fine grid)

        igridf = sendconns(4, sendid)
        igridc = sendconns(3, sendid)

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
        CALL fc%get_ptr(fcptr, igridc)

        CALL pack_single(sendbuf(offset:offset+thismessagelength-1), fcptr, &
            ista, isto, jsta, jsto, ksta, ksto, counter)

        offset = offset + thismessagelength

        IF (offset /= sendcounter + messagelength + thismessagelength + 1) THEN
            WRITE(*, *) "offset:", offset, "  expected:", sendcounter + &
                messagelength + thismessagelength + 1
            CALL errr(__FILE__, __LINE__)
        END IF

        messagelength = messagelength + thismessagelength
    END SUBROUTINE write_buffer


    SUBROUTINE pack_single(buf, fcptr, ista, isto, jsta, jsto, ksta, ksto, &
            counter)
        ! Subroutine arguments
        REAL(realk), INTENT(inout), CONTIGUOUS :: buf(:)
        REAL(realk), INTENT(in), CONTIGUOUS :: fcptr(:, :, :)
        INTEGER(intk), INTENT(in) :: ista, isto, jsta, jsto, ksta, ksto
        INTEGER(int32), INTENT(out) :: counter

        ! Local variables
        INTEGER(intk) :: i, j, k

        counter = 0
        DO i = ista, isto
            DO j = jsta, jsto
                DO k = ksta, ksto
                    counter = counter + 1
                    buf(counter) = fcptr(k, j, i)
                END DO
            END DO
        END DO
    END SUBROUTINE pack_single


    SUBROUTINE post_send(iprocnbr, messagelength, sendcounter)
        ! Subroutine arguments
        INTEGER(int32), INTENT(in) :: iprocnbr
        INTEGER(int32), INTENT(inout) :: messagelength
        INTEGER(int32), INTENT(inout) :: sendcounter

        ! Local variables
        ! none...

        nsend = nsend + 1

        CALL MPI_Isend(sendbuf(sendcounter + 1), messagelength, &
            mglet_mpi_real, iprocnbr, 1, MPI_COMM_WORLD, sendreqs(nsend))

        sendcounter = sendcounter + messagelength
        messagelength = 0

    END SUBROUTINE post_send


    SUBROUTINE process_bufs(ff)
        ! Subroutine arguments
        TYPE(field_t), INTENT(inout) :: ff

        ! Local variables
        TYPE(MPI_Status) :: recvstatus
        INTEGER(intk) :: i
        INTEGER(int32) :: idx, recvmessagelen, unpacklen

        ! - 1: Rank of sending process
        ! - 2: Rank of receiving process
        ! - 3: ID of sending grid (coarse grid)
        ! - 4: ID of receiving grid (fine grid)

        ! recvidxlist(1, i) = iprocnbr
        ! recvidxlist(2, i) = ncells
        ! recvidxlist(3, i) = recvcounter + messagelength
        ! messagelength = messagelength + ncells


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
                        CALL read_buffer(i, ff)
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

    END SUBROUTINE process_bufs


    SUBROUTINE read_buffer(recvid, ff)
        ! Subroutine arguments
        INTEGER(intk), INTENT(in) :: recvid
        TYPE(field_t), INTENT(inout) :: ff

        ! Local variables
        REAL(realk), POINTER, CONTIGUOUS :: ffptr(:, :, :)
        INTEGER(intk) :: igridf, kk, jj, ii, offset, thismessagelength

        igridf = recvconns(3, recvid)
        offset = recvidxlist(3, recvid) + 1

        CALL ff%get_ptr(ffptr, igridf)
        CALL get_mgdims(kk, jj, ii, igridf)
        thismessagelength = kk*jj*ii/8

        CALL unpack_single(recvbuf(offset:offset+thismessagelength-1), &
            ffptr, kk, jj, ii)

        IF (thismessagelength /= recvidxlist(2, recvid)) THEN
            CALL errr(__FILE__, __LINE__)
        END IF
    END SUBROUTINE read_buffer


    SUBROUTINE unpack_single(buf, ffptr, kk, jj, ii)
        ! Subroutine arguments
        REAL(realk), INTENT(in), CONTIGUOUS :: buf(:)
        REAL(realk), INTENT(inout), CONTIGUOUS :: ffptr(:, :, :)
        INTEGER(intk), INTENT(in) :: kk, jj, ii

        ! Local variables
        INTEGER(intk) :: k, j, i, kc, jc, ic
        INTEGER(intk) :: kkc, jjc

        kkc = kk/2
        jjc = jj/2

        DO i = 1, ii
            DO j = 1, jj
                DO k = 1, kk
                    ic = (i-1)/2 + 1
                    jc = (j-1)/2 + 1
                    kc = (k-1)/2 + 1
                    ffptr(k, j, i) = buf(kc + (jc-1)*kkc + (ic-1)*kkc*jjc)
                END DO
            END DO
        END DO
    END SUBROUTINE unpack_single


    ! SUBROUTINE pack_send(buf, kk, jj, ii, fc, igridc, igridf)
    !     ! Subroutine arguments
    !     REAL(realk), INTENT(inout) :: buf(:)
    !     INTEGER(intk), INTENT(in) :: kk, jj, ii
    !     REAL(realk), INTENT(in) :: fc(kk, jj, ii)
    !     INTEGER(intk), INTENT(in) :: igridc
    !     INTEGER(intk), INTENT(in) :: igridf

    !     ! Local variables
    !     INTEGER(intk) :: i, j, k
    !     INTEGER(intk) :: kkf, jjf, iif
    !     INTEGER(intk) :: counter
    !     INTEGER(intk) :: ista, jsta, ksta, isto, jsto, ksto

    !     ! Compute start- and end-positions in coarse grid
    !     ista = iposition(igridf) - 1
    !     jsta = jposition(igridf) - 1
    !     ksta = kposition(igridf) - 1

    !     CALL get_mgdims(kkf, jjf, iif, igridf)
    !     isto = ista + (iif - 4)/2 + 1
    !     jsto = jsta + (jjf - 4)/2 + 1
    !     ksto = ksta + (kkf - 4)/2 + 1

    !     ! Pack buffer
    !     counter = 0
    !     DO i = ista, isto
    !         DO j = jsta, jsto
    !             DO k = ksta, ksto
    !                 counter = counter + 1
    !                 buf(counter) = fc(k, j, i)
    !             END DO
    !         END DO
    !     END DO

    !     ! Sanity checks
    !     IF (counter /= kkf*jjf*iif/8) THEN
    !         WRITE(*, *) "counter = ", counter
    !         WRITE(*, *) "kkf = ", kkf
    !         WRITE(*, *) "jjf = ", jjf
    !         WRITE(*, *) "iif = ", iif
    !         WRITE(*, *) "ksta = ", ksta
    !         WRITE(*, *) "jsta = ", jsta
    !         WRITE(*, *) "ista = ", ista
    !         WRITE(*, *) "ksto = ", ksto
    !         WRITE(*, *) "jsto = ", jsto
    !         WRITE(*, *) "isto = ", isto
    !         CALL errr(__FILE__, __LINE__)
    !     END IF
    !     IF (counter /= SIZE(buf)) THEN
    !         CALL errr(__FILE__, __LINE__)
    !     END IF
    ! END SUBROUTINE pack_send


    ! ! Finish prolongation
    ! !
    ! ! Wait for communication to finish and clean up
    ! SUBROUTINE ctof_end(ff)
    !     ! Subroutine arguments
    !     REAL(realk), INTENT(inout) :: ff(:)

    !     ! Local variables
    !     INTEGER(int32) :: idx

    !     CALL start_timer(232)

    !     IF (.NOT. in_progress) THEN
    !         CALL errr(__FILE__, __LINE__)
    !     END IF

    !     IF (nrecv > 0) THEN
    !         DO WHILE (.TRUE.)
    !             CALL MPI_Waitany(nrecv, recvreqs, idx, MPI_STATUS_IGNORE)

    !             IF (idx /= MPI_UNDEFINED) THEN
    !                 CALL start_timer(235)
    !                 CALL prolong_finish(ff, recvgrids(idx), recvpos(idx))
    !                 CALL stop_timer(235)
    !             ELSE
    !                 EXIT
    !             END IF
    !         END DO
    !     END IF

    !     CALL MPI_Waitall(nsend, sendreqs, MPI_STATUSES_IGNORE)

    !     in_progress = .FALSE.

    !     CALL stop_timer(232)
    ! END SUBROUTINE ctof_end


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
    SUBROUTINE init_ctof1()

        ! Local variables
        INTEGER(intk) :: i, igrid, iprocc, ipar, maxconns
        INTEGER(int32), ALLOCATABLE :: sendcounts(:), sdispls(:)
        INTEGER(int32), ALLOCATABLE :: recvcounts(:), rdispls(:)
        INTEGER(intk), PARAMETER :: ncols = 4
        INTEGER(intk), ALLOCATABLE :: recvtmp(:, :), sendtmp(:, :)

        ! Avoiding repeated initialization
        IF (is_init) CALL errr(__FILE__, __LINE__)

        CALL set_timer(230, "CTOF")
        CALL set_timer(231, "CTOF_BEGIN")
        CALL set_timer(232, "CTOF_END")
        CALL set_timer(235, "CTOF_PROLONG_FINISH")

        ! Setting up the infrastructure for bundling of messages
        maxconns = nmygrids
        ALLOCATE(recvconns(ncols, maxconns))

        ALLOCATE(sendcounts(0:numprocs-1), SOURCE=0)
        ALLOCATE(sdispls(0:numprocs-1), SOURCE=0)
        ALLOCATE(recvcounts(0:numprocs-1), SOURCE=0)
        ALLOCATE(rdispls(0:numprocs-1), SOURCE=0)

        ! Filling the recvconns array by looking up parents
        nrecv = 0
        DO i = 1, nmygrids

            ! Obtaining the grid ID of the fine grid and its parent
            igrid = mygrids(i)
            ipar = iparent(igrid)
            IF (ipar == 0) CYCLE

            ! Obtaining the process ID of the coarse parent grid
            iprocc = idprocofgrd(ipar)

            ! Adding a receive connection
            nrecv = nrecv + 1
            IF (nrecv > maxconns) THEN
                CALL errr(__FILE__, __LINE__)
            END IF

            ! - 1: Rank of sending process
            ! - 2: Rank of receiving process
            ! - 3: ID of sending grid (coarse grid)
            ! - 4: ID of receiving grid (fine grid)

            recvconns(1, nrecv) = iprocc
            recvconns(2, nrecv) = myid
            recvconns(3, nrecv) = ipar
            recvconns(4, nrecv) = igrid

            ! Storing information used later in the context of AllToAll
            recvcounts(iprocc) = recvcounts(iprocc) + ncols
        END DO
        irecv = nrecv

        ! Sort recvconns by process ID (col 1) = sender with the coarse grid
        CALL sort_conns(recvconns(:, 1:nrecv), 1)

        ! Calculate sdispl offset (used in MPI_Alltoallv)
        DO i = 1, numprocs-1
            sdispls(i) = sdispls(i-1) + recvcounts(i-1)
        END DO

        ! In a first step, the NUMBER OF ELEMENTS TO SEND must be computed
        ! by exchanging information about the number of received elements
        CALL MPI_Alltoall(recvcounts, 1, MPI_INTEGER, sendcounts, 1, &
            MPI_INTEGER, MPI_COMM_WORLD)

        ! Array sendcounts is now filled and offsets can be computed
        DO i = 1, numprocs-1
            rdispls(i) = rdispls(i-1) + sendcounts(i-1)
        END DO

        ! The total number of send connections is the sum of all sendcounts
        isend = (rdispls(numprocs-1) + sendcounts(numprocs-1)) / ncols
        ALLOCATE(sendconns(ncols, isend))
        sendconns = 0

        ! Exchange connection information
        CALL MPI_Alltoallv(recvconns, recvcounts, sdispls, MPI_INTEGER, &
            sendconns, sendcounts, rdispls, MPI_INTEGER, MPI_COMM_WORLD)

        ! Both recvconns and sendconns should be fully populated and ordered
        DO i = 1, isend-1
            IF (sendconns(2, i) > sendconns(2, i+1)) THEN
                WRITE(*, *) "Sendconns not sorted by receiving rank"
                CALL errr(__FILE__, __LINE__)
            END IF
        END DO

        DO i = 1, irecv-1
            IF (recvconns(1, i) > recvconns(1, i+1)) THEN
                WRITE(*, *) "Recvconns not sorted by sending rank"
                CALL errr(__FILE__, __LINE__)
            END IF
        END DO

        ! Deallocate arrays which were only needed to operate MPI_Alltoallv
        DEALLOCATE(rdispls)
        DEALLOCATE(recvcounts)
        DEALLOCATE(sdispls)
        DEALLOCATE(sendcounts)

        ! That should now be isend and irecv (--> TO DO)
        ALLOCATE(sendreqs(isend))
        ALLOCATE(recvreqs(irecv))
        ALLOCATE(recvlist(irecv))
        ALLOCATE(recvidxlist(3, irecv))

        ! Reallocating sendconns and recvconns to the actual size
        ALLOCATE(sendtmp(ncols, isend), SOURCE=sendconns(:, 1:isend))
        CALL move_alloc(sendtmp, sendconns)
        ALLOCATE(recvtmp(ncols, irecv), SOURCE=recvconns(:, 1:irecv))
        CALL move_alloc(recvtmp, recvconns)

        nrecv = 0
        nsend = 0
        is_init = .TRUE.

    END SUBROUTINE init_ctof1


    SUBROUTINE finish_ctof1()
        IF (.NOT. is_init) THEN
            RETURN
        END IF

        is_init = .FALSE.

        ! Deallocation of infrastructure arrays
        DEALLOCATE(sendconns)
        DEALLOCATE(recvconns)
        DEALLOCATE(sendreqs)
        DEALLOCATE(recvreqs)
        DEALLOCATE(recvlist)
        DEALLOCATE(recvidxlist)

    END SUBROUTINE finish_ctof1

END MODULE ctof1_mod
