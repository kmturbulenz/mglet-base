MODULE ctof_mod

    USE core_mod
    USE MPI_f08

    IMPLICIT NONE (type, external)
    PRIVATE

    ! The information in the first dimension is sorted as follows:
    !   Field 1: Rank of receiving process
    !   Field 2: Rank of sending process
    !   Field 3: ID of receiving grid (fine grid)
    !   Field 4: ID of sending grid (coarse grid)
    INTEGER(intk), ALLOCATABLE :: sendconns(:, :), recvconns(:, :)

    ! Lists that hold the send and receive request arrays
    TYPE(MPI_Request), ALLOCATABLE :: sendreqs(:), recvreqs(:)

    ! Lists that hold the messages that are ACTUALLY sent and received
    INTEGER(intk) :: nsend, nrecv
    INTEGER(int32), ALLOCATABLE :: recvlist(:)
    INTEGER(intk), ALLOCATABLE :: recvidxlist(:, :)

    ! Number of send and receive connections
    INTEGER(intk) :: isend = 0, irecv = 0

    ! Variable to indicate if the required data structures have been created
    LOGICAL :: is_init = .FALSE.

    PUBLIC :: ctof, init_ctof, finish_ctof

CONTAINS

    SUBROUTINE ctof(ilevel, ff, fc)
        ! Subroutine arguments
        INTEGER(intk), INTENT(in) :: ilevel  ! Level of the *fine* side
        TYPE(field_t), INTENT(inout) :: ff
        TYPE(field_t), INTENT(in) :: fc

        ! Local variables
        ! none...

        CALL start_timer(230)

        IF (.NOT. is_init) CALL errr(__FILE__, __LINE__)

        ! Posting non-blocking receives and preparing lists for unpacking
        CALL recv_all(ilevel)

        ! Packing data into send buffer and posting non-blocking sends
        CALL send_all(ilevel, fc)

        ! Querying MPI for completed receives and unpacking the recv buffer
        CALL process_bufs(ff)

        CALL stop_timer(230)
    END SUBROUTINE ctof


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
        recvlist = 0

        ! Iteration over all receive connections
        DO i = 1, irecv
            igridf = recvconns(3, i)
            IF (ilevel == level(igridf)) THEN
                iprocnbr = recvconns(2, i) ! The sender process (coarse side)

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
                    ! Posting at final receive connection
                    CALL post_recv(iprocnbr, messagelength, recvcounter)
                ELSE IF (recvconns(2, i + 1) /= iprocnbr) THEN
                    ! Posting at final receive connections from current sender
                    CALL post_recv(iprocnbr, messagelength, recvcounter)
                END IF
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

        ! Iteration over all send connections
        DO i = 1, isend
            igridf = sendconns(3, i)
            iprocnbr = sendconns(1, i)
            IF (ilevel == level(igridf)) THEN
                ! Extract data from the coarse grid and compact it into buffer
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
        CALL fc%get_ptr(fcptr, igridc)
        CALL pack_single(sendbuf(offset:offset+thismessagelength-1), fcptr, &
            ista, isto, jsta, jsto, ksta, ksto, counter)

        IF (counter /= thismessagelength) THEN
            WRITE(*, *) "counter:", counter, "expected:", thismessagelength
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


    SUBROUTINE process_bufs(ff)
        ! Subroutine arguments
        TYPE(field_t), INTENT(inout) :: ff

        ! Local variables
        TYPE(MPI_Status) :: recvstatus
        INTEGER(intk) :: i
        INTEGER(int32) :: idx, recvmessagelen, unpacklen

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
        INTEGER(intk) :: kkc, jjc, k, j, i, kc, jc, ic, idx

        kkc = kk/2
        jjc = jj/2

        DO i = 1, ii
            DO j = 1, jj
                DO k = 1, kk
                    ic = (i-1)/2 + 1
                    jc = (j-1)/2 + 1
                    kc = (k-1)/2 + 1
                    idx = kc + (jc-1)*kkc + (ic-1)*kkc*jjc
                    ffptr(k, j, i) = buf(idx)
                END DO
            END DO
        END DO
    END SUBROUTINE unpack_single


    ! Initialize arrays and data types
    SUBROUTINE init_ctof()
        ! Local variables
        INTEGER :: i, igrid, iprocc, ipar, maxconns
        INTEGER(int32), ALLOCATABLE :: sendcounts(:), sdispls(:)
        INTEGER(int32), ALLOCATABLE :: recvcounts(:), rdispls(:)
        INTEGER(intk), PARAMETER :: ncols = 4
        INTEGER(intk), ALLOCATABLE :: recvtmp(:, :)

        CALL set_timer(230, "CTOF")

        IF (is_init) CALL errr(__FILE__, __LINE__)

        ! Fine grids only need information from one coarse grid
        maxconns = nmygrids
        ALLOCATE(recvconns(ncols, maxconns))

        ALLOCATE(recvcounts(0:numprocs-1), source=0_int32)
        ALLOCATE(rdispls(0:numprocs-1), source=0_int32)
        ALLOCATE(sendcounts(0:numprocs-1), source=0_int32)
        ALLOCATE(sdispls(0:numprocs-1), source=0_int32)

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

            recvconns(1, nrecv) = myid
            recvconns(2, nrecv) = iprocc
            recvconns(3, nrecv) = igrid
            recvconns(4, nrecv) = ipar

            ! Storing information used later in the context of AllToAll
            recvcounts(iprocc) = recvcounts(iprocc) + ncols
        END DO
        irecv = nrecv

        ! Sort recvconns by process ID (col 2)
        CALL sort_conns(recvconns(:, 1:nrecv), 2)

        ! Calculate rdispl offset
        DO i = 1, numprocs-1
            rdispls(i) = rdispls(i-1) + recvcounts(i-1)
        END DO

        ! First exchange NUMBER OF ELEMENTS TO SEND, to be able to
        ! calculate sdispls array
        CALL MPI_Alltoall(recvcounts, 1, MPI_INTEGER, sendcounts, 1, &
            MPI_INTEGER, MPI_COMM_WORLD)

        ! Array sendcounts is now filled and offsets can be computed
        DO i = 1, numprocs-1
            sdispls(i) = sdispls(i-1) + sendcounts(i-1)
        END DO

        isend = (sdispls(numprocs-1) + sendcounts(numprocs-1))/ncols
        ALLOCATE(sendconns(ncols, isend))
        sendconns = 0

        ! Exchange connection information
        CALL MPI_Alltoallv(recvconns, recvcounts, rdispls, MPI_INTEGER, &
            sendconns, sendcounts, sdispls, MPI_INTEGER, MPI_COMM_WORLD)

        ! Both recvconns and sendconns should be fully populated and ordered
        DO i = 1, isend-1
            IF (sendconns(1, i) > sendconns(1, i+1)) THEN
                WRITE(*, *) "Sendconns not sorted by receiving rank"
                CALL errr(__FILE__, __LINE__)
            END IF
        END DO

        DO i = 1, irecv-1
            IF (recvconns(2, i) > recvconns(2, i+1)) THEN
                WRITE(*, *) "Recvconns not sorted by sending rank"
                CALL errr(__FILE__, __LINE__)
            END IF
        END DO

        ! Deallocate arrays which were only needed to operate MPI_Alltoallv
        DEALLOCATE(rdispls)
        DEALLOCATE(recvcounts)
        DEALLOCATE(sdispls)
        DEALLOCATE(sendcounts)

        ALLOCATE(sendreqs(isend))
        ALLOCATE(recvreqs(irecv))
        ALLOCATE(recvlist(irecv))
        ALLOCATE(recvidxlist(3, irecv))

        ! Reallocating recvconns to the actual size
        ALLOCATE(recvtmp(ncols, irecv), SOURCE=recvconns(:, 1:irecv))
        CALL move_alloc(recvtmp, recvconns)

        nrecv = 0
        nsend = 0
        is_init = .TRUE.

    END SUBROUTINE init_ctof


    SUBROUTINE finish_ctof()
        IF (.NOT. is_init) RETURN

        is_init = .FALSE.
        isend = 0
        irecv = 0

        ! Deallocation of infrastructure arrays
        DEALLOCATE(sendconns)
        DEALLOCATE(recvconns)
        DEALLOCATE(sendreqs)
        DEALLOCATE(recvreqs)
        DEALLOCATE(recvlist)
        DEALLOCATE(recvidxlist)

    END SUBROUTINE finish_ctof

END MODULE ctof_mod
