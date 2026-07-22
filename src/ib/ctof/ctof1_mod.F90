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
    LOGICAL :: isinit = .FALSE.

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

    PUBLIC :: ctof1, init_ctof1, finish_ctof1

CONTAINS

    ! Nota bene: Routine is thought to be called from the fine side, i.e. the
    ! receiving side of the prolongation. The child grid has parent information
    ! but the parent grid does not have information about its children.

    SUBROUTINE ctof1(ilevel, ff, fc)
        ! Subroutine arguments
        INTEGER(intk), INTENT(in) :: ilevel  ! Level of the *fine* side
        REAL(realk), INTENT(inout) :: ff(:)
        REAL(realk), INTENT(in) :: fc(:)

        CALL ctof_begin(ilevel, fc)
        CALL ctof_end(ff)

    END SUBROUTINE ctof1


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

        IF (.NOT. isInit) THEN
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


    ! Finish prolongation, i.e. distribute the data on the grid
    SUBROUTINE prolong_finish(ff, igridf, pos)
        ! Subroutine arguments
        REAL(realk), INTENT(inout), TARGET :: ff(:)
        INTEGER(intk), INTENT(in) :: igridf
        INTEGER(intk), INTENT(in) :: pos

        ! Local variables
        INTEGER(intk) :: ip3
        INTEGER(intk) :: k, j, i, kc, jc, ic
        INTEGER(intk) :: kk, jj, ii, kkc, jjc, iic
        REAL(realk), POINTER :: fc(:, :, :)
        REAL(realk), POINTER :: fff(:, :, :)

        CALL get_ip3(ip3, igridf)
        CALL get_mgdims(kk, jj, ii, igridf)
        fff(1:kk, 1:jj, 1:ii) => ff(ip3:ip3 + kk*jj*ii - 1)

        ! We map the recvbuf to a 3-D field to make lookup easier
        ! (remember that only 2..kkc-1 are send - so kkc here is not
        ! really the same as kk for the coarse grid)
        kkc = kk/2
        jjc = jj/2
        iic = ii/2
        fc(1:kkc, 1:jjc, 1:iic) => recvbuf(pos:pos + kkc*jjc*iic - 1)

        DO i = 1, ii
            DO j = 1, jj
                DO k = 1, kk
                    ic = (i-1)/2 + 1
                    jc = (j-1)/2 + 1
                    kc = (k-1)/2 + 1
                    fff(k, j, i) = fc(kc, jc, ic)
                END DO
            END DO
        END DO
    END SUBROUTINE prolong_finish


    ! Initialize arrays and data types
    SUBROUTINE init_ctof1()

        ! Local variables
        INTEGER(intk) :: i, igrid, iprocc, ipar, maxconns
        INTEGER(int32), ALLOCATABLE :: sendcounts(:), sdispls(:)
        INTEGER(int32), ALLOCATABLE :: recvcounts(:), rdispls(:)
        INTEGER(intk), PARAMETER :: ncols = 4

        ! Avoiding repeated initialization
        IF (isinit) CALL errr(__FILE__, __LINE__)

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
        ALLOCATE(sendreqs(nmygrids*maxchilds))
        ALLOCATE(recvreqs(nmygrids))
        ALLOCATE(recvgrids(nmygrids))
        ALLOCATE(recvpos(nmygrids))

        nrecv = 0
        nsend = 0

        isinit = .TRUE.
        in_progress = .FALSE.

    END SUBROUTINE init_ctof1


    SUBROUTINE finish_ctof1()
        IF (.NOT. isinit) THEN
            RETURN
        END IF

        isinit = .FALSE.

        DEALLOCATE(sendreqs)
        DEALLOCATE(recvreqs)
        DEALLOCATE(recvgrids)
        DEALLOCATE(recvpos)
    END SUBROUTINE finish_ctof1
END MODULE ctof1_mod
