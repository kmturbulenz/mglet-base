MODULE ctof_core_mod

    USE core_mod
    USE MPI_f08

    IMPLICIT NONE (type, external)
    PRIVATE

    ! The information in the first dimension is sorted as follows:
    ! - field 1: Rank of receiving process
    ! - field 2: Rank of sending process
    ! - field 3: ID of receiving grid (fine grid)
    ! - field 4: ID of sending grid (coarse grid)
    INTEGER(intk), ALLOCATABLE, PROTECTED :: sendconns(:, :), recvconns(:, :)

    ! Number of send and receive connections
    INTEGER(intk), PROTECTED :: isend = 0, irecv = 0

    ! Variable to indicate if the required data structures have been created
    LOGICAL, PROTECTED :: has_infrastructure = .FALSE.

    PUBLIC :: init_core_ctof, finish_core_ctof, sendconns, recvconns, &
        isend, irecv, has_infrastructure

CONTAINS

    ! Initialize the communication infrastructure for coarse-to-fine
    ! communication with bundling of all transfer between two processes.
    !
    SUBROUTINE init_core_ctof()
        ! Local variables
        INTEGER :: i, igrid, iprocc, ipar, maxconns
        INTEGER(int32), ALLOCATABLE :: sendcounts(:), sdispls(:)
        INTEGER(int32), ALLOCATABLE :: recvcounts(:), rdispls(:)
        INTEGER(intk), PARAMETER :: ncols = 4
        INTEGER(intk), ALLOCATABLE :: recvtmp(:, :)

        CALL set_timer(230, "CTOF")

        IF (has_infrastructure) THEN
            WRITE(*, *) "CTOF infrastructure already initialized"
            CALL errr(__FILE__, __LINE__)
        END IF

        ! Fine grids only need information from one coarse grid
        maxconns = nmygrids
        ALLOCATE(recvconns(ncols, maxconns))

        ALLOCATE(recvcounts(0:numprocs-1), source=0_int32)
        ALLOCATE(rdispls(0:numprocs-1), source=0_int32)
        ALLOCATE(sendcounts(0:numprocs-1), source=0_int32)
        ALLOCATE(sdispls(0:numprocs-1), source=0_int32)

        irecv = 0
        DO i = 1, nmygrids

            ! Obtaining the grid ID of the fine grid and its parent
            igrid = mygrids(i)
            ipar = iparent(igrid)
            IF (ipar == 0) CYCLE

            ! Obtaining the process ID of the coarse parent grid
            iprocc = idprocofgrd(ipar)

            ! Adding a receive connection
            irecv = irecv + 1
            IF (irecv > maxconns) THEN
                CALL errr(__FILE__, __LINE__)
            END IF

            recvconns(1, irecv) = myid
            recvconns(2, irecv) = iprocc
            recvconns(3, irecv) = igrid
            recvconns(4, irecv) = ipar

            ! Storing information used later in the context of AllToAll
            recvcounts(iprocc) = recvcounts(iprocc) + ncols
        END DO
        ! irecv is now already correctly set

        ! Sort recvconns by process ID (col 2)
        CALL sort_conns(recvconns(:, 1:irecv), 2)

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

        ! >>> in the actual implementation, the following arrays are not used, but they are
        ! ALLOCATE(sendreqs(isend))
        ! ALLOCATE(recvreqs(irecv))
        ! ALLOCATE(recvlist(irecv))
        ! ALLOCATE(recvidxlist(3, irecv))

        ! Reallocating recvconns to the actual size
        ALLOCATE(recvtmp(ncols, irecv), SOURCE=recvconns(:, 1:irecv))
        CALL move_alloc(recvtmp, recvconns)

        has_infrastructure = .TRUE.

    END SUBROUTINE init_core_ctof


    ! Clearing the allocated communicaiton infrastructure.
    !
    SUBROUTINE finish_core_ctof()

        IF (.NOT. has_infrastructure) RETURN

        has_infrastructure = .FALSE.
        isend = -1
        irecv = -1

        ! Deallocation of infrastructure arrays
        DEALLOCATE(sendconns)
        DEALLOCATE(recvconns)
        DEALLOCATE(sendreqs)
        DEALLOCATE(recvreqs)
        DEALLOCATE(recvlist)
        DEALLOCATE(recvidxlist)

    END SUBROUTINE finish_core_ctof

END MODULE ctof_core_mod
