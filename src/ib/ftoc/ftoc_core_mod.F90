MODULE ftoc_core_mod

    ! Module provides the shared communication infrastructure for the
    ! implementations of FTOC (Fine to Coarse) operations

    USE precision_mod
    USE MPI_f08

    IMPLICIT NONE (type, external)
    PRIVATE

    ! The information in the first dimension is sorted as follows:
    !   Field 1: Rank of sending process
    !   Field 2: Rank of receiving process
    !   Field 3: ID of sending grid (fine grid)
    !   Field 4: ID of receiving grid (coarse grid)
    INTEGER(intk), ALLOCATABLE, PROTECTED :: sendconns(:, :), recvconns(:, :)

    ! Number of send and receive connections
    INTEGER(intk), PROTECTED :: isend = 0, irecv = 0

    ! Variable to indicate if the required data structures have been created
    LOGICAL, PROTECTED :: is_ctof_core_init = .FALSE.


    ! contained functions
    PUBLIC :: init_ftoc_core, finish_ftoc_core, sendconns, recvconns, &
        isend, irecv, is_ctof_core_init


CONTAINS

    SUBROUTINE init_ftoc_core()
        ! Local variables
        INTEGER(intk) :: i, igrid, iprocc, ipar, maxconns
        INTEGER(int32), ALLOCATABLE :: sendcounts(:), sdispls(:)
        INTEGER(int32), ALLOCATABLE :: recvcounts(:), rdispls(:)
        INTEGER(intk), ALLOCATABLE :: sendtmp(:, :)
        INTEGER(intk), PARAMETER :: ncols = 4

        ! Check if the module has already been initialized
        IF (is_ctof_core_init) CALL errr(__FILE__, __LINE__)

        ! Initialize the send connections to maximal size (shorten later)
        maxconns = nmygrids*8
        ALLOCATE(sendconns(ncols, maxconns))

        ! Initialize the send and receive counts and displacements arrays.
        ! These arrays are used in the MPI_Alltoallv call.
        ALLOCATE(sendcounts(0:numprocs-1), SOURCE=0)
        ALLOCATE(sdispls(0:numprocs-1), SOURCE=0)
        ALLOCATE(recvcounts(0:numprocs-1), SOURCE=0)
        ALLOCATE(rdispls(0:numprocs-1), SOURCE=0)

        isend = 0
        DO i = 1, nmygrids
            igrid = mygrids(i)
            ipar = iparent(igrid)
            IF (ipar == 0) CYCLE

            iprocc = idprocofgrd(ipar)

            isend = isend + 1
            IF (isend > maxconns) THEN
                CALL errr(__FILE__, __LINE__)
            END IF

            sendconns(1, isend) = myid
            sendconns(2, isend) = iprocc
            sendconns(3, isend) = igrid
            sendconns(4, isend) = ipar

            sendcounts(iprocc) = sendcounts(iprocc) + ncols
        END DO
        ! Thus, isend is already correctly set

        ! Sort sendconns by process ID (col 2)
        CALL sort_conns(sendconns(:, 1:isend), 2)

        ! Calculate sdispl offset
        DO i = 1, numprocs-1
            sdispls(i) = sdispls(i-1) + sendcounts(i-1)
        END DO

        ! First exchange NUMBER OF ELEMENTS TO RECEIVE, to be able to
        ! calculate rdispls array
        CALL MPI_Alltoall(sendcounts, 1, MPI_INTEGER, recvcounts, 1, &
            MPI_INTEGER, MPI_COMM_WORLD)

        ! Calculate rdispl offset
        DO i = 1, numprocs-1
            rdispls(i) = rdispls(i-1) + recvcounts(i-1)
        END DO

        ! Compute the number of receive connections
        irecv = (rdispls(numprocs-1) + recvcounts(numprocs-1))/ncols
        ! Allocation of recvconns in correct size
        ALLOCATE(recvconns(ncols, irecv))
        recvconns = 0

        ! Exchange connection information
        CALL MPI_Alltoallv(sendconns, sendcounts, sdispls, MPI_INTEGER, &
            recvconns, recvcounts, rdispls, MPI_INTEGER, &
            MPI_COMM_WORLD)

        ! Deallocate temporary arrays needed to operate MPI_Alltoallv
        DEALLOCATE(rdispls)
        DEALLOCATE(recvcounts)
        DEALLOCATE(sdispls)
        DEALLOCATE(sendcounts)

        ! While recvconn has exact extend, sendconn has been allocated to max
        ALLOCATE(sendtmp(ncols, isend), SOURCE=sendconns(:, 1:isend))
        CALL move_alloc(sendtmp, sendconns)

        is_ctof_core_init = .TRUE.

    END SUBROUTINE init_ftoc_core


    SUBROUTINE finish_ftoc_core()
        ! Function to deallocate arrays.
        ! After successful execution, the module varaible "is_init" is set
        ! to false.
        IF (.NOT. is_ctof_core_init) RETURN

        is_ctof_core_init = .FALSE.
        isend = 0
        irecv = 0
        DEALLOCATE(sendconns)
        DEALLOCATE(recvconns)

    END SUBROUTINE finish_ftoc_core

END MODULE ftoc_core_mod
