MODULE ftoc_core_mod
    USE core_mod
    USE MPI_f08
    USE restrict_mod, ONLY: message_length

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
    LOGICAL, PROTECTED :: is_ftoc_core_init = .FALSE.

    ! Public read-only variables
    PUBLIC :: sendconns, recvconns, isend, irecv, is_ftoc_core_init

    ! Public subroutines
    PUBLIC :: init_ftoc_core, finish_ftoc_core, count_ncells
CONTAINS
    SUBROUTINE init_ftoc_core()
        ! Subroutine arguments
        ! none...

        ! Local variables
        INTEGER(intk) :: i, igrid, iprocc, ipar, maxconns, nsend
        INTEGER(int32), ALLOCATABLE :: sendcounts(:), sdispls(:)
        INTEGER(int32), ALLOCATABLE :: recvcounts(:), rdispls(:)
        INTEGER(intk), PARAMETER :: ncols = 4

        IF (is_ftoc_core_init) CALL errr(__FILE__, __LINE__)

        CALL set_timer(220, "FTOC")

        maxconns = nmygrids*8
        ALLOCATE(sendconns(ncols, maxconns))

        ALLOCATE(sendcounts(0:numprocs-1), source=0)
        ALLOCATE(sdispls(0:numprocs-1), source=0)
        ALLOCATE(recvcounts(0:numprocs-1), source=0)
        ALLOCATE(rdispls(0:numprocs-1), source=0)

        nsend = 0
        DO i = 1, nmygrids
            igrid = mygrids(i)
            ipar = iparent(igrid)
            IF (ipar == 0) CYCLE

            iprocc = idprocofgrd(ipar)

            nsend = nsend + 1
            IF (nsend > maxconns) THEN
                CALL errr(__FILE__, __LINE__)
            END IF

            sendconns(1, nsend) = myid
            sendconns(2, nsend) = iprocc
            sendconns(3, nsend) = igrid
            sendconns(4, nsend) = ipar

            sendcounts(iprocc) = sendcounts(iprocc) + ncols
        END DO
        isend = nsend

        ! Sort sendconns by process ID (col 2)
        CALL sort_conns(sendconns(:, 1:nsend), 2)

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

        irecv = (rdispls(numprocs-1) + recvcounts(numprocs-1))/ncols
        ALLOCATE(recvconns(ncols, irecv))
        recvconns = 0

        ! Exchange connection information
        CALL MPI_Alltoallv(sendconns, sendcounts, sdispls, MPI_INTEGER, &
            recvconns, recvcounts, rdispls, MPI_INTEGER, &
            MPI_COMM_WORLD)

        DEALLOCATE(rdispls)
        DEALLOCATE(recvcounts)
        DEALLOCATE(sdispls)
        DEALLOCATE(sendcounts)

        is_ftoc_core_init = .TRUE.
        nsend = 0
    END SUBROUTINE init_ftoc_core


    SUBROUTINE count_ncells(ncells, flag, igrid, v1, v2, v3, s1, s2, s3)
        ! Subroutine arguments
        INTEGER(int32), INTENT(out) :: ncells
        CHARACTER(len=1), INTENT(in) :: flag
        INTEGER(intk), INTENT(in) :: igrid
        TYPE(field_t), INTENT(in), OPTIONAL :: v1, v2, v3, s1, s2, s3

        ! Local variables
        INTEGER(int32) :: messagelength

        ncells = 0

        IF (flag == '*') THEN
            IF (PRESENT(v1)) THEN
                CALL message_length(messagelength, 'U', igrid)
                ncells = ncells + messagelength
            END IF
            IF (PRESENT(v2)) THEN
                CALL message_length(messagelength, 'V', igrid)
                ncells = ncells + messagelength
            END IF
            IF (PRESENT(v3)) THEN
                CALL message_length(messagelength, 'W', igrid)
                ncells = ncells + messagelength
            END IF
            IF (PRESENT(s1)) THEN
                CALL message_length(messagelength, 'P', igrid)
                ncells = ncells + messagelength
            END IF
            IF (PRESENT(s2)) THEN
                CALL message_length(messagelength, 'P', igrid)
                ncells = ncells + messagelength
            END IF
            IF (PRESENT(s3)) THEN
                CALL message_length(messagelength, 'P', igrid)
                ncells = ncells + messagelength
            END IF
        ELSE
            CALL message_length(messagelength, flag, igrid)
            ncells = ncells + messagelength
        END IF
    END SUBROUTINE count_ncells


    SUBROUTINE finish_ftoc_core()
        ! Subroutine arguments
        ! none...

        ! Local variables
        ! none...

        IF (.NOT. is_ftoc_core_init) CALL errr(__FILE__, __LINE__)

        DEALLOCATE(sendconns)
        DEALLOCATE(recvconns)
        isend = 0
        irecv = 0

        is_ftoc_core_init = .FALSE.
    END SUBROUTINE finish_ftoc_core
END MODULE ftoc_core_mod
