
MODULE commbuf_test_mod

    USE precision_mod
    USE commbuf_mod, ONLY: sendbuf, recvbuf, device_sendbuf, device_recvbuf, &
        idim_mg_bufs
    USE comms_mod, ONLY: myid, numprocs
    USE err_mod, ONLY: errr
    USE MPI_f08

    IMPLICIT NONE (type, external)
    PRIVATE

    PUBLIC :: test_commbuf

CONTAINS

    SUBROUTINE test_commbuf()

        INTEGER(intk) :: i, n, ierr, peer_send, peer_recv, mpierr, size
        TYPE(MPI_Request) :: reqs(2)
        TYPE(MPI_Status) :: statuses(2)
        REAL(realk), ALLOCATABLE :: field(:)
        REAL(realk) :: exptd

        ALLOCATE(field(idim_mg_bufs))
        field = 0.0_realk
        !$omp target enter data map(to: field)

        ! Setting up round-robin communication pattern for testing
        peer_send = MOD(myid + 1, numprocs)
        peer_recv = MOD(myid - 1 + numprocs, numprocs)
        size = INT(idim_mg_bufs, intk)

        !$omp target teams distribute parallel do private(i)
        DO i = 1, idim_mg_bufs
            sendbuf(i) = REAL(i, realk) + REAL(myid, realk) * 100.0_realk
        END DO
        !$omp end target teams distribute parallel do

        CALL MPI_Irecv(device_recvbuf, size, mglet_mpi_real, &
            peer_recv, 42, MPI_COMM_WORLD, reqs(1), mpierr)
        IF (mpierr /= MPI_SUCCESS) CALL MPI_Abort(MPI_COMM_WORLD, 3, mpierr)

        CALL MPI_Isend(device_sendbuf, size, mglet_mpi_real, &
            peer_send, 42, MPI_COMM_WORLD, reqs(2), mpierr)
        IF (mpierr /= MPI_SUCCESS) CALL MPI_Abort(MPI_COMM_WORLD, 3, mpierr)

        CALL MPI_Waitall(2, reqs, statuses, mpierr)
        IF (mpierr /= MPI_SUCCESS) CALL MPI_Abort(MPI_COMM_WORLD, 3, mpierr)

        !$omp target teams distribute parallel do private(i)
        DO i = 1, idim_mg_bufs
            field(i) = recvbuf(i)
        END DO
        !$omp end target teams distribute parallel do

        !$omp target update from(field)

        DO i = 1, idim_mg_bufs
            exptd = REAL(i, realk) + REAL(peer_recv, realk) * 100.0_realk
            IF (ABS(field(i) - exptd) > 0.000001) THEN
                WRITE(*, *) "Commbufs failing"
                CALL errr(__FILE__, __LINE__)
            END IF
        END DO

        WRITE(*, *) "Commbufs test passed on rank ", myid

        !$omp target exit data map(release: field)
        DEALLOCATE(field)

    END SUBROUTINE test_commbuf

END MODULE commbuf_test_mod