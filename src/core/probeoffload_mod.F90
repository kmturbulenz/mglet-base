MODULE probeoffload_mod
    USE comms_mod, ONLY: myid, numprocs
    USE precision_mod, ONLY: intk

    IMPLICIT NONE(type, external)
    PRIVATE

    PUBLIC :: init_offload
CONTAINS
    SUBROUTINE init_offload()
        USE MPI_f08

        ! Local variables
        LOGICAL, ALLOCATABLE :: offload_ranks(:)
        INTEGER(intk), ALLOCATABLE :: failed_ranks(:)
        INTEGER(intk) :: i, nfail
        LOGICAL :: available

        available = probe_target()

        IF (myid == 0) THEN
            ALLOCATE(offload_ranks(numprocs))
        ELSE
            ALLOCATE(offload_ranks(0))
        END IF

        CALL MPI_Gather(available, 1, MPI_LOGICAL, offload_ranks, 1, &
            MPI_LOGICAL, 0, MPI_COMM_WORLD)

        IF (myid == 0) THEN
            nfail = 0
            DO i = 1, numprocs
                IF (.NOT. offload_ranks(i)) nfail = nfail + 1
            END DO

            WRITE(*, '("OPENMP OFFLOADING:")')
            IF (nfail == 0) THEN
                WRITE(*, '("    Target probe successful on all ranks")')
            ELSE
                ALLOCATE(failed_ranks(nfail), source=0_intk)
                nfail = 0
                DO i = 1, numprocs
                    IF (.NOT. offload_ranks(i)) THEN
                        nfail = nfail + 1
                        failed_ranks(nfail) = i - 1
                    END IF
                END DO
                WRITE(*, '("    Target probe FAILED on rank(s):", *(1X, I0))') &
                    failed_ranks
                DEALLOCATE(failed_ranks)
            END IF
            WRITE(*, '()')
        END IF

        DEALLOCATE(offload_ranks)
    END SUBROUTINE init_offload

    FUNCTION probe_target() RESULT(available)
#ifdef _MGLET_OFFLOAD_
        USE omp_lib
#endif
        ! Function arguments
        LOGICAL :: available

        ! Local variables
        LOGICAL :: is_initial_device

        is_initial_device = .TRUE.

#ifdef _MGLET_OFFLOAD_
        !$omp target map(from: is_initial_device)
        is_initial_device = omp_is_initial_device()
        !$omp end target
#endif

        available = .NOT. is_initial_device
    END FUNCTION probe_target
END MODULE probeoffload_mod
