MODULE offload_mod
    USE comms_mod, ONLY: myid
    USE precision_mod, ONLY: intk

    IMPLICIT NONE(type, external)
    PRIVATE

    ! Public subroutines
    PUBLIC :: init_offload
CONTAINS
    SUBROUTINE init_offload()
#ifdef _MGLET_OFFLOAD_
        USE omp_lib
        ! Local variables
        INTEGER(intk) :: num_devices
        CHARACTER(len=3) :: message
        LOGICAL :: is_initial_device

        !$omp target map(from: is_initial_device)
        is_initial_device = omp_is_initial_device()
        !$omp end target

        IF (is_initial_device) THEN
            message = "NO "
        ELSE
            message = "YES"
        END IF

        IF (myid == 0) THEN
            WRITE(*, '("OPENMP ENABLED")')
            WRITE(*, '("    Offload available:   ", A)') message
            IF (.NOT. is_initial_device) THEN
                num_devices = omp_get_num_devices()
                WRITE(*, '("    Visible devices (rank 0): ", I0)') num_devices
            END IF
            WRITE(*, '()')
        END IF
#endif
    END SUBROUTINE init_offload
END MODULE offload_mod
