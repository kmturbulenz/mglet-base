MODULE offload_mod
    USE comms_mod, ONLY: myid
    USE precision_mod, ONLY: intk

    IMPLICIT NONE(type, external)
    PRIVATE

    ! Public subroutines
    PUBLIC :: init_offload, is_on_device
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

    LOGICAL FUNCTION is_on_device(cptr, device_num) RESULT(res)
        USE ISO_C_BINDING, ONLY: c_ptr
#ifdef _MGLET_OFFLOAD_
        USE omp_lib
#endif
        ! Subroutine arguments
        TYPE(c_ptr), INTENT(in) :: cptr
        INTEGER(intk), OPTIONAL :: device_num
#ifdef _MGLET_OFFLOAD_
        ! Local variables
        INTEGER(intk) :: has_device_num

        IF (PRESENT(device_num)) THEN
            has_device_num = device_num
        ELSE
            has_device_num = omp_get_default_device()
        END IF

        res = (omp_target_is_present(cptr, has_device_num) /= 0)
#else
        res = .FALSE.
#endif
    END FUNCTION is_on_device
END MODULE offload_mod
