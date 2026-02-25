MODULE offload_mod
    USE comms_mod, ONLY: myid
    USE precision_mod, ONLY: intk, realk, c_intk
    USE err_mod, ONLY: errr

    IMPLICIT NONE(type, external)
    PRIVATE

    INTERFACE
        ! void get_num_devices(int* n_dev)
        SUBROUTINE get_num_devices(n_dev) BIND(C)
            IMPORT :: c_intk
            INTEGER(kind=c_intk), INTENT(INOUT) :: n_dev
        END SUBROUTINE get_num_devices
    END INTERFACE

    INTERFACE
        ! void get_default_device(int* def_dev)
        SUBROUTINE get_default_device(def_dev) BIND(C)
            IMPORT :: c_intk
            INTEGER(kind=c_intk), INTENT(INOUT) :: def_dev
        END SUBROUTINE get_default_device
    END INTERFACE

    INTERFACE
        ! void get_launched_on_gpuint* l)
        SUBROUTINE get_launched_on_gpu(l) BIND(C)
            IMPORT :: c_intk
            INTEGER(kind=c_intk), INTENT(INOUT) :: l
            ! 1 for GPU, 0 for CPU
        END SUBROUTINE get_launched_on_gpu
    END INTERFACE

    ! Public subroutines
    PUBLIC :: init_offload

CONTAINS

    SUBROUTINE init_offload()

#ifdef _MGLET_OFFLOAD_

        ! Local variables
        INTEGER(intk) :: num_devices
        INTEGER(intk) :: gpu_launched
        CHARACTER(len=3) :: message
        REAL(realk) :: a = 1.0

        ! Querying number of devices (via C function)
        CALL get_num_devices(num_devices)

        ! Check if a kernel can launch on the GPU (via C function)
        gpu_launched = -1
        CALL get_launched_on_gpu(gpu_launched)

        IF (gpu_launched == 1) THEN
            message = "YES"
        ELSE IF (gpu_launched == 0) THEN
            message = "NO "
        ELSE
            CALL errr(__FILE__, __LINE__)
        END IF

        IF (myid == 0) THEN
            WRITE(*, '("OPENMP ENABLED")')
            WRITE(*, '("    Offload available:   ", A)') message
            IF (gpu_launched == 1) THEN
                WRITE(*, '("    Visible devices (rank 0): ", I0)') num_devices
            END IF
            WRITE(*, '()')
        END IF
#endif

    END SUBROUTINE init_offload

END MODULE offload_mod
