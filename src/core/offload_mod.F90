MODULE offload_mod
    USE comms_mod, ONLY: myid
    USE precision_mod, ONLY: intk, realk, c_intk, c_realk
    USE err_mod, ONLY: errr

    IMPLICIT NONE(type, external)
    PRIVATE

    ! Experimental interfaces to OpenMP runtime functions in C

    ! INTERFACE
    !     FUNCTION omp_get_num_devices_c() BIND(C) RESULT(res)
    !         USE iso_c_binding
    !         IMPORT :: c_intk
    !         INTEGER(kind=c_intk) :: res
    !     END FUNCTION omp_get_num_devices_c
    ! END INTERFACE

    INTERFACE
        FUNCTION omp_get_num_devices_c() BIND(C)
            IMPORT :: c_intk
            INTEGER(kind=c_intk) :: omp_get_num_devices_c
        END FUNCTION

        FUNCTION omp_is_initial_device_c() BIND(C)
            IMPORT :: c_intk
            INTEGER(kind=c_intk) :: omp_is_initial_device_c
        END FUNCTION

        FUNCTION omp_get_default_device_c() BIND(C)
            IMPORT :: c_intk
            INTEGER(kind=c_intk) :: omp_get_default_device_c
        END FUNCTION

        FUNCTION omp_get_initial_device_c() BIND(C)
            IMPORT :: c_intk
            INTEGER(kind=c_intk) :: omp_get_initial_device_c
        END FUNCTION

        FUNCTION omp_get_num_teams_c() BIND(C)
            IMPORT :: c_intk
            INTEGER(kind=c_intk) :: omp_get_num_teams_c
        END FUNCTION

        FUNCTION omp_get_team_num_c() BIND(C)
            IMPORT :: c_intk
            INTEGER(kind=c_intk) :: omp_get_team_num_c
        END FUNCTION

        FUNCTION omp_target_is_present_c(ptr, device) BIND(C)
            USE, INTRINSIC :: iso_c_binding, only : C_PTR
            USE precision_mod, ONLY: c_intk
            TYPE(C_PTR), VALUE :: ptr
            INTEGER(kind=c_intk), VALUE :: device
            INTEGER(kind=c_intk) :: omp_target_is_present_c
        END FUNCTION
    END INTERFACE

    ! Experimental interfaces to C math functions

    INTERFACE
        FUNCTION exp_c(a) BIND(C) RESULT(res)
            USE iso_c_binding
            IMPORT :: c_realk
            REAL(kind=c_realk), VALUE :: a
            REAL(kind=c_realk) :: res
        END FUNCTION exp_c

        FUNCTION sin_c(a) bind(C) RESULT(res)
            USE iso_c_binding
            IMPORT :: c_realk
            REAL(kind=c_realk), VALUE :: a
            REAL(kind=c_realk) :: res
        END FUNCTION sin_c

        FUNCTION cos_c(a) bind(C) RESULT(res)
            USE iso_c_binding
            IMPORT :: c_realk
            REAL(kind=c_realk), VALUE :: a
            REAL(kind=c_realk) :: res
        END FUNCTION cos_c

        FUNCTION asin_c(a) bind(C) RESULT(res)
            USE iso_c_binding
            IMPORT :: c_realk
            REAL(kind=c_realk), VALUE :: a
            REAL(kind=c_realk) :: res
        END FUNCTION asin_c

        FUNCTION acos_c(a) bind(C) RESULT(res)
            USE iso_c_binding
            IMPORT :: c_realk
            REAL(kind=c_realk), VALUE :: a
            REAL(kind=c_realk) :: res
        END FUNCTION acos_c

        FUNCTION sqrt_c(a) bind(C) RESULT(res)
            USE iso_c_binding
            IMPORT :: c_realk
            REAL(kind=c_realk), VALUE :: a
            REAL(kind=c_realk) :: res
        END FUNCTION sqrt_c

        FUNCTION cbrt_c(a) bind(C) RESULT(res)
            USE iso_c_binding
            IMPORT :: c_realk
            REAL(kind=c_realk), VALUE :: a
            REAL(kind=c_realk) :: res
        END FUNCTION cbrt_c
    END INTERFACE

    ! Public subroutines
    PUBLIC :: init_offload

CONTAINS

    SUBROUTINE init_offload()

#ifdef _MGLET_OFFLOAD_

        ! Local variables
        INTEGER(intk) :: num_devices
        INTEGER(intk) :: is_initial_device_cpu
        INTEGER(intk) :: is_initial_device_gpu
        LOGICAL :: has_gpu_support
        CHARACTER(len=3) :: message
        REAL(realk) :: a = 1.0
        REAL(realk) :: res = 0.0

        ! Querying number of devices (via C function)
        num_devices = omp_get_num_devices_c()

        ! Check if a kernel can launch on the CPU and GPU (via C function)
        is_initial_device_cpu = omp_is_initial_device_c()
        !$omp target map(is_initial_device_gpu)
        is_initial_device_gpu = omp_is_initial_device_c()
        !$omp end target

        IF (is_initial_device_cpu == 1 .AND. &
                is_initial_device_gpu == 0) THEN
            message = "YES"
            has_gpu_support = .TRUE.
        ELSE IF (is_initial_device_cpu == 1 .AND. &
                is_initial_device_gpu == 1) THEN
            message = "NO "
            has_gpu_support = .FALSE.
        ELSE
            CALL errr(__FILE__, __LINE__)
        END IF

        IF (myid == 0) THEN
            WRITE(*, '("OPENMP ENABLED")')
            WRITE(*, '("    Offload available:   ", A)') message
            IF (has_gpu_support) THEN
                WRITE(*, '("    Visible devices (rank 0): ", I0)') num_devices
            END IF
            WRITE(*, '()')
        END IF

        ! Experimental: test offloaded function

        !$omp target map(tofrom: res) map(to: a)
        res = exp_c(a)
        !$omp end target
        WRITE(*, '("Result of offloaded exp_c: ", F6.3)') res
#endif

    END SUBROUTINE init_offload

END MODULE offload_mod
