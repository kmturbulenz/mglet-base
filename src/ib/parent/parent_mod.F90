MODULE parent_mod
    USE MPI_f08
    USE core_mod
    USE parent_core_mod
    USE parent_cpu_mod
    USE parent_gpu_mod

    IMPLICIT NONE (type, external)
    PRIVATE

    ! No state in this module
    ! Acts only as a dispatcher for parent_cpu and parent_gpu

    PUBLIC :: parent, init_parent, finish_parent
CONTAINS
    ! Offload builds default to GPU execution; device=.FALSE. selects CPU.
    ! Non-offload builds always use the CPU implementation.
    SUBROUTINE parent(ilevel, v1, v2, v3, s1, s2, s3, normal, device)
        ! Subroutine arguments
        INTEGER(intk), INTENT(in) :: ilevel
        TYPE(field_t), OPTIONAL, TARGET, INTENT(inout) :: v1, v2, v3, s1, s2, s3
        LOGICAL, OPTIONAL, INTENT(in) :: normal
        LOGICAL, OPTIONAL, INTENT(in) :: device

#ifdef _MGLET_OFFLOAD_
        ! Local variables
        LOGICAL :: device2
#endif

        ! The coarsest level cannot receive values from a parent grid
        IF (ilevel == minlevel) RETURN

        IF (.NOT. is_parent_core_init) CALL errr(__FILE__, __LINE__)

        CALL start_timer(210)

#ifdef _MGLET_OFFLOAD_
        device2 = .TRUE.
        IF (PRESENT(device)) THEN
            device2 = device
        END IF

        IF (device2) THEN
            CALL parent_gpu(ilevel, v1, v2, v3, s1, s2, s3, normal)
        ELSE
            CALL parent_cpu(ilevel, v1, v2, v3, s1, s2, s3, normal)
        END IF
#else
        CALL parent_cpu(ilevel, v1, v2, v3, s1, s2, s3, normal)
#endif

        CALL stop_timer(210)
    END SUBROUTINE parent


    SUBROUTINE init_parent()
        ! Subroutine arguments
        ! none...

        ! Local variables
        ! none...

        CALL init_parent_core()
        CALL init_parent_cpu()
#ifdef _MGLET_OFFLOAD_
        CALL init_parent_gpu()
#endif
    END SUBROUTINE init_parent


    SUBROUTINE finish_parent()
        ! Subroutine arguments
        ! none...

        ! Local variables
        ! none...

        CALL finish_parent_core()
        CALL finish_parent_cpu()
#ifdef _MGLET_OFFLOAD_
        CALL finish_parent_gpu()
#endif
    END SUBROUTINE finish_parent
END MODULE parent_mod
