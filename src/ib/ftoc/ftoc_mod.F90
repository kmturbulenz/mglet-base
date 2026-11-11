MODULE ftoc_mod
    USE MPI_f08
    USE core_mod
    USE ftoc_core_mod
    USE ftoc_cpu_mod
    USE ftoc_gpu_mod

    IMPLICIT NONE (type, external)
    PRIVATE

    INTERFACE ftoc
        MODULE PROCEDURE :: ftoc_one, ftoc_multiple
    END INTERFACE ftoc

    ! No state in this module
    ! Acts only as a dispatcher for ftoc_cpu and ftoc_gpu
    ! Offload builds default to GPU execution; device=.FALSE. selects CPU.
    ! Non-offload builds always use the CPU implementation.

    PUBLIC :: ftoc, init_ftoc, finish_ftoc
CONTAINS
    SUBROUTINE ftoc_one(ilevel, ff, fc, flag, device)
        ! Subroutine arguments
        INTEGER(intk), INTENT(in) :: ilevel
        TYPE(field_t), INTENT(in) :: ff
        TYPE(field_t), INTENT(inout) :: fc
        CHARACTER(len=1), INTENT(in) :: flag
        LOGICAL, OPTIONAL, INTENT(in) :: device

#ifdef _MGLET_OFFLOAD_
        ! Local variables
        LOGICAL :: device2
#endif

        IF (.NOT. is_ftoc_core_init) CALL errr(__FILE__, __LINE__)

        ! The coarsest level cannot restrict values to a coarser level
        IF (ilevel == minlevel) RETURN

        CALL start_timer(220)

#ifdef _MGLET_OFFLOAD_
        device2 = .TRUE.
        IF (PRESENT(device)) THEN
            device2 = device
        END IF

        IF (device2) THEN
            CALL ftoc_gpu(ilevel, ff, fc, flag)
        ELSE
            CALL ftoc_cpu(ilevel, ff, fc, flag)
        END IF
#else
        CALL ftoc_cpu(ilevel, ff, fc, flag)
#endif

        CALL stop_timer(220)
    END SUBROUTINE ftoc_one


    SUBROUTINE ftoc_multiple(ilevel, v1, v2, v3, s1, s2, s3, device)
        ! Subroutine arguments
        INTEGER(intk), INTENT(in) :: ilevel
        TYPE(field_t), INTENT(inout), OPTIONAL :: v1, v2, v3, s1, s2, s3
        LOGICAL, OPTIONAL, INTENT(in) :: device

#ifdef _MGLET_OFFLOAD_
        ! Local variables
        LOGICAL :: device2
#endif

        IF (.NOT. is_ftoc_core_init) CALL errr(__FILE__, __LINE__)

        ! The coarsest level cannot restrict values to a coarser level
        IF (ilevel == minlevel) RETURN

        CALL start_timer(220)

#ifdef _MGLET_OFFLOAD_
        device2 = .TRUE.
        IF (PRESENT(device)) THEN
            device2 = device
        END IF

        IF (device2) THEN
            CALL ftoc_gpu(ilevel, v1, v2, v3, s1, s2, s3)
        ELSE
            CALL ftoc_cpu(ilevel, v1, v2, v3, s1, s2, s3)
        END IF
#else
        CALL ftoc_cpu(ilevel, v1, v2, v3, s1, s2, s3)
#endif

        CALL stop_timer(220)
    END SUBROUTINE ftoc_multiple


    SUBROUTINE init_ftoc()
        ! Subroutine arguments
        ! none...

        ! Local variables
        ! none...

        CALL init_ftoc_core()
        CALL init_ftoc_cpu()
#ifdef _MGLET_OFFLOAD_
        CALL init_ftoc_gpu()
#endif
    END SUBROUTINE init_ftoc


    SUBROUTINE finish_ftoc()
        ! Subroutine arguments
        ! none...

        ! Local variables
        ! none...

        CALL finish_ftoc_core()
        CALL finish_ftoc_cpu()
#ifdef _MGLET_OFFLOAD_
        CALL finish_ftoc_gpu()
#endif
    END SUBROUTINE finish_ftoc
END MODULE ftoc_mod
