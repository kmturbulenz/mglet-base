MODULE ctof_mod
    USE core_mod

    USE ctof_cpu_mod
    USE ctof_gpu_mod
    USE ctof_core_mod

    IMPLICIT NONE (type, external)
    PRIVATE

    PUBLIC :: ctof, init_ctof, finish_ctof
CONTAINS
    ! Offload builds default to GPU execution; device=.FALSE. selects CPU.
    ! Non-offload builds always use the CPU implementation.
    SUBROUTINE ctof(ilevel, ff_f, fc_f, device)
        ! Subroutine arguments
        INTEGER(intk), INTENT(in) :: ilevel
        TYPE(field_t), INTENT(inout) :: ff_f
        TYPE(field_t), INTENT(in) :: fc_f
        LOGICAL, OPTIONAL, INTENT(in) :: device

#ifdef _MGLET_OFFLOAD_
        ! Local variables
        LOGICAL :: device2
#endif

        ! The coarsest level cannot obtain values from a coarser level
        IF (ilevel == minlevel) RETURN

        IF (.NOT. is_ctof_core_init) CALL errr(__FILE__, __LINE__)

        CALL start_timer(230)

#ifdef _MGLET_OFFLOAD_
        device2 = .TRUE.
        IF (PRESENT(device)) THEN
            device2 = device
        END IF

        IF (device2) THEN
            CALL ctof_gpu(ilevel, ff_f, fc_f)
        ELSE
            CALL ctof_cpu(ilevel, ff_f, fc_f)
        END IF
#else
        CALL ctof_cpu(ilevel, ff_f, fc_f)
#endif

        CALL stop_timer(230)
    END SUBROUTINE ctof


    SUBROUTINE init_ctof()
        CALL init_ctof_core()
        CALL init_ctof_cpu()
#ifdef _MGLET_OFFLOAD_
        CALL init_ctof_gpu()
#endif
    END SUBROUTINE init_ctof


    SUBROUTINE finish_ctof()
        CALL finish_ctof_cpu()
#ifdef _MGLET_OFFLOAD_
        CALL finish_ctof_gpu()
#endif
        CALL finish_ctof_core()
    END SUBROUTINE finish_ctof
END MODULE ctof_mod
