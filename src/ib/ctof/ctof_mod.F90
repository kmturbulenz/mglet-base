MODULE ctof_mod

    USE core_mod

    USE ctof1_mod
    USE ctof2_mod
    USE ctof_core_mod

    PRIVATE

        ! Indicator for operation on GPU
        LOGICAL :: use_device = .FALSE.

    PUBLIC :: ctof, init_ctof, finish_ctof

CONTAINS

    SUBROUTINE ctof(ilevel, ff_f, fc_f, device)

        ! This subroutine acts as an interface for ctof (coarse-to-fine) calls.
        ! Dependent on the configuration, it will call either ctof1_mod
        ! or ctof2_mod, which are relevant for CPU and GPU, respectively.

        ! Subroutine arguments
        INTEGER(intk), INTENT(in) :: ilevel
        TYPE(field_t), INTENT(inout) :: ff_f
        TYPE(field_t), INTENT(in) :: fc_f
        LOGICAL, OPTIONAL, INTENT(in) :: device

        ! Check prerequisites for ctof operation
        IF (.NOT. has_infrastructure) THEN
            WRITE(*, *) "ctof_core_mod not initialized"
            CALL errr(__FILE__, __LINE__)
        END IF

        ! In case of GPU support, overwrite with provided device argument
        IF (PRESENT(device)) THEN
#ifdef _MGLET_OFFLOAD_
            use_device = device
#else
            use_device = .FALSE.
#endif
        END IF

        ! Both versions must be initialized to be ready for operation

        IF (use_device) THEN
            ! Use the GPU version
            CALL ctof2(ilevel, ff_f, fc_f)
        ELSE
            ! Use the CPU version
            CALL ctof1(ilevel, ff_f, fc_f)
        END IF

    END SUBROUTINE ctof


    SUBROUTINE init_ctof()

        CALL set_timer(230, "CTOF")
        CALL set_timer(231, "CTOF_BEGIN")
        CALL set_timer(232, "CTOF_END")
        CALL set_timer(235, "CTOF_PROLONG_FINISH")

        IF (has_infrastructure) THEN
            WRITE(*, *) "ctof_core_mod already initialized"
            CALL errr(__FILE__, __LINE__)
        END IF

        ! Initialize the communication infrastructure
        CALL init_core_ctof()

        ! Initialize the CPU version
        CALL init_ctof1()

        ! Initialize the GPU version
        CALL init_ctof2()

    END SUBROUTINE init_ctof


    SUBROUTINE finish_ctof()

        ! Finish the CPU version
        CALL finish_ctof1()

        ! Finish the GPU version
        CALL finish_ctof2()

        ! Check if ctof_core_mod still provides infrastructure and finish
        IF (has_infrastructure) THEN
            CALL finish_core_ctof()
        END IF

    END SUBROUTINE finish_ctof

END MODULE ctof_mod

