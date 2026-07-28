MODULE ctof_mod
    USE core_mod

    USE ctof1_mod
    USE ctof2_mod
    USE ctof_core_mod

    IMPLICIT NONE (type, external)
    PRIVATE

    PUBLIC :: ctof, init_ctof, finish_ctof
CONTAINS
    ! This subroutine acts as an interface for ctof (coarse-to-fine) calls.
    ! Dependent on the configuration, it will call either ctof1_mod
    ! or ctof2_mod, which are relevant for CPU and GPU, respectively.
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
        ! If offloading is enabled, the optional device argument determines
        ! whether device or host version of parent is called.
        IF (PRESENT(device)) THEN
            device2 = device
        ELSE
            device2 = .FALSE.
        END IF

        IF (device2) THEN
            CALL ctof2(ilevel, ff_f, fc_f)
        ELSE
            CALL ctof1(ilevel, ff_f, fc_f)
        END IF
#else
        ! Always call ctof1 if offloading is not enabled
        CALL ctof1(ilevel, ff_f, fc_f)
#endif

        CALL stop_timer(230)
    END SUBROUTINE ctof


    SUBROUTINE init_ctof()
        CALL init_ctof_core()
        CALL init_ctof1()
        CALL init_ctof2()
    END SUBROUTINE init_ctof


    SUBROUTINE finish_ctof()
        CALL finish_ctof1()
        CALL finish_ctof2()
        CALL finish_ctof_core()
    END SUBROUTINE finish_ctof
END MODULE ctof_mod

