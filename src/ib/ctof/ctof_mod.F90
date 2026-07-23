MODULE ctof_mod

    USE core_mod
    USE ctof1_mod
    USE ctof2_mod

    PUBLIC :: ctof, init_ctof, finish_ctof

CONTAINS

    SUBROUTINE ctof(ilevel, ff_f, fc_f)

        ! This subroutine acts as an interface for ctof (coarse-to-fine) calls.
        ! Dependent on the configuration, it will call either ctof1_mod
        ! or ctof2_mod, which are relevant for CPU and GPU, respectively.

        ! Subroutine arguments
        INTEGER(intk), INTENT(in) :: ilevel
        TYPE(field_t), INTENT(inout) :: ff_f
        TYPE(field_t), INTENT(in) :: fc_f

        CALL start_timer(230)
#ifdef _MGLET_OFFLOAD_
        CALL ctof2(ilevel, ff_f, fc_f)
#else
        CALL ctof1(ilevel, ff_f, fc_f)
#endif
        CALL stop_timer(230)
    END SUBROUTINE ctof


    SUBROUTINE init_ctof()
#ifdef _MGLET_OFFLOAD_
        CALL init_ctof2()
#else
        CALL init_ctof1()
#endif
    END SUBROUTINE init_ctof


    SUBROUTINE finish_ctof()
#ifdef _MGLET_OFFLOAD_
        CALL finish_ctof2()
#else
        CALL finish_ctof1()
#endif
    END SUBROUTINE finish_ctof

END MODULE ctof_mod

