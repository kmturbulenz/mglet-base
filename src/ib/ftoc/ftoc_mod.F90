MODULE ftoc_mod

    USE precision_mod
    USE timer_mod, ONLY: start_timer, stop_timer
    USE field_mod
    USE ftoc_core_mod
    USE ftoc1_mod

    INTERFACE ftoc
        MODULE PROCEDURE :: ftoc_one, ftoc_multiple
    END INTERFACE ftoc

    PUBLIC :: ftoc, init_ftoc, finish_ftoc

CONTAINS

    SUBROUTINE ftoc_one(ilevel, ff, fc, flag)
        ! Subroutine arguments
        INTEGER(intk), INTENT(in) :: ilevel
        TYPE(field_t), INTENT(in) :: ff
        TYPE(field_t), INTENT(inout) :: fc
        CHARACTER(len=1), INTENT(in) :: flag

#ifdef _MGLET_OFFLOAD_
        CALL ftoc2(ilevel, ff, fc, flag)
#else
        CALL ftoc1(ilevel, ff, fc, flag)
#endif

    END SUBROUTINE ftoc_one


    SUBROUTINE ftoc_multiple(ilevel, v1, v2, v3, s1, s2, s3)
        ! Subroutine arguments
        INTEGER(intk), INTENT(in) :: ilevel
        TYPE(field_t), INTENT(inout), OPTIONAL :: v1, v2, v3, s1, s2, s3

#ifdef _MGLET_OFFLOAD_
        CALL ftoc2(ilevel, v1, v2, v3, s1, s2, s3)
#else
        CALL ftoc1(ilevel, v1, v2, v3, s1, s2, s3)
#endif

    END SUBROUTINE ftoc_multiple


    SUBROUTINE init_ftoc()

        ! Initialize the core module first (provides infrastructure)
        CALL init_ftoc_core()

#ifdef _MGLET_OFFLOAD_
        CALL init_ftoc2()
#else
        CALL init_ftoc1()
#endif
    END SUBROUTINE init_ftoc


    SUBROUTINE finish_ftoc()

#ifdef _MGLET_OFFLOAD_
        CALL finish_ftoc2()
#else
        CALL finish_ftoc1()
#endif

        ! Finish the core module last
        CALL finish_ftoc_core()

    END SUBROUTINE finish_ftoc

END MODULE ftoc_mod
