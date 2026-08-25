MODULE par_ftoc_mod
    USE core_mod
    USE par_ftoc1_mod
    USE par_ftoc2_mod
    USE par_ftoc_core_mod

    IMPLICIT NONE (type, external)
    PRIVATE

    PUBLIC :: par_ftoc_norm, init_par_ftoc, finish_par_ftoc
CONTAINS
    SUBROUTINE par_ftoc_norm(ilevel, v1, v2, v3, sum, device)
        ! Subroutine arguments
        INTEGER(intk), INTENT(in) :: ilevel
        TYPE(field_t), INTENT(inout) :: v1, v2, v3
        LOGICAL, OPTIONAL, INTENT(in) :: sum, device

        ! Local variables
#ifdef _MGLET_OFFLOAD_
        LOGICAL :: device2
#endif

        ! The coarsest level cannot restrict values to a coarser level
        IF (ilevel == minlevel) RETURN

        IF (.NOT. is_par_ftoc_core_init) CALL errr(__FILE__, __LINE__)

        CALL start_timer(211)

#ifdef _MGLET_OFFLOAD_
        IF (PRESENT(device)) THEN
            device2 = device
        ELSE
            device2 = .FALSE.
        END IF

        IF (device2) THEN
            CALL par_ftoc2(ilevel, v1, v2, v3, sum)
        ELSE
            CALL par_ftoc1(ilevel, v1, v2, v3, sum)
        END IF
#else
        CALL par_ftoc1(ilevel, v1, v2, v3, sum)
#endif

        CALL stop_timer(211)
    END SUBROUTINE par_ftoc_norm


    SUBROUTINE init_par_ftoc()
        ! Subroutine arguments
        ! none...

        ! Local variables
        ! none...

        CALL init_par_ftoc_core()
        CALL init_par_ftoc1()
#ifdef _MGLET_OFFLOAD_
        CALL init_par_ftoc2()
#endif
    END SUBROUTINE init_par_ftoc


    SUBROUTINE finish_par_ftoc()
        ! Subroutine arguments
        ! none...

        ! Local variables
        ! none...

        CALL finish_par_ftoc1()
#ifdef _MGLET_OFFLOAD_
        CALL finish_par_ftoc2()
#endif
        CALL finish_par_ftoc_core()
    END SUBROUTINE finish_par_ftoc
END MODULE par_ftoc_mod
