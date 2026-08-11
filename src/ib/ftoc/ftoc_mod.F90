MODULE ftoc_mod
    USE MPI_f08
    USE core_mod
    USE ftoc_core_mod
    USE ftoc1_mod
    USE ftoc2_mod

    IMPLICIT NONE (type, external)
    PRIVATE

    INTERFACE ftoc
        MODULE PROCEDURE :: ftoc_one, ftoc_multiple
    END INTERFACE ftoc

    ! No state in this module
    ! Acts only as a dispatcher for ftoc1 and ftoc2

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

        ! The coarsest level cannot restrict values to a coarser level
        IF (ilevel == minlevel) RETURN

        IF (.NOT. is_ftoc_core_init) CALL errr(__FILE__, __LINE__)

        CALL start_timer(220)

#ifdef _MGLET_OFFLOAD_
        IF (PRESENT(device)) THEN
            device2 = device
        ELSE
            device2 = .FALSE.
        END IF

        IF (device2) THEN
            CALL ftoc2(ilevel, ff, fc, flag)
        ELSE
            CALL ftoc1(ilevel, ff, fc, flag)
        END IF
#else
        CALL ftoc1(ilevel, ff, fc, flag)
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

        CALL start_timer(220)

        IF (.NOT. is_ftoc_core_init) CALL errr(__FILE__, __LINE__)

#ifdef _MGLET_OFFLOAD_
        IF (PRESENT(device)) THEN
            device2 = device
        ELSE
            device2 = .FALSE.
        END IF

        IF (device2) THEN
            CALL ftoc2(ilevel, v1, v2, v3, s1, s2, s3)
        ELSE
            CALL ftoc1(ilevel, v1, v2, v3, s1, s2, s3)
        END IF
#else
        CALL ftoc1(ilevel, v1, v2, v3, s1, s2, s3)
#endif

        CALL stop_timer(220)
    END SUBROUTINE ftoc_multiple


    SUBROUTINE init_ftoc()
        ! Subroutine arguments
        ! none...

        ! Local variables
        ! none...

        CALL init_ftoc_core()
        CALL init_ftoc1()
#ifdef _MGLET_OFFLOAD_
        CALL init_ftoc2()
#endif
    END SUBROUTINE init_ftoc


    SUBROUTINE finish_ftoc()
        ! Subroutine arguments
        ! none...

        ! Local variables
        ! none...

        CALL finish_ftoc_core()
        CALL finish_ftoc1()
#ifdef _MGLET_OFFLOAD_
        CALL finish_ftoc2()
#endif
    END SUBROUTINE finish_ftoc
END MODULE ftoc_mod
