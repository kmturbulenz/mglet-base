MODULE parent_mod
    USE MPI_f08
    USE core_mod
    USE parent_core_mod
    USE parent1_mod
    USE parent2_mod

    IMPLICIT NONE (type, external)
    PRIVATE

    ! No state in this module
    ! Acts only as a dispatcher for parent1 and parent2

    PUBLIC :: parent, init_parent, finish_parent
CONTAINS
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
        ! If offloading is enabled, the optional device argument determines
        ! whether device or host version of parent is called.
        IF (PRESENT(device)) THEN
            device2 = device
        ELSE
            device2 = .FALSE.
        END IF

        IF (device2) THEN
            CALL parent2(ilevel, v1, v2, v3, s1, s2, s3, normal)
        ELSE
            CALL parent1(ilevel, v1, v2, v3, s1, s2, s3, normal)
        END IF
#else
        ! Always call parent1 if offloading is not enabled
        CALL parent1(ilevel, v1, v2, v3, s1, s2, s3, normal)
#endif

        CALL stop_timer(210)
    END SUBROUTINE parent


    SUBROUTINE init_parent()
        ! Subroutine arguments
        ! none...

        ! Local variables
        ! none...

        CALL init_parent_core()
        CALL init_parent1()
#ifdef _MGLET_OFFLOAD_
        CALL init_parent2()
#endif
    END SUBROUTINE init_parent


    SUBROUTINE finish_parent()
        ! Subroutine arguments
        ! none...

        ! Local variables
        ! none...

        CALL finish_parent_core()
        CALL finish_parent1()
#ifdef _MGLET_OFFLOAD_
        CALL finish_parent2()
#endif
    END SUBROUTINE finish_parent
END MODULE parent_mod
