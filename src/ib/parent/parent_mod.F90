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
        TYPE(field_t), OPTIONAL, INTENT(inout) :: v1, v2, v3, s1, s2, s3
        LOGICAL, OPTIONAL, INTENT(in) :: normal
        LOGICAL, OPTIONAL, INTENT(in) :: device

        ! Local variables
        LOGICAL :: device2

        IF (.NOT. is_parent_core_init) CALL errr(__FILE__, __LINE__)

        IF (PRESENT(device)) THEN
            device2 = device
        ELSE
            device2 = .FALSE.
        END IF

        CALL start_timer(210)

        IF (device2) THEN
            CALL parent2(ilevel, v1, v2, v3, s1, s2, s3, normal)
        ELSE
            CALL parent1(ilevel, v1, v2, v3, s1, s2, s3, normal)
        END IF

        CALL stop_timer(210)
    END SUBROUTINE parent


    SUBROUTINE init_parent()
        ! Subroutine arguments
        ! none...

        ! Local variables
        ! none...

        CALL init_parent_core()
        CALL init_parent1()
        CALL init_parent2()
    END SUBROUTINE init_parent


    SUBROUTINE finish_parent()
        ! Subroutine arguments
        ! none...

        ! Local variables
        ! none...

        CALL finish_parent_core()
        CALL finish_parent1()
        CALL finish_parent2()
    END SUBROUTINE finish_parent
END MODULE parent_mod
