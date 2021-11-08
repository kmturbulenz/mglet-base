MODULE scalar_mod
    ! Do not 'USE' core_mod or ib_mod here, because they will get exported
    ! (this module had no 'PRIVATE' default)

    USE scacore_mod
    USE timeintegrate_scalar_mod

    IMPLICIT NONE(type, external)

    ! PRIVATE ::

CONTAINS
    SUBROUTINE init_scalar()
        CALL init_scacore()
    END SUBROUTINE init_scalar


    SUBROUTINE finish_scalar
        CALL finish_scacore()
    END SUBROUTINE finish_scalar
END MODULE scalar_mod
