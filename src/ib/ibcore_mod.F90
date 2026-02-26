MODULE ibcore_mod
    USE core_mod, ONLY: fort7, config_t, errr, intk, realk
    USE ibmodel_mod, ONLY: ibmodel_t

    IMPLICIT NONE(type, external)
    PRIVATE

    CLASS(ibmodel_t), ALLOCATABLE :: ib

    REAL(realk), PROTECTED :: openaccur

    PUBLIC :: ib, init_ibcore, finish_ibcore, openaccur


CONTAINS
    SUBROUTINE init_ibcore()
        ! Local variables
        TYPE(config_t) :: ibconf

        ! Required values
        CALL fort7%get(ibconf, "/ib")
        CALL ibconf%get_value("/openaccur", openaccur, 0.125_realk)
    END SUBROUTINE init_ibcore


    SUBROUTINE finish_ibcore()
        DEALLOCATE(ib)
    END SUBROUTINE finish_ibcore
END MODULE ibcore_mod
