MODULE ibcore_mod
    USE core_mod, ONLY: fort7, config_t, errr, intk, realk
    USE ibmodel_mod, ONLY: ibmodel_t

    IMPLICIT NONE(type, external)
    PRIVATE

    ABSTRACT INTERFACE
        FUNCTION constructor_i() RESULT(ib)
            IMPORT :: ibmodel_t
            CLASS(ibmodel_t), ALLOCATABLE :: ib
        END FUNCTION constructor_i
    END INTERFACE

    ! For storing IB type plugins
    TYPE :: ib_ext_t
        CHARACTER(len=16) :: name = REPEAT(" ", 16)
        PROCEDURE(constructor_i), POINTER, NOPASS :: constructor => NULL()
    END TYPE ib_ext_t
    TYPE(ib_ext_t) :: ib_ext(8)
    INTEGER(intk) :: n_ib_ext = 0

    CLASS(ibmodel_t), ALLOCATABLE :: ib

    REAL(realk), PROTECTED :: openaccur

    PUBLIC :: ib, init_ibcore, finish_ibcore, register_ib, openaccur

CONTAINS
    SUBROUTINE init_ibcore()
        ! Local variables
        TYPE(config_t) :: ibconf
        CHARACTER(len=16) :: ctyp
        INTEGER(intk) :: i
        LOGICAL :: found_ib

        ! Required values
        CALL fort7%get(ibconf, "/ib")
        CALL ibconf%get_value("/type", ctyp)
        CALL ibconf%get_value("/openaccur", openaccur, 0.125_realk)

        found_ib = .FALSE.
        DO i = 1, n_ib_ext
            IF (ctyp == TRIM(ib_ext(i)%name)) THEN
                ib = ib_ext(i)%constructor()
                found_ib = .TRUE.
            END IF
        END DO
        IF (.NOT. found_ib) THEN
            WRITE(*, '("Could not find IB type: ", A)') ctyp
            CALL errr(__FILE__, __LINE__)
        END IF
    END SUBROUTINE init_ibcore


    SUBROUTINE finish_ibcore()
        INTEGER(intk) :: i

        DO i = 1, n_ib_ext
            ib_ext(i)%name = REPEAT(" ", 16)
        END DO
        n_ib_ext = 0
        DEALLOCATE(ib)
    END SUBROUTINE finish_ibcore


    SUBROUTINE register_ib(name, constructor)
        CHARACTER(len=*), INTENT(in) :: name
        PROCEDURE(constructor_i) :: constructor

        IF (n_ib_ext > SIZE(ib_ext)) ERROR STOP
        IF (LEN_TRIM(name) > 16) ERROR STOP
        n_ib_ext = n_ib_ext + 1

        ib_ext(n_ib_ext)%name = TRIM(name)
        ib_ext(n_ib_ext)%constructor => constructor
    END SUBROUTINE register_ib
END MODULE ibcore_mod
