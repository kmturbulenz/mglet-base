MODULE mgletmath_mod
    USE precision_mod, ONLY: intk, realk
    USE field_mod, ONLY: field_t
    USE grids_mod, ONLY: mygrids, nmygrids

    IMPLICIT NONE(type, external)
    PRIVATE

    PUBLIC :: zero_ghostlayers
CONTAINS
    SUBROUTINE zero_ghostlayers(field)
        ! Subroutine arguments
        TYPE(field_t), INTENT(inout) :: field

        ! Local variables
        INTEGER(intk) :: i, igrid
        REAL(realk), POINTER, CONTIGUOUS :: f(:, :, :)

        DO i = 1, nmygrids
            igrid = mygrids(i)
            CALL field%get_ptr(f, igrid)
            CALL zero_ghostlayers_grid(f)
        END DO
    END SUBROUTINE zero_ghostlayers


    PURE SUBROUTINE zero_ghostlayers_grid(field)
        ! Subroutine arguments
        REAL(realk), INTENT(inout) :: field(:, :, :)

        ! Local variables
        INTEGER(intk) :: kk, jj, ii

        kk = SIZE(field, 1)
        jj = SIZE(field, 2)
        ii = SIZE(field, 3)

        field(:, :, 1:2) = 0.0
        field(:, :, ii-1:ii) = 0.0
        field(:, 1:2, :) = 0.0
        field(:, jj-1:jj, :) = 0.0
        field(1:2, :, :) = 0.0
        field(kk-1:kk, :, :) = 0.0
    END SUBROUTINE zero_ghostlayers_grid
END MODULE mgletmath_mod
