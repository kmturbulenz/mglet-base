MODULE timeintegrate_scalar_mod
    USE core_mod

    IMPLICIT NONE(type, external)
    PRIVATE

    PUBLIC :: timeintegrate_scalar

CONTAINS
    SUBROUTINE timeintegrate_scalar(itstep, ittot, timeph, dt, irk, rkscheme)
        ! Subroutine arguments
        INTEGER(intk), INTENT(in) :: itstep
        INTEGER(intk), INTENT(in) :: ittot
        REAL(realk), INTENT(in) :: timeph
        REAL(realk), INTENT(in) :: dt
        INTEGER(intk), INTENT(in) :: irk
        TYPE(rk_2n_t), INTENT(in) :: rkscheme

        ! Local variables
        ! None...

        CONTINUE
    END SUBROUTINE timeintegrate_scalar
END MODULE timeintegrate_scalar_mod
