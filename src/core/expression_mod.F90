MODULE expression_mod
    USE, INTRINSIC :: ISO_C_BINDING, ONLY: c_char, c_null_char
    USE precision_mod, ONLY: realk, c_realk

    IMPLICIT NONE(type, external)
    PRIVATE

    INTERFACE
        SUBROUTINE eval_real_expr(res, name, expr, rho, gmol, tu_level, &
                timeph, x, y, z, dx, dy, dz, ddx, ddy, ddz, &
                xstag, ystag, zstag) BIND(C)
            IMPORT :: c_realk, c_char
            REAL(c_realk), INTENT(INOUT) :: res(:, :, :)
            CHARACTER(c_char), INTENT(IN) :: name(*)
            CHARACTER(c_char), INTENT(IN) :: expr(*)
            REAL(c_realk), INTENT(IN), VALUE :: rho
            REAL(c_realk), INTENT(IN), VALUE :: gmol
            REAL(c_realk), INTENT(IN), VALUE :: tu_level
            REAL(c_realk), INTENT(IN), VALUE :: timeph
            REAL(c_realk), INTENT(IN) :: x(:)
            REAL(c_realk), INTENT(IN) :: y(:)
            REAL(c_realk), INTENT(IN) :: z(:)
            REAL(c_realk), INTENT(IN) :: dx(:)
            REAL(c_realk), INTENT(IN) :: dy(:)
            REAL(c_realk), INTENT(IN) :: dz(:)
            REAL(c_realk), INTENT(IN) :: ddx(:)
            REAL(c_realk), INTENT(IN) :: ddy(:)
            REAL(c_realk), INTENT(IN) :: ddz(:)
            REAL(c_realk), INTENT(IN) :: xstag(:)
            REAL(c_realk), INTENT(IN) :: ystag(:)
            REAL(c_realk), INTENT(IN) :: zstag(:)
        END SUBROUTINE
    END INTERFACE

    PUBLIC :: initial_condition

CONTAINS

    SUBROUTINE initial_condition(res, name, expr, rho, gmol, tu_level, &
            timeph, x, y, z, dx, dy, dz, ddx, ddy, ddz, xstag, ystag, zstag)
        ! Subroutine arguments
        REAL(realk), INTENT(INOUT) :: res(:, :, :)
        CHARACTER(LEN=*), INTENT(IN) :: name
        CHARACTER(LEN=*), INTENT(IN) :: expr
        REAL(realk), INTENT(IN) :: rho, gmol, tu_level, timeph
        REAL(realk), INTENT(IN) :: x(:)
        REAL(realk), INTENT(IN) :: y(:)
        REAL(realk), INTENT(IN) :: z(:)
        REAL(realk), INTENT(IN) :: dx(:)
        REAL(realk), INTENT(IN) :: dy(:)
        REAL(realk), INTENT(IN) :: dz(:)
        REAL(realk), INTENT(IN) :: ddx(:)
        REAL(realk), INTENT(IN) :: ddy(:)
        REAL(realk), INTENT(IN) :: ddz(:)
        REAL(realk), INTENT(IN) :: xstag(:)
        REAL(realk), INTENT(IN) :: ystag(:)
        REAL(realk), INTENT(IN) :: zstag(:)

        ! Local variables
        CHARACTER(C_CHAR), DIMENSION(LEN(expr)+1) :: c_expr
        CHARACTER(C_CHAR), DIMENSION(LEN(expr)+1) :: c_name

        ! Add trailing C_NULL_CHAR to expr and name
        c_expr = TRANSFER(expr, c_expr)
        c_expr(LEN_TRIM(expr)+1) = C_NULL_CHAR

        c_name = TRANSFER(name, c_name)
        c_name(LEN_TRIM(name)+1) = C_NULL_CHAR

        CALL eval_real_expr(res, c_name, c_expr, rho, gmol, tu_level, &
            timeph, x, y, z, dx, dy, dz, ddx, ddy, ddz, xstag, ystag, zstag)
    END SUBROUTINE initial_condition

END MODULE expression_mod
