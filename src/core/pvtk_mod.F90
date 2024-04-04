MODULE pvtk_mod
    USE, INTRINSIC :: ISO_C_BINDING, ONLY: c_char, c_null_char

    IMPLICIT NONE(type, external)
    PRIVATE

    INTERFACE
        ! void create_pvtk(const char *path, const char *prefix)
        SUBROUTINE create_pvtk(path, prefix) BIND(C)
            IMPORT :: c_char
            CHARACTER(kind=c_char), INTENT(IN) :: path(*)
            CHARACTER(kind=c_char), INTENT(IN) :: prefix(*)
        END SUBROUTINE create_pvtk
    END INTERFACE

    PUBLIC :: pvtk_directory

CONTAINS
    SUBROUTINE pvtk_directory(path, prefix)
        ! Subroutine arguments
        CHARACTER(len=*), INTENT(in) :: path, prefix

        ! Local variables
        CHARACTER(c_char), DIMENSION(LEN(path)+1) :: c_path
        CHARACTER(c_char), DIMENSION(LEN(prefix)+1) :: c_prefix

        ! Add trailing C_NULL_CHAR to key
        c_path(1:LEN(path)) = TRANSFER(path, c_path)
        c_path(LEN_TRIM(path)+1) = C_NULL_CHAR

        c_prefix(1:LEN(prefix)) = TRANSFER(prefix, c_prefix)
        c_prefix(LEN_TRIM(prefix)+1) = C_NULL_CHAR

        CALL create_pvtk(c_path, c_prefix)
    END SUBROUTINE pvtk_directory
END MODULE pvtk_mod
