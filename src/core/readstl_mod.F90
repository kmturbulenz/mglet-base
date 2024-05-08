! Module for wrapping around various STL reading functions
MODULE readstl_mod
    USE precision_mod, ONLY: realk, c_realk
    USE, INTRINSIC :: ISO_C_BINDING, ONLY: c_char, c_int, c_null_char
    USE err_mod, ONLY: errr

    IMPLICIT NONE(type, external)
    PRIVATE

    INTERFACE
        SUBROUTINE readstl(data, filename, ierr) BIND(C)
            IMPORT :: c_char, c_int, c_realk
            REAL(c_realk), ALLOCATABLE, INTENT(inout) :: data(:, :, :)
            CHARACTER(c_char), DIMENSION(*), INTENT(in) :: filename
            INTEGER(c_int), INTENT(out) :: ierr
        END SUBROUTINE readstl
    END INTERFACE

    ! Public data items
    PUBLIC :: stl_read

CONTAINS
    SUBROUTINE stl_read(vertices, filename)
        ! Subroutine arguments
        REAL(realk), ALLOCATABLE, INTENT(inout) :: vertices(:, :, :)
        CHARACTER(len=*), INTENT(in) :: filename

        ! Local variables
        INTEGER(c_int) :: ierr
        CHARACTER(c_char), DIMENSION(LEN(filename)+1) :: c_filename

        ! Convert the filename to a C string
        c_filename(1:LEN(filename)) = TRANSFER(filename, c_filename)
        c_filename(LEN_TRIM(filename)+1) = c_null_char

        CALL readstl(vertices, c_filename, ierr)
        IF (ierr /= 0) CALL errr(__FILE__, __LINE__)
    END SUBROUTINE stl_read
END MODULE readstl_mod
