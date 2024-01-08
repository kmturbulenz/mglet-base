message(STATUS "Checking if Fortran compiler supports C data pointers")
check_fortran_source_runs("
PROGRAM main

    USE, INTRINSIC :: iso_c_binding, only : &
        c_loc, c_ptr, c_f_pointer, c_null_ptr

    IMPLICIT NONE

    TYPE(c_ptr) :: c_array_ptr
    REAL, POINTER, CONTIGUOUS :: f_array_ptr(:,:)
    REAL, TARGET :: arr(3,3) = 0.0
    INTEGER :: i, j

    ! clean preparation
    NULLIFY(f_array_ptr)
    c_array_ptr = c_null_ptr
    arr = 3.0

    ! conversion to a C pointer
    f_array_ptr => arr
    c_array_ptr = c_loc(f_array_ptr)

    ! re-conversion
    CALL c_f_pointer(c_array_ptr, f_array_ptr, [3,3] )

    DO i = 1, 3
        DO j = 1, 3
            IF ( f_array_ptr(j,i) /= 3.0 ) STOP 9
        END DO
    END DO

END PROGRAM main" DATA_CPTR_TEST_OK SRC_EXT "F90")
if(NOT DATA_CPTR_TEST_OK)
    message(FATAL_ERROR "Compiler fails at C data pointers")
endif()
