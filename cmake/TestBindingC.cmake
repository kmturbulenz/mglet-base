message(STATUS "Checking if Fortran compiler supports C binding")
check_fortran_source_runs("
MODULE c_binding_mod

    USE, INTRINSIC :: ISO_C_BINDING, ONLY: C_PTR, C_LOC, C_LONG
    IMPLICIT NONE(type, external)
    PRIVATE

    INTERFACE
        INTEGER(c_long) FUNCTION testfun(crc, buf, len) BIND(C)
            USE, INTRINSIC :: ISO_C_BINDING, ONLY: C_LONG, C_PTR
            INTEGER(c_long), INTENT(in), VALUE :: crc
            TYPE(C_PTR), INTENT(in), VALUE :: buf
            INTEGER(c_long), INTENT(in), VALUE :: len
        END FUNCTION testfun
    END INTERFACE

END MODULE c_binding_mod


PROGRAM main

    USE c_binding_mod

END PROGRAM main" CBIND_TEST_OK SRC_EXT "F90")
if(NOT CBIND_TEST_OK)
    message(FATAL_ERROR "Compiler fails to use C bindings")
endif()
