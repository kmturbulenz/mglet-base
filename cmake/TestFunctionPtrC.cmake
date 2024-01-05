message(STATUS "Checking if Fortran compiler supports C function pointers")
check_fortran_source_runs("
MODULE ptr_mod

    ! defining an abstract function interface
    abstract interface
        subroutine add_type(a, b, c)
            implicit none
            integer, intent(inout) :: a, b, c
        end subroutine add_type
    end interface

contains

    ! function implementation A
    subroutine add_a(a, b, c)
        implicit none
        integer, intent(inout) :: a, b, c
        c = a + b
    end subroutine add_a

    ! function implementation B
    subroutine add_b(a, b, c)
        implicit none
        integer, intent(inout) :: a, b, c
        c = a + b + 1
    end subroutine add_b

END MODULE ptr_mod


PROGRAM main

    USE, INTRINSIC :: iso_c_binding, only : &
        c_funloc, c_funptr, c_f_procpointer, c_null_funptr
    USE ptr_mod

    PROCEDURE(add_type), POINTER :: add_f_fptr
    TYPE(c_funptr) :: add_c_fptr
    INTEGER :: a, b, res

    ! clean preparation
    NULLIFY(add_f_fptr)
    add_c_fptr = c_null_funptr
    a = 1
    b = 2

    ! conversion to a C function pointer
    add_f_fptr => add_a
    add_c_fptr = c_funloc(add_f_fptr)

    ! testing Fortran function pointer
    add_f_fptr => add_b
    CALL add_f_fptr(a,b,res)
    IF ( res /= 4 ) stop 9

    ! recovering from C function pointer
    CALL c_f_procpointer(add_c_fptr, add_f_fptr)
    CALL add_f_fptr(a,b,res)
    IF ( res /= 3 ) stop 9

END PROGRAM main" FUNCTION_CPTR_TEST_OK SRC_EXT "F90")
if(NOT FUNCTION_CPTR_TEST_OK)
    message(FATAL_ERROR "Compiler fails at C function pointers")
endif()
