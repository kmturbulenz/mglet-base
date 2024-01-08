message(STATUS "Checking if Fortran compiler supports (pure) elemental")
check_fortran_source_runs("
MODULE elemental_mod

    IMPLICIT NONE(type, external)
    PUBLIC :: af, bf

CONTAINS

    ELEMENTAL SUBROUTINE af(a, b)
        REAL, INTENT(inout) :: a, b
        a = a + b
    END SUBROUTINE af

    PURE ELEMENTAL REAL FUNCTION bf(a, b) RESULT(res)
        REAL, INTENT(in) :: a, b
        res = a + b
    END FUNCTION bf

END MODULE elemental_mod

PROGRAM main

    USE elemental_mod, ONLY: af, bf
    IMPLICIT NONE

    REAL :: a(10) = 1.0
    REAL :: b(10) = 2.0
    REAL :: c(10) = -1.0
    c = bf( a, b )
    IF ( ABS( MAXVAL(c) - 3.0 ) > 0.00001 ) STOP 9

END PROGRAM main" ELEMENTAL_TEST_OK SRC_EXT "F90")
if(NOT ELEMENTAL_TEST_OK)
    message(FATAL_ERROR "Compiler fails at 'elemental' or / and 'pure elemental'")
endif()
