message(STATUS "Checking if compiler supports basic Fortran-2018 features")
check_fortran_source_runs("
PROGRAM main

    USE, INTRINSIC :: iso_fortran_env, only : int8, int16

    IMPLICIT NONE ( type, external )

    INTEGER(kind=int8) :: m = -1
    INTEGER(kind=int16) :: n = -1

    n = n + m

END PROGRAM main" F2018_TEST_OK SRC_EXT "F90")
if(NOT F2018_TEST_OK)
    message(FATAL_ERROR "Compiler fails at 'implicit none(type,external)'
    (Fortran-2018 feature)")
endif()
