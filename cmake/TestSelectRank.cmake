message(STATUS "Checking if your Fortran compiler support 'SELECT RANK'")
check_fortran_source_runs("
MODULE test_m
    IMPLICIT NONE(type, external)
CONTAINS
    SUBROUTINE test(arr)
        INTEGER :: arr(..)
        SELECT RANK (arr)
        RANK (1)
            WRITE(*, *) 'Rank is 1'
        RANK DEFAULT
            WRITE(*, *) 'Rank is something else'
        END SELECT
    END SUBROUTINE test
END MODULE test_m
PROGRAM main
    USE test_m
    INTEGER, PARAMETER :: testarr(5) = [1, 2, 3, 4, 5]
    CALL test(testarr)
END PROGRAM" ASSUMED_RANK_OK SRC_EXT "F90")
if(NOT ASSUMED_RANK_OK)
    message(WARNING "Your Fortran compiler does not support 'SELECT RANK'")
endif()
