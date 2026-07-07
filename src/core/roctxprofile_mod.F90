MODULE roctxprofile_mod
    USE err_mod, ONLY: errr

    IMPLICIT NONE (type, external)
    PRIVATE

    INTERFACE
        SUBROUTINE roctxrangepush_c(message) BIND(C, NAME="roctxRangePushA")
            USE ISO_C_BINDING, only: C_CHAR
            CHARACTER(C_CHAR) :: message(*)
        END SUBROUTINE roctxrangepush_c

        SUBROUTINE roctxrangepop_c() BIND(C, NAME="roctxRangePop")
            ! No arguments
        END SUBROUTINE roctxrangepop_c

        SUBROUTINE roctxmark_c(message) BIND(C, NAME="roctxMarkA")
            USE ISO_C_BINDING, only: C_CHAR
            CHARACTER(C_CHAR) :: message(*)
        END SUBROUTINE roctxmark_c
    END INTERFACE

    PUBLIC :: roctxrangepush, roctxrangepop, roctxmark
CONTAINS
    SUBROUTINE roctxrangepush(name)
        USE ISO_C_BINDING, only: C_NULL_CHAR
        ! Subroutine arguments
        CHARACTER(*), INTENT(in) :: name

        ! Local variables
        ! none...

#ifdef _MGLET_ROCTX_PROFILING_
        CALL roctxrangepush_c(name // C_NULL_CHAR)
#endif
    END SUBROUTINE roctxrangepush


    SUBROUTINE roctxrangepop()
        ! Subroutine arguments
        ! none...

        ! Local variables
        ! none...

#ifdef _MGLET_ROCTX_PROFILING_
        CALL roctxrangepop_c()
#endif
    END SUBROUTINE roctxrangepop


    SUBROUTINE roctxmark(name)
        USE ISO_C_BINDING, only: C_NULL_CHAR
        ! Subroutine arguments
        CHARACTER(*), INTENT(in) :: name

        ! Local variables
        ! none...

#ifdef _MGLET_ROCTX_PROFILING_
        CALL roctxmark_c(name // C_NULL_CHAR)
#endif
    END SUBROUTINE roctxmark
END MODULE roctxprofile_mod
