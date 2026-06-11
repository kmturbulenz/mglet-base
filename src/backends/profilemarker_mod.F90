MODULE profilemarker_mod
    USE, INTRINSIC :: ISO_C_BINDING, ONLY: c_char, c_null_char, c_int

    IMPLICIT NONE (type, external)
    PRIVATE

    INTERFACE
        SUBROUTINE range_push_c(name) &
            BIND(C, name="gpu_backend_range_push")
            IMPORT :: c_char
            CHARACTER(kind=c_char), INTENT(in) :: name(*)
        END SUBROUTINE range_push_c

        SUBROUTINE range_pop_c() &
            BIND(C, name="gpu_backend_range_pop")
             ! Subroutine arguments
             ! none...
        END SUBROUTINE range_pop_c

        SUBROUTINE mark_c(name) &
            BIND(C, name="gpu_backend_mark")
            IMPORT :: c_char
            CHARACTER(kind=c_char), INTENT(in) :: name(*)
        END SUBROUTINE mark_c
    END INTERFACE

    PUBLIC :: profile_marker_push, profile_marker_pop, profile_marker_set
CONTAINS
    SUBROUTINE profile_marker_push(name)
        ! Subroutine arguments
        CHARACTER(*), INTENT(in) :: name

        ! Local variables
        CHARACTER(len=LEN_TRIM(name)+1) :: trimmed_name

        trimmed_name = TRIM(name) // c_null_char
        CALL range_push_c(trimmed_name)
    END SUBROUTINE profile_marker_push

    SUBROUTINE profile_marker_pop()
        ! Subroutine arguments
        ! none...

        ! Local variables
        ! none...

        CALL range_pop_c()
    END SUBROUTINE profile_marker_pop

    SUBROUTINE profile_marker_set(name)
        ! Subroutine arguments
        CHARACTER(*), INTENT(in) :: name

        ! Local variables
        CHARACTER(len=LEN_TRIM(name)+1) :: trimmed_name

        trimmed_name = TRIM(name) // c_null_char
        CALL mark_c(trimmed_name)
    END SUBROUTINE profile_marker_set
END MODULE profilemarker_mod
