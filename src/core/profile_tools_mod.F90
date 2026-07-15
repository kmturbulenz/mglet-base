MODULE profile_tools_mod
    USE, INTRINSIC :: ISO_C_BINDING, ONLY: c_char, c_null_char

    IMPLICIT NONE (type, external)
    PRIVATE

    INTERFACE
        SUBROUTINE profile_range_push_c(name) &
            BIND(C, name="profile_range_push_impl")
            IMPORT :: c_char
            CHARACTER(kind=c_char), INTENT(in) :: name(*)
        END SUBROUTINE profile_range_push_c

        SUBROUTINE profile_range_pop_c() BIND(C, name="profile_range_pop_impl")
             ! Subroutine arguments
             ! none...
        END SUBROUTINE profile_range_pop_c

        SUBROUTINE profile_mark_c(name) BIND(C, name="profile_mark_impl")
            IMPORT :: c_char
            CHARACTER(kind=c_char), INTENT(in) :: name(*)
        END SUBROUTINE profile_mark_c
    END INTERFACE

    PUBLIC :: profile_range_push, profile_range_pop, profile_mark

CONTAINS

    SUBROUTINE profile_range_push(name)
        ! Subroutine arguments
        CHARACTER(*), INTENT(in) :: name

        ! Local variables
        CHARACTER(c_char), DIMENSION(LEN(name)+1) :: c_name

        c_name(1:LEN(name)) = TRANSFER(name, c_name)
        c_name(LEN_TRIM(name)+1) = c_null_char
        CALL profile_range_push_c(c_name)
    END SUBROUTINE profile_range_push

    SUBROUTINE profile_range_pop()
        ! Subroutine arguments
        ! none...

        ! Local variables
        ! none...

        CALL profile_range_pop_c()
    END SUBROUTINE profile_range_pop

    SUBROUTINE profile_mark(name)
        ! Subroutine arguments
        CHARACTER(*), INTENT(in) :: name

        ! Local variables
        CHARACTER(c_char), DIMENSION(LEN(name)+1) :: c_name

        c_name(1:LEN(name)) = TRANSFER(name, c_name)
        c_name(LEN_TRIM(name)+1) = c_null_char
        CALL profile_mark_c(c_name)
    END SUBROUTINE profile_mark
END MODULE profile_tools_mod
