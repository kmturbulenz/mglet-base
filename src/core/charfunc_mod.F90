MODULE charfunc_mod
    IMPLICIT NONE(type, external)
    PRIVATE

    CHARACTER(26), PARAMETER :: cap = 'ABCDEFGHIJKLMNOPQRSTUVWXYZ'
    CHARACTER(26), PARAMETER :: low = 'abcdefghijklmnopqrstuvwxyz'

    PUBLIC :: upper, lower
CONTAINS
    PURE FUNCTION upper(input) RESULT(output)
        ! Subroutine arguments
        CHARACTER(len=*), INTENT(IN) :: input
        CHARACTER(LEN(input)) :: output

        ! Local variables
        INTEGER :: i, ic

        ! First copy entire input to output
        output = input

        ! Then iterate over the result
        DO i = 1, LEN_TRIM(output)
            ic = INDEX(low, output(i:i))
            IF (ic > 0) output(i:i) = cap(ic:ic)
        END DO
    END FUNCTION upper


    PURE FUNCTION lower(input) RESULT(output)
        ! Subroutine arguments
        CHARACTER(len=*), INTENT(IN) :: input
        CHARACTER(LEN(input)) :: output

        ! Local variables
        INTEGER :: i, ic

        ! First copy entire input to output
        output = input

        ! Then iterate over the result
        DO i = 1, LEN_TRIM(output)
            ic = INDEX(cap, output(i:i))
            IF (ic > 0) output(i:i) = low(ic:ic)
        END DO
    END FUNCTION lower
END MODULE charfunc_mod
