MODULE sym_name_mod
    USE err_mod, ONLY: errr
    USE charfunc_mod, ONLY: lower

    IMPLICIT NONE(type, external)
    PRIVATE

    PUBLIC :: sym_name1, sym_name2
CONTAINS
    SUBROUTINE sym_name1(mod, fun, sym)
        ! Subroutine arguments
        CHARACTER(len=*), INTENT(in) :: mod
        CHARACTER(len=*), INTENT(in) :: fun
        CHARACTER(len=*), INTENT(out) :: sym

        ! Local variables
        CHARACTER(len=LEN(mod)) :: mod_l
        CHARACTER(len=LEN(fun)) :: fun_l

        mod_l = lower(mod)
        fun_l = lower(fun)

        ! Mangling:
        !    __
        !    mod_l
        !    _MOD_
        !    fun_l
        ! i.e. 7 chars + fun and mod name length
        IF (LEN(sym) < LEN(mod_l) + LEN(fun_l) + 7) THEN
            CALL errr(__FILE__, __LINE__)
        END IF

        sym = REPEAT(" ", LEN(sym))
        WRITE(sym, '("__", A, "_MOD_", A)') mod_l, fun_l
    END SUBROUTINE sym_name1


    FUNCTION sym_name2(mod, fun) RESULT(sym)
        ! Function arguments
        CHARACTER(len=*), INTENT(in) :: mod
        CHARACTER(len=*), INTENT(in) :: fun
        CHARACTER(len=LEN(mod) + LEN(fun) + 7) :: sym

        CALL sym_name1(mod, fun, sym)
    END FUNCTION sym_name2
END MODULE sym_name_mod
