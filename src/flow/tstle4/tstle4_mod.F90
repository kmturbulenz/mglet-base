MODULE tstle4_mod

    USE core_mod
    USE tstle1_mod
    USE tstle2_mod

    IMPLICIT NONE(type, external)
    PRIVATE

    PUBLIC :: tstle4, tstle4_init, tstle4_finish

CONTAINS

    SUBROUTINE tstle4_init()
#ifdef __llvm__
        CALL tstle2_init()
#endif
    END SUBROUTINE tstle4_init


    SUBROUTINE tstle4(uo_f, vo_f, wo_f, u_f, v_f, w_f, ut_f, vt_f, wt_f, &
            p_f, g_f)
        ! Subroutine arguments
        TYPE(field_t), INTENT(inout) :: uo_f, vo_f, wo_f
        TYPE(field_t), INTENT(in) :: u_f, v_f, w_f
        TYPE(field_t), INTENT(in) :: ut_f, vt_f, wt_f, p_f, g_f
#ifdef __llvm__
        ! Uses the C-based execution kernels
        CALL tstle2(uo_f, vo_f, wo_f, u_f, v_f, w_f, ut_f, vt_f, wt_f, &
            p_f, g_f)
#else
        ! Uses the Fortran-based execution kernels (initialization not needed)
        CALL tstle1(uo_f, vo_f, wo_f, u_f, v_f, w_f, ut_f, vt_f, wt_f, &
            p_f, g_f)
#endif
    END SUBROUTINE tstle4


    SUBROUTINE tstle4_finish()
#ifdef __llvm__
        CALL tstle2_finish()
#endif
    END SUBROUTINE tstle4_finish

END MODULE tstle4_mod
