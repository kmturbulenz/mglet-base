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

#ifndef _MGLET_OFFLOAD_

        ! Uses the Fortran-based implementation on CPU
        CALL tstle1(uo_f, vo_f, wo_f, u_f, v_f, w_f, ut_f, vt_f, wt_f, &
            p_f, g_f)
#else

#if defined(__flang__) &&defined(__amdflang__)
        ! Uses the Fortran-based device execution kernels with AMD compiler
        CALL tstle1(uo_f, vo_f, wo_f, u_f, v_f, w_f, ut_f, vt_f, wt_f, &
            p_f, g_f)
#else
        ! Uses the C-based execution kernels e.g. for upstream LLVM compiler
        CALL tstle2(uo_f, vo_f, wo_f, u_f, v_f, w_f, ut_f, vt_f, wt_f, &
            p_f, g_f)
#endif

#endif

    END SUBROUTINE tstle4


    SUBROUTINE tstle4_finish()
#ifdef __llvm__
        CALL tstle2_finish()
#endif
    END SUBROUTINE tstle4_finish

END MODULE tstle4_mod
