#ifndef __EXPRTK_WRAPPER_H__
#define __EXPRTK_WRAPPER_H__

#include "mglet_precision.h"

extern "C" {
#include <ISO_Fortran_binding.h>
}

template <typename T>
void eval_expr(CFI_cdesc_t* res, const char* name, const char* expr,
    T rho, T gmol, T tu_level, T timeph,
    CFI_cdesc_t* x, CFI_cdesc_t* y, CFI_cdesc_t* z,
    CFI_cdesc_t* dx, CFI_cdesc_t* dy, CFI_cdesc_t* dz,
    CFI_cdesc_t* ddx, CFI_cdesc_t* ddy, CFI_cdesc_t* ddz,
    CFI_cdesc_t* xstag, CFI_cdesc_t* ystag, CFI_cdesc_t* zstag,
    int* ierr);

extern "C" {
    void eval_real_expr(CFI_cdesc_t* res, const char* name, const char* expr,
        mgletreal rho, mgletreal gmol, mgletreal tu_level, mgletreal timeph,
        CFI_cdesc_t* x, CFI_cdesc_t* y, CFI_cdesc_t* z,
        CFI_cdesc_t* dx, CFI_cdesc_t* dy, CFI_cdesc_t* dz,
        CFI_cdesc_t* ddx, CFI_cdesc_t* ddy, CFI_cdesc_t* ddz, int* ierr);
}

#endif /* __EXPRTK_WRAPPER_H__ */
