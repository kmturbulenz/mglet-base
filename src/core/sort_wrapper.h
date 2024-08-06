#ifndef __SORT_WRAPPER_H__
#define __SORT_WRAPPER_H__

extern "C" {
#include <ISO_Fortran_binding.h>
}

template <typename T>
void sort_fortran_array(CFI_cdesc_t* idx, CFI_cdesc_t* data);

extern "C" {
    void sort_int(CFI_cdesc_t* data,  CFI_cdesc_t* idx);
    void sort_long(CFI_cdesc_t* data,  CFI_cdesc_t* idx);
    void sort_float(CFI_cdesc_t* data,  CFI_cdesc_t* idx);
    void sort_double(CFI_cdesc_t* data,  CFI_cdesc_t* idx);
}

#endif /* __SORT_WRAPPER_H__ */
