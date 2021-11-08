#ifndef __MGLET_PRECISION_H__
#define __MGLET_PRECISION_H__

#include <ISO_Fortran_binding.h>

#ifdef _MGLET_DOUBLE_PRECISION_
typedef double mgletreal;
#define CFI_type_mgletreal CFI_type_double
#else
typedef float mgletreal;
#define CFI_type_mgletreal CFI_type_float
#endif

#ifdef _MGLET_INT64_
typedef long long mgletint;
#define CFI_type_mgletint CFI_type_long_long
#else
typedef int mgletint;
#define CFI_type_mgletint CFI_type_int
#endif

#endif /* __MGLET_PRECISION_H__ */
