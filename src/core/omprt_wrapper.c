
#include "omprt_wrapper.h"

#include "mglet_precision.h"
#include <stdio.h>
#include <math.h>
#include <omp.h>

// OpenMP runtime functions

/* Return number of OpenMP devices */
mgletint omp_get_num_devices_c(void) {
    return (mgletint) omp_get_num_devices();
}

/* Returns true if current code is executing on the initial device */
#pragma omp declare target
mgletint omp_is_initial_device_c(void) {
    return (mgletint) omp_is_initial_device();
}
#pragma omp end declare target

/* Returns the default device number */
mgletint omp_get_default_device_c(void) {
    return (mgletint) omp_get_default_device();
}

/* Returns the initial device number */
mgletint omp_get_initial_device_c(void) {
    return (mgletint) omp_get_initial_device();
}

/* Returns the number of teams in the current parallel region */
mgletint omp_get_num_teams_c(void) {
    return (mgletint) omp_get_num_teams();
}

/* Returns the team number of the calling thread */
mgletint omp_get_team_num_c(void) {
    return (mgletint) omp_get_team_num();
}

/* Returns true if the pointer is present on the target device */
mgletint omp_target_is_present_c(void *ptr, mgletint device) {
    return (mgletint) omp_target_is_present(ptr, (int)device);
}


// Math functions

#define MGLET_EXP(x)   _Generic((x), float: expf,   double: exp)(x)
#define MGLET_SIN(x)   _Generic((x), float: sinf,   double: sin)(x)
#define MGLET_COS(x)   _Generic((x), float: cosf,   double: cos)(x)
#define MGLET_ASIN(x)  _Generic((x), float: asinf,  double: asin)(x)
#define MGLET_ACOS(x)  _Generic((x), float: acosf,  double: acos)(x)
#define MGLET_SQRT(x)  _Generic((x), float: sqrtf,  double: sqrt)(x)
#define MGLET_CBRT(x)  _Generic((x), float: cbrtf,  double: cbrt)(x)


#pragma omp declare target

mgletreal exp_c(const mgletreal a)   { return MGLET_EXP(a); }
mgletreal sin_c(const mgletreal a)   { return MGLET_SIN(a); }
mgletreal cos_c(const mgletreal a)   { return MGLET_COS(a); }
mgletreal asin_c(const mgletreal a)  { return MGLET_ASIN(a); }
mgletreal acos_c(const mgletreal a)  { return MGLET_ACOS(a); }
mgletreal sqrt_c(const mgletreal a)  { return MGLET_SQRT(a); }
mgletreal cbrt_c(const mgletreal a)  { return MGLET_CBRT(a); }

#pragma omp end declare target



// Sinus function

