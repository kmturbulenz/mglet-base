#ifndef __OMPRT_WRAPPER_H__
#define __OMPRT_WRAPPER_H__

#include "mglet_precision.h"


mgletint omp_get_num_devices_c(void);
mgletint omp_is_initial_device_c(void);
mgletint omp_get_default_device_c(void);
mgletint omp_get_initial_device_c(void);
mgletint omp_get_num_teams_c(void);
mgletint omp_get_team_num_c(void);
mgletint omp_target_is_present_c(void *ptr, mgletint device);

mgletreal exp_c(const mgletreal);
mgletreal sin_c(const mgletreal);
mgletreal cos_c(const mgletreal);
mgletreal asin_c(const mgletreal);
mgletreal acos_c(const mgletreal);
mgletreal sqrt_c(const mgletreal);
mgletreal cbrt_c(const mgletreal);


#endif /* __OMPRT_WRAPPER_H__ */