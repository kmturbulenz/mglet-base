#ifndef __OMPRT_WRAPPER_H__
#define __OMPRT_WRAPPER_H__

#include "mglet_precision.h"

void get_num_devices(int*);
void get_default_device(int*);
void get_launched_on_gpu(int*);

void get_exp(mgletreal*, const mgletreal*);

#endif /* __OMPRT_WRAPPER_H__ */