
#include "omprt_wrapper.h"

#include <stdio.h>
#include <math.h>
#include <omp.h>


void get_num_devices(int* n_dev) {
    *n_dev = omp_get_num_devices();
    return;
}


void get_default_device(int* def_dev) {
    *def_dev = omp_get_default_device();
    return;
}


void get_launched_on_gpu(int* l) {
    #pragma omp target map(l)
    {
        if (omp_is_initial_device() == 1) {
            printf("omp_is_initial_device() = 1");
            *l = 0;     // running on host
        }
        else {
            printf("omp_is_initial_device() = 0");
            *l = 1;     // running on device
        }
    }
    return;
}

#pragma omp declare target
void get_exp(float* b, float* a){
    *b = expf(*a);
    return;
}
#pragma omp end declare target