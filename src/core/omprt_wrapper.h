#ifndef __OMPRT_WRAPPER_H__
#define __OMPRT_WRAPPER_H__

void get_num_devices(int*);
void get_default_device(int*);
void get_launched_on_gpu(int*);
void get_exp(float*, float*);

#endif /* __OMPRT_WRAPPER_H__ */