#pragma once

#ifdef __cplusplus
extern "C" {
#endif

void gpu_backend_range_push(const char* name);
void gpu_backend_range_pop();
void gpu_backend_mark(const char* name);

#ifdef __cplusplus
}
#endif
