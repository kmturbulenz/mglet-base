#include "gpu_backend_runtime_api.h"
#include "hip/gpu_backend_hip_tools.h"

#include <hip/hip_runtime.h>

void gpu_backend_runtime_init() {
    HIP_CHECK(hipFree(nullptr));
}

void gpu_backend_runtime_finalize() {
    HIP_CHECK(hipDeviceSynchronize());
}
