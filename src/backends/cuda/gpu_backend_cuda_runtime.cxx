#include "gpu_backend_runtime_api.h"
#include "cuda/gpu_backend_cuda_tools.h"

#include <cuda_runtime.h>

void gpu_backend_runtime_init() {
    CUDA_CHECK(cudaFree(nullptr));
}

void gpu_backend_runtime_finalize() {
    CUDA_CHECK(cudaDeviceSynchronize());
}
