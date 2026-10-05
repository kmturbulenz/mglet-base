#pragma once

#if defined(_MGLET_CUDA_)
#include <cuda_runtime.h>
#elif defined(_MGLET_HIP_)
#include <hip/hip_runtime.h>
#endif
