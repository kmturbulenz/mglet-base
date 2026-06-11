#pragma once
#include <iostream>

#include <hip/hip_runtime.h>

#ifndef HIP_CHECK
// https://rocm.docs.amd.com/projects/HIP/en/docs-7.2.4/how-to/hip_runtime_api/error_handling.html
#define HIP_CHECK(expression)                  \
{                                              \
    const hipError_t status = expression;      \
    if(status != hipSuccess){                  \
        std::cerr << "HIP error "              \
                  << status << ": "            \
                  << hipGetErrorString(status) \
                  << " at " << __FILE__ << ":" \
                  << __LINE__ << std::endl;    \
        std::exit(1);                          \
    }                                          \
}
#endif // HIP_CHECK
