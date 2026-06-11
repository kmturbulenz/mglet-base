#include "gpu_backend_marker_api.h"

#include <rocprofiler-sdk-roctx/roctx.h>

void gpu_backend_range_push(const char* name) {
    roctxRangePush(name);
}

void gpu_backend_range_pop() {
    roctxRangePop();
}

void gpu_backend_mark(const char* name) {
    roctxMark(name);
}
