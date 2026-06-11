#include "gpu_backend_marker_api.h"

#include <nvtx3/nvToolsExt.h>

void gpu_backend_range_push(const char* name) {
    nvtxRangePushA(name);
}

void gpu_backend_range_pop() {
    nvtxRangePop();
}

void gpu_backend_mark(const char* name) {
    nvtxMarkA(name);
}
