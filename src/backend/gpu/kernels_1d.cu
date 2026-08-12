#include <cstdint>

#include <ISO_Fortran_binding.h>

#include "errr.h"
#include "gpu_check.h"
#include "gpu_include.h"
#include "mapped_arr_view.h"

namespace mglet::gpu {

namespace {

__global__ void set_field_arr_realk_kernel(mgletreal* __restrict__ arr,
                                           mgletint n, mgletreal val) {
    const auto i = blockIdx.x * blockDim.x + threadIdx.x;

    if (i >= n) {
        return;
    }

    arr[i] = val;
}

__global__ void set_field_arr_ifk_kernel(mgletifk* __restrict__ arr, mgletint n,
                                         mgletifk val) {
    const auto i = blockIdx.x * blockDim.x + threadIdx.x;

    if (i >= n) {
        return;
    }

    arr[i] = val;
}

__global__ void add_field_arr_realk_kernel(mgletreal* __restrict__ lhs,
                                           const mgletreal* __restrict__ rhs,
                                           mgletint n) {
    const auto i = blockIdx.x * blockDim.x + threadIdx.x;

    if (i >= n) {
        return;
    }

    lhs[i] = lhs[i] + rhs[i];
}

} // namespace

void set_field_arr_realk(MappedArrView<mgletreal> farr, mgletreal val) {
    const auto n = farr.flat_size();

    if (n == 0) {
        return;
    }

    const unsigned block_size = 256;
    const unsigned grid_size = (n + block_size - 1) / block_size;

    set_field_arr_realk_kernel<<<grid_size, block_size>>>(farr.device_ptr(), n,
                                                          val);

    GPU_CHECK(gpuGetLastError());
    GPU_CHECK(gpuDeviceSynchronize());
}

void set_field_arr_ifk(MappedArrView<mgletifk> farr, mgletifk val) {
    const auto n = farr.flat_size();

    if (n == 0) {
        return;
    }

    const unsigned block_size = 256;
    const unsigned grid_size = (n + block_size - 1) / block_size;

    set_field_arr_ifk_kernel<<<grid_size, block_size>>>(farr.device_ptr(), n,
                                                        val);

    GPU_CHECK(gpuGetLastError());
    GPU_CHECK(gpuDeviceSynchronize());
}

void add_field_arr_realk(MappedArrView<mgletreal> lhs,
                         MappedArrView<const mgletreal> rhs) {
    const auto n_lhs = lhs.flat_size();
    const auto n_rhs = rhs.flat_size();

    if (n_lhs != n_rhs) {
        MGLET_ERRR();
    }

    if (n_lhs == 0) {
        return;
    }

    const unsigned threads = 256;
    const unsigned blocks = (n_lhs + threads - 1) / threads;

    add_field_arr_realk_kernel<<<blocks, threads>>>(lhs.device_ptr(),
                                                    rhs.device_ptr(), n_lhs);

    GPU_CHECK(gpuGetLastError());
    GPU_CHECK(gpuDeviceSynchronize());
}

} // namespace mglet::gpu

#ifdef _MGLET_USE_BACKEND_

extern "C" void set_field_arr_realk_c(CFI_cdesc_t* farr, mgletreal val) {
    mglet::gpu::set_field_arr_realk(mglet::gpu::MappedArrView<mgletreal>(farr),
                                    val);
}

extern "C" void set_field_arr_ifk_c(CFI_cdesc_t* farr, mgletifk val) {
    mglet::gpu::set_field_arr_ifk(mglet::gpu::MappedArrView<mgletifk>(farr),
                                  val);
}

extern "C" void add_field_arr_realk_c(CFI_cdesc_t* lhs, CFI_cdesc_t* rhs) {
    mglet::gpu::add_field_arr_realk(
        mglet::gpu::MappedArrView<mgletreal>(lhs),
        mglet::gpu::MappedArrView<const mgletreal>(rhs));
}

#endif // _MGLET_USE_BACKEND_
