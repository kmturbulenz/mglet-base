#include <type_traits>

#include "errr.h"
#include "gpu_check.h"
#include "gpu_include.h"
#include "mapped_arr_view.h"

namespace mglet::gpu {

namespace {

// Only interior residuals need updating. Halos stay in global memory.
// CUDA's portable shared-memory limit is 48 KiB.
#if defined(_MGLET_HIP_)
constexpr mgletint sip_res_capacity = 64 * 1024 / sizeof(mgletreal);
#else
constexpr mgletint sip_res_capacity = 48 * 1024 / sizeof(mgletreal);
#endif

struct sip_launch_config {
    std::size_t shared_bytes = 0;
    mgletint cube_size = -1;
    mgletint mapping_grid = 0;
    std::size_t grid_count = 0;
};

sip_launch_config sip_config(
    std::size_t count, const mgletint* grids, const gridinfo_t* info)
{
    sip_launch_config config;
    std::size_t largest = 0;
    for (std::size_t g = 0; g < count; ++g) {
        if (grids[g] == 0) continue;
        ++config.grid_count;
        const auto& grid = info[grids[g] - 1];
        const auto cube = grid.ii == grid.jj && grid.jj == grid.kk ? grid.ii : 0;
        if (config.cube_size == -1) {
            config.cube_size = cube;
            config.mapping_grid = grids[g] - 1;
        }
        else if (config.cube_size != cube) config.cube_size = 0;
        if (grid.ii <= 4 || grid.jj <= 4 || grid.kk <= 4) continue;
        const auto size = std::size_t(grid.ii - 4) * (grid.jj - 4) * (grid.kk - 4);
        if (size <= sip_res_capacity && size > largest) largest = size;
    }
    config.shared_bytes = largest * sizeof(mgletreal);
    return config;
}

template <class Launch>
void sip_dispatch(mgletint cube_size, Launch launch) {
    switch (cube_size) {
    case 20: launch(std::integral_constant<mgletint, 20>{}); break;
    case 24: launch(std::integral_constant<mgletint, 24>{}); break;
    case 28: launch(std::integral_constant<mgletint, 28>{}); break;
    default: launch(std::integral_constant<mgletint, 0>{}); break;
    }
}

template <bool Backward>
unsigned sip_block_size(const sip_launch_config& config) {
#if defined(_MGLET_HIP_) && !defined(_MGLET_DOUBLE_PRECISION_)
    // Measured on gfx1201; larger 20^3 forward batches favor 320 threads.
    switch (config.cube_size) {
    case 20: return !Backward && config.grid_count >= 256 ? 320 : 256;
    case 24: return 384;
    case 28: return 1024;
    }
#endif
    return 512;
}

template <bool Backward, mgletint CubeSize = 0>
__global__ void sip_hyperplane_level_kernel(
    mgletreal* __restrict__ res,
    mgletreal* __restrict__ dp,
    const mgletreal* __restrict__ rhs,
    const mgletreal* __restrict__ ci,
    const mgletreal* __restrict__ cj,
    const mgletreal* __restrict__ ck,
    const mgletreal* __restrict__ lpr,
    const mgletifk* __restrict__ mip,
    const mgletifk* __restrict__ idxsip,
    mgletint nmygridsonlvl,
    const mgletint* __restrict__ mygridsonlvl,
    const gridinfo_t* __restrict__ gridinfo,
    const mgletint* __restrict__ ip3d,
    mgletint mapping_grid)
{
    const auto block_idx = blockIdx.x;
    if (block_idx >= nmygridsonlvl) return;

    // Level slices are zero-padded to the largest level's grid count.
    const auto grid_id = mygridsonlvl[block_idx];
    if (grid_id == 0) return;

    const auto igrid = grid_id - 1;

    const auto kk = CubeSize ? CubeSize : gridinfo[igrid].kk;
    const auto jj = CubeSize ? CubeSize : gridinfo[igrid].jj;
    const auto ii = CubeSize ? CubeSize : gridinfo[igrid].ii;
    const auto ip3 = ip3d[igrid] - 1;
    // Equal-sized grids have identical traversal tables; share their cache lines.
    const auto map_ip3 = CubeSize ? ip3d[mapping_grid] - 1 : ip3;
    extern __shared__ mgletreal cached_res[];
    const auto nk = kk - 4, nj = jj - 4, ni = ii - 4;
    if (nk <= 0 || nj <= 0 || ni <= 0) return;
    const auto interior_size = ni * nj * nk;
    // Specialize the address space so cached accesses use LDS instructions.
    auto sweep = [&](auto cache) {
        constexpr bool use_lds = decltype(cache)::value;
        const auto n3dmin = 3 + 3 + 3;
        const auto n3dmax = (ii - 2) + (jj - 2) + (kk - 2);
        constexpr mgletint step = Backward ? 1 : -1;
        for (mgletint m = Backward ? n3dmax : n3dmin;
             Backward ? m >= n3dmin : m <= n3dmax; m -= step) {
            const auto lm = static_cast<mgletint>(mip[map_ip3 + m - 1]);
            const auto lp = static_cast<mgletint>(mip[map_ip3 + m]) - lm;

            for (mgletint ipp = 1 + threadIdx.x; ipp <= lp; ipp += blockDim.x) {
                const auto iacc = lm + ipp;
                const auto sip_idx = ip3 + iacc - 1;
                const auto local =
                    static_cast<mgletint>(idxsip[map_ip3 + iacc - 1]) - 1;
                const auto idx = ip3 + local;
                mgletreal val = res[idx];
                if constexpr (!Backward) val = (rhs[idx] + val) * lpr[sip_idx];
                if constexpr (use_lds) {
                    const auto k = local % kk - 2;
                    const auto j = (local / kk) % jj - 2;
                    const auto i = local / (kk * jj) - 2;
                    const auto c = k + nk * (j + nj * i);
                    val -= ck[sip_idx] * ((Backward ? k < nk - 1 : k > 0)
                        ? cached_res[c + step] : res[idx + step]);
                    val -= cj[sip_idx] * ((Backward ? j < nj - 1 : j > 0)
                        ? cached_res[c + step * nk] : res[idx + step * kk]);
                    val -= ci[sip_idx] * ((Backward ? i < ni - 1 : i > 0)
                        ? cached_res[c + step * nk * nj] : res[idx + step * kk * jj]);
                    cached_res[c] = val;
                } else {
                    val -= ck[sip_idx] * res[idx + step];
                    val -= cj[sip_idx] * res[idx + step * kk];
                    val -= ci[sip_idx] * res[idx + step * kk * jj];
                    res[idx] = val;
                }
            }

            __syncthreads();
        }

        if constexpr (use_lds || Backward) {
            for (mgletint lin = threadIdx.x; lin < interior_size; lin += blockDim.x) {
                const auto k = 2 + lin % nk;
                const auto j = 2 + (lin / nk) % nj;
                const auto i = 2 + lin / (nk * nj);
                const auto idx = ip3 + k + j * kk + i * kk * jj;
                const auto val = use_lds ? cached_res[lin] : res[idx];
                if constexpr (use_lds) res[idx] = val;
                if constexpr (Backward) dp[idx] += val;
            }
        }
    };
    if (interior_size <= sip_res_capacity) sweep(std::true_type{});
    else sweep(std::false_type{});
}

void sipiter1_hyperplane_level(
    MappedArrView<mgletreal>(res),
    MappedArrView<const mgletreal> rhs,
    MappedArrView<const mgletreal> siplw,
    MappedArrView<const mgletreal> sipls,
    MappedArrView<const mgletreal> siplb,
    MappedArrView<const mgletreal> siplpr,
    MappedArrView<const mgletifk> miphp,
    MappedArrView<const mgletifk> idxhp,
    MappedArrView<const mgletint> mygridsonlvl,
    MappedArrView<const gridinfo_t> gridinfo,
    MappedArrView<const mgletint> ip3d)
{
    const auto nmygridsonlvl = mygridsonlvl.flat_size();

    const auto blocks = ::dim3{static_cast<unsigned>(nmygridsonlvl)};
    const auto config = sip_config(
        nmygridsonlvl, mygridsonlvl.host_ptr(), gridinfo.host_ptr());
    const auto threads = ::dim3{sip_block_size<false>(config)};

    sip_dispatch(config.cube_size, [&](auto size) {
        sip_hyperplane_level_kernel<false, decltype(size)::value>
            <<<blocks, threads, config.shared_bytes>>>(
            res.device_ptr(),
            nullptr,
            rhs.device_ptr(),
            siplw.device_ptr(),
            sipls.device_ptr(),
            siplb.device_ptr(),
            siplpr.device_ptr(),
            miphp.device_ptr(),
            idxhp.device_ptr(),
            nmygridsonlvl,
            mygridsonlvl.device_ptr(),
            gridinfo.device_ptr(),
            ip3d.device_ptr(),
            config.mapping_grid);
    });

    GPU_CHECK(gpuGetLastError());
    GPU_CHECK(gpuDeviceSynchronize());
}

void sipiter2_hyperplane_level(
    MappedArrView<mgletreal> dp,
    MappedArrView<mgletreal> res,
    MappedArrView<const mgletreal> sipue,
    MappedArrView<const mgletreal> sipun,
    MappedArrView<const mgletreal> siput,
    MappedArrView<const mgletifk> miphp,
    MappedArrView<const mgletifk> idxhp,
    MappedArrView<const mgletint> mygridsonlvl,
    MappedArrView<const gridinfo_t> gridinfo,
    MappedArrView<const mgletint> ip3d)
{
    const auto nmygridsonlvl = mygridsonlvl.flat_size();

    const auto blocks = ::dim3{static_cast<unsigned>(nmygridsonlvl)};
    const auto config = sip_config(
        nmygridsonlvl, mygridsonlvl.host_ptr(), gridinfo.host_ptr());
    const auto threads = ::dim3{sip_block_size<true>(config)};

    sip_dispatch(config.cube_size, [&](auto size) {
        sip_hyperplane_level_kernel<true, decltype(size)::value>
            <<<blocks, threads, config.shared_bytes>>>(
            res.device_ptr(),
            dp.device_ptr(),
            nullptr,
            sipue.device_ptr(),
            sipun.device_ptr(),
            siput.device_ptr(),
            nullptr,
            miphp.device_ptr(),
            idxhp.device_ptr(),
            nmygridsonlvl,
            mygridsonlvl.device_ptr(),
            gridinfo.device_ptr(),
            ip3d.device_ptr(),
            config.mapping_grid);
    });

    GPU_CHECK(gpuGetLastError());
    GPU_CHECK(gpuDeviceSynchronize());
}

} // namespace

extern "C" void sipiter1_hyperplane_level_c(
    CFI_cdesc_t* res,
    CFI_cdesc_t* rhs,
    CFI_cdesc_t* siplw,
    CFI_cdesc_t* sipls,
    CFI_cdesc_t* siplb,
    CFI_cdesc_t* siplpr,
    CFI_cdesc_t* miphp,
    CFI_cdesc_t* idxhp,
    CFI_cdesc_t* mygridsonlvl,
    CFI_cdesc_t* gridinfo,
    CFI_cdesc_t* ip3d)
{
    mglet::gpu::sipiter1_hyperplane_level(
        mglet::gpu::MappedArrView<mgletreal>(res),
        mglet::gpu::MappedArrView<const mgletreal>(rhs),
        mglet::gpu::MappedArrView<const mgletreal>(siplw),
        mglet::gpu::MappedArrView<const mgletreal>(sipls),
        mglet::gpu::MappedArrView<const mgletreal>(siplb),
        mglet::gpu::MappedArrView<const mgletreal>(siplpr),
        mglet::gpu::MappedArrView<const mgletifk>(miphp),
        mglet::gpu::MappedArrView<const mgletifk>(idxhp),
        mglet::gpu::MappedArrView<const mgletint>(mygridsonlvl),
        mglet::gpu::MappedArrView<const gridinfo_t>(gridinfo),
        mglet::gpu::MappedArrView<const mgletint>(ip3d));
}

extern "C" void sipiter2_hyperplane_level_c(
    CFI_cdesc_t* dp,
    CFI_cdesc_t* res,
    CFI_cdesc_t* sipue,
    CFI_cdesc_t* sipun,
    CFI_cdesc_t* siput,
    CFI_cdesc_t* miphp,
    CFI_cdesc_t* idxhp,
    CFI_cdesc_t* mygridsonlvl,
    CFI_cdesc_t* gridinfo,
    CFI_cdesc_t* ip3d)
{
    mglet::gpu::sipiter2_hyperplane_level(
        mglet::gpu::MappedArrView<mgletreal>(dp),
        mglet::gpu::MappedArrView<mgletreal>(res),
        mglet::gpu::MappedArrView<const mgletreal>(sipue),
        mglet::gpu::MappedArrView<const mgletreal>(sipun),
        mglet::gpu::MappedArrView<const mgletreal>(siput),
        mglet::gpu::MappedArrView<const mgletifk>(miphp),
        mglet::gpu::MappedArrView<const mgletifk>(idxhp),
        mglet::gpu::MappedArrView<const mgletint>(mygridsonlvl),
        mglet::gpu::MappedArrView<const gridinfo_t>(gridinfo),
        mglet::gpu::MappedArrView<const mgletint>(ip3d));
}

} // namespace mglet::gpu
