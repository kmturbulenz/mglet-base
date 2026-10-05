#include <cstdint>

#include "errr.h"
#include "gpu_check.h"
#include "gpu_include.h"
#include "mapped_arr_view.h"

namespace mglet::gpu {

namespace {

__global__ void sipiter1_hyperplane_level_kernel(
    mgletreal* __restrict__ res,
    const mgletreal* __restrict__ rhs,
    const mgletreal* __restrict__ lw,
    const mgletreal* __restrict__ ls,
    const mgletreal* __restrict__ lb,
    const mgletreal* __restrict__ lpr,
    const mgletifk* __restrict__ mip,
    const mgletifk* __restrict__ idxsip,
    mgletint nmygridsonlvl,
    const mgletint* __restrict__ mygridsonlvl,
    const gridinfo_t* __restrict__ gridinfo,
    const mgletint* __restrict__ ip3d)
{
    const auto block_idx = blockIdx.x;
    if (block_idx >= nmygridsonlvl)
    {
        return;
    }

    // mygridsonlvl is a per-level slice of a rectangular array sized to the
    // level with the most grids (see mygridslvl in grids_mod.F90), so
    // shorter levels are zero-padded at the end. Skip the padding.
    const auto grid_id = mygridsonlvl[block_idx];
    if (grid_id == 0)
    {
        return;
    }

    const auto igrid = grid_id - 1;

    const auto kk = gridinfo[igrid].kk;
    const auto jj = gridinfo[igrid].jj;
    const auto ii = gridinfo[igrid].ii;
    const auto ip3 = ip3d[igrid] - 1;

    const auto n3dmin = 3 + 3 + 3;
    const auto n3dmax = (ii - 2) + (jj - 2) + (kk - 2);

    for (mgletint m = n3dmin; m <= n3dmax; ++m)
    {
        const mgletint lm = static_cast<mgletint>(mip[ip3 + m - 1]);
        const mgletint lp = static_cast<mgletint>(mip[ip3 + m]) - lm;

        for (mgletint ipp = 1 + threadIdx.x; ipp <= lp; ipp += blockDim.x)
        {
            const mgletint iacc = lm + ipp;

            const mgletint sip_idx = ip3 + iacc - 1;

            const mgletint local_idx = static_cast<mgletint>(idxsip[sip_idx]);
            const mgletint idx = ip3 + local_idx - 1;
            const mgletint idx_km = idx - 1;
            const mgletint idx_jm = idx - kk;
            const mgletint idx_im = idx - kk * jj;

            mgletreal val = (rhs[idx] + res[idx]) * lpr[sip_idx];
            val -= lb[sip_idx] * res[idx_km];
            val -= ls[sip_idx] * res[idx_jm];
            val -= lw[sip_idx] * res[idx_im];

            res[idx] = val;
        }

        __syncthreads();
    }
}

__global__ void sipiter2_hyperplane_level_kernel(
    mgletreal* __restrict__ dp,
    mgletreal* __restrict__ res,
    const mgletreal* __restrict__ ue,
    const mgletreal* __restrict__ un,
    const mgletreal* __restrict__ ut,
    const mgletifk* __restrict__ mip,
    const mgletifk* __restrict__ idxsip,
    const mgletint* __restrict__ mygridsonlvl,
    mgletint nmygridsonlvl,
    const gridinfo_t* __restrict__ gridinfo,
    const mgletint* __restrict__ ip3d)
{
    const auto block_idx = blockIdx.x;
    if (block_idx >= nmygridsonlvl) return;

    // mygridsonlvl is zero-padded past this level's true grid count, see
    // mygridslvl in grids_mod.F90. Skip the padding.
    const auto grid_id = mygridsonlvl[block_idx];
    if (grid_id == 0) return;

    const auto igrid = grid_id - 1;

    const auto kk = gridinfo[igrid].kk;
    const auto jj = gridinfo[igrid].jj;
    const auto ii = gridinfo[igrid].ii;
    const auto ip3 = ip3d[igrid] - 1;

    const mgletint n3dmin = 3 + 3 + 3;
    const mgletint n3dmax = (ii - 2) + (jj - 2) + (kk - 2);

    for (mgletint m = n3dmax; m >= n3dmin; --m) {
        const mgletint lm = static_cast<mgletint>(mip[ip3 + m - 1]);
        const mgletint lp = static_cast<mgletint>(mip[ip3 + m]) - lm;

        for (mgletint ipp = 1 + threadIdx.x; ipp <= lp; ipp += blockDim.x) {
            const mgletint iacc = lm + ipp;
            const mgletint sip_idx = ip3 + iacc - 1;

            const mgletint local_idx = static_cast<mgletint>(idxsip[sip_idx]);
            const mgletint idx    = ip3 + local_idx - 1;
            const mgletint idx_kp = idx + 1;
            const mgletint idx_jp = idx + kk;
            const mgletint idx_ip = idx + kk * jj;

            mgletreal val = res[idx];
            val -= ut[sip_idx] * res[idx_kp];
            val -= un[sip_idx] * res[idx_jp];
            val -= ue[sip_idx] * res[idx_ip];

            res[idx] = val;
        }

        __syncthreads();
    }

    const auto ni = ii - 4;
    const auto nj = jj - 4;
    const auto nk = kk - 4;
    if (ni > 0 && nj > 0 && nk > 0) {
        const auto n = ni * nj * nk;
        for (mgletint lin = threadIdx.x; lin < n; lin += blockDim.x) {
            const int k = 2 + lin % nk;
            const int j = 2 + (lin / nk) % nj;
            const int i = 2 + lin / ((std::int64_t)nk * nj);

            const auto idx = ip3 + k + j * kk + i * kk * jj;
            dp[idx] = dp[idx] + res[idx];
        }
    }
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

    const auto threads = ::dim3{64};
    const auto blocks = ::dim3{static_cast<unsigned>(nmygridsonlvl)};

    sipiter1_hyperplane_level_kernel<<<blocks, threads>>>(
        res.device_ptr(),
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
        ip3d.device_ptr());

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

    const auto threads = ::dim3{64};
    const auto blocks = ::dim3{static_cast<unsigned>(nmygridsonlvl)};

    sipiter2_hyperplane_level_kernel<<<blocks, threads>>>(
        dp.device_ptr(),
        res.device_ptr(),
        sipue.device_ptr(),
        sipun.device_ptr(),
        siput.device_ptr(),
        miphp.device_ptr(),
        idxhp.device_ptr(),
        mygridsonlvl.device_ptr(),
        nmygridsonlvl,
        gridinfo.device_ptr(),
        ip3d.device_ptr());

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
