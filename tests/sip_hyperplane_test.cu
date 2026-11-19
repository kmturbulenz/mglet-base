#include <algorithm>
#include <cmath>
#include <iostream>
#include <limits>
#include <string_view>
#include <vector>

// Run with --benchmark for GPU-event timings at 200-1000 grids/level.
#include "../src/backend/kernels.cu"

#if defined(_MGLET_HIP_)
#define GPU_API(name) hip##name
#else
#define GPU_API(name) cuda##name
#endif

constexpr unsigned block_sizes[] = {
    64, 128, 192, 256, 320, 384, 448, 512, 768, 1024
};

template <class T>
struct device_array {
    T* ptr = nullptr;
    std::size_t size;

    explicit device_array(const std::vector<T>& values) : size(values.size()) {
        GPU_CHECK(GPU_API(Malloc)(reinterpret_cast<void**>(&ptr), size * sizeof(T)));
        GPU_CHECK(GPU_API(Memcpy)(ptr, values.data(), size * sizeof(T),
                                 GPU_API(MemcpyHostToDevice)));
    }

    ~device_array() { GPU_CHECK(GPU_API(Free)(ptr)); }

    std::vector<T> download() const {
        std::vector<T> values(size);
        GPU_CHECK(GPU_API(Memcpy)(values.data(), ptr, size * sizeof(T),
                                 GPU_API(MemcpyDeviceToHost)));
        return values;
    }
};

int run_case(mgletint cube_size = 0, bool benchmark = false) {
    // Include nonzero offsets, rectangular/empty interiors, and a grid too
    // large for LDS. Padding and an unselected grid must remain untouched.
    std::vector<gridinfo_t> grids(benchmark ? 1000 : (cube_size ? 4 : 8));
    const mgletint dims[][3] = {
        {12, 12, 12}, {9, 15, 11}, {5, 5, 5},
        {4, 8, 7}, {8, 72, 76}, {10, 10, 10}
    };
    std::vector<mgletint> offsets;
    mgletint total = 17;
    for (std::size_t g = 0; g < grids.size(); ++g) {
        grids[g].ii = cube_size ? cube_size : (g < 6 ? dims[g][0] : 24);
        grids[g].jj = cube_size ? cube_size : (g < 6 ? dims[g][1] : 24);
        grids[g].kk = cube_size ? cube_size : (g < 6 ? dims[g][2] : 24);
        if (!cube_size && (g == 6 || g == 7)) {
            grids[g].ii = g == 6 ? 20 : 12;
            grids[g].jj = grids[g].kk = 36;
        }
        offsets.push_back(total + 1);
        total += grids[g].ii * grids[g].jj * grids[g].kk;
    }
    std::vector<mgletint> selected;
    if (!cube_size) selected = {1, 0, 5, 2, 4, 3, 0};
    for (mgletint g = cube_size ? 1 : 7;
         g <= static_cast<mgletint>(grids.size()); ++g)
        selected.push_back(g);
    if (cube_size && !benchmark) {
        std::reverse(selected.begin(), selected.end());
        selected.insert(selected.begin(), 0);
        selected.push_back(0);
    }
    std::vector<mgletifk> mip(total, 0), indices(total, -1);
    for (std::size_t g = 0; g < grids.size(); ++g) {
        const auto ii = grids[g].ii, jj = grids[g].jj, kk = grids[g].kk;
        const auto base = offsets[g] - 1;
        mgletint count = 0;
        for (mgletint m = 9; m <= ii + jj + kk - 6; ++m) {
            mip[base + m - 1] = count;
            std::vector<mgletifk> plane;
            for (mgletint i = 2; i < ii - 2; ++i) {
                for (mgletint j = 2; j < jj - 2; ++j) {
                    const auto k = m - 3 - i - j;
                    if (k >= 2 && k < kk - 2)
                        plane.push_back(1 + k + j * kk + i * kk * jj);
                }
            }
            std::sort(plane.begin(), plane.end());
            for (auto index : plane) indices[base + count++] = index;
            mip[base + m] = count;
        }
    }

    std::vector<mgletreal> res(total), dp(total), rhs(total);
    std::vector<mgletreal> lb(total), ls(total), lw(total), lpr(total);
    for (mgletint n = 0; n < total; ++n) {
        res[n] = mgletreal((n % 37) - 18) / 37;
        dp[n] = mgletreal(n % 23) / 23;
        rhs[n] = mgletreal((n % 29) - 14) / 29;
        lb[n] = mgletreal(1 + n % 3) / 100;
        ls[n] = mgletreal(1 + n % 5) / 100;
        lw[n] = mgletreal(1 + n % 7) / 100;
        lpr[n] = mgletreal(0.7);
    }
    device_array d_res(res), d_dp(dp), d_rhs(rhs);
    device_array d_lb(lb), d_ls(ls), d_lw(lw), d_lpr(lpr);
    device_array d_mip(mip), d_indices(indices);
    device_array d_selected(selected), d_offsets(offsets);
    device_array d_grids(grids);

    unsigned threads = 512;
    mgletint active_grids = selected.size();
    const auto config = mglet::gpu::sip_config(
        selected.size(), selected.data(), grids.data());
    auto forward = [&] {
        mglet::gpu::sip_dispatch(config.cube_size, [&](auto size) {
            mglet::gpu::sip_hyperplane_level_kernel<false, decltype(size)::value>
                <<<active_grids, threads, config.shared_bytes>>>(
                d_res.ptr, nullptr, d_rhs.ptr, d_lw.ptr, d_ls.ptr, d_lb.ptr, d_lpr.ptr,
                d_mip.ptr, d_indices.ptr, active_grids, d_selected.ptr,
                d_grids.ptr, d_offsets.ptr, config.mapping_grid);
        });
        GPU_CHECK(gpuGetLastError());
    };
    auto backward = [&] {
        mglet::gpu::sip_dispatch(config.cube_size, [&](auto size) {
            mglet::gpu::sip_hyperplane_level_kernel<true, decltype(size)::value>
                <<<active_grids, threads, config.shared_bytes>>>(
                d_res.ptr, d_dp.ptr, nullptr, d_lw.ptr, d_ls.ptr, d_lb.ptr, nullptr,
                d_mip.ptr, d_indices.ptr, active_grids, d_selected.ptr,
                d_grids.ptr, d_offsets.ptr, config.mapping_grid);
        });
        GPU_CHECK(gpuGetLastError());
    };
    auto reference = [&](bool reverse) {
        for (auto id : selected) {
            if (id == 0) continue;
            const auto& grid = grids[id - 1];
            const auto kk = grid.kk, jj = grid.jj, ii = grid.ii;
            const auto base = offsets[id - 1] - 1;
            const auto last = ii + jj + kk - 6;
            for (auto m = reverse ? last : 9;
                 reverse ? m >= 9 : m <= last; m += reverse ? -1 : 1) {
                for (auto p = mip[base + m - 1]; p < mip[base + m]; ++p) {
                    const auto sip = base + p;
                    const auto idx = base + indices[sip] - 1;
                    auto val = reverse ? res[idx] : (rhs[idx] + res[idx]) * lpr[sip];
                    const auto step = reverse ? 1 : -1;
                    val -= lb[sip] * res[idx + step];
                    val -= ls[sip] * res[idx + step * kk];
                    val -= lw[sip] * res[idx + step * kk * jj];
                    res[idx] = val;
                    if (reverse) dp[idx] += val;
                }
            }
        }
    };
    auto matches = [](const auto& expected, const auto& actual) {
        for (std::size_t n = 0; n < expected.size(); ++n) {
            const auto tolerance = mgletreal(32) *
                std::numeric_limits<mgletreal>::epsilon() *
                std::max(mgletreal(1), std::abs(expected[n]));
            if (!std::isfinite(actual[n]) ||
                std::abs(expected[n] - actual[n]) > tolerance) {
                std::cerr << "Mismatch at " << n << ": " << expected[n]
                          << " != " << actual[n] << '\n';
                return false;
            }
        }
        return true;
    };

    for (unsigned block_size : block_sizes) {
        threads = block_size;
        forward();
        reference(false);
        if (!matches(res, d_res.download())) return 1;
        backward();
        reference(true);
        if (!matches(res, d_res.download()) || !matches(dp, d_dp.download())) return 1;
        if (benchmark) break;
    }
    if (!benchmark) return 0;

    GPU_API(Event_t) start, stop;
    GPU_CHECK(GPU_API(EventCreate)(&start));
    GPU_CHECK(GPU_API(EventCreate)(&stop));
    auto measure = [&](bool reverse) {
        if (reverse) backward();
        else forward();
        std::vector<double> timings;
        constexpr int repeats = 20;
        for (int sample = 0; sample < 5; ++sample) {
            GPU_CHECK(GPU_API(EventRecord)(start));
            for (int repeat = 0; repeat < repeats; ++repeat) {
                if (reverse) backward();
                else forward();
            }
            GPU_CHECK(GPU_API(EventRecord)(stop));
            GPU_CHECK(GPU_API(EventSynchronize)(stop));
            float milliseconds;
            GPU_CHECK(GPU_API(EventElapsedTime)(&milliseconds, start, stop));
            timings.push_back(milliseconds * 1000.0 / repeats);
        }
        std::sort(timings.begin(), timings.end());
        return timings[timings.size() / 2];
    };
    for (mgletint count : {200, 256, 257, 300, 400, 500, 750, 1000}) {
        active_grids = count;
        for (unsigned block_size : block_sizes) {
            threads = block_size;
            int forward_blocks, backward_blocks;
            mglet::gpu::sip_dispatch(config.cube_size, [&](auto size) {
                GPU_CHECK(GPU_API(OccupancyMaxActiveBlocksPerMultiprocessor)(
                    &forward_blocks,
                    mglet::gpu::sip_hyperplane_level_kernel<false, decltype(size)::value>,
                    threads, config.shared_bytes));
                GPU_CHECK(GPU_API(OccupancyMaxActiveBlocksPerMultiprocessor)(
                    &backward_blocks,
                    mglet::gpu::sip_hyperplane_level_kernel<true, decltype(size)::value>,
                    threads, config.shared_bytes));
            });
            const auto forward_us = measure(false);
            const auto backward_us = measure(true);
            std::cout << cube_size << ',' << count << ',' << threads << ','
                      << forward_us << ',' << backward_us << ','
                      << forward_blocks << ',' << backward_blocks << std::endl;
        }
    }
    GPU_CHECK(GPU_API(EventDestroy)(start));
    GPU_CHECK(GPU_API(EventDestroy)(stop));
    return 0;
}

int main(int argc, char** argv) {
    if (argc == 1) {
        if (run_case() != 0) return 1;
        for (mgletint size : {20, 24, 28}) {
            if (run_case(size) != 0) return 1;
        }
        std::cout << "All SIP sweeps match CPU reference, including halos.\n";
        return 0;
    }
    if (argc != 2 || std::string_view(argv[1]) != "--benchmark") {
        std::cerr << "Usage: sip_hyperplane_test [--benchmark]\n";
        return 1;
    }
    std::cout << "cube_with_halos,grids,threads,forward_us,backward_us,"
                 "forward_blocks_per_cu,backward_blocks_per_cu\n";
    for (mgletint size : {20, 24, 28}) {
        if (run_case(size, true) != 0) return 1;
    }
}
