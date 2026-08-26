#include "mglet_precision.h"


extern "C" void setzero_real_c(size_t n, mgletreal* arr) {
    if (n == 0) return;

    #pragma omp target teams loop
    for (size_t i = 0; i < n; ++i) {
        arr[i] = 0.0;
    }
}


extern "C" void setzero_ifk_c(size_t n, mgletifk* arr) {
    if (n == 0) return;

    #pragma omp target teams loop
    for (size_t i = 0; i < n; ++i) {
        arr[i] = 0;
    }
}


extern "C" void copyarr_real_c(size_t n, mgletreal* dest,
        const mgletreal* src) {
    if (n == 0) return;

    #pragma omp target teams loop
    for (size_t i = 0; i < n; ++i) {
        dest[i] = src[i];
    }
}


extern "C" void rkstep_c(const size_t n, mgletreal* p, mgletreal* dp,
        const mgletreal* rhsp, mgletreal frhs, mgletreal dtfu) {
    if (n == 0) return;

    #pragma omp target teams loop
    for (size_t i = 0; i < n; ++i) {
        dp[i] = frhs*dp[i] + rhsp[i];
        p[i] = p[i] + dtfu*dp[i];
    }
}


extern "C" void accumulate_pcorr_c(size_t n, mgletreal* dp,
        const mgletreal* hilf) {
    if (n == 0) return;

    #pragma omp target teams loop
    for (size_t i = 0; i < n; ++i) {
        dp[i] = dp[i] + hilf[i];
    }
}
