message(STATUS "Checking if your C library supports 'qsort_r'")
check_c_source_runs("
#define _GNU_SOURCE
#include <assert.h>
#include <stdio.h>
#include <stdlib.h>

typedef struct
{
    size_t n;
    void* ptr;
} ctx_t;

int cmp_r(const void *a, const void *b, void *arg) {
    int ia = *(const int *)a;
    int ib = *(const int *)b;

    ctx_t* ctx = (ctx_t*)arg;
    int n = ctx->n;
    float *data = ctx->ptr;

    assert(ia >= 0 && ia < n);
    assert(ib >= 0 && ib < n);

    if (data[ia] < data[ib]) {
        return (-1);
    }

    else if (data[ia] > data[ib]) {
        return (1);
    }
    return (0);
}

int main() {
    float data[] = {2.5, 1.2, 3.7, 0.4, 1.2, 8.3};
    size_t n = sizeof(data)/sizeof(data[0]);
    int idx[n];
    ctx_t ctx;

    ctx.n = n;
    ctx.ptr = data;

    for (size_t i = 0; i < n; ++i) {
        idx[i] = i;
    }

    qsort_r(idx, n, sizeof(idx[0]), cmp_r, &ctx);

    for (size_t i = 1; i < n; ++i) {
        assert(data[idx[i]] >= data[idx[i-1]]);
    }
    for (size_t i = 0; i < n; ++i) {
        printf(\" %d %f%c\", idx[i], data[idx[i]], 10);
    }

    return (0);
}" QSORT_OK)

if(NOT QSORT_OK)
    message(FATAL_ERROR "Your C library does not support 'qsort_r'")
endif()
