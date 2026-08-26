#include <float.h>
#include <stddef.h>

#ifdef _MGLET_DOUBLE_PRECISION_
typedef double realk_c;
#else
typedef float realk_c;
#endif

#ifdef _MGLET_INT64_
typedef long long intk_c;
#else
typedef int intk_c;
#endif

typedef struct {
    intk_c kk;
    intk_c jj;
    intk_c ii;
    intk_c iface;
    intk_c ityp;
    intk_c ip3;
    intk_c ipbb;
} boundg_task_t;

static inline size_t idx3(const intk_c kk, const intk_c jj, const intk_c k,
    const intk_c j, const intk_c i)
{
    return (size_t)(k - 1) +
        (size_t)kk * ((size_t)(j - 1) + (size_t)jj * (size_t)(i - 1));
}

static inline size_t idx2(const intk_c n1, const intk_c i1, const intk_c i2)
{
    return (size_t)(i1 - 1) + (size_t)n1 * (size_t)(i2 - 1);
}

static inline realk_c eps_one_c(void)
{
#ifdef _MGLET_DOUBLE_PRECISION_
    return (realk_c)DBL_EPSILON;
#else
    return (realk_c)FLT_EPSILON;
#endif
}

#if defined(_OPENMP)
#pragma omp begin declare target
#endif

static void boundg_nobuffer_device_c(const intk_c kk, const intk_c jj,
    const intk_c ii, const intk_c iface, const intk_c ityp, realk_c *g,
    const realk_c gmol)
{
    intk_c i2, i3, j2, j3, k2, k3;

    switch (iface) {
    case 1:
    case 2:
        if (iface == 1) {
            i2 = 2;
            i3 = 3;
        } else {
            i2 = ii - 1;
            i3 = ii - 2;
        }
#if defined(_OPENMP)
#pragma omp for collapse(2)
#endif
        for (intk_c j = 1; j <= jj; ++j) {
            for (intk_c k = 1; k <= kk; ++k) {
                const size_t ig = idx3(kk, jj, k, j, i2);
                if (ityp == 1) {
                    g[ig] = -eps_one_c() * gmol;
                } else {
                    g[ig] = g[idx3(kk, jj, k, j, i3)];
                }
            }
        }
        break;
    case 3:
    case 4:
        if (iface == 3) {
            j2 = 2;
            j3 = 3;
        } else {
            j2 = jj - 1;
            j3 = jj - 2;
        }
#if defined(_OPENMP)
#pragma omp for collapse(2)
#endif
        for (intk_c i = 1; i <= ii; ++i) {
            for (intk_c k = 1; k <= kk; ++k) {
                const size_t ig = idx3(kk, jj, k, j2, i);
                if (ityp == 1) {
                    g[ig] = -eps_one_c() * gmol;
                } else {
                    g[ig] = g[idx3(kk, jj, k, j3, i)];
                }
            }
        }
        break;
    case 5:
    case 6:
        if (iface == 5) {
            k2 = 2;
            k3 = 3;
        } else {
            k2 = kk - 1;
            k3 = kk - 2;
        }
#if defined(_OPENMP)
#pragma omp for collapse(2)
#endif
        for (intk_c i = 1; i <= ii; ++i) {
            for (intk_c j = 1; j <= jj; ++j) {
                const size_t ig = idx3(kk, jj, k2, j, i);
                if (ityp == 1) {
                    g[ig] = -eps_one_c() * gmol;
                } else {
                    g[ig] = g[idx3(kk, jj, k3, j, i)];
                }
            }
        }
        break;
    default:
        break;
    }
}

static void bfront_device_c(const intk_c kk, const intk_c jj, const intk_c ii,
    const intk_c i2, const intk_c i3, const intk_c ityp,
    const realk_c *buffer, realk_c *g, const realk_c *bp, const realk_c gmol)
{
    (void)ii;

    switch (ityp) {
    case 1:
#if defined(_OPENMP)
#pragma omp for collapse(2)
#endif
        for (intk_c j = 1; j <= jj; ++j) {
            for (intk_c k = 1; k <= kk; ++k) {
                g[idx3(kk, jj, k, j, i2)] = -eps_one_c() * gmol;
            }
        }
        break;
    case 2:
#if defined(_OPENMP)
#pragma omp for collapse(2)
#endif
        for (intk_c j = 1; j <= jj; ++j) {
            for (intk_c k = 1; k <= kk; ++k) {
                g[idx3(kk, jj, k, j, i2)] = g[idx3(kk, jj, k, j, i3)];
            }
        }
        break;
    case 3:
#if defined(_OPENMP)
#pragma omp for collapse(2)
#endif
        for (intk_c j = 2; j <= jj - 1; ++j) {
            for (intk_c k = 2; k <= kk - 1; ++k) {
                const size_t ig = idx3(kk, jj, k, j, i2);
                const realk_c sbp = bp[ig];
                g[ig] = buffer[idx2(kk, k, j)] * sbp + ((realk_c)1.0 - sbp) * g[ig];
            }
        }
        break;
    default:
        break;
    }
}

static void bright_device_c(const intk_c kk, const intk_c jj, const intk_c ii,
    const intk_c j2, const intk_c j3, const intk_c ityp,
    const realk_c *buffer, realk_c *g, const realk_c *bp, const realk_c gmol)
{
    switch (ityp) {
    case 1:
#if defined(_OPENMP)
#pragma omp for collapse(2)
#endif
        for (intk_c i = 1; i <= ii; ++i) {
            for (intk_c k = 1; k <= kk; ++k) {
                g[idx3(kk, jj, k, j2, i)] = -eps_one_c() * gmol;
            }
        }
        break;
    case 2:
#if defined(_OPENMP)
#pragma omp for collapse(2)
#endif
        for (intk_c i = 1; i <= ii; ++i) {
            for (intk_c k = 1; k <= kk; ++k) {
                g[idx3(kk, jj, k, j2, i)] = g[idx3(kk, jj, k, j3, i)];
            }
        }
        break;
    case 3:
#if defined(_OPENMP)
#pragma omp for collapse(2)
#endif
        for (intk_c i = 2; i <= ii - 1; ++i) {
            for (intk_c k = 2; k <= kk - 1; ++k) {
                const size_t ig = idx3(kk, jj, k, j2, i);
                const realk_c sbp = bp[ig];
                g[ig] = buffer[idx2(kk, k, i)] * sbp + ((realk_c)1.0 - sbp) * g[ig];
            }
        }
        break;
    default:
        break;
    }
}

static void bbottom_device_c(const intk_c kk, const intk_c jj,
    const intk_c ii, const intk_c k2, const intk_c k3, const intk_c ityp,
    const realk_c *buffer, realk_c *g, const realk_c *bp, const realk_c gmol)
{
    (void)kk;

    switch (ityp) {
    case 1:
#if defined(_OPENMP)
#pragma omp for collapse(2)
#endif
        for (intk_c i = 1; i <= ii; ++i) {
            for (intk_c j = 1; j <= jj; ++j) {
                g[idx3(kk, jj, k2, j, i)] = -eps_one_c() * gmol;
            }
        }
        break;
    case 2:
#if defined(_OPENMP)
#pragma omp for collapse(2)
#endif
        for (intk_c i = 1; i <= ii; ++i) {
            for (intk_c j = 1; j <= jj; ++j) {
                g[idx3(kk, jj, k2, j, i)] = g[idx3(kk, jj, k3, j, i)];
            }
        }
        break;
    case 3:
#if defined(_OPENMP)
#pragma omp for collapse(2)
#endif
        for (intk_c i = 2; i <= ii - 1; ++i) {
            for (intk_c j = 2; j <= jj - 1; ++j) {
                const size_t ig = idx3(kk, jj, k2, j, i);
                const realk_c sbp = bp[ig];
                g[ig] = buffer[idx2(jj, j, i)] * sbp + ((realk_c)1.0 - sbp) * g[ig];
            }
        }
        break;
    default:
        break;
    }
}

#if defined(_OPENMP)
#pragma omp end declare target
#endif


void apply_boundg_impl_c(const intk_c ntasks, realk_c *g, realk_c *gbuffer,
    const realk_c *bp, const boundg_task_t *tasks, const realk_c gmol_in)
{
#pragma omp target teams distribute
    for (intk_c itask = 0; itask < ntasks; ++itask) {
        const boundg_task_t *const task = &tasks[itask];
        realk_c *const g_g = &g[(size_t)(task->ip3 - 1)];
        const realk_c *const bp_g = &bp[(size_t)(task->ip3 - 1)];

        if (task->ityp != 3) {
#pragma omp parallel
            {
                boundg_nobuffer_device_c(task->kk, task->jj, task->ii,
                    task->iface, task->ityp, g_g, gmol_in);
            }
            continue;
        }

        const realk_c *const buffer_g = &gbuffer[(size_t)(task->ipbb - 1)];
#pragma omp parallel
        {
            switch (task->iface) {
            case 1:
                bfront_device_c(task->kk, task->jj, task->ii, 2, 3,
                    task->ityp, buffer_g, g_g, bp_g, gmol_in);
                break;
            case 2:
                bfront_device_c(task->kk, task->jj, task->ii,
                    task->ii - 1, task->ii - 2, task->ityp, buffer_g, g_g,
                    bp_g, gmol_in);
                break;
            case 3:
                bright_device_c(task->kk, task->jj, task->ii, 2, 3,
                    task->ityp, buffer_g, g_g, bp_g, gmol_in);
                break;
            case 4:
                bright_device_c(task->kk, task->jj, task->ii,
                    task->jj - 1, task->jj - 2, task->ityp, buffer_g, g_g,
                    bp_g, gmol_in);
                break;
            case 5:
                bbottom_device_c(task->kk, task->jj, task->ii, 2, 3,
                    task->ityp, buffer_g, g_g, bp_g, gmol_in);
                break;
            case 6:
                bbottom_device_c(task->kk, task->jj, task->ii,
                    task->kk - 1, task->kk - 2, task->ityp, buffer_g, g_g,
                    bp_g, gmol_in);
                break;
            default:
                break;
            }
        }
    }
}
