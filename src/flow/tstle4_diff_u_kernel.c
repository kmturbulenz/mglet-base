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

static inline size_t idx3(const intk_c kk, const intk_c jj, const intk_c k,
        const intk_c j, const intk_c i)
{
    return (size_t)(k - 1) + (size_t)kk * ((size_t)(j - 1) + (size_t)jj * (size_t)(i - 1));
}

static inline size_t idx1(const intk_c i)
{
    return (size_t)(i - 1);
}

static inline realk_c max2(const realk_c a, const realk_c b)
{
    return (a > b) ? a : b;
}

#if defined(_OPENMP)
#pragma omp declare target
#endif
void tstle4_diff_u_c(const intk_c kk, const intk_c jj, const intk_c ii,
        realk_c *uo,
        const realk_c *u, const realk_c *v, const realk_c *w,
        const realk_c *g,
        const realk_c *dx, const realk_c *dy, const realk_c *dz,
        const realk_c *ddx, const realk_c *ddy, const realk_c *ddz,
        const realk_c *rdx, const realk_c *rdy, const realk_c *rdz,
        const realk_c *rddx, const realk_c *rddy, const realk_c *rddz,
        const intk_c nfro, const intk_c nbac,
        const intk_c ilesmodel, const realk_c gmol, const realk_c rho)
{
    intk_c nfu = 0;
    intk_c nbu = 0;
    intk_c iles = 1;

    (void)dz;
    (void)ddx;

    if (nbac == 7) nbu = 1;
    if (nfro == 3) nfu = 1;
    if (nbac == 3) nbu = 1;

    if (ilesmodel == 0) iles = 0;

#if defined(_OPENMP)
#pragma omp for collapse(3)
#endif
    for (intk_c i = 3 - nfu; i <= ii - 3 + nbu; ++i) {
        for (intk_c j = 3; j <= jj - 2; ++j) {
            for (intk_c k = 3; k <= kk - 2; ++k) {
                const realk_c ax = ddy[idx1(j)] * ddz[idx1(k)];
                const realk_c ay = dx[idx1(i)] * ddz[idx1(k)];
                const realk_c az = dx[idx1(i)] * ddy[idx1(j)];

                const realk_c g_kji = g[idx3(kk, jj, k, j, i)];
                const realk_c ge = g[idx3(kk, jj, k, j, i + 1)];
                const realk_c gw = g_kji;

                const realk_c g_kj1i = g[idx3(kk, jj, k, j + 1, i)];
                const realk_c g_kj1ip1 = g[idx3(kk, jj, k, j + 1, i + 1)];
                const realk_c gn =
                    g_kji * g_kj1i / max2(g_kji + g_kj1i, gmol) +
                    ge * g_kj1ip1 / max2(ge + g_kj1ip1, gmol);

                const realk_c g_kjm1i = g[idx3(kk, jj, k, j - 1, i)];
                const realk_c g_kjm1ip1 = g[idx3(kk, jj, k, j - 1, i + 1)];
                const realk_c gs =
                    g_kjm1i * g_kji / max2(g_kjm1i + g_kji, gmol) +
                    g_kjm1ip1 * ge / max2(g_kjm1ip1 + ge, gmol);

                const realk_c g_kp1ji = g[idx3(kk, jj, k + 1, j, i)];
                const realk_c g_kp1jip1 = g[idx3(kk, jj, k + 1, j, i + 1)];
                const realk_c gt =
                    g_kji * g_kp1ji / max2(g_kji + g_kp1ji, gmol) +
                    ge * g_kp1jip1 / max2(ge + g_kp1jip1, gmol);

                const realk_c g_km1ji = g[idx3(kk, jj, k - 1, j, i)];
                const realk_c g_km1jip1 = g[idx3(kk, jj, k - 1, j, i + 1)];
                const realk_c gb =
                    g_km1ji * g_kji / max2(g_km1ji + g_kji, gmol) +
                    g_km1jip1 * ge / max2(g_km1jip1 + ge, gmol);

                const realk_c u_kji = u[idx3(kk, jj, k, j, i)];
                const realk_c qe = -ge * ax * rddx[idx1(i + 1)] * (u[idx3(kk, jj, k, j, i + 1)] - u_kji);
                const realk_c qw = -gw * ax * rddx[idx1(i)] * (u_kji - u[idx3(kk, jj, k, j, i - 1)]);
                const realk_c qn = -gn * ay * rdy[idx1(j)] * (u[idx3(kk, jj, k, j + 1, i)] - u_kji);
                const realk_c qs = -gs * ay * rdy[idx1(j - 1)] * (u_kji - u[idx3(kk, jj, k, j - 1, i)]);
                const realk_c qt = -gt * az * rdz[idx1(k)] * (u[idx3(kk, jj, k + 1, j, i)] - u_kji);
                const realk_c qb = -gb * az * rdz[idx1(k - 1)] * (u_kji - u[idx3(kk, jj, k - 1, j, i)]);

                const realk_c st =
                    ((ge * (u[idx3(kk, jj, k, j, i + 1)] - u_kji) * rddx[idx1(i + 1)]) -
                    (gw * (u_kji - u[idx3(kk, jj, k, j, i - 1)]) * rddx[idx1(i)])) * ax +
                    ((gn * (v[idx3(kk, jj, k, j, i + 1)] - v[idx3(kk, jj, k, j, i)])) -
                    (gs * (v[idx3(kk, jj, k, j - 1, i + 1)] - v[idx3(kk, jj, k, j - 1, i)]))) * ddz[idx1(k)] +
                    ((gt * (w[idx3(kk, jj, k, j, i + 1)] - w[idx3(kk, jj, k, j, i)])) -
                    (gb * (w[idx3(kk, jj, k - 1, j, i + 1)] - w[idx3(kk, jj, k - 1, j, i)]))) * ddy[idx1(j)];
                const realk_c qc = st * (realk_c)iles;

                const realk_c fak = ((realk_c)1.0 / rho) * rddy[idx1(j)] * rdx[idx1(i)] * rddz[idx1(k)];
                uo[idx3(kk, jj, k, j, i)] -= fak * (qe - qw + qn - qs + qt - qb - qc);
            }
        }
    }
}

void tstle4_diff_v_c(const intk_c kk, const intk_c jj, const intk_c ii,
        realk_c *vo,
        const realk_c *u, const realk_c *v, const realk_c *w,
        const realk_c *g,
        const realk_c *dx, const realk_c *dy, const realk_c *dz,
        const realk_c *ddx, const realk_c *ddy, const realk_c *ddz,
        const realk_c *rdx, const realk_c *rdy, const realk_c *rdz,
        const realk_c *rddx, const realk_c *rddy, const realk_c *rddz,
        const intk_c nrgt, const intk_c nlft,
        const intk_c ilesmodel, const realk_c gmol, const realk_c rho)
{
    intk_c nrv = 0;
    intk_c nlv = 0;
    intk_c iles = 1;

    (void)dx;
    (void)dz;

    if (nlft == 7) nlv = 1;
    if (nrgt == 3) nrv = 1;
    if (nlft == 3) nlv = 1;

    if (ilesmodel == 0) iles = 0;

#if defined(_OPENMP)
#pragma omp for collapse(3)
#endif
    for (intk_c i = 3; i <= ii - 2; ++i) {
        for (intk_c j = 3 - nrv; j <= jj - 3 + nlv; ++j) {
            for (intk_c k = 3; k <= kk - 2; ++k) {
                const realk_c ax = dy[idx1(j)] * ddz[idx1(k)];
                const realk_c ay = ddx[idx1(i)] * ddz[idx1(k)];
                const realk_c az = ddx[idx1(i)] * dy[idx1(j)];

                const realk_c g_kji = g[idx3(kk, jj, k, j, i)];
                const realk_c g_kjip1 = g[idx3(kk, jj, k, j, i + 1)];
                const realk_c g_kj1i = g[idx3(kk, jj, k, j + 1, i)];
                const realk_c g_kj1ip1 = g[idx3(kk, jj, k, j + 1, i + 1)];
                const realk_c ge =
                    g_kji * g_kjip1 / max2(g_kji + g_kjip1, gmol) +
                    g_kj1i * g_kj1ip1 / max2(g_kj1i + g_kj1ip1, gmol);

                const realk_c g_kjim1 = g[idx3(kk, jj, k, j, i - 1)];
                const realk_c g_kj1im1 = g[idx3(kk, jj, k, j + 1, i - 1)];
                const realk_c gw =
                    g_kjim1 * g_kji / max2(g_kjim1 + g_kji, gmol) +
                    g_kj1im1 * g_kj1i / max2(g_kj1im1 + g_kj1i, gmol);

                const realk_c gn = g_kj1i;
                const realk_c gs = g_kji;

                const realk_c g_kp1ji = g[idx3(kk, jj, k + 1, j, i)];
                const realk_c g_kp1j1i = g[idx3(kk, jj, k + 1, j + 1, i)];
                const realk_c gt =
                    g_kji * g_kp1ji / max2(g_kji + g_kp1ji, gmol) +
                    g_kj1i * g_kp1j1i / max2(g_kj1i + g_kp1j1i, gmol);

                const realk_c g_km1ji = g[idx3(kk, jj, k - 1, j, i)];
                const realk_c g_km1j1i = g[idx3(kk, jj, k - 1, j + 1, i)];
                const realk_c gb =
                    g_km1ji * g_kji / max2(g_km1ji + g_kji, gmol) +
                    g_km1j1i * g_kj1i / max2(g_km1j1i + g_kj1i, gmol);

                const realk_c v_kji = v[idx3(kk, jj, k, j, i)];
                const realk_c qe = -ge * ax * rdx[idx1(i)] * (v[idx3(kk, jj, k, j, i + 1)] - v_kji);
                const realk_c qw = -gw * ax * rdx[idx1(i - 1)] * (v_kji - v[idx3(kk, jj, k, j, i - 1)]);
                const realk_c qn = -gn * ay * rddy[idx1(j + 1)] * (v[idx3(kk, jj, k, j + 1, i)] - v_kji);
                const realk_c qs = -gs * ay * rddy[idx1(j)] * (v_kji - v[idx3(kk, jj, k, j - 1, i)]);
                const realk_c qt = -gt * az * rdz[idx1(k)] * (v[idx3(kk, jj, k + 1, j, i)] - v_kji);
                const realk_c qb = -gb * az * rdz[idx1(k - 1)] * (v_kji - v[idx3(kk, jj, k - 1, j, i)]);

                const realk_c st =
                    ((ge * (u[idx3(kk, jj, k, j + 1, i)] - u[idx3(kk, jj, k, j, i)])) -
                    (gw * (u[idx3(kk, jj, k, j + 1, i - 1)] - u[idx3(kk, jj, k, j, i - 1)]))) * ddz[idx1(k)] +
                    ((gn * (v[idx3(kk, jj, k, j + 1, i)] - v_kji) * rddy[idx1(j + 1)]) -
                    (gs * (v_kji - v[idx3(kk, jj, k, j - 1, i)]) * rddy[idx1(j)])) * ay +
                    ((gt * (w[idx3(kk, jj, k, j + 1, i)] - w[idx3(kk, jj, k, j, i)])) -
                    (gb * (w[idx3(kk, jj, k - 1, j + 1, i)] - w[idx3(kk, jj, k - 1, j, i)]))) * ddx[idx1(i)];
                const realk_c qc = st * (realk_c)iles;

                const realk_c fak = ((realk_c)1.0 / rho) * rddx[idx1(i)] * rdy[idx1(j)] * rddz[idx1(k)];
                vo[idx3(kk, jj, k, j, i)] -= fak * (qe - qw + qn - qs + qt - qb - qc);
            }
        }
    }
}

void tstle4_diff_w_c(const intk_c kk, const intk_c jj, const intk_c ii,
        realk_c *wo,
        const realk_c *u, const realk_c *v, const realk_c *w,
        const realk_c *g,
        const realk_c *dx, const realk_c *dy, const realk_c *dz,
        const realk_c *ddx, const realk_c *ddy, const realk_c *ddz,
        const realk_c *rdx, const realk_c *rdy, const realk_c *rdz,
        const realk_c *rddx, const realk_c *rddy, const realk_c *rddz,
        const intk_c nbot, const intk_c ntop,
        const intk_c ilesmodel, const realk_c gmol, const realk_c rho)
{
    intk_c nbw = 0;
    intk_c ntw = 0;
    intk_c iles = 1;

    (void)dx;
    (void)dy;
    (void)ddz;

    if (ntop == 7) ntw = 1;
    if (nbot == 3) nbw = 1;
    if (ntop == 3) ntw = 1;

    if (ilesmodel == 0) iles = 0;

#if defined(_OPENMP)
    #pragma omp for collapse(3)
#endif
    for (intk_c i = 3; i <= ii - 2; ++i) {
        for (intk_c j = 3; j <= jj - 2; ++j) {
            for (intk_c k = 3 - nbw; k <= kk - 3 + ntw; ++k) {
                const realk_c ax = ddy[idx1(j)] * dz[idx1(k)];
                const realk_c ay = ddx[idx1(i)] * dz[idx1(k)];
                const realk_c az = ddx[idx1(i)] * ddy[idx1(j)];

                const realk_c g_kji = g[idx3(kk, jj, k, j, i)];
                const realk_c g_kjip1 = g[idx3(kk, jj, k, j, i + 1)];
                const realk_c ge =
                    g_kji * g_kjip1 / max2(g_kji + g_kjip1, gmol) +
                    g[idx3(kk, jj, k + 1, j, i)] * g[idx3(kk, jj, k + 1, j, i + 1)] /
                    max2(g[idx3(kk, jj, k + 1, j, i)] + g[idx3(kk, jj, k + 1, j, i + 1)], gmol);

                const realk_c gw =
                    g[idx3(kk, jj, k, j, i - 1)] * g_kji /
                    max2(g[idx3(kk, jj, k, j, i - 1)] + g_kji, gmol) +
                    g[idx3(kk, jj, k + 1, j, i - 1)] * g[idx3(kk, jj, k + 1, j, i)] /
                    max2(g[idx3(kk, jj, k + 1, j, i - 1)] + g[idx3(kk, jj, k + 1, j, i)], gmol);

                const realk_c gn =
                    g_kji * g[idx3(kk, jj, k, j + 1, i)] /
                    max2(g_kji + g[idx3(kk, jj, k, j + 1, i)], gmol) +
                    g[idx3(kk, jj, k + 1, j, i)] * g[idx3(kk, jj, k + 1, j + 1, i)] /
                    max2(g[idx3(kk, jj, k + 1, j, i)] + g[idx3(kk, jj, k + 1, j + 1, i)], gmol);

                const realk_c gs =
                    g[idx3(kk, jj, k, j - 1, i)] * g_kji /
                    max2(g[idx3(kk, jj, k, j - 1, i)] + g_kji, gmol) +
                    g[idx3(kk, jj, k + 1, j - 1, i)] * g[idx3(kk, jj, k + 1, j, i)] /
                    max2(g[idx3(kk, jj, k + 1, j - 1, i)] + g[idx3(kk, jj, k + 1, j, i)], gmol);

                const realk_c gt = g[idx3(kk, jj, k + 1, j, i)];
                const realk_c gb = g_kji;

                const realk_c w_kji = w[idx3(kk, jj, k, j, i)];
                const realk_c qe = -ge * ax * rdx[idx1(i)] * (w[idx3(kk, jj, k, j, i + 1)] - w_kji);
                const realk_c qw = -gw * ax * rdx[idx1(i - 1)] * (w_kji - w[idx3(kk, jj, k, j, i - 1)]);
                const realk_c qn = -gn * ay * rdy[idx1(j)] * (w[idx3(kk, jj, k, j + 1, i)] - w_kji);
                const realk_c qs = -gs * ay * rdy[idx1(j - 1)] * (w_kji - w[idx3(kk, jj, k, j - 1, i)]);
                const realk_c qt = -gt * az * rddz[idx1(k + 1)] * (w[idx3(kk, jj, k + 1, j, i)] - w_kji);
                const realk_c qb = -gb * az * rddz[idx1(k)] * (w_kji - w[idx3(kk, jj, k - 1, j, i)]);

                const realk_c st =
                    ((ge * (u[idx3(kk, jj, k + 1, j, i)] - u[idx3(kk, jj, k, j, i)])) -
                    (gw * (u[idx3(kk, jj, k + 1, j, i - 1)] - u[idx3(kk, jj, k, j, i - 1)]))) * ddy[idx1(j)] +
                    ((gn * (v[idx3(kk, jj, k + 1, j, i)] - v[idx3(kk, jj, k, j, i)])) -
                    (gs * (v[idx3(kk, jj, k + 1, j - 1, i)] - v[idx3(kk, jj, k, j - 1, i)]))) * ddx[idx1(i)] +
                    ((gt * (w[idx3(kk, jj, k + 1, j, i)] - w_kji) * rddz[idx1(k + 1)]) -
                    (gb * (w_kji - w[idx3(kk, jj, k - 1, j, i)]) * rddz[idx1(k)])) * az;
                const realk_c qc = st * (realk_c)iles;

                const realk_c fak = ((realk_c)1.0 / rho) * rddx[idx1(i)] * rddy[idx1(j)] * rdz[idx1(k)];
                wo[idx3(kk, jj, k, j, i)] -= fak * (qe - qw + qn - qs + qt - qb - qc);
            }
        }
    }
}
#if defined(_OPENMP)
    #pragma omp end declare target
#endif
