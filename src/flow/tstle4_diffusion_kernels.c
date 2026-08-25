#include "tstle4_kernels_shared.h"

#pragma omp begin declare target


void tstle4_diff_u_c(const intk_c kk, const intk_c jj, const intk_c ii,
    realk_c *uo, const realk_c *u, const realk_c *v, const realk_c *w,
    const realk_c *g, const realk_c *dx, const realk_c *dy, const realk_c *dz,
    const realk_c *ddx, const realk_c *ddy, const realk_c *ddz,
    const intk_c nfro, const intk_c nbac, const intk_c ilesmodel,
    const realk_c gmol, const realk_c rho) {
    intk_c nfu = 0;
    intk_c nbu = 0;
    intk_c iles = 1;

    (void)dz;
    (void)ddx;

    if (nbac == 7)
        nbu = 1;
    if (nfro == 3)
        nfu = 1;
    if (nbac == 3)
        nbu = 1;

    if (ilesmodel == 0)
        iles = 0;

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
                const realk_c gn = g_kji * g_kj1i / max2(g_kji + g_kj1i, gmol) +
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
                const realk_c qe = -ge * ax * ((realk_c)1.0 / ddx[idx1(i + 1)]) *
                    (u[idx3(kk, jj, k, j, i + 1)] - u_kji);
                const realk_c qw = -gw * ax * ((realk_c)1.0 / ddx[idx1(i)]) *
                    (u_kji - u[idx3(kk, jj, k, j, i - 1)]);
                const realk_c qn = -gn * ay * ((realk_c)1.0 / dy[idx1(j)]) *
                    (u[idx3(kk, jj, k, j + 1, i)] - u_kji);
                const realk_c qs = -gs * ay * ((realk_c)1.0 / dy[idx1(j - 1)]) *
                    (u_kji - u[idx3(kk, jj, k, j - 1, i)]);
                const realk_c qt = -gt * az * ((realk_c)1.0 / dz[idx1(k)]) *
                    (u[idx3(kk, jj, k + 1, j, i)] - u_kji);
                const realk_c qb = -gb * az * ((realk_c)1.0 / dz[idx1(k - 1)]) *
                    (u_kji - u[idx3(kk, jj, k - 1, j, i)]);

                const realk_c st =
                    ((ge * (u[idx3(kk, jj, k, j, i + 1)] - u_kji) *
                        ((realk_c)1.0 / ddx[idx1(i + 1)])) -
                        (gw * (u_kji - u[idx3(kk, jj, k, j, i - 1)]) *
                            ((realk_c)1.0 / ddx[idx1(i)]))) *
                        ax +
                    ((gn *
                        (v[idx3(kk, jj, k, j, i + 1)] -
                            v[idx3(kk, jj, k, j, i)])) -
                        (gs *
                            (v[idx3(kk, jj, k, j - 1, i + 1)] -
                                v[idx3(kk, jj, k, j - 1, i)]))) *
                        ddz[idx1(k)] +
                    ((gt *
                        (w[idx3(kk, jj, k, j, i + 1)] -
                            w[idx3(kk, jj, k, j, i)])) -
                        (gb *
                            (w[idx3(kk, jj, k - 1, j, i + 1)] -
                                w[idx3(kk, jj, k - 1, j, i)]))) *
                        ddy[idx1(j)];
                const realk_c qc = st * (realk_c)iles;

                const realk_c fak = ((realk_c)1.0 / rho) *
                    ((realk_c)1.0 / ddy[idx1(j)]) *
                    ((realk_c)1.0 / dx[idx1(i)]) *
                    ((realk_c)1.0 / ddz[idx1(k)]);
                uo[idx3(kk, jj, k, j, i)] -=
                    fak * (qe - qw + qn - qs + qt - qb - qc);
            }
        }
    }
}

// Diffusive contribution for v-momentum with optional LES correction.
void tstle4_diff_v_c(const intk_c kk, const intk_c jj, const intk_c ii,
    realk_c *vo, const realk_c *u, const realk_c *v, const realk_c *w,
    const realk_c *g, const realk_c *dx, const realk_c *dy, const realk_c *dz,
    const realk_c *ddx, const realk_c *ddy, const realk_c *ddz,
    const intk_c nrgt, const intk_c nlft, const intk_c ilesmodel,
    const realk_c gmol, const realk_c rho) {
    intk_c nrv = 0;
    intk_c nlv = 0;
    intk_c iles = 1;

    (void)dx;
    (void)dz;

    if (nlft == 7)
        nlv = 1;
    if (nrgt == 3)
        nrv = 1;
    if (nlft == 3)
        nlv = 1;

    if (ilesmodel == 0)
        iles = 0;

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
                const realk_c qe = -ge * ax * ((realk_c)1.0 / dx[idx1(i)]) *
                    (v[idx3(kk, jj, k, j, i + 1)] - v_kji);
                const realk_c qw = -gw * ax * ((realk_c)1.0 / dx[idx1(i - 1)]) *
                    (v_kji - v[idx3(kk, jj, k, j, i - 1)]);
                const realk_c qn = -gn * ay * ((realk_c)1.0 / ddy[idx1(j + 1)]) *
                    (v[idx3(kk, jj, k, j + 1, i)] - v_kji);
                const realk_c qs = -gs * ay * ((realk_c)1.0 / ddy[idx1(j)]) *
                    (v_kji - v[idx3(kk, jj, k, j - 1, i)]);
                const realk_c qt = -gt * az * ((realk_c)1.0 / dz[idx1(k)]) *
                    (v[idx3(kk, jj, k + 1, j, i)] - v_kji);
                const realk_c qb = -gb * az * ((realk_c)1.0 / dz[idx1(k - 1)]) *
                    (v_kji - v[idx3(kk, jj, k - 1, j, i)]);

                const realk_c st = ((ge *
                                        (u[idx3(kk, jj, k, j + 1, i)] -
                                            u[idx3(kk, jj, k, j, i)])) -
                                    (gw *
                                        (u[idx3(kk, jj, k, j + 1, i - 1)] -
                                            u[idx3(kk, jj, k, j, i - 1)]))) *
                        ddz[idx1(k)] +
                    ((gn * (v[idx3(kk, jj, k, j + 1, i)] - v_kji) *
                        ((realk_c)1.0 / ddy[idx1(j + 1)])) -
                        (gs * (v_kji - v[idx3(kk, jj, k, j - 1, i)]) *
                            ((realk_c)1.0 / ddy[idx1(j)]))) *
                        ay +
                    ((gt *
                        (w[idx3(kk, jj, k, j + 1, i)] -
                            w[idx3(kk, jj, k, j, i)])) -
                        (gb *
                            (w[idx3(kk, jj, k - 1, j + 1, i)] -
                                w[idx3(kk, jj, k - 1, j, i)]))) *
                        ddx[idx1(i)];
                const realk_c qc = st * (realk_c)iles;

                const realk_c fak = ((realk_c)1.0 / rho) *
                    ((realk_c)1.0 / ddx[idx1(i)]) *
                    ((realk_c)1.0 / dy[idx1(j)]) *
                    ((realk_c)1.0 / ddz[idx1(k)]);
                vo[idx3(kk, jj, k, j, i)] -=
                    fak * (qe - qw + qn - qs + qt - qb - qc);
            }
        }
    }
}

// Diffusive contribution for w-momentum with optional LES correction.
void tstle4_diff_w_c(const intk_c kk, const intk_c jj, const intk_c ii,
    realk_c *wo, const realk_c *u, const realk_c *v, const realk_c *w,
    const realk_c *g, const realk_c *dx, const realk_c *dy, const realk_c *dz,
    const realk_c *ddx, const realk_c *ddy, const realk_c *ddz,
    const intk_c nbot, const intk_c ntop, const intk_c ilesmodel,
    const realk_c gmol, const realk_c rho) {
    intk_c nbw = 0;
    intk_c ntw = 0;
    intk_c iles = 1;

    (void)dx;
    (void)dy;
    (void)ddz;

    if (ntop == 7)
        ntw = 1;
    if (nbot == 3)
        nbw = 1;
    if (ntop == 3)
        ntw = 1;

    if (ilesmodel == 0)
        iles = 0;

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
                    g[idx3(kk, jj, k + 1, j, i)] *
                        g[idx3(kk, jj, k + 1, j, i + 1)] /
                        max2(g[idx3(kk, jj, k + 1, j, i)] +
                                g[idx3(kk, jj, k + 1, j, i + 1)],
                            gmol);

                const realk_c gw = g[idx3(kk, jj, k, j, i - 1)] * g_kji /
                        max2(g[idx3(kk, jj, k, j, i - 1)] + g_kji, gmol) +
                    g[idx3(kk, jj, k + 1, j, i - 1)] *
                        g[idx3(kk, jj, k + 1, j, i)] /
                        max2(g[idx3(kk, jj, k + 1, j, i - 1)] +
                                g[idx3(kk, jj, k + 1, j, i)],
                            gmol);

                const realk_c gn = g_kji * g[idx3(kk, jj, k, j + 1, i)] /
                        max2(g_kji + g[idx3(kk, jj, k, j + 1, i)], gmol) +
                    g[idx3(kk, jj, k + 1, j, i)] *
                        g[idx3(kk, jj, k + 1, j + 1, i)] /
                        max2(g[idx3(kk, jj, k + 1, j, i)] +
                                g[idx3(kk, jj, k + 1, j + 1, i)],
                            gmol);

                const realk_c gs = g[idx3(kk, jj, k, j - 1, i)] * g_kji /
                        max2(g[idx3(kk, jj, k, j - 1, i)] + g_kji, gmol) +
                    g[idx3(kk, jj, k + 1, j - 1, i)] *
                        g[idx3(kk, jj, k + 1, j, i)] /
                        max2(g[idx3(kk, jj, k + 1, j - 1, i)] +
                                g[idx3(kk, jj, k + 1, j, i)],
                            gmol);

                const realk_c gt = g[idx3(kk, jj, k + 1, j, i)];
                const realk_c gb = g_kji;

                const realk_c w_kji = w[idx3(kk, jj, k, j, i)];
                const realk_c qe = -ge * ax * ((realk_c)1.0 / dx[idx1(i)]) *
                    (w[idx3(kk, jj, k, j, i + 1)] - w_kji);
                const realk_c qw = -gw * ax * ((realk_c)1.0 / dx[idx1(i - 1)]) *
                    (w_kji - w[idx3(kk, jj, k, j, i - 1)]);
                const realk_c qn = -gn * ay * ((realk_c)1.0 / dy[idx1(j)]) *
                    (w[idx3(kk, jj, k, j + 1, i)] - w_kji);
                const realk_c qs = -gs * ay * ((realk_c)1.0 / dy[idx1(j - 1)]) *
                    (w_kji - w[idx3(kk, jj, k, j - 1, i)]);
                const realk_c qt = -gt * az * ((realk_c)1.0 / ddz[idx1(k + 1)]) *
                    (w[idx3(kk, jj, k + 1, j, i)] - w_kji);
                const realk_c qb = -gb * az * ((realk_c)1.0 / ddz[idx1(k)]) *
                    (w_kji - w[idx3(kk, jj, k - 1, j, i)]);

                const realk_c st = ((ge *
                                        (u[idx3(kk, jj, k + 1, j, i)] -
                                            u[idx3(kk, jj, k, j, i)])) -
                                    (gw *
                                        (u[idx3(kk, jj, k + 1, j, i - 1)] -
                                            u[idx3(kk, jj, k, j, i - 1)]))) *
                        ddy[idx1(j)] +
                    ((gn *
                        (v[idx3(kk, jj, k + 1, j, i)] -
                            v[idx3(kk, jj, k, j, i)])) -
                        (gs *
                            (v[idx3(kk, jj, k + 1, j - 1, i)] -
                                v[idx3(kk, jj, k, j - 1, i)]))) *
                        ddx[idx1(i)] +
                    ((gt * (w[idx3(kk, jj, k + 1, j, i)] - w_kji) *
                        ((realk_c)1.0 / ddz[idx1(k + 1)])) -
                        (gb * (w_kji - w[idx3(kk, jj, k - 1, j, i)]) *
                            ((realk_c)1.0 / ddz[idx1(k)]))) *
                        az;
                const realk_c qc = st * (realk_c)iles;

                const realk_c fak = ((realk_c)1.0 / rho) *
                    ((realk_c)1.0 / ddx[idx1(i)]) *
                    ((realk_c)1.0 / ddy[idx1(j)]) *
                    ((realk_c)1.0 / dz[idx1(k)]);
                wo[idx3(kk, jj, k, j, i)] -=
                    fak * (qe - qw + qn - qs + qt - qb - qc);
            }
        }
    }
}

// Slip-wall correction on boundary planes using the local wall model.
void tstle4_diff_swc_c(const intk_c kk, const intk_c jj, const intk_c ii,
    realk_c *uo, realk_c *vo, realk_c *wo, const realk_c *u, const realk_c *v,
    const realk_c *w, const realk_c *ddx, const realk_c *ddy,
    const realk_c *ddz, const intk_c nfro, const intk_c nbac, const intk_c nrgt,
    const intk_c nlft, const intk_c nbot, const intk_c ntop, const realk_c gmol,
    const realk_c rho) {
    if (nfro == 5) {
        const intk_c i = 3;
#if defined(_OPENMP)
#pragma omp for collapse(2)
#endif
        for (intk_c j = 2; j <= jj - 1; ++j) {
            for (intk_c k = 2; k <= kk - 1; ++k) {
                const realk_c d = ddx[idx1(i)];
                vo[idx3(kk, jj, k, j, i)] -=
                    swcle3d_one_c(d, v[idx3(kk, jj, k, j, i)], gmol, rho);
                wo[idx3(kk, jj, k, j, i)] -=
                    swcle3d_one_c(d, w[idx3(kk, jj, k, j, i)], gmol, rho);
            }
        }
#if defined(_OPENMP)
#pragma omp barrier
#endif
    }

    if (nbac == 5) {
        const intk_c i = ii - 2;
#if defined(_OPENMP)
#pragma omp for collapse(2)
#endif
        for (intk_c j = 2; j <= jj - 1; ++j) {
            for (intk_c k = 2; k <= kk - 1; ++k) {
                const realk_c d = ddx[idx1(i)];
                vo[idx3(kk, jj, k, j, i)] -=
                    swcle3d_one_c(d, v[idx3(kk, jj, k, j, i)], gmol, rho);
                wo[idx3(kk, jj, k, j, i)] -=
                    swcle3d_one_c(d, w[idx3(kk, jj, k, j, i)], gmol, rho);
            }
        }
#if defined(_OPENMP)
#pragma omp barrier
#endif
    }

    if (nrgt == 5) {
        const intk_c j = 3;
#if defined(_OPENMP)
#pragma omp for collapse(2)
#endif
        for (intk_c i = 2; i <= ii - 1; ++i) {
            for (intk_c k = 2; k <= kk - 1; ++k) {
                const realk_c d = ddy[idx1(j)];
                uo[idx3(kk, jj, k, j, i)] -=
                    swcle3d_one_c(d, u[idx3(kk, jj, k, j, i)], gmol, rho);
                wo[idx3(kk, jj, k, j, i)] -=
                    swcle3d_one_c(d, w[idx3(kk, jj, k, j, i)], gmol, rho);
            }
        }
#if defined(_OPENMP)
#pragma omp barrier
#endif
    }

    if (nlft == 5) {
        const intk_c j = jj - 2;
#if defined(_OPENMP)
#pragma omp for collapse(2)
#endif
        for (intk_c i = 2; i <= ii - 1; ++i) {
            for (intk_c k = 2; k <= kk - 1; ++k) {
                const realk_c d = ddy[idx1(j)];
                uo[idx3(kk, jj, k, j, i)] -=
                    swcle3d_one_c(d, u[idx3(kk, jj, k, j, i)], gmol, rho);
                wo[idx3(kk, jj, k, j, i)] -=
                    swcle3d_one_c(d, w[idx3(kk, jj, k, j, i)], gmol, rho);
            }
        }
#if defined(_OPENMP)
#pragma omp barrier
#endif
    }

    if (nbot == 5) {
        const intk_c k = 3;
#if defined(_OPENMP)
#pragma omp for collapse(2)
#endif
        for (intk_c i = 2; i <= ii - 1; ++i) {
            for (intk_c j = 2; j <= jj - 1; ++j) {
                const realk_c d = ddz[idx1(k)];
                uo[idx3(kk, jj, k, j, i)] -=
                    swcle3d_one_c(d, u[idx3(kk, jj, k, j, i)], gmol, rho);
                vo[idx3(kk, jj, k, j, i)] -=
                    swcle3d_one_c(d, v[idx3(kk, jj, k, j, i)], gmol, rho);
            }
        }
#if defined(_OPENMP)
#pragma omp barrier
#endif
    }

    if (ntop == 5) {
        const intk_c k = kk - 2;
#if defined(_OPENMP)
#pragma omp for collapse(2)
#endif
        for (intk_c i = 2; i <= ii - 1; ++i) {
            for (intk_c j = 2; j <= jj - 1; ++j) {
                const realk_c d = ddz[idx1(k)];
                uo[idx3(kk, jj, k, j, i)] -=
                    swcle3d_one_c(d, u[idx3(kk, jj, k, j, i)], gmol, rho);
                vo[idx3(kk, jj, k, j, i)] -=
                    swcle3d_one_c(d, v[idx3(kk, jj, k, j, i)], gmol, rho);
            }
        }
    }
}

#pragma omp end declare target



// Full diffusive implementation: iterate all grids in C using Fortran metadata.
void tstle4_diff_impl_c(realk_c *uo, realk_c *vo, realk_c *wo,
    const realk_c *u, const realk_c *v, const realk_c *w, const realk_c *g,
    const realk_c *dx, const realk_c *dy, const realk_c *dz,
    const realk_c *ddx, const realk_c *ddy, const realk_c *ddz,
    const intk_c ilesmodel_in, const realk_c gmol_in, const realk_c rho_in,
    const intk_c nmygrids_in, const tstle4_grid_t *grids)
{
    #pragma omp target teams distribute
    for (intk_c igrid = 0; igrid < nmygrids_in; ++igrid)
    {
        const tstle4_grid_t *const gr = &grids[igrid];

        realk_c *const uo_g = &uo[(size_t)(gr->ip3 - 1)];
        realk_c *const vo_g = &vo[(size_t)(gr->ip3 - 1)];
        realk_c *const wo_g = &wo[(size_t)(gr->ip3 - 1)];

        const realk_c *const u_g = &u[(size_t)(gr->ip3 - 1)];
        const realk_c *const v_g = &v[(size_t)(gr->ip3 - 1)];
        const realk_c *const w_g = &w[(size_t)(gr->ip3 - 1)];
        const realk_c *const g_g = &g[(size_t)(gr->ip3 - 1)];

        const realk_c *const dx_g = &dx[(size_t)(gr->ipx - 1)];
        const realk_c *const dy_g = &dy[(size_t)(gr->ipy - 1)];
        const realk_c *const dz_g = &dz[(size_t)(gr->ipz - 1)];
        const realk_c *const ddx_g = &ddx[(size_t)(gr->ipx - 1)];
        const realk_c *const ddy_g = &ddy[(size_t)(gr->ipy - 1)];
        const realk_c *const ddz_g = &ddz[(size_t)(gr->ipz - 1)];

        #pragma omp parallel
        {
            tstle4_diff_swc_c(gr->kk, gr->jj, gr->ii, uo_g, vo_g, wo_g,
                u_g, v_g, w_g, ddx_g, ddy_g, ddz_g, gr->nfro, gr->nbac,
                gr->nrgt, gr->nlft, gr->nbot, gr->ntop, gmol_in, rho_in);

            tstle4_diff_u_c(gr->kk, gr->jj, gr->ii, uo_g, u_g, v_g, w_g, g_g,
                dx_g, dy_g, dz_g, ddx_g, ddy_g, ddz_g,
                gr->nfro, gr->nbac, ilesmodel_in, gmol_in, rho_in);

            tstle4_diff_v_c(gr->kk, gr->jj, gr->ii, vo_g, u_g, v_g, w_g, g_g,
                dx_g, dy_g, dz_g, ddx_g, ddy_g, ddz_g,
                gr->nrgt, gr->nlft, ilesmodel_in, gmol_in, rho_in);

            tstle4_diff_w_c(gr->kk, gr->jj, gr->ii, wo_g, u_g, v_g, w_g, g_g,
                dx_g, dy_g, dz_g, ddx_g, ddy_g, ddz_g,
                gr->nbot, gr->ntop, ilesmodel_in, gmol_in, rho_in);
        }
    }
}
