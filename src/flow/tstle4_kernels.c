#include <math.h>
#include <stddef.h>
#include <stdio.h>

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

// Data structure to hold grid metadata for TSTLE4 computations.
typedef struct {
    // Grid dimensions
    intk_c ii; intk_c jj; intk_c kk;
    // Grid offsets in the global arrays
    intk_c ip3; intk_c ipx; intk_c ipy; intk_c ipz;
    // Grid boundary conditions
    intk_c nfro; intk_c nbac; intk_c nrgt; intk_c nlft; intk_c nbot; intk_c ntop;
    // Grid pressure gradient source terms
    realk_c gpx; realk_c gpy; realk_c gpz;
} tstle4_grid_t;

// Convert Fortran-style 1-based (k,j,i) coordinates to a flat C index.
static inline size_t idx3(const intk_c kk, const intk_c jj, const intk_c k,
    const intk_c j, const intk_c i) {
    return (size_t)(k - 1) +
        (size_t)kk * ((size_t)(j - 1) + (size_t)jj * (size_t)(i - 1));
}

static inline size_t idx1(const intk_c i) { return (size_t)(i - 1); }

static inline realk_c max2(const realk_c a, const realk_c b) {
    return (a > b) ? a : b;
}

static inline realk_c abs_c(const realk_c x) {
#ifdef _MGLET_DOUBLE_PRECISION_
    return fabs(x);
#else
    return fabsf(x);
#endif
}

static inline realk_c pow_c(const realk_c x, const realk_c p) {
#ifdef _MGLET_DOUBLE_PRECISION_
    return pow(x, p);
#else
    return powf(x, p);
#endif
}

// Werner-Wengle wall model for one velocity component at one boundary cell.
static inline realk_c swcle3d_one_c(
    const realk_c ddz, const realk_c u, const realk_c gmol, const realk_c rho) {
    const realk_c cwa = (realk_c)8.3;
    const realk_c cwb = (realk_c)(1.0 / 7.0);
    const realk_c cpo1 = (realk_c)1.0 - cwb;
    const realk_c cpo2 = (realk_c)1.0 + cwb;
    const realk_c nu = gmol / rho;
    const realk_c cpo4 = cpo2 / cwa * pow_c(nu, cwb);
    const realk_c cpo5 =
        (realk_c)0.5 * cpo1 * pow_c(cwa, cpo2 / cpo1) * pow_c(nu, cpo2);
    const realk_c cpo6 = (realk_c)0.5 * nu * pow_c(cwa, (realk_c)2.0 / cpo1);
    const realk_c cpo8 = (realk_c)2.0 / cpo2;

    const realk_c vz = (u < (realk_c)0.0) ? (realk_c)-1.0 : (realk_c)1.0;
    const realk_c uquern = abs_c(u);

    if (uquern >= cpo6 / ddz) {
        const realk_c ddsb = pow_c(ddz, -cwb);
        return vz * pow_c(ddsb * (cpo4 * uquern + cpo5 / ddz), cpo8) / ddz;
    }

    return vz * (realk_c)2.0 * gmol * uquern / (rho * ddz * ddz);
}

#pragma omp begin declare target

// Convective contribution for u-momentum (skew-symmetric discretization).
void tstle4_kon_u_c(const intk_c kk, const intk_c jj, const intk_c ii,
    realk_c *uo, const realk_c *u, const realk_c *v, const realk_c *w,
    const realk_c *ut, const realk_c *vt, const realk_c *wt, const realk_c *dx,
    const realk_c *dy, const realk_c *dz, const realk_c *ddx,
    const realk_c *ddy, const realk_c *ddz, const realk_c *rdx,
    const realk_c *rdy, const realk_c *rdz, const realk_c *rddx,
    const realk_c *rddy, const realk_c *rddz, const intk_c nfro,
    const intk_c nbac) {
    intk_c nfu = 0;
    intk_c nbu = 0;

    (void)dy;
    (void)dz;
    (void)rdy;
    (void)rdz;
    (void)rddx;

    if (nbac == 7)
        nbu = 1;
    if (nfro == 3)
        nfu = 1;
    if (nbac == 3)
        nbu = 1;

    for (intk_c i = 3 - nfu; i <= ii - 3 + nbu; ++i) {
        for (intk_c j = 3; j <= jj - 2; ++j) {
            for (intk_c k = 3; k <= kk - 2; ++k) {
                const realk_c ax = ddy[idx1(j)] * ddz[idx1(k)];
                const realk_c ay = dx[idx1(i)] * ddz[idx1(k)];
                const realk_c az = dx[idx1(i)] * ddy[idx1(j)];

                const realk_c fe = ax *
                    (ut[idx3(kk, jj, k, j, i)] +
                        (ut[idx3(kk, jj, k, j, i + 1)] -
                            ut[idx3(kk, jj, k, j, i)]) *
                            (realk_c)0.5 * dx[idx1(i)] / ddx[idx1(i + 1)]);
                const realk_c fw = ax *
                    (ut[idx3(kk, jj, k, j, i - 1)] +
                        (ut[idx3(kk, jj, k, j, i)] -
                            ut[idx3(kk, jj, k, j, i - 1)]) *
                            (realk_c)0.5 * dx[idx1(i - 1)] / ddx[idx1(i)]);
                const realk_c fn = ay *
                    (vt[idx3(kk, jj, k, j, i)] +
                        vt[idx3(kk, jj, k, j, i + 1)]) *
                    (realk_c)0.5;
                const realk_c fs = ay *
                    (vt[idx3(kk, jj, k, j - 1, i)] +
                        vt[idx3(kk, jj, k, j - 1, i + 1)]) *
                    (realk_c)0.5;
                const realk_c ft = az *
                    (wt[idx3(kk, jj, k, j, i)] +
                        wt[idx3(kk, jj, k, j, i + 1)]) *
                    (realk_c)0.5;
                const realk_c fb = az *
                    (wt[idx3(kk, jj, k - 1, j, i)] +
                        wt[idx3(kk, jj, k - 1, j, i + 1)]) *
                    (realk_c)0.5;

                const realk_c qe = (realk_c)0.5 * fe *
                    (u[idx3(kk, jj, k, j, i)] + u[idx3(kk, jj, k, j, i + 1)]);
                const realk_c qw = (realk_c)0.5 * fw *
                    (u[idx3(kk, jj, k, j, i - 1)] + u[idx3(kk, jj, k, j, i)]);
                const realk_c qn = (realk_c)0.5 * fn *
                    (u[idx3(kk, jj, k, j, i)] + u[idx3(kk, jj, k, j + 1, i)]);
                const realk_c qs = (realk_c)0.5 * fs *
                    (u[idx3(kk, jj, k, j - 1, i)] + u[idx3(kk, jj, k, j, i)]);
                const realk_c qt = (realk_c)0.5 * ft *
                    (u[idx3(kk, jj, k, j, i)] + u[idx3(kk, jj, k + 1, j, i)]);
                const realk_c qb = (realk_c)0.5 * fb *
                    (u[idx3(kk, jj, k - 1, j, i)] + u[idx3(kk, jj, k, j, i)]);

                uo[idx3(kk, jj, k, j, i)] = -(qe - qw + qn - qs + qt - qb) *
                    rdx[idx1(i)] * rddy[idx1(j)] * rddz[idx1(k)];
            }
        }
    }
}

// Pressure-gradient source term with boundary-aware loop extents.
void tstle4_gradp_c(const intk_c kk, const intk_c jj, const intk_c ii,
    realk_c *uo, realk_c *vo, realk_c *wo, const realk_c *p, const realk_c *dx,
    const realk_c *dy, const realk_c *dz, const intk_c nfro, const intk_c nbac,
    const intk_c nrgt, const intk_c nlft, const intk_c nbot, const intk_c ntop,
    const realk_c gpx_in, const realk_c gpy_in, const realk_c gpz_in,
    const realk_c rho_in) {

    intk_c nfu = 0;
    intk_c nbu = 0;
    intk_c nrv = 0;
    intk_c nlv = 0;
    intk_c nbw = 0;
    intk_c ntw = 0;
    const realk_c inv_rho = (realk_c)1.0 / rho_in;

    if (nbac == 7)
        nbu = 1;
    if (nlft == 7)
        nlv = 1;
    if (ntop == 7)
        ntw = 1;

    if (nfro == 3)
        nfu = 1;
    if (nbac == 3)
        nbu = 1;
    if (nrgt == 3)
        nrv = 1;
    if (nlft == 3)
        nlv = 1;
    if (nbot == 3)
        nbw = 1;
    if (ntop == 3)
        ntw = 1;

#pragma omp for collapse(3)
    for (intk_c i = 3 - nfu; i <= ii - 3 + nbu; ++i) {
        for (intk_c j = 3; j <= jj - 2; ++j) {
            for (intk_c k = 3; k <= kk - 2; ++k) {
                const realk_c dxi = dx[idx1(i)];
                uo[idx3(kk, jj, k, j, i)] = uo[idx3(kk, jj, k, j, i)] -
                    inv_rho / dxi *
                        (p[idx3(kk, jj, k, j, i + 1)] -
                            p[idx3(kk, jj, k, j, i)] + gpx_in * dxi);
            }
        }
    }

// #pragma omp loop bind(parallel) collapse(3)
    for (intk_c i = 3; i <= ii - 2; ++i) {
        for (intk_c j = 3 - nrv; j <= jj - 3 + nlv; ++j) {
            for (intk_c k = 3; k <= kk - 2; ++k) {
                const realk_c dyj = dy[idx1(j)];
                vo[idx3(kk, jj, k, j, i)] = vo[idx3(kk, jj, k, j, i)] -
                    inv_rho / dyj *
                        (p[idx3(kk, jj, k, j + 1, i)] -
                            p[idx3(kk, jj, k, j, i)] + gpy_in * dyj);
            }
        }
    }

// #pragma omp loop bind(parallel) collapse(3)
    for (intk_c i = 3; i <= ii - 2; ++i) {
        for (intk_c j = 3; j <= jj - 2; ++j) {
            for (intk_c k = 3 - nbw; k <= kk - 3 + ntw; ++k) {
                const realk_c dzk = dz[idx1(k)];
                wo[idx3(kk, jj, k, j, i)] = wo[idx3(kk, jj, k, j, i)] -
                    inv_rho / dzk *
                        (p[idx3(kk, jj, k + 1, j, i)] -
                            p[idx3(kk, jj, k, j, i)] + gpz_in * dzk);
            }
        }
    }
}

// Full gradp implementation: iterate all grids in C using Fortran-provided metadata.
void tstle4_gradp_impl_c(const intk_c nmygrids_in, realk_c *uo, realk_c *vo,
    realk_c *wo, const realk_c *p, const realk_c *dx, const realk_c *dy,
    const realk_c *dz, const intk_c *kk_all, const intk_c *jj_all,
    const intk_c *ii_all, const intk_c *ip3_all, const intk_c *ipx_all,
    const intk_c *ipy_all, const intk_c *ipz_all, const intk_c *nfro_all,
    const intk_c *nbac_all, const intk_c *nrgt_all, const intk_c *nlft_all,
    const intk_c *nbot_all, const intk_c *ntop_all, const realk_c *gpx_all,
    const realk_c *gpy_all, const realk_c *gpz_all, const realk_c rho_in) {
    const realk_c inv_rho = (realk_c)1.0 / rho_in;

    #pragma omp target teams distribute
    for (intk_c igrid = 0; igrid < nmygrids_in; ++igrid) {

        const intk_c kk = kk_all[igrid];
        const intk_c jj = jj_all[igrid];
        const intk_c ii = ii_all[igrid];
        const intk_c nfro = nfro_all[igrid];
        const intk_c nbac = nbac_all[igrid];
        const intk_c nrgt = nrgt_all[igrid];
        const intk_c nlft = nlft_all[igrid];
        const intk_c nbot = nbot_all[igrid];
        const intk_c ntop = ntop_all[igrid];
        const realk_c gpx_in = gpx_all[igrid];
        const realk_c gpy_in = gpy_all[igrid];
        const realk_c gpz_in = gpz_all[igrid];

        realk_c *const uo_g = &uo[(size_t)(ip3_all[igrid] - 1)];
        realk_c *const vo_g = &vo[(size_t)(ip3_all[igrid] - 1)];
        realk_c *const wo_g = &wo[(size_t)(ip3_all[igrid] - 1)];
        const realk_c *const p_g = &p[(size_t)(ip3_all[igrid] - 1)];
        const realk_c *const dx_g = &dx[(size_t)(ipx_all[igrid] - 1)];
        const realk_c *const dy_g = &dy[(size_t)(ipy_all[igrid] - 1)];
        const realk_c *const dz_g = &dz[(size_t)(ipz_all[igrid] - 1)];

        intk_c nfu = 0;
        intk_c nbu = 0;
        intk_c nrv = 0;
        intk_c nlv = 0;
        intk_c nbw = 0;
        intk_c ntw = 0;

        if (nbac == 7)
            nbu = 1;
        if (nlft == 7)
            nlv = 1;
        if (ntop == 7)
            ntw = 1;

        if (nfro == 3)
            nfu = 1;
        if (nbac == 3)
            nbu = 1;
        if (nrgt == 3)
            nrv = 1;
        if (nlft == 3)
            nlv = 1;
        if (nbot == 3)
            nbw = 1;
        if (ntop == 3)
            ntw = 1;

        #pragma omp parallel
        {
            #pragma omp for collapse(3)
            for (intk_c i = 3 - nfu; i <= ii - 3 + nbu; ++i) {
                for (intk_c j = 3; j <= jj - 2; ++j) {
                    for (intk_c k = 3; k <= kk - 2; ++k) {
                        const realk_c dxi = dx_g[idx1(i)];
                        uo_g[idx3(kk, jj, k, j, i)] = uo_g[idx3(kk, jj, k, j, i)] -
                            inv_rho / dxi *
                                (p_g[idx3(kk, jj, k, j, i + 1)] -
                                    p_g[idx3(kk, jj, k, j, i)] + gpx_in * dxi);
                    }
                }
            }

            #pragma omp for collapse(3)
            for (intk_c i = 3; i <= ii - 2; ++i) {
                for (intk_c j = 3 - nrv; j <= jj - 3 + nlv; ++j) {
                    for (intk_c k = 3; k <= kk - 2; ++k) {
                        const realk_c dyj = dy_g[idx1(j)];
                        vo_g[idx3(kk, jj, k, j, i)] = vo_g[idx3(kk, jj, k, j, i)] -
                            inv_rho / dyj *
                                (p_g[idx3(kk, jj, k, j + 1, i)] -
                                    p_g[idx3(kk, jj, k, j, i)] + gpy_in * dyj);
                    }
                }
            }

            #pragma omp for collapse(3)
            for (intk_c i = 3; i <= ii - 2; ++i) {
                for (intk_c j = 3; j <= jj - 2; ++j) {
                    for (intk_c k = 3 - nbw; k <= kk - 3 + ntw; ++k) {
                        const realk_c dzk = dz_g[idx1(k)];
                        wo_g[idx3(kk, jj, k, j, i)] = wo_g[idx3(kk, jj, k, j, i)] -
                            inv_rho / dzk *
                                (p_g[idx3(kk, jj, k + 1, j, i)] -
                                    p_g[idx3(kk, jj, k, j, i)] + gpz_in * dzk);
                    }
                }
            }
        }
    }
}

// PAR boundary treatment and momentum-conserving flux redistribution.
void tstle4_par_c(const intk_c kk, const intk_c jj, const intk_c ii,
    realk_c *uo, realk_c *vo, realk_c *wo, const realk_c *u, const realk_c *v,
    const realk_c *w, const realk_c *ut, const realk_c *vt, const realk_c *wt,
    const realk_c *dx, const realk_c *dy, const realk_c *dz, const realk_c *ddx,
    const realk_c *ddy, const realk_c *ddz, const realk_c *rdx,
    const realk_c *rdy, const realk_c *rdz, const realk_c *rddx,
    const realk_c *rddy, const realk_c *rddz, realk_c *wcu, realk_c *wcv,
    realk_c *wcw, const intk_c nfro, const intk_c nbac, const intk_c nrgt,
    const intk_c nlft, const intk_c nbot, const intk_c ntop) {
    const realk_c wkon = (realk_c)1.0;
    intk_c i, j, k;

    if (nfro == 8) {
        i = 4;
#if defined(_OPENMP)
#pragma omp for
#endif
        for (intk_c j = 3; j <= jj - 2; ++j) {
            realk_c dyj = dy[idx1(j)];
            realk_c ddyj = ddy[idx1(j)];
            realk_c fkdtu = -1.0 * rddy[idx1(j)] * rdx[idx1(i)] * wkon;
            realk_c fkdtv = -1.0 * rddx[idx1(i)] * rdy[idx1(j)] * wkon;
            realk_c fkdtw = -1.0 * rddx[idx1(i)] * rddy[idx1(j)] * wkon;
            for (intk_c k = 3; k <= kk - 2; ++k) {
                realk_c dzk = dz[idx1(k)];
                realk_c rdzk = rdz[idx1(k)];
                realk_c ddzk = ddz[idx1(k)];
                realk_c rddzk = rddz[idx1(k)];
                realk_c avx = dyj * ddzk;
                realk_c awx = ddyj * dzk;
                realk_c fvw = avx *
                    (ut[idx3(kk, jj, k, j, i - 1)] +
                        ut[idx3(kk, jj, k, j + 1, i - 1)]) *
                    0.5;
                realk_c qkvwadd = 0.5 * (fvw - abs_c(fvw)) *
                    (0.5 * v[idx3(kk, jj, k, j, i)] -
                        1.5 * v[idx3(kk, jj, k, j, i - 1)] +
                        v[idx3(kk, jj, k, j, i - 2)]) *
                    0.5 * 0.5;
                vo[idx3(kk, jj, k, j, i)] =
                    vo[idx3(kk, jj, k, j, i)] + fkdtv * rddzk * (-qkvwadd);
                vo[idx3(kk, jj, k, j, i - 1)] =
                    vo[idx3(kk, jj, k, j, i - 1)] + fkdtv * rddzk * (+qkvwadd);
                realk_c fww = awx *
                    (ut[idx3(kk, jj, k, j, i - 1)] +
                        ut[idx3(kk, jj, k + 1, j, i - 1)]) *
                    0.5;
                realk_c qkwwadd = 0.5 * (fww - abs_c(fww)) *
                    (0.5 * w[idx3(kk, jj, k, j, i)] -
                        1.5 * w[idx3(kk, jj, k, j, i - 1)] +
                        w[idx3(kk, jj, k, j, i - 2)]) *
                    0.5 * 0.5;
                wo[idx3(kk, jj, k, j, i)] =
                    wo[idx3(kk, jj, k, j, i)] + fkdtw * rdzk * (-qkwwadd);
                wo[idx3(kk, jj, k, j, i - 1)] =
                    wo[idx3(kk, jj, k, j, i - 1)] + fkdtw * rdzk * (+qkwwadd);
            }
        }
#if defined(_OPENMP)
#pragma omp barrier
#endif
    }
    if (nbac == 8) {
        i = ii - 2;
#if defined(_OPENMP)
#pragma omp for
#endif
        for (intk_c j = 3; j <= jj - 2; ++j) {
            realk_c dyj = dy[idx1(j)];
            realk_c ddyj = ddy[idx1(j)];
            realk_c fkdtu = -1.0 * rddy[idx1(j)] * rdx[idx1(i)] * wkon;
            realk_c fkdtv = -1.0 * rddx[idx1(i)] * rdy[idx1(j)] * wkon;
            realk_c fkdtw = -1.0 * rddx[idx1(i)] * rddy[idx1(j)] * wkon;
            for (intk_c k = 3; k <= kk - 2; ++k) {
                realk_c dzk = dz[idx1(k)];
                realk_c rdzk = rdz[idx1(k)];
                realk_c ddzk = ddz[idx1(k)];
                realk_c rddzk = rddz[idx1(k)];
                realk_c avx = dyj * ddzk;
                realk_c awx = ddyj * dzk;
                realk_c fvw = avx *
                    (ut[idx3(kk, jj, k, j, i - 1)] +
                        ut[idx3(kk, jj, k, j + 1, i - 1)]) *
                    0.5;
                realk_c qkvwadd = 0.5 * (fvw + abs_c(fvw)) *
                    (0.5 * v[idx3(kk, jj, k, j, i - 1)] -
                        1.5 * v[idx3(kk, jj, k, j, i)] +
                        v[idx3(kk, jj, k, j, i + 1)]) *
                    0.5 * 0.5;
                vo[idx3(kk, jj, k, j, i)] =
                    vo[idx3(kk, jj, k, j, i)] + fkdtv * rddzk * (-qkvwadd);
                vo[idx3(kk, jj, k, j, i - 1)] =
                    vo[idx3(kk, jj, k, j, i - 1)] + fkdtv * rddzk * (+qkvwadd);
                realk_c fww = awx *
                    (ut[idx3(kk, jj, k, j, i - 1)] +
                        ut[idx3(kk, jj, k + 1, j, i - 1)]) *
                    0.5;
                realk_c qkwwadd = 0.5 * (fww + abs_c(fww)) *
                    (0.5 * w[idx3(kk, jj, k, j, i - 1)] -
                        1.5 * w[idx3(kk, jj, k, j, i)] +
                        w[idx3(kk, jj, k, j, i + 1)]) *
                    0.5 * 0.5;
                wo[idx3(kk, jj, k, j, i)] =
                    wo[idx3(kk, jj, k, j, i)] + fkdtw * rdzk * (-qkwwadd);
                wo[idx3(kk, jj, k, j, i - 1)] =
                    wo[idx3(kk, jj, k, j, i - 1)] + fkdtw * rdzk * (+qkwwadd);
            }
        }
#if defined(_OPENMP)
#pragma omp barrier
#endif
    }
    if (nrgt == 8) {
#if defined(_OPENMP)
#pragma omp for
#endif
        for (intk_c i = 3; i <= ii - 2; ++i) {
            j = 4;
            realk_c dxi = dx[idx1(i)];
            realk_c ddxi = ddx[idx1(i)];
            realk_c fkdtu = -1.0 * rddy[idx1(j)] * rdx[idx1(i)] * wkon;
            realk_c fkdtv = -1.0 * rddx[idx1(i)] * rdy[idx1(j)] * wkon;
            realk_c fkdtw = -1.0 * rddx[idx1(i)] * rddy[idx1(j)] * wkon;
            for (intk_c k = 3; k <= kk - 2; ++k) {
                realk_c dzk = dz[idx1(k)];
                realk_c rdzk = rdz[idx1(k)];
                realk_c ddzk = ddz[idx1(k)];
                realk_c rddzk = rddz[idx1(k)];
                realk_c auy = dxi * ddzk;
                realk_c awy = ddxi * dzk;
                realk_c fus = auy *
                    (vt[idx3(kk, jj, k, j - 1, i)] +
                        vt[idx3(kk, jj, k, j - 1, i + 1)]) *
                    0.5;
                realk_c qkusadd = 0.5 * (fus - abs_c(fus)) *
                    (0.5 * u[idx3(kk, jj, k, j, i)] -
                        1.5 * u[idx3(kk, jj, k, j - 1, i)] +
                        u[idx3(kk, jj, k, j - 2, i)]) *
                    0.5 * 0.5;
                uo[idx3(kk, jj, k, j, i)] =
                    uo[idx3(kk, jj, k, j, i)] + fkdtu * rddzk * (-qkusadd);
                uo[idx3(kk, jj, k, j - 1, i)] =
                    uo[idx3(kk, jj, k, j - 1, i)] + fkdtu * rddzk * (+qkusadd);
                realk_c fws = awy *
                    (vt[idx3(kk, jj, k, j - 1, i)] +
                        vt[idx3(kk, jj, k + 1, j - 1, i)]) *
                    0.5;
                realk_c qkwsadd = 0.5 * (fws - abs_c(fws)) *
                    (0.5 * w[idx3(kk, jj, k, j, i)] -
                        1.5 * w[idx3(kk, jj, k, j - 1, i)] +
                        w[idx3(kk, jj, k, j - 2, i)]) *
                    0.5 * 0.5;
                wo[idx3(kk, jj, k, j, i)] =
                    wo[idx3(kk, jj, k, j, i)] + fkdtw * rdzk * (-qkwsadd);
                wo[idx3(kk, jj, k, j - 1, i)] =
                    wo[idx3(kk, jj, k, j - 1, i)] + fkdtw * rdzk * (+qkwsadd);
            }
        }
#if defined(_OPENMP)
#pragma omp barrier
#endif
    }
    if (nlft == 8) {
#if defined(_OPENMP)
#pragma omp for
#endif
        for (intk_c i = 3; i <= ii - 2; ++i) {
            j = jj - 2;
            realk_c dxi = dx[idx1(i)];
            realk_c ddxi = ddx[idx1(i)];
            realk_c fkdtu = -1.0 * rddy[idx1(j)] * rdx[idx1(i)] * wkon;
            realk_c fkdtv = -1.0 * rddx[idx1(i)] * rdy[idx1(j)] * wkon;
            realk_c fkdtw = -1.0 * rddx[idx1(i)] * rddy[idx1(j)] * wkon;
            for (intk_c k = 3; k <= kk - 2; ++k) {
                realk_c dzk = dz[idx1(k)];
                realk_c rdzk = rdz[idx1(k)];
                realk_c ddzk = ddz[idx1(k)];
                realk_c rddzk = rddz[idx1(k)];
                realk_c auy = dxi * ddzk;
                realk_c awy = ddxi * dzk;
                realk_c fus = auy *
                    (vt[idx3(kk, jj, k, j - 1, i)] +
                        vt[idx3(kk, jj, k, j - 1, i + 1)]) *
                    0.5;
                realk_c qkusadd = 0.5 * (fus + abs_c(fus)) *
                    (0.5 * u[idx3(kk, jj, k, j - 1, i)] -
                        1.5 * u[idx3(kk, jj, k, j, i)] +
                        u[idx3(kk, jj, k, j + 1, i)]) *
                    0.5 * 0.5;
                uo[idx3(kk, jj, k, j, i)] =
                    uo[idx3(kk, jj, k, j, i)] + fkdtu * rddzk * (-qkusadd);
                uo[idx3(kk, jj, k, j - 1, i)] =
                    uo[idx3(kk, jj, k, j - 1, i)] + fkdtu * rddzk * (+qkusadd);
                realk_c fws = awy *
                    (vt[idx3(kk, jj, k, j - 1, i)] +
                        vt[idx3(kk, jj, k + 1, j - 1, i)]) *
                    0.5;
                realk_c qkwsadd = 0.5 * (fws + abs_c(fws)) *
                    (0.5 * w[idx3(kk, jj, k, j - 1, i)] -
                        1.5 * w[idx3(kk, jj, k, j, i)] +
                        w[idx3(kk, jj, k, j + 1, i)]) *
                    0.5 * 0.5;
                wo[idx3(kk, jj, k, j, i)] =
                    wo[idx3(kk, jj, k, j, i)] + fkdtw * rdzk * (-qkwsadd);
                wo[idx3(kk, jj, k, j - 1, i)] =
                    wo[idx3(kk, jj, k, j - 1, i)] + fkdtw * rdzk * (+qkwsadd);
            }
        }
#if defined(_OPENMP)
#pragma omp barrier
#endif
    }
    if (nbot == 8) {
#if defined(_OPENMP)
#pragma omp for
#endif
        for (intk_c i = 3; i <= ii - 2; ++i) {
            realk_c dxi = dx[idx1(i)];
            realk_c ddxi = ddx[idx1(i)];
            for (intk_c j = 3; j <= jj - 2; ++j) {
                k = 4;
                realk_c fkdtu = -1.0 * rddy[idx1(j)] * rdx[idx1(i)] * wkon;
                realk_c fkdtv = -1.0 * rddx[idx1(i)] * rdy[idx1(j)] * wkon;
                realk_c fkdtw = -1.0 * rddx[idx1(i)] * rddy[idx1(j)] * wkon;
                realk_c dyj = dy[idx1(j)];
                realk_c ddyj = ddy[idx1(j)];
                realk_c auz = dxi * ddyj;
                realk_c avz = ddxi * dyj;
                realk_c rdzk = rdz[idx1(k)];
                realk_c rddzk = rddz[idx1(k)];
                realk_c fub = auz *
                    (wt[idx3(kk, jj, k - 1, j, i)] +
                        wt[idx3(kk, jj, k - 1, j, i + 1)]) *
                    0.5;
                realk_c qkubadd = 0.5 * (fub - abs_c(fub)) *
                    (0.5 * u[idx3(kk, jj, k, j, i)] -
                        1.5 * u[idx3(kk, jj, k - 1, j, i)] +
                        u[idx3(kk, jj, k - 2, j, i)]) *
                    0.5 * 0.5;
                uo[idx3(kk, jj, k, j, i)] =
                    uo[idx3(kk, jj, k, j, i)] + fkdtu * rddzk * (-qkubadd);
                uo[idx3(kk, jj, k - 1, j, i)] =
                    uo[idx3(kk, jj, k - 1, j, i)] + fkdtu * rddzk * (+qkubadd);
                realk_c fvb = avz *
                    (wt[idx3(kk, jj, k - 1, j, i)] +
                        wt[idx3(kk, jj, k - 1, j + 1, i)]) *
                    0.5;
                realk_c qkvbadd = 0.5 * (fvb - abs_c(fvb)) *
                    (0.5 * v[idx3(kk, jj, k, j, i)] -
                        1.5 * v[idx3(kk, jj, k - 1, j, i)] +
                        v[idx3(kk, jj, k - 2, j, i)]) *
                    0.5 * 0.5;
                vo[idx3(kk, jj, k, j, i)] =
                    vo[idx3(kk, jj, k, j, i)] + fkdtv * rddzk * (-qkvbadd);
                vo[idx3(kk, jj, k - 1, j, i)] =
                    vo[idx3(kk, jj, k - 1, j, i)] + fkdtv * rddzk * (+qkvbadd);
            }
        }
#if defined(_OPENMP)
#pragma omp barrier
#endif
    }
    if (ntop == 8) {
#if defined(_OPENMP)
#pragma omp for
#endif
        for (intk_c i = 3; i <= ii - 2; ++i) {
            realk_c dxi = dx[idx1(i)];
            realk_c ddxi = ddx[idx1(i)];
            for (intk_c j = 3; j <= jj - 2; ++j) {
                k = kk - 2;
                realk_c fkdtu = -1.0 * rddy[idx1(j)] * rdx[idx1(i)] * wkon;
                realk_c fkdtv = -1.0 * rddx[idx1(i)] * rdy[idx1(j)] * wkon;
                realk_c fkdtw = -1.0 * rddx[idx1(i)] * rddy[idx1(j)] * wkon;
                realk_c dyj = dy[idx1(j)];
                realk_c ddyj = ddy[idx1(j)];
                realk_c auz = dxi * ddyj;
                realk_c avz = ddxi * dyj;
                realk_c rdzk = rdz[idx1(k)];
                realk_c rddzk = rddz[idx1(k)];
                realk_c fub = auz *
                    (wt[idx3(kk, jj, k - 1, j, i)] +
                        wt[idx3(kk, jj, k - 1, j, i + 1)]) *
                    0.5;
                realk_c qkubadd = 0.5 * (fub + abs_c(fub)) *
                    (0.5 * u[idx3(kk, jj, k - 1, j, i)] -
                        1.5 * u[idx3(kk, jj, k, j, i)] +
                        u[idx3(kk, jj, k + 1, j, i)]) *
                    0.5 * 0.5;
                uo[idx3(kk, jj, k, j, i)] =
                    uo[idx3(kk, jj, k, j, i)] + fkdtu * rddzk * (-qkubadd);
                uo[idx3(kk, jj, k - 1, j, i)] =
                    uo[idx3(kk, jj, k - 1, j, i)] + fkdtu * rddzk * (+qkubadd);
                realk_c fvb = avz *
                    (wt[idx3(kk, jj, k - 1, j, i)] +
                        wt[idx3(kk, jj, k - 1, j + 1, i)]) *
                    0.5;
                realk_c qkvbadd = 0.5 * (fvb + abs_c(fvb)) *
                    (0.5 * v[idx3(kk, jj, k - 1, j, i)] -
                        1.5 * v[idx3(kk, jj, k, j, i)] +
                        v[idx3(kk, jj, k + 1, j, i)]) *
                    0.5 * 0.5;
                vo[idx3(kk, jj, k, j, i)] =
                    vo[idx3(kk, jj, k, j, i)] + fkdtv * rddzk * (-qkvbadd);
                vo[idx3(kk, jj, k - 1, j, i)] =
                    vo[idx3(kk, jj, k - 1, j, i)] + fkdtv * rddzk * (+qkvbadd);
            }
        }
#if defined(_OPENMP)
#pragma omp barrier
#endif
    }
    if (nbac == 8) {
        i = ii - 2;
#if defined(_OPENMP)
#pragma omp for
#endif
        for (intk_c j = 3; j <= jj - 2; ++j) {
            realk_c ddyj = ddy[idx1(j)];
            for (intk_c k = 2; k <= kk - 2; ++k) {
                realk_c dzk = dz[idx1(k)];
                realk_c awx = ddyj * dzk;
                realk_c fwe = awx *
                    (ut[idx3(kk, jj, k, j, i)] +
                        ut[idx3(kk, jj, k + 1, j, i)]) *
                    0.5;
                realk_c qkwe = 0.5 * fwe *
                    (w[idx3(kk, jj, k, j, i)] + w[idx3(kk, jj, k, j, i + 1)]);
                wcw[idx3(kk, jj, k, j, i + 1)] = qkwe;
            }
        }
#if defined(_OPENMP)
#pragma omp barrier
#endif
        i = ii - 2;
#if defined(_OPENMP)
#pragma omp for
#endif
        for (intk_c j = 3; j <= jj - 2; j += 2) {
            realk_c ddyj = ddy[idx1(j)];
            k = 2;
            realk_c dzk = dz[idx1(k)];
            realk_c awx = ddyj * dzk;
            realk_c fwe = awx *
                (2 * ut[idx3(kk, jj, k, j, i)] +
                    2 * ut[idx3(kk, jj, k, j + 1, i)] +
                    ut[idx3(kk, jj, k + 1, j, i)] +
                    ut[idx3(kk, jj, k + 1, j + 1, i)] +
                    ut[idx3(kk, jj, k + 2, j, i)] +
                    ut[idx3(kk, jj, k + 2, j + 1, i)]) *
                0.125;
            realk_c qkwe = 0.5 * fwe *
                (w[idx3(kk, jj, k, j, i)] + w[idx3(kk, jj, k, j, i + 1)]);
            wcw[idx3(kk, jj, k, j, i)] = qkwe;
            qkwe = 0.5 * fwe *
                (w[idx3(kk, jj, k, j + 1, i)] +
                    w[idx3(kk, jj, k, j + 1, i + 1)]);
            wcw[idx3(kk, jj, k, j + 1, i)] = qkwe;
            for (intk_c k = 4; k <= kk - 4; k += 2) {
                dzk = dz[idx1(k)];
                awx = ddyj * dzk;
                fwe = awx *
                    (ut[idx3(kk, jj, k - 1, j, i)] +
                        ut[idx3(kk, jj, k - 1, j + 1, i)] +
                        ut[idx3(kk, jj, k, j, i)] +
                        ut[idx3(kk, jj, k, j + 1, i)] +
                        ut[idx3(kk, jj, k + 1, j, i)] +
                        ut[idx3(kk, jj, k + 1, j + 1, i)] +
                        ut[idx3(kk, jj, k + 2, j, i)] +
                        ut[idx3(kk, jj, k + 2, j + 1, i)]) *
                    0.125;
                qkwe = 0.5 * fwe *
                    (w[idx3(kk, jj, k, j, i)] + w[idx3(kk, jj, k, j, i + 1)]);
                wcw[idx3(kk, jj, k, j, i)] = qkwe;
                qkwe = 0.5 * fwe *
                    (w[idx3(kk, jj, k, j + 1, i)] +
                        w[idx3(kk, jj, k, j + 1, i + 1)]);
                wcw[idx3(kk, jj, k, j + 1, i)] = qkwe;
            }
            k = kk - 2;
            dzk = dz[idx1(k)];
            awx = ddyj * dzk;
            fwe = awx *
                (ut[idx3(kk, jj, k - 1, j, i)] +
                    ut[idx3(kk, jj, k - 1, j + 1, i)] +
                    ut[idx3(kk, jj, k, j, i)] + ut[idx3(kk, jj, k, j + 1, i)] +
                    2 * ut[idx3(kk, jj, k + 1, j, i)] +
                    2 * ut[idx3(kk, jj, k + 1, j + 1, i)]) *
                0.125;
            qkwe = 0.5 * fwe *
                (w[idx3(kk, jj, k, j, i)] + w[idx3(kk, jj, k, j, i + 1)]);
            wcw[idx3(kk, jj, k, j, i)] = qkwe;
            qkwe = 0.5 * fwe *
                (w[idx3(kk, jj, k, j + 1, i)] +
                    w[idx3(kk, jj, k, j + 1, i + 1)]);
            wcw[idx3(kk, jj, k, j + 1, i)] = qkwe;
        }
#if defined(_OPENMP)
#pragma omp barrier
#endif
        i = ii - 2;
#if defined(_OPENMP)
#pragma omp for collapse(2)
#endif
        for (intk_c j = 3; j <= jj - 2; ++j) {
            for (intk_c k = 2; k <= kk - 4; k += 2) {
                wcw[idx3(kk, jj, k + 1, j, i)] = 0.5 *
                    (wcw[idx3(kk, jj, k, j, i)] +
                        wcw[idx3(kk, jj, k + 2, j, i)]);
            }
        }
#if defined(_OPENMP)
#pragma omp barrier
#endif
#if defined(_OPENMP)
#pragma omp for
#endif
        for (intk_c j = 3; j <= jj - 2; ++j) {
            realk_c fkdtw = -1.0 * rddx[idx1(i)] * rddy[idx1(j)] * wkon;
            for (intk_c k = 2; k <= kk - 2; ++k) {
                realk_c rdzk = rdz[idx1(k)];
                wo[idx3(kk, jj, k, j, i)] = wo[idx3(kk, jj, k, j, i)] +
                    fkdtw * rdzk *
                        (-wcw[idx3(kk, jj, k, j, i + 1)] +
                            wcw[idx3(kk, jj, k, j, i)]);
            }
        }
#if defined(_OPENMP)
#pragma omp barrier
#endif
        i = ii - 2;
#if defined(_OPENMP)
#pragma omp for
#endif
        for (intk_c j = 2; j <= jj - 2; ++j) {
            realk_c dyj = dy[idx1(j)];
            for (intk_c k = 3; k <= kk - 2; ++k) {
                realk_c ddzk = ddz[idx1(k)];
                realk_c avx = ddzk * dyj;
                realk_c fve = avx *
                    (ut[idx3(kk, jj, k, j, i)] +
                        ut[idx3(kk, jj, k, j + 1, i)]) *
                    0.5;
                realk_c qkve = 0.5 * fve *
                    (v[idx3(kk, jj, k, j, i)] + v[idx3(kk, jj, k, j, i + 1)]);
                wcv[idx3(kk, jj, k, j, i + 1)] = qkve;
            }
        }
#if defined(_OPENMP)
#pragma omp barrier
#endif
        i = ii - 2;
        j = 2;
        realk_c dyj = dy[idx1(j)];
#if defined(_OPENMP)
#pragma omp for
#endif
        for (intk_c k = 3; k <= kk - 2; k += 2) {
            realk_c ddzk = ddz[idx1(k)];
            realk_c avx = ddzk * dyj;
            realk_c fve = avx *
                (2 * ut[idx3(kk, jj, k, j, i)] +
                    2 * ut[idx3(kk, jj, k + 1, j, i)] +
                    ut[idx3(kk, jj, k, j + 1, i)] +
                    ut[idx3(kk, jj, k + 1, j + 1, i)] +
                    ut[idx3(kk, jj, k, j + 2, i)] +
                    ut[idx3(kk, jj, k + 1, j + 2, i)]) *
                0.125;
            realk_c qkve = 0.5 * fve *
                (v[idx3(kk, jj, k, j, i)] + v[idx3(kk, jj, k, j, i + 1)]);
            wcv[idx3(kk, jj, k, j, i)] = qkve;
            qkve = 0.5 * fve *
                (v[idx3(kk, jj, k + 1, j, i)] +
                    v[idx3(kk, jj, k + 1, j, i + 1)]);
            wcv[idx3(kk, jj, k + 1, j, i)] = qkve;
        }
#if defined(_OPENMP)
#pragma omp barrier
#endif
#if defined(_OPENMP)
#pragma omp for
#endif
        for (intk_c j = 4; j <= jj - 4; j += 2) {
            dyj = dy[idx1(j)];
            for (intk_c k = 3; k <= kk - 2; k += 2) {
                realk_c ddzk = ddz[idx1(k)];
                realk_c avx = ddzk * dyj;
                realk_c fve = avx *
                    (ut[idx3(kk, jj, k, j - 1, i)] +
                        ut[idx3(kk, jj, k + 1, j - 1, i)] +
                        ut[idx3(kk, jj, k, j, i)] +
                        ut[idx3(kk, jj, k + 1, j, i)] +
                        ut[idx3(kk, jj, k, j + 1, i)] +
                        ut[idx3(kk, jj, k + 1, j + 1, i)] +
                        ut[idx3(kk, jj, k, j + 2, i)] +
                        ut[idx3(kk, jj, k + 1, j + 2, i)]) *
                    0.125;
                realk_c qkve = 0.5 * fve *
                    (v[idx3(kk, jj, k, j, i)] + v[idx3(kk, jj, k, j, i + 1)]);
                wcv[idx3(kk, jj, k, j, i)] = qkve;
                qkve = 0.5 * fve *
                    (v[idx3(kk, jj, k + 1, j, i)] +
                        v[idx3(kk, jj, k + 1, j, i + 1)]);
                wcv[idx3(kk, jj, k + 1, j, i)] = qkve;
            }
        }
#if defined(_OPENMP)
#pragma omp barrier
#endif
        j = jj - 2;
        dyj = dy[idx1(j)];
#if defined(_OPENMP)
#pragma omp for
#endif
        for (intk_c k = 3; k <= kk - 2; k += 2) {
            realk_c ddzk = ddz[idx1(k)];
            realk_c avx = ddzk * dyj;
            realk_c fve = avx *
                (ut[idx3(kk, jj, k, j - 1, i)] +
                    ut[idx3(kk, jj, k + 1, j - 1, i)] +
                    ut[idx3(kk, jj, k, j, i)] + ut[idx3(kk, jj, k + 1, j, i)] +
                    2 * ut[idx3(kk, jj, k, j + 1, i)] +
                    2 * ut[idx3(kk, jj, k + 1, j + 1, i)]) *
                0.125;
            realk_c qkve = 0.5 * fve *
                (v[idx3(kk, jj, k, j, i)] + v[idx3(kk, jj, k, j, i + 1)]);
            wcv[idx3(kk, jj, k, j, i)] = qkve;
            qkve = 0.5 * fve *
                (v[idx3(kk, jj, k + 1, j, i)] +
                    v[idx3(kk, jj, k + 1, j, i + 1)]);
            wcv[idx3(kk, jj, k + 1, j, i)] = qkve;
        }
#if defined(_OPENMP)
#pragma omp barrier
#endif
        i = ii - 2;
#if defined(_OPENMP)
#pragma omp for collapse(2)
#endif
        for (intk_c j = 2; j <= jj - 4; j += 2) {
            for (intk_c k = 3; k <= kk - 2; ++k) {
                wcv[idx3(kk, jj, k, j + 1, i)] = 0.5 *
                    (wcv[idx3(kk, jj, k, j, i)] +
                        wcv[idx3(kk, jj, k, j + 2, i)]);
            }
        }
#if defined(_OPENMP)
#pragma omp barrier
#endif
#if defined(_OPENMP)
#pragma omp for
#endif
        for (intk_c j = 2; j <= jj - 2; ++j) {
            realk_c fkdtv = -1.0 * rddx[idx1(i)] * rdy[idx1(j)] * wkon;
            for (intk_c k = 3; k <= kk - 2; ++k) {
                realk_c rddzk = rddz[idx1(k)];
                vo[idx3(kk, jj, k, j, i)] = vo[idx3(kk, jj, k, j, i)] +
                    fkdtv * rddzk *
                        (wcv[idx3(kk, jj, k, j, i)] -
                            wcv[idx3(kk, jj, k, j, i + 1)]);
            }
        }
#if defined(_OPENMP)
#pragma omp barrier
#endif
    }
    if (nfro == 8) {
        i = 3;
#if defined(_OPENMP)
#pragma omp for
#endif
        for (intk_c j = 3; j <= jj - 2; ++j) {
            realk_c ddyj = ddy[idx1(j)];
            for (intk_c k = 2; k <= kk - 2; ++k) {
                realk_c dzk = dz[idx1(k)];
                realk_c awx = ddyj * dzk;
                realk_c fww = awx *
                    (ut[idx3(kk, jj, k, j, i - 1)] +
                        ut[idx3(kk, jj, k + 1, j, i - 1)]) *
                    0.5;
                realk_c qkww = 0.5 * fww *
                    (w[idx3(kk, jj, k, j, i)] + w[idx3(kk, jj, k, j, i - 1)]);
                wcw[idx3(kk, jj, k, j, i - 1)] = qkww;
            }
        }
#if defined(_OPENMP)
#pragma omp barrier
#endif
        i = 3;
#if defined(_OPENMP)
#pragma omp for
#endif
        for (intk_c j = 3; j <= jj - 2; j += 2) {
            realk_c ddyj = ddy[idx1(j)];
            k = 2;
            realk_c dzk = dz[idx1(k)];
            realk_c awx = ddyj * dzk;
            realk_c fww = awx *
                (2 * ut[idx3(kk, jj, k, j, i - 1)] +
                    2 * ut[idx3(kk, jj, k, j + 1, i - 1)] +
                    ut[idx3(kk, jj, k + 1, j, i - 1)] +
                    ut[idx3(kk, jj, k + 1, j + 1, i - 1)] +
                    ut[idx3(kk, jj, k + 2, j, i - 1)] +
                    ut[idx3(kk, jj, k + 2, j + 1, i - 1)]) *
                0.125;
            realk_c qkww = 0.5 * fww *
                (w[idx3(kk, jj, k, j, i)] + w[idx3(kk, jj, k, j, i - 1)]);
            wcw[idx3(kk, jj, k, j, i)] = qkww;
            qkww = 0.5 * fww *
                (w[idx3(kk, jj, k, j + 1, i)] +
                    w[idx3(kk, jj, k, j + 1, i - 1)]);
            wcw[idx3(kk, jj, k, j + 1, i)] = qkww;
            for (intk_c k = 4; k <= kk - 4; k += 2) {
                dzk = dz[idx1(k)];
                awx = ddyj * dzk;
                fww = awx *
                    (ut[idx3(kk, jj, k - 1, j, i - 1)] +
                        ut[idx3(kk, jj, k - 1, j + 1, i - 1)] +
                        ut[idx3(kk, jj, k, j, i - 1)] +
                        ut[idx3(kk, jj, k, j + 1, i - 1)] +
                        ut[idx3(kk, jj, k + 1, j, i - 1)] +
                        ut[idx3(kk, jj, k + 1, j + 1, i - 1)] +
                        ut[idx3(kk, jj, k + 2, j, i - 1)] +
                        ut[idx3(kk, jj, k + 2, j + 1, i - 1)]) *
                    0.125;
                qkww = 0.5 * fww *
                    (w[idx3(kk, jj, k, j, i)] + w[idx3(kk, jj, k, j, i - 1)]);
                wcw[idx3(kk, jj, k, j, i)] = qkww;
                qkww = 0.5 * fww *
                    (w[idx3(kk, jj, k, j + 1, i)] +
                        w[idx3(kk, jj, k, j + 1, i - 1)]);
                wcw[idx3(kk, jj, k, j + 1, i)] = qkww;
            }
            k = kk - 2;
            dzk = dz[idx1(k)];
            awx = ddyj * dzk;
            fww = awx *
                (ut[idx3(kk, jj, k - 1, j, i - 1)] +
                    ut[idx3(kk, jj, k - 1, j + 1, i - 1)] +
                    ut[idx3(kk, jj, k, j, i - 1)] +
                    ut[idx3(kk, jj, k, j + 1, i - 1)] +
                    2 * ut[idx3(kk, jj, k + 1, j, i - 1)] +
                    2 * ut[idx3(kk, jj, k + 1, j + 1, i - 1)]) *
                0.125;
            qkww = 0.5 * fww *
                (w[idx3(kk, jj, k, j, i)] + w[idx3(kk, jj, k, j, i - 1)]);
            wcw[idx3(kk, jj, k, j, i)] = qkww;
            qkww = 0.5 * fww *
                (w[idx3(kk, jj, k, j + 1, i)] +
                    w[idx3(kk, jj, k, j + 1, i - 1)]);
            wcw[idx3(kk, jj, k, j + 1, i)] = qkww;
        }
#if defined(_OPENMP)
#pragma omp barrier
#endif
        i = 3;
#if defined(_OPENMP)
#pragma omp for collapse(2)
#endif
        for (intk_c j = 3; j <= jj - 2; ++j) {
            for (intk_c k = 2; k <= kk - 4; k += 2) {
                wcw[idx3(kk, jj, k + 1, j, i)] = 0.5 *
                    (wcw[idx3(kk, jj, k, j, i)] +
                        wcw[idx3(kk, jj, k + 2, j, i)]);
            }
        }
#if defined(_OPENMP)
#pragma omp barrier
#endif
#if defined(_OPENMP)
#pragma omp for
#endif
        for (intk_c j = 3; j <= jj - 2; ++j) {
            realk_c fkdtw = -1.0 * rddx[idx1(i)] * rddy[idx1(j)] * wkon;
            for (intk_c k = 2; k <= kk - 2; ++k) {
                realk_c rdzk = rdz[idx1(k)];
                wo[idx3(kk, jj, k, j, i)] = wo[idx3(kk, jj, k, j, i)] +
                    fkdtw * rdzk *
                        (wcw[idx3(kk, jj, k, j, i - 1)] -
                            wcw[idx3(kk, jj, k, j, i)]);
            }
        }
#if defined(_OPENMP)
#pragma omp barrier
#endif
        i = 3;
#if defined(_OPENMP)
#pragma omp for
#endif
        for (intk_c j = 2; j <= jj - 2; ++j) {
            realk_c dyj = dy[idx1(j)];
            for (intk_c k = 3; k <= kk - 2; ++k) {
                realk_c ddzk = ddz[idx1(k)];
                realk_c avx = ddzk * dyj;
                realk_c fvw = avx *
                    (ut[idx3(kk, jj, k, j, i - 1)] +
                        ut[idx3(kk, jj, k, j + 1, i - 1)]) *
                    0.5;
                realk_c qkvw = 0.5 * fvw *
                    (v[idx3(kk, jj, k, j, i)] + v[idx3(kk, jj, k, j, i - 1)]);
                wcv[idx3(kk, jj, k, j, i - 1)] = qkvw;
            }
        }
#if defined(_OPENMP)
#pragma omp barrier
#endif
        i = 3;
        j = 2;
        realk_c dyj = dy[idx1(j)];
#if defined(_OPENMP)
#pragma omp for
#endif
        for (intk_c k = 3; k <= kk - 2; k += 2) {
            realk_c ddzk = ddz[idx1(k)];
            realk_c avx = ddzk * dyj;
            realk_c fvw = avx *
                (2 * ut[idx3(kk, jj, k, j, i - 1)] +
                    2 * ut[idx3(kk, jj, k + 1, j, i - 1)] +
                    ut[idx3(kk, jj, k, j + 1, i - 1)] +
                    ut[idx3(kk, jj, k + 1, j + 1, i - 1)] +
                    ut[idx3(kk, jj, k, j + 2, i - 1)] +
                    ut[idx3(kk, jj, k + 1, j + 2, i - 1)]) *
                0.125;
            realk_c qkvw = 0.5 * fvw *
                (v[idx3(kk, jj, k, j, i)] + v[idx3(kk, jj, k, j, i - 1)]);
            wcv[idx3(kk, jj, k, j, i)] = qkvw;
            qkvw = 0.5 * fvw *
                (v[idx3(kk, jj, k + 1, j, i)] +
                    v[idx3(kk, jj, k + 1, j, i - 1)]);
            wcv[idx3(kk, jj, k + 1, j, i)] = qkvw;
        }
#if defined(_OPENMP)
#pragma omp barrier
#endif
#if defined(_OPENMP)
#pragma omp for
#endif
        for (intk_c j = 4; j <= jj - 4; j += 2) {
            dyj = dy[idx1(j)];
            for (intk_c k = 3; k <= kk - 2; k += 2) {
                realk_c ddzk = ddz[idx1(k)];
                realk_c avx = ddzk * dyj;
                realk_c fvw = avx *
                    (ut[idx3(kk, jj, k, j - 1, i - 1)] +
                        ut[idx3(kk, jj, k + 1, j - 1, i - 1)] +
                        ut[idx3(kk, jj, k, j, i - 1)] +
                        ut[idx3(kk, jj, k + 1, j, i - 1)] +
                        ut[idx3(kk, jj, k, j + 1, i - 1)] +
                        ut[idx3(kk, jj, k + 1, j + 1, i - 1)] +
                        ut[idx3(kk, jj, k, j + 2, i - 1)] +
                        ut[idx3(kk, jj, k + 1, j + 2, i - 1)]) *
                    0.125;
                realk_c qkvw = 0.5 * fvw *
                    (v[idx3(kk, jj, k, j, i)] + v[idx3(kk, jj, k, j, i - 1)]);
                wcv[idx3(kk, jj, k, j, i)] = qkvw;
                qkvw = 0.5 * fvw *
                    (v[idx3(kk, jj, k + 1, j, i)] +
                        v[idx3(kk, jj, k + 1, j, i - 1)]);
                wcv[idx3(kk, jj, k + 1, j, i)] = qkvw;
            }
        }
#if defined(_OPENMP)
#pragma omp barrier
#endif
        j = jj - 2;
        dyj = dy[idx1(j)];
#if defined(_OPENMP)
#pragma omp for
#endif
        for (intk_c k = 3; k <= kk - 2; k += 2) {
            realk_c ddzk = ddz[idx1(k)];
            realk_c avx = ddzk * dyj;
            realk_c fvw = avx *
                (ut[idx3(kk, jj, k, j - 1, i - 1)] +
                    ut[idx3(kk, jj, k + 1, j - 1, i - 1)] +
                    ut[idx3(kk, jj, k, j, i - 1)] +
                    ut[idx3(kk, jj, k + 1, j, i - 1)] +
                    2 * ut[idx3(kk, jj, k, j + 1, i - 1)] +
                    2 * ut[idx3(kk, jj, k + 1, j + 1, i - 1)]) *
                0.125;
            realk_c qkvw = 0.5 * fvw *
                (v[idx3(kk, jj, k, j, i)] + v[idx3(kk, jj, k, j, i - 1)]);
            wcv[idx3(kk, jj, k, j, i)] = qkvw;
            qkvw = 0.5 * fvw *
                (v[idx3(kk, jj, k + 1, j, i)] +
                    v[idx3(kk, jj, k + 1, j, i - 1)]);
            wcv[idx3(kk, jj, k + 1, j, i)] = qkvw;
        }
#if defined(_OPENMP)
#pragma omp barrier
#endif
        i = 3;
#if defined(_OPENMP)
#pragma omp for collapse(2)
#endif
        for (intk_c j = 2; j <= jj - 4; j += 2) {
            for (intk_c k = 3; k <= kk - 2; ++k) {
                wcv[idx3(kk, jj, k, j + 1, i)] = 0.5 *
                    (wcv[idx3(kk, jj, k, j, i)] +
                        wcv[idx3(kk, jj, k, j + 2, i)]);
            }
        }
#if defined(_OPENMP)
#pragma omp barrier
#endif
#if defined(_OPENMP)
#pragma omp for
#endif
        for (intk_c j = 2; j <= jj - 2; ++j) {
            realk_c fkdtv = -1.0 * rddx[idx1(i)] * rdy[idx1(j)] * wkon;
            for (intk_c k = 3; k <= kk - 2; ++k) {
                realk_c rddzk = rddz[idx1(k)];
                vo[idx3(kk, jj, k, j, i)] = vo[idx3(kk, jj, k, j, i)] +
                    fkdtv * rddzk *
                        (-wcv[idx3(kk, jj, k, j, i)] +
                            wcv[idx3(kk, jj, k, j, i - 1)]);
            }
        }
#if defined(_OPENMP)
#pragma omp barrier
#endif
    }
    if (ntop == 8) {
        k = kk - 2;
#if defined(_OPENMP)
#pragma omp for
#endif
        for (intk_c i = 2; i <= ii - 2; ++i) {
            realk_c dxi = dx[idx1(i)];
            for (intk_c j = 3; j <= jj - 2; ++j) {
                realk_c ddyj = ddy[idx1(j)];
                realk_c auz = dxi * ddyj;
                realk_c fut = auz *
                    (wt[idx3(kk, jj, k, j, i)] +
                        wt[idx3(kk, jj, k, j, i + 1)]) *
                    0.5;
                realk_c qkut = 0.5 * fut *
                    (u[idx3(kk, jj, k, j, i)] + u[idx3(kk, jj, k + 1, j, i)]);
                wcu[idx3(kk, jj, k + 1, j, i)] = qkut;
            }
        }
#if defined(_OPENMP)
#pragma omp barrier
#endif
        k = kk - 2;
        i = 2;
        realk_c dxi = dx[idx1(i)];
#if defined(_OPENMP)
#pragma omp for
#endif
        for (intk_c j = 3; j <= jj - 2; j += 2) {
            realk_c ddyj = ddy[idx1(j)];
            realk_c auz = dxi * ddyj;
            realk_c fut = auz *
                (2 * wt[idx3(kk, jj, k, j, i)] +
                    2 * wt[idx3(kk, jj, k, j + 1, i)] +
                    wt[idx3(kk, jj, k, j, i + 1)] +
                    wt[idx3(kk, jj, k, j + 1, i + 1)] +
                    wt[idx3(kk, jj, k, j, i + 2)] +
                    wt[idx3(kk, jj, k, j + 1, i + 2)]) *
                0.125;
            realk_c qkut = 0.5 * fut *
                (u[idx3(kk, jj, k, j, i)] + u[idx3(kk, jj, k + 1, j, i)]);
            wcu[idx3(kk, jj, k, j, i)] = qkut;
            qkut = 0.5 * fut *
                (u[idx3(kk, jj, k, j + 1, i)] +
                    u[idx3(kk, jj, k + 1, j + 1, i)]);
            wcu[idx3(kk, jj, k, j + 1, i)] = qkut;
        }
#if defined(_OPENMP)
#pragma omp barrier
#endif
#if defined(_OPENMP)
#pragma omp for
#endif
        for (intk_c i = 4; i <= ii - 4; i += 2) {
            dxi = dx[idx1(i)];
            for (intk_c j = 3; j <= jj - 2; j += 2) {
                realk_c ddyj = ddy[idx1(j)];
                realk_c auz = dxi * ddyj;
                realk_c fut = auz *
                    (wt[idx3(kk, jj, k, j, i - 1)] +
                        wt[idx3(kk, jj, k, j + 1, i - 1)] +
                        wt[idx3(kk, jj, k, j, i)] +
                        wt[idx3(kk, jj, k, j + 1, i)] +
                        wt[idx3(kk, jj, k, j, i + 1)] +
                        wt[idx3(kk, jj, k, j + 1, i + 1)] +
                        wt[idx3(kk, jj, k, j, i + 2)] +
                        wt[idx3(kk, jj, k, j + 1, i + 2)]) *
                    0.125;
                realk_c qkut = 0.5 * fut *
                    (u[idx3(kk, jj, k, j, i)] + u[idx3(kk, jj, k + 1, j, i)]);
                wcu[idx3(kk, jj, k, j, i)] = qkut;
                qkut = 0.5 * fut *
                    (u[idx3(kk, jj, k, j + 1, i)] +
                        u[idx3(kk, jj, k + 1, j + 1, i)]);
                wcu[idx3(kk, jj, k, j + 1, i)] = qkut;
            }
        }
#if defined(_OPENMP)
#pragma omp barrier
#endif
        i = ii - 2;
        dxi = dx[idx1(i)];
#if defined(_OPENMP)
#pragma omp for
#endif
        for (intk_c j = 3; j <= jj - 2; j += 2) {
            realk_c ddyj = ddy[idx1(j)];
            realk_c auz = dxi * ddyj;
            realk_c fut = auz *
                (wt[idx3(kk, jj, k, j, i - 1)] +
                    wt[idx3(kk, jj, k, j + 1, i - 1)] +
                    wt[idx3(kk, jj, k, j, i)] + wt[idx3(kk, jj, k, j + 1, i)] +
                    2 * wt[idx3(kk, jj, k, j, i + 1)] +
                    2 * wt[idx3(kk, jj, k, j + 1, i + 1)]) *
                0.125;
            realk_c qkut = 0.5 * fut *
                (u[idx3(kk, jj, k, j, i)] + u[idx3(kk, jj, k + 1, j, i)]);
            wcu[idx3(kk, jj, k, j, i)] = qkut;
            qkut = 0.5 * fut *
                (u[idx3(kk, jj, k, j + 1, i)] +
                    u[idx3(kk, jj, k + 1, j + 1, i)]);
            wcu[idx3(kk, jj, k, j + 1, i)] = qkut;
        }
#if defined(_OPENMP)
#pragma omp barrier
#endif
        k = kk - 2;
#if defined(_OPENMP)
#pragma omp for collapse(1)
#endif
        for (intk_c i = 2; i <= ii - 4; i += 2) {
            for (intk_c j = 3; j <= jj - 2; ++j) {
                wcu[idx3(kk, jj, k, j, i + 1)] = 0.5 *
                    (wcu[idx3(kk, jj, k, j, i)] +
                        wcu[idx3(kk, jj, k, j, i + 2)]);
            }
        }
#if defined(_OPENMP)
#pragma omp barrier
#endif
        realk_c rddzk = rddz[idx1(k)];
#if defined(_OPENMP)
#pragma omp for collapse(2)
#endif
        for (intk_c i = 2; i <= ii - 2; ++i) {
            for (intk_c j = 3; j <= jj - 2; ++j) {
                realk_c fkdtu = -1.0 * rddy[idx1(j)] * rdx[idx1(i)] * wkon;
                uo[idx3(kk, jj, k, j, i)] = uo[idx3(kk, jj, k, j, i)] +
                    fkdtu * rddzk *
                        (wcu[idx3(kk, jj, k, j, i)] -
                            wcu[idx3(kk, jj, k + 1, j, i)]);
            }
        }
#if defined(_OPENMP)
#pragma omp barrier
#endif
        k = kk - 2;
#if defined(_OPENMP)
#pragma omp for
#endif
        for (intk_c i = 3; i <= ii - 2; ++i) {
            realk_c ddxi = ddx[idx1(i)];
            for (intk_c j = 2; j <= jj - 2; ++j) {
                realk_c dyj = dy[idx1(j)];
                realk_c avz = ddxi * dyj;
                realk_c fvt = avz *
                    (wt[idx3(kk, jj, k, j, i)] +
                        wt[idx3(kk, jj, k, j + 1, i)]) *
                    0.5;
                realk_c qkvt = 0.5 * fvt *
                    (v[idx3(kk, jj, k, j, i)] + v[idx3(kk, jj, k + 1, j, i)]);
                wcv[idx3(kk, jj, k + 1, j, i)] = qkvt;
            }
        }
#if defined(_OPENMP)
#pragma omp barrier
#endif
        k = kk - 2;
#if defined(_OPENMP)
#pragma omp for
#endif
        for (intk_c i = 3; i <= ii - 2; i += 2) {
            realk_c ddxi = ddx[idx1(i)];
            j = 2;
            realk_c dyj = dy[idx1(j)];
            realk_c avz = ddxi * dyj;
            realk_c fvt = avz *
                (2 * wt[idx3(kk, jj, k, j, i)] +
                    2 * wt[idx3(kk, jj, k, j, i + 1)] +
                    wt[idx3(kk, jj, k, j + 1, i)] +
                    wt[idx3(kk, jj, k, j + 1, i + 1)] +
                    wt[idx3(kk, jj, k, j + 2, i)] +
                    wt[idx3(kk, jj, k, j + 2, i + 1)]) *
                0.125;
            realk_c qkvt = 0.5 * fvt *
                (v[idx3(kk, jj, k, j, i)] + v[idx3(kk, jj, k + 1, j, i)]);
            wcv[idx3(kk, jj, k, j, i)] = qkvt;
            qkvt = 0.5 * fvt *
                (v[idx3(kk, jj, k, j, i + 1)] +
                    v[idx3(kk, jj, k + 1, j, i + 1)]);
            wcv[idx3(kk, jj, k, j, i + 1)] = qkvt;
            for (intk_c j = 4; j <= jj - 4; j += 2) {
                dyj = dy[idx1(j)];
                avz = ddxi * dyj;
                fvt = avz *
                    (wt[idx3(kk, jj, k, j - 1, i)] +
                        wt[idx3(kk, jj, k, j - 1, i + 1)] +
                        wt[idx3(kk, jj, k, j, i)] +
                        wt[idx3(kk, jj, k, j, i + 1)] +
                        wt[idx3(kk, jj, k, j + 1, i)] +
                        wt[idx3(kk, jj, k, j + 1, i + 1)] +
                        wt[idx3(kk, jj, k, j + 2, i)] +
                        wt[idx3(kk, jj, k, j + 2, i + 1)]) *
                    0.125;
                qkvt = 0.5 * fvt *
                    (v[idx3(kk, jj, k, j, i)] + v[idx3(kk, jj, k + 1, j, i)]);
                wcv[idx3(kk, jj, k, j, i)] = qkvt;
                qkvt = 0.5 * fvt *
                    (v[idx3(kk, jj, k, j, i + 1)] +
                        v[idx3(kk, jj, k + 1, j, i + 1)]);
                wcv[idx3(kk, jj, k, j, i + 1)] = qkvt;
            }
            j = jj - 2;
            dyj = dy[idx1(j)];
            avz = ddxi * dyj;
            fvt = avz *
                (wt[idx3(kk, jj, k, j - 1, i)] +
                    wt[idx3(kk, jj, k, j - 1, i + 1)] +
                    wt[idx3(kk, jj, k, j, i)] + wt[idx3(kk, jj, k, j, i + 1)] +
                    2 * wt[idx3(kk, jj, k, j + 1, i)] +
                    2 * wt[idx3(kk, jj, k, j + 1, i + 1)]) *
                0.125;
            qkvt = 0.5 * fvt *
                (v[idx3(kk, jj, k, j, i)] + v[idx3(kk, jj, k + 1, j, i)]);
            wcv[idx3(kk, jj, k, j, i)] = qkvt;
            qkvt = 0.5 * fvt *
                (v[idx3(kk, jj, k, j, i + 1)] +
                    v[idx3(kk, jj, k + 1, j, i + 1)]);
            wcv[idx3(kk, jj, k, j, i + 1)] = qkvt;
        }
#if defined(_OPENMP)
#pragma omp barrier
#endif
        k = kk - 2;
#if defined(_OPENMP)
#pragma omp for collapse(2)
#endif
        for (intk_c i = 3; i <= ii - 2; ++i) {
            for (intk_c j = 2; j <= jj - 4; j += 2) {
                wcv[idx3(kk, jj, k, j + 1, i)] = 0.5 *
                    (wcv[idx3(kk, jj, k, j, i)] +
                        wcv[idx3(kk, jj, k, j + 2, i)]);
            }
        }
#if defined(_OPENMP)
#pragma omp barrier
#endif
        rddzk = rddz[idx1(k)];
#if defined(_OPENMP)
#pragma omp for collapse(2)
#endif
        for (intk_c i = 3; i <= ii - 2; ++i) {
            for (intk_c j = 2; j <= jj - 2; ++j) {
                realk_c fkdtv = -1.0 * rddx[idx1(i)] * rdy[idx1(j)] * wkon;
                vo[idx3(kk, jj, k, j, i)] = vo[idx3(kk, jj, k, j, i)] +
                    fkdtv * rddzk *
                        (wcv[idx3(kk, jj, k, j, i)] -
                            wcv[idx3(kk, jj, k + 1, j, i)]);
            }
        }
#if defined(_OPENMP)
#pragma omp barrier
#endif
    }
    if (nbot == 8) {
        k = 3;
#if defined(_OPENMP)
#pragma omp for
#endif
        for (intk_c i = 2; i <= ii - 2; ++i) {
            realk_c dxi = dx[idx1(i)];
            for (intk_c j = 3; j <= jj - 2; ++j) {
                realk_c ddyj = ddy[idx1(j)];
                realk_c auz = dxi * ddyj;
                realk_c fub = auz *
                    (wt[idx3(kk, jj, k - 1, j, i)] +
                        wt[idx3(kk, jj, k - 1, j, i + 1)]) *
                    0.5;
                realk_c qkub = 0.5 * fub *
                    (u[idx3(kk, jj, k, j, i)] + u[idx3(kk, jj, k - 1, j, i)]);
                wcu[idx3(kk, jj, k - 1, j, i)] = qkub;
            }
        }
#if defined(_OPENMP)
#pragma omp barrier
#endif
        k = 3;
        i = 2;
        realk_c dxi = dx[idx1(i)];
#if defined(_OPENMP)
#pragma omp for
#endif
        for (intk_c j = 3; j <= jj - 2; j += 2) {
            realk_c ddyj = ddy[idx1(j)];
            realk_c auz = dxi * ddyj;
            realk_c fub = auz *
                (2 * wt[idx3(kk, jj, k - 1, j, i)] +
                    2 * wt[idx3(kk, jj, k - 1, j + 1, i)] +
                    wt[idx3(kk, jj, k - 1, j, i + 1)] +
                    wt[idx3(kk, jj, k - 1, j + 1, i + 1)] +
                    wt[idx3(kk, jj, k - 1, j, i + 2)] +
                    wt[idx3(kk, jj, k - 1, j + 1, i + 2)]) *
                0.125;
            realk_c qkub = 0.5 * fub *
                (u[idx3(kk, jj, k, j, i)] + u[idx3(kk, jj, k - 1, j, i)]);
            wcu[idx3(kk, jj, k, j, i)] = qkub;
            qkub = 0.5 * fub *
                (u[idx3(kk, jj, k, j + 1, i)] +
                    u[idx3(kk, jj, k - 1, j + 1, i)]);
            wcu[idx3(kk, jj, k, j + 1, i)] = qkub;
        }
#if defined(_OPENMP)
#pragma omp barrier
#endif
#if defined(_OPENMP)
#pragma omp for
#endif
        for (intk_c i = 4; i <= ii - 4; i += 2) {
            dxi = dx[idx1(i)];
            for (intk_c j = 3; j <= jj - 2; j += 2) {
                realk_c ddyj = ddy[idx1(j)];
                realk_c auz = dxi * ddyj;
                realk_c fub = auz *
                    (wt[idx3(kk, jj, k - 1, j, i - 1)] +
                        wt[idx3(kk, jj, k - 1, j + 1, i - 1)] +
                        wt[idx3(kk, jj, k - 1, j, i)] +
                        wt[idx3(kk, jj, k - 1, j + 1, i)] +
                        wt[idx3(kk, jj, k - 1, j, i + 1)] +
                        wt[idx3(kk, jj, k - 1, j + 1, i + 1)] +
                        wt[idx3(kk, jj, k - 1, j, i + 2)] +
                        wt[idx3(kk, jj, k - 1, j + 1, i + 2)]) *
                    0.125;
                realk_c qkub = 0.5 * fub *
                    (u[idx3(kk, jj, k, j, i)] + u[idx3(kk, jj, k - 1, j, i)]);
                wcu[idx3(kk, jj, k, j, i)] = qkub;
                qkub = 0.5 * fub *
                    (u[idx3(kk, jj, k, j + 1, i)] +
                        u[idx3(kk, jj, k - 1, j + 1, i)]);
                wcu[idx3(kk, jj, k, j + 1, i)] = qkub;
            }
        }
#if defined(_OPENMP)
#pragma omp barrier
#endif
        i = ii - 2;
        dxi = dx[idx1(i)];
#if defined(_OPENMP)
#pragma omp for
#endif
        for (intk_c j = 3; j <= jj - 2; j += 2) {
            realk_c ddyj = ddy[idx1(j)];
            realk_c auz = dxi * ddyj;
            realk_c fub = auz *
                (wt[idx3(kk, jj, k - 1, j, i - 1)] +
                    wt[idx3(kk, jj, k - 1, j + 1, i - 1)] +
                    wt[idx3(kk, jj, k - 1, j, i)] +
                    wt[idx3(kk, jj, k - 1, j + 1, i)] +
                    2 * wt[idx3(kk, jj, k - 1, j, i + 1)] +
                    2 * wt[idx3(kk, jj, k - 1, j + 1, i + 1)]) *
                0.125;
            realk_c qkub = 0.5 * fub *
                (u[idx3(kk, jj, k, j, i)] + u[idx3(kk, jj, k - 1, j, i)]);
            wcu[idx3(kk, jj, k, j, i)] = qkub;
            qkub = 0.5 * fub *
                (u[idx3(kk, jj, k, j + 1, i)] +
                    u[idx3(kk, jj, k - 1, j + 1, i)]);
            wcu[idx3(kk, jj, k, j + 1, i)] = qkub;
        }
#if defined(_OPENMP)
#pragma omp barrier
#endif
        k = 3;
#if defined(_OPENMP)
#pragma omp for collapse(2)
#endif
        for (intk_c i = 2; i <= ii - 4; i += 2) {
            for (intk_c j = 3; j <= jj - 2; ++j) {
                wcu[idx3(kk, jj, k, j, i + 1)] = 0.5 *
                    (wcu[idx3(kk, jj, k, j, i)] +
                        wcu[idx3(kk, jj, k, j, i + 2)]);
            }
        }
#if defined(_OPENMP)
#pragma omp barrier
#endif
        realk_c rddzk = rddz[idx1(k)];
#if defined(_OPENMP)
#pragma omp for collapse(2)
#endif
        for (intk_c i = 2; i <= ii - 2; ++i) {
            for (intk_c j = 3; j <= jj - 2; ++j) {
                realk_c fkdtu = -1.0 * rddy[idx1(j)] * rdx[idx1(i)] * wkon;
                uo[idx3(kk, jj, k, j, i)] = uo[idx3(kk, jj, k, j, i)] +
                    fkdtu * rddzk *
                        (-wcu[idx3(kk, jj, k, j, i)] +
                            wcu[idx3(kk, jj, k - 1, j, i)]);
            }
        }
#if defined(_OPENMP)
#pragma omp barrier
#endif
        k = 3;
#if defined(_OPENMP)
#pragma omp for
#endif
        for (intk_c i = 3; i <= ii - 2; ++i) {
            realk_c ddxi = ddx[idx1(i)];
            for (intk_c j = 2; j <= jj - 2; ++j) {
                realk_c dyj = dy[idx1(j)];
                realk_c avz = ddxi * dyj;
                realk_c fvb = avz *
                    (wt[idx3(kk, jj, k - 1, j, i)] +
                        wt[idx3(kk, jj, k - 1, j + 1, i)]) *
                    0.5;
                realk_c qkvb = 0.5 * fvb *
                    (v[idx3(kk, jj, k, j, i)] + v[idx3(kk, jj, k - 1, j, i)]);
                wcv[idx3(kk, jj, k - 1, j, i)] = qkvb;
            }
        }
#if defined(_OPENMP)
#pragma omp barrier
#endif
        k = 3;
#if defined(_OPENMP)
#pragma omp for
#endif
        for (intk_c i = 3; i <= ii - 2; i += 2) {
            realk_c ddxi = ddx[idx1(i)];
            j = 2;
            realk_c dyj = dy[idx1(j)];
            realk_c avz = ddxi * dyj;
            realk_c fvb = avz *
                (2 * wt[idx3(kk, jj, k - 1, j, i)] +
                    2 * wt[idx3(kk, jj, k - 1, j, i + 1)] +
                    wt[idx3(kk, jj, k - 1, j + 1, i)] +
                    wt[idx3(kk, jj, k - 1, j + 1, i + 1)] +
                    wt[idx3(kk, jj, k - 1, j + 2, i)] +
                    wt[idx3(kk, jj, k - 1, j + 2, i + 1)]) *
                0.125;
            realk_c qkvb = 0.5 * fvb *
                (v[idx3(kk, jj, k, j, i)] + v[idx3(kk, jj, k - 1, j, i)]);
            wcv[idx3(kk, jj, k, j, i)] = qkvb;
            qkvb = 0.5 * fvb *
                (v[idx3(kk, jj, k, j, i + 1)] +
                    v[idx3(kk, jj, k - 1, j, i + 1)]);
            wcv[idx3(kk, jj, k, j, i + 1)] = qkvb;
            for (intk_c j = 4; j <= jj - 4; j += 2) {
                dyj = dy[idx1(j)];
                avz = ddxi * dyj;
                fvb = avz *
                    (wt[idx3(kk, jj, k - 1, j - 1, i)] +
                        wt[idx3(kk, jj, k - 1, j - 1, i + 1)] +
                        wt[idx3(kk, jj, k - 1, j, i)] +
                        wt[idx3(kk, jj, k - 1, j, i + 1)] +
                        wt[idx3(kk, jj, k - 1, j + 1, i)] +
                        wt[idx3(kk, jj, k - 1, j + 1, i + 1)] +
                        wt[idx3(kk, jj, k - 1, j + 2, i)] +
                        wt[idx3(kk, jj, k - 1, j + 2, i + 1)]) *
                    0.125;
                qkvb = 0.5 * fvb *
                    (v[idx3(kk, jj, k, j, i)] + v[idx3(kk, jj, k - 1, j, i)]);
                wcv[idx3(kk, jj, k, j, i)] = qkvb;
                qkvb = 0.5 * fvb *
                    (v[idx3(kk, jj, k, j, i + 1)] +
                        v[idx3(kk, jj, k - 1, j, i + 1)]);
                wcv[idx3(kk, jj, k, j, i + 1)] = qkvb;
            }
            j = jj - 2;
            dyj = dy[idx1(j)];
            avz = ddxi * dyj;
            fvb = avz *
                (wt[idx3(kk, jj, k - 1, j - 1, i)] +
                    wt[idx3(kk, jj, k - 1, j - 1, i + 1)] +
                    wt[idx3(kk, jj, k - 1, j, i)] +
                    wt[idx3(kk, jj, k - 1, j, i + 1)] +
                    2 * wt[idx3(kk, jj, k - 1, j + 1, i)] +
                    2 * wt[idx3(kk, jj, k - 1, j + 1, i + 1)]) *
                0.125;
            qkvb = 0.5 * fvb *
                (v[idx3(kk, jj, k, j, i)] + v[idx3(kk, jj, k - 1, j, i)]);
            wcv[idx3(kk, jj, k, j, i)] = qkvb;
            qkvb = 0.5 * fvb *
                (v[idx3(kk, jj, k, j, i + 1)] +
                    v[idx3(kk, jj, k - 1, j, i + 1)]);
            wcv[idx3(kk, jj, k, j, i + 1)] = qkvb;
        }
#if defined(_OPENMP)
#pragma omp barrier
#endif
        k = 3;
#if defined(_OPENMP)
#pragma omp for collapse(2)
#endif
        for (intk_c i = 3; i <= ii - 2; ++i) {
            for (intk_c j = 2; j <= jj - 4; j += 2) {
                wcv[idx3(kk, jj, k, j + 1, i)] = 0.5 *
                    (wcv[idx3(kk, jj, k, j, i)] +
                        wcv[idx3(kk, jj, k, j + 2, i)]);
            }
        }
#if defined(_OPENMP)
#pragma omp barrier
#endif
        rddzk = rddz[idx1(k)];
#if defined(_OPENMP)
#pragma omp for collapse(2)
#endif
        for (intk_c i = 3; i <= ii - 2; ++i) {
            for (intk_c j = 2; j <= jj - 2; ++j) {
                realk_c fkdtv = -1.0 * rddx[idx1(i)] * rdy[idx1(j)] * wkon;
                vo[idx3(kk, jj, k, j, i)] = vo[idx3(kk, jj, k, j, i)] +
                    fkdtv * rddzk *
                        (-wcv[idx3(kk, jj, k, j, i)] +
                            wcv[idx3(kk, jj, k - 1, j, i)]);
            }
        }
#if defined(_OPENMP)
#pragma omp barrier
#endif
    }
    if (nlft == 8) {
        j = jj - 2;
#if defined(_OPENMP)
#pragma omp for
#endif
        for (intk_c i = 2; i <= ii - 2; ++i) {
            realk_c dxi = dx[idx1(i)];
            for (intk_c k = 3; k <= kk - 2; ++k) {
                realk_c ddzk = ddz[idx1(k)];
                realk_c auy = dxi * ddzk;
                realk_c fun = auy *
                    (vt[idx3(kk, jj, k, j, i)] +
                        vt[idx3(kk, jj, k, j, i + 1)]) *
                    0.5;
                realk_c qkun = 0.5 * fun *
                    (u[idx3(kk, jj, k, j, i)] + u[idx3(kk, jj, k, j + 1, i)]);
                wcu[idx3(kk, jj, k, j + 1, i)] = qkun;
            }
        }
#if defined(_OPENMP)
#pragma omp barrier
#endif
        j = jj - 2;
        i = 2;
        realk_c dxi = dx[idx1(i)];
#if defined(_OPENMP)
#pragma omp for
#endif
        for (intk_c k = 3; k <= kk - 2; k += 2) {
            realk_c ddzk = ddz[idx1(k)];
            realk_c auy = dxi * ddzk;
            realk_c fun = auy *
                (2 * vt[idx3(kk, jj, k, j, i)] +
                    2 * vt[idx3(kk, jj, k + 1, j, i)] +
                    vt[idx3(kk, jj, k, j, i + 1)] +
                    vt[idx3(kk, jj, k + 1, j, i + 1)] +
                    vt[idx3(kk, jj, k, j, i + 2)] +
                    vt[idx3(kk, jj, k + 1, j, i + 2)]) *
                0.125;
            realk_c qkun = 0.5 * fun *
                (u[idx3(kk, jj, k, j, i)] + u[idx3(kk, jj, k, j + 1, i)]);
            wcu[idx3(kk, jj, k, j, i)] = qkun;
            qkun = 0.5 * fun *
                (u[idx3(kk, jj, k + 1, j, i)] +
                    u[idx3(kk, jj, k + 1, j + 1, i)]);
            wcu[idx3(kk, jj, k + 1, j, i)] = qkun;
        }
#if defined(_OPENMP)
#pragma omp barrier
#endif
#if defined(_OPENMP)
#pragma omp for
#endif
        for (intk_c i = 4; i <= ii - 4; i += 2) {
            dxi = dx[idx1(i)];
            for (intk_c k = 3; k <= kk - 2; k += 2) {
                realk_c ddzk = ddz[idx1(k)];
                realk_c auy = dxi * ddzk;
                realk_c fun = auy *
                    (vt[idx3(kk, jj, k, j, i - 1)] +
                        vt[idx3(kk, jj, k + 1, j, i - 1)] +
                        vt[idx3(kk, jj, k, j, i)] +
                        vt[idx3(kk, jj, k + 1, j, i)] +
                        vt[idx3(kk, jj, k, j, i + 1)] +
                        vt[idx3(kk, jj, k + 1, j, i + 1)] +
                        vt[idx3(kk, jj, k, j, i + 2)] +
                        vt[idx3(kk, jj, k + 1, j, i + 2)]) *
                    0.125;
                realk_c qkun = 0.5 * fun *
                    (u[idx3(kk, jj, k, j, i)] + u[idx3(kk, jj, k, j + 1, i)]);
                wcu[idx3(kk, jj, k, j, i)] = qkun;
                qkun = 0.5 * fun *
                    (u[idx3(kk, jj, k + 1, j, i)] +
                        u[idx3(kk, jj, k + 1, j + 1, i)]);
                wcu[idx3(kk, jj, k + 1, j, i)] = qkun;
            }
        }
#if defined(_OPENMP)
#pragma omp barrier
#endif
        i = ii - 2;
        dxi = dx[idx1(i)];
#if defined(_OPENMP)
#pragma omp for
#endif
        for (intk_c k = 3; k <= kk - 2; k += 2) {
            realk_c ddzk = ddz[idx1(k)];
            realk_c auy = dxi * ddzk;
            realk_c fun = auy *
                (vt[idx3(kk, jj, k, j, i - 1)] +
                    vt[idx3(kk, jj, k + 1, j, i - 1)] +
                    vt[idx3(kk, jj, k, j, i)] + vt[idx3(kk, jj, k + 1, j, i)] +
                    2 * vt[idx3(kk, jj, k, j, i + 1)] +
                    2 * vt[idx3(kk, jj, k + 1, j, i + 1)]) *
                0.125;
            realk_c qkun = 0.5 * fun *
                (u[idx3(kk, jj, k, j, i)] + u[idx3(kk, jj, k, j + 1, i)]);
            wcu[idx3(kk, jj, k, j, i)] = qkun;
            qkun = 0.5 * fun *
                (u[idx3(kk, jj, k + 1, j, i)] +
                    u[idx3(kk, jj, k + 1, j + 1, i)]);
            wcu[idx3(kk, jj, k + 1, j, i)] = qkun;
        }
#if defined(_OPENMP)
#pragma omp barrier
#endif
        j = jj - 2;
#if defined(_OPENMP)
#pragma omp for collapse(2)
#endif
        for (intk_c i = 2; i <= ii - 4; i += 2) {
            for (intk_c k = 3; k <= kk - 2; ++k) {
                wcu[idx3(kk, jj, k, j, i + 1)] = 0.5 *
                    (wcu[idx3(kk, jj, k, j, i)] +
                        wcu[idx3(kk, jj, k, j, i + 2)]);
            }
        }
#if defined(_OPENMP)
#pragma omp barrier
#endif
        j = jj - 2;
#if defined(_OPENMP)
#pragma omp for
#endif
        for (intk_c i = 2; i <= ii - 2; ++i) {
            realk_c fkdtu = -1.0 * rddy[idx1(j)] * rdx[idx1(i)] * wkon;
            for (intk_c k = 3; k <= kk - 2; ++k) {
                realk_c rddzk = rddz[idx1(k)];
                uo[idx3(kk, jj, k, j, i)] = uo[idx3(kk, jj, k, j, i)] +
                    fkdtu * rddzk *
                        (wcu[idx3(kk, jj, k, j, i)] -
                            wcu[idx3(kk, jj, k, j + 1, i)]);
            }
        }
#if defined(_OPENMP)
#pragma omp barrier
#endif
        j = jj - 2;
#if defined(_OPENMP)
#pragma omp for
#endif
        for (intk_c i = 3; i <= ii - 2; ++i) {
            realk_c ddxi = ddx[idx1(i)];
            for (intk_c k = 2; k <= kk - 2; ++k) {
                realk_c dzk = dz[idx1(k)];
                realk_c awy = ddxi * dzk;
                realk_c fwn = awy *
                    (vt[idx3(kk, jj, k, j, i)] +
                        vt[idx3(kk, jj, k + 1, j, i)]) *
                    0.5;
                realk_c qkwn = 0.5 * fwn *
                    (w[idx3(kk, jj, k, j, i)] + w[idx3(kk, jj, k, j + 1, i)]);
                wcw[idx3(kk, jj, k, j + 1, i)] = qkwn;
            }
        }
#if defined(_OPENMP)
#pragma omp barrier
#endif
        j = jj - 2;
#if defined(_OPENMP)
#pragma omp for
#endif
        for (intk_c i = 3; i <= ii - 2; i += 2) {
            realk_c ddxi = ddx[idx1(i)];
            k = 2;
            realk_c dzk = dz[idx1(k)];
            realk_c awy = ddxi * dzk;
            realk_c fwn = awy *
                (2 * vt[idx3(kk, jj, k, j, i)] +
                    2 * vt[idx3(kk, jj, k, j, i + 1)] +
                    vt[idx3(kk, jj, k + 1, j, i)] +
                    vt[idx3(kk, jj, k + 1, j, i + 1)] +
                    vt[idx3(kk, jj, k + 2, j, i)] +
                    vt[idx3(kk, jj, k + 2, j, i + 1)]) *
                0.125;
            realk_c qkwn = 0.5 * fwn *
                (w[idx3(kk, jj, k, j, i)] + w[idx3(kk, jj, k, j + 1, i)]);
            wcw[idx3(kk, jj, k, j, i)] = qkwn;
            qkwn = 0.5 * fwn *
                (w[idx3(kk, jj, k, j, i + 1)] +
                    w[idx3(kk, jj, k, j + 1, i + 1)]);
            wcw[idx3(kk, jj, k, j, i + 1)] = qkwn;
            for (intk_c k = 4; k <= kk - 4; ++k) {
                dzk = dz[idx1(k)];
                awy = ddxi * dzk;
                fwn = awy *
                    (vt[idx3(kk, jj, k - 1, j, i)] +
                        vt[idx3(kk, jj, k - 1, j, i + 1)] +
                        vt[idx3(kk, jj, k, j, i)] +
                        vt[idx3(kk, jj, k, j, i + 1)] +
                        vt[idx3(kk, jj, k + 1, j, i)] +
                        vt[idx3(kk, jj, k + 1, j, i + 1)] +
                        vt[idx3(kk, jj, k + 2, j, i)] +
                        vt[idx3(kk, jj, k + 2, j, i + 1)]) *
                    0.125;
                qkwn = 0.5 * fwn *
                    (w[idx3(kk, jj, k, j, i)] + w[idx3(kk, jj, k, j + 1, i)]);
                wcw[idx3(kk, jj, k, j, i)] = qkwn;
                qkwn = 0.5 * fwn *
                    (w[idx3(kk, jj, k, j, i + 1)] +
                        w[idx3(kk, jj, k, j + 1, i + 1)]);
                wcw[idx3(kk, jj, k, j, i + 1)] = qkwn;
            }
            k = kk - 2;
            dzk = dz[idx1(k)];
            awy = ddxi * dzk;
            fwn = awy *
                (vt[idx3(kk, jj, k - 1, j, i)] +
                    vt[idx3(kk, jj, k - 1, j, i + 1)] +
                    vt[idx3(kk, jj, k, j, i)] + vt[idx3(kk, jj, k, j, i + 1)] +
                    2 * vt[idx3(kk, jj, k + 1, j, i)] +
                    2 * vt[idx3(kk, jj, k + 1, j, i + 1)]) *
                0.125;
            qkwn = 0.5 * fwn *
                (w[idx3(kk, jj, k, j, i)] + w[idx3(kk, jj, k, j + 1, i)]);
            wcw[idx3(kk, jj, k, j, i)] = qkwn;
            qkwn = 0.5 * fwn *
                (w[idx3(kk, jj, k, j, i + 1)] +
                    w[idx3(kk, jj, k, j + 1, i + 1)]);
            wcw[idx3(kk, jj, k, j, i + 1)] = qkwn;
        }
#if defined(_OPENMP)
#pragma omp barrier
#endif
        j = jj - 2;
#if defined(_OPENMP)
#pragma omp for collapse(2)
#endif
        for (intk_c i = 3; i <= ii - 2; ++i) {
            for (intk_c k = 2; k <= kk - 4; k += 2) {
                wcw[idx3(kk, jj, k + 1, j, i)] = 0.5 *
                    (wcw[idx3(kk, jj, k, j, i)] +
                        wcw[idx3(kk, jj, k + 2, j, i)]);
            }
        }
#if defined(_OPENMP)
#pragma omp barrier
#endif
        j = jj - 2;
#if defined(_OPENMP)
#pragma omp for
#endif
        for (intk_c i = 3; i <= ii - 2; ++i) {
            realk_c fkdtw = -1.0 * rddx[idx1(i)] * rddy[idx1(j)] * wkon;
            for (intk_c k = 2; k <= kk - 2; ++k) {
                realk_c rdzk = rdz[idx1(k)];
                wo[idx3(kk, jj, k, j, i)] = wo[idx3(kk, jj, k, j, i)] +
                    fkdtw * rdzk *
                        (-wcw[idx3(kk, jj, k, j + 1, i)] +
                            wcw[idx3(kk, jj, k, j, i)]);
            }
        }
#if defined(_OPENMP)
#pragma omp barrier
#endif
    }
    if (nrgt == 8) {
        j = 3;
#if defined(_OPENMP)
#pragma omp for
#endif
        for (intk_c i = 2; i <= ii - 2; ++i) {
            realk_c dxi = dx[idx1(i)];
            for (intk_c k = 3; k <= kk - 2; ++k) {
                realk_c ddzk = ddz[idx1(k)];
                realk_c auy = dxi * ddzk;
                realk_c fus = auy *
                    (vt[idx3(kk, jj, k, j - 1, i)] +
                        vt[idx3(kk, jj, k, j - 1, i + 1)]) *
                    0.5;
                realk_c qkus = 0.5 * fus *
                    (u[idx3(kk, jj, k, j, i)] + u[idx3(kk, jj, k, j - 1, i)]);
                wcu[idx3(kk, jj, k, j - 1, i)] = qkus;
            }
        }
#if defined(_OPENMP)
#pragma omp barrier
#endif
        j = 3;
        i = 2;
        realk_c dxi = dx[idx1(i)];
#if defined(_OPENMP)
#pragma omp for
#endif
        for (intk_c k = 3; k <= kk - 2; k += 2) {
            realk_c ddzk = ddz[idx1(k)];
            realk_c auy = dxi * ddzk;
            realk_c fus = auy *
                (2 * vt[idx3(kk, jj, k, j - 1, i)] +
                    2 * vt[idx3(kk, jj, k + 1, j - 1, i)] +
                    vt[idx3(kk, jj, k, j - 1, i + 1)] +
                    vt[idx3(kk, jj, k + 1, j - 1, i + 1)] +
                    vt[idx3(kk, jj, k, j - 1, i + 2)] +
                    vt[idx3(kk, jj, k + 1, j - 1, i + 2)]) *
                0.125;
            realk_c qkus = 0.5 * fus *
                (u[idx3(kk, jj, k, j, i)] + u[idx3(kk, jj, k, j - 1, i)]);
            wcu[idx3(kk, jj, k, j, i)] = qkus;
            qkus = 0.5 * fus *
                (u[idx3(kk, jj, k + 1, j, i)] +
                    u[idx3(kk, jj, k + 1, j - 1, i)]);
            wcu[idx3(kk, jj, k + 1, j, i)] = qkus;
        }
#if defined(_OPENMP)
#pragma omp barrier
#endif
#if defined(_OPENMP)
#pragma omp for
#endif
        for (intk_c i = 4; i <= ii - 4; i += 2) {
            dxi = dx[idx1(i)];
            for (intk_c k = 3; k <= kk - 2; k += 2) {
                realk_c ddzk = ddz[idx1(k)];
                realk_c auy = dxi * ddzk;
                realk_c fus = auy *
                    (vt[idx3(kk, jj, k, j - 1, i - 1)] +
                        vt[idx3(kk, jj, k + 1, j - 1, i - 1)] +
                        vt[idx3(kk, jj, k, j - 1, i)] +
                        vt[idx3(kk, jj, k + 1, j - 1, i)] +
                        vt[idx3(kk, jj, k, j - 1, i + 1)] +
                        vt[idx3(kk, jj, k + 1, j - 1, i + 1)] +
                        vt[idx3(kk, jj, k, j - 1, i + 2)] +
                        vt[idx3(kk, jj, k + 1, j - 1, i + 2)]) *
                    0.125;
                realk_c qkus = 0.5 * fus *
                    (u[idx3(kk, jj, k, j, i)] + u[idx3(kk, jj, k, j - 1, i)]);
                wcu[idx3(kk, jj, k, j, i)] = qkus;
                qkus = 0.5 * fus *
                    (u[idx3(kk, jj, k + 1, j, i)] +
                        u[idx3(kk, jj, k + 1, j - 1, i)]);
                wcu[idx3(kk, jj, k + 1, j, i)] = qkus;
            }
        }
#if defined(_OPENMP)
#pragma omp barrier
#endif
        i = ii - 2;
        dxi = dx[idx1(i)];
#if defined(_OPENMP)
#pragma omp for
#endif
        for (intk_c k = 3; k <= kk - 2; k += 2) {
            realk_c ddzk = ddz[idx1(k)];
            realk_c auy = dxi * ddzk;
            realk_c fus = auy *
                (vt[idx3(kk, jj, k, j - 1, i - 1)] +
                    vt[idx3(kk, jj, k + 1, j - 1, i - 1)] +
                    vt[idx3(kk, jj, k, j - 1, i)] +
                    vt[idx3(kk, jj, k + 1, j - 1, i)] +
                    2 * vt[idx3(kk, jj, k, j - 1, i + 1)] +
                    2 * vt[idx3(kk, jj, k + 1, j - 1, i + 1)]) *
                0.125;
            realk_c qkus = 0.5 * fus *
                (u[idx3(kk, jj, k, j, i)] + u[idx3(kk, jj, k, j - 1, i)]);
            wcu[idx3(kk, jj, k, j, i)] = qkus;
            qkus = 0.5 * fus *
                (u[idx3(kk, jj, k + 1, j, i)] +
                    u[idx3(kk, jj, k + 1, j - 1, i)]);
            wcu[idx3(kk, jj, k + 1, j, i)] = qkus;
        }
#if defined(_OPENMP)
#pragma omp barrier
#endif
        j = 3;
#if defined(_OPENMP)
#pragma omp for collapse(2)
#endif
        for (intk_c i = 2; i <= ii - 4; i += 2) {
            for (intk_c k = 3; k <= kk - 2; ++k) {
                wcu[idx3(kk, jj, k, j, i + 1)] = 0.5 *
                    (wcu[idx3(kk, jj, k, j, i)] +
                        wcu[idx3(kk, jj, k, j, i + 2)]);
            }
        }
#if defined(_OPENMP)
#pragma omp barrier
#endif
        j = 3;
#if defined(_OPENMP)
#pragma omp for
#endif
        for (intk_c i = 2; i <= ii - 2; ++i) {
            realk_c fkdtu = -1.0 * rddy[idx1(j)] * rdx[idx1(i)] * wkon;
            for (intk_c k = 3; k <= kk - 2; ++k) {
                realk_c rddzk = rddz[idx1(k)];
                uo[idx3(kk, jj, k, j, i)] = uo[idx3(kk, jj, k, j, i)] +
                    fkdtu * rddzk *
                        (-wcu[idx3(kk, jj, k, j, i)] +
                            wcu[idx3(kk, jj, k, j - 1, i)]);
            }
        }
#if defined(_OPENMP)
#pragma omp barrier
#endif
        j = 3;
#if defined(_OPENMP)
#pragma omp for
#endif
        for (intk_c i = 3; i <= ii - 2; ++i) {
            realk_c ddxi = ddx[idx1(i)];
            for (intk_c k = 2; k <= kk - 2; ++k) {
                realk_c dzk = dz[idx1(k)];
                realk_c awy = ddxi * dzk;
                realk_c fws = awy *
                    (vt[idx3(kk, jj, k, j - 1, i)] +
                        vt[idx3(kk, jj, k + 1, j - 1, i)]) *
                    0.5;
                realk_c qkws = 0.5 * fws *
                    (w[idx3(kk, jj, k, j, i)] + w[idx3(kk, jj, k, j - 1, i)]);
                wcw[idx3(kk, jj, k, j - 1, i)] = qkws;
            }
        }
#if defined(_OPENMP)
#pragma omp barrier
#endif
        j = 3;
#if defined(_OPENMP)
#pragma omp for
#endif
        for (intk_c i = 3; i <= ii - 2; i += 2) {
            realk_c ddxi = ddx[idx1(i)];
            k = 2;
            realk_c dzk = dz[idx1(k)];
            realk_c awy = ddxi * dzk;
            realk_c fws = awy *
                (2 * vt[idx3(kk, jj, k, j - 1, i)] +
                    2 * vt[idx3(kk, jj, k, j - 1, i + 1)] +
                    vt[idx3(kk, jj, k + 1, j - 1, i)] +
                    vt[idx3(kk, jj, k + 1, j - 1, i + 1)] +
                    vt[idx3(kk, jj, k + 2, j - 1, i)] +
                    vt[idx3(kk, jj, k + 2, j - 1, i + 1)]) *
                0.125;
            realk_c qkws = 0.5 * fws *
                (w[idx3(kk, jj, k, j, i)] + w[idx3(kk, jj, k, j - 1, i)]);
            wcw[idx3(kk, jj, k, j, i)] = qkws;
            qkws = 0.5 * fws *
                (w[idx3(kk, jj, k, j, i + 1)] +
                    w[idx3(kk, jj, k, j - 1, i + 1)]);
            wcw[idx3(kk, jj, k, j, i + 1)] = qkws;
            for (intk_c k = 4; k <= kk - 4; ++k) {
                dzk = dz[idx1(k)];
                awy = ddxi * dzk;
                fws = awy *
                    (vt[idx3(kk, jj, k - 1, j - 1, i)] +
                        vt[idx3(kk, jj, k - 1, j - 1, i + 1)] +
                        vt[idx3(kk, jj, k, j - 1, i)] +
                        vt[idx3(kk, jj, k, j - 1, i + 1)] +
                        vt[idx3(kk, jj, k + 1, j - 1, i)] +
                        vt[idx3(kk, jj, k + 1, j - 1, i + 1)] +
                        vt[idx3(kk, jj, k + 2, j - 1, i)] +
                        vt[idx3(kk, jj, k + 2, j - 1, i + 1)]) *
                    0.125;
                qkws = 0.5 * fws *
                    (w[idx3(kk, jj, k, j, i)] + w[idx3(kk, jj, k, j - 1, i)]);
                wcw[idx3(kk, jj, k, j, i)] = qkws;
                qkws = 0.5 * fws *
                    (w[idx3(kk, jj, k, j, i + 1)] +
                        w[idx3(kk, jj, k, j - 1, i + 1)]);
                wcw[idx3(kk, jj, k, j, i + 1)] = qkws;
            }
            k = kk - 2;
            dzk = dz[idx1(k)];
            awy = ddxi * dzk;
            fws = awy *
                (vt[idx3(kk, jj, k - 1, j - 1, i)] +
                    vt[idx3(kk, jj, k - 1, j - 1, i + 1)] +
                    vt[idx3(kk, jj, k, j - 1, i)] +
                    vt[idx3(kk, jj, k, j - 1, i + 1)] +
                    2 * vt[idx3(kk, jj, k + 1, j - 1, i)] +
                    2 * vt[idx3(kk, jj, k + 1, j - 1, i + 1)]) *
                0.125;
            qkws = 0.5 * fws *
                (w[idx3(kk, jj, k, j, i)] + w[idx3(kk, jj, k, j - 1, i)]);
            wcw[idx3(kk, jj, k, j, i)] = qkws;
            qkws = 0.5 * fws *
                (w[idx3(kk, jj, k, j, i + 1)] +
                    w[idx3(kk, jj, k, j - 1, i + 1)]);
            wcw[idx3(kk, jj, k, j, i + 1)] = qkws;
        }
#if defined(_OPENMP)
#pragma omp barrier
#endif
        j = 3;
#if defined(_OPENMP)
#pragma omp for collapse(2)
#endif
        for (intk_c i = 3; i <= ii - 2; ++i) {
            for (intk_c k = 2; k <= kk - 4; k += 2) {
                wcw[idx3(kk, jj, k + 1, j, i)] = 0.5 *
                    (wcw[idx3(kk, jj, k, j, i)] +
                        wcw[idx3(kk, jj, k + 2, j, i)]);
            }
        }
#if defined(_OPENMP)
#pragma omp barrier
#endif
        j = 3;
#if defined(_OPENMP)
#pragma omp for
#endif
        for (intk_c i = 3; i <= ii - 2; ++i) {
            realk_c fkdtw = -1.0 * rddx[idx1(i)] * rddy[idx1(j)] * wkon;
            for (intk_c k = 2; k <= kk - 2; ++k) {
                realk_c rdzk = rdz[idx1(k)];
                wo[idx3(kk, jj, k, j, i)] = wo[idx3(kk, jj, k, j, i)] +
                    fkdtw * rdzk *
                        (wcw[idx3(kk, jj, k, j - 1, i)] -
                            wcw[idx3(kk, jj, k, j, i)]);
            }
        }
    }
}

// Convective contribution for v-momentum.
void tstle4_kon_v_c(const intk_c kk, const intk_c jj, const intk_c ii,
    realk_c *vo, const realk_c *u, const realk_c *v, const realk_c *w,
    const realk_c *ut, const realk_c *vt, const realk_c *wt, const realk_c *dx,
    const realk_c *dy, const realk_c *dz, const realk_c *ddx,
    const realk_c *ddy, const realk_c *ddz, const realk_c *rdx,
    const realk_c *rdy, const realk_c *rdz, const realk_c *rddx,
    const realk_c *rddy, const realk_c *rddz, const intk_c nrgt,
    const intk_c nlft) {
    intk_c nrv = 0;
    intk_c nlv = 0;

    (void)dx;
    (void)dz;
    (void)rdx;
    (void)rdz;
    (void)rddy;

    if (nlft == 7)
        nlv = 1;
    if (nrgt == 3)
        nrv = 1;
    if (nlft == 3)
        nlv = 1;

#if defined(_OPENMP)
#pragma omp for collapse(3)
#endif
    for (intk_c i = 3; i <= ii - 2; ++i) {
        for (intk_c j = 3 - nrv; j <= jj - 3 + nlv; ++j) {
            for (intk_c k = 3; k <= kk - 2; ++k) {
                const realk_c ax = dy[idx1(j)] * ddz[idx1(k)];
                const realk_c ay = ddx[idx1(i)] * ddz[idx1(k)];
                const realk_c az = ddx[idx1(i)] * dy[idx1(j)];

                const realk_c fe = ax *
                    (ut[idx3(kk, jj, k, j, i)] +
                        ut[idx3(kk, jj, k, j + 1, i)]) *
                    (realk_c)0.5;
                const realk_c fw = ax *
                    (ut[idx3(kk, jj, k, j, i - 1)] +
                        ut[idx3(kk, jj, k, j + 1, i - 1)]) *
                    (realk_c)0.5;
                const realk_c fn = ay *
                    (vt[idx3(kk, jj, k, j, i)] +
                        (vt[idx3(kk, jj, k, j + 1, i)] -
                            vt[idx3(kk, jj, k, j, i)]) *
                            (realk_c)0.5 * dy[idx1(j)] / ddy[idx1(j + 1)]);
                const realk_c fs = ay *
                    (vt[idx3(kk, jj, k, j - 1, i)] +
                        (vt[idx3(kk, jj, k, j, i)] -
                            vt[idx3(kk, jj, k, j - 1, i)]) *
                            (realk_c)0.5 * dy[idx1(j - 1)] / ddy[idx1(j)]);
                const realk_c ft = az *
                    (wt[idx3(kk, jj, k, j, i)] +
                        wt[idx3(kk, jj, k, j + 1, i)]) *
                    (realk_c)0.5;
                const realk_c fb = az *
                    (wt[idx3(kk, jj, k - 1, j, i)] +
                        wt[idx3(kk, jj, k - 1, j + 1, i)]) *
                    (realk_c)0.5;

                const realk_c qe = (realk_c)0.5 * fe *
                    (v[idx3(kk, jj, k, j, i)] + v[idx3(kk, jj, k, j, i + 1)]);
                const realk_c qw = (realk_c)0.5 * fw *
                    (v[idx3(kk, jj, k, j, i - 1)] + v[idx3(kk, jj, k, j, i)]);
                const realk_c qn = (realk_c)0.5 * fn *
                    (v[idx3(kk, jj, k, j, i)] + v[idx3(kk, jj, k, j + 1, i)]);
                const realk_c qs = (realk_c)0.5 * fs *
                    (v[idx3(kk, jj, k, j - 1, i)] + v[idx3(kk, jj, k, j, i)]);
                const realk_c qt = (realk_c)0.5 * ft *
                    (v[idx3(kk, jj, k, j, i)] + v[idx3(kk, jj, k + 1, j, i)]);
                const realk_c qb = (realk_c)0.5 * fb *
                    (v[idx3(kk, jj, k - 1, j, i)] + v[idx3(kk, jj, k, j, i)]);

                vo[idx3(kk, jj, k, j, i)] = -(qe - qw + qn - qs + qt - qb) *
                    rddx[idx1(i)] * rdy[idx1(j)] * rddz[idx1(k)];
            }
        }
    }
}

// Convective contribution for w-momentum.
void tstle4_kon_w_c(const intk_c kk, const intk_c jj, const intk_c ii,
    realk_c *wo, const realk_c *u, const realk_c *v, const realk_c *w,
    const realk_c *ut, const realk_c *vt, const realk_c *wt, const realk_c *dx,
    const realk_c *dy, const realk_c *dz, const realk_c *ddx,
    const realk_c *ddy, const realk_c *ddz, const realk_c *rdx,
    const realk_c *rdy, const realk_c *rdz, const realk_c *rddx,
    const realk_c *rddy, const realk_c *rddz, const intk_c nbot,
    const intk_c ntop) {
    intk_c nbw = 0;
    intk_c ntw = 0;

    (void)dx;
    (void)dy;
    (void)rdx;
    (void)rdy;
    (void)rddz;

    if (ntop == 7)
        ntw = 1;
    if (nbot == 3)
        nbw = 1;
    if (ntop == 3)
        ntw = 1;

#if defined(_OPENMP)
#pragma omp for collapse(3)
#endif
    for (intk_c i = 3; i <= ii - 2; ++i) {
        for (intk_c j = 3; j <= jj - 2; ++j) {
            for (intk_c k = 3 - nbw; k <= kk - 3 + ntw; ++k) {
                const realk_c ax = ddy[idx1(j)] * dz[idx1(k)];
                const realk_c ay = ddx[idx1(i)] * dz[idx1(k)];
                const realk_c az = ddx[idx1(i)] * ddy[idx1(j)];

                const realk_c fe = ax *
                    (ut[idx3(kk, jj, k, j, i)] +
                        ut[idx3(kk, jj, k + 1, j, i)]) *
                    (realk_c)0.5;
                const realk_c fw = ax *
                    (ut[idx3(kk, jj, k, j, i - 1)] +
                        ut[idx3(kk, jj, k + 1, j, i - 1)]) *
                    (realk_c)0.5;
                const realk_c fn = ay *
                    (vt[idx3(kk, jj, k, j, i)] +
                        vt[idx3(kk, jj, k + 1, j, i)]) *
                    (realk_c)0.5;
                const realk_c fs = ay *
                    (vt[idx3(kk, jj, k, j - 1, i)] +
                        vt[idx3(kk, jj, k + 1, j - 1, i)]) *
                    (realk_c)0.5;
                const realk_c ft = az *
                    (wt[idx3(kk, jj, k, j, i)] +
                        (wt[idx3(kk, jj, k + 1, j, i)] -
                            wt[idx3(kk, jj, k, j, i)]) *
                            (realk_c)0.5 * dz[idx1(k)] / ddz[idx1(k + 1)]);
                const realk_c fb = az *
                    (wt[idx3(kk, jj, k - 1, j, i)] +
                        (wt[idx3(kk, jj, k, j, i)] -
                            wt[idx3(kk, jj, k - 1, j, i)]) *
                            (realk_c)0.5 * dz[idx1(k - 1)] / ddz[idx1(k)]);

                const realk_c qe = (realk_c)0.5 * fe *
                    (w[idx3(kk, jj, k, j, i)] + w[idx3(kk, jj, k, j, i + 1)]);
                const realk_c qw = (realk_c)0.5 * fw *
                    (w[idx3(kk, jj, k, j, i - 1)] + w[idx3(kk, jj, k, j, i)]);
                const realk_c qn = (realk_c)0.5 * fn *
                    (w[idx3(kk, jj, k, j, i)] + w[idx3(kk, jj, k, j + 1, i)]);
                const realk_c qs = (realk_c)0.5 * fs *
                    (w[idx3(kk, jj, k, j - 1, i)] + w[idx3(kk, jj, k, j, i)]);
                const realk_c qt = (realk_c)0.5 * ft *
                    (w[idx3(kk, jj, k, j, i)] + w[idx3(kk, jj, k + 1, j, i)]);
                const realk_c qb = (realk_c)0.5 * fb *
                    (w[idx3(kk, jj, k - 1, j, i)] + w[idx3(kk, jj, k, j, i)]);

                wo[idx3(kk, jj, k, j, i)] = -(qe - qw + qn - qs + qt - qb) *
                    rddx[idx1(i)] * rddy[idx1(j)] * rdz[idx1(k)];
            }
        }
    }
}

// Diffusive contribution for u-momentum with optional LES correction.
void tstle4_diff_u_c(const intk_c kk, const intk_c jj, const intk_c ii,
    realk_c *uo, const realk_c *u, const realk_c *v, const realk_c *w,
    const realk_c *g, const realk_c *dx, const realk_c *dy, const realk_c *dz,
    const realk_c *ddx, const realk_c *ddy, const realk_c *ddz,
    const realk_c *rdx, const realk_c *rdy, const realk_c *rdz,
    const realk_c *rddx, const realk_c *rddy, const realk_c *rddz,
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
                const realk_c qe = -ge * ax * rddx[idx1(i + 1)] *
                    (u[idx3(kk, jj, k, j, i + 1)] - u_kji);
                const realk_c qw = -gw * ax * rddx[idx1(i)] *
                    (u_kji - u[idx3(kk, jj, k, j, i - 1)]);
                const realk_c qn = -gn * ay * rdy[idx1(j)] *
                    (u[idx3(kk, jj, k, j + 1, i)] - u_kji);
                const realk_c qs = -gs * ay * rdy[idx1(j - 1)] *
                    (u_kji - u[idx3(kk, jj, k, j - 1, i)]);
                const realk_c qt = -gt * az * rdz[idx1(k)] *
                    (u[idx3(kk, jj, k + 1, j, i)] - u_kji);
                const realk_c qb = -gb * az * rdz[idx1(k - 1)] *
                    (u_kji - u[idx3(kk, jj, k - 1, j, i)]);

                const realk_c st =
                    ((ge * (u[idx3(kk, jj, k, j, i + 1)] - u_kji) *
                        rddx[idx1(i + 1)]) -
                        (gw * (u_kji - u[idx3(kk, jj, k, j, i - 1)]) *
                            rddx[idx1(i)])) *
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

                const realk_c fak = ((realk_c)1.0 / rho) * rddy[idx1(j)] *
                    rdx[idx1(i)] * rddz[idx1(k)];
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
    const realk_c *rdx, const realk_c *rdy, const realk_c *rdz,
    const realk_c *rddx, const realk_c *rddy, const realk_c *rddz,
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
                const realk_c qe = -ge * ax * rdx[idx1(i)] *
                    (v[idx3(kk, jj, k, j, i + 1)] - v_kji);
                const realk_c qw = -gw * ax * rdx[idx1(i - 1)] *
                    (v_kji - v[idx3(kk, jj, k, j, i - 1)]);
                const realk_c qn = -gn * ay * rddy[idx1(j + 1)] *
                    (v[idx3(kk, jj, k, j + 1, i)] - v_kji);
                const realk_c qs = -gs * ay * rddy[idx1(j)] *
                    (v_kji - v[idx3(kk, jj, k, j - 1, i)]);
                const realk_c qt = -gt * az * rdz[idx1(k)] *
                    (v[idx3(kk, jj, k + 1, j, i)] - v_kji);
                const realk_c qb = -gb * az * rdz[idx1(k - 1)] *
                    (v_kji - v[idx3(kk, jj, k - 1, j, i)]);

                const realk_c st = ((ge *
                                        (u[idx3(kk, jj, k, j + 1, i)] -
                                            u[idx3(kk, jj, k, j, i)])) -
                                    (gw *
                                        (u[idx3(kk, jj, k, j + 1, i - 1)] -
                                            u[idx3(kk, jj, k, j, i - 1)]))) *
                        ddz[idx1(k)] +
                    ((gn * (v[idx3(kk, jj, k, j + 1, i)] - v_kji) *
                        rddy[idx1(j + 1)]) -
                        (gs * (v_kji - v[idx3(kk, jj, k, j - 1, i)]) *
                            rddy[idx1(j)])) *
                        ay +
                    ((gt *
                        (w[idx3(kk, jj, k, j + 1, i)] -
                            w[idx3(kk, jj, k, j, i)])) -
                        (gb *
                            (w[idx3(kk, jj, k - 1, j + 1, i)] -
                                w[idx3(kk, jj, k - 1, j, i)]))) *
                        ddx[idx1(i)];
                const realk_c qc = st * (realk_c)iles;

                const realk_c fak = ((realk_c)1.0 / rho) * rddx[idx1(i)] *
                    rdy[idx1(j)] * rddz[idx1(k)];
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
    const realk_c *rdx, const realk_c *rdy, const realk_c *rdz,
    const realk_c *rddx, const realk_c *rddy, const realk_c *rddz,
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
                const realk_c qe = -ge * ax * rdx[idx1(i)] *
                    (w[idx3(kk, jj, k, j, i + 1)] - w_kji);
                const realk_c qw = -gw * ax * rdx[idx1(i - 1)] *
                    (w_kji - w[idx3(kk, jj, k, j, i - 1)]);
                const realk_c qn = -gn * ay * rdy[idx1(j)] *
                    (w[idx3(kk, jj, k, j + 1, i)] - w_kji);
                const realk_c qs = -gs * ay * rdy[idx1(j - 1)] *
                    (w_kji - w[idx3(kk, jj, k, j - 1, i)]);
                const realk_c qt = -gt * az * rddz[idx1(k + 1)] *
                    (w[idx3(kk, jj, k + 1, j, i)] - w_kji);
                const realk_c qb = -gb * az * rddz[idx1(k)] *
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
                        rddz[idx1(k + 1)]) -
                        (gb * (w_kji - w[idx3(kk, jj, k - 1, j, i)]) *
                            rddz[idx1(k)])) *
                        az;
                const realk_c qc = st * (realk_c)iles;

                const realk_c fak = ((realk_c)1.0 / rho) * rddx[idx1(i)] *
                    rddy[idx1(j)] * rdz[idx1(k)];
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