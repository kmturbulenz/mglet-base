#include "tstle4_kernels_shared.h"

#pragma omp begin declare target

// Convective contribution for u-momentum (skew-symmetric discretization).
void tstle4_kon_u_c(const intk_c kk, const intk_c jj, const intk_c ii,
    realk_c *uo, const realk_c *u, const realk_c *v, const realk_c *w,
    const realk_c *ut, const realk_c *vt, const realk_c *wt,
    const realk_c *dx, const realk_c *dy, const realk_c *dz,
    const realk_c *ddx, const realk_c *ddy, const realk_c *ddz,
    const intk_c nfro, const intk_c nbac)
    {

    intk_c nfu = 0; intk_c nbu = 0;
    (void)dy; (void)dz;

    if (nbac == 7) nbu = 1;
    if (nfro == 3) nfu = 1;
    if (nbac == 3) nbu = 1;

    #pragma omp for collapse(3)
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
                    ((realk_c)1.0 / dx[idx1(i)]) *
                    ((realk_c)1.0 / ddy[idx1(j)]) *
                    ((realk_c)1.0 / ddz[idx1(k)]);
            }
        }
    }
}

// Convective contribution for v-momentum.
void tstle4_kon_v_c(const intk_c kk, const intk_c jj, const intk_c ii,
    realk_c *vo, const realk_c *u, const realk_c *v, const realk_c *w,
    const realk_c *ut, const realk_c *vt, const realk_c *wt,
    const realk_c *dx, const realk_c *dy, const realk_c *dz,
     const realk_c *ddx, const realk_c *ddy, const realk_c *ddz,
    const intk_c nrgt, const intk_c nlft)
    {

    intk_c nrv = 0; intk_c nlv = 0;

    (void)dx; (void)dz;

    if (nlft == 7) nlv = 1;
    if (nrgt == 3) nrv = 1;
    if (nlft == 3) nlv = 1;

    #pragma omp for collapse(3)
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
                    ((realk_c)1.0 / ddx[idx1(i)]) *
                    ((realk_c)1.0 / dy[idx1(j)]) *
                    ((realk_c)1.0 / ddz[idx1(k)]);
            }
        }
    }
}

// Convective contribution for w-momentum.
void tstle4_kon_w_c(const intk_c kk, const intk_c jj, const intk_c ii,
    realk_c *wo, const realk_c *u, const realk_c *v, const realk_c *w,
    const realk_c *ut, const realk_c *vt, const realk_c *wt,
    const realk_c *dx, const realk_c *dy, const realk_c *dz,
    const realk_c *ddx, const realk_c *ddy, const realk_c *ddz,
    const intk_c nbot, const intk_c ntop)
    {

    intk_c nbw = 0; intk_c ntw = 0;

    (void)dx; (void)dy;

    if (ntop == 7) ntw = 1;
    if (nbot == 3) nbw = 1;
    if (ntop == 3) ntw = 1;

    #pragma omp for collapse(3)
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
                    ((realk_c)1.0 / ddx[idx1(i)]) *
                    ((realk_c)1.0 / ddy[idx1(j)]) *
                    ((realk_c)1.0 / dz[idx1(k)]);
            }
        }
    }
}

#pragma omp end declare target


// Full convective implementation: iterate all grids in C using Fortran metadata.
void tstle4_kon_impl_c(realk_c *uo, realk_c *vo, realk_c *wo,
    const realk_c *u, const realk_c *v, const realk_c *w,
    const realk_c *ut, const realk_c *vt, const realk_c *wt,
    const realk_c *dx, const realk_c *dy, const realk_c *dz,
    const realk_c *ddx, const realk_c *ddy, const realk_c *ddz,
    const intk_c nmygrids_in, const tstle4_grid_t *grids)
{
    #pragma omp target teams distribute
    for (intk_c igrid = 0; igrid < nmygrids_in; ++igrid)
    {
        // Obtain memory views for the current grid.
        const tstle4_grid_t *const g = &grids[igrid];
        realk_c *const uo_g = &uo[(size_t)(g->ip3 - 1)];
        realk_c *const vo_g = &vo[(size_t)(g->ip3 - 1)];
        realk_c *const wo_g = &wo[(size_t)(g->ip3 - 1)];
        const realk_c *const u_g = &u[(size_t)(g->ip3 - 1)];
        const realk_c *const v_g = &v[(size_t)(g->ip3 - 1)];
        const realk_c *const w_g = &w[(size_t)(g->ip3 - 1)];
        const realk_c *const ut_g = &ut[(size_t)(g->ip3 - 1)];
        const realk_c *const vt_g = &vt[(size_t)(g->ip3 - 1)];
        const realk_c *const wt_g = &wt[(size_t)(g->ip3 - 1)];
        const realk_c *const dx_g = &dx[(size_t)(g->ipx - 1)];
        const realk_c *const dy_g = &dy[(size_t)(g->ipy - 1)];
        const realk_c *const dz_g = &dz[(size_t)(g->ipz - 1)];
        const realk_c *const ddx_g = &ddx[(size_t)(g->ipx - 1)];
        const realk_c *const ddy_g = &ddy[(size_t)(g->ipy - 1)];
        const realk_c *const ddz_g = &ddz[(size_t)(g->ipz - 1)];

        #pragma omp parallel
        {
            // Compution of uo (momentum update in x-direction)
            tstle4_kon_u_c(g->kk, g->jj, g->ii, uo_g, u_g, v_g, w_g,
                ut_g, vt_g, wt_g, dx_g, dy_g, dz_g, ddx_g, ddy_g, ddz_g,
                g->nfro, g->nbac);

            // Compution of vo (momentum update in y-direction)
            tstle4_kon_v_c(g->kk, g->jj, g->ii, vo_g, u_g, v_g, w_g,
                ut_g, vt_g, wt_g, dx_g, dy_g, dz_g, ddx_g, ddy_g, ddz_g,
                g->nrgt, g->nlft);

            // Compution of wo (momentum update in z-direction)
            tstle4_kon_w_c(g->kk, g->jj, g->ii, wo_g, u_g, v_g, w_g,
                ut_g, vt_g, wt_g, dx_g, dy_g, dz_g, ddx_g, ddy_g, ddz_g,
                g->nbot, g->ntop);
        }
    }
}
