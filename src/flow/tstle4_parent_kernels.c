#include "tstle4_kernels_shared.h"

#pragma omp begin declare target

void tstle4_par_c(const intk_c kk, const intk_c jj, const intk_c ii,
    realk_c *uo, realk_c *vo, realk_c *wo, const realk_c *u, const realk_c *v,
    const realk_c *w, const realk_c *ut, const realk_c *vt, const realk_c *wt,
    const realk_c *dx, const realk_c *dy, const realk_c *dz, const realk_c *ddx,
    const realk_c *ddy, const realk_c *ddz, realk_c *wcu, realk_c *wcv,
    realk_c *wcw, const intk_c nfro, const intk_c nbac, const intk_c nrgt,
    const intk_c nlft, const intk_c nbot, const intk_c ntop) {

// Definition to abbreviate the inverse of the grid spacing
#define INV_DX(i) ((realk_c)1.0 / dx[idx1(i)])
#define INV_DY(j) ((realk_c)1.0 / dy[idx1(j)])
#define INV_DZ(k) ((realk_c)1.0 / dz[idx1(k)])
#define INV_DDX(i) ((realk_c)1.0 / ddx[idx1(i)])
#define INV_DDY(j) ((realk_c)1.0 / ddy[idx1(j)])
#define INV_DDZ(k) ((realk_c)1.0 / ddz[idx1(k)])

    const realk_c wkon = (realk_c)1.0;
    intk_c i, j, k;

    // Local upwind correction at front parent boundary.
    if (nfro == 8) {
        i = 4;
#if defined(_OPENMP)
#pragma omp for
#endif
        for (intk_c j = 3; j <= jj - 2; ++j) {
            realk_c dyj = dy[idx1(j)];
            realk_c ddyj = ddy[idx1(j)];
            realk_c fkdtu = -1.0 * INV_DDY(j) * INV_DX(i) * wkon;
            realk_c fkdtv = -1.0 * INV_DDX(i) * INV_DY(j) * wkon;
            realk_c fkdtw = -1.0 * INV_DDX(i) * INV_DDY(j) * wkon;
            for (intk_c k = 3; k <= kk - 2; ++k) {
                realk_c dzk = dz[idx1(k)];
                realk_c rdzk = INV_DZ(k);
                realk_c ddzk = ddz[idx1(k)];
                realk_c rddzk = INV_DDZ(k);
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
    // Local upwind correction at back parent boundary.
    if (nbac == 8) {
        i = ii - 2;
#if defined(_OPENMP)
#pragma omp for
#endif
        for (intk_c j = 3; j <= jj - 2; ++j) {
            realk_c dyj = dy[idx1(j)];
            realk_c ddyj = ddy[idx1(j)];
            realk_c fkdtu = -1.0 * INV_DDY(j) * INV_DX(i) * wkon;
            realk_c fkdtv = -1.0 * INV_DDX(i) * INV_DY(j) * wkon;
            realk_c fkdtw = -1.0 * INV_DDX(i) * INV_DDY(j) * wkon;
            for (intk_c k = 3; k <= kk - 2; ++k) {
                realk_c dzk = dz[idx1(k)];
                realk_c rdzk = INV_DZ(k);
                realk_c ddzk = ddz[idx1(k)];
                realk_c rddzk = INV_DDZ(k);
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
    // Local upwind correction at right parent boundary.
    if (nrgt == 8) {
#if defined(_OPENMP)
#pragma omp for
#endif
        for (intk_c i = 3; i <= ii - 2; ++i) {
            j = 4;
            realk_c dxi = dx[idx1(i)];
            realk_c ddxi = ddx[idx1(i)];
            realk_c fkdtu = -1.0 * INV_DDY(j) * INV_DX(i) * wkon;
            realk_c fkdtv = -1.0 * INV_DDX(i) * INV_DY(j) * wkon;
            realk_c fkdtw = -1.0 * INV_DDX(i) * INV_DDY(j) * wkon;
            for (intk_c k = 3; k <= kk - 2; ++k) {
                realk_c dzk = dz[idx1(k)];
                realk_c rdzk = INV_DZ(k);
                realk_c ddzk = ddz[idx1(k)];
                realk_c rddzk = INV_DDZ(k);
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
    // Local upwind correction at left parent boundary.
    if (nlft == 8) {
#if defined(_OPENMP)
#pragma omp for
#endif
        for (intk_c i = 3; i <= ii - 2; ++i) {
            j = jj - 2;
            realk_c dxi = dx[idx1(i)];
            realk_c ddxi = ddx[idx1(i)];
            realk_c fkdtu = -1.0 * INV_DDY(j) * INV_DX(i) * wkon;
            realk_c fkdtv = -1.0 * INV_DDX(i) * INV_DY(j) * wkon;
            realk_c fkdtw = -1.0 * INV_DDX(i) * INV_DDY(j) * wkon;
            for (intk_c k = 3; k <= kk - 2; ++k) {
                realk_c dzk = dz[idx1(k)];
                realk_c rdzk = INV_DZ(k);
                realk_c ddzk = ddz[idx1(k)];
                realk_c rddzk = INV_DDZ(k);
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
    // Local upwind correction at bottom parent boundary.
    if (nbot == 8) {
#if defined(_OPENMP)
#pragma omp for
#endif
        for (intk_c i = 3; i <= ii - 2; ++i) {
            realk_c dxi = dx[idx1(i)];
            realk_c ddxi = ddx[idx1(i)];
            for (intk_c j = 3; j <= jj - 2; ++j) {
                k = 4;
                realk_c fkdtu = -1.0 * INV_DDY(j) * INV_DX(i) * wkon;
                realk_c fkdtv = -1.0 * INV_DDX(i) * INV_DY(j) * wkon;
                realk_c fkdtw = -1.0 * INV_DDX(i) * INV_DDY(j) * wkon;
                realk_c dyj = dy[idx1(j)];
                realk_c ddyj = ddy[idx1(j)];
                realk_c auz = dxi * ddyj;
                realk_c avz = ddxi * dyj;
                realk_c rdzk = INV_DZ(k);
                realk_c rddzk = INV_DDZ(k);
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
    // Local upwind correction at top parent boundary.
    if (ntop == 8) {
#if defined(_OPENMP)
#pragma omp for
#endif
        for (intk_c i = 3; i <= ii - 2; ++i) {
            realk_c dxi = dx[idx1(i)];
            realk_c ddxi = ddx[idx1(i)];
            for (intk_c j = 3; j <= jj - 2; ++j) {
                k = kk - 2;
                realk_c fkdtu = -1.0 * INV_DDY(j) * INV_DX(i) * wkon;
                realk_c fkdtv = -1.0 * INV_DDX(i) * INV_DY(j) * wkon;
                realk_c fkdtw = -1.0 * INV_DDX(i) * INV_DDY(j) * wkon;
                realk_c dyj = dy[idx1(j)];
                realk_c ddyj = ddy[idx1(j)];
                realk_c auz = dxi * ddyj;
                realk_c avz = ddxi * dyj;
                realk_c rdzk = INV_DZ(k);
                realk_c rddzk = INV_DDZ(k);
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

    // Conservative parent-boundary reconstruction on back boundary:
    // build wcw/wcv fluxes, distribute to intermediate cells, write to wo/vo.
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
            realk_c fkdtw = -1.0 * INV_DDX(i) * INV_DDY(j) * wkon;
            for (intk_c k = 2; k <= kk - 2; ++k) {
                realk_c rdzk = INV_DZ(k);
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
            realk_c fkdtv = -1.0 * INV_DDX(i) * INV_DY(j) * wkon;
            for (intk_c k = 3; k <= kk - 2; ++k) {
                realk_c rddzk = INV_DDZ(k);
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

    // Conservative parent-boundary reconstruction on front boundary:
    // build wcw/wcv fluxes, distribute to intermediate cells, write to wo/vo.
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
            realk_c fkdtw = -1.0 * INV_DDX(i) * INV_DDY(j) * wkon;
            for (intk_c k = 2; k <= kk - 2; ++k) {
                realk_c rdzk = INV_DZ(k);
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
            realk_c fkdtv = -1.0 * INV_DDX(i) * INV_DY(j) * wkon;
            for (intk_c k = 3; k <= kk - 2; ++k) {
                realk_c rddzk = INV_DDZ(k);
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

    // Conservative parent-boundary reconstruction on top boundary:
    // build wcu/wcv fluxes, distribute to intermediate cells, write to uo/vo.
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
        realk_c rddzk = INV_DDZ(k);
#if defined(_OPENMP)
#pragma omp for collapse(2)
#endif
        for (intk_c i = 2; i <= ii - 2; ++i) {
            for (intk_c j = 3; j <= jj - 2; ++j) {
                realk_c fkdtu = -1.0 * INV_DDY(j) * INV_DX(i) * wkon;
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
        rddzk = INV_DDZ(k);
#if defined(_OPENMP)
#pragma omp for collapse(2)
#endif
        for (intk_c i = 3; i <= ii - 2; ++i) {
            for (intk_c j = 2; j <= jj - 2; ++j) {
                realk_c fkdtv = -1.0 * INV_DDX(i) * INV_DY(j) * wkon;
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

    // Conservative parent-boundary reconstruction on bottom boundary:
    // build wcu/wcv fluxes, distribute to intermediate cells, write to uo/vo.
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
        realk_c rddzk = INV_DDZ(k);
#if defined(_OPENMP)
#pragma omp for collapse(2)
#endif
        for (intk_c i = 2; i <= ii - 2; ++i) {
            for (intk_c j = 3; j <= jj - 2; ++j) {
                realk_c fkdtu = -1.0 * INV_DDY(j) * INV_DX(i) * wkon;
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
        rddzk = INV_DDZ(k);
#if defined(_OPENMP)
#pragma omp for collapse(2)
#endif
        for (intk_c i = 3; i <= ii - 2; ++i) {
            for (intk_c j = 2; j <= jj - 2; ++j) {
                realk_c fkdtv = -1.0 * INV_DDX(i) * INV_DY(j) * wkon;
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

    // Conservative parent-boundary reconstruction on left boundary:
    // build wcu/wcw fluxes, distribute to intermediate cells, write to uo/wo.
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
            realk_c fkdtu = -1.0 * INV_DDY(j) * INV_DX(i) * wkon;
            for (intk_c k = 3; k <= kk - 2; ++k) {
                realk_c rddzk = INV_DDZ(k);
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
            realk_c fkdtw = -1.0 * INV_DDX(i) * INV_DDY(j) * wkon;
            for (intk_c k = 2; k <= kk - 2; ++k) {
                realk_c rdzk = INV_DZ(k);
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

    // Conservative parent-boundary reconstruction on right boundary:
    // build wcu/wcw fluxes, distribute to intermediate cells, write to uo/wo.
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
            realk_c fkdtu = -1.0 * INV_DDY(j) * INV_DX(i) * wkon;
            for (intk_c k = 3; k <= kk - 2; ++k) {
                realk_c rddzk = INV_DDZ(k);
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
            realk_c fkdtw = -1.0 * INV_DDX(i) * INV_DDY(j) * wkon;
            for (intk_c k = 2; k <= kk - 2; ++k) {
                realk_c rdzk = INV_DZ(k);
                wo[idx3(kk, jj, k, j, i)] = wo[idx3(kk, jj, k, j, i)] +
                    fkdtw * rdzk *
                        (wcw[idx3(kk, jj, k, j - 1, i)] -
                            wcw[idx3(kk, jj, k, j, i)]);
            }
        }
    }
    #undef INV_DX
    #undef INV_DY
    #undef INV_DZ
    #undef INV_DDX
    #undef INV_DDY
    #undef INV_DDZ
}

#pragma omp end declare target


void tstle4_par_impl_c(realk_c *uo, realk_c *vo, realk_c *wo,
    const realk_c *u, const realk_c *v, const realk_c *w,
    const realk_c *ut, const realk_c *vt, const realk_c *wt,
    const realk_c *dx, const realk_c *dy, const realk_c *dz,
    const realk_c *ddx, const realk_c *ddy, const realk_c *ddz,
    realk_c *wcu, realk_c *wcv, realk_c *wcw,
    const intk_c nmygrids_in, const tstle4_grid_t *grids)
{
    // Execute one team per local grid and apply parent-boundary correction.
    #pragma omp target teams distribute
    for (intk_c igrid = 0; igrid < nmygrids_in; ++igrid)
    {
        const tstle4_grid_t *const g = &grids[igrid];

        // Grid-local views into the global contiguous field storage.
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

        realk_c *const wcu_g = &wcu[(size_t)(g->ip3 - 1)];
        realk_c *const wcv_g = &wcv[(size_t)(g->ip3 - 1)];
        realk_c *const wcw_g = &wcw[(size_t)(g->ip3 - 1)];

        // Threads collaborate on the correction loops for this grid.
        #pragma omp parallel
        {
            tstle4_par_c(g->kk, g->jj, g->ii,
                uo_g, vo_g, wo_g, u_g, v_g, w_g, ut_g, vt_g, wt_g,
                dx_g, dy_g, dz_g, ddx_g, ddy_g, ddz_g, wcu_g, wcv_g, wcw_g,
                g->nfro, g->nbac, g->nrgt, g->nlft, g->nbot, g->ntop);
        }
    }
}
