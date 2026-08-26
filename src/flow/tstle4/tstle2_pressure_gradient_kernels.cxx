#include "tstle4_kernels_shared.h"

#pragma omp begin declare target

// Pressure-gradient source term with boundary-aware loop extents.
void tstle4_gradp_kernel_c(
    const intk_c kk, const intk_c jj, const intk_c ii,
    realk_c *uo, realk_c *vo, realk_c *wo, const realk_c *p,
    const realk_c *dx, const realk_c *dy, const realk_c *dz,
    const intk_c nfro, const intk_c nbac, const intk_c nrgt,
    const intk_c nlft, const intk_c nbot, const intk_c ntop,
    const realk_c gpx_in, const realk_c gpy_in, const realk_c gpz_in,
    const realk_c inv_rho)
    {

    intk_c nfu = 0; intk_c nbu = 0; intk_c nrv = 0;
    intk_c nlv = 0; intk_c nbw = 0; intk_c ntw = 0;

    if (nbac == 7) nbu = 1;
    if (nlft == 7) nlv = 1;
    if (ntop == 7) ntw = 1;

    if (nfro == 3) nfu = 1;
    if (nbac == 3) nbu = 1;
    if (nrgt == 3) nrv = 1;
    if (nlft == 3) nlv = 1;
    if (nbot == 3) nbw = 1;
    if (ntop == 3) ntw = 1;

    #pragma omp for collapse(3)
    for (intk_c i = 3 - nfu; i <= ii - 3 + nbu; ++i) {
        for (intk_c j = 3; j <= jj - 2; ++j) {
            for (intk_c k = 3; k <= kk - 2; ++k) {
                const realk_c dxi = dx[idx1(i)];
                uo[idx3(kk, jj, k, j, i)] = uo[idx3(kk, jj, k, j, i)] -
                    inv_rho / dxi * (p[idx3(kk, jj, k, j, i + 1)] -
                    p[idx3(kk, jj, k, j, i)] + gpx_in * dxi);
            }
        }
    }

    #pragma omp for collapse(3)
    for (intk_c i = 3; i <= ii - 2; ++i) {
        for (intk_c j = 3 - nrv; j <= jj - 3 + nlv; ++j) {
            for (intk_c k = 3; k <= kk - 2; ++k) {
                const realk_c dyj = dy[idx1(j)];
                vo[idx3(kk, jj, k, j, i)] = vo[idx3(kk, jj, k, j, i)] -
                    inv_rho / dyj * (p[idx3(kk, jj, k, j + 1, i)] -
                    p[idx3(kk, jj, k, j, i)] + gpy_in * dyj);
            }
        }
    }

    #pragma omp for collapse(3)
    for (intk_c i = 3; i <= ii - 2; ++i) {
        for (intk_c j = 3; j <= jj - 2; ++j) {
            for (intk_c k = 3 - nbw; k <= kk - 3 + ntw; ++k) {
                const realk_c dzk = dz[idx1(k)];
                wo[idx3(kk, jj, k, j, i)] = wo[idx3(kk, jj, k, j, i)] -
                    inv_rho / dzk * (p[idx3(kk, jj, k + 1, j, i)] -
                    p[idx3(kk, jj, k, j, i)] + gpz_in * dzk);
            }
        }
    }

}

#pragma omp end declare target


// Full gradp implementation: iterate all grids in C using Fortran-provided metadata.
extern "C" void tstle4_gradp_impl_c(realk_c *uo, realk_c *vo, realk_c *wo,
    const realk_c *p, const realk_c *dx, const realk_c *dy, const realk_c *dz,
    const realk_c rho_in, const intk_c nmygrids_in, const tstle4_grid_t *grids)
{

    const realk_c inv_rho = (realk_c) 1.0 / rho_in;

    #pragma omp target teams distribute
    for (intk_c igrid = 0; igrid < nmygrids_in; ++igrid)
    {
        // Obtaining the memory locations for the current grid
        const tstle4_grid_t *const g = &grids[igrid];
        realk_c *const uo_g = &uo[(size_t)(g->ip3 - 1)];
        realk_c *const vo_g = &vo[(size_t)(g->ip3 - 1)];
        realk_c *const wo_g = &wo[(size_t)(g->ip3 - 1)];
        const realk_c *const p_g = &p[(size_t)(g->ip3 - 1)];
        const realk_c *const dx_g = &dx[(size_t)(g->ipx - 1)];
        const realk_c *const dy_g = &dy[(size_t)(g->ipy - 1)];
        const realk_c *const dz_g = &dz[(size_t)(g->ipz - 1)];

        #pragma omp parallel
        { // Start of parallel region

        // Calling the gradp kernel for each grid with its specific metadata
        tstle4_gradp_kernel_c(g->kk, g->jj, g->ii, uo_g, vo_g, wo_g,
            p_g, dx_g, dy_g, dz_g, g->nfro, g->nbac, g->nrgt, g->nlft,
            g->nbot, g->ntop, g->gpx, g->gpy, g->gpz, inv_rho);

        } // End of parallel region
    }
}
