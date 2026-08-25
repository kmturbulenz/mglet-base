#ifndef TSTLE4_KERNELS_SHARED_H
#define TSTLE4_KERNELS_SHARED_H

#include <math.h>
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


static inline size_t idx3(const intk_c kk, const intk_c jj, const intk_c k,
    const intk_c j, const intk_c i)
{
    return (size_t) (k - 1) + (size_t)kk * ((size_t)(j - 1) + (size_t)jj * (size_t)(i - 1));
}

static inline size_t idx1(const intk_c i)
{
    return (size_t)(i - 1);
}

static inline realk_c max2(const realk_c a, const realk_c b)
{
    return (a > b) ? a : b;
}

static inline realk_c abs_c(const realk_c x)
{
#ifdef _MGLET_DOUBLE_PRECISION_
    return fabs(x);
#else
    return fabsf(x);
#endif
}

static inline realk_c pow_c(const realk_c x, const realk_c p)
{
#ifdef _MGLET_DOUBLE_PRECISION_
    return pow(x, p);
#else
    return powf(x, p);
#endif
}


// TO BE MOVED
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

#endif
