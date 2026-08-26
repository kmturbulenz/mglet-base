#include <math.h>
#include <stddef.h>
#include <float.h>

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
    intk_c ii;
    intk_c jj;
    intk_c kk;
    intk_c ip3;
    intk_c ipx;
    intk_c ipy;
    intk_c ipz;
    intk_c nfro;
    intk_c nbac;
    intk_c nrgt;
    intk_c nlft;
    intk_c nbot;
    intk_c ntop;
} lesmodel_grid_t;

static inline size_t idx3(const intk_c kk, const intk_c jj, const intk_c k,
    const intk_c j, const intk_c i)
{
    return (size_t)(k - 1) +
        (size_t)kk * ((size_t)(j - 1) + (size_t)jj * (size_t)(i - 1));
}

static inline size_t idx1(const intk_c i)
{
    return (size_t)(i - 1);
}

static inline realk_c abs_c(const realk_c x)
{
#ifdef _MGLET_DOUBLE_PRECISION_
    return fabs(x);
#else
    return fabsf(x);
#endif
}

static inline realk_c sqrt_c(const realk_c x)
{
#ifdef _MGLET_DOUBLE_PRECISION_
    return sqrt(x);
#else
    return sqrtf(x);
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

static inline realk_c acos_c(const realk_c x)
{
#ifdef _MGLET_DOUBLE_PRECISION_
    return acos(x);
#else
    return acosf(x);
#endif
}

static inline realk_c cos_c(const realk_c x)
{
#ifdef _MGLET_DOUBLE_PRECISION_
    return cos(x);
#else
    return cosf(x);
#endif
}

static inline realk_c max2(const realk_c a, const realk_c b)
{
    return (a > b) ? a : b;
}

static inline realk_c divide0_c(const realk_c a, const realk_c b)
{
    return (b == (realk_c)0.0) ? (realk_c)0.0 : a / b;
}

static inline realk_c cube_root_c(const realk_c x)
{
    return pow_c(x, (realk_c)(1.0 / 3.0));
}

static inline void sort3_desc(realk_c *a, realk_c *b, realk_c *c)
{
    realk_c t;
    if (*a < *b) {
        t = *a;
        *a = *b;
        *b = t;
    }
    if (*b < *c) {
        t = *b;
        *b = *c;
        *c = t;
    }
    if (*a < *b) {
        t = *a;
        *a = *b;
        *b = t;
    }
}

#if defined(_OPENMP)
#pragma omp begin declare target
#endif

static inline realk_c gradp2_c(
    const realk_c uquer, const realk_c dds, const realk_c gmol,
    const realk_c rho)
{
    const realk_c cwa = (realk_c)8.3;
    const realk_c cwb = (realk_c)(1.0 / 7.0);
    const realk_c cpo1 = (realk_c)1.0 - cwb;
    const realk_c cpo2 = (realk_c)1.0 + cwb;
    const realk_c cpo6 =
        (realk_c)0.5 * gmol / rho * pow_c(cwa, (realk_c)2.0 / cpo1);
    const realk_c cpo11 =
        (realk_c)2.0 * cpo6 * cwb * cpo1 / pow_c((realk_c)2.0, cwb);
    const realk_c cpo12 =
        cwb * cpo2 / pow_c((realk_c)2.0, cwb - (realk_c)1.0);

    const realk_c uquern = abs_c(uquer);
    const realk_c vz = (uquer < (realk_c)0.0) ? (realk_c)-1.0 : (realk_c)1.0;
    const realk_c rdds = (realk_c)1.0 / dds;

    if (uquern >= cpo6 * rdds) {
        return rdds * (cpo12 * uquern + cpo11 * rdds) * vz;
    }

    return (realk_c)2.0 * uquern / dds * vz;
}

static inline realk_c sabs_c(const realk_c dudx, const realk_c dudy,
    const realk_c dudz, const realk_c dvdx, const realk_c dvdy,
    const realk_c dvdz, const realk_c dwdx, const realk_c dwdy,
    const realk_c dwdz)
{
    return sqrt_c(dudx * dudx + dvdy * dvdy + dwdz * dwdz +
        (realk_c)0.5 * (dvdx + dudy) * (dvdx + dudy) +
        (realk_c)0.5 * (dwdx + dudz) * (dwdx + dudz) +
        (realk_c)0.5 * (dwdy + dvdz) * (dwdy + dvdz));
}

static inline realk_c smagorinsky_c(const realk_c cm, const realk_c dudx,
    const realk_c dudy, const realk_c dudz, const realk_c dvdx,
    const realk_c dvdy, const realk_c dvdz, const realk_c dwdx,
    const realk_c dwdy, const realk_c dwdz)
{
    const realk_c root2 = sqrt_c((realk_c)2.0);
    const realk_c sabs =
        sabs_c(dudx, dudy, dudz, dvdx, dvdy, dvdz, dwdx, dwdy, dwdz);
    return cm * cm * root2 * sabs;
}

static inline realk_c wale_c(const realk_c cm, const realk_c dudx,
    const realk_c dudy, const realk_c dudz, const realk_c dvdx,
    const realk_c dvdy, const realk_c dvdz, const realk_c dwdx,
    const realk_c dwdy, const realk_c dwdz)
{
    const realk_c tr = dudx * dudx + dudy * dvdx + dudz * dwdx + dvdx * dudy +
        dvdy * dvdy + dvdz * dwdy + dwdx * dudz + dwdy * dvdz + dwdz * dwdz;
    const realk_c sd00 =
        (realk_c)0.5 * (dudx * dudx + dudy * dvdx + dudz * dwdx + dudx * dudx +
                           dudy * dvdx + dudz * dwdx) -
        (realk_c)(1.0 / 3.0) * tr;
    const realk_c sd01 =
        (realk_c)0.5 * (dudx * dudy + dudy * dvdy + dudz * dwdy + dvdx * dudx +
                           dvdy * dvdx + dvdz * dwdx);
    const realk_c sd02 =
        (realk_c)0.5 * (dudx * dudz + dudy * dvdz + dudz * dwdz + dwdx * dudx +
                           dwdy * dvdx + dwdz * dwdx);
    const realk_c sd11 =
        (realk_c)0.5 * (dvdx * dudy + dvdy * dvdy + dvdz * dwdy + dvdx * dudy +
                           dvdy * dvdy + dvdz * dwdy) -
        (realk_c)(1.0 / 3.0) * tr;
    const realk_c sd12 =
        (realk_c)0.5 * (dvdx * dudz + dvdy * dvdz + dvdz * dwdz + dwdx * dudy +
                           dwdy * dvdy + dwdz * dwdy);
    const realk_c sd22 =
        (realk_c)0.5 * (dwdx * dudz + dwdy * dvdz + dwdz * dwdz + dwdx * dudz +
                           dwdy * dvdz + dwdz * dwdz) -
        (realk_c)(1.0 / 3.0) * tr;

    const realk_c sd_abs =
        sqrt_c(sd00 * sd00 + (realk_c)2.0 * sd01 * sd01 +
            (realk_c)2.0 * sd02 * sd02 + sd11 * sd11 +
            (realk_c)2.0 * sd12 * sd12 + sd22 * sd22);
    const realk_c s_abs =
        sabs_c(dudx, dudy, dudz, dvdx, dvdy, dvdz, dwdx, dwdy, dwdz);

#ifdef _MGLET_DOUBLE_PRECISION_
    const realk_c huge_val = (realk_c)DBL_MAX;
#else
    const realk_c huge_val = (realk_c)FLT_MAX;
#endif
    const realk_c sabs_tol =
        (realk_c)1.0 / pow_c(huge_val, (realk_c)(1.0 / 6.0));
    const realk_c sd_abs_root = sqrt_c(sd_abs);

    if (sd_abs_root < sabs_tol * s_abs) {
        return (realk_c)0.0;
    }

    return cm * cm * sd_abs_root /
        (pow_c(divide0_c(s_abs, sd_abs_root), (realk_c)5.0) + (realk_c)1.0);
}

static inline void sigma_eigs_c(const realk_c dudx, const realk_c dudy,
    const realk_c dudz, const realk_c dvdx, const realk_c dvdy,
    const realk_c dvdz, const realk_c dwdx, const realk_c dwdy,
    const realk_c dwdz, realk_c *eig1, realk_c *eig2, realk_c *eig3)
{
    const realk_c a11 = dudx * dudx + dvdx * dvdx + dwdx * dwdx;
    const realk_c a22 = dudy * dudy + dvdy * dvdy + dwdy * dwdy;
    const realk_c a33 = dudz * dudz + dvdz * dvdz + dwdz * dwdz;
    const realk_c a12 = dudx * dudy + dvdx * dvdy + dwdx * dwdy;
    const realk_c a13 = dudx * dudz + dvdx * dvdz + dwdx * dwdz;
    const realk_c a23 = dudy * dudz + dvdy * dvdz + dwdy * dwdz;

    const realk_c p1 = a12 * a12 + a13 * a13 + a23 * a23;

    if (p1 == (realk_c)0.0) {
        *eig1 = a11;
        *eig2 = a22;
        *eig3 = a33;
        sort3_desc(eig1, eig2, eig3);
        return;
    }

    const realk_c q = (a11 + a22 + a33) / (realk_c)3.0;
    const realk_c b11 = a11 - q;
    const realk_c b22 = a22 - q;
    const realk_c b33 = a33 - q;
    const realk_c p2 = b11 * b11 + b22 * b22 + b33 * b33 + (realk_c)2.0 * p1;
    const realk_c p = sqrt_c(p2 / (realk_c)6.0);

    const realk_c inv_p = (realk_c)1.0 / p;
    const realk_c c11 = b11 * inv_p;
    const realk_c c22 = b22 * inv_p;
    const realk_c c33 = b33 * inv_p;
    const realk_c c12 = a12 * inv_p;
    const realk_c c13 = a13 * inv_p;
    const realk_c c23 = a23 * inv_p;

    const realk_c detc = c11 * (c22 * c33 - c23 * c23) -
        c12 * (c12 * c33 - c23 * c13) + c13 * (c12 * c23 - c22 * c13);
    realk_c r = detc * (realk_c)0.5;

    if (r <= (realk_c)-1.0) {
        r = (realk_c)-1.0;
    } else if (r >= (realk_c)1.0) {
        r = (realk_c)1.0;
    }

    const realk_c pi = (realk_c)3.1415926535897932384626433832795;
    const realk_c phi = acos_c(r) / (realk_c)3.0;

    *eig1 = q + (realk_c)2.0 * p * cos_c(phi);
    *eig3 = q + (realk_c)2.0 * p * cos_c(phi + (realk_c)(2.0 * pi / 3.0));
    *eig2 = (realk_c)3.0 * q - *eig1 - *eig3;

    sort3_desc(eig1, eig2, eig3);
}

static inline realk_c sigma_c(const realk_c cm, const realk_c dudx,
    const realk_c dudy, const realk_c dudz, const realk_c dvdx,
    const realk_c dvdy, const realk_c dvdz, const realk_c dwdx,
    const realk_c dwdy, const realk_c dwdz)
{
    realk_c eig1, eig2, eig3;
    sigma_eigs_c(
        dudx, dudy, dudz, dvdx, dvdy, dvdz, dwdx, dwdy, dwdz, &eig1, &eig2,
        &eig3);

    eig1 = max2(eig1, (realk_c)0.0);
    eig2 = max2(eig2, (realk_c)0.0);
    eig3 = max2(eig3, (realk_c)0.0);

    const realk_c sigma1 = sqrt_c(eig1);
    const realk_c sigma2 = sqrt_c(eig2);
    const realk_c sigma3 = sqrt_c(eig3);

    return cm * cm * divide0_c(
        sigma3 * (sigma1 - sigma2) * (sigma2 - sigma3), sigma1 * sigma1);
}

void lesmodel_efvisc_gc_c(const intk_c kk, const intk_c jj, const intk_c ii,
    const intk_c nfro, const intk_c nbac, const intk_c nrgt,
    const intk_c nlft, const intk_c nbot, const intk_c ntop,
    const realk_c *dx, const realk_c *dy, const realk_c *dz,
    const realk_c *ddx, const realk_c *ddy, const realk_c *ddz,
    const realk_c *rddx, const realk_c *rddy, const realk_c *rddz,
    const realk_c *u, const realk_c *v, const realk_c *w, const realk_c *bp,
    realk_c *g, const intk_c ilesmodel_in, const realk_c cm,
    const realk_c gmol, const realk_c rho)
{
    #pragma omp for collapse(3)
    for (intk_c i = 3; i <= ii - 2; ++i) {
        for (intk_c j = 3; j <= jj - 2; ++j) {
            for (intk_c k = 3; k <= kk - 2; ++k) {
                const realk_c dxf = (realk_c)0.25 * dx[idx1(i - 1)] * rddx[idx1(i)];
                const realk_c dudx =
                    (u[idx3(kk, jj, k, j, i)] - u[idx3(kk, jj, k, j, i - 1)]) *
                    rddx[idx1(i)];
                realk_c dudy =
                    ((u[idx3(kk, jj, k, j + 1, i)] -
                         u[idx3(kk, jj, k, j - 1, i)]) *
                        dxf +
                        (u[idx3(kk, jj, k, j + 1, i - 1)] -
                            u[idx3(kk, jj, k, j - 1, i - 1)]) *
                            ((realk_c)0.5 - dxf)) *
                    rddy[idx1(j)];
                realk_c dudz =
                    ((u[idx3(kk, jj, k + 1, j, i)] -
                         u[idx3(kk, jj, k - 1, j, i)]) *
                        dxf +
                        (u[idx3(kk, jj, k + 1, j, i - 1)] -
                            u[idx3(kk, jj, k - 1, j, i - 1)]) *
                            ((realk_c)0.5 - dxf)) *
                    rddz[idx1(k)];

                const realk_c dyf = (realk_c)0.25 * dy[idx1(j - 1)] * rddy[idx1(j)];
                realk_c dvdx =
                    ((v[idx3(kk, jj, k, j, i + 1)] -
                         v[idx3(kk, jj, k, j, i - 1)]) *
                        dyf +
                        (v[idx3(kk, jj, k, j - 1, i + 1)] -
                            v[idx3(kk, jj, k, j - 1, i - 1)]) *
                            ((realk_c)0.5 - dyf)) *
                    rddx[idx1(i)];
                const realk_c dvdy =
                    (v[idx3(kk, jj, k, j, i)] - v[idx3(kk, jj, k, j - 1, i)]) *
                    rddy[idx1(j)];
                realk_c dvdz =
                    ((v[idx3(kk, jj, k + 1, j, i)] -
                         v[idx3(kk, jj, k - 1, j, i)]) *
                        dyf +
                        (v[idx3(kk, jj, k + 1, j - 1, i)] -
                            v[idx3(kk, jj, k - 1, j - 1, i)]) *
                            ((realk_c)0.5 - dyf)) *
                    rddz[idx1(k)];

                const realk_c dzf = (realk_c)0.25 * dz[idx1(k - 1)] * rddz[idx1(k)];
                realk_c dwdx =
                    ((w[idx3(kk, jj, k, j, i + 1)] -
                         w[idx3(kk, jj, k, j, i - 1)]) *
                        dzf +
                        (w[idx3(kk, jj, k - 1, j, i + 1)] -
                            w[idx3(kk, jj, k - 1, j, i - 1)]) *
                            ((realk_c)0.5 - dzf)) *
                    rddx[idx1(i)];
                realk_c dwdy =
                    ((w[idx3(kk, jj, k, j + 1, i)] -
                         w[idx3(kk, jj, k, j - 1, i)]) *
                        dzf +
                        (w[idx3(kk, jj, k - 1, j + 1, i)] -
                            w[idx3(kk, jj, k - 1, j - 1, i)]) *
                            ((realk_c)0.5 - dzf)) *
                    rddy[idx1(j)];
                const realk_c dwdz =
                    (w[idx3(kk, jj, k, j, i)] - w[idx3(kk, jj, k - 1, j, i)]) *
                    rddz[idx1(k)];

                const realk_c delta =
                    cube_root_c(ddx[idx1(i)] * ddy[idx1(j)] * ddz[idx1(k)]) *
                    bp[idx3(kk, jj, k, j, i)];

                if ((nfro == 5) && (i == 3)) {
                    const realk_c dyf2 =
                        (realk_c)0.5 * dy[idx1(j - 1)] * rddy[idx1(j)];
                    dvdx = gradp2_c(
                        dyf2 * v[idx3(kk, jj, k, j, i)] +
                            ((realk_c)1.0 - dyf2) * v[idx3(kk, jj, k, j - 1, i)],
                        ddx[idx1(i)], gmol, rho);

                    const realk_c dzf2 =
                        (realk_c)0.5 * dz[idx1(k - 1)] * rddz[idx1(k)];
                    dwdx = gradp2_c(
                        dzf2 * w[idx3(kk, jj, k, j, i)] +
                            ((realk_c)1.0 - dzf2) * w[idx3(kk, jj, k - 1, j, i)],
                        ddx[idx1(i)], gmol, rho);
                }
                if ((nbac == 5) && (i == ii - 2)) {
                    const realk_c dyf2 =
                        (realk_c)0.5 * dy[idx1(j - 1)] * rddy[idx1(j)];
                    dvdx = -gradp2_c(
                        dyf2 * v[idx3(kk, jj, k, j, i)] +
                            ((realk_c)1.0 - dyf2) * v[idx3(kk, jj, k, j - 1, i)],
                        ddx[idx1(i)], gmol, rho);

                    const realk_c dzf2 =
                        (realk_c)0.5 * dz[idx1(k - 1)] * rddz[idx1(k)];
                    dwdx = -gradp2_c(
                        dzf2 * w[idx3(kk, jj, k, j, i)] +
                            ((realk_c)1.0 - dzf2) * w[idx3(kk, jj, k - 1, j, i)],
                        ddx[idx1(i)], gmol, rho);
                }
                if ((nrgt == 5) && (j == 3)) {
                    const realk_c dxf2 =
                        (realk_c)0.5 * dx[idx1(i - 1)] * rddx[idx1(i)];
                    dudy = gradp2_c(
                        dxf2 * u[idx3(kk, jj, k, j, i)] +
                            ((realk_c)1.0 - dxf2) * u[idx3(kk, jj, k, j, i - 1)],
                        ddy[idx1(j)], gmol, rho);

                    const realk_c dzf2 =
                        (realk_c)0.5 * dz[idx1(k - 1)] * rddz[idx1(k)];
                    dwdy = gradp2_c(
                        dzf2 * w[idx3(kk, jj, k, j, i)] +
                            ((realk_c)1.0 - dzf2) * w[idx3(kk, jj, k - 1, j, i)],
                        ddy[idx1(j)], gmol, rho);
                }
                if ((nlft == 5) && (j == jj - 2)) {
                    const realk_c dxf2 =
                        (realk_c)0.5 * dx[idx1(i - 1)] * rddx[idx1(i)];
                    dudy = -gradp2_c(
                        dxf2 * u[idx3(kk, jj, k, j, i)] +
                            ((realk_c)1.0 - dxf2) * u[idx3(kk, jj, k, j, i - 1)],
                        ddy[idx1(j)], gmol, rho);

                    const realk_c dzf2 =
                        (realk_c)0.5 * dz[idx1(k - 1)] * rddz[idx1(k)];
                    dwdy = -gradp2_c(
                        dzf2 * w[idx3(kk, jj, k, j, i)] +
                            ((realk_c)1.0 - dzf2) * w[idx3(kk, jj, k - 1, j, i)],
                        ddy[idx1(j)], gmol, rho);
                }
                if ((nbot == 5) && (k == 3)) {
                    const realk_c dxf2 =
                        (realk_c)0.5 * dx[idx1(i - 1)] * rddx[idx1(i)];
                    dudz = gradp2_c(
                        dxf2 * u[idx3(kk, jj, k, j, i)] +
                            ((realk_c)1.0 - dxf2) * u[idx3(kk, jj, k, j, i - 1)],
                        ddz[idx1(k)], gmol, rho);

                    const realk_c dyf2 =
                        (realk_c)0.5 * dy[idx1(j - 1)] * rddy[idx1(j)];
                    dvdz = gradp2_c(
                        dyf2 * v[idx3(kk, jj, k, j, i)] +
                            ((realk_c)1.0 - dyf2) * v[idx3(kk, jj, k, j - 1, i)],
                        ddz[idx1(k)], gmol, rho);
                }
                if ((ntop == 5) && (k == kk - 2)) {
                    const realk_c dxf2 =
                        (realk_c)0.5 * dx[idx1(i - 1)] * rddx[idx1(i)];
                    dudz = -gradp2_c(
                        dxf2 * u[idx3(kk, jj, k, j, i)] +
                            ((realk_c)1.0 - dxf2) * u[idx3(kk, jj, k, j, i - 1)],
                        ddz[idx1(k)], gmol, rho);

                    const realk_c dyf2 =
                        (realk_c)0.5 * dy[idx1(j - 1)] * rddy[idx1(j)];
                    dvdz = -gradp2_c(
                        dyf2 * v[idx3(kk, jj, k, j, i)] +
                            ((realk_c)1.0 - dyf2) * v[idx3(kk, jj, k, j - 1, i)],
                        ddz[idx1(k)], gmol, rho);
                }

                realk_c dm = (realk_c)0.0;
                switch (ilesmodel_in) {
                case 1:
                    dm = smagorinsky_c(
                        cm, dudx, dudy, dudz, dvdx, dvdy, dvdz, dwdx, dwdy, dwdz);
                    break;
                case 2:
                    dm = wale_c(
                        cm, dudx, dudy, dudz, dvdx, dvdy, dvdz, dwdx, dwdy, dwdz);
                    break;
                case 5:
                    dm = sigma_c(
                        cm, dudx, dudy, dudz, dvdx, dvdy, dvdz, dwdx, dwdy, dwdz);
                    break;
                default:
                    dm = (realk_c)0.0;
                    break;
                }

                g[idx3(kk, jj, k, j, i)] = rho * delta * delta * dm + gmol;
            }
        }
    }
}

#if defined(_OPENMP)
#pragma omp end declare target
#endif

void lesmodel_gc_impl_c(realk_c *g, const realk_c *u, const realk_c *v,
    const realk_c *w, const realk_c *bp, const realk_c *dx,
    const realk_c *dy, const realk_c *dz, const realk_c *ddx,
    const realk_c *ddy, const realk_c *ddz, const realk_c *rddx,
    const realk_c *rddy, const realk_c *rddz, const intk_c ilesmodel_in,
    const realk_c cm_in, const realk_c gmol_in, const realk_c rho_in,
    const intk_c nmygrids_in, const lesmodel_grid_t *grids)
{

    #pragma omp target teams distribute
    for (intk_c igrid = 0; igrid < nmygrids_in; ++igrid) {
        const lesmodel_grid_t *const ginfo = &grids[igrid];

        realk_c *const g_g = &g[(size_t)(ginfo->ip3 - 1)];
        const realk_c *const u_g = &u[(size_t)(ginfo->ip3 - 1)];
        const realk_c *const v_g = &v[(size_t)(ginfo->ip3 - 1)];
        const realk_c *const w_g = &w[(size_t)(ginfo->ip3 - 1)];
        const realk_c *const bp_g = &bp[(size_t)(ginfo->ip3 - 1)];
        const realk_c *const dx_g = &dx[(size_t)(ginfo->ipx - 1)];
        const realk_c *const dy_g = &dy[(size_t)(ginfo->ipy - 1)];
        const realk_c *const dz_g = &dz[(size_t)(ginfo->ipz - 1)];
        const realk_c *const ddx_g = &ddx[(size_t)(ginfo->ipx - 1)];
        const realk_c *const ddy_g = &ddy[(size_t)(ginfo->ipy - 1)];
        const realk_c *const ddz_g = &ddz[(size_t)(ginfo->ipz - 1)];
        const realk_c *const rddx_g = &rddx[(size_t)(ginfo->ipx - 1)];
        const realk_c *const rddy_g = &rddy[(size_t)(ginfo->ipy - 1)];
        const realk_c *const rddz_g = &rddz[(size_t)(ginfo->ipz - 1)];

        #pragma omp parallel
        {
            lesmodel_efvisc_gc_c(ginfo->kk, ginfo->jj, ginfo->ii, ginfo->nfro,
                ginfo->nbac, ginfo->nrgt, ginfo->nlft, ginfo->nbot,
                ginfo->ntop, dx_g, dy_g, dz_g, ddx_g, ddy_g, ddz_g, rddx_g,
                rddy_g, rddz_g, u_g, v_g, w_g, bp_g, g_g, ilesmodel_in, cm_in,
                gmol_in, rho_in);
        }
    }
}
