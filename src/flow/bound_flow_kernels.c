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

typedef struct {
    intk_c kk;
    intk_c jj;
    intk_c ii;
    intk_c iface;
    intk_c ityp;
    intk_c ip3;
    intk_c ipx;
    intk_c ipy;
    intk_c ipz;
    intk_c ipbb;
    realk_c pinf;
} bound_flow_task_t;

/* Flatten (k,j,i) Fortran indexing into contiguous C storage. */
static inline size_t idx3(const intk_c kk, const intk_c jj, const intk_c k,
                          const intk_c j, const intk_c i) {
    return (size_t)(k - 1) +
           (size_t)kk * ((size_t)(j - 1) + (size_t)jj * (size_t)(i - 1));
}

static inline size_t idx2(const intk_c n1, const intk_c a, const intk_c b) {
    return (size_t)(a - 1) + (size_t)n1 * (size_t)(b - 1);
}

static inline realk_c abs_c(const realk_c x) {
#ifdef _MGLET_DOUBLE_PRECISION_
    return fabs(x);
#else
    return fabsf(x);
#endif
}

static inline realk_c min2(const realk_c a, const realk_c b) {
    return (a < b) ? a : b;
}

static inline realk_c sign_c(const realk_c a, const realk_c b) {
    return (b < (realk_c)0.0) ? -abs_c(a) : abs_c(a);
}

#if defined(_OPENMP)
#pragma omp begin declare target
#endif

/* No-buffer boundary path used for slip/no-slip variants. */
void bound_nobuffer_device_c(const intk_c kk, const intk_c jj, const intk_c ii,
                             const intk_c iface, const intk_c ityp, realk_c *u,
                             realk_c *v, realk_c *w, realk_c *p) {
    intk_c i2, i3, j2, j3, k2, k3, istag2;

    switch (iface) {
    case 1:
    case 2:
        if (iface == 1) {
            i2 = 2;
            i3 = 3;
            istag2 = 2;
        } else {
            i2 = ii - 1;
            i3 = ii - 2;
            istag2 = ii - 2;
        }
#if defined(_OPENMP)
#pragma omp for collapse(2)
#endif
        for (intk_c j = 1; j <= jj; ++j) {
            for (intk_c k = 1; k <= kk; ++k) {
                u[idx3(kk, jj, k, j, istag2)] = (realk_c)0.0;
                if (ityp == 3) {
                    v[idx3(kk, jj, k, j, i2)] = (realk_c)0.0;
                    w[idx3(kk, jj, k, j, i2)] = (realk_c)0.0;
                } else {
                    v[idx3(kk, jj, k, j, i2)] = v[idx3(kk, jj, k, j, i3)];
                    w[idx3(kk, jj, k, j, i2)] = w[idx3(kk, jj, k, j, i3)];
                }
                p[idx3(kk, jj, k, j, i2)] = p[idx3(kk, jj, k, j, i3)];
            }
        }
        break;

    case 3:
    case 4:
        if (iface == 3) {
            j2 = 2;
            j3 = 3;
            istag2 = 2;
        } else {
            j2 = jj - 1;
            j3 = jj - 2;
            istag2 = jj - 2;
        }
#if defined(_OPENMP)
#pragma omp for collapse(2)
#endif
        for (intk_c i = 1; i <= ii; ++i) {
            for (intk_c k = 1; k <= kk; ++k) {
                v[idx3(kk, jj, k, istag2, i)] = (realk_c)0.0;
                if (ityp == 3) {
                    u[idx3(kk, jj, k, j2, i)] = (realk_c)0.0;
                    w[idx3(kk, jj, k, j2, i)] = (realk_c)0.0;
                } else {
                    u[idx3(kk, jj, k, j2, i)] = u[idx3(kk, jj, k, j3, i)];
                    w[idx3(kk, jj, k, j2, i)] = w[idx3(kk, jj, k, j3, i)];
                }
                p[idx3(kk, jj, k, j2, i)] = p[idx3(kk, jj, k, j3, i)];
            }
        }
        break;

    case 5:
    case 6:
        if (iface == 5) {
            k2 = 2;
            k3 = 3;
            istag2 = 2;
        } else {
            k2 = kk - 1;
            k3 = kk - 2;
            istag2 = kk - 2;
        }
#if defined(_OPENMP)
#pragma omp for collapse(2)
#endif
        for (intk_c i = 1; i <= ii; ++i) {
            for (intk_c j = 1; j <= jj; ++j) {
                w[idx3(kk, jj, istag2, j, i)] = (realk_c)0.0;
                if (ityp == 3) {
                    u[idx3(kk, jj, k2, j, i)] = (realk_c)0.0;
                    v[idx3(kk, jj, k2, j, i)] = (realk_c)0.0;
                } else {
                    u[idx3(kk, jj, k2, j, i)] = u[idx3(kk, jj, k3, j, i)];
                    v[idx3(kk, jj, k2, j, i)] = v[idx3(kk, jj, k3, j, i)];
                }
                p[idx3(kk, jj, k2, j, i)] = p[idx3(kk, jj, k3, j, i)];
            }
        }
        break;

    default:
        break;
    }
}

/* X-normal boundary kernel (front/back faces). */
void bfront_device_c(const intk_c kk, const intk_c jj, const intk_c ii,
                     const intk_c i2, const intk_c i3, const intk_c i4,
                     const intk_c istag1, const intk_c istag2, const intk_c dir,
                     const intk_c ityp, const realk_c pinf, realk_c *u,
                     realk_c *v, realk_c *w, realk_c *p, const realk_c *bp,
                     const realk_c *ubuf, const realk_c *vbuf,
                     const realk_c *wbuf, const realk_c *ddy,
                     const realk_c *ddz, const realk_c rho) {
    switch (ityp) {
    case 1:
#if defined(_OPENMP)
#pragma omp for collapse(2)
#endif
        for (intk_c j = 1; j <= jj; ++j) {
            for (intk_c k = 1; k <= kk; ++k) {
                u[idx3(kk, jj, k, j, istag2)] = ubuf[idx2(kk, k, j)];
                v[idx3(kk, jj, k, j, i2)] =
                    (realk_c)2.0 * vbuf[idx2(kk, k, j)] -
                    v[idx3(kk, jj, k, j, i3)];
                w[idx3(kk, jj, k, j, i2)] =
                    (realk_c)2.0 * wbuf[idx2(kk, k, j)] -
                    w[idx3(kk, jj, k, j, i3)];
            }
        }
        break;

    case 2:
#if defined(_OPENMP)
#pragma omp for collapse(2)
#endif
        for (intk_c j = 1; j <= jj - 1; ++j) {
            for (intk_c k = 1; k <= kk - 1; ++k) {
                realk_c flag;

                u[idx3(kk, jj, k, j, istag1)] = u[idx3(kk, jj, k, j, istag2)];
                flag =
                    sign_c((realk_c)1.0,
                           (realk_c)dir * (u[idx3(kk, jj, k, j, istag2)] +
                                           u[idx3(kk, jj, k, j + 1, istag2)]));
                flag = (realk_c)0.5 * (flag + (realk_c)1.0);
                v[idx3(kk, jj, k, j, i2)] =
                    flag * v[idx3(kk, jj, k, j, i3)] +
                    ((realk_c)1.0 - flag) *
                        ((realk_c)2.0 * vbuf[idx2(kk, k, j)] -
                         v[idx3(kk, jj, k, j, i3)]);

                flag =
                    sign_c((realk_c)1.0,
                           (realk_c)dir * (u[idx3(kk, jj, k, j, istag2)] +
                                           u[idx3(kk, jj, k + 1, j, istag2)]));
                flag = (realk_c)0.5 * (flag + (realk_c)1.0);
                w[idx3(kk, jj, k, j, i2)] =
                    flag * w[idx3(kk, jj, k, j, i3)] +
                    ((realk_c)1.0 - flag) *
                        ((realk_c)2.0 * wbuf[idx2(kk, k, j)] -
                         w[idx3(kk, jj, k, j, i3)]);
            }
        }
        break;

    case 3:
#if defined(_OPENMP)
#pragma omp for collapse(2)
#endif
        for (intk_c j = 1; j <= jj; ++j) {
            for (intk_c k = 1; k <= kk; ++k) {
                u[idx3(kk, jj, k, j, istag2)] = (realk_c)0.0;
                v[idx3(kk, jj, k, j, i2)] = (realk_c)0.0;
                w[idx3(kk, jj, k, j, i2)] = (realk_c)0.0;
            }
        }
        break;

    case 4:
#if defined(_OPENMP)
#pragma omp for collapse(2)
#endif
        for (intk_c j = 1; j <= jj; ++j) {
            for (intk_c k = 1; k <= kk; ++k) {
                u[idx3(kk, jj, k, j, istag2)] = (realk_c)0.0;
                v[idx3(kk, jj, k, j, i2)] = v[idx3(kk, jj, k, j, i3)];
                w[idx3(kk, jj, k, j, i2)] = w[idx3(kk, jj, k, j, i3)];
            }
        }
        break;

    case 5:
#if defined(_OPENMP)
#pragma omp for collapse(2)
#endif
        for (intk_c j = 2; j <= jj - 1; ++j) {
            for (intk_c k = 2; k <= kk - 1; ++k) {
                const realk_c sbv =
                    bp[idx3(kk, jj, k, j, i2)] * bp[idx3(kk, jj, k, j + 1, i2)];
                const realk_c sbw =
                    bp[idx3(kk, jj, k, j, i2)] * bp[idx3(kk, jj, k + 1, j, i2)];
                v[idx3(kk, jj, k, j, i2)] =
                    (vbuf[idx2(kk, k, j)] -
                     (realk_c)0.5 * (v[idx3(kk, jj, k, j, i3)] -
                                     v[idx3(kk, jj, k, j, i4)])) *
                        sbv +
                    ((realk_c)1.0 - sbv) * v[idx3(kk, jj, k, j, i2)];
                w[idx3(kk, jj, k, j, i2)] =
                    (wbuf[idx2(kk, k, j)] -
                     (realk_c)0.5 * (w[idx3(kk, jj, k, j, i3)] -
                                     w[idx3(kk, jj, k, j, i4)])) *
                        sbw +
                    ((realk_c)1.0 - sbw) * w[idx3(kk, jj, k, j, i2)];
            }
        }

#if defined(_OPENMP)
#pragma omp for collapse(2)
#endif
        for (intk_c j = 3; j <= jj - 2; j += 2) {
            for (intk_c k = 3; k <= kk - 2; k += 2) {
                const realk_c sbu = bp[idx3(kk, jj, k, j, i2)];
                realk_c sb11;
                realk_c sb12;
                realk_c sb13;
                realk_c sb14;
                realk_c fak;
                realk_c ubfine;

                if (sbu > (realk_c)0.5) {
                    sb11 =
                        bp[idx3(kk, jj, k, j, i2)] * bp[idx3(kk, jj, k, j, i3)];
                    sb12 = bp[idx3(kk, jj, k, j + 1, i2)] *
                           bp[idx3(kk, jj, k, j + 1, i3)];
                    sb13 = bp[idx3(kk, jj, k + 1, j, i2)] *
                           bp[idx3(kk, jj, k + 1, j, i3)];
                    sb14 = bp[idx3(kk, jj, k + 1, j + 1, i2)] *
                           bp[idx3(kk, jj, k + 1, j + 1, i3)];
                } else {
                    sb11 = min2(bp[idx3(kk, jj, k, j, i3)] +
                                    bp[idx3(kk, jj, k, j, i2)] +
                                    bp[idx3(kk, jj, k, j, i4)] +
                                    bp[idx3(kk, jj, k, j - 1, i3)] +
                                    bp[idx3(kk, jj, k, j + 1, i3)] +
                                    bp[idx3(kk, jj, k - 1, j, i3)] +
                                    bp[idx3(kk, jj, k + 1, j, i3)],
                                (realk_c)1.0);
                    sb12 = min2(bp[idx3(kk, jj, k, j + 1, i3)] +
                                    bp[idx3(kk, jj, k, j + 1, i2)] +
                                    bp[idx3(kk, jj, k, j + 1, i4)] +
                                    bp[idx3(kk, jj, k, j, i3)] +
                                    bp[idx3(kk, jj, k, j + 2, i3)] +
                                    bp[idx3(kk, jj, k - 1, j + 1, i3)] +
                                    bp[idx3(kk, jj, k + 1, j + 1, i3)],
                                (realk_c)1.0);
                    sb13 = min2(bp[idx3(kk, jj, k + 1, j, i3)] +
                                    bp[idx3(kk, jj, k + 1, j, i2)] +
                                    bp[idx3(kk, jj, k + 1, j, i4)] +
                                    bp[idx3(kk, jj, k + 1, j - 1, i3)] +
                                    bp[idx3(kk, jj, k + 1, j + 1, i3)] +
                                    bp[idx3(kk, jj, k, j, i3)] +
                                    bp[idx3(kk, jj, k + 2, j, i3)],
                                (realk_c)1.0);
                    sb14 = min2(bp[idx3(kk, jj, k + 1, j + 1, i3)] +
                                    bp[idx3(kk, jj, k + 1, j + 1, i2)] +
                                    bp[idx3(kk, jj, k + 1, j + 1, i4)] +
                                    bp[idx3(kk, jj, k + 1, j, i3)] +
                                    bp[idx3(kk, jj, k + 1, j + 2, i3)] +
                                    bp[idx3(kk, jj, k, j + 1, i3)] +
                                    bp[idx3(kk, jj, k + 2, j + 1, i3)],
                                (realk_c)1.0);
                }

                fak = (realk_c)1.0 /
                      ((ddy[j - 1] + ddy[j]) * (ddz[k - 1] + ddz[k])) *
                      (sb11 * ddy[j - 1] * ddz[k - 1] +
                       sb12 * ddy[j] * ddz[k - 1] + sb13 * ddy[j - 1] * ddz[k] +
                       sb14 * ddy[j] * ddz[k]);
                if (fak < (realk_c)0.1) {
                    fak = (realk_c)1.0;
                }
                fak = (realk_c)1.0 / fak;

                ubfine = (u[idx3(kk, jj, k, j, istag2)] *
                              ((realk_c)1.0 - sb11) * ddy[j - 1] * ddz[k - 1] +
                          u[idx3(kk, jj, k, j + 1, istag2)] *
                              ((realk_c)1.0 - sb12) * ddy[j] * ddz[k - 1] +
                          u[idx3(kk, jj, k + 1, j, istag2)] *
                              ((realk_c)1.0 - sb13) * ddy[j - 1] * ddz[k] +
                          u[idx3(kk, jj, k + 1, j + 1, istag2)] *
                              ((realk_c)1.0 - sb14) * ddy[j] * ddz[k]) /
                         ((ddy[j - 1] + ddy[j]) * (ddz[k - 1] + ddz[k]));

                u[idx3(kk, jj, k, j, istag2)] =
                    (ubuf[idx2(kk, k, j)] - ubfine) * fak * sb11 +
                    ((realk_c)1.0 - sb11) * u[idx3(kk, jj, k, j, istag2)] * sbu;
                u[idx3(kk, jj, k, j + 1, istag2)] =
                    (ubuf[idx2(kk, k, j + 1)] - ubfine) * fak * sb12 +
                    ((realk_c)1.0 - sb12) * u[idx3(kk, jj, k, j + 1, istag2)] *
                        sbu;
                u[idx3(kk, jj, k + 1, j, istag2)] =
                    (ubuf[idx2(kk, k + 1, j)] - ubfine) * fak * sb13 +
                    ((realk_c)1.0 - sb13) * u[idx3(kk, jj, k + 1, j, istag2)] *
                        sbu;
                u[idx3(kk, jj, k + 1, j + 1, istag2)] =
                    (ubuf[idx2(kk, k + 1, j + 1)] - ubfine) * fak * sb14 +
                    ((realk_c)1.0 - sb14) *
                        u[idx3(kk, jj, k + 1, j + 1, istag2)] * sbu;
            }
        }
        break;

    default:
        break;
    }

    if (ityp == 2) {
#if defined(_OPENMP)
#pragma omp for collapse(2)
#endif
        for (intk_c j = 2; j <= jj; ++j) {
            for (intk_c k = 2; k <= kk; ++k) {
                const realk_c vi =
                    (realk_c)0.25 *
                    (v[idx3(kk, jj, k, j - 1, i2)] + v[idx3(kk, jj, k, j, i2)] +
                     v[idx3(kk, jj, k, j - 1, i3)] + v[idx3(kk, jj, k, j, i3)]);
                const realk_c wi =
                    (realk_c)0.25 *
                    (w[idx3(kk, jj, k - 1, j, i2)] + w[idx3(kk, jj, k, j, i2)] +
                     w[idx3(kk, jj, k - 1, j, i3)] + w[idx3(kk, jj, k, j, i3)]);
                const realk_c umagsqr = u[idx3(kk, jj, k, j, istag2)] *
                                            u[idx3(kk, jj, k, j, istag2)] +
                                        vi * vi + wi * wi;
                p[idx3(kk, jj, k, j, i2)] =
                    pinf + min2((realk_c)0.0,
                                (realk_c)0.5 * rho *
                                    sign_c(umagsqr,
                                           (realk_c)dir *
                                               u[idx3(kk, jj, k, j, istag2)]));
            }
        }
    } else if (ityp == 1 || ityp == 3 || ityp == 4) {
#if defined(_OPENMP)
#pragma omp for collapse(2)
#endif
        for (intk_c j = 1; j <= jj; ++j) {
            for (intk_c k = 1; k <= kk; ++k) {
                p[idx3(kk, jj, k, j, i2)] = p[idx3(kk, jj, k, j, i3)];
            }
        }
    }
}

/* Y-normal boundary kernel (right/left faces). */
void bright_device_c(const intk_c kk, const intk_c jj, const intk_c ii,
                     const intk_c j2, const intk_c j3, const intk_c j4,
                     const intk_c jstag1, const intk_c jstag2, const intk_c dir,
                     const intk_c ityp, const realk_c pinf, realk_c *u,
                     realk_c *v, realk_c *w, realk_c *p, const realk_c *bp,
                     const realk_c *ubuf, const realk_c *vbuf,
                     const realk_c *wbuf, const realk_c *ddx,
                     const realk_c *ddz, const realk_c rho) {
    switch (ityp) {
    case 1:
#if defined(_OPENMP)
#pragma omp for collapse(2)
#endif
        for (intk_c i = 1; i <= ii; ++i) {
            for (intk_c k = 1; k <= kk; ++k) {
                u[idx3(kk, jj, k, j2, i)] =
                    (realk_c)2.0 * ubuf[idx2(kk, k, i)] -
                    u[idx3(kk, jj, k, j3, i)];
                v[idx3(kk, jj, k, jstag2, i)] = vbuf[idx2(kk, k, i)];
                w[idx3(kk, jj, k, j2, i)] =
                    (realk_c)2.0 * wbuf[idx2(kk, k, i)] -
                    w[idx3(kk, jj, k, j3, i)];
            }
        }
        break;

    case 2:
#if defined(_OPENMP)
#pragma omp for collapse(2)
#endif
        for (intk_c i = 1; i <= ii - 1; ++i) {
            for (intk_c k = 1; k <= kk - 1; ++k) {
                realk_c flag;

                flag =
                    sign_c((realk_c)1.0,
                           (realk_c)dir * (v[idx3(kk, jj, k, jstag2, i)] +
                                           v[idx3(kk, jj, k, jstag2, i + 1)]));
                flag = (realk_c)0.5 * (flag + (realk_c)1.0);
                u[idx3(kk, jj, k, j2, i)] =
                    flag * u[idx3(kk, jj, k, j3, i)] +
                    ((realk_c)1.0 - flag) *
                        ((realk_c)2.0 * ubuf[idx2(kk, k, i)] -
                         u[idx3(kk, jj, k, j3, i)]);

                v[idx3(kk, jj, k, jstag1, i)] = v[idx3(kk, jj, k, jstag2, i)];

                flag =
                    sign_c((realk_c)1.0,
                           (realk_c)dir * (v[idx3(kk, jj, k, jstag2, i)] +
                                           v[idx3(kk, jj, k + 1, jstag2, i)]));
                flag = (realk_c)0.5 * (flag + (realk_c)1.0);
                w[idx3(kk, jj, k, j2, i)] =
                    flag * w[idx3(kk, jj, k, j3, i)] +
                    ((realk_c)1.0 - flag) *
                        ((realk_c)2.0 * wbuf[idx2(kk, k, i)] -
                         w[idx3(kk, jj, k, j3, i)]);
            }
        }
        break;

    case 3:
#if defined(_OPENMP)
#pragma omp for collapse(2)
#endif
        for (intk_c i = 1; i <= ii; ++i) {
            for (intk_c k = 1; k <= kk; ++k) {
                u[idx3(kk, jj, k, j2, i)] = (realk_c)0.0;
                v[idx3(kk, jj, k, jstag2, i)] = (realk_c)0.0;
                w[idx3(kk, jj, k, j2, i)] = (realk_c)0.0;
            }
        }
        break;

    case 4:
#if defined(_OPENMP)
#pragma omp for collapse(2)
#endif
        for (intk_c i = 1; i <= ii; ++i) {
            for (intk_c k = 1; k <= kk; ++k) {
                u[idx3(kk, jj, k, j2, i)] = u[idx3(kk, jj, k, j3, i)];
                v[idx3(kk, jj, k, jstag2, i)] = (realk_c)0.0;
                w[idx3(kk, jj, k, j2, i)] = w[idx3(kk, jj, k, j3, i)];
            }
        }
        break;

    case 5:
#if defined(_OPENMP)
#pragma omp for collapse(2)
#endif
        for (intk_c i = 2; i <= ii - 1; ++i) {
            for (intk_c k = 2; k <= kk - 1; ++k) {
                const realk_c sbu =
                    bp[idx3(kk, jj, k, j2, i)] * bp[idx3(kk, jj, k, j2, i + 1)];
                const realk_c sbw =
                    bp[idx3(kk, jj, k, j2, i)] * bp[idx3(kk, jj, k + 1, j2, i)];
                u[idx3(kk, jj, k, j2, i)] =
                    (ubuf[idx2(kk, k, i)] -
                     (realk_c)0.5 * (u[idx3(kk, jj, k, j3, i)] -
                                     u[idx3(kk, jj, k, j4, i)])) *
                        sbu +
                    ((realk_c)1.0 - sbu) * u[idx3(kk, jj, k, j2, i)];
                w[idx3(kk, jj, k, j2, i)] =
                    (wbuf[idx2(kk, k, i)] -
                     (realk_c)0.5 * (w[idx3(kk, jj, k, j2, i)] -
                                     w[idx3(kk, jj, k, j4, i)])) *
                        sbw +
                    ((realk_c)1.0 - sbw) * w[idx3(kk, jj, k, j2, i)];
            }
        }

#if defined(_OPENMP)
#pragma omp for collapse(2)
#endif
        for (intk_c i = 3; i <= ii - 2; i += 2) {
            for (intk_c k = 3; k <= kk - 2; k += 2) {
                const realk_c sbv = bp[idx3(kk, jj, k, j2, i)];
                realk_c sb11;
                realk_c sb12;
                realk_c sb13;
                realk_c sb14;
                realk_c fak;
                realk_c vbfine;

                if (sbv > (realk_c)0.5) {
                    sb11 =
                        bp[idx3(kk, jj, k, j2, i)] * bp[idx3(kk, jj, k, j3, i)];
                    sb12 = bp[idx3(kk, jj, k, j2, i + 1)] *
                           bp[idx3(kk, jj, k, j3, i + 1)];
                    sb13 = bp[idx3(kk, jj, k + 1, j2, i)] *
                           bp[idx3(kk, jj, k + 1, j3, i)];
                    sb14 = bp[idx3(kk, jj, k + 1, j2, i + 1)] *
                           bp[idx3(kk, jj, k + 1, j3, i + 1)];
                } else {
                    sb11 = min2(bp[idx3(kk, jj, k, j3, i)] +
                                    bp[idx3(kk, jj, k, j2, i)] +
                                    bp[idx3(kk, jj, k, j4, i)] +
                                    bp[idx3(kk, jj, k, j3, i - 1)] +
                                    bp[idx3(kk, jj, k, j3, i + 1)] +
                                    bp[idx3(kk, jj, k - 1, j3, i)] +
                                    bp[idx3(kk, jj, k + 1, j3, i)],
                                (realk_c)1.0);
                    sb12 = min2(bp[idx3(kk, jj, k, j3, i + 1)] +
                                    bp[idx3(kk, jj, k, j2, i + 1)] +
                                    bp[idx3(kk, jj, k, j4, i + 1)] +
                                    bp[idx3(kk, jj, k, j3, i)] +
                                    bp[idx3(kk, jj, k, j3, i + 2)] +
                                    bp[idx3(kk, jj, k - 1, j3, i + 1)] +
                                    bp[idx3(kk, jj, k + 1, j3, i + 1)],
                                (realk_c)1.0);
                    sb13 = min2(bp[idx3(kk, jj, k + 1, j3, i)] +
                                    bp[idx3(kk, jj, k + 1, j2, i)] +
                                    bp[idx3(kk, jj, k + 1, j4, i)] +
                                    bp[idx3(kk, jj, k + 1, j3, i - 1)] +
                                    bp[idx3(kk, jj, k + 1, j3, i + 1)] +
                                    bp[idx3(kk, jj, k, j3, i)] +
                                    bp[idx3(kk, jj, k + 2, j3, i)],
                                (realk_c)1.0);
                    sb14 = min2(bp[idx3(kk, jj, k + 1, j3, i + 1)] +
                                    bp[idx3(kk, jj, k + 1, j2, i + 1)] +
                                    bp[idx3(kk, jj, k + 1, j4, i + 1)] +
                                    bp[idx3(kk, jj, k + 1, j3, i)] +
                                    bp[idx3(kk, jj, k + 1, j3, i + 2)] +
                                    bp[idx3(kk, jj, k, j3, i + 1)] +
                                    bp[idx3(kk, jj, k + 2, j3, i + 1)],
                                (realk_c)1.0);
                }

                fak = (realk_c)1.0 /
                      ((ddz[k - 1] + ddz[k]) * (ddx[i - 1] + ddx[i])) *
                      (sb11 * ddz[k - 1] * ddx[i - 1] +
                       sb12 * ddz[k - 1] * ddx[i] + sb13 * ddz[k] * ddx[i - 1] +
                       sb14 * ddz[k] * ddx[i]);
                if (fak < (realk_c)0.1) {
                    fak = (realk_c)1.0;
                }
                fak = (realk_c)1.0 / fak;

                vbfine = (v[idx3(kk, jj, k, jstag2, i)] *
                              ((realk_c)1.0 - sb11) * ddz[k - 1] * ddx[i - 1] +
                          v[idx3(kk, jj, k, jstag2, i + 1)] *
                              ((realk_c)1.0 - sb12) * ddz[k - 1] * ddx[i] +
                          v[idx3(kk, jj, k + 1, jstag2, i)] *
                              ((realk_c)1.0 - sb13) * ddz[k] * ddx[i - 1] +
                          v[idx3(kk, jj, k + 1, jstag2, i + 1)] *
                              ((realk_c)1.0 - sb14) * ddz[k] * ddx[i]) /
                         ((ddz[k - 1] + ddz[k]) * (ddx[i - 1] + ddx[i]));

                v[idx3(kk, jj, k, jstag2, i)] =
                    (vbuf[idx2(kk, k, i)] - vbfine) * fak * sb11 +
                    v[idx3(kk, jj, k, jstag2, i)] * ((realk_c)1.0 - sb11) * sbv;
                v[idx3(kk, jj, k, jstag2, i + 1)] =
                    (vbuf[idx2(kk, k, i + 1)] - vbfine) * fak * sb12 +
                    v[idx3(kk, jj, k, jstag2, i + 1)] * ((realk_c)1.0 - sb12) *
                        sbv;
                v[idx3(kk, jj, k + 1, jstag2, i)] =
                    (vbuf[idx2(kk, k + 1, i)] - vbfine) * fak * sb13 +
                    v[idx3(kk, jj, k + 1, jstag2, i)] * ((realk_c)1.0 - sb13) *
                        sbv;
                v[idx3(kk, jj, k + 1, jstag2, i + 1)] =
                    (vbuf[idx2(kk, k + 1, i + 1)] - vbfine) * fak * sb14 +
                    v[idx3(kk, jj, k + 1, jstag2, i + 1)] *
                        ((realk_c)1.0 - sb14) * sbv;
            }
        }
        break;

    default:
        break;
    }

    if (ityp == 2) {
#if defined(_OPENMP)
#pragma omp for collapse(2)
#endif
        for (intk_c i = 2; i <= ii; ++i) {
            for (intk_c k = 2; k <= kk; ++k) {
                const realk_c ui =
                    (realk_c)0.25 *
                    (u[idx3(kk, jj, k, j2, i - 1)] + u[idx3(kk, jj, k, j2, i)] +
                     u[idx3(kk, jj, k, j3, i - 1)] + u[idx3(kk, jj, k, j3, i)]);
                const realk_c wi =
                    (realk_c)0.25 *
                    (w[idx3(kk, jj, k - 1, j2, i)] + w[idx3(kk, jj, k, j2, i)] +
                     w[idx3(kk, jj, k - 1, j3, i)] + w[idx3(kk, jj, k, j3, i)]);
                const realk_c umagsqr = ui * ui +
                                        v[idx3(kk, jj, k, jstag2, i)] *
                                            v[idx3(kk, jj, k, jstag2, i)] +
                                        wi * wi;
                p[idx3(kk, jj, k, j2, i)] =
                    pinf + min2((realk_c)0.0,
                                (realk_c)0.5 * rho *
                                    sign_c(umagsqr,
                                           (realk_c)dir *
                                               v[idx3(kk, jj, k, jstag2, i)]));
            }
        }
    } else if (ityp == 1 || ityp == 3 || ityp == 4) {
#if defined(_OPENMP)
#pragma omp for collapse(2)
#endif
        for (intk_c i = 1; i <= ii; ++i) {
            for (intk_c k = 1; k <= kk; ++k) {
                p[idx3(kk, jj, k, j2, i)] = p[idx3(kk, jj, k, j3, i)];
            }
        }
    }
}

/* Z-normal boundary kernel (bottom/top faces). */
void bbottom_device_c(const intk_c kk, const intk_c jj, const intk_c ii,
                      const intk_c k2, const intk_c k3, const intk_c k4,
                      const intk_c kstag1, const intk_c kstag2,
                      const intk_c dir, const intk_c ityp, const realk_c pinf,
                      realk_c *u, realk_c *v, realk_c *w, realk_c *p,
                      const realk_c *bp, const realk_c *ubuf,
                      const realk_c *vbuf, const realk_c *wbuf,
                      const realk_c *ddx, const realk_c *ddy,
                      const realk_c rho) {
    switch (ityp) {
    case 1:
#if defined(_OPENMP)
#pragma omp for collapse(2)
#endif
        for (intk_c i = 1; i <= ii; ++i) {
            for (intk_c j = 1; j <= jj; ++j) {
                u[idx3(kk, jj, k2, j, i)] =
                    (realk_c)2.0 * ubuf[idx2(jj, j, i)] -
                    u[idx3(kk, jj, k3, j, i)];
                v[idx3(kk, jj, k2, j, i)] =
                    (realk_c)2.0 * vbuf[idx2(jj, j, i)] -
                    v[idx3(kk, jj, k3, j, i)];
                w[idx3(kk, jj, kstag2, j, i)] = wbuf[idx2(jj, j, i)];
            }
        }
        break;

    case 2:
#if defined(_OPENMP)
#pragma omp for collapse(2)
#endif
        for (intk_c i = 1; i <= ii - 1; ++i) {
            for (intk_c j = 1; j <= jj - 1; ++j) {
                realk_c flag;

                flag =
                    sign_c((realk_c)1.0,
                           (realk_c)dir * (w[idx3(kk, jj, kstag2, j, i)] +
                                           w[idx3(kk, jj, kstag2, j, i + 1)]));
                flag = (realk_c)0.5 * (flag + (realk_c)1.0);
                u[idx3(kk, jj, k2, j, i)] =
                    flag * u[idx3(kk, jj, k3, j, i)] +
                    ((realk_c)1.0 - flag) *
                        ((realk_c)2.0 * ubuf[idx2(jj, j, i)] -
                         u[idx3(kk, jj, k3, j, i)]);

                flag =
                    sign_c((realk_c)1.0,
                           (realk_c)dir * (w[idx3(kk, jj, kstag2, j, i)] +
                                           w[idx3(kk, jj, kstag2, j + 1, i)]));
                flag = (realk_c)0.5 * (flag + (realk_c)1.0);
                v[idx3(kk, jj, k2, j, i)] =
                    flag * v[idx3(kk, jj, k3, j, i)] +
                    ((realk_c)1.0 - flag) *
                        ((realk_c)2.0 * vbuf[idx2(jj, j, i)] -
                         v[idx3(kk, jj, k3, j, i)]);

                w[idx3(kk, jj, kstag1, j, i)] = w[idx3(kk, jj, kstag2, j, i)];
            }
        }
        break;

    case 3:
#if defined(_OPENMP)
#pragma omp for collapse(2)
#endif
        for (intk_c i = 1; i <= ii; ++i) {
            for (intk_c j = 1; j <= jj; ++j) {
                u[idx3(kk, jj, k2, j, i)] = (realk_c)0.0;
                v[idx3(kk, jj, k2, j, i)] = (realk_c)0.0;
                w[idx3(kk, jj, kstag2, j, i)] = (realk_c)0.0;
            }
        }
        break;

    case 4:
#if defined(_OPENMP)
#pragma omp for collapse(2)
#endif
        for (intk_c i = 1; i <= ii; ++i) {
            for (intk_c j = 1; j <= jj; ++j) {
                u[idx3(kk, jj, k2, j, i)] = u[idx3(kk, jj, k3, j, i)];
                v[idx3(kk, jj, k2, j, i)] = v[idx3(kk, jj, k3, j, i)];
                w[idx3(kk, jj, kstag2, j, i)] = (realk_c)0.0;
            }
        }
        break;

    case 5:
#if defined(_OPENMP)
#pragma omp for collapse(2)
#endif
        for (intk_c i = 2; i <= ii - 1; ++i) {
            for (intk_c j = 2; j <= jj - 1; ++j) {
                const realk_c sbu =
                    bp[idx3(kk, jj, k2, j, i)] * bp[idx3(kk, jj, k2, j, i + 1)];
                const realk_c sbv =
                    bp[idx3(kk, jj, k2, j, i)] * bp[idx3(kk, jj, k2, j + 1, i)];
                u[idx3(kk, jj, k2, j, i)] =
                    (ubuf[idx2(jj, j, i)] -
                     (realk_c)0.5 * (u[idx3(kk, jj, k3, j, i)] -
                                     u[idx3(kk, jj, k4, j, i)])) *
                        sbu +
                    ((realk_c)1.0 - sbu) * u[idx3(kk, jj, k2, j, i)];
                v[idx3(kk, jj, k2, j, i)] =
                    (vbuf[idx2(jj, j, i)] -
                     (realk_c)0.5 * (v[idx3(kk, jj, k3, j, i)] -
                                     v[idx3(kk, jj, k4, j, i)])) *
                        sbv +
                    ((realk_c)1.0 - sbv) * v[idx3(kk, jj, k2, j, i)];
            }
        }

#if defined(_OPENMP)
#pragma omp for collapse(2)
#endif
        for (intk_c i = 3; i <= ii - 2; i += 2) {
            for (intk_c j = 3; j <= jj - 2; j += 2) {
                const realk_c sbw = bp[idx3(kk, jj, k2, j, i)];
                realk_c sb11;
                realk_c sb12;
                realk_c sb13;
                realk_c sb14;
                realk_c fak;
                realk_c wbfine;

                if (sbw > (realk_c)0.5) {
                    sb11 =
                        bp[idx3(kk, jj, k2, j, i)] * bp[idx3(kk, jj, k3, j, i)];
                    sb12 = bp[idx3(kk, jj, k2, j, i + 1)] *
                           bp[idx3(kk, jj, k3, j, i + 1)];
                    sb13 = bp[idx3(kk, jj, k2, j + 1, i)] *
                           bp[idx3(kk, jj, k3, j + 1, i)];
                    sb14 = bp[idx3(kk, jj, k2, j + 1, i + 1)] *
                           bp[idx3(kk, jj, k3, j + 1, i + 1)];
                } else {
                    sb11 = min2(bp[idx3(kk, jj, k3, j, i)] +
                                    bp[idx3(kk, jj, k2, j, i)] +
                                    bp[idx3(kk, jj, k4, j, i)] +
                                    bp[idx3(kk, jj, k3, j, i - 1)] +
                                    bp[idx3(kk, jj, k3, j, i + 1)] +
                                    bp[idx3(kk, jj, k3, j - 1, i)] +
                                    bp[idx3(kk, jj, k3, j + 1, i)],
                                (realk_c)1.0);
                    sb12 = min2(bp[idx3(kk, jj, k3, j, i + 1)] +
                                    bp[idx3(kk, jj, k2, j, i + 1)] +
                                    bp[idx3(kk, jj, k4, j, i + 1)] +
                                    bp[idx3(kk, jj, k3, j, i)] +
                                    bp[idx3(kk, jj, k3, j, i + 2)] +
                                    bp[idx3(kk, jj, k3, j - 1, i + 1)] +
                                    bp[idx3(kk, jj, k3, j + 1, i + 1)],
                                (realk_c)1.0);
                    sb13 = min2(bp[idx3(kk, jj, k3, j + 1, i)] +
                                    bp[idx3(kk, jj, k2, j + 1, i)] +
                                    bp[idx3(kk, jj, k4, j + 1, i)] +
                                    bp[idx3(kk, jj, k3, j + 1, i - 1)] +
                                    bp[idx3(kk, jj, k3, j + 1, i + 1)] +
                                    bp[idx3(kk, jj, k3, j, i)] +
                                    bp[idx3(kk, jj, k3, j + 2, i)],
                                (realk_c)1.0);
                    sb14 = min2(bp[idx3(kk, jj, k3, j + 1, i + 1)] +
                                    bp[idx3(kk, jj, k2, j + 1, i + 1)] +
                                    bp[idx3(kk, jj, k4, j + 1, i + 1)] +
                                    bp[idx3(kk, jj, k3, j + 1, i)] +
                                    bp[idx3(kk, jj, k3, j + 1, i + 2)] +
                                    bp[idx3(kk, jj, k3, j, i + 1)] +
                                    bp[idx3(kk, jj, k3, j + 2, i + 1)],
                                (realk_c)1.0);
                }

                fak = (realk_c)1.0 /
                      ((ddy[j - 1] + ddy[j]) * (ddx[i - 1] + ddx[i])) *
                      (sb11 * ddy[j - 1] * ddx[i - 1] +
                       sb12 * ddy[j - 1] * ddx[i] + sb13 * ddy[j] * ddx[i - 1] +
                       sb14 * ddy[j] * ddx[i]);
                if (fak < (realk_c)0.1) {
                    fak = (realk_c)1.0;
                }
                fak = (realk_c)1.0 / fak;

                wbfine = (w[idx3(kk, jj, kstag2, j, i)] *
                              ((realk_c)1.0 - sb11) * ddy[j - 1] * ddx[i - 1] +
                          w[idx3(kk, jj, kstag2, j, i + 1)] *
                              ((realk_c)1.0 - sb12) * ddy[j - 1] * ddx[i] +
                          w[idx3(kk, jj, kstag2, j + 1, i)] *
                              ((realk_c)1.0 - sb13) * ddy[j] * ddx[i - 1] +
                          w[idx3(kk, jj, kstag2, j + 1, i + 1)] *
                              ((realk_c)1.0 - sb14) * ddy[j] * ddx[i]) /
                         ((ddy[j - 1] + ddy[j]) * (ddx[i - 1] + ddx[i]));

                w[idx3(kk, jj, kstag2, j, i)] =
                    (wbuf[idx2(jj, j, i)] - wbfine) * fak * sb11 +
                    w[idx3(kk, jj, kstag2, j, i)] * ((realk_c)1.0 - sb11) * sbw;
                w[idx3(kk, jj, kstag2, j, i + 1)] =
                    (wbuf[idx2(jj, j, i + 1)] - wbfine) * fak * sb12 +
                    w[idx3(kk, jj, kstag2, j, i + 1)] * ((realk_c)1.0 - sb12) *
                        sbw;
                w[idx3(kk, jj, kstag2, j + 1, i)] =
                    (wbuf[idx2(jj, j + 1, i)] - wbfine) * fak * sb13 +
                    w[idx3(kk, jj, kstag2, j + 1, i)] * ((realk_c)1.0 - sb13) *
                        sbw;
                w[idx3(kk, jj, kstag2, j + 1, i + 1)] =
                    (wbuf[idx2(jj, j + 1, i + 1)] - wbfine) * fak * sb14 +
                    w[idx3(kk, jj, kstag2, j + 1, i + 1)] *
                        ((realk_c)1.0 - sb14) * sbw;
            }
        }
        break;

    default:
        break;
    }

    if (ityp == 2) {
#if defined(_OPENMP)
#pragma omp for collapse(2)
#endif
        for (intk_c i = 2; i <= ii; ++i) {
            for (intk_c j = 2; j <= jj; ++j) {
                const realk_c ui =
                    (realk_c)0.25 *
                    (u[idx3(kk, jj, k2, j, i - 1)] + u[idx3(kk, jj, k2, j, i)] +
                     u[idx3(kk, jj, k3, j, i - 1)] + u[idx3(kk, jj, k3, j, i)]);
                const realk_c vi =
                    (realk_c)0.25 *
                    (v[idx3(kk, jj, k2, j - 1, i)] + v[idx3(kk, jj, k2, j, i)] +
                     v[idx3(kk, jj, k3, j - 1, i)] + v[idx3(kk, jj, k3, j, i)]);
                const realk_c umagsqr = ui * ui + vi * vi +
                                        w[idx3(kk, jj, kstag2, j, i)] *
                                            w[idx3(kk, jj, kstag2, j, i)];
                p[idx3(kk, jj, k2, j, i)] =
                    pinf + min2((realk_c)0.0,
                                (realk_c)0.5 * rho *
                                    sign_c(umagsqr,
                                           (realk_c)dir *
                                               w[idx3(kk, jj, kstag2, j, i)]));
            }
        }
    } else if (ityp == 1 || ityp == 3 || ityp == 4) {
#if defined(_OPENMP)
#pragma omp for collapse(2)
#endif
        for (intk_c i = 1; i <= ii; ++i) {
            for (intk_c j = 1; j <= jj; ++j) {
                p[idx3(kk, jj, k2, j, i)] = p[idx3(kk, jj, k3, j, i)];
            }
        }
    }
}

/* Unified boundary kernel: face and boundary-type selection happen in C. */
void bound_flow_device_c(const intk_c kk, const intk_c jj, const intk_c ii,
                         const intk_c iface, const intk_c ityp,
                         const realk_c pinf, realk_c *u, realk_c *v,
                         realk_c *w, realk_c *p, const realk_c *bp,
                         realk_c *ubuf, realk_c *vbuf, realk_c *wbuf,
                         const realk_c *ddx, const realk_c *ddy,
                         const realk_c *ddz, const realk_c rho) {

    /* NOS/SLI do not rely on boundary buffers, keep the no-buffer path. */
    if (ityp == 3 || ityp == 4) {
        bound_nobuffer_device_c(kk, jj, ii, iface, ityp, u, v, w, p);
        return;
    }

    switch (iface) {
    case 1:
        bfront_device_c(kk, jj, ii, 2, 3, 4, 1, 2, -1, ityp, pinf, u, v, w, p,
                        bp, ubuf, vbuf, wbuf, ddy, ddz, rho);
        break;

    case 2:
        bfront_device_c(kk, jj, ii, ii - 1, ii - 2, ii - 3, ii - 1, ii - 2, 1,
                        ityp, pinf, u, v, w, p, bp, ubuf, vbuf, wbuf, ddy, ddz,
                        rho);
        break;

    case 3:
        bright_device_c(kk, jj, ii, 2, 3, 4, 1, 2, -1, ityp, pinf, u, v, w, p,
                        bp, ubuf, vbuf, wbuf, ddx, ddz, rho);
        break;

    case 4:
        bright_device_c(kk, jj, ii, jj - 1, jj - 2, jj - 3, jj - 1, jj - 2, 1,
                        ityp, pinf, u, v, w, p, bp, ubuf, vbuf, wbuf, ddx, ddz,
                        rho);
        break;

    case 5:
        bbottom_device_c(kk, jj, ii, 2, 3, 4, 1, 2, -1, ityp, pinf, u, v, w, p,
                         bp, ubuf, vbuf, wbuf, ddx, ddy, rho);
        break;

    case 6:
        bbottom_device_c(kk, jj, ii, kk - 1, kk - 2, kk - 3, kk - 1, kk - 2,
                         1, ityp, pinf, u, v, w, p, bp, ubuf, vbuf, wbuf, ddx,
                         ddy, rho);
        break;

    default:
        break;
    }
}

#if defined(_OPENMP)
#pragma omp end declare target
#endif

void apply_bound_flow_impl_c(const intk_c ntasks, realk_c *u, realk_c *v,
                             realk_c *w, realk_c *p, realk_c *ubuffer,
                             realk_c *vbuffer, realk_c *wbuffer,
                             const realk_c *bp, const realk_c *ddx,
                             const realk_c *ddy, const realk_c *ddz,
                             const bound_flow_task_t *tasks,
                             const realk_c rho_in) {
#pragma omp target teams distribute
    for (intk_c itask = 0; itask < ntasks; ++itask) {
        const bound_flow_task_t *const t = &tasks[itask];

        realk_c *const u_g = &u[(size_t)(t->ip3 - 1)];
        realk_c *const v_g = &v[(size_t)(t->ip3 - 1)];
        realk_c *const w_g = &w[(size_t)(t->ip3 - 1)];
        realk_c *const p_g = &p[(size_t)(t->ip3 - 1)];
        const realk_c *const bp_g = &bp[(size_t)(t->ip3 - 1)];

        realk_c *const ubuf_g = &ubuffer[(size_t)(t->ipbb - 1)];
        realk_c *const vbuf_g = &vbuffer[(size_t)(t->ipbb - 1)];
        realk_c *const wbuf_g = &wbuffer[(size_t)(t->ipbb - 1)];

        const realk_c *const ddx_g = &ddx[(size_t)(t->ipx - 1)];
        const realk_c *const ddy_g = &ddy[(size_t)(t->ipy - 1)];
        const realk_c *const ddz_g = &ddz[(size_t)(t->ipz - 1)];

#pragma omp parallel
        {
            bound_flow_device_c(t->kk, t->jj, t->ii, t->iface, t->ityp,
                                t->pinf, u_g, v_g, w_g, p_g, bp_g, ubuf_g,
                                vbuf_g, wbuf_g, ddx_g, ddy_g, ddz_g, rho_in);
        }
    }
}
