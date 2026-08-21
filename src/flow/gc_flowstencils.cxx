#include "mglet_precision.h"

#include <cstdlib>
#include <limits>

constexpr mgletint nvelpts = 6;

struct flowstencil_t {
    mgletint icell;
    mgletint npts;
    mgletint pts[nvelpts];
    mgletreal coeff[nvelpts];
    mgletreal acoeff;
    mgletreal oldsol;
};

struct fcorrstencil_t {
    mgletint pts[6];
    mgletreal area[6];
    mgletreal darea[3];
    mgletreal dvol;
    mgletreal acoeff;
};

extern "C" void copy_field_values_c(const mgletint narr,
        mgletreal* dstarr, const mgletreal* srcarr,
        const mgletint nbuffers, mgletreal* dstbuffers,
        const mgletreal* srcbuffers) {
    #pragma omp target teams loop
    for (mgletint i = 0; i < narr; ++i) {
        dstarr[i] = srcarr[i];
    }

    #pragma omp target teams loop
    for (mgletint i = 0; i < nbuffers; ++i) {
        dstbuffers[i] = srcbuffers[i];
    }
}

static void apply_flowstencils(const mgletint cmp, const char ityp,
    const mgletint first, const mgletint last, flowstencil_t* stencils,
        mgletreal* u, mgletreal* v, mgletreal* w) {
    mgletreal* var;
    if (cmp == 1) {
        var = u;
    } else if (cmp == 2) {
        var = v;
    } else if (cmp == 3) {
        var = w;
    } else {
        std::abort();
    }

    if (ityp == 'F' || ityp == 'X') {
        #pragma omp target teams loop
        for (mgletint i = first; i <= last; ++i) {
            auto flux = mgletreal{0.0};
            for (mgletint n = 0; n < stencils[i].npts; ++n) {
                const auto pt = stencils[i].pts[n] - 1;
                flux += var[pt] * stencils[i].coeff[n];
            }
            const auto icell = stencils[i].icell - 1;
            var[icell] = flux + stencils[i].acoeff;
        }
    } else if (ityp == 'R' || ityp == 'Y') {
        #pragma omp target teams loop
        for (mgletint i = first; i <= last; ++i) {
            const auto icell = stencils[i].icell - 1;
            var[icell] = stencils[i].oldsol;
        }
    } else if (ityp == 'S' || ityp == 'Z') {
        #pragma omp target teams loop
        for (mgletint i = first; i <= last; ++i) {
            const auto icell = stencils[i].icell - 1;
            stencils[i].oldsol = var[icell];
        }
    } else {
        std::abort();
    }
}

extern "C" void wmxpolquad_c(const mgletint cmp, const char ityp,
        const mgletint first, const mgletint last, flowstencil_t* stencils,
        mgletreal* u, mgletreal* v, mgletreal* w) {
    apply_flowstencils(cmp, ityp, first, last, stencils, u, v, w);
}

extern "C" void wmxpolquadvel_c(const mgletint cmp, const char ityp,
        const mgletint first, const mgletint last, flowstencil_t* stencils,
        mgletreal* u, mgletreal* v, mgletreal* w) {
    apply_flowstencils(cmp, ityp, first, last, stencils, u, v, w);
}

extern "C" void wmxpolquadfcorr_c(const mgletint first,
        const mgletint last, const fcorrstencil_t* stencils,
        mgletreal* u, mgletreal* v, mgletreal* w) {
    #pragma omp target teams loop
    for (mgletint i = first; i <= last; ++i) {
        const auto ax1 = stencils[i].area[0];
        const auto ax2 = stencils[i].area[1];
        const auto ay1 = stencils[i].area[2];
        const auto ay2 = stencils[i].area[3];
        const auto az1 = stencils[i].area[4];
        const auto az2 = stencils[i].area[5];

        const auto u1 = stencils[i].pts[0] - 1;
        const auto u2 = stencils[i].pts[1] - 1;
        const auto v1 = stencils[i].pts[2] - 1;
        const auto v2 = stencils[i].pts[3] - 1;
        const auto w1 = stencils[i].pts[4] - 1;
        const auto w2 = stencils[i].pts[5] - 1;

        auto div = stencils[i].darea[0] * (u[u1] - u[u2]);
        div += stencils[i].darea[1] * (v[v1] - v[v2]);
        div += stencils[i].darea[2] * (w[w1] - w[w2]);
        div += stencils[i].acoeff;

        const auto sarea = ax1 + ax2 + ay1 + ay2 + az1 + az2;
        if (sarea < std::numeric_limits<mgletreal>::min()) {
            std::abort();
        }
        div /= sarea;

        u[u1] -= ax1 * div / stencils[i].darea[0];
        u[u2] += ax2 * div / stencils[i].darea[0];
        v[v1] -= ay1 * div / stencils[i].darea[1];
        v[v2] += ay2 * div / stencils[i].darea[1];
        w[w1] -= az1 * div / stencils[i].darea[2];
        w[w2] += az2 * div / stencils[i].darea[2];
    }
}
