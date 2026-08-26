#include "mglet_precision.h"

#include <cstdlib>
#include <limits>

struct fcorrstencil_t {
    mgletint pts[6];
    mgletreal area[6];
    mgletreal darea[3];
    mgletreal dvol;
    mgletreal acoeff;
};

extern "C" void wmxpolquadfcorr_c(size_t first, size_t count,
        const fcorrstencil_t* stencils, mgletreal* u, mgletreal* v,
        mgletreal* w) {
    if (count == 0) {
        return;
    }

    #pragma omp target teams loop
    for (size_t i = first; i < first + count; ++i) {
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