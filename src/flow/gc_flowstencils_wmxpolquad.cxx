#include "mglet_precision.h"

#include <cstdlib>

constexpr mgletint nvelpts = 6;

struct flowstencil_t {
    mgletint icell;
    mgletint npts;
    mgletint pts[nvelpts];
    mgletreal coeff[nvelpts];
    mgletreal acoeff;
    mgletreal oldsol;
};

extern "C" void wmxpolquad_c(int iop, size_t first, size_t count,
        flowstencil_t* stencils, mgletreal* var) {
    if (count == 0) {
        return;
    }

    if (iop == 'F' || iop == 'X') {
        #pragma omp target teams loop
        for (size_t i = first; i < first + count; ++i) {
            auto flux = mgletreal{0.0};
            for (mgletint j = 0; j < stencils[i].npts; ++j) {
                const auto pt = stencils[i].pts[j] - 1;
                flux += var[pt] * stencils[i].coeff[j];
            }
            const auto icell = stencils[i].icell - 1;
            var[icell] = flux + stencils[i].acoeff;
        }
    } else if (iop == 'R' || iop == 'Y') {
        #pragma omp target teams loop
        for (size_t i = first; i < first + count; ++i) {
            const auto icell = stencils[i].icell - 1;
            var[icell] = stencils[i].oldsol;
        }
    } else if (iop == 'S' || iop == 'Z') {
        #pragma omp target teams loop
        for (size_t i = first; i < first + count; ++i) {
            const auto icell = stencils[i].icell - 1;
            stencils[i].oldsol = var[icell];
        }
    } else {
        std::abort();
    }
}