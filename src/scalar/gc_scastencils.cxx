#include "mglet_precision.h"

constexpr mgletint nperi = 6;

struct scalar_bc_t {
    mgletint flag;
    mgletreal value;
};

struct scastencil_t {
    mgletint icell;
    mgletint body;
    mgletint npts;
    mgletint pts[nperi];
    mgletreal areabyvol;
    mgletreal acoeff;
    mgletreal coeff[nperi];
    mgletreal acoeffp;
    mgletreal coeffp[nperi];
};

extern "C" void set_scastencils_qtt_c(const mgletint nscastencils, 
        const scastencil_t* scastencils, const scalar_bc_t* geometries, 
        mgletreal* qtt) {
    #pragma omp target teams loop
    for (mgletint istencil = 0; istencil < nscastencils; ++istencil) {
        // 1-based Fortran index, but the 0th index in the Fortran array is the
        // default boundary condition. No need to subtract 1 from the index.
        const auto body = scastencils[istencil].body;
        const auto bctype = geometries[body].flag;

        if (bctype == 1) {
            const auto bcval = geometries[body].value;
            
            // 1-based Fortran index
            const auto icell = scastencils[istencil].icell - 1;
            qtt[icell] = qtt[icell] + scastencils[istencil].areabyvol * bcval;
        }
    }
}

extern "C" void set_scastencils_t_c(const char ctyp, mgletint nscastencils,
        const scastencil_t* scastencils, const scalar_bc_t* geometries,
        mgletreal* t) {
    if (ctyp == 'C') {
        #pragma omp target teams loop
        for (mgletint istencil = 0; istencil < nscastencils; ++istencil) {
            // 1-based Fortran index, but the 0th index in the Fortran array is 
            // the default boundary condition. No need to subtract 1 from the 
            // index.
            const auto body = scastencils[istencil].body;
            const auto bctype = geometries[body].flag;

            auto val = mgletreal{0.0};
            for (mgletint n = 0; n < scastencils[istencil].npts; ++n) {
                // 1-based Fortran index
                const auto pt = scastencils[istencil].pts[n] - 1;
                val += t[pt] * scastencils[istencil].coeff[n];
            }

            if (bctype == 0) {
                const auto bcval = geometries[body].value;
                // 1-based Fortran index
                const auto icell = scastencils[istencil].icell - 1;
                t[icell] = val + scastencils[istencil].acoeff * bcval;
            }
        }
    } else if (ctyp == 'P') {
        #pragma omp target teams loop
        for (mgletint istencil = 0; istencil < nscastencils; ++istencil) {
            // 1-based Fortran index, but the 0th index in the Fortran array is 
            // the default boundary condition. No need to subtract 1 from the 
            // index.
            const auto body = scastencils[istencil].body;
            const auto bctype = geometries[body].flag;

            auto val = mgletreal{0.0};
            for (mgletint n = 0; n < scastencils[istencil].npts; ++n) {
                // 1-based Fortran index
                const auto pt = scastencils[istencil].pts[n] - 1;
                val += t[pt] * scastencils[istencil].coeffp[n];
            }

            if (bctype == 0) {
                const auto bcval = geometries[body].value;
                // 1-based Fortran index
                const auto icell = scastencils[istencil].icell - 1;
                t[icell] = val + scastencils[istencil].acoeffp * bcval;
            }
        }
    }
}
