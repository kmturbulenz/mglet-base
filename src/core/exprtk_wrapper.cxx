#include "mglet_precision.h"
#include "exprtk_wrapper.h"

#include "exprtk.hpp"

#include <cassert>
#include <cstdio>
#include <random>
#include <string>

#include <ISO_Fortran_binding.h>


template <typename T>
void eval_expr(CFI_cdesc_t* res, const char* name, const char* expr,
               T rho, T gmol, T tu_level, T timeph,
               CFI_cdesc_t* x, CFI_cdesc_t* y, CFI_cdesc_t* z,
               CFI_cdesc_t* dx, CFI_cdesc_t* dy, CFI_cdesc_t* dz,
               CFI_cdesc_t* ddx, CFI_cdesc_t* ddy, CFI_cdesc_t* ddz,
               CFI_cdesc_t* xstag, CFI_cdesc_t* ystag, CFI_cdesc_t* zstag)
{
    T x_p, y_p, z_p;
    T dx_p, dy_p, dz_p;
    T ddx_p, ddy_p, ddz_p;
    T xstag_p, ystag_p, zstag_p;
    T res_v;
    T randval_p;

    // Initialize RNG with a non-deterministic random number as seed, define
    // a real distribution that lies in the interval [0.0, 1.0)
    std::random_device rd;
    std::mt19937_64 rng(rd());
    std::uniform_real_distribution<T> dis(0.0, 1.0);

    typedef exprtk::symbol_table<T> symbol_table_t;
    typedef exprtk::expression<T> expression_t;
    typedef exprtk::parser<T> parser_t;

    symbol_table_t symbol_table;
    symbol_table.add_constants();  // pi, epsilon and inf
    symbol_table.add_constant("rho", rho);
    symbol_table.add_constant("gmol", gmol);
    symbol_table.add_constant("tu_level", tu_level);
    symbol_table.add_constant("timeph", timeph);

    symbol_table.add_variable("x", x_p);
    symbol_table.add_variable("y", y_p);
    symbol_table.add_variable("z", z_p);

    symbol_table.add_variable("dx", dx_p);
    symbol_table.add_variable("dy", dy_p);
    symbol_table.add_variable("dz", dz_p);

    symbol_table.add_variable("ddx", ddx_p);
    symbol_table.add_variable("ddy", ddy_p);
    symbol_table.add_variable("ddz", ddz_p);

    symbol_table.add_variable("xstag", xstag_p);
    symbol_table.add_variable("ystag", ystag_p);
    symbol_table.add_variable("zstag", zstag_p);

    symbol_table.add_variable(name, res_v);
    symbol_table.add_variable("rand", randval_p);

    expression_t expression;
    expression.register_symbol_table(symbol_table);

    parser_t parser;
    parser.compile(expr, expression);

    // Sanity checks - check that ranks are as expected
    assert (res->rank == 3);
    assert (x->rank == 1);
    assert (y->rank == 1);
    assert (z->rank == 1);

    assert (dx->rank == 1);
    assert (dy->rank == 1);
    assert (dz->rank == 1);

    assert (ddx->rank == 1);
    assert (ddy->rank == 1);
    assert (ddz->rank == 1);

    assert (xstag->rank == 1);
    assert (ystag->rank == 1);
    assert (zstag->rank == 1);

    // Set grid dimensions
    CFI_index_t ii = x->dim[0].extent;
    CFI_index_t jj = y->dim[0].extent;
    CFI_index_t kk = z->dim[0].extent;

    // Check that all grid dimensions are as expected (boooooring...)
    assert (res->dim[0].extent == kk);
    assert (res->dim[1].extent == jj);
    assert (res->dim[2].extent == ii);

    assert (dx->dim[0].extent == ii);
    assert (dy->dim[0].extent == jj);
    assert (dz->dim[0].extent == kk);

    assert (ddx->dim[0].extent == ii);
    assert (ddy->dim[0].extent == jj);
    assert (ddz->dim[0].extent == kk);

    assert (xstag->dim[0].extent == ii);
    assert (ystag->dim[0].extent == jj);
    assert (zstag->dim[0].extent == kk);

    assert (res->dim[0].lower_bound == 0);
    assert (res->dim[1].lower_bound == 0);
    assert (res->dim[2].lower_bound == 0);

    for (CFI_index_t i = 0; i < ii; i++)
    {
        for (CFI_index_t j = 0; j < jj; j++)
        {
            for (CFI_index_t k = 0; k < kk; k++)
            {
                // The randdom value is updated for every cell and lies in the
                // interval [0, 1.0)
                randval_p = dis(rng);

                CFI_index_t sub_x[1] = {i};
                x_p = *((T *) CFI_address(x, sub_x));
                dx_p = *((T *) CFI_address(dx, sub_x));
                ddx_p = *((T *) CFI_address(ddx, sub_x));
                xstag_p = *((T *) CFI_address(xstag, sub_x));

                CFI_index_t sub_y[1] = {j};
                y_p = *((T *) CFI_address(y, sub_y));
                dy_p = *((T *) CFI_address(dy, sub_y));
                ddy_p = *((T *) CFI_address(ddy, sub_y));
                ystag_p = *((T *) CFI_address(ystag, sub_y));

                CFI_index_t sub_z[1] = {k};
                z_p = *((T *) CFI_address(z, sub_z));
                dz_p = *((T *) CFI_address(dz, sub_z));
                ddz_p = *((T *) CFI_address(ddz, sub_z));
                zstag_p = *((T *) CFI_address(zstag, sub_z));

                CFI_index_t sub3[3] = {k, j, i};
                res_v = *((T *) CFI_address(res, sub3));

                // Evaluate expression and write in result array
                res_v = expression.value();
                T* res_p = (T *) CFI_address(res, sub3);
                *res_p = res_v;
            }
        }
    }
}


extern "C" {
    void eval_real_expr(CFI_cdesc_t* res, const char* name, const char* expr,
        mgletreal rho, mgletreal gmol, mgletreal tu_level, mgletreal timeph,
        CFI_cdesc_t* x, CFI_cdesc_t* y, CFI_cdesc_t* z,
        CFI_cdesc_t* dx, CFI_cdesc_t* dy, CFI_cdesc_t* dz,
        CFI_cdesc_t* ddx, CFI_cdesc_t* ddy, CFI_cdesc_t* ddz,
        CFI_cdesc_t* xstag, CFI_cdesc_t* ystag, CFI_cdesc_t* zstag)
    {
        // Should probably check x, y, z etc. as well...
        assert (res->type == CFI_type_mgletreal);
        eval_expr<mgletreal>(res, name, expr, rho, gmol, tu_level, timeph,
            x, y, z, dx, dy, dz, ddx, ddy, ddz, xstag, ystag, zstag);
    }
}

