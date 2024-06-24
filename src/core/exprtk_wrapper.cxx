#include "mglet_precision.h"
#include "exprtk_wrapper.h"

#include "exprtk.hpp"

#include <cassert>
#include <cstdio>
#include <random>
#include <string>

#include <ISO_Fortran_binding.h>


template <typename T>
inline T ramp0(T timeph, T time1, T time2) {
    T result;
    if (timeph < time1) {
        result = 0.0;
    }
    else if (timeph < time2) {
        result = (timeph - time1)/(time2 - time1);
    }
    else {
        result = 1.0;
    }

    return result;
}


template <typename T>
struct ramp final : public exprtk::ifunction<T> {
    using exprtk::ifunction<T>::operator();

    ramp() : exprtk::ifunction<T>(3) {
        exprtk::disable_has_side_effects(*this);
    }

    inline T operator()(const T& timeph, const T& time1, const T& time2) {
        // ramp_time is the time it takers to ramp from 0 to 1, that is the
        // time from time1 to time2
        T ramp_time = time2 - time1;

        // Time for transition, this is fixed to 10% of the ramp time
        T trans_time = 0.1*ramp_time;

        // Compute result without smoothing
        T result = ramp0(timeph, time1, time2);

        // Smooth transition from initial flat to slope is only done if
        // time1 is >= than trans_time/2.0
        T ts1 = time1 - trans_time/2.0;
        T ts2 = time1 + trans_time/2.0;
        if (timeph > ts1 && timeph < ts2 && time1 > trans_time/2.0) {
            T delta = (timeph - ts1)/trans_time;
            T end = ramp0(ts2, time1, time2);
            result = end*delta*delta;
        }

        // Smooth transition from slope to flat again - always done
        ts1 = time2 - trans_time/2.0;
        ts2 = time2 + trans_time/2.0;
        if (timeph > ts1 && timeph < ts2) {
            T delta = (timeph - ts1)/trans_time - 1.0;
            T start = 1.0 - ramp0(ts1, time1, time2);
            result = 1.0 - start*delta*delta;
        }

        return result;
    }
};


template <typename T>
struct ramp_inf final : public exprtk::ifunction<T> {
    using exprtk::ifunction<T>::operator();

    ramp_inf() : exprtk::ifunction<T>(2) {
        exprtk::disable_has_side_effects(*this);
    }

    inline T operator()(const T& timeph, const T& time1) {
        // Compute result without smoothing
        T result = timeph < time1 ? 0.0 : (timeph - time1);

        // No smooth transition
        if (time1 <= 0.0) {
            return result;
        }

        // Time for transition, this is fixed to 10% of the flat plateou time
        T trans_time = 0.1*time1;

        // Smooth transition from initial flat to slope
        T ts1 = time1 - trans_time/2.0;
        T ts2 = time1 + trans_time/2.0;
        if (timeph > ts1 && timeph < ts2) {
            T delta = (timeph - ts1)/trans_time;
            T end = ts2 - time1;
            result = end*delta*delta;
        }

        return result;
    }
};


template <typename T>
void eval_expr(CFI_cdesc_t* res, const char* name, const char* expr,
               T rho, T gmol, T tu_level, T timeph,
               CFI_cdesc_t* x, CFI_cdesc_t* y, CFI_cdesc_t* z,
               CFI_cdesc_t* dx, CFI_cdesc_t* dy, CFI_cdesc_t* dz,
               CFI_cdesc_t* ddx, CFI_cdesc_t* ddy, CFI_cdesc_t* ddz,
               int* ierr)
{
    T x_p, y_p, z_p;
    T dx_p, dy_p, dz_p;
    T ddx_p, ddy_p, ddz_p;
    T res_v;
    T randval_p;

    *ierr = 0;

    // Initialize RNG with a non-deterministic random number as seed, define
    // a real distribution that lies in the interval [0.0, 1.0)
    std::random_device rd;
    std::mt19937_64 rng(rd());
    std::uniform_real_distribution<T> dis(0.0, 1.0);

    typedef exprtk::symbol_table<T> symbol_table_t;
    typedef exprtk::expression<T> expression_t;
    typedef exprtk::parser<T> parser_t;
    typedef exprtk::parser_error::type error_t;

    symbol_table_t symbol_table;
    symbol_table.add_constants();  // pi, epsilon and inf
    symbol_table.add_constant("rho", rho);
    symbol_table.add_constant("gmol", gmol);
    symbol_table.add_constant("tu_level", tu_level);
    symbol_table.add_constant("timeph", timeph);

    ramp<T> ramp_func;
    symbol_table.add_function("ramp", ramp_func);

    ramp_inf<T> ramp_inf_func;
    symbol_table.add_function("ramp_inf", ramp_inf_func);

    symbol_table.add_variable("x", x_p);
    symbol_table.add_variable("y", y_p);
    symbol_table.add_variable("z", z_p);

    symbol_table.add_variable("dx", dx_p);
    symbol_table.add_variable("dy", dy_p);
    symbol_table.add_variable("dz", dz_p);

    symbol_table.add_variable("ddx", ddx_p);
    symbol_table.add_variable("ddy", ddy_p);
    symbol_table.add_variable("ddz", ddz_p);

    symbol_table.add_variable(name, res_v);
    symbol_table.add_variable("rand", randval_p);

    // Give possibility to print to terminal (debugging)
    exprtk::rtl::io::package<T> io_package;
    symbol_table.add_package(io_package);

    expression_t expression;
    expression.register_symbol_table(symbol_table);

    parser_t parser;
    if (!parser.compile(expr, expression)) {
        // See exprtk_simple_example_08.cpp for error handling
        printf("Error: %s\tExpression: %s\n", parser.error().c_str(), expr);

        for (std::size_t i = 0; i < parser.error_count(); ++i) {
            error_t error = parser.get_error(i);
            exprtk::parser_error::update_error(error, expr);

            printf("Error: %02d  Line: %04d Col: %02d Type: [%14s] Msg: %s\n",
                static_cast<unsigned int>(i),
                static_cast<unsigned int>(error.line_no),
                static_cast<unsigned int>(error.column_no),
                exprtk::parser_error::to_str(error.mode).c_str(),
                error.diagnostic.c_str());
        }

        *ierr = 1;
        return;
    }

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

    assert (res->dim[0].lower_bound == 0);
    assert (res->dim[1].lower_bound == 0);
    assert (res->dim[2].lower_bound == 0);

    for (CFI_index_t i = 0; i < ii; i++) {
        for (CFI_index_t j = 0; j < jj; j++) {
            for (CFI_index_t k = 0; k < kk; k++) {
                // The random value is updated for every cell and lies in the
                // interval [0, 1.0)
                randval_p = dis(rng);

                CFI_index_t sub_x[1] = {i};
                x_p = *((T *) CFI_address(x, sub_x));
                dx_p = *((T *) CFI_address(dx, sub_x));
                ddx_p = *((T *) CFI_address(ddx, sub_x));

                CFI_index_t sub_y[1] = {j};
                y_p = *((T *) CFI_address(y, sub_y));
                dy_p = *((T *) CFI_address(dy, sub_y));
                ddy_p = *((T *) CFI_address(ddy, sub_y));

                CFI_index_t sub_z[1] = {k};
                z_p = *((T *) CFI_address(z, sub_z));
                dz_p = *((T *) CFI_address(dz, sub_z));
                ddz_p = *((T *) CFI_address(ddz, sub_z));

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
        int* ierr)
    {
        // Should probably check x, y, z etc. as well...
        assert (res->type == CFI_type_mgletreal);
        eval_expr<mgletreal>(res, name, expr, rho, gmol, tu_level, timeph,
            x, y, z, dx, dy, dz, ddx, ddy, ddz, ierr);
    }
}

