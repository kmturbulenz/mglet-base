#include "mglet_precision.h"
#include "sort_wrapper.h"

#include <algorithm>    // std::stable_sort
#include <cassert>      // assert

extern "C" {
#include <ISO_Fortran_binding.h>
}

template <typename T>
void sort_fortran_array(CFI_cdesc_t* idx, CFI_cdesc_t* data)
{
    // Assure that idx and data are of rank 1
    assert (idx->rank == 1);
    assert (data->rank == 1);

    // Assure that they are the same size
    assert (idx->dim[0].extent == data->dim[0].extent);

    // the idx array must be contiguous, the data not!!
    assert (CFI_is_contiguous(idx));

    // Pointer to the *contiguous* idx-array!! This is why idx must
    // be contiguous!!
    mgletint* index = (mgletint*) idx->base_addr;
    CFI_index_t nelems = idx->dim[0].extent;

    // Lower bound of data array
    CFI_index_t lb = data->dim[0].lower_bound;

    // Initialize index array Fortran-style, this is the values that are
    // actually returned back, therefore we must use n+1 and not n here
    for (CFI_index_t n = 0; n < nelems; ++n)
    {
        index[n] = n + 1;
    }

    // Defines a lambda-function 'compare' for use in sorting the array later
    auto compare = [&data, lb](const mgletint &i, const mgletint &j)
    {
        CFI_index_t subi[1] = {lb+i-1};
        T a = *((T *) CFI_address(data, subi));

        CFI_index_t subj[1] = {lb+j-1};
        T b = *((T *) CFI_address(data, subj));

        return (a < b);
    };
    std::stable_sort(index, index+nelems, compare);
}


extern "C" {
    void sort_int(CFI_cdesc_t* data,  CFI_cdesc_t* idx)
    {
        assert (data->type == CFI_type_int || data->type == CFI_type_int32_t);
        assert (idx->type == CFI_type_mgletint);

        sort_fortran_array<int>(idx, data);
    }

    void sort_long(CFI_cdesc_t* data,  CFI_cdesc_t* idx)
    {
        assert (data->type == CFI_type_long_long || data->type == CFI_type_int64_t);
        assert (idx->type == CFI_type_mgletint);

        sort_fortran_array<long long>(idx, data);
    }

    void sort_float(CFI_cdesc_t* data,  CFI_cdesc_t* idx)
    {
        assert (data->type == CFI_type_float);
        assert (idx->type == CFI_type_mgletint);

        sort_fortran_array<float>(idx, data);
    }

    void sort_double(CFI_cdesc_t* data,  CFI_cdesc_t* idx)
    {
        assert (data->type == CFI_type_double);
        assert (idx->type == CFI_type_mgletint);

        sort_fortran_array<double>(idx, data);
    }
}
