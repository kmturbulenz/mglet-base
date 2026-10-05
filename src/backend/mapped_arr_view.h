#pragma once

#include <ISO_Fortran_binding.h>
#include <omp.h>

#include "errr.h"
#include "mglet_types.h"

namespace mglet::gpu {

template <class T>
class MappedArrView {
  public:
    MappedArrView(CFI_cdesc_t* arr, const char* file = __builtin_FILE(),
                  int line = __builtin_LINE()) {
        if (arr == nullptr) {
            print_location("CFI descriptor is nullptr.", file, line);
            MGLET_ERRR();
        }
        if (arr->base_addr == nullptr) {
            print_location("CFI base_addr is nullptr.", file, line);
            MGLET_ERRR();
        }

        host_ptr_ = static_cast<T*>(arr->base_addr);

        const auto device_num_ = omp_get_default_device();
        device_ptr_ =
            static_cast<T*>(omp_get_mapped_ptr(host_ptr_, device_num_));
        if (device_ptr_ == nullptr) {
            print_location("CFI array is not mapped by OpenMP.", file, line);
            MGLET_ERRR();
        }

        if (arr->rank <= 0) {
            print_location("CFI array does not have rank>=0.", file, line);
            MGLET_ERRR();
        }
        flat_size_ = 1;
        for (std::size_t i = 0; i < arr->rank; ++i) {
            flat_size_ *= static_cast<std::size_t>(arr->dim[i].extent);
        }
    }

    T* host_ptr() const { return host_ptr_; }
    T* device_ptr() const { return device_ptr_; }
    std::size_t flat_size() const { return flat_size_; }

  private:
    T* host_ptr_;
    T* device_ptr_;
    std::size_t flat_size_;
};

} // namespace mglet::gpu
