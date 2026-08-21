MODULE blasbind_mod
    USE, INTRINSIC :: ISO_C_BINDING, ONLY: c_int, c_loc, c_null_ptr, &
        c_ptr, c_size_t, c_sizeof
    USE precision_mod, ONLY: ifk, realk

    IMPLICIT NONE (type, external)
    PRIVATE

    INTERFACE memset
        MODULE PROCEDURE memset_realk
        MODULE PROCEDURE memset_ifk
    END INTERFACE memset

#ifdef _MGLET_OFFLOAD_
    TYPE(c_ptr) :: handle = c_null_ptr
    LOGICAL :: initialized = .FALSE.

    INTERFACE
#ifdef _MGLET_ROCM_
        FUNCTION backend_memset(dst, value, num_bytes) &
                BIND(c, name = 'hipMemset') RESULT(status)
#else
        FUNCTION backend_memset(dst, value, num_bytes) &
                BIND(c, name = 'cudaMemset') RESULT(status)
#endif
            IMPORT :: c_int, c_ptr, c_size_t
            TYPE(c_ptr), VALUE :: dst
            INTEGER(c_int), VALUE :: value
            INTEGER(c_size_t), VALUE :: num_bytes
            INTEGER(c_int) :: status
        END FUNCTION backend_memset

#ifdef _MGLET_ROCM_
        FUNCTION backend_synchronize() BIND(c, name = 'hipDeviceSynchronize') &
                RESULT(status)
#else
        FUNCTION backend_synchronize() &
                BIND(c, name = 'cudaDeviceSynchronize') RESULT(status)
#endif
            IMPORT :: c_int
            INTEGER(c_int) :: status
        END FUNCTION backend_synchronize

#ifdef _MGLET_ROCM_
        FUNCTION create_handle(handle) BIND(c, name = 'rocblas_create_handle') &
                RESULT(status)
#else
        FUNCTION create_handle(handle) BIND(c, name = 'cublasCreate_v2') &
                RESULT(status)
#endif
            IMPORT :: c_int, c_ptr
            TYPE(c_ptr) :: handle
            INTEGER(c_int) :: status
        END FUNCTION create_handle

#ifdef _MGLET_ROCM_
        FUNCTION destroy_handle(handle) &
                BIND(c, name = 'rocblas_destroy_handle') RESULT(status)
#else
        FUNCTION destroy_handle(handle) BIND(c, name = 'cublasDestroy_v2') &
                RESULT(status)
#endif
            IMPORT :: c_int, c_ptr
            TYPE(c_ptr), VALUE :: handle
            INTEGER(c_int) :: status
        END FUNCTION destroy_handle

#ifdef _MGLET_DOUBLE_PRECISION_
#ifdef _MGLET_ROCM_
        FUNCTION backend_axpy(handle, n, alpha, x, incx, y, incy) &
                BIND(c, name = 'rocblas_daxpy') RESULT(status)
#else
        FUNCTION backend_axpy(handle, n, alpha, x, incx, y, incy) &
                BIND(c, name = 'cublasDaxpy_v2') RESULT(status)
#endif
#else
#ifdef _MGLET_ROCM_
        FUNCTION backend_axpy(handle, n, alpha, x, incx, y, incy) &
                BIND(c, name = 'rocblas_saxpy') RESULT(status)
#else
        FUNCTION backend_axpy(handle, n, alpha, x, incx, y, incy) &
                BIND(c, name = 'cublasSaxpy_v2') RESULT(status)
#endif
#endif
            IMPORT :: c_int, c_ptr
            TYPE(c_ptr), VALUE :: handle
            INTEGER(c_int), VALUE :: n
            TYPE(c_ptr), VALUE :: alpha, x, y
            INTEGER(c_int), VALUE :: incx, incy
            INTEGER(c_int) :: status
        END FUNCTION backend_axpy

#ifdef _MGLET_DOUBLE_PRECISION_
#ifdef _MGLET_ROCM_
        FUNCTION backend_scal(handle, n, alpha, x, incx) &
                BIND(c, name = 'rocblas_dscal') RESULT(status)
#else
        FUNCTION backend_scal(handle, n, alpha, x, incx) &
                BIND(c, name = 'cublasDscal_v2') RESULT(status)
#endif
#else
#ifdef _MGLET_ROCM_
        FUNCTION backend_scal(handle, n, alpha, x, incx) &
                BIND(c, name = 'rocblas_sscal') RESULT(status)
#else
        FUNCTION backend_scal(handle, n, alpha, x, incx) &
                BIND(c, name = 'cublasSscal_v2') RESULT(status)
#endif
#endif
            IMPORT :: c_int, c_ptr
            TYPE(c_ptr), VALUE :: handle
            INTEGER(c_int), VALUE :: n
            TYPE(c_ptr), VALUE :: alpha, x
            INTEGER(c_int), VALUE :: incx
            INTEGER(c_int) :: status
        END FUNCTION backend_scal
    END INTERFACE
#endif

    PUBLIC :: init_blasbind, finish_blasbind, axpy, scal, memset, &
            synchronize

CONTAINS
    SUBROUTINE init_blasbind()
#ifdef _MGLET_OFFLOAD_
        INTEGER(c_int) :: status

        IF (initialized) ERROR STOP "BLAS bindings already initialized"
        status = create_handle(handle)
        IF (status /= 0_c_int) ERROR STOP "BLAS handle creation failed"
        initialized = .TRUE.
#endif
    END SUBROUTINE init_blasbind


    SUBROUTINE finish_blasbind()
#ifdef _MGLET_OFFLOAD_
        INTEGER(c_int) :: status

        IF (.NOT. initialized) ERROR STOP "BLAS bindings not initialized"
        status = destroy_handle(handle)
        IF (status /= 0_c_int) ERROR STOP "BLAS handle destruction failed"
        handle = c_null_ptr
        initialized = .FALSE.
#endif
    END SUBROUTINE finish_blasbind


    SUBROUTINE synchronize()
#ifdef _MGLET_OFFLOAD_
        INTEGER(c_int) :: status

        status = backend_synchronize()
        IF (status /= 0_c_int) ERROR STOP "Device synchronization failed"
#endif
    END SUBROUTINE synchronize


    SUBROUTINE memset_realk(arr)
        REAL(realk), CONTIGUOUS, TARGET, INTENT(inout) :: arr(:)

#ifdef _MGLET_OFFLOAD_
        INTEGER(c_int) :: status
        INTEGER(c_size_t) :: num_bytes

        IF (SIZE(arr) == 0) RETURN
        num_bytes = SIZE(arr, kind=c_size_t)*c_sizeof(arr(1))
        !$omp target data use_device_addr(arr)
        status = backend_memset(c_loc(arr), 0_c_int, num_bytes)
        !$omp end target data
        IF (status /= 0_c_int) ERROR STOP "Device memset failed"
#else
        INTEGER :: i

        DO i = 1, SIZE(arr)
            arr(i) = 0.0_realk
        END DO
#endif
    END SUBROUTINE memset_realk


    SUBROUTINE memset_ifk(arr)
        INTEGER(ifk), CONTIGUOUS, TARGET, INTENT(inout) :: arr(:)

#ifdef _MGLET_OFFLOAD_
        INTEGER(c_int) :: status
        INTEGER(c_size_t) :: num_bytes

        IF (SIZE(arr) == 0) RETURN
        num_bytes = SIZE(arr, kind=c_size_t)*c_sizeof(arr(1))
        !$omp target data use_device_addr(arr)
        status = backend_memset(c_loc(arr), 0_c_int, num_bytes)
        !$omp end target data
        IF (status /= 0_c_int) ERROR STOP "Device memset failed"
#else
        INTEGER :: i

        DO i = 1, SIZE(arr)
            arr(i) = 0_ifk
        END DO
#endif
    END SUBROUTINE memset_ifk


    SUBROUTINE axpy(y, alpha, x)
        REAL(realk), CONTIGUOUS, TARGET, INTENT(inout) :: y(:)
        REAL(realk), TARGET, INTENT(in) :: alpha
        REAL(realk), CONTIGUOUS, TARGET, INTENT(in) :: x(:)

#ifdef _MGLET_OFFLOAD_
        INTEGER(c_int) :: status

        IF (SIZE(x) /= SIZE(y)) ERROR STOP "AXPY array sizes differ"
        IF (SIZE(x) == 0) RETURN
        IF (SIZE(x, kind=c_size_t) > HUGE(0_c_int)) &
            ERROR STOP "AXPY array is too large"
        IF (.NOT. initialized) ERROR STOP "BLAS bindings not initialized"

        !$omp target data use_device_addr(x, y)
        status = backend_axpy(handle, INT(SIZE(x), c_int), c_loc(alpha), &
            c_loc(x), 1_c_int, c_loc(y), 1_c_int)
        !$omp end target data
        IF (status /= 0_c_int) ERROR STOP "BLAS AXPY failed"
#else
        INTEGER :: i

        IF (SIZE(x) /= SIZE(y)) ERROR STOP "AXPY array sizes differ"
        DO i = 1, SIZE(x)
            y(i) = y(i) + alpha*x(i)
        END DO
#endif
    END SUBROUTINE axpy


    SUBROUTINE scal(x, alpha)
        REAL(realk), CONTIGUOUS, TARGET, INTENT(inout) :: x(:)
        REAL(realk), TARGET, INTENT(in) :: alpha

#ifdef _MGLET_OFFLOAD_
        INTEGER(c_int) :: status

        IF (SIZE(x) == 0) RETURN
        IF (SIZE(x, kind=c_size_t) > HUGE(0_c_int)) &
                ERROR STOP "SCAL array is too large"
        IF (.NOT. initialized) ERROR STOP "BLAS bindings not initialized"

        !$omp target data use_device_addr(x)
        status = backend_scal(handle, INT(SIZE(x), c_int), &
            c_loc(alpha), c_loc(x), 1_c_int)
        !$omp end target data
        IF (status /= 0_c_int) ERROR STOP "BLAS SCAL failed"
#else
        INTEGER :: i

        DO i = 1, SIZE(x)
            x(i) = alpha*x(i)
        END DO
#endif
    END SUBROUTINE scal
END MODULE blasbind_mod
