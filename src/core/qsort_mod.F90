# if 1

! New sorting module, using C++ stable_sort instead of qsort_r. This guarantee
! that equal elements keep their original order, which is a great advantage.
!
! This implementation also does not involve creating any temporary arrays
! when sorting
MODULE qsort_mod
    USE, INTRINSIC :: ISO_C_BINDING, ONLY: c_int, c_long_long, c_float, c_double

    USE err_mod, ONLY: errr
    USE precision_mod, ONLY: intk, c_intk, int32, int64, real32, real64

    IMPLICIT NONE(type, external)
    PRIVATE

    INTERFACE

        ! Functions from sort_wrapper.cxx
        SUBROUTINE sort_int(data, idx) BIND(C)
            IMPORT :: c_intk, c_int

            INTEGER(c_int), INTENT(in) :: data(:)
            INTEGER(c_intk), INTENT(out) :: idx(:)
        END SUBROUTINE sort_int

        SUBROUTINE sort_long(data, idx) BIND(C)
            IMPORT :: c_intk, c_long_long

            INTEGER(c_long_long), INTENT(in) :: data(:)
            INTEGER(c_intk), INTENT(out) :: idx(:)
        END SUBROUTINE sort_long

        SUBROUTINE sort_float(data, idx) BIND(C)
            IMPORT :: c_intk, c_float

            REAL(c_float), INTENT(in) :: data(:)
            INTEGER(c_intk), INTENT(out) :: idx(:)
        END SUBROUTINE sort_float

        SUBROUTINE sort_double(data, idx) BIND(C)
            IMPORT :: c_intk, c_double

            ! Subroutine arguments
            REAL(c_double), INTENT(in) :: data(:)
            INTEGER(c_intk), INTENT(out) :: idx(:)
        END SUBROUTINE sort_double
    END INTERFACE

    INTERFACE sortix
        MODULE PROCEDURE :: sortix_int, sortix_long
    END INTERFACE sortix

    INTERFACE sortrx
        MODULE PROCEDURE :: sortrx_sp, sortrx_dp
    END INTERFACE sortrx

    PUBLIC :: sortix, sortrx

CONTAINS

    SUBROUTINE sortix_int(n, data, idx)
        ! Interface identical to sortix from qsort.F in old MGLET repo

        ! Subroutine arguments
        INTEGER(intk), INTENT(in) :: n
        INTEGER(int32), INTENT(in) :: data(:)
        INTEGER(intk), INTENT(out) :: idx(:)

        ! Local variables
        ! none...

        ! Sanity check
        IF (SIZE(data) /= n) CALL errr(__FILE__, __LINE__)
        IF (SIZE(idx) /= n) CALL errr(__FILE__, __LINE__)
        IF (n == 0) CALL errr(__FILE__, __LINE__)

        CALL sort_int(data, idx)

#ifdef _MGLET_DEBUG_
        BLOCK
            INTEGER(intk) :: i
            IF (n > 1) THEN
                DO i = 2, n
                    IF (data(idx(i)) < data(idx(i-1))) THEN
                        CALL errr(__FILE__, __LINE__)
                    END IF
                END DO
            END IF
        END BLOCK
#endif
    END SUBROUTINE sortix_int


    SUBROUTINE sortix_long(n, data, idx)
        ! Interface identical to sortix from qsort.F in old MGLET repo

        ! Subroutine arguments
        INTEGER(intk), INTENT(in) :: n
        INTEGER(int64), INTENT(in) :: data(:)
        INTEGER(intk), INTENT(out) :: idx(:)

        ! Local variables
        ! none...

        ! Sanity check
        IF (SIZE(data) /= n) CALL errr(__FILE__, __LINE__)
        IF (SIZE(idx) /= n) CALL errr(__FILE__, __LINE__)
        IF (n == 0) CALL errr(__FILE__, __LINE__)

        CALL sort_long(data, idx)

#ifdef _MGLET_DEBUG_
        BLOCK
            INTEGER(intk) :: i
            IF (n > 1) THEN
                DO i = 2, n
                    IF (data(idx(i)) < data(idx(i-1))) THEN
                        CALL errr(__FILE__, __LINE__)
                    END IF
                END DO
            END IF
        END BLOCK
#endif
    END SUBROUTINE sortix_long


    SUBROUTINE sortrx_sp(n, data, idx)
        ! Interface identical to sortix from qsort.F in old MGLET repo

        ! Subroutine arguments
        INTEGER(intk), INTENT(in) :: n
        REAL(real32), INTENT(in) :: data(:)
        INTEGER(intk), INTENT(out) :: idx(:)

        ! Local variables
        ! none...

        ! Sanity check
        IF (SIZE(data) /= n) CALL errr(__FILE__, __LINE__)
        IF (SIZE(idx) /= n) CALL errr(__FILE__, __LINE__)
        IF (n == 0) CALL errr(__FILE__, __LINE__)

        CALL sort_float(data, idx)

#ifdef _MGLET_DEBUG_
        BLOCK
            INTEGER(intk) :: i
            IF (n > 1) THEN
                DO i = 2, n
                    IF (data(idx(i)) < data(idx(i-1))) THEN
                        CALL errr(__FILE__, __LINE__)
                    END IF
                END DO
            END IF
        END BLOCK
#endif
    END SUBROUTINE sortrx_sp


    SUBROUTINE sortrx_dp(n, data, idx)
        ! Interface identical to sortix from qsort.F in old MGLET repo

        ! Subroutine arguments
        INTEGER(intk), INTENT(in) :: n
        REAL(real64), INTENT(in) :: data(:)
        INTEGER(intk), INTENT(out) :: idx(:)

        ! Local variables
        ! none...

        ! Sanity check
        IF (SIZE(data) /= n) CALL errr(__FILE__, __LINE__)
        IF (SIZE(idx) /= n) CALL errr(__FILE__, __LINE__)
        IF (n == 0) CALL errr(__FILE__, __LINE__)

        CALL sort_double(data, idx)

#ifdef _MGLET_DEBUG_
        BLOCK
            INTEGER(intk) :: i
            IF (n > 1) THEN
                DO i = 2, n
                    IF (data(idx(i)) < data(idx(i-1))) THEN
                        CALL errr(__FILE__, __LINE__)
                    END IF
                END DO
            END IF
        END BLOCK
#endif
    END SUBROUTINE sortrx_dp
END MODULE qsort_mod

#else

! Old qsort implementation
MODULE qsort_mod
    USE, INTRINSIC :: ISO_C_BINDING, ONLY: C_PTR, C_LOC, C_FUNPTR, C_INT, &
        C_F_POINTER, C_SIZE_T, C_FUNLOC, C_SIZEOF

    USE err_mod, ONLY: errr
    USE precision_mod, ONLY: intk, realk, int_bytes, real_bytes

    IMPLICIT NONE(type, external)
    PRIVATE

    TYPE, BIND(C) :: ctx_t
        INTEGER(c_size_t) :: n
        TYPE(c_ptr) :: ptr
    END TYPE ctx_t

    ! Interfaces to qsort_s and qsort_r from stdlib.h
    INTERFACE
        ! GNU qsort_r extension:
        ! https://man7.org/linux/man-pages/man3/qsort_r.3.html
        !
        ! void qsort_r(void *ptr, size_t count, size_t size,
        !              int (*comp)(const void *, const void *, void *),
        !              void *context)
        !
        ! Warning:
        ! https://sourceware.org/legacy-ml/libc-alpha/2008-12/msg00003.html
        ! Does not work on BSD (different interface)!
        SUBROUTINE qsort_r(ptr, count, size, comp, context) BIND(C)
            USE, INTRINSIC :: ISO_C_BINDING, ONLY: c_ptr, c_int, c_size_t, &
                c_funptr
            TYPE(C_PTR), VALUE :: ptr
            INTEGER(c_size_t), INTENT(in), VALUE :: count
            INTEGER(c_size_t), INTENT(in), VALUE :: size
            TYPE(C_FUNPTR), VALUE :: comp
            TYPE(C_PTR), VALUE :: context
        END SUBROUTINE qsort_r

        ! C11 qsort_s:
        ! https://en.cppreference.com/w/c/algorithm/qsort
        !
        ! errno_t qsort_s(void *ptr, rsize_t count, rsize_t size,
        !                 int (*comp)(const void *, const void *, void *),
        !                 void *context)
        !
        ! Can be used as drop-in replacement if available
        INTEGER(c_int) FUNCTION qsort_s(ptr, count, size, comp, context) BIND(C)
            USE, INTRINSIC :: ISO_C_BINDING, ONLY: c_ptr, c_int, c_size_t, &
                c_funptr
            TYPE(C_PTR), VALUE :: ptr
            INTEGER(c_size_t), INTENT(in), VALUE :: count
            INTEGER(c_size_t), INTENT(in), VALUE :: size
            TYPE(C_FUNPTR), VALUE :: comp
            TYPE(C_PTR), VALUE :: context
        END FUNCTION qsort_s
    END INTERFACE

    PUBLIC :: sortix, sortrx

CONTAINS
    SUBROUTINE sortix(n, data, idx)
        ! Interface identical to sortix from qsort.F in old MGLET repo

        ! Subroutine arguments
        INTEGER(intk), INTENT(in) :: n
        INTEGER(intk), INTENT(in), TARGET :: data(n)
        INTEGER(intk), INTENT(out), TARGET :: idx(n)

        ! Local variables
        TYPE(ctx_t), TARGET :: ctx
        INTEGER(c_size_t) :: count_c
        INTEGER(c_size_t) :: size_c
        INTEGER(intk) :: i

        ! Sanity check
        IF (SIZE(data) /= n) CALL errr(__FILE__, __LINE__)
        IF (SIZE(idx) /= n) CALL errr(__FILE__, __LINE__)
        IF (n == 0) CALL errr(__FILE__, __LINE__)

        ! Context is the original data array, this does not change
        ctx%n = INT(n, c_size_t)
        ctx%ptr = C_LOC(data)

        ! Initialize index array - this is re-organized as the output of the
        ! sorting
        DO i = 1, n
            idx(i) = i
        END DO

        ! Size and count parameters to sorting algorithm
        count_c = INT(n, c_size_t)
        size_c = C_SIZEOF(idx(1))

        ! Sort array
        CALL qsort_r(C_LOC(idx), count_c, size_c, C_FUNLOC(cmp_i), C_LOC(ctx))

#ifdef _MGLET_DEBUG_
        IF (n > 1) THEN
            DO i = 2, n
                IF (data(idx(i)) < data(idx(i-1))) THEN
                    CALL errr(__FILE__, __LINE__)
                END IF
            END DO
        END IF
#endif
    END SUBROUTINE sortix


    INTEGER(c_int) FUNCTION cmp_i(e1, e2, context) BIND(c)
        ! Function arguments
        TYPE(C_PTR), VALUE, INTENT(IN) :: e1
        TYPE(C_PTR), VALUE, INTENT(IN) :: e2
        TYPE(C_PTR), VALUE, INTENT(IN) :: context

        ! Local variables
        INTEGER(intk), POINTER :: i1, i2
        TYPE(ctx_t), POINTER :: ctx
        INTEGER(intk), POINTER, CONTIGUOUS :: data(:)

        CALL C_F_POINTER(e1, i1)
        CALL C_F_POINTER(e2, i2)
        CALL C_F_POINTER(context, ctx)
        CALL C_F_POINTER(ctx%ptr, data, shape=[ctx%n])

        IF (data(i1) < data(i2)) THEN
            cmp_i = -1
        ELSE IF (data(i1) > data(i2)) THEN
            cmp_i = 1
        ELSE
            cmp_i = 0
        END IF
    END FUNCTION cmp_i


    SUBROUTINE sortrx(n, data, idx)
        ! Interface identical to sortix from qsort.F in old MGLET repo

        ! Subroutine arguments
        INTEGER(intk), INTENT(in) :: n
        REAL(realk), INTENT(in), TARGET :: data(n)
        INTEGER(intk), INTENT(out), TARGET :: idx(n)

        ! Local variables
        TYPE(ctx_t), TARGET :: ctx
        INTEGER(c_size_t) :: count_c
        INTEGER(c_size_t) :: size_c
        INTEGER(intk) :: i

        ! Sanity check
        IF (SIZE(data) /= n) CALL errr(__FILE__, __LINE__)
        IF (SIZE(idx) /= n) CALL errr(__FILE__, __LINE__)
        IF (n == 0) CALL errr(__FILE__, __LINE__)

        ! Context is the original data array, this does not change
        ctx%n = INT(n, c_size_t)
        ctx%ptr = C_LOC(data)

        ! Initialize index array - this is re-organized as the output of the
        ! sorting
        DO i = 1, n
            idx(i) = i
        END DO

        ! Size and count parameters to sorting algorithm
        count_c = INT(n, c_size_t)
        size_c = C_SIZEOF(idx(1))

        ! Sort array
        CALL qsort_r(C_LOC(idx), count_c, size_c, C_FUNLOC(cmp_r), C_LOC(ctx))

#ifdef _MGLET_DEBUG_
        IF (n > 1) THEN
            DO i = 2, n
                IF (data(idx(i)) < data(idx(i-1))) THEN
                    CALL errr(__FILE__, __LINE__)
                END IF
            END DO
        END IF
#endif
    END SUBROUTINE sortrx


    INTEGER(c_int) FUNCTION cmp_r(e1, e2, context) BIND(c)
        ! Function arguments
        TYPE(C_PTR), VALUE, INTENT(IN) :: e1
        TYPE(C_PTR), VALUE, INTENT(IN) :: e2
        TYPE(C_PTR), VALUE, INTENT(IN) :: context

        ! Local variables
        INTEGER(intk), POINTER :: i1, i2
        TYPE(ctx_t), POINTER :: ctx
        REAL(realk), POINTER, CONTIGUOUS :: data(:)

        CALL C_F_POINTER(e1, i1)
        CALL C_F_POINTER(e2, i2)
        CALL C_F_POINTER(context, ctx)
        CALL C_F_POINTER(ctx%ptr, data, shape=[ctx%n])

        IF (data(i1) < data(i2)) THEN
            cmp_r = -1
        ELSE IF (data(i1) > data(i2)) THEN
            cmp_r = 1
        ELSE
            cmp_r = 0
        END IF
    END FUNCTION cmp_r
END MODULE qsort_mod
#endif
