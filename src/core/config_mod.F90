MODULE config_mod
    USE, INTRINSIC :: ISO_C_BINDING, ONLY: c_char, c_ptr, c_null_ptr, c_int, &
        c_loc, c_null_char, c_associated, c_size_t, c_int64_t, &
        c_double, c_float, c_bool
    USE, INTRINSIC :: ISO_FORTRAN_ENV, ONLY: int32, int64, real32, real64
    USE err_mod, ONLY: errr

    IMPLICIT NONE(type, external)
    PRIVATE

    TYPE :: config_t
        TYPE(c_ptr), PRIVATE :: handle = C_NULL_PTR
    CONTAINS
        PROCEDURE :: read
        PROCEDURE :: get
        PROCEDURE :: dump
        PROCEDURE :: print

        PROCEDURE, PRIVATE :: get_int32
        PROCEDURE, PRIVATE :: get_int64
        PROCEDURE, PRIVATE :: get_real32
        PROCEDURE, PRIVATE :: get_real64
        PROCEDURE, PRIVATE :: get_logical
        PROCEDURE, PRIVATE :: get_char

        PROCEDURE, PRIVATE :: set_int32
        PROCEDURE, PRIVATE :: set_int64
        PROCEDURE, PRIVATE :: set_real32
        PROCEDURE, PRIVATE :: set_real64
        PROCEDURE, PRIVATE :: set_logical
        PROCEDURE, PRIVATE :: set_char

        PROCEDURE, PRIVATE :: get_int32_arr
        PROCEDURE, PRIVATE :: get_int64_arr
        PROCEDURE, PRIVATE :: get_real32_arr
        PROCEDURE, PRIVATE :: get_real64_arr

        PROCEDURE, PRIVATE :: is_type

        PROCEDURE :: is_int
        PROCEDURE :: is_real
        PROCEDURE :: is_logical
        PROCEDURE :: is_char
        PROCEDURE :: is_object
        PROCEDURE :: is_array

        PROCEDURE :: get_value
        PROCEDURE :: get_array
        PROCEDURE :: get_size

        PROCEDURE :: set_value

        PROCEDURE :: exists
        PROCEDURE :: finish

        FINAL :: destructor
    END TYPE config_t

    ! Interfaces to nlohmann's json library
    INTERFACE
        ! jsoncppc_t* json_from_file(const char *filename)
        TYPE(C_PTR) FUNCTION json_from_file(filename, ierr) BIND(C)
            IMPORT :: c_char, c_ptr, c_int
            CHARACTER(C_CHAR), INTENT(IN) :: filename(*)
            INTEGER(C_INT), INTENT(OUT) :: ierr
        END FUNCTION json_from_file

        ! jsoncppc_t* json_from_json(jsoncppc_t*, const char*, int*)
        TYPE(C_PTR) FUNCTION json_from_json(handle, key, ierr) BIND(C)
            IMPORT :: c_char, c_ptr, c_int
            TYPE(C_PTR), VALUE :: handle
            CHARACTER(C_CHAR), INTENT(IN) :: key(*)
            INTEGER(C_INT), INTENT(OUT) :: ierr
        END FUNCTION json_from_json

        ! void json_free(jsoncppc_t* jsonc, int* ierr)
        SUBROUTINE json_free(handle, ierr) BIND(C)
            IMPORT :: c_ptr, c_int
            TYPE(C_PTR), VALUE :: handle
            INTEGER(C_INT), INTENT(OUT) :: ierr
        END SUBROUTINE json_free

        ! void json_dump(jsoncppc_t* jsonc, CFI_cdesc_t* res, int* ierr)
        SUBROUTINE json_dump(handle, res, ierr) BIND(C)
            IMPORT :: c_ptr, c_int, c_char
            TYPE(C_PTR), VALUE :: handle
#if __GNUC__ < 12
            CHARACTER(len=1, kind=c_char), ALLOCATABLE :: res(:)
#else
            CHARACTER(len=:, kind=c_char), ALLOCATABLE :: res(:)
#endif
            INTEGER(C_INT), INTENT(OUT) :: ierr
        END SUBROUTINE json_dump

        ! void json_get_int(jsoncppc_t* jsonc, const char* key, int* val, int* ierr)
        SUBROUTINE json_get_int(handle, key, val, ierr) BIND(C)
            IMPORT :: c_ptr, c_char, c_int, c_int64_t
            TYPE(C_PTR), VALUE :: handle
            CHARACTER(C_CHAR), INTENT(IN) :: key(*)
            INTEGER(C_INT), INTENT(INOUT) :: val
            INTEGER(C_INT), INTENT(OUT) :: ierr
        END SUBROUTINE json_get_int

        ! void json_set_int(jsoncppc_t* jsonc, const char* key, const int* val, int* ierr)
        SUBROUTINE json_set_int(handle, key, val, ierr) BIND(C)
            IMPORT :: c_ptr, c_char, c_int, c_int64_t
            TYPE(C_PTR), VALUE :: handle
            CHARACTER(C_CHAR), INTENT(IN) :: key(*)
            INTEGER(C_INT), INTENT(IN) :: val
            INTEGER(C_INT), INTENT(OUT) :: ierr
        END SUBROUTINE json_set_int

        ! void json_get_int64(jsoncppc_t* jsonc, const char* key, int64_t* val, int* ierr)
        SUBROUTINE json_get_int64(handle, key, val, ierr) BIND(C)
            IMPORT :: c_ptr, c_char, c_int, c_int64_t
            TYPE(C_PTR), VALUE :: handle
            CHARACTER(C_CHAR), INTENT(IN) :: key(*)
            INTEGER(C_INT64_T), INTENT(INOUT) :: val
            INTEGER(C_INT), INTENT(OUT) :: ierr
        END SUBROUTINE json_get_int64

        ! void json_set_int64(jsoncppc_t* jsonc, const char* key, const int64_t* val, int* ierr)
        SUBROUTINE json_set_int64(handle, key, val, ierr) BIND(C)
            IMPORT :: c_ptr, c_char, c_int, c_int64_t
            TYPE(C_PTR), VALUE :: handle
            CHARACTER(C_CHAR), INTENT(IN) :: key(*)
            INTEGER(C_INT64_T), INTENT(IN) :: val
            INTEGER(C_INT), INTENT(OUT) :: ierr
        END SUBROUTINE json_set_int64

        ! void json_get_float(jsoncppc_t* jsonc, const char* key, float* val, int* ierr)
        SUBROUTINE json_get_float(handle, key, val, ierr) BIND(C)
            IMPORT :: c_ptr, c_char, c_int, c_float
            TYPE(C_PTR), VALUE :: handle
            CHARACTER(C_CHAR), INTENT(IN) :: key(*)
            REAL(C_FLOAT), INTENT(INOUT) :: val
            INTEGER(C_INT), INTENT(OUT) :: ierr
        END SUBROUTINE json_get_float

        ! void json_set_float(jsoncppc_t* jsonc, const char* key, const float* val, int* ierr)
        SUBROUTINE json_set_float(handle, key, val, ierr) BIND(C)
            IMPORT :: c_ptr, c_char, c_int, c_float
            TYPE(C_PTR), VALUE :: handle
            CHARACTER(C_CHAR), INTENT(IN) :: key(*)
            REAL(C_FLOAT), INTENT(IN) :: val
            INTEGER(C_INT), INTENT(OUT) :: ierr
        END SUBROUTINE json_set_float

        ! void json_get_double(jsoncppc_t* jsonc, const char* key, double* val, int* ierr)
        SUBROUTINE json_get_double(handle, key, val, ierr) BIND(C)
            IMPORT :: c_ptr, c_char, c_int, c_double
            TYPE(C_PTR), VALUE :: handle
            CHARACTER(C_CHAR), INTENT(IN) :: key(*)
            REAL(C_DOUBLE), INTENT(INOUT) :: val
            INTEGER(C_INT), INTENT(OUT) :: ierr
        END SUBROUTINE json_get_double

        ! void json_set_double(jsoncppc_t* jsonc, const char* key, const double* val, int* ierr)
        SUBROUTINE json_set_double(handle, key, val, ierr) BIND(C)
            IMPORT :: c_ptr, c_char, c_int, c_double
            TYPE(C_PTR), VALUE :: handle
            CHARACTER(C_CHAR), INTENT(IN) :: key(*)
            REAL(C_DOUBLE), INTENT(IN) :: val
            INTEGER(C_INT), INTENT(OUT) :: ierr
        END SUBROUTINE json_set_double

        ! void json_get_bool(jsoncppc_t* jsonc, const char* key, bool* val, int* ierr)
        SUBROUTINE json_get_bool(handle, key, val, ierr) BIND(C)
            IMPORT :: c_ptr, c_char, c_int, c_bool
            TYPE(C_PTR), VALUE :: handle
            CHARACTER(C_CHAR), INTENT(IN) :: key(*)
            LOGICAL(C_BOOL), INTENT(INOUT) :: val
            INTEGER(C_INT), INTENT(OUT) :: ierr
        END SUBROUTINE json_get_bool

        ! void json_set_bool(jsoncppc_t* jsonc, const char* key, const bool* val, int* ierr)
        SUBROUTINE json_set_bool(handle, key, val, ierr) BIND(C)
            IMPORT :: c_ptr, c_char, c_int, c_bool
            TYPE(C_PTR), VALUE :: handle
            CHARACTER(C_CHAR), INTENT(IN) :: key(*)
            LOGICAL(C_BOOL), INTENT(IN) :: val
            INTEGER(C_INT), INTENT(OUT) :: ierr
        END SUBROUTINE json_set_bool

        ! void json_get_char(jsoncppc_t* jsonc, const char* key, char* val, const size_t maxlen, int* ierr)
        SUBROUTINE json_get_char(handle, key, cval, maxlen, length, ierr) &
                BIND(C)
            IMPORT :: c_ptr, c_char, c_int, c_size_t
            TYPE(C_PTR), VALUE :: handle
            CHARACTER(C_CHAR), INTENT(IN) :: key(*)
            CHARACTER(C_CHAR), INTENT(INOUT) :: cval(*)
            INTEGER(C_SIZE_T), INTENT(IN), VALUE :: maxlen
            INTEGER(C_SIZE_T), INTENT(OUT) :: length
            INTEGER(C_INT), INTENT(OUT) :: ierr
        END SUBROUTINE json_get_char

        ! void json_set_char(jsoncppc_t* jsonc, const char* key, const char* val, int* ierr)
        SUBROUTINE json_set_char(handle, key, cval, ierr) BIND(C)
            IMPORT :: c_ptr, c_char, c_int
            TYPE(C_PTR), VALUE :: handle
            CHARACTER(C_CHAR), INTENT(IN) :: key(*)
            CHARACTER(C_CHAR), INTENT(INOUT) :: cval(*)
            INTEGER(C_INT), INTENT(OUT) :: ierr
        END SUBROUTINE json_set_char

        ! void json_get_int_arr(jsoncppc_t* jsonc, const char* key, int* arr, const size_t length, int* ierr)
        SUBROUTINE json_get_int_arr(handle, key, arr, length, ierr) BIND(C)
            IMPORT :: c_ptr, c_char, c_int, c_float, c_size_t
            TYPE(C_PTR), VALUE :: handle
            CHARACTER(C_CHAR), INTENT(IN) :: key(*)
            !REAL(C_FLOAT), INTENT(INOUT) :: arr(:)
            TYPE(C_PTR), VALUE :: arr
            INTEGER(C_SIZE_T), INTENT(IN), VALUE :: length
            INTEGER(C_INT), INTENT(OUT) :: ierr
        END SUBROUTINE json_get_int_arr

        ! void json_get_int64_arr(jsoncppc_t* jsonc, const char* key, int64_t* arr, const size_t length, int* ierr)
        SUBROUTINE json_get_int64_arr(handle, key, arr, length, ierr) BIND(C)
            IMPORT :: c_ptr, c_char, c_int, c_float, c_size_t
            TYPE(C_PTR), VALUE :: handle
            CHARACTER(C_CHAR), INTENT(IN) :: key(*)
            !REAL(C_INT64_T), INTENT(INOUT) :: arr(:)
            TYPE(C_PTR), VALUE :: arr
            INTEGER(C_SIZE_T), INTENT(IN), VALUE :: length
            INTEGER(C_INT), INTENT(OUT) :: ierr
        END SUBROUTINE json_get_int64_arr

        ! void json_get_double_arr(jsoncppc_t* jsonc, const char* key, float* arr, const size_t length, int* ierr)
        SUBROUTINE json_get_float_arr(handle, key, arr, length, ierr) BIND(C)
            IMPORT :: c_ptr, c_char, c_int, c_float, c_size_t
            TYPE(C_PTR), VALUE :: handle
            CHARACTER(C_CHAR), INTENT(IN) :: key(*)
            !REAL(C_FLOAT), INTENT(INOUT) :: arr(:)
            TYPE(C_PTR), VALUE :: arr
            INTEGER(C_SIZE_T), INTENT(IN), VALUE :: length
            INTEGER(C_INT), INTENT(OUT) :: ierr
        END SUBROUTINE json_get_float_arr

        ! void json_get_double_arr(jsoncppc_t* jsonc, const char* key, double* arr, const size_t length, int* ierr)
        SUBROUTINE json_get_double_arr(handle, key, arr, length, ierr) BIND(C)
            IMPORT :: c_ptr, c_char, c_int, c_double, c_size_t
            TYPE(C_PTR), VALUE :: handle
            CHARACTER(C_CHAR), INTENT(IN) :: key(*)
            !REAL(C_DOUBLE), INTENT(INOUT) :: arr(:)
            TYPE(C_PTR), VALUE :: arr
            INTEGER(C_SIZE_T), INTENT(IN), VALUE :: length
            INTEGER(C_INT), INTENT(OUT) :: ierr
        END SUBROUTINE json_get_double_arr

        ! void json_get_size(jsoncppc_t*, const char*, size_t* size, int* ierr)
        SUBROUTINE json_get_size(handle, key, size, ierr) BIND(C)
            IMPORT :: c_ptr, c_char, c_int, c_size_t
            TYPE(C_PTR), VALUE :: handle
            CHARACTER(C_CHAR), INTENT(IN) :: key(*)
            INTEGER(C_SIZE_T), INTENT(INOUT) :: size
            INTEGER(C_INT), INTENT(OUT) :: ierr
        END SUBROUTINE json_get_size

        ! void json_exists(jsoncppc_t*, const char*, bool*, int*)
        SUBROUTINE json_exists(handle, key, exists, type, ierr) BIND(C)
            IMPORT :: c_ptr, c_char, c_int, c_bool
            TYPE(C_PTR), VALUE :: handle
            CHARACTER(C_CHAR), INTENT(IN) :: key(*)
            LOGICAL(C_BOOL), INTENT(INOUT) :: exists
            INTEGER(C_INT), INTENT(OUT) :: type
            INTEGER(C_INT), INTENT(OUT) :: ierr
        END SUBROUTINE json_exists
    END INTERFACE

    PUBLIC :: config_t

CONTAINS
    SUBROUTINE read(this, filename)
        ! Subroutine arguments
        CLASS(config_t), INTENT(inout) :: this
        CHARACTER(len=*), INTENT(in) :: filename

        ! Local variables
        INTEGER(c_int) :: ierr
        CHARACTER(C_CHAR), DIMENSION(LEN(filename)+1) :: c_filename
        c_filename = TRANSFER(filename, c_filename)
        c_filename(LEN_TRIM(filename)+1) = C_NULL_CHAR

        this%handle = C_NULL_PTR
        this%handle = json_from_file(c_filename, ierr)
        IF (ierr /= 0) CALL errr(__FILE__, __LINE__)
    END SUBROUTINE read


    SUBROUTINE get(this, new, key)
        ! Function arguments
        CLASS(config_t), INTENT(inout) :: this
        CLASS(config_t), INTENT(out) :: new
        CHARACTER(len=*), INTENT(in) :: key

        ! Local variables
        INTEGER(c_int) :: ierr
        CHARACTER(C_CHAR), DIMENSION(LEN(key)+1) :: c_key

        ! Add trailing C_NULL_CHAR to key
        c_key = TRANSFER(key, c_key)
        c_key(LEN_TRIM(key)+1) = C_NULL_CHAR

        new%handle = C_NULL_PTR
        new%handle = json_from_json(this%handle, c_key, ierr)
        IF (ierr /= 0) CALL errr(__FILE__, __LINE__)
    END SUBROUTINE get


    SUBROUTINE dump(this, result)
        ! Function arguments
        CLASS(config_t), INTENT(inout) :: this
#if __GNUC__ < 12
        CHARACTER(kind=C_CHAR, len=1), ALLOCATABLE, INTENT(out) :: result(:)
#else
        CHARACTER(kind=C_CHAR, len=:), ALLOCATABLE, INTENT(out) :: result(:)
#endif

        ! Local variables
        INTEGER(c_int) :: ierr

        ! Get JSON dump
        CALL json_dump(this%handle, result, ierr)
        IF (ierr /= 0) CALL errr(__FILE__, __LINE__)

        ! Sanity check that json_dump did not add a NULL char at end
        IF (result(SIZE(result)) == C_NULL_CHAR) THEN
            CALL errr(__FILE__, __LINE__)
        END IF
    END SUBROUTINE dump


    SUBROUTINE print(this)
        ! Function arguments
        CLASS(config_t), INTENT(inout) :: this

        ! Local variables
#if __GNUC__ < 12
        CHARACTER(kind=C_CHAR, len=1), ALLOCATABLE :: jsondump(:)
#else
        CHARACTER(kind=C_CHAR, len=:), ALLOCATABLE :: jsondump(:)
#endif
        INTEGER(c_int) :: i, ierr

        ! Get JSON dump
        CALL json_dump(this%handle, jsondump, ierr)
        IF (ierr /= 0) CALL errr(__FILE__, __LINE__)

        ! Probably not the most efficient way of writing this to the terminal,
        ! but this is just for debugging purposes...
        DO i = 1, SIZE(jsondump)
            WRITE(*, '(A)', advance="no") jsondump(i)
        END DO
        WRITE(*, '()')

        DEALLOCATE(jsondump)
    END SUBROUTINE print


    SUBROUTINE get_int32(this, key, val, found)
        ! Subroutine arguments
        CLASS(config_t), INTENT(inout) :: this
        CHARACTER(len=*), INTENT(in) :: key
        INTEGER(int32), INTENT(inout) :: val
        LOGICAL, INTENT(out), OPTIONAL :: found

        ! Local variables
        INTEGER(c_int) :: ierr
        CHARACTER(C_CHAR), DIMENSION(LEN(key)+1) :: c_key

        ! Add trailing C_NULL_CHAR to key
        c_key = TRANSFER(key, c_key)
        c_key(LEN_TRIM(key)+1) = C_NULL_CHAR

        IF (PRESENT(found)) found = .FALSE.
        CALL json_get_int(this%handle, c_key, val, ierr)
        IF (ierr > 0) CALL errr(__FILE__, __LINE__)
        IF (ierr < 0) RETURN
        IF (PRESENT(found)) found = .TRUE.
    END SUBROUTINE get_int32


    SUBROUTINE get_int64(this, key, val, found)
        ! Subroutine arguments
        CLASS(config_t), INTENT(inout) :: this
        CHARACTER(len=*), INTENT(in) :: key
        INTEGER(int64), INTENT(inout) :: val
        LOGICAL, INTENT(out), OPTIONAL :: found

        ! Local variables
        INTEGER(c_int) :: ierr
        CHARACTER(C_CHAR), DIMENSION(LEN(key)+1) :: c_key

        ! Add trailing C_NULL_CHAR to key
        c_key = TRANSFER(key, c_key)
        c_key(LEN_TRIM(key)+1) = C_NULL_CHAR

        IF (PRESENT(found)) found = .FALSE.
        CALL json_get_int64(this%handle, c_key, val, ierr)
        IF (ierr > 0) CALL errr(__FILE__, __LINE__)
        IF (ierr < 0) RETURN
        IF (PRESENT(found)) found = .TRUE.
    END SUBROUTINE get_int64


    SUBROUTINE get_real32(this, key, val, found)
        ! Subroutine arguments
        CLASS(config_t), INTENT(inout) :: this
        CHARACTER(len=*), INTENT(in) :: key
        REAL(real32), INTENT(inout) :: val
        LOGICAL, INTENT(out), OPTIONAL :: found

        ! Local variables
        INTEGER(c_int) :: ierr
        CHARACTER(C_CHAR), DIMENSION(LEN(key)+1) :: c_key

        ! Add trailing C_NULL_CHAR to key
        c_key = TRANSFER(key, c_key)
        c_key(LEN_TRIM(key)+1) = C_NULL_CHAR

        IF (PRESENT(found)) found = .FALSE.
        CALL json_get_float(this%handle, c_key, val, ierr)
        IF (ierr > 0) CALL errr(__FILE__, __LINE__)
        IF (ierr < 0) RETURN
        IF (PRESENT(found)) found = .TRUE.
    END SUBROUTINE get_real32


    SUBROUTINE get_real64(this, key, val, found)
        ! Subroutine arguments
        CLASS(config_t), INTENT(inout) :: this
        CHARACTER(len=*), INTENT(in) :: key
        REAL(real64), INTENT(inout) :: val
        LOGICAL, INTENT(out), OPTIONAL :: found

        ! Local variables
        INTEGER(c_int) :: ierr
        CHARACTER(C_CHAR), DIMENSION(LEN(key)+1) :: c_key

        ! Add trailing C_NULL_CHAR to key
        c_key = TRANSFER(key, c_key)
        c_key(LEN_TRIM(key)+1) = C_NULL_CHAR

        IF (PRESENT(found)) found = .FALSE.
        CALL json_get_double(this%handle, c_key, val, ierr)
        IF (ierr > 0) CALL errr(__FILE__, __LINE__)
        IF (ierr < 0) RETURN
        IF (PRESENT(found)) found = .TRUE.
    END SUBROUTINE get_real64


    SUBROUTINE get_logical(this, key, val, found)
        ! Subroutine arguments
        CLASS(config_t), INTENT(inout) :: this
        CHARACTER(len=*), INTENT(in) :: key
        LOGICAL, INTENT(inout) :: val
        LOGICAL, INTENT(out), OPTIONAL :: found

        ! Local variables
        LOGICAL(c_bool) :: c_bool_val
        INTEGER(c_int) :: ierr
        CHARACTER(C_CHAR), DIMENSION(LEN(key)+1) :: c_key

        ! Add trailing C_NULL_CHAR to key
        c_key = TRANSFER(key, c_key)
        c_key(LEN_TRIM(key)+1) = C_NULL_CHAR

        IF (PRESENT(found)) found = .FALSE.
        CALL json_get_bool(this%handle, c_key, c_bool_val, ierr)
        IF (ierr > 0) CALL errr(__FILE__, __LINE__)
        IF (ierr < 0) RETURN
        IF (PRESENT(found)) found = .TRUE.
        val = LOGICAL(c_bool_val)
    END SUBROUTINE get_logical


    SUBROUTINE get_char(this, key, cval, found)
        ! Subroutine arguments
        CLASS(config_t), INTENT(inout) :: this
        CHARACTER(len=*), INTENT(in) :: key
        CHARACTER(len=*), INTENT(inout) :: cval
        LOGICAL, INTENT(out), OPTIONAL :: found

        ! Local variables
        CHARACTER(c_char), DIMENSION(LEN(cval)+1) :: c_cval
        INTEGER(c_size_t) :: maxlen, length, i
        INTEGER(c_int) :: ierr
        CHARACTER(c_char), DIMENSION(LEN(key)+1) :: c_key

        ! Add trailing C_NULL_CHAR to key
        c_key = TRANSFER(key, c_key)
        c_key(LEN_TRIM(key)+1) = C_NULL_CHAR

        maxlen = LEN(cval)
        IF (PRESENT(found)) found = .FALSE.
        CALL json_get_char(this%handle, c_key, c_cval, maxlen, length, ierr)
        IF (ierr > 0) CALL errr(__FILE__, __LINE__)
        IF (ierr < 0) RETURN
        IF (PRESENT(found)) found = .TRUE.

        ! Copy 'length' characters read from file
        DO i = 1, length
            cval(i:i) = c_cval(i)
        END DO
        ! Spacepad rest of string
        DO i = length+1, maxlen
            cval(i:i) = " "
        END DO
    END SUBROUTINE get_char


    SUBROUTINE set_int32(this, key, val)
        ! Subroutine arguments
        CLASS(config_t), INTENT(inout) :: this
        CHARACTER(len=*), INTENT(in) :: key
        INTEGER(int32), INTENT(in) :: val

        ! Local variables
        INTEGER(c_int) :: ierr
        CHARACTER(C_CHAR), DIMENSION(LEN(key)+1) :: c_key

        ! Add trailing C_NULL_CHAR to key
        c_key = TRANSFER(key, c_key)
        c_key(LEN_TRIM(key)+1) = C_NULL_CHAR

        CALL json_set_int(this%handle, c_key, val, ierr)
        IF (ierr /= 0) CALL errr(__FILE__, __LINE__)
    END SUBROUTINE set_int32


    SUBROUTINE set_int64(this, key, val)
        ! Subroutine arguments
        CLASS(config_t), INTENT(inout) :: this
        CHARACTER(len=*), INTENT(in) :: key
        INTEGER(int64), INTENT(in) :: val

        ! Local variables
        INTEGER(c_int) :: ierr
        CHARACTER(C_CHAR), DIMENSION(LEN(key)+1) :: c_key

        ! Add trailing C_NULL_CHAR to key
        c_key = TRANSFER(key, c_key)
        c_key(LEN_TRIM(key)+1) = C_NULL_CHAR

        CALL json_set_int64(this%handle, c_key, val, ierr)
        IF (ierr /= 0) CALL errr(__FILE__, __LINE__)
    END SUBROUTINE set_int64


    SUBROUTINE set_real32(this, key, val)
        ! Subroutine arguments
        CLASS(config_t), INTENT(inout) :: this
        CHARACTER(len=*), INTENT(in) :: key
        REAL(real32), INTENT(in) :: val

        ! Local variables
        INTEGER(c_int) :: ierr
        CHARACTER(C_CHAR), DIMENSION(LEN(key)+1) :: c_key

        ! Add trailing C_NULL_CHAR to key
        c_key = TRANSFER(key, c_key)
        c_key(LEN_TRIM(key)+1) = C_NULL_CHAR

        CALL json_set_float(this%handle, c_key, val, ierr)
        IF (ierr /= 0) CALL errr(__FILE__, __LINE__)
    END SUBROUTINE set_real32


    SUBROUTINE set_real64(this, key, val)
        ! Subroutine arguments
        CLASS(config_t), INTENT(inout) :: this
        CHARACTER(len=*), INTENT(in) :: key
        REAL(real64), INTENT(in) :: val

        ! Local variables
        INTEGER(c_int) :: ierr
        CHARACTER(C_CHAR), DIMENSION(LEN(key)+1) :: c_key

        ! Add trailing C_NULL_CHAR to key
        c_key = TRANSFER(key, c_key)
        c_key(LEN_TRIM(key)+1) = C_NULL_CHAR

        CALL json_set_double(this%handle, c_key, val, ierr)
        IF (ierr /= 0) CALL errr(__FILE__, __LINE__)
    END SUBROUTINE set_real64


    SUBROUTINE set_logical(this, key, val)
        ! Subroutine arguments
        CLASS(config_t), INTENT(inout) :: this
        CHARACTER(len=*), INTENT(in) :: key
        LOGICAL, INTENT(in) :: val

        ! Local variables
        INTEGER(c_int) :: ierr
        CHARACTER(C_CHAR), DIMENSION(LEN(key)+1) :: c_key
        LOGICAL(c_bool) :: c_bool_val

        ! Add trailing C_NULL_CHAR to key
        c_key = TRANSFER(key, c_key)
        c_key(LEN_TRIM(key)+1) = C_NULL_CHAR

        c_bool_val = LOGICAL(val, kind=c_bool)
        CALL json_set_bool(this%handle, c_key, c_bool_val, ierr)
        IF (ierr /= 0) CALL errr(__FILE__, __LINE__)
    END SUBROUTINE set_logical


    SUBROUTINE set_char(this, key, cval)
        ! Subroutine arguments
        CLASS(config_t), INTENT(inout) :: this
        CHARACTER(len=*), INTENT(in) :: key
        CHARACTER(len=*), INTENT(in) :: cval

        ! Local variables
        CHARACTER(c_char), DIMENSION(LEN(cval)+1) :: c_cval
        INTEGER(c_int) :: ierr
        CHARACTER(c_char), DIMENSION(LEN(key)+1) :: c_key

        ! Add trailing C_NULL_CHAR to key
        c_key = TRANSFER(key, c_key)
        c_key(LEN_TRIM(key)+1) = C_NULL_CHAR

        ! Add trailing C_NULL_CHAR to cval
        c_cval = TRANSFER(cval, c_cval)
        c_cval(LEN_TRIM(cval)+1) = C_NULL_CHAR

        CALL json_set_char(this%handle, c_key, c_cval, ierr)
        IF (ierr /= 0) CALL errr(__FILE__, __LINE__)
    END SUBROUTINE set_char


    SUBROUTINE get_int32_arr(this, key, arr, required)
        ! Subroutine arguments
        CLASS(config_t), INTENT(inout) :: this
        CHARACTER(len=*), INTENT(in) :: key
        INTEGER(int32), INTENT(inout), TARGET :: arr(:)
        LOGICAL, INTENT(in), OPTIONAL :: required

        ! Local variables
        INTEGER(c_size_t) :: length
        INTEGER(c_int) :: ierr
        CHARACTER(C_CHAR), DIMENSION(LEN(key)+1) :: c_key

        ! Add trailing C_NULL_CHAR to key
        c_key = TRANSFER(key, c_key)
        c_key(LEN_TRIM(key)+1) = C_NULL_CHAR

        length = SIZE(arr)
        CALL json_get_int_arr(this%handle, c_key, C_LOC(arr), length, ierr)
        IF (ierr /= 0) CALL errr(__FILE__, __LINE__)

        ! TODO: handle required
    END SUBROUTINE get_int32_arr


    SUBROUTINE get_int64_arr(this, key, arr, required)
        ! Subroutine arguments
        CLASS(config_t), INTENT(inout) :: this
        CHARACTER(len=*), INTENT(in) :: key
        INTEGER(int64), INTENT(inout), TARGET :: arr(:)
        LOGICAL, INTENT(in), OPTIONAL :: required

        ! Local variables
        INTEGER(c_size_t) :: length
        INTEGER(c_int) :: ierr
        CHARACTER(C_CHAR), DIMENSION(LEN(key)+1) :: c_key

        ! Add trailing C_NULL_CHAR to key
        c_key = TRANSFER(key, c_key)
        c_key(LEN_TRIM(key)+1) = C_NULL_CHAR

        length = SIZE(arr)
        CALL json_get_int64_arr(this%handle, c_key, C_LOC(arr), length, ierr)
        IF (ierr /= 0) CALL errr(__FILE__, __LINE__)

        ! TODO: handle required
    END SUBROUTINE get_int64_arr


    SUBROUTINE get_real32_arr(this, key, arr, required)
        ! Subroutine arguments
        CLASS(config_t), INTENT(inout) :: this
        CHARACTER(len=*), INTENT(in) :: key
        REAL(real32), INTENT(inout), TARGET :: arr(:)
        LOGICAL, INTENT(in), OPTIONAL :: required

        ! Local variables
        INTEGER(c_size_t) :: length
        INTEGER(c_int) :: ierr
        CHARACTER(C_CHAR), DIMENSION(LEN(key)+1) :: c_key

        ! Add trailing C_NULL_CHAR to key
        c_key = TRANSFER(key, c_key)
        c_key(LEN_TRIM(key)+1) = C_NULL_CHAR

        length = SIZE(arr)
        CALL json_get_float_arr(this%handle, c_key, C_LOC(arr), length, ierr)
        IF (ierr /= 0) CALL errr(__FILE__, __LINE__)

        ! TODO: handle required
    END SUBROUTINE get_real32_arr


    SUBROUTINE get_real64_arr(this, key, arr, required)
        ! Subroutine arguments
        CLASS(config_t), INTENT(inout) :: this
        CHARACTER(len=*), INTENT(in) :: key
        REAL(real64), INTENT(inout), TARGET :: arr(:)
        LOGICAL, INTENT(in), OPTIONAL :: required

        ! Local variables
        INTEGER(c_size_t) :: length
        INTEGER(c_int) :: ierr
        CHARACTER(C_CHAR), DIMENSION(LEN(key)+1) :: c_key

        ! Add trailing C_NULL_CHAR to key
        c_key = TRANSFER(key, c_key)
        c_key(LEN_TRIM(key)+1) = C_NULL_CHAR

        length = SIZE(arr)
        CALL json_get_double_arr(this%handle, c_key, C_LOC(arr), length, ierr)
        IF (ierr /= 0) CALL errr(__FILE__, __LINE__)

        ! TODO: handle required
    END SUBROUTINE get_real64_arr


    LOGICAL FUNCTION is_type(this, key, type)
        ! Subroutine arguments
        CLASS(config_t), INTENT(inout) :: this
        CHARACTER(len=*), INTENT(in) :: key
        INTEGER(c_int), INTENT(in) :: type

        ! Local variables
        INTEGER(c_int) :: ierr, actual_type
        LOGICAL(c_bool) :: c_exists
        CHARACTER(C_CHAR), DIMENSION(LEN(key)+1) :: c_key

        ! Add trailing C_NULL_CHAR to key
        c_key = TRANSFER(key, c_key)
        c_key(LEN_TRIM(key)+1) = C_NULL_CHAR

        c_exists = .FALSE.
        is_type = .FALSE.
        actual_type = -1
        CALL json_exists(this%handle, c_key, c_exists, actual_type, ierr)
        IF (ierr /= 0) CALL errr(__FILE__, __LINE__)

        ! https://github.com/nlohmann/json/blob/develop/include/nlohmann/detail/value_t.hpp
        IF (c_exists .AND. actual_type == type) is_type = .TRUE.
    END FUNCTION is_type


    LOGICAL FUNCTION is_int(this, key)
        ! Subroutine arguments
        CLASS(config_t), INTENT(inout) :: this
        CHARACTER(len=*), INTENT(in) :: key

        ! https://github.com/nlohmann/json/blob/develop/include/nlohmann/detail/value_t.hpp
        is_int = this%is_type(key, 6)
    END FUNCTION is_int


    LOGICAL FUNCTION is_real(this, key)
        ! Subroutine arguments
        CLASS(config_t), INTENT(inout) :: this
        CHARACTER(len=*), INTENT(in) :: key

        ! https://github.com/nlohmann/json/blob/develop/include/nlohmann/detail/value_t.hpp
        is_real = this%is_type(key, 7)
    END FUNCTION is_real


    LOGICAL FUNCTION is_logical(this, key)
        ! Subroutine arguments
        CLASS(config_t), INTENT(inout) :: this
        CHARACTER(len=*), INTENT(in) :: key

        ! https://github.com/nlohmann/json/blob/develop/include/nlohmann/detail/value_t.hpp
        is_logical = this%is_type(key, 4)
    END FUNCTION is_logical


    LOGICAL FUNCTION is_char(this, key)
        ! Subroutine arguments
        CLASS(config_t), INTENT(inout) :: this
        CHARACTER(len=*), INTENT(in) :: key

        ! https://github.com/nlohmann/json/blob/develop/include/nlohmann/detail/value_t.hpp
        is_char = this%is_type(key, 3)
    END FUNCTION is_char


    LOGICAL FUNCTION is_object(this, key)
        ! Subroutine arguments
        CLASS(config_t), INTENT(inout) :: this
        CHARACTER(len=*), INTENT(in) :: key

        ! https://github.com/nlohmann/json/blob/develop/include/nlohmann/detail/value_t.hpp
        is_object = this%is_type(key, 1)
    END FUNCTION is_object


    LOGICAL FUNCTION is_array(this, key)
        ! Subroutine arguments
        CLASS(config_t), INTENT(inout) :: this
        CHARACTER(len=*), INTENT(in) :: key

        ! https://github.com/nlohmann/json/blob/develop/include/nlohmann/detail/value_t.hpp
        is_array = this%is_type(key, 2)
    END FUNCTION is_array


    SUBROUTINE get_value(this, key, value, default_value, found)
        ! Subroutine arguments
        CLASS(config_t), INTENT(inout) :: this
        CHARACTER(len=*), INTENT(in) :: key
        CLASS(*), INTENT(inout) :: value
        CLASS(*), INTENT(in), OPTIONAL :: default_value
        LOGICAL, INTENT(out), OPTIONAL :: found

        LOGICAL :: is_found
        is_found = .FALSE.

        SELECT TYPE (value)
        TYPE IS (INTEGER(int32))
            CALL this%get_int32(key, value, is_found)
            IF (.NOT. is_found) THEN
                IF (PRESENT(default_value)) THEN
                    SELECT TYPE (default_value)
                    TYPE IS (INTEGER(int32))
                        value = default_value
                        CALL this%set_int32(key, value)
                    CLASS DEFAULT
                        CALL errr(__FILE__, __LINE__)
                    END SELECT
                ELSE
                    WRITE(*, *) "Key not found: ", key
                    CALL errr(__FILE__, __LINE__)
                END IF
            END IF
        TYPE IS (INTEGER(int64))
            CALL this%get_int64(key, value, is_found)
            IF (.NOT. is_found) THEN
                IF (PRESENT(default_value)) THEN
                    SELECT TYPE (default_value)
                    TYPE IS (INTEGER(int64))
                        value = default_value
                        CALL this%set_int64(key, value)
                    TYPE IS (INTEGER(int32))
                        value = INT(default_value, int32)
                        CALL this%set_int64(key, value)
                    CLASS DEFAULT
                        CALL errr(__FILE__, __LINE__)
                    END SELECT
                ELSE
                    WRITE(*, *) "Key not found: ", key
                    CALL errr(__FILE__, __LINE__)
                END IF
            END IF
        TYPE IS (REAL(real32))
            CALL this%get_real32(key, value, is_found)
            IF (.NOT. is_found) THEN
                IF (PRESENT(default_value)) THEN
                    SELECT TYPE (default_value)
                    TYPE IS (REAL(real32))
                        value = default_value
                        CALL this%set_real32(key, value)
                    CLASS DEFAULT
                        CALL errr(__FILE__, __LINE__)
                    END SELECT
                ELSE
                    WRITE(*, *) "Key not found: ", key
                    CALL errr(__FILE__, __LINE__)
                END IF
            END IF
        TYPE IS (REAL(real64))
            CALL this%get_real64(key, value, is_found)
            IF (.NOT. is_found) THEN
                IF (PRESENT(default_value)) THEN
                    SELECT TYPE (default_value)
                    TYPE IS (REAL(real64))
                        value = default_value
                        CALL this%set_real64(key, value)
                    TYPE IS (REAL(real32))
                        value = REAL(default_value, real64)
                        CALL this%set_real64(key, value)
                    CLASS DEFAULT
                        CALL errr(__FILE__, __LINE__)
                    END SELECT
                ELSE
                    WRITE(*, *) "Key not found: ", key
                    CALL errr(__FILE__, __LINE__)
                END IF
            END IF
        TYPE IS (LOGICAL)
            CALL this%get_logical(key, value, is_found)
            IF (.NOT. is_found) THEN
                IF (PRESENT(default_value)) THEN
                    SELECT TYPE (default_value)
                    TYPE IS (LOGICAL)
                        value = default_value
                        CALL this%set_logical(key, value)
                    CLASS DEFAULT
                        CALL errr(__FILE__, __LINE__)
                    END SELECT
                ELSE
                    WRITE(*, *) "Key not found: ", key
                    CALL errr(__FILE__, __LINE__)
                END IF
            END IF
        TYPE IS (CHARACTER(len=*))
            CALL this%get_char(key, value, is_found)
            IF (.NOT. is_found) THEN
                IF (PRESENT(default_value)) THEN
                    SELECT TYPE (default_value)
                    TYPE IS (CHARACTER(len=*))
                        value = default_value
                        CALL this%set_char(key, value)
                    CLASS DEFAULT
                        CALL errr(__FILE__, __LINE__)
                    END SELECT
                ELSE
                    WRITE(*, *) "Key not found: ", key
                    CALL errr(__FILE__, __LINE__)
                END IF
            END IF
        CLASS DEFAULT
            CALL errr(__FILE__, __LINE__)
        END SELECT

        IF (PRESENT(found)) found = is_found
    END SUBROUTINE get_value


    SUBROUTINE get_array(this, key, array, required)
        ! Subroutine arguments
        CLASS(config_t), INTENT(inout) :: this
        CHARACTER(len=*), INTENT(in) :: key
        CLASS(*), INTENT(inout) :: array(:)
        LOGICAL, INTENT(in), OPTIONAL :: required

        SELECT TYPE (array)
        TYPE IS (INTEGER(int32))
            CALL this%get_int32_arr(key, array, required)
        TYPE IS (INTEGER(int64))
            CALL this%get_int64_arr(key, array, required)
        TYPE IS (REAL(real32))
            CALL this%get_real32_arr(key, array, required)
        TYPE IS (REAL(real64))
            CALL this%get_real64_arr(key, array, required)
        CLASS DEFAULT
            CALL errr(__FILE__, __LINE__)
        END SELECT
    END SUBROUTINE get_array


    SUBROUTINE get_size(this, key, size)
        ! Subroutine arguments
        CLASS(config_t), INTENT(inout) :: this
        CHARACTER(len=*), INTENT(in) :: key
        CLASS(*), INTENT(out) :: size

        ! Local variables
        INTEGER(c_size_t) :: c_size
        INTEGER(c_int) :: ierr
        CHARACTER(C_CHAR), DIMENSION(LEN(key)+1) :: c_key

        ! Add trailing C_NULL_CHAR to key
        c_key = TRANSFER(key, c_key)
        c_key(LEN_TRIM(key)+1) = C_NULL_CHAR

        CALL json_get_size(this%handle, c_key, c_size, ierr)
        IF (ierr /= 0) CALL errr(__FILE__, __LINE__)

        SELECT TYPE (size)
        TYPE IS (INTEGER(int32))
            IF (c_size > HUGE(1_int32)) THEN
                CALL errr(__FILE__, __LINE__)
            END IF
            size = INT(c_size, int32)
        TYPE IS (INTEGER(int64))
            size = INT(c_size, int64)
        CLASS DEFAULT
            CALL errr(__FILE__, __LINE__)
        END SELECT
    END SUBROUTINE get_size


    SUBROUTINE set_value(this, key, value)
        ! Subroutine arguments
        CLASS(config_t), INTENT(inout) :: this
        CHARACTER(len=*), INTENT(in) :: key
        CLASS(*), INTENT(in) :: value

        SELECT TYPE (value)
        TYPE IS (INTEGER(int32))
            CALL this%set_int32(key, value)
        TYPE IS (INTEGER(int64))
            CALL this%set_int64(key, value)
        TYPE IS (REAL(real32))
            CALL this%set_real32(key, value)
        TYPE IS (REAL(real64))
            CALL this%set_real64(key, value)
        TYPE IS (LOGICAL)
            CALL this%set_logical(key, value)
        TYPE IS (CHARACTER(len=*))
            CALL this%set_char(key, value)
        CLASS DEFAULT
            CALL errr(__FILE__, __LINE__)
        END SELECT
    END SUBROUTINE set_value


    LOGICAL FUNCTION exists(this, key)
        ! Subroutine arguments
        CLASS(config_t), INTENT(inout) :: this
        CHARACTER(len=*), INTENT(in) :: key

        ! Local variables
        INTEGER(c_int) :: type, ierr
        LOGICAL(c_bool) :: c_exists
        CHARACTER(C_CHAR), DIMENSION(LEN(key)+1) :: c_key

        ! Add trailing C_NULL_CHAR to key
        c_key = TRANSFER(key, c_key)
        c_key(LEN_TRIM(key)+1) = C_NULL_CHAR

        exists = .FALSE.
        c_exists = .FALSE.
        CALL json_exists(this%handle, c_key, c_exists, type, ierr)
        IF (ierr /= 0) CALL errr(__FILE__, __LINE__)

        exists = LOGICAL(c_exists)
    END FUNCTION exists


    SUBROUTINE finish(this)
        ! Subroutine arguments
        CLASS(config_t), INTENT(inout) :: this

        ! Local variables
        INTEGER(c_int) :: ierr

        IF (.NOT. C_ASSOCIATED(this%handle)) RETURN
        CALL json_free(this%handle, ierr)
        IF (ierr /= 0) CALL errr(__FILE__, __LINE__)
        this%handle = C_NULL_PTR
    END SUBROUTINE finish


    SUBROUTINE destructor(this)
        ! Subroutine arguments
        TYPE(config_t), INTENT(inout) :: this
        CALL this%finish()
    END SUBROUTINE destructor
END MODULE config_mod
