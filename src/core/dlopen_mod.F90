MODULE dlopen_mod
    USE,  INTRINSIC :: ISO_C_BINDING, ONLY: C_PTR, C_CHAR, C_INT, C_NULL_PTR, &
        C_NULL_FUNPTR, C_NULL_CHAR, C_FUNPTR, C_ASSOCIATED, C_F_POINTER
    USE precision_mod
    USE fort7_mod
    USE comms_mod
    USE err_mod, ONLY: errr

    IMPLICIT NONE(type, external)
    PRIVATE

    INTEGER(intk), PARAMETER :: shlib_max = 16
    INTEGER(intk), PARAMETER :: shlib_name_maxlen = 64
    INTEGER(intk) :: shlib_idx = 0

    TYPE(C_PTR) :: handle(shlib_max) = C_NULL_PTR

    INTEGER, PARAMETER :: RTLD_LAZY=1, RTLD_NOW=2, RTLD_GLOBAL=256, RTLD_LOCAL=0
    INTERFACE
        ! void * dlopen(const char *filename, int flag);
        FUNCTION dlopen(filename, mode) RESULT(handle) BIND(C, NAME='dlopen')
            IMPORT :: c_char, c_int, c_ptr
            CHARACTER(C_CHAR), DIMENSION(*), INTENT(IN) :: filename
            INTEGER(C_INT), VALUE :: mode
            TYPE(C_PTR) :: handle
        END FUNCTION dlopen

        ! void *dlsym(void *handle, char *symbol);
        FUNCTION dlsym(handle, symbol) RESULT(funptr) BIND(C, NAME='dlsym')
            IMPORT :: c_char, c_funptr, c_ptr
            TYPE(C_PTR), VALUE :: handle
            CHARACTER(C_CHAR), DIMENSION(*), INTENT(IN) :: symbol
            TYPE(C_FUNPTR) :: funptr
        END FUNCTION dlsym

        ! void *dlsym(void *handle, char *symbol);
        FUNCTION dlsym_ptr(handle, symbol) RESULT(ptr) BIND(C, NAME='dlsym')
            ! Returns a C_PTR instead of a C_FUNPTR
            IMPORT :: c_char, c_ptr
            TYPE(C_PTR), VALUE :: handle
            CHARACTER(C_CHAR), DIMENSION(*), INTENT(IN) :: symbol
            TYPE(C_PTR) :: ptr
        END FUNCTION dlsym_ptr

        ! int dlclose(void *handle);
        FUNCTION dlclose(handle) RESULT(ierror) BIND(C, NAME='dlclose')
            IMPORT :: c_int, c_ptr
            TYPE(C_PTR), VALUE :: handle
            INTEGER(C_INT) :: ierror
        END FUNCTION dlclose

        ! char *dlerror(void);
        FUNCTION dlerror() RESULT(error) BIND(C, NAME='dlerror')
            IMPORT :: c_ptr
            TYPE(C_PTR) :: error
        END FUNCTION dlerror
    END INTERFACE

    PUBLIC :: shlib_load, finish_dlopen, shlib_get_fun, shlib_get_ptr, &
        init_dlopen

CONTAINS
    SUBROUTINE init_dlopen()
        INTEGER(intk) :: i, nlibs
        CHARACTER(len=mglet_filename_max) :: soname, jsonptr

        IF (.NOT. fort7%exists("/libs")) THEN
            RETURN
        END IF
        CALL fort7%get_size("/libs", nlibs)
        IF (nlibs <= 0) RETURN

        IF (myid == 0) WRITE(*, '("LOADING LIBRARIES:")')
        DO i = 1, nlibs
            WRITE(jsonptr, '("/libs/", I0)') i-1
            CALL fort7%get_value(jsonptr, soname)
            CALL shlib_load(soname)
        END DO
        IF (myid == 0) WRITE(*, '()')
    END SUBROUTINE init_dlopen


    ! Load shared libraries
    SUBROUTINE shlib_load(soname, required)
        ! Subroutine arguments
        CHARACTER(len=*), INTENT(in) :: soname
        LOGICAL, INTENT(in), OPTIONAL :: required

        ! Local variables
        TYPE(C_PTR) :: c_string
        CHARACTER(len=1024), POINTER :: f_string
        INTEGER(intk) :: str_len
        LOGICAL :: fileexist

        IF (shlib_idx >= shlib_max) THEN
            CALL errr(__FILE__, __LINE__)
        END IF

        ! Check if library exist, fail if required library does not exist
        INQUIRE(file=soname, EXIST=fileexist)
        IF (PRESENT(required)) THEN
            IF (required .AND. .NOT. fileexist) THEN
                CALL errr(__FILE__, __LINE__)
            END IF
        END IF
        IF (.NOT. fileexist) THEN
            RETURN
        END IF

        ! Load library
        handle(shlib_idx + 1) = dlopen(TRIM(soname)//C_NULL_CHAR, &
            INT(IOR(RTLD_GLOBAL, RTLD_LAZY), KIND=C_INT))

        IF (C_ASSOCIATED(handle(shlib_idx + 1))) THEN
            IF (myid == 0) WRITE(*, '(4X, A)') TRIM(soname)
            shlib_idx = shlib_idx + 1
        ELSE
            c_string = dlerror()
            CALL C_F_POINTER(c_string, f_string)
            str_len = MIN(INDEX(f_string, C_NULL_CHAR), 1024)

            WRITE(*,*) 'shlib_load: Error loading ', soname
            WRITE(*,*) 'shlib_load: dlerror said: ', f_string(1:str_len)
            CALL errr(__FILE__, __LINE__)
        END IF
    END SUBROUTINE shlib_load


    SUBROUTINE shlib_get_fun(name, funptr, required)
        ! Subroutine arguments
        CHARACTER(len=*), INTENT(IN) :: name
        TYPE(C_FUNPTR), INTENT(OUT) :: funptr
        LOGICAL, OPTIONAL :: required

        ! Local variables
        INTEGER :: i
        LOGICAL :: success

        success = .FALSE.
        funptr = C_NULL_FUNPTR
        DO i = 1, shlib_idx
            funptr = dlsym(handle(i), trim(name)//C_NULL_CHAR)
            IF (C_ASSOCIATED(funptr)) THEN
                success = .TRUE.
                EXIT
            END IF
        END DO

        IF (.NOT. success .AND. PRESENT(required)) THEN
            IF (required) THEN
                WRITE(*,*) "Could not find: ", name
                CALL errr(__FILE__, __LINE__)
            END IF
        END IF
    END SUBROUTINE shlib_get_fun


    SUBROUTINE shlib_get_ptr(name, ptr, required)
        ! Subroutine arguments
        CHARACTER(len=*), INTENT(IN) :: name
        TYPE(C_PTR), INTENT(OUT) :: ptr
        LOGICAL, OPTIONAL :: required

        ! Local variables
        INTEGER :: i
        LOGICAL :: success

        success = .FALSE.
        ptr = C_NULL_PTR
        DO i = 1, shlib_idx
            ptr = dlsym_ptr(handle(i), trim(name)//C_NULL_CHAR)
            IF (C_ASSOCIATED(ptr)) THEN
                success = .TRUE.
                EXIT
            END IF
        END DO

        IF (.NOT. success .AND. PRESENT(required)) THEN
            IF (required) THEN
                WRITE(*,*) "Could not find: ", name
                CALL errr(__FILE__, __LINE__)
            END IF
        END IF
    END SUBROUTINE shlib_get_ptr


    SUBROUTINE finish_dlopen()
        ! Local variables
        INTEGER :: i
        INTEGER(C_INT) :: ierror

        DO i = 1, shlib_idx
            ierror = dlclose(handle(i))
            IF (ierror /= 0) then
                WRITE(*,*) 'shlib_close: Error closing library idx = ', i
            END IF
        END DO
    END SUBROUTINE finish_dlopen
END MODULE dlopen_mod
