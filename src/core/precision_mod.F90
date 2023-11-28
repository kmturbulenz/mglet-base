MODULE precision_mod
    USE, INTRINSIC :: ISO_FORTRAN_ENV
    USE HDF5

    USE MPI_f08

    IMPLICIT NONE (type, external)
    PRIVATE

    ! Kind for REAL's
#ifdef _MGLET_DOUBLE_PRECISION_
    INTEGER(int32), PARAMETER :: realk = real64
    INTEGER(int32), PARAMETER :: c_realk = c_double
    INTEGER(int32), PARAMETER :: real_bytes = 8
#else
    INTEGER(int32), PARAMETER :: realk = real32
    INTEGER(int32), PARAMETER :: c_realk = c_float
    INTEGER(int32), PARAMETER :: real_bytes = 4
#endif

    ! Kind for INTEGERS's
#ifdef _MGLET_INT64_
    INTEGER(int32), PARAMETER :: intk = int64
    INTEGER(int32), PARAMETER :: c_intk = c_long_long
    INTEGER(int32), PARAMETER :: int_bytes = 8
#else
    INTEGER(int32), PARAMETER :: intk = int32
    INTEGER(int32), PARAMETER :: c_intk = c_int
    INTEGER(int32), PARAMETER :: int_bytes = 4
#endif

    ! Special kind for integer fields
#ifdef _MGLET_IFK64_
    INTEGER(int32), PARAMETER :: ifk = int64
    INTEGER(int32), PARAMETER :: ifk_bytes = 8
#else
    INTEGER(int32), PARAMETER :: ifk = int32
    INTEGER(int32), PARAMETER :: ifk_bytes = 4
#endif

    ! MPI data types for REAL and INTEGER
    TYPE(MPI_Datatype), PROTECTED :: mglet_mpi_real
    TYPE(MPI_Datatype), PROTECTED :: mglet_mpi_int
    TYPE(MPI_Datatype), PROTECTED :: mglet_mpi_ifk

    ! HDF5 type for REAL and INTEGER
    INTEGER(HID_T), PROTECTED :: mglet_hdf5_real
    INTEGER(HID_T), PROTECTED :: mglet_hdf5_int
    INTEGER(HID_T), PROTECTED :: mglet_hdf5_ifk

    ! MPI datatype to communicate INTEGER(hsize_t)
    TYPE(MPI_Datatype) :: mglet_mpi_hsize_t

    ! Numerical constants
    REAL(realk), PARAMETER :: pi = 3.14159265358979_realk
    REAL(realk), PARAMETER :: eps = EPSILON(1.0_realk)
#ifdef _MGLET_DOUBLE_PRECISION_
    REAL(realk), PARAMETER :: neps = 1e-10_realk
#else
    REAL(realk), PARAMETER :: neps = 1e-4_realk
#endif


    ! Maximum number of characters in a filename
    INTEGER(intk), PARAMETER :: mglet_filename_max = 256

    ! Public data items
    PUBLIC :: init_precision, finish_precision, realk, intk, c_intk, &
        c_realk, mglet_hdf5_real, mglet_hdf5_int, pi, eps, int8, int16, &
        int32, int64, real32, real64, real_bytes, int_bytes, &
        mglet_filename_max, neps, ifk, ifk_bytes, mglet_mpi_ifk, &
        mglet_hdf5_ifk

    PUBLIC :: mglet_mpi_real, mglet_mpi_int, mglet_mpi_hsize_t

CONTAINS
    SUBROUTINE init_precision()
        USE, INTRINSIC :: ISO_FORTRAN_ENV, ONLY: error_unit

        ! Local variables
        INTEGER(hsize_t), PARAMETER :: dummy_hsize_t = 0
        INTEGER(int32) :: size

        ! Create HDF5 types
        mglet_hdf5_real = h5kind_to_type(realk, H5_REAL_KIND)
        mglet_hdf5_int = h5kind_to_type(intk, H5_INTEGER_KIND)
        mglet_hdf5_ifk = h5kind_to_type(ifk, H5_INTEGER_KIND)

        ! Set mglet_mpi_hsize_t
        size = STORAGE_SIZE(dummy_hsize_t)/8
        IF (size == 8) THEN
            mglet_mpi_hsize_t = MPI_INTEGER8
        ELSE IF (size == 4) THEN
            mglet_mpi_hsize_t = MPI_INTEGER
        ELSE
            WRITE(error_unit, '("precision_mod, unknown size: ", I0)') size
            STOP 255
        END IF

        ! Set MPI data types
#ifdef _MGLET_DOUBLE_PRECISION_
        mglet_mpi_real = MPI_DOUBLE_PRECISION
#else
        mglet_mpi_real = MPI_REAL
#endif

#ifdef _MGLET_INT64_
        mglet_mpi_int = MPI_INTEGER8
#else
        mglet_mpi_int = MPI_INTEGER
#endif

#ifdef _MGLET_IFK64_
        mglet_mpi_ifk = MPI_INTEGER8
#else
        mglet_mpi_ifk = MPI_INTEGER
#endif

    END SUBROUTINE init_precision


    SUBROUTINE finish_precision()
        CONTINUE
    END SUBROUTINE finish_precision
END MODULE precision_mod
