MODULE corefields_mod
    USE comms_mod
    USE err_mod, ONLY: errr
    USE field_mod
    USE fieldio2_mod
    USE fields_mod
    USE fort7_mod
    USE grids_mod
    USE hdf5common_mod
    USE pointers_mod
    USE precision_mod

    IMPLICIT NONE(type, external)
    PRIVATE

    ! For convenience in this routine
    TYPE(field_t), POINTER :: rddx, rddy, rddz
    TYPE(field_t), POINTER :: rdx, rdy, rdz
    TYPE(field_t), POINTER :: ddx, ddy, ddz
    TYPE(field_t), POINTER :: dx, dy, dz
    TYPE(field_t), POINTER :: x, y, z
    TYPE(field_t), POINTER :: xstag, ystag, zstag

    PUBLIC :: init_corefields, finish_corefields

CONTAINS
    SUBROUTINE init_corefields()
        ! If we set dread=TRUE in set_field, they will be read in from
        ! fields.h5. We do not want that. Therefore, set no reading
        ! flag on the gridspacing arrays

        CALL set_field("RDDX", ndim=1, get_len=get_ii)
        CALL set_field("RDDY", ndim=1, get_len=get_jj)
        CALL set_field("RDDZ", ndim=1, get_len=get_kk)

        CALL set_field("RDX", ndim=1, get_len=get_ii)
        CALL set_field("RDY", ndim=1, get_len=get_jj)
        CALL set_field("RDZ", ndim=1, get_len=get_kk)

        CALL set_field("DDX", ndim=1, get_len=get_ii)
        CALL set_field("DDY", ndim=1, get_len=get_jj)
        CALL set_field("DDZ", ndim=1, get_len=get_kk)

        CALL set_field("DX", ndim=1, get_len=get_ii)
        CALL set_field("DY", ndim=1, get_len=get_jj)
        CALL set_field("DZ", ndim=1, get_len=get_kk)

        CALL set_field("X", ndim=1, get_len=get_ii)
        CALL set_field("Y", ndim=1, get_len=get_jj)
        CALL set_field("Z", ndim=1, get_len=get_kk)

        CALL set_field("XSTAG", ndim=1, get_len=get_ii)
        CALL set_field("YSTAG", ndim=1, get_len=get_jj)
        CALL set_field("ZSTAG", ndim=1, get_len=get_kk)

        CALL get_field(rddx, "RDDX")
        CALL get_field(rddy, "RDDY")
        CALL get_field(rddz, "RDDZ")

        CALL get_field(rdx, "RDX")
        CALL get_field(rdy, "RDY")
        CALL get_field(rdz, "RDZ")

        CALL get_field(ddx, "DDX")
        CALL get_field(ddy, "DDY")
        CALL get_field(ddz, "DDZ")

        CALL get_field(dx, "DX")
        CALL get_field(dy, "DY")
        CALL get_field(dz, "DZ")

        CALL get_field(x, "X")
        CALL get_field(y, "Y")
        CALL get_field(z, "Z")

        CALL get_field(xstag, "XSTAG")
        CALL get_field(ystag, "YSTAG")
        CALL get_field(zstag, "ZSTAG")

        CALL read_gridspacing()
        CALL calc_reciprocals()
    END SUBROUTINE init_corefields


    SUBROUTINE finish_corefields()

    END SUBROUTINE finish_corefields


    SUBROUTINE read_gridspacing()
        USE HDF5
        USE MPI_f08, ONLY: MPI_Wtime

        ! Local variables
        INTEGER(HID_T) :: file_id
        REAL(real64) :: tic, toc
        CHARACTER(len=mglet_filename_max) :: filename

        filename = REPEAT(" ", mglet_filename_max)
        CALL fort7%get_value("/io/grids", filename)

        tic = MPI_Wtime()

        CALL hdf5common_open(filename, "r", file_id)

        CALL fieldio_read(file_id, x)
        CALL fieldio_read(file_id, dx)
        CALL fieldio_read(file_id, ddx)

        CALL fieldio_read(file_id, y)
        CALL fieldio_read(file_id, dy)
        CALL fieldio_read(file_id, ddy)

        CALL fieldio_read(file_id, z)
        CALL fieldio_read(file_id, dz)
        CALL fieldio_read(file_id, ddz)

        CALL hdf5common_close(file_id)

        toc = MPI_Wtime()
        IF (myid == 0) THEN
            WRITE(*, '("Read gridspacing arrays in ", F6.3, " seconds")') &
                toc - tic
            WRITE(*, '()')
        END IF
    END SUBROUTINE read_gridspacing


#if 0
    SUBROUTINE write_gridspacing(gridid)
        USE HDF5

        ! Subroutine arguments
        INTEGER(intk), INTENT(in) :: gridid(:)

        ! Local variables
        INTEGER(HID_T) :: file_id
        CHARACTER(len=mglet_filename_max) :: filename

        filename = REPEAT(" ", mglet_filename_max)
        CALL fort7%get_value("/io/grids", filename)

        CALL hdf5common_open(filename, "w", file_id)

        CALL fieldio_write(file_id, x)
        CALL fieldio_write(file_id, dx)
        CALL fieldio_write(file_id, ddx)

        CALL fieldio_write(file_id, y)
        CALL fieldio_write(file_id, dy)
        CALL fieldio_write(file_id, ddy)

        CALL fieldio_write(file_id, z)
        CALL fieldio_write(file_id, dz)
        CALL fieldio_write(file_id, ddz)

        CALL hdf5common_close(file_id)
    END SUBROUTINE write_gridspacing
#endif


    SUBROUTINE calc_reciprocals()
        CALL calc_reciprocal(rddx%arr, ddx%arr)
        CALL calc_reciprocal(rddy%arr, ddy%arr)
        CALL calc_reciprocal(rddz%arr, ddz%arr)

        CALL calc_reciprocal(rdx%arr, dx%arr)
        CALL calc_reciprocal(rdy%arr, dy%arr)
        CALL calc_reciprocal(rdz%arr, dz%arr)

        CALL calc_stag(xstag%arr, x%arr, dx%arr)
        CALL calc_stag(ystag%arr, y%arr, dy%arr)
        CALL calc_stag(zstag%arr, z%arr, dz%arr)
    END SUBROUTINE calc_reciprocals


    SUBROUTINE calc_reciprocal(rezip, dx)
        USE simdfunctions_mod, ONLY: divide0

        ! Subroutine arguments
        REAL(realk), INTENT(out) :: rezip(:)
        REAL(realk), INTENT(in) :: dx(:)

        ! Local variables
        INTEGER(intk) :: i, ii

        ii = SIZE(rezip)
        IF (SIZE(dx) /= ii) CALL errr(__FILE__, __LINE__)

        !$omp simd
        DO i = 1, ii
            rezip(i) = divide0(1.0_realk, dx(i))
        END DO
    END SUBROUTINE calc_reciprocal


    SUBROUTINE calc_stag(xstag, x, dx)
        ! Subroutine arguments
        REAL(realk), INTENT(out) :: xstag(:)
        REAL(realk), INTENT(in) :: x(:)
        REAL(realk), INTENT(in) :: dx(:)

        ! Local variables
        INTEGER(intk) :: i, ii

        ii = SIZE(xstag)
        IF (SIZE(x) /= ii) CALL errr(__FILE__, __LINE__)
        IF (SIZE(dx) /= ii) CALL errr(__FILE__, __LINE__)

        !$omp simd
        DO i = 1, ii
            xstag(i) = x(i) + 0.5*dx(i)
        END DO
    END SUBROUTINE calc_stag
END MODULE corefields_mod
