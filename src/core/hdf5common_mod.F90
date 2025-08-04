MODULE hdf5common_mod
    USE HDF5
    USE MPI_f08
    USE precision_mod, ONLY: int32, intk, realk, c_intk, c_realk, &
        mglet_hdf5_int, mglet_hdf5_real, mglet_mpi_hsize_t, &
        mglet_mpi_real, mglet_mpi_int
    USE err_mod, ONLY: errr
    USE comms_mod, ONLY: iocomm, ioproc, myid, iogrcomm

    IMPLICIT NONE (type, external)
    PRIVATE

    ! Maximum number of dimensions for HDF5 dataset
    ! HDF5 library permits up to 32...
    INTEGER(HSIZE_T), PARAMETER :: maxrank = 8

    ! H5FD_MPIO_INDEPENDENT_F is a default-kind integer and has the value 0.
    ! it is not a compile-time constant and cannot be used in an
    ! initialization expression. Therefore setting 0 here.
    INTEGER :: hdf5_io_mode = 0

    ! Compound HDF5 type for time storage in various sampling
    ! features, such as pload, probes, afsi etc.
    TYPE, BIND(C) :: timeinfo_t
        INTEGER(c_intk) :: ittot
        REAL(c_realk) :: time
    END TYPE timeinfo_t
    INTEGER(hid_t) :: timeinfo_h5t

    PUBLIC :: init_hdf5common, finish_hdf5common, &
        hdf5common_open, hdf5common_close, &
        hdf5common_group_open, hdf5common_group_close, &
        hdf5common_dataset_create, hdf5common_dataset_open, &
        hdf5common_dataset_close, hdf5common_dataset_extend, &
        hdf5common_dataset_exists, &
        hdf5common_attr_write, hdf5common_attr_write_arr, &
        hdf5common_attr_exists, &
        hdf5common_attr_read, hdf5common_attr_read_arr, &
        hdf5common_attr_read_collective, hdf5common_attr_read_arr_collective, &
        hdf5_io_mode, timeinfo_h5t, timeinfo_t

CONTAINS
    SUBROUTINE init_hdf5common()
        USE charfunc_mod, ONLY: upper
        USE envvars_mod, ONLY: getenv_char

        ! Local variables
        CHARACTER(len=32) :: iomodestr
        LOGICAL :: is_collective

        ! Default value if no environment variable is set
        is_collective = .TRUE.

        ! There is a Bcast on this later so I don't use _coll here...
        CALL getenv_char(iomodestr, "MGLET_IO_MODE", "")
        IF (upper(TRIM(iomodestr)) == "INDEPENDENT") THEN
            is_collective = .FALSE.
        END IF

        CALL MPI_Bcast(is_collective, 1, MPI_LOGICAL, 0, MPI_COMM_WORLD)

        IF (myid == 0) THEN
            WRITE(*, '("HDF5 collective IO: ", L1)') is_collective
            WRITE(*, '()')
        END IF

        IF (is_collective) THEN
            hdf5_io_mode = H5FD_MPIO_COLLECTIVE_F
        ELSE
            hdf5_io_mode = H5FD_MPIO_INDEPENDENT_F
        END IF

        ! Create some useful compound types
        CALL create_timeinfo_h5t()
    END SUBROUTINE init_hdf5common


    SUBROUTINE create_timeinfo_h5t()
        ! Local variables
        TYPE(timeinfo_t), TARGET :: foo
        INTEGER(int32) :: hdferr

        CALL h5tcreate_f(H5T_COMPOUND_F, C_SIZEOF(foo), timeinfo_h5t, hdferr)
        IF (hdferr /= 0) CALL errr(__FILE__, __LINE__)

        CALL h5tinsert_f(timeinfo_h5t, "ITTOT", &
            H5OFFSETOF(C_LOC(foo), C_LOC(foo%ittot)), mglet_hdf5_int, hdferr)
        IF (hdferr /= 0) CALL errr(__FILE__, __LINE__)

        CALL h5tinsert_f(timeinfo_h5t, "TIME", &
            H5OFFSETOF(C_LOC(foo), C_LOC(foo%time)), mglet_hdf5_real, hdferr)
        IF (hdferr /= 0) CALL errr(__FILE__, __LINE__)
    END SUBROUTINE create_timeinfo_h5t


    SUBROUTINE finish_hdf5common()
        ! Local variables
        INTEGER(int32) :: hdferr

        CALL h5tclose_f(timeinfo_h5t, hdferr)
        IF (hdferr < 0) CALL errr(__FILE__, __LINE__)
    END SUBROUTINE finish_hdf5common


    SUBROUTINE hdf5common_open(filename, mode, file_id)
        ! Subroutine arguments
        CHARACTER(LEN=*), INTENT(IN) :: filename
        CHARACTER(LEN=*), INTENT(IN) :: mode
        INTEGER(HID_T), INTENT(OUT) :: file_id

        ! Local variables
        LOGICAL :: status, fileexist
        INTEGER(kind=int32) :: ierr

        INTEGER(HID_T) :: plist_id

        file_id = -1
        IF (.NOT. ioProc) THEN
            file_id = -1
            RETURN
        END IF

        IF (LEN(mode) < 1 .OR. LEN(mode) > 2) THEN
            WRITE(*, *) "Invalid mode: ", mode
            CALL errr(__FILE__, __LINE__)
        END IF

        ! setup file access property list with parallel I/O access.
        CALL h5pcreate_f(H5P_FILE_ACCESS_F, plist_id, ierr)
        IF (ierr < 0) CALL errr(__FILE__, __LINE__)

        ! this is to return error if objects are still open upon file closure
        CALL h5pset_fclose_degree_f(plist_id, H5F_CLOSE_SEMI_F, ierr)
        IF (ierr < 0) CALL errr(__FILE__, __LINE__)

        ! We want to write 1.10 compatible files as of now
        CALL h5pset_libver_bounds_f(plist_id, H5F_LIBVER_V110_F, &
            H5F_LIBVER_V110_F, ierr)
        IF (ierr < 0) CALL errr(__FILE__, __LINE__)

        ! TODO: Check if this give advantages - according to documentation
        ! it give performance benefits when the amount of metadata is large
        ! this is maybe not the case for MGLET?
        !
        ! Set collective metadata write (require HDF5 ver. > 1.10)
        CALL h5pset_coll_metadata_write_f(plist_id, .TRUE., ierr)
        IF (ierr < 0) CALL errr(__FILE__, __LINE__)

        CALL h5pset_fapl_mpio_f(plist_id, iocomm%MPI_val, &
            MPI_INFO_NULL%MPI_val, ierr)
        IF (ierr /= 0) CALL errr(__FILE__, __LINE__)

        ! Determine if the in/output file already exists
        INQUIRE(file=filename, EXIST=fileexist)
        IF (fileexist) THEN
            CALL h5fis_hdf5_f(filename, status, ierr)

            IF (TRIM(mode) == 'r' .OR. TRIM(mode) == 'r+' .OR. &
                    TRIM(mode) == 'a') THEN
                IF (ierr < 0 .OR. .NOT. status) THEN
                    WRITE(*, *) "Opening file ", TRIM(filename), &
                        " with mode ", mode
                    WRITE(*, *) "Not a valid HDF5 file."
                    WRITE(*, *) "    ierr: ", ierr
                    WRITE(*, *) "  status: ", status
                    CALL errr(__FILE__, __LINE__)
                END IF
            ELSEIF (TRIM(mode) == "wx") THEN
                WRITE(*, *) "Opening file ", TRIM(filename), &
                    " with mode ", mode
                WRITE(*, *) "File already exist."
                CALL errr(__FILE__, __LINE__)
            END IF
        END IF

        ! Read-only mode, file must exist
        IF (TRIM(mode) == "r") THEN
            CALL h5fopen_f(filename, H5F_ACC_RDONLY_F, file_id, ierr, &
                access_prp=plist_id)
            IF (ierr < 0) CALL errr(__FILE__, __LINE__)

        ! Read-write mode, file must exist
        ELSEIF (TRIM(mode) == "r+") THEN
            CALL h5fopen_f(filename, H5F_ACC_RDWR_F, file_id, ierr, &
                access_prp=plist_id)
            IF (ierr < 0) CALL errr(__FILE__, __LINE__)

        ! Open for reading/writing - if file does not exist, create it,
        ! otherwise open file without altering contents
        ELSE IF (TRIM(mode) == "a") THEN
            IF (fileexist) THEN
                CALL h5fopen_f(filename, H5F_ACC_RDWR_F, file_id, ierr, &
                    access_prp=plist_id)
                IF (ierr < 0) CALL errr(__FILE__, __LINE__)
            ELSE
                CALL h5fcreate_f(filename, H5F_ACC_EXCL_F, file_id, ierr, &
                    access_prp=plist_id)
                IF (ierr < 0) CALL errr(__FILE__, __LINE__)
            END IF

        ! Write mode, truncate file if existing
        ELSE IF (TRIM(mode) == "w") THEN
            CALL h5fcreate_f(filename, H5F_ACC_TRUNC_F, file_id, ierr, &
                access_prp=plist_id)
            IF (ierr < 0) CALL errr(__FILE__, __LINE__)

        ! Write mode, fail if existing
        ELSE IF (TRIM(mode) == "wx") THEN
            CALL h5fcreate_f(filename, H5F_ACC_EXCL_F, file_id, ierr, &
                access_prp=plist_id)
            IF (ierr < 0) CALL errr(__FILE__, __LINE__)

        ELSE
            WRITE(*, *) "Invalid mode: ", mode
            CALL errr(__FILE__, __LINE__)
        END IF

        ! close property list
        CALL h5pclose_f(plist_id, ierr)
        IF (ierr < 0) CALL errr(__FILE__, __LINE__)
    END SUBROUTINE hdf5common_open


    SUBROUTINE hdf5common_close(file_id)
        ! Subroutine arguments
        INTEGER(HID_T), INTENT(IN) :: file_id

        ! Local variables
        INTEGER(kind=int32) :: ierr

        IF (.NOT. ioProc) RETURN

        CALL h5fclose_f(file_id, ierr)
        IF (ierr < 0) CALL errr(__FILE__, __LINE__)
    END SUBROUTINE hdf5common_close


    SUBROUTINE hdf5common_group_open(name, parent_id, group_id, track_index)
        ! Subroutine arguments
        CHARACTER(LEN=*), INTENT(IN) :: name
        INTEGER(HID_T), INTENT(IN) :: parent_id
        INTEGER(HID_T), INTENT(OUT) :: group_id
        LOGICAL, INTENT(IN), OPTIONAL :: track_index

        ! Local variables
        INTEGER(HID_T) :: gcpl_id
        INTEGER(kind=int32) :: ierr
        LOGICAL :: link_exists
        LOGICAL :: track

        IF (.NOT. ioProc) THEN
            group_id = -1
            RETURN
        END IF

        track = .FALSE.
        IF (PRESENT(track_index)) THEN
            IF (track_index) THEN
                track = .TRUE.
            END IF
        END IF

        CALL h5lexists_f(parent_id, name, link_exists, ierr)
        IF (ierr /= 0) CALL errr(__FILE__, __LINE__)

        IF (.NOT. link_exists) THEN
            CALL h5pcreate_f(H5P_GROUP_CREATE_F, gcpl_id, ierr)
            IF (ierr /= 0) CALL errr(__FILE__, __LINE__)

            IF (track) THEN
                ! Set order tracking and indexing, useful for frequency information
                CALL h5pset_link_creation_order_f(gcpl_id, &
                    IOR(H5P_CRT_ORDER_TRACKED_F, H5P_CRT_ORDER_INDEXED_F), &
                    ierr)
                IF (ierr /= 0) CALL errr(__FILE__, __LINE__)
            END IF

            CALL h5gcreate_f(parent_id, name, group_id, ierr, gcpl_id=gcpl_id)
            IF (ierr /= 0) THEN
                WRITE(*, *) "Group creation failed: ", name
                CALL errr(__FILE__, __LINE__)
            END IF

            CALL h5pclose_f(gcpl_id, ierr)
            IF (ierr < 0) CALL errr(__FILE__, __LINE__)
        ELSE
            CALL h5gopen_f(parent_id, name, group_id, ierr)
            IF (ierr /= 0) CALL errr(__FILE__, __LINE__)
        END IF
    END SUBROUTINE hdf5common_group_open


    SUBROUTINE hdf5common_group_close(group_id)
        ! Subroutine arguments
        INTEGER(HID_T), INTENT(IN) :: group_id

        ! Local variables
        INTEGER(kind=int32) :: ierr

        IF (.NOT. ioProc) RETURN

        CALL h5gclose_f(group_id, ierr)
        IF (ierr /= 0) CALL errr(__FILE__, __LINE__)
    END SUBROUTINE hdf5common_group_close


    ! This subroutine create or open a dataset with given dimensions.
    ! Intented mainly for writing out data.
    !
    ! If you want to read in a dataset that should exist before, it is better
    ! to use the analogue "read" routine, this assure that any dataset
    ! is not created.
    SUBROUTINE hdf5common_dataset_create(name, shape, dtype, parent_id, &
            dset_id, filespace, maxdims, chunksize)
        ! Subroutine arguments
        CHARACTER(LEN=*), INTENT(IN) :: name
        INTEGER(HSIZE_T), INTENT(IN) :: shape(:)
        INTEGER(HID_T), INTENT(IN) :: dtype
        INTEGER(HID_T), INTENT(IN) :: parent_id
        INTEGER(HID_T), INTENT(OUT) :: dset_id
        INTEGER(HID_T), INTENT(OUT) :: filespace
        INTEGER(HSIZE_T), INTENT(IN), OPTIONAL :: maxdims(:)
        INTEGER(HSIZE_T), INTENT(IN), OPTIONAL :: chunksize(:)

        ! Local variables
        INTEGER(kind=int32) :: ierr
        LOGICAL :: link_exists

        INTEGER :: i, rank
        INTEGER(HSIZE_T) :: dims(maxrank)
        INTEGER(HSIZE_T) :: maxdims2(maxrank)

        INTEGER(HID_T) :: filetype, dcpl
        LOGICAL :: is_equal

        IF (.NOT. ioProc) THEN
            dset_id = -1
            filespace = -1
            RETURN
        END IF

        ! Sanity check
        IF (PRESENT(maxdims)) THEN
            IF (SIZE(maxdims) /= SIZE(shape)) THEN
                CALL errr(__FILE__, __LINE__)
            END IF
        END IF
        IF (PRESENT(chunksize)) THEN
            IF (SIZE(chunksize) /= SIZE(shape)) THEN
                CALL errr(__FILE__, __LINE__)
            END IF
        END IF

        ! If field does not exist, create it
        CALL h5lexists_f(parent_id, name, link_exists, ierr)
        IF (ierr < 0) CALL errr(__FILE__, __LINE__)

        IF (.NOT. link_exists) THEN
            rank = SIZE(shape)
            IF (rank > maxrank) CALL errr(__FILE__, __LINE__)
            DO i = 1, rank
                dims(i) = shape(i)
            END DO

            ! Create dataspace in file
            CALL h5screate_simple_f(rank, dims, filespace, ierr, &
                maxdims=maxdims)
            IF (ierr < 0) CALL errr(__FILE__, __LINE__)

            ! Dataset creation property list for optional chunking
            CALL h5pcreate_f(H5P_DATASET_CREATE_F, dcpl, ierr)
            IF (PRESENT(chunksize)) THEN
                CALL h5pset_chunk_f(dcpl, rank, chunksize, ierr)
                IF (ierr < 0) CALL errr(__FILE__, __LINE__)
            END IF

            ! Create dataset in file with default properties
            CALL h5dcreate_f(parent_id, name, dtype, filespace, &
                dset_id, ierr, dcpl_id=dcpl)
            IF (ierr < 0) CALL errr(__FILE__, __LINE__)

            CALL h5pclose_f(dcpl, ierr)
            IF (ierr < 0) CALL errr(__FILE__, __LINE__)
        ELSE
            CALL h5dopen_f(parent_id, name, dset_id, ierr)
            IF (ierr < 0) CALL errr(__FILE__, __LINE__)

            ! get dataspace from opened dataset
            CALL h5dget_space_f(dset_id, filespace, ierr)
            IF (ierr < 0) CALL errr(__FILE__, __LINE__)

            ! check dataspace dimension
            CALL h5sget_simple_extent_dims_f(filespace, dims, maxdims2, ierr)
            ! ierror = rank when successful, -1 when failure
            IF (ierr < 0) CALL errr(__FILE__, __LINE__)

            ! Check that dataspace is same dims. as requested
            IF (ierr > maxrank) THEN
                WRITE(*, *) "Dataset rank in file: ", ierr
                WRITE(*, *) "Maximum rank: ", maxrank
                CALL errr(__FILE__, __LINE__)
            END IF
            IF (ierr /= SIZE(shape)) THEN
                WRITE(*, *) "Error opening dataset: ", name
                WRITE(*, *) "Griven shape parameter: ", shape
                WRITE(*, *) "Dataset rank in file: ", ierr
                CALL errr(__FILE__, __LINE__)
            END IF
            DO i = 1, SIZE(shape)
                IF (dims(i) /= shape(i)) THEN
                    WRITE(*, *) "Error opening dataset: ", name
                    WRITE(*, *) "Griven shape parameter: ", shape
                    WRITE(*, *) "Dataset shape in file: ", dims(1:ierr)
                    CALL errr(__FILE__, __LINE__)
                END IF
            END DO

            ! Get datatype from dataset and check if the type of the dataset
            ! is the same as requested
            CALL h5dget_type_f(dset_id, filetype, ierr)
            IF (ierr < 0) CALL errr(__FILE__, __LINE__)

            CALL h5tequal_f(filetype, dtype, is_equal, ierr)
            IF (ierr < 0) CALL errr(__FILE__, __LINE__)
            IF (.NOT. is_equal) THEN
                WRITE(*, *) "Fatal error: datatype in opened dataset is not"
                WRITE(*, *) "equal to requested datatype."
                WRITE(*, *) "Dataset name: ", name
                CALL errr(__FILE__, __LINE__)
            END IF

            CALL h5tclose_f(filetype, ierr)
            IF (ierr < 0) CALL errr(__FILE__, __LINE__)
        END IF
    END SUBROUTINE hdf5common_dataset_create


    ! This subroutine open a dataset and return dimensions.
    ! Intented mainly for reading in out data. If the dataset
    ! is not present it will fail.
    !
    ! The "shape" parameter must have the same size as the expected shape of
    ! the HDF5 data array
    SUBROUTINE hdf5common_dataset_open(name, shape, parent_id, dset_id, &
            filespace)
        ! Subroutine arguments
        CHARACTER(LEN=*), INTENT(IN) :: name
        INTEGER(HSIZE_T), INTENT(INOUT) :: shape(:)
        INTEGER(HID_T), INTENT(IN) :: parent_id
        INTEGER(HID_T), INTENT(OUT) :: dset_id
        INTEGER(HID_T), INTENT(OUT) :: filespace

        ! Local variables
        INTEGER(kind=int32) :: ierr
        LOGICAL :: link_exists

        INTEGER :: i
        INTEGER(HSIZE_T) :: dims(maxrank)
        INTEGER(HSIZE_T) :: maxdims(maxrank)

        INTEGER(SIZE_T), PARAMETER :: charlen = 64
        INTEGER(SIZE_T) :: glen, flen
        CHARACTER(len=charlen) :: filename
        CHARACTER(len=charlen) :: groupname

        IF (.NOT. ioProc) THEN
            dset_id = -1
            filespace = -1
            shape = -1
            RETURN
        END IF

        ! If field does not exist, fail.
        CALL h5lexists_f(parent_id, name, link_exists, ierr)
        IF (ierr /= 0) CALL errr(__FILE__, __LINE__)

        IF (.NOT. link_exists) THEN
            CALL h5iget_name_f(parent_id, groupname, charlen, glen, ierr)
            IF (ierr < 0) CALL errr(__FILE__, __LINE__)

            CALL h5fget_name_f(parent_id, filename, flen, ierr)
            IF (ierr < 0) CALL errr(__FILE__, __LINE__)

            WRITE(*, *) "Fatal error in opening dataset: ", name
            WRITE(*, *) "  group: ", groupname(1:glen)
            WRITE(*, *) "   file: ", filename(1:flen)
            WRITE(*, *) "Dataset does not exist!"
            CALL errr(__FILE__, __LINE__)
        END IF

        CALL h5dopen_f(parent_id, name, dset_id, ierr)
        IF (ierr < 0) CALL errr(__FILE__, __LINE__)

        ! get dataspace from opened dataset
        CALL h5dget_space_f(dset_id, filespace, ierr)
        IF (ierr < 0) CALL errr(__FILE__, __LINE__)

        ! check dataspace dimension
        CALL h5sget_simple_extent_dims_f(filespace, dims, maxdims, ierr)
        ! ierror = rank when successful, -1 when failure
        IF (ierr < 0) CALL errr(__FILE__, __LINE__)

        ! Check that dataspace is same dims. as requested
        IF (ierr > maxrank) THEN
            WRITE(*, *) "Dataset rank in file: ", ierr
            WRITE(*, *) "Maximum rank: ", maxrank
            CALL errr(__FILE__, __LINE__)
        END IF
        IF (ierr /= SIZE(shape)) THEN
            WRITE(*, *) "Error opening dataset: ", name
            WRITE(*, *) "Size of 'shape' parameter: ", SIZE(shape)
            WRITE(*, *) "Dataset rank in file: ", ierr
            CALL errr(__FILE__, __LINE__)
        END IF
        DO i = 1, SIZE(shape)
            shape(i) = dims(i)
        END DO
    END SUBROUTINE hdf5common_dataset_open


    ! By default this will broadcast the result to all processes  and all
    ! processes miust call this routine. This is how it is usually used
    ! from for instacne probes and fieldio2. However, there are also
    ! occasions where this is used from only IO processes and then
    ! the optional bcast parameter can be set to false to prevent a deadlock.
    SUBROUTINE hdf5common_dataset_exists(name, parent_id, link_exists, shape, &
            bcast)
        ! Subroutine arguments
        CHARACTER(LEN=*), INTENT(IN) :: name
        INTEGER(HID_T), INTENT(IN) :: parent_id
        LOGICAL, INTENT(OUT) :: link_exists
        INTEGER(HSIZE_T), INTENT(OUT), OPTIONAL :: shape(:)
        LOGICAL, INTENT(IN), OPTIONAL :: bcast

        ! Local variables
        INTEGER(HID_T) :: dset_id, filespace
        INTEGER(HSIZE_T) :: dims(maxrank)
        INTEGER(HSIZE_T) :: maxdims(maxrank)
        INTEGER(intk) :: i, ndims
        INTEGER(int32) :: ierr
        LOGICAL :: bcast2

        IF (ioproc) THEN
            CALL h5lexists_f(parent_id, name, link_exists, ierr)
            IF (ierr < 0) CALL errr(__FILE__, __LINE__)
        END IF

        bcast2 = .TRUE.
        IF (PRESENT(bcast)) THEN
            bcast2 = bcast
        END IF
        IF (bcast2) THEN
            CALL MPI_Bcast(link_exists, 1, MPI_LOGICAL, 0, iogrcomm)
        END IF

        IF (link_exists .AND. PRESENT(shape)) THEN
            IF (ioproc) THEN
                CALL h5dopen_f(parent_id, name, dset_id, ierr)
                IF (ierr < 0) CALL errr(__FILE__, __LINE__)

                ! get dataspace from opened dataset
                CALL h5dget_space_f(dset_id, filespace, ierr)
                IF (ierr < 0) CALL errr(__FILE__, __LINE__)

                ! check dataspace dimension
                CALL h5sget_simple_extent_dims_f(filespace, dims, maxdims, ierr)
                ! ierror = rank when successful, -1 when failure
                IF (ierr < 0) CALL errr(__FILE__, __LINE__)

                ! Check that dataspace is same dims. as requested
                IF (ierr > maxrank) THEN
                    WRITE(*, *) "Dataset rank in file: ", ierr
                    WRITE(*, *) "Maximum rank: ", maxrank
                    CALL errr(__FILE__, __LINE__)
                END IF
                IF (ierr /= SIZE(shape)) THEN
                    WRITE(*, *) "Error opening dataset: ", name
                    WRITE(*, *) "Size of 'shape' parameter: ", SIZE(shape)
                    WRITE(*, *) "Dataset rank in file: ", ierr
                    CALL errr(__FILE__, __LINE__)
                END IF
                ! Copy to output array
                ndims = SIZE(shape)
                DO i = 1, ndims
                    shape(i) = dims(i)
                END DO

                ! close dataspace and dataset
                CALL h5sclose_f(filespace, ierr)
                IF (ierr < 0) CALL errr(__FILE__, __LINE__)
                CALL h5dclose_f(dset_id, ierr)
                IF (ierr < 0) CALL errr(__FILE__, __LINE__)
            END IF

            ! IO processes now has the shape of the dataset, broadcast this
            IF (bcast2) THEN
                CALL MPI_Bcast(shape, INT(ndims, int32), mglet_mpi_hsize_t, &
                    0, iogrcomm)
            END IF
        END IF
    END SUBROUTINE hdf5common_dataset_exists


    ! Takes an already open dset_id and filepace, extends it to the desired
    ! shape
    SUBROUTINE hdf5common_dataset_extend(dset_id, filespace, shape, &
            offset, count)

        ! Subroutine arguments
        INTEGER(HID_T), INTENT(IN) :: dset_id
        INTEGER(HID_T), INTENT(INOUT) :: filespace
        INTEGER(HSIZE_T), INTENT(IN) :: shape(:)
        INTEGER(HSIZE_T), INTENT(IN), OPTIONAL :: offset(:)
        INTEGER(HSIZE_T), INTENT(IN), OPTIONAL :: count(:)

        ! Local variables
        INTEGER(HSIZE_T) :: dims(maxrank)
        INTEGER(HSIZE_T) :: maxdims(maxrank)
        INTEGER(intk) :: i, ndims
        INTEGER(int32) :: ierr

        IF (.NOT. ioproc) RETURN

        ! check dataspace dimension
        CALL h5sget_simple_extent_dims_f(filespace, dims, maxdims, ierr)
        ! ierror = rank when successful, -1 when failure
        IF (ierr < 0) CALL errr(__FILE__, __LINE__)
        ndims = ierr

        IF (ndims /= SIZE(shape)) THEN
            CALL errr(__FILE__, __LINE__)
        END IF

        DO i = 1, ndims
            IF (shape(i) > maxdims(i) .AND. maxdims(i) /= H5S_UNLIMITED_F) THEN
                WRITE(*, *) "i:          ", i
                WRITE(*, *) "shape(i):   ", shape(i)
                WRITE(*, *) "maxdims(i): ", maxdims(i)
                CALL errr(__FILE__, __LINE__)
            END IF
        END DO

        CALL h5sclose_f(filespace, ierr)
        IF (ierr /= 0) CALL errr(__FILE__, __LINE__)

        CALL h5dset_extent_f(dset_id, shape, ierr)
        IF (ierr /= 0) CALL errr(__FILE__, __LINE__)

        CALL h5dget_space_f(dset_id, filespace, ierr)
        IF (ierr /= 0) CALL errr(__FILE__, __LINE__)

        IF (PRESENT(offset) .AND. PRESENT(count)) THEN
            CALL h5sselect_hyperslab_f(filespace, H5S_SELECT_SET_F, offset, &
                count, ierr)
            IF (ierr /= 0) CALL errr(__FILE__, __LINE__)
        END IF
    END SUBROUTINE hdf5common_dataset_extend


    SUBROUTINE hdf5common_dataset_close(dset_id, filespace)
        ! Subroutine arguments
        INTEGER(HID_T), INTENT(IN) :: dset_id
        INTEGER(HID_T), INTENT(IN) :: filespace

        ! Local variables
        INTEGER(kind=int32) :: ierr

        IF (.NOT. ioProc) RETURN

        ! close dataspace
        CALL h5sclose_f(filespace, ierr)
        IF (ierr < 0) CALL errr(__FILE__, __LINE__)

        CALL h5dclose_f(dset_id, ierr)
        IF (ierr < 0) CALL errr(__FILE__, __LINE__)
    END SUBROUTINE hdf5common_dataset_close


    ! This routine write a scalar attribute to an HDF5 object
    ! The contents of the "value" must be identical across all IO processes,
    ! and is completely ignored at non-io processes
    SUBROUTINE hdf5common_attr_write(name, value, parent_id, parent_name)
        ! Subroutine arguments
        CHARACTER(LEN=*), INTENT(IN) :: name
        CLASS(*), INTENT(in), TARGET :: value
        INTEGER(HID_T), INTENT(IN) :: parent_id
        CHARACTER(LEN=*), INTENT(IN), OPTIONAL :: parent_name

        ! Local variables
        INTEGER(kind=int32) :: ierr
        LOGICAL :: attr_exists, is_equal
        INTEGER(HID_T) :: hdf5_dtype, type_id
        INTEGER(HID_T) :: aspace_id, attr_id, obj_id
        INTEGER(HSIZE_T) :: dims(maxrank)
        INTEGER(HSIZE_T) :: maxdims(maxrank)
        TYPE(C_PTR) :: buf

        IF (.NOT. ioProc) RETURN

        SELECT TYPE (value)
            TYPE IS (REAL(realk))
                hdf5_dtype = mglet_hdf5_real
                buf = C_LOC(value)
            TYPE IS (INTEGER(intk))
                hdf5_dtype = mglet_hdf5_int
                buf = C_LOC(value)
            TYPE IS (CHARACTER(len=*))
                CALL h5tcopy_f(H5T_NATIVE_CHARACTER, hdf5_dtype, ierr)
                IF (ierr < 0) CALL errr(__FILE__, __LINE__)
                CALL h5tset_size_f(hdf5_dtype, INT(LEN(value), size_t), ierr)
                IF (ierr < 0) CALL errr(__FILE__, __LINE__)
                buf = C_LOC(value)
            CLASS DEFAULT
                CALL errr(__FILE__, __LINE__)
        END SELECT

        ! Check whether attribute exists
        CALL hdf5common_attr_exists(name, attr_exists, parent_id, parent_name)

        ! If not exist - create
        IF (.NOT. attr_exists) THEN
            CALL h5screate_f(H5S_SCALAR_F, aspace_id, ierr)
            IF (ierr < 0) CALL errr(__FILE__, __LINE__)

            IF (PRESENT(parent_name)) THEN
                CALL h5acreate_by_name_f(parent_id, parent_name, name, &
                    hdf5_dtype, aspace_id, attr_id, ierr)
                IF (ierr < 0) CALL errr(__FILE__, __LINE__)
            ELSE
                CALL h5acreate_f(parent_id, name, hdf5_dtype, aspace_id, &
                    attr_id, ierr)
                IF (ierr < 0) CALL errr(__FILE__, __LINE__)
            END IF
        ELSE
            IF (PRESENT(parent_name)) THEN
                ! This should in theory do the same thing, but does not work
                ! HDF5 library bug?
                ! CALL h5aopen_by_name_f(parent_id, parent_name, name, &
                !    attr_id, ierr)
                CALL h5oopen_f(parent_id, parent_name, obj_id, ierr)
                IF (ierr < 0) CALL errr(__FILE__, __LINE__)

                CALL h5aopen_f(obj_id, name, attr_id, ierr)
                IF (ierr < 0) CALL errr(__FILE__, __LINE__)

                CALL h5oclose_f(obj_id, ierr)
                IF (ierr < 0) CALL errr(__FILE__, __LINE__)
            ELSE
                CALL h5aopen_f(parent_id, name, attr_id, ierr)
                IF (ierr < 0) CALL errr(__FILE__, __LINE__)
            END IF

            CALL h5aget_space_f(attr_id, aspace_id, ierr)
            IF (ierr < 0) CALL errr(__FILE__, __LINE__)

            CALL h5sget_simple_extent_dims_f(aspace_id, dims, maxdims, ierr)
            ! ierror = rank when successful, -1 when failure
            IF (ierr < 0) CALL errr(__FILE__, __LINE__)

            ! Check that dataspace is scalar
            IF (ierr /= 0) THEN
                WRITE(*, *) "Error opening attribute: ", name
                WRITE(*, *) "Attribute rank in file: ", ierr
                CALL errr(__FILE__, __LINE__)
            END IF

            ! Get type and check that type is the same
            CALL h5aget_type_f(attr_id, type_id, ierr)
            CALL h5tequal_f(hdf5_dtype, type_id, is_equal, ierr)
            IF (ierr < 0) CALL errr(__FILE__, __LINE__)
            IF (.NOT. is_equal) THEN
                WRITE(*, *) "Fatal error: datatype in opened attribute is not"
                WRITE(*, *) "equal to requested datatype."
                WRITE(*, *) "Attribute name: ", name
                CALL errr(__FILE__, __LINE__)
            END IF

            CALL h5tclose_f(type_id, ierr)
            IF (ierr < 0) CALL errr(__FILE__, __LINE__)
        END IF

        CALL h5sclose_f(aspace_id, ierr)
        IF (ierr < 0) CALL errr(__FILE__, __LINE__)

        CALL h5awrite_f(attr_id, hdf5_dtype, buf, ierr)
        IF (ierr < 0) CALL errr(__FILE__, __LINE__)

        CALL h5aclose_f(attr_id, ierr)
        IF (ierr < 0) CALL errr(__FILE__, __LINE__)

        SELECT TYPE (value)
            TYPE IS (CHARACTER(len=*))
                CALL h5tclose_f(hdf5_dtype, ierr)
                IF (ierr < 0) CALL errr(__FILE__, __LINE__)
        END SELECT
    END SUBROUTINE hdf5common_attr_write


    ! This routine write a scalar attribute to an HDF5 object
    ! The contents of the "value" must be identical across all IO processes,
    ! and is completely ignored at non-io processes
    SUBROUTINE hdf5common_attr_write_arr(name, value, parent_id, parent_name)
        ! Subroutine arguments
        CHARACTER(LEN=*), INTENT(IN) :: name
        CLASS(*), INTENT(in), TARGET :: value(:)
        INTEGER(HID_T), INTENT(IN) :: parent_id
        CHARACTER(LEN=*), INTENT(IN), OPTIONAL :: parent_name

        ! Local variables
        INTEGER(kind=int32) :: ierr
        LOGICAL :: attr_exists, is_equal
        INTEGER(HID_T) :: hdf5_dtype, type_id
        INTEGER(HID_T) :: aspace_id, attr_id, obj_id
        INTEGER(HSIZE_T) :: dims(maxrank)
        INTEGER(HSIZE_T) :: maxdims(maxrank)
        TYPE(C_PTR) :: buf
        INTEGER(intk) :: nelems

        IF (.NOT. ioProc) RETURN

        SELECT TYPE (value)
            TYPE IS (REAL(realk))
                hdf5_dtype = mglet_hdf5_real
                buf = C_LOC(value)
                nelems = SIZE(value)
            TYPE IS (INTEGER(intk))
                hdf5_dtype = mglet_hdf5_int
                buf = C_LOC(value)
                nelems = SIZE(value)
            TYPE IS (CHARACTER(*))
                CALL h5tcopy_f(H5T_NATIVE_CHARACTER, hdf5_dtype, ierr)
                IF (ierr < 0) CALL errr(__FILE__, __LINE__)
                CALL h5tset_size_f(hdf5_dtype, LEN(value, size_t), ierr)
                IF (ierr < 0) CALL errr(__FILE__, __LINE__)
                buf = C_LOC(value)
                nelems = SIZE(value)
            CLASS DEFAULT
                nelems = 0  ! Too ensure compiler see nelems as uninitialized
                CALL errr(__FILE__, __LINE__)
        END SELECT

        ! Check whether attribute exists
        CALL hdf5common_attr_exists(name, attr_exists, parent_id, parent_name)

        ! If not exist - create
        IF (.NOT. attr_exists) THEN
            dims(1) = nelems
            CALL h5screate_simple_f(1, dims, aspace_id, ierr)
            IF (ierr < 0) CALL errr(__FILE__, __LINE__)

            IF (PRESENT(parent_name)) THEN
                CALL h5acreate_by_name_f(parent_id, parent_name, name, &
                    hdf5_dtype, aspace_id, attr_id, ierr)
                IF (ierr < 0) CALL errr(__FILE__, __LINE__)
            ELSE
                CALL h5acreate_f(parent_id, name, hdf5_dtype, aspace_id, &
                    attr_id, ierr)
                IF (ierr < 0) CALL errr(__FILE__, __LINE__)
            END IF
        ELSE
            IF (PRESENT(parent_name)) THEN
                ! This should in theory do the same thing, but does not work
                ! HDF5 library bug?
                ! CALL h5aopen_by_name_f(parent_id, parent_name, name, &
                !    attr_id, ierr)
                CALL h5oopen_f(parent_id, parent_name, obj_id, ierr)
                IF (ierr < 0) CALL errr(__FILE__, __LINE__)

                CALL h5aopen_f(obj_id, name, attr_id, ierr)
                IF (ierr < 0) CALL errr(__FILE__, __LINE__)

                CALL h5oclose_f(obj_id, ierr)
                IF (ierr < 0) CALL errr(__FILE__, __LINE__)
            ELSE
                CALL h5aopen_f(parent_id, name, attr_id, ierr)
                IF (ierr < 0) CALL errr(__FILE__, __LINE__)
            END IF

            CALL h5aget_space_f(attr_id, aspace_id, ierr)
            IF (ierr < 0) CALL errr(__FILE__, __LINE__)

            CALL h5sget_simple_extent_dims_f(aspace_id, dims, maxdims, ierr)
            ! ierror = rank when successful, -1 when failure
            IF (ierr < 0) CALL errr(__FILE__, __LINE__)

            ! Check that dataspace is rank-1
            IF (ierr /= 1) THEN
                WRITE(*, *) "Error opening attribute: ", name
                WRITE(*, *) "Attribute rank in file: ", ierr
                CALL errr(__FILE__, __LINE__)
            END IF

            ! Get type and check that type is the same
            CALL h5aget_type_f(attr_id, type_id, ierr)
            CALL h5tequal_f(hdf5_dtype, type_id, is_equal, ierr)
            IF (ierr < 0) CALL errr(__FILE__, __LINE__)
            IF (.NOT. is_equal) THEN
                WRITE(*, *) "Fatal error: datatype in opened attribute is not"
                WRITE(*, *) "equal to requested datatype."
                WRITE(*, *) "Attribute name: ", name
                CALL errr(__FILE__, __LINE__)
            END IF

            CALL h5tclose_f(type_id, ierr)
            IF (ierr < 0) CALL errr(__FILE__, __LINE__)
        END IF

        CALL h5sclose_f(aspace_id, ierr)
        IF (ierr < 0) CALL errr(__FILE__, __LINE__)

        CALL h5awrite_f(attr_id, hdf5_dtype, buf, ierr)
        IF (ierr < 0) CALL errr(__FILE__, __LINE__)

        CALL h5aclose_f(attr_id, ierr)
        IF (ierr < 0) CALL errr(__FILE__, __LINE__)
    END SUBROUTINE hdf5common_attr_write_arr


    SUBROUTINE hdf5common_attr_exists(name, attr_exists, parent_id, &
            parent_name, shape)
        ! Subroutine arguments
        CHARACTER(LEN=*), INTENT(IN) :: name
        LOGICAL, INTENT(out) :: attr_exists
        INTEGER(HID_T), INTENT(IN) :: parent_id
        CHARACTER(LEN=*), INTENT(IN), OPTIONAL :: parent_name
        INTEGER(HSIZE_T), INTENT(OUT), OPTIONAL :: shape(:)

        ! Local variables
        INTEGER(kind=int32) :: ierr
        INTEGER(intk) :: i
        INTEGER(HID_T) :: obj_id, attr_id, aspace_id
        INTEGER(HSIZE_T) :: dims(maxrank), maxdims(maxrank)

        attr_exists = .FALSE.
        IF (PRESENT(shape)) THEN
            shape = 0
        END IF

        IF (.NOT. ioProc) RETURN

        ! Check whether attribute exists
        IF (PRESENT(parent_name)) THEN
            CALL h5aexists_by_name_f(parent_id, parent_name, name, &
                attr_exists, ierr)
            IF (ierr < 0) CALL errr(__FILE__, __LINE__)
        ELSE
            CALL h5aexists_f(parent_id, name, attr_exists, ierr)
            IF (ierr < 0) CALL errr(__FILE__, __LINE__)
        END IF

        IF (attr_exists .AND. PRESENT(shape)) THEN
            IF (PRESENT(parent_name)) THEN
                ! This should in theory do the same thing, but does not work
                ! HDF5 library bug?
                ! CALL h5aopen_by_name_f(parent_id, parent_name, name, &
                !    attr_id, ierr)
                CALL h5oopen_f(parent_id, parent_name, obj_id, ierr)
                IF (ierr < 0) CALL errr(__FILE__, __LINE__)

                CALL h5aopen_f(obj_id, name, attr_id, ierr)
                IF (ierr < 0) CALL errr(__FILE__, __LINE__)

                CALL h5oclose_f(obj_id, ierr)
                IF (ierr < 0) CALL errr(__FILE__, __LINE__)
            ELSE
                CALL h5aopen_f(parent_id, name, attr_id, ierr)
                IF (ierr < 0) CALL errr(__FILE__, __LINE__)
            END IF

            CALL h5aget_space_f(attr_id, aspace_id, ierr)
            IF (ierr < 0) CALL errr(__FILE__, __LINE__)

            CALL h5sget_simple_extent_dims_f(aspace_id, dims, maxdims, ierr)
            ! ierror = rank when successful, -1 when failure
            IF (ierr < 0) CALL errr(__FILE__, __LINE__)

            IF (ierr > maxrank) THEN
                WRITE(*, *) "Dataset rank in file: ", ierr
                WRITE(*, *) "Maximum rank: ", maxrank
                CALL errr(__FILE__, __LINE__)
            END IF

            IF (ierr /= SIZE(shape)) THEN
                WRITE(*, *) "Attribute rank in file: ", ierr
                WRITE(*, *) "Size of shape array: ",  SIZE(shape)
                CALL errr(__FILE__, __LINE__)
            END IF

            DO i = 1, SIZE(shape)
                shape(i) = dims(i)
            END DO

            CALL h5aclose_f(attr_id, ierr)
        END IF
    END SUBROUTINE hdf5common_attr_exists


    ! This routine read in an attribute from a HDF5 object.
    ! The value is only set at IO processes, and no changes to the
    ! "value" happens at non-IO processes.
    ! If the value is to be distrinbuted to all MPI ranks, see routine
    ! hdf5_attr_read_collective
    SUBROUTINE hdf5common_attr_read(name, value, parent_id, parent_name)
        ! Subroutine arguments
        CHARACTER(LEN=*), INTENT(IN) :: name
        CLASS(*), INTENT(INOUT), TARGET :: value
        INTEGER(HID_T), INTENT(IN) :: parent_id
        CHARACTER(LEN=*), INTENT(IN), OPTIONAL :: parent_name

        ! Local variables
        INTEGER(kind=int32) :: ierr
        LOGICAL :: attr_exists
        INTEGER(HID_T) :: hdf5_dtype
        INTEGER(HID_T) :: aspace_id, attr_id
        INTEGER(HSIZE_T) :: dims(maxrank)
        INTEGER(HSIZE_T) :: maxdims(maxrank)
        TYPE(C_PTR) :: buf

        IF (.NOT. ioProc) RETURN

        SELECT TYPE (value)
            TYPE IS (REAL(realk))
                hdf5_dtype = mglet_hdf5_real
                buf = C_LOC(value)
            TYPE IS (INTEGER(intk))
                hdf5_dtype = mglet_hdf5_int
                buf = C_LOC(value)
            TYPE IS (CHARACTER(len=*))
                CALL h5tcopy_f(H5T_NATIVE_CHARACTER, hdf5_dtype, ierr)
                IF (ierr < 0) CALL errr(__FILE__, __LINE__)
                CALL h5tset_size_f(hdf5_dtype, INT(LEN(value), size_t), ierr)
                IF (ierr < 0) CALL errr(__FILE__, __LINE__)
                buf = C_LOC(value)
            CLASS DEFAULT
                CALL errr(__FILE__, __LINE__)
        END SELECT

        ! Check whether attribute exists
        CALL hdf5common_attr_exists(name, attr_exists, parent_id, parent_name)

        IF (.NOT. attr_exists) THEN
            CALL errr(__FILE__, __LINE__)
        END IF

        IF (PRESENT(parent_name)) THEN
            CALL h5aopen_by_name_f(parent_id, parent_name, name, attr_id, ierr)
            IF (ierr < 0) CALL errr(__FILE__, __LINE__)
        ELSE
            CALL h5aopen_f(parent_id, name, attr_id, ierr)
            IF (ierr < 0) CALL errr(__FILE__, __LINE__)
        END IF

        CALL h5aget_space_f(attr_id, aspace_id, ierr)
        IF (ierr < 0) CALL errr(__FILE__, __LINE__)

        CALL h5sget_simple_extent_dims_f(aspace_id, dims, maxdims, ierr)
        ! ierror = rank when successful, -1 when failure
        IF (ierr < 0) CALL errr(__FILE__, __LINE__)

        ! Check that dataspace is scalar or 1-D
        IF (ierr > 1) THEN
            WRITE(*, *) "Error opening attribute: ", name
            WRITE(*, *) "Attribute rank in file: ", ierr
            CALL errr(__FILE__, __LINE__)
        ! In case 1-D the dims must be 1
        ELSE IF (ierr == 1) THEN
            IF (dims(1) /= 1) THEN
                WRITE(*, *) "Error opening attribute: ", name
                WRITE(*, *) "Attribute rank in file: ", ierr
                WRITE(*, *) "Attribute shape: ", dims(1:ierr)
                CALL errr(__FILE__, __LINE__)
            END IF
        END IF

        CALL h5sclose_f(aspace_id, ierr)
        IF (ierr < 0) CALL errr(__FILE__, __LINE__)

        CALL h5aread_f(attr_id, hdf5_dtype, buf, ierr)
        IF (ierr < 0) CALL errr(__FILE__, __LINE__)

        CALL h5aclose_f(attr_id, ierr)
        IF (ierr < 0) CALL errr(__FILE__, __LINE__)
    END SUBROUTINE hdf5common_attr_read


    ! This routine read in an attribute from a HDF5 object.
    ! The value is only set at IO processes, and no changes to the
    ! "value" happens at non-IO processes.
    ! If the value is to be distrinbuted to all MPI ranks, see routine
    ! hdf5_attr_read_arr_collective
    SUBROUTINE hdf5common_attr_read_arr(name, value, parent_id, parent_name)
        ! Subroutine arguments
        CHARACTER(LEN=*), INTENT(IN) :: name
        CLASS(*), INTENT(INOUT), TARGET :: value(:)
        INTEGER(HID_T), INTENT(IN) :: parent_id
        CHARACTER(LEN=*), INTENT(IN), OPTIONAL :: parent_name

        ! Local variables
        INTEGER(kind=int32) :: ierr
        LOGICAL :: attr_exists
        INTEGER(HID_T) :: hdf5_dtype, type_id
        INTEGER(HID_T) :: aspace_id, attr_id
        INTEGER(HSIZE_T) :: dims(maxrank)
        INTEGER(HSIZE_T) :: maxdims(maxrank)
        TYPE(C_PTR) :: buf
        INTEGER(intk) :: nelems

        ! NB: Default integer, not an error!
        ! Ref. h5tget_class_f API documentation
        INTEGER :: class_file, class_value

        IF (.NOT. ioProc) RETURN

        SELECT TYPE (value)
            TYPE IS (REAL(realk))
                hdf5_dtype = mglet_hdf5_real
                buf = C_LOC(value)
                nelems = SIZE(value)
            TYPE IS (INTEGER(intk))
                hdf5_dtype = mglet_hdf5_int
                buf = C_LOC(value)
                nelems = SIZE(value)
            CLASS DEFAULT
                CALL errr(__FILE__, __LINE__)
        END SELECT

        ! Check whether attribute exists
        CALL hdf5common_attr_exists(name, attr_exists, parent_id, parent_name)

        IF (.NOT. attr_exists) THEN
            CALL errr(__FILE__, __LINE__)
        END IF

        IF (PRESENT(parent_name)) THEN
            CALL h5aopen_by_name_f(parent_id, parent_name, name, attr_id, ierr)
            IF (ierr < 0) CALL errr(__FILE__, __LINE__)
        ELSE
            CALL h5aopen_f(parent_id, name, attr_id, ierr)
            IF (ierr < 0) CALL errr(__FILE__, __LINE__)
        END IF

        CALL h5aget_space_f(attr_id, aspace_id, ierr)
        IF (ierr < 0) CALL errr(__FILE__, __LINE__)

        CALL h5sget_simple_extent_dims_f(aspace_id, dims, maxdims, ierr)
        ! ierror = rank when successful, -1 when failure
        IF (ierr < 0) CALL errr(__FILE__, __LINE__)

        ! Check that dataspace is rank-1
        IF (ierr /= 1) THEN
            WRITE(*, *) "Error opening attribute: ", name
            WRITE(*, *) "Attribute rank in file: ", ierr
            CALL errr(__FILE__, __LINE__)
        END IF

        ! Check that dataset to be read is same as memory buffer size
        IF (dims(1) /= nelems) THEN
            WRITE(*, *) "Error opening attribute: ", name
            WRITE(*, *) "Attribute size in file: ", dims(1)
            WRITE(*, *) "Memory buffer size: ", nelems
            CALL errr(__FILE__, __LINE__)
        END IF

        CALL h5sclose_f(aspace_id, ierr)
        IF (ierr < 0) CALL errr(__FILE__, __LINE__)

        ! Get type and check that type is the same
        CALL h5aget_type_f(attr_id, type_id, ierr)
        IF (ierr < 0) CALL errr(__FILE__, __LINE__)

        CALL h5tget_class_f(type_id, class_file, ierr)
        IF (ierr < 0) CALL errr(__FILE__, __LINE__)

        CALL h5tclose_f(type_id, ierr)
        IF (ierr < 0) CALL errr(__FILE__, __LINE__)

        CALL h5tget_class_f(hdf5_dtype, class_value, ierr)
        IF (ierr < 0) CALL errr(__FILE__, __LINE__)
        IF (class_file /= class_value) THEN
            WRITE(*, *) "Fatal error: datatype in opened attribute is not"
            WRITE(*, *) "equal to requested datatype."
            WRITE(*, *) "Attribute name: ", name
            CALL errr(__FILE__, __LINE__)
        END IF

        CALL h5aread_f(attr_id, hdf5_dtype, buf, ierr)
        IF (ierr < 0) CALL errr(__FILE__, __LINE__)

        CALL h5aclose_f(attr_id, ierr)
        IF (ierr < 0) CALL errr(__FILE__, __LINE__)
    END SUBROUTINE hdf5common_attr_read_arr


    ! This routine read in an attribute from a HDF5 object.
    ! The hdf5common_attr_read is used to read the value at the IO processes,
    ! and then MPI is used to distribute the values to all processes
    SUBROUTINE hdf5common_attr_read_collective(name, value, parent_id, &
            parent_name)

        ! Subroutine arguments
        CHARACTER(LEN=*), INTENT(IN) :: name
        CLASS(*), INTENT(out) :: value
        INTEGER(HID_T), INTENT(IN) :: parent_id
        CHARACTER(LEN=*), INTENT(IN), OPTIONAL :: parent_name

        ! Local variables
        TYPE(MPI_Datatype) :: mpi_dtype

        CALL hdf5common_attr_read(name, value, parent_id, parent_name)

        SELECT TYPE (value)
            TYPE IS (REAL(realk))
                mpi_dtype = mglet_mpi_real
                CALL MPI_Bcast(value, 1, mpi_dtype, 0, iogrcomm)
            TYPE IS (INTEGER(intk))
                mpi_dtype = mglet_mpi_int
                CALL MPI_Bcast(value, 1, mpi_dtype, 0, iogrcomm)
            CLASS DEFAULT
                CALL errr(__FILE__, __LINE__)
        END SELECT
    END SUBROUTINE hdf5common_attr_read_collective


    ! This routine read in an attribute from a HDF5 object.
    ! The hdf5common_attr_read is used to read the value at the IO processes,
    ! and then MPI is used to distribute the values to all processes
    SUBROUTINE hdf5common_attr_read_arr_collective(name, value, parent_id, &
            parent_name)

        ! Subroutine arguments
        CHARACTER(LEN=*), INTENT(IN) :: name
        CLASS(*), INTENT(out) :: value(:)
        INTEGER(HID_T), INTENT(IN) :: parent_id
        CHARACTER(LEN=*), INTENT(IN), OPTIONAL :: parent_name

        ! Local variables
        TYPE(MPI_Datatype) :: mpi_dtype

        CALL hdf5common_attr_read_arr(name, value, parent_id, parent_name)

        SELECT TYPE (value)
            TYPE IS (REAL(realk))
                mpi_dtype = mglet_mpi_real
                CALL MPI_Bcast(value, SIZE(value), mpi_dtype, 0, iogrcomm)
            TYPE IS (INTEGER(intk))
                mpi_dtype = mglet_mpi_int
                CALL MPI_Bcast(value, SIZE(value), mpi_dtype, 0, iogrcomm)
            CLASS DEFAULT
                CALL errr(__FILE__, __LINE__)
        END SELECT
    END SUBROUTINE hdf5common_attr_read_arr_collective

END MODULE hdf5common_mod
