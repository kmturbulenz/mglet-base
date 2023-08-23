MODULE runinfo_mod
    USE, INTRINSIC :: ISO_C_BINDING, ONLY: c_char, c_double, c_ptr, c_long_long
    USE HDF5
    USE MPI_f08

    USE core_mod

    IMPLICIT NONE(type, external)
    PRIVATE

    INTEGER(intk), PARAMETER :: nchar_programname = 64

    ! derive data type to be used in runinfo table
    TYPE, BIND(C) :: runinfo_t
        ! basic information
        CHARACTER(kind=C_CHAR) :: programname(nchar_programname) = " "
        INTEGER(c_intk) :: starttime(8) = 0
        REAL(c_double) :: walltime = 0.0
        INTEGER(c_intk) :: numprocs = 0

        ! Basics about grid - useful for performance numbers
        INTEGER(c_intk) :: ngrids = 0
        INTEGER(c_long_long) :: ncells = 0

        ! Time
        INTEGER(c_intk) :: itstep = 0
        INTEGER(c_intk) :: ittot = 0
        REAL(c_realk) :: timeph = 0.0
        REAL(c_realk) :: dt = 0.0
        REAL(c_realk) :: targetcfl = 0.0

        INTEGER(c_intk) :: realprec = 0
        INTEGER(c_intk) :: intprec = 0
    END TYPE runinfo_t

    PUBLIC :: read_runinfo, write_runinfo

CONTAINS
    SUBROUTINE read_runinfo(ittot, timeph, dt)
        ! Subroutine arguments
        INTEGER(intk), INTENT(out) :: ittot
        REAL(realk), INTENT(out) :: timeph
        REAL(realk), INTENT(out) :: dt

        ! Local variables
        INTEGER(hid_t) :: dtype, file_id, dset_id, filespace, memspace
        INTEGER(hsize_t) :: shape1(1), offset1(1), count1(1)
        INTEGER(int32) :: hdferr
        TYPE(runinfo_t), TARGET :: wdata
        TYPE(C_PTR) :: c_ptr

        IF (ioproc) THEN
            CALL create_runinfo_h5type(dtype)
            CALL get_fileh(file_id)

            ! Menory space
            shape1(1) = 1
            CALL h5screate_simple_f(1, shape1, memspace, hdferr)
            IF (hdferr /= 0) CALL errr(__FILE__, __LINE__)
            IF (ioid == 0) THEN
                CALL h5sselect_all_f(memspace, hdferr)
                IF (hdferr /= 0) CALL errr(__FILE__, __LINE__)
            ELSE
                CALL h5sselect_none_f(memspace, hdferr)
                IF (hdferr /= 0) CALL errr(__FILE__, __LINE__)
            END IF

            CALL hdf5common_dataset_open("RUNINFO", shape1, file_id, dset_id, &
                filespace)

            ! Only rank 0 read data from file
            IF (ioid == 0) THEN
                offset1(1) = shape1(1) - 1
                count1(1) = 1
                CALL h5sselect_hyperslab_f(filespace, H5S_SELECT_SET_F, &
                    offset1, count1, hdferr)
                IF (hdferr /= 0) CALL errr(__FILE__, __LINE__)
            ELSE
                CALL h5sselect_none_f(filespace, hdferr)
                IF (hdferr /= 0) CALL errr(__FILE__, __LINE__)
            END IF

            c_ptr = C_LOC(wdata)
            CALL h5dread_f(dset_id, dtype, c_ptr, hdferr, &
                file_space_id=filespace, mem_space_id=memspace)
            IF (hdferr /= 0) CALL errr(__FILE__, __LINE__)

            CALL hdf5common_dataset_close(dset_id, filespace)

            CALL h5sclose_f(memspace, hdferr)
            IF (hdferr /= 0) CALL errr(__FILE__, __LINE__)

            CALL h5tclose_f(dtype, hdferr)
            IF (hdferr /= 0) CALL errr(__FILE__, __LINE__)

            ittot = wdata%ittot
            timeph = wdata%timeph
            dt = wdata%dt
        END IF

        ! Broadcast to all processes
        CALL MPI_Bcast(ittot, 1, mglet_mpi_int, 0, MPI_COMM_WORLD)
        CALL MPI_Bcast(timeph, 1, mglet_mpi_real, 0, MPI_COMM_WORLD)
        CALL MPI_Bcast(dt, 1, mglet_mpi_real, 0, MPI_COMM_WORLD)
    END SUBROUTINE read_runinfo


    SUBROUTINE write_runinfo(itstep, ittot, timeph, dt, targetcflmax, deltawt, &
            ncells)
        ! Subroutine arguments
        INTEGER(intk), INTENT(in) :: itstep
        INTEGER(intk), INTENT(in) :: ittot
        REAL(realk), INTENT(in) :: timeph
        REAL(realk), INTENT(in) :: dt
        REAL(realk), INTENT(in) :: targetcflmax
        REAL(real64), INTENT(in) :: deltawt
        INTEGER(int64), INTENT(in) :: ncells

        ! Local variables
        INTEGER(hid_t) :: parent_id, dset_id, dtype, crp_list, &
            filespace, memspace
        INTEGER(hsize_t) :: dims1(1), maxdims1(1), offset1(1), count1(1)
        INTEGER(intk) :: starttime(8)
        TYPE(c_ptr) :: f_ptr
        INTEGER(int32) :: hdferr
        LOGICAL :: link_exists
        TYPE(runinfo_t), TARGET :: wdata

        INTEGER :: status
        CHARACTER(len=nchar_programname) :: programname
        CHARACTER(len=nchar_programname), PARAMETER :: default_programname = &
            "mglet"

        IF (.NOT. ioproc) RETURN

        ! Fill runinfo with data for current run
        CALL get_starttime(starttime)
        CALL GET_COMMAND_ARGUMENT(0, programname, status=status)

        ! If status is positive, an error occured
        ! if it is negative, a the 'programname' buffer was too short and the
        ! returnded value was truncated
        ! if it is zero, everything is good
        IF (status <= 0) THEN
            wdata%programname = TRANSFER(programname, wdata%programname)
        ELSE
            wdata%programname = TRANSFER(default_programname, wdata%programname)
        END IF

        wdata%starttime = starttime
        wdata%walltime = deltawt
        wdata%numprocs = numprocs
        wdata%ngrids = ngrid
        wdata%ncells = ncells
        wdata%itstep = itstep
        wdata%ittot = ittot
        wdata%timeph = timeph
        wdata%dt = dt
        wdata%targetcfl = targetcflmax
        wdata%realprec = real_bytes
        wdata%intprec = int_bytes

        CALL create_runinfo_h5type(dtype)

        CALL get_fileh(parent_id)
        CALL h5lexists_f(parent_id, "RUNINFO", link_exists, hdferr)
        IF (hdferr < 0) CALL errr(__FILE__, __LINE__)

        IF (link_exists) THEN
            CALL h5dopen_f(parent_id, "RUNINFO", dset_id, hdferr)
            IF (hdferr /= 0) CALL errr(__FILE__, __LINE__)

            CALL h5dget_space_f(dset_id, filespace, hdferr)
            IF (hdferr /= 0) CALL errr(__FILE__, __LINE__)

            CALL h5sget_simple_extent_dims_f(filespace, dims1, maxdims1, hdferr)
            IF (hdferr < 0) CALL errr(__FILE__, __LINE__)

            CALL h5sclose_f(filespace, hdferr)
            IF (hdferr /= 0) CALL errr(__FILE__, __LINE__)

            ! Extend dataset and re-open for writing
            offset1(1) = dims1(1)
            dims1(1) = dims1(1) + 1
            CALL h5dset_extent_f(dset_id, dims1, hdferr)
            IF (hdferr /= 0) CALL errr(__FILE__, __LINE__)

            CALL h5dget_space_f(dset_id, filespace, hdferr)
            IF (hdferr /= 0) CALL errr(__FILE__, __LINE__)
        ELSE
            maxdims1(1) = H5S_UNLIMITED_F
            dims1(1) = 1
            CALL h5screate_simple_f(1, dims1, filespace, hdferr, maxdims1)
            IF (hdferr /= 0) CALL errr(__FILE__, __LINE__)

            CALL h5pcreate_f(H5P_DATASET_CREATE_F, crp_list, hdferr)
            IF (hdferr /= 0) CALL errr(__FILE__, __LINE__)

            dims1(1) = 1
            CALL h5pset_chunk_f(crp_list, 1, dims1, hdferr)
            IF (hdferr /= 0) CALL errr(__FILE__, __LINE__)

            CALL h5dcreate_f(parent_id, "RUNINFO", dtype, filespace, &
                dset_id, hdferr, dcpl_id=crp_list)
            IF (hdferr /= 0) CALL errr(__FILE__, __LINE__)

            CALL h5pclose_f(crp_list, hdferr)
            IF (hdferr /= 0) CALL errr(__FILE__, __LINE__)

            offset1(1) = 0
        END IF

        ! Memspace
        dims1(1) = 1
        CALL h5screate_simple_f(1, dims1, memspace, hdferr)
        IF (hdferr /= 0) CALL errr(__FILE__, __LINE__)

        ! Only rank 0 write data
        IF (ioid == 0) THEN
            ! Offset is set earlier
            count1(1) = 1
            CALL h5sselect_hyperslab_f(filespace, H5S_SELECT_SET_F, offset1, &
                count1, hdferr)
            IF (hdferr /= 0) CALL errr(__FILE__, __LINE__)

            offset1(1) = 0
            count1(1) = 1
            CALL h5sselect_hyperslab_f(memspace, H5S_SELECT_SET_F, offset1, &
                count1, hdferr)
            IF (hdferr /= 0) CALL errr(__FILE__, __LINE__)
        ELSE
            CALL h5sselect_none_f(filespace, hdferr)
            IF (hdferr /= 0) CALL errr(__FILE__, __LINE__)

            CALL h5sselect_none_f(memspace, hdferr)
            IF (hdferr /= 0) CALL errr(__FILE__, __LINE__)
        END IF

        f_ptr = C_LOC(wdata)
        CALL h5dwrite_f(dset_id, dtype, f_ptr, hdferr, &
            file_space_id=filespace, mem_space_id=memspace)
        IF (hdferr /= 0) CALL errr(__FILE__, __LINE__)

        CALL h5dclose_f(dset_id, hdferr)
        IF (hdferr /= 0) CALL errr(__FILE__, __LINE__)

        CALL h5sclose_f(filespace, hdferr)
        IF (hdferr /= 0) CALL errr(__FILE__, __LINE__)

        CALL h5tclose_f(dtype, hdferr)
        IF (hdferr /= 0) CALL errr(__FILE__, __LINE__)
    END SUBROUTINE write_runinfo


    SUBROUTINE create_runinfo_h5type(dtype)
        ! Subrouitine arguments
        INTEGER(HID_T), INTENT(out) :: dtype

        ! Local variables
        TYPE(runinfo_t), TARGET :: foo
        INTEGER(HID_T) :: strtype, timetype
        INTEGER(int32) :: hdferr
        INTEGER(hsize_t) :: dims1(1)

        ! programname(10)
        CALL h5tcopy_f(H5T_NATIVE_CHARACTER, strtype, hdferr)
        IF (hdferr /= 0) CALL errr(__FILE__, __LINE__)
        CALL h5tset_size_f(strtype, INT(nchar_programname, size_t), hdferr)
        IF (hdferr /= 0) CALL errr(__FILE__, __LINE__)

        ! starttime(8)
        dims1 = 8
        CALL h5tarray_create_f(mglet_hdf5_int, 1, dims1, &
            timetype, hdferr)
        IF (hdferr /= 0) CALL errr(__FILE__, __LINE__)

        ! Compound datatype
        CALL h5tcreate_f(H5T_COMPOUND_F, C_SIZEOF(foo), dtype, hdferr)
        IF (hdferr /= 0) CALL errr(__FILE__, __LINE__)

        CALL h5tinsert_f(dtype, "PROGRAMNAME", &
            H5OFFSETOF(C_LOC(foo), C_LOC(foo%programname)), strtype, hdferr)
        IF (hdferr /= 0) CALL errr(__FILE__, __LINE__)

        CALL h5tinsert_f(dtype, "STARTTIME", &
            H5OFFSETOF(C_LOC(foo), C_LOC(foo%starttime)), timetype, hdferr)
        IF (hdferr /= 0) CALL errr(__FILE__, __LINE__)

        CALL h5tinsert_f(dtype, "WALLTIME", &
            H5OFFSETOF(C_LOC(foo), C_LOC(foo%walltime)), &
            h5kind_to_type(c_double, H5_REAL_KIND), hdferr)
        IF (hdferr /= 0) CALL errr(__FILE__, __LINE__)

        CALL h5tinsert_f(dtype, "NUMPROCS", &
            H5OFFSETOF(C_LOC(foo), C_LOC(foo%numprocs)), &
            mglet_hdf5_int, hdferr)
        IF (hdferr /= 0) CALL errr(__FILE__, __LINE__)

        CALL h5tinsert_f(dtype, "NGRIDS", &
            H5OFFSETOF(C_LOC(foo), C_LOC(foo%ngrids)), &
            mglet_hdf5_int, hdferr)
        IF (hdferr /= 0) CALL errr(__FILE__, __LINE__)

        CALL h5tinsert_f(dtype, "NCELLS", &
            H5OFFSETOF(C_LOC(foo), C_LOC(foo%ncells)), &
            h5kind_to_type(c_long_long, H5_INTEGER_KIND), hdferr)
        IF (hdferr /= 0) CALL errr(__FILE__, __LINE__)

        CALL h5tinsert_f(dtype, "ITSTEP", &
            H5OFFSETOF(C_LOC(foo), C_LOC(foo%itstep)), &
            mglet_hdf5_int, hdferr)
        IF (hdferr /= 0) CALL errr(__FILE__, __LINE__)

        CALL h5tinsert_f(dtype, "ITTOT", &
            H5OFFSETOF(C_LOC(foo), C_LOC(foo%ittot)), &
            mglet_hdf5_int, hdferr)
        IF (hdferr /= 0) CALL errr(__FILE__, __LINE__)

        CALL h5tinsert_f(dtype, "TIMEPH", &
            H5OFFSETOF(C_LOC(foo), C_LOC(foo%timeph)), &
            mglet_hdf5_real, hdferr)
        IF (hdferr /= 0) CALL errr(__FILE__, __LINE__)

        CALL h5tinsert_f(dtype, "DT", &
            H5OFFSETOF(C_LOC(foo), C_LOC(foo%dt)), &
            mglet_hdf5_real, hdferr)
        IF (hdferr /= 0) CALL errr(__FILE__, __LINE__)

        CALL h5tinsert_f(dtype, "TARGETCFL", &
            H5OFFSETOF(C_LOC(foo), C_LOC(foo%targetcfl)), &
            mglet_hdf5_real, hdferr)
        IF (hdferr /= 0) CALL errr(__FILE__, __LINE__)

        CALL h5tinsert_f(dtype, "REALPREC", &
            H5OFFSETOF(C_LOC(foo), C_LOC(foo%realprec)), &
            mglet_hdf5_int, hdferr)
        IF (hdferr /= 0) CALL errr(__FILE__, __LINE__)

        CALL h5tinsert_f(dtype, "INTPREC", &
            H5OFFSETOF(C_LOC(foo), C_LOC(foo%intprec)), &
            mglet_hdf5_int, hdferr)
        IF (hdferr /= 0) CALL errr(__FILE__, __LINE__)

        ! close datatypes no longer required
        CALL h5tclose_f(strtype, hdferr)
        IF (hdferr /= 0) CALL errr(__FILE__, __LINE__)

        CALL h5tclose_f(timetype, hdferr)
        IF (hdferr /= 0) CALL errr(__FILE__, __LINE__)
    END SUBROUTINE create_runinfo_h5type

END MODULE runinfo_mod
