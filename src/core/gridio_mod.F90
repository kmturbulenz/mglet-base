MODULE gridio_mod
    USE MPI_f08
    USE HDF5

    USE comms_mod
    USE err_mod
    USE hdf5common_mod
    USE precision_mod

    IMPLICIT NONE (type, external)
    PRIVATE

    ! Gridinfo table type
    INTEGER(c_intk), PARAMETER :: gridinfo_elems = 13
    TYPE, BIND(C) :: gridinfo_t
        INTEGER(c_intk) :: igrid

        INTEGER(c_intk) :: ii
        INTEGER(c_intk) :: jj
        INTEGER(c_intk) :: kk

        INTEGER(c_intk) :: level
        INTEGER(c_intk) :: iparent
        INTEGER(c_intk) :: iposition
        INTEGER(c_intk) :: jposition
        INTEGER(c_intk) :: kposition

        ! General purpose "flags"
        INTEGER(c_intk) :: flags

        ! Position in level-datasets
        INTEGER(c_intk) :: igridlvl

        ! Bounding box of grid (x0, x1, y0, y1, z0, z1)
        REAL(c_realk) :: bbox(6)

        ! Neighbouring grid
        INTEGER(c_intk) :: nbrgrid(26)
    END TYPE gridinfo_t

    ! Boundary conditions table type
    INTEGER(c_intk), PARAMETER :: maxboconds = 2
    INTEGER(c_intk), PARAMETER :: nchar = 8
    INTEGER(c_intk), PARAMETER :: bcond_elems = 4
    TYPE, BIND(C) :: bcond_t
        ! Actual maximum number of boundary conditions set for any single face
        INTEGER(c_intk) :: nbocd

        ! Character sequence storing boundary condition identifier
        CHARACTER(kind=C_CHAR) :: type(nchar, maxboconds)

        ! Offset and number of integer and real parameters to boundary condition
        INTEGER(c_intk) :: intprm(2, maxboconds)
        INTEGER(c_intk) :: realprm(2, maxboconds)
    END TYPE bcond_t

    PUBLIC :: read_gridinfo, read_bcondinfo, maxboconds, gridinfo_t, bcond_t, &
        write_gridinfo, write_bcondinfo

CONTAINS
    SUBROUTINE read_gridinfo(parent_id, gridinfo, realprms, intprms, ngrid)
        USE, INTRINSIC :: ISO_C_BINDING, ONLY: C_LOC, C_PTR, C_NULL_PTR

        ! Subroutine arguments
        INTEGER(HID_T), INTENT(in) :: parent_id
        TYPE(gridinfo_t), ALLOCATABLE, INTENT(inout), TARGET :: gridinfo(:)
        REAL(realk), ALLOCATABLE, INTENT(inout) :: realprms(:)
        INTEGER(intk), ALLOCATABLE, INTENT(inout) :: intprms(:)
        INTEGER(intk), INTENT(out) :: ngrid

        ! Local variables
        INTEGER(HID_T) :: dset_id, filespace_id, memspace_id, gridinfo_h5type
        INTEGER(HSIZE_T) :: shape(1), dims(1)
        INTEGER(int32) :: ierr
        TYPE(MPI_datatype) :: gridinfo_mpitype
        TYPE(C_PTR) :: ptr
        LOGICAL :: attribute_exists
        INTEGER(intk) :: nrealprm, nintprm

        ! First we must figure out how many grids there are, i.e. the shape
        ! of the 'GRIDINFO'-table
        !
        ! This opens the dataset and return the filespace along with the shape
        CALL hdf5common_dataset_open('GRIDINFO', shape, parent_id, dset_id, &
            filespace_id)
        ngrid = INT(shape(1), intk)

        ! Only IO processes opened the dataset and knows the size, bcast this
        ! to everyone
        CALL MPI_Bcast(ngrid, 1, mglet_mpi_int, 0, MPI_COMM_WORLD)

        ! Every process ends up with a copy of the 'GRIDINFO' table, so
        ! everyone allocate this
        ALLOCATE(gridinfo(ngrid))

        ! Read GRIDINFO
        IF (ioproc) THEN
            CALL create_gridinfo_h5type(gridinfo_h5type)

            dims(1) = ngrid
            CALL h5screate_simple_f(1, dims, memspace_id, ierr)
            IF (ierr /= 0) CALL errr(__FILE__, __LINE__)

            ! Let rank 0 do all the reading
            ptr = C_LOC(gridinfo(1))
            IF (myid /= 0) THEN
                ptr = C_NULL_PTR

                CALL h5sselect_none_f(memspace_id, ierr)
                IF (ierr /= 0) CALL errr(__FILE__, __LINE__)

                CALL h5sselect_none_f(filespace_id, ierr)
                IF (ierr /= 0) CALL errr(__FILE__, __LINE__)
            END IF

            CALL h5dread_f(dset_id, gridinfo_h5type, ptr, ierr, &
                mem_space_id=memspace_id, file_space_id=filespace_id)
            IF (ierr /= 0) CALL errr(__FILE__, __LINE__)

            CALL h5sclose_f(memspace_id, ierr)
            IF (ierr /= 0) CALL errr(__FILE__, __LINE__)

            CALL h5tclose_f(gridinfo_h5type, ierr)
            IF (ierr /= 0) CALL errr(__FILE__, __LINE__)
        END IF

        ! Reading attributes "REALPRM" and "INTPRM" of GRIDINFO table
        nrealprm = 0
        attribute_exists = .FALSE.
        CALL hdf5common_attr_exists("REALPRMS", attribute_exists, dset_id, &
            shape=shape)
        IF (attribute_exists) THEN
            nrealprm = INT(shape(1), intk)
        END IF
        CALL MPI_Bcast(nrealprm, 1, mglet_mpi_int, 0, MPI_COMM_WORLD)
        ALLOCATE(realprms(nrealprm))

        IF (nrealprm > 0) THEN
            CALL hdf5common_attr_read_arr("REALPRMS", realprms, dset_id)
            CALL MPI_Bcast(realprms, nrealprm, mglet_mpi_real, 0, &
                MPI_COMM_WORLD)
        END IF

        nintprm = 0
        attribute_exists = .FALSE.
        CALL hdf5common_attr_exists("INTPRMS", attribute_exists, dset_id, &
            shape=shape)
        IF (attribute_exists) THEN
            nintprm = INT(shape(1), intk)
        END IF
        CALL MPI_Bcast(nintprm, 1, mglet_mpi_int, 0, MPI_COMM_WORLD)
        ALLOCATE(intprms(nintprm))

        IF (nintprm > 0) THEN
            CALL hdf5common_attr_read_arr("INTPRMS", intprms, dset_id)
            CALL MPI_Bcast(intprms, nintprm, mglet_mpi_int, 0, &
                MPI_COMM_WORLD)
        END IF

        ! Close dataset and filespace
        CALL hdf5common_dataset_close(dset_id, filespace_id)

        ! Broadcast the resulting gridinfo-table
        CALL create_gridinfo_mpitype(gridinfo_mpitype)
        CALL MPI_Bcast(gridinfo, INT(ngrid, int32), gridinfo_mpitype, 0, &
            MPI_COMM_WORLD)
        CALL MPI_Type_free(gridinfo_mpitype)
    END SUBROUTINE read_gridinfo


    SUBROUTINE read_bcondinfo(parent_id, face, bcond_arr)
        USE, INTRINSIC :: ISO_C_BINDING, ONLY: C_LOC, C_PTR, C_NULL_PTR

        ! Subroutine arguments
        INTEGER(HID_T), INTENT(in) :: parent_id
        CHARACTER(len=*), INTENT(in) :: face
        TYPE(bcond_t), INTENT(inout), TARGET :: bcond_arr(:)

        ! Local variables
        INTEGER(HID_T) :: dset_id, filespace_id, memspace_id, bcond_h5type
        INTEGER(HSIZE_T) :: dims(1), shape(1)
        INTEGER(intk) :: ngrid
        INTEGER(int32) :: ierr
        TYPE(MPI_datatype) :: bcond_mpitype
        TYPE(C_PTR) :: ptr

        ! Ngrid is set from size of array
        ngrid = SIZE(bcond_arr)

        ! This opens the dataset and return the filespace along with the shape
        CALL hdf5common_dataset_open(face, shape, parent_id, dset_id, &
            filespace_id)

        ! Sanity check
        IF (ioproc .AND. shape(1) /= ngrid) THEN
            WRITE(*, *) "face: ", face
            WRITE(*, *) "shape(1): ", shape(1)
            WRITE(*, *) "ngrid: ", ngrid
            CALL errr(__FILE__, __LINE__)
        END IF

        ! Read GRIDINFO
        IF (ioproc) THEN
            CALL create_bcond_h5type(bcond_h5type)

            dims(1) = ngrid
            CALL h5screate_simple_f(1, dims, memspace_id, ierr)
            IF (ierr /= 0) CALL errr(__FILE__, __LINE__)

            ! Let rank 0 do all the reading
            ptr = C_LOC(bcond_arr(1))
            IF (myid /= 0) THEN
                ptr = C_NULL_PTR

                CALL h5sselect_none_f(memspace_id, ierr)
                IF (ierr /= 0) CALL errr(__FILE__, __LINE__)

                CALL h5sselect_none_f(filespace_id, ierr)
                IF (ierr /= 0) CALL errr(__FILE__, __LINE__)
            END IF

            CALL h5dread_f(dset_id, bcond_h5type, ptr, ierr, &
                mem_space_id=memspace_id, file_space_id=filespace_id)
            IF (ierr /= 0) CALL errr(__FILE__, __LINE__)

            CALL h5dclose_f(dset_id, ierr)
            IF (ierr /= 0) CALL errr(__FILE__, __LINE__)

            CALL h5sclose_f(memspace_id, ierr)
            IF (ierr /= 0) CALL errr(__FILE__, __LINE__)

            CALL h5sclose_f(filespace_id, ierr)
            IF (ierr /= 0) CALL errr(__FILE__, __LINE__)

            CALL h5tclose_f(bcond_h5type, ierr)
            IF (ierr /= 0) CALL errr(__FILE__, __LINE__)
        END IF

        ! Broadcast the resulting gridinfo-table
        CALL create_bcond_mpitype(bcond_mpitype)
        CALL MPI_Bcast(bcond_arr, INT(ngrid, int32), bcond_mpitype, 0, &
            MPI_COMM_WORLD)
        CALL MPI_Type_free(bcond_mpitype)
    END SUBROUTINE read_bcondinfo


    SUBROUTINE create_gridinfo_h5type(dtype)
        ! Subrouitine arguments
        INTEGER(HID_T), INTENT(out) :: dtype

        ! Local variables
        TYPE(gridinfo_t), TARGET :: wdata
        INTEGER(HID_T) :: bboxtype, nbrtype
        INTEGER(hsize_t) :: dims1(1)
        INTEGER(int32) :: ierr

        ! BBOX datatype
        dims1 = 6
        CALL h5tarray_create_f(mglet_hdf5_real, 1, dims1, bboxtype, ierr)
        IF (ierr /= 0) CALL errr(__FILE__, __LINE__)

        ! Neighbour datatype
        dims1 = 26
        CALL h5tarray_create_f(mglet_hdf5_int, 1, dims1, nbrtype, ierr)
        IF (ierr /= 0) CALL errr(__FILE__, __LINE__)

        ! Gridinfo-table as compound data type
        CALL h5tcreate_f(H5T_COMPOUND_F, C_SIZEOF(wdata), dtype, ierr)
        IF (ierr /= 0) CALL errr(__FILE__, __LINE__)

        CALL h5tinsert_f(dtype, "IGRID", &
            H5OFFSETOF(C_LOC(wdata), C_LOC(wdata%igrid)), &
            mglet_hdf5_int, ierr)
        IF (ierr /= 0) CALL errr(__FILE__, __LINE__)

        CALL h5tinsert_f(dtype, "II", &
            H5OFFSETOF(C_LOC(wdata), C_LOC(wdata%ii)), &
            mglet_hdf5_int, ierr)
        IF (ierr /= 0) CALL errr(__FILE__, __LINE__)

        CALL h5tinsert_f(dtype, "JJ", &
            H5OFFSETOF(C_LOC(wdata), C_LOC(wdata%jj)), &
            mglet_hdf5_int, ierr)
        IF (ierr /= 0) CALL errr(__FILE__, __LINE__)

        CALL h5tinsert_f(dtype, "KK", &
            H5OFFSETOF(C_LOC(wdata), C_LOC(wdata%kk)), &
            mglet_hdf5_int, ierr)
        IF (ierr /= 0) CALL errr(__FILE__, __LINE__)

        CALL h5tinsert_f(dtype, "LEVEL", &
            H5OFFSETOF(C_LOC(wdata), C_LOC(wdata%level)), &
            mglet_hdf5_int, ierr)
        IF (ierr /= 0) CALL errr(__FILE__, __LINE__)

        CALL h5tinsert_f(dtype, "IPARENT", &
            H5OFFSETOF(C_LOC(wdata), C_LOC(wdata%iparent)), &
            mglet_hdf5_int, ierr)
        IF (ierr /= 0) CALL errr(__FILE__, __LINE__)

        CALL h5tinsert_f(dtype, "IPOSITION", &
            H5OFFSETOF(C_LOC(wdata), C_LOC(wdata%iposition)), &
            mglet_hdf5_int, ierr)
        IF (ierr /= 0) CALL errr(__FILE__, __LINE__)

        CALL h5tinsert_f(dtype, "JPOSITION", &
            H5OFFSETOF(C_LOC(wdata), C_LOC(wdata%jposition)), &
            mglet_hdf5_int, ierr)
        IF (ierr /= 0) CALL errr(__FILE__, __LINE__)

        CALL h5tinsert_f(dtype, "KPOSITION", &
            H5OFFSETOF(C_LOC(wdata), C_LOC(wdata%kposition)), &
            mglet_hdf5_int, ierr)
        IF (ierr /= 0) CALL errr(__FILE__, __LINE__)

        CALL h5tinsert_f(dtype, "FLAGS", &
            H5OFFSETOF(C_LOC(wdata), C_LOC(wdata%flags)), &
            mglet_hdf5_int, ierr)
        IF (ierr /= 0) CALL errr(__FILE__, __LINE__)

        CALL h5tinsert_f(dtype, "IGRIDLVL", &
            H5OFFSETOF(C_LOC(wdata), C_LOC(wdata%igridlvl)), &
            mglet_hdf5_int, ierr)
        IF (ierr /= 0) CALL errr(__FILE__, __LINE__)

        CALL h5tinsert_f(dtype, "BBOX", &
            H5OFFSETOF(C_LOC(wdata), C_LOC(wdata%bbox)), &
            bboxtype, ierr)
        IF (ierr /= 0) CALL errr(__FILE__, __LINE__)

        CALL h5tinsert_f(dtype, "NBRGRID", &
            H5OFFSETOF(C_LOC(wdata), C_LOC(wdata%nbrgrid)), &
            nbrtype, ierr)
        IF (ierr /= 0) CALL errr(__FILE__, __LINE__)

        ! Close types not needed any more
        CALL h5tclose_f(bboxtype, ierr)
        IF (ierr /= 0) CALL errr(__FILE__, __LINE__)
        CALL h5tclose_f(nbrtype, ierr)
        IF (ierr /= 0) CALL errr(__FILE__, __LINE__)
    END SUBROUTINE create_gridinfo_h5type


    SUBROUTINE create_gridinfo_mpitype(dtype)
        ! Subrouitine arguments
        TYPE(MPI_Datatype), INTENT(out) :: dtype

        ! Local variables
        INTEGER(intk) :: i
        TYPE(gridinfo_t) :: foo
        INTEGER(MPI_ADDRESS_KIND) :: base, disp(gridinfo_elems)
        INTEGER(int32) :: blocklen(gridinfo_elems)
        TYPE(MPI_Datatype) :: types(gridinfo_elems)
        TYPE(MPI_Datatype) :: bboxtype, nbrtype

        CALL MPI_Type_contiguous(6, mglet_mpi_real, bboxtype)
        CALL MPI_Type_contiguous(26, mglet_mpi_int, nbrtype)

        CALL MPI_Get_address(foo%igrid, disp(1))
        CALL MPI_Get_address(foo%ii, disp(2))
        CALL MPI_Get_address(foo%jj, disp(3))
        CALL MPI_Get_address(foo%kk, disp(4))
        CALL MPI_Get_address(foo%level, disp(5))
        CALL MPI_Get_address(foo%iparent, disp(6))
        CALL MPI_Get_address(foo%iposition, disp(7))
        CALL MPI_Get_address(foo%jposition, disp(8))
        CALL MPI_Get_address(foo%kposition, disp(9))
        CALL MPI_Get_address(foo%flags, disp(10))
        CALL MPI_Get_address(foo%igridlvl, disp(11))
        CALL MPI_Get_address(foo%bbox, disp(12))
        CALL MPI_Get_address(foo%nbrgrid, disp(13))

        types(1) = mglet_mpi_int    ! igrid
        types(2) = mglet_mpi_int    ! ii
        types(3) = mglet_mpi_int    ! jj
        types(4) = mglet_mpi_int    ! kk
        types(5) = mglet_mpi_int    ! level
        types(6) = mglet_mpi_int    ! iparent
        types(7) = mglet_mpi_int    ! iposition
        types(8) = mglet_mpi_int    ! jposition
        types(9) = mglet_mpi_int    ! kposition
        types(10) = mglet_mpi_int   ! flags
        types(11) = mglet_mpi_int   ! igridlvl
        types(12) = bboxtype        ! bbox
        types(13) = nbrtype         ! nbrgrid

        base = disp(1)
        DO i = 1, gridinfo_elems
            disp(i) = disp(i) - base
        END DO

        blocklen = 1
        CALL MPI_Type_create_struct(gridinfo_elems, blocklen, disp, types, &
            dtype)
        CALL MPI_Type_commit(dtype)

        ! Free types used temporarily
        !   OK according to:
        !   https://www.open-mpi.org/doc/v4.0/man3/MPI_Type_free.3.php
        CALL MPI_Type_free(bboxtype)
        CALL MPI_Type_free(nbrtype)
    END SUBROUTINE create_gridinfo_mpitype


    SUBROUTINE create_bcond_h5type(dtype)
        ! Subrouitine arguments
        INTEGER(HID_T), INTENT(out) :: dtype

        ! Local variables
        TYPE(bcond_t), TARGET :: wdata
        INTEGER(HID_T) :: str_t, strarr_t, intarr_t
        INTEGER(int32) :: ierr
        INTEGER(hsize_t) :: dims1(1), dims2(2)

        ! Create type for CHARACTER(C_CHAR) :: type(nchar, maxboconds)
        CALL h5tcopy_f(H5T_NATIVE_CHARACTER, str_t, ierr)
        IF (ierr /= 0) CALL errr(__FILE__, __LINE__)
        CALL h5tset_size_f(str_t, INT(nchar, size_t), ierr)
        IF (ierr /= 0) CALL errr(__FILE__, __LINE__)

        dims1 = maxboconds
        CALL h5tarray_create_f(str_t, 1, dims1, strarr_t, ierr)
        IF (ierr /= 0) CALL errr(__FILE__, __LINE__)
        CALL h5tclose_f(str_t, ierr)
        IF (ierr /= 0) CALL errr(__FILE__, __LINE__)

        ! Create arrays for boundary condition parameters
        dims2(1) = 2
        dims2(2) = maxboconds
        CALL h5tarray_create_f(mglet_hdf5_int, 2, dims2, intarr_t, ierr)
        IF (ierr /= 0) CALL errr(__FILE__, __LINE__)

        ! Bcond-table as compound data type
        CALL h5tcreate_f(H5T_COMPOUND_F, C_SIZEOF(wdata), dtype, ierr)
        IF (ierr /= 0) CALL errr(__FILE__, __LINE__)

        CALL h5tinsert_f(dtype, "NBOCD", &
            H5OFFSETOF(C_LOC(wdata), C_LOC(wdata%nbocd)), &
            mglet_hdf5_int, ierr)
        IF (ierr /= 0) CALL errr(__FILE__, __LINE__)

        CALL h5tinsert_f(dtype, "TYPE", &
            H5OFFSETOF(C_LOC(wdata), C_LOC(wdata%type)), &
            strarr_t, ierr)
        IF (ierr /= 0) CALL errr(__FILE__, __LINE__)

        CALL h5tinsert_f(dtype, "INTPRM", &
            H5OFFSETOF(C_LOC(wdata), C_LOC(wdata%intprm)), &
            intarr_t, ierr)
        IF (ierr /= 0) CALL errr(__FILE__, __LINE__)

        CALL h5tinsert_f(dtype, "REALPRM", &
            H5OFFSETOF(C_LOC(wdata), C_LOC(wdata%realprm)), &
            intarr_t, ierr)
        IF (ierr /= 0) CALL errr(__FILE__, __LINE__)

        ! Close temporary datatypes
        CALL h5tclose_f(strarr_t, ierr)
        IF (ierr /= 0) CALL errr(__FILE__, __LINE__)
        CALL h5tclose_f(intarr_t, ierr)
        IF (ierr /= 0) CALL errr(__FILE__, __LINE__)
    END SUBROUTINE create_bcond_h5type


    SUBROUTINE create_bcond_mpitype(dtype)
        ! Subrouitine arguments
        TYPE(MPI_Datatype), INTENT(out) :: dtype

        ! Local variables
        INTEGER(intk) :: i
        TYPE(bcond_t) :: foo
        INTEGER(MPI_ADDRESS_KIND) :: base, disp(bcond_elems)
        INTEGER(int32) :: blocklen(bcond_elems)
        TYPE(MPI_Datatype) :: types(bcond_elems)
        TYPE(MPI_Datatype) :: typetype, prmtype

        CALL MPI_Type_contiguous(nchar*maxboconds, mpi_byte, typetype)
        CALL MPI_Type_contiguous(2*maxboconds, mglet_mpi_int, prmtype)

        CALL MPI_Get_address(foo%nbocd, disp(1))
        CALL MPI_Get_address(foo%type(1, 1), disp(2))
        CALL MPI_Get_address(foo%intprm(1, 1), disp(3))
        CALL MPI_Get_address(foo%realprm(1, 1), disp(4))

        types(1) = mglet_mpi_int    ! nbocd
        types(2) = typetype         ! type
        types(3) = prmtype          ! intprm
        types(4) = prmtype          ! realprm

        base = disp(1)
        DO i = 1, bcond_elems
            disp(i) = disp(i) - base
        END DO

        blocklen = 1
        CALL MPI_Type_create_struct(bcond_elems, blocklen, disp, types, &
            dtype)
        CALL MPI_Type_commit(dtype)

        ! Free types used temporarily
        !   OK according to:
        !   https://www.open-mpi.org/doc/v4.0/man3/MPI_Type_free.3.php
        CALL MPI_Type_free(typetype)
        CALL MPI_Type_free(prmtype)
    END SUBROUTINE create_bcond_mpitype


    SUBROUTINE write_gridinfo(parent_id, gridinfo, realprms, intprms, new_ngrid)
        ! Subroutine arguments
        INTEGER(HID_T), INTENT(inout) :: parent_id
        TYPE(gridinfo_t), INTENT(in), TARGET :: gridinfo(:)
        REAL(realk), ALLOCATABLE, INTENT(in) :: realprms(:)
        INTEGER(intk), ALLOCATABLE, INTENT(in) :: intprms(:)
        INTEGER(intk), INTENT(in) :: new_ngrid

        ! Local variables
        INTEGER(int32) :: ierr
        TYPE(C_PTR) :: ptr
        INTEGER(HID_T) :: dset_id, dspace_id, gridinfo_h5type
        INTEGER(HSIZE_T) :: dims(1)

        CALL create_gridinfo_h5type(gridinfo_h5type)

        dims(1) = new_ngrid
        CALL h5screate_simple_f(1, dims, dspace_id, ierr)
        IF (ierr /= 0) CALL errr(__FILE__, __LINE__)

        CALL h5dcreate_f(parent_id, "GRIDINFO", gridinfo_h5type, &
            dspace_id, dset_id, ierr)
        IF (ierr /= 0) CALL errr(__FILE__, __LINE__)

        ! Only rank 0 have valid data to write
        IF (myid == 0) THEN
            ptr = C_LOC(gridinfo)
        ELSE
            ptr = C_NULL_PTR
            CALL h5sselect_none_f(dspace_id, ierr)
            IF (ierr /= 0) CALL errr(__FILE__, __LINE__)
        END IF

        CALL h5dwrite_f(dset_id, gridinfo_h5type, ptr, ierr, &
            file_space_id=dspace_id)
        IF (ierr /= 0) CALL errr(__FILE__, __LINE__)

        CALL h5dclose_f(dset_id, ierr)
        IF (ierr /= 0) CALL errr(__FILE__, __LINE__)

        CALL h5sclose_f(dspace_id, ierr)
        IF (ierr /= 0) CALL errr(__FILE__, __LINE__)

        CALL h5tclose_f(gridinfo_h5type, ierr)
        IF (ierr /= 0) CALL errr(__FILE__, __LINE__)

        ! Writing attributes "REALPRM" and "INTPRM" of GRIDINFO table
        ! All processes has identical copies of these for simplicity...
        IF (SIZE(realprms) > 0) THEN
            CALL hdf5common_attr_write_arr("REALPRMS", realprms, &
                parent_id, "GRIDINFO")
        END IF

        IF (SIZE(intprms) > 0) THEN
            CALL hdf5common_attr_write_arr("INTPRMS", intprms, &
                parent_id, "GRIDINFO")
        END IF
    END SUBROUTINE write_gridinfo


    SUBROUTINE write_bcondinfo(parent_id, bcond_arr, face, new_ngrid)
        ! Subroutine arguments
        INTEGER(HID_T), INTENT(inout) :: parent_id
        TYPE(bcond_t), INTENT(IN), TARGET :: bcond_arr(:)
        CHARACTER(LEN=*), INTENT(IN) :: face
        INTEGER(intk), INTENT(in) :: new_ngrid

        ! Local variables
        INTEGER(int32) :: ierr
        TYPE(C_PTR) :: ptr
        INTEGER(HID_T) :: dset_id, dspace_id, bcond_h5type
        INTEGER(HSIZE_T) :: dims(1)

        CALL create_bcond_h5type(bcond_h5type)

        dims(1) = new_ngrid
        CALL h5screate_simple_f(1, dims, dspace_id, ierr)
        IF (ierr /= 0) CALL errr(__FILE__, __LINE__)

        CALL h5dcreate_f(parent_id, face, bcond_h5type, &
            dspace_id, dset_id, ierr)
        IF (ierr /= 0) CALL errr(__FILE__, __LINE__)

        ! Only rank 0 have valid data to write
        IF (myid == 0) THEN
            ptr = C_LOC(bcond_arr)
        ELSE
            ptr = C_NULL_PTR
            CALL h5sselect_none_f(dspace_id, ierr)
            IF (ierr /= 0) CALL errr(__FILE__, __LINE__)
        END IF

        CALL h5dwrite_f(dset_id, bcond_h5type, ptr, ierr, &
            file_space_id=dspace_id)
        IF (ierr /= 0) CALL errr(__FILE__, __LINE__)

        CALL h5dclose_f(dset_id, ierr)
        IF (ierr /= 0) CALL errr(__FILE__, __LINE__)

        CALL h5sclose_f(dspace_id, ierr)
        IF (ierr /= 0) CALL errr(__FILE__, __LINE__)

        CALL h5tclose_f(bcond_h5type, ierr)
        IF (ierr /= 0) CALL errr(__FILE__, __LINE__)
    END SUBROUTINE write_bcondinfo
END MODULE gridio_mod
