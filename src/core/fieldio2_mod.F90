MODULE fieldio2_mod
    USE MPI_f08
    USE HDF5

    USE commbuf_mod
    USE comms_mod
    USE connect2_mod
    USE err_mod, ONLY: errr
    USE field_mod
    USE grids_mod
    USE hdf5common_mod
    USE precision_mod
    USE qsort_mod
    USE stencilio_mod
    USE timer_mod

    USE, INTRINSIC :: ISO_C_BINDING

    IMPLICIT NONE (type, external)
    PRIVATE

    LOGICAL :: is_init = .FALSE.

    ! Number of columns in the 'mygridinfo' table. Need to know this globally
    ! for the counts, displs that are pre-computed for communication
    INTEGER(intk), PARAMETER :: info_dim = 2

    ! Properties for the IO group (iogrcomm comm.) - this is all processes that
    ! are associated with a particular IO process. First stage gathering of
    ! data over MPI is done with these structures
    INTEGER(intk) :: niogrgrids
    INTEGER(int32), ALLOCATABLE :: iogr_counts(:)
    INTEGER(int32), ALLOCATABLE :: iogr_displs(:)

    ! Properties of the IO communicator (iocomm)
    ! All the IO processes that do IO use these structures to correctly
    ! communicate and write data to disk.
    INTEGER(int32), ALLOCATABLE :: io_counts(:)
    INTEGER(int32), ALLOCATABLE :: io_displs(:)

    ! Type containing an offset,length pair used to look up grids
    TYPE, BIND(C) :: offset_t
        INTEGER(c_long_long) :: offset
        INTEGER(c_long_long) :: length
    END TYPE offset_t

    ! Types for reading in attributes for restart-fields
    TYPE, BIND(C) :: rattr_t
        CHARACTER(kind=c_char, len=1) :: name(nchar_name)
        REAL(c_realk) :: value
    END TYPE rattr_t

    TYPE, BIND(C) :: iattr_t
        CHARACTER(kind=c_char, len=1) :: name(nchar_name)
        INTEGER(c_intk) :: value
    END TYPE iattr_t

    INTEGER(hid_t) :: offset_h5t

    TYPE(MPI_Datatype) :: offset_mpit
    TYPE(MPI_Datatype) :: rattr_mpit
    TYPE(MPI_Datatype) :: iattr_mpit

    ! module interfaces
    PUBLIC :: init_fieldio, finish_fieldio, fieldio_write, fieldio_read

CONTAINS
    ! Initialize common data structures
    SUBROUTINE init_fieldio()
        ! Local variables
        INTEGER(intk) :: i

        CALL set_timer(100, "FIELDIO2")
        CALL set_timer(101, "FIELDIO2_GATHER")
        CALL set_timer(102, "FIELDIO2_WRITE")
        CALL set_timer(103, "FIELDIO2_READ")
        CALL set_timer(104, "FIELDIO2_SCATTER")

        IF (ioproc) THEN
            ALLOCATE(iogr_counts(iogrprocs))
            ALLOCATE(iogr_displs(iogrprocs))

            ALLOCATE(io_counts(ioprocs))
            ALLOCATE(io_displs(ioprocs))
        ELSE
            ALLOCATE(iogr_counts(0))
            ALLOCATE(iogr_displs(0))

            ALLOCATE(io_counts(0))
            ALLOCATE(io_displs(0))
        END IF

        ! Establish properties for IO group communicator
        BLOCK
            INTEGER(intk), ALLOCATABLE :: niogridsproc(:)

            IF (ioproc) THEN
                ALLOCATE(niogridsproc(iogrprocs))
            ELSE
                ALLOCATE(niogridsproc(0))
            END IF

            CALL MPI_Gather(nmygrids, 1, mglet_mpi_int, &
                niogridsproc, 1, mglet_mpi_int, 0, iogrcomm)

            iogr_counts = 0
            iogr_displs = 0
            IF (ioproc) THEN
                DO i = 1, iogrprocs
                    iogr_counts(i) = info_dim*niogridsproc(i)
                END DO
                DO i = 2, iogrprocs
                    iogr_displs(i) = iogr_displs(i-1) + iogr_counts(i-1)
                END DO
            END IF

            niogrgrids = SUM(niogridsproc)
            DEALLOCATE(niogridsproc)
        END BLOCK

        ! Establish properties for IO communicator
        init_iocomm: BLOCK
            INTEGER(intk) :: niogrids_tot
            INTEGER(intk), ALLOCATABLE :: niogrids_all(:)

            ! Non-io processes do not need to execute code in this block
            IF (.NOT. ioproc) EXIT init_iocomm

            ALLOCATE(niogrids_all(ioprocs))
            CALL MPI_Allgather(niogrgrids, 1, mglet_mpi_int, &
                niogrids_all, 1, mglet_mpi_int, iocomm)

            ! At this stage we do a sanity check that niogrids_tot == ngrid
            ! Even for fields that are not defined on all levels, we should have
            ! ngrid == niogrids_tot. Grids on levels where no data are present
            ! just have length == 0
            niogrids_tot = SUM(niogrids_all)
            IF (ngrid /= niogrids_tot) CALL errr(__FILE__, __LINE__)

            io_counts = 0
            io_displs = 0
            DO i = 1, ioprocs
                io_counts(i) = info_dim*niogrids_all(i)
            END DO
            DO i = 2, ioprocs
                io_displs(i) = io_displs(i-1) + io_counts(i-1)
            END DO

            DEALLOCATE(niogrids_all)
        END BLOCK init_iocomm

        ! Create HDF5 type for offset,count
        BLOCK
            ! Local variables
            TYPE(offset_t), TARGET :: foo
            INTEGER(int32) :: hdferr
            INTEGER(hid_t) :: int64_h5t

            int64_h5t = h5kind_to_type(c_long_long, H5_INTEGER_KIND)

            CALL h5tcreate_f(H5T_COMPOUND_F, C_SIZEOF(foo), offset_h5t, hdferr)
            IF (hdferr /= 0) CALL errr(__FILE__, __LINE__)

            CALL h5tinsert_f(offset_h5t, "OFFSET", &
                 H5OFFSETOF(C_LOC(foo), C_LOC(foo%offset)), int64_h5t, hdferr)
            IF (hdferr /= 0) CALL errr(__FILE__, __LINE__)

            CALL h5tinsert_f(offset_h5t, "LENGTH", &
                 H5OFFSETOF(C_LOC(foo), C_LOC(foo%length)), int64_h5t, hdferr)
            IF (hdferr /= 0) CALL errr(__FILE__, __LINE__)
        END BLOCK

        ! Create MPI type for offset,count
        BLOCK
            ! Local variables
            TYPE(offset_t), TARGET :: foo
            INTEGER(int32) :: blocklen(2)
            TYPE(MPI_Datatype) :: types(2)
            INTEGER(mpi_address_kind) :: base, disp(2)

            blocklen = 1

            CALL MPI_Get_address(foo%offset, disp(1))
            CALL MPI_Get_address(foo%length, disp(2))

            base = disp(1)
            DO i = 1, SIZE(disp)
                disp(i) = disp(i) - base
            END DO

            types(1) = MPI_INTEGER8
            types(2) = MPI_INTEGER8

            CALL MPI_Type_create_struct(2, blocklen, disp, types, &
                offset_mpit)
            CALL MPI_Type_commit(offset_mpit)
        END BLOCK

        ! Create MPI type for reading attributes
        BLOCK
            ! Local variables
            TYPE(rattr_t), TARGET :: foo_r
            TYPE(iattr_t), TARGET :: foo_i
            INTEGER(int32) :: blocklen(2)
            TYPE(MPI_Datatype) :: types(2)
            INTEGER(mpi_address_kind) :: base, disp(2)

            ! Number of elements in each entry - common for both
            blocklen(1) = nchar_name
            blocklen(2) = 1

            ! Real
            CALL MPI_Get_address(foo_r%name, disp(1))
            CALL MPI_Get_address(foo_r%value, disp(2))

            base = disp(1)
            DO i = 1, SIZE(disp)
                disp(i) = disp(i) - base
            END DO

            types(1) = MPI_CHARACTER
            types(2) = mglet_mpi_real

            CALL MPI_Type_create_struct(2, blocklen, disp, types, &
                rattr_mpit)
            CALL MPI_Type_commit(rattr_mpit)

            ! Integer
            CALL MPI_Get_address(foo_i%name, disp(1))
            CALL MPI_Get_address(foo_i%value, disp(2))

            base = disp(1)
            DO i = 1, SIZE(disp)
                disp(i) = disp(i) - base
            END DO

            types(1) = MPI_CHARACTER
            types(2) = mglet_mpi_int

            CALL MPI_Type_create_struct(2, blocklen, disp, types, &
                iattr_mpit)
            CALL MPI_Type_commit(iattr_mpit)
        END BLOCK

        is_init = .TRUE.
    END SUBROUTINE init_fieldio


    SUBROUTINE finish_fieldio()
        ! Local variables
        INTEGER(int32) :: hdferr

        is_init = .FALSE.

        CALL h5tclose_f(offset_h5t, hdferr)
        IF (hdferr < 0) CALL errr(__FILE__, __LINE__)

        CALL MPI_Type_free(offset_mpit)
        CALL MPI_Type_free(rattr_mpit)
        CALL MPI_Type_free(iattr_mpit)

        IF (ALLOCATED(iogr_counts)) DEALLOCATE(iogr_counts)
        IF (ALLOCATED(iogr_displs)) DEALLOCATE(iogr_displs)

        IF (ALLOCATED(io_counts)) DEALLOCATE(io_counts)
        IF (ALLOCATED(io_displs)) DEALLOCATE(io_displs)
    END SUBROUTINE finish_fieldio


    SUBROUTINE fieldio_write(parent_id, field)
        ! Subroutine arguments
        INTEGER(HID_T), INTENT(in) :: parent_id
        CLASS(basefield_t), INTENT(inout) :: field

        ! Local variables
        INTEGER(HID_T) :: group_id
        INTEGER(intk) :: i, ilevel

        ! Sanity check
        IF (.NOT. is_init) CALL errr(__FILE__, __LINE__)
        CALL start_timer(100)

        CALL hdf5common_group_open(field%name, parent_id, group_id)

        DO ilevel = minlevel, maxlevel
            ! If no grids are to be written at this level, cycle
            IF (.NOT. field%active_level(ilevel)) THEN
                CYCLE
            END IF

            ! 3D fields are "connected" to ensure that all ghost cell locations
            ! are updated before writing data to file to ease postprocessing
            IF (field%ndim == 3) THEN
                CALL connect(ilevel, 2, s1=field, corners=.TRUE.)
            END IF
        END DO

        ! Gather the field on an IO process and write it to disk
        CALL write_data(group_id, field)

        CALL hdf5common_attr_write("ISTAG", field%istag, group_id)
        CALL hdf5common_attr_write("JSTAG", field%jstag, group_id)
        CALL hdf5common_attr_write("KSTAG", field%kstag, group_id)

        CALL hdf5common_attr_write("FIELDNAME", field%description, group_id)
        CALL hdf5common_attr_write_arr("UNITS", field%units, group_id)

        DO i = 1, field%n_iattr
            CALL hdf5common_attr_write(field%iattr_key(i), field%iattr_val(i), &
                group_id)
        END DO

        DO i = 1, field%n_rattr
            CALL hdf5common_attr_write(field%rattr_key(i), field%rattr_val(i), &
                group_id)
        END DO
        CALL hdf5common_group_close(group_id)

        CALL stop_timer(100)
    END SUBROUTINE fieldio_write


    SUBROUTINE fieldio_read(parent_id, field, required)
        ! Subroutine arguments
        INTEGER(hid_t), INTENT(in) :: parent_id
        CLASS(basefield_t), INTENT(inout) :: field
        LOGICAL, INTENT(in), OPTIONAL :: required

        ! Local variables
        INTEGER(hid_t) :: group_id
        INTEGER(int32) :: ierr
        LOGICAL :: link_exists

        ! Check if field is present
        IF (PRESENT(required)) THEN
            IF (ioproc) THEN
                CALL h5lexists_f(parent_id, field%name, link_exists, ierr)
                IF (ierr /= 0) CALL errr(__FILE__, __LINE__)
            END IF
            CALL MPI_Bcast(link_exists, 1, MPI_LOGICAL, 0, MPI_COMM_WORLD)

            IF ((.NOT. required) .AND. (.NOT. link_exists)) THEN
                RETURN
            END IF
        END IF

        CALL start_timer(100)
        CALL hdf5common_group_open(field%name, parent_id, group_id)

        ! Read actual 3D field data
        CALL read_data(group_id, field)

        ! Read attributes (TSAMP etc.)
        CALL read_attrs(group_id, field)

        CALL hdf5common_group_close(group_id)
        CALL stop_timer(100)
    END SUBROUTINE fieldio_read


    SUBROUTINE read_attrs(parent_id, field)
        ! Subroutine arguments
        INTEGER(hid_t), INTENT(in) :: parent_id
        CLASS(basefield_t), INTENT(inout) :: field

        ! Local variables
        INTEGER(hid_t) :: attr_id, type_id
        INTEGER(hsize_t) :: nchar, i
        INTEGER(int32) :: hdferr
        TYPE(h5o_info_t), TARGET :: grp_info
        CHARACTER(len=nchar_name) :: attrname
        INTEGER(intk), TARGET :: intattr
        REAL(realk), TARGET :: realattr

        TYPE(rattr_t), ALLOCATABLE :: iattrs(:)
        TYPE(rattr_t), ALLOCATABLE :: rattrs(:)
        INTEGER(intk) :: n_iattrs, n_rattrs, nattr_max

        ! NB: Default integer, not an error!
        ! Ref. h5tget_class_f API documentation
        INTEGER :: type_class

        ! There is a fixed maximum number of attributes a field can have
        nattr_max = SIZE(field%iattr_val)
        ALLOCATE(iattrs(nattr_max))
        ALLOCATE(rattrs(nattr_max))
        n_iattrs = 0
        n_rattrs = 0

        ! IO processes read attributes form file
        IF (ioproc) THEN
            CALL h5oget_info_f(parent_id, grp_info, hdferr)
            IF (hdferr /= 0) CALL errr(__FILE__, __LINE__)

            DO i = 0, grp_info%num_attrs - 1
                ! https://forum.hdfgroup.org/t/error-accessing-attributes-fortran-h5aopen-by-idx/9094
                CALL H5Aopen_by_idx_f(parent_id, ".", H5_INDEX_CRT_ORDER_F, &
                    H5_ITER_INC_F, i, attr_id, hdferr)
                IF (hdferr /= 0) CALL errr(__FILE__, __LINE__)

                nchar = nchar_name
                CALL h5aget_name_f(attr_id, nchar, attrname, hdferr)
                ! Returns the length of the attribute name if successful
                IF (hdferr < 0) CALL errr(__FILE__, __LINE__)
                IF (hdferr > nchar_name) CALL errr(__FILE__, __LINE__)

                SELECT CASE(TRIM(attrname))
                CASE ("MINLEVEL", "MAXLEVEL", "ISTAG", "JSTAG", "KSTAG", &
                        "FIELDNAME", "UNITS")
                    CONTINUE
                CASE DEFAULT
                    CALL h5aget_type_f(attr_id, type_id, hdferr)
                    IF (hdferr /= 0) CALL errr(__FILE__, __LINE__)

                    CALL h5tget_class_f(type_id, type_class, hdferr)
                    IF (hdferr < 0) CALL errr(__FILE__, __LINE__)

                    IF (type_class == H5T_INTEGER_F) THEN
                        n_iattrs = n_iattrs + 1
                        IF (n_iattrs > nattr_max) CALL errr(__FILE__, __LINE__)
                        CALL read_attr(intattr, attr_id)

                        iattrs(n_iattrs)%name = &
                            TRANSFER(attrname, iattrs(n_iattrs)%name)
                        iattrs(n_iattrs)%value = intattr
                    ELSE IF (type_class == H5T_FLOAT_F) THEN
                        n_rattrs = n_rattrs + 1
                        IF (n_rattrs > nattr_max) CALL errr(__FILE__, __LINE__)
                        CALL read_attr(realattr, attr_id)

                        rattrs(n_rattrs)%name = &
                            TRANSFER(attrname, rattrs(n_rattrs)%name)
                        rattrs(n_rattrs)%value = realattr
                    END IF

                    CALL h5tclose_f(type_id, hdferr)
                    IF (hdferr /= 0) CALL errr(__FILE__, __LINE__)
                END SELECT

                CALL h5aclose_f(attr_id, hdferr)
                IF (hdferr /= 0) CALL errr(__FILE__, __LINE__)
            END DO
        END IF

        ! Distribute attributes
        CALL MPI_Bcast(n_iattrs, 1, mglet_mpi_int, 0, iogrcomm)
        CALL MPI_Bcast(iattrs, n_iattrs, iattr_mpit, 0, iogrcomm)

        CALL MPI_Bcast(n_rattrs, 1, mglet_mpi_int, 0, iogrcomm)
        CALL MPI_Bcast(rattrs, n_rattrs, rattr_mpit, 0, iogrcomm)

        ! Set attributes (phew!)
        DO i = 1, n_iattrs
            attrname = TRANSFER(iattrs(i)%name, attrname)
            CALL field%set_attr(INT(iattrs(i)%value, intk), attrname)
        END DO
        DO i = 1, n_rattrs
            attrname = TRANSFER(rattrs(i)%name, attrname)
            CALL field%set_attr(REAL(rattrs(i)%value, realk), attrname)
        END DO

        DEALLOCATE(iattrs)
        DEALLOCATE(rattrs)
    END SUBROUTINE read_attrs


    SUBROUTINE read_attr(val, attr_id)
        ! Subroutine arguments
        CLASS(*), INTENT(out), TARGET :: val
        INTEGER(hid_t), INTENT(in) :: attr_id

        ! Local variables
        INTEGER(hid_t) :: aspace_id, hdf5_dtype
        INTEGER(int32) :: hdferr
        INTEGER(hsize_t) :: npoints
        TYPE(c_ptr) :: ptr

        SELECT TYPE (val)
        TYPE IS (REAL(realk))
            hdf5_dtype = mglet_hdf5_real
            ptr = C_LOC(val)
        TYPE IS (INTEGER(intk))
            hdf5_dtype = mglet_hdf5_int
            ptr = C_LOC(val)
        CLASS DEFAULT
            CALL errr(__FILE__, __LINE__)
        END SELECT

        CALL h5aget_space_f(attr_id, aspace_id, hdferr)
        IF (hdferr < 0) CALL errr(__FILE__, __LINE__)

        CALL h5sget_simple_extent_npoints_f(aspace_id, npoints, hdferr)
        IF (hdferr < 0) CALL errr(__FILE__, __LINE__)
        IF (npoints /= 1) THEN
            WRITE(*,*) "Error opening attribute"
            WRITE(*,*) "Number of elements in file: ", npoints
            CALL errr(__FILE__, __LINE__)
        END IF

        CALL h5sclose_f(aspace_id, hdferr)
        IF (hdferr < 0) CALL errr(__FILE__, __LINE__)

        CALL h5aread_f(attr_id, hdf5_dtype, ptr, hdferr)
        IF (hdferr < 0) CALL errr(__FILE__, __LINE__)
    END SUBROUTINE read_attr


    SUBROUTINE write_data(parent_id, field)
        ! Subroutine arguments
        INTEGER(HID_T), INTENT(in) :: parent_id
        CLASS(basefield_t), INTENT(inout) :: field

        ! Local variables
        INTEGER(int64) :: my_dim, ioproc_dim
        INTEGER(intk), ALLOCATABLE :: iogridinfo(:, :)

        ! Determine the number of elements the IO process will need to
        ! receive
        my_dim = INT(field%idim, kind=int64)
        ioproc_dim = 0
        CALL MPI_Reduce(my_dim, ioproc_dim, 1, MPI_INTEGER8, MPI_SUM, 0, &
            iogrcomm)

        ! In case bigbuf is not big enough - if it is already big enough
        ! nothing happens. Only IO processes has ioproc_dim /= 0 so OK
        ! to call this collectively
        SELECT TYPE(field)
        TYPE IS (field_t)
            CALL increase_bigbuf(ioproc_dim)
        TYPE IS (intfield_t)
            CALL increase_ifkbuf(ioproc_dim)
        END SELECT

        ! Gather the data in the IO processes' bigbuf
        CALL create_iogridinfo(iogridinfo, field)

        SELECT TYPE(field)
        TYPE IS (field_t)
            CALL gather_grids(bigbuf, iogridinfo, field)
            CALL write_grids(parent_id, bigbuf, iogridinfo)
        TYPE IS (intfield_t)
            CALL gather_grids(ifkbuf, iogridinfo, field)
            CALL write_grids(parent_id, ifkbuf, iogridinfo)
        END SELECT
        CALL write_levels_vds(parent_id, field%hdf5_dtype, iogridinfo)

        DEALLOCATE(iogridinfo)
    END SUBROUTINE write_data


    ! Gather grids on an IO process before writing
    !
    ! This routine gather all data for a specific field on an IO process.
    ! The iogridinfo is allocated and filled with the nessecary metadata
    ! for constructing the HDF5 data structures. All processes
    ! need to participate in this operation.
    SUBROUTINE gather_grids(buffer, iogridinfo, field)
        ! Subroutine arguments
        CLASS(*), INTENT(out) :: buffer(:)
        INTEGER(intk), ALLOCATABLE, INTENT(in) :: iogridinfo(:, :)
        CLASS(basefield_t), INTENT(in), TARGET :: field

        ! Local variables
        INTEGER(int64) :: bufptr
        INTEGER(intk) :: i, igrid, iproc, nelems, ptr, nsend, nrecv
        INTEGER(intk) :: kk, jj, ii
        TYPE(MPI_Request), ALLOCATABLE :: sendreq(:), recvreq(:)
        CLASS(*), POINTER :: transposed(:)

        CALL start_timer(101)

        ! IO processes call MPI_Irecv first
        IF (ioproc) THEN
            ALLOCATE(recvreq(niogrgrids))

            bufptr = 1
            nrecv = 0
            DO i = 1, niogrgrids
                igrid = iogridinfo(1, i)
                iproc = idprocofgrd(igrid)
                nelems = iogridinfo(2, i)

                ! Grids with no data are not communicated
                IF (nelems == 0) CYCLE

                nrecv = nrecv + 1
                ! The SELECT TYPE here is really not neccesary, but exist as
                ! a workaround for an Intel Fortran compiler bug:
                ! https://community.intel.com/t5/Intel-Fortran-Compiler/Calling-a-function-with-TYPE-with-a-CLASS-argument/m-p/1497024
                ! When this is resolved, all that is needed is to pass
                ! the correct MPI type along
                SELECT TYPE (buffer)
                TYPE IS (REAL(realk))
                    CALL MPI_Irecv(buffer(bufptr:bufptr+nelems-1), nelems, &
                        field%mpi_dtype, iproc, igrid, iogrcomm, &
                        recvreq(nrecv))
                TYPE IS (INTEGER(ifk))
                    CALL MPI_Irecv(buffer(bufptr:bufptr+nelems-1), nelems, &
                        field%mpi_dtype, iproc, igrid, iogrcomm, &
                        recvreq(nrecv))
                END SELECT
                bufptr = bufptr + nelems
            END DO
        END IF

        ! All ranks transpose data for all grids, initiate communication,
        ! wait for communication to finish and then de-allocate the transposed
        ! data.
        !
        ! If someone thinks the (ab)use of the 'transpose' here is ugly - they
        ! are perfectly right. Feel free to come up with a more elegant
        ! solution!
        IF (field%ndim == 3) THEN
            SELECT TYPE (field)
            TYPE IS (field_t)
                ALLOCATE(transposed, mold=field%arr)
            TYPE IS (intfield_t)
                ALLOCATE(transposed, mold=field%arr)
            END SELECT
            DO i = 1, nmygrids
                igrid = mygrids(i)
                IF (.NOT. field%active_level(level(igrid))) CYCLE

                CALL get_mgdims(kk, jj, ii, igrid)
                CALL field%get_ip(ptr, igrid)
                SELECT TYPE (field)
                TYPE IS (field_t)
                    CALL fieldio_transpose(transposed(ptr:ptr+kk*jj*ii-1), &
                        field%arr(ptr:ptr+kk*jj*ii-1), kk, jj, ii)
                TYPE IS (intfield_t)
                    CALL fieldio_transpose(transposed(ptr:ptr+kk*jj*ii-1), &
                        field%arr(ptr:ptr+kk*jj*ii-1), kk, jj, ii)
                END SELECT
            END DO
        ELSE
            SELECT TYPE (field)
            TYPE IS (field_t)
                transposed => field%arr
            TYPE IS (intfield_t)
                transposed => field%arr
            END SELECT
        END IF

        ALLOCATE(sendreq(nmygrids))
        nsend = 0
        DO i = 1, nmygrids
            igrid = mygrids(i)
            IF (.NOT. field%active_level(level(igrid))) CYCLE
            nsend = nsend + 1

            CALL field%get_ip(ptr, igrid)
            CALL field%get_len(nelems, igrid)

            ! The SELECT TYPE here is really not neccesary, but exist as
            ! a workaround for an Intel Fortran compiler bug:
            ! https://community.intel.com/t5/Intel-Fortran-Compiler/Calling-a-function-with-TYPE-with-a-CLASS-argument/m-p/1497024
            ! When this is resolved, all that is needed is to pass
            ! the correct MPI type along
            SELECT TYPE (transposed)
            TYPE IS (REAL(realk))
                CALL MPI_Isend(transposed(ptr:ptr+nelems-1), nelems, &
                    field%mpi_dtype, 0, igrid, iogrcomm, sendreq(nsend))
            TYPE IS (INTEGER(ifk))
                CALL MPI_Isend(transposed(ptr:ptr+nelems-1), nelems, &
                    field%mpi_dtype, 0, igrid, iogrcomm, sendreq(nsend))
            END SELECT
        END DO
        CALL MPI_Waitall(nsend, sendreq, MPI_STATUSES_IGNORE)
        DEALLOCATE(sendreq)

        IF (field%ndim == 3) THEN
            DEALLOCATE(transposed)
        ELSE
            NULLIFY(transposed)
        END IF

        IF (ioProc) THEN
            CALL MPI_Waitall(nrecv, recvreq, MPI_STATUSES_IGNORE)
            DEALLOCATE(recvreq)
        END IF

        CALL stop_timer(101)
    END SUBROUTINE gather_grids


    ! Write grids from buffer to disk
    SUBROUTINE write_grids(parent_id, buffer, iogridinfo)
        ! Subroutine arguments
        INTEGER(HID_T), INTENT(in) :: parent_id
        CLASS(*), INTENT(in), TARGET :: buffer(:)
        INTEGER(intk), ALLOCATABLE, INTENT(in) :: iogridinfo(:, :)

        ! Local variables
        INTEGER(intk) :: i
        INTEGER(intk), ALLOCATABLE :: iogridinfo_all(:, :)
        INTEGER(int32) :: ierr
        INTEGER(hid_t) :: dset_id, filespace, memspace, plist_id, memtype
        INTEGER(hsize_t) :: shape1(1)
        INTEGER(hssize_t) :: npoints
        TYPE(offset_t), ALLOCATABLE, TARGET :: offset(:)
        TYPE(c_ptr) :: cptr

        IF (.NOT. ioproc) RETURN
        CALL start_timer(102)

        SELECT TYPE (buffer)
        TYPE IS (REAL(realk))
            cptr = C_LOC(buffer)
            memtype = mglet_hdf5_real
        TYPE IS (INTEGER(ifk))
            cptr = C_LOC(buffer)
            memtype = mglet_hdf5_ifk
        CLASS DEFAULT
            CALL errr(__FILE__, __LINE__)
        END SELECT

        CALL create_iogridinfo_all(iogridinfo_all, iogridinfo)

        ! Compute global offset and length information
        ALLOCATE(offset(ngrid))
        offset(1)%offset = 1
        DO i = 1, ngrid
            offset(i)%length = iogridinfo_all(2, i)
        END DO
        DO i = 2, ngrid
            offset(i)%offset = offset(i-1)%offset + offset(i-1)%length
        END DO

        ! Create/open dataset
        shape1(1) = offset(ngrid)%offset + offset(ngrid)%length - 1
        CALL hdf5common_dataset_create("DATA", shape1, memtype, &
            parent_id, dset_id, filespace)

        ! This makes the filespace hyperslab selection, returns the number
        ! of points selected
        CALL select_hyperslab_data(filespace, npoints, iogridinfo, offset)

        ! Create memspace and select entire memspace
        shape1(1) = npoints
        CALL h5screate_simple_f(1, shape1, memspace, ierr)
        IF (ierr /= 0) CALL errr(__FILE__, __LINE__)
        CALL h5sselect_all_f(memspace, ierr)
        IF (ierr /= 0) CALL errr(__FILE__, __LINE__)

        ! Property list collectively for collective dataset write
        CALL h5pcreate_f(H5P_DATASET_XFER_F, plist_id, ierr)
        IF (ierr /= 0) CALL errr(__FILE__, __LINE__)

        CALL h5pset_dxpl_mpio_f(plist_id, hdf5_io_mode, ierr)
        IF (ierr /= 0) CALL errr(__FILE__, __LINE__)

        CALL h5dwrite_f(dset_id, memtype, cptr, ierr, &
            file_space_id=filespace, mem_space_id=memspace, &
            xfer_prp=plist_id)
        IF (ierr /= 0) CALL errr(__FILE__, __LINE__)

        ! Close handles
        CALL h5pclose_f(plist_id, ierr)
        IF (ierr /= 0) CALL errr(__FILE__, __LINE__)

        CALL h5sclose_f(memspace, ierr)
        IF (ierr /= 0) CALL errr(__FILE__, __LINE__)

        CALL hdf5common_dataset_close(dset_id, filespace)

        ! Write offset, length table
        cptr = C_LOC(offset)
        shape1(1) = ngrid
        CALL stencilio_write_master_cptr(parent_id, "OFFSET", cptr, &
            shape1, offset_h5t)

        CALL stop_timer(102)
    END SUBROUTINE write_grids


    ! Write virtraul datasets "LEVEL" that point into the "DATA" array
    ! TODO: legacy - remove this feature
    SUBROUTINE write_levels_vds(parent_id, hdf5_dtype, iogridinfo)
        ! Subroutine arguments
        INTEGER(HID_T), INTENT(in) :: parent_id
        INTEGER(HID_T), INTENT(in) :: hdf5_dtype
        INTEGER(intk), ALLOCATABLE, INTENT(in) :: iogridinfo(:, :)

        ! Local variables
        LOGICAL :: link_exists
        INTEGER(int64) :: bufptr
        INTEGER(intk) :: i, ilevel, nelems
        INTEGER(intk) :: this_minlevel, this_maxlevel
        INTEGER(intk), ALLOCATABLE :: iogridinfo_all(:, :)
        INTEGER(intk), ALLOCATABLE :: igridlvl(:), nofthislevel(:), &
            maxofthislevel(:)
        INTEGER(int32) :: ierr
        INTEGER(hid_t) :: dset_id, dcpl, vspace, src_space
        INTEGER(hsize_t) :: shape1(1), shape2(2), offset1(1), dims1(1)
        TYPE(offset_t), ALLOCATABLE, TARGET :: offset(:)
        CHARACTER(LEN=32) :: dsetname
        CHARACTER(LEN=32) :: grpname
        INTEGER(size_t) :: name_size, buf_size

        IF (.NOT. ioproc) RETURN

        ! If array "IGRIDLVL" exists already we can skip this function
        ! because it would not write any new data only virtual datasets
        CALL h5lexists_f(parent_id, "IGRIDLVL", link_exists, ierr)
        IF (ierr < 0) CALL errr(__FILE__, __LINE__)
        IF (link_exists) RETURN

        CALL create_iogridinfo_all(iogridinfo_all, iogridinfo)

        ! Compute global offset and length information
        ALLOCATE(offset(ngrid))
        offset(1)%offset = 1
        DO i = 1, ngrid
            offset(i)%length = iogridinfo_all(2, i)
        END DO
        DO i = 2, ngrid
            offset(i)%offset = offset(i-1)%offset + offset(i-1)%length
        END DO

        ! Create igridlvl
        this_minlevel = maxlevel
        this_maxlevel = minlevel

        ALLOCATE(igridlvl(ngrid))
        ALLOCATE(nofthislevel(minlevel:maxlevel))
        ALLOCATE(maxofthislevel(minlevel:maxlevel))

        nofthislevel = 0
        maxofthislevel = 0
        DO i = 1, ngrid
            ilevel = level(i)
            nofthislevel(ilevel) = nofthislevel(ilevel) + 1
            igridlvl(i) = nofthislevel(ilevel)

            nelems = iogridinfo_all(2, i)
            maxofthislevel(ilevel) = MAX(maxofthislevel(ilevel), nelems)

            IF (nelems > 0) THEN
                this_minlevel = MIN(this_minlevel, ilevel)
                this_maxlevel = MAX(this_maxlevel, ilevel)
            END IF
        END DO

        ! Construct path to DATA dataset
        buf_size = LEN(grpname)
        CALL h5iget_name_f(parent_id, grpname, buf_size, name_size, ierr)
        IF (ierr /= 0) CALL errr(__FILE__, __LINE__)

        ! Create datasets
        bufptr = 0
        DO ilevel = this_minlevel, this_maxlevel
            ! Dataspace for the virtual dataset ("LEVEL1")
            shape2(1) = maxofthislevel(ilevel)
            shape2(2) = nofthislevel(ilevel)
            CALL h5screate_simple_f(2, shape2, vspace, ierr)
            IF (ierr /= 0) CALL errr(__FILE__, __LINE__)

            ! Dataspace for source data ("DATA")
            shape1(1) = offset(ngrid)%offset + offset(ngrid)%length - 1
            CALL h5screate_simple_f(1, shape1, src_space, ierr)
            IF (ierr /= 0) CALL errr(__FILE__, __LINE__)

            ! Only when all grids at the level has the same shape
            offset1(1) = bufptr
            dims1(1) = nofthislevel(ilevel)*maxofthislevel(ilevel)
            bufptr = bufptr + dims1(1)
            CALL h5sselect_hyperslab_f(src_space, H5S_SELECT_SET_F, offset1, &
                dims1, ierr)
            IF (ierr /= 0) CALL errr(__FILE__, __LINE__)

            ! Dataset creation property list and set virtual properties
            CALL h5pcreate_f(H5P_DATASET_CREATE_F, dcpl, ierr)
            IF (ierr /= 0) CALL errr(__FILE__, __LINE__)
            CALL h5pset_virtual_f(dcpl, vspace, ".", TRIM(grpname)//"/DATA", &
                src_space, ierr)
            IF (ierr /= 0) CALL errr(__FILE__, __LINE__)

            ! Create dataset
            WRITE(dsetname, '("LEVEL", i0)') ilevel
            CALL h5dcreate_f(parent_id, dsetname, hdf5_dtype, vspace, &
                dset_id, ierr, dcpl_id=dcpl)
            IF (ierr /= 0) CALL errr(__FILE__, __LINE__)

            ! Close handles
            CALL h5dclose_f(dset_id, ierr)
            CALL h5pclose_f(dcpl, ierr)
            CALL h5sclose_f(src_space, ierr)
            CALL h5sclose_f(vspace, ierr)
        END DO

        ! Write IGRIDLVL
        CALL stencilio_write_master(parent_id, "IGRIDLVL", igridlvl)

        ! Write MINLEVEL and MAXLEVEL
        CALL hdf5common_attr_write("MINLEVEL", this_minlevel, parent_id)
        CALL hdf5common_attr_write("MAXLEVEL", this_maxlevel, parent_id)

        DEALLOCATE(maxofthislevel)
        DEALLOCATE(nofthislevel)
        DEALLOCATE(igridlvl)
    END SUBROUTINE write_levels_vds


    SUBROUTINE read_data(parent_id, field)
        ! Subroutine arguments
        INTEGER(HID_T), INTENT(in) :: parent_id
        CLASS(basefield_t), INTENT(inout) :: field

        ! Local variables
        LOGICAL :: link_exists
        INTEGER(int64) :: my_dim, ioproc_dim
        INTEGER(intk), ALLOCATABLE :: iogridinfo(:, :)

        ! Determine the number of elements the IO process will need to
        ! receive
        my_dim = INT(field%idim, kind=int64)
        ioproc_dim = 0
        CALL MPI_Reduce(my_dim, ioproc_dim, 1, MPI_INTEGER8, MPI_SUM, 0, &
            iogrcomm)

        ! In case bigbuf is not big enough - if it is already big enough
        ! nothing happens. Only IO processes has ioproc_dim /= 0 so OK
        ! to call this collectively
        SELECT TYPE(field)
        TYPE IS (field_t)
            CALL increase_bigbuf(ioproc_dim)
        TYPE IS (intfield_t)
            CALL increase_ifkbuf(ioproc_dim)
        END SELECT

        CALL create_iogridinfo(iogridinfo, field)

        ! Variant A: read "DATA", which is the 1-D conracted way of writing
        ! out fields - how MGLET does it presently
        !
        ! Variant B: read "LEVEL0" etc. - legacy file structure - will be
        ! deprecated in the future
        CALL hdf5common_dataset_exists("DATA", parent_id, link_exists)

        SELECT TYPE(field)
        TYPE IS (field_t)
            IF (link_exists) THEN
                CALL read_grids_data(parent_id, bigbuf, iogridinfo)
            ELSE
                CALL read_grids_levels(parent_id, bigbuf, iogridinfo)
            END IF
            CALL scatter_grids(bigbuf, iogridinfo, field)
        TYPE IS (intfield_t)
            IF (link_exists) THEN
                CALL read_grids_data(parent_id, ifkbuf, iogridinfo)
            ELSE
                CALL read_grids_levels(parent_id, ifkbuf, iogridinfo)
            END IF
            CALL scatter_grids(ifkbuf, iogridinfo, field)
        END SELECT

        DEALLOCATE(iogridinfo)
    END SUBROUTINE read_data


    ! Read data level by level (LEVEL1, LEVEL2 etc.)
    ! TODO: legacy - remove this feature
    SUBROUTINE read_grids_levels(parent_id, buffer, iogridinfo)
        ! Subroutine arguments
        INTEGER(HID_T), INTENT(in) :: parent_id
        CLASS(*), INTENT(inout), TARGET :: buffer(:)
        INTEGER(intk), ALLOCATABLE, INTENT(in) :: iogridinfo(:, :)

        ! Local variables
        LOGICAL :: finishblock, link_exists
        INTEGER(intk) :: i, igrid, ilevel
        INTEGER(intk) :: this_minlevel, this_maxlevel
        INTEGER(intk), ALLOCATABLE :: igridlvl(:), nofthislevel(:)
        INTEGER(int32) :: ierr
        INTEGER(hid_t) :: dset_id, filespace, memspace, plist_id, memtype
        INTEGER(hsize_t) :: count1(1), shape2(2), count2(2), offset2(2)
        INTEGER(hssize_t) :: npoints, count_m, bufptr
        TYPE(c_ptr) :: cptr
        CHARACTER(LEN=32) :: dsetname

        IF (.NOT. ioproc) RETURN
        CALL start_timer(103)

        ! If 'IGRIDLVL' is present - read it - if not - make it
        ALLOCATE(igridlvl(ngrid))
        igridlvl = 0

        ! CALL h5lexists_f(parent_id, "IGRIDLVL", link_exists, ierr)
        ! IF (ierr < 0) CALL errr(__FILE__, __LINE__)

        ! "Legacy" MGLET writes out an igridlevel that can have more
        ! entries than grids during CHECKBLOCK - this causes a problem here
        ! when it is allocated to igridlevel length and there is not
        ! enough space in memory to read in the entire thig. Therefore
        ! it is maybe better to always create it.
        link_exists = .FALSE.
        IF (link_exists) THEN
            CALL stencilio_read_master(parent_id, "IGRIDLVL", igridlvl)
        ELSE
            ALLOCATE(nofthislevel(minlevel:maxlevel))
            nofthislevel = 0
            DO i = 1, ngrid
                ilevel = level(i)
                nofthislevel(ilevel) = nofthislevel(ilevel) + 1
                igridlvl(i) = nofthislevel(ilevel)
            END DO
            DEALLOCATE(nofthislevel)
        END IF
        CALL MPI_Bcast(igridlvl, ngrid, mglet_mpi_int, 0, iocomm)

        ! Sanity check on dataset levels
        CALL hdf5common_attr_read("MINLEVEL", this_minlevel, parent_id)
        CALL hdf5common_attr_read("MAXLEVEL", this_maxlevel, parent_id)

        ! Level by level read the datasets
        bufptr = 1
        DO ilevel = this_minlevel, this_maxlevel
            WRITE(dsetname, '("LEVEL", i0)') ilevel
            CALL hdf5common_dataset_open(dsetname, shape2, parent_id, dset_id, &
                filespace)
            CALL h5sselect_none_f(filespace, ierr)
            IF (ierr /= 0) CALL errr(__FILE__, __LINE__)

            ! set up hyperslab dimension arrays in file
            offset2 = 0
            count2 = 0
            count_m = 0
            finishblock = .FALSE.
            DO i = 1, niogrgrids
                igrid = iogridinfo(1, i)
                IF (level(igrid) /= ilevel) CYCLE

                ! start new block?
                IF (count2(2) == 0) THEN
                    offset2(1) = 0
                    offset2(2) = igridlvl(igrid) - 1
                    count2(1) = iogridinfo(2, i)
                    count2(2) = 0
                END IF

                count2(2) = count2(2) + 1

                ! Finish current block?
                IF (i < niogrgrids) THEN
                    ! If it's non-contigous in filespace - stop
                    IF (igridlvl(iogridinfo(1, i+1)) /= &
                            igridlvl(igrid) + 1) THEN
                        finishblock = .TRUE.
                    END IF

                    ! If length of next grid is different - stop
                    IF (iogridinfo(2, i+1) /= iogridinfo(2, i)) THEN
                        finishblock = .TRUE.
                    END IF

                    ! Next grid is at a different level - stop
                    IF (level(iogridinfo(1, i+1)) /= ilevel) THEN
                        finishblock = .TRUE.
                    END IF
                ELSE
                    ! we reach to the end grid of this level
                    finishblock = .TRUE.
                END IF

                IF (finishblock) THEN
                    CALL h5sselect_hyperslab_f(filespace, H5S_SELECT_OR_F, &
                        offset2, count2, ierr)
                    IF (ierr /= 0) CALL errr(__FILE__, __LINE__)

                    ! reset for next iteration
                    count_m = count_m + count2(1)*count2(2)
                    offset2 = 0
                    count2 = 0
                    finishblock = .FALSE.
                END IF
            END DO

            ! Sanity check
            CALL h5sget_select_npoints_f(filespace, npoints, ierr)
            IF (ierr /= 0) CALL errr(__FILE__, __LINE__)

            IF (npoints /= count_m) THEN
                WRITE(*,*) "npoints, count_m", npoints, count_m
                CALL errr(__FILE__, __LINE__)
            END IF

            ! Create memspace and select entire memspace
            count1(1) = count_m
            CALL h5screate_simple_f(1, count1, memspace, ierr)
            IF (ierr /= 0) CALL errr(__FILE__, __LINE__)
            CALL h5sselect_all_f(memspace, ierr)
            IF (ierr /= 0) CALL errr(__FILE__, __LINE__)

            ! Property list collectively for collective dataset write
            CALL h5pcreate_f(H5P_DATASET_XFER_F, plist_id, ierr)
            IF (ierr /= 0) CALL errr(__FILE__, __LINE__)

            CALL h5pset_dxpl_mpio_f(plist_id, hdf5_io_mode, ierr)
            IF (ierr /= 0) CALL errr(__FILE__, __LINE__)

            ! Sanity check on buffer size
            IF (bufptr + count_m - 1 > SIZE(buffer)) THEN
                CALL errr(__FILE__, __LINE__)
            END IF

            SELECT TYPE (buffer)
            TYPE IS (REAL(realk))
                cptr = C_LOC(buffer(bufptr))
                memtype = mglet_hdf5_real
            TYPE IS (INTEGER(ifk))
                cptr = C_LOC(buffer(bufptr))
                memtype = mglet_hdf5_ifk
            CLASS DEFAULT
                CALL errr(__FILE__, __LINE__)
            END SELECT
            CALL h5dread_f(dset_id, memtype, cptr, ierr, &
                 file_space_id=filespace, mem_space_id=memspace, &
                 xfer_prp=plist_id)
            IF (ierr /= 0) CALL errr(__FILE__, __LINE__)
            bufptr = bufptr + count_m

            ! Close handles
            CALL h5pclose_f(plist_id, ierr)
            IF (ierr /= 0) CALL errr(__FILE__, __LINE__)

            CALL h5sclose_f(memspace, ierr)
            IF (ierr /= 0) CALL errr(__FILE__, __LINE__)

            CALL hdf5common_dataset_close(dset_id, filespace)
        END DO

        CALL stop_timer(103)
    END SUBROUTINE read_grids_levels


    ! Read data all levels at once (DATA array)
    SUBROUTINE read_grids_data(parent_id, buffer, iogridinfo)
        ! Subroutine arguments
        INTEGER(HID_T), INTENT(in) :: parent_id
        CLASS(*), INTENT(inout), TARGET :: buffer(:)
        INTEGER(intk), ALLOCATABLE, INTENT(in) :: iogridinfo(:, :)

        ! Local variables
        INTEGER(int32) :: ierr
        INTEGER(hid_t) :: dset_id, filespace, memspace, plist_id, memtype
        INTEGER(hsize_t) :: shape1(1)
        INTEGER(hssize_t) :: npoints
        TYPE(offset_t), ALLOCATABLE, TARGET :: offset(:)
        TYPE(c_ptr) :: cptr

        IF (.NOT. ioproc) RETURN
        CALL start_timer(103)

        ! Read 'OFFSET' and broadcast result to all IO processes
        ALLOCATE(offset(ngrid))
        shape1(1) = ngrid
        cptr = C_LOC(offset)
        CALL stencilio_read_master_cptr(parent_id, "OFFSET", cptr, shape1, &
            offset_h5t)
        CALL MPI_Bcast(offset, ngrid, offset_mpit, 0, iocomm)

        ! Open dataset
        CALL hdf5common_dataset_open("DATA", shape1, parent_id, dset_id, &
            filespace)

        ! This makes the filespace hyperslab selection, returns the number
        ! of points selected
        CALL select_hyperslab_data(filespace, npoints, iogridinfo, offset)

        ! Create memspace and select entire memspace
        shape1(1) = npoints
        CALL h5screate_simple_f(1, shape1, memspace, ierr)
        IF (ierr /= 0) CALL errr(__FILE__, __LINE__)
        CALL h5sselect_all_f(memspace, ierr)
        IF (ierr /= 0) CALL errr(__FILE__, __LINE__)

        ! Property list collectively for collective dataset write
        CALL h5pcreate_f(H5P_DATASET_XFER_F, plist_id, ierr)
        IF (ierr /= 0) CALL errr(__FILE__, __LINE__)

        CALL h5pset_dxpl_mpio_f(plist_id, hdf5_io_mode, ierr)
        IF (ierr /= 0) CALL errr(__FILE__, __LINE__)

        SELECT TYPE (buffer)
        TYPE IS (REAL(realk))
            cptr = C_LOC(buffer)
            memtype = mglet_hdf5_real
        TYPE IS (INTEGER(ifk))
            cptr = C_LOC(buffer)
            memtype = mglet_hdf5_ifk
        CLASS DEFAULT
            CALL errr(__FILE__, __LINE__)
        END SELECT
        CALL h5dread_f(dset_id, memtype, cptr, ierr, &
            file_space_id=filespace, mem_space_id=memspace, &
            xfer_prp=plist_id)
        IF (ierr /= 0) CALL errr(__FILE__, __LINE__)

        ! Close handles
        CALL h5pclose_f(plist_id, ierr)
        IF (ierr /= 0) CALL errr(__FILE__, __LINE__)

        CALL h5sclose_f(memspace, ierr)
        IF (ierr /= 0) CALL errr(__FILE__, __LINE__)

        CALL hdf5common_dataset_close(dset_id, filespace)
        CALL stop_timer(103)
    END SUBROUTINE read_grids_data


    ! Scatter data from buffer to individual processes
    SUBROUTINE scatter_grids(buffer, iogridinfo, field)
        ! Subroutine arguments
        CLASS(*), INTENT(in) :: buffer(:)
        INTEGER(intk), ALLOCATABLE, INTENT(in) :: iogridinfo(:, :)
        CLASS(basefield_t), INTENT(inout), TARGET :: field

        ! Local variables
        INTEGER(int64) :: bufptr
        INTEGER(intk) :: i, igrid, iproc, nelems, ptr, nsend, nrecv
        INTEGER(intk) :: kk, jj, ii
        TYPE(MPI_Request), ALLOCATABLE :: sendreq(:), recvreq(:)

        CALL start_timer(104)

        ! All processes call Recv first
        ALLOCATE(recvreq(nmygrids))
        nrecv = 0
        DO i = 1, nmygrids
            igrid = mygrids(i)
            IF (.NOT. field%active_level(level(igrid))) CYCLE
            nrecv = nrecv + 1

            CALL field%get_ip(ptr, igrid)
            CALL field%get_len(nelems, igrid)
            SELECT TYPE (field)
            TYPE IS (field_t)
                CALL MPI_Irecv(field%arr(ptr:ptr+nelems-1), nelems, &
                    field%mpi_dtype, 0, igrid, iogrcomm, recvreq(nrecv))
            TYPE IS (intfield_t)
                CALL MPI_Irecv(field%arr(ptr:ptr+nelems-1), nelems, &
                    field%mpi_dtype, 0, igrid, iogrcomm, recvreq(nrecv))
            END SELECT
        END DO

        ! IO processes send data
        IF (ioproc) THEN
            ALLOCATE(sendreq(niogrgrids))

            bufptr = 1
            nsend = 0
            DO i = 1, niogrgrids
                igrid = iogridinfo(1, i)
                iproc = idprocofgrd(igrid)
                nelems = iogridinfo(2, i)

                ! Grids with no data are not communicated
                IF (nelems == 0) CYCLE

                nsend = nsend + 1
                SELECT TYPE (buffer)
                TYPE IS (REAL(realk))
                    CALL MPI_Isend(buffer(bufptr:bufptr+nelems-1), nelems, &
                        field%mpi_dtype, iproc, igrid, iogrcomm, &
                        sendreq(nsend))
                TYPE IS (INTEGER(ifk))
                    CALL MPI_Isend(buffer(bufptr:bufptr+nelems-1), nelems, &
                        field%mpi_dtype, iproc, igrid, iogrcomm, &
                        sendreq(nsend))
                END SELECT
                bufptr = bufptr + nelems
            END DO
        END IF

        ! Wait for communication to complete
        CALL MPI_Waitall(nrecv, recvreq, MPI_STATUSES_IGNORE)
        DEALLOCATE(recvreq)

        IF (ioproc) THEN
            CALL MPI_Waitall(nsend, sendreq, MPI_STATUSES_IGNORE)
            DEALLOCATE(sendreq)
        END IF

        ! Transpose data in-place
        ! IMPORTANT: Deliberate use of "wrong" ordering of kk, jj, ii to
        ! transpose backwards!
        IF (field%ndim == 3) THEN
            DO i = 1, nmygrids
                igrid = mygrids(i)
                IF (.NOT. field%active_level(level(igrid))) CYCLE

                CALL get_mgdims(kk, jj, ii, igrid)
                CALL field%get_ip(ptr, igrid)

                SELECT TYPE (field)
                TYPE IS (field_t)
                    CALL fieldio_transpose_inplace( &
                        field%arr(ptr:ptr+kk*jj*ii-1), ii, jj, kk)
                TYPE IS (intfield_t)
                    CALL fieldio_transpose_inplace( &
                        field%arr(ptr:ptr+kk*jj*ii-1), ii, jj, kk)
                END SELECT
            END DO
        END IF

        CALL stop_timer(104)
    END SUBROUTINE scatter_grids


    SUBROUTINE create_iogridinfo(iogridinfo, field)
        ! Subroutine arguments
        INTEGER(intk), ALLOCATABLE, INTENT(out) :: iogridinfo(:, :)
        CLASS(basefield_t), INTENT(in), TARGET :: field

        ! Local variables
        INTEGER(intk) :: i, igrid
        INTEGER(intk), ALLOCATABLE :: mygridinfo(:, :), iogridinfo_tmp(:, :)
        INTEGER(intk), ALLOCATABLE :: sortidx(:)

        ! Fill a table with igrid and length for every grid. Not just the grids
        ! to be written out, but every grid. Grids with length == 0 are nether
        ! communicated nor written out, but this serve as the basis for
        ! lookup-tables when reading the file back in, either in MGLET or
        ! in various postprocessing tools. Therefore it is useful to
        ! know which grids are absent.
        ALLOCATE(mygridinfo(info_dim, nmygrids))
        mygridinfo = 0
        DO i = 1, nmygrids
            igrid = mygrids(i)
            mygridinfo(1, i) = igrid

            IF (.NOT. field%active_level(level(igrid))) CYCLE
            CALL field%get_len(mygridinfo(2, i), igrid)
        END DO

        IF (ioproc) THEN
            ALLOCATE(iogridinfo(info_dim, niogrgrids))
        ELSE
            ALLOCATE(iogridinfo(info_dim, 0))
        END IF
        CALL MPI_Gatherv(mygridinfo, info_dim*nmygrids, mglet_mpi_int, &
            iogridinfo, iogr_counts, iogr_displs, mglet_mpi_int, &
            0, iogrcomm)
        DEALLOCATE(mygridinfo)

        ! IO processes call MPI_Irecv first
        IF (ioproc) THEN
            ! Re-arrange iogrdinfo to get it sorted by igrid
            ALLOCATE(sortidx(niogrgrids))
            CALL sortix(niogrgrids, iogridinfo(1, :), sortidx)
            ALLOCATE(iogridinfo_tmp, MOLD=iogridinfo)
            DO i = 1, niogrgrids
                iogridinfo_tmp(:, i) = iogridinfo(:, sortidx(i))
            END DO
            DEALLOCATE(sortidx)
            CALL MOVE_ALLOC(iogridinfo_tmp, iogridinfo)
        END IF
    END SUBROUTINE create_iogridinfo


    SUBROUTINE create_iogridinfo_all(iogridinfo_all, iogridinfo)
        ! Subroutine arguments
        INTEGER(intk), ALLOCATABLE, INTENT(out) :: iogridinfo_all(:, :)
        INTEGER(intk), INTENT(in) :: iogridinfo(:, :)

        ! Local variables
        INTEGER(intk) :: i
        INTEGER(intk), ALLOCATABLE :: sortidx(:)
        INTEGER(intk), ALLOCATABLE :: iogridinfo_tmp(:, :)

        ALLOCATE(iogridinfo_all(info_dim, ngrid))
        CALL MPI_Allgatherv(iogridinfo, info_dim*niogrgrids, &
            mglet_mpi_int, iogridinfo_all, io_counts, io_displs, &
            mglet_mpi_int, iocomm)

        ! iogridinfo_all needs sorting by grid id
        ALLOCATE(sortidx(ngrid))
        CALL sortix(ngrid, iogridinfo_all(1, :), sortidx)

        ALLOCATE(iogridinfo_tmp, mold=iogridinfo_all)
        DO i = 1, ngrid
            iogridinfo_tmp(:, i) = iogridinfo_all(:, sortidx(i))
        END DO

        ! Sanity check to assure that we have received information
        ! for all grids by checking that the grid-id is sequential and
        ! no grids are missing
        ! TODO: encapsulate in _MGLET_DEBUG_
        DO i = 1, ngrid-1
            IF (iogridinfo_tmp(1, i+1) /= iogridinfo_tmp(1, i) + 1) THEN
                CALL errr(__FILE__, __LINE__)
            END IF
        END DO

        DEALLOCATE(sortidx)
        CALL MOVE_ALLOC(iogridinfo_tmp, iogridinfo_all)
    END SUBROUTINE create_iogridinfo_all


    SUBROUTINE select_hyperslab_data(filespace, npoints, iogridinfo, offset)
        ! Subroutine arguments
        INTEGER(hid_t), INTENT(inout) :: filespace
        INTEGER(hssize_t), INTENT(out) :: npoints
        INTEGER(intk), INTENT(in) :: iogridinfo(:, :)
        TYPE(offset_t), INTENT(in) :: offset(:)

        ! Local variables
        LOGICAL :: finishblock
        INTEGER(intk) :: i, igrid
        INTEGER(hsize_t) :: count1(1), offset1(1)
        INTEGER(hssize_t) :: count_m
        INTEGER(int32) :: ierr

        CALL h5sselect_none_f(filespace, ierr)
        IF (ierr /= 0) CALL errr(__FILE__, __LINE__)

        ! set up hyperslab dimension arrays in file
        offset1 = 0
        count1 = 0
        count_m = 0

        finishblock = .FALSE.
        DO i = 1, niogrgrids
            igrid = iogridinfo(1, i)

            ! start new block?
            IF (count1(1) == 0) THEN
                offset1(1) = offset(igrid)%offset - 1
            END IF

            count1(1) = count1(1) + iogridinfo(2, i)

            ! Finish current block?
            IF (i < niogrgrids) THEN
                ! If it's non-contigous in filespace, stop block here
                IF (iogridinfo(1, i+1) /= igrid + 1) THEN
                    finishblock = .TRUE.
                END IF
            ELSE
                ! we reach to the end grid of this level
                finishblock = .TRUE.
            END IF

            IF (finishblock) THEN
                CALL h5sselect_hyperslab_f(filespace, H5S_SELECT_OR_F, &
                    offset1, count1, ierr)
                IF (ierr /= 0) CALL errr(__FILE__, __LINE__)

                ! reset for next iteration
                count_m = count_m + count1(1)
                offset1 = 0
                count1 = 0
                finishblock = .FALSE.
            END IF
        END DO

        ! Sanity check
        npoints = 0
        CALL h5sget_select_npoints_f(filespace, npoints, ierr)
        IF (ierr /= 0) CALL errr(__FILE__, __LINE__)

        IF (npoints /= count_m) THEN
            WRITE(*,*) "npoints, count_m", npoints, count_m
            CALL errr(__FILE__, __LINE__)
        END IF
    END SUBROUTINE select_hyperslab_data


    SUBROUTINE fieldio_transpose(outfield, field, kk, jj, ii)
        ! Subroutine arguments
        CLASS(*), INTENT(out), CONTIGUOUS :: outfield(:)
        CLASS(*), INTENT(in), CONTIGUOUS :: field(:)
        INTEGER(intk), INTENT(in) :: kk, jj, ii

        ! Local variables
        INTEGER(intk) :: offset_in, offset_out
        INTEGER(intk) :: k, j, i

        ! Check input types
        SELECT TYPE (field)
        TYPE IS (REAL(realk))
            SELECT TYPE (outfield)
            TYPE IS (REAL(realk))
                DO i = 1, ii
                    DO j = 1, jj
                        DO k = 1, kk
                            ! Regular MGLET ordering offset
                            offset_in  = k + (j-1)*kk + (i-1)*kk*jj
                            ! Offset in transposed array
                            offset_out = i + (j-1)*ii + (k-1)*ii*jj
                            ! Assignment of value
                            outfield(offset_out) = field(offset_in)
                        END DO
                    END DO
                END DO
            END SELECT
        TYPE IS (INTEGER(ifk))
            SELECT TYPE (outfield)
            TYPE IS (INTEGER(ifk))
                DO i = 1, ii
                    DO j = 1, jj
                        DO k = 1, kk
                            ! Regular MGLET ordering offset
                            offset_in  = k + (j-1)*kk + (i-1)*kk*jj
                            ! Offset in transposed array
                            offset_out = i + (j-1)*ii + (k-1)*ii*jj
                            ! Assignment of value
                            outfield(offset_out) = field(offset_in)
                        END DO
                    END DO
                END DO
            END SELECT
        CLASS DEFAULT
            CALL errr(__FILE__, __LINE__)
        END SELECT
    END SUBROUTINE fieldio_transpose


    SUBROUTINE fieldio_transpose_inplace(field, kk, jj, ii)
        ! Subroutine arguments
        CLASS(*), INTENT(inout), CONTIGUOUS :: field(:)
        INTEGER(intk), INTENT(in) :: kk, jj, ii

        ! Local variables
        INTEGER(intk) :: offset_in, offset_out
        INTEGER(intk) :: k, j, i

        CLASS(*), ALLOCATABLE :: transposed(:)

        ! Check input types
        SELECT TYPE (field)
        TYPE IS (REAL(realk))
            ALLOCATE(REAL(realk) :: transposed(kk*jj*ii))
            SELECT TYPE (transposed)
            TYPE IS (REAL(realk))
                DO i = 1, ii
                    DO j = 1, jj
                        DO k = 1, kk
                            ! Regular MGLET ordering offset
                            offset_in  = k + (j-1)*kk + (i-1)*kk*jj
                            ! Offset in transposed array
                            offset_out = i + (j-1)*ii + (k-1)*ii*jj
                            ! Assignment of value
                            transposed(offset_out) = field(offset_in)
                        END DO
                    END DO
                END DO
                field = transposed
            END SELECT
            DEALLOCATE(transposed)
        TYPE IS (INTEGER(ifk))
            ALLOCATE(INTEGER(ifk) :: transposed(kk*jj*ii))
            SELECT TYPE (transposed)
            TYPE IS (INTEGER(ifk))
                DO i = 1, ii
                    DO j = 1, jj
                        DO k = 1, kk
                            ! Regular MGLET ordering offset
                            offset_in  = k + (j-1)*kk + (i-1)*kk*jj
                            ! Offset in transposed array
                            offset_out = i + (j-1)*ii + (k-1)*ii*jj
                            ! Assignment of value
                            transposed(offset_out) = field(offset_in)
                        END DO
                    END DO
                END DO
                field = transposed
            END SELECT
            DEALLOCATE(transposed)
        CLASS DEFAULT
            CALL errr(__FILE__, __LINE__)
        END SELECT
    END SUBROUTINE fieldio_transpose_inplace

END MODULE fieldio2_mod
