MODULE stencilio_mod
    USE, INTRINSIC ::  ISO_C_BINDING, ONLY: c_ptr, c_loc

    ! Intel compilers have a problem compiling this file if using MPI_F08:
    ! https://community.intel.com/t5/Intel-Fortran-Compiler/Calling-a-function-with-TYPE-with-a-CLASS-argument/m-p/1497024
    USE MPI
    USE HDF5
    USE comms_mod
    USE err_mod, ONLY: errr
    USE hdf5common_mod
    USE precision_mod
    USE qsort_mod, ONLY: sortix

    IMPLICIT NONE
    PRIVATE

    TYPE :: int_stencils_t
        INTEGER(intk), ALLOCATABLE :: arr(:)
    CONTAINS
        PRIVATE
        FINAL :: int_stencils_destructor
    END TYPE int_stencils_t

    TYPE :: real_stencils_t
        REAL(realk), ALLOCATABLE :: arr(:)
    CONTAINS
        PRIVATE
        FINAL :: real_stencils_destructor
    END TYPE real_stencils_t

    INTEGER(intk), PARAMETER :: maxrank = 4

    PUBLIC :: stencilio_write_list, stencilio_write, stencilio_read_list, &
        stencilio_read, stencilio_read_master, stencilio_read_master_cptr, &
        stencilio_write_master, stencilio_write_master_cptr, &
        stencilio_append_master_cptr, int_stencils_t, real_stencils_t

CONTAINS

    PURE ELEMENTAL SUBROUTINE int_stencils_destructor(this)
        TYPE(int_stencils_t), INTENT(inout) :: this
        IF (ALLOCATED(this%arr)) DEALLOCATE(this%arr)
    END SUBROUTINE int_stencils_destructor


    PURE ELEMENTAL SUBROUTINE real_stencils_destructor(this)
        TYPE(real_stencils_t), INTENT(inout) :: this
        IF (ALLOCATED(this%arr)) DEALLOCATE(this%arr)
    END SUBROUTINE real_stencils_destructor


    SUBROUTINE stencilio_write(group_id, name, stencils, get_ptr, get_len, &
            indexlist, extend, ncmp, same_kind)
        USE grids_mod, ONLY: mygrids, nmygrids
        USE commbuf_mod, ONLY: bigbuf, intbuf, increase_bigbuf, increase_intbuf

        ! Subroutine arguments
        INTEGER(HID_T), INTENT(IN) :: group_id
        CHARACTER(LEN=*), INTENT(IN) :: name
        CLASS(*), INTENT(in) :: stencils(:)

        INTERFACE
            SUBROUTINE get_ptr(ptr, igrid)
                USE precision_mod, ONLY: intk
                INTEGER(intk), INTENT(OUT) :: ptr
                INTEGER(intk), INTENT(IN) :: igrid
            END SUBROUTINE get_ptr

            SUBROUTINE get_len(len, igrid)
                USE precision_mod, ONLY: intk
                INTEGER(intk), INTENT(OUT) :: len
                INTEGER(intk), INTENT(IN) :: igrid
            END SUBROUTINE get_len
        END INTERFACE
        OPTIONAL :: get_ptr, get_len

        LOGICAL, OPTIONAL :: indexlist
        LOGICAL, OPTIONAL :: extend
        INTEGER(intk), OPTIONAL :: ncmp
        LOGICAL, OPTIONAL :: same_kind

        ! Local variables
        INTEGER(kind=intk), ALLOCATABLE :: stencilInfo(:, :), &
            grpStencilInfo(:,:), totStencilInfo(:,:)
        INTEGER(kind=int64), ALLOCATABLE :: grdOffsets(:)
        INTEGER(kind=int64), ALLOCATABLE :: idx64(:)
        INTEGER(kind=intk), ALLOCATABLE :: idx(:)
        INTEGER(kind=intk) :: i, j, igrid
        INTEGER(kind=intk) :: nGrpGrids, nTotGrids
        INTEGER(kind=int64) :: my_len, ioproc_len, total_len

        INTEGER(kind=int32), ALLOCATABLE :: array_of_blocklengths(:), &
            array_of_displacements(:), nblocks(:)

        ! MPI_08 types
        ! TYPE(MPI_Datatype), ALLOCATABLE :: recvTypes(:)
        ! TYPE(MPI_Request), ALLOCATABLE :: recvReqs(:)
        ! TYPE(MPI_Datatype) :: mpi_dtype
        INTEGER(int32), ALLOCATABLE :: recvTypes(:)
        INTEGER(int32), ALLOCATABLE :: recvReqs(:)
        INTEGER(int32) :: mpi_dtype

        INTEGER(intk) :: len, ptr
        INTEGER(intk) :: rank, r, ncmp2
        INTEGER(intk) :: maxarr
        INTEGER(HSIZE_T) :: shape1(1), offset1(1), count1(1)
        INTEGER(HSIZE_T) :: shape2(2), curr_shape2(2)
        INTEGER(HSIZE_T) :: count2(2), offset2(2)
        INTEGER(HSIZE_T) :: maxdims2(2), chunksize2(2)
        INTEGER(HID_T) :: dset_id, plist_id, filespace, memspace
        INTEGER(HID_T) :: hdf5_memtype, hdf5_filetype
        TYPE(C_PTR) :: cptr

        ! Local per-node data buffer
        INTEGER(kind=intk), ALLOCATABLE, TARGET :: localbuf_i(:)
        REAL(kind=realk), ALLOCATABLE, TARGET :: localbuf_r(:)
        CLASS(*), POINTER :: localbuf(:)

        LOGICAL :: finishBlock, write_indexlist, append, link_exists
        INTEGER(kind=int32) :: ierr

        SELECT TYPE (stencils)
        TYPE IS (int_stencils_t)
            hdf5_memtype = mglet_hdf5_int
            mpi_dtype = mglet_mpi_int%MPI_val

            ! Determine max int value
            maxarr = 0
            DO i = 1, nMyGrids
                IF (ALLOCATED(stencils(i)%arr)) THEN
                    IF (SIZE(stencils(i)%arr) > 0) THEN
                        maxarr = MAX(maxarr, MAXVAL(ABS(stencils(i)%arr)))
                    END IF
                END IF
            END DO
            CALL MPI_Allreduce(MPI_IN_PLACE, maxarr, 1, mpi_dtype, &
                MPI_MAX, MPI_COMM_WORLD, ierr)
            IF (maxarr <= HUGE(1_int16)) THEN
                IF (PRESENT(same_kind)) THEN
                    IF (same_kind) THEN
                        hdf5_filetype = mglet_hdf5_int
                    ELSE
                        hdf5_filetype = h5kind_to_type(int16, H5_INTEGER_KIND)
                    END IF
                ELSE
                    hdf5_filetype = h5kind_to_type(int16, H5_INTEGER_KIND)
                END IF
            ELSE
                hdf5_filetype = mglet_hdf5_int
            END IF
        TYPE IS (real_stencils_t)
            hdf5_memtype = mglet_hdf5_real
            hdf5_filetype = mglet_hdf5_real
            mpi_dtype = mglet_mpi_real%MPI_val
        TYPE IS (INTEGER(intk))
            IF (.NOT. PRESENT(get_len)) THEN
                CALL errr(__FILE__, __LINE__)
            END IF
            IF (.NOT. PRESENT(get_ptr)) THEN
                CALL errr(__FILE__, __LINE__)
            END IF

            hdf5_memtype = mglet_hdf5_int
            mpi_dtype = mglet_mpi_int%MPI_val

            ! Determine max int value
            maxarr = 0.0
            DO i = 1, nMyGrids
                igrid = mygrids(i)
                CALL get_len(len, igrid)
                IF (len > 0) THEN
                    CALL get_ptr(ptr, igrid)
                    maxarr = MAX(maxarr, &
                        MAXVAL(ABS(stencils(ptr:ptr+len-1))))
                END IF
            END DO
            CALL MPI_Allreduce(MPI_IN_PLACE, maxarr, 1, mpi_dtype, &
                MPI_MAX, MPI_COMM_WORLD, ierr)
            IF (maxarr <= HUGE(1_int16)) THEN
                IF (PRESENT(same_kind)) THEN
                    IF (same_kind) THEN
                        hdf5_filetype = mglet_hdf5_int
                    ELSE
                        hdf5_filetype = h5kind_to_type(int16, H5_INTEGER_KIND)
                    END IF
                ELSE
                    hdf5_filetype = h5kind_to_type(int16, H5_INTEGER_KIND)
                END IF
            ELSE
                hdf5_filetype = mglet_hdf5_int
            END IF
        TYPE IS (REAL(realk))
            IF (.NOT. PRESENT(get_len)) THEN
                CALL errr(__FILE__, __LINE__)
            END IF
            IF (.NOT. PRESENT(get_ptr)) THEN
                CALL errr(__FILE__, __LINE__)
            END IF

            hdf5_memtype = mglet_hdf5_real
            hdf5_filetype = mglet_hdf5_real
            mpi_dtype = mglet_mpi_real%MPI_val
        CLASS DEFAULT
            CALL errr(__FILE__, __LINE__)
        END SELECT

        ! Number of grids in IO group
        nGrpGrids = 0
        CALL MPI_Allreduce(nMyGrids, nGrpGrids, 1, mglet_mpi_int%MPI_val, &
            MPI_SUM, iogrcomm%MPI_val, ierr)

        ! Allocate memory
        ALLOCATE(stencilInfo(3, nMyGrids))
        ALLOCATE(grpStencilInfo(3, nGrpGrids))

        ! These arrays are used for several purposes, allocate to the longest
        ! length and check to avoid the need to deallocate/allocate between
        ! usage scenarios
        ALLOCATE(array_of_blocklengths(MAX(nGrpGrids, iogrprocs, ioprocs)))
        ALLOCATE(array_of_displacements(MAX(nGrpGrids, iogrprocs, ioprocs)))

        ! Gather and compute offset
        array_of_blocklengths = 0
        CALL MPI_Gather(nMyGrids*3, 1, mglet_mpi_int%MPI_val, &
            array_of_blocklengths, 1, mglet_mpi_int%MPI_val, 0, &
            iogrcomm%MPI_val, ierr)

        array_of_displacements = 0
        DO i = 2, iogrprocs
            array_of_displacements(i) = &
                array_of_displacements(i-1) + array_of_blocklengths(i-1)
        END DO

        IF (MINVAL(array_of_displacements) < 0) THEN
            ! check the preceding summation loop (possible source of overflow)
            WRITE(*,*) "Probably integer overflow in array_of_displacements(:)"
            CALL errr(__FILE__, __LINE__)
        END IF

        ! Each rank allocate and fill in the stencil-info-table with the
        ! following:
        !   1: iogroup-comm rank
        !   2: grid ID
        !   3: number of stencil entries
        ! in total 3*nMygrids entries
        DO i = 1, nMyGrids
            igrid = mygrids(i)
            stencilInfo(1, i) = iogrid
            stencilInfo(2, i) = igrid
            SELECT TYPE (stencils)
            TYPE IS (int_stencils_t)
                IF (ALLOCATED(stencils(i)%arr)) THEN
                    stencilInfo(3, i) = SIZE(stencils(i)%arr)
                ELSE
                    stencilInfo(3, i) = 0
                END IF
            TYPE IS (real_stencils_t)
                IF (ALLOCATED(stencils(i)%arr)) THEN
                    stencilInfo(3, i) = SIZE(stencils(i)%arr)
                ELSE
                    stencilInfo(3, i) = 0
                END IF
            CLASS DEFAULT
                CALL get_len(len, igrid)
                stencilInfo(3, i) = len
            END SELECT
        END DO

        ! Write Index table
        write_indexlist = .TRUE.
        IF (PRESENT(indexlist)) THEN
            IF (.NOT. indexlist) THEN
                write_indexlist = .FALSE.
            END IF
        END IF
        IF (write_indexlist) THEN
            ALLOCATE(idx64(nMyGrids))
            idx64 = stencilInfo(3, :)
            CALL stencilio_write_list(group_id, TRIM(name)//"Index", &
                idx64, index=.TRUE.)
            DEALLOCATE(idx64)
        END IF

        ! Collect list
        CALL MPI_Gatherv(stencilInfo, INT(3*nMyGrids, int32), &
            mglet_mpi_int%MPI_val, grpStencilInfo, array_of_blocklengths, &
            array_of_displacements, mglet_mpi_int%MPI_val, 0, &
            iogrcomm%MPI_val, ierr)

        ! Each node determine own stencil data length, allocate data and copy
        my_len = 0
        DO i = 1, nMyGrids
            my_len = my_len + stencilInfo(3, i)
        END DO

        SELECT TYPE (stencils)
        TYPE IS (int_stencils_t)
            ALLOCATE(localbuf_i(my_len))
            localbuf => localbuf_i
        TYPE IS (real_stencils_t)
            ALLOCATE(localbuf_r(my_len))
            localbuf => localbuf_r
        TYPE IS (INTEGER(intk))
            ALLOCATE(localbuf_i(my_len))
            localbuf => localbuf_i
        TYPE IS (REAL(realk))
            ALLOCATE(localbuf_r(my_len))
            localbuf => localbuf_r
        END SELECT

        my_len = 0
        DO i = 1, nMyGrids
            igrid = mygrids(i)
            SELECT TYPE (stencils)
            TYPE IS (int_stencils_t)
                IF (ALLOCATED(stencils(i)%arr)) THEN
                    localbuf_i(my_len+1:my_len+SIZE(stencils(i)%arr)) = &
                        stencils(i)%arr(:)
                    my_len = my_len + SIZE(stencils(i)%arr)
                END IF
            TYPE IS (real_stencils_t)
                IF (ALLOCATED(stencils(i)%arr)) THEN
                    localbuf_r(my_len+1:my_len+SIZE(stencils(i)%arr)) = &
                        stencils(i)%arr(:)
                    my_len = my_len + SIZE(stencils(i)%arr)
                END IF
            TYPE IS (INTEGER(intk))
                CALL get_len(len, igrid)
                CALL get_ptr(ptr, igrid)
                localbuf_i(my_len+1:my_len+len) = stencils(ptr:ptr+len-1)
                my_len = my_len + len
            TYPE IS (REAL(realk))
                CALL get_len(len, igrid)
                CALL get_ptr(ptr, igrid)
                localbuf_r(my_len+1:my_len+len) = stencils(ptr:ptr+len-1)
                my_len = my_len + len
            END SELECT
        END DO

        IF (my_len < 0) THEN
            ! check the preceding summation loop (possible source of overflow)
            WRITE(*,*) "Probably integer overflow in my_len"
            CALL errr(__FILE__, __LINE__)
        END IF

        ! The IO rank determine the length needed to fetch the data from memory
        IF (ioProc) THEN
            ioproc_len = 0
            DO i = 1, nGrpGrids
                ioproc_len = ioproc_len + grpStencilInfo(3, i)
            END DO
            SELECT TYPE (stencils)
            TYPE IS (int_stencils_t)
                IF (ioproc_len > SIZE(intbuf, KIND=int64)) THEN
                    CALL increase_intbuf(ioproc_len)
                END IF
                cptr = C_LOC(intbuf)
            TYPE IS (INTEGER(intk))
                IF (ioproc_len > SIZE(intbuf, KIND=int64)) THEN
                    CALL increase_intbuf(ioproc_len)
                END IF
                cptr = C_LOC(intbuf)
            TYPE IS (real_stencils_t)
                IF (ioproc_len > SIZE(bigbuf, KIND=int64)) THEN
                    CALL increase_bigbuf(ioproc_len)
                END IF
                cptr = C_LOC(bigbuf)
            TYPE IS (REAL(realk))
                IF (ioproc_len > SIZE(bigbuf, KIND=int64)) THEN
                    CALL increase_bigbuf(ioproc_len)
                END IF
                cptr = C_LOC(bigbuf)
            END SELECT
        END IF

        ! The IO rank create derived MPI types for the receive operations
        ! and call MPI_Irecv
        IF (ioProc) THEN
            ! Determine sorted order of grids
            ALLOCATE(idx(nGrpGrids))
            CALL sortix(nGrpGrids, grpStencilInfo(2, :), idx)

            ! Determine offset of each grid
            !
            ! grdOffsets is checked to not overflow a 32-bit integer, since it
            ! later is downcasted to 32-bit in the array_of_displacements.
            ! If the overall stencil length per IO node exceed what can be
            ! put in a 32-bit integer, then add more IO nodes...
            ALLOCATE(grdOffsets(nGrpGrids))
            offset2(1) = 0
            grdOffsets = 0
            DO i = 1, nGrpGrids
                j = idx(i)
                grdOffsets(j) = offset2(1)
                offset2(1) = offset2(1) + grpStencilInfo(3, j)

                ! Check that grdOffsets does not overflow in next loop
                ! iteration
                !
                ! offset(1) is of kind HSIZE_T, which on all MGLET-relevant
                ! platforms is 64-bit
                IF (offset2(1) > HUGE(1_intk)) THEN
                    WRITE(*,*) "Integer overflow:"
                    WRITE(*,*) "Offset of stencils: ", offset2(1)
                    WRITE(*,*) "Exceed: ", HUGE(1_int32)
                    CALL errr(__FILE__, __LINE__)
                END IF
            END DO

            ALLOCATE(recvTypes(0:iogrprocs-1))
            ALLOCATE(nblocks(0:iogrprocs-1))

            nblocks = 0
            DO i = 0, iogrprocs-1
                nblocks(i) = 0
                DO j = 1, nGrpGrids
                    ! Only create a block if length > 0
                    IF (i == grpStencilInfo(1, j) .AND. &
                        grpStencilInfo(3, j) > 0) THEN
                        nblocks(i) = nblocks(i) + 1
                        array_of_blocklengths(nblocks(i)) = grpStencilInfo(3, j)
                        array_of_displacements(nblocks(i)) = &
                            INT(grdOffsets(j), int32)
                    END IF
                END DO
                IF (nblocks(i) > 0) THEN
                    CALL MPI_Type_indexed(nblocks(i), array_of_blocklengths, &
                        array_of_displacements, mpi_dtype, recvTypes(i), ierr)
                    CALL MPI_Type_commit(recvTypes(i), ierr)
                END IF
            END DO

            ALLOCATE(recvReqs(0:iogrprocs-1))
            recvReqs = MPI_REQUEST_NULL
            DO i = 0, iogrprocs-1
                IF (nblocks(i) > 0) THEN
                    SELECT TYPE (stencils)
                    TYPE IS (int_stencils_t)
                        CALL MPI_Irecv(intbuf, 1, recvTypes(i), i, 1, &
                            iogrcomm%MPI_val, recvReqs(i), ierr)
                    TYPE IS (INTEGER(intk))
                        CALL MPI_Irecv(intbuf, 1, recvTypes(i), i, 1, &
                            iogrcomm%MPI_val, recvReqs(i), ierr)
                    TYPE IS (real_stencils_t)
                        CALL MPI_Irecv(bigbuf, 1, recvTypes(i), i, 1, &
                            iogrcomm%MPI_val, recvReqs(i), ierr)
                    TYPE IS (REAL(realk))
                        CALL MPI_Irecv(bigbuf, 1, recvTypes(i), i, 1, &
                            iogrcomm%MPI_val, recvReqs(i), ierr)
                    END SELECT
                END IF
            END DO

        END IF

        ! All processes call MPI_Send if they have data to send
        ! No need to use non-blocing comm here, since there is only one
        ! send-call to do and nothing else to do afterwards.
        IF (my_len > 0) THEN
            IF (my_len > HUGE(1_int32)) THEN
                write(*,*) "my_len: ", my_len
                write(*,*) "exceed: ", HUGE(1_int32)
                CALL errr(__FILE__, __LINE__)
            END IF
            CALL MPI_Send(localbuf, INT(my_len, int32), mpi_dtype, 0, 1, &
                iogrcomm%MPI_val, ierr)
        END IF

        ! All processes can now deallocate some memory
        NULLIFY(localbuf)
        SELECT TYPE (stencils)
        TYPE IS (int_stencils_t)
            DEALLOCATE(localbuf_i)
        TYPE IS (real_stencils_t)
            DEALLOCATE(localbuf_r)
        TYPE IS (INTEGER(intk))
            DEALLOCATE(localbuf_i)
        TYPE IS (REAL(realk))
            DEALLOCATE(localbuf_r)
        END SELECT
        DEALLOCATE(stencilInfo)

        ! In order to assure that all processes have deallocated localbuf
        ! before proceeding, let's have an MPI_Barrier here.
        CALL MPI_Barrier(iogrcomm%MPI_val, ierr)

        ! Non-io processes can now return
        IF (.NOT. ioproc) THEN
            DEALLOCATE(grpStencilInfo)
            DEALLOCATE(array_of_blocklengths)
            DEALLOCATE(array_of_displacements)
            RETURN
        END IF

        ! --------------
        ! From here on only IO processes work to write out the data
        ! --------------
        ! Wait for the receive buffer to be ready
        CALL MPI_Waitall(iogrprocs, recvReqs, MPI_STATUSES_IGNORE, ierr)

        ! IO processes can free MPI datatypes
        DO i = 0, iogrprocs-1
            IF (nblocks(i) > 0) THEN
                CALL MPI_Type_free(recvTypes(i), ierr)
            END IF
        END DO
        IF (ioProc) THEN
            DEALLOCATE(nblocks)
            DEALLOCATE(recvTypes)
            DEALLOCATE(recvReqs)
        END IF

        ! Bigbuf now contain the nessecary data to be written out to disk.
        ! prepare HDF5 datatypes and write data out

        ! Total number of elements to be written
        CALL MPI_Allreduce(ioproc_len, total_len, 1, MPI_INTEGER8, MPI_SUM, &
            iocomm%MPI_val, ierr)

        rank = 1
        shape2(1) = total_len
        shape2(2) = 0

        ! Optionally support organizing the data in 2-D datasets
        ncmp2 = 1
        IF (PRESENT(ncmp)) THEN
            IF (MOD(total_len, INT(ncmp, int64)) > 0) THEN
                WRITE(*,*) total_len, ncmp
                CALL errr(__FILE__, __LINE__)
            END IF
            IF (ncmp > 1) THEN
                rank = 2
                shape2(1) = ncmp
                shape2(2) = total_len/ncmp
            END IF
            ncmp2 = ncmp
        END IF

        append = .FALSE.
        IF (PRESENT(extend)) THEN
            IF (extend) THEN
                append = .TRUE.
            END IF
        END IF

        ! Open/create dataset and filespace
        curr_shape2 = 0
        IF (append) THEN
            CALL hdf5common_dataset_exists(name, group_id, link_exists, &
                curr_shape2(1:rank))
            IF (link_exists) THEN
                CALL hdf5common_dataset_open(name, curr_shape2(1:rank), &
                    group_id, dset_id, filespace)
                shape2(rank) = shape2(rank) + curr_shape2(rank)
                CALL hdf5common_dataset_extend(dset_id, filespace, &
                    shape2(1:rank))
            ELSE
                maxdims2 = shape2
                maxdims2(rank) = H5S_UNLIMITED_F
                chunksize2 = shape2
                CALL hdf5common_dataset_create(name, shape2(1:rank), &
                    hdf5_filetype, group_id, dset_id, filespace, &
                    maxdims2(1:rank), chunksize2(1:rank))
            END IF
        ELSE
            CALL hdf5common_dataset_create(name, shape2(1:rank), &
                hdf5_filetype, group_id, dset_id, filespace)
        END IF

        ! This is to mitigate an assert-error in HDF5 when the dataset have no
        ! data. There is no need to continue in this case.
        IF (total_len == 0) THEN
            CALL hdf5common_dataset_close(dset_id, filespace)
            RETURN
        END IF

        ! Create memspace and select entire memspace as one hyperslab
        shape1 = ioproc_len
        CALL h5screate_simple_f(1, shape1, memspace, ierr)
        IF (ierr /= 0) CALL errr(__FILE__, __LINE__)

        offset1 = 0
        count1 = ioproc_len
        CALL h5sselect_hyperslab_f(memspace, H5S_SELECT_SET_F, offset1, &
            count1, ierr)
        IF (ierr /= 0) CALL errr(__FILE__, __LINE__)

        ! select hyperslab in filespace (already exists)
        CALL h5sselect_none_f(filespace, ierr)
        IF (ierr /= 0) CALL errr(__FILE__, __LINE__)

        CALL MPI_Allreduce(nGrpGrids, nTotGrids, 1, mglet_mpi_int%MPI_val, &
            MPI_SUM, iocomm%MPI_val, ierr)

        ! Set rank of IO process in IO comm in group information
        DO i = 1, nGrpGrids
            grpStencilInfo(1, i) = ioid
        END DO

        ! Blocklength and displacements for totStencilInfo array
        array_of_blocklengths = 0
        CALL MPI_Allgather(nGrpGrids*3, 1, mglet_mpi_int%MPI_val, &
            array_of_blocklengths, 1, mglet_mpi_int%MPI_val, iocomm%MPI_val, &
            ierr)
        array_of_displacements = 0
        DO i = 2, ioprocs
            array_of_displacements(i) = &
                array_of_displacements(i-1) + array_of_blocklengths(i-1)
        END DO

        ALLOCATE(totStencilInfo(3, nTotGrids))
        CALL MPI_Allgatherv(grpStencilInfo, INT(3*nGrpGrids, int32), &
            mglet_mpi_int%MPI_val, totStencilInfo, array_of_blocklengths, &
            array_of_displacements, mglet_mpi_int%MPI_val, iocomm%MPI_val, ierr)

        DEALLOCATE(grdOffsets)
        DEALLOCATE(idx)

        ALLOCATE(grdOffsets(nTotGrids))
        ALLOCATE(idx(nTotGrids))

        ! Determine offset for each individual grid
        CALL sortix(nTotGrids, totStencilInfo(2, :), idx)

        offset2(1) = 0
        grdOffsets = 0
        DO i = 1, nTotGrids
            j = idx(i)
            igrid = totStencilInfo(2, j)
            IF (igrid > nTotGrids) CALL errr(__FILE__, __LINE__)
            grdOffsets(igrid) = offset2(1)
            offset2(1) = offset2(1) + totStencilInfo(3, j)
        END DO

        IF (MINVAL(grdOffsets) < 0) THEN
            ! check the preceding summation loop (possible source of overflow)
            WRITE(*,*) "Probably integer overflow in grdOffsets(:)"
            CALL errr(__FILE__, __LINE__)
        END IF

        finishBlock = .FALSE.
        offset2 = 0
        count2 = 0
        IF (rank == 1) THEN
            r = 1
        ELSE
            r = 2
            count2(1) = ncmp2
        END IF
        DO i = 1, nGrpGrids
            igrid = grpStencilInfo(2, i)

            ! Start new block?
            IF (count2(r) == 0) THEN
                offset2(r) = curr_shape2(r) + grdOffsets(igrid)/ncmp2
                count2(r) = 0
            END IF

            count2(r) = count2(r) + grpStencilInfo(3, i)/ncmp2

            ! Finish current block
            IF (i < nGrpGrids) THEN
                IF (grpStencilInfo(2, i+1) /= igrid + 1) THEN
                    finishBlock = .TRUE.
                END IF
            ELSE
                finishBlock = .TRUE.
            END IF

            ! Define hyperslab
            IF ( finishBlock ) THEN
                IF ( count2(r) < 0 ) THEN
                    WRITE(*,*) "Probably integer overflow in count(1)"
                    CALL errr(__FILE__, __LINE__)
                END IF
                IF ( offset2(r) < 0 ) THEN
                    WRITE(*,*) "Probably integer overflow in offset(1)"
                    CALL errr(__FILE__, __LINE__)
                END IF
                IF (count2(r) > 0) THEN
                    CALL h5sselect_hyperslab_f(filespace, H5S_SELECT_OR_F, &
                        offset2, count2, ierr)
                    IF (ierr /= 0) CALL errr(__FILE__, __LINE__)
                END IF

                offset2(r) = 0
                count2(r) = 0
                finishBlock = .FALSE.
            END IF
        END DO

        ! Property list for collective dataset write
        CALL h5pcreate_f(H5P_DATASET_XFER_F, plist_id, ierr)
        IF (ierr /= 0) CALL errr(__FILE__, __LINE__)

        CALL h5pset_dxpl_mpio_f(plist_id, hdf5_io_mode, ierr)
        IF (ierr /= 0) CALL errr(__FILE__, __LINE__)

        ! write dataset collectively
        CALL h5dwrite_f(dset_id, hdf5_memtype, cptr, ierr, &
            file_space_id=filespace, mem_space_id=memspace, &
            xfer_prp=plist_id)
        IF (ierr /= 0) CALL errr(__FILE__, __LINE__)

        ! close property list for collective dataset write
        CALL h5pclose_f(plist_id, ierr)
        IF (ierr /= 0) CALL errr(__FILE__, __LINE__)

        ! Close memspace
        CALL h5sclose_f(memspace, ierr)
        IF (ierr /= 0) CALL errr(__FILE__, __LINE__)

        ! Close filespace and dataset
        CALL hdf5common_dataset_close(dset_id, filespace)

        DEALLOCATE(grpStencilInfo)
        DEALLOCATE(totStencilInfo)
        DEALLOCATE(grdOffsets)
        DEALLOCATE(idx)
        DEALLOCATE(array_of_blocklengths)
        DEALLOCATE(array_of_displacements)
    END SUBROUTINE stencilio_write


    SUBROUTINE stencilio_read(group_id, name, stencils, get_ptr, get_len)
        USE hdf5common_mod, ONLY: hdf5common_dataset_open, &
            hdf5common_dataset_close
        USE grids_mod, ONLY: mygrids, nmygrids
        USE commbuf_mod, ONLY: bigbuf, intbuf, increase_bigbuf, increase_intbuf

        ! Subroutine arguments
        INTEGER(HID_T), INTENT(IN) :: group_id
        CHARACTER(LEN=*), INTENT(IN) :: name
        CLASS(*), INTENT(inout) :: stencils(:)

        INTERFACE
            SUBROUTINE get_ptr(ptr, igrid)
                USE precision_mod, ONLY: intk
                INTEGER(intk), INTENT(OUT) :: ptr
                INTEGER(intk), INTENT(IN) :: igrid
            END SUBROUTINE get_ptr

            SUBROUTINE get_len(len, igrid)
                USE precision_mod, ONLY: intk
                INTEGER(intk), INTENT(OUT) :: len
                INTEGER(intk), INTENT(IN) :: igrid
            END SUBROUTINE get_len
        END INTERFACE
        OPTIONAL :: get_ptr, get_len

        ! Local variables
        INTEGER(kind=intk), ALLOCATABLE :: stencilInfo(:, :), &
            grpStencilInfo(:,:)
        INTEGER(kind=int64), ALLOCATABLE :: file_offset(:)
        INTEGER(kind=intk), ALLOCATABLE :: file_length(:)
        INTEGER(kind=intk), ALLOCATABLE :: idx(:), nGridsProc(:), displ(:), &
            gridid(:)
        INTEGER(kind=intk) :: i, j, igrid, iproc
        INTEGER(kind=intk) :: nGrpGrids, nTotGrids
        INTEGER(kind=int32) :: offset32, length
        INTEGER(kind=int64) :: ioproc_len

        INTEGER(kind=int32), ALLOCATABLE :: displs(:), sendcounts(:)
        ! TYPE(MPI_Request), ALLOCATABLE :: recvReqs(:), sendReqs(:)
        ! TYPE(MPI_Datatype) :: mpi_dtype
        INTEGER(int32), ALLOCATABLE :: recvReqs(:), sendReqs(:)
        INTEGER(int32) :: mpi_dtype

        INTEGER(intk) :: len, ptr
        INTEGER :: rank
        INTEGER(HSIZE_T) :: shape(1)
        INTEGER(HSIZE_T) :: count(1), offset(1)
        INTEGER(HID_T) :: dset_id, plist_id, filespace, memspace

        INTEGER(HID_T) :: hdf5_dtype

        LOGICAL :: finishBlock
        INTEGER(kind=int32) :: ierr

        TYPE(C_PTR) :: cptr

        SELECT TYPE (stencils)
        TYPE IS (int_stencils_t)
            hdf5_dtype = mglet_hdf5_int
            mpi_dtype = mglet_mpi_int%MPI_val
        TYPE IS (real_stencils_t)
            hdf5_dtype = mglet_hdf5_real
            mpi_dtype = mglet_mpi_real%MPI_val
        TYPE IS (INTEGER(intk))
            IF (.NOT. PRESENT(get_len)) THEN
                CALL errr(__FILE__, __LINE__)
            END IF
            IF (.NOT. PRESENT(get_ptr)) THEN
                CALL errr(__FILE__, __LINE__)
            END IF
            hdf5_dtype = mglet_hdf5_int
            mpi_dtype = mglet_mpi_int%MPI_val
        TYPE IS (REAL(realk))
            IF (.NOT. PRESENT(get_len)) THEN
                CALL errr(__FILE__, __LINE__)
            END IF
            IF (.NOT. PRESENT(get_ptr)) THEN
                CALL errr(__FILE__, __LINE__)
            END IF
            hdf5_dtype = mglet_hdf5_real
            mpi_dtype = mglet_mpi_real%MPI_val
        CLASS DEFAULT
            CALL errr(__FILE__, __LINE__)
        END SELECT

        ! Number of grids in IO group
        nGrpGrids = 0
        CALL MPI_Allreduce(nMyGrids, nGrpGrids, 1, mglet_mpi_int%MPI_val, &
            MPI_SUM, iogrcomm%MPI_val, ierr)

        nTotGrids = 0
        IF (ioproc) THEN
            ! Overall number of grids reduced across all IO processes
            CALL MPI_Allreduce(nGrpGrids, nTotGrids, 1, mglet_mpi_int%MPI_val, &
                MPI_SUM, iocomm%MPI_val, ierr)

            ! Allocate offset-and-length arrays and read them in from Index table
            ALLOCATE(file_offset(nTotGrids))
            ALLOCATE(file_length(nTotGrids))
            CALL stencilio_read_indexlist(group_id, TRIM(name)//"Index", &
                file_offset, file_length)
        END IF

        ! On each IO process, gather list of grid ID's
        ALLOCATE(nGridsProc(iogrprocs))
        ALLOCATE(displ(iogrprocs))
        ALLOCATE(gridid(nGrpGrids))
        ALLOCATE(idx(nGrpGrids))
        nGridsProc = 0
        displ = 0
        gridid = -1

        CALL MPI_Gather(nMyGrids, 1, mglet_mpi_int%MPI_val, nGridsProc, 1, &
            mglet_mpi_int%MPI_val, 0, iogrcomm%MPI_val, ierr)

        DO i = 2, iogrprocs
            displ(i) = displ(i-1) + nGridsProc(i-1)
        END DO

        CALL MPI_Gatherv(myGrids, INT(nMyGrids, int32), mglet_mpi_int%MPI_val, &
            gridid, nGridsProc, displ, mglet_mpi_int%MPI_val, 0, &
            iogrcomm%MPI_val, ierr)

        DEALLOCATE(displ)

        CALL sortix(nGrpGrids, gridid, idx)

        ! Each IO process now have a list over the grids belonging to the IO
        ! group, and the corresponding offset- and count of these grids
        ! in the stencil list.
        !
        ! Start preparing data read
        IF (ioProc) THEN
            ! Open dataset and filespace
            CALL hdf5common_dataset_open(name, shape, group_id, dset_id, &
                filespace)

            ! This is to mitigate an assert-error in HDF5 when the dataset have no
            ! data. There is no need to continue in this case.
            IF (shape(1) == 0) THEN
                CALL hdf5common_dataset_close(dset_id, filespace)
            ELSE

                ! Create memspace and select entire memspace as one hyperslab
                ioproc_len = 0
                DO i = 1, nGrpGrids
                    ioproc_len = ioproc_len + file_length(gridid(idx(i)))
                END DO

                rank = 1
                shape(1) = ioproc_len
                CALL h5screate_simple_f(rank, shape, memspace, ierr)
                IF (ierr /= 0) CALL errr(__FILE__, __LINE__)

                offset(1) = 0
                count = ioproc_len
                CALL h5sselect_hyperslab_f(memspace, H5S_SELECT_SET_F, &
                    offset, count, ierr)
                IF (ierr /= 0) CALL errr(__FILE__, __LINE__)

                ! select hyperslab in filespace (already exists)
                CALL h5sselect_none_f(filespace, ierr)
                IF (ierr /= 0) CALL errr(__FILE__, __LINE__)

                finishBlock = .FALSE.
                offset = 0
                count = 0
                DO i = 1, nGrpGrids
                    igrid = gridid(idx(i))

                    ! Start new block?
                    IF (count(1) == 0) THEN
                        offset(1) = file_offset(igrid)
                        count(1) = 0
                    END IF

                    count(1) = count(1) + file_length(igrid)

                    ! Finish current block
                    IF (i < nGrpGrids) THEN
                        IF (gridid(idx(i+1)) /= igrid + 1) THEN
                            finishBlock = .TRUE.
                        END IF
                    ELSE
                        finishBlock = .TRUE.
                    END IF

                    ! Define hyperslab
                    IF (finishBlock) THEN
                        IF (count(1) > 0) THEN
                            CALL h5sselect_hyperslab_f(filespace, &
                                H5S_SELECT_OR_F, offset, count, ierr)
                            IF (ierr /= 0) CALL errr(__FILE__, __LINE__)
                        END IF

                        offset = 0
                        count = 0
                        finishBlock = .FALSE.
                    END IF
                END DO

                ! Property list for collective dataset write
                CALL h5pcreate_f(H5P_DATASET_XFER_F, plist_id, ierr)
                IF (ierr /= 0) CALL errr(__FILE__, __LINE__)

                CALL h5pset_dxpl_mpio_f(plist_id, hdf5_io_mode, ierr)
                IF (ierr /= 0) CALL errr(__FILE__, __LINE__)

                ! Make sure bigfuf is large enough
                SELECT TYPE (stencils)
                    TYPE IS (int_stencils_t)
                        IF (ioproc_len > SIZE(intbuf, KIND=int64)) THEN
                            CALL increase_intbuf(ioproc_len)
                        END IF
                        cptr = C_LOC(intbuf)
                    TYPE IS (INTEGER(intk))
                        IF (ioproc_len > SIZE(intbuf, KIND=int64)) THEN
                            CALL increase_intbuf(ioproc_len)
                        END IF
                        cptr = C_LOC(intbuf)
                    TYPE IS (real_stencils_t)
                        IF (ioproc_len > SIZE(bigbuf, KIND=int64)) THEN
                            CALL increase_bigbuf(ioproc_len)
                        END IF
                        cptr = C_LOC(bigbuf)
                    TYPE IS (REAL(realk))
                        IF (ioproc_len > SIZE(bigbuf, KIND=int64)) THEN
                            CALL increase_bigbuf(ioproc_len)
                        END IF
                        cptr = C_LOC(bigbuf)
                END SELECT

                ! Read dataset collectively
                CALL h5dread_f(dset_id, hdf5_dtype, cptr, ierr, &
                    file_space_id=filespace, mem_space_id=memspace, &
                    xfer_prp=plist_id)
                IF (ierr /= 0) CALL errr(__FILE__, __LINE__)

                ! Close filespace and dataset
                CALL hdf5common_dataset_close(dset_id, filespace)

                CALL h5sclose_f(memspace, ierr)
                IF (ierr /= 0) CALL errr(__FILE__, __LINE__)

                CALL h5pclose_f(plist_id, ierr)
                IF (ierr /= 0) CALL errr(__FILE__, __LINE__)
            END IF

            DEALLOCATE(file_offset)

        END IF

        ALLOCATE(stencilInfo(3, nMyGrids))
        ALLOCATE(grpStencilInfo(3, nGrpGrids))
        ALLOCATE(displs(iogrprocs))
        ALLOCATE(sendcounts(iogrprocs))
        stencilInfo = 0
        grpStencilInfo = 0
        displs = -1
        sendcounts = -1

        IF (ioproc) THEN
            j = 1
            DO i = 1, iogrprocs
                displs(i) = (j-1)*3
                sendcounts(i) = nGridsProc(i)*3
                DO j = j, j + nGridsProc(i) - 1
                    grpStencilInfo(1, j) = i-1
                    grpStencilInfo(2, j) = gridid(j)
                    grpStencilInfo(3, j) = file_length(gridid(j))
                END DO
            END DO
            DEALLOCATE(file_length)
        END IF

        DEALLOCATE(nGridsProc)
        DEALLOCATE(gridid)
        ! Scatter this list to IO group members
        CALL MPI_Scatterv(grpStencilInfo, sendcounts, displs, &
            mglet_mpi_int%MPI_val, stencilInfo, INT(nMyGrids*3, int32), &
            mglet_mpi_int%MPI_val, 0, iogrcomm%MPI_val, ierr)

        DEALLOCATE(displs)
        DEALLOCATE(sendcounts)

        ! Each rank allocate stencil array list and call MPI_Recv
        ALLOCATE(recvReqs(nMyGrids))
        ALLOCATE(sendReqs(nGrpGrids))
        recvReqs = MPI_REQUEST_NULL
        sendReqs = MPI_REQUEST_NULL

        DO i = 1, nMyGrids
            igrid = stencilInfo(2, i)
            length = stencilInfo(3, i)
            SELECT TYPE (stencils)
            TYPE IS (int_stencils_t)
                ALLOCATE(stencils(i)%arr(length))
                stencils(i)%arr = 0
                IF (length > 0) THEN
                    CALL MPI_Irecv(stencils(i)%arr, length, mpi_dtype, &
                        0, 0, iogrcomm%MPI_val, recvReqs(i), ierr)
                END IF
            TYPE IS (real_stencils_t)
                ALLOCATE(stencils(i)%arr(length))
                stencils(i)%arr = 0.0
                IF (length > 0) THEN
                    CALL MPI_Irecv(stencils(i)%arr, length, mpi_dtype, &
                        0, 0, iogrcomm%MPI_val, recvReqs(i), ierr)
                END IF
            TYPE IS (INTEGER(intk))
                CALL get_len(len, igrid)
                CALL get_ptr(ptr, igrid)

                IF (length /= len) THEN
                    CALL errr(__FILE__, __LINE__)
                END IF
                IF (length > 0) THEN
                    CALL MPI_Irecv(stencils(ptr), length, mpi_dtype, &
                        0, 0, iogrcomm%MPI_val, recvReqs(i), ierr)
                END IF
            TYPE IS (REAL(realk))
                CALL get_len(len, igrid)
                CALL get_ptr(ptr, igrid)

                IF (length /= len) THEN
                    CALL errr(__FILE__, __LINE__)
                END IF
                IF (length > 0) THEN
                    CALL MPI_Irecv(stencils(ptr), length, mpi_dtype, &
                        0, 0, iogrcomm%MPI_val, recvReqs(i), ierr)
                END IF
            END SELECT
        END DO

        ! IO processes send data
        IF (ioproc) THEN
            CALL sortix(nGrpGrids, grpStencilInfo(2, :), idx)
            offset32 = 1
            DO i = 1, nGrpGrids
                iproc = grpStencilInfo(1, idx(i))
                igrid = grpStencilInfo(2, idx(i))
                length = grpStencilInfo(3, idx(i))

                IF (length > 0) THEN
                    SELECT TYPE (stencils)
                    TYPE IS (int_stencils_t)
                        CALL MPI_Isend(intbuf(offset32), length, &
                            mpi_dtype, iproc, 0, iogrcomm%MPI_val, &
                            sendReqs(idx(i)), ierr)
                    TYPE IS (INTEGER(intk))
                        CALL MPI_Isend(intbuf(offset32), length, &
                            mpi_dtype, iproc, 0, iogrcomm%MPI_val, &
                            sendReqs(idx(i)), ierr)
                    TYPE IS (real_stencils_t)
                        CALL MPI_Isend(bigbuf(offset32), length, &
                            mpi_dtype, iproc, 0, iogrcomm%MPI_val, &
                            sendReqs(idx(i)), ierr)
                    TYPE IS (REAL(realk))
                        CALL MPI_Isend(bigbuf(offset32), length, &
                            mpi_dtype, iproc, 0, iogrcomm%MPI_val, &
                            sendReqs(idx(i)), ierr)
                    END SELECT
                END IF
                offset32 = offset32 + length
            END DO
        END IF

        DEALLOCATE(idx)
        DEALLOCATE(StencilInfo)
        DEALLOCATE(grpStencilInfo)

        ! Wait for Recv to finish
        CALL MPI_Waitall(INT(nMyGrids, int32), recvReqs, MPI_STATUSES_IGNORE, &
            ierr)
        IF (ioproc) THEN
            CALL MPI_Waitall(INT(nGrpGrids, int32), sendReqs, &
                MPI_STATUSES_IGNORE, ierr)
        END IF

        DEALLOCATE(sendreqs)
        DEALLOCATE(recvreqs)
    END SUBROUTINE stencilio_read


    SUBROUTINE stencilio_write_list(group_id, name, list, index)
        USE grids_mod, ONLY: mygrids, nmygrids

        ! Subroutine arguments
        INTEGER(HID_T), INTENT(IN) :: group_id
        CHARACTER(LEN=*), INTENT(IN) :: name
        CLASS(*), INTENT(in) :: list(:)
        LOGICAL, OPTIONAL, INTENT(in) :: index

        ! Local variables
        INTEGER(kind=intk), ALLOCATABLE :: idx(:), &
            nElemsProc(:), displ(:), gridid(:)
        INTEGER(kind=intk) :: i, nelemstot

        CLASS(*), ALLOCATABLE, TARGET :: datalist(:), sorted(:)

        INTEGER(kind=intk) :: rank
        INTEGER(kind=int64) :: indexOffset
        INTEGER(HSIZE_T) :: shape(1)
        INTEGER(HSIZE_T) :: count(1), offset(1)
        INTEGER(HID_T) :: dset_id, filespace, memspace, plist_id

        TYPE(C_PTR) :: cptr

        INTEGER(HID_T) :: hdf5_dtype
        ! TYPE(MPI_Datatype) :: mpi_dtype
        INTEGER(int32) :: mpi_dtype

        INTEGER(kind=int32) :: ierr

        SELECT TYPE (list)
        TYPE IS (INTEGER(kind=intk))
            hdf5_dtype = mglet_hdf5_int
            mpi_dtype = mglet_mpi_int%MPI_val
        TYPE IS (INTEGER(kind=int64))
            hdf5_dtype = h5kind_to_type(int64, H5_INTEGER_KIND)
            mpi_dtype = MPI_INTEGER8
        TYPE IS (REAL(kind=realk))
            hdf5_dtype = mglet_hdf5_real
            mpi_dtype = mglet_mpi_real%MPI_val
        CLASS DEFAULT
            CALL errr(__FILE__, __LINE__)
        END SELECT

        ! Each process communicate how many list elemnents they have
        ALLOCATE(nElemsProc(numprocs))
        ALLOCATE(displ(numprocs))
        nElemsProc = 0
        displ = 0

        ! Do an Allgather here since this data is used to compute nElemsTot,
        ! in this way everyone gets nElemsTot in the end
        CALL MPI_Allgather(nMyGrids, 1, mglet_mpi_int%MPI_val, &
            nElemsProc, 1, mglet_mpi_int%MPI_val, MPI_COMM_WORLD, ierr)

        ! Allocate more....
        nElemsTot = SUM(nElemsProc)
        ALLOCATE(gridid(nElemsTot))
        gridid = -1

        ! Then gather grid number on global rank 0
        ! Only rank 0 need this because only rank 0 receive the data in the
        ! end and do the sorting alone
        DO i = 2, numprocs
            displ(i) = displ(i-1) + nElemsProc(i-1)
        END DO
        CALL MPI_Gatherv(myGrids, INT(nMyGrids, int32), mglet_mpi_int%MPI_val, &
            gridid, nElemsProc, displ, mglet_mpi_int%MPI_val, 0, &
            MPI_COMM_WORLD, ierr)

        ! Allocate data
        SELECT TYPE (list)
        TYPE IS (INTEGER(kind=intk))
            ALLOCATE(INTEGER(kind=intk) :: datalist(nElemsTot))
            ALLOCATE(INTEGER(kind=intk) :: sorted(nElemsTot))
        TYPE IS (INTEGER(kind=int64))
            ALLOCATE(INTEGER(kind=int64) :: datalist(nElemsTot))
            ALLOCATE(INTEGER(kind=int64) :: sorted(nElemsTot))
        TYPE IS (REAL(kind=realk))
            ALLOCATE(REAL(kind=realk) :: datalist(nElemsTot))
            ALLOCATE(REAL(kind=realk) :: sorted(nElemsTot))
        END SELECT

        ! Collect data and store into datalist (make sure that )
        CALL MPI_Gatherv(list, nMyGrids, mpi_dtype, &
            datalist, nElemsProc, displ, mpi_dtype, 0, MPI_COMM_WORLD, ierr)

        ! Sort by grid ID
        IF (ioid == 0) THEN
            ALLOCATE(idx(nElemsTot))
            CALL sortix(nElemsTot, gridid, idx)
            SELECT TYPE (sorted)
            TYPE IS (INTEGER(kind=int64))
                SELECT TYPE (datalist)
                TYPE IS (INTEGER(kind=int64))
                    DO i = 1, nElemsTot
                        sorted(i) = datalist(idx(i))
                    END DO
                END SELECT
                cptr = C_LOC(sorted)

            TYPE IS (INTEGER(kind=intk))
                SELECT TYPE (datalist)
                TYPE IS (INTEGER(kind=intk))
                    DO i = 1, nElemsTot
                        sorted(i) = datalist(idx(i))
                    END DO
                END SELECT
                cptr = C_LOC(sorted)

            TYPE IS (REAL(kind=realk))
                SELECT TYPE (datalist)
                TYPE IS (REAL(kind=realk))
                    DO i = 1, nElemsTot
                        sorted(i) = datalist(idx(i))
                    END DO
                END SELECT
                cptr = C_LOC(sorted)
            END SELECT
            DEALLOCATE(idx)
        END IF

        DEALLOCATE(datalist)
        DEALLOCATE(gridid)
        DEALLOCATE(nelemsproc)
        DEALLOCATE(displ)

        ! Non-io processes can now return
        IF (.NOT. ioproc) THEN
            DEALLOCATE(sorted)
            RETURN
        END IF

        ! Create an index array that skips the first index which is always
        ! 0. That way, the last element is the size of the array. Element i
        ! therefore starts at integer_buffer(i-1) and has the length
        ! integer_buffer(i)-integer_buffer(i-1).
        IF (ioid == 0 .AND. PRESENT(index)) THEN
            IF (index) THEN
                SELECT TYPE (sorted)
                TYPE IS (INTEGER(kind=int64))
                    indexOffset = 0
                    DO i = 1, nElemsTot
                        indexOffset = indexOffset + sorted(i)
                        sorted(i) = indexOffset
                    END DO
                CLASS DEFAULT
                    CALL errr(__FILE__, __LINE__)
                END SELECT
            END IF
        END IF

        ! --------------
        ! From here on only IO processes work to write out the data
        ! --------------
        ! Open/create dataset and filespace
        shape(1) = nElemsTot
        CALL hdf5common_dataset_create(name, shape, hdf5_dtype, group_id, &
            dset_id, filespace)

        ! Create memspace
        rank = 1
        shape(1) = nElemsTot
        CALL h5screate_simple_f(rank, shape, memspace, ierr)
        IF (ierr /= 0) CALL errr(__FILE__, __LINE__)

        ! Only ioid = 0 do the actual write
        offset(1) = 0
        count = nElemsTot
        IF (ioid == 0) THEN
            CALL h5sselect_hyperslab_f(memspace, H5S_SELECT_SET_F, &
                offset, count, ierr)
            IF (ierr /= 0) CALL errr(__FILE__, __LINE__)

            CALL h5sselect_hyperslab_f(filespace, H5S_SELECT_SET_F, &
                offset, count, ierr)
            IF (ierr /= 0) CALL errr(__FILE__, __LINE__)
        ELSE
            CALL h5sselect_none_f(memspace, ierr)
            IF (ierr /= 0) CALL errr(__FILE__, __LINE__)

            CALL h5sselect_none_f(filespace, ierr)
            IF (ierr /= 0) CALL errr(__FILE__, __LINE__)
        END IF

        ! Property list for collective dataset write
        CALL h5pcreate_f(H5P_DATASET_XFER_F, plist_id, ierr)
        IF (ierr /= 0) CALL errr(__FILE__, __LINE__)

        ! Only the IO master rank 0 reads any data. It makes sense to
        ! always use independent IO here, and not use the hdf5_io_mode
        ! switch
        CALL h5pset_dxpl_mpio_f(plist_id, H5FD_MPIO_INDEPENDENT_F, ierr)
        IF (ierr /= 0) CALL errr(__FILE__, __LINE__)

        ! write dataset collectively
        CALL h5dwrite_f(dset_id, hdf5_dtype, cptr, ierr, &
            file_space_id=filespace, mem_space_id=memspace, &
            xfer_prp=plist_id)
        IF (ierr /= 0) CALL errr(__FILE__, __LINE__)

        ! close property list for collective dataset write
        CALL h5pclose_f(plist_id, ierr)
        IF (ierr /= 0) CALL errr(__FILE__, __LINE__)

        ! Close memspace
        CALL h5sclose_f(memspace, ierr)
        IF (ierr /= 0) CALL errr(__FILE__, __LINE__)

        ! Close filespace and dataset
        CALL hdf5common_dataset_close(dset_id, filespace)

        DEALLOCATE(sorted)
    END SUBROUTINE stencilio_write_list


    ! Get index list, collective on all processes of an IO group.
    ! All processes call this routine, and all processes are returned
    ! the offset and count of a grid.
    !
    ! The arrays offset and length must be allocated a priori.
    SUBROUTINE stencilio_read_list(group_id, name, list)
        USE grids_mod, ONLY: mygrids, nmygrids

        ! Subroutine arguments
        INTEGER(HID_T), INTENT(IN) :: group_id
        CHARACTER(LEN=*), INTENT(IN) :: name
        CLASS(*), INTENT(OUT) :: list(:)

        ! Local variables
        INTEGER(kind=intk), ALLOCATABLE :: gridid(:)
        INTEGER(kind=intk) :: i, nelemstot
        CLASS(*), ALLOCATABLE :: datalist(:), sorted(:)
        INTEGER(intk), ALLOCATABLE :: nElemsProc(:), displ(:)
        ! TYPE(MPI_Datatype) :: mpi_dtype
        INTEGER(int32) :: mpi_dtype
        INTEGER(int32) :: ierr

        SELECT TYPE (list)
        TYPE IS (INTEGER(kind=intk))
            mpi_dtype = mglet_mpi_int%MPI_val
        TYPE IS (REAL(kind=realk))
            mpi_dtype = mglet_mpi_real%MPI_val
        CLASS DEFAULT
            CALL errr(__FILE__, __LINE__)
        END SELECT

        ! Each process communicate how many list elemnents they have
        ALLOCATE(nElemsProc(numprocs))
        ALLOCATE(displ(numprocs))
        nElemsProc = 0
        displ = 0

        ! All processes need to know the length of the dataset
        CALL MPI_Allgather(nMyGrids, 1, mglet_mpi_int%MPI_val, &
            nElemsProc, 1, mglet_mpi_int%MPI_val, MPI_COMM_WORLD, ierr)

        ! Allocate more....
        nElemsTot = SUM(nElemsProc)
        ALLOCATE(gridid(nElemsTot))
        gridid = -1

        ! Then gather grid number
        DO i = 2, numprocs
            displ(i) = displ(i-1) + nElemsProc(i-1)
        END DO
        CALL MPI_Gatherv(myGrids, INT(nMyGrids, int32), mglet_mpi_int%MPI_val, &
            gridid, nElemsProc, displ, mglet_mpi_int%MPI_val, 0, &
            MPI_COMM_WORLD, ierr)

        ! Allocate data
        SELECT TYPE (list)
        TYPE IS (INTEGER(kind=intk))
            ALLOCATE(INTEGER(kind=intk) :: datalist(nElemsTot))
            ALLOCATE(INTEGER(kind=intk) :: sorted(nElemsTot))
        TYPE IS (REAL(kind=realk))
            ALLOCATE(REAL(kind=realk) :: datalist(nElemsTot))
            ALLOCATE(REAL(kind=realk) :: sorted(nElemsTot))
        END SELECT

        ! Read data - all IO processes now hold entire dataset
        CALL stencilio_read_master(group_id, name, datalist)

        ! Sort by grid ID
        IF (ioid == 0) THEN
            SELECT TYPE (sorted)
            TYPE IS (INTEGER(kind=intk))
                SELECT TYPE (datalist)
                TYPE IS (INTEGER(kind=intk))
                    DO i = 1, nElemsTot
                        sorted(i) = datalist(gridid(i))
                    END DO
                END SELECT

            TYPE IS (REAL(kind=realk))
                SELECT TYPE (datalist)
                TYPE IS (REAL(kind=realk))
                    DO i = 1, nElemsTot
                        sorted(i) = datalist(gridid(i))
                    END DO
                END SELECT
            END SELECT
        END IF

        ! Deallicate
        DEALLOCATE(datalist)
        DEALLOCATE(gridid)

        ! Scatter data
        SELECT TYPE (sorted)
        TYPE IS (INTEGER(kind=intk))
            SELECT TYPE (list)
            TYPE IS (INTEGER(kind=intk))
                CALL MPI_Scatterv(sorted, nElemsProc, displ, mpi_dtype, list, &
                    INT(nMyGrids, int32), mpi_dtype, 0, MPI_COMM_WORLD, ierr)
            END SELECT

        TYPE IS (REAL(kind=realk))
            SELECT TYPE (list)
                TYPE IS (REAL(kind=realk))
                CALL MPI_Scatterv(sorted, nElemsProc, displ, mpi_dtype, list, &
                    INT(nMyGrids, int32), mpi_dtype, 0, MPI_COMM_WORLD, ierr)
            END SELECT
        END SELECT

        DEALLOCATE(nElemsProc)
        DEALLOCATE(displ)
    END SUBROUTINE stencilio_read_list


    ! Get index list, collective on all processes of an IO group.
    ! All processes call this routine, and all processes are returned
    ! the offset and count of a grid.
    !
    ! The arrays offset and length must be allocated a priori.
    SUBROUTINE stencilio_read_indexlist(group_id, name, file_offset, file_length)
        ! Subroutine arguments
        INTEGER(HID_T), INTENT(IN) :: group_id
        CHARACTER(LEN=*), INTENT(IN) :: name
        INTEGER(kind=int64), INTENT(OUT) :: file_offset(:)
        INTEGER(kind=intk), INTENT(OUT) :: file_length(:)

        ! Local variables
        INTEGER(kind=intk) :: i, nelemstot
        INTEGER(kind=int64), ALLOCATABLE :: indexlist(:)

        IF (.NOT. ioProc) THEN
            RETURN
        END IF

        nElemsTot = SIZE(file_offset)
        ALLOCATE(indexlist(nElemsTot))
        ! in this step, all offset data is stored into indexlist
        CALL stencilio_read_collective(group_id, name, indexlist)

        ! Return offset and length
        file_offset = 0
        file_length = 0

        ! Length of first item - offset always 0
        file_offset(1) = INT(0, int64)
        file_length(1) = INT(indexlist(1), intk)
        DO i = 2, nElemsTot
            file_offset(i) = INT(indexlist(i-1), int64)
            file_length(i) = INT(indexlist(i)-indexlist(i-1), intk)
        END DO

        DEALLOCATE(indexlist)
    END SUBROUTINE stencilio_read_indexlist


    ! All IO processes read the same dataset collectively and all have the
    ! entire dataset in the end
    SUBROUTINE stencilio_read_collective(group_id, name, list)
        ! Subroutine arguments
        INTEGER(HID_T), INTENT(IN) :: group_id
        CHARACTER(LEN=*), INTENT(IN) :: name
        CLASS(*), INTENT(OUT) :: list(:)

        ! Local variables
        INTEGER(intk) :: rank, nelemstot
        INTEGER(HSIZE_T) :: shape(1)
        INTEGER(HSIZE_T) :: count(1), offset(1)
        INTEGER(HID_T) :: dset_id, filespace, memspace, plist_id
        INTEGER(HID_T) :: hdf5_dtype

        INTEGER(kind=int32) :: ierr

        IF (.NOT. ioProc) THEN
            RETURN
        END IF

        SELECT TYPE (list)
        TYPE IS (INTEGER(kind=int64))
            hdf5_dtype = h5kind_to_type(int64, H5_INTEGER_KIND)
            nElemsTot = SIZE(list)
        TYPE IS (REAL(kind=realk))
            hdf5_dtype = mglet_hdf5_real
            nElemsTot = SIZE(list)
        CLASS DEFAULT
            CALL errr(__FILE__, __LINE__)
        END SELECT

        ! Open dataset and filespace
        CALL hdf5common_dataset_open(name, shape, group_id, dset_id, filespace)

        IF (shape(1) /= nElemsTot) THEN
            WRITE(*,*) "Reading list, got wrong dimensions"
            WRITE(*,*) "  myid = ", myid
            WRITE(*,*) "  ioid = ", ioid
            WRITE(*,*) "  Length of dataset in file:", shape(1)
            WRITE(*,*) "  Length of allocated memory:", nElemsTot
            CALL errr(__FILE__, __LINE__)
        END IF

        ! Create memspace
        rank = 1
        shape(1) = nElemsTot
        CALL h5screate_simple_f(rank, shape, memspace, ierr)
        IF (ierr /= 0) CALL errr(__FILE__, __LINE__)

        ! All IO processes read the same data
        offset(1) = 0
        count = nElemsTot
        CALL h5sselect_hyperslab_f(memspace, H5S_SELECT_SET_F, &
            offset, count, ierr)
        IF (ierr /= 0) CALL errr(__FILE__, __LINE__)

        CALL h5sselect_hyperslab_f(filespace, H5S_SELECT_SET_F, &
            offset, count, ierr)
        IF (ierr /= 0) CALL errr(__FILE__, __LINE__)

        ! Property list for collective dataset read
        CALL h5pcreate_f(H5P_DATASET_XFER_F, plist_id, ierr)
        IF (ierr /= 0) CALL errr(__FILE__, __LINE__)

        CALL h5pset_dxpl_mpio_f(plist_id, hdf5_io_mode, ierr)
        IF (ierr /= 0) CALL errr(__FILE__, __LINE__)

        ! Read dataset collectively (data types are already set correctly)
        shape(1) = nElemsTot
        SELECT TYPE (list)
        TYPE IS (INTEGER(kind=int64))
            CALL h5dread_f(dset_id, hdf5_dtype, list, shape, ierr, &
                file_space_id=filespace, mem_space_id=memspace, &
                xfer_prp=plist_id)
            IF (ierr /= 0) CALL errr(__FILE__, __LINE__)
        TYPE IS (REAL(kind=realk))
            CALL h5dread_f(dset_id, hdf5_dtype, list, shape, ierr, &
                file_space_id=filespace, mem_space_id=memspace, &
                xfer_prp=plist_id)
            IF (ierr /= 0) CALL errr(__FILE__, __LINE__)
        END SELECT

        ! Close handlers
        CALL h5pclose_f(plist_id, ierr)
        IF (ierr /= 0) CALL errr(__FILE__, __LINE__)

        CALL h5sclose_f(memspace, ierr)
        IF (ierr /= 0) CALL errr(__FILE__, __LINE__)

        CALL hdf5common_dataset_close(dset_id, filespace)
    END SUBROUTINE stencilio_read_collective


    ! The IO master rank read a list. No other processes get any data.
    SUBROUTINE stencilio_read_master(group_id, name, list)
        ! Subroutine arguments
        INTEGER(HID_T), INTENT(IN) :: group_id
        CHARACTER(LEN=*), INTENT(IN) :: name
        CLASS(*), INTENT(INOUT), TARGET :: list(:)

        ! Local variables
        INTEGER(HSIZE_T) :: shape(1)
        INTEGER(HID_T) :: dtype
        TYPE(C_PTR) :: cptr

        IF (.NOT. ioProc) THEN
            RETURN
        END IF

        SELECT TYPE (list)
        TYPE IS (INTEGER(kind=intk))
            cptr = C_LOC(list)
            shape(1) = SIZE(list)
            dtype = mglet_hdf5_int
        TYPE IS (REAL(kind=realk))
            cptr = C_LOC(list)
            shape(1) = SIZE(list)
            dtype = mglet_hdf5_real
        CLASS DEFAULT
            CALL errr(__FILE__, __LINE__)
        END SELECT
        CALL stencilio_read_master_cptr(group_id, name, cptr, shape, dtype)
    END SUBROUTINE stencilio_read_master


    ! The IO master rank read a list. No other processes get any data.
    SUBROUTINE stencilio_read_master_cptr(group_id, name, cptr, shape, dtype)
        ! Subroutine arguments
        INTEGER(HID_T), INTENT(IN) :: group_id
        CHARACTER(LEN=*), INTENT(IN) :: name
        TYPE(C_PTR), INTENT(INOUT) :: cptr
        INTEGER(HSIZE_T), INTENT(IN) :: shape(:)
        INTEGER(HID_T), INTENT(IN) :: dtype

        ! Local variables
        INTEGER(intk) :: i, rank
        INTEGER(HSIZE_T) :: count(maxrank), offset(maxrank), dsetshape(maxrank)
        INTEGER(HID_T) :: dset_id, filespace, memspace, plist_id
        INTEGER(kind=int32) :: ierr

        IF (.NOT. ioProc) THEN
            RETURN
        END IF

        rank = SIZE(shape)
        IF (rank > maxrank .AND. ioid == 0) THEN
            WRITE(*,*) "rank, maxrank", rank, maxrank
            CALL errr(__FILE__, __LINE__)
        END IF

        ! Open dataset and filespace
        CALL hdf5common_dataset_open(name, dsetshape(1:rank), group_id, &
            dset_id, filespace)

        DO i = 1, rank
            IF (dsetshape(i) /= shape(i) .AND. ioid == 0) THEN
                WRITE(*,*) "i, dsetshape(i), shape(i)", i, &
                    dsetshape(i), shape(i)
                CALL errr(__FILE__, __LINE__)
            END IF
        END DO

        ! Create memspace
        CALL h5screate_simple_f(rank, shape, memspace, ierr)
        IF (ierr /= 0) CALL errr(__FILE__, __LINE__)

        ! Only ioid = 0 read data
        offset = 0
        count = 0
        count(1:rank) = shape(1:rank)
        IF (ioid == 0) THEN
            CALL h5sselect_hyperslab_f(memspace, H5S_SELECT_SET_F, &
                offset, count, ierr)
            IF (ierr /= 0) CALL errr(__FILE__, __LINE__)

            CALL h5sselect_hyperslab_f(filespace, H5S_SELECT_SET_F, &
                offset, count, ierr)
            IF (ierr /= 0) CALL errr(__FILE__, __LINE__)
        ELSE
            CALL h5sselect_none_f(memspace, ierr)
            IF (ierr /= 0) CALL errr(__FILE__, __LINE__)

            CALL h5sselect_none_f(filespace, ierr)
            IF (ierr /= 0) CALL errr(__FILE__, __LINE__)
        END IF

        ! Property list for collective dataset read
        CALL h5pcreate_f(H5P_DATASET_XFER_F, plist_id, ierr)
        IF (ierr /= 0) CALL errr(__FILE__, __LINE__)

        ! Only the IO master rank 0 reads any data. It makes sense to
        ! always use independent IO here, and not use the hdf5_io_mode
        ! switch
        CALL h5pset_dxpl_mpio_f(plist_id, H5FD_MPIO_INDEPENDENT_F, ierr)
        IF (ierr /= 0) CALL errr(__FILE__, __LINE__)

        ! Read dataset collectively
        CALL h5dread_f(dset_id, dtype, cptr, ierr, &
            file_space_id=filespace, mem_space_id=memspace, &
            xfer_prp=plist_id)
        IF (ierr /= 0) CALL errr(__FILE__, __LINE__)

        ! Close handlers
        CALL h5pclose_f(plist_id, ierr)
        IF (ierr /= 0) CALL errr(__FILE__, __LINE__)

        CALL h5sclose_f(memspace, ierr)
        IF (ierr /= 0) CALL errr(__FILE__, __LINE__)

        CALL hdf5common_dataset_close(dset_id, filespace)
    END SUBROUTINE stencilio_read_master_cptr


    ! This routine will write from master IO rank an array. If the array is
    ! given directly it must be a 1-D array, higher dimensions are supported
    ! by giving the C-ptr, shape and dtype of the array explicitly.
    !
    ! All ranks must give the same datatype and shape, the pointer at other
    ! ranks are not used.
    !
    ! Maximum number of dims are 4
    SUBROUTINE stencilio_write_master(group_id, name, list)
        ! Subroutine arguments
        INTEGER(HID_T), INTENT(IN) :: group_id
        CHARACTER(LEN=*), INTENT(IN) :: name
        CLASS(*), INTENT(in), TARGET :: list(:)

        ! Local variables
        INTEGER(HSIZE_T) :: shape(1)
        INTEGER(HID_T) :: dtype
        TYPE(C_PTR) :: cptr

        ! Non-io processes can now return
        IF (.NOT. ioproc) THEN
            RETURN
        END IF

        SELECT TYPE (list)
        TYPE IS (INTEGER(kind=intk))
            cptr = C_LOC(list)
            shape(1) = SIZE(list)
            dtype = mglet_hdf5_int
        TYPE IS (REAL(kind=realk))
            cptr = C_LOC(list)
            shape(1) = SIZE(list)
            dtype = mglet_hdf5_real
        CLASS DEFAULT
            CALL errr(__FILE__, __LINE__)
        END SELECT
        CALL stencilio_write_master_cptr(group_id, name, cptr, shape, dtype)
    END SUBROUTINE stencilio_write_master


    ! This routine will write from master IO rank an arbitrarily shaped array
    ! passed in as C-ptr with shape and dtype
    SUBROUTINE stencilio_write_master_cptr(group_id, name, cptr, shape, dtype)
        ! Subroutine arguments
        INTEGER(HID_T), INTENT(IN) :: group_id
        CHARACTER(LEN=*), INTENT(IN) :: name
        TYPE(C_PTR), INTENT(in) :: cptr
        INTEGER(HSIZE_T), INTENT(IN) :: shape(:)
        INTEGER(HID_T), INTENT(IN) :: dtype

        ! Local variables
        INTEGER(intk) :: rank
        INTEGER(HSIZE_T) :: count(maxrank), offset(maxrank)
        INTEGER(HID_T) :: dset_id, filespace, memspace, plist_id
        INTEGER(int32) :: ierr

        ! Non-io processes can now return
        IF (.NOT. ioproc) THEN
            RETURN
        END IF

        rank = SIZE(shape)
        IF (rank > maxrank) THEN
            WRITE(*,*) "rank, maxrank", rank, maxrank
            CALL errr(__FILE__, __LINE__)
        END IF

        ! Bacst shape into count array such that all IO processes have the
        ! same dataset size
        count = 0
        count(1:rank) = shape
        CALL MPI_Bcast(count, maxrank, mglet_mpi_hsize_t%MPI_Val, 0, &
            iocomm%MPI_val, ierr)

        ! Open/create dataset and filespace
        ! Count is synchronized among all IO processes and thus being the same
        CALL hdf5common_dataset_create(name, count(1:rank), dtype, group_id, &
            dset_id, filespace)

        ! Create memspace
        ! Each individual process use it's own 'shape' (not count) - because
        ! the memspace might differ
        CALL h5screate_simple_f(rank, shape, memspace, ierr)
        IF (ierr /= 0) CALL errr(__FILE__, __LINE__)

        ! Only ioid = 0 do the actual write
        offset = 0
        IF (ioid == 0) THEN
            CALL h5sselect_hyperslab_f(memspace, H5S_SELECT_SET_F, &
                offset, count, ierr)
            IF (ierr /= 0) CALL errr(__FILE__, __LINE__)

            CALL h5sselect_hyperslab_f(filespace, H5S_SELECT_SET_F, &
                offset, count, ierr)
            IF (ierr /= 0) CALL errr(__FILE__, __LINE__)
        ELSE
            CALL h5sselect_none_f(memspace, ierr)
            IF (ierr /= 0) CALL errr(__FILE__, __LINE__)

            CALL h5sselect_none_f(filespace, ierr)
            IF (ierr /= 0) CALL errr(__FILE__, __LINE__)
        END IF

        ! Property list for collective dataset write
        CALL h5pcreate_f(H5P_DATASET_XFER_F, plist_id, ierr)
        IF (ierr /= 0) CALL errr(__FILE__, __LINE__)

        ! Only the IO master rank 0 reads any data. It makes sense to
        ! always use independent IO here, and not use the hdf5_io_mode
        ! switch
        CALL h5pset_dxpl_mpio_f(plist_id, H5FD_MPIO_INDEPENDENT_F, ierr)
        IF (ierr /= 0) CALL errr(__FILE__, __LINE__)

        ! Write dataset collectively
        CALL h5dwrite_f(dset_id, dtype, cptr, ierr, &
            file_space_id=filespace, mem_space_id=memspace, &
            xfer_prp=plist_id)
        IF (ierr /= 0) CALL errr(__FILE__, __LINE__)

        ! Close handles
        CALL h5pclose_f(plist_id, ierr)
        IF (ierr /= 0) CALL errr(__FILE__, __LINE__)

        CALL h5sclose_f(memspace, ierr)
        IF (ierr /= 0) CALL errr(__FILE__, __LINE__)

        CALL hdf5common_dataset_close(dset_id, filespace)
    END SUBROUTINE stencilio_write_master_cptr


    ! Appends a list to the end of a dataset, appending happens at last
    ! dimension
    SUBROUTINE stencilio_append_master_cptr(group_id, name, cptr, shape, dtype)

        ! Subroutine arguments
        INTEGER(hid_t), INTENT(IN) :: group_id
        CHARACTER(LEN=*), INTENT(IN) :: name
        TYPE(C_PTR), INTENT(in) :: cptr
        INTEGER(hsize_t), INTENT(IN) :: shape(:)
        INTEGER(hid_t), INTENT(IN) :: dtype

        ! Local variables
        INTEGER(hsize_t) :: existing_shape(SIZE(shape))
        INTEGER(hsize_t) :: offset(SIZE(shape))
        INTEGER(hsize_t) :: count(SIZE(shape))
        INTEGER(hsize_t) :: new_shape(SIZE(shape))
        INTEGER(hid_t) :: dset_id, filespace, memspace
        INTEGER(intk) :: rank
        INTEGER(int32) :: hdferr

        IF (.NOT. ioproc) RETURN

        rank = SIZE(shape)

        ! Open dataset, extend, make filespace selection
        CALL hdf5common_dataset_open(name, existing_shape, group_id, dset_id, &
            filespace)

        offset = 0
        offset(rank) = existing_shape(rank)
        count = shape
        new_shape = existing_shape
        new_shape(rank) = existing_shape(rank) + shape(rank)

        CALL hdf5common_dataset_extend(dset_id, filespace, new_shape, &
            offset, count)
        IF (myid /= 0) THEN
            CALL h5sselect_none_f(filespace, hdferr)
            IF (hdferr /= 0) CALL errr(__FILE__, __LINE__)
        END IF

        ! Memspace
        CALL h5screate_simple_f(rank, shape, memspace, hdferr)
        IF (hdferr /= 0) CALL errr(__FILE__, __LINE__)

        IF (myid == 0) THEN
            CALL h5sselect_all_f(memspace, hdferr)
            IF (hdferr /= 0) CALL errr(__FILE__, __LINE__)
        ELSE
            CALL h5sselect_none_f(memspace, hdferr)
            IF (hdferr /= 0) CALL errr(__FILE__, __LINE__)
        END IF

        CALL h5dwrite_f(dset_id, dtype, cptr, hdferr, &
            file_space_id=filespace, mem_space_id=memspace)
        IF (hdferr /= 0) CALL errr(__FILE__, __LINE__)

        CALL hdf5common_dataset_close(dset_id, filespace)

        CALL h5sclose_f(memspace, hdferr)
        IF (hdferr /= 0) CALL errr(__FILE__, __LINE__)
    END SUBROUTINE stencilio_append_master_cptr

END MODULE stencilio_mod
