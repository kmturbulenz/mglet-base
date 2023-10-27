MODULE stencils_mod
    USE HDF5
    USE MPI_f08
    USE core_mod
    IMPLICIT NONE(type, external)
    PRIVATE

    ! Common stencils: geometry, IB
    TYPE, ABSTRACT :: stencils_t
        CHARACTER(len=mglet_filename_max) :: file
        CHARACTER(len=mglet_filename_max), ALLOCATABLE :: stlnames(:)

        TYPE(real_stencils_t), ALLOCATABLE :: points(:)
        TYPE(int_stencils_t), ALLOCATABLE :: cells(:)
        TYPE(int_stencils_t), ALLOCATABLE :: cellind(:)
        INTEGER(intk), ALLOCATABLE :: ncells(:)
        TYPE(real_stencils_t), ALLOCATABLE :: area(:)

        TYPE(int_stencils_t), ALLOCATABLE :: bpind(:)

        TYPE(int_stencils_t), ALLOCATABLE :: icellind(:)
        TYPE(int_stencils_t), ALLOCATABLE :: bodyid(:)
        TYPE(real_stencils_t), ALLOCATABLE :: sxsysz(:)
        TYPE(real_stencils_t), ALLOCATABLE :: ucell(:)

        TYPE(int_stencils_t), ALLOCATABLE :: auind(:)
        TYPE(int_stencils_t), ALLOCATABLE :: avind(:)
        TYPE(int_stencils_t), ALLOCATABLE :: awind(:)

        TYPE(real_stencils_t), ALLOCATABLE :: auvalue(:)
        TYPE(real_stencils_t), ALLOCATABLE :: avvalue(:)
        TYPE(real_stencils_t), ALLOCATABLE :: awvalue(:)
    CONTAINS
        PROCEDURE :: finish
        PROCEDURE :: set_geometry
        PROCEDURE :: set_bodyid
        PROCEDURE :: set_intersected
        PROCEDURE :: set_stlnames
        PROCEDURE :: get_bp
        PROCEDURE :: get_icells
        PROCEDURE :: get_intersected
        PROCEDURE :: read_stencils
        PROCEDURE :: write_stencils
        PROCEDURE, PRIVATE :: write_stlnames
        PROCEDURE, PRIVATE :: read_stlnames
        PROCEDURE(write_i), DEFERRED :: read
        PROCEDURE(write_i), DEFERRED :: write
    END TYPE stencils_t

    ABSTRACT INTERFACE
        SUBROUTINE write_i(this)
            IMPORT :: stencils_t
            CLASS(stencils_t), INTENT(inout) :: this
        END SUBROUTINE write_i
    END INTERFACE

    PUBLIC :: stencils_t
CONTAINS
    SUBROUTINE finish(this)
        ! Subroutine arguments
        CLASS(stencils_t), INTENT(inout) :: this

        IF (ALLOCATED(this%points)) DEALLOCATE(this%points)
        IF (ALLOCATED(this%cells)) DEALLOCATE(this%cells)
        IF (ALLOCATED(this%cellind)) DEALLOCATE(this%cellind)
        IF (ALLOCATED(this%ncells)) DEALLOCATE(this%ncells)
        IF (ALLOCATED(this%area)) DEALLOCATE(this%area)
        IF (ALLOCATED(this%bodyid)) DEALLOCATE(this%bodyid)

        IF (ALLOCATED(this%bpind)) DEALLOCATE(this%bpind)

        IF (ALLOCATED(this%icellind)) DEALLOCATE(this%icellind)
        IF (ALLOCATED(this%bodyid)) DEALLOCATE(this%bodyid)
        IF (ALLOCATED(this%sxsysz)) DEALLOCATE(this%sxsysz)
        IF (ALLOCATED(this%ucell)) DEALLOCATE(this%ucell)

        IF (ALLOCATED(this%auind)) DEALLOCATE(this%auind)
        IF (ALLOCATED(this%avind)) DEALLOCATE(this%avind)
        IF (ALLOCATED(this%awind)) DEALLOCATE(this%awind)

        IF (ALLOCATED(this%auvalue)) DEALLOCATE(this%auvalue)
        IF (ALLOCATED(this%avvalue)) DEALLOCATE(this%avvalue)
        IF (ALLOCATED(this%awvalue)) DEALLOCATE(this%awvalue)
    END SUBROUTINE finish


    SUBROUTINE set_geometry(this, icells, icellsall, icelllist, ncon, &
            connect, xx, cellind_p, arealist, igrid)

        ! Subroutine arguments
        CLASS(stencils_t), INTENT(inout) :: this
        INTEGER(intk), INTENT(in) :: icells, icellsall
        INTEGER(intk), INTENT(in) :: icelllist(icells)
        INTEGER(intk), INTENT(in) :: ncon(icellsall*7), connect(6, icellsall*7)
        REAL(realk), INTENT(in) :: xx(3, icellsall*20)
        INTEGER(intk), INTENT(in) :: cellind_p(icellsall)
        REAL(realk), INTENT(in) :: arealist(icellsall)
        INTEGER(intk), INTENT(in) :: igrid

        ! Local variables
        INTEGER(intk) :: i, idx, ika, counter, ipoints
        INTEGER(intk) :: npointswrite, ncellsdata

        IF (.NOT. ALLOCATED(this%points)) ALLOCATE(this%points(nmygrids))
        IF (.NOT. ALLOCATED(this%cells)) ALLOCATE(this%cells(nmygrids))
        IF (.NOT. ALLOCATED(this%ncells)) THEN
           ALLOCATE(this%ncells(nmygrids))
           this%ncells = 0
        END IF
        IF (.NOT. ALLOCATED(this%cellind)) ALLOCATE(this%cellind(nmygrids))
        IF (.NOT. ALLOCATED(this%area)) ALLOCATE(this%area(nmygrids))

        npointswrite = 0
        DO idx = 1, icells
            npointswrite = npointswrite + ncon(icelllist(idx))
        END DO
        ncellsdata = npointswrite + 2*icells

        ! i is the index in the mygrids list
        CALL get_imygrid(i, igrid)

        this%ncells(i) = icells
        ALLOCATE(this%points(i)%arr(3*npointswrite))
        ALLOCATE(this%cells(i)%arr(ncellsdata))
        ALLOCATE(this%cellind(i)%arr(icells))
        ALLOCATE(this%area(i)%arr(icells))

        counter = 1
        DO idx = 1, icells
            DO ika = 1, ncon(icelllist(idx))
                this%points(i)%arr(counter) = &
                    xx(1, connect(ika, icelllist(idx)) + 1)
                this%points(i)%arr(counter+1) = &
                    xx(2, connect(ika, icelllist(idx)) + 1)
                this%points(i)%arr(counter+2) = &
                    xx(3, connect(ika, icelllist(idx)) + 1)
                counter = counter + 3
           END DO
        END DO

        counter = 1
        ipoints = 0
        DO idx = 1, icells
            this%cells(i)%arr(counter) = 3
            this%cells(i)%arr(counter + 1) = ncon(icelllist(idx))
            counter = counter + 2
            DO ika = 1, ncon(icelllist(idx))
                this%cells(i)%arr(counter) = ipoints + ika - 1
                counter = counter + 1
            END DO
            ipoints = ipoints + ncon(icelllist(idx))
        END DO

        counter = 0
        DO idx = 1, icells
            counter = counter + 1
            this%cellind(i)%arr(counter) = cellind_p(icelllist(idx))
        END DO

        counter = 0
        DO idx = 1, icells
            counter = counter + 1
            this%area(i)%arr(counter) = arealist(icelllist(idx))
        END DO
    END SUBROUTINE set_geometry


    SUBROUTINE set_bodyid(this, igrid, icellsall, bodyid_p, icells, icelllist)
        ! Subroutine arguments
        CLASS(stencils_t), INTENT(inout) :: this
        INTEGER(intk), INTENT(in) :: igrid
        INTEGER(intk), INTENT(in) :: icellsall
        INTEGER(intk), INTENT(in) :: bodyid_p(icellsall)
        INTEGER(intk), INTENT(in) :: icells
        INTEGER(intk), INTENT(in), OPTIONAL :: icelllist(icells)

        INTEGER(intk) :: i, idx, counter

        IF (.NOT. ALLOCATED(this%bodyid)) ALLOCATE(this%bodyid(nmygrids))

        ! i is the index in the mygrids list
        CALL get_imygrid(i, igrid)

        IF (PRESENT(icelllist)) THEN
            ALLOCATE(this%bodyid(i)%arr(icells))
            counter = 0
            DO idx = 1, icells
                counter = counter + 1
                this%bodyid(i)%arr(counter) = bodyid_p(icelllist(idx))
            END DO
        ELSE
            ALLOCATE(this%bodyid(i)%arr(icellsall))
            this%bodyid(i)%arr = bodyid_p
        END IF
    END SUBROUTINE set_bodyid


    SUBROUTINE set_intersected(this, icells, icellspointer, &
            bodyid, sxsysz, ucell, bzelltyp_p)
        ! Subroutine arguments
        CLASS(stencils_t), INTENT(inout) :: this
        INTEGER(intk), INTENT(in) :: icells(:)
        INTEGER(intk), INTENT(in) :: icellspointer(:)
        INTEGER(intk), INTENT(in) :: bodyid(:)
        REAL(realk), INTENT(in) :: sxsysz(:, :)
        REAL(realk), INTENT(in) :: ucell(:, :)
        INTEGER(intk), TARGET, INTENT(in) :: bzelltyp_p(idim3d)

        INTEGER(intk) :: imygrid, igrid, ncells, ipp, ind, ip3
        INTEGER(intk) :: kk, jj, ii, k, j, i
        INTEGER(intk), POINTER, CONTIGUOUS :: bzelltyp(:, :, :)

        IF (.NOT. ALLOCATED(this%icellind)) ALLOCATE(this%icellind(nmygrids))
        IF (.NOT. ALLOCATED(this%bodyid)) ALLOCATE(this%bodyid(nmygrids))
        IF (.NOT. ALLOCATED(this%sxsysz)) ALLOCATE(this%sxsysz(nmygrids))
        IF (.NOT. ALLOCATED(this%ucell)) ALLOCATE(this%ucell(nmygrids))

        all_grids: DO imygrid = 1, nmygrids
            igrid = mygrids(imygrid)

            ncells = icells(igrid)
            ipp = icellspointer(igrid)

            ALLOCATE(this%bodyid(imygrid)%arr(ncells))
            ALLOCATE(this%sxsysz(imygrid)%arr(3*ncells))
            ALLOCATE(this%ucell(imygrid)%arr(3*ncells))

            this%bodyid(imygrid)%arr = bodyid(ipp:ipp+ncells-1)
            this%sxsysz(imygrid)%arr = RESHAPE(sxsysz(:, ipp:ipp+ncells-1), &
                [3*ncells])
            this%ucell(imygrid)%arr = RESHAPE(ucell(:, ipp:ipp+ncells-1), &
                [3*ncells])

            CALL get_mgdims(kk, jj, ii, igrid)
            CALL get_ip3(ip3, igrid)
            bzelltyp(1:kk, 1:jj, 1:ii) => bzelltyp_p(ip3:ip3+kk*jj*ii-1)

            ALLOCATE(this%icellind(imygrid)%arr(ncells))
            ! Counter used for sanity check
            ! Loop indices must match calcnormals_mod.F90, main loop in
            ! calcnormals_grid
            ncells = 0
            DO i = 2, ii
                DO j = 2, jj
                    DO k = 2, kk
                        IF (bzelltyp(k, j, i) >= 0) CYCLE
                        ncells = ncells + 1
                        CALL sub2ind(ind, k, j, i, kk, jj, ii)
                        this%icellind(imygrid)%arr(ncells) = ind
                    END DO
                END DO
            END DO
            ! Sanity check on indices
            IF (ncells /= icells(igrid)) THEN
                WRITE(*, '("igrid: ", I0)') igrid
                WRITE(*, '("ncells: ", I0)') ncells
                WRITE(*, '("icells: ", I0)') icells(igrid)
                CALL errr(__FILE__, __LINE__)
            END IF
        END DO all_grids
    END SUBROUTINE set_intersected


    SUBROUTINE set_stlnames(this, stlnames)
        ! Subroutine arguments
        CLASS(stencils_t), INTENT(inout) :: this
        CHARACTER(len=mglet_filename_max), INTENT(in) :: stlnames(:)

        INTEGER(intk) :: i, nstl, idx

        nstl = SIZE(stlnames)
        IF (.NOT. ALLOCATED(this%stlnames)) ALLOCATE(this%stlnames(nstl))

        DO i = 1, nstl
            ! Gets the last occurrence of "/" when present in filename -
            ! this is to strip off paths in the stored stlnames
            idx = INDEX(stlnames(i), "/", .TRUE., intk)
            this%stlnames(i) = stlnames(i)(idx+1:)
        END DO
    END SUBROUTINE set_stlnames


    SUBROUTINE get_bp(this, bp_f)
        ! Subroutine arguments
        CLASS(stencils_t), INTENT(inout) :: this
        TYPE(field_t), INTENT(inout) :: bp_f

        ! Local variables
        INTEGER(intk) :: imygrid, idx, igrid, ind, bpcells, ip3
        INTEGER(intk) :: kk, jj, ii, k, j, i
        REAL(realk), POINTER, CONTIGUOUS :: bp(:, :, :)

        ! Set BP to 1.0 initially
        bp_f%arr = 1.0

        ! Set all marked cells in the body to 0.0
        all_grids: DO imygrid = 1, nmygrids
            igrid = mygrids(imygrid)
            CALL get_mgdims(kk, jj, ii, igrid)
            CALL get_ip3(ip3, igrid)

            CALL bp_f%get_ptr(bp, igrid)
            bpcells = SIZE(this%bpind(imygrid)%arr)
            DO idx = 1, bpcells
                ind = this%bpind(imygrid)%arr(idx)
                CALL ind2sub(ind, k, j, i, kk, jj, ii)
                bp(k, j, i) = 0.0
            END DO
        END DO all_grids
    END SUBROUTINE get_bp


    SUBROUTINE get_icells(this, icells)
        ! Subroutine arguments
        CLASS(stencils_t), INTENT(inout) :: this
        INTEGER(intk), INTENT(out) :: icells(:)

        ! Local variables
        INTEGER(intk) :: imygrid, igrid

        icells = 0
        all_grids: DO imygrid = 1, nmygrids
            igrid = mygrids(imygrid)
            icells(igrid) = SIZE(this%icellind(imygrid)%arr)
        END DO all_grids
    END SUBROUTINE get_icells


    SUBROUTINE get_intersected(this, icells, icellspointer, &
            bodyid, sxsysz, ucell, bzelltyp_p)
        ! Subroutine arguments
        CLASS(stencils_t), INTENT(inout) :: this
        INTEGER(intk), INTENT(in) :: icells(:)
        INTEGER(intk), INTENT(in) :: icellspointer(:)
        INTEGER(intk), INTENT(out) :: bodyid(:)
        REAL(realk), INTENT(out) :: sxsysz(:, :)
        REAL(realk), INTENT(out) :: ucell(:, :)
        INTEGER(intk), TARGET, INTENT(out) :: bzelltyp_p(idim3d)

        ! Local variables
        INTEGER(intk) :: imygrid, igrid, ncells, ipp, ip3, idx, ind
        INTEGER(intk) :: kk, jj, ii, k, j, i
        INTEGER(intk), POINTER, CONTIGUOUS :: bzelltyp(:, :, :)

        bzelltyp_p = 0

        all_grids: DO imygrid = 1, nmygrids
            igrid = mygrids(imygrid)
            CALL get_mgdims(kk, jj, ii, igrid)
            CALL get_ip3(ip3, igrid)

            ncells = icells(igrid)
            ipp = icellspointer(igrid)

            ncells = SIZE(this%icellind(imygrid)%arr)
            IF (ncells /= icells(igrid)) CALL errr(__FILE__, __LINE__)

            bodyid(ipp:ipp+ncells-1) = this%bodyid(imygrid)%arr
            sxsysz(:, ipp:ipp+ncells-1) = RESHAPE(this%sxsysz(imygrid)%arr, &
                [3, ncells])
            ucell(:, ipp:ipp+ncells-1) = RESHAPE(this%ucell(imygrid)%arr, &
                [3, ncells])

            bzelltyp(1:kk, 1:jj, 1:ii) => bzelltyp_p(ip3:ip3+kk*jj*ii-1)
            DO idx = 1, ncells
                ind = this%icellind(imygrid)%arr(idx)
                CALL ind2sub(ind, k, j, i, kk, jj, ii)
                bzelltyp(k, j, i) = -idx
            END DO
        END DO all_grids
    END SUBROUTINE get_intersected


    SUBROUTINE read_stencils(this, file_id)
        ! Subroutine arguments
        CLASS(stencils_t), INTENT(inout) :: this
        INTEGER(hid_t), INTENT(in) :: file_id

        ALLOCATE(this%icellind(nmygrids))
        ALLOCATE(this%bodyid(nmygrids))
        ALLOCATE(this%sxsysz(nmygrids))
        ALLOCATE(this%ucell(nmygrids))

        ALLOCATE(this%auind(nmygrids))
        ALLOCATE(this%avind(nmygrids))
        ALLOCATE(this%awind(nmygrids))

        ALLOCATE(this%auvalue(nmygrids))
        ALLOCATE(this%avvalue(nmygrids))
        ALLOCATE(this%awvalue(nmygrids))

        ALLOCATE(this%bpind(nmygrids))

        ! Read stlnames
        CALL this%read_stlnames(file_id)

        CALL stencilio_read(file_id, 'icellind', this%icellind)
        CALL stencilio_read(file_id, 'bodyid', this%bodyid)
        CALL stencilio_read(file_id, 'sxsysz', this%sxsysz)
        CALL stencilio_read(file_id, 'ucell', this%ucell)

        CALL stencilio_read(file_id, 'auind', this%auind)
        CALL stencilio_read(file_id, 'avind', this%avind)
        CALL stencilio_read(file_id, 'awind', this%awind)

        CALL stencilio_read(file_id, 'auvalue', this%auvalue)
        CALL stencilio_read(file_id, 'avvalue', this%avvalue)
        CALL stencilio_read(file_id, 'awvalue', this%awvalue)

        CALL stencilio_read(file_id, 'bpind', this%bpind)
    END SUBROUTINE read_stencils


    SUBROUTINE write_stencils(this, file_id)
        ! Subroutine arguments
        CLASS(stencils_t), INTENT(inout) :: this
        INTEGER(hid_t), INTENT(in) :: file_id

        CALL stencilio_write(file_id, 'icellind', this%icellind)
        CALL stencilio_write(file_id, 'bodyid', this%bodyid)
        CALL stencilio_write(file_id, 'sxsysz', this%sxsysz)
        CALL stencilio_write(file_id, 'ucell', this%ucell)

        CALL stencilio_write(file_id, 'auind', this%auind)
        CALL stencilio_write(file_id, 'avind', this%avind)
        CALL stencilio_write(file_id, 'awind', this%awind)

        CALL stencilio_write(file_id, 'auvalue', this%auvalue)
        CALL stencilio_write(file_id, 'avvalue', this%avvalue)
        CALL stencilio_write(file_id, 'awvalue', this%awvalue)

        CALL stencilio_write(file_id, 'bpind', this%bpind)

        ! Write out STL names
        CALL this%write_stlnames(file_id)

        IF (ALLOCATED(this%points)) THEN
            CALL stencilio_write(file_id, 'points', this%points)
            CALL stencilio_write(file_id, 'cells', this%cells)
            CALL stencilio_write(file_id, 'cellind', this%cellind)
            CALL stencilio_write_list(file_id, 'nCells', this%ncells)
            CALL stencilio_write(file_id, 'area', this%area)
        END IF
    END SUBROUTINE write_stencils


    SUBROUTINE write_stlnames(this, file_id)
        USE, INTRINSIC :: ISO_C_BINDING, ONLY: C_LOC

        ! Subroutine arguments
        CLASS(stencils_t), INTENT(inout), TARGET :: this
        INTEGER(hid_t), INTENT(in) :: file_id

        ! Local variables
        INTEGER(hid_t) :: str_t
        INTEGER(hsize_t) :: shape(1)
        TYPE(C_PTR) :: cptr
        INTEGER(int32) :: ierr

        shape(1) = SIZE(this%stlnames, dim=1, kind=hsize_t)
        IF (ioproc) THEN
            ! Create type for character array
            CALL h5tcopy_f(H5T_NATIVE_CHARACTER, str_t, ierr)
            IF (ierr /= 0) CALL errr(__FILE__, __LINE__)
            CALL h5tset_size_f(str_t, LEN(this%stlnames, kind=size_t), ierr)
            IF (ierr /= 0) CALL errr(__FILE__, __LINE__)

            cptr = C_LOC(this%stlnames)
            CALL stencilio_write_master_cptr(file_id, 'stlnames', &
                                             cptr, shape, str_t)

            CALL h5tclose_f(str_t, ierr)
            IF (ierr /= 0) CALL errr(__FILE__, __LINE__)
        END IF
    END SUBROUTINE write_stlnames


    SUBROUTINE read_stlnames(this, file_id)
        USE, INTRINSIC :: ISO_C_BINDING, ONLY: C_LOC

        ! Subroutine arguments
        CLASS(stencils_t), INTENT(inout), TARGET :: this
        INTEGER(hid_t), INTENT(in) :: file_id

        ! Local variables
        INTEGER(hid_t) :: str_t, dset_id, filespace
        INTEGER(hsize_t) :: shape(1)
        TYPE(C_PTR) :: cptr
        INTEGER(int32) :: nchar, ierr

        ! Open dataset to get number of STL's to read
        CALL hdf5common_dataset_open('stlnames', shape, file_id, &
            dset_id, filespace)
        CALL hdf5common_dataset_close(dset_id, filespace)

        ! Broadcast number of STL's and allocate memory on all ranks
        CALL MPI_Bcast(shape, 1, mglet_mpi_hsize_t, 0, MPI_COMM_WORLD)

        IF (ALLOCATED(this%stlnames)) THEN
            IF (SIZE(this%stlnames) /= shape(1)) THEN
                CALL errr(__FILE__, __LINE__)
            END IF
        ELSE
            ALLOCATE(this%stlnames(shape(1)))
        END IF

        ! Read actual data
        IF (ioproc .AND. shape(1) > 0) THEN
            ! Create type for character array
            CALL h5tcopy_f(H5T_NATIVE_CHARACTER, str_t, ierr)
            IF (ierr /= 0) CALL errr(__FILE__, __LINE__)
            CALL h5tset_size_f(str_t, LEN(this%stlnames, kind=size_t), ierr)
            IF (ierr /= 0) CALL errr(__FILE__, __LINE__)

            cptr = C_LOC(this%stlnames)
            CALL stencilio_read_master_cptr(file_id, 'stlnames', &
                                            cptr, shape, str_t)

            CALL h5tclose_f(str_t, ierr)
            IF (ierr /= 0) CALL errr(__FILE__, __LINE__)
        END IF

        ! Broadcast STL names to all processes
        nchar = INT(shape(1), int32)*LEN(this%stlnames)
        CALL MPI_Bcast(this%stlnames, nchar, MPI_CHARACTER, 0, &
            MPI_COMM_WORLD)
    END SUBROUTINE read_stlnames
END MODULE stencils_mod
