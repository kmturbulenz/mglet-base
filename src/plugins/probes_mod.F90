MODULE probes_mod
    USE core_mod
    USE MPI_f08
    USE HDF5

    IMPLICIT NONE(type, external)
    PRIVATE

    TYPE :: probearr_t
        CHARACTER(len=mglet_filename_max) :: name
        INTEGER(intk) :: nvars
        INTEGER(intk) :: npts
        INTEGER(intk) :: nmypnts
        INTEGER(intk) :: grouppnts
        INTEGER(intk) :: notfound
        CHARACTER(len=nchar_name), ALLOCATABLE :: variables(:)
        REAL(realk), ALLOCATABLE :: coordinates(:, :)
        INTEGER(intk), ALLOCATABLE :: id(:)
        INTEGER(intk), ALLOCATABLE :: groupid(:)
        INTEGER(intk), ALLOCATABLE :: grid(:)
        REAL(realk), ALLOCATABLE :: buffer(:, :, :)

        ! Data transfer types and helpers
        TYPE(MPI_Datatype), ALLOCATABLE :: mpitype(:)
    CONTAINS
        FINAL :: probearr_destructor
    END TYPE probearr_t

    LOGICAL :: has_probes = .FALSE.
    LOGICAL :: isinit_buffers = .FALSE.

    INTEGER(intk) :: itsamp
    INTEGER(intk) :: navg
    INTEGER(hsize_t) :: bufloc
    INTEGER(hsize_t) :: buflen
    INTEGER(hsize_t) :: fileoffset
    REAL(realk) :: tstart

    ! For initializing HDF5 file and buffers
    INTEGER(intk) :: mtstep_probes
    INTEGER(intk) :: itint_probes
    REAL(realk) :: tend_probes

    ! Storage of the individual probe arrays
    INTEGER(intk) :: narrays
    TYPE(probearr_t), ALLOCATABLE, TARGET :: arr(:)

    ! Storage of sampled time
    TYPE(timeinfo_t), ALLOCATABLE, TARGET :: timeinfo(:)

    ! Filename of probes output file
    CHARACTER(len=mglet_filename_max) :: outfile

    ! Chunksize of 1D datasets (time)
    INTEGER(hsize_t), PARAMETER :: chunksize11 = 16384

    ! Chunksize of 2D datasets in 1st and 2nd dimension
    INTEGER(hsize_t), PARAMETER :: chunksize21 = 16
    INTEGER(hsize_t), PARAMETER :: chunksize22 = 4096

    ! Max. buffered size per IO group = 1 GB
    INTEGER(int64), PARAMETER :: maxbytes = 1073741824

    PUBLIC :: init_probes, sample_probes, finish_probes, &
        can_checkpoint_probes, checkpoint_probes

CONTAINS
    SUBROUTINE init_probes(ittot, mtstep, itint, timeph, dt, tend)
        ! Subroutine arguments
        INTEGER(intk), INTENT(in) :: ittot
        INTEGER(intk), INTENT(in) :: mtstep
        INTEGER(intk), INTENT(in) :: itint
        REAL(realk), INTENT(in) :: timeph
        REAL(realk), INTENT(in) :: dt
        REAL(realk), INTENT(in) :: tend

        ! Local variables
        TYPE(config_t) :: probesconf, array
        CHARACTER(len=64) :: jsonptr
        INTEGER(intk) :: i, j, ivar
        TYPE(field_t), POINTER :: field

        has_probes = .FALSE.
        IF (.NOT. fort7%exists("/probes")) THEN
            RETURN
        END IF
        has_probes = .TRUE.

        CALL fort7%get(probesconf, "/probes")

        ! Arrays - if no, then there are no probes
        CALL probesconf%get_size("/arrays", narrays)
        IF (narrays < 1) THEN
            has_probes = .FALSE.
            CALL probesconf%finish()
            RETURN
        END IF

        CALL probesconf%get_value("/itsamp", itsamp)
        IF (itsamp <= 0) THEN
            WRITE(*, '("probes: itsamp must be positive: ", I7)') itsamp
            CALL errr(__FILE__, __LINE__)
        END IF

        CALL probesconf%get_value("/tstart", tstart, 0.0)
        IF (tstart < 0.0) THEN
            WRITE(*, '("probes: tstart must be >= 0.0: ", F15.7)') tstart
            CALL errr(__FILE__, __LINE__)
        END IF

        CALL probesconf%get_value("/file", outfile, "probes.h5")

        IF (myid == 0) THEN
            WRITE(*, '("READING PROBES:")')
            WRITE(*, '(4X, "itsamp: ", I15)') itsamp
            WRITE(*, '(4X, "tstart: ", F15.7)') tstart
        END IF

        ALLOCATE(arr(narrays))

        DO i = 1, narrays
            WRITE(jsonptr, '("/arrays/", I0)') i-1
            CALL probesconf%get(array, jsonptr)

            ! Array name
            CALL array%get_value("/name", arr(i)%name)

            ! Which variables to sample (U, V, P, etc)
            CALL array%get_size("/variables", arr(i)%nvars)
            ALLOCATE(arr(i)%variables(arr(i)%nvars))
            DO j = 1, arr(i)%nvars
                WRITE(jsonptr, '("/variables/", I0)') j-1
                CALL array%get_value(jsonptr, arr(i)%variables(j))
            END DO

            ! Probe positions from file or directly from JSON
            CALL read_positions(arr(i), array)
            CALL array%finish()
        END DO
        CALL probesconf%finish()

        ! Check that all variables exist
        DO i = 1, narrays
            DO ivar = 1, arr(i)%nvars
                ! This will fail if a requested field is not available
                ! Not the most beautiful error message, though...
                CALL get_field(field, arr(i)%variables(ivar))
            END DO
        END DO

        ! The HDF5 file and communication buffers are initialized at the first
        ! timestep where probes are sampled. The initialization process needs
        ! these values, but they are not passed over in the subroutine
        ! arguments. Therfore "keep" them in a module variable until later...
        mtstep_probes = mtstep
        itint_probes = itint
        tend_probes = tend

        ! Write more statistics to terminal
        IF (myid == 0) THEN
            DO i = 1, narrays
                WRITE(*, '(4X, A, ": nvars: ", I0, ", npts: ", I0, ' &
                    //'", notfound: ", I0)') TRIM(arr(i)%name), arr(i)%nvars, &
                    arr(i)%npts, arr(i)%notfound
            END DO
        END IF

        IF (myid == 0) WRITE(*, '()')
    END SUBROUTINE init_probes


    SUBROUTINE finish_probes()
        ! Subroutine arguments
        ! none...

        ! Local variables
        INTEGER(intk) :: iarr, i, lb, ub

        IF (.NOT. has_probes) RETURN
        IF (.NOT. isinit_buffers) RETURN

        ! Flush leftover data to disk
        CALL write_probes()

        ! Free MPI types
        DO iarr = 1, narrays
            IF (ALLOCATED(arr(iarr)%mpitype)) THEN
                lb = LBOUND(arr(iarr)%mpitype, dim=1)
                ub = UBOUND(arr(iarr)%mpitype, dim=1)
                DO i = lb, ub
                    CALL MPI_Type_free(arr(iarr)%mpitype(i))
                END DO
                DEALLOCATE(arr(iarr)%mpitype)
            END IF
        END DO

        ! Deallocate module-wide data
        DEALLOCATE(arr)
        narrays = 0

        IF (ALLOCATED(timeinfo)) DEALLOCATE(timeinfo)
    END SUBROUTINE finish_probes


    SUBROUTINE can_checkpoint_probes(can_checkpoint)
        ! Subroutine arguments
        LOGICAL, INTENT(out) :: can_checkpoint

        ! Local variables
        ! none...

        ! This is the default assumption
        can_checkpoint = .TRUE.

        ! When no probes are present or when nothing is initialized yet
        IF (.NOT. has_probes) RETURN
        IF (.NOT. isinit_buffers) RETURN

        ! Only when a sample is not finished a checkpoint cannot be made
        IF (navg > 0) THEN
            can_checkpoint = .FALSE.
        END IF
    END SUBROUTINE can_checkpoint_probes


    SUBROUTINE checkpoint_probes()
        ! Subroutine arguments
        ! none...

        ! Local variables
        ! none...

        IF (.NOT. has_probes) RETURN
        IF (.NOT. isinit_buffers) RETURN
        CALL write_probes()
    END SUBROUTINE checkpoint_probes


    SUBROUTINE read_positions(arr, arrayconf)
        ! Read and distribute positions among ranks

        ! Subroutine arguments
        TYPE(probearr_t), INTENT(inout) :: arr
        TYPE(config_t), INTENT(inout) :: arrayconf

        ! Local variables
        CHARACTER(len=mglet_filename_max) :: filename
        CHARACTER(len=64) :: jsonptr
        INTEGER(intk) :: i, igrid, ilevel, ipoint
        INTEGER(intk), ALLOCATABLE :: probesgrids(:)
        INTEGER(intk) :: kk, jj, ii
        INTEGER(intk) :: nfro, nbac, nrgt, nlft, nbot, ntop
        REAL(realk), ALLOCATABLE :: tmpcoords(:, :)
        REAL(realk) :: points(3)
        REAL(realk) :: minx, maxx, miny, maxy, minz, maxz
        REAL(realk) :: dx, dy, dz

        ! Rank 0 reads the coordinates into 'tmpcoords'
        IF (myid == 0) THEN
            IF (arrayconf%is_char("/positions")) THEN
                CALL arrayconf%get_value("/positions", filename)
                CALL read_datfile(tmpcoords, filename)
            ELSE
                CALL arrayconf%get_size("/positions", arr%npts)
                ALLOCATE(tmpcoords(3, arr%npts))
                DO i = 1, arr%npts
                    WRITE(jsonptr, '("/positions/", I0)') i-1
                    CALL arrayconf%get_array(jsonptr, tmpcoords(:, i))
                END DO
            END IF
            arr%npts = SIZE(tmpcoords, dim=2)
        END IF

        ! Broadcast the global number of points
        CALL MPI_Bcast(arr%npts, 1, mglet_mpi_int, 0, MPI_COMM_WORLD)

        ! Other ranks allocate memory and receive coordinates
        IF (myid /= 0) THEN
            ALLOCATE(tmpcoords(3, arr%npts))
        END IF
        CALL MPI_Bcast(tmpcoords, 3*arr%npts, mglet_mpi_real, 0, &
            MPI_COMM_WORLD)

        ! Array to store which grid the probe belongs to
        ALLOCATE(probesgrids(arr%npts))
        probesgrids = 0

        ! Loop over the different grid levels, from finest to coarsest
        ! and see if we can find points that belong to this grid
        DO ilevel = maxlevel, minlevel, -1
            DO i = 1, nmygridslvl(ilevel)
                igrid = mygridslvl(i, ilevel)

                ! Find grid "bounding box"
                CALL get_bbox(minx, maxx, miny, maxy, minz, maxz, igrid)
                CALL get_mgdims(kk, jj, ii, igrid)
                dx = (maxx-minx)/REAL(ii-4, kind=realk)
                dy = (maxy-miny)/REAL(jj-4, kind=realk)
                dz = (maxz-minz)/REAL(kk-4, kind=realk)

                ! In case of a probe in vicinity of a PAR, interpolation
                ! is difficult becuase the PAR buffers are not suited for
                ! interpolation. Account for this by artificially shrink
                ! the gridbox and let the probe belong to the coarse
                ! grid instead
                CALL get_mgbasb(nfro, nbac, nrgt, nlft, nbot, ntop, igrid)
                IF (nfro == 8) minx = minx + dx/2.0
                IF (nbac == 8) maxx = maxx - dx/2.0
                IF (nrgt == 8) miny = miny + dy/2.0
                IF (nlft == 8) maxy = maxy - dy/2.0
                IF (nbot == 8) minz = minz + dz/2.0
                IF (ntop == 8) maxz = maxz - dz/2.0

                ! Loop over all points and see if we can find it on this grid
                DO ipoint = 1, arr%npts
                    ! Point is already found in another grid
                    IF (probesgrids(ipoint) > 0) CYCLE

                    points(:) = tmpcoords(:, ipoint)

                    ! If the current point is outside of the grid, skip the
                    ! rest of the loop
                    IF (points(1) <= minx .OR. points(1) > maxx .OR. &
                        points(2) <= miny .OR. points(2) > maxy .OR. &
                        points(3) <= minz .OR. points(3) > maxz) CYCLE

                    ! Point is found on this grid
                    probesgrids(ipoint) = igrid
                END DO
            END DO
        END DO

        ! Do an allreduce with MPI_MAX to find the finest grid level with this
        ! grid within
        CALL MPI_Allreduce(MPI_IN_PLACE, probesgrids, arr%npts, &
            mglet_mpi_int, MPI_MAX, MPI_COMM_WORLD)

        ! Number of points this process owns
        arr%nmypnts = 0
        DO i = 1, arr%npts
            igrid = probesgrids(i)
            IF (igrid < 1) CYCLE  ! For probes not assigned to any grid
            IF (myid == idprocofgrd(igrid)) arr%nmypnts = arr%nmypnts + 1
        END DO

        ! Allocate and store coordinates for points this process owns
        ALLOCATE(arr%coordinates(3, arr%nmypnts))
        ALLOCATE(arr%id(arr%nmypnts))
        ALLOCATE(arr%grid(arr%nmypnts))

        ipoint = 0
        DO i = 1, arr%npts
            igrid = probesgrids(i)
            IF (igrid < 1) CYCLE  ! For probes not assigned to any grid
            IF (myid == idprocofgrd(igrid)) THEN
                ipoint = ipoint + 1
                arr%coordinates(:, ipoint) = tmpcoords(:, i)
                arr%id(ipoint) = i
                arr%grid(ipoint) = probesgrids(i)
            END IF
        END DO

        ! Sum all not found points, that is points no process claims
        ! ownership to...
        arr%notfound = COUNT(probesgrids == 0)

        DEALLOCATE(probesgrids)
        DEALLOCATE(tmpcoords)
    END SUBROUTINE read_positions


    PURE ELEMENTAL SUBROUTINE probearr_destructor(this)
        ! Subroutine arguments
        TYPE(probearr_t), INTENT(inout) :: this

        ! Local variables
        ! none...

        IF (ALLOCATED(this%variables)) DEALLOCATE(this%variables)
        IF (ALLOCATED(this%coordinates)) DEALLOCATE(this%coordinates)
        IF (ALLOCATED(this%id)) DEALLOCATE(this%id)
        IF (ALLOCATED(this%groupid)) DEALLOCATE(this%groupid)
        IF (ALLOCATED(this%grid)) DEALLOCATE(this%grid)
        IF (ALLOCATED(this%buffer)) DEALLOCATE(this%buffer)
        IF (ALLOCATED(this%mpitype)) DEALLOCATE(this%mpitype)
    END SUBROUTINE probearr_destructor


    SUBROUTINE sample_probes(itstep, ittot, timeph, dt)
        ! Subroutine arguments
        INTEGER(intk), INTENT(in) :: itstep
        INTEGER(intk), INTENT(in) :: ittot
        REAL(realk), INTENT(in) :: timeph
        REAL(realk), INTENT(in) :: dt

        ! Local variables
        INTEGER(intk) :: iarr, ivar, i, igrid
        TYPE(field_t), POINTER :: field

        IF (.NOT. has_probes) RETURN

        ! Start sampling at tstart
        IF (timeph < tstart) RETURN

        ! File and buffers are initialized at first sampled timestep
        IF (.NOT. isinit_buffers) THEN
            CALL init_file_and_buffers(ittot, mtstep_probes, itint_probes, &
                timeph, dt, tend_probes)
        END IF

        ! Do actual sampling
        ! bufloc is set to 1 when buffers are allocated first time
        DO iarr = 1, narrays
            DO ivar = 1, arr(iarr)%nvars
                CALL get_field(field, arr(iarr)%variables(ivar))
                DO i = 1, nmygrids
                    igrid = mygrids(i)
                    CALL sample_variable(arr(iarr)%buffer(:, bufloc, ivar), &
                        field, igrid, iarr)
                END DO
            END DO
        END DO

        ! Sample time
        IF (ioid == 0) THEN
            timeinfo(bufloc)%ittot = ittot
            timeinfo(bufloc)%time = timeph
        END IF

        ! Increment time counter for running averages
        navg = navg + 1
        IF (navg >= itsamp) THEN
            bufloc = bufloc + 1
            navg = 0
        END IF

        ! When buffer(s) are full - write to disk
        IF (bufloc > buflen) THEN
            CALL write_probes()
        END IF
    END SUBROUTINE sample_probes


    SUBROUTINE sample_variable(buffer, field, igrid, iarr)
        ! Subroutine arguments
        REAL(realk), INTENT(inout), CONTIGUOUS :: buffer(:)
        TYPE(field_t), INTENT(in) :: field
        INTEGER(intk), INTENT(in) :: igrid
        INTEGER(intk), INTENT(in) :: iarr

        ! Local variables
        REAL(realk), CONTIGUOUS, POINTER :: field3d(:, :, :), bp(:, :, :)
        REAL(realk), CONTIGUOUS, POINTER :: x(:), y(:), z(:)
        INTEGER(intk) :: iprobe
        REAL(realk) :: samplevar, fac1, fac2

        IF (field%istag == 0) THEN
            CALL get_fieldptr(x, "X", igrid)
        ELSE
            CALL get_fieldptr(x, "XSTAG", igrid)
        END IF

        IF (field%jstag == 0) THEN
            CALL get_fieldptr(y, "Y", igrid)
        ELSE
            CALL get_fieldptr(y, "YSTAG", igrid)
        END IF

        IF (field%kstag == 0) THEN
            CALL get_fieldptr(z, "Z", igrid)
        ELSE
            CALL get_fieldptr(z, "ZSTAG", igrid)
        END IF

        CALL field%get_ptr(field3d, igrid)

        CALL get_fieldptr(bp, "BP", igrid)

        ! Factors for running averages
        fac1 = REAL(navg)/REAL(navg + 1)
        fac2 = 1.0/REAL(navg + 1)

        DO iprobe = 1, arr(iarr)%nmypnts
            IF (arr(iarr)%grid(iprobe) /= igrid) CYCLE

            IF (field%istag == 0 .AND. field%jstag == 0 .AND. &
                    field%kstag == 0) THEN
                CALL solintxyzgrad0(samplevar, &
                    arr(iarr)%coordinates(:, iprobe), field3d, x, y, z, bp)
            ELSE
                CALL solintxyz(samplevar, &
                    arr(iarr)%coordinates(:, iprobe), field3d, x, y, z)
            END IF

            buffer(iprobe) = fac1*buffer(iprobe) + fac2*samplevar
        END DO
    END SUBROUTINE sample_variable


    PURE SUBROUTINE solintxyz(sol, coord, phi, x, y, z)
        ! Subroutine arguments
        REAL(realk), INTENT(out) :: sol
        REAL(realk), INTENT(in) :: coord(3)
        REAL(realk), INTENT(in), CONTIGUOUS :: phi(:, :, :)
        REAL(realk), INTENT(in), CONTIGUOUS :: x(:)
        REAL(realk), INTENT(in), CONTIGUOUS :: y(:)
        REAL(realk), INTENT(in), CONTIGUOUS :: z(:)

        ! Local variables
        INTEGER(intk) :: kp, jp, ip
        INTEGER(intk) :: i
        REAL(realk) :: facx, facy, facz
        REAL(realk) :: psi(8), field(8)

        CALL get_idx(kp, coord(3), z)
        CALL get_idx(jp, coord(2), y)
        CALL get_idx(ip, coord(1), x)

        field(1) = phi(kp, jp, ip)
        field(2) = phi(kp, jp, ip+1)
        field(3) = phi(kp, jp+1, ip+1)
        field(4) = phi(kp, jp+1, ip)
        field(5) = phi(kp+1, jp, ip)
        field(6) = phi(kp+1, jp, ip+1)
        field(7) = phi(kp+1, jp+1, ip+1)
        field(8) = phi(kp+1, jp+1, ip)

        facx = (coord(1) - x(ip))/(x(ip+1) - x(ip))
        facy = (coord(2) - y(jp))/(y(jp+1) - y(jp))
        facz = (coord(3) - z(kp))/(z(kp+1) - z(kp))

        psi(1) = (1.0-facx)*(1.0-facy)*(1.0-facz)
        psi(2) = facx*(1.0-facy)*(1.0-facz)
        psi(3) = facx*facy*(1.0-facz)
        psi(4) = (1.0-facx)*facy*(1.0-facz)
        psi(5) = (1.0-facx)*(1.0-facy)*facz
        psi(6) = facx*(1.0-facy)*facz
        psi(7) = facx*facy*facz
        psi(8) = (1.0-facx)*facy*facz

        sol = 0.0
        DO i = 1, 8
            sol = sol + psi(i)*field(i)
        END DO
    END SUBROUTINE solintxyz


    PURE SUBROUTINE solintxyzgrad0(sol, coord, phi, x, y, z, bp)
        ! Subroutine arguments
        REAL(realk), INTENT(out) :: sol
        REAL(realk), INTENT(in) :: coord(3)
        REAL(realk), INTENT(in), CONTIGUOUS :: phi(:, :, :)
        REAL(realk), INTENT(in), CONTIGUOUS :: x(:)
        REAL(realk), INTENT(in), CONTIGUOUS :: y(:)
        REAL(realk), INTENT(in), CONTIGUOUS :: z(:)
        REAL(realk), INTENT(in), CONTIGUOUS :: bp(:, :, :)

        ! Local variables
        INTEGER(intk) :: kk, jj, ii
        INTEGER(intk) :: kp, jp, ip
        INTEGER(intk) :: k, j, i
        REAL(realk) :: facx, facy, facz, sumfac
        REAL(realk) :: psi(8), field(8), blocked(8)
        REAL(realk) :: lastdistance, distance

        CALL get_idx(kp, coord(3), z)
        CALL get_idx(jp, coord(2), y)
        CALL get_idx(ip, coord(1), x)

        field(1) = phi(kp, jp, ip)
        field(2) = phi(kp, jp, ip+1)
        field(3) = phi(kp, jp+1, ip+1)
        field(4) = phi(kp, jp+1, ip)
        field(5) = phi(kp+1, jp, ip)
        field(6) = phi(kp+1, jp, ip+1)
        field(7) = phi(kp+1, jp+1, ip+1)
        field(8) = phi(kp+1, jp+1, ip)

        blocked(1) = bp(kp, jp, ip)
        blocked(2) = bp(kp, jp, ip+1)
        blocked(3) = bp(kp, jp+1, ip+1)
        blocked(4) = bp(kp, jp+1, ip)
        blocked(5) = bp(kp+1, jp, ip)
        blocked(6) = bp(kp+1, jp, ip+1)
        blocked(7) = bp(kp+1, jp+1, ip+1)
        blocked(8) = bp(kp+1, jp+1, ip)

        facx = (coord(1) - x(ip))/(x(ip+1) - x(ip))
        facy = (coord(2) - y(jp))/(y(jp+1) - y(jp))
        facz = (coord(3) - z(kp))/(z(kp+1) - z(kp))

        psi(1) = (1.0-facx)*(1.0-facy)*(1.0-facz)
        psi(2) = facx*(1.0-facy)*(1.0-facz)
        psi(3) = facx*facy*(1.0-facz)
        psi(4) = (1.0-facx)*facy*(1.0-facz)
        psi(5) = (1.0-facx)*(1.0-facy)*facz
        psi(6) = facx*(1.0-facy)*facz
        psi(7) = facx*facy*facz
        psi(8) = (1.0-facx)*facy*facz

        sumfac = 0.0
        DO i = 1, 8
           sumfac = sumfac + psi(i)*blocked(i)
        END DO

        IF (sumfac > 0.02) THEN
            ! linear interpolation
            ! bp=0 cells are excluded
            sol = 0.0
            DO i = 1, 8
                sol = sol + psi(i)*field(i)
            END DO
        ELSE
            ! Pick nearest point with bp=1
            kk = SIZE(z)
            jj = SIZE(y)
            ii = SIZE(x)

            sol = 0.0
            lastdistance = (x(ii)-x(1))**2 + (y(jj)-y(1))**2 + (z(kk)-z(1))**2

            ! Probe is in between ip and ip+1, therefore the indices are
            ! extended one cell in each direction
            DO i = ip-1, ip+2
                DO j = jp-1, jp+2
                    DO k = kp-1, kp+2
                        IF (bp(k, j, i) < 0.5) CYCLE

                        distance = (x(i)-coord(1))**2 + &
                            (y(j)-coord(2))**2 + (z(k)-coord(3))**2

                        IF (distance < lastdistance) THEN
                            sol = phi(k, j, i)
                            lastdistance = distance
                        END IF
                    END DO
                END DO
            END DO
        END IF
    END SUBROUTINE solintxyzgrad0


    SUBROUTINE init_file_and_buffers(ittot, mtstep, itint, timeph, dt, tend)
        ! Subroutine arguments
        INTEGER(intk), INTENT(in) :: ittot
        INTEGER(intk), INTENT(in) :: mtstep
        INTEGER(intk), INTENT(in) :: itint
        REAL(realk), INTENT(in) :: timeph
        REAL(realk), INTENT(in) :: dt
        REAL(realk), INTENT(in) :: tend

        ! Local variables
        INTEGER(HID_T)  :: fileh
        INTEGER(int32) :: ierr
        INTEGER(intk) :: iarr
        LOGICAL :: link_exists

        IF (isinit_buffers) CALL errr(__FILE__, __LINE__)

        ! Initialize data buffers and transfer patterns - all processes do this
        CALL init_buffers(ittot, mtstep, itint, timeph, dt, tend)

        DO iarr = 1, narrays
            CALL init_data_transfer(arr(iarr))
        END DO

        IF (ioproc) THEN
            ! Open existing file if it extst, create new if not
            IF (dcont) THEN
                CALL hdf5common_open(outfile, 'a', fileh)
            ELSE
                CALL hdf5common_open(outfile, 'w', fileh)
            END IF

            ! Initialize file if required
            CALL h5lexists_f(fileh, "PROBES", link_exists, ierr)
            IF (ierr /= 0) CALL errr(__FILE__, __LINE__)

            IF (.NOT. link_exists) THEN
                CALL init_file(fileh)
            ELSE
                CALL read_file(fileh, ittot)
            END IF

            CALL hdf5common_close(fileh)
        END IF

        isinit_buffers = .TRUE.
        ! Write more statistics to terminal
        IF (myid == 0) THEN
            WRITE(*, '("INITIALIZING PROBES:")')
            WRITE(*, '(4X, "outfile: ", A)') TRIM(outfile)
            WRITE(*, '(4X, "buflen:  ", I15)') buflen
            WRITE(*, '(4X, "offset:  ", I15)') fileoffset
            WRITE(*, '()')
        END IF
    END SUBROUTINE init_file_and_buffers


    SUBROUTINE init_file(fileh)
        ! Subroutine arguments
        INTEGER(HID_T), INTENT(inout) :: fileh

        ! Local variables
        INTEGER(intk) :: iarr, ivar
        INTEGER(hid_t)  :: probes_grouph, arr_grouph, filespace, dset_id
        INTEGER(hsize_t) :: shape1(1), maxdims1(1), chunksize1(1)
        INTEGER(hsize_t) :: shape2(2), maxdims2(2), chunksize2(2)
        INTEGER(hsize_t) :: npts
        INTEGER(int32) :: hdferr

        ! Greate group(s)
        CALL h5gcreate_f(fileh, "PROBES", probes_grouph, hdferr)
        IF (hdferr /= 0) CALL errr(__FILE__, __LINE__)

        ! Initialize arrays
        DO iarr = 1, narrays
            ! Group for individual array
            CALL h5gcreate_f(probes_grouph, TRIM(arr(iarr)%name), arr_grouph, &
                hdferr)
            IF (hdferr /= 0) CALL errr(__FILE__, __LINE__)

            ! Variables legend (P, PC etc.)
            CALL hdf5common_attr_write_arr("VARIABLES", arr(iarr)%variables, &
                arr_grouph)

            ! Probes coordinates and grid id's
            CALL stencilio_write_master_cptr(arr_grouph, 'coordinates', &
                C_LOC(arr(iarr)%coordinates), &
                INT([3, arr(iarr)%npts], HSIZE_T), &
                mglet_hdf5_real)
            CALL stencilio_write_master_cptr(arr_grouph, 'igrid', &
                C_LOC(arr(iarr)%grid), INT([arr(iarr)%npts], HSIZE_T), &
                mglet_hdf5_int)

            DO ivar = 1, arr(iarr)%nvars
                npts = arr(iarr)%npts
                shape2 = [npts, 0_hsize_t]
                maxdims2 = [npts, H5S_UNLIMITED_F]
                chunksize2 = [MIN(chunksize21, npts), chunksize22]
                CALL hdf5common_dataset_create(arr(iarr)%variables(ivar), &
                    shape2, mglet_hdf5_real, arr_grouph, dset_id, &
                    filespace, maxdims2, chunksize2)
                CALL hdf5common_dataset_close(dset_id, filespace)
            END DO

            CALL h5gclose_f(arr_grouph, hdferr)
            IF (hdferr /= 0) CALL errr(__FILE__, __LINE__)
        END DO

        ! Initialize time
        shape1 = [0]
        maxdims1 = [H5S_UNLIMITED_F]
        chunksize1 = [chunksize11]
        CALL hdf5common_dataset_create("time", shape1, timeinfo_h5t, &
            probes_grouph, dset_id, filespace, maxdims1, chunksize1)
        CALL hdf5common_dataset_close(dset_id, filespace)

        CALL h5gclose_f(probes_grouph, hdferr)
        IF (hdferr /= 0) CALL errr(__FILE__, __LINE__)
    END SUBROUTINE init_file


    SUBROUTINE read_file(file_id, ittot)
        ! This subroutine read the time information from the
        ! probes file, to determine the correct offset for restarts

        ! Subroutine arguments
        INTEGER(HID_T), INTENT(inout) :: file_id
        INTEGER(intk), INTENT(in) :: ittot

        ! Local variables
        INTEGER(hid_t)  :: probes_id, dset_id, filespace
        INTEGER(hsize_t) :: i, dims1(1)
        TYPE(timeinfo_t), ALLOCATABLE, TARGET :: timeinfo_old(:)
        TYPE(C_PTR) :: cptr

        CALL hdf5common_group_open("PROBES", file_id, probes_id)

        ! Open and close to get dims - stencilio_read_master_cptr opens
        ! internally again.
        CALL hdf5common_dataset_open("time", dims1, probes_id, dset_id, &
            filespace)
        CALL hdf5common_dataset_close(dset_id, filespace)

        ! When it contains data
        IF (dims1(1) > 0) THEN
            IF (ioid == 0) THEN
                ALLOCATE(timeinfo_old(dims1(1)))
                cptr = C_LOC(timeinfo_old)
                CALL stencilio_read_master_cptr(probes_id, "time", cptr, &
                    dims1, timeinfo_h5t)

                ! The ittot given presently is the ittotbefore the first
                ! timestep of the present simulation is executed. First sampeld
                ! variable in the present simulation will be with ittot + 1.
                fileoffset = 0
                DO i = 1, dims1(1)
                    fileoffset = fileoffset + 1
                    IF (timeinfo_old(i)%ittot > ittot) EXIT
                END DO
                DEALLOCATE(timeinfo_old)
            END IF
            CALL MPI_Bcast(fileoffset, 1, mglet_mpi_hsize_t, 0, iocomm)
        END IF

        CALL hdf5common_group_close(probes_id)
    END SUBROUTINE read_file


    SUBROUTINE init_buffers(ittot, mtstep, itint, timeph, dt, tend)
        ! Subroutine arguments
        INTEGER(intk), INTENT(in) :: ittot
        INTEGER(intk), INTENT(in) :: mtstep
        INTEGER(intk), INTENT(in) :: itint
        REAL(realk), INTENT(in) :: timeph
        REAL(realk), INTENT(in) :: dt
        REAL(realk), INTENT(in) :: tend

        ! Subroutine arguments
        ! none...

        ! Local variables
        INTEGER(intk) :: iarr
        INTEGER(hsize_t) :: buflen_max
        REAL(realk) :: buflen_r
        INTEGER(int64) :: nvalues_max, nbytes

        ! Per-process number of values sampled
        nvalues_max = 0
        DO iarr = 1, narrays
            nvalues_max = nvalues_max + arr(iarr)%nvars*arr(iarr)%npts
        END DO

        ! Sum over IO group
        CALL MPI_Allreduce(MPI_IN_PLACE, nvalues_max, 1, MPI_INTEGER8, &
            MPI_SUM, iogrcomm)

        ! Global max
        CALL MPI_Allreduce(MPI_IN_PLACE, nvalues_max, 1, MPI_INTEGER8, &
            MPI_MAX, MPI_COMM_WORLD)

        ! Compute buffer length
        nbytes = real_bytes*nvalues_max
        buflen_r = REAL(maxbytes, kind=realk)/REAL(nbytes, kind=realk)

        ! Round downwards, only if you end up with a buffer length of 0
        ! increase to 1
        buflen = INT(buflen_r, kind=intk)
        IF (buflen == 0) buflen = 1

        ! Do not have a longer buffer than simulation length
        buflen_max = estimate_buflen(ittot, mtstep, itint, timeph, dt, tend)
        buflen = MIN(buflen, buflen_max)

        ! Nothing is sampled yet...
        bufloc = 1
        navg = 0
        fileoffset = 0

        ! Fileoffset might later in read_file be set in the case of restarts

        ! Allocate memory in buffers
        DO iarr = 1, narrays
            ALLOCATE(arr(iarr)%buffer(arr(iarr)%nmypnts, buflen, &
                arr(iarr)%nvars))
            arr(iarr)%buffer = 0.0
        END DO

        ! Allocate time
        IF (ioid == 0) ALLOCATE(timeinfo(buflen))
    END SUBROUTINE init_buffers


    FUNCTION estimate_buflen(ittot, mtstep, itint, timeph, dt, tend) &
            RESULT(buflen)
        ! Function arguments
        INTEGER(intk) :: buflen
        INTEGER(intk), INTENT(in) :: ittot
        INTEGER(intk), INTENT(in) :: mtstep
        INTEGER(intk), INTENT(in) :: itint
        REAL(realk), INTENT(in) :: timeph
        REAL(realk), INTENT(in) :: dt
        REAL(realk), INTENT(in) :: tend

        ! Local variables
        INTEGER(intk) :: buflen1
        REAL(realk) :: tstartprobes

        ! Stopping by mtstep
        buflen = mtstep

        ! Stopping by tend
        IF (tend > 0.0) THEN
            tstartprobes = MAX(tstart, timeph)
            buflen1 = INT((tend-tstartprobes)/dt, kind=intk)
            buflen = MIN(buflen, buflen1)
        END IF

        ! When no tend is present limit buffer by probes tstart
        IF (tend <= 0.0) THEN
            tstartprobes = MAX(tstart, timeph)
            buflen1 = INT((timeph+mtstep*dt-tstartprobes)/dt)
            buflen = MIN(buflen, buflen1)
        END IF

        ! Stopping by itint (rare use - here no tstart is considered)
        IF (itint > 0) THEN
            buflen1 = itint - ittot
            buflen = MIN(buflen, buflen1)
        END IF

        ! Always add one to avoid flushing buffers one timestep before
        ! simulation end
        buflen = buflen/itsamp + 1
    END FUNCTION estimate_buflen


    SUBROUTINE init_data_transfer(array)
        ! Subroutine arguments
        TYPE(probearr_t), INTENT(inout) :: array

        ! Local variables
        INTEGER(intk) :: i
        INTEGER(intk), ALLOCATABLE :: ngrpnts(:), displs(:), tmpgroupid(:), &
            tmppointowner(:), pointowner(:)

        ! Number of probes in IO group
        CALL MPI_Allreduce(array%nmypnts, array%grouppnts, 1, &
            mglet_mpi_int, MPI_SUM, iogrcomm)

        IF (ioproc) THEN
            ALLOCATE(ngrpnts(iogrprocs))
            ALLOCATE(displs(iogrprocs))
            ALLOCATE(tmpgroupid(array%grouppnts))
            ALLOCATE(tmppointowner(array%grouppnts))
        ELSE
            ALLOCATE(ngrpnts(0))
            ALLOCATE(displs(0))
            ALLOCATE(tmpgroupid(0))
            ALLOCATE(tmppointowner(0))
        END IF

        ! ngrpnts is a list of how many probes each MPI rank in the IO group
        ! has
        CALL MPI_Gather(array%nmypnts, 1, mglet_mpi_int, ngrpnts, 1, &
            mglet_mpi_int, 0, iogrcomm)

        IF (ioproc) THEN
            displs = 0
            DO i = 2, iogrprocs
                displs(i) = displs(i-1) + ngrpnts(i-1)
            END DO
        END IF

        ! tmpgroupid is the ID of each point in the IO group
        CALL MPI_Gatherv(array%id, array%nmypnts, mglet_mpi_int, &
            tmpgroupid, ngrpnts, displs, mglet_mpi_int, 0, iogrcomm)

        ! Local owner of points
        BLOCK
            INTEGER(intk), ALLOCATABLE :: localowner(:)

            ALLOCATE(localowner(array%nmypnts))
            localowner = iogrid
            CALL MPI_Gatherv(localowner, array%nmypnts, mglet_mpi_int, &
                tmppointowner, ngrpnts, displs, mglet_mpi_int, 0, &
                iogrcomm)
            DEALLOCATE(localowner)
        END BLOCK

        ! Non-IO processes can return here - no more work to do for them
        IF (.NOT. ioproc) THEN
            DEALLOCATE(tmppointowner)
            DEALLOCATE(tmpgroupid)
            DEALLOCATE(displs)
            DEALLOCATE(ngrpnts)
            RETURN
        END IF

        ! No need to check IF (ioproc) here any more
        DEALLOCATE(displs)

        ! Still allocated:
        !   - ngrpnts
        !   - tmpgroupid
        !   - tmppointowner

        ! IO process sort list of points
        ALLOCATE(array%groupid(array%grouppnts))
        BLOCK
            INTEGER(intk), ALLOCATABLE :: sortid(:)

            ALLOCATE(sortid(array%grouppnts))
            ALLOCATE(pointowner(array%grouppnts))

            CALL sortix(array%grouppnts, tmpgroupid, sortid)

            DO i = 1, array%grouppnts
                array%groupid(i) = tmpgroupid(sortid(i))
                pointowner(i) = tmppointowner(sortid(i))
            END DO

            DEALLOCATE(sortid)
        END BLOCK

        DEALLOCATE(tmppointowner)
        DEALLOCATE(tmpgroupid)
        DEALLOCATE(ngrpnts)

        ! Still allocated:
        !   - pointowner

        BLOCK
            INTEGER(intk) :: nblocks, proc
            INTEGER(intk), ALLOCATABLE :: blocklengths(:), displacements(:)
            INTEGER(KIND=mpi_address_kind) :: lb, extent
            TYPE(MPI_Datatype) :: tmptype

            ALLOCATE(blocklengths(array%grouppnts))
            ALLOCATE(displacements(array%grouppnts))
            ALLOCATE(array%mpitype(0:iogrprocs-1))

            DO proc = 0, iogrprocs-1
                nblocks = 0
                blocklengths = 0
                displacements = 0

                DO i = 1, array%grouppnts
                    IF (pointowner(i) == proc) THEN
                        ! Start new block?
                        IF (blocklengths(nblocks + 1) == 0) THEN
                            displacements(nblocks + 1) = (i - 1)
                        END IF

                        blocklengths(nblocks + 1) = blocklengths(nblocks + 1) &
                            + 1

                        ! Finish current block?
                        IF (i < array%grouppnts) THEN
                            IF (pointowner(i+1) /= proc) THEN
                                nblocks = nblocks + 1
                            END IF
                        ELSE
                            nblocks = nblocks + 1
                        END IF
                    END IF
                END DO

                ! Create temporary datatype
                CALL MPI_Type_indexed(nblocks, blocklengths, displacements, &
                    mglet_mpi_real, tmptype)

                ! Create new type with correct extent
                CALL MPI_Type_get_extent(tmptype, lb, extent)
                extent = array%grouppnts*real_bytes
                CALL MPI_Type_create_resized(tmptype, lb, extent, &
                    array%mpitype(proc))

                ! Commit type and free temporary type
                CALL MPI_Type_commit(array%mpitype(proc))
                CALL MPI_Type_free(tmptype)
            END DO

            DEALLOCATE(blocklengths)
            DEALLOCATE(displacements)
        END BLOCK

        ! Deallocate last array
        DEALLOCATE(pointowner)
    END SUBROUTINE init_data_transfer


    SUBROUTINE write_probes()
        ! Subroutine arguments
        ! none...

        ! Local variables
        INTEGER(hid_t) :: file_id, probes_id
        INTEGER(intk) :: iarr

        ! If navg == 0 then the bufloc is not a valid sample, decrement by one
        IF (navg == 0) THEN
            bufloc = bufloc - 1
        END IF

        ! Only write if there is data in the buffer
        ! (bufloc is only 0 in the case that it was decremented above, so reset
        ! it to 1 here before returning)
        IF (bufloc == 0) THEN
            bufloc = 1
            navg = 0
            RETURN
        END IF

        IF (myid == 0) THEN
            WRITE(*, "('Writing probes to disk')")
            WRITE(*, "()")
        END IF

        CALL hdf5common_open(outfile, "r+", file_id)
        CALL hdf5common_group_open("PROBES", file_id, probes_id)

        DO iarr = 1, narrays
            CALL write_array(probes_id, arr(iarr))
        END DO

        CALL write_time(probes_id)

        CALL hdf5common_group_close(probes_id)
        CALL hdf5common_close(file_id)

        ! Store file offset for next use, reset counters
        fileoffset = fileoffset + bufloc
        bufloc = 1
        navg = 0
    END SUBROUTINE write_probes


    SUBROUTINE write_array(probes_id, array)
        ! Subroutine arguments
        INTEGER(hid_t), INTENT(inout) :: probes_id
        TYPE(probearr_t), INTENT(inout) :: array

        ! Local variables
        INTEGER(hid_t) :: arr_id
        INTEGER(intk) :: ivar

        CALL hdf5common_group_open(array%name, probes_id, arr_id)

        DO ivar = 1, array%nvars
            CALL write_variable(arr_id, array, ivar)
        END DO

        CALL hdf5common_group_close(arr_id)
    END SUBROUTINE write_array


    SUBROUTINE write_variable(arr_id, array, ivar)
        USE ISO_C_BINDING, ONLY: C_PTR, C_LOC

        ! Subroutine arguments
        INTEGER(hid_t), INTENT(inout) :: arr_id
        TYPE(probearr_t), INTENT(inout) :: array
        INTEGER(intk), INTENT(in) :: ivar

        ! Local variables
        INTEGER(hid_t)  :: dset_id, filespace, memspace, dxpl
        TYPE(MPI_Request), ALLOCATABLE :: recvreq(:)
        TYPE(MPI_Request) :: sendreq
        TYPE(MPI_Status), ALLOCATABLE :: recvstatus(:)
        TYPE(MPI_Status) :: sendstatus
        REAL(realk), POINTER, CONTIGUOUS :: outputbuf(:, :)
        INTEGER(intk) :: i
        INTEGER(int32) :: hdferr, nelems
        INTEGER(hsize_t) :: npts, dims2(2), offset2(2), count2(2)
        TYPE(C_PTR) :: cptr
        INTEGER :: op
        LOGICAL :: finish_block

        ! Gather data on IO process
        IF (ioproc) THEN
            ! This will not do anything if buffer is already bif enough
            CALL increase_bigbuf(array%grouppnts*bufloc)

            ! Convenience for having a 2D buffer
            outputbuf(1:array%grouppnts, 1:bufloc) => &
                bigbuf(1:array%grouppnts*bufloc)

            ALLOCATE(recvreq(0:iogrprocs-1))
            DO i = 0, iogrprocs-1
                CALL MPI_Irecv(outputbuf, INT(bufloc, int32), &
                    array%mpitype(i), INT(i, int32), 0, iogrcomm, &
                    recvreq(i))
            END DO
        END IF

        ! All processes send data to IO process
        nelems = INT(array%nmypnts*bufloc, int32)
        CALL MPI_Isend(array%buffer(:, :, ivar), nelems, mglet_mpi_real, &
            0, 0, iogrcomm, sendreq)

        ! Wait for communication to finish
        IF (ioproc) THEN
            ALLOCATE(recvstatus(iogrprocs))
            CALL MPI_Waitall(iogrprocs, recvreq, recvstatus)
            DEALLOCATE(recvstatus)
            DEALLOCATE(recvreq)
        END IF
        CALL MPI_Wait(sendreq, sendstatus)

        ! Non-IO processes return here, they are finished with the job
        IF (.NOT. ioproc) RETURN

        CALL hdf5common_dataset_open(TRIM(array%variables(ivar)), dims2, &
            arr_id, dset_id, filespace)

        ! Extend dataset
        npts = INT(array%npts, hsize_t)
        dims2 = [npts, fileoffset + bufloc]
        CALL hdf5common_dataset_extend(dset_id, filespace, dims2)

        ! Select hyperslab in filespace, the dimension are as follows:
        !    1: probe id
        !    2: time
        offset2 = [0_hsize_t, fileoffset]
        count2 = [0_hsize_t, bufloc]
        finish_block = .FALSE.
        op = H5S_SELECT_SET_F
        DO i = 1, array%grouppnts
            ! Start new block?
            IF (count2(1) == 0) THEN
                offset2(1) = array%groupid(i) - 1
            END IF

            count2(1) = count2(1) + 1

            ! Finish current block?
            IF (i < array%grouppnts) THEN
                IF (array%groupid(i+1) /= array%groupid(i) + 1) THEN
                    finish_block = .TRUE.
                END IF
            ELSE
                finish_block = .TRUE.
            END IF

            IF (finish_block) THEN
                CALL h5sselect_hyperslab_f(filespace, op, offset2, count2, &
                    hdferr)
                IF (hdferr /= 0) CALL errr(__FILE__, __LINE__)
                offset2(1) = 0
                count2(1) = 0
                op = H5S_SELECT_OR_F
                finish_block = .FALSE.
            END IF
        END DO

        ! Create and select hyperslab in memspace
        dims2 = [INT(array%grouppnts, hsize_t), bufloc]
        offset2 = 0
        count2 = [INT(array%grouppnts, hsize_t), bufloc]
        CALL h5screate_simple_f(2, dims2, memspace, hdferr)
        IF (hdferr /= 0) CALL errr(__FILE__, __LINE__)
        CALL h5sselect_hyperslab_f(memspace, H5S_SELECT_SET_F, offset2, &
            count2, hdferr)
        IF (hdferr /= 0) CALL errr(__FILE__, __LINE__)

        ! If a process has no points, select a none as filespace and memspace
        ! (special)
        IF (array%grouppnts == 0) THEN
            CALL h5sselect_none_f(filespace, hdferr)
            IF (hdferr /= 0) CALL errr(__FILE__, __LINE__)
            CALL h5sselect_none_f(memspace, hdferr)
            IF (hdferr /= 0) CALL errr(__FILE__, __LINE__)
        END IF

        ! Write dataset
        CALL h5pcreate_f(H5P_DATASET_XFER_F, dxpl, hdferr)
        IF (hdferr /= 0) CALL errr(__FILE__, __LINE__)
        CALL h5pset_dxpl_mpio_f(dxpl, hdf5_io_mode, hdferr)
        IF (hdferr /= 0) CALL errr(__FILE__, __LINE__)

        cptr = C_LOC(outputbuf)
        CALL h5dwrite_f(dset_id, mglet_hdf5_real, cptr, hdferr, &
            file_space_id=filespace, mem_space_id=memspace, &
            xfer_prp=dxpl)
        IF (hdferr /= 0) CALL errr(__FILE__, __LINE__)

        CALL h5pclose_f(dxpl, hdferr)
        IF (hdferr /= 0) CALL errr(__FILE__, __LINE__)

        ! Close leftover handles
        CALL h5sclose_f(memspace, hdferr)
        IF (hdferr /= 0) CALL errr(__FILE__, __LINE__)
        CALL hdf5common_dataset_close(dset_id, filespace)

        NULLIFY(outputbuf)
    END SUBROUTINE write_variable


    SUBROUTINE write_time(probes_id)
        ! Subroutine arguments
        INTEGER(hid_t), INTENT(inout) :: probes_id

        ! Local variables
        INTEGER(hid_t) :: dset_id, memspace, filespace
        INTEGER(hsize_t) :: dims1(1), offset1(1), count1(1)
        INTEGER(int32) :: hdferr

        IF (.NOT. ioproc) RETURN

        ! Memspace
        dims1 = [bufloc]
        CALL h5screate_simple_f(1, dims1, memspace, hdferr)
        IF (hdferr /= 0) CALL errr(__FILE__, __LINE__)

        IF (ioid == 0) THEN
            CALL h5sselect_all_f(memspace, hdferr)
            IF (hdferr /= 0) CALL errr(__FILE__, __LINE__)
        ELSE
            CALL h5sselect_none_f(memspace, hdferr)
            IF (hdferr /= 0) CALL errr(__FILE__, __LINE__)
        END IF

        CALL hdf5common_dataset_open("time", dims1, probes_id, dset_id, &
            filespace)

        ! Filespace
        offset1 = [fileoffset]
        count1 = [bufloc]
        dims1 = [fileoffset + bufloc]
        CALL hdf5common_dataset_extend(dset_id, filespace, dims1)
        IF (ioid == 0) THEN
            CALL h5sselect_hyperslab_f(filespace, H5S_SELECT_SET_F, offset1, &
                count1, hdferr)
            IF (hdferr /= 0) CALL errr(__FILE__, __LINE__)
        ELSE
            CALL h5sselect_none_f(filespace, hdferr)
            IF (hdferr /= 0) CALL errr(__FILE__, __LINE__)
        END IF

        ! Write actual data
        CALL h5dwrite_f(dset_id, timeinfo_h5t, C_LOC(timeinfo), hdferr, &
            file_space_id=filespace, mem_space_id=memspace)
        IF (hdferr /= 0) CALL errr(__FILE__, __LINE__)

        CALL hdf5common_dataset_close(dset_id, filespace)

        CALL h5sclose_f(memspace, hdferr)
        IF (hdferr /= 0) CALL errr(__FILE__, __LINE__)

    END SUBROUTINE write_time
END MODULE probes_mod
