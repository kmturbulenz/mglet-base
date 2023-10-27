MODULE fillfollower_mod
    USE MPI_f08

    USE core_mod, ONLY: realk, intk, int32, idim3d, get_mgdims, get_ip3, &
        nmygrids, mygrids, field_t, mglet_mpi_int, get_field, errr, myid
    USE filling_mod, ONLY: blockquad_search

    IMPLICIT NONE (type, external)
    PRIVATE

    PUBLIC :: fillfollower

CONTAINS
    SUBROUTINE fillfollower(knoten, nofluidpoints)
        ! Subroutine arguments
        TYPE(field_t), INTENT(inout) :: knoten
        REAL(realk), INTENT(in) :: nofluidpoints(:, :)

        ! Local variables
        INTEGER(intk) :: i, igrid, kk, jj, ii, ip3, idx
        INTEGER(intk) :: nnofluidpoints, nfound, ninfluid, ninfluid_tot
        INTEGER(intk), ALLOCATABLE :: inofluidpoints(:, :), idxinfluid(:)
        INTEGER(intk), ALLOCATABLE :: foundinfluid(:)
        TYPE(field_t), POINTER :: x_f, y_f, z_f
        REAL(realk), POINTER, CONTIGUOUS :: x(:), y(:), z(:)

        ! Get gridspacing
        CALL get_field(x_f, "X")
        CALL get_field(y_f, "Y")
        CALL get_field(z_f, "Z")

        nnofluidpoints = SIZE(nofluidpoints, 2)
        ALLOCATE(inofluidpoints(3, nnofluidpoints))
        inofluidpoints = 0

        ALLOCATE(idxinfluid(nnofluidpoints))
        idxinfluid = 0

        ALLOCATE(foundinfluid(nnofluidpoints))
        foundinfluid = 0

        ! Identify points in fluid
        DO i = 1, nmygrids
            igrid = mygrids(i)

            CALL get_mgdims(kk, jj, ii, igrid)
            CALL get_ip3(ip3, igrid)
            CALL x_f%get_ptr(x, igrid)
            CALL y_f%get_ptr(y, igrid)
            CALL z_f%get_ptr(z, igrid)

            CALL blockquad_search(nnofluidpoints, nfound, inofluidpoints, &
                nofluidpoints, kk, jj, ii, x, y, z)

            ! In contrast to filling algorithm, there is no need to continue
            ! if none of the nofluidpoints are in this grid
            IF (nfound == 0) CYCLE

            ! ninfluid is the number of points in fluid
            ! nfound is the number of found seeding points in this grid
            CALL checknopoint(kk, jj, ii, knoten%arr(ip3), nfound, &
                inofluidpoints, ninfluid, idxinfluid)

            ! Each point can be found in different grids and by different
            ! ranks. Mark as found here and then allreduce afterwards
            DO idx = 1, ninfluid
                foundinfluid(idxinfluid(idx)) = 1
            END DO
        END DO

        ! Count overall number of points found
        ! Must be allreduce, in that way all processes can count if some points
        ! were found or not
        CALL MPI_Allreduce(MPI_IN_PLACE, foundinfluid, &
            INT(nnofluidpoints, int32), mglet_mpi_int, &
            MPI_MAX, MPI_COMM_WORLD)

        ! Write out diagnostic information
        ninfluid_tot = SUM(foundinfluid)
        IF (myid == 0) THEN
            WRITE(*, '(" Found ", I0, " points in fluid")') ninfluid_tot
        END IF

        ! If no points are locaed in the fluid, return and continue blocking
        IF (ninfluid_tot == 0) RETURN

        ! Write out index and coordinates of points found in fluid
        IF (myid == 0) THEN
            WRITE(*, '(" The following points are in the fluid:")')
            DO idx = 1, nnofluidpoints
                IF (foundinfluid(idx) > 0) THEN
                    WRITE(*, '(4X, I4, F16.3, F16.3, F16.3)') idx, &
                        nofluidpoints(1, idx), nofluidpoints(2, idx), &
                        nofluidpoints(3, idx)
                END IF
            END DO
        END IF

        ! Now identify the fluid paths from the nofluidpoint to the filling
        ! start
        ! TODO: Implement

        ! TODO: Make a graceful and nice stop instead???
        CALL errr(__FILE__, __LINE__)
    END SUBROUTINE fillfollower


    PURE SUBROUTINE checknopoint(kk, jj, ii, knoten, npts, i0, &
            ninfluid, idxinfluid)
        ! Subroutine arguments
        INTEGER(intk), INTENT(in) :: kk, jj, ii
        REAL(realk), INTENT(in) :: knoten(kk, jj, ii)
        INTEGER(intk), INTENT(in) :: npts
        INTEGER(intk), INTENT(in) :: i0(3, npts)
        INTEGER(intk), INTENT(out) :: ninfluid
        INTEGER(intk), INTENT(out) :: idxinfluid(npts)

        ! Local variables
        INTEGER(intk) :: idx, k, j, i

        ninfluid = 0
        idxinfluid = 0
        DO idx = 1, npts
            i = i0(1, idx)
            j = i0(2, idx)
            k = i0(3, idx)

            IF (knoten(k, j, i) > 0) THEN
                ninfluid = ninfluid + 1
                idxinfluid(ninfluid) = idx
            END IF
        END DO
    END SUBROUTINE checknopoint
END MODULE fillfollower_mod
