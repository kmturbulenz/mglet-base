MODULE snapshots_mod
    USE HDF5

    USE core_mod

    IMPLICIT NONE(type, external)
    PRIVATE

    CHARACTER(len=mglet_filename_max) :: outfile
    INTEGER(intk) :: itsamp
    REAL(realk) :: tstart
    INTEGER(intk) :: nfields
    CHARACTER(len=nchar_name), ALLOCATABLE :: fieldlist(:)

    PUBLIC :: init_snapshots, sample_snapshots, finish_snapshots

CONTAINS
    SUBROUTINE init_snapshots(ittot, mtstep, itint, timeph, dt, tend)
        ! Subroutine arguments
        INTEGER(intk), INTENT(in) :: ittot
        INTEGER(intk), INTENT(in) :: mtstep
        INTEGER(intk), INTENT(in) :: itint
        REAL(realk), INTENT(in) :: timeph
        REAL(realk), INTENT(in) :: dt
        REAL(realk), INTENT(in) :: tend

        ! Local variables
        TYPE(field_t), POINTER :: field
        CHARACTER(len=64) :: jsonptr
        INTEGER(intk) :: i
        INTEGER(hid_t) :: file_id, group_id

        nfields = 0
        IF (.NOT. fort7%exists("/snapshots")) THEN
            RETURN
        END IF

        CALL fort7%get_size("/snapshots/fields", nfields)
        IF (nfields <= 0) RETURN

        CALL fort7%get_value("/snapshots/file", outfile, "snapshots.h5")

        CALL fort7%get_value("/snapshots/itsamp", itsamp)
        IF (itsamp <= 0) THEN
            WRITE(*, '("snapshots: itsamp must be positive: ", I7)') itsamp
            CALL errr(__FILE__, __LINE__)
        END IF

        CALL fort7%get_value("/tstart", tstart, 0.0)
        IF (tstart < 0.0) THEN
            WRITE(*, '("snapshots: tstart must be >= 0.0: ", F15.7)') tstart
            CALL errr(__FILE__, __LINE__)
        END IF

        ALLOCATE(fieldlist(nfields))
        DO i = 1, nfields
            WRITE(jsonptr, '("/snapshots/fields/", I0)') i-1
            CALL fort7%get_value(jsonptr, fieldlist(i))

            ! Sanity check if the field does not exist
            ! Not the most beautiful error message, though...
            CALL get_field(field, fieldlist(i))
        END DO

        IF (myid == 0) THEN
            WRITE(*, '("ENABLING SNAPSHOTS:")')
            WRITE(*, '(4X, "file:   ", A)') TRIM(outfile)
            WRITE(*, '(4X, "itsamp: ", I15)') itsamp
            WRITE(*, '(4X, "tstart: ", F15.7)') tstart
            WRITE(*, '(4X, "nfields: ", I15)') nfields
            DO i = 1, nfields
                WRITE(*, '(8X, A)') fieldlist(i)
            END DO
            WRITE(*, '()')
        END IF

        ! Prepare file if not continuing from another run
        IF (.NOT. dcont) THEN
            CALL hdf5common_open(outfile, "w", file_id)
            CALL hdf5common_group_open("SNAPSHOTS", file_id, group_id, &
                track_index=.TRUE.)

            CALL hdf5common_group_close(group_id)
            CALL hdf5common_close(file_id)
        END IF
    END SUBROUTINE init_snapshots


    SUBROUTINE finish_snapshots()
        IF (nfields <= 0) RETURN

        DEALLOCATE(fieldlist)
    END SUBROUTINE finish_snapshots


    SUBROUTINE sample_snapshots(itstep, ittot, timeph, dt)
        ! Subroutine arguments
        INTEGER(intk), INTENT(in) :: itstep
        INTEGER(intk), INTENT(in) :: ittot
        REAL(realk), INTENT(in) :: timeph
        REAL(realk), INTENT(in) :: dt

        ! Local variables
        TYPE(field_t), POINTER :: field
        CHARACTER(len=64) :: grpname
        INTEGER(intk) :: i
        INTEGER(hid_t) :: file_id, group1_id, group2_id

        IF (nfields <= 0) RETURN
        IF (MOD(ittot, itsamp) /= 0) RETURN
        IF (myid == 0) THEN
            WRITE(*, '("Writing snapshot at ittot = ", I0)') ittot
            WRITE(*, '()')
        END IF

        ! File was already prepared in init_snapshots
        CALL hdf5common_open(outfile, "a", file_id)
        CALL hdf5common_group_open("SNAPSHOTS", file_id, group1_id, &
            track_index=.TRUE.)

        WRITE(grpname, '("ITTOT-", I0)') ittot
        CALL hdf5common_group_open(grpname, group1_id, group2_id, &
            track_index=.TRUE.)

        CALL hdf5common_attr_write("ITSTEP", itstep, group2_id)
        CALL hdf5common_attr_write("ITTOT", ittot, group2_id)
        CALL hdf5common_attr_write("TIMEPH", timeph, group2_id)
        CALL hdf5common_attr_write("DT", dt, group2_id)

        CALL init_fieldio()
        DO i = 1, nfields
            CALL get_field(field, fieldlist(i))
            CALL fieldio_write(group2_id, field)
        END DO
        CALL finish_fieldio()

        CALL hdf5common_group_close(group2_id)
        CALL hdf5common_group_close(group1_id)
        CALL hdf5common_close(file_id)
    END SUBROUTINE sample_snapshots
END MODULE snapshots_mod
