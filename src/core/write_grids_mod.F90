MODULE write_grids_mod
    USE HDF5

    USE precision_mod, ONLY: intk, realk
    USE field_mod, ONLY: field_t
    USE fields_mod, ONLY: get_field
    USE gridio_mod, ONLY: gridinfo_t, bcond_t, write_gridinfo, write_bcondinfo
    USE grids_mod, ONLY: rewrite_grids
    USE comms_mod, ONLY: myid, ioproc
    USE fieldio2_mod, ONLY: init_fieldio, fieldio_write, finish_fieldio
    USE hdf5common_mod, ONLY: hdf5common_open, hdf5common_close

    IMPLICIT NONE (type, external)
    PRIVATE

    ! Public subroutines
    PUBLIC :: write_grids

CONTAINS

    SUBROUTINE write_grids(nfluidcells, outfile)
        ! Subroutine arguments
        INTEGER(intk), INTENT(in) :: nfluidcells(:)
        CHARACTER(len=*), INTENT(in) :: outfile

        ! Local variables
        INTEGER(intk) :: new_ngrid
        INTEGER(HID_T) :: file_id
        TYPE(gridinfo_t), ALLOCATABLE :: gridinfo(:)
        TYPE(bcond_t), ALLOCATABLE :: front(:)
        TYPE(bcond_t), ALLOCATABLE :: back(:)
        TYPE(bcond_t), ALLOCATABLE :: right(:)
        TYPE(bcond_t), ALLOCATABLE :: left(:)
        TYPE(bcond_t), ALLOCATABLE :: bottom(:)
        TYPE(bcond_t), ALLOCATABLE :: top(:)
        REAL(realk), ALLOCATABLE :: realprms(:)
        INTEGER(intk), ALLOCATABLE :: intprms(:)
        INTEGER(intk), ALLOCATABLE :: new_gridids(:)

        CALL rewrite_grids(gridinfo, front, back, right, left, bottom, top, &
            realprms, intprms, new_gridids, new_ngrid, nfluidcells)

        IF (myid == 0) THEN
            WRITE(*, '(1X, "Checkblock writing ", I0, " grids into ", A)') &
                new_ngrid, TRIM(outfile)
        END IF

        CALL hdf5common_open(outfile, "w", file_id)

        ! Write gridinfo and boundary condition tables
        IF (ioproc) THEN
            CALL write_gridinfo(file_id, gridinfo, realprms, intprms, new_ngrid)
            CALL write_bcondinfo(file_id, front, "FRONT", new_ngrid)
            CALL write_bcondinfo(file_id, back, "BACK", new_ngrid)
            CALL write_bcondinfo(file_id, right, "RIGHT", new_ngrid)
            CALL write_bcondinfo(file_id, left, "LEFT", new_ngrid)
            CALL write_bcondinfo(file_id, bottom, "BOTTOM", new_ngrid)
            CALL write_bcondinfo(file_id, top, "TOP", new_ngrid)
        END IF

        IF (myid == 0) THEN
            DEALLOCATE(gridinfo)
            DEALLOCATE(front)
            DEALLOCATE(back)
            DEALLOCATE(right)
            DEALLOCATE(left)
            DEALLOCATE(bottom)
            DEALLOCATE(top)
        END IF

        ! Write gridspacing arrays
        CALL init_fieldio(new_gridids)

        CALL write_gridarray(file_id, "X")
        CALL write_gridarray(file_id, "DX")
        CALL write_gridarray(file_id, "DDX")

        CALL write_gridarray(file_id, "Y")
        CALL write_gridarray(file_id, "DY")
        CALL write_gridarray(file_id, "DDY")

        CALL write_gridarray(file_id, "Z")
        CALL write_gridarray(file_id, "DZ")
        CALL write_gridarray(file_id, "DDZ")

        CALL finish_fieldio()

        DEALLOCATE(realprms)
        DEALLOCATE(intprms)
        DEALLOCATE(new_gridids)

        CALL hdf5common_close(file_id)
    END SUBROUTINE write_grids


    SUBROUTINE write_gridarray(parent_id, name)
        ! Subroutine arguments
        INTEGER(HID_T), INTENT(inout) :: parent_id
        CHARACTER(len=*), INTENT(in) :: name

        ! Local variables
        TYPE(field_t), POINTER :: field

        CALL get_field(field, name)
        CALL fieldio_write(parent_id, field)
    END SUBROUTINE write_gridarray

END MODULE write_grids_mod
