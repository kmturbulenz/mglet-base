MODULE fields_mod
    USE HDF5
    USE comms_mod, ONLY: myid
    USE err_mod, ONLY: errr
    USE fieldio2_mod, ONLY: fieldio_read, fieldio_write
    USE fort7_mod
    USE hdf5common_mod, ONLY: hdf5common_open, hdf5common_close, &
        hdf5common_group_open, hdf5common_group_close
    USE precision_mod
    USE field_mod, ONLY: field_t, get_len_i, nchar_name

    IMPLICIT NONE(type, external)
    PRIVATE

    ! Fields storage array
    INTEGER(intk), PARAMETER :: nfields_max = 1000
    INTEGER(intk) :: nfields = 0
    TYPE(field_t), TARGET :: fields(nfields_max)

    ! Field input/output storage
    CHARACTER(len=mglet_filename_max) :: filename = &
        REPEAT(" ", mglet_filename_max)
    INTEGER(hid_t) :: file_id = 0
    INTEGER(hid_t) :: group_id = 0

    ! Mode for writing result files
    ! Flipepd to "a" if:
    !  - Any time when dcont = .T.
    !  - _After_ 1st write operation otherwise
    ! this allows for checkpoints to function properly...
    CHARACTER(len=1) :: writemode = "w"

    INTERFACE get_fieldptr
        MODULE PROCEDURE :: get_fieldptr_grid1
        MODULE PROCEDURE :: get_fieldptr_grid3
        MODULE PROCEDURE :: get_fieldptr
    END INTERFACE get_fieldptr

    INTERFACE get_field
        MODULE PROCEDURE :: get_field_A
        MODULE PROCEDURE :: get_field_B
    END INTERFACE get_field

    PUBLIC :: init_fields, finish_fields, set_field, get_field, get_fieldptr, &
        fields_begin_read, fields_begin_write, fields_write, &
        fields_end_rw, get_fileh

CONTAINS
    SUBROUTINE init_fields()
        CONTINUE
    END SUBROUTINE init_fields


    SUBROUTINE finish_fields()
        ! Cleans up memory. Strictly speaking not neccesary, but it is ugly
        ! and triggers lots of messages in valgrind when memory is not
        ! explicitly deallocated in the end of executing a program.

        ! Local variables
        INTEGER(intk) :: i

        DO i = 1, nfields
            CALL fields(i)%finish()
        END DO
    END SUBROUTINE finish_fields


    SUBROUTINE fields_begin_read()
        IF (.NOT. dread) RETURN

        CALL fort7%get_value("/io/infile", filename, "fields.h5")

        IF (myid == 0) THEN
            WRITE(*, '("Opening file ", A, " for reading")') TRIM(filename)
        END IF

        CALL hdf5common_open(filename, 'r', file_id)
        CALL hdf5common_group_open("VOLUMEFIELDS", file_id, group_id, &
            track_index=.TRUE.)
    END SUBROUTINE fields_begin_read


    SUBROUTINE fields_begin_write()
        ! Local variables
        ! none...

        IF (.NOT. dwrite) RETURN

        CALL fort7%get_value("/io/outfile", filename, "fields.h5")

        IF (dcont) THEN
            writemode = "a"
        END IF

        IF (myid == 0) THEN
            WRITE(*, '("Opening file ", A, " with mode ", A)') &
                TRIM(filename), writemode
        END IF

        CALL hdf5common_open(filename, writemode, file_id)
        CALL hdf5common_group_open("VOLUMEFIELDS", file_id, group_id, &
            track_index=.TRUE.)

        ! Next time opened for writing the mode is "a"
        writemode = "a"
    END SUBROUTINE fields_begin_write


    SUBROUTINE fields_write()
        USE MPI_f08, ONLY: MPI_Wtime

        ! Local variables
        INTEGER(intk) :: i
        REAL(real64) :: tic, toc

        IF (.NOT. dwrite) RETURN

        tic = MPI_Wtime()

        ! Iterate over all fields and write out
        DO i = 1, nfields
            IF (fields(i)%dwrite) THEN
                CALL fieldio_write(group_id, fields(i))
            END IF
        END DO

        toc = MPI_Wtime()
        IF (myid == 0) THEN
            WRITE(*, '("Wrote all fields in ", F6.3, " seconds")') &
                toc - tic
            WRITE(*, '()')
        END IF
    END SUBROUTINE fields_write


    SUBROUTINE fields_end_rw()
        IF (file_id == 0) RETURN

        CALL hdf5common_group_close(group_id)
        CALL hdf5common_close(file_id)

        filename = REPEAT(" ", mglet_filename_max)
        file_id = 0
        group_id = 0
    END SUBROUTINE fields_end_rw


    SUBROUTINE set_field(name, description, ndim, istag, jstag, kstag, &
            units, dread, required, dwrite, buffers, active_level, get_len)

        ! Subroutine arguments
        CHARACTER(len=*), INTENT(in) :: name
        CHARACTER(len=*), INTENT(in), OPTIONAL :: description
        INTEGER(intk), INTENT(in), OPTIONAL :: ndim
        INTEGER(intk), INTENT(in), OPTIONAL :: istag
        INTEGER(intk), INTENT(in), OPTIONAL :: jstag
        INTEGER(intk), INTENT(in), OPTIONAL :: kstag
        INTEGER(intk), INTENT(in), OPTIONAL :: units(7)
        LOGICAL, INTENT(in), OPTIONAL :: dread
        LOGICAL, INTENT(in), OPTIONAL :: required
        LOGICAL, INTENT(in), OPTIONAL :: dwrite
        LOGICAL, INTENT(in), OPTIONAL :: buffers
        LOGICAL, INTENT(in), OPTIONAL :: active_level(:)
        PROCEDURE(get_len_i), OPTIONAL :: get_len

        ! Local variables
        TYPE(field_t), POINTER :: dummy
        LOGICAL :: exists

        ! Check that field does not exist before
        CALL get_field(dummy, name, exists)
        IF (exists) THEN
            WRITE(*, '("Field ", A, " exists already!")') name
            CALL errr(__FILE__, __LINE__)
        END IF

        IF (nfields + 1 > nfields_max) CALL errr(__FILE__, __LINE__)
        nfields = nfields + 1

        ! Initialize fields
        CALL fields(nfields)%init(name, description, ndim, istag, &
            jstag, kstag, units, dread, required, dwrite, &
            active_level, get_len)

        ! Allocate buffers
        IF (PRESENT(buffers)) THEN
            IF (buffers) CALL fields(nfields)%init_buffers()
        END IF

        ! Read if neccesary
        IF (fields(nfields)%dread) THEN
            IF (group_id /= 0) THEN
                CALL fieldio_read(group_id, fields(nfields), required)
            ELSE
                CALL errr(__FILE__, __LINE__)
            END IF
        END IF
    END SUBROUTINE set_field


    SUBROUTINE get_field_A(field, name, found)
        ! TODO: Implement dictionary lookup
        TYPE(field_t), POINTER, INTENT(inout) :: field
        CHARACTER(len=*), INTENT(in) :: name
        LOGICAL, OPTIONAL, INTENT(out) :: found

        INTEGER(intk) :: i
        LOGICAL :: thisfound
        CHARACTER(len=nchar_name) :: name2

        IF (LEN_TRIM(name) > nchar_name) CALL errr(__FILE__, __LINE__)

        ! Avoid calling TRIM in a loop, copy input fieldname to a
        ! variable with fixed length
        name2 = name

        NULLIFY(field)
        thisfound = .FALSE.
        DO i = 1, nfields
            IF (fields(i)%name == name2) THEN
                field => fields(i)
                thisfound = .TRUE.
            END IF
        END DO

        IF (PRESENT(found)) THEN
            found = thisfound
        ELSE
            IF (.NOT. thisfound) THEN
                WRITE(*,*) "Could not find field: ", TRIM(name)
                CALL errr(__FILE__, __LINE__)
            END IF
        END IF
    END SUBROUTINE get_field_A


    SUBROUTINE get_field_B(fieldarr, name, found)
        ! TODO: Implement dictionary lookup
        REAL(realk), POINTER, INTENT(inout) :: fieldarr(:)
        CHARACTER(len=*), INTENT(in) :: name
        LOGICAL, OPTIONAL, INTENT(out) :: found

        INTEGER(intk) :: i
        LOGICAL :: thisfound
        CHARACTER(len=nchar_name) :: name2

        IF (LEN_TRIM(name) > nchar_name) CALL errr(__FILE__, __LINE__)

        ! Avoid calling TRIM in a loop, copy input fieldname to a
        ! variable with fixed length
        name2 = name

        NULLIFY(fieldarr)
        thisfound = .FALSE.
        DO i = 1, nfields
            IF (fields(i)%name == name2) THEN
                fieldarr => fields(i)%arr
                thisfound = .TRUE.
            END IF
        END DO

        IF (PRESENT(found)) THEN
            found = thisfound
        ELSE
            IF (.NOT. thisfound) THEN
                WRITE(*,*) "Could not find field: ", TRIM(name)
                CALL errr(__FILE__, __LINE__)
            END IF
        END IF
    END SUBROUTINE get_field_B


    SUBROUTINE get_fieldptr_grid1(ptr, name, igrid)
        ! Subroutine arguments
        REAL(realk), INTENT(out), POINTER, CONTIGUOUS :: ptr(:)
        CHARACTER(len=*), INTENT(in) :: name
        INTEGER(intk), INTENT(in) :: igrid

        ! Local variables
        TYPE(field_t), POINTER :: field
        LOGICAL :: found

        CALL get_field(field, name, found)

        IF (.NOT. found) THEN
            WRITE(*, '("Field ", A, " does not exist!")') name
            CALL errr(__FILE__, __LINE__)
        END IF

        CALL field%get_ptr(ptr, igrid)
    END SUBROUTINE get_fieldptr_grid1


    SUBROUTINE get_fieldptr_grid3(ptr, name, igrid)
        ! Subroutine arguments
        REAL(realk), INTENT(out), POINTER, CONTIGUOUS :: ptr(:, :, :)
        CHARACTER(len=*), INTENT(in) :: name
        INTEGER(intk), INTENT(in) :: igrid

        ! Local variables
        TYPE(field_t), POINTER :: field
        LOGICAL :: found

        CALL get_field(field, name, found)

        IF (.NOT. found) THEN
            WRITE(*, '("Field ", A, " does not exist!")') name
            CALL errr(__FILE__, __LINE__)
        END IF

        CALL field%get_ptr(ptr, igrid)
    END SUBROUTINE get_fieldptr_grid3


    SUBROUTINE get_fieldptr(ptr, name)
        ! Subroutine arguments
        REAL(realk), INTENT(out), POINTER, CONTIGUOUS :: ptr(:)
        CHARACTER(len=*), INTENT(in) :: name

        ! Local variables
        TYPE(field_t), POINTER :: field
        LOGICAL :: found

        CALL get_field(field, name, found)

        IF (.NOT. found) THEN
            WRITE(*, '("Field ", A, " does not exist!")') name
            CALL errr(__FILE__, __LINE__)
        END IF

        ptr(1:field%idim) => field%arr
    END SUBROUTINE get_fieldptr


    SUBROUTINE get_fileh(fileh)
        INTEGER(hid_t), INTENT(out) :: fileh

        IF (file_id == 0) CALL errr(__FILE__, __LINE__)
        fileh = file_id
    END SUBROUTINE get_fileh

END MODULE fields_mod
