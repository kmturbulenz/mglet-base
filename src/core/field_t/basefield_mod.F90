MODULE basefield_mod
    USE HDF5
    USE MPI_f08
    USE err_mod, ONLY: errr
    USE grids_mod, ONLY: mygrids, nmygrids, minlevel, maxlevel, get_imygrid, &
        level
    USE pointers_mod, ONLY: get_len3, idim2d, get_ibb
    USE precision_mod, ONLY: intk, realk


    IMPLICIT NONE(type, external)
    PRIVATE

    INTEGER(intk), PARAMETER :: nchar_name = 16
    INTEGER(intk), PARAMETER :: nchar_desc = 32
    INTEGER(intk), PARAMETER :: nattr_max = 8

    TYPE, ABSTRACT :: basefield_t
        LOGICAL :: is_init = .FALSE.

        CHARACTER(len=nchar_name) :: name = REPEAT(" ", nchar_name)
        CHARACTER(len=nchar_desc) :: description = REPEAT(" ", nchar_desc)
        INTEGER(intk) :: ndim = 3

        INTEGER(intk) :: istag = 0
        INTEGER(intk) :: jstag = 0
        INTEGER(intk) :: kstag = 0
        INTEGER(intk) :: units(7) = [0, 0, 0, 0, 0, 0, 0]

        ! Try to read field
        LOGICAL :: dread = .FALSE.

        ! Requierd field: if dread = .TRUE. and not found, kill
        LOGICAL :: required = .TRUE.

        ! Write field
        LOGICAL :: dwrite = .FALSE.

        ! Allocated to minlevel:maxlevel and indicate which levels the field
        ! is active on
        LOGICAL, ALLOCATABLE :: active_level(:)

        ! Pointers to procedures for getting the pointer and length
        PROCEDURE(get_len_i), POINTER, NOPASS :: get_len => NULL()
        INTEGER(intk), ALLOCATABLE :: ptr(:)

        ! Length and actual data array can be allocated to anything
        INTEGER(intk) :: idim = 0

        ! HDF5 datatype for IO operations
        INTEGER(hid_t) :: hdf5_dtype

        ! MPI datatype for communication
        TYPE(MPI_Datatype) :: mpi_dtype

        ! Attributes
        INTEGER(intk) :: n_iattr = 0
        CHARACTER(len=nchar_name) :: iattr_key(nattr_max) = &
            REPEAT(" ", nchar_name)
        INTEGER(intk) :: iattr_val(nattr_max) = 0

        INTEGER(intk) :: n_rattr = 0
        CHARACTER(len=nchar_name) :: rattr_key(nattr_max) = &
            REPEAT(" ", nchar_name)
        REAL(realk) :: rattr_val(nattr_max) = 0.0
    CONTAINS
        PROCEDURE :: init_corefield

        GENERIC, PUBLIC :: get_attr => get_rattr, get_iattr
        PROCEDURE, PRIVATE :: get_rattr, get_iattr

        GENERIC, PUBLIC :: set_attr => set_rattr, set_iattr
        PROCEDURE, PRIVATE :: set_rattr, set_iattr

        PROCEDURE :: get_ip

        PROCEDURE :: finish_corefield
    END TYPE basefield_t

    ABSTRACT INTERFACE
        SUBROUTINE get_len_i(len, igrid)
            IMPORT :: intk
            INTEGER(intk), INTENT(out) :: len
            INTEGER(intk), INTENT(in) :: igrid
        END SUBROUTINE get_len_i
    END INTERFACE

    PUBLIC :: basefield_t, get_len_i, nchar_name, nchar_desc, nattr_max

CONTAINS
    SUBROUTINE init_corefield(this, name, description, ndim, istag, &
            jstag, kstag, units, dread, required, dwrite, active_level, get_len)
        ! Subroutine arguments
        CLASS(basefield_t), INTENT(out) :: this
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
        LOGICAL, INTENT(in), OPTIONAL :: active_level(:)
        PROCEDURE(get_len_i), OPTIONAL :: get_len

        ! Local variables
        INTEGER(intk) :: i, igrid, length

        this%is_init = .TRUE.

        IF (LEN_TRIM(name) > nchar_name) CALL errr(__FILE__, __LINE__)
        this%name = name

        this%description = name
        IF (PRESENT(description)) THEN
            IF (LEN_TRIM(description) > nchar_desc) THEN
                CALL errr(__FILE__, __LINE__)
            END IF
            this%description = description
        END IF

        IF (PRESENT(ndim)) THEN
            this%ndim = ndim
        END IF

        IF (PRESENT(istag)) THEN
            this%istag = istag
        END IF

        IF (PRESENT(jstag)) THEN
            this%jstag = jstag
        END IF

        IF (PRESENT(kstag)) THEN
            this%kstag = kstag
        END IF

        IF (PRESENT(units)) THEN
            this%units = units
        END IF

        IF (PRESENT(dread)) THEN
            this%dread = dread
        END IF

        IF (PRESENT(required)) THEN
            this%required = required
        END IF

        IF (PRESENT(dwrite)) THEN
            this%dwrite = dwrite
        END IF

        ALLOCATE(this%active_level(minlevel:maxlevel))
        this%active_level = .TRUE.
        IF (PRESENT(active_level)) THEN
            IF (SIZE(active_level) /= maxlevel-minlevel+1) THEN
                CALL errr(__FILE__, __LINE__)
            END IF
            this%active_level(:) = active_level(:)
        END IF

        ! Only sets default pointers for 3D fields
        IF (this%ndim == 3) THEN
            this%get_len => get_len3
        END IF
        IF (PRESENT(get_len)) THEN
            this%get_len => get_len
        END IF
        IF (.NOT. ASSOCIATED(this%get_len)) CALL errr(__FILE__, __LINE__)

        ! Set pointers
        ALLOCATE(this%ptr(nmygrids))
        this%ptr = 0

        this%idim = 0
        DO i = 1, nmygrids
            igrid = mygrids(i)
            IF (.NOT. this%active_level(level(igrid))) CYCLE

            CALL this%get_len(length, igrid)
            this%ptr(i) = this%idim + 1
            this%idim = this%idim + length
        END DO
    END SUBROUTINE init_corefield


    ELEMENTAL SUBROUTINE finish_corefield(this)
        CLASS(basefield_t), INTENT(inout) :: this

        INTEGER(intk) :: i

        IF (.NOT. this%is_init) RETURN
        this%is_init = .FALSE.

        this%name = REPEAT(" ", nchar_name)
        this%description = REPEAT(" ", nchar_desc)
        this%idim = 0
        this%ndim = 3
        this%istag = 0
        this%jstag = 0
        this%kstag = 0
        this%units = 0
        this%dread = .FALSE.
        this%required = .FALSE.
        this%dwrite = .FALSE.

        DO i = 1, nattr_max
            this%iattr_key = REPEAT(" ", nchar_name)
            this%iattr_val = 0
        END DO
        this%n_iattr = 0

        DO i = 1, nattr_max
            this%rattr_key = REPEAT(" ", nchar_name)
            this%rattr_val = 0.0
        END DO
        this%n_rattr = 0

        DEALLOCATE(this%active_level)
        DEALLOCATE(this%ptr)
    END SUBROUTINE finish_corefield


    SUBROUTINE get_ip(this, ip, igrid)
        CLASS(basefield_t), INTENT(in) :: this
        INTEGER(intk), INTENT(out) :: ip
        INTEGER(intk), INTENT(in) :: igrid

        INTEGER(intk) :: imygrid
        ip = 0
        CALL get_imygrid(imygrid, igrid)
        ip = this%ptr(imygrid)
    END SUBROUTINE get_ip


    SUBROUTINE get_rattr(this, val, key)
        ! Subroutine arguments
        CLASS(basefield_t), INTENT(inout) :: this
        REAL(realk), INTENT(inout) :: val
        CHARACTER(len=*), INTENT(in) :: key

        ! Local variables
        INTEGER(intk) :: i
        LOGICAL :: thisfound

        IF (LEN_TRIM(key) > nchar_name) CALL errr(__FILE__, __LINE__)

        thisfound = .FALSE.
        DO i = 1, this%n_rattr
            IF (TRIM(this%rattr_key(i)) == TRIM(key)) THEN
                val = this%rattr_val(i)
                thisfound = .TRUE.
                EXIT
            END IF
        END DO

        IF (.NOT. thisfound) THEN
            WRITE(*, '("Attribute ", A, " does not exist!")') key
            CALL errr(__FILE__, __LINE__)
        END IF
    END SUBROUTINE get_rattr


    SUBROUTINE get_iattr(this, val, key)
        ! Subroutine arguments
        CLASS(basefield_t), INTENT(inout) :: this
        INTEGER(intk), INTENT(inout) :: val
        CHARACTER(len=*), INTENT(in) :: key

        ! Local variables
        INTEGER(intk) :: i
        LOGICAL :: thisfound

        IF (LEN_TRIM(key) > nchar_name) CALL errr(__FILE__, __LINE__)

        thisfound = .FALSE.
        DO i = 1, this%n_iattr
            IF (TRIM(this%iattr_key(i)) == TRIM(key)) THEN
                val = this%iattr_val(i)
                thisfound = .TRUE.
                EXIT
            END IF
        END DO

        IF (.NOT. thisfound) THEN
            WRITE(*, '("Attribute ", A, " does not exist!")') key
            CALL errr(__FILE__, __LINE__)
        END IF
    END SUBROUTINE get_iattr


    SUBROUTINE set_rattr(this, val, key)
        ! Subroutine arguments
        CLASS(basefield_t), INTENT(inout) :: this
        REAL(realk), INTENT(in) :: val
        CHARACTER(len=*), INTENT(in) :: key

        ! Local variables
        INTEGER(intk) :: i

        IF (LEN_TRIM(key) > nchar_name) CALL errr(__FILE__, __LINE__)

        ! First check if it exists already
        DO i = 1, this%n_rattr
            IF (TRIM(this%rattr_key(i)) == TRIM(key)) THEN
                this%rattr_val(i) = val
                RETURN
            END IF
        END DO

        ! Not existing, is there space for another attribute?
        IF (this%n_rattr + 1 > nattr_max) THEN
            WRITE(*, '("Field ", A, " no space for attributes")') this%name
            CALL errr(__FILE__, __LINE__)
        END IF

        ! Insert attribute
        this%n_rattr = this%n_rattr + 1
        this%rattr_key(this%n_rattr) = TRIM(key)
        this%rattr_val(this%n_rattr) = val
    END SUBROUTINE set_rattr


    SUBROUTINE set_iattr(this, val, key)
        ! Subroutine arguments
        CLASS(basefield_t), INTENT(inout) :: this
        INTEGER(intk), INTENT(in) :: val
        CHARACTER(len=*), INTENT(in) :: key

        ! Local variables
        INTEGER(intk) :: i

        IF (LEN_TRIM(key) > nchar_name) CALL errr(__FILE__, __LINE__)

        ! First check if it exists already
        DO i = 1, this%n_iattr
            IF (TRIM(this%iattr_key(i)) == TRIM(key)) THEN
                this%iattr_val(i) = val
                RETURN
            END IF
        END DO

        ! Not existing, is there space for another attribute?
        IF (this%n_rattr + 1 > nattr_max) THEN
            WRITE(*, '("Field ", A, " no space for attributes")') this%name
            CALL errr(__FILE__, __LINE__)
        END IF

        ! Insert attribute
        this%n_iattr = this%n_iattr + 1
        this%iattr_key(this%n_iattr) = TRIM(key)
        this%iattr_val(this%n_iattr) = val
    END SUBROUTINE set_iattr
END MODULE basefield_mod
