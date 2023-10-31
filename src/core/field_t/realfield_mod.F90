MODULE realfield_mod
    USE err_mod, ONLY: errr
    USE grids_mod, ONLY: get_mgdims, mygrids, nmygrids, level
    USE pointers_mod, ONLY: idim2d, get_ibb
    USE precision_mod, ONLY: intk, realk, mglet_hdf5_real, mglet_mpi_real
    USE utils_mod, ONLY: get_stag_shift
    USE basefield_mod
    USE commbuf_mod, ONLY: bigbuf

    IMPLICIT NONE(type, external)
    PRIVATE

    TYPE :: buffer_t
        LOGICAL :: is_init = .FALSE.
        REAL(realk), ALLOCATABLE :: fr(:)
        REAL(realk), ALLOCATABLE :: ba(:)

        REAL(realk), ALLOCATABLE :: ri(:)
        REAL(realk), ALLOCATABLE :: le(:)

        REAL(realk), ALLOCATABLE :: bo(:)
        REAL(realk), ALLOCATABLE :: to(:)
    CONTAINS
        FINAL :: buffer_destructor
        PROCEDURE :: finish => finish_buffer
        PROCEDURE :: get_buffer
        PROCEDURE :: init => init_buffer
        PROCEDURE, PRIVATE :: copy_buffer
        GENERIC :: ASSIGNMENT(=) => copy_buffer
    END TYPE buffer_t

    TYPE, EXTENDS(basefield_t) :: field_t
        REAL(realk), ALLOCATABLE :: arr(:)
        TYPE(buffer_t) :: buffers
    CONTAINS
        PROCEDURE :: init

        GENERIC, PUBLIC :: get_ptr => get_grid1, get_grid3
        PROCEDURE, PRIVATE :: get_grid1, get_grid3

        PROCEDURE :: copy_from
        PROCEDURE :: multiply
        PROCEDURE :: shift
        PROCEDURE :: get_buffers
        PROCEDURE :: init_buffers
        PROCEDURE :: finish
        FINAL :: destructor

        PROCEDURE, PRIVATE :: set_const
        GENERIC :: ASSIGNMENT(=) => set_const
    END TYPE field_t

    PUBLIC :: field_t, buffer_t

CONTAINS
    SUBROUTINE init(this, name, description, ndim, istag, jstag, kstag, &
            units, dread, required, dwrite, active_level, get_len)
        ! Subroutine arguments
        CLASS(field_t), INTENT(out) :: this
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
        ! none...

        CALL this%init_corefield(name, description, ndim, istag, jstag, kstag, &
            units, dread, required, dwrite, active_level, get_len)

        this%hdf5_dtype = mglet_hdf5_real
        this%mpi_dtype = mglet_mpi_real

        ALLOCATE(this%arr(this%idim))
        this%arr = 0.0
    END SUBROUTINE init


    ELEMENTAL SUBROUTINE finish(this)
        CLASS(field_t), INTENT(inout) :: this

        IF (.NOT. this%is_init) RETURN
        CALL this%finish_corefield()

        DEALLOCATE(this%arr)
        IF (this%buffers%is_init) CALL this%buffers%finish()
    END SUBROUTINE finish


    ELEMENTAL SUBROUTINE destructor(this)
        TYPE(field_t), INTENT(inout) :: this

        CALL this%finish()
    END SUBROUTINE destructor


    SUBROUTINE get_grid1(this, ptr, igrid)
        ! Subroutine arguments
        CLASS(field_t), INTENT(in), TARGET :: this
        REAL(realk), POINTER, CONTIGUOUS, INTENT(out) :: ptr(:)
        INTEGER(intk), INTENT(in) :: igrid

        ! Local variables
        INTEGER(intk) :: ip, len

        CALL this%get_ip(ip, igrid)
        CALL this%get_len(len, igrid)
        IF (len <= 0) CALL errr(__FILE__, __LINE__)

        ptr(1:len) => this%arr(ip:ip+len-1)
    END SUBROUTINE get_grid1


    SUBROUTINE get_grid3(this, ptr, igrid)
        ! Subroutine arguments
        CLASS(field_t), INTENT(in), TARGET :: this
        REAL(realk), POINTER, CONTIGUOUS, INTENT(out) :: ptr(:, :, :)
        INTEGER(intk), INTENT(in) :: igrid

        ! Local variables
        INTEGER(intk) :: kk, jj, ii, ip, len

        IF (.NOT. this%ndim == 3) THEN
            WRITE(*, '("Field ", A, " is not 3D!")') TRIM(this%name)
            CALL errr(__FILE__, __LINE__)
        END IF

        CALL this%get_ip(ip, igrid)
        CALL this%get_len(len, igrid)
        IF (len <= 0) CALL errr(__FILE__, __LINE__)

        CALL get_mgdims(kk, jj, ii, igrid)
        IF (len /= kk*jj*ii) CALL errr(__FILE__, __LINE__)

        ptr(1:kk, 1:jj, 1:ii) => this%arr(ip:ip+kk*jj*ii-1)
    END SUBROUTINE get_grid3


    ! Make a deep copy of another field
    SUBROUTINE copy_from(this, that)
        ! Subroutine arguments
        CLASS(field_t), INTENT(inout) :: this
        CLASS(field_t), INTENT(in) :: that

        this%is_init = that%is_init
        this%name = that%name
        this%description = that%description
        this%ndim = that%ndim
        this%istag = that%istag
        this%jstag = that%jstag
        this%kstag = that%kstag
        this%units = that%units
        this%dread = that%dread
        this%required = that%required
        this%dwrite = that%dwrite

        IF (.NOT. ALLOCATED(this%active_level)) THEN
            ALLOCATE(this%active_level, mold=that%active_level)
        END IF
        this%active_level = that%active_level

        this%get_len => that%get_len

        IF (.NOT. ALLOCATED(this%ptr)) THEN
            ALLOCATE(this%ptr, mold=that%ptr)
        END IF
        this%ptr = that%ptr

        this%idim = that%idim

        IF (ALLOCATED(this%arr)) THEN
            IF (SIZE(this%arr) /= SIZE(that%arr)) THEN
                DEALLOCATE(this%arr)
            END IF
        END IF
        IF (.NOT. ALLOCATED(this%arr)) THEN
            ALLOCATE(this%arr, mold=that%arr)
        END IF
        this%arr = that%arr

        ! TODO: buffers???

        this%n_iattr = that%n_iattr
        this%iattr_key = that%iattr_key
        this%iattr_val = that%iattr_val

        this%n_rattr = that%n_rattr
        this%rattr_key = that%rattr_key
        this%rattr_val = that%rattr_val
    END SUBROUTINE copy_from


    ! Multiply two fields, takes staggering of variables into account. The
    ! destination field istag, jstag, kstag determine the interpolation of the
    ! source fields
    !
    ! Prior content in the field is discarded, but the staggering of the
    ! destination determine the interpolation of the source fields.
    SUBROUTINE multiply(this, a, b)
        ! Subroutine arguments
        CLASS(field_t), INTENT(inout) :: this
        CLASS(field_t), INTENT(in) :: a
        CLASS(field_t), INTENT(in) :: b

        ! Local variables
        INTEGER(intk) :: igr, ilevel, igrid
        INTEGER(intk) :: kk, jj, ii
        INTEGER(intk) :: k, j, i
        INTEGER(intk) :: kstart, jstart, istart
        INTEGER(intk) :: kstop, jstop, istop
        INTEGER(intk) :: k1, j1, i1
        INTEGER(intk) :: k2, j2, i2
        REAL(realk), POINTER, CONTIGUOUS :: out(:, :, :), phi1(:, :, :), &
            phi2(:, :, :)

        DO igr = 1, nmygrids
            igrid = mygrids(igr)
            ilevel = level(igrid)

            ! All fields must be defiend on the same levels
            IF (.NOT. this%active_level(ilevel)) CYCLE
            IF (a%active_level(ilevel) .EQV. .FALSE. .OR. &
                    b%active_level(ilevel) .EQV. .FALSE.) THEN
                CALL errr(__FILE__, __LINE__)
            END IF

            CALL get_mgdims(kk, jj, ii, igrid)

            istart = 1
            istop = ii
            CALL get_stag_shift(i1, istart, istop, ii, this%istag, a%istag)
            CALL get_stag_shift(i2, istart, istop, ii, this%istag, b%istag)

            jstart = 1
            jstop = jj
            CALL get_stag_shift(j1, jstart, jstop, jj, this%jstag, a%jstag)
            CALL get_stag_shift(j2, jstart, jstop, jj, this%jstag, b%jstag)

            kstart = 1
            kstop = kk
            CALL get_stag_shift(k1, kstart, kstop, kk, this%kstag, a%kstag)
            CALL get_stag_shift(k2, kstart, kstop, kk, this%kstag, b%kstag)

            CALL this%get_ptr(out, igrid)
            CALL a%get_ptr(phi1, igrid)
            CALL b%get_ptr(phi2, igrid)

            DO i = istart, istop
                DO j = jstart, jstop
                    DO k = kstart, kstop
                        out(k, j, i) = &
                            0.25*(phi1(k, j, i) + phi1(k+k1, j+j1, i+i1)) &
                            *(phi2(k, j, i) + phi2(k+k2, j+j2, i+i2))
                    END DO
                END DO
            END DO
        END DO
    END SUBROUTINE multiply


    ! Shift position (staggering) of field. The destination field istag, jstag,
    ! kstag determine the interpolation of the source fields
    !
    ! Prior content in the field is discarded, but the staggering of the
    ! destination determine the interpolation of the source fields.
    SUBROUTINE shift(this, that)
        ! Subroutine arguments
        CLASS(field_t), INTENT(inout) :: this
        CLASS(field_t), INTENT(in) :: that

        ! Local variables
        INTEGER(intk) :: igr, ilevel, igrid
        INTEGER(intk) :: kk, jj, ii
        INTEGER(intk) :: k, j, i
        INTEGER(intk) :: kstart, jstart, istart
        INTEGER(intk) :: kstop, jstop, istop
        INTEGER(intk) :: k1, j1, i1
        REAL(realk), POINTER, CONTIGUOUS :: out(:, :, :), in(:, :, :)

        DO igr = 1, nmygrids
            igrid = mygrids(igr)
            ilevel = level(igrid)

            ! All fields must be defiend on the same levels
            IF (.NOT. this%active_level(ilevel)) CYCLE
            IF (that%active_level(ilevel) .EQV. .FALSE.) THEN
                CALL errr(__FILE__, __LINE__)
            END IF

            CALL get_mgdims(kk, jj, ii, igrid)
            CALL this%get_ptr(out, igrid)
            CALL that%get_ptr(in, igrid)

            istart = 1
            istop = ii
            CALL get_stag_shift(i1, istart, istop, ii, this%istag, that%istag)

            jstart = 1
            jstop = jj
            CALL get_stag_shift(j1, jstart, jstop, jj, this%jstag, that%jstag)

            kstart = 1
            kstop = kk
            CALL get_stag_shift(k1, kstart, kstop, kk, this%kstag, that%kstag)

            DO i = istart, istop
                DO j = jstart, jstop
                    DO k = kstart, kstop
                        out(k, j, i) = 0.125*(in(k, j, i) + in(k+k1, j, i) &
                            + in(k, j+j1, i) + in(k, j, i+i1) &
                            + in(k+k1, j+j1, i) + in(k+k1, j, i+i1) &
                            + in(k, j+j1, i+i1) + in(k+k1, j+j1, i+i1))
                    END DO
                END DO
            END DO
        END DO
    END SUBROUTINE shift


    SUBROUTINE init_buffers(this)
        ! Subroutine arguments
        CLASS(field_t), TARGET, INTENT(inout) :: this

        IF (.NOT. this%buffers%is_init) THEN
            CALL this%buffers%init()
        END IF
    END SUBROUTINE init_buffers


    SUBROUTINE get_buffers(this, buffers)
        ! Allocate and (optionally) fetch buffers

        ! Subroutine arguments
        CLASS(field_t), TARGET, INTENT(inout) :: this
        TYPE(buffer_t), POINTER, INTENT(inout) :: buffers

        IF (.NOT. this%buffers%is_init) THEN
            CALL errr(__FILE__, __LINE__)
        END IF
        buffers => this%buffers
    END SUBROUTINE get_buffers


    SUBROUTINE set_const(this, val)
        ! Set a field to a constant scalar value

        ! Subroutine arguments
        CLASS(field_t), INTENT(inout) :: this
        REAL(realk), INTENT(in) :: val

        this%arr = val
    END SUBROUTINE set_const


    SUBROUTINE init_buffer(this)
        ! Subroutine arguments
        CLASS(buffer_t), INTENT(inout) :: this

        ! Local variables
        ! none...

        IF (this%is_init) CALL errr(__FILE__, __LINE__)

        ALLOCATE(this%fr(2*idim2d))
        ALLOCATE(this%ba(2*idim2d))
        ALLOCATE(this%ri(2*idim2d))
        ALLOCATE(this%le(2*idim2d))
        ALLOCATE(this%bo(2*idim2d))
        ALLOCATE(this%to(2*idim2d))

        this%fr = 0.0
        this%ba = 0.0
        this%ri = 0.0
        this%le = 0.0
        this%bo = 0.0
        this%to = 0.0

        this%is_init = .TRUE.
    END SUBROUTINE init_buffer


    ELEMENTAL SUBROUTINE buffer_destructor(this)
        TYPE(buffer_t), INTENT(inout) :: this
        CALL this%finish()
    END SUBROUTINE buffer_destructor


    PURE SUBROUTINE finish_buffer(this)
        CLASS(buffer_t), INTENT(inout) :: this
        this%is_init = .FALSE.
        IF (ALLOCATED(this%fr)) DEALLOCATE(this%fr)
        IF (ALLOCATED(this%ba)) DEALLOCATE(this%ba)
        IF (ALLOCATED(this%ri)) DEALLOCATE(this%ri)
        IF (ALLOCATED(this%le)) DEALLOCATE(this%le)
        IF (ALLOCATED(this%bo)) DEALLOCATE(this%bo)
        IF (ALLOCATED(this%to)) DEALLOCATE(this%to)
    END SUBROUTINE finish_buffer


    SUBROUTINE get_buffer(this, ptr, igrid, iface)
        ! Subroutine arguments
        CLASS(buffer_t), INTENT(inout), TARGET :: this
        REAL(realk), INTENT(out), POINTER, CONTIGUOUS :: ptr(:, :, :)
        INTEGER(intk), INTENT(in) :: igrid
        INTEGER(intk), INTENT(in) :: iface

        ! Local variables
        INTEGER(intk) :: kk, jj, ii, ibb

        IF (.NOT. this%is_init) THEN
            WRITE(*,*) "Buffers not initialized"
            CALL errr(__FILE__, __LINE__)
        END IF

        CALL get_mgdims(kk, jj, ii, igrid)
        CALL get_ibb(ibb, igrid)

        SELECT CASE (iface)
        CASE (1)
            ptr(1:kk, 1:jj, 1:2) => this%fr(ibb:ibb+2*kk*jj-1)
        CASE (2)
            ptr(1:kk, 1:jj, 1:2) => this%ba(ibb:ibb+2*kk*jj-1)
        CASE (3)
            ptr(1:kk, 1:ii, 1:2) => this%ri(ibb:ibb+2*kk*ii-1)
        CASE (4)
            ptr(1:kk, 1:ii, 1:2) => this%le(ibb:ibb+2*kk*ii-1)
        CASE (5)
            ptr(1:jj, 1:ii, 1:2) => this%bo(ibb:ibb+2*jj*ii-1)
        CASE (6)
            ptr(1:jj, 1:ii, 1:2) => this%to(ibb:ibb+2*jj*ii-1)
        CASE DEFAULT
            WRITE(*, '("Invalid face: ", I0)') iface
            CALL errr(__FILE__, __LINE__)
        END SELECT
    END SUBROUTINE get_buffer


    SUBROUTINE copy_buffer(this, that)
        ! Copy buffers (from one field to antother)
        ! Example usage:
        !     pwu%buffers = u%buffers
        !     pwv%buffers = v%buffers
        !     pww%buffers = w%buffers
        ! from setpointvalues

        ! Subroutine arguments
        CLASS(buffer_t), INTENT(inout) :: this
        CLASS(buffer_t), INTENT(in) :: that

        this%fr(:) = that%fr(:)
        this%ba(:) = that%ba(:)
        this%ri(:) = that%ri(:)
        this%le(:) = that%le(:)
        this%bo(:) = that%bo(:)
        this%to(:) = that%to(:)
    END SUBROUTINE copy_buffer

END MODULE realfield_mod
