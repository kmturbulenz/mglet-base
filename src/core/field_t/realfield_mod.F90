MODULE realfield_mod
    USE err_mod, ONLY: errr
    USE grids_mod, ONLY: get_mgdims, mygrids, nmygrids, level, get_imygrid, globalgrids
    USE pointers_mod, ONLY: idimbb, get_ibb
    USE precision_mod, ONLY: intk, realk, mglet_hdf5_real, mglet_mpi_real
    USE utils_mod, ONLY: get_stag_shift
    USE basefield_mod

    IMPLICIT NONE(type, external)
    PRIVATE

    TYPE, EXTENDS(basefield_t) :: field_t
        REAL(realk), ALLOCATABLE :: arr(:)
        REAL(realk), ALLOCATABLE :: buffers(:)
    CONTAINS
        PROCEDURE :: init

        GENERIC, PUBLIC :: get_ptr => get_grid1, get_grid3
        GENERIC, PUBLIC :: multiply => multiply2, multiply3

        PROCEDURE, PRIVATE, NON_OVERRIDABLE :: get_grid1, get_grid3
        PROCEDURE, PRIVATE, NON_OVERRIDABLE :: multiply2, multiply3

        PROCEDURE, NON_OVERRIDABLE :: get_value
        PROCEDURE, NON_OVERRIDABLE :: get_buffer
        PROCEDURE :: copy_from
        PROCEDURE :: shift
        PROCEDURE :: init_buffers
        PROCEDURE :: finish
        FINAL :: destructor
    END TYPE field_t

    PUBLIC :: field_t

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
        IF (ALLOCATED(this%buffers)) THEN
            DEALLOCATE(this%buffers)
        END IF
    END SUBROUTINE finish


    ELEMENTAL SUBROUTINE destructor(this)
        TYPE(field_t), INTENT(inout) :: this

        CALL this%finish()
    END SUBROUTINE destructor


    ! SUBROUTINE get_grid1(this, ptr, igrid)
    ! $omp declare target

        ! Subroutine arguments
        ! CLASS(field_t), INTENT(in), TARGET :: this
        ! REAL(realk), POINTER, CONTIGUOUS, INTENT(out) :: ptr(:)
        ! INTEGER(intk), INTENT(in) :: igrid

        ! Local variables
        ! INTEGER(intk) :: ip, len

        ! WARNUNG: Hier Roituinen des Basistypes
        ! CALL this%get_ip(ip, igrid)
        ! CALL this%get_len(len, igrid)
        ! IF (len <= 0) CALL errr(__FILE__, __LINE__)

        ! ptr(1:len) => this%arr(ip:ip+len-1)

    ! END SUBROUTINE get_grid1


    SUBROUTINE get_grid1(this, ptr, igrid, lin, info)
    !$omp declare target

        ! Subroutine arguments
        CLASS(field_t), INTENT(in), TARGET :: this
        REAL(realk), POINTER, CONTIGUOUS, INTENT(out) :: ptr(:)
        INTEGER(intk), INTENT(in) :: igrid
        LOGICAL, OPTIONAL :: lin
        LOGICAL, OPTIONAL :: info

        ! Local variables
        INTEGER(intk) :: i, len, ip, ii, jj, kk
        LOGICAL :: linearize

        ! Setting for linearization
        linearize = .FALSE.
        IF (PRESENT(lin)) THEN
            IF (lin) THEN
                linearize = .TRUE.
            END IF
        END IF

        IF (.NOT. this%ndim == 1 .AND. .NOT. linearize) THEN
            WRITE(*, '("Field is not 1D!")')
            ! CALL errr(__FILE__, __LINE__)
        END IF

        IF (.NOT. this%ndim == 3 .AND. linearize) THEN
            WRITE(*, '("Field is not 3D! No linearization")')
            ! CALL errr(__FILE__, __LINE__)
        END IF


        i = globalgrids(igrid)
        ip = this%ptr(i)
        len = this%length(i)

        IF (PRESENT(info)) THEN
            IF (info) THEN
               WRITE(*, *) "get_grid1", igrid, i, ip, len, linearize
            END IF
        END IF

        IF (.NOT. linearize) THEN
            IF (len <= 0) WRITE(*, *) "ERROR"
            ptr(1:len) => this%arr(ip:ip+len-1)
        ELSE
            CALL get_mgdims(kk, jj, ii, igrid)
            IF (len /= kk*jj*ii) WRITE(*, *) "ERROR"
            ptr(1:kk*jj*ii) => this%arr(ip:ip+kk*jj*ii-1)
        END IF

    END SUBROUTINE get_grid1




    SUBROUTINE get_grid3(this, ptr, igrid, lin, info)
    !$omp declare target

        ! Subroutine arguments
        CLASS(field_t), INTENT(in), TARGET :: this
        REAL(realk), POINTER, CONTIGUOUS, INTENT(out) :: ptr(:, :, :)
        INTEGER(intk), INTENT(in) :: igrid
        LOGICAL, OPTIONAL :: lin
        LOGICAL, OPTIONAL :: info

        ! Local variables
        INTEGER(intk) :: kk, jj, ii, ip, len, i
        LOGICAL :: linearize

        IF (PRESENT(lin)) THEN
            IF (lin) THEN
                WRITE(*, *) "ERROR: No linearization for 3D fields"
            END IF
        END IF

        IF (.NOT. this%ndim == 3) THEN
            WRITE(*, '("Field ", A, " is not 3D!")') TRIM(this%name)
            ! CALL errr(__FILE__, __LINE__)
        END IF

        CALL get_imygrid(i, igrid)
        ip = this%ptr(i)
        len = this%length(i)

        CALL get_mgdims(kk, jj, ii, igrid)
        IF (len /= kk*jj*ii) WRITE(*, *) "Error"
        IF (PRESENT(info)) THEN
            IF (info) THEN
                WRITE(*, *) "get_grid3", igrid, i, ip, len, linearize
            END IF
        END IF

        ptr(1:kk, 1:jj, 1:ii) => this%arr(ip:ip+kk*jj*ii-1)
    END SUBROUTINE get_grid3


    REAL(realk) FUNCTION get_value(this, k, j, i, igrid)
        ! Function arguments
        CLASS(field_t), INTENT(in) :: this
        INTEGER(intk), INTENT(in) :: k
        INTEGER(intk), INTENT(in) :: j
        INTEGER(intk), INTENT(in) :: i
        INTEGER(intk), INTENT(in) :: igrid

        ! Local variables
        INTEGER(intk) :: ip
        INTEGER(intk) :: kk, jj, ii

        IF (.NOT. this%ndim == 3) THEN
            WRITE(*, '("Field ", A, " is not 3D!")') TRIM(this%name)
            CALL errr(__FILE__, __LINE__)
        END IF

        CALL this%get_ip(ip, igrid)
        CALL get_mgdims(kk, jj, ii, igrid)

        IF (k < 1 .OR. k > kk) THEN
            WRITE(*, *) "k out of bounds:", k, kk
            CALL errr(__FILE__, __LINE__)
        END IF

        IF (j < 1 .OR. j > jj) THEN
            WRITE(*, *) "j out of bounds:", j, jj
            CALL errr(__FILE__, __LINE__)
        END IF

        IF (i < 1 .OR. i > ii) THEN
            WRITE(*, *) "i out of bounds:", i, ii
            CALL errr(__FILE__, __LINE__)
        END IF

        get_value = this%arr(ip + (k-1) + (j-1)*kk + (i-1)*kk*jj)
    END FUNCTION get_value


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

        IF (.NOT. ALLOCATED(this%length)) THEN
            ALLOCATE(this%length, mold=that%length)
        END IF
        this%length = that%length

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
    SUBROUTINE multiply2(this, a, b)
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
    END SUBROUTINE multiply2


    SUBROUTINE multiply3(this, a, b, c)
        ! Subroutine arguments
        CLASS(field_t), INTENT(inout) :: this
        CLASS(field_t), INTENT(in) :: a
        CLASS(field_t), INTENT(in) :: b
        CLASS(field_t), INTENT(in) :: c

        ! Local variables
        INTEGER(intk) :: igr, ilevel, igrid
        INTEGER(intk) :: kk, jj, ii
        INTEGER(intk) :: k, j, i
        INTEGER(intk) :: kstart, jstart, istart
        INTEGER(intk) :: kstop, jstop, istop
        INTEGER(intk) :: k1, j1, i1
        INTEGER(intk) :: k2, j2, i2
        INTEGER(intk) :: k3, j3, i3
        REAL(realk), POINTER, CONTIGUOUS :: out(:, :, :), phi1(:, :, :), &
            phi2(:, :, :), phi3(:, :, :)

        DO igr = 1, nmygrids
            igrid = mygrids(igr)
            ilevel = level(igrid)

            ! All fields must be defiend on the same levels
            IF (.NOT. this%active_level(ilevel)) CYCLE
            IF (a%active_level(ilevel) .EQV. .FALSE. .OR. &
                    b%active_level(ilevel) .EQV. .FALSE. .OR. &
                    c%active_level(ilevel) .EQV. .FALSE.) THEN
                CALL errr(__FILE__, __LINE__)
            END IF

            CALL get_mgdims(kk, jj, ii, igrid)

            istart = 1
            istop = ii
            CALL get_stag_shift(i1, istart, istop, ii, this%istag, a%istag)
            CALL get_stag_shift(i2, istart, istop, ii, this%istag, b%istag)
            CALL get_stag_shift(i3, istart, istop, ii, this%istag, c%istag)

            jstart = 1
            jstop = jj
            CALL get_stag_shift(j1, jstart, jstop, jj, this%jstag, a%jstag)
            CALL get_stag_shift(j2, jstart, jstop, jj, this%jstag, b%jstag)
            CALL get_stag_shift(j3, jstart, jstop, jj, this%jstag, c%jstag)

            kstart = 1
            kstop = kk
            CALL get_stag_shift(k1, kstart, kstop, kk, this%kstag, a%kstag)
            CALL get_stag_shift(k2, kstart, kstop, kk, this%kstag, b%kstag)
            CALL get_stag_shift(k3, kstart, kstop, kk, this%kstag, c%kstag)

            CALL this%get_ptr(out, igrid)
            CALL a%get_ptr(phi1, igrid)
            CALL b%get_ptr(phi2, igrid)
            CALL c%get_ptr(phi3, igrid)

            DO i = istart, istop
                DO j = jstart, jstop
                    DO k = kstart, kstop
                        out(k, j, i) = &
                            0.125*(phi1(k, j, i) + phi1(k+k1, j+j1, i+i1)) &
                            *(phi2(k, j, i) + phi2(k+k2, j+j2, i+i2)) &
                            *(phi3(k, j, i) + phi3(k+k3, j+j3, i+i3))
                    END DO
                END DO
            END DO
        END DO
    END SUBROUTINE multiply3

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

        IF (ALLOCATED(this%buffers)) THEN
            RETURN
        END IF
        ALLOCATE(this%buffers(idimbb))

#ifdef _MGLET_DEBUG_
        BLOCK
            USE, INTRINSIC :: IEEE_ARITHMETIC
            USE, INTRINSIC :: IEEE_EXCEPTIONS

            LOGICAL :: saved_fpe_mode(SIZE(ieee_all))
            REAL(realk) :: nan

            ! Make sure we do not trigger floating point exceptions when
            ! setting the array to NaN
            CALL IEEE_GET_HALTING_MODE(IEEE_ALL, saved_fpe_mode)
            CALL IEEE_SET_HALTING_MODE(IEEE_ALL, .FALSE.)

            ! Define NaN and set that value in the array
            nan = IEEE_VALUE(0.0_realk, IEEE_SIGNALING_NAN)
            this%buffers = nan

            ! Restore the previous floating point exception mode
            CALL IEEE_SET_FLAG(IEEE_ALL, .FALSE.)
            CALL IEEE_SET_HALTING_MODE(IEEE_ALL, saved_fpe_mode)
        END BLOCK
#endif
    END SUBROUTINE init_buffers


    SUBROUTINE get_buffer(this, ptr, igrid, iface)
        ! Subroutine arguments
        CLASS(field_t), INTENT(inout), TARGET :: this
        REAL(realk), INTENT(out), POINTER, CONTIGUOUS :: ptr(:, :)
        INTEGER(intk), INTENT(in) :: igrid
        INTEGER(intk), INTENT(in) :: iface

        ! Local variables
        INTEGER(intk) :: kk, jj, ii, ibb

        IF (.NOT. ALLOCATED(this%buffers)) THEN
            WRITE(*, *) "Buffers not initialized"
            CALL errr(__FILE__, __LINE__)
        END IF

        CALL get_mgdims(kk, jj, ii, igrid)
        CALL get_ibb(ibb, iface, igrid)

        ! Buffers are only allocated on FIX, OP1 and PAR boundaries. If the
        ! returned ibb is zero, this means that get_buffer was called on
        ! another boundary condition which does not have a buffer
        IF (ibb == 0) THEN
            WRITE(*, *) "Buffer not allocated for this boundary condition"
            WRITE(*, *) "  iface: ", iface, " igrid: ", igrid
            CALL errr(__FILE__, __LINE__)
        END IF

        SELECT CASE (iface)
        CASE (1, 2)
            ptr(1:kk, 1:jj) => this%buffers(ibb:ibb+kk*jj-1)
        CASE (3, 4)
            ptr(1:kk, 1:ii) => this%buffers(ibb:ibb+kk*ii-1)
        CASE (5, 6)
            ptr(1:jj, 1:ii) => this%buffers(ibb:ibb+jj*ii-1)
        CASE DEFAULT
            WRITE(*, '("Invalid face: ", I0)') iface
            CALL errr(__FILE__, __LINE__)
        END SELECT
    END SUBROUTINE get_buffer
END MODULE realfield_mod
