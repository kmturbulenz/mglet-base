MODULE intfield_mod
    USE err_mod, ONLY: errr
    USE grids_mod, ONLY: get_mgdims
    USE precision_mod, ONLY: intk, realk, ifk, mglet_hdf5_ifk, mglet_mpi_ifk
    USE basefield_mod

    IMPLICIT NONE(type, external)
    PRIVATE

    TYPE, EXTENDS(basefield_t) :: intfield_t
        INTEGER(ifk), ALLOCATABLE :: arr(:)
    CONTAINS
        PROCEDURE :: init

        GENERIC, PUBLIC :: get_ptr => get_grid1, get_grid3
        PROCEDURE, PRIVATE :: get_grid1, get_grid3

        PROCEDURE :: finish
        FINAL :: destructor
    END TYPE intfield_t

    PUBLIC :: intfield_t

CONTAINS
    SUBROUTINE init(this, name, description, ndim, istag, jstag, kstag, &
            units, dread, required, dwrite, active_level, get_len)
        ! Subroutine arguments
        CLASS(intfield_t), INTENT(out) :: this
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

        this%hdf5_dtype = mglet_hdf5_ifk
        this%mpi_dtype = mglet_mpi_ifk

        ALLOCATE(this%arr(this%idim))
        this%arr = 0.0
    END SUBROUTINE init


    ELEMENTAL SUBROUTINE finish(this)
        CLASS(intfield_t), INTENT(inout) :: this

        IF (.NOT. this%is_init) RETURN
        CALL this%finish_corefield()

        DEALLOCATE(this%arr)
    END SUBROUTINE finish


    ELEMENTAL SUBROUTINE destructor(this)
        TYPE(intfield_t), INTENT(inout) :: this

        CALL this%finish()
    END SUBROUTINE destructor


    SUBROUTINE get_grid1(this, ptr, igrid)
        ! Subroutine arguments
        CLASS(intfield_t), INTENT(in), TARGET :: this
        INTEGER(ifk), POINTER, CONTIGUOUS, INTENT(out) :: ptr(:)
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
        CLASS(intfield_t), INTENT(in), TARGET :: this
        INTEGER(ifk), POINTER, CONTIGUOUS, INTENT(out) :: ptr(:, :, :)
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

END MODULE intfield_mod
