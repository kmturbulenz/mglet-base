MODULE intfield_mod
    USE basefield_mod
    USE err_mod, ONLY: errr
    USE grids_mod, ONLY: get_mgdims, get_imygrid, mygridslvl, nmygrids, &
        nmygridslvl
    USE precision_mod, ONLY: intk, ifk, mglet_hdf5_ifk, mglet_mpi_ifk
    USE profile_tools_mod, ONLY: profile_range_push, profile_range_pop

    IMPLICIT NONE(type, external)
    PRIVATE

    TYPE, EXTENDS(basefield_t) :: intfield_t
        INTEGER(ifk), ALLOCATABLE :: arr(:)
    CONTAINS
        PROCEDURE :: init

        ! Support device offloading and default to device
        PROCEDURE :: zero
        PROCEDURE :: copy_arr_to_device
        PROCEDURE :: copy_arr_from_device

        GENERIC, PUBLIC :: get_ptr => get_grid1, get_grid3
        PROCEDURE, PRIVATE, NON_OVERRIDABLE :: get_grid1, get_grid3

        PROCEDURE :: finish
        FINAL :: destructor
    END TYPE intfield_t

    PUBLIC :: intfield_t

CONTAINS
    SUBROUTINE init(this, name, description, ndim, istag, jstag, kstag, &
            units, dread, required, dwrite, active_level, get_len, zero, device)
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
        LOGICAL, INTENT(in), OPTIONAL :: zero
        LOGICAL, INTENT(in), OPTIONAL :: device

        ! Local variables
        LOGICAL :: zero2

        CALL this%init_corefield(name, description, ndim, istag, jstag, kstag, &
            units, dread, required, dwrite, active_level, get_len)

        this%hdf5_dtype = mglet_hdf5_ifk
        this%mpi_dtype = mglet_mpi_ifk

        IF (PRESENT(zero)) THEN
            zero2 = zero
        ELSE
            zero2 = .TRUE.
        END IF

        this%device = .TRUE.
        IF (PRESENT(device)) THEN
            this%device = device
        END IF

        ALLOCATE(this%arr(this%idim))
        IF (this%device) THEN
            CALL device_alloc(this%arr)
        END IF

        IF (zero2) THEN
            this%arr = 0_ifk
            IF (this%device) THEN
                CALL this%copy_arr_to_device()
            END IF
        END IF
    END SUBROUTINE init


    SUBROUTINE finish(this)
        ! Subroutine arguments
        CLASS(intfield_t), INTENT(inout) :: this

        ! Local variables
        ! none...

        IF (.NOT. this%is_init) RETURN
        CALL this%finish_corefield()

        IF (this%device) THEN
            CALL device_free(this%arr)
        END IF
        DEALLOCATE(this%arr)
    END SUBROUTINE finish


    IMPURE ELEMENTAL SUBROUTINE destructor(this)
        ! Subroutine arguments
        TYPE(intfield_t), INTENT(inout) :: this

        ! Local variables
        ! none...

        CALL this%finish()
    END SUBROUTINE destructor


    SUBROUTINE copy_arr_to_device(this)
        ! Subroutine arguments
        CLASS(intfield_t), INTENT(inout) :: this

        ! Local variables
        ! none...

        IF (.NOT. this%device) THEN
            WRITE(*, *) "Intfield was not device initialized:", TRIM(this%name)
            CALL errr(__FILE__, __LINE__)
        END IF

        CALL copy_to_device(this%arr)
    END SUBROUTINE copy_arr_to_device


    SUBROUTINE copy_arr_from_device(this)
        ! Subroutine arguments
        CLASS(intfield_t), INTENT(inout) :: this

        ! Local variables
        ! none...

        IF (.NOT. this%device) THEN
            WRITE(*, *) "Intfield was not device initialized:", TRIM(this%name)
            CALL errr(__FILE__, __LINE__)
        END IF

        CALL copy_from_device(this%arr)
    END SUBROUTINE copy_arr_from_device


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


    SUBROUTINE zero(this, ilevel, device)
        ! Subroutine arguments
        CLASS(intfield_t), INTENT(inout) :: this
        INTEGER(intk), OPTIONAL, INTENT(in) :: ilevel
        LOGICAL, OPTIONAL, INTENT(in) :: device

        ! Local variables
        LOGICAL :: device2

        device2 = .TRUE.
        IF (PRESENT(device)) THEN
            device2 = device
        END IF

        IF (PRESENT(ilevel)) THEN
            CALL zero_lvl_impl(this%arr, this%ptr, this%length, ilevel, device2)
        ELSE
            CALL zero_impl(this%arr, SIZE(this%arr), device2)
        END IF
    END SUBROUTINE zero


    SUBROUTINE zero_impl(arr, n, device)
        ! Subroutine arguments
        INTEGER(ifk), INTENT(inout) :: arr(*)
        INTEGER(intk), INTENT(in) :: n
        LOGICAL, INTENT(in) :: device

        ! Local variables
        INTEGER(intk) :: i

        CALL profile_range_push("intfield_zero_arr")

        !$omp target teams loop if(device)
        DO i = 1, n
            arr(i) = 0_ifk
        END DO
        !$omp end target teams loop

        CALL profile_range_pop()
    END SUBROUTINE zero_impl


    SUBROUTINE zero_lvl_impl(arr, ptr, length, ilevel, device)
        ! Subroutine arguments
        INTEGER(ifk), INTENT(inout) :: arr(*)
        INTEGER(intk), INTENT(in) :: ptr(*), length(*), ilevel
        LOGICAL, INTENT(in) :: device

        ! Local variables
        INTEGER(intk) :: i, il, igrid, imygrid, ip, n, ngridslvl

        ngridslvl = nmygridslvl(ilevel)
        CALL profile_range_push("intfield_zero_arr_lvl")

        !$omp target teams distribute &
        !$omp& if(device) &
        !$omp& private(il, igrid, imygrid, ip, n) &
        !$omp& shared(ilevel) &
        !$omp& map(to: ptr(1:nmygrids), length(1:nmygrids))
        DO il = 1, ngridslvl
            igrid = mygridslvl(il, ilevel)
            CALL get_imygrid(imygrid, igrid)
            ip = ptr(imygrid)
            n = length(imygrid)

            !$omp parallel do
            DO i = 1, n
                arr(ip + i - 1) = 0_ifk
            END DO
            !$omp end parallel do
        END DO
        !$omp end target teams distribute

        CALL profile_range_pop()
    END SUBROUTINE zero_lvl_impl


    SUBROUTINE device_alloc(data)
        ! Subroutine arguments
        INTEGER(ifk), CONTIGUOUS, INTENT(in) :: data(:)

        ! Local variables
        ! none...

        !$omp target enter data map(alloc: data)
    END SUBROUTINE device_alloc


    SUBROUTINE device_free(data)
        ! Subroutine arguments
        INTEGER(ifk), CONTIGUOUS, INTENT(inout) :: data(:)

        ! Local variables
        ! none...

        !$omp target exit data map(delete: data)
    END SUBROUTINE device_free


    SUBROUTINE copy_to_device(data)
        ! Subroutine arguments
        INTEGER(ifk), CONTIGUOUS, INTENT(in) :: data(:)

        ! Local variables
        ! none...

        !$omp target update to(data)
    END SUBROUTINE copy_to_device


    SUBROUTINE copy_from_device(data)
        ! Subroutine arguments
        INTEGER(ifk), CONTIGUOUS, INTENT(inout) :: data(:)

        ! Local variables
        ! none...

        !$omp target update from(data)
    END SUBROUTINE copy_from_device
END MODULE intfield_mod
