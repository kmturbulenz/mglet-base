! Shared memory module
!
! Wrapper around MPI's shared memory module for easier allocation/deallocation
MODULE shmem_mod
    USE, INTRINSIC :: ISO_C_BINDING, ONLY: C_PTR, C_F_POINTER, C_NULL_PTR

    USE MPI_f08

    USE comms_mod
    USE err_mod, ONLY: errr
    USE precision_mod

    IMPLICIT NONE (type, external)
    PRIVATE

    TYPE, ABSTRACT :: shmem_arr
        TYPE(MPI_Win) :: win
        INTEGER(MPI_ADDRESS_KIND) :: winsize
        TYPE(C_PTR) :: baseptr
    CONTAINS
        PROCEDURE :: barrier
        PROCEDURE :: allocate
        PROCEDURE :: free
    END TYPE shmem_arr

    TYPE, EXTENDS(shmem_arr) :: real_shmem_arr
        REAL(realk), POINTER :: arr(:)
    END TYPE real_shmem_arr

    TYPE, EXTENDS(shmem_arr) :: int_shmem_arr
        INTEGER(intk), POINTER :: arr(:)
    END TYPE int_shmem_arr

    PUBLIC :: int_shmem_arr, real_shmem_arr

CONTAINS
    ! Allocate shared memory segment
    !
    ! array: array to be allocated
    ! length: number of elements (integers, reals) to allocate per process
    !         overall length of shared array is the sum of these on all
    !         processes.
    SUBROUTINE allocate(this, length)
        ! Subroutine arguments
        CLASS(shmem_arr), INTENT(inout) :: this
        INTEGER(int64), INTENT(in) :: length

        ! Local variables
        INTEGER(int32) :: disp_unit
        INTEGER(MPI_ADDRESS_KIND) :: nelems

        SELECT TYPE (this)
        TYPE IS (int_shmem_arr)
            this%winsize = INT(int_bytes, int64)*length
        TYPE IS (real_shmem_arr)
            this%winsize = INT(real_bytes, int64)*length
        CLASS DEFAULT
            CALL errr(__FILE__, __LINE__)
        END SELECT
        CALL MPI_Win_allocate_shared(this%winsize, 1, MPI_INFO_NULL, &
            shmcomm, this%baseptr, this%win)

        CALL MPI_Allreduce(this%winsize, nelems, 1, MPI_Aint, MPI_SUM, &
            shmcomm)

        CALL MPI_Win_shared_query(this%win, MPI_PROC_NULL, &
            this%winsize, disp_unit, this%baseptr)

        ! Associate the fortran-pointer with the C-style-pointer from MPI
        ! This allows us to use the array as a regular, local Fortran
        ! array.
        SELECT TYPE (this)
        TYPE IS (int_shmem_arr)
            nelems = nelems/int_bytes
            CALL C_F_POINTER(this%baseptr, this%arr, [nelems])
        TYPE IS (real_shmem_arr)
            nelems = nelems/real_bytes
            CALL C_F_POINTER(this%baseptr, this%arr, [nelems])
        CLASS DEFAULT
            CALL errr(__FILE__, __LINE__)
        END SELECT
        CALL this%barrier()

        ! Set array to zero
        SELECT TYPE (this)
        TYPE IS (int_shmem_arr)
            this%arr = 0
        TYPE IS (real_shmem_arr)
            this%arr = 0.0
        CLASS DEFAULT
            CALL errr(__FILE__, __LINE__)
        END SELECT
        CALL this%barrier()
    END SUBROUTINE allocate


    SUBROUTINE free(this)
        ! Subroutine arguments
        CLASS(shmem_arr), INTENT(inout) :: this

        SELECT TYPE (this)
        TYPE IS (int_shmem_arr)
            NULLIFY(this%arr)
        TYPE IS (real_shmem_arr)
            NULLIFY(this%arr)
        CLASS DEFAULT
            CALL errr(__FILE__, __LINE__)
        END SELECT

        CALL this%barrier()

        CALL MPI_Win_free(this%win)

        this%win = MPI_WIN_NULL
        this%winsize = 0
        this%baseptr = C_NULL_PTR
    END SUBROUTINE free


    SUBROUTINE barrier(this)
        ! Subroutine arguments
        CLASS(shmem_arr), INTENT(inout) :: this

        ! Local arguments
        ! none...

        CALL MPI_Win_fence(0, this%win)
    END SUBROUTINE barrier

END MODULE shmem_mod
