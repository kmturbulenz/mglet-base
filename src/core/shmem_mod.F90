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
    END TYPE shmem_arr

    TYPE, EXTENDS(shmem_arr) :: real_shmem_arr
        REAL(realk), POINTER :: arr(:)
    END TYPE real_shmem_arr

    TYPE, EXTENDS(shmem_arr) :: int_shmem_arr
        INTEGER(intk), POINTER :: arr(:)
    END TYPE int_shmem_arr

    PUBLIC :: int_shmem_arr, real_shmem_arr, shmem_alloc, shmem_dealloc

CONTAINS
    ! Allocate shared memory segment
    !
    ! array: array to be allocated
    ! length: number of elements (integers, reals) to allocate per process
    !         overall length of shared array is the sum of these on all
    !         processes.
    SUBROUTINE shmem_alloc(array, length)
        ! Subroutine arguments
        CLASS(shmem_arr), INTENT(inout) :: array
        INTEGER(int64), INTENT(in) :: length

        ! Local variables
        INTEGER(int32) :: disp_unit
        INTEGER(MPI_ADDRESS_KIND) :: nelems

        SELECT TYPE (array)
        TYPE IS (int_shmem_arr)
            array%winsize = INT(int_bytes, int64)*length
        TYPE IS (real_shmem_arr)
            array%winsize = INT(real_bytes, int64)*length
        CLASS DEFAULT
            CALL errr(__FILE__, __LINE__)
        END SELECT
        CALL MPI_Win_allocate_shared(array%winsize, 1, MPI_INFO_NULL, &
            shmcomm, array%baseptr, array%win)

        CALL MPI_Allreduce(array%winsize, nelems, 1, MPI_Aint, MPI_SUM, &
            shmcomm)

        CALL MPI_Win_shared_query(array%win, MPI_PROC_NULL, &
            array%winsize, disp_unit, array%baseptr)

        ! Associate the fortran-pointer with the C-style-pointer from MPI
        ! This allows us to use the array as a regular, local Fortran
        ! array.
        SELECT TYPE (array)
        TYPE IS (int_shmem_arr)
            nelems = nelems/int_bytes
            CALL C_F_POINTER(array%baseptr, array%arr, [nelems])
        TYPE IS (real_shmem_arr)
            nelems = nelems/real_bytes
            CALL C_F_POINTER(array%baseptr, array%arr, [nelems])
        CLASS DEFAULT
            CALL errr(__FILE__, __LINE__)
        END SELECT
        CALL array%barrier()

        ! Set array to zero
        SELECT TYPE (array)
        TYPE IS (int_shmem_arr)
            array%arr = 0
        TYPE IS (real_shmem_arr)
            array%arr = 0.0
        CLASS DEFAULT
            CALL errr(__FILE__, __LINE__)
        END SELECT
        CALL array%barrier()
    END SUBROUTINE shmem_alloc


    SUBROUTINE shmem_dealloc(array)
        ! Subroutine arguments
        CLASS(shmem_arr), INTENT(inout) :: array

        SELECT TYPE (array)
        TYPE IS (int_shmem_arr)
            NULLIFY(array%arr)
        TYPE IS (real_shmem_arr)
            NULLIFY(array%arr)
        CLASS DEFAULT
            CALL errr(__FILE__, __LINE__)
        END SELECT

        CALL array%barrier()

        CALL MPI_Win_free(array%win)

        array%win = MPI_WIN_NULL
        array%winsize = 0
        array%baseptr = C_NULL_PTR
    END SUBROUTINE shmem_dealloc


    SUBROUTINE barrier(this)
        ! Subroutine arguments
        CLASS(shmem_arr), INTENT(inout) :: this

        ! Local arguments
        ! none...

        CALL MPI_Win_fence(0, this%win)
    END SUBROUTINE barrier

END MODULE shmem_mod
