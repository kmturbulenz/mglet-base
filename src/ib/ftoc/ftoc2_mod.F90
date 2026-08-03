MODULE ftoc2_mod
    USE MPI_f08
    USE core_mod
    USE ftoc_core_mod
    USE ibcore_mod, ONLY: ib
    USE restrict_mod, ONLY: restrict, message_length, start_and_stop

    IMPLICIT NONE (type, external)
    PRIVATE

    ! Lists that hold the send and receive request arrays
    TYPE(MPI_Request), ALLOCATABLE :: sendreqs(:), recvreqs(:)

    ! Lists that hold the messages that are ACTUALLY sendt and received
    INTEGER(intk) :: nsend, nrecv
    INTEGER(int32), ALLOCATABLE :: recvlist(:)
    INTEGER(intk), ALLOCATABLE :: recvidxlist(:, :)

    ! Variable to indicate if the required data structures have been created
    LOGICAL :: is_init = .FALSE.

    INTERFACE ftoc2
        MODULE PROCEDURE :: ftoc2_one, ftoc2_multiple
    END INTERFACE ftoc2

    ! contained functions
    PUBLIC :: ftoc2, init_ftoc2, finish_ftoc2
CONTAINS
    SUBROUTINE ftoc2_one(ilevel, ff, fc, flag)
        ! Subroutine arguments
        INTEGER(intk), INTENT(in) :: ilevel
        TYPE(field_t), INTENT(in) :: ff
        TYPE(field_t), INTENT(inout) :: fc
        CHARACTER(len=1), INTENT(in) :: flag

        ! Local variables
        ! none...


        IF (.NOT. is_init) CALL errr(__FILE__, __LINE__)

        ! TODO(offload): implement

        nrecv = 0
        nsend = 0
    END SUBROUTINE ftoc2_one


    SUBROUTINE ftoc2_multiple(ilevel, v1, v2, v3, s1, s2, s3)
        ! Subroutine arguments
        INTEGER(intk), INTENT(in) :: ilevel
        TYPE(field_t), INTENT(inout), OPTIONAL :: v1, v2, v3, s1, s2, s3

        ! Local variables
        CHARACTER(len=*), PARAMETER :: flag = '*'

        IF (.NOT. is_init) CALL errr(__FILE__, __LINE__)

        ! TODO(offload): implement

        nrecv = 0
        nsend = 0
    END SUBROUTINE ftoc2_multiple


    SUBROUTINE init_ftoc2()
        ! Subroutine arguments
        ! none...

        ! Local variables
        ! none...

        IF (is_init) CALL errr(__FILE__, __LINE__)

        ! Check if ctof_core_mod provides necessary infrastructure
        IF (.NOT. is_ftoc_core_init) CALL errr(__FILE__, __LINE__)

        ALLOCATE(recvidxlist(3, irecv), source=0_intk)
        ALLOCATE(recvlist(irecv), source=0_int32)
        ALLOCATE(sendreqs(numprocs), source=MPI_REQUEST_NULL)
        ALLOCATE(recvreqs(numprocs), source=MPI_REQUEST_NULL)
        nsend = 0
        nrecv = 0

        is_init = .TRUE.
    END SUBROUTINE init_ftoc2


    SUBROUTINE finish_ftoc2()
        ! Subroutine arguments
        ! none...

        ! Local variables
        ! none...

        IF (.NOT. is_init) CALL errr(__FILE__, __LINE__)

        DEALLOCATE(recvreqs)
        DEALLOCATE(sendreqs)
        DEALLOCATE(recvlist)
        DEALLOCATE(recvidxlist)
        nsend = 0
        nrecv = 0

        is_init = .FALSE.
    END SUBROUTINE finish_ftoc2
END MODULE ftoc2_mod
