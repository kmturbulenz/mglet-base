! Communication module
!
! This module creates various communicators for different purposes:
!   - shmcomm:  A node-local communicator, possible to create shared memory
!               allocations
!   - iocomm:   A communicator for processes doing IO
!   - iogrcomm: Each process belong to a group of processes that hold
!               exactly one IO process.
!
! The IO-related communicators are currently created based on the shared-memory
! communicator, i.e. one process per compute node do IO, and all processes on
! this node participate in an IO group.

MODULE comms_mod
    USE precision_mod, ONLY: int32
    USE MPI_f08

    IMPLICIT NONE (type, external)
    PRIVATE

    ! Communicator over shared memory space, i.e. one compute node
    TYPE(MPI_Comm), PROTECTED :: shmcomm

    ! Communicator with lowest ranks at each SHM segment (SHM masters)
    TYPE(MPI_Comm), PROTECTED :: shm_masters_comm

    ! Communicator gropuping toghether processes associated with an IO process
    TYPE(MPI_Comm), PROTECTED :: iogrcomm

    ! Communicator for processes participating in IO operations
    TYPE(MPI_Comm), PROTECTED :: iocomm

    ! Rank and number of processes in node-local SHM communicator (shmcomm)
    INTEGER(int32), PROTECTED :: shmid = -1, shmprocs = -1

    ! Rank and number of processes in IO group (iogrcomm)
    INTEGER(int32), PROTECTED :: iogrid = -1, iogrprocs = -1

    ! Rank and number of processes in IO communicator (iocomm)
    INTEGER(int32), PROTECTED :: ioid = -1, ioprocs = -1

    ! Rank and number of processes in global communicator (MPI_COMM_WORLD)
    INTEGER(int32), PROTECTED :: myid = -1, numprocs = -1

    ! Number of compute nodes (i.e. number of ranks in iocomm)
    INTEGER(int32), PROTECTED :: numnodes = -1

    ! Id of current compute node (i.e. rank of master in iogrcomm)
    INTEGER(int32), PROTECTED :: nodeid = -1

    ! Flag to indicate if this process is to do any IO work
    LOGICAL :: ioProc = .FALSE.

    ! Flag to indicate if this process is to do any IO work
    LOGICAL :: iogrcomm_is_world = .FALSE.

    ! MPI error handler
    TYPE(MPI_Errhandler) :: errh

    ! Public data items
    PUBLIC :: shmId, shmProcs, iogrId, iogrProcs, ioId, ioProcs, myId, &
        numProcs, ioProc, init_comms, numnodes, nodeid, finish_comms
    PUBLIC :: shmcomm, shm_masters_comm, iogrcomm, iocomm


CONTAINS
    ! Initialize communicator(s) and module constants
    SUBROUTINE init_comms()
        ! Local variables
        INTEGER(int32) :: color = 1
        CHARACTER(len=32) :: serial
        INTEGER :: length, status
        CHARACTER(LEN=MPI_MAX_LIBRARY_VERSION_STRING) :: version

        ! Initialization for deterministic behaviour
        shmcomm = MPI_COMM_NULL
        shm_masters_comm = MPI_COMM_NULL
        iogrcomm = MPI_COMM_NULL
        iocomm = MPI_COMM_NULL

        CALL MPI_Comm_rank(MPI_COMM_WORLD, myid)
        CALL MPI_Comm_size(MPI_COMM_WORLD, numprocs)

        IF (myid == 0) THEN
            WRITE(*, '()')
            WRITE(*, '("M G L E T - Multi Grid Large Eddy Turbulence")')
            WRITE(*, '()')
        END IF

        ! Check if serial IO should be used
        ! TODO: re-work such that we can use getenv_bool_coll here - that
        ! is currently not possible because we create a circular dependency
        ! between comms_mod and envvars_mod
        iogrcomm_is_world = .FALSE.
        CALL get_environment_variable("KMT_MGLET_SERIAL_IO", serial, length, &
            status)
        IF (status < 1) THEN
            IF (serial == "YES" .OR. serial == "1") THEN
                iogrcomm_is_world = .TRUE.
            END IF
        END IF
        CALL MPI_Bcast(iogrcomm_is_world, 1, MPI_LOGICAL, 0, MPI_COMM_WORLD)
        IF (myid == 0) THEN
            WRITE(*, '("MPI INFORMATION:")')
            WRITE(*, '("    Processes:     ", I0)') numprocs
            WRITE(*, '("    Serial IO:     ", L1)') iogrcomm_is_world

            CALL MPI_Get_library_version(version, length)
            ! Strip away newline char in the end if present
            IF (version(length:length) == NEW_LINE('a')) THEN
                length = length - 1
            END IF
            WRITE(*, '("    MPI version:   ", A)') version(1:length)
            WRITE(*, '()')
        END IF

        CALL MPI_Comm_split_type(MPI_COMM_WORLD, MPI_COMM_TYPE_SHARED, 0, &
            MPI_INFO_NULL, shmcomm)

        CALL MPI_comm_rank(shmcomm, shmid)
        CALL MPI_comm_size(shmcomm, shmprocs)

        ! Create IO group communicator
        ! Currently just a copy of shmcomm, but can be changed later
        IF (iogrcomm_is_world) THEN
            CALL MPI_Comm_dup(MPI_COMM_WORLD, iogrcomm)
        ELSE
            CALL MPI_Comm_dup(shmcomm, iogrcomm)
        END IF
        CALL MPI_comm_rank(iogrcomm, iogrid)
        CALL MPI_comm_size(iogrcomm, iogrprocs)

        ! Create I/O process communicator
        ! All the processes with rank 0 in iogrcomm actually do IO
        ioproc = .FALSE.
        IF (iogrid == 0) ioproc = .TRUE.

        color = 1
        IF (iogrid /= 0) color = MPI_UNDEFINED
        CALL MPI_Comm_split(MPI_COMM_WORLD, color, 0, iocomm)

        IF (ioproc) THEN
            CALL MPI_comm_rank(iocomm, ioid)
            CALL MPI_comm_size(iocomm, ioprocs)
        END IF

        ! SHM masters communicator, one rank (the lowest) per SHM segment
        ! If parallel IO is used with one IO rank per SHM segment, this is
        ! equal to iocomm, otherwise they might differ
        color = 1
        IF (shmid /= 0) color = MPI_UNDEFINED
        CALL MPI_Comm_split(MPI_COMM_WORLD, color, 0, shm_masters_comm)

        ! Set number of compute nodes
        ! (rank 0 will always be an IO process)
        numnodes = ioprocs
        nodeid = ioid
        CALL MPI_Bcast(numnodes, 1, MPI_INTEGER, 0, MPI_COMM_WORLD)
        CALL MPI_Bcast(nodeid, 1, MPI_INTEGER, 0, iogrcomm)

        CALL init_errhandler()
    END SUBROUTINE init_comms


    SUBROUTINE finish_comms
        CALL MPI_Comm_free(shmcomm)

        IF (ioproc) THEN
            CALL MPI_Comm_free(iocomm)
        END IF

        CALL MPI_Comm_free(iogrcomm)

        IF (shmid == 0) THEN
            CALL MPI_Comm_free(shm_masters_comm)
        END IF

        CALL MPI_Errhandler_free(errh)

        IF (myid == 0) THEN
            WRITE(*, '("MGLET FINISHED SUCCESSFULLY!")')
            WRITE(*, '()')
        END IF
    END SUBROUTINE finish_comms


    SUBROUTINE init_errhandler()
        ! Set custom error handler
        CALL MPI_Comm_create_errhandler(mglet_errhandler_function, errh)
        CALL MPI_Comm_set_errhandler(MPI_COMM_WORLD, errh)
    END SUBROUTINE init_errhandler


    SUBROUTINE mglet_errhandler_function(comm, error_code)
        USE, INTRINSIC :: ISO_FORTRAN_ENV, ONLY: unit => ERROR_UNIT

#if defined __INTEL_COMPILER
        USE IFCORE
#endif

        ! Subroutine arguments
        TYPE(MPI_Comm) :: comm
        INTEGER(int32) :: error_code

        INTEGER(int32) :: error_class, msglen
        CHARACTER(LEN=MPI_MAX_ERROR_STRING) :: msg

        WRITE(*,*) "MGLET MPI ERROR occured in rank = ", myid
        WRITE(unit,*) "MGLET MPI ERROR"
        WRITE(unit,*) "MPI rank = ", myid
        WRITE(unit,*) ""

        CALL MPI_Error_class(error_code, error_class)
        WRITE(unit,*) "error_class = ", error_class

        CALL MPI_Error_string(error_class, msg, msglen)

        WRITE(unit,*) "Error class type: ", msg(1:msglen)
        WRITE(unit,*) ""
        WRITE(unit,*) "error_code = ", error_code
        WRITE(unit,*) ""

        CALL MPI_Error_string(error_code, msg, msglen)
        WRITE(unit,*) "Full error message:"
        WRITE(unit,*) msg(1:msglen)
        WRITE(unit,*) ""

        WRITE(unit,*) "Full backtrace of process:"
#   if defined __INTEL_COMPILER
        CALL TRACEBACKQQ(user_exit_code=-1)
#   endif
#   if defined __GFORTRAN__
        CALL BACKTRACE()
#   endif

        CALL MPI_Abort(comm, 99)
    END SUBROUTINE mglet_errhandler_function
END MODULE comms_mod
