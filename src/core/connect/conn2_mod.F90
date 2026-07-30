MODULE conn2_mod

    USE MPI_f08
    USE precision_mod
    USE commbuf_mod, ONLY: sendbuf, recvbuf
    USE err_mod, ONLY: errr
    USE grids_mod, ONLY: mygrids, nmygrids, level, idprocofgrd, itypboconds, &
        maxlevel, minlevel, get_neighbours, get_mgdims
    USE comms_mod, ONLY: myid, numprocs
    USE field_mod
    USE connect_core_mod
    USE fieldhelper_mod
    USE fieldpool_mod
    USE pointers_mod, ONLY: get_ip3
    USE profile_tools_mod

    IMPLICIT NONE (type, external)
    PRIVATE

    ! Lists that hold the send and receive request arrays
    TYPE(MPI_Request), ALLOCATABLE :: sendreqs(:), recvreqs(:)
    INTEGER(int32), ALLOCATABLE :: sendlist(:), recvlist(:)
    INTEGER(intk), ALLOCATABLE :: recvidxlist(:, :)
    INTEGER(intk) :: nsend, nrecv

    ! Task extents
    INTEGER(intk), PARAMETER :: buffertasksize = 9
    INTEGER(intk), PARAMETER :: selftasksize = 15
    INTEGER(intk), PARAMETER :: mpitasksize = 3
    INTEGER(intk) :: maxsendtasks, maxrecvtasks
    INTEGER(intk) :: maxselftasks
    INTEGER(intk) :: maxmpisendtasks, maxmpirecvtasks

    ! Workpackages containing individual tasks for packing / unpacking
    INTEGER(intk), ALLOCATABLE :: sendtasks(:, :), recvtasks(:, :)
    INTEGER(intk), ALLOCATABLE :: selftasks(:, :)
    INTEGER(intk), ALLOCATABLE :: mpisendtasks(:, :), mpirecvtasks(:, :)
    !$omp declare target(sendtasks, recvtasks, selftasks)

    ! Type to hold condensed task arrays to execute a certain type of conn
    TYPE :: work_t
        LOGICAL :: is_init = .FALSE.
        INTEGER(intk), ALLOCATABLE :: sendtasks(:, :)
        INTEGER(intk), ALLOCATABLE :: recvtasks(:, :)
        INTEGER(intk), ALLOCATABLE :: selftasks(:, :)
        INTEGER(intk), ALLOCATABLE :: mpisendtasks(:, :)
        INTEGER(intk), ALLOCATABLE :: mpirecvtasks(:, :)
    END TYPE work_t

    ! Multidimensional array to store work_t for different conn types
    TYPE(work_t), ALLOCATABLE, TARGET :: workrecords(:, :, :, :, :, :)
    LOGICAL :: is_recording = .FALSE.

    ! Pointers to fields passed to conn2, point to dummy field if not present
    TYPE(field_t), POINTER :: f1 => NULL(), f2 => NULL(), f3 => NULL(), &
        f4 => NULL(), f5 => NULL(), f6 => NULL()

    LOGICAL :: is_init = .FALSE.

    PUBLIC :: conn2, init_conn2, finish_conn2
CONTAINS
    ! conn2 makes conn1 available for efficient execution on GPUs.
    ! The implementation provides the same functionality as conn1.
    SUBROUTINE conn2(ilevel, layers, v1, v2, v3, s1, s2, s3, corners, normal, &
            forward, ityp)
        ! Subroutine arguments
        INTEGER(intk), INTENT(in), OPTIONAL :: ilevel, layers
        TYPE(field_t), OPTIONAL, TARGET, INTENT(inout) :: v1, v2, v3, s1, s2, s3
        LOGICAL, OPTIONAL, INTENT(in) :: corners, normal
        INTEGER(intk), OPTIONAL, INTENT(in) :: forward
        CHARACTER(len=1), OPTIONAL, INTENT(in) :: ityp

        ! Local variables
        TYPE(field_t), POINTER :: dummy
        INTEGER(intk) :: minconlvl, maxconlvl, nplane, fwd, nvars
        LOGICAL :: vertices, normal2
        CHARACTER(len=1) :: flag

        ! Indices for the position within the workrecords array
        INTEGER(intk) :: idx_ilevel, idx_layers, idx_args, idx_corners, &
            idx_normal, idx_forward

        IF (.NOT. is_init) CALL errr(__FILE__, __LINE__)

        CALL push_field(dummy, "CONN2_DUMMY", zero=.FALSE.)

        f1 => dummy
        f2 => dummy
        f3 => dummy
        f4 => dummy
        f5 => dummy
        f6 => dummy

        IF (PRESENT(ilevel)) THEN
            IF (ilevel < minlevel .OR. ilevel > maxlevel) THEN
                CALL errr(__FILE__, __LINE__)
            END IF
            minconlvl = ilevel
            maxconlvl = ilevel
            idx_ilevel = ilevel
        ELSE
            minconlvl = minlevel
            maxconlvl = maxlevel
            idx_ilevel = maxlevel+1
        END IF

        nplane = 1
        idx_layers = 1
        IF (PRESENT(layers)) THEN
            nplane = layers
            idx_layers = layers
        END IF
        IF (nplane < 1 .OR. nplane > 2) THEN
            CALL errr(__FILE__, __LINE__)
        END IF

        vertices = .FALSE.
        idx_corners = 0
        IF (PRESENT(corners)) THEN
            vertices = corners
            IF (corners) THEN
                idx_corners = 1
            END IF
        END IF

        fwd = 0
        idx_forward = 0
        IF (PRESENT(forward)) THEN
            IF (forward < -1 .OR. forward > 1) THEN
                CALL errr(__FILE__, __LINE__)
            END IF
            fwd = forward
            idx_forward = forward
        END IF

        IF (fwd /= 0 .AND. vertices) THEN
            CALL errr(__FILE__, __LINE__)
        END IF

        flag = ' '
        IF (PRESENT(ityp)) THEN
            ! Radical approach
            CALL errr(__FILE__, __LINE__)
        END IF

        idx_args = 0
        nvars = 0
        IF (PRESENT(v1)) THEN
            f1 => v1
            nvars = nvars + 1
            idx_args = idx_args + 1 * 2**(0)
        END IF
        IF (PRESENT(v2)) THEN
            f2 => v2
            nvars = nvars + 1
            idx_args = idx_args + 1 * 2**(1)
        END IF
        IF (PRESENT(v3)) THEN
            f3 => v3
            nvars = nvars + 1
            idx_args = idx_args + 1 * 2**(2)
        END IF
        IF (PRESENT(s1)) THEN
            f4 => s1
            nvars = nvars + 1
            idx_args = idx_args + 1 * 2**(3)
        END IF
        IF (PRESENT(s2)) THEN
            f5 => s2
            nvars = nvars + 1
            idx_args = idx_args + 1 * 2**(4)
        END IF
        IF (PRESENT(s3)) THEN
            f6 => s3
            nvars = nvars + 1
            idx_args = idx_args + 1 * 2**(5)
        END IF
        IF (nvars == 0) THEN
            CALL errr(__FILE__, __LINE__)
        END IF

        normal2 = .FALSE.
        idx_normal = 0
        IF (PRESENT(normal)) THEN
            normal2 = normal
            IF (normal) THEN
                idx_normal = 1
            END IF
        END IF
        IF (normal2) THEN
            IF (.NOT. PRESENT(v1) .OR. .NOT. PRESENT(v2) .OR. &
                    .NOT. PRESENT(v3)) THEN
                CALL errr(__FILE__, __LINE__)
            END IF
            IF (vertices) THEN
                CALL errr(__FILE__, __LINE__)
            END IF
            nvars = nvars - 2
        END IF

        ! Looking up the workpackage in the workrecords array
        ASSOCIATE(wptr => workrecords(idx_ilevel, idx_layers, idx_args, &
            idx_corners, idx_normal, idx_forward))

        IF (is_recording) THEN
            ! During the recording phase, the workpackage is not yet initialized
            ! and needs to be created. Tasks arrays are offloaded for later
            ! efficient task execution on device
            CALL recording_pass(wptr, minconlvl, maxconlvl, nplane, vertices, &
                normal2, fwd, flag, nvars, v1, v2, v3, s1, s2, s3)
        ELSE
            IF (wptr%is_init) THEN
                ! Using recorded workpackage with tasks arrays on device
                CALL recorded_conn(wptr)
            ELSE
                ! Using just-in-time (jit) variant, where a non-recorded
                ! workpackage is assembled on CPU and offloaded to device
                ! right before execution. If overheads are too high, the jit
                ! variant can be avoided by recording the relevant workpackage.
                CALL jit_conn(minconlvl, maxconlvl, nplane, vertices, normal2, &
                    fwd, flag, nvars, v1, v2, v3, s1, s2, s3)
            END IF
        END IF

        END ASSOCIATE

        f1 => NULL()
        f2 => NULL()
        f3 => NULL()
        f4 => NULL()
        f5 => NULL()
        f6 => NULL()

        CALL pop_field(dummy)
    END SUBROUTINE conn2


    SUBROUTINE jit_conn(minconlvl, maxconlvl, nplane, vertices, normal2, fwd, &
            flag, nvars, v1, v2, v3, s1, s2, s3)
        ! Subroutine arguments
        INTEGER(intk), INTENT(in) :: minconlvl, maxconlvl, nplane
        LOGICAL, INTENT(in) :: vertices, normal2
        INTEGER(intk), INTENT(in) :: fwd
        CHARACTER(len=1), INTENT(in) :: flag
        INTEGER(intk), INTENT(in) :: nvars
        TYPE(field_t), OPTIONAL, INTENT(inout) :: v1, v2, v3, s1, s2, s3

        ! Local variables
        INTEGER(intk) :: nsendtasks, nrecvtasks
        INTEGER(intk) :: nselftasks
        INTEGER(intk) :: nmpisendtasks

        ! Manufacturing the workpackage just in-time and execute
        CALL prepare_tasks_all(sendtasks, nsendtasks, &
            selftasks, nselftasks, mpisendtasks, nmpisendtasks, &
            minconlvl, maxconlvl, nplane, vertices, &
            normal2, fwd, flag, nvars, v1, v2, v3, s1, s2, s3)

        !$omp target update to( &
        !$omp&  sendtasks(1:buffertasksize, 1:nsendtasks+1), &
        !$omp&  selftasks(1:selftasksize, 1:nselftasks+1)) nowait

        CALL recv_mpi_all(minconlvl, maxconlvl, nplane, vertices, &
            normal2, fwd, flag, nvars)

        !$omp taskwait

        CALL process_sendtasks(nsendtasks, sendtasks)
        CALL process_mpisend(nmpisendtasks, mpisendtasks)
        CALL process_selftasks(nselftasks, selftasks)

        CALL prepare_recvtasks_all(recvtasks, nrecvtasks, &
            nplane, normal2, flag, v1, v2, v3, s1, s2, s3)

        !$omp target update to(recvtasks(1:buffertasksize, 1:nrecvtasks+1))

        CALL process_recvtasks(nrecvtasks, recvtasks)
    END SUBROUTINE jit_conn


    SUBROUTINE recorded_conn(wptr)
        ! Subroutine arguments
        TYPE(work_t), INTENT(in) :: wptr

        ! Local variables
        INTEGER(intk) :: nsendtasks, nrecvtasks
        INTEGER(intk) :: nselftasks
        INTEGER(intk) :: nmpisendtasks, nmpirecvtasks

        ! Using the (offloaded) workpackage from the workrecords array
        nmpirecvtasks = SIZE(wptr%mpirecvtasks, 2) - 1
        nmpisendtasks = SIZE(wptr%mpisendtasks, 2) - 1
        nsendtasks = SIZE(wptr%sendtasks, 2) - 1
        nselftasks = SIZE(wptr%selftasks, 2) - 1
        nrecvtasks = SIZE(wptr%recvtasks, 2) - 1

        CALL process_mpirecv(nmpirecvtasks, wptr%mpirecvtasks)
        CALL process_sendtasks(nsendtasks, wptr%sendtasks)
        CALL process_mpisend(nmpisendtasks, wptr%mpisendtasks)
        CALL process_selftasks(nselftasks, wptr%selftasks)
        CALL MPI_Waitall(nrecv, recvreqs, MPI_STATUSES_IGNORE)
        CALL process_recvtasks(nrecvtasks, wptr%recvtasks)
        CALL MPI_Waitall(nsend, sendreqs, MPI_STATUSES_IGNORE)
    END SUBROUTINE recorded_conn


    SUBROUTINE recording_pass(wptr, minconlvl, maxconlvl, nplane, vertices, &
            normal2, fwd, flag, nvars, v1, v2, v3, s1, s2, s3)
        ! Subroutine arguments
        TYPE(work_t), INTENT(inout) :: wptr
        INTEGER(intk), INTENT(in) :: minconlvl, maxconlvl, nplane
        LOGICAL, INTENT(in) :: vertices, normal2
        INTEGER(intk), INTENT(in) :: fwd
        CHARACTER(len=1), INTENT(in) :: flag
        INTEGER(intk), INTENT(in) :: nvars
        TYPE(field_t), OPTIONAL, INTENT(inout) :: v1, v2, v3, s1, s2, s3

        ! Local variables
        INTEGER(intk) :: nsendtasks, nrecvtasks
        INTEGER(intk) :: nselftasks
        INTEGER(intk) :: nmpisendtasks, nmpirecvtasks

        IF (wptr%is_init) THEN
            WRITE(*, *) "Combination already recorded."
            CALL errr(__FILE__, __LINE__)
        END IF

        ! It is necessary to execute one cycle with communication
        ! as otherwise many valuable checks are not possible

        CALL prepare_mpirecvtasks(mpirecvtasks, nmpirecvtasks, &
            minconlvl, maxconlvl, nplane, vertices, normal2, fwd, &
            flag, nvars)

        CALL prepare_tasks_all(sendtasks, nsendtasks, &
            selftasks, nselftasks, mpisendtasks, nmpisendtasks, &
            minconlvl, maxconlvl, nplane, vertices, &
            normal2, fwd, flag, nvars, v1, v2, v3, s1, s2, s3)

        !$omp target update to( &
        !$omp&  sendtasks(1:buffertasksize, 1:nsendtasks+1), &
        !$omp&  selftasks(1:selftasksize, 1:nselftasks+1))

        CALL process_mpirecv(nmpirecvtasks, mpirecvtasks)
        CALL process_sendtasks(nsendtasks, sendtasks)
        CALL process_mpisend(nmpisendtasks, mpisendtasks)
        CALL process_selftasks(nselftasks, selftasks)
        CALL prepare_recvtasks_all(recvtasks, nrecvtasks, &
            nplane, normal2, flag, v1, v2, v3, s1, s2, s3)

        !$omp target update to(recvtasks(1:buffertasksize, 1:nrecvtasks+1))
        CALL process_recvtasks(nrecvtasks, recvtasks)

        ! Allocate the workpackage arrays in the exact sizes
        ALLOCATE(wptr%sendtasks(buffertasksize, nsendtasks+1))
        ALLOCATE(wptr%recvtasks(buffertasksize, nrecvtasks+1))
        ALLOCATE(wptr%selftasks(selftasksize, nselftasks+1))
        ALLOCATE(wptr%mpisendtasks(mpitasksize, nmpisendtasks+1))
        ALLOCATE(wptr%mpirecvtasks(mpitasksize, nmpirecvtasks+1))

        ! Store the workpackage in the workrecords array and map
        wptr%sendtasks = sendtasks(:, 1:nsendtasks+1)
        wptr%recvtasks = recvtasks(:, 1:nrecvtasks+1)
        wptr%selftasks = selftasks(:, 1:nselftasks+1)
        wptr%mpisendtasks = mpisendtasks(:, 1:nmpisendtasks+1)
        wptr%mpirecvtasks = mpirecvtasks(:, 1:nmpirecvtasks+1)

        !$omp target enter data map(to: &
        !$omp&  wptr%sendtasks(1:buffertasksize, 1:nsendtasks+1), &
        !$omp&  wptr%recvtasks(1:buffertasksize, 1:nrecvtasks+1), &
        !$omp&  wptr%selftasks(1:selftasksize, 1:nselftasks+1))

        ! Mark the workpackage as initialized
        wptr%is_init = .TRUE.
    END SUBROUTINE recording_pass


    SUBROUTINE init_conn2()
        ! Subroutine arguments
        ! none...

        ! Local variables
        INTEGER(intk) :: nselfsend, nselfrecv

        IF (is_init) CALL errr(__FILE__, __LINE__)

        ! The maximum number of concurrent communications are the number
        ! of processes
        ALLOCATE(recvidxlist(3, SIZE(recvconns, 2)), source=0_intk)
        ALLOCATE(sendlist(numprocs), source=0_int32)
        ALLOCATE(recvlist(numprocs), source=0_int32)
        ALLOCATE(sendreqs(numprocs))
        ALLOCATE(recvreqs(numprocs))
        nrecv = 0
        nsend = 0

        ! Allocating the workpackage arrays
        ! Always add 1 extra dummy task for error checking purposes

        ! 6 Variables may be exchanged per send/receive connection
        ! Upper bound: 6*irecv+1 or 6*isend+1
        maxsendtasks = 6 * isend + 1
        maxrecvtasks = 6 * irecv + 1
        ALLOCATE(sendtasks(buffertasksize, maxsendtasks))
        ALLOCATE(recvtasks(buffertasksize, maxrecvtasks))
        !$omp target enter data map(always, to: sendtasks, recvtasks)

        ! One grid has up to 26 neighbors that may live on the same rank.
        ! Data may be exchanged in forward and backward direction.
        ! In total, up to 6 variables may be exchanged.
        ! Upper bound: nmygrids*26*2*6+1 (that is a lot of combinations...)
        ! Instead, count the number of selftasks upfront and take the upper
        ! bound for 6 fields
        CALL count_selftasks(nselfsend, nselfrecv)
        maxselftasks = 6 * (nselfsend + nselfrecv) + 1
        ALLOCATE(selftasks(selftasksize, maxselftasks))
        !$omp target enter data map(always, to: selftasks)

        ! MPI tasks are only used on the host and do not exceed isend/irecv+1
        maxmpisendtasks = isend + 1
        maxmpirecvtasks = irecv + 1
        ALLOCATE(mpisendtasks(mpitasksize, maxmpisendtasks))
        ALLOCATE(mpirecvtasks(mpitasksize, maxmpirecvtasks))

        ! Allocate the workrecords array for all possible types of conn
        ! dimension 1 = ilevel (including maxlevel+1 for "all levels")
        ! dimension 2 = layers (1 or 2)
        ! dimension 3 = argument combination (6 bits = 2^6-1 = 63)
        ! dimension 4 = corners (0 or 1)
        ! dimension 5 = normal (0 or 1)
        ! dimension 6 = forward (-1, 0, 1)
        ALLOCATE(workrecords(minlevel:maxlevel+1, 1:2, 1:63, 0:1, 0:1, -1:1))

        is_init = .TRUE.

        ! Recording relevant conn calls for later efficient reusage on device
        CALL run_recording_pass()
    END SUBROUTINE init_conn2


    SUBROUTINE count_selftasks(nselfsend, nselfrecv)
        ! Subroutine arguments
        INTEGER(intk), INTENT(out) :: nselfsend, nselfrecv

        ! Local variables
        INTEGER(intk) :: i, iprocnbr

        nselfsend = 0
        DO i = 1, isend
            iprocnbr = sendconns(1, i)
            IF (iprocnbr == myid) THEN
                nselfsend = nselfsend + 1
            END IF
        END DO

        nselfrecv = 0
        DO i = 1, irecv
            iprocnbr = recvconns(2, i)
            IF (iprocnbr == myid) THEN
                nselfrecv = nselfrecv + 1
            END IF
        END DO

        IF (nselfsend /= nselfrecv) THEN
            CALL errr(__FILE__, __LINE__)
        END IF
    END SUBROUTINE count_selftasks


    SUBROUTINE run_recording_pass()
        ! Subroutine arguments
        ! none...

        ! Local variables
        TYPE(field_t) :: pdummy, udummy, vdummy, wdummy
        INTEGER(intk) :: ilevel

        is_recording = .TRUE.

        CALL pdummy%init("DUMMY")
        CALL udummy%init("DUMMY", istag=1)
        CALL vdummy%init("DUMMY", jstag=1)
        CALL wdummy%init("DUMMY", kstag=1)

        !$omp target enter data map(to: &
        !$omp&  udummy%arr, vdummy%arr, wdummy%arr, pdummy%arr)

        ! START -- This section defines the recored variants of conn2 ---

        DO ilevel = minlevel, maxlevel
            ! Inner pressuresolver iterations
            CALL conn2(ilevel, layers=1, s1=pdummy)
            CALL conn2(ilevel, layers=1, s1=pdummy, forward=-1)
        END DO

        ! Outer pressuresolver iterations
        CALL conn2(layers=1, s1=pdummy)

        ! END -- This section defines the recored variants of conn2 ---

        !$omp target exit data map(delete: &
        !$omp&  udummy%arr, vdummy%arr, wdummy%arr, pdummy%arr)

        CALL pdummy%finish()
        CALL udummy%finish()
        CALL vdummy%finish()
        CALL wdummy%finish()

        is_recording = .FALSE.
    END SUBROUTINE run_recording_pass


    SUBROUTINE finish_conn2()
        ! Subroutine arguments
        ! none...

        ! Local variables
        ! none...

        IF (.NOT. is_init) CALL errr(__FILE__, __LINE__)

        DEALLOCATE(recvidxlist)
        DEALLOCATE(sendlist)
        DEALLOCATE(recvlist)
        DEALLOCATE(sendreqs)
        DEALLOCATE(recvreqs)

        !$omp target exit data map(always, delete: sendtasks, recvtasks, &
        !$omp& selftasks)
        DEALLOCATE(sendtasks)
        DEALLOCATE(recvtasks)
        DEALLOCATE(selftasks)

        DEALLOCATE(mpisendtasks)
        DEALLOCATE(mpirecvtasks)

        CALL purge_workrecords()
        DEALLOCATE(workrecords)

        is_init = .FALSE.
    END SUBROUTINE finish_conn2


    SUBROUTINE purge_workrecords()
        ! Subroutine arguments
        ! none...

        ! Local variables
        TYPE(work_t), POINTER :: wrptr(:)
        INTEGER(intk) :: i

        wrptr(1:SIZE(workrecords)) => workrecords

        DO i = 1, SIZE(workrecords)
            IF (.NOT. wrptr(i)%is_init) CYCLE

            !$omp target exit data map(delete: wrptr(i)%sendtasks, &
            !$omp& wrptr(i)%recvtasks, wrptr(i)%selftasks)
            DEALLOCATE(wrptr(i)%sendtasks)
            DEALLOCATE(wrptr(i)%recvtasks)
            DEALLOCATE(wrptr(i)%selftasks)
            DEALLOCATE(wrptr(i)%mpisendtasks)
            DEALLOCATE(wrptr(i)%mpirecvtasks)
            wrptr(i)%is_init = .FALSE.
        END DO
    END SUBROUTINE purge_workrecords


    ! Host routine for preparing the sendtasks, selftasks and mpisendtasks
    ! Task parameters are gathered in arrays but tasks are not executed.
    !
    SUBROUTINE prepare_tasks_all(stasks, nstasks, &
            etasks, netasks, mpistasks, nmpistasks, &
            minconlvl, maxconlvl, nplane, vertices, &
            normal, fwd, flag, nvars, v1, v2, v3, s1, s2, s3)

        ! Subroutine arguments
        INTEGER(intk), INTENT(inout) :: stasks(buffertasksize, maxsendtasks)
        INTEGER(intk), INTENT(out) :: nstasks
        INTEGER(intk), INTENT(inout) :: etasks(selftasksize, maxselftasks)
        INTEGER(intk), INTENT(out) :: netasks
        INTEGER(intk), INTENT(inout) :: mpistasks(mpitasksize, maxmpisendtasks)
        INTEGER(intk), INTENT(out) :: nmpistasks

        INTEGER(intk), INTENT(in) :: minconlvl, maxconlvl, nplane
        LOGICAL, INTENT(in) :: vertices, normal
        INTEGER(intk), INTENT(in) :: fwd
        CHARACTER(len=1), INTENT(in) :: flag
        INTEGER(intk), INTENT(in) :: nvars
        TYPE(field_t), OPTIONAL, INTENT(inout) :: v1, v2, v3, s1, s2, s3

        ! Local variables
        INTEGER(intk) :: i, iprocnbr, igrid, ifacerecv, facearea
        INTEGER(intk) :: istask, ietask, impistasks
        LOGICAL :: exchange, geometry
        INTEGER(int32) :: sendcounter, messagelength

#ifdef _MGLET_PROFILE_ANNOTATIONS_
        CALL profile_range_push("prepare_tasks_all")
#endif

        geometry = .FALSE.

        ! Pack all buffers and send data
        sendcounter = 0
        messagelength = 0
        nsend = 0

        ! Initializing the task counters to zero
        istask = 0
        ietask = 0
        impistasks = 0

        DO i = 1, SIZE(sendconns, 2)

            ! Just-in-time decision whether to exchange data or not
            exchange = decide(i, sendconns, geometry, vertices, fwd, &
                minconlvl, maxconlvl)
            iprocnbr = sendconns(1, i)

            ! Intra-rank communication with direct copy between local grids
            IF (iprocnbr == myid .AND. exchange) THEN

                ! >>> adding entries to selftasks (= "connect_self")
                CALL prepare_selftasks(etasks, ietask, i, nplane, &
                    normal, flag, v1, v2, v3, s1, s2, s3)
                CYCLE

            END IF

            ! Message is send via MPI and tasks for buffer filling as added
            IF (exchange) THEN

                igrid = sendconns(3, i)
                ifacerecv = sendconns(5, i)
                facearea = face_area(igrid, ifacerecv, nplane, flag)

                ! >>> adding entries to sendtasks (= "write_buffer")
                CALL prepare_sendtasks(stasks, istask, i, messagelength, &
                    sendcounter, nplane, normal, flag, nvars, v1, v2, v3, &
                    s1, s2, s3)

                messagelength = messagelength + nvars*facearea

            END IF

            ! Check if we need to record an MPI task
            ! >>> adding entries to mpisendtasks (= "post_send")
            IF (messagelength > 0) THEN
                IF (i == SIZE(sendconns, 2)) THEN
                    CALL add_mpi_task(mpistasks, impistasks, &
                        iprocnbr, messagelength, sendcounter, "send")
                ELSE IF (sendconns(1, i + 1) /= iprocnbr) THEN
                    CALL add_mpi_task(mpistasks, impistasks, &
                        iprocnbr, messagelength, sendcounter, "send")
                END IF
            END IF

        END DO

        ! Set the output task counters
        nstasks = istask
        netasks = ietask
        nmpistasks = impistasks

        ! Add a harmful dummy task at (ntasks+1) for checking
        stasks(:, nstasks+1) = -1
        etasks(:, netasks+1) = -1
        mpistasks(:, nmpistasks+1) = -1

#ifdef _MGLET_PROFILE_ANNOTATIONS_
        CALL profile_range_pop()
#endif
    END SUBROUTINE prepare_tasks_all


    SUBROUTINE prepare_sendtasks(stasks, istask, sendid, messagelength, &
            sendcounter, nplane, normal, flag, nvars, v1, v2, v3, s1, s2, s3)

        ! Subroutine arguments
        INTEGER(intk), INTENT(inout) :: stasks(buffertasksize, maxsendtasks)
        INTEGER(intk), INTENT(inout) :: istask
        INTEGER(intk), INTENT(in) :: sendid
        INTEGER(intk), INTENT(inout) :: messagelength
        INTEGER(intk), INTENT(inout) :: sendcounter
        INTEGER(intk), INTENT(in) :: nplane
        LOGICAL, INTENT(in) :: normal
        CHARACTER(len=1), INTENT(in) :: flag
        INTEGER(intk), INTENT(in) :: nvars
        TYPE(field_t), OPTIONAL, INTENT(inout) :: v1, v2, v3, s1, s2, s3

        ! Indices of start- and stop of iteration over boundary face
        INTEGER(intk) :: istart, istop, jstart, jstop, kstart, kstop

        ! Grid to send from
        ! Must be intk because it intreface with MGLET
        INTEGER(intk) :: igrid, ifacerecv, ifacesend, fieldid

        ! Message sizes
        ! Must be int32 because it iterface with MPI
        INTEGER(int32) :: thismessagelength, facearea
        INTEGER(int32) :: icount, offset

        ! Flags to indicate exchange of U, V, W
        LOGICAL :: exU, exV, exW, exp1

        ! Set variables from send table
        igrid = sendconns(4, sendid)
        ifacerecv = sendconns(5, sendid)
        ifacesend = sendconns(6, sendid)

        ! Get start- and stop indices of grid
        CALL start_and_stop(igrid, facenbr(ifacerecv), istart, istop, &
            jstart, jstop, kstart, kstop, nplane, flag, nghost=1)

        CALL corr_start_stop(igrid, ifacesend, ifacerecv, &
            istart, istop, jstart, jstop, kstart, kstop, nplane, flag)

        facearea = (istop-istart+1)*(jstop-jstart+1)*(kstop-kstart+1)
        thismessagelength = nvars*facearea

        ! Check that buffer does not overflow
        IF (sendcounter + messagelength + thismessagelength &
                > SIZE(sendbuf)) THEN
            CALL errr(__FILE__, __LINE__)
        END IF

        ! Reset message size counter
        offset = sendcounter + messagelength
        icount = offset

        ! Fill buffers
        IF (flag == 'W') THEN
            exU = (ifacerecv == 1)
        ELSE
            exU = (normal .AND. ifacerecv < 3) .OR. (.NOT. normal)
        END IF
        IF (PRESENT(v1) .AND. exU) THEN
            fieldid = 1   ! (for v1)
            istask = istask + 1
            CALL add_single_task(stasks(:, istask), fieldid, icount, &
                igrid, istart, istop, jstart, jstop, kstart, kstop)
        END IF

        IF (flag == 'W') THEN
            exV = (ifacerecv == 3)
        ELSE
            exV = (normal .AND. (ifacerecv > 2 .AND. ifacerecv < 5)) .OR. &
                (.NOT. normal)
        END IF
        IF (PRESENT(v2) .AND. exV) THEN
            fieldid = 2   ! (for v2)
            istask = istask + 1
            CALL add_single_task(stasks(:, istask), fieldid, icount, &
                igrid, istart, istop, jstart, jstop, kstart, kstop)
        END IF

        IF (flag == 'W') THEN
            exW = (ifacerecv == 5)
        ELSE
            exW = (normal .AND. ifacerecv > 4) .OR. (.NOT. normal)
        END IF
        IF (PRESENT(v3) .AND. exW) THEN
            fieldid = 3   ! (for v3)
            istask = istask + 1
            CALL add_single_task(stasks(:, istask), fieldid, icount, &
                igrid, istart, istop, jstart, jstop, kstart, kstop)
        END IF

        IF (flag == 'W') THEN
            exp1 = (ifacerecv == 2) .OR. (ifacerecv == 4) .OR. &
                (ifacerecv == 6)
        ELSE
            exp1 = .TRUE.
        END IF
        IF (PRESENT(s1) .AND. exp1) THEN
            fieldid = 4   ! (for s1)
            istask = istask + 1
            CALL add_single_task(stasks(:, istask), fieldid, icount, &
                igrid, istart, istop, jstart, jstop, kstart, kstop)
        END IF

        IF (PRESENT(s2)) THEN
            fieldid = 5   ! (for s2)
            istask = istask + 1
            CALL add_single_task(stasks(:, istask), fieldid, icount, &
                igrid, istart, istop, jstart, jstop, kstart, kstop)
        END IF

        IF (PRESENT(s3)) THEN
            fieldid = 6   ! (for s3)
            istask = istask + 1
            CALL add_single_task(stasks(:, istask), fieldid, icount, &
                igrid, istart, istop, jstart, jstop, kstart, kstop)
        END IF

        ! Check that message length was calculated correctly
        IF (thismessagelength /= (icount - offset)) THEN
            write(*, *) "thismessagelength:", thismessagelength, &
                "icount:", icount
            CALL errr(__FILE__, __LINE__)
        END IF

    END SUBROUTINE prepare_sendtasks


    ! Routine to prepare the workpackage for intra-rank self connect
    !
    SUBROUTINE prepare_selftasks(etasks, ietask, sendid, nplane, &
            normal, flag, v1, v2, v3, s1, s2, s3)

        ! Subroutine arguments
        INTEGER(intk), INTENT(inout) :: etasks(selftasksize, maxselftasks)
        INTEGER(intk), INTENT(inout) :: ietask
        INTEGER(int32), INTENT(in) :: sendid
        INTEGER(int32), INTENT(in) :: nplane
        LOGICAL, INTENT(in) :: normal
        CHARACTER(len=1), INTENT(in) :: flag
        TYPE(field_t), INTENT(inout), OPTIONAL :: v1, v2, v3, s1, s2, s3

        ! Indices of start- and stop of iteration over boundary face
        ! Source face
        INTEGER(intk) :: istart, istop, jstart, jstop, kstart, kstop

        ! Indices of start- and stop of iteration over boundary face
        ! Destination face
        INTEGER(intk) :: istart_d, istop_d, jstart_d, jstop_d, kstart_d, kstop_d

        ! Grid to send from
        ! Must be intk because it intreface with MGLET
        INTEGER(intk) :: igrid, igrid_d, ifacerecv, ifacesend, fieldid

        ! Message sizes
        ! Must be int32 because it iterface with MPI
        INTEGER(int32) :: dest_size, source_size

        ! Flags to indicate exchange of U, V, W
        LOGICAL :: exU, exV, exW, exp1

        ! Set variables from send table
        igrid_d = sendconns(3, sendid)
        igrid = sendconns(4, sendid)
        ifacerecv = sendconns(5, sendid)
        ifacesend = sendconns(6, sendid)

        ! Get start- and stop indices of source grid
        CALL start_and_stop(igrid, facenbr(ifacerecv), &
            istart, istop, jstart, jstop, kstart, kstop, nplane, flag, nghost=1)
        CALL corr_start_stop(igrid, ifacesend, ifacerecv, &
            istart, istop, jstart, jstop, kstart, kstop, nplane, flag)

        ! Get start- and stop indices of destination grid
        CALL start_and_stop(igrid_d, ifacerecv, &
            istart_d, istop_d, jstart_d, jstop_d, kstart_d, kstop_d, &
            nplane, flag)

        ! Sanity check of message length
        source_size = (istop-istart+1)*(jstop-jstart+1)*(kstop-kstart+1)
        dest_size = (istop_d-istart_d+1) &
            *(jstop_d-jstart_d+1)*(kstop_d-kstart_d+1)
        IF (source_size /= dest_size) THEN
            CALL errr(__FILE__, __LINE__)
        END IF

        IF (flag == 'W') THEN
            exU = (ifacerecv == 1)
        ELSE
            exU = (normal .AND. ifacerecv < 3) .OR. (.NOT. normal)
        END IF
        IF (PRESENT(v1) .AND. exU) THEN
            fieldid = 1   ! (for v1)
            ietask = ietask + 1
            CALL add_self_task(etasks(:, ietask), fieldid, igrid, &
                igrid_d, istart, istop, jstart, jstop, kstart, kstop, &
                istart_d, istop_d, jstart_d, jstop_d, kstart_d, kstop_d)
        END IF

        IF (flag == 'W') THEN
            exV = (ifacerecv == 3)
        ELSE
            exV = (normal .AND. (ifacerecv > 2 .AND. ifacerecv < 5)) .OR. &
                (.NOT. normal)
        END IF
        IF (PRESENT(v2) .AND. exV) THEN
            fieldid = 2   ! (for v2)
            ietask = ietask + 1
            CALL add_self_task(etasks(:, ietask), fieldid, igrid, &
                igrid_d, istart, istop, jstart, jstop, kstart, kstop, &
                istart_d, istop_d, jstart_d, jstop_d, kstart_d, kstop_d)
        END IF

        IF (flag == 'W') THEN
            exW = (ifacerecv == 5)
        ELSE
            exW = (normal .AND. ifacerecv > 4) .OR. (.NOT. normal)
        END IF
        IF (PRESENT(v3) .AND. exW) THEN
            fieldid = 3   ! (for v3)
            ietask = ietask + 1
            CALL add_self_task(etasks(:, ietask), fieldid, igrid, &
                igrid_d, istart, istop, jstart, jstop, kstart, kstop, &
                istart_d, istop_d, jstart_d, jstop_d, kstart_d, kstop_d)
        END IF

        IF (flag == 'W') THEN
            exp1 = (ifacerecv == 2) .OR. (ifacerecv == 4) .OR. &
                (ifacerecv == 6)
        ELSE
            exp1 = .TRUE.
        END IF
        IF (PRESENT(s1) .AND. exp1) THEN
            fieldid = 4   ! (for s1)
            ietask = ietask + 1
            CALL add_self_task(etasks(:, ietask), fieldid, igrid, &
                igrid_d, istart, istop, jstart, jstop, kstart, kstop, &
                istart_d, istop_d, jstart_d, jstop_d, kstart_d, kstop_d)
        END IF

        IF (PRESENT(s2)) THEN
            fieldid = 5   ! (for s2)
            ietask = ietask + 1
            CALL add_self_task(etasks(:, ietask), fieldid, igrid, &
                igrid_d, istart, istop, jstart, jstop, kstart, kstop, &
                istart_d, istop_d, jstart_d, jstop_d, kstart_d, kstop_d)
        END IF

        IF (PRESENT(s3)) THEN
            fieldid = 6   ! (for s3)
            ietask = ietask + 1
            CALL add_self_task(etasks(:, ietask), fieldid, igrid, &
                igrid_d, istart, istop, jstart, jstop, kstart, kstop, &
                istart_d, istop_d, jstart_d, jstop_d, kstart_d, kstop_d)
        END IF
    END SUBROUTINE prepare_selftasks


    ! Store the paremters for a single non-blocking MPI send call
    !
    SUBROUTINE add_mpi_task(mpitasks, impitask, iprocnbr, &
        messagelength, counter, type)

        ! Subroutine arguments
        INTEGER(intk), INTENT(inout) :: mpitasks(mpitasksize, *)
        INTEGER(intk), INTENT(inout) :: impitask
        INTEGER(intk), INTENT(in) :: iprocnbr
        INTEGER(intk), INTENT(inout) :: messagelength
        INTEGER(intk), INTENT(inout) :: counter
        CHARACTER(len=4), INTENT(in) :: type

        ! Local variables
        ! none...

        IF (type == 'send') THEN
            nsend = nsend + 1
            sendlist(nsend) = iprocnbr
        ELSE IF (type == 'recv') THEN
            nrecv = nrecv + 1
            recvlist(nrecv) = iprocnbr
        ELSE
            CALL errr(__FILE__, __LINE__)
        END IF

        ! Add the MPI send task to the mpitasks array
        impitask = impitask + 1

        mpitasks(1, impitask) = iprocnbr
        mpitasks(2, impitask) = messagelength
        mpitasks(3, impitask) = counter

        counter = counter + messagelength
        messagelength = 0

    END SUBROUTINE add_mpi_task


    ! Routine to prepare all workpackages for processing the receive buffer
    !
    SUBROUTINE prepare_recvtasks_all(rtasks, nrtasks, nplane, normal, &
            flag, v1, v2, v3, s1, s2, s3)

        ! Subroutine arguments
        INTEGER(intk), INTENT(inout) :: rtasks(buffertasksize, maxrecvtasks)
        INTEGER(intk), INTENT(out) :: nrtasks
        INTEGER(int32), INTENT(in) :: nplane
        LOGICAL, INTENT(in) :: normal
        CHARACTER(len=1), INTENT(in) :: flag
        TYPE(field_t), OPTIONAL, INTENT(inout) :: &
            v1, v2, v3, s1, s2, s3

        ! Local variables
        INTEGER(int32) :: idx, i
        TYPE(MPI_Status) :: recvstatus
        INTEGER(int32) :: recvmessagelen
        INTEGER(int32) :: unpacklen
        INTEGER(intk) :: irtask

#ifdef _MGLET_PROFILE_ANNOTATIONS_
        CALL profile_range_push("prepare_recvtasks_all")
#endif

        irtask = 0

        DO WHILE (.TRUE.)
            IF (nrecv == 0) EXIT

            ! Get index of an already arrived message
            CALL MPI_Waitany(nrecv, recvreqs, idx, recvstatus)

            IF (idx /= MPI_UNDEFINED) THEN
                CALL MPI_Get_count(recvstatus, mglet_mpi_real, &
                    recvmessagelen)

                unpacklen = 0
                DO i = 1, SIZE(recvidxlist, 2)
                    IF (recvidxlist(1, i) == recvlist(idx) &
                            .AND. recvidxlist(2, i) > 0) THEN

                        ! >>> adding entries to recvtasks ("read_buffer")
                        CALL prepare_recvtasks(rtasks, irtask, &
                            i, nplane, normal, flag, v1, v2, v3, s1, s2, s3)

                        unpacklen = unpacklen + recvidxlist(2, i)
                    END IF
                END DO

                ! Security check -- not possible before message arrival
                IF (recvmessagelen /= unpacklen) THEN
                    CALL errr(__FILE__, __LINE__)
                END IF
            ELSE
                EXIT
            END IF
        END DO

        ! Enfore the completion of all send requests
        CALL MPI_Waitall(nsend, sendreqs, MPI_STATUSES_IGNORE)

        ! Setting the number of recvtasks to the actual number of tasks
        nrtasks = irtask

        ! Add a harmful dummy task at (ntasks+1) to detect execution overshoot
        rtasks(:, nrtasks+1) = -1

#ifdef _MGLET_PROFILE_ANNOTATIONS_
        CALL profile_range_pop()
#endif
    END SUBROUTINE prepare_recvtasks_all


    ! Routine to prepare a workpackage for reading received values from buffer
    !
    SUBROUTINE prepare_recvtasks(rtasks, irtask, recvid, nplane, &
            normal, flag, v1, v2, v3, s1, s2, s3)

        ! Subroutine arguments
        INTEGER(intk), INTENT(inout) :: rtasks(buffertasksize, maxrecvtasks)
        INTEGER(intk), INTENT(inout) :: irtask
        INTEGER(intk), INTENT(in) :: recvid
        INTEGER(int32), INTENT(in) :: nplane
        LOGICAL, INTENT(in) :: normal
        CHARACTER(len=1), INTENT(in) :: flag
        TYPE(field_t), OPTIONAL, INTENT(inout) :: v1, v2, v3, s1, s2, s3

        ! Indices of start- and stop of iteration over boundary face
        INTEGER(intk) :: istart, istop, jstart, jstop, kstart, kstop

        ! Grid to send from
        ! Must be intk because it intreface with MGLET
        INTEGER(intk) :: igrid, ifacerecv, fieldid

        ! Message sizes
        ! Must be int32 because it iterface with MPI
        INTEGER(int32) :: offset, icount

        ! Flags to indicate exchange of U, V, W
        LOGICAL :: exU, exV, exW, exp1


        ! Array for the workpackage:
        ! -> wp(10, SIZE(recvconns, 2))

        ! Set variables from send table
        igrid = recvconns(3, recvid)
        ifacerecv = recvconns(5, recvid)

        ! ioperation = 0

        ! Get start- and stop indices of grid
        CALL start_and_stop(igrid, ifacerecv, &
            istart, istop, jstart, jstop, kstart, kstop, nplane, flag)

        ! Zeroise message size counter
        offset = recvidxlist(3, recvid)
        icount = offset

        IF (flag == 'W') THEN
            exU = (ifacerecv == 1)
        ELSE
            exU = (normal .AND. ifacerecv < 3) .OR. (.NOT. normal)
        END IF
        IF (PRESENT(v1) .AND. exU) THEN
            irtask = irtask + 1
            fieldid = 1   ! (for v1)
            CALL add_single_task(rtasks(:, irtask), fieldid, icount, &
                igrid, istart, istop, jstart, jstop, kstart, kstop)
        END IF

        IF (flag == 'W') THEN
            exV = (ifacerecv == 3)
        ELSE
            exV = (normal .AND. (ifacerecv > 2 .AND. ifacerecv < 5)) .OR. &
                (.NOT. normal)
        END IF
        IF (PRESENT(v2) .AND. exV) THEN
            irtask = irtask + 1
            fieldid = 2   ! (for v2)
            CALL add_single_task(rtasks(:, irtask), fieldid, icount, &
                igrid, istart, istop, jstart, jstop, kstart, kstop)
        END IF

        IF (flag == 'W') THEN
            exW = (ifacerecv == 5)
        ELSE
            exW = (normal .AND. ifacerecv > 4) .OR. (.NOT. normal)
        END IF
        IF (PRESENT(v3) .AND. exW) THEN
            irtask = irtask + 1
            fieldid = 3   ! (for v3)
            CALL add_single_task(rtasks(:, irtask), fieldid, icount, &
                igrid, istart, istop, jstart, jstop, kstart, kstop)
        END IF

        IF (flag == 'W') THEN
            exp1 = (ifacerecv == 2) .OR. (ifacerecv == 4) .OR. &
                (ifacerecv == 6)
        ELSE
            exp1 = .TRUE.
        END IF
        IF (PRESENT(s1) .AND. exp1) THEN
            irtask = irtask + 1
            fieldid = 4   ! (for s1)
            CALL add_single_task(rtasks(:, irtask), fieldid, icount, &
                igrid, istart, istop, jstart, jstop, kstart, kstop)
        END IF

        IF (PRESENT(s2)) THEN
            irtask = irtask + 1
            fieldid = 5   ! (for s2)
            CALL add_single_task(rtasks(:, irtask), fieldid, icount, &
                igrid, istart, istop, jstart, jstop, kstart, kstop)
        END IF

        IF (PRESENT(s3)) THEN
            irtask = irtask + 1
            fieldid = 6   ! (for s3)
            CALL add_single_task(rtasks(:, irtask), fieldid, icount, &
                igrid, istart, istop, jstart, jstop, kstart, kstop)
        END IF

        ! Check that message length is calculated correctly
        IF ((icount - offset) /= recvidxlist(2, recvid)) THEN
            WRITE(*, *) "icount:", icount, &
                "recvidxlist(2, recvid):", recvidxlist(2, recvid)
            CALL errr(__FILE__, __LINE__)
        END IF

    END SUBROUTINE prepare_recvtasks




    ! Core routine to add a single self-copy task to the workpackage
    !
    PURE SUBROUTINE add_self_task(selftask, fieldid, igrid, igrid_d, &
            istart, istop, jstart, jstop, kstart, kstop, istart_d, &
            istop_d, jstart_d, jstop_d, kstart_d, kstop_d)

        ! Subroutine arguments
        INTEGER(intk), INTENT(inout) :: selftask(selftasksize)
        INTEGER(intk), INTENT(in) :: fieldid, igrid, igrid_d, istart, istop, &
            jstart, jstop, kstart, kstop, istart_d, istop_d, jstart_d, &
            jstop_d, kstart_d, kstop_d

        ! Filling the task
        selftask(1) = fieldid
        selftask(2) = igrid
        selftask(3) = igrid_d
        ! for the source grid (igrid) = for packing
        selftask(4) = istart
        selftask(5) = istop
        selftask(6) = jstart
        selftask(7) = jstop
        selftask(8) = kstart
        selftask(9) = kstop
        ! for the destination grid (igrid_d) = for unpacking
        selftask(10) = istart_d
        selftask(11) = istop_d
        selftask(12) = jstart_d
        selftask(13) = jstop_d
        selftask(14) = kstart_d
        selftask(15) = kstop_d

    END SUBROUTINE add_self_task


    ! Core routine to add a single buffer-related task to the workpackage
    !
    PURE SUBROUTINE add_single_task(task, fieldid, icount, igrid, &
            istart, istop, jstart, jstop, kstart, kstop)

        ! Subroutine arguments
        INTEGER(intk), INTENT(out) :: task(9)
        INTEGER(intk), INTENT(in) :: fieldid
        INTEGER(intk), INTENT(inout) :: icount
        INTEGER(intk), INTENT(in) :: igrid, istart, istop, jstart, jstop, &
            kstart, kstop

        ! Local variables
        INTEGER(intk) :: tasksize

        ! Filling the task
        task(1) = fieldid
        task(2) = icount
        task(3) = igrid
        task(4) = istart
        task(5) = istop
        task(6) = jstart
        task(7) = jstop
        task(8) = kstart
        task(9) = kstop

        ! Increment the icount by the number of elements in the task
        tasksize = (istop-istart+1)*(jstop-jstart+1)*(kstop-kstart+1)
        icount = icount + tasksize

    END SUBROUTINE add_single_task




    ! Now the routines with offloaded kernels

    ! Routine with offloaded kernel to process all sendtasks on the device
    !
    SUBROUTINE process_sendtasks(nstasks, stasks)

        ! Subroutine arguments
        INTEGER(intk), INTENT(in) :: nstasks
        INTEGER(intk), INTENT(in) :: stasks(buffertasksize, nstasks)

        ! Local variables
        INTEGER(intk) :: itask, fieldid, icount, igrid, istart, istop, &
            jstart, jstop, kstart, kstop, ii, jj, kk, ip3

        IF (nstasks == 0) THEN
            RETURN
        END IF

        ! At all cost, avoid pointers within the kernel or, even worse,
        ! pointer operations within the kernel!
        ASSOCIATE(a1 => f1%arr, a2 => f2%arr, a3 => f3%arr, &
                  a4 => f4%arr, a5 => f5%arr, a6 => f6%arr)

#ifdef _MGLET_PROFILE_ANNOTATIONS_
        CALL profile_range_push("process_sendtasks")
#endif
        !$omp target teams distribute private(itask, fieldid, icount, &
        !$omp&  igrid, istart, istop, jstart, jstop, kstart, kstop, &
        !$omp&  ii, jj, kk, ip3)
        DO itask = 1, nstasks

            ! Set variables from sendtasks workpackage
            fieldid = stasks(1, itask)
            icount = stasks(2, itask)
            igrid = stasks(3, itask)
            istart = stasks(4, itask)
            istop = stasks(5, itask)
            jstart = stasks(6, itask)
            jstop = stasks(7, itask)
            kstart = stasks(8, itask)
            kstop = stasks(9, itask)

            ! The following replaces "read_single_buffer"
            CALL get_ip3(ip3, igrid)
            CALL get_mgdims(kk, jj, ii, igrid)

            ! Assign the correct field pointer based on ifield
            SELECT CASE (fieldid)
            CASE (1)
                CALL arr_to_buf(kk, jj, ii, a1(ip3), istart, istop, &
                    jstart, jstop, kstart, kstop, icount)
            CASE (2)
                CALL arr_to_buf(kk, jj, ii, a2(ip3), istart, istop, &
                    jstart, jstop, kstart, kstop, icount)
            CASE (3)
                CALL arr_to_buf(kk, jj, ii, a3(ip3), istart, istop, &
                    jstart, jstop, kstart, kstop, icount)
            CASE (4)
                CALL arr_to_buf(kk, jj, ii, a4(ip3), istart, istop, &
                    jstart, jstop, kstart, kstop, icount)
            CASE (5)
                CALL arr_to_buf(kk, jj, ii, a5(ip3), istart, istop, &
                    jstart, jstop, kstart, kstop, icount)
            CASE (6)
                CALL arr_to_buf(kk, jj, ii, a6(ip3), istart, istop, &
                    jstart, jstop, kstart, kstop, icount)
            CASE DEFAULT
                CALL errr(__FILE__, __LINE__)
            END SELECT

        END DO
        !$omp end target teams distribute
#ifdef _MGLET_PROFILE_ANNOTATIONS_
        CALL profile_range_pop()
#endif
        END ASSOCIATE
    END SUBROUTINE process_sendtasks


    SUBROUTINE arr_to_buf(kk, jj, ii, arr, istart, istop, &
        jstart, jstop, kstart, kstop, icount)
        !$omp declare target
        ! Subroutine arguments
        INTEGER(intk), INTENT(in) :: kk, jj, ii
        REAL(realk), INTENT(in) :: arr(kk, jj, ii)
        INTEGER(intk), INTENT(in) :: istart, istop, jstart, jstop, kstart, &
            kstop, icount
        ! Local variables
        INTEGER(intk) :: i, j, k, idx_b, kkl, jjl

        kkl = kstop - kstart + 1
        jjl = jstop - jstart + 1

        !$omp parallel do collapse(3) private(i, j, k, idx_b)
        DO i = istart, istop
            DO j = jstart, jstop
                DO k = kstart, kstop
                    idx_b = 1 + (k-kstart) + (j-jstart)*kkl + &
                        (i-istart)*jjl*kkl + icount
                    sendbuf(idx_b) = arr(k, j, i)
                END DO
            END DO
        END DO
        !$omp end parallel do

    END SUBROUTINE arr_to_buf


    ! Routine with offloaded kernel to process all receive tasks on the device
    !
    SUBROUTINE process_recvtasks(nrtasks, rtasks)

        ! Subroutine arguments
        INTEGER(intk), INTENT(in) :: nrtasks
        INTEGER(intk), INTENT(in) :: rtasks(buffertasksize, nrtasks)

        ! Local variables
        INTEGER(intk) :: itask, fieldid, icount, igrid, istart, istop, &
            jstart, jstop, kstart, kstop, ii, jj, kk, ip3

        IF (nrtasks == 0) THEN
            RETURN
        END IF

        ! At all cost, avoid pointers within the kernel or, even worse,
        ! pointer operations within the kernel!
        ASSOCIATE(a1 => f1%arr, a2 => f2%arr, a3 => f3%arr, &
                  a4 => f4%arr, a5 => f5%arr, a6 => f6%arr)


#ifdef _MGLET_PROFILE_ANNOTATIONS_
        CALL profile_range_push("process_recvtasks")
#endif

        !$omp target teams distribute private(itask, fieldid, icount, &
        !$omp&  igrid, istart, istop, jstart, jstop, kstart, kstop, &
        !$omp&  ii, jj, kk, ip3)
        DO itask = 1, nrtasks

            ! Set variables from recvtasks workpackage
            fieldid = rtasks(1, itask)
            icount  = rtasks(2, itask)
            igrid   = rtasks(3, itask)
            istart  = rtasks(4, itask)
            istop   = rtasks(5, itask)
            jstart  = rtasks(6, itask)
            jstop   = rtasks(7, itask)
            kstart  = rtasks(8, itask)
            kstop   = rtasks(9, itask)

            ! The following replaces "read_single_buffer"
            CALL get_ip3(ip3, igrid)
            CALL get_mgdims(kk, jj, ii, igrid)

            ! Assign the correct field pointer
            SELECT CASE (fieldid)
            CASE (1)
                CALL buf_to_arr(kk, jj, ii, a1(ip3), istart, istop, &
                    jstart, jstop, kstart, kstop, icount)
            CASE (2)
                CALL buf_to_arr(kk, jj, ii, a2(ip3), istart, istop, &
                    jstart, jstop, kstart, kstop, icount)
            CASE (3)
                CALL buf_to_arr(kk, jj, ii, a3(ip3), istart, istop, &
                    jstart, jstop, kstart, kstop, icount)
            CASE (4)
                CALL buf_to_arr(kk, jj, ii, a4(ip3), istart, istop, &
                    jstart, jstop, kstart, kstop, icount)
            CASE (5)
                CALL buf_to_arr(kk, jj, ii, a5(ip3), istart, istop, &
                    jstart, jstop, kstart, kstop, icount)
            CASE (6)
                CALL buf_to_arr(kk, jj, ii, a6(ip3), istart, istop, &
                    jstart, jstop, kstart, kstop, icount)
            CASE DEFAULT
                CALL errr(__FILE__, __LINE__)
            END SELECT

        END DO
        !$omp end target teams distribute
#ifdef _MGLET_PROFILE_ANNOTATIONS_
        CALL profile_range_pop()
#endif
        END ASSOCIATE
    END SUBROUTINE process_recvtasks


    SUBROUTINE buf_to_arr(kk, jj, ii, arr, istart, istop, &
        jstart, jstop, kstart, kstop, icount)
        !$omp declare target
        ! Subroutine arguments
        INTEGER(intk), INTENT(in) :: kk, jj, ii
        REAL(realk), INTENT(inout) :: arr(kk, jj, ii)
        INTEGER(intk), INTENT(in) :: istart, istop, jstart, jstop, kstart, &
            kstop, icount
        ! Local variables
        INTEGER(intk) :: i, j, k, idx_b, kkl, jjl

        kkl = kstop - kstart + 1
        jjl = jstop - jstart + 1

        !$omp parallel do collapse(3) private(i, j, k, idx_b)
        DO i = istart, istop
            DO j = jstart, jstop
                DO k = kstart, kstop
                    idx_b = 1 + (k-kstart) + (j-jstart)*kkl + &
                        (i-istart)*jjl*kkl + icount
                    arr(k, j, i) = recvbuf(idx_b)
                END DO
            END DO
        END DO
        !$omp end parallel do

    END SUBROUTINE buf_to_arr


    ! Routine with offloaded kernel to process all selftasks on the device
    !
    SUBROUTINE process_selftasks(nstasks, stasks)

        ! Subroutine arguments
        INTEGER(intk), INTENT(in) :: nstasks
        INTEGER(intk), INTENT(in) :: stasks(selftasksize, nstasks)

        ! Local variables
        INTEGER(intk) :: itask, fieldid, igrid, igrid_d, ip3, ip3_d, &
            kk, jj, ii, istart, istop, jstart, jstop, kstart, kstop, &
            istart_d, istop_d, jstart_d, jstop_d, kstart_d, kstop_d

        ! Precheck to avoid kernel launch overhead
        IF (nstasks == 0) THEN
            RETURN
        END IF

        ! At all cost, avoid pointers within the kernel or, even worse,
        ! pointer operations within the kernel!
        ASSOCIATE(a1 => f1%arr, a2 => f2%arr, a3 => f3%arr, &
                  a4 => f4%arr, a5 => f5%arr, a6 => f6%arr)

#ifdef _MGLET_PROFILE_ANNOTATIONS_
        CALL profile_range_push("process_selftasks")
#endif

        !$omp target teams distribute private(itask, fieldid, igrid, igrid_d, &
        !$omp&  istart, istop, jstart, jstop, kstart, kstop, &
        !$omp&  istart_d, istop_d, jstart_d, jstop_d, kstart_d, kstop_d, &
        !$omp&  ip3, ip3_d, kk, jj, ii)
        DO itask = 1, nstasks

            ! Set variables from selftasks workpackage
            fieldid  = stasks(1, itask)
            igrid    = stasks(2, itask)
            igrid_d  = stasks(3, itask)
            istart   = stasks(4, itask)
            istop    = stasks(5, itask)
            jstart   = stasks(6, itask)
            jstop    = stasks(7, itask)
            kstart   = stasks(8, itask)
            kstop    = stasks(9, itask)
            istart_d = stasks(10, itask)
            istop_d  = stasks(11, itask)
            jstart_d = stasks(12, itask)
            jstop_d  = stasks(13, itask)
            kstart_d = stasks(14, itask)
            kstop_d  = stasks(15, itask)

            ! The following replaces "read_single_buffer"
            CALL get_ip3(ip3, igrid)
            CALL get_ip3(ip3_d, igrid_d)
            CALL get_mgdims(kk, jj, ii, igrid)

            SELECT CASE (fieldid)
            CASE (1)
                CALL arr_to_arr(kk, jj, ii, a1(ip3_d), a1(ip3), &
                    istart, istop, jstart, jstop, kstart, kstop, &
                    istart_d, istop_d, jstart_d, jstop_d, kstart_d, kstop_d)
            CASE (2)
                CALL arr_to_arr(kk, jj, ii, a2(ip3_d), a2(ip3), &
                    istart, istop, jstart, jstop, kstart, kstop, &
                    istart_d, istop_d, jstart_d, jstop_d, kstart_d, kstop_d)
            CASE (3)
                CALL arr_to_arr(kk, jj, ii, a3(ip3_d), a3(ip3), &
                    istart, istop, jstart, jstop, kstart, kstop, &
                    istart_d, istop_d, jstart_d, jstop_d, kstart_d, kstop_d)
            CASE (4)
                CALL arr_to_arr(kk, jj, ii, a4(ip3_d), a4(ip3), &
                    istart, istop, jstart, jstop, kstart, kstop, &
                    istart_d, istop_d, jstart_d, jstop_d, kstart_d, kstop_d)
            CASE (5)
                CALL arr_to_arr(kk, jj, ii, a5(ip3_d), a5(ip3), &
                    istart, istop, jstart, jstop, kstart, kstop, &
                    istart_d, istop_d, jstart_d, jstop_d, kstart_d, kstop_d)
            CASE (6)
                CALL arr_to_arr(kk, jj, ii, a6(ip3_d), a6(ip3), &
                    istart, istop, jstart, jstop, kstart, kstop, &
                    istart_d, istop_d, jstart_d, jstop_d, kstart_d, kstop_d)
            CASE DEFAULT
                CALL errr(__FILE__, __LINE__)
            END SELECT

        END DO
        !$omp end target teams distribute
#ifdef _MGLET_PROFILE_ANNOTATIONS_
        CALL profile_range_pop()
#endif
        END ASSOCIATE
    END SUBROUTINE process_selftasks


    PURE SUBROUTINE arr_to_arr(kk, jj, ii, dst_rarr, src_rarr, &
            istart, istop, jstart, jstop, kstart, kstop, &
            istart_d, istop_d, jstart_d, jstop_d, kstart_d, kstop_d)
        !$omp declare target
        ! Subroutine arguments
        INTEGER(intk), INTENT(in) :: kk, jj, ii
        REAL(realk), INTENT(inout) :: dst_rarr(kk, jj, ii)
        REAL(realk), INTENT(in) :: src_rarr(kk, jj, ii)
        INTEGER(intk), INTENT(in) :: istart, istop, jstart, jstop, kstart, &
             kstop, istart_d, istop_d, jstart_d, jstop_d, kstart_d, kstop_d
        ! Local variables
        INTEGER(intk) :: koff, joff, ioff, i, j, k

        koff = kstart - kstart_d
        joff = jstart - jstart_d
        ioff = istart - istart_d

        !$omp parallel do collapse(3) private(i, j, k)
        DO i = istart_d, istop_d
            DO j = jstart_d, jstop_d
                DO k = kstart_d, kstop_d
                    dst_rarr(k, j, i) = &
                        src_rarr(k + koff, j + joff, i + ioff)
                END DO
            END DO
        END DO
        !$omp end parallel do

    END SUBROUTINE arr_to_arr


    ! Now the routines which which launch non-blocking MPI calls

    ! Perform a single MPI Recv
    !
    SUBROUTINE post_recv(iprocnbr, messagelength, recvcounter)
        ! Subroutine arguments
        INTEGER(int32), INTENT(in) :: iprocnbr
        INTEGER(int32), INTENT(inout) :: messagelength
        INTEGER(int32), INTENT(inout) :: recvcounter

        nrecv = nrecv + 1
        recvlist(nrecv) = iprocnbr

        IF (recvcounter + messagelength > SIZE(recvbuf)) THEN
            CALL errr(__FILE__, __LINE__)
        END IF

        !$omp target data use_device_addr(recvbuf)
        CALL MPI_Irecv(recvbuf(recvcounter+1), messagelength, &
            mglet_mpi_real, iprocnbr, 1, MPI_COMM_WORLD, recvreqs(nrecv))
        !$omp end target data

        recvcounter = recvcounter + messagelength
        messagelength = 0
    END SUBROUTINE post_recv


    ! Perform all non-blocking MPI Send calls based on the mpisendtasks
    !
    SUBROUTINE process_mpisend(nmpistasks, mpistasks)

        ! Subroutine arguments
        INTEGER(intk), INTENT(in) :: nmpistasks
        INTEGER(intk), INTENT(in) :: mpistasks(mpitasksize, nmpistasks+1)

        ! Local variables
        INTEGER(int32) :: itask, iprocnbr, messagelength, sendcounter

#ifdef _MGLET_PROFILE_ANNOTATIONS_
        CALL profile_range_push("process_mpisend")
#endif

        ! Iterate over task and post all non-blocking MPI send calls
        DO itask = 1, nmpistasks

            iprocnbr      = INT(mpistasks(1, itask), int32)
            messagelength = INT(mpistasks(2, itask), int32)
            sendcounter   = INT(mpistasks(3, itask), int32)

            !$omp target data use_device_addr(sendbuf)
            CALL MPI_Isend(sendbuf(sendcounter + 1), messagelength, &
                mglet_mpi_real, iprocnbr, 1, MPI_COMM_WORLD, &
                sendreqs(itask))
            !$omp end target data
        END DO

        nsend = nmpistasks

        ! Safety check based on final dummy entry
        IF (nmpistasks < maxmpisendtasks) THEN
            IF (.NOT. ALL(mpistasks(:, nmpistasks+1) == -1)) THEN
                WRITE(*, *) "Did not encounter the expected dummy task."
                CALL errr(__FILE__, __LINE__)
            END IF
        END IF

#ifdef _MGLET_PROFILE_ANNOTATIONS_
        CALL profile_range_pop()
#endif
    END SUBROUTINE process_mpisend


    ! Perform all non-blocking MPI Send calls based on the mpisendtasks
    !
    SUBROUTINE process_mpirecv(nmpirtasks, mpirtasks)

        ! Subroutine arguments
        INTEGER(intk), INTENT(in) :: nmpirtasks
        INTEGER(intk), INTENT(in) :: mpirtasks(mpitasksize, nmpirtasks+1)

        ! Local variables
        INTEGER(int32) :: itask, iprocnbr, messagelength, recvcounter

#ifdef _MGLET_PROFILE_ANNOTATIONS_
        CALL profile_range_push("process_mpirecv")
#endif

        ! Iterate over task and post all non-blocking MPI send calls
        DO itask = 1, nmpirtasks

            iprocnbr      = INT(mpirtasks(1, itask), int32)
            messagelength = INT(mpirtasks(2, itask), int32)
            recvcounter   = INT(mpirtasks(3, itask), int32)

            !$omp target data use_device_addr(recvbuf)
            CALL MPI_Irecv(recvbuf(recvcounter+1), messagelength, &
                mglet_mpi_real, iprocnbr, 1, MPI_COMM_WORLD, recvreqs(itask))
            !$omp end target data

        END DO

        nrecv = nmpirtasks

        ! Safety check based on final dummy entry
        IF (nmpirtasks < maxmpirecvtasks) THEN
            IF (.NOT. ALL(mpirtasks(:, nmpirtasks+1) == -1)) THEN
                WRITE(*, *) "Did not encounter the expected dummy task."
                CALL errr(__FILE__, __LINE__)
            END IF
        END IF

#ifdef _MGLET_PROFILE_ANNOTATIONS_
        CALL profile_range_pop()
#endif
    END SUBROUTINE process_mpirecv


    ! Perform all non-blocking MPI Recv calls
    !
    SUBROUTINE prepare_mpirecvtasks(mpirtasks, nmpirtasks, &
            minconlvl, maxconlvl, nplane, vertices, &
            normal, fwd, flag, nvars)

        ! Subroutine arguments
        INTEGER(intk), INTENT(inout) :: mpirtasks(mpitasksize, maxmpirecvtasks)
        INTEGER(intk), INTENT(out) :: nmpirtasks
        INTEGER(intk), INTENT(in) :: minconlvl, maxconlvl, nplane
        LOGICAL, INTENT(in) :: vertices, normal
        INTEGER(intk), INTENT(in) :: fwd
        CHARACTER(len=1), INTENT(in) :: flag
        INTEGER(intk), INTENT(in) :: nvars


        ! Local variables
        INTEGER(intk) :: i, iprocnbr, igrid, iface, facearea, impirtask
        LOGICAL :: exchange, geometry
        INTEGER(int32) :: recvcounter, messagelength

        geometry = .FALSE.

        ! Post all receive calls
        recvcounter = 0
        messagelength = 0
        nrecv = 0
        recvidxlist = 0
        impirtask = 0

        DO i = 1, SIZE(recvconns, 2)
            exchange = decide(i, recvconns, geometry, vertices, fwd, &
                minconlvl, maxconlvl)
            iprocnbr = recvconns(2, i)

            ! Communication with self is handled specially in
            ! send_all - nothing to do here
            IF (iprocnbr == myid) THEN
                CYCLE
            END IF

            IF (exchange) THEN
                igrid = recvconns(3, i)
                iface = recvconns(5, i)

                facearea = face_area(igrid, iface, nplane, flag)
                recvidxlist(1, i) = iprocnbr
                recvidxlist(2, i) = nvars*facearea
                recvidxlist(3, i) = recvcounter + messagelength
                messagelength = messagelength + nvars*facearea
            END IF

            IF (messagelength > 0) THEN
                IF (i == SIZE(recvconns, 2)) THEN
                    CALL add_mpi_task(mpirtasks, impirtask, &
                        iprocnbr, messagelength, recvcounter, "recv")
                ELSE IF (recvconns(2, i + 1) /= iprocnbr) THEN
                    CALL add_mpi_task(mpirtasks, impirtask, &
                        iprocnbr, messagelength, recvcounter, "recv")
                END IF
            END IF
        END DO

        ! Set the output task counters
        nmpirtasks = impirtask

        ! Add a harmful dummy task at (ntasks+1) for checking
        mpirtasks(:, nmpirtasks+1) = -1

    END SUBROUTINE prepare_mpirecvtasks


    ! Perform all non-blocking MPI Recv calls
    !
    SUBROUTINE recv_mpi_all(minconlvl, maxconlvl, nplane, vertices,  &
            normal, fwd, flag, nvars)

        ! Subroutine arguments
        INTEGER(intk), INTENT(in) :: minconlvl, maxconlvl, nplane
        LOGICAL, INTENT(in) :: vertices, normal
        INTEGER(intk), INTENT(in) :: fwd
        CHARACTER(len=1), INTENT(in) :: flag
        INTEGER(intk), INTENT(in) :: nvars

        ! Local variables
        INTEGER(intk) :: i, iprocnbr, igrid, iface, facearea
        LOGICAL :: exchange, geometry
        INTEGER(int32) :: recvcounter, messagelength

        geometry = .FALSE.

        ! Post all receive calls
        recvcounter = 0
        messagelength = 0
        nrecv = 0
        recvidxlist = 0

        DO i = 1, SIZE(recvconns, 2)
            exchange = decide(i, recvconns, geometry, vertices, fwd, &
                minconlvl, maxconlvl)
            iprocnbr = recvconns(2, i)

            ! Communication with self is handled specially in
            ! send_all - nothing to do here
            IF (iprocnbr == myid) THEN
                CYCLE
            END IF

            IF (exchange) THEN
                igrid = recvconns(3, i)
                iface = recvconns(5, i)

                facearea = face_area(igrid, iface, nplane, flag)
                recvidxlist(1, i) = iprocnbr
                recvidxlist(2, i) = nvars*facearea
                recvidxlist(3, i) = recvcounter + messagelength
                messagelength = messagelength + nvars*facearea
            END IF

            IF (messagelength > 0) THEN
                IF (i == SIZE(recvconns, 2)) THEN
                    CALL post_recv(iprocnbr, messagelength, recvcounter)
                ELSE IF (recvconns(2, i + 1) /= iprocnbr) THEN
                    CALL post_recv(iprocnbr, messagelength, recvcounter)
                END IF
            END IF
        END DO
    END SUBROUTINE recv_mpi_all

END MODULE conn2_mod
