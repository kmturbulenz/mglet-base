MODULE ftoc2_mod
    USE MPI_f08
    USE core_mod
    USE ftoc1_mod, ONLY: ftoc1
    USE ftoc_core_mod
    USE ibcore_mod, ONLY: ib
    USE gc_restrict_mod
    USE noib_restrict_mod
    USE restrict_mod

    ! FTOC2 uses the sendbuf from commbuf_mod simultaneously while packing for
    ! MPI sends and for self communication. Since the self communication must
    ! never overlap with packing, and to avoid allocating a temporary buffer
    ! for self communication, packing and self communication use different
    ! parts of the same sendbuf. Due to the nature of fine to coarse, where
    ! the restricted grid (1/8 of the fine grid) is sent to the coarser grid,
    ! the size should not be exceeded at any point.
    ! The sendbuf is used as following:
    ! Packing: start to max size for worst-case ftoc2 call
    ! Self communication: max size for worst-case ftoc2 call to end

    IMPLICIT NONE (type, external)
    PRIVATE

    ! Lists that hold the send and receive request arrays
    TYPE(MPI_Request), ALLOCATABLE :: sendreqs(:), recvreqs(:)
    INTEGER(int32), ALLOCATABLE :: recvlist(:)
    INTEGER(intk), ALLOCATABLE :: recvidxlist(:, :)
    INTEGER(intk) :: nsend, nrecv

    ! Task extents
    INTEGER(intk), PARAMETER :: sendtasksize = 11
    INTEGER(intk), PARAMETER :: recvtasksize = 15
    INTEGER(intk), PARAMETER :: selftasksize = 15
    INTEGER(intk), PARAMETER :: mpitasksize = 3
    INTEGER(intk) :: maxsendtasks, maxrecvtasks
    INTEGER(intk) :: maxselftasks
    INTEGER(intk) :: maxmpisendtasks, maxmpirecvtasks

    ! Pointers to fields passed to ftoc2, point to dummy field if not present
    TYPE(field_t), POINTER :: f1 => NULL(), f2 => NULL(), f3 => NULL(), &
        f4 => NULL(), f5 => NULL(), f6 => NULL()

    ! Workpackages containing individual tasks for packing / unpacking
    INTEGER(intk), ALLOCATABLE :: sendtasks(:, :), recvtasks(:, :)
    INTEGER(intk), ALLOCATABLE :: selftasks(:, :)
    INTEGER(intk), ALLOCATABLE :: mpisendtasks(:, :), mpirecvtasks(:, :)
    !$omp declare target(sendtasks, recvtasks, selftasks)

    TYPE :: work_t
        LOGICAL :: is_init = .FALSE.
        INTEGER(intk), ALLOCATABLE :: sendtasks(:, :)
        INTEGER(intk), ALLOCATABLE :: recvtasks(:, :)
        INTEGER(intk), ALLOCATABLE :: selftasks(:, :)
        INTEGER(intk), ALLOCATABLE :: mpisendtasks(:, :)
        INTEGER(intk), ALLOCATABLE :: mpirecvtasks(:, :)
    END TYPE work_t

    ! Array to store instructions for different argument combinations
    TYPE(work_t), ALLOCATABLE, TARGET :: workrecords(:, :, :)
    LOGICAL :: is_recording = .FALSE.

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
        TYPE(field_t), TARGET, INTENT(in) :: ff
        TYPE(field_t), TARGET, INTENT(inout) :: fc
        CHARACTER(len=1), INTENT(in) :: flag

        ! Local variables
        ! none...

        IF (.NOT. is_init) CALL errr(__FILE__, __LINE__)

        IF (TRIM(ff%name) /= TRIM(fc%name)) THEN
            WRITE(*, *) "ftoc2_one only supports equal ff and fc fields."
            WRITE(*, *) "Fine to coarse from a fine field to a different coarse"
            WRITE(*, *) "field is not supported on the device. Use a temporary"
            WRITE(*, *) "field instead and copy result over."
            WRITE(*, *) "ff: ", TRIM(ff%name)
            WRITE(*, *) "fc: ", TRIM(fc%name)
            CALL errr(__FILE__, __LINE__)
        END IF

        ! ff and fc are the same field
        CALL ftoc2_impl(ilevel, v1=fc, flag=flag)

        nrecv = 0
        nsend = 0
    END SUBROUTINE ftoc2_one


    SUBROUTINE ftoc2_multiple(ilevel, v1, v2, v3, s1, s2, s3)
        ! Subroutine arguments
        INTEGER(intk), INTENT(in) :: ilevel
        TYPE(field_t), TARGET, INTENT(inout), OPTIONAL :: v1, v2, v3, s1, s2, s3

        ! Local variables
        ! none...

        IF (.NOT. is_init) CALL errr(__FILE__, __LINE__)

        CALL ftoc2_impl(ilevel, v1, v2, v3, s1, s2, s3, flag='*')

        nrecv = 0
        nsend = 0
    END SUBROUTINE ftoc2_multiple


    SUBROUTINE ftoc2_impl(ilevel, v1, v2, v3, s1, s2, s3, flag)
        ! Subroutine arguments
        INTEGER(intk), INTENT(in) :: ilevel
        TYPE(field_t), TARGET, OPTIONAL, INTENT(inout) :: v1, v2, v3, s1, s2, s3
        CHARACTER(len=1), INTENT(in) :: flag

        ! Local variables
        TYPE(field_t), POINTER :: dummy
        INTEGER(intk) :: nvars, idx_args, idx_flag

        CALL push_field(dummy, "FTOC2_DUMMY")

        f1 => dummy
        f2 => dummy
        f3 => dummy
        f4 => dummy
        f5 => dummy
        f6 => dummy

        nvars = 0
        idx_args = 0
        IF (PRESENT(v1)) THEN
            nvars = nvars + 1
            f1 => v1
            idx_args = idx_args + 1 * 2**(0)
        END IF
        IF (PRESENT(v2)) THEN
            nvars = nvars + 1
            f2 => v2
            idx_args = idx_args + 1 * 2**(1)
        END IF
        IF (PRESENT(v3)) THEN
            nvars = nvars + 1
            f3 => v3
            idx_args = idx_args + 1 * 2**(2)
        END IF
        IF (PRESENT(s1)) THEN
            nvars = nvars + 1
            f4 => s1
            idx_args = idx_args + 1 * 2**(3)
        END IF
        IF (PRESENT(s2)) THEN
            nvars = nvars + 1
            f5 => s2
            idx_args = idx_args + 1 * 2**(4)
        END IF
        IF (PRESENT(s3)) THEN
            nvars = nvars + 1
            f6 => s3
            idx_args = idx_args + 1 * 2**(5)
        END IF

        IF (idx_args == 0) THEN
            WRITE(*, *) "No fields present."
            CALL errr(__FILE__, __LINE__)
        END IF

        idx_flag = IACHAR(flag)

        ASSOCIATE(wptr => workrecords(ilevel, idx_args, idx_flag))

        IF (is_recording) THEN
            ! During the recording phase, the workpackage is not yet initialized
            ! and needs to be created. Tasks arrays are offloaded for later
            ! efficient task execution on device
            CALL recording_pass(wptr, ilevel, flag, nvars, &
                v1, v2, v3, s1, s2, s3)
        ELSE
            IF (wptr%is_init) THEN
                CALL recorded_ftoc2(wptr)
            ELSE
                CALL jit_ftoc2(ilevel, flag, v1, v2, v3, s1, s2, s3)
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
    END SUBROUTINE ftoc2_impl


    SUBROUTINE jit_ftoc2(ilevel, flag, v1, v2, v3, s1, s2, s3)
        ! Subroutine arguments
        INTEGER(intk), INTENT(in) :: ilevel
        CHARACTER(len=1), INTENT(in) :: flag
        TYPE(field_t), OPTIONAL, INTENT(inout) :: v1, v2, v3, s1, s2, s3

        ! Local variables
        INTEGER(intk) :: nsendtasks, nrecvtasks
        INTEGER(intk) :: nselftasks
        INTEGER(intk) :: nmpisendtasks

        ! Manufacturing the workpackage just in-time and execute
        CALL prepare_tasks_all(sendtasks, nsendtasks, selftasks, nselftasks, &
            mpisendtasks, nmpisendtasks, ilevel, flag, v1, v2, v3, s1, s2, s3)

        !$omp target update to( &
        !$omp& sendtasks(1:sendtasksize, 1:nsendtasks+1), &
        !$omp& selftasks(1:selftasksize, 1:nselftasks+1)) nowait

        CALL recv_mpi_all(ilevel, flag, v1, v2, v3, s1, s2, s3)

        !$omp taskwait

        CALL process_sendtasks(nsendtasks, sendtasks)
        CALL process_mpisend(nmpisendtasks, mpisendtasks)
        CALL process_selftasks(nselftasks, selftasks)

        CALL prepare_recvtasks_all(recvtasks, nrecvtasks, flag, v1, v2, v3, &
            s1, s2, s3)

        !$omp target update to(recvtasks(1:recvtasksize, 1:nrecvtasks+1))

        CALL process_recvtasks(nrecvtasks, recvtasks)
    END SUBROUTINE jit_ftoc2


    SUBROUTINE recorded_ftoc2(wptr)
        ! Subroutine arguments
        TYPE(work_t), INTENT(inout) :: wptr

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
    END SUBROUTINE recorded_ftoc2


    SUBROUTINE recv_mpi_all(ilevel, flag, v1, v2, v3, s1, s2, s3)
        ! Subroutine arguments
        INTEGER(intk), INTENT(in) :: ilevel
        CHARACTER(len=1), INTENT(in) :: flag
        TYPE(field_t), INTENT(in), OPTIONAL :: v1, v2, v3, s1, s2, s3

        ! Local variables
        INTEGER(intk) :: i, iprocnbr, igridf
        INTEGER(int32) :: recvcounter, messagelength, ncells

        recvcounter = 0
        messagelength = 0
        nrecv = 0
        recvidxlist = -HUGE(1_intk)
        recvlist = 0

        DO i = 1, irecv
            igridf = recvconns(3, i)
            iprocnbr = recvconns(1, i)  ! The sender process (fine side)
            IF (iprocnbr == myid) CYCLE

            IF (ilevel == level(igridf)) THEN
                CALL count_ncells(ncells, flag, igridf, v1, v2, v3, s1, s2, s3)
                recvidxlist(1, i) = iprocnbr
                recvidxlist(2, i) = ncells
                recvidxlist(3, i) = recvcounter + messagelength
                messagelength = messagelength + ncells

                IF (recvcounter + messagelength > idim_mg_bufs) THEN
                    CALL errr(__FILE__, __LINE__)
                END IF
            END IF

            IF (messagelength > 0) THEN
                IF (i == irecv) THEN
                    CALL post_recv(iprocnbr, messagelength, recvcounter)
                ELSE IF (recvconns(1, i + 1) /= iprocnbr) THEN
                    CALL post_recv(iprocnbr, messagelength, recvcounter)
                END IF
            END IF
        END DO
    END SUBROUTINE recv_mpi_all


    SUBROUTINE post_recv(iprocnbr, messagelength, recvcounter)
        ! Subroutine arguments
        INTEGER(intk), INTENT(in) :: iprocnbr
        INTEGER(int32), INTENT(inout) :: messagelength
        INTEGER(int32), INTENT(inout) :: recvcounter

        ! Local variables
        ! none...

        nrecv = nrecv + 1
        recvlist(nrecv) = iprocnbr

        !$omp target data use_device_addr(recvbuf)
        CALL MPI_Irecv(recvbuf(recvcounter+1), messagelength, &
            mglet_mpi_real, iprocnbr, 1, MPI_COMM_WORLD, recvreqs(nrecv))
        !$omp end target data

        recvcounter = recvcounter + messagelength
        messagelength = 0
    END SUBROUTINE post_recv


    SUBROUTINE recording_pass(wptr, ilevel, flag, nvars, v1, v2, v3, s1, s2, s3)
        ! Subroutine arguments
        TYPE(work_t), INTENT(inout) :: wptr
        INTEGER(intk), INTENT(in) :: ilevel
        CHARACTER(len=1), INTENT(in) :: flag
        INTEGER(intk), INTENT(in) :: nvars
        TYPE(field_t), OPTIONAL, INTENT(inout) :: v1, v2, v3, s1, s2, s3

        ! Local variables
        INTEGER(intk) :: nsendtasks, nrecvtasks, nselftasks
        INTEGER(intk) :: nmpisendtasks, nmpirecvtasks

        IF (wptr%is_init) THEN
            WRITE(*, *) "Combination already recorded."
            CALL errr(__FILE__, __LINE__)
        END IF

        ! It is necessary to execute one cycle with communication
        ! as otherwise many valuable checks are not possible

        CALL prepare_mpirecvtasks(mpirecvtasks, nmpirecvtasks, ilevel, flag, &
            v1, v2, v3, s1, s2, s3)
        CALL prepare_tasks_all(sendtasks, nsendtasks, selftasks, nselftasks, &
            mpisendtasks, nmpisendtasks, ilevel, flag, v1, v2, v3, s1, s2, s3)

        !$omp target update to( &
        !$omp& sendtasks(1:sendtasksize, 1:nsendtasks+1), &
        !$omp& selftasks(1:selftasksize, 1:nselftasks+1))

        CALL process_mpirecv(nmpirecvtasks, mpirecvtasks)
        CALL process_sendtasks(nsendtasks, sendtasks)
        CALL process_mpisend(nmpisendtasks, mpisendtasks)
        CALL process_selftasks(nselftasks, selftasks)
        CALL prepare_recvtasks_all(recvtasks, nrecvtasks, flag, &
            v1, v2, v3, s1, s2, s3)

        !$omp target update to(recvtasks(1:recvtasksize, 1:nrecvtasks+1))
        CALL process_recvtasks(nrecvtasks, recvtasks)

        ! Allocate the workpackage arrays in the exact sizes
        ALLOCATE(wptr%sendtasks(sendtasksize, nsendtasks+1))
        ALLOCATE(wptr%recvtasks(recvtasksize, nrecvtasks+1))
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
        !$omp&  wptr%sendtasks(1:sendtasksize, 1:nsendtasks+1), &
        !$omp&  wptr%recvtasks(1:recvtasksize, 1:nrecvtasks+1), &
        !$omp&  wptr%selftasks(1:selftasksize, 1:nselftasks+1))

        ! Mark the workpackage as initialized
        wptr%is_init = .TRUE.
    END SUBROUTINE recording_pass


    SUBROUTINE process_recvtasks(nrtasks, rtasks)
        ! Subroutine arguments
        INTEGER(intk), INTENT(in) :: nrtasks
        INTEGER(intk), INTENT(in) :: rtasks(recvtasksize, nrtasks+1)

        ! Local variables
        ! none...

        IF (nrtasks == 0) RETURN

#ifdef _MGLET_PROFILE_ANNOTATIONS_
        CALL profile_range_push("process_recvtasks")
#endif

        CALL process_recvtasks_impl(nrtasks, rtasks, f1%arr, f2%arr, &
            f3%arr, f4%arr, f5%arr, f6%arr)

#ifdef _MGLET_PROFILE_ANNOTATIONS_
        CALL profile_range_pop()
#endif

        CALL check_recvtasks_dummy(nrtasks, rtasks)
        END SUBROUTINE process_recvtasks


        SUBROUTINE process_recvtasks_impl(nrtasks, rtasks, a1, a2, a3, a4, &
            a5, a6)
        ! Subroutine arguments
        INTEGER(intk), INTENT(in) :: nrtasks
        INTEGER(intk), INTENT(in) :: rtasks(recvtasksize, nrtasks+1)
        REAL(realk), INTENT(inout) :: a1(*), a2(*), a3(*), a4(*), a5(*), a6(*)

        ! Local variables
        INTEGER(intk) :: itask, fieldid, tasksize, recvidx, igridc
        CHARACTER(len=1) :: flag
        INTEGER(intk) :: kk, jj, ii, ip3
        INTEGER(intk) :: istart, istop, jstart, jstop, kstart, kstop
        INTEGER(intk) :: ipos, jpos, kpos

        !$omp target teams distribute private(itask, fieldid, flag, &
        !$omp& tasksize, recvidx, igridc, istart, istop, jstart, jstop, &
        !$omp& kstart, kstop, ipos, jpos, kpos, kk, jj, ii, ip3)
        DO itask = 1, nrtasks
            fieldid = rtasks(1, itask)
            flag = ACHAR(rtasks(2, itask))
            tasksize = rtasks(3, itask)
            recvidx = rtasks(4, itask)
            igridc = rtasks(6, itask)
            istart = rtasks(7, itask)
            istop = rtasks(8, itask)
            jstart = rtasks(9, itask)
            jstop = rtasks(10, itask)
            kstart = rtasks(11, itask)
            kstop = rtasks(12, itask)
            ipos = rtasks(13, itask)
            jpos = rtasks(14, itask)
            kpos = rtasks(15, itask)

            CALL get_mgdims(kk, jj, ii, igridc)
            CALL get_ip3(ip3, igridc)

            !$omp parallel
            SELECT CASE(fieldid)
            CASE (1)
                CALL unpack_restricted_buffer(flag, kk, jj, ii, a1(ip3), &
                    recvbuf(recvidx:recvidx+tasksize-1), tasksize, istart, &
                    istop, jstart, jstop, kstart, kstop, ipos, jpos, kpos)
            CASE (2)
                CALL unpack_restricted_buffer(flag, kk, jj, ii, a2(ip3), &
                    recvbuf(recvidx:recvidx+tasksize-1), tasksize, istart, &
                    istop, jstart, jstop, kstart, kstop, ipos, jpos, kpos)
            CASE (3)
                CALL unpack_restricted_buffer(flag, kk, jj, ii, a3(ip3), &
                    recvbuf(recvidx:recvidx+tasksize-1), tasksize, istart, &
                    istop, jstart, jstop, kstart, kstop, ipos, jpos, kpos)
            CASE (4)
                CALL unpack_restricted_buffer(flag, kk, jj, ii, a4(ip3), &
                    recvbuf(recvidx:recvidx+tasksize-1), tasksize, istart, &
                    istop, jstart, jstop, kstart, kstop, ipos, jpos, kpos)
            CASE (5)
                CALL unpack_restricted_buffer(flag, kk, jj, ii, a5(ip3), &
                    recvbuf(recvidx:recvidx+tasksize-1), tasksize, istart, &
                    istop, jstart, jstop, kstart, kstop, ipos, jpos, kpos)
            CASE (6)
                CALL unpack_restricted_buffer(flag, kk, jj, ii, a6(ip3), &
                    recvbuf(recvidx:recvidx+tasksize-1), tasksize, istart, &
                    istop, jstart, jstop, kstart, kstop, ipos, jpos, kpos)
#ifdef _MGLET_DEBUG_
            CASE DEFAULT
                CALL errr(__FILE__, __LINE__)
#endif
            END SELECT
            !$omp end parallel
        END DO
        !$omp end target teams distribute
    END SUBROUTINE process_recvtasks_impl


    SUBROUTINE check_recvtasks_dummy(nrtasks, rtasks)
        ! Subroutine arguments
        INTEGER(intk), INTENT(in) :: nrtasks
        INTEGER(intk), INTENT(in) :: rtasks(recvtasksize, nrtasks+1)

        ! Local variables
        ! none...

        IF (.NOT. ALL(rtasks(:, nrtasks+1) == -1)) THEN
            WRITE(*, *) "Did not encounter the expected dummy task."
            CALL errr(__FILE__, __LINE__)
        END IF
    END SUBROUTINE check_recvtasks_dummy


    SUBROUTINE prepare_recvtasks_all(rtasks, nrtasks, flag, &
            v1, v2, v3, s1, s2, s3)
        ! Subroutine arguments
        INTEGER(intk), INTENT(inout) :: rtasks(recvtasksize, maxrecvtasks)
        INTEGER(intk), INTENT(out) :: nrtasks
        CHARACTER(len=1), INTENT(in) :: flag
        TYPE(field_t), OPTIONAL, INTENT(inout) :: v1, v2, v3, s1, s2, s3

        ! Local variables
        INTEGER(intk) :: irtask
        TYPE(MPI_Status) :: recvstatus
        INTEGER(int32) :: idx, i, recvmessagelen, unpacklen

#ifdef _MGLET_PROFILE_ANNOTATIONS_
        CALL profile_range_push("prepare_recvtasks_all")
#endif

        irtask = 0
        DO WHILE(.TRUE.)
            IF (nrecv == 0) EXIT
            CALL MPI_Waitany(nrecv, recvreqs, idx, recvstatus)

            IF (idx /= MPI_UNDEFINED) THEN
                CALL MPI_Get_count(recvstatus, mglet_mpi_real, recvmessagelen)

                unpacklen = 0
                DO i = 1, irecv
                    IF (recvidxlist(1, i) == recvlist(idx) &
                            .AND. recvidxlist(2, i) > 0) THEN
                        CALL prepare_recvtask(rtasks, irtask, i, flag, &
                            v1, v2, v3, s1, s2, s3)
                        unpacklen = unpacklen + recvidxlist(2, i)
                    END IF
                END DO

                IF (recvmessagelen /= unpacklen) THEN
                    CALL errr(__FILE__, __LINE__)
                END IF
            ELSE
                EXIT
            END IF
        END DO

        CALL MPI_Waitall(nsend, sendreqs, MPI_STATUSES_IGNORE)
        nrtasks = irtask
        rtasks(:, nrtasks+1) = -1

#ifdef _MGLET_PROFILE_ANNOTATIONS_
        CALL profile_range_pop()
#endif
    END SUBROUTINE prepare_recvtasks_all


    SUBROUTINE prepare_recvtask(rtasks, irtask, recvid, flag, &
            v1, v2, v3, s1, s2, s3)
        ! Subroutine arguments
        INTEGER(intk), INTENT(inout) :: rtasks(recvtasksize, maxrecvtasks)
        INTEGER(intk), INTENT(inout) :: irtask
        INTEGER(intk), INTENT(in) :: recvid
        CHARACTER(len=1), INTENT(in) :: flag
        TYPE(field_t), OPTIONAL, INTENT(inout) :: v1, v2, v3, s1, s2, s3

        ! Local variables
        INTEGER(intk) :: igridf, igridc
        INTEGER(int32) :: ncells, offset, icount
        CHARACTER(len=1) :: flag_u

        igridf = recvconns(3, recvid)
        igridc = recvconns(4, recvid)

        CALL count_ncells(ncells, flag, igridf, &
            v1, v2, v3, s1, s2, s3)

        offset = INT(recvidxlist(3, recvid), kind=int32) + 1
        icount = offset

        IF (PRESENT(v1)) THEN
            IF (flag == '*') THEN
                flag_u = 'U'
            ELSE
                flag_u = flag
            END IF

            irtask = irtask + 1
            CALL add_single_recvtask(rtasks(:, irtask), 1, flag_u, icount, &
                igridf, igridc)
        END IF

        IF (PRESENT(v2)) THEN
            irtask = irtask + 1
            CALL add_single_recvtask(rtasks(:, irtask), 2, 'V', icount, &
                igridf, igridc)
        END IF

        IF (PRESENT(v3)) THEN
            irtask = irtask + 1
            CALL add_single_recvtask(rtasks(:, irtask), 3, 'W', icount, &
                igridf, igridc)
        END IF

        IF (PRESENT(s1)) THEN
            irtask = irtask + 1
            CALL add_single_recvtask(rtasks(:, irtask), 4, 'P', icount, &
                igridf, igridc)
        END IF

        IF (PRESENT(s2)) THEN
            irtask = irtask + 1
            CALL add_single_recvtask(rtasks(:, irtask), 5, 'P', icount, &
                igridf, igridc)
        END IF

        IF (PRESENT(s3)) THEN
            irtask = irtask + 1
            CALL add_single_recvtask(rtasks(:, irtask), 6, 'P', icount, &
                igridf, igridc)
        END IF

        IF (ncells /= (icount - offset)) THEN
            WRITE(*, *) "ncells:", ncells, "recv_len:", icount - offset
            CALL errr(__FILE__, __LINE__)
        END IF
    END SUBROUTINE prepare_recvtask


    SUBROUTINE process_selftasks(netasks, etasks)
        ! Subroutine arguments
        INTEGER(intk), INTENT(in) :: netasks
        INTEGER(intk), INTENT(in) :: etasks(selftasksize, netasks+1)

        ! Local variables
        TYPE(field_t), POINTER :: ddx, ddy, ddz, bp, bt
        LOGICAL :: foundbt

        IF (netasks == 0) RETURN

        CALL get_field(ddx, "DDX")
        CALL get_field(ddy, "DDY")
        CALL get_field(ddz, "DDZ")

        SELECT CASE (ib%type)
        CASE ("NONE")
            CALL process_selftasks_noib(netasks, etasks, ddx, ddy, ddz)
        CASE ("GHOSTCELL")
            CALL get_field(bp, "BP")
            CALL get_field(bt, "BT", foundbt)
            IF (.NOT. foundbt) THEN
                ! For non-scalar cases, point bt to bp. It is unused anyways.
                bt => bp
            END IF
            CALL process_selftasks_gc(netasks, etasks, ddx, ddy, ddz, bp, bt)
        CASE DEFAULT
            CALL errr(__FILE__, __LINE__)
        END SELECT

        ! Safety check based on final dummy entry
        IF (.NOT. ALL(etasks(:, netasks+1) == -1)) THEN
            WRITE(*, *) "Did not encounter the expected dummy task."
            CALL errr(__FILE__, __LINE__)
        END IF
    END SUBROUTINE process_selftasks


    SUBROUTINE process_selftasks_noib(netasks, etasks, ddx, ddy, ddz)
        ! Subroutine arguments
        INTEGER(intk), INTENT(in) :: netasks
        INTEGER(intk), INTENT(in) :: etasks(selftasksize, netasks+1)
        TYPE(field_t), INTENT(in) :: ddx, ddy, ddz

        ! Local variables
        ! none...

#ifdef _MGLET_PROFILE_ANNOTATIONS_
        CALL profile_range_push("process_selftasks_noib")
#endif

        CALL process_selftasks_noib_impl(netasks, etasks, f1%arr, f2%arr, &
            f3%arr, f4%arr, f5%arr, f6%arr, ddx%arr, ddy%arr, ddz%arr)

#ifdef _MGLET_PROFILE_ANNOTATIONS_
        CALL profile_range_pop()
#endif
    END SUBROUTINE process_selftasks_noib


    SUBROUTINE process_selftasks_noib_impl(netasks, etasks, a1, a2, a3, a4, &
            a5, a6, ddx, ddy, ddz)
        ! Subroutine arguments
        INTEGER(intk), INTENT(in) :: netasks
        INTEGER(intk), INTENT(in) :: etasks(selftasksize, netasks+1)
        REAL(realk), INTENT(inout) :: a1(*), a2(*), a3(*), a4(*), a5(*), a6(*)
        REAL(realk), INTENT(in) :: ddx(*), ddy(*), ddz(*)

        ! Local variables
        INTEGER(intk) :: itask, kkf, jjf, iif, kkc, jjc, iic
        INTEGER(intk) :: ip3f, ip3c, ipx, ipy, ipz
        INTEGER(intk) :: fieldid, igridf, igridc, istart, istop, jstart, &
            jstop, kstart, kstop, ipos, jpos, kpos
        INTEGER(int32) :: tasksize, scratchidx
        CHARACTER(len=1) :: flag

        !$omp target teams distribute private(itask, fieldid, flag, &
        !$omp& tasksize, scratchidx, igridf, igridc, istart, istop, jstart, &
        !$omp& jstop, kstart, kstop, ipos, jpos, kpos, kkf, jjf, iif, kkc, &
        !$omp& jjc, iic, ip3f, ip3c, ipx, ipy, ipz)
        DO itask = 1, netasks
            fieldid = etasks(1, itask)
            flag = ACHAR(etasks(2, itask))
            tasksize = INT(etasks(3, itask), kind=int32)
            scratchidx = INT(etasks(4, itask), kind=int32)
            igridf = etasks(5, itask)
            igridc = etasks(6, itask)
            istart = etasks(7, itask)
            istop = etasks(8, itask)
            jstart = etasks(9, itask)
            jstop = etasks(10, itask)
            kstart = etasks(11, itask)
            kstop = etasks(12, itask)
            ipos = etasks(13, itask)
            jpos = etasks(14, itask)
            kpos = etasks(15, itask)

            CALL get_mgdims(kkf, jjf, iif, igridf)
            CALL get_mgdims(kkc, jjc, iic, igridc)
            CALL get_ip3(ip3f, igridf)
            CALL get_ip3(ip3c, igridc)
            CALL get_ip1x(ipx, igridf)
            CALL get_ip1y(ipy, igridf)
            CALL get_ip1z(ipz, igridf)

            !$omp parallel
            SELECT CASE(fieldid)
            CASE (1)
                CALL restrict_noib_flag(flag, kkf, jjf, iif, &
                    a1(ip3f), ddx(ipx), ddy(ipy), ddz(ipz), &
                    scratchidx, tasksize, istart, istop, jstart, jstop, &
                    kstart, kstop)
                CALL unpack_restricted_buffer(flag, kkc, jjc, iic, &
                    a1(ip3c), sendbuf(scratchidx:scratchidx+tasksize-1), &
                    tasksize, istart, istop, jstart, jstop, kstart, kstop, &
                    ipos, jpos, kpos)
            CASE (2)
                CALL restrict_noib_flag(flag, kkf, jjf, iif, &
                    a2(ip3f), ddx(ipx), ddy(ipy), ddz(ipz), &
                    scratchidx, tasksize, istart, istop, jstart, jstop, &
                    kstart, kstop)
                CALL unpack_restricted_buffer(flag, kkc, jjc, iic, &
                    a2(ip3c), sendbuf(scratchidx:scratchidx+tasksize-1), &
                    tasksize, istart, istop, jstart, jstop, kstart, kstop, &
                    ipos, jpos, kpos)
            CASE (3)
                CALL restrict_noib_flag(flag, kkf, jjf, iif, &
                    a3(ip3f), ddx(ipx), ddy(ipy), ddz(ipz), &
                    scratchidx, tasksize, istart, istop, jstart, jstop, &
                    kstart, kstop)
                CALL unpack_restricted_buffer(flag, kkc, jjc, iic, &
                    a3(ip3c), sendbuf(scratchidx:scratchidx+tasksize-1), &
                    tasksize, istart, istop, jstart, jstop, kstart, kstop, &
                    ipos, jpos, kpos)
            CASE (4)
                CALL restrict_noib_flag(flag, kkf, jjf, iif, &
                    a4(ip3f), ddx(ipx), ddy(ipy), ddz(ipz), &
                    scratchidx, tasksize, istart, istop, jstart, jstop, &
                    kstart, kstop)
                CALL unpack_restricted_buffer(flag, kkc, jjc, iic, &
                    a4(ip3c), sendbuf(scratchidx:scratchidx+tasksize-1), &
                    tasksize, istart, istop, jstart, jstop, kstart, kstop, &
                    ipos, jpos, kpos)
            CASE (5)
                CALL restrict_noib_flag(flag, kkf, jjf, iif, &
                    a5(ip3f), ddx(ipx), ddy(ipy), ddz(ipz), &
                    scratchidx, tasksize, istart, istop, jstart, jstop, &
                    kstart, kstop)
                CALL unpack_restricted_buffer(flag, kkc, jjc, iic, &
                    a5(ip3c), sendbuf(scratchidx:scratchidx+tasksize-1), &
                    tasksize, istart, istop, jstart, jstop, kstart, kstop, &
                    ipos, jpos, kpos)
            CASE (6)
                CALL restrict_noib_flag(flag, kkf, jjf, iif, &
                    a6(ip3f), ddx(ipx), ddy(ipy), ddz(ipz), &
                    scratchidx, tasksize, istart, istop, jstart, jstop, &
                    kstart, kstop)
                CALL unpack_restricted_buffer(flag, kkc, jjc, iic, &
                    a6(ip3c), sendbuf(scratchidx:scratchidx+tasksize-1), &
                    tasksize, istart, istop, jstart, jstop, kstart, kstop, &
                    ipos, jpos, kpos)
#ifdef _MGLET_DEBUG_
            CASE DEFAULT
                CALL errr(__FILE__, __LINE__)
#endif
            END SELECT
            !$omp end parallel
        END DO
        !$omp end target teams distribute
    END SUBROUTINE process_selftasks_noib_impl


    SUBROUTINE process_selftasks_gc(netasks, etasks, ddx_f, ddy_f, ddz_f, &
            bp_f, bt_f)
        ! Subroutine arguments
        INTEGER(intk), INTENT(in) :: netasks
        INTEGER(intk), INTENT(in) :: etasks(selftasksize, netasks+1)
        TYPE(field_t), INTENT(in) :: ddx_f, ddy_f, ddz_f, bp_f, bt_f

        ! Local variables
        ! none...

#ifdef _MGLET_PROFILE_ANNOTATIONS_
        CALL profile_range_push("process_selftasks_gc")
#endif

        CALL process_selftasks_gc_impl(netasks, etasks, f1%arr, f2%arr, &
            f3%arr, f4%arr, f5%arr, f6%arr, ddx_f%arr, ddy_f%arr, &
            ddz_f%arr, bp_f%arr, bt_f%arr)

#ifdef _MGLET_PROFILE_ANNOTATIONS_
        CALL profile_range_pop()
#endif
    END SUBROUTINE process_selftasks_gc


    SUBROUTINE process_selftasks_gc_impl(netasks, etasks, a1, a2, a3, a4, &
            a5, a6, ddx, ddy, ddz, bp, bt)
        ! Subroutine arguments
        INTEGER(intk), INTENT(in) :: netasks
        INTEGER(intk), INTENT(in) :: etasks(selftasksize, netasks+1)
        REAL(realk), INTENT(inout) :: a1(*), a2(*), a3(*), a4(*), a5(*), a6(*)
        REAL(realk), INTENT(in) :: ddx(*), ddy(*), ddz(*), bp(*), bt(*)

        ! Local variables
        INTEGER(intk) :: itask, kkf, jjf, iif, kkc, jjc, iic
        INTEGER(intk) :: ip3f, ip3c, ipx, ipy, ipz
        INTEGER(intk) :: fieldid, igridf, igridc, istart, istop, jstart, &
            jstop, kstart, kstop, ipos, jpos, kpos
        INTEGER(int32) :: tasksize, scratchidx
        CHARACTER(len=1) :: flag

        !$omp target teams distribute private(itask, fieldid, flag, &
        !$omp& tasksize, scratchidx, igridf, igridc, istart, istop, jstart, &
        !$omp& jstop, kstart, kstop, ipos, jpos, kpos, kkf, jjf, iif, kkc, &
        !$omp& jjc, iic, ip3f, ip3c, ipx, ipy, ipz)
        DO itask = 1, netasks
            fieldid = etasks(1, itask)
            flag = ACHAR(etasks(2, itask))
            tasksize = INT(etasks(3, itask), kind=int32)
            scratchidx = INT(etasks(4, itask), kind=int32)
            igridf = etasks(5, itask)
            igridc = etasks(6, itask)
            istart = etasks(7, itask)
            istop = etasks(8, itask)
            jstart = etasks(9, itask)
            jstop = etasks(10, itask)
            kstart = etasks(11, itask)
            kstop = etasks(12, itask)
            ipos = etasks(13, itask)
            jpos = etasks(14, itask)
            kpos = etasks(15, itask)

            CALL get_mgdims(kkf, jjf, iif, igridf)
            CALL get_mgdims(kkc, jjc, iic, igridc)
            CALL get_ip3(ip3f, igridf)
            CALL get_ip3(ip3c, igridc)
            CALL get_ip1x(ipx, igridf)
            CALL get_ip1y(ipy, igridf)
            CALL get_ip1z(ipz, igridf)

            !$omp parallel
            SELECT CASE(fieldid)
            CASE (1)
                CALL restrict_gc_flag(flag, kkf, jjf, iif, a1(ip3f), &
                    ddx(ipx), ddy(ipy), ddz(ipz), bp(ip3f), bt(ip3f), &
                    scratchidx, tasksize, istart, istop, jstart, jstop, &
                    kstart, kstop)
                !$omp barrier
                CALL unpack_restricted_buffer(flag, kkc, jjc, iic, &
                    a1(ip3c), sendbuf(scratchidx:scratchidx+tasksize-1), &
                    tasksize, istart, istop, jstart, jstop, kstart, kstop, &
                    ipos, jpos, kpos)
                !$omp barrier
            CASE (2)
                CALL restrict_gc_flag(flag, kkf, jjf, iif, a2(ip3f), &
                    ddx(ipx), ddy(ipy), ddz(ipz), bp(ip3f), bt(ip3f), &
                    scratchidx, tasksize, istart, istop, jstart, jstop, &
                    kstart, kstop)
                !$omp barrier
                CALL unpack_restricted_buffer(flag, kkc, jjc, iic, &
                    a2(ip3c), sendbuf(scratchidx:scratchidx+tasksize-1), &
                    tasksize, istart, istop, jstart, jstop, kstart, kstop, &
                    ipos, jpos, kpos)
                !$omp barrier
            CASE (3)
                CALL restrict_gc_flag(flag, kkf, jjf, iif, a3(ip3f), &
                    ddx(ipx), ddy(ipy), ddz(ipz), bp(ip3f), bt(ip3f), &
                    scratchidx, tasksize, istart, istop, jstart, jstop, &
                    kstart, kstop)
                !$omp barrier
                CALL unpack_restricted_buffer(flag, kkc, jjc, iic, &
                    a3(ip3c), sendbuf(scratchidx:scratchidx+tasksize-1), &
                    tasksize, istart, istop, jstart, jstop, kstart, kstop, &
                    ipos, jpos, kpos)
                !$omp barrier
            CASE (4)
                CALL restrict_gc_flag(flag, kkf, jjf, iif, a4(ip3f), &
                    ddx(ipx), ddy(ipy), ddz(ipz), bp(ip3f), bt(ip3f), &
                    scratchidx, tasksize, istart, istop, jstart, jstop, &
                    kstart, kstop)
                !$omp barrier
                CALL unpack_restricted_buffer(flag, kkc, jjc, iic, &
                    a4(ip3c), sendbuf(scratchidx:scratchidx+tasksize-1), &
                    tasksize, istart, istop, jstart, jstop, kstart, kstop, &
                    ipos, jpos, kpos)
                !$omp barrier
            CASE (5)
                CALL restrict_gc_flag(flag, kkf, jjf, iif, a5(ip3f), &
                    ddx(ipx), ddy(ipy), ddz(ipz), bp(ip3f), bt(ip3f), &
                    scratchidx, tasksize, istart, istop, jstart, jstop, &
                    kstart, kstop)
                !$omp barrier
                CALL unpack_restricted_buffer(flag, kkc, jjc, iic, &
                    a5(ip3c), sendbuf(scratchidx:scratchidx+tasksize-1), &
                    tasksize, istart, istop, jstart, jstop, kstart, kstop, &
                    ipos, jpos, kpos)
                !$omp barrier
            CASE (6)
                CALL restrict_gc_flag(flag, kkf, jjf, iif, a6(ip3f), &
                    ddx(ipx), ddy(ipy), ddz(ipz), bp(ip3f), bt(ip3f), &
                    scratchidx, tasksize, istart, istop, jstart, jstop, &
                    kstart, kstop)
                !$omp barrier
                CALL unpack_restricted_buffer(flag, kkc, jjc, iic, &
                    a6(ip3c), sendbuf(scratchidx:scratchidx+tasksize-1), &
                    tasksize, istart, istop, jstart, jstop, kstart, kstop, &
                    ipos, jpos, kpos)
                !$omp barrier
            CASE DEFAULT
                CALL errr(__FILE__, __LINE__)
            END SELECT
            !$omp end parallel
        END DO
        !$omp end target teams distribute
    END SUBROUTINE process_selftasks_gc_impl


    SUBROUTINE unpack_restricted_buffer(flag, kk, jj, ii, fc, buffer, &
            tasksize, istart, istop, jstart, jstop, kstart, kstop, ipos, &
            jpos, kpos)
        !$omp declare target
        ! Subroutine arguments
        CHARACTER(len=1), INTENT(in) :: flag
        INTEGER(intk), INTENT(in) :: kk, jj, ii, tasksize
        REAL(realk), INTENT(inout) :: fc(kk, jj, ii)
        REAL(realk), INTENT(in) :: buffer(tasksize)
        INTEGER(intk), INTENT(in) :: istart, istop, jstart, jstop, kstart, kstop
        INTEGER(intk), INTENT(in) :: ipos, jpos, kpos

        IF (flag == 'N') THEN
            CALL unpack_restricted_buffer_n(kk, jj, ii, fc, buffer, tasksize, &
                istart, istop, jstart, jstop, kstart, kstop, ipos, jpos, kpos)
        ELSE
            CALL unpack_restricted_buffer_open(flag, kk, jj, ii, fc, &
                buffer, tasksize, istart, istop, jstart, jstop, kstart, &
                kstop, ipos, jpos, kpos)
        END IF
    END SUBROUTINE unpack_restricted_buffer


    SUBROUTINE unpack_restricted_buffer_open(flag, kk, jj, ii, fc, buffer, &
            tasksize, istart, istop, jstart, jstop, kstart, kstop, ipos, &
            jpos, kpos)
        !$omp declare target
        ! Subroutine arguments
        CHARACTER(len=1), INTENT(in) :: flag
        INTEGER(intk), INTENT(in) :: kk, jj, ii, tasksize
        REAL(realk), INTENT(inout) :: fc(kk, jj, ii)
        REAL(realk), INTENT(in) :: buffer(tasksize)
        INTEGER(intk), INTENT(in) :: istart, istop, jstart, jstop, kstart, kstop
        INTEGER(intk), INTENT(in) :: ipos, jpos, kpos

        ! Local variables
        INTEGER(intk) :: i, j, k, ic, jc, kc, ic0, jc0, kc0, icount
        INTEGER(intk) :: ni, nj, nk

        ic0 = 0
        jc0 = 0
        kc0 = 0
        SELECT CASE (flag)
        CASE ('U')
            ic0 = 1
        CASE ('V')
            jc0 = 1
        CASE ('W')
            kc0 = 1
        END SELECT

        ni = (istop - istart)/2 + 1
        nj = (jstop - jstart)/2 + 1
        nk = (kstop - kstart)/2 + 1

        !$omp do collapse(3) private(ic, jc, kc, icount)
        DO i = istart, istop, 2
            DO j = jstart, jstop, 2
                DO k = kstart, kstop, 2
                    ic = ipos + (i-3+ic0)/2 - ic0
                    jc = jpos + (j-3+jc0)/2 - jc0
                    kc = kpos + (k-3+kc0)/2 - kc0
                    icount = ((i-istart)/2*nj + (j-jstart)/2)*nk + &
                        (k-kstart)/2 + 1
                    fc(kc, jc, ic) = buffer(icount)
                END DO
            END DO
        END DO
        !$omp end do

#ifdef _MGLET_DEBUG_
        IF (ni*nj*nk /= tasksize) CALL errr(__FILE__, __LINE__)
#endif
    END SUBROUTINE unpack_restricted_buffer_open


    SUBROUTINE unpack_restricted_buffer_n(kk, jj, ii, fc, buffer, tasksize, &
            istart, istop, jstart, jstop, kstart, kstop, ipos, jpos, kpos)
        !$omp declare target
        ! Subroutine arguments
        INTEGER(intk), INTENT(in) :: kk, jj, ii, tasksize
        REAL(realk), INTENT(inout) :: fc(kk, jj, ii)
        REAL(realk), INTENT(in) :: buffer(tasksize)
        INTEGER(intk), INTENT(in) :: istart, istop, jstart, jstop, kstart, kstop
        INTEGER(intk), INTENT(in) :: ipos, jpos, kpos

        ! Local variables
        INTEGER(intk) :: i, j, k, ic, jc, kc, icount
        INTEGER(intk) :: ni, nj, nk

        ni = (istop - istart)/2 + 1
        nj = (jstop - jstart)/2 + 1
        nk = (kstop - kstart)/2 + 1

        !$omp do collapse(3) private(ic, jc, kc, icount)
        DO i = istart, istop, 2
            DO j = jstart, jstop, 2
                DO k = kstart, kstop, 2
                    ic = ipos + (i-2)/2 - 1
                    jc = jpos + (j-2)/2 - 1
                    kc = kpos + (k-2)/2 - 1
                    icount = ((i-istart)/2*nj + (j-jstart)/2)*nk + &
                        (k-kstart)/2 + 1
                    IF (ABS(fc(kc, jc, ic)) < TINY(1.0_realk) .AND. &
                            buffer(icount) >= 1.0_realk) THEN
                        fc(kc, jc, ic) = buffer(icount)
                    END IF
                END DO
            END DO
        END DO
        !$omp end do

#ifdef _MGLET_DEBUG_
        IF (ni*nj*nk /= tasksize) CALL errr(__FILE__, __LINE__)
#endif
    END SUBROUTINE unpack_restricted_buffer_n


    SUBROUTINE process_mpisend(nmpistasks, mpistasks)
        ! Subroutine arguments
        INTEGER(intk), INTENT(in) :: nmpistasks
        INTEGER(intk), INTENT(in) :: mpistasks(mpitasksize, nmpistasks+1)

        ! Local variables
        INTEGER(intk) :: itask
        INTEGER(int32) :: iprocnbr, messagelength, sendcounter

#ifdef _MGLET_PROFILE_ANNOTATIONS_
        CALL profile_range_push("process_mpisend")
#endif

        DO itask = 1, nmpistasks
            iprocnbr = INT(mpistasks(1, itask), kind=int32)
            messagelength = INT(mpistasks(2, itask), kind=int32)
            sendcounter = INT(mpistasks(3, itask), kind=int32)

            !$omp target data use_device_addr(sendbuf)
            CALL MPI_Isend(sendbuf(sendcounter+1), messagelength, &
                mglet_mpi_real, iprocnbr, 1, MPI_COMM_WORLD, sendreqs(itask))
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


    SUBROUTINE process_sendtasks(nstasks, stasks)
        ! Subroutine arguments
        INTEGER(intk), INTENT(in) :: nstasks
        INTEGER(intk), INTENT(in) :: stasks(sendtasksize, nstasks+1)

        ! Local variables
        TYPE(field_t), POINTER :: ddx, ddy, ddz, bp, bt
        LOGICAL :: foundbt

        IF (nstasks == 0) RETURN

        CALL get_field(ddx, "DDX")
        CALL get_field(ddy, "DDY")
        CALL get_field(ddz, "DDZ")

        SELECT CASE (ib%type)
        CASE ("NONE")
            CALL process_sendtasks_noib(nstasks, stasks, ddx, ddy, ddz)
        CASE ("GHOSTCELL")
            CALL get_field(bp, "BP")
            CALL get_field(bt, "BT", foundbt)
            IF (.NOT. foundbt) THEN
                ! For non-scalar cases, point bt to bp. It is unused anyways.
                bt => bp
            END IF
            CALL process_sendtasks_gc(nstasks, stasks, ddx, ddy, ddz, bp, bt)
        CASE DEFAULT
            CALL errr(__FILE__, __LINE__)
        END SELECT

        ! Safety check based on final dummy entry
        IF (.NOT. ALL(stasks(:, nstasks+1) == -1)) THEN
            WRITE(*, *) "Did not encounter the expected dummy task."
            CALL errr(__FILE__, __LINE__)
        END IF
    END SUBROUTINE process_sendtasks


    SUBROUTINE process_sendtasks_noib(nstasks, stasks, ddx_f, ddy_f, ddz_f)
        ! Subroutine arguments
        INTEGER(intk), INTENT(in) :: nstasks
        INTEGER(intk), INTENT(in) :: stasks(sendtasksize, nstasks+1)
        TYPE(field_t), INTENT(in) :: ddx_f, ddy_f, ddz_f

        ! Local variables
        ! none...

#ifdef _MGLET_PROFILE_ANNOTATIONS_
        CALL profile_range_push("process_sendtasks_noib")
#endif

        CALL process_sendtasks_noib_impl(nstasks, stasks, f1%arr, f2%arr, &
            f3%arr, f4%arr, f5%arr, f6%arr, ddx_f%arr, ddy_f%arr, &
            ddz_f%arr)

#ifdef _MGLET_PROFILE_ANNOTATIONS_
        CALL profile_range_pop()
#endif
    END SUBROUTINE process_sendtasks_noib


    SUBROUTINE process_sendtasks_noib_impl(nstasks, stasks, a1, a2, a3, a4, &
            a5, a6, ddx, ddy, ddz)
        ! Subroutine arguments
        INTEGER(intk), INTENT(in) :: nstasks
        INTEGER(intk), INTENT(in) :: stasks(sendtasksize, nstasks+1)
        REAL(realk), INTENT(in) :: a1(*), a2(*), a3(*), a4(*), a5(*), a6(*)
        REAL(realk), INTENT(in) :: ddx(*), ddy(*), ddz(*)

        ! Local variables
        INTEGER(intk) :: itask, kk, jj, ii, ip3, ipx, ipy, ipz
        INTEGER(intk) :: fieldid, igrid, istart, istop, jstart, jstop, kstart, &
            kstop
        INTEGER(int32) :: icount, tasksize
        CHARACTER(len=1) :: flag

        !$omp target teams distribute private(itask, fieldid, flag, icount, &
        !$omp& tasksize, igrid, istart, istop, jstart, jstop, kstart, kstop, &
        !$omp& kk, jj, ii, ip3, ipx, ipy, ipz)
        DO itask = 1, nstasks
            fieldid = stasks(1, itask)
            flag = ACHAR(stasks(2, itask))
            icount = INT(stasks(3, itask), kind=int32)
            tasksize = INT(stasks(4, itask), kind=int32)
            igrid = stasks(5, itask)
            istart = stasks(6, itask)
            istop = stasks(7, itask)
            jstart = stasks(8, itask)
            jstop = stasks(9, itask)
            kstart = stasks(10, itask)
            kstop = stasks(11, itask)

            CALL get_mgdims(kk, jj, ii, igrid)
            CALL get_ip3(ip3, igrid)
            CALL get_ip1x(ipx, igrid)
            CALL get_ip1y(ipy, igrid)
            CALL get_ip1z(ipz, igrid)

            !$omp parallel
            SELECT CASE(fieldid)
            CASE (1)
                CALL restrict_noib_flag(flag, kk, jj, ii, a1(ip3), &
                    ddx(ipx), ddy(ipy), ddz(ipz), icount, tasksize, &
                    istart, istop, jstart, jstop, kstart, kstop)
            CASE (2)
                CALL restrict_noib_flag(flag, kk, jj, ii, a2(ip3), &
                    ddx(ipx), ddy(ipy), ddz(ipz), icount, tasksize, &
                    istart, istop, jstart, jstop, kstart, kstop)
            CASE (3)
                CALL restrict_noib_flag(flag, kk, jj, ii, a3(ip3), &
                    ddx(ipx), ddy(ipy), ddz(ipz), icount, tasksize, &
                    istart, istop, jstart, jstop, kstart, kstop)
            CASE (4)
                CALL restrict_noib_flag(flag, kk, jj, ii, a4(ip3), &
                    ddx(ipx), ddy(ipy), ddz(ipz), icount, tasksize, &
                    istart, istop, jstart, jstop, kstart, kstop)
            CASE (5)
                CALL restrict_noib_flag(flag, kk, jj, ii, a5(ip3), &
                    ddx(ipx), ddy(ipy), ddz(ipz), icount, tasksize, &
                    istart, istop, jstart, jstop, kstart, kstop)
            CASE (6)
                CALL restrict_noib_flag(flag, kk, jj, ii, a6(ip3), &
                    ddx(ipx), ddy(ipy), ddz(ipz), icount, tasksize, &
                    istart, istop, jstart, jstop, kstart, kstop)
#ifdef _MGLET_DEBUG_
            CASE DEFAULT
                CALL errr(__FILE__, __LINE__)
#endif
            END SELECT
            !$omp end parallel
        END DO
        !$omp end target teams distribute
    END SUBROUTINE process_sendtasks_noib_impl


    SUBROUTINE process_sendtasks_gc(nstasks, stasks, ddx_f, ddy_f, ddz_f, &
            bp_f, bt_f)
        ! Subroutine arguments
        INTEGER(intk), INTENT(in) :: nstasks
        INTEGER(intk), INTENT(in) :: stasks(sendtasksize, nstasks+1)
        TYPE(field_t), INTENT(in) :: ddx_f, ddy_f, ddz_f, bp_f, bt_f

        ! Local variables
        ! none...

#ifdef _MGLET_PROFILE_ANNOTATIONS_
        CALL profile_range_push("process_sendtasks_gc")
#endif

        CALL process_sendtasks_gc_impl(nstasks, stasks, f1%arr, f2%arr, &
            f3%arr, f4%arr, f5%arr, f6%arr, ddx_f%arr, ddy_f%arr, &
            ddz_f%arr, bp_f%arr, bt_f%arr)

#ifdef _MGLET_PROFILE_ANNOTATIONS_
        CALL profile_range_pop()
#endif
    END SUBROUTINE process_sendtasks_gc


    SUBROUTINE process_sendtasks_gc_impl(nstasks, stasks, a1, a2, a3, a4, &
            a5, a6, ddx, ddy, ddz, bp, bt)
        ! Subroutine arguments
        INTEGER(intk), INTENT(in) :: nstasks
        INTEGER(intk), INTENT(in) :: stasks(sendtasksize, nstasks+1)
        REAL(realk), INTENT(in) :: a1(*), a2(*), a3(*), a4(*), a5(*), a6(*)
        REAL(realk), INTENT(in) :: ddx(*), ddy(*), ddz(*), bp(*), bt(*)

        ! Local variables
        INTEGER(intk) :: itask, kk, jj, ii, ip3, ipx, ipy, ipz
        INTEGER(intk) :: fieldid, igrid, istart, istop, jstart, jstop, kstart, &
            kstop
        INTEGER(int32) :: icount, tasksize
        CHARACTER(len=1) :: flag

        !$omp target teams distribute private(itask, fieldid, flag, icount, &
        !$omp& tasksize, igrid, istart, istop, jstart, jstop, kstart, kstop, &
        !$omp& kk, jj, ii, ip3, ipx, ipy, ipz)
        DO itask = 1, nstasks
            fieldid = stasks(1, itask)
            flag = ACHAR(stasks(2, itask))
            icount = INT(stasks(3, itask), kind=int32)
            tasksize = INT(stasks(4, itask), kind=int32)
            igrid = stasks(5, itask)
            istart = stasks(6, itask)
            istop = stasks(7, itask)
            jstart = stasks(8, itask)
            jstop = stasks(9, itask)
            kstart = stasks(10, itask)
            kstop = stasks(11, itask)

            CALL get_mgdims(kk, jj, ii, igrid)
            CALL get_ip3(ip3, igrid)
            CALL get_ip1x(ipx, igrid)
            CALL get_ip1y(ipy, igrid)
            CALL get_ip1z(ipz, igrid)

            !$omp parallel
            SELECT CASE(fieldid)
            CASE (1)
                CALL restrict_gc_flag(flag, kk, jj, ii, a1(ip3), &
                    ddx(ipx), ddy(ipy), ddz(ipz), bp(ip3), bt(ip3), icount, &
                    tasksize, istart, istop, jstart, jstop, kstart, kstop)
            CASE (2)
                CALL restrict_gc_flag(flag, kk, jj, ii, a2(ip3), &
                    ddx(ipx), ddy(ipy), ddz(ipz), bp(ip3), bt(ip3), icount, &
                    tasksize, istart, istop, jstart, jstop, kstart, kstop)
            CASE (3)
                CALL restrict_gc_flag(flag, kk, jj, ii, a3(ip3), &
                    ddx(ipx), ddy(ipy), ddz(ipz), bp(ip3), bt(ip3), icount, &
                    tasksize, istart, istop, jstart, jstop, kstart, kstop)
            CASE (4)
                CALL restrict_gc_flag(flag, kk, jj, ii, a4(ip3), &
                    ddx(ipx), ddy(ipy), ddz(ipz), bp(ip3), bt(ip3), icount, &
                    tasksize, istart, istop, jstart, jstop, kstart, kstop)
            CASE (5)
                CALL restrict_gc_flag(flag, kk, jj, ii, a5(ip3), &
                    ddx(ipx), ddy(ipy), ddz(ipz), bp(ip3), bt(ip3), icount, &
                    tasksize, istart, istop, jstart, jstop, kstart, kstop)
            CASE (6)
                CALL restrict_gc_flag(flag, kk, jj, ii, a6(ip3), &
                    ddx(ipx), ddy(ipy), ddz(ipz), bp(ip3), bt(ip3), icount, &
                    tasksize, istart, istop, jstart, jstop, kstart, kstop)
#ifdef _MGLET_DEBUG_
            CASE DEFAULT
                CALL errr(__FILE__, __LINE__)
#endif
            END SELECT
            !$omp end parallel
        END DO
        !$omp end target teams distribute
    END SUBROUTINE process_sendtasks_gc_impl


    SUBROUTINE process_mpirecv(nmpirtasks, mpirtasks)
        ! Subroutine arguments
        INTEGER(intk), INTENT(in) :: nmpirtasks
        INTEGER(intk), INTENT(in) :: mpirtasks(mpitasksize, nmpirtasks+1)

        ! Local variables
        INTEGER(intk) :: itask
        INTEGER(int32) :: iprocnbr, messagelength, recvcounter

#ifdef _MGLET_PROFILE_ANNOTATIONS_
        CALL profile_range_push("process_mpirecv")
#endif

        DO itask = 1, nmpirtasks
            iprocnbr = INT(mpirtasks(1, itask), kind=int32)
            messagelength = INT(mpirtasks(2, itask), kind=int32)
            recvcounter = INT(mpirtasks(3, itask), kind=int32)

            !$omp target data use_device_addr(recvbuf)
            CALL MPI_Irecv(recvbuf(recvcounter+1), messagelength, &
                mglet_mpi_real, iprocnbr, 1, MPI_COMM_WORLD, recvreqs(itask))
            !$omp end target data
        END DO

        nrecv = nmpirtasks

        ! Safety check based on final dummy entry
        IF (.NOT. ALL(mpirtasks(:, nmpirtasks+1) == -1)) THEN
            WRITE(*, *) "Did not encounter the expected dummy task."
            CALL errr(__FILE__, __LINE__)
        END IF

#ifdef _MGLET_PROFILE_ANNOTATIONS_
        CALL profile_range_pop()
#endif
    END SUBROUTINE process_mpirecv


    SUBROUTINE prepare_mpirecvtasks(mpirtasks, nmpirtasks, ilevel, flag, &
            v1, v2, v3, s1, s2, s3)
        ! Subroutine arguments
        INTEGER(intk), INTENT(inout) :: mpirtasks(mpitasksize, maxmpirecvtasks)
        INTEGER(intk), INTENT(out) :: nmpirtasks
        INTEGER(intk), INTENT(in) :: ilevel
        CHARACTER(len=1), INTENT(in) :: flag
        TYPE(field_t), OPTIONAL, INTENT(inout) :: v1, v2, v3, s1, s2, s3

        ! Local variables
        INTEGER(intk) :: i, iprocnbr, igridf, impirtask
        INTEGER(int32) :: ncells, messagelength, recvcounter

        recvcounter = 0
        messagelength = 0
        nrecv = 0
        recvidxlist = -HUGE(1_intk)
        recvlist = 0
        impirtask = 0

        DO i = 1, irecv
            iprocnbr = recvconns(1, i)  ! The sender process (fine side)
            IF (iprocnbr == myid) CYCLE

            igridf = recvconns(3, i)
            IF (ilevel == level(igridf)) THEN
                CALL count_ncells(ncells, flag, igridf, v1, v2, v3, s1, s2, s3)
                recvidxlist(1, i) = iprocnbr
                recvidxlist(2, i) = ncells
                recvidxlist(3, i) = recvcounter + messagelength
                messagelength = messagelength + ncells

                IF (recvcounter + messagelength > idim_mg_bufs) THEN
                    CALL errr(__FILE__, __LINE__)
                END IF
            END IF

            IF (messagelength > 0) THEN
                IF (i == irecv) THEN
                    CALL add_mpi_task(mpirtasks, impirtask, iprocnbr, &
                        messagelength, recvcounter, 'recv')
                ELSE IF (recvconns(1, i + 1) /= iprocnbr) THEN
                    CALL add_mpi_task(mpirtasks, impirtask, iprocnbr, &
                        messagelength, recvcounter, 'recv')
                END IF
            END IF
        END DO

        nmpirtasks = impirtask
        ! Add a harmful dummy task at (ntasks+1) for checking
        mpirtasks(:, nmpirtasks+1) = -1
    END SUBROUTINE prepare_mpirecvtasks


    SUBROUTINE prepare_tasks_all(stasks, nstasks, etasks, netasks, mpistasks, &
        nmpistasks, ilevel, flag, v1, v2, v3, s1, s2, s3)
        ! Subroutine arguments
        INTEGER(intk), INTENT(inout) :: stasks(sendtasksize, maxsendtasks)
        INTEGER(intk), INTENT(out) :: nstasks
        INTEGER(intk), INTENT(inout) :: etasks(selftasksize, maxselftasks)
        INTEGER(intk), INTENT(out) :: netasks
        INTEGER(intk), INTENT(inout) :: mpistasks(mpitasksize, maxmpisendtasks)
        INTEGER(intk), INTENT(out) :: nmpistasks
        INTEGER(intk), INTENT(in) :: ilevel
        CHARACTER(len=1), INTENT(in) :: flag
        TYPE(field_t), OPTIONAL, INTENT(inout) :: v1, v2, v3, s1, s2, s3

        ! Local variables
        INTEGER(intk) :: i, iprocnbr, igridf, istask, ietask, impistask
        INTEGER(int32) :: ncells, ncells_total, messagelength, sendcounter, &
            selfcounter

#ifdef _MGLET_PROFILE_ANNOTATIONS_
        CALL profile_range_push("prepare_tasks_all")
#endif

        ! Pack all buffers and send data
        sendcounter = 0
        CALL count_sendcells(ncells_total, ilevel, flag, v1, v2, v3, &
            s1, s2, s3)
        selfcounter = 1 + ncells_total
        messagelength = 0
        nsend = 0

        ! Initialize the task counters to zero
        istask = 0
        ietask = 0
        impistask = 0

        DO i = 1, isend
            iprocnbr = sendconns(2, i)
            igridf = sendconns(3, i)

            IF (ilevel == level(igridf)) THEN
                IF (iprocnbr == myid) THEN
                    CALL prepare_selftask(etasks, ietask, i, flag, &
                        selfcounter, v1, v2, v3, s1, s2, s3)
                ELSE
                    CALL prepare_sendtask(stasks, istask, i, &
                        messagelength, sendcounter, flag, &
                        v1, v2, v3, s1, s2, s3)
                    CALL count_ncells(ncells, flag, igridf, &
                        v1, v2, v3, s1, s2, s3)
                    ! ncells already sums over all args (nvars not needed)
                    messagelength = messagelength + ncells
                END IF
            END IF

            IF (messagelength > 0) THEN
                IF (i == isend) THEN
                    CALL add_mpi_task(mpistasks, impistask, iprocnbr, &
                        messagelength, sendcounter, 'send')
                ELSE IF (sendconns(2, i + 1) /= iprocnbr) THEN
                    CALL add_mpi_task(mpistasks, impistask, iprocnbr, &
                        messagelength, sendcounter, 'send')
                END IF
            END IF
        END DO

        ! Set the output task counters
        nstasks = istask
        netasks = ietask
        nmpistasks = impistask

        IF (sendcounter >= selfcounter) CALL errr(__FILE__, __LINE__)
        IF (selfcounter > idim_mg_bufs + 1) CALL errr(__FILE__, __LINE__)

        ! Add a harmful dummy task at (ntasks+1) for checking
        stasks(:, nstasks+1) = -1
        etasks(:, netasks+1) = -1
        mpistasks(:, nmpistasks+1) = -1

#ifdef _MGLET_PROFILE_ANNOTATIONS_
        CALL profile_range_pop()
#endif
    END SUBROUTINE prepare_tasks_all


    SUBROUTINE prepare_selftask(etasks, ietask, sendid, flag, selfcounter, &
            v1, v2, v3, s1, s2, s3)
        ! Subroutine arguments
        INTEGER(intk), INTENT(inout) :: etasks(selftasksize, maxselftasks)
        INTEGER(intk), INTENT(inout) :: ietask
        INTEGER(intk), INTENT(in) :: sendid
        CHARACTER(len=1), INTENT(in) :: flag
        INTEGER(int32), INTENT(inout) :: selfcounter
        TYPE(field_t), OPTIONAL, INTENT(inout) :: v1, v2, v3, s1, s2, s3

        ! Local variables
        INTEGER(intk) :: igridf, igridc
        INTEGER(int32) :: ncells, icount
        CHARACTER(len=1) :: flag_u

        igridf = sendconns(3, sendid)
        igridc = sendconns(4, sendid)

        CALL count_ncells(ncells, flag, igridf, &
            v1, v2, v3, s1, s2, s3)

        icount = 0

        IF (PRESENT(v1)) THEN
            IF (flag == '*') THEN
                flag_u = 'U'
            ELSE
                flag_u = flag
            END IF

            ietask = ietask + 1
            CALL add_single_selftask(etasks(:, ietask), 1, flag_u, icount, &
                selfcounter, igridf, igridc)
        END IF

        IF (PRESENT(v2)) THEN
            ietask = ietask + 1
            CALL add_single_selftask(etasks(:, ietask), 2, 'V', icount, &
                selfcounter, igridf, igridc)
        END IF

        IF (PRESENT(v3)) THEN
            ietask = ietask + 1
            CALL add_single_selftask(etasks(:, ietask), 3, 'W', icount, &
                selfcounter, igridf, igridc)
        END IF

        IF (PRESENT(s1)) THEN
            ietask = ietask + 1
            CALL add_single_selftask(etasks(:, ietask), 4, 'P', icount, &
                selfcounter, igridf, igridc)
        END IF

        IF (PRESENT(s2)) THEN
            ietask = ietask + 1
            CALL add_single_selftask(etasks(:, ietask), 5, 'P', icount, &
                selfcounter, igridf, igridc)
        END IF

        IF (PRESENT(s3)) THEN
            ietask = ietask + 1
            CALL add_single_selftask(etasks(:, ietask), 6, 'P', icount, &
                selfcounter, igridf, igridc)
        END IF

        IF (ncells /= icount) THEN
            WRITE(*, *) "ncells:", ncells, "self_len:", icount
            CALL errr(__FILE__, __LINE__)
        END IF
    END SUBROUTINE prepare_selftask


    SUBROUTINE prepare_sendtask(stasks, istask, sendid, messagelength, &
            sendcounter, flag, v1, v2, v3, s1, s2, s3)
        ! Subroutine arguments
        INTEGER(intk), INTENT(inout) :: stasks(sendtasksize, maxsendtasks)
        INTEGER(intk), INTENT(inout) :: istask
        INTEGER(intk), INTENT(in) :: sendid
        INTEGER(int32), INTENT(in) :: messagelength
        INTEGER(int32), INTENT(inout) :: sendcounter
        CHARACTER(len=1), INTENT(in) :: flag
        TYPE(field_t), OPTIONAL, INTENT(inout) :: v1, v2, v3, s1, s2, s3

        ! Local variables
        INTEGER(intk) :: igrid
        INTEGER(int32) :: ncells, offset, icount
        CHARACTER(len=1) :: flag_u

        igrid = sendconns(3, sendid)

        CALL count_ncells(ncells, flag, igrid, &
            v1, v2, v3, s1, s2, s3)

        ! Check that buffer does not overflow
        IF (sendcounter + messagelength + ncells > idim_mg_bufs) THEN
            CALL errr(__FILE__, __LINE__)
        END IF

        ! Reset message size counter
        offset = sendcounter + messagelength + 1
        icount = offset

        IF (PRESENT(v1)) THEN
            IF (flag == '*') THEN
                flag_u = 'U'
            ELSE
                flag_u = flag
            END IF

            istask = istask + 1
            CALL add_single_sendtask(stasks(:, istask), 1, flag_u, icount, &
                igrid)
        END IF

        IF (PRESENT(v2)) THEN
            istask = istask + 1
            CALL add_single_sendtask(stasks(:, istask), 2, 'V', icount, igrid)
        END IF

        IF (PRESENT(v3)) THEN
            istask = istask + 1
            CALL add_single_sendtask(stasks(:, istask), 3, 'W', icount, igrid)
        END IF

        IF (PRESENT(s1)) THEN
            istask = istask + 1
            CALL add_single_sendtask(stasks(:, istask), 4, 'P', icount, igrid)
        END IF

        IF (PRESENT(s2)) THEN
            istask = istask + 1
            CALL add_single_sendtask(stasks(:, istask), 5, 'P', icount, igrid)
        END IF

        IF (PRESENT(s3)) THEN
            istask = istask + 1
            CALL add_single_sendtask(stasks(:, istask), 6, 'P', icount, igrid)
        END IF

        IF (ncells /= (icount - offset)) THEN
            WRITE(*, *) "ncells:", ncells, "packed_len:", icount - offset
            CALL errr(__FILE__, __LINE__)
        END IF
    END SUBROUTINE prepare_sendtask


    SUBROUTINE add_single_selftask(task, fieldid, flag, icount, selfcounter, &
            igridf, igridc)
        ! Subroutine arguments
        INTEGER(intk), INTENT(out) :: task(selftasksize)
        INTEGER(intk), INTENT(in) :: fieldid
        CHARACTER(len=1), INTENT(in) :: flag
        INTEGER(int32), INTENT(inout) :: icount
        INTEGER(int32), INTENT(inout) :: selfcounter
        INTEGER(intk), INTENT(in) :: igridf, igridc

        ! Local variables
        INTEGER(intk) :: flagidx, istart, istop, jstart, jstop, kstart, kstop
        INTEGER(int32) :: tasksize

        flagidx = IACHAR(flag)

        CALL message_length(tasksize, flag, igridf)
        CALL start_and_stop(istart, istop, jstart, jstop, kstart, kstop, &
            flag, igridf)

        IF (selfcounter + tasksize - 1 > idim_mg_bufs) THEN
            CALL errr(__FILE__, __LINE__)
        END IF

        task(1) = fieldid
        task(2) = flagidx
        task(3) = INT(tasksize, kind=intk)
        task(4) = INT(selfcounter, kind=intk)
        task(5) = igridf
        task(6) = igridc
        task(7) = istart
        task(8) = istop
        task(9) = jstart
        task(10) = jstop
        task(11) = kstart
        task(12) = kstop
        task(13) = iposition(igridf)
        task(14) = jposition(igridf)
        task(15) = kposition(igridf)

        icount = icount + tasksize
        selfcounter = selfcounter + tasksize
    END SUBROUTINE add_single_selftask


    SUBROUTINE count_sendcells(ncells_total, ilevel, flag, v1, v2, v3, &
            s1, s2, s3)
        ! Subroutine arguments
        INTEGER(int32), INTENT(out) :: ncells_total
        INTEGER(intk), INTENT(in) :: ilevel
        CHARACTER(len=1), INTENT(in) :: flag
        TYPE(field_t), OPTIONAL, INTENT(inout) :: v1, v2, v3, s1, s2, s3

        ! Local variables
        INTEGER(intk) :: i, iprocnbr, igridf
        INTEGER(int32) :: ncells

        ncells_total = 0
        DO i = 1, isend
            iprocnbr = sendconns(2, i)
            IF (iprocnbr == myid) CYCLE

            igridf = sendconns(3, i)
            IF (ilevel == level(igridf)) THEN
                CALL count_ncells(ncells, flag, igridf, &
                    v1, v2, v3, s1, s2, s3)
                ncells_total = ncells_total + ncells
            END IF
        END DO
    END SUBROUTINE count_sendcells


    SUBROUTINE add_single_sendtask(task, fieldid, flag, icount, igrid)
        ! Subroutine arguments
        INTEGER(intk), INTENT(out) :: task(sendtasksize)
        INTEGER(intk), INTENT(in) :: fieldid
        CHARACTER(len=1), INTENT(in) :: flag
        INTEGER(int32), INTENT(inout) :: icount
        INTEGER(intk), INTENT(in) :: igrid

        ! Local variables
        INTEGER(intk) :: flagidx, istart, istop, jstart, jstop, kstart, kstop
        INTEGER(int32) :: tasksize

        flagidx = IACHAR(flag)

        CALL message_length(tasksize, flag, igrid)
        CALL start_and_stop(istart, istop, jstart, jstop, kstart, kstop, flag, &
            igrid)

        task(1) = fieldid
        task(2) = flagidx
        task(3) = INT(icount, kind=intk)
        task(4) = INT(tasksize, kind=intk)
        task(5) = igrid
        task(6) = istart
        task(7) = istop
        task(8) = jstart
        task(9) = jstop
        task(10) = kstart
        task(11) = kstop

        icount = icount + tasksize
    END SUBROUTINE add_single_sendtask


    SUBROUTINE add_single_recvtask(task, fieldid, flag, icount, igridf, igridc)
        ! Subroutine arguments
        INTEGER(intk), INTENT(out) :: task(recvtasksize)
        INTEGER(intk), INTENT(in) :: fieldid
        CHARACTER(len=1), INTENT(in) :: flag
        INTEGER(int32), INTENT(inout) :: icount
        INTEGER(intk), INTENT(in) :: igridf, igridc

        ! Local variables
        INTEGER(intk) :: flagidx, istart, istop, jstart, jstop, kstart, kstop
        INTEGER(int32) :: tasksize

        flagidx = IACHAR(flag)

        CALL message_length(tasksize, flag, igridf)
        CALL start_and_stop(istart, istop, jstart, jstop, kstart, kstop, &
            flag, igridf)

        task(1) = fieldid
        task(2) = flagidx
        task(3) = INT(tasksize, kind=intk)
        task(4) = INT(icount, kind=intk)
        task(5) = igridf
        task(6) = igridc
        task(7) = istart
        task(8) = istop
        task(9) = jstart
        task(10) = jstop
        task(11) = kstart
        task(12) = kstop
        task(13) = iposition(igridf)
        task(14) = jposition(igridf)
        task(15) = kposition(igridf)

        icount = icount + tasksize
    END SUBROUTINE add_single_recvtask


    SUBROUTINE add_mpi_task(mpitasks, impitask, iprocnbr, messagelength, &
            counter, type)
        ! Subroutine arguments
        INTEGER(intk), INTENT(inout) :: mpitasks(mpitasksize, *)
        INTEGER(intk), INTENT(inout) :: impitask
        INTEGER(intk), INTENT(in) :: iprocnbr
        INTEGER(int32), INTENT(inout) :: messagelength
        INTEGER(int32), INTENT(inout) :: counter
        CHARACTER(len=4), INTENT(in) :: type

        ! Local variables
        ! none...

        IF (type == 'send') THEN
            nsend = nsend + 1
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


    SUBROUTINE init_ftoc2()
        ! Subroutine arguments
        ! none...

        ! Local variables
        INTEGER(intk) :: nselfsend, nselfrecv

        IF (is_init) CALL errr(__FILE__, __LINE__)

        ! Check if ctof_core_mod provides necessary infrastructure
        IF (.NOT. is_ftoc_core_init) CALL errr(__FILE__, __LINE__)

        ALLOCATE(recvidxlist(3, irecv), source=0_intk)
        ALLOCATE(recvlist(irecv), source=0_int32)
        ALLOCATE(sendreqs(numprocs), source=MPI_REQUEST_NULL)
        ALLOCATE(recvreqs(numprocs), source=MPI_REQUEST_NULL)
        nsend = 0
        nrecv = 0

        ! Allocating the workpackage arrays
        ! Always add 1 extra dummy task for error checking purposes
        maxsendtasks = 6 * isend + 1
        maxrecvtasks = 6 * irecv + 1
        ALLOCATE(sendtasks(sendtasksize, maxsendtasks))
        ALLOCATE(recvtasks(recvtasksize, maxrecvtasks))
        !$omp target enter data map(always, to: sendtasks, recvtasks)

        CALL count_selftasks(nselfsend, nselfrecv)
        maxselftasks = 6 * (nselfsend + nselfrecv) + 1
        ALLOCATE(selftasks(selftasksize, maxselftasks))
        !$omp target enter data map(always, to: selftasks)

        ! MPI tasks are only used on the host and do not exceed isend/irecv+1
        maxmpisendtasks = isend + 1
        maxmpirecvtasks = irecv + 1
        ALLOCATE(mpisendtasks(mpitasksize, maxmpisendtasks))
        ALLOCATE(mpirecvtasks(mpitasksize, maxmpirecvtasks))

        ! Allocate the workrecords array for all possible types of parent
        ! dimension 1 = ilevel
        ! dimension 2 = argument combination (6 bits = 2^6-1 = 63)
        ! dimension 3 = flag (any ASCII char from '*'=42 to 'Z'=90)
        ALLOCATE(workrecords(minlevel:maxlevel, 1:63, IACHAR('*'):IACHAR('Z')))

        CALL check_sendbuf_capacity()

        is_init = .TRUE.

        ! Record relevant calls for later efficient reuse on device
        ! (includes a parity check with the CPU version)
        CALL run_recording_pass()

    END SUBROUTINE init_ftoc2


    SUBROUTINE check_sendbuf_capacity()
        ! Subroutine arguments
        ! none...

        ! Local variables
        INTEGER(intk) :: i, iprocnbr, igridf
        INTEGER(int32) :: sendcells, selfcells, ncells, messagelength

        sendcells = 0
        selfcells = 0

        DO i = 1, isend
            iprocnbr = sendconns(2, i)
            igridf = sendconns(3, i)

            ncells = 0
            CALL message_length(messagelength, 'U', igridf)
            ncells = ncells + messagelength
            CALL message_length(messagelength, 'V', igridf)
            ncells = ncells + messagelength
            CALL message_length(messagelength, 'W', igridf)
            ncells = ncells + messagelength
            CALL message_length(messagelength, 'P', igridf)
            ncells = ncells + 3*messagelength

            IF (iprocnbr == myid) THEN
                selfcells = selfcells + ncells
            ELSE
                sendcells = sendcells + ncells
            END IF
        END DO

        IF (sendcells + selfcells > idim_mg_bufs) THEN
            WRITE(*, *) "ftoc2 sendbuf too small for send+self upper bound."
            WRITE(*, *) "send cells:", sendcells, " self cells:", selfcells
            WRITE(*, *) "required:", sendcells+selfcells, &
                " available:", idim_mg_bufs
            CALL errr(__FILE__, __LINE__)
        END IF
    END SUBROUTINE check_sendbuf_capacity


    SUBROUTINE run_recording_pass()
        ! Subroutine arguments
        ! none...

        ! Local variables
        TYPE(field_t) :: udummy, vdummy, wdummy, pdummy
        INTEGER(intk) :: ilevel

        is_recording = .TRUE.

        CALL udummy%init("U_DUMMY", istag=1)
        CALL vdummy%init("V_DUMMY", jstag=1)
        CALL wdummy%init("W_DUMMY", kstag=1)
        CALL pdummy%init("P_DUMMY")

        CALL init_dummy_fields_cpu(udummy, vdummy, wdummy, pdummy)

        CALL connect(layers=2, v1=udummy, v2=vdummy, v3=wdummy, s1=pdummy, &
            corners=.TRUE.)

        !$omp target enter data map(to: udummy%arr, vdummy%arr, wdummy%arr, &
        !$omp& pdummy%arr)

        ! START -- This section defines the record variants of ftoc2 ---

        DO ilevel = maxlevel, minlevel, -1
            ! Outer pressuresolver iterations
            CALL ftoc2(ilevel, pdummy, pdummy, flag='P')
            ! Pre pressuresolver iterations in mgpoisl
            CALL ftoc2(ilevel, pdummy, pdummy, flag='R')
            ! Post pressuresolver iterations in mgpoisl
            CALL ftoc2(ilevel, udummy, vdummy, wdummy, pdummy)
        END DO

        ! Repeating the same calls with the CPU function !
        ! ("ftoc1" is already initialized at this point)

        DO ilevel = maxlevel, minlevel, -1
            ! Outer pressuresolver iterations
            CALL ftoc1(ilevel, pdummy, pdummy, flag='P')
            ! Pre pressuresolver iterations in mgpoisl
            CALL ftoc1(ilevel, pdummy, pdummy, flag='R')
            ! Post pressuresolver iterations in mgpoisl
            CALL ftoc1(ilevel, udummy, vdummy, wdummy, pdummy)
        END DO

        ! END -- This section defines the record variants of ftoc2 ---

        CALL assert_same_field(udummy)
        CALL assert_same_field(vdummy)
        CALL assert_same_field(wdummy)
        CALL assert_same_field(pdummy)

        !$omp target exit data map(delete: udummy%arr, vdummy%arr, wdummy%arr, &
        !$omp& pdummy%arr)

        CALL udummy%finish()
        CALL vdummy%finish()
        CALL wdummy%finish()
        CALL pdummy%finish()

        is_recording = .FALSE.
    END SUBROUTINE run_recording_pass


    SUBROUTINE init_dummy_fields_cpu(u_f, v_f, w_f, p_f)
        ! Subroutine arguments
        TYPE(field_t), INTENT(inout) :: u_f, v_f, w_f, p_f
        ! Local variables
        INTEGER(intk) :: i
        DO i = 1, SIZE(u_f%arr)
            u_f%arr(i) = REAL(MOD(17_intk*i + 11_intk, 997_intk), realk) * &
                1.0e-3_realk
        END DO
        DO i = 1, SIZE(v_f%arr)
            v_f%arr(i) = REAL(MOD(23_intk*i + 7_intk, 991_intk), realk) * &
                1.0e-3_realk
        END DO
        DO i = 1, SIZE(w_f%arr)
            w_f%arr(i) = REAL(MOD(31_intk*i + 5_intk, 983_intk), realk) * &
                1.0e-3_realk
        END DO
        DO i = 1, SIZE(p_f%arr)
            p_f%arr(i) = REAL(MOD(47_intk*i + 3_intk, 977_intk), realk) * &
                1.0e-3_realk
        END DO
    END SUBROUTINE init_dummy_fields_cpu


    SUBROUTINE assert_same_field(field_f)
        ! Subroutine arguments
        TYPE(field_t), INTENT(in) :: field_f
        ! Local variables
        REAL(realk), ALLOCATABLE :: arr(:)
        REAL(realk) :: diff_local, diff_global
        INTEGER(intk) :: ierr
        LOGICAL :: is_same = .TRUE.

        ! Store the field on the host into antoher array
        ALLOCATE(arr(SIZE(field_f%arr)))
        arr = field_f%arr

        ! Update the field from the device to the host
        !$omp target update from(field_f%arr)

        ! Compute the maximum absolute difference between the array
        diff_local = MAXVAL(ABS(arr - field_f%arr))

        CALL MPI_Allreduce(diff_local, diff_global, 1, mglet_mpi_real, &
            MPI_MAX, MPI_COMM_WORLD, ierr)

        IF (diff_global > 3.0 * eps) THEN
            IF (myid == 0) THEN
                WRITE(*, *) "field name: ", TRIM(field_f%name), &
                    "  with maxdiff = ", diff_global
            END IF
            is_same = .FALSE.
        END IF
        IF (.NOT. is_same) CALL errr(__FILE__, __LINE__)

        DEALLOCATE(arr)

    END SUBROUTINE assert_same_field


    SUBROUTINE count_selftasks(nselfsend, nselfrecv)
        ! Subroutine arguments
        INTEGER(intk), INTENT(out) :: nselfsend, nselfrecv

        ! Local variables
        INTEGER(intk) :: i, iprocnbr

        nselfsend = 0
        DO i = 1, isend
            iprocnbr = sendconns(2, i)
            IF (iprocnbr == myid) THEN
                nselfsend = nselfsend + 1
            END IF
        END DO

        nselfrecv = 0
        DO i = 1, irecv
            iprocnbr = recvconns(1, i)
            IF (iprocnbr == myid) THEN
                nselfrecv = nselfrecv + 1
            END IF
        END DO

        IF (nselfsend /= nselfrecv) THEN
            CALL errr(__FILE__, __LINE__)
        END IF
    END SUBROUTINE count_selftasks


    SUBROUTINE finish_ftoc2()
        ! Subroutine arguments
        ! none...

        ! Local variables
        ! none...

        IF (.NOT. is_init) CALL errr(__FILE__, __LINE__)

        nsend = 0
        nrecv = 0
        DEALLOCATE(recvreqs)
        DEALLOCATE(sendreqs)
        DEALLOCATE(recvlist)
        DEALLOCATE(recvidxlist)

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
    END SUBROUTINE finish_ftoc2


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
END MODULE ftoc2_mod
