MODULE parent2_mod
    USE MPI_f08
    USE core_mod
    USE parent1_mod, ONLY: parent1
    USE parent_core_mod

    IMPLICIT NONE (type, external)
    PRIVATE

    ! Lists that hold the send and receive request arrays
    TYPE(MPI_Request), ALLOCATABLE :: sendreqs(:), recvreqs(:)
    INTEGER(int32), ALLOCATABLE :: sendlist(:), recvlist(:)
    INTEGER(intk), ALLOCATABLE :: recvidxlist(:, :)
    INTEGER(intk) :: nsend, nrecv

    ! Task extents
    INTEGER(intk), PARAMETER :: sendtasksize = 9
    INTEGER(intk), PARAMETER :: recvtasksize = 8
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

    ! Pointers to fields passed to parent2, point to dummy field if not present
    TYPE(field_t), POINTER :: f1 => NULL(), f2 => NULL(), f3 => NULL(), &
        f4 => NULL(), f5 => NULL(), f6 => NULL()

    LOGICAL :: is_init = .FALSE.

    PUBLIC :: parent2, init_parent2, finish_parent2
CONTAINS
    SUBROUTINE parent2(ilevel, v1, v2, v3, s1, s2, s3, normal)
        ! Subroutine arguments
        INTEGER(intk), INTENT(in) :: ilevel
        TYPE(field_t), OPTIONAL, TARGET, INTENT(inout) :: v1, v2, v3, s1, s2, s3
        LOGICAL, OPTIONAL, INTENT(in) :: normal

        ! Local variables
        TYPE(field_t), POINTER :: dummy
        INTEGER(intk) :: nvars
        INTEGER(intk) :: idx_args, idx_normal
        LOGICAL :: normal2

        IF (.NOT. is_init) CALL errr(__FILE__, __LINE__)

        ! Count number of variables to send/receive
        IF ((PRESENT(v1) .NEQV. PRESENT(v2)) .OR. &
                (PRESENT(v1) .NEQV. PRESENT(v3))) THEN
            CALL errr(__FILE__, __LINE__)
        END IF

        IF (ilevel < minlevel .OR. ilevel > maxlevel) THEN
            CALL errr(__FILE__, __LINE__)
        END IF

        normal2 = .FALSE.
        idx_normal = 0
        IF (PRESENT(normal)) THEN
            normal2 = normal
            IF (normal2) idx_normal = 1
        END IF

        CALL push_field(dummy, "PARENT2_DUMMY")

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

        ! If one vector component is present, all must be present
        IF (normal2 .AND. PRESENT(v1)) nvars = nvars - 2

        ! idx_level = ilevel in parent

        ! Looking up the workpackage in the workrecords array
        ASSOCIATE(wptr => workrecords(ilevel, idx_args, idx_normal))

        IF (is_recording) THEN
            ! During the recording phase, the workpackage is not yet initialized
            ! and needs to be created. Tasks arrays are offloaded for later
            ! efficient task execution on device
            CALL recording_pass(wptr, ilevel, normal2, nvars, v1, v2, v3, s1, &
                s2, s3)
        ELSE
            IF (wptr%is_init) THEN
                CALL recorded_parent(wptr)
            ELSE
                CALL jit_parent(ilevel, normal2, nvars, v1, v2, v3, s1, s2, s3)
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
    END SUBROUTINE parent2


    SUBROUTINE jit_parent(ilevel, normal, nvars, v1, v2, v3, s1, s2, s3)
        ! Subroutine arguments
        INTEGER(intk), INTENT(in) :: ilevel
        LOGICAL, INTENT(in) :: normal
        INTEGER(intk), INTENT(in) :: nvars
        TYPE(field_t), OPTIONAL, INTENT(inout) :: v1, v2, v3, s1, s2, s3

        ! Local variables
        INTEGER(intk) :: nsendtasks, nrecvtasks
        INTEGER(intk) :: nselftasks
        INTEGER(intk) :: nmpisendtasks

        ! Manufacturing the workpackage just in-time and execute
        CALL prepare_tasks_all(sendtasks, nsendtasks, selftasks, nselftasks, &
            mpisendtasks, nmpisendtasks, ilevel, normal, nvars, v1, v2, v3, &
            s1, s2, s3)

        !$omp target update to( &
        !$omp& sendtasks(1:sendtasksize, 1:nsendtasks+1), &
        !$omp& selftasks(1:selftasksize, 1:nselftasks+1)) nowait

        CALL recv_mpi_all(ilevel, nvars)

        !$omp taskwait

        CALL process_sendtasks(nsendtasks, sendtasks)
        CALL process_mpisend(nmpisendtasks, mpisendtasks)
        CALL process_selftasks(nselftasks, selftasks)

        CALL prepare_recvtasks_all(recvtasks, nrecvtasks, normal, v1, v2, v3, &
            s1, s2, s3)

        !$omp target update to(recvtasks(1:recvtasksize, 1:nrecvtasks+1))

        CALL process_recvtasks(nrecvtasks, recvtasks)
    END SUBROUTINE jit_parent


    SUBROUTINE recv_mpi_all(ilevel, nvars)
        ! Subroutine arguments
        INTEGER(intk), INTENT(in) :: ilevel
        INTEGER(intk), INTENT(in) :: nvars

        ! Local variables
        INTEGER(intk) :: i, iprocnbr, igrid, iface, ilevelgrid, facearea
        INTEGER(int32) :: recvcounter, messagelength

        ! Post all receive calls
        recvcounter = 0
        messagelength = 0
        nrecv = 0
        recvidxlist = 0

        DO i = 1, irecv
            iprocnbr = recvconns(2, i)
            igrid = recvconns(3, i)
            iface = recvconns(5, i)
            ilevelgrid = level(igrid)

            IF (iprocnbr == myid) CYCLE

            IF (ilevel == ilevelgrid) THEN
                facearea = face_area(igrid, iface)
                recvidxlist(1, i) = iprocnbr
                recvidxlist(2, i) = nvars*facearea
                recvidxlist(3, i) = recvcounter + messagelength
                messagelength = messagelength + nvars*facearea

                IF (recvcounter + messagelength > idim_mg_bufs) THEN
                    CALL errr(__FILE__, __LINE__)
                END IF
            END IF

            IF (messagelength > 0) THEN
                IF (i == irecv) THEN
                    CALL post_recv(iprocnbr, messagelength, recvcounter)
                ELSE IF (recvconns(2, i + 1) /= iprocnbr) THEN
                    CALL post_recv(iprocnbr, messagelength, recvcounter)
                END IF
            END IF
        END DO
    END SUBROUTINE recv_mpi_all


    SUBROUTINE post_recv(iprocnbr, messagelength, recvcounter)
        ! Identifier of receive connection
        INTEGER(int32), INTENT(in) :: iprocnbr
        INTEGER(int32), INTENT(inout) :: messagelength
        INTEGER(int32), INTENT(inout) :: recvcounter

        ! Local variables (for convenience)
        ! none...

        nrecv = nrecv + 1
        recvlist(nrecv) = iprocnbr

        CALL MPI_Irecv(device_recvbuf(recvcounter+1), messagelength, &
            mglet_mpi_real, iprocnbr, 1, MPI_COMM_WORLD, recvreqs(nrecv))

        recvcounter = recvcounter + messagelength
        messagelength = 0
    END SUBROUTINE post_recv


    SUBROUTINE recorded_parent(wptr)
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
    END SUBROUTINE recorded_parent


    SUBROUTINE recording_pass(wptr, ilevel, normal, nvars, v1, v2, v3, s1, s2, &
            s3)
        ! Subroutine arguments
        TYPE(work_t), INTENT(inout) :: wptr
        INTEGER(intk), INTENT(in) :: ilevel
        LOGICAL, INTENT(in) :: normal
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

        CALL prepare_mpirecvtasks(mpirecvtasks, nmpirecvtasks, ilevel, &
            normal, nvars)

        CALL prepare_tasks_all(sendtasks, nsendtasks, selftasks, nselftasks, &
            mpisendtasks, nmpisendtasks, ilevel, normal, nvars, v1, v2, v3, &
            s1, s2, s3)

        !$omp target update to( &
        !$omp& sendtasks(1:sendtasksize, 1:nsendtasks+1), &
        !$omp& selftasks(1:selftasksize, 1:nselftasks+1))

        CALL process_mpirecv(nmpirecvtasks, mpirecvtasks)
        CALL process_sendtasks(nsendtasks, sendtasks)
        CALL process_mpisend(nmpisendtasks, mpisendtasks)
        CALL process_selftasks(nselftasks, selftasks)
        CALL prepare_recvtasks_all(recvtasks, nrecvtasks, normal, v1, v2, v3, &
            s1, s2, s3)

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


    SUBROUTINE process_mpirecv(nmpirtasks, mpirtasks)
        ! Subroutine arguments
        INTEGER(intk), INTENT(in) :: nmpirtasks
        INTEGER(intk), INTENT(in) :: mpirtasks(mpitasksize, nmpirtasks+1)

        ! Local variables
        INTEGER(int32) :: itask, iprocnbr, messagelength, recvcounter

#ifdef _MGLET_PROFILE_ANNOTATIONS_
        CALL profile_range_push("process_mpirecv")
#endif

        DO itask = 1, nmpirtasks
            iprocnbr = INT(mpirtasks(1, itask), kind=int32)
            messagelength = INT(mpirtasks(2, itask), kind=int32)
            recvcounter = INT(mpirtasks(3, itask), kind=int32)

            CALL MPI_Irecv(device_recvbuf(recvcounter+1), messagelength, &
                mglet_mpi_real, iprocnbr, 1, MPI_COMM_WORLD, recvreqs(itask))
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


    SUBROUTINE process_sendtasks(nstasks, stasks)
        ! Subroutine arguments
        INTEGER(intk), INTENT(in) :: nstasks
        INTEGER(intk), INTENT(in) :: stasks(sendtasksize, nstasks+1)

        ! Local variables
        ! none...

        IF (nstasks == 0) RETURN

#ifdef _MGLET_PROFILE_ANNOTATIONS_
        CALL profile_range_push("process_sendtasks")
#endif

        CALL process_sendtasks_impl(sendbuf, nstasks, stasks, f1%arr, f2%arr, &
            f3%arr, f4%arr, f5%arr, f6%arr)

#ifdef _MGLET_PROFILE_ANNOTATIONS_
        CALL profile_range_pop()
#endif
    END SUBROUTINE process_sendtasks


    SUBROUTINE process_sendtasks_impl(sbuf, nstasks, stasks, a1, a2, a3, a4, &
        a5, a6)
        ! Subroutine arguments
        REAL(realk), INTENT(inout) :: sbuf(*)
        INTEGER(intk), INTENT(in) :: nstasks
        INTEGER(intk), INTENT(in) :: stasks(sendtasksize, nstasks+1)
        REAL(realk), INTENT(in) :: a1(*), a2(*), a3(*), a4(*), a5(*), a6(*)


        ! Local variables
        INTEGER(intk) :: itask, fieldid, icount, igrid, istart, istop, &
            jstart, jstop, kstart, kstop, ii, jj, kk, ip3

        !$omp target teams distribute private(itask, fieldid, icount, igrid, &
        !$omp& istart, istop, jstart, jstop, kstart, kstop, ii, jj, kk, ip3)
        DO itask = 1, nstasks
            fieldid = stasks(1, itask)
            icount = stasks(2, itask)
            igrid = stasks(3, itask)
            istart = stasks(4, itask)
            istop = stasks(5, itask)
            jstart = stasks(6, itask)
            jstop = stasks(7, itask)
            kstart = stasks(8, itask)
            kstop = stasks(9, itask)

            CALL get_ip3(ip3, igrid)
            CALL get_mgdims(kk, jj, ii, igrid)

            !$omp parallel
            SELECT CASE(fieldid)
            CASE (1)
                CALL arr_to_sendbuf(sbuf(icount), kk, jj, ii, a1(ip3), &
                    istart, istop, jstart, jstop, kstart, kstop)
            CASE (2)
                CALL arr_to_sendbuf(sbuf(icount), kk, jj, ii, a2(ip3), &
                    istart, istop, jstart, jstop, kstart, kstop)
            CASE (3)
                CALL arr_to_sendbuf(sbuf(icount), kk, jj, ii, a3(ip3), &
                    istart, istop, jstart, jstop, kstart, kstop)
            CASE (4)
                CALL arr_to_sendbuf(sbuf(icount), kk, jj, ii, a4(ip3), &
                    istart, istop, jstart, jstop, kstart, kstop)
            CASE (5)
                CALL arr_to_sendbuf(sbuf(icount), kk, jj, ii, a5(ip3), &
                    istart, istop, jstart, jstop, kstart, kstop)
            CASE (6)
                CALL arr_to_sendbuf(sbuf(icount), kk, jj, ii, a6(ip3), &
                    istart, istop, jstart, jstop, kstart, kstop)
#ifdef _MGLET_DEBUG_
            CASE DEFAULT
                CALL errr(__FILE__, __LINE__)
#endif
            END SELECT
            !$omp end parallel
        END DO
        !$omp end target teams distribute
    END SUBROUTINE process_sendtasks_impl


    SUBROUTINE arr_to_sendbuf(sbuf, kk, jj, ii, arr, istart, istop, &
            jstart, jstop, kstart, kstop)
        !$omp declare target
        ! Subroutine arguments
        REAL(realk), INTENT(inout) :: sbuf(*)
        INTEGER(intk), INTENT(in) :: kk, jj, ii
        REAL(realk), INTENT(in) :: arr(kk, jj, ii)
        INTEGER(intk), INTENT(in) :: istart, istop, jstart, jstop, kstart, kstop

        ! Local variables
        INTEGER(intk) :: i, j, k, idx, kkl, jjl

        kkl = kstop - kstart + 1
        jjl = jstop - jstart + 1

        !$omp do collapse(3) private(i, j, k, idx)
        DO i = istart, istop
            DO j = jstart, jstop
                DO k = kstart, kstop
                    idx = (k-kstart)+(j-jstart)*kkl+(i-istart)*jjl*kkl+1
                    sbuf(idx) = arr(k, j, i)
                END DO
            END DO
        END DO
        !$omp end do
    END SUBROUTINE arr_to_sendbuf


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

        CALL process_recvtasks_impl(nrtasks, rtasks, f1%buffers, &
            f2%buffers, f3%buffers, f4%buffers, f5%buffers, f6%buffers, &
            recvbuf)

#ifdef _MGLET_PROFILE_ANNOTATIONS_
        CALL profile_range_pop()
#endif
    END SUBROUTINE process_recvtasks


    SUBROUTINE process_recvtasks_impl(nrtasks, rtasks, a1, a2, a3, a4, &
            a5, a6, rbuf)
        ! Subroutine arguments
        INTEGER(intk), INTENT(in) :: nrtasks
        INTEGER(intk), INTENT(in) :: rtasks(recvtasksize, nrtasks+1)
        REAL(realk), INTENT(inout) :: a1(*), a2(*), a3(*), a4(*), a5(*), a6(*)
        REAL(realk), INTENT(in) :: rbuf(*)

        ! Local variables
        INTEGER(intk) :: itask, fieldid
        INTEGER(intk) :: jjc2d, ii2d, jj2d, ibb, stag1, stag2
        INTEGER(int32) :: icount

        !$omp target teams distribute private(itask, fieldid, icount, ibb, &
        !$omp& jjc2d, jj2d, ii2d, stag1, stag2)
        DO itask = 1, nrtasks
            fieldid = rtasks(1, itask)
            icount = INT(rtasks(2, itask), kind=int32)
            ibb = rtasks(3, itask)
            jjc2d = rtasks(4, itask)
            jj2d = rtasks(5, itask)
            ii2d = rtasks(6, itask)
            stag1 = rtasks(7, itask)
            stag2 = rtasks(8, itask)

            !$omp parallel
            SELECT CASE (fieldid)
            CASE (1)
                CALL recvbuf_to_buffers(a1(ibb), rbuf(icount), jjc2d, jj2d, &
                    ii2d, stag1, stag2)
            CASE (2)
                CALL recvbuf_to_buffers(a2(ibb), rbuf(icount), jjc2d, jj2d, &
                    ii2d, stag1, stag2)
            CASE (3)
                CALL recvbuf_to_buffers(a3(ibb), rbuf(icount), jjc2d, jj2d, &
                    ii2d, stag1, stag2)
            CASE (4)
                CALL recvbuf_to_buffers(a4(ibb), rbuf(icount), jjc2d, jj2d, &
                    ii2d, stag1, stag2)
            CASE (5)
                CALL recvbuf_to_buffers(a5(ibb), rbuf(icount), jjc2d, jj2d, &
                    ii2d, stag1, stag2)
            CASE (6)
                CALL recvbuf_to_buffers(a6(ibb), rbuf(icount), jjc2d, jj2d, &
                    ii2d, stag1, stag2)
#ifdef _MGLET_DEBUG_
            CASE DEFAULT
                CALL errr(__FILE__, __LINE__)
#endif
            END SELECT
            !$omp end parallel
        END DO
        !$omp end target teams distribute
    END SUBROUTINE process_recvtasks_impl


    SUBROUTINE recvbuf_to_buffers(buffers, rbuf, jjc2d, jj2d, ii2d, &
            stag1, stag2)
        !$omp declare target

        REAL(realk), INTENT(inout) :: buffers(*)
        REAL(realk), INTENT(in) :: rbuf(*)
        INTEGER(intk), INTENT(in) :: jjc2d, jj2d, ii2d
        INTEGER(intk), INTENT(in) :: stag1, stag2

        INTEGER(intk) :: j, i, jc, ic, idst
        REAL(realk) :: val_c, val_jm1, val_i, val_im1, val_out
        LOGICAL :: odd_j, odd_i

        !$omp do collapse(2) private(i, j, jc, ic, idst, val_c, &
        !$omp& val_jm1, val_i, val_im1, val_out, odd_j, odd_i)
        DO i = 1, ii2d
            DO j = 1, jj2d
                odd_j = MOD(j, 2_intk) == 1_intk
                odd_i = MOD(i, 2_intk) == 1_intk

                jc = 2 + (j-1)/2
                ic = 2 + (i-1)/2

                val_c = rbuf(1 + (jc-1) + (ic-1)*jjc2d + 1)
                val_jm1 = val_c
                IF (stag1 == 1 .AND. odd_j) THEN
                    val_jm1 = rbuf(1 + (jc-2) + (ic-1)*jjc2d + 1)
                END IF
                CALL prolong_face_first_direction(val_i, val_c, val_jm1, &
                    stag1, odd_j)

                val_im1 = val_i
                IF (stag2 == 1 .AND. odd_i) THEN
                    val_c = rbuf(1 + (jc-1) + (ic-2)*jjc2d + 1)
                    val_jm1 = val_c
                    IF (stag1 == 1 .AND. odd_j) THEN
                        val_jm1 = rbuf(1 + (jc-2) + (ic-2)*jjc2d + 1)
                    END IF
                    CALL prolong_face_first_direction(val_im1, val_c, val_jm1, &
                        stag1, odd_j)
                END IF

                CALL prolong_face_second_direction(val_out, val_i, val_im1, &
                    stag2, odd_i)

                idst = 1 + (j-1) + (i-1)*jj2d
                buffers(idst) = val_out
            END DO
        END DO
        !$omp end do
    END SUBROUTINE recvbuf_to_buffers


    SUBROUTINE prolong_face_first_direction(val, val_c, val_m1, stag1, odd)
        !$omp declare target
        ! Subroutine arguments
        REAL(realk), INTENT(out) :: val
        REAL(realk), INTENT(in) :: val_c, val_m1
        INTEGER(intk), INTENT(in) :: stag1
        LOGICAL, INTENT(in) :: odd

        ! Local variables
        ! none...

        SELECT CASE (stag1)
        CASE (0)
            val = val_c
        CASE (1)
            IF (odd) THEN
                val = 0.5_realk*(val_c + val_m1)
            ELSE
                val = val_c
            END IF
#ifdef _MGLET_DEBUG_
        CASE DEFAULT
            CALL errr(__FILE__, __LINE__)
#endif
        END SELECT
    END SUBROUTINE prolong_face_first_direction


    SUBROUTINE prolong_face_second_direction(val, val_c, val_m1, stag2, odd)
        !$omp declare target
        ! Subroutine arguments
        REAL(realk), INTENT(out) :: val
        REAL(realk), INTENT(in) :: val_c, val_m1
        INTEGER(intk), INTENT(in) :: stag2
        LOGICAL, INTENT(in) :: odd

        ! Local variables
        ! none...

        SELECT CASE (stag2)
        CASE (0)
            val = val_c
        CASE (1)
            IF (odd) THEN
                val = 0.5_realk*(val_c + val_m1)
            ELSE
                val = val_c
            END IF
#ifdef _MGLET_DEBUG_
        CASE DEFAULT
            CALL errr(__FILE__, __LINE__)
#endif
        END SELECT
    END SUBROUTINE prolong_face_second_direction


    SUBROUTINE process_selftasks(nstasks, stasks)
        ! Subroutine arguments
        INTEGER(intk), INTENT(in) :: nstasks
        INTEGER(intk), INTENT(in) :: stasks(selftasksize, nstasks+1)

        ! Local variables
        ! none...

        IF (nstasks == 0) RETURN

#ifdef _MGLET_PROFILE_ANNOTATIONS_
        CALL profile_range_push("process_selftasks")
#endif

        CALL process_selftasks_impl(nstasks, stasks, f1%arr, f2%arr, &
            f3%arr, f4%arr, f5%arr, f6%arr, f1%buffers, f2%buffers, &
            f3%buffers, f4%buffers, f5%buffers, f6%buffers)

#ifdef _MGLET_PROFILE_ANNOTATIONS_
        CALL profile_range_pop()
#endif
        END SUBROUTINE process_selftasks


        SUBROUTINE process_selftasks_impl(nstasks, stasks, a1, a2, a3, a4, a5, &
            a6, b1, b2, b3, b4, b5, b6)
        ! Subroutine arguments
        INTEGER(intk), INTENT(in) :: nstasks
        INTEGER(intk), INTENT(in) :: stasks(selftasksize, nstasks+1)
        REAL(realk), INTENT(in) :: a1(*), a2(*), a3(*), a4(*), a5(*), a6(*)
        REAL(realk), INTENT(inout) :: b1(*), b2(*), b3(*), b4(*), b5(*), b6(*)

        ! Local variables
        INTEGER(intk) :: itask, fieldid, igridc, ibb, istart, istop, jstart
        INTEGER(intk) :: jstop, kstart, kstop, jj2d, ii2d
        INTEGER(intk) :: stag1, stag2, kk, jj, ii, ip3

        !$omp target teams distribute private(itask, fieldid, igridc, ibb, &
        !$omp& istart, istop, jstart, jstop, kstart, kstop, jj2d, ii2d, &
        !$omp& stag1, stag2, kk, jj, ii, ip3)
        DO itask = 1, nstasks
            fieldid = stasks(1, itask)
            igridc = stasks(2, itask)
            ibb = stasks(3, itask)
            istart = stasks(4, itask)
            istop = stasks(5, itask)
            jstart = stasks(6, itask)
            jstop = stasks(7, itask)
            kstart = stasks(8, itask)
            kstop = stasks(9, itask)
            jj2d = stasks(12, itask)
            ii2d = stasks(13, itask)
            stag1 = stasks(14, itask)
            stag2 = stasks(15, itask)

            CALL get_ip3(ip3, igridc)
            CALL get_mgdims(kk, jj, ii, igridc)

            !$omp parallel
            SELECT CASE(fieldid)
            CASE (1)
                CALL arr_to_buffers(kk, jj, ii, a1(ip3), &
                    b1, ibb, istart, istop, jstart, jstop, kstart, kstop, &
                    jj2d, ii2d, stag1, stag2)
            CASE (2)
                CALL arr_to_buffers(kk, jj, ii, a2(ip3), &
                    b2, ibb, istart, istop, jstart, jstop, kstart, kstop, &
                    jj2d, ii2d, stag1, stag2)
            CASE (3)
                CALL arr_to_buffers(kk, jj, ii, a3(ip3), &
                    b3, ibb, istart, istop, jstart, jstop, kstart, kstop, &
                    jj2d, ii2d, stag1, stag2)
            CASE (4)
                CALL arr_to_buffers(kk, jj, ii, a4(ip3), &
                    b4, ibb, istart, istop, jstart, jstop, kstart, kstop, &
                    jj2d, ii2d, stag1, stag2)
            CASE (5)
                CALL arr_to_buffers(kk, jj, ii, a5(ip3), &
                    b5, ibb, istart, istop, jstart, jstop, kstart, kstop, &
                    jj2d, ii2d, stag1, stag2)
            CASE (6)
                CALL arr_to_buffers(kk, jj, ii, a6(ip3), &
                    b6, ibb, istart, istop, jstart, jstop, kstart, kstop, &
                    jj2d, ii2d, stag1, stag2)
#ifdef _MGLET_DEBUG_
            CASE DEFAULT
                CALL errr(__FILE__, __LINE__)
#endif
            END SELECT
            !$omp end parallel
        END DO
        !$omp end target teams distribute
    END SUBROUTINE process_selftasks_impl


    SUBROUTINE arr_to_buffers(kk, jj, ii, arr, buffers, ibb, &
            istart, istop, jstart, jstop, kstart, kstop, jj2d, ii2d, stag1, &
            stag2)
        !$omp declare target
        INTEGER(intk), INTENT(in) :: kk, jj, ii
        REAL(realk), INTENT(in) :: arr(kk, jj, ii)
        REAL(realk), INTENT(inout) :: buffers(*)
        INTEGER(intk), INTENT(in) :: ibb, istart, istop, jstart, jstop
        INTEGER(intk), INTENT(in) :: kstart, kstop, jj2d, ii2d, stag1, stag2

        INTEGER(intk) :: j, i, jc, ic, idst
        REAL(realk) :: val_c, val_jm1, val_i, val_im1, val_out
        LOGICAL :: odd_j, odd_i

        !$omp do collapse(2) private(i, j, jc, ic, idst, val_c, &
        !$omp& val_jm1, val_i, val_im1, val_out, odd_j, odd_i)
        DO i = 1, ii2d
            DO j = 1, jj2d
                ! Fine face indices are doubled relative to the parent face.
                ! Odd fine indices lie between two parent values if staggered.
                odd_j = MOD(j, 2_intk) == 1_intk
                odd_i = MOD(i, 2_intk) == 1_intk

                ! Map fine face indices to coarse parent-face indices. Two fine
                ! points share one parent point; parent index 1 is available for
                ! staggered interpolation from the parent-side neighbor.
                jc = 2_intk + (j-1_intk) / 2_intk
                ic = 2_intk + (i-1_intk) / 2_intk

                ! First interpolate along the face-local j direction.
                CALL get_parent_face_value(val_c, kk, jj, ii, arr, istart, &
                    istop, jstart, jstop, kstart, kstop, jc, ic)
                val_jm1 = val_c
                IF (stag1 == 1 .AND. odd_j) THEN
                    CALL get_parent_face_value(val_jm1, kk, jj, ii, arr, &
                        istart, istop, jstart, jstop, kstart, kstop, jc-1, ic)
                END IF
                CALL prolong_face_first_direction(val_i, val_c, val_jm1, &
                    stag1, odd_j)

                ! Then, if needed, interpolate the neighboring parent i line
                ! the same way before combining along the face-local i
                ! direction.
                val_im1 = val_i
                IF (stag2 == 1 .AND. odd_i) THEN
                    CALL get_parent_face_value(val_c, kk, jj, ii, arr, &
                        istart, istop, jstart, jstop, kstart, kstop, jc, ic-1)
                    val_jm1 = val_c
                    IF (stag1 == 1 .AND. odd_j) THEN
                        CALL get_parent_face_value(val_jm1, kk, jj, ii, arr, &
                            istart, istop, jstart, jstop, kstart, kstop, &
                            jc-1, ic-1)
                    END IF
                    CALL prolong_face_first_direction(val_im1, val_c, val_jm1, &
                        stag1, odd_j)
                END IF

                CALL prolong_face_second_direction(val_out, val_i, val_im1, &
                    stag2, odd_i)

                idst = ibb + (j-1) + (i-1)*jj2d
                buffers(idst) = val_out
            END DO
        END DO
        !$omp end do
    END SUBROUTINE arr_to_buffers


    SUBROUTINE get_parent_face_value(value, kk, jj, ii, arr, istart, istop, &
            jstart, jstop, kstart, kstop, jc, ic)
        !$omp declare target
        ! Subroutine arguments
        REAL(realk), INTENT(out) :: value
        INTEGER(intk), INTENT(in) :: kk, jj, ii
        REAL(realk), INTENT(in) :: arr(kk, jj, ii)
        INTEGER(intk), INTENT(in) :: istart, istop, jstart, jstop, kstart, kstop
        INTEGER(intk), INTENT(in) :: jc, ic

        ! Local variables
        ! none...

        ! The constant start/stop pair identifies which 3D direction is normal
        ! to the face; jc and ic address the two directions within that face.
        IF (istart == istop) THEN
            value = arr(kstart+jc-1, jstart+ic-1, istart)
        ELSE IF (jstart == jstop) THEN
            value = arr(kstart+jc-1, jstart, istart+ic-1)
        ELSE IF (kstart == kstop) THEN
            value = arr(kstart, jstart+jc-1, istart+ic-1)
#ifdef _MGLET_DEBUG_
        ELSE
            CALL errr(__FILE__, __LINE__)
#endif
        END IF
    END SUBROUTINE get_parent_face_value


    SUBROUTINE process_mpisend(nmpistasks, mpistasks)
        ! Subroutine arguments
        INTEGER(intk), INTENT(in) :: nmpistasks
        INTEGER(intk), INTENT(in) :: mpistasks(mpitasksize, nmpistasks+1)

        ! Local variables
        INTEGER(int32) :: itask, iprocnbr, messagelength, sendcounter

#ifdef _MGLET_PROFILE_ANNOTATIONS_
        CALL profile_range_push("process_mpisend")
#endif

        DO itask = 1, nmpistasks
            iprocnbr  = INT(mpistasks(1, itask), kind=int32)
            messagelength = INT(mpistasks(2, itask), kind=int32)
            sendcounter = INT(mpistasks(3, itask), kind=int32)

            CALL MPI_Isend(device_sendbuf(sendcounter+1), messagelength, &
                mglet_mpi_real, iprocnbr, 1, MPI_COMM_WORLD, sendreqs(itask))
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


    SUBROUTINE prepare_recvtasks_all(rtasks, nrtasks, normal, v1, v2, v3, &
            s1, s2, s3)
        ! Subroutine arguments
        INTEGER(intk), INTENT(inout) :: rtasks(recvtasksize, maxrecvtasks)
        INTEGER(intk), INTENT(out) :: nrtasks
        LOGICAL, INTENT(in) :: normal
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
                CALL MPI_Get_count(recvstatus, mglet_mpi_real, &
                    recvmessagelen)

                unpacklen = 0
                DO i = 1, irecv
                    IF (recvidxlist(1, i) == recvlist(idx) &
                            .AND. recvidxlist(2, i) > 0) THEN
                        CALL prepare_recvtask(rtasks, irtask, i, normal, &
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


    SUBROUTINE prepare_recvtask(rtask, irtask, recvid, normal, v1, v2, v3, s1, &
            s2, s3)
        ! Subroutine arguments
        INTEGER(intk), INTENT(inout) :: rtask(recvtasksize, maxrecvtasks)
        INTEGER(intk), INTENT(inout) :: irtask
        INTEGER(intk), INTENT(in) :: recvid
        LOGICAL, INTENT(in) :: normal
        TYPE(field_t), OPTIONAL, INTENT(inout) :: v1, v2, v3, s1, s2, s3

        ! Local variables
        INTEGER(intk) :: iface, igrid, facearea, ibb
        INTEGER(intk) :: kk, jj, ii, kkc, jjc, iic, jj2d, ii2d, jjc2d, iic2d
        INTEGER(intk) :: ustag1, ustag2, vstag1, vstag2, wstag1, wstag2
        INTEGER(int32) :: offset, icount
        LOGICAL :: exU, exV, exW


        igrid = recvconns(3, recvid)
        iface = recvconns(5, recvid)
        facearea = face_area(igrid, iface)

        offset = recvidxlist(3, recvid)
        icount = offset

        CALL get_mgdims(kk, jj, ii, igrid)
        CALL idx2d(kk, jj, ii, iface, kkc, jjc, iic, jj2d, ii2d, jjc2d, iic2d)
        CALL stag(iface, ustag1, ustag2, vstag1, vstag2, wstag1, wstag2)
        CALL get_ipbb(ibb, iface, igrid)

        exU = (normal .AND. iface < 3) .OR. (.NOT. normal)
        exV = (normal .AND. (iface > 2 .AND. iface < 5)) .OR. (.NOT. normal)
        exW = (normal .AND. iface > 4) .OR. (.NOT. normal)

        IF (PRESENT(v1) .AND. exU) THEN
            irtask = irtask + 1
            CALL add_single_recvtask(rtask(:, irtask), 1, icount, ibb, jjc2d, &
                iic2d, jj2d, ii2d, ustag1, ustag2, facearea)
        END IF
        IF (PRESENT(v2) .AND. exV) THEN
            irtask = irtask + 1
            CALL add_single_recvtask(rtask(:, irtask), 2, icount, ibb, jjc2d, &
                iic2d, jj2d, ii2d, vstag1, vstag2, facearea)
        END IF
        IF (PRESENT(v3) .AND. exW) THEN
            irtask = irtask + 1
            CALL add_single_recvtask(rtask(:, irtask), 3, icount, ibb, jjc2d, &
                iic2d, jj2d, ii2d, wstag1, wstag2, facearea)
        END IF
        IF (PRESENT(s1)) THEN
            irtask = irtask + 1
            CALL add_single_recvtask(rtask(:, irtask), 4, icount, ibb, jjc2d, &
                iic2d, jj2d, ii2d, 0, 0, facearea)
        END IF
        IF (PRESENT(s2)) THEN
            irtask = irtask + 1
            CALL add_single_recvtask(rtask(:, irtask), 5, icount, ibb, jjc2d, &
                iic2d, jj2d, ii2d, 0, 0, facearea)
        END IF
        IF (PRESENT(s3)) THEN
            irtask = irtask + 1
            CALL add_single_recvtask(rtask(:, irtask), 6, icount, ibb, jjc2d, &
                iic2d, jj2d, ii2d, 0, 0, facearea)
        END IF

        ! Check that message length is calculated correctly
        IF ((icount - offset) /= recvidxlist(2, recvid)) THEN
            WRITE(*, *) "icount:", icount, &
                "recvidxlist(2, recvid):", recvidxlist(2, recvid)
            CALL errr(__FILE__, __LINE__)
        END IF
    END SUBROUTINE prepare_recvtask


    PURE SUBROUTINE add_single_recvtask(task, fieldid, icount, ibb, jjc2d, &
            iic2d, jj2d, ii2d, stag1, stag2, unpacklen)
        ! Subroutine arguments
        INTEGER(intk), INTENT(out) :: task(recvtasksize)
        INTEGER(intk), INTENT(in) :: fieldid
        INTEGER(intk), INTENT(inout) :: icount
        INTEGER(intk), INTENT(in) :: ibb, jjc2d, iic2d, jj2d, ii2d
        INTEGER(intk), INTENT(in) :: stag1, stag2
        INTEGER(intk), INTENT(in) :: unpacklen

        ! Local variables
        ! none...

        task(1) = fieldid
        task(2) = icount
        task(3) = ibb
        task(4) = jjc2d
        task(5) = jj2d
        task(6) = ii2d
        task(7) = stag1
        task(8) = stag2

        icount = icount + unpacklen
    END SUBROUTINE add_single_recvtask


    PURE SUBROUTINE add_single_selftask(task, fieldid, igridc, ibb, istart, &
            istop, jstart, jstop, kstart, kstop, jjc2d, iic2d, jj2d, ii2d, &
            stag1, stag2)
        ! Subroutine arguments
        INTEGER(intk), INTENT(out) :: task(selftasksize)
        INTEGER(intk), INTENT(in) :: fieldid, igridc, ibb
        INTEGER(intk), INTENT(in) :: istart, istop, jstart, jstop, kstart, kstop
        INTEGER(intk), INTENT(in) :: jjc2d, iic2d, jj2d, ii2d, stag1, stag2

        task(1) = fieldid
        task(2) = igridc
        task(3) = ibb
        task(4) = istart
        task(5) = istop
        task(6) = jstart
        task(7) = jstop
        task(8) = kstart
        task(9) = kstop
        task(10) = jjc2d
        task(11) = iic2d
        task(12) = jj2d
        task(13) = ii2d
        task(14) = stag1
        task(15) = stag2
    END SUBROUTINE add_single_selftask


    SUBROUTINE prepare_mpirecvtasks(mpirtasks, nmpirtasks, ilevel, normal, &
            nvars)
        ! Subroutine arguments
        INTEGER(intk), INTENT(inout) :: mpirtasks(mpitasksize, maxmpirecvtasks)
        INTEGER(intk), INTENT(out) :: nmpirtasks
        INTEGER(intk), INTENT(in) :: ilevel
        LOGICAL, INTENT(in) :: normal
        INTEGER(intk), INTENT(in) :: nvars

        ! Local variables
        INTEGER(intk) :: i, iprocnbr, igrid, iface, ilevelgrid, facearea
        INTEGER(intk) :: impirtask
        INTEGER(int32) :: recvcounter, messagelength

        ! Post all receive calls
        recvcounter = 0
        messagelength = 0
        nrecv = 0
        recvidxlist = 0

        impirtask = 0

        DO i = 1, irecv
            iprocnbr = recvconns(2, i)
            igrid = recvconns(3, i)
            iface = recvconns(5, i)
            ilevelgrid = level(igrid)

            IF (iprocnbr == myid) CYCLE

            IF (ilevel == ilevelgrid) THEN
                facearea = face_area(igrid, iface)
                recvidxlist(1, i) = iprocnbr
                recvidxlist(2, i) = nvars*facearea
                recvidxlist(3, i) = recvcounter + messagelength
                messagelength = messagelength + nvars*facearea

                IF (recvcounter + messagelength > idim_mg_bufs) THEN
                    CALL errr(__FILE__, __LINE__)
                END IF
            END IF

            IF (messagelength > 0) THEN
                IF (i == irecv) THEN
                    CALL add_mpi_task(mpirtasks, impirtask, iprocnbr, &
                        messagelength, recvcounter, 'recv')
                ELSE IF (recvconns(2, i + 1) /= iprocnbr) THEN
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
            nmpisendtasks, ilevel, normal, nvars, v1, v2, v3, s1, s2, s3)
        ! Subroutine arguments
        INTEGER(intk), INTENT(inout) :: stasks(sendtasksize, maxsendtasks)
        INTEGER(intk), INTENT(out) :: nstasks
        INTEGER(intk), INTENT(inout) :: etasks(selftasksize, maxselftasks)
        INTEGER(intk), INTENT(out) :: netasks
        INTEGER(intk), INTENT(inout) :: mpistasks(mpitasksize, maxmpisendtasks)
        INTEGER(intk), INTENT(out) :: nmpisendtasks

        INTEGER(intk), INTENT(in) :: ilevel
        LOGICAL, INTENT(in) :: normal
        INTEGER(intk), INTENT(in) :: nvars
        TYPE(field_t), OPTIONAL, INTENT(inout) :: v1, v2, v3, s1, s2, s3

        ! Local variables
        INTEGER(intk) :: i, iprocnbr, igrid, iface, ilevelgrid, facearea
        INTEGER(intk) :: istask, ietask, impitasks
        INTEGER(int32) :: sendcounter, messagelength

#ifdef _MGLET_PROFILE_ANNOTATIONS_
        CALL profile_range_push("prepare_tasks_all")
#endif

        ! Pack all buffers and send data
        sendcounter = 0
        messagelength = 0
        nsend = 0

        ! Initializing the task counters to zero
        istask = 0
        ietask = 0
        impitasks = 0

        DO i = 1, isend
            iprocnbr = sendconns(1, i)
            igrid = sendconns(3, i)
            iface = sendconns(5, i)
            ilevelgrid = level(igrid)

            IF (ilevel == ilevelgrid) THEN
                IF (iprocnbr == myid) THEN
                    CALL prepare_selftasks(etasks, ietask, i, normal, &
                        v1, v2, v3, s1, s2, s3)
                    CYCLE
                ELSE
                    CALL prepare_sendtask(stasks, istask, i, messagelength, &
                        sendcounter, normal, nvars, v1, v2, v3, s1, s2, s3)
                    facearea = face_area(igrid, iface)

                    messagelength = messagelength + nvars*facearea
                END IF
            END IF

            IF (messagelength > 0) THEN
                IF (i == isend) THEN
                    CALL add_mpi_task(mpistasks, impitasks, iprocnbr, &
                        messagelength, sendcounter, 'send')
                ELSE IF (sendconns(1, i + 1) /= iprocnbr) THEN
                    CALL add_mpi_task(mpistasks, impitasks, iprocnbr, &
                        messagelength, sendcounter, 'send')
                END IF
            END IF
        END DO

        ! Set the output task counters
        nstasks = istask
        netasks = ietask
        nmpisendtasks = impitasks

        ! Add a harmful dummy task at (ntasks+1) for checking
        stasks(:, nstasks+1) = -1
        etasks(:, netasks+1) = -1
        mpistasks(:, nmpisendtasks+1) = -1

#ifdef _MGLET_PROFILE_ANNOTATIONS_
        CALL profile_range_pop()
#endif
    END SUBROUTINE prepare_tasks_all


    SUBROUTINE prepare_selftasks(etasks, ietask, sendid, normal, v1, v2, v3, &
            s1, s2, s3)
        ! Subroutine arguments
        INTEGER(intk), INTENT(inout) :: etasks(selftasksize, maxselftasks)
        INTEGER(intk), INTENT(inout) :: ietask
        INTEGER(intk), INTENT(in) :: sendid
        LOGICAL, INTENT(in) :: normal
        TYPE(field_t), OPTIONAL, INTENT(inout) :: v1, v2, v3, s1, s2, s3

        ! Local variables
        INTEGER(intk) :: igrid, igridc, iface, ibb
        INTEGER(intk) :: kk, jj, ii, kkc, jjc, iic, jj2d, ii2d, jjc2d, iic2d
        INTEGER(intk) :: istart, istop, jstart, jstop, kstart, kstop
        INTEGER(intk) :: ustag1, ustag2, vstag1, vstag2, wstag1, wstag2
        LOGICAL :: exU, exV, exW

        igrid = sendconns(3, sendid)
        igridc = sendconns(4, sendid)
        iface = sendconns(5, sendid)

        CALL get_mgdims(kk, jj, ii, igrid)
        CALL idx2d(kk, jj, ii, iface, kkc, jjc, iic, jj2d, ii2d, jjc2d, iic2d)
        CALL stag(iface, ustag1, ustag2, vstag1, vstag2, wstag1, wstag2)
        CALL get_ipbb(ibb, iface, igrid)

        exU = (normal .AND. iface < 3) .OR. (.NOT. normal)
        exV = (normal .AND. (iface > 2 .AND. iface < 5)) .OR. (.NOT. normal)
        exW = (normal .AND. iface > 4) .OR. (.NOT. normal)

        IF (PRESENT(v1) .AND. exU) THEN
            ietask = ietask + 1
            CALL get_senddims(istart, istop, jstart, jstop, kstart, kstop, &
                igrid, iface, v1)
            CALL add_single_selftask(etasks(:, ietask), 1, igridc, ibb, &
                istart, istop, jstart, jstop, kstart, kstop, jjc2d, iic2d, &
                jj2d, ii2d, ustag1, ustag2)
        END IF
        IF (PRESENT(v2) .AND. exV) THEN
            ietask = ietask + 1
            CALL get_senddims(istart, istop, jstart, jstop, kstart, kstop, &
                igrid, iface, v2)
            CALL add_single_selftask(etasks(:, ietask), 2, igridc, ibb, &
                istart, istop, jstart, jstop, kstart, kstop, jjc2d, iic2d, &
                jj2d, ii2d, vstag1, vstag2)
        END IF
        IF (PRESENT(v3) .AND. exW) THEN
            ietask = ietask + 1
            CALL get_senddims(istart, istop, jstart, jstop, kstart, kstop, &
                igrid, iface, v3)
            CALL add_single_selftask(etasks(:, ietask), 3, igridc, ibb, &
                istart, istop, jstart, jstop, kstart, kstop, jjc2d, iic2d, &
                jj2d, ii2d, wstag1, wstag2)
        END IF
        IF (PRESENT(s1)) THEN
            ietask = ietask + 1
            CALL get_senddims(istart, istop, jstart, jstop, kstart, kstop, &
                igrid, iface, s1)
            CALL add_single_selftask(etasks(:, ietask), 4, igridc, ibb, &
                istart, istop, jstart, jstop, kstart, kstop, jjc2d, iic2d, &
                jj2d, ii2d, 0, 0)
        END IF
        IF (PRESENT(s2)) THEN
            ietask = ietask + 1
            CALL get_senddims(istart, istop, jstart, jstop, kstart, kstop, &
                igrid, iface, s2)
            CALL add_single_selftask(etasks(:, ietask), 5, igridc, ibb, &
                istart, istop, jstart, jstop, kstart, kstop, jjc2d, iic2d, &
                jj2d, ii2d, 0, 0)
        END IF
        IF (PRESENT(s3)) THEN
            ietask = ietask + 1
            CALL get_senddims(istart, istop, jstart, jstop, kstart, kstop, &
                igrid, iface, s3)
            CALL add_single_selftask(etasks(:, ietask), 6, igridc, ibb, &
                istart, istop, jstart, jstop, kstart, kstop, jjc2d, iic2d, &
                jj2d, ii2d, 0, 0)
        END IF
    END SUBROUTINE prepare_selftasks


    SUBROUTINE prepare_sendtask(stasks, istask, sendid, messagelength, &
        sendcounter, normal, nvars, v1, v2, v3, s1, s2, s3)
        ! Subroutine arguments
        INTEGER(intk), INTENT(inout) :: stasks(sendtasksize, maxsendtasks)
        INTEGER(intk), INTENT(inout) :: istask
        INTEGER(intk), INTENT(in) :: sendid
        INTEGER(int32), INTENT(in) :: messagelength
        INTEGER(int32), INTENT(inout) :: sendcounter
        LOGICAL, INTENT(in) :: normal
        INTEGER(intk), INTENT(in) :: nvars
        TYPE(field_t), OPTIONAL, INTENT(inout) :: v1, v2, v3, s1, s2, s3


        ! Local variables
        INTEGER(intk) :: igrid, igridc, iface
        INTEGER(intk) :: istart, istop, jstart, jstop, kstart, kstop
        INTEGER(int32) :: thismessagelength, facearea, offset, icount
        LOGICAL :: exU, exV, exW

        ! Set variables from send table - *fine* grid and face
        igrid = sendconns(3, sendid)
        igridc = sendconns(4, sendid)
        iface = sendconns(5, sendid)

        ! Which vectors to exchange
        exU = (normal .AND. iface < 3) .OR. (.NOT. normal)
        exV = (normal .AND. (iface > 2 .AND. iface < 5)) .OR. (.NOT. normal)
        exW = (normal .AND. iface > 4) .OR. (.NOT. normal)

        ! Face area
        facearea = face_area(igrid, iface)
        thismessagelength = nvars*facearea

        ! Check that buffer does not overflow
        IF (sendcounter + messagelength + thismessagelength > idim_mg_bufs) THEN
            CALL errr(__FILE__, __LINE__)
        END IF

        ! Reset message size counter
        offset = sendcounter + messagelength + 1
        icount = offset

        ! Fill buffers
        IF (PRESENT(v1) .AND. exU) THEN
            istask = istask + 1
            CALL get_senddims(istart, istop, jstart, jstop, kstart, kstop, &
                igrid, iface, v1)
            CALL add_single_sendtask(stasks(:, istask), 1, icount, igridc, &
                istart, istop, jstart, jstop, kstart, kstop)
        END IF

        IF (PRESENT(v2) .AND. exV) THEN
            istask = istask + 1
            CALL get_senddims(istart, istop, jstart, jstop, kstart, kstop, &
                igrid, iface, v2)
            CALL add_single_sendtask(stasks(:, istask), 2, icount, igridc, &
                istart, istop, jstart, jstop, kstart, kstop)
        END IF

        IF (PRESENT(v3) .AND. exW) THEN
            istask = istask + 1
            CALL get_senddims(istart, istop, jstart, jstop, kstart, kstop, &
                igrid, iface, v3)
            CALL add_single_sendtask(stasks(:, istask), 3, icount, igridc, &
                istart, istop, jstart, jstop, kstart, kstop)
        END IF

        IF (PRESENT(s1)) THEN
            istask = istask + 1
            CALL get_senddims(istart, istop, jstart, jstop, kstart, kstop, &
                igrid, iface, s1)
            CALL add_single_sendtask(stasks(:, istask), 4, icount, igridc, &
                istart, istop, jstart, jstop, kstart, kstop)
        END IF

        IF (PRESENT(s2)) THEN
            istask = istask + 1
            CALL get_senddims(istart, istop, jstart, jstop, kstart, kstop, &
                igrid, iface, s2)
            CALL add_single_sendtask(stasks(:, istask), 5, icount, igridc, &
                istart, istop, jstart, jstop, kstart, kstop)
        END IF

        IF (PRESENT(s3)) THEN
            istask = istask + 1
            CALL get_senddims(istart, istop, jstart, jstop, kstart, kstop, &
                igrid, iface, s3)
            CALL add_single_sendtask(stasks(:, istask), 6, icount, igridc, &
                istart, istop, jstart, jstop, kstart, kstop)
        END IF

        IF (thismessagelength /= (icount - offset)) THEN
            WRITE(*, *) "thismessagelength:", thismessagelength, &
                "icount:", icount
            CALL errr(__FILE__, __LINE__)
        END IF
    END SUBROUTINE prepare_sendtask


    SUBROUTINE get_senddims(istart, istop, jstart, jstop, kstart, kstop, &
            igrid, iface, field)
        ! Subroutine arguments
        INTEGER(intk), INTENT(out) :: istart, istop, jstart, jstop, kstart, &
            kstop
        INTEGER(intk), INTENT(in) :: igrid, iface
        TYPE(field_t), INTENT(in) :: field

        ! Local variables
        INTEGER(intk) :: icomp

        icomp = 0
        SELECT CASE(iface)
        CASE (1, 2)
            IF (field%istag == 1) icomp = 1
        CASE (3, 4)
            IF (field%jstag == 1) icomp = 2
        CASE (5, 6)
            IF (field%kstag == 1) icomp = 3
        END SELECT

        CALL start_and_stop(igrid, iface, icomp, istart, istop, jstart, jstop, &
            kstart, kstop)
    END SUBROUTINE get_senddims


    PURE SUBROUTINE add_single_sendtask(task, fieldid, icount, igrid, &
            istart, istop, jstart, jstop, kstart, kstop)
        ! Subroutine arguments
        INTEGER(intk), INTENT(out) :: task(sendtasksize)
        INTEGER(intk), INTENT(in) :: fieldid
        INTEGER(int32), INTENT(inout) :: icount
        INTEGER(intk), INTENT(in) :: igrid
        INTEGER(intk), INTENT(in) :: istart, istop, jstart, jstop, kstart, kstop

        ! Local variables
        INTEGER(intk) :: tasksize

        task(1) = fieldid
        task(2) = INT(icount, kind=intk)
        task(3) = igrid
        task(4) = istart
        task(5) = istop
        task(6) = jstart
        task(7) = jstop
        task(8) = kstart
        task(9) = kstop

        tasksize = (istop-istart+1)*(jstop-jstart+1)*(kstop-kstart+1)
        icount = icount + tasksize
    END SUBROUTINE add_single_sendtask


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


    SUBROUTINE init_parent2()
        ! Subroutine arguments
        ! none...

        ! Local variables
        INTEGER(intk) :: nselfsend, nselfrecv

        IF (is_init) CALL errr(__FILE__, __LINE__)

        ! Check if parent_core_mod provides necessary infrastructure
        IF (.NOT. is_parent_core_init) CALL errr(__FILE__, __LINE__)

        ALLOCATE(recvidxlist(3, irecv))
        ALLOCATE(sendlist(isend))
        ALLOCATE(recvlist(irecv))
        ALLOCATE(sendreqs(isend))
        ALLOCATE(recvreqs(irecv))
        nrecv = 0
        nsend = 0

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
        ! dimension 3 = normal (0 or 1)
        ALLOCATE(workrecords(minlevel:maxlevel, 1:63, 0:1))

        is_init = .TRUE.

        ! Record relevant conn calls for later efficient reuse on device
        CALL run_recording_pass()
    END SUBROUTINE init_parent2


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
        CALL udummy%init_buffers()
        CALL vdummy%init_buffers()
        CALL wdummy%init_buffers()
        CALL pdummy%init_buffers()

        !$omp target enter data map(to: &
        !$omp&  udummy%arr, vdummy%arr, wdummy%arr, pdummy%arr, &
        !$omp&  udummy%buffers, vdummy%buffers, wdummy%buffers, pdummy%buffers)

        ! START -- This section defines the record variants of parent2 ---

        WRITE(*, *) minlevel, maxlevel
        DO ilevel = minlevel, maxlevel
            ! Inner pressuresolver iterations

            CALL init_dummy_fields_cpu(udummy, vdummy, wdummy, pdummy)
            !$omp target update to(udummy%arr, vdummy%arr, wdummy%arr, &
            !$omp& pdummy%arr, udummy%buffers, vdummy%buffers, &
            !$omp& wdummy%buffers, pdummy%buffers)

            CALL parent2(ilevel, s1=pdummy)
            CALL parent1(ilevel, s1=pdummy)

            CALL assert_same_field(pdummy, ilevel)
        END DO

        ! Outer pressuresolver iterations
        CALL init_dummy_fields_cpu(udummy, vdummy, wdummy, pdummy)
        !$omp target update to(udummy%arr, vdummy%arr, wdummy%arr, &
        !$omp& pdummy%arr, udummy%buffers, vdummy%buffers, &
        !$omp& wdummy%buffers, pdummy%buffers)

        CALL parent2(minlevel, v1=udummy, v2=vdummy, v3=wdummy, s1=pdummy)
        CALL parent1(minlevel, v1=udummy, v2=vdummy, v3=wdummy, s1=pdummy)

        CALL assert_same_field(udummy)
        CALL assert_same_field(vdummy)
        CALL assert_same_field(wdummy)
        CALL assert_same_field(pdummy)

        ! END -- This section defines the record variants of parent2 ---

        !$omp target exit data map(delete: &
        !$omp&  udummy%arr, vdummy%arr, wdummy%arr, pdummy%arr, &
        !$omp&  udummy%buffers, vdummy%buffers, wdummy%buffers, pdummy%buffers)

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

        DO i = 1, SIZE(u_f%buffers)
            u_f%buffers(i) = REAL(MOD(13_intk*i + 19_intk, 941_intk), realk) * &
                1.0e-3_realk
        END DO
        DO i = 1, SIZE(v_f%buffers)
            v_f%buffers(i) = REAL(MOD(29_intk*i + 17_intk, 937_intk), realk) * &
                1.0e-3_realk
        END DO
        DO i = 1, SIZE(w_f%buffers)
            w_f%buffers(i) = REAL(MOD(37_intk*i + 13_intk, 929_intk), realk) * &
                1.0e-3_realk
        END DO
        DO i = 1, SIZE(p_f%buffers)
            p_f%buffers(i) = REAL(MOD(43_intk*i + 9_intk, 919_intk), realk) * &
                1.0e-3_realk
        END DO
    END SUBROUTINE init_dummy_fields_cpu


    SUBROUTINE assert_same_field(field_f, ilevel)
        ! Subroutine arguments
        TYPE(field_t), INTENT(in) :: field_f
        INTEGER(intk), INTENT(in), OPTIONAL :: ilevel

        ! Local variables
        REAL(realk), ALLOCATABLE :: arr(:)
        REAL(realk), ALLOCATABLE :: buffers(:)
        REAL(realk) :: diff_local, diff_global
        INTEGER(intk) :: ierr
        LOGICAL :: do_buffers


        ALLOCATE(arr(SIZE(field_f%arr)))
        arr = field_f%arr

        IF (.NOT. ALLOCATED(field_f%buffers)) CALL errr(__FILE__, __LINE__)
        ALLOCATE(buffers(SIZE(field_f%buffers)))
        buffers = field_f%buffers

        !$omp target update from(field_f%arr)
        !$omp target update from(field_f%buffers)


        ! Looking at the buffers
        diff_local = MAXVAL(ABS(buffers - field_f%buffers))

        CALL MPI_Allreduce(diff_local, diff_global, 1, mglet_mpi_real, &
            MPI_MAX, MPI_COMM_WORLD, ierr)

        IF (diff_global > 3.0_realk*eps) THEN
            IF (myid == 0) THEN
                WRITE(*, *) "Buffers are not the same!"
                WRITE(*, *) "Max difference:", diff_global
                IF (PRESENT(ilevel)) THEN
                    WRITE(*, *) "Parent called from level:", ilevel
                END IF
            END IF
            ! CALL errr(__FILE__, __LINE__)
        END IF

        ! Looking at the arrays
        diff_local = MAXVAL(ABS(arr - field_f%arr))

        CALL MPI_Allreduce(diff_local, diff_global, 1, mglet_mpi_real, &
            MPI_MAX, MPI_COMM_WORLD, ierr)

        IF (diff_global > 3.0_realk*eps) THEN
            IF (myid == 0) THEN
                WRITE(*, *) "Arrays are not the same!"
                WRITE(*, *) "Max difference:", diff_global
                IF (PRESENT(ilevel)) THEN
                    WRITE(*, *) "Parent called from level:", ilevel
                END IF
            END IF
            ! CALL errr(__FILE__, __LINE__)
        END IF

        DEALLOCATE(buffers)
        DEALLOCATE(arr)

    END SUBROUTINE assert_same_field


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


    SUBROUTINE finish_parent2()
        ! Subroutine arguments
        ! none...

        ! Local variables
        ! none...

        IF (.NOT. is_init) CALL errr(__FILE__, __LINE__)

        nrecv = 0
        nsend = 0
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
    END SUBROUTINE finish_parent2


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
END MODULE parent2_mod
