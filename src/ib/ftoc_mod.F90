MODULE ftoc_mod
    USE core_mod
    USE ibcore_mod, ONLY: ib
    USE MPI_f08

    IMPLICIT NONE (type, external)
    PRIVATE

    ! The information in the first dimension is sorted as follows:
    !   Field 1: Rank of sending process
    !   Field 2: Rank of receiving process
    !   Field 3: ID of sending grid (fine grid)
    !   Field 4: ID of receiving grid (coarse grid)
    INTEGER(intk), ALLOCATABLE :: sendconns(:, :), recvconns(:, :)

    ! Lists that hold the send and receive request arrays
    TYPE(MPI_Request), ALLOCATABLE :: sendreqs(:), recvreqs(:)

    ! Lists that hold the messages that are ACTUALLY sendt and received
    INTEGER(intk) :: nsend, nrecv
    INTEGER(int32), ALLOCATABLE :: recvlist(:)
    INTEGER(intk), ALLOCATABLE :: recvidxlist(:, :)

    ! Number of send and receive connections
    INTEGER(intk) :: isend = 0, irecv = 0

    ! Variable to indicate if the required data structures have been created
    LOGICAL :: is_init = .FALSE.

    INTERFACE ftoc
        MODULE PROCEDURE :: ftoc_one, ftoc_multiple
    END INTERFACE ftoc

    ! contained functions
    PUBLIC :: ftoc, init_ftoc, finish_ftoc


CONTAINS
    SUBROUTINE ftoc_one(ilevel, ff, fc, flag)
        ! Subroutine arguments
        INTEGER(intk), INTENT(in) :: ilevel
        TYPE(field_t), INTENT(in) :: ff
        TYPE(field_t), INTENT(inout) :: fc
        CHARACTER(len=1), INTENT(in) :: flag

        ! Local variables
        ! none...

        CALL start_timer(220)

        IF (.NOT. is_init) CALL errr(__FILE__, __LINE__)

        CALL recv_all(ilevel, flag, fc)
        CALL send_all(ilevel, flag, ff)
        CALL process_bufs(flag, fc)

        nrecv = 0
        nsend = 0

        CALL stop_timer(220)
    END SUBROUTINE ftoc_one


    SUBROUTINE ftoc_multiple(ilevel, v1, v2, v3, s1, s2, s3)
        ! Subroutine arguments
        INTEGER(intk), INTENT(in) :: ilevel
        TYPE(field_t), INTENT(inout), OPTIONAL :: v1, v2, v3, s1, s2, s3

        ! Local variables
        CHARACTER(len=*), PARAMETER :: flag = '*'

        CALL start_timer(220)

        IF (.NOT. is_init) CALL errr(__FILE__, __LINE__)

        CALL recv_all(ilevel, flag, v1, v2, v3, s1, s2, s3)
        CALL send_all(ilevel, flag, v1, v2, v3, s1, s2, s3)
        CALL process_bufs(flag, v1, v2, v3, s1, s2, s3)

        nrecv = 0
        nsend = 0

        CALL stop_timer(220)
    END SUBROUTINE ftoc_multiple


    SUBROUTINE recv_all(ilevel, flag, v1, v2, v3, s1, s2, s3)
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
            ! Receiving grid
            igridf = recvconns(3, i)
            IF (ilevel == level(igridf)) THEN
                iprocnbr = recvconns(1, i)  ! The sender process (fine side)

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
    END SUBROUTINE recv_all


    ! Perform a single Recv
    SUBROUTINE post_recv(iprocnbr, messagelength, recvcounter)
        ! Identifier of receive connection
        INTEGER(int32), INTENT(in) :: iprocnbr
        INTEGER(int32), INTENT(inout) :: messagelength
        INTEGER(int32), INTENT(inout) :: recvcounter

        ! Local variables (for convenience)
        ! none...

        nrecv = nrecv + 1
        recvlist(nrecv) = iprocnbr
        CALL MPI_Irecv(recvbuf(recvcounter+1), messagelength, &
            mglet_mpi_real, iprocnbr, 1, MPI_COMM_WORLD, recvreqs(nrecv))

        recvcounter = recvcounter + messagelength
        messagelength = 0
    END SUBROUTINE post_recv


    SUBROUTINE count_ncells(ncells, flag, igrid, v1, v2, v3, s1, s2, s3)
        ! Subroutine arguments
        INTEGER(int32), INTENT(out) :: ncells
        CHARACTER(len=1), INTENT(in) :: flag
        INTEGER(intk), INTENT(in) :: igrid
        TYPE(field_t), INTENT(in), OPTIONAL :: v1, v2, v3, s1, s2, s3

        ! Local variables
        INTEGER(int32) :: n

        ncells = 0

        IF (flag == '*') THEN
            IF (PRESENT(v1)) THEN
                n = ib%restrict_op%message_length('U', igrid)
                ncells = ncells + n
            END IF
            IF (PRESENT(v2)) THEN
                n = ib%restrict_op%message_length('V', igrid)
                ncells = ncells + n
            END IF
            IF (PRESENT(v3)) THEN
                n = ib%restrict_op%message_length('W', igrid)
                ncells = ncells + n
            END IF
            IF (PRESENT(s1)) THEN
                n = ib%restrict_op%message_length('P', igrid)
                ncells = ncells + n
            END IF
            IF (PRESENT(s2)) THEN
                n = ib%restrict_op%message_length('P', igrid)
                ncells = ncells + n
            END IF
            IF (PRESENT(s3)) THEN
                n = ib%restrict_op%message_length('P', igrid)
                ncells = ncells + n
            END IF
        ELSE
            ncells = ib%restrict_op%message_length(flag, igrid)
        END IF
    END SUBROUTINE count_ncells


    ! Perform all send calls
    SUBROUTINE send_all(ilevel, flag, v1, v2, v3, s1, s2, s3)
        ! Subroutine arguments
        INTEGER(intk), INTENT(in) :: ilevel
        CHARACTER(len=1), INTENT(in) :: flag
        TYPE(field_t), OPTIONAL, INTENT(in) :: v1, v2, v3, s1, s2, s3

        ! Local variables
        INTEGER(intk) :: i, iprocnbr, igridf
        INTEGER(int32) :: sendcounter, messagelength

        ! Pack all buffers and send data
        sendcounter = 0
        messagelength = 0
        nsend = 0

        DO i = 1, isend
            igridf = sendconns(3, i)
            iprocnbr = sendconns(2, i)
            IF (ilevel == level(igridf)) THEN
                CALL write_buffer(i, messagelength, sendcounter, flag, &
                    v1, v2, v3, s1, s2, s3)
            END IF

            IF (messagelength > 0) THEN
                IF (i == isend) THEN
                    CALL post_send(iprocnbr, messagelength, sendcounter)
                ELSE IF (sendconns(2, i + 1) /= iprocnbr) THEN
                    CALL post_send(iprocnbr, messagelength, sendcounter)
                END IF
            END IF
        END DO
    END SUBROUTINE send_all


    ! Perform a single send call
    SUBROUTINE post_send(iprocnbr, messagelength, sendcounter)
        ! Subroutine arguments
        INTEGER(int32), INTENT(in) :: iprocnbr
        INTEGER(int32), INTENT(inout) :: messagelength
        INTEGER(int32), INTENT(inout) :: sendcounter

        ! Local variables (for convenience)
        ! none...

        nsend = nsend + 1
        CALL MPI_Isend(sendbuf(sendcounter + 1), messagelength, &
            mglet_mpi_real, iprocnbr, 1, MPI_COMM_WORLD, sendreqs(nsend))

        sendcounter = sendcounter + messagelength
        messagelength = 0
    END SUBROUTINE post_send


    SUBROUTINE write_buffer(sendid, messagelength, sendcounter, flag, &
            v1, v2, v3, s1, s2, s3)
        ! Subroutine arguments
        INTEGER(int32), INTENT(in) :: sendid
        INTEGER(int32), INTENT(inout) :: messagelength
        INTEGER(int32), INTENT(in) :: sendcounter
        CHARACTER(len=1), INTENT(in) :: flag
        TYPE(field_t), OPTIONAL, INTENT(in) :: &
            v1, v2, v3, s1, s2, s3

        ! Local variables
        INTEGER(intk) :: kk, jj, ii
        INTEGER(intk) :: igrid
        INTEGER(int32) :: thismessagelength, offset, ncells
        REAL(realk), POINTER, CONTIGUOUS :: field(:, :, :)
        CHARACTER(len=1) :: flag_u

        igrid = sendconns(3, sendid)
        CALL get_mgdims(kk, jj, ii, igrid)

        CALL count_ncells(thismessagelength, flag, igrid, &
            v1, v2, v3, s1, s2, s3)

        ! Check that buffer does not overflow
        IF (sendcounter + messagelength + thismessagelength > idim_mg_bufs) THEN
            CALL errr(__FILE__, __LINE__)
        END IF

        ! Reset message size counter
        offset = sendcounter + messagelength + 1

        ! Fill buffers
        IF (PRESENT(v1)) THEN
            IF (flag == '*') THEN
                flag_u = 'U'
            ELSE
                flag_u = flag
            END IF
            ncells = ib%restrict_op%message_length(flag_u, igrid)
            CALL v1%get_ptr(field, igrid)
            CALL ib%restrict_op%restrict(kk, jj, ii, field, &
                sendbuf(offset:offset+ncells-1), &
                flag_u, igrid)
            offset = offset + ncells
        END IF

        IF (PRESENT(v2)) THEN
            ncells = ib%restrict_op%message_length('V', igrid)
            CALL v2%get_ptr(field, igrid)
            CALL ib%restrict_op%restrict(kk, jj, ii, field, &
                sendbuf(offset:offset+ncells-1), &
                'V', igrid)
            offset = offset + ncells
        END IF

        IF (PRESENT(v3)) THEN
            ncells = ib%restrict_op%message_length('W', igrid)
            CALL v3%get_ptr(field, igrid)
            CALL ib%restrict_op%restrict(kk, jj, ii, field, &
                sendbuf(offset:offset+ncells-1), &
                'W', igrid)
            offset = offset + ncells
        END IF

        IF (PRESENT(s1)) THEN
            ncells = ib%restrict_op%message_length('P', igrid)
            CALL s1%get_ptr(field, igrid)
            CALL ib%restrict_op%restrict(kk, jj, ii, field, &
                sendbuf(offset:offset+ncells-1), &
                'P', igrid)
            offset = offset + ncells
        END IF

        IF (PRESENT(s2)) THEN
            ncells = ib%restrict_op%message_length('P', igrid)
            CALL s2%get_ptr(field, igrid)
            CALL ib%restrict_op%restrict(kk, jj, ii, field, &
                sendbuf(offset:offset+ncells-1), &
                'P', igrid)
            offset = offset + ncells
        END IF

        IF (PRESENT(s3)) THEN
            ncells = ib%restrict_op%message_length('P', igrid)
            CALL s3%get_ptr(field, igrid)
            CALL ib%restrict_op%restrict(kk, jj, ii, field, &
                sendbuf(offset:offset+ncells-1), &
                'P', igrid)
            offset = offset + ncells
        END IF

        IF (offset /= sendcounter + messagelength + thismessagelength + 1) THEN
            WRITE(*, *) "offset:", offset, &
                "expected:", sendcounter + messagelength + thismessagelength + 1
            CALL errr(__FILE__, __LINE__)
        END IF

        messagelength = messagelength + thismessagelength
    END SUBROUTINE write_buffer


    ! Process receive buffers as they arrive, wait for send
    ! buffers to be free
    SUBROUTINE process_bufs(flag, v1, v2, v3, s1, s2, s3)
        ! Subroutine arguments
        CHARACTER(len=1), INTENT(in) :: flag
        TYPE(field_t), OPTIONAL, INTENT(inout) :: &
            v1, v2, v3, s1, s2, s3

        INTEGER(int32) :: idx, i
        TYPE(MPI_Status) :: recvstatus
        INTEGER(int32) :: recvmessagelen
        INTEGER(int32) :: unpacklen

        DO WHILE (.TRUE.)
            IF (nrecv == 0) EXIT
            CALL MPI_Waitany(nrecv, recvreqs, idx, recvstatus)

            IF (idx /= MPI_UNDEFINED) THEN
                CALL MPI_Get_count(recvstatus, mglet_mpi_real, &
                    recvmessagelen)

                unpacklen = 0
                DO i = 1, irecv
                    IF (recvidxlist(1, i) == recvlist(idx) &
                            .AND. recvidxlist(2, i) > 0) THEN
                        CALL read_buffer(i, flag, v1, v2, v3, s1, s2, s3)
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
    END SUBROUTINE process_bufs


    ! Read Receive buffers
    !
    ! Write the contents of the receive buffers back in their
    ! matching fields
    SUBROUTINE read_buffer(recvid, flag, v1, v2, v3, s1, s2, s3)
        ! Subroutine arguments
        INTEGER(intk), INTENT(in) :: recvid
        CHARACTER(len=1), INTENT(in) :: flag
        TYPE(field_t), OPTIONAL, INTENT(inout) :: &
            v1, v2, v3, s1, s2, s3

        ! Grid dimensions and pointers
        INTEGER(intk) :: kk, jj, ii
        INTEGER(intk) :: igrid, igridc
        INTEGER(int32) :: offset, ncells
        REAL(realk), POINTER, CONTIGUOUS :: field(:, :, :)
        CHARACTER(len=1) :: flag_u

        igrid = recvconns(3, recvid)
        igridc = recvconns(4, recvid)
        CALL get_mgdims(kk, jj, ii, igridc)

        offset = recvidxlist(3, recvid) + 1

        IF (PRESENT(v1)) THEN
            IF (flag == '*') THEN
                flag_u = 'U'
            ELSE
                flag_u = flag
            END IF
            ncells = ib%restrict_op%message_length(flag_u, igrid)
            CALL v1%get_ptr(field, igridc)

            IF (flag_u == 'N') THEN
                CALL restrict_recieve_open_n(kk, jj, ii, field, &
                    recvbuf(offset:offset+ncells-1), igrid, flag_u)
            ELSE
                CALL restrict_recieve_open(kk, jj, ii, field, &
                    recvbuf(offset:offset+ncells-1), igrid, flag_u)
            END IF
            offset = offset + ncells
        END IF

        IF (PRESENT(v2)) THEN
            ncells = ib%restrict_op%message_length('V', igrid)
            CALL v2%get_ptr(field, igridc)
            CALL restrict_recieve_open(kk, jj, ii, field, &
                recvbuf(offset:offset+ncells-1), igrid, 'V')
            offset = offset + ncells
        END IF

        IF (PRESENT(v3)) THEN
            ncells = ib%restrict_op%message_length('W', igrid)
            CALL v3%get_ptr(field, igridc)
            CALL restrict_recieve_open(kk, jj, ii, field, &
                recvbuf(offset:offset+ncells-1), igrid, 'W')
            offset = offset + ncells
        END IF

        IF (PRESENT(s1)) THEN
            ncells = ib%restrict_op%message_length('P', igrid)
            CALL s1%get_ptr(field, igridc)
            CALL restrict_recieve_open(kk, jj, ii, field, &
                recvbuf(offset:offset+ncells-1), igrid, 'P')
            offset = offset + ncells
        END IF

        IF (PRESENT(s2)) THEN
            ncells = ib%restrict_op%message_length('P', igrid)
            CALL s2%get_ptr(field, igridc)
            CALL restrict_recieve_open(kk, jj, ii, field, &
                recvbuf(offset:offset+ncells-1), igrid, 'P')
            offset = offset + ncells
        END IF

        IF (PRESENT(s3)) THEN
            ncells = ib%restrict_op%message_length('P', igrid)
            CALL s3%get_ptr(field, igridc)
            CALL restrict_recieve_open(kk, jj, ii, field, &
                recvbuf(offset:offset+ncells-1), igrid, 'P')
            offset = offset + ncells
        END IF

        IF (offset - recvidxlist(3, recvid) - 1 /= recvidxlist(2, recvid)) THEN
            CALL errr(__FILE__, __LINE__)
        END IF
    END SUBROUTINE read_buffer


    SUBROUTINE init_ftoc()
        ! Local variables
        INTEGER(intk) :: i, igrid, iprocc, ipar, maxconns
        INTEGER(int32), ALLOCATABLE :: sendcounts(:), sdispls(:)
        INTEGER(int32), ALLOCATABLE :: recvcounts(:), rdispls(:)
        INTEGER(intk), PARAMETER :: ncols = 4

        CALL set_timer(220, "FTOC")

        IF (is_init) CALL errr(__FILE__, __LINE__)

        maxconns = nmygrids*8
        ALLOCATE(sendconns(ncols, maxconns))

        ALLOCATE(sendcounts(0:numprocs-1), SOURCE=0)
        ALLOCATE(sdispls(0:numprocs-1), SOURCE=0)
        ALLOCATE(recvcounts(0:numprocs-1), SOURCE=0)
        ALLOCATE(rdispls(0:numprocs-1), SOURCE=0)

        nsend = 0
        DO i = 1, nmygrids
            igrid = mygrids(i)
            ipar = iparent(igrid)
            IF (ipar == 0) CYCLE

            iprocc = idprocofgrd(ipar)

            nsend = nsend + 1
            IF (nsend > maxconns) THEN
                CALL errr(__FILE__, __LINE__)
            END IF

            sendconns(1, nsend) = myid
            sendconns(2, nsend) = iprocc
            sendconns(3, nsend) = igrid
            sendconns(4, nsend) = ipar

            sendcounts(iprocc) = sendcounts(iprocc) + ncols
        END DO
        isend = nsend

        ! Sort sendconns by process ID (col 2)
        CALL sort_conns(sendconns(:, 1:nsend), 2)

        ! Calculate sdispl offset
        DO i = 1, numprocs-1
            sdispls(i) = sdispls(i-1) + sendcounts(i-1)
        END DO

        ! First exchange NUMBER OF ELEMENTS TO RECEIVE, to be able to
        ! calculate rdispls array
        CALL MPI_Alltoall(sendcounts, 1, MPI_INTEGER, recvcounts, 1, &
            MPI_INTEGER, MPI_COMM_WORLD)

        ! Calculate rdispl offset
        DO i = 1, numprocs-1
            rdispls(i) = rdispls(i-1) + recvcounts(i-1)
        END DO

        irecv = (rdispls(numprocs-1) + recvcounts(numprocs-1))/ncols
        ALLOCATE(recvconns(ncols, irecv))
        recvconns = 0

        ! Exchange connection information
        CALL MPI_Alltoallv(sendconns, sendcounts, sdispls, MPI_INTEGER, &
            recvconns, recvcounts, rdispls, MPI_INTEGER, &
            MPI_COMM_WORLD)

        DEALLOCATE(rdispls)
        DEALLOCATE(recvcounts)
        DEALLOCATE(sdispls)
        DEALLOCATE(sendcounts)

        ALLOCATE(sendreqs(numprocs))
        ALLOCATE(recvreqs(numprocs))
        ALLOCATE(recvlist(irecv))
        ALLOCATE(recvidxlist(3, irecv))

        is_init = .TRUE.
        nsend = 0
        nrecv = 0
    END SUBROUTINE init_ftoc


    SUBROUTINE finish_ftoc()
        ! Function to deallocate arrays.
        ! After successful execution, the module varaible "is_init" is set
        ! to false.
        IF (.NOT. is_init) RETURN

        is_init = .FALSE.
        isend = 0
        irecv = 0
        DEALLOCATE(sendconns)
        DEALLOCATE(recvconns)
        DEALLOCATE(sendreqs)
        DEALLOCATE(recvreqs)
        DEALLOCATE(recvlist)
        DEALLOCATE(recvidxlist)
    END SUBROUTINE finish_ftoc


    SUBROUTINE restrict_recieve_open(kk, jj, ii, fc, buffer, igridf, flag)
        ! Subroutine arguments
        INTEGER(intk), INTENT(in) :: kk, jj, ii
        REAL(realk), INTENT(inout) :: fc(kk, jj, ii)
        REAL(realk), CONTIGUOUS, INTENT(in) :: buffer(:)
        INTEGER(intk), INTENT(in) :: igridf
        CHARACTER(len=1), INTENT(in) :: flag

        ! Local variables
        INTEGER(intk) :: istart, istop, jstart, jstop, kstart, kstop
        INTEGER(intk) :: ipos, jpos, kpos
        INTEGER(intk) :: i, j, k, ic, jc, kc, ic0, jc0, kc0, icount

        CALL ib%restrict_op%start_and_stop(istart, istop, jstart, jstop, &
            kstart, kstop, flag, igridf)

        ic0 = 0
        jc0 = 0
        kc0 = 0
        SELECT CASE (flag)
        CASE ("U")
            ic0 = 1
        CASE ("V")
            jc0 = 1
        CASE ("W")
            kc0 = 1
        END SELECT

        ipos = iposition(igridf)
        jpos = jposition(igridf)
        kpos = kposition(igridf)

        icount = 0
        DO i = istart, istop, 2
            ic = ipos + (i-3+ic0)/2 - ic0
            DO j = jstart, jstop, 2
                jc = jpos + (j-3+jc0)/2 - jc0
                DO k = kstart, kstop, 2
                    kc = kpos + (k-3+kc0)/2 - kc0

                    icount = icount + 1
                    fc(kc, jc, ic) = buffer(icount)
                END DO
            END DO
        END DO
    END SUBROUTINE restrict_recieve_open


    SUBROUTINE restrict_recieve_open_n(kk, jj, ii, fc, buffer, igridf, flag)
        ! Subroutine arguments
        INTEGER(intk), INTENT(in) :: kk, jj, ii
        REAL(realk), INTENT(inout) :: fc(kk, jj, ii)
        REAL(realk), CONTIGUOUS, INTENT(in) :: buffer(:)
        INTEGER(intk), INTENT(in) :: igridf
        CHARACTER(len=1), INTENT(in) :: flag

        ! Local variables
        INTEGER(intk) :: istart, istop, jstart, jstop, kstart, kstop
        INTEGER(intk) :: ipos, jpos, kpos
        INTEGER(intk) :: i, j, k, ic, jc, kc, icount

        CALL ib%restrict_op%start_and_stop(istart, istop, jstart, jstop, &
            kstart, kstop, flag, igridf)

        ipos = iposition(igridf)
        jpos = jposition(igridf)
        kpos = kposition(igridf)

        ! Warning: ic, jc, kc defined differently compared to
        ! restrict_receive_open
        ! TODO: check!
        icount = 0
        DO i = istart, istop, 2
            ic = ipos + (i-2)/2 - 1
            DO j = jstart, jstop, 2
                jc = jpos + (j-2)/2 - 1
                DO k = kstart, kstop, 2
                    kc = kpos + (k-2)/2 - 1

                    icount = icount + 1
                    IF (ABS(fc(kc, jc, ic)) < TINY(1.0_realk) .AND. &
                            buffer(icount) >= 1.0_realk) THEN
                        fc(kc, jc, ic) = buffer(icount)
                    END IF
                END DO
            END DO
        END DO
    END SUBROUTINE restrict_recieve_open_n

END MODULE ftoc_mod
