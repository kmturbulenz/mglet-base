MODULE ftoc_mod
    ! Module is responsible for the communication from fine to coarse grids.
    ! An injection method is used for the restriction.
    ! Under consideration of the shift of the momentum cells, the values of every second
    ! node on the finer level are used to define a value on the coarser grid.
    USE core_mod

    USE ibcore_mod, ONLY: ib
    USE MPI_f08

    IMPLICIT NONE (type, external)
    PRIVATE

    ! Variable to indicate if the required data structures have been created
    LOGICAL :: is_init = .FALSE.

    ! If .TRUE., a restriction is already in process and you cannot start
    ! another one
    LOGICAL :: in_progress = .FALSE.

    ! Maximum allowed number of childs per parent (i.e. maximum number of
    ! fine grids inside any coarse grid). Does not need to be exact, some grids
    ! can have many more.
    INTEGER(intk), PARAMETER :: max_childs = 8

    ! Lists that hold the send and receive request arrays
    TYPE(MPI_Request), ALLOCATABLE :: sendReqs(:), recvReqs(:)

    ! Actual number of messages that are sendt and received in one
    ! "round" of operations
    INTEGER(intk) :: nSend, nRecv

    ! Counters for send- and receive operations (for locations in
    ! send and receive buffers)
    INTEGER(int32) :: send_counter, recv_counter

    ! Message offsets, i.e. location in Send/Recv buffer.
    ! This list is filled each time 'ftoc' is called.
    ! Columns as follows:
    !
    !   Column 1: Offset in the send/recv buffer, first element has offset 0
    !   Column 2: Lengt of message, i.e. number of elements in the send/recv
    !             buffer
    INTEGER(int32), ALLOCATABLE :: recvOffsetLength(:,:)

    ! Flag to indicate behaviour (A, B, C, D, E, F, L, N, P, R, S, U, V, W)
    CHARACTER(len=1) :: flag

    ! List of *fine grids* to receive data *from*.
    INTEGER(intk), ALLOCATABLE :: recvGrids(:)

    ! pointers to fine and coarse field used in restriction
    REAL(realk), POINTER, CONTIGUOUS  :: ffg(:), fcg(:)

    ! contained functions
    PUBLIC :: ftoc, init_ftoc, finish_ftoc


CONTAINS
    SUBROUTINE ftoc(ilevel, ff_p, fc_p, flag_p)
        INTEGER(intk), INTENT(in) :: ilevel

        ! Fields to be interpolated
        REAL(realk), CONTIGUOUS, INTENT(inout) :: ff_p(:)

        ! field on the fine grid
        REAL(realk), CONTIGUOUS, INTENT(inout) :: fc_p(:)

        ! flag to specify the mode
        CHARACTER(len=1), INTENT(in) :: flag_p

        CALL start_timer(220)
        CALL ftoc_begin(ff_p, fc_p, flag_p, noflevel(ilevel), &
            igrdoflevel(1, ilevel))
        CALL ftoc_end()
        CALL stop_timer(220)
    END SUBROUTINE ftoc


    SUBROUTINE ftoc_begin(ff_p, fc_p, flag_p, ngrids, lofgrids)
        ! Function initiates the restriction.
        ! The module variable "in_progress" is set to true.
        ! Interpolate the results from a fine level to a coarse level
        ! This function initiate the process, and the ftoc_finish
        ! must be called afterwards to clean up.

        ! Fields to be interpolated
        REAL(realk), TARGET, CONTIGUOUS, INTENT(inout) :: ff_p(:)
        ! field on the fine grid
        REAL(realk), TARGET, CONTIGUOUS, INTENT(inout) :: fc_p(:)

        ! flag to specify the mode
        CHARACTER(len=1), INTENT(in) :: flag_p

        ! number of grids on the level
        INTEGER(intk), INTENT(in) :: ngrids

        ! grid ID's on level
        INTEGER(intk), INTENT(in) :: lofgrids(ngrids)

        ! Local variables
        INTEGER(intk) :: i, igrid, iprocc, iprocf, ipar

        CALL start_timer(221)

        IF (.NOT. is_init) THEN
            CALL errr(__FILE__, __LINE__)
        END IF

        IF (in_progress) THEN
            CALL errr(__FILE__, __LINE__)
        END IF

        ! reminder: The scope of the following variables is module wide
        in_progress = .TRUE.
        ! reminder: scope of the pointers is module-wide
        ffg => ff_p
        fcg => fc_p
        flag = flag_p

        nRecv = 0
        nSend = 0

        send_counter = 0
        recv_counter = 0

        ! Make all Recv-calls
        DO i = 1, ngrids
            igrid = lofgrids(i)
            ipar = iparent(igrid)
            IF (ipar /= 0) THEN
                iprocf = idprocofgrd(igrid)
                iprocc = idprocofgrd(ipar)

                IF (myid == iprocc) THEN
                    nRecv = nRecv + 1
                    CALL restrict_recv(igrid)
                    recvGrids(nRecv) = igrid
                END IF
            END IF
        END DO

        ! Make all Send-calls
        DO i = 1, ngrids
            igrid = lofgrids(i)
            ipar = iparent(igrid)
            IF (ipar /= 0 ) THEN
                iprocf = idprocofgrd(igrid)
                iprocc = idprocofgrd(ipar)

                IF (myid == iprocf) THEN
                    nSend = nSend + 1
                    CALL restrict_send(igrid, ipar)
                END IF
            END IF
        END DO

        CALL stop_timer(221)
    END SUBROUTINE ftoc_begin


    SUBROUTINE ftoc_end()
        ! Function finishes the restriction.
        ! The module variable "in_progress" is set to false.
        ! This means to wait for communication to finish and to clean up
        ! after restriction.

        INTEGER(int32) :: idx

        CALL start_timer(222)

        IF (.NOT. in_progress) THEN
            CALL errr(__FILE__, __LINE__)
        END IF

        IF (nRecv > 0) THEN
            DO WHILE (.TRUE.)
                CALL MPI_Waitany(nRecv, recvReqs, idx, MPI_STATUS_IGNORE)

                IF (idx /= MPI_UNDEFINED) THEN
                    CALL restrict_finish(recvGrids(idx))
                ELSE
                    EXIT
                END IF
            END DO
        END IF

        CALL MPI_Waitall(nSend, sendReqs, MPI_STATUSES_IGNORE)

        ! pointers are voided again
        NULLIFY(ffg)
        NULLIFY(fcg)

        flag = ' '
        in_progress = .FALSE.

        CALL stop_timer(222)
    END SUBROUTINE ftoc_end


    SUBROUTINE restrict_recv(igridf)
        INTEGER, INTENT(in) :: igridf

        INTEGER(intk) :: iprocf
        INTEGER(intk) :: message_length

        iprocf = idprocofgrd(igridf)
        message_length = ib%restrict_op%message_length(flag, igridf)

        ! Check that there is sufficient space in the receive buffer
        IF (recv_counter + message_length > idim_mg_bufs) THEN
            CALL errr(__FILE__, __LINE__)
        END IF

        CALL MPI_Irecv(recvbuf(recv_counter + 1), message_length, &
            mglet_mpi_real, iprocf, igridf, MPI_COMM_WORLD, recvReqs(nRecv))

        recvOffsetLength(1, igridf) = recv_counter
        recvOffsetLength(2, igridf) = message_length
        recv_counter = recv_counter + message_length
    END SUBROUTINE restrict_recv


    SUBROUTINE restrict_send(igridf, igridc)
        ! Function to pack to information into a message and to send it.
        ! In a pattern, values are extracted from the fine fields and
        ! stored within a 1-dimensional buffer array for sending
        ! (see: subroutine restrict_send_packaging)
        ! MPI_Isend is called.

        INTEGER(intk), INTENT(in) :: igridf
        INTEGER(intk), INTENT(in) :: igridc

        INTEGER(intk) :: iprocc
        INTEGER(intk) :: ip3
        INTEGER(intk) :: ii, jj, kk
        INTEGER(intk) :: message_length

        message_length = ib%restrict_op%message_length(flag, igridf)

        IF (send_counter + message_length > idim_mg_bufs) THEN
            CALL errr(__FILE__, __LINE__)
        END IF

        CALL get_mgdims(kk, jj, ii, igridf)
        CALL get_ip3(ip3, igridf)
        CALL ib%restrict_op%restrict(kk, jj, ii, ffg(ip3:ip3+kk*jj*ii-1), &
            sendbuf(send_counter+1:send_counter+message_length), &
            flag, igridf)

        iprocc = idprocofgrd(igridc)
        CALL MPI_Isend(sendbuf(send_counter + 1), message_length, &
            mglet_mpi_real, iprocc, igridf,  MPI_COMM_WORLD, &
            sendReqs(nSend))

        send_counter = send_counter + message_length
    END SUBROUTINE restrict_send


    SUBROUTINE restrict_finish(igridf)
        ! Subroutine arguments
        INTEGER(intk), INTENT(in) :: igridf

        ! Local variables
        INTEGER(intk) :: igridc, position, length
        INTEGER(intk) :: kk, jj, ii, ip3

        position = recvOffsetLength(1, igridf)
        length = recvOffsetLength(2, igridf)

        igridc = iparent(igridf)
        CALL get_mgdims(kk, jj, ii, igridc)
        CALL get_ip3(ip3, igridc)

        IF (flag == 'N') THEN
            CALL restrict_recieve_open_n(kk, jj, ii, fcg(ip3:ip3+kk*jj*ii-1), &
                recvbuf(position+1:position+length), igridf)
        ELSE
            CALL restrict_recieve_open(kk, jj, ii, fcg(ip3:ip3+kk*jj*ii-1), &
                recvbuf(position+1:position+length), igridf)
        END IF
    END SUBROUTINE restrict_finish


    SUBROUTINE init_ftoc()
        ! Function to initialize arrays and data types.
        ! After successful execution, the module varaible "is_init" is set to
        ! true
        INTEGER :: igrid, iprocc, ipar

        CALL set_timer(220, "FTOC")
        CALL set_timer(221, "FTOC_BEGIN")
        CALL set_timer(222, "FTOC_END")

        IF (.NOT. is_init) THEN
            ALLOCATE(sendReqs(nMyGrids))
            ALLOCATE(recvReqs(nMyGrids*max_childs))
            ALLOCATE(recvGrids(nMyGrids*max_childs))
            ALLOCATE(recvOffsetLength(2, ngrid))
        END IF

        nRecv = 0
        nSend = 0

        ! Make all Send- and Recv-types
        DO igrid = 1, ngrid
            ipar = iparent(igrid)
            IF (ipar /= 0 ) THEN
                iprocc = idprocofgrd(ipar)

                IF (myid == iprocc) THEN
                    nRecv = nRecv + 1
                    IF (nRecv > nMyGrids*max_childs) THEN
                        CALL errr(__FILE__, __LINE__)
                    END IF
                END IF
            END IF
        END DO

        is_init = .TRUE.
        in_progress = .FALSE.

        ! Nullify pointers
        nullify(ffg)
        nullify(fcg)
    END SUBROUTINE init_ftoc


    SUBROUTINE finish_ftoc()
        ! Function to deallocate arrays.
        ! After successful execution, the module varaible "is_init" is set
        ! to false.
        IF (is_init .NEQV. .TRUE.) THEN
            RETURN
        END IF

        is_init = .FALSE.

        DEALLOCATE(sendReqs)
        DEALLOCATE(recvReqs)
        DEALLOCATE(recvGrids)
        DEALLOCATE(recvOffsetLength)
    END SUBROUTINE finish_ftoc


    SUBROUTINE restrict_recieve_open(kk, jj, ii, fc, buffer, igridf)
        ! Subroutine arguments
        INTEGER(intk), INTENT(in) :: kk, jj, ii
        REAL(realk), INTENT(inout) :: fc(kk, jj, ii)
        REAL(realk), CONTIGUOUS, INTENT(in) :: buffer(:)
        INTEGER(intk), INTENT(in) :: igridf

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


    SUBROUTINE restrict_recieve_open_n(kk, jj, ii, fc, buffer, igridf)
        ! Subroutine arguments
        INTEGER(intk), INTENT(in) :: kk, jj, ii
        REAL(realk), INTENT(inout) :: fc(kk, jj, ii)
        REAL(realk), CONTIGUOUS, INTENT(in) :: buffer(:)
        INTEGER(intk), INTENT(in) :: igridf

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
