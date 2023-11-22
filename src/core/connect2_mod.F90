MODULE connect2_mod
    USE precision_mod
    USE MPI_f08
    USE commbuf_mod, ONLY: sendbuf, recvbuf, isendbuf, irecvbuf
    USE err_mod, ONLY: errr
    USE timer_mod, ONLY: start_timer, set_timer, stop_timer
    USE pointers_mod, ONLY: get_ip3
    USE grids_mod, ONLY: mygrids, nmygrids, level, idprocofgrd, itypboconds, &
        maxlevel, minlevel, get_neighbours, get_mgdims
    USE comms_mod, ONLY: myid, numprocs
    USE field_mod

    IMPLICIT NONE (type, external)
    PRIVATE

    ! Maximum number of connections on one single process, either
    ! outgoing or incomming, on any single grid level
    INTEGER(intk) :: maxConns

    ! Lists that hold the Send and Recv connections per grid level
    ! This list must be pre-compiled before the first call to 'connect'
    ! is being made. The reason for this is because it is expensive
    ! to compute every single time a connect is being made.
    !
    ! Dimensions contain:
    !   Dim 1: Information about a specific connection
    !   Dim 2: The different connections
    !
    ! The information in the first dimension is sorted as follows:
    !   Field 1: Rank of receiving process
    !   Field 2: Rank of sending process
    !   Field 3: ID of receiving grid
    !   Field 4: ID of sending grid
    !   Field 5: Which face (1..26) to receive
    !   Field 6: Which face (1..26) to send
    !   Field 7: Message tag (for MPI)
    !   Field 8: Geometry exchange flag
    INTEGER(intk), ALLOCATABLE :: sendConns(:,:), recvConns(:,:)

    ! Lists that hold the send and receive request arrays
    TYPE(MPI_Request), ALLOCATABLE :: sendReqs(:), recvReqs(:)

    ! Lists that hold the messages that are ACTUALLY sendt and received
    INTEGER(intk) :: nSend, nRecv, nRecvFaces
    INTEGER(int32), ALLOCATABLE :: sendList(:), recvList(:)
    INTEGER(intk), ALLOCATABLE :: recvIdxList(:,:)

    ! Number of send and receive connections
    INTEGER(intk) :: iSend = 0, iRecv = 0

    ! Counters for send- and receive operations (for locations in
    ! send and receive buffers)
    INTEGER(intk) :: sendCounter, recvCounter

    ! Number of variables per cell to exchange
    ! Number of planes to exchange
    INTEGER(intk) :: nVars, nplane

    ! Exchange faces with geom only
    LOGICAL :: geometry

    ! Forward-flag. If 1, only send data forward (in positive
    ! coordinate direction) and if -1 only send data backward
    ! NOt to be used with corners
    INTEGER(intk) :: fwd

    ! If true, exchange only surface normal component of
    ! vector field
    LOGICAL :: sn

    ! Exchange corner lines/points as well
    LOGICAL :: vertices

    ! Flag to indicate very special behaviour
    CHARACTER(len=1) :: flag

    ! Accumulated message length
    INTEGER(int32) :: messageLength

    ! Minimum and maximul level to connect
    INTEGER(intk) :: minConLvl, maxConLvl

    ! Variable to indicate if the connection information has
    ! been created.
    LOGICAL :: isInit = .FALSE.

    ! A connect can be either of either REAL or INTEGER, never both in the
    ! same call. When connect_integer is .TRUE. we use isendbuf, irecvbuf
    ! instead of sendbuf and recvbuf
    LOGICAL :: connect_integer = .FALSE.

    ! Fields
    REAL(realk), POINTER, CONTIGUOUS :: u(:), v(:), w(:), &
        p1(:), p2(:), p3(:)
    INTEGER(ifk), POINTER, CONTIGUOUS :: ui(:), vi(:), wi(:), &
        pi1(:), pi2(:), pi3(:)

    INTEGER(intk), PARAMETER :: facelist(4,26) = RESHAPE((/ &
        1, 1, 0, 0, &
        1, 2, 0, 0, &
        1, 3, 0, 0, &
        1, 4, 0, 0, &
        1, 5, 0, 0, &
        1, 6, 0, 0, &
        2, 1, 3, 0, &
        2, 1, 4, 0, &
        2, 1, 5, 0, &
        2, 1, 6, 0, &
        2, 2, 3, 0, &
        2, 2, 4, 0, &
        2, 2, 5, 0, &
        2, 2, 6, 0, &
        2, 3, 5, 0, &
        2, 3, 6, 0, &
        2, 4, 5, 0, &
        2, 4, 6, 0, &
        3, 1, 3, 5, &
        3, 1, 3, 6, &
        3, 1, 4, 5, &
        3, 1, 4, 6, &
        3, 2, 3, 5, &
        3, 2, 3, 6, &
        3, 2, 4, 5, &
        3, 2, 4, 6 /), SHAPE(facelist))

    INTEGER(intk), PARAMETER :: facenbr(26) = (/ &
        2, &
        1, &
        4, &
        3, &
        6, &
        5, &
        12, &
        11, &
        14, &
        13, &
        8, &
        7, &
        10, &
        9, &
        18, &
        17, &
        16, &
        15, &
        26, &
        25, &
        24, &
        23, &
        22, &
        21, &
        20, &
        19 /)

    ! These patterns come from a Python program, however, it was discovered,
    ! that the GC fcorr-stencils need in inner corners a special rescue
    ! neighbor. fcorr-stencils exist from 3..ii-1 and so on. In the inner
    ! corners in a three-grid configuration, in the front-left (8),
    ! front-top (10) and right-top (16), these need data to come from
    ! the grids that treat the face-normal velocity also face-normal. This is
    ! because the face-normal and face-tangential velocities are
    ! treated diffrently.
    !
    ! The original un-altered order of these faces are in the comment
    ! behind the adapted ones.
    INTEGER(intk), PARAMETER :: rescue_dir(7, 26) = RESHAPE((/ &
        1, 0, 0, 0, 0, 0, 0, &
        2, 0, 0, 0, 0, 0, 0, &
        3, 0, 0, 0, 0, 0, 0, &
        4, 0, 0, 0, 0, 0, 0, &
        5, 0, 0, 0, 0, 0, 0, &
        6, 0, 0, 0, 0, 0, 0, &
        7, 1, 3, 0, 0, 0, 0, &
        8, 4, 1, 0, 0, 0, 0, &  ! 8, 1, 4, 0, 0, 0, 0, &
        9, 1, 5, 0, 0, 0, 0, &
        10, 6, 1, 0, 0, 0, 0, &  ! 10, 1, 6, 0, 0, 0, 0, &
        11, 2, 3, 0, 0, 0, 0, &
        12, 2, 4, 0, 0, 0, 0, &
        13, 2, 5, 0, 0, 0, 0, &
        14, 2, 6, 0, 0, 0, 0, &
        15, 3, 5, 0, 0, 0, 0, &
        16, 6, 3, 0, 0, 0, 0, &  ! 16, 3, 6, 0, 0, 0, 0, &
        17, 4, 5, 0, 0, 0, 0, &
        18, 4, 6, 0, 0, 0, 0, &
        19, 7, 9, 15, 1, 3, 5, &
        20, 7, 10, 16, 1, 3, 6, &
        21, 8, 9, 17, 1, 4, 5, &
        22, 8, 10, 18, 1, 4, 6, &
        23, 11, 13, 15, 2, 3, 5, &
        24, 11, 14, 16, 2, 3, 6, &
        25, 12, 13, 17, 2, 4, 5, &
        26, 12, 14, 18, 2, 4, 6 /), SHAPE(rescue_dir))

    INTEGER(intk), PARAMETER :: rescue_nbr(7, 26) = RESHAPE((/ &
        2, 0, 0, 0, 0, 0, 0, &
        1, 0, 0, 0, 0, 0, 0, &
        4, 0, 0, 0, 0, 0, 0, &
        3, 0, 0, 0, 0, 0, 0, &
        6, 0, 0, 0, 0, 0, 0, &
        5, 0, 0, 0, 0, 0, 0, &
        12, 11, 8, 0, 0, 0, 0, &
        11, 7, 12, 0, 0, 0, 0, &  ! 11, 12, 7, 0, 0, 0, 0, &
        14, 13, 10, 0, 0, 0, 0, &
        13, 9, 14, 0, 0, 0, 0, &  ! 13, 14, 9, 0, 0, 0, 0, &
        8, 7, 12, 0, 0, 0, 0, &
        7, 8, 11, 0, 0, 0, 0, &
        10, 9, 14, 0, 0, 0, 0, &
        9, 10, 13, 0, 0, 0, 0, &
        18, 17, 16, 0, 0, 0, 0, &
        17, 15, 18, 0, 0, 0, 0, &  ! 17, 18, 15, 0, 0, 0, 0, &
        16, 15, 18, 0, 0, 0, 0, &
        15, 16, 17, 0, 0, 0, 0, &
        26, 25, 24, 22, 23, 21, 20, &
        25, 26, 23, 21, 24, 22, 19, &
        24, 23, 26, 20, 25, 19, 22, &
        23, 24, 25, 19, 26, 20, 21, &
        22, 21, 20, 26, 19, 25, 24, &
        21, 22, 19, 25, 20, 26, 23, &
        20, 19, 22, 24, 21, 23, 26, &
        19, 20, 21, 23, 22, 24, 25 /), SHAPE(rescue_nbr))

    PUBLIC :: connect, connect_real, init_connect2, finish_connect2

CONTAINS

    ! Main connect functions
    SUBROUTINE connect(ilevel, layers, v1, v2, v3, &
            s1, s2, s3, geom, corners, normal, forward, ityp, minlvl, maxlvl)

        ! Subroutine arguments
        INTEGER(intk), INTENT(in), OPTIONAL :: ilevel, layers
        CLASS(basefield_t), TARGET, OPTIONAL, INTENT(inout) :: &
            v1, v2, v3, s1, s2, s3
        LOGICAL, OPTIONAL, INTENT(in) :: geom, corners, normal
        INTEGER(intk), OPTIONAL, INTENT(in) :: forward
        CHARACTER(len=1), OPTIONAL, INTENT(in) :: ityp
        INTEGER(intk), INTENT(in), OPTIONAL :: minlvl, maxlvl

        ! Local variables
        LOGICAL :: has_v1, has_v2, has_v3, has_s1, has_s2, has_s3

        has_v1 = .FALSE.
        has_v2 = .FALSE.
        has_v3 = .FALSE.
        has_s1 = .FALSE.
        has_s2 = .FALSE.
        has_s3 = .FALSE.

        IF (PRESENT(v1)) THEN
            has_v1 = .TRUE.
            SELECT TYPE (v1)
            TYPE IS (field_t)
                u => v1%arr
            TYPE IS (intfield_t)
                ui => v1%arr
            END SELECT
        END IF

        IF (PRESENT(v2)) THEN
            has_v2 = .TRUE.
            SELECT TYPE (v2)
            TYPE IS (field_t)
                v => v2%arr
            TYPE IS (intfield_t)
                vi => v2%arr
            END SELECT
        END IF

        IF (PRESENT(v3)) THEN
            has_v3 = .TRUE.
            SELECT TYPE (v3)
            TYPE IS (field_t)
                w => v3%arr
            TYPE IS (intfield_t)
                wi => v3%arr
            END SELECT
        END IF

        IF (PRESENT(s1)) THEN
            has_s1 = .TRUE.
            SELECT TYPE (s1)
            TYPE IS (field_t)
                p1 => s1%arr
            TYPE IS (intfield_t)
                pi1 => s1%arr
            END SELECT
        END IF

        IF (PRESENT(s2)) THEN
            has_s2 = .TRUE.
            SELECT TYPE (s2)
            TYPE IS (field_t)
                p2 => s2%arr
            TYPE IS (intfield_t)
                pi2 => s2%arr
            END SELECT
        END IF

        IF (PRESENT(s3)) THEN
            has_s3 = .TRUE.
            SELECT TYPE (s3)
            TYPE IS (field_t)
                p3 => s3%arr
            TYPE IS (intfield_t)
                pi3 => s3%arr
            END SELECT
        END IF

        CALL connect_impl(ilevel, layers, has_v1, has_v2, has_v3, &
            has_s1, has_s2, has_s3, geom, corners, normal, forward, ityp, &
            minlvl, maxlvl)
    END SUBROUTINE connect


    ! Main connect functions
    SUBROUTINE connect_real(ilevel, layers, v1, v2, v3, &
        s1, s2, s3, geom, corners, normal, forward, ityp)

        ! Subroutine arguments
        INTEGER(intk), INTENT(in), OPTIONAL :: ilevel, layers
        REAL(realk), TARGET, OPTIONAL, CONTIGUOUS, INTENT(inout) :: &
            v1(:), v2(:), v3(:), s1(:), s2(:), s3(:)
        LOGICAL, OPTIONAL, INTENT(in) :: geom, corners, normal
        INTEGER(intk), OPTIONAL, INTENT(in) :: forward
        CHARACTER(len=1), OPTIONAL, INTENT(in) :: ityp

        ! Local variables
        LOGICAL :: has_v1, has_v2, has_v3, has_s1, has_s2, has_s3

        has_v1 = .FALSE.
        has_v2 = .FALSE.
        has_v3 = .FALSE.
        has_s1 = .FALSE.
        has_s2 = .FALSE.
        has_s3 = .FALSE.

        IF (PRESENT(v1)) THEN
            has_v1 = .TRUE.
            u => v1
        END IF
        IF (PRESENT(v2)) THEN
            has_v2 = .TRUE.
            v => v2
        END IF
        IF (PRESENT(v3)) THEN
            has_v3 = .TRUE.
            w => v3
        END IF
        IF (PRESENT(s1)) THEN
            has_s1 = .TRUE.
            p1 => s1
        END IF
        IF (PRESENT(s2)) THEN
            has_s2 = .TRUE.
            p2 => s2
        END IF
        IF (PRESENT(s3)) THEN
            has_s3 = .TRUE.
            p3 => s3
        END IF

        CALL connect_impl(ilevel, layers, has_v1, has_v2, has_v3, &
            has_s1, has_s2, has_s3, geom, corners, normal, forward, ityp)
    END SUBROUTINE connect_real


    ! Main connect function
    SUBROUTINE connect_impl(ilevel, layers, has_v1, has_v2, has_v3, &
            has_s1, has_s2, has_s3, geom, corners, normal, forward, ityp, &
            minlvl, maxlvl)

        ! The ilevel and nplane parameters are the only required
        ! parameters
        INTEGER(intk), INTENT(in), OPTIONAL :: ilevel, layers

        ! I am surprised this is allowed as these are not optional...
        LOGICAL :: has_v1, has_v2, has_v3, has_s1, has_s2, has_s3

        ! Optional parameters to control special behaviour
        LOGICAL, OPTIONAL, INTENT(in) :: geom, corners, normal
        INTEGER(intk), OPTIONAL, INTENT(in) :: forward
        CHARACTER(len=1), OPTIONAL, INTENT(in) :: ityp
        INTEGER(intk), INTENT(in), OPTIONAL :: minlvl, maxlvl

        ! Local variables
        LOGICAL :: has_real_arg
        LOGICAL :: has_int_arg

        CALL start_timer(150)

        ! Check if the connection information has been created
        IF (isInit .EQV. .FALSE.) THEN
            WRITE(*,*) "'connect' not initialized"
            CALL errr(__FILE__, __LINE__)
        END IF

        ! Check that no other transfers are in progress
        IF (nSend > 0 .OR. nRecv > 0) THEN
            WRITE(*,*) "Other transfer in progress."
            CALL errr(__FILE__, __LINE__)
        END IF

        ! Check that number of connect layers are either 1 or 2
        IF (PRESENT(layers)) THEN
            IF ( .NOT. (layers  ==  1 .OR. layers  ==  2)) THEN
                WRITE(*,*) "Invalid layers=", layers
                CALL errr(__FILE__, __LINE__)
            END IF
            nplane = layers
        ELSE
            nplane = 1
        END IF

        ! If one vector argument is given, check that all three is present.
        nVars = 0
        IF (has_v1) THEN
            IF (.NOT. (has_v2 .AND. has_v3)) THEN
                WRITE(*,*) "If one vector arg is present, all three must be present."
                CALL errr(__FILE__, __LINE__)
            END IF
            nVars = nVars + 3
        ELSE IF (has_v2 .OR. has_v3) THEN
            WRITE(*,*) "If one vector arg is present, all three must be present."
            CALL errr(__FILE__, __LINE__)
        END IF

        ! Set pointers to scalars
        IF (has_s1) THEN
            nVars = nVars + 1
        END IF
        IF (has_s2) THEN
            nVars = nVars + 1
        END IF
        IF (has_s3) THEN
            nVars = nVars + 1
        END IF

        geometry = .FALSE.
        IF (PRESENT(geom)) THEN
            IF (geom .EQV. .TRUE.) THEN
                geometry = .TRUE.
            END IF
        END IF

        vertices = .FALSE.
        IF (PRESENT(corners)) THEN
            IF (corners .EQV. .TRUE.) THEN
                vertices = .TRUE.
            END IF
        END IF

        fwd = 0
        IF (PRESENT(forward)) THEN
            IF (vertices) THEN
                CALL errr(__FILE__, __LINE__)
            END IF
            IF (forward > 0) THEN
                fwd = 1
            ELSE IF (forward < 0) THEN
                fwd = -1
            END IF
        END IF

        ! If exchanging surface normal component, check that a vector is
        ! given.
        sn = .FALSE.
        IF (PRESENT(normal)) THEN
            IF (vertices) THEN
                CALL errr(__FILE__, __LINE__)
            END IF
            IF (normal .eqv. .true.) THEN
                IF (has_v1) THEN
                    nVars = nVars - 2
                    sn = .TRUE.
                ELSE
                    WRITE(*,*) "normal=.TRUE. require a vector v1, v2, v3."
                    CALL errr(__FILE__, __LINE__)
                END IF
            END IF
        END IF

        ! Check special flags
        flag = ' '
        IF (PRESENT(ityp)) THEN
            IF ((ityp == 'Y')) THEN
                flag = ityp
            ELSE
                write(*,*) "Invalid ITYP=", ityp
                CALL errr(__FILE__, __LINE__)
            END IF
        END IF

        ! Not specifying any fields would be very strange
        IF (nVars == 0) THEN
            WRITE(*,*) "You have not specified any fields to exchange."
            CALL errr(__FILE__, __LINE__)
        END IF
        ! TODO: Check why? Maybe because of buffer limitations?
        ! In the meantime, comment out
        ! IF (nVars > maxnVars) THEN
        !     WRITE(*,*) "You specified too many fields to exchange."
        !     CALL errr(__FILE__, __LINE__)
        ! END IF

        ! Either ilevel or both minlvl and maxlvl
        IF (PRESENT(ilevel) .AND. (PRESENT(minlvl) .OR. PRESENT(maxlvl))) THEN
            CALL errr(__FILE__, __LINE__)
        END IF
        IF (PRESENT(minlvl) .NEQV. PRESENT(maxlvl)) THEN
            CALL errr(__FILE__, __LINE__)
        END IF

        minconlvl = minlevel
        maxconlvl = maxlevel
        IF (PRESENT(ilevel)) THEN
            minconlvl = ilevel
            maxconlvl = ilevel
        ELSE IF (PRESENT(minlvl) .AND. PRESENT(maxlvl)) THEN
            minconlvl = minlvl
            maxconlvl = maxlvl
        END IF
        IF ((minconlvl < minlevel) .OR. (maxconlvl > maxlevel)) THEN
            WRITE(*, *) "Invalid levels: ", minconlvl, maxconlvl
            CALL errr(__FILE__, __LINE__)
        END IF

        ! Sanity check if integers/realk
        ! can only connect either INTEGER _or_ REAL
        has_real_arg = .FALSE.
        has_int_arg = .FALSE.
        IF (ASSOCIATED(u) .OR. ASSOCIATED(v) .OR. ASSOCIATED(w) .OR. &
                ASSOCIATED(p1) .OR. ASSOCIATED(p2) .OR. ASSOCIATED(p3)) THEN
            has_real_arg = .TRUE.
        END IF
        IF (ASSOCIATED(ui) .OR. ASSOCIATED(vi) .OR. ASSOCIATED(wi) .OR. &
                ASSOCIATED(pi1) .OR. ASSOCIATED(pi2) .OR. ASSOCIATED(pi3)) THEN
            has_int_arg = .TRUE.
        END IF

        IF (has_real_arg .EQV. has_int_arg) THEN
            CALL errr(__FILE__, __LINE__)
        END IF

        connect_integer = .FALSE.
        IF (has_int_arg) connect_integer = .TRUE.

        ! Exchnage data
        CALL recv_all()
        CALL send_all()
        CALL process_bufs()

        ! Clear counters and unset pointers. This is important to avoid
        ! problems at next function call
        nRecv = 0
        nSend = 0
        nplane = 0
        nVars = 0

        sn = .FALSE.
        flag = ' '
        fwd = 0
        vertices = .FALSE.
        geometry = .FALSE.

        NULLIFY(u)
        NULLIFY(v)
        NULLIFY(w)
        NULLIFY(p1)
        NULLIFY(p2)
        NULLIFY(p3)

        NULLIFY(ui)
        NULLIFY(vi)
        NULLIFY(wi)
        NULLIFY(pi1)
        NULLIFY(pi2)
        NULLIFY(pi3)

        CALL stop_timer(150)
    END SUBROUTINE connect_impl


    FUNCTION decide(i, list) RESULT(exchange)
        INTEGER(intk), INTENT(in) :: i
        INTEGER(intk), INTENT(in) :: list(:,:)

        INTEGER(intk) :: ifacerecv, ilevel
        LOGICAL :: exchange

        exchange = .TRUE.
        ifacerecv = list(5, i)

        IF (i > SIZE(list, 2)) THEN
            exchange = .FALSE.
            RETURN
        END IF

        ! Levels
        ilevel = level(list(3, i))
        IF (ilevel > maxConLvl .OR. ilevel < minConLvl) THEN
            exchange = .FALSE.
        END IF

        ! Only exchange geometry
#ifndef _IB_CUTCELL_MOVING_
        IF (geometry .AND. list(8, i) == 0) THEN
            exchange = .FALSE.
        END IF
#endif

        ! Corners
        IF ((.NOT. vertices) .AND. (ifacerecv > 6)) THEN
            exchange = .FALSE.
        END IF

        ! Forward
        IF (fwd == 1 .AND. &
            (ifacerecv == 2 .OR. ifacerecv == 4 .OR. ifacerecv == 6)) THEN
            exchange = .FALSE.
        ELSE IF (fwd == -1 .AND. &
            (ifacerecv == 1 .OR. ifacerecv == 3 .OR. ifacerecv == 5)) THEN
            exchange = .FALSE.
        END IF

        RETURN
    END FUNCTION decide


    ! Perform all Recv-calls
    SUBROUTINE recv_all()

        INTEGER(intk) :: i, iprocnbr, igrid, iface, faceArea
        LOGICAL :: exchange

        ! Post all receive calls
        recvCounter = 0
        messageLength = 0
        nRecv = 0
        recvIdxList = 0

        DO i = 1, iRecv
            exchange = decide(i, recvConns)
            iprocnbr = recvConns(2, i)

            ! Communication with self is handled specially in
            ! send_all - nothing to do here
            IF (iprocnbr == myid) THEN
                CYCLE
            END IF

            IF (exchange) THEN
                igrid         = recvConns(3, i)
                iface         = recvConns(5, i)

                faceArea      = face_area(igrid, iface)
                nRecvFaces = nRecvFaces + 1
                recvIdxList(1, i) = iprocnbr
                recvIdxList(2, i) = nVars*faceArea
                recvIdxList(3, i) = recvCounter + messageLength
                messageLength = messageLength + nVars*faceArea
            END IF

            IF (messageLength > 0) THEN
                IF (i == iRecv) THEN
                    CALL post_recv(iprocnbr)
                ELSE IF (recvConns(2, i + 1) /= iprocnbr) THEN
                    CALL post_recv(iprocnbr)
                END IF
            END IF
        END DO

    END SUBROUTINE recv_all


    ! Perform a single Recv
    SUBROUTINE post_recv(iprocnbr)
        ! Subroutine arguments
        INTEGER(int32), INTENT(in) :: iprocnbr

        nRecv = nRecv + 1
        recvList(nRecv) = iprocnbr

        IF (connect_integer) THEN
            IF (recvcounter + messagelength > SIZE(irecvbuf)) THEN
                CALL errr(__FILE__, __LINE__)
            END IF

            CALL MPI_Irecv(irecvbuf(recvcounter+1) , messagelength, &
                mglet_mpi_ifk, iprocnbr, 1, MPI_COMM_WORLD, recvreqs(nrecv))
        ELSE
            IF (recvcounter + messagelength > SIZE(recvbuf)) THEN
                CALL errr(__FILE__, __LINE__)
            END IF

            CALL MPI_Irecv(recvbuf(recvcounter+1) , messagelength, &
                mglet_mpi_real, iprocnbr, 1, MPI_COMM_WORLD, recvreqs(nrecv))
        END IF

        recvcounter = recvcounter + messagelength
        messagelength = 0
    END SUBROUTINE post_recv


    ! Perform all send calls
    SUBROUTINE send_all()
        INTEGER(intk) :: i, iprocnbr, igrid, ifacerecv, faceArea
        LOGICAL :: exchange

        ! Pack all buffers and send data
        sendCounter = 0
        messageLength = 0
        nSend = 0

        DO i = 1, iSend
            exchange = decide(i, sendConns)
            iprocnbr = sendConns(1, i)

            ! Communication with self copies directly from source to
            ! destination grid - then skip the rest
            IF (iprocnbr == myid .AND. exchange) THEN
                CALL connect_self(i)
                CYCLE
            END IF

            IF (exchange) THEN
                igrid         = sendConns(3, i)
                ifacerecv     = sendConns(5, i)
                faceArea      = face_area(igrid, ifacerecv)

                CALL write_buffer(i)

                messageLength = messageLength + nVars*faceArea
            END IF

            IF (messageLength > 0) THEN
                IF (i == iSend) THEN
                    CALL post_send(iprocnbr)
                ELSE IF (sendConns(1, i + 1) /= iprocnbr) THEN
                    CALL post_send(iprocnbr)
                END IF
            END IF
        END DO
    END SUBROUTINE send_all


    ! Perform a single send call
    SUBROUTINE post_send(iprocnbr)
        ! Subroutine arguments
        INTEGER(int32), INTENT(in) :: iprocnbr

        ! Local variables
        ! none...

        nSend = nSend + 1
        sendList(nSend) = iprocnbr

        IF (connect_integer) THEN
            CALL MPI_Isend(isendbuf(sendCounter + 1), messageLength, &
                mglet_mpi_ifk, iprocnbr, 1, MPI_COMM_WORLD, sendReqs(nSend))
        ELSE
            CALL MPI_Isend(sendbuf(sendCounter + 1), messageLength, &
                mglet_mpi_real, iprocnbr, 1, MPI_COMM_WORLD, sendReqs(nSend))
        END IF

        sendCounter = sendCounter + messageLength
        messageLength = 0
    END SUBROUTINE post_send


    ! Write Send buffers
    !
    ! Write the relevant fields into the send buffers
    SUBROUTINE write_buffer(sendId)

        ! Input parameter
        INTEGER(int32), INTENT(in) :: sendId

        ! Grid dimensions and pointers
        INTEGER(intk) :: kk, jj, ii
        INTEGER(intk) :: ip3

        ! Indices of start- and stop of iteration over boundary face
        INTEGER(intk) :: istart, istop, jstart, jstop, kstart, kstop

        ! Grid to send from
        ! Must be intk because it intreface with MGLET
        INTEGER(intk) :: igrid, ifacerecv, ifacesend

        ! Message sizes
        ! Must be int32 because it iterface with MPI
        INTEGER(int32) :: thisMessageLength, faceArea
        INTEGER(int32) :: iCount, pos, offset

        ! Iterators
        INTEGER(intk) :: i, j, k

        ! Flags to indicate exchange of U, V, W
        LOGICAL :: exU, exV, exW

        ! Set variables from send table
        igrid         = sendConns(4, sendId)
        ifacerecv     = sendConns(5, sendId)
        ifacesend     = sendConns(6, sendId)

        ! Get grid dimentsions and pointers
        CALL get_mgdims(kk, jj, ii, igrid)
        CALL get_ip3(ip3, igrid)

        ! Get start- and stop indices of grid
        CALL start_and_stop(igrid, facenbr(ifacerecv), &
            istart, istop, jstart, jstop, kstart, kstop, nghost=1)
        CALL corr_start_stop(igrid, ifacesend, ifacerecv, &
            istart, istop, jstart, jstop, kstart, kstop)

        faceArea = (istop-istart+1)*(jstop-jstart+1)*(kstop-kstart+1)
        thisMessageLength = nVars*faceArea

        ! Check that buffer does not overflow
        IF (sendcounter + messagelength + thismessagelength > SIZE(sendbuf)) THEN
            CALL errr(__FILE__, __LINE__)
        END IF

        ! Reset message size counter
        offset = sendCounter + messageLength
        iCount = 0

        ! Fill buffers
        IF (connect_integer) THEN
            ! INTEGER
            exU = (sn .AND. ifacerecv < 3) .OR. (.NOT. sn)
            IF (ASSOCIATED(ui) .AND. exU) THEN
                DO i = istart, istop
                    DO j = jstart, jstop
                        DO k = kstart, kstop
                            iCount = iCount + 1
                            pos = jj*kk*(i-1) + kk*(j-1) + (k-1)
                            isendbuf(offset + iCount) = ui(ip3 + pos)
                        END DO
                    END DO
                END DO
            END IF
            exV = (sn .AND. (ifacerecv > 2 .AND. ifacerecv < 5)) .OR. (.NOT. sn)
            IF (ASSOCIATED(vi) .AND. exV) THEN
                DO i = istart, istop
                    DO j = jstart, jstop
                        DO k = kstart, kstop
                            iCount = iCount + 1
                            pos = jj*kk*(i-1) + kk*(j-1) + (k-1)
                            isendbuf(offset + iCount) = vi(ip3 + pos)
                        END DO
                    END DO
                END DO
            END IF
            exW = (sn .AND. ifacerecv > 4) .OR. (.NOT. sn)
            IF (ASSOCIATED(wi) .AND. exW) THEN
                DO i = istart, istop
                    DO j = jstart, jstop
                        DO k = kstart, kstop
                            iCount = iCount + 1
                            pos = jj*kk*(i-1) + kk*(j-1) + (k-1)
                            isendbuf(offset + iCount) = wi(ip3 + pos)
                        END DO
                    END DO
                END DO
            END IF
            IF (ASSOCIATED(pi1)) THEN
                DO i = istart, istop
                    DO j = jstart, jstop
                        DO k = kstart, kstop
                            iCount = iCount + 1
                            pos = jj*kk*(i-1) + kk*(j-1) + (k-1)
                            isendbuf(offset + iCount) = pi1(ip3 + pos)
                        END DO
                    END DO
                END DO
            END IF
            IF (ASSOCIATED(pi2)) THEN
                DO i = istart, istop
                    DO j = jstart, jstop
                        DO k = kstart, kstop
                            iCount = iCount + 1
                            pos = jj*kk*(i-1) + kk*(j-1) + (k-1)
                            isendbuf(offset + iCount) = pi2(ip3 + pos)
                        END DO
                    END DO
                END DO
            END IF
            IF (ASSOCIATED(pi3)) THEN
                DO i = istart, istop
                    DO j = jstart, jstop
                        DO k = kstart, kstop
                            iCount = iCount + 1
                            pos = jj*kk*(i-1) + kk*(j-1) + (k-1)
                            isendbuf(offset + iCount) = pi3(ip3 + pos)
                        END DO
                    END DO
                END DO
            END IF
        ELSE
            ! REAL
            exU = (sn .AND. ifacerecv < 3) .OR. (.NOT. sn)
            IF (ASSOCIATED(u) .AND. exU) THEN
                DO i = istart, istop
                    DO j = jstart, jstop
                        DO k = kstart, kstop
                            iCount = iCount + 1
                            pos = jj*kk*(i-1) + kk*(j-1) + (k-1)
                            sendbuf(offset + iCount) = u(ip3 + pos)
                        END DO
                    END DO
                END DO
            END IF
            exV = (sn .AND. (ifacerecv > 2 .AND. ifacerecv < 5)) .OR. (.NOT. sn)
            IF (ASSOCIATED(v) .AND. exV) THEN
                DO i = istart, istop
                    DO j = jstart, jstop
                        DO k = kstart, kstop
                            iCount = iCount + 1
                            pos = jj*kk*(i-1) + kk*(j-1) + (k-1)
                            sendbuf(offset + iCount) = v(ip3 + pos)
                        END DO
                    END DO
                END DO
            END IF
            exW = (sn .AND. ifacerecv > 4) .OR. (.NOT. sn)
            IF (ASSOCIATED(w) .AND. exW) THEN
                DO i = istart, istop
                    DO j = jstart, jstop
                        DO k = kstart, kstop
                            iCount = iCount + 1
                            pos = jj*kk*(i-1) + kk*(j-1) + (k-1)
                            sendbuf(offset + iCount) = w(ip3 + pos)
                        END DO
                    END DO
                END DO
            END IF
            IF (ASSOCIATED(p1)) THEN
                DO i = istart, istop
                    DO j = jstart, jstop
                        DO k = kstart, kstop
                            iCount = iCount + 1
                            pos = jj*kk*(i-1) + kk*(j-1) + (k-1)
                            sendbuf(offset + iCount) = p1(ip3 + pos)
                        END DO
                    END DO
                END DO
            END IF
            IF (ASSOCIATED(p2)) THEN
                DO i = istart, istop
                    DO j = jstart, jstop
                        DO k = kstart, kstop
                            iCount = iCount + 1
                            pos = jj*kk*(i-1) + kk*(j-1) + (k-1)
                            sendbuf(offset + iCount) = p2(ip3 + pos)
                        END DO
                    END DO
                END DO
            END IF
            IF (ASSOCIATED(p3)) THEN
                DO i = istart, istop
                    DO j = jstart, jstop
                        DO k = kstart, kstop
                            iCount = iCount + 1
                            pos = jj*kk*(i-1) + kk*(j-1) + (k-1)
                            sendbuf(offset + iCount) = p3(ip3 + pos)
                        END DO
                    END DO
                END DO
            END IF
        END IF

        ! Check that message length was calculated correctly
        IF (thisMessageLength /= iCount) THEN
            write(*,*) "thisMessageLength:", thisMessageLength, "iCount:", iCount
            CALL errr(__FILE__, __LINE__)
        END IF
    END SUBROUTINE write_buffer


    SUBROUTINE corr_start_stop(igrid, ifacesend, ifacerecv, &
        istart, istop, jstart, jstop, kstart, kstop)

        ! Subroutine arguments
        INTEGER(intk), INTENT(in) :: igrid, ifacesend, ifacerecv
        INTEGER(intk), INTENT(inout) :: istart, istop, jstart, jstop, &
            kstart, kstop

        ! Local variables
        INTEGER(intk) :: i, j
        INTEGER(intk) :: iface

        IF (ifacerecv <= 6) RETURN

        ! Check if send- and recv- face have elementary faces in common
        DO i = 1, facelist(1, ifacesend)
            DO j = 1, facelist(1, ifacerecv)
                IF (facelist(i + 1, ifacesend) == facelist(j + 1, ifacerecv)) THEN
                    iface = facelist(i + 1, ifacesend)
                    CALL start_and_stop_face(igrid, iface, istart, istop, &
                        jstart, jstop, kstart, kstop)
                END IF
            END DO
        END DO
    END SUBROUTINE corr_start_stop


    ! Read Receive buffers
    !
    ! Write the contents of the receive buffers back in their
    ! matching fields
    SUBROUTINE read_buffer(recvId)
        ! Input parameter
        INTEGER(intk), INTENT(in) :: recvId

        ! Grid dimensions and pointers
        INTEGER(intk) :: kk, jj, ii
        INTEGER(intk) :: ip3

        ! Indices of start- and stop of iteration over boundary face
        INTEGER(intk) :: istart, istop, jstart, jstop, kstart, kstop

        ! Grid to send from
        ! Must be intk because it intreface with MGLET
        INTEGER(intk) :: igrid, ifacerecv

        ! Message sizes
        ! Must be int32 because it iterface with MPI
        INTEGER(int32) :: offset, iCount, pos

        ! Iterators
        INTEGER(intk) :: i, j, k

        ! Flags to indicate exchange of U, V, W
        LOGICAL :: exU, exV, exW

        ! Set variables from send table
        igrid         = recvConns(3, recvId)
        ifacerecv     = recvConns(5, recvId)

        ! Get grid dimentsions and pointers
        CALL get_mgdims(kk, jj, ii, igrid)
        CALL get_ip3(ip3, igrid)

        ! Get start- and stop indices of grid
        CALL start_and_stop(igrid, ifacerecv, &
            istart, istop, jstart, jstop, kstart, kstop)

        ! Zeroise message size counter
        iCount = 0
        offset = recvIdxList(3, recvId)

        ! Fill velicity field from buffer
        IF (connect_integer) THEN
            exU = (sn .AND. ifacerecv < 3) .OR. (.NOT. sn)
            IF (ASSOCIATED(ui) .AND. exU) THEN
                DO i = istart, istop
                    DO j = jstart, jstop
                        DO k = kstart, kstop
                            iCount = iCount + 1
                            pos = jj*kk*(i-1) + kk*(j-1) + (k-1)
                            ui(ip3 + pos) = irecvbuf(offset + iCount)
                        END DO
                    END DO
                END DO
            END IF
            exV = (sn .AND. (ifacerecv > 2 .AND. ifacerecv < 5)) .OR. (.NOT. sn)
            IF (ASSOCIATED(vi) .AND. exV) THEN
                DO i = istart, istop
                    DO j = jstart, jstop
                        DO k = kstart, kstop
                            iCount = iCount + 1
                            pos = jj*kk*(i-1) + kk*(j-1) + (k-1)
                            vi(ip3 + pos) = irecvbuf(offset + iCount)
                        END DO
                    END DO
                END DO
            END IF
            exW = (sn .AND. ifacerecv > 4) .OR. (.NOT. sn)
            IF (ASSOCIATED(wi) .AND. exW) THEN
                DO i = istart, istop
                    DO j = jstart, jstop
                        DO k = kstart, kstop
                            iCount = iCount + 1
                            pos = jj*kk*(i-1) + kk*(j-1) + (k-1)
                            wi(ip3 + pos) = irecvbuf(offset + iCount)
                        END DO
                    END DO
                END DO
            END IF
            IF (ASSOCIATED(pi1)) THEN
                DO i = istart, istop
                    DO j = jstart, jstop
                        DO k = kstart, kstop
                            iCount = iCount + 1
                            pos = jj*kk*(i-1) + kk*(j-1) + (k-1)
                            pi1(ip3 + pos) = irecvbuf(offset + iCount)
                        END DO
                    END DO
                END DO
            END IF
            IF (ASSOCIATED(pi2)) THEN
                DO i = istart, istop
                    DO j = jstart, jstop
                        DO k = kstart, kstop
                            iCount = iCount + 1
                            pos = jj*kk*(i-1) + kk*(j-1) + (k-1)
                            pi2(ip3 + pos) = irecvbuf(offset + iCount)
                        END DO
                    END DO
                END DO
            END IF
            IF (ASSOCIATED(pi3)) THEN
                DO i = istart, istop
                    DO j = jstart, jstop
                        DO k = kstart, kstop
                            iCount = iCount + 1
                            pos = jj*kk*(i-1) + kk*(j-1) + (k-1)
                            pi3(ip3 + pos) = irecvbuf(offset + iCount)
                        END DO
                    END DO
                END DO
            END IF
        ELSE
            exU = (sn .AND. ifacerecv < 3) .OR. (.NOT. sn)
            IF (ASSOCIATED(u) .AND. exU) THEN
                DO i = istart, istop
                    DO j = jstart, jstop
                        DO k = kstart, kstop
                            iCount = iCount + 1
                            pos = jj*kk*(i-1) + kk*(j-1) + (k-1)
                            u(ip3 + pos) = recvbuf(offset + iCount)
                        END DO
                    END DO
                END DO
            END IF
            exV = (sn .AND. (ifacerecv > 2 .AND. ifacerecv < 5)) .OR. (.NOT. sn)
            IF (ASSOCIATED(v) .AND. exV) THEN
                DO i = istart, istop
                    DO j = jstart, jstop
                        DO k = kstart, kstop
                            iCount = iCount + 1
                            pos = jj*kk*(i-1) + kk*(j-1) + (k-1)
                            v(ip3 + pos) = recvbuf(offset + iCount)
                        END DO
                    END DO
                END DO
            END IF
            exW = (sn .AND. ifacerecv > 4) .OR. (.NOT. sn)
            IF (ASSOCIATED(w) .AND. exW) THEN
                DO i = istart, istop
                    DO j = jstart, jstop
                        DO k = kstart, kstop
                            iCount = iCount + 1
                            pos = jj*kk*(i-1) + kk*(j-1) + (k-1)
                            w(ip3 + pos) = recvbuf(offset + iCount)
                        END DO
                    END DO
                END DO
            END IF
            IF (ASSOCIATED(p1)) THEN
                DO i = istart, istop
                    DO j = jstart, jstop
                        DO k = kstart, kstop
                            iCount = iCount + 1
                            pos = jj*kk*(i-1) + kk*(j-1) + (k-1)
                            p1(ip3 + pos) = recvbuf(offset + iCount)
                        END DO
                    END DO
                END DO
            END IF
            IF (ASSOCIATED(p2)) THEN
                DO i = istart, istop
                    DO j = jstart, jstop
                        DO k = kstart, kstop
                            iCount = iCount + 1
                            pos = jj*kk*(i-1) + kk*(j-1) + (k-1)
                            p2(ip3 + pos) = recvbuf(offset + iCount)
                        END DO
                    END DO
                END DO
            END IF
            IF (ASSOCIATED(p3)) THEN
                DO i = istart, istop
                    DO j = jstart, jstop
                        DO k = kstart, kstop
                            iCount = iCount + 1
                            pos = jj*kk*(i-1) + kk*(j-1) + (k-1)
                            p3(ip3 + pos) = recvbuf(offset + iCount)
                        END DO
                    END DO
                END DO
            END IF
        END IF

        ! Check that message length is calculated correctly
        IF (iCount /= recvIdxList(2, recvId)) THEN
            WRITE(*,*) "iCount:", iCount, &
                "recvIdxList(2, recvId):", recvIdxList(2, recvId)
            CALL errr(__FILE__, __LINE__)
        END IF
    END SUBROUTINE read_buffer


    ! Connect to self
    !
    ! Directly copes data from source to destination buffer
    SUBROUTINE connect_self(sendId)
        ! Input parameter
        INTEGER(int32), INTENT(in) :: sendId

        ! Grid dimensions and pointers
        INTEGER(intk) :: kk, jj, ii
        INTEGER(intk) :: ip3

        ! Grid dimensions and pointers, destination face
        INTEGER(intk) :: kk_d, jj_d, ii_d
        INTEGER(intk) :: ip3_d

        ! Indices of start- and stop of iteration over boundary face
        ! Source face
        INTEGER(intk) :: istart, istop, jstart, jstop, kstart, kstop

        ! Indices of start- and stop of iteration over boundary face
        ! Destination face
        INTEGER(intk) :: istart_d, istop_d, jstart_d, jstop_d, kstart_d, kstop_d

        ! Grid to send from
        ! Must be intk because it intreface with MGLET
        INTEGER(intk) :: igrid, igrid_d, ifacerecv, ifacesend

        ! Message sizes
        ! Must be int32 because it iterface with MPI
        INTEGER(int32) :: dest_size, source_size

        ! Flags to indicate exchange of U, V, W
        LOGICAL :: exU, exV, exW

        ! Pointers to fields
        REAL(realk), POINTER, CONTIGUOUS :: &
            src_field(:, :, :), dst_field(:, :, :)
        INTEGER(ifk), POINTER, CONTIGUOUS :: &
            src_ifield(:, :, :), dst_ifield(:, :, :)

        ! Set variables from send table
        igrid_d       = sendConns(3, sendId)
        igrid         = sendConns(4, sendId)
        ifacerecv     = sendConns(5, sendId)
        ifacesend     = sendConns(6, sendId)

        ! Get grid dimentsions and pointers
        CALL get_mgdims(kk, jj, ii, igrid)
        CALL get_ip3(ip3, igrid)

        CALL get_mgdims(kk_d, jj_d, ii_d, igrid_d)
        CALL get_ip3(ip3_d, igrid_d)

        ! Get start- and stop indices of source grid
        CALL start_and_stop(igrid, facenbr(ifacerecv), &
            istart, istop, jstart, jstop, kstart, kstop, nghost=1)
        CALL corr_start_stop(igrid, ifacesend, ifacerecv, &
            istart, istop, jstart, jstop, kstart, kstop)

        ! Get start- and stop indices of destination grid
        CALL start_and_stop(igrid_d, ifacerecv, &
            istart_d, istop_d, jstart_d, jstop_d, kstart_d, kstop_d)

        ! Sanity check of message length
        source_size = (istop-istart+1)*(jstop-jstart+1)*(kstop-kstart+1)
        dest_size = (istop_d-istart_d+1)*(jstop_d-jstart_d+1)*(kstop_d-kstart_d+1)
        IF (source_size /= dest_size) THEN
            CALL errr(__FILE__, __LINE__)
        END IF

        IF (connect_integer) THEN
            exU = (sn .AND. ifacerecv < 3) .OR. (.NOT. sn)
            IF (ASSOCIATED(ui) .AND. exU) THEN
                src_ifield(1:kk, 1:jj, 1:ii) => ui(ip3:ip3+kk*jj*ii-1)
                dst_ifield(1:kk_d, 1:jj_d, 1:ii_d) => ui(ip3_d:ip3_d+kk_d*jj_d*ii_d-1)
                dst_ifield(kstart_d:kstop_d, jstart_d:jstop_d, istart_d:istop_d) = &
                    src_ifield(kstart:kstop, jstart:jstop, istart:istop)
                NULLIFY(dst_ifield)
                NULLIFY(src_ifield)
            END IF
            exV = (sn .AND. (ifacerecv > 2 .AND. ifacerecv < 5)) .OR. (.NOT. sn)
            IF (ASSOCIATED(vi) .AND. exV) THEN
                src_ifield(1:kk, 1:jj, 1:ii) => vi(ip3:ip3+kk*jj*ii-1)
                dst_ifield(1:kk_d, 1:jj_d, 1:ii_d) => vi(ip3_d:ip3_d+kk_d*jj_d*ii_d-1)
                dst_ifield(kstart_d:kstop_d, jstart_d:jstop_d, istart_d:istop_d) = &
                    src_ifield(kstart:kstop, jstart:jstop, istart:istop)
                NULLIFY(dst_ifield)
                NULLIFY(src_ifield)
            END IF
            exW = (sn .AND. ifacerecv > 4) .OR. (.NOT. sn)
            IF (ASSOCIATED(wi) .AND. exW) THEN
                src_ifield(1:kk, 1:jj, 1:ii) => wi(ip3:ip3+kk*jj*ii-1)
                dst_ifield(1:kk_d, 1:jj_d, 1:ii_d) => wi(ip3_d:ip3_d+kk_d*jj_d*ii_d-1)
                dst_ifield(kstart_d:kstop_d, jstart_d:jstop_d, istart_d:istop_d) = &
                    src_ifield(kstart:kstop, jstart:jstop, istart:istop)
                NULLIFY(dst_ifield)
                NULLIFY(src_ifield)
            END IF
            IF (ASSOCIATED(pi1)) THEN
                src_ifield(1:kk, 1:jj, 1:ii) => pi1(ip3:ip3+kk*jj*ii-1)
                dst_ifield(1:kk_d, 1:jj_d, 1:ii_d) => pi1(ip3_d:ip3_d+kk_d*jj_d*ii_d-1)
                dst_ifield(kstart_d:kstop_d, jstart_d:jstop_d, istart_d:istop_d) = &
                    src_ifield(kstart:kstop, jstart:jstop, istart:istop)
                NULLIFY(dst_ifield)
                NULLIFY(src_ifield)
            END IF
            IF (ASSOCIATED(pi2)) THEN
                src_ifield(1:kk, 1:jj, 1:ii) => pi2(ip3:ip3+kk*jj*ii-1)
                dst_ifield(1:kk_d, 1:jj_d, 1:ii_d) => pi2(ip3_d:ip3_d+kk_d*jj_d*ii_d-1)
                dst_ifield(kstart_d:kstop_d, jstart_d:jstop_d, istart_d:istop_d) = &
                    src_ifield(kstart:kstop, jstart:jstop, istart:istop)
                NULLIFY(dst_ifield)
                NULLIFY(src_ifield)
            END IF
            IF (ASSOCIATED(pi3)) THEN
                src_ifield(1:kk, 1:jj, 1:ii) => pi3(ip3:ip3+kk*jj*ii-1)
                dst_ifield(1:kk_d, 1:jj_d, 1:ii_d) => pi3(ip3_d:ip3_d+kk_d*jj_d*ii_d-1)
                dst_ifield(kstart_d:kstop_d, jstart_d:jstop_d, istart_d:istop_d) = &
                    src_ifield(kstart:kstop, jstart:jstop, istart:istop)
                NULLIFY(dst_ifield)
                NULLIFY(src_ifield)
            END IF
        ELSE
            exU = (sn .AND. ifacerecv < 3) .OR. (.NOT. sn)
            IF (ASSOCIATED(u) .AND. exU) THEN
                src_field(1:kk, 1:jj, 1:ii) => u(ip3:ip3+kk*jj*ii-1)
                dst_field(1:kk_d, 1:jj_d, 1:ii_d) => u(ip3_d:ip3_d+kk_d*jj_d*ii_d-1)
                dst_field(kstart_d:kstop_d, jstart_d:jstop_d, istart_d:istop_d) = &
                    src_field(kstart:kstop, jstart:jstop, istart:istop)
                NULLIFY(dst_field)
                NULLIFY(src_field)
            END IF
            exV = (sn .AND. (ifacerecv > 2 .AND. ifacerecv < 5)) .OR. (.NOT. sn)
            IF (ASSOCIATED(v) .AND. exV) THEN
                src_field(1:kk, 1:jj, 1:ii) => v(ip3:ip3+kk*jj*ii-1)
                dst_field(1:kk_d, 1:jj_d, 1:ii_d) => v(ip3_d:ip3_d+kk_d*jj_d*ii_d-1)
                dst_field(kstart_d:kstop_d, jstart_d:jstop_d, istart_d:istop_d) = &
                    src_field(kstart:kstop, jstart:jstop, istart:istop)
                NULLIFY(dst_field)
                NULLIFY(src_field)
            END IF
            exW = (sn .AND. ifacerecv > 4) .OR. (.NOT. sn)
            IF (ASSOCIATED(w) .AND. exW) THEN
                src_field(1:kk, 1:jj, 1:ii) => w(ip3:ip3+kk*jj*ii-1)
                dst_field(1:kk_d, 1:jj_d, 1:ii_d) => w(ip3_d:ip3_d+kk_d*jj_d*ii_d-1)
                dst_field(kstart_d:kstop_d, jstart_d:jstop_d, istart_d:istop_d) = &
                    src_field(kstart:kstop, jstart:jstop, istart:istop)
                NULLIFY(dst_field)
                NULLIFY(src_field)
            END IF
            IF (ASSOCIATED(p1)) THEN
                src_field(1:kk, 1:jj, 1:ii) => p1(ip3:ip3+kk*jj*ii-1)
                dst_field(1:kk_d, 1:jj_d, 1:ii_d) => p1(ip3_d:ip3_d+kk_d*jj_d*ii_d-1)
                dst_field(kstart_d:kstop_d, jstart_d:jstop_d, istart_d:istop_d) = &
                    src_field(kstart:kstop, jstart:jstop, istart:istop)
                NULLIFY(dst_field)
                NULLIFY(src_field)
            END IF
            IF (ASSOCIATED(p2)) THEN
                src_field(1:kk, 1:jj, 1:ii) => p2(ip3:ip3+kk*jj*ii-1)
                dst_field(1:kk_d, 1:jj_d, 1:ii_d) => p2(ip3_d:ip3_d+kk_d*jj_d*ii_d-1)
                dst_field(kstart_d:kstop_d, jstart_d:jstop_d, istart_d:istop_d) = &
                    src_field(kstart:kstop, jstart:jstop, istart:istop)
                NULLIFY(dst_field)
                NULLIFY(src_field)
            END IF
            IF (ASSOCIATED(p3)) THEN
                src_field(1:kk, 1:jj, 1:ii) => p3(ip3:ip3+kk*jj*ii-1)
                dst_field(1:kk_d, 1:jj_d, 1:ii_d) => p3(ip3_d:ip3_d+kk_d*jj_d*ii_d-1)
                dst_field(kstart_d:kstop_d, jstart_d:jstop_d, istart_d:istop_d) = &
                    src_field(kstart:kstop, jstart:jstop, istart:istop)
                NULLIFY(dst_field)
                NULLIFY(src_field)
            END IF
        END IF
    END SUBROUTINE connect_self


    ! Process receive buffers as they arrive, wait for send
    ! buffers to be free
    SUBROUTINE process_bufs()

        INTEGER(int32) :: idx, i
        TYPE(MPI_Status) :: recvstatus
        INTEGER(int32) :: recvmessagelen
        INTEGER(int32) :: unpacklen

        DO WHILE (.TRUE.)
            CALL MPI_Waitany(nRecv, recvReqs, idx, recvstatus)

            IF (idx /= MPI_UNDEFINED) THEN
                IF (connect_integer) THEN
                    CALL MPI_Get_count(recvstatus, mglet_mpi_ifk, &
                        recvmessagelen)
                ELSE
                    CALL MPI_Get_count(recvstatus, mglet_mpi_real, &
                        recvmessagelen)
                END IF

                unpackLen = 0
                DO i = 1, iRecv
                    IF (recvIdxList(1, i) == recvList(idx) .AND. recvIdxList(2, i) > 0) THEN
                        CALL read_buffer(i)
                        unpacklen = unpacklen + recvIdxList(2, i)
                    END IF
                END DO

                IF (recvmessagelen /= unpacklen) THEN
                    CALL errr(__FILE__, __LINE__)
                END IF
            ELSE
                EXIT
            END IF
        END DO
        CALL MPI_Waitall(nSend, sendReqs, MPI_STATUSES_IGNORE)
    END SUBROUTINE process_bufs


    SUBROUTINE init_connect2()
        INTEGER(intk) :: i, iface, igrid
        INTEGER(intk) :: iface1, iface2, iface3
        INTEGER(intk) :: itypbc1, itypbc2, itypbc3
        INTEGER(intk) :: iprocnbr, itypbc, inbrface, inbrgrid

        INTEGER(int32), ALLOCATABLE :: maxTag(:)
        INTEGER(int32), ALLOCATABLE :: sendcounts(:), sdispls(:)
        INTEGER(int32), ALLOCATABLE :: recvcounts(:), rdispls(:)

        INTEGER(intk) :: nFaceTot
        INTEGER(intk) :: nFaceGeom
        INTEGER(intk) :: nLineTot
        INTEGER(intk) :: nLineGeom
        INTEGER(intk) :: nCornerTot
        INTEGER(intk) :: nCornerGeom

        INTEGER(intk) :: neighbours(26)

        LOGICAL :: exchange
        INTEGER :: iexchange

        nFaceTot = 0
        nFaceGeom = 0
        nLineTot = 0
        nLineGeom = 0
        nCornerTot = 0
        nCornerGeom = 0

        CALL set_timer(150, "CONNECT2")

        ! Maximum number of connections for "simple" cases is number
        ! of grids*26. However, due to the possible prescence of
        ! precursors etc, we add a few more.
        maxConns = INT((nMyGrids+1)*26.0*1.2, intk)
        ALLOCATE(sendConns(8, maxConns))
        ALLOCATE(recvConns(8, maxConns))
        sendConns = 0
        recvConns = 0

        ! The maximum number of concurrent communications are the number
        ! of processes
        ALLOCATE(recvIdxList(3, maxConns))
        ALLOCATE(sendList(numprocs))
        ALLOCATE(recvList(numprocs))
        ALLOCATE(sendReqs(numprocs))
        ALLOCATE(recvReqs(numprocs))
        recvIdxList = 0
        sendList = 0
        recvList = 0

        ALLOCATE(maxTag(0:numprocs-1))
        ALLOCATE(sendcounts(0:numprocs-1))
        ALLOCATE(sdispls(0:numprocs-1))
        ALLOCATE(recvcounts(0:numprocs-1))
        ALLOCATE(rdispls(0:numprocs-1))
        maxTag = 0
        sendcounts = 0
        sdispls = 0
        recvcounts = 0
        rdispls = 0

        ! It is really important that nplane = 2 also in preconnect
        nplane = 2
        nRecv = 0

#if 0
#ifdef _IB_CUTCELL_
        bp = bpbp
#ifndef _IB_CUTCELL_MOVING_
        DO i = 1, nMyGrids
            igrid = myGrids(i)
            CALL get_ip3(ip3, igrid)
            CALL get_mgdims(kk, jj, ii, igrid)
            CALL mgdpcc(ip1cc, ip3cc, ip4cc, ip6cc, ip8cc, ip24cc, igrid)
            CALL create_ccbp(kk, jj, ii, bp(ip3), icellsgrids(igrid), &
                xpsInd(ip3cc), smallP(ip1cc), masterP(ip8cc), nmasterp(igrid), &
                masterSmallP(ip8cc))
        END DO
#endif
#endif
#endif

        DO i = 1, nMyGrids
            igrid = myGrids(i)
            CALL get_neighbours(neighbours, igrid)

            ! Loop over the boundary faces 1..6
            !
            ! Meaning of iterator:
            !   1 : FRONT  ( low X)
            !   2 : BACK   (high X)
            !   3 : RIGHT  ( low Y)
            !   4 : LEFT   (high Y)
            !   5 : BOTTOM ( low Z)
            !   6 : TOP    (high Z)
            !
            ! See also setcobone.src
            DO iface = 1, 6
                ! Get type of BC (assuming CON is ibocond = 1)
                itypbc = itypboconds(1, iface, igrid)

                ! If the process that holds this grid is to send something
                ! to another process
                !
                ! Meaning of 'itypbc':
                !    7: CON - grids connect directly in two directions
                !   19: CO1 - one-way grid coupling (only front do this)
                !
                ! See also setcobone.F
                IF (itypbc == 7 .OR. itypbc == 19) THEN
                    CALL get_nbrs(iface, neighbours, inbrgrid, inbrface)
                    IF (inbrgrid == 0) THEN
                        CYCLE
                    END IF
                    iprocnbr = idprocofgrd(inbrgrid)
                    nFaceTot = nFaceTot + 1

                    ! Check face
                    iexchange = 0
#if 0
#ifdef _IB_
                    CALL check_bp(exchange, igrid, iface, kk, jj, ii, bp(ip3))
                    ! One-sided connects need to be marked as if they contain
                    ! geometry. There is a conenct in mgpoisit_cc:
                    !     CALL connect2(ilevel, 1, v1=GSAW, v2=GSAS, v3=GSAB,
                    !         geom=.TRUE., normal=.TRUE., forward=1)
                    ! where all CO1 faces need to be connected.
                    ! Instead of penalizing the performance in every case
                    ! we choose to flag all CO1-faces as geometry.
                    IF (itypbc == 19) THEN
                        exchange = .TRUE.
                    END IF
#else
                    exchange = .FALSE.
#endif
#else
                    exchange = .TRUE.
#endif
                    IF (exchange) THEN
                        iexchange = 1
                        nFaceGeom = nFaceGeom + 1
                    END IF

                    nRecv = nRecv + 1
                    maxTag(iprocnbr) = maxTag(iprocnbr) + 1

                    recvConns(1, nRecv) = myid      ! Receiving process (this process)
                    recvConns(2, nRecv) = iprocnbr  ! Sending process (neighbour process)
                    recvConns(3, nRecv) = igrid     ! Receiving grid (on current process)
                    recvConns(4, nRecv) = inbrgrid  ! Sending grid (on neighbour process)
                    recvConns(5, nRecv) = iface     ! Which face receive (1..26)
                    recvConns(6, nRecv) = inbrface  ! Which face receive from (sending face) (1..26)
                    recvConns(7, nRecv) = maxTag(iprocnbr)  ! Message tag
                    recvConns(8, nRecv) = iexchange  ! Geometry exchange flag

                    sendcounts(iprocnbr) = sendcounts(iprocnbr) + SIZE(recvConns, 1)
                END IF
            END DO

            ! Check lines
            DO iface = 7, 18
                iface1 = facelist(2, iface)
                iface2 = facelist(3, iface)

                ! Get type of BC
                itypbc1 = itypboconds(1, iface1, igrid)
                itypbc2 = itypboconds(1, iface2, igrid)

                IF (itypbc1 == 7 .OR. itypbc2 == 7 .OR. &
                    itypbc1 == 19 .OR. itypbc2 == 19) THEN
                    CALL get_nbrs(iface, neighbours, inbrgrid, inbrface)
                    IF (inbrgrid == 0) THEN
                        CYCLE
                    END IF
                    iprocnbr = idprocofgrd(inbrgrid)
                    nLineTot = nLineTot + 1
                    iexchange = 0
#if 0
#ifdef _IB_
                    CALL check_bp(exchange, igrid, iface, kk, jj, ii, bp(ip3))
                    IF (itypbc1 == 19 .OR. itypbc2 == 19) THEN
                        exchange = .TRUE.
                    END IF
#else
                    exchange = .FALSE.
#endif
#else
                    exchange = .TRUE.
#endif
                    IF (exchange) THEN
                        iexchange = 1
                        nLineGeom = nLineGeom + 1
                    END IF

                    nRecv = nRecv + 1
                    maxTag(iprocnbr) = maxTag(iprocnbr) + 1

                    recvConns(1, nRecv) = myid      ! Receiving process (this process)
                    recvConns(2, nRecv) = iprocnbr  ! Sending process (neighbour process)
                    recvConns(3, nRecv) = igrid     ! Receiving grid (on current process)
                    recvConns(4, nRecv) = inbrgrid  ! Sending grid (on neighbour process)
                    recvConns(5, nRecv) = iface     ! Which face receive (1..26)
                    recvConns(6, nRecv) = inbrface  ! Which face receive from (sending face) (1..26)
                    recvConns(7, nRecv) = maxTag(iprocnbr)  ! Message tag
                    recvConns(8, nRecv) = iexchange  ! Geometry exchange flag

                    sendcounts(iprocnbr) = sendcounts(iprocnbr) + SIZE(recvConns, 1)
                END IF
            END DO

            ! Check corners
            DO iface = 19, 26
                iface1 = facelist(2, iface)
                iface2 = facelist(3, iface)
                iface3 = facelist(4, iface)

                ! Get type of BC
                itypbc1 = itypboconds(1, iface1, igrid)
                itypbc2 = itypboconds(1, iface2, igrid)
                itypbc3 = itypboconds(1, iface3, igrid)

                IF (itypbc1 == 7 .OR. itypbc2 == 7 .OR. itypbc3 == 7 .OR. &
                    itypbc1 == 19 .OR. itypbc2 == 19 .OR. itypbc3 == 19) THEN
                    CALL get_nbrs(iface, neighbours, inbrgrid, inbrface)
                    IF (inbrgrid == 0) THEN
                        CYCLE
                    END IF
                    iprocnbr = idprocofgrd(inbrgrid)
                    nCornerTot = nCornerTot + 1
                    iexchange = 0
#if 0
#ifdef _IB_
                    CALL check_bp(exchange, igrid, iface, kk, jj, ii, bp(ip3))
                    IF (itypbc1 == 19 .OR. itypbc2 == 19 .OR. &
                        itypbc3 == 19) THEN
                        exchange = .TRUE.
                    END IF
#else
                    exchange = .FALSE.
#endif
#else
                    exchange = .TRUE.
#endif
                    IF (exchange) THEN
                        iexchange = 1
                        nCornerGeom = nCornerGeom + 1
                    END IF

                    nRecv = nRecv + 1
                    maxTag(iprocnbr) = maxTag(iprocnbr) + 1

                    recvConns(1, nRecv) = myid      ! Receiving process (this process)
                    recvConns(2, nRecv) = iprocnbr  ! Sending process (neighbour process)
                    recvConns(3, nRecv) = igrid     ! Receiving grid (on current process)
                    recvConns(4, nRecv) = inbrgrid  ! Sending grid (on neighbour process)
                    recvConns(5, nRecv) = iface     ! Which face receive (1..26)
                    recvConns(6, nRecv) = inbrface  ! Which face receive from (sending face) (1..26)
                    recvConns(7, nRecv) = maxTag(iprocnbr)  ! Message tag
                    recvConns(8, nRecv) = iexchange  ! Geometry exchange flag

                    sendcounts(iprocnbr) = sendcounts(iprocnbr) + SIZE(recvConns, 1)
                END IF
            END DO
        END DO

        iRecv = nRecv

        ! Sort recvConns by process ID
        CALL sort_conns(recvConns(:,1:nRecv))

        ! Calculate sdispl offset
        DO i=1,numprocs-1
            sdispls(i) = sdispls(i-1) + sendcounts(i-1)
        END DO

        ! First exchange NUMBER OF ELEMENTS TO RECEIVE, to be able to
        ! calculate rdispls array
        CALL MPI_Alltoall(sendcounts, 1, MPI_INTEGER, recvcounts, 1, &
            MPI_INTEGER, MPI_COMM_WORLD)

        ! Calculate rdispl offset
        DO i=1,numprocs-1
            rdispls(i) = rdispls(i-1) + recvcounts(i-1)
        END DO

        ! Check that number of connections fit in array
        iSend = (rdispls(numprocs-1) + recvcounts(numprocs-1))/SIZE(sendConns, 1)
        IF (iSend > maxConns) THEN
            write(*,*) "Number of connections exceeded on process ", myid
            write(*,*) "maxConns =", maxConns, "nMyGrids =", nMyGrids, &
                "iSend = ", iSend
            CALL errr(__FILE__, __LINE__)
        END IF

        ! Exchange connection information
        CALL MPI_Alltoallv(recvConns(1, 1), sendcounts, sdispls, MPI_INTEGER, &
            sendConns(1, 1), recvcounts, rdispls, MPI_INTEGER, &
            MPI_COMM_WORLD)

        isInit = .TRUE.

        ! Agglomerate statistics
        CALL MPI_Allreduce(MPI_IN_PLACE, nFaceTot, 1, mglet_mpi_int, &
            MPI_SUM, MPI_COMM_WORLD)
        CALL MPI_Allreduce(MPI_IN_PLACE, nFaceGeom, 1, mglet_mpi_int, &
            MPI_SUM, MPI_COMM_WORLD)
        CALL MPI_Allreduce(MPI_IN_PLACE, nLineTot, 1, mglet_mpi_int, &
            MPI_SUM, MPI_COMM_WORLD)
        CALL MPI_Allreduce(MPI_IN_PLACE, nLineGeom, 1, mglet_mpi_int, &
            MPI_SUM, MPI_COMM_WORLD)
        CALL MPI_Allreduce(MPI_IN_PLACE, nCornerTot, 1, mglet_mpi_int, &
            MPI_SUM, MPI_COMM_WORLD)
        CALL MPI_Allreduce(MPI_IN_PLACE, nCornerGeom, 1, mglet_mpi_int, &
            MPI_SUM, MPI_COMM_WORLD)

        IF (myid == 0) THEN
            WRITE(*, '("CONNECT STATISTICS:")')
            WRITE(*, '(4X, "Faces:         ", I7, 4X, "with geometry: ", I7)') &
                nFaceTot, nFaceGeom
            WRITE(*, '(4X, "Lines:         ", I7, 4X, "with geometry: ", I7)') &
                nLineTot, nLineGeom
            WRITE(*, '(4X, "Corners:       ", I7, 4X, "with geometry: ", I7)') &
                nCornerTot, nCornerGeom
            WRITE(*, '()')
        END IF

        nRecv = 0

        ! Nullify pointers
        nullify(u)
        nullify(v)
        nullify(w)
        nullify(p1)
        nullify(p2)
        nullify(p3)
    END SUBROUTINE init_connect2


    SUBROUTINE get_nbrs(iface, neighbours, nbrgrid, nbrface)
        INTEGER(intk), INTENT(IN) :: iface
        INTEGER(intk), INTENT(IN) :: neighbours(26)
        INTEGER(intk), INTENT(OUT) :: nbrgrid
        INTEGER(intk), INTENT(OUT) :: nbrface

        INTEGER(intk) :: n_rescue, i, dir
        INTEGER(intk) :: iface1, iface2, iface3
        INTEGER(intk) :: itypbc1, itypbc2, itypbc3

        ! Should be 7...
        n_rescue = SIZE(rescue_nbr, 1)

        ! 0 means no connect
        nbrgrid = 0
        nbrface = 0
        DO i = 1, n_rescue
            dir = rescue_dir(i, iface)

            ! rescue_dir is ordered and when a 0 is encountered there is
            ! nothing more to do...
            IF (dir == 0) THEN
                EXIT
            END IF

            ! If there is a neighbour in this position, use this
            IF (neighbours(dir) > 0) THEN
                nbrgrid = neighbours(dir)
                nbrface = rescue_nbr(i, iface)

                ! Check if this is suited for a connect (symmetry req.)
                ! This require knowledge of the global grid structure -
                ! currently this is OK.
                IF (nbrface > 18) THEN
                    ! Get adjacent primary faces
                    iface1 = facelist(2, nbrface)
                    iface2 = facelist(3, nbrface)
                    iface3 = facelist(4, nbrface)

                    ! Get type of BC on these
                    itypbc1 = itypboconds(1, iface1, nbrgrid)
                    itypbc2 = itypboconds(1, iface2, nbrgrid)
                    itypbc3 = itypboconds(1, iface3, nbrgrid)

                    ! If none of the neighboring faces are CON or CO1, the connect
                    ! should not be carried out - check next neighbour
                    IF ((.NOT. (itypbc1 == 7 .OR. itypbc1 == 19)) .AND. &
                            (.NOT. (itypbc2 == 7 .OR. itypbc2 == 19)) .AND. &
                            (.NOT. (itypbc3 == 7 .OR. itypbc3 == 19))) THEN
                        ! Reset neighbour information and cycle loop
                        nbrgrid = 0
                        nbrface = 0
                        CYCLE
                    END IF
                END IF

                ! If sofar, all is good!
                EXIT
            END IF
        END DO
    END SUBROUTINE get_nbrs


    SUBROUTINE finish_connect2()
        isInit = .FALSE.

        DEALLOCATE(sendConns)
        DEALLOCATE(recvConns)

        DEALLOCATE(recvIdxList)
        DEALLOCATE(sendList)
        DEALLOCATE(recvList)
        DEALLOCATE(sendReqs)
        DEALLOCATE(recvReqs)
    END SUBROUTINE finish_connect2


    SUBROUTINE sort_conns(list)
        ! Input array to be sorted
        INTEGER(int32), INTENT(inout) :: list(:,:)

        INTEGER(intk) :: i,j

        ! Temporary storage
        INTEGER(int32) :: temp(8)

        IF (SIZE(list, 1) /= SIZE(temp)) THEN
            CALL errr(__FILE__, __LINE__)
        END IF

        ! Sort by sending processor number (field 2)
        DO i = 2, SIZE(list, 2)
            j = i - 1
            temp(:) = list(:,i)
            DO WHILE (j >= 1)
                IF (list(2,j) > temp(2)) THEN
                    list(:,j+1) = list(:,j)
                    j = j - 1
                ELSE
                    EXIT
                END IF
            END DO
            list(:,j+1) = temp(:)
        END DO

    END SUBROUTINE sort_conns


    SUBROUTINE start_and_stop(igrid, iface, istart, &
        istop, jstart, jstop, kstart, kstop, nghost)

        ! Subroutine arguments
        INTEGER(intk), INTENT(in) :: igrid, iface
        INTEGER(intk), INTENT(out) :: istart, istop, jstart, jstop, &
            kstart, kstop
        INTEGER(intk), INTENT(in), OPTIONAL :: nghost

        ! Local variables
        INTEGER(intk) :: kk, jj, ii, i

        ! Get diensions of grid
        CALL get_mgdims(kk, jj, ii, igrid)

        istart = 3
        istop = ii - 2

        jstart = 3
        jstop = jj - 2

        kstart = 3
        kstop = kk - 2

        DO i = 1, facelist(1, iface)
            CALL start_and_stop_face(igrid, facelist(i + 1, iface), &
                istart, istop, jstart, jstop, kstart, kstop, nghost)
        END DO
    END SUBROUTINE start_and_stop


    SUBROUTINE start_and_stop_face(igrid, iface, istart, &
        istop, jstart, jstop, kstart, kstop, nghost)

        ! Subroutine arguments
        INTEGER(intk), INTENT(in) :: igrid, iface
        INTEGER(intk), INTENT(inout) :: istart, istop, jstart, jstop, &
            kstart, kstop
        INTEGER(intk), INTENT(in), OPTIONAL :: nghost

        ! Internal variables
        INTEGER(intk) :: kk, jj, ii

        ! Get diensions of grid
        CALL get_mgdims(kk, jj, ii, igrid)

        ! Sanity check
        IF (nplane > 2 .OR. nplane < 1) THEN
            CALL errr(__FILE__, __LINE__)
        END IF

        SELECT CASE(iface)
            ! Front
            CASE (1)
                istart = 3 - nplane
                istop = 2
            !   Back
            CASE (2)
                istart = ii - 1
                istop = ii - (2 - nplane)
            !   Right
            CASE (3)
                jstart = 3 - nplane
                jstop = 2
            !   Left
            CASE (4)
                jstart = jj - 1
                jstop = jj - (2 - nplane)
            !   Bottom
            CASE (5)
                kstart = 3 - nplane
                kstop = 2
            !   Top
            CASE (6)
                kstart = kk - 1
                kstop = kk - (2 - nplane)
            CASE DEFAULT
                CALL errr(__FILE__, __LINE__)
        END SELECT

        ! Indices for the non-ghost region, i.e. the region inside
        ! the computational domain
        IF ((PRESENT(nghost) .AND. flag /= 'Y') .OR. &
                (.NOT. PRESENT(nghost) .AND. flag == 'Y')) THEN
            SELECT CASE(iface)
                ! Front
                CASE (1)
                    istart = istart + nplane
                    istop = istop + nplane
                !   Back
                CASE (2)
                    istart = istart - nplane
                    istop = istop - nplane
                !   Right
                CASE (3)
                    jstart = jstart + nplane
                    jstop = jstop + nplane
                !   Left
                CASE (4)
                    jstart = jstart - nplane
                    jstop = jstop - nplane
                !   Bottom
                CASE (5)
                    kstart = kstart + nplane
                    kstop = kstop + nplane
                !   Top
                CASE (6)
                    kstart = kstart - nplane
                    kstop = kstop - nplane
                CASE DEFAULT
                    CALL errr(__FILE__, __LINE__)
            END SELECT
        END IF
    END SUBROUTINE start_and_stop_face


    ! Calculate the area of a boundary face, i.e. the number of cells
    ! to be exchanged (all planes). Used to calculate message lengths.
    FUNCTION face_area(igrid, iface, nghost) RESULT(area)

        ! Result = length of message to be passed
        INTEGER(intk) :: area

        ! Input parameters, grid and boundary face information
        ! Must be intk because they are intreface to other MGLET functions
        INTEGER(intk), INTENT(in) :: igrid, iface

        ! Indicate that we want the non-ghost region of the grid,
        ! i.e. the part in the computational domain that correspond to the
        ! neighbour ghsot region.
        !
        ! It does not matter which value it has, it's prescence act as a
        ! switch.
        INTEGER(intk), INTENT(in), OPTIONAL :: nghost

        ! Indices of start- and stop of iteration over boundary face
        ! Must be intk because they are intreface to other MGLET functions
        INTEGER(intk) :: istart, istop, jstart, jstop, kstart, kstop

        CALL start_and_stop(igrid, iface, istart, istop, &
            jstart, jstop, kstart, kstop, nghost)

        area = (istop-istart+1)*(jstop-jstart+1)*(kstop-kstart+1)

        RETURN
    END FUNCTION face_area
END MODULE connect2_mod
