MODULE connect_core_mod
    USE MPI_f08
    USE precision_mod
    USE err_mod, ONLY: errr
    USE timer_mod, ONLY: set_timer
    USE grids_mod, ONLY: mygrids, nmygrids, level, idprocofgrd, itypboconds, &
        maxlevel, minlevel, get_neighbours, get_mgdims
    USE comms_mod, ONLY: myid, numprocs
    USE qsort_mod, ONLY: sort_conns

    IMPLICIT NONE (type, external)
    PRIVATE

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
    INTEGER(intk), ALLOCATABLE, PROTECTED :: sendconns(:, :), recvconns(:, :)

    ! Number of send and receive connections
    ! (not used from conn1_mod so not declare target)
    INTEGER(intk), PROTECTED :: isend = 0, irecv = 0

    INTEGER(intk), PARAMETER :: facelist(4, 26) = RESHAPE([ &
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
        3, 2, 4, 6], SHAPE(facelist))

    INTEGER(intk), PARAMETER :: facenbr(26) = [ &
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
        19]

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
    INTEGER(intk), PARAMETER :: rescue_dir(7, 26) = RESHAPE([ &
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
        26, 12, 14, 18, 2, 4, 6], SHAPE(rescue_dir))

    INTEGER(intk), PARAMETER :: rescue_nbr(7, 26) = RESHAPE([ &
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
        19, 20, 21, 23, 22, 24, 25], SHAPE(rescue_nbr))

    PUBLIC :: init_connect_core, finish_connect_core, sendconns, recvconns, &
        isend, irecv, start_and_stop, start_and_stop_face, face_area, &
        facelist, facenbr, corr_start_stop, decide

CONTAINS


    SUBROUTINE init_connect_core()
        ! Local variables
        INTEGER(intk) :: i, iface, igrid
        INTEGER(intk) :: iface1, iface2, iface3
        INTEGER(intk) :: itypbc1, itypbc2, itypbc3
        INTEGER(intk) :: iprocnbr, itypbc, inbrface, inbrgrid
        INTEGER(intk) :: nrecv, maxconns

        INTEGER(int32), ALLOCATABLE :: maxtag(:)
        INTEGER(int32), ALLOCATABLE :: sendcounts(:), sdispls(:)
        INTEGER(int32), ALLOCATABLE :: recvcounts(:), rdispls(:)
        INTEGER(intk), ALLOCATABLE :: recvconns2(:, :)
        INTEGER(intk) :: neighbours(26)
        INTEGER :: iexchange

        CALL set_timer(150, "CONNECT")

        ! Maximum number of connections for "simple" cases is number
        ! of grids*26. However, due to the possible prescence of
        ! precursors etc, we add a few more.
        maxconns = INT((nmygrids+1)*26.0*1.2, intk)
        ALLOCATE(recvconns2(8, maxconns))
        recvconns2 = 0

        ALLOCATE(maxtag(0:numprocs-1))
        ALLOCATE(sendcounts(0:numprocs-1))
        ALLOCATE(sdispls(0:numprocs-1))
        ALLOCATE(recvcounts(0:numprocs-1))
        ALLOCATE(rdispls(0:numprocs-1))
        maxtag = 0
        sendcounts = 0
        sdispls = 0
        recvcounts = 0
        rdispls = 0

        nrecv = 0

        DO i = 1, nmygrids
            igrid = mygrids(i)
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
                    iexchange = 1

                    nrecv = nrecv + 1
                    maxtag(iprocnbr) = maxtag(iprocnbr) + 1

                    recvconns2(1, nrecv) = myid      ! Receiving process (this process)
                    recvconns2(2, nrecv) = iprocnbr  ! Sending process (neighbour process)
                    recvconns2(3, nrecv) = igrid     ! Receiving grid (on current process)
                    recvconns2(4, nrecv) = inbrgrid  ! Sending grid (on neighbour process)
                    recvconns2(5, nrecv) = iface     ! Which face receive (1..26)
                    recvconns2(6, nrecv) = inbrface  ! Which face receive from (sending face) (1..26)
                    recvconns2(7, nrecv) = maxtag(iprocnbr)  ! Message tag
                    recvconns2(8, nrecv) = iexchange  ! Geometry exchange flag

                    sendcounts(iprocnbr) = sendcounts(iprocnbr) &
                        + SIZE(recvconns2, 1)
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
                    iexchange = 1

                    nrecv = nrecv + 1
                    maxtag(iprocnbr) = maxtag(iprocnbr) + 1

                    recvconns2(1, nrecv) = myid      ! Receiving process (this process)
                    recvconns2(2, nrecv) = iprocnbr  ! Sending process (neighbour process)
                    recvconns2(3, nrecv) = igrid     ! Receiving grid (on current process)
                    recvconns2(4, nrecv) = inbrgrid  ! Sending grid (on neighbour process)
                    recvconns2(5, nrecv) = iface     ! Which face receive (1..26)
                    recvconns2(6, nrecv) = inbrface  ! Which face receive from (sending face) (1..26)
                    recvconns2(7, nrecv) = maxtag(iprocnbr)  ! Message tag
                    recvconns2(8, nrecv) = iexchange  ! Geometry exchange flag

                    sendcounts(iprocnbr) = sendcounts(iprocnbr) &
                        + SIZE(recvconns2, 1)
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
                    iexchange = 1

                    nrecv = nrecv + 1
                    maxtag(iprocnbr) = maxtag(iprocnbr) + 1

                    recvconns2(1, nrecv) = myid      ! Receiving process (this process)
                    recvconns2(2, nrecv) = iprocnbr  ! Sending process (neighbour process)
                    recvconns2(3, nrecv) = igrid     ! Receiving grid (on current process)
                    recvconns2(4, nrecv) = inbrgrid  ! Sending grid (on neighbour process)
                    recvconns2(5, nrecv) = iface     ! Which face receive (1..26)
                    recvconns2(6, nrecv) = inbrface  ! Which face receive from (sending face) (1..26)
                    recvconns2(7, nrecv) = maxtag(iprocnbr)  ! Message tag
                    recvconns2(8, nrecv) = iexchange  ! Geometry exchange flag

                    sendcounts(iprocnbr) = sendcounts(iprocnbr) &
                        + SIZE(recvconns2, 1)
                END IF
            END DO
        END DO

        irecv = nrecv

        ! Sort recvconns by process ID
        CALL sort_conns(recvconns2(:, 1:nrecv), 2)

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

        ! Allocate sendconns to the correct size
        isend = (rdispls(numprocs-1) + recvcounts(numprocs-1)) &
            /SIZE(recvconns2, 1)
        ALLOCATE(sendconns(SIZE(recvconns2, 1), isend), SOURCE=0)

        ! Exchange connection information
        CALL MPI_Alltoallv(recvconns2, sendcounts, sdispls, MPI_INTEGER, &
            sendconns, recvcounts, rdispls, MPI_INTEGER, &
            MPI_COMM_WORLD)

        ! Allocate the recvconns array to the correct size and copy the data
        ! from recvconns2
        ALLOCATE(recvconns(8, irecv), SOURCE=recvconns2(:, 1:irecv))
        DEALLOCATE(recvconns2)

        !$omp target enter data map(always, to: sendconns, recvconns)
    END SUBROUTINE init_connect_core


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

                    ! If none of the neighboring faces are CON or CO1, the
                    ! connect should not be carried out - check next neighbour
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


    SUBROUTINE finish_connect_core()
        !$omp target exit data map(always, delete: sendconns, recvconns)

        isend = 0
        irecv = 0
        DEALLOCATE(sendconns)
        DEALLOCATE(recvconns)
    END SUBROUTINE finish_connect_core


    SUBROUTINE start_and_stop(igrid, iface, istart, istop, jstart, jstop, &
            kstart, kstop, nplane, flag, nghost)

        ! Subroutine arguments
        INTEGER(intk), INTENT(in) :: igrid, iface
        INTEGER(intk), INTENT(out) :: istart, istop, jstart, jstop, &
            kstart, kstop
        INTEGER(intk), INTENT(in) :: nplane
        CHARACTER(len=*), INTENT(in) :: flag
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
                istart, istop, jstart, jstop, kstart, kstop, nplane, flag, &
                nghost)
        END DO
    END SUBROUTINE start_and_stop


    SUBROUTINE start_and_stop_face(igrid, iface, istart, &
        istop, jstart, jstop, kstart, kstop, nplane, flag, nghost)

        ! Subroutine arguments
        INTEGER(intk), INTENT(in) :: igrid, iface
        INTEGER(intk), INTENT(inout) :: istart, istop, jstart, jstop, &
            kstart, kstop
        INTEGER(intk), INTENT(in) :: nplane
        CHARACTER(len=*), INTENT(in) :: flag
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


    SUBROUTINE corr_start_stop(igrid, ifacesend, ifacerecv, &
        istart, istop, jstart, jstop, kstart, kstop, nplane, flag)

        ! Subroutine arguments
        INTEGER(intk), INTENT(in) :: igrid, ifacesend, ifacerecv
        INTEGER(intk), INTENT(inout) :: istart, istop, jstart, jstop, &
            kstart, kstop
        INTEGER(intk), INTENT(in) :: nplane
        CHARACTER(len=*), INTENT(in) :: flag

        ! Local variables
        INTEGER(intk) :: i, j
        INTEGER(intk) :: iface

        IF (ifacerecv <= 6) RETURN

        ! Check if send- and recv- face have elementary faces in common
        DO i = 1, facelist(1, ifacesend)
            DO j = 1, facelist(1, ifacerecv)
                IF (facelist(i + 1, ifacesend) &
                        == facelist(j + 1, ifacerecv)) THEN
                    iface = facelist(i + 1, ifacesend)
                    CALL start_and_stop_face(igrid, iface, istart, istop, &
                        jstart, jstop, kstart, kstop, nplane, flag)
                END IF
            END DO
        END DO
    END SUBROUTINE corr_start_stop


    ! Calculate the area of a boundary face, i.e. the number of cells
    ! to be exchanged (all planes). Used to calculate message lengths.
    FUNCTION face_area(igrid, iface, nplane, flag, nghost) RESULT(area)

        ! Result = length of message to be passed
        INTEGER(intk) :: area

        ! Input parameters, grid and boundary face information
        ! Must be intk because they are intreface to other MGLET functions
        INTEGER(intk), INTENT(in) :: igrid, iface

        INTEGER(intk), INTENT(in) :: nplane
        CHARACTER(len=*), INTENT(in) :: flag

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
            jstart, jstop, kstart, kstop, nplane, flag, nghost)

        area = (istop-istart+1)*(jstop-jstart+1)*(kstop-kstart+1)

        RETURN
    END FUNCTION face_area


    FUNCTION decide(i, list, geometry, vertices, fwd, minconlvl, maxconlvl) &
            RESULT(exchange)
        INTEGER(intk), INTENT(in) :: i
        INTEGER(intk), INTENT(in) :: list(:, :)
        LOGICAL, INTENT(in) :: geometry, vertices
        INTEGER(intk), INTENT(in) :: fwd, minconlvl, maxconlvl

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
        IF (ilevel > maxconlvl .OR. ilevel < minconlvl) THEN
            exchange = .FALSE.
        END IF

        ! Only exchange geometry
        IF (geometry .AND. list(8, i) == 0) THEN
            exchange = .FALSE.
        END IF

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
    END FUNCTION decide

END MODULE connect_core_mod
