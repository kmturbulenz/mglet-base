MODULE parent_core_mod
    USE core_mod
    USE MPI_f08

    IMPLICIT NONE
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
    INTEGER(intk), ALLOCATABLE, PROTECTED :: sendconns(:, :), recvconns(:, :)

    ! Number of send and receive connections
    INTEGER(intk), PROTECTED :: isend = 0, irecv = 0

    LOGICAL :: is_parent_core_init = .FALSE.

    ! Public read-only variables
    PUBLIC :: sendconns, recvconns, isend, irecv, is_parent_core_init

    ! Public subroutines
    PUBLIC :: init_parent_core, finish_parent_core, idx2d, stag, &
        start_and_stop, face_area
CONTAINS
    SUBROUTINE init_parent_core()
        ! Subroutine arguments
        ! none...

        ! Local variables
        INTEGER(intk) :: i, iface, igrid, inbr, iprocnbr, itypbc, maxconns, &
            nrecv
        INTEGER(int32), ALLOCATABLE :: sendcounts(:), sdispls(:)
        INTEGER(int32), ALLOCATABLE :: recvcounts(:), rdispls(:)
        INTEGER(int32) :: ierr

        IF (is_parent_core_init) CALL errr(__FILE__, __LINE__)

        CALL set_timer(210, "PARENT")

        ! Maximum number of parents for "simple" cases is number
        ! of grids*6. However, due to the possible prescence of
        ! strange grid structures we add a few more
        maxconns = INT((nmygrids+1.0)*6.0*1.2)
        ALLOCATE(recvconns(5, maxconns), source=0_intk)

        ALLOCATE(sendcounts(0:numprocs-1), source=0_int32)
        ALLOCATE(sdispls(0:numprocs-1), source=0_int32)
        ALLOCATE(recvcounts(0:numprocs-1), source=0_int32)
        ALLOCATE(rdispls(0:numprocs-1), source=0_int32)

        nrecv = 0

        DO i = 1, nmygrids
            igrid = myGrids(i)

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
            ! See also setcobone.F
            DO iface = 1, 6
                ! Get type of BC (assuming PAR is ibocond = 1)
                itypbc = itypboconds(1, iface, igrid)

                ! See also setcobone.F
                IF (itypbc == 8) THEN
                    ! Get neighbouring grid and process

                    inbr = iparent(igrid)
                    iprocnbr = idprocofgrd(inbr)

                    nrecv = nrecv + 1

                    IF (nrecv > maxconns) THEN
                        write(*, *) "Number of PAR's exceeded on process ", myid
                        write(*, *) "maxconns =", maxconns, &
                            "nmygrids =", nmygrids, "nrecv = ", nrecv
                        CALL errr(__FILE__, __LINE__)
                    END IF

                    ! Receiving process (this process)
                    recvconns(1, nrecv) = myid

                    ! Sending process (neighbour process)
                    recvconns(2, nrecv) = iprocnbr

                    ! Receiving grid (on current process)
                    recvconns(3, nrecv) = igrid

                    ! Sending grid (on neighbour process)
                    recvconns(4, nrecv) = inbr

                    ! Which face receive (1..6)
                    recvconns(5, nrecv) = iface

                    sendcounts(iprocnbr) = sendcounts(iprocnbr) &
                        + SIZE(recvconns, 1)
                END IF
            END DO
        END DO

        irecv = nrecv

        ! Sort recvconns by process ID
        CALL sort_conns(recvconns(:, 1:nrecv), 2)

        ! Calculate sdispl offset
        DO i = 1, numprocs-1
            sdispls(i) = sdispls(i-1) + sendcounts(i-1)
        END DO

        ! First exchange NUMBER OF ELEMENTS TO RECEIVE, to be able to
        ! calculate rdispls array
        CALL MPI_Alltoall(sendcounts, 1, MPI_INTEGER, recvcounts, 1, &
            MPI_INTEGER, MPI_COMM_WORLD, ierr)

        ! Calculate rdispl offset
        DO i=1, numprocs-1
            rdispls(i) = rdispls(i-1) + recvcounts(i-1)
        END DO

        ! Allocate sendconns array
        isend = (rdispls(numprocs-1) &
            + recvcounts(numprocs-1))/SIZE(recvconns, 1)
        ALLOCATE(sendconns(5, isend), source=0_intk)

        ! Exchange connection information
        CALL MPI_Alltoallv(recvconns, sendcounts, sdispls, MPI_INTEGER, &
            sendconns, recvcounts, rdispls, MPI_INTEGER, &
            MPI_COMM_WORLD)

        is_parent_core_init = .TRUE.
    END SUBROUTINE init_parent_core


    SUBROUTINE idx2d(kk, jj, ii, iface, kkc, jjc, iic, jj2d, ii2d, jjc2d, iic2d)
        ! Subroutine arguments
        INTEGER(intk), INTENT(in) :: kk, jj, ii, iface
        INTEGER(intk), INTENT(out) :: kkc, jjc, iic
        INTEGER(intk), INTENT(out) :: jj2d, ii2d, jjc2d, iic2d

        kkc = kk/2 + 2
        jjc = jj/2 + 2
        iic = ii/2 + 2

        IF (iface == 1 .OR. iface == 2) THEN
            ii2d = jj
            jj2d = kk
            iic2d = jjc
            jjc2d = kkc
        ELSE IF (iface == 3 .OR. iface == 4) THEN
            ii2d = ii
            jj2d = kk
            iic2d = iic
            jjc2d = kkc
        ELSE IF (iface == 5 .OR. iface == 6) THEN
            ii2d = ii
            jj2d = jj
            iic2d = iic
            jjc2d = jjc
        ELSE
            CALL errr(__FILE__, __LINE__)
        END IF
    END SUBROUTINE idx2d


    SUBROUTINE stag(iface, ustag1, ustag2, vstag1, vstag2, wstag1, wstag2)
        ! Subroutine arguments
        INTEGER(intk), INTENT(in) :: iface
        INTEGER(intk), INTENT(out) :: ustag1, ustag2, vstag1, vstag2, &
            wstag1, wstag2

        ustag1 = 0
        ustag2 = 0
        vstag1 = 0
        vstag2 = 0
        wstag1 = 0
        wstag2 = 0

        IF (iface == 1 .OR. iface == 2) THEN
            vstag2 = 1
            wstag1 = 1
        ELSE IF (iface == 3 .OR. iface == 4) THEN
            ustag2 = 1
            wstag1 = 1
        ELSE IF (iface == 5 .OR. iface == 6) THEN
            ustag2 = 1
            vstag1 = 1
        ELSE
            CALL errr(__FILE__, __LINE__)
        END IF
    END SUBROUTINE stag


    SUBROUTINE start_and_stop(igrid, iface, icomp, &
        istart, istop, jstart, jstop, kstart, kstop)

        ! Subroutine arguments
        ! Igrid and iface of *fine* grid is to be given
        INTEGER(intk), INTENT(in) :: igrid, iface, icomp

        ! Returns start- and stop indices of the corresponding *coarse* grid
        ! region to be copied
        INTEGER(intk), INTENT(out) :: istart, istop, jstart, jstop, &
            kstart, kstop

        ! Local variables
        INTEGER(intk) :: kk, jj, ii
        INTEGER(intk) :: kpos, jpos, ipos

        ! Get diensions of (fine)grid
        CALL get_mgdims(kk, jj, ii, igrid)
        kpos = kposition(igrid)
        jpos = jposition(igrid)
        ipos = iposition(igrid)

        ! Select entire slab in coarse grid
        kstart = kpos - 2
        jstart = jpos - 2
        istart = ipos - 2

        kstop = kstart + kk/2 + 1
        jstop = jstart + jj/2 + 1
        istop = istart + ii/2 + 1

        ! Reduce selection to only contain face with 1 or 2 slices
        IF (iface == 1) THEN
            istart = ipos - 1
            istop = ipos - 1
        ELSE IF (iface == 2) THEN
            istart = ipos + (ii-4)/2
            IF (icomp == 1) istart = istart - 1
            istop = istart
        ELSE IF (iface == 3) THEN
            jstart = jpos - 1
            jstop = jpos - 1
        ELSE IF (iface == 4) THEN
            jstart = jpos + (jj-4)/2
            IF (icomp == 2) jstart = jstart - 1
            jstop = jstart
        ELSE IF (iface == 5) THEN
            kstart = kpos - 1
            kstop = kpos - 1
        ELSE IF (iface == 6) THEN
            kstart = kpos + (kk-4)/2
            IF (icomp == 3) kstart = kstart - 1
            kstop = kstart
        END IF
    END SUBROUTINE start_and_stop


    ! Calculate the area of a boundary face, i.e. the number of cells
    ! to be exchanged (all planes). Used to calculate message lengths.
    ! Needs to be given *fine grid* igrid and iface
    FUNCTION face_area(igrid, iface) RESULT(area)

        ! Result = length of message to be passed
        INTEGER(intk) :: area

        ! Input parameters, grid and boundary face information
        ! Must be intk because they are intreface to other MGLET functions
        INTEGER(intk), INTENT(in) :: igrid, iface

        ! Indices of start- and stop of iteration over boundary face
        ! Must be intk because they are intreface to other MGLET functions
        INTEGER(intk) :: istart, istop, jstart, jstop, kstart, kstop

        CALL start_and_stop(igrid, iface, 0, istart, istop, &
            jstart, jstop, kstart, kstop)

        area = (istop-istart+1)*(jstop-jstart+1)*(kstop-kstart+1)

        RETURN
    END FUNCTION face_area


    SUBROUTINE finish_parent_core()
        ! Subroutine arguments
        ! none...

        ! Local variables
        ! none...

        IF (.NOT. is_parent_core_init) CALL errr(__FILE__, __LINE__)

        DEALLOCATE(sendconns)
        DEALLOCATE(recvconns)
        isend = 0
        irecv = 0
        is_parent_core_init = .FALSE.
    END SUBROUTINE finish_parent_core
END MODULE parent_core_mod
