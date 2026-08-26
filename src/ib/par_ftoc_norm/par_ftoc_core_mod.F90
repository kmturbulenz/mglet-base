MODULE par_ftoc_core_mod
    USE MPI_f08
    USE core_mod

    IMPLICIT NONE (type, external)
    PRIVATE

    ! The information in the first dimension is sorted as follows:
    !   Field 1: Rank of sending process
    !   Field 2: Rank of receiving process
    !   Field 3: ID of sending grid
    !   Field 4: ID of receiving grid
    !   Field 5: Which face (1..6) to send
    !   Field 6: Which face (1..6) to receive, -1 for internal
    INTEGER(intk), ALLOCATABLE, PROTECTED :: sendconns(:, :), recvconns(:, :)

    ! Number of send and receive connections
    INTEGER(intk), PROTECTED :: isend = 0, irecv = 0

    LOGICAL :: is_par_ftoc_core_init = .FALSE.

    PUBLIC :: sendconns, recvconns, isend, irecv, &
        is_par_ftoc_core_init, init_par_ftoc_core, &
        finish_par_ftoc_core, face_area, start_and_stop
CONTAINS
    SUBROUTINE init_par_ftoc_core()
        ! Subroutine arguments
        ! none...

        ! Local variables
        LOGICAL :: one_connect
        INTEGER(intk) :: i, iface, ifacerecv, igrid, inbr, iprocnbr, pos
        INTEGER(intk) :: maxconns, nsend
        INTEGER(intk) :: neighbours(26)
        INTEGER(int32) :: ncols
        CHARACTER(len=8) :: ctyp
        INTEGER(int32), ALLOCATABLE :: sendcounts(:), sdispls(:)
        INTEGER(int32), ALLOCATABLE :: recvcounts(:), rdispls(:)

        IF (is_par_ftoc_core_init) CALL errr(__FILE__, __LINE__)

        CALL set_timer(211, "PAR_FTOC")

        ncols = 6
        maxconns = INT((nmygrids+1.0)*6.0*1.2)
        ALLOCATE(sendconns(ncols, maxconns))
        sendconns = 0

        ALLOCATE(sendcounts(0:numprocs-1))
        ALLOCATE(sdispls(0:numprocs-1))
        ALLOCATE(recvcounts(0:numprocs-1))
        ALLOCATE(rdispls(0:numprocs-1))
        sendcounts = 0
        sdispls = 0
        recvcounts = 0
        rdispls = 0

        nsend = 0

        ! On the coarsest level there are not allowed to be any PAR BC's.
        DO i = 1, nmygrids
            igrid = mygrids(i)

            ! Loop over the boundary faces 1..6
            !
            ! Meaning of iterator:
            !   1 : FRONT  ( low X)
            !   2 : BACK   (high X)
            !   3 : RIGHT  ( low Y)
            !   4 : LEFT   (high Y)
            !   5 : BOTTOM ( low Z)
            !   6 : TOP    (high Z)
            DO iface = 1, 6
                ! We assume the PAR is always first/primary BC for that face
                CALL get_bc_ctyp(ctyp, 1, iface, igrid)
                IF (ctyp(1:3) /= "PAR") CYCLE

                ! The PAR BC needs always to transfer data to the parent grid,
                ! independent on the position within this grid.
                !
                ! Additionally, when the PAR is on top of a CON, we transport
                ! data to the neighbour grid as well.

                ! First - always send data to parent grid!
                inbr = iparent(igrid)
                iprocnbr = idprocofgrd(inbr)

                ! Position of fine grid within coarse grid
                SELECT CASE (iface)
                CASE (1, 2)
                    pos = iposition(igrid)
                CASE (3, 4)
                    pos = jposition(igrid)
                CASE (5, 6)
                    pos = kposition(igrid)
                END SELECT

                ! We allow a fine grid to completely refine a coarse grid
                ! in this case both PAR's of the fine grid can sit on top of
                ! a CON.
                !
                ! There are just two valid configurations:
                !   - iif == iic: This means the fine grid refine exactly half
                !     of the coarse grid in that direction. One PAR sits on
                !     top of a CON and another PAR is completely immersed in
                !     the coarse grid.
                !   - iif > iic: Actually, iif = (iic-4)*2 + 4, precisely.
                !     In this case the fine grid cover the entire extent of
                !     the coarse grid in the i-direction. Both PAR BC's sit on
                !     top of a CON in the coarse grid.
                one_connect = .TRUE.
                BLOCK
                    INTEGER :: iif, jjf, kkf, iic, jjc, kkc
                    CALL get_mgdims(kkf, jjf, iif, igrid)
                    CALL get_mgdims(kkc, jjc, iic, inbr)
                    SELECT CASE (iface)
                    CASE (1, 2)
                        IF (iif > iic) one_connect = .FALSE.
                    CASE (3, 4)
                        IF (jjf > jjc) one_connect = .FALSE.
                    CASE (5, 6)
                        IF (kkf > kkc) one_connect = .FALSE.
                    END SELECT
                END BLOCK

                ! In this "first" message, front send to front, back to back
                ! or it goes somewhere internal to the grid.
                ifacerecv = iface
                IF (pos == 3 .AND. one_connect) THEN
                    IF (iface == 2 .OR. iface == 4 .OR. iface == 6) &
                        ifacerecv = -1
                ELSE IF (pos > 3) THEN
                    IF (iface == 1 .OR. iface == 3 .OR. iface == 5) &
                        ifacerecv = -1
                END IF

                nsend = nsend + 1
                IF (nsend > maxconns) THEN
                    write(*, *) "Number of PAR's exceeded on process ", myid
                    write(*, *) "maxconns =", maxconns, &
                        "nmygrids =", nmygrids, "nsend = ", nsend
                    CALL errr(__FILE__, __LINE__)
                END IF

                sendconns(1, nsend) = myid
                sendconns(2, nsend) = iprocnbr
                sendconns(3, nsend) = igrid
                sendconns(4, nsend) = inbr
                sendconns(5, nsend) = iface
                sendconns(6, nsend) = ifacerecv

                sendcounts(iprocnbr) = sendcounts(iprocnbr) + ncols

                ! The same PAR face data can be sent a second time to a
                ! different grid and process.
                IF (pos == 3 .AND. one_connect) THEN
                    IF (iface == 2 .OR. iface == 4 .OR. iface == 6) CYCLE
                ELSE IF (pos > 3) THEN
                    IF (iface == 1 .OR. iface == 3 .OR. iface == 5) CYCLE
                END IF

                ! Second message/transfer follows.
                CALL get_neighbours(neighbours, inbr)
                inbr = neighbours(iface)
                iprocnbr = idprocofgrd(inbr)

                ! Get receiving face: front receive from back etc.
                SELECT CASE (iface)
                CASE (1)
                    ifacerecv = 2
                CASE (2)
                    ifacerecv = 1
                CASE (3)
                    ifacerecv = 4
                CASE (4)
                    ifacerecv = 3
                CASE (5)
                    ifacerecv = 6
                CASE (6)
                    ifacerecv = 5
                END SELECT

                nsend = nsend + 1
                IF (nsend > maxconns) THEN
                    write(*, *) "Number of PAR's exceeded on process ", myid
                    write(*, *) "maxconns =", maxconns, &
                        "nmygrids =", nmygrids, "nsend = ", nsend
                    CALL errr(__FILE__, __LINE__)
                END IF

                sendconns(1, nsend) = myid
                sendconns(2, nsend) = iprocnbr
                sendconns(3, nsend) = igrid
                sendconns(4, nsend) = inbr
                sendconns(5, nsend) = iface
                sendconns(6, nsend) = ifacerecv

                sendcounts(iprocnbr) = sendcounts(iprocnbr) + ncols
            END DO
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

        ! Allocate recvconns array
        irecv = (rdispls(numprocs-1) + recvcounts(numprocs-1))/ncols
        ALLOCATE(recvconns(ncols, irecv))
        recvconns = 0

        ! Exchange connection information
        CALL MPI_Alltoallv(sendconns, sendcounts, sdispls, MPI_INTEGER, &
            recvconns, recvcounts, rdispls, MPI_INTEGER, &
            MPI_COMM_WORLD)

        is_par_ftoc_core_init = .TRUE.

        DEALLOCATE(rdispls)
        DEALLOCATE(recvcounts)
        DEALLOCATE(sdispls)
        DEALLOCATE(sendcounts)
    END SUBROUTINE init_par_ftoc_core


    ! Calculate the number of values in the packed fine-grid face buffer.
    FUNCTION face_area(igrid, iface) RESULT(area)
        ! Function arguments
        INTEGER(intk) :: area
        INTEGER(intk), INTENT(in) :: igrid, iface

        ! Local variables
        INTEGER(intk) :: kstart, kstop, jstart, jstop, istart, istop

        CALL start_and_stop(kstart, kstop, jstart, jstop, istart, istop, &
            igrid, iface)

        area = ((istop-istart)/2+1)*((jstop-jstart)/2+1) &
            *((kstop-kstart)/2+1)
    END FUNCTION face_area


    ! Indices from which to send data from the fine grid
    SUBROUTINE start_and_stop(kstart, kstop, jstart, jstop, istart, istop, &
            igrid, iface)
        ! Subroutine arguments
        INTEGER(intk), INTENT(out) :: kstart, kstop, jstart, jstop, &
            istart, istop
        INTEGER(intk), INTENT(in) :: igrid, iface

        ! Local variables
        INTEGER(intk) :: kk, jj, ii

        CALL get_mgdims(kk, jj, ii, igrid)

        ! Select entire grid - later reduce to only contain face
        kstart = 3
        jstart = 3
        istart = 3

        kstop = kk-2
        jstop = jj-2
        istop = ii-2

        ! Reduce selection to only contain face
        ! The velocities here are always solved at the fine grid - so at the
        ! front the fine velocities at the PAR is stored in pos. 2, at the
        ! back in pos. ii-3 and so on for the other faces.
        SELECT CASE (iface)
        CASE (1)
            istart = 2
            istop = 2
        CASE (2)
            istart = ii-2
            istop = ii-2
        CASE (3)
            jstart = 2
            jstop = 2
        CASE (4)
            jstart = jj-2
            jstop = jj-2
        CASE (5)
            kstart = 2
            kstop = 2
        CASE (6)
            kstart = kk-2
            kstop = kk-2
        END SELECT
    END SUBROUTINE start_and_stop


    SUBROUTINE finish_par_ftoc_core()
        ! Subroutine arguments
        ! none...

        ! Local variables
        ! none...

        IF (.NOT. is_par_ftoc_core_init) CALL errr(__FILE__, __LINE__)

        isend = 0
        irecv = 0
        DEALLOCATE(sendconns)
        DEALLOCATE(recvconns)

        is_par_ftoc_core_init = .FALSE.
    END SUBROUTINE finish_par_ftoc_core
END MODULE par_ftoc_core_mod
