MODULE ctof_mod
    USE core_mod
    USE MPI_f08

    IMPLICIT NONE (type, external)
    PRIVATE

    ! Variable to indicate if the required data structures and MPI-types
    ! have been created
    LOGICAL :: isInit = .FALSE.

    ! If .TRUE., a prolongation is already in process and you cannot start
    ! another one
    LOGICAL :: in_progress = .FALSE.

    ! Maximum allowed number of childs per parent (i.e. maximum number of
    ! send-conenctions per grid)
    INTEGER(intk), PARAMETER :: maxChilds = 8

    ! Lists that hold the send and receive request arrays
    TYPE(MPI_Request), ALLOCATABLE :: sendReqs(:), recvReqs(:)

    ! Actual number of messages that are to be sendt and received in one
    ! "round" of operations
    INTEGER(intk) :: nSend, nRecv

    ! List of grids to receive data on
    INTEGER(intk), ALLOCATABLE :: recvGrids(:)

    ! Datatypes for send- and receive-operations
    TYPE(MPI_Datatype), ALLOCATABLE :: sendTypes(:)
    TYPE(MPI_Datatype), ALLOCATABLE :: recvTypes(:)

    ! Field to prolongate
    REAL(realk), POINTER :: ff(:), fc(:)

    PUBLIC :: ctof, init_ctof, finish_ctof

CONTAINS
    SUBROUTINE ctof(ilevel, ff_p, fc_p)
        INTEGER(intk), INTENT(in) :: ilevel

        ! Fields to be interpolated
        REAL(realk), TARGET, INTENT(inout) :: ff_p(:)

        ! field on the fine grid
        REAL(realk), TARGET, INTENT(inout) :: fc_p(:)

        CALL start_timer(230)
        CALL ctof_begin(ff_p, fc_p, noflevel(ilevel), &
            igrdoflevel(1, ilevel))
        CALL ctof_end()
        CALL stop_timer(230)
    END SUBROUTINE ctof


    ! Initiate prolongation
    !
    ! Interpolate the results from a coarse level to a fine level
    ! This function initiate the process, and the ctof_finish
    ! must be called afterwards to clean up.
    SUBROUTINE ctof_begin(ff_p, fc_p, ngrids, lofgrids)
        ! Fields to be interpolated
        REAL(realk), TARGET, INTENT(inout) :: ff_p(:)
        REAL(realk), TARGET, INTENT(inout) :: fc_p(:)

        ! Number of grids ol level
        INTEGER(intk), INTENT(in) :: ngrids

        ! Grid ID's on level
        INTEGER(intk), INTENT(in) :: lofgrids(ngrids)

        ! Local variables
        INTEGER(intk) :: i, igrid, iprocc, iprocf, ipar

        CALL start_timer(231)

        IF (.NOT. isInit) THEN
            CALL errr(__FILE__, __LINE__)
        END IF

        IF (in_progress) THEN
            CALL errr(__FILE__, __LINE__)
        END IF

        in_progress = .TRUE.
        ff => ff_p
        fc => fc_p

        nRecv = 0
        nSend = 0

        ! Make all Recv-calls
        DO i = 1, ngrids
            igrid = lofgrids(i)
            ipar = iparent(igrid)
            IF (ipar /= 0 ) THEN
                iprocf = idprocofgrd(igrid)
                iprocc = idprocofgrd(ipar)

                IF (myid == iprocf) THEN
                    nRecv = nRecv + 1
                    CALL prolong_recv(igrid, ipar)
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

                IF (myid == iprocc) THEN
                    nSend = nSend + 1
                    CALL prolong_send(igrid, ipar)
                END IF
            END IF
        END DO

        CALL stop_timer(231)
    END SUBROUTINE ctof_begin


    ! Finish prolongation
    !
    ! Wait for communication to finish and clean up after prolongation
    !
    ! If the optional argument 'pbufs' is present the result from the
    ! prolongation will also be copied into the buffers PFR, PBA, PRI
    ! on the faces with 'parent' boundary conditions. This eliminate the
    ! need to update these ones manually afterward with a call to BPARMG,
    ! thus saving one communication stage.
    SUBROUTINE ctof_end(pbufs)
        ! USE copy_pbufs_mod, ONLY: copy_pbufs

        ! Subroutine arguments
        LOGICAL, OPTIONAL :: pbufs

        ! Local variables
        INTEGER(int32) :: idx

        CALL start_timer(232)

        IF (.NOT. in_progress) THEN
            CALL errr(__FILE__, __LINE__)
        END IF

        IF (nRecv > 0) THEN
            DO WHILE (.TRUE.)
                CALL MPI_Waitany(nRecv, recvReqs, idx, MPI_STATUS_IGNORE)

                IF (idx /= MPI_UNDEFINED) THEN
                    CALL prolong_finish(recvGrids(idx))

                    IF (PRESENT(pbufs)) THEN
                        ! TODO: Implement or remove
                        ! IF (pbufs) THEN
                        !     CALL copy_pbufs(recvGrids(idx), SIZE(ff), ff)
                        ! END IF
                        CALL errr(__FILE__, __LINE__)
                    END IF
                ELSE
                    EXIT
                END IF
            END DO
        END IF

        CALL MPI_Waitall(nSend, sendReqs, MPI_STATUSES_IGNORE)

        NULLIFY(ff)
        NULLIFY(fc)
        in_progress = .FALSE.

        CALL stop_timer(232)
    END SUBROUTINE ctof_end


    ! Initiate communication on receiver side
    SUBROUTINE prolong_recv(igridf, igridc)
        ! Suibroutine arguments
        INTEGER, INTENT(in) :: igridf
        INTEGER, INTENT(in) :: igridc

        ! Local variables
        INTEGER(intk) :: iprocc
        INTEGER(intk) :: ip3

        iprocc = idprocofgrd(igridc)

        CALL get_ip3(ip3, igridf)
        CALL MPI_Irecv(ff(ip3), 1, recvTypes(igridf), iprocc, igridf, &
            MPI_COMM_WORLD, recvReqs(nRecv))
    END SUBROUTINE prolong_recv


    ! Initiate communication on sender side
    SUBROUTINE prolong_send(igridf, igridc)
        ! Subroutine arguments
        INTEGER, INTENT(in) :: igridf
        INTEGER, INTENT(in) :: igridc

        ! Local variables
        INTEGER(intk) :: iprocf, iprocc
        INTEGER(intk) :: ip3

        iprocf = idprocofgrd(igridf)
        iprocc = idprocofgrd(igridc)

        CALL get_ip3(ip3, igridc)
        CALL MPI_Isend(fc(ip3), 1, sendTypes(igridf), iprocf, igridf, &
            MPI_COMM_WORLD, sendReqs(nSend))
    END SUBROUTINE prolong_send


    ! Finish prolongation, i.e. distribute the data on the grid
    SUBROUTINE prolong_finish(igridf)
        INTEGER, INTENT(in) :: igridf

        INTEGER(intk) :: ip3, pos, posc
        INTEGER(intk) :: k, j, i
        INTEGER(intk) :: kk, jj, ii
        INTEGER(intk) :: kc, jc, ic

        CALL get_ip3(ip3, igridf)
        CALL get_mgdims(kk, jj, ii, igridf)

        !$omp simd private(pos, posc, ic, jc, kc)
        DO i = 1, ii
            DO j = 1, jj
                DO k = 1, kk
                    pos = jj*kk*(i-1) + kk*(j-1) + (k-1)

                    ! Position in fine grid where the values from the
                    ! coarse grid is located
                    ic = i - (1 - MOD(i, 2))
                    jc = j - (1 - MOD(j, 2))
                    kc = k - (1 - MOD(k, 2))
                    posc = jj*kk*(ic-1) + kk*(jc-1) + (kc-1)
                    ff(ip3 + pos) = ff(ip3 + posc)
                END DO
            END DO
        END DO
    END SUBROUTINE prolong_finish


    ! Initialize arrays and data types
    SUBROUTINE init_ctof()
        ! Local variables
        INTEGER :: igrid, iprocc, iprocf, ipar

        CALL set_timer(230, "CTOF")
        CALL set_timer(231, "CTOF_BEGIN")
        CALL set_timer(232, "CTOF_END")

        IF (.NOT. isInit) THEN
            ALLOCATE(sendReqs(nMyGrids*maxChilds))
            ALLOCATE(recvReqs(nMyGrids))
            ALLOCATE(recvGrids(nMyGrids))
            ALLOCATE(sendTypes(ngrid))
            ALLOCATE(recvTypes(ngrid))

            sendTypes = MPI_DATATYPE_NULL
            recvTypes = MPI_DATATYPE_NULL
        END IF

        nRecv = 0
        nSend = 0

        ! Make all Send- and Recv-types
        DO igrid = 1, ngrid
            ipar = iparent(igrid)
            IF (ipar /= 0 ) THEN
                iprocf = idprocofgrd(igrid)
                iprocc = idprocofgrd(ipar)

                IF (myid == iprocf) THEN
                    nRecv = nRecv + 1
                    CALL create_recvtype(igrid)
                END IF
                IF (myid == iprocc) THEN
                    nSend = nSend + 1

                    IF (nSend > nMyGrids*maxChilds) THEN
                        CALL errr(__FILE__, __LINE__)
                    END IF

                    CALL create_sendtype(igrid, ipar)
                END IF
            END IF
        END DO

        isInit = .TRUE.
        in_progress = .FALSE.

        ! Nullify pointers
        nullify(ff)
        nullify(fc)
    END SUBROUTINE init_ctof


    ! Create MPI-datatypes for send-operation from a grid
    SUBROUTINE create_sendtype(igridf, igridc)
        ! Subroutine arguments
        INTEGER(intk), INTENT(in) :: igridf, igridc

        ! Local variables
        INTEGER(intk) :: kkc, jjc, iic
        INTEGER(intk) :: kkf, jjf, iif
        INTEGER(intk) :: ipos, jpos, kpos
        INTEGER(int32) :: array_of_sizes(3)
        INTEGER(int32) :: array_of_subsizes(3)
        INTEGER(int32) :: array_of_starts(3)

        ipos = iposition(igridf) - 1
        jpos = jposition(igridf) - 1
        kpos = kposition(igridf) - 1

        CALL get_mgdims(kkf, jjf, iif, igridf)
        CALL get_mgdims(kkc, jjc, iic, igridc)

        array_of_sizes(1) = kkc
        array_of_sizes(2) = jjc
        array_of_sizes(3) = iic

        array_of_subsizes(1) = kkf/2
        array_of_subsizes(2) = jjf/2
        array_of_subsizes(3) = iif/2

        ! NB: The starting value must be given in a C-index style, i.e.
        ! first element is 0, second element is 1 etc. Starting at Fortran
        ! index 2 means C index 1.
        array_of_starts(1) = kpos - 1
        array_of_starts(2) = jpos - 1
        array_of_starts(3) = ipos - 1

        CALL MPI_Type_create_subarray(3, array_of_sizes, array_of_subsizes, &
            array_of_starts, MPI_ORDER_FORTRAN,  mglet_mpi_real, &
            sendTypes(igridf))
        CALL MPI_Type_commit(sendTypes(igridf))
    END SUBROUTINE create_sendtype


    ! Create MPI-datatypes for receive operation onto a grid
    SUBROUTINE create_recvtype(igridf)
        ! Subroutine arguments
        INTEGER(intk), INTENT(in) :: igridf

        ! Local variables
        INTEGER(intk) :: kk, jj, ii
        TYPE(MPI_Datatype) :: ktype, jtype
        INTEGER(mpi_address_kind) :: stride

        CALL get_mgdims(kk, jj, ii, igridf)

        ! Temporary type for jumping in z-direction
        CALL MPI_Type_vector(kk/2, 1, 2, mglet_mpi_real, ktype)

        ! Temporary type for jumping in y-direction
        stride = 2*kk*real_bytes
        CALL MPI_Type_create_hvector(jj/2, 1, stride, ktype, jtype)

        ! Final type including jumps in x-direction
        stride = 2*kk*jj*real_bytes
        CALL MPI_Type_create_hvector(ii/2, 1, stride, jtype, recvTypes(igridf))
        CALL MPI_Type_commit(recvTypes(igridf))

        ! Free temporary types
        CALL MPI_Type_free(ktype)
        CALL MPI_Type_free(jtype)
    END SUBROUTINE create_recvtype


    SUBROUTINE finish_ctof()
        INTEGER(intk) :: i

        IF (isInit .NEQV. .TRUE.) THEN
            RETURN
        END IF

        isInit = .FALSE.

        ! Free datatypes
        DO i = 1, SIZE(sendTypes)
            IF (sendTypes(i) /= MPI_DATATYPE_NULL) THEN
                CALL MPI_Type_free(sendTypes(i))
            END IF
        END DO
        DO i = 1, SIZE(recvTypes)
            IF (recvTypes(i) /= MPI_DATATYPE_NULL) THEN
                CALL MPI_Type_free(recvTypes(i))
            END IF
        END DO

        DEALLOCATE(sendReqs)
        DEALLOCATE(recvReqs)
        DEALLOCATE(recvGrids)
        DEALLOCATE(sendTypes)
        DEALLOCATE(recvTypes)
    END SUBROUTINE finish_ctof
END MODULE ctof_mod
