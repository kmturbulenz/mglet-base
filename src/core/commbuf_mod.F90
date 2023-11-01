MODULE commbuf_mod
    USE, INTRINSIC :: ISO_FORTRAN_ENV, ONLY: int8
    USE, INTRINSIC :: ISO_C_BINDING, ONLY: C_PTR, C_F_POINTER, C_LOC
    USE MPI_f08

    USE precision_mod, ONLY: int64, intk, realk, int_bytes, real_bytes, &
        ifk, ifk_bytes
    USE pointers_mod, ONLY: idim2d, idim3d

    IMPLICIT NONE (type, external)
    PRIVATE

    ! Maximum number of variables in any connect-call
    INTEGER(intk) :: maxnvars = 6

    INTEGER(int64), PROTECTED :: idim_mg_bufs = 0
    INTEGER(int64), PROTECTED :: idim_mg_big = 0
    INTEGER(int64), PROTECTED :: idim_mg_intbuf = 0

    ! A 1-byte integer data buffer as a core for simplicity
    INTEGER(int8), ALLOCATABLE, TARGET :: buffer(:)

    ! Various buffers that all point to the same core buffer
    REAL(realk), POINTER, CONTIGUOUS :: sendbuf(:)
    REAL(realk), POINTER, CONTIGUOUS :: recvbuf(:)
    REAL(realk), POINTER, CONTIGUOUS :: bigbuf(:)
    INTEGER(intk), POINTER, CONTIGUOUS :: intbuf(:)

    INTEGER(ifk), POINTER, CONTIGUOUS :: ifkbuf(:)
    INTEGER(ifk), POINTER, CONTIGUOUS :: isendbuf(:)
    INTEGER(ifk), POINTER, CONTIGUOUS :: irecvbuf(:)

    PUBLIC :: sendbuf, recvbuf, bigbuf, intbuf, &
        idim_mg_bufs, idim_mg_big, idim_mg_intbuf, &
        increase_bigbuf, increase_intbuf, increase_ifkbuf, &
        init_commbuf, finish_commbuf, ifkbuf, isendbuf, irecvbuf

CONTAINS
    SUBROUTINE init_commbuf()
        ! Local variables
        INTEGER(int64) :: commbuflen, bigbuflen

        ! Max memory requirement for sendbuf/recvbuf in connect:
        !   - One 2D slice of each grid (idim2d)
        !   - 2 planes per face
        !   - 6 faces
        !   - times the maximum number of variables
        ! This calculation overestimate the buffer if many levels
        ! are used, since the idim2d variable include grids on all
        ! levels, and connect only work on one level each time.
        !
        ! original formulation: 1.2*idim2d*2*6*maxnvars
        ! only INT's: 1.2*2*6 = 14.4 round up to 15
        commbuflen = 15*idim2d*maxnvars

        ! Max memory requirement for bigbuf: All 3D fields at once,
        ! rounded up to a number dividable by two
        bigbuflen = idim3d + MOD(idim3d,2) + 2

        ! All processes allocate the same buffer
        bigbuflen = MAX(bigbuflen, 2*commbuflen)
        CALL MPI_Allreduce(MPI_IN_PLACE, bigbuflen, 1, MPI_INTEGER8, &
            mpi_max, MPI_COMM_WORLD)

        CALL increase_bigbuf(bigbuflen)
    END SUBROUTINE init_commbuf


    SUBROUTINE finish_commbuf()
        idim_mg_bufs = 0
        idim_mg_big = 0
        idim_mg_intbuf = 0
        NULLIFY(sendbuf)
        NULLIFY(recvbuf)
        NULLIFY(intbuf)
        NULLIFY(recvbuf)
        DEALLOCATE(buffer)
    END SUBROUTINE finish_commbuf


    SUBROUTINE increase_bigbuf(length)
        ! Length of bigbuf in number of elements
        INTEGER(kind=int64), INTENT(in) :: length
        IF (length > idim_mg_big) THEN
            CALL allocate_buffer(length*real_bytes)
        END IF
    END SUBROUTINE increase_bigbuf


    SUBROUTINE increase_intbuf(length)
        ! Length of bigbuf in number of elements
        INTEGER(kind=int64), INTENT(in) :: length

        IF (length > idim_mg_big) THEN
            CALL allocate_buffer(length*int_bytes)
        END IF
    END SUBROUTINE increase_intbuf


    SUBROUTINE increase_ifkbuf(length)
        ! Length of bigbuf in number of elements
        INTEGER(kind=int64), INTENT(in) :: length

        IF (length > SIZE(ifkbuf)) THEN
            CALL allocate_buffer(length*ifk_bytes)
        END IF
    END SUBROUTINE increase_ifkbuf


    SUBROUTINE allocate_buffer(length)
        ! Length in *bytes*
        INTEGER(int64), INTENT(in) :: length

        ! Local variables
        INTEGER(int64) :: corrlength
        TYPE(C_PTR) :: cptr
        INTEGER(int64) :: ifklength

        ! We correct the length to be a multiple of 32, to have a size that is
        ! dividable by two of quad prec reals
        corrlength = length + (32_int64 - MOD(length, 32_int64))

        ! Nullify all associated pointers, deallocate storage buffer and
        ! re-allocate the new length
        IF (ASSOCIATED(sendbuf)) NULLIFY(sendbuf)
        IF (ASSOCIATED(recvbuf)) NULLIFY(recvbuf)
        IF (ASSOCIATED(bigbuf)) NULLIFY(bigbuf)
        IF (ASSOCIATED(intbuf)) NULLIFY(intbuf)
        IF (ASSOCIATED(ifkbuf)) NULLIFY(ifkbuf)
        IF (ASSOCIATED(isendbuf)) NULLIFY(isendbuf)
        IF (ASSOCIATED(irecvbuf)) NULLIFY(irecvbuf)

        IF (ALLOCATED(buffer)) THEN
            DEALLOCATE(buffer)
        END IF
        ALLOCATE(buffer(corrlength))

        idim_mg_big = corrlength/real_bytes
        idim_mg_bufs = idim_mg_big/2
        idim_mg_intbuf = corrlength/int_bytes
        cptr = C_LOC(buffer)

        CALL C_F_POINTER(cptr, bigbuf, [idim_mg_big])
        sendbuf => bigbuf(1:idim_mg_bufs)
        recvbuf => bigbuf(idim_mg_bufs+1:2*idim_mg_bufs)

        CALL C_F_POINTER(cptr, intbuf, [idim_mg_intbuf])

        ifklength = corrlength/ifk_bytes
        CALL C_F_POINTER(cptr, ifkbuf, [ifklength])
        isendbuf => ifkbuf(1:ifklength/2)
        irecvbuf => ifkbuf(ifklength/2+1:2*(ifklength/2))
    END SUBROUTINE allocate_buffer

END MODULE commbuf_mod
