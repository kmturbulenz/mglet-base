MODULE commbuf_mod
    USE, INTRINSIC :: ISO_FORTRAN_ENV, ONLY: int8
    USE, INTRINSIC :: ISO_C_BINDING, ONLY: C_PTR, C_F_POINTER, C_LOC
    USE MPI_f08
#ifdef _MGLET_DEVICE_BUFFER_
    USE omp_lib
#endif

    USE precision_mod, ONLY: int64, intk, realk, int_bytes, real_bytes, &
        ifk, ifk_bytes
    USE pointers_mod, ONLY: idim3d

    IMPLICIT NONE (type, external)
    PRIVATE

    INTEGER(int64), PROTECTED :: idim_mg_bufs = 0
    INTEGER(int64), PROTECTED :: idim_mg_big = 0
    INTEGER(int64), PROTECTED :: idim_mg_intbuf = 0

    ! A 1-byte integer data buffer as a core for simplicity
#ifdef _MGLET_DEVICE_BUFFER_
    TYPE(C_PTR) :: buffer
    INTEGER(intk) :: idev
#else
    INTEGER(int8), ALLOCATABLE, TARGET :: buffer(:)
#endif

    ! Various buffers that all point to the same core buffer
    REAL(realk), POINTER, CONTIGUOUS :: sendbuf(:) => NULL()
    REAL(realk), POINTER, CONTIGUOUS :: recvbuf(:) => NULL()
    REAL(realk), POINTER, CONTIGUOUS :: bigbuf(:) => NULL()
    INTEGER(intk), POINTER, CONTIGUOUS :: intbuf(:) => NULL()

    INTEGER(ifk), POINTER, CONTIGUOUS :: ifkbuf(:) => NULL()
    INTEGER(ifk), POINTER, CONTIGUOUS :: isendbuf(:) => NULL()
    INTEGER(ifk), POINTER, CONTIGUOUS :: irecvbuf(:) => NULL()

    PUBLIC :: sendbuf, recvbuf, bigbuf, intbuf, &
        idim_mg_bufs, idim_mg_big, idim_mg_intbuf, &
        increase_bigbuf, increase_intbuf, increase_ifkbuf, &
        init_commbuf, finish_commbuf, ifkbuf, isendbuf, irecvbuf

CONTAINS
    SUBROUTINE init_commbuf()
        ! Local variables
        INTEGER(int64) :: bigbuflen

        ! All processes allocate the same buffer
        bigbuflen = 6*idim3d
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
#ifdef _MGLET_DEVICE_BUFFER_
        CALL omp_target_free(buffer, idev)
        idev = -1
#else
        DEALLOCATE(buffer)
#endif
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

#ifdef _MGLET_DEVICE_BUFFER_
        IF (omp_target_is_present(buffer, idev) /= 0) THEN
            CALL omp_target_free(buffer, idev)
            idev = -1
        END IF
        idev = omp_get_default_device()
        buffer = omp_target_alloc(corrlength, idev)
        cptr = buffer
#else
        IF (ALLOCATED(buffer)) THEN
            DEALLOCATE(buffer)
        END IF
        ALLOCATE(buffer(corrlength))
        cptr = C_LOC(buffer)
#endif

        idim_mg_big = corrlength/real_bytes
        idim_mg_bufs = idim_mg_big/2
        idim_mg_intbuf = corrlength/int_bytes

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
