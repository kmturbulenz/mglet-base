
MODULE commbuf_mod

    USE, INTRINSIC :: ISO_FORTRAN_ENV, ONLY: int8
    USE, INTRINSIC :: ISO_C_BINDING, ONLY: C_PTR, C_F_POINTER, C_LOC, &
        C_ASSOCIATED, C_NULL_PTR, c_size_t, c_int

    USE MPI_f08
    USE omp_lib, ONLY: omp_get_default_device, omp_target_alloc, &
        omp_target_associate_ptr, omp_target_disassociate_ptr, omp_target_free

    USE precision_mod, ONLY: int64, intk, realk, int_bytes, real_bytes, &
        ifk, ifk_bytes
    USE pointers_mod, ONLY: idim3d
    USE err_mod, ONLY: errr

    IMPLICIT NONE (type, external)
    PRIVATE

    INTEGER(int64), PROTECTED :: idim_mg_bufs = 0
    INTEGER(int64), PROTECTED :: idim_mg_big = 0
    INTEGER(int64), PROTECTED :: idim_mg_intbuf = 0

    ! A 1-byte integer data buffer as a core for simplicity
    INTEGER(int8), ALLOCATABLE, TARGET :: buffer(:)
    TYPE(C_PTR) :: device_buffer = C_NULL_PTR

    ! Various buffers that all point to the same core buffer
    REAL(realk), POINTER, CONTIGUOUS :: sendbuf(:) => NULL()
    REAL(realk), POINTER, CONTIGUOUS :: recvbuf(:) => NULL()
    REAL(realk), POINTER, CONTIGUOUS :: bigbuf(:) => NULL()
    INTEGER(intk), POINTER, CONTIGUOUS :: intbuf(:) => NULL()
    !$omp declare target(sendbuf, recvbuf, bigbuf, intbuf)

    ! Device-native aliases used only as MPI buffer arguments
    REAL(realk), POINTER, CONTIGUOUS :: device_bigbuf(:) => NULL()
    REAL(realk), POINTER, CONTIGUOUS :: device_sendbuf(:) => NULL()
    REAL(realk), POINTER, CONTIGUOUS :: device_recvbuf(:) => NULL()

    INTEGER(ifk), POINTER, CONTIGUOUS :: ifkbuf(:) => NULL()
    INTEGER(ifk), POINTER, CONTIGUOUS :: isendbuf(:) => NULL()
    INTEGER(ifk), POINTER, CONTIGUOUS :: irecvbuf(:) => NULL()
    !$omp declare target(ifkbuf, isendbuf, irecvbuf)

    PUBLIC :: sendbuf, recvbuf, bigbuf, intbuf, device_sendbuf, device_recvbuf, &
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
        INTEGER(c_int) :: errorcode
        INTEGER :: device

        idim_mg_bufs = 0
        idim_mg_big = 0
        idim_mg_intbuf = 0
        NULLIFY(ifkbuf)
        NULLIFY(isendbuf)
        NULLIFY(irecvbuf)
        NULLIFY(sendbuf)
        NULLIFY(recvbuf)
        NULLIFY(bigbuf)
        NULLIFY(intbuf)
        NULLIFY(device_sendbuf)
        NULLIFY(device_recvbuf)
        NULLIFY(device_bigbuf)
        !$omp target update to(ifkbuf, isendbuf, irecvbuf, sendbuf, recvbuf, &
        !$omp& bigbuf, intbuf)

        device = omp_get_default_device()
        errorcode = omp_target_disassociate_ptr(C_LOC(buffer(1)), device)
        IF (errorcode /= 0) CALL errr(__FILE__, __LINE__)
        CALL omp_target_free(device_buffer, device)

        DEALLOCATE(buffer)
        device_buffer = C_NULL_PTR

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
        INTEGER(c_int) :: errorcode
        INTEGER :: device

        ! We correct the length to be a multiple of 32, to have a size that is
        ! dividable by two of quad prec reals
        corrlength = length + (32_int64 - MOD(length, 32_int64))

        ! Nullify all associated pointers, deallocate storage buffer and
        ! re-allocate the new length
        IF (ASSOCIATED(sendbuf)) NULLIFY(sendbuf)
        IF (ASSOCIATED(recvbuf)) NULLIFY(recvbuf)
        IF (ASSOCIATED(bigbuf)) NULLIFY(bigbuf)
        IF (ASSOCIATED(intbuf)) NULLIFY(intbuf)
        IF (ASSOCIATED(device_bigbuf)) NULLIFY(device_bigbuf)
        IF (ASSOCIATED(device_sendbuf)) NULLIFY(device_sendbuf)
        IF (ASSOCIATED(device_recvbuf)) NULLIFY(device_recvbuf)
        IF (ASSOCIATED(ifkbuf)) NULLIFY(ifkbuf)
        IF (ASSOCIATED(isendbuf)) NULLIFY(isendbuf)
        IF (ASSOCIATED(irecvbuf)) NULLIFY(irecvbuf)

        IF (ALLOCATED(buffer)) THEN
            device = omp_get_default_device()
            errorcode = omp_target_disassociate_ptr(C_LOC(buffer(1)), device)
            IF (errorcode /= 0) CALL errr(__FILE__, __LINE__)
            CALL omp_target_free(device_buffer, device)
            DEALLOCATE(buffer)
            device_buffer = C_NULL_PTR
        END IF

        ALLOCATE(buffer(corrlength))

        ! Allocating on device and checking device pointer association
        device = omp_get_default_device()
        device_buffer = omp_target_alloc(INT(corrlength, c_size_t), device)
        IF (.NOT. C_ASSOCIATED(device_buffer)) CALL errr(__FILE__, __LINE__)

        ! Telling OpenMP about the connections (H: buffer = D: device_buffer)
        errorcode = omp_target_associate_ptr(C_LOC(buffer(1)), device_buffer, &
            INT(corrlength, c_size_t), 0_c_size_t, device)
        IF (errorcode /= 0) CALL errr(__FILE__, __LINE__)

        idim_mg_big = corrlength/real_bytes
        idim_mg_bufs = idim_mg_big/2
        idim_mg_intbuf = corrlength/int_bytes
        cptr = C_LOC(buffer)

        CALL C_F_POINTER(cptr, bigbuf, [idim_mg_big])
        !$omp target update to(bigbuf)
        sendbuf => bigbuf(1:idim_mg_bufs)
        recvbuf => bigbuf(idim_mg_bufs+1:2*idim_mg_bufs)
        !$omp target update to(sendbuf, recvbuf)

        ! Conversion of device pointer to Fortran pointer
        CALL C_F_POINTER(device_buffer, device_bigbuf, [idim_mg_big])
        ! Setting the Fortran pointers to relevant parts of the device buffer
        device_sendbuf => device_bigbuf(1:idim_mg_bufs)
        device_recvbuf => device_bigbuf(idim_mg_bufs+1:2*idim_mg_bufs)

        CALL C_F_POINTER(cptr, intbuf, [idim_mg_intbuf])
        !$omp target update to(intbuf)

        ifklength = corrlength/ifk_bytes
        CALL C_F_POINTER(cptr, ifkbuf, [ifklength])
        !$omp target update to(ifkbuf)
        isendbuf => ifkbuf(1:ifklength/2)
        irecvbuf => ifkbuf(ifklength/2+1:2*(ifklength/2))
        !$omp target update to(isendbuf, irecvbuf)
    END SUBROUTINE allocate_buffer

END MODULE commbuf_mod
