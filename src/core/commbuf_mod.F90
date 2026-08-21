
MODULE commbuf_mod

    USE, INTRINSIC :: ISO_FORTRAN_ENV, ONLY: int8
    USE, INTRINSIC :: ISO_C_BINDING, ONLY: C_PTR, C_F_POINTER, C_LOC, &
        C_ASSOCIATED, C_NULL_PTR, c_size_t, c_int

    USE MPI_f08

#ifdef _MGLET_OFFLOAD_
    USE omp_lib, ONLY: omp_get_default_device, omp_get_initial_device, &
        omp_target_alloc, &
        omp_target_associate_ptr, omp_target_disassociate_ptr, omp_target_free
#endif

    USE precision_mod, ONLY: int64, intk, realk, int_bytes, real_bytes, &
        ifk, ifk_bytes, mglet_mpi_real
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
    !$omp declare target(intbuf)

    INTEGER(ifk), POINTER, CONTIGUOUS :: ifkbuf(:) => NULL()
    INTEGER(ifk), POINTER, CONTIGUOUS :: isendbuf(:) => NULL()
    INTEGER(ifk), POINTER, CONTIGUOUS :: irecvbuf(:) => NULL()
    !$omp declare target(ifkbuf, isendbuf, irecvbuf)

    ! Device-native aliases used only as MPI buffer arguments
    REAL(realk), POINTER, CONTIGUOUS :: device_bigbuf(:) => NULL()
    REAL(realk), POINTER, CONTIGUOUS :: device_sendbuf(:) => NULL()
    REAL(realk), POINTER, CONTIGUOUS :: device_recvbuf(:) => NULL()

#ifdef _MGLET_OFFLOAD_
    LOGICAL :: commbuf_uses_device = .FALSE.
#endif

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

#ifdef _MGLET_OFFLOAD_
    IF (commbuf_uses_device .AND. C_ASSOCIATED(device_buffer)) THEN
        device = omp_get_default_device()
        errorcode = omp_target_disassociate_ptr(C_LOC(buffer(1)), device)
        IF (errorcode /= 0) CALL errr(__FILE__, __LINE__)
        CALL omp_target_free(device_buffer, device)
    END IF
    commbuf_uses_device = .FALSE.
#endif

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
        INTEGER :: device, initial_device
        INTEGER :: nprobe
        INTEGER :: myrank, nranks

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
#ifdef _MGLET_OFFLOAD_
            IF (commbuf_uses_device .AND. C_ASSOCIATED(device_buffer)) THEN
                device = omp_get_default_device()
                errorcode = omp_target_disassociate_ptr(C_LOC(buffer(1)), device)
                IF (errorcode /= 0) CALL errr(__FILE__, __LINE__)
                CALL omp_target_free(device_buffer, device)
            END IF
            commbuf_uses_device = .FALSE.
#endif
            DEALLOCATE(buffer)
            device_buffer = C_NULL_PTR
        END IF

        ! Allocating the complete buffer and getting a C pointer to it
        ALLOCATE(buffer(corrlength))
        cptr = C_LOC(buffer)

        ! Computign dimensions for the various buffers
        idim_mg_big = corrlength/real_bytes
        idim_mg_bufs = idim_mg_big/2
        idim_mg_intbuf = corrlength/int_bytes

        CALL C_F_POINTER(cptr, intbuf, [idim_mg_intbuf])
        !$omp target update to(intbuf)

        ifklength = corrlength/ifk_bytes
        CALL C_F_POINTER(cptr, ifkbuf, [ifklength])
        !$omp target update to(ifkbuf)
        isendbuf => ifkbuf(1:ifklength/2)
        irecvbuf => ifkbuf(ifklength/2+1:2*(ifklength/2))
        !$omp target update to(isendbuf, irecvbuf)


#ifdef _MGLET_OFFLOAD_
        ! Runtime fallback for host-only mode: skip device association.
        initial_device = omp_get_initial_device()
        device = omp_get_default_device()
        commbuf_uses_device = .FALSE.
        device_buffer = C_NULL_PTR

        ! The complete buffer is allocated on the device, and the complete
        ! host buffer is associated with it.

        IF (device /= initial_device) THEN
            device_buffer = omp_target_alloc(INT(corrlength, c_size_t), device)
            IF (C_ASSOCIATED(device_buffer)) THEN
                errorcode = omp_target_associate_ptr(C_LOC(buffer(1)), &
                    device_buffer, INT(corrlength, c_size_t), 0_c_size_t, &
                    device)
                IF (errorcode == 0) THEN
                    commbuf_uses_device = .TRUE.
                ELSE
                    CALL omp_target_free(device_buffer, device)
                    device_buffer = C_NULL_PTR
                END IF
            END IF

            !$omp target update to(buffer)

        END IF


        CALL MPI_Barrier(MPI_COMM_WORLD)
        WRITE(*,*) "Test passed! (1)"
        CALL MPI_Barrier(MPI_COMM_WORLD)

        IF (commbuf_uses_device) THEN

            CALL C_F_POINTER(cptr, bigbuf, [idim_mg_big])
            sendbuf => bigbuf(1:idim_mg_bufs)
            recvbuf => bigbuf(idim_mg_bufs+1:2*idim_mg_bufs)

            CALL C_F_POINTER(device_buffer, device_bigbuf, [idim_mg_big])
            device_sendbuf => device_bigbuf(1:idim_mg_bufs)
            device_recvbuf => device_bigbuf(idim_mg_bufs+1:2*idim_mg_bufs)
        ELSE

            CALL C_F_POINTER(cptr, bigbuf, [idim_mg_big])
            sendbuf => bigbuf(1:idim_mg_bufs)
            recvbuf => bigbuf(idim_mg_bufs+1:2*idim_mg_bufs)

            device_sendbuf => sendbuf
            device_recvbuf => recvbuf
        END IF

        sendbuf(:) = -1.0
        recvbuf(:) = +1.0

        !$omp target update to(sendbuf, recvbuf)

        CALL MPI_Comm_rank(MPI_COMM_WORLD, myrank)
        CALL MPI_Comm_size(MPI_COMM_WORLD, nranks)

        CALL MPI_Barrier(MPI_COMM_WORLD)
        WRITE(*,*) "Test passed! (2)"
        CALL MPI_Barrier(MPI_COMM_WORLD)


        sendbuf(:) = 0.0_realk
        recvbuf(:) = 0.0_realk
        !$omp target update to(sendbuf, recvbuf)

        CALL MPI_Barrier(MPI_COMM_WORLD)
        WRITE(*,*) "Test passed! (3)"
        CALL MPI_Barrier(MPI_COMM_WORLD)

#else
        ! Setting these points to the host buffer
        device_sendbuf => sendbuf
        device_recvbuf => recvbuf
#endif

    END SUBROUTINE allocate_buffer

END MODULE commbuf_mod
