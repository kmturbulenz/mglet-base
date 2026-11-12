MODULE conn_mod
    USE precision_mod
    USE timer_mod, ONLY: start_timer, stop_timer
    USE field_mod
    USE connect2_mod, ONLY: connect, init_connect2, finish_connect2
    USE conn_gpu_mod

    IMPLICIT NONE (type, external)
    PRIVATE

    PUBLIC :: conn, init_connect, finish_connect
CONTAINS
    ! Offload builds default to GPU execution; device=.FALSE. selects connect.
    ! Non-offload builds always use the CPU implementation.
    SUBROUTINE conn(ilevel, layers, v1, v2, v3, s1, s2, s3, corners, normal, &
            forward, ityp, device)
        ! Subroutine arguments
        INTEGER(intk), INTENT(in), OPTIONAL :: ilevel, layers
        TYPE(field_t), OPTIONAL, TARGET, INTENT(inout) :: v1, v2, v3, s1, s2, s3
        LOGICAL, OPTIONAL, INTENT(in) :: corners, normal
        INTEGER(intk), OPTIONAL, INTENT(in) :: forward
        CHARACTER(len=1), OPTIONAL, INTENT(in) :: ityp
        LOGICAL, OPTIONAL, INTENT(in) :: device

#ifdef _MGLET_OFFLOAD_
        ! Local variables
        LOGICAL :: device2

        device2 = .TRUE.
        IF (PRESENT(device)) device2 = device

        IF (device2) THEN
            CALL start_timer(150)
            CALL conn_gpu(ilevel, layers, v1, v2, v3, s1, s2, s3, corners, &
                normal, forward, ityp)
            CALL stop_timer(150)
            RETURN
        END IF
#endif
        CALL connect(ilevel, layers, v1, v2, v3, s1, s2, s3, &
            corners=corners, normal=normal, forward=forward, ityp=ityp)
    END SUBROUTINE conn


    SUBROUTINE init_connect()
        CALL init_connect2()
#ifdef _MGLET_OFFLOAD_
        CALL init_conn_gpu()
#endif
    END SUBROUTINE init_connect


    SUBROUTINE finish_connect()
#ifdef _MGLET_OFFLOAD_
        CALL finish_conn_gpu()
#endif
        CALL finish_connect2()
    END SUBROUTINE finish_connect
END MODULE conn_mod
