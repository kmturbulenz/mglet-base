MODULE conn_mod

    USE precision_mod
    USE timer_mod, ONLY: start_timer, stop_timer
    USE field_mod
    USE conn1_mod
    USE conn2_mod
    USE conn3_mod

    PUBLIC :: conn, init_conn, finish_conn

CONTAINS

    SUBROUTINE conn(ilevel, layers, v1, v2, v3, s1, s2, s3, corners, normal, &
            forward, ityp, record)

        ! This subroutine acts as an interface for connect calls.
        ! Dependent on the configuration, it will call either conn1_mod
        ! or conn2_mod, which are relevant for CPU and GPU, respectively.

        ! Subroutine arguments
        INTEGER(intk), INTENT(in), OPTIONAL :: ilevel, layers
        TYPE(field_t), OPTIONAL, TARGET, INTENT(inout) :: v1, v2, v3, s1, s2, s3
        LOGICAL, OPTIONAL, INTENT(in) :: corners, normal
        INTEGER(intk), OPTIONAL, INTENT(in) :: forward
        CHARACTER(len=1), OPTIONAL, INTENT(in) :: ityp
        LOGICAL, OPTIONAL, INTENT(in) :: record

        CALL start_timer(150)
#ifdef _MGLET_OFFLOAD_
        CALL conn3(ilevel, layers, v1, v2, v3, s1, s2, s3, corners, &
            normal, forward, ityp, record)
#else
        ! CALL conn1(ilevel, layers, v1, v2, v3, s1, s2, s3, corners, &
        !     normal, forward, ityp)
        CALL conn3(ilevel, layers, v1, v2, v3, s1, s2, s3, corners, &
            normal, forward, ityp, record)
#endif
        CALL stop_timer(150)
    END SUBROUTINE conn


    SUBROUTINE init_conn()
#ifdef _MGLET_OFFLOAD_
        CALL init_conn3()
#else
        CALL init_conn3()
#endif
    END SUBROUTINE init_conn


    SUBROUTINE finish_conn()
#ifdef _MGLET_OFFLOAD_
        CALL finish_conn3()
#else
        CALL finish_conn3()
#endif
    END SUBROUTINE finish_conn

END MODULE conn_mod
