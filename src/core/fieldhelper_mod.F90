
MODULE fieldhelper_mod
    USE blasbind_mod, ONLY: memset
    USE err_mod, ONLY: errr
    USE field_mod, ONLY: field_t, intfield_t
    USE grids_mod, ONLY: get_mgdims, mygrids, nmygrids, level, get_imygrid
    USE pointers_mod, ONLY: get_ipbb, get_ip3
    USE precision_mod, ONLY: realk, intk, ifk
    USE profile_tools_mod, ONLY: profile_range_push, profile_range_pop

    IMPLICIT NONE(type, external)
    PRIVATE

    INTERFACE zero_field_arr
        PROCEDURE zero_field_arr_realk
        PROCEDURE zero_field_arr_ifk
    END INTERFACE zero_field_arr

    PUBLIC :: zero_field_arr, map_arr_to_device, map_arr_from_device, &
        map_buffers_to_device, map_buffers_from_device
CONTAINS
    SUBROUTINE zero_field_arr_realk(field, device)
        ! Subroutine arguments
        TYPE(field_t), INTENT(inout) :: field
        LOGICAL, OPTIONAL, INTENT(in) :: device

        ! Local variables
        LOGICAL :: device2

        IF (PRESENT(device)) THEN
            device2 = device
        ELSE
            device2 = .FALSE.
        END IF

#ifdef _MGLET_PROFILE_ANNOTATIONS_
            CALL profile_range_push("zero_field_arr_realk")
#endif

        IF (device2) THEN
#ifdef _MGLET_WORKAROUNDS_
            CALL memset(field%arr)
#else
            BLOCK
                INTEGER(intk) :: i
                ASSOCIATE(arr => field%arr)
                    !$omp target teams loop
                    DO i = 1, SIZE(arr)
                        arr(i) = 0.0_realk
                    END DO
                    !$omp end target teams loop
                END ASSOCIATE
            END BLOCK
#endif
        ELSE
            field%arr = 0.0_realk
        END IF

#ifdef _MGLET_PROFILE_ANNOTATIONS_
            CALL profile_range_pop()
#endif
    END SUBROUTINE zero_field_arr_realk


    SUBROUTINE zero_field_arr_ifk(field, device)
        ! Subroutine arguments
        TYPE(intfield_t), INTENT(inout) :: field
        LOGICAL, OPTIONAL, INTENT(in) :: device

        ! Local variables
        LOGICAL :: device2

        IF (PRESENT(device)) THEN
            device2 = device
        ELSE
            device2 = .FALSE.
        END IF

#ifdef _MGLET_PROFILE_ANNOTATIONS_
        CALL profile_range_push("zero_field_arr_ifk")
#endif

        IF (device2) THEN
#ifdef _MGLET_WORKAROUNDS_
            CALL memset(field%arr)
#else
            BLOCK
                INTEGER(intk) :: i
                ASSOCIATE(arr => field%arr)
                    !$omp target teams loop
                    DO i = 1, SIZE(arr)
                        arr(i) = 0_intk
                    END DO
                    !$omp end target teams loop
                END ASSOCIATE
            END BLOCK
#endif
        ELSE
            field%arr = 0_intk
        END IF

#ifdef _MGLET_PROFILE_ANNOTATIONS_
            CALL profile_range_pop()
#endif
    END SUBROUTINE zero_field_arr_ifk


    SUBROUTINE map_arr_to_device(f1, f2, f3, f4, f5, f6, f7, message)
        ! Subroutine arguments
        TYPE(field_t), INTENT(in) :: f1
        TYPE(field_t), INTENT(in), OPTIONAL :: f2, f3, f4, f5, f6, f7
        CHARACTER(*), INTENT(in), OPTIONAL :: message

        ! Local variables
        LOGICAL :: has_message

        has_message = PRESENT(message)

        IF (has_message) THEN
        ! Can not wrap the IF statement or the compiler will complain about
        ! an unused variable...
#if defined(_MGLET_PROFILE_ANNOTATIONS_) && defined(_MGLET_OFFLOAD_)
            CALL profile_range_push(message)
#endif
        END IF

        !$omp target update to(f1%arr)

        IF (PRESENT(f2)) THEN
            !$omp target update to(f2%arr)
        END IF
        IF (PRESENT(f3)) THEN
            !$omp target update to(f3%arr)
        END IF
        IF (PRESENT(f4)) THEN
            !$omp target update to(f4%arr)
        END IF
        IF (PRESENT(f5)) THEN
            !$omp target update to(f5%arr)
        END IF
        IF (PRESENT(f6)) THEN
            !$omp target update to(f6%arr)
        END IF
        IF (PRESENT(f7)) THEN
            !$omp target update to(f7%arr)
        END IF

#if defined(_MGLET_PROFILE_ANNOTATIONS_) && defined(_MGLET_OFFLOAD_)
        IF (has_message) THEN
            CALL profile_range_pop()
        END IF
#endif
    END SUBROUTINE map_arr_to_device


    SUBROUTINE map_buffers_to_device(f1, f2, f3, message)
        ! Subroutine arguments
        TYPE(field_t), INTENT(in) :: f1
        TYPE(field_t), INTENT(in), OPTIONAL :: f2, f3
        CHARACTER(*), INTENT(in), OPTIONAL :: message

        ! Local variables
        LOGICAL :: has_message

        has_message = PRESENT(message)

        IF (has_message) THEN
        ! Can not wrap the IF statement or the compiler will complain about
        ! an unused variable...
#if defined(_MGLET_PROFILE_ANNOTATIONS_) && defined(_MGLET_OFFLOAD_)
            CALL profile_range_push(message)
#endif
        END IF

        !$omp target update to(f1%buffers)

        IF (PRESENT(f2)) THEN
            !$omp target update to(f2%buffers)
        END IF
        IF (PRESENT(f3)) THEN
            !$omp target update to(f3%buffers)
        END IF

#if defined(_MGLET_PROFILE_ANNOTATIONS_) && defined(_MGLET_OFFLOAD_)
        IF (has_message) THEN
            CALL profile_range_pop()
        END IF
#endif
    END SUBROUTINE map_buffers_to_device


    SUBROUTINE map_arr_from_device(f1, f2, f3, f4, f5, f6, f7, message)
        ! Subroutine arguments
        TYPE(field_t), INTENT(in) :: f1
        TYPE(field_t), INTENT(in), OPTIONAL :: f2, f3, f4, f5, f6, f7
        CHARACTER(*), INTENT(in), OPTIONAL :: message

        ! Local variables
        LOGICAL :: has_message

        has_message = PRESENT(message)

        IF (has_message) THEN
        ! Can not wrap the IF statement or the compiler will complain about
        ! an unused variable...
#if defined(_MGLET_PROFILE_ANNOTATIONS_) && defined(_MGLET_OFFLOAD_)
            CALL profile_range_push(message)
#endif
        END IF

        !$omp target update from(f1%arr)

        IF (PRESENT(f2)) THEN
            !$omp target update from(f2%arr)
        END IF
        IF (PRESENT(f3)) THEN
            !$omp target update from(f3%arr)
        END IF
        IF (PRESENT(f4)) THEN
            !$omp target update from(f4%arr)
        END IF
        IF (PRESENT(f5)) THEN
            !$omp target update from(f5%arr)
        END IF
        IF (PRESENT(f6)) THEN
            !$omp target update from(f6%arr)
        END IF
        IF (PRESENT(f7)) THEN
            !$omp target update from(f7%arr)
        END IF

#if defined(_MGLET_PROFILE_ANNOTATIONS_) && defined(_MGLET_OFFLOAD_)
        IF (has_message) THEN
            CALL profile_range_pop()
        END IF
#endif
    END SUBROUTINE map_arr_from_device


    SUBROUTINE map_buffers_from_device(f1, f2, f3, message)
        ! Subroutine arguments
        TYPE(field_t), INTENT(in) :: f1
        TYPE(field_t), INTENT(in), OPTIONAL :: f2, f3
        CHARACTER(*), INTENT(in), OPTIONAL :: message

        ! Local variables
        LOGICAL :: has_message

        has_message = PRESENT(message)

        IF (has_message) THEN
        ! Can not wrap the IF statement or the compiler will complain about
        ! an unused variable...
#if defined(_MGLET_PROFILE_ANNOTATIONS_) && defined(_MGLET_OFFLOAD_)
            CALL profile_range_push(message)
#endif
        END IF

        !$omp target update from(f1%buffers)

        IF (PRESENT(f2)) THEN
            !$omp target update from(f2%buffers)
        END IF
        IF (PRESENT(f3)) THEN
            !$omp target update from(f3%buffers)
        END IF

#if defined(_MGLET_PROFILE_ANNOTATIONS_) && defined(_MGLET_OFFLOAD_)
        IF (has_message) THEN
            CALL profile_range_pop()
        END IF
#endif
    END SUBROUTINE map_buffers_from_device
END MODULE fieldhelper_mod
