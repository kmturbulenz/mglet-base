
MODULE fieldhelper_mod
    USE err_mod, ONLY: errr
    USE field_mod, ONLY: field_t, intfield_t
    USE grids_mod, ONLY: get_mgdims, mygrids, nmygrids, level, get_imygrid
    USE pointers_mod, ONLY: get_ipbb
    USE precision_mod, ONLY: realk, intk, ifk
    USE profile_tools_mod, ONLY: profile_range_push, profile_range_pop

    IMPLICIT NONE(type, external)
    PRIVATE

    INTERFACE
        SUBROUTINE set_field_arr_realk_backend(arr, val) &
                BIND(C, name="set_field_arr_realk_c")
            USE, INTRINSIC :: iso_c_binding, ONLY: c_int64_t
            USE precision_mod, ONLY: c_realk
            REAL(c_realk), INTENT(inout) :: arr(:)
            REAL(c_realk), VALUE, INTENT(in) :: val
        END SUBROUTINE set_field_arr_realk_backend

        SUBROUTINE set_field_arr_ifk_backend(arr, val) &
                BIND(C, name="set_field_arr_ifk_c")
            USE, INTRINSIC :: iso_c_binding, ONLY: c_int64_t
            USE precision_mod, ONLY: c_ifk
            INTEGER(c_ifk), INTENT(inout) :: arr(:)
            INTEGER(c_ifk), VALUE, INTENT(in) :: val
        END SUBROUTINE set_field_arr_ifk_backend
    END INTERFACE

    INTERFACE set_field_arr
        PROCEDURE set_field_arr_realk
        PROCEDURE set_field_arr_ifk
    END INTERFACE set_field_arr

    PUBLIC :: set_field_arr, map_arr_to_device, map_arr_from_device, &
        map_buffers_to_device, map_buffers_from_device
CONTAINS
    SUBROUTINE set_field_arr_realk(field, val, device)
        ! Subroutine arguments
        TYPE(field_t), INTENT(inout) :: field
        REAL(realk), INTENT(in) :: val
        LOGICAL, OPTIONAL, INTENT(in) :: device

        ! Local variables
        INTEGER(intk) :: i, n
        LOGICAL :: device2

        IF (PRESENT(device)) THEN
            device2 = device
        ELSE
            device2 = .FALSE.
        END IF

        IF (device2) THEN
#ifdef _MGLET_USE_BACKEND_
            CALL set_field_arr_realk_backend(field%arr, val)
#else
            n = SIZE(field%arr)
            ASSOCIATE(arr => field%arr)
#ifdef _MGLET_PROFILE_ANNOTATIONS_
            CALL profile_range_push("set_field_arr_realk")
#endif
            !$omp target teams loop
            DO i = 1, n
                arr(i) = val
            END DO
            !$omp end target teams loop
#ifdef _MGLET_PROFILE_ANNOTATIONS_
            CALL profile_range_pop()
#endif
            END ASSOCIATE
#endif
        ELSE
            field%arr = val
        END IF
    END SUBROUTINE set_field_arr_realk


    SUBROUTINE set_field_arr_ifk(field, val, device)
        ! Subroutine arguments
        TYPE(intfield_t), INTENT(inout) :: field
        INTEGER(ifk), INTENT(in) :: val
        LOGICAL, OPTIONAL, INTENT(in) :: device

        ! Local variables
        INTEGER(intk) :: i, n
        LOGICAL :: device2

        IF (PRESENT(device)) THEN
            device2 = device
        ELSE
            device2 = .FALSE.
        END IF

        IF (device2) THEN
#ifdef _MGLET_USE_BACKEND_
            CALL set_field_arr_ifk_backend(field%arr, val)
#else
            n = SIZE(field%arr)
            ASSOCIATE(arr => field%arr)
#ifdef _MGLET_PROFILE_ANNOTATIONS_
            CALL profile_range_push("set_field_arr_ifk")
#endif
            !$omp target teams loop
            DO i = 1, n
                arr(i) = val
            END DO
            !$omp end target teams loop
#ifdef _MGLET_PROFILE_ANNOTATIONS_
            CALL profile_range_pop()
#endif
            END ASSOCIATE
#endif
        ELSE
            field%arr = val
        END IF
    END SUBROUTINE set_field_arr_ifk


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
