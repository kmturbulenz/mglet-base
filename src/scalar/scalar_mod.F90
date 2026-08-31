MODULE scalar_mod
    ! Do not 'USE' core_mod or ib_mod here, because they will get exported
    ! (this module had no 'PRIVATE' default)
    USE timeintegrate_scalar_mod
    USE scacore_mod

    IMPLICIT NONE(type, external)

    ! PRIVATE ::

CONTAINS
    SUBROUTINE init_scalar()
        USE core_mod, ONLY: dcont, set_timer
        USE bound_scalar_mod, ONLY: init_bound_scalar
        USE itinfo_scalar_mod, ONLY: init_itinfo_scalar
        USE scastat_mod, ONLY: init_scastat
        USE gc_scastencils_mod
        USE ib_mod

        CALL init_scacore()
        IF (.NOT. has_scalar) RETURN

        CALL set_timer(400, "SCALAR")
        CALL set_timer(401, "SCA_PRE_TIME")
        CALL set_timer(402, "SCA_TIME")
        CALL set_timer(410, "SCA_TSTSCA4")
        CALL set_timer(411, "SCA_FLUXBALANCE")
        CALL set_timer(412, "SCA_STENCILS")
        CALL set_timer(413, "SCA_SOURCES")
        CALL set_timer(420, "SCA_ITINFO")

        CALL init_itinfo_scalar(dcont)
        CALL init_tfield()
        CALL init_bound_scalar()
        CALL init_scastat()

        ! Need to call this here - cannot be in scacore because that
        ! create a circular dependency
        SELECT TYPE(ib)
        TYPE IS (gc_t)
            IF (solve_scalar) THEN
                CALL init_scastencils(ib)
            END IF
        END SELECT
    END SUBROUTINE init_scalar


    SUBROUTINE finish_scalar
        USE bound_scalar_mod, ONLY: finish_bound_scalar
        USE itinfo_scalar_mod, ONLY: finish_itinfo_scalar
        USE scastat_mod, ONLY: finish_scastat
        USE ib_mod, ONLY: ib, gc_t
        USE gc_scastencils_mod, ONLY: finish_scastencils

        IF (.NOT. has_scalar) RETURN

        SELECT TYPE(ib)
        TYPE IS (gc_t)
            IF (solve_scalar) THEN
                CALL finish_scastencils()
            END IF
        END SELECT

        CALL finish_scastat()
        CALL finish_itinfo_scalar()
        CALL finish_bound_scalar()
        CALL finish_scacore()
    END SUBROUTINE finish_scalar


    SUBROUTINE init_tfield()
        USE core_mod, ONLY: conn, copy_arr, get_field, field_t, &
            map_arr_to_device, map_buffers_to_device, minlevel, maxlevel, &
            zero_ghostlayers
        USE ib_mod
        USE setboundarybuffers_scalar_mod, ONLY: setboundarybuffers_scalar

        TYPE(field_t), POINTER :: t, t_old
        INTEGER(intk) :: l, ilevel

        DO l = 1, nsca
            CALL get_field(t, scalar(l)%name)

            ! Write boundary conditions into buffers
            DO ilevel = minlevel, maxlevel
                CALL setboundarybuffers_scalar%bound(ilevel, t)
            END DO

            ! TODO: set initial condition when not dread

            CALL zero_ghostlayers(t)
            CALL map_arr_to_device(t, message="to:init_scalar")
            CALL map_buffers_to_device(t, message="to:init_scalar_buffers")

            DO ilevel = minlevel, maxlevel
                CALL conn(ilevel, 2, s1=t, corners=.TRUE.)
            END DO

            DO ilevel = maxlevel, minlevel+1, -1
                CALL ftoc(ilevel, t, t, 'T', device=.TRUE.)
            END DO

            CALL get_field(t_old, TRIM(scalar(l)%name)//"_OLD")
            CALL copy_arr(t_old%arr, t%arr)
        END DO
    END SUBROUTINE init_tfield
END MODULE scalar_mod
