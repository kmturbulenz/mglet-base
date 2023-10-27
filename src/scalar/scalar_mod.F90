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
        CALL set_timer(420, "SCA_ITINFO")

        CALL init_itinfo_scalar(dcont)
        CALL init_tfield()
        CALL init_scastat()

        ! Need to call this here - cannot be in scacore because that
        ! create a circular dependency
        SELECT TYPE(ib)
        TYPE IS (gc_t)
            IF (solve_scalar) THEN
                CALL create_scastencils(ib)
            END IF
        END SELECT
    END SUBROUTINE init_scalar


    SUBROUTINE finish_scalar
        USE itinfo_scalar_mod, ONLY: finish_itinfo_scalar
        USE scastat_mod, ONLY: finish_scastat

        IF (.NOT. has_scalar) RETURN

        CALL finish_scastat()
        CALL finish_itinfo_scalar()
        CALL finish_scacore()
    END SUBROUTINE finish_scalar


    SUBROUTINE init_tfield()
        USE core_mod, ONLY: connect, get_field, field_t, &
            minlevel, maxlevel, zero_ghostlayers
        USE ib_mod
        USE setboundarybuffers_scalar_mod, ONLY: setboundarybuffers_scalar

        TYPE(field_t), POINTER :: t
        INTEGER(intk) :: l, ilevel

        DO l = 1, nsca
            CALL get_field(t, scalar(l)%name)

            ! Write boundary conditions into buffers
            DO ilevel = minlevel, maxlevel
                CALL setboundarybuffers_scalar%bound(ilevel, t)
            END DO

            ! TODO: set initial condition when not dread

            CALL zero_ghostlayers(t)

            DO ilevel = minlevel, maxlevel
                CALL connect(ilevel, 2, s1=t, corners=.TRUE.)
            END DO

            DO ilevel = maxlevel, minlevel+1, -1
                CALL ftoc(ilevel, t%arr, t%arr, 'T')
            END DO
        END DO
    END SUBROUTINE init_tfield
END MODULE scalar_mod
