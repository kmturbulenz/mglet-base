MODULE blockbp_mod
    USE MPI_f08, ONLY: MPI_Wtime
    USE core_mod, ONLY: realk, intk, real64, config_t, mglet_filename_max, &
        errr, fort7, myid
    USE topol_mod, ONLY: topol_t
    IMPLICIT NONE(type, external)
    PRIVATE

    ! Common routines for blocking that IB models can inherit from
    TYPE, ABSTRACT :: blockbp_t
        LOGICAL :: do_blocking
        TYPE(config_t) :: blockconf

        INTEGER(intk) :: nfluidpoints
        REAL(realk), ALLOCATABLE :: fluidpoints(:, :)

        INTEGER(intk) :: nnofluidpoints
        REAL(realk), ALLOCATABLE :: nofluidpoints(:, :)

        CHARACTER(len=mglet_filename_max) :: blocking_type

        TYPE(topol_t) :: topol
        REAL(real64) :: time0
    CONTAINS
        PROCEDURE :: init
        PROCEDURE :: finish
        PROCEDURE :: read_velocity
        PROCEDURE, NOPASS :: set_icellspointer
    END TYPE blockbp_t

    PUBLIC :: blockbp_t
CONTAINS
    SUBROUTINE init(this)
        ! Subroutine arguments
        CLASS(blockbp_t), INTENT(inout) :: this

        ! Local variables
        INTEGER(intk) :: i
        LOGICAL :: has_nofluidpoints
        CHARACTER(len=64) :: jsonptr

        this%do_blocking = .FALSE.
        CALL fort7%get(this%blockconf, "/ib")

        this%blocking_type = REPEAT(" ", LEN(this%blocking_type))
        CALL this%blockconf%get_value("/blocking", this%blocking_type)
        IF (TRIM(this%blocking_type) == "newblock") THEN
            this%do_blocking = .TRUE.
        ELSE IF (TRIM(this%blocking_type) == "useblock") THEN
            RETURN
        ELSE
            IF (myid == 0) THEN
                WRITE(*, '("Invalid blocking type: ", A)') &
                    TRIM(this%blocking_type)
            END IF
            CALL errr(__FILE__, __LINE__)
        END IF

        this%time0 = MPI_Wtime()

        ! Get fluidpoints (required to be present if blocking)
        CALL this%blockconf%get_size("/fluidpoints", this%nfluidpoints)
        ALLOCATE(this%fluidpoints(3, this%nfluidpoints))
        DO i = 1, this%nfluidpoints
            WRITE(jsonptr, '("/fluidpoints/", I0)') i-1
            CALL this%blockconf%get_array(jsonptr, this%fluidpoints(:, i))
        END DO

        ! Get nofluidpoints (optional)
        has_nofluidpoints = this%blockconf%exists("/nofluidpoints")
        this%nnofluidpoints = 0
        IF (has_nofluidpoints) THEN
            CALL this%blockconf%get_size("/nofluidpoints", this%nnofluidpoints)
            ALLOCATE(this%nofluidpoints(3, this%nnofluidpoints))
            DO i = 1, this%nnofluidpoints
                WRITE(jsonptr, '("/nofluidpoints/", I0)') i-1
                CALL this%blockconf%get_array(jsonptr, &
                    this%nofluidpoints(:, i))
            END DO
        END IF

        IF (myid == 0) THEN
           WRITE (*,'("BLOCKING: ", A40, ", WTIME:", F16.3)') &
               'READING GEOMETRIES', MPI_Wtime() - this%time0
        END IF

        ! Read STL's
        CALL this%topol%init(this%blockconf)
    END SUBROUTINE init


    SUBROUTINE finish(this)
        ! Subroutine arguments
        CLASS(blockbp_t), INTENT(inout) :: this

        this%do_blocking = .FALSE.
        this%blocking_type = REPEAT(" ", LEN(this%blocking_type))
        this%nfluidpoints = 0
        this%nnofluidpoints = 0

        IF (ALLOCATED(this%fluidpoints)) DEALLOCATE(this%fluidpoints)
        IF (ALLOCATED(this%nofluidpoints)) DEALLOCATE(this%nofluidpoints)

        CALL this%topol%finish()
        CALL this%blockconf%finish()
    END SUBROUTINE finish


    SUBROUTINE read_velocity(this, velocity)
        ! Subroutine arguments
        CLASS(blockbp_t), INTENT(inout) :: this
        REAL(realk), ALLOCATABLE, INTENT(out) :: velocity(:, :)

        ! Local variables
        TYPE(config_t) :: geometry
        INTEGER(intk) :: i, ngeom, id, id_max
        LOGICAL :: has_velocity
        CHARACTER(len=64) :: jsonptr

        id_max = MAXVAL(this%topol%bodyid)
        ALLOCATE(velocity(3, id_max))
        velocity = 0.0

        CALL this%blockconf%get_size("/geometries", ngeom)
        DO i = 1, ngeom
            WRITE(jsonptr, '("/geometries/", I0)') i-1
            CALL this%blockconf%get(geometry, jsonptr)
            CALL geometry%get_value("/id", id)

            has_velocity = geometry%exists("/velocity")
            IF (has_velocity) THEN
                CALL geometry%get_array("/velocity", velocity(:, id))
            ELSE
                velocity(:, id) = 0.0
            END IF

            ! TODO: Set velocities that are read back into the JSON
            ! structure again

            CALL geometry%finish()
        END DO
    END SUBROUTINE read_velocity


    SUBROUTINE set_icellspointer(icells, icellspointer)
        ! Subroutine arguments
        INTEGER(intk), INTENT(in) :: icells(:)
        INTEGER(intk), INTENT(out) :: icellspointer(:)

        ! Local variables
        INTEGER(intk) :: i, ncells

        IF (SIZE(icells) /= SIZE(icellspointer)) THEN
            CALL errr(__FILE__, __LINE__)
        END IF

        icellspointer(1) = 1
        DO i = 2, SIZE(icells)
            icellspointer(i) = icellspointer(i-1) + icells(i-1)
        END DO

        ! All pointers that are out-of-bounds are set to the last element,
        ! and only for grids without intersected cells
        ncells = SUM(icells)
        DO i = 1, SIZE(icells)
            IF (icellspointer(i) > ncells) THEN
                ! Sanity check - grids with cells should not be out of
                ! bounds in the first place
                IF (icells(i) == 0) THEN
                    icellspointer(i) = ncells
                ELSE
                    CALL errr(__FILE__, __LINE__)
                END IF
            END IF
        END DO
    END SUBROUTINE set_icellspointer

END MODULE blockbp_mod
