MODULE scacore_mod
    USE core_mod
    USE ib_mod, ONLY: ftoc, parent, ib

    IMPLICIT NONE(type, external)
    PRIVATE

    ! Control parameters
    LOGICAL, PROTECTED :: has_sca = .FALSE.
    LOGICAL, PROTECTED :: solve_sca = .FALSE.
    LOGICAL, PROTECTED :: dread_sca = .FALSE.
    LOGICAL, PROTECTED :: dwrite_sca = .FALSE.
    LOGICAL, PROTECTED :: dcont_sca = .FALSE.

    ! Scalar/physical paramters
    INTEGER(intk), PARAMETER :: nsca_max = 8
    INTEGER(intk), PROTECTED :: nsca
    REAL(realk), PROTECTED :: D(nsca_max)  ! TODO: rename appropriately

    PUBLIC :: init_scacore, finish_scacore, has_sca, dread_sca, dcont_sca, &
        solve_sca

CONTAINS
    SUBROUTINE init_scacore()
        ! Subroutine arguments
        ! None...

        ! Local variables
        TYPE(config_t) :: scaconf
        INTEGER(intk) :: t
        CHARACTER(len=8) :: fieldname

        ! Read configuration values - if not exists no timeintegration is
        ! performed
        has_sca = .FALSE.
        IF (.NOT. fort7%exists("/scalar")) THEN
            IF (myid == 0) THEN
                WRITE(*, '("NO SCALAR")')
                WRITE(*, '()')
            END IF
            RETURN
        END IF
        has_sca = .TRUE.

        ! Required values
        scaconf = fort7%get("/scalar")
        CALL scaconf%get_value("/nsca", nsca)
        IF (nsca > nsca_max) THEN
            WRITE(*,*) "nsca: ", nsca
            WRITE(*,*) "nsca_max: ", nsca_max
            CALL errr(__FILE__, __LINE__)
        END IF
        CALL scaconf%get_array("/D", D(1:nsca))

        ! Optional values
        CALL scaconf%get_value("/solve", solve_sca, .TRUE.)

        ! If not reading anything, dcont is always false
        IF (.NOT. dread_sca) THEN
            dcont_sca = .FALSE.
        END IF

        DO t = 1, nsca
            WRITE(fieldname, '("T", I0)') t
            CALL set_field(fieldname, dread=dread_sca, required=dread_sca, &
                dwrite=dwrite_sca)

            ! For RK time integration
            WRITE(fieldname, '("DT", I0)') t
            CALL set_field(fieldname)
        END DO
    END SUBROUTINE init_scacore


    SUBROUTINE finish_scacore
        CONTINUE
    END SUBROUTINE finish_scacore

END MODULE scacore_mod
