MODULE scacore_mod
    USE core_mod
    USE ib_mod, ONLY: ib, gc_t
    USE flow_mod, ONLY: has_flow

    IMPLICIT NONE(type, external)
    PRIVATE

    TYPE :: scalar_bc_t
        ! Integer flag corresponding to type of BC
        !   flag = 0: Fixed value (Diriclet) boundary condition ("value")
        !   flag = 1: Fixed flux (Neumann) boundary condition ("flux")
        INTEGER(intk) :: flag

        ! Actual value or flux
        REAL(realk) :: value
    END TYPE scalar_bc_t

    TYPE :: scalar_source_t
        ! Actual value or sourceterm, either as absolute value or
        ! proportionality factor
        REAL(realk) :: value

        ! Fieldname for proportionality factor
        CHARACTER(len=nchar_name) :: field
    END TYPE scalar_source_t

    TYPE :: scalar_t
        CHARACTER(len=nchar_name) :: name
        REAL(realk) :: prmol
        INTEGER(intk) :: units(7)
        INTEGER(intk) :: kayscrawford
        TYPE(scalar_bc_t), ALLOCATABLE :: geometries(:)
        TYPE(scalar_source_t), ALLOCATABLE :: sources(:)
    CONTAINS
        PROCEDURE :: prt
    END TYPE scalar_t

    ! Control parameters
    LOGICAL, PROTECTED :: has_scalar = .FALSE.
    LOGICAL, PROTECTED :: solve_scalar = .FALSE.
    REAL(realk) :: prturb

    ! Scalar/physical paramters
    INTEGER(intk), PROTECTED :: nsca
    TYPE(scalar_t), ALLOCATABLE, PROTECTED :: scalar(:)

    PUBLIC :: init_scacore, finish_scacore, scalar_t, scalar, nsca, prturb, &
        has_scalar, solve_scalar, maskbt, scalar_source_t

CONTAINS
    SUBROUTINE init_scacore()
        USE blockbt_mod
        ! Subroutine arguments
        ! None...

        ! Local variables
        TYPE(config_t) :: scaconf, sc, source
        TYPE(field_t), POINTER :: t, t_old, bt
        INTEGER(intk) :: l, n, nstl, itype, nsource
        CHARACTER(len=mglet_filename_max + 64) :: jsonptr
        CHARACTER(len=16) :: type
        CHARACTER(len=10240) :: initial_expr
        LOGICAL :: kayscrawford
        REAL(realk) :: rvalue

        ! Read configuration values - if not exists no timeintegration is
        ! performed
        has_scalar = .FALSE.
        IF (.NOT. fort7%exists("/scalar")) THEN
            IF (myid == 0) THEN
                WRITE(*, '("NO SCALAR")')
                WRITE(*, '()')
            END IF
            RETURN
        END IF
        has_scalar = .TRUE.

        ! Scalar is dependent on flow to de defined (but not neccesarily solved)
        ! to have variables like gmol and rho defined
        IF (.NOT. has_flow) THEN
            WRITE(*, *) "Scalar needs flow!"
            CALL errr(__FILE__, __LINE__)
        END IF

        ! Required values
        CALL fort7%get(scaconf, "/scalar")
        CALL scaconf%get_size("/scalars", nsca)
        ALLOCATE(scalar(nsca))

        ! Read parameters per scalar
        DO l = 1, nsca
            WRITE(jsonptr, '("/scalars/", I0)') l-1
            CALL scaconf%get(sc, jsonptr)

            CALL sc%get_value("/name", scalar(l)%name)
            ! To allow for field "{name}_AVG" to fit in nchar_name characters
            IF (LEN_TRIM(scalar(l)%name) + 4 > nchar_name) THEN
                CALL errr(__FILE__, __LINE__)
            END IF

            CALL sc%get_value("/prmol", scalar(l)%prmol)

            ! Retrieving the scalar units (noisy output for test)
            IF (sc%exists("/units")) THEN
                CALL sc%get_array("/units", scalar(l)%units)
            ELSE
                scalar(l)%units = 0
            END IF

            ! The parameters.json should have a logical, but it is easier
            ! to integrate computations with an integer
            CALL sc%get_value("/kayscrawford", kayscrawford, .FALSE.)
            scalar(l)%kayscrawford = l_to_i(kayscrawford)

            ! Read boundary conditions for IB geometries
            SELECT TYPE(ib)
            TYPE IS (gc_t)
                nstl = SIZE(ib%stlnames)
                ALLOCATE(scalar(l)%geometries(nstl))
                DO n = 1, nstl
                    ! Read BC type
                    jsonptr = "/geometries/" // TRIM(ib%stlnames(n)) // "/type"
                    CALL sc%get_value(jsonptr, type, "flux")
                    SELECT CASE(TRIM(type))
                    CASE("value")
                        itype = 0
                    CASE("flux")
                        itype = 1
                    CASE DEFAULT
                        WRITE(*, *) "Invalid type: ", type
                        CALL errr(__FILE__, __LINE__)
                    END SELECT
                    scalar(l)%geometries(n)%flag = itype

                    ! Read BC value
                    jsonptr = "/geometries/" // TRIM(ib%stlnames(n)) // "/value"
                    CALL sc%get_value(jsonptr, rvalue, 0.0)
                    scalar(l)%geometries(n)%value = rvalue
                END DO
            END SELECT

            ! Read source terms
            IF (sc%exists("/sources")) THEN
                CALL sc%get_size("/sources", nsource)
            ELSE
                nsource = 0
            END IF
            ALLOCATE(scalar(l)%sources(nsource))
            DO n = 1, nsource
                WRITE(jsonptr, '("/sources/", I0)') n-1
                CALL sc%get(source, jsonptr)

                CALL source%get_value("/value", scalar(l)%sources(n)%value)
                CALL source%get_value("/field", scalar(l)%sources(n)%field, "")

                CALL source%finish()
            END DO

            CALL sc%finish()
        END DO

        ! Optional values
        CALL scaconf%get_value("/solve", solve_scalar, .TRUE.)
        CALL scaconf%get_value("/prturb", prturb, 1.0)

        IF (myid == 0) THEN
            WRITE(*, '("SCALAR TRANSPORT:")')
            WRITE(*, '(2X, "solve:                ", L1)') solve_scalar
            WRITE(*, '(2X, "nsca:                 ", I0)') nsca
        END IF

        ! Declare fields
        DO l = 1, nsca
            CALL set_field(scalar(l)%name, units=scalar(l)%units, dread=dread, &
                required=dcont, dwrite=dwrite, buffers=.TRUE.)

            ! Get field, set PRMOL and scalar index as attribute
            CALL get_field(t, scalar(l)%name)
            CALL t%set_attr(scalar(l)%prmol, "PRMOL")
            CALL t%set_attr(l, "SCAIDX")

            ! For RK time integration
            CALL set_field("D"//TRIM(scalar(l)%name))

            ! To store previous time value
            CALL set_field(TRIM(scalar(l)%name)//"_OLD")

            ! Scalar initial condition
            WRITE(jsonptr, '("/scalars/", I0, "/initial_value")') l-1
            IF (scaconf%exists(jsonptr)) THEN
                IF (scaconf%is_char(jsonptr)) THEN
                    CALL scaconf%get_value(jsonptr, initial_expr)
                    CALL set_initial_expr(t, initial_expr)
                ELSE
                    CALL scaconf%get_value(jsonptr, rvalue)
                    CALL set_initial_value(t, rvalue)
                END IF

                CALL get_field(t_old, TRIM(scalar(l)%name)//"_OLD")
                t_old%arr = t%arr
            ELSE
                CALL scaconf%set_value(jsonptr, 0.0_realk)
            END IF

            IF (myid == 0) THEN
                WRITE(*, '(2X, "Scalar:               ", A, " prmol: ", G0, ' &
                    //'"   unit: [", 7I3, " ]")') &
                    scalar(l)%name, scalar(l)%prmol, scalar(l)%units
            END IF
        END DO

        IF (myid == 0) THEN
            WRITE(*, '()')
        END IF

        ! Compute BT field from BP
        CALL set_field("BT", dwrite=.TRUE.)
        CALL get_field(bt, "BT")
        CALL blockbt(bt)
    END SUBROUTINE init_scacore


    SUBROUTINE finish_scacore
        IF (ALLOCATED(scalar)) DEALLOCATE(scalar)
    END SUBROUTINE finish_scacore


    PURE ELEMENTAL REAL(realk) FUNCTION prt(this, gtgmol)
        ! Calculation of local turbulent Prandtl number
        !
        ! limitations of Kays/Crawford: 0.5 < prmol < 7,
        !                               any Re
        !                               any dp/dx

        ! Kays/Crawford provides a relatively high prturb
        ! near the wall (in the sublayer) but approaches
        ! 0.85 as y+ increases into the log layer

        ! function arguments
        CLASS(scalar_t), INTENT(in) :: this
        REAL(realk), INTENT(IN) :: gtgmol

        ! Local variables
        REAL(realk) :: kayscrawford

        ! Kays/Crawford solution
        IF (this%kayscrawford == 0) THEN
            prt = prturb
        ELSE
            IF (gtgmol > 0.0) THEN
                kayscrawford = 0.5882 + 0.228*gtgmol &
                    - 0.0441*gtgmol**2*(1.0 - exp(-5.165/gtgmol))
            ELSE
                kayscrawford = this%prmol
            END IF
            prt = kayscrawford
        END IF
    END FUNCTION prt


    SUBROUTINE maskbt(t_f)
        ! Subroutine arguments
        TYPE(field_t), INTENT(inout) :: t_f

        ! Local variables
        INTEGER(intk) :: i, igrid
        INTEGER(intk) :: kk, jj, ii
        TYPE(field_t), POINTER :: bt_f
        REAL(realk), POINTER, CONTIGUOUS :: t(:, :, :), bt(:, :, :)

        CALL get_field(bt_f, "BT")

        DO i = 1, nmygrids
            igrid = mygrids(i)
            CALL get_mgdims(kk, jj, ii, igrid)

            CALL bt_f%get_ptr(bt, igrid)
            CALL t_f%get_ptr(t, igrid)

            CALL maskbt_grid(kk, jj, ii, t, bt)
        END DO
    END SUBROUTINE maskbt


    PURE SUBROUTINE maskbt_grid(kk, jj, ii, t, bt)
        ! Subroutine arguments
        INTEGER(intk), INTENT(in) :: kk, jj, ii
        REAL(realk), INTENT(inout) :: t(kk, jj, ii)
        REAL(realk), INTENT(in) :: bt(kk, jj, ii)

        ! Local variables
        INTEGER :: k, j, i

        ! TODO: Indices?
        DO i = 3, ii-2
            DO j = 3, jj-2
                DO k = 3, kk-2
                    t(k, j, i) = t(k, j, i)*bt(k, j, i)
                END DO
            END DO
        END DO
    END SUBROUTINE maskbt_grid


    SUBROUTINE set_initial_value(sca_f, val)
        ! Subroutine arguments
        TYPE(field_t), INTENT(inout) :: sca_f
        REAL(realk), INTENT(in) :: val

        ! Local variables
        INTEGER(intk) :: igrid
        REAL(realk), POINTER, CONTIGUOUS :: t(:, :, :)

        DO igrid = 1, nmygrids
            CALL sca_f%get_ptr(t, igrid)
            t = val
        END DO
    END SUBROUTINE set_initial_value


    SUBROUTINE set_initial_expr(sca_f, expr)
        ! Subroutine arguments
        TYPE(field_t), INTENT(inout) :: sca_f
        CHARACTER(LEN=*), INTENT(IN) :: expr

        ! Local variables
        INTEGER(intk) :: i, igrid
        REAL(realk), PARAMETER :: zero = 0.0_realk
        REAL(realk), POINTER, CONTIGUOUS :: t(:, :, :)
        REAL(realk), POINTER, CONTIGUOUS :: x(:), y(:), z(:)
        REAL(realk), POINTER, CONTIGUOUS :: dx(:), dy(:), dz(:)
        REAL(realk), POINTER, CONTIGUOUS :: ddx(:), ddy(:), ddz(:)

        DO i = 1, nmygrids
            igrid = mygrids(i)

            CALL sca_f%get_ptr(t, igrid)

            CALL get_fieldptr(x, "X", igrid)
            CALL get_fieldptr(y, "Y", igrid)
            CALL get_fieldptr(z, "Z", igrid)

            CALL get_fieldptr(dx, "DX", igrid)
            CALL get_fieldptr(dy, "DY", igrid)
            CALL get_fieldptr(dz, "DZ", igrid)

            CALL get_fieldptr(ddx, "DDX", igrid)
            CALL get_fieldptr(ddy, "DDY", igrid)
            CALL get_fieldptr(ddz, "DDZ", igrid)

            ! rho, gmol, tu_level, timeph are zero
            ! This allows the "t" field to be set as an expression of spatial
            ! coordinates only - not other fields.
            CALL initial_condition(t, "t", expr, zero, zero, zero, zero, &
                x, y, z, dx, dy, dz, ddx, ddy, ddz)
        END DO
    END SUBROUTINE set_initial_expr

END MODULE scacore_mod
