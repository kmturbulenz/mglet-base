MODULE write3d_mod
    USE precision_mod
    USE grids_mod
    USE field_mod
    USE fields_mod

    IMPLICIT NONE(type, external)
    PRIVATE

    PUBLIC :: writevtk

CONTAINS
    SUBROUTINE writevtk(field1, field2, field3, field4, prefix)
        ! Subroutine arguments
        TYPE(field_t), INTENT(in) :: field1
        TYPE(field_t), INTENT(in), OPTIONAL :: field2
        TYPE(field_t), INTENT(in), OPTIONAL :: field3
        TYPE(field_t), INTENT(in), OPTIONAL :: field4
        CHARACTER(len=*), INTENT(in), OPTIONAL :: prefix

        ! Local variables
        INTEGER :: unit
        CHARACTER(len=mglet_filename_max) :: filename
        INTEGER(intk) :: igr, igrid
        INTEGER(intk) :: kk, jj, ii
        INTEGER(intk) :: k, j, i
        REAL(realk), POINTER, CONTIGUOUS :: x(:), y(:), z(:)
        REAL(realk), POINTER, CONTIGUOUS :: dx(:), dy(:), dz(:)
        REAL(realk), POINTER, CONTIGUOUS :: phi(:, :, :)
        REAL(realk), ALLOCATABLE :: xstag(:), ystag(:), zstag(:)

        DO igr = 1, nmygrids
            igrid = mygrids(igr)
            CALL get_mgdims(kk, jj, ii, igrid)

            CALL get_fieldptr(x, "X", igrid)
            CALL get_fieldptr(y, "Y", igrid)
            CALL get_fieldptr(z, "Z", igrid)

            CALL get_fieldptr(dx, "DX", igrid)
            CALL get_fieldptr(dy, "DY", igrid)
            CALL get_fieldptr(dz, "DZ", igrid)

            ALLOCATE(xstag(ii+1))
            xstag(1) = x(1) - dx(1)/2.0
            DO i = 1, ii
                xstag(i+1) = x(i) + 0.5*dx(i)
            END DO

            ALLOCATE(ystag(jj+1))
            ystag(1) = y(1) - dy(1)/2.0
            DO j = 1, jj
                ystag(j+1) = y(j) + 0.5*dy(j)
            END DO

            ALLOCATE(zstag(kk+1))
            zstag(1) = z(1) - dz(1)/2.0
            DO k = 1, kk
                zstag(k+1) = z(k) + 0.5*dz(k)
            END DO

            IF (PRESENT(prefix)) THEN
                WRITE(filename, '(A, "-igrid-", I0, ".vtk")') &
                    TRIM(prefix), igrid
            ELSE
                WRITE(filename, '("igrid-", I0, ".vtk")') igrid
            END IF
            OPEN(newunit=unit, file=filename)

            WRITE(unit, '("# vtk DataFile Version 3.0")')
            WRITE(unit, '("MGLET SUBROUTINE writevtk")')
            WRITE(unit, '("ASCII")')

            WRITE(unit, '("DATASET RECTILINEAR_GRID")')
            WRITE(unit, '("DIMENSIONS ", I0, 1X, I0, 1X, I0)') ii+1, jj+1, kk+1
            WRITE(unit, '("X_COORDINATES ", i0, " float")') ii+1
            WRITE(unit, '(E14.7)') (xstag(i), i = 1, ii+1)
            WRITE(unit, '("Y_COORDINATES ", I0, " float")') jj+1
            WRITE(unit, '(E14.7)') (ystag(j), j = 1, jj+1)
            WRITE(unit, '("Z_COORDINATES ", I0, " float")') kk+1
            WRITE(unit, '(E14.7)') (zstag(k), k = 1, kk+1)

            WRITE(unit, '("CELL_DATA ", I0)') ii*jj*kk

            CALL field1%get_ptr(phi, igrid)
            WRITE(unit, '("SCALARS ", A, " float 1")') TRIM(field1%name)
            WRITE(unit, '("LOOKUP_TABLE default")')
            WRITE(unit, '(E14.7)') &
                (((phi(k, j, i), i = 1, ii), j = 1, jj), k = 1, kk)

            IF (PRESENT(field2)) THEN
                CALL field2%get_ptr(phi, igrid)
                WRITE(unit, '("SCALARS ", A, " float 1")') TRIM(field2%name)
                WRITE(unit, '("LOOKUP_TABLE default")')
                WRITE(unit, '(E14.7)') &
                    (((phi(k, j, i), i = 1, ii), j = 1, jj), k = 1, kk)
            END IF

            IF (PRESENT(field3)) THEN
                CALL field3%get_ptr(phi, igrid)
                WRITE(unit, '("SCALARS ", A, " float 1")') TRIM(field3%name)
                WRITE(unit, '("LOOKUP_TABLE default")')
                WRITE(unit, '(E14.7)') &
                    (((phi(k, j, i), i = 1, ii), j = 1, jj), k = 1, kk)
            END IF

            IF (PRESENT(field4)) THEN
                CALL field4%get_ptr(phi, igrid)
                WRITE(unit, '("SCALARS ", A, " float 1")') TRIM(field4%name)
                WRITE(unit, '("LOOKUP_TABLE default")')
                WRITE(unit, '(E14.7)') &
                    (((phi(k, j, i), i = 1, ii), j = 1, jj), k = 1, kk)
            END IF

            CLOSE(unit)

            DEALLOCATE(xstag)
            DEALLOCATE(ystag)
            DEALLOCATE(zstag)
        END DO
    END SUBROUTINE writevtk

END MODULE write3d_mod
