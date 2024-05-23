MODULE checkblock_mod
    USE MPI_f08
    USE HDF5

    USE core_mod

    IMPLICIT NONE(type, external)
    PRIVATE

    PUBLIC :: checkblock

CONTAINS
    SUBROUTINE checkblock(bp_f, outfile)
        TYPE(field_t), INTENT(in) :: bp_f
        CHARACTER(len=*), INTENT(in) :: outfile

        ! Local variables
        REAL(realk), POINTER, CONTIGUOUS :: bp(:, :, :)
        INTEGER(intk) :: n1
        INTEGER(intk) :: igr, igrid
        INTEGER(intk) :: i, j, k, ii, jj, kk
        INTEGER(intk), ALLOCATABLE :: nfluidcells(:)

        ALLOCATE(nfluidcells(ngrid))
        nfluidcells = 0

        DO igr = 1, nmygrids
            igrid = mygrids(igr)

            CALL bp_f%get_ptr(bp, igrid)
            CALL get_mgdims(kk, jj, ii, igrid)

            ! 1: Grid contains fluid?
            n1 = 0
            DO i = 1, ii
                DO j = 1, jj
                    DO k = 1, kk
                        n1 = n1 + NINT(bp(k, j, i), intk)
                    END DO
                END DO
            END DO
            IF (n1 > 0) nfluidcells(igrid) = IBSET(nfluidcells(igrid), 0)

            ! 2: Front contains fluid?
            n1 = 0
            DO i = 2, 3
                DO j = 1, jj
                    DO k = 1, kk
                        n1 = n1 + NINT(bp(k, j, i), intk)
                    END DO
                END DO
            END DO
            IF (n1 > 0) nfluidcells(igrid) = IBSET(nfluidcells(igrid), 1)

            ! 4: Back contains fluid?
            n1 = 0
            DO i = ii-2, ii-1
                DO j = 1, jj
                    DO k = 1, kk
                        n1 = n1 + NINT(bp(k, j, i), intk)
                    END DO
                END DO
            END DO
            IF (n1 > 0) nfluidcells(igrid) = IBSET(nfluidcells(igrid), 2)

            ! 8: Right contains fluid?
            n1 = 0
            DO i = 1, ii
                DO j = 2, 3
                    DO k = 1, kk
                        n1 = n1 + NINT(bp(k, j, i), intk)
                    END DO
                END DO
            END DO
            IF (n1 > 0) nfluidcells(igrid) = IBSET(nfluidcells(igrid), 3)

            ! 16: Left contains fluid?
            n1 = 0
            DO i = 1, ii
                DO j = jj-2, jj-1
                    DO k = 1, kk
                        n1 = n1 + NINT(bp(k, j, i), intk)
                    END DO
                END DO
            END DO
            IF (n1 > 0) nfluidcells(igrid) = IBSET(nfluidcells(igrid), 4)

            ! 32: Bottom contains fluid?
            n1 = 0
            DO i = 1, ii
                DO j = 1, jj
                    DO k = 2, 3
                        n1 = n1 + NINT(bp(k, j, i), intk)
                    END DO
                END DO
            END DO
            IF (n1 > 0) nfluidcells(igrid) = IBSET(nfluidcells(igrid), 5)

            ! 64: Top contains fluid?
            n1 = 0
            DO i = 1, ii
                DO j = 1, jj
                    DO k = kk-2, kk-1
                        n1 = n1 + NINT(bp(k, j, i), intk)
                    END DO
                END DO
            END DO
            IF (n1 > 0) nfluidcells(igrid) = IBSET(nfluidcells(igrid), 6)
        END DO

        ! All ranks must know which grids to write out, even though only rank 0
        ! actually write the data
        CALL MPI_Allreduce(MPI_IN_PLACE, nfluidcells, ngrid, mglet_mpi_int, &
            mpi_sum, MPI_COMM_WORLD)

        CALL write_grids(nfluidcells, outfile)
    END SUBROUTINE checkblock
END MODULE checkblock_mod
