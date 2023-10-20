MODULE itinfo_scalar_mod
    USE, INTRINSIC :: ISO_C_BINDING, ONLY: C_F_POINTER, C_PTR
    USE core_mod
    USE MPI_f08
    USE scacore_mod

    IMPLICIT NONE
    PRIVATE

    INTEGER(intk), PARAMETER :: itinfo_nsca_max = 16

    TYPE, BIND(C) :: itinfo_scalar_t
        INTEGER(c_intk) :: level
        REAL(c_realk) :: volume
        INTEGER(c_intk) :: exploded

        ! Scalar values
        REAL(c_realk) :: tmean(itinfo_nsca_max)
        REAL(c_realk) :: tmeansqr(itinfo_nsca_max)
    END TYPE itinfo_scalar_t

    INTEGER(intk), PARAMETER :: itinfo_elems = 5

    TYPE(MPI_Datatype) :: mpitype
    TYPE(MPI_Op) :: mpiop

    TYPE(itinfo_scalar_t), ALLOCATABLE :: itinfo_grd(:)
    TYPE(itinfo_scalar_t), ALLOCATABLE :: itinfo_level(:)

    CHARACTER(len=*), PARAMETER :: logfile = logdir//"/scalar.log"

    PUBLIC :: init_itinfo_scalar, finish_itinfo_scalar, itinfo_scalar_sample, &
        itinfo_scalar_print

CONTAINS
    ! Tag the grids, create MPI datatype and reduction operator
    SUBROUTINE init_itinfo_scalar(dcont)
        ! Subroutine arguments
        LOGICAL, INTENT(in) :: dcont

        ! Local variables
        INTEGER(intk) :: i
        INTEGER(int32) :: blocklen(itinfo_elems)
        TYPE(MPI_Datatype) :: types(itinfo_elems)
        TYPE(MPI_Datatype) :: nscaltype
        INTEGER(MPI_ADDRESS_KIND) :: base, disp(itinfo_elems)
        TYPE(itinfo_scalar_t) :: foo

        IF (nsca > itinfo_nsca_max) CALL errr(__FILE__, __LINE__)

        ALLOCATE(itinfo_grd(nmygrids))
        ALLOCATE(itinfo_level(minlevel:maxlevel))

        CALL zero_itinfo(itinfo_grd)
        CALL zero_itinfo(itinfo_level)

        ! Temporary type for conting. lists
        CALL MPI_Type_contiguous(itinfo_nsca_max, mglet_mpi_real, nscaltype)

        ! Create user-defined MPI reduction and data types
        blocklen = 1

        CALL MPI_Get_address(foo%level, disp(1))
        CALL MPI_Get_address(foo%volume, disp(2))
        CALL MPI_Get_address(foo%exploded, disp(3))
        CALL MPI_Get_address(foo%tmean, disp(4))
        CALL MPI_Get_address(foo%tmeansqr, disp(5))

        base = disp(1)
        DO i = 1, itinfo_elems
            disp(i) = disp(i) - base
        END DO

        types(1) = mglet_mpi_int
        types(2) = mglet_mpi_real
        types(3) = mglet_mpi_int
        types(4) = nscaltype
        types(5) = nscaltype

        CALL MPI_Type_create_struct(itinfo_elems, blocklen, disp, types, &
            mpitype)
        CALL MPI_Type_commit(mpitype)

        ! Free types used temporarily
        !   OK according to:
        !   https://www.open-mpi.org/doc/v4.0/man3/MPI_Type_free.3.php
        CALL MPI_Type_free(nscaltype)

        CALL MPI_Op_create(itinfo_reduce_mpi, .TRUE., mpiop)

        ! Initialize logs
        CALL itinfo_init_logs(dcont)
    END SUBROUTINE init_itinfo_scalar


    SUBROUTINE finish_itinfo_scalar()
        DEALLOCATE(itinfo_grd)
        DEALLOCATE(itinfo_level)

        CALL MPI_Op_free(mpiop)

        CALL MPI_Type_free(mpitype)
    END SUBROUTINE finish_itinfo_scalar


    ! Initialize log files
    SUBROUTINE itinfo_init_logs(dcont)
        ! Subroutine arguments
        LOGICAL, INTENT(in) :: dcont

        ! Local variables
        LOGICAL :: exists
        INTEGER :: logunit

        IF (myid == 0) THEN
            INQUIRE(FILE=logfile, EXIST=exists)
            IF (.NOT. exists .OR. .NOT. dcont) THEN
                OPEN(NEWUNIT=logunit, FILE=logfile)
                WRITE(logunit, '(A9, A18, A6, A6, 2A14)') &
                    "#  ITTOT", "TIME", "LEVEL", "T", "TMEAN", "TMEANSQR"
                CLOSE(logunit)
            END IF
        END IF
    END SUBROUTINE itinfo_init_logs


    SUBROUTINE itinfo_reduce(ivec, iovec, length, datatype)
        ! Subroutine arguments
        INTEGER(int32), INTENT(in) :: length
        TYPE(itinfo_scalar_t), INTENT(in) :: ivec(length)
        TYPE(itinfo_scalar_t), INTENT(inout) :: iovec(length)
        TYPE(MPI_Datatype), INTENT(in) :: datatype

        ! Local variables
        INTEGER(intk) :: i

        DO i = 1, length
            IF (iovec(i)%level /= ivec(i)%level &
                    .AND. iovec(i)%level /= HUGE(1_intk)) THEN
                CALL errr(__FILE__, __LINE__)
            END IF

            iovec(i)%exploded = IOR(iovec(i)%exploded, ivec(i)%exploded)

            IF (iovec(i)%volume + ivec(i)%volume > 0.0) THEN
                iovec(i)%tmean = (iovec(i)%tmean*iovec(i)%volume + &
                    ivec(i)%tmean*ivec(i)%volume) &
                    /(iovec(i)%volume + ivec(i)%volume)

                iovec(i)%tmeansqr = (iovec(i)%tmeansqr*iovec(i)%volume + &
                    ivec(i)%tmeansqr*ivec(i)%volume) &
                    /(iovec(i)%volume + ivec(i)%volume)
            ELSE
                iovec(i)%tmean = 0.0
                iovec(i)%tmeansqr = 0.0
            END IF

            ! It's very important that this update last for the
            ! average calculations to be correct
            iovec(i)%volume = iovec(i)%volume + ivec(i)%volume
        END DO
    END SUBROUTINE itinfo_reduce


    SUBROUTINE itinfo_reduce_mpi(invec, inoutvec, length, datatype)
        ! Subroutine arguments
        TYPE(C_PTR), VALUE :: invec
        TYPE(C_PTR), VALUE :: inoutvec
        INTEGER(int32) :: length        ! Declaring INTENT makes it incompatible
        TYPE(MPI_Datatype) :: datatype  ! with MPI_f08 MPI_User_function

        ! Local variables
        TYPE(itinfo_scalar_t), POINTER :: idata(:)
        TYPE(itinfo_scalar_t), POINTER :: iodata(:)

        CALL C_F_POINTER(invec, idata, [length])
        CALL C_F_POINTER(inoutvec, iodata, [length])

        CALL itinfo_reduce(idata, iodata, length, datatype)
    END SUBROUTINE itinfo_reduce_mpi


    SUBROUTINE itinfo_scalar_sample(igrid, tmean, tmeansqr)

        USE, INTRINSIC :: IEEE_ARITHMETIC, ONLY: IEEE_IS_FINITE

        ! Subroutine arguments
        INTEGER(intk), INTENT(in) :: igrid
        REAL(realk), INTENT(in) :: tmean(:)
        REAL(realk), INTENT(in) :: tmeansqr(:)

        ! Local variables
        INTEGER(intk) :: i, l
        REAL(realk) :: volume

        ! Sanity checks
        IF (SIZE(tmean) /= nsca) CALL errr(__FILE__, __LINE__)
        IF (SIZE(tmeansqr) /= nsca) CALL errr(__FILE__, __LINE__)

        CALL get_imygrid(i, igrid)
        CALL get_gridvolume(volume, igrid)

        itinfo_grd(i)%level = level(igrid)
        itinfo_grd(i)%volume = volume

        itinfo_grd(i)%tmean(1:nsca) = tmean(:)
        itinfo_grd(i)%tmeansqr(1:nsca) = tmeansqr(:)

        itinfo_grd(i)%exploded = 0

        DO l = 1, nsca
            IF (.NOT. IEEE_IS_FINITE(tmeansqr(l))) THEN
                itinfo_grd(i)%exploded = 1
            END IF
        END DO
    END SUBROUTINE itinfo_scalar_sample


    SUBROUTINE itinfo_scalar_print(itstep, ittot, timeph, exploded)
        ! Subroutine arguments
        INTEGER(intk), INTENT(in) :: itstep
        INTEGER(intk), INTENT(in) :: ittot
        REAL(realk), INTENT(in) :: timeph
        INTEGER(intk), INTENT(inout) :: exploded

        ! Local variables
        INTEGER(intk) :: i, l, igrid, nlevels
        INTEGER(int32) :: ierr
        INTEGER :: logunit
        TYPE(itinfo_scalar_t), ALLOCATABLE :: itinfo_level_tmp(:)

        ALLOCATE(itinfo_level_tmp(minlevel:maxlevel))
        CALL zero_itinfo(itinfo_level_tmp)
        CALL zero_itinfo(itinfo_level)

        ! First do a local reduction
        DO i = 1, nmygrids
            igrid = mygrids(i)
            CALL itinfo_reduce(itinfo_grd(i), itinfo_level_tmp(level(igrid)), &
                1, mpitype)
        END DO

        ! Then reduce across all processes
        nlevels = maxlevel - minlevel + 1
        CALL MPI_Allreduce(itinfo_level_tmp, itinfo_level, nlevels, &
            mpitype, mpiop, MPI_COMM_WORLD, ierr)

        ! Check if solution is exploded
        DO i = minlevel, maxlevel
            exploded = MAX(exploded, itinfo_level(i)%exploded)
        END DO

        IF (myid == 0) THEN
            OPEN(NEWUNIT=logunit, FILE=logfile, POSITION="APPEND")
            DO i = minlevel, maxlevel
                DO l = 1, nsca
                    WRITE(logunit, '(I9, ES18.10, I6, I6, 2ES14.5)') &
                        ittot, timeph, i, l, itinfo_level(i)%tmean(l), &
                        itinfo_level(i)%tmeansqr(l)
                END DO
            END DO
            CLOSE(logunit)
        END IF

        ! For good measure...
        DEALLOCATE(itinfo_level_tmp)
        CALL zero_itinfo(itinfo_grd)
    END SUBROUTINE itinfo_scalar_print


    PURE SUBROUTINE zero_itinfo(itinfo)
        ! Subroutine arguments
        TYPE(itinfo_scalar_t), INTENT(out) :: itinfo(:)

        ! Local variables
        INTEGER(intk) :: i

        DO i = LBOUND(itinfo, 1), UBOUND(itinfo, 1)
            itinfo(i)%level = HUGE(1_intk)
            itinfo(i)%volume = 0.0
            itinfo(i)%exploded = 0

            itinfo(i)%tmean = 0.0
            itinfo(i)%tmeansqr = 0.0
        END DO
    END SUBROUTINE zero_itinfo
END MODULE itinfo_scalar_mod
