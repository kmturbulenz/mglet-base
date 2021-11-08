MODULE itinfo_mod
    USE, INTRINSIC :: ISO_C_BINDING, ONLY: C_F_POINTER, C_PTR
    USE core_mod
    USE MPI_f08

    IMPLICIT NONE
    PRIVATE

    TYPE, BIND(C) :: itinfo_t
        INTEGER(c_intk) :: level
        INTEGER(c_intk) :: ngrids
        REAL(c_realk) :: volume
        INTEGER(c_intk) :: ipc

        ! Information to be propagated to all processes, currently this
        ! is first if the solution is exploded and secondly a custom
        ! message (unused).
        INTEGER(c_intk) :: exploded
        INTEGER(c_intk) :: message

        ! Maximum divergence and its location
        REAL(c_realk) :: divmax
        REAL(c_realk) :: divmax_x
        REAL(c_realk) :: divmax_y
        REAL(c_realk) :: divmax_z
        INTEGER(c_intk) :: divmax_grid

        ! Maximum CFL number and its location
        REAL(c_realk) :: cflmax
        REAL(c_realk) :: cflmax_x
        REAL(c_realk) :: cflmax_y
        REAL(c_realk) :: cflmax_z
        INTEGER(c_intk) :: cflmax_grid

        ! Kinetic energy
        REAL(c_realk) :: esumg
        REAL(c_realk) :: esums
    END TYPE itinfo_t

    INTEGER(intk), PARAMETER :: itinfo_elems = 18

    TYPE(MPI_Datatype) :: mpitype
    TYPE(MPI_Op) :: mpiop

    TYPE(itinfo_t), ALLOCATABLE :: itinfo_grd(:)
    TYPE(itinfo_t), ALLOCATABLE :: itinfo_level(:)

    CHARACTER(len=*), PARAMETER :: generalfile = logdir//"/general.log"

    PUBLIC :: init_itinfo, finish_itinfo, itinfo_sample, &
        itinfo_print, itinfo_get_cflmax

CONTAINS
    ! Tag the grids, create MPI datatype and reduction operator
    SUBROUTINE init_itinfo(dcont)
        ! Subroutine arguments
        LOGICAL, INTENT(in) :: dcont

        ! Local variables
        INTEGER(intk) :: i
        INTEGER(int32) :: blocklen(itinfo_elems)
        TYPE(MPI_Datatype) :: types(itinfo_elems)
        INTEGER(MPI_ADDRESS_KIND) :: base, disp(itinfo_elems)
        TYPE(itinfo_t) :: foo

        ALLOCATE(itinfo_grd(nmygrids))
        ALLOCATE(itinfo_level(minlevel:maxlevel))

        CALL zero_itinfo(itinfo_grd)
        CALL zero_itinfo(itinfo_level)

        ! Create user-defined MPI reduction and data types
        blocklen = 1

        CALL MPI_Get_address(foo%level, disp(1))
        CALL MPI_Get_address(foo%ngrids, disp(2))
        CALL MPI_Get_address(foo%volume, disp(3))
        CALL MPI_Get_address(foo%ipc, disp(4))

        CALL MPI_Get_address(foo%exploded, disp(5))
        CALL MPI_Get_address(foo%message, disp(6))

        CALL MPI_Get_address(foo%divmax, disp(7))
        CALL MPI_Get_address(foo%divmax_x, disp(8))
        CALL MPI_Get_address(foo%divmax_y, disp(9))
        CALL MPI_Get_address(foo%divmax_z, disp(10))
        CALL MPI_Get_address(foo%divmax_grid, disp(11))

        CALL MPI_Get_address(foo%cflmax, disp(12))
        CALL MPI_Get_address(foo%cflmax_x, disp(13))
        CALL MPI_Get_address(foo%cflmax_y, disp(14))
        CALL MPI_Get_address(foo%cflmax_z, disp(15))
        CALL MPI_Get_address(foo%cflmax_grid, disp(16))

        CALL MPI_Get_address(foo%esumg, disp(17))
        CALL MPI_Get_address(foo%esums, disp(18))

        base = disp(1)
        DO i = 1, itinfo_elems
            disp(i) = disp(i) - base
        END DO

        types(1) = mglet_mpi_int
        types(2) = mglet_mpi_int
        types(3) = mglet_mpi_real
        types(4) = mglet_mpi_int

        types(5) = mglet_mpi_int
        types(6) = mglet_mpi_int

        types(7) = mglet_mpi_real
        types(8) = mglet_mpi_real
        types(9) = mglet_mpi_real
        types(10) = mglet_mpi_real
        types(11) = mglet_mpi_int

        types(12) = mglet_mpi_real
        types(13) = mglet_mpi_real
        types(14) = mglet_mpi_real
        types(15) = mglet_mpi_real
        types(16) = mglet_mpi_int

        types(17) = mglet_mpi_real
        types(18) = mglet_mpi_real

        CALL MPI_Type_create_struct(itinfo_elems, blocklen, disp, types, &
            mpitype)
        CALL MPI_Type_commit(mpitype)

        CALL MPI_Op_create(itinfo_reduce_mpi, .TRUE., mpiop)

        ! Initialize logs
        CALL itinfo_init_logs(dcont)
    END SUBROUTINE init_itinfo


    SUBROUTINE finish_itinfo()
        DEALLOCATE(itinfo_grd)
        DEALLOCATE(itinfo_level)

        CALL MPI_Op_free(mpiop)

        CALL MPI_Type_free(mpitype)
    END SUBROUTINE finish_itinfo


    ! Initialize log files
    SUBROUTINE itinfo_init_logs(dcont)
        ! Subroutine arguments
        LOGICAL, INTENT(in) :: dcont

        ! Local variables
        LOGICAL :: exists
        INTEGER :: logunit

        IF (myid == 0) THEN
            INQUIRE(FILE=generalfile, EXIST=exists)
            IF (.NOT. exists .OR. .NOT. dcont) THEN
                OPEN(NEWUNIT=logunit, FILE=generalfile)
                WRITE(logunit, '(A9, A18, A6, A4, 14A14)') &
                    "#  ITTOT", "TIME", "LEVEL", "IPC", &
                    "DIVMAX", "DIVMAX_X", "DIVMAX_Y", "DIVMAX_Z", &
                    "DIVMAX_GRID", &
                    "CFLMAX", "CFLMAX_X", "CFLMAX_Y", "CFLMAX_Z", &
                    "CFLMAX_GRID", &
                    "ESUMG", "ESUMS"
                CLOSE(logunit)
            END IF
        END IF
    END SUBROUTINE itinfo_init_logs


    SUBROUTINE itinfo_reduce(ivec, iovec, length, datatype)
        ! Subroutine arguments
        INTEGER(int32), INTENT(in) :: length
        TYPE(itinfo_t), INTENT(in) :: ivec(length)
        TYPE(itinfo_t), INTENT(inout) :: iovec(length)
        TYPE(MPI_Datatype), INTENT(in) :: datatype

        ! Local variables
        INTEGER(intk) :: i

        DO i = 1, length
            IF (iovec(i)%level /= ivec(i)%level &
                    .AND. iovec(i)%level /= HUGE(1_intk)) THEN
                CALL errr(__FILE__, __LINE__)
            END IF

            iovec(i)%ipc = MAX(iovec(i)%ipc, ivec(i)%ipc)

            iovec(i)%exploded = MAX(iovec(i)%exploded, ivec(i)%exploded)
            iovec(i)%message = MAX(iovec(i)%message, ivec(i)%message)

            IF (ABS(iovec(i)%divmax) < ABS(ivec(i)%divmax)) THEN
                iovec(i)%divmax = ivec(i)%divmax
                iovec(i)%divmax_x = ivec(i)%divmax_x
                iovec(i)%divmax_y = ivec(i)%divmax_y
                iovec(i)%divmax_z = ivec(i)%divmax_z
                iovec(i)%divmax_grid = ivec(i)%divmax_grid
            END IF

            IF (ABS(iovec(i)%cflmax) < ABS(ivec(i)%cflmax)) THEN
                iovec(i)%cflmax = ivec(i)%cflmax
                iovec(i)%cflmax_x = ivec(i)%cflmax_x
                iovec(i)%cflmax_y = ivec(i)%cflmax_y
                iovec(i)%cflmax_z = ivec(i)%cflmax_z
                iovec(i)%cflmax_grid = ivec(i)%cflmax_grid
            END IF

            IF (iovec(i)%volume + ivec(i)%volume > 0.0) THEN
                iovec(i)%esumg = (iovec(i)%esumg*iovec(i)%volume + &
                    ivec(i)%esumg*ivec(i)%volume) &
                    /(iovec(i)%volume + ivec(i)%volume)

                iovec(i)%esums = (iovec(i)%esums*iovec(i)%volume + &
                    ivec(i)%esums*ivec(i)%volume) &
                    /(iovec(i)%volume + ivec(i)%volume)
            ELSE
                iovec(i)%esumg = 0.0
                iovec(i)%esums = 0.0
            END IF

            ! It's very important that this update last for the
            ! average calculations to be correct
            iovec(i)%ngrids = iovec(i)%ngrids + ivec(i)%ngrids
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
        TYPE(itinfo_t), POINTER :: idata(:)
        TYPE(itinfo_t), POINTER :: iodata(:)

        CALL C_F_POINTER(invec, idata, [length])
        CALL C_F_POINTER(inoutvec, iodata, [length])

        CALL itinfo_reduce(idata, iodata, length, datatype)
    END SUBROUTINE itinfo_reduce_mpi


    SUBROUTINE itinfo_sample(igrid, divmax, divmax_pos, cflmax, cflmax_pos, &
            esumg, esums, ipc)

        USE, INTRINSIC :: IEEE_ARITHMETIC, ONLY: IEEE_IS_FINITE

        ! Subroutine arguments
        INTEGER(intk), INTENT(in) :: igrid
        REAL(realk), INTENT(in), OPTIONAL :: divmax
        REAL(realk), INTENT(in), OPTIONAL :: divmax_pos(3)
        REAL(realk), INTENT(in), OPTIONAL :: cflmax
        REAL(realk), INTENT(in), OPTIONAL :: cflmax_pos(3)
        REAL(realk), INTENT(in), OPTIONAL :: esumg
        REAL(realk), INTENT(in), OPTIONAL :: esums
        INTEGER(intk), INTENT(in), OPTIONAL :: ipc

        ! Local variables
        INTEGER(intk) :: i
        REAL(realk) :: volume

        ! Sanity checks of optional arguments
        IF (PRESENT(divmax) .NEQV. PRESENT(divmax_pos)) THEN
            CALL errr(__FILE__, __LINE__)
        END IF
        IF (PRESENT(cflmax) .NEQV. PRESENT(cflmax_pos)) THEN
            CALL errr(__FILE__, __LINE__)
        END IF
        IF (PRESENT(esumg) .NEQV. PRESENT(esums)) THEN
            CALL errr(__FILE__, __LINE__)
        END IF

        CALL get_imygrid(i, igrid)
        CALL get_gridvolume(volume, igrid)

        itinfo_grd(i)%level = level(igrid)
        itinfo_grd(i)%ngrids = 1
        itinfo_grd(i)%volume = volume

        IF (PRESENT(ipc)) THEN
            itinfo_grd(i)%ipc = ipc
        END IF

        IF (PRESENT(divmax)) THEN
            itinfo_grd(i)%divmax = divmax
            itinfo_grd(i)%divmax_x = divmax_pos(1)
            itinfo_grd(i)%divmax_y = divmax_pos(2)
            itinfo_grd(i)%divmax_z = divmax_pos(3)
            itinfo_grd(i)%divmax_grid = igrid

            IF (.NOT. IEEE_IS_FINITE(divmax)) THEN
                itinfo_grd(i)%exploded = 1
            END IF
        END IF

        IF (PRESENT(cflmax)) THEN
            itinfo_grd(i)%cflmax = cflmax
            itinfo_grd(i)%cflmax_x = cflmax_pos(1)
            itinfo_grd(i)%cflmax_y = cflmax_pos(2)
            itinfo_grd(i)%cflmax_z = cflmax_pos(3)
            itinfo_grd(i)%cflmax_grid = igrid

            IF (.NOT. IEEE_IS_FINITE(cflmax)) THEN
                itinfo_grd(i)%exploded = 1
            END IF
        END IF

        IF (PRESENT(esumg)) THEN
            itinfo_grd(i)%esumg = esumg
            itinfo_grd(i)%esums = esums

            IF (.NOT. IEEE_IS_FINITE(esumg)) THEN
                itinfo_grd(i)%exploded = 1
            END IF
        END IF
    END SUBROUTINE itinfo_sample


    SUBROUTINE itinfo_print(itstep, ittot, timeph, globalcflmax, exploded)
        ! Subroutine arguments
        INTEGER(intk), INTENT(in) :: itstep
        INTEGER(intk), INTENT(in) :: ittot
        REAL(realk), INTENT(in) :: timeph
        REAL(realk), INTENT(out) :: globalcflmax
        INTEGER(intk), INTENT(out) :: exploded

        ! Local variables
        INTEGER(intk) :: i, igrid, nlevels
        INTEGER(int32) :: ierr
        INTEGER :: logunit
        TYPE(itinfo_t), ALLOCATABLE :: itinfo_level_tmp(:)

        ! Initialize INTENT(out)
        globalcflmax = 0.0
        exploded = 0

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
            globalcflmax = MAX(globalcflmax, itinfo_level(i)%cflmax)
            exploded = MAX(exploded, itinfo_level(i)%exploded)
        END DO

        IF (myid == 0) THEN
            ! Print summary to screen
            DO i = minlevel, maxlevel
                WRITE(*, '("LEVEL:", I3, " IPC= ", I3, " DIVMAX= ", ' // &
                    'ES12.4, " CFLMAX= ", F9.5)') &
                    i, itinfo_level(i)%ipc, itinfo_level(i)%divmax, &
                    itinfo_level(i)%cflmax
            END DO

            OPEN(NEWUNIT=logunit, FILE=generalfile, POSITION="APPEND")
            DO i = minlevel, maxlevel
                WRITE(logunit, '(I9, ES18.10, I6, I4, ES14.5, F14.5, ' // &
                    'F14.5, F14.5, I14, ES14.5, F14.5, F14.5, F14.5, ' // &
                    'I14, ES14.5, ES14.5)') &
                    ittot, timeph, i, itinfo_level(i)%ipc, &
                    itinfo_level(i)%divmax, itinfo_level(i)%divmax_x, &
                    itinfo_level(i)%divmax_y, itinfo_level(i)%divmax_z, &
                    itinfo_level(i)%divmax_grid, &
                    itinfo_level(i)%cflmax, itinfo_level(i)%cflmax_x, &
                    itinfo_level(i)%cflmax_y, itinfo_level(i)%cflmax_z, &
                    itinfo_level(i)%cflmax_grid, &
                    itinfo_level(i)%esumg, itinfo_level(i)%esums
            END DO
            CLOSE(logunit)
        END IF

        ! For good measure...
        DEALLOCATE(itinfo_level_tmp)
        CALL zero_itinfo(itinfo_grd)
    END SUBROUTINE itinfo_print


    PURE SUBROUTINE zero_itinfo(itinfo)
        ! Subroutine arguments
        TYPE(itinfo_t), INTENT(out) :: itinfo(:)

        ! Local variables
        INTEGER(intk) :: i

        DO i = LBOUND(itinfo, 1), UBOUND(itinfo, 1)
            itinfo(i)%level = HUGE(1_intk)
            itinfo(i)%ngrids = 0
            itinfo(i)%volume = 0.0
            itinfo(i)%ipc = 0

            itinfo(i)%exploded = 0
            itinfo(i)%message = 0

            itinfo(i)%divmax = 0.0
            itinfo(i)%divmax_x = 0.0
            itinfo(i)%divmax_y = 0.0
            itinfo(i)%divmax_z = 0.0
            itinfo(i)%divmax_grid = 0

            itinfo(i)%cflmax = 0.0
            itinfo(i)%cflmax_x = 0.0
            itinfo(i)%cflmax_y = 0.0
            itinfo(i)%cflmax_z = 0.0
            itinfo(i)%cflmax_grid = 0

            itinfo(i)%esumg = 0.0
            itinfo(i)%esums = 0.0
        END DO
    END SUBROUTINE zero_itinfo


    SUBROUTINE itinfo_get_cflmax(cflmax)
        ! Subroutine arguments
        REAL(realk), INTENT(out):: cflmax

        ! Local variables
        INTEGER(intk) :: i

        cflmax = 0.0

        DO i = minlevel, maxlevel
           cflmax = MAX(cflmax, itinfo_level(i)%cflmax)
        END DO
    END SUBROUTINE itinfo_get_cflmax
END MODULE itinfo_mod
