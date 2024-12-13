MODULE core_mod
    ! Should absolutely NOT 'USE' external modules like MPI and HDF5 here,
    ! because then they will get re-exported from this module!!!
    USE bound_mod
    USE buildinfo_mod
    USE charfunc_mod
    USE checksum_mod
    USE commbuf_mod
    USE comms_mod
    USE config_mod
    USE connect2_mod
    USE corefields_mod
    USE dlopen_mod
    USE envvars_mod
    USE err_mod
    USE expression_mod
    USE fieldio2_mod
    USE field_mod
    USE fields_mod
    USE fort7_mod
    USE grids_mod
    ! Should not be neccesary to export gridio_mod - since this is only used
    ! in core functions
    USE hdf5common_mod
    USE mgletmath_mod
    USE plugins_mod
    USE pointers_mod
    USE precision_mod
    USE pvtk_mod
    USE qsort_mod
    USE readstl_mod
    USE rungekutta_mod
    USE shmem_mod
    USE simdfunctions_mod
    USE statistics_mod
    USE stencilio_mod
    USE tensormath_mod
    USE timekeeper_mod
    USE timer_mod
    USE utils_mod
    USE write3d_mod
    USE write_grids_mod

    IMPLICIT NONE (type, external)
    ! In this module everything is public. It exports all the public
    ! definitions from the other modules as well, which makes it easier
    ! to use them, you do not need to "USE" so many modules each place.

    CHARACTER(len=*), PARAMETER :: logdir = "LOGS"

    ! These are not to be called from other places in the code and are PRIVATE
    PRIVATE :: set_fpe_traps, set_underflow_mode, print_ieee_modes

CONTAINS
    ! Call all initializer-functions in the core
    SUBROUTINE init_core()
        USE HDF5, ONLY: h5open_f
        USE MPI_f08, ONLY: MPI_Init
        USE, INTRINSIC :: IEEE_EXCEPTIONS

        ! Local variables
        INTEGER(int32) :: ierr
        LOGICAL :: saved_fpe_mode(SIZE(ieee_all))

        ! Fetch the IEEE halting modes and disable all of them temporarily for
        ! the MPI_Init and h5open_f calls. This is to avoid the MPI and HDF5
        ! libraries to trigger floating point exceptions (sometimes they do)
        CALL IEEE_GET_HALTING_MODE(IEEE_ALL, saved_fpe_mode)
        CALL IEEE_SET_HALTING_MODE(IEEE_ALL, .FALSE.)

        ! Initialize MPI and HDF5
        CALL MPI_Init()
        CALL h5open_f(ierr)

        ! Restore the IEEE halting modes that was selected when the application
        ! was compiled
        CALL IEEE_SET_FLAG(IEEE_ALL, .FALSE.)
        CALL IEEE_SET_HALTING_MODE(IEEE_ALL, saved_fpe_mode)

        ! Need to set the communicators first before set_fp_traps or
        ! set_underflow_mode since they depend on them
        CALL init_comms()

        ! Allow the user to set desired halting and underflow modes via
        ! environment variables
        CALL set_fpe_traps()
        CALL set_underflow_mode()
        CALL print_ieee_modes()

        ! Continue MGLET startup

        CALL init_precision()
        CALL init_buildinfo()
        CALL init_timer()

        CALL set_timer(1, "MGLET")
        CALL set_timer(100, "FIELDIO2")
        CALL set_timer(101, "FIELDIO2_GATHER")
        CALL set_timer(102, "FIELDIO2_WRITE")
        CALL set_timer(103, "FIELDIO2_READ")
        CALL set_timer(104, "FIELDIO2_SCATTER")

        CALL start_timer(1)

        CALL init_hdf5common()
        CALL init_fort7()
        CALL init_dlopen()
        CALL init_grids()
        CALL init_pointers()
        CALL init_commbuf()
        CALL init_fields()
        CALL init_corefields()
        CALL init_connect2()

        ! Create log directory
        CALL create_directory(logdir)
    END SUBROUTINE init_core


    ! Cleanup core
    SUBROUTINE finish_core
        ! The calls in here should be in the reverse order of the init_core
        ! init-calls - first initialized is last closed and so on.
        USE HDF5, ONLY: h5close_f
        USE MPI_f08, ONLY: MPI_Finalize

        ! Local variables
        INTEGER(int32) :: ierr

        CALL finish_connect2()
        CALL finish_corefields()
        CALL finish_fields()
        CALL finish_commbuf()
        CALL finish_pointers()
        CALL finish_grids()
        CALL finish_dlopen()
        CALL finish_fort7()
        CALL finish_hdf5common()

        CALL stop_timer(1)
        CALL finish_timer()
        CALL finish_buildinfo()
        CALL finish_precision()
        CALL finish_comms()

        CALL h5close_f(ierr)
        CALL MPI_Finalize()
    END SUBROUTINE finish_core


    ! Set floating point exceptions handling
    ! When the environment variable MGLET_FPE_TRAP is set to a string
    ! containing one or more of the following words, the corresponding
    ! floating point exception will be enabled:
    ! - divide_by_zero
    ! - inexact
    ! - invalid
    ! - overflow
    ! - underflow
    ! The default is to do nothing (i.e. the method chosen by the compiler is
    ! used). An emptry string (or absence of any of the keywords above) will
    ! disable all floating point exceptions.
    SUBROUTINE set_fpe_traps()
        USE, INTRINSIC :: IEEE_EXCEPTIONS

        ! Local variables
        CHARACTER(len=1024) :: mglet_fpe_trap
        LOGICAL :: found

        ! Get the FPE mode from the environment
        CALL getenv_char_coll(mglet_fpe_trap, "MGLET_FPE_TRAP", found=found)

        ! When the environment variable is not set, do nothing
        IF (.NOT. found) RETURN

        ! Start by setting all flags to zero
        CALL IEEE_SET_HALTING_MODE(IEEE_ALL, .FALSE.)

        ! Convert to lower case
        mglet_fpe_trap = LOWER(mglet_fpe_trap)

        ! Set the flags according to the environment variable
        IF (INDEX(mglet_fpe_trap, "divide_by_zero") > 0) THEN
            CALL IEEE_SET_HALTING_MODE(IEEE_DIVIDE_BY_ZERO, .TRUE.)
        END IF
        IF (INDEX(mglet_fpe_trap, "inexact") > 0) THEN
            CALL IEEE_SET_HALTING_MODE(IEEE_INEXACT, .TRUE.)
        END IF
        IF (INDEX(mglet_fpe_trap, "invalid") > 0) THEN
            CALL IEEE_SET_HALTING_MODE(IEEE_INVALID, .TRUE.)
        END IF
        IF (INDEX(mglet_fpe_trap, "overflow") > 0) THEN
            CALL IEEE_SET_HALTING_MODE(IEEE_OVERFLOW, .TRUE.)
        END IF
        IF (INDEX(mglet_fpe_trap, "underflow") > 0) THEN
            CALL IEEE_SET_HALTING_MODE(IEEE_UNDERFLOW, .TRUE.)
        END IF
    END SUBROUTINE set_fpe_traps


    ! Set floating point underflow mode
    ! When the environment variable MGLET_UNDERFLOW is set to "gradual", the
    ! underflow mode is set to gradual. When it is set to "abrupt", the
    ! underflow mode is set to abrupt. The default is to do nothing (i.e. the
    ! method chosen by the compiler is used).
    SUBROUTINE set_underflow_mode()
        USE, INTRINSIC :: IEEE_ARITHMETIC

        ! Local variables
        CHARACTER(len=1024) :: mglet_underflow
        LOGICAL :: found

        ! Get the underflow mode from the environment
        CALL getenv_char_coll(mglet_underflow, "MGLET_UNDERFLOW", found=found)

        ! When the environment variable is not set, do nothing
        IF (.NOT. found) RETURN

        ! Convert to lower case
        mglet_underflow = LOWER(mglet_underflow)

        ! Set the flags according to the environment variable
        SELECT CASE(TRIM(mglet_underflow))
! TODO: monitor when flang implements IEEE_SET_UNDERFLOW_MODE
#ifndef __flang__
        CASE ("gradual")
            CALL IEEE_SET_UNDERFLOW_MODE(.TRUE.)
        CASE ("abrupt")
            CALL IEEE_SET_UNDERFLOW_MODE(.FALSE.)
#endif
        CASE DEFAULT
            WRITE(*, *) "Invalid value for MGLET_UNDERFLOW: ", &
                TRIM(mglet_underflow)
            CALL errr(__FILE__, __LINE__)
        END SELECT
    END SUBROUTINE set_underflow_mode


    SUBROUTINE print_ieee_modes()
        USE, INTRINSIC :: IEEE_ARITHMETIC

        ! Local variables
        LOGICAL :: flag
        CHARACTER(len=1024) :: mglet_fpe_trap

        IF (myid /= 0) RETURN
        mglet_fpe_trap = ""

        WRITE(*, '("IEEE aritmetic flags:")')

! TODO: monitor when flang implements IEEE_GET_UNDERFLOW_MODE
#ifndef __flang__
        CALL IEEE_GET_UNDERFLOW_MODE(flag)
        IF (flag) THEN
            WRITE(*, '("    Underflow:     gradual")')
        ELSE
            WRITE(*, '("    Underflow:     abrupt")')
        END IF
#endif

        CALL IEEE_GET_HALTING_MODE(IEEE_DIVIDE_BY_ZERO, flag)
        IF (flag) THEN
            mglet_fpe_trap = "divide_by_zero"
        END IF

        CALL IEEE_GET_HALTING_MODE(IEEE_INEXACT, flag)
        IF (flag) THEN
            mglet_fpe_trap = TRIM(mglet_fpe_trap) // " inexact"
        END IF

        CALL IEEE_GET_HALTING_MODE(IEEE_INVALID, flag)
        IF (flag) THEN
            mglet_fpe_trap = TRIM(mglet_fpe_trap) // " invalid"
        END IF

        CALL IEEE_GET_HALTING_MODE(IEEE_OVERFLOW, flag)
        IF (flag) THEN
            mglet_fpe_trap = TRIM(mglet_fpe_trap) // " overflow"
        END IF

        CALL IEEE_GET_HALTING_MODE(IEEE_UNDERFLOW, flag)
        IF (flag) THEN
            mglet_fpe_trap = TRIM(mglet_fpe_trap) // " underflow"
        END IF

        IF (LEN_TRIM(mglet_fpe_trap) > 0) THEN
            WRITE(*, '("    FPE traps:     ", A)') TRIM(mglet_fpe_trap)
        ELSE
            WRITE(*, '("    FPE traps:     none")')
        END IF

        WRITE(*, '()')
    END SUBROUTINE print_ieee_modes

END MODULE core_mod
