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
    USE qsort_mod
    USE readstl_mod
    USE rungekutta_mod
    USE shmem_mod
    USE simdfunctions_mod
    USE statistics_mod
    USE stencilio_mod
    USE sym_name_mod
    USE tensormath_mod
    USE timekeeper_mod
    USE timer_mod
    USE utils_mod
    USE write3d_mod

    IMPLICIT NONE
    ! In this module everything is public. It exports all the public
    ! definitions from the other modules as well, which makes it easier
    ! to use them, you do not need to "USE" so many modules each place.

    CHARACTER(len=*), PARAMETER :: logdir = "LOGS"

CONTAINS
    ! Call all initializer-functions in the core
    SUBROUTINE init_core()
        USE HDF5, ONLY: h5open_f
        USE MPI_f08, ONLY: MPI_Init

        ! Local variables
        INTEGER(int32) :: ierr

        CALL MPI_Init()
        CALL h5open_f(ierr)

        CALL init_comms()
        CALL init_precision()
        CALL init_buildinfo()
        CALL init_timer()

        CALL set_timer(1, "MGLET")
        CALL start_timer(1)

        CALL init_hdf5common()
        CALL init_fort7()
        CALL init_dlopen()
        CALL init_grids()
        CALL init_pointers()
        CALL init_commbuf()
        CALL init_fieldio()
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
        CALL finish_fieldio()
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

END MODULE core_mod
