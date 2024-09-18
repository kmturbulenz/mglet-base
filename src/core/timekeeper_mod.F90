! Timekeeper module
!
! The purpose of this module is to keep exact track of time, even in single
! precision versions of the code, in which roundoff-errors become very large
! when running huge number of timesteps. Running many timesteps and
! keep track of time by the "traditional" 'time = time + dt' is not ideal,
! since the ratio between 'time' and 'dt' is equal to the number
! of timesteps. When running for example one million timesteps, the ratio
! dt/time = 1e-6, which is (too) close to single precsision accuracy.
!
! This module keeps track of the time in double precision, independently
! on the presion of the rest of MGLET. In addition to that, to really
! reduce roundoff errors, the module use Kahan's summation algorithm
! to make sure that the roundoff errors are bounded and independent on
! the number of timesteps run.
MODULE timekeeper_mod
    USE precision_mod, ONLY: real32, real64

    IMPLICIT NONE(type, external)
    PRIVATE

    TYPE :: timekeeper_t
        ! Time value in double precision
        REAL(real64), PRIVATE :: time = 0.0_real64

        ! Kahan summation coefficients
        REAL(real64), PRIVATE :: c = 0.0_real64
        REAL(real64), PRIVATE :: y = 0.0_real64
        REAL(real64), PRIVATE :: t = 0.0_real64
    CONTAINS
        GENERIC, PUBLIC :: init => init_sp, init_dp
        PROCEDURE, PRIVATE :: init_sp, init_dp

        GENERIC, PUBLIC :: get_time => get_time_sp, get_time_dp
        PROCEDURE, PRIVATE :: get_time_sp, get_time_dp

        GENERIC, PUBLIC :: add_to_time => add_to_time_sp, add_to_time_dp
        PROCEDURE, PRIVATE :: add_to_time_sp, add_to_time_dp

        GENERIC, PUBLIC :: tmpadd_to_time => tmpadd_to_time_sp, &
            tmpadd_to_time_dp
        PROCEDURE, PRIVATE :: tmpadd_to_time_sp, tmpadd_to_time_dp
    END TYPE timekeeper_t

    ! Public data items
    PUBLIC :: timekeeper_t

CONTAINS
    SUBROUTINE init_sp(this, starttime)
        CLASS(timekeeper_t), INTENT(out) :: this
        REAL(real32), INTENT(in) :: starttime

        this%time = starttime
        this%c = 0.0
        this%y = 0.0
        this%t = 0.0
    END SUBROUTINE init_sp


    SUBROUTINE init_dp(this, starttime)
        CLASS(timekeeper_t), INTENT(out) :: this
        REAL(real64), INTENT(in) :: starttime

        this%time = starttime
        this%c = 0.0
        this%y = 0.0
        this%t = 0.0
    END SUBROUTINE init_dp


    ! Return current time value
    SUBROUTINE get_time_sp(this, timeph)
        CLASS(timekeeper_t), INTENT(inout) :: this
        REAL(real32), INTENT(out) :: timeph
        timeph = REAL(this%time, real32)
    END SUBROUTINE get_time_sp


    ! Return current time value
    SUBROUTINE get_time_dp(this, timeph)
        CLASS(timekeeper_t), INTENT(inout) :: this
        REAL(real64), INTENT(out) :: timeph
        timeph = this%time
    END SUBROUTINE get_time_dp


    ! Adds a small value dt to time and update global time
    SUBROUTINE add_to_time_dp(this, dt)
        CLASS(timekeeper_t), INTENT(inout) :: this
        REAL(real64), INTENT(in) :: dt

        ! Kahan's summation algorithm
        this%y = dt - this%c
        this%t = this%time + this%y
        this%c = (this%t - this%time) - this%y
        this%time = this%t
    END SUBROUTINE add_to_time_dp


    ! Adds a small value dt to time and update global time
    SUBROUTINE add_to_time_sp(this, dt)
        CLASS(timekeeper_t), INTENT(inout) :: this
        REAL(real32), INTENT(in) :: dt

        CALL this%add_to_time(REAL(dt, real64))
    END SUBROUTINE add_to_time_sp


    ! Adds a small value dt to time, without updating global time
    SUBROUTINE tmpadd_to_time_dp(this, dt, timerk)
        CLASS(timekeeper_t), INTENT(inout) :: this
        REAL(real64), INTENT(in) :: dt
        REAL(real64), INTENT(out) :: timerk

        ! Since we are only updating the time value temporarily, and
        ! not store the result in the module, we make local copied of
        ! of the Kahan's summation coefficients here in this function
        REAL(real64) :: ly
        REAL(real64) :: lt
        REAL(real64) :: lc
        REAL(real64) :: ltime

        ly = this%y
        lt = this%t
        lc = this%c
        ltime = this%time

        ! Kahan's summation algorithm
        ly = dt - lc
        lt = ltime + ly
        lc = (lt - ltime) - ly  ! Strictly not neeccesary since we do
                                ! not keep this lc variable
        ltime = lt

        ! Return time
        timerk = ltime
    END SUBROUTINE tmpadd_to_time_dp


    SUBROUTINE tmpadd_to_time_sp(this, dt, timerk)
        CLASS(timekeeper_t), INTENT(inout) :: this
        REAL(real32), INTENT(in) :: dt
        REAL(real32), INTENT(out) :: timerk

        REAL(real64) :: timerk_dp

        CALL this%tmpadd_to_time(REAL(dt, real64), timerk_dp)
        timerk = REAL(timerk_dp, real32)
    END SUBROUTINE tmpadd_to_time_sp
END MODULE timekeeper_mod
