PROGRAM main
    USE builtin_plugins_mod, ONLY: register_builtin_plugins
    USE core_mod, ONLY: init_core, finish_core, fields_begin_read, &
        fields_begin_write, fields_write, fields_end_rw, init_plugins, &
        finish_plugins
    USE flow_mod, ONLY: init_flow, finish_flow
    USE ib_mod, ONLY: init_ib, finish_ib, ib
    USE timeloop_mod, ONLY: init_timeloop, finish_timeloop, timeloop
    USE scalar_mod, ONLY: init_scalar, finish_scalar

    ! Initialization of core data structures
    CALL init_core()

    ! Open result file for reading
    CALL fields_begin_read()

    ! Immersed boundary
    CALL init_ib()
    CALL ib%blockbp()
    CALL ib%read_stencils()
    CALL ib%giteig()

    ! Initialize plugins
    CALL register_builtin_plugins()
    CALL init_plugins()

    ! Initialize builtin physical models
    CALL init_flow()
    CALL init_scalar()

    ! This initialize the time loop. Reads the RUNINFO table in case of DCONT.
    CALL init_timeloop()

    ! After this position no data is allowed to be read from fields.h5 any more
    CALL fields_end_rw()

    ! Run time loop
    CALL timeloop()

    ! Writes RUNINFO-table
    CALL fields_begin_write()
    CALL finish_timeloop()

    ! Finish plugins
    CALL finish_plugins()

    ! Finish physical models
    CALL finish_scalar()
    CALL finish_flow()
    CALL finish_ib()

    ! Write out fields and close fields.h5
    CALL fields_write()
    CALL fields_end_rw()

    ! Deallocate core consructs
    CALL finish_core()
END PROGRAM main
