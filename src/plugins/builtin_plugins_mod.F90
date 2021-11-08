MODULE builtin_plugins_mod
    USE core_mod

    IMPLICIT NONE(type, external)
    PRIVATE

    PUBLIC :: register_builtin_plugins

CONTAINS
    SUBROUTINE register_builtin_plugins()
        USE derivfields_mod
        USE probes_mod
        USE snapshots_mod
        USE uvwbulk_mod

        CALL register_plugin("DERIVFIELDS", init_late=init_derivfields, &
            postprocess=calc_derivfields, finish=finish_derivfields)
        CALL register_plugin("PROBES", init_late=init_probes, &
            postprocess=sample_probes, finish=finish_probes)
        CALL register_plugin("SNAPSHOTS", init_late=init_snapshots, &
            postprocess=sample_snapshots, finish=finish_snapshots)
        CALL register_plugin("UVWBULK", init_late=init_uvwbulk, &
            itinfo=itinfo_uvwbulk, finish=finish_uvwbulk)
    END SUBROUTINE register_builtin_plugins
END MODULE builtin_plugins_mod
