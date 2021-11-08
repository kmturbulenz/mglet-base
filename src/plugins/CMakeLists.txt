add_library(plugins STATIC
    builtin_plugins_mod.F90
    derivfields_mod.F90
    probes_mod.F90
    snapshots_mod.F90
    uvwbulk_mod.F90
)

target_include_directories(plugins
    PRIVATE ${MPI_Fortran_INCLUDE_PATH}
    PRIVATE ${HDF5_Fortran_INCLUDE_DIRS}
)

target_link_libraries(plugins
    PRIVATE core
    PRIVATE ib
    PRIVATE ${MPI_Fortran_LIBRARIES}
    PRIVATE ${HDF5_Fortran_LIBRARIES}
)
