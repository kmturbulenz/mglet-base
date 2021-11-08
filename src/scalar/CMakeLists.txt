add_library(scalar STATIC
    scacore_mod.F90
    scalar_mod.F90
    timeintegrate_scalar_mod.F90
)

target_include_directories(scalar
    PRIVATE ${MPI_Fortran_INCLUDE_PATH}
)

target_link_libraries(scalar
    PRIVATE core
    PRIVATE ib
    PRIVATE ${MPI_Fortran_LIBRARIES}
)
