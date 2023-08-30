add_library(scalar STATIC
    blockbt_mod.F90
    bound_scalar_mod.F90
    gc_scastencils_mod.F90
    itinfo_scalar_mod.F90
    scacore_mod.F90
    scalar_mod.F90
    scastat_mod.F90
    setboundarybuffers_scalar_mod.F90
    timeintegrate_scalar_mod.F90
)

target_include_directories(scalar
    PRIVATE ${MPI_Fortran_INCLUDE_PATH}
)

target_link_libraries(scalar
    PRIVATE core
    PRIVATE ib
    PRIVATE flow
    PRIVATE ${MPI_Fortran_LIBRARIES}
)
