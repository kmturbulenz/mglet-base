add_library(scalar STATIC
    blockbt_mod.F90
    bound_scalar_mod.F90
    gc_scastencils_mod.F90
    gc_scastencils.cxx
    itinfo_scalar_mod.F90
    scacore_mod.F90
    scalar_mod.F90
    scastat_mod.F90
    setboundarybuffers_scalar_mod.F90
    timeintegrate_scalar_mod.F90
)

target_include_directories(scalar
    PRIVATE ${MPI_Fortran_INCLUDE_PATH}
    PRIVATE ${HDF5_Fortran_INCLUDE_DIRS}
)

target_link_libraries(scalar
    PRIVATE core
    PRIVATE ib
    PRIVATE flow
    PRIVATE ${MPI_Fortran_LIBRARIES}
    PRIVATE $<$<BOOL:${MGLET_OFFLOAD}>:OpenMP::OpenMP_Fortran>
    PRIVATE $<$<BOOL:${MGLET_OFFLOAD}>:OpenMP::OpenMP_C>
    PRIVATE $<$<BOOL:${MGLET_OFFLOAD}>:OpenMP::OpenMP_CXX>
)

target_compile_options(scalar
    PRIVATE "$<$<AND:$<BOOL:${MGLET_OFFLOAD}>,$<COMPILE_LANGUAGE:C>>:${MGLET_OFFLOAD_COMPILE_FLAGS}>"
    PRIVATE "$<$<AND:$<BOOL:${MGLET_OFFLOAD}>,$<COMPILE_LANGUAGE:C>>:${MGLET_OFFLOAD_ARCH_FLAGS}>"
    PRIVATE "$<$<AND:$<BOOL:${MGLET_OFFLOAD}>,$<COMPILE_LANGUAGE:CXX>>:${MGLET_OFFLOAD_COMPILE_FLAGS}>"
    PRIVATE "$<$<AND:$<BOOL:${MGLET_OFFLOAD}>,$<COMPILE_LANGUAGE:CXX>>:${MGLET_OFFLOAD_ARCH_FLAGS}>"
)

# GCC complains about unknown pragmas if we compile without -fopenmp
set_source_files_properties(gc_scastencils.cxx PROPERTIES COMPILE_FLAGS
    $<$<STREQUAL:"${CMAKE_CXX_COMPILER_ID}","GNU">:-Wno-unknown-pragmas>
)
