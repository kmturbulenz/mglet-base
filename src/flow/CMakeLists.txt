add_library(flow STATIC
    bound_flow_mod.F90
    bound_pressure_mod.F90
    boussinesqterm_mod.F90
    coriolisterm_mod.F90
    flowcore_mod.F90
    flow_mod.F90
    flowstat_mod.F90
    gc_compbodyforce_mod.F90
    gc_flowstencils_mod.F90
    sip_classic_mod.F90
    sip_hyperplane_mod.F90
    tstle4_kernels.c
    laplacephi_mod.F90
    sor_mod.F90
    itinfo_mod.F90
    lesmodel_mod.F90
    plog_mod.F90
    pressuresolver_mod.F90
    setboundarybuffers_mod.F90
    timeintegration_mod.F90
    tstle4_mod.F90
    wernerwengle_mod.F90
)

target_include_directories(flow
    PRIVATE ${MPI_Fortran_INCLUDE_PATH}
    PRIVATE ${HDF5_Fortran_INCLUDE_DIRS}
)

target_link_libraries(flow
    PRIVATE core
    PRIVATE ib
    PRIVATE ${MPI_Fortran_LIBRARIES}
    PRIVATE $<$<BOOL:${MGLET_OFFLOAD}>:OpenMP::OpenMP_Fortran>
    PRIVATE $<$<BOOL:${MGLET_OFFLOAD}>:OpenMP::OpenMP_C>
)

target_compile_options(flow
    PRIVATE "$<$<AND:$<BOOL:${MGLET_OFFLOAD}>,$<COMPILE_LANGUAGE:C>>:${MGLET_OFFLOAD_COMPILE_FLAGS}>"
    PRIVATE "$<$<AND:$<BOOL:${MGLET_OFFLOAD}>,$<COMPILE_LANGUAGE:C>>:${MGLET_OFFLOAD_ARCH_FLAGS}>"
)

set_source_files_properties(tstle4_kernels.c PROPERTIES COMPILE_FLAGS
    $<$<STREQUAL:"${CMAKE_C_COMPILER_ID}","GNU">:-Wno-unknown-pragmas>
)

set_source_files_properties(pressuresolver_mod.F90 PROPERTIES COMPILE_FLAGS
    $<$<STREQUAL:"${CMAKE_Fortran_COMPILER_ID}","GNU">:-Wno-unused-function>
)
