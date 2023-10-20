add_library(flow STATIC
    bound_flow_mod.F90
    boussinesqterm_mod.F90
    flowcore_mod.F90
    flow_mod.F90
    flowstat_mod.F90
    gc_flowstencils_mod.F90
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
)

target_link_libraries(flow
    PRIVATE core
    PRIVATE ib
    PRIVATE ${MPI_Fortran_LIBRARIES}
)

set_source_files_properties(pressuresolver_mod.F90 PROPERTIES COMPILE_FLAGS
    $<$<STREQUAL:"${CMAKE_Fortran_COMPILER_ID}","GNU">:-Wno-unused-function>
)
