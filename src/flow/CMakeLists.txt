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
    gc_flowstencils_wmxpolquad.cxx
    gc_flowstencils_wmxpolquadfcorr.cxx
    sip_classic_mod.F90
    sip_hyperplane_mod.F90
    bound_flow_kernels.c
    laplacephi_mod.F90
    sor_mod.F90
    itinfo_mod.F90
    lesmodel_mod.F90
    plog_mod.F90
    pressuresolver_mod.F90
    setboundarybuffers_mod.F90
    timeintegration_mod.F90
    tstle4/tstle1_mod.F90
    tstle4/tstle2_mod.F90
    tstle4/tstle2_convection_kernels.cxx
    tstle4/tstle2_diffusion_kernels.cxx
    tstle4/tstle2_parent_kernels.cxx
    tstle4/tstle2_pressure_gradient_kernels.cxx
    tstle4/tstle4_mod.F90
    wernerwengle_mod.F90
)

target_include_directories(flow
    PRIVATE ${PROJECT_SOURCE_DIR}/src/core
    PRIVATE ${MPI_Fortran_INCLUDE_PATH}
    PRIVATE ${HDF5_Fortran_INCLUDE_DIRS}
)

target_link_libraries(flow
    PRIVATE core
    PRIVATE ib
    PRIVATE ${MPI_Fortran_LIBRARIES}
    PRIVATE $<$<BOOL:${MGLET_OFFLOAD}>:mglet_openmp>
)

set_source_files_properties(pressuresolver_mod.F90 PROPERTIES COMPILE_FLAGS
    $<$<STREQUAL:"${CMAKE_Fortran_COMPILER_ID}","GNU">:-Wno-unused-function>
)

