add_library(core STATIC
    bound_mod.F90
    buildinfo_mod.F90
    charfunc_mod.F90
    checksum_mod.F90
    commbuf_mod.F90
    comms_mod.F90
    config_mod.F90
    connect2_mod.F90
    corefields_mod.F90
    core_mod.F90
    create_directory.c
    dlopen_mod.F90
    envvars_mod.F90
    err_mod.F90
    expression_mod.F90
    exprtk_wrapper.cxx
    fieldio2_mod.F90
    field_t/basefield_mod.F90
    field_t/field_mod.F90
    field_t/intfield_mod.F90
    field_t/realfield_mod.F90
    fields_mod.F90
    fort7_mod.F90
    gridio_mod.F90
    grids_mod.F90
    hdf5common_mod.F90
    jsoncpp_wrapper.cxx
    mgletmath_mod.F90
    plugins_mod.F90
    pointers_mod.F90
    precision_mod.F90
    pvtk.c
    pvtk_mod.F90
    qsort_mod.F90
    readstl_mod.F90
    rungekutta_mod.F90
    shmem_mod.F90
    simdfunctions_mod.F90
    sort_wrapper.cxx
    statistics_mod.F90
    stencilio_mod.F90
    stl_is_binary.c
    sym_name_mod.F90
    tensormath_mod.F90
    timekeeper_mod.F90
    timer_mod.F90
    utils_mod.F90
    write3d_mod.F90
)

target_include_directories(core
    PRIVATE ${MPI_Fortran_INCLUDE_PATH}
    PRIVATE ${HDF5_Fortran_INCLUDE_DIRS}
)

target_link_libraries(core
    PRIVATE ${MPI_Fortran_LIBRARIES}
    PRIVATE ${HDF5_Fortran_LIBRARIES}
    PRIVATE ZLIB::ZLIB
)

# Set include directories for specific files
set_source_files_properties(exprtk_wrapper.cxx PROPERTIES INCLUDE_DIRECTORIES
    ${exprtk_SOURCE_DIR}
)
set_source_files_properties(jsoncpp_wrapper.cxx PROPERTIES INCLUDE_DIRECTORIES
    ${json_SOURCE_DIR}/include
)

set_property(TARGET core PROPERTY CXX_STANDARD 11)

# Other neccesary compile flags
set_source_files_properties(comms_mod.F90 PROPERTIES COMPILE_FLAGS
    $<$<STREQUAL:"${CMAKE_Fortran_COMPILER_ID}","GNU">:-fall-intrinsics>
)

set_source_files_properties(connect2_mod.F90 PROPERTIES COMPILE_FLAGS
    $<$<STREQUAL:"${CMAKE_Fortran_COMPILER_ID}","GNU">:-Wno-error=function-elimination>
)

set_source_files_properties(err_mod.F90 PROPERTIES COMPILE_FLAGS
    $<$<STREQUAL:"${CMAKE_Fortran_COMPILER_ID}","GNU">:-fall-intrinsics>
)

set_source_files_properties(simdfunctions_mod.F90 PROPERTIES COMPILE_FLAGS
    $<$<STREQUAL:"${CMAKE_Fortran_COMPILER_ID}","GNU">:-Wno-compare-reals>
)

set_source_files_properties(stl_is_binary.c PROPERTIES COMPILE_FLAGS
    $<$<STREQUAL:"${CMAKE_C_COMPILER_ID}","GNU">:-Wno-unused-result>
)

set_source_files_properties(utils_mod.F90 PROPERTIES COMPILE_FLAGS
    $<$<STREQUAL:"${CMAKE_Fortran_COMPILER_ID}","IntelLLVM">:"-nostand -warn nostderrors">
)
