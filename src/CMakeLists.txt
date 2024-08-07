# Add separate libraries/subdirectories
add_subdirectory("core")
add_subdirectory("flow")
add_subdirectory("ib")
add_subdirectory("plugins")
add_subdirectory("scalar")

add_library(mgletlib STATIC
    runinfo_mod.F90
    timeloop_mod.F90
)
set_source_files_properties(timeloop_mod.F90 PROPERTIES COMPILE_FLAGS
    $<$<STREQUAL:"${CMAKE_Fortran_COMPILER_ID}","GNU">:-fall-intrinsics>
)
target_include_directories(mgletlib
    PRIVATE ${MPI_Fortran_INCLUDE_PATH}
    PRIVATE ${HDF5_Fortran_INCLUDE_DIRS}
)
target_link_libraries(mgletlib
    PUBLIC core
    PUBLIC flow
    PUBLIC ib
    PUBLIC plugins
    PUBLIC scalar
    PRIVATE ${MPI_Fortran_LIBRARIES}
)

# Executable
add_executable(mglet
    mglet.F90
)
target_include_directories(mglet
    PRIVATE ${MPI_Fortran_INCLUDE_PATH}
    PRIVATE ${HDF5_Fortran_INCLUDE_DIRS}
)
target_link_libraries(mglet
    PUBLIC mgletlib
    PRIVATE ${HDF5_Fortran_LIBRARIES}
    PRIVATE ${MPI_Fortran_LIBRARIES}
)

set_property(TARGET mglet PROPERTY LINKER_LANGUAGE Fortran)
