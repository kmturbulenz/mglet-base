add_library(backend STATIC)

set_property(TARGET backend PROPERTY CXX_STANDARD 17)

set(GPU_SOURCES
    gpu/kernels_1d.cu
)

if(MGLET_OFFLOAD_BACKEND STREQUAL "cuda")
    set_source_files_properties(${GPU_SOURCES} PROPERTIES LANGUAGE CUDA)
    target_sources(backend PRIVATE ${GPU_SOURCES})
elseif(MGLET_OFFLOAD_BACKEND STREQUAL "hip")
    set_source_files_properties(${GPU_SOURCES} PROPERTIES LANGUAGE HIP)
    target_sources(backend PRIVATE ${GPU_SOURCES})
else()
    target_sources(backend PRIVATE empty.cpp)
endif()

target_include_directories(backend
    PUBLIC ${CMAKE_CURRENT_SOURCE_DIR}
)

target_link_libraries(backend
    PRIVATE $<$<BOOL:${MGLET_OFFLOAD}>:OpenMP::OpenMP_Fortran>
    PRIVATE $<$<STREQUAL:${MGLET_PROFILE_ANNOTATIONS},roctx>:rocprofiler-sdk-roctx::rocprofiler-sdk-roctx>
    PRIVATE $<$<STREQUAL:${MGLET_PROFILE_ANNOTATIONS},nvtx>:CUDA::nvtx3>
    PRIVATE $<$<STREQUAL:${MGLET_OFFLOAD_BACKEND},cuda>:CUDA::cudart>
    PRIVATE $<$<STREQUAL:${MGLET_OFFLOAD_BACKEND},hip>:hip::device>
)
