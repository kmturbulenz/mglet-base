add_library(backends STATIC
    gpuruntime_mod.F90
    profilemarker_mod.F90
)
target_include_directories(backends PRIVATE ${CMAKE_CURRENT_SOURCE_DIR})

if(MGLET_GPU_BACKEND STREQUAL "cuda")
    find_package(CUDAToolkit REQUIRED)
    if(MGLET_ENABLE_GPU_RUNTIME)
        enable_language(CUDA)
        target_link_libraries(backends PRIVATE CUDA::cudart)
        target_sources(backends PRIVATE
            cuda/gpu_backend_cuda_runtime.cxx
        )
    endif()
    if(MGLET_PROFILE_MARKERS)
        target_link_libraries(backends PRIVATE CUDA::nvtx3)
        target_sources(backends PRIVATE cuda/gpu_backend_cuda_marker.cxx)
    endif()
elseif(MGLET_GPU_BACKEND STREQUAL "rocm")
    if(MGLET_ENABLE_GPU_RUNTIME)
        find_package(hip REQUIRED)
        enable_language(HIP)
        target_link_libraries(backends PRIVATE hip::host)
        target_sources(backends PRIVATE
            hip/gpu_backend_hip_runtime.cxx
        )
    endif()
    if(MGLET_PROFILE_MARKERS)
        find_package(rocprofiler-sdk-roctx REQUIRED)
        target_link_libraries(backends PRIVATE rocprofiler-sdk-roctx::rocprofiler-sdk-roctx)
        target_sources(backends PRIVATE hip/gpu_backend_hip_marker.cxx)
    endif()
elseif(MGLET_GPU_BACKEND STREQUAL "none")
    target_sources(backends PRIVATE
        none/gpu_backend_none_marker.cxx
        none/gpu_backend_none_runtime.cxx
    )
    if(MGLET_ENABLE_GPU_RUNTIME)
        message(WARNING "GPU runtime 'none' selected but MGLET_ENABLE_GPU_RUNTIME is ON. No GPU runtime will be used in the backend.")
    endif()
    if(MGLET_PROFILE_MARKERS)
        message(WARNING "Profiling markers enabled but no GPU backend selected. Markers will not be functional.")
    endif()
else()
    message(FATAL_ERROR "Invalid GPU backend selected: '${MGLET_GPU_BACKEND}'.")
endif()
