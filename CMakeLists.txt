cmake_minimum_required(VERSION 3.21)
if(POLICY CMP0126)
    cmake_policy(SET CMP0126 NEW)
endif()
if(POLICY CMP0135)
	cmake_policy(SET CMP0135 NEW)
endif()

# Set the project name
project(mglet-base LANGUAGES C CXX Fortran)

include(FortranCInterface)
FortranCInterface_VERIFY(CXX)

find_package(ZLIB REQUIRED)

find_package(MPI COMPONENTS C Fortran REQUIRED)
if(MPI_Fortran_HAVE_F08_MODULE)
    message(STATUS "MPI Library supports MPI_f08")
else()
    message(FATAL_ERROR "MPI Library does not support MPI_f08")
endif()

set(HDF5_PREFER_PARALLEL TRUE)
set(HDF5_USE_STATIC_LIBRARIES FALSE)
find_package(HDF5 COMPONENTS C Fortran REQUIRED)
if(HDF5_IS_PARALLEL)
    message(STATUS "HDF5 Library supports parallel I/O")
else()
    message(WARNING "HDF5 does not support parallel I/O")
endif()

# CMake build options
option(MGLET_REAL64 "Compile MGLET with double precision support" OFF)
if(MGLET_REAL64)
    message(STATUS "Compiling MGLET in double precision mode")
    add_compile_definitions(_MGLET_DOUBLE_PRECISION_=1)
else()
    message(STATUS "Compiling MGLET in single precision mode")
endif()

option(MGLET_IFK64 "Compile MGLET with 64 bit integer fields" ON)
if(MGLET_IFK64)
    message(STATUS "Compiling MGLET with 64 ifk fields")
    add_compile_definitions(_MGLET_IFK64_=1)
else()
    message(STATUS "Compiling MGLET with 32 bit ifk fields")
endif()

option(MGLET_INT64 "Compile MGLET with 64 bit integer support (DOES NOT WORK)" OFF)
if(MGLET_INT64)
    message(STATUS "Compiling MGLET with 64 bit integers")
    message(WARNING "64-bit integer support does not work ATM!")
    add_compile_definitions(_MGLET_INT64_=1)
else()
    message(STATUS "Compiling MGLET with 32 bit integers")
endif()

option(MGLET_OPENMP "Enable OpenMP (experimental)" OFF)
if(MGLET_OPENMP)
    message(STATUS "Enabling OpenMP")
    find_package(OpenMP REQUIRED)
endif()

# CMake build options
set(MGLET_ENVPREFIX "" CACHE STRING "Prefix for all MGLET environment variables")
if(NOT MGLET_ENVPREFIX STREQUAL "")
    add_compile_definitions(_MGLET_ENVPREFIX_="${MGLET_ENVPREFIX}")
endif ()

# All fortran modules are placed here
set(CMAKE_Fortran_MODULE_DIRECTORY ${CMAKE_BINARY_DIR}/modules)

set(CMAKE_POSITION_INDEPENDENT_CODE ON)

# Debug define directive
add_compile_definitions("$<$<CONFIG:Debug>:_MGLET_DEBUG_=1>")

set(MGLET_C_FLAGS "" CACHE STRING "MGLET C flags")
set(MGLET_CXX_FLAGS "" CACHE STRING "MGLET CXX flags")
set(MGLET_Fortran_FLAGS "" CACHE STRING "MGLET Fortran flags")

set(MGLET_C_FLAGS_RELEASE "" CACHE STRING "MGLET C flags (release)")
set(MGLET_CXX_FLAGS_RELEASE "" CACHE STRING "MGLET CXX flags (release)")
set(MGLET_Fortran_FLAGS_RELEASE "" CACHE STRING "MGLET Fortran flags (release)")

set(MGLET_C_FLAGS_DEBUG "" CACHE STRING "MGLET C flags (debug)")
set(MGLET_CXX_FLAGS_DEBUG "" CACHE STRING "MGLET CXX flags (debug)")
set(MGLET_Fortran_FLAGS_DEBUG "" CACHE STRING "MGLET Fortran flags (debug)")

add_compile_options("$<$<COMPILE_LANGUAGE:C>:${MGLET_C_FLAGS}>"
    "$<$<COMPILE_LANGUAGE:CXX>:${MGLET_CXX_FLAGS}>"
    "$<$<COMPILE_LANGUAGE:Fortran>:${MGLET_Fortran_FLAGS}>"

    "$<$<AND:$<COMPILE_LANGUAGE:C>,$<CONFIG:Release>>:${MGLET_C_FLAGS_RELEASE}>"
    "$<$<AND:$<COMPILE_LANGUAGE:CXX>,$<CONFIG:Release>>:${MGLET_CXX_FLAGS_RELEASE}>"
    "$<$<AND:$<COMPILE_LANGUAGE:Fortran>,$<CONFIG:Release>>:${MGLET_Fortran_FLAGS_RELEASE}>"

    "$<$<AND:$<COMPILE_LANGUAGE:C>,$<CONFIG:Debug>>:${MGLET_C_FLAGS_DEBUG}>"
    "$<$<AND:$<COMPILE_LANGUAGE:CXX>,$<CONFIG:Debug>>:${MGLET_CXX_FLAGS_DEBUG}>"
    "$<$<AND:$<COMPILE_LANGUAGE:Fortran>,$<CONFIG:Debug>>:${MGLET_Fortran_FLAGS_DEBUG}>"

    # The preprocessor macro __FILE__ is often used in MGLET. This can be a
    # long path. Up until Fortran 2018 there was a max line length of
    # 132 characters. A filename can easily overflow this length.
    # GCC allows to increase this length to infinity with this flag.
    # Fortran 2023 will have the maximum line length increased to 10000
    # (ten thousand) characters, so this is a problem that will disappear
    # as soon as a Fortran 2023 compiler is available.
    "$<$<AND:$<COMPILE_LANGUAGE:Fortran>,$<STREQUAL:${CMAKE_Fortran_COMPILER_ID},GNU>>:-ffree-line-length-none>"
    "$<$<AND:$<COMPILE_LANGUAGE:Fortran>,$<STREQUAL:${CMAKE_Fortran_COMPILER_ID},IntelLLVM>>:-diag-disable 5268>"
)

if (MGLET_OPENMP)
    add_compile_options("$<$<COMPILE_LANGUAGE:C>:${OpenMP_C_FLAGS}>"
        "$<$<COMPILE_LANGUAGE:CXX>:${OpenMP_CXX_FLAGS}>"
        "$<$<COMPILE_LANGUAGE:Fortran>:${OpenMP_Fortran_FLAGS}>"
    )
endif()

# nlohman JSON parser library
include(FetchContent)
FetchContent_Declare(json
    URL https://github.com/nlohmann/json/releases/download/v3.11.3/json.tar.xz
    URL_HASH SHA256=d6c65aca6b1ed68e7a182f4757257b107ae403032760ed6ef121c9d55e81757d
)
FetchContent_MakeAvailable(json)

# Exprtk expression parser
FetchContent_Declare(exprtk
    GIT_REPOSITORY https://github.com/kmturbulenz/exprtk
    GIT_TAG 48390febb53664cc3577d2329b63b45231037cb3
)
FetchContent_MakeAvailable(exprtk)

# Perform various compiler checks
add_subdirectory("cmake")
add_subdirectory("src")
add_subdirectory("tests")
