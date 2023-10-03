# https://discourse.cmake.org/t/help-with-cmake-intel-mpi-and-gnu-fortran/377/2
# https://gitlab.kitware.com/cmake/cmake/-/merge_requests/4153/diffs?commit_id=5861c6d450937a38cd3c5ab9aa878427a50b069d
cmake_minimum_required(VERSION 3.18)

# set the project name
project(mglet-base LANGUAGES C CXX Fortran)

include(FortranCInterface)
FortranCInterface_VERIFY(CXX)

find_package(ZLIB REQUIRED)

find_package(MPI COMPONENTS C Fortran REQUIRED)
if ( MPI_Fortran_HAVE_F08_MODULE )
    message( STATUS "MPI Library supports MPI_f08" )
else ( MPI_Fortran_HAVE_F08_MODULE )
    message( FATAL_ERROR "MPI Library does not support MPI_f08" )
endif ( MPI_Fortran_HAVE_F08_MODULE )

set(HDF5_PREFER_PARALLEL TRUE)
set(HDF5_USE_STATIC_LIBRARIES FALSE)
find_package( HDF5 COMPONENTS C Fortran REQUIRED )
if ( HDF5_IS_PARALLEL )
    message( STATUS "HDF5 Library supports parallel I/O" )
else ( HDF5_IS_PARALLEL )
    message( WARNING "HDF5 does not support parallel I/O" )
endif ( HDF5_IS_PARALLEL )

# CMake build options
option(REAL64 "Compile MGLET with double precision support" OFF)
if ( REAL64 )
    message( STATUS "Compiling MGLET in double precision mode" )
    add_compile_definitions(_MGLET_DOUBLE_PRECISION_=1)
else ( REAL64 )
    message( STATUS "Compiling MGLET in single precision mode" )
endif ( REAL64 )

# CMake build options
option(INT64 "Compile MGLET with 64 bit integer support (DOES NOT WORK)" OFF)
if ( INT64 )
    message( STATUS "Compiling MGLET with 64 bit integers" )
    message( WARNING "64-bit integer support does not work ATM!" )
    add_compile_definitions(_MGLET_INT64_=1)
else ( INT64 )
    message( STATUS "Compiling MGLET with 32 bit integers" )
endif ( INT64 )

# All fortran modules are placed here
set(CMAKE_Fortran_MODULE_DIRECTORY ${CMAKE_BINARY_DIR}/modules)

set(CMAKE_POSITION_INDEPENDENT_CODE ON)

# Flags for builds
if (CMAKE_Fortran_COMPILER_ID STREQUAL "GNU")
    # Common flags for all builds and languages
    set(MGLET_COMMON_FLAGS -Wall -Wextra -Wpedantic -Werror -Wno-maybe-uninitialized)

    # Flags for Release (=optimized builds), all languages
    set(MGLET_RELEASE_FLAGS -O3 -fopenmp-simd)

    # Debug flags for all languages
    set(MGLET_DEBUG_FLAGS -Og -g -D_MGLET_DEBUG_)

    # Common Fortran flags
    #
    # Regarding: -Wno-uninitialized
    # GFortran have several problems with the detection of uninitialized
    # variables, ref:
    # https://fortran-lang.discourse.group/t/gfortran-uninitialized-warnings/3838
    # https://gcc.gnu.org/bugzilla/show_bug.cgi?id=56670
    # https://gcc.gnu.org/bugzilla/show_bug.cgi?id=106089
    # https://gcc.gnu.org/bugzilla/show_bug.cgi?id=106108
    # therefore setting this flag to turn off this feature.
    #
    # Regarding: -ffree-line-length-none
    # The preprocessor macro __FILE__ is often used in MGLET. This can be a
    # long path. Up until Fortran 2018 there was a max line length of
    # 132 characters. A filename can easily overflow this length.
    # GCC allows to increase this length to infinity with this flag.
    # Fortran 2023 will have the maximum line length increased to 10000
    # (ten thousand) characters, so this is a problem that will disappear
    # as soon as a Fortran 2023 compiler is available.
    set(MGLET_COMMON_FORTRAN_FLAGS -std=f2018 -fimplicit-none -ffpe-trap=invalid,zero,overflow -Wno-unused-dummy-argument -Wno-uninitialized -ffree-line-length-none)

    # Fortran-specific debug flags
    set(MGLET_DEBUG_FORTRAN_FLAGS -fcheck=all)
elseif (CMAKE_Fortran_COMPILER_ID STREQUAL "IntelLLVM")
    # Flags for Release (=optimized builds), all languages
    set(MGLET_RELEASE_FLAGS -O3 -qopenmp-simd)

    # Debug flags for all languages
    set(MGLET_DEBUG_FLAGS -O1 -g -D_MGLET_DEBUG_)

    # Common Fortran flags
    set(MGLET_COMMON_FORTRAN_FLAGS -stand f18 -auto -warn all,noexternals,nounused,errors,stderrors)

    # Fortran-specific debug flags
    # TODO: add -init=snan,huge
    # https://community.intel.com/t5/Intel-Fortran-Compiler/ICE-with-locals-declared-in-a-block-construct-IFX/m-p/1489269#M166473
    # TODO: add -check shape
    set(MGLET_DEBUG_FORTRAN_FLAGS -traceback -ftrapuv -check arg_temp_created,bounds,contiguous,pointers)
endif()


# https://stackoverflow.com/questions/48895687/cmake-add-compile-options-error-with-release-debug-config
add_compile_options("${MGLET_COMMON_FLAGS}"
    "$<$<COMPILE_LANGUAGE:Fortran>:${MGLET_COMMON_FORTRAN_FLAGS}>"
    "$<$<CONFIG:Release,RelWithDebInfo>:${MGLET_RELEASE_FLAGS}>"
    "$<$<CONFIG:Debug>:${MGLET_DEBUG_FLAGS}>"
    "$<$<AND:$<COMPILE_LANGUAGE:Fortran>,$<CONFIG:Debug>>:${MGLET_DEBUG_FORTRAN_FLAGS}>"
)

if(POLICY CMP0126)
    cmake_policy(SET CMP0126 NEW)
endif()

# nlohman JSON parser library
include(FetchContent)
FetchContent_Declare(json
    GIT_REPOSITORY https://github.com/ArthurSonzogni/nlohmann_json_cmake_fetchcontent
    GIT_TAG v3.10.4
)

set(JSON_ImplicitConversions OFF)
FetchContent_GetProperties(json)
if(NOT json_POPULATED)
  FetchContent_Populate(json)
  set(JSON_BuildTests OFF CACHE INTERNAL "")
  add_subdirectory(${json_SOURCE_DIR} ${json_BINARY_DIR} EXCLUDE_FROM_ALL)
endif()

# Exprtk expression parser
FetchContent_Declare(exprtk
    GIT_REPOSITORY https://github.com/ArashPartow/exprtk
    GIT_TAG f46bffcd6966d38a09023fb37ba9335214c9b959
)
FetchContent_MakeAvailable(exprtk)

# Perform various compiler checks
add_subdirectory("cmake")

# Add code subdirectory
add_subdirectory("src")

# Add all extensions in the end (they depend on stuff compiled in "src")
MACRO(SUBDIRLIST result curdir)
  FILE(GLOB children RELATIVE ${curdir} ${curdir}/*)
  SET(dirlist "")
  FOREACH(child ${children})
    IF(IS_DIRECTORY ${curdir}/${child} AND NOT ${child} STREQUAL ".git")
      LIST(APPEND dirlist ${child})
    ENDIF()
  ENDFOREACH()
  SET(${result} ${dirlist})
ENDMACRO()
SUBDIRLIST(MGLET_EXTENSIONS ${CMAKE_CURRENT_SOURCE_DIR}/extensions)
FOREACH(ext ${MGLET_EXTENSIONS})
  ADD_SUBDIRECTORY(${CMAKE_CURRENT_SOURCE_DIR}/extensions/${ext})
ENDFOREACH()

add_subdirectory("tests")
