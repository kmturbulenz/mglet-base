# MGLET Offloading

## Table of Contents

- [Build Instructions](#build-instructions)
   - [Supported Compilers](#supported-compilers)
   - [Creating CMake Presets](#creating-cmake-presets-for-offloading)
   - [Enabling Profiling Annotations](#enabling-profiling-annotations)
- [Run Instructions](#run-instructions)
- [Development Guidelines](#development-guidelines)
   - [Custom field mappers](#custom-field-mappers)
   - [Common Blocks](#common-blocks)
   - [Field Subroutines](#field-subroutines)
   - [Connect Communication](#connect-communication)
   - [using Profiling Annotations](#using-profiling-annotations)

# Build Instructions

## Supported Compilers

`CMakePresets.json` provides build configurations for tested and supported target architectures and compilers.

Currently, only the clang-based [AMD compilers](https://github.com/amd/InfinityHub-CI/blob/main/fortran/README.md) support MGLET offloading. While a wide range of architectures are available, MGLET is tested on `gfx942` and `gfx1201`. It can not be guaranteed that other architectures run without problems.

## Creating CMake Presets

To create a new preset, make sure to define the following CMake variables:
| Variable | Description |
| --- | --- |
| `MGLET_OFFLOAD` | Toggles (`ON` or `OFF`) OpenMP support, which enables the use of offload-specific implementations in MGLET. |
| `MGLET_OFFLOAD_COMPILE_FLAGS` | Compiler-specific flags for OpenMP offloading (e.g. `-fopenmp-version=...` for clang-based compilers). |
| `MGLET_OFFLOAD_ARCH_FLAGS` | Compiler-specific flag to specify the offload-device's architecture (e.g. `--offload-arch=...` for clang-based compilers). |

> [!CAUTION]  
> MGLET does not support CPU-parallelization with OpenMP. Performance degradation of hybrid MPI+OpenMP shall be avoided and OpenMP directives are layed out to benefit running in target regions. `OMP_NUM_THREADS=1` might work, but is generally not tested or supported in MGLET.

## Enabling Profiling Annotations

MGLET supports profiling annotations, such as markers and ranges from `roctx` and `nvtx`. To build with annotation support enabled, set `MGLET_PROFILE_ANNOTATIONS={roctx|nvtx}`. Linking against the underlying libraries requires that either `rocm` or `cuda` respectively is loaded in the build environment.

When enabled, MGLET timers are annotated as profiling ranges by default for both CPU-only and offloaded runs. For offload builds specifically, data copies using the [relevant helper functions](#custom-field-mappers) are also annotated.

# Run Instructions

## Simple Profiling

To collect metrics and profiling annotations during runtime, different tools can be used.

### rocm-based

### cuda-based

`nsys`: collect traces 

# Development Guidelines

## Custom field mappers

Custom mappers for derived field types are provided in `fieldmapper_mod`. Using this module enables the usage of custom mappers to map fields between the host and device. Note that transitive uses of this module are not supported.

> [!CAUTION]  
> The default mapper of derived field types are overridden in `fieldmapper_mod`. When entering a GPU kernel, make sure that `fieldmapper_mod` is not visible in its scope. An overridden default mapper causes significant kernel launch overhead due to implicit copies of descriptors.

To copy cell or buffer data between the host and device, use the helper functions provided in `field_helper_mod`.

Rule of thumb: Never USE `fieldmapper_mod` unless absolutely necessary. Limit its usage to the innermost core implementations.

## Common Blocks

Non-scalar variables in common-blocks marked with `declare target` must be entered to the device using the `always` map-type modifier to enforce performing an actual data transfer (see [AMD HPC Training Examples](https://github.com/amd/HPCTrainingExamples/tree/main/Pragma_Examples/OpenMP/Fortran/Common_blocks_on_device)).

Examples of this paradigm can for example be found in `corefields_mod`, `grids_mod` and `pointers_mod`.

## Field Subroutines

Subroutines, such as `get_ptr` and `get_buffer`, must not be called on the device due to compilers having trouble with [polymorphics](https://fortran-lang.org/learn/oop_features_in_fortran/object_oriented_programming_techniques/#polymorphism) (use of `CLASS`). Instead, use the helper subroutines `get_grid3_real` or `get_grid1_real` etc. provided in `field_helper_mod`.

## Connect Communication

To make use of device to device communication in connect, use the `conn` subroutine. Depending on whether offloading is enabled in CMake, `conn` will use device or host buffers and communication.

## Using Profiling Annotations

To annotate the code with profiling ranges or markers, use `profile_tools_mod`.

> [!IMPORTANT]
> Profiling annotations are not meant to be used in production. They incur a small overhead and shall be used in development or for debugging. Thus, make sure not to leave calls to profiling tools that are not guarded with preprocessor macros when merging to the master branch.
