# MGLET Offloading

## Build instructions

`CMakePresets.json` provides build configurations for tested and supported target architectures and compilers.

To create a new preset, make sure to define the following CMake variables:
- `MGLET_OFFLOAD` (`ON` or `OFF`): Toggles OpenMP support and enables the `_MGLET_OFFLOAD_` compile definition to use offload-specific implementations.
- `MGLET_OFFLOAD_COMPILE_FLAGS`: Compiler-specific flags for OpenMP offloading
- `MGLET_OFFLOAD_ARCH_FLAGS`: Compiler-specific flag specifying the offload-device's architecture

> [!WARNING]  
> Do NOT toggle `MGLET_OFFLOAD` without specifying an offload architecture. CPU-only OpenMP is not supported in MGLET and may cause bad performance if ill-configured.

## Implementation Guidelines

### Grid-Cell Loops

Grid-cell loops can be offloaded by assigning one grid to each execution unit (team), while performing grid-internal work on threads within. Make sure that variables that differ per grid are marked as `private`.

```Fortran
!$omp target teams distribute private(kk, jj, ii, u, v, w)
DO igrid = 1, nmygrids
    CALL get_mgdims(kk, jj, ii, igrid)

    CALL get_grid3_real(u, u_f, igrid)
    CALL get_grid3_real(v, v_f, igrid)
    CALL get_grid3_real(w, w_f, igrid)

    ! ...

    !$omp parallel
    CALL do_work_grid(kk, jj, ii, u, v, w, ...)
    !$omp end parallel
END DO
!$omp end target teams distribute

SUBROUTINE do_work_grid(kk, jj, ii, u, v, w)
    ! ...

    !$omp do collapse(3)
    DO i = 3, ii-2
        DO j = 3, jj-2
            DO k = 3, kk-2
                u(k, j, i) = ...
            END DO
        END DO
    END DO
    !$omp end do
END SUBROUTINE do_work_grid
```

### Custom field mappers

Custom mappers for the derived field types are provided in `fieldmapper_mod`. Using this module enables the usage of the provided custom mappers. Note that transitive uses are not supported. Rule of thumb: if a custom mapper is needed, use `fieldmapper_mod` in the tightest scope possible.

> [!WARNING]  
> The default mapper of derived field types are overridden in `fieldmapper_mod`. When entering a GPU kernel, make sure that the use of `fieldmapper_mod` is not visible in its scope. An overridden default mapper causes significant kernel launch overhead due to redundant implicit copies of descriptors.

### Common Blocks

Non-scalar variables in common-blocks declared as target with `declare target` must be mapped using the `always` map-type modifier to enforce performing an actual data transfer (see [AMD HPC Training Examples](https://github.com/amd/HPCTrainingExamples/tree/main/Pragma_Examples/OpenMP/Fortran/Common_blocks_on_device)).

```Fortran
MODULE common_blocks_mod
    IMPLICIT NONE (type, external)

    REAL, ALLOCATABLE :: common_block_variable(:)
    !$omp declare target(common_block_variable)
CONTAINS
    SUBROUTINE init()
        ALLOCATE(common_block_variable(n))
        !$omp target enter data map(always, to: common_block_variable)
    END SUBROUTINE init
END MODULE common_blocks_mod
```
