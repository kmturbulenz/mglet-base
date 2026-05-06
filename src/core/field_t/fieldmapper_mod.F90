MODULE fieldmapper_mod
    USE realfield_mod, ONLY: field_t

    IMPLICIT NONE(type, external)
    ! Not PRIVATE'ized

    PRIVATE :: field_t

    ! Use this module directly. Mapper declarations are not visible through
    ! transitive use statements if any intermediate module is PRIVATE'ized.
    ! Mangle custom mapper-ids to avoid conflicts with variable names.
    !$omp declare mapper(field_t :: t) map(t%arr, t%buffers, t%ptr, t%length)
    !$omp declare mapper(field_t__map_arr: field_t :: t) map(t%arr)
    !$omp declare mapper(field_t__map_buffers: field_t :: t) map(t%buffers)
    !$omp declare mapper(field_t__map_deferred: field_t :: t) &
    !$omp& map(alloc: t%arr, t%buffers) map(to: t%ptr, t%length)
END MODULE fieldmapper_mod
