MODULE fieldmapper_mod
    USE field_mod, ONLY: field_t, intfield_t

    IMPLICIT NONE(type, external)
    ! Not PRIVATE'ized

    PRIVATE :: field_t, intfield_t

    ! Use this module directly. Mapper declarations are not visible through
    ! transitive use statements if any intermediate module is PRIVATE'ized.
    ! Mangle custom mapper-ids to avoid conflicts with variable names.
    !$omp declare mapper(field_t :: t) map(t%arr, t%buffers, t%ptr, t%length)
    !$omp declare mapper(field_t__map_arr: field_t :: t) map(t%arr)
    !$omp declare mapper(field_t__map_buffers: field_t :: t) map(t%buffers)

    !$omp declare mapper(intfield_t :: t) map(t%arr, t%ptr, t%length)
    !$omp declare mapper(intfield_t__map_arr: intfield_t :: t) map(t%arr)
END MODULE fieldmapper_mod
