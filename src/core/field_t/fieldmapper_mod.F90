MODULE fieldmapper_mod
    USE realfield_mod, ONLY: field_t
    IMPLICIT NONE
    ! NOT private by default

    PRIVATE :: field_t

    ! Why override the default mapper?
    ! (Intel) Only maps the underlying allocations of field_t if the default
    !         mapper is overridden and the allocatable members are listed
    ! (LLVM)  Maps the underlying allocations of field_t by default


    ! (Intel) Why is t%name mapped: Intel currently does not support updating
    !         from and to the host with custom mappers, the compiler segfaults
    !         during compilation. Thus, all "target update" directives must use
    !         field%arr directly. However, still it seems the actual data being
    !         mapped is not ONLY arr but also other members. Thus, t%name gets
    !         overridden when mapped back. This is not a sustainable fix, but
    !         just something to get it working on Intel.
    !         To workaround for Intel, override the default mapper:
    ! !$omp declare mapper(field_t :: t) &
    ! !$omp& map(t%arr, t%buffers, t%ptr, t%length, t%name)
    !         and add t%name to mapfielddefer.

    !$omp declare mapper(field_t :: t) &
    !$omp& map(t%arr, t%buffers, t%ptr, t%length)
    !$omp declare mapper(maparr: field_t :: t) map(t%arr)
    !$omp declare mapper(mapbuffers: field_t :: t) map(t%buffers)
    !$omp declare mapper(mapfielddefer: field_t :: t) &
    !$omp& map(alloc: t%arr, t%buffers) map(to: t%ptr, t%length)
END MODULE fieldmapper_mod
