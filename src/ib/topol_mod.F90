MODULE topol_mod
    USE MPI_f08
    USE core_mod
    USE, INTRINSIC :: ISO_C_BINDING, ONLY: C_F_POINTER

    IMPLICIT NONE (type, external)
    PRIVATE

    TYPE :: topol_t
        INTEGER(intk) :: n = 0
        REAL(realk), POINTER, CONTIGUOUS :: topol(:)
        INTEGER(intk), POINTER, CONTIGUOUS :: bodyid(:)

        INTEGER(intk) :: nbody
        INTEGER(intk), ALLOCATABLE :: ids(:)
        CHARACTER(len=mglet_filename_max), ALLOCATABLE :: geometries(:)

        LOGICAL :: shmem_is_allocated = .FALSE.
        TYPE(real_shmem_arr), PRIVATE :: shmtopol
        TYPE(int_shmem_arr), PRIVATE :: shmbodyid
    CONTAINS
        PROCEDURE :: init
        PROCEDURE :: finish
    END TYPE topol_t

    PUBLIC :: topol_t, add_ids

CONTAINS

    SUBROUTINE init(this, blockconf, id)
        ! Reads the geometries specified in the blockconf. If the (optional)
        ! argument 'id' is present, then only geometries with that tag is read.
        !
        ! Geometries without an ID is assigned a unique ID automatically, in
        ! increasing order.

        ! Subroutine arguments
        CLASS(topol_t), INTENT(inout) :: this
        TYPE(config_t), INTENT(inout) :: blockconf
        INTEGER(intk), INTENT(in), OPTIONAL :: id

        ! Local variables
        TYPE(config_t) :: geometries, geometry
        INTEGER(intk) :: i, j, offset, ngeom, this_id, only_id

        REAL(realk), ALLOCATABLE :: vertices(:, :, :)
        REAL(realk), POINTER :: topol3d(:, :, :)
        INTEGER(intk), ALLOCATABLE :: ntri(:)
        LOGICAL, ALLOCATABLE :: is_binary(:)

        CHARACTER(len=64) :: jsonptr

        ! Check
        IF (this%shmem_is_allocated) CALL errr(__FILE__, __LINE__)

        ! Read which STL's to read
        this%nbody = 0
        CALL blockconf%get_size("/geometries", ngeom)

        ! Add ID's to all geometries
        CALL blockconf%get(geometries, "/geometries")
        CALL add_ids(geometries)

        only_id = 0
        IF (PRESENT(id)) THEN
            only_id = id
        END IF

        ! Count *number of* STL's to read
        DO i = 1, ngeom
            WRITE(jsonptr, '("/geometries/", I0)') i-1
            CALL blockconf%get(geometry, jsonptr)

            ! All geometries have either a user-specified id or have been
            ! assigned one
            CALL geometry%get_value("/id", this_id)

            IF (only_id == 0 .OR. this_id == only_id) THEN
                this%nbody = this%nbody + 1
            END IF

            ! Hmm. Check why this is neccesary. It shouln't be.. Maybe
            ! compiler bug? If this is not present, a memory leak appear
            CALL geometry%finish()
        END DO

        ! Read actual names to read
        ALLOCATE(this%geometries(this%nbody))
        ALLOCATE(this%ids(this%nbody))
        this%ids = 0

        j = 0
        DO i = 1, ngeom
            WRITE(jsonptr, '("/geometries/", I0)') i-1
            CALL blockconf%get(geometry, jsonptr)

            ! All geometries have either a user-specified id or have been
            ! assigned one
            CALL geometry%get_value("/id", this_id)

            IF (only_id == 0 .OR. this_id == only_id) THEN
                j = j + 1
                this%geometries(j) = &
                    REPEAT(" ", LEN(this%geometries(j)))
                CALL geometry%get_value("/file", this%geometries(j))
                this%ids(j) = this_id
            END IF

            ! Hmm. Check why this is neccesary. It shouln't be.. Maybe
            ! compiler bug? If this is not present, a memory leak appear
            CALL geometry%finish()
        END DO

        ! Read STL's
        ALLOCATE(ntri(this%nbody))
        ntri = 0
        ALLOCATE(is_binary(this%nbody))
        is_binary = .FALSE.

        IF (myid == 0) THEN
            DO i = 1, this%nbody
                CALL stl_info(this%geometries(i), is_binary(i), ntri(i))
            END DO
        END IF

        CALL MPI_Bcast(ntri, INT(this%nbody, int32), mglet_mpi_int, &
            0, MPI_COMM_WORLD)

        CALL MPI_Bcast(is_binary, INT(this%nbody, int32), MPI_LOGICAL, &
            0, MPI_COMM_WORLD)

        ! Overall number of triangles to be read in
        ! Ranks with smhid > 0 shall not allocate storage space, so keep
        ! this zero at these as of now
        !
        ! Also, C_F_POINTER does not like to associate zero size allocations,
        ! so always let each SHM master have at least space for 1 triangle
        ! allocated, otherwise "topol" is not correctly associated in
        ! blockbp/blockbpcc.
        this%n = 0
        IF (shmid == 0) this%n = MAX(SUM(ntri), 1)

        ! Allocate storage - if shmid > 0 then ntopol == 0 and no space is
        ! allocated
        !
        ! Wrapper shmem_alloc want the number of elements to allocate,
        ! not bytes.
        CALL shmem_alloc(this%shmtopol, &
            3*3*INT(this%n, int64))

        ! Only one rank per compute node/SHM group do this to spare memory
        IF (shmid == 0) THEN
            CALL C_F_POINTER(this%shmtopol%baseptr, &
                topol3d, [3, 3, this%n])
            offset = 1
            DO i = 1, this%nbody
                ! If the topol have zero triangles, do nothing
                IF (ntri(i) <= 0) CYCLE

                ALLOCATE(vertices(3, 3, ntri(i)))

                ! Global rank 0 read the STL
                IF (myid == 0) THEN
                    CALL stl_read(this%geometries(i), &
                        is_binary(i), vertices)
                END IF

                ! Rank 0 broadcast the vertices to other shm master ranks
                CALL MPI_Bcast(vertices, 3*3*ntri(i), &
                    mglet_mpi_real, 0, shm_masters_comm)

                ! Now all ranks with shmid == 0 have the vertices data
                DO j = 1, ntri(i)
                    IF (offset > this%n) THEN
                        CALL errr(__FILE__, __LINE__)
                    END IF

                    topol3d(1, :, offset) = vertices(:, 1, j) ! x1, y1, z1
                    topol3d(2, :, offset) = vertices(:, 2, j) ! x2, y2, z2
                    topol3d(3, :, offset) = vertices(:, 3, j) ! x3, y3, z3

                    offset = offset + 1
                END DO
                DEALLOCATE(vertices)
            END DO
            NULLIFY(topol3d)
        END IF
        CALL this%shmtopol%barrier()

        ! Allocate and set bodyid
        CALL shmem_alloc(this%shmbodyid, INT(this%n, int64))
        IF (shmid == 0) THEN
            offset = 1
            DO i = 1, this%nbody
                IF (ntri(i) <= 0) CYCLE
                this%shmbodyid%arr(offset:offset+ntri(i)-1) = &
                    this%ids(i)
                offset = offset + ntri(i)
            END DO
        END IF
        CALL this%shmbodyid%barrier()

        ! In the end ntopol is set correct at all processes
        ! Also in the case of zero trinagles, ntopol is zero from now on..
        this%n = SUM(ntri)

        DEALLOCATE(is_binary, ntri)

        ! Set convenience pointers
        CALL C_F_POINTER(this%shmtopol%baseptr, this%topol, &
            [SIZE(this%shmtopol%arr)])
        CALL C_F_POINTER(this%shmbodyid%baseptr, this%bodyid, &
            [SIZE(this%shmbodyid%arr)])

        this%shmem_is_allocated = .TRUE.
    END SUBROUTINE init


    IMPURE ELEMENTAL SUBROUTINE finish(this)
        CLASS(topol_t), INTENT(inout) :: this

        this%nbody = 0
        this%n = 0
        IF (ALLOCATED(this%ids)) DEALLOCATE(this%ids)
        IF (ALLOCATED(this%geometries)) DEALLOCATE(this%geometries)

        NULLIFY(this%topol)
        NULLIFY(this%bodyid)

        IF (this%shmem_is_allocated) THEN
            CALL shmem_dealloc(this%shmtopol)
            CALL shmem_dealloc(this%shmbodyid)
        END IF
        this%shmem_is_allocated = .FALSE.
    END SUBROUTINE finish


    SUBROUTINE add_ids(geometries)
        ! Subroutine that adds ID tag to geometries that are not already
        ! equipped with this

        ! Subroutine arguments
        TYPE(config_t) :: geometries

        ! Local variables
        INTEGER(intk) :: ngeom, last_id, this_id, i
        INTEGER(intk), ALLOCATABLE :: ids(:)
        TYPE(config_t) :: geometry
        LOGICAL :: has_id
        CHARACTER(len=64) :: jsonptr

        CALL geometries%get_size("", ngeom)
        ALLOCATE(ids(ngeom))
        ids = 0

        ! This algorithm is a two-stage pass. First all geomeetries *with* id
        ! are checked and assigned. Thereafter geometries without ID are
        ! assigned free ID's in the range above the predefined range
        DO i = 1, ngeom
            WRITE(jsonptr, '("/", I0)') i-1
            CALL geometries%get(geometry, jsonptr)

            has_id = geometry%exists("/id")
            IF (has_id) THEN
                CALL geometry%get_value("/id", this_id)
                IF (this_id < 1) THEN
                    WRITE(*, *) "Invalid id:", this_id
                    CALL errr(__FILE__, __LINE__)
                END IF
                ids(i) = this_id
            END IF

            ! If this is not present, a memory leak appear
            CALL geometry%finish()
        END DO

        ! Now the ids-array are willed with predetermined ID's, assign
        ! is to the rest
        last_id = MAXVAL(ids)
        DO i = 1, ngeom
            WRITE(jsonptr, '("/", I0)') i-1
            CALL geometries%get(geometry, jsonptr)

            IF (ids(i) < 1) THEN
                last_id = last_id + 1
                CALL geometry%set_value("/id", last_id)
            END IF

            ! If this is not present, a memory leak appear
            CALL geometry%finish()
        END DO
    END SUBROUTINE add_ids
END MODULE topol_mod
