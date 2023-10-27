MODULE gc_stencils_mod
    USE HDF5
    USE core_mod, ONLY: realk, intk, eps, errr, int_stencils_t, &
        real_stencils_t, nmygrids, get_imygrid, sub2ind, field_t, &
        mygrids, nmygrids, get_mgdims, hdf5common_attr_read, &
        hdf5common_open, hdf5common_close, hdf5common_attr_write, get_ip3, &
        idim3d, myid, ind2sub
    USE stencils_mod, ONLY: stencils_t
    IMPLICIT NONE(type, external)
    PRIVATE

    TYPE, EXTENDS(stencils_t) :: gc_stencils_t
    CONTAINS
        PROCEDURE :: set_bp
        PROCEDURE :: set_auavaw
        PROCEDURE :: get_auavaw
        PROCEDURE :: read
        PROCEDURE :: write
        FINAL :: destructor
    END TYPE gc_stencils_t

    PUBLIC :: gc_stencils_t
CONTAINS

    IMPURE ELEMENTAL SUBROUTINE destructor(this)
        TYPE(gc_stencils_t), INTENT(inout) :: this
        CALL this%finish()
    END SUBROUTINE destructor


    SUBROUTINE set_bp(this, bp_f)
        ! Subroutine arguments
        CLASS(gc_stencils_t), INTENT(inout) :: this
        TYPE(field_t), INTENT(in) :: bp_f

        ! Local variables
        INTEGER(intk) :: imygrid, igrid, ind, bpcells
        INTEGER(intk) :: kk, jj, ii, k, j, i
        REAL(realk), POINTER, CONTIGUOUS :: bp(:, :, :)

        IF (.NOT. ALLOCATED(this%bpind)) ALLOCATE(this%bpind(nmygrids))

        all_grids: DO imygrid = 1, nmygrids
            igrid = mygrids(imygrid)
            CALL get_mgdims(kk, jj, ii, igrid)
            CALL bp_f%get_ptr(bp, igrid)

            ! Count blocked cells
            bpcells = 0
            DO i = 1, ii
                DO j = 1, jj
                    DO k = 1, kk
                        IF (bp(k, j, i) < 0.5) bpcells = bpcells + 1
                    END DO
                END DO
            END DO

            ! Set blocked cells
            ALLOCATE(this%bpind(imygrid)%arr(bpcells))
            bpcells = 0
            DO i = 1, ii
                DO j = 1, jj
                    DO k = 1, kk
                        IF (bp(k, j, i) < 0.5) THEN
                            CALL sub2ind(ind, k, j, i, kk, jj, ii)
                            bpcells = bpcells + 1
                            this%bpind(imygrid)%arr(bpcells) = ind
                        END IF
                    END DO
                END DO
            END DO

        END DO all_grids
    END SUBROUTINE set_bp


    SUBROUTINE set_auavaw(this, bu_f, bv_f, bw_f, au_f, av_f, aw_f)
        ! Subroutine arguments
        CLASS(gc_stencils_t), INTENT(inout) :: this
        TYPE(field_t), INTENT(in) :: bu_f
        TYPE(field_t), INTENT(in) :: bv_f
        TYPE(field_t), INTENT(in) :: bw_f
        TYPE(field_t), INTENT(in) :: au_f
        TYPE(field_t), INTENT(in) :: av_f
        TYPE(field_t), INTENT(in) :: aw_f

        ! Local variables
        INTEGER(intk) :: imygrid, igrid, ind, ip3
        INTEGER(intk) :: aucells
        INTEGER(intk) :: kk, jj, ii, k, j, i
        REAL(realk), POINTER, CONTIGUOUS :: bu(:, :, :), au(:, :, :)

        IF (.NOT. ALLOCATED(this%auind)) ALLOCATE(this%auind(nmygrids))
        IF (.NOT. ALLOCATED(this%avind)) ALLOCATE(this%avind(nmygrids))
        IF (.NOT. ALLOCATED(this%awind)) ALLOCATE(this%awind(nmygrids))

        IF (.NOT. ALLOCATED(this%auvalue)) ALLOCATE(this%auvalue(nmygrids))
        IF (.NOT. ALLOCATED(this%avvalue)) ALLOCATE(this%avvalue(nmygrids))
        IF (.NOT. ALLOCATED(this%awvalue)) ALLOCATE(this%awvalue(nmygrids))

        all_grids: DO imygrid = 1, nmygrids
            igrid = mygrids(imygrid)
            CALL get_mgdims(kk, jj, ii, igrid)
            CALL get_ip3(ip3, igrid)

            ! Count blocked cells aucells
            CALL bu_f%get_ptr(bu, igrid)
            CALL au_f%get_ptr(au, igrid)
            aucells = 0
            DO i = 1, ii
                DO j = 1, jj
                    DO k = 1, kk
                        IF (ABS(bu(k, j, i) - au(k, j, i)) >= eps) THEN
                            aucells = aucells + 1
                        END IF
                    END DO
                END DO
            END DO
            ALLOCATE(this%auind(imygrid)%arr(aucells))
            ALLOCATE(this%auvalue(imygrid)%arr(aucells))

            ! Set aucells
            aucells = 0
            DO i = 1, ii
                DO j = 1, jj
                    DO k = 1, kk
                        IF (ABS(bu(k, j, i) - au(k, j, i)) >= eps) THEN
                            CALL sub2ind(ind, k, j, i, kk, jj, ii)
                            aucells = aucells + 1
                            this%auind(imygrid)%arr(aucells) = ind
                            this%auvalue(imygrid)%arr(aucells) = &
                                au(k, j, i)
                        END IF
                    END DO
                END DO
            END DO

            ! Notice re-using pointer
            CALL bv_f%get_ptr(bu, igrid)
            CALL av_f%get_ptr(au, igrid)
            aucells = 0
            DO i = 1, ii
                DO j = 1, jj
                    DO k = 1, kk
                        IF (ABS(bu(k, j, i) - au(k, j, i)) >= eps) THEN
                            aucells = aucells + 1
                        END IF
                    END DO
                END DO
            END DO
            ALLOCATE(this%avind(imygrid)%arr(aucells))
            ALLOCATE(this%avvalue(imygrid)%arr(aucells))

            aucells = 0
            DO i = 1, ii
                DO j = 1, jj
                    DO k = 1, kk
                        IF (ABS(bu(k, j, i) - au(k, j, i)) >= eps) THEN
                            CALL sub2ind(ind, k, j, i, kk, jj, ii)
                            aucells = aucells + 1
                            this%avind(imygrid)%arr(aucells) = ind
                            this%avvalue(imygrid)%arr(aucells) = &
                                au(k, j, i)
                        END IF
                    END DO
                END DO
            END DO

            CALL bw_f%get_ptr(bu, igrid)
            CALL aw_f%get_ptr(au, igrid)
            aucells = 0
            DO i = 1, ii
                DO j = 1, jj
                    DO k = 1, kk
                        IF (ABS(bu(k, j, i) - au(k, j, i)) >= eps) THEN
                            aucells = aucells + 1
                        END IF
                    END DO
                END DO
            END DO
            ALLOCATE(this%awind(imygrid)%arr(aucells))
            ALLOCATE(this%awvalue(imygrid)%arr(aucells))

            aucells = 0
            DO i = 1, ii
                DO j = 1, jj
                    DO k = 1, kk
                        IF (ABS(bu(k, j, i) - au(k, j, i)) >= eps) THEN
                            CALL sub2ind(ind, k, j, i, kk, jj, ii)
                            aucells = aucells + 1
                            this%awind(imygrid)%arr(aucells) = ind
                            this%awvalue(imygrid)%arr(aucells) = &
                                au(k, j, i)
                        END IF
                    END DO
                END DO
            END DO

        END DO all_grids
    END SUBROUTINE set_auavaw


    SUBROUTINE get_auavaw(this, bu, bv, bw, au_f, av_f, aw_f)
        ! Subroutine arguments
        CLASS(gc_stencils_t), INTENT(inout) :: this
        TYPE(field_t), INTENT(in) :: bu
        TYPE(field_t), INTENT(in) :: bv
        TYPE(field_t), INTENT(in) :: bw
        TYPE(field_t), INTENT(inout) :: au_f
        TYPE(field_t), INTENT(inout) :: av_f
        TYPE(field_t), INTENT(inout) :: aw_f

        ! Local variables
        INTEGER(intk) :: imygrid, idx, igrid, ind, ip3
        INTEGER(intk) :: aucells
        INTEGER(intk) :: kk, jj, ii, k, j, i
        REAL(realk), POINTER, CONTIGUOUS :: au(:, :, :)

        ! First set AU to value of BU
        au_f%arr = bu%arr
        av_f%arr = bv%arr
        aw_f%arr = bw%arr

        ! Now set AU, AV, AW
        all_grids: DO imygrid = 1, nmygrids
            igrid = mygrids(imygrid)
            CALL get_mgdims(kk, jj, ii, igrid)
            CALL get_ip3(ip3, igrid)

            CALL au_f%get_ptr(au, igrid)
            aucells = SIZE(this%auind(imygrid)%arr)
            DO idx = 1, aucells
                ind = this%auind(imygrid)%arr(idx)
                CALL ind2sub(ind, k, j, i, kk, jj, ii)
                au(k, j, i) = this%auvalue(imygrid)%arr(idx)
            END DO

            CALL av_f%get_ptr(au, igrid)
            aucells = SIZE(this%avind(imygrid)%arr)
            DO idx = 1, aucells
                ind = this%avind(imygrid)%arr(idx)
                CALL ind2sub(ind, k, j, i, kk, jj, ii)
                au(k, j, i) = this%avvalue(imygrid)%arr(idx)
            END DO

            CALL aw_f%get_ptr(au, igrid)
            aucells = SIZE(this%awind(imygrid)%arr)
            DO idx = 1, aucells
                ind = this%awind(imygrid)%arr(idx)
                CALL ind2sub(ind, k, j, i, kk, jj, ii)
                au(k, j, i) = this%awvalue(imygrid)%arr(idx)
            END DO

        END DO all_grids
    END SUBROUTINE get_auavaw


    SUBROUTINE read(this)
        ! Subroutine arguments
        CLASS(gc_stencils_t), INTENT(inout) :: this

        ! Local variables
        INTEGER(HID_T) :: file_id
        CHARACTER(len=16) :: ibtype

        CALL hdf5common_open(this%file, 'r', file_id)

        ibtype = REPEAT(" ", LEN(ibtype))
        CALL hdf5common_attr_read("IBTYPE", ibtype, file_id)
        IF (myid == 0) THEN
            IF (TRIM(ibtype) /= 'GHOSTCELL') THEN
                WRITE(*, *) "Wrong ibtype, expected 'GHOSTCELL', got: ", &
                    TRIM(ibtype)
                CALL errr(__FILE__, __LINE__)
            END IF
        END IF
        CALL this%read_stencils(file_id)
        CALL hdf5common_close(file_id)
    END SUBROUTINE read


    SUBROUTINE write(this)
        ! Subroutine arguments
        CLASS(gc_stencils_t), INTENT(inout) :: this

        ! Local variables
        INTEGER(HID_T) :: file_id

        CALL hdf5common_open(this%file, 'w', file_id)
        CALL this%write_stencils(file_id)
        CALL hdf5common_attr_write("IBTYPE", "GHOSTCELL", file_id)
        CALL hdf5common_close(file_id)
    END SUBROUTINE write
END MODULE gc_stencils_mod
