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
        TYPE(field_t) :: bp_f

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


    SUBROUTINE set_auavaw(this, bu_p, bv_p, bw_p, au_p, av_p, aw_p)
        ! Subroutine arguments
        CLASS(gc_stencils_t), INTENT(inout) :: this
        REAL(realk), TARGET, CONTIGUOUS, INTENT(in) :: bu_p(:)
        REAL(realk), TARGET, CONTIGUOUS, INTENT(in) :: bv_p(:)
        REAL(realk), TARGET, CONTIGUOUS, INTENT(in) :: bw_p(:)
        REAL(realk), TARGET, CONTIGUOUS, INTENT(in) :: au_p(:)
        REAL(realk), TARGET, CONTIGUOUS, INTENT(in) :: av_p(:)
        REAL(realk), TARGET, CONTIGUOUS, INTENT(in) :: aw_p(:)

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
            bu(1:kk, 1:jj, 1:ii) => bu_p(ip3:ip3+kk*jj*ii-1)
            au(1:kk, 1:jj, 1:ii) => au_p(ip3:ip3+kk*jj*ii-1)
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
            bu(1:kk, 1:jj, 1:ii) => bv_p(ip3:ip3+kk*jj*ii-1)
            au(1:kk, 1:jj, 1:ii) => av_p(ip3:ip3+kk*jj*ii-1)
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

            bu(1:kk, 1:jj, 1:ii) => bw_p(ip3:ip3+kk*jj*ii-1)
            au(1:kk, 1:jj, 1:ii) => aw_p(ip3:ip3+kk*jj*ii-1)
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


    SUBROUTINE get_auavaw(this, bu, bv, bw, au_p, av_p, aw_p)
        ! Subroutine arguments
        CLASS(gc_stencils_t), INTENT(inout) :: this
        REAL(realk), INTENT(in), CONTIGUOUS :: bu(:)
        REAL(realk), INTENT(in), CONTIGUOUS :: bv(:)
        REAL(realk), INTENT(in), CONTIGUOUS :: bw(:)
        REAL(realk), TARGET, INTENT(out), CONTIGUOUS :: au_p(:)
        REAL(realk), TARGET, INTENT(out), CONTIGUOUS :: av_p(:)
        REAL(realk), TARGET, INTENT(out), CONTIGUOUS :: aw_p(:)

        ! Local variables
        INTEGER(intk) :: imygrid, idx, igrid, ind, ip3
        INTEGER(intk) :: aucells
        INTEGER(intk) :: kk, jj, ii, k, j, i
        REAL(realk), POINTER, CONTIGUOUS :: au(:, :, :)

        ! First set AU to value of BU
        au_p = bu
        av_p = bv
        aw_p = bw

        ! Now set AU, AV, AW
        all_grids: DO imygrid = 1, nmygrids
            igrid = mygrids(imygrid)
            CALL get_mgdims(kk, jj, ii, igrid)
            CALL get_ip3(ip3, igrid)

            au(1:kk, 1:jj, 1:ii) => au_p(ip3:ip3+kk*jj*ii-1)
            aucells = SIZE(this%auind(imygrid)%arr)
            DO idx = 1, aucells
                ind = this%auind(imygrid)%arr(idx)
                CALL ind2sub(ind, k, j, i, kk, jj, ii)
                au(k, j, i) = this%auvalue(imygrid)%arr(idx)
            END DO

            au(1:kk, 1:jj, 1:ii) => av_p(ip3:ip3+kk*jj*ii-1)
            aucells = SIZE(this%avind(imygrid)%arr)
            DO idx = 1, aucells
                ind = this%avind(imygrid)%arr(idx)
                CALL ind2sub(ind, k, j, i, kk, jj, ii)
                au(k, j, i) = this%avvalue(imygrid)%arr(idx)
            END DO

            au(1:kk, 1:jj, 1:ii) => aw_p(ip3:ip3+kk*jj*ii-1)
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
