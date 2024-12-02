MODULE gc_stencils_mod
    USE HDF5
    USE core_mod, ONLY: realk, intk, eps, errr, nmygrids, sub2ind, field_t, &
        mygrids, nmygrids, get_mgdims, hdf5common_attr_read, hdf5common_open, &
        hdf5common_close, hdf5common_attr_write, get_ip3, myid, ind2sub, &
        get_fieldptr, get_field
    USE stencils_mod, ONLY: stencils_t
    IMPLICIT NONE(type, external)
    PRIVATE

    TYPE, EXTENDS(stencils_t) :: gc_stencils_t
    CONTAINS
        PROCEDURE :: set_bp
        PROCEDURE :: set_auavaw
        PROCEDURE :: get_auavaw
        PROCEDURE :: calc_volp
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


    SUBROUTINE get_auavaw(this, bu_f, bv_f, bw_f, areau_f, areav_f, areaw_f)
        ! Subroutine arguments
        CLASS(gc_stencils_t), INTENT(inout) :: this
        TYPE(field_t), INTENT(in) :: bu_f
        TYPE(field_t), INTENT(in) :: bv_f
        TYPE(field_t), INTENT(in) :: bw_f
        TYPE(field_t), INTENT(inout) :: areau_f
        TYPE(field_t), INTENT(inout) :: areav_f
        TYPE(field_t), INTENT(inout) :: areaw_f

        ! Local variables
        TYPE(field_t), POINTER :: ddx_f, ddy_f, ddz_f
        INTEGER(intk) :: imygrid, idx, igrid, ind
        INTEGER(intk) :: aucells
        INTEGER(intk) :: kk, jj, ii, k, j, i
        REAL(realk) :: area
        REAL(realk), POINTER, CONTIGUOUS :: areau(:, :, :)
        REAL(realk), POINTER, CONTIGUOUS :: bu(:, :, :)
        REAL(realk), POINTER, CONTIGUOUS :: ddx(:), ddy(:), ddz(:)

        CALL get_field(ddx_f, "DDX")
        CALL get_field(ddy_f, "DDY")
        CALL get_field(ddz_f, "DDZ")

        ! Now set AU, AV, AW ( absolute areas in unit [L^2] )
        all_grids: DO imygrid = 1, nmygrids

            igrid = mygrids(imygrid)
            CALL get_mgdims(kk, jj, ii, igrid)

            ! Getting the pressure cell side lengths
            CALL ddx_f%get_ptr(ddx, igrid)
            CALL ddy_f%get_ptr(ddy, igrid)
            CALL ddz_f%get_ptr(ddz, igrid)

            ! Area AU from BU
            CALL areau_f%get_ptr(areau, igrid)
            CALL bu_f%get_ptr(bu, igrid)
            DO i = 1, ii
                DO j = 1, jj
                    DO k = 1, kk
                        area = ddz(k) * ddy(j)
                        areau(k, j, i) = bu(k, j, i) * area
                    END DO
                END DO
            END DO
            ! Refine for intersected cells
            aucells = SIZE(this%auind(imygrid)%arr)
            DO idx = 1, aucells
                ind = this%auind(imygrid)%arr(idx)
                CALL ind2sub(ind, k, j, i, kk, jj, ii)
                area = ddz(k) * ddy(j)
                areau(k, j, i) = this%auvalue(imygrid)%arr(idx) * area
            END DO

            ! Area AV from BV
            CALL areav_f%get_ptr(areau, igrid)
            CALL bv_f%get_ptr(bu, igrid)
            DO i = 1, ii
                DO j = 1, jj
                    DO k = 1, kk
                        area = ddx(i) * ddz(k)
                        areau(k, j, i) = bu(k, j, i) * area
                    END DO
                END DO
            END DO
            ! Refine for intersected cells
            aucells = SIZE(this%avind(imygrid)%arr)
            DO idx = 1, aucells
                ind = this%avind(imygrid)%arr(idx)
                CALL ind2sub(ind, k, j, i, kk, jj, ii)
                area = ddx(i) * ddz(k)
                areau(k, j, i) = this%avvalue(imygrid)%arr(idx) * area
            END DO

            ! Area AW from BW
            CALL areaw_f%get_ptr(areau, igrid)
            CALL bw_f%get_ptr(bu, igrid)
            DO i = 1, ii
                DO j = 1, jj
                    DO k = 1, kk
                        area = ddx(i) * ddy(j)
                        areau(k, j, i) = bu(k, j, i) * area
                    END DO
                END DO
            END DO
            ! Refine for intersected cells
            aucells = SIZE(this%awind(imygrid)%arr)
            DO idx = 1, aucells
                ind = this%awind(imygrid)%arr(idx)
                CALL ind2sub(ind, k, j, i, kk, jj, ii)
                area = ddx(i) * ddy(j)
                areau(k, j, i) = this%awvalue(imygrid)%arr(idx) * area
            END DO

        END DO all_grids
    END SUBROUTINE get_auavaw



    SUBROUTINE calc_volp(this, bzelltyp, areau_f, areav_f, areaw_f, icells, &
        icellspointer, xpsw, volp_f)
        ! Subroutine arguments
        CLASS(gc_stencils_t), INTENT(inout) :: this
        INTEGER(intk), INTENT(in) :: bzelltyp(*)
        TYPE(field_t), INTENT(in) :: areau_f
        TYPE(field_t), INTENT(in) :: areav_f
        TYPE(field_t), INTENT(in) :: areaw_f
        INTEGER(intk), INTENT(in) :: icells(:)
        INTEGER(intk), INTENT(in) :: icellspointer(:)
        REAL(realk), INTENT(in) :: xpsw(:, :)
        TYPE(field_t), INTENT(inout) :: volp_f

        ! Local variables
        INTEGER(intk) :: i, igrid, kk, jj, ii, ip3, ipp, ncells
        REAL(realk), POINTER, CONTIGUOUS :: xstag(:), ystag(:), zstag(:)
        REAL(realk), POINTER, CONTIGUOUS :: ddx(:), ddy(:), ddz(:)
        REAL(realk), POINTER, CONTIGUOUS :: areau(:, :, :)
        REAL(realk), POINTER, CONTIGUOUS :: areav(:, :, :)
        REAL(realk), POINTER, CONTIGUOUS :: areaw(:, :, :)
        REAL(realk), POINTER, CONTIGUOUS :: volp(:, :, :)

        DO i = 1, nmygrids
            igrid = mygrids(i)

            ! Index pointers and dimensions
            CALL get_mgdims(kk, jj, ii, igrid)
            CALL get_ip3(ip3, igrid)
            ipp = icellspointer(igrid)
            ncells = icells(igrid)

            ! Cell side lengths of pressure cells
            CALL get_fieldptr(ddx, "DDX", igrid)
            CALL get_fieldptr(ddy, "DDY", igrid)
            CALL get_fieldptr(ddz, "DDZ", igrid)
            ! Get staggered coordinate positions
            CALL get_fieldptr(xstag, "XSTAG", igrid)
            CALL get_fieldptr(ystag, "YSTAG", igrid)
            CALL get_fieldptr(zstag, "ZSTAG", igrid)

            ! Areas including info on intersected cells
            CALL areau_f%get_ptr(areau, igrid)
            CALL areav_f%get_ptr(areav, igrid)
            CALL areaw_f%get_ptr(areaw, igrid)
            CALL volp_f%get_ptr(volp, igrid)

            CALL calc_volp_grid(kk, jj, ii, ddx, ddy, ddz, &
                xstag, ystag, zstag, bzelltyp(ip3), areau, areav, areaw, &
                xpsw(:, ipp:ipp+ncells-1), volp)
        END DO
    END SUBROUTINE calc_volp


    SUBROUTINE calc_volp_grid(kk, jj, ii, ddx, ddy, ddz, &
            xstag, ystag, zstag, bzelltyp, areau, areav, areaw, xpsw, volp)
        ! Subroutine arguments
        INTEGER(intk), INTENT(in) :: kk, jj, ii
        REAL(realk), INTENT(in) :: ddx(ii), ddy(jj), ddz(kk)
        REAL(realk), INTENT(in) :: xstag(ii), ystag(jj), zstag(kk)
        INTEGER(intk), INTENT(in) :: bzelltyp(kk, jj, ii)
        REAL(realk), INTENT(in) :: areau(kk, jj, ii)
        REAL(realk), INTENT(in) :: areav(kk, jj, ii)
        REAL(realk), INTENT(in) :: areaw(kk, jj, ii)
        REAL(realk), INTENT(in) :: xpsw(:, :)
        REAL(realk), INTENT(inout) :: volp(kk, jj, ii)

        ! Local variables
        INTEGER(intk) :: k, j, i, icell
        REAL(realk) :: swx, swy, swz, areax, areay, areaz
        REAL(realk) :: vpx, vpy, vpz

        ! Now set VOLP ( absolute volumes in unit [L^3] )
        DO i = 2, ii
            DO j = 2, jj
                DO k = 2, kk

                    IF (bzelltyp(k, j, i) == 0 .OR. &
                        bzelltyp(k, j, i) == 1) THEN

                        ! Cell is fully open or fully blocked
                        volp(k, j, i) = REAL(bzelltyp(k, j, i), realk) * &
                            ddx(i) * ddy(j) * ddz(k)

                    ELSE IF (bzelltyp(k, j, i) < 0) THEN

                        ! Cell is an interface cell and intersected
                        icell = -bzelltyp(k, j, i)

                        swx = xpsw(1, icell)
                        swy = xpsw(2, icell)
                        swz = xpsw(3, icell)

                        areax = areau(k, j, i) - areau(k, j, i-1)
                        areay = areav(k, j, i) - areav(k, j-1, i)
                        areaz = areaw(k, j, i) - areaw(k-1, j, i)

                        vpx = (areau(k, j, i-1) * ddx(i)) &
                            + areax * (xstag(i) - swx)

                        vpy = (areav(k, j-1, i) * ddy(j)) &
                            + areay * (ystag(j) - swy)

                        vpz = (areaw(k-1, j, i) * ddz(k)) &
                            + areaz * (zstag(k) - swz)

                        volp(k, j, i) = 1.0 / 3.0 * (vpx + vpy + vpz)

                    ELSE

                        ! Unexpected type indicator
                        WRITE(*, *) "bzelltype > 1 not expected"
                        CALL errr(__FILE__, __LINE__)

                    END IF

                END DO
            END DO
        END DO
    END SUBROUTINE calc_volp_grid


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
