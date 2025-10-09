MODULE gc_mod
    USE core_mod
    USE ibmodel_mod, ONLY: ibmodel_t
    USE ibconst_mod, ONLY: maccur
    USE noib_mod, ONLY: noib_t
    USE bubvbw_mod, ONLY: bubvbw
    USE gc_restrict_mod, ONLY: gc_restrict_t
    USE gc_blockbp_mod, ONLY: gc_blockbp_t
    USE gc_stencils_mod, ONLY: gc_stencils_t

    IMPLICIT NONE(type, external)
    PRIVATE

    TYPE, EXTENDS(noib_t) :: gc_t
        INTEGER(intk) :: ncells
        INTEGER(intk), ALLOCATABLE :: icells(:)
        INTEGER(intk), ALLOCATABLE :: icellspointer(:)

        REAL(realk), ALLOCATABLE :: xpsw(:, :), nvecs(:, :), ucell(:, :)
        REAL(realk), ALLOCATABLE :: arealist(:)
        INTEGER(intk), ALLOCATABLE :: bodyid(:)
        INTEGER(intk), ALLOCATABLE :: bzelltyp(:)

        TYPE(gc_stencils_t) :: stencils
        CHARACTER(len=mglet_filename_max) :: stencilfile
        CHARACTER(len=mglet_filename_max), ALLOCATABLE :: stlnames(:)
    CONTAINS
        PROCEDURE :: blockbp
        PROCEDURE :: read_stencils
        PROCEDURE :: divcal
        PROCEDURE, PRIVATE :: calc_nvecs
        PROCEDURE, NOPASS, PRIVATE :: calc_nvecs_grid
        FINAL :: destructor
    END TYPE gc_t

    PUBLIC :: gc_t, constructor

CONTAINS
    SUBROUTINE constructor(ib)
        ! Subroutine arguments
        CLASS(ibmodel_t), ALLOCATABLE, INTENT(out) :: ib

        ! Local variables
        REAL(realk), POINTER, CONTIGUOUS :: bp(:), bu(:), bv(:), bw(:)
        REAL(realk), POINTER, CONTIGUOUS :: areau(:), areav(:), areaw(:)
        REAL(realk), POINTER, CONTIGUOUS :: volp(:)

        IF (myid == 0) THEN
            WRITE(*, '(A)') "Using 'ghostcell' immersed boundary method"
            WRITE(*, '()')
        END IF

        ALLOCATE(gc_t :: ib)
        ALLOCATE(gc_restrict_t :: ib%restrict_op)

        ib%type = "GHOSTCELL"

        CALL set_field("BP", dwrite=.TRUE.)
        CALL set_field("BU", istag=1 )
        CALL set_field("BV", jstag=1 )
        CALL set_field("BW", kstag=1 )

        CALL set_field("AREAU", istag=1, units=[0, 2, 0, 0, 0, 0, 0])
        CALL set_field("AREAV", jstag=1, units=[0, 2, 0, 0, 0, 0, 0])
        CALL set_field("AREAW", kstag=1, units=[0, 2, 0, 0, 0, 0, 0])
        CALL set_field("VOLP", units=[0, 3, 0, 0, 0, 0, 0], dwrite=.TRUE.)

        CALL set_field("SDIV")

        CALL get_fieldptr(bp, "BP")
        CALL get_fieldptr(bu, "BU")
        CALL get_fieldptr(bv, "BV")
        CALL get_fieldptr(bw, "BW")

        CALL get_fieldptr(areau, "AREAU")
        CALL get_fieldptr(areav, "AREAV")
        CALL get_fieldptr(areaw, "AREAW")
        CALL get_fieldptr(volp, "VOLP")

        bp = 1.0
        bu = 1.0
        bv = 1.0
        bw = 1.0

        ! Negative values to check that values are overwritten later
        areau = -1.0
        areav = -1.0
        areaw = -1.0
        volp = -1.0

        SELECT TYPE(ib)
        TYPE IS (gc_t)
            CALL fort7%get_value("/ib/stencilfile", ib%stencils%file, &
                "ib_stencils.h5")

            ALLOCATE(ib%icells(ngrid))
            ALLOCATE(ib%icellspointer(ngrid))
            ib%icells = 0
            ib%icellspointer = 0
        CLASS DEFAULT
            CALL errr(__FILE__, __LINE__)
        END SELECT
    END SUBROUTINE constructor


    SUBROUTINE destructor(this)
        TYPE(gc_t), INTENT(inout) :: this

        IF (ALLOCATED(this%restrict_op)) THEN
            DEALLOCATE(this%restrict_op)
        END IF
    END SUBROUTINE destructor


    SUBROUTINE blockbp(this, stop_now)
        ! Subroutine arguments
        CLASS(gc_t), INTENT(inout) :: this
        LOGICAL, INTENT(out) :: stop_now

        ! Local variables
        TYPE(gc_blockbp_t) :: blockbp_op

        CALL blockbp_op%blockbp(this%icells, this%icellspointer, &
            this%stencils, stop_now)
    END SUBROUTINE blockbp


    SUBROUTINE read_stencils(this)
        ! Subroutine arguments
        CLASS(gc_t), INTENT(inout) :: this

        ! Local variables
        TYPE(field_t), POINTER :: bp, bu, bv, bw
        TYPE(field_t), POINTER :: areau, areav, areaw
        TYPE(field_t), POINTER :: volp
        TYPE(gc_blockbp_t) :: blockbp_op

        CALL this%stencils%read()

        CALL get_field(bp, "BP")
        CALL this%stencils%get_bp(bp)

        CALL get_field(bu, "BU")
        CALL get_field(bv, "BV")
        CALL get_field(bw, "BW")
        CALL bubvbw(bp, bu, bv, bw)

        CALL get_field(areau, "AREAU")
        CALL get_field(areav, "AREAV")
        CALL get_field(areaw, "AREAW")
        CALL this%stencils%get_auavaw(bu, bv, bw, areau, areav, areaw)

        CALL this%stencils%get_icells(this%icells)
        CALL blockbp_op%set_icellspointer(this%icells, this%icellspointer)

        this%ncells = SUM(this%icells)
        ALLOCATE(this%xpsw(3, this%ncells))
        ALLOCATE(this%ucell(3, this%ncells))
        ALLOCATE(this%bodyid(this%ncells))
        ALLOCATE(this%arealist(this%ncells))
        ALLOCATE(this%bzelltyp(idim3d))

        CALL this%stencils%get_intersected(this%icells, this%icellspointer, &
            this%bodyid, this%xpsw, this%ucell, this%bzelltyp)

        CALL get_field(volp, "VOLP")
        CALL this%stencils%calc_volp(this%bzelltyp, areau, areav, areaw, &
            this%icells, this%icellspointer, this%xpsw, bp, volp)

        ALLOCATE(this%nvecs(4, this%ncells))
        CALL this%calc_nvecs(this%bzelltyp, areau, areav, areaw, this%icells, &
            this%icellspointer, this%xpsw, this%nvecs, this%arealist)

        ALLOCATE(this%stlnames, SOURCE=this%stencils%stlnames)

        CALL this%stencils%finish()
    END SUBROUTINE read_stencils


    SUBROUTINE calc_nvecs(this, bzelltyp, areau, areav, areaw, &
            icells, icellspointer, xpsw, nvecs, arealist)
        ! Subroutine arguments
        CLASS(gc_t), INTENT(inout) :: this
        INTEGER(intk), INTENT(in) :: bzelltyp(*)
        TYPE(field_t), INTENT(in) :: areau
        TYPE(field_t), INTENT(in) :: areav
        TYPE(field_t), INTENT(in) :: areaw
        INTEGER(intk), INTENT(in) :: icells(:)
        INTEGER(intk), INTENT(in) :: icellspointer(:)
        REAL(realk), INTENT(in) :: xpsw(:, :)
        REAL(realk), INTENT(out) :: nvecs(:, :)
        REAL(realk), INTENT(out) :: arealist(:)

        ! Local variables
        INTEGER(intk) :: i, igrid, kk, jj, ii, ip3, ipp, ncells
        REAL(realk), POINTER, CONTIGUOUS :: ddx(:), ddy(:), ddz(:)

        DO i = 1, nmygrids
            igrid = mygrids(i)
            CALL get_mgdims(kk, jj, ii, igrid)
            CALL get_ip3(ip3, igrid)
            ipp = icellspointer(igrid)
            ncells = icells(igrid)

            CALL get_fieldptr(ddx, "DDX", igrid)
            CALL get_fieldptr(ddy, "DDY", igrid)
            CALL get_fieldptr(ddz, "DDZ", igrid)

            CALL this%calc_nvecs_grid(kk, jj, ii, ddx, ddy, ddz, &
                bzelltyp(ip3), areau%arr(ip3), areav%arr(ip3), areaw%arr(ip3), &
                xpsw(:, ipp:ipp+ncells-1), nvecs(:, ipp:ipp+ncells-1), &
                arealist(ipp:ipp+ncells-1))
        END DO
    END SUBROUTINE calc_nvecs


    SUBROUTINE calc_nvecs_grid(kk, jj, ii, ddx, ddy, ddz, bzelltyp, &
            areau, areav, areaw, xpsw, nvecs, arealist)
        ! Subroutine arguments
        INTEGER(intk), INTENT(in) :: kk, jj, ii
        REAL(realk), INTENT(in) :: ddx(ii), ddy(jj), ddz(kk)
        INTEGER(intk), INTENT(in) :: bzelltyp(kk, jj, ii)
        REAL(realk), INTENT(in) :: areau(kk, jj, ii)
        REAL(realk), INTENT(in) :: areav(kk, jj, ii)
        REAL(realk), INTENT(in) :: areaw(kk, jj, ii)
        REAL(realk), INTENT(in) :: xpsw(:, :)
        REAL(realk), INTENT(out) :: nvecs(:, :)
        REAL(realk), INTENT(out) :: arealist(:)

        ! Local variables
        INTEGER(intk) :: k, j, i, icell
        REAL(realk) :: nx, ny, nz, area, cartarea

        ! Corresponding to indices in calcauavaw_mod.F90 where this
        ! code initially came from
        DO i = 2, ii
            DO j = 2, jj
                DO k = 2, kk
                    IF (bzelltyp(k, j, i) >= 0) CYCLE
                    icell = -bzelltyp(k, j, i)

                    area = SQRT((areau(k, j, i)-areau(k, j, i-1))**2 &
                              + (areav(k, j, i)-areav(k, j-1, i))**2 &
                              + (areaw(k, j, i)-areaw(k-1, j, i))**2)

                    ! From calcauavaw
                    ! Referenzflaeche: Mittel der kartesischen
                    ! Flaechen aus allen drei Richtungen.
                    cartarea = (1.0/3.0) * (ddx(i)*ddy(j) + &
                        ddx(i)*ddz(k) + ddz(k)*ddy(j))

                    ! Kann vorkommen, dass Flaeche zu klein,
                    ! dann muss Notfallnormale berechnet werden.
                    IF (area/cartarea < maccur**2) THEN
                        nx = 0.0
                        ny = 0.0
                        nz = 0.0
                    ELSE
                        nx = (areau(k, j, i)-areau(k, j, i-1)) / area
                        ny = (areav(k, j, i)-areav(k, j-1, i)) / area
                        nz = (areaw(k, j, i)-areaw(k-1, j, i)) / area
                    END IF

                    arealist(icell) = area
                    nvecs(1, icell) = nx
                    nvecs(2, icell) = ny
                    nvecs(3, icell) = nz
                    nvecs(4, icell) = nx*xpsw(1, icell) &
                        + ny*xpsw(2, icell) &
                        + nz*xpsw(3, icell)
                END DO
            END DO
        END DO
    END SUBROUTINE calc_nvecs_grid


    SUBROUTINE divcal(this, div, u, v, w, fak, ctyp)
        ! Subroutine arguments
        CLASS(gc_t), INTENT(inout) :: this
        TYPE(field_t), INTENT(inout) :: div
        TYPE(field_t), INTENT(in) :: u
        TYPE(field_t), INTENT(in) :: v
        TYPE(field_t), INTENT(in) :: w
        REAL(realk), INTENT(in) :: fak
        CHARACTER(len=1), INTENT(in), OPTIONAL :: ctyp

        ! Local variables
        INTEGER(intk) :: i, igrid, ilevel
        INTEGER(intk) :: kk, jj, ii
        LOGICAL :: use_sdiv
        TYPE(field_t), POINTER :: rddx_f, rddy_f, rddz_f
        TYPE(field_t), POINTER :: sdiv_f, bp_f
        REAL(realk), CONTIGUOUS, POINTER :: rddx(:), rddy(:), rddz(:)
        REAL(realk), CONTIGUOUS, POINTER :: sdiv(:, :, :), bp(:, :, :), &
            div_p(:, :, :), u_p(:, :, :), v_p(:, :, :), w_p(:, :, :)

        CALL start_timer(240)

        ! For safety
        NULLIFY(bp_f)
        NULLIFY(sdiv_f)
        NULLIFY(sdiv)

        CALL get_field(rddx_f, "RDDX")
        CALL get_field(rddy_f, "RDDY")
        CALL get_field(rddz_f, "RDDZ")
        CALL get_field(bp_f, "BP")

        use_sdiv = .TRUE.
        IF (PRESENT(ctyp)) THEN
            IF (ctyp == "W") THEN
                use_sdiv = .FALSE.
            END IF
        END IF

        IF (use_sdiv) THEN
            CALL get_field(sdiv_f, "SDIV")
        END IF

        DO ilevel = minlevel, maxlevel
            ! Assume that U, V, W and DIV are defined on the same levels!!!
            IF (.NOT. u%active_level(ilevel)) CYCLE

            DO i = 1, nmygridslvl(ilevel)
                igrid = mygridslvl(i, ilevel)
                CALL get_mgdims(kk, jj, ii, igrid)

                CALL rddx_f%get_ptr(rddx, igrid)
                CALL rddy_f%get_ptr(rddy, igrid)
                CALL rddz_f%get_ptr(rddz, igrid)
                CALL bp_f%get_ptr(bp, igrid)

                CALL div%get_ptr(div_p, igrid)
                CALL u%get_ptr(u_p, igrid)
                CALL v%get_ptr(v_p, igrid)
                CALL w%get_ptr(w_p, igrid)

                IF (use_sdiv) CALL sdiv_f%get_ptr(sdiv, igrid)

                CALL this%divcal_grid(kk, jj, ii, fak, div_p, &
                    u_p, v_p, w_p, rddx, rddy, rddz, bp, sdiv)
            END DO
        END DO

        CALL stop_timer(240)
    END SUBROUTINE divcal

END MODULE gc_mod
