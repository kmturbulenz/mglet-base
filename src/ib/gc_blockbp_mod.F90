MODULE gc_blockbp_mod
    USE MPI_f08, ONLY: MPI_Wtime

    USE core_mod
    USE calcauavaw_mod, ONLY: calcauavaw
    USE checkzelle_mod, ONLY: checkzelle
    USE cutcorner_mod, ONLY: cutcorner
    USE blocknodes_mod, ONLY: blocknodes
    USE blockcheck_mod, ONLY: blockcheck
    USE blockclosetoboundary_mod, ONLY: blockclosetoboundary
    USE blockluecken_mod, ONLY: blockluecken, blockluecken_grid
    USE filling_mod, ONLY: fillfluid, blockparentboundary_p
    USE fillfollower_mod, ONLY: fillfollower
    USE freepressure_mod, ONLY: freepressure
    USE freekante_mod, ONLY: freekante
    USE gc_zelltyp_mod, ONLY: zelltyp
    USE gc_blockface_mod, ONLY: blockface
    USE gc_blockbpfeld_mod, ONLY: blockbpfeld
    USE openbubvbw_mod, ONLY: openbubvbw
    USE topol_mod, ONLY: topol_t
    USE blockbp_mod, ONLY: blockbp_t
    USE parent_mod, ONLY: parent
    USE ftoc_mod, ONLY: ftoc
    USE gc_totwasser_mod, ONLY: totwasser
    USE bubvbw_mod, ONLY: bubvbw
    USE gc_finishknotenbezelltyp_mod, ONLY: finishknotenbezelltyp
    USE calcnormals_mod, ONLY: calcnormals
    USE gc_stencils_mod, ONLY: gc_stencils_t

    IMPLICIT NONE (type, external)
    PRIVATE

    TYPE, EXTENDS(blockbp_t) :: gc_blockbp_t
    CONTAINS
        PROCEDURE :: blockbp
    END TYPE gc_blockbp_t

    PUBLIC :: gc_blockbp_t, list_to_field, blockluecken_closetoboundary

CONTAINS
    SUBROUTINE blockbp(this, icells, icellspointer, stencils)
        ! Subroutine arguments
        CLASS(gc_blockbp_t), INTENT(inout) :: this
        INTEGER(intk), INTENT(inout) :: icells(:)
        INTEGER(intk), INTENT(inout) :: icellspointer(:)
        TYPE(gc_stencils_t), INTENT(inout) :: stencils

        ! Local variables
        INTEGER(intk), PARAMETER :: ntrimax = 2
        INTEGER(intk) :: ncells

        INTEGER(intk), ALLOCATABLE :: triau(:), triav(:), triaw(:), bzelltyp(:)
        REAL(realk), ALLOCATABLE :: xpsw(:, :), ucell(:, :)
        INTEGER(intk), ALLOCATABLE :: bodyid(:)
        REAL(realk), ALLOCATABLE :: velocity(:, :)

        TYPE(field_t), POINTER :: bp, bu, bv, bw
        TYPE(field_t) :: knoten
        TYPE(field_t) :: au, av, aw
        TYPE(field_t) :: kanteu, kantev, kantew

        ! Read configuration etc.
        CALL this%init()

        IF (.NOT. this%do_blocking) THEN
            RETURN
        END IF

        icells = 0
        icellspointer = 0

        ALLOCATE(triau(ntrimax*idim3d))
        ALLOCATE(triav(ntrimax*idim3d))
        ALLOCATE(triaw(ntrimax*idim3d))
        ALLOCATE(bzelltyp(idim3d))

        CALL kanteu%init("KANTEU", jstag=1, kstag=1)
        CALL kantev%init("KANTEV", istag=1, kstag=1)
        CALL kantew%init("KANTEW", istag=1, jstag=1)

        CALL au%init("AU", istag=1)
        CALL av%init("AV", jstag=1)
        CALL aw%init("AW", kstag=1)

        ! Important with proper staggering for parent to work
        CALL knoten%init("KNOTEN", istag=1, jstag=1, kstag=1)

        IF (myid == 0) THEN
           WRITE (*,'("BLOCKING: ", A40, ", WTIME:", F16.3)') &
               'FIND CELLS INTERSECTED BY TRIANGLES', MPI_Wtime() - this%time0
        END IF

        CALL cutcorner(this%topol, ntrimax, kanteu, kantev, &
            kantew, triau, triav, triaw)
        CALL blocknodes(kanteu, kantev, kantew, knoten)

        IF (myid == 0) THEN
           WRITE (*,'("BLOCKING: ", A40, ", WTIME:", F16.3)') &
               'FILLING OF FLUID CELLS', MPI_Wtime() - this%time0
        END IF

        CALL fillfluid(this%fluidpoints, knoten)

        IF (this%nnofluidpoints > 0) THEN
            IF (myid == 0) THEN
               WRITE (*,'("BLOCKING: ", A40, ", WTIME:", F16.3)') &
                   'FILLFOLLOWER', MPI_Wtime() - this%time0
            END IF
            CALL fillfollower(knoten, this%nofluidpoints)
        END IF

        IF (myid == 0) THEN
           WRITE (*,'("BLOCKING: ", A40, ", WTIME:", F16.3)') &
               'FINISHING INTERSECTED CORNERS', MPI_Wtime() - this%time0
        END IF

        CALL freepressure(kanteu, kantev, kantew, knoten)

        CALL freekante(this%topol, ntrimax, knoten, kanteu, kantev, &
            kantew, triau, triav, triaw)

        CALL blockcheck(kanteu, kantev, kantew, knoten, 0)

        CALL checkzelle(this%topol, ntrimax, knoten, kanteu, kantev, &
            kantew, triau, triav, triaw)

        CALL zelltyp(knoten, bzelltyp, icells)
        CALL this%set_icellspointer(icells, icellspointer)

        CALL get_field(bu, "BU")
        CALL get_field(bv, "BV")
        CALL get_field(bw, "BW")
        CALL blockface(knoten, bu, bv, bw)

        ncells = SUM(icells)
        ALLOCATE(xpsw(3, ncells))
        CALL calcauavaw(this%topol, ntrimax, triau, triav, triaw, &
            knoten, kanteu, kantev, kantew, bzelltyp, au, av, aw, &
            icells, icellspointer, xpsw)

        CALL openbubvbw(au, av, aw, bu, bv, bw)

        CALL get_field(bp, "BP")
        CALL blockbpfeld(bu, bv, bw, bp)

        CALL blockluecken_closetoboundary(bp)

        IF (myid == 0) THEN
           WRITE (*,'("BLOCKING: ", A40, ", WTIME:", F16.3)') &
               'REMOVING DEAD WATER ZONES', MPI_Wtime() - this%time0
        END IF

        CALL totwasser(this%fluidpoints, bp)
        CALL bubvbw(bp, bu, bv, bw)

        ! BP field is now finished
        IF (myid == 0) THEN
           WRITE (*,'("BLOCKING: ", A40, ", WTIME:", F16.3)') &
               'WRITE GEOMETRY DATA, CREATE STENCILS', MPI_Wtime() - this%time0
        END IF

        CALL finishknotenbezelltyp(this%topol, ntrimax, triau, triav, triaw, &
            bzelltyp, bp, knoten, kanteu, kantev, kantew, icells)
        CALL this%set_icellspointer(icells, icellspointer)

        ! Re-allocate storage
        ncells = SUM(icells)
        DEALLOCATE(xpsw)
        ALLOCATE(xpsw(3, ncells))

        CALL calcauavaw(this%topol, ntrimax, triau, triav, triaw, &
            knoten, kanteu, kantev, kantew, bzelltyp, au, av, aw, &
            icells, icellspointer, xpsw)

        CALL this%read_velocity(velocity)
        ALLOCATE(ucell(3, ncells))
        ALLOCATE(bodyid(ncells))

        CALL calcnormals(velocity, this%topol, ntrimax, triau, triav, triaw, &
            knoten, kanteu, kantev, kantew, bzelltyp, au, av, aw, icells, &
            icellspointer, ncells, bodyid, ucell, stencils, .TRUE.)

        ! Deallocate fields that are no longer needed
        CALL kanteu%finish()
        CALL kantev%finish()
        CALL kantew%finish()
        DEALLOCATE(triau)
        DEALLOCATE(triav)
        DEALLOCATE(triaw)
        CALL knoten%finish()

        CALL stencils%set_stlnames(this%topol%geometries)
        CALL this%topol%finish()

        CALL update_bodyid(icells, icellspointer, bzelltyp, bodyid)
        CALL stencils%set_bp(bp)
        CALL stencils%set_auavaw(bu, bv, bw, au, av, aw)
        CALL stencils%set_intersected(icells, icellspointer, &
            bodyid, xpsw, ucell, bzelltyp)

        IF (myid == 0) THEN
            WRITE (*,'("BLOCKING: ", A40, ", WTIME:", F16.3)') &
                'WRITE STENCILS', MPI_Wtime() - this%time0
        END IF

        ! Write out stencilfile
        CALL stencils%write()
        CALL stencils%finish()

        DEALLOCATE(xpsw)
        DEALLOCATE(ucell)
        DEALLOCATE(bodyid)
        DEALLOCATE(bzelltyp)

        CALL au%finish()
        CALL av%finish()
        CALL aw%finish()

        IF (myid == 0) THEN
           WRITE (*,'("BLOCKING: ", A40, ", WTIME:", F16.3)') &
               'BLOCKING COMPLETE', MPI_Wtime() - this%time0
        END IF

        CALL this%finish()

        IF (myid == 0)  WRITE(*, '()')
    END SUBROUTINE blockbp


    SUBROUTINE blockluecken_closetoboundary(bp)
        ! Subroutine arguments
        TYPE(field_t), INTENT(inout) :: bp

        ! Local variables
        INTEGER(intk) :: iloop, ilevel, igrid, i, kk, jj, ii, ip3
        INTEGER(intk) :: nfro, nbac, nrgt, nlft, nbot, ntop

        ! TODO: Maybe move to  gc_mod.F90 where BP is initialized
        CALL bp%init_buffers()

        CALL blockluecken(.FALSE., bp%arr)

        ! TODO: change - variable loop iterations and check for convergence
        DO iloop = 1, 3
            DO ilevel = minlevel, maxlevel
                CALL connect(ilevel, 2, s1=bp)
                CALL parent(ilevel, s1=bp)

                DO i = 1, nmygridslvl(ilevel)
                    igrid = mygridslvl(i, ilevel)

                    CALL get_mgdims(kk, jj, ii, igrid)
                    CALL get_ip3(ip3, igrid)
                    CALL get_mgbasb(nfro, nbac, nrgt, nlft, nbot, ntop, igrid)

                    IF (ilevel > minlevel) THEN
                        CALL blockparentboundary_p(bp, igrid)
                        CALL blockclosetoboundary(kk, jj, ii, nfro, nbac, &
                            nrgt, nlft, nbot, ntop, .FALSE., bp%arr(ip3))
                    END IF
                    CALL blockluecken_grid(kk, jj, ii, nfro, nbac, &
                        nrgt, nlft, nbot, ntop, bp%arr(ip3))
                END DO
            END DO

            DO ilevel = maxlevel, minlevel+1, -1
                CALL ftoc(ilevel, bp%arr, bp%arr, 'E')
            END DO
        END DO
    END SUBROUTINE blockluecken_closetoboundary


    SUBROUTINE update_bodyid(icells, icellspointer, bzelltyp, bodyid)
        ! Subroutine arguments
        INTEGER(intk), INTENT(in) :: icells(:)
        INTEGER(intk), INTENT(in) :: icellspointer(:)
        INTEGER(intk), INTENT(in) :: bzelltyp(*)
        INTEGER(intk), INTENT(inout) :: bodyid(:)

        ! Local variables
        INTEGER(intk) :: ilevel, igrid, i, kk, jj, ii, ip3, ipp, ncells
        INTEGER(intk) :: nfro, nbac, nrgt, nlft, nbot, ntop
        REAL(realk), POINTER, CONTIGUOUS :: bodyid_3d(:, :, :)
        TYPE(field_t) :: bodyid_3d_f

        CALL bodyid_3d_f%init("BODYID")
        CALL bodyid_3d_f%init_buffers()

        ! Copy from list to 3D field
        DO ilevel = minlevel, maxlevel
            DO i = 1, nmygridslvl(ilevel)
                igrid = mygridslvl(i, ilevel)

                CALL get_mgdims(kk, jj, ii, igrid)
                CALL get_ip3(ip3, igrid)
                CALL bodyid_3d_f%get_ptr(bodyid_3d, igrid)

                ipp = icellspointer(igrid)
                ncells = icells(igrid)

                CALL list_to_field(kk, jj, ii, bzelltyp(ip3), &
                    bodyid(ipp:ipp+ncells-1), bodyid_3d)
            END DO
        END DO

        ! Update parent boundaries in field
        DO ilevel = minlevel, maxlevel
            CALL parent(ilevel, s1=bodyid_3d_f)

            DO i = 1, nmygridslvl(ilevel)
                igrid = mygridslvl(i, ilevel)

                CALL get_mgdims(kk, jj, ii, igrid)
                CALL get_ip3(ip3, igrid)
                CALL get_mgbasb(nfro, nbac, nrgt, nlft, nbot, ntop, igrid)

                CALL blockparentboundary_p(bodyid_3d_f, igrid)
            END DO
        END DO

        ! Fine-to-coarse
        DO ilevel = maxlevel, minlevel+1, -1
            CALL ftoc(ilevel, bodyid_3d_f%arr, bodyid_3d_f%arr, 'I')
        END DO

        ! Connect
        DO ilevel = minlevel, maxlevel
            CALL connect(ilevel, 2, s1=bodyid_3d_f, corners=.TRUE.)
        END DO

        ! Copy back to list
        DO ilevel = minlevel, maxlevel
            DO i = 1, nmygridslvl(ilevel)
                igrid = mygridslvl(i, ilevel)

                CALL get_mgdims(kk, jj, ii, igrid)
                CALL get_ip3(ip3, igrid)
                CALL bodyid_3d_f%get_ptr(bodyid_3d, igrid)

                ipp = icellspointer(igrid)
                ncells = icells(igrid)

                CALL field_to_list(kk, jj, ii, bzelltyp(ip3), &
                    bodyid(ipp:ipp+ncells-1), bodyid_3d)
            END DO
        END DO

        CALL bodyid_3d_f%finish()
    END SUBROUTINE update_bodyid


    SUBROUTINE list_to_field(kk, jj, ii, bzelltyp, list, field)
        ! Subroutine arguments
        INTEGER(intk), INTENT(in) :: kk, jj, ii
        INTEGER(intk), INTENT(in) :: bzelltyp(kk, jj, ii)
        INTEGER(intk), INTENT(in) :: list(:)
        REAL(realk), INTENT(out) :: field(kk, jj, ii)

        ! Local variables
        INTEGER(intk) :: k, j, i, icell

        field = 0.0

        DO i = 2, ii
            DO j = 2, jj
                DO k = 2, kk
                    IF (bzelltyp(k, j, i) < 0) THEN
                        icell = -bzelltyp(k, j, i)
                        field(k, j, i) = REAL(list(icell), realk)
                    END IF
                END DO
            END DO
        END DO
    END SUBROUTINE list_to_field


    SUBROUTINE field_to_list(kk, jj, ii, bzelltyp, list, field)
        ! Subroutine arguments
        INTEGER(intk), INTENT(in) :: kk, jj, ii
        INTEGER(intk), INTENT(in) :: bzelltyp(kk, jj, ii)
        INTEGER(intk), INTENT(out) :: list(:)
        REAL(realk), INTENT(in) :: field(kk, jj, ii)

        ! Local variables
        INTEGER(intk) :: k, j, i, icell

        list = 0.0

        DO i = 2, ii
            DO j = 2, jj
                DO k = 2, kk
                    IF (bzelltyp(k, j, i) < 0) THEN
                        icell = -bzelltyp(k, j, i)
                        list(icell) = NINT(field(k, j, i), intk)
                    END IF
                END DO
            END DO
        END DO
    END SUBROUTINE field_to_list
END MODULE gc_blockbp_mod
