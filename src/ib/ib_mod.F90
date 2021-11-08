MODULE ib_mod
    USE blocknodes_mod
    USE ctof_mod
    USE cutcorner_mod
    USE ftoc_mod
    USE ibconst_mod
    USE ibcore_mod
    USE parent_mod
    USE topol_mod
    USE gc_mod, gc_constructor => constructor
    USE noib_mod, noib_constructor => constructor

    IMPLICIT NONE(type, external)
    ! In this module everything is public. It exports all the public
    ! definitions from the other modules as well, which makes it easier
    ! to use them, you do not need to "USE" so many modules each place.

CONTAINS
    SUBROUTINE init_ib()
        CALL register_ib("ghostcell", gc_constructor)
        CALL register_ib("noib", noib_constructor)

        CALL init_ibcore()
        CALL init_ctof()
        CALL init_ftoc()
        CALL init_parent()

        CALL set_finecell()
    END SUBROUTINE init_ib


    SUBROUTINE finish_ib()
        CALL finish_parent()
        CALL finish_ftoc()
        CALL finish_ctof()
        CALL finish_ibcore()
    END SUBROUTINE finish_ib


    SUBROUTINE set_finecell()
        USE core_mod, ONLY: set_field, get_fieldptr, maxlevel, minlevel, &
            connect, idim3d

        ! Local variables
        INTEGER :: ilevel
        REAL(realk), ALLOCATABLE :: hilf(:)
        REAL(realk), POINTER, CONTIGUOUS :: finecell(:)

        CALL set_field("FINECELL")
        ALLOCATE(hilf(idim3d))
        hilf = 0.0

        CALL get_fieldptr(finecell, "FINECELL")
        finecell = 1.0

        DO ilevel = maxlevel, minlevel, -1
            CALL ftoc(ilevel, hilf, finecell, 'P')
        END DO

        DO ilevel = maxlevel, minlevel, -1
            CALL connect(ilevel, 2, s1=finecell, corners=.TRUE.)
        END DO

        DEALLOCATE(hilf)
    END SUBROUTINE set_finecell
END MODULE ib_mod
