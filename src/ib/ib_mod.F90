MODULE ib_mod
    USE blockcheck_mod, ONLY: blockcheck
    USE blocknodes_mod
    USE calcauavaw_mod, ONLY: calcauavaw
    USE checkzelle_mod, ONLY: checkzelle
    USE ctof_mod
    USE cutcorner_mod
    USE filling_mod, ONLY: fillfluid
    USE flzelle_mod
    USE freekante_mod, ONLY: freekante
    USE freepressure_mod, ONLY: freepressure
    USE ftoc_mod
    USE gc_blockbpfeld_mod, ONLY: blockbpfeld
    USE gc_blockbp_mod, ONLY: blockluecken_closetoboundary
    USE gc_blockface_mod, ONLY: blockface
    USE gc_totwasser_mod, ONLY: totwasser
    USE gc_zelltyp_mod, ONLY: zelltyp
    USE ibconst_mod
    USE ibcore_mod
    USE knotenundkanten_mod
    USE openbubvbw_mod, ONLY: openbubvbw
    USE parent_mod
    USE punktekoordinaten_mod
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
        USE core_mod

        ! Local variables
        INTEGER :: ilevel
        REAL(realk), ALLOCATABLE :: hilf(:)
        TYPE(field_t), POINTER :: finecell

        CALL set_field("FINECELL")
        CALL get_field(finecell, "FINECELL")
        finecell%arr = 1.0

        ALLOCATE(hilf(SIZE(finecell%arr)))
        hilf = 0.0

        DO ilevel = maxlevel, minlevel, -1
            CALL ftoc(ilevel, hilf, finecell%arr, 'P')
        END DO

        DO ilevel = maxlevel, minlevel, -1
            CALL connect(ilevel, 2, s1=finecell, corners=.TRUE.)
        END DO

        DEALLOCATE(hilf)
    END SUBROUTINE set_finecell
END MODULE ib_mod
