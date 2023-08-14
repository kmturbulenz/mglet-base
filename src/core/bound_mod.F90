MODULE bound_mod
    USE precision_mod, ONLY: realk, intk
    USE field_mod, ONLY: field_t
    USE grids_mod, ONLY: nmygridslvl, mygridslvl, nboconds, get_bc_ctyp

    IMPLICIT NONE(type, external)
    PRIVATE

    TYPE, ABSTRACT :: bound_t
    CONTAINS
        PROCEDURE :: bound
        PROCEDURE(bound_face_i), DEFERRED, NOPASS :: front
        PROCEDURE(bound_face_i), DEFERRED, NOPASS :: back
        PROCEDURE(bound_face_i), DEFERRED, NOPASS :: right
        PROCEDURE(bound_face_i), DEFERRED, NOPASS :: left
        PROCEDURE(bound_face_i), DEFERRED, NOPASS :: bottom
        PROCEDURE(bound_face_i), DEFERRED, NOPASS :: top
    END TYPE bound_t

    ABSTRACT INTERFACE
        SUBROUTINE bound_face_i(igrid, iface, ibocd, ctyp, f1, f2, f3, f4, timeph)
            IMPORT :: bound_t, intk, realk, field_t
            INTEGER(intk), INTENT(in) :: igrid, iface, ibocd
            CHARACTER(len=*), INTENT(in) :: ctyp
            TYPE(field_t), INTENT(inout) :: f1
            TYPE(field_t), INTENT(inout), OPTIONAL :: f2, f3, f4
            REAL(realk), INTENT(in), OPTIONAL :: timeph
        END SUBROUTINE bound_face_i
    END INTERFACE

    PUBLIC :: bound_t
CONTAINS

    SUBROUTINE bound(this, ilevel, f1, f2, f3, f4, timeph)
        CLASS(bound_t) :: this
        INTEGER(intk), INTENT(in) :: ilevel
        TYPE(field_t), INTENT(inout) :: f1
        TYPE(field_t), INTENT(inout), OPTIONAL :: f2, f3, f4
        REAL(realk), INTENT(in), OPTIONAL :: timeph

        INTEGER(intk) :: i, igrid, iface, nbocd, ibocd
        CHARACTER(len=8) :: ctyp

        DO i = 1, nmygridslvl(ilevel)
            igrid = mygridslvl(i, ilevel)

            DO iface = 1, 6
                nbocd = nboconds(iface, igrid)
                DO ibocd = 1, nbocd
                    CALL get_bc_ctyp(ctyp, ibocd, iface, igrid)
                    SELECT CASE(iface)
                    CASE(1)
                        CALL this%front(igrid, iface, ibocd, ctyp, &
                            f1, f2, f3, f4, timeph)
                    CASE(2)
                        CALL this%back(igrid, iface, ibocd, ctyp, &
                            f1, f2, f3, f4, timeph)
                    CASE(3)
                        CALL this%right(igrid, iface, ibocd, ctyp, &
                            f1, f2, f3, f4, timeph)
                    CASE(4)
                        CALL this%left(igrid, iface, ibocd, ctyp, &
                            f1, f2, f3, f4, timeph)
                    CASE(5)
                        CALL this%bottom(igrid, iface, ibocd, ctyp, &
                            f1, f2, f3, f4, timeph)
                    CASE(6)
                        CALL this%top(igrid, iface, ibocd, ctyp, &
                            f1, f2, f3, f4, timeph)
                    END SELECT
                END DO
            END DO
        END DO
    END SUBROUTINE bound
END MODULE bound_mod
