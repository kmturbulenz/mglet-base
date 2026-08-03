MODULE ibmodel_mod
    USE core_mod, ONLY: realk, intk, field_t

    IMPLICIT NONE(type, external)
    PRIVATE

    TYPE, ABSTRACT :: ibmodel_t
        CHARACTER(len=16) :: type
    CONTAINS
        PROCEDURE(blockbp_i), DEFERRED :: blockbp
        PROCEDURE(read_stencils_i), DEFERRED :: read_stencils
        PROCEDURE(giteig_i), NOPASS, DEFERRED :: giteig
        PROCEDURE(divcal_i), DEFERRED :: divcal
    END TYPE ibmodel_t

    ABSTRACT INTERFACE
        SUBROUTINE start_and_stop_i(ista, isto, jsta, jsto, &
                ksta, ksto, ctyp, igrid)
            IMPORT :: intk
            INTEGER(intk), INTENT(out) :: ista, isto, jsta, jsto, ksta, ksto
            CHARACTER(len=1), INTENT(in) :: ctyp
            INTEGER(intk), INTENT(in) :: igrid
        END SUBROUTINE start_and_stop_i

        SUBROUTINE blockbp_i(this, stop_now)
            IMPORT :: ibmodel_t
            CLASS(ibmodel_t), INTENT(inout) :: this
            LOGICAL, INTENT(out) :: stop_now
        END SUBROUTINE blockbp_i

        SUBROUTINE read_stencils_i(this)
            IMPORT :: ibmodel_t
            CLASS(ibmodel_t), INTENT(inout) :: this
        END SUBROUTINE read_stencils_i

        SUBROUTINE giteig_i()
        END SUBROUTINE giteig_i

        SUBROUTINE divcal_i(this, div, u, v, w, fak, ctyp)
            IMPORT :: field_t, ibmodel_t, realk
            CLASS(ibmodel_t), INTENT(inout) :: this
            TYPE(field_t), INTENT(inout) :: div
            TYPE(field_t), INTENT(in) :: u
            TYPE(field_t), INTENT(in) :: v
            TYPE(field_t), INTENT(in) :: w
            REAL(realk), INTENT(in) :: fak
            CHARACTER(len=1), INTENT(in), OPTIONAL :: ctyp
        END SUBROUTINE divcal_i
    END INTERFACE

    PUBLIC :: ibmodel_t
END MODULE ibmodel_mod
