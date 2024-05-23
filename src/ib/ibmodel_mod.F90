MODULE ibmodel_mod
    USE core_mod, ONLY: realk, intk, field_t
    USE blockbp_mod, ONLY: blockbp_t

    IMPLICIT NONE(type, external)
    PRIVATE

    TYPE, ABSTRACT :: restrict_t
    CONTAINS
        PROCEDURE :: message_length
        PROCEDURE(start_and_stop_i), NOPASS, DEFERRED :: start_and_stop
        PROCEDURE(restrict_i), DEFERRED :: restrict
    END TYPE restrict_t

    TYPE, ABSTRACT :: ibmodel_t
        CHARACTER(len=16) :: type
        CLASS(restrict_t), ALLOCATABLE :: restrict_op
        CLASS(blockbp_t), ALLOCATABLE :: blockbp_op
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

        SUBROUTINE restrict_i(this, kk, jj, ii, ff, sendbuf, ctyp, igrid)
            IMPORT :: restrict_t, intk, realk
            CLASS(restrict_t), INTENT(inout) :: this
            INTEGER(intk), INTENT(IN) :: kk, jj, ii
            REAL(realk), INTENT(IN) :: ff(kk, jj, ii)
            REAL(realk), CONTIGUOUS, INTENT(INOUT) :: sendbuf(:)
            CHARACTER(len=1), INTENT(in) :: ctyp
            INTEGER(intk), INTENT(IN) :: igrid
        END SUBROUTINE restrict_i

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

    PUBLIC :: ibmodel_t, restrict_t

CONTAINS
    INTEGER(intk) FUNCTION message_length(this, ctyp, igrid)
        ! Subroutine arguments
        CLASS(restrict_t), INTENT(inout) :: this
        CHARACTER(len=1), INTENT(in) :: ctyp
        INTEGER(intk), INTENT(IN) :: igrid

        ! Local variables
        INTEGER(intk) :: istart, istop, jstart, jstop, kstart, kstop

        CALL this%start_and_stop(istart, istop, jstart, jstop, &
            kstart, kstop, ctyp, igrid)

        message_length = ((istop-istart)/2+1) &
            *((jstop-jstart)/2+1) &
            *((kstop-kstart)/2+1)
    END FUNCTION message_length
END MODULE ibmodel_mod
