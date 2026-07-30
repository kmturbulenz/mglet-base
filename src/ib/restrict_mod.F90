MODULE restrict_mod
    USE core_mod
    USE ibmodel_mod, ONLY: ibmodel_t
    USE noib_restrict_mod
    USE gc_restrict_mod

    IMPLICIT NONE (type, external)
    PRIVATE

    INTERFACE message_length
        MODULE PROCEDURE :: message_length_A
        MODULE PROCEDURE :: message_length_B
    END INTERFACE message_length

    PUBLIC :: restrict, message_length, start_and_stop
CONTAINS
    SUBROUTINE restrict(ib, ctyp, field, igrid, offset)
        ! Subroutine arguments
        CLASS(ibmodel_t), INTENT(in) :: ib
        CHARACTER(len=1), INTENT(in) :: ctyp
        TYPE(field_t), INTENT(in) :: field
        INTEGER(intk), INTENT(in) :: igrid
        INTEGER(int32), INTENT(inout) :: offset

        SELECT CASE (ib%type)
        CASE ("GHOSTCELL")
            CALL restrict_gc(ctyp, field, igrid, offset)
        CASE ("NONE")
            CALL restrict_noib(ctyp, field, igrid, offset)
        CASE DEFAULT
            CALL errr(__FILE__, __LINE__)
        END SELECT
    END SUBROUTINE restrict


    SUBROUTINE restrict_noib(ctyp, f_t, igrid, offset)
        ! Subroutine arguments
        CHARACTER(len=1), INTENT(in) :: ctyp
        TYPE(field_t), INTENT(in) :: f_t
        INTEGER(intk), INTENT(in) :: igrid
        INTEGER(int32), INTENT(inout) :: offset

        ! Local variables
        TYPE(field_t), POINTER :: ddx_f, ddy_f, ddz_f
        REAL(realk), POINTER, CONTIGUOUS, DIMENSION(:, :, :) :: f
        REAL(realk), POINTER, CONTIGUOUS, DIMENSION(:) :: ddx, ddy, ddz

        INTEGER(intk) :: kk, jj, ii
        INTEGER(intk) :: istart, istop, jstart, jstop, kstart, kstop
        INTEGER(int32) :: messagelength, targetidx

        CALL get_mgdims(kk, jj, ii, igrid)
        CALL start_and_stop(istart, istop, jstart, jstop, kstart, kstop, ctyp, &
            igrid)
        CALL message_length(messagelength, istart, istop, jstart, jstop, &
            kstart, kstop)
        targetidx = offset + messagelength - 1

        CALL get_field(ddx_f, "DDX")
        CALL get_field(ddy_f, "DDY")
        CALL get_field(ddz_f, "DDZ")

        CALL ddx_f%get_ptr(ddx, igrid)
        CALL ddy_f%get_ptr(ddy, igrid)
        CALL ddz_f%get_ptr(ddz, igrid)

        CALL f_t%get_ptr(f, igrid)

        SELECT CASE (ctyp)
        CASE ("U")
            CALL restrict_noib_u(kk, jj, ii, f, ddx, ddy, ddz, &
                messagelength, sendbuf(offset:targetidx), &
                istart, istop, jstart, jstop, kstart, kstop)
        CASE ("V")
            CALL restrict_noib_v(kk, jj, ii, f, ddx, ddy, ddz, &
                messagelength, sendbuf(offset:targetidx), &
                istart, istop, jstart, jstop, kstart, kstop)
        CASE ("W")
            CALL restrict_noib_w(kk, jj, ii, f, ddx, ddy, ddz, &
                messagelength, sendbuf(offset:targetidx), &
                istart, istop, jstart, jstop, kstart, kstop)
        CASE ("P", "R", "S", "T")
            CALL restrict_noib_s(kk, jj, ii, f, ddx, ddy, ddz, &
                messagelength, sendbuf(offset:targetidx), &
                istart, istop, jstart, jstop, kstart, kstop)
        CASE DEFAULT
            CALL errr(__FILE__, __LINE__)
        END SELECT

        offset = offset + messagelength
    END SUBROUTINE restrict_noib


    SUBROUTINE restrict_gc(ctyp, f_t, igrid, offset)
        ! Subroutine arguments
        CHARACTER(len=1), INTENT(in) :: ctyp
        TYPE(field_t), INTENT(in) :: f_t
        INTEGER(intk), INTENT(in) :: igrid
        INTEGER(int32), INTENT(inout) :: offset

        ! Local variables
        TYPE(field_t), POINTER :: ddx_f, ddy_f, ddz_f, b_f
        REAL(realk), POINTER, CONTIGUOUS, DIMENSION(:, :, :) :: f, b
        REAL(realk), POINTER, CONTIGUOUS, DIMENSION(:) :: ddx, ddy, ddz

        INTEGER(intk) :: kk, jj, ii
        INTEGER(intk) :: istart, istop, jstart, jstop, kstart, kstop
        INTEGER(int32) :: messagelength, targetidx

        CALL get_mgdims(kk, jj, ii, igrid)
        CALL start_and_stop(istart, istop, jstart, jstop, kstart, kstop, ctyp, &
            igrid)
        CALL message_length(messagelength, istart, istop, jstart, jstop, &
            kstart, kstop)
        targetidx = offset + messagelength - 1

        CALL get_field(ddx_f, "DDX")
        CALL get_field(ddy_f, "DDY")
        CALL get_field(ddz_f, "DDZ")

        CALL ddx_f%get_ptr(ddx, igrid)
        CALL ddy_f%get_ptr(ddy, igrid)
        CALL ddz_f%get_ptr(ddz, igrid)

        CALL f_t%get_ptr(f, igrid)

        SELECT CASE (ctyp)
        CASE ("U")
            CALL restrict_noib_u(kk, jj, ii, f, ddx, ddy, ddz, &
                messagelength, sendbuf(offset:targetidx), &
                istart, istop, jstart, jstop, kstart, kstop)
        CASE ("V")
            CALL restrict_noib_v(kk, jj, ii, f, ddx, ddy, ddz, &
                messagelength, sendbuf(offset:targetidx), &
                istart, istop, jstart, jstop, kstart, kstop)
        CASE ("W")
            CALL restrict_noib_w(kk, jj, ii, f, ddx, ddy, ddz, &
                messagelength, sendbuf(offset:targetidx), &
                istart, istop, jstart, jstop, kstart, kstop)
        CASE ("P")
            CALL get_field(b_f, "BP")
            CALL b_f%get_ptr(b, igrid)
            CALL restrict_gc_p_t(kk, jj, ii, f, ddx, ddy, ddz, b, &
                messagelength, sendbuf(offset:targetidx), &
                istart, istop, jstart, jstop, kstart, kstop)
        CASE ("T")
            CALL get_field(b_f, "BT")
            CALL b_f%get_ptr(b, igrid)
            CALL restrict_gc_p_t(kk, jj, ii, f, ddx, ddy, ddz, b, &
                messagelength, sendbuf(offset:targetidx), &
                istart, istop, jstart, jstop, kstart, kstop)
        CASE ("R")
            CALL get_field(b_f, "BP")
            CALL b_f%get_ptr(b, igrid)
            CALL restrict_gc_r(kk, jj, ii, f, ddx, ddy, ddz, b, &
                messagelength, sendbuf(offset:targetidx), &
                istart, istop, jstart, jstop, kstart, kstop)
        CASE ("S")
            CALL restrict_noib_s(kk, jj, ii, f, ddx, ddy, ddz, &
                messagelength, sendbuf(offset:targetidx), &
                istart, istop, jstart, jstop, kstart, kstop)
        CASE ("N")
            CALL restrict_gc_n(kk, jj, ii, f, messagelength, &
                sendbuf(offset:targetidx), &
                istart, istop, jstart, jstop, kstart, kstop)
        CASE ("E")
            CALL restrict_gc_e(kk, jj, ii, f, messagelength, &
                sendbuf(offset:targetidx), &
                istart, istop, jstart, jstop, kstart, kstop)
        CASE ("F")
            CALL restrict_gc_f(kk, jj, ii, f, messagelength, &
                sendbuf(offset:targetidx), &
                istart, istop, jstart, jstop, kstart, kstop)
        CASE ("I")
            CALL restrict_gc_i(kk, jj, ii, f, messagelength, &
                sendbuf(offset:targetidx), &
                istart, istop, jstart, jstop, kstart, kstop)
        CASE DEFAULT
            CALL errr(__FILE__, __LINE__)
        END SELECT

        offset = offset + messagelength
    END SUBROUTINE restrict_gc


    SUBROUTINE message_length_A(messagelength, ctyp, igrid)
        ! Subroutine arguments
        INTEGER(int32), INTENT(out) :: messagelength
        CHARACTER(len=1), INTENT(in) :: ctyp
        INTEGER(intk), INTENT(in) :: igrid

        ! Local variables
        INTEGER(intk) :: istart, istop, jstart, jstop, kstart, kstop

        CALL start_and_stop(istart, istop, jstart, jstop, kstart, kstop, ctyp, &
            igrid)
        CALL message_length_impl(messagelength, istart, istop, jstart, jstop, &
            kstart, kstop)
    END SUBROUTINE message_length_A


    SUBROUTINE message_length_B(messagelength, istart, istop, jstart, &
            jstop, kstart, kstop)
        ! Subroutine arguments
        INTEGER(int32), INTENT(out) :: messagelength
        INTEGER(intk), INTENT(in) :: istart, istop, jstart, jstop, kstart, kstop

        CALL message_length_impl(messagelength, istart, istop, jstart, jstop, &
            kstart, kstop)
    END SUBROUTINE message_length_B


    SUBROUTINE message_length_impl(messagelength, istart, istop, jstart, &
            jstop, kstart, kstop)
        ! Function arguments
        INTEGER(int32), INTENT(out) :: messagelength
        INTEGER(intk), INTENT(in) :: istart, istop, jstart, jstop, kstart, kstop

        messagelength = INT(((istop-istart)/2+1) &
            *((jstop-jstart)/2+1) &
            *((kstop-kstart)/2+1), kind=int32)
    END SUBROUTINE message_length_impl


    SUBROUTINE start_and_stop(ista, isto, jsta, jsto, &
            ksta, ksto, ctyp, igrid)
        USE core_mod, ONLY: get_mgdims
        ! Subroutine arguments
        INTEGER(intk), INTENT(out) :: ista, isto, jsta, jsto, ksta, ksto
        CHARACTER(len=1), INTENT(in) :: ctyp
        INTEGER(intk), INTENT(in) :: igrid

        ! Local variables
        INTEGER(intk) :: ii, jj, kk
        INTEGER(intk) :: istart, iend, jstart, jend, kstart, kend

        ! Loop ranges are DO i = istart, ii-iend
        ! returned as DO i = ista, isto, 2
        SELECT CASE (ctyp)
        CASE ("U")
            istart = 2
            iend = 2
            jstart = 3
            jend = 3
            kstart = 3
            kend = 3
        CASE ("V")
            istart = 3
            iend = 3
            jstart = 2
            jend = 2
            kstart = 3
            kend = 3
        CASE ("W")
            istart = 3
            iend = 3
            jstart = 3
            jend = 3
            kstart = 2
            kend = 2
        CASE ("E", "F", "P", "R", "S", 'I', 'T')
            istart = 3
            iend = 3
            jstart = 3
            jend = 3
            kstart = 3
            kend = 3
        CASE ("N")
            istart = 2
            iend = 2
            jstart = 2
            jend = 2
            kstart = 2
            kend = 2
        CASE DEFAULT
            CALL errr(__FILE__, __LINE__)
        END SELECT

        CALL get_mgdims(kk, jj, ii, igrid)

        ista = istart
        isto = ii - iend
        jsta = jstart
        jsto = jj - jend
        ksta = kstart
        ksto = kk - kend

        ! Since the loop has strides of 2, do a sanity check that
        ! the upper index is correct. The stop-index MUST be the actual
        ! index of the last iteration for the message size to be computed
        ! correctly
        !
        ! E.g. The loop:
        !   DO i = 3, kk-3, 2
        ! with kk = 28 will iterate like:
        !   i = 3, 5, 7, 9, 11, 13, 15, 17, 19, 21, 23, 25
        ! with 12 iterations.
        !
        ! If the loop had been:
        !   DO i = 3, kk-2, 2
        ! the iterations would have been the same (the last one being i = 25)
        ! but the message size would be wrong.
        !
        ! Therefore require that MOD(isto-ista, 2) == 0
        IF (MOD(isto-ista, 2) > 0) CALL errr(__FILE__, __LINE__)
        IF (MOD(jsto-jsta, 2) > 0) CALL errr(__FILE__, __LINE__)
        IF (MOD(ksto-ksta, 2) > 0) CALL errr(__FILE__, __LINE__)
    END SUBROUTINE start_and_stop
END MODULE restrict_mod
