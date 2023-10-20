MODULE setboundarybuffers_scalar_mod
    USE core_mod
    USE scacore_mod

    IMPLICIT NONE(type, external)
    PRIVATE

    ! Bound operation 'T' operate on U, V, W, P
    TYPE, EXTENDS(bound_t) :: setboundarybuffers_scalar_t
    CONTAINS
        PROCEDURE, NOPASS :: front => bfront
        PROCEDURE, NOPASS :: back => bfront
        PROCEDURE, NOPASS :: right => bright
        PROCEDURE, NOPASS :: left => bright
        PROCEDURE, NOPASS :: bottom => bbottom
        PROCEDURE, NOPASS :: top => bbottom
    END TYPE setboundarybuffers_scalar_t
    TYPE(setboundarybuffers_scalar_t) :: setboundarybuffers_scalar

    PUBLIC :: setboundarybuffers_scalar

CONTAINS
    SUBROUTINE bfront(igrid, iface, ibocd, ctyp, f1, f2, f3, f4, timeph)
        ! Subroutine arguments
        INTEGER(intk), INTENT(in) :: igrid, iface, ibocd
        CHARACTER(len=*), INTENT(in) :: ctyp
        TYPE(field_t), INTENT(inout) :: f1
        TYPE(field_t), INTENT(inout), OPTIONAL :: f2, f3, f4
        REAL(realk), INTENT(in), OPTIONAL :: timeph

        ! Local variables
        INTEGER(intk) :: idx
        REAL(realk), POINTER, CONTIGUOUS :: tbuf(:, :, :)
        REAL(realk) :: scbvalue(nsca)

        ! Only works on scalar SWA and SIO boundaries, should do nothing
        ! otherwise
        SELECT CASE (ctyp)
        CASE ("SIO", "SWA")
            CONTINUE
        CASE DEFAULT
            RETURN
        END SELECT

        ! Assure that required fields are present
        IF (PRESENT(f2) .OR. PRESENT(f3) .OR. PRESENT(f4)) THEN
            CALL errr(__FILE__, __LINE__)
        END IF

        ! Get value of buffer (we do not need the type here...)
        CALL f1%get_attr(idx, "SCAIDX")
        CALL get_bcprms(scbvalue, igrid, iface, ibocd)

        ! Apply value to buffer
        CALL f1%buffers%get_buffer(tbuf, igrid, iface)
        tbuf = scbvalue(idx)
    END SUBROUTINE bfront


    SUBROUTINE bright(igrid, iface, ibocd, ctyp, f1, f2, f3, f4, timeph)
        ! Subroutine arguments
        INTEGER(intk), INTENT(in) :: igrid, iface, ibocd
        CHARACTER(len=*), INTENT(in) :: ctyp
        TYPE(field_t), INTENT(inout) :: f1
        TYPE(field_t), INTENT(inout), OPTIONAL :: f2, f3, f4
        REAL(realk), INTENT(in), OPTIONAL :: timeph

        ! Local variables
        INTEGER(intk) :: idx
        REAL(realk), POINTER, CONTIGUOUS :: tbuf(:, :, :)
        REAL(realk) :: scbvalue(nsca)

        ! Only works on scalar SWA and SIO boundaries, should do nothing
        ! otherwise
        SELECT CASE (ctyp)
        CASE ("SIO", "SWA")
            CONTINUE
        CASE DEFAULT
            RETURN
        END SELECT

        ! Assure that required fields are present
        IF (PRESENT(f2) .OR. PRESENT(f3) .OR. PRESENT(f4)) THEN
            CALL errr(__FILE__, __LINE__)
        END IF

        ! Get value of buffer (we do not need the type here...)
        CALL f1%get_attr(idx, "SCAIDX")
        CALL get_bcprms(scbvalue, igrid, iface, ibocd)

        ! Apply value to buffer
        CALL f1%buffers%get_buffer(tbuf, igrid, iface)
        tbuf = scbvalue(idx)
    END SUBROUTINE bright


    SUBROUTINE bbottom(igrid, iface, ibocd, ctyp, f1, f2, f3, f4, timeph)
        ! Subroutine arguments
        INTEGER(intk), INTENT(in) :: igrid, iface, ibocd
        CHARACTER(len=*), INTENT(in) :: ctyp
        TYPE(field_t), INTENT(inout) :: f1
        TYPE(field_t), INTENT(inout), OPTIONAL :: f2, f3, f4
        REAL(realk), INTENT(in), OPTIONAL :: timeph

        ! Local variables
        INTEGER(intk) :: idx
        REAL(realk), POINTER, CONTIGUOUS :: tbuf(:, :, :)
        REAL(realk) :: scbvalue(nsca)

        ! Only works on scalar SWA and SIO boundaries, should do nothing
        ! otherwise
        SELECT CASE (ctyp)
        CASE ("SIO", "SWA")
            CONTINUE
        CASE DEFAULT
            RETURN
        END SELECT

        ! Assure that required fields are present
        IF (PRESENT(f2) .OR. PRESENT(f3) .OR. PRESENT(f4)) THEN
            CALL errr(__FILE__, __LINE__)
        END IF

        ! Get value of buffer (we do not need the type here...)
        CALL f1%get_attr(idx, "SCAIDX")
        CALL get_bcprms(scbvalue, igrid, iface, ibocd)

        ! Apply value to buffer
        CALL f1%buffers%get_buffer(tbuf, igrid, iface)
        tbuf = scbvalue(idx)
    END SUBROUTINE bbottom
END MODULE setboundarybuffers_scalar_mod
