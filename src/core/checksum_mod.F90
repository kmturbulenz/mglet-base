MODULE checksum_mod
    USE precision_mod, ONLY: intk, realk, int_bytes, real_bytes
    USE grids_mod, ONLY: mygridslvl, nmygridslvl, get_mgdims, &
        minlevel, maxlevel
    USE pointers_mod, ONLY: get_ip3, idim3d
    USE, INTRINSIC :: ISO_C_BINDING, ONLY: C_PTR, C_LOC, C_LONG

    IMPLICIT NONE(type, external)
    PRIVATE

    INTERFACE
        ! Interface to zlib crc32
        ! unsigned long crc32(unsigned long crc, const unsigned char * buf,
        !                     unsigned int len);
        !
        ! It sucks that Fortran does not have unsigned datatypes, I just
        ! have to declare unsigned long as long...
        INTEGER(c_long) FUNCTION crc32(crc, buf, len) BIND(C)
            USE, INTRINSIC :: ISO_C_BINDING, ONLY: C_LONG, C_PTR
            INTEGER(c_long), INTENT(in), VALUE :: crc
            TYPE(C_PTR), INTENT(in), VALUE :: buf
            INTEGER(c_long), INTENT(in), VALUE :: len
        END FUNCTION crc32
    END INTERFACE

    INTERFACE checksum
        MODULE PROCEDURE :: checksum_r1, checksum_r3, checksum_i1, checksum_i3
    END INTERFACE checksum

    INTERFACE checksum_field
        MODULE PROCEDURE :: checksum_rfield, checksum_ifield
    END INTERFACE checksum_field

    INTERFACE checksum_level
        MODULE PROCEDURE :: checksum_rlevel, checksum_ilevel
    END INTERFACE checksum_level

    INTERFACE checksum_grid
        MODULE PROCEDURE :: checksum_rgrid, checksum_igrid
    END INTERFACE checksum_grid

    PUBLIC :: checksum, checksum_field, checksum_level, checksum_grid

CONTAINS
    INTEGER(c_long) FUNCTION checksum_r1(field)
        ! Function arguments
        REAL(realk), INTENT(in), TARGET :: field(:)

        ! Local variables
        INTEGER(c_long) :: len, crc
        TYPE(C_PTR) :: ptr

        len = real_bytes*SIZE(field)
        ptr = C_LOC(field)
        crc = 0

        checksum_r1 = crc32(crc, ptr, len)
    END FUNCTION checksum_r1


    INTEGER(c_long) FUNCTION checksum_r3(field)
        ! Function arguments
        REAL(realk), INTENT(in), TARGET :: field(:, :, :)

        ! Local variables
        INTEGER(c_long) :: len, crc
        TYPE(C_PTR) :: ptr

        len = real_bytes*SIZE(field)
        ptr = C_LOC(field)
        crc = 0

        checksum_r3 = crc32(crc, ptr, len)
    END FUNCTION checksum_r3


    INTEGER(c_long) FUNCTION checksum_i1(field)
        ! Function arguments
        INTEGER(intk), INTENT(in), TARGET :: field(:)

        ! Local variables
        INTEGER(c_long) :: len, crc
        TYPE(C_PTR) :: ptr

        len = int_bytes*SIZE(field)
        ptr = C_LOC(field)
        crc = 0

        checksum_i1 = crc32(crc, ptr, len)
    END FUNCTION checksum_i1


    INTEGER(c_long) FUNCTION checksum_i3(field)
        ! Function arguments
        INTEGER(intk), INTENT(in), TARGET :: field(:, :, :)

        ! Local variables
        INTEGER(c_long) :: len, crc
        TYPE(C_PTR) :: ptr

        len = int_bytes*SIZE(field)
        ptr = C_LOC(field)
        crc = 0

        checksum_i3 = crc32(crc, ptr, len)
    END FUNCTION checksum_i3


    SUBROUTINE checksum_rfield(msg, field)
        ! Subroutine arguments
        CHARACTER(len=*), INTENT(in) :: msg
        REAL(realk), INTENT(in) :: field(idim3d)

        ! Local variables
        INTEGER(intk) :: ilevel

        DO ilevel = minlevel, maxlevel
            CALL checksum_level(ilevel, msg, field)
        END DO
    END SUBROUTINE checksum_rfield


    SUBROUTINE checksum_rlevel(ilevel, msg, field)
        ! Subroutine arguments
        INTEGER(intk), INTENT(in) :: ilevel
        CHARACTER(len=*), INTENT(in) :: msg
        REAL(realk), INTENT(in) :: field(idim3d)

        ! Local variables
        INTEGER(intk) :: i, igrid, kk, jj, ii, ip3

        DO i = 1, nmygridslvl(ilevel)
            igrid = mygridslvl(i, ilevel)

            CALL get_mgdims(kk, jj, ii, igrid)
            CALL get_ip3(ip3, igrid)
            CALL checksum_rgrid(kk, jj, ii, igrid, msg, field(ip3))
        END DO
    END SUBROUTINE checksum_rlevel


    SUBROUTINE checksum_rgrid(kk, jj, ii, igrid, msg, field)
        ! Subroutine arguments
        INTEGER(intk), INTENT(in) :: kk, jj, ii, igrid
        CHARACTER(len=*), INTENT(in) :: msg
        REAL(realk), INTENT(in) :: field(kk, jj, ii)

        WRITE(*, *) msg, " igrid, checksum: ", igrid, checksum(field)
    END SUBROUTINE checksum_rgrid


    SUBROUTINE checksum_ifield(msg, field)
        ! Subroutine arguments
        CHARACTER(len=*), INTENT(in) :: msg
        INTEGER(intk), INTENT(in) :: field(idim3d)

        ! Local variables
        INTEGER(intk) :: ilevel

        DO ilevel = minlevel, maxlevel
            CALL checksum_level(ilevel, msg, field)
        END DO
    END SUBROUTINE checksum_ifield


    SUBROUTINE checksum_ilevel(ilevel, msg, field)
        ! Subroutine arguments
        INTEGER(intk), INTENT(in) :: ilevel
        CHARACTER(len=*), INTENT(in) :: msg
        INTEGER(intk), INTENT(in) :: field(idim3d)

        ! Local variables
        INTEGER(intk) :: i, igrid, kk, jj, ii, ip3

        DO i = 1, nmygridslvl(ilevel)
            igrid = mygridslvl(i, ilevel)

            CALL get_mgdims(kk, jj, ii, igrid)
            CALL get_ip3(ip3, igrid)
            CALL checksum_igrid(kk, jj, ii, igrid, msg, field(ip3))
        END DO
    END SUBROUTINE checksum_ilevel


    SUBROUTINE checksum_igrid(kk, jj, ii, igrid, msg, field)
        ! Subroutine arguments
        INTEGER(intk), INTENT(in) :: kk, jj, ii, igrid
        CHARACTER(len=*), INTENT(in) :: msg
        INTEGER(intk), INTENT(in) :: field(kk, jj, ii)

        WRITE(*, *) msg, " igrid, checksum: ", igrid, checksum(field)
    END SUBROUTINE checksum_igrid
END MODULE checksum_mod
