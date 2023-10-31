! Module for wrapping around various STL reading functions
MODULE readstl_mod
    USE precision_mod, ONLY: intk, realk, int32, real32
    USE err_mod, ONLY: errr

    IMPLICIT NONE(type, external)
    PRIVATE

    INTERFACE
        FUNCTION stl_is_binary(filename) RESULT(ntri) BIND(C)
            USE, INTRINSIC :: ISO_C_BINDING, ONLY: C_CHAR, C_INT
            CHARACTER(C_CHAR), DIMENSION(*), INTENT(IN) :: filename
            INTEGER(C_INT) :: ntri
        END FUNCTION stl_is_binary
    END INTERFACE

    INTERFACE stl_read
        PROCEDURE stl_read_1d, stl_read_3d
    END INTERFACE

    ! Public data items
    PUBLIC :: stl_info, stl_read

CONTAINS
    SUBROUTINE stl_info(filename, binary, ntri)
        USE, INTRINSIC :: ISO_C_BINDING, ONLY: C_NULL_CHAR

        ! Subroutine arguments
        CHARACTER(len=*), INTENT(in) :: filename
        LOGICAL, INTENT(out) :: binary
        INTEGER(intk), INTENT(out) :: ntri

        binary = .TRUE.
        ntri = stl_is_binary(TRIM(filename) // C_NULL_CHAR)

        ! The stl_is_binary return either the number of trinagles in case the
        ! file is binary, or -1 if the file is ASCII
        IF (ntri < 0) THEN
            binary = .FALSE.
            CALL read_ascii(filename, ntri)
        END IF

        IF (binary) THEN
            WRITE(*, "(' STL: ', A, ' is binary, containing ', I0, " &
                //"' triangles')") TRIM(filename), ntri
        ELSE
            WRITE(*, "(' STL: ', A, ' is ASCII, containing ', I0, " &
                //"' triangles')") TRIM(filename), ntri
        END IF
    END SUBROUTINE stl_info


    SUBROUTINE stl_read_1d(filename, binary, vertices)
        ! Subroutine arguments
        CHARACTER(len=*), INTENT(in) :: filename
        LOGICAL, INTENT(in) :: binary
        REAL(realk), INTENT(inout), TARGET :: vertices(:)

        ! Local variables
        INTEGER(intk) :: ntri
        REAL(realk), POINTER :: vertices3d(:, :, :)

        ntri = SIZE(vertices)/9
        IF (ntri*9 /= SIZE(vertices)) CALL errr(__FILE__, __LINE__)

        vertices3d(1:3, 1:3, 1:ntri) => vertices(1:3*3*ntri)
        CALL stl_read_3d(filename, binary, vertices3d)
        NULLIFY(vertices3d)
    END SUBROUTINE stl_read_1d


    SUBROUTINE stl_read_3d(filename, binary, vertices)
        ! Subroutine arguments
        CHARACTER(len=*), INTENT(in) :: filename
        LOGICAL, INTENT(in) :: binary
        REAL(realk), INTENT(inout) :: vertices(:, :, :)

        ! Local variables
        INTEGER(intk) :: ntri

        IF (binary) THEN
            CALL read_binary(filename, ntri, vertices)
        ELSE
            CALL read_ascii(filename, ntri, vertices)
        END IF
    END SUBROUTINE stl_read_3d


    SUBROUTINE read_ascii(filename, itri, vertices)
        ! Subroutine arguments
        CHARACTER(len=*), INTENT(in) :: filename
        INTEGER(intk), INTENT(out) :: itri
        REAL(realk), INTENT(inout), OPTIONAL :: vertices(:, :, :)

        ! Local variables
        INTEGER(intk) :: fu, ntri, ios, idir
        CHARACTER(len=64) :: cdummy

        ! Sanity check
        IF (PRESENT(vertices)) THEN
            IF (SIZE(vertices, 1) /= 3) CALL errr(__FILE__, __LINE__)
            IF (SIZE(vertices, 2) /= 3) CALL errr(__FILE__, __LINE__)
            ntri = SIZE(vertices, 3)
        END IF

        OPEN(newunit=fu, file=filename, action="READ")
        itri = 0
        DO
            ! reading line into dummy variable
            READ(fu, *, IOSTAT=ios) cdummy

            ! Exit loop if no further lines left to read
            IF (ios /= 0) EXIT

            IF (cdummy(1:5) == 'outer') THEN
                itri = itri + 1

                IF (PRESENT(vertices)) THEN
                    IF (itri > ntri) THEN
                        CALL errr(__FILE__, __LINE__)
                    END IF

                    READ(fu, *) cdummy, (vertices(idir, 1, itri), idir=1,3)
                    READ(fu, *) cdummy, (vertices(idir, 2, itri), idir=1,3)
                    READ(fu, *) cdummy, (vertices(idir, 3, itri), idir=1,3)
                ELSE
                    READ(fu, *)
                    READ(fu, *)
                    READ(fu, *)
                END IF
            END IF
        END DO
        CLOSE(fu)

        IF (PRESENT(vertices)) THEN
            IF (itri /= ntri) THEN
                CALL errr(__FILE__, __LINE__)
            END IF
        END IF
    END SUBROUTINE read_ascii


    SUBROUTINE read_binary(filename, itri, vertices)
        USE, INTRINSIC :: ISO_FORTRAN_ENV, ONLY: int16

        ! Subroutine arguments
        CHARACTER(len=*), INTENT(in) :: filename
        INTEGER(intk), INTENT(out) :: itri
        REAL(realk), INTENT(inout) :: vertices(:, :, :)

        ! Local variables
        INTEGER(intk) :: fu, ntri, ios
        INTEGER(int32) :: itri32
        CHARACTER(len=80) :: cline
        REAL(real32) :: tridata(12)
        INTEGER(int16) :: attr

        ! Sanity check
        IF (SIZE(vertices, 1) /= 3) CALL errr(__FILE__, __LINE__)
        IF (SIZE(vertices, 2) /= 3) CALL errr(__FILE__, __LINE__)
        ntri = SIZE(vertices, 3)

        OPEN(newunit=fu, file=filename, form='UNFORMATTED', access='STREAM', &
            action="READ")
        READ(fu) cline
        READ(fu) itri32
        IF (itri32 /= ntri) THEN
            CALL errr(__FILE__, __LINE__)
        END IF

        DO itri = 1, ntri
            READ(fu, IOSTAT=ios) tridata, attr

            ! In binary mode we do not read beyond the end of file deliberately,
            ! so an non-successful read is a failure.
            IF (ios /= 0) THEN
                CALL errr(__FILE__, __LINE__)
            END IF

            vertices(1, 1, itri) = tridata(4)
            vertices(2, 1, itri) = tridata(5)
            vertices(3, 1, itri) = tridata(6)
            vertices(1, 2, itri) = tridata(7)
            vertices(2, 2, itri) = tridata(8)
            vertices(3, 2, itri) = tridata(9)
            vertices(1, 3, itri) = tridata(10)
            vertices(2, 3, itri) = tridata(11)
            vertices(3, 3, itri) = tridata(12)
        END DO
        CLOSE(fu)
    END SUBROUTINE read_binary
END MODULE readstl_mod
