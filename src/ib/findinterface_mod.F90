MODULE findinterface_mod
    USE core_mod, ONLY: realk, intk, errr

    IMPLICIT NONE(type, external)
    PRIVATE

    INTERFACE findinterface2
        MODULE PROCEDURE :: findinterface2_a, findinterface2_b
    END INTERFACE findinterface2

    PUBLIC :: findinterface2, findinterface3

CONTAINS
    PURE SUBROUTINE findinterface2_a(k, j, i, kk, jj, ii, blocked, found, foundx1, &
            foundx2, foundy1, foundy2, foundz1, foundz2, foundnr)
        ! Subroutine arguments
        INTEGER(intk), INTENT(in) :: k, j, i
        INTEGER(intk), INTENT(in) :: kk, jj, ii
        REAL(realk), INTENT(in) :: blocked(kk, jj, ii)
        INTEGER(intk), INTENT(out) :: found, foundx1, foundx2, &
            foundy1, foundy2, foundz1, foundz2, foundnr

        found = 0
        foundx1 = 0
        foundx2 = 0
        foundy1 = 0
        foundy2 = 0
        foundz1 = 0
        foundz2 = 0
        foundnr = 0

        IF (blocked(k, j ,i) <= 0.5) THEN
            IF (blocked(k, j, i+1) >= 0.5) THEN
                foundnr = foundnr + 1
                found = 1
                foundx1 = 1
            END IF
            IF (blocked(k, j, i-1) >= 0.5) THEN
                foundnr = foundnr + 1
                found = 1
                foundx2 = 1
            END IF
            IF (blocked(k, j+1, i) >= 0.5) THEN
                foundnr = foundnr + 1
                found = 1
                foundy1 = 1
            END IF
            IF (blocked(k, j-1, i) >= 0.5) THEN
                foundnr = foundnr + 1
                found = 1
                foundy2 = 1
            END IF
            IF (blocked(k+1, j, i) >= 0.5) THEN
                foundnr = foundnr + 1
                found = 1
                foundz1 = 1
            END IF
            IF (blocked(k-1, j, i) >= 0.5) THEN
                foundnr = foundnr + 1
                found = 1
                foundz2 = 1
            END IF
        END IF
    END SUBROUTINE findinterface2_a


    PURE SUBROUTINE findinterface2_b(k, j, i, kk, jj, ii, blocked, found)
        ! Subroutine arguments
        INTEGER(intk), INTENT(in) :: k, j, i
        INTEGER(intk), INTENT(in) :: kk, jj, ii
        REAL(realk), INTENT(in) :: blocked(kk, jj, ii)
        INTEGER(intk), INTENT(out) :: found

        ! Local variables
        INTEGER(intk) :: foundx1, foundx2, foundy1, foundy2, foundz1, &
            foundz2, foundnr

        CALL findinterface2_a(k, j, i, kk, jj, ii, blocked, found, foundx1, &
            foundx2, foundy1, foundy2, foundz1, foundz2, foundnr)
    END SUBROUTINE findinterface2_b


    SUBROUTINE findinterface3(k, j, i, kk, jj, ii, compon, &
            bp, bu, bv, bw, found)

        ! Subroutine arguments
        INTEGER(intk), INTENT(in) :: k, j, i
        INTEGER(intk), INTENT(in) :: kk, jj, ii
        INTEGER(intk), INTENT(in) :: compon
        REAL(realk), INTENT(in) :: bp(kk, jj, ii), bu(kk, jj, ii), &
            bv(kk, jj, ii), bw(kk, jj, ii)
        INTEGER(intk), INTENT(out) :: found

        found = 0

        SELECT CASE(compon)
        CASE (1)
            IF (NINT(bu(k,j,i)) == 0) THEN
                IF (NINT(bp(k, j, i)) == 1) found = 1
                IF (NINT(bp(k, j, i+1)) == 1) found = 1
            END IF
        CASE (2)
            IF (NINT(bv(k, j, i)) == 0) then
                IF (NINT(bp(k, j, i )) == 1) found = 1
                IF (NINT(bp(k, j+1, i)) == 1) found = 1
            END IF
        CASE (3)
            IF (NINT(bw(k, j, i)) == 0) THEN
                IF (NINT(bp(k, j, i)) == 1) found = 1
                IF (NINT(bp(k+1, j, i)) == 1) found = 1
            END IF
        CASE DEFAULT
            CALL errr(__FILE__, __LINE__)
        END SELECT
    END SUBROUTINE findinterface3

END MODULE findinterface_mod
