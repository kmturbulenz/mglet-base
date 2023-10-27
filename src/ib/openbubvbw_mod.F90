MODULE openbubvbw_mod
    USE core_mod, ONLY: realk, intk, errr, minlevel, maxlevel, nmygridslvl, &
        mygridslvl, get_ip3, get_mgdims, field_t
    USE ibcore_mod, ONLY: openaccur

    IMPLICIT NONE(type, external)
    PRIVATE

    PUBLIC :: openbubvbw

CONTAINS
    SUBROUTINE openbubvbw(au, av, aw, bu, bv, bw)
        ! Subroutine arguments
        TYPE(field_t), INTENT(in) :: au, av, aw
        TYPE(field_t), INTENT(inout) :: bu, bv, bw

        ! Local variables
        INTEGER(intk) :: ilevel

        DO ilevel = minlevel, maxlevel
            CALL openbubvbw_level(ilevel, au, av, aw, bu, bv, bw)
        END DO
    END SUBROUTINE openbubvbw


    SUBROUTINE openbubvbw_level(ilevel, au, av, aw, bu, bv, bw)
        ! Subroutine arguments
        INTEGER(intk), INTENT(in) :: ilevel
        TYPE(field_t), INTENT(in) :: au, av, aw
        TYPE(field_t), INTENT(inout) :: bu, bv, bw

        ! Local variables
        INTEGER(intk) :: i, igrid, kk, jj, ii, ip3

        DO i = 1, nmygridslvl(ilevel)
            igrid = mygridslvl(i, ilevel)

            CALL get_mgdims(kk, jj, ii, igrid)
            CALL get_ip3(ip3, igrid)
            CALL openbubvbw_grid(kk, jj, ii, au%arr(ip3), av%arr(ip3), &
                aw%arr(ip3), bu%arr(ip3), bv%arr(ip3), bw%arr(ip3))
        END DO
    END SUBROUTINE openbubvbw_level


    PURE SUBROUTINE openbubvbw_grid(kk, jj, ii, au, av, aw, bu, bv, bw)
        ! Oeffnet bu,bv,bw, wenn mehr als (1-openaccur) der Flaeche offen ist

        ! Subroutine arguments
        INTEGER(intk), INTENT(in) :: kk, jj, ii
        REAL(realk), INTENT(in) :: au(kk, jj, ii), av(kk, jj, ii), &
            aw(kk, jj, ii)
        REAL(realk), INTENT(inout) :: bu(kk, jj, ii), bv(kk, jj, ii), &
            bw(kk, jj, ii)

        ! Local variables
        INTEGER(intk) :: i, j, k

        DO i = 2, ii-1
            DO j = 2, jj-1
                DO k = 2, kk-1
                    IF (au(k, j, i) >= (1.0-openaccur) &
                            .AND. NINT(bu(k, j, i)) == 0) THEN
                        bu(k,j,i) = 1.0
                    END IF
                    IF (av(k, j, i) >= (1.0-openaccur) &
                            .AND. NINT(bv(k, j, i)) == 0) THEN
                        bv(k,j,i) = 1.0
                    END IF
                    IF (aw(k, j, i) >= (1.0-openaccur) &
                            .AND. NINT(bw(k, j, i)) == 0) THEN
                        bw(k,j,i) = 1.0
                    END IF
                END DO
            END DO
        END DO
    END SUBROUTINE openbubvbw_grid
END MODULE openbubvbw_mod
