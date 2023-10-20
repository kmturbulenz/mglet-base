MODULE scastat_mod
    USE core_mod
    USE scacore_mod

    PUBLIC :: init_scastat, finish_scastat

CONTAINS
    SUBROUTINE init_scastat()
        ! Subroutine arguments
        ! none...

        ! Local variables
        INTEGER(intk) :: l

        DO l = 1, nsca
            IF (LEN_TRIM(scalar(l)%name) + 4 > nchar_name) THEN
                CALL errr(__FILE__, __LINE__)
            END IF
            CALL register_statfield(TRIM(scalar(l)%name)//"_AVG", comp_avg)
        END DO
    END SUBROUTINE init_scastat


    SUBROUTINE finish_scastat
        CONTINUE
    END SUBROUTINE finish_scastat
END MODULE scastat_mod
