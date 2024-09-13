MODULE scastat_mod
    USE core_mod
    USE scacore_mod

    IMPLICIT NONE (type, external)
    PRIVATE

    PUBLIC :: init_scastat, finish_scastat

CONTAINS

    SUBROUTINE init_scastat()
        ! Subroutine arguments
        ! none...

        ! Local variables
        INTEGER(intk) :: l

        DO l = 1, nsca
            ! checking length of scalar name string
            IF (LEN_TRIM(scalar(l)%name) + 4 > nchar_name) THEN
                CALL errr(__FILE__, __LINE__)
            END IF

            ! T_AVG
            CALL register_statfield(TRIM(scalar(l)%name)&
                //"_AVG", comp_avg)

            ! TT_AVG
            CALL register_statfield(TRIM(scalar(l)%name)&
                //TRIM(scalar(l)%name)&
                //"_AVG", comp_sqr_avg)

            ! TxTx_AVG
            CALL register_statfield(TRIM(scalar(l)%name)&
                //"x"//TRIM(scalar(l)%name)&
                //"x"//"_AVG", comp_txtx_avg)
            ! TyTy_AVG
            CALL register_statfield(TRIM(scalar(l)%name)&
                //"y"//TRIM(scalar(l)%name)&
                //"y"//"_AVG", comp_txtx_avg)
            ! TzTz_AVG
            CALL register_statfield(TRIM(scalar(l)%name)&
                //"z"//TRIM(scalar(l)%name)&
                //"z"//"_AVG", comp_txtx_avg)

            ! UT_AVG
            CALL register_statfield("U"//TRIM(scalar(l)%name)//&
                "_AVG", comp_ut_avg)
            ! VT_AVG
            CALL register_statfield("V"//TRIM(scalar(l)%name)//&
                "_AVG", comp_ut_avg)
            ! WT_AVG
            CALL register_statfield("W"//TRIM(scalar(l)%name)//&
                "_AVG", comp_ut_avg)

            !UTT_AVG
            CALL register_statfield("U"//TRIM(scalar(l)%name)&
            //TRIM(scalar(l)%name)//"_AVG", comp_utt_avg)
            !VTT_AVG
            CALL register_statfield("V"//TRIM(scalar(l)%name)&
            //TRIM(scalar(l)%name)//"_AVG", comp_utt_avg)
            !WTT_AVG
            CALL register_statfield("W"//TRIM(scalar(l)%name)&
            //TRIM(scalar(l)%name)//"_AVG", comp_utt_avg)
        END DO

    END SUBROUTINE init_scastat


    SUBROUTINE finish_scastat
        CONTINUE
    END SUBROUTINE finish_scastat


    ! Routine to compute the TxTx_AVG, TyTy_AVG and TzTz_AVG fields
    ! Computation of spatial derivatives with central difference
    !
    SUBROUTINE  comp_txtx_avg(field, name, dt)
        ! Subroutine arguments
        TYPE(field_t), INTENT(inout) :: field
        CHARACTER(len=*), INTENT(in) :: name
        REAL(realk), INTENT(in) :: dt

        !Local Variables
        TYPE(field_t) :: tx_f ! derivative of a scalar in any direction
        ! The units of the scalar are asummed to be °K (Temperature)
        INTEGER(intk), PARAMETER :: units(*) = [0, -2, 0, 2, 0, 0, 0]
        INTEGER(intk), PARAMETER :: units_tx(*) = [0, -1, 0, 1, 0, 0, 0]
        TYPE(field_t), POINTER :: t_f
        INTEGER(intk) :: istag, jstag, kstag
        CHARACTER(len=3) :: ivar
        CHARACTER(len=nchar_name) :: name_tx
        INTEGER(intk) :: nchar
        INTEGER(intk) :: sca_name_length
        CHARACTER(len=nchar_name) :: sca_name

        nchar = LEN_TRIM(name)
        sca_name_length = (nchar-6)/2
        sca_name = name(1:sca_name_length)
        name_tx = name(1:sca_name_length+1)

        CALL get_field(t_f, sca_name)

        istag = 0
        jstag = 0
        kstag = 0

        SELECT CASE (TRIM(name(nchar-4:nchar)))
        CASE ("x_AVG")
            ivar = "DXT"
        CASE ("y_AVG")
            ivar = "DYT"
        CASE ("z_AVG")
            ivar = "DZT"
        CASE DEFAULT
            CALL errr(__FILE__, __LINE__)
        END SELECT

        CALL field%init(name, istag=istag, jstag=jstag, kstag=kstag, &
            units=units)
        CALL tx_f%init('tmp', istag=istag, jstag=jstag, kstag=kstag, &
            units=units_tx)

        ! central difference on scalar (= no staggering)
        CALL differentiate(tx_f, t_f, ivar)

        field%arr = tx_f%arr(:)**2
        CALL tx_f%finish()

    END SUBROUTINE


    ! Routine to compute the UT_AVG, VT_AVG, and WT_AVG fields
    ! Computation at staggered positions of velocity
    !
    SUBROUTINE comp_ut_avg(field, name, dt)
        !Subroutine arguments
        TYPE(field_t), INTENT(inout) :: field
        CHARACTER(len=*), INTENT(in) :: name
        REAL(realk), INTENT(in) :: dt

        ! Local variables
        INTEGER(intk), PARAMETER :: units(*) = [0, 1, -1, 1, 0, 0, 0]
        TYPE(field_t), POINTER :: u_f, t_f
        INTEGER(intk) :: istag, jstag, kstag
        INTEGER(intk) :: sca_name_length
        INTEGER(intk) :: nchar
        CHARACTER(len=nchar_name) :: sca_name

        nchar = LEN_TRIM(name)
        sca_name_length = nchar-5
        sca_name = name(2:sca_name_length+1)

        SELECT CASE (name(1:1))
        CASE ("U")
            CALL get_field(u_f, "U")
            istag = 1
            jstag = 0
            kstag = 0
        CASE ("V")
            CALL get_field(u_f, "V")
            istag = 0
            jstag = 1
            kstag = 0
        CASE ("W")
            CALL get_field(u_f, "W")
            istag = 0
            jstag = 0
            kstag = 1
        CASE DEFAULT
            CALL errr(__FILE__, __LINE__)
        END SELECT

        CALL get_field(t_f, sca_name)

        CALL field%init(name, istag=istag, jstag=jstag, kstag=kstag, &
            units=units)

        ! multiplication at staggered positions
        CALL field%multiply(u_f, t_f)

    END SUBROUTINE comp_ut_avg


    ! Routine to compute the UTT_AVG, VTT_AVG, and WTT_AVG fields
    ! Computation at staggered positions of velocity
    !
    SUBROUTINE comp_utt_avg(field, name, dt)
       ! Subroutine arguments
        TYPE(field_t), INTENT(inout) :: field
        CHARACTER(len=*), INTENT(in) :: name
        REAL(realk), INTENT(in) :: dt

        !Local variables
        INTEGER(intk), PARAMETER :: units(*) = [0, 1, -1, 2, 0, 0, 0]
        TYPE(field_t), POINTER :: u_f, t_f
        INTEGER(intk) :: istag, jstag, kstag
        INTEGER(intk) :: sca_name_length
        CHARACTER(len=nchar_name) :: sca_name
        INTEGER(intk) :: nchar

        nchar = LEN_TRIM(name)
        IF ( MODULO((nchar-5),2) /= 0 ) THEN
            CALL errr(__FILE__, __LINE__)
        END IF
        sca_name_length = (nchar-5) / 2
        sca_name = name(2:sca_name_length+1)

        SELECT CASE (name(1:1))
        CASE ("U")
            CALL get_field(u_f, "U")
            istag = 1
            jstag = 0
            kstag = 0
        CASE ("V")
            CALL get_field(u_f, "V")
            istag = 0
            jstag = 1
            kstag = 0
        CASE ("W")
            CALL get_field(u_f, "W")
            istag = 0
            jstag = 0
            kstag = 1
        CASE DEFAULT
            CALL errr(__FILE__, __LINE__)
        END SELECT

        CALL get_field(t_f, sca_name)

        CALL field%init(name, istag=istag, jstag=jstag, kstag=kstag, &
            units=units)

        ! multiplication at staggered positions
        CALL field%multiply(u_f, t_f, t_f)

    END SUBROUTINE

END MODULE scastat_mod
