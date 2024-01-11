MODULE flowstat_mod
    USE core_mod

    IMPLICIT NONE(type, external)
    PRIVATE

    PUBLIC :: init_flowstat, finish_flowstat

CONTAINS
    SUBROUTINE init_flowstat()
        ! Subroutine arguments
        ! none...

        ! Local variables
        ! none....

        CALL register_statfield("U_AVG", comp_avg)
        CALL register_statfield("V_AVG", comp_avg)
        CALL register_statfield("W_AVG", comp_avg)
        CALL register_statfield("P_AVG", comp_avg)
        CALL register_statfield("G_AVG", comp_avg)

        CALL register_statfield("UU_AVG", comp_sqr_avg)
        CALL register_statfield("VV_AVG", comp_sqr_avg)
        CALL register_statfield("WW_AVG", comp_sqr_avg)
        CALL register_statfield("PP_AVG", comp_sqr_avg)
        CALL register_statfield("GG_AVG", comp_sqr_avg)

        CALL register_statfield("UV_AVG", comp_uv_avg)
        CALL register_statfield("UW_AVG", comp_uv_avg)
        CALL register_statfield("VW_AVG", comp_uv_avg)

        CALL register_statfield("laplaceP_AVG", comp_laplacep_avg)
        CALL register_statfield("laplaceP_SQR_AVG", comp_laplacep_sqr_avg)
    END SUBROUTINE init_flowstat


    SUBROUTINE finish_flowstat
        CONTINUE
    END SUBROUTINE finish_flowstat


    ! Routine to compute the UV_AVG, UW_AVG and VW_AVG fields
    SUBROUTINE comp_uv_avg(field, name, dt)
        ! Subroutine arguments
        TYPE(field_t), INTENT(inout) :: field
        CHARACTER(len=*), INTENT(in) :: name
        REAL(realk), INTENT(in) :: dt

        ! Local variables
        INTEGER(intk), PARAMETER :: units(*) = [0, 2, -2, 0, 0, 0, 0]
        TYPE(field_t), POINTER :: in1, in2
        INTEGER(intk) :: istag, jstag, kstag

        SELECT CASE (TRIM(name))
        CASE ("UV_AVG")
            CALL get_field(in1, "U")
            CALL get_field(in2, "V")
            istag = 1
            jstag = 1
            kstag = 0
        CASE ("UW_AVG")
            CALL get_field(in1, "U")
            CALL get_field(in2, "W")
            istag = 1
            jstag = 0
            kstag = 1
        CASE ("VW_AVG")
            CALL get_field(in1, "V")
            CALL get_field(in2, "W")
            istag = 0
            jstag = 1
            kstag = 1
        CASE DEFAULT
            CALL errr(__FILE__, __LINE__)
        END SELECT

        CALL field%init(name, istag=istag, jstag=jstag, kstag=kstag, &
            units=units)
        CALL field%multiply(in1, in2)
    END SUBROUTINE comp_uv_avg


    SUBROUTINE comp_laplacep_avg(field, name, dt)
        ! Subroutine arguments
        TYPE(field_t), INTENT(inout) :: field
        CHARACTER(len=*), INTENT(in) :: name
        REAL(realk), INTENT(in) :: dt

        ! Local variables
        TYPE(field_t), POINTER :: p_f, bp_f, dx_f, dy_f, dz_f
        REAL(realk), CONTIGUOUS, POINTER :: p(:, :, :), bp(:, :, :), &
            lpp(:, :, :), dx(:), dy(:), dz(:)
        INTEGER(intk), PARAMETER :: units(*) = [1, -3, -2, 0, 0, 0, 0]
        INTEGER(intk) :: i, igrid
        INTEGER(intk) :: nfro, nbac, nlft, nrgt, ntop, nbot
        INTEGER(intk) :: kk, jj, ii

        IF (name /= "laplaceP_AVG") CALL errr(__FILE__, __LINE__)

        CALL get_field(p_f, "P")
        CALL get_field(bp_f, "BP")
        CALL get_field(dx_f, "DX")
        CALL get_field(dy_f, "DY")
        CALL get_field(dz_f, "DZ")

        ! Create a field to store the result, this creates an empty field
        CALL field%init(name, units=units)

        ! Compute laplaceP
        DO i = 1, nmygrids
            igrid = mygrids(i)

            CALL field%get_ptr(lpp, igrid)
            CALL p_f%get_ptr(p, igrid)
            CALL bp_f%get_ptr(bp, igrid)
            CALL dx_f%get_ptr(dx, igrid)
            CALL dy_f%get_ptr(dy, igrid)
            CALL dz_f%get_ptr(dz, igrid)
            CALL get_mgbasb(nfro, nbac, nrgt, nlft, nbot, ntop, igrid)
            CALL get_mgdims(kk, jj, ii, igrid)

            CALL calclpp_grid(kk, jj, ii, lpp, p, bp, dx, dy, dz, &
                nfro, nbac, nrgt, nlft, nbot, ntop)
        END DO
    END SUBROUTINE comp_laplacep_avg


    SUBROUTINE comp_laplacep_sqr_avg(field, name, dt)
        ! Subroutine arguments
        TYPE(field_t), INTENT(inout) :: field
        CHARACTER(len=*), INTENT(in) :: name
        REAL(realk), INTENT(in) :: dt

        ! Local variables
        ! nobe...

        IF (name /= "laplaceP_SQR_AVG") CALL errr(__FILE__, __LINE__)

        CALL comp_laplacep_avg(field, "laplaceP_AVG", dt)
        field%arr = field%arr(:)**2
        field%name = "laplaceP_SQR_AVG"
        field%units = field%units*2
    END SUBROUTINE comp_laplacep_sqr_avg


    SUBROUTINE calclpp_grid(kk, jj, ii, lpp, p, bp, dx, dy, dz, &
            nfro, nbac, nrgt, nlft, nbot, ntop)
        ! Subroutine arguments
        INTEGER(intk), INTENT(in) :: kk, jj, ii
        REAL(realk), INTENT(inout) :: lpp(kk, jj, ii)
        REAL(realk), INTENT(in) :: p(kk, jj, ii)
        REAL(realk), INTENT(in) :: bp(kk, jj, ii)
        REAL(realk), INTENT(in) :: dx(ii), dy(jj), dz(kk)
        INTEGER(intk), INTENT(in) :: nfro, nbac, nrgt, nlft, nbot, ntop

        ! Local variables
        INTEGER(intk) :: k, j, i
        INTEGER(intk) :: ifr, iba, jri, jle, kbo, kto
        REAL(realk) :: flfr, flba, flri, flle, flto, flbo
        REAL(realk) :: d2pdx2, d2pdy2, d2pdz2

        ! 2 = FIX, 5 = NOS, 6 = SLI, 19 = CO1
        IF (nfro == 2 .OR. nfro == 5 .OR. nfro == 6 .OR. nfro == 19) ifr = 1
        IF (nrgt == 2 .OR. nrgt == 5 .OR. nrgt == 6 .OR. nrgt == 19) jri = 1
        IF (nbot == 2 .OR. nbot == 5 .OR. nbot == 6 .OR. nbot == 19) kbo = 1
        IF (nbac == 2 .OR. nbac == 5 .OR. nbac == 6 .OR. nbac == 19) iba = 1
        IF (nlft == 2 .OR. nlft == 5 .OR. nlft == 6 .OR. nlft == 19) jle = 1
        IF (ntop == 2 .OR. ntop == 5 .OR. ntop == 6 .OR. ntop == 19) kto = 1

        ! lpp is declared INTENT(inout), therefore there are no need to set
        ! buffers (where no lpp is computed) to zero. If it were INTENT(out)
        ! the entire field should have been completely defined inside this
        ! routine

        DO i = 3, ii-2
            IF (ifr == 1 .AND. i == 3) THEN
                flfr = 1.0
            ELSE
                flfr = 0.0
            ENDIF
            IF (iba == 1 .AND. i == ii-2) THEN
                flba = 1.0
            ELSE
                flba = 0.0
            ENDIF

            DO j = 3, jj-2
                IF (jri == 1 .AND. j == 3) THEN
                    flri = 1.0
                ELSE
                    flri = 0.0
                ENDIF
                IF (jle == 1 .AND. j == jj-2) THEN
                    flle = 1.0
                ELSE
                    flle = 0.0
                ENDIF

                DO k = 3, kk-2
                    IF (kbo == 1 .AND. k == 3) THEN
                        flbo = 1.0
                    ELSE
                        flbo = 0.0
                    ENDIF
                    IF(kto == 1 .AND. k == kk-2) THEN
                        flto = 1.0
                    ELSE
                        flto = 0.0
                    ENDIF

                    d2pdx2 = ((1.0-flba)*bp(k, j, i+1)*p(k, j, i+1) &
                        + (-bp(k, j, i+1)-bp(k, j, i-1)+flfr+flba)*p(k, j, i) &
                        + (1.0-flfr)*bp(k, j, i-1)*p(k, j, i-1)) &
                        /(dx(i)*dx(i-1))

                    d2pdy2 = ((1.0-flle)*bp(k, j+1, i)*p(k, j+1, i) &
                        + (-bp(k, j+1, i)-bp(k, j-1, i)+flri+flle)*p(k, j, i) &
                        + (1.0-flri)*bp(k, j-1, i)*p(k, j-1, i)) &
                        /(dy(j)*dy(j-1))

                    d2pdz2 = ((1.0-flto)*bp(k+1, j, i)*p(k+1, j, i) &
                        + (-bp(k+1, j, i)-bp(k-1, j, i)+flto+flbo)*p(k, j, i) &
                        + (1.0-flbo)*bp(k-1, j, i)*p(k-1, j, i)) &
                        /(dz(k)*dz(k-1))

                    lpp(k, j, i) = bp(k, j, i)*(d2pdx2 + d2pdy2 + d2pdz2)
                END DO
            END DO
        END DO
    END SUBROUTINE calclpp_grid

END MODULE flowstat_mod
