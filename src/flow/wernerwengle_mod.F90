MODULE wernerwengle_mod
    USE core_mod
    USE flowcore_mod, ONLY: rho, gmol

    IMPLICIT NONE(type, external)
    PRIVATE

    REAL(realk), PROTECTED :: cwa, cwb
    real(realk), PROTECTED :: cpo1, cpo2, cpo3, cpo4, cpo5, cpo6, &
        cpo7, cpo8, cpo9, cpo10, cpo11, cpo12

    PUBLIC :: init_wernerwengle, finish_wernerwengle, gradp2, tauwin, qwallfix
CONTAINS
    SUBROUTINE init_wernerwengle()
        cwa = 8.3_realk
        cwb = 1.0_realk/7.0_realk

        cpo1 = 1.0_realk - cwb
        cpo2 = 1.0_realk + cwb
        cpo3 = 1.0_realk/cpo2
        cpo4 = cpo2/cwa*(gmol/rho)**cwb
        cpo5 = 0.5_realk*cpo1*(cwa**(cpo2/cpo1))*(gmol/rho)**cpo2
        cpo6 = 0.5_realk*gmol/rho*cwa**(2.0_realk/cpo1)
        cpo7 = cwa*cwb/((gmol/rho)**cwb)
        cpo8 = 2.0_realk/cpo2
        cpo9 = -cpo1*cpo6/cpo2
        cpo10 = gmol/rho*cwa**(1.0_realk/cpo1)
        cpo11 = 2.0_realk*cpo6*cwb*cpo1/(2.0_realk** cwb)
        cpo12 = cwb*cpo2/(2.0_realk**(cwb-1.0_realk))
    END SUBROUTINE init_wernerwengle


    SUBROUTINE finish_wernerwengle()
        CONTINUE
    END SUBROUTINE finish_wernerwengle


    PURE ELEMENTAL REAL(realk) FUNCTION gradp2(uquer, dds)
        !$omp declare simd(gradp2)

        ! Computes the gradient dds/2 away from the wall using the WW wall
        ! model

        ! Function arguments
        REAL(realk), INTENT(in) :: uquer, dds

        ! Local variables
        REAL(realk) :: uquern, vz, rdds

        uquern = ABS(uquer)
        vz = SIGN(1.0_realk, uquer)
        rdds = 1.0_realk/dds
        IF (uquern >= cpo6*rdds) THEN
            ! Fully turbulent regime
            gradp2 = rdds*(cpo12*uquern + cpo11*rdds)*vz
        ELSE
            ! Laminar sublayer
            gradp2 = 2.0_realk*uquern/dds*vz
        END IF
    END FUNCTION gradp2


    PURE ELEMENTAL REAL(realk) FUNCTION tauwin(uquer, dds)
        !$omp declare simd(tauwin)

        ! Function arguments
        REAL(realk), INTENT(IN) :: uquer, dds

        ! Local variables
        REAL(realk) :: vz, uquern
        REAL(realk) :: ddsb

        vz = SIGN(1.0_realk, uquer)
        uquern = ABS(uquer)

        IF (uquern >= cpo6/dds) THEN
            ! Fully turbulent regime
            ddsb = dds**(-cwb)
            tauwin  = vz*rho*(ddsb*(cpo4*uquern + cpo5/dds))**cpo8
        ELSE
            ! Laminar sublayer
            tauwin  = vz*2.0*gmol*uquern/dds
        END IF
    END FUNCTION tauwin


    ! Experimental and untested. Should be tested and verified before actual
    ! usage!
    PURE ELEMENTAL REAL(realk) FUNCTION qwallfix(tbound, tfluid, uquer, dds, &
            prmol)
        !$omp declare simd(qwallfix)

        ! Function arguments
        REAL(realk), INTENT(IN) :: tbound, tfluid, uquer, dds, prmol

        ! Local variables
        REAL(realk) :: uquern, ddsb, tauwin, sm, fsm

        uquern = ABS(uquer)
        IF (uquern >= cpo6/dds) THEN
            ! Fully turbulent regime
            ddsb = dds**(-cwb)
            tauwin = rho*(ddsb*(cpo4*uquern + cpo5/dds))**cpo8
            ! cwa**(1.0/(1.0-cwb)) = 11.8 => intersection between linear and
            ! exponential in plus-units sm is the "real" wall distance [m]
            ! of the switching point
            sm = gmol / tauwin**0.5 * cwa**(1.0/(1.0-cwb))
            fsm = 1/dds**2 *( 0.5*sm**2.0*prmol + (tauwin/gmol)**(-cpo1) &
                * cwa/cpo2 * (dds**cpo2 - sm**cpo2) )

            qwallfix = gmol / rho * (tbound - tfluid) / dds / fsm
        ELSE
            ! Laminar sublayer
            qwallfix = gmol/rho/prmol*(tbound - tfluid)*2.0/dds
        END IF
    END FUNCTION qwallfix
END MODULE wernerwengle_mod
