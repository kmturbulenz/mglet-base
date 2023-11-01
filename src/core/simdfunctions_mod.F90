MODULE simdfunctions_mod
    USE, INTRINSIC :: ISO_FORTRAN_ENV, ONLY: real32, real64
    USE precision_mod, ONLY: intk, realk

    IMPLICIT NONE(type, external)
    PRIVATE

    INTERFACE cube_root
        MODULE PROCEDURE :: cube_root_sp, cube_root_dp
    END INTERFACE cube_root

    INTERFACE cross
        MODULE PROCEDURE :: cross_sp, cross_dp
    END INTERFACE cross

    PUBLIC :: divide00, divide0, l_to_i, i_to_l, lcm, gcd, cube_root, cross
CONTAINS

    PURE ELEMENTAL REAL(realk) FUNCTION divide00(a, b, bp) RESULT(res)
    !$omp declare simd(divide00)
        REAL(realk), INTENT(in) :: a, b, bp

        IF (bp == 0.0_realk) THEN
            res = 0.0_realk
        ELSE
            res = a/b
        END IF
    END FUNCTION divide00


    PURE ELEMENTAL REAL(realk) FUNCTION divide0(a, b) RESULT(res)
    !$omp declare simd(divide0)
        REAL(realk), INTENT(in) :: a, b

        IF (b == 0.0_realk) THEN
            res = 0.0
        ELSE
            res = a/b
        END IF
    END FUNCTION divide0


    PURE ELEMENTAL INTEGER(intk) FUNCTION l_to_i(l) RESULT(i)
    !$omp declare simd(l_to_i)
        LOGICAL, INTENT(in) :: l

        IF (l) THEN
            i = 1
        ELSE
            i = 0
        END IF
    END FUNCTION l_to_i


    PURE ELEMENTAL LOGICAL FUNCTION i_to_l(i) RESULT(l)
    !$omp declare simd(i_to_l)
        INTEGER(intk), INTENT(in) :: i

        IF (i == 0) THEN
            l = .FALSE.
        ELSE
            l = .TRUE.
        END IF
    END FUNCTION i_to_l


    PURE ELEMENTAL INTEGER(intk) FUNCTION lcm(a, b)
    !$omp declare simd(lcm)
        INTEGER(intk), INTENT(in) :: a, b
        lcm = a*b/gcd(a, b)
    END FUNCTION lcm


    PURE ELEMENTAL INTEGER(intk) FUNCTION gcd(a, b)
    !$omp declare simd(gcd)
        INTEGER(intk), INTENT(in) :: a, b
        INTEGER(intk) :: aa, bb, t

        aa = a
        bb = b
        DO WHILE (bb /= 0)
            t = bb
            bb = MOD(aa, bb)
            aa = t
        END DO

        gcd = ABS(aa)
    END FUNCTION gcd


    PURE ELEMENTAL REAL(real32) FUNCTION cube_root_sp(a)
    !$omp declare simd(cube_root_sp)
        REAL(real32), INTENT(in) :: a

        INTERFACE
            ! float cbrtf(float arg )
            PURE FUNCTION cbrtf(arg) RESULT(res) BIND(C, NAME='cbrtf')
                USE, INTRINSIC :: ISO_C_BINDING, ONLY: C_FLOAT
                REAL(C_FLOAT), INTENT(IN), VALUE :: arg
                REAL(C_FLOAT) :: res
            END FUNCTION cbrtf
        END INTERFACE

        cube_root_sp = cbrtf(a)
    END FUNCTION cube_root_sp


    PURE ELEMENTAL REAL(real64) FUNCTION cube_root_dp(a)
    !$omp declare simd(cube_root_dp)
        REAL(real64), INTENT(in) :: a

        INTERFACE
            ! double cbrt(double arg )
            PURE FUNCTION cbrt(arg) RESULT(res) BIND(C, NAME='cbrt')
                USE, INTRINSIC :: ISO_C_BINDING, ONLY: C_DOUBLE
                REAL(C_DOUBLE), INTENT(IN), VALUE :: arg
                REAL(C_DOUBLE) :: res
            END FUNCTION cbrt
        END INTERFACE

        cube_root_dp = cbrt(a)
    END FUNCTION cube_root_dp


    PURE SUBROUTINE cross_sp(c, a, b)
        REAL(real32), INTENT(out) :: c(3)
        REAL(real32), INTENT(in) :: a(3), b(3)

        c(1) = a(2)*b(3) - a(3)*b(2)
        c(2) = a(3)*b(1) - a(1)*b(3)
        c(3) = a(1)*b(2) - a(2)*b(1)
    END SUBROUTINE cross_sp


    PURE SUBROUTINE cross_dp(c, a, b)
        REAL(real64), INTENT(out) :: c(3)
        REAL(real64), INTENT(in) :: a(3), b(3)

        c(1) = a(2)*b(3) - a(3)*b(2)
        c(2) = a(3)*b(1) - a(1)*b(3)
        c(3) = a(1)*b(2) - a(2)*b(1)
    END SUBROUTINE cross_dp

END MODULE simdfunctions_mod
