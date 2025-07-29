MODULE tensormath_mod
    USE precision_mod, ONLY: realk, intk

    IMPLICIT NONE (type, external)
    PRIVATE

    TYPE :: tensor_t
        REAL(realk) :: g(3, 3)
    CONTAINS
        PRIVATE

        PROCEDURE, PUBLIC, NON_OVERRIDABLE :: abst
        PROCEDURE, PUBLIC, NON_OVERRIDABLE :: t
        PROCEDURE, PUBLIC, NON_OVERRIDABLE :: sqr
        PROCEDURE, PUBLIC, NON_OVERRIDABLE :: trace
        PROCEDURE, PUBLIC, NON_OVERRIDABLE :: det
        PROCEDURE, PUBLIC, NON_OVERRIDABLE :: eig_a, eig_b

#if defined __GNUC__
        PROCEDURE :: add
        PROCEDURE :: subtract
        PROCEDURE :: multiply
#else
        PROCEDURE, NON_OVERRIDABLE :: add
        PROCEDURE, NON_OVERRIDABLE :: subtract
        PROCEDURE, NON_OVERRIDABLE :: multiply
#endif

        GENERIC, PUBLIC :: OPERATOR(+) => add
        GENERIC, PUBLIC :: OPERATOR(-) => subtract
        GENERIC, PUBLIC :: OPERATOR(*) => multiply
    END TYPE tensor_t

    INTERFACE OPERATOR(*)
        MODULE PROCEDURE scamul
    END INTERFACE

    TYPE(tensor_t), PARAMETER :: ident = tensor_t(RESHAPE([1.0, 0.0, 0.0, &
        0.0, 1.0, 0.0, 0.0, 0.0, 1.0], [3, 3]))

    PUBLIC :: tensor_t, OPERATOR(*)

CONTAINS

    PURE REAL(realk) FUNCTION abst(this)
        ! Computes:
        ! abst = |T_ij| = sqrt(Tij*Tij)

        ! Function arguments
        CLASS(tensor_t), INTENT(in) :: this

        ! Local variables
        REAL(realk) :: sumsqr
        INTEGER(intk) :: i, j

        sumsqr = 0.0
        DO i = 1, 3
            DO j = 1, 3
                sumsqr = sumsqr + this%g(j, i)**2
            END DO
        END DO
        abst = SQRT(sumsqr)
    END FUNCTION abst


    PURE TYPE(tensor_t) FUNCTION t(this)
        ! Transpose a tensor

        ! Function arguments
        CLASS(tensor_t), INTENT(in) :: this

        t%g = TRANSPOSE(this%g)
    END FUNCTION t


    PURE TYPE(tensor_t) FUNCTION sqr(this)
        ! Multiply a tensor with itself (standard matrix multiplication)
        ! T_ij^2 = T_ik*T_kj

        ! Function arguments
        CLASS(tensor_t), INTENT(in) :: this

        INTEGER(intk) :: i, j

        DO i = 1, 3
            DO j = 1, 3
                sqr%g(j, i) = this%g(1, i)*this%g(j, 1) &
                    + this%g(2, i)*this%g(j, 2) + this%g(3, i)*this%g(j, 3)
            END DO
        END DO
    END FUNCTION sqr


    PURE REAL(realk) FUNCTION trace(this)
        ! Trace of a tensor

        ! Function arguments
        CLASS(tensor_t), INTENT(in) :: this

        trace = this%g(1, 1) + this%g(2, 2) + this%g(3, 3)
    END FUNCTION trace


    PURE TYPE(tensor_t) FUNCTION add(this, that)
        ! Add two tensors together, element by element

        ! Function arguments
        CLASS(tensor_t), INTENT(in) :: this
        CLASS(tensor_t), INTENT(in) :: that

        add%g = this%g + that%g
    END FUNCTION add


    PURE TYPE(tensor_t) FUNCTION subtract(this, that)
        ! Subtract one tensor from another, element by element

        ! Function arguments
        CLASS(tensor_t), INTENT(in) :: this
        CLASS(tensor_t), INTENT(in) :: that

        subtract%g = this%g - that%g
    END FUNCTION subtract


    PURE TYPE(tensor_t) FUNCTION multiply(this, that)
        ! Multiply two tensors, standard matrix multiplication rules
        ! R_ij = T_ik*G_kj

        ! Function arguments
        CLASS(tensor_t), INTENT(in) :: this
        CLASS(tensor_t), INTENT(in) :: that

        ! Local variables
        INTEGER(intk) :: i, j, k

        ! g(j, i) -> j=column, i=row
        multiply%g = 0.0
        DO i = 1, 3
            DO j = 1, 3
                DO k = 1, 3
                    multiply%g(j, i) = multiply%g(j, i) + &
                        this%g(k, i)*that%g(j, k)
                END DO
            END DO
        END DO
    END FUNCTION multiply


    PURE REAL(realk) FUNCTION det(this)
        ! Computes determinant

        ! Function arguments
        CLASS(tensor_t), INTENT(in) :: this

        ! g(j, i) -> j=column, i=row
        det = this%g(1, 1)*(this%g(2, 2)*this%g(3, 3) &
            - this%g(3, 2)*this%g(2, 3)) &
            - this%g(2, 1)*(this%g(1, 2)*this%g(3, 3) &
            - this%g(3, 2)*this%g(1, 3)) &
            + this%g(3, 1)*(this%g(1, 2)*this%g(2, 3) &
            - this%g(2, 2)*this%g(1, 3))
    END FUNCTION det


    PURE SUBROUTINE eig_a(this, eig1, eig2, eig3)
        ! Eigenvalue calculation
        ! from: https://en.wikipedia.org/wiki/Eigenvalue_algorithm
        ! For SYMMETRIC tensors/matrices
        !
        ! Warning: There are no sanity check if the symmetry requirements
        ! is fulfilled!

        ! Subroutine arguments
        CLASS(tensor_t), INTENT(in) :: this
        REAL(realk), INTENT(out) :: eig1, eig2, eig3

        ! Local variables
        REAL(realk), PARAMETER :: pi = 4.0_realk*ATAN(1.0_realk)
        REAL(realk) :: e1, e2, e3
        REAL(realk) :: p, p1, p2, q, r, phi
        TYPE(tensor_t) :: B

        p1 = this%g(1, 2)**2 + this%g(1, 3)**2 + this%g(2, 3)**2
        IF (p1 < TINY(1.0_realk)) THEN
            ! Tensor is diagonal
            e1 = this%g(1, 1)
            e2 = this%g(2, 2)
            e3 = this%g(3, 3)

            eig1 = MAX(e1, e2, e3)
            eig3 = MIN(e1, e2, e3)
            eig2 = e1 + e2 + e3 - eig1 - eig3
        ELSE
            q = (1.0/3.0)*this%trace()
            p2 = (this%g(1, 1) - q)**2 &
                + (this%g(2, 2) - q)**2 &
                + (this%g(3, 3) - q)**2 + 2*p1
            p = SQRT(p2/6.0)
            B = (1.0/p)*(this - q*ident)
            r = 0.5*B%det()

            ! In exact arithmetic for a symmetric matrix  -1 <= r <= 1
            ! but computation error can leave it slightly outside this range.
            IF (r <= -1.0) THEN
                phi = pi/3.0
            ELSE IF (r >= 1.0) THEN
                phi = 0.0
            ELSE
                phi = ACOS(r)/3.0
            END IF

            ! The eigenvalues satisfy eig3 <= eig2 <= eig1
            eig1 = q + 2.0*p*COS(phi)
            eig3 = q + 2.0*p*COS(phi + (2.0*pi/3.0))
            eig2 = 3.0*q - eig1 - eig3  ! since trace(A) = eig1 + eig2 + eig3
        END IF
    END SUBROUTINE eig_a


    PURE SUBROUTINE eig_b(this, eig1, eig2, eig3)
        ! Eigenvalue calculation
        ! K. Hasan, P. Basser, D. Parker, and A. Alexander, "Analytical
        ! computation of the eigenvalues and eigenvectors in dt-mri,"
        ! J. Magn. Reson. 152, 41 (2001).
        !
        ! Warning: There are no sanity check if the symmetry requirements
        ! is fulfilled!

        USE simdfunctions_mod, ONLY: divide0

        ! Function arguments
        CLASS(tensor_t), INTENT(in) :: this
        REAL(realk), INTENT(out) :: eig1, eig2, eig3

        ! Local variables
        REAL(realk), PARAMETER :: pi = 4.0_realk*ATAN(1.0_realk)
        REAL(realk) :: i1, i2, i3
        REAL(realk) :: a1, a2, a3
        REAL(realk) :: arg
        TYPE(tensor_t) :: gg

        gg = this%sqr()

        ! Invariants
        i1 = this%trace()
        i2 = 0.5*(i1**2 - gg%trace())
        i3 = this%det()

        ! Angles
        a1 = i1**2/9.0_realk - i2/3.0_realk
        a2 = i1**3/27.0_realk - (i1*i2)/6.0_realk + i3/2.0_realk

        ! a3 = (1.0/3.0)*ACOS(a2/a1**(3.0/2.0))
        ! Due to numerical inaccuracies the argument to ACOS can be slightly
        ! outside the range -1.0 to 1.0, clip to this range
        arg = divide0(a2, a1**(3.0_realk/2.0_realk))
        arg = MAX(MIN(arg, 1.0_realk), -1.0_realk)
        a3 = (1.0_realk/3.0_realk)*ACOS(arg)

        ! Eigenvalues
        eig1 = i1/3.0_realk + 2.0_realk*SQRT(a1)*COS(a3)
        eig2 = i1/3.0_realk - 2.0_realk*SQRT(a1)*COS(pi/3.0_realk + a3)
        eig3 = i1/3.0_realk - 2.0_realk*SQRT(a1)*COS(pi/3.0_realk - a3)
    END SUBROUTINE eig_b


    PURE TYPE(tensor_t) FUNCTION scamul(a, that)
        ! Scale a tensor with a scalar:
        ! R_ij = a*T_ij

        ! Function arguments
        REAL(realk), INTENT(in) :: a
        CLASS(tensor_t), INTENT(in) :: that

        scamul%g = a*that%g
    END FUNCTION scamul

END MODULE tensormath_mod
