
MODULE sip_classic_mod

    USE core_mod
    use laplacephi_mod

    IMPLICIT NONE(type, external)

    PRIVATE

    LOGICAL :: is_init = .FALSE.

    PUBLIC :: sip_classic_init, sip_classic_finish, sipiter1, sipiter2

CONTAINS

    SUBROUTINE sip_classic_init()

        ! Subroutine arguments
        ! none...

        ! Local variables
        INTEGER(intk) :: igr, igrid
        INTEGER(intk) :: k, j, i
        INTEGER(intk) :: kk, jj, ii
        REAL(realk), PARAMETER :: alfa = 0.92
        REAL(realk) :: p1, p2, p3
        REAL(realk), POINTER, CONTIGUOUS :: bp(:, :, :)
        REAL(realk), POINTER, CONTIGUOUS :: lb(:, :, :), lw(:, :, :), &
            ls(:, :, :), ue(:, :, :), un(:, :, :), ut(:, :, :), lpr(:, :, :)
        REAL(realk), POINTER, CONTIGUOUS :: ae(:),  aw(:),  an(:), as(:), &
            at(:), ab(:), ap(:, :, :)


        ! Conducting the incomplete LU decomposition for all local grids
        DO igr = 1, nmygrids
            igrid = mygrids(igr)
            CALL get_mgdims(kk, jj, ii, igrid)

            CALL get_fieldptr(lw, "SIPLW", igrid)
            CALL get_fieldptr(ls, "SIPLS", igrid)
            CALL get_fieldptr(lb, "SIPLB", igrid)
            CALL get_fieldptr(ue, "SIPUE", igrid)
            CALL get_fieldptr(un, "SIPUN", igrid)
            CALL get_fieldptr(ut, "SIPUT", igrid)
            CALL get_fieldptr(lpr, "SIPLPR", igrid)

            CALL get_fieldptr(aw, "GSAW", igrid)
            CALL get_fieldptr(ae, "GSAE", igrid)
            CALL get_fieldptr(as, "GSAS", igrid)
            CALL get_fieldptr(an, "GSAN", igrid)
            CALL get_fieldptr(ab, "GSAB", igrid)
            CALL get_fieldptr(at, "GSAT", igrid)

            CALL get_fieldptr(ap, "GSAP", igrid)

            CALL get_fieldptr(bp, "BP", igrid)

            DO i = 3, ii-2
                DO j = 3, jj-2
                    DO k = 3, kk-2
                        lw(k, j, i) = aw(i)*bu(k, j, i-1) &
                            /(1.0 + alfa*(un(k, j, i-1) + ut(k, j, i-1)))

                        ls(k, j, i) = as(j)*bv(k, j-1, i) &
                            /(1.0 + alfa*(ue(k, j-1, i) + ut(k, j-1, i)))

                        lb(k, j, i) = ab(k)*bw(k-1, j, i) &
                            /(1.0 + alfa*(un(k-1, j, i) + ue(k-1, j, i)))

                        p1 = alfa*(lb(k, j, i)*ue(k-1, j, i) &
                            + ls(k, j, i)*ue(k, j-1, i))
                        p2 = alfa*(lb(k, j, i)*un(k-1, j, i) &
                            + lw(k, j, i)*un(k, j, i-1))
                        p3 = alfa*(lw(k, j, i)*ut(k, j, i-1) &
                            + ls(k, j, i)*ut(k, j-1, i))

                        lpr(k, j, i) = 1.0/(ap(k, j, i) + p1 + p2 + p3 &
                            - lb(k, j, i)*ut(k-1, j, i) &
                            - lw(k, j, i)*ue(k, j, i-1) &
                            - ls(k, j, i)*un(k, j-1, i) &
                            + 1.0e-20)

                        ue(k, j, i) = (ae(i)*bu(k, j, i) - p1)*lpr(k, j, i)
                        un(k, j, i) = (an(j)*bv(k, j, i) - p2)*lpr(k, j, i)
                        ut(k, j, i) = (at(k)*bw(k, j, i) - p3)*lpr(k, j, i)
                    END DO
                END DO
            END DO

            DO i = 3, ii-2
                DO j = 3, jj-2
                    DO k = 3, kk-2
                        lw(k, j, i) = lw(k, j, i)*lpr(k, j, i)
                    END DO
                    DO k = 3, kk-2
                        ls(k, j, i) = ls(k, j, i)*lpr(k, j, i)
                    END DO
                    DO k = 3, kk-2
                        lb(k, j, i) = lb(k, j, i)*lpr(k, j, i)
                    END DO
                END DO
            END DO
        END DO

        is_init = .TRUE.

    CONTAINS

        ! These routines are not protected against out-of-bounds and not valid
        ! for the edges ii, jj, kk - but in the scope above that is OK
        PURE REAL(realk) FUNCTION bu(k, j, i)
            INTEGER(intk), INTENT(in) :: k, j, i
            bu = bp(k, j, i)*bp(k, j, i+1)
        END FUNCTION bu

        PURE REAL(realk) FUNCTION bv(k, j, i)
            INTEGER(intk), INTENT(in) :: k, j, i
            bv = bp(k, j, i)*bp(k, j+1, i)
        END FUNCTION bv

        PURE REAL(realk) FUNCTION bw(k, j, i)
            INTEGER(intk), INTENT(in) :: k, j, i
            bw = bp(k, j, i)*bp(k+1, j, i)
        END FUNCTION bw

    END SUBROUTINE sip_classic_init


    SUBROUTINE sip_classic_finish()

        ! Subroutine arguments
        ! none...

        ! Local variables
        ! none...

        is_init = .FALSE.

    END SUBROUTINE sip_classic_finish


    PURE SUBROUTINE sipiter1(kk, jj, ii, rhs, res, lw, ls, lb, lpr)
        ! Subroutine arguments
        INTEGER(intk), INTENT(in) :: kk, jj, ii
        REAL(realk), INTENT(inout) :: res(kk, jj, ii)
        REAL(realk), INTENT(in) :: rhs(kk, jj, ii)
        REAL(realk), INTENT(in) :: lw(kk, jj, ii), ls(kk, jj, ii), &
            lb(kk, jj, ii), lpr(kk, jj, ii)

        ! Local variables
        INTEGER(intk) :: k, j, i

        DO i = 3, ii-2
            DO j = 3, jj-2
                DO k = 3, kk-2
                    res(k, j, i) = (rhs(k, j, i) + res(k, j, i))*lpr(k, j, i)
                    res(k, j, i) = res(k, j, i) - lw(k, j, i)*res(k, j, i-1) &
                        - ls(k, j, i)*res(k, j-1, i)
                END DO
                DO k = 3, kk-2
                    res(k, j, i) = res(k, j, i) - lb(k, j, i)*res(k-1, j, i)
                END DO
            END DO
        END DO
    END SUBROUTINE sipiter1


    PURE SUBROUTINE sipiter2(kk, jj, ii, phi, res, ue, un, ut)
        ! Subroutine arguments
        INTEGER(intk), INTENT(in) :: kk, jj, ii
        REAL(realk), INTENT(inout) :: phi(kk, jj, ii), res(kk, jj, ii)
        REAL(realk), INTENT(in) :: ue(kk, jj, ii), un(kk, jj, ii), &
            ut(kk, jj, ii)

        ! Local variables
        INTEGER(intk) :: k, j, i

        DO i = ii-2, 3, -1
            DO j = jj-2, 3, -1
                DO k = 3, kk-2
                    res(k, j, i) = res(k, j, i) - un(k, j, i)*res(k, j+1, i) &
                        - ue(k, j, i)*res(k, j, i+1)
                END DO
                DO k = kk-2, 3, -1
                    res(k, j, i) = res(k, j, i) - ut(k, j, i)*res(k+1, j, i)
                END DO
            END DO
        END DO

        DO i = 3, ii-2
            DO j = 3, jj-2
                DO k = 3, kk-2
                    phi(k, j, i) = phi(k, j, i) + res(k, j, i)
                END DO
            END DO
        END DO
    END SUBROUTINE sipiter2

END MODULE sip_classic_mod
