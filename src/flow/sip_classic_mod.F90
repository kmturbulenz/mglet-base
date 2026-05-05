
MODULE sip_classic_mod

    USE core_mod
    use laplacephi_mod

    IMPLICIT NONE(type, external)

    PRIVATE

    LOGICAL :: is_init = .FALSE.

    ! A-versions: simple versions not considering the BP field
    ! B-versions: IB versions using the BP field
    INTERFACE sipiter1
        MODULE PROCEDURE :: sipiter1_A, sipiter1_B
    END INTERFACE sipiter1

    PUBLIC :: sip_classic_init, sip_classic_finish, sip_classic, &
        sipiter0, sipiter1, sipiter2

        ! sipiter{0,1,2} decleared public to avoid "unused fucntion" warnings

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


    SUBROUTINE sip_classic(ilevel, iloop, ninner, dp, res, rhs, &
            siplw, sipls, siplb, siplpr, sipue, sipun, siput, bp)

        ! Subroutine arguments
        INTEGER(intk), INTENT(in) :: ilevel
        INTEGER(intk), INTENT(in) :: iloop
        INTEGER(intk), INTENT(in) :: ninner
        TYPE(field_t), INTENT(inout) :: dp
        TYPE(field_t), INTENT(inout) :: res
        TYPE(field_t), INTENT(in) :: rhs
        TYPE(field_t), INTENT(in) :: siplw
        TYPE(field_t), INTENT(in) :: sipls
        TYPE(field_t), INTENT(in) :: siplb
        TYPE(field_t), INTENT(in) :: siplpr
        TYPE(field_t), INTENT(in) :: sipue
        TYPE(field_t), INTENT(in) :: sipun
        TYPE(field_t), INTENT(in) :: siput
        TYPE(field_t), INTENT(in), OPTIONAL :: bp

        ! Local variables
        INTEGER(intk) :: i, igrid
        INTEGER(intk) :: kk, jj, ii
        REAL(realk), POINTER, CONTIGUOUS :: lw(:, :, :), ls(:, :, :), &
            lb(:, :, :), ue(:, :, :), un(:, :, :), ut(:, :, :), &
            lpr(:, :, :)
        REAL(realk), POINTER, CONTIGUOUS :: dp_p(:, :, :), res_p(:, :, :), &
            rhs_p(:, :, :)

        IF (.NOT. is_init) CALL errr(__FILE__, __LINE__)

        CALL laplacephi_level(ilevel, res, dp, bp)

        DO i = 1, nmygridslvl(ilevel)
            igrid = mygridslvl(i, ilevel)
            CALL get_mgdims(kk, jj, ii, igrid)

            CALL res%get_ptr(res_p, igrid)
            CALL rhs%get_ptr(rhs_p, igrid)

            CALL siplw%get_ptr(lw, igrid)
            CALL sipls%get_ptr(ls, igrid)
            CALL siplb%get_ptr(lb, igrid)
            CALL siplpr%get_ptr(lpr, igrid)

            CALL sipiter1(kk, jj, ii, rhs_p, res_p, lw, ls, lb, lpr)
            ! >>> corresponds to sipiter1_B
        END DO

        IF (iloop < ninner) THEN
            CALL connect(ilevel, 1, s1=res)
        ELSE
            CALL connect(ilevel, 1, s1=res, forward=-1)
        END IF

        DO i = 1, nmygridslvl(ilevel)
            igrid = mygridslvl(i, ilevel)
            CALL get_mgdims(kk, jj, ii, igrid)

            CALL dp%get_ptr(dp_p, igrid)
            CALL res%get_ptr(res_p, igrid)

            CALL sipue%get_ptr(ue, igrid)
            CALL sipun%get_ptr(un, igrid)
            CALL siput%get_ptr(ut, igrid)

            CALL sipiter2(kk, jj, ii, dp_p, res_p, ue, un, ut)
        END DO
    END SUBROUTINE sip_classic


    PURE SUBROUTINE sipiter0(kk, jj, ii, rhs, res, lpr)
        ! Subroutine arguments
        INTEGER(intk), INTENT(in) :: kk, jj, ii
        REAL(realk), INTENT(in) :: rhs(kk, jj, ii)
        REAL(realk), INTENT(inout) :: res(kk, jj, ii)
        REAL(realk), INTENT(in) :: lpr(kk, jj, ii)

        ! Local variables
        INTEGER(intk) :: k, j, i

        DO i = 3, ii-2
            DO j = 3, jj-2
                DO k = 3, kk-2
                    res(k, j, i) = (rhs(k, j, i) + res(k, j, i))*lpr(k, j, i)
                END DO
            END DO
        END DO
    END SUBROUTINE sipiter0


    PURE SUBROUTINE sipiter1_A(kk, jj, ii, rhs, res, lw, ls, lb)
        ! Subroutine arguments
        INTEGER(intk), INTENT(in) :: kk, jj, ii
        REAL(realk), INTENT(inout) :: res(kk, jj, ii)
        REAL(realk), INTENT(in) :: rhs(kk, jj, ii)
        REAL(realk), INTENT(in) :: lw(kk, jj, ii), ls(kk, jj, ii), &
            lb(kk, jj, ii)

        ! Local variables
        INTEGER(intk) :: k, j, i

        DO i = 3, ii-2
            DO j = 3, jj-2
                DO k = 3, kk-2
                    res(k, j, i) = res(k, j, i) - lw(k, j, i)*res(k, j, i-1) &
                        - ls(k, j, i)*res(k, j-1, i)
                END DO
                DO k = 3, kk-2
                    res(k, j, i) = res(k, j, i) - lb(k, j, i)*res(k-1, j, i)
                END DO
            END DO
        END DO
    END SUBROUTINE sipiter1_A


    PURE SUBROUTINE sipiter1_B(kk, jj, ii, rhs, res, lw, ls, lb, lpr)
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
    END SUBROUTINE sipiter1_B


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