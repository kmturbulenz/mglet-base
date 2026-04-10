    SUBROUTINE sip(ilevel, iloop, dp, res, rhs, siplw, sipls, siplb, &
            sipue, sipun, siput, siplpr, bp)

        USE hyperplane_mod, ONLY: mip_sip_list, idx_sip_list

        ! Subroutine arguments
        INTEGER(intk), INTENT(in) :: ilevel
        INTEGER(intk), INTENT(in) :: iloop
        TYPE(field_t), INTENT(inout) :: dp
        TYPE(field_t), INTENT(inout) :: res
        TYPE(field_t), INTENT(in) :: rhs
        TYPE(field_t), INTENT(in) :: siplw
        TYPE(field_t), INTENT(in) :: sipls
        TYPE(field_t), INTENT(in) :: siplb
        TYPE(field_t), INTENT(in) :: sipue
        TYPE(field_t), INTENT(in) :: sipun
        TYPE(field_t), INTENT(in) :: siput
        TYPE(field_t), INTENT(in) :: siplpr
        TYPE(field_t), INTENT(in), OPTIONAL :: bp

        ! Local variables
        INTEGER(intk) :: i, il, igrid, ip3
        INTEGER(intk) :: kk, jj, ii, n3dmin, n3dmax, m, lm, lp, ip, iacc, &
            idx, idx_km, idx_jm, idx_im, idx_kp, idx_jp, idx_ip, i2, j2, k2
        REAL(realk), POINTER, CONTIGUOUS :: dp_1d_p(:), res_1d_p(:), &
            rhs_1d_p(:)

        CALL laplacephi_level(ilevel, res, dp, bp)

        ! private(igrid, i, kk, jj, ii, ip3, n3dmin, n3dmax)
        ! private(iacc, idx, idx_km, idx_jm, idx_im)


        WRITE(*, *) "A"

        !$omp target teams loop bind(teams) private(igrid, i, kk, jj, ii, ip3, n3dmin, n3dmax) map(to: dp%arr, res%arr) map(to: rhs%arr, mip_sip_list, idx_sip_list)
        DO il = 1, nmygridslvl(ilevel)

            igrid = mygridslvl(il, ilevel)
            kk = gridinfo(igrid)%kk
            jj = gridinfo(igrid)%jj
            ii = gridinfo(igrid)%ii
            i = igrid
            ip3 = kk*jj*ii*(i-1)

            n3dmin = 3 + 3 + 3
            n3dmax = (ii-2) + (jj-2) + (kk-2)

            WRITE(*, *) igrid, ii, jj, kk, ip3, n3dmin, n3dmax

            !$omp parallel private(lm, lp, ip, iacc, idx, idx_km, idx_jm, idx_im) num_threads(1000)

            DO m = n3dmin, n3dmax

                lm = mip_sip_list(i)%arr(m)
                lp = mip_sip_list(i)%arr(m+1) - lm

                WRITE(*, *) "A"

                ip = omp_get_thread_num() + 1

                IF (ip <= lp) THEN

                    ! Computing the contiguous access index
                    iacc = lm + ip

                    ! Computing the required indices
                    idx = idx_sip_list(i)%arr(iacc) + ip3
                    idx_km = idx - 1
                    idx_jm = idx - kk
                    idx_im = idx - (kk*jj)

                    ! Accounting for RHS
                    res%arr(idx) = (rhs%arr(idx) + &
                        res%arr(idx)) * lpr_sip_list(i)%arr(iacc)

                    ! Performing the forward substitution
                    res%arr(idx) = res%arr(idx) - &
                        lb_sip_list(i)%arr(iacc) * res%arr(idx_km) - &
                        ls_sip_list(i)%arr(iacc) * res%arr(idx_jm) - &
                        lw_sip_list(i)%arr(iacc) * res%arr(idx_im)

                END IF

                !$omp barrier

            END DO

            !$omp end parallel

        END DO
        !$omp end target teams loop


        WRITE(*, '("Finished forward substitution")')
        IF (iloop < ninner) THEN
            CALL connect(ilevel, 1, s1=res)
        ELSE
            CALL connect(ilevel, 1, s1=res, forward=-1)
        END IF


        !$omp target teams loop bind(teams) private(igrid, i, kk, jj, ii, ip3, n3dmin, n3dmax)
        DO il = 1, nmygridslvl(ilevel)

            igrid = mygridslvl(il, ilevel)
            kk = gridinfo(igrid)%kk
            jj = gridinfo(igrid)%jj
            ii = gridinfo(igrid)%ii
            i = igrid
            ip3 = kk*jj*ii*(i-1)   ! hack

            n3dmin = 3 + 3 + 3
            n3dmax = (ii-2) + (jj-2) + (kk-2)

            !$omp parallel private(lm, lp, ip, iacc, idx, idx_km, idx_jm, idx_im) num_threads(900)

            DO m = n3dmax, n3dmin, -1

                lm = mip_sip_list(i)%arr(m)
                lp = mip_sip_list(i)%arr(m+1) - lm

                ip = omp_get_thread_num() + 1

                IF (ip <= lp) THEN

                    ! Computing the contiguous access index
                    iacc = lm + ip

                    ! Computing the required indices
                    idx = idx_sip_list(i)%arr(iacc) + ip3
                    idx_kp = idx + 1
                    idx_jp = idx + kk
                    idx_ip = idx + (kk*jj)

                    ! Performing the backward substitution
                    res%arr(idx) = res%arr(idx) - &
                        ut_sip_list(i)%arr(iacc) * res%arr(idx_kp) - &
                        un_sip_list(i)%arr(iacc) * res%arr(idx_jp) - &
                        ue_sip_list(i)%arr(iacc) * res%arr(idx_ip)

                END IF

                !$omp barrier

            END DO

            !$omp end parallel

            !$omp parallel do private(idx) collapse(3)
            DO i2 = 3, ii-2
                DO j2 = 3, jj-2
                    DO k2 = 3, kk-2
                        idx = k2 + (j2-1)*kk + (i2-1)*kk*jj + ip3
                        dp%arr(idx) = dp%arr(idx) + res%arr(idx)
                    END DO
                END DO
            END DO
            !$omp end parallel do


        END DO
        !$omp end target teams loop

    END SUBROUTINE sip