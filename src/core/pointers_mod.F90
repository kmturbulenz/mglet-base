MODULE pointers_mod
    USE precision_mod, ONLY: intk
    USE err_mod, ONLY: errr
    USE comms_mod, ONLY: myid
    USE grids_mod, ONLY: ngrid, get_mgdims, mygrids, nmygrids, idprocofgrd

    IMPLICIT NONE(type, external)
    PRIVATE

    INTEGER(intk), PROTECTED :: idim3d, idim2d, idim1d
    INTEGER(intk), ALLOCATABLE, PROTECTED :: ip3d(:), ip2d(:), ip1d(:)

    PUBLIC :: init_pointers, finish_pointers, get_ip1, &
        get_ip3, get_ip3n, get_ibb, get_ibbn, &
        idim3d, idim2d, get_len3

CONTAINS
    SUBROUTINE init_pointers()
        idim3d = 0
        idim2d = 0
        idim1d = 0

        ALLOCATE(ip3d(ngrid))
        ALLOCATE(ip2d(ngrid))
        ALLOCATE(ip1d(ngrid))

        ip3d = 0
        ip2d = 0
        ip1d = 0

        BLOCK
            ! Initialize and set pointers for all grids this process owns
            INTEGER(intk) :: i, igrid, kk, jj, ii
            INTEGER(intk) :: nsize1d, nsize2d, nsize3d

            DO i = 1, nmygrids
                igrid = mygrids(i)

                CALL get_mgdims(kk, jj, ii, igrid)

                nsize3d = ii*jj*kk
                nsize2d = MAX(ii*jj, ii*kk, jj*kk)
                nsize1d = MAX(ii, jj, kk)

                ip3d(igrid) = idim3d + 1
                ip2d(igrid) = idim2d + 1
                ip1d(igrid) = idim1d + 1

                idim3d = idim3d + nsize3d
                idim2d = idim2d + nsize2d
                idim1d = idim1d + nsize1d
            END DO
        END BLOCK

        IF (myid == 0) THEN
            WRITE(*, '("ARRAY DIMENSIONS:")')
            WRITE(*, '("    idim3d:        ", I0)') idim3d
            WRITE(*, '("    idim2d:        ", I0)') idim2d
            WRITE(*, '("    idim1d:        ", I0)') idim1d
            WRITE(*, '()')
        END IF
    END SUBROUTINE init_pointers


    SUBROUTINE finish_pointers()
        idim3d = 0
        idim2d = 0
        idim1d = 0

        DEALLOCATE(ip3d)
        DEALLOCATE(ip2d)
        DEALLOCATE(ip1d)
    END SUBROUTINE finish_pointers


    SUBROUTINE get_ip3(ip3, igrid)
        INTEGER(intk), INTENT(in) :: igrid
        INTEGER(intk), INTENT(out) :: ip3

#ifdef _MGLET_DEBUG_
        IF (myid /= idprocofgrd(igrid)) THEN
            CALL errr(__FILE__, __LINE__)
        END IF
#endif

        ip3 = ip3d(igrid)
    END SUBROUTINE get_ip3


    SUBROUTINE get_ip3n(ip3n, ncomp, igrid)
        INTEGER(intk), INTENT(in) :: igrid
        INTEGER(intk), INTENT(in) :: ncomp
        INTEGER(intk), INTENT(out) :: ip3n

        IF (myid /= idprocofgrd(igrid)) THEN
            CALL errr(__FILE__, __LINE__)
        END IF

        ip3n = ncomp*ip3d(igrid) - (ncomp-1)
    END SUBROUTINE get_ip3n


    SUBROUTINE get_ibb(ibb, igrid)
        INTEGER(intk), INTENT(in) :: igrid
        INTEGER(intk), INTENT(out) :: ibb

        IF (myid /= idprocofgrd(igrid)) THEN
            CALL errr(__FILE__, __LINE__)
        END IF

        ibb = 2*ip2d(igrid) - 1
    END SUBROUTINE get_ibb


    SUBROUTINE get_ibbn(ibbt, ncomp, igrid)
        INTEGER(intk), INTENT(in) :: igrid, ncomp
        INTEGER(intk), INTENT(out) :: ibbt

        IF (myid /= idprocofgrd(igrid)) THEN
            CALL errr(__FILE__, __LINE__)
        END IF

        ibbt = 2*ncomp*ip2d(igrid) - (2*ncomp-1)
    END SUBROUTINE get_ibbn


    SUBROUTINE get_ip1(ip1, igrid)
        INTEGER(intk), INTENT(in) :: igrid
        INTEGER(intk), INTENT(out) :: ip1

        IF (myid /= idprocofgrd(igrid)) THEN
            CALL errr(__FILE__, __LINE__)
        END IF

        ip1 = ip1d(igrid)
    END SUBROUTINE get_ip1


    SUBROUTINE get_len3(len, igrid)
        INTEGER(intk), INTENT(in) :: igrid
        INTEGER(intk), INTENT(out) :: len

        INTEGER(intk) :: kk, jj, ii

        IF (myid /= idprocofgrd(igrid)) THEN
            CALL errr(__FILE__, __LINE__)
        END IF

        CALL get_mgdims(kk, jj, ii, igrid)
        len = kk*jj*ii
    END SUBROUTINE get_len3
END MODULE pointers_mod
