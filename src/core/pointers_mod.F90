MODULE pointers_mod
    USE precision_mod, ONLY: intk
    USE err_mod, ONLY: errr
    USE comms_mod, ONLY: myid
    USE grids_mod, ONLY: ngrid, get_mgdims, mygrids, nmygrids, idprocofgrd, &
        get_bc_ctyp, nboconds

    IMPLICIT NONE(type, external)
    PRIVATE

    INTEGER(intk), PROTECTED :: idim3d, idimbb, idim1dx, idim1dy, idim1dz
    INTEGER(intk), ALLOCATABLE, PROTECTED :: ip3d(:)
    INTEGER(intk), ALLOCATABLE, PROTECTED :: ipbb(:, :)
    INTEGER(intk), ALLOCATABLE, PROTECTED :: ip1dx(:)
    INTEGER(intk), ALLOCATABLE, PROTECTED :: ip1dy(:)
    INTEGER(intk), ALLOCATABLE, PROTECTED :: ip1dz(:)
    !$omp declare target(ip3d, ip1dx, ip1dy, ip1dz)

    PUBLIC :: init_pointers, finish_pointers, get_ip3, get_ip3n, get_ibb, &
        idim3d, idimbb, get_len3, get_ipx, get_ipy, get_ipz, ip3d, ip1dx, &
        ip1dy, ip1dz

CONTAINS
    SUBROUTINE init_pointers()
        idim3d = 0
        idimbb = 0
        idim1dx = 0
        idim1dy = 0
        idim1dz = 0

        ALLOCATE(ip3d(ngrid))
        ALLOCATE(ipbb(6, ngrid))
        ALLOCATE(ip1dx(ngrid))
        ALLOCATE(ip1dy(ngrid))
        ALLOCATE(ip1dz(ngrid))

        ip3d = 0
        ipbb = 0
        ip1dx = 0
        ip1dy = 0
        ip1dz = 0

        BLOCK
            ! Initialize and set pointers for all grids this process owns
            INTEGER(intk) :: i, igrid, iface, ibocd, nbocd, kk, jj, ii
            INTEGER(intk) :: nsizebb, nsize3d
            CHARACTER(len=8) :: ctyp

            DO i = 1, nmygrids
                igrid = mygrids(i)
                CALL get_mgdims(kk, jj, ii, igrid)

                ! 3-D pointers are allocated for the entire grid
                nsize3d = ii*jj*kk
                ip3d(igrid) = idim3d + 1
                idim3d = idim3d + nsize3d

                ! 1-D pointers
                ip1dx(igrid) = idim1dx + 1
                idim1dx = idim1dx + ii
                ip1dy(igrid) = idim1dy + 1
                idim1dy = idim1dy + jj
                ip1dz(igrid) = idim1dz + 1
                idim1dz = idim1dz + kk

                ! BB (boundary buffer) pointers are allocated for every FIX,
                ! OP1, PAR, SIO and SWA buffer for every face. If a face has
                ! for example both FIX and SIO, only one buffer is
                ! allocated for that face.
                DO iface = 1, 6
                    nsizebb = 0
                    nbocd = nboconds(iface, igrid)
                    DO ibocd = 1, nbocd
                        CALL get_bc_ctyp(ctyp, ibocd, iface, igrid)
                        SELECT CASE (ctyp)
                        CASE ("FIX", "OP1", "PAR", "SIO", "SWA")
                            SELECT CASE (iface)
                            CASE (1, 2)
                                nsizebb = kk*jj
                            CASE (3, 4)
                                nsizebb = kk*ii
                            CASE (5, 6)
                                nsizebb = jj*ii
                            END SELECT
                        END SELECT
                    END DO

                    ! ipbb needs to be 0 for unallocated buffers
                    IF (nsizebb > 0) THEN
                        ipbb(iface, igrid) = idimbb + 1
                        idimbb = idimbb + nsizebb
                    END IF
                END DO
            END DO
        END BLOCK
        !$omp target enter data map(always, to: ip3d, ip1dx, ip1dy, ip1dz)

        IF (myid == 0) THEN
            WRITE(*, '("ARRAY DIMENSIONS:")')
            WRITE(*, '("    idim3d:        ", I0)') idim3d
            WRITE(*, '("    idimbb:        ", I0)') idimbb
            WRITE(*, '()')
        END IF
    END SUBROUTINE init_pointers


    SUBROUTINE finish_pointers()
        idim3d = 0
        idimbb = 0
        idim1dx = 0
        idim1dy = 0
        idim1dz = 0

        DEALLOCATE(ip3d)
        DEALLOCATE(ipbb)
        DEALLOCATE(ip1dx)
        DEALLOCATE(ip1dy)
        DEALLOCATE(ip1dz)
    END SUBROUTINE finish_pointers


    SUBROUTINE get_ip3(ip3, igrid)
        !$omp declare target
        INTEGER(intk), INTENT(in) :: igrid
        INTEGER(intk), INTENT(out) :: ip3

! #ifdef _MGLET_DEBUG_
!         IF (myid /= idprocofgrd(igrid)) THEN
!             CALL errr(__FILE__, __LINE__)
!         END IF
! #endif

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


    SUBROUTINE get_ipx(ipx, igrid)
        !$omp declare target
        INTEGER(intk), INTENT(in) :: igrid
        INTEGER(intk), INTENT(out) :: ipx

        ipx = ip1dx(igrid)
    END SUBROUTINE get_ipx


    SUBROUTINE get_ipy(ipy, igrid)
        !$omp declare target
        INTEGER(intk), INTENT(in) :: igrid
        INTEGER(intk), INTENT(out) :: ipy

        ipy = ip1dy(igrid)
    END SUBROUTINE get_ipy


    SUBROUTINE get_ipz(ipz, igrid)
        !$omp declare target
        INTEGER(intk), INTENT(in) :: igrid
        INTEGER(intk), INTENT(out) :: ipz

        ipz = ip1dz(igrid)
    END SUBROUTINE get_ipz


    SUBROUTINE get_ibb(ibb, iface, igrid)
        INTEGER(intk), INTENT(out) :: ibb
        INTEGER(intk), INTENT(in) :: iface
        INTEGER(intk), INTENT(in) :: igrid

#ifdef _MGLET_DEBUG_
        IF (myid /= idprocofgrd(igrid)) THEN
            CALL errr(__FILE__, __LINE__)
        END IF
        IF (iface < 1 .OR. iface > 6) THEN
            CALL errr(__FILE__, __LINE__)
        END IF
#endif

        ibb = ipbb(iface, igrid)
    END SUBROUTINE get_ibb


    SUBROUTINE get_len3(len, igrid)
    !$omp declare target

        INTEGER(intk), INTENT(in) :: igrid
        INTEGER(intk), INTENT(out) :: len

        INTEGER(intk) :: kk, jj, ii

        ! IF (myid /= idprocofgrd(igrid)) THEN
        !     CALL errr(__FILE__, __LINE__)
        ! END IF

        CALL get_mgdims(kk, jj, ii, igrid)
        len = kk*jj*ii
    END SUBROUTINE get_len3
END MODULE pointers_mod
