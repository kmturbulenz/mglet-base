MODULE noib_mod
    USE core_mod
    USE ibmodel_mod, ONLY: ibmodel_t
    USE noib_restrict_mod, ONLY: noib_restrict_t

    IMPLICIT NONE(type, external)
    PRIVATE

    TYPE, EXTENDS(ibmodel_t) :: noib_t
    CONTAINS
        PROCEDURE :: blockbp
        PROCEDURE :: read_stencils
        PROCEDURE, NOPASS :: giteig
        PROCEDURE :: divcal
        PROCEDURE, NOPASS :: divcal_grid
        FINAL :: destructor
    END TYPE noib_t

    INTERFACE noib_t
        MODULE PROCEDURE :: constructor
    END INTERFACE noib_t

    PUBLIC :: noib_t, constructor

CONTAINS
    FUNCTION constructor() RESULT(ib)
        ! Subroutine arguments
        CLASS(ibmodel_t), ALLOCATABLE :: ib

        ! Local variables
        REAL(realk), POINTER, CONTIGUOUS :: bp(:), sdiv(:)

        IF (myid == 0) THEN
            WRITE(*, '(A)') "Using 'noib' immersed boundary method"
            WRITE(*, '()')
        END IF

        ALLOCATE(noib_t :: ib)
        ALLOCATE(noib_restrict_t :: ib%restrict_op)

        ib%type = "NONE"

        ! "NONE" also has a BP field, always defined to 1.0. This saves a lot
        ! of conditional programming in places where performance is not
        ! critical, for instance initialization of pressure solver variables
        ! etc.
        CALL set_field("BP")
        CALL get_fieldptr(bp, "BP")
        bp = 1.0

        CALL set_field("SDIV")
        CALL get_fieldptr(sdiv, "SDIV")
        sdiv = 0.0
    END FUNCTION constructor


    SUBROUTINE destructor(this)
        TYPE(noib_t), INTENT(inout) :: this
        IF (ALLOCATED(this%restrict_op)) THEN
            DEALLOCATE(this%restrict_op)
        END IF
    END SUBROUTINE destructor


    SUBROUTINE blockbp(this)
        CLASS(noib_t), INTENT(inout) :: this

        ! Not doing anything
        CONTINUE
    END SUBROUTINE blockbp


    SUBROUTINE read_stencils(this)
        CLASS(noib_t), INTENT(inout) :: this

        ! Not doing anything
        CONTINUE
    END SUBROUTINE read_stencils


    SUBROUTINE giteig()
        ! Subroutine arguments
        ! none...

        ! Local variables
        INTEGER(intk) :: k, j, i
        INTEGER(intk) :: kk, jj, ii
        INTEGER(intk) :: igr, igrid
        INTEGER(intk) :: nfro, nbac, nrgt, nlft, nbot, ntop

        REAL(realk), POINTER, CONTIGUOUS :: dx(:), dy(:), dz(:)
        REAL(realk), POINTER, CONTIGUOUS :: ae(:),  aw(:),  an(:), as(:), &
            at(:), ab(:)
        REAL(realk), POINTER, CONTIGUOUS :: ap(:, :, :)
        REAL(realk), POINTER, CONTIGUOUS :: bp(:, :, :)

        ! 1-D fields used in the pressure solver
        CALL set_field("GSAW", ndim=1, get_len=get_ii)
        CALL set_field("GSAE", ndim=1, get_len=get_ii)
        CALL set_field("GSAS", ndim=1, get_len=get_jj)
        CALL set_field("GSAN", ndim=1, get_len=get_jj)
        CALL set_field("GSAB", ndim=1, get_len=get_kk)
        CALL set_field("GSAT", ndim=1, get_len=get_kk)

        ! 3-D fields in the pressure solver
        CALL set_field("GSAP")

        DO igr = 1, nmygrids
            igrid = mygrids(igr)

            CALL get_mgdims(kk, jj, ii, igrid)

            CALL get_fieldptr(dx, "DX", igrid)
            CALL get_fieldptr(dy, "DY", igrid)
            CALL get_fieldptr(dz, "DZ", igrid)
            CALL get_fieldptr(bp, "BP", igrid)

            CALL get_fieldptr(aw, "GSAW", igrid)
            CALL get_fieldptr(ae, "GSAE", igrid)
            CALL get_fieldptr(as, "GSAS", igrid)
            CALL get_fieldptr(an, "GSAN", igrid)
            CALL get_fieldptr(ab, "GSAB", igrid)
            CALL get_fieldptr(at, "GSAT", igrid)

            CALL get_fieldptr(ap, "GSAP", igrid)

            DO i = 3, ii-2
                ae(i) = 2.0/((dx(i-1)+dx(i))*dx(i))
                aw(i) = 2.0/((dx(i-1)+dx(i))*dx(i-1))
            END DO
            DO j = 3, jj-2
                an(j) = 2.0/((dy(j-1)+dy(j))*dy(j))
                as(j) = 2.0/((dy(j-1)+dy(j))*dy(j-1))
            END DO
            DO k = 3, kk-2
                at(k) = 2.0/((dz(k-1)+dz(k))*dz(k))
                ab(k) = 2.0/((dz(k-1)+dz(k))*dz(k-1))
            END DO

            DO i = 3, ii-2
                DO j = 3, jj-2
                    DO k = 3, kk-2
                        ap(k, j, i) = -2.0/(dx(i-1)*dx(i)) &
                            -2.0/(dy(j-1)*dy(j)) &
                            -2.0/(dz(k-1)*dz(k))
                    END DO
                END DO
            END DO

            DO i = 3, ii-2
                DO j = 3, jj-2
                    DO k = 3, kk-2
                        ap(k, j, i) = ap(k, j, i) &
                            + aw(i)*(1.0-bp(k, j, i-1)*bp(k, j, i)) &
                            + ae(i)*(1.0-bp(k, j, i  )*bp(k, j, i+1)) &
                            + as(j)*(1.0-bp(k, j-1, i)*bp(k, j  , i)) &
                            + an(j)*(1.0-bp(k, j  , i)*bp(k, j+1, i)) &
                            + ab(k)*(1.0-bp(k-1, j, i)*bp(k  , j, i)) &
                            + at(k)*(1.0-bp(k  , j, i)*bp(k+1, j, i))
                    END DO
                END DO
            END DO

            CALL get_mgbasb(nfro, nbac, nrgt, nlft, nbot, ntop, igrid)

            ! Front/West
            IF (nfro == 2 .OR. nfro == 5 .OR. nfro == 6 .OR. nfro == 19) THEN
                DO j = 3, jj-2
                    DO k = 3, kk-2
                        ap(k, j, 3) = ap(k, j, 3) &
                            + aw(3)*(bp(k, j, 2)*bp(k, j, 3))
                    END DO
                END  DO
                aw(3) = 0.0
            END IF

            ! Back/East
            IF (nbac == 2 .OR. nbac == 5 .OR. nbac == 6) THEN
                DO j = 3, jj-2
                    DO k = 3, kk-2
                        ap(k, j, ii-2) = ap(k, j, ii-2) &
                            + ae(ii-2)*(bp(k, j, ii-2)*bp(k, j, ii-1))
                    END DO
                END  DO
                ae(ii-2) = 0.0
            END IF

            ! Right/South
            IF (nrgt == 2 .OR. nrgt == 5 .OR. nrgt == 6 .OR. nrgt == 19) THEN
                DO i = 3, ii-2
                    DO k = 3, kk-2
                        ap(k, 3, i) = ap(k, 3, i) &
                            + as(3)*(bp(k, 2, i)*bp(k, 3, i))
                    END DO
                END  DO
                as(3) = 0.0
            END IF

            ! Left/North
            IF (nlft == 2 .OR. nlft == 5 .OR. nlft == 6) THEN
                DO i = 3, ii-2
                    DO k = 3, kk-2
                        ap(k, jj-2, i) = ap(k, jj-2, i) &
                            + an(jj-2)*(bp(k, jj-2, i)*bp(k, jj-1, i))
                    END DO
                END  DO
                an(jj-2) = 0.0
            END IF

            ! Bottom
            IF (nbot == 2 .OR. nbot == 5 .OR. nbot == 6 .OR. nbot == 19) THEN
                DO i = 3, ii-2
                    DO j = 3, jj-2
                        ap(3, j, i) = ap(3, j, i) &
                            + ab(3)*(bp(2, j, i)*bp(3, j, i))
                    END DO
                END  DO
                ab(3) = 0.0
            END IF

            ! Top
            IF (ntop == 2 .OR. ntop == 5 .OR. ntop == 6) THEN
                DO i = 3, ii-2
                    DO j = 3, jj-2
                        ap(kk-2, j, i) = ap(kk-2, j, i) &
                            + at(kk-2)*(bp(kk-2, j, i)*bp(kk-1, j, i))
                    END DO
                END  DO
                at(kk-2) = 0.0
            END IF
        END DO
    END SUBROUTINE giteig


    SUBROUTINE divcal(this, div, u, v, w, fak, ctyp)
        ! Subroutine arguments
        CLASS(noib_t), INTENT(inout) :: this
        TYPE(field_t), INTENT(inout) :: div
        TYPE(field_t), INTENT(in) :: u
        TYPE(field_t), INTENT(in) :: v
        TYPE(field_t), INTENT(in) :: w
        REAL(realk), INTENT(in) :: fak
        CHARACTER(len=1), INTENT(in), OPTIONAL :: ctyp

        ! Local variables
        INTEGER(intk) :: i, igrid, ip3
        INTEGER(intk) :: kk, jj, ii
        TYPE(field_t), POINTER :: rddx_f, rddy_f, rddz_f
        REAL(realk), CONTIGUOUS, POINTER :: rddx(:), rddy(:), rddz(:)

        CALL get_field(rddx_f, "RDDX")
        CALL get_field(rddy_f, "RDDY")
        CALL get_field(rddz_f, "RDDZ")

        DO i = 1, nmygrids
            igrid = mygrids(i)
            CALL get_mgdims(kk, jj, ii, igrid)
            CALL get_ip3(ip3, igrid)
            CALL rddx_f%get_ptr(rddx, igrid)
            CALL rddy_f%get_ptr(rddy, igrid)
            CALL rddz_f%get_ptr(rddz, igrid)

            CALL divcal_grid(kk, jj, ii, fak, div%arr(ip3), u%arr(ip3), &
                v%arr(ip3), w%arr(ip3), rddx, rddy, rddz)
        END DO
    END SUBROUTINE divcal


    PURE SUBROUTINE divcal_grid(kk, jj, ii, fak, div, u, v, w, rddx, rddy, &
            rddz, bp, sdiv)
        ! Subroutine arguments
        INTEGER(intk), INTENT(in) :: kk, jj, ii
        REAL(realk), INTENT(in) :: fak
        REAL(realk), INTENT(inout) :: div(kk, jj, ii)
        REAL(realk), INTENT(in) :: u(kk, jj, ii)
        REAL(realk), INTENT(in) :: v(kk, jj, ii)
        REAL(realk), INTENT(in) :: w(kk, jj, ii)
        REAL(realk), INTENT(in) :: rddx(ii)
        REAL(realk), INTENT(in) :: rddy(jj)
        REAL(realk), INTENT(in) :: rddz(kk)
        REAL(realk), INTENT(in), OPTIONAL :: bp(kk, jj, ii)
        REAL(realk), INTENT(in), OPTIONAL :: sdiv(kk, jj, ii)

        ! Local variables
        INTEGER(intk) :: k, j, i

        DO i = 3, ii-2
            DO j = 3, jj-2
                DO k = 3, kk-2
                    div(k, j, i) = fak*((u(k, j, i) - u(k, j, i-1))*rddx(i) &
                        + (v(k, j, i) - v(k, j-1, i))*rddy(j) &
                        + (w(k, j, i) - w(k-1, j, i) )*rddz(k))
                END DO
            END DO
        END DO

        ! TODO: If SDIV is properly masked in the ghost layers, the indices
        ! here could be from 1 to ii etc. That will lead to ever so slightly
        ! better performance. I leave it as is until SDIV is properly checked.
        IF (PRESENT(sdiv)) THEN
            DO i = 3, ii-2
                DO j = 3, jj-2
                    DO k = 3, kk-2
                        div(k, j, i) = div(k, j, i) + fak*sdiv(k, j, i)
                    END DO
                END DO
            END DO
        END IF

        IF (PRESENT(bp)) THEN
            DO i = 1, ii
                DO j = 1, jj
                    DO k = 1, kk
                        div(k, j, i) = bp(k, j, i)*div(k, j, i)
                    END DO
                END DO
            END DO
        END IF
    END SUBROUTINE divcal_grid
END MODULE noib_mod
