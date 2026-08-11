MODULE noib_mod
    USE core_mod
    USE ibmodel_mod, ONLY: ibmodel_t
    USE fieldhelper_mod, ONLY: map_arr_to_device

    IMPLICIT NONE(type, external)
    PRIVATE

    TYPE, EXTENDS(ibmodel_t) :: noib_t
    CONTAINS
        PROCEDURE :: blockbp
        PROCEDURE :: read_stencils
        PROCEDURE, NOPASS :: giteig
        FINAL :: destructor
    END TYPE noib_t

    PUBLIC :: noib_t, constructor

CONTAINS
    SUBROUTINE constructor(ib)
        ! Subroutine arguments
        CLASS(ibmodel_t), ALLOCATABLE, INTENT(out) :: ib

        ! Local variables
        TYPE(field_t), POINTER :: bp, sdiv

        IF (myid == 0) THEN
            WRITE(*, '(A)') "Using 'noib' immersed boundary method"
            WRITE(*, '()')
        END IF

        ALLOCATE(noib_t :: ib)

        ib%type = "NONE"

        ! "NONE" also has a BP field, always defined to 1.0. This saves a lot
        ! of conditional programming in places where performance is not
        ! critical, for instance initialization of pressure solver variables
        ! etc.
        CALL set_field("BP")
        CALL get_field(bp, "BP")
        bp%arr = 1.0
        CALL map_arr_to_device(bp, message="to:bp%arr")

        CALL set_field("SDIV")
        CALL get_field(sdiv, "SDIV")
        sdiv%arr = 0.0
        CALL map_arr_to_device(sdiv, message="to:sdiv%arr")
    END SUBROUTINE constructor


    SUBROUTINE destructor(this)
        TYPE(noib_t), INTENT(inout) :: this

        ! Nothing to do
    END SUBROUTINE destructor


    SUBROUTINE blockbp(this, stop_now)
        CLASS(noib_t), INTENT(inout) :: this
        LOGICAL, INTENT(out) :: stop_now

        stop_now = .FALSE.
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
        TYPE(field_t), POINTER :: gsaw, gsae, gsas, gsan, gsab, gsat, gsap
        TYPE(field_t), POINTER :: dx_f, dy_f, dz_f, bp_f

        ! 1-D fields used in the pressure solver
        CALL set_field("GSAW", ndim=1, get_len=get_ii)
        CALL set_field("GSAE", ndim=1, get_len=get_ii)
        CALL set_field("GSAS", ndim=1, get_len=get_jj)
        CALL set_field("GSAN", ndim=1, get_len=get_jj)
        CALL set_field("GSAB", ndim=1, get_len=get_kk)
        CALL set_field("GSAT", ndim=1, get_len=get_kk)
        CALL get_field(gsaw, "GSAW")
        CALL get_field(gsae, "GSAE")
        CALL get_field(gsas, "GSAS")
        CALL get_field(gsan, "GSAN")
        CALL get_field(gsab, "GSAB")
        CALL get_field(gsat, "GSAT")

        ! 3-D fields in the pressure solver
        CALL set_field("GSAP")
        CALL get_field(gsap, "GSAP")

        CALL get_field(dx_f, "DX")
        CALL get_field(dy_f, "DY")
        CALL get_field(dz_f, "DZ")
        CALL get_field(bp_f, "BP")

        DO igr = 1, nmygrids
            igrid = mygrids(igr)

            CALL get_mgdims(kk, jj, ii, igrid)

            CALL dx_f%get_ptr(dx, igrid)
            CALL dy_f%get_ptr(dy, igrid)
            CALL dz_f%get_ptr(dz, igrid)
            CALL bp_f%get_ptr(bp, igrid)

            CALL gsaw%get_ptr(aw, igrid)
            CALL gsae%get_ptr(ae, igrid)
            CALL gsas%get_ptr(as, igrid)
            CALL gsan%get_ptr(an, igrid)
            CALL gsab%get_ptr(ab, igrid)
            CALL gsat%get_ptr(at, igrid)

            CALL gsap%get_ptr(ap, igrid)

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
                            + ae(i)*(1.0-bp(k, j, i)*bp(k, j, i+1)) &
                            + as(j)*(1.0-bp(k, j-1, i)*bp(k, j, i)) &
                            + an(j)*(1.0-bp(k, j, i)*bp(k, j+1, i)) &
                            + ab(k)*(1.0-bp(k-1, j, i)*bp(k, j, i)) &
                            + at(k)*(1.0-bp(k, j, i)*bp(k+1, j, i))
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

        CALL map_arr_to_device(gsaw, gsae, gsas, gsan, gsab, gsat, &
            message="to:gs-neighbors")
        CALL map_arr_to_device(gsap, message="to:gs-ap")
    END SUBROUTINE giteig
END MODULE noib_mod
