MODULE gc_createstencils_mod
    USE MPI_f08, ONLY: MPI_Wtime

    USE core_mod, ONLY: intk, realk, errr

    ! Commom routines shared among flow, scalar etc. stencil creation
    ! TODO: Fix subroutine names in this file...

    IMPLICIT NONE (type, external)
    PRIVATE

    PUBLIC :: wmcheckneighbor, choosestencil, wmmultimatrix, wmaddcoefflist, &
        wmindexlistn

CONTAINS

    PURE SUBROUTINE wmcheckneighbor(k, j, i, kk, jj, ii, bp, bconds, ncells, &
            inlst, jnlst, knlst, istc, jstc, kstc, nstc)

        ! Subroutine arguments
        INTEGER(intk), INTENT(in) :: k, j, i
        INTEGER(intk), INTENT(in) :: kk, jj, ii
        REAL(realk), INTENT(in) :: bp(kk, jj, ii)
        INTEGER(intk), INTENT(in) :: bconds(6)
        INTEGER(intk), INTENT(in) :: ncells

        INTEGER(intk), CONTIGUOUS, INTENT(in) :: inlst(:), jnlst(:), knlst(:)
        INTEGER(intk), CONTIGUOUS, INTENT(inout) :: istc(:), jstc(:), kstc(:)
        INTEGER(intk), INTENT(out) :: nstc

        ! Local variables
        INTEGER(intk) :: irfro, irbac, irrgt, irlft, irbot, irtop
        INTEGER(intk) :: n, noneigh, it, jt, kt

        nstc = 0

        irfro = 2
        irbac = 2
        irrgt = 2
        irlft = 2
        irbot = 2
        irtop = 2

        ! OP1
        IF (bconds(1) == 3) irfro = 0
        IF (bconds(2) == 3) irbac = 0
        IF (bconds(3) == 3) irrgt = 0
        IF (bconds(4) == 3) irlft = 0
        IF (bconds(5) == 3) irbot = 0
        IF (bconds(6) == 3) irtop = 0

        ! OP2
        IF (bconds(1) == 4) irfro = 0
        IF (bconds(2) == 4) irbac = 0
        IF (bconds(3) == 4) irrgt = 0
        IF (bconds(4) == 4) irlft = 0
        IF (bconds(5) == 4) irbot = 0
        IF (bconds(6) == 4) irtop = 0

        ! PAR
        IF (bconds(1) == 8) irfro = 1
        IF (bconds(2) == 8) irbac = 1
        IF (bconds(3) == 8) irrgt = 1
        IF (bconds(4) == 8) irlft = 1
        IF (bconds(5) == 8) irbot = 1
        IF (bconds(6) == 8) irtop = 1

        ! CON
        IF (bconds(1) == 7) irfro = 1
        IF (bconds(2) == 7) irbac = 1
        IF (bconds(3) == 7) irrgt = 1
        IF (bconds(4) == 7) irlft = 1
        IF (bconds(5) == 7) irbot = 1
        IF (bconds(6) == 7) irtop = 1

        ! SLI
        IF (bconds(1) == 6) irfro = 0
        IF (bconds(2) == 6) irbac = 0
        IF (bconds(3) == 6) irrgt = 0
        IF (bconds(4) == 6) irlft = 0
        IF (bconds(5) == 6) irbot = 0
        IF (bconds(6) == 6) irtop = 0


        DO n = 1, ncells
            noneigh = 0

            it = i + inlst(n)
            IF (it < 3 - irfro) noneigh = 1
            IF (it > ii-2+irbac) noneigh = 1

            jt = j + jnlst(n)
            IF (jt < 3 - irrgt) noneigh = 1
            IF (jt > jj-2+irlft) noneigh = 1

            kt = k + knlst(n)
            IF (kt < 3 - irbot) noneigh = 1
            IF (kt > kk-2+irtop) noneigh = 1

            IF (noneigh == 0) THEN
                IF (NINT(bp(kt, jt, it)) == 1) THEN
                    nstc = nstc + 1
                    istc(nstc) = it
                    jstc(nstc) = jt
                    kstc(nstc) = kt
                END IF
            END IF
        END DO
    END SUBROUTINE wmcheckneighbor


    SUBROUTINE choosestencil(k, j, i, nstc, istc, jstc, kstc, velpts)
        ! Subroutine arguments
        INTEGER(intk), INTENT(in) :: k, j, i, nstc
        INTEGER(intk), CONTIGUOUS, INTENT(inout) :: istc(:), jstc(:), kstc(:)
        INTEGER(intk), INTENT(out) :: velpts

        ! Local variables
        INTEGER(intk) :: n, id, jd, kd, idges, jdges, kdges
        INTEGER(intk) :: istctmp(SIZE(istc)), jstctmp(SIZE(jstc)), &
            kstctmp(SIZE(kstc))
        INTEGER, PARAMETER :: output = 0

        istctmp = 0
        jstctmp = 0
        kstctmp = 0

        velpts = 0
        DO n = 1, nstc
            id = ABS(istc(n) - i)
            jd = ABS(jstc(n) - j)
            kd = ABS(kstc(n) - k)

            IF (id + jd + kd == 1) THEN
                velpts = velpts + 1
                istctmp(velpts) = istc(n)
                jstctmp(velpts) = jstc(n)
                kstctmp(velpts) = kstc(n)
            END IF
        END DO

        IF (velpts == 1) THEN
            idges = (istctmp(1) - i)*2 + i
            jdges = (jstctmp(1) - j)*2 + j
            kdges = (kstctmp(1) - k)*2 + k

            DO n = 1, nstc
                IF (idges == istc(n) .AND. jdges == jstc(n) .AND. &
                        kdges == kstc(n)) THEN
                    velpts = velpts + 1
                    istctmp(velpts) = istc(n)
                    jstctmp(velpts) = jstc(n)
                    kstctmp(velpts) = kstc(n)
                END IF
            END DO
        END IF

        IF (velpts == 1) THEN
            IF (output == 1) THEN
                WRITE(*, *) 'choosestencil i, j, k', i, j, k
                CALL errr(__FILE__, __LINE__)
            END IF
            velpts = 0
        END IF

        IF (output == 1) THEN
            IF (velpts == 0) THEN
                WRITE(*, *) 'choosestencil i, j, k', i, j, k
                CALL errr(__FILE__, __LINE__)
            END IF
        END IF

        istc = istctmp
        jstc = jstctmp
        kstc = kstctmp
    END SUBROUTINE choosestencil


    SUBROUTINE wmmultimatrix(compon, k, j, i, kt, jt, it, kk, jj, ii, &
            x, y, z, xstag, ystag, zstag, bzelltyp, istc, jstc, kstc, &
            velpts, nvecs, ucell, coeffstc, coeffstcvel, &
            acoeffstc, acoeffstcvel, found)

        ! Subroutine arguments
        INTEGER(intk), INTENT(in) :: compon
        INTEGER(intk), INTENT(in) :: k, j, i
        INTEGER(intk), INTENT(in) :: kt, jt, it
        INTEGER(intk), INTENT(in) :: kk, jj, ii
        REAL(realk), INTENT(in) :: x(ii), y(jj), z(kk)
        REAL(realk), INTENT(in) :: xstag(ii), ystag(jj), zstag(kk)
        INTEGER(intk), INTENT(in) :: bzelltyp(kk, jj, ii)
        INTEGER(intk), CONTIGUOUS, INTENT(in) :: istc(:), jstc(:), kstc(:)
        INTEGER(intk), INTENT(in) :: velpts
        REAL(realk), CONTIGUOUS, INTENT(in) :: nvecs(:, :), ucell(:, :)
        REAL(realk), CONTIGUOUS, INTENT(out) :: coeffstc(:), coeffstcvel(:)
        REAL(realk), INTENT(out) :: acoeffstc, acoeffstcvel
        INTEGER(intk), INTENT(out) :: found

        ! Local variables
        INTEGER(intk) :: counter
        REAL(realk) :: a, b, c, dpl, velu, velv, velw

        ! Parameter for the accuracy of the numerical integration
        INTEGER(intk), PARAMETER :: parts = 20

        counter = -bzelltyp(kt, jt, it)
        a = nvecs(1, counter)
        b = nvecs(2, counter)
        c = nvecs(3, counter)
        dpl = nvecs(4, counter)
        velu = ucell(1, counter)
        velv = ucell(2, counter)
        velw = ucell(3, counter)

        coeffstc = 0.0
        acoeffstc = 0.0

        coeffstcvel = 0.0
        acoeffstcvel = 0.0

        ! Warning: Calling calcfluxu in three different ways depending on the
        ! direction, changing order of arguments. Be very careful when/if
        ! changing interface!
        SELECT CASE (compon)
        CASE (1)
            CALL calcfluxu(k, j, i, kk, jj, ii, x, ystag, zstag, xstag, y, z, &
                a, b, c, dpl, velu, velpts, parts, istc, jstc, kstc, &
                coeffstc, coeffstcvel, acoeffstc, acoeffstcvel, found)
        CASE (2)
            CALL calcfluxu(k, i, j, kk, ii, jj, y, xstag, zstag, ystag, x, z, &
                b, a, c, dpl, velv, velpts, parts, jstc, istc, kstc, &
                coeffstc, coeffstcvel, acoeffstc, acoeffstcvel, found)
        CASE (3)
            CALL calcfluxu(i, j, k, ii, jj, kk, z, ystag, xstag, zstag, y, x, &
                c, b, a, dpl, velw, velpts, parts, kstc, jstc, istc, &
                coeffstc, coeffstcvel, acoeffstc, acoeffstcvel, found)
        CASE DEFAULT
            CALL errr(__FILE__, __LINE__)
        END SELECT
    END SUBROUTINE wmmultimatrix


    SUBROUTINE calcfluxu(k, j, i, kk, jj, ii, x, y, z, xstag, ystag, zstag, &
            a, b, c, dpl, avel, velpts, parts, istc, jstc, kstc, &
            coeffstc, coeffstcvel, acoeffstc, acoeffstcvel, found)

        ! Subroutine arguments
        ! NB: Strange call of this subroutine from wmmultimatrix, variables
        ! (k, j, i, xstag, x etc) are deliberately interchanged to compute for
        ! different directiosn
        INTEGER(intk), INTENT(in) :: k, j, i
        INTEGER(intk), INTENT(in) :: kk, jj, ii
        REAL(realk), INTENT(in) :: x(ii), y(jj), z(kk)
        REAL(realk), INTENT(in) :: xstag(ii), ystag(jj), zstag(kk)
        ! Components of nvec 1..4
        REAL(realk), INTENT(in) :: a, b, c, dpl
        REAL(realk), INTENT(in) :: avel
        INTEGER(intk), INTENT(in) :: velpts
        INTEGER(intk), INTENT(in) :: parts
        INTEGER(intk), CONTIGUOUS, INTENT(in) :: istc(:), jstc(:), kstc(:)
        REAL(realk), CONTIGUOUS, INTENT(inout) :: coeffstc(:), coeffstcvel(:)
        REAL(realk), INTENT(inout) :: acoeffstc, acoeffstcvel
        INTEGER(intk), INTENT(out) :: found

        ! Local variables
        INTEGER(intk), PARAMETER :: output = 0
        INTEGER(intk) :: it, jt, kt
        INTEGER(intk) :: n
        INTEGER(intk) :: carea, carea2
        REAL(realk) :: dx, dy, dz
        REAL(realk) :: xs, ys, zs
        REAL(realk) :: ddx, ddy, ddz, xm, ym, zm, dt, area, dir
        REAL(realk) :: dist(velpts), ce(velpts), cesum, dall, dpoint
        REAL(realk) :: sum

        found = 0

        ! Abstand fuer Mitte zu interpolierender Flaeche
        dpoint = a*xstag(i) + b*ystag(j) + c*zstag(k) - dpl

        dall = 0.0
        dist = 0.0
        DO n = 1, velpts
            it = istc(n)
            jt = jstc(n)
            kt = kstc(n)
            dist(n) = a*xstag(it) + b*ystag(jt) + c*zstag(kt) - dpl

            ! Stencils auf neg. Seite der Ebene tragen nichts bei
            !     dist(n) = max(dist(n), 0.0)
            ! und Stencils, die naeher an der Wand sind als dpoint
            !    if (dist(n) < abs(dpoint)) then
            !    JK 23.2.2009: um 10% naeher ist erlaubt
            IF (dist(n) < abs(dpoint)*0.9) THEN
                dist(n) = 0.0
            END IF

            dall = dall + dist(n)**2
        END DO

        ! Wenn a: Normalenvektor = 0
        ! oder b: alle Stencil-Punkte auf der falschen Seite der Ebene
        ! erzeugen keinen stencil
        IF (ABS(dall) <= TINY(dall)) THEN
            IF (output == 1) THEN
                CALL errr(__FILE__, __LINE__)
            END IF
            ce = 0.0
            RETURN
        END IF

        cesum = 0.0
        DO n = 1, velpts
            ce(n) = dist(n)/dall
            cesum = cesum + ce(n)
        END DO

        dx = x(i) - x(i-1)
        dy = y(j) - y(j-1)
        dz = z(k) - z(k-1)

        ddx = dx/parts
        ddy = dy/parts
        ddz = dz/parts

        xs = x(i-1) - ddx/2.0
        ys = y(j-1) - ddy/2.0
        zs = z(k-1) - ddz/2.0

        carea = 0
        carea2 = 0
        area = 0.0
        xm = xstag(i)
        dt = -dpl + a*xm

        DO jt = 1, parts
            ym = ys + ddy*REAL(jt)
            DO kt = 1, parts
                zm = zs + ddz*REAL(kt)
                dir = dt + b*ym + c*zm
                IF (dir > 0.0) THEN
                    dpoint =  a*xm + b*ym + c*zm - dpl
                    carea = carea + 1
                    DO n = 1, velpts
                        coeffstc(n) = coeffstc(n) + dpoint*ce(n)
                    END DO
                    acoeffstc = acoeffstc + 1.0-dpoint*cesum
                ELSE
                    carea2 = carea2 + 1
                END IF
            END DO
        END DO

        area = ddy*ddz

        ! Flux stencil
        DO n = 1, velpts
            coeffstc(n) = coeffstc(n)*area/(dy*dz)
        END DO

        ! Additional flux stencil 1: in the fluid
        acoeffstc = acoeffstc*avel*area/(dy*dz)
        ! Additional flux stencil 2: in the body
        acoeffstc = acoeffstc + avel*(carea2*area)/(dy*dz)

        ! Velocity stencil... ATTENTION XSTAG or X?
        dpoint =  a*xstag(i) + b*ystag(j) + c*zstag(k) - dpl
        DO n = 1, velpts
            coeffstcvel(n) = dpoint*ce(n)
        END DO
        acoeffstcvel = avel * (1.0 - dpoint*cesum)

        ! Stencil erzeugt
        found = 1

        sum = 0.0
        DO n = 1, velpts
            sum = sum + coeffstc(n)
        END DO

        ! jk 7.2.2008
        ! Grenze auf 3/5: Wert fuer Wuerfelgitter, wenn Ebene durch Ecke geht
        IF (sum > 0.61) THEN
            found = -1
            IF (output == 1) THEN
                CALL errr(__FILE__, __LINE__)
            END IF
        END IF
    END SUBROUTINE calcfluxu


    PURE SUBROUTINE wmaddcoefflist(velpts, xpolr, xpolrvel, &
            pntxpolr, acoeffstc, acoeffstcvel, coeffstc, coeffstcvel)

        ! Subroutine arguments
        INTEGER(intk), INTENT(in) :: velpts
        REAL(realk), CONTIGUOUS, INTENT(inout) :: xpolr(:)
        REAL(realk), CONTIGUOUS, INTENT(inout) :: xpolrvel(:)
        INTEGER(intk), INTENT(inout) :: pntxpolr
        REAL(realk), INTENT(in) :: acoeffstc, acoeffstcvel
        REAL(realk), CONTIGUOUS, INTENT(in) :: coeffstc(:)
        REAL(realk), CONTIGUOUS, INTENT(in) :: coeffstcvel(:)

        ! Local variables
        INTEGER(intk) :: n

        pntxpolr = pntxpolr - (velpts+1)

        DO n = 1, velpts
            pntxpolr = pntxpolr + 1
            xpolr(pntxpolr) = 0.5*xpolr(pntxpolr) + 0.5*coeffstc(n)
            xpolrvel(pntxpolr) = 0.5*xpolrvel(pntxpolr) + 0.5*coeffstcvel(n)
        END DO

        pntxpolr = pntxpolr + 1
        xpolr(pntxpolr) = 0.5*xpolr(pntxpolr) + 0.5*acoeffstc
        xpolrvel(pntxpolr) = 0.5*xpolrvel(pntxpolr) + 0.5*acoeffstcvel
    END SUBROUTINE wmaddcoefflist


    PURE SUBROUTINE wmindexlistn(inlst, jnlst, knlst, nnlst)
        ! Subroutine arguments
        INTEGER(intk), CONTIGUOUS, INTENT(out) :: inlst(:), jnlst(:), knlst(:)
        INTEGER(intk), INTENT(out) :: nnlst

        ! Local variables
        INTEGER(intk) :: k, j, i, kk, jj, ii, n
        INTEGER(intk), PARAMETER :: boxsize = 7
        INTEGER(intk), PARAMETER :: radius = (boxsize-1)/2

        ! Initialize INTENT(out)
        inlst = 0
        jnlst = 0
        knlst = 0

        n = 0
        DO i = 1, boxsize
            DO j = 1, boxsize
                DO k = 1, boxsize
                    ii = i - (INT(boxsize/2.0)+1)
                    jj = j - (INT(boxsize/2.0)+1)
                    kk = k - (INT(boxsize/2.0)+1)

                    IF ((ii /= 0) .OR. (jj /= 0) .OR. (kk /= 0)) THEN
                        IF ((ABS(ii)**2 + ABS(jj)**2 + ABS(kk)**2)**0.5 <= &
                                radius) THEN
                            n = n + 1
                            inlst(n) = ii
                            jnlst(n) = jj
                            knlst(n) = kk
                        END IF
                    END IF
                END DO
            END DO
        END DO

        nnlst = n
    END SUBROUTINE wmindexlistn

END MODULE gc_createstencils_mod
