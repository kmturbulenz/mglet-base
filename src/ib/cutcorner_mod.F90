MODULE cutcorner_mod
    USE core_mod, ONLY: realk, intk, sortrx, field_t, mygridslvl, nmygridslvl, &
        minlevel, maxlevel, get_fieldptr, get_mgdims, get_ip3, get_ip3n
    USE ibconst_mod, ONLY: maccur
    USE topol_mod, ONLY: topol_t

    IMPLICIT NONE (type, external)
    PRIVATE

    PUBLIC :: cutcorner, cutcorner_level, cutcorner_grid, calcint

CONTAINS
    SUBROUTINE cutcorner(topol, ntrimax, kanteu, kantev, &
        kantew, triau, triav, triaw)

        ! Subroutine arguments
        TYPE(topol_t), INTENT(in) :: topol
        INTEGER(intk), INTENT(in) :: ntrimax
        TYPE(field_t), INTENT(inout) :: kanteu
        TYPE(field_t), INTENT(inout) :: kantev
        TYPE(field_t), INTENT(inout) :: kantew
        INTEGER(intk), INTENT(out) :: triau(*)
        INTEGER(intk), INTENT(out) :: triav(*)
        INTEGER(intk), INTENT(out) :: triaw(*)

        ! Local variables
        INTEGER(intk) :: ilevel

        DO ilevel = minlevel, maxlevel
            CALL cutcorner_level(ilevel, topol, ntrimax, kanteu, kantev, &
                kantew, triau, triav, triaw)
        END DO
    END SUBROUTINE cutcorner


    SUBROUTINE cutcorner_level(ilevel, topol, ntrimax, kanteu, kantev, &
            kantew, triau, triav, triaw)
        ! Calls cutcorner on all grids on the same level

        ! Subroutine arguments
        INTEGER(intk), INTENT(in) :: ilevel
        TYPE(topol_t), INTENT(in) :: topol
        INTEGER(intk), INTENT(in) :: ntrimax
        TYPE(field_t), INTENT(inout) :: kanteu
        TYPE(field_t), INTENT(inout) :: kantev
        TYPE(field_t), INTENT(inout) :: kantew
        INTEGER(intk), INTENT(out) :: triau(*)
        INTEGER(intk), INTENT(out) :: triav(*)
        INTEGER(intk), INTENT(out) :: triaw(*)

        ! Local variables
        INTEGER(intk) :: i, igrid, kk, jj, ii, ip3, ip3n
        REAL(realk), POINTER, CONTIGUOUS :: xstag(:), ystag(:), zstag(:), &
            ddx(:), ddy(:), ddz(:)

        DO i = 1, nmygridslvl(ilevel)
            igrid = mygridslvl(i, ilevel)

            CALL get_fieldptr(xstag, "XSTAG", igrid)
            CALL get_fieldptr(ystag, "YSTAG", igrid)
            CALL get_fieldptr(zstag, "ZSTAG", igrid)

            CALL get_fieldptr(ddx, "DDX", igrid)
            CALL get_fieldptr(ddy, "DDY", igrid)
            CALL get_fieldptr(ddz, "DDZ", igrid)

            CALL get_mgdims(kk, jj, ii, igrid)
            CALL get_ip3(ip3, igrid)
            CALL get_ip3n(ip3n, ntrimax, igrid)

            CALL cutcorner_grid(kk, jj, ii, topol%n, ntrimax, xstag, &
                ystag, zstag, ddx, ddy, ddz, topol%topol, kanteu%arr(ip3), &
                kantev%arr(ip3), kantew%arr(ip3), triau(ip3n), triav(ip3n), &
                triaw(ip3n))
        END DO
    END SUBROUTINE cutcorner_level


    SUBROUTINE cutcorner_grid(kk, jj, ii, ntopol, ntrimax, xstag, ystag, &
        zstag, ddx, ddy, ddz, topol, kanteu, kantev, kantew, triau, triav, &
        triaw, exactedge)

        ! Subroutine arguments
        INTEGER(intk), INTENT(in) :: kk, jj, ii, ntopol, ntrimax
        REAL(realk), INTENT(in) :: xstag(ii), ystag(jj), zstag(kk), &
            ddx(ii), ddy(jj), ddz(kk), topol(3, 3, ntopol)

        REAL(realk), INTENT(out) :: kanteu(kk, jj, ii), kantev(kk, jj, ii), &
            kantew(kk, jj, ii)
        INTEGER(intk), INTENT(out) :: triau(ntrimax, kk, jj, ii), &
            triav(ntrimax, kk, jj, ii), triaw(ntrimax, kk, jj, ii)
        LOGICAL, OPTIONAL, INTENT(in) :: exactedge

        ! Local variables
        INTEGER(intk) :: found, foundx, foundy, foundz
        INTEGER(intk) :: i, j, k, itri
        INTEGER(intk) :: ntrilim(3, 2)
        REAL(realk) :: xp, yp, zp, a, b, c, x1, y1, z1
        REAL(realk) :: x2, y2, z2, x3, y3, z3, betrag
        REAL(realk) :: length, eps
        REAL(realk) :: xmin, xmax, ymin, ymax, zmin, zmax
        LOGICAL :: exactedge2

        ! Initializing variables
        kanteu = 0.0
        kantev = 0.0
        kantew = 0.0

        triau = 0
        triav = 0
        triaw = 0

        IF (PRESENT(exactedge)) THEN
            exactedge2 = exactedge
        ELSE
            exactedge2 = .FALSE.
        END IF

        ! Loop over all triangles and compute the intersections with the grid
        triangles: DO itri = 1, ntopol

            ! Triangle corners
            x1 = topol(1, 1, itri)
            x2 = topol(2, 1, itri)
            x3 = topol(3, 1, itri)
            y1 = topol(1, 2, itri)
            y2 = topol(2, 2, itri)
            y3 = topol(3, 2, itri)
            z1 = topol(1, 3, itri)
            z2 = topol(2, 3, itri)
            z3 = topol(3, 3, itri)

            ! Triangle bounding box
            xmin = MIN(x1, x2, x3)
            xmax = MAX(x1, x2, x3)
            ymin = MIN(y1, y2, y3)
            ymax = MAX(y1, y2, y3)
            zmin = MIN(z1, z2, z3)
            zmax = MAX(z1, z2, z3)

            ! If triangle is outisde grid, skip triangle
            IF (xmax < xstag(1) - (2*maccur+1)*ddx(1)) CYCLE
            IF (xmin > xstag(ii) + (2*maccur)*ddx(ii)) CYCLE
            IF (ymax < ystag(1) - (2*maccur+1)*ddy(1)) CYCLE
            IF (ymin > ystag(jj) + (2*maccur)*ddy(jj)) CYCLE
            IF (zmax < zstag(1) - (2*maccur+1)*ddz(1)) CYCLE
            IF (zmin > zstag(kk) + (2*maccur)*ddz(kk)) CYCLE

            ! Triangle normal vector
            a = (y2-y1)*(z3-z1) - (z2-z1)*(y3-y1)
            b = (z2-z1)*(x3-x1) - (x2-x1)*(z3-z1)
            c = (x2-x1)*(y3-y1) - (y2-y1)*(x3-x1)

            ! If triangle has zero area, skip triangle
            betrag = sqrt(a**2 + b**2 + c**2)
            IF (betrag < TINY(1.0_realk)) CYCLE
            a = a/betrag
            b = b/betrag
            c = c/betrag

            ! Ausdehnung des Dreiecks im kartesischen Gitter
            ! -->  Indices, berechnet aus Koordinaten
            DO i = 1, 3
                ntrilim(i, 1) = 999999
                ntrilim(i, 2) = -1
            END DO

            DO i = 1, ii
                IF (xmin <= xstag(i) + 2*maccur*ddx(i)) THEN
                    ntrilim(1, 1) = i
                    EXIT
                END IF
            END DO
            DO i = ii, 1, -1
                IF (xmax >= xstag(i) - 2*maccur*ddx(i)) THEN
                    ntrilim(1, 2) = i+1
                    EXIT
                END IF
            END DO
            DO j = 1, jj
                IF (ymin <= ystag(j) + 2*maccur*ddy(j)) THEN
                    ntrilim(2, 1) = j
                    EXIT
                END IF
            END DO
            DO j = jj, 1, -1
                IF (ymax >= ystag(j) - 2*maccur*ddy(j)) THEN
                    ntrilim(2, 2) = j+1
                    EXIT
                END IF
            END DO
            DO k = 1, kk
                IF (zmin <= zstag(k) + 2*maccur*ddz(k)) THEN
                    ntrilim(3, 1) = k
                    EXIT
                END IF
            END DO
            DO k = kk, 1, -1
                IF (zmax >= zstag(k) - 2*maccur*ddz(k)) THEN
                    ntrilim(3, 2) = k+1
                    EXIT
                END IF
            END DO

            ntrilim(1, 2) = MIN(ntrilim(1, 2), ii)
            ntrilim(2, 2) = MIN(ntrilim(2, 2), jj)
            ntrilim(3, 2) = MIN(ntrilim(3, 2), kk)

            ! Zellen suchen, die das Dreieck schneiden
            ! Innerhalb der zuvor bestimmten Index-Grenzen
            ! wird nach Schnittpunkten gesucht
            DO i = ntrilim(1,1), ntrilim(1,2)
                DO j = ntrilim(2,1), ntrilim(2,2)
                    DO k = ntrilim(3,1), ntrilim(3,2)

                        ! The length is the middle length when the cell has
                        ! different edge lengths.
                        length = (ddx(i) + ddy(j) + ddz(k)) &
                            - MIN(ddx(i), ddy(j), ddz(k)) &
                            - MAX(ddx(i), ddy(j), ddz(k))
                        eps = maccur*length

                        ! TODO: Cleanup usage of xp, yp, zp - remove??
                        xp = xstag(i)
                        yp = ystag(j)
                        zp = zstag(k)

                        CALL intinface(kk, jj, ii, xstag, ystag, zstag, &
                            ddx, ddy, ddz, &
                            i, j, k, a, b, c, xp, yp, zp, &
                            x1, y1, z1, x2, y2, z2, x3, y3, z3, &
                            eps, xmin, xmax, ymin, ymax, zmin, zmax, &
                            found, foundx, foundy, foundz, .FALSE.)

                        IF (foundx == 1) THEN
                            CALL insertcutpoint(kk, jj, ii, ntopol, ntrimax, &
                                xstag, ystag, zstag, ddx, ddy, ddz, topol, &
                                kanteu, triau, i, j, k, 1, itri)
                        END IF

                        IF (foundy == 1) THEN
                            CALL insertcutpoint(kk, jj, ii, ntopol, ntrimax, &
                                xstag, ystag, zstag, ddx, ddy, ddz, topol, &
                                kantev, triav, i, j, k, 2, itri)
                        END IF

                        IF (foundz == 1) THEN
                            CALL insertcutpoint(kk, jj, ii, ntopol, ntrimax, &
                                xstag, ystag, zstag, ddx, ddy, ddz, topol, &
                                kantew, triaw, i, j, k, 3, itri)
                        END IF

                    END DO
                END DO
            END DO
        END DO triangles

        DO i = 1, ii
            DO j = 1, jj
                DO k = 1, kk
                    IF (ABS(kanteu(k, j, i)) < TINY(1.0_realk)) THEN
                        kanteu(k, j, i) = 1.0_realk
                    END IF
                    IF (ABS(kantev(k, j, i)) < TINY(1.0_realk)) THEN
                        kantev(k, j, i) = 1.0_realk
                    END IF
                    IF (ABS(kantew(k, j, i)) < TINY(1.0_realk)) THEN
                        kantew(k, j, i) = 1.0_realk
                    END IF
                END DO
            END DO
        END DO

    END SUBROUTINE cutcorner_grid


    SUBROUTINE insertcutpoint(kk, jj, ii, ntopol, ntrimax, xstag, ystag, &
            zstag, ddx, ddy, ddz, topol, kanteu, triau, &
            i, j, k, direction, itrineu)

        ! Subroutine arguments
        INTEGER(intk), INTENT(in) :: kk, jj, ii, ntopol, ntrimax
        REAL(realk), INTENT(in) :: xstag(ii), ystag(jj), zstag(kk), &
            ddx(ii), ddy(jj), ddz(kk), topol(3, 3, ntopol)
        REAL(realk), INTENT(inout) :: kanteu(kk, jj, ii)
        INTEGER(intk), INTENT(inout) :: triau(ntrimax, kk, jj, ii)
        INTEGER(intk), INTENT(in) :: i, j, k, direction, itrineu

        ! Local variables
        INTEGER(intk) :: itri
        REAL(realk) :: xs, ys, zs
        REAL(realk) :: xp, yp, zp
        REAL(realk) :: x(3)
        INTEGER(intk) :: ind(3)

        ! ntrimax=2 in usual blocking, =1 in srep and calcvolumesexact
        IF (-kanteu(k,j,i) < ntrimax) THEN
            kanteu(k, j, i) = kanteu(k, j, i) - 1.0_realk
            triau(-INT(kanteu(k, j, i)),k, j, i) = itrineu
        ELSE IF (ntrimax == 2) THEN
            ! In case of ntrimax=2, maximal two intersection points per corner
            ! are stored. The following code chooses, which ones to store.
            ! (the ones mit the min and max coordinate in directino of the
            ! corner)

            xp = xstag(i)
            yp = ystag(j)
            zp = zstag(k)

            ! TODO: Is this numerically safe? Use xstag(i-1) instead like in
            ! intinface?
            IF (direction == 1) xp = xp - ddx(i)
            IF (direction == 2) yp = yp - ddy(j)
            IF (direction == 3) zp = zp - ddz(k)

            ! Schnittpunkt erstes Dreieck
            itri = triau(1, k, j, i)
            CALL calcint(xp, yp, zp, xs, ys, zs, &
                itri, ntopol, direction, topol)
            IF (direction == 1) x(1) = xs
            IF (direction == 2) x(1) = ys
            IF (direction == 3) x(1) = zs

            ! Schnittpunkt zweites Dreieck
            itri = triau(2, k, j, i)
            CALL calcint(xp, yp, zp, xs, ys, zs, &
                itri, ntopol, direction, topol)
            IF (direction == 1) x(2) = xs
            IF (direction == 2) x(2) = ys
            IF (direction == 3) x(2) = zs

            ! Schnittpunkt neues Dreieck
            itri = itrineu
            CALL calcint(xp, yp, zp, xs, ys, zs, &
                itri, ntopol, direction, topol)
            IF (direction == 1) x(3) = xs
            IF (direction == 2) x(3) = ys
            IF (direction == 3) x(3) = zs

            CALL sortrx(3, x, ind)

            IF (ind(2) == 1) THEN
                ! First intersection point between second and the new
                triau(1, k, j, i) = itrineu
            ELSE IF (ind(2) == 2) THEN
                ! Second intersection point between second and the new
                triau(2, k, j, i) = itrineu
            END IF
        END IF
    END SUBROUTINE insertcutpoint


    PURE SUBROUTINE calcint(xp, yp, zp, xx, yy, zz, itri, &
            ntopol, direction, topol)

        ! Subroutine arguments
        REAL(realk), INTENT(in) ::  xp, yp, zp
        REAL(realk), INTENT(out) :: xx, yy, zz
        INTEGER(intk), INTENT(in) :: itri
        INTEGER(intk), INTENT(in) :: ntopol, direction
        REAL(realk), INTENT(in) :: topol(3, 3, ntopol)

        ! Local variables
        REAL(realk) :: a, b, c, x1, y1, z1, x2, y2, z2, x3, y3, z3
        REAL(realk) :: betrag

        ! Corner of triangle
        x1 = topol(1, 1, itri)
        x2 = topol(2, 1, itri)
        x3 = topol(3, 1, itri)
        y1 = topol(1, 2, itri)
        y2 = topol(2, 2, itri)
        y3 = topol(3, 2, itri)
        z1 = topol(1, 3, itri)
        z2 = topol(2, 3, itri)
        z3 = topol(3, 3, itri)

        ! Normal vector
        a = (y2-y1)*(z3-z1) - (z2-z1)*(y3-y1)
        b = (z2-z1)*(x3-x1) - (x2-x1)*(z3-z1)
        c = (x2-x1)*(y3-y1) - (y2-y1)*(x3-x1)

        ! Normalization
        betrag = sqrt(a**2 + b**2 + c**2)
        a = a/betrag
        b = b/betrag
        c = c/betrag

        SELECT CASE(direction)
        CASE (1)
            xx = (-c*(zp-z1)-b*(yp-y1))/a + x1
            yy = yp
            zz = zp
        CASE (2)
            xx = xp
            yy = (-c*(zp-z1)-a*(xp-x1))/b + y1
            zz = zp
        CASE (3)
            xx = xp
            yy = yp
            zz = (-b*(yp-y1)-a*(xp-x1))/c + z1
        END SELECT
    END SUBROUTINE calcint


    PURE SUBROUTINE intinface(kk, jj, ii, xstag, ystag, zstag, ddx, ddy, ddz, &
            i, j, k, a, b, c, xp, yp, zp, &
            x1, y1, z1, x2, y2, z2, x3, y3, z3, &
            eps, xmin, xmax, ymin, ymax, zmin, zmax,&
            found, foundx, foundy, foundz, exactedge)

        ! Subroutine arguments
        INTEGER(intk), INTENT(in) :: kk, jj, ii
        REAL(realk), INTENT(in) :: xstag(ii), ystag(jj), zstag(kk), &
            ddx(ii), ddy(jj), ddz(kk)
        INTEGER(intk), INTENT(in) :: i, j, k
        REAL(realk), INTENT(in) :: a, b, c
        REAL(realk), INTENT(in) :: xp, yp, zp
        REAL(realk), INTENT(in) :: x1, y1, z1, x2, y2, z2, x3, y3, z3
        REAL(realk), INTENT(in) :: eps
        REAL(realk), INTENT(in) :: xmin, xmax, ymin, ymax, zmin, zmax
        INTEGER(intk), INTENT(out) :: found, foundx, foundy, foundz
        LOGICAL, INTENT(in) :: exactedge

        ! Local variables
        REAL(realk) :: betrag
        REAL(realk) :: px, py, pz, abx, aby, abz, s1, s2, s3
        REAL(realk) :: xx, yy, zz
        REAL(realk) :: xfront, yright, zbottom
        REAL(realk) :: epsedge

        ! Initialize return variables to default value
        found = 0
        foundx = 0
        foundy = 0
        foundz = 0

        ! xstag(i)-ddx(i) is numerically different from xstag(i-1) - make sure
        ! both cells sharing a face "sees" the same value! Small errors still
        ! done on front, right, bottom faces of a grid!
        IF (i == 1) THEN
            xfront = xstag(1) - ddx(1) - eps
        ELSE
            xfront = xstag(i-1)
        END IF

        IF (j == 1) THEN
            yright = ystag(1) - ddy(1) - eps
        ELSE
            yright = ystag(j-1)
        END IF

        IF (k == 1) THEN
            zbottom = zstag(1) - ddz(1) - eps
        ELSE
            zbottom = zstag(k-1)
        END IF

        ! Controls what happens when an intersection is in between
        ! two edges edge. The exact edge intersection will only cut the
        ! lower side edge, otherwise bot are marked as intersected.
        IF (exactedge) THEN
            epsedge = 0.0
        ELSE
            epsedge = eps
        END IF

        IF (ABS(a) > TINY(1.0_realk)) THEN
            ! x-Kante
            xx = (-c*(zp-z1)-b*(yp-y1))/a + x1
            yy = yp
            zz = zp

            ! This checks if the intersection point is somewhere on the edge
            ! itself
            IF ((xx > xfront - epsedge) .AND. (xx <= xstag(i) + epsedge)) THEN
                ! d1 Vektor berechnen = Gerade von Dreiecksseitenmittelpkt.
                ! zu Schnittpunkt xx,yy,zz
                px = xx-0.5*(x2+x1)
                py = yy-0.5*(y2+y1)
                pz = zz-0.5*(z2+z1)

                ! Hier nur ein paar Abkürzungen
                abx = (x2-x1)
                aby = (y2-y1)
                abz = (z2-z1)

                ! normieren
                betrag = SQRT(abx**2 + aby**2 + abz**2)
                abx = abx/betrag
                aby = aby/betrag
                abz = abz/betrag

                ! ab=[abx,aby,abz] Einheitsvektor von Eckpunkt 1 zu 2
                ! n=[a,b,c]        Normalenvektor des Dreiecks
                ! p=[px,py,pz]     von Dreieckseitenmitelpunkt zum Schniittpunkt

                ! Jetzt p * ( ab x n )
                ! Vorzeichen vertauscht, um lt.eps zu schreiben
                s1 = px*(aby*c-abz*b) + py*(abz*a-abx*c) + pz*(abx*b-aby*a)

                px = xx-0.5*(x3+x2)
                py = yy-0.5*(y3+y2)
                pz = zz-0.5*(z3+z2)

                abx = (x3-x2)
                aby = (y3-y2)
                abz = (z3-z2)

                ! normieren
                betrag = SQRT(abx**2 + aby**2 + abz**2)
                abx = abx/betrag
                aby = aby/betrag
                abz = abz/betrag

                s2 = px*(aby*c-abz*b) + py*(abz*a-abx*c) + pz*(abx*b-aby*a)

                px = xx-0.5*(x1+x3)
                py = yy-0.5*(y1+y3)
                pz = zz-0.5*(z1+z3)

                abx = (x1-x3)
                aby = (y1-y3)
                abz = (z1-z3)

                ! normieren
                betrag = sqrt(abx**2 + aby**2 + abz**2)
                abx = abx/betrag
                aby = aby/betrag
                abz = abz/betrag

                s3 = px*(aby*c-abz*b) + py*(abz*a-abx*c) + pz*(abx*b-aby*a)

                ! Cheks that the intersection point is within the triangle.
                ! Please refer to Nikolaus Peller' theis, "Numerische
                ! Simulation turbulenter Strömungen mit Immersed Boundaries"
                ! page 44-46.
                IF ((s1 <= eps) .AND. (s2 <= eps) .AND. (s3 <= eps)) THEN
                    IF ((xx >= xmin-eps) .AND. (xx <= xmax+eps) .AND. &
                            (yy >= ymin-eps) .AND. (yy <= ymax+eps) .AND. &
                            (zz >= zmin-eps) .AND. (zz <= zmax+eps) ) THEN
                        found  = 1
                        foundx = 1
                    END IF
                END IF
            END IF
        END IF

        IF (ABS(b) > TINY(1.0_realk)) THEN
            xx = xp
            yy = (-c*(zp-z1)-a*(xp-x1))/b + y1
            zz = zp

            ! This checks if the intersection point is somewhere on the edge
            ! itself
            IF ((yy > yright - epsedge) .AND. (yy <= ystag(j) + epsedge)) THEN
                px = xx-0.5*(x2+x1)
                py = yy-0.5*(y2+y1)
                pz = zz-0.5*(z2+z1)

                abx = (x2-x1)
                aby = (y2-y1)
                abz = (z2-z1)

                ! normieren
                betrag = SQRT(abx**2 + aby**2 + abz**2)
                abx = abx/betrag
                aby = aby/betrag
                abz = abz/betrag

                s1 = px*(aby*c-abz*b) + py*(abz*a-abx*c) + pz*(abx*b-aby*a)
                px = xx-0.5*(x3+x2)
                py = yy-0.5*(y3+y2)
                pz = zz-0.5*(z3+z2)

                abx = (x3-x2)
                aby = (y3-y2)
                abz = (z3-z2)

                ! normieren
                betrag = SQRT(abx**2 + aby**2 + abz**2)
                abx = abx/betrag
                aby = aby/betrag
                abz = abz/betrag

                s2 = px*(aby*c-abz*b) + py*(abz*a-abx*c) + pz*(abx*b-aby*a)

                px = xx-0.5*(x1+x3)
                py = yy-0.5*(y1+y3)
                pz = zz-0.5*(z1+z3)

                abx = (x1-x3)
                aby = (y1-y3)
                abz = (z1-z3)

                ! normieren
                betrag = sqrt(abx**2 + aby**2 + abz**2)
                abx = abx/betrag
                aby = aby/betrag
                abz = abz/betrag

                s3 = px*(aby*c-abz*b) + py*(abz*a-abx*c) + pz*(abx*b-aby*a)

                IF ((s1 <= eps) .AND. (s2 <= eps) .AND. (s3 <= eps)) THEN
                    IF ((xx >= xmin-eps) .AND. (xx <= xmax+eps) .AND. &
                            (yy >= ymin-eps) .AND. (yy <= ymax+eps).AND.&
                            (zz >= zmin-eps) .AND. (zz <= zmax+eps) ) THEN
                        found  = 1
                        foundy = 1
                    END IF
                END IF
            END IF
        END IF

        IF (ABS(c) > TINY(1.0_realk)) THEN
            xx = xp
            yy = yp
            zz = (-b*(yp-y1)-a*(xp-x1))/c + z1

            ! This checks if the intersection point is somewhere on the edge
            ! itself
            IF ((zz > zbottom - epsedge).AND.(zz <= zstag(k) + epsedge)) THEN
                px = xx-0.5*(x2+x1)
                py = yy-0.5*(y2+y1)
                pz = zz-0.5*(z2+z1)

                abx = (x2-x1)
                aby = (y2-y1)
                abz = (z2-z1)
                ! normieren
                betrag = sqrt(abx**2 + aby**2 + abz**2)
                abx = abx/betrag
                aby = aby/betrag
                abz = abz/betrag

                s1 = px*(aby*c-abz*b) + py*(abz*a-abx*c) + pz*(abx*b-aby*a)
                px = xx-0.5*(x3+x2)
                py = yy-0.5*(y3+y2)
                pz = zz-0.5*(z3+z2)

                abx = (x3-x2)
                aby = (y3-y2)
                abz = (z3-z2)
                ! normieren
                betrag = SQRT(abx**2 + aby**2 + abz**2)
                abx = abx/betrag
                aby = aby/betrag
                abz = abz/betrag

                s2 = px*(aby*c-abz*b) + py*(abz*a-abx*c) + pz*(abx*b-aby*a)
                px = xx-0.5*(x1+x3)
                py = yy-0.5*(y1+y3)
                pz = zz-0.5*(z1+z3)

                abx = (x1-x3)
                aby = (y1-y3)
                abz = (z1-z3)
                ! normieren
                betrag = SQRT(abx**2 + aby**2 + abz**2)
                abx = abx/betrag
                aby = aby/betrag
                abz = abz/betrag

                s3 = px*(aby*c-abz*b) + py*(abz*a-abx*c) + pz*(abx*b-aby*a)

                IF ((s1 <= eps) .AND. (s2 <= eps) .AND. (s3 <= eps)) THEN
                    IF ((xx >= xmin-eps) .AND. (xx <= xmax+eps) .AND. &
                            (yy >= ymin-eps) .AND. (yy <= ymax+eps) .AND. &
                            (zz >= zmin-eps) .AND. (zz <= zmax+eps) ) THEN
                        found  = 1
                        foundz = 1
                    END IF
                END IF
            END IF
        END IF
    END SUBROUTINE intinface

END MODULE cutcorner_mod
