MODULE calcfacearea_mod
    USE core_mod, ONLY: realk, intk, divide0
    USE ibconst_mod, ONLY: maccur

    IMPLICIT NONE(type, external)
    PRIVATE

    PUBLIC :: calcfacedata, calcwallfacecenter, calcwallfacecenterrescue

CONTAINS
    SUBROUTINE calcfacedata(k, j, i, dir, nptsface, face, npoint, xpoints, &
            area, s1, s2)
        ! Berechnet die flaeche eines konvexen polygons und dessen schwerpunkt

        ! Subroutine arguments
        INTEGER(intk), INTENT(in) :: k, j, i
        INTEGER(intk), INTENT(in) :: dir
        INTEGER(intk), INTENT(in) :: nptsface
        INTEGER(intk), INTENT(in) :: face(nptsface)
        INTEGER(intk), INTENT(in) :: npoint
        REAL(realk), INTENT(in) :: xpoints(3, npoint)
        REAL(realk), INTENT(out) :: area, s1, s2

        ! Local variables
        INTEGER(intk) :: icorner, idir
        REAL(realk) :: x0(3), x1(3), x2(3), x3, x4, x5
        REAL(realk) :: xmin, xmax, ymin, ymax

        area = 0.0
        s1 = 0.0
        s2 = 0.0

        IF (nptsface <= 2) THEN
            RETURN
        END IF

        x0(1) = xpoints(1, face(1))
        x0(2) = xpoints(2, face(1))
        x0(3) = xpoints(3, face(1))

        xmin = x0(MOD(dir+1-1, 3) + 1)
        xmax = x0(MOD(dir+1-1, 3) + 1)
        ymin = x0(MOD(dir+2-1, 3) + 1)
        ymax = x0(MOD(dir+2-1, 3) + 1)

        DO icorner = 1, nptsface
            DO idir = 1, 3
                x1(idir) = xpoints(idir, face(icorner))
                x2(idir) = xpoints(idir, face(MOD(icorner, nptsface)+1))
            END DO

            x3 = ABS((x1(MOD(dir+1-1, 3) + 1) - x0(MOD(dir+1-1, 3) + 1)) &
                * (x2(MOD(dir+2-1, 3) + 1) - x0(MOD(dir+2-1, 3) + 1)) &
                - (x1(MOD(dir+2-1, 3) + 1) - x0(MOD(dir+2-1, 3) + 1)) &
                * (x2(MOD(dir+1-1, 3) + 1) - x0(MOD(dir+1-1, 3) + 1)))

            x4 = (x1(MOD(dir+1-1, 3) + 1) &
                + x2(MOD(dir+1-1, 3) + 1) &
                - 2*x0(MOD(dir+1-1, 3) + 1))

            x5 = (x1(MOD(dir+2-1, 3) + 1) &
                + x2(MOD(dir+2-1, 3) + 1) &
                - 2*x0(MOD(dir+2-1, 3) + 1))

            area = area + 0.5*x3
            s1 = s1 + x4*x3
            s2 = s2 + x5*x3

            xmin = MIN(xmin, x2(MOD(dir+1-1, 3) + 1))
            xmax = MAX(xmax, x2(MOD(dir+1-1, 3) + 1))
            ymin = MIN(ymin, x2(MOD(dir+2-1, 3) + 1))
            ymax = MAX(ymax, x2(MOD(dir+2-1, 3) + 1))
        END DO

        s1 = x0(MOD(dir+1-1, 3) + 1) + divide0(s1, 6*area)
        s2 = x0(MOD(dir+2-1, 3) + 1) + divide0(s2, 6*area)

        IF (s1 < xmin .OR. s1 > xmax) THEN
            WRITE(*, *) 'warn calcfacedata k, j, i', k, j, i
            WRITE(*, *) 'warn calcfacedata s1, xmin, xmax', s1, xmin, xmax
        END IF

        IF (s2 < ymin .OR. s2 > ymax) THEN
            WRITE(*, *) 'warn calcfacedata k,j,i', k, j, i
            WRITE(*, *) 'warn calcfacedata s2, ymin, ymax', s2, ymin, ymax
        END IF
    END SUBROUTINE calcfacedata


    ! TODO: Why is area part of the argument list at all? This routine is only
    ! called from calcauavaw and area is not used there... Remove!
    SUBROUTINE calcwallfacecenter(k, j, i, nwa, npoints, connectwa, xpoints, &
            cartarea, area, s1, s2, s3)
        ! Berechnet den mittelpunkt der koerperschnittflaeche

        ! Subroutine arguments
        INTEGER(intk), INTENT(in) :: k, j, i
        INTEGER(intk), INTENT(in) :: nwa
        INTEGER(intk), INTENT(in) :: npoints
        INTEGER(intk), INTENT(in) :: connectwa(7)
        REAL(realk), INTENT(in) :: xpoints(3, npoints)
        REAL(realk), INTENT(in) :: cartarea
        REAL(realk), INTENT(out) :: area, s1, s2, s3

        ! Local variables
        INTEGER(intk) :: ipoint, idir, iipoint
        REAL(realk) :: cross(3), x1(3), x2(3), x0, y0, z0
        REAL(realk) :: abscross
        REAL(realk) :: s1one, s2one, s3one, areaone

        s1 = 0.0
        s2 = 0.0
        s3 = 0.0
        area = 0.0

        IF (nwa <= 2) THEN
            RETURN
        END IF

        DO iipoint = 1, nwa-2
            ! for each starting point compute cog (s) and area
            x0 = xpoints(1, npoints-20+8+connectwa(iipoint))
            y0 = xpoints(2, npoints-20+8+connectwa(iipoint))
            z0 = xpoints(3, npoints-20+8+connectwa(iipoint))

            s1one = 0.0
            s2one = 0.0
            s3one = 0.0
            areaone = 0.0

            ! divide polygon into triangles
            DO ipoint = 2, nwa-1
                DO idir = 1, 3
                    x1(idir) = xpoints(idir, npoints-20+8+connectwa( &
                        MOD(iipoint-1+ipoint-1, nwa)+1))
                    x2(idir) = xpoints(idir, npoints-20+8+connectwa( &
                        MOD(iipoint-1+ipoint+1-1, nwa)+1))
                END DO

                cross(1) = ((x1(2)-y0)*(x2(3)-z0) - (x1(3)-z0)*(x2(2)-y0))
                cross(2) = ((x1(3)-z0)*(x2(1)-x0) - (x1(1)-x0)*(x2(3)-z0))
                cross(3) = ((x1(1)-x0)*(x2(2)-y0) - (x1(2)-y0)*(x2(1)-x0))
                abscross = SQRT(cross(1)**2 + cross(2)**2 + cross(3)**2)

                s1one = s1one + abscross*(x1(1) + x2(1) - 2*x0)
                s2one = s2one + abscross*(x1(2) + x2(2) - 2*y0)
                s3one = s3one + abscross*(x1(3) + x2(3) - 2*z0)
                areaone = areaone + 0.5*abscross
            END DO

            s1one = x0 + divide0(s1one, 6.0*areaone)
            s2one = y0 + divide0(s2one, 6.0*areaone)
            s3one = z0 + divide0(s3one, 6.0*areaone)

            s1 = s1 + s1one
            s2 = s2 + s2one
            s3 = s3 + s3one
            area = area + areaone
        END DO

        s1 = s1/(nwa-2)
        s2 = s2/(nwa-2)
        s3 = s3/(nwa-2)
        area = area/(nwa-2)
    END SUBROUTINE calcwallfacecenter


    SUBROUTINE calcwallfacecenterrescue(nwa, npoints, connectwa, xpoints, &
            s1, s2, s3)
        ! Subroutine arguments
        INTEGER(intk), INTENT(in) :: nwa, npoints
        INTEGER(intk), INTENT(in) :: connectwa(7)
        REAL(realk), INTENT(in) :: xpoints(3, npoints)
        REAL(realk), INTENT(out) :: s1, s2, s3

        ! Local variables
        INTEGER(intk) :: ipoint, idir
        REAL(realk) :: x1(3)

        s1 = 0.0
        s2 = 0.0
        s3 = 0.0

        DO ipoint = 1, nwa
            DO idir = 1, 3
                x1(idir) = xpoints(idir, npoints-20+8+connectwa(ipoint))
            END DO
            s1 = s1 + x1(1)
            s2 = s2 + x1(2)
            s3 = s3 + x1(3)
        END DO

        s1 = s1/nwa
        s2 = s2/nwa
        s3 = s3/nwa
    END SUBROUTINE calcwallfacecenterrescue
END MODULE calcfacearea_mod
