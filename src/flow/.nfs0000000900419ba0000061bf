module Grid_Area

   IMPLICIT NONE
   PRIVATE

   !public
   PUBLIC :: set_geometry
   public :: print_geometry
   PUBLIC :: init_Grid_Area, finish_Grid_Area
   REAL(8), ALLOCATABLE, PUBLIC :: show(:, :)
   PUBLIC :: Area_Calculation

   !  private
   LOGICAL, PROTECTED :: has_Grid_Area = .FALSE.
   REAL(8), private :: inR
   REAL(8), private :: outR
   REAL(8), private :: BladeLength
   REAL(8), private :: TankSize
   REAL(8), private :: w
   REAL(8), private :: area_dt
   REAL(8), private :: rotate_dt
   REAL(8), private :: tmax
   REAL(8), private :: centrePoint(2)
   REAL(8), private :: grid_dx
   REAL(8), private :: grid_dy
   REAL(8), private :: previous_startangle1
   REAL(8), private :: previous_startangle2
   REAL(8), private :: angle1
   REAL(8), private :: angle2
   !grid
   PRIVATE :: generate_mesh
   REAL(8), ALLOCATABLE, PRIVATE :: Xgrid(:, :), Ygrid(:, :)
contains
   SUBROUTINE init_Grid_Area()

      ! leaving inactive if no parameters specified
      has_Grid_Area = .FALSE.
      ! IF (.NOT. fort7%exists("/flow/Grid_Area")) RETURN   


      ! retrieving rotation rate vector from parameters.json

      ! display obtained parameters
      ! IF (myid == 0) THEN
      !    WRITE (*, '("grid TERM:")')
      ! END IF

      ! set active
      has_Grid_Area = .TRUE.

   END SUBROUTINE init_Grid_Area

   SUBROUTINE finish_Grid_Area

      ! revoking activity
      has_Grid_Area = .FALSE.

      RETURN

   END SUBROUTINE finish_Grid_Area

   SUBROUTINE generate_mesh()
      IMPLICIT NONE
      INTEGER :: nx, ny, i, j
      REAL(8), ALLOCATABLE :: XX(:), YY(:)
      REAL(8) :: dx_step, dy_step

      nx = INT(TankSize/grid_dx + 0.5)
      ny = INT(TankSize/grid_dy + 0.5)

      dx_step = TankSize/REAL(nx, 8)
      dy_step = TankSize/REAL(ny, 8)

      ALLOCATE (XX(nx + 1), YY(ny + 1))
      ALLOCATE (Xgrid(ny + 1, nx + 1), Ygrid(ny + 1, nx + 1))

      DO i = 1, nx + 1
         XX(i) = (i - 1)*dx_step
      END DO

      DO j = 1, ny + 1
         YY(j) = (j - 1)*dy_step
      END DO

      DO j = 1, ny + 1
         DO i = 1, nx + 1
            Xgrid(j, i) = XX(i)
            Ygrid(j, i) = YY(j)
         END DO
      END DO

      ALLOCATE (show(ny, nx))
      show = 0.0D0
   END SUBROUTINE generate_mesh

   SUBROUTINE set_geometry(rinR, routR, rBlade, rTank, rw, rAreaDT, rRotateDT, rtmax, rCentre, rdx, rdy, rAngle1, rAngle2)
      REAL(8), INTENT(IN) :: rinR, routR, rBlade, rTank, rw, rAreaDT, rRotateDT, rtmax
      REAL(8), INTENT(IN) :: rCentre(2), rdx, rdy, rAngle1, rAngle2

      inR = rinR
      outR = routR
      BladeLength = rBlade
      TankSize = rTank
      w = rw
      area_dt = rAreaDT
      rotate_dt = rRotateDT
      tmax = rtmax
      centrePoint = rCentre
      grid_dx = rdx
      grid_dy = rdy
      previous_startangle1 = rAngle1
      previous_startangle2 = rAngle2
      angle1 = rAngle1 + rAreaDT*rw
      angle2 = rAngle2 + rAreaDT*rw
      CALL generate_mesh()
   END SUBROUTINE set_geometry

   SUBROUTINE print_geometry()
      PRINT *, "==== Current Geometry Parameters ===="
      PRINT *, "Inner Radius (inR): ", inR
      PRINT *, "Outer Radius (outR): ", outR
      PRINT *, "Blade Length: ", BladeLength
      PRINT *, "Tank Size: ", TankSize
      PRINT *, "Angular Velocity (w): ", w
      PRINT *, "Area Time Step (area_dt): ", area_dt
      PRINT *, "Rotation Time Step (rotate_dt): ", rotate_dt
      PRINT *, "Total Rotation Time (tmax): ", tmax
      PRINT *, "Center Position: ", centrePoint
      PRINT *, "Grid dx: ", grid_dx
      PRINT *, "Grid dy: ", grid_dy
      PRINT *, "Blade1 Initial Angle: ", previous_startangle1
      PRINT *, "Blade2 Initial Angle: ", previous_startangle2
      PRINT *, "==== Current Mesh Grid ====", Xgrid
   END SUBROUTINE print_geometry

   function trapezoid_area(points) result(area)
      implicit none
      real(8), intent(in) :: points(4, 2)
      real(8) :: area
      real(8) :: sorted_points(4, 2)
      real(8) :: top_points(2, 2), bottom_points(2, 2)
      real(8) :: top_base, bottom_base, height
      integer :: i, j, max_idx, k
      real(8) :: y_values(4)

      do i = 1, 4
         do j = 1, 2
            sorted_points(i, j) = points(i, j)
         end do
         y_values(i) = points(i, 2)
      end do

      do i = 1, 3
         max_idx = i
         do j = i + 1, 4
            if (sorted_points(j, 2) > sorted_points(max_idx, 2)) then
               max_idx = j
            end if
         end do
         if (max_idx /= i) then
            do k = 1, 2
               call swap(sorted_points(i, k), sorted_points(max_idx, k))
            end do
         end if
      end do

      top_points = sorted_points(1:2, :)
      bottom_points = sorted_points(3:4, :)

      top_base = abs(top_points(2, 1) - top_points(1, 1))
      bottom_base = abs(bottom_points(2, 1) - bottom_points(1, 1))

      height = (top_points(1, 2) + top_points(2, 2))/2.0d0 - (bottom_points(1, 2) + bottom_points(2, 2))/2.0d0

      area = 0.5d0*(top_base + bottom_base)*height

   contains
      subroutine swap(a, b)
         real(8), intent(inout) :: a, b
         real(8) :: temp
         temp = a
         a = b
         b = temp
      end subroutine swap
   end function trapezoid_area

   function find_adjacent(p1, p2) result(adjacent_points)
      implicit none
      integer, intent(in) :: p1, p2
      integer :: grid(2, 2) = reshape([1, 2, 3, 4], [2, 2])
      integer :: adjacent_points(2)
      integer :: i, j, idx

      idx = 1

      do i = 1, 2
         do j = 1, 2
            if (grid(i, j) /= p1 .and. grid(i, j) /= p2) then
               adjacent_points(idx) = grid(i, j)
               idx = idx + 1
            end if
         end do
      end do

   end function find_adjacent

   function triangle_area(point1, point2, point3) result(areaa)
      implicit none
      real(8), intent(in) :: point1(2), point2(2), point3(2)
      real(8) :: areaa, base, height

      base = sqrt((point2(1) - point1(1))**2 + (point2(2) - point1(2))**2)

      height = sqrt((point2(1) - point3(1))**2 + (point2(2) - point3(2))**2)

      areaa = 0.5d0*base*height
   end function triangle_area

   function line_circle_intersection(x1, y1, x2, y2, a, b, r) result(points)
      use, intrinsic :: ieee_arithmetic
      implicit none

      real(8), intent(in) :: x1, y1, x2, y2, a, b, r
      real(8), allocatable :: points(:, :)
      real(8) :: dx, dy, coef_a, coef_b, coef_c, disc, sqrt_disc, t1, t2
      real(8), allocatable :: temp_points(:, :)
      integer :: count
      real(8), parameter :: tol = 1.0d-10

      dx = x2 - x1
      dy = y2 - y1

      coef_a = dx*dx + dy*dy

      ! 线段退化为点
      if (coef_a < tol) then
         if (abs((x1 - a)**2 + (y1 - b)**2 - r**2) < tol) then
            allocate (points(1, 2))
            points(1, 1) = x1
            points(1, 2) = y1
         else
            allocate (points(1, 2))
            points = reshape([ieee_value(0.0_8, ieee_quiet_nan), ieee_value(0.0_8, ieee_quiet_nan)], [1, 2])

         end if
         return
      end if

      coef_b = 2.0d0*(dx*(x1 - a) + dy*(y1 - b))
      coef_c = (x1 - a)**2 + (y1 - b)**2 - r**2
      disc = coef_b**2 - 4.0d0*coef_a*coef_c

      if (disc < -tol) then
         allocate (points(1, 2))
         points = reshape([ieee_value(0.0_8, ieee_quiet_nan), ieee_value(0.0_8, ieee_quiet_nan)], [1, 2])

         return
      end if

      allocate (temp_points(2, 2))
      count = 0

      if (disc >= -tol) then
         sqrt_disc = sqrt(max(disc, 0.0d0))
         t1 = (-coef_b - sqrt_disc)/(2.0d0*coef_a)
         t2 = (-coef_b + sqrt_disc)/(2.0d0*coef_a)

         if (t1 >= -tol .and. t1 <= 1.0d0 + tol) then
            count = count + 1
            temp_points(count, 1) = x1 + t1*dx
            temp_points(count, 2) = y1 + t1*dy
         end if
         if (abs(t2 - t1) > tol .and. t2 >= -tol .and. t2 <= 1.0d0 + tol) then
            count = count + 1
            temp_points(count, 1) = x1 + t2*dx
            temp_points(count, 2) = y1 + t2*dy
         end if
      end if

      if (count == 0) then
         allocate (points(1, 2))
         points = reshape([ieee_value(0.0_8, ieee_quiet_nan), ieee_value(0.0_8, ieee_quiet_nan)], [1, 2])

      else
         allocate (points(count, 2))
         points = temp_points(1:count, :)
      end if

   end function line_circle_intersection

   function lineSegmentIntersection(x1, y1, x2, y2, x3, y3, x4, y4) result(intersection)
      use, intrinsic :: ieee_arithmetic
      implicit none
      real(8), intent(in) :: x1, y1, x2, y2, x3, y3, x4, y4
      real(8) :: dx1, dy1, dx2, dy2, denom, t, u, px, py
      real(8), allocatable :: intersection(:)

      dx1 = x2 - x1
      dy1 = y2 - y1
      dx2 = x4 - x3
      dy2 = y4 - y3

      denom = dx1*dy2 - dy1*dx2

      if (abs(denom) > 1e-10) then  ! Avoid division by zero
         t = ((x3 - x1)*dy2 - (y3 - y1)*dx2)/denom
         u = ((x3 - x1)*dy1 - (y3 - y1)*dx1)/denom

         if (t >= 0.0d0 .and. t <= 1.0d0 .and. u >= 0.0d0 .and. u <= 1.0d0) then
            px = x1 + t*dx1
            py = y1 + t*dy1
            allocate (intersection(2))
            intersection(1) = px
            intersection(2) = py
         else
            allocate (intersection(2))
            intersection = [ieee_value(0.0_8, ieee_quiet_nan), ieee_value(0.0_8, ieee_quiet_nan)]
         end if
      else
         allocate (intersection(2))
         intersection = [ieee_value(0.0_8, ieee_quiet_nan), ieee_value(0.0_8, ieee_quiet_nan)]
      end if
   end function lineSegmentIntersection

   logical function PointInSector(px, py, centre, innerR, outerR, startAngle, endAngle)
      implicit none
      real(8), intent(in) :: px, py
      real(8), intent(in) :: centre(2)
      real(8), intent(in) :: innerR, outerR, startAngle, endAngle
      real(8) :: distance, angle
      real(8) :: modStartAngle, modEndAngle

      modStartAngle = mod(startAngle, 2.0d0*3.141592653589793d0)
      modEndAngle = mod(endAngle, 2.0d0*3.141592653589793d0)

      distance = sqrt((px - centre(1))**2 + (py - centre(2))**2)

      angle = atan2(py - centre(2), px - centre(1))

      if (angle < 0.0d0) then
         angle = angle + 2.0d0*3.141592653589793d0
      end if

      if (distance >= innerR .and. distance <= outerR .and. &
          ((modStartAngle <= modEndAngle .and. angle >= modStartAngle .and. angle <= modEndAngle) .or. &
           (modStartAngle > modEndAngle .and. (angle >= modStartAngle .or. angle <= modEndAngle)))) then
         PointInSector = .true.
      else
         PointInSector = .false.
      end if

   end function PointInSector

   subroutine Area_Calculation(step)
      use, intrinsic :: ieee_arithmetic
      real(8)::sector1_startAngle
      real(8)::sector1_endAngle
      real(8)::sector2_startAngle
      real(8)::sector2_endAngle
      real(8)::add_angle
      real(8)::t
      REAL(8) :: innerStartX1, innerStartY1
      REAL(8) :: innerEndX1, innerEndY1
      REAL(8) :: OuterStartX1, OuterStartY1
      REAL(8) :: OuterEndX1, OuterEndY1
      REAL(8) :: innerStartX2, innerStartY2
      REAL(8) :: innerEndX2, innerEndY2
      REAL(8) :: OuterStartX2, OuterStartY2
      REAL(8) :: OuterEndX2, OuterEndY2
      integer::step
      INTEGER :: i, j
      INTEGER :: ni, nj
      real(8) :: x1, y1, x2, y2, x3, y3, x4, y4
      real(8) :: xpoints(4), ypoints(4)
      integer :: q(4)
      integer :: mapping(4, 2)
      integer :: k
      integer :: idx, m
      integer :: New_m, m_count
      integer :: New_idx(2)
      integer :: adjacent_points(2)
      real(8) :: orthoPoint(2)
      integer :: calcupoint(2)
      integer :: New_New_m, New_New_idx
      real(8) :: New_orthoPoint(2)
      integer :: New_calcupoint(2)
      real(8) :: coordinateA(2), coordinateB(2), coordinateC(2), coordinateD(2), coordinateE(2), coordinateF(2)
      real(8) :: coordinate1(2), coordinate2(2), coordinate3(2), coordinate4(2)
      real(8) :: coordinate7(2), coordinate8(2), coordinate9(2), coordinate10(2)
      real(8) :: coordinate13(2), coordinate14(2), coordinate15(2), coordinate16(2)
      real(8) :: coordinate19(2), coordinate20(2), coordinate21(2), coordinate22(2)
      real(8) :: coordinate25(2), coordinate26(2), coordinate27(2), coordinate28(2)
      real(8) :: coordinate31(2), coordinate32(2), coordinate33(2), coordinate34(2)
      real(8), allocatable :: coordinate5(:, :), coordinate6(:, :)
      real(8), allocatable :: coordinate11(:, :), coordinate12(:, :)
      real(8), allocatable :: coordinate17(:, :), coordinate18(:, :)
      real(8), allocatable :: coordinate23(:, :), coordinate24(:, :)
      real(8), allocatable :: coordinate29(:, :), coordinate30(:, :)
      real(8), allocatable :: coordinate35(:, :), coordinate36(:, :)
      real(8) ::area, New_area
      real(8), dimension(4, 2) :: trapezoid_Point
      real(8) :: trapezoid_Area_calculated
      t = step*rotate_dt
      add_angle = w*t
      show = 0.0d0
      !sectror1
      sector1_startAngle = previous_startangle1 + add_angle
      sector1_endAngle = angle1 + add_angle
      !sector2
      sector2_startAngle = previous_startangle2 + add_angle
      sector2_endAngle = angle2 + add_angle
      !endPoints-sector1
      innerStartX1 = centrePoint(1) + inR*cos(sector1_startAngle); 
      innerStartY1 = centrePoint(2) + inR*sin(sector1_startAngle); 
      innerEndX1 = centrePoint(1) + inR*cos(sector1_endAngle); 
      innerEndY1 = centrePoint(2) + inR*sin(sector1_endAngle); 
      OuterStartX1 = centrePoint(1) + outR*cos(sector1_startAngle); 
      OuterStartY1 = centrePoint(2) + outR*sin(sector1_startAngle); 
      OuterEndX1 = centrePoint(1) + outR*cos(sector1_endAngle); 
      OuterEndY1 = centrePoint(2) + outR*sin(sector1_endAngle); 
      !endPoints-sector2
      innerStartX2 = centrePoint(1) + inR*cos(sector2_startAngle); 
      innerStartY2 = centrePoint(2) + inR*sin(sector2_startAngle); 
      innerEndX2 = centrePoint(1) + inR*cos(sector2_endAngle); 
      innerEndY2 = centrePoint(2) + inR*sin(sector2_endAngle); 
      OuterStartX2 = centrePoint(1) + outR*cos(sector2_startAngle); 
      OuterStartY2 = centrePoint(2) + outR*sin(sector2_startAngle); 
      OuterEndX2 = centrePoint(1) + outR*cos(sector2_endAngle); 
      OuterEndY2 = centrePoint(2) + outR*sin(sector2_endAngle); 
      !loop every grid
      ni = size(Xgrid, 2)
      nj = size(Xgrid, 1)
      DO i = 1, nj - 1
         DO j = 1, ni - 1
            x1 = Xgrid(i, j); y1 = Ygrid(i, j)
            x2 = Xgrid(i, j + 1); y2 = Ygrid(i, j + 1)
            x3 = Xgrid(i + 1, j); y3 = Ygrid(i + 1, j)
            x4 = Xgrid(i + 1, j + 1); y4 = Ygrid(i + 1, j + 1)

            xpoints = [x1, x2, x3, x4]
            ypoints = [y1, y2, y3, y4]

            q = [0, 0, 0, 0]

            mapping(1, :) = [2, 3]
            mapping(2, :) = [1, 4]
            mapping(3, :) = [1, 4]
            mapping(4, :) = [2, 3]
            Do k = 1, 4
               if (PointInSector(xpoints(k), ypoints(k), centrePoint, inR, outR, sector1_startAngle, sector1_endAngle) .or. &
                   PointInSector(xpoints(k), ypoints(k), centrePoint, inR, outR, sector2_startAngle, sector2_endAngle)) then
                  q(k) = 1
               else
                  q(k) = 0
               end if
            end do
            !4 points
            if (sum(q) == 4) then
               show(i, j) = 1
               !3 points
            elseif (sum(q) == 3) then
               do m = 1, 4
                  if (q(m) == 0) then
                     idx = m
                     exit
                  end if
               end do
               orthoPoint(1) = xpoints(idx)
               orthoPoint(2) = ypoints(idx)
               calcupoint = mapping(idx, :)
               !Calculate A Point
               coordinate1=lineSegmentIntersection(xpoints(idx),ypoints(idx),xpoints(calcupoint(1)),ypoints(calcupoint(1)),innerStartX1,innerStartY1,OuterStartX1,OuterStartY1)
               if (.not. ieee_is_nan(coordinate1(1))) then
                  coordinateA = coordinate1
               else
                  coordinate2=lineSegmentIntersection(xpoints(idx),ypoints(idx),xpoints(calcupoint(1)),ypoints(calcupoint(1)),innerEndX1,innerEndY1,OuterEndX1,OuterEndY1)
                  if (.not. ieee_is_nan(coordinate2(1))) then
                     coordinateA = coordinate2
                  else
                     coordinate3=lineSegmentIntersection(xpoints(idx),ypoints(idx),xpoints(calcupoint(1)),ypoints(calcupoint(1)),innerStartX2,innerStartY2,OuterStartX2,OuterStartY2)
                     if (.not. ieee_is_nan(coordinate3(1))) then
                        coordinateA = coordinate3
                     else
                      coordinate4=lineSegmentIntersection(xpoints(idx),ypoints(idx),xpoints(calcupoint(1)),ypoints(calcupoint(1)),innerEndX2,innerEndY2,OuterEndX2,OuterEndY2)
                        if (.not. ieee_is_nan(coordinate4(1))) then
                           coordinateA = coordinate4
                        else
                           coordinate5 = line_circle_intersection(xpoints(calcupoint(1)), ypoints(calcupoint(1)), xpoints(idx), ypoints(idx), centrePoint(1), centrePoint(2), inR)
                           if (.not. ieee_is_nan(coordinate5(1, 1))) then
                              coordinateA = coordinate5(1, :)
                           else
                              coordinate6 = line_circle_intersection(xpoints(calcupoint(1)), ypoints(calcupoint(1)), xpoints(idx), ypoints(idx), centrePoint(1), centrePoint(2), outR)
                              if (.not. ieee_is_nan(coordinate6(1, 1))) then
                                 coordinateA = coordinate6(1, :)
                              end if
                           end if
                        end if
                     end if
                  end if
               end if
               !Calculate B Point
               coordinate7=lineSegmentIntersection(xpoints(idx),ypoints(idx),xpoints(calcupoint(2)),ypoints(calcupoint(2)),innerStartX1,innerStartY1,OuterStartX1,OuterStartY1)
               if (.not. ieee_is_nan(coordinate7(1))) then
                  coordinateB = coordinate7
               else
                  coordinate8=lineSegmentIntersection(xpoints(idx),ypoints(idx),xpoints(calcupoint(2)),ypoints(calcupoint(2)),innerEndX1,innerEndY1,OuterEndX1,OuterEndY1)
                  if (.not. ieee_is_nan(coordinate8(1))) then
                     coordinateB = coordinate8
                  else
                     coordinate9=lineSegmentIntersection(xpoints(idx),ypoints(idx),xpoints(calcupoint(2)),ypoints(calcupoint(2)),innerStartX2,innerStartY2,OuterStartX2,OuterStartY2)
                     if (.not. ieee_is_nan(coordinate9(1))) then
                        coordinateB = coordinate9
                     else
                      coordinate10=lineSegmentIntersection(xpoints(idx),ypoints(idx),xpoints(calcupoint(2)),ypoints(calcupoint(2)),innerEndX2,innerEndY2,OuterEndX2,OuterEndY2)
                        if (.not. ieee_is_nan(coordinate10(1))) then
                           coordinateB = coordinate10
                        else
                           coordinate11 = line_circle_intersection(xpoints(calcupoint(2)), ypoints(calcupoint(2)), xpoints(idx), ypoints(idx), centrePoint(1), centrePoint(2), inR)
                           if (.not. ieee_is_nan(coordinate11(1, 1))) then
                              coordinateB = coordinate11(1, :)
                           else
                              coordinate12 = line_circle_intersection(xpoints(calcupoint(2)), ypoints(calcupoint(2)), xpoints(idx), ypoints(idx), centrePoint(1), centrePoint(2), outR)
                              if (.not. ieee_is_nan(coordinate12(1, 1))) then
                                 coordinateB = coordinate12(1, :)
                              end if
                           end if
                        end if
                     end if
                  end if
               end if
               !calculate cut area
               area = triangle_area(coordinateA, orthoPoint, coordinateB)
               show(i, j) = ((grid_dx*grid_dy) - area)/(grid_dx*grid_dy)
               ! 2 points
            elseif (sum(q) == 2) then
               m_count = 0
               do New_m = 1, 4
                  if (q(New_m) == 0) then
                     m_count = m_count + 1
                     New_idx(m_count) = New_m
                     if (m_count == 2) exit
                  end if
               end do
               adjacent_points = find_adjacent(New_idx(1), New_idx(2))
               !Calculate C Point
               coordinate13=lineSegmentIntersection(xpoints(New_idx(1)),ypoints(New_idx(1)),xpoints(adjacent_points(1)),ypoints(adjacent_points(1)),innerStartX1,innerStartY1,OuterStartX1,OuterStartY1)
               if (.not. ieee_is_nan(coordinate13(1))) then
                  coordinateC = coordinate13
               else
                  coordinate14=lineSegmentIntersection(xpoints(New_idx(1)),ypoints(New_idx(1)),xpoints(adjacent_points(1)),ypoints(adjacent_points(1)),innerEndX1,innerEndY1,OuterEndX1,OuterEndY1)
                  if (.not. ieee_is_nan(coordinate14(1))) then
                     coordinateC = coordinate14
                  else
                  coordinate15=lineSegmentIntersection(xpoints(New_idx(1)),ypoints(New_idx(1)),xpoints(adjacent_points(1)),ypoints(adjacent_points(1)),innerStartX2,innerStartY2,OuterStartX2,OuterStartY2)
                     if (.not. ieee_is_nan(coordinate15(1))) then
                        coordinateC = coordinate15
                     else
                  coordinate16=lineSegmentIntersection(xpoints(New_idx(1)),ypoints(New_idx(1)),xpoints(adjacent_points(1)),ypoints(adjacent_points(1)),innerEndX2,innerEndY2,OuterEndX2,OuterEndY2)
                        if (.not. ieee_is_nan(coordinate16(1))) then
                           coordinateC = coordinate16
                        else
                     coordinate17=line_circle_intersection(xpoints(adjacent_points(1)),ypoints(adjacent_points(1)),xpoints(New_idx(1)),ypoints(New_idx(1)),centrePoint(1), centrePoint(2), inR)
                           if (.not. ieee_is_nan(coordinate17(1, 1))) then
                              coordinateC = coordinate17(1, :)
                           else
                         coordinate18=line_circle_intersection(xpoints(adjacent_points(1)),ypoints(adjacent_points(1)),xpoints(New_idx(1)),ypoints(New_idx(1)),centrePoint(1), centrePoint(2), outR)
                              if (.not. ieee_is_nan(coordinate18(1, 1))) then
                                 coordinateC = coordinate18(1, :)
                              end if
                           end if
                        end if
                     end if
                  end if
               end if
               !Calculate D Point
               coordinate19=lineSegmentIntersection(xpoints(New_idx(2)),ypoints(New_idx(2)),xpoints(adjacent_points(2)),ypoints(adjacent_points(2)),innerStartX1,innerStartY1,OuterStartX1,OuterStartY1)
               if (.not. ieee_is_nan(coordinate19(1))) then
                  coordinateD = coordinate19
               else
                  coordinate20=lineSegmentIntersection(xpoints(New_idx(2)),ypoints(New_idx(2)),xpoints(adjacent_points(2)),ypoints(adjacent_points(2)),innerEndX1,innerEndY1,OuterEndX1,OuterEndY1)
                  if (.not. ieee_is_nan(coordinate20(1))) then
                     coordinateD = coordinate20
                  else
                  coordinate21=lineSegmentIntersection(xpoints(New_idx(2)),ypoints(New_idx(2)),xpoints(adjacent_points(2)),ypoints(adjacent_points(2)),innerStartX2,innerStartY2,OuterStartX2,OuterStartY2)
                     if (.not. ieee_is_nan(coordinate21(1))) then
                        coordinateD = coordinate21
                     else
                  coordinate22=lineSegmentIntersection(xpoints(New_idx(2)),ypoints(New_idx(2)),xpoints(adjacent_points(2)),ypoints(adjacent_points(2)),innerEndX2,innerEndY2,OuterEndX2,OuterEndY2)
                        if (.not. ieee_is_nan(coordinate22(1))) then
                           coordinateD = coordinate22
                        else
                     coordinate23=line_circle_intersection(xpoints(adjacent_points(2)),ypoints(adjacent_points(2)),xpoints(New_idx(2)),ypoints(New_idx(2)),centrePoint(1), centrePoint(2), inR)
                           if (.not. ieee_is_nan(coordinate23(1, 1))) then
                              coordinateD = coordinate23(1, :)
                           else
                         coordinate24=line_circle_intersection(xpoints(adjacent_points(2)),ypoints(adjacent_points(2)),xpoints(New_idx(2)),ypoints(New_idx(2)),centrePoint(1), centrePoint(2), outR)
                              if (.not. ieee_is_nan(coordinate24(1, 1))) then
                                 coordinateD = coordinate24(1, :)
                              end if
                           end if
                        end if
                     end if
                  end if
               end if
               !Calculate trapzoid area
               trapezoid_Point(1, 1) = xpoints(New_idx(1))
               trapezoid_Point(1, 2) = ypoints(New_idx(1))
               trapezoid_Point(2, 1) = xpoints(New_idx(2))
               trapezoid_Point(2, 2) = ypoints(New_idx(2))
               trapezoid_Point(3, 1) = coordinateC(1)
               trapezoid_Point(3, 2) = coordinateC(2)
               trapezoid_Point(4, 1) = coordinateD(1)
               trapezoid_Point(4, 2) = coordinateD(2)
               trapezoid_Area_calculated = trapezoid_area(trapezoid_Point)
               show(i, j) = ((grid_dx*grid_dy) - trapezoid_Area_calculated)/(grid_dx*grid_dy)
            elseif (sum(q) == 1) then
               do New_New_m = 1, 4
                  if (q(New_New_m) == 1) then
                     New_New_idx = New_New_m
                     exit
                  end if
               end do
               New_orthoPoint(1) = xpoints(New_New_idx)
               New_orthoPoint(2) = ypoints(New_New_idx)
               New_calcupoint = mapping(New_New_idx, :)
               !Calculate E point
               coordinate25=lineSegmentIntersection(xpoints(New_New_idx),ypoints(New_New_idx),xpoints(New_calcupoint(1)),ypoints(New_calcupoint(1)),innerStartX1,innerStartY1,OuterStartX1,OuterStartY1)
               if (.not. ieee_is_nan(coordinate25(1))) then
                  coordinateE = coordinate25
               else
                   coordinate26=lineSegmentIntersection(xpoints(New_New_idx),ypoints(New_New_idx),xpoints(New_calcupoint(1)),ypoints(New_calcupoint(1)),innerEndX1,innerEndY1,OuterEndX1,OuterEndY1)
                  if (.not. ieee_is_nan(coordinate26(1))) then
                     coordinateE = coordinate26
                  else
                  coordinate27=lineSegmentIntersection(xpoints(New_New_idx),ypoints(New_New_idx),xpoints(New_calcupoint(1)),ypoints(New_calcupoint(1)),innerStartX2,innerStartY2,OuterStartX2,OuterStartY2)
                     if (.not. ieee_is_nan(coordinate27(1))) then
                        coordinateE = coordinate27
                     else
                  coordinate28=lineSegmentIntersection(xpoints(New_New_idx),ypoints(New_New_idx),xpoints(New_calcupoint(1)),ypoints(New_calcupoint(1)),innerEndX2,innerEndY2,OuterEndX2,OuterEndY2)
                        if (.not. ieee_is_nan(coordinate28(1))) then
                           coordinateE = coordinate28
                        else
                 coordinate29=line_circle_intersection(xpoints(New_calcupoint(1)),ypoints(New_calcupoint(1)),xpoints(New_New_idx),ypoints(New_New_idx), centrePoint(1), centrePoint(2), inR) 
                           if (.not. ieee_is_nan(coordinate29(1, 1))) then
                              coordinateE = coordinate29(1, :)
                           else
                               coordinate30=line_circle_intersection(xpoints(New_calcupoint(1)),ypoints(New_calcupoint(1)), xpoints(New_New_idx),ypoints(New_New_idx),centrePoint(1), centrePoint(2), outR)
                              if (.not. ieee_is_nan(coordinate30(1, 1))) then
                                 coordinateE = coordinate30(1, :)
                              end if
                           end if
                        end if
                     end if
                  end if
               end if
               !Calculate F point
               coordinate31=lineSegmentIntersection(xpoints(New_New_idx),ypoints(New_New_idx),xpoints(New_calcupoint(2)),ypoints(New_calcupoint(2)),innerStartX1,innerStartY1,OuterStartX1,OuterStartY1)
               if (.not. ieee_is_nan(coordinate31(1))) then
                  coordinateF = coordinate31
               else
                   coordinate32=lineSegmentIntersection(xpoints(New_New_idx),ypoints(New_New_idx),xpoints(New_calcupoint(2)),ypoints(New_calcupoint(2)),innerEndX1,innerEndY1,OuterEndX1,OuterEndY1)
                  if (.not. ieee_is_nan(coordinate32(1))) then
                     coordinateF = coordinate32
                  else
                  coordinate33=lineSegmentIntersection(xpoints(New_New_idx),ypoints(New_New_idx),xpoints(New_calcupoint(2)),ypoints(New_calcupoint(2)),innerStartX2,innerStartY2,OuterStartX2,OuterStartY2)
                     if (.not. ieee_is_nan(coordinate33(1))) then
                        coordinateF = coordinate33
                     else
                  coordinate34=lineSegmentIntersection(xpoints(New_New_idx),ypoints(New_New_idx),xpoints(New_calcupoint(2)),ypoints(New_calcupoint(2)),innerEndX2,innerEndY2,OuterEndX2,OuterEndY2)
                        if (.not. ieee_is_nan(coordinate34(1))) then
                           coordinateF = coordinate34
                        else
                 coordinate35=line_circle_intersection(xpoints(New_calcupoint(2)),ypoints(New_calcupoint(2)),xpoints(New_New_idx),ypoints(New_New_idx), centrePoint(1), centrePoint(2), inR) 
                           if (.not. ieee_is_nan(coordinate35(1, 1))) then
                              coordinateF = coordinate35(1, :)
                           else
                               coordinate36=line_circle_intersection(xpoints(New_calcupoint(2)),ypoints(New_calcupoint(2)), xpoints(New_New_idx),ypoints(New_New_idx),centrePoint(1), centrePoint(2), outR)
                              if (.not. ieee_is_nan(coordinate36(1, 1))) then
                                 coordinateF = coordinate36(1, :)
                              end if
                           end if
                        end if
                     end if
                  end if
               end if
               New_area = triangle_area(coordinateE, New_orthoPoint, coordinateF)
               show(i, j) = (New_area)/(grid_dx*grid_dy)
            end if

         end do
      end do

   end subroutine Area_Calculation

end module Grid_Area

PROGRAM learn
   USE Grid_Area
   IMPLICIT NONE

   REAL(8) :: rinR, routR, rBlade, rTank, rw
   REAL(8) :: rAreaDT, rRotateDT, rtmax
   REAL(8) :: rCentre(2), rdx, rdy, rAngle1, rAngle2
   integer ::L, nsteps
   ! user input
   PRINT *, "Enter inner radius:"
   READ (*, *) rinR
   PRINT *, "Enter outer radius:"
   READ (*, *) routR
   PRINT *, "Enter blade length:"
   READ (*, *) rBlade
   PRINT *, "Enter tank size:"
   READ (*, *) rTank
   PRINT *, "Enter angular velocity:"
   READ (*, *) rw
   PRINT *, "Enter area time step:"
   READ (*, *) rAreaDT
   PRINT *, "Enter rotate time step:"
   READ (*, *) rRotateDT
   PRINT *, "Enter total rotate time:"
   READ (*, *) rtmax
   PRINT *, "Enter center point coordinates (x y):"
   READ (*, *) rCentre(1), rCentre(2)
   PRINT *, "Enter grid dx:"
   READ (*, *) rdx
   PRINT *, "Enter grid dy:"
   READ (*, *) rdy
   PRINT *, "Enter blade 1 initial angle:"
   READ (*, *) rAngle1
   PRINT *, "Enter blade 2 initial angle:"
   READ (*, *) rAngle2

   CALL set_geometry(rinR, routR, rBlade, rTank, rw, rAreaDT, rRotateDT, rtmax, &
                     rCentre, rdx, rdy, rAngle1, rAngle2)

   CALL print_geometry
   !calculate grid
   call Area_Calculation(15)

   nsteps = INT(aint(rtmax/rRotateDT))
   do L = 0, nsteps
      call Area_Calculation(L)
   end do
   print *, show
END PROGRAM learn

