module Target_Region2d_Poly_mod
!
!  Geometry operations on regions defined by polygons.
!
!  The polygons are passed by the geom structure.
!
!  Used Modules:
!  ============
!
   use Const_mod,                only: INT4, INT8, REAL8,  &
                                       ZERO, HALF, BIG,  &
                                       ENUM

   use Geom_Polygon_mod,         only: isInsidePolygon, distPolygon, distRectPolygon

   use Math_IntMath_mod,         only: index1to2

   use Param_Geom_mod,           only: geom

   use Target_Surface2d_SDF_mod, only: grid
!
!  Default Type and Access:
!  =======================
!
   implicit none
   private
   save
!
!  Public Procedures:
!  =================
!
   public :: getDist2FromRegion2dPoly        ! determine square of distance from region
!            ========================
   public :: getDist2RectRegion2dPoly        ! determine square of distance from region
!            ========================
   public :: isOutsideRegion2dPoly           ! determine if point is outside given region
!            =====================
!
!=======================================================================
!
contains
!
!=======================================================================
!
   subroutine getDist2FromRegion2dPoly( x, z, ir, dist2 )
!             ========================  i  i  i     o
!
!     Determine the square of the distance of the point [x,z] from region ir
!
!  Arguments:
!
      real(REAL8), intent(in) :: x              ! x coordinate of point
      real(REAL8), intent(in) :: z              ! z coordinate of point
      integer(INT4), intent(in) :: ir           ! index of region
      real(REAL8), intent(out) :: dist2         ! square of the distance from region
!
!  Local Variables:
!
!     None
!
!  Start of Executable Code:
!  ------------------------
!
      call distPolygon( dist2, [x, z], geom%region(ir) )
!          ===========    o      i           i

   end subroutine getDist2FromRegion2dPoly
!
!=======================================================================
!
    subroutine getDist2RectRegion2dPoly( ic, ir, cinside, dist2 )
!             =========================  i   i      o       o
!
!     Determine the minimum distance of any point inside a rectangular cell from any point
!     inside region ir.
!
!  Arguments:
!
      integer(INT8), intent(in) :: ic           ! cell index
      integer(INT4), intent(in) :: ir           ! index of region
      logical, intent(out) :: cinside           ! flag indicating rectangle center is inside
                                                !    region
      real(REAL8), intent(out) :: dist2         ! square of distance betweeen ractangle and
                                                !    region
!
!  Local Variables:
!
      integer(INT4) :: ii, jj, i1, i2, j1, j2, ix, iz
      real(REAL8) :: x(2), z(2), xx(2), zz(2), xc, zc, dist2_tmp, dummy
!
!  Start of Executable Code:
!  ------------------------
!
      call index1to2( ic, grid%nx, ix, iz )
!          =========  i      i     o   o

      x(:) = grid%x(ix-1:ix)
      z(:) = grid%z(iz-1:iz)
!
!  Check cell center (this will often cause return from subroutine)
!
      xc = HALF*sum( x(:) )
      zc = HALF*sum( z(:) )
      call distRectPolygon( dist2, dummy, cinside,  &
!          ===============    o      o       o
                            x(1), x(2), z(1), z(2), geom%region(ir) )
!                            i     i     i     i           i

      if( cinside ) then
         dist2 = ZERO
         return
      end if
!
!  Check cell corners, considering periodic images of region
!
      if( geom%bndtypex(1) == ENUM%bndtype%periodic ) then
         i1 = -1
      else
         i1 = 0
      end if
      if( geom%bndtypex(2) == ENUM%bndtype%periodic ) then
         i2 = 1
      else
         i2 = 0
      end if
      if( geom%bndtypez(1) == ENUM%bndtype%periodic ) then
         j1 = -1
      else
         j1 = 0
      end if
      if( geom%bndtypez(2) == ENUM%bndtype%periodic ) then
         j2 = 1
      else
         j2 = 0
      end if

      dist2 = BIG

      do ii = i1, i2
         do jj = j1, j2
            xx(:) = x(:) + ii * (geom%xlim(2) - geom%xlim(1))
            zz(:) = z(:) + jj * (geom%zlim(2) - geom%zlim(1))
            call distRectPolygon( dist2_tmp, dummy, cinside,  &
!                ===============      o        o       o
                                  xx(1), xx(2), zz(1), zz(2), geom%region(ir) )
!                                   i      i      i      i           i
            dist2 = min( dist2_tmp, dist2 )
         end do
      end do

   end subroutine getDist2RectRegion2dPoly
!
!=======================================================================
!
   subroutine isOutsideRegion2dPoly( x, z, ir, outside )
!             =====================  i  i  i      o
!
!     Determine the region the point [x,z] is located outside region ir
!
!     "outside" here refers to the position being outside the region polygon. This is not
!     exactly the same as sdf < 0.
!
!  Arguments:
!
      real(REAL8), intent(in) :: x              ! x coordinate of point
      real(REAL8), intent(in) :: z              ! z coordinate of point
      integer(INT4), intent(in) :: ir           ! index of region
      logical, intent(out) :: outside           ! flag indicating point is outside region
!
!  Local Variables:
!
      integer(INT4) :: iloc
!
!  Start of Executable Code:
!  ------------------------
!
      iloc = isInsidePolygon( [x, z], geom%region(ir) )
!            ===============    i           i
      outside = ( iloc < 0 )

   end subroutine isOutsideRegion2dPoly
!
!=======================================================================
!
end module Target_Region2d_Poly_mod
