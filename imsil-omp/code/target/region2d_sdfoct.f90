module Target_Region2d_SDFOct_mod
!
!  Geometry operations on 2D regions using signed distance functions (SDFs) on a hyperoctree.
!
!  Reads in the SDFs from files.
!
!  Used Modules:
!  ============
!
   use Const_mod,                      only: INT4, INT8, REAL8,  &
                                             ZERO, TWO, QUART, BIG,  &
                                             ENUM

   use Target_Surface2d_SDFOct_mod,    only: octree
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
   public :: getDist2RectRegion2dSDFOct         ! determine square of distance between grid
!            ==========================         !    cell and region
   public :: getRegionIndex2dSDFOct             ! determine the region a given point is located
!            ======================             !    in (returns 0 if outside target)
!
!=======================================================================
!
contains
!
!=======================================================================
!
   subroutine getDist2RectRegion2dSDFOct( ic, ir, cinside, dist2 )
!             ==========================  i   i      i       o
!
!     Determine the square of the minimum distance between any point inside cell ic and
!     region ir
!
!  Arguments:
!
      integer(INT8), intent(in) :: ic           ! cell index
      integer(INT4), intent(in) :: ir           ! index of region
      logical, intent(out) :: cinside           ! flag indicating rectangle center is inside
                                                !    region (.false. may also indicate
      real(REAL8), intent(out) :: dist2         ! square of distance betweeen ractangle and
                                                !    region
!
!  Local Variables:
!
      integer(INT4) :: i, j, jr, level
      real(REAL8) :: x(4), z(4), sdf(4), xc, zc
      real(REAL8) :: dx2, dz2, dxmin, dzmin, wx0, wz0, wx, wz
!
!  Start of Executable Code:
!  ------------------------
!
!  Determine corner SDF values
!
      call octree%getCornersSDFValuesAtLeafNodeId( ic, ir, level, x, z, sdf )
!          ======================================  i   i     o    o  o   o
!
!  Determine whether cell center is inside region ir
!
      xc = QUART * sum( x(:) )
      zc = QUART * sum( z(:) )
      call octree%getRegionIndex( xc, zc, jr )
!          =====================  i   i   o
      cinside = ( jr == ir )
!
!  Determine minimum distance between cell and region
!
      if( any( sdf(:) > ZERO ) ) then
         dist2 = ZERO
      else if( any( sdf(:) == -BIG ) ) then
         dist2 = -BIG
      else
         dist2 = minval( sdf(:)**2 )  ! corner distances

         if( level > 0 ) then                ! face distances
            call octree%getGridDelta( wx0, wz0 )
!                ===================   o    o
            wx = 2**level * wx0
            wz = 2**level * wz0
            dxmin = ( 2**level-2 ) * wx0     ! min. distance from next level=0 cell
            dzmin = ( 2**level-2 ) * wz0
            do i = 1, 4
               do j = i+1, 4
                  if( x(i) == x(j) ) then
                     if( abs( sdf(i) ) + abs( sdf(j) ) < wz ) then
                        dx2 = ZERO
                     else
                        dx2 = sdf(i)**2 - ( (sdf(i)**2-sdf(j)**2+wz**2)/(TWO*wz) )**2
                     end if
                     dist2 = min( dist2, max( dx2, dxmin**2 ) )
                  else if( z(i) == z(j) ) then
                     if( abs( sdf(i) ) + abs( sdf(j) ) < wx ) then
                        dz2 = ZERO
                     else
                        dz2 = sdf(i)**2 - ( (sdf(i)**2-sdf(j)**2+wx**2)/(TWO*wx) )**2
                     end if
                     dist2 = min( dist2, max( dz2, dzmin**2 ) )
                  end if
               end do
            end do
         end if
      end if

   end subroutine getDist2RectRegion2dSDFOct
!
!=======================================================================
!
   subroutine getRegionIndex2dSDFOct( x, z, ir_, outside_ )
!             ======================  i  i   o      o
!
!     Determine the region the point [x,z] is located in (if pos is inside the target)
!     or closest to (if point is outside the target)
!
!  Arguments:
!
      real(REAL8), intent(in) :: x                 ! x coordinate of position to be checked
      real(REAL8), intent(in) :: z                 ! z coordinate of position to be checked
      integer(INT4), optional, intent(out) :: ir_  ! index of region
      logical, optional, intent(out) :: outside_   ! flag indicating whether pos is outside any
!                                                  !    region
!  Local Variables:
!
      integer(INT4) :: ir
      logical :: outside
!
!  Start of Executable Code:
!  ------------------------
!
!  Just ask the octree
!
      call octree%getRegionIndex( x, z, ir, outside )
!          =====================  i  i  o      o

      if( present( ir_ ) )  ir_ = ir

      if( present( outside_ ) )  outside_ = outside

   end subroutine getRegionIndex2dSDFOct
!
!=======================================================================
!
end module Target_Region2d_SDFOct_mod
