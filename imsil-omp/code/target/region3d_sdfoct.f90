module Target_Region3d_SDFOct_mod
!
!  Geometry operations on 3D regions using signed distance functions (SDFs) on a hyperoctree.
!
!  Reads in the SDFs from files.
!
!  Used Modules:
!  ============
!
   use Const_mod,                      only: INT4, INT8, REAL8,  &
                                             ZERO, TWO, EIGHTH, BIG,  &
                                             ENUM

   use Target_Surface3d_SDFOct_mod,    only: octree
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
   public :: getDist2CuboidRegion3dSDFOct       ! determine square of distance between grid
!            ============================       !    cell and region
   public :: getRegionIndex3dSDFOct             ! determine the region a given point is located
!            ======================             !    in (returns 0 if outside target)
!
!=======================================================================
!
contains
!
!=======================================================================
!
   subroutine getDist2CuboidRegion3dSDFOct( ic, ir, cinside, dist2 )
!             ============================  i   i      i       o
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
      integer(INT4) :: i, j, k, jr, level
      real(REAL8) :: x(8), y(8), z(8), sdf(8), xc, yc, zc
      real(REAL8) :: dx2, dy2, dz2, dxmin, dymin, dzmin, wx0, wy0, wz0, wx, wy, wz
!
!  Start of Executable Code:
!  ------------------------
!
!  Determine corner SDF values
!
      call octree%getCornersSDFValuesAtLeafNodeId( ic, ir, level, x, y, z, sdf )
!          ======================================  i   i     o    o  o  o   o
!
!  Determine whether cell center is inside region ir
!
      xc = EIGHTH * sum( x(:) )
      yc = EIGHTH * sum( y(:) )
      zc = EIGHTH * sum( z(:) )
      call octree%getRegionIndex( xc, yc, zc, jr )
!          =====================  i   i   i   o
      cinside = ( jr == ir )
!
!  Determine minimum distance between cell and region
!
      if( any( sdf(:) > ZERO ) ) then
         dist2 = ZERO
      else if( any( sdf(:) == -BIG ) ) then
         dist2 = BIG
      else
         dist2 = minval( sdf(:)**2 )  ! corner distances

         if( level > 0 ) then                ! face distances
            ! This branch should never be entered
            error stop 'getDist2CuboidRegion3dSDFOct: level>0 without any sdf>0 or =-BIG'

            call octree%getGridDelta( wx0, wy0, wz0 )
!                ===================   o    o    o
            wx = 2**level * wx0
            wy = 2**level * wy0
            wz = 2**level * wz0
            dxmin = ( 2**level-2 ) * wx0     ! min. distance from next level=0 cell
            dymin = ( 2**level-2 ) * wy0
            dzmin = ( 2**level-2 ) * wz0

            ! Extending this to 3D is possible but complex (see private notes 20.3.2025)
            do i = 1, 8
               do j = i+1, 8
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

   end subroutine getDist2CuboidRegion3dSDFOct
!
!=======================================================================
!
   subroutine getRegionIndex3dSDFOct( x, y, z, ir_, outside_ )
!             ======================  i  i  i   o      o
!
!     Determine the region the point [x,z] is located in (if pos is inside the target)
!     or closest to (if point is outside the target)
!
!  Arguments:
!
      real(REAL8), intent(in) :: x                 ! x coordinate of position to be checked
      real(REAL8), intent(in) :: y                 ! y coordinate of position to be checked
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
      call octree%getRegionIndex( x, y, z, ir, outside )
!          =====================  i  i  i  o      o

      if( present( ir_ ) )  ir_ = ir

      if( present( outside_ ) )  outside_ = outside

   end subroutine getRegionIndex3dSDFOct
!
!=======================================================================
!
end module Target_Region3d_SDFOct_mod
