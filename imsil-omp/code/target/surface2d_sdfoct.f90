module Target_Surface2d_SDFOct_mod
!
!  Define the surface by a distance function (SDF) on a octree grid
!
!  Used Modules:
!  ============
!
   use Const_mod,                      only: INT4, INT8, REAL8,  &
                                             ONE, TWO, EPSMB,  &
                                             ENUM

   use Cinterface_LSHyperoctree2d_mod, only: LSHyperoctree2d_type

   use Param_Setup_mod,                only: nr, lsput
   use Param_Dam_mod,                  only: dam
   use Param_Geom_mod,                 only: geom
   use Param_Out_mod,                  only: out
   use Param_Sn_mod,                   only: sn

   use Stat_His2d_mod,                 only: His2d_type, writeHis2d
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
   public :: initSurface2dSDFOct             ! initialize the octree representing the surface
!            ===================
   public :: getNcellSurface2dSDFOct         ! get the number of cells
!            =======================
   public :: getCellIndex2dSDFOct            ! get cell index at given position
!            ====================
   public :: getDistSurface2dSDFOct          ! get signed distance from surface
!            ======================
!
!  Public Data:
!  ===========
!
   type(LSHyperoctree2d_type), public :: octree
!
!=======================================================================
!
contains
!
!=======================================================================
!
   subroutine initSurface2dSDFOct( xlim, zlim )
!             ===================   o     o
!
!     Initialize the simulation domain and the signed distance function from the surface
!     defined on an octree
!
!  Arguments:
!
      real(REAL8), intent(out) :: xlim(2)       ! minimum and maximum x values to be queried
      real(REAL8), intent(out) :: zlim(2)       ! minimum and maximum z values to be queried
!
!  Local Variables:
!
      integer(INT4) :: bndtypex(2), bndtypez(2)
      real(REAL8) :: denseDistance, wvacuum

      logical :: outside
      integer(INT4) :: nboxx, nboxz, i, j, ir
      real(REAL8) :: wx, wz, x, z
      type(His2d_type) :: his
!
!  Start of Executable Code:
!  ------------------------
!
!  Initialize the octree
!
      denseDistance = max( sqrt( sn%ffpmax**2 + sn%pmax**2 ), geom%wrough )
      if( lsput )  denseDistance = max( dam%dnear, denseDistance )

      octree = LSHyperoctree2d_type( geom%sdffile, nr, denseDistance )
!              ====================       i        i         i
!
!  Get bounding box
!
      call octree%getBoundaryLocationMin( xlim(1), zlim(1) )
!          =============================     o        o
      call octree%getBoundaryLocationMax( xlim(2), zlim(2) )
!          =============================     o        o
!
!  Check boundary types
!
      call octree%getBoundaryConditionMin( bndtypex(1), bndtypez(1) )
!          ==============================       o            o
      call octree%getBoundaryConditionMin( bndtypex(2), bndtypez(2) )
!          ==============================       o            o
      if( ( any( bndtypex(:) /= geom%bndtypex(:) ) .and..not.  &
            ( all( bndtypex(:) == ENUM%bndtype%reflect ) .and.  &
              any( geom%bndtypex(:) == ENUM%bndtype%reflect ) .and.  &
              any( geom%bndtypex(:) == ENUM%bndtype%periodic ) ) ) .or.  &
          any( bndtypez(:) /= geom%bndtypez(:) ) ) then
         error stop 'Boundary types on octree inconsistent with those specified in INP'
      end if
!
!  Adjust bounding box according to boundary types
!
      wvacuum = sn%pmax * (ONE+TWO*EPSMB)
      if( geom%bndtypex(1) == ENUM%bndtype%free )  xlim(1) = xlim(1) - wvacuum
      if( geom%bndtypex(2) == ENUM%bndtype%free )  xlim(2) = xlim(2) + wvacuum
      if( geom%bndtypez(1) == ENUM%bndtype%free )  zlim(1) = zlim(1) - wvacuum
      if( geom%bndtypez(2) == ENUM%bndtype%free )  zlim(2) = zlim(2) + wvacuum
!
!  Write signed distance function to file (for debugging)
!
      if( out%ldist ) then

         wx = 1
         wz = 1
         nboxx = (xlim(2) - xlim(1) + 40) / wx
         nboxz = (zlim(2) - zlim(1) + 40) / wz
         his%hbeg(:) = [ xlim(1)-20 + 0.5*wx, zlim(1)-20 + 0.5*wz ]

!         nboxx = 1000
!         nboxz = 10000
!         wx = (xlim(2)-xlim(1))/nboxx
!         wz = (grid%zmax-zlim(1))/nboxz
!         his%hbeg(:) = [ xlim(1) + 0.5*wx, zlim(1) + 0.5*wz ]

!         nboxx = 1000
!         nboxz = 1000
!         wx = 2
!         wz = 2
!         his%hbeg(:) = [ -1000+0.5*wx, -1000+0.5*wz ]

         his%nbox(:) = [ nboxx, nboxz ]
         his%imin(:) = [ 1, 1 ]
         his%imax(:) = [ nboxx, nboxz ]
         his%wbox(:) = [ wx, wz ]
         allocate( his%indhis( 3 ) )
         allocate( his%hist( 3, nboxx, nboxz ) )

         open( 33, file='../DIST' )
         write( 33, * )  3, nboxz
         write( 33, * )  0, 0, 0

         do i = 1, nboxx
            x = his%hbeg(1)+(i-1)*wx
!            x = 0
!            i = 1
            do j = 1, nboxz
               z = his%hbeg(2)+(j-1)*wz
               call octree%getDistSurface( x, z, his%hist(3,i,j), isOutside=outside )
!                   =====================  i  i       o              o
               if( outside ) then
                  his%hist(1,i,j) = 1
               else
                  his%hist(1,i,j) = 0
               end if
               call octree%getRegionIndex(x, z, ir)
!                   ===================== i  i  o
               his%hist(2,i,j) = ir
               write( 33, * )  z, ir, his%hist(3,i,j)
            end do
         end do
!            close( 33 )
!            error stop

         call writeHis2d( '../DIST2', 'write', his )
!             ==========     i        i      i
         error stop

      end if

   end subroutine initSurface2dSDFOct
!
!=======================================================================
!
   subroutine getNcellSurface2dSDFOct( ncell )
!             =======================    o
!
!     Get the number of cells
!
!  Arguments:
!
      integer(INT8), intent(out) :: ncell       ! number of grid cells
!
!  Local Variables:
!
!     None
!
!  Start of Executable Code:
!  ------------------------
!
      ncell = octree%getNumberOfLeafNodes( )
!             ===========================

   end subroutine getNcellSurface2dSDFOct
!
!=======================================================================
!
   subroutine getCellIndex2dSDFOct( x, z, ic )
!             ====================  i  i  o
!
!     Get the cell index at position [x,z]
!
!  Arguments:
!
      real(REAL8), intent(in) :: x               ! x coordinate of point
      real(REAL8), intent(in) :: z               ! z coordinate of point
      integer(INT8), intent(out) :: ic          ! cell index
!
!  Local Variables:
!
!     None
!
!  Start of Executable Code:
!  ------------------------
!
      ic = octree%getLeafNodeIdAtCoords( x, z )
!          ============================  i  i

   end subroutine getCellIndex2dSDFOct
!
!=======================================================================
!
   subroutine getDistSurface2dSDFOct( x, z, sdf, dirx_, dirz_ )
!             ======================  i  i   o     o      o
!
!     Get the signed distance of a point [x,z] from the surface by interpolation
!     ( + inside, - outside ) and
!     the direction towards/from the surface (pointing outwards the target)
!
!     NOTE: It is assumed that the point is inside the grid
!
!  Arguments:
!
      real(REAL8), intent(in) :: x                 ! x coordinate of point
      real(REAL8), intent(in) :: z                 ! x coordinate of point
      real(REAL8), intent(out) :: sdf              ! signed distance: (+) inside, (-) outside
      real(REAL8), optional, intent(out) :: dirx_  ! x/z component of direction towards/from
      real(REAL8), optional, intent(out) :: dirz_  !    surface
!
!  Local Variables:
!
!     None
!
!  Start of Executable Code:
!  ------------------------
!
      call octree%getDistSurface( x, z, sdf, dirx_, dirz_ )
!          =====================  i  i   o     o      o

   end subroutine getDistSurface2dSDFOct
!
!=======================================================================
!
end module Target_Surface2d_SDFOct_mod
