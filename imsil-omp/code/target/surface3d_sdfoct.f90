module Target_Surface3d_SDFOct_mod
!
!  Define the surface by a distance function (SDF) on a octree grid
!
!  Used Modules:
!  ============
!
   use Const_mod,                      only: INT4, INT8, REAL8,  &
                                             ONE, TWO, EPSMB,  &
                                             ENUM

   use Cinterface_LSHyperoctree3d_mod, only: LSHyperoctree3d_type

   use Param_Setup_mod,                only: nr, lsput
   use Param_Dam_mod,                  only: dam
   use Param_Geom_mod,                 only: geom
   use Param_Out_mod,                  only: out
   use Param_Sn_mod,                   only: sn

   use Stat_His3d_mod,                 only: His3d_type, writeHis3d
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
   public :: initSurface3dSDFOct             ! initialize the octree representing the surface
!            ===================
   public :: getNcellSurface3dSDFOct         ! get the number of cells
!            =======================
   public :: getCellIndex3dSDFOct            ! get cell index at given position
!            ====================
   public :: getDistSurface3dSDFOct          ! get signed distance from surface
!            ======================
!
!  Public Data:
!  ===========
!
   type(LSHyperoctree3d_type), public :: octree
!
!=======================================================================
!
contains
!
!=======================================================================
!
   subroutine initSurface3dSDFOct( xlim, ylim, zlim )
!             ===================   o     o     o
!
!     Initialize the simulation domain and the signed distance function from the surface
!     defined on an octree
!
!  Arguments:
!
      real(REAL8), intent(out) :: xlim(2)       ! minimum and maximum x values to be queried
      real(REAL8), intent(out) :: ylim(2)       ! minimum and maximum y values to be queried
      real(REAL8), intent(out) :: zlim(2)       ! minimum and maximum z values to be queried
!
!  Local Variables:
!
      integer(INT4) :: bndtypex(2), bndtypey(2), bndtypez(2)
      real(REAL8) :: denseDistance, wvacuum

      logical :: outside
      integer(INT4) :: nboxx, nboxy, nboxz, i, j, k, ir
      real(REAL8) :: wx, wy, wz, x, y, z
      type(His3d_type) :: his
!
!  Start of Executable Code:
!  ------------------------
!
!  Initialize the octree
!
      denseDistance = max( sqrt( sn%ffpmax**2 + sn%pmax**2 ), geom%wrough )
      if( lsput )  denseDistance = max( dam%dnear, denseDistance )

      octree = LSHyperoctree3d_type( geom%sdffile, nr, denseDistance )
!              ====================       i        i         i
!
!  Get number of cells, bounding box and boundary types
!
      call octree%getBoundaryLocationMin( xlim(1), ylim(1), zlim(1) )
!          =============================     o        o        o
      call octree%getBoundaryLocationMax( xlim(2), ylim(2), zlim(2) )
!          =============================     o        o        o
      call octree%getBoundaryConditionMin( bndtypex(1), bndtypey(1), bndtypez(1) )
!          ==============================       o            o            o
      call octree%getBoundaryConditionMin( bndtypex(2), bndtypey(2), bndtypez(2) )
!          ==============================       o            o            o
      if( ( any( bndtypex(:) /= geom%bndtypex(:) ) .and..not.  &
            ( all( bndtypex(:) == ENUM%bndtype%reflect ) .and.  &
              any( geom%bndtypex(:) == ENUM%bndtype%reflect ) .and.  &
              any( geom%bndtypex(:) == ENUM%bndtype%periodic ) ) ) .or.  &
          ( any( bndtypey(:) /= geom%bndtypey(:) ) .and..not.  &
            ( all( bndtypey(:) == ENUM%bndtype%reflect ) .and.  &
              any( geom%bndtypey(:) == ENUM%bndtype%reflect ) .and.  &
              any( geom%bndtypey(:) == ENUM%bndtype%periodic ) ) ) .or.  &
          any( bndtypez(:) /= geom%bndtypez(:) ) ) then
         error stop 'Boundary types on octree inconsistent with those specified in INP'
      end if
!
!  Adjust bounding box according to boundary types
!
      wvacuum = sn%pmax * (ONE+TWO*EPSMB)
      if( geom%bndtypex(1) == ENUM%bndtype%free )  xlim(1) = xlim(1) - wvacuum
      if( geom%bndtypex(2) == ENUM%bndtype%free )  xlim(2) = xlim(2) + wvacuum
      if( geom%bndtypey(1) == ENUM%bndtype%free )  ylim(1) = ylim(1) - wvacuum
      if( geom%bndtypey(2) == ENUM%bndtype%free )  ylim(2) = ylim(2) + wvacuum
      if( geom%bndtypez(1) == ENUM%bndtype%free )  zlim(1) = zlim(1) - wvacuum
      if( geom%bndtypez(2) == ENUM%bndtype%free )  zlim(2) = zlim(2) + wvacuum
!
!  Write signed distance function to file (for debugging)
!
      if( out%ldist ) then

         wx = 10
         wy = 10
         wz = 10
         nboxx = (xlim(2) - xlim(1) + 40) / wx
         nboxy = (ylim(2) - ylim(1) + 40) / wy
         nboxz = (zlim(2) - zlim(1) + 40) / wz
         print *, 'nboxx, nboxy, nboxz = ', nboxx, nboxy, nboxz
         his%hbeg(:) = [ xlim(1)-20 + 0.5*wx, ylim(1)-20 + 0.5*wy, zlim(1)-20 + 0.5*wz ]

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

         his%nbox(:) = [ nboxx, nboxy, nboxz ]
         his%imin(:) = [ 1, 1, 1 ]
         his%imax(:) = [ nboxx, nboxy, nboxz ]
         his%wbox(:) = [ wx, wy, wz ]
         allocate( his%indhis( 3 ) )
         allocate( his%hist( 3, nboxx, nboxy, nboxz ) )

!         open( 33, file='../DIST' )
!         write( 33, * )  3, nboxz
!         write( 33, * )  0, 0, 0

         do i = 1, nboxx
            x = his%hbeg(1)+(i-1)*wx
!            x = 0
!            i = 1
            do j = 1, nboxy
               y = his%hbeg(2)+(i-1)*wy
!               y = 0
!               j = 1
               do k = 1, nboxz
                  z = his%hbeg(3)+(k-1)*wz
                  call octree%getDistSurface( x, y, z, his%hist(3,i,j,k), isOutside=outside )
!                      =====================  i  i  i         o                   o
                  if( outside ) then
                     his%hist(1,i,j,k) = 1
                  else
                     his%hist(1,i,j,k) = 0
                  end if
                  call octree%getRegionIndex(x, y, z, ir)
!                      ===================== i  i  i  o
                  his%hist(2,i,j,k) = ir
!                  write( 33, * )  z, ir, his%hist(3,i,j,k)
               end do
            end do
         end do
!            close( 33 )
!            error stop

         call writeHis3d( '../DIST3.his3', 'write', his )
!             ==========       i         i      i
         error stop

      end if

   end subroutine initSurface3dSDFOct
!
!=======================================================================
!
   subroutine getNcellSurface3dSDFOct( ncell )
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

   end subroutine getNcellSurface3dSDFOct
!
!=======================================================================
!
   subroutine getCellIndex3dSDFOct( x, y, z, ic )
!             ====================  i  i  i  o
!
!     Get the cell index at position [x,z]
!
!  Arguments:
!
      real(REAL8), intent(in) :: x               ! x coordinate of point
      real(REAL8), intent(in) :: y               ! y coordinate of point
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
      ic = octree%getLeafNodeIdAtCoords( x, y, z )
!          ============================  i  i  i

   end subroutine getCellIndex3dSDFOct
!
!=======================================================================
!
   subroutine getDistSurface3dSDFOct( x, y, z, sdf, dirx_, diry_, dirz_ )
!             ======================  i  i  i   o     o      o      o
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
      real(REAL8), intent(in) :: y                 ! y coordinate of point
      real(REAL8), intent(in) :: z                 ! z coordinate of point
      real(REAL8), intent(out) :: sdf              ! signed distance: (+) inside, (-) outside
      real(REAL8), optional, intent(out) :: dirx_  ! x/z component of direction towards/from
      real(REAL8), optional, intent(out) :: diry_  !    surface
      real(REAL8), optional, intent(out) :: dirz_  !
!
!  Local Variables:
!
!     None
!
!  Start of Executable Code:
!  ------------------------
!
      call octree%getDistSurface( x, y, z, sdf, dirx_, diry_, dirz_ )
!          =====================  i  i  i   o     o      o      o

   end subroutine getDistSurface3dSDFOct
!
!=======================================================================
!
end module Target_Surface3d_SDFOct_mod
