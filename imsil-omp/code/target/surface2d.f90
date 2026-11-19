module Target_Surface2d_mod
!
!  Define the surface and provide surface related operations
!
!  Used Modules:
!  ============
!
   use Const_mod,                   only: INT4, INT8, REAL8,  &
                                          ZERO, ONE, TWO, SMALL, BIG, EPSMB,  &
                                          ENUM
   use Undefine_mod,                only: defined
   use Types_mod,                   only: Proj_type

   use Param_Geom_mod,              only: geom
   use Param_Sn_mod,                only: sn
   use Param_Out_mod,               only: out

   use Geom_BndCond_mod,            only: applyBndCond2d
   use Geom_CurviCoord_mod,         only: toCurviCoord2d, fromCurviZDir, fromGeomNormDir2d
   use Geom_Polygon_Old_mod,        only: writePolygons

   use Target_Surface2d_SDF_mod,    only: initSurface2dSDF, initSurface2dSDFfromSegments,  &
                                          getNcellSurface2dSDF, getDistSurface2dSDF,  &
                                          getCellIndex2dSDF
#ifdef GTS
   use Target_Surface2d_SDFOct_mod, only: initSurface2dSDFOct, getNcellSDFOct,  &
                                          getDistSurface2dSDFOct, getCellIndex2dSDFOct
#endif
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
   public :: initSurface2d                   ! initialize surface
!            =============
   public :: initSurface2dFromSegments       ! initialize surface from segments
!            =========================
   public :: getNcellSurface2d               ! get the number of cells
!            =================
   public :: getCellIndexSurface2d           ! get cell index at given position
!            ==============
   public :: getDistSurface2d                ! get the signed distance of a point
!            ================                !    from the surface
   public :: getEntranceGrid2d               ! get the entrance point into the grid
!            =================               !    covering the target
!
!  Private Procedures:
!  ==================
!
   private :: setDomainBoundaries            ! set domain boundaries
!             ===================
!
!  Private Data:
!  ============
!
   type :: Domain_type                       ! the domain where particles can move
!          ===========
      real(REAL8) :: xentr(2)                ! x coordinate of possible entrance point
      real(REAL8) :: zentr(2)                ! z coordinate of possible entrance point
      real(REAL8) :: rentr2(2)               ! squares of radii of possible entrance points
   end type Domain_type

   type(Domain_type) :: domain
!
!=======================================================================
!
contains
!
!=======================================================================
!
   subroutine initSurface2d( )
!             =============
!
!     Initialize the simulation domain and the signed distance from the surface.
!
!  Arguments:
!
!     None
!
!  Local Variables:
!
!     None
!
!  Start of Executable Code:
!  ------------------------
!
!  Read SDF from file
!
      if( geom%usehyperoctree ) then

#ifdef GTS
         call initSurface2dSDFOct( geom%xlim(:), geom%zlim(:) )
!             ===================       o             o
#else
         error stop 'initSurface2d: usehyperoctree=t not supported'
#endif
!
!  Read SDF from file or initialize from regions
!
      else

         call initSurface2dSDF( )
!             ================

      end if
!
!  Adjust geom%xlim(:) and geom%zlim(:) at reflective boundaries (this is done only at init,
!  not at an update using initSurface2dFromSegments, because for ldyn=t reflective boundaries
!  have been eliminated; in contrast, at init the geometry is reflected at a reflective
!  boundary)
!
      if( geom%bndtypex(1) == ENUM%bndtype%reflect ) then
         geom%xlim(3) = geom%xlim(1)
         geom%xlim(1) = TWO*geom%xlim(1) - geom%xlim(2)
      else if( geom%bndtypex(2) == ENUM%bndtype%reflect ) then
         geom%xlim(3) = geom%xlim(2)
         geom%xlim(2) = TWO*geom%xlim(2) - geom%xlim(1)
      end if

      if( geom%bndtypez(1) == ENUM%bndtype%reflect ) then
         geom%zlim(3) = geom%zlim(1)
         geom%zlim(1) = TWO*geom%zlim(1) - geom%zlim(2)
      else if( geom%bndtypez(2) == ENUM%bndtype%reflect ) then
         geom%zlim(3) = geom%zlim(2)
         geom%zlim(2) = TWO*geom%zlim(2) - geom%zlim(1)
      end if
!
!  Calculate domain boundaries
!
      call setDomainBoundaries( )
!          ===================

   end subroutine initSurface2d
!
!=======================================================================
!
   subroutine initSurface2dFromSegments( x, z, ffgeom, segments )
!             =========================  i  i    i        i
!
!     Initialize the signed distance function from the cell fill fractions and a set of
!     surface segments. The segments need not be connected. The fill fractions are needed
!     to determine the sign of the SDF.
!
!     The rectilinear grid must cover all SDF > -sn%pmax*(1+2*epsB)
!
      real(REAL8), intent(in) :: x(0:)             ! (0:nx) x coordinates of rectilinear grid
                                                   !    for distance values
      real(REAL8), intent(in) :: z(0:)             ! (0:nz) z coordinates of rectilinear grid
                                                   !    for distance values
      real(REAL8), intent(in) :: ffgeom(:,:)       ! (nx,nz) fill fractions
      real(REAL8), target, intent(in) :: segments(:,:,:) ! (2,2,ns) surface segments
!
!  Local Variables
!
!     None
!
!  Start of Executable Code:
!  ------------------------
!
      call initSurface2dSDFfromSegments( x(:), z(:), ffgeom(:,:), segments(:,:,:) )
!          ============================   i     i        i              i
!
!  Calculate domain boundaries
!
      call setDomainBoundaries( )
!          ===================

   end subroutine initSurface2dFromSegments
!
!=======================================================================
!
   subroutine setDomainBoundaries( )
!             ===================
!
!     Set the boundaries of the simulation domain depending on boundary types
!
!  Arguments:
!
!     None
!
!  Local Variables:
!
!     None
!
!  Start of Executable Code:
!  ------------------------
!
!  Slightly inside boundaries of simulation domain
!
      if( geom%bndtypex(1) == ENUM%bndtype%free ) then
         domain%xentr(1) = geom%xlim(1) + sn%pmax*EPSMB
      else if( geom%bndtypex(1) == ENUM%bndtype%reflect .and.  &
               geom%bndtypex(2) == ENUM%bndtype%free ) then
         domain%xentr(1) = TWO*geom%xlim(1) - geom%xlim(2) + sn%pmax*EPSMB
      else
         domain%xentr(1) = -BIG
      end if

      if( geom%bndtypex(2) == ENUM%bndtype%free ) then
         domain%xentr(2) = geom%xlim(2) - sn%pmax*EPSMB
      else if( geom%bndtypex(2) == ENUM%bndtype%reflect .and.  &
               geom%bndtypex(1) == ENUM%bndtype%free ) then
         domain%xentr(2) = TWO*geom%xlim(2) - geom%xlim(1) - sn%pmax*EPSMB
      else
         domain%xentr(2) = BIG
      end if

      if( geom%bndtypez(1) == ENUM%bndtype%free ) then
         domain%zentr(1) = geom%zlim(1) + sn%pmax*EPSMB
      else if( geom%bndtypez(1) == ENUM%bndtype%reflect .and.  &
               geom%bndtypez(2) == ENUM%bndtype%free ) then
         domain%zentr(1) = TWO*geom%zlim(1) - geom%zlim(2) + sn%pmax*EPSMB
      else
         domain%zentr(1) = -BIG
      end if

      if( geom%bndtypez(2) == ENUM%bndtype%free ) then
         domain%zentr(2) = geom%zlim(2) - sn%pmax*EPSMB
      else if( geom%bndtypez(2) == ENUM%bndtype%reflect .and.  &
               geom%bndtypez(1) == ENUM%bndtype%free ) then
         domain%zentr(2) = TWO*geom%zlim(2) - geom%zlim(1) - sn%pmax*EPSMB
      else
         domain%zentr(2) = BIG
      end if

      if( .not.geom%lcart ) then
         domain%rentr2(1) = max( domain%zentr(1), ZERO )**2
         domain%rentr2(2) = domain%zentr(2)**2
      end if

   end subroutine setDomainBoundaries
!
!=======================================================================
!
   subroutine getNcellSurface2d( ncell )
!             =================    o
!
!     Get the number of cells
!
!  Arguments:
!
      integer(INT8), intent(out) :: ncell       ! number of cells
!
!  Local Variables:
!
!     None
!
!  Start of Executable Code:
!  ------------------------
!
      if( geom%usehyperoctree ) then

#ifdef GTS
         call getNcellSurface2dSDFOct( ncell )
!             =======================    o
#else
         error stop 'getNcellSurface2d: usehyperoctree=t not supported'
#endif

      else

         call getNcellSurface2dSDF( ncell )
!             ====================    o

      end if

   end subroutine getNcellSurface2d
!
!=======================================================================
!
   subroutine getCellIndexSurface2d( x, z, ic, outside, what_ )
!             =====================  i  i  o      o       i
!
!     Determine the cell index of point [x,z]. x and z are the coordinates in the geometry
!     coordinate system
!
!     If point is outside the grid, the nearest cell is chosen
!
!  Arguments:
!
      real(REAL8), intent(in) :: x              ! x coordinate of point
      real(REAL8), intent(in) :: z              ! z coordinate of point
      integer(INT8), intent(out) :: ic          ! cell index
      logical, intent(out) :: outside           ! indicates point is outside grid
      character(*), intent(in), optional :: what_  ! indication what point is checked
!
!  Local Variables:
!
!     None
!
!  Start of Executable Code:
!  ------------------------
!
      if( geom%usehyperoctree ) then

#ifdef GTS
         call getCellIndex2dSDFOct( x, z, ic )
!             ====================  i  i  o
#else
         error stop 'getCellIndexSurface2d: usehyperoctree=t not supported'
#endif
         outside = ( x < geom%xlim(1) .or. x > geom%xlim(2) .or.  &
                     z < geom%zlim(1) .or. z > geom%zlim(2) )

         return

      else

         call getCellIndex2dSDF( x, z, ic, outside, what_ )
!             =================  i  i  o      o       i

      end if

   end subroutine getCellIndexSurface2d
!
!=======================================================================
!
   subroutine getDistSurface2d( pos, sdf, outside_, dir_ )
!             ================   i    o      o       o
!
!     Get
!     - the signed distance of a point from the surface by interpolation ( + inside, - outside )
!     - the direction towards/from the surface (pointing outwards the target)
!
!     NOTE: If the point is outside the grid, dist is set to -BIG, and dir_ is undefined.
!
!  Arguments:
!
      real(REAL8), intent(in) :: pos(3)            ! position of point
      real(REAL8), intent(out) :: sdf              ! signed distance: (+) inside, (-) outside
      logical, optional, intent(out) :: outside_   ! indicates point is outside grid
      real(REAL8), optional, intent(out) :: dir_(3)! direction towards/from surface
!
!  Local Variables:
!
      real(REAL8) :: x, z, dirx, dirz, length, zdir(3)
!
!  Start of Executable Code:
!  ------------------------
!
      if( geom%lcart ) then
         x = pos(1)
         z = pos(3)
      else
         call toCurviCoord2d( pos(:), x, z )
!             ==============    i     o  o
      end if
      call applyBndCond2d( x, z )
!          ==============
!
!  Check if point is outside grid
!
      if( x < geom%xlim(1) .or. x > geom%xlim(2) .or.  &
          z < geom%zlim(1) .or. z > geom%zlim(2) ) then
         sdf = -BIG
         if( present( outside_ ) )  outside_ = .true.
         return
      end if

      if( present( outside_ ) )  outside_ = .false.
!
!  Get SDF and direction towards/from surface
!
      if( geom%usehyperoctree ) then

#ifdef GTS
         if( present( dir_ ) ) then
            call getDistSurface2dSDFOct( x, z, sdf, dirx, dirz )
!                ======================  i  i   o    o     o
         else
            call getDistSurface2dSDFOct( x, z, sdf )
!                ======================  i  i   o
         end if
#else
         error stop 'getDistSurface2d: usehyperoctree=t not supported'
#endif

      else

         if( present( dir_ ) ) then
            call getDistSurface2dSDF( x, z, sdf, dirx, dirz )
!                ===================  i  i   o    o     o
         else
            call getDistSurface2dSDF( x, z, sdf )
!                ===================  i  i   o
         end if

      end if
!
!  Normalize direction towards/from surface
!
      if( present( dir_ ) ) then

         length = sqrt( dirx**2 + dirz**2 )
         if( length < SMALL ) then
            dirx = ZERO
            dirz = - ONE   ! *sign(dirz) !
         else
            dirx = - dirx / length              ! "-" in order that direction points outwards
            dirz = - dirz / length
         end if

         if( geom%lcart ) then
            dir_(:) = [ dirx, ZERO, dirz ]
         else
            call fromCurviZDir( pos(:), zdir(:) )
!                =============    i       o
            call fromGeomNormDir2d( zdir(:), dirx, dirz, dir_(:) )
!                =================    i       i     i      o
         end if

      end if

   end subroutine getDistSurface2d
!
!=======================================================================
!
   subroutine getEntranceGrid2d( proj, found, onlyforward )
!             =================   io     o         i
!
!     Determine the entrance point of a projectile into the grid covering the target
!     The entrance point may also lie backwards of the projectile if onlyforward=f.
!     The entrance point is identical to the current position if the position is inside the
!     grid and onlyforward=t.
!
!  Arguments:
!
      type(Proj_type), intent(inout) :: proj ! projectile
      logical, intent(out) :: found          ! flag indicating whether proj%pos contains
                                             !    a valid position
      logical, intent(in) :: onlyforward     ! flag indicating that entrance point may only
                                             !    lie in the forward direction
!
!  Local Variables:
!
      integer(INT4) :: i, j
      real(REAL8) :: x, y, z, dirx, diry, dirz
      real(REAL8) :: r2, plusminus(2)
      real(REAL8) :: a, b, c, discr, sqrtd, xx, yy, zz, lambda, lambdamin
!
!  Start of Executable Code:
!  ------------------------
!
!  Determine relative coordinates
!  In case of cylindrical symmetry x is in the direction of the axis, y and z are perpendicular
!  to it, so that sqrt(y**2+z**2) is the distance from the axis (=radius)
!
      if( geom%lcart ) then
         x = proj%pos(1)
         z = proj%pos(3)
         dirx = proj%dir(1)
         dirz = proj%dir(3)
      else if( geom%lcylx ) then
         x = proj%pos(1) - geom%xcen
         y = proj%pos(2) - geom%ycen
         z = proj%pos(3) - geom%zcen
         dirx = proj%dir(1)
         diry = proj%dir(2)
         dirz = proj%dir(3)
      else if( geom%lcyly ) then
         x = proj%pos(2) - geom%ycen
         y = proj%pos(1) - geom%xcen
         z = proj%pos(3) - geom%zcen
         dirx = proj%dir(2)
         diry = proj%dir(1)
         dirz = proj%dir(3)
      else if( geom%lcylz ) then
         x = proj%pos(3) - geom%zcen
         y = proj%pos(1) - geom%xcen
         z = proj%pos(2) - geom%ycen
         dirx = proj%dir(3)
         diry = proj%dir(1)
         dirz = proj%dir(2)
      else
         error stop 'getEntranceGrid2d: not all of lcart, lcylx, lcyly, lcylz must be false'
      end if
!
!  If onlyforward=t, check if position is inside target. If so, return
!
      if( onlyforward ) then
         if( geom%lcart ) then
            found = ( x >= domain%xentr(1) .and. x <= domain%xentr(2) .and.  &
                      z >= domain%zentr(1) .and. z <= domain%zentr(2) )
         else
            r2 = y**2 + z**2
            found = ( x >= domain%xentr(1) .and. x <= domain%xentr(2) .and.  &
                      r2 >= domain%rentr2(1) .and. r2 <= domain%rentr2(2) )
         end if
         if( found )  return
      end if
!
!  Look for entrance point
!
      lambdamin = BIG

      if( geom%lcart ) then

         if( .not.any( geom%bndtypex(:) == ENUM%bndtype%free ) ) then
                                                ! more efficient version
            if( dirz > SMALL ) then
               lambdamin = (domain%zentr(1) - z) / dirz
            else if( dirz < -SMALL ) then
               lambdamin = (domain%zentr(2) - z) / dirz
            end if

         else

            if( abs( dirz ) > SMALL ) then      ! top/bottom boundary
               do i = 1, 2
                  lambda = (domain%zentr(i) - z) / dirz
                  xx = x + lambda * dirx
                  if( xx >= domain%xentr(1) .and. xx <= domain%xentr(2) .and.  &
                      (.not.onlyforward .or. lambda > ZERO) )  &
                        lambdamin = min( lambdamin, lambda )
               end do
            end if

            if( abs( dirx ) > SMALL ) then      ! left/right boundary
               do i = 1, 2
                  lambda = (domain%xentr(i) - x) / dirx
                  zz = z + lambda * dirz
                  if( zz >= domain%zentr(1) .and. zz <= domain%zentr(2) .and.  &
                      (.not.onlyforward .or. lambda > ZERO) )  &
                        lambdamin = min( lambdamin, lambda )
               end do
            end if

         end if

      else

         if( abs( dirx ) > SMALL ) then         ! left/right boundary
            do i = 1, 2
               lambda = (domain%xentr(i) - x) / dirx
               yy = y + lambda * diry
               zz = z + lambda * dirz
               r2 = yy**2 + zz**2
               if( r2 >= domain%rentr2(1) .and. r2 <= domain%rentr2(2) .and.  &
                   (.not.onlyforward .or. lambda > ZERO) )  &
                     lambdamin = min( lambdamin, lambda )
            end do
         end if

         a = diry**2 + dirz**2                  ! inner and outer cylinder surface
         b = 2 * (y*diry + z*dirz)
         do i = 1, 2
            c = y**2 + z**2 - domain%rentr2(i)
            discr = b**2 - 4*a*c
            if( discr >= ZERO ) then
               sqrtd = sqrt(discr)
               plusminus(:) = [ +ONE, -ONE ]    ! first and second solution
               do j = 1, 2
                  lambda = (-b + plusminus(i)*sqrtd) / (2*a)
                  xx = x + lambda * dirx
                  if( xx >= domain%xentr(1) .and. xx <= domain%xentr(2) .and.  &
                      (.not.onlyforward .or. lambda > ZERO) )  &
                        lambdamin = min( lambdamin, lambda )
               end do
            end if
         end do

      end if
!
!  Set new position
!
      if( lambdamin == BIG ) then
         found = .false.
      else
         found = .true.
         proj%posold(:) = proj%pos(:)
         proj%pos(:) = proj%pos(:) + lambdamin * proj%dir(:)
      end if
      !print *, 'Entrance position =', proj%pos(:)

   end subroutine getEntranceGrid2d
!
!=======================================================================
!
end module Target_Surface2d_mod
