module Target_Surface3d_mod
!
!  Define the surface and provide surface related operations
!
!  Used Modules:
!  ============
!
   use Const_mod,                   only: INT4, INT8, REAL8,  &
                                          ZERO, HALF, ONE, TWO, SMALL, BIG, EPSMB,  &
                                          ENUM
   use Types_mod,                   only: Proj_type
   use Undefine_mod,                only: defined

   use Param_Geom_mod,              only: geom
   use Param_Sn_mod,                only: sn

   use Geom_BndCond_mod,            only: applyBndCond3d

   use Target_Surface3d_SDF_mod,    only: initSurface3dSDF, getNcellSurface3dSDF,  &
                                          getDistSurface3dSDF, getCellIndex3dSDF
#ifdef GTS
   use Target_Surface3d_SDFOct_mod, only: initSurface3dSDFOct, getNcellSurface3dSDFOct,  &
                                          getDistSurface3dSDFOct, getCellIndex3dSDFOct
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
   public :: initSurface3d                   ! initialize surface
!            =============
   public :: getNcellSurface3d               ! get the number of cells
!            =================
   public :: getCellIndexSurface3d                  ! get cell index at given position
!            ==============
   public :: getDistSurface3d                ! get the signed distance of a point
!            ================                !    from the surface
   public :: getEntranceGrid3d               ! get the entrance point into the grid
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
   type :: Domain_type                       ! the domain on which the grid is defined
!          ===========
      real(REAL8) :: xentr(2)                ! x coordinate of possible entrance point
      real(REAL8) :: yentr(2)                ! y coordinate of possible entrance point
      real(REAL8) :: zentr(2)                ! z coordinate of possible entrance point
   end type Domain_type

   type(Domain_type) :: domain
!
!=======================================================================
!
contains
!
!=======================================================================
!
   subroutine initSurface3d( )
!             =============
!
!     Initialize the simulation domain and the signed distance function from the surface.
!     This subroutine sets/modifies geom%xlim, geom%ylim, geom%zlim
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
!  Initialize SDF if SDF is given on files
!
      if( geom%usehyperoctree ) then

#ifdef GTS
         call initSurface3dSDFOct( geom%xlim(:), geom%ylim(:), geom%zlim(:) )
!             ===================       o             o             o
#else
         error stop 'initSurface3d: usehyperoctree=t not supported'
#endif

      else if( geom%usesdf ) then

         call initSurface3dSDF( geom%xlim(:), geom%ylim(:), geom%zlim(:) )
!             ================       o             o             o

      end if
!
!  Adjust geom%xlim(:), geom%ylim(:), and geom%zlim(:) at reflective boundaries
!
      if( geom%bndtypex(1) == ENUM%bndtype%reflect ) then
         geom%xlim(3) = geom%xlim(1)
         geom%xlim(1) = TWO*geom%xlim(1) - geom%xlim(2)
      else if( geom%bndtypex(2) == ENUM%bndtype%reflect ) then
         geom%xlim(3) = geom%xlim(2)
         geom%xlim(2) = TWO*geom%xlim(2) - geom%xlim(1)
      end if

      if( geom%bndtypey(1) == ENUM%bndtype%reflect ) then
         geom%ylim(3) = geom%ylim(1)
         geom%ylim(1) = TWO*geom%ylim(1) - geom%ylim(2)
      else if( geom%bndtypey(2) == ENUM%bndtype%reflect ) then
         geom%ylim(3) = geom%ylim(2)
         geom%ylim(2) = TWO*geom%ylim(2) - geom%ylim(1)
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

   end subroutine initSurface3d
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

      if( geom%bndtypey(1) == ENUM%bndtype%free ) then
         domain%yentr(1) = geom%ylim(1) + sn%pmax*EPSMB
      else if( geom%bndtypey(1) == ENUM%bndtype%reflect .and.  &
               geom%bndtypey(2) == ENUM%bndtype%free ) then
         domain%yentr(1) = TWO*geom%ylim(1) - geom%ylim(2) + sn%pmax*EPSMB
      else
         domain%yentr(1) = -BIG
      end if

      if( geom%bndtypey(2) == ENUM%bndtype%free ) then
         domain%yentr(2) = geom%ylim(2) - sn%pmax*EPSMB
      else if( geom%bndtypey(2) == ENUM%bndtype%reflect .and.  &
               geom%bndtypey(1) == ENUM%bndtype%free ) then
         domain%yentr(2) = TWO*geom%ylim(2) - geom%ylim(1) - sn%pmax*EPSMB
      else
         domain%yentr(2) = BIG
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

   end subroutine setDomainBoundaries
!
!=======================================================================
!
   subroutine getNcellSurface3d( ncell )
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
         call getNcellSurface3dSDFOct( ncell )
!             =======================    o
#else
         error stop 'getNcellSurface3d: usehyperoctree=t not supported'
#endif

      else

         call getNcellSurface3dSDF( ncell )
!             ====================    o

      end if

   end subroutine getNcellSurface3d
!
!=======================================================================
!
   subroutine getCellIndexSurface3d( x, y, z, ic, outside, what_ )
!             =====================  i  i  i  o      o       i
!
!     Determine the cell index of point [x,z]. x and z are the coordinates in the geometry
!     coordinate system
!
!     If point is outside the grid, the nearest cell is chosen
!
!  Arguments:
!
      real(REAL8), intent(in) :: x              ! x coordinate of point
      real(REAL8), intent(in) :: y              ! y coordinate of point
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
         call getCellIndex3dSDFOct( x, y, z, ic )
!             ====================  i  i  i  o
#else
         error stop 'getCellIndexSurface3d: usehyperoctree=t not supported'
#endif
         outside = ( x < geom%xlim(1) .or. x > geom%xlim(2) .or.  &
                     y < geom%ylim(1) .or. y > geom%ylim(2) .or.  &
                     z < geom%zlim(1) .or. z > geom%zlim(2) )

         return

      else

         call getCellIndex3dSDF( x, y, z, ic, outside, what_ )
!             =================  i  i  i  o      o       i

      end if

   end subroutine getCellIndexSurface3d
!
!=======================================================================
!
   subroutine getDistSurface3d( pos, sdf, outside_, dir_ )
!             ================   i    o      o       o
!
!     Get the signed distance of a point pos=(x,y,z) from the surface by interpolation
!     ( + inside, - outside ) and
!     the direction towards/from the surface (pointing outwards the target)
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
      real(REAL8) :: x, y, z, dirx, diry, dirz, length
!
!  Start of Executable Code:
!  ------------------------
!
      x = pos(1)
      y = pos(2)
      z = pos(3)
      call applyBndCond3d( x, y, z )
!          ============== io io io
!
!  Check if point is outside grid
!
      if( x < geom%xlim(1) .or. x > geom%xlim(2) .or.  &
          y < geom%ylim(1) .or. y > geom%ylim(2) .or.  &
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
            call getDistSurface3dSDFOct( x, y, z, sdf, dirx, diry, dirz )
!                ======================  i  i  i   o    o     o     o
         else
            call getDistSurface3dSDFOct( x, y, z, sdf )
!                ======================  i  i  i   o
         end if
#else
         error stop 'getDistSurface3d: usehyperoctree=t not supported'
#endif

      else

         if( present( dir_ ) ) then
            call getDistSurface3dSDF( x, y, z, sdf, dirx, diry, dirz )
!                ===================  i  i  i   o    o     o     o
         else
            call getDistSurface3dSDF( x, y, z, sdf )
!                ===================  i  i  i   o
         end if

      end if
!
!  Normalize direction towards/from surface ("-" in order that direction points outwards)
!
      if( present( dir_ ) ) then

         length = sqrt( dirx**2 + diry**2 + dirz**2 )
         dir_(:) = [ -dirx, -diry, -dirz ] / length

      end if

   end subroutine getDistSurface3d
!
!=======================================================================
!
   subroutine getEntranceGrid3d( proj, found, onlyforward )
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
      integer(INT4) :: i
      real(REAL8) :: x, y, z, dirx, diry, dirz
      real(REAL8) :: xx, yy, zz, rlam, lambda
!
!  Start of Executable Code:
!  ------------------------
!
!  Determine relative coordinates
!
      x = proj%pos(1)
      y = proj%pos(2)
      z = proj%pos(3)
      dirx = proj%dir(1)
      diry = proj%dir(2)
      dirz = proj%dir(3)
!
!  If onlyforward=t, check if position is inside grid. If so, return
!
      if( onlyforward ) then
         found = ( x >= domain%xentr(1) .and. x <= domain%xentr(2) .and.  &
                   y >= domain%yentr(1) .and. y <= domain%yentr(2) .and.  &
                   z >= domain%zentr(1) .and. z <= domain%zentr(2) )
         if( found )  return
      end if
!
!  Look for entrance point
!
      lambda = BIG

      if( .not.any( geom%bndtypex(:) == ENUM%bndtype%free ) .and.  &
          .not.any( geom%bndtypey(:) == ENUM%bndtype%free ) ) then
                                                ! more efficient version
         if( dirz > SMALL ) then
            lambda = (domain%zentr(1) - z) / dirz
         else if( dirz < -SMALL ) then
            lambda = (domain%zentr(2) - z) / dirz
         end if

      else

         if( abs( dirz ) > SMALL ) then      ! top/bottom boundary
            do i = 1, 2
               rlam = (domain%zentr(i) - z) / dirz
               xx = x + rlam * dirx
               yy = y + rlam * diry
               if( xx >= domain%xentr(1) .and. xx <= domain%xentr(2) .and.  &
                   yy >= domain%yentr(1) .and. yy <= domain%yentr(2) .and.  &
                   (.not.onlyforward .or. rlam > ZERO) )  lambda = min( lambda, rlam )
            end do
         end if

         if( abs( diry ) > SMALL ) then      ! front/back boundary
            do i = 1, 2
               rlam = (domain%yentr(i) - y) / diry
               xx = x + rlam * dirx
               zz = z + rlam * dirz
               if( xx >= domain%xentr(1) .and. xx <= domain%xentr(2) .and.  &
                   zz >= domain%zentr(1) .and. zz <= domain%zentr(2) .and.  &
                   (.not.onlyforward .or. rlam > ZERO) )  lambda = min( lambda, rlam )
            end do
         end if

         if( abs( dirx ) > SMALL ) then      ! left/right boundary
            do i = 1, 2
               rlam = (domain%xentr(i) - x) / dirx
               yy = y + rlam * diry
               zz = z + rlam * dirz
               if( yy >= domain%yentr(1) .and. yy <= domain%yentr(2) .and.  &
                   zz >= domain%zentr(1) .and. zz <= domain%zentr(2) .and.  &
                   (.not.onlyforward .or. rlam > ZERO) )  lambda = min( lambda, rlam )
            end do
         end if

      end if
!
!  Set new position
!
      if( lambda == BIG ) then
         found = .false.
      else
         found = .true.
         proj%posold(:) = proj%pos(:)
         proj%pos(:) = proj%pos(:) + lambda * proj%dir(:)
      end if
      !print *, 'Entrance position =', proj%pos(:)

   end subroutine getEntranceGrid3d
!
!=======================================================================
!
end module Target_Surface3d_mod
