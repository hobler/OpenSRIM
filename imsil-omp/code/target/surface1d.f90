module Target_Surface1d_mod
!
!  Define the surface by two points along the axis
!
!  Used Modules:
!  ============
!
   use Const_mod,             only: INT4, REAL8,  &
                                    ZERO, HALF, ONE, TWO, SMALL, BIG, EPSMB,  &
                                    ENUM
   use Types_mod,             only: Proj_type

   use Param_Setup_mod,       only: ldyn, usecell, cellrec, cellfile, nelem, nelemIon
   use Param_Geom_mod,        only: geom
   use Param_Sn_mod,          only: sn

   use Geom_BndCond_mod,      only: applyBndCond1d
   use Geom_CurviCoord_mod,   only: toCurviCoord1d, fromCurviZDir
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
   public :: initSurface1d                   ! initialize the surface
!            =============                   !
   public :: updateSurface1d                 ! update the surface using new surface positions
!            ===============                 !
   public :: getDistSurface1d                ! get the signed distance of a point
!            ================                !    from the surface
   public :: getEntranceGrid1d               ! get the entrance point into the grid
!            =================               !    covering the target
   public :: getExitPoint1d                  ! get the exit point from the target (at a 
!            ==============                  !    distance sn%pmax from the surface)
!
!  Private Procedures:
!  ==================
!
   private :: isInsideGrid                   ! determine if position is inside the grid
!             ============
   private :: intersect                      ! intersect a line with a z-plane, a cylinder
!             =========                      !    or a sphere
!
!  Private Data:
!  ============
!
   real(REAL8) :: zsurf(2)                   ! surface positions
                                             ! NOTE: by increasing the dimension, one could
                                             !       define holes in the target
   real(REAL8) :: zrefract(2)                ! positions of refraction planes (sn%pmax outside
                                             !    zsurf)
   real(REAL8) :: zrefract1(2)               ! a little inside zrefract
   real(REAL8) :: zrefract2(2)               ! a little outside zrefract
   ! note that geom%zlim(:) is a little outside zrefract2
!
!=======================================================================
!
contains
!
!=======================================================================
!
   subroutine initSurface1d( )
!             =============
!
!     Initialize the surfaces
!
!     geom%zlim(:) must be the minimum and maximum coordinates of the target,
!     when the procedure is entered. It is modified at reflective and free
!     boundaries.
!
!  Arguments:
!
!
!  Local Variables:
!
!     None
!
!  Start of Executable Code:
!  ------------------------
!
!  Adjust geom%zlim(:) at reflective boundaries (this is done only at init, not at update,
!  because for ldyn=t reflective boundaries have been eliminated; in contrast, at init the
!  geometry is reflected at a reflective boundary)
!
      if( geom%bndtypez(1) == ENUM%bndtype%reflect ) then
         geom%zlim(3) = geom%zlim(1)
         geom%zlim(1) = TWO*geom%zlim(1) - geom%zlim(2)
      else if( geom%bndtypez(2) == ENUM%bndtype%reflect ) then
         geom%zlim(3) = geom%zlim(2)
         geom%zlim(2) = TWO*geom%zlim(2) - geom%zlim(1)
      end if
!
!  Define the surfaces
!
      call updateSurface1d( geom%zlim(1), geom%zlim(2) )
!          ===============       i             i

   end subroutine initSurface1d
!
!=======================================================================
!
   subroutine updateSurface1d( zmin, zmax )
!             ===============   i     i
!
!     Initialize the signed distance function from the positions of the surface.
!
!     Adds vacuum layers to free boundaries
!
!  Arguments:
!
      real(REAL8), intent(in) :: zmin           ! position of top surface
      real(REAL8), intent(in) :: zmax           ! position of bottom surface
!
!  Local Variables:
!
      real(REAL8) :: dz(2)
!
!  Start of Executable Code:
!  ------------------------
!
      zsurf(:) = [ zmin, zmax ]

      ! the following will be used only at free boundaries:
      zrefract(:) = [ zmin-sn%pmax, zmax+sn%pmax ]
      dz(:) = [ - max( abs(zmin), abs(sn%pmax) ), max( abs(zmax), abs(sn%pmax) ) ] * EPSMB
      zrefract1(:) = zrefract(:) - dz(:)
      zrefract2(:) = zrefract(:) + dz(:)

      where( geom%bndtypez(:) == ENUM%bndtype%free )
         geom%zlim(1:2) = zrefract(:) + TWO*dz(:)
      else where
         geom%zlim(1:2) = zsurf(:)
      end where

   end subroutine updateSurface1d
!
!=======================================================================
!
   subroutine getDistSurface1d( pos, dist, outside_, dir_ )
!             ================   i    o       o       o
!
!     Get the signed distance of a point from the surface and orientation of outward pointing
!     direction
!
!     If pos is outside (zlim(1),zlim(2)), then dist=-BIG is returned
!
!  Arguments:
!
      real(REAL8), intent(in) :: pos(3)            ! position of point
      real(REAL8), intent(out) :: dist             ! signed distance: (+) inside, (-) outside
      logical, optional, intent(out) :: outside_   ! indicates point is outside grid
      real(REAL8), optional, intent(out) :: dir_(3)! outwards direction
!
!  Local Variables:
!
      real(REAL8) :: z, dirfac, zdir(3)
!
!  Start of Executable Code:
!  ------------------------
!
      if( geom%lcart ) then
         z = pos(3)
      else
         call toCurviCoord1d( pos(:), z )
!             ==============    i     o
      end if
      call applyBndCond1d( z )
!          ============== io

      if( z < HALF * ( zsurf(1) + zsurf(2) ) ) then
         dist = z - zsurf(1)
         if( present( dir_ ) )  dirfac = -1
      else
         dist = zsurf(2) - z
         if( present( dir_ ) )  dirfac = +1
      end if

      if( z < geom%zlim(1) .or. z > geom%zlim(2) ) then
         dist = -BIG
         if( present( outside_ ) )  outside_ = .true.
      else
         if( present( outside_ ) )  outside_ = .false.
      end if

      if( present( dir_ ) ) then
         if( geom%lcart ) then
            dir_(:) = dirfac * [ ZERO, ZERO, ONE ]
         else
            call fromCurviZDir( pos(:), zdir(:) )
!                =============    i       o
            dir_(:) = dirfac * zdir(:)
         end if
      end if

   end subroutine getDistSurface1d
!
!=======================================================================
!
   subroutine getEntranceGrid1d( proj, found, onlyforward )
!             =================   io     o         i
!
!     Determine the entrance point of a projectile into the grid covering the target
!     The entrance point may also lie backwards of the projectile if onlyforward=f.
!
!     NOTE: The projectile is put slightly inside the grid
!
!  Arguments:
!
      type(Proj_type), intent(inout) :: proj ! projectile
      logical, intent(out) :: found          ! flag indicating whether entrance point
                                             !    has been found
      logical, intent(in) :: onlyforward     ! flag indicating that entrance point may only
                                             !    lie in the forward direction
!
!  Local Variables:
!
      integer(INT4) :: isurf
      real(REAL8) :: rlam, lambda
!
!  Start of Executable Code:
!  ------------------------
!
!  If projectile is already inside the grid and onlyforward=t, nothing to do
!
      if( isInsideGrid( proj%pos(:) ) .and. onlyforward ) then
!         ============      i
         found = .true.
         return
      end if
!
!  Check intersections with top and bottom surface
!
      lambda = BIG
!
      do isurf = 1, 2
         call intersect( proj%pos(:), proj%dir(:), zrefract2(isurf), onlyforward, rlam, found )
!             =========      i            i             i              i        o      o
         if( found )  lambda = min( lambda, rlam )
      end do
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

   end subroutine getEntranceGrid1d
!
!=======================================================================
!
   subroutine getExitPoint1d( pos1, pos2, found )
!             ==============   io    io     o
!
!     Determine the position on the line pos1-pos2 where the distance from the target
!     is sn%pmax. pos1 is assumed to be closer to the target than sn%pmax, pos2 is
!     further away from the target than sn%pmax
!
!     On exit pos1 and pos2 contain approximations to a position sn%pmax away from the 
!     target (they bracket the exact solution)
!
!  Arguments:
!
      real(REAL8), intent(inout):: pos1(3)      ! position closer than sn%pmax
      real(REAL8), intent(inout):: pos2(3)      ! position further away than sn%pmax
      logical, intent(out) :: found             ! flag indicating success
!
!  Local Variables:
!
      integer(INT4) :: isurf
      real(REAL8) :: z1, z2, rlam, rlam1, rlam2, vec(3), pos(3)
!
!  Start of Executable Code:
!  ------------------------
!
      if( geom%lcart ) then
         z1 = pos1(3)
         z2 = pos2(3)
      else
         call toCurviCoord1d( pos1(:), z1 )
!             ==============     i     o
         call toCurviCoord1d( pos2(:), z2 )
!             ==============     i     o
      end if
      call applyBndCond1d( z1 )
!          ==============
      call applyBndCond1d( z2 )
!          ==============
!
!  Determine surface
!
      if( z1 >= zrefract(1) .and. z2 <= zrefract(1) ) then
         isurf = 1                              ! top/inner surface
      else if( z1 <= zrefract(2) .and. z2 >= zrefract(2) ) then
         isurf = 2                              ! bottom/outer surface
      else
         found = .false.
         return                                 ! erroneous pos1, pos2
      end if
!
!  Get exit point bracket
!
      vec(:) = pos2(:)-pos1(:)
      call intersect( pos1(:), vec(:), zrefract(isurf), .true., rlam, found )
!          =========    i        i         i           i       o     o
      if( .not. found )  return
      pos(:) = pos1(:) + rlam * vec(:)

      call intersect( pos(:), -vec(:), zrefract1(isurf), .true., rlam1, found )
!          =========    i        i          i            i      o      o
      if( .not. found )  return
      call intersect( pos(:), vec(:), zrefract2(isurf), .true., rlam2, found )
!          =========    i       i          i            i      o      o
      if( .not. found )  return
      pos1(:) = pos(:) - rlam1 * vec(:)
      pos2(:) = pos(:) + rlam2 * vec(:)

   end subroutine getExitPoint1d
!
!=======================================================================
!
   logical function isInsideGrid( pos )
!                   ============   i
!
!     Determine if position pos is inside the grid.
!
!  Arguments:
!
      real(REAL8), intent(in) :: pos(3)        ! position on the line
!
!  Local Variables:
!
      real(REAL8) :: z
!
!  Start of Executable Code:
!  ------------------------
!
      if( geom%lcart ) then
         z = pos(3)
      else
         call toCurviCoord1d( pos(:), z )
!             ==============    i     o
      end if
      call applyBndCond1d( z )
!          ==============

      isInsideGrid = (z > geom%zlim(1)) .and. (z < geom%zlim(2))

   end function isInsideGrid
!
!=======================================================================
!
   subroutine intersect( pos, vec, z, onlyforward, rlam, found )
!             =========   i    i   i       i        o      o
!
!     Determine rlam such that pos(:)+rlam*vec(:) is on the plane/cylinder/sphere defined
!     by z.
!
!     For lcylx=t, lcyly=t, lcylz=t, or lspher=t the smaller of the two rlam is returned
!     unless onlyforward=t and the smaller rlam is negative.
!
!  Arguments:
!
      real(REAL8), intent(in) :: pos(3)         ! position on the line
      real(REAL8), intent(in) :: vec(3)         ! vector definining the direction of the line
      real(REAL8), intent(in) :: z              ! position of the plane or radius of the
                                                !    cylinder/sphere
      logical, intent(in) :: onlyforward        ! search only in forward direction
      real(REAL8), intent(out) :: rlam          ! the parameter to be determined
      logical, intent(out) :: found             ! flag indicating success
!
!  Local Variables:
!
      real(REAL8) :: a, b, d, vec2
!
!  Start of Executable Code:
!  ------------------------
!
!  Determine whether there is an intersection point of the unlimited line
!
      if( geom%lcart ) then

         found = ( abs(vec(3)) >= SMALL )

      else

         if( z <= ZERO ) then
            found = .false.
         else
            if( geom%lcylx ) then
               a = ( geom%ycen - pos(2) ) * vec(2) + ( geom%zcen - pos(3) ) * vec(3)
               b = ( geom%ycen - pos(2) )**2 + ( geom%zcen - pos(3) )**2 - z**2
               vec2 = vec(2)**2 + vec(3)**2
            else if( geom%lcyly ) then
               a = ( geom%xcen - pos(1) ) * vec(1) + ( geom%zcen - pos(3) ) * vec(3)
               b = ( geom%xcen - pos(1) )**2 + ( geom%zcen - pos(3) )**2 - z**2
               vec2 = vec(1)**2 + vec(3)**2
            else if( geom%lcylz ) then
               a = ( geom%xcen - pos(1) ) * vec(1) + ( geom%ycen - pos(2) ) * vec(2)
               b = ( geom%xcen - pos(1) )**2 + ( geom%ycen - pos(2) )**2 - z**2
               vec2 = vec(1)**2 + vec(2)**2
            else if( geom%lspher ) then
               a = ( geom%xcen - pos(1) ) * vec(1) + ( geom%ycen - pos(2) ) * vec(2) +  &
                   ( geom%zcen - pos(3) ) * vec(3)
               b = ( geom%xcen - pos(1) )**2 + ( geom%ycen - pos(2) )**2 +  &
                   ( geom%zcen - pos(3) )**2 - z**2
               vec2 = vec(1)**2 + vec(2)**2 + vec(3)**2
            else
               error stop 'intersect: Inconsistent geom%lcart'
            end if
            d = a**2 - b*vec2
            found = ( d >= ZERO )
         end if

      end if

      if( .not. found )  return
!
!  Determine intersection point
!
      if( geom%lcart ) then

         rlam = ( z - pos(3) ) / vec(3)
         if( onlyforward )  found = ( rlam > ZERO )

      else

         if( a > ZERO ) then                    ! avoid numerical error
            rlam = b / ( a + sqrt(d) )
         else
            rlam = ( a - sqrt(d) ) / vec2
         end if
         if( onlyforward ) then
            found = ( rlam > ZERO )
            if( .not.found ) then
               if( a > ZERO ) then
                  rlam = ( a + sqrt(d) ) / vec2
               else
                  rlam = b / ( a - sqrt(d) )
               end if
               found = ( rlam > ZERO )
            end if
         end if

      end if

   end subroutine intersect
!
!=======================================================================
!
end module Target_Surface1d_mod
