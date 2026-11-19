module Target_Surface_mod
!
!  Provide procedures for ONLY geometry related operations on the target
!
!  Used Modules:
!  ============
!
   use Const_mod,             only: INT4, REAL8,  &
                                    ONE, HALF, EPSMB, SMALL

   use Types_mod,             only: Proj_type

   use Param_Setup_mod,       only: ndim
   use Param_Sn_mod,          only: sn

   use Target_Surface1d_mod,  only: initSurface1d, getEntranceGrid1d, getDistSurface1d,  &
                                    getExitPoint1d
   use Target_Surface2d_mod,  only: initSurface2d, getEntranceGrid2d, getDistSurface2d
   use Target_Surface3d_mod,  only: initSurface3d, getEntranceGrid3d, getDistSurface3d
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
   public :: initSurface               ! initialize the surface
!            ===========
   public :: getDistSurface            ! get the signed distance from the surface
!            ==============            !    (-) outside, (+) inside
   public :: getEntranceGrid           ! get the entrance point into the grid covering the
!            ===============           !    target
   public :: getEntrancePoint          ! get the entrance point into the target (to be more
!            ================          !    precise, a distance pmax from the target)
   public :: getExitPoint              ! get the exit point from the target (to be more
!            ============              !    precise, a distance pmax from the target)
!
contains
!
!=======================================================================
!
   subroutine initSurface( )
!             ===========
!
!     Initialize the surface
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
      if( ndim == 1 ) then
         call initSurface1d( )
!             =============
      else if( ndim == 2 ) then
         call initSurface2d( )
!             =============
      else
         call initSurface3d( )
!             =============
      end if

      end subroutine initSurface
!
!=======================================================================
!
   subroutine getDistSurface( pos, dist, outside_, dir_ )
!             ==============   i    o       o       o
!  
!     Determine the distance of the point from the target surface
!
!  Arguments:
!
      real(REAL8), intent(in) :: pos(3)            ! position to be checked
      real(REAL8), intent(out) :: dist             ! signed distance: (+) inside, (-) outside
                                                   !    set to -BIG if position is outside grid
                                                   !    in the 2d case
      logical, optional, intent(out) :: outside_   ! indicates point is outside grid
      real(REAL8), optional, intent(out) :: dir_(3)! outwards direction
!
!  Local Variables:
!
!     None
!
!  Start of Executable Code:
!  ------------------------
!
      if( ndim == 1 ) then
         call getDistSurface1d( pos(:), dist, outside_, dir_ )
!             ================    i      o       o       o
      else if( ndim == 2 ) then
         call getDistSurface2d( pos(:), dist, outside_, dir_ )
!             ================    i      o       o       o
      else
         call getDistSurface3d( pos(:), dist, outside_, dir_ )
!             ================    i      o       o       o
      end if

   end subroutine getDistSurface
!
!=======================================================================
!
   subroutine getEntranceGrid( proj, found, onlyforward_ )
!             ===============   io     o         i
!  
!     Determine the entrance point of a projectile into the grid covering the target
!     The entrance point may also lie backwards of the projectile if onlyforward=f.
!     The entrance point is identical to the current position if the position is inside the
!     grid and onlyforward=t.
!
!  Arguments:
!
      type(Proj_type), intent(inout) :: proj       ! projectile
      logical, intent(out) :: found                ! flag indicating whether proj%pos contains
                                                   !    a valid position
      logical, optional, intent(in) :: onlyforward_! flag indicating that entrance point may
                                                   !    only lie in the forward direction
                                                   !    (default: false)
!
!  Local Variables:
!
      logical :: onlyforward
!
!  Start of Executable Code:
!  ------------------------
!
      if( present( onlyforward_ ) ) then
         onlyforward = onlyforward_
      else
         onlyforward = .false.
      end if

      if( ndim == 1 ) then
         call getEntranceGrid1d( proj, found, onlyforward )
!             =================   io     o         i
      else if( ndim == 2 ) then
         call getEntranceGrid2d( proj, found, onlyforward )
!             =================   io     o         i
      else
         call getEntranceGrid3d( proj, found, onlyforward )
!             =================   io     o         i
      end if

   end subroutine getEntranceGrid
!
!=======================================================================
!
   subroutine getEntrancePoint( pos1, pos2, found )
!             ================   io    io     o
!  
!     Determine the position on the line pos1-pos2 where the distance from the target
!     is sn%pmax. pos1 is assumed to be further away from the target than sn%pmax,
!     pos2 is closer to the target than sn%pmax
!
!     On exit pos1 and pos2 contain approximations to a position sn%pmax away from the 
!     target (they bracket the exact solution)
!
!  Arguments:
!
      real(REAL8), intent(inout):: pos1(3)      ! position further away than sn%pmax
      real(REAL8), intent(inout):: pos2(3)      ! position closer than sn%pmax
      logical, intent(out) :: found             ! flag indicating success
!
!  Local Variables:
!
!     None
!
!  Start of Executable Code:
!  ------------------------
!
      call getExitPoint( pos2(:), pos1(:), found )
!          ============    io       io       o
!
   end subroutine getEntrancePoint
!
!=======================================================================
!
   subroutine getExitPoint( pos1, pos2, found )
!             ============   io    io     o
!  
!     Determine the position on the line pos1-pos2 where the distance from the target
!     is sn%pmax. pos1 is assumed to be closer to the target than sn%pmax, pos2 is
!     further away from the target than sn%pmax
!
!     On exit pos1 and pos2 contain approximations to a position sn%pmax away from the 
!     target, pos1 inside and pos2 outside -sn%pmax (they bracket the exact solution)
!
!  Arguments:
!
      real(REAL8), intent(inout):: pos1(3)      ! position closer than sn%pmax
      real(REAL8), intent(inout):: pos2(3)      ! position further away than sn%pmax
      logical, intent(out) :: found             ! flag indicating success
!
!  Local Variables:
!
      real(REAL8) :: dist, dist1, dist2, pos(3), pos1_tmp(3)
!
!  Start of Executable Code:
!  ------------------------
!
!  More efficient treatment of 1D case
!
      if( ndim == 1 ) then
         call getExitPoint1d( pos1(:), pos2(:), found )
!             ==============    io       io       o
         return
      end if
!
!  Initialize iteration
!
      call getDistSurface( pos1(:), dist1 )
!          ==============     i       o
      call getDistSurface( pos2(:), dist2 )
!          ==============     i       o
!
      if( abs( dist2 - dist1 ) < SMALL .or. dist2 > -sn%pmax .or. dist1 < -sn%pmax ) then
         found = .false.
         return
      end if

      found = .true.
!
!  Find two positions near -sn%pmax by bisection, exit before |dist2-dist1| < EPSMB*sn%pmax
!  so we can still calculate the derivative
!
      do
         pos(:) = HALF * ( pos1(:) + pos2(:) )

         call getDistSurface( pos(:), dist )
!             ==============    i      o

         if( dist > -sn%pmax ) then             ! pos will replace pos1
            if( abs( dist - dist2 ) < EPSMB * sn%pmax )  exit
            pos1(:) = pos(:)
            dist1 = dist
         else                                   ! pos will replace pos2
            if( abs( dist - dist1 ) < EPSMB * sn%pmax )  exit
            pos2(:) = pos(:)
            dist2 = dist
         end if
      end do
!
!  Along the line connecting the two positions find the points with dist = -sn%pmax*(1+/-EPSMB)
!  assuming linear dependence of dist on pos
!
      pos1_tmp(:) = pos1(:)
      pos1(:) = pos1_tmp(:) + ( -sn%pmax*(ONE-EPSMB) - dist1 ) / (dist2 - dist1 ) *  &
                              ( pos2(:) - pos1_tmp(:) )
      pos2(:) = pos1_tmp(:) + ( -sn%pmax*(ONE+EPSMB) - dist1 ) / (dist2 - dist1 ) *  &
                              ( pos2(:) - pos1_tmp(:) )

   end subroutine getExitPoint
!
!=======================================================================
!
end module Target_Surface_mod
