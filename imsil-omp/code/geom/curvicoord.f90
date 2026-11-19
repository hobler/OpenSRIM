module Geom_CurviCoord_mod
!
!  Manages the transformation between the Cartesian target coordinate system and the
!  the curvilinear (cylindrical or spherical) coordinate system.
!
!  Used Modules:
!  ============
!
   use Const_mod,       only: REAL8, INT4,  &
                              ZERO, ONE, TWO,  &
                              ENUM
   use Undefine_mod,    only: undefine

   use Param_Geom_mod,  only: geom
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
   public :: toCurviCoord1d                  ! given a position in the target coordinate
!            =============                   !    system, get the 1d cylindrical or spherical
                                             !    coordinate
   public :: toCurviCoord2d                  ! given a position in the target coordinate
!            =============                   !    system, get the 2d cylindrical coordinates
   public :: fromCurviZDir                   ! get the direction (defined in the target
!            ============                    !    coordinate system) corresponding to the
                                             !    positive z direction of the curvilinear
                                             !    coordinate system
   public :: fromGeomNormDir2d               ! convert surface normal direction from curvilinear
!            =================               !    coordinate system to the target coordinate
!                                            !    system
!
contains
!
!=======================================================================
!
   subroutine toCurviCoord1d( pos, z )
!             ==============   i   o
!
!     Get the 1d cylindrical or spherical coordinate z ("radius") for a given position pos.
!
!  Arguments:
!
      real(REAL8), intent(in) :: pos(3)      ! position given in target coordinates
      real(REAL8), intent(out) :: z          ! For cylinder coordinate system:
                                             !    distance from axis
                                             ! For spherical coordinate system:
                                             !    distance from center
!
!  Local Variables:
!
!     None
!
!  Start of Executable Code:
!  ------------------------
!
      if( geom%lcylx ) then
         z = sqrt( (pos(2)-geom%ycen)**2 + (pos(3)-geom%zcen)**2 )
      else if( geom%lcyly ) then
         z = sqrt( (pos(1)-geom%xcen)**2 + (pos(3)-geom%zcen)**2 )
      else if( geom%lcylz ) then
         z = sqrt( (pos(1)-geom%xcen)**2 + (pos(2)-geom%ycen)**2 )
      else if( geom%lspher ) then
         z = sqrt( (pos(1)-geom%xcen)**2 + (pos(2)-geom%ycen)**2 + (pos(3)-geom%zcen)**2 )
      else
         error stop 'toCurviCoord1d: must be called for cylindrical or spherical coordinates'
      end if

   end subroutine toCurviCoord1d
!
!=======================================================================
!
   subroutine toCurviCoord2d( pos, x, z )
!             ==============   i   o  o
!
!     Get the cylindrical coordinates x (along the axis) and z ("radius") for a given position
!     pos.
!
!  Arguments:
!
      real(REAL8), intent(in) :: pos(3)      ! position given in target coordinates
      real(REAL8), intent(out) :: x          ! Coordinate along axis measured from "center"
      real(REAL8), intent(out) :: z          ! Distance from axis (radius)
!
!  Local Variables:
!
!     None
!
!  Start of Executable Code:
!  ------------------------
!
      if( geom%lcylx ) then
         x = pos(1) - geom%xcen
         z = sqrt( (pos(2)-geom%ycen)**2 + (pos(3)-geom%zcen)**2 )
      else if( geom%lcyly ) then
         x = pos(2) - geom%ycen
         z = sqrt( (pos(1)-geom%xcen)**2 + (pos(3)-geom%zcen)**2 )
      else if( geom%lcylz ) then
         x = pos(3) - geom%zcen
         z = sqrt( (pos(1)-geom%xcen)**2 + (pos(2)-geom%ycen)**2 )
      else
         error stop 'toCurviCoord2d: must be called for cylindrical coordinates only'
      end if

   end subroutine toCurviCoord2d
!
!=======================================================================
!
   subroutine fromCurviZDir( pos, dir )
!             =============   i    o
!
!     Get the direction in the target coordinate system corresponding to the
!     positive z direction of the geometry coordinate system
!
!  Arguments:
!
      real(REAL8), intent(in) :: pos(3)      ! position given in target coordinates
      real(REAL8), intent(out) :: dir(3)     ! direction of the z coordinate of the geometry
                                             !    coordinate system measured in the target
                                             !    coordinate system
!
!  Local Variables:
!
      real(REAL8) :: z
!
!  Start of Executable Code:
!  ------------------------
!
      if( geom%lcylx ) then
         z = sqrt( (pos(2)-geom%ycen)**2 + (pos(3)-geom%zcen)**2 )
         dir(1) = ZERO
         dir(2) = (pos(2)-geom%ycen) / z
         dir(3) = (pos(3)-geom%zcen) / z
      else if( geom%lcyly ) then
         z = sqrt( (pos(1)-geom%xcen)**2 + (pos(3)-geom%zcen)**2 )
         dir(1) = (pos(1)-geom%xcen) / z
         dir(2) = ZERO
         dir(3) = (pos(3)-geom%zcen) / z
      else if( geom%lcylz ) then
         z = sqrt( (pos(1)-geom%xcen)**2 + (pos(2)-geom%ycen)**2 )
         dir(1) = (pos(1)-geom%xcen) / z
         dir(2) = (pos(2)-geom%ycen) / z
         dir(3) = ZERO
      else if( geom%lspher ) then
         z = sqrt( (pos(1)-geom%xcen)**2 + (pos(2)-geom%ycen)**2 + (pos(3)-geom%zcen)**2 )
         dir(1) = (pos(1)-geom%xcen) / z
         dir(2) = (pos(2)-geom%ycen) / z
         dir(3) = (pos(3)-geom%zcen) / z
      else
         error stop 'fromCurviZDir: must be called for cylindrical or spherical coordinates'
      end if

   end subroutine fromCurviZDir
!
!=======================================================================
!
   subroutine fromGeomNormDir2d( zdir, dirx, dirz, dir )
!             =================   i     i     i     o
!
!     Convert the surface normal direction (dirx, dirz) in the 2d geometry coordinate system
!     to the surface normal direction in the target coordinate system.
!
!  Arguments:
!
      real(REAL8), intent(in) :: zdir(3)     ! direction of the geometry z coordinate
      real(REAL8), intent(in) :: dirx        ! x component of the surface normal direction in
                                             !    the geometry coordinate system
      real(REAL8), intent(in) :: dirz        ! z component of the surface normal direction in
                                             !    the geometry coordinate system
      real(REAL8), intent(out) :: dir(3)     ! surface normal direction in the target
                                             !    coordinate system
!
!  Local Variables:
!
!     None
!
!  Start of Executable Code:
!  ------------------------
!
      dir(:) = zdir(:) * dirz
      if( geom%lcylx ) then
         dir(1) = dirx
      else if( geom%lcyly ) then
         dir(2) = dirx
      else if( geom%lcylz ) then
         dir(3) = dirx
      else
         error stop 'fromGeomNormDir2d: must be called for cylindrical coordinates only'
      end if

   end subroutine fromGeomNormDir2d
!
!=======================================================================
!
end module Geom_CurviCoord_mod
