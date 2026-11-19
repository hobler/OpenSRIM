module Geom_CoordSys_mod
!
!  Manages the transformation between the target coordinate system (K) and a second
!  coordinate system (K'). Both coordinate systems are orthogonal. The target coordinate
!  system is Cartesian and is considered the primary coordinate system. The second coordinate
!  system may have scaled axes, which allows to represents a crystallographic lattice by
!  fractional coordinates.
!
!     K  ... target coordinate system (units: A)
!     K' ... second/lattice coordinate system (units: A or stretch(1:3))
!
!  The transformation of the coordinates (x,y,z) of a point P in K to its coordinates (x',y',z')
!  in K' is given by
!
!              |x'|           |x|            | 1/stretch(1)     0         0     |
!     K -> K': |y'| = A * T * |y|   with A = |     0     1/stretch(2)     0     |
!              |z'|           |z|            |     0         0     1/stretch(3) |
!
!                                           |ex(1:3)|
!                                   and T = |ey(1:3)| = |ex'(1:3) ey'(1:3) ez'(1:3)|
!                                           |ez(1:3)|
!
!  where ex, ey, and ez are the unit vectors of the K' system in coordinates of K,
!  thus T defines a rotation matrix. ex', ey', and ez' are the unit vectors of the K system
!  in coordinates of K'. The former are stored as properties of the coordinate system
!
!  "to" coordinate system means to calculate (x',y',z') from (x,y,z)
!  "from" coordinate system means to calculate (x,y,z) from (x',y',z')
!
!  Note: If we ever need rotation about an arbitrary axis, google for "Rodriguez formula".
!
!  Used Modules:
!  ============
!
   use Const_mod,    only: REAL8, INT4,  &
                           ZERO, ONE, EPSMA, PI
   use Vector_mod,   only: operator(.dot.), operator(.cross.), abs, abs2
   use Random_mod,   only: ranrot
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
   public :: initCoordSys                    ! define transformation from K to K'
!            ============
   public :: toCoordSys                      ! determine coordinates in K' from
!            ==========                      !    coordinates in K
   public :: fromCoordSys                    ! determine coordinates in K from
!            ============                    !    coordinates in K'
   public :: getVolume                       ! calculate the volume of the unit cell
!            =========
   public :: ranrotCoordSys                  ! generate randomly oriented coordinate system K'
!            ==============
   public :: invertCoordSys                  ! exchange the roles of K and K'
!            ==============
!
!  Private Procedures:
!  ==================
!
   private :: rotXCoordSys                   ! rotate K' about x' axis
!             ============
   private :: rotYCoordSys                   ! rotate K' about y' axis
!             ============
   private :: rotZCoordSys                   ! rotate K' about z' axis
!             ============
!
!  Public Data Types:
!  =================
!
   type, public :: CoordSys_type             ! data type of coordinate system (K')
!                  =============
      !private
      logical :: rotated                     ! flag indicating T is not the identity matrix
      logical :: scaled                      ! flag indicating stretch(:) not equal [ONE,ONE,ONE]
      real(REAL8) :: ex(3)                   ! x', y', and z' axis unit vectors in coordinates
      real(REAL8) :: ey(3)                   !    of K (corresponding to ex, ey, ez in the
      real(REAL8) :: ez(3)                   !    comment above)
      real(REAL8) :: stretch(3)              ! scaling factors = lengths of axes of K'
                                             !    measured in units of K
   end type CoordSys_type

contains
!
!=======================================================================
!
   subroutine initCoordSys( coord, dirz_, dirx_, stretch_, rotate_, tiltx_, tilty_ )
!             ============    o     i     i       i       i       i      i
!
!     Define transformation from K to K'
!
!     If only coord is present, K = K'.
!     If neither dirz_ nor dirx_ is present, the axes of K and K' are parallel.
!     If dirz_ is present but not dirx_, the y axes of K and K' are assumed to be identical, i.e.,
!     K is rotated about its y axis to arrive at K'.
!     If stretch_ is not present, stretch_=[1,1,1], i.e., the axes are unscaled.
!     If rotate_, tiltx_, and/or tilty_ are not specified, they are assumed zero

!     If rotate_, tiltx_, and/or tilty_ are specified, they are applied in the order
!
!        rotate_ (about z axis) -> tiltx_ -> tilty_
!
!  Arguments:
!
      type(CoordSys_type), intent(out) :: coord       ! coordinate system K'
      real(REAL8), optional, intent(in) :: dirz_(3)   ! direction of z' axis in coordinates of K
      real(REAL8), optional, intent(in) :: dirx_(3)   ! direction of x' axis in coordinates of K
                                                      !    Note: dirz_ and dirx_ need not be
                                                      !    scaled to 1
      real(REAL8), optional, intent(in) :: stretch_(3)! unit lengths for the x', y', and z' axis
      real(REAL8), optional, intent(in) :: rotate_    ! rotation angle about z axis (deg)
      real(REAL8), optional, intent(in) :: tiltx_     ! tilt angle about x axis (deg)
      real(REAL8), optional, intent(in) :: tilty_     ! tilt angle about y axis (deg)
!
!  Local Variables:
!
!     None
!
!  Start of Executable Code:
!  ------------------------
!
!  Define unit vectors in z' direction
!
      if( present( dirz_ ) ) then
         coord%ez(:) = dirz_(:) / abs( dirz_(:) )
!                                 ===
      else
         coord%ez(:) = [ ZERO, ZERO, ONE ]
      end if
!
!  Define unit vectors in x' and y' direction
!
      if( present( dirx_ ) ) then

         coord%ex(:) = dirx_(:) / abs( dirx_(:) )
!                                 ===
         if( abs(coord%ex(:).dot.coord%ez(:)) > 10*EPSMA )  &
!            ===
               error stop 'initCoordSys: x and z axes not perpendicular'

         coord%ey(:) = coord%ez(:) .cross. coord%ex(:)

      else

         coord%ey(:) = [ ZERO, ZERO, ONE ] .cross. coord%ez(:)

         if( abs2( coord%ey(:) ) > EPSMA ) then
!            ====
            coord%ey(:) = coord%ey(:) / abs( coord%ey(:) )
!                                       ===
         else
            coord%ey(:) = [ ZERO, ONE, ZERO ]
         end if

         coord%ex(:) = coord%ey(:) .cross. coord%ez(:)

      end if
!
!  Rotate and tilt
!
      if( present( rotate_ ) )  call rotZCoordSys( coord, rotate_ )
!                                    ============   io       i
      if( present( tiltx_ ) )  call rotXCoordSys( coord, tiltx_ )
!                                   ============   io      i
      if( present( tilty_ ) )  call rotYCoordSys( coord, tilty_ )
!                                   ============   io      i
!
!  Is there a rotation at all?
!
      coord%rotated = abs2( coord%ex(:) - [ONE,ZERO,ZERO] ) > EPSMA  .or.  &
!                     ====
                      abs2( coord%ez(:) - [ZERO,ZERO,ONE] ) > EPSMA
!                     ====
!
!  Scaling length
!
      if( present( stretch_ ) ) then
         coord%stretch(:) = stretch_(:)
      else
         coord%stretch(:) = ONE
      end if
!
!  Is there a scaling at all?
!
      coord%scaled = abs2( coord%stretch(:) - [ONE,ONE,ONE] ) > EPSMA
!                    ====
   end subroutine initCoordSys
!
!=======================================================================
!
   subroutine toCoordSys( coord, v1, v2 )
!             ==========    i    i   o
!
!     Determine coordinates v2 in K' from coordinates v1 in K.
!
!  Arguments:
!
      type(CoordSys_type), intent(in) :: coord  ! coordinate system
      real(REAL8), intent(in) :: v1(3)          ! coordinates in K
      real(REAL8), intent(out) :: v2(3)         ! coordinates in K'
!
!  Local Variables:
!
!     None
!
!  Start of Executable Code:
!  ------------------------
!
!  Rotate
!
      if( coord%rotated ) then
         v2(1) = coord%ex(:) .dot. v1(:)
         v2(2) = coord%ey(:) .dot. v1(:)
         v2(3) = coord%ez(:) .dot. v1(:)
      else
         v2(:) = v1(:)
      end if
!
!  Scale
!
      if( coord%scaled )  v2(:) = v2(:) / coord%stretch(:)

   end subroutine toCoordSys
!
!=======================================================================
!
   subroutine fromCoordSys( coord, v1, v2 )
!             ============    i    i   o
!
!     Determine coordinates v2 in K from coordinates v1 in K".
!
!  Arguments:
!
      type(CoordSys_type), intent(in) :: coord  ! coordinate system
      real(REAL8), intent(in) :: v1(3)          ! coordinates in K
      real(REAL8), intent(out) :: v2(3)         ! coordinates in K'(K")
!
!  Local Variables:
!
      real(REAL8) :: ex(3), ey(3), ez(3)        ! unit vectors of K in K'
      real(REAL8) :: v1_unscal(3)
!
!  Start of Executable Code:
!  ------------------------
!
!  Unscale
!
      if( coord%scaled ) then
         v1_unscal(:) = coord%stretch(:) * v1(:)
      else
         v1_unscal(:) = v1(:)
      end if
!
!  Rotate
!
      if( coord%rotated ) then
         ex(:) = [ coord%ex(1) , coord%ey(1) , coord%ez(1) ]
         ey(:) = [ coord%ex(2) , coord%ey(2) , coord%ez(2) ]
         ez(:) = [ coord%ex(3) , coord%ey(3) , coord%ez(3) ]

         v2(1) = ex(:) .dot. v1_unscal(:)
         v2(2) = ey(:) .dot. v1_unscal(:)
         v2(3) = ez(:) .dot. v1_unscal(:)
      else
         v2(:) = v1_unscal(:)
      end if

   end subroutine fromCoordSys
!
!=======================================================================
!
   subroutine getVolume( coord, vol )
!             =========    i     o
!
!     calculate the volume of the unit cell (assumes orthogonal axes!)
!
!  Arguments:
!
      type(CoordSys_type), intent(in) :: coord  ! coordinate system
      real(REAL8) :: vol                        ! volume of unit cell
!
!  Local Variables:
!
!     None
!
!  Start of Executable Code:
!  ------------------------
!
      vol = product( coord%stretch(:) )

   end subroutine getVolume
!
!=======================================================================
!
   subroutine ranrotCoordSys( coord )
!             ==============   io
!
!     Generate a randomly oriented coordinate system K'
!
!  Arguments:
!
      type(CoordSys_type), intent(out) :: coord ! coordinate system K'
!
!  Local Variables:
!
      real(REAL8) :: dirx(3), diry(3), dirz(3)
!
!  Start of Executable Code:
!  ------------------------
!
      call ranrot( dirx, diry, dirz )
!          ======   o     o     o

      coord%ex(:) = dirx(:)
      coord%ey(:) = diry(:)
      coord%ez(:) = dirz(:)

      coord%rotated = .true.

   end subroutine ranrotCoordSys
!
!=======================================================================
!
   subroutine invertCoordSys( coord )
!             ==============   io
!
!     exchange the role of K and K'
!
!  Arguments:
!
      type(CoordSys_type), intent(inout) :: coord  ! K' defined in K
!
!  Local Variables:
!
      real(REAL8) :: ex(3), ey(3), ez(3)        ! unit vectors of new K'
!
!  Start of Executable Code:
!  ------------------------
!
      if( coord%rotated ) then
         ex(:) = [ coord%ex(1) , coord%ey(1) , coord%ez(1) ]
         ey(:) = [ coord%ex(2) , coord%ey(2) , coord%ez(2) ]
         ez(:) = [ coord%ex(3) , coord%ey(3) , coord%ez(3) ]

         coord%ex(:) = ex(:)
         coord%ey(:) = ey(:)
         coord%ez(:) = ez(:)
      end if

      if( coord%scaled )  coord%stretch(:) = ONE / coord%stretch(:)

   end subroutine invertCoordSys
!
!=======================================================================
!
   subroutine rotXCoordSys( coord, angle )
!             ============   io      i
!
!     Rotate the coordinate system K' about its x' axis
!
!  Arguments:
!
      type(CoordSys_type), intent(inout) :: coord  ! coordinate system K'
      real(REAL8), intent(in) :: angle             ! rotation angle (deg)
!
!  Local Variables:
!
      real(REAL8) :: diry(3), dirz(3)
      real(REAL8) :: sina, cosa
!
!  Start of Executable Code:
!  ------------------------
!
      cosa = cos(angle*PI/180)
      sina = sin(angle*PI/180)

      diry(:) = cosa * coord%ey(:) + sina * coord%ez(:)
      dirz(:) = - sina * coord%ey(:) + cosa * coord%ez(:)

      coord%ey(:) = diry(:)
      coord%ez(:) = dirz(:)

      coord%rotated = abs2( coord%ex(:) - [ONE,ZERO,ZERO] ) > EPSMA  .or.  &
!                     ====
                      abs2( coord%ez(:) - [ZERO,ZERO,ONE] ) > EPSMA
!                     ====

   end subroutine rotXCoordSys
!
!=======================================================================
!
   subroutine rotYCoordSys( coord, angle )
!             ============   io      i
!
!     Rotate the coordinate system K' about its y' axis
!
!  Arguments:
!
      type(CoordSys_type), intent(inout) :: coord  ! coordinate system K'
      real(REAL8), intent(in) :: angle             ! rotation angle (deg)
!
!  Local Variables:
!
      real(REAL8) :: dirz(3), dirx(3)
      real(REAL8) :: sina, cosa
!
!  Start of Executable Code:
!  ------------------------
!
      cosa = cos(angle*PI/180)
      sina = sin(angle*PI/180)

      dirx(:) = cosa * coord%ex(:) - sina * coord%ez(:)
      dirz(:) = sina * coord%ex(:) + cosa * coord%ez(:)

      coord%ex(:) = dirx(:)
      coord%ez(:) = dirz(:)

      coord%rotated = abs2( coord%ex(:) - [ONE,ZERO,ZERO] ) > EPSMA  .or.  &
!                     ====
                      abs2( coord%ez(:) - [ZERO,ZERO,ONE] ) > EPSMA
!                     ====

   end subroutine rotYCoordSys
!
!=======================================================================
!
   subroutine rotZCoordSys( coord, angle )
!             ============   io      i
!
!     Rotate the coordinate system K' about its z' axis
!
!  Arguments:
!
      type(CoordSys_type), intent(inout) :: coord  ! coordinate system K'
      real(REAL8), intent(in) :: angle             ! rotation angle (deg)
!
!  Local Variables:
!
      real(REAL8) :: dirx(3), diry(3)
      real(REAL8) :: sina, cosa
!
!  Start of Executable Code:
!  ------------------------
!
      cosa = cos(angle*PI/180)
      sina = sin(angle*PI/180)

      dirx(:) = cosa * coord%ex(:) + sina * coord%ey(:)
      diry(:) = - sina * coord%ex(:) + cosa * coord%ey(:)

      coord%ex(:) = dirx(:)
      coord%ey(:) = diry(:)

      coord%rotated = abs2( coord%ex(:) - [ONE,ZERO,ZERO] ) > EPSMA  .or.  &
!                     ====
                      abs2( coord%ez(:) - [ZERO,ZERO,ONE] ) > EPSMA
!                     ====

   end subroutine rotZCoordSys
!
!=======================================================================
!
end module Geom_CoordSys_mod
