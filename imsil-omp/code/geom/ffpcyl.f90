module Geom_FFPcyl_mod
!
!  Provide geometric operations on the projection of free flight path cylinder to the 
!  xz-plane, which consist of two semi-ellipses and two connecting straight lines.
!
!  Used Modules:
!  ============
!
   use Const_mod,       only: INT4, REAL8,  &
                              ZERO

   use Geom_Line_mod,   only: mcisct
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
   public :: isInsideFFPcyl                  ! determine whether point is inside ffp 
!            ==============                  !    cylinder
   public :: intersectFFPcyl                 ! determine intersection of straight line
!            ===============                 !    with ffp cylinder

contains
!
!=======================================================================
!
   logical function isInsideFFPcyl( pos, pos0, dir, ffp, pmax, b_a )  result( isInside )
!                   ==============   i    i     i    i     i    i                o
!
!     Determine whether point pos=(x,z) is inside the free flight path cylinder
!
!  Arguments:
!
      real(REAL8), intent(in) :: pos(2)         ! point to be investigated
      real(REAL8), intent(in) :: pos0(2)        ! starting point of ffp
      real(REAL8), intent(in) :: dir(2)         ! direction of ffp (dir(1)**2+dir(2)**2=1)
      real(REAL8), intent(in) :: ffp            ! projection of ffp to xz-plane = of height of
                                                !    ffp cylinder
      real(REAL8), intent(in) :: pmax           ! maximum impact parameter = radius of ffp 
                                                !    cylinder
      real(REAL8), intent(in) :: b_a            ! ratio of axis lengths of ellipe (<=1)
!
!  Local Variables:
!
      real(REAL8) :: xi, zeta
!
!  Start of Executable Code:
!  ------------------------
!
!  Coordinates in the (ffp,pmax)-coordinate system
!
      xi   = dir(1) * ( pos(2) - pos0(2) ) - dir(2) * ( pos(1) - pos0(1) )
      zeta = dir(1) * ( pos(1) - pos0(1) ) + dir(2) * ( pos(2) - pos0(2) )
!
!  Check position
!
      if( abs(xi) > pmax ) then
         isInside =.false.                            ! outside infinite cylinder
      else if( zeta >= ZERO .and. zeta <= ffp ) then
         isInside =.true.                             ! inside rectangle
      else if( b_a**2 * xi**2 + zeta**2 <= b_a**2 * pmax**2 ) then
         isInside =.true.                             ! inside semiellipse around ion
      else if( b_a**2 * xi**2 + ( zeta - ffp )**2 <= b_a**2 * pmax**2 ) then
         isInside =.true.                             ! inside other semiellipse
      else
         isInside =.false.
      end if

   end function isInsideFFPcyl
!
!=======================================================================
!
   logical function intersectFFPcyl( pos1, pos2, pos0, dir, ffp, pmax, b_a ) result( cut )
!                   ===============   i     i     i     i    i    i     i             o
!
!     Determine whether a straight line pos1-pos2 intersects the boundary of the free
!     flight path cylinder. If there is an intersection point return its position along
!     pos1-pos2.
!
!     NOTE: an intersection point is also reported if the line intersects the semi-ellipses
!           inside the ffp cylinder. This could be avoided by considering that no intersection
!           point may exist if both points are insided the ffp cylinder (using isInsideFFPcyl)
!
!  Arguments:
!
      real(REAL8), intent(in) :: pos1(2)        ! starting point of straight line
      real(REAL8), intent(in) :: pos2(2)        ! end point of straight line
      real(REAL8), intent(in) :: pos0(2)        ! starting point of ffp
      real(REAL8), intent(in) :: dir(2)         ! direction of ffp (dir_p(1)**2+dir_p(2)**2=1)
      real(REAL8), intent(in) :: ffp            ! projection of ffp to xz-plane = of height of
                                                !    ffp cylinder
      real(REAL8), intent(in) :: pmax           ! maximum impact parameter = radius of ffp 
                                                !    cylinder
      real(REAL8), intent(in) :: b_a            ! ratio of axis lengths of ellipe (<=1)
!
!  Local Variables:
!
      logical :: cut1
      real(REAL8) :: xi1, xi2, zeta1, zeta2
      real(REAL8) :: rlam, rmu
      real(REAL8) :: a, b, c
!
!  Start of Executable Code:
!  ------------------------
!
      cut = .false.
!
!  Coordinates in the (ffp,pmax)-coordinate system
!
      xi1   = dir(1) * ( pos1(2) - pos0(2) ) - dir(2) * ( pos1(1) - pos0(1) )
      xi2   = dir(1) * ( pos2(2) - pos0(2) ) - dir(2) * ( pos2(1) - pos0(1) )
      zeta1 = dir(1) * ( pos1(1) - pos0(1) ) + dir(2) * ( pos1(2) - pos0(2) )
      zeta2 = dir(1) * ( pos2(1) - pos0(1) ) + dir(2) * ( pos2(2) - pos0(2) )
!
!  Check straight lines connecting the semi-ellipses
!
      call mcisct( rlam, rmu, cut1, cut, xi1, zeta1, xi2, zeta2, pmax, ZERO, pmax, ffp )
!          ======   o     o    o     o    i     i     i     i     i     i     i     i
      if( cut ) return
!
      call mcisct( rlam, rmu, cut1, cut, xi1, zeta1, xi2, zeta2, -pmax, ZERO, -pmax, ffp )
!          ======   o     o    o     o    i     i     i     i      i     i      i     i
      if( cut ) return
!
!  Check ellipse around pos0
!
      a = b_a**2 * ( xi2 - xi1 )**2 + ( zeta2 - zeta1 )**2
      b = b_a**2 * xi1 * ( xi2 - xi1 ) + zeta1 * ( zeta2 - zeta1 )
      c = b_a**2 * ( xi1**2 - pmax**2 ) + zeta1**2

      if( a /= 0 ) then

         b = b/a
         c = c/a

         if( c <= b**2 ) then
            rlam = -b + sqrt( b**2 - c )
            rmu = -b - sqrt( b**2 - c )
!!!            print *, '1st ellipse: lambda,mu=', rlam, rmu
         end if

         if( b < -1 ) then
            cut = ( c >= 0 .and. c <= -(2*b+1) )
         else if( b < -0.5 ) then
            cut = ( c >= 0 .and. c <= b**2 )
         else if( b < 0 ) then
            cut = ( c >= -(2*b+1) .and. c <= b**2 )
         else
            cut = ( c >= -(2*b+1) .and. c <= 0 )
         endif
         if( cut ) return

      end if
!
!  Check ellipse around pos0+ffp*dir
!
      zeta1 = zeta1 - ffp
      zeta2 = zeta2 - ffp

      a = b_a**2 * ( xi2 - xi1 )**2 + ( zeta2 - zeta1 )**2
      b = b_a**2 * xi1 * ( xi2 - xi1 ) + zeta1 * ( zeta2 - zeta1 )
      c = b_a**2 * ( xi1**2 - pmax**2 ) + zeta1**2

      if( a /= 0 ) then

         b = b/a
         c = c/a

         if( c <= b**2 ) then
            rlam = -b + sqrt( b**2 - c )
            rmu = -b - sqrt( b**2 - c )
!!!            print *, '2nd ellipse: lambda,mu=', rlam, rmu
         end if

         if( b < -1 ) then
            cut = ( c >= 0 .and. c <= -(2*b+1) )
         else if( b < -0.5 ) then
            cut = ( c >= 0 .and. c <= b**2 )
         else if( b < 0 ) then
            cut = ( c >= -(2*b+1) .and. c <= b**2 )
         else
            cut = ( c >= -(2*b+1) .and. c <= 0 )
         endif
         if( cut ) return

      end if

   end function intersectFFPcyl
!
!=======================================================================
!
end module Geom_FFPcyl_mod
