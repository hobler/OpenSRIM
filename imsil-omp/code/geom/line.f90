module Geom_Line_mod
!
!  Provide a collection of procedures operating on straight lines
!
!  Used Modules:
!  ============
!
   use Const_mod,       only: INT4, REAL8,  &
                              ZERO, ONE, TWO
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
   public :: mcisct                 ! intersection of straight lines
!            ======
   public :: mccut                  ! intersection of straight lines
!            =====                  !    (one parallel to z axis)
   public :: mcdist                 ! distance between point and line
!            ======                 !
   public :: mccutc                 ! intersection of straight line and
!            ======                 !    circle
   public :: mccuti                 ! intersection of straight line with
!            ======                 !    interface
!
contains
!
!=======================================================================
!
   subroutine mcisct( rlam, rmu, cut1, cut2, x1, z1, x2, z2, px1, pz1, px2, pz2, eps_ )
!             ======   o     o    o     o    i   i   i   i    i    i    i    i    i
!
!     Determines if and where the two straight lines (x1,z1)-(x2,z2) and
!     (px1,pz1)-(px2,pz2) are intersecting
!
!     1) uses Cramer method.
!
!     2) as the solution is not so save near 1 as it is near 0, the
!        corresponding equations for (1-rlam) and (1-rmu) are also
!        solved.
!
!     3) it is important to take (as far as possible) the same
!        variables in the present step to determine whether rlam>0
!        as have been taken in the previous step to determine whether
!        rlam<1.
!
!  Arguments:
!
      real(REAL8), intent(out) :: rlam          ! coordinate of intersection point in
                                                !    (x1,z1)-(x2,z2) coordinate system. Only
                                                !    set if there is an intersection point
      real(REAL8), intent(out) :: rmu           ! coordinate of intersection point in
                                                !    (px1,pz1)-(px2,pz2) coordinate system.
                                                !    Only set if there is an intersection point
      logical, intent(out) :: cut1              ! flag indicating there is an intersection point
                                                !    with rlam>=0
      logical, intent(out) :: cut2              ! flag indicating there is an intersection point
                                                !    with 0<=rlam<=1
      real(REAL8), intent(in) :: x1             ! x-coordinate of starting point of first line
      real(REAL8), intent(in) :: z1             ! z-coordinate of starting point of first line
      real(REAL8), intent(in) :: x2             ! x-coordinate of terminating point of first
                                                !    line
      real(REAL8), intent(in) :: z2             ! z-coordinate of terminating point of first
                                                !    line
      real(REAL8), intent(in) :: px1            ! x-coordinate of starting point of second line
      real(REAL8), intent(in) :: pz1            ! z-coordinate of starting point of second line
      real(REAL8), intent(in) :: px2            ! x-coordinate of terminating point of second
                                                !    line
      real(REAL8), intent(in) :: pz2            ! z-coordinate of terminating point of second
                                                !    line
      real(REAL8), optional, intent(in) :: eps_ ! tolerance of rmu, i.e. cut=.true. for
                                                !    -eps_ < rmu <= 1+eps_
!
!  Local Variables:
!
      real(REAL8) :: a11, a21, a12, a22, det, b1, b2, c1, c2, d1, d2
      real(REAL8) :: rlam1, rmu1, eps
!
!  Start of Executable Code:
!  ------------------------
!
! the equation set for rlam and rmu is:
!        a11*rlam + a21*rmu = b1
!        a12*rlam + a22*rmu = b2
! the equation set for rlam and (1-rmu) is:
!        a11*rlam + a21*(rmu-1) = c1
!        a12*rlam + a22*(rmu-1) = c2
! the equation set for (1-rlam) and rmu is analogous:
!        a11*(rlam-1) + a21*rmu = d1
!        a12*(rlam-1) + a22*rmu = d2
!
      a11 = x2 - x1
      a21 = px1 - px2
      a12 = z2 - z1
      a22 = pz1 - pz2

      cut1 = .false.
      cut2 = .false.

      det = a11*a22 - a12*a21

      if( det == ZERO )  return                 ! exit, if the lines are parallel

      b1 = px1 - x1
      b2 = pz1 - z1

      rlam = (b1*a22 - b2*a21) / det
      if( rlam <= ZERO )  return

      if( present( eps_ ) ) then
         eps = eps_
      else
         eps = ZERO
      end if

      rmu = (a11*b2 - a12*b1) / det
      if( rmu <= -eps )  return

      c1 = px2 - x1
      c2 = pz2 - z1

      rmu1 = (a11*c2 - a12*c1) / det            ! rmu1 = rmu - 1
      if( rmu1 > eps )  return

      if( rmu > 0.5 )  rmu = ONE + rmu1
      cut1 = .true.

      d1 = px1 - x2
      d2 = pz1 - z2
      rlam1 = (d1*a22 - d2*a21) / det           ! rlam1 = rlam -1
      if( rlam1 > ZERO )  return

      if( rlam > 0.5)  rlam = ONE + rlam1
      cut2 = .true.

   end subroutine mcisct
!
!=======================================================================
!
   logical function mccut( x1, x2, z, px1, pz1, px2, pz2 )
!                   =====  i   i   i   i    i    i    i
!
!     Determines if two straight lines (one of them parallel to the z axis)
!     are intersecting
!
!  Arguments:
!
      real(REAL8), intent(in) :: x1             ! x-coordinate of starting point of first line
      real(REAL8), intent(in) :: x2             ! x-coordinate of terminating point of first
                                                !    line
      real(REAL8), intent(in) :: z              ! z-coordinate of starting and terminating point
                                                !    of first line
      real(REAL8), intent(in) :: px1            ! x-coordinate of starting point of second line
      real(REAL8), intent(in) :: pz1            ! z-coordinate of starting point of second line
      real(REAL8), intent(in) :: px2            ! x-coordinate of terminating point of second
                                                !    line
      real(REAL8), intent(in) :: pz2            ! z-coordinate of terminating point of second
                                                !    line
!
!  Local Variables:
!
      real(REAL8) :: rmu, rlam
!
!  Start of Executable Code:
!  ------------------------
!
      if( pz1 /= pz2 ) then
         rmu = (z-pz1) / (pz2-pz1)
         rlam = ((px1-x1) + rmu*(px2-px1)) / (x2-x1)
         if( rmu >= ZERO .and. rmu <= ONE .and. rlam >= ZERO .and. rlam <= ONE ) then
            mccut = .true.
         else
            mccut = .false.
         end if
      else
         if( z == pz1 .and. max( px1, px2 ) >= x1 .and. min( px1, px2 ) <= x2 ) then
            mccut = .true.
         else
            mccut = .false.
         end if
      end if

   end function mccut
!
!=======================================================================
!
   subroutine mcdist( distsq, x, z, px1, pz1, px2, pz2, distx_, distz_ )
!             ======    o     i  i   i    i    i    i     o      o
!
!     Calculate the square of the minimum distance between a point (x,z)
!     and a straight line (px1,pz1)-(px2,pz2)
!
!  Arguments:
!
      real(REAL8), intent(out) :: distsq        ! square of the minimum distance between point
                                                !    and line
      real(REAL8), intent(in) :: x              ! x-coordinate of point
      real(REAL8), intent(in) :: z              ! z-coordinate of point
      real(REAL8), intent(in) :: px1            ! x-coordinate of starting point of line
      real(REAL8), intent(in) :: pz1            ! z-coordinate of starting point of line
      real(REAL8), intent(in) :: px2            ! x-coordinate of terminating point of line
      real(REAL8), intent(in) :: pz2            ! z-coordinate of terminating point of line
      real(REAL8), intent(out), optional :: distx_, distz_ ! distance vector
!
!  Local Variables:
!
      real(REAL8) :: px12, pz12, x1, z1, x2, z2
!
!  Start of Executable Code:
!  ------------------------
!
      if( any( [ present( distx_ ), present( distz_ ) ] ) .and.  &
          .not.all( [ present( distx_ ), present( distz_ ) ] ) )  &
            error stop "mcdist: both or none of distx, distz must be present"

      px12 = px1 - px2
      pz12 = pz1 - pz2
      x1 = x - px1
      z1 = z - pz1
      x2 = x - px2
      z2 = z - pz2

      if( -x1*px12 - z1*pz12 <= ZERO ) then
         distsq = x1**2 + z1**2
         if( present( distx_ ) ) then
            distx_ = -x1
            distz_ = -z1
         end if
      else if( px12*x2 + pz12*z2 <= ZERO ) then
         distsq = x2**2 + z2**2
         if( present( distx_ ) ) then
            distx_ = -x2
            distz_ = -z2
         end if
      else
         distsq = (pz12*x1 - px12*z1)**2 / (px12**2 + pz12**2)
         if( present( distx_ ) ) then
            distx_ = px12*(x1*px12+z1*pz12)/(px12**2 + pz12**2)-x1
            distz_ = pz12*(x1*px12+z1*pz12)/(px12**2 + pz12**2)-z1
         end if
      end if

   end subroutine mcdist
!
!=======================================================================
!
   subroutine mccutc( xn, zn, x0, z0, xm, zm, rcirc, cutcir, rlam )
!             ======  i   i   i   i   i   i     i      o      o
!
!     Determines if and where the straight line (xn,zn)-(x0,z0)
!     intersects the circle with center (xm,zm) and radius rcirc
!
!     If there are two intersection points, that one closer to (xn,zn)
!     is taken (smaller rlam).
!
!  Arguments:
!
      real(REAL8), intent(in) :: xn             ! x-coordinate of starting point of line
      real(REAL8), intent(in) :: zn             ! z-coordinate of starting point of line
      real(REAL8), intent(in) :: x0             ! x-coordinate of terminating point of line
      real(REAL8), intent(in) :: z0             ! z-coordinate of terminating point of line
      real(REAL8), intent(in) :: xm             ! x-coordinate of circle center
      real(REAL8), intent(in) :: zm             ! z-coordinate of circle center
      real(REAL8), intent(in) :: rcirc          ! radius of circle
      logical, intent(out) :: cutcir            ! flag indicating if circle intersects line
      real(REAL8), intent(out) :: rlam          ! coordinate of intersection point in
                                                !    (xn,zn)-(x0,z0) coordinate system.
!
!  Local Variables:
!
      real(REAL8) :: xtr, ztr, dirx, dirz, a, b, c, delta, rlam1, rlam2
      real(REAL8) :: dsqrr
!
!  Start of Executable Code:
!  ------------------------
!
!
!  Intersecting a circle with a straight line means solving an equation
!  of second order for rlam:
!
!    ((x0-xm)+rlam*(xn-x0))**2+((x0-xm)+rlam*(xn-x0))**2=rcirc**2
!
!  Initialize the parameter of the equation
!
      cutcir = .false.

      xtr = x0 - xm
      ztr = z0 - zm
      dirx = xn - x0
      dirz = zn - z0
!
!  Employing these notations, the equation is transformed to the form
!  a*rlam**2 + b*rlam + c = 0
!
!  Initialize the parameters of this equation
!
      a = dirx**2 + dirz**2
      b = 2*(xtr*dirx + ztr*dirz)
      c = xtr**2 + ztr**2 - rcirc**2
!
!  Determine the discriminant of the equation
!
      delta = b**2 - 4*a*c
!
!  Different cases for delta to be taken into account
!
      if( delta < ZERO ) then                    ! no intersection point
         cutcir = .false.
         return
      else if( delta == ZERO ) then              ! one intersection point
         rlam1 = - b / (2*a)
         if( rlam1 >= ZERO .and. rlam1 <= ONE ) then
            cutcir = .true.                     ! intersection found
            rlam = rlam1
         end if
         return
      else                                      ! two intersection points
         dsqrr = sqrt(delta)
         rlam1 = (-b+dsqrr) / (2*a)
         rlam2 = (-b-dsqrr) / (2*a)

         if( rlam1 < ZERO .or. rlam1 > ONE ) then
            if( rlam2 < ZERO .or. rlam2 > ONE ) then
               return
            else
               cutcir = .true.
               rlam = rlam2
               return
            end if
         else
            if( rlam2 < ZERO .or. rlam2 > ONE ) then
               cutcir = .true.
               rlam = rlam1
               return
            else
               cutcir = .true.
               if( rlam1 < rlam2 ) then
                  rlam = rlam1
               else
                  rlam = rlam2
               end if
               return
            end if
         end if
      end if

   end subroutine mccutc
!
!=======================================================================
!
   subroutine mccuti( xn, zn, x0, z0, px1, pz1, px2, pz2, dthici, cutint, rlam )
!             ======  i   i   i   i    i    i    i    i      i      o      o
!
!     Determines if and where the straight line (xn,zn)-(x0,z0) intersects
!     the interface layer associated with the line (px1,pz1)-(px2,pz2). The
!     interface layer is the region with distance less than dthici from the
!     line.
!
!  Arguments:
!
      real(REAL8), intent(in) :: xn             ! x-coordinate of starting point of
                                                !    line (xn,zn)-(x0,z0)
      real(REAL8), intent(in) :: zn             ! z-coordinate of starting point of
                                                !    line (xn,zn)-(x0,z0)
      real(REAL8), intent(in) :: x0             ! x-coordinate of terminating point of
                                                !    line (xn,zn)-(x0,z0)
      real(REAL8), intent(in) :: z0             ! z-coordinate of terminating point of
                                                !    line (xn,zn)-(x0,z0)
      real(REAL8), intent(in) :: px1            ! x-coordinate of starting point of
                                                !    line (px1,pz1)-(px2,pz2)
      real(REAL8), intent(in) :: pz1            ! z-coordinate of starting point of
                                                !    line (px1,pz1)-(px2,pz2)
      real(REAL8), intent(in) :: px2            ! x-coordinate of terminating point of
                                                !    line (px1,pz1)-(px2,pz2)
      real(REAL8), intent(in) :: pz2            ! z-coordinate of terminating point of
                                                !    line (px1,pz1)-(px2,pz2)
      real(REAL8), intent(in) :: dthici         ! thickness of the interface layer
      logical, intent(out) :: cutint            ! flag indicating intersection
      real(REAL8), intent(out) :: rlam          ! coordinate of intersection point in the
                                                !    (xn,zn)-(x0,z0) coordinate system.
                                                !    Only set if cutint=.true.
!
!  Local Variables:
!
      logical :: cut1, cut2, cutcir
      integer(INT4) :: ifound
      real(REAL8) :: rlamh, rlamh1, rmu
      real(REAL8) :: dirbx, dirbz, dirnbx, dirnbz, disnor
      real(REAL8) :: pbnd1x, pbnd1z, pbnd2x, pbnd2z
!
!  Start of Executable Code:
!  ------------------------
!
      cut1 = .false.
      cut2 = .false.
      cutcir = .false.

      cutint = .false.
      ifound = 0
!
!  Check, whether the line intersects the circles surrounding the
!  interface layer at the margins.
!
      call mccutc( xn, zn, x0, z0, px1, pz1, dthici, cutcir, rlamh1 )
!          ======  i   i   i   i    i    i     i       o       o

      if( cutcir ) then                          ! first circle intersected by line
         ifound = ifound + 1
         rlamh = rlamh1
         cutcir = .false.
      end if

      call mccutc( xn, zn, x0, z0, px2, pz2, dthici, cutcir, rlamh1 )
!          ======  i   i   i   i    i    i     i       o       o

      if( cutcir ) then                          ! second circle intersected by line
         if( ifound == 0 ) then                  ! first intersection found
            ifound = ifound + 1
            rlamh = rlamh1
         else                                   ! an intersection has previously been found.
            ifound = ifound + 1                 !    Take smaller rlam
            if( rlamh1 < rlamh )  rlamh = rlamh1
         end if
         cutcir = .false.
      end if
!
!  Construct a normal vector of length unity to p1,p2
!
      dirbx = px2 - px1
      dirbz = pz2 - pz1

      dirnbx = - dirbz
      dirnbz = dirbx

      disnor = sqrt( dirnbx**2 + dirnbz**2 )    ! scale to unity length
      dirnbx = dirnbx/disnor
      dirnbz = dirnbz/disnor
!
!  Construct the first boundary line
!
      pbnd1x = px1 + dthici*dirnbx
      pbnd1z = pz1 + dthici*dirnbz
      pbnd2x = px2 + dthici*dirnbx
      pbnd2z = pz2 + dthici*dirnbz
!
!  Cut with line pn,p0
!
      call mcisct( rlamh1, rmu, cut1, cut2, x0, z0, xn, zn, pbnd1x, pbnd1z, pbnd2x, pbnd2z )
!          ======    o      o    o     o    i   i   i   i     i       i       i       i

      if( cut2 ) then                            ! intersect boundary of interface
         if( ifound == 0 ) then                  ! first intersection
            ifound = ifound + 1
            rlamh = rlamh1
         else                                   ! an intersection has previously been found.
            ifound = ifound + 1                 !    Take smaller rlam
            if( rlamh1 < rlamh )  rlamh = rlamh1
         end if
         cut1 = .false.
         cut2 = .false.
      end if
!
!  Construct the second boundary line
!
      pbnd1x = px1 - dthici*dirnbx
      pbnd1z = pz1 - dthici*dirnbz
      pbnd2x = px2 - dthici*dirnbx
      pbnd2z = pz2 - dthici*dirnbz
!
!  Cut with line pn,p0
!
      call mcisct( rlamh1, rmu, cut1, cut2, x0, z0, xn, zn, pbnd1x, pbnd1z, pbnd2x, pbnd2z )
!          ======    o      o    o     o    i   i   i   i     i       i       i       i
      if( cut2 ) then                            ! intersect boundary of interface
         if( ifound == 0 ) then                  ! first intersection
            ifound = ifound + 1
            rlamh = rlamh1
         else                                   ! an intersection has previously been found.
            ifound = ifound + 1                 !    Take smaller rlam
            if( rlamh1 < rlamh )  rlamh = rlamh1
         end if
         cut1 = .false.
         cut2 = .false.
      end if
!
! Check whether there is an intersection and set output parameters
!
      if( ifound >= 1 ) then
         rlam = rlamh
         cutint = .true.
      endif

   end subroutine mccuti
!
!=======================================================================
!
end module Geom_Line_mod
