module Geom_Polygon_Old_mod
!
!  Provide a collection of procedures operating on a closed polygon
!
!  Used Modules:
!  ============
!
   use Const_mod,       only: INT4, REAL8,  &
                              ZERO, ONE, TWO, HALF, PI, BIG, EPSMA, EPSMB
   use Open_mod,        only: openFile

   use Geom_Line_mod,   only: mccut, mcdist, mcisct
   use Math_Sort_mod,   only: quickSort
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
   public :: mcauxl                 ! calculate direction of auxiliary line
!            ======
   public :: mcloc                  ! location of a point relative to polygon
!            =====
   public :: mccut4                 ! intersection between rectangle (parallel to axes)
!            ======                 !    and polygon
   public :: mcdis4                 ! minimum distance between rectangle (parallel to axes)
!            ======                 !    and polygon
   public :: getDistPolygon         ! get the minimum distance of a point from a polygon
!            ==============
   public :: getEntrancePolygon     ! get the entrance point into a polygon given an
!            ==================     !    initial position and direction
   public :: getPolyArea            ! calculate the area of the polygon
!            ===========
   public :: relPosPolygon          ! determine the relative position of a point with respect
!            =============          !    to a polygon
   public :: pointInsidePolygon     ! return point safely inside polygon
!            ==================
   public :: getDistPolyPoly        ! maximum and minimum of the minimum distance of any point
!            ===============        !    of one polygon from all points of another polygon
   public :: writePolygons          ! write the surface polygon to the output file
!            =============          !    as passed to the subroutine
!
!  Private Data:
!  ============
!
!  None
!
contains
!
!=======================================================================
!
   subroutine mcauxl( coseta, sineta, px, pz )
!             ======    o       o     i   i
!
!     Calculates the direction of the auxiliary line for a given polygon.
!     It is calculated as the direction which is most different from any
!     direction of the polygon lines.
!
!  Arguments:
!
      real(REAL8), intent(out) :: coseta        ! cosine of direction angle of auxiliary line
      real(REAL8), intent(out) :: sineta        ! sine of direction angle of auxiliary line
      real(REAL8), intent(in) :: px(:)          ! x-coordinates of polygon corners
      real(REAL8), intent(in) :: pz(:)          ! z-coordinates of polygon corners
!
!  Local Variables:
!
      integer(INT4) :: np, ip, jp, iip
      real(REAL8) :: dpx, dpz, etamin, h, detamx, deta, etaaux
      real(REAL8) :: eta(size( px ))

!
!  Start of Executable Code:
!  ------------------------
!
      np = size( px ) - 1
!
!  Calculate direction angles of the boundary lines
!
      do ip = 1, np
         dpx = px(ip+1) - px(ip)
         dpz = pz(ip+1) - pz(ip)
         eta(ip) = atan2( dpz, dpx )
         if( eta(ip) < ZERO )  eta(ip) = eta(ip) + PI
      end do
!
!  Sort direction angles
!
      do ip=1,np-1
         jp = ip
         etamin = eta(ip)
         do iip = ip+1, np
            if( eta(iip) < etamin ) then
               jp = iip
               etamin = eta(iip)
            end if
         end do

         h = eta(ip)
         eta(ip) = eta(jp)
         eta(jp) = h
      end do

      eta(np+1) = eta(1) + PI
!
!  Search maximum difference between adjacent direction angles
!
      detamx = eta(2) - eta(1)
      jp = 1

      do ip = 2, np
         deta = eta(ip+1) - eta(ip)
         if( deta > detamx ) then
            jp = ip
            detamx = deta
         end if
      end do
!
!  The direction angle of the auxiliary line is approximately the mean
!  value of the direction angles of the two boundary lines with the
!  maximum difference in the direction angles (see above)
!
      etaaux = HALF* (eta(jp+1) + eta(jp)) + detamx/99.
      coseta = cos( etaaux )
      sineta = sin( etaaux )

   end subroutine mcauxl
!
!=======================================================================
!
   integer(INT4) function mcloc( x, z, px, pz, coseta, sineta )
!                         =====  i  i  i   i     i       i
!
!     Determine location of point (x,z) relative to a closed polygon:
!
!        mcloc = +1 if the point is inside the polygon
!                 0 if the point is on the polygon
!                -1 if the point is outside or on the polygon
!
!     The information about the location of the point is obtained by
!     intersecting an auxiliary line, which contains the point, with
!     all lines of the polygon, and counting the number of times the
!     auxiliary line is leaving or entering the inner part of the
!     polygon (in the direction off the point). It is assumed the auxiliary
!     line is not nearly parallel to one of the polygon lines.
!
!  Arguments:
!
      real(REAL8), intent(in) :: x              ! x-coordinate of point
      real(REAL8), intent(in) :: z              ! z-coordinate of point
      real(REAL8), intent(in) :: px(:)          ! x-coordinates of polygon corners (first=last)
      real(REAL8), intent(in) :: pz(:)          ! z-coordinates of polygon corners (first=last)
      real(REAL8), intent(in) :: coseta         ! cosine of direction angle of auxiliary line
      real(REAL8), intent(in) :: sineta         ! sine of  direction angle of auxiliary line
!
!  Local Variables:
!
      integer(INT4) :: np, nin, nout, ip, ip0
      real(REAL8) :: dpx, dpz, dx, dz, det, rmu
      real(REAL8) :: dpx0, dpz0, rlam
      real(REAL8) :: dpx1, dpz1, dx1, dz1, det1, rmu1
!
!  Start of Executable Code:
!  ------------------------
!
      np = size( px ) - 1

      nin = 0
      nout = 0

      do ip = 1, np

         dpx = px(ip+1) - px(ip)
         dpz = pz(ip+1) - pz(ip)
         dx = x - px(ip)
         dz = z - pz(ip)

         det =  - dpz*coseta + dpx*sineta
         rmu = - dz*coseta + dx*sineta

         dpx1 = -dpx
         dpz1 = -dpz
         dx1  = x - px(ip+1)
         dz1  = z - pz(ip+1)

         det1 =  - dpz1*coseta + dpx1*sineta
         rmu1 = - dz1*coseta + dx1*sineta
!
!  rmu=0 ... corner point lies on auxiliary line
!    dx=0 .and. dz=0 ... (x,z) is identical with the corner point
!    the next block-if takes the inner product of the interconnection
!      of the corner point and (x,z)  with  the normal vectors of the
!      two lines adjacent to the corner point
!      with that, it is decided whether the auxiliary line is leaving,
!      entering, or grazing the region
!
!  0<=rmu/det<1 ... there is an intersection between the auxiliary line
!                  and the boundary line within the corner points
!    the next block-if takes the inner product of the interconnection
!      of the corner point and (x,z) with the normal vector of the
!      boundary line
!      with that, it is decided whether the auxiliary line is leaving
!      or entering the region, or if( x,z) lies on the straight line
!
!  note: det/=0 is guaranteed, because the auxiliary line is not
!        chosen parallel to a boundary line
!
         if( rmu == ZERO ) then
            if( dx == ZERO .and. dz == ZERO ) then
               mcloc = 0
               return
            end if

            ip0 = ip - 1
            if( ip0 == 0 )  ip0 = np
            dpx0 = px(ip) - px(ip0)
            dpz0 = pz(ip) - pz(ip0)

            if( dpz*dx-dpx*dz > ZERO .and. dpz0*dx-dpx0*dz > ZERO ) then
               nout = nout + 1
            else if( dpz*dx-dpx*dz < 0 .and. dpz0*dx-dpx0*dz < 0 ) then
               nin = nin + 1
            end if
         else if( rmu*det >= ZERO .and. rmu1*det1 > ZERO ) then
            rlam = dpz*dx - dpx*dz
            if( rlam > ZERO ) then
               nout = nout + 1
            else if( rlam < ZERO ) then
               nin = nin + 1
            else
               mcloc = 0
               return
            end if
         end if

      end do
!
!  The difference between nout and nin determines whether (x,z) is
!    inside or outside the (closed) polygon
!  Take the absolute value, because the boundary lines may be both right-handed
!    and left-handed
!
      nout = abs( nout - nin )

      if( nout == 0 ) then
         mcloc = -1
      else if( nout == 2 ) then
         mcloc = +1
      else
         mcloc = 0
         print *, 'x,z=', x,z
         do ip = 1, size( px )
            print *, 'ip,x,y=', ip, px(ip), pz(ip)
         end do
         error stop 'error in function mcloc'
      end if

      end function mcloc
!
!=======================================================================
!
   logical function mccut4( x1, x2, z1, z2, px, pz )
!                   ======  i   i   i   i   i   i
!
!     Determines if a rectangle (parallel to axes) intersects a polygon
!
!     mccut4 = .true. if rectangle intersects polygon
!              .false. if rectangle does  n o t  intersect polygon
!
!  Arguments:
!
      real(REAL8), intent(in) :: x1             ! left boundary of rectangle
      real(REAL8), intent(in) :: x2             ! right boundary of rectangle
      real(REAL8), intent(in) :: z1             ! lower boundary of rectangle
      real(REAL8), intent(in) :: z2             ! upper boundary of rectangle
      real(REAL8), intent(in) :: px(:)          ! x-coordinates of polygon corners
                                                !    (px(np)=px(1))
      real(REAL8), intent(in) :: pz(:)          ! z-coordinates of polygon corners
                                                !    (pz(np)=pz(1))
!
!  Local Variables:
!
      integer(INT4) :: np, ip
!
!  Start of Executable Code:
!  ------------------------
!
      np = size( px ) - 1
!
      mccut4 = .false.

      do ip = 1, np
         mccut4 = mccut4 .or. mccut( x1, x2, z1, px(ip), pz(ip), px(ip+1), pz(ip+1) )
!                             =====  i   i   i    i       i       i         i
         if( mccut4 )  return
      end do

      do ip = 1, np
         mccut4 = mccut4 .or. mccut( x1, x2, z2, px(ip), pz(ip), px(ip+1), pz(ip+1) )
!                             =====  i   i   i    i       i       i         i
         if( mccut4 )  return
      end do

      do ip = 1, np
         mccut4 = mccut4 .or. mccut( z1, z2, x1, pz(ip), px(ip), pz(ip+1), px(ip+1) )
!                             =====  i   i   i    i       i       i         i
         if( mccut4 )  return
      end do

      do ip = 1, np
         mccut4 = mccut4 .or. mccut( z1, z2, x1, pz(ip), px(ip), pz(ip+1), px(ip+1) )
!                             =====  i   i   i    i       i       i         i
         if( mccut4 )  return
      end do

   end function mccut4
!
!=======================================================================
!
   subroutine mcdis4( distminsq, distmaxsq, inside,  &
!             ======      o          o        o
                      x1, x2, z1, z2, rsq, px, pz, reentr, coseta, sineta)
!                     i   i   i   i    i   i   i     i       i       i
!
!     Calculates the square of the minimum and maximum distance of points inside a rectangle
!     (parallel to axes) from a polygon (px,pz).
!
!     If rectangle and polygon are intersecting, or the whole polygon is
!     inside the rectangle, distminsq = 0.
!
!     NOTE: By "distance" between a point and a polygon we mean the minimum distance between
!           the point and any point on the polygon.
!
!  Arguments:
!
      real(REAL8), intent(out) :: distminsq     ! minimum distance of any point inside the
                                                !    rectangle from the polygon
      real(REAL8), intent(out) :: distmaxsq     ! maximum distance of any point inside the
                                                !    rectangle from the polygon
      logical, intent(out) :: inside            ! flag indicating center of rectangle is inside
                                                !    the polygon
      real(REAL8), intent(in) :: x1             ! left boundary of rectangle
      real(REAL8), intent(in) :: x2             ! right boundary of rectangle
      real(REAL8), intent(in) :: z1             ! lower boundary of rectangle
      real(REAL8), intent(in) :: z2             ! upper boundary of rectangle
      real(REAL8), intent(in) :: rsq            ! radius of circle circumscribing rectangle
                                                !     = 0.25 * ((x2-x1)**2 + (z2-z1)**2)
      real(REAL8), intent(in) :: px(:)          ! x-coordinates of polygon corners (first=last)
      real(REAL8), intent(in) :: pz(:)          ! z-coordinates of polygon corners (first=last)
      logical, intent(in) :: reentr(:)          ! flag indicating reentrant corners
      real(REAL8), intent(in) :: coseta         ! cosine of direction angle of auxiliary line
      real(REAL8), intent(in) :: sineta         ! sine of direction angle of auxiliary line
!
!  Local Variables:
!
      logical :: cuting
      integer(INT4) :: np, ip

      real(REAL8) :: dx, dz, xm, zm, dmsqmin, distsq, dsqmin, dsqmax
      real(REAL8) :: dmsq(size( px )-1)         ! distance of a boundary line from the middle
                                                !    of the rectangle
!
!  Start of Executable Code:
!  ------------------------
!
      np = size( px ) - 1
!
!  Calculate the square of the distance between the center of the rectangle and the polygon
!
      dx = x2 - x1
      dz = z2 - z1
      xm = HALF * (x1+x2)
      zm = HALF * (z1+z2)
      call mcdist( dmsq(1), xm, zm, px(1), pz(1), px(2), pz(2) )
!          ======    o      i   i    i      i      i      i
      dmsqmin = dmsq(1)
      do ip = 2, np
         call mcdist( dmsq(ip), xm, zm, px(ip), pz(ip), px(ip+1), pz(ip+1) )
!             ======    o       i   i    i       i       i         i
         dmsqmin = min( dmsqmin, dmsq(ip) )
      end do
!
!  Determine whether the center of the rectangle is inside or outside the polygon
!
      inside = ( mcloc( xm, zm, px(1:np+1), pz(1:np+1), coseta, sineta ) > 0 )
!                =====  i   i      i           i          i       i
!
!  If this distance is greater than half of the length of the diagonal of
!  the rectangle, the rectangle and the polygon may not intersect.
!  If this distance is less than half of the length of the shorter line of
!  the rectangle, the rectangle and the polygon intersect in any case
!
      if( dmsqmin > rsq ) then
         cuting = .false.
      else if( dmsqmin < (HALF*min( dx, dz ))**2 ) then
         cuting = .true.
      else
         cuting = mccut4( x1, x2, z1, z2, px(1:np+1), pz(1:np+1) )
!                 ======  i   i   i   i      i           i
      end if
!
!  Calculate distminsq.
!  If polygon and rectangle are intersecting, the distance is 0.
!  If they are not intersecting, the rectangle may either be totally
!    inside the polygon or totally outside.
!  If it is totally inside, the minimum distance is given by one
!    of the corners of the rectangle and a line of the the polygon
!    or by a reentrant corner of the polygon and a line of the rectangle
!
      ! the next line calcuates a close upper limit for distsq
      ! (a+b)**2 <= 2*a**2 + 2*b**2.
      ! Explanation: all points of the rectangle lie within a
      ! distance of r from the middle of the rectangle (rsq=r**2)
      ! Therefore no point of the rectangle can have a larger
      ! distance from the boundary than dmmn+r (dmsqmin=dmmn**2).
      ! a=dmmn, b=r above.
      dmsqmin = (dmsqmin + rsq) * TWO
      distminsq = dmsqmin
      distmaxsq = dmsqmin

      do ip = 1, np
         if( dmsq(ip) > dmsqmin )  cycle

         call mcdist( distsq, x1, z1, px(ip), pz(ip), px(ip+1), pz(ip+1) )
!             ======    o     i   i    i       i       i         i
         dsqmin = distsq
         dsqmax = distsq
         call mcdist( distsq, x1, z2, px(ip), pz(ip), px(ip+1), pz(ip+1) )
!             ======    o     i   i    i       i       i         i
         dsqmin = min( dsqmin, distsq )
         dsqmax = max( dsqmax, distsq )
         call mcdist( distsq, x2, z1, px(ip), pz(ip), px(ip+1), pz(ip+1) )
!             ======    o     i   i    i       i       i         i
         dsqmin = min( dsqmin, distsq )
         dsqmax = max( dsqmax, distsq )
         call mcdist( distsq, x2, z2, px(ip), pz(ip), px(ip+1), pz(ip+1) )
!             ======    o     i   i    i       i       i         i
         dsqmin = min( dsqmin, distsq )
         dsqmax = max( dsqmax, distsq )

         distminsq = min( distminsq, dsqmin )
         distmaxsq = min( distmaxsq, dsqmax )

         if( inside ) then
            if( .not.reentr(ip) )  cycle
         else
            if( reentr(ip) )  cycle
         end if

         call mcdist( distsq, px(ip), pz(ip), x1, z1, x1, z2 )
!             ======    o      i       i      i   i   i   i
         distminsq = min( distminsq, distsq )
         call mcdist( distsq, px(ip), pz(ip), x1, z1, x2, z1 )
!             ======    o      i       i      i   i   i   i
         distminsq = min( distminsq, distsq )
         call mcdist( distsq, px(ip), pz(ip), x1, z2, x2, z2 )
!             ======    o      i       i      i   i   i   i
         distminsq = min( distminsq, distsq )
         call mcdist( distsq, px(ip), pz(ip), x2, z1, x2, z2)
!             ======    o      i       i      i   i   i   i
         distminsq = min( distminsq, distsq )
      end do

      if( cuting )  distminsq = ZERO

   end subroutine mcdis4
!
!==========================================================================================
!
   subroutine getDistPolygon( distsq, x, z, px, pz, virtual_ )
!             ==============    o     i  i  i   i     i
!
!     Get the distance of point (x,z) from the polygon (px,pz)
!
!  Arguments:
!
      real(REAL8), intent(out) :: distsq           ! square of the distance of point (x,z) from
                                                   !    polygon (px,pz)
      real(REAL8), intent(in) :: x                 ! x-coordinate of point
      real(REAL8), intent(in) :: z                 ! z-coordinate of point
      real(REAL8), intent(in) :: px(:)             ! x-coordinates of polygon corners
                                                   !    (first=last)
      real(REAL8), intent(in) :: pz(:)             ! z-coordinates of polygon corners
                                                   !    (first=last)
      logical, optional, intent(in) :: virtual_(:) ! flag indicating virtual (non-existent)
                                                   !    segments
!
!  Local Variables:
!
      integer(INT4) :: ip, np
      real(REAL8) :: dist2
!
!  Start of Executable Code:
!  ------------------------
!
      np = size( px ) - 1

      distsq = BIG

      do ip = 1, np
         if( present( virtual_ ) ) then
            if( virtual_(ip) )  cycle
         end if
         call mcdist( dist2, x, z, px(ip), pz(ip), px(ip+1), pz(ip+1) )
!             ======    o    i  i   i       i       i         i
         distsq = min( dist2, distsq )
      end do

   end subroutine getDistPolygon
!
!=======================================================================
!
   subroutine getEntrancePolygon( xentr, zentr, cosinc, found,  &
!             ==================    o      o      o       o
                                  x, z, dirx, dirz, px, pz, lfront, rmin_, rmax_ )
!                                 i  i   i     i    i   i     i       i     i
!
!     Get the entrance point into the polygon (px,pz) along the line defined by (x,z) and
!     (dirx,dirz).
!
!  Arguments:
!
      real(REAL8), intent(out) :: xentr         ! x-coordinate of entrance point
      real(REAL8), intent(out) :: zentr         ! z-coordinate of entrance point
      real(REAL8), intent(out) :: cosinc        ! cosine of the incidence angle at entrance
                                                !    point
      logical, intent(out) :: found             ! flag indicating whether an entrance point
                                                !    has been found
      real(REAL8), intent(in) :: x              ! x-coordinate of point on line
      real(REAL8), intent(in) :: z              ! z-coordinate of point on line
      real(REAL8), intent(in) :: dirx           ! x-component of direction vector
      real(REAL8), intent(in) :: dirz           ! z-component of direction vector
                                                !    NOTE: dirx**2+dirz**2=1 not required,
                                                !          but dirx**2+dirz**2/=0.
      real(REAL8), intent(in) :: px(:)          ! x-ccordinates of polygon corner (first=last)
      real(REAL8), intent(in) :: pz(:)          ! z-ccordinates of polygon corner (first=last)
      logical, intent(in) :: lfront             ! flag indicating whether only entrance points
                                                !    in front of (x,z) shall be considered
      real(REAL8), intent(in), optional :: rmin_(2)   ! lower bounds of limiting box
      real(REAL8), intent(in), optional :: rmax_(2)   ! upper bounds of limiting box
!
!  Local Variables:
!
      integer(INT4) :: ip, jp, np
      real(REAL8) :: a11, a12, a21, a22, b1, b2, c1, c2, a22j, a21j
      real(REAL8) :: det, detj
      real(REAL8) :: rlam, rlamj, rmu, rmu1, rmuj
!
!  Start of Executable Code:
!  ------------------------
!
!  Note: This is not actually necessary while the contour can be limited only at
!        one side but for the moment this is the fastest solution
!
      if( any( [ present( rmin_ ), present( rmax_ ) ] ) .and.  &
          .not.all( [ present( rmin_ ), present( rmax_ ) ] ) )  &
            error stop 'getEntrancePolygon: both or none of rmin, rmax must be present'

      np = size( px ) - 1

      jp = 0
      found = .false.

      a11 = dirx
      a12 = dirz

      do ip = 1, np
!
!  Do not consider line segments completely outside the limiting box
!
         if( present( rmin_ ) ) then
            if( ( px(ip)<rmin_(1) .or. px(ip)>rmax_(1) .or.  &
                  pz(ip)<rmin_(2) .or. pz(ip)>rmax_(2) )  &
                .or.  &   ! or maybe .and. but than the entrance point can be outside the target
                 ( px(ip+1)<rmin_(1) .or. px(ip+1)>rmax_(1) .or. pz(ip+1)<rmin_(2)  &
                   .or. pz(ip+1)>rmax_(2) ) )  cycle
         end if

         a21 = - px(ip+1) + px(ip)
         a22 = - pz(ip+1) + pz(ip)
         det = a11*a22 - a12*a21
         if( det == ZERO )  cycle

         b1 = px(ip) - x
         b2 = pz(ip) - z

         rlam = (b1*a22 - b2*a21) / det
         if( lfront .and. rlam <= ZERO )  cycle

         rmu = (a11*b2 - a12*b1) / det
         if( rmu <= ZERO )  cycle

         c1 = px(ip+1) - x
         c2 = pz(ip+1) - z

         rmu1 = (a11*c2 - a12*c1) / det
         if( rmu1 > ZERO )  cycle
         if( rmu > HALF)  rmu = ONE + rmu1
!
! There is an intersection point !
!
         if( .not.found ) then
            found = .true.
         else
            if( rlam >= rlamj )  cycle
         end if

         jp = ip
         rlamj = rlam
         rmuj = rmu
         a22j = a22
         a21j = a21
         detj = det

      end do
!
!  If no intersection point has been found, return
!
      if( .not.found ) return

      xentr = px(jp) + rmuj * ( px(jp+1) - px(jp) )
      zentr = pz(jp) + rmuj * ( pz(jp+1) - pz(jp) )

                                                ! normal vector to px_r,pz_r has direction
                                                !    (a22,-a21). cosinc is cosine of incidence
                                                !    angle with respect to normal vector
      cosinc = abs(detj) / sqrt( a22j**2 + a21j**2 )

   end subroutine getEntrancePolygon
!
!=======================================================================
!
   function getPolyArea( x, z, np )   result( area )
!           ===========  i  i  i               o
!
!     Code converted using TO_F90 by Alan Miller
!     given a sequence of np points (x(1:np),z(1:np)), polyarea computes the area
!     bounded by the closed polygonal curve which passes through the points in
!     the order that they are indexed. The points are assumed to defined a closed
!     polygon. Therefore, the first point need not be listed at the end
!     of x and y. The curve is not required to be simple e.g. It may cross over
!     itself.
!
!     The area is positive if the vertices of the polygon are ordered counter
!     clockwise and negative for clockwise ordering.
!
!  Arguments:
!
      real(REAL8), intent(in) :: x(:)       ! x-coordinates of the points
      real(REAL8), intent(in) :: z(:)       ! z-
      integer(INT4), intent(in) :: np       ! number of points
!
!  Local Variables:
!
      real(REAL8) :: area   ! resulting area
      integer(INT4) :: i, n, nm1
      real(REAL8) :: a
!
!  Start of Executable Code:
!  ------------------------
!
      n = np
      if( x(1) == x(n) .and. z(1) == z(n) )  n = n - 1

      select case( n )
         case( :2 )
            area = ZERO
         case( 3 )
            area = HALF * ( (x(2) - x(1))*(z(3) - z(1)) - (x(3) - x(1))*(z(2) - z(1)) )
         case default
            nm1 = n - 1
            a = x(1)*( z(2) - z(n) ) + x(n)*( z(1) - z(nm1) )
            do  i = 2, nm1
               a = a + x(i)*( z(i+1) - z(i-1) )
            end do
            area = HALF * a
      end select

   end function getPolyArea
!
!=========================================================================================
!
!      subroutine getCellAreas_try( px, pz, ox, oz, wx, wz, area )
!!                ================  i   i   i   i   i   i    o
!!
!!     Obtain the fractional cell area inside the given polygon.
!!
!!  Arguments:
!!
!      real(REAL8), intent(in) :: px(:)          ! x-coordinates of polygon corners
!      real(REAL8), intent(in) :: pz(:)          ! z-coordinates of polygon corners
!      real(REAL8), intent(in) :: ox             ! origin of cell array (smallest x-value)
!      real(REAL8), intent(in) :: oz             ! origin of cell array (smallest z-value)
!      real(REAL8), intent(in) :: wx             ! cell width in x-direction
!      real(REAL8), intent(in) :: wz             ! cell width in x-direction
!      real(REAL8), intent(inout) :: area(:,:)   ! area inside polygon
!                                                !    first index runs in x-direction
!                                                !    second index runs in z-direction
!!
!!  Local Variables:
!!
!      type :: Point_type
!         real(REAL8) :: x
!         real(REAL8) :: z
!      end type
!      type(Point_type) :: pp(size( px )-1)      ! polygon points in scaled coordinates
!!
!      type(Point_type) :: pa(size( px )+3)      ! polygon for area calculation
!!
!      type(Point_type), allocatable :: pi(:)    ! intersection points
!      real(REAL8), allocatable :: di(:)         ! corresponding distances along polygon line
!      real(REAL8), allocatable :: isort(:)      ! index for ordering according to d
!      logical, allocatable :: lx(:)             ! lx=.true.:  x-value is on grid line
!                                                ! lx=.false.: z-value is on grid line
!!
!      integer(INT4) :: ip, np
!      integer(INT4) :: icx, icz, icx1, icz1, icx2, icz2, ncx, ncz
!      integer(INT4) :: ix1, ix2, ixi, iz1, iz2, izi, nxi, nzi, ni, ii, ix, iz
!      integer(INT4) :: ia
!!
!      real(REAL8) :: px1, pz1, px2, pz2, dx, dz, d
!!
!!  Start of Executable Code:
!!  ------------------------
!!
!!  Convert px, pz into polygon points and scale with cell width
!!
!      np = size( px ) - 1
!!
!      do ip = 1, np
!         pp(ip)%x = ( px(ip) - ox ) / wx
!         pp(ip)%z = ( pz(ip) - oz ) / wz
!      end do
!!
!!  Cells are uniquely defined by the maximum index
!!
!      ncx = size( area, 1 )
!      ncz = size( area, 2 )
!!
!      area(:,:) = ZERO
!!
!!  Check that all polygon points are inside grid
!!
!      do ip = 1, np
!         icx = ceiling( pp(ip)%x )
!         icz = ceiling( pp(ip)%z )
!         if( icx < 1 .or. icx > ncx .or. icz < 1 .or. icz > ncz ) then
!            error stop 'getCellAreas: Polygon point outside grid'
!         end if
!      end do
!!
!!  Follow polygon as it is given by its points
!!
!      do ip = 1, np
!!
!!  Cell indices of start and end point and length of line
!!
!         px1 = pp(ip)%x
!         pz1 = pp(ip)%z
!         px2 = pp( inext(ip,np) )%x
!         pz2 = pp( inext(ip,np) )%z
!         icx1 = ceiling( px1 )
!         icz1 = ceiling( pz1 )
!         icx2 = ceiling( px2 )
!         icz2 = ceiling( pz2 )
!!
!         dx = px2 - px1
!         dz = pz2 - pz1
!         d = sqrt( dx**2 + dz**2 )
!!
!!  Determine intersection points with grid
!!
!         if( icx2 >= icx1 ) then
!            ix1 = icx1 + 1
!            ix2 = icx2
!            ixi = 1
!            nxi = icx2 - icx1
!         else
!            ix1 = icx1
!            ix2 = icx2 + 1
!            ixi = -1
!            nxi = icx1 - icx2
!         end if
!!
!         if( icz2 >= icz1 ) then
!            iz1 = icz1 + 1
!            iz2 = icz2
!            izi = 1
!            nzi = icz2 - icz1
!         else
!            iz1 = icz1
!            iz2 = icz2 + 1
!            izi = -1
!            nzi = icz1 - icz2
!         end if
!!
!         ni = nxi + nzi
!!
!         allocate( pi(ni), di(ni), isort(ni), lx(ni) )
!         ii = 0
!!
!         do ix = ix1, ix2, ixi
!            ii = ii + 1
!            lx(ii) = .true.
!            pi(ii)%x = ix
!            pi(ii)%z = pz1 + dz/dx*(ix-px1)
!            di(ii) = (ix-px1)/dx * d
!         end do
!!
!         do iz = iz1, iz2, izi
!            ii = ii + 1
!            lx(ii) = .false.
!            pi(ii)%z = iz
!            pi(ii)%x = px1 + dx/dz*(iz-pz1)
!            di(ii) = (iz-pz1)/dz * d
!         end do
!!
!!  Sort the points according to their distance from the first polygon line point
!!
!         call quickSort( di(:), isort(:), 1, size( di ) )
!!             =========   io      io      i    i
!         pi(:) = pi(isort(:))
!         di(:) = di(isort(:))
!         lx(:) = lx(isort(:))
!!
!!  Determine on which cell border the first point is located, using the following
!!  counting rule
!!                              |  4  |
!!                            --+-----+--      z
!!                              |     |       ^
!!                             1|     |3      |
!!                              |     |       |
!!                            --+-----+--     +---> x
!!                              |  2  |
!!
!         ii = 1    ! to be modified
!         ia = 1    !
!!
!         if( ii < ni ) then
!            ii = ii + 1
!            ia = ia + 1
!            pa(ia) = pi(ii)
!            ! close polygon with cell lines
!            ! calculate fractional cell area
!            ! start new polygon
!         else
!            pa(ia) = pp( inext(.,.) )
!            ! step forward to new polygon line
!            ! calculate new intersection points
!         end if
!
!
!         deallocate( pi, di, isort, lx )
!
!      contains
!!
!!-----------------------------------------------------------------------------------------
!!
!         integer(INT4) function inext( i, n )
!!                               =====  i  i
!!
!!        Provide next index in a cyclic way
!!
!!     Arguments:
!!
!         integer(INT4), intent(in) :: i         ! index to be incremented
!         integer(INT4), intent(in) :: n         ! maximum index
!!
!!     Local Variables:
!!
!!        None
!!
!!     Start of Executable Code:
!!     ------------------------
!!
!         if( i < n ) then
!            inext = i + 1
!         else
!            inext = 1
!         end if
!!
!         end function inext
!!
!!-----------------------------------------------------------------------------------------
!!
!         integer(INT4) function iprev( i, n )
!!                               =====  i  i
!!
!!        Provide previous index in a cyclic way
!!
!!     Arguments:
!!
!         integer(INT4), intent(in) :: i         ! index to be decremented
!         integer(INT4), intent(in) :: n         ! maximum index
!!
!!     Local Variables:
!!
!!        None
!!
!!     Start of Executable Code:
!!     ------------------------
!!
!         if( i > 1 ) then
!            iprev = i - 1
!         else
!            iprev = n
!         end if
!!
!         end function iprev
!!
!!-----------------------------------------------------------------------------------------
!!
!      end subroutine getCellAreas_try
!
!=========================================================================================
!
   subroutine getCellAreas_old( px, pz, ox, oz, wx, wz, area )
!             ================  i   i   i   i   i   i    o
!
!     Obtain the fractional cell area inside the given polygon.
!     Currently resticted to maximum 4-intersection points between the polygon and the
!     the cell. Not more than one polygon point per cell is allowed.
!
!  Arguments:
!
      real(REAL8), intent(in) :: px(:)          ! x-coordinates of polygon corners
      real(REAL8), intent(in) :: pz(:)          ! z-coordinates of polygon corners
      real(REAL8), intent(in) :: ox             ! origin of cell array (smallest x-value)
      real(REAL8), intent(in) :: oz             ! origin of cell array (smallest z-value)
      real(REAL8), intent(in) :: wx             ! cell width in x-direction
      real(REAL8), intent(in) :: wz             ! cell width in x-direction
      real(REAL8), intent(inout) :: area(:,:)   ! area inside polygon
!
!  Local Variables:
!
      type :: cellPts
!             =======
         integer(INT4) :: np     ! number of points in cell
                                 ! (input polygon points + intersection points
         integer(INT4) :: nisct  ! number of intersection points in the cell
                                 ! input polygon intersects the cell edges
         real(REAL8) :: x(6)     ! coordinates of polygon
         real(REAL8) :: z(6)     ! intersection points in the cell
      end type cellPts
!
      type(cellPts), allocatable :: cell(:,:)
      integer(INT4) :: nx, nz, ixp1, ixp2, izp1, izp2, ix, iz, ix1, iz1, ip, ncins, ncp
      integer(INT4) :: dix, diz
      integer(INT4) :: np1, status, icp, ip1, ip2, np, isct, icx1, icx2, icz1, icz2
      real(REAL8) :: x1, x2, z1, z2, rlam1, rlam2, ascale
      real(REAL8) :: xtmp, ztmp, xisct, zisct, coseta, sineta, cx(4), cz(4), cpx(8), cpz(8)
      logical :: cut1, cut2, lsub
!
!  Start of Executable Code:
!  ------------------------
!
      np = size( px )
!
      area(:,:) = ZERO
!
!  Allocate and initialize auxiliary cell array with same shape as area(:,:)
!
      allocate( cell( size( area, 1 ), size( area, 2 ) ), stat=status )
      if( status /= 0 )  error stop 'getCellAreas: cell not allocated'

      cell%np = 0
      cell%nisct = 0
      nx = size( area, 1 )
      nz = size( area, 2 )
      do ix = 1, nx
         do iz = 1, nz
            cell(ix,iz)%x = -BIG
            cell(ix,iz)%z = -BIG
         end do
      end do
!
!  Initialize coseta and sineta. Required for mcloc
!
      call mcauxl( coseta, sineta, px, pz )
!          ======    o        o     i   i

      ascale = wx * wz                      ! area scaling factor
      ix1 = ceiling( ( px(1) - ox ) / wx )  ! Detect the cell where the first
      iz1 = ceiling( ( pz(1) - oz ) / wz )  ! point located is
      isct = 0
!
!  Follow the polygon as it is given by its point.
!
      do ip = 1, np-1
!
!  Check intersection with horizontal (z=const) lines
!
         izp1 = ceiling( ( pz(ip) - oz ) / wz )
         izp2 = ceiling( ( pz(ip+1) - oz ) / wz )
         diz = nint( sign( ONE, real(izp2 - izp1, REAL8) ) )

         if( diz > 0 ) then
            icz1 = izp1
            icz2 = izp2-1
         else
            icz1 = izp1-1
            icz2 = izp2
         end if

         do iz = icz1, icz2, diz
            x1 = ox
            x2 = ox + nx * wx
            z1 = oz + iz * wz

            call mcisct( rlam1, rlam2, cut1, cut2, x1, z1, x2, z1,  &
!                ======    o      o     o     o    i   i   i   i
                         px(ip), pz(ip), px(ip+1), pz(ip+1) )
!                          i       i         i         i

            if( .not.cut2 )  error stop 'getCellAreas: error intersection with horizontal line'

            xisct = x1 + ( x2 - x1 ) * rlam1    ! x-coordinate of intersection point
            ix = ceiling( ( xisct - ox ) / wx )
!                                               ! Add point to cell (ix,iz)
            if( all( [ix,iz] == [ix1,iz1] ) .and. ip == 1 ) then
               isct = isct+1
               xtmp = xisct                     ! Intersection of first polygon segment
               ztmp = z1                        !    with closest line should be added later
            else
               cell(ix,iz)%np = cell(ix,iz)%np + 1
               cell(ix,iz)%nisct = cell(ix,iz)%nisct + 1
               cell(ix,iz)%x(cell(ix,iz)%np) = xisct
               cell(ix,iz)%z(cell(ix,iz)%np) = z1
            end if
!                                               ! Add point to cell (ix,iz+1)
            if( all( [ix,iz+1] == [ix1,iz1] ) .and. ip == 1 ) then
               isct = isct+1
               xtmp = xisct
               ztmp = z1
            else
               cell(ix,iz+1)%np = cell(ix,iz+1)%np + 1
               cell(ix,iz+1)%nisct = cell(ix,iz+1)%nisct + 1
               cell(ix,iz+1)%x(cell(ix,iz+1)%np) = xisct
               cell(ix,iz+1)%z(cell(ix,iz+1)%np) = z1
            end if
         end do
!
!  Check intersection with vertical (x=const) lines
!
         ixp1 = ceiling( (px(ip) - ox) / wx )
         ixp2 = ceiling( (px(ip+1) - ox) / wx )
         dix = nint( sign( ONE, real(ixp2 - ixp1, REAL8) ) )

         if( dix > 0 ) then
            icx1 = ixp1
            icx2 = ixp2-1
         else
            icx1 = ixp1-1
            icx2 = ixp2
         end if

         do ix = icx1, icx2, dix
            z1 = oz
            z2 = oz + nz * wz
            x1 = ox + ix * wx

            call mcisct( rlam1, rlam2, cut1, cut2, x1, z1, x1, z2,  &
!                ======    o      o     o     o    i   i   i   i
                         px(ip), pz(ip), px(ip+1), pz(ip+1) )
!                          i       i         i         i

            if( .not.cut2 )  error stop 'getCellAreas: error intersection with vertical line'

            zisct = z1 + ( z2 - z1 ) * rlam1    ! z-coordinate of intersection point
            iz = ceiling( ( zisct - oz ) / wz )
                                                ! Add point to cell (ix,iz)
            if( all( [ix,iz] == [ix1,iz1] ) .and. ip == 1 ) then
               isct = isct+1
               xtmp = x1                        ! Intersection of first polygon segment
               ztmp = zisct                     !    with closest line should be added later
            else
               cell(ix,iz)%np = cell(ix,iz)%np + 1
               cell(ix,iz)%nisct = cell(ix,iz)%nisct + 1
               cell(ix,iz)%x(cell(ix,iz)%np) = x1
               cell(ix,iz)%z(cell(ix,iz)%np) = zisct
            end if
                                                ! Add point to cell (ix+1,iz)
            if( all( [ix+1,iz] == [ix1,iz1] ) .and. ip == 1 ) then
               isct = isct+1
               xtmp = x1
               ztmp = zisct
            else
               cell(ix+1,iz)%np = cell(ix+1,iz)%np + 1
               cell(ix+1,iz)%nisct = cell(ix+1,iz)%nisct + 1
               cell(ix+1,iz)%x(cell(ix+1,iz)%np) = x1
               cell(ix+1,iz)%z(cell(ix+1,iz)%np) = zisct
            end if
         end do
!
!  Score second point (ip+1) of line segment into the appropriate cell
!
         ix = ceiling( ( px(ip+1) - ox ) / wx )
         iz = ceiling( ( pz(ip+1) - oz ) / wz )
         cell(ix,iz)%np =  cell(ix,iz)%np + 1
         cell(ix,iz)%x(cell(ix,iz)%np) = px(ip+1)
         cell(ix,iz)%z(cell(ix,iz)%np) = pz(ip+1)

      end do   ! loop over points

      !  check whether the cell ix1,iz1 is intersected more than once
      if( isct > 1)  error stop 'cell(ix1,iz1) intersected more than once by line segment 1'
!
!  Add saved intersection point to the ic1, iz1 cell. Must not happened twice
!
      if( isct == 0 )  error stop 'first line does not intersect the first cell'
      cell(ix1,iz1)%np = cell(ix1,iz1)%np + 1
      cell(ix1,iz1)%nisct = cell(ix1,iz1)%nisct + 1
      cell(ix1,iz1)%x(cell(ix1,iz1)%np) = xtmp
      cell(ix1,iz1)%z(cell(ix1,iz1)%np) = ztmp
!
!  Construct the cell polygons and detect the fractional cell area inside the polygon
!
      do ix = 1, nx
         do iz = 1, nz
            if( modulo( cell(ix,iz)%nisct, 2 ) /= 0 )  &
                  error stop 'getCellAreas: error: odd number of intersection points'
            if( cell(ix,iz)%np == 0 ) then ! no polygon(inters.) points in the cell
               !  cell completely inside or outside the polygon
               xtmp = ox + (ix-0.5)*wx   ! enough to check whether the cell center
               ztmp = oz + (iz-0.5)*wz   ! is inside/outside the region

               if( mcloc( xtmp, ztmp, px, pz, coseta, sineta ) > 0 )  area(ix,iz) = ONE
!                  =====    i    i    i    i     i      i                ! inside
            else ! np > 0
               ncins = 0  ! number of cell corners inside the polygon
               do iz1 = -1, 0                   ! loop over all four corners
                  do ix1 = -1, 0                ! Detect how many cell corners are
                     xtmp = ox + (ix+ix1)*wx    ! inside the region (max=4)
                     ztmp = oz + (iz+iz1)*wz    !
                     if( mcloc( xtmp, ztmp, px, pz, coseta, sineta ) > 0 ) then
!                        =====    i    i    i    i     i      i
                        ncins = ncins + 1
                        cx(ncins) = xtmp
                        cz(ncins) = ztmp
                     end if
                  end do
               end do
!
!  Several cases depending on number of intersection points between cell and region
!
               lsub = .false.
               np1 = cell(ix,iz)%np  ! simplify the notation

               select case( cell(ix,iz)%nisct )

                  case( 2 )  ! copy the first portion of the array
                     cpx(1:np1) = cell(ix,iz)%x(1:np1)
                     cpz(1:np1) = cell(ix,iz)%z(1:np1)

                     select case( ncins )   ! number of included corner points
                        case( 0 )              ! just close the polygon
                           ncp = np1 + 1
                        case( 1 )              ! one corner included append it
                           ncp = np1 + 2       ! corner + closing point (first point)
                           cpx(np1+1) = cx(1)  ! first: append the corner point
                           cpz(np1+1) = cz(1)
                        case( 2 )
                           ncp = np1 + 3  ! 2 corners + first point
                           ! append closest corner point to the last scored point
                           xtmp = cpx(np1)
                           ztmp = cpz(np1)
                           if( (xtmp-cx(1))**2 + (ztmp-cz(1))**2 <  &
                               (xtmp-cx(2))**2 + (ztmp-cz(2))**2 ) then
                              cpx(np1+1) = cx(1)
                              cpz(np1+1) = cz(1)
                              cpx(np1+2) = cx(2)
                              cpz(np1+2) = cz(2)
                           else
                              cpx(np1+1) = cx(2)
                              cpz(np1+1) = cz(2)
                              cpx(np1+2) = cx(1)
                              cpz(np1+2) = cz(1)
                           end if
                        case( 3 ) ! three corners are inside
                           ! close the complementary polygon with the outside corner and
                           ! subtract the area
                           ncp = np1 + 2  ! corner + closing point (first point)
                           cpx(np1+1) = cx(1)   ! first: append the corner point
                           cpz(np1+1) = cz(1)
                           lsub = .true.        ! subtract area flag
                        case( 4 )  ! just close the polygon
                           lsub = .true.
                           ncp = np1 + 1
                        case default
                           error stop 'getCellAreas: error(1) nins > 4'
                     end select
                  case( 4 )
                     select case( ncins )
                        case( 0 )  ! just close the polygon
                           ncp = np1 + 1
                           cpx(1:np1) = cell(ix,iz)%x(1:np1)
                           cpz(1:np1) = cell(ix,iz)%z(1:np1)
                        case( 1, 3 )  ! insert corner point such that makes the 90deg
                           ncp = np1 + 2  ! corner + closing point (first point)
                           do icp = 2, np1-1
                              x1=cell(ix,iz)%x(icp)
                              z1=cell(ix,iz)%z(icp)
                              x2=cell(ix,iz)%x(icp+1)
                              z2=cell(ix,iz)%z(icp+1)
                              if( abs( (cx(1)-x1)*(x2-cx(1)) + (cz(1)-z1)*(z2-cz(1)) ) <  &
                                  EPSMB)  exit   ! insert the corner
                           end do
                           cpx(1:icp) = cell(ix,iz)%x(1:icp) !  Copy the points up to ip
                           cpz(1:icp) = cell(ix,iz)%z(1:icp)
                           cpx(icp+1) = cx(1)                  ! insert the corner
                           cpz(icp+1) = cz(1)
                           cpx(icp+2:np1+1) = cell(ix,iz)%x(icp+1:np1) ! insert the rest points
                           cpz(icp+2:np1+1) = cell(ix,iz)%z(icp+1:np1)
                           if( ncins == 3 )  lsub = .true.  ! TAKE CARE ABOUT ORDER OF POINTS
                        case( 2 )   ! we have to insert two corners so they made 90deg with the
                           ncp = np1 + 3  ! pre- and succesive points
                           ! if the corner points are neighboring than append them to the end
                           if( abs( cx(2)-cx(1) ) + abs( cz(2)-cz(1) )-wz < EPSMB .or.  &
                               abs( cx(2)-cx(1) ) + abs( cz(2)-cz(1) )-wx < EPSMB ) then
                              cpx(1:np1) = cell(ix,iz)%x(1:np1)
                              cpz(1:np1) = cell(ix,iz)%z(1:np1)
                              xtmp = cell(ix,iz)%x(np1)
                              ztmp = cell(ix,iz)%z(np1)
                              if( (xtmp-cx(1))**2 + (ztmp-cz(1))**2 <  &
                                  (xtmp-cx(2))**2 + (ztmp-cz(2))**2 ) then  ! corner 1 closer
                                 cpx(np1+1) = cx(1)
                                 cpz(np1+1) = cz(1)
                                 cpx(np1+2) = cx(2)
                                 cpz(np1+2) = cz(2)
                              else
                                 cpx(np1+1) = cx(2)
                                 cpz(np1+1) = cz(2)
                                 cpx(np1+2) = cx(1)
                                 cpz(np1+2) = cz(1)
                              end if
                           else
                              ip = 1   ! first try to insert first corner point
                              do icp = 2, np1-1
                                 x1 = cell(ix,iz)%x(icp)
                                 z1 = cell(ix,iz)%z(icp)
                                 x2 = cell(ix,iz)%x(icp+1)
                                 z2 = cell(ix,iz)%z(icp+1)
                                 if( ip > 2 )  error stop 'getCellAreas: internal error ' //  &
                                                    'nisct=4, nins=2'
                                 if( abs( (cx(ip)-x1)*(x2-cx(ip)) + (cz(ip)-z1)*(z2-cz(ip)) )  &
                                     < EPSMB ) then
                                    if( ip == 1 ) then
                                       ip1 = icp
                                    else
                                       ip2 = icp
                                    end if
                                    exit     ! insert the corner
                                    ip = ip + 1 ! second corner point
                                 end if
                              end do
                              cpx(1:ip1) = cell(ix,iz)%x(1:ip1)   ! Copy the points up to first
                              cpz(1:ip1) = cell(ix,iz)%z(1:ip1)   ! insertion point
                              cpx(ip1+1) = cx(1)                  ! insert the first corner
                              cpz(ip1+1) = cz(1)
                              cpx(ip1+2:ip2+1) = cell(ix,iz)%x(ip1+1:ip2) ! copy the rest portion
                              cpz(ip1+2:ip2+1) = cell(ix,iz)%z(ip1+1:ip2) ! before the second corner
                              cpx(ip2+2) = cx(2)                   ! insert the second corner
                              cpz(ip2+2) = cz(2)
                              cpx(ip2+3:np1+2) = cell(ix,iz)%x(ip2+1:np1) ! insert the rest points
                              cpz(ip2+3:np1+2) = cell(ix,iz)%z(ip2+1:np1)
                           end if
                        case default
                           error stop 'getCellAreas: error(2) cell corner points'
                     end select
                  case default
                     error stop 'getCellAreas: cannot handle the case when a cell is ' //  &
                          'intersected more than 4 times by a region polygon'
               end select
               ! close the polygon with the first point
               cpx(ncp) = cell(ix,iz)%x(1)
               cpz(ncp) = cell(ix,iz)%z(1)
               !  Cell polygon area = fraction of cell inside the polygon
               area(ix,iz) = getPolyArea( cpx(1:ncp), cpz(1:ncp), ncp )
!                            ===========       i           i       i
               area(ix,iz) = abs(area(ix,iz)) / ascale
               if( lsub )  area(ix,iz) = ONE - area(ix,iz)
            end if
         end do
      end do

      deallocate( cell )

   end subroutine getCellAreas_old
!
!=======================================================================
!
   subroutine relPosPolygon( x, z, px, pz, inside, onedge )
!             =============  i  i  i   i     o       o
!
!  Determine whether point (x,z) is inside the polygon with vertices (px,pz). In addition,
!  determine whether point (x,z) is on a polygon edge.
!
!  THIS VERSION IS CURRENTLY RESTRICTED TO QUADRILATERALS.
!
!  Arguments
!
      real(REAL8), intent(in) :: x        ! x-coordinate of point to be checked
      real(REAL8), intent(in) :: z        ! z-coordinate of point to be checked
      real(REAL8), intent(in) :: px(4)    ! x-coordinates of quadrilateral vertices
      real(REAL8), intent(in) :: pz(4)    ! z-coordinates of quadrilateral vertices
      logical, intent(out) :: inside      ! flag indicating point is inside the polygon (not
                                          !    on an edge and not outside)
      logical, intent(out) :: onedge(4)   ! flags indicating whether point is on each of the
                                          !    edges
!
!  Local Variables
!
      real(REAL8) :: area(4)
      logical :: notoutside
      integer(INT4) :: i
!
!  Start of Executable Code
!  ------------------------
!
!  Determine (twice the) signed area of triangle formed by the point and each of the polygon
!  edges. If the area is positive, the point is left of the edge seen in the direction from the
!  first to the second polygon vertex (or the triangle x->p(i)->p(i+1) is counterclockwise)
!
      do i = 1, 4
         area(i) = (px(i)-x)*(pz(next( i ))-pz(i)) - (pz(i)-z)*(px(next( i ))-px(i))
!                                ====  i                           ====  i
      end do
!
!  Look for a reentrant corner (in a quadrilateral there can be at most one)
!
      do i = 1, 4
         if( (px(prev(i))-px(i))*(pz(next( i ))-pz(prev( i ))) -  &
!                                    ====  i       ====  i
             (pz(prev(i))-pz(i))*(px(next( i ))-px(prev( i ))) > 0 )  exit
!                                    ====  i       ====  i
      end do
!
!  Check whether point is inside the polygon
!
      if( i == 5 ) then                            ! no reentrant corner
         inside = all( area(:) > 0 )
      else                                         ! reentrant corner at vertex i
         inside = (area(prev( i )) > 0 .or. area(i) > 0) .and.  &
!                       ====  i
                  area(next( i )) > 0 .and. area(next( next( i ) )) > 0
!                      ====  i                   ====      i
      end if
!
!  Check whether point is on one of the edges
!
      if( inside ) then
         onedge(:) = .false.
      else
         if( i == 5 ) then
            notoutside = all( area(:) >= 0 )
         else
            notoutside = (area(prev( i )) >= 0 .or. area(i) >= 0) .and.  &
!                              ====  i
                         area(next( i )) >= 0 .and. area(next( next( i ) )) >= 0
!                             ====  i                    ====      i
         end if
         onedge(:) = ( area(:) == 0 ) .and. notoutside
      end if

   contains
!
!-----------------------------------------------------------------------
!
      integer(INT4) function next( i )
!                            ====  i
!
!  Get the next index (cyclic)
!
!  Replaces the statement function next(i) = modulo(i,4) + 1
!
         integer(INT4), intent(in) :: i

         next = modulo(i,4) + 1

      end function next
!
!-----------------------------------------------------------------------
!
      integer(INT4) function prev( i )
!                            ====  i
!
!  Get the previous index (cyclic)
!
!  Replaces the statement function prev(i) = modulo(i-2,4) + 1
!
         integer(INT4), intent(in) :: i

         prev = modulo(i-2,4) + 1

      end function prev
!
!-----------------------------------------------------------------------
!
   end subroutine relPosPolygon
!
!=======================================================================
!
   subroutine pointInsidePolygon( px, pz, x, z )
!             ==================  i   i   o  o
!
!  Return a point (x,z) that is safely inside the polygon with vertices (px,pz).
!  The point is at the center of one of the diameters.
!
!  THIS VERSION IS CURRENTLY RESTRICTED TO QUADRILATERALS.
!
!  Arguments
!
      real(REAL8), intent(in) :: px(4)    ! x-coordinates of quadrilateral vertices
      real(REAL8), intent(in) :: pz(4)    ! z-coordinates of quadrilateral vertices
      real(REAL8), intent(out) :: x       ! x-coordinate of point to be checked
      real(REAL8), intent(out) :: z       ! z-coordinate of point to be checked
!
!  Local Variables
!
      real(REAL8) :: a, amax
      integer(INT4) :: i, imax
!
!  Start of Executable Code
!  ------------------------
!
!  Determine the "most-reentrant" (or "least non-reentrant") corner
!
      amax = -BIG
      do i = 1, 4
         a = (px(prev(i))-px(i))*(pz(next(i))-pz(prev(i))) -  &
             (pz(prev(i))-pz(i))*(px(next(i))-px(prev(i)))
         if( a > amax ) then
            amax = a
            imax = i
         end if
      end do
!
!  Determine the center of the diagonal of corner i
!
      x = HALF * ( px(imax) + px(next(next(imax))) )
      z = HALF * ( pz(imax) + pz(next(next(imax))) )

   contains
!
!-----------------------------------------------------------------------
!
      integer(INT4) function next( i )
!                            ====  i
!
!  Get the next index (cyclic)
!
!  Replaces the statement function next(i) = modulo(i,4) + 1
!
         integer(INT4), intent(in) :: i

         next = modulo(i,4) + 1

      end function next
!
!-----------------------------------------------------------------------
!
      integer(INT4) function prev( i )
!                            ====  i
!
!  Get the previous index (cyclic)
!
!  Replaces the statement function prev(i) = modulo(i-2,4) + 1
!
         integer(INT4), intent(in) :: i

         prev = modulo(i-2,4) + 1

      end function prev
!
!-----------------------------------------------------------------------
!
   end subroutine pointInsidePolygon
!
!=======================================================================
!
   subroutine getDistPolyPoly( px1, pz1, px2, pz2, distmin2min, distmin2max )
!             ===============   i    i    i    i        o            o
!
!  Determine the minimum distance (distmin2min) between any two points of two polygons, and the
!  maximum over all points of the first polygon of the minimum distance from any point of the
!  second polygon (distmin2max). By polygon we mean the area enclosed by the polygon.
!
!  THIS VERSION IS CURRENTLY RESTRICTED TO QUADRILATERALS.
!
!  Arguments
!
      real(REAL8), intent(in) :: px1(4)         ! x-coordinate of point to be checked
      real(REAL8), intent(in) :: pz1(4)         ! z-coordinate of point to be checked
      real(REAL8), intent(in) :: px2(4)         ! x-coordinates of quadrilateral vertices
      real(REAL8), intent(in) :: pz2(4)         ! z-coordinates of quadrilateral vertices
      real(REAL8), intent(out) :: distmin2min   ! square of minimum distance
      real(REAL8), intent(out) :: distmin2max   ! square of maximum minimum distance of
                                                !    polygon 1
!
!  Local Variables
!
      real(REAL8) :: dist2, distmin2
      logical :: inside, onedge(4)
      integer(INT4) :: ip1, ip2
!
!  Start of Executable Code
!  ------------------------
!
!  Consider distances of polygon 1 corners from polygon 2
!
      distmin2min = BIG
      distmin2max = ZERO
      do ip1 = 1, 4
         call relPosPolygon( px1(ip1), pz1(ip1), px2(:), pz2(:), inside, onedge(:) )
!             =============     i         i        i       i       o        o
         if( inside .or. any( onedge(:) ) ) then
            distmin2 = ZERO
         else
            distmin2 = BIG
            do ip2 = 1, 4
               call mcdist( dist2, px1(ip1), pz1(ip1),  &
!                   ======    o       i         i
                            px2(ip2), pz2(ip2), px2(next(ip2)), pz2(next(ip2)) )
!                              i         i          i                i
               distmin2 = min( distmin2, dist2 )
            end do
         end if
         distmin2min = min( distmin2min, distmin2 )
         distmin2max = max( distmin2max, distmin2 )
      end do
!
!  Consider distances of polygon 2 corners from polygon 1
!  (note: distmin2max is always from a polygon 1 corner, therefore there are no more updates
!  for distmin2max)
!
      do ip2 = 1, 4
         call relPosPolygon( px2(ip2), pz2(ip2), px1(:), pz1(:), inside, onedge(:) )
!             =============    i        i        i       i         o        o
         if( inside .or. any( onedge(:) ) ) then
            distmin2min = ZERO
            exit
         end if
         do ip1 = 1, 4
            call mcdist( distmin2, px2(ip2), pz2(ip2),  &
!                ======     o         i         i
                         px1(ip1), pz1(ip1), px1(next( ip1 )), pz1(next( ip1 )) )
!                           i         i             i                i
            distmin2min = min( distmin2min, distmin2 )
         end do
      end do

   contains
!
!-----------------------------------------------------------------------
!
      integer(INT4) function next( i )
!                            ====  i
!
!  Get the next index (cyclic)
!
!  Replaces the statement function next(i) = modulo(i,4) + 1
!
         integer(INT4), intent(in) :: i

         next = modulo(i,4) + 1

      end function next
!
!-----------------------------------------------------------------------
!
   end subroutine getDistPolyPoly
!
!=======================================================================
!
   subroutine writePolygons( ipend, px, pz, dose )
!             =============    i    i   i    i
!
!     Write the contour as a set of closed surfaces (connected components) to the SURF file.
!
!  Arguments:
!
      integer(INT4), intent(in) :: ipend(:)     ! last point index of connected component
      real(REAL8), intent(in) :: px(:)          ! x-coordinates of polygon corners
      real(REAL8), intent(in) :: pz(:)          ! z-coordinates of polygon corners
      real(REAL8), intent(in) :: dose           ! dose per ion
!
!  Local Variables:
!
      integer(INT4) :: icc, ncc, ip, ipbeg(size( ipend ))
!
      character(6), save :: fmode = 'write'
!
!  Start of Executable Code:
!  ------------------------
!
      call openFile( 22, 'SURF', 'formatted', fmode )
!          ========  i     i          i         i
      fmode = 'append'
!
!  Write dose to first line
!
      write( 22, '(a,1p,e15.6)' ) '# Dose:', dose
!
!  Number of connected components
!
      ncc = size( ipend )
      write( 22, * )  ncc
!
!  Begin index of each connected component
!
      ipbeg(1) = 1
      ipbeg(2:ncc) = ipend(1:ncc-1)+1
!
!  For each connected component write number of points and points to the file
!
      do icc = 1, ncc

         write( 22, * )  ipend(icc)-ipbeg(icc)+1
         do ip = ipbeg(icc), ipend(icc)
            write( 22, '(1p,2e15.6)' )  px(ip), pz(ip)
         end do

      end do

      close( 22 )

   end subroutine writePolygons
!
!=======================================================================
!
end module Geom_Polygon_Old_mod
