module Geom_Polygon_mod
!
!  Provide a collection of procedures operating on a closed polygon
!
!  Used Modules:
!  ============
!
   use Const_mod,             only: INT4, REAL8,  &
                                    ZERO, ONE, TWO, HALF, QUART, PI, BIG, EPSMA, EPSMB
   use Array_mod,             only: alloc

   use Geom_Line_mod,         only: mccut, mcdist, mcisct
   use Geom_Polygon_old_mod,  only: mcauxl, mcloc, mcdis4, getDistPolygon
   use Geom_FFPcyl_mod,       only: isInsideFFPcyl, intersectFFPcyl
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
   public :: initPolygon            ! initialize derived properties of a polygon
!            ===========
   public :: isInsidePolygon        ! determine whether point is inside polygon
!            ===============
   public :: distPolygon            ! determine distance between point and polygon
!            ===========
   public :: distRectPolygon        ! determine distance between axes-parallel rectangle and
!            ===============        !  polygon
   public :: ovlpFFPcylPolygon      ! determine whether ffp cylinder overlaps polygon
!            =================
!
!  Private Procedures:
!  ==================
!
   private :: mcreentr              ! determine reentrant corners of a polygon
!             ========
!
!  Public Data Type:
!  ================
!
   type, public :: Polygon_type
      integer(INT4) :: np                 ! number of vertices
      real(REAL8), allocatable :: px(:)   ! (np+1) x-values of vertices (last=first)
      real(REAL8), allocatable :: pz(:)   ! (np+1) y-values of vertices (last=first)
      real(REAL8) :: coseta               ! cosine of direction of auxiliary line
      real(REAL8) :: sineta               ! sine of direction of auxiliary line
      logical, allocatable :: reentr(:)   ! (np+1) flag indicating reentrant corners
      logical, allocatable :: virtual(:)  ! (np) flag indicating virtual segments (if .true.,
                                          !  segment does not exist)
   end type Polygon_type

contains
!
!=======================================================================
!
   subroutine initPolygon( poly )
!             ===========   io
!
!     Initialize derived properties of a polygon. px(:) and pz(:) must already exist.
!
!  Arguments:
!
      type(Polygon_type), intent(inout) :: poly ! polygon to be initialized
!
!  Local Variables:
!
!     None
!
!  Start of Executable Code:
!  ------------------------
!
!  Calculate direction of auxiliary line
!
      call mcauxl( poly%coseta, poly%sineta, poly%px(:), poly%pz(:) )
!          ======       o            o           i           i
!
!  Determine reentrant corners
!
      if( allocated( poly%reentr ) )  deallocate( poly%reentr )
      call alloc( poly%reentr, poly%np+1, 'poly%reentr', 'initPolygon' )
!          =====       o           i            i              i
      call mcreentr( poly%reentr(:), poly%px, poly%pz, poly%coseta, poly%sineta )
!          ========       o             i        i          i            i
!
!  Initialize virtual segment flag
!
      if( allocated( poly%virtual ) )  deallocate( poly%virtual )
      call alloc( poly%virtual, poly%np, 'poly%virtual', 'initPolygon' )
!          =====       o          i            i               i
      poly%virtual(:) = .false.

   end subroutine initPolygon
!
!=======================================================================
!
   subroutine mcreentr( reentr, px, pz, coseta, sineta )
!             ========    o     i   i     i       i
!
!     Determines the reentrant corners of the polygon
!
!  Arguments:
!
      logical :: reentr(:)                ! flag indicating reentrant corner
      real(REAL8), intent(in) :: px(:)    ! x-coordinates of polygon corners (first=last)
      real(REAL8), intent(in) :: pz(:)    ! z-coordinates of polygon corners (first=last)
      real(REAL8), intent(in) :: coseta   ! cosine of direction of auxiliary line
      real(REAL8), intent(in) :: sineta   ! sine of direction of auxiliary line
!
!  Local Variables:
!
      integer(INT4) :: np, jp, iloc
      real(REAL8) :: rx(size( reentr )-1), rz(size( reentr )-1)
!
!  Start of Executable Code:
!  ------------------------
!
      np = size( reentr ) - 1

      do jp = 1, np+1

         rx(1:jp-1) = px(1:jp-1)
         rz(1:jp-1) = pz(1:jp-1)
         rx(jp:np) = px(jp+1:np+1)
         rz(jp:np) = pz(jp+1:np+1)

         if( jp == 1 ) then
            rx(np) = rx(1)
            rz(np) = rz(1)
         else if( jp == np+1 ) then
            rx(1) = rx(np)
            rz(1) = rz(np)
         end if

         iloc = mcloc( px(jp), pz(jp), rx(:), rz(:), coseta, sineta )
!               =====    i       i       i      i      i       i

         reentr(jp) = (iloc == 1)

      end do

   end subroutine mcreentr
!
!=======================================================================
!
   integer(INT4) function isInsidePolygon( pos, poly )  result( iloc )
!                         ===============   i    i               o
!
!     Determine whether point pos=(x,z) is inside the polygon poly
!
!        iloc = +1 if the point is inside the polygon
!                0 if the point is on the polygon
!               -1 if the point is outside or on the polygon
!
!  Arguments:
!
      real(REAL8), intent(in) :: pos(2)         ! point to be investigated
      type(Polygon_type), intent(in) :: poly    ! polygon
!
!  Local Variables:
!

!
!  Start of Executable Code:
!  ------------------------
!
      iloc = mcloc( pos(1), pos(2), poly%px(:), poly%pz(:), poly%coseta, poly%sineta )
!            =====   i       i          i           i            i            i

   end function isInsidePolygon
!
!=======================================================================
!
   subroutine distPolygon( distsq, pos, poly )
!             ===========    o      i    i
!
!     Get the distance of position pos from the polygon poly
!
!  Arguments:
!
      real(REAL8), intent(out) :: distsq        ! square of the distance of position pos from
                                                !    polygon poly
      real(REAL8), intent(in) :: pos(2)         ! x- and z-coordinate of point
      type(Polygon_type), intent(in) :: poly    ! polygon
!
!  Start of Executable Code:
!  ------------------------
!
      call getDistPolygon( distsq, pos(1), pos(2), poly%px(:), poly%pz(:), poly%virtual(:) )
!          ==============    o       i       i        i           i               i

   end subroutine distPolygon
!
!==========================================================================================
!
   subroutine distRectPolygon( distminsq, distmaxsq, inside,  &
!             ===============      o          o        o
                               x1, x2, z1, z2, poly )
!                              i   i   i   i    i
!
!     Calculates the square of the minimum and maximum distance of points inside an
!     axis-parallel rectangle from a polygon.
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
                                                !  rectangle from the polygon
      real(REAL8), intent(out) :: distmaxsq     ! maximum distance of any point inside the
                                                !  rectangle from the polygon
                                                !  between distminsq and distmaxsq
      logical, intent(out) :: inside            ! flag indicating center of rectangle is inside
                                                !  the polygon
      real(REAL8), intent(in) :: x1             ! left boundary of rectangle
      real(REAL8), intent(in) :: x2             ! right boundary of rectangle
      real(REAL8), intent(in) :: z1             ! lower boundary of rectangle
      real(REAL8), intent(in) :: z2             ! upper boundary of rectangle
      type(Polygon_type), intent(in) :: poly    ! polygon
!
!  Local variables
!
      real(REAL8) :: rsq                        ! radius of circle circumscribing rectangle
!
!  Start of Executable Code:
!  ------------------------
!
      rsq = QUART * ((x2-x1)**2 + (z2-z1)**2)

      call mcdis4( distminsq, distmaxsq, inside,  &
!          ======      o          o        o
                   x1, x2, z1, z2, rsq, poly%px, poly%pz, poly%reentr, poly%coseta, poly%sineta)
!                  i   i   i   i    i      i        i          i            i            i

   end subroutine distRectPolygon
!
!=======================================================================
!
   logical function ovlpFFPcylPolygon( poly, pos0, dir, ffp, pmax, b_a, chkInside_ )  &
!                   =================   i     i     i    i    i     i       i
                    result( ovlp )
!                            o
!
!     Determine whether the free flight path cylinder overlaps the polygon area.
!     If the ffp cylinder is completely inside the polygon, this is counted as an overlap
!     only if chkInside_ is set.
!
!  Arguments:
!
      type(Polygon_type), intent(in) :: poly       ! polygon
      real(REAL8), intent(in) :: pos0(2)           ! starting point of ffp
      real(REAL8), intent(in) :: dir(2)            ! direction of ffp (dir(1)**2+dir(2)**2=1)
      real(REAL8), intent(in) :: ffp               ! projection of ffp to xz-plane = of height
                                                   !    of ffp cylinder
      real(REAL8), intent(in) :: pmax              ! maximum impact parameter = radius of ffp
                                                   !    cylinder
      real(REAL8), intent(in) :: b_a               ! ratio of axis lengths of ellipse (<=1)
      logical, optional, intent(in) :: chkInside_  ! check if ffp cylinder is completely inside
!                                                  ! polygon
!  Local Variables:
!
      logical :: chkInside
      integer(INT4) :: ip, np
!
!  Start of Executable Code:
!  ------------------------
!
      if( present( chkInside_ ) ) then
         chkInside = chkInside_
      else
         chkInside = .false.
      end if
!
!  Check whether the the ffpcyl and polygon boundaries intersect
!
      ovlp = .false.
      np = poly%np
      do ip = 1, np
         ovlp = intersectFFPcyl( [poly%px(ip),poly%pz(ip)],  &
!               ===============               i
                                 [poly%px(ip+1),poly%pz(ip+1)],  &
!                                               i
                                 pos0, dir, ffp, pmax, b_a )
!                                 i     i    i    i     i
         if( ovlp )  exit
      end do
!
!  Check whether the whole polygon is inside the ffp cylinder. If there is not an intersection
!  point, it is enough to check one point
!
      if( .not.ovlp ) then
         ovlp = isInsideFFPcyl( [poly%px(1),poly%pz(1)], pos0, dir, ffp, pmax, b_a )
!               ==============              i             i     i    i    i     i
      end if
!
!  Check whether the whole ffp cylinder is inside the polygon. If there is not an intersection
!  point, it is enough to check one point
!
      if( chkInside ) then
         if( .not.ovlp ) then
            ovlp = ( isInsidePolygon( pos0, poly ) > 0 )
!                    ===============   i     i
         end if
      end if

end function ovlpFFPcylPolygon
!
!=======================================================================
!
end module Geom_Polygon_mod
