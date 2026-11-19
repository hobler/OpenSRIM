module Geom_CartGrid2d_mod
!
!  Provides operations on a 2D Cartesian grid
!
!  Used Modules:
!  ============
!
   use Const_mod,             only: INT4, REAL8,  &
                                    ONE, HALF, ZERO, EPSMB, EPSMA
   use Array_mod,             only: check_allocation
   use Stack_mod,             only: Int4Stack, Vec2Stack,  &
                                    newStack, sizeOfStack, pushStack, popStack, stack2Array
   use Math_Sort_mod,         only: quickSort

   use Geom_Polygon_Old_mod,  only: mcauxl, mcloc, getPolyArea
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
   public :: findZeroLevelSetCartGrid2d   ! detects the zero level contour of input 2D array
!            ==========================
   public :: getOvlpAreasCartGrid2d       ! obtain the cell areas inside the polygon
!            ======================
   public :: bndConnected                 ! set elements of mask to false which are not
!            ============                 !    connected to the boundary of the array via other
!                                         !    false elements (UNUSED)
!  Private Procedures:
!  ==================
!
   private :: secondPoint                    ! interpolates the second point of the level
!             ===========                    !    set line of a cell
   private :: scorePoint                     ! as secondPoint provided first point is
!             ==========                     !    between corners 1 and 3
   private :: rotatePoint                    ! rotates the point placed at the cell edges
!             ===========
   private :: rotateCor                      ! rotates the four corners of the cell.
!             =========

contains
!
!==========================================================================================
!
   subroutine findZeroLevelSetCartGrid2d( val, ipend, px, pz )
!             ==========================   i     o    o   o
!
!     Returns the zero isovalue line (val=zero). The line may consist of disconnected closed
!     polygons (= connected components).
!
!     We number the cell corners like this:
!
!                             3|       |4
!                       iz+1 --+-------+--      z
!                              |       |       ^
!        edge(ix,iz,2) ------->|       |       |
!       ... "vertical"         |       |       |
!       line segment      iz --+-------+--     +---> x
!                             1|   ^   |2
!                             ix   | ix+1
!                                  |
!                          edge(ix,iz,1) ... "horizontal" line segment
!
!     The x-coordinate of vertical line ix is    x=ix-1
!     The z-coordinate of horizontal line iz is  z=iz-1
!
!     NOTE: Do not change the order of the indices since it is critical for the rest
!           of the code.
!
!  Arguments:
!
      real(REAL8), intent(in) :: val(:,:)                   ! cell values
      integer(INT4), intent(out), allocatable :: ipend(:)   ! last index of connected component
      real(REAL8), intent(out), allocatable :: px(:), pz(:) ! surface points
!
!  Local Variables:
!
      character(256) :: msg
      integer(INT4) :: ix1, iz1, ior, ix, iz, istat, case, invcase, ip, np
      real(REAL8) :: f(4), r(2), p(2), p1(2)
      real(REAL8), allocatable :: p2(:,:)
      logical :: newCont
      type(Vec2Stack) :: pStack
      type(Int4Stack) :: ipendStack
      logical :: checked(size( val, 1 ),size( val, 2 ),2)   ! checked grid lines array
!
!  Start of Executable Code:
!  ------------------------
!
      newCont = .false.

      checked(:,:,:) = .false.                  ! unprocessed grid lines at the start

      call newStack( pStack )                   ! stack of contour points
!          ========    io
      call newStack( ipendStack )               ! stack of end point indices of connected
!          ========      io                     !   components
!
!  Loop over all line segments to find the first point of a new contour
!
      do ix1 = 1, size( val, 1 )-1
         do iz1 = 1, size( val, 2 )-1
            do ior = 1, 2
               newCont = .false.

               f(1:4) = reshape( val(ix1:ix1+1,iz1:iz1+1), [4] )  ! neighboring values

               if( checked(ix1,iz1,ior) )  cycle ! line segment already checked

               checked(ix1,iz1,ior) = .true.     ! label as checked
!
!  Check for vertical/horizontal sign change
!
               if( ior == 2 ) then
                  if( f(1) * f(3) < ZERO ) then
                     ! interpolate the first point at the vertical line
                     r(1) = ZERO
                     r(2) = max( EPSMB, min( ONE-EPSMB, f(1) / ( f(1) - f(3) ) ) )
                     newCont = .true.
                     case = 0       ! 0°
                     invcase = 0    ! 0°
                  end if
               else
                  if( f(1) * f(2) < ZERO ) then
                     ! interpolate the first point at the horizontal line
                     r(1) = max( EPSMB, min( ONE-EPSMB, f(1) / ( f(1) - f(2) ) ) )
                     r(2) = ZERO
                     newCont = .true.
                     case = 4       ! -90°
                     invcase = 1    ! +90°
                  end if
               end if

               if( .not. newCont )  cycle
!
!  We have a new contour. Insert the first point into the list.
!
               p1(:) = [ix1-1,iz1-1] + r(:)
               call pushStack( p1, pStack )
!                   =========   i     io

               ix = ix1   ! separate the indices from those used in contour loop
               iz = iz1

               do                               ! loop over contour points

                  call secondPoint( f, r, case, invcase )  ! find second point in cell
!                      ===========  i  o    i      i

                  p(:) = [ix-1,iz-1] + r(:)
                  call pushStack( p, pStack ) ! insert the point
!                      =========  i    io
                                                ! point == first point ?
                  if( sum( (p(:)-p1(:))**2 ) < EPSMA )  exit
!
                  call getNextGrid( r(:), ix, iz, case, invcase ) ! get the next grid
!                      ===========   i    io  io    o      o

                  f(1:4) = reshape( val(ix:ix+1,iz:iz+1), [4] ) ! neighboring values

               end do
!
               call pushStack( sizeOfStack( pStack ), ipendStack )
!                   =========  ===========    i           io

            end do
         end do
      end do
!
!  Allocate the arrays to store the contour points
!
      call stack2Array( ipendStack, ipend )
!          ===========      io        o
      call stack2Array( pStack, p2 )
!          ===========    io    o
      np = size( p2, 2 )
      allocate( px( np ), pz( np ), stat=istat, errmsg=msg )
      call check_allocation( 'findZeroLevelSetCartGrid2d', 'px, pz', istat, msg )
      do ip = 1, np
         px(ip) = p2(1,ip)
         pz(ip) = p2(2,ip)
      end do
!
!------------------------------------------------------------------------------------------
!
   contains
!
!------------------------------------------------------------------------------------------
!
      subroutine getNextGrid( r, ix, iz, casef, invcase )
!                ===========  i  io  io    o       o
!
!     Returns the indices of the neighboring grid cell that contains the contour point.
!     Detects the rotation to be applied to a cell in order to conform the reference case.
!     REFERENCE CASE: the first point (scored in the previous step) is between f(1)&f(3)
!
!  Arguments:
!
         real(REAL8), intent(in) :: r(2)        ! second point in cell
         integer(INT4), intent(inout) :: ix, iz ! actual grid indices
         integer(INT4), intent(out) :: casef    ! case defining the rotation angle
         integer(INT4), intent(out) :: invcase  ! case(angle) defining inverse rotation
!
!  Local Variables:
!
         integer(INT4) :: ior, igx, igz          ! checked line indices
!
!  Start of Executable Code:
!  ------------------------
!
         if( abs(r(1)) < EPSMA .and. r(2) > EPSMA ) then ! r on the left side of grid(ix,iz)
            igx = ix
            igz = iz
            ix = ix - 1
            iz = iz
            casef = 5      ! rotate -180�
            invcase = 2    ! next grid is left to previous surface point between f2,f4
            ior = 2        ! point is on the vertical line
         elseif( r(1) > EPSMA .and. abs(r(2)) < EPSMA ) then ! right side
            igx = ix
            igz = iz
            ix = ix
            iz = iz-1
            casef = 6      ! rotate -270�
            invcase = 3
            ior = 1
         elseif( abs( r(1)-ONE ) < EPSMA .and. r(2) > ZERO ) then
            igx = ix+1
            igz = iz
            ix = ix+1
            iz = iz
            casef = 0      ! reference case deg=0�
            invcase = 0    ! next grid is left to previous surface point between f1,f3
            ior = 2        ! point is on the horizontal line
         elseif( r(1) > EPSMA .and. abs( r(2)-ONE ) < EPSMA ) then
            igx = ix
            igz = iz + 1
            ix = ix
            iz = iz+1
            casef = 4      ! rotate -90�
            invcase = 1    ! next grid is left to previous surface point between f1,f2
            ior = 1
         else
            error stop "getNextGrid: Error. Point inconsistent with grid indices!"
         end if
!
!  Set the correspoinding line segment to checked
!
         checked(igx,igz,ior) = .true.

      end subroutine getNextGrid
!
!------------------------------------------------------------------------------------------
!
   end subroutine findZeroLevelSetCartGrid2d
!
!==========================================================================================
!
   subroutine secondPoint( f, r, case, invcase )
!             ===========  i  o   i       i
!
!     Assuming that the first point is between f1&f3 corner interpolates the second point.
!
!  Arguments:
!
      real(REAL8), intent(in) :: f(4)           ! fill factors of edge cells
      real(REAL8), intent(out) :: r(2)          ! interpolated point
      integer(INT4), intent(in) :: case         ! case defining the rotation angle
      integer(INT4), intent(in) :: invcase      ! case(angle) defining inverse rotation
!
!  Local Variables:
!
      real(REAL8) :: floc(4)  ! auxiliary array
!
!  Start of Executable Code:
!  ------------------------
!
      floc = f

      if( case /= 0 )  call rotateCor( floc, case )   ! rotate the cell corners
!                           =========   io   i
      call scorePoint( floc, r )
!          ==========    i   o
      if( invcase /= 0 )  call rotatePoint( r, invcase ) ! rotate back the second point
!                              ===========  o     o

   end subroutine secondPoint
!
!==========================================================================================
!
   subroutine scorePoint( f, r )
!             ==========  i  o
!
!     Assuming that the first point in the cell is between f1&f3 corner detects
!     and interpolates the second point.
!
!  Arguments:
!
      real(REAL8), intent(in) :: f(4)           ! fill factors of edge cells
      real(REAL8), intent(out) :: r(2)          ! interpolated point
!
!  Local Variables:
!
!
!  Start of Executable Code:
!  ------------------------
!
      if( f(3) * f(2) > ZERO ) then             ! point between f1 & f2
         r(1) = max( EPSMB, min( ONE-EPSMB, f(1) / ( f(1) - f(2) ) ) )
         r(2) = ZERO
      elseif( f(3) * f(4) > ZERO ) then         ! point between f2 & f4
         r(1) = ONE
         r(2) = max( EPSMB, min( ONE-EPSMB, f(2) / ( f(2) - f(4) ) ) )
      else                                      ! point between f3 & f4
         r(1) = max( EPSMB, min( ONE-EPSMB, f(3) / ( f(3) - f(4) ) ) )
         r(2) = ONE
      end if

   end subroutine scorePoint
!
!==========================================================================================
!
   subroutine rotateCor( f, deg )
!             =========  io  i
!
!     Rotates the four corners of the cell.
!
!  Arguments:
!
      real(REAL8), intent(inout) :: f(4)        ! values to be rotated
      integer(INT4), intent(in) :: deg          ! flag denoting degrees
!
!  Local Variables:
!
      real(REAL8) :: tmp
!
!  Start of Executable Code:
!  ------------------------
!
      select case( deg )
         case( 1, 6 )      ! 90deg, -270deg
            tmp = f(2)
            f(2) = f(1)
            f(1) = f(3)
            f(3) = f(4)
            f(4) = tmp
         case( 2, 5 )      ! 180deg, -180deg
            tmp = f(1)
            f(1) = f(4)
            f(4) = tmp
            tmp = f(2)
            f(2) = f(3)
            f(3) = tmp
         case( 3, 4 )      ! 270deg, -90deg
            tmp = f(1)
            f(1) = f(2)
            f(2) = f(4)
            f(4) = f(3)
            f(3) = tmp
         case default
            error stop "rotate: error deg"
      end select

   end subroutine rotateCor
!
!===========================================================================================
!
   subroutine rotatePoint( r, deg )
!             ===========  io  i
!
!     Rotates the point placed at the cell edges.
!     Coordinates or point r(:) must be in unit coordinates: 0<r(:)<1
!
!  Arguments:
!
      real(REAL8), intent(inout) :: r(2)        ! point to be rotated
      integer(INT4), intent(in) :: deg          ! flag denoting degrees
!
!  Local Variables:
!
      real(REAL8) :: rtmp(2)
!
!  Start of Executable Code:
!  ------------------------
!
      select case( deg )
         case( 1, 6 )      ! 90deg, -270deg
            rtmp(1) = ONE - r(2)
            rtmp(2) = r(1)
            r(:) = rtmp(:)
         case( 2, 5 )      ! 180deg, -180deg
            rtmp(2) = ONE - r(2)
            rtmp(1) = ONE - r(1)
            r(:) = rtmp(:)
         case( 3, 4 )      ! 270deg, -90deg
            rtmp(1) = r(2)
            rtmp(2) = ONE - r(1)
            r(:) = rtmp(:)
         case default
            error stop "rotate: error deg"
      end select

   end subroutine rotatePoint
!
!===========================================================================================
!
   subroutine getOvlpAreasCartGrid2d( px, pz, nx, nz, area, ix_area, iz_area)
!             ======================  i   i   i   i    o       o        o
!
!     Obtain the fractional cell area inside the given polygon assuming a Cartesian grid,
!     i.e., the grid lines are assumed at integer positions:
!
!     positions of grid lines:   x= 0   1   2   3    nx-1  nx     (and similar for z)
!     cell index:               ix= | 1 | 2 | 3 | ... | nx |
!
!     The polygon may extend outside [0,nx] a,d [0,nz], but corresponding areas are ignored.
!
!     The first polygon corner must be repeated at the end of px, pz.
!
!     The subroutine can be used for rectilinear grids if polygon is scaled appropriately
!     before it is passed to this subroutine.
!
!  Arguments:
!
      real(REAL8), intent(in) :: px(:)                   ! x-coordinates of polygon corners
      real(REAL8), intent(in) :: pz(:)                   ! z-coordinates of polygon corners
      integer(INT4), intent(in) :: nx                    ! x-dimension of grid
      integer(INT4), intent(in) :: nz                    ! z-dimension of grid
      real(REAL8), allocatable, intent(out) :: area(:)   ! area inside polygon
                                                         !    first index runs in x-direction
                                                         !    second index runs in z-direction
      integer(INT4), allocatable, intent(out) :: ix_area(:) ! x-indices corresponding to the
                                                            !    areas
      integer(INT4), allocatable, intent(out) :: iz_area(:) ! z-indices corresponding to the
                                                            !    areas
!
!  Local Data Types:
!
      type :: Point_type
!             ==========
         real(REAL8) :: x                       ! x-coordinate
         real(REAL8) :: z                       ! z-coordinate
         integer(INT4) :: ix                    ! cell index in x-direction
         integer(INT4) :: iz                    ! cell index in z-direction
      end type

      type :: Line_type
!             =========
         type(Point_type) :: p1                 ! start point
         type(Point_type) :: p2                 ! end point
         real(REAL8) :: d                       ! length
         real(REAL8) :: dx                      ! x-component of length
         real(REAL8) :: dz                      ! z-component of length
      end type
!
!  Local Variables:
!
      logical :: lend
      character(256) :: errmsg
      integer(INT4) :: ixmin, ixmax, izmin, izmax, iarea, istat
      integer(INT4) :: ip, np, ix, iz, ii, ni, ia, na
      real(REAL8) :: polyarea, coseta, sineta, xc, zc

      type(Point_type) :: pp( size( px )-1 )    ! polygon points in scaled coordinates
      type(Point_type) :: pa( size( pp )+5 )    ! polygon for area calculation
      type(Point_type) :: pi( nx+nz )           ! intersection points (was: ...-2)

      logical :: lxi( size( pi ) )              ! lx*=.true.:  x-value is on grid line
      logical :: lxa( size( pa ) )              ! lx*=.false.: z-value is on grid line
      integer(INT4) :: isort( size( pp ) )

      type(Line_type) :: line
!
!  Start of Executable Code:
!  ------------------------
!
!  Convert px, pz into polygon
!
      np = size( px ) - 1

      do ip = 1, np
         pp(ip)%x = px(ip)
         pp(ip)%z = pz(ip)
      end do
!
!  Determine cell indices of polygon points and range of area indices affected
!
      ixmin = nx
      ixmax = 1
      izmin = nz
      izmax = 1
      do ip = 1, np
         pp(ip)%ix = ceiling( pp(ip)%x )
         pp(ip)%iz = ceiling( pp(ip)%z )
         ixmin = min( ixmin, pp(ip)%ix )
         ixmax = max( ixmax, pp(ip)%ix )
         izmin = min( izmin, pp(ip)%iz )
         izmax = max( izmax, pp(ip)%iz )
      end do
      ixmin = max( 1, ixmin )
      ixmax = min( nx, ixmax )
      izmin = max( 1, izmin )
      izmax = min( nz, izmax )

      allocate( area( (ixmax-ixmin+1)*(izmax-izmin+1) ), stat=istat, errmsg=errmsg )
      call check_allocation( 'getOvlpAreasCartGrid2d', 'area', istat, errmsg )
      area(:) = ZERO

      allocate( ix_area( (ixmax-ixmin+1)*(izmax-izmin+1) ), stat=istat, errmsg=errmsg )
      call check_allocation( 'getOvlpAreasCartGrid2d', 'ix_area', istat, errmsg )

      allocate( iz_area( (ixmax-ixmin+1)*(izmax-izmin+1) ), stat=istat, errmsg=errmsg )
      call check_allocation( 'getOvlpAreasCartGrid2d', 'iz_area', istat, errmsg )

      do iz = izmin, izmax
         do ix = ixmin, ixmax
            iarea = (ixmax-ixmin+1) * (iz-izmin) + (ix-ixmin+1)
            ix_area(iarea) = ix
            iz_area(iarea) = iz
         end do
      end do
!
!  Change orientation so area is positive
!
      polyarea = getPolyArea( pp(:)%x, pp(:)%z, np )
!                ===========     i        i     i

      if( polyarea < ZERO ) then
         isort(:) = [ (ip, ip=np,1,-1) ]      ! Note: pp has size np
         pp(:) = pp( isort(:) )
         polyarea = - polyarea                  ! Remember polyarea in case polygon is
      end if                                    !   completely inside a cell
!
!  Look for first intersection point of polygon with the grid
!
      do ip = 1, np

         call getLine( ip, line )
!             =======  i    o

         call getIntersectionPoints( line, ni, pi, lxi )
!             =====================   i    o   o    o

         if( ni /= 0 ) exit

      end do
!
!  If no intersection point has been found, polygon is completely inside a cell
!
      if( ni == 0 ) then

         ix = pp(1)%ix
         iz = pp(1)%iz
         if( ix >= 1 .and. ix <= nx .and. iz >= 1 .and. iz <= nz ) then
            iarea = (ixmax-ixmin+1) * (iz-izmin) + (ix-ixmin+1)
            area(iarea) = polyarea
         end if
         return

      end if
!
!  Initialize intersection point index and end flag
!
      ii = 1
      lend = .false.
!
!  Loop over cells
!  ---------------
!
      do
!
!  Start a new area polygon
!
         ia = 1                                 ! second old intersection point is
         pa(ia) = pi(ii)                        !    first point of new area polygon
         lxa(ia) = lxi(ii)
!
!  If this is the last intersection point of current polygon line, follow polygon until next
!  intersection point is found. All polygon corners passed belong to the area polygon
!
         if( ii == ni ) then

            do

               ia = ia + 1                      ! assign polygon corner to area polygon
               pa(ia) = line%p2

               if( ip == np ) then
                  ip = 1
                  lend = .true.                 ! having reached first line again, we have to
                                                ! run the loop body to the end and then exit
               else
                  ip = ip + 1
               end if

               call getLine( ip, line )         ! look for intersection points on the next line
!                   =======  i    o

               call getIntersectionPoints( line, ni, pi, lxi )
!                   =====================   i    o   o    o

               if( ni /= 0 )  exit              ! exit when intersection point has been found

            end do

            ii = 1                              ! reset intersection point index

         else

            ii = ii + 1                         ! iterate intersection point index

         end if
!
!  Assign intersection point to area polygon
!
         ia = ia + 1
         pa(ia) = pi(ii)
         lxa(ia) = lxi(ii)
!
!  Calculate cell index and polygon area
!
         na = ia
         call closePolygon( ix, iz, na, pa, lxa )
!             ============  o   o   io  io   i

         if( ix >= 1 .and. ix <= nx .and. iz >= 1 .and. iz <= nz ) then

            polyarea = getPolyArea( pa(:)%x, pa(:)%z, na )
!                      ===========     i        i     i      ! AREA SHOULD NOT BE NEGATIVE
!                                                            ! IF getPolyArea AND closePolygon
!                                                            ! are consistent
            if( polyarea < ZERO )  print *, 'ix,iz,polyarea=',ix,iz,polyarea
            iarea = (ixmax-ixmin+1) * (iz-izmin) + (ix-ixmin+1)
            area(iarea) = area(iarea) + polyarea

         end if

         if( lend )  exit

      end do
!
!  subtract integer such that area is in [0,1)  ! Note: this is essential and works
!
      do iarea = 1, size( area )
         area(iarea) = area(iarea) - floor(area(iarea))
      end do
!
!  set area of cells completely inside the polygon equal to 1
!  cell is completely inside polygon, if polygon does not intersect cell and cell center
!  is inside the polygon
!
      call mcauxl( coseta, sineta, px, pz )
!          ======    o       o     i   i

      do iarea = 1, size( area )
         if( area(iarea) /= ZERO )  cycle    ! exclude cells intersected by polygon
         xc = ix_area(iarea) - HALF
         zc = iz_area(iarea) - HALF
         if( mcloc( xc, zc, px, pz, coseta, sineta ) > 0 )  area(iarea) = ONE
!            =====  i   i   i   i     i       i
      enddo
!
! debug
!      do ix = 1, nx
!         print '(1000i3)',  (int(100*area(ix,iz)), iz=1,nz)
!      end do
! end-debug
!      error stop
!
   contains
!
!-----------------------------------------------------------------------------------------
!
      subroutine getLine( ip, line )
!                =======  i    o
!
!        Get the next line along the polygon pp
!
         integer(INT4), intent(in) :: ip
         type(Line_type), intent(out) :: line

         line%p1 = pp(ip)
         if( ip < np ) then
            line%p2 = pp(ip+1)
         else
            line%p2 = pp(1)
         end if

         line%dx = line%p2%x - line%p1%x
         line%dz = line%p2%z - line%p1%z
         line%d = sqrt( line%dx**2 + line%dz**2 )

      end subroutine getLine
!
!-----------------------------------------------------------------------------------------
!
      subroutine getIntersectionPoints( line, ni, pi, lx )
!                =====================   i    o   o   o
!
!        Determine the intersection points of line with integer-values grid lines
!        and sort them according to their distance from the first point
!
         type(Line_type), intent(in) :: line
         integer(INT4), intent(out) :: ni
         type(Point_type), intent(out) :: pi(:)
         logical, intent(out) :: lx(:)

         integer(INT4) :: ix1, ix2, iz1, iz2, ix, ii
         integer(INT4) :: isort(size( pi ))     ! index for ordering according to di
         real(REAL8) :: di(size( pi ))          ! distance from start point of line

         ii = 0
!
!  Determine intersection points with vertical lines (x=const)
!
         ix1 = min( line%p1%ix, line%p2%ix )
         ix2 = max( line%p1%ix, line%p2%ix ) - 1

         do ix = ix1, ix2
            ii = ii + 1
            lx(ii) = .true.
            pi(ii)%x = ix
            pi(ii)%z = line%p1%z + line%dz / line%dx * (pi(ii)%x-line%p1%x)
            pi(ii)%ix = ix
            pi(ii)%iz = ceiling( pi(ii)%z )
            di(ii) = (pi(ii)%x - line%p1%x) / line%dx * line%d
         end do
!
!  Determine intersection points with horizontal lines (z=const)
!
         iz1 = min( line%p1%iz, line%p2%iz )
         iz2 = max( line%p1%iz, line%p2%iz ) - 1

         do iz = iz1, iz2
            ii = ii + 1
            lx(ii) = .false.
            pi(ii)%z = iz
            pi(ii)%x = line%p1%x + line%dx / line%dz * (pi(ii)%z-line%p1%z)
            pi(ii)%iz = iz
            pi(ii)%ix = ceiling( pi(ii)%x )
            di(ii) = (pi(ii)%z - line%p1%z) / line%dz * line%d
         end do

         ni = ii
!
!  Sort the points according to their distance from the first polygon line point
!
         isort(1:ni) = [ (ii, ii=1,ni) ]
         call quickSort( di(1:ni), isort(1:ni), 1, ni )
!             =========    io         io        i  i

         pi(1:ni) = pi( isort(1:ni) )
         lx(1:ni) = lx( isort(1:ni) )

      end subroutine getIntersectionPoints
!
!-----------------------------------------------------------------------------------------
!
      subroutine closePolygon( ix, iz, na, pa, lxa )
!                ============  o   o   io  io   i
!
!        Close the polygon counter-clockwise along the edges of the cell
!
         integer(INT4), intent(out) :: ix
         integer(INT4), intent(out) :: iz
         integer(INT4), intent(inout) :: na
         type(Point_type), intent(inout) :: pa(:)
         logical, intent(in) :: lxa(:)

         integer(INT4) :: ie, ie1, ie2
!
!  We number the cell edges like this:
!
!                              |   4   |
!                         iz --+-------+--      z
!                              |       |       ^
!                             1|(ix,iz)|3      |
!                              |       |       |
!                       iz-1 --+-------+--     +---> x
!                              |   2   |
!                             ix-1     ix
!
!  Determine cell indices.
!
!  Note: pa(.)%ix is cell index ix except when pa is on line 1
!        pa(.)%iz is cell index iz except when pa is on line 2
!
!!!         if( lxa(1) .and. lxa(na) .and. pa(1)%ix == pa(na)%ix ) then
!!!            ix = pa(2)%ix                       ! two vertical lines with the same index =>
!!!         else                                   !    there must be a point (ia=2) inside cell
!!!            ix = max( pa(1)%ix, pa(na)%ix )     !    which determines the cell index
!!!         end if
!!!!
!!!         if( .not.lxa(1) .and. .not.lxa(na) .and. pa(1)%iz == pa(na)%iz ) then
!!!            iz = pa(2)%iz                       ! two horizontal lines with the same index =>
!!!         else                                   !    there must be a point (ia=2) inside cell
!!!            iz = max( pa(1)%iz, pa(na)%iz )     !    which determines the cell index
!!!         end if
         ix = maxval(pa(1:na)%ix)
         iz = maxval(pa(1:na)%iz)
!
!  Determine cell edges first and last point are on
!
         if( lxa(1) ) then
            if( pa(1)%ix == ix ) then
               ie1 = 3
            else
               ie1 = 1
            end if
         else
            if( pa(1)%iz == iz ) then
               ie1 = 4
            else
               ie1 = 2
            end if
         end if

         if( lxa(na) ) then
            if( pa(na)%ix == ix ) then
               ie2 = 3
            else
               ie2 = 1
            end if
         else
            if( pa(na)%iz == iz ) then
               ie2 = 4
            else
               ie2 = 2
            end if
         end if
!
!  If both points are on same edge and first point is counter clockwise from last point,
!  no additional points are needed
!
         if( ie1 == ie2 ) then                  ! first and last point on same edge
            select case( ie1 )
               case( 1 )
                  if( pa(1)%z <= pa(na)%z )  return
               case( 2 )
                  if( pa(1)%x >= pa(na)%x )  return
               case( 3 )
                  if( pa(1)%z >= pa(na)%z )  return
               case( 4 )
                  if( pa(1)%x <= pa(na)%x )  return
            end select
         end if
!
!  Starting with edge ie2, add corner points until we reach edge ie1
!
         ie = ie2

         do

            na = na + 1

            select case( ie )
               case( 1 )
                  pa(na)%x = ix - 1
                  pa(na)%z = iz - 1
               case( 2 )
                  pa(na)%x = ix
                  pa(na)%z = iz - 1
               case( 3 )
                  pa(na)%x = ix
                  pa(na)%z = iz
               case( 4 )
                  pa(na)%x = ix - 1
                  pa(na)%z = iz
            end select

            if( ie < 4 ) then
               ie = ie + 1
            else
               ie = 1
            end if

            if( ie == ie1 )  exit

         end do

      end subroutine closePolygon
!
!-----------------------------------------------------------------------------------------
!
   end subroutine getOvlpAreasCartGrid2d
!
!=========================================================================================
!
   subroutine bndConnected( mask )
!             ============   io
!
!     Among the set elements of mask determine those connected to the boundary
!     via other set elements, and unset all other elements
!
!  Arguments:
!
      logical, intent(inout) :: mask(:,:)       ! mask array
!
!  Local Variables:
!
      logical :: isEmpty
      integer(INT4) :: ix, iz, nx, nz, ixn, izn
      type(Int4Stack) :: ixStack, izStack
      logical :: mask_tmp(size( mask, 1 ),size( mask, 2 ) )
!
!  Start of Executable Code:
!  ------------------------
!
      nx = size( mask, 1 )
      nz = size( mask, 2 )

      mask_tmp(:,:) = mask(:,:)
      mask(:,:) = .false.

      call newStack( ixStack )
!          ========     o
      call newStack( izStack )
!          ========     o
!
!  Select those cells on circumference that are .false.
!
      ix = 1
      do iz = 1, nz
         if( .not.mask_tmp(ix,iz) ) cycle
         mask(ix,iz) = .true.
         call pushStack( ix, ixStack )
!             =========  i     io
         call pushStack( iz, izStack )
!             =========  i     io
      end do

      ix = nx
      do iz = 1, nz
         if( .not.mask_tmp(ix,iz) ) cycle
         mask(ix,iz) = .true.
         call pushStack( ix, ixStack )
!             =========  i     io
         call pushStack( iz, izStack )
!             =========  i     io
      end do

      iz = 1
      do ix = 2, nx-1
         if( .not.mask_tmp(ix,iz) ) cycle
         mask(ix,iz) = .true.
         call pushStack( ix, ixStack )
!             =========  i     io
         call pushStack( iz, izStack )
!             =========  i     io
      end do

      iz = nz
      do ix = 2, nx-1
         if( .not.mask_tmp(ix,iz) ) cycle
         mask(ix,iz) = .true.
         call pushStack( ix, ixStack )
!             =========  i     io
         call pushStack( iz, izStack )
!             =========  i     io
      end do
!
!  Look for mask_tmp=false cells that are connected to the mask=.false. cells determined above
!
      do

         call popStack( ixStack, ix, isEmpty )
!             ========    io     o      o
         call popStack( izStack, iz, isEmpty )
!             ========    io     o      o
         if( isEmpty ) exit

         ixn = max( 1, min( nx, ix - 1 ) )
         izn = max( 1, min( nz, iz ) )
         if( .not.mask(ixn,izn) .and. mask_tmp(ixn,izn) ) then
            mask(ixn,izn) = .true.
            call pushStack( ixn, ixStack )
!                =========   i     io
            call pushStack( izn, izStack )
!                =========   i     io
         end if

         ixn = max( 1, min( nx, ix ) )
         izn = max( 1, min( nz, iz - 1 ) )
         if( .not.mask(ixn,izn) .and. mask_tmp(ixn,izn) ) then
            mask(ixn,izn) = .true.
            call pushStack( ixn, ixStack )
!                =========   i     io
            call pushStack( izn, izStack )
!                =========   i     io
         end if

         ixn = max( 1, min( nx, ix + 1 ) )
         izn = max( 1, min( nz, iz ) )
         if( .not.mask(ixn,izn) .and. mask_tmp(ixn,izn) ) then
            mask(ixn,izn) = .true.
            call pushStack( ixn, ixStack )
!                =========   i     io
            call pushStack( izn, izStack )
!                =========   i     io
         end if

         ixn = max( 1, min( nx, ix ) )
         izn = max( 1, min( nz, iz + 1 ) )
         if( .not.mask(ixn,izn) .and. mask_tmp(ixn,izn) ) then
            mask(ixn,izn) = .true.
            call pushStack( ixn, ixStack )
!                =========   i     io
            call pushStack( izn, izStack )
!                =========   i     io
         end if
!
      end do

   end subroutine bndConnected
!
!=======================================================================
!
end module Geom_CartGrid2d_mod
