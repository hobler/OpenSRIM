module Target_Surface2d_SDF_mod
!
!  Define the surface by a signed distance function (SDF) defined on a rectilinear grid
!
!  NOTE: The indices for the distance function grid are defined consistently with
!        findZeroLevelSetCartGrid2d, but different from Target_Grid2d_*_mod
!
!  Used Modules:
!  ============
!
   use Const_mod,             only: INT4, INT8, REAL8,  &
                                    ZERO, ONE, TWO, HALF, QUART, BIG, EPSMB,  &
                                    ENUM
   use Array_mod,             only: alloc
   use Open_mod,              only: openFile

   use Param_Setup_mod,       only: ndim
   use Param_Geom_mod,        only: geom
   use Param_Ion_mod,         only: ion
   use Param_Sn_mod,          only: sn
   use Param_Out_mod,         only: out

   use Geom_CartGrid2d_mod,   only: findZeroLevelSetCartGrid2d
   use Geom_Line_mod,         only: mcdist
   use Geom_Polygon_Old_mod,  only: getDistPolygon, mcauxl, mcloc, writePolygons

   use Math_BinSearch_mod,    only: binSearch
   use Math_IntMath_mod,      only: index2to1
   use Math_Table_mod,        only: Table2d_type, initTable, interpolate

   use Stat_His2d_mod,        only: His2d_type, writeHis2d
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
   public :: initSurface2dSDF                ! initialize the grid and SDF
!            ================
   public :: initSurface2dSDFfromPolygons    ! initialize the surface from the surface
!            ============================    !    polygon(s)
   public :: initSurface2dSDFfromSegments    ! initialize surface from fill fractions and
!            ============================    !    surface segments
   public :: getNcellSurface2dSDF            ! get the number of cells
!            ====================
   public :: getGridSurface2dSDF             ! get the grid
!            ===================
   public :: getCellIndex2dSDF               ! get cell index at given position
!            =================
   public :: getCartSurface2dSDFPolygons     ! get the surface polygon in Cartesian grid
!            ===========================     !    coordinates
   public :: getDistSurface2dSDF             ! get signed distance from surface
!            ===================
   public :: writeSurface2d                  ! write the surface polygon to the output file
!            ==============                  !    as constructed from the distance function
!
!  Private Procedures:
!  ==================
!
   private :: readSurface2dSDF               ! read the SDF from sdffile
!             ================
   private :: initSurface2dSDFfromGeom       ! initialize the surface from regions
!             ========================
   private :: writeDist                      ! write the distance function to file DIST2
!             =========
!
!  Public and Private Data:
!  ============
!
   type :: Grid_type
!          =========
      logical :: regular                     ! indicates constant cell widths
      integer(INT4) :: nx                    ! number of cells in x direction
      integer(INT4) :: nz                    ! number of cells in z direction
      real(REAL8), allocatable :: wx(:)      ! (nx) cell width in x direction
      real(REAL8), allocatable :: wz(:)      ! (nz) cell width in z direction
      real(REAL8), allocatable :: x(:)       ! (0:nx) positions of x interfaces between cells
      real(REAL8), allocatable :: z(:)       ! (0:nz) positions of z interfaces between cells
   end type Grid_type

   type(Grid_type), public :: grid

   type(Table2d_type) :: sdf_tab             ! table of signed distance function
!
!=======================================================================
!
contains
!
!=======================================================================
!
   subroutine initSurface2dSDF( )
!             ================
!
!     Initialize the simulation domain and the signed distance function from the surface
!
!     Currently, reading from a 2D cellfile is not implemented
!
!  Arguments:
!
!     None
!
!  Local Variables:
!
      real(REAL8) :: wvacuum
!
!  Start of Executable Code:
!  ------------------------
!
      if( geom%usesdf ) then

         call readSurface2dSDF( )
!             ================

         ! Bounding box (we rely on free boundaries already containing a vacuum layer)
         geom%xlim(1) = grid%x(0)
         geom%xlim(2) = grid%x(grid%nx)
         geom%zlim(1) = grid%z(0)
         geom%zlim(2) = grid%z(grid%nz)

      else

         ! Add vaccuum layer at free boundaries
         wvacuum = sn%pmax * (ONE+TWO*EPSMB)
         if( geom%bndtypex(1) == ENUM%bndtype%free )  geom%xlim(1) = geom%xlim(1) - wvacuum
         if( geom%bndtypex(2) == ENUM%bndtype%free )  geom%xlim(2) = geom%xlim(2) + wvacuum
         if( geom%bndtypez(1) == ENUM%bndtype%free )  geom%zlim(1) = geom%zlim(1) - wvacuum
         if( geom%bndtypez(2) == ENUM%bndtype%free )  geom%zlim(2) = geom%zlim(2) + wvacuum

         call initSurface2dSDFfromGeom( )
!             ========================

      end if

   end subroutine initSurface2dSDF
!
!=======================================================================
!
   subroutine readSurface2dSDF( )
!             ================
!
!     Initialize the simulation domain and the signed distance function from the surface
!
!     The rectilinear grid read from sdffile(0) must cover all SDF > -sn%pmax*(1+2*epsB)
!
!  Arguments:
!
!     None
!
!  Local Variables:
!
      character(256) :: line
      integer(INT4) :: i, nd, istat, ix, iz
      real(REAL8), allocatable :: sdf(:,:)
!
!  Start of Executable Code:
!  ------------------------
!
      call openFile( 11, geom%sdffile(0), 'formatted', 'read' )
!          ========   i        i               i         i
!
!  Read header
!
      grid%nx = 0
      grid%nz = 0

      do
         read( 11, '(a)', err=900, end=900 )  line
         i = index( line, ':' )
         if( i == 0 )  cycle
         if( line(:i-1) == 'Dimension' ) then
            read( line(i+1:), '(i99)' )  nd
            if( nd /= ndim )  goto 910
         else if( line(:i-1) == 'Unit' ) then
            if( adjustl( line(i+1:) ) /= 'Angstrom' )  goto 911
         else if( line(:i-1) == 'nX' ) then
            read( line(i+1:), '(i99)' )  grid%nx
            grid%nx = grid%nx - 1
            call alloc( grid%wx, grid%nx, 'grid%wx', 'readSurface2dSDF' )
!                =====     o        i         i              i
            grid%wx(:) = ZERO
         else if( line(:i-1) == 'nY' ) then
            read( line(i+1:), '(i99)' )  grid%nz
            grid%nz = grid%nz - 1
            call alloc( grid%wz, grid%nz, 'grid%wz', 'readSurface2dSDF' )
!                =====     o        i         i              i
            grid%wz(:) = ZERO
         else if( line(:i-1) == 'dX' ) then
            if( .not.allocated( grid%wx ) )  goto 912
            grid%regular = ( index( trim( adjustl( line(i+1:) ) ), ' ' ) == 0 .and.  &
                             index( trim( adjustl( line(i+1:) ) ), ',' ) == 0 )
            if( grid%regular ) then
               read( line(i+1:), '(f99.0)' )  grid%wx(1)
               grid%wx(:) = grid%wx(1)
            else
               backspace( 11 )
               read( 11, *, iostat=istat )  line, grid%wx(:)
               if( istat /= 0 )  goto 913
            end if
         else if( line(:i-1) == 'dY' ) then
            if( .not.allocated( grid%wz ) )  goto 914
            grid%regular = ( index( trim( adjustl( line(i+1:) ) ), ' ' ) == 0 .and.  &
                             index( trim( adjustl( line(i+1:) ) ), ',' ) == 0 )
            if( grid%regular ) then
               read( line(i+1:), '(f99.0)' )  grid%wz(1)
               grid%wz(:) = grid%wz(1)
            else
               backspace( 11 )
               read( 11, *, iostat=istat )  line, grid%wz(:)
               if( istat /= 0 )  goto 915
            end if
         else if( line(:i-1) == 'Values' ) then
            exit
         end if
      end do

      if( grid%nx == 0 .or. grid%nz == 0 )  goto 919
      if( any( grid%wx(:) == ZERO ) .or. any( grid%wz(:) == ZERO ) )  goto 919
!
!  Grid positions
!
      call alloc( grid%x, 0, grid%nx, 'grid%x', 'initSurface2dSDF' )
!          =====    o     i     i        i               i
      grid%x(0) = - HALF * sum( grid%wx(:) )
      do ix = 1, grid%nx
         grid%x(ix) = grid%x(ix-1) + grid%wx(ix)
      end do

      call alloc( grid%z, 0, grid%nz,  'grid%z', 'initSurface2dSDF' )
!          =====    o     i     i         i              i
      grid%z(0) = ZERO
      do iz = 1, grid%nz
         grid%z(iz) = grid%z(iz-1) + grid%wz(iz)
      end do
!
!  SDF values
!
      call alloc( sdf, [0, 0], [grid%nx, grid%nz], 'sdf', 'initSurface2dSDF' )
!          =====   o     i             i             i            i
      do ix = 0, grid%nx
         read( 11, * )  ( sdf(ix,iz), iz=0,grid%nz )
      end do

      close( 11 )
!
!  Initialize SDF table
!
      call initTable( sdf_tab, grid%x(:), grid%z(:), sdf(:,:) )
!          =========     o        i          i          i

      if( out%ldist )  call writeDist( )
!                           =========

      return

900   print *, 'readSurface2dSDF: Error reading file ''', geom%sdffile(0), ''''
      error stop

910   error stop 'readSurface2dSDF: nDim=2 required on sdffile'
911   error stop 'readSurface2dSDF: Units=Angstrom required on sdffile'
912   error stop 'readSurface2dSDF: nX must be specified before dX in sdffile'
913   error stop 'readSurface2dSDF: error when reading dX record of sdffile'
914   error stop 'readSurface2dSDF: nY must be specified before dY in sdffile'
915   error stop 'readSurface2dSDF: error when reading dY record of sdffile'
919   error stop 'readSurface2dSDF: Incomplete grid definition'

   end subroutine readSurface2dSDF
!
!=======================================================================
!
   subroutine initSurface2dSDFfromGeom( )
!             ========================
!
!     Initialize the signed distance function from a set of polygons defining the surface.
!     The surface may be disconnected, i.e. composed of a set of connected components.
!
!  Arguments:
!
!     None
!
!  Local Variables:
!
      integer(INT4) :: ix, ixmin, ixmax, nx, iz, izmin, izmax, nz
      real(REAL8) :: pxmin, pxmax, xmin, xmax, pzmin, pzmax, zmin, zmax
      real(REAL8), allocatable :: x(:), z(:)
!
!  Start of Executable Code:
!  ------------------------
!
!  Determine minimum boundaries of grid area (not including reflected targets)
!
      if( geom%bndtypex(1) == ENUM%bndtype%reflect ) then
         xmin = geom%xlim(1)
         xmax = geom%xlim(3)
      else if( geom%bndtypex(2) == ENUM%bndtype%reflect ) then
         xmin = geom%xlim(3)
         xmax = geom%xlim(2)
      else
         xmin = geom%xlim(1)
         xmax = geom%xlim(2)
      end if

      if( geom%bndtypez(1) == ENUM%bndtype%reflect ) then
         zmin = geom%zlim(1)
         zmax = geom%zlim(3)
      else if( geom%bndtypez(2) == ENUM%bndtype%reflect ) then
         zmin = geom%zlim(3)
         zmax = geom%zlim(2)
      else
         zmin = geom%zlim(1)
         zmax = geom%zlim(2)
      end if
!
!  Determine grid from minimum boundaries, geom%wbox2, and target geometry
!
      pxmin = minval( geom%region(0)%px(:) )
      pxmax = maxval( geom%region(0)%px(:) )
      ixmin = - ceiling( (pxmin-xmin) / geom%wbox2(1) )
      ! first term is close to integer, but may be slightly less or more due to rounding error
      ! do not join the two terms!
      ixmax = nint( (pxmax-pxmin) / geom%wbox2(1) ) + ceiling( (xmax-pxmax) / geom%wbox2(1) )
      nx = ixmax - ixmin

      call alloc( x, 0, nx, 'x', 'initSurface2dSDFfromGeom' )
!          =====  o  i  i    i               i
      x(:) = [ ( pxmin+ix*geom%wbox2(1), ix=ixmin,ixmax ) ]

      pzmin = minval( geom%region(0)%pz(:) )
      pzmax = maxval( geom%region(0)%pz(:) )
      izmin = - ceiling( (pzmin-zmin) / geom%wbox2(2) )
      ! first term is close to integer, but may be slightly less or more due to rounding error
      ! do not join the two terms!
      izmax = nint( (pzmax-pzmin) / geom%wbox2(2) ) + ceiling( (zmax-pzmax) / geom%wbox2(2) )
      nz = izmax - izmin

      call alloc( z, 0, nz, 'z', 'initSurface2dSDFfromGeom' )
!          =====  o  i  i    i               i
      z(:) = [ ( pzmin+iz*geom%wbox2(2), iz=izmin,izmax ) ]
!
!  Initialize the SDF from the surface polygon. Lines marked as virtual are ignored
!
      call initSurface2dSDFfromPolygons( x(:), z(:), [geom%region(0)%np+1],  &
!          ============================  i     i              i
                                         geom%region(0)%px(:), geom%region(0)%pz(:),  &
!                                                i                     i
                                         geom%region(0)%virtual(:) )
!                                                   i
   end subroutine initSurface2dSDFfromGeom
!
!=======================================================================
!
   subroutine initSurface2dSDFfromPolygons( x, z, ipend, px, pz, virtual_ )
!             ============================  i  i    i    i   i      i
!
!     Initialize the signed distance function from a set of polygons defining the surface.
!     The surface may be disconnected, i.e. composed of a set of connected components.
!
!  Arguments:
!
      real(REAL8), intent(in) :: x(0:)             ! x-coordinates of rectilinear grid for
                                                   !    distance values
      real(REAL8), intent(in) :: z(0:)             ! z-coordinates of rectilinear grid for
                                                   !    distance values
      integer(INT4), intent(in) :: ipend(:)        ! last point index of connected component
      real(REAL8), target, intent(in) :: px(:)     ! x-coordinates of polygon corners
      real(REAL8), target, intent(in) :: pz(:)     ! z-coordinates of polygon corners
                                                   !    in each polygon last=first point
      logical, target, intent(in), optional :: virtual_(:)  ! flag indicating virtual (non-
                                                            !    existing) segment
!
!  Local Variables:
!
      integer(INT4) :: ii, jj, i1, i2, j1, j2, ix, iz
      integer(INT4) :: icc, ncc, iloc, isign, ipbeg(size( ipend ))
      real(REAL8), allocatable :: sdf(:,:)

      real(REAL8) :: xx, zz, dist, distsq, distsq_tmp
      real(REAL8), pointer :: pccx(:), pccz(:)
      logical, pointer :: virtualcc(:)
      real(REAL8) :: coseta(size( ipend )), sineta(size( ipend ))
!
!  Start of Executable Code:
!  ------------------------
!
      grid%nx = ubound( x, 1 )
      grid%nz = ubound( z, 1 )

      call alloc( grid%x, 0, grid%nx, 'grid%x', 'initSurface2dSDFfromPolygons' )
!          =====    o     i     i        i                    i
      call alloc( grid%z, 0, grid%nz, 'grid%z', 'initSurface2dSDFfromPolygons' )
!          =====    o     i     i        i                    i
      grid%x(:) = x(:)
      grid%z(:) = z(:)

      grid%wx = x(1:) - x(:grid%nx-1)
      grid%wz = z(1:) - z(:grid%nz-1)

      grid%regular = all( abs( grid%wx(:) - grid%wx(1) ) < EPSMB * abs( grid%wx(1) ) ) .and.  &
                     all( abs( grid%wz(:) - grid%wz(1) ) < EPSMB * abs( grid%wz(1) ) )
!
!  Number of connected components and begin index of each connected component
!
      ncc = size( ipend )
      ipbeg(1) = 1
      ipbeg(2:ncc) = ipend(1:ncc-1)+1
!
!  Determine auxiliary direction for each connected component
!
      do icc = 1, ncc
         pccx => px(ipbeg(icc):ipend(icc))
         pccz => pz(ipbeg(icc):ipend(icc))
         call mcauxl( coseta(icc), sineta(icc), pccx, pccz )
!             ======      o            o         i     i
      end do
!
!  Calculate signed distance function
!
      if( geom%bndtypex(1) == ENUM%bndtype%periodic ) then
         i1 = -1
      else
         i1 = 0
      end if
      if( geom%bndtypex(2) == ENUM%bndtype%periodic ) then
         i2 = 1
      else
         i2 = 0
      end if
      if( geom%bndtypez(1) == ENUM%bndtype%periodic ) then
         j1 = -1
      else
         j1 = 0
      end if
      if( geom%bndtypez(2) == ENUM%bndtype%periodic ) then
         j2 = 1
      else
         j2 = 0
      end if

      call alloc( sdf, [0, 0], [grid%nx, grid%nz], 'sdf', 'initSurface2dSDFfromPolygons' )
!          =====   o     i             i             i                  i

      do ix = 0, grid%nx
         do iz = 0, grid%nz

            dist = BIG                          ! distance from surface
            isign = -1                          ! outside target

            do icc = 1, ncc
               pccx => px(ipbeg(icc):ipend(icc))
               pccz => pz(ipbeg(icc):ipend(icc))
!               do i = 1, size( pccx )
!                  print *, 'px, pz =', pccx(i), pccz(i)
!               end do
!               print *, '------'
               distsq = BIG
               do ii = i1, i2
                  do jj = j1, j2
                     xx = x(ix) + ii * (geom%xlim(2) - geom%xlim(1))
                     zz = z(iz) + jj * (geom%zlim(2) - geom%zlim(1))
                     if( present( virtual_ ) ) then
                        virtualcc => virtual_(ipbeg(icc):ipend(icc))
                        call getDistPolygon( distsq_tmp, xx, zz, pccx, pccz, virtualcc )
!                            ==============      o       i   i    i     i        i
                     else
                        call getDistPolygon( distsq_tmp, xx, zz, pccx, pccz )
!                            ==============      o       i   i    i     i
                     end if
                     distsq = min( distsq_tmp, distsq )
                  end do
               end do

               dist = min( dist, sqrt( distsq ) )

               iloc = mcloc( x(ix), z(iz), pccx, pccz, coseta(icc), sineta(icc) )
!                     =====    i      i      i     i       i            i
               if( iloc == 0 )  iloc = 1
               iloc = - iloc                    ! -1 if inside current polygon
               isign = isign * iloc             ! -1 if inside even number of polygons
            end do

            sdf(ix,iz) = isign * dist

         end do
      end do
!
!  Convert signed distance function to table
!
      call initTable( sdf_tab, x(:), z(:), sdf(:,:) )
!          =========     o      i     i       i

      if( out%ldist )  call writeDist( )
!                           =========

   end subroutine initSurface2dSDFfromPolygons
!
!=======================================================================
!
   subroutine initSurface2dSDFfromSegmentsOld( x, z, ffgeom, segments )
!             ============================  i  i    i        i
!
!     Initialize the signed distance function from the cell fill fractions and a set of
!     surface segments. The segments need not be connected. The fill fractions are needed
!     to determine the sign of the SDF.
!
!     The rectilinear grid must cover all SDF > -sn%pmax*(1+2*epsB)
!
      real(REAL8), intent(in) :: x(0:)             ! (0:nx) x coordinates of rectilinear grid
                                                   !    for distance values
      real(REAL8), intent(in) :: z(0:)             ! (0:nz) z coordinates of rectilinear grid
                                                   !    for distance values
      real(REAL8), intent(in) :: ffgeom(:,:)       ! (nx,nz) fill fractions
      real(REAL8), target, intent(in) :: segments(:,:,:) ! (2,2,ns) surface segments
!
!  Local Variables
!
      integer(INT4) :: is, ns, ix, iz, dix, diz, jx, jz
      real(REAL8) :: dist2
      real(REAL8), pointer :: px1(:), px2(:), pz1(:), pz2(:)
      real(REAL8), allocatable :: ffinterp(:,:)
      real(REAL8), allocatable :: sdf(:,:)
!
!  Start of Executable Code:
!  ------------------------
!
!  Set up the grid
!
      grid%nx = ubound( x, 1 )
      grid%nz = ubound( z, 1 )

      grid%x = x(:)
      grid%z = z(:)

      grid%wx = x(1:) - x(:grid%nx-1)
      grid%wz = z(1:) - z(:grid%nz-1)

      grid%regular = all( abs( grid%wx(:) - grid%wx(1) ) < EPSMB * abs( grid%wx(1) ) ) .and.  &
                     all( abs( grid%wz(:) - grid%wz(1) ) < EPSMB * abs( grid%wz(1) ) )
!
!  Interpolate the fill fractions to the grid points
!
      call alloc( ffinterp, [0, 0], [grid%nx, grid%nz], 'ffinterp',  &
!          =====     io       i             i               i
                  'initSurface2dSDFfromSegmentsOld' )
!                               i
      do iz = 0, grid%nz
         do ix = 0, grid%nx
            ffinterp(ix,iz) = ZERO
            do diz = 0, 1
               jz = max( 1, min( grid%nz, iz+diz ) )
               do dix = 0, 1
                  jx = max( 1, min( grid%nx, ix+dix ) )
                  ffinterp(ix,iz) = ffinterp(ix,iz) + ffgeom(jx,jz)
               end do
            end do
            ffinterp(ix,iz) = QUART * ffinterp(ix,iz)
         end do
      end do
!
!  Determine the signed distance function
!
      ns = size( segments, 3 )
      px1 => segments(1,1,:)
      px2 => segments(1,2,:)
      pz1 => segments(2,1,:)
      pz2 => segments(2,2,:)

      call alloc( sdf, [0, 0],  [grid%nx, grid%nz], 'sdf', 'initSurface2dSDFfromSegmentsOld' )
!          =====   o     i              i             i                  i
      do iz = 0, grid%nz
         do ix = 0, grid%nx
            sdf(ix,iz) = BIG
            do is = 1, ns
               call mcdist( dist2, x(ix), z(iz), px1(is), pz1(is), px2(is), pz2(is) )
!                   ======    o      i      i       i        i        i        i
               sdf(ix,iz) = min( sdf(ix,iz), dist2 )
            end do
            if( ffinterp(ix,iz) > HALF ) then
               sdf(ix,iz) = sqrt( sdf(ix,iz) )
            else
               sdf(ix,iz) = - sqrt( sdf(ix,iz) )
            end if
         end do
      end do
!
!  Convert signed distance function to table
!
      call initTable( sdf_tab, x(:), z(:), sdf(:,:) )
!          =========     o      i     i       i

      if( out%ldist )  call writeDist( )
!                           =========

   end subroutine initSurface2dSDFfromSegmentsOld
!
!=======================================================================
!
   subroutine initSurface2dSDFfromSegments( x, z, ffgeom, segments )
!             ============================  i  i    i        i
!
!     Initialize the signed distance function from the cell fill fractions and a set of
!     surface segments. The segments need not be connected, but they must be oriented so
!     that the target is on the left-hand side of the segment. The end points of the
!     segments must lie on the grid lines.
!
      real(REAL8), intent(in) :: x(0:)             ! (0:nx) x coordinates of rectilinear grid
                                                   !    for distance values
      real(REAL8), intent(in) :: z(0:)             ! (0:nz) z coordinates of rectilinear grid
                                                   !    for distance values
      real(REAL8), intent(in) :: ffgeom(:,:)       ! (nx,nz) fill fractions
      real(REAL8), target, intent(in) :: segments(:,:,:) ! (2,2,ns) surface segments
!
!  Local Variables
!
      integer(INT4) :: is, ns, ixc, izc, ix, iz, ix1, iz1, ix2, iz2, jx, jz, jsign
      integer(INT4) :: ii, jj, i1, i2, j1, j2
      integer(INT4) :: ixpos_min, ixneg_min, dix, diz
      integer(INT4), allocatable :: isign(:,:), counter(:,:)
      real(REAL8) :: px1, px2, pz1, pz2, dist2
      real(REAL8), allocatable :: ff(:,:), sdf(:,:)
!
!  Start of Executable Code:
!  ------------------------
!
!  Set up the grid
!
      grid%nx = ubound( x, 1 )
      grid%nz = ubound( z, 1 )

      if( allocated( grid%x ) )  deallocate( grid%x )
      if( allocated( grid%z ) )  deallocate( grid%z )
      call alloc( grid%x, 0, grid%nx, 'grid%x', 'initSurface2dSDFfromSegments' )
!          =====    o     i     i        i                   i
      call alloc( grid%z, 0, grid%nz, 'grid%z', 'initSurface2dSDFfromSegments' )
!          =====    o     i     i        i                   i
      grid%x(:) = x(:)
      grid%z(:) = z(:)

      grid%wx = x(1:) - x(:grid%nx-1)
      grid%wz = z(1:) - z(:grid%nz-1)

      grid%regular = all( abs( grid%wx(:) - grid%wx(1) ) < EPSMB * abs( grid%wx(1) ) ) .and.  &
                     all( abs( grid%wz(:) - grid%wz(1) ) < EPSMB * abs( grid%wz(1) ) )

      call alloc( sdf, [0, 0], [grid%nx, grid%nz], 'sdf', 'initSurface2dSDFfromSegments' )
!          =====   o     i             i             i                  i
      sdf(:,:) = ZERO
!
!  Special case: if there are no segments
!  (could happen if everything is sputtered away or if we start with an infinite medium)
!
      ns = size( segments, 3 )

      if( ns == 0 ) then
         if( any( geom%bndtypex(:) == ENUM%bndtype%free ) .or.  &
             any( geom%bndtypez(:) /= ENUM%bndtype%free ) ) then     ! completely empty
            sdf(:,:) = - BIG
         else                                                        ! completely full
            sdf(:,:) = BIG
         end if
         call initTable( sdf_tab, x(:), z(:), sdf(:,:) )
!             =========     o      i     i       i
         return
      end if
!
!  Interpolate the fill fractions to the grid points
!
      call alloc( ff, [0, 0], [grid%nx, grid%nz], 'ff', 'initSurface2dSDFfromSegments' )
!          =====  o     i             i            i                   i
      do iz = 0, grid%nz
         do ix = 0, grid%nx
            ff(ix,iz) = ZERO
            do diz = 0, 1
               jz = max( 1, min( grid%nz, iz+diz ) )
               do dix = 0, 1
                  jx = max( 1, min( grid%nx, ix+dix ) )
                  ff(ix,iz) = ff(ix,iz) + ffgeom(jx,jz)
               end do
            end do
            ff(ix,iz) = QUART * ff(ix,iz)
         end do
      end do
!
!  Add the signs (+1 for inside, -1 for outside) and SDFs on all grid points that are corners
!  of cells containing segments
!
      call alloc( isign, [0, 0], [grid%nx, grid%nz], 'isign', 'initSurface2dSDFfromSegments' )
!          =====    o      i             i              i                  i
      isign(:,:) = 0
      call alloc( counter, [0, 0], [grid%nx, grid%nz], 'counter',  &
!          =====     o       i             i               i
                  'initSurface2dSDFfromSegments' )
!                                i
      counter(:,:) = 0

      do is = 1, ns
         ! scaled coordinates of segment end points (NOTE: this only works for regular grids)
         px1 = (segments(1,1,is)-grid%x(0)) / grid%wx(1)
         px2 = (segments(1,2,is)-grid%x(0)) / grid%wx(1)
         pz1 = (segments(2,1,is)-grid%z(0)) / grid%wz(1)
         pz2 = (segments(2,2,is)-grid%z(0)) / grid%wz(1)

         ! cell indices of segment
         ixc = HALF*(px1+px2) + 1
         izc = HALF*(pz1+pz2) + 1

         call backCorner( ix1, iz1, px1, pz1, ixc, izc )
!             ----------   o    o    i    i    i    i
         call backCorner( ix2, iz2, px2, pz2, ixc, izc )
!             ----------   o    o    i    i    i    i
         ix = ix1
         iz = iz1
         jsign = +1                                         ! starting inside
         px1 = segments(1,1,is)
         px2 = segments(1,2,is)
         pz1 = segments(2,1,is)
         pz2 = segments(2,2,is)
         do
            call nextCorner( ix, iz, ixc, izc )
!                ----------  io  io   i    i
            call mcdist( dist2, grid%x(ix), grid%z(iz), px1, pz1, px2, pz2 )
!                ======    o        i           i        i    i    i    i
            isign(ix,iz) = isign(ix,iz) + jsign
            counter(ix,iz) = counter(ix,iz) + 1
            sdf(ix,iz) = sdf(ix,iz) + jsign*sqrt( dist2 )
            if( ix == ix2 .and. iz == iz2 )  jsign = -1     ! going outside
            if( ix == ix1 .and. iz == iz1 )  exit
         end do
      end do

      if( any( geom%bndtypex(:) == ENUM%bndtype%periodic ) ) then
         do iz = 0, grid%nz
            counter(0,iz) = counter(0,iz) + counter(grid%nx,iz)
            isign(0,iz) = isign(0,iz) + isign(grid%nx,iz)
            counter(grid%nx,iz) = counter(0,iz)
            isign(grid%nx,iz) = isign(0,iz)
         end do
         i1 = -1
         i2 = +1
      else
         i1 = 0
         i2 = 0
      end if
      if( any( geom%bndtypez(:) == ENUM%bndtype%periodic ) ) then
         do ix = 0, grid%nx
            counter(ix,0) = counter(ix,0) + counter(ix,grid%nz)
            isign(ix,0) = isign(ix,0) + isign(ix,grid%nz)
            counter(ix,grid%nz) = counter(ix,0)
            isign(ix,grid%nz) = isign(ix,0)
         end do
         j1 = -1
         j2 = +1
      else
         j1 = 0
         j2 = 0
      end if

      where( counter(:,:) > 0 )
         isign(:,:) = max( -1, min( 1, isign(:,:) ) )       ! clip to [-1, +1]
         sdf(:,:) = sdf(:,:) / counter(:,:)                 ! convert sum to average
      end where
!
!  Determine sign of corner grid point ix=iz=0
!
      do iz = 0, grid%nz
         if( any( isign(:,iz) /= 0 ) )  exit
      end do
      ixpos_min = findloc( isign(:,iz) > 0, .true., dim=1 )
      ixneg_min = findloc( isign(:,iz) < 0, .true., dim=1 )
      if( ixpos_min == 0 ) then
         isign(0,0) = -1
      else if( ixneg_min == 0 ) then
         isign(0,0) = +1
      else if( ixneg_min < ixpos_min ) then
         isign(0,0) = -1
      else if( ixpos_min < ixneg_min ) then
         isign(0,0) = +1
      end if
!
!  Calculate SDFs as the signed minimum of the distances from all segments, keeping averages
!  calculated above where counter > 0 and isign == 0
!
      do iz = 0, grid%nz

         do ix = 0, grid%nx

            if( counter(ix,iz) == 0 ) then
               if( ff(ix,iz) > HALF ) then               ! take sign from fill fraction
                  isign(ix,iz) = +1
               else
                  isign(ix,iz) = -1
               end if
            else if( isign(ix,iz) == 0 ) then
               cycle                                     ! SDF already calculated
            else
               continue                                  ! sign already calculated
            end if

            sdf(ix,iz) = BIG
            do ii = i1, i2
               do jj = j1, j2
                  do is = 1, ns
                     px1 = segments(1,1,is) + ii * ( geom%xlim(2) - geom%xlim(1) )
                     px2 = segments(1,2,is) + ii * ( geom%xlim(2) - geom%xlim(1) )
                     pz1 = segments(2,1,is) + jj * ( geom%zlim(2) - geom%zlim(1) )
                     pz2 = segments(2,2,is) + jj * ( geom%zlim(2) - geom%zlim(1) )
                     call mcdist( dist2, grid%x(ix), grid%z(iz), px1, pz1, px2, pz2 )
      !                   ======    o        i           i        i    i    i    i
                     sdf(ix,iz) = min( dist2, sdf(ix,iz) )
                  end do
               end do
            end do
            sdf(ix,iz) = isign(ix,iz) * sqrt( sdf(ix,iz) )

         end do

      end do
!
!  Convert signed distance function to table
!
      call initTable( sdf_tab, x(:), z(:), sdf(:,:) )
!          =========     o      i     i       i
!
!  Set domain boundaries
!
      if( geom%bndtypex(1) == ENUM%bndtype%free )  geom%xlim(1) = grid%x(0)
      if( geom%bndtypex(2) == ENUM%bndtype%free )  geom%xlim(2) = grid%x(grid%nx)
      if( geom%bndtypez(1) == ENUM%bndtype%free )  geom%zlim(1) = grid%z(0)
      if( geom%bndtypez(2) == ENUM%bndtype%free )  geom%zlim(2) = grid%z(grid%nz)

   contains
!
!-----------------------------------------------------------------------------
!
      subroutine nextCorner( ix, iz, ixc, izc )
!                ----------  io  io   i   i
!
!     Get indices of next corner in counter-clockwise direction around cell (ixc,izc)
!
         integer(INT4), intent(inout) :: ix, iz    ! corner indices
         integer(INT4), intent(inout) :: ixc, izc  ! cell indices

         if( ix == ixc-1 .and. iz == izc-1 ) then
            iz = izc
         else if( ix == ixc-1 .and. iz == izc ) then
            ix = ixc
         else if( ix == ixc .and. iz == izc ) then
            iz = izc - 1
         else if( ix == ixc .and. iz == izc-1 ) then
            ix = ixc - 1
         else
            error stop 'initSurface2dSDFfromSegments::nextCorner: illegal indices passed'
         end if

      end subroutine nextCorner
!
!-----------------------------------------------------------------------------
!
      subroutine backCorner( ix, iz, x, z, ixc, izc )
!                ----------   o   o  i  i   i    i
!
!     Get indices of corner counter-clockwise in backward direction (i.e., clockwise) of a
!     point on the cell edge
!
         integer(INT4), intent(out) :: ix, iz   ! corner indices
         real(REAL8), intent(in) :: x, z        ! point coordinates
         integer(INT4), intent(in) :: ixc, izc  ! cell indices

         if( x < (ixc-1)*(ONE-EPSMB) .or. x > ixc*(1+EPSMB) .or.  &
             z < (izc-1)*(ONE-EPSMB) .or. z > izc*(1+EPSMB) )  &
               error stop 'initSurface2dSDFfromSegments::backCorner: point not within cell'

         if( abs( x - (ixc-1) ) < EPSMB ) then
            ix = ixc - 1
            iz = izc - 1
         else if( abs( z - izc ) < EPSMB ) then
            ix = ixc - 1
            iz = izc
         else if( abs( x - ixc ) < EPSMB ) then
            ix = ixc
            iz = izc
         else if( abs( z - (izc-1) ) < EPSMB ) then
            ix = ixc
            iz = izc - 1
         else
            error stop  'initSurface2dSDFfromSegments::backCorner: point not on cell edge'
         end if

      end subroutine backCorner
!
!-----------------------------------------------------------------------------
!
   end subroutine initSurface2dSDFfromSegments
!
!=======================================================================
!
   subroutine getNcellSurface2dSDF( ncell )
!             ====================    o
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
      ncell = int( grid%nx, INT8 ) * int( grid%nz, INT8 )

   end subroutine getNcellSurface2dSDF
!
!=======================================================================
!
   subroutine getGridSurface2dSDF( x, z )
!             ===================  o  o
!
!     Get the grid
!
!  Arguments:
!
      real(REAL8), allocatable, intent(out) :: x(:)   ! x grid
      real(REAL8), allocatable, intent(out) :: z(:)   ! z grid
!
!  Local Variables:
!
!     None
!
!  Start of Executable Code:
!  ------------------------
!
      call alloc( x, 0, grid%nx, 'x', 'getGridSurface2dSDF' )
!          =====  o  i     i      i             i
      x = grid%x(:)

      call alloc( z, 0, grid%nz, 'z', 'getGridSurface2dSDF' )
!          =====  o  i     i      i             i
      z = grid%z(:)

   end subroutine getGridSurface2dSDF
!
!=======================================================================
!
   subroutine getCellIndex2dSDF( x, z, ic, outside, what_ )
!             =================  i  i  o      o       i
!
!     Get the cell index at position [x,z]
!
!  Arguments:
!
      real(REAL8), intent(in) :: x              ! x coordinate of point
      real(REAL8), intent(in) :: z              ! z coordinate of point
      integer(INT8), intent(out) :: ic          ! cell index
      logical, intent(out) :: outside           ! flag indicating point is outside grid
      character(*), intent(in), optional :: what_  ! indication what point is checked
!
!  Local Variables:
!
      integer(INT4) :: ix, iz
      real(REAL8) :: xtmp, ztmp
!
!  Start of Executable Code:
!  ------------------------
!
      if( grid%regular ) then
         xtmp = max( grid%x(0)-HALF*grid%wx(1), min( grid%x(grid%nx)+HALF*grid%wx(1), x ) )
         ztmp = max( grid%z(0)-HALF*grid%wz(1), min( grid%z(grid%nz)+HALF*grid%wz(1), z ) )
         ix = ceiling( ( xtmp - grid%x(0) ) / grid%wx(1) )  ! ix can be [0,nx+1]
         iz = ceiling( ( ztmp - grid%z(0) ) / grid%wz(1) )  ! iz can be [0,nz+1]
         outside = ( any( [ix,iz] == 0 ) .or. any( [ix,iz] == [grid%nx,grid%nz]+1 ) )
         if( outside ) then
            ix = max( 1, min( grid%nx, ix ) )
            iz = max( 1, min( grid%nz, iz ) )
         end if
      else
         call binSearch( grid%nx+1, grid%x(0:), x, ix )
!             =========      i          i       i  o
         call binSearch( grid%nz+1, grid%z(0:), z, iz )
!             =========      i          i       i  o
         outside = ( x < grid%x(0) .or. x > grid%x(grid%nx) .or.  &
                     z < grid%z(0) .or. z > grid%z(grid%nz) )
      end if

      call index2to1( ix, iz, grid%nx, ic )
!          =========  i   i      i     o

   end subroutine getCellIndex2dSDF
!
!=======================================================================
!
   subroutine getCartSurface2dSDFPolygons( ipend, px, pz )
!             ===========================    o    o   o
!
!     Get the surface polygon in Cartesian coordinates
!
!  Arguments:
!
      integer(INT4), allocatable :: ipend(:)     ! last point index of connected component
      real(REAL8), allocatable :: px(:)          ! x-coordinates of polygon corners
      real(REAL8), allocatable :: pz(:)          ! z-coordinates of polygon corners
!
!  Local Variables:
!
!     None
!
!  Start of Executable Code:
!  ------------------------
!
      call findZeroLevelSetCartGrid2d( sdf_tab%a00(:,:), ipend, px, pz )
!          ==========================         i            o    o   o

   end subroutine getCartSurface2dSDFPolygons
!
!=======================================================================
!
   subroutine getDistSurface2dSDF( x, z, sdf, dirx_, dirz_ )
!             ===================  i  i   o     o      o
!
!     Get the signed distance of a point [x,z] from the surface by interpolation
!     ( + inside, - outside ) and
!     the direction towards/from the surface (pointing outwards the target)
!
!     NOTE: It is assumed that the point is inside the grid
!
!  Arguments:
!
      real(REAL8), intent(in) :: x                 ! x coordinate of point
      real(REAL8), intent(in) :: z                 ! z coordinate of point
      real(REAL8), intent(out) :: sdf              ! signed distance: (+) inside, (-) outside
      real(REAL8), optional, intent(out) :: dirx_  ! x/z component of direction towards/from
      real(REAL8), optional, intent(out) :: dirz_  !    surface
!
!  Local Variables:
!
!     None
!
!  Start of Executable Code:
!  ------------------------
!
      call interpolate( sdf_tab, x, z, sdf, dirx_, dirz_ )
!          ===========     i     i  i   o    o      o

   end subroutine getDistSurface2dSDF
!
!=======================================================================
!
   subroutine writeSurface2d( iion )
!             ==============    i
!
!     Write the contour as a set of closed surfaces (connected components) to the SURF file.
!     The surface polygon is constructed from the distance function
!
!  Arguments:
!
      integer(INT4), intent(in) :: iion         ! index of ion
!
!  Local Variables:
!
      integer(INT4), allocatable :: ipend(:)    ! last point index of connected component
      real(REAL8), allocatable :: px(:)         ! x-coordinates of polygon corners
      real(REAL8), allocatable :: pz(:)         ! z-coordinates of polygon corners
!
      integer(INT4) :: ix, iz, ip
      real(REAL8) :: dx, dz
!
!  Start of Executable Code:
!  ------------------------
!
      call getCartSurface2dSDFPolygons( ipend, px, pz )
!          ===========================    o    o   o

      do ip = 1, size( px )
         ix = floor( px(ip) )
         iz = floor( pz(ip) )
         dx = px(ip) - ix
         dz = pz(ip) - iz                       ! the following is not nice:
         px(ip) = sdf_tab%x(ix+1) + dx * ( sdf_tab%x(ix+2) - sdf_tab%x(ix+1) )
         pz(ip) = sdf_tab%y(iz+1) + dz * ( sdf_tab%y(iz+2) - sdf_tab%y(iz+1) )
      end do

      call writePolygons( ipend, px, pz, iion*ion%ddose)
!          =============    i    i   i         i

      deallocate( ipend, px, pz )

   end subroutine writeSurface2d
!
!=======================================================================
!
   subroutine writeDist( )
!             =========
!
!     Write the distance from the surface to the DIST2 file (same format as HIS2 file)
!
!  Arguments:
!
!     None
!
!  Local Variables:
!
      integer(INT4) :: nx, nz
      type(His2d_type) :: his
      real(REAL8) :: wx, wz

      character(6), save :: fmode = 'write'
!
!  Start of Executable Code:
!  ------------------------
!
!      print *, 'x =', grid%x
!      print *, 'z =', grid%z
!      print *, 'dist =', grid%dist(grid%nx/2-30,:)
!      error stop

      nx = size( sdf_tab%x )
      nz = size( sdf_tab%y )
      his%nbox(:) = [nx, nz]
      his%imin(:) = 1
      his%imax(:) = [nx, nz]
      wx = sdf_tab%wx
      wz = sdf_tab%wy
      his%hbeg(:) = [grid%x(0)-HALF*grid%wx(1), grid%z(0)-HALF*grid%wz(1)]
      his%wbox(:) = [wx, wz]

      allocate( his%hist( 1, grid%nx+1, grid%nz+1 ) )
      his%hist(1,:,:) = sdf_tab%a00(:,:)

      allocate( his%indhis( 1 ) )
      his%indhis(:) = 0

      call writeHis2d( 'DIST2', fmode, his )
!          ==========     i       i     i

      deallocate( his%hist )
      deallocate( his%indhis )

      fmode = 'append'

   end subroutine writeDist
!
!=======================================================================
!
end module Target_Surface2d_SDF_mod
