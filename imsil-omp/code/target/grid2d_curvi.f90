module Target_Grid2d_Curvi_mod
!
!  Defines a 2d regular grid and provides subroutines operating on it
!
!  Used Modules:
!  ============
!
   use Const_mod,                only: INT4, REAL8,  &
                                       ZERO, QUART, HALF, ONE, TWO, FOUR, SIX,  &
                                       BIG, EPSMB, EPSMA, SMALL, MATHPI => PI,  &
                                       ENUM
   use Array_mod,                only: check_allocation
   use Undefine_mod,             only: defined

   use Open_mod,                 only: openFile

   use Param_Setup_mod,          only: nelem
   use Param_Geom_mod,           only: geom
   use Param_Elem_mod,           only: elem
   use Param_Ion_mod,            only: ion
   use Param_Mat_mod,            only: mat, irxtal
   use Param_Dam_mod,            only: dam
   use Param_Sn_mod,             only: sn
   use Param_Out_mod,            only: out

   use Types_mod,                only: Proj_type
   use Geom_Line_mod,            only: mcisct
   use Geom_Polygon_Old_mod,     only: mcauxl, mcloc, relPosPolygon, getPolyArea,  &
                                       getDistPolyPoly, pointInsidePolygon
   use Geom_FFPcyl_mod,          only: isInsideFFPcyl
   use Math_Sort_mod,            only: quickSort
   use Target_Surface2d_mod,     only: getDistSurface2d
   use Target_Surface2d_SDF_mod, only: initSurface2dSDFfromPolygons, getCartSurface2dSDFPolygons
#ifdef AMAT
   use Target_Mechanics2d_mod,   only: relaxAreas2d_simple, relaxViscoElastic2d
#endif
   use Geom_CartGrid2d_mod,      only: findZeroLevelSetCartGrid2d, getOvlpAreasCartGrid2d
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
   public :: initCurviGrid2d              ! initialize the target
!            ===============
   public :: getCellCurviGrid2d           ! determines the cell a given point is located in
!            ==================
   public :: getOtherNearCellsCurviGrid2d ! determine cells other than nearest cell that
!            ============================ !    intersect with the free flight path cylinder
   public :: readCurviGrid2d              ! read target grid from file
!            ===============
   public :: writeCurviGrid2d             ! write target grid to file
!            ================
   public :: relaxCurviGrid2d             ! relax the cell densities
!            ================
   public :: getCellsCurviGrid2d          ! get the cell infos
!            ===================
   public :: getOvlpAreasCurviGrid2d      ! obtain the cell areas inside the polygon
!            =======================
!
!  Private Procedures:
!  ==================
!
   private :: initOnceCurviGrid2d         ! futher initializations in addition to
!             ===================         !    initCurviGrid2d, but not after relaxCurviGrid2d
   private :: initDerivedCurviGrid2d      ! further initialization of the target after
!             ======================      !    initCurviGrid2d and relaxCurviGrid2d
   private :: icfun                       ! convert x- and z-indices into cell index
!             =====
   private :: ixfun                       ! get x-index from cell index
!             =====
   private :: izfun                       ! get z-index from cell index
!             =====
!
!  Private Data:
!  ============
!
   type :: RegCell_type
!          ============
      integer(INT4), allocatable :: icovlp(:)   ! (curvilinear) grid cells overlapping the
   end type RegCell_type                        !    (regular) grid cell

   type :: RegGrid_type                         ! regular 2d grid
!          ============
      integer(INT4) :: nx                    ! number of cells in x-direction
      integer(INT4) :: nz                    ! number of cells in z-direction
      real(REAL8) :: wx                      ! cell width in x-direction
      real(REAL8) :: wz                      ! cell width in z-direction
      real(REAL8), allocatable :: x(:)       ! (0:nx) x-coordinates of grid positions
      real(REAL8), allocatable :: z(:)       ! (0:nz) z-coordinates of grid positions
      real(REAL8) :: xcmin                   ! minimum cell center position in x-direction
      real(REAL8) :: xcmax                   ! maximum cell center position in x-direction
      real(REAL8) :: zcmin                   ! minimum cell center position in z-direction
      real(REAL8) :: zcmax                   ! maximum cell center position in z-direction
      type(RegCell_type), allocatable :: cell(:,:) ! cell information
   end type RegGrid_type

   type(RegGrid_type) :: reggrid             ! auxiliary grid for accelerating index lookup

   type :: Cell_type
!          =========
      sequence
      real(REAL8), allocatable :: xv(:)         ! x-coordinates of cell vertices (closed
                                                !    polygon)
      real(REAL8), allocatable :: zv(:)         ! z-coordinates of cell vertices (closed
                                                !    polygon)
      integer(INT4), allocatable :: icnear(:)   ! indices of near nonempty cells, ordered
                                                !    according to increasing distmin2
      real(REAL8), allocatable :: prob(:)       ! approximate probability that icnear is
                                                !    nearest nonempty cell to a point inside
                                                !    the current cell
   end type Cell_type

   type :: CurviGrid_type                       ! curvilinear 2d grid
!          ==============
      integer(INT4) :: nx                       ! number of cells in x-direction
      integer(INT4) :: nz                       ! number of cells in z-direction
      integer(INT4) :: nfrozenx(2)              ! number of frozen cells at the left/right
                                                !    boundary
      integer(INT4) :: nfrozenz(2)              ! number of frozen cells at the top/bottom
                                                !    boundary
      real(REAL8), allocatable :: x(:,:)        ! (0:nx,0:nz) x-coordinates of grid positions
      real(REAL8), allocatable :: z(:,:)        ! (0:nx,0:nz) z-coordinates of grid positions
      type(Cell_type), allocatable :: cell(:,:) ! (nx,nz) cell information
   end type CurviGrid_type

   type(CurviGrid_type) :: grid              ! main grid holding the cells (= regular grid
                                             !    except near to the surface where the nodes are
                                             !    moved to the surface)
   type(CurviGrid_type) :: relaxedgrid       ! grid after a relaxation step


   type :: Surface_type
!          ============
      integer(INT4), allocatable :: ix(:)       ! x-index of surface point
      integer(INT4), allocatable :: iz(:)       ! z-index of surface point
      integer(INT4), allocatable :: ipend(:)    ! last point of closed surface
   end type Surface_type

   type(Surface_type) :: surface             ! surface
!
!   integer(INT4), allocatable :: icell(:,:)     ! Conversion table from Cartesian indices to
                                                !    cell indices
!
!==========================================================================================
!
contains
!
!==========================================================================================
!
   subroutine initCurviGrid2d( )
!             ===============
!
!     Initialize the grid
!
!  Arguments:
!
!     None
!
!  Local Variables:
!
      character(256) :: errmsg
      integer(INT4) :: ix, iz, ixmin, izmin, ixmax, izmax, istat, np
      integer(INT4) :: ixmin1, izmin1, ixmax1, izmax1, ixmin2, izmin2, ixmax2, izmax2
      real(REAL8) :: offset(2), rixmin, rizmin, rixmax, rizmax
      real(REAL8), allocatable :: xc(:), zc(:)
!
!  Start of Executable Code:
!  ------------------------
!
!  Regular grid
!
      reggrid%wx = out%wbox2(1)
      reggrid%wz = out%wbox2(2)
!
!  number of cells
!                                               ! avoid cell boundaries at lattice planes
      if( irxtal > 0 ) then
         offset(:) = -100 * EPSMA * [ maxval(geom%region(0)%px(:)),  &
                                       maxval(geom%region(0)%pz(:)) ]
      else
         offset(:) = ZERO
      end if

      rixmin = ( minval(geom%region(0)%px(:)) - offset(1) ) / reggrid%wx
      rizmin = ( minval(geom%region(0)%pz(:)) - offset(2) ) / reggrid%wz
      rixmax = ( maxval(geom%region(0)%px(:)) - offset(1) ) / reggrid%wx
      rizmax = ( maxval(geom%region(0)%pz(:)) - offset(2) ) / reggrid%wz

      ixmin1 = ceiling( rixmin - HALF )
      if( geom%freex(1) )  ixmin1 = ixmin1 - 1
      ixmin2 = ceiling( rixmin - sn%pmax*(ONE+TWO*EPSMB)/reggrid%wx )
      ixmin = min( ixmin1, ixmin2 )
      grid%nfrozenx(1) = ixmin1 - ixmin
      if( .not. geom%freex(1) )  grid%nfrozenx(1) = grid%nfrozenx(1) + 1

      ixmax1 = ceiling( rixmax + HALF )
      if( geom%freex(2) )  ixmax1 = ixmax1 + 1
      ixmax2 = ceiling( rixmax + sn%pmax*(ONE+TWO*EPSMB)/reggrid%wx )
      ixmax = max( ixmax1, ixmax2 )
      grid%nfrozenx(2) = ixmax - ixmax1
      if( .not. geom%freex(2) )  grid%nfrozenx(2) = grid%nfrozenx(2) + 1

      izmin1 = ceiling( rizmin - HALF )
      if( geom%freez(1) )  izmin1 = izmin1 - 1
      izmin2 = ceiling( rizmin - sn%pmax*(ONE+TWO*EPSMB)/reggrid%wz )
      izmin = min( izmin1, izmin2 )
      grid%nfrozenz(1) = izmin1 - izmin
      if( .not. geom%freez(1) )  grid%nfrozenz(1) = grid%nfrozenz(1) + 1

      izmax1 = ceiling( rizmax + HALF )
      if( geom%freez(2) )  izmax1 = izmax1 + 1
      izmax2 = ceiling( rizmax + sn%pmax*(ONE+TWO*EPSMB)/reggrid%wz )
      izmax = max( izmax1, izmax2 )
      grid%nfrozenz(2) = izmax - izmax1
      if( .not. geom%freez(2) )  grid%nfrozenz(2) = grid%nfrozenz(2) + 1
!old1
!      if( geom%freex(1) ) then                  ! margin consistent with movable points
!         rixmin = rixmin - max( ONE+HALF, sn%pmax/reggrid%wx )
!      else
!         rixmin = rixmin - HALF
!      end if
!      if( geom%freez(1) ) then
!         rizmin = rizmin - max( ONE+HALF, sn%pmax/reggrid%wz )
!      else
!         rizmin = rizmin - HALF
!      end if
!      if( geom%freex(2) ) then
!        rixmax = rixmax + max( ONE+HALF, sn%pmax/reggrid%wx )
!      else
!         rixmax = rixmax + HALF
!      end if
!      if( geom%freez(2) ) then
!         rizmax = rizmax + max( ONE+HALF, sn%pmax/reggrid%wz )
!      else
!         rizmax = rizmax + HALF
!      end if
!     ixmin = ceiling( rixmin )
!      izmin = ceiling( rizmin )
!      ixmax = ceiling( rixmax )
!      izmax = ceiling( rizmax )
!old2
!      ixmin = ceiling( ( minval( geom%region(0)%px(:) ) - offset(1) ) / reggrid%wx )
!      izmin = ceiling( ( minval( geom%region(0)%pz(:) ) - offset(2) ) / reggrid%wz )
!      ixmax = ceiling( ( maxval( geom%region(0)%px(:) ) - offset(1) ) / reggrid%wx )
!      izmax = ceiling( ( maxval( geom%region(0)%pz(:) ) - offset(2) ) / reggrid%wz )
!!                                               ! safety border of pmax, min 1 cell
!      ixmin = ixmin - ceiling( sn%pmax / reggrid%wx )
!      izmin = izmin - ceiling( sn%pmax / reggrid%wz )
!      ixmax = ixmax + ceiling( sn%pmax / reggrid%wx )
!      izmax = izmax + ceiling( sn%pmax / reggrid%wz )
!
      reggrid%nx = ixmax - ixmin + 1
      reggrid%nz = izmax - izmin + 1
!
!  position of grid lines
!
      allocate( reggrid%x( 0:reggrid%nx ) )
      allocate( reggrid%z( 0:reggrid%nz ) )
      reggrid%x = [ (offset(1) + ix * reggrid%wx, ix=ixmin-1,ixmax) ]
      reggrid%z = [ (offset(2) + iz * reggrid%wz, iz=izmin-1,izmax) ]
!
!  Initialize the surface
!
      allocate( xc( 0:reggrid%nx+1 ), stat=istat, errmsg=errmsg )
      call check_allocation( 'initCurviGrid2d', 'xc', istat, errmsg )

      allocate( zc( 0:reggrid%nz+1 ), stat=istat, errmsg=errmsg )
      call check_allocation( 'initCurviGrid2d', 'zc', istat, errmsg )
                                                !  distance grid at midpoints
      xc(1:reggrid%nx) = HALF * ( reggrid%x(0:reggrid%nx-1) + reggrid%x(1:reggrid%nx) )
      zc(1:reggrid%nz) = HALF * ( reggrid%z(0:reggrid%nz-1) + reggrid%z(1:reggrid%nz) )
      xc(0) = TWO * xc(1) - xc(2)
      zc(0) = TWO * zc(1) - zc(2)
      xc(reggrid%nx+1) = TWO * xc(reggrid%nx) - xc(reggrid%nx-1)
      zc(reggrid%nz+1) = TWO * zc(reggrid%nz) - zc(reggrid%nz-1)

      np = geom%region(0)%np + 1
      call initSurface2dSDFfromPolygons( xc, zc, [ np ],  &
!          ============================  i   i      i
                                         geom%region(0)%px(1:np), geom%region(0)%pz(1:np) )
!                                                  i                        i
      deallocate( xc, zc )
!
!  Initializations that are also done at the end of readCurviGrid2d
!
      call initOnceCurviGrid2d( )
!          ===================
!
!  Initializations that are also done at the end of relaxCurviGrid2d
!
      call initDerivedCurviGrid2d( )
!          ======================

   end subroutine initCurviGrid2d
!
!===================================================================================
!
   subroutine getCellCurviGrid2d( x, z, ic, outside_, lx_, lz_ )
!             ==================  i  i  o      o       o    o
!
!     Determine cell the point pos is located in
!
!  Arguments:
!
      real(REAL8), intent(in) :: x              ! x-coordinate of position to be checked
      real(REAL8), intent(in) :: z              ! z-coordinate of position to be checked
      integer(INT4), intent(out) :: ic          ! index of cell (region)
      logical, optional, intent(out) :: outside_! indicates whether pos is outside all cells
      logical, optional, intent(out) :: lx_     ! flag indicating point is on "horizontal"
                                                !    grid line iz
      logical, optional, intent(out) :: lz_     ! flag indicating point is on "vertical"
                                                !    grid line ix
!
!  Local Variables:
!
      logical :: inside, onedge(4)
      integer(INT4) :: ix, iz, iovlp, novlp, ixx, izz, ixv(4), izv(4), iv
      real(REAL8) :: zz, xx, xv(4), zv(4)
!
!  Start of Executable Code:
!  ------------------------
!
!  Determine regular cell pos is located in (positions outside grid map to nearest cell)
!
      if( present( outside_ ) ) then
         outside_ = x < reggrid%x(0) .or. x > reggrid%x(reggrid%nx) .or.  &
                    z < reggrid%z(0) .or. z > reggrid%z(reggrid%nx)
      end if

      xx = max( reggrid%xcmin, min( reggrid%xcmax, x ) ) ! don't use min/max with the indices
      zz = max( reggrid%zcmin, min( reggrid%zcmax, z ) ) !   to avoid integer overflow

      ix = ceiling( ( xx - reggrid%x(0) ) / reggrid%wx )
      iz = ceiling( ( zz - reggrid%z(0) ) / reggrid%wz )
!
!  Determine among the overlapped curvilinear cells the one pos is actually in
!  (this code should be moved to getCellIndexGrid2d)
!
      novlp = size( reggrid%cell(ix,iz)%icovlp )
      if( novlp == 1 ) then
         ic = reggrid%cell(ix,iz)%icovlp(1)
      else
         do iovlp = 1, novlp
            ic = reggrid%cell(ix,iz)%icovlp(iovlp)
            ixx = ixfun(ic)
            izz = izfun(ic)
            ixv(:) = [ixx-1,ixx,ixx,ixx-1]
            izv(:) = [izz-1,izz-1,izz,izz]
            do iv = 1, 4
               xv(iv) = grid%x(ixv(iv),izv(iv))
               zv(iv) = grid%z(ixv(iv),izv(iv))
            end do

            call relPosPolygon( xx, zz, xv(:), zv(:), inside, onedge )
!                =============  i   i   i      i        o       o
            if( inside .or. onedge(2) .or. onedge(3) ) exit
         end do
      end if

      if( present( lx_ ) ) then
         lx_ = onedge(3)
         lz_ = onedge(2)
      end if

   end subroutine getCellCurviGrid2d
!
!==========================================================================================
!
   subroutine getOtherNearCellsCurviGrid2d( proj, posc, ffp, pmax, ncn, icnear )
!             ============================   i     i     i    i     o     o
!
!     Determine cells that are within the free flight path (ffp) cylinder (= cylinder with
!     height ffp and radius pmax).
!     Only cells in the range (1:nx,1:nz) are returned.
!
!  Arguments:
!
      type(Proj_type), intent(in) :: proj       ! projectile
      real(REAL8), intent(in) :: posc(3)        ! position of end of free flight path
      real(REAL8), intent(in) :: ffp            ! free flight path (= height of cylinder)
      real(REAL8), intent(in) :: pmax           ! maximum impact parameter (= radius of cyl.)
      integer(INT4), intent(out) :: ncn         ! number of regions found
      integer(INT4), allocatable, intent(inout) :: icnear(:)   ! index of regions found
!
!  Local Variables:
!
      character(256) :: errmsg
      integer(INT4) :: ixmin, izmin, ixmax, izmax
      integer(INT4) :: ix, iz, ixmin_ab, izmin_ab, ixmax_ab, izmax_ab
      integer(INT4) :: ic, ixx, izz, iovlp, novlp, istat
      real(REAL8) :: z, dz, z1, z2, zmin, zmax, zisct
      real(REAL8) :: x, dx, x1, x2, xmin, xmax, xisct, ffp_p
      real(REAL8) :: ax, az, bx, bz
      real(REAL8) :: rlam1, rlam2, dirx, diry, dirz, dirxy, dirxz, diryz, dirx_p, dirz_p
      logical, allocatable :: regoverlap(:,:)   ! 2D mask indicating regular cells overlaping
                                                !    with the ffp cylinder
      logical, allocatable :: overlap(:,:)      ! same for main cells

      logical :: isInside, cut1, cut2
!
!  Start of Executable Code:
!  ------------------------
!
      x = proj%pos(1)
      z = proj%pos(3)
      dirx = proj%dir(1)
      diry = proj%dir(2)
      dirz = proj%dir(3)
      dirxy = sqrt( dirx**2 + diry**2 )
      dirxz = sqrt( dirx**2 + dirz**2 )
      diryz = sqrt( diry**2 + dirz**2 )
!
!  Determine coordinates of the points A and B (outmost points of ffp cylinder
!  on the projection of the direction of motion)
!
      if( abs(dirxz) < EPSMB ) then
         dirx_p = ZERO
         dirz_p = ONE                           ! arbitrarily set direction parallel to z
      else
         dirx_p = dirx / dirxz
         dirz_p = dirz / dirxz
      end if

      ax = x - dirx_p * pmax * abs(diry)
      az = z - dirz_p * pmax * abs(diry)
      bx = posc(1) + dirx_p * pmax * abs(diry)
      bz = posc(3) + dirz_p * pmax * abs(diry)
!
!  Define the rectangle that contains the projection projection of the ffp cylinder and
!  the corresponding range of cell indices
!
      dx = pmax * diryz
      dz = pmax * dirxy

      x1 = x - sign( dx, dirx )
      x2 = posc(1) + sign( dx, dirx )
      z1 = z - sign( dz, dirz )
      z2 = posc(3) + sign( dz, dirz )

      xmin = min( x1, x2 )
      xmax = max( x1, x2 )
      zmin = min( z1, z2 )
      zmax = max( z1, z2 )

      xmin = max( reggrid%xcmin, min( reggrid%xcmax, xmin ) )     ! don't apply bounds
      xmax = max( reggrid%xcmin, min( reggrid%xcmax, xmax ) )     !    to ixmin,ixmax
      zmin = max( reggrid%zcmin, min( reggrid%zcmax, zmin ) )     !    because of possible
      zmax = max( reggrid%zcmin, min( reggrid%zcmax, zmax ) )     !    integer overflow

      ixmin = ceiling( ( xmin - reggrid%x(0) ) / reggrid%wx )  ! indices are valid
      ixmax = ceiling( ( xmax - reggrid%x(0) ) / reggrid%wx )  !    cell indices in
      izmin = ceiling( ( zmin - reggrid%z(0) ) / reggrid%wz )  !    the range
      izmax = ceiling( ( zmax - reggrid%z(0) ) / reggrid%wz )  !    [1..nxz(:)]
!
!  Allocate mask
!
      allocate( regoverlap( ixmin:ixmax, izmin:izmax ) )
      regoverlap(:,:) = .false.
!
!  Find the range of cell indices that contain the line A-B
!
      xmin = min( ax, bx )
      xmax = max( ax, bx )
      zmin = min( az, bz )
      zmax = max( az, bz )

      xmin = max( reggrid%xcmin, min( reggrid%xcmax, xmin ) )
      xmax = max( reggrid%xcmin, min( reggrid%xcmax, xmax ) )
      zmin = max( reggrid%zcmin, min( reggrid%zcmax, zmin ) )
      zmax = max( reggrid%zcmin, min( reggrid%zcmax, zmax ) )

      ixmin_ab = ceiling( ( xmin - reggrid%x(0) ) / reggrid%wx )  ! indices must be in the
      ixmax_ab = ceiling( ( xmax - reggrid%x(0) ) / reggrid%wx )  !    range [ixmin(:)...
      izmin_ab = ceiling( ( zmin - reggrid%z(0) ) / reggrid%wz )  !    ixmax(:)]
      izmax_ab = ceiling( ( zmax - reggrid%z(0) ) / reggrid%wz )  !
!
!  Check for cell corners inside the ffp cylinder. All four adjoining cells overlap the
!  ffp cylinder
!
      ffp_p = ffp * dirxz                          ! p = projected to xz-plane

      do iz = izmin, izmax-1
         do ix = ixmin, ixmax-1

            isInside = .false.

            x1 = reggrid%x(ix)
            z1 = reggrid%z(iz)

            isInside = isInsideFFPcyl( [x1,z1], [x,z], [dirx_p,dirz_p], ffp_p, pmax,  &
!                      ==============      i         i             i            i     i
                                       abs(diry) )
!                                          i

            if( isInside )  regoverlap(ix:ix+1,iz:iz+1) = .true.

         end do
      end do
!
!  It is possible that the ffp cylinder passes between all cell corners. Therefore
!  check in addition for cell edges intersecting the extended ffp (line AB)
!
!  Check intersection with horizontal (z=const) lines
!
      do iz = izmin_ab, izmax_ab-1

         x1 = reggrid%x(0)
         x2 = reggrid%x(reggrid%nx)
         z1 = reggrid%z(iz)
         z2 = reggrid%z(iz)

         call mcisct( rlam1, rlam2, cut1, cut2, x1, z1, x2, z2, ax, az, bx, bz )
!             ======    o      o     o     o    i   i   i   i   i   i   i   i

         if( .not. cut2 ) cycle ! no valid intersection only if A^B outside the polygon
                                ! (can this happen?)

         xisct = x1 + ( x2 - x1 ) * rlam1       ! x-coordinate of intersection point

         ix = ceiling( ( xisct - reggrid%x(0) ) / reggrid%wx )
         ix = max( ixmin, min( ixmax, ix ) )    ! just in case AB intersects at
                                                !    cell area boundary
         regoverlap(ix,iz:iz+1) = .true.

      end do
!
!  Check intersection with the vertical (x=const) lines
!
      do ix = ixmin_ab, ixmax_ab-1

         x1 = reggrid%x(ix)
         x2 = reggrid%x(ix)
         z1 = reggrid%z(0)
         z2 = reggrid%z(reggrid%nz)

         call mcisct( rlam1, rlam2, cut1, cut2, x1, z1, x2, z2, ax, az, bx, bz )
!             ======    o      o     o     o    i   i   i   i   i   i   i   i

         if( .not. cut2 )  cycle ! (?)

         zisct = z1 + ( z2 - z1 ) * rlam1       ! z-coordinate of intersection point

         iz = ceiling( ( zisct - reggrid%z(0) ) / reggrid%wz )
         iz = max( izmin, min( izmax, iz ) )    ! just in case AB intersects at
                                                !    cell area boundary
         regoverlap(ix:ix+1,iz) = .true.

      end do
!
!  Use icovlp array to determine (curvilinear) cells that possibly overlap with the ffp
!  cylinder
!
      allocate( overlap( ixmin-1:ixmax+1, izmin-1:izmax+1 ), stat=istat, errmsg=errmsg )
      call check_allocation( 'getOtherNearCellsCurviGrid2d', 'overlap', istat, errmsg )
      overlap(:,:) = .false.

      do iz = izmin, izmax
         do ix = ixmin, ixmax
            if( regoverlap(ix,iz) ) then
               novlp = size( reggrid%cell(ix,iz)%icovlp )
               do iovlp = 1, novlp
                  ic = reggrid%cell(ix,iz)%icovlp(iovlp)
                  if( ic == icfun(ix,iz) ) cycle   ! do not include the cell itself
                  ixx = ixfun(ic)
                  izz = izfun(ic)
                  if( ixx < ixmin-1 .or. ixx > ixmax+1 .or.  &
                      izz < izmin-1 .or. izz > izmax+1 ) then
                     print *,  'getOtherNearCellsCurviGrid2d: ixx/izz out of range'
                     cycle
                  end if
                  overlap(ixx,izz) = .true.
               end do
            end if
         end do
      end do
!
!  Convert the logical overlap array into a list of cell indices
!
      deallocate( icnear )
      allocate( icnear( (ixmax-ixmin+3)*(izmax-izmin+3) ) ) ! maximum possible size

      ncn = 0

      do iz = izmin-1, izmax+1
         do ix = ixmin-1, ixmax+1
            if( .not.overlap(ix,iz) )  cycle
            ncn = ncn + 1
            icnear(ncn) = icfun(ix,iz)
         end do
      end do
!
!  the version below works, but is slower
!
!      call newStack( icnearStack )
!!          ========       o
!!
!      do iz = izmin, izmax
!         do ix = ixmin, ixmax
!            if( overlap(ix,iz) )  call pushStack( icfun(ix,iz), icnearStack )
!!                                      =========       i            io
!         end do
!      end do
!!
!      deallocate( icnear )
!      call stack2Array( icnearStack, icnear )
!!          ===========      io         o
!      ncn = size( icnear )
!
      deallocate( overlap )
!
   end subroutine getOtherNearCellsCurviGrid2d
!
!=======================================================================
!
   subroutine getCellsCurviGrid2d( cell )
!             ===================   o
!
!     Get the cell infos
!
!     NOTE: We have to do it the slow way, since reshape does not seem to work with
!           allocatable arrays containing allocatable arrays
!
!  Arguments:
!
      type(Cell_type), allocatable, intent(out) :: cell(:)  ! cells
!
!  Local Variables:
!
      character(80) :: text
      character(256) :: errmsg
      integer(INT4) :: ncell, istat, ix, iz, ic
!
!  Start of Executable Code:
!  ------------------------
!
      ncell = grid%nx * grid%nz                 !

      allocate( cell( ncell ), stat=istat, errmsg=errmsg )
      call check_allocation( 'getCellsRegGrid2d', 'cell', istat, errmsg )

      do iz = 1, grid%nz
         do ix = 1, grid%nx                  ! contrary to its name, allocate also copies
            ic = icfun(ix,iz)                !    source to destination

            allocate( cell(ic)%xv, source=grid%cell(ix,iz)%xv, stat=istat, errmsg=errmsg )
            write( text, '(a,i0,a)' )  'cell(', ic, ')%xv'
            call check_allocation( 'getCellsCurviGrid2d', text, istat, errmsg )

            allocate( cell(ic)%zv, source=grid%cell(ix,iz)%zv, stat=istat, errmsg=errmsg )
            write( text, '(a,i0,a)' )  'cell(', ic, ')%zv'
            call check_allocation( 'getCellsCurviGrid2d', text, istat, errmsg )

            allocate( cell(ic)%icnear, source=grid%cell(ix,iz)%icnear, stat=istat, &
                      errmsg=errmsg )
            write( text, '(a,i0,a)' )  'cell(', ic, ')%icnear'
            call check_allocation( 'getCellsCurviGrid2d', text, istat, errmsg )

            allocate( cell(ic)%prob, source=grid%cell(ix,iz)%prob, stat=istat, errmsg=errmsg )
            write( text, '(a,i0,a)' )  'cell(', ic, ')%prob'
            call check_allocation( 'getCellsCurviGrid2d', text, istat, errmsg )

         end do
      end do

   end subroutine getCellsCurviGrid2d
!
!=======================================================================
!
   subroutine getOvlpAreasCurviGrid2d( px, pz, area, ic_area )
!             =======================  i   i    o       o
!
!     Obtain the fractional cell areas of a curvilinear grid inside the given polygon.
!
!     The polygon may extend outside [0,nx] a,d [0,nz], but corresponding areas are ignored.
!
!     px, pz must contain the repeated first point as the last point.
!
!  Arguments:
!
      real(REAL8), intent(in) :: px(:)                   ! x-coordinates of polygon corners
      real(REAL8), intent(in) :: pz(:)                   ! z-coordinates of polygon corners
      real(REAL8), allocatable, intent(out) :: area(:)   ! area inside polygon
                                                         !    first index runs in x-direction
                                                         !    second index runs in z-direction
      integer(INT4), allocatable, intent(out) :: ic_area(:) ! cell indices corresponding to the
                                                            !    areas
!
!  Local Data Types:
!
      type :: Point_type
!             ==========
         real(REAL8) :: x                    ! x-coordinate
         real(REAL8) :: z                    ! z-coordinate
         integer(INT4) :: ix                 ! cell index in x-direction
         integer(INT4) :: iz                 ! cell index in z-direction
         logical :: lx                       ! flag indicating point is on a "horizontal" line
         logical :: lz                       ! flag indicating point is on a "vertical" line
      end type

      type :: Line_type
!             =========
         type(Point_type) :: p1                 ! start point
         type(Point_type) :: p2                 ! end point
      end type
!
!  Local Variables:
!
      logical :: lend, inside, onedge(4), found, lx, lz
      character(256) :: msg
      integer(INT4) :: i, ip, np, ix0, iz0, ix, iz, ii, ni, ia, na, iarea, narea, istat
      integer(INT4) :: izmin, izmax
      integer(INT4), allocatable :: ixmin(:), ixmax(:)
      real(REAL8) :: polyarea, area_scal, coseta, sineta
      real(REAL8) :: xv(4), zv(4), cellarea, x, z
      real(REAL8), allocatable :: area2d(:,:)
      real(REAL8), allocatable :: pxx(:), pzz(:)

      type(Point_type), allocatable :: pp(:)    ! polygon points in scaled coordinates
      type(Point_type), allocatable :: pa(:)    ! polygon for area calculation
      type(Point_type), allocatable :: pi(:)    ! intersection points

      integer(INT4), allocatable :: isort(:)

      type(Line_type) :: line
!
!  Start of Executable Code:
!  ------------------------
!
      allocate( ixmin( 0:grid%nz+1 ), stat=istat, errmsg=msg )
      call check_allocation( 'getOvlpAreasCurviGrid2d', 'ixmin', istat, msg )
      allocate( ixmax( 0:grid%nz+1 ), stat=istat, errmsg=msg )
      call check_allocation( 'getOvlpAreasCurviGrid2d', 'ixmax', istat, msg )
      allocate( area2d( grid%nx, grid%nz ), stat=istat, errmsg=msg )
      call check_allocation( 'getOvlpAreasCurviGrid2d', 'area2d', istat, msg )
      allocate( pp( size( px ) ), stat=istat, errmsg=msg )
      call check_allocation( 'getOvlpAreasCurviGrid2d', 'pp', istat, msg )
      allocate( pa( size( pp )+5 ), stat=istat, errmsg=msg )
      call check_allocation( 'getOvlpAreasCurviGrid2d', 'pa', istat, msg )
      allocate( pi( grid%nx+grid%nz ), stat=istat, errmsg=msg )
      call check_allocation( 'getOvlpAreasCurviGrid2d', 'pi', istat, msg )
      allocate( isort( size( pp ) ), stat=istat, errmsg=msg )
      call check_allocation( 'getOvlpAreasCurviGrid2d', 'isort', istat, msg )
!
!  Convert px, pz into polygon points
!
      np = size( px )

      do ip = 1, np
         pp(ip)%x = px(ip)                      ! set coordinates
         pp(ip)%z = pz(ip)
                                                ! search for cell indices in concentric squares
         ix0 = ( pp(ip)%x - reggrid%x(0) ) / reggrid%wx
         iz0 = ( pp(ip)%z - reggrid%z(0) ) / reggrid%wz
         found = .false.
         do i = 0, max( ix0, grid%nx-ix0, iz0, grid%nz-iz0 )
            do iz = iz0-i, iz0+i
               if( iz < 1 .or. iz > grid%nz )  cycle
               do ix = ix0-i, ix0+i
                  if( ix < 1 .or. ix > grid%nx )  cycle
                  if( abs( iz-iz0 ) < i .and. abs( ix-ix0 ) < i )  cycle

                  xv = [ grid%x(ix-1,iz-1), grid%x(ix,iz-1), grid%x(ix,iz), grid%x(ix-1,iz) ]
                  zv = [ grid%z(ix-1,iz-1), grid%z(ix,iz-1), grid%z(ix,iz), grid%z(ix-1,iz) ]
                  call relPosPolygon( pp(ip)%x, pp(ip)%z, xv(:), zv(:), inside, onedge )
!                      =============     i         i        i      i      o       o
                  found = inside .or. any( onedge(:) )
                  if( found ) exit
               end do
               if( found ) exit
            end do
            if( found ) exit
         end do

         if( .not. found ) then
            print *,  'ix0,iz0=', ix0, iz0
            print *,  'ip,x,z=', ip, pp(ip)%x, pp(ip)%z
            do iz = iz0-1, iz0+1
               if( iz < 0 .or. iz > grid%nz )  cycle
               do ix = ix0-1, ix0+1
                  if( ix < 0 .or. ix > grid%nx )  cycle
                  print *,  'ix,iz,x,z=',ix,iz,grid%x(ix,iz),grid%z(ix,iz)
               end do
            end do
            error stop 'getOvlpAreasCurviGrid2d: point outside grid'
         end if

         if( onedge(1) ) then
            iz = iz - 1
            onedge(1) = .false.
            onedge(3) = .true.
         end if
         if( onedge(4) ) then
            ix = ix - 1
            onedge(4) = .false.
            onedge(2) = .true.
         end if
         pp(ip)%ix = ix                         ! set cell indices
         pp(ip)%iz = iz
         pp(ip)%lx = onedge(3)                  ! set edge flags
         pp(ip)%lz = onedge(2)
      end do
!
!  Initialize the arrays to be used to collect the areas
!
      ixmin(:) = grid%nx + 2
      ixmax(:) = -1
      izmin = grid%nz + 2
      izmax = -1
      area2d(:,:) = ZERO
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

         call getIntersectionPoints( line, ni, pi )
!             =====================   i    o   o

         if( ni /= 0 ) exit

      end do
!
!  If no intersection point has been found, polygon is completely inside a cell
!
      if( ni == 0 ) then

         ix = pp(1)%ix
         iz = pp(1)%iz
         lx = pp(1)%lx
         lz = pp(1)%lz

         if( ( ix >= 1 .and. ix <= grid%nx .and. iz >= 1 .and. iz <= grid%nz ) .or.  &
             ( ix == 0 .and. lz ) .or. ( iz == 0 .and. lx ) ) then
            allocate( area( 1 ), stat=istat, errmsg=msg )
            call check_allocation( 'getOvlpAreasCurviGrid2d', 'area', istat, msg )
            area(1) = polyarea
            allocate( ic_area( 1 ), stat=istat, errmsg=msg )
            call check_allocation( 'getOvlpAreasCurviGrid2d', 'ic_area', istat, msg )
            ic_area(1) = icfun(max(ix,1),max(iz,1))
         else
            allocate( area( 0 ), stat=istat, errmsg=msg )
            call check_allocation( 'getOvlpAreasCurviGrid2d', 'area', istat, msg )
            allocate( ic_area( 0 ), stat=istat, errmsg=msg )
            call check_allocation( 'getOvlpAreasCurviGrid2d', 'ic_area', istat, msg )
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

               call getIntersectionPoints( line, ni, pi )
!                   =====================   i    o   o

               if( ni /= 0 )  exit              ! exit when intersection point has been found

            end do

            ii = 1                              ! reset intersection point index

         else

            ii = ii + 1                         ! iterate intersection point index
            if( ii > size( pi ) ) then
               do ii = 1, size( pi )
                  print *, 'ii,pi=', ii, pi(ii)%x, pi(ii)%z, pi(ii)%ix, pi(ii)%iz,  &
                           pi(ii)%lx, pi(ii)%lz
               end do
               ix0 = pi(size( pi ))%ix
               iz0 = pi(size( pi ))%iz
               do iz = iz0-1, iz0+1
                  if( iz < 0 .or. iz > grid%nz )  cycle
                  do ix = ix0-1, ix0+1
                     if( ix < 0 .or. ix > grid%nx )  cycle
                     write( *, * )  'ix,iz,grid%x,grid%z=',ix,iz,grid%x(ix,iz),grid%z(ix,iz)
                  end do
               end do
               error stop 'getOvlpAreasCurviGrid2d: array pi too small'
            end if

         end if
!
!  Assign intersection point to area polygon
!
         ia = ia + 1
         pa(ia) = pi(ii)
!
!  Calculate cell index and polygon area
!
         na = ia
         call closePolygon( ix, iz, na, pa )
!             ============  o   o   io  io

         if( ix >= 1 .and. ix <= grid%nx .and. iz >= 1 .and. iz <= grid%nz ) then

            polyarea = getPolyArea( pa(:)%x, pa(:)%z, na )
!                      ===========     i        i     i
                                                ! AREA SHOULD NOT BE NEGATIVE IF getPolyArea
                                                ! AND closePolygon are consistent
            if( polyarea < ZERO )  print *, 'ix,iz,polyarea=',ix,iz,polyarea
            area2d(ix,iz) = area2d(ix,iz) + polyarea
         end if

         izmin = min( izmin, iz )
         izmax = max( izmax, iz )
         ixmin(iz) = min( ixmin(iz), ix )
         ixmax(iz) = max( ixmax(iz), ix )

         if( lend )  exit

      end do
!
!  End of loop over cells
!  ----------------------
!
!  Clip ixmin, ixmax (izmin, izmax)
!
      izmin = max( 1, min( grid%nz, izmin ) )
      izmax = max( 1, min( grid%nz, izmax ) )
      ixmin(:) = max( 1, min( grid%nx, ixmin(:) ) )
      ixmax(:) = max( 1, min( grid%nx, ixmax(:) ) )
!
!  Finalize area calculation: subtract multiples of the cell area as needed
!
      pxx = [ px(:), px(1) ]
      pzz = [ pz(:), pz(1) ]
      call mcauxl( coseta, sineta, pxx, pzz )
!          ======    o       o     i   i

      narea = 0
      do iz = izmin, izmax
         do ix = ixmin(iz), ixmax(iz)
            xv(:) = [ grid%x(ix-1,iz-1), grid%x(ix,iz-1), grid%x(ix,iz), grid%x(ix-1,iz) ]
            zv(:) = [ grid%z(ix-1,iz-1), grid%z(ix,iz-1), grid%z(ix,iz), grid%z(ix-1,iz) ]
            cellarea = getPolyArea( xv(:), zv(:), 4 )
!                      ===========   i      i     i
            !
            !  Subtract integer multiple of cell area such that area is in [0,cell_area)
            !  Note: this is essential and works (not sure about round-off errors!?!)
            !
            area_scal = area2d(ix,iz) / cellarea
            area_scal = area_scal - floor(area_scal)

            if( area_scal < EPSMB .or. area_scal > ONE-EPSMB ) then
                                                ! check whether cell "center" is inside polygon
               call pointInsidePolygon( xv(:), zv(:), x, z )
!                   ==================   i      i     o  o
               if( mcloc( x, z, pxx, pzz, coseta, sineta ) > 0 ) then  ! cell inside polygon
!                  =====  i  i  i   i     i       i
                  area2d(ix,iz) = cellarea
                  narea = narea + 1
               else                                                    ! cell outside polygon
                  area2d(ix,iz) = ZERO
               endif
            else
               area2d(ix,iz) = cellarea * area_scal
               narea = narea + 1
            end if
         end do
      end do
!
!  Convert 2d area array into 1d array
!
      allocate( area( narea ), stat=istat, errmsg=msg )
      call check_allocation( 'getOvlpAreasCurviGrid2d', 'area', istat, msg )
      allocate( ic_area( narea ), stat=istat, errmsg=msg )
      call check_allocation( 'getOvlpAreasCurviGrid2d', 'ic_area', istat, msg )

      iarea = 0
      do iz = izmin, izmax
         do ix = ixmin(iz), ixmax(iz)
            if( area2d(ix,iz) > ZERO ) then
               iarea = iarea + 1
               area(iarea) = area2d(ix,iz)
               ic_area(iarea) = icfun(ix,iz)
            end if
         end do
      end do
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

      end subroutine getLine
!
!-----------------------------------------------------------------------------------------
!
      subroutine getIntersectionPoints( line, ni, pi )
!                =====================   i    o   o
!
!        Determine the intersection points of line with the grid lines.
!        This is done iteratively starting from line%p1.
!
         type(Line_type), intent(in) :: line
         integer(INT4), intent(out) :: ni
         type(Point_type), intent(out) :: pi(:)

         logical :: cut1, cut2, found
         logical, allocatable :: mask(:)
         integer(INT4) :: ix, iz, ii, i, imin
         integer(INT4), allocatable :: ix1(:), ix2(:), iz1(:), iz2(:)
         integer(INT4) :: isort(size( pi ))     ! index for ordering according to di
         real(REAL8) :: rlam, rmu, rlammin, rmumin
         real(REAL8) :: di(size( pi ))          ! distance from start point of line
         real(REAL8) :: px(9), pz(9)
         type(Point_type) :: p1                 ! starting point of remaining line

         ii = 0

         p1 = line%p1

         do
!
!  Determine the possibly intersecting grid lines
!
            ix = p1%ix
            iz = p1%iz
            if( p1%lx .and. p1%lz ) then        ! starting point is a grid point
               ix1 = [ ix+1, ix+1, ix  , ix-1, ix-1, ix-1, ix  , ix+1 ]
               iz1 = [ iz  , iz+1, iz+1, iz+1, iz  , iz-1, iz-1, iz-1 ]
            else if( p1%lx ) then               ! starting point is on a horizontal grid line
               ix1 = [ ix, ix  , ix-1, ix-1, ix-1, ix   ]
               iz1 = [ iz, iz+1, iz+1, iz  , iz-1, iz-1 ]
            else if( p1%lz ) then               ! starting point is on a vertical grid line
               ix1 = [ ix, ix+1, ix+1, ix  , ix-1, ix-1 ]
               iz1 = [ iz, iz  , iz-1, iz-1, iz-1, iz   ]
            else                                ! starting point is inside a grid cell
               ix1 = [ ix  , ix, ix-1, ix-1 ]
               iz1 = [ iz-1, iz, iz  , iz-1 ]
            end if

            ix2 = cshift( ix1(:), -1 )          ! end point indices
            iz2 = cshift( iz1(:), -1 )
                                                ! exclude non-existing grid lines
            mask = ix1(:) >= 0 .and. ix1(:) <= grid%nx .and.  &
                   ix2(:) >= 0 .and. ix2(:) <= grid%nx .and.  &
                   iz1(:) >= 0 .and. iz1(:) <= grid%nz .and.  &
                   iz2(:) >= 0 .and. iz2(:) <= grid%nz
            ix1 = pack( ix1, mask )
            ix2 = pack( ix2, mask )
            iz1 = pack( iz1, mask )
            iz2 = pack( iz2, mask )
!
!  Determine intersection point of line with grid lines
!
            rlammin = BIG
            found = .false.
            do i = 1, size( ix1 )
               call mcisct( rlam, rmu, cut1, cut2, p1%x, p1%z, line%p2%x, line%p2%z,  &
!                   ======   o     o    o     o     i     i        i          i
                            grid%x(ix1(i),iz1(i)), grid%z(ix1(i),iz1(i)),  &
!                                     i                      i
                            grid%x(ix2(i),iz2(i)), grid%z(ix2(i),iz2(i)), EPSMB )
!                                  i                      i                 i
               if( cut2 .and. rlam < rlammin) then
                  rlammin = rlam
                  rmumin = rmu
                  imin = i
                  found = .true.
               end if
            end do

            if( .not. found )  exit
!
!  Save intersection point
!
            ii = ii + 1
            if( ii > size( pi ) ) then
               do ii = 1, size( pi )
                  print *, 'ii,pi=', ii, pi(ii)%x, pi(ii)%z, pi(ii)%ix, pi(ii)%iz,  &
                           pi(ii)%lx, pi(ii)%lz
               end do
               ix0 = pi(size( pi ))%ix
               iz0 = pi(size( pi ))%iz
               do iz = iz0-1, iz0+1
                  if( iz < 0 .or. iz > grid%nz )  cycle
                  do ix = ix0-1, ix0+1
                     if( ix < 0 .or. ix > grid%nx )  cycle
                     print *,  'ix,iz,grid%x,grid%z=',ix,iz,grid%x(ix,iz),grid%z(ix,iz)
                  end do
               end do
               error stop 'getIntersectionPoints: pi too small'
            end if
            i = imin
            rmu = rmumin
            if( abs(rmu) < EPSMB ) then       ! intersection point is a grid point
               pi(ii)%x = grid%x(ix1(i),iz1(i))
               pi(ii)%z = grid%z(ix1(i),iz1(i))
               pi(ii)%ix = ix1(i)
               pi(ii)%iz = iz1(i)
               pi(ii)%lx = .true.
               pi(ii)%lz = .true.
            else if( rmu > ONE-EPSMB ) then     ! intersection point is a grid point
               pi(ii)%x = grid%x(ix2(i),iz2(i))
               pi(ii)%z = grid%z(ix2(i),iz2(i))
               pi(ii)%ix = ix2(i)
               pi(ii)%iz = iz2(i)
               pi(ii)%lx = .true.
               pi(ii)%lz = .true.
            else                          ! intersection point is on a grid line
               pi(ii)%x = (ONE-rmu) * grid%x(ix1(i),iz1(i)) + rmu * grid%x(ix2(i),iz2(i))
               pi(ii)%z = (ONE-rmu) * grid%z(ix1(i),iz1(i)) + rmu * grid%z(ix2(i),iz2(i))
               pi(ii)%ix = max( ix1(i), ix2(i) )
               pi(ii)%iz = max( iz1(i), iz2(i) )
               pi(ii)%lx = ( iz1(i) == iz2(i) )
               pi(ii)%lz = ( ix1(i) == ix2(i) )
            end if
!
!  Prepare next intersection
!
            p1 = pi(ii)

         end do

         ni = ii

      end subroutine getIntersectionPoints
!
!-----------------------------------------------------------------------------------------
!
      subroutine closePolygon( ix, iz, na, pa )
!                ============  o   o   io  io
!
!        Close the polygon counter-clockwise along the edges of the cell
!
         integer(INT4), intent(out) :: ix
         integer(INT4), intent(out) :: iz
         integer(INT4), intent(inout) :: na
         type(Point_type), intent(inout) :: pa(:)

         integer(INT4) :: ie, ie1, ie2
         real(REAL8) :: xv(4), zv(4)
!
!  We number the cell edges and corners like this:
!
!                            4 |   4   | 3
!                         iz --+-------+--      z
!                              |       |       ^
!                             1|(ix,iz)|3      |
!                              |       |       |
!                       iz-1 --+-------+--     +---> x
!                            1 |   2   | 2
!                             ix-1     ix
!
!  The cell (ix,iz) is between the grid lines ix-1 and ix, and iz-1 and iz.
!  Edge 1 belongs to cell (ix-1,iz), edge 2 to (ix,iz-1), and edges 3 and 4 to (ix,iz).
!  The corners belong to the cell with the same indices.
!
!  Determine the cell index the area polygon belongs to
!
         ix = maxval(pa(1:na)%ix)
         iz = maxval(pa(1:na)%iz)

         if( ix == 0 .or. iz == 0 ) then        ! cell does not exist
            na = 0
            return
         end if
!
!  The corner points
!
         xv(:) = [ grid%x(ix-1,iz-1), grid%x(ix,iz-1), grid%x(ix,iz), grid%x(ix-1,iz) ]
         zv(:) = [ grid%z(ix-1,iz-1), grid%z(ix,iz-1), grid%z(ix,iz), grid%z(ix-1,iz) ]

!
!  The cell edges the first and last point are on
!
         if( pa(1)%lz ) then
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

         if( pa(na)%lz ) then
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
            if( (pa(1)%x-xv(ie1))**2 + (pa(1)%z-zv(ie1))**2 <=  &
                (pa(na)%x-xv(ie1))**2 + (pa(na)%z-zv(ie1))**2 )  return
         end if
!
!  Starting with edge ie2, add corner points until we reach edge ie1
!
         ie = ie2

         do

            na = na + 1

            pa(na)%x = xv(ie)
            pa(na)%z = zv(ie)

            ie = modulo( ie, 4 ) + 1

            if( ie == ie1 )  exit

         end do

      end subroutine closePolygon
!
!-----------------------------------------------------------------------------------------
!
   end subroutine getOvlpAreasCurviGrid2d
!
!=========================================================================================
!
   subroutine relaxCurviGrid2d( ff0, ff, dpa, stress, ovlp )
!             ================   i   i    i     io     o
!
!     Redistribute the cell densities so that stress relaxes. This is done in two steps:
!     First, depending on dam%relaxmodel, the geometry of the grid is adjusted so that the cell
!     densities (ratio of cell content and cell size) are identical to the equilibrium densities
!     for the respective target composition, or they are relaxed using a viscoelastic model.
!     Then the cell contents are mapped back to an equidistant grid.
!
!     This subroutine should not be called too often in order to avoid artificial diffusion,
!     but often enough to maintain realistic densities.
!
!  Arguments:
!
      type :: Ovlp_type
!             =========
         sequence
         integer(INT4), allocatable :: ic(:)    ! cell indices of overlapped cells
         real(REAL8), allocatable :: area(:)    ! overlap area
         real(REAL8) :: emptyarea               ! total overlap area with empty cells
      end type Ovlp_type

      real(REAL8), intent(in) :: ff0(:)         ! (ncell) fill fractions before BC time step
      real(REAL8), intent(in) :: ff(:)          ! (ncell) fill fractions after BC time step if
                                                !    there were no relaxation
      real(REAL8), intent(in) :: dpa(:)         ! (ncell) displacements per atom during BC step
      real(REAL8), intent(inout) :: stress(:,:) ! (13,ncell) pressure and deviatoric stress
                                                !    tensor elements
      type(Ovlp_type), intent(out) :: ovlp(:)   ! (ncell) overlap information
!
!  Local Variables:
!
      logical, allocatable :: mask(:,:)
      logical, allocatable :: movable(:,:)
      character(256) :: msg
      integer(INT4) :: ix, iz, ic, ixx, izz
      integer(INT4) :: isurf, istat, iovlp, novlp, ixo, izo, ico, in, io
      real(REAL8) :: x4(5), z4(5)
      real(REAL8), allocatable :: area0(:,:), area(:,:), xc(:), zc(:)
      real(REAL8), allocatable :: xsurf(:), zsurf(:)
      real(REAL8), allocatable :: ff0_2(:,:), ff_2(:,:)
      real(REAL8), allocatable :: dpa_2(:,:), stress_2(:,:,:)

      integer(INT4), allocatable :: icovlp(:)
      real(REAL8), allocatable :: areaovlp(:)
!
!  Start of Executable Code:
!  ------------------------
!
      allocate( mask( reggrid%nx, reggrid%nz ), stat=istat, errmsg=msg )
      call check_allocation( 'relaxCurviGrid2d', 'mask', istat, msg )
      allocate( movable( 0:reggrid%nx, 0:reggrid%nz ), stat=istat, errmsg=msg )
      call check_allocation( 'relaxCurviGrid2d', 'movable', istat, msg )
      allocate( area0( grid%nx, grid%nz ), stat=istat, errmsg=msg )
      call check_allocation( 'relaxCurviGrid2d', 'area0', istat, msg )
      allocate( area( grid%nx, grid%nz ), stat=istat, errmsg=msg )
      call check_allocation( 'relaxCurviGrid2d', 'area', istat, msg )
      allocate( xc( 0:reggrid%nx+1 ), stat=istat, errmsg=msg )
      call check_allocation( 'relaxCurviGrid2d', 'xc', istat, msg )
      allocate( zc( 0:reggrid%nz+1 ), stat=istat, errmsg=msg )
      call check_allocation( 'relaxCurviGrid2d', 'zc', istat, msg )
      allocate( xsurf( size( surface%ix(:) ) ), stat=istat, errmsg=msg )
      call check_allocation( 'relaxCurviGrid2d', 'xsurf', istat, msg )
      allocate( zsurf( size( surface%ix(:) ) ), stat=istat, errmsg=msg )
      call check_allocation( 'relaxCurviGrid2d', 'zsurf', istat, msg )
      allocate( ff0_2( reggrid%nx, reggrid%nz ), stat=istat, errmsg=msg )
      call check_allocation( 'relaxCurviGrid2d', 'ff0_2', istat, msg )
      allocate( ff_2( reggrid%nx, reggrid%nz ), stat=istat, errmsg=msg )
      call check_allocation( 'relaxCurviGrid2d', 'ff_2', istat, msg )
      allocate( dpa_2( reggrid%nx, reggrid%nz ), stat=istat, errmsg=msg )
      call check_allocation( 'relaxCurviGrid2d', 'dpa_2', istat, msg )
      allocate( stress_2( 13, reggrid%nx, reggrid%nz ), stat=istat, errmsg=msg )
      call check_allocation( 'relaxCurviGrid2d', 'stress_2', istat, msg )
!
!  Calculate relaxed cell sizes and densities
!
      do iz = 1, grid%nz
         do ix = 1, grid%nx
            area0(ix,iz) = getPolyArea( grid%cell(ix,iz)%xv(:), grid%cell(ix,iz)%zv(:), 5 )
!                          ===========           i                       i              i
         end do
      end do
      area(:,:) = area0(:,:) * reshape( ff, [ grid%nx, grid%nz ] )
!!
!      call bndConnected( mask )            ! empty cells connected to boundary
!!          ============   io
!      mask(:,:) = .not.mask(:,:)       ! full cells and unconnected empty cells
!
!  Set cells to be relaxed
!
      mask(:,:) = .false.                   ! possibly do not relax the outermost cells
      do iz = 1+grid%nfrozenz(1), grid%nz-grid%nfrozenz(2)
         do ix = 1+grid%nfrozenx(1), grid%nx-grid%nfrozenx(2)
            mask(ix,iz) = ( grid%cell(ix,iz)%icnear(1) == icfun(ix,iz) )
         end do
      end do

      movable(:,:) = .false.                    ! do not move the outermost grid points
      do iz = 1+grid%nfrozenz(1), grid%nz-1-grid%nfrozenz(2)
         do ix = 1+grid%nfrozenx(1), grid%nx-1-grid%nfrozenx(2)
            do izz = max( iz, 1 ), min( iz+1, grid%nz )
               do ixx = max( ix, 1 ), min( ix+1, grid%nx)
                  !if( mask(ixx,izz) .and. dpa(icfun(ixx,izz)) > dam%dparelax )  &
                     movable(ix,iz) = .true.
               end do
            end do
         end do
      end do

      print *, count(movable(:,:)), 'movable points'
!
!  Convert ff0, ff, dpa to 2d arrays
!
      ff0_2 = reshape(ff0, [ grid%nx, grid%nz ] )
      ff_2 = reshape(ff, [ grid%nx, grid%nz ] )
      dpa_2 = reshape(dpa, [ grid%nx, grid%nz ] )
      stress_2 = reshape( stress, [ 13, grid%nx, grid%nz ] )
!
!  Copy grid to relaxedgrid
!
      relaxedgrid%nx = grid%nx
      relaxedgrid%nz = grid%nz

      if( .not. allocated( relaxedgrid%x ) ) then
         allocate( relaxedgrid%x( relaxedgrid%nx, relaxedgrid%nz ), stat=istat, errmsg=msg )
         call check_allocation( 'relaxCurviGrid2d', 'relaxedgrid%x', istat, msg )
         allocate( relaxedgrid%z( relaxedgrid%nx, relaxedgrid%nz ), stat=istat, errmsg=msg )
         call check_allocation( 'relaxCurviGrid2d', 'relaxedgrid%z', istat, msg )
         !allocate( grid%cell( relaxedgrid%nx, relaxedgrid%nz ), stat=istat, errmsg=msg )
         !                                                        if( istat /= 0 ) goto 900
      end if

      relaxedgrid%x(:,:) = grid%x(:,:)
      relaxedgrid%z(:,:) = grid%z(:,:)          ! We don't need the other relaxedgrid variables
!
!  Relax
!
      if( dam%relaxmodel == ENUM%relaxmodel%none ) then
         area(:,:) = area0(:,:)
      else if( dam%relaxmodel == ENUM%relaxmodel%full ) then
#ifdef AMAT
         call relaxAreas2d_simple( relaxedgrid%x, relaxedgrid%z, area, mask, movable )
!             ===================       io             io         io     i      i
#else
         error stop 'relaxCurviGrid2d: relaxAreas2d_simple not available'
#endif
      else if( dam%relaxmodel == ENUM%relaxmodel%viscoelastic ) then
#ifdef AMAT
         call relaxViscoElastic2d( relaxedgrid%x, relaxedgrid%z, mask, movable, ff0_2, ff_2,  &
!             ===================       io             io         i       i       i     i
                                   dpa_2, stress_2 )
!                                    i      io
#else
         error stop 'relaxViscoElastic2d: relaxViscoElastic2d not available'
#endif
         stress = reshape( stress_2, [ 13, grid%nx*grid%nz ] )
      else
         error stop 'relaxCurviGrid2d: Invalid dam%relaxmodel'
      end if
!
!  Initialize the surface
!                                               !  distance grid at midpoints
      xc(1:reggrid%nx) = HALF * ( reggrid%x(0:reggrid%nx-1) + reggrid%x(1:reggrid%nx) )
      zc(1:reggrid%nz) = HALF * ( reggrid%z(0:reggrid%nz-1) + reggrid%z(1:reggrid%nz) )
      xc(0) = TWO * xc(1) - xc(2)
      zc(0) = TWO * zc(1) - zc(2)
      xc(reggrid%nx+1) = TWO * xc(reggrid%nx) - xc(reggrid%nx-1)
      zc(reggrid%nz+1) = TWO * zc(reggrid%nz) - zc(reggrid%nz-1)

      do isurf = 1, size( surface%ix(:) )
         xsurf(isurf) = relaxedgrid%x(surface%ix(isurf),surface%iz(isurf))
         zsurf(isurf) = relaxedgrid%z(surface%ix(isurf),surface%iz(isurf))
      end do
      call initSurface2dSDFfromPolygons( xc, zc, surface%ipend(:), xsurf(:), zsurf(:) )
!          ============================  i   i          i             i         i
!
!  Initialize derived quantities
!
      call initDerivedCurviGrid2d( )
!          ======================
!
!  Determine overlap areas with original grid   !!! can grid move out of original grid?
!
      do iz = 1, relaxedgrid%nz                 ! loop over relaxed cells
         do ix = 1, relaxedgrid%nx
            ic = icfun(ix,iz)

            ovlp(ic)%emptyarea = ZERO

            x4(1:4) = [ relaxedgrid%x(ix-1,iz-1), relaxedgrid%x(ix,iz-1),  &
                         relaxedgrid%x(ix,iz), relaxedgrid%x(ix-1,iz) ]
            z4(1:4) = [ relaxedgrid%z(ix-1,iz-1), relaxedgrid%z(ix,iz-1),  &
                         relaxedgrid%z(ix,iz), relaxedgrid%z(ix-1,iz) ]
            x4(5) = x4(1)
            z4(5) = z4(1)

            call getOvlpAreasCurviGrid2d( x4(:), z4(:), areaovlp, icovlp )
!                =======================    i      i       o        o

            if( .true. ) then                   ! do not redistribute overlap area with empty
                                                !    cells (does not stricly conserve mass)
               novlp = 0
               do io = 1, size( icovlp(:) )
                  ico = icovlp(io)
                  ixo = ixfun( ico )
                  izo = izfun( ico )
                  if( grid%cell(ixo,izo)%icnear(1) == ico )  novlp = novlp + 1
               end do

               allocate( ovlp(ic)%area( novlp ), stat=istat, errmsg=msg )
               call check_allocation( 'relaxCurviGrid2d', 'ovlp%area', istat, msg )
               allocate( ovlp(ic)%ic( novlp ), stat=istat, errmsg=msg )
               call check_allocation( 'relaxCurviGrid2d', 'ovlp%ic', istat, msg )

               iovlp = 0
               do io = 1, size( icovlp(:) )
                  ico = icovlp(io)
                  ixo = ixfun( ico )
                  izo = izfun( ico )
                  if( grid%cell(ixo,izo)%icnear(1) == ico ) then
                     iovlp = iovlp + 1
                     ovlp(ic)%ic(iovlp) = ico
                     ovlp(ic)%area(iovlp) = areaovlp(io)
                  else
                     ovlp(ic)%emptyarea = ovlp(ic)%emptyarea + areaovlp(io)
                  end if
               end do

            else                                ! redistribute overlap area with empty cells
                                                !    to near cells (conserves mass but causes
                                                !    problems during relaxation)

               novlp = 0                        ! for empty overlapped cells we need to divide
               do io = 1, size( icovlp(:) )     !    overlap area among icnear
                  ixo = ixfun( icovlp(io) )     ! non-empty cells are included since icnear=ic
                  izo = izfun( icovlp(io) )
                  novlp = novlp + size( grid%cell(ixo,izo)%icnear(:) )
               end do

               allocate( ovlp(ic)%area( novlp ), stat=istat, errmsg=msg )
               call check_allocation( 'relaxCurviGrid2d', 'ovlp%area', istat, msg )
               allocate( ovlp(ic)%ic( novlp ), stat=istat, errmsg=msg )
               call check_allocation( 'relaxCurviGrid2d', 'ovlp%ic', istat, msg )

               iovlp = 0
               do io = 1, size( icovlp(:) )
                  ixo = ixfun( icovlp(io) )
                  izo = izfun( icovlp(io) )
                  do in = 1, size( grid%cell(ixo,izo)%icnear(:) )
                     iovlp = iovlp + 1
                     ovlp(ic)%ic(iovlp) = grid%cell(ixo,izo)%icnear(in)
                     ovlp(ic)%area(iovlp) = areaovlp(io) * grid%cell(ixo,izo)%prob(in)
                  end do
               end do

            end if

         end do
      end do

      return

   end subroutine relaxCurviGrid2d
!
!=======================================================================
!
   subroutine writeCurviGrid2d( nout, relaxed_ )
!             ================   i       i
!
!     Write the grid to the file with unit number nout
!
!  Arguments:
!
      integer(INT4), intent(in) :: nout         ! unit number
      logical, optional, intent(in) :: relaxed_ ! indicates relaxed grid to be written
!
!  Local Variables:
!
      integer(INT4) :: ix, iz
!
!  Start of Executable Code:
!  ------------------------
!
      if( present( relaxed_ ) ) then
         write( nout, '(a)' )  '# Grid: curvilinear'
         write( nout, '(2i15)' )  relaxedgrid%nx, relaxedgrid%nz
         do iz = 0, relaxedgrid%nz
            do ix = 0, relaxedgrid%nx
               write( nout, '(2f15.6)' )  relaxedgrid%x(ix,iz), relaxedgrid%z(ix,iz)
            end do
         end do
      else
         write( nout, '(a)' )  '# Grid: curvilinear'
         write( nout, '(2i15)' )  grid%nx, grid%nz
         do iz = 0, grid%nz
            do ix = 0, grid%nx
               write( nout, '(2f15.6)' )  grid%x(ix,iz), grid%z(ix,iz)
            end do
         end do
      end if

   end subroutine writeCurviGrid2d
!
!=======================================================================
!
   subroutine readCurviGrid2d( nin, valid )
!             ===============   i     o
!
!     Read the grid from the CELL2 file
!
!     UNTESTED CODE !
!
!  Arguments:
!
      integer(INT4), intent(in) :: nin          ! unit number of CELL2 file
      logical, intent(out) :: valid             ! indicates valid grid has been read
!
!  Local Variables:
!
      character(80) :: line
      integer(INT4) :: ix, iz
      real(REAL8) :: xdummy, zdummy
!
!  Start of Executable Code:
!  ------------------------
!
!  Read grid info
!
      write( nin, '(a)' )  line
      if( line == '# Grid: curvilinear' ) then
         read( nin, '(2i15)' )  grid%nx, grid%nz
         do iz = 0, grid%nz
            do ix = 0, grid%nx
               read( nin, *, err=900, end=900 )  xdummy, zdummy
            end do
         end do
         valid = .false.
      else if( line == '# Grid: regular' ) then
         read( nin, *, err=900, end=900 )  reggrid%nx, reggrid%nz
         if( allocated( reggrid%x ) )  deallocate( reggrid%x, reggrid%z )
         allocate( reggrid%x( 0:reggrid%nx ) )
         allocate( reggrid%z( 0:reggrid%nz ) )
         read( nin, *, err=900, end=900 )  reggrid%x(0), reggrid%z(0)
         read( nin, *, err=900, end=900 )  reggrid%wx, reggrid%wz
         reggrid%x(1:reggrid%nx) = [ (reggrid%x(0) + ix * reggrid%wx, ix = 1, reggrid%nx) ]
         reggrid%z(1:reggrid%nz) = [ (reggrid%z(0) + iz * reggrid%wz, iz = 1, reggrid%nz) ]
         valid = .true.
      else
         goto 901
      end if
!
!  Shouldn't initOnceCurviGrid2d be called here?
!
      return

900   error stop 'readCurviGrid2d: error during reading CELL2 file'
901   error stop 'readCurviGrid2d: unknown grid type'

   end subroutine readCurviGrid2d
!
!=======================================================================
!
   subroutine initOnceCurviGrid2d( )
!             ===================
!
!     Further initializations of the grid in addition to initCurviGrid2d
!
!  Arguments:
!
!     None
!
!  Local Variables:
!
      character(80) :: text
      character(256) :: errmsg
      integer(INT4) :: istat, ix, iz
      integer(INT4) :: iv, ixv(5), izv(5)
!
!  Start of Executable Code:
!  ------------------------
!
!  Minimum and maximum cell center positions
!
      reggrid%xcmin = HALF * ( reggrid%x(0) + reggrid%x(1) )
      reggrid%xcmax = HALF * ( reggrid%x(reggrid%nx-1) + reggrid%x(reggrid%nx) )
      reggrid%zcmin = HALF * ( reggrid%z(0) + reggrid%z(1) )
      reggrid%zcmax = HALF * ( reggrid%z(reggrid%nz-1) + reggrid%z(reggrid%nz) )
!
!  Allocate main grid
!
      grid%nx = reggrid%nx
      grid%nz = reggrid%nz

      allocate( grid%x( 0:grid%nx, 0:grid%nz ), stat=istat, errmsg=errmsg )
      call check_allocation( 'initOnceCurviGrid2d', 'grid%x', istat, errmsg )

      allocate( grid%z( 0:grid%nx, 0:grid%nz ), stat=istat, errmsg=errmsg )
      call check_allocation( 'initOnceCurviGrid2d', 'grid%z', istat, errmsg )

      grid%x(:,:) = spread( reggrid%x(:), dim=2, ncopies=grid%nz+1 )
      grid%z(:,:) = spread( reggrid%z(:), dim=1, ncopies=grid%nx+1 )
!
!  Allocate cells and initialize cell vertices
!
      allocate( reggrid%cell( reggrid%nx, reggrid%nz ), stat=istat, errmsg=errmsg )
      call check_allocation( 'initOnceCurviGrid2d', 'reggrid%cell', istat, errmsg )

      allocate( grid%cell( grid%nx, grid%nz ), stat=istat, errmsg=errmsg )
      call check_allocation( 'initOnceCurviGrid2d', 'grid%cell', istat, errmsg )

      do iz = 1, grid%nz
         do ix = 1, grid%nx

            allocate( grid%cell(ix,iz)%xv( 5 ), stat=istat, errmsg=errmsg )
            write( text, '(a,i0,a,i0,a)' )  'grid%cell(', ix, ',', iz, ')%xv'
            call check_allocation( 'initOnceCurviGrid2d', text, istat, errmsg )

            allocate( grid%cell(ix,iz)%zv( 5 ), stat=istat, errmsg=errmsg )
            write( text, '(a,i0,a,i0,a)' )  'grid%cell(', ix, ',', iz, ')%zv'
            call check_allocation( 'initOnceCurviGrid2d', text, istat, errmsg )

            ixv(:) = [ix-1,ix,ix,ix-1,ix-1]
            izv(:) = [iz-1,iz-1,iz,iz,iz-1]
            do iv = 1, 5
               grid%cell(ix,iz)%xv(iv) = grid%x( ixv(iv), izv(iv) )
               grid%cell(ix,iz)%zv(iv) = grid%z( ixv(iv), izv(iv) )
            end do

            ! dummy allocations
            allocate( grid%cell(ix,iz)%icnear( 1 ), stat=istat, errmsg=errmsg )
            write( text, '(a,i0,a,i0,a)' )  'grid%cell(', ix, ',', iz, ')%icnear'
            call check_allocation( 'initOnceCurviGrid2d', text, istat, errmsg )

            allocate( grid%cell(ix,iz)%prob( 1 ), stat=istat, errmsg=errmsg )
            write( text, '(a,i0,a,i0,a)' )  'grid%cell(', ix, ',', iz, ')%prob'
            call check_allocation( 'initOnceCurviGrid2d', text, istat, errmsg )
         end do
      end do
!
!  Curvilinear grid set to a Cartesian/regular grid
!
      if( .not. allocated( relaxedgrid%x ) ) then
         relaxedgrid%nx = reggrid%nx
         relaxedgrid%nz = reggrid%nz

         allocate( relaxedgrid%x( 0:relaxedgrid%nx, 0:relaxedgrid%nz ), stat=istat,  &
                   errmsg=errmsg )
         call check_allocation( 'initOnceCurviGrid2d', 'relaxedgrid%x', istat, errmsg )

         allocate( relaxedgrid%z( 0:relaxedgrid%nx, 0:relaxedgrid%nz ), stat=istat,  &
                   errmsg=errmsg )
         call check_allocation( 'initOnceCurviGrid2d', 'relaxedgrid%z', istat, errmsg )

      end if

      relaxedgrid%x(:,:) = spread( reggrid%x(:), dim=2, ncopies=relaxedgrid%nz+1 )
      relaxedgrid%z(:,:) = spread( reggrid%z(:), dim=1, ncopies=relaxedgrid%nx+1 )

   end subroutine initOnceCurviGrid2d
!
!=======================================================================
!
   subroutine initDerivedCurviGrid2d( )
!             ======================
!
!     Initialize derived quantities of the cell structure of the 2d target
!
!  Arguments:
!
!     None
!
!  Local Variables:
!
      character(256) :: msg
      integer(INT4) :: ix, iz, ixx, izz, i, n, ncn, istat, inear, nnear, stride, ic, icc, ip
      integer(INT4) :: ix0, ix1, ix2, iz0, iz1, iz2, ixv(5), izv(5), iv, iv1, iv2
!      real(REAL8) :: distmin2, distmin2_1, distmin2_2, distmax2, distmax2_1, distmax2_2
      real(REAL8) :: distmin2, distmin2min, distmax2, distmax2min, diameter, diameter2, ddist
      real(REAL8) :: dir(3), dirx, dirz
      real(REAL8) :: xv(4), zv(4)

      real(REAL8), allocatable :: xc(:), zc(:), ff2(:,:)

      logical, allocatable :: empty(:,:)
      real(REAL8), allocatable :: dist(:,:)
      integer(INT4), allocatable :: ic_area(:), isort(:)
      real(REAL8), allocatable :: area(:)
      real(REAL8), allocatable :: distmin2_tmp(:)

      integer(INT4), allocatable :: ipbeg(:)    ! first point index of connected component
      integer(INT4), allocatable :: ipend(:)    ! last point index of connected component
      real(REAL8), allocatable :: px(:)         ! x-coordinates of polygon corners
      real(REAL8), allocatable :: pz(:)         ! z-coordinates of polygon corners

      type :: NearCell_type                     ! cell with distmin2 < minval(distmax2)
!             =============
         integer(INT4) :: ix                    ! cell index in x-direction
         integer(INT4) :: iz                    ! cell index in z-direction
         real(REAL8) :: distmin2                ! square of minimum minimum distance of
                                                !    point inside the cell from boundary of
                                                !    region ir
         real(REAL8) :: distmax2                ! square of maximum minimum distance of
                                                !    point inside the cell from boundary of
                                                !    region ir
      end type NearCell_type

      type(Nearcell_type), allocatable :: nearcell(:)
!
!  Start of Executable Code:
!  ------------------------
!
      allocate( xc( 0:reggrid%nx+1 ), stat=istat, errmsg=msg )
      call check_allocation( 'initDerivedCurviGrid2d', 'xc', istat, msg )
      allocate( zc( 0:reggrid%nz+1 ), stat=istat, errmsg=msg )
      call check_allocation( 'initDerivedCurviGrid2d', 'zc', istat, msg )
      allocate( ff2( 0:reggrid%nx+1, 0:reggrid%nz+1 ), stat=istat, errmsg=msg )
      call check_allocation( 'initDerivedCurviGrid2d', 'ff2', istat, msg )
      allocate( empty( grid%nx, grid%nz ), stat=istat, errmsg=msg )
      call check_allocation( 'initDerivedCurviGrid2d', 'empty', istat, msg )
      allocate( dist( 0:grid%nx, 0:grid%nz ), stat=istat, errmsg=msg )
      call check_allocation( 'initDerivedCurviGrid2d', 'dist', istat, msg )
      allocate( nearcell( reggrid%nx*reggrid%nz ), stat=istat, errmsg=msg )
      call check_allocation( 'initDerivedCurviGrid2d', 'nearcell', istat, msg )
!
!  Set up a new main grid
!
      grid%x(:,:) = spread( reggrid%x(:), dim=2, ncopies=grid%nz+1 )
      grid%z(:,:) = spread( reggrid%z(:), dim=1, ncopies=grid%nx+1 )
!
!  Modify grid near the surface so there are nodes on the surface
!
!  initialize the distance array
!
      do iz = 0, grid%nz
         do ix = 0, grid%nx
            call getDistSurface2d( [grid%x(ix,iz), ZERO, grid%z(ix,iz)], dist(ix,iz) )
!                ================                   i                         o
            if( dist(ix,iz) == -BIG )  error stop 'initCurviGrid: node outside grid'
         end do
      end do
!
!  get the surface polygon in Cartesian grid coordinates (each vertex either has an
!  integer+0.5 x- or integer+0.5 z-coordinate)
!
      call getCartSurface2dSDFPolygons( ipend, px, pz )   ! px, pz relative to xc
!          ===========================    o    o   o
      allocate( ipbeg, source=ipend, stat=istat, errmsg=msg )
      call check_allocation( 'initDerivedCurviGrid2d', 'ipbeg', istat, msg )
      ipbeg(:) = cshift( ipbeg(:) + 1, -1 )
      ipbeg(1) = 1
      px(:) = px(:) - HALF                            ! shift so px are relative to grid%x(0)
      pz(:) = pz(:) - HALF                            ! shift so pz are relative to grid%z(0)
!
!  move nodes near the boundary to the boundary and remember their indices (they are the
!  surface nodes)
!
      if( allocated( surface%ix ) )  deallocate( surface%ix, surface%iz, surface%ipend )
      allocate( surface%ix( size( px ) ), stat=istat, errmsg=msg )
      call check_allocation( 'initDerivedCurviGrid2d', 'surface%ix', istat, msg )
      allocate( surface%iz( size( px ) ), stat=istat, errmsg=msg )
      call check_allocation( 'initDerivedCurviGrid2d', 'surface%iz', istat, msg )
      allocate( surface%ipend, source=ipend, stat=istat, errmsg=msg )
      call check_allocation( 'initDerivedCurviGrid2d', 'surface%ipend', istat, msg )

      i = 1
      do ip = 1, size( px )
         if( ip == ipend(i) ) then              ! ipend(last_index)=size( px )!
            surface%ix(ip) = surface%ix(ipbeg(i))
            surface%iz(ip) = surface%iz(ipbeg(i))
            i = i + 1
            cycle
         end if
         ix = nint( HALF*(px(ip)+px(ip+1)) )
         iz = nint( HALF*(pz(ip)+pz(ip+1)) )
         call getDistSurface2d( [grid%x(ix,iz), ZERO, grid%z(ix,iz)], dist(ix,iz), dir_=dir )
!             ================                   i                        o           o
         dirx = dir(1)
         dirz = dir(3)
!         grid%x(ix,iz) = grid%x(ix,iz) + dist(ix,iz) * dirx
!         grid%z(ix,iz) = grid%z(ix,iz) + dist(ix,iz) * dirz
         grid%x(ix,iz) = reggrid%x(0) + HALF*(px(ip)+px(ip+1)) * reggrid%wx
         grid%z(ix,iz) = reggrid%z(0) + HALF*(pz(ip)+pz(ip+1)) * reggrid%wz
         surface%ix(ip) = ix
         surface%iz(ip) = iz
      end do
!
!  update the grid cell vertices
!
      do iz = 1, grid%nz
         do ix = 1, grid%nx
            ixv(:) = [ix-1,ix,ix,ix-1,ix-1]
            izv(:) = [iz-1,iz-1,iz,iz,iz-1]
            do iv = 1, 5
               grid%cell(ix,iz)%xv(iv) = grid%x( ixv(iv), izv(iv) )
               grid%cell(ix,iz)%zv(iv) = grid%z( ixv(iv), izv(iv) )
            end do
         end do
      end do
!
!  Determine for each reggrid cell the grid cells that overlap the grid cell
!
      do iz = 1, reggrid%nz
         do ix = 1, reggrid%nx
            xv(:) = [ reggrid%x(ix-1), reggrid%x(ix), reggrid%x(ix), reggrid%x(ix-1) ]
            zv(:) = [ reggrid%z(iz-1), reggrid%z(iz-1), reggrid%z(iz), reggrid%z(iz) ]
            call getOvlpAreasCurviGrid2d( xv, zv, area, ic_area )
!                =======================  i   i    o      o
                                                ! sort ic_area according to decreasing area
            n = size( area )
            isort = [ (i, i=1,n) ]
            call quickSort( area(:), isort(:), 1, n )
!                =========    io       io      i  i
            reggrid%cell(ix,iz)%icovlp = ic_area( isort(n:1:-1) )
         end do
      end do
!
!  Determine non-empty cells that possibly are closest to points inside each grid cell
!
!  mark as empty, cells that have all four corners with negative distance from the surface,
!  or that are left of the surface polygon (surface is clockwise)
!
      empty(:,:) = dist(0:grid%nx-1,0:grid%nz-1) < 0 .and. dist(0:grid%nx-1,1:grid%nz) < 0  &
                   .and. dist(1:grid%nx,0:grid%nz-1) < 0 .and. dist(1:grid%nx,1:grid%nz) < 0

      i = 1
      ix1 = surface%ix(ipend(1)-1)
      iz1 = surface%iz(ipend(1)-1)
      ix2 = surface%ix(1)
      iz2 = surface%iz(1)
      do ip = 2, size( surface%ix )
         if( ip == ipend(i)+1 ) then
            i = i + 1
            ix1 = surface%ix(ipend(i)-1)
            iz1 = surface%iz(ipend(i)-1)
            ix2 = surface%ix(ip)
            iz2 = surface%iz(ip)
            cycle
         end if
         ix0 = ix1
         iz0 = iz1
         ix1 = ix2
         iz1 = iz2
         ix2 = surface%ix(ip)
         iz2 = surface%iz(ip)
         if( ix2 == ix1 .and. iz2 == iz1 + 1 ) then
            empty(ix1,iz1+1) = .true.
            if( ix0 == ix1 + 1 .and. iz0 == iz1 )  empty(ix1,iz1) = .true.
         else if( ix2 == ix1 + 1 .and. iz2 == iz1 ) then
            empty(ix1+1,iz1+1) = .true.
            if( ix0 == ix1 .and. iz0 == iz1 - 1 )  empty(ix1,iz1+1) = .true.
         else if( ix2 == ix1 .and. iz2 == iz1 - 1 ) then
            empty(ix1+1,iz1) = .true.
            if( ix0 == ix1 - 1 .and. iz0 == iz1 )  empty(ix1+1,iz1+1) = .true.
         else if( ix2 == ix1 - 1 .and. iz2 == iz1 ) then
            empty(ix1,iz1) = .true.
            if( ix0 == ix1 .and. iz0 == iz1 + 1 )  empty(ix1+1,iz1) = .true.
         else
            print *,  'Surface points:'
            do i = 1, size( surface%ix )
               print *,  'x,z=', grid%x(surface%ix(i),surface%iz(i)),  &
                                 grid%z(surface%ix(i),surface%iz(i))
            end do
            print *,  'initCurviGrid2d: Surface points are not neighbors:'
            print *,  'Surface point 1 indices:', ix1, iz1
            print *,  'Surface point 2 indices:', ix2, iz2
            error stop
         end if
      end do

      deallocate( ipbeg, ipend, px, pz )
!
!  Loop over cells
!  ---------------
!
      do iz = 1, grid%nz
         do ix = 1, grid%nx
            ic = icfun(ix,iz)

            deallocate( grid%cell(ix,iz)%icnear )
            deallocate( grid%cell(ix,iz)%prob )

            if( .not. empty(ix,iz) ) then       ! The single near cell of a non-empty cell is
                                                !    the cell itself
               allocate( grid%cell(ix,iz)%icnear( 1 ), stat=istat, errmsg=msg )
               call check_allocation( 'initDerivedCurviGrid2d', 'grid%cell%icnear', istat, msg )
               grid%cell(ix,iz)%icnear(1) = ic
               allocate( grid%cell(ix,iz)%prob( 1 ), stat=istat, errmsg=msg )
               call check_allocation( 'initDerivedCurviGrid2d', 'grid%cell%prob', istat, msg )
               grid%cell(ix,iz)%prob(1) = ONE
               cycle
            end if

            nearcell(:)%distmin2 = BIG          ! distmin2, distmax2 from the current cell
            nearcell(:)%distmax2 = BIG          !    (ix,iz)
!            distmin2_1 = BIG                    ! smallest distmin2
!            distmin2_2 = BIG                    ! second-smallest distmin2
!            distmax2_1 = BIG                    ! distmax2 corresponding to smallest distmin2
!            distmax2_2 = BIG                    ! distmax2 corresponding to second-smallest
                                                 !    distmin2
            distmax2min = BIG

            ncn = 0
            !
            !  loop over shells (should always exit via exit)
            !
            do n = 1, max( grid%nx, grid%nz )
!               if( (n-1)*min(reggrid%wx,reggrid%wz) > distmax2_2 )  exit
               distmin2min = BIG
               !
               !  loop over cells in shell
               !
               do ixx = ix-n, ix+n              ! skip cells outside the grid
                  if( ixx < 1 .or. ixx > grid%nx )  cycle
                  if( ixx == ix-n .or. ixx == ix+n )  then
                     stride = 1
                  else
                     stride = 2*n
                  end if
                  do izz = iz-n, iz+n, stride   ! skip cells outside the grid
                     if( izz < 1 .or. izz > grid%nz )  cycle
                     icc = icfun(ixx,izz)
                                                ! skip empty cells
                     if( empty(ixx,izz) )  cycle

                     call getDistPolyPoly( grid%cell(ix,iz)%xv(:), grid%cell(ix,iz)%zv(:),  &
!                         ===============              i                       i
                                           grid%cell(ixx,izz)%xv(:), grid%cell(ixx,izz)%zv(:), &
!                                                      i                          i
                                           distmin2, distmax2 )
!                                             o         o
                     !
                     !  maintain distances from second-nearest cell
                     !
!                     if( distmin2 < distmin2_1 ) then
!                        distmin2_2 = distmin2_1
!                        distmin2_1 = distmin2
!                        distmax2_2 = distmax2_1
!                        distmax2_1 = distmax2
!                     else if( distmin2 < distmin2_2 ) then
!                        distmin2_2 = distmin2
!                        distmax2_2 = distmax2
!                     end if

                     distmax2min = min( distmax2min, distmax2 )   ! overall minimum of distmax2
                     distmin2min = min( distmin2min, distmax2 )   ! shell minimum of distmin2

                     ncn = ncn + 1

                     nearcell(ncn)%ix = ixx
                     nearcell(ncn)%iz = izz
                     nearcell(ncn)%distmin2 = distmin2
                     nearcell(ncn)%distmax2 = distmax2

                  end do
               end do
               !
               !  end of loop over cells in shell
               !
               if( distmin2min == BIG )  cycle           ! only empty cells in this shell
               if( distmin2min >= distmax2min )  exit    ! no more relevant cells

            end do
            !
            !  Keep only near cells that obey distmin2 < distmax2_2 and sort according to
            !  distmin2
            !
!            nnear = count( nearcell(1:ncn)%distmin2 <= distmax2_2 )
!            nearcell(1:nnear) = pack( nearcell(1:ncn), nearcell(1:ncn)%distmin2 <= distmax2_2 )
            nnear = count( nearcell(1:ncn)%distmin2 <= distmax2min )
            nearcell(1:nnear) = pack( nearcell(1:ncn), nearcell(1:ncn)%distmin2 <= distmax2min )

            isort = [ (i, i=1,nnear) ]
            distmin2_tmp = nearcell(1:nnear)%distmin2
            call quickSort( distmin2_tmp(:), isort(:), 1, nnear )
!                =========        io            io     i    i

            allocate( grid%cell(ix,iz)%icnear( nnear ), stat=istat, errmsg=msg )
            call check_allocation( 'initDerivedCurviGrid2d', 'grid%cell%icnear', istat, msg )
            allocate( grid%cell(ix,iz)%prob( nnear ), stat=istat, errmsg = msg )
            call check_allocation( 'initDerivedCurviGrid2d', 'grid%cell%prob', istat, msg )

            diameter2 = ZERO                    ! for the calculation of the probabilities we
            do iv1 = 1, 3                       ! need the cell diameter
               do iv2 = iv1+1, 4
                  diameter2 = max( diameter2,  &
                                   (grid%cell(ix,iz)%xv(iv2)-grid%cell(ix,iz)%xv(iv1))**2 +  &
                                   (grid%cell(ix,iz)%zv(iv2)-grid%cell(ix,iz)%zv(iv1))**2 )
               end do
            end do
            diameter = sqrt(diameter2)

            do inear = 1, nnear
               ixx = nearcell( isort(inear) )%ix
               izz = nearcell( isort(inear) )%iz
               icc = icfun(ixx,izz)
               grid%cell(ix,iz)%icnear(inear) = icc
               ddist = sqrt( nearcell(isort(inear))%distmin2 ) -  &
                       sqrt( nearcell(isort(1))%distmin2 )
               if( FOUR/MATHPI*ddist/diameter < ONE ) then  ! rough approximation
                  grid%cell(ix,iz)%prob(inear) = (ONE - FOUR/MATHPI*ddist/diameter)**2
               else
                  grid%cell(ix,iz)%prob(inear) = ZERO
               end if
            end do
            grid%cell(ix,iz)%prob(:) = grid%cell(ix,iz)%prob(:) /  &
                                       sum( grid%cell(ix,iz)%prob(:) )
         end do
      end do

   end subroutine initDerivedCurviGrid2d
!
!==========================================================================================
!
   integer(INT4) function icfun( ix, iz )
!                         =====  i   i
!
!     Convert x- and z-index into cell index
!
!  Arguments:
!
      integer(INT4), intent(in) :: ix
      integer(INT4), intent(in) :: iz
!
!  Local Variables:
!
!     None
!
!  Start of Executable Code:
!  ------------------------
!
      icfun = ( iz - 1 ) * grid%nx + ix

   end function icfun
!
!==========================================================================================
!
   integer(INT4) function ixfun( ic )
!                         =====  i
!
!     Get x-index from cell index
!
!  Arguments:
!
      integer(INT4), intent(in) :: ic
!
!  Local Variables:
!
!     None
!
!  Start of Executable Code:
!  ------------------------
!
      ixfun = mod( ic-1, grid%nx ) + 1

   end function ixfun
!
!==========================================================================================
!
   integer(INT4) function izfun( ic )
!                         =====  i
!
!     Get z-index from cell index
!
!  Arguments:
!
      integer(INT4), intent(in) :: ic
!
!  Local Variables:
!
!     None
!
!  Start of Executable Code:
!  ------------------------
!
      izfun = (ic-1) / grid%nx + 1

   end function izfun
!
!==========================================================================================
!
end module Target_Grid2d_Curvi_mod
