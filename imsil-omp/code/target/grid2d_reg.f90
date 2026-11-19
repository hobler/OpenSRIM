module Target_Grid2d_Reg_mod
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
   use Param_Mat_mod,            only: mat, irxtal
   use Param_Sn_mod,             only: sn
   use Param_Dam_mod,            only: dam
   use Param_Sn_mod,             only: sn
   use Param_Out_mod,            only: out

   use Types_mod,                only: Proj_type
   use Geom_Line_mod,            only: mcisct
   use Geom_Polygon_Old_mod,     only: mcauxl, mcdis4
   use Geom_FFPcyl_mod,          only: isInsideFFPcyl
   use Math_Sort_mod,            only: quickSort
   use Target_Surface2d_mod,     only: getDistSurface2d
   use Target_Surface2d_SDF_mod, only: initSurface2dSDFfromPolygons
#ifdef AMAT
   use Target_Mechanics2d_mod,   only: relaxAreas2d
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
   public :: initRegGrid2d                ! initialize the target
!            =============
   public :: initPostCellRegGrid2d        ! grid initializations after cell initializations
!            =====================
   public :: getCellRegGrid2d             ! determines the cell a given point is located in
!            ================
   public :: getOtherNearCellsRegGrid2d   ! determine cells other than nearest cell that
!            ==========================   !    intersect with the free flight path cylinder
   public :: getCellsRegGrid2d            ! get the cell infos
!            =================
   public :: getOvlpAreasRegGrid2d        ! obtain the cell areas inside the polygon
!            =====================
   public :: relaxRegGrid2d               ! relax the cell densities
!            ==============
   public :: writeRegGrid2d               ! write target grid to file
!            ==============
   public :: readRegGrid2d                ! read target grid from file
!            =============
!
!  Private Procedures:
!  ==================
!
   private :: initOnceRegGrid2d           ! initializations in initRegGrid2d and readRegGrid2d,
!             =================           !    but not after relaxRegGrid2d
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

   type :: RegGrid_type                         ! regular 2d grid
!          ============
      integer(INT4) :: nx                       ! number of cells in x-direction
      integer(INT4) :: nz                       ! number of cells in z-direction
      real(REAL8) :: wx                         ! cell width in x-direction
      real(REAL8) :: wz                         ! cell width in z-direction
      real(REAL8), allocatable :: x(:)          ! (0:nx) x-coordinates of grid positions
      real(REAL8), allocatable :: z(:)          ! (0:nz) z-coordinates of grid positions
      real(REAL8) :: xcmin                      ! minimum cell center position in x-direction
      real(REAL8) :: xcmax                      ! maximum cell center position in x-direction
      real(REAL8) :: zcmin                      ! minimum cell center position in z-direction
      real(REAL8) :: zcmax                      ! maximum cell center position in z-direction
      type(Cell_type), allocatable :: cell(:,:) ! cell info
   end type RegGrid_type

   type(RegGrid_type) :: grid

   type :: CurviGrid_type                    ! curvilinear 2d grid
!          ==============
      integer(INT4) :: nx                    ! number of cells in x-direction
      integer(INT4) :: nz                    ! number of cells in z-direction
      real(REAL8), allocatable :: x(:,:)     ! (0:nx,0:nz) x-coordinates of grid positions
      real(REAL8), allocatable :: z(:,:)     ! (0:nx,0:nz) z-coordinates of grid positions
   end type CurviGrid_type

   type(CurviGrid_type) :: relaxedgrid

!   integer(INT4), allocatable :: icell(:,:)     ! Conversion table from Cartesian indices to
!                                                !    cell indices
!
!==========================================================================================
!
contains
!
!==========================================================================================
!
   subroutine initRegGrid2d( )
!             =============
!
!     Initialize the grid
!
!  Arguments:
!
!     None
!
!  Local Variables:
!
      integer(INT4) :: ix, iz, ixmin, izmin, ixmax, izmax
      real(REAL8) :: dx, offset(2)
!
!  Start of Executable Code:
!  ------------------------
!
      grid%wx = out%wbox2(1)
      grid%wz = out%wbox2(2)
!
!  Number of cells
!                                                     ! avoid cell boundaries at lattice planes
      if( irxtal > 0 ) then
         offset(:) = -100 * EPSMA * [ maxval(geom%region(0)%px(:)),  &
                                       maxval(geom%region(0)%pz(:)) ]
      else
         offset(:) = ZERO
      end if

      dx = sn%pmax*(ONE+TWO*EPSMB)
      ixmin = ceiling( ( minval(geom%region(0)%px(:)) - dx - offset(1) ) / grid%wx )
      izmin = ceiling( ( minval(geom%region(0)%pz(:)) - dx - offset(2) ) / grid%wz )
      ixmax = ceiling( ( maxval(geom%region(0)%px(:)) + dx - offset(1) ) / grid%wx )
      izmax = ceiling( ( maxval(geom%region(0)%pz(:)) + dx - offset(2) ) / grid%wz )
                                                ! adding a safety border of pmax is not
      grid%nx = ixmax - ixmin + 1               !    strictly logical. It should rather be
      grid%nz = izmax - izmin + 1               !    the maximum displacement by the time
!                                               !    integral
!  Position of cell interfaces
!
      allocate( grid%x( 0:grid%nx ) )
      allocate( grid%z( 0:grid%nz ) )
      grid%x = [ (offset(1) + ix * grid%wx, ix=ixmin-1,ixmax) ]
      grid%z = [ (offset(2) + iz * grid%wz, iz=izmin-1,izmax) ]

      call initOnceRegGrid2d( )
!          =================

   end subroutine initRegGrid2d
!
!===================================================================================
!
   subroutine initPostCellRegGrid2d( ff, empty )
!             =====================  i     i
!
!     Initialize derived quantities of the cell structure of the 2d target
!
!  Arguments:
!
      real(REAL8), intent(in) :: ff(:)             ! fill fractions
      logical, intent(in) :: empty(:)              ! indicate empty cells
!
!  Local Variables:
!
      character(80) :: text
      character(256) :: errmsg
      logical :: inside, reentr(4)
      integer(INT4) :: ix, iz, ixx, izz, i, n, ncn, istat, inear, nnear, stride, ic, icc
      integer(INT4) :: iv1, iv2
      real(REAL8) :: x1, x2, z1, z2, rsqi, cosetac, sinetac
      real(REAL8) :: distmin2, distmin2_1, distmin2_2, distmax2, distmax2_1, distmax2_2
      real(REAL8) :: diameter2, diameter, ddist
      real(REAL8) :: px_c(5), pz_c(5)

      real(REAL8) :: xc(0:grid%nx+1), zc(0:grid%nz+1)
      real(REAL8) :: ff2(0:grid%nx+1,0:grid%nz+1)

      integer(INT4), allocatable :: ipend(:)
      real(REAL8), allocatable :: px(:), pz(:)

      integer(INT4) :: isort(grid%nx*grid%nz)
      real(REAL8) :: distmin2_tmp(grid%nx*grid%nz)

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

      type(Nearcell_type) :: nearcell(grid%nx*grid%nz)
!
!  Start of Executable Code:
!  ------------------------
!
!  Curvilinear grid set to a Cartesian/regular grid
!
      if( .not. allocated( relaxedgrid%x ) ) then
         relaxedgrid%nx = grid%nx
         relaxedgrid%nz = grid%nz

         allocate( relaxedgrid%x( 0:relaxedgrid%nx, 0:relaxedgrid%nz ), stat=istat,  &
                   errmsg=errmsg )
         call check_allocation( 'initPostCellRegGrid2d', 'relaxedgrid%x', istat, errmsg )

         allocate( relaxedgrid%z( 0:relaxedgrid%nx, 0:relaxedgrid%nz ), stat=istat,  &
                   errmsg=errmsg )
         call check_allocation( 'initPostCellRegGrid2d', 'relaxedgrid%z', istat, errmsg )

      end if

      relaxedgrid%x(:,:) = spread( grid%x(:), dim=2, ncopies=relaxedgrid%nz+1 )
      relaxedgrid%z(:,:) = spread( grid%z(:), dim=1, ncopies=relaxedgrid%nx+1 )
!
!  Initialize near cells
!
      rsqi = QUART * ( grid%wx**2 + grid%wz**2 )    ! cell radius ^2

      px_c(:) = [ZERO,grid%wx,grid%wx,ZERO,ZERO]  ! auxiliary line for rectangles
      pz_c(:) = [ZERO,ZERO,grid%wz,grid%wz,ZERO]

      call mcauxl( cosetac, sinetac, px_c, pz_c )
!          ======    o          o      i     i

      reentr(:) = .false.                                   ! no reentrant cell corners

      do iz = 1, grid%nz
         z1 = grid%z(iz-1)
         z2 = grid%z(iz)
         do ix = 1, grid%nx
            x1 = grid%x(ix-1)
            x2 = grid%x(ix)

            ic = icfun(ix,iz)
            deallocate( grid%cell(ix,iz)%icnear )
            deallocate( grid%cell(ix,iz)%prob )

            if( .not. empty(ic) ) then
               ! The single near cell of a non-empty cell is the cell itself
               allocate( grid%cell(ix,iz)%icnear( 1 ), stat=istat, errmsg=errmsg )
               write( text, '(a,i0,a,i0,a)' )  'grid%cell(', ix, ',', iz, ')%icnear'
               call check_allocation( 'initPostCellRegGrid2d', text, istat, errmsg )
               grid%cell(ix,iz)%icnear(1) = ic
               !    empty cell is the cell itself
               allocate( grid%cell(ix,iz)%prob( 1 ), stat=istat, errmsg=errmsg )
               write( text, '(a,i0,a,i0,a)' )  'grid%cell(', ix, ',', iz, ')%prob'
               call check_allocation( 'initPostCellRegGrid2d', text, istat, errmsg )
               grid%cell(ix,iz)%prob(1) = ONE

               cycle

            end if

            nearcell(:)%distmin2 = BIG          ! distmin2, distmax2 from the current cell
            nearcell(:)%distmax2 = BIG          !    (ix,iz)
            distmin2_1 = BIG                    ! smallest distmin2
            distmin2_2 = BIG                    ! second-smallest distmin2
            distmax2_1 = BIG                    ! distmax2 corresponding to smallest distmin2
            distmax2_2 = BIG                    ! distmax2 corresponding to second-smallest
                                                !    distmin2
            ncn = 0
            !
            !  loop over shells (should always exit via exit)
            !
            do n = 1, max( grid%nx, grid%nz )
               if( (n-1) * min( grid%wx, grid%wz ) > distmax2_2 )  exit
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
                     if( empty(icc) )  cycle

                     call mcdis4( distmin2, distmax2, inside, x1, x2, z1, z2, rsqi,  &
!                         ======     o         o        o     i   i   i   i    i
                                  grid%cell(ixx,izz)%xv(:), grid%cell(ixx,izz)%zv(:),  &
!                                         i                    i
                                  reentr, cosetac, sinetac )
!                                   i        i        i

                     if( inside )  &
                        error stop 'initPostCellRegGrid2d: cell center cannot be within another cell'
                     !
                     !  maintain distances from second-nearest cell
                     !
                     if( distmin2 < distmin2_1 ) then
                        distmin2_2 = distmin2_1
                        distmin2_1 = distmin2
                        distmax2_2 = distmax2_1
                        distmax2_1 = distmax2
                     else if( distmin2 < distmin2_2 ) then
                        distmin2_2 = distmin2
                        distmax2_2 = distmax2
                     end if

                     ncn = ncn + 1

                     nearcell(ncn)%ix = ixx
                     nearcell(ncn)%iz = izz
                     nearcell(ncn)%distmin2 = distmin2
                     nearcell(ncn)%distmax2 = distmax2

                  end do
               end do
            end do
            !
            !  Keep only near cells that obey distmin2 < distmax2_2 and sort according to
            !  distmin2
            !
            nnear = count( nearcell(1:ncn)%distmin2 <= distmax2_2 )
            nearcell(1:nnear) = pack( nearcell(1:ncn), nearcell(1:ncn)%distmin2 <= distmax2_2 )

            isort(1:nnear) = [ (i, i=1,nnear) ]
            distmin2_tmp(1:nnear) = nearcell(1:nnear)%distmin2
            call quickSort( distmin2_tmp(1:nnear), isort(1:nnear), 1, nnear )
!                =========          io                  io         i    i

            allocate( grid%cell(ix,iz)%icnear( nnear ), stat=istat, errmsg=errmsg )
            write( text, '(a,i0,a,i0,a)' )  'grid%cell(', ix, ',', iz, ')%icnear'
            call check_allocation( 'initPostCellRegGrid2d', text, istat, errmsg )

            allocate( grid%cell(ix,iz)%prob( nnear ), stat=istat, errmsg=errmsg )
            write( text, '(a,i0,a,i0,a)' )  'grid%cell(', ix, ',', iz, ')%prob'
            call check_allocation( 'initPostCellRegGrid2d', text, istat, errmsg )

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
!
!  Calculate grid for fill fractions (extended by 1 row on each side)
!
      xc(1:grid%nx) = HALF * ( grid%x(0:grid%nx-1) + grid%x(1:grid%nx) )
      zc(1:grid%nz) = HALF * ( grid%z(0:grid%nz-1) + grid%z(1:grid%nz) )
      xc(0) = TWO * xc(1) - xc(2)
      zc(0) = TWO * zc(1) - zc(2)
      xc(grid%nx+1) = TWO * xc(grid%nx) - xc(grid%nx-1)
      zc(grid%nz+1) = TWO * zc(grid%nz) - zc(grid%nz-1)
!
!  Expand fill fraction array
!
      ff2(:,:) = ZERO                           ! outer cells will remain 0
      ff2(1:grid%nx,1:grid%nz) = max( ZERO, min( ONE, reshape( ff, [grid%nx,grid%nz] ) ) )

      if( .not. geom%freex(1) )  ff2(0,:) = ff2(1,:)                 ! outer cells, taking into
      if( .not. geom%freex(2) )  ff2(grid%nx+1,:) = ff2(grid%nx,:)   !    account the real/
      if( .not. geom%freez(1) )  ff2(:,0) = ff2(:,1)                 !    artificial boundary
      if( .not. geom%freez(2) )  ff2(:,grid%nz+1) = ff2(:,grid%nz)
!
!  Initialize the surface from the fill fractions
!
      call findZeroLevelSetCartGrid2d( ff2-HALF, ipend, px, pz )
!          ==========================     i       o    o   o

      px(:) = xc(0) + grid%wx * px(:)
      pz(:) = zc(0) + grid%wz * pz(:)

      call initSurface2dSDFfromPolygons( xc, zc, ipend, px, pz )
!          ============================  i   i     i    i   i

      deallocate( ipend, px, pz )

   end subroutine initPostCellRegGrid2d
!
!==========================================================================================
!
   subroutine getCellRegGrid2d( x, z, ic, outside_ )
!             ================  i  i  o      o
!
!     Determine cell the point pos is located in
!
!  Arguments:
!
      real(REAL8), intent(in) :: x              ! x-coordinate of position to be checked
      real(REAL8), intent(in) :: z              ! z-coordinate of position to be checked
      integer(INT4), intent(out) :: ic          ! index of cell (region)
      logical, optional, intent(out) :: outside_! indicates whether pos is outside all cells
!
!  Local Variables:
!
      integer(INT4) :: ix, iz
      real(REAL8) :: xx, zz
!
!  Start of Executable Code:
!  ------------------------
!
!  Determine cell pos is located (positions outside grid map to nearest cell)
!
      if( present( outside_ ) ) then
         outside_ = x < grid%x(0) .or. x > grid%x(grid%nx) .or.  &
                    z < grid%z(0) .or. z > grid%z(grid%nx)
      end if

      xx = max( grid%xcmin, min( grid%xcmax, x ) )   ! don't use min/max with the indices
      zz = max( grid%zcmin, min( grid%zcmax, z ) )   !   to avoid integer overflow

      ix = ceiling( ( xx - grid%x(0) ) / grid%wx )
      iz = ceiling( ( zz - grid%z(0) ) / grid%wz )

      ic = icfun(ix,iz)

   end subroutine getCellRegGrid2d
!
!==========================================================================================
!
   subroutine getOtherNearCellsRegGrid2d( proj, posc, ffp, pmax, ncn, icnear )
!             ==========================   i     i     i    i     o     o
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
      integer(INT4) :: ixmin, izmin, ixmax, izmax
      integer(INT4) :: ix, iz, ixmin_ab, izmin_ab, ixmax_ab, izmax_ab
      real(REAL8) :: z, dz, z1, z2, zmin, zmax, zisct
      real(REAL8) :: x, dx, x1, x2, xmin, xmax, xisct, ffp_p
      real(REAL8) :: ax, az, bx, bz
      real(REAL8) :: rlam1, rlam2, dirx, diry, dirz, dirxy, dirxz, diryz, dirx_p, dirz_p
      logical, allocatable :: overlap(:,:)      ! 2D mask indicating cells overlaping with the
                                                !    ffp cylinder
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
!  Define the rectangle that contains the projection of the ffp cylinder and
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

      xmin = max( grid%xcmin, min( grid%xcmax, xmin ) )     ! don't apply bounds
      xmax = max( grid%xcmin, min( grid%xcmax, xmax ) )     !    to ixmin,ixmax
      zmin = max( grid%zcmin, min( grid%zcmax, zmin ) )     !    because of possible
      zmax = max( grid%zcmin, min( grid%zcmax, zmax ) )     !    integer overflow

      ixmin = ceiling( ( xmin - grid%x(0) ) / grid%wx )  ! indices are valid
      ixmax = ceiling( ( xmax - grid%x(0) ) / grid%wx )  !    cell indices in
      izmin = ceiling( ( zmin - grid%z(0) ) / grid%wz )  !    the range
      izmax = ceiling( ( zmax - grid%z(0) ) / grid%wz )  !    [1..nxz(:)]
!
!  Allocate mask
!
      allocate( overlap( ixmin:ixmax, izmin:izmax ) )
      overlap(:,:) = .false.
!
!  Find the range of cell indices that contain the line A-B
!
      xmin = min( ax, bx )
      xmax = max( ax, bx )
      zmin = min( az, bz )
      zmax = max( az, bz )

      xmin = max( grid%xcmin, min( grid%xcmax, xmin ) )
      xmax = max( grid%xcmin, min( grid%xcmax, xmax ) )
      zmin = max( grid%zcmin, min( grid%zcmax, zmin ) )
      zmax = max( grid%zcmin, min( grid%zcmax, zmax ) )

      ixmin_ab = ceiling( ( xmin - grid%x(0) ) / grid%wx )  ! indices must be in the
      ixmax_ab = ceiling( ( xmax - grid%x(0) ) / grid%wx )  !    range [ixmin(:)...
      izmin_ab = ceiling( ( zmin - grid%z(0) ) / grid%wz )  !    ixmax(:)]
      izmax_ab = ceiling( ( zmax - grid%z(0) ) / grid%wz )  !
!
!  Check for cell corners inside the ffp cylinder. All four adjoining cells overlap the
!  ffp cylinder
!
      ffp_p = ffp * dirxz                          ! p = projected to xz-plane

      do iz = izmin, izmax-1
         do ix = ixmin, ixmax-1

            isInside = .false.

            x1 = grid%x(ix)
            z1 = grid%z(iz)

            isInside = isInsideFFPcyl( [x1,z1], [x,z], [dirx_p,dirz_p], ffp_p, pmax,  &
!                      ==============      i         i             i            i     i
                                       abs(diry) )
!                                          i

            if( isInside )  overlap(ix:ix+1,iz:iz+1) = .true.

         end do
      end do
!
!  It is possible that the ffp cylinder passes between all cell corners. Therefore
!  check in addition for cell edges intersecting the extended ffp (line AB)
!
!  Check intersection with horizontal (z=const) lines
!
      do iz = izmin_ab, izmax_ab-1

         x1 = grid%x(0)
         x2 = grid%x(grid%nx)
         z1 = grid%z(iz)
         z2 = grid%z(iz)

         call mcisct( rlam1, rlam2, cut1, cut2, x1, z1, x2, z2, ax, az, bx, bz )
!             ======    o      o     o     o    i   i   i   i   i   i   i   i

         if( .not. cut2 ) cycle ! no valid intersection only if A^B outside the polygon
                                ! (can this happen?)

         xisct = x1 + ( x2 - x1 ) * rlam1       ! x-coordinate of intersection point

         ix = ceiling( ( xisct - grid%x(0) ) / grid%wx )
         ix = max( ixmin, min( ixmax, ix ) )    ! just in case AB intersects at
                                                !    cell area boundary
         overlap(ix,iz:iz+1) = .true.

      end do
!
!  Check intersection with the vertical (x=const) lines
!
      do ix = ixmin_ab, ixmax_ab-1

         x1 = grid%x(ix)
         x2 = grid%x(ix)
         z1 = grid%z(0)
         z2 = grid%z(grid%nz)

         call mcisct( rlam1, rlam2, cut1, cut2, x1, z1, x2, z2, ax, az, bx, bz )
!             ======    o      o     o     o    i   i   i   i   i   i   i   i

         if( .not. cut2 )  cycle ! (?)

         zisct = z1 + ( z2 - z1 ) * rlam1       ! z-coordinate of intersection point

         iz = ceiling( ( zisct - grid%z(0) ) / grid%wz )
         iz = max( izmin, min( izmax, iz ) )    ! just in case AB intersects at
                                                !    cell area boundary
         overlap(ix:ix+1,iz) = .true.

      end do
!
!  Make sure icnear is big enough
!
      deallocate( icnear )
      allocate( icnear( (ixmax-ixmin+1)*(izmax-izmin+1) ) )
!
!  Include in the near cell list the indices of non-empty cells overlaping with the ffp
!  cylinder other than the present cell
!
      ncn = 0

      do iz = izmin, izmax
         do ix = ixmin, ixmax
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

   end subroutine getOtherNearCellsRegGrid2d
!
!=======================================================================
!
   subroutine getCellsRegGrid2d( cell )
!             =================   o
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
      ncell = grid%nx * grid%nz

      allocate( cell( ncell ), stat=istat, errmsg=errmsg )
      call check_allocation( 'getCellsRegGrid2d', 'cell', istat, errmsg )

      do iz = 1, grid%nz
         do ix = 1, grid%nx                  ! contrary to its name, allocate also copies
            ic = icfun(ix,iz)                !    source to destination

            allocate( cell(ic)%xv, source=grid%cell(ix,iz)%xv, stat=istat, errmsg=errmsg )
            write( text, '(a,i0,a)' )  'cell(', ic, ')%xv'
            call check_allocation( 'getCellsRegGrid2d', text, istat, errmsg )

            allocate( cell(ic)%zv, source=grid%cell(ix,iz)%zv, stat=istat, errmsg=errmsg )
            write( text, '(a,i0,a)' )  'cell(', ic, ')%zv'
            call check_allocation( 'getCellsRegGrid2d', text, istat, errmsg )

            allocate( cell(ic)%icnear, source=grid%cell(ix,iz)%icnear, stat=istat,  &
                      errmsg=errmsg )
            write( text, '(a,i0,a)' )  'cell(', ic, ')%icnear'
            call check_allocation( 'getCellsRegGrid2d', text, istat, errmsg )

            allocate( cell(ic)%prob, source=grid%cell(ix,iz)%prob, stat=istat, errmsg=errmsg )
            write( text, '(a,i0,a)' )  'cell(', ic, ')%prob'
            call check_allocation( 'getCellsRegGrid2d', text, istat, errmsg )

         end do
      end do

   end subroutine getCellsRegGrid2d
!
!=======================================================================
!
   subroutine getOvlpAreasRegGrid2d( px, pz, area, ic_area)
!             =====================  i   i    o      o
!
!     Obtain the overlap areas of the given polygon with the regular grid.
!
!  Arguments:
!
      real(REAL8), intent(in) :: px(:)                   ! x-coordinates of polygon corners
      real(REAL8), intent(in) :: pz(:)                   ! z-coordinates of polygon corners
      real(REAL8), allocatable, intent(out) :: area(:)   ! area inside polygon
                                                         !    first index runs in x-direction
                                                         !    second index runs in z-direction
      integer(INT4), allocatable, intent(out) :: ic_area(:) ! cell indices corresponding to
                                                         !    the areas
!
!  Local Variables:
!
      character(256) :: errmsg
      integer(INT4) :: iarea, narea, istat
      integer(INT4), allocatable :: ix_area(:), iz_area(:)
!
!  Start of Executable Code:
!  ------------------------
!
      call getOvlpAreasCartGrid2d( (px(:)-grid%x(0))/grid%wx, (pz(:)-grid%z(0))/grid%wz,  &
!          ======================            i                          i
                                   grid%nx, grid%nz, area, ix_area, iz_area )
!                                     i        i      o       o        o
      area(:) = area(:)*grid%wx*grid%wz

      narea = size( area )
      allocate( ic_area( narea ), stat=istat, errmsg=errmsg )
      call check_allocation( 'getOvlpAreasRegGrid2d', 'ic_area', istat, errmsg )

      do iarea = 1, narea
         ic_area(iarea) = icfun(ix_area(iarea),iz_area(iarea))
      end do

      deallocate( ix_area, iz_area )

   end subroutine getOvlpAreasRegGrid2d
!
!=======================================================================
!
   subroutine relaxRegGrid2d( ff, dpa, ovlp )
!             ==============  i    i    o
!
!     Redistribute the cell densities so that stress relaxes. This is done in two steps:
!     First the geometry of the grid is adjusted so that the cell densities (given by the
!     cell content and the cell size) are identical to the equilibrium densities for the
!     respective target composition. Then the cell contents are mapped to an equidistant
!     grid again.
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

      real(REAL8), intent(in) :: ff(:)          ! (ncell) fill fractions
      real(REAL8), intent(in) :: dpa(:)         ! (ncell) displacements per atom during BC step
      type(Ovlp_type), intent(out) :: ovlp(:)   ! (ncell) overlap information
!
!  Local Variables:
!
      logical :: mask(grid%nx,grid%nz)
      logical :: movable(0:grid%nx,0:grid%nz)
      integer(INT4) :: ix, iz, ic
      integer(INT4) :: ix1, ix2, iz1, iz2
      real(REAL8) :: x4(5), z4(5)
      real(REAL8) :: dist
      real(REAL8) :: area(grid%nx,grid%nz)
!
!  Start of Executable Code:
!  ------------------------
!
!  Calculate relaxed cell sizes and densities
!
      area = reshape( ff, [ grid%nx, grid%nz ] ) * grid%wx * grid%wz
!!
!      call bndConnected( mask )                 ! empty cells connected to boundary
!!          ============   io
!      mask(:,:) = .not.mask(:,:)                ! full cells and unconnected empty cells
!
!  Relax cells completely inside target
!
      mask(:,:) = .true.        ! areas that may be relaxed
      movable(:,:) = .true.     ! cell corners that may be moved

      do ix = 0, grid%nx
         do iz = 0, grid%nz
            call getDistSurface2d( [grid%x(ix), ZERO, grid%z(iz)], dist )
!                ================                 i                 o
            ix1 = max( ix, 1 )
            ix2 = min( ix+1, grid%nx )
            iz1 = max( iz, 1 )
            iz2 = min( iz+1, grid%nz )
            mask(ix1:ix2,iz1:iz2) = mask(ix1:ix2,iz1:iz2) .and. (dist > ZERO)
!            movable(ix,iz) = all(  &
!               dpa(reshape( icfun(ix1:ix2,iz1:iz2), [ (ix2-ix1+1)*(iz2-iz1+1) ] )) >  &
!               dam%dparelax )
            movable(ix,iz) = (dpa(icfun(ix1,iz1)) > dam%dparelax ) .and.  &
                             (dpa(icfun(ix1,iz2)) > dam%dparelax ) .and.  &
                             (dpa(icfun(ix2,iz1)) > dam%dparelax ) .and.  &
                             (dpa(icfun(ix2,iz2)) > dam%dparelax )
         end do
      end do

      print *, count(movable(:,:)), 'movable points'
      if( dam%relaxmodel == ENUM%relaxmodel%none ) then
         area(:,:) = grid%wx * grid%wz
      else if( dam%relaxmodel == ENUM%relaxmodel%full ) then
#ifdef AMAT
         call relaxAreas2d( relaxedgrid%x, relaxedgrid%z, area, mask, movable )
!             ============       io             io         io    i       i
#else
         error stop 'relaxRegGrid2d: relaxAreas2d not available'
#endif
      else
         error stop 'relaxRegGrid2d: Invalid dam%relaxmodel'
      end if
!
!  Determine overlap areas with original grid   !!! can grid move out of original grid?
!
      do iz = 1, grid%nz                    ! loop over relaxed cells
         do ix = 1, grid%nx
            ic = icfun(ix,iz)

            x4(1:4) = [ relaxedgrid%x(ix-1,iz-1), relaxedgrid%x(ix,iz-1),  &
                         relaxedgrid%x(ix,iz), relaxedgrid%x(ix-1,iz) ]
            z4(1:4) = [ relaxedgrid%z(ix-1,iz-1), relaxedgrid%z(ix,iz-1),  &
                         relaxedgrid%z(ix,iz), relaxedgrid%z(ix-1,iz) ]
            x4(5) = x4(1)
            z4(5) = z4(1)

            call getOvlpAreasRegGrid2d( x4(:), z4(:), ovlp(ic)%area, ovlp(ic)%ic )
!                =====================    i      i           o              o
            ovlp(ic)%emptyarea = ZERO

         end do
      end do

   end subroutine relaxRegGrid2d
!
!=======================================================================
!
   subroutine writeRegGrid2d( nout, relaxed_ )
!             ==============   i       i
!
!     Write the grid to the file with unit number nout
!
!  Arguments:
!
      integer(INT4), intent(in) :: nout         ! unit number
      logical, optional :: relaxed_             ! indicates relaxed grid to be written
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
         write( nout, '(2i15)' )  grid%nx, grid%nz
         do iz = 0, grid%nz
            do ix = 0, grid%nx
               write( nout, '(2f15.6)' )  relaxedgrid%x(ix,iz), relaxedgrid%z(ix,iz)
            end do
         end do
      else
         write( nout, '(a)' )  '# Grid: regular'
         write( nout, '(2i15)' )  grid%nx, grid%nz
         write( nout, '(2f15.6)' )  grid%x(0), grid%z(0)
         write( nout, '(2f15.6)' )  grid%wx, grid%wz
      end if

   end subroutine writeRegGrid2d
!
!=======================================================================
!
   subroutine readRegGrid2d( nin, valid )
!             =============   i     o
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
         read( nin, *, err=900, end=900 )  grid%nx, grid%nz
         if( allocated( grid%x ) )  deallocate( grid%x, grid%z )
         allocate( grid%x( 0:grid%nx ) )
         allocate( grid%z( 0:grid%nz ) )
         read( nin, *, err=900, end=900 )  grid%x(0), grid%z(0)
         read( nin, *, err=900, end=900 )  grid%wx, grid%wz
         grid%x(1:grid%nx) = [ (grid%x(0) + ix * grid%wx, ix = 1, grid%nx) ]
         grid%z(1:grid%nz) = [ (grid%z(0) + iz * grid%wz, iz = 1, grid%nz) ]
         valid = .true.
      else
         goto 901
      end if
!
!  Other initializations
!
      call initOnceRegGrid2d( )
!          =================

      return

900   error stop 'readRegGrid2d: error during reading CELL2 file'
901   error stop 'readRegGrid2d: unknown grid type'

   end subroutine readRegGrid2d
!
!=======================================================================
!
   subroutine initOnceRegGrid2d( )
!             =================
!
!     Initializations of the grid common to initRegGrid2d and readRegGrid2d
!
!  Arguments:
!
!     None
!
!  Local Variables:
!
      character(80) :: text
      character(256) :: errmsg
      integer(INT4) :: istat, ix, iz, ncell
!
!  Start of Executable Code:
!  ------------------------
!
!  Number of cells
!
      ncell = grid%nx * grid%nz
!
!  Minimum and maximum cell center positions
!
      grid%xcmin = HALF * ( grid%x(0) + grid%x(1) )
      grid%xcmax = HALF * ( grid%x(grid%nx-1) + grid%x(grid%nx) )
      grid%zcmin = HALF * ( grid%z(0) + grid%z(1) )
      grid%zcmax = HALF * ( grid%z(grid%nz-1) + grid%z(grid%nz) )
!
!  Allocate cells and initialize cell vertices
!
      allocate( grid%cell( grid%nx, grid%nz ), stat=istat, errmsg=errmsg )
      call check_allocation( 'initOnceRegGrid2d', 'grid%cell', istat, errmsg )

      do iz = 1, grid%nz
         do ix = 1, grid%nx

            allocate( grid%cell(ix,iz)%xv( 5 ), stat=istat, errmsg=errmsg )
            write( text, '(a,i0,a,i0,a)' )  'grid%cell(', ix, ',', iz, ')%xv'
            call check_allocation( 'initOnceRegGrid2d', text, istat, errmsg )

            allocate( grid%cell(ix,iz)%zv( 5 ), stat=istat, errmsg=errmsg )
            write( text, '(a,i0,a,i0,a)' )  'grid%cell(', ix, ',', iz, ')%zv'
            call check_allocation( 'initOnceRegGrid2d', text, istat, errmsg )

            grid%cell(ix,iz)%xv(1:5) = grid%x( [ix-1,ix,ix,ix-1,ix-1] )
            grid%cell(ix,iz)%zv(1:5) = grid%z( [iz-1,iz-1,iz,iz,iz-1] )
            ! dummy allocations
            allocate( grid%cell(ix,iz)%icnear( 1 ), stat=istat, errmsg=errmsg )
            write( text, '(a,i0,a,i0,a)' )  'grid%cell(', ix, ',', iz, ')%icnear'
            call check_allocation( 'initOnceRegGrid2d', text, istat, errmsg )

            allocate( grid%cell(ix,iz)%prob( 1 ), stat=istat, errmsg=errmsg )
            write( text, '(a,i0,a,i0,a)' )  'grid%cell(', ix, ',', iz, ')%prob'
            call check_allocation( 'initOnceRegGrid2d', text, istat, errmsg )

         end do
      end do

   end subroutine initOnceRegGrid2d
!
!=======================================================================
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
end module Target_Grid2d_Reg_mod
