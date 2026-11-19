module Target_Cell2d_Reg_mod
!
!  Used Modules:
!  ============
!
   use Const_mod,                only: INT4, INT8, REAL8,  &
                                       ZERO, HALF, ONE, TWO, BIG, EPSMA, EPSMB,  &
                                       ENUM
   use Array_mod,                only: alloc, changeShape, addReflect
   use String_mod,               only: str
   use Types_mod,                only: Proj_type

   use Param_Setup_mod,          only: nr, nelem, usecell
   use Param_Elem_mod,           only: elem
   use Param_Geom_mod,           only: geom
   use Param_Mat_mod,            only: mat, irxtal
   use Param_Sn_mod,             only: sn
   use Param_Dam_mod,            only: dam
   use Param_Out_mod,            only: out

   use Geom_FFPcyl_mod,          only: isInsideFFPcyl
   use Geom_CartGrid2d_mod,      only: getOvlpAreasCartGrid2d
   use Geom_Line_mod,            only: mcisct
   use Geom_Polygon_mod,         only: initPolygon, isInsidePolygon
   use Math_IntMath_mod,         only: index1to2, index2to1

   use Stat_His2d_mod,           only: His2d_type, writeHis2d
   use Target_Surface2d_SDF_mod, only: getGridSurface2dSDF
#ifdef AMAT
   use Cinterface_Mechanics2d_mod,  only: relaxViscoElasticEuler2d, initSurface2dFromFF
#endif
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
   public :: initCell2dReg                   ! initialize grid and cells
!            =============
   public :: getNearestCell2dReg             ! get nearest non-empty cell index
!            ===================
   public :: getOtherNearCellIndices2dReg    ! determine cells other than nearest cell that
!            ============================    !    intersect with the free flight path cylinder
   public :: fixNegDensCell2dReg             ! fix negative densities by drawing material from
!            ===================             !    neighboring cells
   public :: relaxCell2dReg                  ! relax the cell densities
!            ==============
   public :: writeCell2dReg                  ! write target cell contents to CELL2 file
!            ==============
!
!  Private Procedures:
!  ==================
!
   private :: initCell2dRegFromGeom          ! initialize grid and cell from region geometry
!             =====================
!
!=======================================================================
!
contains
!
!=======================================================================
!
   subroutine initCell2dReg( x, z, dens, ffgeom, lxtal, pressure, stress )
!             ============= io io  io     io     io       io       io
!
!     Read the contents of the target cells from the cellfile file (same format as HIS file)
!
!     At the target boundaries, the densities are averaged over the cell.
!     Geometric fill fractions (ffgeom) are determined
!
!     TODO: eliminate shiftto and reflect or consider it in initSurface2d
!
!  Arguments:
!
      real(REAL8), allocatable, intent(inout) :: x(:)             ! x positions of grid points
      real(REAL8), allocatable, intent(inout) :: z(:)             ! z positions of grid points
      real(REAL8), allocatable, intent(inout) :: dens(:,:,:)      ! (nelem,nc) atom densities
      real(REAL8), allocatable, intent(inout) :: ffgeom(:,:)      ! geometric fill fraction
      logical, allocatable, intent(inout) :: lxtal(:,:)           ! crystallinity flag
      real(REAL8), allocatable, intent(inout) :: pressure(:,:)    ! pressure
      real(REAL8), allocatable, intent(inout) :: stress(:,:,:,:)  ! 2d stress tensor
!
!  Local Variables:
!
      logical :: true(4) = [.true., .true., .true., .true.]
      integer(INT4) :: direction(2)
!
!  Start of Executable Code:
!  ------------------------
!
!  Initialize the arrays
!
      if( usecell ) then
!         call readCell2dReg( z, dens, ffgeom, lxtal, pressure, stress )
!!             =============  io  io     io     io       io       io
         error stop 'initCell2dReg: reading from 2d cell file not yet implemented'
      else
         call initCell2dRegFromGeom( x, z, dens, ffgeom, lxtal, pressure, stress )
!             ===================== io io   io     io     io       io       io
      end if
!
!  Apply reflective boundary conditions (if any)
!
      if( geom%bndtypex(1) == ENUM%bndtype%reflect ) then
         direction(1) = -1
         geom%bndtypex(1) = geom%bndtypex(2)
      else if( geom%bndtypex(2) == ENUM%bndtype%reflect ) then
         direction(1) = +1
         geom%bndtypex(2) = geom%bndtypex(1)
      else
         direction(1) = 0
      end if

      if( geom%bndtypez(1) == ENUM%bndtype%reflect ) then
         direction(2) = -1
         geom%bndtypez(1) = geom%bndtypez(2)
      else if( geom%bndtypez(2) == ENUM%bndtype%reflect ) then
         direction(2) = +1
         geom%bndtypez(2) = geom%bndtypez(1)
      else
         direction(2) = 0
      end if

      call addReflect( x, direction(1), .false., 'x', 'initCell2dReg' )
!          ========== io       i           i      i          i
      call addReflect( z, direction(2), .false., 'z', 'initCell2dReg' )
!          ========== io       i           i      i          i
      call addReflect( dens, [0, direction], true(1:3), 'dens', 'initCell2dReg' )
!          ==========   io         i             i        i            i
      call addReflect( ffgeom, direction, true(1:2), 'ffgeom', 'initCell2dReg' )
!          ==========    io        i          i         i             i
      call addReflect( lxtal, direction, true(1:2), 'lxtal', 'initCell2dReg' )
!          ==========   io        i          i         i            i
      call addReflect( pressure, direction, true(1:2), 'pressure', 'initCell2dReg' )
!          ==========     io         i          i          i              i
      call addReflect( stress, [0, 0, direction], true(1:4), 'stress', 'initCell2dReg' )
!          ==========    io          i                i         i             i

   end subroutine initCell2dReg
!
!=======================================================================
!
   subroutine initCell2dRegFromGeom( x, z,  &
!             =====================  o  o
                                     dens, ffgeom, lxtal, pressure, stress )
!                                     io     io     io       io       io
!
!     Initialize the cell structure of the target from a 2d region definition by polygons.
!
!     The target extends from minval(geom%region(0)%px) to maxval(geom%region(0)%px) and from
!     minval(geom%region(0)%pz) to maxval(geom%region(0)%pz) in x and z direction, respectively.
!     On entry, geom%xlim(:) and geom%zlim(:) have already been initialized to include a vacuum
!     surface layer at free boundaries and the reflected region at reflected boundaries.
!
!     This procedure does NOT perform the reflection at a reflective boundary.
!
!  Arguments:
!
      real(REAL8), allocatable, intent(out) :: x(:)            ! x coordinates of grid points
      real(REAL8), allocatable, intent(out) :: z(:)            ! z coordinates of grid points
      real(REAL8), allocatable, intent(inout) :: dens(:,:,:)   ! (nelem,nx,nz) densities
      real(REAL8), allocatable, intent(inout) :: ffgeom(:,:)   ! (nx,nz) geometric fill fraction
      logical, allocatable, intent(inout) :: lxtal(:,:)        ! crystallinity flag
      real(REAL8), allocatable, intent(inout) :: pressure(:,:) ! pressure
      real(REAL8), allocatable, intent(inout) :: stress(:,:,:,:)  ! components of stress tensor
!
!  Local Variables:
!
      integer(INT4) :: ix, iz, iloc, iffgeomc, nx, nz
      integer(INT4) :: ir, ielemr, ielem
      integer(INT4), allocatable :: ix_ffgeomc(:), iz_ffgeomc(:)
      real(REAL8) :: xc, zc, ddens
      real(REAL8), allocatable :: ffgeomc(:)
!
!  Start of Executable Code:
!  ------------------------
!
!  Get grid from the surface SDF
!
      call getGridSurface2dSDF( x, z )
!          ===================  o  o
      nx = ubound( x, 1 )
      nz = ubound( z, 1 )
!
!  Determine fill fractions
!
      call alloc( ffgeom, [nx, nz], 'ffgeom', 'initCell2dRegFromGeom' )
!          =====    io       i         i                 i
      ffgeom(:,:) = ZERO

      call getOvlpAreasCartGrid2d( (geom%region(0)%px(:)-x(0)) / geom%wbox2(1),  &
!          ======================                   i
                                   (geom%region(0)%pz(:)-z(0)) / geom%wbox2(2),  &
!                                                   i
                                   nx, nz, ffgeomc, ix_ffgeomc, iz_ffgeomc )
!                                  i   i      o           o         o
      do iffgeomc = 1, size( ffgeomc )
         ffgeom(ix_ffgeomc(iffgeomc),iz_ffgeomc(iffgeomc)) = ffgeomc(iffgeomc)
      end do
!
!  Cell densities
!
      call alloc( dens, [nelem, nx, nz], 'dens', 'initCell2dRegFromGeom' )
!          =====   io             i               i                i
      dens(:,:,:) = ZERO

      do ir = 1, nr
         do ielemr = 1, mat(ir)%nelemr
            ielem = mat(ir)%ielem(ielemr)
            do iz = 1, nz
               do ix = 1, nx
                  ddens = mat(ir)%dens * mat(ir)%prob(ielemr) * ffgeom(ix,iz)
                  dens(ielem,ix,iz) = dens(ielem,ix,iz) + ddens
               end do
            end do
         end do
      end do
!
!  Pressure and stress
!
      if( dam%relaxmodel == ENUM%relaxmodel%viscoelastic ) then
         call alloc( pressure, [nx, nz], 'pressure', ' initCell2dRegFromGeom' )
!             =====     io        i          i                   i
         pressure(:,:) = ZERO
         call alloc( stress, [2, 2, nx, nz], 'stress', ' initCell2dRegFromGeom' )
!             =====    io           i           i                  i
         stress(:,:,:,:) = ZERO
      end if
!
!  Crystal flag
!
      call alloc( lxtal, [nx, nz], 'lxtal', ' initCell2dRegFromGeom' )
!          =====   io       i         i                  i
      lxtal(:,:) = .false.
      if( irxtal > 0 ) then
         call initPolygon( geom%region(irxtal) )
!             ===========         io
         do iz = 1, nz
            zc = HALF * (z(iz) - z(iz-1))
            do ix = 1, nx
               xc = HALF * (x(ix) - x(ix-1))
               iloc = isInsidePolygon( [xc, zc ] , geom%region(irxtal) )  ! iloc=1 if inside
!                     ===============      i                 i
               lxtal(ix,iz) = ( iloc > 0 )
            end do
         end do
      end if

   end subroutine initCell2dRegFromGeom
!
!=======================================================================
!
   subroutine getNearestCell2dReg( ic, nx, nz, wx, wz, dens, jc )
!             ===================  i   i   i   i   i    i    o
!
!     Determine the nearest cell with dens /= 0
!
!  Arguments:
!
      integer(INT4), intent(in) :: ic           ! index of cell to which the nearest non-empty
                                                !    cell is sought
      integer(INT4), intent(in) :: nx           ! extents of cell arrray in x direction
      integer(INT4), intent(in) :: nz           ! extents of cell arrray in z direction
      real(REAL8), intent(in) :: wx             ! cell width in x direction
      real(REAL8), intent(in) :: wz             ! cell width in z direction
      real(REAL8), intent(in) :: dens(:)        ! densities
      integer(INT4), intent(out) :: jc          ! index of non-empty cell nearest to ic
!
!  Local Variables:
!
      integer(INT4) :: ix, iz, ishell, jx, jz, jcc, stride
      real(REAL8) :: dist2, dist2min
!
!  Start of Executable Code:
!  ------------------------
!
      call index1to2( ic, nx, ix, iz )
!          =========  i   i   o   o

      dist2min = BIG

      do ishell = 0, max( nx, nz )

         if( ishell*min( wx, wz )**2 > dist2min )  exit  ! no cell in the remaining layers
         !                                               !    can be nearer that dist2min
         !  loop over cells in shell
         !
         do jx = ix-ishell, ix+ishell
            if( jx < 1 .or. jx > nx )  cycle          ! skip cells outside the grid
            if( jx == ix-ishell .or. jx == ix+ishell )  then
               stride = 1
            else
               stride = 2*ishell
            end if
            do jz = iz-ishell, iz+ishell, stride
               if( jz < 1 .or. jz > nz )  cycle       ! skip cells outside the grid
               call index2to1( jx, jz, nx, jc )
!                   =========  i   i   i    o
               if( dens(jc) == ZERO )  cycle          ! skip empty cells

               dist2 = ((jx-ix)*wx)**2 + ((jz-iz)*wz)**2
               if( dist2 < dist2min ) then
                  dist2min = dist2
                  jcc = jc
               end if
            end do
         end do

      end do

      jc = jcc

   end subroutine getNearestCell2dReg
!
!=======================================================================
!
   subroutine getOtherNearCellIndices2dReg( x, xmin2, xmax1, z, zmin2, zmax1,  &
!             ============================  i    i      i    i    i      i
                                            proj, posc, ffp, pmax, ncn, icnear )
!                                            i     i     i    i     o     o
!
!     Determine cells that are within the free flight path (ffp) cylinder (= cylinder with
!     height ffp and radius pmax).
!
!     Only cells in the range (1:nx,1:nz) are returned.
!
!     Note: this is very similar to getOtherNearCellsRegGrid2d
!
!  Arguments:
!
      real(REAL8), intent(in) :: x(0:)          ! positions of grid points
      real(REAL8), intent(in) :: xmin2          ! smallest x larger than x(0)
      real(REAL8), intent(in) :: xmax1          ! largest x smaller than x(nx)
      real(REAL8), intent(in) :: z(0:)          ! positions of grid points
      real(REAL8), intent(in) :: zmin2          ! smallest z larger than z(0)
      real(REAL8), intent(in) :: zmax1          ! largest z smaller than z(nz)
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
      integer(INT4) :: nx, nz, ix, iz, ixmin_ab, izmin_ab, ixmax_ab, izmax_ab
      real(REAL8) :: wx, dx, x1, x2, xmin, xmax, xisct, ffp_p
      real(REAL8) :: wz, dz, z1, z2, zmin, zmax, zisct
      real(REAL8) :: ax, az, bx, bz
      real(REAL8) :: rlam1, rlam2, dirx, diry, dirz, dirxy, dirxz, diryz, dirx_p, dirz_p
      logical, allocatable :: overlap(:,:)      ! 2D mask indicating cells overlaping with the
                                                !    ffp cylinder
      logical :: isInside, cut1, cut2
!
!  Start of Executable Code:
!  ------------------------
!
      nx = ubound( x, 1 )
      wx = x(1) - x(0)
      nz = ubound( z, 1 )
      wz = z(1) - z(0)

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

      ax = proj%pos(1) - dirx_p * pmax * abs(diry)
      az = proj%pos(3) - dirz_p * pmax * abs(diry)
      bx = posc(1) + dirx_p * pmax * abs(diry)
      bz = posc(3) + dirz_p * pmax * abs(diry)
!
!  Define the rectangle that contains the projection of the ffp cylinder and
!  the corresponding range of cell indices
!
      dx = pmax * diryz
      dz = pmax * dirxy

      x1 = proj%pos(1) - sign( dx, dirx )
      x2 = posc(1) + sign( dx, dirx )
      z1 = proj%pos(3) - sign( dz, dirz )
      z2 = posc(3) + sign( dz, dirz )

      xmin = min( x1, x2 )
      xmax = max( x1, x2 )
      zmin = min( z1, z2 )
      zmax = max( z1, z2 )

      xmin = max( xmin2, min( xmax1, xmin ) )     ! don't apply bounds
      xmax = max( xmin2, min( xmax1, xmax ) )     !    to ixmin,ixmax
      zmin = max( zmin2, min( zmax1, zmin ) )     !    because of possible
      zmax = max( zmin2, min( zmax1, zmax ) )     !    integer overflow

      ixmin = ceiling( ( xmin - x(0) ) / wx )  ! indices are valid cell indices in
      ixmax = ceiling( ( xmax - x(0) ) / wx )  !    the ranges [1..nx], [1..nz]
      izmin = ceiling( ( zmin - z(0) ) / wz )  !
      izmax = ceiling( ( zmax - z(0) ) / wz )  !
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

      xmin = max( xmin2, min( xmax1, xmin ) )
      xmax = max( xmin2, min( xmax1, xmax ) )
      zmin = max( zmin2, min( zmax1, zmin ) )
      zmax = max( zmin2, min( zmax1, zmax ) )

      ixmin_ab = ceiling( ( xmin - x(0) ) / wx )  ! indices must be in the ranges
      ixmax_ab = ceiling( ( xmax - x(0) ) / wx )  !    [ixmin...ixmax], [izmin...izmax]
      izmin_ab = ceiling( ( zmin - z(0) ) / wz )  !
      izmax_ab = ceiling( ( zmax - z(0) ) / wz )  !
!
!  Check for cell corners inside the ffp cylinder. All four adjoining cells overlap the
!  ffp cylinder
!
      ffp_p = ffp * dirxz                          ! p = projected to xz-plane

      do iz = izmin, izmax-1
         do ix = ixmin, ixmax-1

            isInside = .false.

            x1 = x(ix)
            z1 = z(iz)

            isInside = isInsideFFPcyl( [x1,z1], [proj%pos(1),proj%pos(3)], [dirx_p,dirz_p],  &
!                      ==============     i                 i                     i
                                       ffp_p, pmax, abs(diry) )
!                                        i     i        i

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

         x1 = x(0)
         x2 = x(nx)
         z1 = z(iz)
         z2 = z(iz)

         call mcisct( rlam1, rlam2, cut1, cut2, x1, z1, x2, z2, ax, az, bx, bz )
!             ======    o      o     o     o    i   i   i   i   i   i   i   i

         if( .not. cut2 ) cycle ! no valid intersection only if A^B outside the polygon
                                ! (can this happen?)

         xisct = x1 + ( x2 - x1 ) * rlam1       ! x-coordinate of intersection point

         ix = ceiling( ( xisct - x(0) ) / wx )
         ix = max( ixmin, min( ixmax, ix ) )    ! just in case AB intersects at
                                                !    cell area boundary
         overlap(ix,iz:iz+1) = .true.

      end do
!
!  Check intersection with the vertical (x=const) lines
!
      do ix = ixmin_ab, ixmax_ab-1

         x1 = x(ix)
         x2 = x(ix)
         z1 = z(0)
         z2 = z(nz)

         call mcisct( rlam1, rlam2, cut1, cut2, x1, z1, x2, z2, ax, az, bx, bz )
!             ======    o      o     o     o    i   i   i   i   i   i   i   i

         if( .not. cut2 )  cycle ! (?)

         zisct = z1 + ( z2 - z1 ) * rlam1       ! z-coordinate of intersection point

         iz = ceiling( ( zisct - z(0) ) / wz )
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
            call index2to1( ix, iz, nx, icnear(ncn) )
!                =========  i   i   i       o
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
!      ncn = size(icnear)
!

   end subroutine getOtherNearCellIndices2dReg
!
!=======================================================================
!
   subroutine fixNegDensCell2dReg( dens )
!             ===================   io
!
!     Eliminate negative cell densities by drawing material from neighboring cells.
!     Upon exit, densities less than zero on entry, are set to zero.
!
!     OpenMP: This subroutine must be called from the master thread
!
!  Arguments:
!
      real(REAL8), target, intent(inout) :: dens(:,:,:)  ! (nelem,nx, nz) densities of atoms
!
!  Local Variables:
!
      logical, allocatable :: mask(:)
      integer(INT4) :: ielem, nelem_, ix, nx, iz, nz, i, k, n
      integer(INT4), allocatable :: jx(:), jz(:)
      real(REAL8) :: sum_dens
!
!  Start of Executable Code:
!  ------------------------
!
      nelem_ = size( dens, 1 )
      nx = size( dens, 2 )
      nz = size( dens, 3 )
      if( nelem_ /= nelem )  error stop 'fixNegDensCell2dReg: Inconsistent nelem'

      do iz = 1, nz
         do ix = 1, nx
            do ielem = 1, nelem
               if( dens(ielem,ix,iz) >= ZERO )  cycle
               print '(a)', 'fixNegDensCell2dReg: dens(' // str( ielem, 'i3' ) // ',' //  &
                            str( ix, 'i9' ) // ',' // str( iz, 'i9' ) // ')=' //  &
                            str( dens(ielem,ix,iz), '1es10.3' )
               do k = 1, max( nx, nz )          ! go outwards from (ix,iz) in shells
                  jx = [ (i, i=ix-k,ix+k-1), (ix+k, i=1,2*k),  &
                         (i, i=ix+k,ix-k+1,-1), (ix-k, i=1,2*k) ]
                  jz = [ (iz+k, i=1,2*k), (i, i=iz+k,iz-k+1,-1),  &
                         (iz-k, i=1,2*k), (i, i=iz-k,iz+k-1) ]
                  mask = ( jx(:) >= 1 .and. jx(:) <= nx .and. jz(:) >= 1 .and. jz(:) <= nz )
                  if( .not.any( mask(:) ) ) then
                     print *, 'WARNING: Negative density could not completely be eliminated'
                     exit
                  end if

                  n = size( mask )
                  sum_dens = ZERO
                  do i = 1, n
                     ! avoid redistributing negative densities
                     if( mask(i) )  mask(i) = ( dens(ielem,jx(i),jz(i)) > ZERO )
                     if( mask(i) )  sum_dens = sum_dens + dens(ielem,jx(i),jz(i))
                  end do

                  if( sum_dens > -dens(ielem,ix,iz) ) then
                     do i = 1, n
                        if( mask(i) )  dens(ielem,jx(i),jz(i)) =  &
                              dens(ielem,jx(i),jz(i)) * (ONE + dens(ielem,ix,iz)/sum_dens)
                     end do
                     dens(ielem,ix,iz) = ZERO
                     exit
                  else
                     do i = 1, n
                        if( mask(i) )  dens(ielem,jx(i),jz(i)) = ZERO
                     end do
                     dens(ielem,ix,iz) = dens(ielem,ix,iz) + sum_dens
                  end if
               end do
            end do
         end do
      end do

   end subroutine fixNegDensCell2dReg
!
!=======================================================================
!
   subroutine relaxCell2dReg( x, z,  &
!             ==============  io io
                              dens, damdens, dens0, ffgeom, lxtal, dpa, pressure, stress,  &
!                              io     io      io      io     io    io      io       io
                              istep )
!                               i
!
!     Redistribute the cell densities so that stress relaxes.
!     A viscoelastic relaxation model is employed.
!
!     This subroutine should not be called too often in order to avoid artificial diffusion,
!     but often enough to keep the surface sufficiently accurate.
!
!     OpenMP: This subroutine must be called from the master thread
!
!  Arguments:
!
      real(REAL8), allocatable, intent(inout) :: x(:)          ! x positions of grid points
      real(REAL8), allocatable, intent(inout) :: z(:)          ! z positions of grid points
      real(REAL8), allocatable, intent(inout) :: dens(:,:,:)   ! (nelem,nx,nz) densities
                                                               !    on entry: before relax
                                                               !    on exit: after  relax
      real(REAL8), allocatable, intent(inout) :: damdens(:,:,:)! (nelemIn,nx,nz) damage dens.
      real(REAL8), allocatable, intent(inout) :: dens0(:,:,:)  ! (nelem,nx,nz) densities
                                                               !    after previous relax
      real(REAL8), allocatable, intent(inout) :: ffgeom(:,:)   ! geometric fill fraction:
                                                               !    on entry: after last relax
                                                               !    on exit: after this relax
      logical, allocatable, intent(inout) :: lxtal(:,:)        ! crystallinity flag
      real(REAL8), allocatable, intent(inout) :: dpa(:,:)      ! displacements per atom
      real(REAL8), allocatable, intent(inout) :: pressure(:,:) ! pressure
      real(REAL8), allocatable, intent(inout) :: stress(:,:,:,:)  ! stress tensor
      integer(INT4), intent(in) :: istep                       ! step counter
!
!  Local Variables:
!
      integer(INT4) :: ix, nx, iz, nz, nelemIn
      integer(INT4) :: naddx(2), naddz(2), nxnew, nznew, shift(2)
      real(REAL8) :: wx, wz
      real(REAL8), allocatable :: xnew(:), znew(:)
!
!  Start of Executable Code:
!  ------------------------
!
      nelemIn = 2*nelem
!
!  Add layer of empty cells at free boundaries to allow swelling
!
      nx = ubound( x, 1 )
      wx = x(1) - x(0)
      nz = ubound( z, 1 )
      wz = z(1) - z(0)

      where( geom%bndtypex(:) == ENUM%bndtype%free )
         naddx(:) = 1
      elsewhere
         naddx(:) = 0
      end where
      nxnew = nx + sum( naddx(:) )

      where( geom%bndtypez(:) == ENUM%bndtype%free )
         naddz(:) = 1
      elsewhere
         naddz(:) = 0
      end where
      nznew = nz + sum( naddz(:) )

      call alloc( xnew, 0, nxnew, 'xnew', 'relaxCell2dReg' )
!          =====   io   i   i       i             i
      xnew(naddx(1):nxnew-naddx(2)) = x(:)
      if( naddx(1) == 1 )  xnew(0) = xnew(1) - wx
      if( naddx(2) == 1 )  xnew(nxnew) = xnew(nxnew-1) + wx
      call move_alloc( xnew, x )

      call alloc( znew, 0, nznew, 'znew', 'relaxCell2dReg' )
!          =====   io   i   i       i             i
      znew(naddz(1):nznew-naddz(2)) = z(:)
      if( naddz(1) == 1 )  znew(0) = znew(1) - wz
      if( naddz(2) == 1 )  znew(nznew) = znew(nznew-1) + wz
      call move_alloc( znew, z )
!
!  Reallocate and shift dependent variables
!
      shift = [naddx(1), naddz(1)]
      call changeShape( ffgeom, [nxnew, nznew], 'ffgeom', 'relaxCell2dReg', shift_=shift )
!          ===========    io          i            i              i              i
      call changeShape( lxtal, [nxnew, nznew], 'lxtal', 'relaxCell2dReg', shift_=shift )
!          ===========   io          i            i             i              i
      call changeShape( dpa, [nxnew, nznew], 'dpa', 'relaxCell2dReg', shift_=shift )
!          ===========  io         i           i            i              i
      if( dam%relaxmodel == ENUM%relaxmodel%viscoelastic ) then
         call changeShape( pressure, [nxnew, nznew], 'pressure', 'relaxCell2dReg',  &
!             ===========     io           i             i               i
                           shift_=shift )
!                               i
         call changeShape( stress, [2, 2, nxnew, nznew], 'stress', 'relaxCell2dReg',  &
!             ===========    io          i             i              i
                           shift_=[0, 0, shift] )
!                                    i
      end if
      call changeShape( dens, [nelem, nxnew, nznew], 'dens', 'relaxCell2dReg',  &
!          ===========   io            i               i             i
                        shift_=[0, shift] )
!                               i
      if( out%lhis2 ) then
         call changeShape( damdens, [nelemIn, nxnew, nznew], 'damdens', 'relaxCell2dReg', &
!             ===========     io               i                 i              i
                           shift_=[0, shift] )
!                                  i
      end if
      call changeShape( dens0, [nelem, nxnew, nznew], 'dens0', 'relaxCell2dReg',  &
!          ===========   io             i                i             i
                        shift_=[0, shift] )
!                               i
!
!  Relax the fill fractions
!
#ifdef AMAT
      call relaxViscoElasticEuler2d( x(:), z(:), dens(:,:,:), damdens(:,:,:),  &
!          ========================   io    io       io           io
                                     dens0(:,:,:), ffgeom(:,:), dpa(:,:),  &
!                                        i             io          i
                                     pressure(:,:), stress(:,:,:,:), istep )
!                                         io              io           i
#else
      error stop 'ViscoElasticEuler2d model requires separate license'
#endif
!
!  Remove empty cells on both sides and add cells at free boundaries to accomodate vacuum layer
!
      nx = ubound( x, 1 )
      nz = ubound( z, 1 )

      do ix = 1, nx
         if( any( dens(:,ix,:) > ZERO ) )  exit
      end do
      naddx(1) = - (ix-1)
      if( geom%bndtypex(1) == ENUM%bndtype%free )  &
            naddx(1) = naddx(1) + ceiling( sn%pmax*(ONE+TWO*EPSMB) / wx )

      do ix = nx, 1, -1
         if( any( dens(:,ix,:) > ZERO ) )  exit
      end do
      naddx(2) = - (nx-ix)
      if( geom%bndtypex(2) == ENUM%bndtype%free )  &
            naddx(2) = naddx(2) + ceiling( sn%pmax*(ONE+TWO*EPSMB) / wx )
      nxnew = nx + sum( naddx(:) )

      call alloc( xnew, 0, nxnew, 'xnew', 'relaxCell2dReg' )
!          =====   io   i   i       i             i
      xnew(:) = [( x(0)+ix*wx, ix=-naddx(1),nx+naddx(2) )]
      call move_alloc( xnew, x )

      do iz = 1, nz
         if( any( dens(:,:,iz) > ZERO ) )  exit
      end do
      naddz(1) = - (iz-1)
      if( geom%bndtypez(1) == ENUM%bndtype%free )  &
            naddz(1) = naddz(1) + ceiling( sn%pmax*(ONE+TWO*EPSMB) / wz )

      do iz = nz, 1, -1
         if( any( dens(:,:,iz) > ZERO ) )  exit
      end do
      naddz(2) = - (nz-iz)
      if( geom%bndtypez(2) == ENUM%bndtype%free )  &
            naddz(2) = naddz(2) + ceiling( sn%pmax*(ONE+TWO*EPSMB) / wz )
      nznew = nz + sum( naddz(:) )

      call alloc( znew, 0, nznew, 'znew', 'relaxCell2dReg' )
!          =====   io   i   i       i             i
      znew(:) = [( z(0)+iz*wz, iz=-naddz(1),nz+naddz(2) )]
      call move_alloc( znew, z )
!
!  Reallocate and shift dependent variables (same code as above)
!
      shift = [naddx(1), naddz(1)]
      call changeShape( ffgeom, [nxnew, nznew], 'ffgeom', 'relaxCell2dReg', shift_=shift )
!          ===========    io          i            i              i              i
      call changeShape( lxtal, [nxnew, nznew], 'lxtal', 'relaxCell2dReg', shift_=shift )
!          ===========   io          i            i             i              i
      call changeShape( dpa, [nxnew, nznew], 'dpa', 'relaxCell2dReg', shift_=shift )
!          ===========  io         i           i            i              i
      if( dam%relaxmodel == ENUM%relaxmodel%viscoelastic ) then
         call changeShape( pressure, [nxnew, nznew], 'pressure', 'relaxCell2dReg',  &
!             ===========     io           i             i               i
                           shift_=shift )
!                               i
         call changeShape( stress, [2, 2, nxnew, nznew], 'stress', 'relaxCell2dReg',  &
!             ===========    io          i             i              i
                           shift_=[0, 0, shift] )
!                                    i
      end if
      call changeShape( dens, [nelem, nxnew, nznew], 'dens', 'relaxCell2dReg',  &
!          ===========   io            i               i             i
                        shift_=[0, shift] )
!                               i
      if( out%lhis2 ) then
         call changeShape( damdens, [nelemIn, nxnew, nznew], 'damdens', 'relaxCell2dReg', &
!             ===========     io               i                 i              i
                           shift_=[0, shift] )
!                                  i
      end if
      call changeShape( dens0, [nelem, nxnew, nznew], 'dens0', 'relaxCell2dReg',  &
!          ===========   io             i                i             i
                        shift_=[0, shift] )
!                               i
!                                    i
!  Update surface
!
#ifdef AMAT
      call initSurface2dFromFF( x(:), z(:), ffgeom(:,:), dens(:,:,:), istep )
!          ===================   i     i         i           i          i
#else
         error stop 'ViscoElasticEuler2d model requires separate license'
#endif

   end subroutine relaxCell2dReg
!
!=======================================================================
!
   subroutine writeCell2dReg( x, z, dens, pressure, stress )
!             ==============  i  i   i       i        i
!
!     Write the contents of the target cells to the CELL2 file (same format as HIS2 file)
!
!  Arguments:
!
      real(REAL8), intent(in) :: x(0:)          ! (0:nx) x coordinates of grid points
      real(REAL8), intent(in) :: z(0:)          ! (0:nz) y coordinates of grid points
      real(REAL8), intent(in) :: dens(:,:,:)    ! (nelem,nx,nz) densities
      real(REAL8), intent(in) :: pressure(:,:)  ! (nx,nz) pressure
      real(REAL8), intent(in) :: stress(:,:,:,:)! (2,2,nx,nz) zz component of stress tensor
!
!  Local Variables:
!
      type(His2d_type) :: his
      integer(INT4) :: nx, nz, iz
      real(REAL8) :: wx, wz
!
      character(6), save :: fmode = 'write'
!
!  Start of Executable Code:
!  ------------------------
!
      nx = ubound( x, 1 )
      wx = x(1) - x(0)
      nz = ubound( z, 1 )
      wz = z(1) - z(0)

      his%nbox(:) = [nx,nz]
      his%imin(:) = [1,1]
      his%imax(:) = [nx,nz]
      his%hbeg(:) = [x(0),z(0)]
      his%wbox(:) = [wx,wz]

      if( dam%relaxmodel == ENUM%relaxmodel%viscoelastic ) then
         call alloc( his%hist, [nelem+5,nx,nz], 'his%hist', 'writeCell1dReg' )
!             =====     io            i             i               i
         do iz = 1, nz
            his%hist(1:nelem,:,:) = dens(:,:,:)
            his%hist(nelem+1,:,:) = pressure(:,:) / 1e24     ! compensate for scaling
            his%hist(nelem+2,:,:) = stress(1,1,:,:) / 1e24   ! compensate for scaling
            his%hist(nelem+3,:,:) = stress(2,1,:,:) / 1e24   ! compensate for scaling
            his%hist(nelem+4,:,:) = stress(1,2,:,:) / 1e24   ! compensate for scaling
            his%hist(nelem+5,:,:) = stress(2,2,:,:) / 1e24   ! compensate for scaling
         end do

         call alloc( his%indhis, nelem+5, 'his%indhis', 'writeCell1dReg' )
!             =====      io         i           i               i
         his%indhis(1:nelem) = elem(1:nelem)%z
         his%indhis(nelem+1:nelem+5) = 0
      else
         call alloc( his%hist, [nelem,nx,nz], 'his%hist', 'writeCell1dReg' )
!             =====     io           i            i               i
         his%hist(:,:,:) = dens(:,:,:)

         call alloc( his%indhis, nelem, 'his%indhis', 'writeCell1dReg' )
!             =====      io        i          i               i
         his%indhis(:) = elem(1:nelem)%z
      end if
!
      call writeHis2d( 'CELL2', fmode, his, -BIG, factor_=1e24_REAL8 )
!          ==========     i       i     i     i          i

      fmode = 'append'

   end subroutine writeCell2dReg
!
!=======================================================================
!
   end module Target_Cell2d_Reg_mod
