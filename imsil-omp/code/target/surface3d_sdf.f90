module Target_Surface3d_SDF_mod
!
!  Define the surface by a signed distance function (SDF) defined on a rectilinear grid
!
!  Used Modules:
!  ============
!
   use Const_mod,             only: INT4, INT8, REAL8,  &
                                    ZERO, HALF, ONE, TWO, SMALL, BIG, EPSMB
   use Array_mod,             only: alloc
   use Open_mod,              only: openFile

   use Math_BinSearch_mod,    only: binSearch
   use Math_IntMath_mod,      only: index3to1
   use Math_Table_mod,        only: Table3d_type, initTable, interpolate

   use Param_Setup_mod,       only: ndim
   use Param_Geom_mod,        only: geom
   use Param_Ion_mod,         only: ion
   use Param_Out_mod,         only: out

   use Stat_His3d_mod,        only: His3d_type, writeHis3d
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
   public :: initSurface3dSDF                ! initialize grid and SDF
!            ================
   public :: getNcellSurface3dSDF            ! get the number of cells
!            ====================
   public :: getCellIndex3dSDF               ! get cell index at given position
!            =================
   public :: getDistSurface3dSDF             ! get the signed distance of a point
!            ===================             !    from the surface
   public :: writeSurface3d                  ! write the surface polygon to the output file
!            ==============                  !    as constructed from the distance function
!
!  Private Procedures:
!  ==================
!
   private :: writeDist                      ! write the distance function to file DIST3
!             =========
!
!  Private Data:
!  ============
!
   type :: Grid_type                         ! distance function grid
!          =========
      logical :: regular                     ! indicates constant cell widths
      integer(INT4) :: nx                    ! number of cells in x direction
      integer(INT4) :: ny                    ! number of cells in y direction
      integer(INT4) :: nz                    ! number of cells in z direction
      real(REAL8), allocatable :: wx(:)      ! (nx) cell width in x direction
      real(REAL8), allocatable :: wy(:)      ! (nx) cell width in y direction
      real(REAL8), allocatable :: wz(:)      ! (nz) cell width in z direction
      real(REAL8), allocatable :: x(:)       ! (0:nx) positions of x interfaces between cells
      real(REAL8), allocatable :: y(:)       ! (0:nx) positions of y interfaces between cells
      real(REAL8), allocatable :: z(:)       ! (0:nz) positions of z interfaces between cells
   end type Grid_type

   type(Grid_type), public :: grid

   type(Table3d_type) :: sdf_tab             ! table of signed distance function
!
!=======================================================================
!
contains
!
!=======================================================================
!
   subroutine initSurface3dSDF( xlim, ylim, zlim )
!             ================   o     o     o
!
!     Initialize the simulation domain and the signed distance function from the surface
!
!     The rectilinear grid must cover all SDF > -sn%pmax*(1+2*epsB)
!
!  Arguments:
!
      real(REAL8), intent(out) :: xlim(2)       ! minimum and maximum x values to be queried
      real(REAL8), intent(out) :: ylim(2)       ! minimum and maximum y values to be queried
      real(REAL8), intent(out) :: zlim(2)       ! minimum and maximum z values to be queried
!
!  Local Variables:
!
      character(256) :: line
      integer(INT4) :: i, nd, istat, ix, iy, iz
      real(REAL8), allocatable :: sdf(:,:,:)
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
      grid%ny = 0
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
            call alloc( grid%wx, grid%nx, 'grid%wx', 'initSurface3dSDF' )
!                =====     o        i         i               i
            grid%wx(:) = ZERO
         else if( line(:i-1) == 'nY' ) then
            read( line(i+1:), '(i99)' )  grid%ny
            grid%ny = grid%ny - 1
            call alloc( grid%wy, grid%ny, 'grid%wy', 'initSurface3dSDF' )
!                =====     o        i         i               i
            grid%wy(:) = ZERO
         else if( line(:i-1) == 'nZ' ) then
            read( line(i+1:), '(i99)' )  grid%nz
            grid%nz = grid%nz - 1
            call alloc( grid%wz, grid%nz, 'grid%wz', 'initSurface3dSDF' )
!                =====     o        i         i               i
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
            if( .not.allocated( grid%wy ) )  goto 912
            grid%regular = ( index( trim( adjustl( line(i+1:) ) ), ' ' ) == 0 .and.  &
                             index( trim( adjustl( line(i+1:) ) ), ',' ) == 0 )
            if( grid%regular ) then
               read( line(i+1:), '(f99.0)' )  grid%wy(1)
               grid%wy(:) = grid%wy(1)
            else
               backspace( 11 )
               read( 11, *, iostat=istat )  line, grid%wy(:)
               if( istat /= 0 )  goto 913
            end if
         else if( line(:i-1) == 'dZ' ) then
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

      if( grid%nx == 0 .or. grid%ny == 0 .or. grid%nz == 0 ) goto 919
      if( any( grid%wx(:) == ZERO ) .or. any( grid%wy(:) == ZERO ) .or.  &
          any( grid%wz(:) == ZERO ) )  goto 919
!
!  Grid positions
!
      call alloc( grid%x, 0, grid%nx, 'grid%nx', 'initSurface3dSDF' )
!          =====    o     i     i         i              i
      grid%x(0) = - HALF * sum( grid%wx(:) )
      do ix = 1, grid%nx
         grid%x(ix) = grid%x(ix-1) + grid%wx(ix)
      end do

      call alloc( grid%y, 0, grid%ny, 'grid%ny', 'initSurface3dSDF' )
!          =====    o     i     i         i              i
      grid%y(0) = - HALF * sum( grid%wy(:) )
      do iy = 1, grid%ny
         grid%y(iy) = grid%y(iy-1) + grid%wy(iy)
      end do

      call alloc( grid%z, 0, grid%nz, 'grid%nz', 'initSurface3dSDF' )
!          =====    o     i     i         i              i
      grid%z(0) = - HALF * sum( grid%wz(:) )
      do iz = 1, grid%nz
         grid%z(iz) = grid%z(iz-1) + grid%wz(iz)
      end do
!
!  SDF values
!
      call alloc( sdf, [0, 0, 0], [grid%nx, grid%ny, grid%nz], 'sdf', 'initSurface3dSDF' )
!          =====   o       i                   i                 i            i
      do ix = 0, grid%nx
         do iy = 0, grid%ny
            read( 11, * )  ( sdf(ix,iy,iz), iz=0,grid%nz )
         end do
      end do

      close( 11 )
!
!  Initialize SDF table
!
      call initTable( sdf_tab, grid%x(:), grid%y(:), grid%z(:), sdf(:,:,:) )
!          =========     o        i          i          i           i
!
!  Bounding box (we rely on free boundaries already containing a vacuum layer)
!
      xlim(1) = grid%x(0)
      xlim(2) = grid%x(grid%nx)
      ylim(1) = grid%y(0)
      ylim(2) = grid%y(grid%ny)
      zlim(1) = grid%z(0)
      zlim(2) = grid%z(grid%nz)

      if( out%ldist )  call writeDist( )
!                           =========

      return

900   print *, 'initSurface3dSDF: Error reading file ''', geom%sdffile(0), ''''
      error stop

910   error stop 'initSurface3dSDF: nDim=2 required on sdffile'
911   error stop 'initSurface3dSDF: Units=Angstrom required on sdffile'
912   error stop 'initSurface3dSDF: nX must be specified before dX in sdffile'
913   error stop 'initSurface3dSDF: error when reading dX record of sdffile'
914   error stop 'initSurface3dSDF: nY must be specified before dY in sdffile'
915   error stop 'initSurface3dSDF: error when reading dY record of sdffile'
919   error stop 'initSurface3dSDF: Incomplete grid definition'

   end subroutine initSurface3dSDF
!
!=======================================================================
!
   subroutine getNcellSurface3dSDF( ncell )
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
      ncell = int( grid%nx, INT8 ) * int( grid%ny, INT8 ) * int( grid%nz, INT8 )

   end subroutine getNcellSurface3dSDF
!
!=======================================================================
!
   subroutine getCellIndex3dSDF( x, y, z, ic, outside, what_ )
!             =================  i  i  i  o      o       i
!
!     Get the cell index at position [x,z]
!
!  Arguments:
!
      real(REAL8), intent(in) :: x              ! x coordinate of point
      real(REAL8), intent(in) :: y              ! y coordinate of point
      real(REAL8), intent(in) :: z              ! z coordinate of point
      integer(INT8), intent(out) :: ic          ! cell index
      logical, intent(out) :: outside           ! flag indicating point is outside grid
      character(*), intent(in), optional :: what_  ! indication what point is checked
!
!  Local Variables:
!
      integer(INT4) :: ix, iy, iz
      real(REAL8) :: xtmp, ytmp, ztmp
!
!  Start of Executable Code:
!  ------------------------
!
      if( grid%regular ) then
         xtmp = max( grid%x(0)-HALF*grid%wx(1), min( grid%x(grid%nx)+HALF*grid%wx(1), x ) )
         ytmp = max( grid%y(0)-HALF*grid%wy(1), min( grid%y(grid%nx)+HALF*grid%wy(1), y ) )
         ztmp = max( grid%z(0)-HALF*grid%wz(1), min( grid%z(grid%nz)+HALF*grid%wz(1), z ) )
         ix = ceiling( ( xtmp - grid%x(0) ) / grid%wx(1) )  ! ix can be [0,nx+1]
         iy = ceiling( ( ytmp - grid%y(0) ) / grid%wy(1) )  ! iy can be [0,ny+1]
         iz = ceiling( ( ztmp - grid%z(0) ) / grid%wz(1) )  ! iz can be [0,nz+1]
         outside = ( any( [ix,iy,iz] == 0 ) .or.  &
                     any( [ix,iy,iz] == [grid%nx,grid%ny,grid%nz]+1 ) )
         if( outside ) then
            ix = max( 1, min( grid%nx, ix ) )
            iy = max( 1, min( grid%ny, iy ) )
            iz = max( 1, min( grid%nz, iz ) )
         end if
      else
         call binSearch( grid%nx+1, grid%x(0:), x, ix )
!             =========      i          i       i  o
         call binSearch( grid%ny+1, grid%y(0:), y, iy )
!             =========      i          i       i  o
         call binSearch( grid%nz+1, grid%z(0:), z, iz )
!             =========      i          i       i  o
         outside = ( x < grid%x(0) .or. x > grid%x(grid%nx) .or.  &
                     y < grid%y(0) .or. y > grid%y(grid%ny) .or.  &
                     z < grid%z(0) .or. z > grid%z(grid%nz) )
      end if

      call index3to1( ix, iy, iz, grid%nx, grid%ny, ic )
!          =========  i   i   i      i        i     o

   end subroutine getCellIndex3dSDF
!
!=======================================================================
!
   subroutine getDistSurface3dSDF( x, y, z, sdf, dirx_, diry_, dirz_ )
!             ===================  i  i  i   o     o      o      o
!
!     Get the signed distance of a point pos=(x,y,z) from the surface by interpolation
!     ( + inside, - outside ) and
!     the direction towards/from the surface (pointing outwards the target)
!
!     NOTE: If the point is outside the grid, dist is set to -BIG, and dir_ is undefined.
!
!  Arguments:
!
      real(REAL8), intent(in) :: x                 ! x coordinate of point
      real(REAL8), intent(in) :: y                 ! y coordinate of point
      real(REAL8), intent(in) :: z                 ! z coordinate of point
      real(REAL8), intent(out) :: sdf              ! signed distance: (+) inside, (-) outside
      real(REAL8), optional, intent(out) :: dirx_  ! x/y/z component of direction towards/from
      real(REAL8), optional, intent(out) :: diry_  !    surface
      real(REAL8), optional, intent(out) :: dirz_  !
!
!  Local Variables:
!
!     None
!
!  Start of Executable Code:
!  ------------------------
!
!  Bilinear interpolation
!
      call interpolate( sdf_tab, x, y, z, sdf, dirx_, diry_, dirz_ )
!          ===========     i     i  i  i   o     o      o      o

   end subroutine getDistSurface3dSDF
!
!=======================================================================
!
   subroutine writeSurface3d( iion )
!             ==============    i
!
!     Write the contour as a set of closed surfaces (connected components) to the SURF file.
!
!  Arguments:
!
      integer(INT4), intent(in) :: iion         ! index of ion
!
!  Local Variables:
!
!     None
!
!  Start of Executable Code:
!  ------------------------
!
      error stop 'writeSurface3d: Not yet implemented'

   end subroutine writeSurface3d
!
!=======================================================================
!
   subroutine writeDist( )
!             =========
!
!     Write the distance from the surface to the DIST3 file (same format as HIS3 file)
!
!  Arguments:
!
!     None
!
!  Local Variables:
!
      integer(INT4) :: nx, ny, nz
      type(His3d_type) :: his
      real(REAL8) :: wx, wy, wz

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
      ny = size( sdf_tab%y )
      nz = size( sdf_tab%z )
      his%nbox(:) = [nx, ny, nz]
      his%imin(:) = 1
      his%imax(:) = [nx, ny, nz]
      wx = sdf_tab%wx
      wy = sdf_tab%wy
      wz = sdf_tab%wz
      his%hbeg(:) = [grid%x(0)-HALF*grid%wx(1), grid%y(0)-HALF*grid%wy(1),  &
                     grid%z(0)-HALF*grid%wz(1)]
      his%wbox(:) = [wx, wy, wz]

      allocate( his%hist( 1, grid%nx+1, grid%ny+1, grid%nz+1 ) )
      his%hist(1,:,:,:) = sdf_tab%a000(:,:,:)

      allocate( his%indhis( 1 ) )
      his%indhis(:) = 0

      call writeHis3d( 'DIST3', fmode, his )
!          ==========     i       i     i

      deallocate( his%hist )
      deallocate( his%indhis )

      fmode = 'append'

   end subroutine writeDist
!
!=======================================================================
!
end module Target_Surface3d_SDF_mod
