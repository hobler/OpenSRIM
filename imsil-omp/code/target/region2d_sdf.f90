module Target_Region2d_SDF_mod
!
!  Geometry operations on 2D regions defined by signed distance functions (SDFs).
!
!  Reads in the SDFs from files.
!
!  Used Modules:
!  ============
!
   use Const_mod,                only: INT4, INT8, REAL8,  &
                                       ZERO, ONE, TWO, HALF, QUART,  &
                                       EPSMA, EPSMB, SMALL, BIG,  &
                                       ENUM
   use Array_mod,                only: alloc, check_allocation
   use Open_mod,                 only: openFile
   use Types_mod,                only: Proj_type

   use Math_BinSearch_mod,       only: binSearch
   use Math_IntMath_mod,         only: index1to2
   use Math_Table_mod,           only: Table2d_type, initTable, interpolate

   use Param_Setup_mod,          only: nr, ndim, nelem
   use Param_Geom_mod,           only: geom

   use Target_Surface2d_SDF_mod, only: grid
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
   public :: initRegion2dSDF                 ! read the signed distance functions
!            ===============
   public :: getDist2FromRegion2dSDF         ! determine square of distance from region
!            =======================
   public :: getDist2RectRegion2dSDF         ! determine distance between rectangle and
!            =======================         !    region
   public :: isOutsideRegion2dSDF            ! determine if point is outside given region
!            ====================
!
!  Private Data:
!  ============
!
!  Signed distance function
!
   type(Table2d_type), allocatable :: sdf_tab(:) ! (nr) signed distance function
!
!=======================================================================
!
contains
!
!=======================================================================
!
   subroutine initRegion2dSDF( )
!             ===============
!
!     Read SDFs and initialize the signed distance functions
!
!  Arguments:
!
!     None
!
!  Local Variables:
!
      logical :: regular
      character(256) :: line, errmsg
      integer(INT4) :: i, ir, nd, istat, ix, iz, nx, nz
      real(REAL8), allocatable :: wx(:), wz(:)
      real(REAL8), allocatable :: sdf(:,:)
!
!  Start of Executable Code:
!  ------------------------
!
      allocate( sdf_tab( nr ), stat=istat, errmsg=errmsg )
      call check_allocation( 'initRegion2dSDF', 'sdf_tab', istat, errmsg )

      do ir = 1, nr

         call openFile( 11, geom%sdffile(ir), 'formatted', 'read' )
!             ========   i         i               i         i
!
!  Read header
!
         nx = 0
         nz = 0

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
               read( line(i+1:), '(i99)' )  nx
               nx = nx - 1
               call alloc( wx, nx, 'wx', 'initRegion2dSDF')
!                   =====  o   i    i            i
               wx(:) = ZERO
            else if( line(:i-1) == 'nY' ) then
               read( line(i+1:), '(i99)' )  nz
               nz = nz - 1
               call alloc( wz, nz, 'wz', 'initRegion2dSDF' )
!                   =====  o   i    i            i
               wz(:) = ZERO
            else if( line(:i-1) == 'dX' ) then
               if( .not.allocated( wx ) )  goto 912
               regular = ( index( trim( adjustl( line(i+1:) ) ), ' ' ) == 0 .and.  &
                           index( trim( adjustl( line(i+1:) ) ), ',' ) == 0 )
               if( regular ) then
                  read( line(i+1:), '(f99.0)' )  wx(1)
                  wx(:) = wx(1)
               else
                  backspace( 11 )
                  read( 11, *, iostat=istat )  line, wx(:)
                  if( istat /= 0 )  goto 913
               end if
            else if( line(:i-1) == 'dY' ) then
               if( .not.allocated( wz ) )  goto 914
               regular = ( index( trim( adjustl( line(i+1:) ) ), ' ' ) == 0 .and.  &
                           index( trim( adjustl( line(i+1:) ) ), ',' ) == 0 )
               if( regular ) then
                  read( line(i+1:), '(f99.0)' )  wz(1)
                  wz(:) = wz(1)
               else
                  backspace( 11 )
                  read( 11, *, iostat=istat )  line, wz(:)
                  if( istat /= 0 )  goto 915
               end if
            else if( line(:i-1) == 'Values' ) then
               exit
            end if
         end do

         if( nx == 0 .or. nz == 0 )  goto 919
         if( any( wx(:) == ZERO ) .or. any( wz(:) == ZERO ) )  goto 919
!
!  Check grid positions
!
         if( nx /= grid%nx .or. nz /= grid%nz .or.  &
             any( abs( wx(:) - grid%wx(:) ) > EPSMB * abs( grid%wx(:) ) ) .or.  &
             any( abs( wz(:) - grid%wz(:) ) > EPSMB * abs( grid%wz(:) ) ) )  goto 930
         deallocate( wx, wz )
!
!  SDF values
!
         call alloc( sdf, [0, 0], [grid%nx, grid%nz], 'sdf', 'initRegion2dSDF' )
!             =====   o     i             i             i            i
         do ix = 0, grid%nx
            read( 11, * )  ( sdf(ix,iz), iz=0,grid%nz )
         end do

         close( 11 )

         call initTable( sdf_tab(ir), grid%x(:), grid%z(:), sdf(:,:) )
!             =========      o            i          i          i

         deallocate( sdf )

      end do

      return

900   print *, 'initRegion2dSDF: Error reading file ''', geom%sdffile(ir), ''''
      error stop

910   error stop 'initRegion2dSDF: nDim=2 required on sdffile'
911   error stop 'initRegion2dSDF: Units=Angstrom required on sdffile'
912   error stop 'initRegion2dSDF: nX must be specified before dX in sdffile'
913   error stop 'initRegion2dSDF: error when reading dX record of sdffile'
914   error stop 'initRegion2dSDF: nY must be specified before dY in sdffile'
915   error stop 'initRegion2dSDF: error when reading dY record of sdffile'
919   error stop 'initRegion2dSDF: Incomplete grid definition'

930   error stop 'initRegion2dSDF: Identical grids on all sdffiles required'

   end subroutine initRegion2dSDF
!
!=======================================================================
!
   subroutine getDist2FromRegion2dSDF( x, z, ir, dist2 )
!             =======================  i  i  i     o
!
!     Determine the square of the distance of the point [x,z] from region ir
!
!  Arguments:
!
      real(REAL8), intent(in) :: x               ! x coordinate of point
      real(REAL8), intent(in) :: z               ! z coordinate of point
      integer(INT4), intent(in) :: ir           ! index of region
      real(REAL8), intent(out) :: dist2         ! square of the distance from region
!
!  Local Variables:
!
      real(REAL8) :: sdf_val
!
!  Start of Executable Code:
!  ------------------------
!
      call interpolate( sdf_tab(ir), x, z, sdf_val )
!          ===========      i        i  i     o
      dist2 = sdf_val**2

   end subroutine getDist2FromRegion2dSDF
!
!=======================================================================
!
   subroutine getDist2RectRegion2dSDF( ic, ir, cinside, dist2 )
!             =======================  i   i      i       o
!
!     Determine the minimum distance of any point inside a rectangle from any point inside
!     region ir.
!
!  Arguments:
!
      integer(INT8), intent(in) :: ic           ! cell index
      integer(INT4), intent(in) :: ir           ! index of region
      logical, intent(out) :: cinside           ! flag indicating rectangle center is inside
                                                !    region
      real(REAL8), intent(out) :: dist2         ! square of distance betweeen ractangle and
                                                !    region
!
!  Local Variables:
!
      logical :: all_outside
      integer(INT4) :: i, j, ii, jj, i1, i2, j1, j2, ix, iz
      real(REAL8) :: x(2), z(2), xx(2), zz(2), xc, zc, sdf_val
!
!  Start of Executable Code:
!  ------------------------
!
      call index1to2( ic, grid%nx, ix, iz )
!          =========  i      i     o   o

      x(:) = grid%x(ix-1:ix)
      z(:) = grid%z(iz-1:iz)
!
!  Check cell center (this will often cause return from subroutine)
!
      xc = HALF*sum( x(:) )
      zc = HALF*sum( z(:) )
      call interpolate( sdf_tab(ir), xc, zc, sdf_val )
!          ===========      i        i   i      o
      cinside = (sdf_val > ZERO)

      if( cinside ) then
         dist2 = ZERO
         return
      end if
!
!  Check cell corners considering periodic images of region
!  (Note: extrema of sdf are in the corners with linear interpolation)
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

      dist2 = BIG

      do ii = i1, i2
         do jj = j1, j2
            xx(:) = x(:) + ii * (geom%xlim(2) - geom%xlim(1))
            zz(:) = z(:) + jj * (geom%zlim(2) - geom%zlim(1))
            all_outside = .true.
            do i = 1, 2
               do j = 1, 2
                  call interpolate( sdf_tab(ir), xx(i), zz(j), sdf_val )
   !                   ===========      i          i      i       o
                  all_outside = (all_outside .and. sdf_val < ZERO )
                  dist2 = min( dist2, sdf_val**2 )
               end do
            end do
                                          ! set distmin2 to zero if cell intersects
                                          ! boundary (Note: corners may not be all inside)
            if( .not.all_outside ) then
               dist2 = ZERO
               return
            end if
         end do
      end do

   end subroutine getDist2RectRegion2dSDF
!
!=======================================================================
!
   subroutine isOutsideRegion2dSDF( x, z, ir, outside )
!             ====================  i  i  i      o
!
!     Determine the region the point [x,z] is located outside region ir
!
!  Arguments:
!
      real(REAL8), intent(in) :: x               ! x coordinate of point
      real(REAL8), intent(in) :: z               ! z coordinate of point
      integer(INT4), intent(in) :: ir           ! index of region
      logical, intent(out) :: outside           ! flag indicating point is outside region
!
!  Local Variables:
!
      real(REAL8) :: sdf_val
!
!  Start of Executable Code:
!  ------------------------
!
      call interpolate( sdf_tab(ir), x, z, sdf_val )
!          ===========      i        i  i     o
      outside = ( sdf_val < ZERO )

   end subroutine isOutsideRegion2dSDF
!
!=======================================================================
!
end module Target_Region2d_SDF_mod
