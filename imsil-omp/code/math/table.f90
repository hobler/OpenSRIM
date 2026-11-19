module Math_Table_mod
!
!  Defines 1D, 2D, and 3D tables and provides interpolation subroutines
!
!  Used Modules:
!  ============
!
   use Const_mod,          only: INT4, REAL8,  &
                                 ZERO, EPSMA
   use Array_mod,          only: check_allocation
   use Math_BinSearch_mod, only: binSearch
   use Math_Solve_mod,     only: tridiag
!
!  Default Type and Access:
!  =======================
!
   implicit none
   private
   save
!
!  Public Procedures
!  =================
!
   public :: initTable                       ! initialize the table
!            =========
   interface initTable
      module procedure initTable1d
      module procedure initTable2d
      module procedure initTable3d
   end interface

   public :: interpolate                     ! interpolate in the table
!            ===========
   interface interpolate
      module procedure interpolate1d
      module procedure interpolate2d
      module procedure interpolate3d
   end interface

   public :: getDerivatives                ! get the derivatives of 1D tabulated data
!            ==============
!
!  Public Data Types
!  =================
!
   type, public :: Table1d_type
!                  ============
      logical :: regular                     ! indicates equidistant abscissa
      real(REAL8), allocatable :: x(:)       ! abscissa
      real(REAL8) :: wx                      ! x(2)-x(1) (only used if regular=t)
      real(REAL8), allocatable :: a0(:)      ! table values
      real(REAL8), allocatable :: a1(:)      ! derivatives times intervals
      real(REAL8), allocatable :: a2(:)      ! 2nd order coefficient for interpolation
      real(REAL8), allocatable :: a3(:)      ! 3rd order coefficient for interpolation
   end type Table1d_type

   type, public :: Table2d_type
!                  ============
      logical :: regular                     ! indicates equidistant x and y
      real(REAL8), allocatable :: x(:)       ! independent variable
      real(REAL8), allocatable :: y(:)       ! independent variable
      real(REAL8) :: wx                      ! x(2)-x(1) (only used if regular=t)
      real(REAL8) :: wy                      ! y(2)-y(1) (only used if regular=t)
      real(REAL8), allocatable :: a00(:,:)   ! table values
      real(REAL8), allocatable :: a10(:,:)   ! 1st order coefficient for interpolation
      real(REAL8), allocatable :: a01(:,:)   ! 1st order coefficient for interpolation
      real(REAL8), allocatable :: a20(:,:)   ! 2nd order coefficient for interpolation
      real(REAL8), allocatable :: a11(:,:)   ! 2nd order coefficient for interpolation
      real(REAL8), allocatable :: a02(:,:)   ! 2nd order coefficient for interpolation
      real(REAL8), allocatable :: a30(:,:)   ! 3rd order coefficient for interpolation
      real(REAL8), allocatable :: a21(:,:)   ! 3rd order coefficient for interpolation
      real(REAL8), allocatable :: a12(:,:)   ! 3rd order coefficient for interpolation
      real(REAL8), allocatable :: a03(:,:)   ! 3rd order coefficient for interpolation
      real(REAL8), allocatable :: a31(:,:)   ! 4th order coefficient for interpolation
      real(REAL8), allocatable :: a13(:,:)   ! 4th order coefficient for interpolation
   end type Table2d_type

   type, public :: Table3d_type
!                  ============
      logical :: regular                        ! indicates equidistant x and y
      real(REAL8), allocatable :: x(:)          ! independent variable
      real(REAL8), allocatable :: y(:)          ! independent variable
      real(REAL8), allocatable :: z(:)          ! independent variable
      real(REAL8) :: wx                         ! x(2)-x(1) (used only if regular=t)
      real(REAL8) :: wy                         ! y(2)-y(1) (used only if regular=t)
      real(REAL8) :: wz                         ! z(2)-z(1) (used only if regular=t)
      real(REAL8), allocatable :: a000(:,:,:)   ! table values
      real(REAL8), allocatable :: a100(:,:,:)   ! 1st order coefficient for interpolation
      real(REAL8), allocatable :: a010(:,:,:)   ! 1st order coefficient for interpolation
      real(REAL8), allocatable :: a001(:,:,:)   ! 1st order coefficient for interpolation
      real(REAL8), allocatable :: a110(:,:,:)   ! 2nd order coefficient for interpolation
      real(REAL8), allocatable :: a101(:,:,:)   ! 2nd order coefficient for interpolation
      real(REAL8), allocatable :: a011(:,:,:)   ! 2nd order coefficient for interpolation
      real(REAL8), allocatable :: a111(:,:,:)   ! 3rd order coefficient for interpolation
   end type Table3d_type
!
!=======================================================================
!
contains
!
!=======================================================================
!
   subroutine initTable1d( table, x, val, dvaldx_, cubic_ )
!             ===========    o    i   i      i       i
!
!     Initialize the table from values and optionally derivatives
!
!     Depending on cubic_, linear or cubic interpolation is prepared. If cubic_=.true. and
!     derivatives are not specified, they are calculated using getDerivatives.
!
!  Arguments:
!
      type(Table1d_type), intent(out) :: table        ! table
      real(REAL8), intent(in) :: x(:)                 ! x values
      real(REAL8), intent(in) :: val(:)               ! values
      real(REAL8), intent(in), optional :: dvaldx_(:) ! d(val)/dx
      logical, intent(in), optional :: cubic_         ! type of interpolation (default: .true.)
!
!  Local Variables:
!
      logical :: cubic
      character(256) :: errmsg
      integer(INT4) :: nx, istat
      real(REAL8), allocatable :: dx(:)
      real(REAL8), allocatable :: f0(:), f1(:)
      real(REAL8), allocatable :: fx0(:), fx1(:)

      real(REAL8), allocatable :: dvaldx_auto(:)
      target :: dvaldx_, dvaldx_auto
      real(REAL8), pointer :: dvaldx_ptr(:) => null()
!
!  Start of Executable Code:
!  ------------------------
!
!  Check the input variables
!
      !print *, 'Entering initTable1d...'
      !print *, 'x =', x
      !print *, 'val =', val
      !print *, 'dvaldx =', dvaldx_

      nx = size( x )

      if( size( val ) /= nx )  &
            error stop 'initTable1d: Shape of val is inconsistent with size of axes'

      if( present( dvaldx_ ) ) then
         if( size( dvaldx_ ) /= nx )  &
               error stop 'initTable1d: Size of dvaldx must match size of val'
      end if

      if( present( cubic_ ) ) then
         cubic = cubic_
      else
         cubic = .true.
      end if
!
!  Automatic calculation of derivatives
!
      if( cubic ) then

         if( present( dvaldx_ ) ) then

            dvaldx_ptr => dvaldx_

         else

            allocate( dvaldx_auto( nx ), stat=istat, errmsg=errmsg )
            call check_allocation( 'getDerivatives', 'dvaldx_auto', istat, errmsg )

            call getDerivatives( x, val, dvaldx_auto )
!                ==============  i   i        o

            dvaldx_ptr => dvaldx_auto

         end if

      else

         nullify( dvaldx_ptr )

      end if
!
!  Table coefficient arrays will be 1 smaller than value array (except a0)
!
      nx = nx - 1
!
!  Deallocate arrays if already allocated
!
      if( allocated( table%x ) )  deallocate( table%x )
      if( allocated( table%a0 ) )  deallocate( table%a0 )
      if( allocated( table%a1 ) )  deallocate( table%a1 )
      if( allocated( table%a2 ) )  deallocate( table%a2 )
      if( allocated( table%a3 ) )  deallocate( table%a3 )
!
!  Copy x and val into the table
!
      allocate( table%x, source=x, stat=istat, errmsg=errmsg )
      call check_allocation( 'initTable1d', 'table%x', istat, errmsg )

      allocate( table%a0( nx+1 ), stat=istat, errmsg=errmsg )
      call check_allocation( 'initTable1d', 'table%a0', istat, errmsg )
      table%a0(:) = val(:)
!
!  Allocate table coefficients
!
      allocate( table%a1( nx ), stat=istat, errmsg=errmsg )
      call check_allocation( 'initTable1d', 'table%a1', istat, errmsg )

      if( associated(dvaldx_ptr) ) then
         allocate( table%a2( nx ), stat=istat, errmsg=errmsg )
         call check_allocation( 'initTable1d', 'table%a2', istat, errmsg )

         allocate( table%a3( nx ), stat=istat, errmsg=errmsg )
         call check_allocation( 'initTable1d', 'table%a3', istat, errmsg )
      end if
!
!  Some auxiliary arrays
!
      allocate( dx( nx ), stat=istat, errmsg=errmsg )
      call check_allocation( 'initTable1d', 'dx', istat, errmsg )
      dx = x(2:nx+1) - x(1:nx)

      allocate( f0( nx ), stat=istat, errmsg=errmsg )
      call check_allocation( 'initTable1d', 'f0', istat, errmsg )
      f0 = val(1:nx)

      allocate( f1( nx ), stat=istat, errmsg=errmsg )
      call check_allocation( 'initTable1d', 'f1', istat, errmsg )
      f1 = val(2:nx+1)

      if( associated(dvaldx_ptr) ) then
         allocate( fx0( nx ), stat=istat, errmsg=errmsg )
         call check_allocation( 'initTable1d', 'fx0', istat, errmsg )
         fx0 = dvaldx_ptr(1:nx) * dx(:)

         allocate( fx1( nx ), stat=istat, errmsg=errmsg )
         call check_allocation( 'initTable1d', 'fx1', istat, errmsg )
         fx1 = dvaldx_ptr(2:nx+1) * dx(:)
      end if
!
!  Set wx and determine if grid is regular
!
      table%wx = dx(1)
      table%regular = ( all( dx(:) - table%wx <= EPSMA*max( abs( dx(:) ), abs( table%wx ) ) ) )
!
!  Calculate table coefficients
!
      if( associated(dvaldx_ptr) ) then      ! coefficients for cubic interpolation
         table%a1 = fx0
         table%a2 = - 3*f0 + 3*f1 - 2*fx0 - fx1
         table%a3 = fx1 - 2*f1 + 2*f0 + fx0
      else                                   ! coefficients for linear interpolation
         table%a1 = f1 - f0
      end if

      if( cubic .and..not.present( dvaldx_ ) )  deallocate( dvaldx_auto )

   end subroutine initTable1d
!
!=======================================================================
!
   subroutine initTable2d( table, x, y, val, dvaldx_, dvaldy_ )
!             ===========    o    i  i   i      i        i
!
!     Initialize the table from values and optionally derivatives
!
!     Depending on whether derivatives are defined for the table, bilinear
!     or bicubic interpolation will be performed in interpolate2d
!
!  Arguments:
!
      type(Table2d_type), intent(out) :: table           ! table
      real(REAL8), intent(in) :: x(:)                    ! x values
      real(REAL8), intent(in) :: y(:)                    ! y values
      real(REAL8), intent(in) :: val(:,:)                ! values
      real(REAL8), intent(in), optional :: dvaldx_(:,:)  ! d(val)/dx
      real(REAL8), intent(in), optional :: dvaldy_(:,:)  ! d(val)/dy
!
!  Local Variables:
!
      character(256) :: errmsg
      integer(INT4) :: ix, iy, nx, ny, istat
      real(REAL8), allocatable :: dx(:,:), dy(:,:)
      real(REAL8), allocatable :: f00(:,:), f01(:,:), f10(:,:), f11(:,:)
      real(REAL8), allocatable :: fx00(:,:), fx01(:,:), fx10(:,:), fx11(:,:)
      real(REAL8), allocatable :: fy00(:,:), fy01(:,:), fy10(:,:), fy11(:,:)
!
!  Start of Executable Code:
!  ------------------------
!
!  Check the input variables
!
      nx = size( x )
      ny = size( y )

      if( size( val, 1 ) /= nx .or. size( val, 2 ) /= ny )  &
            error stop 'initTable2d: Shape of val is inconsistent with size of axes'

      if( present( dvaldx_ ) ) then
         if( .not.present( dvaldy_ ) )  &
               error stop 'initTable2d: Either none or both of dvaldx and dvaldy must be ' //  &
                          'present'
         if( .not.all( shape( dvaldx_ ) == shape( val ) ) )  &
               error stop 'initTable2d: Shape of dvaldx must match shape of val'
         if( .not.all( shape( dvaldy_ ) == shape( val ) ) )  &
               error stop 'initTable2d: Shape of dvaldy must match shape of val'
      else
         if( present( dvaldy_ ) )  &
               error stop 'initTable2d: Either none or both of dvaldx and dvaldy must be ' //  &
                          'present'
      end if
!
!  Table coefficient arrays will be 1 smaller than value array (except a00)
!
      nx = nx - 1
      ny = ny - 1
!
!  Deallocate arrays if already allocated
!
      if( allocated( table%x ) )  deallocate( table%x )
      if( allocated( table%y ) )  deallocate( table%y )
      if( allocated( table%a00 ) )  deallocate( table%a00 )
      if( allocated( table%a10 ) )  deallocate( table%a10 )
      if( allocated( table%a20 ) )  deallocate( table%a20 )
      if( allocated( table%a30 ) )  deallocate( table%a30 )
      if( allocated( table%a01 ) )  deallocate( table%a01 )
      if( allocated( table%a02 ) )  deallocate( table%a02 )
      if( allocated( table%a03 ) )  deallocate( table%a03 )
      if( allocated( table%a11 ) )  deallocate( table%a11 )
      if( allocated( table%a21 ) )  deallocate( table%a21 )
      if( allocated( table%a31 ) )  deallocate( table%a31 )
      if( allocated( table%a12 ) )  deallocate( table%a12 )
      if( allocated( table%a13 ) )  deallocate( table%a13 )
!
!  Copy x, y, and val into the table
!
      allocate( table%x, source=x, stat=istat, errmsg=errmsg )
      call check_allocation( 'initTable2d', 'table%x', istat, errmsg )

      allocate( table%y, source=y, stat=istat, errmsg=errmsg )
      call check_allocation( 'initTable2d', 'table%y', istat, errmsg )

      allocate( table%a00( nx+1,ny+1 ), stat=istat, errmsg=errmsg )
      call check_allocation( 'initTable2d', 'table%a00', istat, errmsg )
      table%a00(:,:) = val(:,:)
!
!  Allocate table coefficients
!
      allocate( table%a10( nx,ny ), stat=istat, errmsg=errmsg )
      call check_allocation( 'initTable2d', 'table%a10', istat, errmsg )

      allocate( table%a01( nx,ny ), stat=istat, errmsg=errmsg )
      call check_allocation( 'initTable2d', 'table%a01', istat, errmsg )

      allocate( table%a11( nx,ny ), stat=istat, errmsg=errmsg )
      call check_allocation( 'initTable2d', 'table%a11', istat, errmsg )

      if( present( dvaldx_ ) ) then
         allocate( table%a20( nx,ny ), stat=istat, errmsg=errmsg )
         call check_allocation( 'initTable2d', 'table%a20', istat, errmsg )

         allocate( table%a02( nx,ny ), stat=istat, errmsg=errmsg )
         call check_allocation( 'initTable2d', 'table%a02', istat, errmsg )

         allocate( table%a30( nx,ny ), stat=istat, errmsg=errmsg )
         call check_allocation( 'initTable2d', 'table%a30', istat, errmsg )

         allocate( table%a03( nx,ny ), stat=istat, errmsg=errmsg )
         call check_allocation( 'initTable2d', 'table%a03', istat, errmsg )

         allocate( table%a21( nx,ny ), stat=istat, errmsg=errmsg )
         call check_allocation( 'initTable2d', 'table%a21', istat, errmsg )

         allocate( table%a12( nx,ny ), stat=istat, errmsg=errmsg )
         call check_allocation( 'initTable2d', 'table%a12', istat, errmsg )

         allocate( table%a31( nx,ny ), stat=istat, errmsg=errmsg )
         call check_allocation( 'initTable2d', 'table%a31', istat, errmsg )

         allocate( table%a13( nx,ny ), stat=istat, errmsg=errmsg )
         call check_allocation( 'initTable2d', 'table%a13', istat, errmsg )
      end if
!
!  Some auxiliary arrays
!
      allocate( dx( nx,ny ), stat=istat, errmsg=errmsg )
      call check_allocation( 'initTable2d', 'dx', istat, errmsg )
      do iy = 1, ny
         dx(:,iy) = x(2:nx+1) - x(1:nx)
      end do

      allocate( dy( nx,ny ), stat=istat, errmsg=errmsg )
      call check_allocation( 'initTable2d', 'dy', istat, errmsg )
      do ix = 1, nx
         dy(ix,:) = y(2:ny+1) - y(1:ny)
      end do

      f00 = val(1:nx,1:ny)
      f01 = val(1:nx,2:ny+1)
      f10 = val(2:nx+1,1:ny)
      f11 = val(2:nx+1,2:ny+1)

      if( present( dvaldx_ ) ) then
         fx00 = dvaldx_(1:nx,1:ny) * dx(:,:)
         fx01 = dvaldx_(1:nx,2:ny+1) * dx(:,:)
         fx10 = dvaldx_(2:nx+1,1:ny) * dx(:,:)
         fx11 = dvaldx_(2:nx+1,2:ny+1) * dx(:,:)
         fy00 = dvaldy_(1:nx,1:ny) * dy(:,:)
         fy01 = dvaldy_(1:nx,2:ny+1) * dy(:,:)
         fy10 = dvaldy_(2:nx+1,1:ny) * dy(:,:)
         fy11 = dvaldy_(2:nx+1,2:ny+1) * dy(:,:)
      end if
!
!  Set wx, wy and determine if grid is regular
!
      table%wx = dx(1,1)
      table%wy = dy(1,1)
      table%regular = ( all( abs( dx(:,1) - table%wx )  &
                             <= EPSMA*max( abs( dx(:,1) ), abs( table%wx ) ) ) .and.  &
                        all( abs( dy(1,:) - table%wy )  &
                             <= EPSMA*max( abs( dy(1,:) ), abs( table%wy ) ) ) )
!
!  Calculate table coefficients
!
      if( present( dvaldx_ ) ) then              ! coefficients for bicubic interpolation
         table%a10 = fx00
         table%a01 = fy00

         table%a11 = (f10 + f01 - f11 - f00) + fx01 + fy10 - fx00 - fy00

         table%a20 = - 3*f00 + 3*f10 - 2*fx00 - fx10
         table%a02 = - 3*f00 + 3*f01 - 2*fy00 - fy01

         table%a30 = fx10 - 2*f10 + 2*f00 + fx00
         table%a03 = fy01 - 2*f01 + 2*f00 + fy00

         table%a21 = - 3*(f10 + f01 - f11 - f00) - 2*fx01 + 2*fx00 - fx11 + fx10
         table%a12 = - 3*(f10 + f01 - f11 - f00) - 2*fy10 + 2*fy00 - fy11 + fy01

         table%a31 = 2*(f10 + f01 - f11 - f00) + fx11 - fx10 + fx01 - fx00
         table%a13 = 2*(f10 + f01 - f11 - f00) + fy11 - fy01 + fy10 - fy00
      else                                      ! coefficients for bilinear interpolation
         table%a10 = f10 - f00
         table%a01 = f01 - f00
         table%a11 = f11 + f00 - f10 - f01
      end if

   end subroutine initTable2d
!
!=======================================================================
!
   subroutine initTable3d( table, x, y, z, val )
!             ===========    o    i  i  i   i
!
!     Initialize the table from values
!
!     Note: cubic interpolation has not been implemented
!
!  Arguments:
!
      type(Table3d_type), intent(out) :: table           ! table
      real(REAL8), intent(in) :: x(:)                    ! x values
      real(REAL8), intent(in) :: y(:)                    ! y values
      real(REAL8), intent(in) :: z(:)                    ! y values
      real(REAL8), intent(in) :: val(:,:,:)              ! values
!
!  Local Variables:
!
      character(256) :: errmsg
      integer(INT4) :: nx, ny, nz, istat
      real(REAL8), allocatable :: dx(:), dy(:), dz(:)
      real(REAL8), allocatable :: f000(:,:,:), f100(:,:,:), f010(:,:,:), f110(:,:,:)
      real(REAL8), allocatable :: f001(:,:,:), f101(:,:,:), f011(:,:,:), f111(:,:,:)
!
!  Start of Executable Code:
!  ------------------------
!
!  Check the input variables
!
      nx = size( x )
      ny = size( y )
      nz = size( z )

      if( size( val, 1 ) /= nx .or. size( val, 2 ) /= ny .or. size( val, 3 ) /= nz )  &
            error stop 'initTable3d: Shape of val is inconsistent with size of axes'
!
!  Table coefficient arrays will be 1 smaller than value array (except a000)
!
      nx = nx - 1
      ny = ny - 1
      nz = nz - 1
!
!  Deallocate arrays if already allocated
!
      if( allocated( table%x ) )  deallocate( table%x )
      if( allocated( table%y ) )  deallocate( table%y )
      if( allocated( table%z ) )  deallocate( table%z )
      if( allocated( table%a000 ) )  deallocate( table%a000 )
      if( allocated( table%a100 ) )  deallocate( table%a100 )
      if( allocated( table%a010 ) )  deallocate( table%a010 )
      if( allocated( table%a110 ) )  deallocate( table%a110 )
      if( allocated( table%a001 ) )  deallocate( table%a001 )
      if( allocated( table%a101 ) )  deallocate( table%a101 )
      if( allocated( table%a011 ) )  deallocate( table%a011 )
      if( allocated( table%a111 ) )  deallocate( table%a111 )
!
!  Copy x, y, z, and val into the table
!
      allocate( table%x, source=x, stat=istat, errmsg=errmsg )
      call check_allocation( 'initTable3d', 'table%x', istat, errmsg )

      allocate( table%y, source=y, stat=istat, errmsg=errmsg )
      call check_allocation( 'initTable3d', 'table%y', istat, errmsg )

      allocate( table%z, source=z, stat=istat, errmsg=errmsg )
      call check_allocation( 'initTable3d', 'table%z', istat, errmsg )

      allocate( table%a000, source=val, stat=istat, errmsg=errmsg )
      call check_allocation( 'initTable3d', 'table%a000', istat, errmsg )
!
!  Allocate table coefficients
!
      allocate( table%a100( nx,ny,nz ), stat=istat, errmsg=errmsg )
      call check_allocation( 'initTable3d', 'table%a100', istat, errmsg )

      allocate( table%a010( nx,ny,nz ), stat=istat, errmsg=errmsg )
      call check_allocation( 'initTable3d', 'table%a010', istat, errmsg )

      allocate( table%a110( nx,ny,nz ), stat=istat, errmsg=errmsg )
      call check_allocation( 'initTable3d', 'table%a110', istat, errmsg )

      allocate( table%a001( nx,ny,nz ), stat=istat, errmsg=errmsg )
      call check_allocation( 'initTable3d', 'table%a001', istat, errmsg )

      allocate( table%a101( nx,ny,nz ), stat=istat, errmsg=errmsg )
      call check_allocation( 'initTable3d', 'table%a101', istat, errmsg )

      allocate( table%a011( nx,ny,nz ), stat=istat, errmsg=errmsg )
      call check_allocation( 'initTable3d', 'table%a011', istat, errmsg )

      allocate( table%a111( nx,ny,nz ), stat=istat, errmsg=errmsg )
      call check_allocation( 'initTable3d', 'table%a111', istat, errmsg )
!
!  Some auxiliary arrays
!
      allocate( dx( nx ), stat=istat, errmsg=errmsg )
      call check_allocation( 'initTable3d', 'dx', istat, errmsg )
      dx(:) = x(2:nx+1) - x(1:nx)

      allocate( dy( ny ), stat=istat, errmsg=errmsg )
      call check_allocation( 'initTable3d', 'dy', istat, errmsg )
      dy(:) = y(2:ny+1) - y(1:ny)

      allocate( dz( nz ), stat=istat, errmsg=errmsg )
      call check_allocation( 'initTable3d', 'dz', istat, errmsg )
      dz(:) = z(2:nz+1) - z(1:nz)

      f000 = val(1:nx,1:ny,1:nz)
      f010 = val(1:nx,2:ny+1,1:nz)
      f100 = val(2:nx+1,1:ny,1:nz)
      f110 = val(2:nx+1,2:ny+1,1:nz)
      f001 = val(1:nx,1:ny,2:nz+1)
      f011 = val(1:nx,2:ny+1,2:nz+1)
      f101 = val(2:nx+1,1:ny,2:nz+1)
      f111 = val(2:nx+1,2:ny+1,2:nz+1)
!
!  Set wx, wy, wz, and determine if grid is regular
!
      table%wx = dx(1)
      table%wy = dy(1)
      table%wz = dz(1)
      table%regular = ( all( abs( dx(:) - table%wx )  &
                             <= EPSMA*max( abs( dx(:) ), abs( table%wx ) ) ) .and.  &
                        all( abs( dy(:) - table%wy )  &
                             <= EPSMA*max( abs( dy(:) ), abs( table%wy ) ) ) .and.  &
                        all( abs( dz(:) - table%wz )  &
                             <= EPSMA*max( abs( dz(:) ), abs( table%wz ) ) ) )
!
!  Calculate table coefficients
!
      table%a100 = f100 - f000
      table%a010 = f010 - f000
      table%a001 = f001 - f000
      table%a110 = f000 - f100 - f010 + f110
      table%a101 = f000 - f100 - f001 + f101
      table%a011 = f000 - f010 - f001 + f011
      table%a111 = -f000 + f100 + f010 + f001 - f110 - f101 - f011 + f111

   end subroutine initTable3d
!
!=======================================================================
!
   subroutine interpolate1d( table, x, val, dvaldx_ )
!             =============    i    i   o      o
!
!     Interpolate in the 1D table
!
!     Depending on whether derivatives are defined for the table, linear
!     or cubic interpolation is performed
!
!  Arguments:
!
      type(Table1d_type), intent(in) :: table      ! table
      real(REAL8), intent(in) :: x                 ! x value
      real(REAL8), intent(out) :: val              ! interpolated value
      real(REAL8), intent(out), optional :: dvaldx_! derivative of val wrt x
!
!  Local Variables:
!
      integer(INT4) :: ix
      real(REAL8) :: wx, xbar
!
!  Start of Executable Code:
!  ------------------------
!
      if( table%regular ) then
         ix = max( 1, min( size( table%x )-1, ceiling( (x - table%x(1))/table%wx ) ) )
         wx = table%wx
      else
         call binSearch( size( table%x ), table%x, x, ix )
!             =========         i            i     i  o
         wx = table%x(ix+1) - table%x(ix)
      end if

      xbar = (x - table%x(ix)) / wx

      if( allocated( table%a2 ) ) then
         val = ( (table%a3(ix)*xbar + table%a2(ix)) * xbar + table%a1(ix) ) * xbar  &
               + table%a0(ix)
      else
         val = table%a1(ix) * xbar + table%a0(ix)
      end if

      if( present( dvaldx_ ) ) then
         if( allocated( table%a2 ) ) then
            dvaldx_ = (3*table%a3(ix)*xbar + 2*table%a2(ix)) * xbar + table%a1(ix)
            dvaldx_ = dvaldx_ / wx
         else
            dvaldx_ = table%a1(ix) / wx
         end if
      end if

   end subroutine interpolate1d
!
!=======================================================================
!
   subroutine interpolate2d( table, x, y, val, dvaldx_, dvaldy_ )
!             =============    i    i  i   o      o        o
!
!     Interpolate in the 2D table
!
!     Depending on whether derivatives are defined for the table, linear
!     or cubic interpolation is performed
!
!  Arguments:
!
      type(Table2d_type), intent(in) :: table      ! table
      real(REAL8), intent(in) :: x                 ! x value
      real(REAL8), intent(in) :: y                 ! y value
      real(REAL8), intent(out) :: val              ! interpolated value
      real(REAL8), intent(out), optional :: dvaldx_! derivative of val wrt x
      real(REAL8), intent(out), optional :: dvaldy_! derivative of val wrt y
!
!  Local Variables:
!
      integer(INT4) :: ix, iy
      real(REAL8) :: wx, wy, xbar, ybar
!
!  Start of Executable Code:
!  ------------------------
!
      if( table%regular ) then
         ix = max( 1, min( size( table%x )-1, ceiling( (x - table%x(1))/table%wx ) ) )
         iy = max( 1, min( size( table%y )-1, ceiling( (y - table%y(1))/table%wy ) ) )
         wx = table%wx
         wy = table%wy
      else
         call binSearch( size(table%x), table%x, x, ix )
!             =========        i           i     i  o
         call binSearch( size(table%y), table%y, y, iy )
!             =========        i           i     i  o
         wx = table%x(ix+1) - table%x(ix)
         wy = table%y(iy+1) - table%y(iy)
      end if

      xbar = (x - table%x(ix)) / wx
      ybar = (y - table%y(iy)) / wy

      if( allocated( table%a20 ) ) then
         val = ( ( (table%a13(ix,iy)*xbar+table%a03(ix,iy)) * ybar  &
                   + (table%a12(ix,iy)*xbar+table%a02(ix,iy)) ) * ybar  &
                 + ((table%a31(ix,iy)*xbar+table%a21(ix,iy))*xbar+table%a11(ix,iy))*xbar  &
                   +table%a01(ix,iy) ) * ybar  &
               + ((table%a30(ix,iy)*xbar+table%a20(ix,iy))*xbar+table%a10(ix,iy))*xbar  &
                 +table%a00(ix,iy)
         if( present( dvaldx_ ) ) then  ! untested
            dvaldx_ = ( ( table%a13(ix,iy) * ybar + table%a12(ix,iy) ) * ybar  &
                        + ((3*table%a31(ix,iy)*xbar+2*table%a21(ix,iy))*xbar+table%a11(ix,iy)) &
                          ) * ybar  &
                      + ((3*table%a30(ix,iy)*xbar+2*table%a20(ix,iy))*xbar+table%a10(ix,iy))
            dvaldx_ = dvaldx_ / wx
         end if
         if( present( dvaldy_ ) ) then  ! untested
            dvaldy_ = ( 3*(table%a13(ix,iy)*xbar+table%a03(ix,iy)) * ybar  &
                       + 2*(table%a12(ix,iy)*xbar+table%a02(ix,iy)) ) * ybar  &
                     + ((table%a31(ix,iy)*xbar+table%a21(ix,iy))*xbar+table%a11(ix,iy))*xbar  &
                       +table%a01(ix,iy)
            dvaldy_ = dvaldy_ / wy
         end if
      else
         val = ( table%a11(ix,iy)*xbar + table%a01(ix,iy) ) * ybar  &
               + table%a10(ix,iy)*xbar + table%a00(ix,iy)
         if( present( dvaldx_ ) )  &
               dvaldx_ = ( table%a11(ix,iy) * ybar + table%a10(ix,iy) ) / wx
         if( present( dvaldy_ ) )  &
               dvaldy_ = ( table%a11(ix,iy) * xbar + table%a01(ix,iy) ) / wy
      end if

   end subroutine interpolate2d
!
!=======================================================================
!
   subroutine interpolate3d( table, x, y, z, val, dvaldx_, dvaldy_, dvaldz_ )
!             =============    i    i  i  i   o      o        o        o
!
!     Interpolate in the 3D table
!
!  Arguments:
!
      type(Table3d_type), intent(in) :: table      ! table
      real(REAL8), intent(in) :: x                 ! x value
      real(REAL8), intent(in) :: y                 ! y value
      real(REAL8), intent(in) :: z                 ! y value
      real(REAL8), intent(out) :: val              ! interpolated value
      real(REAL8), intent(out), optional :: dvaldx_! derivative of val wrt x
      real(REAL8), intent(out), optional :: dvaldy_! derivative of val wrt y
      real(REAL8), intent(out), optional :: dvaldz_! derivative of val wrt z
!
!  Local Variables:
!
      integer(INT4) :: ix, iy, iz
      real(REAL8) :: wx, wy, wz, xbar, ybar, zbar
!
!  Start of Executable Code:
!  ------------------------
!
      if( table%regular ) then
         ix = max( 1, min( size( table%x )-1, ceiling( (x - table%x(1))/table%wx ) ) )
         iy = max( 1, min( size( table%y )-1, ceiling( (y - table%y(1))/table%wy ) ) )
         iz = max( 1, min( size( table%z )-1, ceiling( (z - table%z(1))/table%wz ) ) )
         wx = table%wx
         wy = table%wy
         wz = table%wz
      else
         call binSearch( size( table%x ), table%x, x, ix )
!             =========         i            i     i  o
         call binSearch( size( table%y ), table%y, y, iy )
!             =========         i            i     i  o
         call binSearch( size( table%z ), table%z, z, iz )
!             =========         i            i     i  o
         wx = table%x(ix+1) - table%x(ix)
         wy = table%y(iy+1) - table%y(iy)
         wz = table%z(iz+1) - table%z(iz)
      end if

      xbar = (x - table%x(ix)) / wx
      ybar = (y - table%y(iy)) / wy
      zbar = (z - table%z(iz)) / wz

      val = ( ( table%a111(ix,iy,iz)*xbar + table%a011(ix,iy,iz) ) * ybar  &
              + ( table%a101(ix,iy,iz)*xbar + table%a001(ix,iy,iz) ) ) * zbar  &
            + ( table%a110(ix,iy,iz)*xbar + table%a010(ix,iy,iz) ) * ybar  &
            + ( table%a100(ix,iy,iz)*xbar + table%a000(ix,iy,iz) )
      if( present( dvaldx_ ) )  &
            dvaldx_ = ( ( table%a111(ix,iy,iz) * ybar + table%a101(ix,iy,iz) ) * zbar  &
                        + ( table%a110(ix,iy,iz) * ybar + table%a100(ix,iy,iz) ) ) / wx
      if( present( dvaldy_ ) )  &
            dvaldy_ = ( ( table%a111(ix,iy,iz)*xbar + table%a011(ix,iy,iz) ) * zbar  &
                        + ( table%a110(ix,iy,iz)*xbar + table%a010(ix,iy,iz) ) ) / wy
      if( present( dvaldz_ ) )  &
            dvaldz_ = ( ( table%a111(ix,iy,iz)*xbar + table%a011(ix,iy,iz) ) * ybar  &
                        + ( table%a101(ix,iy,iz)*xbar + table%a001(ix,iy,iz) ) ) / wz

   end subroutine interpolate3d
!
!=======================================================================
!
   subroutine getDerivatives( x, y, dydx )
!             ==============  i  i   o
!
!     Determine the derivatives dydx of tabulated data of the function y(x) at abscissae x.
!
!     Derivatives are determined as the derivatives of a cubic spline through the data with
!     not-a-knot boundary conditions. The algorithm is taken from the SciPy library class
!     CubicSpline.
!
!  Arguments:
!
      real(REAL8), intent(in) :: x(:)       ! independent variable
      real(REAL8), intent(in) :: y(:)       ! dependent variable
      real(REAL8), intent(out) :: dydx(:)   ! derivatives
!
!  Local Variables:
!
      character(256) :: errmsg
      integer(INT4) :: nx, istat
      real(REAL8), allocatable :: dx(:), dy(:), am(:), au(:), al(:), b(:)
!
!  Start of Executable Code:
!  ------------------------
!
!  Allocate arrays
!
      nx = size( x )

      allocate( dx( nx-1 ), dy( nx-1 ), stat=istat, errmsg=errmsg )
      call check_allocation( 'getDerivatives', 'dx, dy', istat, errmsg )

      allocate( am( nx ), au( nx ), al( nx ), b( nx ), stat=istat, errmsg=errmsg )
      call check_allocation( 'getDerivatives', 'am, au, al, b', istat, errmsg )
!
!  Set up tridiagonal system of equations
!
      dx(1:nx-1) = x(2:nx) - x(1:nx-1)
      dy(1:nx-1) = y(2:nx) - y(1:nx-1)

      am(1) = dx(2)
      am(2:nx-1) = 2 * (dx(1:nx-2) + dx(2:nx-1))
      am(nx) = dx(nx-2)

      au(1) = dx(1) + dx(2)
      au(2:nx-1) = dx(1:nx-2)

      al(2:nx-1) = dx(2:nx-1)
      al(nx) = dx(nx-2) + dx(nx-1)

      b(1) = ( (3*dx(1) + 2*dx(2)) * dx(2) * dy(1)/dx(1) + dx(1)**2 * dy(2)/dx(2) )  &
             / (dx(1) + dx(2))
      b(2:nx-1) = 3 * (dx(2:nx-1) * dy(1:nx-2)/dx(1:nx-2) + dx(1:nx-2) * dy(2:nx-1)/dx(2:nx-1))
      b(nx) = ( dx(nx-1)**2 * dy(nx-2)/dx(nx-2)  &
                + (3*dx(nx-1) + 2*dx(nx-2)) * dx(nx-2) * dy(nx-1)/dx(nx-1) )  &
              / (dx(nx-2) + dx(nx-1))
!
!  Solve the system of equations
!
      call tridiag( al, am, au, b, dydx )
!          =======  i   i   i   i   o

   end subroutine getDerivatives
!
!=======================================================================
!
end module Math_Table_mod
