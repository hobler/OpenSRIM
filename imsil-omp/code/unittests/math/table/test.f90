module Test_Table_mod
!
!     Tests the Table_mod module
!
!  Used Modules:
!  ============
!
   use Const_mod,          only: INT4, REAL8,  &
                                 ZERO
   use open_mod,           only: openFile
   use Math_Table_mod
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
   public :: testPolynomial                     ! test derivatives of polynomial
!            ==============
   public :: testTable2d                        ! test interpolation in a table2d
!            ===========
   public :: testDerivatives                    ! test derivative calculation
!            ===============
   public :: testTable1d                        ! test interpolation in a table1d
!            ===========
!
!=======================================================================
!
contains
!
!=======================================================================
!
   subroutine polynomial( coeff, x, y, val, dvaldx, dvaldy )
!             ==========    i    i  i   o     o       o
!
!     Evaluate a bicubic polynomial and its derivatives given its coefficients at
!     a given point (x,y).
!
!  Arguments:
!
      real(REAL8), intent(in) :: coeff(12)      ! coefficients of the polynomial
      real(REAL8), intent(in) :: x              ! x value for evaluation of polynomial
      real(REAL8), intent(in) :: y              ! y value for evaluation of polynomial
      real(REAL8), intent(out) :: val           ! value of polynomial at (x,y)
      real(REAL8), intent(out) :: dvaldx        ! d(val)/dx at (x,y)
      real(REAL8), intent(out) :: dvaldy        ! d(val)/dy at (x,y)
!
!  Local Variables:
!
!     None
!
!  Start of Executable Code:
!  ------------------------
!
      val = coeff(1) + coeff(2)*x + coeff(3)*x**2 + coeff(4)*x**3  &
            + coeff(5)*y + coeff(6)*x*y + coeff(7)*x**2*y + coeff(8)*x**3*y  &
            + coeff(9)*y**2 + coeff(10)*x*y**2  &
            + coeff(11)*y**3 + coeff(12)*x*y**3

      dvaldx = coeff(2) + 2*coeff(3)*x + 3*coeff(4)*x**2  &
               + coeff(6)*y + 2*coeff(7)*x*y + 3*coeff(8)*x**2*y  &
               + coeff(10)*y**2  &
               + coeff(12)*y**3

      dvaldy = coeff(5) + coeff(6)*x + coeff(7)*x**2 + coeff(8)*x**3  &
               + 2*coeff(9)*y + 2*coeff(10)*x*y  &
               + 3*coeff(11)*y**2 + 3*coeff(12)*x*y**2

   end subroutine polynomial
!
!=======================================================================
!
   subroutine testPolynomial()
!
!     Test the derivatives of the polynomial
!
!  Arguments:
!
!     None
!
!  Local Variables:
!
      real(REAL8) :: x, y, val00, val10, val01, dvaldx, dvaldy, dummy
      real(REAL8) :: coeff(12)
!
!  Start of Executable Code:
!  ------------------------
!
      do
         call random_number( coeff )
!             =============    o
         call random_number( x )
!             =============  o
         call random_number( y )
!             =============  o
         print *, 'x=', x, 'y=', y
         call polynomial( coeff, x, y, val00, dvaldx, dvaldy )
!             ==========    i    i  i    o      o        o
         call polynomial( coeff, x+0.0001, y, val10, dummy, dummy )
!             ==========    i      i     i    o      o      o
         call polynomial( coeff, x, y+0.0001, val01, dummy, dummy )
!             ==========    i    i    i       o      o      o
         print *, 'val=', val00
         print *, 'dvaldx(der)=', dvaldx, 'dvaldx(diff)=', (val10-val00)/0.0001
         print *, 'dvaldy(der)=', dvaldy, 'dvaldy(diff)=', (val01-val00)/0.0001
         print *, ' '
         print *, 'Enter <CR> to continue'
         read *
      end do

   end subroutine testPolynomial
!
!=======================================================================
!
   subroutine testTable2d()
!
!     Test the interpolation in a table
!
!  Arguments:
!
!     None
!
!  Local Variables:
!
      integer(INT4), parameter :: nx=3, ny=4
      integer(INT4) :: ix, iy
      real(REAL8) :: xval, yval, dummy, valexact, valint
      real(REAL8) :: x(nx), y(ny), val(nx,ny), dvaldx(nx,ny), dvaldy(nx,ny)
      real(REAL8) :: coeff(12)
      type(Table2d_type) :: table
!
!  Start of Executable Code:
!  ------------------------
!
      x(:) = [ 0., 0.3, 1. ]
      y(:) = [ 0., 0.6, 0.8, 1. ]
      do
         call random_number( coeff )
!             =============    o
         do iy = 1, ny
            do ix = 1, nx
               call polynomial( coeff, x(ix), y(iy), val(ix,iy), dvaldx(ix,iy), dvaldy(ix,iy) )
!                   ==========    i      i      i        o             o              o
            end do
         end do

         call initTable( table, x, y, val, dvaldx, dvaldy )
!             =========    o    i  i   i     i       i

         call random_number( xval )
!             =============   o
         call random_number( yval )
!             =============   o

         !xval = 0.5
         !yval = 1

         call interpolate( table, xval, yval, valint )
!             ===========    i     i     i      o
         call polynomial( coeff, xval, yval, valexact, dummy, dummy )
!             ==========    i     i     i       o        o      o

         print *, 'x=', xval, 'y=', yval
         print *, 'value=', valexact, 'interpolated=', valint
         print *, ' '
         print *, 'Enter <CR> to continue'
         read *
      end do

   end subroutine testTable2d
!
!=======================================================================
!
   subroutine testDerivatives()
!
!     Test the subroutine for taking the derivatives
!
!  Arguments:
!
!     None
!
!  Local Variables:
!
      integer(INT4), parameter :: nx=100
      integer(INT4) :: ix
      real(REAL8) :: x(0:nx), y(0:nx), dydx(0:nx)
!
!  Start of Executable Code:
!  ------------------------
!
      !x(:) = [ ( ix/20., ix=0,nx ) ]
      x(:) = [ ( (ix/20.)**2 / 5, ix=0,nx ) ]
      y(:) = exp(-x(:)) - exp(-10*x(:))

      call getDerivatives( x, y, dydx )
!          ==============  i  i   o

      call openFile( 21, 'deriv.dat', 'formatted', 'write' )
!          ========  i        i            i          i
      do ix = 0, nx
         write( 21, * )  x(ix), dydx(ix)
      end do

      close( 21 )
!
!  Write and run python script to display errors
!
      call openFile( 12, 'plot_deriv.py', 'formatted', 'write' )
!          ========  i          i               i          i
      write( 12, '(a)' )  'import numpy as np'
      write( 12, '(a)' )  'from scipy.interpolate import CubicSpline'
      write( 12, '(a)' )  'import matplotlib.pyplot as plt'
      write( 12, '(a)' )  ' '
      write( 12, '(a)' )  'x = np.linspace(0, 5, 101)'
      write( 12, '(a)' )  'y = np.exp(-x) - np.exp(-10*x)'
      write( 12, '(a)' )  'ppoly = CubicSpline(x, y)'
      write( 12, '(a)' )  'dydx = ppoly.derivative()(x)'
      write( 12, '(a)' )  ' '
      write( 12, '(a)' )  'plt.plot(x, dydx, "x-", label="SciPy")'
      write( 12, '(a)' )  ' '
      write( 12, '(a)' )  'x, dydx = np.loadtxt("deriv.dat", unpack=True)'
      write( 12, '(a)' )  'plt.plot(x, dydx, ".--", label="IMSIL")'
      write( 12, '(a)' )  ' '
      write( 12, '(a)' )  'plt.legend(loc="best")'
      write( 12, '(a)' )  'plt.xlabel("x")'
      write( 12, '(a)' )  'plt.ylabel("dy/dx")'
      write( 12, '(a)' )  'plt.tight_layout()'
      write( 12, '(a)' )  'plt.show()'
      close( 12 )
      call execute_command_line('python plot_deriv.py')

   end subroutine testDerivatives
!
!=======================================================================
!
   subroutine testTable1d()
!
!     Test the interpolation in a 1d table
!
!  Arguments:
!
!     None
!
!  Local Variables:
!
      integer(INT4), parameter :: nx=1000
      integer(INT4) :: ix
      real(REAL8) :: x(0:nx), y(0:nx)
      type(Table1d_type) :: table
!
!  Start of Executable Code:
!  ------------------------
!
      x(:) = [ ( (ix/200.)**2 / 5, ix=0,nx ) ]
      y(:) = exp(-x(:)) - exp(-10*x(:))

      call initTable( table, x(::10), y(::10) )!, cubic=.true. )
!          =========    o    i  i       i

      call openFile( 21, 'table1d.dat', 'formatted', 'write' )
!          ========  i        i            i          i
      do ix = 0, nx
         call interpolate( table, x(ix), y(ix) )
!             ===========    i      i      o
         write( 21, * )  x(ix), y(ix)
      end do

      close( 21 )
!
!  Write and run python script to display errors
!
      call openFile( 12, 'plot_table1d.py', 'formatted', 'write' )
!          ========  i          i               i          i
      write( 12, '(a)' )  'import numpy as np'
      write( 12, '(a)' )  'from scipy.interpolate import CubicSpline'
      write( 12, '(a)' )  'import matplotlib.pyplot as plt'
      write( 12, '(a)' )  ' '
      write( 12, '(a)' )  'x = np.linspace(0, 5, 101)'
      write( 12, '(a)' )  'y = np.exp(-x) - np.exp(-10*x)'
      write( 12, '(a)' )  'plt.plot(x, y, "x", label="SciPy")'
      write( 12, '(a)' )  'ppoly = CubicSpline(x, y)'
      write( 12, '(a)' )  'x = np.linspace(0, 5, 1001)'
      write( 12, '(a)' )  'y = ppoly(x)'
      write( 12, '(a)' )  ' '
      write( 12, '(a)' )  'plt.plot(x, y, "-")'
      write( 12, '(a)' )  ' '
      write( 12, '(a)' )  'x, y = np.loadtxt("table1d.dat", unpack=True)'
      write( 12, '(a)' )  'plt.plot(x, y, ".--", label="IMSIL")'
      write( 12, '(a)' )  ' '
      write( 12, '(a)' )  'plt.legend(loc="best")'
      write( 12, '(a)' )  'plt.xlabel("x")'
      write( 12, '(a)' )  'plt.ylabel("y")'
      write( 12, '(a)' )  'plt.tight_layout()'
      write( 12, '(a)' )  'plt.show()'
      close( 12 )
      call execute_command_line('python plot_table1d.py')

   end subroutine testTable1d
!
!=======================================================================
!
end module Test_Table_mod
!
!#######################################################################
!
program test

   use const_mod, only: INT4
   use Test_Table_mod

   integer(INT4) :: choice

   do
      print *, ' '
      print *, 'Choose among the following tests:'
      print *, ' '
      print *, '   1 ... test derivatives of 2d test polynomial'
      print *, '   2 ... test interpolation using the 2d test polynomial'
      print *, '   3 ... test automatic derivatives in 1d interpolation'
      print *, '   4 ... test 1d interpolation'
      print *, ' '
      read *, choice

      select case( choice )
         case( 1 )
            call testPolynomial()
         case( 2 )
            call testTable2d()
         case( 3 )
            call testDerivatives()
         case( 4 )
            call testTable1d()
         case default
            print *, 'Invalid answer.'
      end select
   end do

end program test
