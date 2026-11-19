!
!##############################################################################################
!
module Test_Quadrature_mod
!
!  Test the Math_Quadrature_mod module
!
!  Used Modules:
!  ============
!
   use Const_mod,          only: INT4, REAL8, PI
   use Math_Quadrature_mod
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
   public :: testLegendre                       ! test Gauss-Legendre quadrature
!            ============
   public :: testMehler                         ! test Gauss-Mehler quadrature
!            ==========
!
!  Private Procedures
!  ==================
!
   private :: square                            ! test function
!             ======
!
!==============================================================================================
!
contains
!
!==============================================================================================
!
   real(REAL8) function square( x )
!
!     Take the square of the argument
!
!  Arguments:
!
      real(REAL8), intent(in) :: x              ! number to be squared
!
!  Start of Executable Code
!  ------------------------
!
      square = x**2

   end function square
!
!==============================================================================================
!
   subroutine testLegendre()
!
!     Test the convergence of the Gauss-Legendre quadrature
!
!  Arguments:
!
!     None
!
!  Local Variables:
!
      integer(INT4), parameter :: nmax = 8
      integer(INT4) :: n
      real(REAL8) :: result
!
!  Start of Executable Code:
!  ------------------------
!
      print *, 'Int_0^1 x^2 dx = 1/3'
      do n = 1, nmax
         result = quadLegendre( square, n )
!                 ============    p     i
         print '(a,i2,a,g15.7)', 'n=', n, ', integral=', result
      end do

   end subroutine testLegendre
!
!==============================================================================================
!
   subroutine testMehler()
!
!     Test the convergence of the Gauss-Mehler quadratures
!
!  Arguments:
!
!     None
!
!  Local Variables:
!
      integer(INT4), parameter :: nmax = 8
      integer(INT4) :: n
      real(REAL8) :: resultSymm, resultLin
!
!  Start of Executable Code:
!  ------------------------
!
      print *, 'Int_0^1 x^2/sqrt(1-x^2) dx = PI/4 =', PI/4
      do n = 1, nmax
         resultSymm = quadMehlerSymm( square, n )
!                     ==============    p     i
         resultLin = quadMehlerLin( square, n )
!                    =============    p     i
         print '(a,i2,2(a,g15.7))', 'n=', n, ', integral(symm)=', resultSymm,  &
                                             ', integral(lin)=', resultLin
      end do

   end subroutine testMehler
!
!==============================================================================================
!
end module Test_Quadrature_mod
!
!##############################################################################################
!
program test

   use const_mod, only: INT4
   use Test_Quadrature_mod

   integer(INT4) :: choice

   do
      print *, ' '
      print *, 'Choose among the following tests:'
      print *, ' '
      print *, '   1 ... test Gauss-Legendre quadrature'
      print *, '   2 ... test Gauss-Mehler quadrature'
      print *, ' '
      read *, choice

      select case( choice )
         case( 1 )
            call testLegendre( )
         case( 2 )
            call testMehler( )
         case default
            print *, 'Invalid answer.'
      end select
   end do

end program test
