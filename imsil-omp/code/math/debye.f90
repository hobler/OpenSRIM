module Math_Debye_mod
!
!  Calculate Debye function of order one
!
!  Used modules:
!  ============
!
   use Const_mod,  only: REAL8, INT4,  &
                         ZERO, HALF, EIGHT, FOUR, ONE, QUART, TWO 
!
!  Default type of access:
!  ======================
!
   implicit none
   private 
   save
!
!  Public procedures:
!  =================
!
   public :: debye1                        ! calculate Debye function of order 1
!            ======
!
!  Private procedures
!  ==================
!
   private :: cheval                       ! evaluates a Chebyshev series
!             ======

contains
!
!=========================================================================================
!
   function debye1( xvalue ) 
!           ======    i
!
!     DEBYE1 calculates the Debye function of order 1.
!
!    The function is defined by:
!
!      DEBYE1(x) = 1 / x * Integral ( 0 <= t <= x ) t / ( exp ( t ) - 1 ) dt
!
!    The code uses Chebyshev series whose coefficients are given to 20 decimal places.
!
!  Author:
!
!    Allan McLeod,
!    Department of Mathematics and Statistics,
!    Paisley University, High Street, Paisley, Scotland, PA12BE
!    macl_ms0@paisley.ac.uk
!
!  Reference:
!
!    Allan Mcleod,
!    Algorithm 757, MISCFUN: A software package to compute uncommon 
!      special functions,
!    ACM Transactions on Mathematical Software,
!    Volume 22, Number 3, September 1996, pages 288-301.
!
!  Arguments:
!
      real(REAL8), intent(in) :: xvalue   ! input value
      real(REAL8) :: debye1                 ! returning value
!
!  Local variables
!
      real(REAL8) :: expmx, rk, sum1, t, xk, x
      real(REAL8), parameter :: debinf = 0.60792710185402662866d0
      real(REAL8), parameter :: nine = 9.0d0, thirt6=36.0d0, xlim = 708.39642d0
      real(REAL8), parameter :: xlow = 0.298023d-7, xupper = 35.35051d0
      real(REAL8), parameter :: adeb1(0:18) = [ 2.40065971903814101941d0, &
           0.19372130421893600885d0, -0.623291245548957703d-2, 0.35111747702064800d-3, & 
             -0.2282224667012310d-4,     0.158054678750300d-5,   -0.11353781970719d-6, &
                  0.835833611875d-8,        -0.62644247872d-9,       0.4760334890d-10, &
                   -0.365741540d-11,           0.28354310d-12,         -0.2214729d-13, &
                       0.174092d-14,             -0.13759d-15,             0.1093d-16, & 
                          -0.87d-18,                  0.7d-19,               -0.1d-19 ]

      integer, parameter :: nterms = 15
      integer(INT4) :: i, nexp
!
!  Start of executable code
!  ------------------------
!
      x = xvalue

      if( x < ZERO ) then
         print '(a)',  ' '
         print '(a)',  'DEBYE1 - Fatal error!'
         print '(a)',  '  Argument X < 0.'
         debye1 = ZERO
      else if( x < xlow ) then
         debye1 = ( ( x - nine ) * x + thirt6 ) / thirt6
      else if( x <= FOUR ) then
         t = ( ( x * x / eight ) - HALF ) - HALF
         debye1 = cheval( nterms, adeb1, t ) - QUART * x
!                 ======    i       i    i
      else

         debye1 = ONE / ( x * debinf )
         if( x < xlim ) then
            expmx = exp ( -x )
            if( xupper < x ) then
               debye1 = debye1 - expmx * ( ONE + ONE / x )
            else
               sum1 = ZERO
               rk = aint ( xlim / x )
               nexp = int ( rk )
               xk = rk * x
               do i = nexp, 1, -1
                  t = ( ONE + ONE / xk ) / rk
                  sum1 = sum1 * expmx + t
                  rk = rk - ONE
                  xk = xk - x
               end do
               debye1 = debye1 - sum1 * expmx
            end if
         end if
      end if
        
      return

   end function debye1
!
!==========================================================================================
!
   function cheval( n, a, t )
!                   i  i  i
! 
!    Evaluates a Chebyshev series.
!
!    This function evaluates a Chebyshev series, using the Clenshaw method with Reinsch 
!    modification, as analysed in the paper by Oliver.
!
!  Author:
!
!    Allan McLeod,
!    Department of Mathematics and Statistics,
!    Paisley University, High Street, Paisley, Scotland, PA12BE
!    macl_ms0@paisley.ac.uk
!
!  Reference:
!
!    Allan Mcleod,
!    Algorithm 757, MISCFUN: A software package to compute uncommon  special functions,
!    ACM Transactions on Mathematical Software, 
!    Volume 22, Number 3, September 1996, pages 288-301.
!
!    J Oliver,
!    "An error analysis of the modified Clenshaw method for
!    evaluating Chebyshev and Fourier series",
!    J.I.M.A., Volume 20, 1977, pages 379-391.
!
!  Arguments:
!
      integer(INT4), intent(in) :: n
      real(REAL8), intent(in) :: a(0:n)
      real(REAL8), intent(in) :: t
      real(REAL8) :: cheval                 ! returning value
!
!  Local variables
!
      real(REAL8) :: d1, d2, tt, u0, u1, u2
      integer(INT4) ::  i
      real(REAL8), parameter :: test = 0.6D+00
!
!  Start of executable code
!  ------------------------
!
      u1 = ZERO
!
!  T <= -0.6, Reinsch modification.
!
      if( t <= -test ) then

         d1 = ZERO
         tt = ( t + HALF ) + HALF
         tt = tt + tt
           
         do i = n, 0, -1
            d2 = d1
            u2 = u1
            d1 = tt * u2 + a(i) - d2
            u1 = d1 - u2
         end do
           
         cheval = ( d1 - d2 ) / TWO
!
!  -0.6 < T < 0.6, Standard Clenshaw method.
!
      else if( t < test ) then

         u0 = ZERO
         tt = t + t
          
         do i = n, 0, -1
            u2 = u1
            u1 = u0
            u0 = tt * u1 + a(i) - u2
         end do
           
         cheval = ( u0 - u2 ) / TWO
!
!  0.6 <= T, Reinsch modification.
!
      else
           
         d1 = ZERO
         tt = ( t - HALF ) - HALF
         tt = tt + tt
           
         do i = n, 0, -1
            d2 = d1
            u2 = u1
            d1 = tt * u2 + a(i) + d2
            u1 = d1 + u2
         end do
          
         cheval = ( d1 + d2 ) / TWO
           
      end if
        
      return

   end function cheval
!
!===========================================================================================
!
end module Math_Debye_mod
