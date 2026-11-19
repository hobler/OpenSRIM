module Sn_Cmangm_mod
!
!  Compute the scattering angle in the center-of-mass system and the time integral 
!  using Gauss-Mehler quadrature 
!
!  Used Modules:
!  ============
!
   use Const_mod,       only: INT4, REAL8,  &
                              EPSMA, ZERO, ONE, TWO, FOUR, PI

   use Param_Sn_mod,    only: sn

   use Sn_Screen_mod,   only: screen, get_rmax
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
   public :: initCmangm             ! initialize weights for quadrature
!            ==========
   public :: cmangm                 ! calculate scattering angle
!            ======
!
!  Private Data:
!  ============
!
   real(REAL8), allocatable :: b(:)          ! (nmehl) weights for quadrature
!
!=======================================================================
!
contains
!
!=======================================================================
!
   ! Allocate the quadrature nodes for sn%nmehl once before calling cmangm.
   ! This routine does not deallocate an existing rule; comparison drivers use
   ! separate processes when changing orders to respect that lifecycle.
   subroutine initCmangm( )
!             ==========
!
!     Initialize weights for Gauss Mehler quadrature
!
!  Arguments:
!
!     None
!
!  Local Variables:
!
      integer(INT4) :: i
!
!  Start of Executable Code:
!  ------------------------
!
      allocate( b( sn%nmehl ) )

      do i = 1, sn%nmehl
         b(i) = cos( (TWO*i-ONE) / (FOUR*sn%nmehl) * PI )
      end do

   end subroutine initCmangm
!
!=======================================================================
!
   ! Return theta and tau using the initialized native Gauss-Mehler rule.
   ! R0 is an input estimate: refinement changes a local copy only. Closely
   ! spaced nodes near the turning point need a more accurate root than the
   ! general apsis stopping criterion. Nonconvergence terminates the native run.
   subroutine cmangm( theta, tau, e, p, r0, ielem1, ielem2 )
!             ======    o     o   i  i  i     i       i
!
!     Calculate the scattering angle in the center-of-mass system and the
!     time integral using Gauss-Mehler quadrature
!
!  Arguments:
!
      real(REAL8), intent(out) :: theta         ! deflection angle in CM system
      real(REAL8), intent(out) :: tau           ! time integral [rnorm]
      real(REAL8), intent(in) :: e              ! ion energy [enorm]
      real(REAL8), intent(in) :: p              ! impact parameter [rnorm]
      real(REAL8), intent(in) :: r0             ! scaled distance of closest approach [rnorm]
      integer(INT4), intent(in) :: ielem1       ! index of projectile
      integer(INT4), intent(in) :: ielem2       ! index of recoil
!
!  Local Variables:
!
      integer(INT4) :: i, iteration
      real(REAL8) :: sumth, sumtau, r, phi, phi1, g, gg, rturn, correction
!
!  Start of Executable Code:
!  ------------------------
!
      !
      !new
      if( p >= get_rmax(ielem1, ielem2) ) then
         theta = ZERO
         tau = ZERO
         return
      end if

      ! Mehler nodes approach the turning point closely. Refine the supplied
      ! apsis before evaluating the square-root singularity.
      rturn = r0
      do iteration = 1, 20
         call screen(rturn, ielem1, ielem2, phi, phi1)
         !
         !new
         correction = ((rturn-p)*(ONE+p/rturn) - phi/e) / &
                      (ONE + (p/rturn)**2 - phi1/e)
         rturn = rturn - correction

         !
         !new
         if( abs(correction) <= 1e-12_REAL8*rturn ) exit
      end do
      if( iteration > 20 ) error stop 'cmangm: apsis refinement did not converge'

      sumth = ZERO
      sumtau = ZERO
      do i = 1, sn%nmehl
         r = rturn / b(i)**2
         call screen( r, ielem1, ielem2, phi, phi1 )
!             ======  i    i       io    o
         g = sqrt (ONE - phi/(r*e) - (p/r)**2)
         gg = sqrt (ONE - (p/r)**2)
         sumth = sumth + b(i)*b(sn%nmehl-i+1) / g
         sumtau = sumtau + b(sn%nmehl-i+1)/b(i)**3 * (1/g - 1/gg)
      end do

      theta = PI - 2*p/rturn * PI/sn%nmehl * sumth
                                                ! avoid unphysical values due to numerical
      theta = max( EPSMA, min( PI, theta) )     !    errors

      tau = sqrt( (rturn-p)*(rturn+p) ) - rturn * PI/sn%nmehl * sumtau

   end subroutine cmangm
!
!=======================================================================
!
end module Sn_Cmangm_mod
