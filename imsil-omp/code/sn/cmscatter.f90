module Sn_CMscatter_mod
!
!  Calculate the scattering angle in the center-of-mass system and the time integral
!
!  Used Modules:
!  ============
!
   use Const_mod,             only: INT4, REAL8,  &
                                    ZERO, ONE, HALF, TWO, FOUR, SQRT2, PI, BIG

   use Param_Sn_mod,          only: sn
   use Sn_Apsis_mod,          only: apsis
   use Sn_Screen_mod,         only: screen, get_rmax
   use Math_Quadrature_mod,   only: quadLegendre, quadMehlerSymm, quadMehlerLin
!
!  Default Type and Access:
!  =======================
!
   implicit none
   private
   save
!
!  Public Data Types
!  =================
!
!     None
!
!  Public Procedures
!  =================
!
   public :: cmScatterLegendre                  ! LegendreChiCorr angle and Legendre time integral
   public :: calcPhiChi                         ! screening and chi used by integrands
   public :: cmAngle                            ! calculate CM angle
!            =======
   public :: cmang_magic                        ! calculate CM angle using Biersack's Magic
!            ===========                        !    formula
!
!  Private Data
!  ============
!
   integer(INT4) :: ielem1_save = 0             ! index of projectile
   integer(INT4) :: ielem2_save = 0             ! index of recoil
   integer(INT4) :: nabsc_save = 0              ! number of abscissas
   real(REAL8) :: e_save = ZERO                 ! energy [enorm]
   real(REAL8) :: p_save = ZERO                 ! impact parameter [rnorm]
   real(REAL8) :: r0_save                       ! apsis [rnorm]
   real(REAL8) :: phi0_save                     ! screening function at r0
   real(REAL8) :: umax_save                     ! upper integration boundary (depending on rmax)
   real(REAL8) :: pinum_legendre
   real(REAL8) :: arcsin_legendre
   real(REAL8) :: arcsin_save
   real(REAL8) :: arccos_save
!$omp threadprivate( ielem1_save, ielem2_save, nabsc_save, e_save, p_save )
!$omp threadprivate( r0_save, phi0_save, umax_save, arcsin_save, arccos_save )
!$omp threadprivate( pinum_legendre, arcsin_legendre )
!
!==============================================================================================
!
contains
!
!==============================================================================================
!
   subroutine calcPhiChi(u, r0, ielem1, ielem2, phi, chi)
      ! Production chi shared by quadrature and callers inspecting the integrand.
      real(REAL8), intent(in) :: u, r0
      integer(INT4), intent(in) :: ielem1, ielem2
      real(REAL8), intent(out) :: phi, chi
      real(REAL8) :: phi0, dphi0
      call screen(r0, ielem1, ielem2, phi0, dphi0)
      call screen(r0/(ONE-u**2), ielem1, ielem2, phi)
      if(u < 3e-4_REAL8) then
         chi = phi0-r0*dphi0
      else
         chi = (phi0-phi*(ONE-u**2))/u**2
      end if
   end subroutine calcPhiChi

   real(REAL8) function f_pinum_legendre( u )
!                       ================  i
!
!     PI calculated by Legendre quadrature
!
!  Arguments:
!
      real(REAL8), intent(in) :: u                 ! integration variable (0<u<1)
!
!  Local variables
!
!     None
!
!  Start of Executable Code
!  ------------------------
!
      f_pinum_legendre = FOUR / sqrt(TWO-u**2)

   end function f_pinum_legendre
!
!==============================================================================================
!
   real(REAL8) function f_theta_chi_legendre( u )
!                       ====================  i
!
!     Integrand of Legendre quadrature for the calculation of PI minus the scattering angle
!     in the center-of-mass-system (inaccurate for large impact parameters, but can be
!     corrected) - Version using the chi(u) function.
!
!     Energy and impact parameter are passed by module variables
!
!  Arguments:
!
      real(REAL8), intent(in) :: u                 ! integration variable (0<u<1)
!
!  Local variables
!
      real(REAL8) :: rho, phi, chi
!
!  Start of Executable Code
!  ------------------------
!
      rho = e_save * p_save**2 / r0_save
      call calcPhiChi(u, r0_save, ielem1_save, ielem2_save, phi, chi)

      f_theta_chi_legendre = FOUR * sqrt(rho/(rho*(TWO-u**2)+chi))

   end function f_theta_chi_legendre
!
!==============================================================================================
!
   ! Evaluate four times the free-minus-interacting angle integrand.
   ! Uses the saved collision state set by cmAngle; this is not a standalone API.
   ! Rationalization retains weak deflections near the finite-range cutoff.
   real(REAL8) function f_theta_deflection_legendre(u)
      ! Rationalized free-minus-interacting integrand. The small-u limit
      ! avoids cancellation of screening values close to the apsis.
      real(REAL8), intent(in) :: u
      real(REAL8) :: phi, chi, base, h, interaction
      call calcPhiChi(u, r0_save, ielem1_save, ielem2_save, phi, chi)
      base = sqrt(TWO-u**2)
      interaction = r0_save/e_save*chi
      h = sqrt(interaction+p_save**2*base**2)
      f_theta_deflection_legendre = FOUR*interaction/(base*h*(h+p_save*base))
   end function f_theta_deflection_legendre

   real(REAL8) function f_tau_legendre(u)
      ! Time-delay integrand in reduced length units, with r=r0/(1-u**2).
      real(REAL8), intent(in) :: u
      real(REAL8) :: phi, chi, h
      call calcPhiChi(u, r0_save, ielem1_save, ielem2_save, phi, chi)
      h = sqrt(r0_save/e_save*chi + p_save**2*(TWO-u**2))
      f_tau_legendre = (p_save**2 + r0_save/e_save*phi/(ONE-u**2)) / &
                       (h*(ONE+u/r0_save*h))
   end function f_tau_legendre

   subroutine cmScatterLegendre(e, p, ielem1, ielem2, theta, tau, nabsc)
      ! Public joint API: angle uses LegendreChiCorr, tau uses the same native
      ! Gauss-Legendre rule and apsis. No Mehler or scattering tables are used.
      ! At p=0 adopt the standalone Python API's tau=0 convention.
      real(REAL8), intent(in) :: e, p
      integer(INT4), intent(in) :: ielem1, ielem2, nabsc
      real(REAL8), intent(out) :: theta, tau
      real(REAL8) :: rmax, sp
      if(e <= ZERO .or. p < ZERO .or. nabsc < 1) error stop 'cmScatterLegendre: invalid input'
      rmax = get_rmax(ielem1, ielem2)
      if(p >= rmax) then
         theta=ZERO
         tau=ZERO
         return
      end if
      call cmAngle(e, p, ielem1, ielem2, theta, method='LegendreChiCorr', nabsc=nabsc)
      if(p == ZERO) then
         tau=ZERO
         return
      end if
      tau = r0_save
      if(rmax /= BIG) then
         sp = sqrt((rmax-p)*(rmax+p))
         ! Rationalized finite-range tail, including collapsed intervals.
         tau = tau-p**2/(rmax+sp)
      end if
      if(umax_save > ZERO) tau = tau-TWO*quadLegendre(f_tau_legendre, nabsc, umax_save)
   end subroutine cmScatterLegendre

   real(REAL8) function f_theta_g_legendre( u )
!                       ==================  i
!
!     Integrand of Legendre quadrature for the calculation of PI minus the scattering angle
!     in the center-of-mass-system (inaccurate for large impact parameters, but can be
!     corrected) - Version using the g(u) function
!
!     Energy and impact parameter are passed by module variables
!
!  Arguments:
!
      real(REAL8), intent(in) :: u                 ! integration variable (0<u<1)
!
!  Local variables
!
      real(REAL8) :: phi, r, g
!
!  Start of Executable Code
!  ------------------------
!
      r = r0_save / (ONE-u**2)
      call screen( r, ielem1_save, ielem2_save, phi )
!          ======  i       i            i        o
      g = sqrt(ONE - phi/(e_save*r) - (p_save/r)**2)

      f_theta_g_legendre = FOUR * p_save/r0_save * u / g

   end function f_theta_g_legendre
!
!==============================================================================================
! below are old integrands
!==============================================================================================
!
   real(REAL8) function f_theta_legendre( u )
!                       ================  i
!
!     Integrand of Legendre quadrature for the calculation of PI minus the scattering angle
!     in the center-of-mass-system (inaccurate for large impact parameters)
!
!     Energy and impact parameter are passed by module variables
!
!  Arguments:
!
      real(REAL8), intent(in) :: u                 ! integration variable (0<u<1)
!
!  Local variables
!
      real(REAL8) :: rho, phi, chi
!
!  Start of Executable Code
!  ------------------------
!
      rho = e_save * p_save**2 / r0_save
      call calcPhiChi(u, r0_save, ielem1_save, ielem2_save, phi, chi)

      f_theta_legendre = FOUR * (1 / sqrt(TWO-u**2) - sqrt(rho/(rho*(TWO-u**2)+chi)))

   end function f_theta_legendre
!
!==============================================================================================
!
   real(REAL8) function f_theta_legendre_acc( u )
!                       ====================  i
!
!     Integrand of Legendre quadrature for the calculation of PI minus the scattering angle
!     in the center-of-mass-system (inaccurate for large impact parameters)
!
!     Energy and impact parameter are passed by module variables
!
!  Arguments:
!
      real(REAL8), intent(in) :: u                 ! integration variable (0<u<1)
!
!  Local variables
!
      real(REAL8) :: rho, phi, chi, s1, s2
!
!  Start of Executable Code
!  ------------------------
!
      rho = e_save * p_save**2 / r0_save
      call calcPhiChi(u, r0_save, ielem1_save, ielem2_save, phi, chi)
      s1 = sqrt(rho*(TWO-u**2))
      s2 = sqrt(rho*(TWO-u**2)+chi)

      f_theta_legendre_acc = FOUR * (chi / (sqrt(TWO-u**2) * (s1+s2) * s2) )

   end function f_theta_legendre_acc
!
!==============================================================================================
!
   real(REAL8) function f_pi_theta_legendre( u )
!                       ===================  i
!
!     Integrand of Legendre quadrature for the calculation of PI minus the scattering angle
!     in the center-of-mass-system (inaccurate for large impact parameters)
!
!     Energy and impact parameter are passed by module variables
!
!  Arguments:
!
      real(REAL8), intent(in) :: u                 ! integration variable (0<u<1)
!
!  Local variables
!
      real(REAL8) :: rho, phi, chi
!
!  Start of Executable Code
!  ------------------------
!
      rho = e_save * p_save**2 / r0_save
      call calcPhiChi(u, r0_save, ielem1_save, ielem2_save, phi, chi)

      f_pi_theta_legendre = FOUR * sqrt(rho/(rho*(TWO-u**2)+chi))

   end function f_pi_theta_legendre
!
!==============================================================================================
!
   real(REAL8) function f_pi_theta_mehler( u )
!                       =================  i
!
!     Integrand of Gauss-Mehler quadrature for the calculation of PI minus the scattering
!     angle in the center-of-mass-system (inaccurate for large impact parameters)
!
!     Energy and impact parameter are passed by module variables
!
!  Arguments:
!
      real(REAL8), intent(in) :: u                 ! integration variable (0<u<1)
!
!  Local variables
!
      real(REAL8) :: rho, phi, varphi
!
!  Start of Executable Code
!  ------------------------
!
      rho = e_save * p_save**2 / r0_save
      call screen( r0_save/u, ielem1_save, ielem2_save, phi )
!          ======      i           i            i        o
      varphi = (phi0_save - phi*u) / (1-u**2)

      f_pi_theta_mehler = TWO * sqrt(rho/(rho+varphi))

   end function f_pi_theta_mehler
!
!==============================================================================================
!
   real(REAL8) function f_theta_mehler( u )
!                       ==============  i
!
!     Integrand of Gauss-Mehler quadrature for the calculation of the scattering angle in the
!     center-of-mass-system (less accurate version)
!
!     Energy and impact parameter are passed by module variables
!
!  Arguments:
!
      real(REAL8), intent(in) :: u                 ! integration variable (0<u<1)
!
!  Local variables
!
      real(REAL8) :: rho, phi, varphi
!
!  Start of Executable Code
!  ------------------------
!
      rho = e_save * p_save**2 / r0_save
      call screen( r0_save/u, ielem1_save, ielem2_save, phi )
!          ======      i           i            i        o
      varphi = (phi0_save - phi*u) / (1-u**2)

      f_theta_mehler = TWO * (ONE - sqrt(rho/(rho+varphi)))

   end function f_theta_mehler
!
!==============================================================================================
!
   real(REAL8) function f_theta_mehler_acc( u )
!                       ==================  i
!
!     Integrand of Gauss-Mehler quadrature for the calculation of the scattering angle in the
!     center-of-mass-system (accurate version)
!
!     Energy and impact parameter are passed by module variables
!
!  Arguments:
!
      real(REAL8), intent(in) :: u                 ! integration variable (0<u<1)
!
!  Local variables
!
      real(REAL8) :: rho, phi, varphi
!
!  Start of Executable Code
!  ------------------------
!
      rho = e_save * p_save**2 / r0_save
      call screen( r0_save/u, ielem1_save, ielem2_save, phi )
!          ======      i           i            i        o
      varphi = (phi0_save - phi*u) / (1-u**2)

      f_theta_mehler_acc = TWO * varphi / (rho+varphi + sqrt(rho*(rho+varphi)))

   end function f_theta_mehler_acc
!
!==============================================================================================
!
   ! Return the deflection angle in radians (not its complement).
   ! METHOD and NABSC select the native quadrature; finite-range p >= rmax
   ! returns zero before apsis evaluation. A collapsed integration interval
   ! uses its geometric limit before any quadrature normalization is divided.
   subroutine cmAngle( e, p, ielem1, ielem2, theta, method, nabsc )
!             =======  i  i    i       i       o      i       i
!
!     Calculate the scattering angle in the center-of-mass system
!
!  Arguments:
!
      real(REAL8), intent(in) :: e                 ! energy [enorm]
      real(REAL8), intent(in) :: p                 ! impact parameter [rnorm]
      integer(INT4), intent(in) :: ielem1          ! index of projectile
      integer(INT4), intent(in) :: ielem2          ! index of recoil
      real(REAL8), intent(out) :: theta            ! scattering angle [rad]
      character(*), intent(in) :: method ! integration method
      integer(INT4), intent(in) :: nabsc ! number of abscissa for integration
!
!  Local Variables:
!
      logical :: updated
      real(REAL8) :: rmax, sp, s0, geometric
!
!  Start of Executable Code:
!  ------------------------
!
      rmax = get_rmax( ielem1, ielem2 )
!            ========    i       i
      if( p >= rmax ) then
         theta = ZERO
         return
      end if
!
!  Calculate apsis and value of the screening function at the apsis
!
      updated = .false.
      if( e /= e_save .or. p /= p_save .or.  &
          ielem1 /= ielem1_save .or. ielem2 /= ielem2_save) then
         updated = .true.
         e_save = e
         p_save = p
         ielem1_save = ielem1
         ielem2_save = ielem2

         call apsis( e, p, ielem1, ielem2, r0_save )
!             =====  i  i    i       i        o
         call screen( r0_save, ielem1, ielem2, phi0_save )
!             ======     i       i       i         o

         if( rmax == BIG ) then                ! meaning rmax=infinity
            umax_save = ONE
            arcsin_save = PI
            arccos_save = PI
         else
            umax_save = sqrt((rmax-r0_save)/rmax)
            s0 = sqrt((rmax-r0_save)*(rmax+r0_save))
            sp = sqrt((rmax-p_save)*(rmax+p_save))
            arcsin_save = TWO*atan2(s0,r0_save)
            arccos_save = TWO*atan2(sp,p_save)
         end if
      end if
!
!  Calculate numerical value of pi
!
      if( umax_save == ZERO ) then
         theta = arccos_save
         return
      end if

      if( method(1:8) == 'Legendre' .and. (nabsc /= nabsc_save .or. updated) ) then
         nabsc_save = nabsc
         if( method == 'LegendreAcc' ) then
            pinum_legendre = quadLegendre( f_pinum_legendre, nabsc, ONE )
!                            ============          p           i     i
         else
            arcsin_legendre = quadLegendre( f_pinum_legendre, nabsc, umax_save )
!                             ============          p           i        i
         end if
      end if
!
!  Perform integration
!
      select case( method )

         case( 'MehlerSymm' )             ! recommended for calculating theta alone
            theta = quadMehlerSymm( f_theta_mehler, nabsc )
!                   ==============         p              i

         case( 'LegendreChiCorr' )        ! needs one iteration more to reach same accuracy
            if(rmax == BIG) then
               theta = quadLegendre(f_theta_chi_legendre, nabsc, umax_save)
               theta = arccos_save-arcsin_save/arcsin_legendre*theta
            else
               sp = sqrt((rmax-p_save)*(rmax+p_save))
               s0 = sqrt((rmax-r0_save)*(rmax+r0_save))
               geometric = TWO*asin((r0_save-p_save)*(r0_save+p_save)/(r0_save*sp+p_save*s0))
               theta = quadLegendre(f_theta_deflection_legendre, nabsc, umax_save)
               theta = geometric+arcsin_save/arcsin_legendre*theta
            end if

         case( 'MehlerSymmInacc' )        ! almost same accuracy and expense as 'MehlerSymm'
            theta = quadMehlerSymm( f_pi_theta_mehler, nabsc )
!                   ==============          p            i
            theta = PI - theta

         case( 'MehlerSymmAcc' )          ! almost same accuracy and expense as 'MehlerSymm'
            theta = quadMehlerSymm( f_theta_mehler_acc, nabsc )
!                   ==============         p              i

         case( 'LegendreInacc' )          ! almost same as 'LegendreChiCorr'
            theta = quadLegendre( f_pi_theta_legendre, nabsc, umax_save )
!                   ============           p             i        i
            theta = arccos_save - arcsin_save/arcsin_legendre * theta

         case( 'LegendreAcc' )            ! almost same as 'LegendreChiCorr'
                                          ! Not sure if this is correct for umax/=1
            theta = quadLegendre( f_theta_legendre_acc, nabsc, umax_save ) * PI/pinum_legendre
!                   ============            p             i        i

         case( 'LegendreChi' )            ! less accurate than 'LegendreChiCorr' at large p
                                          ! Not sure if this is correct for umax/=1
            theta = quadLegendre( f_theta_chi_legendre, nabsc, umax_save )
!                   ============            p             i        i
            theta = PI - theta

         case( 'Legendre' )               ! less accurate than 'LegendreChiCorr' at small p
                                          ! Not sure if this is correct for umax/=1
            theta = quadLegendre( f_theta_legendre, nabsc, umax_save )
!                   ============          p           i        i

         case( 'LegendreGCorr' )          ! less accurate than 'LegendreChiCorr' for large nabsc
                                          ! WARNING: This may sometimes return NaN
            theta = quadLegendre( f_theta_g_legendre, nabsc, umax_save )
!                   ============          p             i        i
            theta = arccos_save - arcsin_save/arcsin_legendre * theta

         case( 'MehlerLinAcc' )           ! very inaccurate
            theta = quadMehlerLin( f_theta_mehler_acc, nabsc )
!                   =============         p              i

         case( 'MehlerLin' )              ! almost same as 'MehlerLinAcc'
            theta = quadMehlerLin( f_theta_mehler, nabsc )
!                   =============         p              i

         case( 'MehlerLinInacc' )         ! (Eckstein) worst; VERY inaccurate for large p
            theta = quadMehlerLin( f_pi_theta_mehler, nabsc )
!                   =============          p            i
            theta = PI - theta

         case default
            error stop 'cmAngle: unknown integration method'

      end select

   end subroutine cmAngle
!
!==============================================================================================
!
   subroutine cmang_magic( e, p, ielem1, ielem2, theta )
!             ===========  i  i    i       i       o
!
!     Calculate the scattering angle in the CM system according to Biersack's Magic formula
!
!     The c values are from the [ZBL85] book
!
!  Arguments:
!
      real(REAL8), intent(in) :: e                 ! reduced energy (enorm)
      real(REAL8), intent(in) :: p                 ! reduced impact parameter (rnorm)
      integer(INT4), intent(in) :: ielem1          ! index of projectile
      integer(INT4), intent(in) :: ielem2          ! index of recoil
      real(REAL8), intent(out) :: theta            ! scattering angle in CM system (rad)
!
!  Local Variables:
!
      real(REAL8) :: r0, phi, phi1, rho, sqrte
      real(REAL8) :: c(5), alpha, beta, gamma, a, g, delta, costhh
!
!  Start of Executable Code:
!  ------------------------
!
      select case( sn%scrfun(ielem2,ielem1) )
         case( 'KrC' )
            c(:) = [ 0.7887, 0.01166, 0.006913, 17.16, 10.79 ]   ! different values in [Eck91]
         case( 'Moliere' )
            c(:) = [ 0.6743, 0.009611, 0.005175, 10., 6.314 ]
         case default
            c(:) = [ 0.99229, 0.011615, 0.007122, 14.813, 9.3066 ]
      end select
!
!  Calculate distance of closest approach, its derivatives, and the screening function at r0
!
      call apsis( e, p, ielem1_save, ielem2_save, r0 )
!          =====  i  i       i            i       o
      call screen( r0, ielem1_save, ielem2_save, phi, phi1 )
!          ======  i        i            i        o    o
!
!  Evaluate Magic formula
!
      rho = 2 * r0 * (e*r0 - phi) / (phi - r0*phi1)

      sqrte = sqrt(e)
      alpha = 1 + c(1)/sqrte
      beta = (c(2)+sqrte) / (c(3) + sqrte)
      gamma = (c(4)+e) / (c(5)+e)

      a = 2 * alpha * e * p**beta
      g = gamma / (sqrt(1+a**2) - a)

      delta = a * (r0-p) / (1+g)

      costhh = (p+rho+delta) / (r0+rho)

      theta = 2 * acos(costhh)

   end subroutine cmang_magic
!
!=======================================================================
!
end module Sn_CMscatter_mod
