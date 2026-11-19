module Math_Bessel_mod
!
!  Calculate Bessel function
!
!  Used Modules:
!  ============
!
   use Const_mod, only: INT4, REAL8, &
                        ZERO, ONE, TWO, HALF
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
   public :: mbes1                  ! calculates the modified Bessel function of first
!            =====                  !   order

contains
!
!=======================================================================
!
   subroutine mbes1( arg, value )
!             =====   i     o
!
!     Calculates the modified Bessel function of order 1.
!
!     arg has to be >0.
!
!     formulae for mbes1 and bes1 have been taken from Abramowitz Stegun    
!     (formulae 9.8.7, 9.8.8 and 9.8.3, respectively).
!
!  Arguments:
!
      real(REAL8), intent(in) :: arg            ! argument of the function
      real(REAL8), intent(out) :: value         ! result of the function
!
!  Local Variables:
!
      real(REAL8) :: argb, argb2, bes1
      real(REAL8) :: argmb, argmb2, besm1
!
!  Start of Executable Code
!  ------------------------
!
!  Check validity of the function's argument
!
      if( arg <= ZERO )  goto 900                   
!
!  Calculate modified Bessel function
!
      if( arg <= TWO ) then

         argb = arg / 3.75                      ! calculate Bessel function of order 1 of arg
         argb2 = argb**2                        !    using Horner's scheme
         bes1 = ZERO
         bes1 = (0.00032411*argb2 + 0.00301532) * argb2
         bes1 = ((bes1+0.02658733)*argb2 + 0.15084934) * argb2
         bes1 = ((bes1+0.51498869)*argb2 + 0.87890594) * argb2 + HALF
         bes1 = bes1 * arg 

         argmb = arg / TWO                      ! calculate modified Bessel function of order 1
         argmb2 = argmb**2                      !    of arg using Horner's scheme
         besm1 = (-0.00004686*argmb2 - 0.00110404) * argmb2
         besm1 = ((besm1-0.01919402)*argmb2 - 0.18156897) * argmb2
         besm1 = ((besm1-0.67278579)*argmb2 + 0.15443144) * argmb2 + ONE
         besm1 = besm1/arg + bes1 * log(argmb)

         value = besm1                    

      else 
                                                ! for arg>2 Bessel function is not needed
         argmb = TWO / arg
         besm1 = ((-0.00068245*argmb + 0.00325614) * argmb - 0.00780353) * argmb
         besm1 = ((besm1+0.01504268)*argmb - 0.03655620) * argmb
         besm1 = (besm1+0.23498619)*argmb + 1.25331414
         besm1 = besm1 / (exp(arg)*sqrt(arg))

         value = besm1          

      end if      

      return

900   error stop 'mbes1: argument must be >0' 

   end subroutine mbes1
!
!==============================================================================
!
end module Math_Bessel_mod
