module Se_OR_mod
!
!  Handle the impact parameter dependent part of electronic stopping according to Oen and
!  Robinson
!
!  Used Modules:
!  ============
!
   use Const_mod,       only: INT4, REAL8,  &
                              ZERO, ONE, PI

   use Param_Setup_mod, only: nelem
   use Param_Elem_mod,  only: elem
   use Param_Se_mod,    only: se
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
   public :: evalSeOR               ! evaluate local part of electronic energy loss according
!            ========               !    to Oen and Robinson
   public :: evalCorXnlOR           ! evaluate xnl correction due to neglect of local part
!            ============           !    of originating from p > pmax
!
!=======================================================================
!
contains
!
!=======================================================================
!
   subroutine evalSeOR( fac, p, ascr )
!             ========   o   i   i
!
!     Evaluate the local electronic energy loss factor according to Oen and Robinson.
!     The local electronic energy loss factor has to be multiplied with the stopping power to
!     obtain the local electronic energy loss.
!
!  Arguments:
!
      real(REAL8), intent(out) :: fac           ! local electronic energy loss factor
      real(REAL8), intent(in) :: p              ! impact parameter
      real(REAL8), intent(in) :: ascr           ! screening length
!
!  Local variables:
!
!     None
!
!  Start of Executable Code:
!  ------------------------
!
      fac = ONE / ( 2 * PI * ascr**2 ) * exp( - p / ascr )

   end subroutine evalSeOR
!
!=======================================================================
!
   subroutine evalCorXnlOR( corxnl, pmax, ascr )
!             ============    o      i     i
!
!     Evaluate the correction to the nonlocal electronic energy loss due to the maximum impact
!     parameter
!
!  Arguments:
!
      real(REAL8), intent(out) :: corxnl        ! correction factor for nonlocal fraction
      real(REAL8), intent(in) :: pmax           ! maximum impact parameter used when
                                                !    calling the subroutine
      real(REAL8), intent(in) :: ascr           ! screening length
!
!  Local variables:
!
      real(REAL8) :: pm
!
!  Start of Executable Code:
!  ------------------------
!
      pm = pmax / ascr
      corxnl = (ONE + pm) * exp(-pm)

   end subroutine evalCorXnlOR
!
!=======================================================================
!
end module Se_OR_mod
