module Sn_Screen_old_mod
!
!  Calculate the screening function and its derivative using the ZBL potential
!
!  Used Modules:
!  ============
!
   use Const_mod,    only: INT4, REAL8,  &
                           ZERO

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
   public :: screen                 ! calculate screening function and
!            ======                 ! derivative
!
!  Private Data:
!  ============
!
   integer(INT4), parameter :: NCOEF = 4
   real(REAL8), dimension(NCOEF), parameter :: A = [ 0.1818, 0.50983, 0.2802, 0.02817 ]
   real(REAL8), dimension(NCOEF), parameter :: B = [ 3.2   , 0.9423 , 0.4028, 0.2016  ]
!
!=======================================================================
!
contains
!
!=======================================================================
!
   subroutine screen( r, phi, phi1)
!             ======  i   o    o
!
!     Calculate the screening function and its derivative using the ZBL 
!     potential
!
!  Arguments:
!
      real(REAL8), intent(in) :: r              ! interatomic distance [rnorm]
      real(REAL8), intent(out) :: phi           ! screening function
      real(REAL8), intent(out) :: phi1          ! derivative of screening function
!
!  Local Variables:
!
      integer(INT4) :: i
      real(REAL8) :: fexp
!
!  Start of Executable Code:
!  ------------------------
!
      phi = ZERO
      phi1 = ZERO
      do i = 1, NCOEF
         fexp = exp( -B(i)*r )
         phi = phi + A(i)*fexp
         phi1 = phi1 - A(i)*B(i)*fexp
      end do

   end subroutine screen
!
!=======================================================================
!
end module Sn_Screen_mod
