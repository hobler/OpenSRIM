module Sn_Ndira_mod
!
!  Compute the new direction after a deflection by a given angle with random azimuthal angle
!
!  Used Modules:
!  ============
!
   use Const_mod,       only: INT4, REAL8,  &
                              ZERO, ONE, TWO, THREE, HALF, PI
   use Vector_mod,      only: abs

   use Random_mod,      only: rannum
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
   public :: ndira
!            =====
!
!=======================================================================
!
contains
!
!=======================================================================
!
   subroutine ndira( dir, cospsi, sinpsi, dirn )
!             =====   i     i       i      o
!
!     compute the new direction after a deflection by a given angle with 
!     random azimuthal angle
!
!     the azimuthal angle is chosen randomly within the routine
!
!     dirxn**2 + diryn**2 + dirzn**2 = 1 is enforced by proper scaling
!
!  Arguments:
!
      real(REAL8), intent(in) :: dir(3)         ! direction of motion
      real(REAL8), intent(in) :: cospsi         ! cosine of deflection angle
      real(REAL8), intent(in) :: sinpsi         ! sine of deflection angle
      real(REAL8), intent(out) :: dirn(3)       ! new direction of motion
!
!  Local Variables:
!
!      real(REAL8), parameter :: AMIN = 1.E-10

      integer(INT4) :: i, j, k
      real(REAL8) :: rand, fi, cosfi, sinfi, cosalf, sinalf 
      real(REAL8) :: cosphi, sinphi
      real(REAL8) :: h0, h1, h2, dirnn
!
!  Start of Executable Code:
!  ------------------------
!
!  Select the azimuthal angle randomly
!
      call rannum( rand )
!          ======   o
      fi = TWO*PI * rand
      cosfi = cos( fi )
      sinfi = sin( fi )
!
!  Convert direction vector to polar angles
!
      !i = 3   ! the original version
      !j = 2
      !k = 1
      k = minloc( dir(:), dim=1 )   ! make k point to the smallest dir(:) so sinalf > sqrt(2/3)
      i = mod( k, 3 ) + 1
      j = mod( i, 3 ) + 1
      cosalf = dir(k)
      sinalf = sqrt( dir(i)**2 + dir(j)**2 )
      !if( sinalf > AMIN ) then
         cosphi = dir(i) / sinalf
         sinphi = dir(j) / sinalf
      !else
      !   cosphi = ONE
      !   sinphi = ZERO
      !end if
!
!  Calculate new direction
!
      h0 = sinpsi*cosfi
      h1 = sinalf*cospsi + h0*cosalf
      h2 = sinpsi*sinfi
      dirn(i) = h1*cosphi - h2*sinphi
      dirn(j) = h1*sinphi + h2*cosphi
      dirn(k) = cosalf*cospsi - h0*sinalf

      !dirn(:) = TWO/(ONE+abs2(dirn(:))) * dirn(:)  ! for unknown reason below is faster
!                        ====
      dirnn = abs( dirn(:) )                    ! scale direction vector to length 1
      dirn(:) = dirn(:) / dirnn

   end subroutine ndira
!
!=======================================================================
!
end module Sn_Ndira_mod
