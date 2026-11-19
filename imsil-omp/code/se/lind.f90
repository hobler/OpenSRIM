module Se_Lind_mod
!
!  Handle electronic stopping power according to Lindhard, interpolated with Bethe stopping 
!  at high energies
!
!  Used Modules:
!  ============
!
   use Const_mod,       only: INT4, REAL8,  &
                              ZERO, ONE, FIVE, HALF, BIG

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
   public :: initSeLind             ! initialize parameters
!            ==========
   public :: evalSeLind             ! evaluate the energy loss
!            ==========
!
!  Private Data:
!  ============
!
   logical :: init = .false.

   real(REAL8), allocatable :: faclin(:,:)   ! (NELEM,NELEM)
   real(REAL8), allocatable :: facbet(:,:)   ! (NELEM,NELEM)
   real(REAL8), allocatable :: ebethe(:,:)   ! (NELEM,NELEM)
   real(REAL8), allocatable :: elinmx(:,:)   ! (NELEM,NELEM)
!
!=======================================================================
!
contains
!
!=======================================================================
!
   subroutine initSeLind( )
!             ==========
!  
!     Initialize the parameters of Lindhard/Bethe stopping
!
!  Arguments:
!
!     None
!
!  Local Variables:
!
      integer(INT4) :: ielem1, ielem2
      real(REAL8) :: eexcit
!
!  Start of Executable Code:
!  ------------------------
!
      init = .true.
!
!  Allocate memory
!
      allocate( faclin( nelem, nelem ) )

      allocate( facbet( nelem, nelem ) )
      allocate( ebethe( nelem, nelem ) )

      allocate( elinmx( nelem, nelem ) )
!
!  Calculate parameters for the Lindhard stopping power
!
      do ielem1 = 1, nelem
         do ielem2 = 1, nelem
            faclin(ielem2,ielem1) = ( 1.212*real(elem(ielem1)%z)**(7./6.) * elem(ielem2)%z /  &
                 (real(elem(ielem1)%z)**(2./3.) + real(elem(ielem2)%z)**(2./3.))**1.5 /  &
                  sqrt(elem(ielem1)%m) )
            faclin(ielem2,ielem1) = faclin(ielem2,ielem1) * se%corlin(ielem2,ielem1)
         end do
      end do
!
!  Calculate parameters for the modified Bethe stopping power
!
      do ielem1 = 1, nelem
         do ielem2 = 1, nelem
            facbet(ielem2,ielem1) =  &
                    ( 2.375e6 * elem(ielem1)%z**2 * elem(ielem2)%z * elem(ielem1)%m )
!     1          ( 5210.34*elem(ielem1)%z**2*elem(ielem2)%z*elem(ielem1)%m )  

            if( elem(ielem2)%z < 13 ) then
               eexcit = 12.0 + 7.0 * real(elem(ielem2)%z)**(-ONE)
            else
               eexcit = 9.76 + 58.5 * real(elem(ielem2)%z)**(-1.19)
            end if
            ebethe(ielem2,ielem1) =  &
                    ( 455.7 * elem(ielem1)%m * elem(ielem2)%z * eexcit )
!                   ( 0.25 * elem(ielem1)%m * elem(ielem2)%z * eexcit)
         end do
      end do
!
!  Calculate the limit below which only Lindhard stopping has to be taken into
!  account
!
      do ielem1=1,nelem
         do ielem2 =1, nelem
            if( faclin(ielem2,ielem1) == ZERO ) then
               elinmx(ielem2,ielem1) = BIG
            else
               elinmx(ielem2,ielem1) = ( 0.01 * facbet(ielem2,ielem1) / faclin(ielem2,ielem1) ) &
                    ** ( ONE / (se%powlin(ielem1) + ONE) )
            end if
         end do
      end do

      return

   end subroutine initSeLind
!
!=======================================================================
!
   subroutine evalSeLind( e, ielem1, ielem2, se_val )
!             ==========  i    i        i      o
!
!     Evaluate energy loss according to Lindhard,
!     interpolated with Bethe stopping at high energies
!
!  Arguments:
!
      real(REAL8), intent(in) :: e              ! ion energy [eV]
      integer(INT4), intent(in) :: ielem1       ! index of the ion
      integer(INT4), intent(in) :: ielem2       ! index of the atom
      real(REAL8), intent(out) :: se_val        ! stopping power [eV A^2]
!
!  Local Variables:
!
      real(REAL8) :: seb, eb, seb1, se1
!
!  Start of Executable Code:
!  ------------------------
!
      if( .not.init )  goto 901

      if( faclin(ielem2,ielem1) == ZERO ) then      ! no electronic stopping if faclin=0
         se_val = ZERO
         return
      endif

      if( se%powlin(ielem1) == HALF ) then
         se_val = faclin(ielem2,ielem1) * sqrt(e)
      else
         se_val = faclin(ielem2,ielem1) * e**se%powlin(ielem1)
      endif

      if( e > elinmx(ielem2,ielem1) ) then
         eb = e/ebethe(ielem2,ielem1)
         seb = facbet(ielem2,ielem1)/e * &
               log( eb + se%c0bethe(ielem2,ielem1) + se%c1bethe(ielem2,ielem1)/eb )
         if( se%powint(ielem2,ielem1) == ONE ) then
            se_val = ONE/(ONE/se_val+ONE/seb)
         else
            se1 = (ONE/se_val)**se%powint(ielem2,ielem1)
            seb1 = (ONE/seb)**se%powint(ielem2,ielem1)
            se_val = (ONE/(se1+seb1))**(ONE/se%powint(ielem2,ielem1))
         endif
      endif

      return

 901  error stop 'evalSeLind: not initialized'

   end subroutine evalSeLind
!
!=======================================================================
!
end module Se_Lind_mod
