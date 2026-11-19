module Se_KKK_mod
!
!  Handles the electronic stopping power according to Konac, Klatt, and Kalbitzer 
!  (NIM B146, 106, 1998) for all ions in C (not implemented) and Si
!
!  Used Modules:
!  ============
!
   use Const_mod,       only: INT4, REAL8,  &
                              ZERO, TEN, EULER, BIG

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
   public :: initSeKKK              ! initialize parameters
!            =========
   public :: evalSeKKK              ! evaluate the energy loss
!            =========
!
!  Private Data:
!  ============
!
   logical :: init = .false.

   real(REAL8), allocatable :: s(:,:)        ! (NELEM,NELEM)
   real(REAL8), allocatable :: alpha(:,:,:)  ! (0:3,NELEM,NELEM)
   real(REAL8), allocatable :: beta(:,:)     ! (NELEM,NELEM)
!
!=======================================================================
!
contains
!
!=======================================================================
!
   subroutine initSeKKK ()
!             =========
!
!     Initializes the namelist variables corresponding to KKK stopping
!
!  Arguments:
!
!     None
!
!  Local Variables:
!
      integer(INT4) :: ielem1, ielem2
!
!  Start of Executable Code:
!  ------------------------
!
      init = .true.
!
!  Allocate memory
!
      allocate( s( nelem, nelem ) )
      allocate( alpha( 0:3, nelem, nelem ) )
      allocate( beta( nelem, nelem ) )
!
!  Calculate parameters
!
      do ielem1 = 1, nelem
         do ielem2 = 1, nelem
            if( elem(ielem2)%z == 14 ) then
               s(ielem2,ielem1) = (0.18/real(elem(ielem1)%z) + 0.51)
               alpha(0,ielem2,ielem1) = 0.0035*real(elem(ielem1)%z)**(-0.47)
               alpha(1,ielem2,ielem1) = (-0.0122)*real(elem(ielem1)%z)**(-1.63)
               alpha(2,ielem2,ielem1) = 0.0574*real(elem(ielem1)%z)**(-1.24)
               alpha(3,ielem2,ielem1) = 0.333*real(elem(ielem1)%z)**(-2.18)
               beta(ielem2,ielem1) = 15.7
            else
               if( se%semodel(ielem2,ielem1) == 'kkk' )  &
                     stop '&separ: KKK model only allowed for Si target atoms'
            end if
            if( alpha(0,ielem2,ielem1) == ZERO )  alpha(0,ielem2,ielem1) = BIG
         end do
      end do

      return

   end subroutine initSeKKK
!
!=======================================================================
!
   subroutine evalSeKKK( e, ielem1, ielem2, se )
!             =========  i    i       i     o
!
!     evaluate energy loss according to Konac, Klatt, and Kalbitzer
!
!  Arguments:
!
      real(REAL8), intent (in) :: e             ! ion energy [eV]
      integer(INT4), intent (in) :: ielem1      ! index of the ion
      integer(INT4), intent (in) :: ielem2      ! index of the atom
      real(REAL8), intent (out) :: se           ! stopping power [eV A^2]
!
!  Local Variables:
!
      real(REAL8) :: em, es
!
!  Start of Executable Code
!  ------------------------
!
      if( .not.init )  goto 901
!
!  energy in MeV/amu
! 
      em = e / (1e6*elem(ielem1)%m)

      es = em**s(ielem2,ielem1)

      se = es * log( EULER + beta(ielem2,ielem1)*em ) / &
           ( alpha(0,ielem2,ielem1) + alpha(1,ielem2,ielem1) * sqrt(em) + &
             (alpha(2,ielem2,ielem1) + alpha(3,ielem2,ielem1) * es) * em )
!
!  Convert from eVcm^2/1e15 to eV A^2
!
      se = se * TEN

      return

 901  error stop 'evalSeKKK: not initialized'

   end subroutine evalSeKKK
!
!=======================================================================
!
end module Se_KKK_mod
