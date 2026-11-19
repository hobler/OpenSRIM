module Se_Straggle_mod
!
!  Calculate the electronic energy loss straggling figure
!
!  Used Modules:
!  ============
!
   use Const_mod,          only: INT4, REAL8,  &
                                 ZERO, ONE

   use Param_Setup_mod,    only: nelem
   use Param_Elem_mod,     only: elem
   use Param_Se_mod,       only: se
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
   public :: initWe                          ! initialize parameters
!            ======
   public :: evalWe                          ! calculate straggling figure
!            ======
!
!  Private Data:
!  ============
!
   logical :: init = .false.
   real(REAL8), allocatable :: webohr(:,:)   ! (NELEM,NELEM) Bohr straggling
   logical, allocatable :: lstrag(:,:)       ! (NELEM,NELEM) straggling flag
   real(REAL8), allocatable :: estrag(:,:)   ! (NELEM,NELEM) transition energy
!
!=======================================================================
!
contains
!
!=======================================================================
!
   subroutine initWe( )
!             ======
!
!     Initialize the electronic stopping parameters
!
!     The Bohr straggling formula has been checked on 2018-05-25 using the formula given
!     by E. Bonderup and P. Hvelplund (PRA 4, 562, 1971).
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

      allocate( lstrag( nelem, nelem ) )
      allocate( webohr( nelem, nelem ) )
      allocate( estrag( nelem, nelem ) )

      lstrag(:,:) = .false.
      webohr(:,:) = ZERO
      estrag(:,:) = ZERO

      do ielem1 = 1, nelem
         do ielem2 = 1, nelem
            lstrag(ielem2,ielem1) = se%lstrag(ielem2,ielem1)
            if( se%lstrag(ielem2,ielem1) ) then
               webohr(ielem2,ielem1) = 2599 * elem(ielem1)%z**2 * elem(ielem2)%z
               estrag(ielem2,ielem1) = se%estrag(ielem2,ielem1)
            end if
         end do
      end do

   end subroutine initWe
!
!=======================================================================
!
   subroutine evalWe( we, e, ielem1, ielem2 )
!             ======  o   i     i       i
!
!     Evaluate electronic energy loss straggling
!
!     The standard deviation of the electronic energy loss is given by sqrt(we*dens*dx)
!     where we is the electronic straggling figure computed here, dens is the atomic
!     density, and dx is the flight path.
!
!     The model is according to Konac et al. (NIM-B 136-138, 159, 1998). estrag can be used
!     as a fit parameter.
!
!     An alternative model with parameters for all (?) Z1, Z2 has been proposed by Q. Yang
!     et al. (NIM-B 61, 149, 1991).
!
!  Arguments:
!
      real(REAL8), intent(out) :: we            ! electronic straggling figure [eV^2 A^2]
      real(REAL8), intent(in) :: e              ! ion energy [eV]
      integer(INT4), intent(in) :: ielem1       ! chemical element index of the ion
      integer(INT4), intent(in) :: ielem2       ! index of the region
!
!  Local Variables:
!
!     None
!
!  Start of Executable Code:
!  ------------------------
!
      if( .not.init )  error stop 'evalWe: not initialized'

      if( lstrag(ielem2,ielem1) ) then
         we = webohr(ielem2,ielem1) * ( ONE - exp( -e/estrag(ielem2,ielem1) ) )
      else
         we = ZERO
      end if

   end subroutine evalWe
!
!=======================================================================
!
end module Se_Straggle_mod
