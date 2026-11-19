module Se_Firs_Save_mod
!
!  Handles inelastic energy loss according to Firsov (Sov. Phys. JETP 9 (1959) 1076)
!
!  Used Modules:
!  ============
!
   use Const_mod

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
   public :: initSeFirs             ! initializes the table of the loss function
!            ==========             !   
   public :: evalSeFirs             ! evaluates the energy loss according to Firsov
!            ==========             !   
!
!  Private Data:
!  ============
!
   logical :: init = .false.

   real(REAL8), allocatable :: facfir(:,:)   ! prefactor (sqrt(eV))
   real(REAL8), allocatable :: ascr(:,:)     ! screening length (A)
!
!=======================================================================
!
contains
!
!=======================================================================
!
   subroutine initSeFirs( )
!             ==========
!
!     Initializes the parameters of Firsov stopping
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
!  Factors in the energy loss formula
!
      allocate( facfir( nelem, nelem ) )
      allocate( ascr( nelem, nelem ) )

      do ielem1 = 1, nelem
         do ielem2 = 1, nelem
            if( .true. ) then                   ! Firsov
               facfir(ielem2,ielem1) = se%corfir(ielem2,ielem1) * 5.973e-2  &
                     * real( elem(ielem1)%z  + elem(ielem2)%z ) ** (FIVE*THIRD)  &
                     / sqrt(elem(ielem1)%m)
               ascr(ielem1,ielem2) = 3.29 / real( elem(ielem1)%z  + elem(ielem2)%z ) ** THIRD
            else                                ! Kishinevskii
               facfir(ielem2,ielem1) = se%corfir(ielem2,ielem1) * 4.17e-2  &
                     * max( elem(ielem1)%z, elem(ielem2)%z )  &
                     * ( sqrt( real(elem(ielem1)%z) ) + sqrt( real(elem(ielem2)%z) ) )  &
                     * ( elem(ielem1)%z ** SIXTH + elem(ielem2)%z ** SIXTH )  &
                     / sqrt(elem(ielem1)%m)
               ascr(ielem1,ielem2) = 0.701  &
                     * ( elem(ielem1)%z ** SIXTH + elem(ielem2)%z ** SIXTH )  &
                     / sqrt( real( max( elem(ielem1)%z, elem(ielem2)%z ) ) )
            end if
         enddo
      enddo

   end subroutine initSeFirs
!
!=======================================================================
!
   subroutine evalSeFirs( e, p, ielem1, ielem2, de )
!             ==========  i  i    i       i     o
!
!     Calculates Firsov's electronic energy loss
!
!  Arguments:
!
      real(REAL8), intent (in) :: e             ! energy of the ion [eV]
      real(REAL8), intent (in) :: p             ! impact parameter [A]
      integer(INT4), intent (in) :: ielem1      ! index of the ion
      integer(INT4), intent (in) :: ielem2      ! index of the target atom
      real(REAL8), intent (out) :: de           ! energy loss [eV]
!
!  Local Variables:
!
!     None
!
!  Start of Executable Code:
!  ------------------------
!
      if( .not.init )  goto 901

      if( .true. ) then                         ! Firsov
         de = facfir(ielem2,ielem1) / (ONE + p/ascr(ielem2,ielem1))**5 * sqrt(e)
      else                                      ! Kishinevskii
         de = facfir(ielem2,ielem1) / (ONE + p/ascr(ielem2,ielem1))**3 * sqrt(e)
      end if

      return

901   error stop 'evalSeFirs: not initialized'

   end subroutine evalSeFirs
!
!=======================================================================
!
end module Se_Firs_mod
