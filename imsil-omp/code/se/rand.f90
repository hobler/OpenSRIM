module Se_Rand_mod
!
!  Handle random electronic stopping power.
!
!  Used Modules:
!  ============
!
   use Const_mod,    only: INT4, REAL8

   use Param_Se_mod, only: se

   use Se_Lind_mod,  only: initSeLind, evalSeLind
   use Se_ZBL_mod,   only: initSeZBL, evalSeZBL
   use Se_KKK_mod,   only: initSeKKK, evalSeKKK
   use Se_Table_mod, only: initSeTable, evalSeTable
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
   public :: initSeRand             ! initialize electronic stopping
!            ==========             !   parameters
   public :: evalSeRand             ! evaluate electronic stopping
!            ==========
!
!  Private Data:
!  ============
!
   logical :: init = .false.
!
!=======================================================================
!
contains
!
!=======================================================================
!
   subroutine initSeRand( )
!             ==========
!
!     Initialize random electronic energy loss
!
!  Arguments:
!
!     None
!
!  Start of Executable Code:
!  ------------------------
!
      init=.true.

      if( any( se%semodel(:,:) == 'lindhard' ) )  call initSeLind( )
!                                                      ==========
      if( any( se%semodel(:,:) == 'zbl' ) )  call initSeZBL( )
!                                                 =========
      if( any( se%semodel(:,:) == 'kkk' ) )  call initSeKKK( )
!                                                 =========
      if( any( se%semodel(:,:) == 'table' .or. se%semodel(:,:) == 'srim' ) )  &
         call initSeTable( )
!             ===========

   end subroutine initSeRand
!
!=======================================================================
!
   subroutine evalSeRand( e, ielem1, ielem2, se_val )
!             ==========  i     i        i    o
!
!     Evaluate the random stopping power by calling specific subroutines
!     according the semodel
!
!  Arguments:
!
      real(REAL8), intent(in) :: e              ! ion energy [eV]
      integer(INT4), intent(in) :: ielem1       ! chemical element index of ion
      integer(INT4), intent(in) :: ielem2       ! chemical element index of target atom
      real(REAL8), intent(out) :: se_val        ! electronic stopping power [eV A^2]
!
!  Start of Executable Code:
!  ------------------------
!
      if( .not.init )  goto 901

      select case( se%semodel(ielem2, ielem1) )

         case( 'lindhard' )                      ! Lindhard+modified Bethe stopping
            call evalSeLind( e, ielem1, ielem2, se_val )
!                ==========  i    i        i       o
         case( 'zbl' )                           ! ZBL85 stopping power
            call evalSeZBL( e, ielem1, ielem2, se_val )
!                =========  i    i        i     o
         case( 'kkk' )                           ! KKK model
            call evalSeKKK( e, ielem1, ielem2, se_val )
!                =========  i    i        i      o
         case( 'table', 'srim' )                 ! SE to be interpolated in table
            call evalSeTable( e, ielem1, ielem2, se_val )
!                ===========  i    i        i      o
         case( 'none')
         case default
            error stop 'evalSeRand: unknown semodel!'

      end select

      return

 901  error stop 'evalSeRand: not initialized'

   end subroutine evalSeRand
!
!=======================================================================
!
end module Se_rand_mod
