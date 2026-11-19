module Se_Loc_mod
!
!  Handle the impact parameter dependent part of electronic stopping
!
!  Used Modules:
!  ============
!
   use Const_mod,       only: INT4, REAL8,  &
                              ZERO, ONE, PI

   use Param_Setup_mod, only: nelem
   use Param_Elem_mod,  only: elem
   use Param_Se_mod,    only: se

   use Se_Firs_mod,     only: evalSeFirs, evalCorXnlFirs
   use Se_OR_mod,       only: evalSeOR, evalCorXnlOR
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
   public :: initSeLoc              ! initialize parameters for local electronic stopping
!            =========
   public :: setSeScreen            ! set screening length for impact parameter dependent
!            ===========            !    electronic stopping
   public :: evalSeLoc              ! evaluate local part of electronic energy loss
!            =========
   public :: evalCorXnl             ! evaluate xnl correction due to neglect of local part
!            ==========             !    of originating from p > pmax
!
!  Private Data:
!  ============
!
   real(REAL8), allocatable :: aor(:,:)         ! (nelem,nelem)
   real(REAL8), allocatable :: ascr_save(:,:)   ! (nelem,nelem)
   real(REAL8), allocatable :: ascr
!$omp threadprivate( ascr )
!
!=======================================================================
!
contains
!
!=======================================================================
!
   subroutine initSeLoc( )
!             =========
!
!     Initializes paramerters for local electronic stopping calculation
!
!     This subroutine must be called from a single thread
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
!  Allocate memory
!
      allocate( aor( nelem, nelem ) )
      allocate( ascr_save( nelem, nelem ) )
!
!  Calculate screening length
!
      do ielem1 = 1, nelem
         do ielem2 = 1, nelem
            aor(ielem2,ielem1) = 0.8853 * 0.529 /  &
               ( real(elem(ielem1)%z)**0.23 + real(elem(ielem2)%z)**0.23 ) / 0.3
            if( se%nfacscr(ielem2,ielem1) == 1 ) then
               ascr_save(ielem2,ielem1) = se%facscr(1,ielem2,ielem1) * aor(ielem2,ielem1)
            end if
         end do
      end do

   end subroutine initSeLoc
!
!=======================================================================
!
   subroutine setSeScreen( e, xnl, ielem1, ielem2 )
!             ===========  i   i     i       i
!
!     Set the screening length of local electronic stopping
!
!  Arguments:
!
      real(REAL8), intent(in) :: e              ! projectile energy (eV)
      real(REAL8), intent(in) :: xnl            ! nonlocal fraction
      integer(INT4), intent(in) :: ielem1       ! element index of projectile
      integer(INT4), intent(in) :: ielem2       ! element index of recoil
!
!  Local variables:
!
      integer(INT4) :: i
      real(REAL8) :: facscr
!
!  Start of Executable Code:
!  ------------------------
!
      if( se%nfacscr(ielem2,ielem1) == 1 ) then
         if( se%lfac2scr(ielem1) ) then
            facscr = se%facscr(1,ielem2,ielem1) +  &
                     ( se%fac2scr(ielem2,ielem1) - se%facscr(1,ielem2,ielem1) ) * xnl
         else
            ascr = ascr_save(ielem2,ielem1)
            return
         end if
      else
         if( e <= se%enl(1,ielem2,ielem1) ) then
            facscr = se%facscr(1,ielem2,ielem1)
         else if( e >= se%enl(se%nfacscr(ielem2,ielem1),ielem2,ielem1) ) then
            facscr = se%facscr(se%nfacscr(ielem2,ielem1),ielem2,ielem1)
         else
            do i = 2, se%nfacscr(ielem2,ielem1)-1
               if( e <= se%enl(i,ielem2,ielem1) )  exit
            end do
            facscr = ( se%facscr(i-1,ielem2,ielem1) * (se%enl(i,ielem2,ielem1) - e) +  &
                       se%facscr(i,ielem2,ielem1) * (e - se%enl(i-1,ielem2,ielem1)) ) /  &
                     ( se%enl(i,ielem2,ielem1) - se%enl(i-1,ielem2,ielem1) )
         end if
      end if

      ascr = facscr * aor(ielem2,ielem1)

   end subroutine setSeScreen
!
!=======================================================================
!
   subroutine evalSeLoc( fac, p, ielem1, ielem2 )
!             =========   o   i    i       i
!
!     Evaluate the local electronic energy loss factor.
!     The local electronic energy loss factor has to be multiplied with the stopping power to
!     obtain the local electronic energy loss.
!
!     NOTE: Before this subroutine is called, the screening length ascr has to be set using
!           the subroutine setSeScreen.
!
!  Arguments:
!
      real(REAL8), intent(out) :: fac           ! local electronic energy loss factor
      real(REAL8), intent(in) :: p              ! impact parameter
      integer(INT4), intent(in) :: ielem1       ! element index of projectile
      integer(INT4), intent(in) :: ielem2       ! element index of recoil
!
!  Local variables:
!
!     None
!
!  Start of Executable Code:
!  ------------------------
!
      if( se%locmodel(ielem2,ielem1) == 'or' ) then
         call evalSeOR( fac, p, ascr )
!             ========   o   i   i
      else if( se%locmodel(ielem2,ielem1) == 'firsov' ) then
         call evalSeFirs( fac, p, ascr, ielem1, ielem2 )
!             ==========   o   i   i      i       i
      else
         fac = ZERO
      end if

   end subroutine evalSeLoc
!
!=======================================================================
!
   subroutine evalCorXnl( corxnl, pmax, ielem1, ielem2 )
!             ==========    o      i      i       i
!
!     Evaluate the correction to the nonlocal electronic energy loss due to the maximum impact
!     parameter
!
!     NOTE: Before this subroutine is called, the screening length ascr has to be set using
!           the subroutine setSeScreen.
!
!  Arguments:
!
      real(REAL8), intent(out) :: corxnl        ! correction factor for nonlocal fraction
      real(REAL8), intent(in) :: pmax           ! maximum impact parameter used when
                                                !    calling the subroutine
      integer(INT4), intent(in) :: ielem1       ! element index of projectile
      integer(INT4), intent(in) :: ielem2       ! element index of recoil
!
!  Local variables:
!
!     None
!
!  Start of Executable Code:
!  ------------------------
!
      if( se%locmodel(ielem2,ielem1) == 'or' ) then
         call evalCorXnlOR( corxnl, pmax, ascr )
!             ============    o      i     i
      else if( se%locmodel(ielem2,ielem1) == 'firsov' ) then
         call evalCorXnlFirs( corxnl, pmax, ascr, ielem1, ielem2 )
!             ==============    o      i     i      i       i
      else
         corxnl = ZERO
      end if

   end subroutine evalCorXnl
!
!=======================================================================
!
end module Se_Loc_mod
