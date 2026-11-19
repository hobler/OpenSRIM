module Se_Firs_mod
!
!  Handle the impact parameter dependent part of electronic stopping according to Firsov
!
!  Used Modules:
!  ============
!
   use Const_mod,       only: INT4, REAL8,  &
                              ZERO, ONE, TWO, HALF, PI

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
   public :: evalSeFirs             ! evaluate local part of electronic energy loss according
!            ==========             !    to Oen and Robinson
   public :: evalCorXnlFirs         ! evaluate xnl correction due to neglect of local part
!            ==============         !    of originating from p > pmax
!
!  Private Procedures:
!  ==================
!
   private:: fun                    ! helper function for implementation of cutoff impact
!            ===                    !    parameter
!
!  Private Data:
!  ============
!
   logical :: usepmaxf = .true.     ! flag activating cutoff impact parameter
   real(REAL8) :: pmaxf = 3.38      ! cutoff impact parameter of electronic stopping
                                    !    (should be made an input parameter if successful)
!
!=======================================================================
!
contains
!
!=======================================================================
!
   subroutine fun( f, p, c, pmax )
!             ===  o  i  i    i
!
!     Helper function to evaluate the local electronic energy loss factor according to Firsov
!     with a cutoff.
!
!  Arguments:
!
      real(REAL8), intent(out) :: f             ! function value
      real(REAL8), intent(in) :: p              ! impact parameter / screening length
      real(REAL8), intent(in) :: c              ! power
      real(REAL8), intent(in) :: pmax           ! cutoff impact parameter / screening length
!
!  Local variables:
!
!     None
!
!  Start of Executable Code:
!  ------------------------
!
      f = (ONE+(c-ONE)*p) * (ONE+p)**(-c+ONE)
!      if( p == pmax ) then
!         f = (ONE+c*p*(ONE+HALF*(c-ONE)*p)) * (ONE+p)**(-c)
!      else
!         f = (ONE+(c-ONE)*p) * (ONE+p)**(-c+ONE) + HALF*(c-ONE)*(c-TWO)*p**2 * (ONE+pmax)**(-c)
!      end if

   end subroutine fun
!
!=======================================================================
!
   subroutine evalSeFirs( fac, p, ascr, ielem1, ielem2 )
!             ==========   o   i   i      i       i
!
!     Evaluate the local electronic energy loss factor according to Firsov.
!     The local electronic energy loss factor has to be multiplied with the stopping power to
!     obtain the local electronic energy loss.
!
!  Arguments:
!
      real(REAL8), intent(out) :: fac           ! local electronic energy loss factor
      real(REAL8), intent(in) :: p              ! impact parameter
      real(REAL8), intent(in) :: ascr           ! screening length
      integer(INT4), intent(in) :: ielem1       ! element index of projectile
      integer(INT4), intent(in) :: ielem2       ! element index of recoil
!
!  Local variables:
!
      real(REAL8) :: c, ff, k
!
!  Start of Executable Code:
!  ------------------------
!
      c = se%powfirs(ielem2,ielem1)

      if( usepmaxf ) then
         if( p >= pmaxf ) then
            fac = ZERO
         else
            call fun( ff, pmaxf/ascr, c, pmaxf/ascr )
!                ===  o       i       i      i
            k = (c-ONE)*(c-TWO) / (TWO*PI*ascr**2 * (ONE-ff))
            fac = k * (ONE+p/ascr)**(-c)
!            fac = k * ((ONE+p/ascr)**(-c) - (ONE+pmaxf/ascr)**(-c))
         end if
      else
         k = (c-ONE)*(c-TWO) / (TWO*PI*ascr**2)
         fac = k * (ONE+p/ascr)**(-c)
      end if

   end subroutine evalSeFirs
!
!=======================================================================
!
   subroutine evalCorXnlFirs( corxnl, pmax, ascr, ielem1, ielem2 )
!             ==============    o      i     i      i       i
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
      integer(INT4), intent(in) :: ielem1       ! element index of projectile
      integer(INT4), intent(in) :: ielem2       ! element index of recoil
!
!  Local variables:
!
      real(REAL8) :: pm, c, fc, ff
!
!  Start of Executable Code:
!  ------------------------
!
      c = se%powfirs(ielem2,ielem1)

      if( usepmaxf ) then
         if( pmax >= pmaxf ) then
            corxnl = ZERO
         else
            call fun( fc, pmax/ascr, c, pmaxf/ascr )
!                ===  o       i      i      i
            call fun( ff, pmaxf/ascr, c, pmaxf/ascr )
!                ===  o       fi      i      i
            corxnl = (fc-ff) / (ONE-ff)
         end if
      else
         pm = pmax / ascr
         corxnl = (ONE+(c-ONE)*pm) * (1+pm)**(-c+ONE)
      end if

   end subroutine evalCorXnlFirs
!
!=======================================================================
!
end module Se_Firs_mod
