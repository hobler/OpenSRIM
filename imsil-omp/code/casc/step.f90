module Casc_Step_mod
!
!  Proceed one step in crystalline material
!
!  Used Modules:
!  ============
!
   use Const_mod,          only: INT4, REAL8, &
                                 ONE
   use Types_mod,          only: Proj_type

   use Param_Xtal_mod,     only: xtal
   use Param_Mat_mod,      only: irxtal
   use Param_Dam_mod,      only: dam
   use Param_Sn_mod,       only: sn
   use Param_Out_mod,      only: out

   use Casc_Selncp_mod,    only: selncp
   use Casc_Branch_mod,    only: branch
   use Sn_Scatter_mod,     only: scatter
   use Sn_FFPcyl_mod,      only: newFFPcyl
   use Xtal_OldSites_mod,  only: updateOldSites
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
   public :: step                            ! perform one step in crystalline material
!            ====
!
contains
!
!=======================================================================
!
   subroutine step( proj, ig )
!             ====   io   io
!
!     Proceed one step in any material, consisting of a free flight 
!     path and a collision.
!
!  Arguments:
!
      type(Proj_type), allocatable, intent(inout) :: proj(:)   ! projectile stack
      integer(INT4), intent(inout) :: ig                       ! index of recoil generation
!
!  Local Variables:
!
      logical :: lrecoil
      type(Proj_type) :: recoil
      real(REAL8) :: posc(3)

      real(REAL8) ::  ffp, p
      real(REAL8) :: pmaxse, dee, cospsi
!
!  Start of Executable Code:
!  ------------------------
!
!  Determine collision partner
!
      call selncp( proj(ig), posc, ffp, pmaxse, lrecoil, recoil, p )
!          ======     i       o     o     o        o       o     o
!
!  Treat scattering event (including electronic energy loss)
!
      call scatter( proj(ig), posc, ffp, pmaxse, dee, cospsi, lrecoil, recoil, p )
!          =======    io       i    io     i      o     o       io       io    i
!
!  Update old site list and ffp cylinder
!
      if( lrecoil ) then
         if( irxtal /= 0 )  call updateOldSites( ig, ffp, cospsi, recoil )
!                                ==============  i    i     i       i
         call newFFPcyl( proj(ig)%pos(:), ig, posc(:) )
!             =========        i          i     i
      else
         if( irxtal /= 0 )  call updateOldSites( ig, ffp )
!                                ==============  i    i
      end if
!
!  Damage statistics and branching of the trajectory
!
      if( lrecoil ) then
         call branch( proj, ig, dee, recoil, posc(:) )
!             ======   io   io   i     io      i
      else
         call branch( proj, ig, dee )
!             ======   io   io   i
      end if

   end subroutine step
!
!=======================================================================
!
end module Casc_Step_mod
