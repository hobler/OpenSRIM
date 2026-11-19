module Eval_mod
!
!  Evaluates the results of the cascade calculations
!
!  Used Modules:
!  ============
!
!$ use omp_lib,                  only: omp_get_thread_num

   use Const_mod,                only: INT4,REAL8,  &
                                       ZERO, ONE
   use Types_mod,                only: Proj_type
   use Undefine_mod,             only: undefine
   use Vector_mod,               only: abs, abs2, operator( .cross. )

   use Param_Setup_mod,          only: nelem, ldyn, ndim, nionbk, niongeom, nionhis, nionout
   use Param_Geom_mod,           only: geom
   use Param_Elem_mod,           only: elem
   use Param_Dam_mod,            only: dam
   use Param_Out_mod,            only: out

   use Geom_BndCond_mod,         only: applyBndCond3d

   use Target_Region_mod,        only: sputter
   use Target_Cell_mod,          only: scoreCell, writeCell, writeDamCell, mergeCell, relaxCell
   use Target_Surface2d_SDF_mod, only: writeSurface2d
   use Xtal_Ipos_mod,            only: getPosI, operator( == )
   use Eval_Casc_mod,            only: casc10, casc11, casc21, casc23, casc26
   use Eval_Flow_mod,            only: flow10, flow11, flow20, flow21, flow23, flow26
   use Eval_Rbs_mod,             only: rbs10, rbs22, rbs30, rbs40
   use Eval_Stat_mod,            only: stat10, stat11, stat12,  &
                                       stat22, stat2x, stat24, stat25,  &
                                       stat26, stat27, stat2e, stat30, stat40
   use Eval_Tra_mod,             only: tra10, tra12, tra2x, tra30
   use Eval_Dam_mod,             only: dam20, dam21, dam23
   use Eval_Det_mod,             only: det23
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
   public :: eval10                 ! Initialize at begin of simulation
!            ======                  
   public :: eval11                 ! Initialize at begin of molecule cascade
!            ======
   public :: eval12                 ! Initialize at beginning of projectile trajectory
!            ======
   public :: eval20                 ! Evaluate when no trajectory is followed
!            ======
   public :: eval21                 ! Evaluate at beginning of trajectory
!            ======
   public :: eval22                 ! Evaluate along trajectory
!            ======
   public :: eval23                 ! Evaluate at end of trajectory
!            ======
   public :: eval24                 ! Evaluate at end of trajectory of a backscattered or
!            ======                 !    transmitted projectile
   public :: eval26                 ! Actions at end of each cascade
!            ======
   public :: eval27                 ! Actions after a certain number of cascades
!            ======
   public :: eval2e                 ! Evaluate at end of trajectory if projectile is lost
!            ======
   public :: eval30                 ! Evaluate at end of simulation
!            ======
!
contains
!
!=======================================================================
!
   subroutine eval10( )
!             ======
!
!     Initialize variables at the beginning of the whole simulation
!
!  Arguments: 
!
!     None
!
!  Local Variables:
!
!     None
!
!  Start of Executable Code:
!  ------------------------
!
!  Additional recoil flow
!
      if( dam%lflow )  call flow10( )
!                           ======
!
!  Cascade analysis
!
      if( dam%lcasc )  call casc10( )
!                           ======
!
!  Statistics
!
      call stat10( )
!          ======
!
!  Trajectory recording
!
      if( any( out%ltra(:) ) )  call tra10( )
!                                    =====
!
!  RBS statistics
!
      if( out%lrbs )  call rbs10( )
!                          =====

   end subroutine eval10
!
!=======================================================================
!
   subroutine eval11( )
!             ======
!
!     Initialize variables at the beginning of the simulation of a recoil
!     cascade
!
!  Arguments: 
!
!     None
!
!  Local Variables:
!
!     None
!
!  Remark:
!
!     Stack should always be clear at the beginning of a recoil cascade, if
!     the previous cascade has been handled properly. But just to make sure
!     we clear it again.
!
!  Start of Executable Code:
!  ------------------------
!
!  Additional recoil flow
!
      if( dam%lflow )  call flow11( )
!                           ======
!
!  Cascade analysis
!
      if( dam%lcasc )  call casc11( )
!                           ======
!
!  Statistics
!
      call stat11( )
!          ======
!
!  Trajectory recording
!
!     None
!
   end subroutine eval11
!
!=======================================================================
!
   subroutine eval12( proj )
!             ======   i
!
!     Evaluations at the beginning of the simulation of a new ion atom
!
!  Arguments:
!
      type(Proj_type), intent(in) :: proj       ! projectile at start of trajectory
!
!  Local Variables:
!
!     None
!
!  Start of Executable Code:
!  ------------------------
!
!  Cascade analysis
!
!     None
!
!  Statistics
!
      call stat12( proj%pos )
!          ======    i
!
!  Trajectory recording
!
      if( any( out%ltra(:) ) )  call tra12( proj )
!                                    =====   i

   end subroutine eval12
!
!=======================================================================
!
   subroutine eval20( proj )
!             ======   i
!
!     Evaluate when no trajectory is followed
!
!  Arguments:
!
      type(Proj_type), intent(in) :: proj       ! projectile
!
!  Local Variables:
!
      logical :: allow_score
      real(REAL8) :: ddam, edam, x, y, z
!
!  Start of Executable Code:
!  ------------------------
!
!  Additional recoil flow
!
      if( dam%lflow ) then

         call flow20( proj )
!             ======   i
!
!  Cascade analysis
!
      else if( dam%lcasc .and. proj%origin%internal ) then

         call casc23( proj, sub=.true. )
!             ======   i        i

      else if( dam%ldet ) then

         continue                               ! do nothing
!
!  Statistics
!
      else

         call dam20( proj, ddam, edam )
!             =====   i     o     o

         call stat2x( proj, ddam, ddam, edam, 0 )
!             ======   i     i     i     i    i

         if( ldyn ) then
            x = proj%pos(1)
            y = proj%pos(2)
            z = proj%pos(3)
            call applyBndCond3d( x, y, z, .false., allow_score )
!                ============== io io io     i          o
            if( allow_score ) then           ! proj%id should already be the correct cell index
               call scoreCell( proj%ielem, proj%id, ZERO, ddam, ddam, edam )
!                   =========       i         i       i    i     i     i
            end if
         end if

      end if

   end subroutine eval20
!
!=======================================================================
!
   subroutine eval21( proj )
!             ======   i
!
!     Evaluate at the beginning of a trajectory
!
!  Arguments:
!
      type(Proj_type), intent(in) :: proj       ! vacancy
!
!  Local Variables:
!
      logical :: allow_score
      real(REAL8) :: w, ddami, ddamv, x, y, z
!
!  Start of Executable Code:
!  ------------------------
!
      if( dam%lflow ) then
!
!  Additional recoil flow
!
         call flow21( proj )
!             ======   i

      else if( dam%lcasc ) then
!
!  Cascade analysis
!
         call casc21( proj )                    ! store vacancy
!             ======   i

      elseif( dam%ldet ) then
!                                               ! score vacancy into simulation cell
         call det23( proj, -1 )
!             =====   i    i

      else                                      ! increment damage histograms
!
!  Damage and cell contents
!
         call dam21( proj, ddami, ddamv )
!             =====   i      o      o

         call stat2x( proj, ddami, ddamv, ZERO, 1 )
!             ======   i      i      i     i    i

         if( ldyn ) then
            if( elem(proj%ielem)%volatile ) then
               w = ZERO
            else
               w = -ONE
            end if

            x = proj%pos(1)
            y = proj%pos(2)
            z = proj%pos(3)
            call applyBndCond3d( x, y, z, .false., allow_score )
!                ============== io io io     i          o
            if( allow_score ) then           ! proj%id should already be the correct cell index
               call scoreCell( proj%ielem, proj%id, w, ddami, ddamv, ZERO )
!                   =========       i         i     i    i      i     i
            end if
         end if

      end if
!
   end subroutine eval21
!
!==========================================================================================
!
   subroutine eval22( proj, dee_, recoil_, posc_, lstart_ )
!             ======   i     i       i       i       i
!
!     Evaluate along the trajectory
!
!     All passed variables are those after ffp and collision(s)
!
!  Arguments:
!
      type(Proj_type), intent(in) :: proj             ! projectile
      real(REAL8), optional, intent(in) :: dee_       ! electronic energy loss
      type(Proj_type), optional, intent(in) :: recoil_! recoils
      real(REAL8), optional, intent(in) :: posc_(:)   ! collision point(s) (where impact
                                                      !    parameter is perpendicular to 
                                                      !    direction of motion)
      logical, optional, intent(in) :: lstart_        ! flag indicating recoil starts trajectory
!
!  Local Variables:
!
      logical :: lstart
!
!  Start of Executable Code:
!  ------------------------
!
!  Cascade analysis
!
!     None
!
!  Statistics (call it also when lcasc=t; statRBS !)
!
      call stat22( proj, dee_, recoil_ )
!          ======   i     i       i
!
!  Trajectory recording
!
      if( present( lstart_ ) ) then
         lstart = lstart_
      else
         lstart =.true.
      end if

      if( out%ltra(2) .and. (.not.out%shorttra .or. lstart) )  call tra2x( proj, 2 )
!                                                                    =====   i    i
!
!  RBS statistics
!
      if( present( recoil_ ) .and. proj%ielem == 1 .and. out%lrbs )  then
!
         call rbs22( proj, abs2( recoil_%pos(:)-posc_(:) ), recoil_%ielem )
!             =====    i   ====             i                     i
      end if

   end subroutine eval22
!
!=======================================================================
!
   subroutine eval23( proj )
!             ======   i
!
!     Evaluate at the end of a trajectory
!
!  Arguments:
!
      type(Proj_type), intent(in) :: proj    ! interstitial
!
!  Local Variables:
!
      logical :: allow_score
      real(REAL8) :: w, ddami, ddamv, edam, x, y, z
!
!  Start of Executable Code:
!  ------------------------
!
!  Additional recoil flow
!
      if( dam%lflow ) then

         call flow23( proj )
!             ======   i
!
!  Update cell contents of dynamic target cell structure
!
      else if( dam%lcasc .and. proj%origin%internal ) then
!
!  Cascade analysis
!
         call casc23( proj, sub=.false. )
!             ======   i         i

      else if( dam%ldet ) then
!
!  Score recoil into simulation cell
!
         call det23( proj, +1 )
!             =====   i    i

      else

         call dam23( proj, ddami, ddamv, edam )
!             =====   i      o      o     o

         call stat2x( proj, ddami, ddamv, edam, 3 )
!             ======    i     i      i     i    i

         if( ldyn ) then
            if( elem(proj%ielem)%volatile ) then
               w = ZERO
            else
               w = ONE
            end if
            x = proj%pos(1)
            y = proj%pos(2)
            z = proj%pos(3)
            call applyBndCond3d( x, y, z, .false., allow_score )
!                ============== io io io     i          o
            if( allow_score ) then           ! proj%id should already be the correct cell index
               call scoreCell( proj%ielem, proj%id, w, ddami, ddamv, proj%e )
!                   =========       i         i     i    i      i      i
            end if
         end if

      end if

   end subroutine eval23
!
!=======================================================================
!
   subroutine eval24( proj )
!             ======   i
!
!     Evaluate at the end of the trajectory of a backscattered or transmitted ion
!
!  Arguments:
!
      type(Proj_type), intent(in) :: proj       ! projectile
!
!  Local Variables:
!
!     None
!
!  Start of Executable Code:
!  ------------------------
!
      if( proj%dir(3) <= ZERO ) then            !    ... backscattered
!
!  Cascade analysis
!
!     None
!
!  Statistics
!
         call stat24( proj )
!             ======   i
!
!  Trajectory recording
!
         if( out%ltra(4) )  call tra2x( proj, 4 )
!                                =====   i    i

      else                                      !    ... transmitted
!
!  Cascade analysis
!
!     None
!
!  Statistics
!
         call stat25( proj )
!             ======   i
!
!  Trajectory recording
!
         if( out%ltra(5) )  call tra2x( proj, 5 )
!                                =====   i    i

      end if

   end subroutine eval24
!
!=======================================================================
!
   subroutine eval26( )
!             ======
!
!     Perform actions after each ion molecule has been implanted, working on the
!     own image only
!
!  Arguments: 
!
!     None
!
!  Local Variables:
!
!     None
!
!  Start of Executable Code:
!  ------------------------
!
!  Cascade analysis
!
      if( dam%lcasc )  call casc26( )
!                           ======
!
!  Statistics
!
      call stat26( )
!          ======

   end subroutine eval26
!
!=======================================================================
!
   subroutine eval27( iion, final_ )
!             ======   i      i
!
!     Perform actions after a certain number of ion chunks, working with or updating
!     global quantities.
!
!     The final evaluation (final argument present) serves to gather remaining data between
!     regular updates and the end of the simulation. For the final evaluation iion must
!     equal nion.
!
!     OpneMP: This subroutine must be called from a parallel region after a barrier.
!
!  Arguments:
!
      integer(INT4), intent(in) :: iion         ! index of ion
      logical, intent(in), optional :: final_   ! indicates final call of stat27
!
!  Local Variables:
!
      logical :: final, lbk, lmerge, lrelax
      integer(INT4) :: ithread
!
!  Start of Executable Code:
!  ------------------------
!
      if( present( final_ ) ) then
         final = final_
      else
         final = .false.
      end if
!
!  Sputtering of the target (requires passing of ielem to subroutine)
!
!      call sputter( ielem )
!           =======    i
!
!  Additional recoil flow
!
!      if( dam%lflow )  call flow26( )
!                           ======
!
!  Statistics
!
      call stat27( iion, final_ )
!          ======   i      i
!
!  Trajectory recording
!
!     None
!
!  Make a backup after nionbk ions
!
      ithread = 1
!$    ithread = 1 + omp_get_thread_num( )
      if( nionbk == 0 .or. ithread > 1 ) then
         lbk = .false.
      else
         lbk = ( mod( iion, nionbk ) == 0 )
         if( final )  lbk = .not.lbk
      end if

      if( lbk ) then                            ! will not be entered for ithread > 1
         call stat40( iion )
!             ======    i
         if( out%lrbs )  call rbs40( )
!                             =====
      end if
!
!  Cell relaxation after each niongeom molecules
!
      if( ldyn ) then

         if( nionhis <= 0 ) then
            lmerge = .false.
         else
            lmerge = ( mod( iion, nionhis ) == 0 )
            if( final )  lmerge = .not.lmerge
         end if

         if( lmerge )  call mergeCell( )
!                           =========

         if( niongeom == 0 ) then               ! we assume mod(niongeom,nionhis) == 0
            lrelax = .false.                    ! so there is also a merge before a relax
         else
            lrelax = ( mod( iion, niongeom ) == 0 )
            if( final )  lrelax = .not.lrelax
         end if

         if( lrelax ) then
!$omp barrier
            call relaxCell( iion )
!                =========   i
!$omp master
            if( mod( iion, nionout ) == 0 ) then
               if( out%lcell )  call writeCell( iion )
!                                    =========    i
               if( ndim == 2 .and. out%lsurf )  call writeSurface2d( iion )
!                                                    ==============   i
            end if
!$omp end master
!$omp barrier
         end if

      end if

   end subroutine eval27
!
!=======================================================================
!
   subroutine eval2e( proj )
!             ======   i
!
!     Evaluate at the end of the trajectory if ion is lost
!
!  Arguments:
!
      type(Proj_type), intent(in) :: proj       ! projectile
!
!  Local Variables:
!
!     None
!
!  Start of Executable Code:
!  ------------------------
!
!  Cascade analysis
!
!     None
!
!  Statistics
!
      call stat2e( proj )
!          ======   i
!
!  Trajectory recording
!
!     None
!
   end subroutine eval2e
!
!=======================================================================
!
   subroutine eval30( iion )
!             ======   i
!
!     Evaluate at the end of simulation
!
!  Arguments:
!
      integer(INT4), intent(in) :: iion         ! index of ion
!
!  Local Variables:
!
!     None
!
!  Start of Executable Code:
!  ------------------------
!
!  Dynamic cell structure
!
!$omp master
      if( ldyn ) then
         if( mod( iion, nionout ) /= 0 ) then
            if( out%lcell )  call writeCell( iion )
!                                 =========    i
            if( ndim == 2 .and. out%lsurf )  call writeSurface2d( iion )
!                                                 ==============   i
         end if
         if( out%lhis )  call writeDamCell( )
!                             ============
      end if
!$omp end master
!
!  Cascade analysis
!
!     None
!
!  Statistics (currently only correctly calculated at end of simulation)
!
      call stat30( )
!          ======
!
!  Trajectory recording
!
      if( any( out%ltra(:) ) )  call tra30( )
!                                    =====
!
!  RBS statistics
!
      if( out%lrbs )  call rbs30( )
!                          =====
!
   end subroutine eval30
!
!=======================================================================
!
end module Eval_mod
