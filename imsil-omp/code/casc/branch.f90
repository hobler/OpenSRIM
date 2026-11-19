module Casc_Branch_mod
!
!  Manages the branching of the recoil cascade and termination of trajectories. Also calls
!  all evaluations during trajectory simulations and treats free flights through the
!  vacuum if the ion reenters the target.
!
!  Note: this module is complicated by the necessity of virtual projectiles for accurate
!        sputtering simulations
!
!  Used Modules:
!  ============
!
   use Const_mod,          only: INT4, REAL8,  &
                                 ZERO, ONE, SQRT2,  &
                                 STATE
   use Array_mod,          only: changeShape

   use Types_mod,          only: Proj_type, Origin_type
   use String_mod,         only: str

   use Param_Setup_mod,    only: ldyn, nelemIon
   use Param_Mat_mod,      only: irxtal
   use Param_Dam_mod,      only: dam
   use Param_Out_mod,      only: out

   use Eval_mod,           only: eval20, eval21, eval23, eval24, eval2e
   use Eval_Tra_mod,       only: tra2x
   use Xtal_OldSites_mod,  only: copyOldSites, updateOldSites
   use Sn_FFPcyl_mod,      only: newFFPcyl
   use Sn_RestFFP_mod,     only: resetRestFFP

   use Casc_Check_mod,     only: chkRepl, chkLeave, chkStart, chkStop
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
   public :: branch                 ! check recoil cascade on branching, manages
!            ======                 !    the stack and considers damage generation
!
!  Private Procedures:
!  ==================
!
   private :: treatLeave            ! treat leaving or reentering projectile
!             ==========
   private :: treatStop             ! treat stopped projectile
!             =========
   private :: treatSaved            ! treat saved projectile states
!             ==========
   private :: treatRecoil           ! treat the reacoil
!             ===========
!
contains
!
!=======================================================================
!
   recursive subroutine branch( proj, ig, dee, recoil_, posc_ )
!                       ======   io   io   i     io      i
!
!     Check whether the recoil cascade branches or a trajectory terminates, adapt the
!     projectile stack, and take care of evaluations
!
!     Remarks
!
!     1) ig will be zero after exit of the routine if the last
!        projectile of a cascade has stopped
!
!     2) some of the components of recoil may be changed
!
!     3) it is assumed that recoil_%virtual == proj(ig)%virtual on entry
!
!  Arguments:
!
      type(Proj_type), allocatable, intent(inout) :: proj(:)   ! projectile stack
      integer(INT4), intent(inout) :: ig                       ! projectile index
      real(REAL8), intent(in) :: dee                           ! electronic energy loss on the
                                                               !    preceding free flight path
      type(Proj_type), optional, intent(inout) :: recoil_       ! recoils (undefined on exit)
      real(REAL8), optional, intent(in) :: posc_(3)            ! collision point (where impact
                                                               !    parameter is perpendicular
!                                                              !    to direction of motion)
!  Local Variables:
!
      logical :: lstart, virtual, turn_virtual, stopped
!
!  Start of Executable Code:
!  ------------------------
      !print *, 'e, z, lrecoil=', proj(ig)%e, Proj(ig)%pos(3), present(recoil_)
!
!  Check for replacement collisions and whether a recoil is to be followed
!
      if( present( recoil_ ) ) then

         if( dam%lrcoil .and.  &
             ( ( dam%lrcoil1 .and. .not.proj(ig)%reentered ) .or.  &
               ( dam%lrcoil2 .and. proj(ig)%reentered ) ) ) then

            call chkRepl( proj(ig), recoil_ )
!                =======    io        io

            call chkStart( recoil_, lstart, virtual )
!                ========     i       o        o

            recoil_%virtual = ( recoil_%virtual .or. virtual )
                                                ! avoid recoil_%virtual = .true.->.false.
         else

            lstart = .false.
            recoil_%virtual = .false.

         end if

      else

         lstart = .false.

      end if
!
!  Check if projectile leaves/reenters the target (result is in proj(ig)%state)
!  This also writes statistics along the trajectory (iflag=2)
!  If projectile is leaving and/or reentering, the projectile position, direction, and energy
!  are changed
!
      call chkLeave( proj(ig), dee, recoil_, posc_(:), lstart )
!          ========     io      i     i        i        i
!
!  Treat leaving/reentering projectiles
!
      if( proj(ig)%state == STATE%ERROR ) then

         call eval2e( proj(ig) )
!             ======     i

         proj(ig)%state = STATE%LEAVING         ! later treat as if leaving

      else

         if( proj(ig)%state == STATE%LEAVING .or. proj(ig)%state == STATE%REENTERED ) then

            call treatLeave( proj, ig )
!                ==========   io   i

         end if

         if( proj(ig)%state == STATE%MOVING ) then
!
!  Check if projectile stops or turns virtual
!
            call chkStop( proj(ig), stopped, virtual )
!                =======     i         o        o
            turn_virtual = virtual .and. .not. proj(ig)%virtual
!
!  Treat stopped/virtual projectiles (ig is incremented if turn_virtual=.true.)
!
            call treatStop( proj, ig, stopped, turn_virtual )
!                =========   io   io     i          i

         end if

      end if
!
!  For stopped or backscattered/sputtered projectile ...
!
      if( proj(ig)%state == STATE%STOPPED .or. proj(ig)%state == STATE%LEAVING ) then
!
!  ... delete projectile
!
         ig = ig - 1
!
!  ... treat saved states down the projectile stack, unless there is a recoil to be followed.
!
         if( .not.lstart ) then

            call treatSaved( proj(:), ig )
!                ==========    io     io

         end if

      end if

!
!  Treat the recoil
!
      if( present( recoil_ ) ) then

         call treatRecoil( proj, ig, recoil_, lstart )
!             ===========   io   io     i       i

      end if

      !do i=1,ig
      !   write( 998, * )  '***** ig=', str(i), ', proj=', proj(i)
      !end do
      !write( 998, * )  ig, proj(ig)%state, proj(ig)%virtual

   end subroutine branch
!
!=======================================================================
!
   subroutine treatLeave( proj, ig )
!             ==========   io   i
!
!     Adapt projectile stack due to a leaving or reentering projectile
!
!  Arguments:
!
      type(Proj_type), allocatable, intent(inout) :: proj(:)   ! projectile stack
      integer(INT4), intent(in) :: ig                          ! projectile index
!
!  Local Variables:
!
!     None
!
!  Start of Executable Code:
!  ------------------------
!                                               ! Note: virtual projectile always has a parent
      if( proj(ig)%virtual ) then               !       with a saved state

         proj(ig-1)%ededuct = proj(ig-1)%ededuct + proj(ig)%e

         select case( proj(ig-1)%state )
            case( STATE%SUBTHRESHOLD )
               proj(ig-1)%state = STATE%STARTING
            case( STATE%STOPPED )
               if( ldyn .or. proj(ig-1)%ielem > nelemIon )  &
                     proj(ig-1)%state = STATE%SUBTHRESHOLD
            case default
               error stop 'treatLeave: invalid parent projectile state'
         end select

         proj(ig)%virtual = .false.             ! turn projectile real (Note: parent's
                                                !    virtuality remains untouched)
      end if

      select case( proj(ig)%state )
         case( STATE%LEAVING )
            !if( proj(ig)%ielem == 1 .and. .not.proj(ig)%origin%internal ) then
            !   write( 998, * ) 'backscattered'
            !   write( 998, * ) proj(ig)
            !end if
            call eval24( proj(ig) )
!                ======     i
         case( STATE%REENTERED )
            proj(ig)%state = STATE%MOVING
         case default
            error stop 'treatLeave: invalid projectile state'
      end select

   end subroutine treatLeave
!
!=======================================================================
!
   subroutine treatStop( proj, ig, stopped, turn_virtual )
!             =========   io   io     i          i
!
!     Adapt projectile stack due to a stopped projectile or a projectile turning virtual
!
!  Arguments:
!
      type(Proj_type), allocatable, intent(inout) :: proj(:)   ! projectile stack
      integer(INT4), intent(inout) :: ig                       ! projectile index
      logical, intent(in) :: stopped                           ! flag indicating projectile
                                                               !    stops
      logical, intent(in) :: turn_virtual                      ! flag indicating projectile
                                                               !    turns virtual
!
!  Local Variables:
!
!     None
!
!  Start of Executable Code:
!  ------------------------
!
      if( stopped ) then                        ! stopped projectile

         proj(ig)%state = STATE%STOPPED

         if( proj(ig)%virtual ) then
            proj(ig-1)%ededuct = proj(ig-1)%ededuct + proj(ig)%ededuct
         else
            !if( proj(ig)%ielem == 1 .and. .not.proj(ig)%origin%internal ) then
            !   write( 998, * ) 'real stopped'
            !   write( 998, * ) proj(ig)
            !end if
            call eval23( proj(ig) )
!                ======     i
         end if

         if( out%ltra(3) )  call tra2x( proj(ig), 3 )
!                                =====     i      i

      else if( turn_virtual ) then              ! projectile turns virtual

         proj(ig)%state = STATE%STOPPED         ! save state

         ig = ig + 1
         if( ig > size( proj ) )  call changeShape( proj, ig, 'proj', 'treatStop' )
!                                      ===========   io   i     i          i
         proj(ig) = proj(ig-1)
         proj(ig)%ig = ig
         proj(ig)%state = STATE%MOVING
         proj(ig)%virtual = .true.
         proj(ig)%virtualparent = .false.
         proj(ig)%ededuct = ZERO

         if( irxtal /= 0 )  call copyOldSites( ig-1, ig )
!                                =============  i    i

         call newFFPcyl( proj(ig)%pos(:), ig )  ! initialize the first ffp cylinder
!             =========        i          i

         call resetRestFFP( ig )                ! reset old ffp list
!             ============  i

      end if

   end subroutine treatStop
!
!=======================================================================
!
   subroutine treatSaved( proj, ig )
!             ==========   io   io
!
!     Adapt projectile stack due to a stopped or leaving projectile, treating saved states
!     Note that there are only virtual saved states
!
!  Arguments:
!
      type(Proj_type), intent(inout) :: proj(:) ! projectile stack
      integer(INT4), intent(inout) :: ig        ! projectile index
!
!  Local Variables:
!
      logical :: saved, savedparent
!
!  Start of Executable Code:
!  ------------------------
!
      if( ig <= 0 )  return

      saved = ( proj(ig)%state == STATE%SUBTHRESHOLD .or.  &
                proj(ig)%state == STATE%STARTING .or.  &
                proj(ig)%state == STATE%STOPPED )

      do while( saved )

         if( ig == 1 ) then
            savedparent = .false.
         else
            savedparent = ( proj(ig-1)%state == STATE%SUBTHRESHOLD .or.  &
                            proj(ig-1)%state == STATE%STARTING .or.  &
                            proj(ig-1)%state == STATE%STOPPED )
         end if

         if( savedparent ) then

            proj(ig-1)%ededuct = proj(ig-1)%ededuct + proj(ig)%ededuct

            if( proj(ig)%state == STATE%STARTING )  call eval21( proj(ig) )
!                                                        ======      i
         else

            !if( proj(ig)%ededuct > ZERO )  print *, 'ededuct=', proj(ig)%ededuct
            proj(ig)%e = proj(ig)%e - proj(ig)%ededuct

            select case( proj(ig)%state )
               case( STATE%SUBTHRESHOLD )
                  if( .not. proj(ig)%virtualparent )  call eval20( proj(ig) )
!                                                          ======     i
                  !write( 998, * ) 'saved state subthreshold'
               case( STATE%STARTING )
                  call eval21( proj(ig) )
!                      ======      i
                  !write( 998, * ) 'saved state starting'
               case( STATE%STOPPED )
                  !if( proj(ig)%ielem == 1 .and. .not.proj(ig)%origin%internal ) then
                  !write( 998, * ) 'saved state stopped'
                  !   write( 998, * ) proj(ig)
                  !end if
                  call eval23( proj(ig) )
!                      ======      i
            end select

         end if

         ig = ig - 1

         saved = savedparent

      end do

   end subroutine treatSaved
!
!=======================================================================
!
   subroutine treatRecoil( proj, ig, recoil, lstart )
!             ===========   io   io    i       i
!
!     Adapt projectile stack due to a new recoil
!
!  Arguments:
!
      type(Proj_type), allocatable, intent(inout) :: proj(:)   ! projectile stack
      integer(INT4), intent(inout) :: ig                       ! projectile index
      type(Proj_type), intent(in) :: recoil                    ! new recoil
      logical, intent(in) :: lstart                            ! flag indicating recoil starts
                                                               !    trajectory
!
!  Local Variables:
!
!     None
!
!  Start of Executable Code:
!  ------------------------
!
      if( lstart ) then                         ! A new recoil is generated

         ig = ig + 1
         if( ig > size( proj ) )  call changeShape( proj, ig, 'proj', 'treatRecoil' )
!                                      ===========   io   i     i           i
         proj(ig) = recoil                      ! note: no binding energy considered
         proj(ig)%ig = ig

         if( proj(ig)%virtual ) then            ! save state of virtual recoil

            proj(ig)%state = STATE%SUBTHRESHOLD

            ig = ig + 1
            if( ig > size( proj ) )  call changeShape( proj, ig, 'proj', 'treatRecoil' )
!                                         ===========   io   i     i           i
            proj(ig) = proj(ig-1)
            proj(ig)%ig = ig
            proj(ig)%virtualparent = .true.
            proj(ig)%ededuct = ZERO

         else                                   ! evaluate real recoil immediately

            call eval21( proj(ig) )
!                ======     i

         end if

         if( out%ltra(-1) )  call tra2x( proj(ig), -1 )
!                                 =====     i       i
         if( out%ltra(1) )  call tra2x( proj(ig), 1 )
!                                =====     i      i

         proj(ig)%state = STATE%MOVING

                                                ! reset old lattice site list and store recoil's
                                                !    lattice site in old lattice sites list
         if( irxtal /= 0 )  call updateOldSites( ig, proj_=proj(ig) )
!                                ==============  i         i

         call newFFPcyl( proj(ig)%pos(:), ig )  ! initialize the first ffp cylinder
!             =========        i          i

         call resetRestFFP( ig )                ! reset old ffp list
!             ============  i

         call branch( proj, ig, ZERO )          ! results in ig=ig-1 if trajectory terminates
!             ======   io   io   i

      else                                      ! No recoil (subthreshold collision)

         if( .not.recoil%virtualparent )  call eval20( recoil )
!                                              ======    i

         if( out%ltra(0) )  call tra2x( recoil, 0)
!                                =====    i     i

      end if

   end subroutine treatRecoil
!
!=======================================================================
!
end module Casc_Branch_mod
