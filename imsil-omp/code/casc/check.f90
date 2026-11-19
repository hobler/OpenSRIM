module Casc_Check_mod
!
!  Provides subroutines that check whether a projectile stops or can be started
!
!  Used Modules:
!  ============
!
   use Const_mod,          only: INT4, REAL8,  &
                                 ZERO, ONE, HALF, SQRT2, PI,  &
                                 STATE
   use Vector_mod,         only: operator(.dot.), abs
   use Types_mod,          only: Proj_type, Origin_type

   use Param_Setup_mod,    only: nelem, lsput
   use Param_Geom_mod,     only: geom
   use Param_Elem_mod,     only: elem
   use Param_Ion_mod,      only: ion
   use Param_Xtal_mod,     only: xtal
   use Param_Mat_mod,      only: coord_xtal, irxtal
   use Param_Dam_mod,      only: dam
   use Param_Sn_mod,       only: sn
   use Param_Out_mod,      only: out

   use Casc_Surface_mod,   only: enterTarget, exitTarget
   use Eval_mod,           only: eval22
   use Geom_Coordsys_mod,  only: fromCoordSys
   use Target_mod,         only: getDomainProperty
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
   public :: chkRepl                ! check if replacement collision occurs
!            =======
   public :: chkLeave               ! check if projectile leaves the target
!            ========
   public :: chkStart               ! check whether recoil can start
!            ========
   public :: chkStop                ! check whether projectile stops motion
!            =======
!
!  Private Procedures:
!  ==================
!
   private :: chkKP                 ! check whether the Kinchin-Pease model can be
!             =====                 !    applied
!
contains
!
!=======================================================================
!
   subroutine chkRepl( proj, recoil )
!             =======   io     io
!
!     Check whether a replacement collision occurs and exchange projectile and
!     recoil if necessary. The exchange makes possible that the projectile (then
!     recoil) does not initiate a trajectory.
!
!  Arguments:
!
      type(Proj_type), intent(inout) :: proj    ! projectile
      type(Proj_type), intent(inout) :: recoil  ! recoil
!
!  Local Variables:
!
      logical :: lstart, stopped, virtual, virtualp, virtualr, virtualparentp, virtualparentr
      integer(INT4) :: igp, igr, statep, stater
      real(REAL8) :: posoldp(3), posoldr(3)
      type(Proj_type) :: proj_tmp
      type(Origin_type) :: originp, originr
!
!  Start of Executable Code:
!  ------------------------
!
!  Replacement collsion only when demanded by user
!
      if( .not.dam%lrepl )  return
!
!  Replacement collision only when projectile and recoil are the same atom species
!
      if( recoil%ielem /= proj%ielem )  return
!
!  Replacement collision only if projectile could not start as a non-virtual recoil
!
      call chkStart( proj, lstart, virtual )
!          ========   i      o        o
      if( lstart .and. .not.virtual )  return
!
!  Replacement collision only if recoil could continue to move as a projectile
!
      call chkStop( recoil, stopped )
!          =======    i        o
      if( stopped )  return
!
!  Replacement collision only if recoil energy is larger than projectile energy
!
      if( recoil%e <= proj%e )  return
!
!  Exchange recoil and projectile, but not recoil generation ig, old position posold,
!  origin, virtual, and state
!
      igp = proj%ig                             ! save properties not to exchange
      igr = recoil%ig
      posoldp(:) = proj%posold(:)
      posoldr(:) = recoil%posold(:)
      originp = proj%origin
      originr = recoil%origin
      virtualp = proj%virtual
      virtualr = recoil%virtual
      virtualparentp = proj%virtualparent
      virtualparentr = recoil%virtualparent
      statep = proj%state
      stater = recoil%state

      proj_tmp = proj                           ! exchange projectile and recoil
      proj = recoil
      recoil = proj_tmp

      proj%ig = igp                             ! restore properties not to exchange
      recoil%ig = igr
      proj%posold(:) = posoldp(:)
      recoil%posold(:) = posoldr(:)
      proj%origin = originp
      recoil%origin = originr
      proj%virtual = virtualp
      recoil%virtual = virtualr
      proj%virtualparent = virtualparentp
      recoil%virtualparent = virtualparentr
      proj%state = statep
      recoil%state = stater

   end subroutine chkRepl
!
!=======================================================================
!
   subroutine chkLeave( proj, dee_, recoil_, posc_, lstart_ )
!             ========   io     i      i       i       i
!
!     Check whether the projectile leaves the target and, if it leaves, whether it reenters.
!     In case it reenters, move projectile to the entrance point.
!
!     Leaving and reentering, or when an error occurs, is indicated by the projectile state
!     proj%state. It may have the values LEAVING, REENTERED, or ERROR.
!
!  Arguments:
!
      type(Proj_type), intent(inout) :: proj          ! projectile
      real(REAL8), optional, intent(in) :: dee_       ! electronic energy loss
      type(Proj_type), optional, intent(in) :: recoil_! recoil
      real(REAL8), optional, intent(in) :: posc_(3)   ! collision point (where impact
                                                      !    parameter is perpendicular to
                                                      !    direction of motion)
      logical, optional, intent(in) :: lstart_        ! flag indicating recoil start trajectory
!
!  Local Variables:
!
      logical :: leaving, error, found
!
!  Start of Executable Code:
!  ------------------------
!
!  Check whether the projectile leaves the target
!
      if( proj%distsurf < - sn%pmax ) then      ! projectile has reached the surface

         call exitTarget( proj, leaving, error )   ! set projectile to a position where
!             ==========   io      o       o       !    distsurf = -sn%pmax

         if( error ) then
            proj%state = STATE%ERROR
            return                    ! no exit point could be determined
         end if

      else

         leaving = .false.

      end if

      if( leaving )  proj%state = STATE%LEAVING
!
!  Evaluation along the trajectory
!
      call eval22( proj, dee_, recoil_, posc_(:), lstart_ )
!          ======   i     i       i        i         i
!
!  Further actions for leaving projectiles
!
      if( leaving ) then

         if( proj%virtual )  proj%ededuct = proj%ededuct + proj%e

                                                ! check whether projectile reenters the target
         if( ( proj%ig == 1 .and. ion%reenter ) .or. ( proj%ig > 1 .and. dam%reenter ) ) then

            call enterTarget( proj, found )
!                ===========   io     o

            if( found ) then
               call eval22( proj )
!                   ======   i
               proj%state = STATE%REENTERED
               proj%reentered = .true.
            end if

         end if

      end if

   end subroutine chkLeave
!
!=======================================================================
!
   subroutine chkStart( recoil, lstart, virtual_ )
!             ========    i       o        o
!
!     Check whether a recoil can start a trajectory
!
!  Arguments:
!
      type(Proj_type), intent(in) :: recoil     ! recoil
      logical, intent(out) :: lstart            ! flag indicating whether trajectory starts
      logical, optional, intent(out) :: virtual_! flag indicating recoil is virtual if it starts
!
!  Local Variables:
!
      logical :: lkp
      real(REAL8) :: elim, ed(nelem)
!
!  Start of Executable Code:
!  ------------------------
!
!  Displacement energy criterion
!
      call getDomainProperty( recoil%id, ed_=ed(:) )
!          =================      i         o

      if( lsput .and. recoil%distsurf < dam%dnear ) then
         elim = min( dam%esurfmin, ed(recoil%ielem) )
         if( recoil%distsurf > HALF*dam%dnear )  &   ! smooth transition close to dnear
               elim = elim + (recoil%distsurf/(HALF*dam%dnear) - ONE) *  &
                     (ed(recoil%ielem) - elim)
      else
         elim = ed(recoil%ielem)
      end if
      lstart = ( recoil%e > elim )
!
!  Kinchin-Pease criterion
!
      call chkKP( recoil, lkp )
!          =====    i      o
      if( lkp )  lstart = .false.

      if( present( virtual_ ) )  virtual_ = .not.lkp .and. recoil%e < ed(recoil%ielem)

   end subroutine chkStart
!
!=======================================================================
!
   subroutine chkStop( proj, stopped, virtual )
!             =======   i       o        o
!
!     Check whether a projectile stops motion based on various criteria.
!     Motion is labeled virtual if energy is below trajectory cutoff energy (possible
!     near the surface)
!
!  Arguments:
!
      type(Proj_type), intent(in) :: proj       ! projectile
      logical, intent(out) :: stopped           ! flag indicating stopped motion
      logical, intent(out), optional :: virtual ! flag indicating projectile turns virtual
!                                               !    if it continues to move
!  Local Variables:
!
      logical :: lkp
      real(REAL8) :: elim
!
!  Start of Executable Code:
!  ------------------------
!
!  Cutoff energy criterion
!
      if( lsput .and. proj%distsurf < dam%dnear ) then
         elim = min( sn%ef, dam%esurfmin )
         if( proj%distsurf > HALF*dam%dnear )  &   ! smooth transition close to dnear
               elim = elim + (proj%distsurf/(HALF*dam%dnear) - ONE) * (sn%ef - elim)
      else
         elim = sn%ef
      end if
      !print *, 'elim, e=', elim, proj%e
      if( sn%lapsisef ) then
         !print *, 'check ', proj%emin
         stopped = ( proj%emin < elim )         ! cutoff energy criterion
      else
         stopped = ( proj%e < elim )            ! cutoff energy criterion
      end if
!
!  Cutoff angle criterion
!
      if( sn%lpsif ) then                       ! cutoff angle criterion
         if( ion%lpinit ) then
            if( abs( -proj%dir(1)*sin(ion%rotate*PI/180)  &
                     + proj%dir(2)*cos(ion%rotate*PI/180) ) > sn%cospsf )  &
!               ===
                  stopped = .true.
         else
            if( ( proj%dir(1)**2 + proj%dir(2)**2 ) > sn%dirl2f )  stopped = .true.
         end if
      end if
!
!  Surface reflection criterion
!
      if( proj%nreflect > sn%reflectmax )  stopped = .true.
!
!  Kinchin-Pease criterion
!
      call chkKP( proj, lkp )
!          =====   i     o
      if( lkp )  stopped = .true.

      if( present( virtual ) ) then
         if( sn%lapsisef ) then
            virtual = ( proj%emin < sn%ef )
         else
            virtual = ( proj%e < sn%ef )
         end if
      end if

   end subroutine chkStop
!
!=======================================================================
!
   subroutine chkKP( proj, lkp )
!             =====   i     o
!
!     Check whether the Kinchin Pease model can be applied
!
!  Arguments:
!
      type(Proj_type), intent(in) :: proj       ! projectile
      logical, intent(out) :: lkp               ! flag indicating Kinchin Pease model can be
                                                !    applied
!
!  Local variables:
!
      real(REAL8) :: dir1(3), dir2(3), dir3(3), dir4(3), dir5(3), dir6(3)
!
!  Start of executable code
!  ------------------------
!
      if( .not.proj%origin%internal .or. proj%e >= dam%elimkp ) then
         lkp = .false.
         return
      end if

      if( proj%id == irxtal .and. irxtal /= 0 .and. dam%cospsi110 /= ONE ) then
                                             ! do not cut off projectiles channeled
                                             ! along <110> or deeper than zmaxkp
         call fromCoordSys( coord_xtal, [ONE,ONE,ZERO], dir1(:) )
!             ============         i                  i            o
         dir1(:) = dir1(:) / abs( dir1(:) )
!                            ===
         call fromCoordSys( coord_xtal, [ONE,-ONE,ZERO], dir2(:) )
!             ============         i                   i            o
         dir2(:) = dir2(:) / abs( dir2(:) )
!                            ===
         call fromCoordSys( coord_xtal, [ONE,ZERO,ONE], dir3(:) )
!             ============         i                  i            o
         dir3(:) = dir3(:) / abs( dir3(:) )
!                            ===
         call fromCoordSys( coord_xtal, [ONE,ZERO,-ONE], dir4(:) )
!             ============         i                  i            o
         dir4(:) = dir4(:) / abs( dir4(:) )
!                            ===
         call fromCoordSys( coord_xtal, [ZERO,ONE,ONE], dir5(:) )
!             ============         i                  i            o
         dir5(:) = dir5(:) / abs( dir5(:) )
!                            ===
         call fromCoordSys( coord_xtal, [ZERO,ONE,-ONE], dir6(:) )
!             ============         i                  i            o
         dir6(:) = dir6(:) / abs( dir6(:) )
!                            ===

         lkp = max( abs( proj%dir(:).dot.dir1(:) ), abs( proj%dir(:).dot.dir2(:) ),  &
!                   ===                             ===
                    abs( proj%dir(:).dot.dir3(:) ), abs( proj%dir(:).dot.dir4(:) ),  &
!                   ===                             ===
                    abs( proj%dir(:).dot.dir5(:) ), abs( proj%dir(:).dot.dir6(:) ) ) <  &
!                   ===                             ===
               SQRT2 * dam%cospsi110  .and.  proj%pos(3) < dam%zmaxkp

      else

         lkp = .true.

      end if

   end subroutine chkKP
!
!=======================================================================
!
end module Casc_Check_mod
