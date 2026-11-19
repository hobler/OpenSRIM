module Eval_Flow_mod
!
!  Stores and analyzes vacancies and interstitials generated in one ion's recoil cascade
!
!  Note: currently no valid model. Kept for possible later use.
!
!  Used Modules:
!  ============
!
   use Const_mod,          only: INT4,REAL8,  &
                                 ZERO, ONE, TWO, THIRD, BIG,  &
                                 PI
   use Array_mod,          only: check_allocation
   use Vector_mod,         only: abs2
   use Stack_mod,          only: ProjStack,  &
                                 newStack, clearStack, pushStack, stack2Array
   use Types_mod,          only: Proj_type

   use Param_Setup_mod,    only: ldyn
   use Param_Elem_mod,     only: elem
   use Param_Dam_mod,      only: dam

   use Eval_Dam_mod,       only: dam20, dam21, dam23
   use Eval_Stat_mod,      only: stat2x

   use Geom_BndCond_mod,   only: applyBndCond3d

   use Target_Cell_mod,    only: scoreCell
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
   public :: flow10                 ! Initialize at beginning of simulation
!            ======                  
   public :: flow11                 ! Initialize at beginning of molecule cascade
!            ======
   public :: flow20                 ! Evaluate subthreshold recoil
!            ======
   public :: flow21                 ! Evaluate at beginning of trajectory
!            ======
   public :: flow23                 ! Evaluate at end of trajectory of a stopped ion
!            ======
   public :: flow26                 ! Actions at end of cascade
!            ======
!
!  Private Data:
!  ============
!
   type(ProjStack) :: IStack        ! Interstitial stack
!$omp threadprivate( IStack )
!
!=======================================================================
!
contains
!
!=======================================================================
!
   subroutine flow10( )
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
!  Initialize stacks
!
      call newStack( IStack )
!          ========    o
!
   end subroutine flow10
!
!=======================================================================
!
   subroutine flow11( )
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
!  Start of Executable Code:
!  ------------------------
!
!     Stack should always be clear at the beginning of a recoil cascade, if
!     the previous cascade has been handled properly. But just to make sure
!     we clear it again.
!
      call clearStack( IStack )
!          ==========    io

   end subroutine flow11
!
!=======================================================================
!
   subroutine flow20( proj )
!             ======   i
!
!     Evaluate subthreshold collision
!
!  Arguments:
!
      type(Proj_type), intent(in) :: proj       ! vacancy
!
!  Local Variables:
!
      logical :: allow_score
      real(REAL8) :: ddam, edam, x, y, z
!
!  Start of Executable Code:
!  ------------------------
!
      if( proj%e < dam%eminflow ) then          ! recoil does not flow
                                                !    do statistics as if lflow=.false.
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

      else                                      ! treat vacancy and interstitial as if there
                                                !    were a trajectory
         call flow21( proj )
!             ======   i

         call flow23( proj )
!             ======   i

      end if

   end subroutine flow20
!
!=======================================================================
!
   subroutine flow21( proj )
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
      real(REAL8) :: ddami, ddamv, x, y, z
!
!  Start of Executable Code:
!  ------------------------
!
      call dam21( proj, ddami, ddamv )
!          =====   i      o      o

      call stat2x( proj, ddami, ddamv, ZERO, 1 )
!          ======   i      i      i     i    i

      if( ldyn ) then
         x = proj%pos(1)
         y = proj%pos(2)
         z = proj%pos(3)
         call applyBndCond3d( x, y, z, .false., allow_score )
!             ============== io io io     i          o
         if( allow_score ) then           ! proj%id should already be the correct cell index
            call scoreCell( proj%ielem, proj%id, -ONE, ddami, ddamv, ZERO )
!                =========       i         i       i     i      i     i
         end if
      end if

   end subroutine flow21
!
!=======================================================================
!
   subroutine flow23( proj )
!             ======   i
!
!     Evaluate at the end of a trajectory
!
!  Arguments:
!
      type(Proj_type), intent(in) :: proj       ! interstitial
!
!  Local Variables:
!
!     None
!
!  Start of Executable Code:
!  ------------------------
!
      call pushStack( proj, IStack )
!          =========   i      io

   end subroutine flow23
!
!=======================================================================
!
   subroutine flow26( )
!             ======
!
!     Perform actions after simulation of a collision cascade
!
!     - move interstitial according to its energy
!     - move interstitial depending on momentum of nearby interstitials
!
!  Arguments: 
!
!     None
!
!  Local Variables:
!
      logical :: allow_score
      character(256) :: errmsg
      integer(INT4) :: ii, ji, ni, istat
      real(REAL8) :: ddami, ddamv, edam, x, y, z
      real(REAL8) :: r2
      real(REAL8), allocatable :: esum(:), psum(:,:)
      type(Proj_type), allocatable :: I(:)
!
!  Start of Executable Code:
!  ------------------------
!
!  Convert stack to array, get size, and allocate esum
!
      call stack2Array( IStack, I )
!          ===========    io    o

      ni = size( I )
!
!  Increase range of projectiles below cutoff energy
!
      do ii = 1, ni
         I(ii)%pos(:) = I(ii)%pos(:) + dam%drdeflow * (I(ii)%e - dam%eminflow ) * I(ii)%dir(:)
      end do
!
!  Enhance range due to nearby deposited energy
!
      if( dam%lflow .and. dam%facflow > ZERO ) then
         allocate( esum( ni ), stat=istat, errmsg=errmsg )
         call check_allocation( 'flow26', 'esum', istat, errmsg )

         do ii = 1, ni
            esum(ii) = ZERO
            do ji = 1, ni
               if( ji == ii )  cycle
               r2 = abs2( I(ii)%pos(:) - I(ji)%pos(:) )
!                   ====               i
               if( r2 < dam%rcap2 )  esum(ii) = esum(ii) + I(ji)%e
            end do
            !print *, esum(ii)
         end do

         do ii = 1, ni
            I(ii)%pos(:) = I(ii)%pos(:) + dam%facflow * esum(ii) * I(ii)%dir(:)
         end do

         deallocate( esum )
      end if
!
!  Enhance range due to nearby deposited momentum
!
      if( dam%lflow .and. dam%pfacflow > ZERO ) then
         allocate( psum( 3, ni ), stat=istat, errmsg=errmsg )
         call check_allocation( 'flow26', 'psum', istat, errmsg )

         do ii = 1, ni
            psum(:,ii) = ZERO
            do ji = 1, ni
               if( ji == ii )  cycle
               r2 = abs2( I(ii)%pos(:) - I(ji)%pos(:) )
!                   ====               i
               if( r2 < dam%rcap2 )  &
                     psum(:,ii) = psum(:,ii) + I(ji)%dir(:) *  &
                                  sqrt( TWO * elem(I(ji)%ielem)%m * I(ji)%e )
            end do
            !print *, psum(ii)
         end do

         do ii = 1, ni
            I(ii)%pos(:) = I(ii)%pos(:) + dam%pfacflow * psum(:,ii)
         end do

         deallocate( psum )
      end if

!
!  Evaluate statistics
!
      do ii = 1, ni

         call dam23( I(ii), ddami, ddamv, edam )
!             =====    i      o      o     o

         call stat2x( I(ii), ddami, ddamv, edam, 3 )
!             ======    i      i      i     i    i

         if( ldyn ) then
            x = I(ii)%pos(1)
            y = I(ii)%pos(2)
            z = I(ii)%pos(3)
            call applyBndCond3d( x, y, z, .false., allow_score )
!                ============== io io io     i          o
            if( allow_score ) then           ! proj%id should already be the correct cell index
               call scoreCell( I(ii)%ielem, I(ii)%id, ONE, ddami, ddamv, I(ii)%e )
!                   =========       i          i       i     i      i       i
            end if
         end if
      end do
!
!  Deallocate arrays
!
      deallocate( I )
!
!-----------------------------------------------------------------------
!
   end subroutine flow26
!
!=======================================================================
!
end module Eval_Flow_mod
