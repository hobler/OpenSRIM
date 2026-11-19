module Sn_RestFFP_mod
!
!  Update and set the remaining free flight paths for each recoil generation and domain.
!  This is needed when the projectile moves within pmax of more than one domain. Then,
!  after a collision partner in one domain has been selected, rest ffps remain for all
!  the other domains.
!
!  Used Modules:
!  ============
!
   use Const_mod,       only: INT4, REAL8,  &
                              ZERO
   use Param_Setup_mod, only: nr
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
   public :: initRestFFP             ! allocate rest ffp and pmax arrays
!            ===========
   public :: resetRestFFP            ! reset the rest ffp list
!            ============
   public :: updateRestFFP           ! update rest ffp list considering that projectile
!            =============           !    proceeds specified ffp
   public :: setRestFFP              ! set one particular rest ffp
!            ==========
   public :: getRestFFP              ! get one particular rest ffp
!            ==========
!
!  Private Data:
!  ============
!
   integer(INT4), parameter :: NGMAX = 10    ! maximum number of recoil generations to
                                             ! be stored. Higher generations can probably
                                             ! be neglected because they have low energy
                                             ! and thus small ffps
   real(REAL8), allocatable :: restffp(:,:)
   real(REAL8), allocatable :: restpmax(:,:)
!$omp threadprivate( restffp, restpmax )
!
!=======================================================================
!
contains
!
!=======================================================================
!
   subroutine initRestFFP( )
!             ===========
!
!     Allocate the array restffp, restpmax
!
!     NOTE: arrays will be reallocated if necessary
!
!     This subroutine must be called from a parallel region
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
      allocate( restffp( nr, NGMAX ) )
      allocate( restpmax( nr, NGMAX ) )

   end subroutine initRestFFP
!
!=======================================================================
!
   subroutine resetRestFFP( ig )
!             ============  i
!
!     Reset the the ffp list for one recoil generation
!
!  Arguments:
!
      integer(INT4), intent(in) :: ig           ! recoil generation index
!
!  Local Variables:
!
!     None
!
!  Start of Executable Code:
!  ------------------------
!
      if( ig > NGMAX )  return

      restffp(:,ig) = ZERO
      restpmax(:,ig) = ZERO

   end subroutine resetRestFFP
!
!=======================================================================
!
   subroutine updateRestFFP( ig, ffp )
!             =============  i    i
!
!     Update old ffp list by subtracting ffp from each restffp(:,ig)
!
!  Arguments:
!
      integer(INT4), intent(in) :: ig           ! index of recoil generation
      real(REAL8), intent(in) :: ffp           ! free flight path to be subtracted
!
!  Local Variables:
!
!     None
!
!  Start of Executable Code:
!  ------------------------
!
      if( ig > NGMAX )  return
!
!  Reduce all ffp's of recoil generation ig by value ffp
!
      restffp(:,ig) = restffp(:,ig) - ffp
!
   end subroutine updateRestFFP
!
!=======================================================================
!
   subroutine setRestFFP( id, ig, ffp, pmax )
!             ==========  i   i    i    i
!
!     Set ffp(id,ig) and pmax(id,ig) to the passed values
!
!  Arguments:
!
      integer(INT4), intent(in) :: id                 ! domain index
      integer(INT4), intent(in) :: ig                 ! index of recoil generation
      real(REAL8), intent(in) :: ffp                  ! free flight path (id,ig)
      real(REAL8), intent(in) :: pmax                 ! max impact parameter (id,ig)
!
!  Local Variables:
!
      integer(INT4) :: nr_tmp
      real(REAL8), allocatable :: tmp(:,:)
!
!  Start of Executable Code:
!  ------------------------
!
      if( ig > NGMAX )  return
!
!  If id is out of range, reallocate restffp and restpmax
!
      if( id > size( restffp, 1 ) ) then
         nr_tmp = size( restffp, 1 )

         allocate( tmp(nr_tmp,NGMAX) )

         tmp(:,:) = restffp(:,:)                 ! reallocate restffp and copy values

         deallocate( restffp )
         allocate( restffp(id,NGMAX) )

         restffp(1:nr_tmp,:) = tmp(:,:)
         restffp(nr_tmp+1:,:) = ZERO
         
         tmp(:,:) = restpmax(:,:)                ! reallocate restpmax and copy values

         deallocate( restpmax )
         allocate( restpmax(id,NGMAX) )

         restpmax(1:nr_tmp,:) = tmp(:,:)
         restpmax(nr_tmp+1:,:) = ZERO

         deallocate( tmp )
      end if         
!
!  Set new values of restffp(id,ig) and restpmax(id,ig)
!
      restffp(id,ig) = ffp
      restpmax(id,ig) = pmax

   end subroutine setRestFFP
!
!=======================================================================
!
   subroutine getRestFFP( id, ig, ffp, pmax )
!             ==========  i   i    o    o
!
!     Get the ffp values at (id,ig)
!
!  Arguments:
!
      integer(INT4), intent(in) :: id           ! domain index
      integer(INT4), intent(in) :: ig           ! index of recoil generation
      real(REAL8), intent(out) :: ffp           ! free flight path
      real(REAL8), intent(out) :: pmax          ! maximum impact parameter
!
!  Local Variables:
!
!     None
!
!  Start of Executable Code:
!  ------------------------
!
      if( ig > NGMAX .or. id > size( restffp, 1 ) ) then
         ffp = ZERO
         pmax = ZERO
         return
      end if

      ffp = restffp(id,ig)
      pmax = restpmax(id,ig)

   end subroutine getRestFFP
!
!=======================================================================
!
end module Sn_RestFFP_mod
