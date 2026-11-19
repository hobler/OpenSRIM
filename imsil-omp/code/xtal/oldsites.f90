module Xtal_OldSites_mod
!
!  Update and set old site list
!
!  Used Modules:
!  ============
!
   use Const_mod,       only: INT4, REAL8,  &
                              ZERO
   use Array_mod,       only: changeShape
   use Types_mod,       only: Proj_type, Ipos_type

   use Param_Xtal_mod,  only: xtal 

   use Xtal_Ipos_mod,   only: operator( == )
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
   public :: initOldSites           ! set the lattice site the ion originates from
!            ============
   public :: updateOldSites         ! update the old site list
!            ==============
   public :: copyOldSites           ! copy old sites from one recoil generation to another
!            ============
   public :: testOldSites           ! check whether a site is contained in the old site list
!            ============
!
!  Private Procedures:
!  ==================
!
   private :: changeShapeOldSites   ! change shape of the arrays defining old lattice sites
!             ===================
!
!  Private Data:
!  ============
!
   integer(INT4), parameter :: NSOMAX = 20
   integer(INT4), allocatable :: nsold(:)       ! (NGM)
   type(Ipos_type), allocatable :: iposold(:,:) ! (NSOMAX,NGM)
   integer(INT4), allocatable :: isaold(:,:)    ! (NSOMAX,NGM)
   real(REAL8), allocatable :: ffpold(:)        ! (NGM)
                                                ! sum of ffps since last crystalline collision
!$omp threadprivate( nsold, iposold, isaold, ffpold )
!
!=======================================================================
!
contains
!
!=======================================================================
!
   subroutine initOldSites( )
!             ============
!
!     Allocate the arrays to store the old sites
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
      if( allocated( nsold ) ) then
         deallocate( nsold )
         deallocate( iposold )
         deallocate( isaold )
         deallocate( ffpold )
      end if

      allocate( nsold( 1 ) )
      allocate( iposold( xtal%nloldm, 1 ) )
      allocate( isaold( xtal%nloldm, 1 ) )
      allocate( ffpold( 1 ) )

      nsold(1) = 0
      ffpold(1) = ZERO

   end subroutine initOldSites
!
!=======================================================================
!
   subroutine copyOldSites( ig1, ig2 )
!             ============   i    i
!
!     Copy old sites from ig1 to ig2
!
!  Arguments:
!
      integer(INT4), intent(in) :: ig1          ! index of recoil generation of source
      integer(INT4), intent(in) :: ig2          ! index of recoil generation of destination
!
!  Local Variables:
!
!
!  Start of Executable Code:
!  ------------------------
!
!  Extend the old lattice site list if necessary
!
      if( ig2 > size( nsold ) )  call changeShapeOldSites( ig2 )
!                                     ===================   i
      nsold(ig2) = nsold(ig1)
      iposold(:,ig2) = iposold(:,ig1)
      isaold(:,ig2) = isaold(:,ig1)
      ffpold(ig2) = ffpold(ig1)

   end subroutine copyOldSites
!
!=======================================================================
!
   subroutine changeShapeOldSites( ig )
!             ===================  i
!
!     Change the shape of the old site list
!
!  Arguments:
!
      integer(INT4), intent(in) :: ig           ! index of recoil generation
!
!  Local Variables:
!
!
!  Start of Executable Code:
!  ------------------------
!
      call changeShape( nsold, ig, 'nsold', 'changeShapeOldSites' )
!          ===========   io    i      i               i
      call changeShape( iposold, [xtal%nloldm, ig], 'iposold', 'changeShapeOldSites' )
!          ===========    io             i              i                i
      call changeShape( isaold, [xtal%nloldm, ig], 'isaold', 'changeShapeOldSites' )
!          ===========    io            i             i                i
      call changeShape( ffpold, ig, 'ffpold', 'changeShapeOldSites' )
!          ===========    io    i      i                i

   end subroutine changeShapeOldSites
!
!=======================================================================
!
   subroutine updateOldSites( ig, ffp_, cospsi_, proj_ )
!             ==============  i    i       i       i
!
!     Update old site list
!
!     The actions performed depend on the presence of the optional arguments.
!     Note, in particular, that the old lattice site list for recoil generation ig
!     is deleted in the following cases:
!
!     - ffp_ is not passed as an argument.
!     - the accumulated ffp_ exceeds xtal%ffoldm
!     - the scattering angle exceeds acos(xtal%cpsolm)
!
!  Arguments:
!
      integer(INT4), intent(in) :: ig                 ! index of recoil generation
      real(REAL8), optional, intent(in) :: ffp_       ! free flight path
      real(REAL8), optional, intent(in) :: cospsi_    ! cosine of scattering angle
      type(Proj_type), optional, intent(in) :: proj_  ! projectile
!
!  Local Variables:
!
      integer(INT4) :: i
!
!  Start of Executable Code:
!  ------------------------
!
!  Extend the old lattice site list if necessary
!  (Note that changeShape sets new elements to zero)
!
      if( ig > size( nsold ) )  call changeShapeOldSites( ig )
!                                    ===================  i
!
!  Delete old lattice site for this ig, if accumulated ffp_ exceeds ffoldm or ffp_ is not passed
!
      if( present( ffp_ ) ) then
         ffpold(ig) = ffpold(ig) + ffp_
         if( ffpold(ig) > xtal%ffoldm ) then
            nsold(ig) = 0
            ffpold(ig) = ZERO
         end if
      else
         nsold(ig) = 0
         ffpold(ig) = ZERO
      end if
!
!  Delete old lattice site list for this ig, if cos(scrattering angle) exceeds cpsolm
!
      if( present( cospsi_ ) ) then
         if( cospsi_ < xtal%cpsolm ) then
            nsold(ig) = 0
            ffpold(ig) = ZERO
         end if
      end if
!
!  If the proj_ originates from a crystalline region, store it in the old sites list
!
      if( present( proj_ ) ) then
         if( .not.proj_%origin%amorphous ) then
            nsold(ig) = min( xtal%nloldm, nsold(ig)+1 )
            do i = nsold(ig), 2, -1             ! make room at position 1 for new site
               iposold(i,ig) = iposold(i-1,ig)
               isaold(i,ig) = isaold(i-1,ig)
            end do

            iposold(1,ig) = proj_%origin%ipos    ! store new site at position 1
            isaold(1,ig) = proj_%origin%isa

            ffpold(ig) = ZERO
         end if
      end if

   end subroutine updateOldSites
!
!=======================================================================
!
   subroutine testOldSites( ipos, isa, ig, found )
!             ============    i    i   i     o
!
!     Check whether a site is contained in the old site list
!
!  Arguments:
!
      type(Ipos_type), intent(in) :: ipos       ! position of site to be checked
      integer(INT4), intent(in) :: isa          ! atom counter in split interstitial
                                                !    configurations
      integer(INT4), intent(in) :: ig           ! index of recoil generation
      logical, intent(out) :: found             ! flag indicating whether site has been found
                                                !    in list
!
!  Local Variables:
!
      integer(INT4) :: i
!
!  Start of Executable Code:
!  ------------------------
!
      do i = 1, nsold(ig)
         if( i > size( iposold, 1 ) .or. ig>size( iposold, 2 ) ) then
            print *, 'i,ig,shape(iposold)=', i, ig, shape(iposold)
            print *, 'nsold(ig)=', nsold(ig)
         end if
         if( ipos == iposold(i,ig) .and. isa == isaold(i,ig) .and. isa /= -1 ) then
            found = .true.
            return
         end if
      end do

      found = .false.

   end subroutine testOldSites
!
!=======================================================================
!
end module Xtal_OldSites_mod
