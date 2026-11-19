module Sn_FFPcyl_mod
!
!  Update and set the free flight path cylinder since last collision.
!  There is one cylinder axis (=the projectile path) for all nearby domain, but different
!  cylinder radii for each domain (=pmax).
!
!  Used Modules:
!  ============
!
   use Const_mod,       only: INT4, REAL8,  &
                              ZERO, ONE
   use Array_mod,       only: changeShape
   use Vector_mod,      only: operator(.dot.), abs2

   use Param_Setup_mod, only: nr
   use Param_Geom_mod,  only: geom
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
   public :: initFFPcyl           ! allocate arrays
!            ==========
   public :: newFFPcyl            ! copy ffp cylinder to old ffp cylinder and set the begin
!            =========            !    point of a new ffp cylinder
   public :: setFFPcylRadius      ! set the radius of the ffp cylinder
!            ===============
   public :: isInOldFFPcyl        ! determine if point is inside old ffp cylinder
!            =============
!
!  Private Data:
!  ============
!
   real(REAL8), allocatable :: posbeg(:,:)      ! (3,ng) begin of ffp cylinder axis
   real(REAL8), allocatable :: radius(:,:)      ! (nr,ng) ffp cylinder radii
   real(REAL8), allocatable :: oldposbeg(:,:)   ! (3,ng) posbeg of old ffp cylinders
   real(REAL8), allocatable :: oldposend(:,:)   ! (3,ng) posend of old ffp cylinders
   real(REAL8), allocatable :: oldradius(:,:)   ! (nr,ng) radii of old ffp cylinders
!                                               !   ng initially is 1 and grows as necessary
                                                !   nr may be increased in a dynamic simulation
!$omp threadprivate( posbeg, radius, oldposbeg, oldposend, oldradius )
!
!=======================================================================
!
contains
!
!=======================================================================
!
   subroutine initFFPcyl( )
!             ==========
!
!     Allocate the arrays posbeg, posend, radius
!
!     NOTE: arrays will be reallocated as necessary
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
      allocate( posbeg( 3, 1 ) )
      allocate( oldposbeg( 3, 1 ) )
      allocate( oldposend( 3, 1 ) )
      if( geom%lcell ) then
         allocate( radius( 1, 1 ) )
         allocate( oldradius( 1, 1 ) )
      else
         allocate( radius( nr, 1 ) )
         allocate( oldradius( nr, 1 ) )
      end if

   end subroutine initFFPcyl
!
!=======================================================================
!
   subroutine newFFPcyl( posbeg_new, ig, posend_new_ )
!             =========      i       i       i
!
!     Copy ffp cylinder to old ffp cylinder and set the beginning of the ffp cylinder.
!     Reallocate the arrays if necessary.
!
!  Arguments:
!
      real(REAL8), intent(in) :: posbeg_new(3)              ! beginning of new free flight path
                                                            !   cylinder
      integer(INT4), intent(in) :: ig                       ! index of recoil generation
      real(REAL8), optional, intent(in) :: posend_new_(3)   ! end of free flight path cylinder
!
!  Local Variables:
!
      integer(INT4) :: ng, nd
!
!  Start of Executable Code:
!  ------------------------
!
      ng = size( posbeg, 2 )

      if( ig > ng ) then

         nd = size( radius, 1 )
         call changeShape( posbeg, [3, ig], 'posbeg', 'newFFPcyl' )
!             ===========    io       i        i           i
         call changeShape( radius, [nd, ig], 'radius', 'newFFPcyl' )
!             ===========    io       i         i           i
         call changeShape( oldposbeg, [3, ig], 'oldposbeg', 'newFFPcyl' )
!             ===========     io         i          i            i
         call changeShape( oldposend, [3, ig], 'oldposend', 'newFFPcyl' )
!             ===========     io         i          i            i
         call changeShape( oldradius, [nd, ig], 'oldradius', 'newFFPcyl' )
!             ===========     io         i           i            i
      end if
!
!  Copy ffp cylinder to old ffp cylinder, completing the ffp cylinder with its end position,
!  or clear old ffp cylinder
!
      if( present( posend_new_ ) ) then
         oldposbeg(:,ig) = posbeg(:,ig)
         oldposend(:,ig) = posend_new_(:)
         oldradius(:,ig) = radius(:,ig)
      else
         oldposbeg(:,ig) = ZERO
         oldposend(:,ig) = ZERO
         oldradius(:,ig) = ZERO
      end if
!
!  Set beginning of the new ffp cylinder
!
      posbeg(:,ig) = posbeg_new(:)
      radius(:,ig) = ZERO                       ! these must be set separately

   end subroutine newFFPcyl
!
!=======================================================================
!
   subroutine setFFPcylRadius( radiusl, ig, id )
!             ===============     i     i   i
!
!     Set the radius of the ffp cylinder. Reallocate the array if necessary.
!
!  Arguments:
!
      real(REAL8), intent(in) :: radiusl        ! max impact parameter (id,ig)
      integer(INT4), intent(in) :: ig           ! index of recoil generation
      integer(INT4), intent(in) :: id           ! domain index
!
!  Local Variables:
!
      integer(INT4) :: nd, ng
!
!  Start of Executable Code:
!  ------------------------
!
      ng = size( radius, 2 )
      if( ig > ng )  goto 900
!
!  If id is out of range, reallocate radius
!
      nd = size( radius, 1 )

      if( id > nd ) then

         call changeShape( radius, [id, ng], 'radius', 'setFFPcylRadius' )
!             ===========    io       i
         call changeShape( oldradius, [id, ng], 'oldradius', 'setFFPcylRadius' )
!             ===========     io         i
      end if
!
!  Set new value of radius
!
      radius(id,ig) = radiusl

      return

900   error stop 'setFFPcylRadius: ig out of range'

   end subroutine setFFPcylRadius
!
!=======================================================================
!
   logical function isInOldFFPcyl( pos, ig, id )  result( isInside )
!                   =============   i   i   i                 o
!
!     Determine whether pos(:) is inside the old ffp cylinder
!
!  Arguments:
!
      real(REAL8), intent(in) :: pos(3)         ! position to be tested
      integer(INT4), intent(in) :: ig           ! index of recoil generation
      integer(INT4), intent(in) :: id           ! domain index
!
!      logical, intent(out) :: isInside          ! result
!
!  Local Variables:
!
      real(REAL8) :: ffpold(3), ffpold2, dist(3), lambda
!
!  Start of Executable Code:
!  ------------------------
!
      isInside = .false.

      if( ig > size( radius, 2 ) .or. id > size( radius, 1 ) )  return

      ffpold(:) = oldposend(:,ig) - oldposbeg(:,ig)
      ffpold2 = abs2( ffpold(:) )
!               ====     i
      if( ffpold2 /= ZERO ) then
         lambda = ( (pos(:) - oldposbeg(:,ig)) .dot. ffpold(:) ) / ffpold2
         if( lambda >= ZERO .and. lambda <= ONE ) then
            dist(:) = oldposbeg(:,ig) + lambda * ffpold(:) - pos(:)
            isInside = ( abs2( dist(:) ) <= oldradius(id,ig)**2 )
!                        ====    i
         end if
      end if

   end function isInOldFFPcyl
!
!=======================================================================
!
end module Sn_FFPcyl_mod
