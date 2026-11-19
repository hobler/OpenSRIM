module Xtal_Ipos_mod
!
!  Provides a class for lattice site indices
!
!  Used Modules:
!  ============
!
   use Const_mod,          only: INT4, REAL8
   use Types_mod,          only: Ipos_type

   use Param_Mat_mod,      only: coord_xtal
   use Param_Xtal_mod,     only: xtal

   use Geom_CoordSys_mod,  only: fromCoordSys
!
!  Default Type and Access:
!  =======================
!
   implicit none
   private
   save
!
!  Public Procedures and Operators:
!  ===============================
!
   public :: initIpos                        ! set ipos to defined value (all 0)
!            ========
   public :: getPosI                         ! get the position of the lattice site in the
!            ======                          !   crystal coordinate system
   public :: getPos                          ! get the position of the lattice site in the 
!            ======                          !   target coordinate system

   public :: operator( == )

   interface operator( == )
      module procedure equal
   end interface
!
!=======================================================================
!
contains
!
!=======================================================================
!
   subroutine initIpos( ipos )
!             ========   o
!
!     Initialize ipos to defined values (lattice site at the origin)
!
!  Arguments:
!
      type(Ipos_type), intent(out) :: ipos       ! lattice site indices
!
!  Start of Executable Code:
!  ------------------------
!
      ipos%ipos(:) = 0
      ipos%iuc(:) = 0
      ipos%iucs = 1

   end subroutine initIpos
!
!=======================================================================
!
   subroutine getPosI( ipos, posi )
!             =======   i     o
!
!     Convert lattice site indices to position in the crystal coordinate system
!     Coordinates will still be scaled by the (quarter) lattice constant
!
!  Arguments:
!
      type(Ipos_type), intent(in) :: ipos       ! lattice site indices
      real(REAL8), intent(out) :: posi(3)       ! coordinates of the position (crystal
                                                !    coordinate system)
!
!  Local Variables:
!
!     None
!
!  Start of Executable Code:
!  ------------------------
!
      if( xtal%uc%name == ' ' ) then            ! default diamond lattice
         posi(:) = ipos%ipos(:)
      else                                      ! arbitrary lattice
         posi(:) = ipos%iuc(:) + xtal%uc%pos(:,ipos%iucs)
      end if

   end subroutine getPosI
!
!=======================================================================
!
   subroutine getPos( ipos, ir, pos )
!             ======   i    i    o
!
!     Convert lattice site indices to position in the target coordinate system
!     Coordinates will still be scaled by the (quarter) lattice constant
!
!  Arguments:
!
      type(Ipos_type), intent(in) :: ipos       ! lattice site indices
      integer(INT4), intent(in) :: ir           ! index of region
      real(REAL8), intent(out) :: pos(3)        ! coordinates of the position (target
                                                !    coordinate system)
!
!  Local Variables:
!
      real(REAL8) :: posi(3)
!
!  Start of Executable Code:
!  ------------------------
!
      call getPosI( ipos, posi )
!          =======   i     o

      call fromCoordSys( coord_xtal, posi(:), pos(:) )
!          ============        i          i        o

   end subroutine getPos
!
!=======================================================================
!
   logical function equal( ipos1, ipos2 )
!                   =====    i      i
!
!     Test whether ipos1 and ipos2 are identical
!     Assumes that iuc and iucs are always set so it is sufficient to test their values
! 
!  Arguments:
!
      type(Ipos_type), intent(in) :: ipos1      ! atom indices of first lattice site
      type(Ipos_type), intent(in) :: ipos2      ! atom indices of second lattice site
!
!  Start of Executable Code:
!  ------------------------
!
      equal = all( ipos1%iuc(:) == ipos2%iuc(:) ) .and. ipos1%iucs == ipos2%iucs

   end function equal
!
!=======================================================================
!
end module Xtal_Ipos_mod
