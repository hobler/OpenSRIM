module Xtal_Lattice_Diam_mod
!
!  Generate lattice sites around a given position
!
!     -) The crystal is assumed to have a diamond lattice
!     -) ********************************************************
!        * The lattice positions as well as all coordinates and *
!        * lengths are multiples of a quarter lattice constant  *
!        ********************************************************
!     -) Algorithm after Posselt
!     -) The "search cube" is defined as the smallest cube containing
!        all points which are within the search distance sqrt(nrad2+1)
!        from any point of the unit cube
!
!  Used Modules:
!  ============
!
   use Const_mod,       only: INT4, REAL8,  &
                              ZERO, ONE, TWO, HALF, SMALL

   use Types_mod,       only: Ipos_type

   use Param_Xtal_mod,  only: xtal
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
   public :: initLatticeDiam                 ! initialize the lattice site positions
!            ===============
   public :: maxSizeLatticeDiam              ! return the maximum size of the lattice site
!            ==================              !    list
   public :: getSitesDiam                    ! get list of lattice sites for a diamond
!            ============                    !    lattice

!  Private Procedures:
!  ==================
!
   private :: convertIpos        ! convert old to new lattice indices
!             ===========
!
!  Private Data:
!  ============
!
   integer(INT4) :: npos(2) = [ 0, 0 ]     ! number of lattice sites in each search cell
   integer(INT4), allocatable :: ipos(:,:,:) ! (3,maxval(npos),2) position of lattice sites
                                             !    in search cell
!
!=======================================================================
!
contains
!
!=======================================================================
!
   subroutine initLatticeDiam( )
!             ===============
!
!     Extract search cell for given search radius (nrad2)
!     (Version for unit cubes with edge length 2).
!
!     This subroutine must be called from a single thread.
!
!  Arguments:
!
!     None
!
!  Local Variables:
!
      integer(INT4) :: i, j, icube
      
      integer(INT4), save :: nposTab(0:32,2)
      integer(INT4), save :: iposTab(3,184,2)

      data (nposTab(i,1), i = 0, 32) /  &
           5,  5, 17, 17, 29, 29, 29, 29, 41, 47,  &
          47, 71, 75, 75, 75, 75, 87, 87, 99, 99,  &
         123,123,123,123,135,135,159,167,167,167,  &
         167,167,179 /
      data (nposTab(i,2), i = 0, 32) /  &
           4, 10, 10, 18, 30, 30, 30, 30, 42, 42,  &
          66, 66, 70, 70, 70, 70, 82, 82, 82,106,  &
         130,130,130,130,142,148,148,172,172,172,  &
         172,172,184 /
!
!  Full cube, coordinates with respect to the cube center
!
      data ((iposTab(j,i,1), j = 1, 3), i = 1, 179) /  &
         -1,-1,-1, -1,+1,+1, +0,+0,+0, +1,-1,+1, +1,+1,-1,  &
         -2,-2,+0, -2,+0,-2, -2,+0,+2, -2,+2,+0, +0,-2,-2, +0,-2,+2,  &
         +0,+2,-2, +0,+2,+2, +2,-2,+0, +2,+0,-2, +2,+0,+2, +2,+2,+0,  &
         -3,-1,+1, -3,+1,-1, -1,-3,+1, -1,-1,+3, -1,+1,-3, -1,+3,-1,  &
         +1,-3,-1, +1,-1,-3, +1,+1,+3, +1,+3,+1, +3,-1,-1, +3,+1,+1,  &
         -3,-3,-1, -3,-1,-3, -3,+1,+3, -3,+3,+1, -1,-3,-3, -1,+3,+3,  &
         +1,-3,+3, +1,+3,-3, +3,-3,+1, +3,-1,+3, +3,+1,-3, +3,+3,-1,  &
         -4,+0,+0, +0,-4,+0, +0,+0,-4, +0,+0,+4, +0,+4,+0, +4,+0,+0,  &
         -4,-2,-2, -4,-2,+2, -4,+2,-2, -4,+2,+2, -2,-4,-2, -2,-4,+2,  &
         -2,-2,-4, -2,-2,+4, -2,+2,-4, -2,+2,+4, -2,+4,-2, -2,+4,+2,  &
         +2,-4,-2, +2,-4,+2, +2,-2,-4, +2,-2,+4, +2,+2,-4, +2,+2,+4,  &
         +2,+4,-2, +2,+4,+2, +4,-2,-2, +4,-2,+2, +4,+2,-2, +4,+2,+2,  &
         -3,-3,+3, -3,+3,-3, +3,-3,-3, +3,+3,+3,  &
         -5,-1,-1, -5,+1,+1, -1,-5,-1, -1,-1,-5, -1,+1,+5, -1,+5,+1,  &
         +1,-5,+1, +1,-1,+5, +1,+1,-5, +1,+5,-1, +5,-1,+1, +5,+1,-1,  &
         -4,-4,+0, -4,+0,-4, -4,+0,+4, -4,+4,+0, +0,-4,-4, +0,-4,+4,  &
         +0,+4,-4, +0,+4,+4, +4,-4,+0, +4,+0,-4, +4,+0,+4, +4,+4,+0,  &
         -5,-3,+1, -5,-1,+3, -5,+1,-3, -5,+3,-1, -3,-5,+1, -3,-1,+5,  &
         -3,+1,-5, -3,+5,-1, -1,-5,+3, -1,-3,+5, -1,+3,-5, -1,+5,-3,  &
         +1,-5,-3, +1,-3,-5, +1,+3,+5, +1,+5,+3, +3,-5,-1, +3,-1,-5,  &
         +3,+1,+5, +3,+5,+1, +5,-3,-1, +5,-1,-3, +5,+1,+3, +5,+3,+1,  &
         -5,-3,-3, -5,+3,+3, -3,-5,-3, -3,-3,-5, -3,+3,+5, -3,+5,+3,  &
         +3,-5,+3, +3,-3,+5, +3,+3,-5, +3,+5,-3, +5,-3,+3, +5,+3,-3,  &
         -6,-2,+0, -6,+0,-2, -6,+0,+2, -6,+2,+0, -2,-6,+0, -2,+0,-6,  &
         -2,+0,+6, -2,+6,+0, +0,-6,-2, +0,-6,+2, +0,-2,-6, +0,-2,+6,  &
         +0,+2,-6, +0,+2,+6, +0,+6,-2, +0,+6,+2, +2,-6,+0, +2,+0,-6,  &
         +2,+0,+6, +2,+6,+0, +6,-2,+0, +6,+0,-2, +6,+0,+2, +6,+2,+0,  &
         -4,-4,-4, -4,-4,+4, -4,+4,-4, -4,+4,+4, +4,-4,-4, +4,-4,+4,  &
         +4,+4,-4, +4,+4,+4,  &
         -5,-5,-1, -5,-1,-5, -5,+1,+5, -5,+5,+1, -1,-5,-5, -1,+5,+5,  &
         +1,-5,+5, +1,+5,-5, +5,-5,+1, +5,-1,+5, +5,+1,-5, +5,+5,-1 /
!
!  Empty cube, coordinates with respect to the cube center
!
      data ((iposTab(j,i,2), j = 1, 3), i = 1, 184) /  &
         -1,-1,+1, -1,+1,-1, +1,-1,-1, +1,+1,+1,  &
         -2,+0,+0, +0,-2,+0, +0,+0,-2, +0,+0,+2, +0,+2,+0, +2,+0,+0,  &
         -2,-2,-2, -2,-2,+2, -2,+2,-2, -2,+2,+2, +2,-2,-2, +2,-2,+2,  &
         +2,+2,-2, +2,+2,+2,  &
         -3,-1,-1, -3,+1,+1, -1,-3,-1, -1,-1,-3, -1,+1,+3, -1,+3,+1,  &
         +1,-3,+1, +1,-1,+3, +1,+1,-3, +1,+3,-1, +3,-1,+1, +3,+1,-1,  &
         -3,-3,+1, -3,-1,+3, -3,+1,-3, -3,+3,-1, -1,-3,+3, -1,+3,-3,  &
         +1,-3,-3, +1,+3,+3, +3,-3,-1, +3,-1,-3, +3,+1,+3, +3,+3,+1,  &
         -4,-2,+0, -4,+0,-2, -4,+0,+2, -4,+2,+0, -2,-4,+0, -2,+0,-4,  &
         -2,+0,+4, -2,+4,+0, +0,-4,-2, +0,-4,+2, +0,-2,-4, +0,-2,+4,  &
         +0,+2,-4, +0,+2,+4, +0,+4,-2, +0,+4,+2, +2,-4,+0, +2,+0,-4,  &
         +2,+0,+4, +2,+4,+0, +4,-2,+0, +4,+0,-2, +4,+0,+2, +4,+2,+0,  &
         -3,-3,-3, -3,+3,+3, +3,-3,+3, +3,+3,-3,  &
         -5,-1,+1, -5,+1,-1, -1,-5,+1, -1,-1,+5, -1,+1,-5, -1,+5,-1,  &
         +1,-5,-1, +1,-1,-5, +1,+1,+5, +1,+5,+1, +5,-1,-1, +5,+1,+1,  &
         -4,-4,-2, -4,-4,+2, -4,-2,-4, -4,-2,+4, -4,+2,-4, -4,+2,+4,  &
         -4,+4,-2, -4,+4,+2, -2,-4,-4, -2,-4,+4, -2,+4,-4, -2,+4,+4,  &
         +2,-4,-4, +2,-4,+4, +2,+4,-4, +2,+4,+4, +4,-4,-2, +4,-4,+2,  &
         +4,-2,-4, +4,-2,+4, +4,+2,-4, +4,+2,+4, +4,+4,-2, +4,+4,+2,  &
         -5,-3,-1, -5,-1,-3, -5,+1,+3, -5,+3,+1, -3,-5,-1, -3,-1,-5,  &
         -3,+1,+5, -3,+5,+1, -1,-5,-3, -1,-3,-5, -1,+3,+5, -1,+5,+3,  &
         +1,-5,+3, +1,-3,+5, +1,+3,-5, +1,+5,-3, +3,-5,+1, +3,-1,+5,  &
         +3,+1,-5, +3,+5,-1, +5,-3,+1, +5,-1,+3, +5,+1,-3, +5,+3,-1,  &
         -5,-3,+3, -5,+3,-3, -3,-5,+3, -3,-3,+5, -3,+3,-5, -3,+5,-3,  &
         +3,-5,-3, +3,-3,-5, +3,+3,+5, +3,+5,+3, +5,-3,-3, +5,+3,+3,  &
         -6,+0,+0, +0,-6,+0, +0,+0,-6, +0,+0,+6, +0,+6,+0, +6,+0,+0,  &
         -6,-2,-2, -6,-2,+2, -6,+2,-2, -6,+2,+2, -2,-6,-2, -2,-6,+2,  &
         -2,-2,-6, -2,-2,+6, -2,+2,-6, -2,+2,+6, -2,+6,-2, -2,+6,+2,  &
         +2,-6,-2, +2,-6,+2, +2,-2,-6, +2,-2,+6, +2,+2,-6, +2,+2,+6,  &
         +2,+6,-2, +2,+6,+2, +6,-2,-2, +6,-2,+2, +6,+2,-2, +6,+2,+2,  &
         -5,-5,+1, -5,-1,+5, -5,+1,-5, -5,+5,-1, -1,-5,+5, -1,+5,-5,  &
         +1,-5,-5, +1,+5,+5, +5,-5,-1, +5,-1,-5, +5,+1,+5, +5,+5,+1 /
!
!  Start of Executable Code:
!  ------------------------
!
!  Check crystal name
!
      if( xtal%uc%name /= ' ' )  &
            error stop 'initLatticeDiam: entered with arbitrary crystal structure'
!
!  Calculate the positions of the stored lattice sites 
!
      do icube = 1, 2
         npos(icube) = nposTab(xtal%nrad2,icube)
      end do

      allocate( ipos( 3, maxval(npos(:)), 2 ) )
      ipos(1:3,1:npos(1),1) = iposTab(1:3,1:npos(1),1)
      ipos(1:3,1:npos(2),2) = iposTab(1:3,1:npos(2),2)
!
   end subroutine initLatticeDiam
!
!===========================================================================================
!
   integer(INT4) function maxSizeLatticeDiam( ) 
!                         ==================
!
!     Return an upper limit to the size of the lattice site list to be returned by
!     getSitesDiam
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
      maxSizeLatticeDiam = maxval(npos(:))

      if( maxSizeLatticeDiam == 0 )  error stop 'maxSizeLatticeDiam = 0'

   end function maxSizeLatticeDiam
!
!=======================================================================
!
   subroutine getSitesDiam( pos, nl, iposl, cub )
!             ============   i   o     o     o 
!
!     Set up a list of lattice sites around a given position using a lattice 
!     site list around the unit cube. (Version for unit cubes with edge 
!     length 2).
!
!     -) The crystal is assumed to have a diamond lattice
!     -) ********************************************************
!        * The lattice positions as well as all coordinates and *
!        * lengths are multiples of a quarter lattice constant  *
!        ********************************************************
!     -) Algorithm after Posselt
!     -) The "search cube" is defined as the smallest cube containing
!        all points which are within the search distance sqrt(nrad2+1)
!        from any point of the unit cube
!
!  Arguments:
!
      real(REAL8), intent(in) :: pos(3)            ! position
      integer(INT4), intent(out) :: nl             ! number of lattice sites
      type(Ipos_type), intent(inout) :: iposl(:)   ! lattice site
      real(REAL8), intent(out) :: cub(3,2)         ! lower (:,1) and upper (:,2) corner of search cube
!
!  Local Variables:
!
      integer(INT4) :: l, icube
      integer(INT4) :: icen(3), ih(3)
      real(REAL8) :: dx
      real(REAL8) :: posh(3)
!
!  Start of Executable Code:
!  ------------------------
!
!  Determine index of subcube containing x
!
      posh(:) = HALF*pos(:)
      ih(:) = nint( posh(:) - HALF )
!
!  Determine whether full cube (icube=1) or empty cube (icube=2)
!
      icube = abs( mod( sum( ih(:) ), 2 ) ) + 1
!
!  Determine coordinates of the center of the cube
!
      icen(:) = 2 * ih(:) + 1
!
!  Calculate the positions of the stored lattice sites 
!
      nl = npos(icube)
      do l = 1, nl
         iposl(l)%ipos(:) = icen(:) + ipos(:,l,icube)
         call convertIpos( iposl(l) )           ! convert the position to new indices
!             ===========    io
      end do

      dx = xtal%rrad + ONE                      ! using xtal%rrad(nrad2) here limits 
                                                !    reusability of this code
      cub(:,1) = icen(:) - dx
      cub(:,2) = icen(:) + dx

   end subroutine getSitesDiam
!
!===========================================================================================
!
   subroutine convertIpos( ipos )
!             ===========   io
!
!     Convert old to new lattice site indices
!
!  Arguments:
!
      type(Ipos_type), intent(inout) :: ipos    ! lattice site indices
!
!  Local Variables:
!
      integer(INT4) :: l, m, n
!
!  Start of Executable Code:
!  ------------------------
!
      ipos%iuc(:) = ( ipos%ipos(:) - modulo( ipos%ipos(:), 4 ) ) / 4

      l = modulo( ( ipos%ipos(1) - modulo( ipos%ipos(1), 2 ) ) / 2, 2 )
      m = modulo( ( ipos%ipos(2) - modulo( ipos%ipos(2), 2 ) ) / 2, 2 )
      n = modulo( ipos%ipos(1), 2 )
      ipos%iucs = l + 2*m + 4*n

   end subroutine convertIpos
!
!=======================================================================
!
end module Xtal_Lattice_Diam_mod
