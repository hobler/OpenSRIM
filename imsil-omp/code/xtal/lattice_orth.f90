module Xtal_Lattice_Orth_mod
!
!  Generate lattice sites around a given position
!
!     -) The crystal is assumed to have an orthogonal unit cell with arbitrary
!        positions of the lattice sites
!     -) ********************************************************
!        * The lattice positions as well as all coordinates and *
!        * lengths are multiples of the lattice constant        *
!        ********************************************************
!     -) The "search cube" is defined as the smallest cuboid containing
!        all points which are within the search distance rrad
!        from any point of the unit cell
!
!  Used Modules:
!  ============
!
   use Const_mod,       only: INT4, REAL8,  &
                              ZERO, ONE, TWO, HALF, SMALL
   use Stack_mod,       only: Int4Stack,  &
                              newStack, sizeOfStack, pushStack, popStack
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
   public :: initLatticeOrth         ! initialize the lattice site positions
!            ===============
   public :: maxSizeLatticeOrth      ! return the maximum size of the lattice site list
!            ==================
   public :: getSitesOrth   ! get list of lattice sites for a cubic lattice
!            ============
!
!  Private Data:
!  ============
!
   integer(INT4) :: npos = 0                 ! number of lattice sites in search cell
   type(Ipos_type), allocatable :: ipos(:)   ! (npos) position of lattice sites
                                             !    in search cell
   real(REAL8) :: cub(3,2)                   ! search cube
!
!  Program History:
!  ===============
!
!   2-Jul-2003  G.Hobler  creation from subroutines lsunit, lslist, fndnls
!   4-Jan-2004  G.Hobler  conversion to free form
!   9-Jun-2004  G.Hobler  remove lslist, rename lsunit to getLSList
!  25-Jan-2007  G.Hobler  adapt from xtal_lattice_mod; use Ipos_type and allocatable ipos
!  26-Jan-2007  G.Hobler  getSitesOrth, maxSizeLatticeOrth, initLatticeOrth
!  12-Mar-2007  G.Hobler  renamed getLSList to getSites, fndnls to getNearestSite
!
!=======================================================================
!
contains
!
!=======================================================================
!
   subroutine initLatticeOrth( )
!             ===============
!
!     Extract search cell for given search radius (rrad)
!
!     This subroutine must be called from a single thread
!
!  Arguments:
!
!     None
!
!  Local Variables:
!
      logical :: flag, isEmpty
      integer(INT4) :: n, ix, iy, iz, isite
      real(REAL8) :: pos(3), dpos(3), dist2
      type(Int4Stack) :: stack1, stack2, stack3, stack4
!
!  Start of Executable Code:
!  ------------------------
!
!  Check crystal name
!
      if( xtal%uc%name == ' ' )  &
            error stop 'initLatticeOrth: entered with optimized diamond crystal lattice'
!
!  Calculate the positions of the stored lattice sites
!
      call newStack( stack1 )                   ! stacks for position
!          ========    o
      call newStack( stack2 )
!          ========    o
      call newStack( stack3 )
!          ========    o
      call newStack( stack4 )
!          ========    o

      n = -1
      do                                        ! go in shells outwards until no new
         n = n + 1                              !    site is found
         flag = .false.
         do ix = -n, n
            do iy = -n, n
               do iz= -n, n
                  if( abs(ix) < n .and. abs(iy) < n .and. abs(iz) < n )  cycle
                  do isite = 1, xtal%uc%nsite
                     pos(:) = [ ix, iy, iz ] + xtal%uc%pos(:,isite)
                     where( pos > ONE )
                        dpos = pos - ONE
                     elsewhere
                        where( pos < ZERO )
                           dpos = pos
                        elsewhere
                           dpos = ZERO
                        end where
                     end where
                     !dist2 = sum( dpos(:)**2 )  ! only lower limit in case of non-cubic unit
                                                !    cell (assuming rrad given in units of
                                                !    the smallest lattice constant)
                     dist2 = sum( (xtal%uc%alat(:) * dpos(:))**2 )

                     !if( dist2 < xtal%rrad2 ) then
                     if( dist2 < minval(xtal%uc%alat(:))**2 * xtal%rrad2 ) then
                        flag = .true.
                        call pushStack( ix, stack1 )  ! put lattice site indices on stack
!                            =========  i     io
                        call pushStack( iy, stack2 )
!                            =========  i     io
                        call pushStack( iz, stack3 )
!                            =========  i     io
                        call pushStack( isite, stack4 )
!                            =========    i      io
                     end if
                  end do
               end do
            end do
         end do
         if( .not.flag )  exit                  ! exit if no new lattice site has been found
      end do

      npos = sizeOfStack( stack1 )              ! get size of stack (all stacks have same size)
!            ===========    i
      allocate( ipos( npos ) )                  ! allocate array
      
      do isite = 1, npos                        ! copy stacks to array
         call popStack( stack1, ipos(isite)%iuc(1), isEmpty )
!             ========    io           o               o
         if( isEmpty )  error stop 'initLatticeOrth: inconsistent stack size'
         call popStack( stack2, ipos(isite)%iuc(2), isEmpty )
!             ========    io           o               o
         if( isEmpty )  error stop 'initLatticeOrth: inconsistent stack size'
         call popStack( stack3, ipos(isite)%iuc(3), isEmpty )
!             ========    io           o               o
         if( isEmpty )  error stop 'initLatticeOrth: inconsistent stack size'
         call popStack( stack4, ipos(isite)%iucs, isEmpty )
!             ========    io           o             o
         if( isEmpty )  error stop 'initLatticeOrth: inconsistent stack size'
         !print '(4i3)', ipos(isite)%iuc, ipos(isite)%iucs
      end do
      !print *, npos

      cub(:,1) = - xtal%rrad
      cub(:,2) = ONE + xtal%rrad
!
   end subroutine initLatticeOrth
!
!===========================================================================================
!
   integer(INT4) function maxSizeLatticeOrth( ) 
!                         ==================
!
!     Return an upper limit to the size of the lattice site list to be returned by
!     getSitesOrth
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
      maxSizeLatticeOrth = npos

      if( maxSizeLatticeOrth == 0 )  error stop 'maxSizeLatticeOrth = 0'

   end function maxSizeLatticeOrth
!
!=======================================================================
!
   subroutine getSitesOrth( pos, nl, iposl, cubl )
!             ============   i   o     o     o 
!
!     Set up a list of lattice sites around a given position using a lattice 
!     site list around the unit cube. (Version for unit cubes with edge 
!     length 2).
!
!     -) The crystal is assumed to have a cubic lattice
!     -) ********************************************************
!        * The lattice positions as well as all coordinates and *
!        * lengths are multiples of the lattice constant        *
!        ********************************************************
!     -) The "search cube" is defined as the smallest cube containing
!        all points which are within the search distance rrad
!        from any point of the unit cube
!
!  Arguments:
!
      real(REAL8), intent(in) :: pos(3)         ! position
      integer(INT4), intent(out) :: nl          ! number of lattice sites
      type(Ipos_type), intent(out) :: iposl(:)  ! lattice site
      real(REAL8), intent(out) :: cubl(3,2)     ! lower (:,1) and upper (:,2) corner of search
                                                !    cube
!
!  Local Variables:
!
      integer(INT4) :: i, l, icube(3)
!
!  Start of Executable Code:
!  ------------------------
!
!  Determine index of subcube containing x
!
      icube(:) = nint( pos(:) - HALF )
!
!  Calculate the positions of the stored lattice sites 
!
      nl = npos
      do l = 1, nl
         iposl(l)%iuc(:) = icube(:) + ipos(l)%iuc(:)
         iposl(l)%iucs = ipos(l)%iucs
      end do

      do i = 1, 2
         cubl(:,i) = icube(:) + cub(:,i)
      end do

   end subroutine getSitesOrth
!
!===========================================================================================
!
end module Xtal_Lattice_Orth_mod
