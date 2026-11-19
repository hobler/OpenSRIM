module Xtal_Tvib_mod
!
!  Perform lattice vibrations
!
!  Used Modules:
!  ============
!
   use Const_mod,       only: INT4, REAL8,  &
                              ZERO, ONE, QUART

   use Math_Debye_mod,  only: debye1
   use Random_mod,      only: ranno
   use Undefine_mod

   use Param_Setup_mod, only: nelem
   use Param_Elem_mod,  only: elem
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
   public :: tvib2                  ! perform lattice vibration only perpendicular to
!            =====                  !    direction
   public :: tvib3                  ! perform 3-d lattice vibration
!            =====
   public :: initFakeVib            ! initialize fake lattice vibrations
!            ===========
   public :: addFakeVib             ! add fake lattice vibrations
!            ==========
!
!  Private Data:
!  ============
!
   logical, parameter :: POWERLAW = .false.  ! power law vs logarithmic for fake thermal vib
   logical, allocatable :: lfake(:)          ! (nelem) flag indicating whether to consider
   real(REAL8), allocatable :: efake(:)      ! (nelem) scaling energy (eV)
   real(REAL8), allocatable :: qfake(:)      ! (nelem) exponent for energy dependence
   real(REAL8), allocatable :: vibfake(:)    ! (nelem) 1-d fake vibration amplitude at efake
!
!=======================================================================
!
contains
!
!=======================================================================
!
   subroutine tvib2( vib, dir, posr )
!             =====   i    i    io
!
!     Perform thermal lattice vibration of an atom perpendicular to
!     direction dir (Gaussian probability density).
!
!  Arguments:
!
      real(REAL8), intent(in) :: vib            ! thermal vibration amplitude 
      real(REAL8), intent(in) :: dir(3)         ! direction of motion of the ion
      real(REAL8), intent(inout) :: posr(3)     ! position of target atom
!
!  Local Variables:
!
      real(REAL8) cosalf, sinalf, cosphi, sinphi, rand1, rand2, dxx, dyy
!
!  Start of Executable Code:
!  ------------------------
!
      if( vib == ZERO ) return
!
!  Rotate planar displacement into plane perpendicular to the direction of 
!  motion
!
      cosalf = dir(3)
      sinalf = sqrt( dir(1)**2+dir(2)**2 )

      if( sinalf > 1e-10 ) then
         cosphi = dir(1) / sinalf
         sinphi = dir(2) / sinalf
      else
         cosphi = ONE
         sinphi = ZERO
      end if

      call ranno( rand1 )
!          =====    o
      call ranno( rand2 )
!          =====    o
      dxx = vib*rand1
      dyy = vib*rand2

      posr(1) = posr(1) + cosalf*cosphi*dxx - sinphi*dyy
      posr(2) = posr(2) + cosalf*sinphi*dxx + cosphi*dyy
      posr(3) = posr(3) - sinalf*dxx

      end subroutine tvib2
!
!=======================================================================
!
      subroutine tvib3( vib, posr )
!                =====   i    io
!  
!     Perform thermal lattice vibration of an atom in random direction
!     (Gaussian probability density).
!
!  Arguments:
!
      real(REAL8), intent(in) :: vib            ! thermal vibration amplitude 
      real(REAL8), intent(inout) :: posr(3)     ! position of target atom
!
!  Local Variables:
!
      real(REAL8) :: rand1, rand2, rand3
!
!  Start of Executable Code:
!  ------------------------
!
      if( vib == ZERO ) return

      call ranno( rand1 )
!          =====    o
      call ranno( rand2 )
!          =====    o
      call ranno( rand3 )
!          =====    o
      posr(1) = posr(1) + vib*rand1
      posr(2) = posr(2) + vib*rand2
      posr(3) = posr(3) + vib*rand3

   end subroutine tvib3
!
!=======================================================================
!
   subroutine initFakeVib( )
!             ===========
!
!  Initialize fake lattice vibrations
!
!  This subroutine must be called from a single thread
!
!  Arguments:
!
!     None
!
!  Local Variables:
!
      integer(INT4) :: ielem
      real(REAL8) :: x, u1
!
!  Start of Executable Code:
!  ------------------------
!
      allocate( lfake( nelem ) )
      allocate( efake( nelem ) )
      allocate( qfake( nelem ) )
      allocate( vibfake( nelem ) )

      do ielem = 1, nelem
         lfake(ielem) = defined( elem(ielem)%efake(1) )
         if( .not.lfake(ielem) )  cycle

         x = elem(ielem)%tdebye / elem(ielem)%tempfake(1)
         u1 = 12.1 * sqrt( (debye1( x )/x + QUART) / (elem(ielem)%tdebye * elem(ielem)%m) )
!                           ======  i
         vibfake(ielem) = u1

         if( defined( elem(ielem)%efake(2) ) .and.  &
             elem(ielem)%efake(1) /= elem(ielem)%efake(2) ) then
            efake(ielem) = elem(ielem)%efake(1)
            x = elem(ielem)%tdebye / elem(ielem)%tempfake(2)
            u1 = 12.1 * sqrt( (debye1( x )/x + QUART) / (elem(ielem)%tdebye * elem(ielem)%m) )
!                              ======  i
            if( POWERLAW ) then
               qfake(ielem) = log( u1/vibfake(ielem) ) /  &
                              log( elem(ielem)%efake(2)/elem(ielem)%efake(1) )
            else
               qfake(ielem) = (u1 - vibfake(ielem)) /  &
                              log( elem(ielem)%efake(2)/elem(ielem)%efake(1) )
            end if
         else
            qfake(ielem) = ZERO
         end if
      end do

   end subroutine initFakeVib
!
!=======================================================================
!
   subroutine addFakeVib( vib, e, ielem )
!             ==========  io   i    i
!
!     Add the 1-d fake lattice vibration amplitude quadratically to the passed thermal lattice
!     vibration amplitude
!
!  Arguments:
!
      real(REAL8), intent(inout) :: vib         ! lattice vibration amplitude modified by
                                                !    fake lattice vibration amplitude (A)
      real(REAL8), intent(in) :: e              ! ion energy (eV)
      integer(INT4), intent(in) :: ielem        ! index of target atom
!
!  Local Variables:
!
      real(REAL8) :: u1
!      logical, allocatable :: lfake(:)       ! (nelem) flag indicating whether to consider
!      real(REAL8), allocatable :: efake(:)   ! (nelem) scaling energy (eV)
!      real(REAL8), allocatable :: qfake(:)   ! (nelem) exponent for energy dependence
!      real(REAL8), allocatable :: vibfake(:) ! (nelem) 1-d fake vibration amplitude at efake
!
!  Start of Executable Code:
!  ------------------------
!
      if( lfake(ielem) ) then
         if( qfake(ielem) == ZERO ) then
            u1 = vibfake(ielem)
         else
            if( POWERLAW ) then
               u1 = vibfake(ielem) * (e/efake(ielem))**qfake(ielem)
            else
               u1 = max( vibfake(ielem) + qfake(ielem) * log(e/efake(ielem)), ZERO )
            end if
         end if
         !vib = max( vib, u1 )
         vib = sqrt( vib**2 + u1**2 )
      else
         return
      end if

   end subroutine addFakeVib
!
!=======================================================================
!
end module Xtal_Tvib_mod
