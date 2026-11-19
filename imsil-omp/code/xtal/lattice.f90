module Xtal_Lattice_mod
!
!  Generate lattice sites and perform closely related operations
!
!  Used Modules:
!  ============
!
   use Const_mod,             only: INT4, REAL8,  &
                                    ZERO, ONE, TWO, HALF, SMALL
   use Vector_mod,            only: abs2
   use Types_mod,             only: Ipos_type

   use Param_Mat_mod,         only: coord_xtal
   use Param_Xtal_mod,        only: xtal

   use Geom_CoordSys_mod,     only: toCoordSys
   use Xtal_Ipos_mod,         only: getPos
   use Xtal_Lattice_Diam_mod, only: initLatticeDiam, maxSizeLatticeDiam, getSitesDiam
   use Xtal_Lattice_Orth_mod, only: initLatticeOrth, maxSizeLatticeOrth, getSitesOrth
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
   public :: initLattice         ! initialize the stored lattices site search cell
!            ===========
   public :: maxSizeLattice      ! return the maximum size of the lattice site list
!            ==============
   public :: getSites            ! get list of lattice sites
!            ========
   public :: getNearestSite      ! find the lattice site nearest to a given position
!            ==============
!
!=======================================================================
!
contains
!
!=======================================================================
!
   subroutine getNearestSite( pose, ir, iposl )
!             ==============   i    i     o
!
!     Find nearest lattice site
!
!  Arguments:
!
      real(REAL8), intent(in) :: pose(3)        ! input position (target coordinate system)
      integer(INT4), intent(in) :: ir           ! index of region
      type(Ipos_type), intent(out) :: iposl     ! nearest lattice site (indices)
!
!  Local Variables:
!
      type(Ipos_type), save, allocatable :: jposl(:)
      integer(INT4) :: il, nl
      real(REAL8) :: dist2, dist2l
      real(REAL8) :: pos(3), posl(3), cub(3,2)
!$omp threadprivate( jposl )
!
!  Start of Executable Code:
!  ------------------------
!
!  Calculate a list of lattice sites surrounding the ion position
!
      if( .not. allocated( jposl ) )  &
            allocate( jposl( maxsizeLattice() ) )
!                            ==============
      call toCoordSys( coord_xtal, pose(:), pos(:) )
!          ==========        i          i       o
      call getSites( pos(:), nl, jposl(:), cub(:,:) )
!          ========    i     o      o         o
!
!  Select the nearest ion position
!
      call getPos( jposl(1), ir, posl(:) )
!          ======    i       i     o
      dist2 = abs2( posl(:)-pose(:) )
!             ====
      iposl = jposl(1)
      
      do il = 2, nl
         call getPos( jposl(il), ir, posl(:) )
!             ======    i          i    o
         dist2l = abs2( posl(:)-pose(:) )
!                 ====
         if( dist2l < dist2 ) then
            dist2 = dist2l
            iposl = jposl(il)
         end if
      end do

   end subroutine getNearestSite
!
!=======================================================================
!
   subroutine initLattice( )
!             ===========
!
!     Initialize lattice site search cell
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
      if( xtal%uc%name == ' ' ) then
         call initLatticeDiam( )
!             ===============
      else
         call initLatticeOrth( )
!             ===============
      end if

   end subroutine initLattice
!
!=======================================================================
!
   integer(INT4) function maxSizeLattice( ) 
!                         ==============
!
!     Return an upper limit to the size of the lattice site list to be returned by
!     getSites
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
      if( xtal%uc%name == ' ' ) then
         maxSizeLattice = maxSizeLatticeDiam( )
!                         ==================
      else
         maxSizeLattice = maxSizeLatticeOrth( )
!                         ==================
      end if

   end function maxSizeLattice
!
!=======================================================================
!
   subroutine getSites( posi, nl, iposl, cub )
!             ========   i    o     o     o
!
!     Set up a list of lattice sites around a given position
!
!  Arguments:
!
      real(REAL8), intent(in) :: posi(3)           ! position (crystal coordinate system)
      integer(INT4), intent(out) :: nl             ! number of lattice sites
      type(Ipos_type), intent(inout) :: iposl(:)   ! lattice site
      real(REAL8), intent(out) :: cub(3,2)         ! lower (:,1) and upper (:,2) corner of search cube
!
!  Local Variables:
!
!     None
!
!  Start of Executable Code:
!  ------------------------
!
      if( xtal%uc%name == ' ' ) then
         call getSitesDiam( posi, nl, iposl, cub )
!             ============    i   o     o     o
      else
         call getSitesOrth( posi, nl, iposl, cub )
!             ============    i   o     o     o
      end if

   end subroutine getSites
!
!=======================================================================
!
end module Xtal_Lattice_mod
