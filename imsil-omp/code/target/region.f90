module Target_Region_mod
!
!  Holds region structure for static targets and provides subroutines operating on it.
!
!  Used Modules:
!  ============
!
   use Const_mod,                   only: INT4, REAL8

   use Types_mod,                   only: Proj_type

   use Param_Setup_mod,             only: ndim, nr
   use Param_Geom_mod,              only: geom
   use Param_Mat_mod,               only: mat

   use Target_Region1d_mod,         only: getRegionIndex1d,  &
                                          getOtherNearRegionIndices1d, sputter1d
   use Target_Region2d_mod,         only: initRegion2d, getRegionIndex2d,  &
                                          getOtherNearRegionIndices2d
   use Target_Region3d_mod,         only: initRegion3d, getRegionIndex3d,  &
                                          getOtherNearRegionIndices3d
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
   public :: initRegion                ! initialize the region target
!            ==========
   public :: getRegionIndex            ! determine the region nearest to point pos
!            ==============
   public :: getOtherNearRegionIndices ! determine regions other than nearest region that
!            ========================= !    intersect with the free flight path cylinder
   public :: sputter                   ! erode the surface due to sputtering
!            =======
!
contains
!
!=======================================================================
!
   subroutine initRegion( )
!             ==========
!
!     Initialize the region target
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
      if( ndim == 1 ) then
         continue                               ! nothing to initialize
      else if( ndim == 2 ) then
         call initRegion2d( )
!             ============
      else
         call initRegion3d( )
!             ============
      end if

   end subroutine initRegion
!
!=======================================================================
!
   subroutine getRegionIndex( pos, ir_, outside_ )
!             ==============   i   io      o
!
!     Determine the nearest region to the point pos
!
!  Arguments:
!
      real(REAL8), intent(in) :: pos(3)            ! position to be checked
      integer(INT4), optional, intent(inout) :: ir_! region index
      logical, optional, intent(out) :: outside_   ! flag indicating whether pos(:) is outside
                                                   !    all regions (thus ir_ being only the
                                                   !    region nearest to ir_)
!
!  Local Variables:
!
!     None
!
!  Start of Executable Code:
!  ------------------------
!
      if( ndim == 1 ) then
         call getRegionIndex1d( pos, ir_, outside_ )
!             ================   i   io      o
      else if( ndim == 2 ) then
         call getRegionIndex2d( pos, ir_, outside_ )
!             ================   i   io      o
      else
         call getRegionIndex3d( pos, ir_, outside_ )
!             ================   i   io      o
      end if

   end subroutine getRegionIndex
!
!=======================================================================
!
   subroutine getOtherNearRegionIndices( proj, posc, ffp, pmax, nrn, irnear )
!             =========================   i     i     i    i     o     o
!
!     Determine regions other than proj%id that are within the free flight path cylinder
!     (= cylinder with height ffp and radius pmax)
!
!  Arguments:
!
      type(Proj_type), intent(in) :: proj       ! projectile
      real(REAL8), intent(in) :: posc(3)        ! position of end of free flight path
      real(REAL8), intent(in) :: ffp            ! free flight path (= height of cylinder)
      real(REAL8), intent(in) :: pmax           ! maximum impact parameter (= radius of cyl.)
      integer(INT4), intent(out) :: nrn         ! number of regions found
      integer(INT4), allocatable, intent(inout) :: irnear(:)
                                                ! index of regions found
!
!  Local Variables:
!
!     None
!
!  Start of Executable Code:
!  ------------------------
!
      if( ndim == 1 ) then
         call getOtherNearRegionIndices1d( proj, posc, ffp, pmax, nrn, irnear )
!             ===========================   i     i     i    i     o     o
      else if( ndim == 2 ) then
         call getOtherNearRegionIndices2d( proj, posc, ffp, pmax, nrn, irnear )
!             ===========================   i     i     i    i     o     o
      else
         call getOtherNearRegionIndices3d( proj, posc, ffp, pmax, nrn, irnear )
!             ===========================   i     i     i    i     o     o
      end if

   end subroutine getOtherNearRegionIndices
!
!=======================================================================
!
   subroutine sputter( ielem )
!             =======    i
!
!     Erode the surface due to sputtering
!
!  Arguments:
!
      integer(INT4), intent(in) :: ielem        ! index of ion atom
!
!  Local Variables:
!
!     None
!
!  Start of Executable Code:
!  ------------------------
!
      if( geom%lcell )  return

      if( ndim == 1 ) then
         call sputter1d( ielem )
!             =========    i
      else
         if( any( mat(1:nr)%lsput ) )  &
            error stop 'sputter: not implemented for ndim > 1'
      end if

   end subroutine sputter
!
!=======================================================================
!
end module Target_Region_mod
