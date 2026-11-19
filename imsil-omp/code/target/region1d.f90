module Target_Region1d_mod
!
!  Check point location or free flight paths with respect to the regions for 1d geometry.
!
!  Used Modules:
!  ============
!
   use Const_mod,             only: INT4, REAL8,  &
                                    ZERO, HALF
   use Types_mod,             only: Proj_type

   use Param_Setup_mod,       only: nr
   use Param_Geom_mod,        only: geom
   use Param_Ion_mod,         only: ion
   use Param_Mat_mod,         only: mat

   use Geom_BndCond_mod,      only: applyBndCond1d
   use Geom_CurviCoord_mod,   only: toCurviCoord1d

   use Target_Surface1d_mod,  only: updateSurface1d
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
   public :: getRegionIndex1d             ! determines the region a given point is
!            ================             !   located in
   public :: getOtherNearRegionIndices1d  ! determine regions other than nearest region that
!            ===========================  !    intersect with the free flight path cylinder
   public :: sputter1d                    ! erode the surface (1-d)
!            =========
!
!=======================================================================
!
contains
!
!=======================================================================
!
   subroutine getRegionIndex1d( pos, ir_, outside_ )
!             ================   i   io      o
!
!     Determine the region the point pos is located in
!
!  Arguments:
!
      real(REAL8), intent(in) :: pos(3)            ! position to be checked
      integer(INT4), optional, intent(inout) :: ir_! index of region
      logical, optional, intent(out) :: outside_   ! flag indicating whether pos is outside any
                                                   !    region
!
!  Local Variables:
!
      real(REAL8) :: z
!
!  Start of Executable Code:
!  ------------------------
!
      if( geom%lcart ) then
         z = pos(3)
      else
         call toCurviCoord1d( pos(:), z )
!             ==============    i     o
      end if
      call applyBndCond1d( z )
!          ============== io
!
!  Check whether pos is inside any region
!
      if( present( outside_ ) )  &
            outside_ = ( z < geom%posif(0) .or. z > geom%posif(nr) )
!
!  Determine layer z is located in
!
      if( present( ir_ ) ) then
         do
            if( geom%posif(ir_) < z ) then
               if( ir_ == nr )  exit
               ir_ = ir_ + 1
            else if( geom%posif(ir_-1) > z ) then
               if( ir_ == 1 )  exit
               ir_ = ir_ - 1
            else
               exit
            end if
         end do
      end if

   end subroutine getRegionIndex1d
!
!=======================================================================
!
   subroutine getOtherNearRegionIndices1d( proj, posc, ffp, pmax, nrn, irnear )
!             ===========================   i     i     i    i     o     o
!
!     Determine regions other than ir that must be checked for potential collision partners.
!
!     For planar layers, regions are returned that overlap the free flight path cylinder
!     (= cylinder with height ffp and radius pmax). For cylindrical and spherical layers,
!     regions are returned that overlap the sphere circumscribing the free flight path
!     cylinders.
!
!     Only regions in the range (1:nr) are returned. Vacuum regions (0, nr+1) are not
!     returned.
!
!  Arguments:
!
      type(Proj_type), intent(in) :: proj       ! projectile
      real(REAL8), intent(in) :: posc(3)        ! position of end of free flight path
      real(REAL8), intent(in) :: ffp            ! free flight path (= height of cylinder)
      real(REAL8), intent(in) :: pmax           ! maximum impact parameter (= radius of cyl.)
      integer(INT4), intent(out) :: nrn         ! number of regions found
      integer(INT4), intent(out) :: irnear(:)   ! index of regions found except proj%id
!
!  Local Variables:
!
      integer(INT4) :: ir, jr
      real(REAL8) :: z, dz, z1, z2, zmin, zmax, zcen, distother, ffpcen(3)
!
!  Start of Executable Code:
!  ------------------------
!
!  If there is only one region there is no "other near region"
!
      if( nr == 1 ) then
         nrn = 0
         return
      end if

      ir = proj%id                              ! region to be excluded from the search

      if( geom%lcart ) then
         z = proj%pos(3)
      else
         call toCurviCoord1d( proj%pos(:), z )
!             ==============       i       o
      end if
      call applyBndCond1d( z )
!          ============== io
!
!  First do a rough check whether other regions are far from the projectile
!
      if( z < geom%posif(1) ) then           
         distother = geom%posif(1) - z
      else if( z > geom%posif(nr-1) ) then
         distother = z - geom%posif(nr-1)
      else
         distother = min( geom%posif(ir) - z, z - geom%posif(ir-1) )
      end if

      if( distother < ZERO )  error stop 'getOtherNearRegions1d: inconsistent ir'

      if( distother**2 > ffp**2 + pmax**2 ) then
         nrn = 0
         return
      end if
!
!  Determine minimum and maximum z-values of the ffp cylinder or the sphere circumscribing
!  the ffp cylinder
!
      if( geom%lcart ) then
         dz = pmax * sqrt( proj%dir(1)**2 + proj%dir(2)**2 )
         z1 = proj%pos(3) - sign( dz, proj%dir(3) )
         z2 = posc(3) + sign( dz, proj%dir(3) )
         zmin = min( z1, z2 )
         zmax = max( z1, z2 )
      else
         ffpcen(:) = HALF * ( proj%pos(:) + posc(:) )
         call toCurviCoord1d( ffpcen(:), zcen )
!             ==============     i        o
         dz = sqrt( (HALF*ffp)**2 + pmax**2 )
         zmin = zcen - dz
         zmax = zcen + dz
      end if
!
!  Starting from ir, search for regions intersecting the free flight path cylinder
!
      nrn = 0

      do jr = ir, 2, -1
         if( geom%posif(jr-1) < zmin )  exit
         nrn = nrn + 1
         irnear(nrn) = jr - 1
      end do

      do jr = ir, nr-1
         if( geom%posif(jr) > zmax )  exit      
         nrn = nrn + 1
         irnear(nrn) = jr + 1
      end do

   end subroutine getOtherNearRegionIndices1d
!
!=======================================================================
!
   subroutine sputter1d( ielem )
!             =========  i
!
!     Erode the surface due to sputtering (1-d)
!
!  Arguments:
!
      integer(INT4), intent(in) :: ielem        ! index of ion atom
!
!  Local Variables:
!
      integer(INT4) :: ir
      real(REAL8) :: dposif
!
!  Start of Executable Code:
!  ------------------------
!
      if( mat(1)%lsput ) then

         do ir = 1, nr-1
            if( geom%posif(ir-1) /= geom%posif(ir) )  exit
         end do

         dposif = mat(ir)%syield(ielem) * ion%ddose / mat(ir)%dens

         geom%posif(0) = geom%posif(0) + dposif
         do ir = 1, nr
            geom%posif(ir) = max( geom%posif(ir), geom%posif(0) )
         end do

      end if

   end subroutine sputter1d
!
!=======================================================================
!
end module Target_Region1d_mod
