module Target_Region2d_mod
!
!  Geometry operations on regions defined by polygons.
!
!  To increase efficiency, a regular grid is maintained whose cells contain information on
!  nearby regions. The grid is shared with Target_Surface2d_mod. Cell indices are obtained
!  by calling getCellIndexSurface2d from Target_Surface2d_mod.
!
!  Used Modules:
!  ============
!
   use Const_mod,                   only: INT4, INT8, REAL8,  &
                                          ZERO, HALF, BIG
   use Array_mod,                   only: alloc, check_allocation
   use Types_mod,                   only: Proj_type

   use Param_Setup_mod,             only: nr
   use Param_Geom_mod,              only: geom

   use Geom_BndCond_mod,            only: applyBndCond2d
   use Geom_CurviCoord_mod,         only: toCurviCoord2d
   use Math_IntMath_mod,            only: index2to1

   use Target_Region2d_Poly_mod,    only: isOutsideRegion2dPoly, getDist2RectRegion2dPoly,  &
                                          getDist2FromRegion2dPoly
   use Target_Region2d_SDF_mod,     only: initRegion2dSDF, getDist2RectRegion2dSDF,  &
                                          isOutsideRegion2dSDF, getDist2FromRegion2dSDF
#ifdef GTS
   use Target_Region2d_SDFOct_mod,  only: getDist2RectRegion2dSDFOct, getRegionIndex2dSDFOct
#endif
   use Target_Surface2d_mod,        only: initSurface2d, getNcellSurface2d,  &
                                          getCellIndexSurface2d
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
   public :: initRegion2d                 ! do initializations
!            ============
   public :: getRegionIndex2d             ! determine the region a given point is located in
!            ================             !    (returns 0 if outside target)
   public :: getOtherNearRegionIndices2d  ! determine regions other than nearest region that
!            ===========================  !    intersect with the free flight path cylinder
!
!  Private Procedures:
!  ==================
!
   private :: getDist2FromRegion2d        ! determine square of distance from region
!             ====================
   private :: getDist2RectRegion2d        ! determine square of distance of rectangle from
!             ====================        !    region
   private :: isOutsideRegion2d           ! determine if point is outside region
!             =================
!
!  Private Data:
!  ============
!
   type :: Cell_type
!          =========
      integer(INT4) :: nrovlp                ! number of regions the cell overlaps with
      integer(INT4), allocatable :: ir(:)    ! (nr) list of regions ordered according to
                                             !    increasing distmin2
                                             !    ir(1) corresponds to the region the cell
                                             !    center is located in; the cell overlaps with
                                             !    regions ir(1:nrovlp)
      real(REAL8), allocatable :: distmin2(:)! (nr) square of minimum minimum distance of
                                             !    point inside the cell from boundary of
                                             !    region ir
   end type Cell_type

   type(Cell_type), allocatable :: cell(:)   ! (grid%nx*grid%nz) cells
!
!=======================================================================
!
contains
!
!=======================================================================
!
   subroutine initRegion2d( )
!             ============
!
!     Do initializations
!
!     We use the same grid for surface and regions
!
!  Arguments:
!
!     None
!
!  Local Variables:
!
      logical :: cinside
      character(256) :: errmsg
      integer(INT4) :: i, j, ir, ih, i1, irc, istat
      integer(INT8) :: ic, ncell
      real(REAL8) :: dh, distmin2(nr)
!
!  Start of Executable Code:
!  ------------------------
!
!  Initialize the surface
!
      call getNcellSurface2d( ncell )
!          =================    o
!
!  Read the region SDF files
!
      if( geom%usesdf .and. .not.geom%usehyperoctree )  call initRegion2dSDF( )
!                                                            ===============
!  Initialize the regions (can be skipped if nr=1 since cell will not be used in this case)
!
      if( nr == 1 )  return

      allocate( cell( ncell ), stat=istat, errmsg=errmsg )
      call check_allocation( 'initRegion2d', 'cell', istat, errmsg )
      do ic = 1, ncell
         call alloc( cell(ic)%ir, nr, 'cell(ic)%ir', 'initRegion2d' )
!             =====       o       i         i              i
         call alloc( cell(ic)%distmin2, nr, 'cell(ic)%ir', 'initRegion2d' )
!             =====          o          i         i              i
      end do
!
!  Calculate distance of cells from regions
!
      do ic = 1, ncell

         irc = 0
         do ir = 1, nr
            call getDist2RectRegion2d( ic, ir, cinside, cell(ic)%distmin2(ir) )
!                ====================  i   i      i          o
            if( cinside )  irc = ir
         end do
!
!  Put the region the cell center is located in at position 1 and sort remaining regions
!  in ascending order of distmin2
!
         cell(ic)%ir(:) = [ (ir, ir = 1, nr) ]

         if( irc > 0 ) then
            ih = cell(ic)%ir(irc)               ! swap region cell center is in
            cell(ic)%ir(irc) = cell(ic)%ir(1)
            cell(ic)%ir(1) = ih
            dh = cell(ic)%distmin2(irc)
            cell(ic)%distmin2(irc) = cell(ic)%distmin2(1)
            cell(ic)%distmin2(1) = dh
            i1 = 2
         else                                   ! nothing to swap since point is outside
            i1 = 1                              !    or point position is unknown
         end if

         do i = i1, nr-1                        ! sort regions with respect to distmin2
            j = i + minloc( cell(ic)%distmin2(i+1:nr), dim=1 )    ! dim=1 needed to get scalar
            if( cell(ic)%distmin2(j) < cell(ic)%distmin2(i) ) then
               ih = cell(ic)%ir(i)              ! swap
               cell(ic)%ir(i) = cell(ic)%ir(j)
               cell(ic)%ir(j) = ih
               dh = cell(ic)%distmin2(i)
               cell(ic)%distmin2(i) = cell(ic)%distmin2(j)
               cell(ic)%distmin2(j) = dh
            end if
         end do
!
         cell(ic)%nrovlp = count( distmin2(:) == ZERO )

      end do

   end subroutine initRegion2d
!
!=======================================================================
!
   subroutine getRegionIndex2d( pos, ir_, outside_ )
!             ================   i   io      o
!
!     Determine the region the point pos is located in (if pos is inside the target) or closest
!     to (if pos is outside the target)
!
!     "outside" here refers to the position being outside any target region. For polygon-type
!     regions, for numerical reasons, this is not always exactly the same as distsurf < 0.
!
!  Arguments:
!
      real(REAL8), intent(in) :: pos(3)            ! position to be checked
      integer(INT4), optional, intent(inout) :: ir_! index of region (on entry: guess)
      logical, optional, intent(out) :: outside_   ! flag indicating whether pos is outside any
!                                                  !    region
!  Local Variables:
!
      logical :: outside
      integer(INT4) :: jr, irnot, i
      integer(INT8) :: ic
      real(REAL8) :: x, z, dist2, dist2_
!
!  Start of Executable Code:
!  ------------------------
!
      if( geom%lcart ) then
         x = pos(1)
         z = pos(3)
      else
         call toCurviCoord2d( pos(:), x, z )
!             ==============    i     o  o
      end if
      call applyBndCond2d( x, z )
!          ============== io io

      if( geom%usehyperoctree ) then               ! octree does all the work itself
#ifdef GTS
         call getRegionIndex2dSDFOct( x, z, ir_, outside_ )
!             ======================  i  i  io      o
         return
#else
         error stop 'getRegionIndex2d: usehyperoctree=t not supported'
#endif
      end if
!
!  If the target has only one region, return ir_=1, since ir_=1 is either the region pos is
!  located in or nearest to.
!
      if( nr == 1 ) then
         if( present( ir_ ) )  ir_ = 1
         if( present( outside_ ) )   call isOutsideRegion2d( x, z, 1, outside_ )
!                                         =================  i  i  i     o
         return
      end if
!
!  Determine index of cell [x,z] is located in
!
      call getCellIndexSurface2d( x, z, ic, outside, 'x,z(getRegionIndex)' )
!          =====================  i  i  o      o                i
      if( ic < 1 .or. ic > size( cell ) )  error stop 'getRegionIndex2d: ic out of range'
!
!  If cell overlaps with only one region, return this one
!
      if( cell(ic)%nrovlp == 1 ) then
         jr = cell(ic)%ir(1)
         if( present( ir_ ) )  ir_ = jr
         if( present( outside_ ) )  call isOutsideRegion2d( x, z, jr, outside_ )
!                                        =================  i  i  i      o
         return
      end if
!
!  If cell overlaps with more than one region, first check ir_ passed to the subroutine
!
      if( cell(ic)%nrovlp > 1 ) then
         if( present( ir_ ) ) then
            call isOutsideRegion2d( x, z, ir_, outside )
!                =================  i  i  i      o
            if( present( outside_ ) )  outside_ = outside
            if( .not.outside )  return
            irnot = ir_                         ! index that does not need to be checked again
         else
            irnot = -1                          ! all indices need to be checked
         end if
!
!  ... then check the other regions the cell overlaps with
!
         do i = 1, cell(ic)%nrovlp
            jr = cell(ic)%ir(i)
            if( jr == irnot )  cycle
            call isOutsideRegion2d( x, z, jr, outside )
!                =================  i  i  i      o
            if( .not.outside ) then
               if( present( ir_ ) )  ir_ = jr
               if( present( outside_ ) )  outside_ = .false.
               return
            end if
         end do
      end if
!
!  If this point is reached, position is outside the target. Check distance to all regions
!  to determine region with minimum distance
!
      if( present( outside_ ) )  outside_ = .true.

      if( present( ir_ ) ) then
         dist2 = BIG
         do jr = 1, nr
            call getDist2FromRegion2d( x, z, jr, dist2_ )
!                ====================  i  i  i     o
            if( dist2_ < dist2 ) then
               ir_ = jr
               dist2 = dist2_
            end if
         end do
      end if

   end subroutine getRegionIndex2d
!
!=======================================================================
!
   subroutine getOtherNearRegionIndices2d( proj, posc, ffp, pmax, nrn, irnear )
!             ===========================   i     i     i    i     o     o
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
      integer(INT4), intent(out) :: irnear(:)   ! index of regions found
!
!  Local Variables:
!
      logical :: mask(nr), outside
      integer(INT8) :: ic
      integer(INT4) :: irnot, i
      real(REAL8) :: xc, zc, x, z, dirx, diry, dirz
      real(REAL8) :: r2, ffp_p
!
!  Start of Executable Code:
!  ------------------------
!
!  If there is only one region, there is no "other near region"
!
      if( nr == 1 ) then
         nrn = 0
         return
      end if

      irnot = proj%id                              ! region to be excluded from the search
      mask(:) = .true.
      mask(irnot) = .false.

      if( irnot < 1 .or. irnot > nr )  &
            error stop 'getOtherNearRegionIndices2d: irnot must be in [1,nr]'
!
!  First do a rough check using the cell the projectile is in
!
      if( geom%lcart ) then
         x = proj%pos(1)
         z = proj%pos(3)
      else
         call toCurviCoord2d( proj%pos(:), x, z )
!             ==============       i       o  o
      end if
      call applyBndCond2d( x, z )
!          ============== io io

      call getCellIndexSurface2d( x, z, ic, outside, 'x,z' )
!          =====================  i  i  o      o       i
      if( ic < 1 .or. ic > size( cell ) )  &
            error stop 'getOtherNearRegionIndices2d (x): ic out of range'
      if( outside ) then
!         print *, 'x,z=', x, z
         error stop 'getOtherNearRegionIndices2d: position outside grid'
      end if

      if( minval( cell(ic)%distmin2(:), mask=mask(:) ) > ffp**2 + pmax**2 ) then
         nrn = 0
         return
      end if
!
!  Get cell index of center point of ffp
!
      if( geom%lcart ) then
         xc = HALF*(proj%pos(1)+posc(1))
         zc = HALF*(proj%pos(3)+posc(3))
      else
         call toCurviCoord2d( HALF*(proj%pos(:)+posc(:)), xc, zc )
!             ==============               i              o   o
      end if
      call applyBndCond2d( xc, zc )
!          ==============

      call getCellIndexSurface2d( xc, zc, ic, outside, 'xc,zc' )
!          =====================  i   i   o      o        i
      if( ic < 1 .or. ic > size( cell ) )  &
            error stop 'getOtherNearRegionIndices2d (xc): ic out of range'
!
!  If center point is outside the grid, we must consider all regions except irnot
!
      if( outside ) then
         nrn = nr - 1                           ! we must consider all regions except ir
         irnear(1:nrn) = [ ( i, i=1,irnot-1 ), ( i, i=irnot+1,nr ) ]
         return
      end if
!
!  If center point is inside the grid, list all cells within a circle of radius sqrt(r2)
!
      r2 = pmax**2 + (HALF*ffp)**2
      diry = proj%dir(2)
      if( r2*diry**2 > pmax**2 ) then
         dirx = proj%dir(1)
         dirz = proj%dir(3)
         ffp_p = ffp * sqrt( dirx**2 + dirz**2 )
         r2 = ( HALF*ffp_p + pmax*abs( diry ) )**2
      end if

      nrn = 0
      do i = 1, nr
         if( cell(ic)%ir(i) == irnot )  cycle
         if( cell(ic)%distmin2(i) > r2 )  exit  ! We can exit because distmin2 are sorted
         nrn = nrn + 1
         irnear(nrn) = cell(ic)%ir(i)
      end do
!
!  TODO: For long ffp_p it might be advisable to subdivide ffp_p into several segments and use
!        a circle in each segment ...
!
   end subroutine getOtherNearRegionIndices2d
!
!=======================================================================
!
   subroutine getDist2FromRegion2d( x, z, ir, dist2 )
!             ====================  i  i  i     o
!
!     Determine the square of the distance of the point [x,z] from region ir
!
!  Arguments:
!
      real(REAL8), intent(in) :: x               ! x coordinate of point
      real(REAL8), intent(in) :: z               ! z coordinate of point
      integer(INT4), intent(in) :: ir           ! index of region
      real(REAL8), intent(out) :: dist2         ! square of the distance from region
!
!  Local Variables:
!
!     None
!
!  Start of Executable Code:
!  ------------------------
!
      if( geom%usehyperoctree ) then
         error stop 'getDist2FromRegion2dSDFOct not implemented'
      else if( geom%usesdf ) then
         call getDist2FromRegion2dSDF( x, z, ir, dist2 )
!             =======================  i  i  i     o
      else
         call getDist2FromRegion2dPoly( x, z, ir, dist2 )
!             ========================  i  i  i     o
      end if

   end subroutine getDist2FromRegion2d
!
!=======================================================================
!
   subroutine getDist2RectRegion2d( ic, ir, cinside, dist2 )
!             ====================  i   i      i       o
!
!     Determine the square of the minimum distance between any point inside cell ic and
!     region ir
!
!  Arguments:
!
      integer(INT8), intent(in) :: ic           ! cell index
      integer(INT4), intent(in) :: ir           ! index of region
      logical, intent(out) :: cinside           ! flag indicating rectangle center is inside
                                                !    region (.false. may also indicate
      real(REAL8), intent(out) :: dist2         ! square of distance betweeen ractangle and
                                                !    region
!
!  Local Variables:
!
!     None
!
!  Start of Executable Code:
!  ------------------------
!
      if( geom%usehyperoctree ) then
#ifdef GTS
         call getDist2RectRegion2dSDFOct( ic, ir, cinside, dist2 )
!             ==========================  i   i      i       o
#else
         error stop 'getDist2RectRegion2d: getDist2RectRegion2dSDFOct requires separate license'
#endif
         cinside = .false.
      else if( geom%usesdf ) then
         call getDist2RectRegion2dSDF( ic, ir, cinside, dist2 )
!             =======================  i   i      i       o
      else
         call getDist2RectRegion2dPoly( ic, ir, cinside, dist2 )
!             ========================  i   i      i       o
      end if

   end subroutine getDist2RectRegion2d
!
!=======================================================================
!
   subroutine isOutsideRegion2d( x, z, ir, outside )
!             =================  i  i  i      o
!
!     Determine if the point [x,z] is located outside region ir
!
!     "outside" here refers to the position being outside the region polygon. This is not
!     exactly the same as sdf < 0 if the region is defined by a polygon.
!
!  Arguments:
!
      real(REAL8), intent(in) :: x               ! x coordinate of point
      real(REAL8), intent(in) :: z               ! z coordinate of point
      integer(INT4), intent(in) :: ir           ! index of region
      logical, intent(out) :: outside           ! flag indicating point is outside region
!
!  Local Variables:
!
!     None
!
!  Start of Executable Code:
!  ------------------------
!
      if( geom%usehyperoctree ) then
         error stop 'isOutsideRegion2dSDFOct not implemented'
      else if( geom%usesdf ) then
         call isOutsideRegion2dSDF( x, z, ir, outside )
!             ====================  i  i  i      o
      else
         call isOutsideRegion2dPoly( x, z, ir, outside )
!             =====================  i  i  i     o
      end if

   end subroutine isOutsideRegion2d
!
!=======================================================================
!
end module Target_Region2d_mod
