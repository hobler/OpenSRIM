module Cinterface_LSHyperoctree2d_mod
!
!  Provides procedures to create and query a 2D hyperoctree holding signed distance functions
!
!  The procedures are implemented as type-bound procedures of LSHyperoctree2d_type
!
!  Used Modules:
!  ============
!
    use Const_mod,             only: INT4, INT8, REAL8,  &
                                     BIG

    use, intrinsic :: iso_c_binding
!
!  Default Type and Access:
!  =======================
!
   implicit none
   private
!
!  Public Data Types
!  =================
!
   public :: LSHyperoctree2d_type
!            ====================
   type LSHyperoctree2d_type
      private
      type(c_ptr) :: objectpointer
   contains
      ! final :: delete_lshyperoctree2d
      ! procedure :: delete => delete_lshyperoctree2d_polymorph ! Apparently this is for gfortran

      procedure :: getRegionIndex  &                           ! get index of nearest region
                   => lshyperoctree2d_getRegionIndex           !
      procedure :: getDistSurface  &                           ! get signed distance from
                   => lshyperoctree2d_getDistSurface           !    surface
      procedure :: getBoundaryLocationMin  &                   ! get lower boundary positions
                   => lshyperoctree2d_getBoundaryLocationMin   !
      procedure :: getBoundaryLocationMax  &                   ! get upper boundary positions
                   => lshyperoctree2d_getBoundaryLocationMax   !
      procedure :: getBoundaryConditionMin &                   ! get boundary type of lower
                   => lshyperoctree2d_getBoundaryConditionMin  !    boundary
      procedure :: getBoundaryConditionMax  &                  ! get boundary type of upper
                   => lshyperoctree2d_getBoundaryConditionMax  !    boundary
      procedure :: getGridDelta &                              ! get level 0 grid spacings
                   => lshyperoctree2d_getGridDelta             !
      procedure :: getNumberOfLeafNodes  &                     ! get number of cells
                   => lshyperoctree2d_getNumberOfLeafNodes     !
      procedure :: getLeafNodeIdAtCoords  &                    ! get cell index
                   => lshyperoctree2d_getLeafNodeIdAtCoords    !
      procedure :: getCornersSDFValuesAtLeafNodeId  &          ! get signed distance of cell
                   => lshyperoctree2d_getCornersSDFValuesAtLeafNodeId ! corners from region
   end type

   interface LSHyperoctree2d_type
      procedure create_lshyperoctree2d          ! create hyperoctree
   end interface
!
!  Private Procedures:
!  ==================
!
   private :: getNegativeUndefinedDouble        ! get value indicating undefined signed distance
!             ==========================        !    (outside target/region)
   private :: getPositiveUndefinedDouble        ! get value indicating undefined signed distance
!             ==========================        !    (inside target/region)
!
!  Private Data:
!  ============
!
   real(c_double), private :: BIG_c             ! return value of getPositiveUndefinedDouble
   real(c_double), private :: NEGBIG_c          ! return value of getNegativeUndefinedDouble
!
!  C Function Interfaces:
!  =====================
!
   interface
      ! Constructor
      function create_lshyperoctree2d_c( inputVTUFiles, numberOfRegions, denseDistance )  &
                                              bind(C, name="create_IMSILLSHyperoctree2DWrapper")
         use, intrinsic :: iso_c_binding
         implicit none
         type(c_ptr) :: create_lshyperoctree2d_c
         type (c_ptr), dimension(*), intent(in) ::  inputVTUFiles
         integer(c_int32_t), value, intent(in) :: numberOfRegions
         real(c_double), value, intent(in) :: denseDistance
      end function

      ! Destructor
      subroutine delete_lshyperoctree2d_c( hyperoctree )  &
                                              bind(C, name="delete_IMSILLSHyperoctree2DWrapper")
         use, intrinsic :: iso_c_binding
         implicit none
         type(c_ptr), value :: hyperoctree
      end subroutine

      ! Function getRegionIndex
      subroutine lshyperoctree2d_getRegionIndex_c( hyperoctree, x, y, regionIndex, isOutside ) &
                                      bind(C, name="IMSILLSHyperoctree2DWrapper_getRegionIndex")
         use, intrinsic :: iso_c_binding
         implicit none
         type(c_ptr), value, intent(in) :: hyperoctree
         real(c_double), value, intent(in) :: x
         real(c_double), value, intent(in) :: y
         type(c_ptr), value :: regionIndex
         type(c_ptr), value :: isOutside
      end subroutine

      ! Function getDistanceToSurfaceEstimate
      subroutine lshyperoctree2d_getDistanceToSurfaceEstimate_c(  &
            hyperoctree, x, y, isOutsideHyperoctree, signedDistance, isEstimate)  &
                        bind(C, name="IMSILLSHyperoctree2DWrapper_getDistanceToSurfaceEstimate")
         use, intrinsic :: iso_c_binding
         implicit none
         type(c_ptr), value, intent(in) :: hyperoctree
         real(c_double), value, intent(in) :: x
         real(c_double), value, intent(in) :: y
         type(c_ptr), value :: isOutsideHyperoctree
         type(c_ptr), value :: signedDistance
         type(c_ptr), value :: isEstimate
      end subroutine

      ! Function getDistanceToSurfaceAndDirectionEstimate
      subroutine lshyperoctree2d_getDistanceToSurfaceAndDirectionEstimate_c(  &
            hyperoctree, x, y, isOutsideHyperoctree, signedDistance, xDir, yDir, isEstimate)  &
            bind(C, name="IMSILLSHyperoctree2DWrapper_getDistanceToSurfaceAndDirectionEstimate")
         use, intrinsic :: iso_c_binding
         implicit none
         type(c_ptr), value, intent(in) :: hyperoctree
         real(c_double), value, intent(in) :: x
         real(c_double), value, intent(in) :: y
         type(c_ptr), value :: isOutsideHyperoctree
         type(c_ptr), value :: signedDistance
         type(c_ptr), value :: xDir
         type(c_ptr), value :: yDir
         type(c_ptr), value :: isEstimate
      end subroutine

      ! Function getBoundaryLocationMin
      subroutine lshyperoctree2d_getBoundaryLocationMin_c( hyperoctree, x, y )  &
                              bind(C, name="IMSILLSHyperoctree2DWrapper_getBoundaryLocationMin")
         use, intrinsic :: iso_c_binding
         implicit none
         type(c_ptr), value, intent(in) :: hyperoctree
         type(c_ptr), value :: x
         type(c_ptr), value :: y
      end subroutine

      ! Function getBoundaryLocationMax
      subroutine lshyperoctree2d_getBoundaryLocationMax_c( hyperoctree, x, y )  &
                              bind(C, name="IMSILLSHyperoctree2DWrapper_getBoundaryLocationMax")
         use, intrinsic :: iso_c_binding
         implicit none
         type(c_ptr), value, intent(in) :: hyperoctree
         type(c_ptr), value :: x
         type(c_ptr), value :: y
      end subroutine

      ! Function getBoundaryConditionMin
      subroutine lshyperoctree2d_getBoundaryConditionMin_c( hyperoctree, x, y )  &
                             bind(C, name="IMSILLSHyperoctree2DWrapper_getBoundaryConditionMin")
         use, intrinsic :: iso_c_binding
         implicit none
         type(c_ptr), value, intent(in) :: hyperoctree
         type(c_ptr), value :: x
         type(c_ptr), value :: y
      end subroutine

      ! Function getBoundaryConditionMax
      subroutine lshyperoctree2d_getBoundaryConditionMax_c( hyperoctree, x, y )  &
                             bind(C, name="IMSILLSHyperoctree2DWrapper_getBoundaryConditionMax")
         use, intrinsic :: iso_c_binding
         implicit none
         type(c_ptr), value, intent(in) :: hyperoctree
         type(c_ptr), value :: x
         type(c_ptr), value :: y
      end subroutine
        
      ! Function getGridDelta
      subroutine lshyperoctree2d_getGridDelta_c( hyperoctree, x, y )  &
                                        bind(C, name="IMSILLSHyperoctree2DWrapper_getGridDelta")
         use, intrinsic :: iso_c_binding
         implicit none
         type(c_ptr), value, intent(in) :: hyperoctree
         type(c_ptr), value :: x
         type(c_ptr), value :: y
      end subroutine

      ! Function getNumberOfLeafNodes
      function lshyperoctree2d_getNumberOfLeafNodes_c( hyperoctree )  &
                                bind(C, name="IMSILLSHyperoctree2DWrapper_getNumberOfLeafNodes")
         use, intrinsic :: iso_c_binding
         implicit none
         integer(c_int64_t) :: lshyperoctree2d_getNumberOfLeafNodes_c
         type(c_ptr), value, intent(in) :: hyperoctree
      end function

      ! Function getLeafNodeIdAtCoords
      function lshyperoctree2d_getLeafNodeIdAtCoords_c( hyperoctree, x, y )  &
                               bind(C, name="IMSILLSHyperoctree2DWrapper_getLeafNodeIdAtCoords")
         use, intrinsic :: iso_c_binding
         implicit none
         integer(c_int64_t) :: lshyperoctree2d_getLeafNodeIdAtCoords_c
         type(c_ptr), value, intent(in) :: hyperoctree
         real(c_double), value, intent(in) :: x
         real(c_double), value, intent(in) :: y
      end function

      ! Function getCornersSDFValuesAtLeafNodeId
      subroutine lshyperoctree2d_getCornersSDFValuesAtLeafNodeId_c(  &
            hyperoctree, leafNodeId, regionIndex, level, xLocations, yLocations, sdfValues)  &
                     bind(C, name="IMSILLSHyperoctree2DWrapper_getCornersSDFValuesAtLeafNodeId")
         use, intrinsic :: iso_c_binding
         implicit none
         type(c_ptr), value, intent(in) :: hyperoctree
         integer(c_int64_t), value, intent(in) :: leafNodeId
         integer(c_int32_t), value, intent(in) :: regionIndex
         type(c_ptr), value :: level
         type(c_ptr), value :: xLocations
         type(c_ptr), value :: yLocations
         type(c_ptr), value :: sdfValues
      end subroutine

      ! Function getNegativeUndefinedDouble
      function lshyperoctree2d_getNegativeUndefinedDouble_c( )  &
                          bind(C, name="IMSILLSHyperoctree2DWrapper_getNegativeUndefinedDouble")
         use, intrinsic :: iso_c_binding
         implicit none
         real(c_double) :: lshyperoctree2d_getNegativeUndefinedDouble_c
      end function

      ! Function getPositiveUndefinedDouble
      function lshyperoctree2d_getPositiveUndefinedDouble_c( )  &
                          bind(C, name="IMSILLSHyperoctree2DWrapper_getPositiveUndefinedDouble")
         use, intrinsic :: iso_c_binding
         implicit none
         real(c_double) :: lshyperoctree2d_getPositiveUndefinedDouble_c
      end function
   end interface

contains
!
!=======================================================================
!
   function create_lshyperoctree2d( inputVTUFiles, numberOfRegions, denseDistance )
!           ======================        i               i               i
!
!  Create the 2D hyperoctree
!
!  TODO: Cleanup unused inputVTUFiles and/or numberOfRegions
!
      type(LSHyperoctree2d_type) :: create_lshyperoctree2d
      character(80), allocatable, intent(in) :: inputVTUFiles(:)  ! paths to the input files
                                                                  !    describing the geometry
      integer(INT4), intent(in) :: numberOfRegions                ! number of material regions
                                                                  !    = number of input files
      real(REAL8), intent(in) :: denseDistance                    ! distance around surfaces
                                                                  !    to be densely filled
!
!  Local variables
!
      integer(INT4) :: i
      integer(INT4) :: filenameLength

      ! Turn the fortran array of "strings" into the corresponding c-style format
      character(len=81, kind=c_char), allocatable, target :: cFileNames(:)
      type(c_ptr), allocatable :: cFileNamesArrayPointer(:)
!
!  Start of Executable Code:
!  ------------------------
!
      BIG_c = getPositiveUndefinedDouble( )
!             ==========================
      NEGBIG_c = getNegativeUndefinedDouble( )
!                ==========================

      allocate( cFileNames( numberOfRegions ) ) ! We drop SDF encompassing all material regions

      ! Convert names of all single-material regions
      do i = 1, numberOfRegions
         cFileNames(i) = inputVTUFiles(i)
         filenameLength = len(trim(cFileNames(i)))
         cFileNames(i)(filenameLength+1:filenameLength+1) = C_NULL_CHAR
      end do

      allocate( cFileNamesArrayPointer( numberOfRegions ) )
      do i = 1, numberOfRegions
         cFileNamesArrayPointer(i) = c_loc( cFileNames(i) )
      end do

      create_lshyperoctree2d%objectpointer =  &
            create_lshyperoctree2d_c( cFileNamesArrayPointer, numberOfRegions, denseDistance )
!           ========================

      ! Deallocate temporary arrays
      deallocate( cFileNamesArrayPointer )
      deallocate( cFileNames )

   end function create_lshyperoctree2d
!
!=======================================================================
!
   subroutine delete_lshyperoctree2d( this )
!             ======================   io
!
!  Delete the hyperoctree
!
!  This subroutine is currently unused, due to the invailablity of smart pointers.
!  Using it would lead to an immediate deletion of the underlying c-object
!
      type(LSHyperoctree2d_type) :: this
!
!  Start of Executable Code:
!  ------------------------
!
      call delete_lshyperoctree2d_c( this%objectpointer )
!          ========================

      this%objectpointer = c_null_ptr           ! Protect against double delete

   end subroutine delete_lshyperoctree2d
!
!=======================================================================
!
   subroutine delete_lshyperoctree2d_polymorph( this )
!             ================================
!
!  Delete the hyperoctree
!
!  This subroutine is currently unused, due to the invailablity of smart pointers.
!  Using it would lead to an immediate deletion of the underlying c-object
!
      class(LSHyperoctree2d_type) :: this
!
!  Start of Executable Code:
!  ------------------------
!
      call delete_lshyperoctree2d_c( this%objectpointer )
!          ========================

      this%objectpointer = c_null_ptr           ! Protect against double delete

   end subroutine delete_lshyperoctree2d_polymorph
!
!=======================================================================
!
   subroutine lshyperoctree2d_getRegionIndex( this, x, z, regionIndex, isOutside )
!             ==============================   i    i  i       o          (o)
!
!  Get the index of the region (x,z) is in or closest to
!
      class(LSHyperoctree2d_type), intent(in) :: this
      real(REAL8), intent(in) :: x                 ! x coordinate of point to be queried
      real(REAL8), intent(in) :: z                 ! z coordinate of point to be queried
      integer(INT4), intent(out) :: regionIndex    ! index of nearest region
      logical, optional, intent(out) :: isOutside  ! flag indicating whether (x,y) is outside
                                                   !    any region
!
!  Local Variables:
!
      integer(c_int32_t), target :: local_regionIndex
      integer(c_int8_t), target :: local_isOutside
!
!  Start of Executable Code:
!  ------------------------
!
!  Call c function using local variables
!
      call lshyperoctree2d_getRegionIndex_c(  &
!          ================================
               this%objectpointer, x, z, c_loc(local_regionIndex), c_loc(local_isOutside) )

      regionIndex = local_regionIndex + 1    ! +1 to get 1-based from 0-based index
!
!  Copy over optional variables
!
      if( present( isOutside ) ) then
         if (local_isOutside == 0) then
            isOutside = .false.
         else
            isOutside = .true.
         end if
      end if

   end subroutine lshyperoctree2d_getRegionIndex
!
!=======================================================================
!
   subroutine lshyperoctree2d_getDistSurface( this, x, z, signedDistance, xDir, zDir,  &
!             ==============================   i    i  i        o         (o)   (o)
                                              isOutside, isEstimate )
!                                                (o)        (o)
!
!  Get the signed distance of the point (x,z) from the surface
!
!  The signed distance is positive inside the target and negative outside
!
!  The returned value of signedDistance should be accurate within denseDistance (see
!  create_lshyperoctree2d) from the surface. Positive values larger than denseDistance are
!  BIG. Negative signed distances less than -denseDistance are returned as -BIG,
!  indicating that the point is more than denseDistance outside the target.
!
!  xDir and zDir do not contain valid values if abs(signedDistance)=BIG
!
      class(LSHyperoctree2d_type), intent(in) :: this
      real(REAL8), intent(in) :: x                       ! x-coordinate of point to be queried
      real(REAL8), intent(in) :: z                       ! z-coordinate of point to be queried
      real(REAL8), target, intent(out) :: signedDistance ! signed distance from the surface
      real(REAL8), optional, intent(out) :: xDir         ! x-component of direction towards/
                                                         !    from surface
      real(REAL8), optional, intent(out) :: zDir         ! z-component of direction towards/
                                                         !    from surface
      logical, optional, intent(out) :: isOutside        ! flag indicating point is outside the
                                                         !    octree
      logical, optional, intent(out) :: isEstimate       ! flag indicating whether sdf value is
                                                         !    only lower limit
!
!  Local Variables:
!
      integer(c_int8_t), target :: local_isOutside
      real(REAL8), target :: local_xDir
      real(REAL8), target :: local_zDir
      integer(c_int8_t), target :: local_isEstimate
!
!  Start of Executable Code:
!  ------------------------
!
      if ( .not. ( present( xDir ) .or. present( zDir ) ) ) then

         ! Call c function using local variables
         call lshyperoctree2d_getDistanceToSurfaceEstimate_c(  &
!             ==============================================
                  this%objectpointer, x, z, c_loc(local_isOutside), c_loc(signedDistance),  &
                  c_loc(local_isEstimate) )

         if( signedDistance == BIG_c ) then
            signedDistance = -BIG
         else if( signedDistance == NEGBIG_c ) then
            signedDistance = BIG
         else
            signedDistance = -signedDistance    ! IMPORTANT: Invert the sign of the sdf value
         end if

      else

         ! Call c function using local variables
         call lshyperoctree2d_getDistanceToSurfaceAndDirectionEstimate_c(  &
!             ==========================================================
                  this%objectpointer, x, z, c_loc(local_isOutside), c_loc(signedDistance),  &
                  c_loc(local_xDir), c_loc(local_zDir), c_loc(local_isEstimate) )

         if( signedDistance == BIG_c ) then
            signedDistance = -BIG
         else if( signedDistance == NEGBIG_c ) then
            signedDistance = BIG
         else
            signedDistance = -signedDistance    ! IMPORTANT: Invert the sign of the sdf value
         end if
!
!  Copy over optional variables
!
         if( present( xDir ) ) then
            xDir = local_xDir
         end if

         if( present( zDir ) ) then
            zDir = local_zDir
         end if

      end if
!
!  Copy over optional variables
!
      if ( present( isOutside ) ) then
         if ( local_isOutside == 0 ) then
            isOutside = .false.
         else
            isOutside = .true.
         end if
      end if

      if( present( isEstimate ) ) then
         if ( local_isEstimate == 0 ) then
            isEstimate = .false.
         else
            isEstimate = .true.
         end if
      end if

   end subroutine lshyperoctree2d_getDistSurface
!
!=======================================================================
!
   subroutine lshyperoctree2d_getBoundaryLocationMin( this, xmin, zmin )
!             ======================================   i     o     o
!
!  Get the coordinates of the lower boundaries in x and z direction
!
      class(LSHyperoctree2d_type), intent(in) :: this
      real(REAL8), target, intent(out) :: xmin     ! coordinate of lower boundary in x direction
      real(REAL8), target, intent(out) :: zmin     ! coordinate of lower boundary in z direction
!
!  Start of Executable Code:
!  ------------------------
!
      ! Call c function
      call lshyperoctree2d_getBoundaryLocationMin_c(  &
!          ========================================
               this%objectpointer, c_loc(xmin), c_loc(zmin) )
        
   end subroutine lshyperoctree2d_getBoundaryLocationMin
!
!=======================================================================
!
   subroutine lshyperoctree2d_getBoundaryLocationMax( this, xmax, zmax )
!             ======================================   i     o     o
!
!  Get the coordinates of the upper boundaries in x and z direction
!
      class(LSHyperoctree2d_type), intent(in) :: this
      real(REAL8), target, intent(out) :: xmax     ! coordinate of upper boundary in x direction
      real(REAL8), target, intent(out) :: zmax     ! coordinate of upper boundary in z direction
!
!  Start of Executable Code:
!  ------------------------
!
      ! Call c function
      call lshyperoctree2d_getBoundaryLocationMax_c(  &
!          ========================================
               this%objectpointer, c_loc(xmax), c_loc(zmax) )

   end subroutine lshyperoctree2d_getBoundaryLocationMax
!
!=======================================================================
!
   subroutine lshyperoctree2d_getBoundaryConditionMin( this, bndtypex, bndtypez )
!             =======================================   i       o         o
!
!  Get the boundary types of the lower boundaries in x and z direction
!
!  The returned values are integers with the possible values and their meaning
!
!     0 - unknown
!     1 - periodic
!     2 - reflective
!     3 - extend
!     4 - infinite (free)
!
      class(LSHyperoctree2d_type), intent(in) :: this
      integer(INT4), target, intent(out) :: bndtypex  ! lower x boundary type
      integer(INT4), target, intent(out) :: bndtypez  ! lower z boundary type
!
!  Start of Executable Code:
!  ------------------------
!
      ! Call c function
      call lshyperoctree2d_getBoundaryConditionMin_c(  &
!          =========================================
               this%objectpointer, c_loc(bndtypex), c_loc(bndtypez) )

   end subroutine lshyperoctree2d_getBoundaryConditionMin
!
!=======================================================================
!
   subroutine lshyperoctree2d_getBoundaryConditionMax( this, bndtypex, bndtypez )
!             =======================================   i       o         o
!
!  Get the boundary types of the upper boundaries in x and z direction
!
!  The returned values are integers with the possible values and their meaning
!
!     0 - unknown
!     1 - periodic
!     2 - reflective
!     3 - extend
!     4 - infinite (free)
!
      class(LSHyperoctree2d_type), intent(in) :: this
      integer(INT4), target, intent(out) :: bndtypex  ! upper x boundary type
      integer(INT4), target, intent(out) :: bndtypez  ! upper z boundary type
!
!  Start of Executable Code:
!  ------------------------
!
      ! Call c function
      call lshyperoctree2d_getBoundaryConditionMax_c(  &
!          =========================================
               this%objectpointer, c_loc(bndtypex), c_loc(bndtypez) )

   end subroutine lshyperoctree2d_getBoundaryConditionMax
!
!=======================================================================
!
   subroutine lshyperoctree2d_getGridDelta( this, dx, dz)
!             ============================   i    o   o
!
!  Get the widths of level 0  cells in x and z direction
!
      class(LSHyperoctree2d_type), intent(in) :: this
      real(REAL8), target, intent(out) :: dx       ! grid delta in x direction
      real(REAL8), target, intent(out) :: dz       ! grid delta in z direction
!
!  Start of Executable Code:
!  ------------------------
!
      ! Call c function
      call lshyperoctree2d_getGridDelta_c( this%objectpointer, c_loc(dx), c_loc(dz) )
!          ==============================

   end subroutine lshyperoctree2d_getGridDelta
!
!=======================================================================
!
   integer(c_int64_t) function lshyperoctree2d_getNumberOfLeafNodes( this )
!                              ====================================   i
!
!  Get the number of cells of the octree
!
      class(LSHyperoctree2d_type), intent(in) :: this
!
!  Start of Executable Code:
!  ------------------------
!
      ! Call c function
      lshyperoctree2d_getNumberOfLeafNodes =  &
            lshyperoctree2d_getNumberOfLeafNodes_c( this%objectpointer )
!           ======================================

   end function lshyperoctree2d_getNumberOfLeafNodes
!
!=======================================================================
!
   integer(c_int64_t) function lshyperoctree2d_getLeafNodeIdAtCoords( this, x, z )
!                              =====================================   i    i  i
!
!  Get the cell index at point (x,z)
!
      class(LSHyperoctree2d_type), intent(in) :: this
      real(REAL8), intent(in) :: x              ! x-coordinate of point to be queried
      real(REAL8), intent(in) :: z              ! z-coordinate of point to be queried
!
!  Start of Executable Code:
!  ------------------------
!
      ! Call c function plus 0-based => 1-based
      lshyperoctree2d_getLeafNodeIdAtCoords =  &
            lshyperoctree2d_getLeafNodeIdAtCoords_c( this%objectpointer, x, z ) + 1
!           =======================================

   end function lshyperoctree2d_getLeafNodeIdAtCoords
!
!=======================================================================
!
   subroutine lshyperoctree2d_getCornersSDFValuesAtLeafNodeId( this, leafNodeId, regionIndex,  &
!             ===============================================   i         i           i
                                                               level, xLocations, zLocations,  &
!                                                                o         o           o
                                                               sdfValues )
!                                                                  o
!  Get the coordinates and signed distances from a region at the corners of a cell
!
      class(LSHyperoctree2d_type), intent(in) :: this
      integer(INT8), intent(in) :: leafNodeId            ! cell index
      integer(INT4), intent(in) :: regionIndex           ! region index
      integer(INT4), intent(out) :: level                ! level of the queried node
      real(REAL8), target, intent(out) :: xLocations(4)  ! x coordinates of the corner points
      real(REAL8), target, intent(out) :: zLocations(4)  ! z coordinates of the corner points
      real(REAL8), target, intent(out) :: sdfValues(4)   ! SDF values at the corner points
!
!  Local Variables:
!
      integer(c_int64_t) :: local_leafNodeId
      integer(c_int32_t) :: local_regionIndex
      integer(c_int16_t), target :: local_level
      integer(INT4) :: i
!
!  Start of Executable Code:
!  ------------------------
!
      ! 1-based => 0 based plus type conversion
      local_leafNodeId = leafNodeId - 1
      local_regionIndex = regionIndex - 1
        
      ! Call c function
      call lshyperoctree2d_getCornersSDFValuesAtLeafNodeId_c(  &
!          =================================================
               this%objectpointer, local_leafNodeId, local_regionIndex, c_loc(local_level),  &
               c_loc(xLocations), c_loc(zLocations), c_loc(sdfValues) )

      ! Flip sign of SDF values
      do i = 1, 4
         if( sdfValues(i) == BIG_c ) then
            sdfValues(i) = -BIG
         else if( sdfValues(i) == NEGBIG_c ) then
            sdfValues(i) = BIG
         else
            sdfValues(i) = -sdfValues(i)
         end if
      end do

      ! Copy over non-optional variables
      level = local_level

   end subroutine lshyperoctree2d_getCornersSDFValuesAtLeafNodeId
!
!=======================================================================
!
   real(c_double) function getNegativeUndefinedDouble( )
!                          ==========================
!
!  Get the negative number that indicates undefined values
!
!  Start of Executable Code:
!  ------------------------
!
      ! Call c function
      getNegativeUndefinedDouble = lshyperoctree2d_getNegativeUndefinedDouble_c( )
!                                  ============================================

   end function getNegativeUndefinedDouble
!
!=======================================================================
!
   real(c_double) function getPositiveUndefinedDouble( )
!                          ==========================
!
!  Get the positive number that indicates undefined values
!
!  Start of Executable Code:
!  ------------------------
!
      ! Call c function
      getPositiveUndefinedDouble = lshyperoctree2d_getPositiveUndefinedDouble_c( )
!                                  ============================================

   end function getPositiveUndefinedDouble
!
!=======================================================================
!
end module Cinterface_LSHyperoctree2d_mod
