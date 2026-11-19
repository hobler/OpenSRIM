module Target_Grid2d_mod
!
!  Interface layer for the Target_Grid2d_*_mod
!
!  Used Modules:
!  ============
!
   use Const_mod,                only: INT4, REAL8,  &
                                       ZERO, QUART, HALF, ONE, TWO, SIX,  &
                                       BIG, SMALL,  &
                                       ENUM

   use Param_Geom_mod,           only: geom

   use Types_mod,                only: Proj_type
   use Geom_Polygon_Old_mod,     only: getDistPolygon, getPolyArea
   use Target_Grid2d_Reg_mod,    only: getCellRegGrid2d, getOtherNearCellsRegGrid2d,  &
                                       initRegGrid2d, readRegGrid2d, writeRegGrid2d,  &
                                       getCellsRegGrid2d, initPostCellRegGrid2d,  &
                                       relaxRegGrid2d, getOvlpAreasRegGrid2d
   use Target_Grid2d_Curvi_mod,  only: initCurviGrid2d, getCellCurviGrid2d,  &
                                       getCellsCurviGrid2d, getOvlpAreasCurviGrid2d,  &
                                       writeCurviGrid2d, getOtherNearCellsCurviGrid2d,  &
                                       relaxCurviGrid2d
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
   public :: initGrid2d                      ! initialize the target grid
!            ==========
   public :: initPostCellGrid2d              ! grid initializations after cell initializations
!            ==================
   public :: getCellIndexGrid2d              ! determines the cell a given point is located in
!            ==================
   public :: getOtherNearCellIndicesGrid2d   ! determine cells other than nearest cell that
!            =============================   !    intersect with the free flight path cylinder
   public :: getAreasGrid2d                  ! obtain the cell areas
!            ==============
   public :: getOvlpAreasGrid2d              ! obtain the cell areas inside the polygon
!            ==================
   public :: relaxGrid2d                     ! relax the cell densities
!            ===========
   public :: writeGrid2d                     ! write target grid file
!            ===========
   public :: readGrid2d                      ! read target grid from file
!            ==========
!
!  Private Data:
!  ============
!
   type :: Cell_type
!          =========
      sequence
      real(REAL8), allocatable :: xv(:)         ! x-coordinates of cell vertices (closed
                                                !    polygon)
      real(REAL8), allocatable :: zv(:)         ! z-coordinates of cell vertices (closed
                                                !    polygon)
      integer(INT4), allocatable :: icnear(:)   ! indices of near nonempty cells, ordered
                                                !    according to increasing distmin2
      real(REAL8), allocatable :: prob(:)       ! approximate probability that icnear is
                                                !    nearest nonempty cell to a point inside
                                                !    the current cell
   end type Cell_type

   type(Cell_type), allocatable :: cell(:)
!
!==========================================================================================
!
contains
!
!==========================================================================================
!
   subroutine initGrid2d( ncell )
!             ==========    o
!
!     Initialize the grid
!
!  Arguments:
!
      integer(INT4), intent(out) :: ncell    ! number of cells
!
!  Local Variables:
!
!     None
!
!  Start of Executable Code:
!  ------------------------
!
      if( geom%grid == ENUM%grid%regular ) then
         call initRegGrid2d( )
!             =============
         call getCellsRegGrid2d( cell )
!             =================   o
      else if( geom%grid == ENUM%grid%curvilinear ) then
         call initCurviGrid2d( )
!             ===============
         call getCellsCurviGrid2d( cell )
!             ===================   o
      else
         error stop 'initGrid2d: invalid geom%grid'
      end if

      ncell = size( cell )

   end subroutine initGrid2d
!
!===================================================================================
!
   subroutine initPostCellGrid2d( ff, empty )
!             ==================  i     i
!
!     Initialize derived cell properties that do not depend on the grid
!
!  Arguments:
!
      real(REAL8), intent(in) :: ff(:)             ! fill fractions
      logical, intent(in) :: empty(:)              ! indicate empty cells
!
!  Local Variables:
!
!     None
!
!  Start of Executable Code:
!  ------------------------
!
      if( geom%grid == ENUM%grid%regular ) then
         call initPostCellRegGrid2d( ff, empty )
!             =====================  i     i
         call getCellsRegGrid2d( cell )
!             =================   o
      else if( geom%grid == ENUM%grid%curvilinear ) then
         ! nothing to do
      else
         error stop 'initPostCellGrid2d: invalid geom%grid'
      end if

   end subroutine initPostCellGrid2d
!
!==========================================================================================
!
   subroutine getCellIndexGrid2d( x, z, ic, outside_ )
!             ==================  i  i  o      o
!
!     Determine non-empty cell the point pos is located in or closest to.
!
!  Arguments:
!
      real(REAL8), intent(in) :: x              ! x-coordinate of position to be checked
      real(REAL8), intent(in) :: z              ! z-coordinate of position to be checked
      integer(INT4), intent(out) :: ic          ! 1D index of cell (region); if outside all
                                                !    cells, ic=0
      logical, optional, intent(out) :: outside_! flag indication whether pos is outside all
                                                !    cells with content
!
!  Local Variables:
!
     integer(INT4) :: in, icn, icnmin
     real(REAL8) :: dist2, dist2tmp
!
!  Start of Executable Code:
!  ------------------------
!
      if( geom%grid == ENUM%grid%regular ) then
         call getCellRegGrid2d( x, z, ic, outside_ )
!             ================  i  i  o      o
      else if( geom%grid == ENUM%grid%curvilinear ) then
         call getCellCurviGrid2d( x, z, ic, outside_ )
!             ==================  i  i  o      o
      else
         error stop 'getCellIndexGrid2d: invalid geom%grid'
      end if
!
!  If cell is empty, determine nearest non-empty cell
!
      if( cell(ic)%icnear(1) /= ic ) then

         if( present( outside_ ) )  outside_ = .true.

         dist2 = BIG

         do in = 1, size( cell(ic)%icnear(:) )  ! check only cells in icnear(:)
            icn = cell(ic)%icnear(in)

            call getDistPolygon( dist2tmp, x, z, cell(icn)%xv(:), cell(icn)%zv(:) )
!                ==============     o      i  i         i                i

            if( dist2tmp < dist2 ) then
               dist2 = dist2tmp
               icnmin = icn
            end if
         end do

         ic = icnmin                            ! replace current by nearest non-empty cell

      end if
!
   end subroutine getCellIndexGrid2d
!
!==========================================================================================
!
   subroutine getOtherNearCellIndicesGrid2d( proj, posc, ffp, pmax, ncn, icnear )
!             =============================   i     i     i    i     o     o
!
!     Determine cells other than ic that are within the free flight path (ffp) cylinder
!     (= cylinder with height ffp and radius pmax)
!     Only cells in the range (1:ncell) are returned.
!
!  Arguments:
!
      type(Proj_type), intent(in) :: proj       ! projectile
      real(REAL8), intent(in) :: posc(3)        ! position of end of free flight path
      real(REAL8), intent(in) :: ffp            ! free flight path (= height of cylinder)
      real(REAL8), intent(in) :: pmax           ! maximum impact parameter (= radius of cyl.)
      integer(INT4), intent(out) :: ncn         ! number of regions found
      integer(INT4), allocatable, intent(inout) :: icnear(:)
                                                ! index of regions found except proj%ir
!
!  Local Variables:
!
!     None
!
!  Start of Executable Code:
!  ------------------------
!
      if( geom%grid == ENUM%grid%regular ) then
         call getOtherNearCellsRegGrid2d( proj, posc, ffp, pmax, ncn, icnear )
!             ==========================   i     i     i    i     o     o
      else if( geom%grid == ENUM%grid%curvilinear ) then
         call getOtherNearCellsCurviGrid2d( proj, posc, ffp, pmax, ncn, icnear )
!             ============================   i     i     i    i     o     o
      else
         error stop 'getOtherNearCellIndicesGrid2d: invalid geom%grid'
      end if

   end subroutine getOtherNearCellIndicesGrid2d
!
!=======================================================================
!
   subroutine getAreasGrid2d( area )
!             ==============   o
!
!     Get the cell areas
!
!  Arguments:
!
      real(REAL8), intent(out) :: area(:)       ! cell area
!
!  Local Variables:
!
      integer(INT4) :: ic, nv
!
!  Start of Executable Code:
!  ------------------------
!
      do ic = 1, size( cell )
         nv = size( cell(ic)%xv ) - 1
         area(ic) = getPolyArea( cell(ic)%xv(1:nv), cell(ic)%zv(1:nv), nv )
!                   ===========          o                  o          o
      end do

   end subroutine getAreasGrid2d
!
!===================================================================================
!
   subroutine getOvlpAreasGrid2d( px, pz, area, ic_area)
!             ==================  i   i    o      o
!
!     Obtain the fractional cell area inside the given polygon assuming a Cartesian grid,
!     for further details see getCartCellAreas
!
!  Arguments:
!
      real(REAL8), intent(in) :: px(:)                   ! x-coordinates of polygon corners
      real(REAL8), intent(in) :: pz(:)                   ! z-coordinates of polygon corners
      real(REAL8), allocatable, intent(out) :: area(:)   ! area inside polygon
                                                         !    first index runs in x-direction
                                                         !    second index runs in z-direction
      integer(INT4), allocatable, intent(out) :: ic_area(:) ! cell indices corresponding to
                                                         !    the areas
!
!  Local Variables:
!
!     None
!
!  Start of Executable Code:
!  ------------------------
!
      if( geom%grid == ENUM%grid%regular ) then
         call getOvlpAreasRegGrid2d( px(:), pz(:), area, ic_area)
!             =====================   i      i      o       o
      else if( geom%grid == ENUM%grid%curvilinear ) then
         call getOvlpAreasCurviGrid2d( px(:), pz(:), area, ic_area)
!             =======================   i      i      o       o
      else
         error stop 'getOvlpAreasGrid2d: invalid geom%grid'
      end if

   end subroutine getOvlpAreasGrid2d
!
!=======================================================================
!
   subroutine relaxGrid2d( ff0, ff, dpa, stress, ovlp )
!             ===========   i   i    i     io     o
!
!     Redistribute the cell densities so that stress relaxes. This is done in two steps:
!     First the geometry of the grid is adjusted so that the cell densities (given by the
!     cell content and the cell size) are identical to the equilibrium densities for the
!     respective target composition. Then the cell contents are mapped to an equidistant
!     grid again.
!
!     This subroutine should not be called too often in order to avoid artificial diffusion,
!     but often enough to maintain realistic densities.
!
!  Arguments:
!
      type :: Ovlp_type
!             =========
         sequence
         integer(INT4), allocatable :: ic(:)    ! cell indices of overlapped cells
         real(REAL8), allocatable :: area(:)    ! overlap area
         real(REAL8) :: emptyarea               ! total overlap area with empty cells
      end type Ovlp_type

      real(REAL8), intent(in) :: ff0(:)         ! (ncell) fill fractions before BC time step
      real(REAL8), intent(in) :: ff(:)          ! (ncell) fill fractions after BC time step if
                                                !    there were no relaxation
      real(REAL8), intent(in) :: dpa(:)         ! (ncell) displacements per atom during BC step
      real(REAL8), intent(inout) :: stress(:,:) ! (13,ncell) pressure and deviatoric stress
                                                !    tensor elements
      type(Ovlp_type), intent(out) :: ovlp(:)   ! (ncell) overlap information
!
!  Local Variables:
!
!     None
!
!  Start of Executable Code:
!  ------------------------
!
      if( geom%grid == ENUM%grid%regular ) then
         call relaxRegGrid2d( ff(:), dpa(:), ovlp(:) )
!             ==============   i       i       o
      else if( geom%grid == ENUM%grid%curvilinear ) then
         call relaxCurviGrid2d( ff0(:), ff(:), dpa(:), stress(:,:), ovlp(:) )
!             ================    i      i       i        io         o
         call getCellsCurviGrid2d( cell )
!             ===================   o
      else
         error stop 'relaxGrid2d: invalid geom%grid'
      end if

   end subroutine relaxGrid2d
!
!=======================================================================
!
   subroutine writeGrid2d( nout, relaxed_ )
!             ===========   i       i
!
!     Write the contents of the target cells to the CELL file (same format as HIS file)
!
!  Arguments:
!
      integer(INT4), intent(in) :: nout         ! unit number of file
      logical, optional, intent(in) :: relaxed_ ! indicates relaxed grid to be written
!
!  Local Variables:
!
!     None
!
!  Start of Executable Code:
!  ------------------------
!
      if( geom%grid == ENUM%grid%regular ) then
         call writeRegGrid2d( nout, relaxed_ )
!             ==============   i       i
      else if( geom%grid == ENUM%grid%curvilinear ) then
         call writeCurviGrid2d( nout, relaxed_ )
!             ================   i       i
      else
         error stop 'writeGrid2d: invalid geom%grid'
      end if

   end subroutine writeGrid2d
!
!=======================================================================
!
   subroutine readGrid2d( nin, ncell, valid )
!             ==========   i     o      o
!
!     Read the contents of the target cells from the hisfile file (same format as HIS file)
!
!     UNTESTED CODE !
!
!  Arguments:
!
      integer(INT4), intent(in) :: nin          ! unit number of file
      integer(INT4), intent(out) :: ncell       ! number of cells
      logical, intent(out) :: valid             ! indicates valid grid has been read
!
!  Local Variables:
!
!     None
!
!  Start of Executable Code:
!  ------------------------
!
      if( geom%grid == ENUM%grid%regular ) then
         call readRegGrid2d( nin, valid )
!             =============   i     o
         call getCellsRegGrid2d( cell )
!             =================   o
      else
         error stop 'readGrid2d: invalid geom%grid'
      end if

      ncell = size( cell )

   end subroutine readGrid2d
!
!=======================================================================
!
end module Target_Grid2d_mod
