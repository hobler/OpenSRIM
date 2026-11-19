module Target_Cell_Reg_mod
!
!  Holds 1d/2d/3d regular grid and cell structure and provides subroutines operating on it.
!  A regular grid is maintained throughout the simulation.
!
!  A grid derived type variable with extent 1, 2, or 3 will be used to define the grid
!  points/lines/planes. The indices correspond to the directions as follows:
!
!  1d: grid(1)...z     2d: grid(1)...x     3d: grid(1)...x
!                          grid(2)...z         grid(2)...y
!                                              grid(3)...z
!
!  The grid lines are indexed ix(i) in the range (0:nx(i)). A cell is identified by the three
!  indices (ix(1),ix(2),ix(3)) in the range (1:nx(1),1:nx(2),1:nx(3)), where nx(i)=1 for
!  i>ndim. Where possible, a 1D index ic derived from ix(i) is used to address the cell
!  properties.
!
!  Note: avoid ldyn=t with lspher=t
!
!  Used Modules:
!  ============
!
!$ use omp_lib,               only: omp_get_thread_num
   use Const_mod,             only: INT4, REAL8,  &
                                    ZERO, ONE, HALF, TENTH, EPSMA, SMALL,  &
                                    ENUM
   use Array_mod,             only: alloc, check_allocation
   use String_mod,            only: str
   use Undefine_mod,          only: defined
   use Math_IntMath_mod,      only: index2to1, index1to2, index3to1, index1to3

   use Param_Setup_mod,       only: ndim, nelem, nelemIon, nionhis,  &
                                    usecell, cellfile, cellrec
   use Param_Geom_mod,        only: geom
   use Param_Elem_mod,        only: elem
   use Param_Dam_mod,         only: dam
   use Param_Out_mod,         only: out

   use Types_mod,             only: Proj_type

   use Target_Cell1d_Reg_mod, only: readAtomsFromCell1dReg, initCell1dReg, initHis1dReg,  &
                                    getNearestCell1dReg, getOtherNearCellIndices1dReg,  &
                                    fixNegDensCell1dReg, &
                                    diffuseCell1dReg, relaxCell1dReg, relaxCell1dRegFull,  &
                                    writeCell1dReg, writeDamCell1dReg
   use Target_Cell2d_Reg_mod, only: initCell2dReg, getNearestCell2dReg,  &
                                    getOtherNearCellIndices2dReg, fixNegDensCell2dReg, &
                                    relaxCell2dReg, writeCell2dReg

#ifdef AMAT
   use Cinterface_Diffusion2d_mod,  only: diffuseCell2dReg
#endif
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
   public :: initCellReg                     ! initialize the target
!            ===========
   public :: readAtomsFromCellReg            ! read atomic numbers from cell file
!            ====================
   public :: scoreCellReg                    ! score a count in the target cell
!            ============
   public :: getCellIndexReg                 ! determines the cell a given point is located in
!            ===============
   public :: getOtherNearCellIndicesReg      ! determine cells other than nearest cell that
!            ==========================      !    intersect with the free flight path cylinder
   public :: getProbCellReg                  ! get the probabilities to encounter each of the
!            ==============                  !    atoms in a given cell
   public :: getDensCellReg                  ! get the total atomic density in a given cell
!            ==============
   public :: getVibCellReg                   ! get the array of thermal vibration amplitudes
!            =============                   !    in a given cell
   public :: getXtalCellReg                  ! get the crystallinity flag of the cell
!            ==============
   public :: getEdCellReg                    ! get the displacement energy after most recent
!            ============                    !    target cell update
   public :: getEsurfCellReg                 ! get the surface binding energy after most recent
!            ===============                 !    target cell update
   public :: broadcastCellReg                ! broadcast cells from master thread to other
!            ================                !    threads
   public :: mergeCellReg                    ! merge the contents of cell_tmp into cell
!            ============
   public :: relaxCellReg                    ! relax the cell densities
!            ============
   public :: writeCellReg                    ! write target cell contents to output file
!            ============
   public :: writeDamCellReg                 ! write damage to output file
!            ===============
!
!  Private Procedures:
!  ==================
!
   private :: initDerivedCellReg             ! initialize derived quantities
!             ==================
!
!  Private Data:
!  ============
!
   type :: Grid_type
!          =========
      integer(INT4) :: nx                    ! number of cells in the respective direction
      real(REAL8) :: wx                      ! cell width in the respective direction
      real(REAL8), allocatable :: x(:)       ! (0:nx) positions of interfaces between
                                             !    cells
      real(REAL8) :: xmin1                   ! largest position to the left of x(0)
      real(REAL8) :: xmin2                   ! smallest position to the right of x(0)
      real(REAL8) :: xmax1                   ! largest position to the left of x(nx)
      real(REAL8) :: xmax2                   ! smallest position to the right of x(nx)
   end type Grid_type

   type(Grid_type), allocatable :: grid(:)   ! grid
   
   type :: Cell_type
!          =========
      integer(INT4) :: nc                       ! number of cells
      integer(INT4) :: nx(3)                    ! extents of cell array
      real(REAL8) :: size                       ! width/area/volume of a cell
      real(REAL8), allocatable :: dens(:,:)     ! (nelem,nc) density of atom ielem in cell
      real(REAL8), allocatable :: dens0(:,:)    ! (nelem,nc) density of atom ielem in cell when
                                                !    initDerivedCell1dReg was called (i.e, after
                                                !    last relaxation)
      real(REAL8), allocatable :: denstot(:)    ! (nc) total atom density in cell
      real(REAL8), allocatable :: prob(:,:)     ! (nelem,nc) probability to encounter atom
      real(REAL8), allocatable :: ff(:)         ! (nc) fill fraction after last relaxation
      real(REAL8), allocatable :: ffgeom(:)     ! (nc) geometric fill fraction
      real(REAL8), allocatable :: dam(:,:)      ! (nelemIn,nc) damage densities
      real(REAL8), allocatable :: vib(:,:)      ! (nelem,nc) thermal vibration amplitude
      logical, allocatable :: lxtal(:)          ! (nc) flag indicating crystallinity
      real(REAL8), allocatable :: ed(:,:)       ! (nelem,nc) displacement energy after
                                                !    last relaxation
      real(REAL8), allocatable :: esurf(:,:)    ! (nelem,nc) surface binding energy after
                                                !    last relaxation
      real(REAL8), allocatable :: edam(:)       ! (nc) damage energy per cm3 since last
                                                !    relaxation
      real(REAL8), allocatable :: pressure(:)   ! (nc) pressure after last relaxation
      real(REAL8), allocatable :: stress(:,:,:) ! (ndim,ndim,nc) stress tensor after last
                                                !    relaxation
   end type Cell_type

   type(Cell_type) :: cell                      ! cells
   type(Cell_type) :: cell_tmp                  ! same as cell but containing density and
                                                !    damage changes since last update
!$omp threadprivate( cell, cell_tmp )

   integer(INT4) :: nelemIn
!
!=======================================================================
!
contains
!
!=======================================================================
!
   subroutine initCellReg( )
!             ===========
!
!     Initialize the cell structure of the target from either a CELL* file or a region
!     definition
!
!     OpenMP: This subroutine should be called from the master thread
!
!  Arguments:
!
!     None
!
!  Local Variables:
!
      character(256) :: errmsg
      integer(INT4) :: istat
      logical, allocatable :: lxtal2d(:,:)
      real(REAL8), allocatable :: dens2d(:,:,:), ffgeom2d(:,:)
      real(REAL8), allocatable :: pressure2d(:,:), stress2d(:,:,:,:)
!
!  Start of Executable Code:
!  ------------------------
!
      nelemIn = 2 * nelem

      allocate( grid( ndim ), stat=istat, errmsg=errmsg )
      call check_allocation( 'grid', 'initCellReg', istat, errmsg )
!          ================    i           i          i      i
!
!  Initialize grid and cell except for cell%dam and derived cell properties
!
      if( ndim == 1 ) then

         call initCell1dReg( grid(1)%x,  &
!             =============      o
                             cell%dens, cell%ffgeom, cell%lxtal, cell%pressure, cell%stress )
!                                o           o           o              o            o

      else if( ndim == 2 ) then

         call initCell2dReg( grid(1)%x, grid(2)%x,  &
!             =============      o          o
                             dens2d, ffgeom2d, lxtal2d, pressure2d, stress2d )
!                              o        o         o         o          o
         grid(1)%nx = ubound( grid(1)%x, 1 )
         grid(2)%nx = ubound( grid(2)%x, 1 )
         cell%nc = product( grid(1:2)%nx )
         cell%dens = reshape( dens2d, [nelem, cell%nc] )
         cell%ffgeom = reshape( ffgeom2d, [cell%nc] )
         cell%lxtal = reshape( lxtal2d, [cell%nc] )
         if( dam%relaxmodel == ENUM%relaxmodel%viscoelastic ) then
            cell%pressure = reshape( pressure2d, [cell%nc] )
            cell%stress = reshape( stress2d, [2,2,cell%nc] )
         end if

      else

!         call initCell3dReg( grid(1)%x, grid(2)%x, grid(3)%x,  &
!!             =============      o
!                             cell%dens, cell%ffgeom, cell%lxtal, cell%pressure, cell%stress )
!!                                o           o           o              o             o
         error stop 'initCellReg: 3D case not yet implemented'

      end if
!
!  Initialize cell damage
!
      if( out%lhis ) then
         if( ndim == 1 ) then
            call initHis1dReg( grid(1)%x, cell%dam )
!                ============      i         io
         else
            error stop 'initHisReg: Damage not yet implemented for 2D and 3D'
         end if
      end if
!
!  Calculate derived quantities
!
      call initDerivedCellReg( )
!          ==================
!
!  Relax densities, warranting equilibrium densities
!
      if( dam%relaxmodel == ENUM%relaxmodel%full )  call relaxCellReg( init_=.true. )
!                                                        ============       i

   end subroutine initCellReg
!
!=======================================================================
!
   subroutine readAtomsFromCellReg( )
!             ====================
!
!     Read the atomic numbers from record cellrec of the CELL file cellfile.
!     The CELL file has the same format as the HIS file
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
         call readAtomsFromCell1dReg( )
!             ======================
      else if( ndim == 2 ) then
!         call readAtomsFromCell2dReg( )
!!             ======================
         error stop 'readAtomsFromCellReg: 2D not yet implemanted'
      else
         error stop 'readAtomsFromCellReg: 3D not yet implemanted'
      end if

   end subroutine readAtomsFromCellReg
!
!=======================================================================
!
   subroutine scoreCellReg( ielem, ic, w, ddami, ddamv, edam )
!             ============    i    i   i    i      i     i
!
!     Score a count in the 1d cell histogram
!
!     The cell contents may go below zero.
!
!  Arguments:
!
      integer(INT4), intent(in) :: ielem        ! index of atom species
      integer(INT4), intent(in) :: ic           ! cell index
      real(REAL8), intent(in) :: w              ! weight (can also be negative)
      real(REAL8), intent(in) :: ddami          ! interstitial increment
      real(REAL8), intent(in) :: ddamv          ! vacancy increment
      real(REAL8), intent(in) :: edam           ! damage energy (eV)
!
!  Local Variables:
!
      real(REAL8) :: ddens, ddam, dedam
!
!  Start of Executable Code:
!  ------------------------
!
!  Cell densities
!
      if( w /= ZERO ) then
         ddens = w * out%fachis(ndim) / cell%size
         if( nionhis > 0 ) then
            cell_tmp%dens(ielem,ic) = cell_tmp%dens(ielem,ic) + ddens
         else
            cell%dens(ielem,ic) = cell%dens(ielem,ic) + ddens
            cell%denstot(ic) = sum( cell%dens(:,ic) )
            cell%prob(:,ic) = cell%dens(:,ic) / cell%denstot(ic)
         end if
      end if
!
!  Damage densities
!
      if( out%lhis ) then
         if( ddami /= ZERO ) then
            ddam = ddami * out%fachis(ndim) / cell%size
            if( nionhis > 0 ) then
               cell_tmp%dam(ielem,ic) = cell_tmp%dam(ielem,ic) + ddam
            else
               cell%dam(ielem,ic) = cell%dam(ielem,ic) + ddam
               cell%dam(ielem,ic) = max( SMALL, cell%dam(ielem,ic) )
            end if
         end if
         if( ddamv /= ZERO ) then
            ddam = ddamv * out%fachis(ndim) / cell%size
            if( nionhis > 0 ) then
               cell_tmp%dam(ielem+nelem,ic) = cell_tmp%dam(ielem+nelem,ic) + ddam
            else
               cell%dam(ielem+nelem,ic) = cell%dam(ielem+nelem,ic) + ddam
               cell%dam(ielem+nelem,ic) = max( SMALL, cell%dam(ielem+nelem,ic) )
            end if
         end if
      end if
!
!  Damage energy per cell volume
!
      if( edam /= ZERO ) then
         dedam = edam * out%fachis(ndim) / cell%size
         if( nionhis > 0 ) then
            cell_tmp%edam(ic) = cell_tmp%edam(ic) + dedam
         else
            cell%edam(ic) = cell%edam(ic) + dedam
         end if
      end if
!
!  Do not update ed and esurf to avoid unreasonable values when cell gets empty!
!
   end subroutine scoreCellReg
!
!=======================================================================
!
   subroutine getCellIndexReg( pos, ic, outside_ )
!             ===============   i   o      o
!  
!     Determine the cell the point pos is located in or closest to.
!
!  Arguments:
!
      real(REAL8), intent(in) :: pos(3)         ! position to be checked
      integer(INT4), intent(out) :: ic          ! index of cell (region)
      logical, optional, intent(out) :: outside_! flag indicating whether pos is outside any
                                                !    region
!
!  Local Variables:                             !     / 3 1 1 \
!                                               ! j = | 0 3 2 |
      integer(INT4) :: ix(3), i                 !     \ 0 0 3 /
      integer(INT4) :: j(3,3) = reshape( [3, 0, 0, 1, 3, 0, 1, 2, 3], [3, 3] )
      real(REAL8) :: x(3)
!
!  Start of Executable Code:
!  ------------------------
!
      if( present( outside_ ) )  outside_ = .false.

      do i = 1, ndim                   ! don't apply bounds to iz because of integer overflow
         x(i) = max( grid(i)%xmin1, min( grid(i)%xmax2, pos(j(i,ndim)) ) )
         ix(i) = ceiling( ( x(i) - grid(i)%x(0) ) / grid(i)%wx )
         if( present( outside_ ) )  &
               outside_ = outside_ .or. ( ix(i) < 1 .or. ix(i) > grid(i)%nx )
         ix(i) = max( 1, min( grid(i)%nx, ix(i) ) )
      end do

      if( ndim == 1 ) then
         ic = ix(1)
      else if( ndim == 2 ) then
         call index2to1( ix(1), ix(2), grid(1)%nx, ic )
!             =========    i      i         i      o
      else
         call index3to1( ix(1), ix(2), ix(3), grid(1)%nx, grid(2)%nx, ic )
!             =========    i      i      i         i           i      o
      end if

   end subroutine getCellIndexReg
!
!=======================================================================
!
   subroutine getOtherNearCellIndicesReg( proj, posc, ffp, pmax, ncn, icnear )
!             ==========================   i     i     i    i     o     o
!
!     Determine cells other than proj%id that are within the free flight path cylinder
!     (= cylinder with height ffp and radius pmax)
!
!  Arguments:
!
      type(Proj_type), intent(in) :: proj       ! projectile
      real(REAL8), intent(in) :: posc(3)        ! position of end of free flight path
      real(REAL8), intent(in) :: ffp            ! free flight path (= height of cylinder)
      real(REAL8), intent(in) :: pmax           ! maximum impact parameter (= radius of cyl.)
      integer(INT4), intent(out) :: ncn         ! number of regions found
      integer(INT4), allocatable, intent(inout) :: icnear(:)   
                                                ! index of cells found except proj%id
!
!  Local Variables:
!
      logical :: allfar
      integer(INT4) :: ix(3), i
      integer(INT4) :: j(3,3) = reshape( [3, 0, 0, 1, 3, 0, 1, 2, 3], [3, 3] )
      real(REAL8) :: distother
!
!  Start of Executable Code:
!  ------------------------
!
!  First do a rough check whether other cells are far from the projectile
!
      call index1to3( proj%id, cell%nx(:), ix(:) )
!          =========     i         i         o

      allfar = .true.

      do i = 1, ndim
         if( proj%pos(j(i,ndim)) < grid(i)%x(1) ) then
            distother = grid(i)%x(1) - proj%pos(j(i,ndim))
         else if( proj%pos(j(i,ndim)) > grid(i)%x(grid(i)%nx-1) ) then
            distother = proj%pos(j(i,ndim)) - grid(i)%x(grid(i)%nx-1)
         else
            distother = min( grid(i)%x(ix(i)) - proj%pos(j(i,ndim)),  &
                             proj%pos(j(i,ndim)) - grid(i)%x(ix(i)-1) )
         end if

         allfar = allfar .and. ( distother**2 > ffp**2 + pmax**2 )
      end do

      if( allfar ) then
         ncn = 0
         return
      end if
!
!  Call dimension-specific subroutines
!
      if( ndim == 1 ) then
         call getOtherNearCellIndices1dReg( grid(1)%x, grid(1)%xmin2, grid(1)%xmax1,  &
!             ============================      i           i               i
                                            proj, posc, ffp, pmax, ncn, icnear )
!                                            i     i     i    i     o     o
      else if( ndim == 2 ) then
         call getOtherNearCellIndices2dReg( grid(1)%x, grid(1)%xmin2, grid(1)%xmax1,  &
!             ============================      i           i               i
                                            grid(2)%x, grid(2)%xmin2, grid(2)%xmax1,  &
!                                                 i           i               i
                                            proj, posc, ffp, pmax, ncn, icnear )
!                                            i     i     i    i     o     o
      else

         error stop 'getOtherNearCellIndicesReg: 3D not yet implemented'

      end if

   end subroutine getOtherNearCellIndicesReg
!
!=======================================================================
!
   subroutine getProbCellReg( ic, prob )
!             ==============  i    o
!
!     Get the probabilities to encounter each of the atoms in cell ic
!
!  Arguments:
!
      integer(INT4), intent(in) :: ic           ! cell index
      real(REAL8), intent(out) :: prob(:)       ! (nelem) probabilities
!
!  Local Variables:
!
!     None
!
!  Start of Executable Code:
!  ------------------------
!
      prob(:) = cell%prob(:,ic)

   end subroutine getProbCellReg
!
!=======================================================================
!
   subroutine getDensCellReg( ic, dens )
!             ==============  i    o
!
!     Get the total atomic density in cell ic
!
!     Zero is returned instead of a negative value.
!
!     For surface cells, denstot represents the average density of the cell. The actual
!     density outside the surface is zero. Inside the surface, the density has to be
!     divided by ffgeom. Checking on whether the point of interest is inside or outside
!     the target must be done outside this procedure.
!
!     For cells completely outside the surface, dens is set to zero.
!
!  Arguments:
!
      integer(INT4), intent(in) :: ic           ! cell index
      real(REAL8), intent(out) :: dens          ! density
!
!  Local Variables:
!
!     None
!
!  Start of Executable Code:
!  ------------------------
!
      if( cell%ffgeom(ic) == ZERO ) then
         dens = ZERO
      else
         dens = max( cell%denstot(ic), ZERO ) / cell%ffgeom(ic)
      end if

      dens = min( TENTH, dens )                 ! limit to 0.1 A^-3

   end subroutine getDensCellReg
!
!=======================================================================
!
   subroutine getVibCellReg( ic, vib )
!             =============  i    o
!  
!     Get the thermal vibration amplitude in cell ic
!
!  Arguments:
!
      integer(INT4), intent(in) :: ic           ! cell index
      real(REAL8), intent(out) :: vib(:)        ! thermal vibration amplitudes
!
!  Local Variables:
!
!     None
!
!  Start of Executable Code:
!  ------------------------
!
      vib(:) = cell%vib(:,ic)

   end subroutine getVibCellReg
!
!=======================================================================
!
   subroutine getXtalCellReg( ic, lxtal )
!             ==============  i     o
!
!     Get the crystallinity flag of the cell ic
!
!  Arguments:
!
      integer(INT4), intent(in) :: ic           ! cell index
      logical, intent(out) :: lxtal             ! crystallinity flag
!
!  Local Variables:
!
!     None
!
!  Start of Executable Code:
!  ------------------------
!
      lxtal = cell%lxtal(ic)

   end subroutine getXtalCellReg
!
!=======================================================================
!
   subroutine getEdCellReg( ic, ed )
!             ============  i   o
!
!     Get the displacement energy in cell ic
!
!  Arguments:
!
      integer(INT4), intent(in) :: ic           ! cell index
      real(REAL8), intent(out) :: ed(:)         ! displacement energy
!
!  Local Variables:
!
!     None
!
!  Start of Executable Code:
!  ------------------------
!
      ed(:) = cell%ed(:,ic)

   end subroutine getEdCellReg
!
!=======================================================================
!
   subroutine getEsurfCellReg( ic, esurf )
!             ===============  i     o
!
!     Get the surface binding energy in cell ic
!
!  Arguments:
!
      integer(INT4), intent(in) :: ic           ! cell index
      real(REAL8), intent(out) :: esurf(:)      ! surface binding energy
!
!  Local Variables:
!
!     None
!
!  Start of Executable Code:
!  ------------------------
!
      esurf(:) = cell%esurf(:,ic)

   end subroutine getEsurfCellReg
!
!=======================================================================
!
   subroutine broadcastCellReg( )
!             ================
!
!     Copy cells from master thread to other threads
!
!  Arguments:
!
!     None
!
!  Local Variables:
!
      type(Cell_type), save :: cell_shared   ! shared temporary cell
!
!  Start of Executable Code:
!  ------------------------
!
!$omp barrier
!$omp master
      cell_shared = cell
!$omp end master
!$omp barrier
      cell = cell_shared

      if( nionhis > 0 ) then
!$omp barrier
!$omp master
         cell_shared = cell_tmp
!$omp end master
!$omp barrier
         cell_tmp = cell_shared
      end if

   end subroutine broadcastCellReg
!
!=======================================================================
!
   subroutine mergeCellReg( )
!             ============
!
!     Merge contents of cell_tmp into cell and reinitialize cell_tmp. No changes are made to
!     the grid.
!
!     This subroutine must only be called if cell_tmp exists, i.e. when nionhis>0.
!
!     OpenMP: This subroutine must be called from within a parallel region. On exit, cell and
!     cell_tmp are updated on all threads.
!
!  Arguments:
!
!     None
!
!  Local Variables:
!
      integer(INT4) :: ithread, ic, nneg
      real(REAL8), pointer :: dens2d(:,:,:)
      type(Cell_type), target, save :: cell_shared  ! shared temporary cell
!
!  Start of Executable Code:
!  ------------------------
!
!  Sum up all cell_tmp into cell_shared
!
!$omp master
      cell_shared = cell_tmp
!$omp end master

!$omp barrier
!$omp critical( critical_cell1d_mergeCell1d )
      ithread = 1
!$    ithread = 1 + omp_get_thread_num( )
      if( ithread > 1 ) then
         cell_shared%dens(:,:) = cell_shared%dens(:,:) + cell_tmp%dens(:,:)
         if( out%lhis )  cell_shared%dam(:,:) = cell_shared%dam(:,:) + cell_tmp%dam(:,:)
         cell_shared%edam(:) = cell_shared%edam(:) + cell_tmp%edam(:)
      end if
!$omp end critical( critical_cell1d_mergeCell1d )
!$omp barrier
!
!  Add old cell densities (which are the same on all threads)
!
!$omp master
      cell_shared%dens(:,:) = cell_shared%dens(:,:) + cell%dens(:,:)
      if( out%lhis )  cell_shared%dam(:,:) = cell_shared%dam(:,:) + cell%dam(:,:)
      cell_shared%edam(:) = cell_shared%edam(:) + cell%edam(:)
!
!  Fix negative values of dens if applicable
!
      nneg = count( cell_shared%dens(:,:) < ZERO )

      if( nneg > 0 ) then
         print '(a)', 'mergeCellReg: Fixing ' // str( nneg ) // ' negative dens values'
         if( ndim == 1 ) then
            call fixNegDensCell1dReg( cell_shared%dens(:,:) )
!                ===================           io
         else if( ndim == 2 ) then
            dens2d(1:nelem,1:grid(1)%nx,1:grid(2)%nx) => cell_shared%dens
            call fixNegDensCell2dReg( dens2d(:,:,:) )
!                ===================       io
         end if
      end if
!$omp end master
!$omp barrier
!
!  Copy cell_shared to cell
!
      cell%dens(:,:) = cell_shared%dens(:,:)
      if( out%lhis )  cell%dam(:,:) = cell_shared%dam(:,:)
      cell%edam(:) = cell_shared%edam(:)
!
!  Calculate derived quantities
!
      do ic = 1, cell%nc
         cell%denstot(ic) = sum( cell%dens(:,ic) )
         if( cell%denstot(ic) == ZERO ) then
            cell%prob(:,ic) = ZERO
         else
            cell%prob(:,ic) = cell%dens(:,ic) / cell%denstot(ic)
         end if
      end do
!
!  Reinitialize cell_tmp
!
      cell_tmp%dens(:,:) = ZERO
      if( out%lhis )  cell_tmp%dam(:,:) = ZERO
      cell_tmp%edam(:) = ZERO

   end subroutine mergeCellReg
!
!=======================================================================
!
   subroutine relaxCellReg( init_ )
!             ============    i
!
!     Redistribute the cell densities so that stress relaxes.
!
!     This subroutine should not be called too often in order to avoid artificial diffusion,
!     but often enough to keep the surface sufficiently accurate.
!
!     OpenMP: This subroutine must be called from the master thread
!
!  Arguments:
!
      logical, optional :: init_                ! init_=t indicates call directly after init
                                                !    if unspecified, init_=f is assumed
!
!  Local Variables:
!
      logical :: init
      logical, allocatable :: lxtal2d(:,:)
      integer(INT4) :: ic, nx, nz, nc
      real(REAL8) :: fracmax, ed, denstot
      real(REAL8), allocatable :: ff(:), dpa(:)
      real(REAL8), allocatable :: dens2d(:,:,:), dens02d(:,:,:), ffgeom2d(:,:)
      real(REAL8), allocatable :: dam2d(:,:,:), dpa2d(:,:), pressure2d(:,:), stress2d(:,:,:,:)

      logical, save :: first = .true.           ! first call except init
      integer(INT4), save :: istep = 0          ! counter for relaxation (only used for 2D)
!
!  Start of Executable Code:
!  ------------------------
!
      if( present( init_ ) ) then
         init = init_
      else
         init = .false.
      end if

      if( first ) then
         first = .false.
         if( dam%relaxmodel == ENUM%relaxmodel%viscoelastic ) then
            call execute_command_line( "find . -maxdepth 1 -type d " //  &
                                       "-regextype posix-extended -regex './[0-9]+' " //  &
                                       "-exec rm -rf {} +" )
         end if
      end if
!
!  Initialize dpa
!
      call alloc( dpa, cell%nc, 'dpa', 'relaxCellReg' )
!          =====  io      i       i          i
      dpa(:) = ZERO

      if( .not.init ) then
!
!  Calculate dpa's
!
         if( any( dam%diffusivity(:) /= ZERO ) .or.  &
             dam%relaxmodel == ENUM%relaxmodel%viscoelastic ) then

            do ic = 1, cell%nc
               ed = sum( cell%ed(:,ic) * cell%prob(:,ic) )
               denstot = HALF * (sum( cell%dens0(:,ic) ) + cell%denstot(ic))
               if( ed == 0 .or. denstot == ZERO ) then   ! edam/=0, denstot=0 possible because
                  dpa(ic) = ZERO                         !    of stopped ions and recoils
                  !if( ed == ZERO .or. denstot == ZERO )  print *, 'ic,edam=', ic, cell%edam(ic)
               else
                  if( ed == ZERO )  print *, 'WARNING: ed=0'
                  if( denstot == ZERO )  print *, 'WARNING: denstot=0'
                  dpa(ic) = cell%edam(ic) / (2.5 * ed) / denstot
               end if
            end do
         end if
!
!  Increment istep and create step directory
!
         if( dam%relaxmodel == ENUM%relaxmodel%viscoelastic ) then
            istep = istep + 1
            call execute_command_line( 'mkdir ' // str( istep ) )
!                ====================             i
         end if
!
!  Diffuse implanted species
!
         if( any( dam%diffusivity(:) /= ZERO ) ) then
            if( any( dpa(:) > ZERO ) ) then

               if( ndim == 1 ) then
                  call diffuseCell1dReg( grid(1)%x(:), cell%dens(:,:), dpa(:) )
!                      ================       i             io           i
               else if( ndim == 2 ) then
#ifdef AMAT
                  call diffuseCell2dReg( grid(1)%x(:), grid(2)%x(:), cell%ffgeom(:),  &
!                      ================       i             i              i
                                         cell%dens(:,:), dpa(:), istep )
!                                              io          i       i
#else
                  error stop ' relaxCellReg: 2D diffusion requires separate license'
#endif
               else
                  error stop 'relaxCellReg: 3D not yet implemented'
               end if

               do ic = 1, cell%nc
                  cell%denstot(ic) = sum( cell%dens(:,ic) )
                  cell%prob(:,ic) = cell%dens(:,ic) / cell%denstot(ic)
               end do
            end if
         end if
!
!  Check fracupd
!
         call alloc( ff, cell%nc, 'ff', 'relaxCellReg' )
   !          =====  io      i     i           i
         do ic = 1, cell%nc
            ff(ic) = sum( cell%dens(:,ic) / elem(:)%dens )
         end do

         fracmax = maxval( abs( ff(:)/cell%ffgeom(:) - ONE ), mask=(cell%ffgeom(:)>ZERO) )

         if( .not.init .and. fracmax < geom%fracupd )  return

      end if
!
!  Call dimension-specific subroutines
!
      if( dam%relaxmodel == ENUM%relaxmodel%viscoelastic ) then

         if( ndim == 1 ) then
            call relaxCell1dReg( grid(1)%x,  &
!                ==============     io
                                 cell%dens, cell%dam, cell%dens0, cell%ffgeom, cell%lxtal,  &
!                                    io         io        io          io           io
                                 dpa, cell%pressure, cell%stress, istep )
!                                io         io           io         i
         else if( ndim == 2 ) then
            nx = grid(1)%nx
            nz = grid(2)%nx
            dens2d = reshape( cell%dens, [nelem,nx,nz] )
            if( out%lhis )  dam2d = reshape( cell%dam, [nelemIn,nx,nz] )
            dens02d = reshape( cell%dens0, [nelem,nx,nz] )
            ffgeom2d = reshape( cell%ffgeom, [nx,nz] )
            lxtal2d = reshape( cell%lxtal, [nx,nz] )
            dpa2d = reshape( dpa, [nx,nz] )
            pressure2d = reshape( cell%pressure, [nx,nz] )
            stress2d = reshape( cell%stress, [2,2,nx,nz] )

            call relaxCell2dReg( grid(1)%x, grid(2)%x,  &
!                ==============      i          i
                                 dens2d, dam2d, dens02d, ffgeom2d, lxtal2d, dpa2d,  &
!                                  io     io      io        io       io      io
                                 pressure2d, stress2d, istep )
!                                    io         io       i

            nx = ubound( grid(1)%x, 1 )
            nz = ubound( grid(2)%x, 1 )
            nc = nx * nz
            cell%dens = reshape( dens2d, [nelem,nc] )
            if( out%lhis )  cell%dam = reshape( dam2d, [nelemIn,nc] )
            cell%dens0 = reshape( dens02d, [nelem,nc] )
            cell%ffgeom = reshape( ffgeom2d, [nc] )
            cell%lxtal = reshape( lxtal2d, [nc] )
            cell%pressure = reshape( pressure2d, [nc] )
            cell%stress = reshape( stress2d, [2,2,nc] )
         else
            error stop 'relaxCellReg: 3D relaxation not yet implemented'
         end if

      else

         if( ndim /= 1 )  error stop 'relaxCellReg: Full relaxation only in 1D'

         call relaxCell1dRegFull( grid(1)%x, grid(1)%xmin1, grid(1)%xmax2,  &
!             ==================     io            i              i
                                  cell%dens, cell%dam, cell%ffgeom, cell%lxtal )
!                                    io         io         io           io

      end if
!
!  Reinitialize derived quantities
!
      call initDerivedCellReg( )
!          ==================

   end subroutine relaxCellReg
!
!=======================================================================
!
   subroutine writeCellReg( )
!             ============
!
!     Write the contents of the target cells to the CELL* file (same format as HIS file)
!
!  Arguments:
!
!     None
!
!  Local Variables:
!
      integer(INT4) :: nx, nz
      real(REAL8), allocatable :: dens(:,:,:), pressure(:,:), stress(:,:,:,:)
!
!  Start of Executable Code:
!  ------------------------
!
      if( ndim == 1 ) then
         call writeCell1dReg( grid(1)%x, cell%dens, cell%pressure, cell%stress )
!             ==============      i          i            i             i
      else if( ndim == 2 ) then
         nx = grid(1)%nx
         nz = grid(2)%nx
         dens = reshape( cell%dens, [nelem,nx,nz] )
         if( dam%relaxmodel == ENUM%relaxmodel%viscoelastic ) then
            pressure = reshape( cell%pressure, [nx,nz] )
            stress = reshape( cell%stress, [2,2,nx,nz] )
         end if
         call writeCell2dReg( grid(1)%x, grid(2)%x, dens, pressure, stress )
!             ==============      i          i       i       i        i
      else
         error stop 'writeCellReg: 3D not not yet implemented'
      end if

   end subroutine writeCellReg
!
!=======================================================================
!
   subroutine writeDamCellReg( )
!             ===============
!
!     Write the damage of the target cells to the HIS* file
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
         call writeDamCell1dReg( grid(1)%x, cell%dam )
!             =================      i         i
      else if( ndim == 2 ) then
!         call writeDamCell2dReg( grid(1)%x, grid(2)%x, cell%dam )
!!             =================      i          i         i
         error stop 'writeDamCellReg: 2D not not yet implemented'
      else
         error stop 'writeDamCellReg: 3D not not yet implemented'
      end if

   end subroutine writeDamCellReg
!
!=======================================================================
!
   subroutine initDerivedCellReg( )
!             ==================
!
!     Initialize derived quantities of the cell structure of the 1d target:
!     grid%(zmin*,zmax*), cell%(denstot, dens0,prob,vib,ed,esurf,edam)
!
!  Arguments:
!
!     None
!
!  Local Variables:
!
      integer(INT4) :: i, ic, ielem, jc
!
!  Start of Executable Code:
!  ------------------------
!
      do i = 1, ndim

         grid(i)%nx = ubound( grid(i)%x, 1 )
         grid(i)%wx = grid(i)%x(1) - grid(i)%x(0)
!
!  Smallest and largest position safely inside target
!
         if( grid(i)%x(0) == ZERO ) then
            grid(i)%xmin1 = -SMALL
            grid(i)%xmin2 = +SMALL
         else
            grid(i)%xmin1 = grid(i)%x(0) - EPSMA * abs(grid(i)%x(0))
            grid(i)%xmin2 = grid(i)%x(0) + EPSMA * abs(grid(i)%x(0))
         end if

         if( grid(i)%x(grid(i)%nx) == ZERO ) then
            grid(i)%xmax1 = -SMALL
            grid(i)%xmax2 = +SMALL
         else
            grid(i)%xmax1 = grid(i)%x(grid(i)%nx) - EPSMA*abs(grid(i)%x(grid(i)%nx))
            grid(i)%xmax2 = grid(i)%x(grid(i)%nx) + EPSMA*abs(grid(i)%x(grid(i)%nx))
         end if

      end do
!
!  Determine shape of (hypothetical) cell array
!
      cell%nx(:) = [ 1, 1, 1 ]
      cell%nx(1:ndim) = grid(1:ndim)%nx
      cell%size = product( grid(1:ndim)%wx )

      cell%nc = product( cell%nx(:) )
      if( cell%nc /= size( cell%dens, 2 ) )  &
            error stop 'initDerivedCellReg: inconsistent number of cells'
!
!  Derived quantities (note that array sizes might have changed)
!
      if( allocated( cell%denstot ) )  deallocate( cell%denstot )
      call alloc( cell%denstot, cell%nc, 'cell%denstot', 'initDerivedCell1dReg' )
!          =====        io         i           i                   i
      cell%denstot(:) = sum( cell%dens(:,:), dim=1 )

      if( allocated( cell%dens0 ) )  deallocate( cell%dens0 )
      call alloc( cell%dens0, [nelem, cell%nc], 'cell%dens0', 'initDerivedCell1dReg' )
!          =====        io           i                i                  i
      cell%dens0(:,:) = cell%dens(:,:)

      if( allocated( cell%prob ) )  deallocate( cell%prob )
      call alloc( cell%prob, [nelem, cell%nc], 'cell%prob', 'initDerivedCell1dReg' )
!          =====      io           i                i                 i
      do ic = 1, cell%nc
         if( cell%denstot(ic) == ZERO ) then
            cell%prob(:,ic) = ZERO
         else
            cell%prob(:,ic) = cell%dens(:,ic) / cell%denstot(ic)
         end if
      end do

      if( allocated( cell%vib ) )  deallocate( cell%vib )
      call alloc( cell%vib, [nelem, cell%nc], 'cell%vib', 'initDerivedCell1dReg' )
!          =====     io           i               i                 i
      do ic = 1, cell%nc
         cell%vib(:,ic) = elem(:)%vib
      end do

      if( allocated( cell%ed ) )  deallocate( cell%ed )
      call alloc( cell%ed, [nelem, cell%nc], 'cell%ed', 'initDerivedCell1dReg' )
!          =====    io           i               i                i
      if( allocated( cell%esurf ) )  deallocate( cell%esurf )
      call alloc( cell%esurf, [nelem, cell%nc], 'cell%esurf', 'initDerivedCell1dReg' )
!          =====      io            i             i                  i

      do ic = 1, cell%nc
         if( cell%denstot(ic) /= ZERO ) then
            jc = ic
         else
            if( ndim == 1 ) then
               call getNearestCell1dReg( ic, cell%denstot(:), jc )
!                   ===================  i          i         o
            else if( ndim == 2 ) then
               call getNearestCell2dReg( ic, cell%nx(1), cell%nx(2), grid(1)%wx, grid(2)%wx,  &
!                   ===================  i        i          i              i           i
                                         cell%denstot(:), jc )
!                                                i        o
            else
               error stop 'initDerivedCellReg: 3D not yet implemented'
            end if
         end if

         do ielem = 1, nelem
            cell%ed(ielem,ic) = sum( cell%prob(:,jc) * dam%ed(:,ielem) )
         end do

         do ielem = 1, nelem
            cell%esurf(ielem,ic) = sum( cell%prob(:,jc) * dam%esurf(:,ielem) )
         end do
      end do

      if( allocated( cell%edam ) )  deallocate( cell%edam )
      call alloc( cell%edam, cell%nc, 'cell%edam', 'initDerivedCell1dReg' )
!          =====     io          i           i                i
      cell%edam(:) = ZERO
!
!  (Re-)Initialize cell_tmp
!
      if( nionhis > 0 ) then

         if( allocated( cell_tmp%dens ) )  deallocate( cell_tmp%dens )
         call alloc( cell_tmp%dens, [nelem, cell%nc], 'cell_tmp%dens', 'initDerivedCell1dReg' )
!             =====       io              i                  i                    i
         cell_tmp%dens(:,:) = ZERO

         if( out%lhis ) then
            if( allocated( cell_tmp%dam ) )  deallocate( cell_tmp%dam )
            call alloc( cell_tmp%dam, shape( cell%dam ),  &
!                =====       io               i
                        'cell_tmp%dam', 'initDerivedCell1dReg' )
!                              i                     i
            cell_tmp%dam(:,:) = ZERO
         end if

         if( allocated( cell_tmp%edam ) )  deallocate( cell_tmp%edam )
         call alloc( cell_tmp%edam, cell%nc, 'cell_tmp%edam', 'initDerivedCell1dReg' )
!             =====       io           i            i                    i
         cell_tmp%edam(:) = ZERO

      end if

   end subroutine initDerivedCellReg
!
!=======================================================================
!
end module Target_Cell_Reg_mod
