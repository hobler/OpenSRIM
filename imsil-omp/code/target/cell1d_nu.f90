module Target_Cell1d_NU_mod
!
!  Holds 1d cell structure on a non-uniform grid and provides subroutines operating on it
!
!  Note: avoid ldyn=t with lspher=t
!
!  Used Modules:
!  ============
!
!$ use omp_lib,               only: omp_get_thread_num
   use Const_mod,             only: INT4, REAL8,  &
                                    ZERO, ONE, TWO, HALF, BIG, EPSMA, SMALL
   use Array_mod,             only: alloc, check_allocation
   use Open_mod,              only: openFile
   use String_mod,            only: str ! debug
   use Types_mod,             only: Proj_type
   use Undefine_mod,          only: defined

   use Geom_BndCond_mod,      only: applyBndCond1d
   use Geom_CurviCoord_mod,   only: toCurviCoord1d
   use Math_BinSearch_mod,    only: binSearch

   use Param_Setup_mod,       only: nr, nelem, nelemIon, nelemTarg, nionhis,  &
                                    ldyn, usecell, cellfile, cellrec, usehis, hisfile,  &
                                    reflect, shiftto
   use Param_Geom_mod,        only: geom
   use Param_Elem_mod,        only: elem
   use Param_Ion_mod,         only: ion
   use Param_Mat_mod,         only: mat, irxtal
   use Param_Dam_mod,         only: dam
   use Param_Out_mod,         only: out

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
   public :: initCell1dNU                 ! initialize the target
!            ============
   public :: scoreCell1dNU                ! score a count in the target cell
!            =============
   public :: getCellIndex1dNU             ! determines the cell a given point is located in
!            ================
   public :: getOtherNearCellIndices1dNU  ! determine cells other than nearest cell that
!            ===========================  !    intersect with the free flight path cylinder
   public :: getProbCell1dNU              ! get the probabilities to encounter each of the atoms
!            ===============              !    in a given cell
   public :: getDensCell1dNU              ! get the total atomic density in a given cell
!            ===============
   public :: getVibCell1dNU               ! get the array of thermal vibration amplitudes
!            ==============               !    in a given cell
   public :: getXtalCell1dNU              ! get the crystallinity flag of the cell
!            ===============
   public :: getEdCell1dNU                ! get the displacement energy after most recent
!            =============                !    target cell update
   public :: getEsurfCell1dNU             ! get the surface binding energy after most recent
!            ================             !    target cell update
   public :: broadcastCell1dNU            ! broadcast cells from master thread to other threads
!            =================
   public :: mergeCell1dNU                ! merge the contents of cell_tmp into cell
!            =============
   public :: relaxCell1dNU                ! relax the cell densities (target cell update)
!            =============
   public :: writeCell1dNU                ! write target cell contents to output file
!            =============
   public :: writeDamCell1dNU             ! write damage to output file
!            ================
   public :: readAtomsFromCell1dNU        ! read atomic numbers from cell file
!            =====================
   public :: readBndFromCell1dNU          ! obtain the geometry boundaries (surface positions)
!            ===================          !    from cell file
   public :: readCell1dNUraw              ! read target cell contents from file returning arrays
!            ===============
!
!  Private Procedures:
!  ==================
!
   private :: readCell1dNU                ! read target cell contents from file
!             ============
   private :: readHis1dNU                 ! read damage from file
!             ===========
   private :: initCell1dNUFromGeom        ! initialize grid and cell from region geometry
!             ====================
   private :: initDerivedCell1dNU         ! initialize derived quantities
!             ===================
   private :: getCellIndex1dNU_core       ! core part of getCellIndex1dNU working with z only
!             =====================
   private :: diffuseCell1dNU             ! diffuse the implanted species
!             ===============
!
!  Private Data:
!  ============
!
   type :: Grid_type
!          =========
      integer(INT4) :: ncell                 ! number of cells
      real(REAL8), allocatable :: wcell(:)   ! (ncell) cell widths; in case of surface
!                                            !   roughness, after relaxation considering
!                                            !   density reduction due to roughness
      real(REAL8), allocatable :: posif(:)   ! (0:ncell) positions of interfaces between
                                             !   cells
      real(REAL8) :: posmin1                 ! position just left of posif(0)
      real(REAL8) :: posmin2                 ! position just right of posif(0)
      real(REAL8) :: posmax1                 ! position just left of posif(ncell)
      real(REAL8) :: posmax2                 ! position just right of posif(ncell)
   end type Grid_type

   type(Grid_type) :: grid                   ! grid
   
   type :: Cell_type
!          =========
      real(REAL8), allocatable :: dens(:)    ! (nelem) density of atom ielem in cell
      real(REAL8) :: dens_tot                ! total atom density in cell
      real(REAL8), allocatable :: prob(:)    ! (nelem) probability to encounter atom
      real(REAL8), allocatable :: dam(:)     ! (nelemIn) damage densities
      real(REAL8), allocatable :: vib(:)     ! (nelem) thermal vibration amplitude
                                             !    ielem in cell
      logical :: lxtal                       ! flag indicating crystallinity
      real(REAL8), allocatable :: ed(:)      ! (nelem) displacement energy after
                                             !    most recent target cell update
      real(REAL8), allocatable :: esurf(:)   ! (nelem) surface binding energy after
                                             !    most recent target cell update
      real(REAL8) :: edam                    ! damage energy per cm3 since last relax
   end type Cell_type

   type(Cell_type), allocatable :: cell(:)      ! (grid%ncell) cells
   type(Cell_type), allocatable :: cell_tmp(:)  ! same as cell but containing density and
                                                !    damage changes since last update
!$omp threadprivate( cell, cell_tmp )

   integer(INT4) :: nelemIn
   real(REAL8), parameter :: HISMIN = 0.0001
!
!=======================================================================
!
contains
!
!=======================================================================
!
   subroutine initCell1dNU( )
!             ============
!
!     Initialize the cell structure of the 1d target
!
!     OpenMP: This subroutine should be called from the master thread. On exit, cell
!     and cell_tmp will be initialized on the master thread only
!
!  Arguments:
!
!     None
!
!  Local Variables:
!
      character(80) :: text
      character(256) :: msg
      integer(INT4) :: ic, istat
!
!  Start of Executable Code:
!  ------------------------
!
      if( usecell ) then

         call readCell1dNU( )
!             ============
      else

         call initCell1dNUFromGeom( )
!             ====================
      end if
!
!  Cell damage
!
      nelemIn = 2 * nelem        ! maximum atom index (note: ielem>nelem ... vacancies)

      if( out%lhis ) then
         do ic = 1, grid%ncell
            allocate( cell(ic)%dam( nelemIn ), stat=istat, errmsg=msg )
            write( text, '(a,i0,a)' )  'cell(', ic, ')%dam'
            call check_allocation( 'initCell1dNU', text, istat, msg )
            cell(ic)%dam(:) = ZERO
         end do
      end if

      if( usehis )  call readHis1dNU( hisfile )
!                        ===========     i
!
!  Cell damage energy since last relax
!
      cell(:)%edam = ZERO                ! needed in relaxCell1dNU
!
!  Initialize cell_tmp
!
      if( nionhis > 0 ) then

         allocate( cell_tmp( grid%ncell ), stat=istat, errmsg=msg )
         call check_allocation( 'initCell1dNU', 'cell_tmp', istat, msg )

         do ic = 1, grid%ncell
            allocate( cell_tmp(ic)%dens( nelem ), stat=istat, errmsg=msg )
            write( text, '(a,i0,a)' )  'cell_tmp(', ic, ')%dens'
            call check_allocation( 'initCell1dNU', text, istat, msg )
            cell_tmp(ic)%dens(:) = ZERO
         end do

         if( out%lhis ) then
            do ic = 1, grid%ncell
               allocate( cell_tmp(ic)%dam( nelemIn ), stat=istat, errmsg=msg )
               write( text, '(a,i0,a)' )  'cell_tmp(', ic, ')%dam'
               call check_allocation( 'initCell1dNU', text, istat, msg )
               cell_tmp(ic)%dam(:) = ZERO
            end do
         end if

         cell_tmp(:)%edam = ZERO

      end if
!
!  Initialize the surface
!
      call updateSurface1d( grid%posif(0), grid%posif(grid%ncell) )
!          ===============        i                   i
!
!  Calculate derived quantities
!
      call initDerivedCell1dNU( )
!          ===================
!
!  Relax the cells
!
      call relaxCell1dNU( first_=.true. )
!          =============       i

   end subroutine initCell1dNU
!
!=======================================================================
!
   subroutine scoreCell1dNU( ielem, ic, w, ddami, ddamv, edam )
!             =============    i    i   i    i      i     i
!
!     Score a count in the 1d cell histogram
!
!     The cell contents may go below zero.
!
!  Arguments:
!
      integer(INT4), intent(in) :: ielem        ! index of atom species
      integer(INT4), intent(in) :: ic           ! cell index
      real(REAL8), intent(in) :: w              ! atom increment (can also be negative)
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
!      print *, 'ic,ielem=', ic, ielem
!      print *, 'w, ddami, ddamv=', w, ddami, ddamv
!      print *, 'cell(ic)%dens(ielem)=', cell(ic)%dens
!      if( ielem==1 ) print *, 'ic, w, ddami, ddamv=', ic, w, ddami, ddamv
!
!  Cell density
!
      if( w /= ZERO ) then
         ddens = w * ion%ddose / grid%wcell(ic)
         if( nionhis > 0 ) then
            cell_tmp(ic)%dens(ielem) = cell_tmp(ic)%dens(ielem) + ddens
         else
            cell(ic)%dens(ielem) = cell(ic)%dens(ielem) + ddens
            cell(ic)%dens(ielem) = cell(ic)%dens(ielem)
            cell(ic)%dens_tot = sum( cell(ic)%dens(:) )
            cell(ic)%prob(:) = cell(ic)%dens(:) / cell(ic)%dens_tot
!
!         print *, 'cell(ic)%dens(ielem)=', cell(ic)%dens
!
         end if
      end if
!
!  Damage densities
!
      if( out%lhis ) then
         if( ddami /= ZERO ) then
            ddam = ddami * ion%ddose / grid%wcell(ic)
            if( nionhis > 0 ) then
               cell_tmp(ic)%dam(ielem) = cell_tmp(ic)%dam(ielem) + ddam
            else
               cell(ic)%dam(ielem) = cell(ic)%dam(ielem) + ddam
               cell(ic)%dam(ielem) = max( SMALL, cell(ic)%dam(ielem) )
            end if
         end if
         if( ddamv /= ZERO ) then
            ddam = ddamv * ion%ddose / grid%wcell(ic)
            if( nionhis > 0 ) then
               cell_tmp(ic)%dam(ielem+nelem) = cell_tmp(ic)%dam(ielem+nelem) + ddam
            else
               cell(ic)%dam(ielem+nelem) = cell(ic)%dam(ielem+nelem) + ddam
               cell(ic)%dam(ielem+nelem) = max( SMALL, cell(ic)%dam(ielem+nelem) )
            end if
         end if
      end if
!
!  Damage energy per cell volume
!
      if( edam /= ZERO ) then
         dedam = edam * ion%ddose / grid%wcell(ic)
         if( nionhis > 0 ) then
            cell_tmp(ic)%edam = cell_tmp(ic)%edam + dedam
         else
            cell(ic)%edam = cell(ic)%edam + dedam
         end if
      end if
!
!  Do not update ed and esurf to avoid unreasonable values when cell gets empty!
!
   end subroutine scoreCell1dNU
!
!=======================================================================
!
   subroutine getCellIndex1dNU( pos, ic, outside_ )
!             ================   i   o      o
!  
!     Determine the cell the point pos is located in or closest to.
!
!  Arguments:
!
      real(REAL8), intent(in) :: pos(3)         ! position to be checked
      integer(INT4), intent(out) :: ic          ! index of cell (region)
      logical, optional, intent(out) :: outside_! flag indicationg whether pos is outside any
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

      call getCellIndex1dNU_core( z, ic, outside_ )
!          =====================  i   o      o

   end subroutine getCellIndex1dNU
!
!=======================================================================
!
   subroutine getOtherNearCellIndices1dNU( proj, posc, ffp, pmax, ncn, icnear )
!             ===========================   i     i     i    i     o     o
!
!     Determine cells other than proj%id that are within the free flight path cylinder
!     (= cylinder with height ffp and radius pmax)
!
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
                                                ! index of regions found except proj%id
!
!  Local Variables:
!
      integer(INT4) :: i, ic, icmin, icmax
      real(REAL8) :: z, dz, z1, z2, zmin, zmax, distother, zcen, ffpcen(3)
!
!  Start of Executable Code:
!  ------------------------
!
!  If there is only one cell there is no "other" cell
!
      if( grid%ncell == 1 ) then
         ncn = 0
         return
      end if

      ic = proj%id                              ! cell to be excluded from the search

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
      if( z < grid%posif(1) ) then           
         distother = grid%posif(1) - z
      else if( z > grid%posif(grid%ncell-1) ) then
         distother = z - grid%posif(grid%ncell-1)
      else
         distother = min( grid%posif(ic) - z, z - grid%posif(ic-1) )
      end if

      if( distother**2 > ffp**2 + pmax**2 ) then
         ncn = 0
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
!             ==============      i       o
         dz = sqrt( (HALF*ffp)**2 + pmax**2 )
         zmin = zcen - dz
         zmax = zcen + dz
      end if

      zmin = max( grid%posmin2, min( grid%posmax1, zmin ) )      ! don't apply bounds to 
      zmax = max( grid%posmin2, min( grid%posmax1, zmax ) )      !    icmin,icmax because 
                                                                 !    of integer overflow
      call binSearch( grid%ncell+1, grid%posif, zmin, icmin )
!          =========       i             i       i      o        ! icmin,icmax are in
      call binSearch( grid%ncell+1, grid%posif, zmax, icmax )    !    [1,ncell]
!          =========       i             i       i      o
!
!  Make sure icnear is big enough
!
      if( size( icnear ) < grid%ncell-1 ) then
         deallocate( icnear )
         allocate( icnear(grid%ncell-1) )
      end if
!
!  Set icnear array
!
      if( ic < icmin .or. ic > icmax ) then
         print *, 'icmin, icmax, ic = ', icmin, icmax, ic
         error stop 'getOtherNearCellInidices1d: cell index error'
         ncn = icmax - icmin + 1
         icnear(1:ncn) = [ ( i, i=icmin,icmax ) ]
      else
         ncn = icmax - icmin
         icnear(1:ncn) = [ ( i, i=icmin,ic-1 ), ( i, i=ic+1,icmax ) ]
      end if

   end subroutine getOtherNearCellIndices1dNU
!
!=======================================================================
!
   subroutine getProbCell1dNU( ic, prob )
!             ===============  i    o
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
      prob(:) = cell(ic)%prob(:)

   end subroutine getProbCell1dNU
!
!=======================================================================
!
   subroutine getDensCell1dNU( ic, dens )
!             ===============  i    o
!
!     Get the total atomic density in cell ic
!
!     Zero is returned instead of a negative value
!
!     NOTE: In case of surface roughness this is the density before reduction due to roughness
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
      dens = max( cell(ic)%dens_tot, ZERO )

   end subroutine getDensCell1dNU
!
!=======================================================================
!
   subroutine getVibCell1dNU( ic, vib )
!             ==============  i    o
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
      vib(:) = cell(ic)%vib(:)

   end subroutine getVibCell1dNU
!
!=======================================================================
!
   subroutine getXtalCell1dNU( ic, lxtal )
!             =============  i     o
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
      lxtal = cell(ic)%lxtal

   end subroutine getXtalCell1dNU
!
!=======================================================================
!
   subroutine getEdCell1dNU( ic, ed )
!             =============  i   o
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
      ed(:) = cell(ic)%ed(:)

   end subroutine getEdCell1dNU

!=======================================================================
!
   subroutine getEsurfCell1dNU( ic, esurf )
!             ================  i     o
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
      esurf(:) = cell(ic)%esurf(:)

   end subroutine getEsurfCell1dNU
!
!=======================================================================
!
   subroutine broadcastCell1dNU( )
!             ===============
!
!     Copy cells from master thread to other threads
!
!  Arguments:
!
!     None
!
!  Local Variables:
!
      type(Cell_type), allocatable, save :: cell_shared(:)   ! shared temporary cell
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

   end subroutine broadcastCell1dNU
!
!=======================================================================
!
   subroutine mergeCell1dNU( )
!             =============
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
      integer(INT4) :: ithread, ic
      type(Cell_type), allocatable, save :: cell_tmp_shared(:)   ! shared temporary cell
!
!  Start of Executable Code:
!  ------------------------
!
!  Sum up all cell_tmp into cell_tmp_shared
!
!$omp master
      cell_tmp_shared = cell_tmp
!$omp end master

!$omp barrier
!$omp critical( critical_cell1d_mergeCell1d )
      ithread = 1
!$    ithread = 1 + omp_get_thread_num( )
      if( ithread > 1 ) then
         do ic = 1, grid%ncell
            cell_tmp_shared(ic)%dens(:) = cell_tmp_shared(ic)%dens(:) + cell_tmp(ic)%dens(:)
            if( out%lhis )  &
                  cell_tmp_shared(ic)%dam(:) = cell_tmp_shared(ic)%dam(:) + cell_tmp(ic)%dam(:)
            cell_tmp_shared(ic)%edam = cell_tmp_shared(ic)%edam + cell_tmp(ic)%edam
         end do
      end if
!$omp end critical( critical_cell1d_mergeCell1d )
!$omp barrier
!
!  Add cell_tmp_shared to cell
!
      do ic = 1, grid%ncell
         cell(ic)%dens(:) = cell(ic)%dens(:) + cell_tmp_shared(ic)%dens(:)
         if( out%lhis )  cell(ic)%dam(:) = cell(ic)%dam(:) + cell_tmp_shared(ic)%dam(:)
         cell(ic)%edam = cell(ic)%edam + cell_tmp_shared(ic)%edam
      end do
!
!  Calculate derived quantities
!
      do ic = 1, grid%ncell
         cell(ic)%dens_tot = sum( cell(ic)%dens(:) )
         cell(ic)%prob(:) = cell(ic)%dens(:) / cell(ic)%dens_tot
      end do
!
!  Reinitialize cell_tmp
!
      do ic = 1, grid%ncell
         cell_tmp(ic)%dens(:) = ZERO
         if( out%lhis )  cell_tmp(ic)%dam(:) = ZERO
         cell_tmp(ic)%edam = ZERO
      end do

   end subroutine mergeCell1dNU
!
!=======================================================================
!
   subroutine relaxCell1dNU( first_ )
!             =============    i
!
!     Adjust the cells so that stress relaxes. This is done in two steps:
!     First the cell sizes are adjusted so that the equilibrium densities are assumed. Then
!     the cells are merged or split to maintain an approximately constant grid spacing.
!
!     For nionhis>0, this subroutine must be called after merging the cell_tmp into cell.
!
!     This subroutine should not be called too often in order to avoid unnecessary computations,
!     but often enough to keep the cell densities sufficiently accurate.
!
!     OpenMP: This subroutine must be called from within a parallel region. It is sufficient
!     that cell and cell_tmp are up to date on the master thread. On exit, cell and cell_tmp
!     are updated on all threads.
!
!  Arguments:
!
      logical, optional :: first_               ! flag indicating first call
!
!  Local Variables:
!
      logical :: first, lmerge(2:grid%ncell), lsplit(grid%ncell)
      logical, allocatable :: lxtalnew(:)
      character(80) :: text
      character(256) :: msg
      integer(INT4) :: ielem, ic, icfix, icnew, icnext, jc, istat
      integer(INT4) :: nmerge, nsplit
      integer(INT4) :: ncellold, ncellnew
      integer(INT4) :: icfixnew
      real(REAL8) :: ed, ddz, wsum, posfix
      real(REAL8) :: fracmax
      real(REAL8) :: posfixrlx0
      real(REAL8) :: posfixrlx1
      real(REAL8), allocatable :: ff(:), wrlx(:), densrlx(:,:)
      real(REAL8), allocatable :: wnew(:)
      real(REAL8), allocatable :: densnew(:,:)
      real(REAL8), allocatable :: damnew(:,:)
      real(REAL8), allocatable :: edamnew(:)
      real(REAL8), allocatable :: vibnew(:,:)
      real(REAL8), allocatable :: dens(:), dpa(:)
!
!  Start of Executable Code:
!  ------------------------
!
      if( present( first_ ) ) then
         first = first_
      else
         first = .false.
      end if

      ncellold = grid%ncell                     ! remember old number of cells
!
!  Diffuse the implanted species
!
      if( any( dam%diffusivity(:) /= ZERO ) ) then
         allocate( dpa( grid%ncell ), stat=istat, errmsg=msg )
         call check_allocation( 'relaxCell1dNU', 'dpa', istat, msg )
         allocate( dens( grid%ncell ), stat=istat, errmsg=msg )
         call check_allocation( 'relaxCell1dNU', 'dens', istat, msg )

         do ic = 1, grid%ncell
            ed = sum( cell(ic)%ed(:) * cell(ic)%prob(:) )
            dpa(ic) = cell(ic)%edam / (2.5 * ed) / sum( cell(ic)%dens(:) )
         end do

         do ielem = 1, nelemIon
            if( dam%diffusivity(ielem) == ZERO )  cycle
            do ic = 1, grid%ncell
               dens(ic) = cell(ic)%dens(ielem)
            end do
            call diffuseCell1dNU( grid%posif(:), dens(:), dpa(:), ielem )
!                ===============        i          io       i       i
            do ic = 1, grid%ncell
               cell(ic)%dens(ielem) = dens(ic)
            end do
         end do

         deallocate( dpa, dens )
      end if
!
!  Calculate relaxed cell sizes and densities
!
      allocate( ff( grid%ncell ), stat=istat, errmsg=msg )
      call check_allocation( 'relaxCell1dNU', 'ff', istat, msg )
      allocate( wrlx( grid%ncell ), stat=istat, errmsg=msg )
      call check_allocation( 'relaxCell1dNU', 'wrlx', istat, msg )
      allocate( densrlx( nelem, grid%ncell ), stat=istat, errmsg=msg )
      call check_allocation( 'relaxCell1dNU', 'densrlx', istat, msg )

      do ic = 1, grid%ncell
         ff(ic) = sum( cell(ic)%dens(:) / elem(:)%dens ) ! fill fraction = to what fraction the
      end do                                             !   volume is filled with atoms

      fracmax = ZERO
      do ic = 1, grid%ncell
         wrlx(ic) = grid%wcell(ic) * ff(ic)
         fracmax = max( fracmax, abs( wrlx(ic) - grid%wcell(ic) ) / out%wbox )
      end do
      if( fracmax < geom%fracupd )  return

      do ic = 1, grid%ncell
         densrlx(:,ic) = cell(ic)%dens(:) / ff(ic)
      end do
!
!  Determine posfix and the box corresponding to it. We will update the grid starting at this
!  box
!
      if( geom%posfix > grid%posif(0) .and. geom%posfix < grid%posif(grid%ncell) ) then
         posfix = geom%posfix
      else
         posfix = HALF * (grid%posif(0) + grid%posif(grid%ncell)) ! Is this necessary?
      end if

      call getCellIndex1dNU_core( posfix, icfix )
!          =====================    i       o
!
!  Label cells to be split
!
      lsplit(:) = (wrlx(:) > 1.5 * out%wbox)
!
!  Label cells to be merged with their left-side neighbors
!
      wsum = wrlx(icfix)                     ! cell width of icfix

      do ic = icfix+1, grid%ncell            ! cells right of icfix
         if( lsplit(ic-1) ) then
            lmerge(ic) = .false.
            wsum = wrlx(ic)
         else
            wsum = wsum + wrlx(ic)
            lmerge(ic) = (wsum < out%wbox .or.  &
                          wrlx(ic) < HISMIN*ion%ddose/sum(densrlx(:,ic)))
            if( .not.lmerge(ic) )  wsum = wrlx(ic)
         end if
      end do

      wsum = wrlx(icfix)                     ! (already merged) cell width of icfix
      do ic = icfix+1, grid%ncell
         if( .not.lmerge(ic) )  exit
         wsum = wsum + wrlx(ic)
      end do

      do ic = icfix-1, 1, -1                 ! cells left of icfix
         if( lsplit(ic+1) ) then
            lmerge(ic+1) = .false.
            wsum = wrlx(ic)
         else
            wsum = wsum + wrlx(ic)
            lmerge(ic+1) = (wsum < out%wbox .or.  &
                            wrlx(ic) < HISMIN*ion%ddose/sum(densrlx(:,ic)))
            if( .not.lmerge(ic+1) )  wsum = wrlx(ic)
         end if
      end do

      nmerge = count( lmerge(:) )
      nsplit = count( lsplit(:) )
      ncellnew = grid%ncell + nsplit - nmerge
!
!  Merge and split the labeled cells; store results in local arrays
!
      allocate( wnew( ncellnew ) )
      allocate( densnew( nelem, ncellnew ) )
      if( out%lhis )  allocate( damnew( nelemIn, ncellnew ) )
      allocate( edamnew( ncellnew ) )
      allocate( vibnew( nelem, ncellnew ) )
      allocate( lxtalnew( ncellnew ) )

      icnew = 1
      ic = 1
      icfixnew = ncellnew

      do
         if( lsplit(ic) ) then               ! split cell
            if( ic == icfix ) then
               if( posfix < HALF*(grid%posif(ic-1)+grid%posif(ic)) ) then
                  icfixnew = icnew
               else
                  icfixnew = icnew + 1
               end if
            end if

            wnew(icnew:icnew+1) = HALF*wrlx(ic)
            densnew(:,icnew) = cell(ic)%dens(:) * grid%wcell(ic) / wrlx(ic)
            densnew(:,icnew+1) = densnew(:,icnew)
            if( out%lhis ) then
               damnew(:,icnew) = cell(ic)%dam(:) * grid%wcell(ic) / wrlx(ic)
               damnew(:,icnew+1) = damnew(:,icnew)
            end if
            edamnew(icnew) = cell(ic)%edam * grid%wcell(ic) / wrlx(ic)
            edamnew(icnew+1) = edamnew(icnew)
            vibnew(:,icnew) = cell(ic)%vib(:)
            vibnew(:,icnew+1) = cell(ic)%vib(:)
            lxtalnew(icnew:icnew+1) = cell(ic)%lxtal

            icnew = icnew + 2
            ic = ic + 1
         else                                ! copy cell or merge cells
            do icnext = ic+1, grid%ncell
               if( .not.lmerge(icnext) )  exit
            end do

            if( icnext > icfix .and. icfixnew == ncellnew)  icfixnew = icnew

            wnew(icnew) = sum( wrlx(ic:icnext-1) )
            densnew(:,icnew) = ZERO
            do jc = ic, icnext-1
               densnew(:,icnew) = densnew(:,icnew) + cell(jc)%dens(:) * grid%wcell(jc)
            end do
            densnew(:,icnew) = densnew(:,icnew) / wnew(icnew)
            if( out%lhis ) then
               damnew(:,icnew) = ZERO
               do jc = ic, icnext-1
                  damnew(:,icnew) = damnew(:,icnew) + cell(jc)%dam(:) * grid%wcell(jc)
               end do
               damnew(:,icnew) = damnew(:,icnew) / wnew(icnew)
            end if
            edamnew(icnew) = ZERO
            do jc = ic, icnext-1
               edamnew(icnew) = edamnew(icnew) + cell(jc)%edam * grid%wcell(jc)
            end do
            edamnew(icnew) = edamnew(icnew) / wnew(icnew)

                                             ! for vib and lxtal take the values from the
            wsum = ZERO                      !   old cell containing the center of the new
            do jc = ic, icnext-1             !   cell
               wsum = wsum + wrlx(jc)
               if( wsum > HALF*wnew(icnew) ) exit
            enddo
            vibnew(:,icnew) = cell(jc)%vib(:)
            lxtalnew(icnew) = cell(jc)%lxtal

            icnew = icnew + 1
            ic = icnext
         end if

!         if( ic == grid%ncell-1 ) then       ! debug
!            continue
!         end if
         if( ic > grid%ncell )  exit
      end do

!      print *, 'densnew(:,ncellnew) = ', densnew(:,ncellnew)
      if( icnew /= ncellnew + 1 )  error stop 'relaxCell1dNU: Internal error with icnew.'
!
!  Determine grid positions of cell containing posfix
!  (TODO: check if this code breaks if icfix is within surface roughness layer)
!
      ddz = wrlx(icfix) - grid%wcell(icfix)     ! without consideration of merge and split...
      posfixrlx0 = grid%posif(icfix-1) -  &
                   ddz * ( posfix - grid%posif(icfix-1) ) / grid%wcell(icfix)
      posfixrlx1 = grid%posif(icfix) +  &
                   ddz * ( grid%posif(icfix) - posfix ) / grid%wcell(icfix)

      if( lsplit(icfix) ) then                  ! correction if cell is split
         if( posfix < HALF*(grid%posif(icfix-1)+grid%posif(icfix)) ) then
            posfixrlx1 = HALF*(posfixrlx0+posfixrlx1)
         else
            posfixrlx0 = HALF*(posfixrlx0+posfixrlx1)
         end if
      end if

      do ic = icfix+1, grid%ncell               ! correction if cell is merged with neighbors
         if( .not.lmerge(ic) )  exit            !   on the right
         posfixrlx1 = posfixrlx1 + wrlx(ic)
      end do
      do ic = icfix-1, 1, -1                    ! correction if cell is merged with neighbors
         if( .not.lmerge(ic+1) )  exit          !   on the left
         posfixrlx0 = posfixrlx0 - wrlx(ic)
      end do
!
!  Reallocate and set cell widths
!
      grid%ncell = ncellnew
      deallocate( grid%wcell, stat=istat )
      allocate( grid%wcell( grid%ncell ), stat=istat, errmsg=msg )
      call check_allocation( 'relaxCell1dNU', 'grid%wcell', istat, msg )

      grid%wcell(:) = wnew(:)
!
!  Reallocate and set grid positions
!
      deallocate( grid%posif )
      allocate( grid%posif( 0:grid%ncell ), stat=istat, errmsg=msg )
      call check_allocation( 'relaxCell1dNU', 'grid%posif', istat, msg )

      icfix = icfixnew
      grid%posif(icfix) = posfixrlx1
      do ic = icfix+1, grid%ncell
         grid%posif(ic) = grid%posif(ic-1) + grid%wcell(ic)
      end do

      grid%posif(icfix-1) = posfixrlx0
      do ic = icfix-2, 0, -1
         grid%posif(ic) = grid%posif(ic+1) - grid%wcell(ic+1)
      end do

!
!  Deallocate cell
!
      do ic = 1, ncellold                       ! to avoid a memory leak! (?)
         deallocate( cell(ic)%dens )
         if( out%lhis )  deallocate( cell(ic)%dam )
         deallocate( cell(ic)%vib )
      end do
      if( .not.first ) then
         do ic = 1, ncellold                    ! to avoid a memory leak! (?)
            deallocate( cell(ic)%prob )
            deallocate( cell(ic)%ed )
            deallocate( cell(ic)%esurf )
         end do
      end if
      deallocate( cell )

      if( nionhis > 0 ) then
         do ic = 1, ncellold                    ! to avoid a memory leak! (?)
            deallocate( cell_tmp(ic)%dens )
            if( out%lhis )  deallocate( cell_tmp(ic)%dam )
         end do
         deallocate( cell_tmp )
      end if
!
!  Allocate and set cell contents
!
      allocate( cell( grid%ncell ), stat=istat, errmsg=msg )
      call check_allocation( 'relaxCell1dNU', 'cell', istat, msg )
      if( nionhis > 0 ) then
         allocate( cell_tmp( grid%ncell ), stat=istat, errmsg=msg )
         call check_allocation( 'relaxCell1dNU', 'cell_tmp', istat, msg )
      end if

      do ic = 1, grid%ncell
         allocate( cell(ic)%dens( nelem ), stat=istat, errmsg=msg )
         write( text, '(a,i0,a)' )  'cell(', ic, ')%dens'
         call check_allocation( 'relaxCell1dNU', text, istat, msg )
         cell(ic)%dens(:) = densnew(:,ic)
      end do
      if( nionhis > 0 ) then
         do ic = 1, grid%ncell
            allocate( cell_tmp(ic)%dens( nelem ), stat=istat, errmsg=msg )
            write( text, '(a,i0,a)' )  'cell_tmp(', ic, ')%dens'
            call check_allocation( 'relaxCell1dNU', text, istat, msg )
            cell_tmp(ic)%dens(:) = ZERO
         end do
      end if

      if( out%lhis ) then
         do ic = 1, grid%ncell
            allocate( cell(ic)%dam( nelemIn ), stat=istat, errmsg=msg )
            write( text, '(a,i0,a)' )  'cell(', ic, ')%dam'
            call check_allocation( 'relaxCell1dNU', text, istat, msg )
            cell(ic)%dam(:) = damnew(:,ic)
         end do
         if( nionhis > 0 ) then
            do ic = 1, grid%ncell
               allocate( cell_tmp(ic)%dam( nelemIn ), stat=istat, errmsg=msg )
               write( text, '(a,i0,a)' )  'cell_tmp(', ic, ')%dam'
               call check_allocation( 'relaxCell1dNU', text, istat, msg )
               cell_tmp(ic)%dam(:) = ZERO
            end do
         end if
      end if

      do ic = 1, grid%ncell
         cell(ic)%edam = edamnew(ic)
      end do

      cell(:)%lxtal = lxtalnew(:)
!
!  Deallocate arrays no longer needed
!
      deallocate( lxtalnew, wnew, densnew, edamnew, vibnew )
      if( out%lhis )  deallocate( damnew )
!
!  Reinitialize the surface
!
      call updateSurface1d( grid%posif(0), grid%posif(grid%ncell) )
!          ===============        i              i
!
!  Calculate derived quantities
!
      call initDerivedCell1dNU( )
!          ===================

   end subroutine relaxCell1dNU
!
!=======================================================================
!
   subroutine writeCell1dNU( )
!             =============
!
!     Write the contents of the target cells to the CELL file (same format as HIS file)
!
!     NOTE: In case of surface roughness the written cell densities are the densities before
!           their reduction due to roughness
!
!  Arguments:
!
!     none
!
!  Local Variables:
!
      integer(INT4) :: ic, ielem

      character(6), save :: fmode = 'write'
!
!  Start of Executable Code:
!  ------------------------
!
      call openFile( 22, 'CELL', 'formatted', fmode )
!          ========  i     i          i         i
!
!  Write header
!
      write( 22, * )  1+nelem, 2*grid%ncell+2
      write( 22, * )  0, elem(1:nelem)%z
!
!  Write first point
!
      write( 22, '(f11.3, t13, 5es13.5)' )  &
           grid%posif(0), (1e24*HISMIN*ion%ddose/grid%wcell(1), ielem=1,nelem)
!
!  Write density values
!
      do ic = 1, grid%ncell
         write( 22, '(f11.3, t13, 5es13.5)' )  &
            grid%posif(ic-1),  &
            ( 1e24 * max( HISMIN*ion%ddose/grid%wcell(ic), cell(ic)%dens(ielem) ),  &
              ielem=1,nelem )
         write( 22, '(f11.3, t13, 5es13.5)' )  &
            grid%posif(ic),  &
            ( 1e24 * max( HISMIN*ion%ddose/grid%wcell(ic), cell(ic)%dens(ielem) ),  &
              ielem=1,nelem )
      end do
!
!  Write last point
!
      write( 22, '(f11.3, t13, 5es13.5)' )  &
           grid%posif(grid%ncell), (1e24*HISMIN*ion%ddose/grid%wcell(1), ielem=1,nelem)

      close( 22 )

      fmode = 'append'

   end subroutine writeCell1dNU
!
!=======================================================================
!
   subroutine writeDamCell1dNU( )
!             ================
!
!     Write the damage of the target cells to the HIS file
!
!  Arguments:
!
!     none
!
!  Local Variables:
!
      integer(INT4) :: ic, ielem
!
!  Start of Executable Code:
!  ------------------------
!
      call openFile( 22, 'HIS', 'formatted', 'write' )
!          ========  i     i         i          i
!
!  Write header
!
      write( 22, * )  1+nelemIn, 2*grid%ncell+2
      write( 22, * )  0, elem(1:nelem)%z, -elem(1:nelem)%z
!
!  Write first point
!
      write( 22, '(f11.3, t13, 100es13.6)' )  &
           grid%posif(0), (1e24*HISMIN*ion%ddose/grid%wcell(1), ielem=1,nelemIn)
!
!  Write histogram values
!
      do ic = 1, grid%ncell
         write( 22, '(f11.3, t13, 100es13.6)' )  &
            grid%posif(ic-1),  &
            ( 1e24 * max( HISMIN*ion%ddose/grid%wcell(ic), cell(ic)%dam(ielem) ),  &
              ielem=1,nelemIn )
         write( 22, '(f11.3, t13, 100es13.6)' )  &
            grid%posif(ic),  &
            ( 1e24 * max( HISMIN*ion%ddose/grid%wcell(ic), cell(ic)%dam(ielem) ),  &
              ielem=1,nelemIn )
      end do
!
!  Write last point
!
      write( 22, '(f11.3, t13, 100es13.6)' )  &
           grid%posif(grid%ncell), (HISMIN*ion%ddose/grid%wcell(1), ielem=1,nelemIn)

      close( 22 )

   end subroutine writeDamCell1dNU
!
!=======================================================================
!
   subroutine readAtomsFromCell1dNU( )
!             =====================
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
      integer(INT4) :: idummy
      integer(INT4) :: i, ielem
      real(REAL8), allocatable :: z(:)
!
!  Start of Executable Code:
!  ------------------------
!
!  Open file
!
      call openFile( 22, cellfile, 'formatted', 'read' )
!          ========  i      i           i         i
!
!  First line can be ignored
!
      read( 22, *, err=901, end=901 )  ( idummy, i=1,2 )
!
!  Read the atomic numbers and determine the atom masses
!
      if( ldyn ) then
         allocate( z( nelem ) )
         read( 22, *, err=900, end=900 ) idummy, z(:)
         do ielem = 1, nelemIon
            if( z(ielem) /= elem(ielem)%z )  &
               error stop 'readAtomsFromCell1dNU: Inconsistent ion atomic number'
         end do
         do ielem = nelemIon+1, nelem
            elem(ielem)%z = z(ielem)
         end do
         deallocate( z )
      else
         read( 22, *, err=900, end=900 ) idummy, elem(nelemIon+1:nelem)%z
      end if

      close( 22 )

      return

900   error stop 'readAtomsFromCell1dNU: error during reading CELL file'
901   error stop 'readAtomsFromCell1dNU: no data on CELL file'

   end subroutine readAtomsFromCell1dNU
!
!=======================================================================
!
   subroutine readCell1dNU( )
!             ============
!
!     Read the contents of the target cells from the CELL file cellfile.
!     The CELL file has the same format as the HIS file
!
!     TODO: call readCell1dNUraw
!
!  Arguments:
!
!     None
!
!  Local Variables:
!
      character(80) :: text
      character(256) :: msg
      integer(INT4) :: ncol, nrow, idummy, irec, nrec
      integer(INT4) :: ic, istat, i, ielem1
      real(REAL8) :: rdummy
!
!  Start of Executable Code:
!  ------------------------
!
!  Open file
!
      call openFile( 22, cellfile, 'formatted', 'read' )
!          ========  i      i           i         i

      nrec = cellrec
      if( nrec == 0 )  nrec = huge(0_INT4)
      grid%ncell = 0

      do irec = 1, nrec
!
!  Read first line and allocate memory
!
         read( 22, *, err=800, end=800 )  ncol, nrow

         if( ncol-1 /= nelemTarg )  &
               error stop 'readCell1dNU: inconsistent number of columns on file'

         if( (nrow-2)/2 /= grid%ncell ) then

            if( allocated( grid%posif ) )  deallocate( grid%posif )

            if( allocated( cell ) ) then
               do ic = 1, grid%ncell
                  deallocate( cell(ic)%dens )
               end do
               deallocate( cell )
            end if

            grid%ncell = ( nrow - 2 ) / 2

            allocate( grid%posif( 0:grid%ncell ), stat=istat, errmsg=msg )
            call check_allocation( 'readCell1dNU', 'grid%posif', istat, msg )

            allocate( cell( grid%ncell ), stat=istat, errmsg=msg )
            call check_allocation( 'readCell1dNU', 'cell', istat, msg )

            do ic = 1, grid%ncell
               allocate( cell(ic)%dens( nelem ), stat=istat, errmsg=msg )
               write( text, '(a,i0,a)' )  'cell(', ic, ')%dens'
               call check_allocation( 'readCell1dNU', text, istat, msg )
            end do

         end if
!
!  Second and third line can be ignored
!
         read( 22, *, err=900, end=900 )  ( idummy, i=1,ncol )
         read( 22, *, err=900, end=900 )  ( rdummy, i=1,ncol )
!
!  Read cell densities
!
         if( ldyn ) then
            ielem1 = 1
         else
            do ic = 1, grid%ncell
               cell(ic)%dens(1:nelemIon) = ZERO
            end do
            ielem1 = nelemIon+1
         end if

         do ic = 1, grid%ncell
            read( 22, *, err=900, end=900 )  grid%posif(ic-1), cell(ic)%dens(ielem1:nelem)
            read( 22, *, err=900, end=900 )  grid%posif(ic), cell(ic)%dens(ielem1:nelem)
            cell(ic)%dens(ielem1:nelem) = cell(ic)%dens(ielem1:nelem) / 1e24
         end do
!
!  Ignore last line
!
         read( 22, *, err=900, end=900 )  ( rdummy, i=1,ncol )

      end do

800   if( irec == 1 )  goto 901

      close( 22 )
!
!  Shift and reflect the distribution if desired
!
      if( defined( shiftto ) )  &
         grid%posif(:) = shiftto + grid%posif(:) - grid%posif(0)

      if( reflect ) then
         grid%posif(:) = - grid%posif(grid%ncell:0:-1)
         cell(:) = cell(grid%ncell:1:-1)
      end if
!
!  Cell widths
!
      allocate( grid%wcell( grid%ncell ) )
      grid%wcell(:) = grid%posif(1:grid%ncell) - grid%posif(0:grid%ncell-1)
!
!  CELL file does not contain info on crystallinity, therefore set to amorphous
!
      cell(:)%lxtal = .false.

      return
!
900   print *, 'ic=', ic
      error stop 'readCell1dNU: error during reading CELL file'
901   error stop 'readCell1dNU: no data on CELL file'

   end subroutine readCell1dNU
!
!=======================================================================
!
   subroutine readBndFromCell1dNU( zmin, zmax )
!             ===================   o     o
!
!     Read the contents of the target cells from the CELL file cellfile, and determine
!     the geometry boundaries (the positions of the surfaces).
!
!  Arguments:
!
      real(REAL8), intent(out) :: zmin          ! position of top surface
      real(REAL8), intent(out) :: zmax          ! position of bottom surface
!
!  Local Variables:
!
      integer(INT4) :: iz, jz, nz
      integer(INT4), allocatable :: indcell(:)
      real(REAL8) :: ff
      real(REAL8), allocatable :: z(:), dens(:,:)
!
!  Start of Executable Code:
!  ------------------------
!
!  Read file record
!
      call readCell1dNUraw( indcell, z, dens )
!          ===============     o     o   o
!
!  Search for first non-empty cell
!
      nz = ubound( z, 1 )

      do iz = 1, nz
         if( any( dens(1:nelem,iz) > ZERO ) )  exit
      end do
      jz = iz

      ff = sum( dens(1:nelem,jz)/1e24 / elem(:)%dens )
      ff = max( ZERO, min( ONE, ff ) )
      zmin = ff*z(jz-1) + (ONE-ff)*z(jz)
!
!  Search for last non-empty cell
!
      do iz = jz+1, nz
         if( any( dens(1:nelem,iz) == ZERO ) )  exit
      end do
      jz = iz - 1

      ff = sum( dens(1:nelem,jz)/1e24 / elem(:)%dens )
      ff = max( ZERO, min( ONE, ff ) )
      zmax = ff*z(jz) + (ONE-ff)*z(jz-1)

   end subroutine readBndFromCell1dNU
!
!=======================================================================
!
   subroutine readCell1dNUraw( indcell, z, data )
!             ===============     o     o   o
!
!     Read the contents of the target cells from the CELL file cellfile.
!     The CELL file has the same format as the HIS file.
!
!     Version returning arrays
!
!  Arguments:
!
      integer(INT4), allocatable, intent(inout) :: indcell(:)  ! (nspec) indicator of column
                                                               !     meaning
      real(REAL8), allocatable, intent(inout) :: z(:)          ! (0:nz) bin boundaries
      real(REAL8), allocatable, intent(inout) :: data(:,:)     ! (nspec,nz) (mostly) densities
!
!  Local Variables:
!
      integer(INT4) :: ncol, nrow, nspec, nz, idummy, irec, nrec
      integer(INT4) :: iz, i
      real(REAL8) :: rdummy
!
!  Start of Executable Code:
!  ------------------------
!
!  Open file
!
      call openFile( 22, cellfile, 'formatted', 'read' )
!          ========  i      i           i         i

      nrec = cellrec
      if( nrec == 0 )  nrec = huge(0_INT4)
      nz = 0
      nspec = 0

      do irec = 1, nrec
!
!  Read first line and allocate memory
!
         read( 22, *, err=800, end=800 )  ncol, nrow

         if( nspec /= ncol-1 ) then
            call alloc( indcell, ncol-1, 'indcell', 'readCell1dNUraw' )
!                =====     o        i        i              i
         end if

         if( nz /= (nrow-2)/2 ) then
            if( allocated( z ) )  deallocate( z )
            call alloc( z, 0, (nrow-2)/2, 'z', 'readCell1dNUraw' )
!                =====  o  i       i       i           i
         end if

         if( nspec /= ncol-1 .or. nz /= (nrow-2)/2 ) then
            if( allocated( data ) )  deallocate( data )
            call alloc( data, [ ncol-1, (nrow-2)/2 ], 'data', 'readCell1dNUraw' )
!                =====   o         i         i          i             i
         end if

         nspec = ncol - 1
         nz  = ( nrow - 2 ) / 2
!
!  Second has the indicator
!
         read( 22, *, err=900, end=900 )  idummy, indcell(:)
!
!  Ignore third line
!
         read( 22, *, err=900, end=900 )  ( rdummy, i=1,ncol )
!
!  Read cell densities
!
         do iz = 1, nz
            read( 22, *, err=900, end=900 )  z(iz-1), data(:,iz)
            read( 22, *, err=900, end=900 )  z(iz), data(:,iz)
         end do
!
!  Ignore last line
!
         read( 22, *, err=900, end=900 )  ( rdummy, i=1,ncol )

      end do

800   if( irec == 1 )  goto 901

      close( 22 )

      return
!
900   print *, 'iz=', iz
      error stop 'readCell1dNU: error during reading CELL file'
901   error stop 'readCell1dNU: no data on CELL file'

   end subroutine readCell1dNUraw
!
!=======================================================================
!
   subroutine initCell1dNUFromGeom( )
!             ====================
!
!     Initialize the cell structure of the target from a 1d region definition
!
!     The initial grid is regular. At region boundaries, the cell densities are averages over
!     the regions covered by the cell.
!
!  Arguments:
!
!     None
!
!  Local Variables:
!
      character(80) :: text
      character(256) :: msg
      integer(INT4) :: ic, icmin, icmax, ir, istat
      integer(INT4) :: ielem, ielemr

      real(REAL8) :: offset, zminr, zmaxr, wcell
      real(REAL8) :: ddens, shift
!
!  Start of Executable Code:
!  ------------------------
!
!  Determine grid from target geometry (geom%posif)
!
      if( irxtal > 0 ) then
         offset = -100 * EPSMA * geom%posif(nr)  ! avoid cell boundaries at lattice planes
      else
         offset = ZERO
      end if

      if( geom%nbox == 0 ) then

         wcell = geom%wbox

         icmin = floor( ( geom%posif(0) - offset ) / wcell )
         icmax = ceiling( ( geom%posif(nr) - offset ) / wcell )
         grid%ncell = icmax - icmin

         call alloc( grid%posif, 0, grid%ncell, 'grid%posif', 'initCell1dNUFromGeom' )
!             =====      o       i       i           i                  i
         grid%posif(:) = [ (offset + ic * wcell, ic=icmin,icmax) ]

      else if( geom%nbox == 1 ) then             ! probably not useful

         error stop 'initCell1dNUFromGeom: requires geom%nbox /= 1'

      else

         grid%ncell = geom%nbox

         wcell = (geom%posif(nr) - geom%posif(0)) / (grid%ncell - 1)

         allocate( grid%posif( 0:grid%ncell ) , stat=istat, errmsg=msg )
         call check_allocation( 'initCell1dNUFromGeom', 'grid%posif', istat, msg )
         grid%posif(:) = [ (geom%posif(0) + ic * wcell, ic=0,grid%ncell) ]

         shift = modulo( grid%posif(0)-offset, wcell )
         grid%posif(:) = grid%posif(:) - shift

      end if

      grid%posif(0) = geom%posif(0)
      grid%posif(grid%ncell) = geom%posif(nr)
!
!  Cell widths
!
      allocate( grid%wcell( grid%ncell ) )
      grid%wcell(:) = grid%posif(1:grid%ncell) - grid%posif(0:grid%ncell-1)
!
!  Cell densities
!
      allocate( cell( grid%ncell ), stat=istat, errmsg=msg )
      call check_allocation( 'initCell1dNU', 'cell', istat, msg )

      do ic = 1, grid%ncell
         allocate( cell(ic)%dens( nelem ), stat=istat, errmsg=msg )
         write( text, '(a,i0,a)' )  'cell(', ic, ')%dens'
         call check_allocation( 'initCell1dNU', text, istat, msg )
         cell(ic)%dens(:) = ZERO
      end do

      do ir = 1, nr
         zminr = ( geom%posif(ir-1) - grid%posif(0) ) / geom%wbox
         zmaxr = ( geom%posif(ir) - grid%posif(0) ) / geom%wbox

         icmin = max( 1, ceiling( zminr ) )
         icmax = max( 1, ceiling( zmaxr ) )

         do ielemr = 1, mat(ir)%nelemr
            ielem = mat(ir)%ielem(ielemr)
            ddens = mat(ir)%dens * mat(ir)%prob(ielemr)

            do ic = icmin, icmax
               cell(ic)%dens(ielem) = cell(ic)%dens(ielem) + ddens
            end do
!
            cell(icmin)%dens(ielem) = cell(icmin)%dens(ielem) - ddens * ( zminr - icmin + 1 )
            cell(icmax)%dens(ielem) = cell(icmax)%dens(ielem) - ddens * ( icmax - zmaxr )
         end do
      end do
!
!  Crystal flag
!
      cell(:)%lxtal = .true.

      do ir = 1, nr
         if( mat(ir)%lxtal )  cycle
         call getCellIndex1dNU_core( geom%posif(ir-1), icmin )
!             =====================         i            o
         call getCellIndex1dNU_core( geom%posif(ir), icmax )
!             =====================         i          o

         cell(icmin:icmax)%lxtal = .false.
      end do

   end subroutine initCell1dNUFromGeom
!
!=======================================================================
!
   subroutine readHis1dNU( hisfile )
!             ===========    i
!
!     Read the damage concentrations from the hisfile file (HIS file format)
!
!  Arguments:
!
      character(*), intent(in) :: hisfile       ! file to be read from
!
!  Local Variables:
!
      integer(INT4) :: ncol, nrow, idummy, ncell
      integer(INT4) :: ielem, ic, i
      real(REAL8) :: rdummy, posif, posif1, posif2, dam(nelemIn)
!
!  Start of Executable Code:
!  ------------------------
!
!  Open file
!
      call openFile( 22, hisfile, 'formatted', 'read' )
!          ========  i     i         i          i
!
!  Read first line and allocate memory
!
      read( 22, *, err=900, end=900 ) ncol, nrow

      if( ncol-1 /= nelemIn )  error stop 'readHis1dNU: inconsistent number of columns on file'

      ncell = ( nrow - 2 ) / 2
!
!  Second and third line can be ignored
!
      read( 22, *, err=900, end=900 )  ( idummy, i=1,ncol )
      read( 22, *, err=900, end=900 )  ( rdummy, i=1,ncol )
!
!  Read histogram values
!
      do ic = 1, min( ncell, grid%ncell )
         read( 22, *, err=900, end=900 )  posif1, dam(:)
         read( 22, *, err=900, end=900 )  posif2, dam(:)
         posif = HALF * ( posif1 + posif2 )
         dam(:) = dam(:) * (posif2 - posif1) / ion%ddose
         do ielem = 1, nelem
            call scoreCell1dNU( ielem, ic, ZERO, dam(ielem)/1e24, dam(ielem+nelem)/1e24, ZERO )
!                =============    i    i    i          i                    i             i
         end do
      end do
!
!  Ignore last line
!
      close( 22 )
      return

900   error stop 'readHis1dNU: error during reading histogram file'

   end subroutine readHis1dNU
!
!=======================================================================
!
   subroutine initDerivedCell1dNU( )
!             ===================
!
!     Initialize derived quantities of the cell structure of the 1d target
!
!     OpenMP: This suboutine must be called from within a parallel region.
!
!  Arguments:
!
!     none
!
!  Local Variables:
!
      character(80) :: text
      character(256) :: msg
      integer(INT4) :: ic, istat, ielem
!
!  Start of Executable Code:
!  ------------------------
!
!  Smallest and largest position safely inside/outside target
!
      if( grid%posif(0) == ZERO ) then
         grid%posmin1 = -SMALL
         grid%posmin2 = +SMALL
      else
         grid%posmin1 = grid%posif(0) - EPSMA * abs(grid%posif(0))
         grid%posmin2 = grid%posif(0) + EPSMA * abs(grid%posif(0))
      end if

      if( grid%posif(grid%ncell) == ZERO ) then
         grid%posmax1 = -SMALL
         grid%posmax2 = +SMALL
      else
         grid%posmax1 = grid%posif(grid%ncell) - EPSMA*abs(grid%posif(grid%ncell))
         grid%posmax2 = grid%posif(grid%ncell) + EPSMA*abs(grid%posif(grid%ncell))
      end if
!
!  Derived quantities
!
      do ic = 1, grid%ncell
         cell(ic)%dens_tot = sum( cell(ic)%dens(:) )
      end do

      do ic = 1, grid%ncell
         allocate( cell(ic)%prob( nelem ), stat=istat, errmsg=msg )
         write( text, '(a,i0,a)' )  'cell(', ic, ')%prob'
         call check_allocation( 'initDerivedCell1dNU', text, istat, msg )
         if( cell(ic)%dens_tot == ZERO ) then
            cell(ic)%prob(:) = ZERO
         else
            cell(ic)%prob(:) = cell(ic)%dens(:) / cell(ic)%dens_tot
         end if
      end do

      do ic = 1, grid%ncell
         allocate( cell(ic)%vib( nelem ), stat=istat, errmsg=msg )
         write( text, '(a,i0,a)' )  'cell(', ic, ')%vib'
         call check_allocation( 'initDerivedCell1dNU', text, istat, msg )
         cell(ic)%vib(:) = elem(:)%vib
      end do

      do ic = 1, grid%ncell
         allocate( cell(ic)%ed( nelem ), stat=istat, errmsg=msg )
         write( text, '(a,i0,a)' )  'cell(', ic, ')%ed'
         call check_allocation( 'initDerivedCell1dNU', text, istat, msg )
         do ielem = 1, nelem
            cell(ic)%ed(ielem) = sum( cell(ic)%prob(:) * dam%ed(:,ielem) )
         end do
      end do

      do ic = 1, grid%ncell
         allocate( cell(ic)%esurf( nelem ), stat=istat, errmsg=msg )
         write( text, '(a,i0,a)' )  'cell(', ic, ')%esurf'
         call check_allocation( 'initDerivedCell1dNU', text, istat, msg )
         do ielem = 1, nelem
            cell(ic)%esurf(ielem) = sum( cell(ic)%prob(:) * dam%esurf(:,ielem) )
         end do
      end do

      do ic = 1, grid%ncell
         cell(ic)%edam = ZERO
      end do

   end subroutine initDerivedCell1dNU
!
!=======================================================================
!
   subroutine getCellIndex1dNU_core( z, ic, outside_ )
!             =====================  i   o      o
!
!     Determine the cell the point z is located in or closest to.
!
!  Arguments:
!
      real(REAL8), intent(in) :: z              ! position to be checked
      integer(INT4), intent(out) :: ic          ! index of cell (region)
      logical, optional, intent(out) :: outside_! flag indicationg whether pos is outside any
                                                !    region
!
!  Local Variables:
!
!     None
!
!  Start of Executable Code:
!  ------------------------
!
!  Determine cell pos is located in (positions outside target map to ic=1 or
!  ic=ncell)
!
      call binSearch( grid%ncell+1, grid%posif, z, ic )
!          =========       i             i      i  o

      if( present( outside_ ) )  &
            outside_ = ( z < grid%posif(0) .or. z > grid%posif(grid%ncell) )

      ic = max( 1, min( grid%ncell, ic ) )

   end subroutine getCellIndex1dNU_core
!
!=======================================================================
!
   subroutine diffuseCell1dNU( posif, dens, dpa, ielem )
!             ===============    i     io    i     i
!
!     Diffuse dens(posif) with diffusivity induced by dpa
!
!  Arguments:
!
      real(REAL8), intent(in) :: posif(0:)   ! boundaries of histogram boxes
      real(REAL8), intent(inout) :: dens(:)  ! concentrations in histogram
      real(REAL8), intent(in) :: dpa(:)      ! dpa's
      integer(INT4), intent(in) :: ielem     ! element index of diffusing species
!
!  Local Variables:
!
      character(256) :: msg
      integer(INT4) :: ncell, istat, i, niter
      real(REAL8) :: fac
      real(REAL8), allocatable :: dz(:), dz_half(:), dpa_half(:), a(:), b(:), ddens(:)
!
!  Start of Executable Code:
!  ------------------------
!
      ncell = size( dens )
      if( size( posif ) /= ncell+1 .or. size( dpa ) /= ncell )  goto 900
      allocate( dz( ncell ), stat=istat, errmsg=msg )
      call check_allocation( 'diffuseCell1dNU', 'dz', istat, msg )
      allocate( dz_half( 0:ncell ), stat=istat, errmsg=msg )
      call check_allocation( 'diffuseCell1dNU', 'dz_half', istat, msg )
      allocate( dpa_half( 0:ncell ), stat=istat, errmsg=msg )
      call check_allocation( 'diffuseCell1dNU', 'dpa_half', istat, msg )
      allocate( a( ncell ), stat=istat, errmsg=msg )
      call check_allocation( 'diffuseCell1dNU', 'a', istat, msg )
      allocate( b( ncell ), stat=istat, errmsg=msg )
      call check_allocation( 'diffuseCell1dNU', 'b', istat, msg )
      allocate( ddens( ncell ), stat=istat, errmsg=msg )
      call check_allocation( 'diffuseCell1dNU', 'ddens', istat, msg )

      dz(:) = posif(1:ncell) - posif(0:ncell-1)
      dz_half(0) = dz(1)
      dz_half(1:ncell-1) = HALF * (dz(1:ncell-1) + dz(2:ncell))
      dz_half(ncell) = dz(ncell)
      dpa_half(0) = dpa(1)
!      dpa_half(1:ncell-1) = HALF * (dpa(1:ncell-1) + dpa(2:ncell))
      dpa_half(1:ncell-1) = (dpa(1:ncell-1)*dz(1:ncell-1) + dpa(2:ncell)*dz(2:ncell)) /  &
                            (TWO * dz_half(1:ncell-1))
      dpa_half(ncell) = dpa(ncell)

      fac = dam%diffusivity(ielem)
      a(:) = fac * dpa_half(0:ncell-1) / dz_half(0:ncell-1) / dz(:)
      b(:) = fac * dpa_half(1:ncell) / dz_half(1:ncell) / dz(:)

      niter = ceiling( 1.25 * maxval( a(:) + b(:) ) ) ! iteration count to guarantee stability
                                                      !    of explicit time integration
      print *, 'niter=', niter
      print *, 'maxval(dpa)=', maxval(dpa)
      a(:) = a(:) / niter
      b(:) = b(:) / niter

      do i = 1, niter
         ddens(:) = - (a(:) + b(:)) * dens(:)
         ddens(2:ncell) = ddens(2:ncell) + a(2:ncell) * dens(1:ncell-1)
         ddens(1:ncell-1) = ddens(1:ncell-1) + b(1:ncell-1) * dens(2:ncell)

         dens(:) = dens(:) + ddens(:)
      end do

      deallocate( dz, dz_half, dpa_half, a, b, ddens )

      return

900   error stop 'diffuseCell1dNU: incompatible dimensions'

   end subroutine diffuseCell1dNU
!
!=======================================================================
!
end module Target_Cell1d_NU_mod
