module Target_Cell2d_Old_mod
!
!  Holds old 2d cell structure and provides subroutines operating on it.
!  All grid related features are handled by calling subroutines from Target_Grid2d_mod
!
!  Note: avoid ldyn=t with lspher=t
!
!  Used Modules:
!  ============
!
   use Const_mod,                only: INT4, REAL8,  &
                                       ZERO, QUART, HALF, ONE, TWO, SIX,  &
                                       BIG, EPSMA, SMALL,  &
                                       ENUM
   use Array_mod,                only: check_allocation

   use Open_mod,                 only: openFile

   use Param_Setup_mod,          only: nr, nelem, nelemIon, nelemTarg, niongeom, nionout,  &
                                       ldyn, usecell, cellfile, cellrec, nion
   use Param_Geom_mod,           only: geom
   use Param_Elem_mod,           only: elem
   use Param_Ion_mod,            only: ion
   use Param_Mat_mod,            only: mat
   use Param_Dam_mod,            only: dam
   use Param_Out_mod,            only: out

   use Types_mod,                only: Proj_type
   use Target_Surface2d_SDF_mod, only: writeSurface2d
   use Target_Grid2d_mod,        only: getCellIndexGrid2d, getOtherNearCellIndicesGrid2d,  &
                                       initGrid2d, readGrid2d, writeGrid2d,  &
                                       initPostCellGrid2d, relaxGrid2d,  &
                                       getOvlpAreasGrid2d, getAreasGrid2d
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
   public :: initCell2dOld                ! initialize the target
!            =============
   public :: scoreCell2dOld               ! score a count in the target cell
!            ==============
   public :: getCellIndex2dOld            ! determines the cell a given point is located in
!            =================
   public :: getOtherNearCellIndices2dOld ! determine cells other than nearest cell that
!            ============================ !    intersect with the free flight path cylinder
   public :: getProbCell2dOld             ! get the probabilities to encounter each of the atoms
!            ================             !    in a given cell
   public :: getDensCell2dOld             ! get the total atomic density in a given cell
!            ================
   public :: getVibCell2dOld              ! get the array of thermal vibration amplitudes
!            ===============              !    in a given cell
   public :: getXtalCell2dOld             ! get the crystallinity flag of the cell
!            ================
   public :: getEdCell2dOld               ! get the displacement energy after most recent
!            ==============               !    target cell update
   public :: getEsurfCell2dOld            ! get the surface binding energy after most recent
!            =================            !    target cell update
   public :: writeCell2dOld               ! write target cell contents to output file
!            ==============
   public :: relaxCell2dOld               ! relax the cell densities
!            ==============
   public :: readAtomsFromCell2dOld       ! read target cell contents from file
!            ======================
   public :: readCell2dOld                ! read target cell contents from file
!            =============
!
!  Private Procedures:
!  ==================
!
   private :: initDerivedCell2dOld     ! initialize derived cell properties
!             ====================
!
!  Private Data:
!  ============
!
   type :: NearCell_type                     ! cell with distmin2 < minval(distmax2)
!          =============
      integer(INT4) :: ix                    ! cell index in x-direction
      integer(INT4) :: iz                    ! cell index in z-direction
      real(REAL8) :: distmin2                ! square of minimum minimum distance of
                                             !    point inside the cell from boundary of
                                             !    region ir
      real(REAL8) :: distmax2                ! square of maximum minimum distance of
                                             !    point inside the cell from boundary of
                                             !    region ir
   end type NearCell_type

   type :: Cell_type
!          =========
      real(REAL8) :: area                       ! cell area
      real(REAL8), allocatable :: dens(:)       ! (nelem) density of atom ielem in cell
      real(REAL8) :: dens_tot                   ! total atom density in cell
      real(REAL8), allocatable :: prob(:)       ! (nelem) probability to encounter atom ielem
                                                !    in cell
      real(REAL8) :: ff                         ! fill fraction
      logical :: empty                          ! indicates empty cell at last update
      real(REAL8), allocatable :: vib(:)        ! (nelem) thermal vibration amplitude
      logical :: lxtal                          ! flag indicating crystallinity
      real(REAL8), allocatable :: ed(:)         ! (nelem) displacement energy after
                                                !    most recent target cell update
      real(REAL8), allocatable :: esurf(:)      ! (nelem) surface binding energy after
                                                !    most recent target cell update
      integer(INT4), allocatable :: icnear(:)   ! indices of near nonempty cells, ordered
                                                !    according to increasing distmin2
      real(REAL8) :: stress(13)                 ! pressure and deviatoric stress tensor elements
      real(REAL8) :: edam                       ! damage energy per cm3 since last relax
   end type Cell_type

   type(Cell_type), allocatable :: cell(:)   ! (ncell) cells

   integer(INT4) :: ncell                    ! number of cells
!
!==========================================================================================
!
contains
!
!==========================================================================================
!
   subroutine initCell2dOld( )
!             =============
!
!     Initialize the cell structure of the target
!
!  Arguments:
!
!  Local Variables:
!
      character(80) :: text
      character(256) :: errmsg
      integer(INT4) :: ic, istat, ir, np, ielemr, ielem, iarea
      real(REAL8) :: ddens, dvib(nelem)
      integer(INT4), allocatable :: ic_area(:)
      real(REAL8), allocatable :: area(:)
      logical, allocatable :: set(:)
!
!  Start of Executable Code:
!  ------------------------
!
      if( usecell ) then
         call initDerivedCell2dOld( )
!             ====================
         return                                 ! Cells have already been read
      end if
!
!  Initialize the grid
!
      call initGrid2d( ncell )
!          ==========    o
!
!  Allocate and initialize cell contents
!
      allocate( cell( ncell ), stat=istat, errmsg=errmsg )
      call check_allocation( 'initCell2dOld', 'cell', istat, errmsg )

      do ic = 1, ncell
                                                ! Cell densities
         allocate( cell(ic)%dens( nelem ), stat=istat, errmsg=errmsg )
         write( text, '(a,i0,a)' )  'cell(', ic, ')%dens'
         call check_allocation( 'initCell2dOld', text, istat, errmsg )
         cell(ic)%dens(:) = SMALL
                                             ! Thermal vibration amplitudes
         allocate( cell(ic)%vib( nelem ), stat=istat, errmsg=errmsg )
         write( text, '(a,i0,a)' )  'cell(', ic, ')%vib'
         call check_allocation( 'initCell2dOld', text, istat, errmsg )
         cell(ic)%vib(:) = ZERO
                                             ! Crystal flag
         cell(ic)%lxtal = .false.            !    amorphous simulation can always be done

         cell(ic)%stress(:) = ZERO

         cell(ic)%edam = ZERO                ! THIS SEEMS TO BE DONE TWICE

      end do
!
!  Get geometric properties of the cells
!
      call getAreasGrid2d( cell(:)%area )
!          ==============       o
!
!  Set cell contents using the region geometries
!
      allocate( set( ncell ), stat=istat, errmsg=errmsg )
      call check_allocation( 'initCell2dOld', 'set', istat, errmsg )
      set(:) = .false.                             ! remember which cells have set lxtal value

      do ir = 1, nr
         np = geom%region(ir)%np + 1
         call getOvlpAreasGrid2d( geom%region(ir)%px(1:np), geom%region(ir)%pz(1:np), area,  &
!             ==================             i                             i            o
                                  ic_area)
!                                    o

         do ielemr = 1, mat(ir)%nelemr
            ielem = mat(ir)%ielem(ielemr)
            ddens = mat(ir)%dens * mat(ir)%prob(ielemr)
            do iarea = 1, size( area )
               ic = ic_area(iarea)
               cell(ic)%dens(ielem) = cell(ic)%dens(ielem) + area(iarea)/cell(ic)%area * ddens
            end do
         end do

         dvib(:) = mat(ir)%vib(:)
         do iarea = 1, size( area )
            ic = ic_area(iarea)
            cell(ic)%vib(:) = cell(ic)%vib(:) + area(iarea)/cell(ic)%area * dvib(:)
         end do

         if( mat(ir)%lxtal ) then
            do iarea = 1, size( area )
               ic = ic_area(iarea)
               if( area(iarea) > ZERO .and. .not.set(ic) ) then
                  cell(ic)%lxtal = .true.
                  set(ic) = .true.
               end if
            end do
         else
            do iarea = 1, size( area )
               ic = ic_area(iarea)
               if( area(iarea) > ZERO ) then
                  cell(ic)%lxtal = .false.
                  set(ic) = .true.
               end if
            end do
         end if
         deallocate( area, ic_area )
      end do

      deallocate( set )

!      open( 111, file='DUMP_OUTPUT' )
!      do ielem=1,nelem
!         write( 111, '(a,i3)' ) 'ielem=',ielem
!         do iz=1,grid%nz
!            write( 111, '(a,i3)', advance='no' ) 'iz=',iz
!            write( 111, '(100e10.2)' ) (cell(ix,iz)%dens(ielem), ix=1,grid%nx)
!         end do
!      end do
!      close( 111 )
!
!  Calculate derived quantities
!
      call initDerivedCell2dOld( )
!          ====================
!
!  Test
!
!      do ic = 1, ncell
!         cell(ic)%dens(:) = cell(ic)%dens * 0.5
!      end do
!      call relaxCell2dOld( )
!          ==============
!      error stop
!
   end subroutine initCell2dOld
!
!===================================================================================
!
   subroutine scoreCell2dOld( ielem, ic, w, edam )
!             ==============    i    i   i   i
!
!     Score a count in the cell histogram
!
!     Cell contents are forced not get less than zero.
!
!  Arguments:
!
      integer(INT4), intent(in) :: ielem        ! index of atom species
      integer(INT4), intent(in) :: ic           ! cell index
      real(REAL8), intent(in) :: w              ! weight (can also be negative)
      real(REAL8), intent(in) :: edam           ! damage energy (eV)
!
!  Local Variables:
!
!     None
!
!  Start of Executable Code:
!  ------------------------
!
!
!  Cell contents
!
      if( w /= ZERO ) then
         cell(ic)%dens(ielem) = max( SMALL, cell(ic)%dens(ielem) +  &   ! SMALL is used to
                                            w * ion%ddose * ion%dxinit / cell(ic)%area )
                                                                        ! check on empty cells
         cell(ic)%dens_tot = sum( cell(ic)%dens(:) )

         cell(ic)%prob(:) = cell(ic)%dens(:) / cell(ic)%dens_tot
      end if
!
!  Damage energy per cell volume
!
      if( edam /= ZERO ) then
         cell(ic)%edam = cell(ic)%edam +  &
              edam * ion%ddose * ion%dxinit / cell(ic)%area
      end if
!
!  Do not update ed and esurf to avoid unreasonable values when cell gets empty!
!
   end subroutine scoreCell2dOld
!
!=======================================================================
!
   subroutine getCellIndex2dOld( pos, ic, outside_ )
!             =================   i   o      o
!
!     Determine non-empty cell the point pos is located in or closest to.
!
!  Arguments:
!
      real(REAL8), intent(in) :: pos(3)         ! position to be checked
      integer(INT4), intent(out) :: ic          ! 1D index of cell (region); if outside all
                                                !    cells, ic=0
      logical, optional, intent(out) :: outside_! flag indication whether pos is outside all
                                                !    cells with content
!
!  Local Variables:
!
!     None
!
!  Start of Executable Code:
!  ------------------------
!
      call getCellIndexGrid2d( pos(1), pos(3), ic, outside_ )
!          ==================    i       i     o      o

   end subroutine getCellIndex2dOld
!
!==========================================================================================
!
   subroutine getOtherNearCellIndices2dOld( proj, posc, ffp, pmax, ncn, icnear )
!             ============================   i     i     i    i     o     o
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
                                                ! index of regions found except proj%id
!
!  Local Variables:
!
      integer(INT4) :: in, nn, ic
!
!  Start of Executable Code:
!  ------------------------
!
      call getOtherNearCellIndicesGrid2d( proj, posc, ffp, pmax, ncn, icnear )
!          =============================   i     i     i    i     o     o
!
!  Exclude current cell and empty cells
!
      nn = ncn
      ncn = 0
      do in = 1, nn
         ic = icnear(in)
         if( cell(ic)%empty )  cycle            ! exclude empty cells
         if( ic == proj%id )  cycle             ! exclude "present" cell
         ncn = ncn + 1
         icnear(ncn) = ic
      end do

   end subroutine getOtherNearCellIndices2dOld
!
!=======================================================================
!
   subroutine getProbCell2dOld( ic, prob )
!             ================  i    o
!
!     Get the probabilities to encounter each of the atoms in cell ic
!
!  Arguments:
!
      integer(INT4), intent(in) :: ic           ! cell index x
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
!
   end subroutine getProbCell2dOld
!
!=======================================================================
!
   subroutine getDensCell2dOld( ic, dens )
!             ================  i    o
!
!     Get the total atomic density in cell ic
!
!  Arguments:
!
      integer(INT4), intent(in) :: ic           ! cell index x
      real(REAL8), intent(out) :: dens          ! density
!
!  Local Variables:
!
!     None
!
!  Start of Executable Code:
!  ------------------------
!
      dens = cell(ic)%dens_tot

   end subroutine getDensCell2dOld
!
!=======================================================================
!
   subroutine getVibCell2dOld( ic, vib )
!             ===============  i    o
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

   end subroutine getVibCell2dOld
!
!=======================================================================
!
   subroutine getXtalCell2dOld( ic, lxtal )
!             ================  i     o
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

   end subroutine getXtalCell2dOld
!
!==========================================================================================
!
   subroutine getEdCell2dOld( ic, ed )
!             ==============  i   o
!
!     Get the displacement energy in cell ic
!
!  Arguments:
!
      integer(INT4), intent(in) :: ic           ! cell index
      real(REAL8), intent(out) :: ed(:)         ! surface binding energy
!
!  Local Variables:
!
!     None
!
!  Start of Executable Code:
!  ------------------------
!
      ed(:) = cell(ic)%ed(:)

   end subroutine getEdCell2dOld
!
!=======================================================================
!
   subroutine getEsurfCell2dOld( ic, esurf )
!             =================  i     o
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

   end subroutine getEsurfCell2dOld
!
!=======================================================================
!
   subroutine relaxCell2dOld( iion )
!             ==============   i
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
      integer(INT4), intent(in) :: iion         ! index of ion
!
!  Local Variables:
!
      type :: Ovlp_type
!             =========
         sequence
         integer(INT4), allocatable :: ic(:)    ! cell indices of overlapped cells
         real(REAL8), allocatable :: area(:)    ! overlap area
         real(REAL8) :: emptyarea               ! total overlap area with empty cells
      end type Ovlp_type

      integer(INT4) :: ic, iarea, ico, istat
      character(256) :: msg
      real(REAL8) :: dens0(nelem), relaxedarea
      real(REAL8), allocatable :: densrlx(:,:)
      real(REAL8), allocatable :: stress(:,:), ff0(:), dpa(:)
      real(REAL8) :: ed
      type(Ovlp_type), allocatable :: ovlp(:)
!
!  Start of Executable Code:
!  ------------------------
!
      allocate( densrlx( nelem, ncell ), stat=istat, errmsg=msg )
      call check_allocation( 'relaxCell2dOld', 'densrlx', istat, msg )
      allocate( stress( 13, ncell ), stat=istat, errmsg=msg )
      call check_allocation( 'relaxCell2dOld', 'stress', istat, msg )
      allocate( ff0( ncell ), stat=istat, errmsg=msg )
      call check_allocation( 'relaxCell2dOld', 'ff0', istat, msg )
      allocate( dpa( ncell ), stat=istat, errmsg=msg )
      call check_allocation( 'relaxCell2dOld', 'dpa', istat, msg )
      allocate( ovlp( ncell ), stat=istat, errmsg=msg )
      call check_allocation( 'relaxCell2dOld', 'ovlp', istat, msg )
!
!  Control output
!
      dens0(:) = ZERO
      do ic = 1, ncell
         dens0(:) = dens0(:) + cell(ic)%dens(:) * cell(ic)%area
      end do
      print '(a,1p,100e10.2)',  'Material BEFORE relaxation:', dens0(:)
!
!  Write CELL2 and SURF files before relaxation (we need the fill fractions)
!
      do ic = 1, ncell
         ff0(ic) = cell(ic)%ff
         cell(ic)%ff = sum( cell(ic)%dens(:) / elem(:)%dens )
      end do

      print *, 'iion, niongeom, nionout =', iion, niongeom, nionout
      if( nionout > 0 ) then
         if( mod( iion, nionout ) == 0 .or. iion == nion ) then
            if( out%lcell )  call writeCell2dOld( iion )
!                                 ==============    i
            if( out%lsurf )  call writeSurface2d( iion )
!                                 ==============   i
         end if
      end if
!
!  Relaxation
!
      if( dam%relaxmodel == ENUM%relaxmodel%viscoelastic ) then
         do ic = 1, ncell
            stress(:,ic) = cell(ic)%stress(:)
         end do
      end if

      do ic = 1, ncell
         ed = sum( cell(ic)%ed(:) * cell(ic)%prob(:) )
         dpa(ic) = cell(ic)%edam / (2.5 * ed) / sum( cell(ic)%dens(:) )
      end do

      call relaxGrid2d( ff0(:), cell(:)%ff, dpa(:), stress(:,:), ovlp(:) )
!          ===========   i          i         i        io          o

      if( dam%relaxmodel == ENUM%relaxmodel%viscoelastic ) then
         do ic = 1, ncell
            cell(ic)%stress(:) = stress(:,ic)
         end do
      end if
!
!  Calculate relaxed densities
!
!  Note: cell(ic)%area is the area of the unrelaxed old cell
!        sum( ovlp(ic)%area(:) ) is the area of the relaxed old cell
!
      do ic = 1, ncell
         relaxedarea = sum( ovlp(ic)%area(:) ) + ovlp(ic)%emptyarea
         cell(ic)%ff = cell(ic)%ff * cell(ic)%area / relaxedarea
         cell(ic)%dens(:) = cell(ic)%dens(:) * cell(ic)%area / relaxedarea
         if( dam%relaxmodel == ENUM%relaxmodel%viscoelastic )  &
            cell(ic)%stress(:) = cell(ic)%stress(:) * cell(ic)%area / relaxedarea
      end do
!
!  Write CELL2 and SURF file
!
      if( nionout > 0 ) then
         if( mod( iion, nionout ) == 0 .or. iion == nion ) then
            if( out%lcell )  call writeCell2dOld( iion, relaxed_=.true., dpa_=dpa(:) )
!                                 ==============   i           i             i
            if( out%lsurf )  call writeSurface2d( iion )
!                                 ==============   i
         end if
      end if
!
!  Interpolate relaxed cell contents back to original cells
!
      do ic = 1, ncell
         densrlx(:,ic) = cell(ic)%dens(:)
         cell(ic)%dens(:) = SMALL
         if( dam%relaxmodel == ENUM%relaxmodel%viscoelastic ) then
            stress(:,ic) = cell(ic)%stress(:)
            cell(ic)%stress(:) = ZERO
         end if
      end do

      do ic = 1, ncell
         do iarea = 1, size( ovlp(ic)%area(:) )
            ico = ovlp(ic)%ic(iarea)
            cell(ico)%dens(:) = cell(ico)%dens(:) + densrlx(:,ic) * ovlp(ic)%area(iarea)
            if( dam%relaxmodel == ENUM%relaxmodel%viscoelastic )  &
               cell(ico)%stress(:) = cell(ico)%stress(:) + stress(:,ic) * ovlp(ic)%area(iarea)
         end do
      end do
!
!  Calculate areas of new cells and convert cell contents into densities
!
      call getAreasGrid2d( cell(:)%area )
!          ==============        o
!
!      write( 10, '(a)' ) ' '
      do ic = 1, ncell
!         write( 10, '(i5,10(i5,es10.3))' )  &
!            ic, (ovlp(ic)%ic(iarea), ovlp(ic)%area(iarea), iarea=1,size( ovlp(ic)%area ))
!         do iarea = 1, size( ovlp(ic)%area(:) )
!            ico = ovlp(ic)%ic(iarea)
!            cell(ico)%area = cell(ico)%area + ovlp(ic)%area(iarea)
!         end do
         deallocate( ovlp(ic)%area, ovlp(ic)%ic )
      end do
!
!debug
!      write( 10, '(a)' ) ' '
!      write( 10, '((es12.3))' ) cell(:)%area
!end-debug
!
      do ic = 1, ncell
         if( .not. all( cell(ic)%dens(:) == SMALL ) )  &
            cell(ic)%dens(:) = cell(ic)%dens(:) / cell(ic)%area
         if( dam%relaxmodel == ENUM%relaxmodel%viscoelastic )  &
            cell(ic)%stress(:) = cell(ic)%stress(:) / cell(ic)%area
      end do
!
!  Calculate derived quantities
!
      call initDerivedCell2dOld( )
!          ====================
!
!  Control output
!
      dens0(:) = ZERO
      do ic = 1, ncell
         dens0(:) = dens0(:) + cell(ic)%dens(:) * cell(ic)%area
      end do
      print '(a,1p,100e10.2)',  'Material AFTER relaxation:', dens0(:)
!
!  For the moment we do not move lattice vibration amplitude and crystal flag back to
!  original grid. For frequent enough updates there would not be a movement anyway
!
!  Write CELL2 and SURF files after relaxation
!
      if( nionout > 0 ) then
         if( mod( iion, nionout ) == 0 .or. iion == nion ) then
            if( out%lcell )  call writeCell2dOld( iion )
!                                 ==============   i
            if( out%lsurf )  call writeSurface2d( iion )
!                                 ==============   i
         end if
      end if

      return

   end subroutine relaxCell2dOld
!
!=======================================================================
!
   subroutine writeCell2dOld( iion, relaxed_, dpa_ )
!             ==============   i       i       i
!
!     Write the contents of the target cells to the CELL file (same format as HIS file)
!
!  Arguments:
!
      integer(INT4), intent(in) :: iion         ! index of ion
      logical, optional :: relaxed_             ! indicates relaxed grid to be written
      real(REAL8), optional :: dpa_(:)          ! displacements per atom
!
!  Local Variables:
!
      integer(INT4) :: ic
      real(REAL8) :: dose, p

      character(6), save :: fmode = 'write'
!
!  Start of Executable Code:
!  ------------------------
!
!  Open CELL2 file
!
      call openFile( 22, 'CELL2', 'formatted', fmode )
!          ========  i      i          i         i
      fmode = 'append'
!
!  Write dose to first line
!
      dose = iion * ion%ddose
      write( 22, '(a,1p,e15.6)' ) '# Dose:', dose
!
!  Write grid
!
      call writeGrid2d( 22, relaxed_ )
!          ===========  i      i
!
!  Write cell densities
!
      write( 22, '(a)' )  '# Values:'
      if( dam%relaxmodel == ENUM%relaxmodel%none  &
          .or. dam%relaxmodel == ENUM%relaxmodel%full ) then
         write( 22, '(2i15)' )  nelem+1, ncell
         write( 22, '(12i15)' )  0, elem(:)%z
      else if( dam%relaxmodel == ENUM%relaxmodel%viscoelastic ) then
         write( 22, '(2i15)' )  nelem+6, ncell
         write( 22, '(12i15)' )  0, elem(:)%z, 0, 11, 22, 12, 0
      else
         error stop 'writeCell2dOld: Invalid dam%relaxmodel'
      end if

      if( dam%relaxmodel == ENUM%relaxmodel%none  &
          .or. dam%relaxmodel == ENUM%relaxmodel%full ) then
         do ic = 1, ncell
            write( 22, '(1p,12e15.6)' )  cell(ic)%ff, cell(ic)%dens(:)
         end do
      else if( dam%relaxmodel == ENUM%relaxmodel%viscoelastic ) then
         do ic = 1, ncell
            p = cell(ic)%stress(1)
            write( 22, '(1p,12e15.6)', advance='no' )  &
               cell(ic)%ff, 1e24*cell(ic)%dens(:), p,  &
               sum( cell(ic)%stress(2::3) ) / 4 - p,  &
               sum( cell(ic)%stress(3::3) ) / 4 - p,  &
               sum( cell(ic)%stress(4::3) ) / 4
            if( present( dpa_ ) ) then
               write( 22, '(es15.6)' )  dpa_(ic)
            else
               write( 22, '(a)' )  ' 0.'
            end if
         end do
      end if

      close( 22 )

   end subroutine writeCell2dOld
!
!=======================================================================
!
   subroutine readAtomsFromCell2dOld( )
!             ======================
!
!     Read the atomic numbers from record cellrec of the CELL2 file cellfile
!
!     UNTESTED CODE !
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
!  Open file
!
      call openFile( 22, cellfile, 'formatted', 'read' )
!          ========  i      i           i         i

      error stop 'readAtomsFromCell2dOld: Not yet implemented'

   end subroutine readAtomsFromCell2dOld
!
!===================================================================================
!
   subroutine readCell2dOld( )
!             =============
!
!     Read the contents of the target cells from record cellrec of cellfile
!
!     UNTESTED CODE !
!
!  Arguments:
!
!     None
!
!  Local Variables:
!
      logical :: valid
      character(80) :: line, text
      character(256) :: errmsg
      integer(INT4) :: ic, nrec, nrow, ncol, i, ielem, istat, irec
      real(REAL8), allocatable :: z(:)
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

      do irec = 1, nrec
!
!  Ignore first line
!
         read( 22, '(a)', err=800, end=800 )  line
!
!  Deallocate densities
!
         if( allocated( cell ) ) then
            do ic = 1, ncell
               deallocate( cell(ic)%dens )
            end do
            deallocate( cell )
         end if
!
!  Read grid
!
         call readGrid2d( 22, ncell, valid )
!             ==========  i     o      o
!
!  Read cell contents
!
         read( 22, '(a)', err=900, end=900 )  line
         read( 22, *, err=900, end=900 )  ncol, nrow
!         if( dam%relaxmodel == ENUM%relaxmodel%none  &
!             .or. dam%relaxmodel == ENUM%relaxmodel%full ) then
!            if( ncol /= nelem )  goto 902
!         else if( dam%relaxmodel == ENUM%relaxmodel%viscoelastic ) then
!            if( ncol /= nelem+3 )  goto 902
!         else
!            error stop 'readCell2dOld: Invalid dam%relaxmodel'
!         end if                                 ! dam not available at this point
         if( nrow /= ncell )  goto 903
!
!  Read the atomic numbers and determine the atom masses (there are probably bugs in this
!  code paragraph)
!
         if( ldyn ) then
            allocate( z( nelem ) )
            read( 22, *, err=900, end=900 ) i, z(:)
            do ielem = 1, nelemIon
               if( z(ielem) /= elem(ielem)%z )  &
                  error stop 'readCell2dOld: Inconsistent ion atomic number'
            end do
            do ielem = nelemIon+1, nelem
               elem(ielem)%z = z(ielem)
            end do
            deallocate( z )
         else
            read( 22, *, err=900, end=900 ) i, elem(nelemIon+1:nelem)%z
         end if
! old:
!         read( 22, *, err=900, end=900 )  (i, ielem=0,ncol)
!
         allocate( cell( ncell ), stat=istat, errmsg=errmsg )
         call check_allocation( 'readCell2dOld', 'cell', istat, errmsg )

         do ic = 1, ncell
            allocate( cell(ic)%dens( nelem ), stat=istat, errmsg=errmsg )
            write( text, '(a,i0,a)' )  'cell(', ic, ')%dens'
            call check_allocation( 'readCell2dOld', text, istat, errmsg )

            if( ncol == nelem ) then            ! nelem+1 ?
               read( 22, *, err=900, end=900 )  cell(ic)%ff, cell(ic)%dens(:)
               cell(ic)%dens(:) = cell(ic)%dens(:) / 1e24
            else
               read( 22, *, err=900, end=900 )  &
                  cell(ic)%ff, cell(ic)%dens(:), cell(ic)%stress(:)
            end if
         end do

      end do

800   if( irec == 1 )  goto 901
      if( .not. valid )  goto 904

      close( 22 )
!
!  Allocate and initialize further cell properties
!
      do ic = 1, ncell
                                             ! Thermal vibration amplitudes
         allocate( cell(ic)%vib( nelem ), stat=istat, errmsg=errmsg )
         write( text, '(a,i0,a)' )  'cell(', ic, ')%vib'
         call check_allocation( 'readCell2dOld', text, istat, errmsg )

         cell(ic)%vib(:) = ZERO
                                             ! Crystal flag
         cell(ic)%lxtal = .false.            !    amorphous simulation can always be done

         cell(ic)%edam = ZERO

      end do
                                             ! Cell areas
      call getAreasGrid2d( cell(:)%area )
!          ==============       o
!
!  Calculate derived quantities
!
      call initDerivedCell2dOld( )
!          ====================

      return

900   error stop 'readCell2dOld: error during reading CELL2 file'
901   error stop 'readCell2dOld: no data on CELL2 file'
902   error stop 'readCell2dOld: number of elements on CELL2 file inconsistent with input file'
903   error stop 'readCell2dOld: inconsistent dimensions on CELL2 file'
904   error stop 'readCell2dOld: grid type cannot be read'

   end subroutine readCell2dOld
!
!=======================================================================
!
   subroutine initDerivedCell2dOld( )
!             ====================
!
!     Initialize derived cell properties that do not depend on the grid
!
!  Arguments:
!
!     None
!
!  Local Variables:
!
      character(80) :: text
      character(256) :: errmsg
      integer(INT4) :: ielem, istat, ic
!
!  Start of Executable Code:
!  ------------------------
!
!  Allocate unallocated arrays
!
      if( .not. allocated( cell(1)%prob ) ) then
         do ic = 1, ncell
            allocate( cell(ic)%prob( nelem ), stat=istat, errmsg=errmsg )
            write( text, '(a,i0,a)' )  'cell(', ic, ')%prob'
            call check_allocation( 'initDerivedCell2dOld', text, istat, errmsg )

            allocate( cell(ic)%ed( nelem ), stat=istat, errmsg=errmsg )
            write( text, '(a,i0,a)' )  'cell(', ic, ')%ed'
            call check_allocation( 'initDerivedCell2dOld', text, istat, errmsg )

            allocate( cell(ic)%esurf( nelem ), stat=istat, errmsg=errmsg )
            write( text, '(a,i0,a)' )  'cell(', ic, ')%esurf'
            call check_allocation( 'initDerivedCell2dOld', text, istat, errmsg )
         end do
      end if
!
!  Derived quantities
!
      do ic = 1, ncell
                                                ! total density
         cell(ic)%dens_tot = sum( cell(ic)%dens(:) )
                                                ! probabilities
         if( cell(ic)%dens_tot /= ZERO ) then
            cell(ic)%prob(:) = cell(ic)%dens(:) / cell(ic)%dens_tot
         else
            cell(ic)%prob(:) = ZERO
         end if
                                                ! fill fraction
         cell(ic)%ff = sum( cell(ic)%dens(:) / elem(:)%dens )
                                                ! empty flag
         cell(ic)%empty = all( cell(ic)%dens(:) == SMALL )
                                                ! displacement energies
         do ielem = 1, nelem
            cell(ic)%ed(ielem) = sum( cell(ic)%prob(:) * dam%ed(:,ielem) )
         end do
                                                ! surface binding energies
         do ielem = 1, nelem
            cell(ic)%esurf(ielem) = sum( cell(ic)%prob(:) * dam%esurf(:,ielem) )
         end do
                                                ! cumulative damage energy
         cell(ic)%edam = ZERO
      end do
!
!  Call grid dependent derived initializations
!
      call initPostCellGrid2d( cell(:)%ff, cell(:)%empty )
!          ==================      i             i

   end subroutine initDerivedCell2dOld
!
!==========================================================================================
!
end module Target_Cell2d_Old_mod
