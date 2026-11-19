module Eval_Det_mod
!
!  Provide an atom list for given lattice sites. Defects are
!
!     -) deterministic, i.e. their types and positions do not change during the simulation
!     -) of certain types defined by the positions of their atoms relative to the lattice
!        sites
!
!  NOTE: It would be more elegant to introduce a Column_type and a module providing 
!        operations on columns (analogous to histograms in the statistic mode) 
!
!  Used Modules:
!  ============
!
   use Const_mod,             only: INT4, REAL8,  &
                                    ZERO, ONE, PI, HALF, TWO, THIRD
   use Types_mod,             only: Proj_type, Ipos_type, Site_type, Atom_type,  &
                                    LSDefect_type,  NATOMM

   use Open_mod,              only: openFile
   use Random_mod,            only: rannum, ranint

   use Param_Setup_mod,       only: nr, ldamdyn, usehis, usebk, deffile, nelemIon, nelem,  &
                                    atom2, usecell
   use Param_Elem_mod,        only: elem
   use Param_Dam_mod,         only: dam
   use Param_Mat_mod,         only: mat, coord_xtal, irxtal
   use Param_Out_mod,         only: out

   use Geom_CoordSys_mod,     only: toCoordSys, fromCoordSys
   use Stat_His1d_mod,        only: His1d_type, scoreHis1d, getHis1d
   use Target_mod,            only: getDomainDensities
   use Xtal_Defects_Diam_mod, only: getDefnAtoms, getDefnNeigh, getDefType, initDef,  &
                                    getDefStrain, getDefLatPos, getDefnLS, getProbDef
   use Xtal_Lattice_mod,      only: getNearestSite, maxSizeLattice, getSites
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
   public :: det10                  ! initialize deterministic defects structure
!            =====                   
   public :: his10                  ! initialize histograms when usedef=T and statistic 
!            =====                  !    model of damage is to be used
   public :: det23                  ! score a vacancy/recoil into lsdef structure
!            =====
   public :: det30                  ! write defects into file
!            =====
   public :: his30                  ! score damage from lsdef to his1
!            =====
   public :: det60                  ! determine LSDefect at lattice site
!            =====
   public :: getCol                 ! get the boundaries of the column
!            ======
!
!  Private Procedures:
!  ==================
!
   private :: createLattice                  ! create xtal lattice (lsdef array)
!             =============                 
   private :: createGrid                     ! create random grid where the defects 
!             ==========                     !    will be placed
   private :: assocLSDef                     ! associate lsdef with lsdef_ptr
!             ==========  
   private :: strainDef                      ! apply strain to one defect
!             =========
   private :: strain                         ! apply strain to lattice sites 
!             ======
   private :: dethis10                       ! score defects read from the file into 
!             ========                       !    1-D damage histograms
   private :: balanceDef                     ! equalize number of vacancy and interstitials
!             ==========                     !
   private :: getIndices                     ! calculate the indices of the lattice sites
!             ==========                     !    from the lattice site position
   private :: writeLSDRas                    ! write deterministic defects in rasmol format
!             ===========
   private :: writeDefects                   ! write defects into the file 
!             ============                   
   private :: readDefects                    ! read defects from the file 
!             ===========                        
!
!  Private Data Types:
!  ==================
!
   type :: LSDefect_ptr_type                 ! Lattice site defect pointer type
!          =================
      type(LSDefect_type), pointer :: ptr => null( )
                                             ! pointer to lattice site defect
   end type LSDefect_ptr_type
!
!  Private Data:
!  ============
!
   integer(INT4) :: icolmin(3)               ! lower boundaries of column 
   integer(INT4) :: icolmax(3)               ! upper boundaries of column 

   integer(INT4) :: nlsdef                   ! number of lattice site defects
   integer(INT4) :: nlsstr                   ! number of strained lattice sites

   type(LSDefect_type), allocatable, target :: lsdef(:) ! lattice site defects
   type(LSDefect_ptr_type), allocatable :: lsdef_ptr(:,:,:,:,:,:)
                                             ! lattice site defect pointer

contains
!
!===========================================================================================
!
   subroutine det10( his1 )
!             =====   io
!
!     Initialize the structure necessary for deterministic damage
!
!     This subroutine must be called by the master thread only
!
!  Arguments
!
      type(His1d_type), intent(inout) :: his1      ! 1d histogram
!
!  Local Variables:
!
      integer(INT4) :: ilsdef
      logical, parameter :: writeRAS = .false.  ! to write RAS
!
!  Start of Executable Code:
!  ------------------------
!
!  Make sure the wafer surface normal is not perpendicular to the crystal z axis
!
      if( usecell )  error stop 'det10: usecell=t not allowed with deterministic defects'
      if( any( mat(1:nr)%wafer(3) == ZERO ) )  &
            error stop 'det10: wafer(3)=0 not allowed with determinisitic defects'

      call initDef( )
!          =======

      if( usehis .or. usebk ) then

         icolmin(:) = [ 0, 0, 0 ]             ! set the column bounds
         icolmax(:) = icolmin(:) + dam%wcol(:)

         call createLattice( )                  ! lattice in inter. coor. sys.
!             =============   
         if( out%lhis ) then 
            call createGrid( his1 )             ! Create the pseudorandom grid 
!                ==========    i      
            his1%hist(:,:) = ZERO               ! reset histograms to zero
         end if
      else
                                                   ! read defects from the file 
         call readDefects( lsdef, nlsdef, nlsstr, icolmin, icolmax )
!             ===========   io     o         o       o       o
         call assocLSDef( )                        ! lsdef_ptr => lsdef
!             ==========  
      end if
                                                  ! rasmol file for checking
      if( writeRAS )  call writeLSDRas( lsdef(:), count( lsdef(:)%type /= 0 ), 1 )
!                          ===========    i                 i                i
!
!  Apply strain to defects neigboring lattice sites 
!
      do ilsdef = 1, size( lsdef )
         if( lsdef(ilsdef)%type <= 0 ) cycle      ! avoid unnecessary computation
         call strainDef( lsdef(ilsdef) )
!             =========       io
      end do

   end subroutine det10
!
!===========================================================================================
!
   subroutine his10( his1 )
!             =====   io
!
!    Initialize the damage histograms according to damage distribution read
!    from the deterministic defect file. 
!
!  Arguments
!
      type(His1d_type), intent(inout) :: his1   ! 1d histogram
!
!  Local Variables:
!
      integer(INT4) :: status
!
!  Start of Executable Code:
!  ------------------------
!
!  Initialize the strain
!
      call initDef( )
!          =======
                                            ! read defects from the file
      call readDefects( lsdef, nlsdef, nlsstr, icolmin, icolmax )
!          ===========   io     o        o       o        o
!
!  In case usedef=t score distribution of deterministic defects in histogram. 
!
      call dethis10( his1 )
!          ========   io
!
!  Deallocate lsdef since deterministic model will not be used
!
      deallocate( lsdef, stat = status )
      if( status /= 0 )  error stop 'his10: lsdef not deallocated'

   end subroutine his10
!
!===========================================================================================
!
   subroutine det23( proj, aincr )
!             =====   i      i 
!
!     Score vacancy or interstitial into the deterministic defect structure.
!
!     Remark: Coordinates of lattice sites are in crystal coordinate system
!
!  Arguments:
!
      type(Proj_type), intent(in) :: proj
      integer(INT4), intent(in) :: aincr        ! atom increment (-1 vacancy, +1 recoil)
!
!  Local Variables:
!
      logical :: isDefect
      integer(INT4) :: ir, naold, nanew, typold, typnew, orold, ornew
      integer(INT4) :: i, j, k, l, m, n, ipos(3)
      type(Ipos_type) :: jpos
      real(REAL8) :: posi(3), rand, dpos(3)
!
!  Start of Executable Code:
!  ------------------------
!
!  Check whether to store this projectile
!
      if( .not.dam%ldam .and. proj%origin%internal )  return

      ir = proj%id

      if( ir /= irxtal )  return                 ! error stop 'det23: not crystalline region'
!
!  Find projectile nearest lattice site position 
!
      call getNearestSite( proj%pos(:), ir, jpos )
!          ==============      i        i    o
!
!  Convert input to crystal coordinate system
!
      call toCoordSys( coord_xtal, proj%pos(:), posi(:) )
!          ==========        i             i          o
!
      dpos(:) = posi(:) - jpos%ipos(:)
!
!  Map position into the column parallel to the surface
!
      ipos(1:2) = icolmin(1:2) + modulo( jpos%ipos(1:2) - icolmin(1:2), dam%wcol(1:2) )
      ipos(3) = jpos%ipos(3) + sum( nint( (jpos%ipos(1:2) - ipos(1:2)) *  &
                                          mat(ir)%wafer(1:2) / mat(ir)%wafer(3) ) )
!
      call getIndices( ipos(:), i, j, k, l, m, n )
!          ==========    i      o  o  o  o  o  o
!
!  Determine whether a defect is already at this lattice site. It has to be
!
      isDefect = all( [n,m,l,k,j,i] >= lbound( lsdef_ptr ) ) .and.  &
                 all( [n,m,l,k,j,i] <= ubound( lsdef_ptr ) )
      if( isDefect )  isDefect = associated( lsdef_ptr(n,m,l,k,j,i)%ptr )

      if( isDefect ) then
!
!  Defect. Add vacancy to the structure. Remove one atom
!                           
         if( any( lsdef_ptr(n,m,l,k,j,i)%ptr%ipos(:) /= ipos(:) ) ) then 
            print *,  lsdef_ptr(n,m,l,k,j,i)%ptr%ipos(:), ipos(:)
            error stop 'det23: inconsistent lattice site position'
         end if
!
!  See what we already have on this site. Check natoms. 
!
         naold = lsdef_ptr(n,m,l,k,j,i)%ptr%nsa
         orold = lsdef_ptr(n,m,l,k,j,i)%ptr%ior
         typold = lsdef_ptr(n,m,l,k,j,i)%ptr%type

         nanew = naold + aincr
!
!  Separate treatment for arbitrary defects
!
         if( typold == -1 ) then
            lsdef_ptr(n,m,l,k,j,i)%ptr%nsa = nanew
            if( aincr > 0 )  lsdef_ptr(n,m,l,k,j,i)%ptr%dpos(:,nanew) = dpos(:)
            return
         end if
         if( naold /= getDefnAtoms( typold ) )  return !error stop' det23: inconsistent typold '
!                     ============
!
!  Do not score if too many/low atoms at lattice site
!
         if( aincr < 0 .and. naold == 0 ) return
         if( aincr > 0 .and. naold == NATOMM ) return

         call getDefType( typnew, nanew ) ! get new defect type
!             ==========    o       i

         if( nanew /= getDefnAtoms( typnew ) )  error stop 'det23: inconsistent typnew'
!                     ============    i

!  Determinine properties of new defect type at this lattice site
!
         ornew = orold
         if( typold <= 1 .and. typnew > 1 ) then
            call rannum( rand )
!                ======    o
            ornew = int( 24 * rand ) + 1
         end if
!
!  But we have to remove strain that this defect caused 
!
         call strainDef( lsdef_ptr(n,m,l,k,j,i)%ptr, sub_=.true. )
!             =========               io                 i

         lsdef_ptr(n,m,l,k,j,i)%ptr%ior = ornew
         lsdef_ptr(n,m,l,k,j,i)%ptr%nsa = nanew
         lsdef_ptr(n,m,l,k,j,i)%ptr%type = typnew
         if( aincr > 0 )  lsdef_ptr(n,m,l,k,j,i)%ptr%ielem(nanew) = proj%ielem
!
!  Add the strain of new defect 
!
         call strainDef( lsdef_ptr(n,m,l,k,j,i)%ptr, sub_=.false. )
!             =========              io                 i             
!
         if( typnew > 0 .and. typold < 1 )  nlsdef = nlsdef + 1
         if( typnew < 1 .and. typold > 0 )  nlsdef = nlsdef - 1
      end if

   end subroutine det23
!
!============================================================================================
!
   subroutine det30( )
!             =====  
!
!     Calls writeDefects. Just interface to eval_stat_mod.f90 
!
!  Calling arguments:
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
!  If deterministic model is used just write the lsdef into file
!
      call writeDefects( lsdef(:), count( lsdef(:)%type /= 0 ), icolmin(:), icolmax(:) )
!          ============    i                    i                  i           i

   end subroutine det30
!
!============================================================================================
!
   subroutine his30( his1 )
!             =====   io
!
!    Initialize the damage histograms according to damage in simulation cell. 
!
!  Calling arguments
!
      type(His1d_type), intent(inout) :: his1   ! 1d histogram
!
!  Local Variables:
!
!     None
!
!  Start of Executable Code:
!  ------------------------
!
!  In case usedef=t score distribution of deterministic defects in histogram. 
!
      call dethis10( his1 )
!          ========   io

   end subroutine his30
!
!===========================================================================================
!
   subroutine det60( jpos, ir, lsdef_rtn, isDefect )
!             =====    i   i       o         o
!
!     Determine whether there is a lattice site defect at the given lattice position,
!     and return its value
!
!  Arguments:
!
      integer(INT4), intent(in) :: jpos(3)            ! lattice site position
      integer(INT4), intent(in) :: ir                 ! index of region
      type(LSDefect_type), intent(out) :: lsdef_rtn   ! lattice site defect (only set if 
                                                      !    isDefect=.true.)
      logical, intent(out) :: isDefect                ! flag indicating whether there is a defect
!
!  Local variables:
!
      integer(INT4) :: ipos(3), i, j, k, l, m, n
!
!  Start of Executable Code:
!  ------------------------
!
!  Map position into the column parallel to the surface
!
      ipos(1:2) = icolmin(1:2) + modulo( jpos(1:2) - icolmin(1:2), dam%wcol(1:2) )
      ipos(3) = jpos(3) + sum( nint( (jpos(1:2) - ipos(1:2)) *  &
                                     mat(ir)%wafer(1:2) / mat(ir)%wafer(3) ) )

      call getIndices( ipos(:), i, j, k, l, m, n )
!          ==========    i      o  o  o  o  o  o
!
!  Determine whether a defect is at this lattice site
!
      isDefect = .false.

      if( nlsdef == 0 ) return

      if( any( [n,m,l,k,j,i] < lbound( lsdef_ptr ) ) .or.  &
          any( [n,m,l,k,j,i] > ubound( lsdef_ptr ) ) )  return
          
      if( .not. associated( lsdef_ptr(n,m,l,k,j,i)%ptr ) ) return

                                                ! iff dpos(:,:) == 0 consider site as a
      if( all( lsdef_ptr(n,m,l,k,j,i)%ptr%dpos(:,:) == ZERO ) )  return  ! lattice site

      if( .not. ldamdyn )  return               ! considered only if damage to be considered

      isDefect = .true.                         ! all checks passed ... is a defect
      lsdef_rtn = lsdef_ptr(n,m,l,k,j,i)%ptr
!
!  Consistency check
!
      if( any( lsdef_ptr(n,m,l,k,j,i)%ptr%ipos(:) /= ipos(:) ) )  &
            error stop 'det60: inconsistent lattice site position'

   end subroutine det60      
!
!===========================================================================================
!
   subroutine getCol( imin, imax )
!             ======   o     o
!
!     Get the boundaries of the defect column 
!
!  Arguments:
!
      integer(INT4), intent(out) :: imin(1:2)   ! lower boundaries in x and y direction
      integer(INT4), intent(out) :: imax(1:2)   ! upper boundaries in x and y direction
!
!  Local Variables:
!
!     None
!
!  Start of Executable Code:
!  ------------------------
!
      imin(1:2) = icolmin(1:2)
      imax(1:2) = icolmax(1:2)

   end subroutine getCol
!
!===========================================================================================
!
   subroutine createLattice( )
!             =============      
!
!     Create the array of Silicon lattice coordinates with given bounds in inter. coor. sys
!
!  Arguments:
!
!     None
!
!  Local Variables:
!
      integer(INT4) :: cmin(3),  cmax(3)        ! minmax unit cell index 
      integer(INT4) :: ix, iy, iz, i, j, k, l, m, n, ielemr, nelemr, ielem(nelem)
      integer(INT4) :: status, isite, nsite, ilatt
      type(Ipos_type), allocatable, save :: iposl(:)
      real(REAL8) :: c1(3,2), posi(3), prob(nelem)
!
!  Start of Executable Code:
!  ------------------------
!
      if( any( icolmax <= icolmin ) )  error stop 'createLattice: ERROR coulumn size!'
      ! check on full number of unit cells
      if( all( modulo( icolmax(:), 4 ) /= 0 ) )  error stop 'createLattice: error icolmax'
      if( all( modulo( icolmin(:), 4 ) /= 0 ) )  error stop 'createLattice: error icolmin'

      if( .not. allocated( iposl ) )  allocate( iposl( maxSizeLattice() ) )
!                                                      ==============
!
!  Calculate bounds of 6D pointer and number of lattice sites
!
      cmin(:) = icolmin(:) / 4
      cmax(:) = ( icolmax(:) - 1 ) / 4

      allocate( lsdef_ptr( 0:1, 0:1, 0:1, cmin(3):cmax(3), cmin(2):cmax(2), cmin(1):cmax(1) ), &
                stat = status )
      if( status /= 0 )  error stop 'createLattice: latt_ptr not allocated'

      allocate( lsdef( 1:size( lsdef_ptr ) ), stat = status  )
      if( status /= 0 )  error stop 'createLattice: lsdef not allocated'

      ilatt = 0

      do iz = icolmin(3), icolmax(3), 3
         do iy = icolmin(2), icolmax(2), 3
            do ix = icolmin(1), icolmax(1), 3

               posi(:) = [ real(ix), real(iy), real(iz) ]
!
!  Get new lattice site list around posi(:)
!
               call getSites( posi(:), nsite, iposl(:), c1(:,:) )
!                   ========     i       o       o         o

               site: do isite = 1, nsite
                  if( any( iposl(isite)%ipos(:) < icolmin(:) ) )  cycle site
                  if( any( iposl(isite)%ipos(:) >= icolmax(:) ) )  cycle site

                  call getIndices( iposl(isite)%ipos(:), i, j, k, l, m, n )
!                      ==========          i             o  o  o  o  o  o

                  if( .not. associated( lsdef_ptr(n,m,l,k,j,i)%ptr ) ) then
                     ilatt = ilatt + 1
                     lsdef_ptr(n,m,l,k,j,i)%ptr => lsdef(ilatt)
                     lsdef(ilatt)%ipos(:) = iposl(isite)%ipos(:) 
                     lsdef(ilatt)%dpos(:,:) = ZERO ! initial values

                     call getDomainDensities( irxtal, nelemr, ielem(:), prob_=prob(:) )
!                         ==================    i        o       o           o
                     call ranint( prob(1:nelemr), ielemr )
!                         ======       i            o
                     lsdef(ilatt)%ielem(:) = ielem(ielemr)
                  end if
               end do site
            end do ! loop over x
         end do ! loop over y
      end do ! loop over z

      if( ilatt /= size( lsdef ) )  error stop 'createLattice: ERROR size(lsdef)!'

   end subroutine createLattice
!
!===========================================================================================
!
   subroutine createGrid( his1 )
!             ==========    i 
!    
!     Create the grid for population with defects. Everything here is done in units of 
!     quarter of lattice constants 
!
!  Arguments:
!
      type(His1d_type), intent(in) :: his1    ! damage histogram
!
!  Local Variables:
!  
      integer(INT4), allocatable :: nx(:), ny(:), nz(:)
      real(REAL8), allocatable :: d(:)
      real(REAL8) :: hval, posi(3), pos(3), r1, r2, r3, probDefS, probInt, r11, r22, r33, wboxi
      integer(INT4) :: imin(3), imax(3), defType, ncz
      integer(INT4) :: il, ix, iy, iz, i, j, k, l, m, n, is
      integer(INT4) :: dipos(3), ref(3), iint, ivac, ior, iodd, nl
      type(Ipos_type) :: ipos
!
!  Start of Executable Code:
!  ------------------------
!
!  Get probability of defective site
!
      call getDefType( defType )
!          ==========     o
      call getProbDef( defType, probDefS, probInt )
!          ==========     i         o        o

      imin(:) = icolmin(:) - 4                  ! virtually extend the bounds to have better
      imax(:) = icolmax(:) + 4                  ! defect statistics when creating the grid

      wboxi = his1%wbox / mat(irxtal)%xscale
      ncz = int( wboxi / 4. )       ! number of unit cells in one wbox
      nl = int( ( imax(3) - imin(3) ) / wboxi ) + 1 
      allocate( nx( 1:nl ), ny( 1:nl ), nz( 1:nl ), d( 1:nl ) )

      do il = 1, nl
         pos(3) = il * wboxi * mat(irxtal)%xscale - ONE ! close to end of his1 box point

         call getHis1d( atom2, pos(3), his1, hval )
!             ========    i      i      i     o

         d(il) = TWO * ( mat(irxtal)%dens / ( probDefS * hval ) ) ** THIRD
         if( d(il) < ONE )  error stop 'createGrid: distance between defects too small'
         nx(il) = int( ( imax(1) - imin(1) ) / d(il) ) + 1
         ny(il) = int( ( imax(2) - imin(2) ) / d(il) ) + 1
         if( il == nl ) then
            nz(il) = int( ( ( imax(3) - imin(3) ) - ( il - 1 ) * ncz * 4 ) / d(il) ) + 1
         else
            nz( il ) = int( ( ncz * 4 ) / d( il ) ) + 1
         end if
      end do

      iint = 0 
      ivac = 0 
      nlsdef = 0

      do il = 1, nl
         do iz = 1, nz(il)
            do iy = 1, ny(il)
               do ix = 1, nx(il)             

                  call rannum( r1 )
!                      ======   o
                  call rannum( r2 )
!                      ======   o
                  call rannum( r3 )
!                      ======   o
                  call rannum( r11 )
!                      ======   o
                  call rannum( r22 )
!                      ======   o
                  call rannum( r33 )
!                      ======   o

                  posi(:) = [ ix-1, iy-1, iz-1 ] * d(il)
                  posi(3) = posi(3) + ( il - 1 ) * ncz * 4 ! wboxi
                  posi(1) = posi(1) + HALF * r1 * d(il) * ( ONE - TWO * nint(r11) )
                  posi(2) = posi(2) + HALF * r2 * d(il) * ( ONE - TWO * nint(r22) )
                  posi(3) = posi(3) + HALF * r3 * d(il) * ( ONE - TWO * nint(r33) )
                  posi(:) = posi(:) + imin(:)

                  call fromCoordSys( coord_xtal, posi(:), pos(:) )
!                      ============          i            i        o
                  call getNearestSite( pos(:), irxtal, ipos )
!                      ==============    i       i      o

                  if( any( ipos%ipos(:) >= icolmax(:) ) )  cycle
                  if( any( ipos%ipos(:) < icolmin(:) ) )  cycle

                  ref(:) = ipos%ipos(:)
                  iodd = modulo( ipos%ipos(1), 2 )
                  ior = int( 24 * r1 ) + 1

                  if( r1 < ProbInt ) then                 ! interstitial
                     do is = 1, getDefnLS( defType )
!                               =========     i
                        dipos(:) = getDefLatPos( iodd, 1, ior, defType+is-1 )
!                                  ============    i   i   i        i
                        ipos%ipos(:) = ref(:) + dipos(:)

                        call getIndices( ipos%ipos(:), i, j, k, l, m, n )
!                            ==========        i       o  o  o  o  o  o

                        if( any( ipos%ipos(:) >= icolmax(:) ) )  cycle  ! if outside box 
                        if( any( ipos%ipos(:) < icolmin(:) ) )  cycle   ! cycle

                        if( lsdef_ptr(n,m,l,k,j,i)%ptr%type == 0 ) then
                           lsdef_ptr(n,m,l,k,j,i)%ptr%type =  defType + is - 1  
                           iint = iint + getDefnAtoms( defType+is-1 ) - 1
!                                        ============        i
                           nlsdef = nlsdef + 1
                           lsdef_ptr(n,m,l,k,j,i)%ptr%ior = ior
                        end if
                     end do
                  else                                           ! vacancy
                     call getIndices( ref(:), i, j, k, l, m, n )
!                         ==========    i     o  o  o  o  o  o

                     if( lsdef_ptr(n,m,l,k,j,i)%ptr%type == 0 ) then
                        lsdef_ptr(n,m,l,k,j,i)%ptr%type = 1
                        ivac = ivac + 1
                        nlsdef = nlsdef + 1
                     end if
                  end if
               end do
            end do
         end do
      end do
!
!  Adjust nI and nV. No need for something like this
!
      iint = min( ivac, iint )
      iint = iint - mod( iint, ( getDefnAtoms( defType ) - 1 ) ) 
!                                ============     i
      ivac = iint
!
!  Equalize number of vacancies and interstitials if necessary
!
      call balanceDef( iint, ivac ) !   #int = #vac 
!          ==========   i     i  

   end subroutine createGrid
!
!===============================================================================================
!
   subroutine assocLSDef( )
!             ==========  
!
!     Associate lsdef with lsdef_ptr
!
!  Arguments
!
!     None
!
!  Local Variables:
!
      integer(INT4) ::  status, ilsdef, i, j, k, l, m, n
!
!  Start of Executable Code:
!  ------------------------
!
!  Allocate and associate pointers
!
      allocate( lsdef_ptr( 0:1, 0:1, 0:1, icolmin(3):icolmax(3), icolmin(2):icolmax(2),  &
                           icolmin(1):icolmax(1) ), stat = status )
      if( status /= 0 )  error stop 'assocLSDef: lsdef_ptr not allocated'

      do ilsdef = 1, size( lsdef )

         call getIndices( lsdef(ilsdef)%ipos(:), i, j, k, l, m, n )
!             ==========           i             o  o  o  o  o  o

         lsdef_ptr(n,m,l,k,j,i)%ptr => lsdef(ilsdef)
      end do

   end subroutine assocLSDef
!
!===========================================================================================
!
   subroutine strainDef( defls, sub_ )
!             =========   io     i
!
!     Apply or subtract strain to defects neighboring lattice sites.
!
!  Arguments:
!
      type(LSDefect_type), intent(inout) :: defls     ! defect lattice site
      logical, intent(in), optional :: sub_            ! subtract strain
!
!  Local Variables:
!
      integer(INT4) :: iodd, iposl(3), isa
      real(REAL8) :: dposl(3), cc
!
!  Start of Executable Code:
!  ------------------------
!
      cc = ONE

      if( present( sub_ ) ) then
         if( sub_ )  cc = -ONE
      end if

      iodd = modulo( defls%ipos(1), 2 ) ! Even or odd site
      do isa = 1, defls%nsa
         defls%dpos(:,isa) = defls%dpos(:,isa) +  &
                             cc * getDefStrain( iodd, isa, defls%ior, defls%type )
!                                 ============    i    i      i          i
      end do

      do isa = defls%nsa + 1, getDefnNeigh( defls%type )
!                             ============      i

         iposl(:) = getDefLatPos( iodd, isa, defls%ior, defls%type )
!                   ============    i    i       i          i
         iposl(:) = iposl(:) + defls%ipos(:)

         dposl(:) = cc * getDefStrain( iodd, isa, defls%ior, defls%type )
!                        ============    i    i       i          i
         call strain( iposl(:), dposl(:) )
!             ======     i         i
      end do

   end subroutine strainDef
!
!============================================================================================
!
   subroutine strain( ipos, dpos )
!             ======   i     i
!
!     Apply strain (dpos(1:3)) to lattice site placed at ipos(1:3) 
!
!  Arguments:
!
      integer(INT4), intent(in) :: ipos(1:3)    ! relative positions of neighboring atoms
      real(REAL8), intent(in) :: dpos(1:3)      ! displacement vectors
!                                              
!  Local Variables:
!
      integer(INT4) :: isa, ielemr, nelemr, idir
      integer(INT4) :: i, j, k, l, m, n
      integer(INT4) :: ielem(nelem), iposl(3)
      real(REAL8) :: prob(nelem)
!
!  Start of Executable Code:
!  ------------------------
!
      if( irxtal == 0 )  error stop 'strain: no crystalline region'

      iposl(1:3) = ipos(1:3)

      do idir = 1, 2                      ! map position into column
         if( iposl(idir) < icolmin(idir) ) then
            iposl(idir) = iposl(idir) + dam%wcol(idir)
            iposl(3) = iposl(3) - nint( dam%wcol(idir) * mat(irxtal)%wafer(idir) /  &
                                        mat(irxtal)%wafer(3) )
         else if( ipos(idir) >= icolmax(idir) ) then
            iposl(idir) = ipos(idir) - dam%wcol(idir)
            iposl(3) = iposl(3) + nint( dam%wcol(idir) * mat(irxtal)%wafer(idir) /  &
                                        mat(irxtal)%wafer(3) )
         end if
      end do
!
!  determine whether there is a defect already
!
      call getIndices( iposl(:), i, j, k, l, m, n )
!          ==========    i       o  o  o  o  o  o
!
!  neglect positions outside pointer array boundaries
!
      if( any( [n,m,l,k,j,i] < lbound( lsdef_ptr ) ) .or.  &
          any( [n,m,l,k,j,i] > ubound( lsdef_ptr ) ) )  return

      if( associated( lsdef_ptr(n,m,l,k,j,i)%ptr ) ) then
!
!  there is a defect: add displacements
!
         if( all( lsdef_ptr(n,m,l,k,j,i)%ptr%dpos(:,:) == ZERO ) )  nlsstr = nlsstr + 1
         do isa = 1, lsdef_ptr(n,m,l,k,j,i)%ptr%nsa
            lsdef_ptr(n,m,l,k,j,i)%ptr%dpos(:,isa) =  &
               lsdef_ptr(n,m,l,k,j,i)%ptr%dpos(:,isa) + dpos(:) 
         end do
      else
!
!  no defect: create one
!
         nlsstr = nlsstr + 1
         if( nlsstr > size( lsdef ) )  error stop 'strain: lsdef array too small'
         lsdef(nlsstr)%ipos(1:3) = iposl(1:3)
         lsdef(nlsstr)%nsa = 1

         call getDomainDensities( irxtal, nelemr, ielem(:), prob_=prob(:) )
!             ==================    i        o       o           o
         call ranint( prob(1:nelemr), ielemr )
!             ======        i           o
         lsdef(nlsstr)%ielem(1) = ielem(ielemr)

         lsdef(nlsstr)%dpos(1:3,1) = dpos(1:3)
         lsdef_ptr(n,m,l,k,j,i)%ptr => lsdef(nlsstr)
      end if

   end subroutine strain
!
!==========================================================================================
!
   subroutine dethis10( his1 )
!             ========   io
!
!     Initialize damage histograms according to deterministic defect distribution 
!     Used when the defects are read from deterministic defect file
!
!  Remark: Coordinates of lattice sites are in crystal coordinate system
!
!  Arguments:
!
      type(His1d_type), intent(inout) :: his1 ! damage histogram
!
!  Local Variables:
!
      integer(INT4) :: ilsdef, ielem, isa
      real(REAL8) :: hincr, posi(3)
!
!  Start of Executable Code:
!  ------------------------
!
!  Are the vacancies included
!
      hincr = 1e24 / ( dam%wcol(1) * dam%wcol(2) * mat(irxtal)%xscale**2 )

      do ilsdef = 1, size( lsdef )
         if( lsdef(ilsdef)%type == 0 )  cycle   ! lattice atoms do not contribute to histogram
!
! Rotate to target coordinate system. Only depth (posi(3)) is important
!
         posi(:) = lsdef(ilsdef)%ipos(:) * mat(irxtal)%xscale
         if( lsdef(ilsdef)%nsa == 0 ) then      ! vacancy
            ielem = lsdef(ilsdef)%ielem(1)      ! this could be a problem
            if( dam%lrcoil ) then
               call scoreHis1d( 1+nelem, posi(3), his1, hincr )
!                   ==========     i        i      io     i
            end if
         else
            do isa = 1, lsdef(ilsdef)%nsa - 1   ! one atom has no extra contribution
               ielem = lsdef(ilsdef)%ielem(isa)
                 
               if( ielem <= nelemIon ) then     ! score ions anyway
                  call scoreHis1d( ielem, posi(3), his1, hincr )
!                      ==========    i       i      i      i
               else                             ! target atoms only if dam%ldam or dam%lrcoil
                  if( dam%ldam .or. dam%lrcoil ) then
                     call scoreHis1d( ielem, posi(3), his1, hincr )
!                         ==========    i       i      i      i
                  end if
               end if
            end do
         end if
      end do
!
   end subroutine dethis10
!
!============================================================================================
!
   subroutine balanceDef( nI, nV )
!             ==========  i   i
!
!     Equalize nI and nV and initialize nAtoms for defect Types. 
!
!  Arguments:
!
      integer(INT4), intent(in) :: nI        ! number of interstitials
      integer(INT4), intent(in) :: nV        ! number of vacancies
!
!  Local Variables:
!
      integer(INT4) :: i, iint, ivac
!
!  Start of Executable Code:
!  ------------------------
!
      iint = 0 
      ivac = 0 
      nlsdef = 0

      do i = 1, size( lsdef )
         if( lsdef(i)%type == 1 ) then          ! vacancy
            if( ivac < nV ) then                ! confirm vacancy
               ivac = ivac + 1
               nlsdef = nlsdef + 1
               lsdef(i)%nsa = 0
            else                                ! delete vacancy and 
               lsdef(i)%type = 0                ! set site as lattice atom
               lsdef(i)%nsa = 1
               lsdef(i)%ior = 0                 ! not used
            end if
         elseif( lsdef(i)%type > 1 ) then       ! defect
            if( iint < nI ) then                ! confirm interstitial
               lsdef(i)%nsa = getDefnAtoms( lsdef(i)%type )
!                             ============         i
               lsdef(i)%ielem(1:lsdef(i)%nsa) = atom2 ! set atom type
               iint = iint + lsdef(i)%nsa - 1
               nlsdef = nlsdef + 1
            else                                ! delete interstitial
               lsdef(i)%type = 0                ! set lattice atom
               lsdef(i)%nsa = 1
               lsdef(i)%ior = 0                 ! not used
            end if
         end if
      end do

      nlsstr = nlsdef             ! strained atoms = defect atoms 
                                  ! check
      print *, 'nI=', iint, 'nV=', ivac 
      if( iint /= ivac )  error stop 'balanceDef: iint /=ivac'

   end subroutine balanceDef
!
!==============================================================================================
!
   subroutine getIndices( ipos, i, j, k, l, m, n )
!             ==========   i    o  o  o  o  o  o
!
!     Calculate the indices of the lattice site from its position 
!     NOTE: very similar code in convertIpos (xtal_lattice_mod)
!
!  Arguments:
!
      integer(INT4), intent(in) :: ipos(3)      ! position
      integer(INT4), intent(out) :: i           !
      integer(INT4), intent(out) :: j           !
      integer(INT4), intent(out) :: k           ! indices
      integer(INT4), intent(out) :: l           !
      integer(INT4), intent(out) :: m           !
      integer(INT4), intent(out) :: n           !
!
!  Start of Executable Code:
!  ------------------------
!
      i = ( ipos(1) - modulo( ipos(1), 4 ) ) / 4
      j = ( ipos(2) - modulo( ipos(2), 4 ) ) / 4
      k = ( ipos(3) - modulo( ipos(3), 4 ) ) / 4
      l = modulo( ( ipos(1) - modulo( ipos(1), 2 ) ) / 2, 2 )
      m = modulo( ( ipos(2) - modulo( ipos(2), 2 ) ) / 2, 2 )
      n = modulo( ipos(1), 2 )

   end subroutine getIndices
!
!==========================================================================================
!
   subroutine writeLSDRas( lsdef, nlsdef, imod )
!             ===========    i      i     i
!
!     Write the deterministic defect file RAS in rasmol format.
!
!  Arguments:
!
      type(LSDefect_type), intent(in) :: lsdef(:)   ! 
      integer(INT4), intent(in) :: nlsdef           ! #rasmol points
      integer(INT4), intent(in) :: imod             ! what to write
                                          ! imod=1, write only defect lattice sites
                                          ! imod=2, write all latice sites
                                          ! imod=3, write all atom sites
!
!  Local Variables:
!
      integer(INT4) :: ilsdef, isa
!
!  Start of Executable Code:
!  ------------------------
!
      if( imod < 1 .or. imod > 3 )  error stop 'writeAtom: error imod'
      if( nlsdef == 0 )  return

      call openFile( 413, 'RAS', 'formatted', 'write' )
!          ========  i       i            i         i

      write( 413, * ) nlsdef
      write( 413, * ) '' 

      do ilsdef = 1, size( lsdef )
         select case( imod )
            case( 1 )       ! write those only with type 1
               write( 413, * ) 'O', lsdef(ilsdef)%ipos(:)
            case( 2 )    ! write all ipos
               if( lsdef( ilsdef )%type .ne. 0 ) then
                  write( 413, * ) 'O', lsdef(ilsdef)%ipos(:)
               end if
            case default
               do isa = 1, lsdef(ilsdef)%nsa
                  write( 413, * ) 'O', lsdef(ilsdef)%dpos(3,isa) + lsdef(ilsdef)%ipos(:)
               end do
         end select
      end do

      close( 413 )

   end subroutine writeLSDRas
!
!=============================================================================================
!
   subroutine readDefects( lsdef, nlsdef, nlsstr, icolmin, icolmax )
!             ===========   io      o      o        o          o
!
!     Read the defects from the file deffile and initialize the lsdef structure.
!     The units are 1/4 of the lattice constant!
!
!  Arguments:
!
      type(LSDefect_type), allocatable, intent(inout) :: lsdef(:)  
      integer(INT4), intent(out) :: nlsdef             ! # lsdefects
      integer(INT4), intent(out) :: nlsstr             ! # lsdefects
      integer(INT4), intent(out) :: icolmin(3)
      integer(INT4), intent(out) :: icolmax(3)
!
!  Local Variables:
!
      integer(INT4) ::  status, ndef, ilsdef, isa, ielem
      integer(INT4) :: i, j, k, imin, imax, jmin, jmax, kmin, kmax
!
!  Start of Executable Code:
!  ------------------------
!
      call openFile( 41, deffile, 'formatted', 'read' )
!          ========  i      i          i         i
!
!  Read the boundaries of the column
!
      read( 41, * ) icolmin(:), icolmax(:)
      dam%wcol(:) = icolmax(:) - icolmin(:)
!
!  Allocate lsdef to cover all lattice sites in the simulation cell
!
      ndef = product( dam%wcol(:) ) / 8  ! contains all lattice sites in cell 
!
!  Read the lattice site defects
!
      read( 41, * ) nlsdef
      nlsstr = nlsdef

      if( nlsdef == 0 ) return

      allocate( lsdef( 1:ndef ), stat=status )      ! No additional strain in this case
      if( status /= 0 )  error stop 'readDefects: lsdef not allocated'

      do ilsdef = 1, nlsdef
         read( 41, * ) lsdef(ilsdef)%ipos(:)
         read( 41, * ) lsdef(ilsdef)%type
         if( lsdef(ilsdef)%type == -1 ) then
            read( 41, * ) lsdef(ilsdef)%nsa
            do isa = 1, lsdef(ilsdef)%nsa
               read( 41, * ) lsdef(ilsdef)%ielem(isa), lsdef(ilsdef)%dpos(:,isa)
            end do
         elseif( lsdef(ilsdef)%type > 0 ) then
            lsdef(ilsdef)%nsa = getDefnAtoms( lsdef(ilsdef)%type )
!                               ============           i
            read( 41, *, iostat=status ) lsdef(ilsdef)%ior,  &
                                         lsdef(ilsdef)%ielem(1:lsdef(ilsdef)%nsa)
            if( status /= 0 )  error stop 'readDefects: error during reading deffile'
         end if
         if( lsdef(ilsdef)%nsa > NATOMM )  error stop 'readDefects: nsa exceeds limits'
         do isa = 1, lsdef(ilsdef)%nsa
            do ielem = nelem, 1, -1
               if( elem(ielem)%z == lsdef(ilsdef)%ielem(isa) )  exit
            end do
            if( ielem == 0 )  error stop 'readDefects: unknown atom type'
            lsdef(ilsdef)%ielem(isa) = ielem
         end do
                                               ! detect column size (optional)
         i = ( lsdef(ilsdef)%ipos(1) - modulo( lsdef(ilsdef)%ipos(1), 4 ) ) / 4
         j = ( lsdef(ilsdef)%ipos(2) - modulo( lsdef(ilsdef)%ipos(2), 4 ) ) / 4
         k = ( lsdef(ilsdef)%ipos(3) - modulo( lsdef(ilsdef)%ipos(3), 4 ) ) / 4
         if( ilsdef == 1 ) then
            imin = i
            imax = i
            jmin = j
            jmax = j
            kmin = k
            kmax = k
         else
            imin = min( imin, i )
            imax = max( imax, i )
            jmin = min( jmin, j )
            jmax = max( jmax, j )
            kmin = min( kmin, k )
            kmax = max( kmax, k )
         end if
      end do

   end subroutine readDefects
!
!===========================================================================================
!
   subroutine writeDefects( lsdef, nlsdef, icolmin, icolmax )
!             ============    i      i       i        i
!
!     Write defects into the deterministic defect file 'DEF'.
!
!  Arguments:
!
      type(LSDefect_type), intent(in) :: lsdef(:)   ! 
      integer(INT4), intent(in) :: nlsdef           ! #rasmol points
      integer(INT4), intent(in) :: icolmin(3)
      integer(INT4), intent(in) :: icolmax(3)
!
!  Local Variables:
!
      integer(INT4) :: ilsdef, isa
!
!  Start of Executable Code:
!  ------------------------
!
      if( nlsdef == 0 ) return

      call openFile( 43, 'DEF', 'formatted', 'write' )
!          ========  i       i            i         i

      write( 43, '(6I8)' ) icolmin(1:3), icolmax(1:3) !  Bounds of the column
      write( 43, * ) nlsdef

      do ilsdef = 1, size( lsdef )
         if( lsdef(ilsdef)%type == 0 )  cycle         ! skip nondefect sites
         write( 43, * ) lsdef(ilsdef)%ipos(:)
         write( 43, * ) lsdef(ilsdef)%type
         if( lsdef(ilsdef)%type == -1 ) then
            write( 43, * ) lsdef(ilsdef)%nsa
            do isa = 1, lsdef(ilsdef)%nsa
               write( 43, * ) elem(lsdef(ilsdef)%ielem(isa))%z, lsdef(ilsdef)%dpos(:,isa)
            end do
         else
            write( 43, * ) lsdef(ilsdef)%ior, elem(lsdef(ilsdef)%ielem(1:lsdef(ilsdef)%nsa))%z
         end if
      end do

      close( 43 )

   end subroutine writeDefects
!
!=============================================================================================
!
end module Eval_Det_mod
