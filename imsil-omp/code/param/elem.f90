module Param_Elem_mod
!
!  Hold data about atom species
!
!  Used Modules:
!  ============
!
   use Const_mod,          only: REAL8, INT4,  &
                                 ZERO, QUART, TENTH, BIG
   use Array_mod,          only: check_allocation
   use Undefine_mod,       only: undefine, defined
   use Open_mod,           only: openFile

   use Param_Setup_mod,    only: nelem, nelemIon, temp

   use Math_Debye_mod,     only: debye1
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
   public :: readAtoms              ! read the &atoms record
!            =========
   public :: initAtoms              ! initialize the elem record
!            =========
   public :: writeAtoms             ! write the elem record to output file
!            ==========
   public :: anaMolecule            ! analyse molecule name and return atom name, atomic
!            ===========            !    number, and abundance
!
!  Public data:
!  ===========
!
   type :: Elem_type
!          =========
      character(2) :: name          ! element name (chemical abbreviation)
      integer(INT4) :: z            ! atomic number ("Z")
      real(REAL8) :: m              ! atomic mass (atomic mass unit)
      real(REAL8) :: dens           ! atomic density (cm-3)
      real(REAL8) :: tdebye         ! Debye temperature (K) of the pure material
      real(REAL8) :: vib            ! thermal vibration amplitude (A) (currently unused)
      real(REAL8) :: esurf          ! surface binding energy (eV) of the pure material
      real(REAL8) :: efake(2)       ! energies for fake thermal vibrations (eV)
      real(REAL8) :: tempfake(2)    ! fake temperatures (K)
      logical :: volatile           ! flag indicating element disappears when stopped
   end type Elem_type

   type(Elem_type), allocatable, public :: elem(:)   ! (nelem) atom array
!
!  Private data:
!  ============
!
   type :: Default_Elem_type           ! chemical element data type
!          =================
      character(2) :: name
      real(REAL8) :: mass1
      real(REAL8) :: mass2
      real(REAL8) :: dens
      character(12) :: long_name
      real(REAL8) :: tdebye
      real(REAL8) :: esurf
   end type Default_Elem_type

   type(Default_Elem_type) :: default_elem(92)

   data default_elem( 1: 2) / &
      Default_Elem_type( 'H ',  1.0078, 1.008,  4.271e22,'HYDROGEN    ', 110. , 0.   ) , &
      Default_Elem_type( 'HE',  4.0026, 4.0026, 1.894e22,'HELIUM      ',   0. , 0.   ) /
   data default_elem( 3:10) / &
      Default_Elem_type( 'LI',  7.016,  6.939,  4.597e22,'LITHIUM     ', 344. , 1.67 ) , &
      Default_Elem_type( 'BE',  9.012,  9.0122,12.046e22,'BERYLLIUM   ',1440. , 3.38 ) , &
      Default_Elem_type( 'B ', 11.009, 10.811, 13.093e22,'BORON       ',1250. , 5.73 ) , &
      Default_Elem_type( 'C ', 12.   , 12.011, 11.364e22,'CARBON      ',2230. , 7.41 ) , &
      ! density of N adjusted to give the correct density of Si3N4 (real value: 3.481)
      Default_Elem_type( 'N ', 14.003, 14.007, 22.590e22,'NITROGEN    ',   0. , 0.   ) , &
      ! density of O adjusted to give the correct density of SiO2 (real value: 4.302)
      Default_Elem_type( 'O ', 15.995, 15.999,  8.759e22,'OXYGEN      ',   0. , 0.   ) , &
      Default_Elem_type( 'F ', 18.998, 18.998,  3.522e22,'FLUORINE    ',   0. , 0.   ) , &
      Default_Elem_type( 'NE', 19.992, 20.183,  3.585e22,'NEON        ',  75. , 0.   ) /
   data default_elem(11:18) / &
      Default_Elem_type( 'NA', 22.99 , 22.989,  2.541e22,'SODIUM      ', 158. , 1.12 ) , &
      Default_Elem_type( 'MG', 23.985, 24.312,  4.302e22,'MAGNESIUM   ', 400. , 1.54 ) , &
      Default_Elem_type( 'AL', 26.982, 26.981,  6.023e22,'ALUMINUM    ', 428. , 3.36 ) , &
      Default_Elem_type( 'SI', 27.977, 28.086,  4.994e22,'SILICON     ', 493.6, 4.70 ) , &
      Default_Elem_type( 'P ', 30.994, 30.973,  3.542e22,'PHOSPHORUS  ',   0. , 3.27 ) , &
      Default_Elem_type( 'S ', 31.97 , 32.064,  3.885e22,'SULFUR      ',   0. , 2.88 ) , &
      Default_Elem_type( 'CL', 34.969, 35.453,  3.220e22,'CHLORINE    ',   0. , 0.   ) , &
      Default_Elem_type( 'AR', 39.98 , 39.948,  2.488e22,'ARGON       ',  92. , 0.   ) /
   data default_elem(19:36) / &
      Default_Elem_type( 'K ', 38.96 , 39.102,  1.329e22,'POTASSIUM   ',  91. , 0.93 ) , &
      Default_Elem_type( 'CA', 39.96 , 40.08 ,  2.014e22,'CALCIUM     ', 230. , 1.83 ) , &
      Default_Elem_type( 'SC', 44.956, 44.956,  4.015e22,'SCANDIUM    ', 360. , 3.49 ) , &
      Default_Elem_type( 'TI', 47.95 , 47.9  ,  5.682e22,'TITANIUM    ', 420. , 4.89 ) , &
      Default_Elem_type( 'V ', 50.94 , 50.942,  7.213e22,'VANADIUM    ', 380. , 5.33 ) , &
      Default_Elem_type( 'CR', 51.94 , 51.996,  8.330e22,'CHROMIUM    ', 630. , 4.12 ) , &
      Default_Elem_type( 'MN', 54.94 , 54.938,  8.150e22,'MANGANESE   ', 410. , 2.98 ) , &
      Default_Elem_type( 'FE', 55.94 , 55.847,  8.483e22,'IRON        ', 470. , 4.34 ) , &
      Default_Elem_type( 'CO', 58.93 , 58.933,  8.989e22,'COBALT      ', 445. , 4.43 ) , &
      Default_Elem_type( 'NI', 57.94 , 58.71 ,  9.125e22,'NICKEL      ', 450. , 4.46 ) , &
      Default_Elem_type( 'CU', 62.93 , 63.54 ,  8.483e22,'COPPER      ', 343. , 3.52 ) , &
      Default_Elem_type( 'ZN', 63.93 , 65.37 ,  6.546e22,'ZINC        ', 327. , 1.35 ) , &
      Default_Elem_type( 'GA', 68.93 , 69.72 ,  5.104e22,'GALLIUM     ', 320. , 2.82 ) , &
      Default_Elem_type( 'GE', 73.92 , 72.59 ,  4.428e22,'GERMANIUM   ', 374. , 3.88 ) , &
      Default_Elem_type( 'AS', 74.92 , 74.922,  4.597e22,'ARSENIC     ', 282. , 1.26 ) , &
      Default_Elem_type( 'SE', 79.92 , 78.96 ,  3.650e22,'SELENIUM    ',  90. , 2.14 ) , &
      Default_Elem_type( 'BR', 78.92 , 79.909,  2.562e22,'BROMINE     ',   0. , 0.   ) , &
      Default_Elem_type( 'KR', 83.92 , 83.8  ,  1.870e22,'KRYPTON     ',  72. , 0.   ) /
   data default_elem(37:54) / &
      Default_Elem_type( 'RB', 84.91 , 85.47 ,  1.077e22,'RUBIDIUM    ',  56. , 0.86 ) , &
      Default_Elem_type( 'SR', 87.91 , 87.62 ,  1.787e22,'STRONTIUM   ', 147. , 1.70 ) , &
      Default_Elem_type( 'Y ', 89.91 , 88.905,  3.041e22,'YTTRIUM     ', 280. , 4.24 ) , &
      Default_Elem_type( 'ZR', 89.9  , 91.22 ,  4.271e22,'ZIRCONIUM   ', 291. , 6.33 ) , &
      Default_Elem_type( 'NB', 92.91 , 92.906,  5.576e22,'NIOBIUM     ', 275. , 7.59 ) , &
      Default_Elem_type( 'MO', 95.5  , 95.94 ,  6.407e22,'MOLYBDENUM  ', 450. , 6.83 ) , &
      Default_Elem_type( 'TC',  0.   ,  99.  ,  0.      ,'TECHNETIUM  ', 453. , 0.   ) , &
      Default_Elem_type( 'RU',101.9  ,101.07 ,  7.256e22,'RUTHENIUM   ', 600. , 6.69 ) , &
      Default_Elem_type( 'RH',102.9  ,102.91 ,  7.256e22,'RHODIUM     ', 480. , 5.78 ) , &
      Default_Elem_type( 'PD',105.9  ,106.4  ,  6.767e22,'PALLADIUM   ', 274. , 3.91 ) , &
      Default_Elem_type( 'AG',106.9  ,107.87 ,  5.847e22,'SILVER      ', 225. , 2.97 ) , &
      Default_Elem_type( 'CD',113.9  ,112.4  ,  4.597e22,'CADMIUM     ', 209. , 1.16 ) , &
      Default_Elem_type( 'IN',114.9  ,114.82 ,  3.836e22,'INDIUM      ', 108. , 2.49 ) , &
      Default_Elem_type( 'SN',119.9  ,118.69 ,  3.695e22,'TIN         ', 200. , 3.12 ) , &
      Default_Elem_type( 'SB',120.9  ,121.75 ,  3.273e22,'ANTIMONY    ', 211. , 2.72 ) , &
      Default_Elem_type( 'TE',130.   ,127.6  ,  2.938e22,'TELLURIUM   ', 153. , 2.02 ) , &
      Default_Elem_type( 'J ',126.9  ,126.9  ,  2.343e22,'IODINE      ',   0. , 0.   ) , &
      Default_Elem_type( 'XE',129.5  ,131.3  ,  1.403e22,'XENON       ',  64. , 0.   ) /
   data default_elem(55:86) / &
      Default_Elem_type( 'CS',134.   ,132.91 ,  0.860e22,'CAESIUM     ',  38. , 0.81 ) , &
      Default_Elem_type( 'BA',138.   ,137.34 ,  1.544e22,'BARIUM      ', 110. , 1.84 ) , &
      Default_Elem_type( 'LA',139.   ,138.91 ,  2.676e22,'LANTHANUM   ', 142. , 4.42 ) , &
      Default_Elem_type( 'CE',140.   ,140.12 ,  2.868e22,'CERIUM      ',   0. , 4.23 ) , &
      Default_Elem_type( 'PR',141.   ,140.91 ,  2.895e22,'PRASEODYMIUM',   0. , 3.71 ) , &
      Default_Elem_type( 'ND',142.   ,144.24 ,  2.932e22,'NEODYMIUM   ',   0. , 3.28 ) , &
      Default_Elem_type( 'PM',  0.   ,147.   ,  0.      ,'PROMETHIUM  ',   0. , 0.   ) , &
      Default_Elem_type( 'SM',152.   ,150.35 ,  3.026e22,'SAMARIUM    ', 166. , 2.16 ) , &
      Default_Elem_type( 'EU',153.   ,151.96 ,  2.084e22,'EUROPIUM    ',   0. , 1.85 ) , &
      Default_Elem_type( 'GD',157.9  ,157.25 ,  3.026e22,'GADOLINIUM  ', 200. , 3.57 ) , &
      Default_Elem_type( 'TB',159.9  ,158.92 ,  3.136e22,'TERBIUM     ',   0. , 3.81 ) , &
      Default_Elem_type( 'DY',164.   ,162.5  ,  3.170e22,'DYSPROSIUM  ', 210. , 2.89 ) , &
      Default_Elem_type( 'HO',165.   ,164.93 ,  3.220e22,'HOLMIUM     ',   0. , 3.05 ) , &
      Default_Elem_type( 'ER',166.   ,167.26 ,  3.273e22,'ERBIUM      ',   0. , 3.05 ) , &
      Default_Elem_type( 'TM',169.   ,168.93 ,  3.327e22,'THULIUM     ',   0. , 2.52 ) , &
      Default_Elem_type( 'YB',174.   ,173.04 ,  2.428e22,'YTTERBIUM   ', 120. , 1.74 ) , &
      Default_Elem_type( 'LU',175.   ,174.97 ,  3.383e22,'LUTETIUM    ', 210. , 4.29 ) , &
      Default_Elem_type( 'HF',180.   ,178.49 ,  4.428e22,'HAFNIUM     ', 252. , 6.31 ) , &
      Default_Elem_type( 'TA',181.   ,180.95 ,  5.525e22,'TANTALUM    ', 240. , 8.10 ) , &
      Default_Elem_type( 'W ',184.   ,183.85 ,  6.320e22,'TUNGSTEN    ', 400. , 8.68 ) , &
      Default_Elem_type( 'RE',187.   ,186.2  ,  6.805e22,'RHENIUM     ', 430. , 8.09 ) , &
      Default_Elem_type( 'OS',192.   ,190.2  ,  7.144e22,'OSMIUM      ', 500. , 8.13 ) , &
      Default_Elem_type( 'IR',193.   ,192.2  ,  7.052e22,'IRIDIUM     ', 420. , 6.90 ) , &
      Default_Elem_type( 'PT',195.   ,195.09 ,  6.618e22,'PLATINUM    ', 240. , 5.86 ) , &
      Default_Elem_type( 'AU',197.   ,196.97 ,  5.904e22,'GOLD        ', 165. , 3.80 ) , &
      Default_Elem_type( 'HG',202.   ,200.59 ,  4.069e22,'MERCURY     ',  71.9, 0.64 ) , &
      Default_Elem_type( 'TL',205.   ,204.37 ,  3.501e22,'THALLIUM    ',  78.5, 1.88 ) , &
      Default_Elem_type( 'PB',208.   ,207.19 ,  3.291e22,'LEAD        ', 105. , 2.03 ) , &
      Default_Elem_type( 'BI',209.   ,208.98 ,  2.827e22,'BISMUTH     ', 119. , 2.17 ) , &
      Default_Elem_type( 'PO',210.   ,210.   ,  2.653e22,'POLONIUM    ',   0. , 1.50 ) , &
      Default_Elem_type( 'AT',210.   ,210.   ,  0.      ,'ASTATINE    ',   0. , 0.   ) , &
      Default_Elem_type( 'RN',222.   ,222.   ,  0.      ,'RADON       ',   0. , 0.   ) /
   data default_elem(87:92) / &
      Default_Elem_type( 'FR',223.   ,223.   ,  0.      ,'FRANCIUM    ',   0. , 0.   ) , &
      Default_Elem_type( 'RA',226.   ,226.   ,  1.338e22,'RADIUM      ',   0. , 0.   ) , &
      Default_Elem_type( 'AC',227.   ,227.   ,  0.      ,'ACTINIUM    ',   0. , 0.   ) , &
      Default_Elem_type( 'TH',232.   ,232.   ,  3.026e22,'THORIUM     ', 163. , 5.93 ) , &
      Default_Elem_type( 'PA',231.   ,231.   ,  4.015e22,'PROTACTINIUM', 0.   , 0.   ) , &
      Default_Elem_type( 'U ',238.04 ,238.04 ,  4.818e22,'URANIUM     ', 207. , 5.42 ) /
!
!  Program History:
!  ===============
!
!  26-Nov-2002  G.Hobler  creation
!  21-Feb-2017  G.Hobler  read &atoms record
!
!=======================================================================
!
contains
!
!=======================================================================
!
   subroutine readAtoms( )
!             =========
!
!     Read the &atoms record from the input file and initialize
!     the elem structure. Undefined parameters will be set somewhere
!     else.
!
!  Arguments:
!
!     None
!
!  Local Variables:
!
      integer(INT4) :: atom                     ! index parameters

      character(2) :: name

      real(REAL8) :: density                    ! atomic density (cm-3) of the pure material
      real(REAL8) :: efake(2)                   ! energies for fake thermal vibrations (keV)
      real(REAL8) :: esurf                      ! surface binding energy (eV) of the pure
                                                !    material
      real(REAL8) :: mass                       ! atom mass (atomic mass units)
      real(REAL8) :: tdebye                     ! Debye termperature (K) of the pure material
      real(REAL8) :: tempfake(2)                ! fake temperatures (K)
      real(REAL8) :: vib                        ! thermal vibration amplitude (A)

      namelist /atoms/  atom,  &
                        name,  &
                        density, efake, esurf, mass, tdebye, tempfake, vib

      character(256) :: msg, errmsg
      integer(INT4) :: istat, ios, ielem
!
!  Start of Executable Code:
!  ------------------------
!
!  Allocate elem array
!
      allocate( elem( nelem ), stat=istat, errmsg=errmsg )
      call check_allocation( 'readAtoms', 'elem', istat, errmsg )
!
!  Set default values
!
      do ielem = 1, nelem
         call undefine( elem(ielem)%name )
         call undefine( elem(ielem)%z )
         call undefine( elem(ielem)%m )
         call undefine( elem(ielem)%dens )
         call undefine( elem(ielem)%efake(:) )
         call undefine( elem(ielem)%esurf )
         call undefine( elem(ielem)%tdebye )
         call undefine( elem(ielem)%tempfake(:) )
         call undefine( elem(ielem)%vib )
         elem(ielem)%volatile = .false.
      end do
!
!  Loop over reading records
!  -------------------------
!
      call openFile( 11, 'INP', 'formatted', 'read' )
!          ========  i     i         i         i

      do
!
!  Set default values of namelist variables
!
         call undefine( atom )
         call undefine( name )
         call undefine( mass )
         call undefine( density )
         call undefine( efake(:) )
         call undefine( esurf )
         call undefine( tdebye )
         call undefine( tempfake(:) )
         call undefine( vib )
!
!  Read parameters from the snpar namlist
!
         read( 11, nml=atoms, iostat=ios, iomsg=msg )
         if( ios < 0  )  exit
         if( ios > 0 ) then
            print '(/a)',  trim( msg )
            goto 900
         end if
!
!  Check index parameter
!
         if( .not. defined( atom ) )  goto 920
         if( atom < 1 .or. atom > nelem )  goto 921
         ielem = atom
!
!  Copy parameters to elem structure
!
         if( defined( name ) )  elem(ielem)%name = name
         if( defined( mass ) )  elem(ielem)%m = mass
         if( defined( density ) )  elem(ielem)%dens = density / 1e24          ! 1/cm3 -> 1/A^3
         where( defined(efake(:)) )  elem(ielem)%efake(:) = efake(:) * 1000   ! keV -> eV
         if( defined( esurf ) )  elem(ielem)%esurf = esurf
         if( defined( tdebye ) )  elem(ielem)%tdebye = tdebye
         where( defined(tempfake(:)) )  elem(ielem)%tempfake(:) = tempfake(:)
         if( defined( vib ) )  elem(ielem)%vib = vib

      end do

      close( 11 )

      return
!
!  Error messages
!
900   error stop '&atoms: error during reading namelist'

920   error stop '&atoms: atom must be specified'
921   error stop '&atoms: 1<=atom<=natom required'

   end subroutine readAtoms
!
!=======================================================================
!
   subroutine initAtoms( )
!             =========
!
!     Initialize the elem structure
!
!     Note that the elem structure may have defined elements on entry from readAtoms, readIon,
!     and readAtomsFromCell
!
!  Arguments:
!
!     None
!
!  Program History:
!
!      7-Mar-2017  G.Hobler  creation
!
!  Local Variables:
!
      integer(INT4) :: ielem, iz
      real(REAL8) :: x, u1
!
!  Start of Executable Code:
!  ------------------------
!
!  Set default values
!
      do ielem = 1, nelem

         if( .not.defined( elem(ielem)%name ) .and. .not.defined( elem(ielem)%z ) )  goto 922

         if( defined( elem(ielem)%z ) ) then
            iz = elem(ielem)%z
         else
            do iz = 1, 92
               if( elem(ielem)%name == default_elem(iz)%name )  exit
            end do
            if( iz > 92 )  goto 923
            elem(ielem)%z = iz
         end if

         if( .not.defined( elem(ielem)%name ) )  elem(ielem)%name = default_elem(iz)%name

         if( .not.defined( elem(ielem)%m ) ) then
            if( ielem <= nelemIon ) then
               elem(ielem)%m = default_elem(iz)%mass1
            else
               elem(ielem)%m = default_elem(iz)%mass2
            end if
         end if

         if( .not.defined( elem(ielem)%dens ) )  elem(ielem)%dens = default_elem(iz)%dens / 1e24
         if( .not.defined( elem(ielem)%tdebye ) )  elem(ielem)%tdebye = default_elem(iz)%tdebye
         if( .not.defined( elem(ielem)%esurf ) )  &
               elem(ielem)%esurf = max( TENTH, default_elem(iz)%esurf )

      end do
!
!  Check validity of the parameters
!
      do ielem = 1, nelem

         if( elem(ielem)%m <= ZERO )  goto 924
         if( elem(ielem)%dens <= ZERO )  goto 925
         if( elem(ielem)%tdebye < ZERO )  goto 926
         if( elem(ielem)%esurf <= ZERO )  goto 927
         if( defined( elem(ielem)%tempfake(2) ) .and.  &
             .not. defined( elem(ielem)%tempfake(1) ) )  goto 928
         if( all( defined( elem(ielem)%tempfake(:) ) ) .and.  &
             .not. all( defined( elem(ielem)%efake(:) ) ) )  goto 929

      end do
!
!  Initialize derived parameters
!
      do ielem = 1, nelem
         if( .not.defined( elem(ielem)%vib ) ) then
            if( elem(ielem)%tdebye == 0. ) then
               elem(ielem)%vib = 0.083
            else
               x = elem(ielem)%tdebye / temp
               u1 = 12.1 * sqrt( (debye1( x )/x + QUART) /  &
!                                 ======  i
                                 (elem(ielem)%tdebye * elem(ielem)%m) )
               if( u1 > BIG ) then
                  elem(ielem)%vib = 0.083
               else
                 elem(ielem)%vib = u1
               end if
            end if
         end if
      end do

      return
!
!  Error messages
!
922   error stop '&atoms: name or atomic number must be specified'
923   error stop '&atoms: unknown name'
924   error stop '&atoms: mass>0 required'
925   error stop '&atoms: density>0 required'
926   error stop '&atoms: tdebye>=0 required'
927   error stop '&atoms: esurf>0 required'
928   error stop '&atoms: tempfake(1) required if tempfake(2) specified'
929   error stop '&atoms: efake(1:2) required if tempfake(1:2) specified'

   end subroutine initAtoms
!
!=======================================================================
!
   subroutine writeAtoms( )
!             ==========
!
!     Write data derived from &atoms parameters to the output file
!
!  Arguments:
!
!     None
!
!  Program History:
!
!      6-Mar-2017  G.Hobler  creation
!
!  Local Variables:
!
      logical :: anydef
      integer(INT4) :: i, ielemmin, ielemmax, ielem, ifake
!
!  Start of Executable Code:
!  ------------------------
!
!  Open the output file
!
      call openFile( 21, 'OUT', 'formatted', 'append' )
!          ========  i     i         i          i
!
!  Make list of atoms with their properties and abundance
!  (cannot be done in Param_Elem_mod since this would require recursive use)
!
      write( 21, '(t2,a,t13,a)' )  &
         '&atoms:', 'Properties of atoms and pure materials used in the simulation:'

      do i = 0, (nelem-1)/8

         ielemmin = 8*i + 1
         ielemmax = min( 8*i + 8, nelem )

         write( 21, '(t13,a)', advance='no' ) '+------------+'
         write( 21, '(8a)' ) &
               ('------+', ielem=ielemmin,ielemmax)

         write( 21, '(t13,a, t26,a)', advance='no' ) '|', '|'
         write( 21, '(8(a,i2,a))' ) &
               ('Atom', ielem, '|', ielem=ielemmin,ielemmax)

         write( 21, '(t13,a, t26,a)', advance='no' ) '|', '|'
         write( 21, '(8(a,a2,a))' ) &
               (' (', elem(ielem)%name, ') |', ielem=ielemmin,ielemmax)

         write( 21, '(t13,a)', advance='no' ) '+------------+'
         write( 21, '(8a)' ) &
               ('------+', ielem=ielemmin,ielemmax)

         write( 21, '(t13,a)', advance='no' ) '|     Z      |'
         write( 21, '(8(i4,a))' ) &
               (elem(ielem)%z, '  |', ielem=ielemmin,ielemmax)

         write( 21, '(t13,a)', advance='no' ) '|     M      |'
         write( 21, '(8(f6.2,a))' ) &
               (elem(ielem)%m, '|', ielem=ielemmin,ielemmax)

         write( 21, '(t13,a)', advance='no' ) '|density/1e22|'
         write( 21, '(8(f6.3,a))' ) &
               (elem(ielem)%dens*1e24/1e22, '|', ielem=ielemmin,ielemmax)

         write( 21, '(t13,a)', advance='no' ) '| Debye temp |'
         write( 21, '(8(f6.1,a))' ) &
               (elem(ielem)%tdebye, '|', ielem=ielemmin,ielemmax)

         write( 21, '(t13,a)', advance='no' ) '|    esurf   |'
         write( 21, '(8(f6.3,a))' ) &
               (elem(ielem)%esurf, '|', ielem=ielemmin,ielemmax)

         do ifake = 1, 2
            anydef = .false.
            do ielem = ielemmin, ielemmax
               anydef = anydef .or. defined( elem(ielem)%efake(ifake) )
            end do
            if( anydef ) then
               write( 21, '(t13,a,i1,a)', advance='no' ) '|  efake(', ifake, ')  |'
               do ielem = ielemmin, ielemmax
                  write( 21, '(i6,a)', advance='no' )  int(elem(ielem)%efake(ifake)/1000), '|'
               end do
               write( 21, '()' )
            end if
         end do

         do ifake = 1, 2
            anydef = .false.
            do ielem = ielemmin, ielemmax
               anydef = anydef .or. defined( elem(ielem)%tempfake(ifake) )
            end do
            if( anydef ) then
               write( 21, '(t13,a,i1,a)', advance='no' ) '| tempfake(', ifake, ')|'
               do ielem = ielemmin, ielemmax
                  write( 21, '(f6.1,a)', advance='no' )  elem(ielem)%tempfake(ifake), '|'
               end do
               write( 21, '()' )
            end if
         end do

         write( 21, '(t13,a)', advance='no' ) '+------------+'
         write( 21, '(8a)' ) &
               ('------+', ielem=ielemmin,ielemmax)
      end do

      write( 21, '()' )

      close( 21 )

   end subroutine writeAtoms
!
!=======================================================================
!
   subroutine anaMolecule( molec, ele, iz, nabund )
!             ===========   io     o   o     o
!
!     Analyze the molecule name and extract the name, atomic number, and
!     abundance of the first atom referenced. Return the molecule name minus
!     the first atom.
!
!  Arguments:
!
      character(*), intent(inout) :: molec      ! chemical name of the molecule
      character(2), intent(out) :: ele          ! chemical name of the atom
      integer(INT4), intent(out) :: iz          ! atomic number
      integer(INT4), intent(out) :: nabund      ! abundance of atom in molecule
!
!  Local Variables:
!
      integer(INT4) :: k, i

      character(52) :: alf = 'ABCDEFGHIJKLMNOPQRSTUVWXYZabcdefghijklmnopqrstuvwxyz'
      character(10) :: num = '1234567890'
!
!  Start of Executable Code:
!  ------------------------
!
      ele = ' '
      iz = 0
      nabund = 0
!
!  Skip leading blanks
!
      do k = 1, len(molec)
         if( molec(k:k) /= ' '  )  exit
      enddo
      if( k > len(molec) ) return
!
!  Get element name
!
      i = index( alf, molec(k:k) )
      if( i == 0 )  return                        ! no letter
      if( i > 26 )  i = i - 26                    ! convert to uppercase
      ele = alf(i:i)

      k = k + 1                                 ! test second letter
      if( k > len(molec) ) return
      i = index( alf, molec(k:k) )
      if( i > 0 ) then                           ! if there is a second letter, assign it
         if( i > 26)  i = i - 26                 !    to the element name preliminarly
         ele(2:2) = alf(i:i)
         k = k + 1
      endif
!
!  Search for the element name in the list
!
      do i = 1, 92
         if( ele == default_elem(i)%name )  exit
      enddo
                                                ! the second letter could belong to the next
                                                !    element; if element has not been found,
      if( i > 92 .and. ele(2:2) /= ' ' ) then    !    try one-character name
         ele(2:2) = ' '
         k = k - 1
         do i = 1, 92
            if( ele == default_elem(i)%name )  exit
         enddo
      endif

      if( i > 92 )  error stop 'invalid ion species or target material'

      iz = i                                    ! element name found => ele, iz
!
!  Get abundance of the atom in the molecule
!
      if( k > len(molec) ) then                  ! if at end of molecule name ...
         nabund = 1
      else if( index(num,molec(k:k)) == 0 ) then ! if next character is a letter ...
         nabund = 1
      else                                      ! otherwise (next character a digit) ...
         do
            read( molec(k:k),'(i1)' ) i
            nabund = 10 * nabund + i
            k = k + 1
            if( k > len(molec) )  exit
            if( index( num, molec(k:k) ) == 0 )  exit
         end do
      end if

      molec = molec(k:len(molec))

   end subroutine anaMolecule
!
!=======================================================================
!
end module Param_Elem_mod
