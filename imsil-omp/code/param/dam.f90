module Param_Dam_mod
!
!  Holds data on damage models, handles their input and output, and calulates
!  derived quantities
!
!  Used Modules:
!  ============
!
   use Const_mod,          only: INT4, REAL8,  &
                                 ZERO, ONE, HALF, THIRD, PI, BIG, EPSMB,  &
                                 ENUM
   use String_mod,         only: str, to_lower
   use Undefine_mod,       only: undefine, defined
   use Open_mod,           only: openFile
   use Param_WriteTab_mod, only: writeTab

   use Param_Setup_mod,    only: nelem, ndim_setup => ndim, nion, nionhis, nr,  &
                                 nthread,  &
                                 ldamdyn, ldamdyn_set, ldyn, lsput, usebk, useion,  &
                                 usehis, usedef, nelemIon, usecell
   use Param_Geom_mod,     only: geom
   use Param_Elem_mod,     only: elem
   use Param_Ion_mod,      only: ion
   use Param_Xtal_mod,     only: xtal
   use Param_Mat_mod,      only: mat, irxtal
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
   public :: readDam                ! read the damage parameters
!            =======
   public :: initDam                ! initialize the dam structure
!            =======
   public :: writeDam               ! write damage parameters to output file
!            ========
!
!  Public Data:
!  ===========
!
   type :: Dam_type

      logical :: lamzon                      ! amorphous zone model flag
      logical :: lamo2                       ! 2-d displacement of amorphous atoms flag
      logical :: lcap                        ! capture radius model flag
      logical :: lcasc                       ! cascade statistics flag
      logical :: ldam                        ! damage calculation and consideration flag
      logical :: ldet                        ! deterministic defects flag
      logical :: lflow                       ! recoil flow flag
      logical :: libiia                      ! ion beam induced interfacial amorphization flag
      logical :: lpock                       ! amorphous pocket model flag
      logical :: lrcoil                      ! recoil cascade flag
      logical :: lrcoil1                     ! recoil cascade before reentering flag
      logical :: lrcoil2                     ! recoil cascade after reentering flag
      logical :: lrepl                       ! replacement collision flag
      logical :: lstrmd                      ! use strain from MD simulations for deterministic
                                             !   defects
      logical :: reenter                     ! reenter flag for recoils

      integer(INT4) :: ndamdim               ! dimensionality of histogram used for damage
                                             !    consideration
      integer(INT4) :: relaxmethod           ! numerical method for target relaxation in a
                                             !    dynamic simulation
      integer(INT4) :: relaxmodel            ! physical model for target relaxation in a dynamic
                                             !    simulation
      integer(INT4) :: wcol(1:3)             ! column width in x, y and z direction
                                             ! used when creating the defect file

      real(REAL8) :: cospsi110               ! cos(critical angle wrt <110> below which
                                             !    elimkp is not applied)
      real(REAL8) :: damamo                  ! amorphization damage concentration [1/A^3]
      real(REAL8) :: damsat                  ! damage saturation concentration [1/A^3]
                                             !    (undefined if not to be used)
      real(REAL8) :: dens                    ! reference atom density [1/A^3]
      real(REAL8) :: dibiia                  ! spatial interval for energy deposition
      real(REAL8), allocatable :: diffusivity(:)   ! radiation induced diffusivity [cm2/dpa]
      real(REAL8) :: dnear                   ! maximum distance to be considered near to
                                             !    the surface
      real(REAL8) :: dparelax                ! critical number of displacements per atom
                                             !   accumulated between cell relaxations so that
                                             !   relaxation is possible
      real(REAL8) :: drdeflow                ! flow coefficient [A/eV]
      real(REAL8), allocatable :: ed(:,:)    ! (nelem,nelem) displacement energy [eV]
      real(REAL8) :: elimkp                  ! energy limit for Kinchin-Pease model [eV]
      real(REAL8) :: emelt                   ! heat of melting
      real(REAL8) :: eminflow                ! minimum energy for flow [eV]
      real(REAL8) :: epockmin                ! minimum energy of amorphous pocket
      real(REAL8), allocatable :: esurf(:,:) ! (nelem,nelem) surface binding energy [eV]
      real(REAL8) :: esurfmin                ! minval(esurf(:,:))
      real(REAL8) :: facflow                 ! factor determining additional displacement at
                                             !    the end of trajectories [A/eV]
      real(REAL8) :: facsat                  ! quantity derived from damsat and dens
      real(REAL8) :: fpock                   ! factor for pocket recrystalliation
      real(REAL8), allocatable :: fracdest(:)! fraction of diffusivity due to destination
      real(REAL8) :: fracdi                  ! fraction of di-interstitials
      real(REAL8) :: fracfour                ! fraction of four-interstitials
      real(REAL8) :: frachex                 ! fraction of hexagonal-interstitials
      real(REAL8) :: fracrg                  ! fraction of frec-1 considered during recoil
                                             !    generation
      real(REAL8) :: fracspl                 ! fraction of split-<110> interstitials
      real(REAL8) :: fractet                 ! fraction of tetragonal-interstitials
      real(REAL8) :: fractri                 ! fraction of tri-interstitial
      real(REAL8) :: frec                    ! recombination factor within cascade
      real(REAL8) :: gammasurf               ! surface energy density [J/m^2]
      real(REAL8) :: ksurf                   ! nonplanarity parameter of the surface potential
      real(REAL8) :: modbulk                 ! bulk modulus [Pa]
      real(REAL8) :: modshear                ! shear modulus [Pa]
      real(REAL8) :: npockmin                ! minimum number of pocket atoms after
                                             !    recrystallization
      real(REAL8) :: pfacflow                ! factor determining additional displacement at
                                             !    the end of trajectories [A/sqrt(eV)]
      real(REAL8) :: rcap                    ! capture radius
      real(REAL8) :: rcap2                   ! rcap**2
      real(REAL8) :: ribiia                  ! factor for IBIIA effect [A/(eV/A^3)]
      real(REAL8) :: viscosity               ! radiation induced viscosity [Pa*dpa]
      real(REAL8) :: zmaxkp                  ! maximum depth at which elimkp is applied

   end type Dam_type

   type(Dam_type), public :: dam             ! dam structure
!
!=======================================================================
!
contains
!
!=======================================================================
!
   subroutine readDam( )
!             =======
!
!     Read &damage record from the input file INP and initialize the dam structure
!
!  Arguments:
!
!     None
!
!  Local Variables:
!
      integer(INT4) :: atom1                    ! atom index of projectile
      integer(INT4) :: atom2                    ! atom index of target atom

      logical :: lamo2
      logical :: lamzon
      logical :: lcap
      logical :: ldam
      logical :: ldet
      logical :: lflow                          ! recoil flow flag
      logical :: libiia
      logical :: lpock
      logical :: lrcoil
      logical :: lrcoil1
      logical :: lrcoil2
      logical :: lrepl
      logical :: lstrmd
      logical :: reenter

      character(80) :: relaxmethod
      character(80) :: relaxmodel

      integer(INT4) :: ndamdim
      integer(INT4) :: ndim                     ! obsolete
      integer(INT4) :: nionupd
      integer(INT4) :: npockmin
      integer(INT4) :: wcol(1:3)

      real(REAL8) :: damamo
      real(REAL8) :: damsat
      real(REAL8) :: density
      real(REAL8) :: dibiia
      real(REAL8) :: diffusivity
      real(REAL8) :: dnear
      real(REAL8) :: dparelax
      real(REAL8) :: drdeflow                   ! flow coefficient
      real(REAL8) :: ed
      real(REAL8) :: elimkp
      real(REAL8) :: emelt
      real(REAL8) :: eminflow                   ! minimum energy for flow
      real(REAL8) :: esurf
      real(REAL8) :: facflow
      real(REAL8) :: fpock
      real(REAL8) :: fracdest
      real(REAL8) :: fracdi
      real(REAL8) :: fractet
      real(REAL8) :: fracfour
      real(REAL8) :: frachex
      real(REAL8) :: fracrg
      real(REAL8) :: fracspl
      real(REAL8) :: fractri
      real(REAL8) :: frec
      real(REAL8) :: gammasurf
      real(REAL8) :: ksurf
      real(REAL8) :: modbulk
      real(REAL8) :: modshear
      real(REAL8) :: pfacflow
      real(REAL8) :: psi110
      real(REAL8) :: rcap
      real(REAL8) :: ribiia
      real(REAL8) :: viscosity
      real(REAL8) :: zmaxkp

      namelist /damage/ atom1, atom2,  &
                        lamo2, lamzon, lcap, ldam, ldet, lflow, libiia, lpock, lrcoil,  &
                        lrcoil1, lrcoil2, lrepl, lstrmd, reenter,  &
                        relaxmethod, relaxmodel,  &
                        ndamdim, npockmin, wcol,  &
                        damamo, damsat, density, dibiia, diffusivity, dnear, dparelax,  &
                        drdeflow, ed, emelt, eminflow, esurf, facflow, fpock, fracdest,  &
                        fracdi, fracfour, frachex, fracrg, fracspl, fractet, fractri, frec,  &
                        gammasurf, ksurf, elimkp, modbulk, modshear, pfacflow, psi110, rcap,  &
                        ribiia, viscosity, zmaxkp

      namelist /damage/ nionupd, ndim           ! obsolete

      character(256) :: msg
      integer(INT4) :: ielem1, ielem1min, ielem1max
      integer(INT4) :: ielem2, ielem2min, ielem2max, ios
!
!  Start of Executable Code:
!  ------------------------
!
!  Allocate memory for displacement energies
!
      allocate( dam%diffusivity( nelem ) )
      allocate( dam%ed( nelem, nelem ) )
      allocate( dam%esurf( nelem, nelem ) )
      allocate( dam%fracdest( nelem ) )
!
!  Default values
!
      dam%lamo2 = .true.
      dam%lamzon = .false.
      dam%lcap = .false.
      dam%ldam = .not.ldyn
      dam%ldet = .false.
      dam%lflow = .false.
      dam%libiia = .false.
      dam%lpock = .false.
      dam%lrcoil = .false.
      dam%lrcoil1 = .true.
      dam%lrcoil2 = .true.
      dam%lrepl = .true.
      dam%lstrmd = .false.
      dam%reenter = ( .not.geom%lcart .or. ndim_setup > 1 )

      call undefine( dam%ndamdim )
      call undefine( dam%relaxmethod )
      call undefine( dam%relaxmodel )
      call undefine( dam%wcol(1:3) )

      call undefine( dam%cospsi110 )
      call undefine( dam%damamo )
      call undefine( dam%damsat )
      call undefine( dam%dibiia )
      call undefine( dam%dens )
      call undefine( dam%diffusivity(:) )
      call undefine( dam%dnear )
      call undefine( dam%dparelax )
      call undefine( dam%drdeflow )
      call undefine( dam%ed(:,:) )
      call undefine( dam%elimkp )
      call undefine( dam%emelt )
      call undefine( dam%eminflow )
      call undefine( dam%esurf )
      call undefine( dam%facflow )
      call undefine( dam%fpock )
      call undefine( dam%fracdi )
      call undefine( dam%fracdest(:) )
      call undefine( dam%fracfour )
      call undefine( dam%frachex )
      call undefine( dam%fracrg )
      call undefine( dam%fracspl )
      call undefine( dam%fractet )
      call undefine( dam%fractri )
      call undefine( dam%frec )
      call undefine( dam%gammasurf )
      call undefine( dam%ksurf )
      call undefine( dam%modbulk )
      call undefine( dam%modshear )
      call undefine( dam%npockmin )
      call undefine( dam%pfacflow )
      call undefine( dam%rcap )
      call undefine( dam%ribiia )
      call undefine( dam%viscosity )
      call undefine( dam%zmaxkp )
!
!  Loop over reading records
!  -------------------------
!
!  Open the input file
!
      call openFile( 11, 'INP', 'formatted', 'read' )
!          ========  i     i         i         i

      do
!
!  Set default values of namelist variables
!
         call undefine( atom1 )
         call undefine( atom2 )

         lamo2 = dam%lamo2
         lamzon = dam%lamzon
         lcap = dam%lcap
         ldam = dam%ldam
         ldet = dam%ldet
         lflow = dam%lflow
         libiia = dam%libiia
         lpock = dam%lpock
         lrcoil = dam%lrcoil
         lrcoil1 = dam%lrcoil1
         lrcoil2 = dam%lrcoil2
         lrepl = dam%lrepl
         reenter = dam%reenter

         call undefine( relaxmethod )
         call undefine( relaxmodel )

         call undefine( ndamdim )
         call undefine( ndim )
         call undefine( nionupd )
         call undefine( npockmin )
         call undefine( wcol(1:3) )

         call undefine( damamo )
         call undefine( damsat )
         call undefine( density )
         call undefine( dibiia )
         call undefine( diffusivity )
         call undefine( dnear )
         call undefine( dparelax )
         call undefine( drdeflow )
         call undefine( ed )
         call undefine( elimkp )
         call undefine( emelt )
         call undefine( eminflow )
         call undefine( esurf )
         call undefine( facflow )
         call undefine( fpock )
         call undefine( fracdest )
         call undefine( fracdi )
         call undefine( fracfour )
         call undefine( frachex )
         call undefine( fracrg )
         call undefine( fracspl )
         call undefine( fractet )
         call undefine( fractri )
         call undefine( frec )
         call undefine( gammasurf )
         call undefine( ksurf )
         call undefine( modbulk )
         call undefine( modshear )
         call undefine( pfacflow )
         call undefine( psi110 )
         call undefine( rcap )
         call undefine( ribiia )
         call undefine( viscosity )
         call undefine( zmaxkp )
!
!  Read first &damage record from the input file
!
         read( 11, nml=damage, iostat=ios, iomsg=msg )
         if( ios < 0 )  exit
         if( ios > 0 ) then
            print '(/a)',  trim( msg )
            goto 900
         end if
!
!  Check index parameters
!
         if( defined( atom1 ) ) then
            if( atom1 < 1 .or. atom1 > nelem )  goto 906
         end if

         if( defined( atom2 ) ) then
            if( atom2 < 1 .or. atom2 > nelem )  goto 906
         end if
!
!  Check on parameters inconsistent with index parameters
!
         if( defined( atom1 ) .or. defined( atom2 ) ) then       ! parameters without indices

            if( defined( relaxmethod ) )  goto 91093
            if( defined( relaxmodel ) )  goto 91093

            if( defined( ndamdim ) )  goto 9179
            if( defined( npockmin ) )  goto 918
            if( any( defined( wcol ) ) )  goto 91092

            if( defined( damamo ) )  goto 922
            if( defined( damsat ) )  goto 912
            if( defined( density ) )  goto 923
            if( defined( dibiia ) )  goto 9100
            if( defined( dnear ) )  goto 9104
            if( defined( dparelax ) )  goto 91094
            if( defined( drdeflow ) )  goto 9249
            if( defined( elimkp ) )  goto 910
            if( defined( emelt ) )  goto 917
            if( defined( eminflow ) )  goto 9251
            if( defined( facflow ) )  goto 914
            if( defined( fpock ) )  goto 919
            if( defined( fracdi ) )  goto 9106
            if( defined( fracfour ) )  goto 9108
            if( defined( frachex ) )  goto 9109
            if( defined( fracrg ) )  goto 915
            if( defined( fracspl ) )  goto 9105
            if( defined( fractet ) )  goto 91091
            if( defined( fractri ) )  goto 9107
            if( defined( frec ) )  goto 911
            if( defined( gammasurf ) )  goto 913
            if( defined( ksurf ) )  goto 921
            if( defined( modbulk ) )  goto 91096
            if( defined( modshear ) )  goto 91097
            if( defined( pfacflow ) )  goto 9141
            if( defined( psi110 ) )  goto 9102
            if( defined( ribiia ) )  goto 9101
            if( defined( rcap ) )  goto 916
            if( defined( viscosity ) )  goto 91095
            if( defined( zmaxkp ) )  goto 9103

         end if

         if( defined( atom1 ) ) then
            if( defined( diffusivity ) )  goto 91098
            if( defined( fracdest ) )  goto 91099
         end if
!
!  Define range of regions and atoms to which parameters apply
!
         if( defined( atom1 ) ) then
            ielem1min = atom1
            ielem1max = atom1
         else
            ielem1min = 1
            ielem1max = nelem
         endif
         ielem1 = ielem1min

         if( defined( atom2 ) ) then
            ielem2min = atom2
            ielem2max = atom2
         else
            ielem2min = 1
            ielem2max = nelem
         endif
         ielem2 = ielem2min
!
!  Copy parameters to dam structure
!
         dam%lamo2 = lamo2
         dam%lamzon = lamzon
         dam%lcap = lcap
         dam%ldam = ldam
         dam%ldet = ldet
         dam%lflow = lflow
         dam%libiia = libiia
         dam%lpock = lpock
         dam%lrcoil = lrcoil
         dam%lrcoil1 = lrcoil1
         dam%lrcoil2 = lrcoil2
         dam%lrepl = lrepl
         dam%lstrmd = lstrmd
         dam%reenter = reenter

         if( defined( ndamdim ) )  dam%ndamdim = ndamdim
         if( defined( ndim ) )  goto 981
         if( defined( nionupd ) )  goto 980
         if( any( defined( wcol(:) ) ) )  dam%wcol(:) = wcol(:)

         if( defined( relaxmethod ) ) then
            print *, '&damage: relaxmethod is currently ignored'
            if( to_lower( relaxmethod ) == 'euler' ) then
!               ========       i
               dam%relaxmethod = ENUM%relaxmethod%euler
            else if( to_lower( relaxmethod ) == 'lagrange' ) then
!                    ========       i
               dam%relaxmethod = ENUM%relaxmethod%lagrange
            else
               goto 9738
            end if
         end if

         if( defined( relaxmodel ) ) then
            if( to_lower( relaxmodel ) == 'none' ) then
!               ========      i
               dam%relaxmodel = ENUM%relaxmodel%none
            else if( to_lower( relaxmodel ) == 'full' ) then
!                    ========      i
               dam%relaxmodel = ENUM%relaxmodel%full
            else if( to_lower( relaxmodel ) == 'viscoelastic' ) then
!                    ========      i
               dam%relaxmodel = ENUM%relaxmodel%viscoelastic
            else
               goto 9739
            end if
         end if

         if( defined( damamo ) )  dam%damamo = damamo / 1e24   ! 1/cm3 -> 1/A^3
         if( defined( damsat ) )  dam%damsat = damsat / 1e24   ! 1/cm3 -> 1/A^3
         if( defined( density ) )  dam%dens = density / 1e24   ! 1/cm3 -> 1/A^3
         if( defined( dibiia ) )  dam%dibiia = dibiia
         if( defined( diffusivity ) )  &
               dam%diffusivity( ielem2min:ielem2max ) = diffusivity * 1e16 ! cm^2/dpa -> A^2/dpa
         if( defined( dnear ) )  dam%dnear = dnear
         if( defined( dparelax ) )  dam%dparelax = dparelax
         if( defined( drdeflow ) )  dam%drdeflow = drdeflow
         if( defined( ed ) )  dam%ed(ielem2min:ielem2max,ielem1min:ielem1max) = ed
         if( defined( elimkp ) )  dam%elimkp = elimkp
         if( defined( emelt ) )  dam%emelt = emelt
         if( defined( eminflow ) )  dam%eminflow = eminflow
         if( defined( esurf ) )  dam%esurf( ielem2min:ielem2max, ielem1min:ielem1max ) = esurf
         if( defined( facflow ) )  dam%facflow = facflow
         if( defined( fpock ) )  dam%fpock = fpock
         if( defined( fracdest ) )  dam%fracdest( ielem2min:ielem2max ) = fracdest
         if( defined( fracdi ) )  dam%fracdi = fracdi
         if( defined( fracfour ) )  dam%fracfour = fracfour
         if( defined( frachex ) )  dam%frachex = frachex
         if( defined( fracrg ) )  dam%fracrg = fracrg
         if( defined( fracspl ) )  dam%fracspl = fracspl
         if( defined( fractet ) )  dam%fractet = fractet
         if( defined( fractri ) )  dam%fractri = fractri
         if( defined( frec ) )  dam%frec = frec
         if( defined( gammasurf ) )  dam%gammasurf = gammasurf
         if( defined( ksurf ) )  dam%ksurf = ksurf
         if( defined( modbulk ) )  dam%modbulk = modbulk
         if( defined( modshear ) )  dam%modshear = modshear
         if( defined( npockmin ) )  dam%npockmin = npockmin
         if( defined( pfacflow ) )  dam%pfacflow = pfacflow
         if( defined( psi110 ) )  dam%cospsi110 = cos( psi110*PI/180 )
         if( defined( ribiia ) )  dam%ribiia = ribiia
         if( defined( rcap ) )  dam%rcap = rcap
         if( defined( viscosity ) )  dam%viscosity = viscosity
         if( defined( zmaxkp ) )  dam%zmaxkp = zmaxkp
!
!  Checks that cannot be done in initDam
!
         if( defined( psi110 ) .and. psi110 < ZERO .or. psi110 > 180 )  goto 966

      end do
!
!  End of loop over reading records
!  --------------------------------
!
      close( 11 )
!
!  Set parameters in other records
!
      if( .not.ldamdyn_set )  ldamdyn = dam%ldam

      return
!
!  Error messages
!
900   error stop '&damage: error during reading namelist'

906   error stop '&damage: 1<=atom1,atom2<=natom required'

910   error stop '&damage: elimkp and atom1/atom2 must not be on same record'
911   error stop '&damage: frec and atom1/atom2 must not be on same record'
912   error stop '&damage: damsat and atom1/atom2 must not be on same record'
913   error stop '&damage: gammasurf and atom1/atom2 must not be on same record'
914   error stop '&damage: facflow and atom1/atom2 must not be on same record'
9141  error stop '&damage: pfacflow and atom1/atom2 must not be on same record'
915   error stop '&damage: fracrg and atom1/atom2 must not be on same record'
916   error stop '&damage: rcap and atom1/atom2 must not be on same record'
917   error stop '&damage: emelt and atom1/atom2 must not be on same record'
9179  error stop '&damage: ndamdim and atom1/atom2 must not be on same record'
918   error stop '&damage: npockmin and atom1/atom2 must not be on same record'
919   error stop '&damage: fpock and atom1/atom2 must not be on same record'
9100  error stop '&damage: dibiia and atom1/atom2 must not be on same record'
9101  error stop '&damage: ribiia and atom1/atom2 must not be on same record'
9102  error stop '&damage: psi110 and atom1/atom2 must not be on same record'
9103  error stop '&damage: zmaxkp and atom1/atom2 must not be on same record'
9104  error stop '&damage: dnear and atom1/atom2 must not be on same record'
9105  error stop '&damage: fracspl and atom1/atom2 must not be on same record'
9106  error stop '&damage: fracdi and atom1/atom2 must not be on same record'
9107  error stop '&damage: fractri and atom1/atom2 must not be on same record'
9108  error stop '&damage: fracfour and atom1/atom2 must not be on same record'
9109  error stop '&damage: frachex and atom1/atom2 must not be on same record'
91091 error stop '&damage: fractet and atom1/atom2 must not be on same record'
91092 error stop '&damage: wcol and atom1/atom2 must not be on same record'
91093 error stop '&damage: relaxmodel or relaxmethod and atom1/atom2 must not be on same record'
91094 error stop '&damage: dparelax and atom1/atom2 must not be on same record'
91095 error stop '&damage: viscosity and atom1/atom2 must not be on same record'
91096 error stop '&damage: modbulk and atom1/atom2 must not be on same record'
91097 error stop '&damage: modshear and atom1/atom2 must not be on same record'
91098 error stop '&damage: diffusivity and atom1 must not be on same record'
91099 error stop '&damage: fracdest and atom1 must not be on same record'
!
921   error stop '&damage: ksurf and atom1/atom2 must not be on same record'
922   error stop '&damage: damamo and atom1/atom2 must not be on same record'
923   error stop '&damage: density and atom1/atom2 must not be on same record'
9249  error stop '&snpar: atom1 or atom2 and drdeflow must not be on same record'
9251  error stop '&damage: atom1 or atom2 and eminflow must not be on same record'

966   error stop '&damage: 0<=psi110<=180 required'

9738  error stop '&damage: illegal value of relaxmethod'
9739  error stop '&damage: illegal value of relaxmodel'

980   error stop '&damage: nionupd is obsolete. Use nionhis of &setup instead.'
981   error stop '&damage: ndim is obsolete. Use ndamdim instead.'

   end subroutine readDam
!
!=======================================================================
!
   subroutine initDam( )
!             =======
!
!     Initialize the dam structure
!
!  Arguments:
!
!     None
!
!  Local Variables:
!
      integer(INT4) :: ielem1, ielem2, ir
!
!  Start of Executable Code:
!  ------------------------
!
!  Checks that have to be done before assigning default values
!
      if( dam%ldet .and. defined( dam%elimkp ) )  goto 9611
!
!  Set default values
!
      if( .not.defined( dam%ndamdim ) ) then
         if( ion%doseunits == 'cm-2' ) then
            if( ion%dxinit == ZERO ) then
               dam%ndamdim = 1
            else if( ion%dyinit == ZERO ) then
               dam%ndamdim = 2
            else
               dam%ndamdim = 3
            end if
         else if( ion%doseunits == 'cm-1' ) then
            if( ion%dyinit == ZERO ) then
               dam%ndamdim = 2
            else
               dam%ndamdim = 3
            end if
         else
            dam%ndamdim = 3
         end if
      end if

      if( .not.defined( dam%wcol(1) ) )  dam%wcol(1) = 40
      if( .not.defined( dam%wcol(2) ) )  dam%wcol(2) = 40
      if( .not.defined( dam%wcol(3) ) )  dam%wcol(3) = 8000

      ! default value of dam%damamo see after dam%dens
      if( .not.defined( dam%cospsi110 ) )  dam%cospsi110 = ONE

      if( .not.defined( dam%dens ) ) then
         if( irxtal > 0 ) then
            dam%dens = mat(irxtal)%dens
         else if( .not.usecell .and. nr == 1 ) then
            dam%dens = mat(1)%dens
         else
            dam%dens = 4.994e-2     ! = 4.994e22 cm-3 !
         end if
      end if

      if( .not.defined( dam%damamo ) ) then
         if( usecell ) then
            dam%damamo = dam%dens
         else if( any( mat(:)%name == 'si' ) .or. any( mat(:)%name == 'Si' )  &
                 .or. any( mat(:)%name == 'SI' ) ) then
            dam%damamo = 0.027            ! = 2.7e22 cm-3
         else
            dam%damamo = dam%dens
         end if
      end if

      if( .not.defined( dam%damsat ) )  dam%damsat = BIG / 1e24
      if( .not.defined( dam%dibiia ) )  dam%dibiia = 10
      if( .not.defined( dam%dnear ) )  dam%dnear = 30.
      if( .not.defined( dam%dparelax ) )  dam%dparelax = 1e-3
      if( .not.defined( dam%drdeflow ) )  dam%drdeflow = 0.2

      do ielem2 = 1, nelem
         if( .not.defined( dam%diffusivity(ielem2) ) )  dam%diffusivity(ielem2) = ZERO
      end do

      do ielem1 = 1, nelem
         do ielem2 = 1, nelem
            if( .not.defined( dam%ed(ielem2,ielem1) ) )  dam%ed(ielem2,ielem1) = 15.
         end do
      end do

      do ielem1 = 1, nelem
         if( .not.defined( dam%esurf(ielem1,ielem1) ) )  &
               dam%esurf(ielem1,ielem1) = elem(ielem1)%esurf
      end do
      do ielem1 = 1, nelem
         do ielem2 = 1, nelem
            if( ielem2 == ielem1 )  cycle
            if( .not. defined( dam%esurf(ielem2,ielem1) ) )  &
                  dam%esurf(ielem2,ielem1) =  &
                        HALF * ( dam%esurf(ielem1,ielem1) + dam%esurf(ielem2,ielem2) )
         end do
      end do

      do ielem2 = 1, nelem
         if( .not.defined( dam%fracdest(ielem2) ) )  dam%fracdest(ielem2) = ZERO
      end do

      if( .not.defined( dam%elimkp ) )  dam%elimkp = ZERO
      if( .not.defined( dam%emelt ) )  dam%emelt = 0.52
      if( .not.defined( dam%eminflow ) )  dam%eminflow = 5.
      if( .not.defined( dam%facflow ) )  dam%facflow = ZERO
      if( .not.defined( dam%fpock ) )  dam%fpock = ONE
      if( .not.defined( dam%fracdi ) )  dam%fracdi = ZERO
      if( .not.defined( dam%fracfour ) )  dam%fracfour = ZERO
      if( .not.defined( dam%frachex ) )  dam%frachex = ZERO
      if( .not.defined( dam%fracrg ) )  dam%fracrg = ONE
      if( .not.defined( dam%fracspl ) )  dam%fracspl = ONE
      if( .not.defined( dam%fractet ) )  dam%fractet = ZERO
      if( .not.defined( dam%fractri ) )  dam%fractri = ZERO

      if( .not.defined( dam%frec ) ) then
         select case( ion%name )
            case( 'B', 'b' )
               dam%frec = 0.3          ! used to be: frec = 0.125, damsat = 0.004 = 4.e21 cm-3,
            case( 'BF2', 'bf2' )       !             damamo = ZERO
               dam%frec = 1.2
            case( 'SI', 'Si', 'si' )
               dam%frec = 0.8
            case( 'P', 'p' )
               dam%frec = 0.8
            case( 'AS', 'As', 'as' )
               dam%frec = 2.
            case default
               dam%frec = 0.0264*ion%m
         end select
      end if

      if( .not.defined( dam%gammasurf ) )  dam%gammasurf = 1.13      ! value for Si
      if( .not.defined( dam%ksurf ) )  dam%ksurf = ZERO
      if( .not.defined( dam%modbulk ) )  dam%modbulk = 4.8e10        ! value for Si
      if( .not.defined( dam%modshear ) )  dam%modshear = 3.3e10      ! value for Si
      if( .not.defined( dam%npockmin ) )  dam%npockmin = 10.
      if( .not.defined( dam%pfacflow ) )  dam%pfacflow = ZERO
      if( .not.defined( dam%rcap ) )  dam%rcap = (3/(4*PI*2*dam%damsat))**THIRD
      if( .not.defined( dam%relaxmethod ) )  dam%relaxmethod = ENUM%relaxmethod%euler
      if( .not.defined( dam%relaxmodel ) )  dam%relaxmodel = ENUM%relaxmodel%full
      if( .not.defined( dam%ribiia ) )  dam%ribiia = 100
      if( .not.defined( dam%viscosity ) )  dam%viscosity = 0.966e9   ! value for Si
      if( .not.defined( dam%zmaxkp ) )  dam%zmaxkp = BIG
!
!  Check validity of the parameters
!
      if( ldamdyn .and. .not.(dam%ldam .or. dam%lrcoil) )  goto 930
      if( nthread > 1 .and. dam%ldet )  goto 940
      if( dam%lcap .and. dam%lpock )  goto 992
      if( dam%lcap .and. .not.dam%lrcoil )  goto 993
      if( dam%lpock .and. .not.dam%lrcoil )  goto 994
      if( dam%lamzon .and. usecell )  goto 9501
      if( dam%lamzon .and. irxtal == 0 )  goto 951
      if( dam%lamzon .and. dam%ndamdim /= 1 )  goto 952
      if( dam%ldet .and. irxtal == 0 )  goto 955
      if( dam%ldet .and. .not.( usehis .or. usebk .or. usedef ) )  goto 9725
      if( dam%ldet .and. usebk .and. useion )  goto 9726
      if( dam%ldet .and. (dam%lcap .or. dam%lpock) )  goto 9727
      if( dam%ldet .and. dam%ldam .and. usedef .and. .not.dam%lrcoil )  goto 9728
      if( dam%ldet .and. xtal%uc%name /= ' ' )  goto 9729
      if( dam%libiia .and. (.not. dam%ldam .or. dam%ndamdim /= 1) )  goto 954

      if( dam%ndamdim < 1 .or. dam%ndamdim > 3 )  goto 9739
      if( ion%doseunits == 'cm-2' ) then
         if( ion%dxinit == ZERO ) then
            if( dam%ndamdim > 1 )  goto 9744
         else if( ion%dyinit == ZERO ) then
            if( dam%ndamdim > 2 )  goto 9744
         end if
      else if( ion%doseunits == 'cm-1' ) then
         if( ion%dyinit == ZERO ) then
            if( dam%ndamdim > 2 )  goto 9744
         end if
      end if

#ifndef AMAT
!      if( dam%relaxmethod == ENUM%relaxmethod%euler )  goto 9746
#endif
#ifndef __INTELCOMPILER
      if( ndim_setup == 2 .and. dam%relaxmethod == ENUM%relaxmethod%lagrange  &
          .and. dam%relaxmodel == ENUM%relaxmodel%viscoelastic )  goto 9745
#endif
      if( dam%relaxmodel == 3 .and. .not. geom%grid == ENUM%grid%regular )  goto 967

      if( .not. all( defined( dam%wcol(:) ) ) )  goto 9737
      if( all( modulo( dam%wcol(:), 4 ) /= 0 ) )  goto 9738

      if( dam%damamo <= ZERO )  goto 964
      if( dam%damsat <= ZERO )  goto 963
      if( dam%dens <= ZERO )  goto 965
      if( any( dam%diffusivity(:) < ZERO ) )  goto 975
      if( dam%dnear < ZERO )  goto 966
      if( dam%dibiia <= ZERO )  goto 984
      if( dam%dparelax < ZERO ) goto 9740
      if( dam%drdeflow < ZERO )  goto 959
      if( any( dam%ed(:,:) <= ZERO ) )  goto 960
      if( dam%elimkp < ZERO )  goto 961
      if( dam%emelt <= ZERO )  goto 981
      do ielem1 = 1, nelem
         do ielem2 = 1, nelem
            if( ( ldyn .or. ( dam%lrcoil .and. ielem1 > nelemIon ) ) .and.  &
                lsput .and. dam%esurf(ielem2,ielem1) < 0.1 )  &
                  print '(a,i2,a,i2,a,es10.3e2,a)',  &
                        ' (WARNING) &damage: esurf(', ielem2, ',', ielem1, ') =',  &
                        dam%esurf(ielem2,ielem1), ' will result in very slow simulation'
         end do
      end do
      if( any( dam%esurf(:,:) < ZERO ) )  goto 9601
      if( lsput .and. any( dam%esurf(:,:) <= ZERO ) )  goto 9602
      if( dam%facflow < ZERO )  goto 986
      if( dam%fpock < ZERO )  goto 982
      if( any( dam%fracdest(:) < ZERO ) .or. any( dam%fracdest(:) > ONE ) )  goto 976
      if( dam%fracdi < ZERO .or. dam%fracdi > ONE )  goto 9731
      if( dam%fracfour < ZERO .or. dam%fracfour > ONE )  goto 9733
      if( dam%frachex < ZERO .or. dam%frachex > ONE )  goto 9734
      if( dam%fracrg < ZERO .or. dam%fracrg > ONE )  goto 972
      if( dam%fractet < ZERO .or. dam%fractet > ONE )  goto 9735
      if( dam%fracspl < ZERO .or. dam%fracspl > ONE )  goto 9730
      if( dam%fractri < ZERO .or. dam%fractri > ONE )  goto 9732
      if( dam%ldet .and. abs( dam%fracspl + dam%fracdi + dam%fractri + dam%fracfour +  &
                              dam%frachex + dam%fractet - ONE ) > EPSMB )  goto 9736
      if( dam%frec <= 0 )  goto 962
      if( dam%gammasurf < 0 )  goto 968
      if( dam%modbulk <= ZERO )  goto 9742
      if( dam%modshear <= ZERO )  goto 9743
      if( dam%npockmin < 0 )  goto 983
      if( dam%facflow < ZERO )  goto 9861
      if( dam%rcap <= ZERO )  goto 980
      if( dam%ribiia <= ZERO )  goto 985
      if( dam%viscosity <= ZERO )  goto 9741
!
!  Calculate derived quantities
!
      dam%lcasc = (dam%lcap .or. dam%lpock) .and. dam%ldam
                                                ! epockmin
      dam%epockmin = dam%emelt * (dam%npockmin**THIRD + dam%fpock)**3
                                                ! esurfmin(ielem1)
      if( ldyn ) then
         dam%esurfmin = minval(dam%esurf(:,:))
      else
         dam%esurfmin = minval(dam%esurf(nelemIon+1:,nelemIon+1:))
      end if
                                                ! facsat
      dam%facsat = ( ONE - ONE / ( ONE + dam%frec - dam%damsat / dam%dens *  &
            (ONE-dam%fracrg + dam%fracrg*dam%frec) ) ) / dam%damsat
      if( dam%facsat < ZERO )  goto 990
                                                ! rcap2
      dam%rcap2 = dam%rcap**2
                                                ! mat%ed and mat%esurf
      if( .not.usecell ) then
         do ir = 1, nr
            allocate( mat(ir)%ed( nelem ) )
            allocate( mat(ir)%esurf( nelem ) )
            do ielem1 = 1, nelem
               mat(ir)%ed(ielem1) = sum( mat(ir)%prob(1:mat(ir)%nelemr) *  &
                                         dam%ed(mat(ir)%ielem(1:mat(ir)%nelemr),ielem1) )
               mat(ir)%esurf(ielem1) = sum( mat(ir)%prob(1:mat(ir)%nelemr) *  &
                                            dam%esurf(mat(ir)%ielem(1:mat(ir)%nelemr),ielem1) )
            end do
         end do
      end if
!
!  Initialize ion%ldet
!
      ion%ldet = dam%ldet

      return

930   error stop '&damage: ldamdyn=t requires ldam=t or lrcoil=t'

940   error stop '&damage: ldet=t not allowed for parallel execution'

9501  error stop '&damage: lamzon=t requires usecell=f'
951   error stop '&damage: lamzon=t requires a crystalline region'
952   error stop '&damage; lamzon=t requires 1-d damage histogram (dam%ndamdim=1)'
954   error stop '&damage: libiia=t requires ldam=t and 1-d damage histogram (dam%ndamdim=1)'
955   error stop '&damage: ldet=t requires crystalline region'

959   error stop '&damage: drdeflow>=0 required'
960   error stop '&damage: ed>0 required'
9601  error stop '&damage: esurf>=0 required'
9602  error stop '&damage: esurf>0 required'
961   error stop '&damage: elimkp>=0 required'
9611  error stop '&damage: specifying elimkp with ldet=t is not allowed'
962   error stop '&damage: frec>0 required'
963   error stop '&damage: damsat>=0 required'
964   error stop '&damage: damamo>0 required'
965   error stop '&damage: density>=0 required'
966   error stop '&damage: dnear>=0 required'
967   error stop '&damage: relaxmodel="Euler" requires grid=''regular'' in &geom'
968   error stop '&damage: gammasurf>=0 required'

972   error stop '&damage: 0<=fracrg<=1 required'
9725  error stop '&damage: ldet=t requires one of usehis, usebk, usedef in &setup'
9726  error stop '&damage: specifying both useion=t and usebk=t with ldet=t not allowed'
9727  error stop '&damage: specifying lcap=t or lpock=t with ldet=t not allowed'
9728  error stop '&damage: ldet=t and ldam=t and usedef=t require lrcoil=t'
9729  error stop '&damage: ldet=t requires default diamond lattice (name='' '')'
9730  error stop '&damage: 0<=fracspl<=1 required'
9731  error stop '&damage: 0<=fracdi<=1 required'
9732  error stop '&damage: 0<=fractri<=1 required'
9733  error stop '&damage: 0<=fracfour<=1 required'
9734  error stop '&damage: 0<=frachex<=1 required'
9735  error stop '&damage: 0<=fractet<=1 required'
9736  error stop '&damage: fracspl+fracdi+fractri+fracfour+fractet+frachex=1 required'
9737  error stop '&damage: all or no components of wcol must be specified'
9738  error stop '&damage: all components of wcol must be full multiples of 4'
9739  error stop '&damage: 1<=ndamdim<=3 required'
9740  error stop '&damage: dparelax>=0 required'
9741  error stop '&damage: viscosity>0 required'
9742  error stop '&damage: modbulk>0 required'
9743  error stop '&damage: modshear>0 required'
9744  error stop '&damage: ndamdim too large because damage concentration cannot be calculated'
9745  error stop '&damage: ndim=2 and relaxmodel=''Viscoelastic'' requires Intel compiler'
!9746  error stop '&damage: illegal value of relaxmethod'
975   error stop '&damage: diffusivity>=0 required'
976   error stop '&damage: 0<=fracdest<=1 required'

980   error stop '&damage: rcap>0 required'
981   error stop '&damage: emelt>0 required'
982   error stop '&damage: fpock>=0 required'
983   error stop '&damage: npockmin>=0 required'
984   error stop '&damage: dibiia>0 required'
985   error stop '&damage: ribiia>0 required'
986   error stop '&damage: facflow>=0 required'
9861  error stop '&damage: pfacflow>=0 required'

990   error stop '&damage: damsat too big'
991   error stop '&damage: lcap=t requires specification of rcap or damsat'
992   error stop '&damage: specifying both lpock=t and lcap=t not allowed'
993   error stop '&damage: lcap=t requires lrcoil=t'
994   error stop '&damage: lpock=t requires lrcoil=t'

   end subroutine initDam
!
!=======================================================================
!
   subroutine writeDam( )
!             ========
!
!     Write parameters derived from the &damage record to the output file
!
!  Arguments:
!
!     None
!
!  Local Variables:
!
      real(REAL8) :: rmax
!
!  Start of Executable Code:
!  ------------------------
!
!  Open the output file
!
      call openFile( 21, 'OUT', 'formatted', 'append' )
!          ========  i     i         i          i

      write( 21, '(t2,a)', advance='no' ) '&damage:   '
      if( dam%lrcoil ) then
         write( 21, '(a)' )  &
               'lrcoil=t: Recoils are followed'
         if( ndim_setup > 1 .or. .not.geom%lcart ) then
            if( .not.dam%lrcoil1 ) then
               write( 21, '(t16,a)' )  'except before reentering the target (lrcoil1=f)'
            end if
            if( .not.dam%lrcoil2 ) then
               write( 21, '(t16,a)' )  'except after reentering the target (lrcoil2=f)'
            end if
            if( dam%reenter ) then
               write( 21, '(t13,a)' )  &
                     'reenter=t: Recoils are followed through vacuum regions'
            else
               write( 21, '(t13,a)' )  &
                     'reenter=f: Recoils are not followed through vacuum regions'
            end if
         end if
         if( dam%elimkp > ZERO ) then
            write( 21, '(t13,a)', advance='no' )  &
                  'Recoils are followed above elimkp=' // str( nint(dam%elimkp) )
            if( irxtal /= 0 ) then
               write( 21, '(a)' ) ', if they are further from'
               write( 21, '(t13,a)' )  &
                     '   <110> than' // str( 180/PI*acos( dam%cospsi110 ), 'f99.1' ) //  &
                     ' degrees and shallower than z=' // str( dam%zmaxkp, 'f99.1' )
            else
               write( 21, '()' )
            end if
         else
            write( 21, '(t13,a)' )  &
                  'elimkp=0: Recoils are followed until the end of their trajectories'
         end if
      else
         write( 21, '(a)' )  &
               'lrcoil=f: Recoils are not followed'
      end if

      if( dam%ldam ) then
         write( 21, '(t13,a)' )  &
               'ldam=t: Refined damage model'
         if( dam%ldet ) then
            write( 21, '(t13,a)' )  &
                  'ldet=t: Damage is treated deterministically'
            write( 21, '(t13,a)' )  &
                  'Column dimensions wcol=(' // str( dam%wcol(1) ) //  &
                  ',' // str( dam%wcol(2) ) // ')'
            if( dam%lstrmd ) then
               write( 21, '(t13, a)' )  &
                     'lstrmd=t: Strain from MD calculations is used'
            else
               write( 21, '(t13, a)' )  &
                     'lstrmd=f: Strain from ab-initio calculations is used'
            end if
         else
            write( 21, '(t13, a)' )  &
                  'ldet=f: Damage is treated statistically based on '  &
                  // str( dam%ndamdim, 'i1' ) // '-D densities'
         end if
      else
         if( dam%lrcoil ) then
            write( 21, '(t13,a)' )  &
                  'ldam=f: No refined damage model'
         else
            write( 21, '(t13,a)' )  &
                  'ldam=f: Damage is not calculated'
         end if
      end if

      if( irxtal > 0 .and. ldamdyn ) then
         write( 21, '(t13,a)' )  &
               'Damage is considered for subsequent trajectories'
         write( 21, '(t13,a)' )  &
               'ndamdim=' // str( dam%ndamdim, 'i1' ) // ': using ' //  &
               str( dam%ndamdim, 'i1' ) // '-D densities'
      end if

      if( dam%lflow ) then
         if( dam%facflow > ZERO ) then
            write( 21, '(t13,a)' )  &
                  'Flow coefficient facflow=' // str( dam%facflow, 'f99.3' ) //  &
                  ' relating flow displacement'
            write( 21, '(t13,a)' )  &
                   '   to energy deposition within capture radius rcap=' //  &
                   str( dam%rcap, 'f99.3' )
         end if
         if( dam%pfacflow > ZERO ) then
            write( 21, '(t13,a)' )  &
                  'Flow coefficient pfacflow=' // str( dam%pfacflow, 'f99.3' ) //  &
                  ' relating flow displacement'
            write( 21, '(t13,a)' )  &
                   '   to momentum deposition within capture radius rcap=' //  &
                   str( dam%rcap, 'f99.3' )
         end if
      end if

      if( dam%lcap ) then
         write( 21, '(t13,a)' )  &
               'lcap=t: Recombination within one recoil cascade is taken into'
         rmax = 9999.999
         write( 21, '(t13,a)' )  &
               '   account by the capture radius model: rcap=' // str( dam%rcap, 'f99.3' )
      end if

      if( dam%lpock ) then
         write( 21, '(/t13,a)' ) 'lpock=t: Amorphous pocket model: '
         write( 21, '(t13,a)' )  &
               'emelt=' // str( dam%emelt, 'f99.3' ) //  &
               ', fpock=' // str( dam%fpock, 'f99.3' ) //  &
               ', npockmin=' // str( dam%npockmin, 'f99.3' ) //  &
               ', epockmin=' // str( dam%epockmin )
         write( 21, '()' )
      end if

      write( 21, '(/t13,a,f6.3)', advance='no' )  &
            'frec=' // str( dam%frec, 'f99.3' ) //  &
            ', fracrg=' // str( dam%fracrg, 'f99.4' )
      if( defined( dam%damsat ) )  &
         write( 21, '(a)', advance='no' )  &
               ', damsat=' // str( 1e24*dam%damsat, 'es8.2e2' )
      write( 21, '(a)' )  &
            ', damamo=' // str( 1e24*dam%damamo, 'es8.2e2' )

      write( 21, '(t13,a)' )  &
            'density=' // str( 1e24*dam%dens, 'es8.2e2' )

      write( 21, '(t13,a)' )  &
            'lrepl=' // str( dam%lrepl ) //  &
            ', lamzon=' // str( dam%lamzon ) //  &
            ', lamo2=' // str( dam%lamo2 )

      if( dam%libiia ) then
         write( 21, '(/t13,a)' )  &
               'libiia=t: ion beam induced amorphization model'
         write( 21, '(t13,a)' )  &
               'ribiia=' // str( dam%ribiia, 'es10.3' ) //  &
               ', dibiia=' // str( dam%dibiia, 'es10.3' )
      end if

      if( dam%ldet ) then
         write( 21, '(t13,a,6(a,0pf6.4))' )  &
               'Fractions of defects: fracspl=' // str( dam%fracspl, 'f99.4' ) //  &
               ', fracdi=' // str( dam%fracdi, 'f99.4' ) //  &
               ', fractri=' // str( dam%fractri, 'f99.4' ) //  &
               ', fracfour=' // str( dam%fracfour, 'f99.4' ) //  &
               ', frachex=' // str( dam%frachex, 'f99.4' ) //  &
               ', fractet=' // str( dam%fractet, 'f99.4' )
      end if
!
!  Displacement energy
!
      write( 21, '()' )
      call writeTab( dam%ed(:,:), 'ed' )
!          ========            i              i
!
!  Surface binding energy
!
      if( lsput ) then
         call writeTab( dam%esurf(:,:), 'esurf' )
!             ========            i                   i
         write( 21, '(t13,a)' )  &
               'ksurf=' // str( dam%ksurf, 'f99.3' )

         write( 21, '(/t13,a)' )  &
               'Virtual recoils are generated within dnear=' // str( dam%dnear, 'f99.1' ) //  &
               ' from the surface'
         write( 21, '()' )
      end if

      if( ldyn ) then
         if( any( dam%diffusivity(:) /= ZERO ) ) then
            call writeTab( dam%diffusivity(:), 'diffusivity' )
!                ========          i                 i
            if( any( dam%fracdest(:) /= ZERO ) ) then
               call writeTab( dam%fracdest(:), 'diffusivity' )
!                   ========        i                i
            else
               write( 21, '(t13,a)' )  &
                     'fracdest=0'
            end if
         else
            write( 21, '(/t13,a)' )  &
                  'diffusivity=0'
         end if
         if( dam%relaxmodel == ENUM%relaxmodel%full ) then
            write( 21, '(/t13,a)' )  &
                  'Relax model: full'
            if( ndim_setup == 2 ) then
               write( 21, '(t13,a)' )  &
                     '             Cell relaxation where displacements per atom exceed ' //  &
                     str( dam%dparelax, 'es10.3' )
            end if
         else if( dam%relaxmodel == ENUM%relaxmodel%viscoelastic ) then
            write( 21, '(/t13,a)' )  'Relax model: viscoelastic'
            write( 21, '(t16,a)' )  &
                  'viscosity=' // str( dam%viscosity, 'es10.3' ) //  &
                  ', surface energy=' // str( dam%gammasurf, 'es10.3' )
            write( 21, '(t16,a)' )  &
                  'bulk modulus=' // str( dam%modbulk, 'es10.3' ) //  &
                  ', shear modulus=' // str( dam%modshear, 'es10.3' )
         end if
         if( dam%relaxmodel /= ENUM%relaxmodel%none ) then
            !if( dam%relaxmethod == ENUM%relaxmethod%euler ) then
            !   write( 21, '(t16,a)' )  &
            !         'Relax method: Eulerian'
            !else if( dam%relaxmodel == ENUM%relaxmethod%lagrange ) then
            !   write( 21, '(t16,a)' )  &
            !         'Relax method: Lagrangian'
            !end if
            continue
         end if
      end if
      write( 21, '()' )

      if( dam%lflow ) then
         write( 21, '(t13,a)' )  &
               'lflow=t: additional recoil flow at end of trajectories'
         write( 21, '(t13,a)' )  &
               '   eminflow=' // str( dam%eminflow, 'f99.3' ) //  &
               ', drdeflow=' // str( dam%drdeflow, 'f99.3' )
!      else
!         write( 21, '(t13,a)' )  &
!               'lflow=f: no additional recoil flow at end of trajectories'
      end if

      close( 21 )

   end subroutine writeDam
!
!=======================================================================
!
end module Param_Dam_mod
