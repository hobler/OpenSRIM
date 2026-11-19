module Param_Out_mod
!
!  Read &output record from the input file INP and calculate derived quantities
!
!  Used Modules:
!  ============
!
   use Const_mod,       only: INT4, REAL8, &
                              ZERO, ONE, HALF, PI, BIG,  &
                              ENUM
   use String_mod,      only: str
   use Undefine_mod,    only: undefine, defined, definedlength
   use Open_mod,        only: openFile

   use Param_Setup_mod, only: nr, nthread,  &
                              nion, niongeom, nionout, niontra_setup => niontra,  &
                              usehis, ldyn, ndim, nion, usecell
   use Param_Geom_mod,  only: geom
   use Param_Elem_mod,  only: elem
   use Param_Ion_mod,   only: ion
   use Param_Mat_mod,   only: mat, irxtal
   use Param_Sn_mod,    only: sn
   use Param_Dam_mod,   only: dam
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
   public :: readOut                ! read the output parameters
!            =======
   public :: initOut                ! initialize the out structure
!            =======
   public :: writeOut               ! write output parameters to output file
!            ========
!
!  Public Data:
!  ===========
!
   type :: Out_type
!          ========
      logical :: lcell                    ! dynamic cell target output flag
      logical :: lcyl                     ! cylinder coordinate output flag
      logical :: ldef                     ! deterministic defects output flag
      logical :: ldist                    ! distance function output flag
      logical :: lextendx(2)              ! flag indicating that the score domain extends
      logical :: lextendy(2)              !    beyond the geometry domain in that direction
      logical :: lextendz(2)              !
      logical :: lhis                     ! 1D spatial distribution histogram flag
      logical :: lhis2                    ! 2D spatial distribution histogram flag
      logical :: lhis2b                   ! backscattered atoms 2D spatial distribution flag
      logical :: lhis2t                   ! transmitted atoms 2D spatial distribution flag
      logical :: lhis3                    ! 3D spatial distribution histogram flag
      logical :: lhisa2b                  ! 2D backscattered atoms angular distribution flag
      logical :: lhisa2t                  ! 2D transmitted atoms angular distribution flag
      logical :: lhisaab                  ! backscattered atoms 3D angular distribution flag
      logical :: lhisaat                  ! transmitted atoms 3D angular distribution flag
      logical :: lhisab                   ! backscattered atoms 2D angular distribution flag
      logical :: lhisat                   ! transmitted atoms 2D angular distribution flag
      logical :: lhisb                    ! backscattered atoms spatial distribution flag
      logical :: lhiseb                   ! backscattered atoms energy distribution flag
      logical :: lhiset                   ! transmitted atoms energy distribution flag
      logical :: lhisd                    ! displacement histogram flag
      logical :: lhisee                   ! deposited electronic energy histogram flag
      logical :: lhisee2                  ! 2d deposited electronic energy histogram flag
      logical :: lhisee3                  ! 3d deposited electronic energy histogram flag
      logical :: lhisne                   ! deposited nuclear energy histogram flag
      logical :: lhisne2                  ! 2d deposited nuclear energy histogram flag
      logical :: lhisne3                  ! 3d deposited nuclear energy histogram flag
      logical :: lhisiv                   ! I-V pair distribution flag
      logical :: lhism                    ! displacement moment histogram flag
      logical :: lhisp                    ! deposited momentum histogram flag
      logical :: lhist                    ! transmitted atoms spatial distribution flag
      logical :: lmom                     ! global moments flag
      logical :: lmomaab                  ! backscattered atoms angular moments flag (3D)
      logical :: lmomaat                  ! transmitted atoms angular moments flag (3D)
      logical :: lmomab                   ! backscattered atoms angular moments flag (2D)
      logical :: lmomat                   ! transmitted atoms angular moments flag (2D)
      logical :: lmomb                    ! backscattered atoms spatial moments flag
      logical :: lmomeb                   ! backscattered atoms energy moments flag
      logical :: lmomee                   ! electronic energy deposition moments
      logical :: lmomet                   ! transmitted atoms energy moments flag
      logical :: lmomne                   ! nuclear energy deposition moments
      logical :: lmomt                    ! transmitted atoms spatial moments flag
      logical :: lperx                    ! periodic output flag for xdirection
      logical :: lpery                    ! periodic output flag for y direction
      logical :: lrbs                     ! rbs flag
      logical :: lsurf                    ! surface output flag
      logical :: ltra(-1:5)               ! trajectory output flag
      logical :: shorttra                 ! shorten TRA file flag

      integer(INT4) :: dpower             ! power of displacements for the histogram
      integer(INT4) :: epower             ! power of recoil energies for the histogram
      integer(INT4) :: mcoord             ! model for statistics coordinate system
      integer(INT4) :: nbox               ! max. number of 1-d histogram boxes
      integer(INT4) :: nbox2(2)           ! max. number of 2-d histogram boxes
      integer(INT4) :: nbox3(3)           ! max. number of 3-d histogram boxes
      integer(INT4) :: nboxa              ! max. number of anglular histogram boxes
      integer(INT4) :: nboxa2(2)          ! max. number of 2-d anglular histogram boxes
      integer(INT4) :: nboxd              ! max. number of displacement/pair distribution
                                          !    boxes
      integer(INT4) :: nboxe              ! max. number of energy histogram boxes
      integer(INT4) :: nlatmom            ! maximum order of lateral moments
      integer(INT4) :: nmixmom            ! maximum order of mixed moments

      real(REAL8) :: amaxa                ! max. acceptance angle [deg]
      real(REAL8) :: amina                ! min. acceptance angle [deg]
      real(REAL8) :: dira(3)              ! reference direction for acceptance angle
      real(REAL8) :: emaxa                ! max. acceptance energy
      real(REAL8) :: emina                ! min. acceptance energy
      real(REAL8) :: fachis(3)            ! factor to convert counts to 1-d/2-d/3-d histogram
                                          !    entries
      real(REAL8) :: posivmax             ! max. depth for pair distribution
      real(REAL8) :: posivmin             ! min. depth for pair distribution
      real(REAL8) :: wbox                 ! initial 1-d histogram box width
      real(REAL8) :: wbox2(2)             ! initial 2-d histogram box width
      real(REAL8) :: wbox3(3)             ! initial 3-d histogram box width
      real(REAL8) :: wboxa                ! initial anglular histogram box width
      real(REAL8) :: wboxa2(2)            ! initial 2-d anglular histogram box width
      real(REAL8) :: wboxd                ! initial box width for pair distribution
      real(REAL8) :: wboxe                ! initial energy histogram box width
      real(REAL8) :: xcyl                 ! coordinates of axis for cylinder coordinates
      real(REAL8) :: ycyl                 !

   end type Out_type

   type(Out_type), public :: out                  ! out structure
!
!=======================================================================
!
contains
!
!=======================================================================
!
   subroutine readOut( )
!             =======
!
!     Read &output record from the input file
!
!  Arguments:
!
!     None
!
!  Local Variables:
!
      logical :: lbint     ! obsolete
      logical :: lcell
      logical :: lcyl
      logical :: ldef
      logical :: ldist
      logical :: lendp     ! obsolete
      logical :: lextendx(2)
      logical :: lextendy(2)
      logical :: lextendz(2)
      logical :: lhis
      logical :: lhis2
      logical :: lhis2b
      logical :: lhis2t
      logical :: lhis3
      logical :: lhisa2b
      logical :: lhisa2t
      logical :: lhisaab
      logical :: lhisaat
      logical :: lhisab
      logical :: lhisat
      logical :: lhisb
      logical :: lhisd
      logical :: lhise     ! obsolete
      logical :: lhiseb
      logical :: lhisee
      logical :: lhisee2
      logical :: lhisee3
      logical :: lhiset
      logical :: lhisne
      logical :: lhisne2
      logical :: lhisne3
      logical :: lhisiv
      logical :: lhism
      logical :: lhisp
      logical :: lhist
      logical :: lint      ! obsolete
      logical :: llincoll  ! obsolete
      logical :: lmom
      logical :: lmomaab
      logical :: lmomaat
      logical :: lmomab
      logical :: lmomat
      logical :: lmomb
      logical :: lmome     ! obsolete
      logical :: lmomeb
      logical :: lmomee
      logical :: lmomet
      logical :: lmomne
      logical :: lmomt
      logical :: lperx     ! obsolete
      logical :: lpery
      logical :: lrbs
      logical :: lsurf
      logical :: ltra
      logical :: shorttra

      integer(INT4) :: dpower
      integer(INT4) :: epower
      integer(INT4) :: itra(8)
      integer(INT4) :: mcoord
      integer(INT4) :: nbox
      integer(INT4) :: nbox2(2)
      integer(INT4) :: nbox3(3)
      integer(INT4) :: nboxa
      integer(INT4) :: nboxa2(2)
      integer(INT4) :: nboxd
      integer(INT4) :: nboxe
      integer(INT4) :: nionbk
      integer(INT4) :: niontra   ! obsolete
      integer(INT4) :: nlatmom
      integer(INT4) :: nmixmom
      integer(INT4) :: nupdout   ! obsolete

      real(REAL8) :: amaxa
      real(REAL8) :: amina
      real(REAL8) :: emaxa
      real(REAL8) :: emina
      real(REAL8) :: posivmax
      real(REAL8) :: posivmin
      real(REAL8) :: rotata                     ! analyzer rotation angle
      real(REAL8) :: tilta                      ! analyzer tilt angle
      real(REAL8) :: wbox
      real(REAL8) :: wbox2(2)
      real(REAL8) :: wbox3(3)
      real(REAL8) :: wboxa
      real(REAL8) :: wboxa2(2)
      real(REAL8) :: wboxd
      real(REAL8) :: wboxe
      real(REAL8) :: wboxi                      ! same as wbox (kept for compatibility)
      real(REAL8) :: wbox2i(2)                  ! same as wbox2 (kept for compatibility)
      real(REAL8) :: wbox3i(3)                  ! same as wbox3 (kept for compatibility)
      real(REAL8) :: wboxai                     ! same as wboxa (kept for compatibility)
      real(REAL8) :: wboxa2i(2)                 ! same as wboxa2 (kept for compatibility)
      real(REAL8) :: wboxdi                     ! same as wboxd (kept for compatibility)
      real(REAL8) :: wboxei                     ! same as wboxe (kept for compatibility)

      namelist /output/ lcell, lcyl, ldef, ldist, lextendx, lextendy, lextendz,  &
                        lhis, lhis2, lhis2b, lhis2t, lhis3,  &
                        lhisa2b, lhisa2t, lhisaab, lhisaat, lhisab, lhisat, lhisb,  &
                        lhisd, lhiseb, lhisee, lhisee2, lhisee3, lhiset,  &
                        lhisne, lhisne2, lhisne3, lhisiv, lhism, lhisp, lhist,  &
                        lmom, lmomaab, lmomaat, lmomab, lmomat, lmomb,  &
                        lmomeb, lmomee, lmomet, lmomne, lmomt,  &
                        lrbs, lsurf, ltra, shorttra,  &
                        dpower, epower, itra, mcoord, nbox, nbox2, nbox3, nboxa,  &
                        nboxa2, nboxd, nboxe, nlatmom, nmixmom,  &
                        amaxa, amina, emaxa, emina, posivmax, posivmin,  &
                        rotata, tilta, wbox, wbox2, wbox3, wboxa, wboxa2, wboxd, wboxe,  &
                        wboxi, wbox2i, wbox3i, wboxai, wboxa2i, wboxdi, wboxei
      ! obsolete parameters
      namelist /output/ lbint, lendp, lhise, lint, llincoll, lperx, lpery, lmome,  &
                        nionbk, niontra, nupdout

      character(256) :: msg
      integer(INT4) :: ios, i
!
!  Start of Executable Code:
!  ------------------------
!
!  Set default values
!
      lbint = .false.   ! obsolete
      lcell = .false.
      lcyl = .false.
      ldef = .false.
      ldist = .false.
      lendp = .false.   ! obsolete
      lextendx(:) = .false.
      lextendy(:) = .false.
      lextendz(:) = .false.
      lhis = .false.
      lhis2 = .false.
      lhis2b = .false.
      lhis2t = .false.
      lhisa2b = .false.
      lhisa2t = .false.
      lhis3 = .false.
      lhisaab = .false.
      lhisaat = .false.
      lhisab = .false.
      lhisat = .false.
      lhisb = .false.
      lhisd = .false.
      lhise = .false.   ! obsolete
      lhiseb = .false.
      lhisee = .false.
      lhiset = .false.
      lhisee2 = .false.
      lhisee3 = .false.
      lhisne = .false.
      lhisne2 = .false.
      lhisne3 = .false.
      lhisiv = .false.
      lhism = .false.
      lhisp = .false.
      lhist = .false.
      lint = .false.    ! obsolete
      llincoll = .false.! obsolete
      lmom = .true.
      lmomaab = .false.
      lmomaat = .false.
      lmomab = .false.
      lmomat = .false.
      lmomb = .false.
      lmome = .false.   ! obsolete
      lmomeb = .false.
      lmomee = .false.
      lmomet = .false.
      lmomne = .false.
      lmomt = .false.
      lperx = any( geom%bndtypex == ENUM%bndtype%periodic ) ! obsolete
      lpery = any( geom%bndtypey == ENUM%bndtype%periodic ) ! obsolete
      lrbs = .false.
      lsurf = .false.
      ltra = .false.
      out%ltra(:) = .false.
      shorttra = .false.

      call undefine( dpower )
      call undefine( epower )
      call undefine( itra(:) )
      call undefine( mcoord )
      call undefine( nbox )
      call undefine( nbox2(:) )
      call undefine( nbox3(:) )
      call undefine( nboxa )
      call undefine( nboxa2(:) )
      call undefine( nboxd )
      call undefine( nboxe )
      call undefine( nionbk )    ! obsolete
      call undefine( niontra )   ! obsolete
      call undefine( nlatmom )
      call undefine( nmixmom )
      call undefine( nupdout )   ! obsolete

      call undefine( amaxa )
      call undefine( amina )
      call undefine( emaxa )
      call undefine( emina )
      call undefine( posivmax )
      call undefine( posivmin )
      call undefine( rotata )
      call undefine( tilta )
      call undefine( wbox )
      call undefine( wbox2(:) )
      call undefine( wbox3(:) )
      call undefine( wboxa )
      call undefine( wboxa2(:) )
      call undefine( wboxd )
      call undefine( wboxe )
      call undefine( wboxi )
      call undefine( wbox2i(:) )
      call undefine( wbox3i(:) )
      call undefine( wboxai )
      call undefine( wboxa2i(:) )
      call undefine( wboxdi )
      call undefine( wboxei )
!
!  Read &output record(s) from the input file
!
      call openFile( 11, 'INP', 'formatted', 'read' )
!          ========  i     i         i         i

      do
         read( 11, nml=output, iostat=ios, iomsg=msg )
         if( ios < 0 )  exit
         if( ios > 0 ) then
            print '(/a)',  trim( msg )
            goto 900
         end if

         if( defined( itra(size( itra )) ) )  goto 941
         do i = 1, definedlength( itra(:) )     ! this has to be done within the loop so all
            if( itra(i) < -1 .or. itra(i) > 5 )  goto 954   !    records containing itra are
            out%ltra(itra(i)) = .true.                      !    considered
         end do

      end do

      close( 11 )
!
!  Copy the parameters to the out structure
!
      out%lcell = lcell
      out%lcyl = lcyl
      out%ldef = ldef
      out%ldist = ldist
      out%lextendx(:) = lextendx(:)
      out%lextendy(:) = lextendy(:)
      out%lextendz(:) = lextendz(:)
      out%lhis = lhis
      out%lhis2 = lhis2
      out%lhis2b = lhis2b
      out%lhis2t = lhis2t
      out%lhis3 = lhis3
      out%lhisa2b = lhisa2b
      out%lhisa2t = lhisa2t
      out%lhisaab = lhisaab
      out%lhisaat = lhisaat
      out%lhisab = lhisab
      out%lhisat = lhisat
      out%lhisb = lhisb
      out%lhisd = lhisd
      out%lhiseb = lhiseb
      out%lhisee = lhisee
      out%lhisee2 = lhisee2
      out%lhisee3 = lhisee3
      out%lhiset = lhiset
      out%lhisne = lhisne
      out%lhisne2 = lhisne2
      out%lhisne3 = lhisne3
      out%lhisiv = lhisiv
      out%lhism = lhism
      out%lhisp = lhisp
      out%lhist = lhist
      out%lmom = lmom
      out%lmomaab = lmomaab
      out%lmomaat = lmomaat
      out%lmomab = lmomab
      out%lmomat = lmomat
      out%lmomb = lmomb
      out%lmomeb = lmomeb
      out%lmomee = lmomee
      out%lmomet = lmomet
      out%lmomne = lmomne
      out%lmomt = lmomt
      out%lrbs = lrbs
      out%lsurf = lsurf
      out%shorttra = shorttra

      if( ltra .and. .not.any( out%ltra(:) ) )  out%ltra(:) = .true.

      out%dpower = dpower
      out%epower = epower
      out%mcoord = mcoord
      out%nbox = nbox
      out%nbox2(:) = nbox2(:)
      out%nbox3(:) = nbox3(:)
      out%nboxa = nboxa
      out%nboxa2(:) = nboxa2(:)
      out%nboxd = nboxd
      out%nboxe = nboxe
      if( defined( nionbk ) )  goto 980
      out%nlatmom = nlatmom
      out%nmixmom = nmixmom

      out%amaxa = amaxa
      out%amina = amina
      out%emaxa = emaxa
      out%emina = emina
      out%posivmin = posivmin
      out%posivmax = posivmax
      if( .not.defined( rotata ) )  rotata = ZERO
      if( .not.defined( tilta ) )  tilta = ZERO
      out%dira(1) = sin(tilta*PI/180)*cos(rotata*PI/180)
      out%dira(2) = sin(tilta*PI/180)*sin(rotata*PI/180)
      out%dira(3) = cos(tilta*PI/180)
      if( defined( wbox ) ) then
         out%wbox = wbox
      else
         out%wbox = wboxi
      end if
      if( defined( wbox2(1) ) ) then
         out%wbox2(:) = wbox2(:)
      else
         out%wbox2(:) = wbox2i(:)
      end if
      if( defined( wbox3(1) ) ) then
         out%wbox3(:) = wbox3(:)
      else
         out%wbox3(:) = wbox3i(:)
      end if
      if( defined( wboxa ) ) then
         out%wboxa = wboxa
      else
         out%wboxa = wboxai
      end if
      if( defined( wboxa2(1) ) ) then
         out%wboxa2(:) = wboxa2(:)
      else
         out%wboxa2(:) = wboxa2i(:)
      end if
      if( defined( wboxd ) ) then
         out%wboxd = wboxd
      else
         out%wboxd = wboxdi
      end if
      if( defined( wboxe ) ) then
         out%wboxe = wboxe
      else
         out%wboxe = wboxei
      end if
!
!  Signal to &setup whether niontra is to be ignored
!
      if( .not.any(out%ltra(:)) )  call undefine( niontra_setup )
!
!  Some checks that cannot be done in initOut
!
      if( lbint )  goto 950
      if( lendp )  goto 951
      if( lhise )  goto 955
      if( lint )  goto 952
      if( llincoll )  goto 956
      if( lmome )  goto 958
      if( lperx .neqv. any( geom%bndtypex == ENUM%bndtype%periodic ) )  goto 961
      if( lpery .neqv. any( geom%bndtypey == ENUM%bndtype%periodic ) )  goto 962
      if( defined( niontra ) )  goto 959
      if( defined( nupdout ) )  goto 960
      if( tilta < ZERO .or. tilta > 90 )  goto 974
!
!  Set default values of geom%nbox2, geom%wbox2
!
      if( ldyn ) then
         if( .not.defined( geom%nbox ) )  geom%nbox = out%nbox
         if( .not.defined( geom%wbox ) )  geom%wbox = out%wbox
         where( .not.defined( geom%nbox2(:) ) )  geom%nbox2(:) = out%nbox2(:)
         where( .not.defined( geom%wbox2(:) ) )  geom%wbox2(:) = out%wbox2(:)
      end if

      return
!
!  Error messages
!
900   error stop '&output: error during reading namelist'

941   error stop '&output: Too many values of itra'

950   error stop '&output: lbint is obsolete'
951   error stop '&output: lendp is obsolete. Use itra instead'
952   error stop '&output: lint is obsolete. Use itra for partial replacement'
953   error stop '&output: error during reading itra'
954   error stop '&output: -1<=itra<=5 required'
955   error stop '&output: lhise is obsolete. Use lhisne instead'
956   error stop '&output: llincoll is obsolete'
958   error stop '&output: lmome is obsolete. Use lmomne instead'
959   error stop '&output: niontra is obsolete. Use niontra of &setup instead'
960   error stop '&output: nupdout is obsolete. Use nionout of &setup instead'
961   error stop '&output: lperx is obsolete. Use lextendx instead'
962   error stop '&output: lpery is obsolete. Use lextendy instead'

974   error stop '&output: 0<=tilta<=90 required'

980   error stop '&output: nionbk must be specified on &setup record'

   end subroutine readOut
!
!=================================================================================
!
   subroutine initOut( )
!             =======
!
!      Initialize the out structure
!
!  Arguments:
!
!     None
!
!  Local Variables:
!
      integer(INT4) :: ir
      real(REAL8) :: dxinit, dyinit
!
!  Start of Executable Code:
!  ------------------------
!
!  Set default values
!
      if( dam%ldam .and. dam%ndamdim==1 )  out%lhis = .true.
      if( dam%ldam .and. dam%ndamdim==2 )  out%lhis2 = .true.
      if( dam%ldam .and. dam%ndamdim==3 )  out%lhis3 = .true.

      if( .not.defined( out%dpower ) )  out%dpower = 3
      if( .not.defined( out%epower ) )  out%epower = 3
      if( .not.defined( out%mcoord ) )  out%mcoord = 1
      if( ldyn ) then
         if( .not.defined( out%nbox ) )  out%nbox = geom%nbox
         if( .not.defined( out%nbox2(1) ) )  out%nbox2(1) = geom%nbox2(1)
         if( .not.defined( out%nbox2(2) ) )  out%nbox2(2) = geom%nbox2(2)
      else
         if( .not.defined( out%nbox ) )  out%nbox = 32
         if( .not.defined( out%nbox2(1) ) )  out%nbox2(1) = 36
         if( .not.defined( out%nbox2(2) ) )  out%nbox2(2) = 24
      end if
      if( .not.defined( out%nbox3(1) ) )  out%nbox3(1) = 36
      if( .not.defined( out%nbox3(2) ) )  out%nbox3(2) = 36
      if( .not.defined( out%nbox3(3) ) )  out%nbox3(3) = 24
      if( .not.defined( out%nboxa ) )  out%nboxa = 90
      if( .not.defined( out%nboxa2(1) ) )  out%nboxa2(1) = 50
      if( .not.defined( out%nboxa2(2) ) )  out%nboxa2(2) = out%nboxa2(1)
      if( .not.defined( out%nboxd ) )  out%nboxd = 100
      if( .not.defined( out%nboxe ) )  out%nboxe = 100
      if( .not.defined( out%nlatmom ) )  out%nlatmom = 2
      if( .not.defined( out%nmixmom ) )  out%nmixmom = 0

      if( .not.defined( out%amaxa ) )  out%amaxa = 180
      if( .not.defined( out%amina ) )  out%amina = ZERO
      if( .not.defined( out%emaxa ) )  out%emaxa = ion%e
      if( .not.defined( out%emina ) )  out%emina = ZERO
      if( .not.defined( out%posivmax ) )  out%posivmax = BIG
      if( .not.defined( out%posivmin ) )  out%posivmin = ZERO
      if( ldyn ) then
         if( .not.defined( out%wbox ) )  out%wbox = geom%wbox
         if( .not.defined( out%wbox2(1) ) )  out%wbox2(1) = geom%wbox2(1)
         if( .not.defined( out%wbox2(2) ) )  out%wbox2(2) = geom%wbox2(2)
      else
         if( .not.defined( out%wbox ) )  out%wbox = 100./128.
         if( .not.defined( out%wbox2(1) ) )  out%wbox2(1) = out%wbox
         if( .not.defined( out%wbox2(2) ) )  out%wbox2(2) = out%wbox2(1)
      end if

      if( .not.defined( out%wbox3(1) ) )  out%wbox3(1) = out%wbox
      if( .not.defined( out%wbox3(2) ) )  out%wbox3(2) = out%wbox3(1)
      if( .not.defined( out%wbox3(3) ) )  out%wbox3(3) = out%wbox3(1)
      if( .not.defined( out%wboxa ) )  out%wboxa = ONE
      if( .not.defined( out%wboxa2(1) ) )  out%wboxa2(1) = 0.01
      if( .not.defined( out%wboxa2(2) ) )  out%wboxa2(2) = out%wboxa2(1)
      if( .not.defined( out%wboxd ) )  out%wboxd = 0.1
      if( .not.defined( out%wboxe ) ) then
         if( out%nboxe == 0 ) then
            out%wboxe = ion%e**(ONE/out%epower) / 100
         else
            out%wboxe = ion%e**(ONE/out%epower) / out%nboxe
         end if
      end if
!
!  Check validity of the parameters
!
      if( nthread > 1 ) then
         if( any( out%ltra(:) ) )  goto 921
      end if

      if( out%lcell .and. .not.ldyn )  goto 959
      if( out%lcyl ) then
         if( ion%xinit(1) /= ion%xinit(2) .or.  &
             ion%yinit(1) /= ion%yinit(2) )  goto 960
      endif
      if( out%ldef .and. .not. dam%ldet )  goto 958
      if( any( out%lextendx(:) .and. geom%bndtypex(:)==ENUM%bndtype%free ) .or.  &
          any( out%lextendy(:) .and. geom%bndtypey(:)==ENUM%bndtype%free ) .or.  &
          any( out%lextendz(:) .and. geom%bndtypez(:)==ENUM%bndtype%free ) )  goto 910
      if( geom%lcell .and.  &
          ( any( out%lextendx(:) .and. geom%bndtypex(:)==ENUM%bndtype%free ) .or.  &
            any( out%lextendy(:) .and. geom%bndtypey(:)==ENUM%bndtype%free ) .or.  &
            any( out%lextendz(:) .and. geom%bndtypez(:)==ENUM%bndtype%free ) ) ) goto 911
      if( ( any( geom%bndtypex(:) == ENUM%bndtype%reflect .or.  &
                 geom%bndtypex(:) == ENUM%bndtype%periodic ) .and.  &
            out%lextendx(1) .neqv. out%lextendx(2) ) .or.  &
          ( any( geom%bndtypey(:) == ENUM%bndtype%reflect .or.  &
                 geom%bndtypey(:) == ENUM%bndtype%periodic ) .and.  &
            out%lextendy(1) .neqv. out%lextendy(2) ) .or.  &
          ( any( geom%bndtypez(:) == ENUM%bndtype%reflect .or.  &
                 geom%bndtypez(:) == ENUM%bndtype%periodic ) .and.  &
            out%lextendz(1) .neqv. out%lextendz(2) ) ) goto 912
      if( .not. out%lhis .and. usehis  )  goto 940
      if( .not. out%lhis .and. dam%ldam .and. dam%ndamdim == 1 )  goto 941
      if( .not. out%lhis2 .and. dam%ldam .and. dam%ndamdim == 2 )  goto 942
      if( .not. out%lhis3 .and. dam%ldam .and. dam%ndamdim == 3 )  goto 943
      if( out%lhis2 .and. dam%ldet )  goto 930
      if( out%lhis3 .and. dam%ldet )  goto 931
      if( out%lhisiv .and. (.not.dam%lcap.or.ldyn.or.usecell) )  goto 950
      if( out%lmom .and. dam%ldet )  goto 932
      if( out%lrbs ) then
         if( ldyn .or. usecell )  goto 979
         do ir = 1, nr
            if( any( mat(ir)%vib( mat(ir)%ielem(1:mat(ir)%nelemr ) ) <= ZERO) )  goto 9791
         end do
      end if
      if( out%lsurf .and. .not.ldyn )  goto 9591

      if( out%dpower <= 0 )  goto 987
      if( out%epower <= 0 )  goto 985
      if( out%mcoord < -1 .or. out%mcoord > 5 )  goto 961
      if( out%mcoord == 3 .and. ion%ranrot )  goto 962
      if( out%nlatmom < 0 .or. out%nlatmom > 4)  goto 933
      if( out%nmixmom < 0 .or. out%nmixmom > 4)  goto 934

      if( out%amina < ZERO .or. out%amaxa > 180 .or. out%amina > out%amaxa )  goto 975
      if( out%emina < ZERO .or. out%emina > out%emaxa )  goto 976
      if( out%posivmin >= out%posivmax )  goto 972
      if( out%nbox < 0 )  goto 9800
      if( any( out%nbox2(:) < 0 ) )  goto 9801
      if( any( out%nbox3(:) < 0 ) )  goto 9802
      if( out%nboxa < 0 )  goto 9803
      if( any( out%nboxa2(:) < 0 ) )  goto 9806
      if( out%nboxd < 0 )  goto 9804
      if( out%nboxe < 0 )  goto 9805
      if( out%wbox <= ZERO )  goto 980
      if( any( out%wbox2(:) <= ZERO ) )  goto 981
      if( any( out%wbox3(:) <= ZERO ) )  goto 982
      if( out%wboxa <= ZERO )  goto 984
      if( any( out%wboxa2(:) <= ZERO ) )  goto 9841
      if( out%wboxd <= ZERO )  goto 986
      if( out%wboxe <= ZERO )  goto 983
!
!  Modify parameters and determine derived parameters
!
      out%emaxa = out%emaxa * 1000
      out%emina = out%emina * 1000
!
!  If vacancy distribution is to be computed and a crystalline region exists,
!  the wbox* have to be multiples of xscale in order to avoid artificial
!  jitter in the net point defect concentrations. Since we have one histogram
!  grid for all regions, only one crystalline region can be considered
!
      if( dam%lrcoil .and. .not.usecell ) then
         do ir = 1, nr
            if( .not.mat(ir)%lxtal ) cycle
            out%wbox = mat(ir)%xscale * max( 1, nint( out%wbox / mat(ir)%xscale ) )
            out%wbox2(:) = mat(ir)%xscale * max( 1, nint( out%wbox2(:) / mat(ir)%xscale ) )
            out%wbox3(:) = mat(ir)%xscale * max( 1, nint( out%wbox3(:) / mat(ir)%xscale ) )
            exit
         end do
      end if
!
!  fachis
!
      if( ion%dxinit == ZERO ) then
         dxinit = ONE
      else
         dxinit = ion%dxinit
      end if
      if( ion%dyinit == ZERO ) then
         dyinit = ONE
      else
         dyinit = ion%dyinit
      end if
      if( ion%doseunits == 'cm-2' ) then
         out%fachis(1) = ion%ddose
         out%fachis(2) = ion%ddose * dxinit
         out%fachis(3) = ion%ddose * dxinit * dyinit
      else if( ion%doseunits == 'cm-1' ) then
         out%fachis(1) = ion%ddose / dxinit
         out%fachis(2) = ion%ddose
         out%fachis(3) = ion%ddose * dyinit
      else
         out%fachis(1) = ion%ddose / (dxinit * dyinit)
         out%fachis(2) = ion%ddose / dyinit
         out%fachis(3) = ion%ddose
      end if

      out%xcyl = ion%xinit(1)
      out%ycyl = ion%yinit(1)

      return

910   error stop '&output: lextend*=t not allowed at free boundaries'
911   error stop '&output: lextend*=t not allowed at periodic boundaries if usecell=t or ldyn=t'
912   error stop '&output: lextend*(1)=lextend*(2) required at periodic or reflective ' //  &
                 'boundaries'

921   error stop '&output: ltra=t or itra not allowed for parallel execution'

930   error stop '&output: lhis2=t not allowed with deterministic defects'
931   error stop '&output: lhis3=t not allowed with deterministic defects'
932   error stop '&output: lmom=t not allowed with deterministic defects'
933   error stop '&output: 0<=nlatmom<=4 required'
934   error stop '&output: 0<=nmixmom<=4 required'

940   error stop '&output: usehis=t requires lhis=t'
941   error stop '&output: ndamdim=1 on &damage requires lhis if ldam=t'
942   error stop '&output: ndamdim=2 on &damage requires lhis2 if ldam=t'
943   error stop '&output: ndamdim=3 on &damage requires lhis3 if ldam=t'

950   error stop '&output: lhisiv=t requires lcap=t, ldyn=f'
958   error stop '&output: ldef=t requires deterministic model of damage. &damage ldet=t'
959   error stop '&output: lcell=t requires dynamic simulation mode: &setup ldyn=t'
9591  error stop '&output: lsurf=t requires dynamic simulation mode: &setup ldyn=t'

960   error stop '&output: lcyl requires xinit(1)=xinit(2), yinit(1)=yinit(2)'
961   error stop '&output: 1<=mcoord<=5 required'
962   error stop '&output: mcoord=3 must not be specified for ranrot=t'
963   error stop '&output: lextendx=t, lextendy=t, lextendz=t not allowed at free boundaries'

972   error stop '&output: posivmin<posivmax required'
975   error stop '&output: 0<=amina<=amaxa<=180 required'
976   error stop '&output: 0<=emina<=emaxa required'
979   error stop '&output: lrbs=t and ldyn=t or usecell=t not allowed'
9791  error stop '&output: lrbs=t and vib=0 not allowed'

9800  error stop '&output: nbox>=0 required'
9801  error stop '&output: nbox2>=0 required'
9802  error stop '&output: nbox3>=0 required'
9803  error stop '&output: nboxa>=0 required'
9804  error stop '&output: nboxd>=0 required'
9805  error stop '&output: nboxe>=0 required'
9806  error stop '&output: nboxa2>=0 required'
980   error stop '&output: wbox>0 required'
981   error stop '&output: wbox2>0 required'
982   error stop '&output: wbox3>0 required'
983   error stop '&output: wboxe>0 required'
986   error stop '&output: wboxd>0 required'
984   error stop '&output: wboxa>0 required'
9841  error stop '&output: wboxa2>0 required'
985   error stop '&output: epower>0 required'
987   error stop '&output: dpower>0 required'

   end subroutine initOut
!
!=================================================================================
!
   subroutine writeOut( )
!             ========
!
!     Write data derived from the &output record to the output file
!
!  Arguments:
!
!     None
!
!  Local Variables:
!
      integer(INT4) :: i, itra, ngeom
!
!  Start of Executable Code:
!  ------------------------
!
!  Open the output file
!
      call openFile( 21, 'OUT', 'formatted', 'append' )
!          ========  i     i         i          i
!
!  Moments and histogram flags
!
      write( 21, '(/t2,a,t13)', advance='no' )  &
            '&output:'
      if( ldyn ) then
         if( out%lcell ) then
            write( 21, '(a)', advance='no' )  &
                  'lcell=t: CELL'
            if( ndim == 2 )  &
               write( 21, '(a)', advance='no' )  &
                     '2'
            if( niongeom == nion ) then
               write( 21, '(a)', advance='no' )  &
                  ' file is written at the end of the simulation'
            else
               ngeom = nion / niongeom
               write( 21, '(a)', advance='no' )  &
                  ' file is written after every'
               write( 21, '(a)', advance='no' )  &
                     str( ngeom )
               if( mod( ngeom, 10 ) == 1 ) then
                  write( 21, '(a)', advance='no' )  '-st'
               else if( mod( ngeom, 10 ) == 2 ) then
                  write( 21, '(a)', advance='no' )  '-nd'
               else if( mod( ngeom, 10 ) == 3 ) then
                  write( 21, '(a)', advance='no' )  '-rd'
               else
                  write( 21, '(a)', advance='no' )  '-th'
               end if
               write( 21, '(a)', advance='no' )  &
                     ' geometry update'
            end if
         end if
         write( 21, '(/t13)', advance='no' )
      end if
      write( 21, '(a)' )  &
            'lmom=' // str( out%lmom ) // &
            ', lmomee=' // str( out%lmomee ) //  &
            ', lmomne=' // str( out%lmomne ) // &
            ', nlatmom=' // str( out%nlatmom ) // &
            ', nmixmom=' // str( out%nmixmom )
      write( 21, '(t13,a)' )  &
            'lmomaab=' // str( out%lmomaab ) //  &
            ', lmomab=' // str( out%lmomab ) //  &
            ', lmomb=' // str( out%lmomb ) //  &
            ', lmomeb=' // str( out%lmomeb )
      write( 21, '(t13,a)' )  &
            'lmomaat=' // str( out%lmomaat ) //  &
            ', lmomat=' // str( out%lmomat ) //  &
            ', lmomt=' // str( out%lmomt ) //  &
            ', lmomet=' // str( out%lmomet )
      write( 21, '(t13,a)' )  &
            'lhis=' // str( out%lhis ) //  &
            ', lhis2=' // str( out%lhis2 ) //  &
            ', lhis3=' // str( out%lhis3 )
      write( 21, '(t13,a)' )  &
            'lhisee=' // str( out%lhisee ) //  &
            ', lhisee2=' // str( out%lhisee2 ) //  &
            ', lhisee3=' // str( out%lhisee3 ) //  &
            ', lhisne=' // str( out%lhisne ) //  &
            ', lhisne2=' // str( out%lhisne2 ) // &
            ', lhisne3=' // str( out%lhisne3 )
      write( 21, '(t13,a)' )  &
            'lhisaab=' // str( out%lhisaab ) //  &
            ', lhisab=' // str( out%lhisab ) //  &
            ', lhisa2b=' // str( out%lhisa2b ) //  &
            ', lhisb=' // str( out%lhisb ) //  &
            ', lhis2b=' // str( out%lhis2b ) //  &
            ', lhiseb=' // str( out%lhiseb )
      write( 21, '(t13,a)' )  &
            'lhisaat=' // str( out%lhisaat ) //  &
            ', lhisat=' // str( out%lhisat ) //  &
            ', lhisa2t=' // str( out%lhisa2t ) //  &
            ', lhist=' // str( out%lhist ) //  &
            ', lhis2t=' // str( out%lhis2t ) //  &
            ', lhiset=' // str( out%lhiset )
      if( out%lhiseb .or. out%lhiset )  &
         write( 21, '(t13,a)' )  &
               'Abscissas of the energy histogram for target atoms: energy**' //  &
               str( out%epower )
      if( out%lmomaab .or. out%lmomaat .or. out%lhisaab .or. out%lhisaat ) then
         write( 21, '(t13,a)' )  &
               'Reference direction=(' // str( out%dira(1), 'f99.3' ) //  &
               ',' // str( out%dira(2), 'f99.3' ) //  &
               ',' // str( out%dira(3), 'f99.3' ) // ')'
         if( out%amina > ZERO .or. out%amaxa < 180 .or.  &
             out%emina > ZERO .or. out%emaxa < ion%e ) then
            write( 21, '(t13,a)' )  &
                  'Backscattering and transmission statistics contains only atoms with'
            if( out%amina > ZERO .or. out%amaxa < 180 )  &
               write( 21, '(t13,a)' )  &
                     '   angles from the reference direction in the range (' //  &
                     str( out%amina, 'f99.1' ) //  &
                     ',' // str( out%amaxa, 'f99.1' ) // ')'
            if( out%emina > ZERO .or. out%emaxa < ion%e )  &
               write( 21, '(t13,a)' )  &
                     '   energies in the range (' // str( out%emina, 'f99.1' ) //  &
                     ',' // str( out%emaxa, 'f99.1' ) // ')'
         end if
      end if

      write( 21, '(/t13,a)', advance='no' )  &
            'lhisd=' // str( out%lhisd ) //  &
            ', lhisiv=' // str( out%lhisiv )
      if( out%lhisd )  &
         write( 21, '(a)', advance='no' )  &
               ', dpower=' // str( out%dpower )
      write( 21, '()' )

      if( out%lhisiv )  &
            write( 21, '(t13,a)' )  &
                  '   Pair distribution function' //  &
                  ' in z=(' // str( out%posivmin, 'es10.3' ) //  &
                  ',' // str( out%posivmax, 'es10.3' )  //  ')'

      write( 21, '(t13,a)' )  &
            'lhism=' // str( out%lhism ) //  &
            ', lhisp=' // str( out%lhisp )

      write( 21, '(t13,a)', advance='no' )  &
            'lcyl=' // str( out%lcyl )
      if( out%lcyl )  &
         write( 21, '(a)', advance='no' ) &
               ', xcyl=' // str( out%xcyl, 'f99.1' ) //  &
               ', ycyl=' // str( out%ycyl, 'f99.1' )
      write( 21, '(a)' )

      if( any( geom%bndtypex(:) == ENUM%bndtype%reflect ) .or.  &
          any( geom%bndtypex(:) == ENUM%bndtype%periodic ) ) then ! both lextendx(:) the same
         if( out%lextendx(1) ) then
            write( 21, '(a)', advance='no' )  &
                  'lextendx=t: Scores outside geometry domain are allowed'
         else
            if( any( geom%bndtypex(:) == ENUM%bndtype%periodic ) ) then
               write( 21, '(t13,a)' )  &
                     'lextendx=f: Scores are mapped into the geometric period'
            else if( any( geom%bndtypex(:) == ENUM%bndtype%extend ) ) then
               write( 21, '(t13,a)' )  &
                     'lextendx=f: Scores outside geometry domain are ignored'
            end if
         end if
      else
         if( out%lextendx(1) ) then
            write( 21, '(a)', advance='no' )  &
                  'lextendx(1)=t: Scores at x < lower limit of geometry domain are allowed'
         else if( geom%bndtypex(1) == ENUM%bndtype%extend ) then
            write( 21, '(t13,a)' )  &
                  'lextendx(1)=f: Scores at x < lower limit of geometry domain are ignored'
         end if
         if( out%lextendx(2) ) then
            write( 21, '(a)', advance='no' )  &
                  'lextendx(2)=t: Scores at x > upper limit of geometry domain are allowed'
         else if( geom%bndtypex(2) == ENUM%bndtype%extend ) then
            write( 21, '(t13,a)' )  &
                  'lextendx(2)=f: Scores at x < lower limit of geometry domain are ignored'
         end if
      end if

      if( any( geom%bndtypey(:) == ENUM%bndtype%reflect ) .or.  &
          any( geom%bndtypey(:) == ENUM%bndtype%periodic ) ) then ! both lextendy(:) the same
         if( out%lextendy(1) ) then
            write( 21, '(a)', advance='no' )  &
                  'lextendy=t: Scores outside geometry domain are allowed'
         else
            if( any( geom%bndtypey(:) == ENUM%bndtype%periodic ) ) then
               write( 21, '(t13,a)' )  &
                     'lextendy=f: Scores are mapped into the geometric period'
            else if( any( geom%bndtypey(:) == ENUM%bndtype%extend ) ) then
               write( 21, '(t13,a)' )  &
                     'lextendy=f: Scores outside geometry domain are ignored'
            end if
         end if
      else
         if( out%lextendy(1) ) then
            write( 21, '(a)', advance='no' )  &
                  'lextendy(1)=t: Scores at y < lower limit of geometry domain are allowed'
         else if( geom%bndtypey(1) == ENUM%bndtype%extend ) then
            write( 21, '(t13,a)' )  &
                  'lextendy(1)=f: Scores at y < lower limit of geometry domain are ignored'
         end if
         if( out%lextendy(2) ) then
            write( 21, '(a)', advance='no' )  &
                  'lextendy(2)=t: Scores at y > upper limit of geometry domain are allowed'
         else if( geom%bndtypey(2) == ENUM%bndtype%extend ) then
            write( 21, '(t13,a)' )  &
                  'lextendy(2)=f: Scores at y < lower limit of geometry domain are ignored'
         end if
      end if

      if( any( geom%bndtypez(:) == ENUM%bndtype%reflect ) .or.  &
          any( geom%bndtypez(:) == ENUM%bndtype%periodic ) ) then ! both lextendz(:) the same
         if( out%lextendz(1) ) then
            write( 21, '(a)', advance='no' )  &
                  'lextendz=t: Scores outside geometry domain are allowed'
         else
            if( any( geom%bndtypez(:) == ENUM%bndtype%periodic ) ) then
               write( 21, '(t13,a)' )  &
                     'lextendz=f: Scores are mapped into the geometric period'
            else if( any( geom%bndtypez(:) == ENUM%bndtype%extend ) ) then
               write( 21, '(t13,a)' )  &
                     'lextendz=f: Scores outside geometry domain are ignored'
            end if
         end if
      else
         if( out%lextendz(1) ) then
            write( 21, '(a)', advance='no' )  &
                  'lextendz(1)=t: Scores at z < lower limit of geometry domain are allowed'
         else if( geom%bndtypez(1) == ENUM%bndtype%extend ) then
            write( 21, '(t13,a)' )  &
                  'lextendz(1)=f: Scores at z < lower limit of geometry domain are ignored'
         end if
         if( out%lextendz(2) ) then
            write( 21, '(a)', advance='no' )  &
                  'lextendz(2)=t: Scores at z > upper limit of geometry domain are allowed'
         else if( geom%bndtypez(2) == ENUM%bndtype%extend ) then
            write( 21, '(t13,a)' )  &
                  'lextendz(2)=f: Scores at z < lower limit of geometry domain are ignored'
         end if
      end if

      write( 21, '(t13,a)' )  &
            'ldef=' // str( out%ldef )

      if( out%lrbs ) then
         write( 21, '(t13,a)' )  &
               'lrbs=t, RBS detector angle=' // str( 180.-(180/PI)*acos(out%dira(3)), 'f99.2' )
      else
         write( 21, '(t13,a)' ) 'lrbs=f'
      endif

      write( 21, '(t13,a)' )  &
            'lsurf=' // str( out%lsurf ) //  &
            ', ldist=' // str( out%ldist )

      if( any( out%ltra(:) ) ) then
         write( 21, '(t13,a)', advance='no' )  &
               'Events '
         i = 0
         do itra = -1, 5
            if( out%ltra(itra) ) then
               i = i + 1
               if( i > 1 )  &
                  write( 21, '(a)', advance='no' )  &
                        ','
               write( 21, '(a)', advance='no' )  &
                     str( itra )
            end if
         end do
         write( 21, '(a)', advance='no' )  &
               ' are written to the TRA file'
         if( niontra_setup /= 1 ) then
            write( 21, '(a)' )  &
                  ' for every ' // str( niontra_setup ) // 'th ion'
         else
            write( 21, '(a)' )
         end if
         if( out%ltra(2) ) then
            write( 21, '(t13,a)', advance='no' )  &
                  'shorttra=' // str( out%shorttra ) // ': '
            if( out%shorttra ) then
               write( 21, '(a)' )  &
                     'Collision points without recoils are not recorded'
            else
               write( 21, '(a)' )  &
                     'All collision points are recorded'
            end if
         end if
      else
         write( 21, '(t13,a)' )  &
               'ltra=f'
      end if
!
!  Properties of all/most output
!
      if( dam%ldam ) then
         if( dam%lrcoil ) then
            write( 21, '(/t13,a)' )  &
                  'Ions, interstitials and vacancies are recorded'
         else
            write( 21, '(/t13,a)' )  &
                  'Ions and Frenkel pairs are recorded'
         end if
      else
         write( 21, '(/t13,a)' )  &
            'Only ions are recorded'
      end if

      write( 21, '(t13,a)', advance='no' )  &
            'mcoord=' // str( out%mcoord ) // ': '
      select case( out%mcoord )
         case( 1 )
            write( 21, '(a)' ) &
                  'The coordinate system is tied to the geometry'
         case( 2 )
            write( 21, '(a)' ) &
                  'The coordinate system is tied to the actual incidence direction'
         case( 3 )
            write( 21, '(a)' ) &
                  'The coordinate system is tied to the nominal incidence direction'
         case( 4 )
            write( 21, '(a)' ) &
                  'The coordinate system is tied to the analyzer direction'
            write( 21, '(t13,a)' )  &
                  '   Analyzer direction=(' // str( out%dira(1), 'f99.3' ) //  &
                  ',' // str( out%dira(2), 'f99.3' ) //  &
                  ',' // str( out%dira(3), 'f99.3' ) // ')'
      end select
!
!  Histogram descriptions
!
      write( 21, '(t13,a)' ) &
            'Maximum number of boxes and initial box widths for'

      if( out%lhis .or. dam%lamzon .or. out%lhism .or. out%lhisee .or.  &
          out%lhisne .or. out%lhisp .or. out%lhisb .or. out%lhist )  &
         write( 21, '(t13,a)' )  &
               '   1d spatial histograms: nbox=' // str( out%nbox ) //  &
               ', wbox=' // str( out%wbox, 'f99.3' )
      if( out%lhis2 .or. out%lhisee2 .or. out%lhisne2 .or. out%lhis2b .or. out%lhis2t)  &
         write( 21, '(t13,a)' ) &
               '   2d spatial histograms: nbox2=(' // str( out%nbox2(1) ) //  &
               ',' // str( out%nbox2(2) ) //  &
               '), wbox2=(' // str( out%wbox2(1), 'f99.3' ) //  &
               ',' // str( out%wbox2(2), 'f99.3' ) // ')'
      if( out%lhis3 .or. out%lhisee3 .or. out%lhisne3 ) then
         write( 21, '(t13,a)' ) &
               '   3d spatial histograms: nbox3=(' // str( out%nbox3(1) ) //  &
               ',' // str( out%nbox3(2) ) //  &
               ',' // str( out%nbox3(3) ) // ')'
         write( 21, '(t13,a)' ) &
               '                          wbox3=(' // str( out%wbox3(1), 'f99.3' ) //  &
               ',' // str( out%wbox3(2), 'f99.3' ) //  &
               ',' // str( out%wbox3(3), 'f99.3' ) // ')'
      end if
      if( out%lhiseb .or. out%lhiset )  &
         write( 21, '(t13,a)' )  &
               '   1d energy histograms: nboxe=' // str( out%nboxe ) //  &
               ', wboxe=' // str( out%wboxe, 'f99.3' )
      if( out%lhisab .or. out%lhisaab .or. out%lhisat .or. out%lhisaat )  &
         write( 21, '(t13,a)' )  &
               '   1d angular histograms: nboxa=' // str( out%nboxa ) //  &
               ', wboxa=' // str( out%wboxa, 'f99.3' )
      if( out%lhisa2b .or. out%lhisa2t )  &
         write( 21, '(t13,a)' )  &
               '   2d angular histograms: nboxa2=(' // str( out%nboxa2(1) ) //  &
               ',' // str( out%nboxa2(2) ) // ')' //  &
               ', wboxa2=(' // str( out%wboxa2(1), 'f99.3' ) //  &
               ',' // str( out%wboxa2(2), 'f99.3' ) // ')'
      if( out%lhisiv .or. out%lhisd )  &
         write( 21, '(t13,a)' )  &
               '   1d displacement histograms: nboxd=' // str( out%nboxd ) //  &
               ', wboxd=' // str( out%wboxd, 'f99.3' )

      write( 21, '()' )

      close( 21 )

      end subroutine writeOut
!
!=======================================================================
!
      end module Param_Out_mod
