module Param_Setup_mod
!
!  Holds the dimensions of various arrays and manages their input and output.
!
!  Note: The data are not put into a structure since some are used very often.
!
!  Used Modules:
!  ============
!
!$ use omp_lib,            only: omp_get_max_threads

   use Const_mod,          only: INT4, REAL8
   use String_mod,         only: str
   use Undefine_mod,       only: undefine, defined
   use Open_mod,           only: openFile
   use Math_IntMath_mod,   only: ggT
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
   public :: readSetup              ! read &setup record from input file
!            =========
   public :: initSetup              ! initialize parameters that need input from other records
!            =========
   public :: writeSetup             ! write data derived from &setup record
!            ==========             !    output file
!
!  Public Data:
!  ===========
!
   logical, public :: ldamdyn             ! consider damage dynamically
   logical, public :: ldamdyn_set         ! flag indicating ldamdyn has been set in &setup
   logical, public :: ldyn                ! consider dynamic changes of the target
   logical, public :: lsput               ! ffp, ed, and ef considerations for sputtering
   logical, public :: reflect             ! initial condition reflection flag
   logical, public :: usebk               ! use backup file
   logical, public :: useion
   logical, public :: usegeom
   logical, public :: usehis              ! use histogram to create damage profile
   logical, public :: usedef              ! use deterministic defect file
   logical, public :: usecell             ! use target composed of cells

   character(80), public :: bkfile        ! backup file
   character(80), public :: cellfile      ! file containing information about target
                                          !    consisting of cells instead of regions
   character(80), public :: deffile       ! file containing deterministic defects
   character(80), public :: hisfile       ! histogram file
   character(7), public :: rng            ! name of random number generator

   integer(INT4), public :: atom1         ! species index from the histogram file (hisfile),
                                          ! which  concentrations will be used to initialize
                                          ! concentration of atom species atom2 (1,..nelem).
   integer(INT4), public :: atom2         ! atom species index (1...nelem) to which the atom1
                                          ! concentration is applied
   integer(INT4), public :: cellrec       ! record to be used for reading from file
   integer(INT4), public :: iarand        ! seed of random number generator
   integer(INT4), public :: ibrand        ! seed of random number generator
   integer(INT4), public :: irand         ! seed of random number generator
   integer(INT4), public :: mrand         ! seed of random number generator
   integer(INT4), public :: ndim          ! number of dimensions of the geometry
   integer(INT4), public :: nelem         ! number of chemical elements (=atom species) in
                                          !    simulation
   integer(INT4), public :: nelemIon      ! number of chemical elements (=atom species) in ion
                                          !    (will be set in readIon)
                                          !    elem(1:nelemIon) belong to ion
   integer(INT4), public :: nelemTarg     ! number of target atom species (will be set in
                                          !    readIon)
                                          !    ldyn=t: elem(1:nelem) belong to target
                                          !    ldyn=f: elem(nelemIon+1:nelem) belong to target
   integer(INT4), public :: nion          ! number of simulated ions
   integer(INT4), public :: nion1         ! 0 or read from backup file
   integer(INT4), public :: nionbk        ! make backup every nionbk ions
   integer(INT4), public :: nionchunk     ! number of simulated ions per chunk
   integer(INT4), public :: niongeom      ! number of ions for mesh updates
   integer(INT4), public :: nionhis       ! number of ions for damage updates
   integer(INT4), public :: nionout       ! number of ions for cell outputs
   integer(INT4), public :: niontra       ! number of ions for trajectory output
   integer(INT4), public :: nr            ! number of regions
   integer(INT4), public :: nthread       ! number of threads
   integer(INT4), public :: randskip      ! random numbers to be initially skipped

   real(REAL8), public :: shiftto         ! position to shift the surface of a CELL or HIS
                                          !    file to
   real(REAL8), public :: temp            ! temperature of the target
!
!  Private Data
!
   logical, private :: equicell           ! equidistant cell grid flag (obsolete)
!
!=========================================================================================
!
contains
!
!==========================================================================================
!
   subroutine readSetup( )
!             =========
!
!     Read the &setup record from the input file and allocate
!     memory for atoms
!
!  Arguments:
!
!     None
!
!  Namelist:
!
      character(80) :: bkfil                    ! obsolete
      character(80) :: filcell                  ! obsolete
      character(80) :: fildef                   ! obsolete
      character(80) :: filhis                   ! obsolete
      integer(INT4) :: natom

      namelist /setup/ ldamdyn, ldamdyn_set, ldyn, lsput, reflect, usecell, usebk,  &
                       usedef, usehis, useion,  &  ! usegeom
                       bkfile, cellfile, deffile, hisfile, rng,  &
                       atom1, atom2, cellrec, iarand, ibrand, irand, mrand,  &
                       natom, ndim, nion, nionbk, nionchunk, niongeom, nionhis, nionout,  &
                       niontra, nr, nthread,  &
                       randskip, shiftto, temp
      ! obsolete
      namelist /setup/ bkfil, equicell, filcell, fildef, filhis
!
!  Local Variables:
!
      character(256) :: msg
      integer(INT4) :: ios
!
!  Start of Executable Code:
!  ------------------------
!
!  Determine whether ldamdyn is set in the &setup record
!
      ldamdyn_set = .false.

      call openFile( 11, 'INP', 'formatted', 'read' )
!          ========  i     i         i         i

      ldamdyn = .false.
      read( 11, nml=setup, iostat=ios, iomsg=msg )
      if( ios /= 0 ) then
         print '(/a)',  trim( msg )
         if( ios < 0 )  goto 900
         if( ios > 0 )  goto 901
      end if
      do
         read( 11, nml=setup, iostat=ios, iomsg=msg )
         if( ios < 0 )  exit
         if( ios > 0 ) then
            print '(/a)',  trim( msg )
            goto 901
         end if
      end do
      if( ldamdyn )  ldamdyn_set = .true.

      rewind 11
      ldamdyn = .true.
      do
         read( 11, nml=setup, iostat=ios, iomsg=msg )
         if( ios < 0 )  exit
         if( ios > 0 ) then
            print '(/a)',  trim( msg )
            goto 901
         end if
      end do
      if( .not.ldamdyn )  ldamdyn_set = .true.
!
!  Default values
!
      lsput = .false.
      ldamdyn = .false.
      ldyn = .false.
      usebk = .false.
      useion = .true.
!      usegeom = .true.       ! deactivated due to bug
      usegeom = .false.
      usehis = .false.
      usedef = .false.
      usecell = .false.
      equicell = .false.                        ! obsolete
      reflect = .false.

      call undefine( bkfil )                    ! obsolete
      bkfile = 'BCK'
      call undefine( cellfile )
      call undefine( deffile )
      call undefine( filcell )                  ! obsolete
      call undefine( fildef )                   ! obsolete
      call undefine( filhis )                   ! obsolete
      call undefine( hisfile )
      rng = 'Knuth'

      call undefine( atom1 )
      call undefine( atom2 )
      cellrec = 0
      iarand = 15678
      ibrand = 67638
      irand = 5421
      mrand = 6103
      call undefine( natom )
      ndim = 1
      call undefine( nion )
      call undefine( nionbk )
      call undefine( nionchunk )
      call undefine( niongeom )
      call undefine( nionhis )
      call undefine( nionout )
      niontra = 1
      call undefine( nr )
      call undefine( nthread )
      randskip = 0

      call undefine( shiftto )
      temp = 300.
!
!  Read the &setup record from the input file
!
      rewind 11

      do
         read( 11, nml=setup, iostat=ios, iomsg=msg )
         if( ios < 0 )  exit
         if( ios > 0 ) then
            print '(/a)',  trim( msg )
            goto 901
         end if
      end do

      close( 11 )
!
!  Check obsolete parameters
!
      if( equicell )  goto 939
      if( defined( bkfil ) )  goto 940
      if( defined( filcell ) )  goto 942
      if( defined( fildef ) )  goto 943
      if( defined( filhis ) )  goto 941
!
!  Check the validity of the parameters
!
      if( ldyn .and. usedef )  goto 1001
      if( ldyn .and. usebk )  goto 1002
      if( usecell .and. usebk )  goto 1004
      if( usecell .and. usedef )  goto 1006
      if( usedef .and. usebk )  goto 924
      if( usehis .and. usebk )  goto 917
      if( usehis .and. usedef )  goto 923

      if( usecell .and. .not. defined( cellfile ) )  goto 1003
      if( usedef .and. .not. defined( deffile ) )  goto 922
      if( usehis .and. .not. defined( hisfile ) )  goto 916

      if( rng /= 'Haas' .and. rng /= 'Fortran' .and. rng /= 'Knuth' )  goto 927
      if( mrand < 1 .or. mrand > 8000)  goto 9281
      if( iarand < 1 .or. iarand > 70000 )  goto 9282
      if( ibrand < 1 .or. ibrand > 90000 ) goto 9283
      if( irand < 1 .or. irand > 10000 )  goto 9284

      if( cellrec < 0 )  goto 1007

      if( .not.defined( natom ) )  goto 911
      if( natom <= 0 )  goto 912

      if( ndim < 1 .or. ndim > 3 )  goto 910
      if( ndim > 1 .and. defined( shiftto ) )  goto 925
      if( ndim > 1 .and. reflect )  goto 926

      if( usecell ) then
         if( defined( nr ) ) then
            print *,  '(WARNING) &setup: nr ignored when usecell=t'
         else
            nr = 0
         end if
      else
         if( .not.defined( nr ) )  goto 913
         if( nr <= 0 )  goto 914
      end if

      if( .not.defined( nthread ) ) then
         nthread = 1
!$       nthread = omp_get_max_threads( )
      end if
      if( nthread <= 0 )  goto 9141

      if( defined( randskip ) ) then
         if( randskip < 0 )  goto 931
      end if

      if( temp < 0 )  goto 915
!
!  Check parameters that are temporarily not allowed for parallel execution
!
      if( nthread > 1 ) then
         !if( ldyn )  goto 950
         !if( usebk )  goto 951
         !if( usecell )  goto 952
         if( usedef )  goto 953
         !if( usehis )  goto 954
      end if
!
!  Internally we use nelem (number of chemical elements) instead of natom (number of atom 
!  species)
!
      nelem = natom
!
!  Check validity of atom1 and atom2
!
      if( usehis ) then
         if( ldyn ) then
            if( defined( atom1 ) .or. defined( atom2 ) )  goto 1000
         else
            if( .not. defined( atom1 ) )  goto 918
            if( .not. defined( atom2 ) )  goto 919
            if( (atom2 < 1) .or. (atom2 > nelem) )  goto 920
            if( atom1 < 1 )  goto 921
         end if
      end if

      return
!
!  Error messages
!
900   error stop '&setup: not found in input file'
901   error stop '&setup: error during reading namelist'

910   error stop '&setup: ndim=1, 2 or 3 required'
911   error stop '&setup: natom must be specified'
912   error stop '&setup: natom>0 required'
913   error stop '&setup: Either nr or usecell must be specified'
9131  error stop '&setup: nr and usecell=t must not be specified both'
914   error stop '&setup: nr>0 required if usecell=F'
9141  error stop '&setup: nthread>0 required'
915   error stop '&setup: temp>0 required'
916   error stop '&setup: hisfile required if usehis=T'
917   error stop '&setup: specifiying both usehis=T and usebk=T not allowed'
918   error stop '&setup: atom1 required if usehis=T'
919   error stop '&setup: atom2 required if usehis=T'
920   error stop '&setup: 1<=atom2<=natom required'
921   error stop '&setup: 1<=atom1 required'
922   error stop '&setup: deffile required if usedef=T'
923   error stop '&setup: specifiying both usehis=T and usedef=T not allowed'
924   error stop '&setup: specifiying both usedef=T and usebk=T not allowed'
925   error stop '&setup: shiftto must not be specified for ndim>1'
926   error stop '&setup: reflect=t must not be specified for ndim>1'
927   error stop '&setup: rng must be one of ''Haas'', ''Fortran'', ''Knuth'''
928   error stop '&setup: mrand, iarand, ibrand, irand must be positive'
9281  error stop '&setup: 0<mrand<=8000 required'
9282  error stop '&setup: 0<iarand<=70000 required'
9283  error stop '&setup: 0<ibrand<=90000 required'
9284  error stop '&setup: 0<irand<=10000 required'

931   error stop '&setup: randskip>=0 required'

939   error stop '&setup: equicell is obsolete. Use grid in &geom instead'
940   error stop '&setup: bkfil is obsolete. Use bkfile instead'
941   error stop '&setup: filhis is obsolete. Use hisfile instead'
942   error stop '&setup: filcell is obsolete. Use cellfile instead'
943   error stop '&setup: fildef is obsolete. Use deffile instead'

950   error stop '&setup: ldyn=t not allowed for parallel execution'
!951   error stop '&setup: usebk=t not allowed for parallel execution'
!952   error stop '&setup: usecell=t not allowed for parallel execution'
953   error stop '&setup: usedet=t not allowed for parallel execution'
!954   error stop '&setup: usehis=t not allowed for parallel execution'

1000  error stop '&setup: specifiying atom1 or atom2 with ldyn=T, usehis=T not allowed'
1001  error stop '&setup: specifiying both ldyn=T and usedef=T not allowed'
1002  error stop '&setup: specifiying both ldyn=T and usebk=T not allowed'
1003  error stop '&setup: cellfile required if usecell=T'
1004  error stop '&setup: specifiying both usecell=T and usebk=T not allowed'
1006  error stop '&setup: specifiying both usecell=T and usedef=T not allowed'
1007  error stop '&setup: cellrec>=0 required'

   end subroutine readSetup
!
!=======================================================================
!
   subroutine initSetup( )
!             =========
!
!     Make initializations that need data from other records.
!
!  Arguments:
!
!     None
!
!  Program History:
!
!     24-Jun-2023  G.Hobler  creation
!
!  Local Variables:
!
      real(REAL8) :: dummy
!
!  Start of Executable Code:
!  ------------------------
!
!  Check obligatory parameter nion and validity of parameters
!
      if( .not.defined( nion ) )  goto 961

      if( nion <= 0 )  goto 960

      if( defined( nionchunk ) ) then
         if( nionchunk <= 0 )  goto 963
      end if

      if( defined( nionbk ) ) then
         if( nionbk < 0 )  goto 9630
         if( nionbk > 0 .and. ldyn )  goto 9631
      end if

      if( defined( niongeom ) ) then
         if( niongeom <= 0 )  goto 927
      endif

      if( defined( nionhis ) ) then
         if( nthread > 1 ) then
            if( nionhis <= 0 .or. nionhis > nion )  goto 9745
         else
            if( nionhis < 0 .or. nionhis > nion )  goto 9746
         end if
      end if

      if( defined( nionout ) ) then
         if( nionout <= 0 )  goto 929
      endif

      if( defined( niontra ) ) then
         if( niontra <= 0 )  goto 928
      end if
!
!  Temporarily undefine nionhis and niongeom if they will not be used
!
      if( .not.ldyn .and. defined( niongeom ) ) then
         call undefine( niongeom )
         print *, '(WARNING) &setup: niongeom ignored when ldyn=f'
      end if

      if( .not.ldamdyn .and. .not.ldyn .and. defined( nionhis ) ) then
         call undefine( nionhis )
         print *, '(WARNING) &setup: nionhis ignored when ldamdyn=f and ldyn=f'
      end if
!
!  If defined, make nionchunk multiple of nthread
!
      if( defined( nionchunk ) .and. mod( nionchunk, nthread ) /= 0 )  &
            nionchunk = (nionchunk/nthread + 1) * nthread
!
!  Make nionhis multiple of nionchunk, or if undefined of nthread
!
      if( defined( nionhis ) ) then
         if( defined( nionchunk ) ) then
            if( mod( nionhis, nionchunk ) /= 0 )  &
                  nionhis = (nionhis/nionchunk + 1) * nionchunk
         else
            if( mod( nionhis, nthread ) /= 0 )  &
                  nionhis = (nionhis/nthread + 1) * nthread
         end if
      end if
!
!  Make niongeom multiple of nionhis, or if undefined of nionchunk, or if undefined of nthread
!
      if( defined( niongeom ) ) then
         if( defined( nionhis ) .and. nionhis /= 0 ) then
            if( mod( niongeom, nionhis ) /= 0 )  &
                  niongeom = (niongeom/nionhis + 1) * nionhis
         else if( defined( nionchunk ) ) then
            if( mod( niongeom, nionchunk ) /= 0 )  &
                  niongeom = (niongeom/nionchunk + 1) * nionchunk
         else
            if( mod( niongeom, nthread ) /= 0 )  &
                  niongeom = (niongeom/nthread + 1) * nthread
         end if
      end if
!
!  Make nion multiple of niongeom, or if undefined of nionhis, or if undefined of nionchunk,
!  or if undefined of nthread
!
      if( defined( niongeom ) ) then
         if( mod( nion, niongeom ) /= 0 )  &
               nion = (nion/niongeom + 1) * niongeom
      else if( defined( nionhis ) .and. nionhis /= 0 ) then
         if( mod( nion, nionhis ) /= 0 )  &
               nion = (nion/nionhis + 1) * nionhis
      else if( defined( nionchunk ) ) then
         if( mod( nion, nionchunk ) /= 0 )  &
               nion = (nion/nionchunk + 1) * nionchunk
      else
         if( mod( nion, nthread ) /= 0 )  &
               nion = (nion/nthread + 1) * nthread
      end if
!
!  Make niontra multiple of nionchunk if both are defined
!
      if( defined( niontra ) .and. defined( nionchunk ) ) then
         if( mod( niontra, nionchunk ) /= 0 )  &
               niontra = (niontra/nionchunk + 1) * nionchunk
      end if
!
!  Set default value of niongeom
!
      if( .not.defined( niongeom ) ) then
         if( ldyn ) then
            if( defined( nionhis ) .and. nionhis /= 0 ) then
               niongeom = nionhis
            else if( defined( nionchunk ) ) then
               niongeom = nionchunk
            else
               niongeom = nthread
            end if
         else
            niongeom = 0                           ! make niongeom==0 alias of ldyn=f
         end if
      end if
!
!  Set default value of nionhis
!
      if( .not.defined( nionhis ) ) then
         if( ldamdyn .or. ldyn ) then
            if( nthread == 1 ) then
               nionhis = 0
            else
               if( defined( nionchunk ) ) then
                  nionhis = nionchunk
               else
                  nionhis = nthread
               end if
            end if
         else
            nionhis = -1                           ! make nionhis==-1 alias of ldamdyn=f and
         end if                                    !    ldyn=f
      end if
!
!  Set default value of nionchunk
!
      if( .not.defined( nionchunk ) ) then
         if( nionhis > 0 ) then
            nionchunk = nionhis
         else if( niongeom > 0 ) then
            nionchunk = niongeom
         else
            nionchunk = nion
         end if
         if( defined( niontra ) ) then
            nionchunk = ggT( nionchunk, niontra )
         end if
      end if
!
!  Make nionout multiple of niongeom
!
      if( niongeom > 0 .and. defined( nionout ) ) then
         if( mod( nionout, niongeom ) /= 0 )  &
               nionout = (nionout/niongeom + 1) * niongeom
      else
         nionout = niongeom
      end if
!
!  Set default value of nionbk
!
      if( .not.defined( nionbk ) )  nionbk = 0
!
!  If nionbk is defined, make it multiple of nthread and adjust nionchunk so that also nionbk
!  is a multiple of nionchunk
!
      if( nionbk /= 0 ) then
         if( mod( nionbk, nionchunk ) /= 0 )  &
               nionbk = (nionbk/nionchunk + 1) * nionchunk
         if( nionhis > 0 ) then
            if( mod( nionbk, nionhis ) /= 0 )  &
                  print *, '(WARNING) &setup: nionbk is not a multiple of nionhis'
         end if
      end if
!
!  Calculate derived parameters
!
      if( usebk .and. useion ) then
         call openFile( 22, bkfile, 'unformatted', 'read' )
!             ========  i     i           i          i
         !read( 22, err=930, end=930 ) nion1, dummy
         read( 22 ) nion1, dummy
         if( nion < nion1 )  &!goto 920
               error stop '&setup: nion must not be smaller than nion on backup file'
         if( mod( nion1, nthread ) /= 0 )  &!goto 921
               error stop '&setup: nion on backup file must be multiple of nthread'
         rewind 22
         close( 22 )
      else
         nion1 = 0
      end if

      return

!920      error stop '&setup: nion must not be smaller than nion on backup file'
!921      error stop '&setup: nion on backup file must be multiple of nthread'
!930      error stop '&setup: error during reading backup file'

927   error stop '&setup: niongeom>0 required'
928   error stop '&setup: niontra>0 required'
929   error stop '&setup: nionout>0 required'

960   error stop '&setup: nion>0 required'
961   error stop '&setup: nion must be specified'
963   error stop '&setup: nionchunk>0 required'

9630  error stop '&output: nionbk>0 or =0 required'
9631  error stop '&output: nionbk requires ldyn=f'

9745  error stop '&setup: 0<nionhis<=nion required if nthread>1'
9746  error stop '&setup: 0<=nionhis<=nion required if nthread=1'

   end subroutine initSetup
!
!=========================================================================================
!
   subroutine writeSetup( )
!             ==========
!
!     Write information derived from the &setup record to the output file
!
!  Arguments:
!
!     None
!
!  Local Variables:
!
!      integer :: m, ia, ib, irand
!
!  Start of Executable Code:
!  ------------------------
!
!  Open the output file
!
      call openFile( 21, 'OUT', 'formatted', 'append' )
!          ========  i     i         i          i

      write( 21, '(t2,a, t13,a)' )  &
            '&setup:',  &
            'ndim=' // str( ndim ) //  &
            ', nr=' // str( nr ) //  &
            ', natom=' // str( nelem ) //  &
            ', natomIon=' // str( nelemIon ) //  &
            ', natomTarg=' // str( nelemTarg )

      if( ldamdyn ) then
         write( 21, '(t13,a)' )  &
               'ldamdyn=t: Damage is dynamically taken into account'
      else
         write( 21, '(t13,a)' )  &
               'ldamdyn=f: Damage is not taken into account dynamically'
      end if
      if( ldyn ) then
         write( 21, '(t13,a)' )  &
               'ldyn=t: Target composition change is dynamically taken into account'
      else
         write( 21, '(t13,a)' )  &
               'ldyn=f: Target composition changes are not taken into account'
      end if
      if( lsput ) then
         write( 21, '(t13,a)' )  &
               'lsput=t: Enhanced-accuracy mode for sputtering simulation enabled'
      else
         write( 21, '(t13,a)' )  &
               'lsput=f: No enhanced accuracy for sputtering simulation'
      end if            
      if( usedef )  &
         write( 21, '(t13,a)' )  &
               'Deterministic defect file used: ' // trim( deffile )
      if( usehis ) then
         write( 21, '(t13,a)' )  &
               'Histogram file used: ' // trim( hisfile )
         write( 21, '(t13,a)' )  &
               'Atom ' // str( atom2 ) //  &
               ' will be initialized with atom ' // str( atom1 ) //  &
               ' from histogram file'
      end if
      if( usecell ) then
         write( 21, '(t13,a)' )  &
               'usecell=t'
         if( cellrec == 0 ) then
            write( 21, '(t13,a)' )  &
                  'Read last record from cell file: ' // trim( cellfile )
         else
            write( 21, '(t13,a)' )  &
                  'Read record' // str(cellrec) //  &
                  ' from cell file: '
            write( 21, '(t13,a)' )  &
                  '   ' // trim( cellfile )
         end if
      else
         write( 21, '(t13,a)' )  &
               'usecell=f: Geometry and materials defined by regions'
      end if
      if( usehis .or. usecell ) then
         if( defined( shiftto ) )  &
            write( 21, '(t13,a)' )  &
                  'Histogram/cell distribution is shifted so that surface is at z=' //  &
                  str( shiftto, 'f99.1' )
         if( reflect )  &
            write( 21, '(t13,a)' )  &
                  'Histogram/cell distribution is reflected at z=0'
      end if
      write( 21, '(t13,a)' )  &
            'Number of ions to be simulated: nion=' // str( nion )
      if( ldamdyn ) then
         if( nionhis == 0 ) then
            write( 21, '(t13,a)' )  &
                  'Damage histograms are updated immediately (nionhis=0)'
         else
            write( 21, '(t13,a)' )  &
                  'Damage histograms are updated every nionhis=' // str( nionhis, 'i99' ) //  &
                  ' ions'
         end if
      end if
      if( ldyn ) then
         write( 21, '(t13,a)' )  &
               'Cell contents are updated every nionhis=' //  &
               str( nionhis, 'i99' ) // ' ions'
         write( 21, '(t13,a)' )  &
               'Cells are checked for geometry updates every niongeom=' //  &
               str( niongeom, 'i99' ) // ' ions'
         write( 21, '(t13,a)' )  &
               'Cell output is written after every nionout=' //  &
               str( nionout, 'i99' ) // ' ions'
      end if
      if( nionbk>0 ) then
         write( 21, '(t13,a)' )  &
               'Backup file is requested after every nionbk=' // str( nionbk, 'i99' ) // ' ions'
      end if
      write( 21, '(t13,a)' )  &
            'Chunk size: nionchunk=' // str( nionchunk )
      write( 21, '(t13,a)' )  &
            'Number of threads: nthread=' // str( nthread )

      write( 21, '(/t13,a,a)' )  &
            'Random number generator: ', rng
!      call ranget( m, ia, ib, irand )
!!          ======  o  o   o     o
!      write( 21, '(t13,a,i4,a,i5,a,i6,a,i4)' )  &
!            'Initial state: m=' // str( m ) //  &
!            ', ia=' // str( ia ) //  &
!            ', ib=' // str( ib ) //  &
!            ', irand=' // str( irand )

      write( 21, '()' )

      close( 21 )

   end subroutine writeSetup
!
!=========================================================================================
!
end module Param_Setup_mod
