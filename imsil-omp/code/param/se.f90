module Param_Se_mod
!
!  Hold data about the electronic stopping power and manage their input and output.
!
!  Used Modules:
!  ============
!
   use Const_mod,          only: INT4, REAL8,  &
                                 ZERO, ONE, TWO, FIVE, HALF, THIRD, FIVE
   use Undefine_mod,       only: undefine, defined, definedlength
   use Stack_mod,          only: sizeOfStack, Real8Stack,  &
                                 newStack, sizeOfStack, popStack, pushStack, clearStack
   use Open_mod,           only: openFile
   use Param_WriteTab_mod, only: writeTab
   use String_mod,         only: to_lower

   use Param_Setup_mod,    only: nelem, usecell
   use Param_Elem_mod,     only: elem
   use Param_Ion_mod,      only: ion
   use Param_Mat_mod,      only: mat
   use Param_Dam_mod,      only: dam
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
   public :: readSe                 ! read &separ record from the input file and initialize
!            ======                 !    further parameters
   public :: initSe                 ! initialize the se structure
!            ======
   public :: writeSe                ! write data derived from &separ record to output file
!            =======
!
!  Private Procedures:
!  ==================
!
   private :: defaultSeLind         ! provide default values of Lindhard stopping power
!             =============         !    parameters
   private :: defaultXnl            ! provide default values of OR model parameters
!             ==========            !
!
!  Public Data:
!  ===========
!
   type :: Se_type
!          =======

      logical :: lapsis                         ! flag indicating whether the apsis or the
                                                !    impact parameter be used for electronic
                                                !    stopping evaluation
      logical, allocatable :: lfac2scr(:)       ! (nelem) flag indicating that fac2scr
                                                !    differs from facscr
      logical :: lse                            ! flag indicating whether electronic stopping
                                                !    shall be considered at all
      logical, allocatable :: lstrag(:,:)       ! (nelem,nelem) flag indicating electronic
                                                !    straggling to be taken into account
      logical :: ltab                           ! flag indicating that a table of the stopping
                                                !    power as a function of energy be generated

      character(6), allocatable :: locmodel(:,:)! (nelem,nelem) local electronic stopping model
      character(80), allocatable :: sefile(:,:) ! (nelem,nelem) electronic stopping file
      character(8), allocatable :: semodel(:,:) ! (nelem,nelem) electronic stopping model
      character(2), allocatable :: xnlmodel(:,:)! (nelem,nelem) model for the energy dependence
                                                !    of the nonlocal fraction

      integer(INT4), allocatable :: nenl(:,:)   ! (nelem,nelem) number of enl values
      integer(INT4), allocatable :: nfacscr(:,:)! (nelem,nelem) number of facscr values
      integer(INT4), allocatable :: nxnl(:,:)   ! (nelem,nelem) number of xnl values

      real(REAL8), allocatable :: c0bethe(:,:)  ! (nelem,nelem) c0 in Bethe-Bloch formula
      real(REAL8), allocatable :: c1bethe(:,:)  ! (nelem,nelem) c1 in Bethe-Bloch formula
      real(REAL8), allocatable :: corlin(:,:)   ! (nelem,nelem) correction to Lindhard stopping
      real(REAL8), allocatable :: eminse(:,:)   ! (nelem,nelem) energy below which se vanishes.
                                                !    e-eminse instead of e is passed to the
                                                !    subroutines to calculate se
      real(REAL8), allocatable :: enl(:,:,:)    ! (nenlm,nelem,nelem) energies for xnl
      real(REAL8), allocatable :: estrag(:,:)   ! (nelem,nelem) transition energy from linear
                                                !    regime to constant Bohr straggling
      real(REAL8), allocatable :: facscr(:,:,:) ! (nfacscrm,nelem,nelem) factor for screening
                                                !    length
      real(REAL8), allocatable :: fac2scr(:,:)  ! (nelem,nelem) high-energy value of facscr
      real(REAL8), allocatable :: powfirs(:,:)  ! (nelem,nelem) power in Firsov model
      real(REAL8), allocatable :: powint(:,:)   ! (nelem,nelem) power for high-energy
                                                !    interpolation
      real(REAL8), allocatable :: powlin(:)     ! (nelem) power for Lindhard stopping
      real(REAL8) :: ptab                       ! impact parameter [A] to be used for table
                                                !    generation
      real(REAL8), allocatable :: xnl(:,:,:)    ! (nxnlm,nelem,nelem) nonlocal fractions

   end type Se_type

   type(Se_type), public :: se
!
!=======================================================================
!
contains
!
!=======================================================================
!
   subroutine readSe( )
!             ======
!
!     Read the &separ record from the input file
!
!  Arguments:
!
!     None
!
!  Namelist variables:
                                                ! Index parameters
      integer(INT4) :: atom1
      integer(INT4) :: atom2

      logical :: lapsis
      logical :: lse                            ! obsolete
      logical :: lse_save                       ! obsolete
      logical :: ltab

      character(80)  :: filse                   ! obsolete
      character(6) :: locmodel
      character(80)  :: sefile                  ! electronic stopping file
      character(3)   :: firsov
      character(8)   :: semodel
      character(3)   :: straggle
      character(2)   :: xnlmodel

      integer(INT4) :: model                    ! obsolete

      real(REAL8) :: c0bethe
      real(REAL8) :: c1bethe
      real(REAL8) :: corlin
      real(REAL8) :: eminse
      real(REAL8) :: enl(21)
      real(REAL8) :: estrag
      real(REAL8) :: facscr(21)
      real(REAL8) :: fac2scr
      real(REAL8) :: powfirs
      real(REAL8) :: powint
      real(REAL8) :: powlin
      real(REAL8) :: ptab
      real(REAL8) :: xnl(21)

      namelist /separ/  atom1, atom2,  &
                        lapsis, lse, ltab,  &
                        firsov, locmodel, sefile, semodel, straggle, xnlmodel,  &
                        model,  &
                        c0bethe, c1bethe, corlin, eminse, enl, estrag, facscr, fac2scr,  &
                        powfirs, powint, powlin, ptab, xnl
      ! obsolete
      namelist /separ/  filse
!
!  Other local variables:
!
      character(256) :: msg
      logical :: isEmpty
      integer(INT4) :: ielem1, ielem1min, ielem1max
      integer(INT4) :: ielem2, ielem2min, ielem2max
      integer(INT4) :: i, nenl, nfacscr, nxnl, nenlm, nfacscrm, nxnlm, ios
      type(Real8Stack), allocatable :: enlStack(:,:), facscrStack(:,:), xnlStack(:,:)
!
!  Start of Executable Code:
!  ------------------------
!
!  Allocate memory
!                                               ! Parameters defining the model
      allocate( se%lfac2scr( nelem ) )
      allocate( se%lstrag( nelem, nelem ) )

      allocate( se%locmodel( nelem, nelem ) )
      allocate( se%sefile( nelem, nelem ) )
      allocate( se%semodel( nelem, nelem ) )
      allocate( se%xnlmodel( nelem, nelem ) )

      allocate( se%nenl( nelem, nelem ) )
      allocate( se%nfacscr( nelem, nelem ) )
      allocate( se%nxnl( nelem, nelem ) )

      allocate( se%c0bethe( nelem, nelem ) )
      allocate( se%c1bethe( nelem, nelem ) )
      allocate( se%corlin( nelem, nelem ) )
      allocate( se%eminse( nelem, nelem ) )
      allocate( se%estrag( nelem, nelem ) )
      allocate( se%fac2scr( nelem, nelem ) )
      allocate( se%powfirs( nelem, nelem ) )
      allocate( se%powint( nelem, nelem ) )
      allocate( se%powlin( nelem ) )
!
!  Set default values
!
      if( dam%lrcoil .and. dam%elimkp<100 ) then
         se%lapsis = .true.
      else
         se%lapsis = .false.
      end if
      se%lfac2scr(:) = .false.
      lse_save = .true.                         ! obsolete
      se%lstrag(:,:) = .false.                  ! Electronic stopping parameters
      se%ltab = .false.

      call undefine( se%locmodel(:,:) )
      call undefine( se%sefile(:,:) )
      call undefine( se%semodel(:,:) )
      call undefine( se%xnlmodel(:,:) )

      call undefine( se%nenl(:,:) )
      call undefine( se%nfacscr(:,:) )
      call undefine( se%nxnl(:,:) )

      call undefine( se%c0bethe(:,:) )
      call undefine( se%c1bethe(:,:) )
      call undefine( se%corlin(:,:) )
      call undefine( se%eminse(:,:) )
      call undefine( se%estrag(:,:) )
      call undefine( se%fac2scr(:,:) )
      call undefine( se%powfirs(:,:) )
      call undefine( se%powint(:,:) )
      call undefine( se%powlin(:) )
      call undefine( se%ptab )
!
!  Set up stacks for enl, facscr, and xnl
!
      allocate( enlStack( nelem, nelem ) )
      allocate( facscrStack( nelem, nelem ) )
      allocate( xnlStack( nelem, nelem ) )
      do ielem1 = 1, nelem
         do ielem2 = 1, nelem
            call newStack( enlStack(ielem2,ielem1) )
!                ========            o
            call newStack( facscrStack(ielem2,ielem1) )
!                ========             o
            call newStack( xnlStack(ielem2,ielem1) )
!                ========            o
         end do
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
         call undefine( atom1 )
         call undefine( atom2 )

         lapsis = se%lapsis
         lse = lse_save
         ltab = se%ltab

         call undefine( filse )                 ! obsolete
         call undefine( firsov )
         call undefine( locmodel )
         call undefine( sefile )
         call undefine( semodel )
         call undefine( straggle )
         call undefine( xnlmodel )

         call undefine( model )

         call undefine( c0bethe )
         call undefine( c1bethe )
         call undefine( corlin )
         call undefine( eminse )
         call undefine( enl(:) )
         call undefine( estrag )
         call undefine( facscr(:) )
         call undefine( fac2scr )
         call undefine( powfirs )
         call undefine( powint )
         call undefine( powlin )
         call undefine( ptab )
         call undefine( xnl(:) )
!
!  Read &separ record from the input file
!
         read( 11, nml=separ, iostat=ios, iomsg=msg )
         if( ios < 0 )  exit
         if( ios > 0 ) then
            print '(/a)',  trim( msg )
            goto 900
         end if
!
!  Check obsolete parameters
!
         if( defined( filse ) )  goto 910
!
!  Check index parameters
!
         if( defined( atom1 ) ) then
            if( atom1 < 1 .or. atom1 > nelem )  goto 906
         end if

         if( defined( atom2 ) ) then
            if( atom2 < 1 .or. atom2 > nelem )  goto 907
         end if
!
!  Check parameters on consistency with index parameters
!
         if( defined( atom1 ) ) then
            if( defined( ptab ) )  goto 924
         end if

         if( defined( atom2 ) ) then
            if( defined( ptab ) )  goto 933
            if( defined( powlin ) )  goto 935
         end if
!
!  Define range of atoms to which parameters apply
!
         if( defined( atom1 ) ) then
!            =======
            ielem1min = atom1
            ielem1max = atom1
         else
            ielem1min = 1
            ielem1max = nelem
         end if
         ielem1 = ielem1min

         if( defined( atom2 ) ) then
!            =======
            ielem2min = atom2
            ielem2max = atom2
         else
            ielem2min = 1
            ielem2max = nelem
         endif
         ielem2 = ielem2min
!
!  Check other consistencies
!
         if( defined( firsov ) .and. defined( locmodel ) )  goto 9521
         if( (lse.neqv.lse_save) .and. defined( semodel ) )  goto 952
!
!  Automatic setting of parameters
!
         if( defined( estrag ) .and. .not.defined( straggle ) )  straggle = 'on'
!
!  Copy parameters to se structure
!
         se%lapsis = lapsis
         select case( firsov )
            case( 'on', 'On', 'ON' )
               se%locmodel( ielem2min:ielem2max, ielem1min:ielem1max ) = 'firsov'
            case( 'off', 'Off', 'OFF' )
            case default
               if( defined( firsov ) )  goto 951
         end select
         select case( straggle )
               case( 'on', 'On', 'ON' )
                  se%lstrag( ielem2min:ielem2max, ielem1min:ielem1max ) = .true.
               case( 'off', 'Off', 'OFF' )
                  se%lstrag( ielem2min:ielem2max, ielem1min:ielem1max ) = .false.
               case default
                  if( defined( straggle ) )  goto 953
         end select
         se%ltab = ltab

         do ielem1 = ielem1min, ielem1max       ! Annoyingly, we need this explicit loop because
            do ielem2 = ielem2min, ielem2max    !    of an ifx compiler error
               if( defined( locmodel ) )  &
                     se%locmodel(ielem2,ielem1) = to_lower( locmodel )
      !                                           ========     i
               if( defined( sefile ) )  &
                     se%sefile(ielem2,ielem1) = sefile
               if( defined( semodel ) )  &
                     se%semodel(ielem2,ielem1) = to_lower( semodel )
      !                                          ========     i
               if( defined( xnlmodel ) )  &
                     se%xnlmodel(ielem2,ielem1) = to_lower( xnlmodel )
      !                                           ========     i
            end do
         end do

         if( lse_save .and. .not.lse )  &
               se%semodel( ielem2min:ielem2max, ielem1min:ielem1max ) = 'none'
         if( defined( model ) ) then
            select case( model )
               case( 1, 3 )
                  se%semodel( ielem2min:ielem2max, ielem1min:ielem1max ) = 'lindhard'
               case( 2, 4 )
                  se%semodel( ielem2min:ielem2max, ielem1min:ielem1max ) = 'zbl'
               case( 5, 6 )
                  se%semodel( ielem2min:ielem2max, ielem1min:ielem1max ) = 'kkk'
               case( 7 )
                  se%semodel( ielem2min:ielem2max, ielem1min:ielem1max ) = 'table'
               case( 8 )
                  se%semodel( ielem2min:ielem2max, ielem1min:ielem1max ) = 'srim'
               case default
                  goto 950
            end select
            select case( model )
               case( 1, 2, 5 )
                  se%xnlmodel( ielem2min:ielem2max, ielem1min:ielem1max ) = 'e'
               case( 3, 4, 6, 7, 8 )
                  se%xnlmodel( ielem2min:ielem2max, ielem1min:ielem1max ) = 'se'
            end select
         end if

         if( defined( c0bethe ) )  &
               se%c0bethe( ielem2min:ielem2max, ielem1min:ielem1max ) = c0bethe
         if( defined( c1bethe ) )  &
               se%c1bethe( ielem2min:ielem2max, ielem1min:ielem1max ) = c1bethe
         if( defined( corlin ) )  &
               se%corlin( ielem2min:ielem2max, ielem1min:ielem1max ) = corlin
         if( defined( eminse ) )  &
               se%eminse( ielem2min:ielem2max, ielem1min:ielem1max ) = eminse

         if( defined( estrag ) )  &
               se%estrag( ielem2min:ielem2max, ielem1min:ielem1max ) = estrag
         if( defined( powfirs ) )  &
               se%powfirs( ielem2min:ielem2max, ielem1min:ielem1max ) = powfirs
         if( defined( powint ) )  &
               se%powint( ielem2min:ielem2max, ielem1min:ielem1max ) = powint
         if( defined( powlin ) )  &
               se%powlin( ielem1min:ielem1max ) = powlin
         if( defined( ptab ) )  se%ptab = ptab

         if( defined( enl(size( enl )) ) )  goto 941
         nenl = definedlength( enl(:) )
         if( nenl > 0 ) then
            do ielem1 = ielem1min, ielem1max
               do ielem2 = ielem2min, ielem2max
                  call clearStack( enlStack(ielem2,ielem1) )
!                      ==========           io
                  do i = 1, nenl
                     call pushStack( enl(i), enlStack(ielem2,ielem1) )
!                         =========    i              io
                  end do
               end do
            end do
         end if

         if( defined( facscr(size( facscr )) ) )  goto 942
         nfacscr = definedlength( facscr(:) )
         if( nfacscr > 0 ) then
            do ielem1 = ielem1min, ielem1max
               do ielem2 = ielem2min, ielem2max
                  call clearStack( facscrStack(ielem2,ielem1) )
!                      ==========            io
                  do i = 1, nfacscr
                     call pushStack( facscr(i), facscrStack(ielem2,ielem1) )
!                         =========    i              io
                  end do
               end do
            end do
         end if

         if( defined( xnl(size( xnl )) ) )  goto 943
         nxnl = definedlength( xnl(:) )
         if( nxnl > 0 ) then
            do ielem1 = ielem1min, ielem1max
               do ielem2 = ielem2min, ielem2max
                  call clearStack( xnlStack(ielem2,ielem1) )
!                      ==========           io
                  do i = 1, nxnl
                     call pushStack( xnl(i), xnlStack(ielem2,ielem1) )
!                         =========    i              io
                  end do
               end do
            end do
         end if

         if( defined( fac2scr ) ) then
            se%fac2scr( ielem2min:ielem2max, ielem1min:ielem1max ) = fac2scr
            se%lfac2scr( ielem1min:ielem1max ) = .true.
         end if

      end do
!
!  End of loop over reading records
!  --------------------------------
!
      close( 11 )
!
!  Copy stacks to arrays
!
      nenlm = 3   ! must be at least 3 to reserve space for default values
      nfacscrm = 1
      nxnlm = 3   ! must be at least 3 to reserve space for default values
      do ielem1 = 1, nelem
         do ielem2 = 1, nelem
            se%nenl(ielem2,ielem1) = sizeOfStack( enlStack(ielem2,ielem1) )
!                                    ===========             i
            se%nfacscr(ielem2,ielem1) = sizeOfStack( facscrStack(ielem2,ielem1) )
!                                       ===========             i
            se%nxnl(ielem2,ielem1) = sizeOfStack( xnlStack(ielem2,ielem1) )
!                                    ===========             i
            nenlm = max( nenlm, se%nenl(ielem2,ielem1) )
            nfacscrm = max( nfacscrm, se%nfacscr(ielem2,ielem1) )
            nxnlm = max( nxnlm, se%nxnl(ielem2,ielem1) )
         end do
      end do

      allocate( se%enl( nenlm, nelem, nelem) )
      allocate( se%facscr( nfacscrm, nelem, nelem) )
      allocate( se%xnl( nxnlm, nelem, nelem) )
      call undefine( se%enl(:,:,:) )
      call undefine( se%facscr(:,:,:) )
      call undefine( se%xnl(:,:,:) )

      do ielem1 = 1, nelem
         do ielem2 = 1, nelem
            do i = se%nenl(ielem2,ielem1), 1, -1
               call popStack( enlStack(ielem2,ielem1), se%enl(i,ielem2,ielem1), isEmpty )
!                   ========           io                         o                o
            end do
            do i = se%nfacscr(ielem2,ielem1), 1, -1
               call popStack( facscrStack(ielem2,ielem1), se%facscr(i,ielem2,ielem1),  &
!                   ========            io                           o
                              isEmpty )
!                                o
            end do
            do i = se%nxnl(ielem2,ielem1), 1, -1
               call popStack( xnlStack(ielem2,ielem1), se%xnl(i,ielem2,ielem1), isEmpty )
!                   ========           io                          o               o
            end do
         end do
      end do

      return
!
!  Error messages
!
900   error stop '&separ: error during reading namelist'

906   error stop '&separ: 1<=atom1<=natom required'
907   error stop '&separ: 1<=atom2<=natom required'

910   error stop '&separ: filse is obsolete. Use sefile instead'

924   error stop '&separ: atom1 and ptab must not be on same record'

933   error stop '&separ: atom2 and ptab must not be on same record'
935   error stop '&separ: atom2 and powlin must not be on same record'

941   error stop '&separ: Too many values of enl'
942   error stop '&separ: Too many values of facscr'
943   error stop '&separ: Too many values of xnl'

950   error stop '&separ: semodel=''Lindhard'',''ZBL'',''KKK'',''table'', ''SRIM'', or ''none'' required'
951   error stop '&separ: invalid value of firsov'
952   error stop '&separ: lse and semodel must not appear on same record'
9521  error stop '&separ: firsov and locmodel must not appear on same record'
953   error stop '&separ: invalid value of straggle'

970   error stop '&separ: error during reading enl or xnl'
971   error stop '&separ: internal error. Error number = 971'
974   error stop '&separ: internal error. Error number = 974'

980   error stop '&separ: sestep>0 required'
982   error stop '&separ: rs0>0 required'

   end subroutine readSe
!
!=======================================================================
!
   subroutine initSe( )
!             ======
!
!     Initializes the se structure
!
!  Arguments
!
!     None
!
!  Local variables
!
      integer(INT4) :: ielem1, ielem2, i
      real(REAL8) :: corlin, powlin, powint
      real(real8), allocatable :: enl(:), xnl(:)
!
!  Start of Executable Code:
!  ------------------------
!
!  Set default values
!
      do, ielem1 = 1, nelem
         do, ielem2 = 1, nelem

            ! determine whether there are default parameters for the Lindhard model
            call undefine( corlin )
            call undefine( powlin )
            call undefine( powint )
            call defaultSeLind( ielem2, ielem1, corlin, powlin, powint )
!                =============    i       i       o       o       o
!
!  Default Se model and Se parameters
!
            if( .not.defined( se%semodel(ielem2,ielem1) ) ) then
               if( defined( corlin ) ) then
                  se%semodel(ielem2,ielem1) = 'lindhard'
               else
                  se%semodel(ielem2,ielem1) = 'srim'
               end if
            end if

            if( .not.defined( corlin ) )  corlin = ONE
            if( .not.defined( powlin ) )  powlin = HALF
            if( .not.defined( powint ) )  powint = ONE
            if( .not.defined( se%corlin(ielem2,ielem1) ) )  se%corlin(ielem2,ielem1) = corlin
            if( .not.defined( se%powlin(ielem1) ) )  se%powlin(ielem1) = powlin
            if( .not.defined( se%powint(ielem2,ielem1) ) )  se%powint(ielem2,ielem1) = powint

            if( .not.defined( se%c0bethe(ielem2,ielem1) ) ) then
               if( elem(ielem1)%z == 1 ) then
                  se%c0bethe(ielem2,ielem1) = 1.43    ! fit to PSTAR
!                  se%c0bethe(ielem2,ielem1) = 1.7     ! fit to own+PSTAR
               else
                  se%c0bethe(ielem2,ielem1) = ONE
               end if
            end if

            if( .not.defined( se%c1bethe(ielem2,ielem1) ) ) then
               if( elem(ielem1)%z == 2 ) then
                  se%c1bethe(ielem2,ielem1) =  &
                        100. * real(elem(ielem1)%z) / real(elem(ielem2)%z)
               else if( elem(ielem1)%z == 1 ) then    ! H randSe fitted to PSTAR
                  se%c1bethe(ielem2,ielem1) =  &
                        40. * real(elem(ielem1)%z) / real(elem(ielem2)%z)
                  ! = 2.857  fit to PSTAR(exp.)
                  ! mat(ir)%prob(ielemr) * ( 53. * elem(ielem1)%z / elem(ielem2)%z )
                  ! = 3.785  fit to own+PSTAR
               else
                  se%c1bethe(ielem2,ielem1) = FIVE
               end if
            end if

            if( .not.defined( se%eminse(ielem2,ielem1) ) )  se%eminse(ielem2,ielem1) = ZERO

            if( .not.defined( se%estrag(ielem2,ielem1) ) .and. se%lstrag(ielem2,ielem1) ) then
               se%estrag(ielem2,ielem1) = 1e5 * elem(ielem1)%m
                                                   ! Default values according to Konac et al.
               select case( ielem1 )               !    (NIM-B 136-138, 159, 1998)
                  case( 1 )                        ! H
                     select case( ielem2 )
                        case( 6 )                  ! C
                           se%estrag(ielem2,ielem1) = 1e6 / 14.9 * elem(ielem1)%m
                        case( 14 )                 ! Si
                           if( elem(ielem1)%m < 1.5 ) then
                              se%estrag(ielem2,ielem1) = 1e6 / 17.2 * elem(ielem1)%m
                           else
                              se%estrag(ielem2,ielem1) = 1e6 / 16.2 * elem(ielem1)%m
                           end if
                     end select
                  case( 2 )                        ! He
                     select case( ielem2 )
                        case( 6 )                  ! C
                           se%estrag(ielem2,ielem1) = 1e6 / 6.4 * elem(ielem1)%m
                        case( 14 )                 ! Si
                           if( elem(ielem1)%m < 3.5 ) then
                              se%estrag(ielem2,ielem1) = 1e6 / 10.5 * elem(ielem1)%m
                           else
                              se%estrag(ielem2,ielem1) = 1e6 / 10.0 * elem(ielem1)%m
                           end if
                     end select
               end select
            end if

!
!  Default xnl and local model and parameters
!
            if( .not.defined( se%locmodel(ielem2,ielem1) ) )  &
                  se%locmodel(ielem2,ielem1) = 'or'

            if( .not.defined( se%xnlmodel(ielem2,ielem1) ) ) then
               if( ion%e > 10000*elem(ielem1)%m ) then
                  se%xnlmodel(ielem2,ielem1) = 'se'
               else
                  se%xnlmodel(ielem2,ielem1) = 'e'
               end if
            end if

            if( se%nenl(ielem2,ielem1) == 0 ) then
               if( se%nxnl(ielem2,ielem1) == 0 ) then
                  call defaultXnl( ielem2, ielem1, se%nenl(ielem2,ielem1),  &
!                      ==========    i       i              o
                                   se%enl(:,ielem2,ielem1), se%xnl(:,ielem2,ielem1) )
!                                            o                        o
                  se%nxnl(ielem2,ielem1) = se%nenl(ielem2,ielem1)
               else                                      ! allow to specify xnl without enl
                  allocate( xnl( size( se%xnl, 1 ) ) )   ! get default value(s) of enl
                  call defaultXnl( ielem2, ielem1, se%nenl(ielem2,ielem1),  &
!                      ==========    i       i             o
                                   se%enl(:,ielem2,ielem1), xnl(:) )
!                                               o             o
                  deallocate( xnl )
               end if
            else
               if( se%nxnl(ielem2,ielem1) == 0 ) then
                  allocate( enl( size( se%enl, 1 ) ) )      ! get default value(s) of xnl
                  call defaultXnl( ielem2, ielem1, se%nxnl(ielem2,ielem1), enl(:),  &
!                      ==========    i       i             o                 o
                                   se%xnl(:,ielem2,ielem1) )
!                                               o
                  deallocate( enl )
                  if( se%nxnl(ielem2,ielem1) > 1 )  goto 970
                                                ! do not allow to specify enl without xnl
                                                ! except when there is only one default value
               else
                  continue                      ! allow to specify both enl and xnl
               end if
            end if

            if( se%nfacscr(ielem2,ielem1) == 0 ) then
               se%nfacscr(ielem2,ielem1) = 1
               se%facscr(1,ielem2,ielem1) = 1.37 / real(elem(ielem1)%z)**THIRD
            end if
            if( .not.defined( se%fac2scr(ielem2,ielem1) ) )  &
                  se%fac2scr(ielem2,ielem1) = se%facscr(1,ielem2,ielem1)

            if( .not.defined( se%powfirs(ielem2,ielem1) ) )  se%powfirs(ielem2,ielem1) = FIVE

         end do
      end do

      if( .not.defined( se%ptab ) )  se%ptab = ZERO
!
!  Check validity of the parameters
!
      if( se%ltab .and. usecell )  goto 9521

      do ielem1 = 1, nelem
         do ielem2 = 1, nelem

            if( se%locmodel(ielem2,ielem1) /= 'or' .and.  &
                se%locmodel(ielem2,ielem1) /= 'firsov' .and.  &
                se%locmodel(ielem2,ielem1) /= 'none' )  goto 949
            if( se%semodel(ielem2,ielem1) /= 'lindhard' .and.  &
                se%semodel(ielem2,ielem1) /= 'zbl' .and.  &
                se%semodel(ielem2,ielem1) /= 'kkk' .and.  &
                se%semodel(ielem2,ielem1) /= 'table' .and.  &
                se%semodel(ielem2,ielem1) /= 'srim' .and.  &
                se%semodel(ielem2,ielem1) /= 'none' )  goto 950
            if( se%xnlmodel(ielem2,ielem1) /= 'e' .and.  &
                se%xnlmodel(ielem2,ielem1) /= 'se' )  goto 951

            if( se%powfirs(ielem2,ielem1) <= TWO )  goto 963
            if( se%corlin(ielem2,ielem1) < ZERO )  goto 960

            if( se%nfacscr(ielem2,ielem1) > 1 .and.  &
                .not.defined( se%enl(se%nfacscr(ielem2,ielem1),ielem2,ielem1) ) )  goto 9739
            if( se%nxnl(ielem2,ielem1) > 1 .and.  &
                .not.defined( se%enl(se%nxnl(ielem2,ielem1),ielem2,ielem1) ) )  goto 974
            if( se%nenl(ielem2,ielem1) > 0 .and. se%enl(1,ielem2,ielem1) <= ZERO )  goto 975
            do i = 2, se%nenl(ielem2,ielem1)
               if( se%enl(i,ielem2,ielem1) <= se%enl(i-1,ielem2,ielem1) .and.  &
                   se%enl(i,ielem2,ielem1) /= ZERO )  goto 976 ! enl=0 means e in Se maximum
            end do

            if( se%nxnl(ielem2,ielem1) == 1 ) then
               if( se%xnl(1,ielem2,ielem1) < ZERO )  goto 972
            else
               if( any( se%xnl(1:se%nxnl(ielem2,ielem1),ielem2,ielem1) <= ZERO ) )  goto 973
            end if

            if( se%lstrag(ielem2,ielem1) .and. se%estrag(ielem2,ielem1) < ZERO )  goto 954

            if( any( se%facscr(1:se%nfacscr(ielem2,ielem1),ielem2,ielem1) <= ZERO ) )  goto 977
            if( se%fac2scr(ielem2,ielem1) < ZERO )  goto 978
            if( se%powint(ielem2,ielem1) <= ZERO )  goto 962

         end do
         if( se%powlin(ielem1) <= ZERO .or. se%powlin(ielem1) > ONE )  goto 961
      end do

      if( se%ptab < ZERO )  goto 952
!
!  For SRIM model, set path of SRIM file
!
      do ielem1 = 1, nelem
         do ielem2 = 1, nelem
            if( se%semodel(ielem2,ielem1) /= 'srim' )  cycle
            se%sefile(ielem2,ielem1) = 'srimtab/SRIM2013-xx.dat'
            write( se%sefile(ielem2,ielem1)(18:19), '(i2.2)' )  elem(ielem1)%z
         end do
      end do

      return

949   error stop '&separ: locmodel=''OR'',''Firsov'', or ''none'' required'
950   error stop '&separ: semodel=''Lindhard'',''ZBL'',''KKK'',''table'', ''SRIM'', or ''none'' required'
951   error stop '&separ: xnlmodel=''e'' or ''se'' required'
952   error stop '&separ: ptab>=0 required'
9521  error stop '&setab: ltab=T not allowed with usecell=t'
954   error stop '&separ: estrag>=0 required'

960   error stop '&separ: corlin>=0 required'
961   error stop '&separ: 0<powlin<=1 required'
962   error stop '&separ: powint>0 required'
963   error stop '&separ: powfirs>2 required'

970   error stop '&separ: xnl required if enl specified'
972   error stop '&separ: xnl>=0 required'
973   error stop '&separ: xnl>0 required'
9739  error stop '&separ; at least as many enl as facscr values required'
974   error stop '&separ; at least as many enl as xnl values required'
975   error stop '&separ: enl must be positive'
976   error stop '&separ: enl must increase monotonically'
977   error stop '&separ: facscr>0 required'
978   error stop '&separ: fac2scr>=0 required'

   end subroutine initSe
!
!=======================================================================
!
   subroutine writeSe( )
!             =======
!
!     Write parameters to the output file
!
!  Arguments:
!
!     None
!
!  Local Variables:
!
      integer(INT4) :: ielem1, ielem2, i, ielemmin, ielemmax, ienl, ifacscr, ixnl

      character(10), allocatable :: model(:,:)
!
!  Start of Executable Code:
!  ------------------------
!
!  Open the output file
!
      call openFile( 21, 'OUT', 'formatted', 'append' )
!          ========  i     i         i          i

      write( 21, '(/t2,a,t13)', advance='no' ) '&separ:'

      if( all( se%semodel(:,:) == 'none' ) ) then
         write( 21, '(a)' ) 'no electronic stopping considered'
         return
      end if
!
!  Report model names
!
      allocate( model( nelem, nelem ) )
      do ielem1 = 1, nelem
         do ielem2 = 1, nelem
            select case( se%semodel(ielem2,ielem1) )
               case( 'lindhard' )
                  model(ielem2,ielem1) = ' Lindhard '
               case( 'zbl' )
                  model(ielem2,ielem1) = '   ZBL    '
               case( 'kkk' )
                  model(ielem2,ielem1) = '   KKK    '
               case( 'table' )
                  model(ielem2,ielem1) = '  table   '
               case( 'srim' )
                  model(ielem2,ielem1) = '   SRIM   '
               case( 'none' )
                  model(ielem2,ielem1) = '   none   '
            end select
         end do
      end do
      call writeTab( model(:,:), 'Se model', indent_=.false. )
!          ========      i           i             i

      do ielem1 = 1, nelem
         do ielem2 = 1, nelem
            select case( se%xnlmodel(ielem2,ielem1) )
               case( 'e' )
                  model(ielem2,ielem1) = '    E     '
               case( 'se' )
                  model(ielem2,ielem1) = '    Se    '
            end select
         end do
      end do
      call writeTab( model(:,:), 'Energy dependence of nonlocal fraction' )
!          ========      i                          i

      do ielem1 = 1, nelem
         do ielem2 = 1, nelem
            select case( se%locmodel(ielem2,ielem1) )
               case( 'or' )
                  model(ielem2,ielem1) = '    OR    '
               case( 'Firsov' )
                  model(ielem2,ielem1) = '  Firsov  '
               case( 'none' )
                  model(ielem2,ielem1) = '   none   '
            end select
         end do
      end do
      call writeTab( model(:,:), 'Local electronic stopping model' )
!          ========      i                      i

      deallocate( model )
!
!  Report minimum energy for electronic stopping
!
      call writeTab( se%eminse(:,:), 'Minimum energy for electronic stopping (eminse)' )
!          ========       i                                 i
!
!  Lindhard stopping power
!
      if( any( se%semodel(:,:) == 'lindhard' ) ) then
         write( 21, '(t13,a)' ) &
              'Parameters for interpolated Lindhard and Bethe-Bloch stopping power:'

         write( 21, '()' )
         call writeTab( se%corlin(:,:), 'corlin' )
!             ========       i             i

         call writeTab( reshape(se%powlin(:), [1, size( se%powlin(:) )] ), 'powlin' )
!             ========                      i                                 i

         call writeTab( se%powint(:,:), 'powint' )
!             ========        i            i

         call writeTab( se%c0bethe(:,:), 'c0bethe' )
!             ========         i             i

         call writeTab( se%c1bethe(:,:), 'c1bethe' )
!             ========         i             i

      end if
!
!  Stopping power defined in file
!
      if( any( se%semodel(:,:) == 'table' ) ) then
         call writeTab( se%sefile(:,:), 'Stopping power defined in file' )
!             ========       i                        i
      end if
!
!  Oen Robinson model
!
      if( se%lapsis ) then
         write( 21, '(t13,a)' )  &
              'Parameters for mixed nonlocal and apsis dependent '
      else
         write( 21, '(t13,a)' )  &
              'Parameters for mixed nonlocal and impact parameter dependent '
      end if
      write( 21, '(t16,a)' )  &
           'electronic stopping:'
!
!  Report enl
!
      write( 21, '(/t13,a)' ) 'enl:'
      do i = 0, (nelem-1)/5
         ielemmin = 5*i + 1
         ielemmax = min( 5*i+5, nelem )

         write( 21, '(t13,a)', advance='no' ) '+----------+'
         write( 21, '(5a)' ) ('----------+', ielem1=ielemmin,ielemmax)

         write( 21, '(t13,a)', advance='no' ) '|Target\Ion|'
         write( 21, '(5(a,i2,a,a2,a))' ) &
              ('Atom', ielem1, '(', elem(ielem1)%name, ')|', ielem1=ielemmin,ielemmax)

         write( 21, '(t13,a)', advance='no' ) '+----------+'
         write( 21, '(5a)' ) ('----------+', ielem1=ielemmin,ielemmax)

         do ielem2 = 1, nelem
            write( 21, '(t13,a,i2,a,a2,a)', advance='no' ) &
                 '|Atom', ielem2, '(', elem(ielem2)%name, ')|'
            write( 21, '(5(es10.3e2,a))' ) &
                 (se%enl(1,ielem2,ielem1), '|', ielem1=ielemmin,ielemmax)
            do ienl = 2, maxval(se%nenl(ielem2,:))
               write( 21, '(t13,a,t24,a,5(es10.3e2,a))' )  '|', '|',  &
                     (se%enl(min(ienl,se%nenl(ielem2,ielem1)),ielem2,ielem1), '|',  &
                      ielem1=ielemmin,ielemmax)
            end do
         end do

         write( 21, '(t13,a)', advance='no' ) '+----------+'
         write( 21, '(5a)' ) ('----------+', ielem1=ielemmin,ielemmax)

      end do
      write( 21, '()' )
!
!  Report xnl
!
      write( 21, '(t13,a)' ) 'xnl:'
      do i = 0, (nelem-1)/5
         ielemmin = 5*i + 1
         ielemmax = min( 5*i+5, nelem )

         write( 21, '(t13,a)', advance='no' ) '+----------+'
         write( 21, '(5a)' ) ('----------+', ielem1=ielemmin,ielemmax)

         write( 21, '(t13,a)', advance='no' ) '|Target\Ion|'
         write( 21, '(5(a,i2,a,a2,a))' ) &
              ('Atom', ielem1, '(', elem(ielem1)%name, ')|', ielem1=ielemmin,ielemmax)

         write( 21, '(t13,a)', advance='no' ) '+----------+'
         write( 21, '(5a)' ) ('----------+', ielem1=ielemmin,ielemmax)

         do ielem2 = 1, nelem
            write( 21, '(t13,a,i2,a,a2,a)', advance='no' ) &
                 '|Atom', ielem2, '(', elem(ielem2)%name, ')|'
            write( 21, '(5(es10.3e2,a))' )  &
                 (se%xnl(1,ielem2,ielem1), '|', ielem1=ielemmin,ielemmax)
            do ixnl = 2, maxval(se%nxnl(ielem2,:))
               write( 21, '(t13,a,t24,a,5(es10.3e2,a))' ) '|', '|', &
                     (se%xnl(min(ixnl,se%nxnl(ielem2,ielem1)),ielem2,ielem1), '|',  &
                      ielem1=ielemmin,ielemmax)
            end do
         end do

         write( 21, '(t13,a)', advance='no' ) '+----------+'
         write( 21, '(5a)' ) ('----------+', ielem1=ielemmin,ielemmax)

      end do
      write( 21, '()' )
!
!  Report facscr
!
      write( 21, '(t13,a)' ) 'facscr:'
      do i = 0, (nelem-1)/5
         ielemmin = 5*i + 1
         ielemmax = min( 5*i+5, nelem )

         write( 21, '(t13,a)', advance='no' ) '+----------+'
         write( 21, '(5a)' ) ('----------+', ielem1=ielemmin,ielemmax)

         write( 21, '(t13,a)', advance='no' ) '|Target\Ion|'
         write( 21, '(5(a,i2,a,a2,a))' ) &
              ('Atom', ielem1, '(', elem(ielem1)%name, ')|', ielem1=ielemmin,ielemmax)

         write( 21, '(t13,a)', advance='no' ) '+----------+'
         write( 21, '(5a)' ) ('----------+', ielem1=ielemmin,ielemmax)

         do ielem2 = 1, nelem
            write( 21, '(t13,a,i2,a,a2,a)', advance='no' ) &
                 '|Atom', ielem2, '(', elem(ielem2)%name, ')|'
            write( 21, '(5(es10.3e2,a))' ) &
                 (se%facscr(1,ielem2,ielem1), '|', ielem1=ielemmin,ielemmax)
            do ifacscr = 2, maxval(se%nfacscr(ielem2,:))
               write( 21, '(t13,a,t24,a,5(es10.3e2,a))' )  '|', '|',  &
                     (se%facscr(min(ifacscr,se%nfacscr(ielem2,ielem1)),ielem2,ielem1), '|',  &
                      ielem1=ielemmin,ielemmax)
            end do
         end do

         write( 21, '(t13,a)', advance='no' ) '+----------+'
         write( 21, '(5a)' ) ('----------+', ielem1=ielemmin,ielemmax)

      end do
      write( 21, '()' )
!
!  Report fac2scr
!
      if( any( se%lfac2scr(:) ) ) then
         call writeTab( se%fac2scr(:,:), 'fac2scr' )
!             ========        i              i
      end if
!
!  Call Firsov electronic stopping model writeout procedure
!
      if( any( se%locmodel(1:nelem,1:nelem) == 'firsov' ) ) then
         call writeTab( se%powfirs(:,:), 'powfirs' )
!             ========        i              i
      end if
!
!  Report electronic straggling
!
      if( any( se%lstrag(:,:) ) ) then

         write( 21, '(t13,a)' ) 'Electronic straggling is enabled'
         write( 21, '()' )
         call writeTab( se%estrag(:,:), 'estrag', undeftext_='no straggle' )
!             ========        i            i                i

      else

         write( 21, '(t13,a,/)' ) 'Electronic straggling is not considered'

      end if
!
!  Report generation of a table
!
      if( se%ltab ) then
         write( 21, '(t13,a,es10.3e2)' ) &
              'Generating a table of stopping powers for p=', se%ptab
         write( 21, '()' )
      end if

      close( 21 )

   end subroutine writeSe
!
!=======================================================================
!
   subroutine defaultSeLind( ielem2, ielem1, corlin, powlin, powint)
!             =============     i       i       o       o       o
!
!     Initializes the default values of the parameters of the Lindhard
!     stopping power. If there are no default values known, the parameters
!     are returned unchanged.
!
!     The default values provided are from experiments using Si targets.
!
!  Arguments
!
      integer(INT4), intent(in) :: ielem2       ! index of atom
      integer(INT4), intent(in) :: ielem1       ! index of projectile
      real(REAL8), intent(inout) :: corlin      ! correction factor
      real(REAL8), intent(inout) :: powlin      ! power of energy dependence
      real(REAL8), intent(inout) :: powint      ! parameter for interpolation with
                                                !    Bethe-Bloch formula
!
!  Local variables
!
!     None
!
!  Start of Executable Code:
!  ------------------------
!
!  Special default values
!
      select case( elem(ielem1)%z )
         case( 1 )                              ! H
!            corlin = 1.4                        !    adjusted to ZBL and KKK
            corlin = 1.7                        !    adjusted to PSTAR
            powint = 0.86
         case( 2 )                              ! He
!            corlin = 1.35                          ! adjusted to ZBL and KKK
            corlin = 1.63                          ! adjusted to Nte23: PRB 107 (2023) 155145
            powint = 0.68                          ! adjusted to DPASS
         case( 3 )                              ! Li
            corlin = 1.4
            powint = 1.05
         case( 4 )                               ! Be
            corlin = 1.4
            powint = 0.95
         case( 5 )                               ! B (note: powint is not from Jiang
!            if( mat(ir)%name == 'SIO2' .or. &   !    but from previous work)
!                mat(ir)%name == 'SiO2' .or. mat(ir)%name == 'sio2' ) then
!               corlin = 1.1                     ! Stolmaijer, JECS 135(9)2309,1988
!            else
            corlin = 1.5
!            end if
            powint = 0.9
         case( 6 )                               ! C
            corlin = 1.5
            powint = 0.87
         case( 7 )                               ! N
!            corlin = 1.5
            corlin = 1.66                          ! adjusted to Nte23: PRB 107 (2023) 155145
            powint = 0.85                          ! adjusted to DPASS
         case( 8 )                               ! O
            corlin = 1.25
            powint = 0.95
         case( 9 )                               ! F
            corlin = 1.2
            powint = 1.
         case( 10 )                              ! Ne
            corlin = 1.1
            powint = 1.07
         case( 11 )                              ! Na
            corlin = 0.9
            powint = 1.35
         case( 12 )                              ! Mg
            corlin = 0.9
            powint = 1.35
         case( 13 )                              ! Al
            corlin = 1.05
            powint = 1.
         case( 14 )                              ! Si
            corlin = 1.25
            powint = 0.9
         case( 15 )                              ! P
!            if( mat(ir)%name == 'SIO2' .or. mat(ir)%name == 'SiO2' .or. &
!                mat(ir)%name == 'sio2' ) then
!               corlin = 1.0                     ! Stolmaijer, JECS 135(9)2309,1988
!            else
            corlin = 1.23
!            endif
            powint = 0.9
         case( 16 )                              ! S
            corlin = 1.4
            powint = 0.82
         case( 17 )                              ! Cl
            corlin = 1.4
            powint = 0.82
         case( 18 )                              ! Ar
            corlin = 1.3
            powint = 0.87
         case( 19 )                              ! K
            corlin = 1.2
            powint = 0.95
         case( 20 )                              ! Ca
            corlin = 1.1
            powint = 1.
         case( 21 )                              ! Sc
            corlin = 1.05
            powint = 1.1
         case( 22 )                              ! Ti
            corlin = 1.
            powint = 1.1
         case( 23 )                              ! V
            corlin = 1.
            powint = 1.1
         case( 24 )                              ! Cr
            corlin = 1.
            powint = 1.09
         case( 25 )                              ! Mn
            corlin = 1.
            powint = 1.09
         case( 26 )                              ! Fe
            corlin = 1.
            powint = 1.09
         case( 27 )                              ! Co
            corlin = 1.
            powint = 1.08
         case( 28 )                              ! Ni
            corlin = 1.
            powint = 1.08
         case( 29 )                              ! Cu
            corlin = 1.
            powint = 1.07
         case( 30 )                              ! Zn
            corlin = 1.
            powint = 1.07
         case( 31 )                              ! Ga
            corlin = 0.95
            powint = 1.12
         case( 32 )                              ! Ge
            corlin = 1.
            powint = 1.06
         case( 33 )                              ! As
            corlin = 1.
            powint = 1.06
         case( 51 )                              ! Sb
            corlin = 1.25
         case( 68 )                              ! Er
            corlin = 0.8
         case( 79 )                              ! Au
            corlin = 1.2
      end select

   end subroutine defaultSeLind
!
!=======================================================================
!
   subroutine defaultXnl( ielem2, ielem1, nnl, enl, xnl )
!             ==========     i      i      o    o    o
!
!     Sets the default values of enl and xnl
!
!  Arguments:
!
      integer(INT4), intent(in) :: ielem2       ! index of the atom
      integer(INT4), intent(in) :: ielem1       ! index of the projectile
      integer(INT4), intent(out) :: nnl         ! number of (enl,xnl) pairs
      real(REAl8), intent(out) :: enl(:)        ! energy corresponding to xnl [eV]
      real(REAL8), intent(out) :: xnl(:)        ! nonlocal fraction
!
!  Local Variables:
!
!     None
!
!  Start of Executable Code:
!  ------------------------
!
!  Default values if nothing is known
!
      nnl = 1
      enl(:) = 1
      xnl(:) = HALF
!
!  Default values derived from literature data (see IIT'2000 paper)
!
      nnl = 2                                   ! usually 2 values are given
      if( size( enl ) < nnl )  goto 900

      select case( elem(ielem1)%z )
!
!  Refined models
!
         case( 5 )
            enl(1) = 1
            xnl(1) = 0.15
            enl(2) = 1000
            xnl(2) = 0.320694                   ! corresponding to qse=0.11

         case( 15 )
            nnl = 3
            if( size( enl ) < nnl )  goto 900
            enl(1) = 1
            xnl(1) = 0.25
            enl(2) = 1e6
            xnl(2) = 0.50
            enl(3) = 0                          ! 0 indicates stopping power maximum
            xnl(3) = 0.949

         case( 33 )                             ! As
            enl(1) = 50000                      ! modified to fit Simonescu, NIMB 100 (1995) 483
            enl(2) = ZERO                       ! 0 indicates stopping power maximum
            xnl(1) = 0.30
            xnl(2) = 0.95

         case default
!
!  From Eisen and similar data
!
            enl(1) = 11660*elem(ielem1)%m       ! Eisen energy
            enl(2) = ZERO                       ! 0 indicates stopping power maximum

            select case( elem(ielem1)%z )
               case( 1 )                        ! H (taken the same as for He)
                  xnl(1) = 0.40
                  xnl(2) = 0.75
               case( 2 )                        ! He (rough estimate)
                  xnl(1) = 0.40
                  xnl(2) = 0.75
               case( 3 )                        ! Li
                  xnl(1) = 0.35
                  xnl(2) = 0.845
               case( 4 )                        ! Be
                  xnl(1) = 0.35
                  xnl(2) = 0.825
               case( 5 )                        ! B
                  xnl(1) = 0.41
                  xnl(2) = 0.804
               case( 6 )                        ! C
                  xnl(1) = 0.47
                  xnl(2) = 0.746
               case( 7 )                        ! N
                  xnl(1) = 0.31
                  xnl(2) = 0.748
               case( 8 )                        ! O
                  xnl(1) = 0.29
                  xnl(2) = 0.790
               case( 9 )                        ! F
                  xnl(1) = 0.12
                  xnl(2) = 0.83
               case( 10 )                             ! Ne
                  xnl(1) = 0.04
                  xnl(2) = 0.88
               case( 11 )                             ! Na
                  xnl(1) = 0.02
                  xnl(2) = 0.90
               case( 12 )                             ! Mg
                  xnl(1) = 0.09
                  xnl(2) = 0.92
               case( 13 )                             ! Al
                  xnl(1) = 0.17
                  xnl(2) = 0.93
               case( 14 )                             ! Si
                  xnl(1) = 0.29
                  xnl(2) = 0.94
               case( 15 )                             ! P
                  xnl(1) = 0.415
                  xnl(2) = 0.949
               case( 16 )                             ! S
                  xnl(1) = 0.43
                  xnl(2) = 0.949
               case( 17 )                             ! Cl
                  xnl(1) = 0.43
                  xnl(2) = 0.949
               case( 18 )                             ! Ar
                  xnl(1) = 0.39
                  xnl(2) = 0.95
               case( 19 )                             ! K
                  xnl(1) = 0.36
                  xnl(2) = 0.95
               case( 20 )                             ! Ca
                  xnl(1) = 0.52
                  xnl(2) = 0.95
               case( 21 )                             ! Sc
                  xnl(1) = 0.55
                  xnl(2) = 0.95
               case( 22 )                             ! Ti
                  xnl(1) = 0.59
                  xnl(2) = 0.95
               case( 23 )                             ! V
                  xnl(1) = 0.52
                  xnl(2) = 0.95
               case( 24 )                             ! Cr
                  xnl(1) = 0.50
                  xnl(2) = 0.95
               case( 25 )                             ! Mn
                  xnl(1) = 0.40
                  xnl(2) = 0.95
               case( 26 )                             ! Fe
                  xnl(1) = 0.30
                  xnl(2) = 0.95
               case( 27 )                             ! Co
                  xnl(1) = 0.30
                  xnl(2) = 0.95
               case( 28 )                             ! Ni
                  xnl(1) = 0.30
                  xnl(2) = 0.95
               case( 29 )                             ! Cu
                  xnl(1) = 0.30
                  xnl(2) = 0.95
               case( 30 )                             ! Zn
                  xnl(1) = 0.30
                  xnl(2) = 0.95
               case( 31 )                             ! Ga
                  xnl(1) = 0.30
                  xnl(2) = 0.95
               case( 32 )                             ! Ge
                  xnl(1) = 0.30
                  xnl(2) = 0.95
               case( 33 )                             ! As
                  xnl(1) = 0.30
                  xnl(2) = 0.95
            end select

      end select

      return

900   error stop 'defaultXnl: array too small'

   end subroutine defaultXnl
!
!=======================================================================
!
end module Param_Se_mod
