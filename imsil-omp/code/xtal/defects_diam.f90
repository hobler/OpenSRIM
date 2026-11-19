module Xtal_Defects_Diam_mod
!
!  Manages information about defect types
!
!  Used Modules:
!  ============
!
   use Const_mod,           only: REAL8, INT4,  &
                                    ONE, EPSMB, ZERO
   use Param_dam_mod,       only: dam
   use Random_mod,          only: rannum
   use Open_mod,            only: openFile
!
!  Default Type of Access:
!  ======================
!
   implicit none
   private 
   save
!
!  Public Procedures:
!  =================
!
   public :: initDef                         ! initialize strain and cumulative probability
!            =======
   public :: getDefType                      ! get a defect type for given number of atoms
!            ==========
   public :: getProbDef
!            ==========
   public :: getDefnAtoms                    ! get number of atoms for defect type
!            ============
   public :: getDefnNeigh                    ! get number of neighbors for defect type
!            ============                          
   public :: getDefnLS                       ! get number of lattice sites for defect type
!            =========
   public :: getDefStrain                    ! get strain
!            ============
   public :: getDefLatPos                    ! get relative position of defect lattice site
!            ============
!
!  Private Procedures:
!  ==================
!
   private :: symmetry                       ! apply symmetry operation to coordinates
!             ========
   interface symmetry
      module procedure SymmetryReal
      module procedure SymmetryInt
   end interface
   private :: writeDefType                   ! write defect to file
!             ============
   private :: readDefType                    ! read defect from file
!             ===========
!
!  Private data:
!  ============
!
   type :: Defect_type
!          ===========
      character(10) :: name
      integer(INT4) :: nsite                    ! number of lattice sites occupied by cluster
      integer(INT4) :: nAtoms(1:4)              ! number of atoms at lattice site atoms
      integer(INT4) :: nNeigh                   ! number of neighbors (valid only for defType(1)
      integer(INT4) :: dist(1:3,1:4)            ! distance between site atoms in cluster
      integer(INT4), allocatable :: ipos(:,:)   ! distance of neighboring site of defType(1)
      real(REAL8), allocatable :: dpos(:,:)     ! strain for atoms
      real(REAL8) :: defects(1:3,1:2,1:4)  
   end type Defect_type

   type(Defect_type) :: defects( 7 ) 

   integer(INT4), allocatable :: defact(:)      ! array of active defects
   real(REAL8), allocatable :: probdefact(:)    ! cummulative probability of active def.
   integer(INT4), allocatable :: numatom(:)     ! number of atoms for defect types
   integer(INT4), allocatable :: numneig(:)     ! number of neighbors for def. types
   integer(INT4), allocatable :: nsite(:)       ! number of lattice sites occupied by def.
   integer(INT4), allocatable :: twoact(:)      ! two atomic active defects

   integer(INT4), parameter :: or(1:3, 1:24) = reshape( [  &
          1,  2,  3,  &    ! orientation  1: ( X, Y, Z)
         -1, -2,  3,  &    ! orientation  2: (-X,-Y, Z) 
          1, -2, -3,  &    ! orientation  3: ( X,-Y,-Z)
         -1,  2, -3,  &    ! orientation  4: (-X, Y,-Z)
          3,  1,  2,  &    ! orientation  5: ( Z, X, Y)
         -3, -1,  2,  &    ! orientation  6: (-Z,-X, Y)
          3, -1, -2,  &    ! orientation  7: ( Z,-X,-Y)
         -3,  1, -2,  &    ! orientation  8: (-Z, X,-Y)
          2,  3,  1,  &    ! orientation  9: ( Y, Z, X)
         -2, -3,  1,  &    ! orientation 10: (-Y,-Z, X)
          2, -3, -1,  &    ! orientation 11: ( Y,-Z,-X)
         -2,  3, -1,  &    ! orientation 12: (-Y, Z,-X)
          3,  2,  1,  &    ! orientation 13: ( Z, Y, X)
         -3, -2,  1,  &    ! orientation 14: (-Z,-Y, X)
          3, -2, -1,  &    ! orientation 15: ( Z,-Y,-X)
         -3,  2, -1,  &    ! orientation 16: (-Z, Y,-X)
          1,  3,  2,  &    ! orientation 17: ( X, Z, Y)
         -1, -3,  2,  &    ! orientation 18: (-X,-Z, Y)
          1, -3, -2,  &    ! orientation 19: ( X,-Z,-Y)
         -1,  3, -2,  &    ! orientation 20: (-X, Z,-Y)
          2,  1,  3,  &    ! orientation 21: ( Y, X, Z)
         -2, -1,  3,  &    ! orientation 22: (-Y,-X, Z)
          2, -1, -3,  &    ! orientation 23: ( Y,-X,-Z)
         -2,  1, -3 ],  &  ! orientation 24: (-Y, X,-Z)
         [ 3, 24 ] )
!____________________________________________________________________________________________
!||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
!||                                                                                       |||
!||     TABLE OF DEFECTS - NUMERICAL                                                      |||
!||                                                                                       |||
!||     0 - LATTICE ATOM                                                                  |||
!||     1 - VACANCY                                                                       |||
!||     2 - SPLIT-<110> INTERSTITIAL                                                      |||
!||     3 - DI-<110> ; DI-INTERSTITIAL MADE ADDING 1 EXTRA ATOM TO SPL-<110>              |||
!||     4 - TRI-<110>; TRI-INTERSTITIAL MADE ADDING 1 EXTRA ATOM TO DI-<110>              |||
!||     5 - HEXAGONAL INTERSTITIAL                                                        |||
!||     6 - TETRAHEDRAL INTERSTITIAL                                                      |||
!||     7 - FOUR1 ; MADE FROM FOUR SPL-<100> INTERSTITIAL!                                |||
!||     8 - FOUR1 ; SECOND site                                                           |||
!||     9 - FOUR1 ; THIRD site                                                            |||
!||     10 - FOUR1 ; FOURTH site                                                          |||
!||                                                                                       |||
!||_______________________________________________________________________________________|||
!||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
!
!===========================================================================================
!
contains
!
!===========================================================================================
!
   subroutine initDef( )
!             =======
!
!     Initialize strain and cummulative probability for defTypes. 
!     
!  Arguments:
!
!     None
!
!  Local Variables:
!
      integer(INT4) :: iact, nact, i, ntwoa, nunit=43, nunitw=44
      real(REAL8) :: probdef(2:7)
      logical :: dact(2:7), lwr=.false.
!
!  Start of Executable Code:
!  ------------------------
!
      call openFile( nunit, "DEFECTS", 'formatted', 'read' )
!          ========   i        i         i          i
!
!  Init the strain
!
      do i = 1, size( defects )
         call readDefType( defects(i), nunit, i )
!             ===========      io        i    i
      end do

      close( nunit )

      if( lwr ) then
         call openFile( nunitw, "DEFECTS", 'formatted', 'write' )
         do i = 1, size( defects )
            call writeDefType( defects(i), nunitw, i )
!                ============     io          i      i
         end do
         error stop "initDef: defect types written in file DEFECTS"
      end if

      iact = 0
      ntwoa = 0
      if( dam%fracspl > EPSMB ) then
         iact = iact + 1
         dact(2) = .true.
         probdef(2) = dam%fracspl
         ntwoa = ntwoa + 1
      end if
      if( dam%fracdi > EPSMB )  then
         iact = iact + 1
         dact(3) = .true.
         probdef(3) = dam%fracdi
      end if
      if( dam%fractri > EPSMB ) then
         iact = iact + 1
         dact(4) = .true.
         probdef(4) = dam%fractri
      end if
      if( dam%frachex > EPSMB ) then
         iact = iact + 1
         dact(5) = .true.
         probdef(5) = dam%frachex
         ntwoa = ntwoa + 1
      end if
      if( dam%fractet > EPSMB ) then
         iact = iact + 1
         dact(6) = .true.
         probdef(6) = dam%fractet
         ntwoa = ntwoa + 1
      end if
      if( dam%fracfour > EPSMB ) then
         iact = iact + 1
         dact(7) = .true.
         probdef(7) = dam%fracfour
      end if

      nact = iact
      allocate( defact( 1:nact ), probdefact( 1:nact ), twoact( 1:ntwoa ) )
      defact(:) = 0
      probdefact(:) = ZERO         

      iact = 0
      ntwoa = 0
      do i = 2, 7
         if( dact(i) ) then
            iact = iact + 1
            defact(iact) = i
            probdefact(iact) = sum( probdef(2:i) ) 
            if( i == 2 .or. i == 5 .or. i == 6 ) then
               ntwoa = ntwoa + 1
               twoact(ntwoa) = i
            end if
         end if
      end do

      if ( iact /= nact ) error stop 'initDef: error'
!
!  Allocate numatom, numneig, nsite
!
      allocate( numatom( 0:10 ), numneig( 0:10 ), nsite( 0:7 ) )
      numatom(0) = 1
      numatom(1:7) = defects(1:7)%natoms(1)
      numatom(8:10) = defects(7)%natoms(2:4)

      numneig(0) = 1
      numneig(1:7) = defects(1:7)%nNeigh
      numneig(8:10) = 2

      nsite(0) = 1
      nsite(1:7) = defects(1:7)%nsite

   end subroutine initDef
!
!===========================================================================================
!
   subroutine getDefType( defType, natom_ )
!             ==========     o       i
!
!     Get the possible defect type for a given number of atoms at the lattice site. 
!
!  Calling arguments:
!
      integer(INT4), intent(out) :: defType          ! defect type
      integer(INT4), intent(in), optional :: natom_  ! natom
!
!  Local Variables:
!
      integer(INT4) :: i
      real(REAL8) :: rand
!
!  Start of Executable Code:
!  ------------------------
!
      if( present( natom_ ) ) then
         if( natom_ == 0 ) then
            defType = 1                     ! vacancy
         elseif( natom_ == 1 ) then
            defType = 0                     ! lattice atom
         elseif( natom_ == 3 ) then
            defType = 3                     ! di-interstititial
         elseif( natom_ == 4 ) then
            defType = 4                     ! tri-interstitial
         elseif( natom_ == 2 ) then          ! may be spl110, tetra, hex, (four??)
            if( size( twoact ) == 1 ) then
               defType = twoact(1)
            elseif ( size( twoact ) == 1 ) then
                                            ! Choose between active defect types
               call rannum( rand )          ! which have 2 atoms at lattice site
!                   ======   o
               defType = twoact( int( size( twoact ) * rand ) + 1 )
            else                            ! default value is spl-<110> if no 
               defType = 2                  ! other twoatomic defect is active
            end if
         else
            error stop 'getDefType: error natom'
         end if
      else                                  ! natom_ not present obtain defect
         if( size( defact ) == 1 ) then     ! according to defect fraction
            defType = defact(1)
         else
            call rannum( rand )
!                ======   o
            do i = 1, size( defact )
               if( rand <= probdefact(i) )  exit
            end do
            defType = defact(i)
         end if
      end if

   end subroutine getDefType
!
!===========================================================================================
!
   subroutine getProbDef( defType, probDefS, probInt )
!             =========     i         o         o
!
!     Get the probability of defect lattice site for a given defType.
!
!  Calling arguments:
!
      integer(INT4), intent(in) :: defType  ! defect type
      real(REAL8), intent(out) :: probDefS  ! prob. of defect site (def+vac)
      real(REAL8), intent(out) :: probInt   ! prob. of interstitial site
                                                ! for already chosen defect site   
!
!  Start of Executable Code:
!  ------------------------
!
      select case( defType )

         case( 2 )
            probDefS = 2.
            probInt = 1./2.
         case( 3 )
            probDefS = 3./2.
            probInt = 1./3.
         case( 4 )
            probDefS = 4./3.
            probInt = 1./4.
         case( 5 )
            probDefS = 2.
            probInt = 1./2.
         case( 6 )
            probDefS = 2.
            probInt = 1./2.
         case( 7 )
            probDefS = 5./4.
            probInt = 1./5.
         case DEFAULT
            error stop 'getProbDef:ERROR'

      end select

   end subroutine getProbDef
!
!===========================================================================================
!
   function getDefnAtoms( defType ) result( natom )
!           ============     i                o
!
!     Give the number of lattice site atoms for defType
!
!  Calling arguments:
!
      integer(INT4), intent(in) :: defType  ! defect type
!
!  Local variables
!
      integer(INT4) :: natom 
!
!  Start of Executable Code:
!  ------------------------
!
      if( defType < 0 .or. defType > size( numatom ) )  error stop 'getDefnAtoms: error'
      natom = numatom( defType )

   end function getDefnAtoms
!
!===========================================================================================
!
   function getDefnNeigh( defType ) result( neigDef ) 
!           ============     i                 o       
!
!     Get the number of neighbors for active defect type & vacancy  
!
!  Arguments:
!
      integer(INT4), intent(in) :: defType  
!
!  Local variables:
!
      integer(INT4) :: neigDef
!
!  Start of Executable Code:
!  ------------------------
!
      if( defType < 0 .or. defType > size( numneig ) )  error stop 'getDefnNeigh: error'
      neigDef = numneig( defType )

   end function  getDefnNeigh
!
!=============================================================================================
!
   function getDefnLS( defType ) result( NLSdefects ) 
!           =========     i                  o          
!
!     Get the number of lattice sites occupied by the defect!
!
!  Arguments:
!
      integer(INT4), intent(in) :: defType  ! defect type
!
!  Local variables:
!
      integer(INT4) :: NLSdefects 
!
!  Start of Executable Code:
!  ------------------------ 
!
      if( defType < 0 .or. defType > size( nsite ) )  error stop 'getDefnLS: error'
      NLSDefects = nsite( defType )

   end function  getDefnLS
!
!=============================================================================================
!
   function getDefStrain( iodd, ineigh, ior, defTyp )  result( dpos ) 
!           ============    i      i     i      i               o        
!
!     Get the strain for defined defect type, neighbor, position and orientation
!
!  Arguments:
!
      integer(INT4), intent(in) :: iodd     ! odd or even site 
      integer(INT4), intent(in) :: ineigh   ! index of neighboring atom
      integer(INT4), intent(in) :: ior      ! defect orientation
      integer(INT4), intent(in) :: defTyp   ! defect type
!
!  Local variables:
!
      real(REAL8) :: dpos(1:3)              ! displacement vectors
!
!  Start of executable code:
!  ------------------------
!
      if( iodd < 0 .or. iodd > 1 )  error stop 'getDefStrain: wrong parity index!'
      if( ineigh > getDefnNeigh( defTyp ) )  &
!                  ============    i
            error stop 'getDefStrain: not enough data available!'
      select case( defTyp )
         case( 0:7 )
            dpos(:) = defects(defTyp)%dpos(:,ineigh)
         case( 8:10 )
            if( ineigh > 2 )  error stop 'getDefStrain: error defTyp 8,9,10'
            dpos(:) = defects(7)%defects(:,ineigh,defTyp-6)
         case default
            error stop 'ERROR getDefStrain!'
      end select
!
!  Find appropriate orientation. Vacancy has not symmetry orientation
!
      if( defTyp /= 1 )  dpos(1:3) = Symmetry( ior, dpos(1:3) )
!                                    ========   i      io
!
!  Even or odd site?
!
      if( iodd == 1 )  dpos(3) = - dpos(3)

   end function getDefStrain
!
!===========================================================================================
!
   function getDefLatPos( iodd, ineigh, ior, defTyp )  result( ipos ) 
!           ============     i      i     i      i                 o        
! 
!     Get the strain for defined defect type, neighbor, position and orientation
!
!     ATTENTION: In case of defTyp=8,9,10 function must be called with ineigh=1
!
!  Arguments:
!
      integer(INT4), intent(in) :: iodd    ! odd or even site 
      integer(INT4), intent(in) :: ineigh  ! index of neighboring atom
      integer(INT4), intent(in) :: ior     ! defect orientation
      integer(INT4), intent(in) :: defTyp ! defect type
!
!  Local variables:
!
      integer(INT4) :: ipos(1:3)! displacement vectors
!
!  Start of executable code:
!  ------------------------
!
      if( iodd < 0 .or. iodd > 1 )  error stop 'getDefLatPos: wrong parity index!'
      if( ineigh > getDefnNeigh( defTyp ) )  &
!                  ============    i
            error stop 'getDefLatPos: not enough data available!'
      select case( defTyp )
         case( 0:7 )
            ipos(:) = defects(defTyp)%ipos(:,ineigh)
         case( 8:10 )
            if( ineigh /= 1 )  error stop 'getDefLatPos: error defTyp 8,9,10'
            ipos(:) = defects(7)%dist(:,defTyp-6)
         case default
            error stop 'ERROR getDefLatPos!'
      end select
!
!  Apply symmetry-equvivalence operation. Vacancy has not symmetry orientation
!
      if( defTyp /= 1)  ipos(1:3) = Symmetry( ior,ipos(1:3) )
!                                   ========   i     io  
!  Even or odd?
!
      if( iodd == 1 )  ipos(3) = - ipos(3)

   end function getDefLatPos
!
!===========================================================================================
!
   function SymmetryReal( i, rold ) result( rnew )
!           ============  i   i               o
!
!     Apply defined symmetry operations (among 24) for given coordinates
!
!  Arguments:
!
      integer(INT4), intent(in) :: i        ! index of symmetry operations
      real(REAL8), intent(in) :: rold(1:3)  ! input coordinates
      real(REAL8) :: rnew(1:3)              ! symmetry coordinates
!
!  Local variables:
! 
      integer(INT4) :: j
!
!  Start of executable code:
!  ------------------------
!
      if( i > 24 .or. i < 1 )  error stop 'symmetryReal: ERROR: i<1 or i>24'

      do j = 1, 3
         rnew(j) = sign( 1, or(j,i) ) * rold( abs(or(j,i)) )
      end do

   end function SymmetryReal
!
!=========================================================================================
!
   function SymmetryInt( i, rold ) result( rnew )
!           ===========  i    i             o
!
!     Apply defined symmetry operations (among 24) for given coordinates
!
!  Calling arguments:
!
      integer(INT4), intent(IN) :: i          ! index of symmetry operations
      integer(INT4), intent(IN) :: rold(1:3)  ! input coordinates
      integer(INT4) :: rnew(1:3)              ! symmetry coordinates
!
!  Local variables:
! 
      integer(INT4) :: j
!
!  Start of executable code:
!  ------------------------
!
      if ( i > 24 .or. i < 1 ) error stop 'symmetryInt: ERROR i<1 or i>24'

      do j = 1, 3
         rnew(j) = sign( 1, or(j,i) ) * rold( abs(or(j,i)) )
      end do

   end function SymmetryInt
!
!===========================================================================================
!
   subroutine writeDefType( defect, nunit, deftype )
!             ============     i       i      i
!
!     Writes the defect ion file
!
!  Arguments:
!
      type(Defect_type), intent(in) :: defect
      integer(INT4), intent(in) :: nunit, deftype  ! unit number
!
!  Local Variables:
!
      integer(INT4) :: i, j, il
!
!  Start of Executable Code:
!  ------------------------
!
      write(nunit,'(t5,a,t49,a)')  trim(defect%name), "=name"
      write(nunit,'(i8,t49,a)') deftype, "=defType"
      write(nunit,'(i8,t49,a)') defect%nsite, "= #defect lattice Sites"
      do i = 1, defect%nsite
         write(nunit,'(i8,t49,a)')  i, "= LatticeSite"
         write(nunit,'(i8,t49,a)') defect%nAtoms(i), "= nAtoms"
         if( i .eq. 1 ) then
            write(nunit,'(i8,t49,a)') defect%nNeigh, "= nNeighbors"
            do il = 1, defect%nNeigh
               write(nunit,'(3i8, t49,a, i2,a)') defect%ipos(1:3,il), "=ipos(:,",il,")"
               write(nunit,'(3f15.8, t49, a)') defect%dpos(1:3,il), "=dpos"
            end do
         else
            write(nunit,'(3i8, t49,a,i2,a)') defect%dist(1:3,i), "=ipos(:,",i,")"
            do j = 1, defect%nAtoms(i)
               write(nunit,'(3f15.8, t49, a)') defect%defects(1:3,j,i), "=dpos"
            end do
         end if
      end do

   end subroutine writeDefType
!
!=======================================================================
!
   subroutine readDefType( defect, nunit, deftype )
!                ===========    io       i      i
!  
!     Writes the defect ion file
!
!  Arguments:
!
      type(Defect_type), intent(inout) :: defect
      integer(INT4), intent(in) :: nunit, deftype  ! unit number
!
!  Local Variables:
!
      integer(INT4) :: i, j, il, i1, status, dType
!
!  Start of Executable Code:
!  ------------------------
!

      read(nunit,*, err=900, end=900) defect%name
      read(nunit,*, err=900, end=900) dType
      if( dtype /= deftype )  error stop "readdeftype: wrong defecttype number"
      read(nunit,*, err=900, end=900) defect%nsite
      do i = 1, defect%nsite
         read(nunit,*, err=900, end=900) i1
         if( i1 /= i )  error stop "readDefType: ERROR wrong LS index" 
         read(nunit,*, err=900, end=900) defect%nAtoms(i)
         if( i .eq. 1 ) then
            read(nunit,*, err=900, end=900) defect%nNeigh

            allocate( defect%dpos( 1:3, 1:defect%nNeigh ), STAT = status )
            if( status /= 0 )  error stop 'readDefType: dposi not allocated'

            allocate( defect%ipos( 1:3, 1:defect%nNeigh ), STAT = status )
            if( status /= 0 )  error stop 'readDefType: iposi not allocated'



            do il = 1, defect%nNeigh
               read(nunit,*, err=900, end=900) defect%ipos(1:3,il)
               read(nunit,*, err=900, end=900) defect%dpos(1:3,il)
            end do
         else
            read(nunit,*, err=900, end=900) defect%dist(1:3,i)
            do j = 1, defect%nAtoms(i)
               read(nunit,*, err=900, end=900) defect%defects(1:3,j,i)
            end do
         end if
      end do

      return

900   error stop 'readDefType: error during reading histogram file'

   end subroutine readDefType
!
!=======================================================================
!
end module Xtal_Defects_Diam_mod
