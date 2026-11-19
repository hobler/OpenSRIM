module Eval_Stat_mod
!
!  Manage statistical results
!
!
!  Used Modules:
!  ============
!
!$ use omp_lib,            only: omp_get_thread_num
   use Const_mod,          only: INT4, REAL8,  &
                                 BIG, ZERO, ONE, TWO, HALF, THIRD, PI,  &
                                 EPSMA,  &
                                 ENUM
   use Array_mod,          only: check_allocation
   use Vector_mod,         only: operator(.dot.), abs
   use Types_mod,          only: Proj_type
   use Undefine_mod,       only: defined

   use Param_Setup_mod,    only: ndim, nelem, nelemIon, nelemTarg, nion, nion1,   &
                                 nionbk, nionhis, nr, nthread, atom1, atom2,  &
                                 mrand, iarand, ibrand, irand,  &
                                 ldamdyn, ldyn, usedef, useion, usebk, usegeom, bkfile,  &
                                 usehis, hisfile, rng, usecell
   use Param_Geom_mod,     only: geom
   use Param_Elem_mod,     only: elem
   use Param_Ion_mod,      only: ion
   use Param_Xtal_mod,     only: xtal
   use Param_Mat_mod,      only: mat, coord_xtal, irxtal
   use Param_Dam_mod,      only: dam
   use Param_Sn_mod,       only: sn
   use Param_Out_mod,      only: out

   use Geom_BndCond_mod,   only: applyBndCond1d, applyBndCond3d
   use Geom_CoordSys_mod,  only: CoordSys_type, initCoordSys, toCoordSys
   use Geom_Polygon_mod,   only: initPolygon

   use Stat_Mom1d_mod,     only: Mom1d_type, initMom1d, scoreMom1d, updateMom1d, writeMom1d,  &
                                 mergeIntoMasterMom1d
   use Stat_Mom2d_mod,     only: Mom2d_type, initMom2d, incrMom2d, writeMom2d,  &
                                 mergeIntoMasterMom2d
   use Stat_His1d_mod,     only: His1d_type, initHis1d, scoreHis1d, writeHis1d,  &
                                 saveHis1d, loadHis1d, getHis1d, readHis1d,  &
                                 mergeHis1d, mergeIntoMasterHis1d, broadcastHis1d
   use Stat_His2d_mod,     only: His2d_type, initHis2d, scoreHis2d, writeHis2d,  &
                                 saveHis2d, loadHis2d, getHis2d,  &
                                 mergeHis2d, mergeIntoMasterHis2d, broadcastHis2d
   use Stat_His3d_mod,     only: His3d_type, initHis3d, scoreHis3d, writeHis3d,  &
                                 saveHis3d, loadHis3d, getHis3d,  &
                                 mergeHis3d, mergeIntoMasterHis3d, broadcastHis3d
   use Stat_His_Util_mod,  only: getAmoCPDD, deleteInAir
   use Eval_RRobin_mod,    only: rob10, rob25, rob30  ! should go into eval-mod, but 
                                                      !    would complicate things
!
   use Ion_Proj_mod,       only: getXStart
   use Open_mod,           only: openFile
   use Random_mod,         only: ranset, ranget, rannum
   use Time_mod,           only: initCPUtime, getCPUtime, getElapsedTime
   use Eval_Det_mod,       only: det10, det30, his10, his30
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
   public :: stat10                 ! initialization
!            ======
   public :: stat11                 ! score at start of cascade
!            ======
   public :: stat12                 ! score at start of projectile trajectory
!            ======
   public :: stat22                 ! score along the trajectory
!            ======
   public :: stat2x                 ! score at the beginning of a trajectory, when a recoil is
!            ======                 !    not followed, or when a projectile is stopped
   public :: stat24                 ! score at end of trajectory of a backscattered projectile
!            ======
   public :: stat25                 ! score at end of trajectory of a transmitted projectile
!            ======
   public :: stat26                 ! perform actions after each ion molecule
!            ======
   public :: stat27                 ! perform actions after a chunk of ion molecules
!            ======
   public :: stat2e                 ! perform actions after ion has been lost
!            ======
   public :: stat30                 ! output of statistical quantities
!            ======
   public :: stat40                 ! write backup file
!            ======

   public :: scoreIV                ! score for I-V pair distribution function
!            =======
   public :: getDamage              ! return interstitial and vacancy concentrations
!            =========
!
!  Private Procedures:
!  ==================
!
   private :: stat50                ! read from backup file
!             ======
   private :: scoreEED              ! score electronic energy deposition
!             ========
   private :: scoreNED              ! score nuclear energy deposition
!             ========
   private :: updateMoments         ! update moments after a set of dependent scores has been
!             =============         !    completed
   private :: updateDamageHistogram ! update his* histograms from his*_tmp (after nionhis ion
!             ===================== !    molecules if ldamdyn=t and nionhis>0)
   private :: updateAmo             ! update amorphous zones considering IBIIA model
!             =========
   private :: setHisAmzon           ! set histogram of amorphous zones
!             ===========
   private :: mergeIntoMasterStat   ! Merge statistics into master thread
!             ===================   !    (his* only if ldamdyn=f)
!
!  Public Data:
!  ===========
!
!  debug
      public sumdee, sumedam, sumeb, sumet
!  end debug
!
!  Private Data:
!  ============
!
!  Number of elements for the statistics
!  nelemIon ... number of elements in the ion (= 1 for atomic ions)
!  nelemTarg ... number of elements that may become targets
!                = nelem for ldyn=t
!                = nelem-nelemIon for ldyn=f
!  Note that for ldyn=t implanted and recoiled ion atom species are considered separate species
!
   integer(INT4) :: nelemIn      ! number of stopped atom species and vacancies
                                 !    = nelemIon or nelemIon+nelemTarg or nelemIon+2*nelemTarg
                                 !    1..nelemIon ... ion atoms
                                 !    nelemIon+1..nelemIon+nelemTarg ... target atoms
                                 !    nelemIon+nelemTarg+1..nelemIon+2*nelemTarg ... vacancies
   integer(INT4) :: nelemAtom    ! number of ion and target atom species
                                 !    = nelemIon+nelemTarg
                                 !    1..nelemIon ... ion atoms
                                 !    nelemIon+1..nelemIon+nelemTarg ... target atoms
   integer(INT4) :: nelemProj    ! number of projectile species (= atoms whose trajectories are
                                 !    followed
                                 !    = nelemIon or nelemIon+nelemTarg
                                 !    1..nelemIon ... ion atoms
                                 !    nelemIon+1..nelemIon+nelemTarg ... target atoms
   real(REAL8), parameter :: HISMIN = 0.0001

   type(CoordSys_type) :: coord_statis
!$omp threadprivate( coord_statis )
!
!  Particle counting
!
   real(REAL8), allocatable :: rin(:)        ! (nelemIn)
   real(REAL8), allocatable :: redam(:)      ! (nelemAtom)
   integer(INT4), allocatable :: nback(:)    ! (nelemProj)
   integer(INT4), allocatable :: ntrans(:)   ! (nelemProj)
   integer(INT4), allocatable :: nlost(:)    ! (nelemProj)
!$omp threadprivate( rin, redam, nback, ntrans, nlost )
!
!  Energy losses
!
   real(REAL8), allocatable :: sumdee(:)     ! (nelemAtom) sum of electronic energy loss of
                                             !    projectile or non-followed recoil
   real(REAL8), allocatable :: sumden(:,:)   ! (nelemAtom,nelemProj) sum of energy loss of
                                             !    projectile to recoil
   real(REAL8), allocatable :: sumedam(:)    ! (nelemAtom) sum of nuclear energy deposition
                                             !    into recoil (not including recoils that
                                             !    are followed)
   real(REAL8), allocatable :: sumeb(:)      ! (nelemProj) average energy of backscattered
                                             !    atoms
   real(REAL8), allocatable :: sumet(:)      ! (nelemProj) average energy of transmitted
                                             !    atoms
!$omp threadprivate( sumdee, sumden, sumedam, sumeb, sumet )
!
!  Amorphous zones (libiia=t)
!
   integer(INT4) :: namo = 0                    ! Number of amorphous zones
   real(REAL8), allocatable :: zamo_left(:)     ! (nbox+2) left boundary of amorphous zone
   real(REAL8), allocatable :: zamo_right(:)    ! (nbox+2) right boundary of amorphous zone
   real(REAL8), allocatable :: dzamo_left(:)    ! (nbox+2) displacement of zamo_left
   real(REAL8), allocatable :: dzamo_right(:)   ! (nbox+2) displacement of zamo_right
!$omp threadprivate( dzamo_left, dzamo_right )
!
!  Global moments (spatial moments unless otherwise noted)
!
   type(Mom1d_type), allocatable :: momyldin(:) ! (nelemIn) yield of stopped atoms and vacancies
   type(Mom1d_type), allocatable :: momylde(:)  ! (nelemAtom) yield of nuclear energy deposition
   type(Mom1d_type), allocatable :: momyldb(:)  ! (nelemProj) yield of backscattered atoms
   type(Mom1d_type), allocatable :: momyldt(:)  ! (nelemProj) yield of transmitted atoms
!$omp threadprivate( momyldin, momylde, momyldb, momyldt )

   type(Mom1d_type), allocatable :: momz(:)     ! (nelemIn) stopped atoms and vacancies
   type(Mom1d_type), allocatable :: momx(:)     ! (nelemIn) stopped atoms and vacancies
   type(Mom1d_type), allocatable :: momy(:)     ! (nelemIn) stopped atoms and vacancies
!$omp threadprivate( momz, momx, momy )

   type(Mom2d_type), allocatable :: momxy(:)    ! (nelemIn) stopped atoms and vacancies
   type(Mom2d_type), allocatable :: momzx(:)    ! (nelemIn) stopped atoms and vacancies
   type(Mom2d_type), allocatable :: momzy(:)    ! (nelemIn) stopped atoms and vacancies
!$omp threadprivate( momxy, momzx, momzy )

   type(Mom1d_type), allocatable :: momzee(:)   ! (nelemAtom) electronic energy deposition
   type(Mom1d_type), allocatable :: momxee(:)   ! (nelemAtom) electronic energy deposition
   type(Mom1d_type), allocatable :: momyee(:)   ! (nelemAtom) electronic energy deposition
!$omp threadprivate( momzee, momxee, momyee )

   type(Mom2d_type), allocatable :: momxyee(:)  ! (nelemAtom) electronic energy deposition
   type(Mom2d_type), allocatable :: momzxee(:)  ! (nelemAtom) electronic energy deposition
   type(Mom2d_type), allocatable :: momzyee(:)  ! (nelemAtom) electronic energy deposition
!$omp threadprivate( momxyee, momzxee, momzyee )

   type(Mom1d_type), allocatable :: momzne(:)   ! (nelemAtom) nuclear energy deposition
   type(Mom1d_type), allocatable :: momxne(:)   ! (nelemAtom) nuclear energy deposition
   type(Mom1d_type), allocatable :: momyne(:)   ! (nelemAtom) nuclear energy deposition
!$omp threadprivate( momzne, momxne, momyne )

   type(Mom2d_type), allocatable :: momxyne(:)  ! (nelemAtom) nuclear energy deposition
   type(Mom2d_type), allocatable :: momzxne(:)  ! (nelemAtom) nuclear energy deposition
   type(Mom2d_type), allocatable :: momzyne(:)  ! (nelemAtom) nuclear energy deposition
!$omp threadprivate( momxyne, momzxne, momzyne )

   type(Mom1d_type), allocatable :: momab(:)    ! (nelemProj) 2D projection angle of
                                                !    backscattered atoms
   type(Mom1d_type), allocatable :: momaab(:)   ! (nelemProj) angle wrt to analyzer direction
                                                !    of backscattered atoms
   type(Mom1d_type), allocatable :: momeb(:)    ! (nelemProj) energy of backscattered atoms
   type(Mom1d_type), allocatable :: momxb(:)    ! (nelemIn) backscattered atoms and vacancies
   type(Mom1d_type), allocatable :: momyb(:)    ! (nelemIn) backscattered atoms and vacancies
   type(Mom2d_type), allocatable :: momxyb(:)   ! (nelemIn) backscattered atoms and vacancies
!$omp threadprivate( momab, momaab, momeb, momxb, momyb, momxyb )

   type(Mom1d_type), allocatable :: momat(:)    ! (nelemProj) 2D projection angle of
                                                !    transmitted atoms
   type(Mom1d_type), allocatable :: momaat(:)   ! (nelemProj) angle wrt to analyzer direction
                                                !    of transmitted atoms
   type(Mom1d_type), allocatable :: momet(:)    ! (nelemProj) energy of transmitted atoms
   type(Mom1d_type), allocatable :: momxt(:)    ! (nelemIn) transmitted atoms and vacancies
   type(Mom1d_type), allocatable :: momyt(:)    ! (nelemIn) transmitted atoms and vacancies
   type(Mom2d_type), allocatable :: momxyt(:)   ! (nelemIn) transmitted atoms and vacancies
!$omp threadprivate( momat, momaat, momet, momxt, momyt, momxyt )
!
!  Histograms: hisxx1y ... 1d, hisxx2y ... 2d, hisxx3y ... 3d;
!              y=(empty) ... inside, y=b ... backscattered, y=t ... transmitted
!
!  Spatial histograms for quantities inside the target
!
   type(His1d_type) :: his1                     ! (nelemIn,nbox) stopped atoms and vacancies
   type(His2d_type) :: his2                     ! (nelemIn,nbox2(1),nbox2(2)) stopped atoms
                                                !    and vacancies
   type(His3d_type) :: his3                     ! (nelemIn,nbox3(1),nbox3(2),nbox3(3))
   type(His1d_type) :: his1_tmp                 ! same as his*, but containing the
   type(His2d_type) :: his2_tmp                 !    concentrations since the last update
   type(His3d_type) :: his3_tmp                 !    of his* if nionhis>0
!$omp threadprivate( his1, his2, his3, his1_tmp, his2_tmp, his3_tmp )

   type(His1d_type) :: hisee1                   ! (nelemAtom,nbox) electronic energy deposition
   type(His2d_type) :: hisee2                   ! (nelemAtom,nbox2(1),nbox2(2)) electronic
                                                !    energy deposition
   type(His3d_type) :: hisee3                   ! (nelemAtom,nbox3(1),nbox3(2),nbox3(3))
                                                !    electronic energy deposition
!$omp threadprivate( hisee1, hisee2, hisee3 )

   type(His1d_type) :: hisne1                   ! (nelemAtom,nbox) nuclear energy deposition
   type(His2d_type) :: hisne2                   ! (nelemAtom,nbox2(1),nbox2(2)) nuclear energy
                                                !    deposition
   type(His3d_type) :: hisne3                   ! (nelemAtom,nbox3(1),nbox3(2),nbox3(3))
                                                !    nuclear energy deposition
!$omp threadprivate( hisne1, hisne2, hisne3 )

   type(His1d_type) :: hisa1                    ! (nelemIn,nbox) auxiliary histogram for
                                                !     amorphous zones
   type(His1d_type) :: hism1                    ! (3*nelemAtom,nbox) displacement moment
   type(His1d_type) :: hisp1                    ! (3*nelemAtom,nbox) deposited momentum
!$omp threadprivate( hisa1, hism1, hisp1 )
!
!  Displacement and pair distribution histograms
!
   type(His1d_type) :: hisd1                    ! (nelemAtom,nboxd) displacement

   type(His1d_type) :: hisiv1                   ! (nelem,nboxd) if ldyn=t
                                                ! (nelemIon+1:nelem,nboxd) otherwise
                                                !     pair distribution
   real(REAL8), allocatable :: rvaciv(:)        ! (nelem) if ldyn=t
                                                ! (nelemIon+1:nelem) otherwise
                                                !     total number of vacancies
!$omp threadprivate( hisd1, hisiv1, rvaciv )
!
!  Histograms of backscattered and transmitted atoms
!
   type(His1d_type) :: his1b                    ! (nelemIn,nbox) backscattered atoms and
                                                !    vacancies
   type(His2d_type) :: his2b                    ! (nelemIn,nbox2(1),nbox2(2)) backscattered
                                                !    atoms and vacancies
                                                ! 1..nelemIon+nelemTarg: where atom leaves
                                                !    boundary layer
                                                ! > nelemIon+nelemTarg: where recoil originates
                                                !    from
   type(His1d_type) :: hise1b                   ! (nelemProj,nboxe) energy of backscattered
                                                !    atoms
   type(His1d_type) :: hisa1b                   ! (nelemProj,nboxa) 2D projection angle of
                                                !    backscattered atoms
   type(His1d_type) :: hisaa1b                  ! (nelemProj,nboxa) angle wrt analyzer direction
                                                !    of backscattered atoms
   type(His2d_type) :: hisa2b                   ! (nelemProj,nboxa2(1),nboxa2(2)) direction
                                                !    cosines of backscattered atoms
!$omp threadprivate( his1b, his2b, hise1b, hisa1b, hisaa1b, hisa2b )

   type(His1d_type) :: his1t                    ! (nelemIn,nbox) transmitted atoms and vacancies
   type(His2d_type) :: his2t                    ! (nelemIn,nbox2(1),nbox2(2)) transmitted atoms
                                                !    and vacancies
                                                ! 1..nelemIon+nelemTarg: where atom leaves
                                                !    boundary layer
                                                ! > nelemIon+nelemTarg: where recoil originates
                                                !    from
   type(His1d_type) :: hise1t                   ! (nelemProj,nboxe) energy of transmitted atoms
   type(His1d_type) :: hisa1t                   ! (nelemProj,nboxa) 2D projection angle of
                                                !    transmitted atoms
   type(His1d_type) :: hisaa1t                  ! (nelemProj,nboxa) angle wrt analyzer direction
                                                !    of transmitted atoms
   type(His2d_type) :: hisa2t                   ! (nelemProj,nboxa2(1),nboxa2(2)) direction
                                                !    cosines of transmitted atoms
!$omp threadprivate( his1t, his2t, hise1t, hisa1t, hisaa1t, hisa2t )

contains
!
!=======================================================================
!
   subroutine stat10( )
!             ======
!
!     Initialize the statistic quantities
!
!     This subroutine must be called from within a parallel region
!
!  Arguments:       
!
!     None
!
!  Local Variables:
!
      integer(INT4) :: ielem, ielem1, i
      integer(INT4), allocatable :: indhisIn(:), indhisRec(:), indhisRec3(:), indhisProj(:)
      integer(INT4), allocatable :: indhisVac(:)
      real(REAL8) :: hbeg, hbeg2(2), hbeg3(3)
!
!  Start of Executable Code:
!  ------------------------
!
!  Call round robin subroutine
!
      if( geom%lrobin )  call rob10( )
!                             =====
!
!  Set number of species that can stop inside the target (including vacancies)
!  Note the different meaning of nelemTarg depending on ldyn
!
!$omp master
      if( dam%lrcoil ) then
         nelemIn = nelemIon + 2*nelemTarg
      else if( dam%ldam ) then
         nelemIn = nelemIon + nelemTarg
      else
         nelemIn = nelemIon
      end if
!
!  Set number of atom species in the simulation
!
      nelemAtom = nelemIon + nelemTarg
!
!  Set number of species that can be moving projectiles
!
      if( dam%lrcoil ) then
         nelemProj = nelemIon + nelemTarg
      else
         nelemProj = nelemIon
      end if
!$omp end master
!$omp barrier
!
!  Index arrays for histograms
!
      if( ldyn ) then
         ielem1 = 1
      else
         ielem1 = nelemIon + 1
      end if

      allocate( indhisIn( nelemIn ) )
      indhisIn(1:nelemIon) = elem(1:nelemIon)%z
      if( dam%lrcoil ) then
         indhisIn(nelemIon+1:nelemIon+nelemTarg) = elem(ielem1:nelem)%z
         indhisIn(nelemIon+nelemTarg+1:nelemIn) = - elem(ielem1:nelem)%z
      else if( dam%ldam ) then
         indhisIn(nelemIon+1:nelemIon+nelemTarg) = elem(ielem1:nelem)%z
      end if

      allocate( indhisRec( 1:nelemAtom ) )
      indhisRec(1:nelemIon) = elem(1:nelemIon)%z
      indhisRec(nelemIon+1:nelemAtom) = elem(ielem1:nelem)%z

      allocate( indhisRec3( 1:3*nelemAtom ) )
      do i = 1, 3
         indhisRec3(i:3*nelemIon:3) = elem(1:nelemIon)%z
         indhisRec3(3*nelemIon+i:3*nelemAtom:3) = elem(ielem1:nelem)%z
      end do

      allocate( indhisProj( nelemProj ) )
      indhisProj(1:nelemIon) = elem(1:nelemIon)%z
      if( dam%lrcoil ) then
         indhisProj(nelemIon+1:nelemProj) = elem(ielem1:nelem)%z
      end if

      allocate( indhisVac( ielem1:nelem ) )
      indhisVac(:) = - elem(ielem1:nelem)%z
!
!  Initialize reference coordinate system for statistics
!
      select case( out%mcoord )
         case( 3 )
            call initCoordSys( coord_statis, rotate_=ion%rotate, tiltx_=ion%tiltx,  &
!                ============       o                i                  i
                               tilty_=ion%tilty )
!                                      i
         case( 4 )
            call initCoordSys( coord_statis, out%dira)
!                ============       o          i
         case( 5 )
            if( .not. mat(irxtal)%lpoly )  coord_statis = coord_xtal
      end select
!
!  Initialize ions inside, backscattered, and transmitted
!
      allocate( rin( nelemIn ) )
      allocate( redam( nelemAtom ) )
      allocate( nback( nelemProj ) )
      allocate( ntrans( nelemProj ) )
      allocate( nlost( nelemProj ) )

      rin(:) = ZERO
      redam(:) = ZERO
      nback(:) = 0
      ntrans(:) = 0
      nlost(:) = 0
!
!  Initialize energy loss sums
!
      allocate( sumdee( nelemAtom ) )
      allocate( sumden( nelemAtom, nelemProj ) )
      allocate( sumedam( nelemAtom ) )
      allocate( sumeb( nelemProj ) )
      allocate( sumet( nelemProj ) )

      sumdee(:) = ZERO
      sumden(:,:) = ZERO
      sumedam(:) = ZERO
      sumeb(:) = ZERO
      sumet(:) = ZERO
!
!  Initialize the yield moments
!
      allocate( momyldin( nelemIn ) )
      do ielem = 1, nelemIn
         call initMom1d( momyldin(ielem) )
!             =========       o
      end do

      allocate( momylde( nelemAtom ) )
      do ielem = 1, nelemAtom
         call initMom1d( momylde(ielem) )
!             =========       o
      end do

      allocate( momyldb( nelemProj ) )
      do ielem = 1, nelemProj
         call initMom1d( momyldb(ielem) )
!             =========       o
      end do

      allocate( momyldt( nelemProj ) )
      do ielem = 1, nelemProj
         call initMom1d( momyldt(ielem) )
!             =========       o
      end do
!
!  Initialize the spatial moments of stopped atoms
! 
      if( out%lmom ) then

         allocate( momz( nelemIn ) )
         do ielem = 1, nelemIn
            call initMom1d( momz(ielem) )
!                =========     o
         end do

         if( out%nlatmom > 0 ) then
            allocate( momx( nelemIn ) )
            do ielem = 1, nelemIn
               call initMom1d( momx(ielem) )
!                   =========     o
            end do

            allocate( momy( nelemIn ) )
            do ielem = 1, nelemIn
               call initMom1d( momy(ielem) )
!                   =========     o
            end do
         end if

         if( out%nmixmom > 0 ) then
            allocate( momxy( nelemIn ) )
            allocate( momzx( nelemIn ) )
            allocate( momzy( nelemIn ) )
            do ielem = 1, nelemIn
               call initMom2d( momxy(ielem) )
!                   =========      o
               call initMom2d( momzx(ielem) )
!                   =========      o
               call initMom2d( momzy(ielem) )
!                   =========      o
            end do
         end if
      end if
!
!  Initialize the spatial moments of electronic energy deposition
!
      if( out%lmomee ) then

         allocate( momzee( nelemAtom ) )
         do ielem = 1, nelemAtom
            call initMom1d( momzee(ielem) )
!                =========      o
         end do

         if( out%nlatmom > 0 ) then
            allocate( momxee( nelemAtom ) )
            do ielem = 1, nelemAtom
               call initMom1d( momxee(ielem) )
!                   =========      o
            end do

            allocate( momyee( nelemAtom ) )
            do ielem = 1, nelemAtom
               call initMom1d( momyee(ielem) )
!                   =========      o
            end do
         end if

         if( out%nmixmom > 0 ) then
            allocate( momxyne( nelemAtom ) )
            allocate( momzxne( nelemAtom ) )
            allocate( momzyne( nelemAtom ) )
            do ielem = 1, nelemAtom
               call initMom2d( momxyee(ielem) )
!                   =========       o
               call initMom2d( momzxee(ielem) )
!                   =========       o
               call initMom2d( momzyee(ielem) )
!                   =========       o
            end do
         end if
      end if
!
!  Initialize the spatial moments of nuclear energy deposition
!
      if( out%lmomne ) then

         allocate( momzne( nelemAtom ) )
         do ielem = 1, nelemAtom
            call initMom1d( momzne(ielem) )
!                =========      o
         end do

         if( out%nlatmom > 0 ) then
            allocate( momxne( nelemAtom ) )
            do ielem = 1, nelemAtom
               call initMom1d( momxne(ielem) )
!                   =========      o
            end do

            allocate( momyne( nelemAtom ) )
            do ielem = 1, nelemAtom
               call initMom1d( momyne(ielem) )
!                   =========      o
            end do
         end if

         if( out%nmixmom > 0 ) then
            allocate( momxyne( nelemAtom ) )
            allocate( momzxne( nelemAtom ) )
            allocate( momzyne( nelemAtom ) )
            do ielem = 1, nelemAtom
               call initMom2d( momxyne(ielem) )
!                   =========       o
               call initMom2d( momzxne(ielem) )
!                   =========       o
               call initMom2d( momzyne(ielem) )
!                   =========       o
            end do
         end if
      end if
!
!  Initialize the backscattered ions moments
!
      if( out%lmomaab ) then
         allocate( momaab( nelemProj ) )
         do ielem = 1, nelemProj
            call initMom1d( momaab(ielem) )
!                =========    o
         end do
      end if

      if( out%lmomab ) then
         allocate( momab( nelemProj ) )
         do ielem = 1, nelemProj
            call initMom1d( momab(ielem) )
!                =========    o
         end do
      end if

      if( out%lmomeb ) then
         allocate( momeb( nelemProj ) )
         do ielem = 1, nelemProj
            call initMom1d( momeb(ielem) )
!                =========    o
         end do
      end if

      if( out%lmomb ) then
         if( out%nlatmom > 0 ) then
            allocate( momxb( nelemIn ) )
            allocate( momyb( nelemIn ) )
         end if
         if( out%nmixmom > 0 )  allocate( momxyb( nelemIn ) )

         do ielem = 1, nelemIn
            if( out%nlatmom > 0 ) then
               call initMom1d( momxb(ielem) )
!                   =========    o
               call initMom1d( momyb(ielem) )
!                   =========    o
            end if
            if( out%nmixmom > 0 )  call initMom2d( momxyb(ielem) )
!                                       =========    o
         end do
      end if
!
!  Initialize the transmitted ions moments
!
      if( out%lmomaat ) then
         allocate( momaat( nelemProj ) )
         do ielem = 1, nelemProj
            call initMom1d( momaat(ielem) )
!                =========    o
         end do
      end if

      if( out%lmomat ) then
         allocate( momat( nelemProj ) )
         do ielem = 1, nelemProj
            call initMom1d( momat(ielem) )
!                =========    o
         end do
      end if

      if( out%lmomet ) then
         allocate( momet( nelemProj ) )
         do ielem = 1, nelemProj
            call initMom1d( momet(ielem) )
!                =========    o
         end do
      end if

      if( out%lmomt ) then
         if( out%nlatmom > 0 ) then
            allocate( momxt( nelemIn ) )
            allocate( momyt( nelemIn ) )
         end if
         if( out%nmixmom > 0 )  allocate( momxyt( nelemIn ) )

         do ielem = 1, nelemIn
            if( out%nlatmom > 0 ) then
               call initMom1d( momxt(ielem) )
!                   =========   o
               call initMom1d( momyt(ielem) )
!                   =========   o
            end if
            if( out%nmixmom > 0 )  call initMom2d( momxyt(ielem) )
!                                       =========     o
         end do
      end if
!
!  Box boundaries for the histograms
!
      if( irxtal > 0 ) then            ! slightly left of origin
         hbeg = -100 * EPSMA * out%nbox * out%wbox
         hbeg2(:) = -100 * EPSMA * out%nbox2(:) * out%wbox2(:)
         hbeg3(:) = -100 * EPSMA * out%nbox3(:) * out%wbox3(:)
      else
         hbeg = ZERO
         hbeg2(:) = ZERO
         hbeg3(:) = ZERO
      end if

      if( any( geom%bndtypex(:) == ENUM%bndtype%periodic .and. .not.out%lextendx(:) ) ) then
         hbeg2(1) = geom%xlim(1)
         hbeg3(1) = geom%xlim(1)
      end if

      if( any( geom%bndtypey(:) == ENUM%bndtype%periodic .and. .not.out%lextendy(:) ) ) then
         hbeg3(1) = geom%ylim(1)
      end if

      if( any( geom%bndtypez(:) == ENUM%bndtype%periodic .and. .not.out%lextendz(:) ) ) then
         hbeg = geom%zlim(1)
         hbeg2(1) = geom%zlim(1)
         hbeg3(1) = geom%zlim(1)
      end if
!
!  Initialize the 1d histogram
!
      if( out%lhis ) then
         call initHis1d( indhisIn, out%nbox, out%wbox, hbeg, his1 )
!             =========      i        i         i       i     o
         if( nionhis > 0 .and. dam%ndamdim == 1 )  &
               call initHis1d( indhisIn, out%nbox, out%wbox, hbeg, his1_tmp )
!                   =========      i        i         i       i       o
      end if
!
!  Initizalize other 1d histograms for projectiles landing inside the target
!
      if( dam%lamzon )  &
            call initHis1d( indhisIn, out%nbox, out%wbox, hbeg, hisa1 )
!                =========      i        i         i       i      o
      if( out%lhisd )  &
            call initHis1d( indhisRec, out%nboxd, out%wboxd, ZERO, hisd1 )
!                =========      i          i          i       i      o
      if( out%lhisiv ) then
         call initHis1d( indhisVac, out%nboxd, out%wboxd, ZERO, hisiv1 )
!             =========      i          i          i       i      o
         allocate( rvaciv( ielem1:nelem ) )
         rvaciv(:) = ZERO
      end if
      if( out%lhism )  &
            call initHis1d( indhisRec3, out%nbox, out%wbox, hbeg, hism1 )
!                =========      i          i         i       i      o
      if( out%lhisee )  &
            call initHis1d( indhisRec, out%nbox, out%wbox, hbeg, hisee1 )
!                =========      i         i         i       i      o
      if( out%lhisne )  &
            call initHis1d( indhisRec, out%nbox, out%wbox, hbeg, hisne1 )
!                =========      i         i         i       i      o
      if( out%lhisp )  &
            call initHis1d( indhisRec3, out%nbox, out%wbox, hbeg, hisp1 )
!                =========      i          i         i       i      o
!
!  Initialize the 2d histogram 
!
      if( out%lhis2 ) then
         call initHis2d( indhisIn, out%nbox2, out%wbox2, hbeg2, his2 )
!             =========     i          i          i        i     o
         if( nionhis > 0 .and. dam%ndamdim == 2 )  &
               call initHis2d( indhisIn, out%nbox2, out%wbox2, hbeg2, his2_tmp )
!                   =========     i          i          i        i       o
      end if
!
!  Initialize other 2D histograms
!
      if( out%lhisee2 )  &
            call initHis2d( indhisRec, out%nbox2, out%wbox2, hbeg2, hisee2 )
!                =========      i          i          i        i      o
      if( out%lhisne2 )  &
            call initHis2d( indhisRec, out%nbox2, out%wbox2, hbeg2, hisne2 )
!                =========      i          i          i        i      o
!
!  Initialize the 3d histogram
!
      if( out%lhis3 ) then
         call initHis3d( indhisIn, out%nbox3, out%wbox3, hbeg3, his3 )
!             =========     i          i          i        i     o
         if( nionhis > 0 .and. dam%ndamdim == 3 )  &
               call initHis3d( indhisIn, out%nbox3, out%wbox3, hbeg3, his3_tmp )
!                   =========     i          i          i        i       o
      end if
!
!  Initialize other 3d histograms
!
      if( out%lhisee3 )  &
            call initHis3d( indhisRec, out%nbox3, out%wbox3, hbeg3, hisee3 )
!                =========      i          i          i        i      o
      if( out%lhisne3 )  &
            call initHis3d( indhisRec, out%nbox3, out%wbox3, hbeg3, hisne3 )
!                =========      i          i          i        i      o
!
!  Initialize the backscattered ions statistics
!
      if( out%lhisaab )  &
            call initHis1d( indhisProj, out%nboxa, out%wboxa, ZERO, hisaa1b )
!                =========      i           i          i       i       o
      if( out%lhisab )  &
            call initHis1d( indhisProj, 2*out%nboxa, out%wboxa, ZERO, hisa1b )
!                =========      i            i           i       i      o
      if( out%lhisa2b )  &
            call initHis2d( indhisProj, out%nboxa2, out%wboxa2, [ZERO, ZERO], hisa2b )
!                =========      i           i           i            i          o
      if( out%lhiseb )  &
            call initHis1d( indhisProj, out%nboxe, out%wboxe, ZERO, hise1b )
!                =========      i           i          i       i      o
      if( out%lhisb )  &
            call initHis1d( indhisIn, out%nbox, out%wbox, ZERO, his1b )
!                =========     i         i         i       i      o
      if( out%lhis2b )  &
            call initHis2d( indhisIn, out%nbox2, out%wbox2, [ZERO, ZERO], his2b )
!                =========     i          i          i           i          o
!
!  Initialize the transmitted ions statistics
!
      if( out%lhisaat )  &
            call initHis1d( indhisProj, out%nboxa, out%wboxa, ZERO, hisaa1t )
!                =========      i           i          i       i       o
      if( out%lhisat )  &
            call initHis1d( indhisProj, 2*out%nboxa, out%wboxa, ZERO, hisa1t )
!                =========      i            i           i       i      o
      if( out%lhisa2t )  &
            call initHis2d( indhisProj, out%nboxa2, out%wboxa2, [ZERO, ZERO], hisa2t )
!                =========      i           i           i            i          o
      if( out%lhiset )  &
            call initHis1d( indhisProj, out%nboxe, out%wboxe, ZERO, hise1t )
!                =========      i           i          i       i      o
      if( out%lhist )  &
            call initHis1d( indhisIn, out%nbox, out%wbox, ZERO, his1t )
!                =========     i         i         i       i      o
      if( out%lhis2t )  &
            call initHis2d( indhisIn, out%nbox2, out%wbox2, [ZERO, ZERO], his2t )
!                =========     i          i          i           i          o
!
!  Deallocate index arrays
!
      deallocate( indhisIn )
      deallocate( indhisRec )
      deallocate( indhisRec3 )
      deallocate( indhisProj )
!
!  Load the backup or histogram file if requested
!
!$omp master
      if( usebk ) then
         call stat50( )
!             ======
      else if( usehis .and. .not.ldyn ) then
         call readHis1d( hisfile, his1, atom1, atom2, factor_=1e-24_REAL8 )
!             =========     i      o      i      i           i
      end if
!
!  Initialize deterministic defect structure
!
      if( dam%ldet )  call det10( his1 )
!                          =====   io
!
!  Initialize histogram according to defect file if statistic model is to be used
!
      if( .not.dam%ldet .and. usedef .and. out%lhis )  call his10( his1 )
!                                                           =====   io
!$omp end master
!$omp barrier
!
!  In case of dynamic damage consideration, spread histogram to all threads
!
      if( dam%ndamdim == 1 .and. out%lhis ) then
         call broadcastHis1d( his1 )
!             ==============   io
         if( dam%lamzon )  &
               call setHisAmzon( his1, hisa1 )
!                   ===========   i      o
      else if( dam%ndamdim == 2 .and. out%lhis2 ) then
         call broadcastHis2d( his2 )
!             ==============   io
      else if( dam%ndamdim == 3 .and. out%lhis3 ) then
         call broadcastHis3d( his3 )
!             ==============   io
      end if

   end subroutine stat10
!
!=======================================================================
!
   subroutine stat11( )
!             ======
!
!     Initialize statistic quantities at the beginning of the simulation of a 
!     new ion molecule (not including those done at the beginning of the 
!     simulation of each ion atom -> stat12 )
!
!     Currently no action
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
      continue

   end subroutine stat11
!
!=======================================================================
!
   subroutine stat12( pos )
!             ======   i
!
!     Initialize statistic quantities at the beginning of the simulation of a 
!     new ion atom
!
!  Arguments:
!
      real(REAL8),  intent(in) :: pos(3)      ! starting point of trajectory
!
!  Local Variables:
!
      real(REAL8) :: posdum(3)
      real(REAL8) :: dir_s(3)
!
!  Start of Executable Code:
!  ------------------------
!
!  Initialize the statistics coordinate system
!
      select case( out%mcoord )
         case( 2 )
            call getXStart( posdum, dir_s )
!                =========    o       o
            call initCoordSys( coord_statis, dir_s )
!                ============       o          i
         case( 5 )
            if( mat(irxtal)%lpoly )  coord_statis = coord_xtal
      end select

   end subroutine stat12
!
!=======================================================================
!
   subroutine stat22( proj, dee_, recoil_ )
!             ======   i     i       i
!
!     Score statistic quantities along the trajectory
!
!     All passed variables are those after ffp and collision
!
!  Arguments:
!
      type(Proj_type), intent(in) :: proj             ! projectile
      real(REAL8), optional, intent(in) :: dee_       ! electronic energy loss during
                                                      !    preceeding ffp [eV] 
      type(Proj_type), optional, intent(in) :: recoil_! recoil
!
!  Local Variables:
!
      logical :: allow_score
      integer(INT4) :: ielemp, ielema
      real(REAL8) :: x, y, z, pos_s(3)
!
!  Start of Executable Code:
!  ------------------------
!
      ielemp = proj%ielem
      if( ldyn .and. proj%origin%internal )  ielemp = nelemIon + ielemp
!
!  Increase the electronic energy sum
!
      if( present( dee_ ) .and. .not.proj%virtual ) then

         pos_s(:) = proj%pos(:)
         x = pos_s(1)
         y = pos_s(2)
         z = pos_s(3)

         call applyBndCond3d( x, y, z, .false., allow_score )
!             ============== io io io     i          o

         if( .not.allow_score )  return

         call scoreEED( x, y, z, ielemp, dee_ )
!             ========  i  i  i    i      i

      end if
!
!  Increase the nuclear energy sum
! 
      if( present( recoil_ ) ) then
         ielema = recoil_%ielem
         if( ldyn .and. recoil_%origin%internal )  ielema = nelemIon + ielema
         sumden(ielema,ielemp) = sumden(ielema,ielemp) + recoil_%e * ion%ddose
      end if

   end subroutine stat22
!
!=======================================================================
!
   subroutine stat2x( proj, ddami, ddamv, edam, iflag )
!             ======   i      i      i     i      i
!
!     Score statistic quantities at the beginning of a trajectory, when a recoil is not
!     followed, or at the end of a trajectory
!
!  Arguments:
!
      type(Proj_type), intent(in) :: proj       ! projectile
      real(REAL8), intent(in) :: ddami          ! weight of projectile/interstitial
      real(REAL8), intent(in) :: ddamv          ! weight of vacancy
      real(REAL8), intent(in) :: edam           ! damage energy
      integer(INT4), intent(in) :: iflag        ! flag indicating meaning of projectile
                                                !    0 ... recoil without trajectory
                                                !    1 ... initial state of recoil
!                                               !    3 ... final state of projectile
!  Local Variables:
!
      logical :: allow_score, scorei, scorev, scoreee, scorene
      integer(INT4) :: ielemp, ielemv, i
      integer(INT4) :: iamo
      real(REAL8) :: value2(2), value3(3)
      real(REAL8) :: x, y, z, zo, d, p(3), vec(3)
      real(REAL8) :: pos_s(3), poso_s(3)
      real(REAL8) :: z_left, z_right, dz
!
!  Start of Executable Code:
!  ------------------------
!
!  Element index of projectile
!
      ielemp = proj%ielem
      if( ldyn .and. proj%origin%internal )  ielemp = nelemIon + ielemp
      !if( ielemp == 1 )  write( 998, '(a,i1,a,i2,l2)' ) 'stat2', iflag, ', ig,internal=',  &
      !                         proj%ig, proj%origin%internal
!
!  Original and current position of projectile
!
      pos_s(:) = proj%pos(:)
      poso_s(:) = proj%origin%pos(:)
!                                   ! revert shift due to time integral and lattice vibrations
      if( iflag == 0 .or. iflag == 1 )  pos_s(:) = poso_s(:)
!
!  First treat scores which have to be performed at the origin of the trajectory
!
      zo = poso_s(3)

      call applyBndCond1d( zo,  .false., allow_score )
!          ==============  io      i          o

      if( allow_score .and. out%lhism .and. ielemp > nelemIon ) then
         vec(:) = pos_s(:) - poso_s(:)
         do i = 1, 3
            call scoreHis1d( 3*(ielemp-1)+i, zo, hism1, vec(i)*ion%ddose/ion%dose )
!                ==========         i        i    io               i
         end do
      end if
!
!  Now treat scores of current position...
!
      x = pos_s(1)
      y = pos_s(2)
      z = pos_s(3)

      call applyBndCond3d( x, y, z, .false., allow_score )
!          ============== io io io     i          o

      if( .not.allow_score )  return
!
!  Determine what to score
!
      scorei = ( proj%ielem <= nelemIn .and. ddami /= ZERO )
      scorev = ( proj%ielem <= nelemIn .and. ddamv /= ZERO .and. dam%lrcoil )
      scoreee = ( iflag /= 1 .and. .not.proj%virtual )
      scorene = ( edam > 0 )
!
!  Determine where to store
!
      if( scorev ) then
         ielemv = nelemTarg + proj%ielem
         if( ldyn )  ielemv = nelemIon + ielemv
      end if
!
!  Count the ions stopping inside the target
! 
      if( scorei )  rin(ielemp) = rin(ielemp) + ddami
      if( scorev )  rin(ielemv) = rin(ielemv) + ddamv
!
!  Evaluate the global moments
! 
      if( out%lmom ) then
         if( scorei ) then
            call scoreMom1d( momz(ielemp), z, ddami*ion%ddose )
!                ==========       io       i        i
            if( out%nlatmom > 0 ) then
               call scoreMom1d( momx(ielemp), x, ddami*ion%ddose )
!                   ==========       io       i        i
               call scoreMom1d( momy(ielemp), y, ddami*ion%ddose )
!                   ==========       io       i        i
            end if
            if( out%nmixmom > 0 ) then
               call incrMom2d( momxy(ielemp), x, y, ddami*ion%ddose )
!                   =========       io        i  i      i
               call incrMom2d( momzx(ielemp), z, x, ddami*ion%ddose )
!                   =========       io        i  i      i
               call incrMom2d( momzy(ielemp), z, y, ddami*ion%ddose )
!                   =========       io        i  i      i
            end if
         end if
         if( scorev ) then
            call scoreMom1d( momz(ielemv), z, ddamv*ion%ddose )
!                ==========       io       i        i
            if( out%nlatmom > 0 ) then
               call scoreMom1d( momx(ielemv), x, ddamv*ion%ddose )
!                   ==========       io       i        i
               call scoreMom1d( momy(ielemv), y, ddamv*ion%ddose )
!                   ==========       io       i        i
            end if
            if( out%nmixmom > 0 ) then
               call incrMom2d( momxy(ielemv), x, y, ddamv*ion%ddose )
!                   =========       io        i  i        i
               call incrMom2d( momzx(ielemv), z, x, ddamv*ion%ddose )
!                   =========       io        i  i        i
               call incrMom2d( momzy(ielemv), z, y, ddamv*ion%ddose )
!                   =========       io        i  i        i
            end if
         end if
      end if
!
!  Evaluate the 1d atom distribution histogram
! 
      if( out%lhis ) then
         if( nionhis > 0 .and. dam%ndamdim == 1 ) then
            if( scorei )  call scoreHis1d( ielemp, z, his1_tmp, ddami*out%fachis(1) )
!                              ==========    i     i     io           i
            if( scorev )  call scoreHis1d( ielemv, z, his1_tmp, ddamv*out%fachis(1) )
!                              ==========    i     i     io           i
         else
            if( scorei )  call scoreHis1d( ielemp, z, his1, ddami*out%fachis(1) )
!                              ==========    i     i   io         i
            if( scorev )  call scoreHis1d( ielemv, z, his1, ddamv*out%fachis(1) )
!                              ==========    i     i   io         i
         end if
      end if
!
!  Evaluate 2d atom distribution histogram
!
      if( out%lhis2 ) then
         if( out%lcyl ) then
            value2(1) = sqrt( ( x-out%xcyl )**2 + ( y-out%ycyl )**2 )
         else
            value2(1) = x
         end if
         value2(2) = z
         if( nionhis > 0 .and. dam%ndamdim == 2 ) then
            if( scorei )  call scoreHis2d( ielemp, value2, his2_tmp, ddami*out%fachis(2) )
!                              ==========    i       i        io           i
            if( scorev )  call scoreHis2d( ielemv, value2, his2_tmp, ddamv*out%fachis(2) )
!                              ==========    i       i        io           i
         else
            if( scorei )  call scoreHis2d( ielemp, value2, his2, ddami*out%fachis(2) )
!                              ==========    i       i      io         i
            if( scorev )  call scoreHis2d( ielemv, value2, his2, ddamv*out%fachis(2) )
!                              ==========    i       i      io         i
         end if
      end if
!
!  Evaluate 3d atom distribution histogram
!
      if( out%lhis3 ) then
         value3(1:3) = [ x, y, z ]
         if( nionhis > 0 .and. dam%ndamdim == 3 ) then
            if( scorei )  call scoreHis3d( ielemp, value3, his3_tmp, ddami*out%fachis(3) )
!                              ==========    i       i        io           i
            if( scorev )  call scoreHis3d( ielemv, value3, his3_tmp, ddamv*out%fachis(3) )
!                              ==========    i       i        io           i
         else
            if( scorei )  call scoreHis3d( ielemp, value3, his3, ddami*out%fachis(3) )
!                              ==========    i       i      io         i
            if( scorev )  call scoreHis3d( ielemv, value3, his3, ddamv*out%fachis(3) )
!                              ==========    i       i      io         i
         end if
      end if
! 
!  Evaluate electronic energy statistics
!
      if( scoreee ) then
         !sumdee(ielemp) = sumdee(ielemp) + (proj%e - edam) * ion%ddose
         call scoreEED( x, y, z, ielemp, proj%e-edam )
!             ========  i  i  i    i          i
      end if
!
!  Evaluate the displacement histogram
!
      if( out%lhisd ) then
         d = abs( pos_s(:) - poso_s(:) )
!            ===
         if( d > ZERO ) then
            call scoreHis1d( ielemp, d, hisd1, ion%ddose/ion%dose, power_=out%dpower )
!                ==========    i     i   io            i                  i
         end if
      end if

      if( out%lhism .and. ielemp <= nelemIon ) then   ! ielemp > nelemIon already treated
         vec(:) = pos_s(:) - poso_s(:)
         do i = 1, 3
            call scoreHis1d( 3*(ielemp-1)+i, z, hism1, vec(i)*ion%ddose/ion%dose )
!                ==========         i        i   io               i
         end do
      end if

      if( iflag /= 1 .and. out%lhisp ) then
         p = sqrt( TWO * elem(proj%ielem)%m * proj%e ) * proj%dir(:)
         do i = 1, 3
            call scoreHis1d( 3*(ielemp-1)+i, z, hisp1, p(i)*ion%ddose/ion%dose )
!                ==========         i        i   io             i
         end do
      end if
!
!  If there is deposited nuclear energy to be recorded...
!  ------------------------------------------------------
!
      if( scorene ) then

         call scoreNED( x, y, z, ielemp, edam )
!             ========  i  i  i    i      i
!
!  Move the amorphous zone boundaries according to the IBIIA effect
!
         if( dam%libiia ) then

            do iamo = 1, namo
               z_left = max( geom%posif(irxtal-1), zamo_left(iamo) - dam%dibiia/2 )
               z_right = min( geom%posif(irxtal), zamo_left(iamo) + dam%dibiia/2 )
               if( z < z_left .or. z > z_right )  cycle
               dz = dam%ribiia * edam * ion%ddose / (z_right - z_left)
               if( nionhis <= 0 ) then
!$omp atomic update
                  zamo_left(iamo) = zamo_left(iamo) - dz
               else
                  dzamo_left(iamo) = dzamo_left(iamo) - dz
               end if
            end do

            do iamo = 1, namo
               z_left = max( geom%posif(irxtal-1), zamo_right(iamo) - dam%dibiia/2 )
               z_right = min( geom%posif(irxtal), zamo_right(iamo) + dam%dibiia/2 )
               if( z < z_left .or. z > z_right )  cycle
               dz = dam%ribiia * edam * ion%ddose / (z_right - z_left)
               if( nionhis <= 0 ) then
!$omp atomic update
                  zamo_right(iamo) = zamo_right(iamo) + dz
               else
                  dzamo_right(iamo) = dzamo_right(iamo) + dz
               end if
            end do

         end if

      end if

   end subroutine stat2x
!
!=======================================================================
!
   subroutine stat24( proj )
!             ======    i
!
!     Score statistic quantities of a backscattered atom
!
!  Arguments:
!
      type(Proj_type), intent(in) :: proj    ! projectile
!
!  Local Variables:
!
      logical :: allow_score, scorev
      integer(INT4) :: ielem, ielemv
      real(REAL8) :: x,  y,  z,  xo, yo, zo, dirx,  diry,  dirz,  angle, anglea,  e,  fac
      real(REAL8) :: xp,  yp, xmax
      real(REAL8) :: pos_s(3), dir_s(3), poso_s(3)
!
!  Start of Executable Code:
!  ------------------------
!
!  First check if projectile fulfills analyzer requirements
!
      if( out%mcoord == 1 ) then
         dir_s(:) = proj%dir(:)
      else
         call toCoordSys( coord_statis, proj%dir, dir_s )
!             ==========       i           i        o
      end if
      dirx = dir_s(1)
      diry = dir_s(2)
      dirz = dir_s(3)

      angle = 180/PI * atan2( dirx, -dirz )     ! angle of 2D projection in coord_statis

      anglea = 180/PI * acos( -dirz )           ! angle with respect to analyzer axis

      e = proj%e

      if( anglea < out%amina .or. anglea > out%amaxa .or.  &
          max( ZERO, e ) < out%emina .or. e > out%emaxa )  return

      if( geom%lrobin ) then                    ! no boundary conditions for lrobin=t
         fac = ( -proj%pos(3) ) / dirz
         xp = proj%pos(1) + fac*dirx
         yp = proj%pos(2) + fac*diry
         xmax = 13.575
         if( anglea > 80 .or. xp < -xmax .or. xp > xmax .or. yp < -xmax .or. yp > xmax )  return
      end if
!
!  Treat scores of the vacancy position...
!
      scorev = ( proj%origin%internal .and. dam%lrcoil )

      if( scorev ) then
         poso_s(:) = proj%origin%pos(:)
         xo = poso_s(1)
         yo = poso_s(2)
         zo = poso_s(3)
         call applyBndCond3d( xo, yo, zo, .false., allow_score )
   !          ==============  io  io  io     i          o
         scorev = ( scorev .and. allow_score )
      end if

      if( scorev ) then

         ielemv = nelemTarg + proj%ielem
         if( ldyn )  ielemv = nelemIon + ielemv
!
!  Evaluate backscattered atoms vacancy moments
!
         if( out%lmomb ) then
            if( out%nlatmom > 0 ) then
               call scoreMom1d( momxb(ielemv), xo, ion%ddose )
!                   ==========       io        i      i
               call scoreMom1d( momyb(ielemv), yo, ion%ddose )
!                   ==========       io        i      i
            end if
            if( out%nmixmom > 0 )  call incrMom2d( momxyb(ielemv), xo, yo, ion%ddose )
!                                       =========        io        i   i      i
         end if
!
!  Evaluate backscattered atoms vacancy histogram
!
         if( out%lhisb )  &
               call scoreHis1d( ielemv, xo, his1b, ion%ddose/ion%dose )
!                   ==========    i     i    io           i
         if( out%lhis2b )  &
               call scoreHis2d( ielemv, [xo, yo], his2b, ion%ddose/ion%dose )
!                   ==========    i        i       io            i

      end if
!
!  Treat scores at current position...
!
      pos_s(:) = proj%pos(:)

      x = pos_s(1)
      y = pos_s(2)
      z = pos_s(3)

      call applyBndCond3d( x, y, z, .false., allow_score )
!          ============== io io io     i          o

      if( .not.allow_score )  return

      ielem = proj%ielem
      if( ldyn .and. proj%origin%internal )  ielem = nelemIon + ielem
      !if( ielem == 1 )  write( 998, * ) 'stat24, ig, internal=', proj%ig, proj%origin%internal
!
!  Count backscattered atoms
! 
      nback(ielem) = nback(ielem) + 1
!
!  Add energy loss sums
!
      sumeb(ielem) = sumeb(ielem) + e*ion%ddose
! 
!  Evaluate backscattered atoms moments
!
      if( out%lmomaab )  &
            call scoreMom1d( momaab(ielem), anglea, ion%ddose )
!                ==========       io          i         i
      if( out%lmomab )  &
            call scoreMom1d( momab(ielem), angle, ion%ddose )
!                ==========      io          i        i
      if( out%lmomeb )  &
            call scoreMom1d( momeb(ielem), e/1000, ion%ddose )
!                ==========      io          i         i
      if( out%lmomb ) then
         if( out%nlatmom > 0 ) then
            call scoreMom1d( momxb(ielem), x, ion%ddose )
!                ==========     io         i      i
            call scoreMom1d( momyb(ielem), y, ion%ddose )
!                ==========     io         i      i
         end if
         if( out%nmixmom > 0 )  call incrMom2d( momxyb(ielem), x, y, ion%ddose )
!                                    =========      io         i  i      i
      end if
!
!  Evaluate backscattered atoms histogram
!
      if( out%lhisaab )  &
            call scoreHis1d( ielem, anglea, hisaa1b, ion%ddose/ion%dose )
!                ==========    i      i       io            i
      if( out%lhisab )  &
            call scoreHis1d( ielem, angle, hisa1b, ion%ddose/ion%dose )
!                ==========    i      i      io           i
      if( out%lhisa2b )  &
            call scoreHis2d( ielem, [dirx, diry], hisa2b, ion%ddose/ion%dose )
!                ==========    i         i          io           i
      if( out%lhiseb )  &
            call scoreHis1d( ielem, e, hise1b, ion%ddose/ion%dose, power_=out%epower )
!                ==========    i    i    io           i                   i
      if( out%lhisb ) then
         call scoreHis1d( ielem, x, his1b, ion%ddose/ion%dose )
!             ==========    i    i   io            i
      end if
      if( out%lhis2b ) then
         call scoreHis2d( ielem, [x, y], his2b, ion%ddose/ion%dose )
!             ==========    i      i      io            i
      end if

   end subroutine stat24
!
!=======================================================================
!
   subroutine stat25( proj )
!             ======   i
!
!     Score statistic quantities of a transmitted atom
!
!  Arguments:
!
      type(Proj_type), intent(in) :: proj    ! projectile
!
!  Local Variables:
!
      logical :: lnot, allow_score, scorev
      integer(INT4) :: ielem, ielemv
      real(REAL8) :: x, y, z, xo, yo, zo, dirx, diry, dirz, e, angle, anglea
      real(REAL8) :: pos_s(3), dir_s(3), poso_s(3)
!
!  Start of Executable Code:
!  ------------------------
!
!  First check if projectile fulfills analyzer requirements
!
      if( out%mcoord == 1 ) then
         dir_s(:) = proj%dir(:)
      else
         call toCoordSys( coord_statis, proj%dir, dir_s )
!             ==========       i           i        o
      end if
      dirx = dir_s(1)
      diry = dir_s(2)
      dirz = dir_s(3)

      angle = 180/PI * atan2( dirx, dirz )      ! angle of 2D projection in coord_statis

      anglea = 180/PI * acos( dirz )            ! angle with respect to analyzer axis

      e = proj%e

      if( anglea < out%amina .or. anglea > out%amaxa .or.  &
          max( ZERO, e ) < out%emina .or. e > out%emaxa )  return

      if( geom%lrobin ) then                    ! no boundary conditions for lrobin=t
         call rob25( e, proj%pos, dir_s, lnot )
!             =====  i     i        i     o
         if( lnot ) return
      end if
!
!  Treat scores of the vacancy position...
!
      scorev = ( proj%origin%internal .and. dam%lrcoil )

      if( scorev ) then
         poso_s(:) = proj%origin%pos(:)
         xo = poso_s(1)
         yo = poso_s(2)
         zo = poso_s(3)
         call applyBndCond3d( xo, yo, zo, .false., allow_score )
   !          ==============  io  io  io     i          o
         scorev = ( scorev .and. allow_score )
      end if

      if( scorev ) then

         ielemv = nelemTarg + proj%ielem
         if( ldyn )  ielemv = nelemIon + ielemv
!
!  Evaluate transmitted atoms vacancy moments
!
         if( out%lmomt ) then
            if( out%nlatmom > 0 ) then
               call scoreMom1d( momxt(ielemv), xo, ion%ddose )
!                   ==========       io        i      i
               call scoreMom1d( momyt(ielemv), yo, ion%ddose )
!                   ==========       io        i      i
            end if
            if( out%nmixmom > 0 )  call incrMom2d( momxyt(ielemv), xo, yo, ion%ddose )
!                                       =========        io        i   i       i
         end if
!
!  Evaluate transmitted atoms histogram
!
         if( out%lhist )  &
               call scoreHis1d( ielemv, xo, his1t, ion%ddose/ion%dose )
!                   ==========    i     i    io           i
         if( out%lhis2t )  &
               call scoreHis2d( ielemv, [xo, yo], his2t, ion%ddose/ion%dose )
!                   ==========    i        i       io            i

      end if
!
!  Treat scores at current position...
!
      pos_s(:) = proj%pos(:)

      x = pos_s(1)
      y = pos_s(2)
      z = pos_s(3)

      call applyBndCond3d( x, y, z, .false., allow_score )
!          ============== io io io     i          o

      if( .not.allow_score )  return

      ielem = proj%ielem
      if( ldyn .and. proj%origin%internal )  ielem = nelemIon + ielem
! 
!  Count transmitted atom
! 
      ntrans(ielem) = ntrans(ielem) + 1
! 
!  Add energy loss sums
!
      sumet(ielem) = sumet(ielem) + e*ion%ddose
! 
!  Evaluate transmitted atoms moments
!
      if( out%lmomaat )  &
            call scoreMom1d( momaat(ielem), anglea, ion%ddose )
!                ==========       io          i         i
      if( out%lmomat )  &
            call scoreMom1d( momat(ielem), angle, ion%ddose )
!                ==========       io         i        i
      if( out%lhisa2t )  &
            call scoreHis2d( ielem, [dirx, diry], hisa2t, ion%ddose/ion%dose )
!                ===========    i         i         io           i
      if( out%lmomet )  &
            call scoreMom1d( momet(ielem), e/1000, ion%ddose )
!                ==========       io         i         i
      if( out%lmomt ) then
         if( out%nlatmom > 0 ) then
            call scoreMom1d( momxt(ielem), x, ion%ddose )
!                ==========      io       i      i
            call scoreMom1d( momyt(ielem), y, ion%ddose )
!                ==========      io        i      i
         end if
         if( out%nmixmom > 0 )  call incrMom2d( momxyt(ielem), x, y, ion%ddose )
!                                    =========      io         i  i      i
      end if
!
!  Evaluate transmitted atoms histogram
!
      if( out%lhisaat )  &
            call scoreHis1d( ielem, anglea, hisaa1t, ion%ddose/ion%dose )
!                ==========    i      i       io            i
      if( out%lhisat )  &
            call scoreHis1d( ielem, angle, hisa1t, ion%ddose/ion%dose )
!                ==========    i      i      io           i
      if( out%lhiset )  &
            call scoreHis1d( ielem, e, hise1t, ion%ddose/ion%dose, power_=out%epower )
!                ==========    i    i    io           i                    i
      if( out%lhist ) then
         call scoreHis1d( ielem, x, his1t, ion%ddose/ion%dose )
!             ==========    i    i   io           i
         if( scorev )  call scoreHis1d( ielemv, xo, his1t, ion%ddose/ion%dose )
!                           ==========    i     i    io           i
      end if
      if( out%lhis2t )  &
            call scoreHis2d( ielem, [x, y], his2t, ion%ddose/ion%dose )
!                ==========    i      i      io            i

   end subroutine stat25
!
!=======================================================================
!
   subroutine scoreEED( x, y, z, ielemp, dee )
!             ========  i  i  i    i      i
!
!     Score electronic energy deposition
!
!  Arguments:
!
      real(REAL8), intent(in) :: x              ! x coordinate of electronic energy deposition
      real(REAL8), intent(in) :: y              ! y coordinate of electronic energy deposition
      real(REAL8), intent(in) :: z              ! z coordinate of electronic energy deposition
      integer(INT4), intent(in) :: ielemp       ! atom species of the projectile
      real(REAL8), intent(in) :: dee            ! electronic energy deposition
!
!  Local Variables:
!
      real(REAL8) :: value2(2), value3(3)
!
!  Start of Executable Code:
!  ------------------------
!
!
!  Increase energy deposition sum
!
      sumdee(ielemp) = sumdee(ielemp) + dee * ion%ddose
!
!  Evaluate global electronic energy deposition moments
!
      if( out%lmomee ) then
         call scoreMom1d( momzee(ielemp), z, dee*ion%ddose )
!             ==========       io         i        i
         if( out%nlatmom > 0 ) then
            call scoreMom1d( momxee(ielemp), x, dee*ion%ddose )
!                ==========       io         i        i
            call scoreMom1d( momyee(ielemp), y, dee*ion%ddose )
!                ==========       io        i        i
         end if
         if( out%nmixmom > 0 ) then
            call incrMom2d( momxyee(ielemp), x, y, dee*ion%ddose )
!                =========        io         i  i        i
            call incrMom2d( momzxee(ielemp), z, x, dee*ion%ddose )
!                =========        io         i  i        i
            call incrMom2d( momzyee(ielemp), z, y, dee*ion%ddose )
!                =========        io         i  i        i
         end if
      end if
!
!  Evaluate 1d deposited electronic energy histogram
!
      if( out%lhisee )  call scoreHis1d( ielemp, z, hisee1, dee*out%fachis(1) )
!                            ==========    i     i    io          i
!
!  Evaluate 2d deposited electronic energy histogram
!
      if( out%lhisee2 ) then
         if( out%lcyl ) then
            value2(1) = sqrt( ( x-out%xcyl )**2 + ( y-out%ycyl )**2 )
         else
            value2(1) = x
         end if
         value2(2) = z
         call scoreHis2d( ielemp, value2, hisee2, dee*out%fachis(2) )
!             ==========    i       i       io         i
      end if
!
!  Evaluate 3d deposited electronic energy histogram
!
      if( out%lhisee3 ) then
         value3(1:3) = [ x, y, z ]
         call scoreHis3d( ielemp, value3, hisee3, dee*out%fachis(2) )
!             ==========    i       i       io         i
      end if

   end subroutine scoreEED
!
!=======================================================================
!
   subroutine scoreNED( x, y, z, ielemp, edam )
!             ========  i  i  i    i      i
!
!     Score nuclear energy deposition
!
!  Arguments:
!
      real(REAL8), intent(in) :: x              ! x coordinate of nuclear energy deposition
      real(REAL8), intent(in) :: y              ! y coordinate of nuclear energy deposition
      real(REAL8), intent(in) :: z              ! z coordinate of nuclear energy deposition
      integer(INT4), intent(in) :: ielemp       ! atom species of the projectile
      real(REAL8), intent(in) :: edam           ! damage energy (nuclear energy deposition)
!
!  Local Variables:
!
      real(REAL8) :: value2(2), value3(3)
!
!  Start of Executable Code:
!  ------------------------
!
!
!  Increase energy deposition sum
!
      sumedam(ielemp) = sumedam(ielemp) + edam * ion%ddose

      redam(ielemp) = redam(ielemp) + edam
!
!  Evaluate global nuclear energy deposition moments
!
      if( out%lmomne ) then
         call scoreMom1d( momzne(ielemp), z, edam*ion%ddose )
!             ==========       io         i        i
         if( out%nlatmom > 0 ) then
            call scoreMom1d( momxne(ielemp), x, edam*ion%ddose )
!                ==========       io         i        i
            call scoreMom1d( momyne(ielemp), y, edam*ion%ddose )
!                ==========       io        i        i
         end if
         if( out%nmixmom > 0 ) then
            call incrMom2d( momxyne(ielemp), x, y, edam*ion%ddose )
!                =========        io         i  i        i
            call incrMom2d( momzxne(ielemp), z, x, edam*ion%ddose )
!                =========        io         i  i        i
            call incrMom2d( momzyne(ielemp), z, y, edam*ion%ddose )
!                =========        io         i  i        i
         end if
      end if
!
!  Evaluate 1d deposited nuclear energy histogram
!
      if( out%lhisne )  call scoreHis1d( ielemp, z, hisne1, edam*out%fachis(1) )
!                            ==========    i     i     io           i
!
!  Evaluate 2d deposited nuclear energy histogram
!
      if( out%lhisne2 ) then
         if( out%lcyl ) then
            value2(1) = sqrt( ( x-out%xcyl )**2 + ( y-out%ycyl )**2 )
         else
            value2(1) = x
         end if
         value2(2) = z
         call scoreHis2d( ielemp, value2, hisne2, edam*out%fachis(2) )
!             ==========    i       i       io         i
      end if
!
!  Evaluate 3d deposited nuclear energy histogram
!
      if( out%lhisne3 ) then
         value3(1:3) = [ x, y, z ]
         call scoreHis3d( ielemp, value3, hisne3, edam*out%fachis(3) )
!             ==========    i       i       io         i
      end if

   end subroutine scoreNED
!
!=======================================================================
!
   subroutine stat26( )
!             ======
!
!     Perform actions after each ion molecule, which operate only on the own image
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
!  Update moments
!
      call updateMoments( )
!          =============

   end subroutine stat26
!
!=======================================================================
!
   subroutine updateMoments( )
!             =============
!
!     Update moments after a set of dependent scores has been completed
!
!  Arguments:
!
!     None
!
!  Local Variables:
!
      integer(INT4) :: ielem
!
!  Start of Executable Code:
!  ------------------------
!
!  Increment yield moments
!
      do ielem = 1, nelemIn
         call scoreMom1d( momyldin(ielem), rin(ielem), ion%ddose )
!             ==========        io             i           i
         call updateMom1d( momyldin(ielem) )
!             ===========        io
      end do
      do ielem = 1, nelemAtom
         call scoreMom1d( momylde(ielem), redam(ielem), ion%ddose )
!             ==========       io              i            i
         call updateMom1d( momylde(ielem) )
!             ===========       io
      end do
      do ielem = 1, nelemProj
         call scoreMom1d( momyldb(ielem), real(nback(ielem),REAL8), ion%ddose )
!             ==========       io                   i                   i
         call updateMom1d( momyldb(ielem) )
!             ===========       io
      end do
      do ielem = 1, nelemProj
         call scoreMom1d( momyldt(ielem), real(ntrans(ielem),REAL8), ion%ddose )
!             ==========       io                    i                   i
         call updateMom1d( momyldt(ielem) )
!             ===========       io
      end do

      rin(:) = ZERO
      redam(:) = ZERO
      nback(:) = 0
      ntrans(:) = 0
!
!  Evaluate the global moments
!
      if( out%lmom ) then
         do ielem = 1, nelemIn
            call updateMom1d( momz(ielem) )
!                ===========      io
            if( out%nlatmom > 0 ) then
               call updateMom1d( momx(ielem) )
!                   ===========       io
               call updateMom1d( momy(ielem) )
!                   ===========       io
            end if
         end do
      end if
!
!  Evaluate global nuclear energy deposition moments
!
      if( out%lmomee ) then
         do ielem = 1, nelemAtom
            call updateMom1d( momzee(ielem) )
!                ===========       io
            if( out%nlatmom > 0 ) then
               call updateMom1d( momxee(ielem) )
!                   ===========       io
               call updateMom1d( momyee(ielem) )
!                   ===========       io
            end if
         end do
      end if
!
!  Evaluate global nuclear energy deposition moments
!
      if( out%lmomne ) then
         do ielem = 1, nelemAtom
            call updateMom1d( momzne(ielem) )
!                ===========       io
            if( out%nlatmom > 0 ) then
               call updateMom1d( momxne(ielem) )
!                   ===========       io
               call updateMom1d( momyne(ielem) )
!                   ===========       io
            end if
         end do
      end if
!
!  Evaluate backscattered atoms moments
!
      if( out%lmomaab ) then
         do ielem = 1, nelemProj
            call updateMom1d( momaab(ielem) )
!                ===========       io
         end do
      end if
      if( out%lmomab ) then
         do ielem = 1, nelemProj
            call updateMom1d( momab(ielem) )
!                ===========      io
         end do
      end if
      if( out%lmomeb ) then
         do ielem = 1, nelemProj
            call updateMom1d( momeb(ielem) )
!                ===========      io
         end do
      end if
      if( out%lmomb ) then
         if( out%nlatmom > 0 ) then
            do ielem = 1, nelemIn
               call updateMom1d( momxb(ielem) )
!                   ===========     io
               call updateMom1d( momyb(ielem) )
!                   ===========     io
            end do
         end if
      end if
!
!  Evaluate transmitted atoms moments
!
      if( out%lmomaat ) then
         do ielem = 1, nelemProj
            call updateMom1d( momaat(ielem) )
!                ===========       io
         end do
      end if
      if( out%lmomat ) then
         do ielem = 1, nelemProj
            call updateMom1d( momat(ielem) )
!                ===========       io
         end do
      end if
      if( out%lmomet ) then
         do ielem = 1, nelemProj
            call updateMom1d( momet(ielem) )
!                ===========       io
         end do
      end if
      if( out%lmomt ) then
         if( out%nlatmom > 0 ) then
            do ielem = 1, nelemIn
               call updateMom1d( momxt(ielem) )
!                   ===========      io
               call updateMom1d( momyt(ielem) )
!                   ===========      io
            end do
         end if
      end if

   end subroutine updateMoments
!
!=======================================================================
!
   subroutine stat27( iion, final_ )
!             ======   i      i
!
!     Perform actions after a certain number of ion molecules, which operate on
!     coindexed variables.
!
!     The final evaluation (final_ argument present and .true.) serves to gather remaining data
!     between regular updates and the end of the simulation. For the final evaluation iion must
!     equal nion.
!
!     This subroutine must be called from a parallel region.
!
!  Arguments:
!
      integer(INT4), intent(in) :: iion         ! index of ion
      logical, intent(in), optional :: final_   ! indicates final call of stat27
!
!  Local Variables:
!
      logical :: final, lupdate
!
!  Start of Executable Code:
!  ------------------------
!
      if( present( final_ ) ) then
         final = final_
      else
         final = .false.
      end if
!
!  Update damage histograms
!
      if( nionhis <= 0 .or. .not.ldamdyn ) then
         lupdate = .false.
      else
         lupdate = ( mod( iion, nionhis ) == 0 )
         if( final )  lupdate = .not.lupdate
      end if

      if( lupdate ) then
!$omp barrier
         call updateDamageHistogram( )
!             =====================
      end if
!
!  Update amorphous zones (also call this when nionhis<=0)
!
      if( dam%libiia ) then
         lupdate = ( lupdate .or. (nionhis <= 0 .and. .not.final) )
         if( lupdate ) then
!$omp barrier
            call updateAmo( )
!                =========
         end if
      end if

   end subroutine stat27
!
!=======================================================================
!
   subroutine stat2e( proj )
!             ======   i
!
!     Evaluate at the end of the trajectory if ion is lost
!
!  Arguments:
!
      type(Proj_type), intent(in) :: proj       ! projectile
!
!  Local Variables:
!
      integer(INT4) :: ielem
!
!  Start of Executable Code:
!  ------------------------
!
! debug
      if( nthread == 1 ) then
         print *, 'lost atom for ion = ', proj%iion
!$    else
!$       print *, 'lost atom for ion = ', proj%iion, ' in thread ', omp_get_thread_num( )
      end if
!      error stop
! end debug

      ielem = proj%ielem

      nlost(ielem) = nlost(ielem) + 1

   end subroutine stat2e
!
!=======================================================================
!
   subroutine updateDamageHistogram( )
!             =====================
!
!     Update damage histograms his* from temporary histograms his*_tmp,
!     and reinitialize his*_tmp
!
!     This subroutine must be called from within a parallel region
!
!     The threads are "almost synchronized" when returning from this subroutine
!
!  Arguments:
!
!     None
!
!  Local Variables:
!
      integer(INT4), allocatable :: indhis(:)
      real(REAL8) :: hbeg, hbeg2(2), hbeg3(3)
!
!  Start of Executable Code:
!  ------------------------
!
      if( dam%ndamdim == 1 ) then
!
!  Update his1 with his1_tmp
!
         call mergeIntoMasterHis1d( his1_tmp )  ! merge his1_tmp into master
!             ====================     io
!$omp master
         call mergeHis1d( his1_tmp, his1 )      ! update his1 on master
!             ==========     i       io
!$omp end master
!$omp barrier
         call broadcastHis1d( his1 )            ! copy his1 to all threads
!             ==============   io
!
!  Update hisa1
!
         if( dam%lamzon )  &                    ! Note: dam%ndamdim=1 if dam%lamzon=t
               call setHisAmzon( his1, hisa1 )
!                   ===========   i     io
!
!  Reset his1_tmp
!
         indhis = his1%indhis(:)
         ! use his1 (not his1_tmp) to avoid accumulation of roundoff errors)
         hbeg = his1%hbeg - nint(his1%hbeg/his1%wbox) * his1%wbox

         call initHis1d( indhis, out%nbox, out%wbox, hbeg, his1_tmp )
!             =========    i        i         i       i       o

      else if( dam%ndamdim == 2 ) then
!
!  Update his2 with his2_tmp
!
         call mergeIntoMasterHis2d( his2_tmp )  ! merge all his2_tmp into master
!             ====================     io
!$omp master
         call mergeHis2d( his2_tmp, his2 )      ! update his2 on master
!             ==========     i       io
!$omp end master
!$omp barrier
         call broadcastHis2d( his2 )            ! copy his2 to all threads
!             ==============   io
!
!  Reset his2_tmp
!
         indhis = his2%indhis(:)
         ! use his2 (not his2_tmp) to avoid accumulation of roundoff errors)
         hbeg2(:) = his2%hbeg(:) - nint(his2%hbeg(:)/his2%wbox(:)) * his2%wbox(:)

         call initHis2d( indhis, out%nbox2(:), out%wbox2(:), hbeg2(:), his2_tmp )
!             =========    i          i             i           i         o

      else if( dam%ndamdim == 3 ) then
!
!  Update his3 with his3_tmp
!
         call mergeIntoMasterHis3d( his3_tmp )  ! merge all his3_tmp into master
!             ====================     io
!$omp master
         call mergeHis3d( his3_tmp, his3 )      ! update his3 on master
!             ==========     i       io
!$omp end master
!$omp barrier
         call broadcastHis3d( his3 )            ! copy his2 to all threads
!             ==============   io
!
!  Reset his3_tmp
!
         indhis = his3%indhis(:)
         ! use his3 (not his3_tmp) to avoid accumulation of roundoff errors)
         hbeg3(:) = his3%hbeg(:) - nint(his3%hbeg(:)/his3%wbox(:)) * his3%wbox(:)

         call initHis3d( indhis, out%nbox3(:), out%wbox3(:), hbeg3(:), his3_tmp )
!             =========    i          i             i           i         o

      end if

   end subroutine updateDamageHistogram
!
!=======================================================================
!
   subroutine setHisAmzon( his1, hisa1 )
!             ===========   i     io
!
!     Construct the amorphous zone histogram
!
!     To be called from within a parallel region
!
!  Arguments:
!
      type(His1d_type), intent(in) :: his1      ! 1d histogram
      type(His1d_type), intent(inout) :: hisa1  ! 1d histogram with amorphous zones
!
!  Local Variables:
!
      integer(INT4) :: ibox, ielem, ielem1
      real(REAL8) :: histpd, rand
!
!  Start of Executable Code:
!  ------------------------
!
      hisa1 = his1

      if( ldyn ) then
         ielem1 = 1
      else
         ielem1 = nelemIon+1
      end if
      do ibox = hisa1%imin, hisa1%imax
         do ielem = ielem1, nelemIn
            histpd = his1%hist(ielem,ibox)
            call rannum( rand )
!                ======   o
            if( rand < histpd/mat(irxtal)%dens ) then
               hisa1%hist(ielem,ibox) = mat(irxtal)%dens
            else
               hisa1%hist(ielem,ibox) = ZERO
            end if
         end do
      end do

   end subroutine setHisAmzon
!
!=======================================================================
!
   subroutine updateAmo( )
!             =========
!
!     Update the amorphous zones. Must be called only if dam%libiia=t.
!
!     This subroutine relies on his1 being the same on all threads.
!     The threads are "almost synced" on exit. Thread 1 may lag behind a little since it
!     calculates the current amorphous zone boundaries while the other images already exit.
!     Note that nionhis=0 requires nthread=1.
!
!  Arguments:
!
!     None
!
!  Local Variables:
!
      logical :: lamo_left, lamo_right
      character(256) :: errmsg
      integer(INT4) :: istat
      integer(INT4) :: iamo, jamo, kamo, iamo_cpdd, namo_cpdd
      integer(INT4) :: iamo_del, jamo_del, ndel, nbox, nbox_old
      real(REAL8), allocatable :: zamo_tmp(:), dzamo_tmp(:)
      real(REAL8), allocatable :: zamo_cpdd_left(:), zamo_cpdd_right(:)

!
!  Start of Executable Code:
!  ------------------------
!
!  Move existing amorphous zone boundaries due to IBIIA
!
      if( nionhis > 0 ) then
         do iamo = 1, namo
!$omp atomic update
            zamo_left(iamo) = zamo_left(iamo) + dzamo_left(iamo)
!$omp atomic update
            zamo_right(iamo) = zamo_right(iamo) + dzamo_right(iamo)
         end do
         dzamo_left(:) = ZERO
         dzamo_right(:) = ZERO
      end if
!
!  Make sure the arrays have sufficient size (size of his1%hist may have increased)
!  If size of arrays changes, this has to be done for all images and causes implicit
!  synchronization of images
!
      nbox = his1%nbox

      if( allocated(zamo_left) ) then

         nbox_old = size( zamo_left ) - 2

         if( nbox /= nbox_old ) then

!$omp master
            allocate( zamo_tmp( nbox+2 ), stat=istat, errmsg=errmsg )
            call check_allocation( 'updateAmo', 'zamo_tmp', istat, errmsg )
            zamo_tmp(1:namo) = zamo_left(1:namo)
            call move_alloc(zamo_tmp, zamo_left)

            allocate( zamo_tmp( nbox+2 ), stat=istat, errmsg=errmsg )
            call check_allocation( 'updateAmo', 'zamo_tmp', istat, errmsg )
            zamo_tmp(1:namo) = zamo_right(1:namo)
            call move_alloc(zamo_tmp, zamo_right)
!$omp end master

            if( nionhis > 0 ) then
               allocate( dzamo_tmp( nbox+2 ), stat=istat, errmsg=errmsg )
               call check_allocation( 'updateAmo', 'dzamo_tmp', istat, errmsg )
               dzamo_tmp(1:namo) = dzamo_left(1:namo)
               call move_alloc(dzamo_tmp, dzamo_left)

               allocate( dzamo_tmp( nbox+2 ), stat=istat, errmsg=errmsg )
               call check_allocation( 'updateAmo', 'dzamo_tmp', istat, errmsg )
               dzamo_tmp(1:namo) = dzamo_right(1:namo)
               call move_alloc(dzamo_tmp, dzamo_right)
            end if

         end if

      else

!$omp master
         allocate( zamo_left( nbox+2 ), stat=istat, errmsg=errmsg )
         call check_allocation( 'updateAmo', 'zamo_left', istat, errmsg )

         allocate( zamo_right( nbox+2 ), stat=istat, errmsg=errmsg )
         call check_allocation( 'updateAmo', 'zamo_right', istat, errmsg )
!$omp end master

         namo = 0

         if( nionhis > 0 ) then
            allocate( dzamo_left( nbox+2 ), stat=istat, errmsg=errmsg )
            call check_allocation( 'updateAmo', 'dzamo_left', istat, errmsg )
            dzamo_left(:) = ZERO

            allocate( dzamo_right( nbox+2 ), stat=istat, errmsg=errmsg )
            call check_allocation( 'updateAmo', 'dzamo_right', istat, errmsg )
            dzamo_right(:) = ZERO
         end if

      end if
!
!  Since the amorphous zones are currently only used to write them to the OUT file (which
!  is done from image 1), the remainder needs to be done for image 1 only
!
!$omp master
!
!  Allocate temporary arrays
!
      allocate( zamo_cpdd_left( nbox ), stat=istat, errmsg=errmsg )
      call check_allocation( 'updateAmo', 'zamo_cpdd_left', istat, errmsg )
      allocate( zamo_cpdd_right( nbox ), stat=istat, errmsg=errmsg )
      call check_allocation( 'updateAmo', 'zamo_cpdd_right', istat, errmsg )
!
!  Eliminate non-monotonous zamo_right's and zamo_left's
!  (unlikely, but can happen due to discretization errors)
!
      do iamo = 1, namo-1
         if( zamo_left(iamo) > zamo_left(iamo+1) ) then
            zamo_left( [ iamo,iamo+1 ] ) = zamo_left( [ iamo+1,iamo ] )
         end if
         if( zamo_right(iamo) > zamo_right(iamo+1) ) then
            zamo_right( [ iamo,iamo+1 ] ) = zamo_right( [ iamo+1,iamo ] )
         end if
      end do
!
!  Join overlapping amorphous zones
!
      iamo = 1
      do
         if( iamo >= namo )  exit
         if( zamo_right(iamo) > zamo_left(iamo+1) ) then    ! eliminate both
            zamo_left(iamo:namo-1) = zamo_left(iamo+1:namo)
            zamo_right(iamo:namo-1) = zamo_right(iamo+1:namo)
            namo = namo - 1
         else                                               ! check next
            iamo = iamo + 1
         end if
      end do
!
!  Recalculate amorphous zones from damamo and damage histogram
!
      call getAmoCPDD( his1, namo_cpdd, zamo_cpdd_left, zamo_cpdd_right )
!          ==========   i        o            o               o
!
!  Join the amorphous zones according to IBIIA and CPDD
!
      do iamo_cpdd = 1, namo_cpdd
                        ! determine left-boundary ibiia intervals where cpdd boundaries lie
         do iamo = 1, namo
            if( zamo_left(iamo) > zamo_cpdd_left(iamo_cpdd) )  exit
         end do
         iamo = iamo - 1                        ! iamo now points to the largest left ibiia
                                                !    boundary smaller than the left cpdd
                                                !    boundary
         if( iamo == 0 ) then                   ! lamo_left indicates left cpdd boudary is in
            lamo_left = .false.                 !    amorphous ibiia region
         else
            lamo_left = ( zamo_cpdd_left(iamo_cpdd) <= zamo_right(iamo) )
         end if

         do jamo = 1, namo
            if( zamo_left(jamo) > zamo_cpdd_right(iamo_cpdd) )  exit
         end do
         jamo = jamo - 1                        ! jamo now points to the largest left ibiia
                                                !    boundary smaller than the right cpdd
                                                !    boundary
         lamo_right = ( zamo_cpdd_right(iamo_cpdd) <= zamo_right(jamo) )
                                                ! lamo_right indicates right cpdd boundary is
                                                !    in amorphous ibiia region
         if( iamo == jamo .and. .not.lamo_left ) then
                                                ! insert new interval
            if( namo >= out%nbox+2 )  error stop 'damage: too many amorphous zones'
            do kamo = namo, iamo+1, -1
               zamo_left(kamo+1) = zamo_left(kamo)
               zamo_right(kamo+1) = zamo_right(kamo)
            end do
            zamo_left(iamo+1) = zamo_cpdd_left(iamo_cpdd)
            zamo_right(iamo+1) = zamo_cpdd_right(iamo_cpdd)
            namo = namo + 1
         else
                                                ! delete ibiia boundaries not needed any more
            if( lamo_left ) then
               iamo_del = iamo + 1
            else
               iamo_del = iamo + 2              ! iamo_del and jamo_del is the range of
            end if                              !    zamo_left boundaries to be deleted
            jamo_del = jamo
            ndel = jamo_del - iamo_del + 1

            do kamo = jamo_del+1, namo
               zamo_left(kamo-ndel) = zamo_left(kamo)
            end do

            iamo_del = iamo_del - 1             ! iamo_del and jamo_del is the range of
            jamo_del = jamo_del - 1             !    zamo_right boundaries to be deleted

            do kamo = jamo_del+1, namo
               zamo_right(kamo-ndel) = zamo_right(kamo)
            end do

            namo = namo - ndel
                                                ! move boundary points if needed
            if( .not.lamo_left )  zamo_left(iamo+1) = zamo_cpdd_left(iamo_cpdd)

            if( .not.lamo_right )  zamo_right(jamo) = zamo_cpdd_right(iamo_cpdd)
         end if
      end do

!$omp end master

   end subroutine updateAmo
!
!=======================================================================
!
   subroutine stat30( )
!             ======
!
!     Write statistics to output files
!
!     This subroutine must be called from within a parallel region
!
!  Arguments:
!
!     None
!
!  Local Variables:
!
      character(9) :: ctime
      character(80) :: title
      character(256) :: errmsg
      character(12) :: headerIn(nelemIn), headerAtom(nelemAtom), headerProj(nelemProj)
      integer(INT4) :: ielem, ielem1, jelem, ielema
      integer(INT4) :: i, iamo, namo_cpdd, istat
      integer(INT4) :: mrand, iarand, ibrand, irand
      real(REAL8) :: xmin, xmax
      real(REAL8) :: posifm
      real(REAL8), allocatable :: hmax(:)
      real(REAL8), allocatable :: hcrit(:)
      real(REAL8), allocatable :: hmin(:)
      real(REAL8), allocatable :: zamo_cpdd_left(:), zamo_cpdd_right(:)

      character(6), save :: fmode = 'write'
!
!  Start of Executable Code:
!  ------------------------
!
!  Gather statistics of other threads into master thread
!
      call mergeIntoMasterStat( )
!          ===================
!
!$omp master
!
!  Open the output file
!
      call openFile( 21, 'OUT', 'formatted', 'append' )
!          ========  i     i         i          i
!
!  Write header to output file
!
      write( 21, '(/t2,a)' ) 'SIMULATION RESULTS:'
!$omp end master
!
!  Call specialized statistic subroutines
!
      if( geom%lrobin )  call rob30( )
!                             =====
                                       ! output of damage statistics at the 
                                       ! end of this subroutine
!
!  Write deterministic defects 
!
      if( out%ldef )  call det30( )
!                          =====  
      if( dam%ldet .and. out%lhis )  call his30( his1 )
!                                         =====   io
!
!  Column headers
!
!$omp master
      if( ldyn ) then
         ielem1 = 1
      else
         ielem1 = nelemIon + 1
      end if

      do jelem = 1, nelemIon
         write( headerIn(jelem), '( a, i2, a, a, a )' ) &
              '  Ion', jelem, '(', elem(jelem)%name, ')'
      end do
      if( dam%ldam .or. dam%lrcoil ) then
         do ielem = ielem1, nelem
            write( headerIn(jelem), '( a, i2, a, a, a )' ) &
                 ' Atom', ielem, '(', elem(ielem)%name, ')'
            jelem = jelem + 1
         end do
      end if
      if( dam%lrcoil ) then
         do ielem = ielem1, nelem
            write( headerIn(jelem), '( a, i2, a, a, a )' ) &
                 '  Vac', ielem, '(', elem(ielem)%name, ')'
            jelem = jelem + 1
         end do
      end if

      do jelem = 1, nelemIon
         write( headerProj(jelem), '( a, i2, a, a, a )' ) &
              '  Ion', jelem, '(', elem(jelem)%name, ')'
      end do
      if( dam%lrcoil ) then
         do ielem = ielem1, nelem
            write( headerProj(jelem), '( a, i2, a, a, a )' ) &
                 ' Atom', ielem, '(', elem(ielem)%name, ')'
            jelem = jelem + 1
         end do
      end if

      do jelem = 1, nelemIon
         write( headerAtom(jelem), '( a, i2, a, a, a )' ) &
              '  Ion', jelem, '(', elem(jelem)%name, ')'
      end do
      do ielem = ielem1, nelem
         write( headerAtom(jelem), '( a, i2, a, a, a )' ) &
              ' Atom', ielem, '(', elem(ielem)%name, ')'
         jelem = jelem + 1
      end do

      write( 21, '( /t2, a )' ) 'Energy statistics:'
!
!  Write energy deposition statistics
!
      write( 21, '( /t5, a )' ) 'Mean energy per ion (keV) deposited into ... by projectile:'

      write( 21, '( t5, a )', advance='no' ) '+---------------------+'
      write( 21, '( 100a )' ) ('------------+', i=1,nelemAtom)

      write( 21, '( t5, a )', advance='no' ) '|       ...\Proj      |'
      write( 21, '( 100( a, a ) )' )  (headerAtom(i), '|', i=1,nelemAtom)

      write( 21, '( t5, a )', advance='no' ) '+---------------------+'
      write( 21, '( 100( a ) )' ) ('------------+', i=1,nelemAtom)
!
!  sumdee
!
      write( 21, '( t5, a )', advance='no' ) '| Electrons           |'
      write( 21, '( 100( f11.3, a ) )' ) ( sumdee(i)/(1000*ion%dose), ' |', i=1,nelemAtom )
      write( 21, '( t5, a )', advance='no' ) '| Phonons+Potential   |'
      write( 21, '( 100( f11.3, a ) )' ) ( sumedam(i)/(1000*ion%dose), ' |', i=1,nelemAtom )

      write( 21, '( t5, a )', advance='no' ) '| Backscattered atoms |'
      write( 21, '( 100( f11.3, a ) )', advance='no' )  &
              ( sumeb(i)/(1000*ion%dose), ' |', i=1,nelemProj )
      write( 21, '( 100( f11.3, a ) )' ) ( ZERO, ' |', i=nelemProj+1,nelemAtom )

      write( 21, '( t5, a )', advance='no' ) '| Transmitted atoms   |'
      write( 21, '( 100( f11.3, a ) )', advance='no' )  &
              ( sumet(i)/(1000*ion%dose), ' |', i=1,nelemProj )
      write( 21, '( 100( f11.3, a ) )' ) ( ZERO, ' |', i=nelemProj+1,nelemAtom )

      write( 21, '( t5, a )', advance='no' ) '+---------------------+'
      write( 21, '( 100( a ) )' ) ('------------+', i=1,nelemAtom)
!
!  Write energy transfer statistics
!
      write( 21, '( /t5, a )' ) 'Mean energy per ion (keV) transferred to recoil by projectile:'

      write( 21, '( t5, a )', advance='no' ) '+------------+'
      write( 21, '( 100a )' ) ('------------+', i=1,nelemProj)

      write( 21, '( t5, a )', advance='no' ) '| Recoil\Proj|'
      write( 21, '( 100( a, a ) )' )  (headerProj(i), '|', i=1,nelemProj)

      write( 21, '( t5, a )', advance='no' ) '+------------+'
      write( 21, '( 100( a ) )' ) ('------------+', i=1,nelemProj)

      do ielema = 1, nelemAtom
         write( 21, '( ( t5, a, a, a ) )', advance='no' ) '|', headerAtom(ielema), '|'
         write( 21, '( 100( f11.3, a ) )' )  &
                 ( sumden(ielema,i)/( 1000*ion%dose ), ' |', i=1,nelemProj )
      end do

      write( 21, '( t5, a )', advance='no' ) '+------------+'
      write( 21, '( 100( a ) )' ) ('------------+', i=1,nelemProj)
!
!  Report yields
!
      write( 21, '( /t2, a )' ) 'Yield moments:'

      write( title, '(a)' ) 'Stopped atoms and vacancies per ion:'
      call writeMom1d( momyldin, title, headerIn, 4 )
!          ==========     i        i       i

      write( title, '(a)' ) 'Nuclear energy deposition (eV) into recoils per ion:'
      call writeMom1d( momylde, title, headerAtom, 4 )
!          ==========     i        i       i

      write( title, '(a)' ) 'Backscattered atoms per ion:'
      call writeMom1d( momyldb, title, headerProj, 4 )
!          ==========    i       i       i

      write( title, '(a)' ) 'Transmitted atoms per ion:'
      call writeMom1d( momyldt, title, headerProj, 4 )
!          ==========    i       i       i
!
      if( all(nlost(:) == 0) ) then
         write( 21, '( /t5, a )' ) 'No atoms lost'
      else
         write( 21, '( )' )
         do ielem = 1, nelemProj
            if( nlost(ielem) == 0 ) cycle
            write( 21, '( t2, i10, a, a, a )' )  &
                 nlost(ielem), ' atoms', headerProj(ielem)(6:), ' lost'
         end do
      end if
! 
!  Write spatial moments of stopped atoms and vacancies
!
      if( out%lmom ) then

         write( 21, '( /t2, a )' )  'Spatial moments of stopped atoms and vacancies:'

         write( title, '(a)' ) 'Vertical moments (A):'
         call writeMom1d( momz, title, headerIn, 4 )
!             ==========   i      i       i      i
         if( out%nlatmom > 0 ) then
            write( title, '(a)' ) 'Lateral x-moments (A):'
            call writeMom1d( momx, title, headerIn, out%nlatmom )
!                ==========   i      i       i          i
            write( title, '(a)' ) 'Lateral y-moments (A):'
            call writeMom1d( momy, title, headerIn, out%nlatmom )
!                ==========   i      i       i          i
         end if

         if( out%nmixmom > 0 ) then
            write( title, '(a)' ) 'Mixed xy-moments:'
            call writeMom2d( momxy, title, headerIn, out%nmixmom )
!                ==========    i      i       i           i
            write( title, '(a)' ) 'Mixed zx-moments:'
            call writeMom2d( momzx, title, headerIn, out%nmixmom )
!                ==========    i      i       i           i
            write( title, '(a)' ) 'Mixed zy-moments:'
            call writeMom2d( momzy, title, headerIn, out%nmixmom )
!                ==========    i      i       i           i
         end if

      end if
! 
!  Write spatial moments of electronic energy deposition
!
      if( out%lmomee ) then

         write( 21, '( /t2, a )' )  'Spatial moments of electronic energy deposition:'

         write( title, '(a)' ) 'Vertical moments (A):'
         call writeMom1d( momzee, title, headerAtom, 4 )
!             ==========    i       i        i      i
         if( out%nlatmom > 0 ) then
            write( title, '(a)' ) 'Lateral x-moments (A):'
            call writeMom1d( momxee, title, headerAtom, out%nlatmom )
!                ==========    i       i        i           i
            write( title, '(a)' ) 'Lateral y-moments (A):'
            call writeMom1d( momyee, title, headerAtom, out%nlatmom )
!                ==========    i       i        i           i
         end if

         if( out%nmixmom > 0 ) then
            write( title, '(a)' ) 'Mixed xy-moments:'
            call writeMom2d( momxyee, title, headerAtom, out%nmixmom )
!                ==========     i       i        i           i
            write( title, '(a)' ) 'Mixed zx-moments:'
            call writeMom2d( momzxee, title, headerAtom, out%nmixmom )
!                ==========     i       i        i           i
            write( title, '(a)' ) 'Mixed zy-moments:'
            call writeMom2d( momzyee, title, headerAtom, out%nmixmom )
!                ==========     i       i        i           i
         end if

      end if
!
!  Write spatial moments of nuclear energy deposition
!
      if( out%lmomne ) then

         write( 21, '( /t2, a )' )  'Spatial moments of nuclear energy deposition:'

         write( title, '(a)' ) 'Vertical moments (A):'
         call writeMom1d( momzne, title, headerAtom, 4 )
!             ==========    i       i        i      i
         if( out%nlatmom > 0 ) then
            write( title, '(a)' ) 'Lateral x-moments (A):'
            call writeMom1d( momxne, title, headerAtom, out%nlatmom )
!                ==========    i       i        i           i
            write( title, '(a)' ) 'Lateral y-moments (A):'
            call writeMom1d( momyne, title, headerAtom, out%nlatmom )
!                ==========    i       i        i           i
         end if

         if( out%nmixmom > 0 ) then
            write( title, '(a)' ) 'Mixed xy-moments:'
            call writeMom2d( momxyne, title, headerAtom, out%nmixmom )
!                ==========     i       i        i           i
            write( title, '(a)' ) 'Mixed zx-moments:'
            call writeMom2d( momzxne, title, headerAtom, out%nmixmom )
!                ==========     i       i        i           i
            write( title, '(a)' ) 'Mixed zy-moments:'
            call writeMom2d( momzyne, title, headerAtom, out%nmixmom )
!                ==========     i       i        i           i
         end if

      end if
!
!  Write energy and angular moments of backscattered atoms
!
      if( out%lmomaab .or. out%lmomab .or. out%lmomeb .or. out%lmomb )  &
            write( 21, '( /t2, a )' )  'Moments of backscattered atoms:'

      if( out%lmomaab ) then
         write( title, '(a)' ) 'Angular moments wrt analyzer direction (deg):'
         call writeMom1d( momaab, title, headerProj, 4 )
!             ==========    i       i         i
      end if

      if( out%lmomab ) then
         write( title, '(a)' ) 'Angular moments of 2D projection (deg):'
         call writeMom1d( momab, title, headerProj, 4 )
!             ==========    i      i         i
      end if

      if( out%lmomeb ) then
         write( title, '(a)' ) 'Energy moments (keV):'
         call writeMom1d( momeb, title, headerProj, 4 )
!             ==========    i      i         i
      end if

      if( out%lmomb ) then
         if( out%nlatmom > 0 ) then
            write( title, '(a)' ) 'Spatial moments of exit x-positions (A):'
            call writeMom1d( momxb, title, headerIn, out%nlatmom )
!                ==========    i      i        i         i
            write( title, '(a)' ) 'Spatial moments of exit y-positions (A):'
            call writeMom1d( momyb, title, headerIn, out%nlatmom )
!                ==========    i      i       i           i
         end if
         if( out%nmixmom > 0 ) then
            write( title, '(a)' ) 'Mixed spatial moments of exit x- and y-positions (A):'
            call writeMom2d( momxyb, title, headerIn, out%nmixmom )
!                ==========     i      i       i           i
         end if
      end if 
!
!  Write energy and angular moments of transmitted atoms
!
      if( out%lmomaat .or. out%lmomat .or. out%lmomet .or. out%lmomt )  &
            write( 21, '( /t2, a )' )  'Moments of transmitted atoms:'

      if( out%lmomaat ) then
         write( title, '(a)' ) 'Angular moments wrt analyzer direction (deg):'
         call writeMom1d( momaat, title, headerProj, 4 )
!             ==========    i       i         i
      end if

      if( out%lmomat ) then
         write( title, '(a)' ) 'Angular moments of 2D projection (deg):'
         call writeMom1d( momat, title, headerProj, 4 )
!             ==========    i      i         i
      end if

      if( out%lmomet ) then
         write( title, '(a)' ) 'Energy moments (keV):'
         call writeMom1d( momet, title, headerProj, 4 )
!             ==========    i      i         i
      end if

      if( out%lmomt ) then
         if( out%nlatmom > 0 ) then
            write( title, '(a)' ) 'Spatial moments of exit x-positions (A):'
            call writeMom1d( momxt, title, headerIn, out%nlatmom )
!                ==========   i      i        i           i
            write( title, '(a)' ) 'Spatial moments of exit y-positions (A):'
            call writeMom1d( momyt, title, headerIn, out%nlatmom )
!                ==========   i      i        i           i
         end if
         if( out%nmixmom > 0 ) then
            write( title, '(a)' ) 'Mixed spatial moments of exit x- and y-positions (A):'
            call writeMom2d( momxyt, title, headerIn, out%nmixmom )
!                ==========     i      i       i           i
         end if
      end if
!
!  Write information about amorphous zones considering IBIIA
!
      if( dam%libiia ) then
         write( 21,'(/t2,a)' ) 'Amorphous zone(s) (IBIIA + critical point defect density):'
         xmin = -99999.9
         xmax = 99999.9
         do iamo = 1, namo
            if( zamo_left(iamo) == zamo_right(iamo) )  cycle
            write( 21, '(t5,a,f8.1,a,a,f8.1,a)' ) &
               'between ', min( xmax, max( xmin, zamo_left( iamo ) ) ), ' A', &
               ' and ', min( xmax, max( xmin, zamo_right( iamo ) ) ), ' A'
         end do
      end if
!
!  Write information about amorphous zones considering damamo
!
      if( dam%ldam .and. his1%imin <= his1%imax .and. irxtal > 0 .and.  &
          ndim == 1 ) then
         if( dam%damamo < mat(irxtal)%dens ) then

            allocate( zamo_cpdd_left( his1%nbox ), stat=istat, errmsg=errmsg )
            call check_allocation( 'stat30', 'zamo_cpdd_left', istat, errmsg )
            allocate( zamo_cpdd_right( his1%nbox ), stat=istat, errmsg=errmsg )
            call check_allocation( 'stat30', 'zamo_cpdd_right', istat, errmsg )

            call getAmoCPDD( his1, namo_cpdd, zamo_cpdd_left, zamo_cpdd_right )
!                ==========   i        o          o               o

            if( namo_cpdd == 0 ) then
               write( 21, '(/t2,a)' )  &
                     'Amorphous zone(s): none'
            else
               write( 21, '(/t2,a)' )  &
                     'Amorphous zone(s) (exceeding critical point defect density):'
               call openFile( 22, 'AMO', 'formatted', fmode )
               write( 22, * ) 2, 4*namo_cpdd
               write( 22, * ) 0, 0
               xmin = -99999.9
               xmax = 99999.9
            end if

            do iamo = 1, namo_cpdd
               write( 21, '(t5,a,f8.1,a,f8.1,a)' ) &
                  'between ', min( xmax, max( xmin, zamo_cpdd_left(iamo) ) ), &
                  ' A and ', min( xmax, max( xmin, zamo_cpdd_right(iamo) ) ), ' A'
               write( 22, * ) zamo_cpdd_left(iamo), ONE
               write( 22, * ) zamo_cpdd_left(iamo), BIG
               write( 22, * ) zamo_cpdd_right(iamo), BIG
               write( 22, * ) zamo_cpdd_right(iamo), ONE
            end do

            if( namo_cpdd > 0 )  close( 22 )

         end if
      end if

      allocate( hmin( nelemIn ) )
      allocate( hmax( nelemIn ) )
      allocate( hcrit( nelemIn ) )
!
!  Limits of concentrations
!     
      hmin(:) = HISMIN*ion%ddose
      hmax(1:nelemIon) = BIG
      hcrit(1:nelemIon) = BIG
      if( irxtal == 0 ) then
         hmax(nelemIon+1:) = BIG
         hcrit(nelemIon+1:) = BIG
      else
         hmax(nelemIon+1:) = mat(irxtal)%dens
         hcrit(nelemIon+1:) = dam%damamo
      end if
! 
!  Write 1d histogram
!
      if( out%lhis ) then
                                                ! delete concentrations in the air
         if( .not.usecell .and. ndim == 1 .and. mat(1)%lsput )  &
               call deleteInAir( his1 )
!                   ===========   io

         call writeHis1d( 'HIS', fmode, 1e24_REAL8, his1, hmin(:), hmax(:), hcrit(:) )
!             ==========    i      i         i       i       i        i       i
     end if
!
!  Write displacement histogram
!
      if( out%lhisd ) then
         call writeHis1d( 'HISD', fmode, hisd1, power_=out%dpower )
!             ==========    i       i      i            i
      end if
!
!  Write pair distribution function
!
      if( out%lhisiv ) then
         call writeHis1d( 'HISIV', fmode, 1/rvaciv(:), hisiv1, power_=out%dpower )
!             ==========     i       i         i         i            i
      end if
!
!  Write displacement moment histogram
!
      if( out%lhism ) then
         call writeHis1d( 'HISM', fmode, hism1 )
!             ==========    i       i      i
      end if
!
!  Write deposited momentum histogram
!
      if( out%lhisp ) then
         call writeHis1d( 'HISP', fmode, hisp1 )
!             ==========    i       i      i
      end if
!
!  Write deposited electronic energy
!
      if( out%lhisee ) then
         call writeHis1d( 'HISEE', fmode, hisee1, HISMIN*ion%ddose )
!             ==========     i       i      i             i
      end if
!
!  Write deposited nuclear energy
!
      if( out%lhisne ) then
         call writeHis1d( 'HISNE', fmode, hisne1, HISMIN*ion%ddose )
!             ==========     i       i      i             i
      end if
!
!  Write 2d histogram
! 
      if( out%lhis2 ) then
                                                ! delete concentrations in the air
         if( .not.usecell .and. ndim == 1 .and. mat(1)%lsput )  &
               call deleteInAir( his2 )
!                   ===========   io

         call writeHis2d( 'HIS2', fmode, his2, hmin, hmax, hcrit, factor_=1e24_REAL8 )
!             ==========    i       i     i     i     i      i            i

         close( 22 )
!
      end if
!
!  Write 2d histogram of deposited electronic energy
!
      if( out%lhisee2 ) then

         call writeHis2d( 'HISEE2', fmode, hisee2 )
!             ==========     i        i      i
!
         close( 22 )

      end if
!
!  Write 2d histogram of deposited nuclear energy
! 
      if( out%lhisne2 ) then

         call writeHis2d( 'HISNE2', fmode, hisne2 )
!             ==========     i        i      i

         close( 22 )

      end if
!
!  Write 3d histogram
! 
      if( out%lhis3 ) then
                                                ! delete concentrations in the air
         if( .not.usecell .and. ndim == 1 .and. mat(1)%lsput )  &
               call deleteInAir( his3 )
!                   ===========   io

         call writeHis3d( 'HIS3', fmode, his3, hmin, hmax, hcrit, factor_=1e24_REAL8 )
!             ==========    i       i     i     i     i      i            i

      end if
!
!  Write 3d histogram of deposited electronic energy
!
      if( out%lhisee3 ) then

         call writeHis3d( 'HISEE3', fmode, hisee3 )
!             ==========     i        i      i

      end if
!
!  Write 3d histogram of deposited nuclear energy
!
      if( out%lhisne3 ) then

         call writeHis3d( 'HISNE3', fmode, hisne3 )
!             ==========     i        i      i

      end if
!
!  Write backscattered atoms histograms
! 
      if( out%lhisaab )  &
            call writeHis1d( 'HISAAB', fmode, hisaa1b, HISMIN*ion%ddose/ion%dose )
!                ==========    i        i        i               i
      if( out%lhisab )  &
            call writeHis1d( 'HISAB', fmode, hisa1b, HISMIN*ion%ddose/ion%dose )
!                ==========    i        i      i               i
      if( out%lhisa2b )  &
            call writeHis2d( 'HISA2B', fmode, hisa2b, HISMIN*ion%ddose/ion%dose )
!                ==========     i        i       i               i
      if( out%lhiseb )  &
            call writeHis1d( 'HISEB', fmode, hise1b, HISMIN*ion%ddose/ion%dose,  &
!                ==========    i       i       i               i
                             power_=out%epower )
!                                    i
      if( out%lhisb )  &
            call writeHis1d( 'HISB', fmode, his1b, HISMIN*ion%ddose/ion%dose )
!                ==========    i       i      i               i
      if( out%lhis2b )  &
            call writeHis2d( 'HIS2B', fmode, his2b, HISMIN*ion%ddose/ion%dose )
!                ==========     i       i      i               i
!
!  Write transmitted atoms histograms
! 
      if( out%lhisaat )  &
            call writeHis1d( 'HISAAT', fmode, hisaa1t, HISMIN*ion%ddose/ion%dose )
!                ==========    i        i        i                i
      if( out%lhisat )  &
            call writeHis1d( 'HISAT', fmode, hisa1t, HISMIN*ion%ddose/ion%dose )
!                ==========    i        i      i                i
      if( out%lhisa2t )  &
            call writeHis2d( 'HISA2T', fmode, hisa2t, HISMIN*ion%ddose/ion%dose )
!                ==========     i        i       i               i
      if( out%lhiset )  &
            call writeHis1d( 'HISET', fmode, hise1t, HISMIN*ion%ddose/ion%dose,  &
!                ==========    i       i       i               i
                             power_=out%epower )
!                                    i
      if( out%lhist )  &
            call writeHis1d( 'HIST', fmode, his1t, HISMIN*ion%ddose/ion%dose )
!                ==========    i       i      i               i
      if( out%lhis2t )  &
            call writeHis2d( 'HIS2T', fmode, his2t, HISMIN*ion%ddose/ion%dose )
!                ==========     i       i      i               i
!
!  Report final position of the surface
!
      if( ndim == 1 .and. .not.geom%lcell ) then
         if( mat(1)%lsput ) then
            posifm = 99999.9
            write( 21, '(/t2,a,f7.1,a)' ) &
               'Final surface position:',  min( geom%posif(0), posifm ), ' A'
         end if
      end if
!
!  Seed of random number generator and execution time
!
      if( rng /= 'Knuth' ) then
         call ranget( mrand, iarand, ibrand, irand )
!                       o      o       o       o
         write( 21, '(/t2,a, a,i5, a,i5, a,i6, a,i4)' ) &
            'Final state of random number generator: ', &
            'm=',mrand, ',ia=',iarand, ',ib=',ibrand, ',irand=',irand
      end if

      write( 21, '(a)' ) ' '
      call getCPUtime( ctime )
!          ==========    o
      write( 21, '(t2,a,a,a)' ) 'Execution time: ', ctime, ' (hh:mm:ss)'
      call getElapsedTime( ctime )
!          ==============    o
      write( 21, '(t2,a,a,a)' ) 'Elapsed time  : ', ctime, ' (hh:mm:ss)'
!
!  Close output file
!
      close( 21 )
!
!  Set file mode to 'append' (for those using fmode)
!
      fmode = 'append'

!$omp end master

   end subroutine stat30
!
!=======================================================================
!
   subroutine mergeIntoMasterStat( )
!             ===================
!
!     Merge statistics from other streads into master thread
!
!     This subroutine must be called from within a parallel region
!
!  Arguments:
!
!     None
!
!  Local Variables:
!
      integer(INT4) :: ielem
!
!  Start of Executable Code:
!  ------------------------
!
!  Merge particle counts
!
      call mergeIntoMasterReal1d( rin(:) )
!          =====================    io
      call mergeIntoMasterReal1d( redam(:) )
!          =====================     io
      call mergeIntoMasterInt1d( nback(:) )
!          ====================     io
      call mergeIntoMasterInt1d( ntrans(:) )
!          ====================     io
      call mergeIntoMasterInt1d( nlost(:) )
!          ====================     io
!
!  Merge energy loss sums
!
      call mergeIntoMasterReal1d( sumdee(:) )
!          =====================     io
      call mergeIntoMasterReal2d( sumden(:,:) )
!          =====================      io
      call mergeIntoMasterReal1d( sumedam(:) )
!          =====================      io
      call mergeIntoMasterReal1d( sumeb(:) )
!          =====================     io
      call mergeIntoMasterReal1d( sumet(:) )
!          =====================     io
!
!  Merge the yield moments
!
      do ielem = 1, nelemIn
         call mergeIntoMasterMom1d( momyldin(ielem) )
!             ====================        io
      end do

      do ielem = 1, nelemAtom
         call mergeIntoMasterMom1d( momylde(ielem) )
!             ====================       io
      end do

      do ielem = 1, nelemProj
         call mergeIntoMasterMom1d( momyldb(ielem) )
!             ====================       io
      end do

      do ielem = 1, nelemProj
         call mergeIntoMasterMom1d( momyldt(ielem) )
!             ====================       io
      end do
!
!  Merge the spatial moments of stopped atoms
!
      if( out%lmom ) then

         do ielem = 1, nelemIn
            call mergeIntoMasterMom1d( momz(ielem) )
!                ====================     io
         end do

         if( out%nlatmom > 0 ) then
            do ielem = 1, nelemIn
               call mergeIntoMasterMom1d( momx(ielem) )
!                   ====================     io
            end do

            do ielem = 1, nelemIn
               call mergeIntoMasterMom1d( momy(ielem) )
!                   ====================     io
            end do
         end if

         if( out%nmixmom > 0 ) then
            do ielem = 1, nelemIn
               call mergeIntoMasterMom2d( momxy(ielem) )
!                   ====================      io
               call mergeIntoMasterMom2d( momzx(ielem) )
!                   ====================      io
               call mergeIntoMasterMom2d( momzy(ielem) )
!                   ====================      io
            end do
         end if
      end if
!
!  Merge the spatial moments of electronic energy deposition
!
      if( out%lmomee ) then

         do ielem = 1, nelemAtom
            call mergeIntoMasterMom1d( momzee(ielem) )
!                ====================      io
         end do

         if( out%nlatmom > 0 ) then
            do ielem = 1, nelemAtom
               call mergeIntoMasterMom1d( momxee(ielem) )
!                   ====================      io
            end do

            do ielem = 1, nelemAtom
               call mergeIntoMasterMom1d( momyee(ielem) )
!                   ====================      io
            end do
         end if

         if( out%nmixmom > 0 ) then
            do ielem = 1, nelemAtom
               call mergeIntoMasterMom2d( momxyee(ielem) )
!                   ====================       io
               call mergeIntoMasterMom2d( momzxee(ielem) )
!                   ====================       io
               call mergeIntoMasterMom2d( momzyee(ielem) )
!                   ====================       io
            end do
         end if
      end if
!
!  Merge the spatial moments of nuclear energy deposition
!
      if( out%lmomne ) then

         do ielem = 1, nelemAtom
            call mergeIntoMasterMom1d( momzne(ielem) )
!                ====================      io
         end do

         if( out%nlatmom > 0 ) then
            do ielem = 1, nelemAtom
               call mergeIntoMasterMom1d( momxne(ielem) )
!                   ====================      io
            end do

            do ielem = 1, nelemAtom
               call mergeIntoMasterMom1d( momyne(ielem) )
!                   ====================      io
            end do
         end if

         if( out%nmixmom > 0 ) then
            do ielem = 1, nelemAtom
               call mergeIntoMasterMom2d( momxyne(ielem) )
!                   ====================       io
               call mergeIntoMasterMom2d( momzxne(ielem) )
!                   ====================       io
               call mergeIntoMasterMom2d( momzyne(ielem) )
!                   ====================       io
            end do
         end if
      end if
!
!  Merge the backscattered and transmitted ions moments
!
      if( out%lmomaab ) then
         do ielem = 1, nelemProj
            call mergeIntoMasterMom1d( momaab(ielem) )
!                ====================      io
         end do
      end if

      if( out%lmomab ) then
         do ielem = 1, nelemProj
            call mergeIntoMasterMom1d( momab(ielem) )
!                ====================      io
         end do
      end if

      if( out%lmomeb ) then
         do ielem = 1, nelemProj
            call mergeIntoMasterMom1d( momeb(ielem) )
!                ====================      io
         end do
      end if

      if( out%lmomb ) then
         do ielem = 1, nelemIn
            if( out%nlatmom > 0 ) then
               call mergeIntoMasterMom1d( momxb(ielem) )
!                   ====================      io
               call mergeIntoMasterMom1d( momyb(ielem) )
!                   ====================      io
            end if
            if( out%nmixmom > 0 )  call mergeIntoMasterMom2d( momxyb(ielem) )
!                                       ====================      io
         end do
      end if

      if( out%lmomaat ) then
         do ielem = 1, nelemProj
            call mergeIntoMasterMom1d( momaat(ielem) )
!                ====================      io
         end do
      end if

      if( out%lmomat ) then
         do ielem = 1, nelemProj
            call mergeIntoMasterMom1d( momat(ielem) )
!                ====================      io
         end do
      end if

      if( out%lmomet ) then
         do ielem = 1, nelemProj
            call mergeIntoMasterMom1d( momet(ielem) )
!                ====================      io
         end do
      end if

      if( out%lmomt ) then
         do ielem = 1, nelemIn
            if( out%nlatmom > 0 ) then
               call mergeIntoMasterMom1d( momxt(ielem) )
!                   ====================      io
               call mergeIntoMasterMom1d( momyt(ielem) )
!                   ====================      io
            end if
            if( out%nmixmom > 0 )  call mergeIntoMasterMom2d( momxyt(ielem) )
!                                       ====================      io
         end do
      end if
!
!  Merge the 1d histogram
!
      if( out%lhis .and. .not.( nionhis > 0 .and. dam%ndamdim == 1 ) ) then
         call mergeIntoMasterHis1d( his1 )
!             ====================   io
      end if
!
!  Merge other 1D histograms
!
      if( out%lhisd )  &
            call mergeIntoMasterHis1d( hisd1 )
!                ====================   io
      if( out%lhisiv ) then
         call mergeIntoMasterHis1d( hisiv1 )
!             ====================    io
         call mergeIntoMasterReal1d( rvaciv(:) )
!             =====================     io
      end if
      if( out%lhism )  &
            call mergeIntoMasterHis1d( hism1 )
!                ====================   io
      if( out%lhisee )  &
            call mergeIntoMasterHis1d( hisee1 )
!                ====================    io
      if( out%lhisne )  &
            call mergeIntoMasterHis1d( hisne1 )
!                ====================    io
      if( out%lhisp )  &
            call mergeIntoMasterHis1d( hisp1 )
!                ====================   io
!
!  Merge the 2d histogram
!
      if( out%lhis2 .and. .not.( nionhis > 0 .and. dam%ndamdim == 2 ) ) then
         call mergeIntoMasterHis2d( his2 )
!             ====================   io
      end if
!
!  Merge other 2D histograms
!
      if( out%lhisee2 )  &
            call mergeIntoMasterHis2d( hisee2 )
!                ====================    io
      if( out%lhisne2 )  &
            call mergeIntoMasterHis2d( hisne2 )
!                ====================    io
!
!  Merge the 3d histogram
!
      if( out%lhis3 .and. .not.( nionhis > 0 .and. dam%ndamdim == 3 ) ) then
         call mergeIntoMasterHis3d( his3 )
!             ====================   io
      end if
!
!  Merge other 3d histograms
!
      if( out%lhisee3 )  &
            call mergeIntoMasterHis3d( hisee3 )
!                ====================    io
      if( out%lhisne3 )  &
            call mergeIntoMasterHis3d( hisne3 )
!                ====================    io
!
!  Merge the backscattered ions statistics
!
      if( out%lhisaab )  &
            call mergeIntoMasterHis1d( hisaa1b )
!                ====================    io
      if( out%lhisab )  &
            call mergeIntoMasterHis1d( hisa1b )
!                ====================    io
      if( out%lhisa2b )  &
            call mergeIntoMasterHis2d( hisa2b )
!                ====================    io
      if( out%lhiseb )  &
            call mergeIntoMasterHis1d( hise1b )
!                ====================    io
      if( out%lhisb )  &
            call mergeIntoMasterHis1d( his1b )
!                ====================   io
      if( out%lhis2b )  &
            call mergeIntoMasterHis2d( his2b )
!                ====================   io
!
!  Merge the transmitted ions statistics
!
      if( out%lhisaat )  &
            call mergeIntoMasterHis1d( hisaa1t )
!                ====================    io
      if( out%lhisat )  &
            call mergeIntoMasterHis1d( hisa1t )
!                ====================    io
      if( out%lhisa2t )  &
            call mergeIntoMasterHis2d( hisa2t )
!                ====================    io
      if( out%lhiset )  &
            call mergeIntoMasterHis1d( hise1t )
!                ====================    io
      if( out%lhist )  &
            call mergeIntoMasterHis1d( his1t )
!                ====================   io
      if( out%lhis2t )  &
            call mergeIntoMasterHis2d( his2t )
!                ====================   io

   end subroutine mergeIntoMasterStat
!
!=======================================================================
!
   subroutine mergeIntoMasterInt1d( i )
!             ====================  io
!
!     Merge integer array of other threads into array of master thread
!
!     To be called from within a parallel region
!
!  Arguments:
!
      integer(INT4), intent(inout) :: i(:)   ! array to be merged
!
!  Local Variables:
!
      integer(INT4), allocatable, save :: i_tmp(:) ! shared temporary array
!
!  Start of Executable Code:
!  ------------------------
!
!$omp master
      allocate( i_tmp, source=i )
      i_tmp(:) = 0
!$omp end master

!$omp barrier
!$omp critical( critical_stat_mergeIntoMasterInt1d )
      i_tmp(:) = i_tmp(:) + i(:)
!$omp end critical( critical_stat_mergeIntoMasterInt1d )

!$omp barrier
!$omp master
      i(:) = i_tmp(:)
      deallocate( i_tmp )
!$omp end master

   end subroutine mergeIntoMasterInt1d
!
!=======================================================================
!
   subroutine mergeIntoMasterReal1d( r )
!             =====================  io
!
!     Merge real array of other threads into array of master thread
!
!     This subroutine must be called from within a parallel region
!
!  Arguments:
!
      real(REAL8), intent(inout) :: r(:)   ! array to be merged
!
!  Local Variables:
!
      real(REAL8), allocatable, save :: r_tmp(:) ! shared temporary array
!
!  Start of Executable Code:
!  ------------------------
!
!$omp master
      allocate( r_tmp, source=r )
      r_tmp(:) = ZERO
!$omp end master

!$omp barrier
!$omp critical( critical_stat_mergeIntoMasterReal1d )
      r_tmp(:) = r_tmp(:) + r(:)
!$omp end critical( critical_stat_mergeIntoMasterReal1d )

!$omp barrier
!$omp master
      r(:) = r_tmp(:)
      deallocate( r_tmp )
!$omp end master

   end subroutine mergeIntoMasterReal1d
!
!=======================================================================
!
   subroutine mergeIntoMasterReal2d( r )
!             =====================  io
!
!     Merge real array of other threads into array of master thread
!
!     This subroutine must be called from within a parallel region
!
!  Arguments:
!
      real(REAL8), intent(inout) :: r(:,:)      ! array to be merged
!
!  Local Variables:
!
      real(REAL8), allocatable, save :: r_tmp(:,:) ! shared temporary array
!
!  Start of Executable Code:
!  ------------------------
!
!$omp master
      allocate( r_tmp, source=r )
      r_tmp(:,:) = ZERO
!$omp end master

!$omp barrier
!$omp critical( critical_stat_mergeIntoMasterReal2d )
      r_tmp(:,:) = r_tmp(:,:) + r(:,:)
!$omp end critical( critical_stat_mergeIntoMasterReal2d )

!$omp barrier
!$omp master
      r(:,:) = r_tmp(:,:)
      deallocate( r_tmp )
!$omp end master

   end subroutine mergeIntoMasterReal2d
!
!=======================================================================
!
   subroutine stat40( iion )
!             ======   i
!
!     Write statistics to the backup file
!
!     This subroutine should be called from the master threads only after the statistics
!     obtained by the other threads have been collected an distributed.
!
!  Arguments:
!
      integer(INT4), intent(in) :: iion         ! index of ion
!
!  Local Variables:
!
      logical, save :: firstbk = .true.
      character(9) :: ctime
      integer(INT4) :: ir
      integer(INT4) :: mrand, iarand, ibrand, irand
      real(REAL8) :: time
!
!  Start of Executable Code:
!  ------------------------
!
!  Open backup file
!
      call openFile( 22, 'BCK', 'unformatted', 'write' )
!          ========  i     i          i           i
      rewind 22
!
!  Write number of simulated ions and dose
!
      if( iion == nion ) then
         write( 22 )  nion, ion%dose
      else
         write( 22 )  iion, ion%dose1 + ion%ddose*( iion-nion1 )
      end if
!
!  Write seed of random number generator
!
      call ranget( mrand, iarand, ibrand, irand )
!          ======    o      o       o       o
      write( 22 ) mrand, iarand, ibrand, irand
!
!  Write geometry
!
      if( .not.geom%lcell .and. .not.geom%usesdf ) then
         if( ndim == 1 ) then
            write( 22 ) &
               geom%posif( 0:nr )
         else
            write( 22 ) &
               geom%region(0:nr)%np
            do ir = 0, nr
               write( 22 ) &
                  geom%region(ir)%px(:), geom%region(ir)%pz(:)
            end do
         end if
      end if
!
!  Write statistical data
!
      write( 22 ) nlost(:)
      write( 22 ) sumdee(:), sumden(:,:), sumedam(:), sumeb(:), sumet(:)

      write( 22 ) momyldin(:), momylde(:), momyldb(:), momyldt(:)

      if( out%lmom ) then
         if( out%nlatmom > 0 ) then
            if( out%nmixmom > 0 ) then
               write( 22 ) momz(:), momx(:), momy(:), momxy(:), momzx(:), momzy(:)
            else
               write( 22 ) momz(:), momx(:), momy(:)
            end if
         else
            if( out%nmixmom > 0 ) then
               write( 22 ) momz(:), momxy(:), momzx(:), momzy(:)
            else
               write( 22 ) momz(:)
            end if
         end if
      end if

      if( out%lmomee ) then
         if( out%nlatmom > 0 ) then
            if( out%nmixmom > 0 ) then
               write( 22 ) momzee(:), momxee(:), momyee(:), momxyee(:), momzxee(:), momzyee(:)
            else
               write( 22 ) momzee(:), momxee(:), momyee(:)
            end if
         else
            if( out%nmixmom > 0 ) then
               write( 22 ) momzee(:), momxyee(:), momzxee(:), momzyee(:)
            else
               write( 22 ) momzee(:)
            end if
         end if
      end if

      if( out%lmomne ) then
         if( out%nlatmom > 0 ) then
            if( out%nmixmom > 0 ) then
               write( 22 ) momzne(:), momxne(:), momyne(:), momxyne(:), momzxne(:), momzyne(:)
            else
               write( 22 ) momzne(:), momxne(:), momyne(:)
            end if
         else
            if( out%nmixmom > 0 ) then
               write( 22 ) momzne(:), momxyne(:), momzxne(:), momzyne(:)
            else
               write( 22 ) momzne(:)
            end if
         end if
      end if

      if( out%lmomaab )  write( 22 ) momaab(:)
      if( out%lmomab )  write( 22 ) momab(:)
      if( out%lmomeb )  write( 22 ) momeb(:)
      if( out%lmomb ) then
         if( out%nlatmom > 0 ) then
            if( out%nmixmom > 0 ) then
               write( 22 ) momxb(:), momyb(:), momxyb(:)
            else
               write( 22 ) momxb(:), momyb(:)
            end if
         else
            if( out%nmixmom > 0 ) then
               write( 22 ) momxyb(:)
            end if
         end if
      end if

      if( out%lmomaat )  write( 22 ) momaat(:)
      if( out%lmomat )  write( 22 ) momat(:)
      if( out%lmomet )  write( 22 ) momet(:)
      if( out%lmomt ) then
         if( out%nlatmom > 0 ) then
            if( out%nmixmom > 0 ) then
               write( 22 ) momxt(:), momyt(:), momxyt(:)
            else
               write( 22 ) momxt(:), momyt(:)
            end if
         else
            if( out%nmixmom > 0 ) then
               write( 22 ) momxyt(:)
            end if
         end if
      end if

      if( out%lhis )  call saveHis1d( his1 )
!                          =========   i
      if( out%lhisd )  call saveHis1d( hisd1 )
!                           =========    i
      if( out%lhisiv ) then
         write( 22 ) rvaciv(:)
         call saveHis1d( hisiv1 )
!             =========    i
      end if
      if( out%lhism )  call saveHis1d( hism1 )
!                           =========    i
      if( out%lhisp )  call saveHis1d( hisp1 )
!                           =========    i
      if( out%lhisee )  call saveHis1d( hisee1 )
!                            =========    i
      if( out%lhisne )  call saveHis1d( hisne1 )
!                            =========    i
      if( out%lhis2 )  call saveHis2d( his2 )
!                           =========   i
      if( out%lhisee2 )  call saveHis2d( hisee2 )
!                             =========    i
      if( out%lhisne2 )  call saveHis2d( hisne2 )
!                             =========    i
      if( out%lhis3 )  call saveHis3d( his3 )
!                           =========   i
      if( out%lhisee3 )  call saveHis3d( hisee3 )
!                             =========    i
      if( out%lhisne3 )  call saveHis3d( hisne3 )
!                             =========    i
      if( out%lhisaab )  call saveHis1d( hisaa1b )
!                             =========     i
      if( out%lhisab )  call saveHis1d( hisa1b )
!                            =========    i
      if( out%lhisa2b )  call saveHis2d( hisa2b )
!                             =========    i
      if( out%lhiseb )  call saveHis1d( hise1b )
!                            =========    i
      if( out%lhisb )  call saveHis1d( his1b )
!                           =========    i
      if( out%lhis2b )  call saveHis2d( his2b )
!                            =========    i
      if( out%lhisaat )  call saveHis1d( hisaa1t )
!                             =========     i
      if( out%lhisat )  call saveHis1d( hisa1t )
!                            =========    i
      if( out%lhisa2t )  call saveHis2d( hisa2t )
!                             =========    i
      if( out%lhiset )  call saveHis1d( hise1t )
!                            =========    i
      if( out%lhist )  call saveHis1d( his1t )
!                           =========    i
      if( out%lhis2t )  call saveHis2d( his2t )
!                            =========    i
!
!  Write simulation time
!
      call getCPUtime( time )
!          ==========   o
      write( 22 ) time
!
!  Close the backup file
!
      close( 22 )
!
!  Report writing of the backup file to the output file
!
      call openFile( 21, 'OUT', 'formatted', 'append' )
!          ========  i     i         i          i

      if( nionbk /= 0 .and. firstbk ) then
         firstbk = .false.
         write( 21, '(t2,a)' ) 'BACKUPS:'
      end if

      if( iion == nionbk )  write( 21, '(a)' ) ' '
      call getCPUtime( ctime )
!          ==========    o
      write( 21, '(t2, a,i9, a,a)', advance='no' )  &
            'iion=',iion, ', texe=', ctime
      if( rng == 'Knuth' ) then
         write( 21, '()' )
      else
         write( 21, '(a,i5,a,i5,a,i6,a,i4)' )  &
               ', m=',mrand, ', ia=',iarand, ', ib=',ibrand, ', irand=',irand
      end if

      close( 21 )

   end subroutine stat40
!
!=======================================================================
!
   subroutine stat50( )
!             ======
!
!     Read statistics from the backup file
!
!     This subroutine assumes usebk=t
!
!  Arguments:
!
!     None
!
!  Local Variables:
!
      character(9) :: ctime
      integer(INT4) :: ielem, ir
      real(REAL8) :: time
!
!  Start of Executable Code:
!  ------------------------
!
!  Open backup file
!
      call openFile( 22, bkfile, 'unformatted', 'read' )
!          ========  i     i          i          i
!
!  Read and set seed of random number generator
!
      read( 22, err=930, end=930 )              ! first line has already been read
      if( useion ) then
         read( 22, err=930, end=930 )  mrand, iarand, ibrand, irand
      else
         read( 22, err=930, end=930 )           ! ignore random number seed
      end if
!
!  Read geometry
!
      if( .not.geom%lcell .and. .not.geom%usesdf ) then
         if( ndim == 1 ) then
            if( usegeom ) then
               read( 22, err=930, end=930 ) &
                  geom%posif( 0:nr )
            else                                ! skip record
               read( 22 )
            end if
         else
            if( usegeom ) then
               allocate( geom%region( 0:nr ) )
               read( 22, err=930, end=930 )  &
                  geom%region(0:nr)%np
               do ir = 0, nr
                  allocate( geom%region(ir)%px( geom%region(ir)%np+1 ) )
                  allocate( geom%region(ir)%pz( geom%region(ir)%np+1 ) )
                  read( 22, err=930, end=930 )  &
                     geom%region(ir)%px(:), geom%region(ir)%pz(:)
                  call initPolygon( geom%region(ir) )
!                      ===========       io
               end do
            else                                ! skip records
               read( 22 )
               do ir = 0, nr
                  read( 22 )
               end do
            end if
         end if
      end if
!
!  Read statistical ion data
!
      read( 22, err=930, end=930 ) nlost(:)
      read( 22, err=930, end=930 ) sumdee(:), sumden(:,:), sumedam(:), sumeb(:), sumet(:)

      read( 22, err=930, end=930 ) momyldin(:), momylde(:), momyldb(:), momyldt(:)
      if( out%lmom ) then
         if( out%nlatmom > 0 ) then
            if( out%nmixmom > 0 ) then
               read( 22, err=930, end=930 ) momz(:), momx(:), momy(:),  &
                                            momxy(:), momzx(:), momzy(:)
            else
               read( 22, err=930, end=930 ) momz(:), momx(:), momy(:)
            end if
         else
            if( out%nmixmom > 0 ) then
               read( 22, err=930, end=930 ) momz(:), momxy(:), momzx(:), momzy(:)
            else
               read( 22, err=930, end=930 ) momz(:)
            end if
         end if
      end if
      if( out%lmomee ) then
         if( out%nlatmom > 0 ) then
            if( out%nmixmom > 0 ) then
               read( 22, err=930, end=930 ) momzee(:), momxee(:), momyee(:), momxyee(:),  &
                                            momzxee(:), momzyee(:)
            else
               read( 22, err=930, end=930 ) momzee(:), momxee(:), momyee(:)
            end if
         else
            if( out%nmixmom > 0 ) then
               read( 22, err=930, end=930 ) momzee(:), momxyee(:), momzxee(:), momzyee(:)
            else
               read( 22, err=930, end=930 ) momzee(:)
            end if
         end if
      end if
      if( out%lmomne ) then
         if( out%nlatmom > 0 ) then
            if( out%nmixmom > 0 ) then
               read( 22, err=930, end=930 ) momzne(:), momxne(:), momyne(:), momxyne(:),  &
                                            momzxne(:), momzyne(:)
            else
               read( 22, err=930, end=930 ) momzne(:), momxne(:), momyne(:)
            end if
         else
            if( out%nmixmom > 0 ) then
               read( 22, err=930, end=930 ) momzne(:), momxyne(:), momzxne(:), momzyne(:)
            else
               read( 22, err=930, end=930 ) momzne(:)
            end if
         end if
      end if

      if( out%lmomaab)  read( 22, err=930, end=930 )  momaab(:)
      if( out%lmomab)  read( 22, err=930, end=930 )  momab(:)
      if( out%lmomeb)  read( 22, err=930, end=930 )  momeb(:)
      if( out%lmomb ) then
         if( out%nlatmom > 0 ) then
            if( out%nmixmom > 0 ) then
               read( 22, err=930, end=930 )  momxb(:), momyb(:), momxyb(:)
            else
               read( 22, err=930, end=930 )  momxb(:), momyb(:)
            end if
         else
            if( out%nmixmom > 0 ) then
               read( 22, err=930, end=930 )  momxyb(:)
            end if
         end if
      end if
      if( out%lmomaat)  read( 22, err=930, end=930 )  momaat(:)
      if( out%lmomat)  read( 22, err=930, end=930 )  momat(:)
      if( out%lmomet)  read( 22, err=930, end=930 )  momet(:)
      if( out%lmomt ) then
         if( out%nlatmom > 0 ) then
            if( out%nmixmom > 0 ) then
               read( 22, err=930, end=930 )  momxt(:), momyt(:), momxyt(:)
            else
               read( 22, err=930, end=930 )  momxt(:), momyt(:)
            end if
         else
            if( out%nmixmom > 0 ) then
               read( 22, err=930, end=930 )  momxyt(:)
            end if
         end if
      end if

      if( out%lhis )  call loadHis1d( his1 )
!                          =========   io
      if( out%lhisd )  call loadHis1d( hisd1 )
!                           =========   io
      if( out%lhism )  call loadHis1d( hism1 )
!                           =========   io
      if( out%lhisp )  call loadHis1d( hisp1 )
!                           =========   io
      if( out%lhisee )  call loadHis1d( hisee1 )
!                            =========    io
      if( out%lhisne )  call loadHis1d( hisne1 )
!                            =========    io
      if( out%lhis2 )  call loadHis2d( his2 )
!                           =========   io
      if( out%lhisee2 )  call loadHis2d( hisee2 )
!                             =========    io
      if( out%lhisne2 )  call loadHis2d( hisne2 )
!                             =========    io
      if( out%lhis3 )  call loadHis3d( his3 )
!                           =========   io

      if( out%lhisee3 )  call loadHis3d( hisee3 )
!                             =========    io
      if( out%lhisne3 )  call loadHis3d( hisne3 )
!                             =========    io
      if( out%lhisiv ) then
         read( 22, err=930, end=930) rvaciv(:)
         call loadHis1d( hisiv1 )
!             =========    io
      end if

      if( out%lhisaab )  call loadHis1d( hisaa1b )
!                             =========    io
      if( out%lhisab )  call loadHis1d( hisa1b )
!                            =========    io
      if( out%lhisa2b )  call loadHis2d( hisa2b )
!                             =========    io
      if( out%lhiseb )  call loadHis1d( hise1b )
!                            =========    io
      if( out%lhis2b )  call loadHis2d( his2b )
!                            =========   io
      if( out%lhisb )  call loadHis1d( his1b )
!                           =========    io
      if( out%lhisaat )  call loadHis1d( hisaa1t )
!                             =========    io
      if( out%lhisat )  call loadHis1d( hisa1t )
!                            =========    io
      if( out%lhisa2t )  call loadHis2d( hisa2t )
!                             =========    io
      if( out%lhiset )  call loadHis1d( hise1t )
!                            =========    io
      if( out%lhist )  call loadHis1d( his1t )
!                           =========   io
      if( out%lhis2t )  call loadHis2d( his2t )
!                            =========   io

      if( .not. useion ) then
         nlost(1:nelemIon) = 0
         sumdee(1:nelemIon) = ZERO
         sumden(:,1:nelemIon) = ZERO
         sumedam(1:nelemIon) = ZERO
         sumeb(1:nelemProj) = ZERO
         sumet(1:nelemProj) = ZERO
         do ielem = 1, nelemIn
            call initMom1d( momyldin(ielem) )
!                =========       o
         end do
         do ielem = 1, nelemAtom
            call initMom1d( momylde(ielem) )
!                =========       o
         end do
         do ielem = 1, nelemProj
            call initMom1d( momyldb(ielem) )
!                =========       o
         end do
         do ielem = 1, nelemProj
            call initMom1d( momyldt(ielem) )
!                =========       o
         end do
         if( out%lmom ) then
            do ielem = 1, nelemIn
               call initMom1d( momz(ielem) )
!                   =========     o
               if( out%nlatmom > 0 ) then
                  call initMom1d( momx(ielem) )
!                      =========     o
                  call initMom1d( momy(ielem) )
!                      =========     o
               end if
               if( out%nmixmom > 0 ) then
                  call initMom2d( momxy(ielem) )
!                      =========      o
                  call initMom2d( momzx(ielem) )
!                      =========      o
                  call initMom2d( momzy(ielem) )
!                      =========      o
               end if
            end do
         endif
         if( out%lmomee ) then
            do ielem = 1, nelemAtom
               call initMom1d( momzee(ielem) )
!                   =========      o
               if( out%nlatmom > 0 ) then
                  call initMom1d( momxee(ielem) )
!                      =========      o
                  call initMom1d( momyee(ielem) )
!                      =========      o
               end if
               if( out%nmixmom > 0 ) then
                  call initMom2d( momxyee(ielem) )
!                      =========       o
                  call initMom2d( momzxee(ielem) )
!                      =========       o
                  call initMom2d( momzyee(ielem) )
!                      =========       o
               end if
            end do
         endif
         if( out%lmomne ) then
            do ielem = 1, nelemAtom
               call initMom1d( momzne(ielem) )
!                   =========      o
               if( out%nlatmom > 0 ) then
                  call initMom1d( momxne(ielem) )
!                      =========      o
                  call initMom1d( momyne(ielem) )
!                      =========      o
               end if
               if( out%nmixmom > 0 ) then
                  call initMom2d( momxyne(ielem) )
!                      =========       o
                  call initMom2d( momzxne(ielem) )
!                      =========       o
                  call initMom2d( momzyne(ielem) )
!                      =========       o
               end if
            end do
         endif
         if( out%lmomaab ) then
            do ielem = 1, nelemProj
               call initMom1d( momaab(ielem) )
!                   =========     o
            end do
         end if
         if( out%lmomab ) then
            do ielem = 1, nelemProj
               call initMom1d( momab(ielem) )
!                   =========     o
            end do
         end if
         if( out%lmomeb ) then
            do ielem = 1, nelemProj
               call initMom1d( momeb(ielem) )
!                   =========     o
            end do
         end if
         if( out%lmomb ) then
            do ielem = 1, nelemIn
               if( out%nlatmom > 0 ) then
                  call initMom1d( momxb(ielem) )
!                      =========    o
                  call initMom1d( momyb(ielem) )
!                      =========    o
               end if
               if( out%nmixmom > 0 )  call initMom2d( momxyb(ielem) )
!                                          =========    o
            end do
         end if
         if( out%lmomaat ) then
            do ielem = 1, nelemProj
               call initMom1d( momaat(ielem) )
!                   =========     o
            end do
         end if
         if( out%lmomat ) then
            do ielem = 1, nelemProj
               call initMom1d( momat(ielem) )
!                   =========     o
            end do
         end if
         if( out%lmomet ) then
            do ielem = 1, nelemProj
               call initMom1d( momet(ielem) )
!                   =========     o
            end do
         end if
         if( out%lmomt ) then
            do ielem = 1, nelemIn
               if( out%nlatmom > 0 ) then
                  call initMom1d( momxt(ielem) )
!                      =========    o
                  call initMom1d( momyt(ielem) )
!                      =========    o
               end if
               if( out%nmixmom > 0 )  call initMom2d( momxyt(ielem) )
!                                          =========    o
            end do
         end if

         if( out%lhis )  his1%hist(1:nelemIon,:) = ZERO
         if( out%lhisd )  hisd1%hist(1:nelemIon,:) = ZERO
         if( out%lhism )  hism1%hist(1:3*nelemIon,:) = ZERO
         if( out%lhisp )  hisp1%hist(1:3*nelemIon,:) = ZERO
         if( out%lhisee )  hisee1%hist(1:nelemIon,:) = ZERO
         if( out%lhisne )  hisne1%hist(1:nelemIon,:) = ZERO
         if( out%lhis2 )  his2%hist(1:nelemIon,:,:) = ZERO
         if( out%lhisee2 )  hisee2%hist(1:nelemIon,:,:) = ZERO
         if( out%lhisne2 )  hisne2%hist(1:nelemIon,:,:) = ZERO
         if( out%lhis3 )  his3%hist(1:nelemIon,:,:,:) = ZERO
         if( out%lhisee3 )  hisee3%hist(1:nelemIon,:,:,:) = ZERO
         if( out%lhisne3 )  hisne3%hist(1:nelemIon,:,:,:) = ZERO
         if( out%lhisaab )  hisaa1b%hist(1:nelemIon,:) = ZERO
         if( out%lhisab )  hisa1b%hist(1:nelemIon,:) = ZERO
         if( out%lhisa2b )  hisa2b%hist(1:nelemIon,:,:) = ZERO
         if( out%lhiseb )  hise1b%hist(1:nelemIon,:) = ZERO
         if( out%lhisb )  his1b%hist(1:nelemIon,:) = ZERO
         if( out%lhis2b )  his2b%hist(1:nelemIon,:,:) = ZERO
         if( out%lhisaat )  hisaa1t%hist(1:nelemIon,:) = ZERO
         if( out%lhisat )  hisa1t%hist(1:nelemIon,:) = ZERO
         if( out%lhisa2t )  hisa2t%hist(1:nelemIon,:,:) = ZERO
         if( out%lhiset )  hise1t%hist(1:nelemIon,:) = ZERO
         if( out%lhist )  his1t%hist(1:nelemIon,:) = ZERO
         if( out%lhis2t )  his2t%hist(1:nelemIon,:,:) = ZERO
      end if

      read( 22, err=930, end=930 ) time
      call initCPUtime( time )
!          ===========   i
!
!  Close backup file
!
      close( 22 )
! 
!  Report the use of a backup file to the output file 
!
      call openFile( 21, 'OUT', 'formatted', 'append' )
!          ========  i     i         i          i
      write( 21, '(/t2,a,a60)' ) 'Backup file used: ', bkfile
      if( useion ) then
         call getCPUtime( ctime )
!             ==========    o
         write( 21, '(t20,a,i7,a,a)' ) 'iion=', nion1, ', texe=', ctime
      else
         write( 21, '(t2,a)' ) 'Only the damage distribution is used'
      end if
      close( 21 )

      return

930   error stop 'stat50: error during reading backup file'

   end subroutine stat50
!
!=======================================================================
!
   subroutine scoreIV( pos, newvac, dist, ielem)
!             =======   i     i       i     i
!
!     Score for I-V pair distribution function
!
!  Arguments:
!
      real(REAL8), intent(in) :: pos(3)         ! vacancy position
      logical, intent(in) :: newvac             ! flag indicating a new vacancy
      real(REAL8), intent(in) :: dist           ! distance of an interstitial from the vacancy
      integer(INT4), intent(in) :: ielem        ! atom index of the vacancy (= of the
                                                !    interstitial)
!
!  Local Variables:
!
      integer(INT4) :: ielem0
      real(REAL8) :: x, y, z
      real(REAL8) :: pos_s(3)

!
!  Start of Executable Code:
!  ------------------------
!
!  Change the coordinate system
!
      pos_s(:) = pos(:)

      x = pos_s(1)
      y = pos_s(2)
      z = pos_s(3)

      if( z < out%posivmin .or. z > out%posivmax )  return

      if( ldyn ) then
         ielem0 = 0
         if( ielem > nelem )  error stop 'scoreIV: atom index error.'
      else
         ielem0 = nelemIon
         if( ielem <= nelemIon .or. ielem > nelem )  error stop 'scoreIV: atom index error.'
      end if
!
!  1d pair distribution function
!
      if( newvac )  rvaciv(ielem) = rvaciv(ielem) + ion%ddose
      call scoreHis1d( ielem-ielem0, dist, hisiv1, ion%ddose, power_=out%dpower )
!          ==========       i         i      io        i              i

   end subroutine scoreIV
!
!=======================================================================
!
   subroutine getDamage( pos, ielem, histi, histv)
!             =========   i     i      o      o
!
!     Determine the interstitial and vacancy concentration
!
!     For ldam=.false. no damage is pretended to exist
!     For lrcoil=.false. the vacancy concentration is assumed to equal the 
!     interstitial concentration
!
!  Arguments:
!
      real(REAL8), intent(in) :: pos(3)         ! recoil position
      integer(INT4), intent(in) :: ielem        ! atom index of the defect
      real(REAL8), intent(out) :: histi         ! interstitial concentration
      real(REAL8), intent(out) :: histv         ! vacancy concentration
!
!  Local Variables:
!
      integer(INT4) :: ielemi, ielemv
      real(REAL8) :: x, y, z
      real(REAL8) :: pos_s(3)

      !integer(INT4) ::  ncall = 0   ! debug
!
!  Start of Executable Code:
!  ------------------------
!
      !ncall = ncall + 1

      if( .not.ldamdyn ) then
         histi = ZERO                           ! no damage has been recorded
         histv = ZERO
         return
      end if

      if( out%mcoord == 1 ) then
         pos_s(:) = pos(:)
      else
         call toCoordSys( coord_statis, pos, pos_s)
!             ==========       i         i     o
      end if

      x = pos_s(1)
      y = pos_s(2)
      z = pos_s(3)

      if( ldyn ) then
         if( ielem > nelem ) then
            error stop 'getDamage: atom index error.'
         end if
      else
         if( ielem <= nelemIon .or. ielem > nelem ) then
            error stop 'getDamage: atom index error.'
         end if
      end if

      ielemi = ielem
      if( ldyn )  ielemi = nelemIon + ielemi    ! this to work for ldyn=t, usecell=t would
      ielemv = nelemTarg + ielem                ! require to also read the damage histogram
      if( ldyn )  ielemv = nelemIon + ielemv


      if( dam%lamzon ) then
!
!  Amorphous zone model
!
         call getHis1d( ielemi, z, hisa1, histi )
!             ========    i     i    i      o
         if( dam%lrcoil )  call getHis1d( ielemv, z, hisa1, histv )
!                               ========    i     i    i      o

      else
!
!  Point defect model
!
         if( dam%ndamdim == 1 ) then
            call getHis1d( ielemi, z, his1, histi )
!                ========    i     i   i      o
            if( dam%lrcoil )  call getHis1d( ielemv, z, his1, histv )
!                                  ========    i     i   i      o
         else if( dam%ndamdim == 2 ) then
            call getHis2d( ielemi, [ x, z ], his2, histi )
!                ========    i        i       i      o
            if( dam%lrcoil )  call getHis2d( ielemv, [ x, z ], his2, histv )
!                                  ========    i        i       i      o
         else if( dam%ndamdim == 3 ) then
            call getHis3d( ielemi, [ x, y, z ], his3, histi )
!                ========    i          i        i      o
            if( dam%lrcoil )  call getHis3d( ielemv, [ x, y, z ], his3, histv )
!                                  ========    i          i        i      o
         end if
      end if

      if( .not.dam%lrcoil )  histv = histi

      if( irxtal /= 0 ) then
         if( histi > dam%damamo ) then
            histi = mat(irxtal)%dens
            histv = mat(irxtal)%dens
         end if
      end if

      !if( modulo(ncall, 10000) == 0 )  print *, z, histi, histv
      return

   end subroutine getDamage
!
!=======================================================================
!
end module Eval_Stat_mod
