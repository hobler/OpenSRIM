module Eval_Rbs_mod
!
!  Manage statistical results about Rutherford Backscattering
!
!  Used modules:
!  ============
!
   use Const_mod,       only: REAL8, INT4,  &
                              PI, ZERO, ONE, BIG, EPSMB, EPSMA, TWO
   use Undefine_mod,    only: defined
   use Vector_mod,      only: operator(.dot.)
   use Types_mod,       only: Proj_type

   use Param_Setup_mod, only: nelemIon, nr, nelem, useion, ndim, nion, nthread, usebk, usecell
   use Param_Geom_mod,  only: geom
   use Param_Elem_mod,  only: elem
   use Param_Ion_mod,   only: ion
   use Param_Mat_mod,   only: mat
   use Param_Out_mod,   only: out

   use Stat_His1d_mod,  only: His1d_type, saveHis1d, loadHis1d, scoreHis1d,  &
                              initHis1d, mergeIntoMasterHis1d, writeHis1d
   use Se_Rbs_mod,      only: initRevPath, findRbsEnergy
!
!  Default type of access:
!  ======================
!
   implicit none
   private 
   save
!
!  Public procedures:
!  =================
!
   public :: rbs10        ! initialize RBS histograms and outgoing path 
!            =====
   public :: rbs22        ! score event into RBS histograms
!            =====
   public :: rbs30        ! write RBS statistics into rbs file
!            =====
   public :: rbs40        ! write RBS histograms into backup file
!            =====
!
!  Private procedure:
!  =================
!
   private :: rbs50        ! load from backup file
!             =====
!
!   Private Data:
!   ============
!
   integer(INT4), allocatable :: indhis( : ) ! indicator of histogram meaning
   type(His1d_type) :: hise                  ! energy histogram
   type(His1d_type) :: hisz                  ! depth histogram
!$omp threadprivate( hise, hisz )
!
   real(REAL8), allocatable :: rm12(:)       ! (nelem) m1/m2
   real(REAL8), allocatable :: cm12(:)       ! (nelem) (1+m1/m2)^2
   real(REAL8), allocatable :: z12(:)        ! (nelem) z1*z2
   real(REAL8), allocatable :: vib1(:,:)     ! (nr, nelem) 1/(2*PI*u1^2)
   real(REAL8), allocatable :: vib2(:,:)     ! (nr, nelem) 1/(2*u1^2)
   real(REAL8), allocatable :: lec(:)        ! (nelem) light energy corrections to 
                                             ! Rutherford cross section (not used)
!
contains
!
!=========================================================================================
!
   subroutine rbs10( )
!             =====
!
!     Intialize the statistics
!
!  Local variables
!
      integer(INT4) :: ielem1, ielem2, ir!, nbox
      real(REAL8) ::  fi!, wbox
!
!  Start of executable code
!  ------------------------
!
!$omp master
!
!  Check are all flags correctly set
!
!  ndim has to equal 1
!
      if( ndim /= 1 )  error stop 'rbs10: With RBS simulations ndim=1 required'
      if( usecell )  error stop 'rbs10: usecell=t not allowed with RBS simulations'
!
!  nelemIon (number of ions)  has to be equal 1. It is easy to implement and allow more than one
!  ion species but it doesn't have physical meaninig and it is for now forbidden.
!
      if( nelemIon /= 1 )  error stop 'rbs10: With RBS simulations nelemIon = 1 required!'
!
!  Ion must be lighter than targer (for all regions and ions)
!
      do ielem1 = 1, nelemIon
         do ielem2 = nelemIon+1, nelem
            if( elem(ielem1)%m > elem(ielem2)%m )  &
                  error stop 'rbs10: ion mass < target mass required!'
         end do
      end do
!
!  Set the RBS acceptance angle
!
      if( defined( out%dira(3) ) ) then
         fi = PI - acos( out%dira(3) )
      else
         fi = ( 170. / 180. ) * PI              ! 170 deg
      end if
!
!  Initialize the RBS histogram structure
!
      allocate( indhis( nelem-nelemIon ) )
      do ielem1 = nelemIon+1, nelem
         indhis(ielem1-nelemIon) = elem(ielem1)%z
      end do
!$omp end master
!$omp barrier
      call initHis1d( indhis, out%nboxe, out%wboxe, ZERO, hise )
!          =========    i         i          i       i     o
      call initHis1d( indhis, out%nbox, out%wbox, ZERO, hisz)
!          =========    i         i          i      i     o
!$omp master
!
!  Initialize structure necessary for calculation of outgoing path
!
      call initRevPath( nelemIon, fi, ion%e )
!          ===========     i       i    i  
!
!  Allocate and initialize parameters rm12, z12, vib
!
      allocate( rm12( 1:nelem ) )
      allocate( cm12( 1:nelem ) )
      allocate( z12( 1:nelem ) )
      allocate( lec( 1:nelem ) )
      allocate( vib1( 1:nr, 1:nelem ) )
      allocate( vib2( 1:nr, 1:nelem ) )

      do ielem2 = 1, nelem
         rm12(ielem2) = elem(nelemIon)%m / elem(ielem2)%m
         z12(ielem2) = elem(nelemIon)%z * elem(ielem2)%z
         cm12(ielem2) = ( ONE + rm12(ielem2) )**2
         lec(ielem2) = ( ONE - 49. * z12(ielem2)**( 4. / 3. ) ) / ( ion%e / cm12(ielem2) )
      end do

      do ir = 1, nr
         do ielem2 = 1, nelem
            vib1(ir,ielem2) = ONE / ( 2 * PI * mat(ir)%vib(ielem2)**2 )
            vib2(ir,ielem2) = ONE / ( 2 * mat(ir)%vib(ielem2)**2 )
         end do
      end do
!$omp end master
!$omp barrier

      if( usebk ) then
         call rbs50( )
!             =====
      end if

   end subroutine rbs10
!
!=========================================================================================
!
   subroutine rbs22( proj, r2, ielem2 )
!             =====    i   i     i    
!
!     Score statistic quantitions along the trajectory
!
!  Arguments:
!
      type(Proj_type), intent(in) :: proj       ! projectile 
      real(REAL8), intent(in) :: r2             ! square of distance between collision point
                                                !    and recoil position
      integer(INT4) :: ielem2                   ! element index of recoil
!
!  Local variables
!
      integer(INT4) :: ielem1, ir
      real(REAL8) :: psi, eScatt, a, b, c, sinpsi2, cospsi, rcs=ZERO, Pgaus, Prbs, eRbs
      real(REAL8), parameter :: e2 = 1.44D1  !  e^2=[eV*A]
      logical :: valRbs
!
!  Start of executable code
!  ------------------------
!
      ir = proj%id
      ielem1 = proj%ielem
!
!  Forbid ions as a targets
!
      if( proj%pos(3) > geom%posif(nr) ) return
      if( ielem2 < nelemIon+1 .or. ielem2 > nelem )  &
            error stop 'rbs22: ielem2 < nelemIon .or. ielem2 > natom!'
!
!  Detect the real scattering angle which may lead to score in the detector 
!  Take into account actual direction of moving projectile
!
      psi = PI - abs( acos( proj%dir(:).dot.out%dira(:) ) )
      cospsi = cos( psi ) 
      sinpsi2 = ONE - cospsi**2    ! avoid sqrt()

      if( sinpsi2 > ZERO ) then     ! eventual limiting of backscattering angle
!
!  Rutherford differential cross section (Handbook of Modern Ion Beam Analysis)
!
         b = sqrt( ONE - rm12(ielem2)**2 * sinpsi2 )
         a = ( ( (z12(ielem2)*e2) / (TWO*proj%e) )**2 ) * ( ONE/sinpsi2**2 ) 
         c = ( cospsi + b )**2 / b
         rcs = a * c
!
!  Check possible division by zero
!
         if( sinpsi2 < EPSMA .or. b < EPSMB ) then
            write(*,*) 'Cross section->INFINITY'
            rcs = BIG
         end if
!
!  Take into account the thermal vibration of target atom
!
         Pgaus = vib1(ir,ielem2) * exp( - r2 * vib2(ir,ielem2) )
!
!  Energy loss during large angle scattering event. Kinematic factor
!
         eScatt = ( ( b + rm12(ielem2) * cospsi )**2 / cm12(ielem2) ) * proj%e
!
!  Simulate the energy loss in outgoing particle path. Only those atoms which have enough
!  energy to reach the detector are taken into account ( valRbs = T )
!
         call findRbsEnergy( eScatt, ielem2, proj%pos(3), eRbs, valRbs )
!             =============    i        i        i         o       o
!
!  Do not score if eRbs is not in the specified energy window
!
         if( eRbs < out%emina .or. eRbs > out%emaxa ) return

         if( valRbs ) then    ! score only if reach the detector

            Prbs = Pgaus * rcs ! Probability of being detected in RBS detector

            call scoreHis1d( ielem2-nelemIon, eRbs, hise, Prbs )
!                ==========         i          i     io     i
            call scoreHis1d( ielem2-nelemIon, proj%pos(3), hisz, Prbs )
!                ==========         i              i        io    i
         end if

      end if

   end subroutine rbs22
!
!=========================================================================================
!
   subroutine rbs30( )
!             =====   
!
!     Write rbs statistics into rbs file
!
!     This subroutine should be called from a serial region
!
!  Calling arguments:
!
!     None
!
!  Local variables
! 
      real(REAL8) :: norm                       ! normalization factor
!
!  Start of executable code
!  ------------------------
!
!  Merge histograms into master thread
!
      call mergeIntoMasterHis1d( hise )
!          ====================   io
      call mergeIntoMasterHis1d( hisz )
!          ====================   io
!$omp barrier
!$omp master
!
!  hise
!
      norm = 1e16 * ion%dose / nion
      hise%hist(:,:) =  hise%hist(:,:) * norm
                                             ! total signal of all target atom species
      call writeHis1d( 'RBS', 'write', hise, ZERO, lsum_=.true. )
!          ===========   i       i      i     i         i
      if( nelem - nelemIon > 1 )  &          ! separate signals of the target atom species
            call writeHis1d( 'RBS', 'append', hise, ZERO, lsum_=.false. )
!                ==========    i       i       i     i          i
!
!  hisz
!
      hisz%hist(:,:) =  hisz%hist(:,:) * norm
                                             ! total signal of all target atom species
      call writeHis1d( 'RBS', 'append', hisz, ZERO, lsum_=.true. )
!          ===========   i       i       i     i          i
      if( nelem - nelemIon > 1 )  &          ! separate signals of the target atom species
            call writeHis1d( 'RBS', 'append', hisz, ZERO, lsum_=.false. )
!                ==========    i       i       i     i          i
!$omp end master

   end subroutine rbs30
!
!=========================================================================================
!
   subroutine rbs40( )
!             =====   
!
!     Write rbs statistics into backup file
!
!     This subroutine should be called from thread 1 only after the statistics obtained by the
!     other threads have been collected.
!
!  Start of executable code
!  ------------------------
!
      call saveHis1d( hise )
!          =========    i
      call saveHis1d( hisz )
!          =========    i

   end subroutine rbs40
!
!===========================================================================================
!
   subroutine rbs50( )
!             =====   
!
!     Read rbs statistics from backup file
!
!  Start of executable code
!  ------------------------
!
!  Set rbs statistics to zero (on all threads)
!
      hise%hist(:,:) = ZERO
      hise%hist(:,:) = ZERO
!
!  Load the backup file if requested
!
!$omp master
      if( useion ) then
         call loadHis1d( hise )
!             =========   i
         call loadHis1d( hisz )
!             =========   i
      end if
!$omp end master
!$omp barrier

   end subroutine rbs50
!
!===========================================================================================
!
end module Eval_Rbs_mod
