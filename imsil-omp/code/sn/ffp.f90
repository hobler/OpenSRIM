module Sn_FFP_mod
!
!  Calculate (deterministic) free flight path in amorphous targets as a function of
!  ion energy for given ion species and materials parameters.
!
!
!  Used Modules:
!  ============
!
   use Const_mod,       only: INT4, REAL8,  &
                              ZERO, ONE, HALF, THIRD, PI, EPSMB, BIG

   use Types_mod,       only: Proj_type
   use Undefine_mod,    only: defined

   use Param_Setup_mod, only: nelem, nr, lsput
   use Param_Geom_mod,  only: geom
   use Param_Elem_mod,  only: elem
   use Param_Ion_mod,   only: ion
   use Param_Mat_mod,   only: mat
   use Param_Sn_mod,    only: sn
   use Param_Dam_mod,   only: dam

   use Math_Bessel_mod, only: mbes1
   use Target_mod,      only: getDomainDensities, getDomainProperty
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
   public :: initFFP                ! initialize parameters used for ffp
!            =======                !    calculation
   public :: evalFFP                ! computes the free flight path ("MCFFP")
!            =======
   public :: tabPmax                ! writes tables of maximum impact parameters to 
!            =======                !    file 'PMAX'
!
!  Private Procedures:
!  ==================
!
   private :: getPmaxTable          ! calculates a table of pmax values for a table of
!             ============          !    energy values for given elements ielem1, ielem2
   private :: fr1                   ! calculates the negative of the derivate of the 
!             ===                   !    screening function  dependent term of the 
!                                   !    continuum potential
   private :: fr1inv                ! inverse function of fr1
!             ======
!
!  Private Data:
!  ============
!
   integer(INT4), parameter :: NT = 1000     ! number of tabulated e(p) pairs
   real(REAL8) :: et(NT)                     ! tabulated energy values [eV]
!
!  pmax tables for geom%lcell=t; shape=(NT,nelem)
!
   real(REAL8), allocatable, target :: pt_cell_sput_back(:,:)
                                             !    for z<0,       if lsput=t, lback=t
   real(REAL8), allocatable, target :: pt_cell_x_back(:,:)
                                             !    for z<0,       if lsput=f, lback=t
   real(REAL8), allocatable, target :: pt_cell_sput_x(:,:)
                                             !    for 0<z<dnear, if lsput=t, lback=t
                                             !    for z<dnear,   if lsput=t, lback=f
   real(REAL8), allocatable, target :: pt_cell_x_x(:,:)
                                             !    for all other cases
!
!  pmax tables for geom%lcell=f; shape=(NT,nr,nelem)
!
   real(REAL8), allocatable, target :: pt_region_sput_back(:,:,:)
                                             !    for z<0,       if lsput=t, lback=t
   real(REAL8), allocatable, target :: pt_region_x_back(:,:,:)
                                             !    for z<0,       if lsput=f, lback=t
   real(REAL8), allocatable, target :: pt_region_sput_x(:,:,:)
                                             !    for 0<z<dnear, if lsput=t, lback=t
                                             !    for z<dnear,   if lsput=t, lback=f
   real(REAL8), allocatable, target :: pt_region_x_x(:,:,:)
                                             !    for all other cases
!
   real(REAL8) :: facint                     ! factor for interpolation
                                             !     =(NT-1)/log(et(NT)/et(1))
!
!==============================================================================
!
contains
!
!==============================================================================
!
   subroutine initFFP( )
!             =======
!
!     Initialize the the table of maximum impact parameters as a function of energy
!
!     The maximum impact parameter is calculated from the minimum scattering angle and the 
!     minimum energy transfer to be considered.
!     Ref: G.Hobler, A. Simionescu, NIM-B 102, 24, 1995. 
!
!     In addition the impact parameters are forced in the range [sn%pmaxmin,sn%pmaxmax].
!     pmax from the scattering angle criterion is given an upper limit of sn%pmaxpsi.
!     Alternatively, the minimum energy for the application of the scattering angle 
!     criterion may be given as sn%eminpsi. Note that the two criteria are not equivalent,
!     since the relationship between sn%pmaxpsi and sn%eminpsi depends on ielem1 and ielem2.
!
!     In a non-dynamic simulation with fixed material composition in each region, only the
!     elements present in each region need to be considered for the maximum impact parameter
!     selection. Therefore it makes sense to define pmax tables for each region:
!     pt_region_xxx(it,ir,ielem1). In a dynamic simulation, any element can occur in any region,
!     so the worst case has to be considered anywhere, and it makes no sense to make differnt
!     tables in different regions: pt_cell_xxx(it,ielem1).
!
!     If the flag lsput is set, possibly larger maximum impact parameters are chosen
!     in a near-surface layer to guarantee that the minimum energy transfer is not smaller
!     than the smallest surface binding energy. For projectiles traveling outside the target
!     (but usually inside the surface potential barrier) the minimum scattering angle is
!     reduced for grazingly incident ions in order to describe reflection more accurately.
!
!     This subroutine must be called from on thread only.
!
!  Arguments: 
!
!     None
!
!  Local Variables:
!
      logical :: lback
      integer(INT4) :: i, ielem1, ielem2, it, ir
      integer(INT4) :: ielem(nelem), ielemr, nelemr
      real(REAL8) :: e, ef, denffp, psiffp, eminpsi
      real(REAL8), allocatable :: ptl(:,:)      ! (NT,nelem) temporary array holding pmax
                                                !            values (indices: it,ielem2) 
      real(REAL8), pointer :: pt_region(:,:)    ! (NT,nr)
      real(REAL8), pointer :: pt_cell(:)        ! (NT)
!
!  Start of Executable Code:
!  ------------------------
!
      lback = defined( sn%psiffp2 )
!
!  Allocate and initialize arrays
!
      if( lsput ) then
         ef = min( sn%ef, dam%esurfmin )
      else
         ef = sn%ef
      end if
      e = ion%e * (ONE+EPSMB)
      et(:) = [ ( ef * ( e/ef ) ** ( real(it-1)/real(NT-1) ), it=1,NT ) ]
      facint = (NT-1) / log( et(NT)/et(1) )

      if( geom%lcell ) then

         if( lsput ) then

            if( lback ) then

               allocate( pt_cell_x_x( NT, nelem ) );        pt_cell_x_x(:,:) = ZERO
               allocate( pt_cell_sput_x( NT, nelem ) );     pt_cell_sput_x(:,:) = ZERO
               allocate( pt_cell_sput_back( NT, nelem ) );  pt_cell_sput_back(:,:) = ZERO

            else

               allocate( pt_cell_x_x( NT, nelem ) );        pt_cell_x_x(:,:) = ZERO
               allocate( pt_cell_sput_x( NT, nelem ) );     pt_cell_sput_x(:,:) = ZERO

            end if

         else

            if( lback ) then

               allocate( pt_cell_x_x( NT, nelem ) );        pt_cell_x_x(:,:) = ZERO
               allocate( pt_cell_x_back( NT, nelem ) );     pt_cell_x_back(:,:) = ZERO

            else

               allocate( pt_cell_x_x( NT, nelem ) );        pt_cell_x_x(:,:) = ZERO

            end if

         end if

      else

         if( lsput ) then

            if( lback ) then

               allocate( pt_region_x_x( NT, nr, nelem ) );        pt_region_x_x(:,:,:) = ZERO
               allocate( pt_region_sput_x( NT, nr, nelem ) );     pt_region_sput_x(:,:,:) = ZERO
               allocate( pt_region_sput_back( NT, nr, nelem ) )
                                                               pt_region_sput_back(:,:,:) = ZERO

            else

               allocate( pt_region_x_x( NT, nr, nelem ) );        pt_region_x_x(:,:,:) = ZERO
               allocate( pt_region_sput_x( NT, nr, nelem ) );     pt_region_sput_x(:,:,:) = ZERO

            end if

         else

            if( lback ) then

               allocate( pt_region_x_x( NT, nr, nelem ) );        pt_region_x_x(:,:,:) = ZERO
               allocate( pt_region_x_back( NT, nr, nelem ) );     pt_region_x_back(:,:,:) = ZERO

            else

               allocate( pt_region_x_x( NT, nr, nelem ) );        pt_region_x_x(:,:,:) = ZERO

            end if

         end if

      end if

      allocate( ptl( NT, nelem ) )
!
!  Loop over elements
!  ------------------
!
      sn%pmax = ZERO

      do ielem1 = 1, nelem

         do i = 1, 8

            select case( i )

               case( 1 )

                  if( .not.allocated( pt_cell_sput_back ) )  cycle
                  pt_cell => pt_cell_sput_back(:,ielem1)
                  denffp = min( dam%esurfmin, sn%denffp )
                  psiffp = sn%psiffp2
                  eminpsi = sn%eminpsi2

               case( 2 )

                  if( .not.allocated( pt_cell_sput_x ) ) cycle
                  pt_cell => pt_cell_sput_x(:,ielem1)
                  denffp = min( dam%esurfmin, sn%denffp )
                  psiffp = sn%psiffp
                  eminpsi = sn%eminpsi

               case( 3 )

                  if( .not.allocated( pt_cell_x_back ) )  cycle
                  pt_cell => pt_cell_x_back(:,ielem1)
                  denffp = sn%denffp
                  psiffp = sn%psiffp2
                  eminpsi = sn%eminpsi2

               case( 4 )

                  if( .not.allocated( pt_cell_x_x ) )  cycle
                  pt_cell => pt_cell_x_x(:,ielem1)
                  denffp = sn%denffp
                  psiffp = sn%psiffp
                  eminpsi = sn%eminpsi

               case( 5 )

                  if( .not.allocated( pt_region_sput_back ) )  cycle
                  pt_region => pt_region_sput_back(:,:,ielem1)
                  denffp = min( dam%esurfmin, sn%denffp )
                  psiffp = sn%psiffp2
                  eminpsi = sn%eminpsi2

               case( 6 )

                  if( .not.allocated( pt_region_sput_x ) ) cycle
                  pt_region => pt_region_sput_x(:,:,ielem1)
                  denffp = min( dam%esurfmin, sn%denffp )
                  psiffp = sn%psiffp
                  eminpsi = sn%eminpsi

               case( 7 )

                  if( .not.allocated( pt_region_x_back ) )  cycle
                  pt_region => pt_region_x_back(:,:,ielem1)
                  denffp = sn%denffp
                  psiffp = sn%psiffp2
                  eminpsi = sn%eminpsi2

               case( 8 )

                  if( .not.allocated( pt_region_x_x ) )  cycle
                  pt_region => pt_region_x_x(:,:,ielem1)
                  denffp = sn%denffp
                  psiffp = sn%psiffp
                  eminpsi = sn%eminpsi

            end select
!
!  Calculate tables for ion/target pairs
!
            do ielem2 = 1, nelem

               call getPmaxTable( ptl(:,ielem2), et(:), ielem1, ielem2, denffp, psiffp,  &
!                   ============        o          i      i       i       i       i
                                  eminpsi )
!                                    i
            end do
!
!  Calculate quantity of interest
!
            select case( i )

               case( 1, 2, 3, 4 )               ! geom%lcell=t

                  do ielem2 = 1, nelem
                     pt_cell(:) = max( pt_cell(:), ptl(:,ielem2) )
                  end do

               case( 5, 6, 7, 8 )               ! geom%lcell=f

                  do ir = 1, nr
                     call getDomainDensities( ir, nelemr, ielem )
!                         ==================  i     o       o

                     do ielemr = 1, nelemr      ! use only elements in this regions
                        ielem2 = ielem(ielemr)
                        pt_region(:,ir) = max( pt_region(:,ir), ptl(:,ielem2) )
                     end do

                  end do

            end select
!
!  Initialize sn%pmax
!
            select case( i )

               case( 1, 2, 3, 4 )               ! geom%lcell=t

                  sn%pmax = max( sn%pmax, pt_cell(1) )    ! maximum pmax at position 1

               case( 5, 6, 7, 8 )               ! geom%lcell=f

                  do ir = 1, nr
                     if( mat(ir)%lxtal ) then
                        sn%pmax = max( sn%pmax, sn%pmaxc(ir) )
                     else
                        sn%pmax = max( sn%pmax, pt_region(1,ir) )
                     end if
                  end do

            end select

         end do

      end do

      deallocate( ptl )
!
! debug
!      if( .not.geom%lcell )  error stop 'initFFP debug output only for ldyn=t or usecell=t'
!      open( 25, file='PMAX', form='formatted' )
!      write( 25, * ) 1+nelem, NT
!      write( 25, * ) 0, (elem(ielem1)%z, ielem1=1,nelem)
!      do it = 1, NT
!         write( 25, * ) et(it), ptdyn(it,:)
!      end do
!      close( 25 )
!      error stop
! end debug
!      
   end subroutine initFFP
!
!==============================================================================
!
   subroutine evalFFP( proj, id, ffp, pmax )
!             =======   i    i    o    o
!
!     Computes the (deterministic) free flight path as a function of projectile energy.
!     The calculation is done by linear interpolation in the pre-calculated table.
!
!  Arguments:
!
      type(Proj_type), intent(in) :: proj       ! projectile
      integer(INT4), intent(in) :: id           ! domain index
      real(REAL8), intent(out) :: ffp           ! free flight path
      real(REAL8), intent(out) :: pmax          ! corresponding maximum impact parameter
!
!  Local Variables:
!
      integer(INT4) :: ielem1, it
      real(REAL8) :: e, distsurf
      real(REAL8) :: fac, dens
      real(REAL8), pointer :: pt(:)             ! impact parameter table
!
!  Start of Executable Code
!  ------------------------
!
      e = proj%e
      ielem1 = proj%ielem
      distsurf = proj%distsurf
!
!  Make pt point to the right section of the right pmax table array
!
      if( geom%lcell ) then

         if( lsput ) then
            if( distsurf > dam%dnear ) then
               pt => pt_cell_x_x(:,ielem1)
            else if( distsurf > sn%pmax .or. .not.defined( sn%psiffp2 ) ) then
               pt => pt_cell_sput_x(:,ielem1)
            else
               pt => pt_cell_sput_back(:,ielem1)
            end if
         else
            if( distsurf > sn%pmax .or. .not.defined( sn%psiffp2 ) ) then
               pt => pt_cell_x_x(:,ielem1)
            else
               pt => pt_cell_x_back(:,ielem1)
            end if
         end if

      else

         if( lsput ) then
            if( distsurf > dam%dnear ) then
               pt => pt_region_x_x(:,id,ielem1)
            else if( distsurf > sn%pmax .or. .not.defined( sn%psiffp2 ) ) then
               pt => pt_region_sput_x(:,id,ielem1)
            else
               pt => pt_region_sput_back(:,id,ielem1)
            end if
         else
            if( distsurf > sn%pmax .or. .not.defined( sn%psiffp2 ) ) then
               pt => pt_region_x_x(:,id,ielem1)
            else
               pt => pt_region_x_back(:,id,ielem1)
            end if
         end if

      end if
!
!  Calculate pmax by linear interpolation in pre-calculated table
!
      if( e < et(1) ) then
         pmax = pt(1)
      else if( e >= et(NT) ) then
         pmax = pt(NT)
      else
         it = 1 + facint * log( e/et(1) )
         fac = (e-et(it)) / (et(it+1)-et(it))
         pmax = (ONE-fac) * pt(it) + fac * pt(it+1)
      end if
!
!  Calculate ffp
!
      call getDomainProperty( id, dens_tot_=dens )
!          =================  i         o

      if( dens == ZERO ) then
         ffp = BIG
      else
         ffp = ONE / (PI * pmax**2 * dens)
      end if
!
! debug
!      write( 25, * ) e, pmax, ffp
! end-debug
!
   end subroutine evalFFP
!
!==============================================================================
!
   subroutine tabPmax( )
!             =======
!
!     Makes a table of pmax values as a function of energy and writes it to the file PMAX.
!     Depending on whether lsput and psiffp2 are set and how many regions are present,
!     the table has up to 4 columns:
!
!        1) energy
!        2) pmax values to be used deep inside the target
!        3) pmax values to be used near the surface (if lsput=t)
!        4) pmax values to be used outside the target (if psiffp2 is set)
!
!     Such tables are written for all target regions (if geom%lcell=f) and all atom species
!
!  Arguments:
!
!     None
!
!  Local Variables:
!
      integer(INT4) :: ir, ielem1, it, ic, nc
!
!  Start of Executable Code:
!  ------------------------
!
      if( .not.sn%lpmaxtab ) return

      open( 122, file='PMAX',  status='replace', form='formatted' )
!
!  Dynamic simulation
!
      if( geom%lcell ) then

         nc = 1                                 ! determine number of columns
         if( allocated( pt_cell_x_x ) )  nc = nc + 1
         if( allocated( pt_cell_sput_x ) )  nc = nc + 1
         if( allocated( pt_cell_x_back ) )  nc = nc + 1
         if( allocated( pt_cell_sput_back ) )  nc = nc + 1

         do ielem1 = 1, nelem

            write( 122, '(a,a)' )  &           ! write header
               '* Ion: ', elem(ielem1)%name

            write( 122, '(a)', advance='no' )  &
               '* energy [eV], pmax [A] '
            if( allocated( pt_cell_sput_back ) )  &
               write( 122, '(a)', advance='no' )  'if z<0, else'
            if( allocated( pt_cell_x_back ) )  &
               write( 122, '(a)', advance='no' )  'if z<0, else'
            if( allocated( pt_cell_sput_x ) )  &
               write( 122, '(a)', advance='no' )  'if z<dnear, else'
            if( allocated( pt_cell_x_x ) )  &
               write( 122, '(a)', advance='no' )  ' (bulk)'
            write( 122, '()' )

            write( 122, '(a)', advance='no' )  &
               '* (den[eV], psi[deg]) = '
            if( allocated( pt_cell_sput_back ) )  &
               write( 122, '(a,f7.1,a,f7.4,a)', advance='no' )  &
                  '(', min( dam%esurfmin, sn%denffp ), ',', sn%psiffp2, '), '
            if( allocated( pt_cell_x_back ) )  &
               write( 122, '(a,f7.1,a,f7.4,a)', advance='no' )  &
                  '(', sn%denffp, ',', sn%psiffp2, '), '
            if( allocated( pt_cell_sput_x ) )  &
               write( 122, '(a,f7.1,a,f7.4,a)', advance='no' )  &
                  '(', min( dam%esurfmin, sn%denffp ), ',', sn%psiffp, '), '
            if( allocated( pt_cell_x_x ) )  &
               write( 122, '(a,f7.1,a,f7.4,a)', advance='no' )  &
                  '(', sn%denffp, ',', sn%psiffp, ') '
            write( 122, '()' )

            write( 122, '(i15,i10)' ) nc, NT              ! write table
            write( 122, '(i15,4i10)' ) ( 0, ic=1,nc )

            do it = 1, NT

               write( 122, '(f15.3)', advance='no' )  et(it)
               if( allocated( pt_cell_sput_back ) )  &
                  write( 122, '(f10.5)', advance='no' )  pt_cell_sput_back(it,ielem1)
               if( allocated( pt_cell_x_back ) )  &
                  write( 122, '(f10.5)', advance='no' )  pt_cell_x_back(it,ielem1)
               if( allocated( pt_cell_sput_x ) )  &
                  write( 122, '(f10.5)', advance='no' )  pt_cell_sput_x(it,ielem1)
               if( allocated( pt_cell_x_x ) )  &
                  write( 122, '(f10.5)', advance='no' )  pt_cell_x_x(it,ielem1)
               write( 122, '()' )

            end do

         end do

      else
!
!  Non-dynamic simulation
!
         nc = 1                                 ! determine number of columns
         if( allocated( pt_region_x_x ) )  nc = nc + 1
         if( allocated( pt_region_sput_x ) )  nc = nc + 1
         if( allocated( pt_region_x_back ) )  nc = nc + 1
         if( allocated( pt_region_sput_back ) )  nc = nc + 1

         do ielem1 = 1, nelem

            do ir = 1, nr
               write( 122, '(a,a,a,a)' )  &    ! write header
                  '* Ion: ', elem(ielem1)%name, ', Material = ', trim(mat(ir)%name)

               write( 122, '(a)', advance='no' )  &
                  '* energy [eV], pmax '
               if( allocated( pt_region_sput_back ) )  &
                  write( 122, '(a)', advance='no' )  'if z<0, else'
               if( allocated( pt_region_x_back ) )  &
                  write( 122, '(a)', advance='no' )  'if z<0, else'
               if( allocated( pt_region_sput_x ) )  &
                  write( 122, '(a)', advance='no' )  'if z<dnear, else'
               if( allocated( pt_region_x_x ) )  &
                  write( 122, '(a)', advance='no' )  ' (bulk)'
               write( 122, '()' )

               write( 122, '(a)', advance='no' )  &
                  '* (den[eV], psi[deg]) = '
               if( allocated( pt_region_sput_back ) )  &
                  write( 122, '(a,f7.1,a,f7.4,a)', advance='no' )  &
                     '(', min( dam%esurfmin, sn%denffp ), ',', sn%psiffp2, '), '
               if( allocated( pt_region_x_back ) )  &
                  write( 122, '(a,f7.1,a,f7.4,a)', advance='no' )  &
                     '(', sn%denffp, ',', sn%psiffp2, '), '
               if( allocated( pt_region_sput_x ) )  &
                  write( 122, '(a,f7.1,a,f7.4,a)', advance='no' )  &
                     '(', min( dam%esurfmin, sn%denffp ), ',', sn%psiffp, '), '
               if( allocated( pt_region_x_x ) )  &
                  write( 122, '(a,f7.1,a,f7.4,a)', advance='no' )  &
                     '(', sn%denffp, ',', sn%psiffp, ') '
               write( 122, '()' )

               write( 122, '(i15,i10)' ) nc, NT              ! write table
               write( 122, '(i15,4i10)' ) ( 0, ic=1,nc )

               do it = 1, NT

                  write( 122, '(f15.3)', advance='no' )  et(it)
                  if( allocated( pt_region_sput_back ) )  &
                     write( 122, '(f10.5)', advance='no' )  pt_region_sput_back(it,ir,ielem1)
                  if( allocated( pt_region_x_back ) )  &
                     write( 122, '(f10.5)', advance='no' )  pt_region_x_back(it,ir,ielem1)
                  if( allocated( pt_region_sput_x ) )  &
                     write( 122, '(f10.5)', advance='no' )  pt_region_sput_x(it,ir,ielem1)
                  if( allocated( pt_region_x_x ) )  &
                     write( 122, '(f10.5)', advance='no' )  pt_region_x_x(it,ir,ielem1)
                  write( 122, '()' )

               end do

            end do

         end do

      end if

   end subroutine tabPmax
!
!=======================================================================
!
   subroutine getPmaxTable( pt, et, ielem1, ielem2, denffp, psiffp, eminpsi )
!             ============  o   i      i       i      i       i        i
!
!     Calculate a table of pmax values corresponding to a table of energy values for
!     given ion and target atoms species (ielem1, ielem2)
!
!  Arguments:
!
      real(REAL8), intent(out) :: pt(:)         ! impact parameter table
      real(REAL8), intent(in) :: et(:)          ! energy table
      integer(INT4), intent(in) :: ielem1       ! index of ion
      integer(INT4), intent(in) :: ielem2       ! index of target atom
      real(REAL8), intent(in) :: denffp         ! minimum energy transfer
      real(REAL8), intent(in) :: psiffp         ! minimum scattering angle
      real(REAL8), intent(in) :: eminpsi        ! energy down to which psiffp is applied
!
!  Local Variables:
!
      integer(INT4) :: it
      real(REAL8) :: psirad, rnorm, enormf, elim
      real(REAL8) :: pnmin, pnpsi, pnmax, pnlim
      real(REAL8) :: pmaxpsi
      real(REAL8) :: m1, m2
      real(REAL8) :: elim1, elim2, elim3, elim4, elim5, elim6
      real(REAL8) :: fr1val
!
!  Start of Executable Code:
!  ------------------------
!
      psirad = psiffp * PI/180

      rnorm = 0.8853 * 0.529 / ( real(elem(ielem1)%z)**0.23 + real(elem(ielem2)%z)**0.23 )

      enormf = 14.39979 * elem(ielem1)%z * elem(ielem2)%z / rnorm

      pnmax = sn%pmaxmax / rnorm
      pnmin = sn%pmaxmin / rnorm

      m1 = elem(ielem1)%m
      m2 = elem(ielem2)%m

      elim = (denffp / psirad**2) * m2 / m1
      fr1val = psirad * elim / enormf
      pnlim = fr1inv( fr1val, pnmin, pnmax )    ! pnlim is artificially restricted to
!             ======    i       i      i        !    [pnmin,pnmax], but we don't care

      fr1val = fr1( pnmax )
!              ===    i
      elim1 = (m1 / m2) / denffp * (fr1val * enormf) ** 2
      elim2 = (ONE / psirad) * (fr1val * enormf)

      fr1val = fr1( pnmin )
!              ===    i
      elim5 = (m1 / m2) / denffp * (fr1val * enormf) ** 2
      elim6 = (ONE / psirad) * (fr1val * enormf)

      if( defined( eminpsi ) ) then

         elim4 = min( elim6, max( elim2, eminpsi ) )  ! guarrantees pnmin<pnpsi<pnmax

         fr1val = psirad * elim4 / enormf
         pnpsi = fr1inv( fr1val, pnmin, pnmax )
!                ======    i       i      i
         pmaxpsi = pnpsi * rnorm

         elim3 = (m1 / m2) / denffp * (fr1val * enormf) ** 2

      else

         pmaxpsi = sn%pmaxpsi
         pnpsi = pmaxpsi / rnorm

         fr1val = fr1( pnpsi )
!                 ===    i
         elim3 = (m1 / m2) / denffp * (fr1val * enormf) ** 2
         elim4 = (ONE / psirad) * (fr1val * enormf)

      end if
!
!  Loop over table indices
!
      do it = 1, size( et )

         if( elim < elim3 ) then

            if( et(it) < elim1 ) then
               pt(it) = sn%pmaxmax
            else if( et(it) < elim5 ) then
               fr1val = sqrt( et(it) * denffp * m2 / m1 ) / enormf
               pt(it) = rnorm * fr1inv( fr1val, pnmin, pnmax )
!                               ======    i       i      i
            else
               pt(it) = sn%pmaxmin
            end if

         else if( elim < elim5 ) then

            if( et(it) < elim1 ) then
               pt(it) = sn%pmaxmax
            else if( et(it) < elim3 ) then
               fr1val = sqrt( et(it) * denffp * m2 / m1 ) / enormf
               pt(it) = rnorm * fr1inv( fr1val, pnpsi, pnmax )
!                               ======    i       i      i
            else if( et(it) < elim4 ) then
               pt(it) = pmaxpsi
            else if( et(it) < elim ) then
               fr1val = psirad * et(it) / enormf
               pt(it) = rnorm * fr1inv( fr1val, pnlim, pnpsi )
!                               ======    i       i      i
            else if( et(it) < elim5 ) then
               fr1val = sqrt( et(it) * denffp * m2 / m1 ) / enormf
               pt(it) = rnorm * fr1inv( fr1val, pnmin, pnlim )
!                               ======    i       i      i
            else
               pt(it) = sn%pmaxmin
            end if

         else

            if( et(it) < elim1 ) then
               pt(it) = sn%pmaxmax
            else if( et(it) < elim3 ) then
               fr1val = sqrt( et(it) * denffp * m2 / m1 ) / enormf
               pt(it) = rnorm * fr1inv( fr1val, pnpsi, pnmax )
!                               ======    i       i      i
            else if( et(it) < elim4 ) then
               pt(it) = pmaxpsi
            else if( et(it) < elim6 ) then
               fr1val = psirad * et(it) / enormf
               pt(it) = rnorm * fr1inv( fr1val, pnmin, pnpsi )
!                               ======    i       i      i
            else
               pt(it) = sn%pmaxmin
            end if

         end if                                 

      end do

   end subroutine getPmaxTable
!
!=======================================================================
!
   function fr1( arg ) result( value )
!           ===   i              o
!
!     Calculates the derivate of the screening function dependent term of the 
!     continuum potential
!
!  Arguments:
!
      real(REAL8), intent(in) :: arg            ! argument of the function to be calculated
      real(REAL8) :: value                      ! function value
!
!  Local Variables:
! 
      integer(INT4), parameter :: NCOEF = 4
      real(REAL8) :: a(NCOEF) = [ 0.1818, 0.50983, 0.2802, 0.02817 ]
      real(REAL8) :: b(NCOEF) = [ 3.2   , 0.9423 , 0.4028, 0.2016  ]

      integer(INT4) :: i
      real(REAL8) :: argb, valb
!
!  Start of Executable Code
!  ========================
!
      value = ZERO

      do i=1,NCOEF
         argb = b(i) * arg
         call mbes1( argb, valb )
!             =====   i     o
         value = value + a(i) * b(i) * valb
      end do

   end function fr1
!
!==============================================================================
!
   function fr1inv( arg, valmin, valmax ) result( value )
!           ======   i     i       i                o
!
!     Inverse function of fr1. We rely on fr1 being monotonically decreasing. 
!
!  Arguments:
!
      real(REAL8), intent(in) :: arg            ! argument (= fr1 value)
      real(REAL8), intent(in) :: valmin         ! lower bound to the result value
      real(REAL8), intent(in) :: valmax         ! upper bound to the result value
      real(REAL8) :: value                      ! result value (= fr1 argument)
!
!  Local Variables:
!
      real(REAL8) :: v, v1, v2, f, f1, f2
!
!  Start of Executable Code:
!  ------------------------
!
      if( valmax < valmin )  error stop 'fr1inv: valmax < valmin'

      v1 = valmin
      v2 = valmax
      f1 = fr1( v1 )
!          ===   i
      f2 = fr1( v2 )
!          ===   i

      do

         v = HALF * (v1 + v2)
         if( v2 - v1 < abs(v)*EPSMB )  exit

         f = fr1( v )
!            ===  i
         if( f < arg ) then
            v2 = v
            f2 = f
         else
            v1 = v
            f1 = f
         end if
         
      end do

      value = v

   end function fr1inv
!
!=======================================================================
!
end module Sn_FFP_mod
