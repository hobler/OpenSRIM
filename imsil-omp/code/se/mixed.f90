module Se_Mixed_mod
!
!  Handle mixed nonlocal and impact parameter dependent electronic stopping
!
!  Used Modules:
!  ============
!
   use Const_mod,       only: INT4, REAL8,  &
                              ZERO, ONE, TWO, FOUR, HALF

   use Param_Setup_mod, only: nelem
   use Param_Elem_mod,  only: elem
   use Param_Se_mod,    only: se

   use Se_Loc_mod,      only: setSeScreen, evalCorXnl, evalSeLoc, initSeLoc
   use Se_Rand_mod,     only: evalSeRand, initSeRand
   use Se_Straggle_mod, only: evalWe, initWe
   use Random_mod,      only: ranno
   use Target_mod,      only: getDomainDensities
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
   public :: initSeMixed            ! initialize parameters needed for considering local and
!            ===========            !    nonlocal contributions
   public :: evalSeMixed            ! evaluate the electronic energy loss considering local and
!            ===========            !    nonlocal contributions
!
!  Private Procedures:
!  ==================
!
   private :: evalSeMixed1          ! evaluate the electronic energy loss according to a mixed
!             ============          !     local/nonlocal model given the random stopping power
!
!  Private Data:
!  ============
!
   logical :: init = .false.

   real(REAL8), allocatable :: qnl(:,:,:)       ! (nnlm-1,nelem,nelem)
   real(REAL8), allocatable :: senl(:,:,:)      ! (nnlm,nelem,nelem)
!
!=======================================================================
!
contains
!
!=======================================================================
!
   subroutine initSeMixed( )
!             ===========
!
!     Initializes the se structure corresponding to Oen Robinson
!     electronic stopping. Calls initSeRand
!
!  Arguments:
!
!     None
!
!  Local Variables:
!
      logical :: first
      integer(INT4) :: ielem1, ielem2, inl, nnlm
      real(REAL8) :: p, se1m, se1p, se1, se2, e1
!
!  Start of Executable Code:
!  ------------------------
!
      init = .true.
!
!  Call random stopping power initialization routine
!
      call initSeRand( )
!          ==========
!
!  Call electronic straggling initialization routine
!
      call initWe( )
!          ======
!
!  Allocate memory
!
      nnlm = size( se%xnl, dim=1 )
      allocate( qnl( max( 1, nnlm-1 ), nelem, nelem ) )
      allocate( senl( nnlm, nelem, nelem ) )
!
!  Calculate derived quantities
!
      do ielem1 = 1, nelem
         do ielem2 = 1, nelem
!
!  If there is only one value, qnl is assumed zero. This will be treated as
!  a special case in evalSeMixed
!
            if( se%nxnl(ielem2,ielem1) == 1 ) then
               qnl(1,ielem2,ielem1) = ZERO
               cycle
            end if
!
!  proportional to E^q
!
            select case( se%xnlmodel(ielem2,ielem1) )
            case( 'e' )

               if( any( se%enl(1:se%nxnl(ielem2,ielem1)-1,ielem2,ielem1) == ZERO ) )  goto 961

               do inl = 1, se%nxnl(ielem2,ielem1) - 1

                  if( se%enl(inl+1,ielem2,ielem1) == ZERO ) then

                     call evalSeRand( se%enl(inl,ielem2,ielem1), ielem1, ielem2, se1 )
!                         ==========              i                i       i      o
                                                ! calculate position of stopping power maximum
                     if( se1 == ZERO ) then
                                                ! probably Se has been turned off => value
                                                !    cannot be calculated but does not matter
                        se%enl(inl+1,ielem2,ielem1) = 1000 * se%enl(inl,ielem2,ielem1)
                     else
                        call findSeMax( ielem1, ielem2, se%enl(inl+1,ielem2,ielem1) )
!                            =========    i       i                       o
                                                ! calculate stopping powers
                     end if
                     call evalSeRand( se%enl(inl+1,ielem2,ielem1), ielem1, ielem2, se2 )
!                         ==========                i                i       i      o
                                                ! calculate power of Se at lower energy
                     call evalSeRand( HALF*se%enl(inl,ielem2,ielem1), ielem1, ielem2, se1m )
!                         ==========             i                        i       i      o
                     call evalSeRand( TWO*se%enl(inl,ielem2,ielem1 ), ielem1, ielem2, se1p )
!                         ==========               i                      i      i      o
                     if( se1m == ZERO .or. se1p == ZERO ) then
                        qnl(inl,ielem2,ielem1) = &    ! emergency exit
                             log( se%xnl(inl+1,ielem2,ielem1) / se%xnl(inl,ielem2,ielem1) ) /  &
                             log( se%enl(inl+1,ielem2,ielem1 ) / se%enl(inl,ielem2,ielem1 ) )
                        cycle
                     end if
                     p = log( se1p / se1m ) / log( FOUR )
                                                ! calculate qnl
                     qnl(inl,ielem2,ielem1) =  &
                          p * log( se%xnl(inl+1,ielem2,ielem1) / se%xnl(inl,ielem2,ielem1) ) / &
                          log( se2 / se1 )

                  else

                     qnl( inl, ielem2, ielem1 ) =  &
                          log( se%xnl(inl+1,ielem2,ielem1) / se%xnl(inl,ielem2,ielem1) ) /  &
                          log( se%enl(inl+1,ielem2,ielem1) / se%enl(inl,ielem2,ielem1) )

                  end if

               end do
!
!  proportional to Se^q
!
            case( 'se' )
                                                ! replace enl=0 by position of stopping power
               first = .true.                   !    maximum
               do inl = 1, se%nxnl(ielem2,ielem1)
                  if( se%enl(inl,ielem2,ielem1) /= ZERO )  cycle
                  if( .not.first )  goto 962
                  first = .false.
                                                ! check if Se has been turned off
                  if( inl == 1 ) then
                     e1 = 1000
                  else
                     e1 = se%enl(inl-1,ielem2,ielem1)
                  end if
                  call evalSeRand( e1, ielem1, ielem2, se1 )
!                      ==========  i     i       i      o

                  if( se1 == ZERO ) then
                                             ! probably Se has been turned off => stopping
                                             !    power maximum cannot be determined, but does
                                             !    not matter
                     se%enl(inl+1,ielem2,ielem1) = 1000 * e1
                  else
                                             ! calculate position of stopping power maximum
                     call findSeMax( ielem1, ielem2, se%enl(inl,ielem2,ielem1) )
!                         =========     i       i                 o
                  end if
               end do
                                                ! calculate senl
               do inl = 1, se%nxnl(ielem2,ielem1)
                  call evalSeRand( se%enl(inl,ielem2,ielem1), ielem1, ielem2,  &
!                      ==========            i                   i      i
                                   senl(inl, ielem2, ielem1) )
!                                              o
               end do
                                                ! calculate qnl
               do inl = 1, se%nxnl(ielem2,ielem1) - 1
                  qnl(inl,ielem2,ielem1) =  &
                       log( se%xnl(inl+1,ielem2,ielem1) / se%xnl(inl,ielem2,ielem1) ) /  &
                       log( senl(inl+1,ielem2,ielem1) / senl(inl,ielem2,ielem1) )
               end do

            end select

         end do
      end do
!
!  screening length
!
      call initSeLoc()
!          =========

      return

961   error stop 'initSeMixed: only last enl=0 allowed for xnlmodel=''e'''
962   error stop 'initSeMixed: only one enl=0 allowed for xnlmodel=''se'''

   contains
!
!-----------------------------------------------------------------------
!
      subroutine findSeMax( ielem1, ielem2, emax )
!                =========     i     i      o
!
!     Find the energy where the random stopping power is maximum
!
!  Arguments:
!
         integer(INT4), intent(in) :: ielem1, ielem2
         real(REAL8), intent(out) :: emax
!
!  Local Variables:
!
         real(REAL8) :: e1, e2, e3, e4, e5, se1, se2, se3, se4, se5
         integer(INT4) :: i
!
!  Start of Executable Code:
!  ------------------------
!
!  Find three se values such that the second one is the biggest
!
         e2 = 1
         call evalSeRand( e2, ielem1, ielem2, se2 )
!             ==========  i      i     i      o
         e3 = 2
         call evalSeRand( e3, ielem1, ielem2, se3 )
!             ==========  i      i     i      o
         do
            e1 = e2
            se1 = se2
            e2 = e3
            se2 = se3
            e3 = 2*e2
            call evalSeRand( e3, ielem1, ielem2, se3 )
!                ==========  i      i     i      o
            if( se3 < se2 )  exit
         end do
!
!  Find three se values with equidistant e values such that the second
!  se value is the biggest
!
         e4 = e3
         se4 = se3
         e3 = (e2+e4)/2
         call evalSeRand( e3, ielem1, ielem2, se3 )
!             ==========  i      i     i      o
         if( se3 > se2 ) then
            e1 = e2
            se1 = se2
            e2 = e3
            se2 = se3
            e3 = e4
            se3 = se4
         end if
!
!  Subdivide the intervals and discard two of the five values such that the
!  remaining second value is the biggest
!
         do, i = 1, 10
            e5 = e3
            se5 = se3
            e3 = e2
            se3 = se2
            e2 = (e1+e3)/2
            call evalSeRand( e2, ielem1, ielem2, se2 )
!                ==========  i      i     i      o
            e4 = (e3+e5)/2
            call evalSeRand( e4, ielem1, ielem2, se4 )
!                ==========  i      i     i      o
            if( se4 > se3 ) then
               e1 = e3
               se1 = se3
               e2 = e4
               se2 = se4
               e3 = e5
               se3 = se5
            else if( se3 > se2 ) then
               e1 = e2
               se1 = se2
               e2 = e3
               se2 = se3
               e3 = e4
               se3 = se4
            end if
         end do

         emax = e2

      end subroutine findSeMax
!
!-----------------------------------------------------------------------
!
   end subroutine initSeMixed
!
!=======================================================================
!
   subroutine evalSeMixed( dee, dedxnl, deeloc, e, ielem1, id, ffp, pmax, ef, p_, ielem2_ )
!             ===========   o     o       o     i    i     i   io    i    i   i      i
!
!     Evaluate the electronic energy loss according to Oen and Robinson with a 
!     nonlocal component, considering electronic straggling
!
!     It is assumed that ielem2_ is present if and only if p_ is present
!
!  Arguments:
!
      real(REAL8), intent(out) :: dee                 ! electronic energy loss [eV]
      real(REAL8), intent(out) :: dedxnl              ! nonlocal energy loss per path length
                                                      !    [eV/A] (at end of ffp if iterated)
      real(REAL8), intent(out) :: deeloc              ! local energy loss [eV]
      real(REAL8), intent(in) :: e                    ! ion energy [eV]
      integer(INT4), intent(in) :: ielem1             ! element index of the ion
      integer(INT4), intent(in) :: id                 ! domain index
      real(REAL8), intent(inout) :: ffp               ! free flight path [A] (reduced on exit if
                                                      !    projectile stops)
      real(REAL8), intent(in) :: pmax                 ! maximum impact parameter [A]
      real(REAL8), intent(in) :: ef                   ! trajectory cutoff energy [eV]
      real(REAL8), optional, intent(in) :: p_         ! impact parameter [A]
      integer(INT4), optional, intent(in) :: ielem2_  ! element index of recoil
!
!  Local Variables:
!
      logical :: calc_loc
      integer(INT4) :: nelemr, ielemr, ielemr2
      integer(INT4) :: ielem(nelem)
      real(REAL8) :: de, ee, eee, seval, weval, rand, dffp, ffprest, straggle
      real(real8) :: dens(nelem), facstrag(nelem), dedxnlr(nelem)
!
!  Start of Executable Code:
!  ------------------------
!
      dee = ZERO
      dedxnl = ZERO
      deeloc = ZERO

      ee = e               ! ee will be decreased as electronic stopping accumulates
      ffprest = ffp        ! remaining free flight path
      straggle = ZERO      ! electronic energy loss straggling
!
!  Get the property of the target
!
      call getDomainDensities( id, nelemr, ielem, dens_=dens )
!          ==================  i     o       o        o
!
!  Loop over free flight path increments
!
      do
!
!  Sum up the nonlocal electronic stopping weighted with the abundance of the
!  atoms in the target.
!  If p_ is present, calculate local electronic stopping for use once the ffp loop has been
!  completed.
!
         do ielemr = 1, nelemr

            eee = ee - se%eminse(ielem(ielemr),ielem1)    ! energy to be used in Se calc.

            if( eee < ZERO ) then
               seval = ZERO
            else
               call evalSeRand( eee, ielem1, ielem(ielemr), seval )   ! random stopping power
!                   ==========   i     i          i           o
            end if

            if( seval <= ZERO .or. se%semodel(ielem(ielemr),ielem1) == 'none' ) then
               dedxnlr(ielemr) = ZERO
               if( present( p_ ) ) then
                  if( ielem(ielemr) == ielem2_ )  deeloc = ZERO
               end if
               facstrag(ielemr) = ZERO
            else
               if( present( p_ ) ) then
                  calc_loc = ( ielem(ielemr) == ielem2_ )
               else
                  calc_loc = .false.
               end if
               if( calc_loc ) then
                  ielemr2 = ielemr
                  call evalSeMixed1( seval, eee, pmax, dens(ielemr), ielem1, ielem2_,  &
!                      ============    i     i    i         i          i       i
                                     dedxnlr(ielemr), deeloc, p_ )
!                                          o            o     i
               else
                  call evalSeMixed1( seval, eee, pmax, dens(ielemr), ielem1,  &
!                      ============    i     i    i         i          i
                                     ielem(ielemr), dedxnlr(ielemr) )
!                                          i              o
               end if

               call evalWe( weval, eee, ielem1, ielem(ielemr) )
!                   ======    o     i       i           i
               facstrag(ielemr) = weval / seval
            end if

         end do

         dedxnl = sum( dedxnlr(1:nelemr) )
!
!  Decrement energy by nonlocal energy loss; if necessary, reduce free flight path and iterate
!
         de = dedxnl * ffprest

         if( de <= 0.01*ee ) then
            ee = ee - de
            straggle = straggle + sum( facstrag(1:nelemr) * dedxnlr(1:nelemr) ) * ffprest
            exit                                ! end of ffp reached
         else
            de = 0.01*ee
            ee = ee - de

            dffp = de / dedxnl
            ffprest = ffprest - dffp

            straggle = straggle + sum( facstrag(1:nelemr) * dedxnlr(1:nelemr) ) * dffp
         end if

         if( ee < ef ) then
            dee = e - ee
            ffp = ffp - ffprest
            return                     ! don't do straggle, this could cause ee>ef
         end if

      end do
!
!  Consider local electronic energy loss
!
      if( deeloc > ZERO ) then
         ee = ee - deeloc
         straggle = straggle + facstrag(ielemr2) * deeloc
      end if
!
!  Deterministic electronic energy loss
!
      dee = e - ee
!
!  Consider electronic straggling
!
      if( straggle > 0 ) then

         call ranno( rand )
!             =====   o
         de = sqrt(straggle) * rand
         !print *, 'dee, de=', dee, de
         !de = max( de, -dee )                   ! avoid negative dee
         deeloc = deeloc * ( ONE + de/dee )     ! distributing straggle equally among local and
         dedxnl = dedxnl * ( ONE + de/dee )     !    nonlocal stopping is not exact but
         dee = dee + de                         !    probably accurate enough

      end if

   end subroutine evalSeMixed
!
!=======================================================================
!
   subroutine evalSeMixed1( seval, e, pmax, dens, ielem1, ielem2, dedxnl, deeloc_, p_)
!             ============    i     i   i     i      i       i       o       o     i
!
!     Evaluate the mixed local and nonlocal electronic energy loss, given the random electronic
!     stopping power
!
!  Arguments:
!
      real(REAL8), intent(in) :: seval          ! random stopping power [eV A^2]
      real(REAL8), intent(in) :: e              ! projectile energy [eV]
      real(REAL8), intent(in) :: pmax           ! maximum impact parameter used when
                                                !    calling the subroutine
      real(REAL8), intent(in) :: dens           ! density of ielem2
      integer(INT4), intent(in) :: ielem1       ! element index of the projectile
      integer(INT4), intent(in) :: ielem2       ! element index of recoil

      real(REAL8), intent(out) :: dedxnl           ! nonlocal electronic energy loss per path
                                                   !    length [eV/A]
      real(REAL8), optional, intent(out) :: deeloc_! local electronic energy loss [eV]
      real(REAL8), optional, intent(in) :: p_      ! impact parameter [A]
!
!  Local variables:
!
      real(REAL8) :: xnl, xloc, dxnl, fac, corxnl
!
!  Start of Executable Code:
!  ------------------------
!
      if( .not.init )  goto 901
!
!  nonlocal and local fraction
!
      if( se%locmodel(ielem2,ielem1) == 'none' ) then
         xnl = ONE
      else
         call evalXnl( xnl, seval, e, ielem1, ielem2 )
!             =======   o     i    i    i       i
      end if
      xloc = ONE - xnl

      call setSeScreen( e, xnl, ielem1, ielem2 )
!          ===========  i   i     i       i
!
!  pmax correction
!
      if( xnl /= ONE ) then
         call evalCorXnl( corxnl, pmax, ielem1, ielem2 )
!             ==========    o      i      i       i
         dxnl = xloc * corxnl
         xnl = xnl + dxnl
      end if
!
!  calculate nonlocal and local energy loss
!
      dedxnl = xnl * dens * seval

      if( present( p_ ) ) then
         if( xnl == ONE ) then
            deeloc_ = ZERO
         else
            call evalSeLoc( fac, p_, ielem1, ielem2 )
!                =========   o   i     i       i
            deeloc_ = xloc * seval * fac
         end if
      end if

      return

901   error stop 'evalSeMixed1: not initialized'

   end subroutine evalSeMixed1
!
!=======================================================================
!
   subroutine evalXnl( xnl, seval, e, ielem1, ielem2 )
!             =======   o     i    i    i       i
!
!     Evaluate nonlocal and local fraction of electronic energy loss.
!     Note that the sum of the two does not equal unity because of the maximum impact parameter
!     correction.
!
!  Arguments:
!
      real(REAL8), intent(out) :: xnl           ! nonlocal fraction
      real(REAL8), intent(in) :: seval          ! random stopping power [eV A^2]
      real(REAL8), intent(in) :: e              ! ion energy [eV]
      integer(INT4), intent(in) :: ielem1       ! element index of the projectile
      integer(INT4), intent(in) :: ielem2       ! element index of recoil
!
!  Local variables:
!
      integer(INT4) :: inl
      real(REAL8) :: qqnl
!
!  Start of Executable Code:
!  ------------------------
!
      if( .not.init )  goto 901

      do inl = 1, se%nxnl(ielem2,ielem1) - 2
         if( e < se%enl(inl+1,ielem2,ielem1) )  exit
      end do
                                                ! simplify notation
      qqnl = qnl(inl,ielem2,ielem1)

      if( qqnl == ZERO ) then
         xnl = se%xnl(inl,ielem2,ielem1)
      else
         select case( se%xnlmodel(ielem2,ielem1) )
            case( 'e' )
               xnl = se%xnl(inl,ielem2,ielem1) * ( e / se%enl(inl,ielem2,ielem1) ) ** qqnl
            case( 'se' )
               xnl = se%xnl(inl,ielem2,ielem1) * ( seval / senl(inl,ielem2,ielem1) ) ** qqnl
            case default
               error stop 'evalSeMixed: unknown model'
         end select
      end if

      xnl = min( xnl, ONE )

      return

 901  error stop 'evalSeMixed: not initialized'

   end subroutine evalXnl
!
!=======================================================================
!
   subroutine evalSeOR_old( dee, dedxnl, deeloc, e, ielem1, id, ffp, pmax, p_, ielem2_ )
!             ============   o     o       o     i    i     i    i    i    i      i
!
!     Evaluate the electronic energy loss according to Oen and Robinson with a
!     nonlocal component considering electronic straggling
!
!  Arguments:
!
      real(REAL8), intent(out) :: dee                 ! electronic energy loss [eV]
      real(REAL8), intent(out) :: dedxnl              ! nonlocal energy loss per path length
                                                      !    [eV/A]
      real(REAL8), intent(out) :: deeloc              ! local energy loss [eV]
      real(REAL8), intent(in) :: e                    ! ion energy [eV]
      integer(INT4), intent(in) :: ielem1             ! element index of the ion
      integer(INT4), intent(in) :: id                 ! domain index
      real(REAL8), intent(in) :: ffp                  ! free flight path [A]
      real(REAL8), intent(in) :: pmax                 ! maximum impact parameter [A]
      real(REAL8), optional, intent(in) :: p_         ! impact parameter [A]
      integer(INT4), optional, intent(in) :: ielem2_  ! element index of recoil
!
!  Local Variables:
!
      logical :: first
      integer(INT4) :: nelemr, ielemr
      real(REAL8) :: ee, dedxnl1, deeloc1, seval, dee1, dedxnl_sum
      real(REAL8) :: facstrag, facstrag1, weval, deenl, rand, sigma, dee_sum
      real(real8) :: prob(nelem), dens(nelem)
      integer(INT4) :: ielem(nelem)
!
!  Start of Executable Code:
!  ------------------------
!
!  Get the property of the target
!
      call getDomainDensities( id, nelemr, ielem, prob, dens_=dens )
!          ==================  i     o       o     o        o

      ee = e - sum( prob(1:nelemr) * se%eminse(ielem(1:nelemr),ielem1) )
      first = .true.
      dedxnl_sum = ZERO
      dee_sum = ZERO
      deeloc = ZERO

      do                                        ! this loop will be executed at most twice
                                                !    (see comment towards the end of the loop)
         if( ee <= ZERO ) then
                                                ! no electronic stopping for negative energy
            if( first ) then
               dee_sum = ZERO
               dedxnl_sum = ZERO
               deeloc = ZERO
               facstrag = ZERO
            else
               dee_sum = HALF * dee1
               dedxnl_sum = HALF * dedxnl1
               deeloc = HALF * deeloc1
               facstrag = HALF * facstrag1
            end if

            exit

         end if
!
!  Sum up the nonlocal electronic stopping weighted with the abundance of the
!  atoms in the target. If p_ is present, calculate local electronic stopping
!  according to the ion species of the recoil.
!
         do ielemr = 1, nelemr

            call evalSeRand( ee, ielem1, ielem(ielemr), seval )   ! random stopping power
!                ==========   i     i         i           o

            if( seval <= ZERO ) then               ! no electronic energy loss if seval=0
               facstrag = ZERO
               dee_sum = ZERO
               dedxnl_sum = ZERO
            else
               call evalWe( weval, e, ielem1, ielem(ielemr) )
!                   ======    o    i      i           i
               facstrag = weval / seval

               if( present( p_ ) .and. ielem(ielemr) == ielem2_ ) then  ! local component of Se

                  call evalSeMixed1( seval, ee,  &
!                      ============    i    i
!                                     ffp,   ! this was required in the old version
!                                      i
                                     pmax, dens(ielemr), ielem1, ielem2_,  &
!                                     i         i          i        i
!                                     dee,   ! this was returned in the old version
!                                      o
                                     dedxnl, deeloc, p_ )
!                                      o       o     i
                  dee_sum = dee_sum + deeloc
               else
                  call evalSeMixed1( seval, ee,  &
!                      ============    i    i
!                                     ffp,
!                                      i
                                     pmax, dens(ielemr), ielem1, ielem(ielemr),  &
!                                     i         i          i             i
!                                     dee,
!                                      o
                                     dedxnl )
!                                      o
               end if

               dedxnl_sum = dedxnl_sum + prob(ielemr) * dedxnl

            end if

         end do     ! loop over atoms

         dee_sum = dee_sum + dedxnl_sum * ffp

         if( first ) then
            deenl = dedxnl_sum * ffp
            if( deenl > 0.01*e ) then           ! if the energy loss exceeds one percent of
               first = .false.                  !    the particle energy, compute the energy
               facstrag1 = facstrag             !    loss from the mean value of the stopping
               dedxnl1 = dedxnl_sum             !    at the beginning and at the end of the
               deeloc1 = deeloc                 !    free flight path.
               dee1 = dee_sum
               ee = ee - dee_sum
               dee_sum = ZERO
               dedxnl_sum = ZERO
               cycle
            end if
         else
            facstrag = HALF * (facstrag1 + facstrag)
            dedxnl_sum = HALF * (dedxnl1 + dedxnl_sum)
            deeloc = HALF * (deeloc1 + deeloc)
            dee_sum = HALF * (dee1 + dee_sum)
         end if

         exit

      end do

      dedxnl = dedxnl_sum
      dee = dee_sum
      if( abs( dee - ( dedxnl*ffp + deeloc ) ) > 1e-4 * abs(dee) ) then
         print *, dee, dedxnl, deeloc
         error stop 'evalSeMixed: error dee/=loc+nl'
      end if
!
!  Consider electronic straggling
!
      if( facstrag > 0 ) then

         deenl = dedxnl * ffp
         sigma = sqrt( facstrag * deenl )
         call ranno( rand )
!             =====   o
         deenl = deenl + sigma * rand
         dedxnl = deenl / ffp

         sigma = sqrt( facstrag * deeloc )
         call ranno( rand )
!             =====   o
         deeloc = deeloc + sigma * rand

         dee = deenl + deeloc

      end if

   end subroutine evalSeOR_old
!
!=======================================================================
!
end module Se_Mixed_mod
