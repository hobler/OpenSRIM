module Sn_Scatter_mod
!
!  Evaluate scattering event with possibly electronic stopping
!
!  Used Modules:
!  ============
!
   use Const_mod,          only: INT4, REAL8,  &
                                 ZERO, ONE, TWO, FOUR, HALF, QUART,  &
                                 EPSMA, BIG, PI
   use Vector_mod,         only: operator(.dot.), abs

   use Types_mod,          only: Proj_type

   use Param_Setup_mod,    only: nelem, lsput
   use Param_Elem_mod,     only: elem
   use Param_Ion_mod,      only: ion
   use Param_Xtal_mod,     only: xtal
   use Param_Dam_mod,      only: dam
   use Param_Sn_mod,       only: sn
   use Param_Se_mod,       only: se
   use Param_Out_mod,      only: out

   use Se_mod,             only: evalSe
   use Sn_Apsis_mod,       only: initApsis, apsis
   use Sn_Cmangm_mod,      only: initCmangm, cmangm
   use Sn_Cmangl_mod,      only: initCmangl, cmangl
   use Sn_CMscatter_mod,   only: cmAngle
   use Sn_Screen_mod,      only: initScreen, screen
   use Sn_TimeInt_mod,     only: initTimeInt, timeint
   use Target_mod,         only: getDomainIndex
   use Target_Surface_mod, only: getDistSurface
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
   public :: initScatter            ! initialize parameters
!            ===========
   public :: scatter                ! evaluate scattering event
!            =======
!
!  Private Procedures:
!  ==================
!
   private :: apsis_est             ! calculate apsis of collision (unused)
!             =========
   private :: apsis_acc             ! calculate apsis of collision (unused)
!             =========
!
!  Private Data:
!  ============
!
   real(REAL8), allocatable :: rnorm(:,:)    ! (nelem,nelem) scaling length
   real(REAL8), allocatable :: enorm(:,:)    ! (nelem,nelem) scaling energy
   real(REAL8), allocatable :: rm12(:,:)     ! (nelem,nelem) m1/m2
   real(REAL8), allocatable :: denfac(:,:)   ! (nelem,nelem) 
   real(REAL8), allocatable :: dirfac(:,:)   ! (nelem,nelem) 
!
!=======================================================================
!
contains
!
!=======================================================================
!
   subroutine initScatter ()
!             ===========
!
!     Initialize parameters needed for scattering angle evaluation
!
!     This subroutine must be called from a single thread
!
!  Arguments:
!
!     None
!
!  Local Variables:
!
      integer(INT4) :: ielem1, ielem2
!
!  Start of Executable Code:
!  ------------------------
!
      allocate( rnorm( nelem, nelem ) ) 
      allocate( enorm( nelem, nelem ) ) 
      allocate( rm12( nelem, nelem ) )
      allocate( denfac( nelem, nelem ) )
      allocate( dirfac( nelem, nelem ) )

      do ielem1 = 1, nelem
         do ielem2 = 1, nelem
            rm12(ielem2,ielem1) = elem(ielem1)%m / elem(ielem2)%m
            denfac(ielem2,ielem1) = FOUR * rm12(ielem2,ielem1) / (rm12(ielem2,ielem1)+ONE)**2
            dirfac(ielem2,ielem1) = TWO / (rm12(ielem2,ielem1) + ONE)
            rnorm(ielem2,ielem1) = 0.4685 /  &
               (real(elem(ielem1)%z)**0.23 + real(elem(ielem2)%z)**0.23)
            enorm(ielem2,ielem1) = 14.39979 * elem(ielem1)%z * elem(ielem2)%z /  &
               rnorm(ielem2,ielem1) * (ONE + rm12(ielem2,ielem1))
         enddo
      enddo
!
! begin test
      call initTimeInt( )
!          ===========
! end test

      if( any( sn%lquad(:,:) ) ) then
         call initScreen( )
!             ==========
         call initApsis( )
!             =========
      end if

      if( sn%lmehl ) then
         call initCmangm( )
!             ==========
      else
         call initCmangl( )
!             ==========
      end if

   end subroutine initScatter
!
!=======================================================================
!
   subroutine scatter( proj, posc, ffp, pmax, dee, cospsi, lrecoil, recoil, p )
!             =======   io    i    io    i     o     o       io       io    i
!
!     Evaluate scattering event including electronic stopping during the preceeding 
!     free flight path and during collision
!
!     Note: It is possible that on exit proj%e < 0
!
!  Arguments:
!
      type(Proj_type), intent(inout) :: proj    ! projectile
      real(REAL8), intent(in) :: posc(:)        ! collision point (= where p is normal to ffp)
      real(REAL8), intent(inout) :: ffp         ! free flight path:
                                                !    on entry: |posc-proj%pos|
                                                !    on exit:  |proj%pos-proj%posold|
      real(REAL8), intent(in) :: pmax           ! maximum impact parameter used during 
                                                !    selection of collision partners 
      real(REAL8), intent(out) :: dee           ! electronic energy loss
      real(REAL8), intent(out) :: cospsi        ! cosine of scattering angle
      logical, intent(inout) :: lrecoil         ! indicates whether there is a collision
                                                !    (=recoil); may be set to .false. if
                                                !    collision is aborted
      type(Proj_type), intent(inout) :: recoil  ! recoil
      real(REAL8), intent(in) :: p              ! impact parameter: |recoil%pos-posc|
!
!  Local Variables:
!
      logical :: ltime                          ! flag indicating use of time integral

      integer(INT4) :: ielem1, ielem2

      real(REAL8) :: dirn(3), dirp(3), dirrl(3)
      real(REAL8) :: posrn(3)
      real(REAL8) :: r0

      real(REAL8) :: fac, dirrn, sinps, cosps, x12, x1, x2
      real(REAL8) :: dirnn, denom, deeloc, dedxnl, en, rn0, pn, ef, de
      real(REAL8) :: f, sinthh2, costhh, tanps, tfac, theta, e1
      real(REAL8) :: cotthh, tau
!
!  Start of Executable Code:
!  ------------------------
!
      ielem1 = proj%ielem                       ! abbreviation
!
!  Trajectory cutoff energy
!
      if( lsput .and. proj%distsurf < dam%dnear ) then
         ef = min( sn%ef, dam%esurfmin )
         if( proj%distsurf > HALF*dam%dnear )  &   ! smooth transition close to dnear
               ef = ef + (proj%distsurf/(HALF*dam%dnear) - ONE) * (sn%ef - ef)
      else
         ef = sn%ef
      end if
!
!  Electronic energy loss
!
      if( lrecoil ) then

         ielem2 = recoil%ielem

         if( se%lapsis ) then                   ! estimate apsis of collision
            en = proj%e / enorm(ielem2,ielem1)
            pn = p / rnorm(ielem2,ielem1)
            call apsis_est( en, pn, rn0 )
!                =========  i   i    o
            r0 = rn0 * rnorm(ielem2,ielem1)
         else
            r0 = p
         end if

         call evalSe( dee, dedxnl, deeloc, proj, posc, ffp, pmax, ef, r0, ielem2 )
!             ======   o     o       o      i     i    io    i    i   i     i

      else

         call evalSe( dee, dedxnl, deeloc, proj, posc, ffp, pmax, ef )
!             ======   o     o       o      i     i    io    i    i
!
      end if
!
!  Check if energy drops below ef due to electronic energy loss alone
!
      if( proj%e - dee < ef ) then
         lrecoil = .false.                      ! abort collision if cutoff energy is reached
         dee = min( dee, proj%e )
      end if
!
!  If there is no collision partner, there is no deflection and we are done
!  (note: subroutine can be called with lrecoil=.false. already) 
      
      if( .not.lrecoil ) then

         proj%e = proj%e - dee
         proj%posold(:) = proj%pos(:)
         proj%pos(:) = proj%pos(:) + ffp * proj%dir(:)
         call getDomainIndex( proj%pos(:), proj%id )
!             ==============      i           io
         call getDistSurface( proj%pos(:), proj%distsurf )
!             ==============       i             o
         if( sn%lapsisef )  proj%emin = proj%e
         proj%ffpmin = ZERO
         cospsi = ONE

         return

      end if
!
!  Determine deflection and nuclear energy loss
!
      dirn(:) = proj%dir(:)

      if( sn%linel ) then
         e1 = proj%e - (dee - deeloc)
      else
         e1 = proj%e - (dee - HALF*deeloc)      ! not sure if this can be justified
      end if
!
!  Get the deflection angle in the center-of-mass system
!
      en = e1 / enorm(ielem2,ielem1)
      pn = p / rnorm(ielem2,ielem1)

      if( sn%lmehl ) then
         call apsis( en, pn, ielem1, ielem2, rn0 )
!             =====  i   i     i       i      o
         call cmangm( theta, tau, en, pn, rn0, ielem1, ielem2 )
!             ======    o     o   i   i    i     i       i
         if( theta < EPSMA ) then
            cotthh = BIG
         else
            cotthh = tan( HALF*(PI-theta) )
         end if
      else if( sn%lquad(ielem2,ielem1) ) then
         if( sn%quad(ielem2,ielem1) == 'Mehler' ) then
            call cmAngle( en, pn, ielem1, ielem2, theta, method='MehlerSymm',  &
!                =======  i   i     i       i       o              i
                          nabsc=sn%nabsc(ielem2,ielem1) )
!                                     i
         else if( sn%quad(ielem2,ielem1) == 'Legendre' ) then
            call cmAngle( en, pn, ielem1, ielem2, theta, method='LegendreChiCorr',  &
!                =======  i   i     i       i       o              i
                          nabsc=sn%nabsc(ielem2,ielem1) )
!                                     i
         end if
         if( theta < EPSMA ) then
            cotthh = BIG
         else
            cotthh = tan( HALF*(PI-theta) )
         end if
      else
         call cmangl( cotthh, tau, en, pn, ielem1, ielem2 ) 
!             ======    o      o   i   i     i       i
      end if

      tau = tau * rnorm(ielem2,ielem1)
!
!  Calculate nuclear energy loss
!
      sinthh2 = ONE / (ONE + cotthh**2)
      if( sn%linel ) then
         f = sqrt( max( ZERO, ONE - deeloc/e1 * ( ONE + rm12(ielem2,ielem1) ) ) )
         tfac = f*sinthh2 + QUART*(ONE-f)**2
      else
         f = ONE
         tfac = sinthh2
      end if
      recoil%e = denfac(ielem2,ielem1) * e1 * tfac
!
!  Calculate kinetic energy at the apsis of collision
!
      if( sn%lapsisef ) then
         costhh = sqrt( ONE - sinthh2 )
         call apsis_acc( en, pn, ielem1, ielem2, rn0 )
!             =========  i   i     i       i      o
         proj%emin = e1 * ( (ONE/(rm12(ielem2,ielem1)+ONE) * pn/rn0)**2 +  &
                            HALF*denfac(ielem2,ielem1) * pn/rn0 * costhh +  &
                            (rm12(ielem2,ielem1) / (rm12(ielem2,ielem1)+ONE))**2 )
         !print *, 'emin=', proj%emin
      end if
!
!  Get the new directions of the ion and the recoils
!
                                                ! direction of the impact parameter
      if( p == ZERO ) then
         dirp(:) = proj%dir(:)
      else
         dirp(:) = (recoil%pos(:)-posc(:)) / p
      end if
                                                ! relative momentum of the recoil
      if( sn%linel ) then
         tanps = 2*f*cotthh / ((1-f)*cotthh**2+(1+f))
      else
         tanps = cotthh
      end if
      denom = sqrt( ONE + tanps**2 )
      sinps = tanps / denom
      cosps = ONE / denom
      if( sn%linel ) then
         fac = sqrt( tfac )
      else
         fac = cosps
      end if
      dirrl(:) = dirfac(ielem2,ielem1) * fac * (cosps*proj%dir(:) + sinps*dirp(:))
                                                ! change momentum of the ion
      dirn(:) = dirn(:) - dirrl(:)

      if( dam%lrcoil .or. any( out%ltra(:) ) ) then
                                                ! scale recoil direction vector
         dirrn = abs( dirrl(:) ) 
!                ===    i
         if( dirrn == ZERO ) then
            recoil%dir(:) = dirrl(:)
         else
            recoil%dir(:) = dirrl(:) / dirrn
         end if
      end if
!
!  Recoil position and trajectory turning point
!
      ltime = proj%e < sn%etime
      ltime = ltime .and. cotthh /= ZERO .and. f /= ZERO

      if( ltime ) then
         if( .not.sn%lmehl .and.  &             ! make sure time integral is
             ( sn%lquad(ielem2,ielem1) .or.  &  ! calculated from ZBLuniv
              sn%scatfile(ielem2,ielem1) /= 'scattab/ZBLuniv' ) ) then
            call timeint( tau, en, pn )
!                =======   o   i   i
            tau = tau * rnorm(ielem2,ielem1)
         end if
         x12 = p / cotthh
         x2 = (ONE+ONE/f) * elem(ielem1)%m/(elem(ielem1)%m+elem(ielem2)%m) * (x12 - tau)
         posrn(:) = recoil%pos(:) + x2 * proj%dir(:)
         x1 = x2 - x12                          ! x1 corresponds to x1p of the Eckstein book
         if( sn%lnegx1 )  x1 = min( ZERO, x1 )  ! avoid positive x1
         if( sn%lposffp )  x1 = max( -ffp, x1 ) ! avoid negative ffp
         ffp = ffp + x1
      else
         posrn(:) = recoil%pos(:)
         x1 = ZERO
      end if
!
!  Scale projectile direction vector to unit length and calculate scattering angle
!
      dirnn = abs( dirn(:) )
!             ===
      if( dirnn == ZERO ) then
         dirn(:) = proj%dir(:)
         cospsi = ONE
      else
         dirn(:) = dirn(:) / dirnn
         cospsi = proj%dir(:) .dot. dirn(:)
      end if
!
!  New ion position ("turning point" of the trajectory)
!
      proj%posold(:) = proj%pos(:)

      if( ltime ) then
         proj%pos(:) = proj%pos(:) + ffp * proj%dir(:)
         proj%ffpmin = (recoil%pos(:) - proj%pos(:)) .dot. dirn(:)
         proj%ffpmin = max( ZERO, proj%ffpmin )
      else
         proj%pos(:) = posc(:)
         proj%ffpmin = ZERO
      end if
!
!  The following code was used until 2016:
!
!      if( ltime ) then
!         proj%posold(:) = proj%pos(:)
!         proj%pos(:) = proj%pos(:) + ffpl * proj%dir(:)
!         if( abs( cospsi ) > xtal%cpsolm ) then
!            ffp2 = (recoil%pos(:) - proj%pos(:)) .dot. dirn(:)
!            proj%pos(:) = proj%pos(:) + ffp2 * dirn(:)
!            ffpl = abs( proj%pos(:) - proj%posold(:) )
!!                  ===
!         else
!            ffpl = abs( ffpl )
!!                  ===
!         end if
!      else
!         proj%posold(:) = proj%pos(:)
!         proj%pos(:) = posc(:)
!         ffpl = ffp
!      end if
!
      call getDomainIndex( proj%pos(:), proj%id )
!          ==============      i           io
      call getDistSurface( proj%pos(:), proj%distsurf )
!          ==============      i              o
!
!  New energy and direction of the projectile
!
      dee = dee + dedxnl * x1
      de = dee + recoil%e
      if( de > proj%e ) then                    ! avoid negative energies that could occur
         recoil%e = recoil%e * proj%e / de      !     for linel=f
         dee = dee * proj%e / de
         proj%e = ZERO
      else
         proj%e = proj%e - de
      end if

      proj%dir(:) = dirn(:)

      recoil%pos(:) = posrn(:)

      if( ltime ) then                          ! if recoil position has changed due to time
                                                !    integral, recalculate id and distsurf
         call getDomainIndex( recoil%pos(:), recoil%id )
!             ==============        i            io
         call getDistSurface( recoil%pos(:), recoil%distsurf )
!             ==============         i               o
      end if
!
! test
!      if( e > ion%e ) then
!         print *,  'e,dee,den,z=',e,dee,(recoil%e,i=1,nc),z
!         error stop 'e>ion%e in SCATTER'
!      end if
! test-end
!
!  Store recoil energy and direction as properties of the recoil origin
!
      recoil%origin%e = recoil%e
      recoil%origin%dir(:) = recoil%dir(:)

   end subroutine scatter
!
!=======================================================================
!
   subroutine apsis_est( e, p, r0 )
!             =========  i  i  o
!
!     Calculate the distance of closest approach (the apsis) in the collision
!     approximating the ZBL screening function by a linear, a 1/R, and a 1/R^3
!     function. Worst case inaccuracy is 8%.
!
!  Arguments:
!
      real(REAL8), intent(in) :: e              ! scaled energy
      real(REAL8), intent(in) :: p              ! scaled impact parameter
      real(REAL8), intent(out) :: r0            ! scaled apsis
!
!  Local Variables:
!
      real(REAL8), parameter :: k2 = 0.38       ! factor of the 1/R part
      real(REAL8), parameter :: k3 = 7.2        ! factor of the 1/R^3 part
      real(REAL8), parameter :: k1 = ONE/(FOUR*k2)
      real(REAL8), parameter :: r12sq = (TWO*k2)**2
      real(REAL8), parameter :: r23sq = k3 / k2

      real(REAL8) :: psq, psqh, r0sq
!
!  Start of Executable Code:
!  ------------------------
!
      psq = p**2
      psqh = HALF * psq
      r0sq = psqh + sqrt( psqh**2 + k3/e )

      if( r0sq > r23sq ) then
         r0 = sqrt(r0sq)
         return
      end if

      r0sq = psq + k2/e

      if( r0sq > r12sq ) then
         r0 = sqrt(r0sq)
         return
      end if

      r0 = ( ONE + sqrt( ONE + FOUR*e*(e+k1)*psq ) ) / ( TWO*(e+k1) )

   end subroutine apsis_est
!
!=======================================================================
!
   subroutine apsis_acc( e, p, ielem1, ielem2, r0 )
!             =========  i  i    i       i     o
!
!     Calculate the distance of closest approach (the apsis) in the collision
!     using Newton iterations.
!
!  Arguments:
!
      real(REAL8), intent(in) :: e              ! scaled energy
      real(REAL8), intent(in) :: p              ! scaled impact parameter
      integer(INT4), intent(in) :: ielem1       ! index of projectile
      integer(INT4), intent(in) :: ielem2       ! index of recoil
      real(REAL8), intent(out) :: r0            ! scaled apsis
!
!  Local Variables:
!
      real(REAL8), parameter :: c = 17.43       ! value good for all potentials
                                                !    weaker than the Moliere
                                                !    at large distance
      real(REAL8) :: f, f1, phie, dr0, p2h, p2, phi, phi1

!
!  Start of Executable Code:
!  ------------------------
!
!  Initial solution
!
      p2 = p**2
      p2h = HALF * p2
      r0 = sqrt( p2h + sqrt( p2h**2 + c/e ) )
      call screen( p, ielem1, ielem2, phi, phi1 )
!          ======  i    i       i      o    o
      if( phi < c/r0**3 ) then
         phie = phi/(TWO*e)
         r0 = phie + sqrt( phie**2+p2 )
      end if
!
!  Newton iteration loop
!
      do
         call screen( r0, ielem1, ielem2, phi, phi1 )
!             ======  i     i       i      o    o
         f = (r0**2 - r0*phi/e - p2)

         if( f < ZERO ) then
            dr0 = dr0/2                ! avoid negative values
            r0 = r0 + dr0
         else
            f1 = (2*r0 - (phi + r0*phi1)/e)
            dr0 = f / f1
            r0 = r0 - dr0
         end if

         if( abs( dr0/r0 ) < 1.e-6 )  exit
      end do

   end subroutine apsis_acc
!
!=======================================================================
!
end module Sn_Scatter_mod
