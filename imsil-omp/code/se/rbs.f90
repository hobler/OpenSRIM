module Se_Rbs_mod
!
!  Approximate the paths of the backscattered particles by a straight
!  trajectories in a random medium. 
!
!  Used modules:
!  ============
!
   use Const_mod,           only: REAL8, INT4,  &
                                    PI, ZERO, BIG, EPSMB
   use Se_Rand_mod,         only: evalSeRand

   use Param_Setup_mod,     only: nelemIon, nelem
   use Param_Geom_mod,      only: geom
   use Param_Elem_mod,      only: elem
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
   public :: initRevPath 
!            ===========
   public :: findRbsEnergy
!            =============
!
!  Private procedures:
!  ==================
!
   private :: inversePath
!             ===========
   private :: drawLine
!             ========
!
!  Private Data:
!  ============
!
   logical :: init = .false.

   real(REAL8), allocatable :: ener(:,:)  ! energy [eV]
   real(REAL8), allocatable :: depth(:,:) ! depth [A]
   real(REAL8), allocatable :: sepow(:,:) ! random stopping power
   integer(INT4) :: il                    ! low bound of e, se ,rp
   integer(INT4) :: iu                    ! upper bound of e, se, rp
!
!=========================================================================================
!
contains
!
!=========================================================================================
!
   subroutine initRevPath( ielem, fi, eIon )
!             ===========    i    i    i 
!
!     Intialize the random electronic stopping and reverse path structure
!
!  Arguments:
!
      integer(INT4), intent(IN) :: ielem      ! ion atom element index
      real(REAL8), intent(IN) :: fi           ! sc. angle for reverse path initializ.
      real(REAL8), intent(IN) :: eIon         ! initial ion energy
!
!  Local variables
!
      integer(INT4) :: ielem2, status, i
      real(REAL8) :: eres
!
!  Start of executable code
!  ------------------------
!
      init = .true.
!
!  Allocate the structure and initialize the values
!
      il = 0          ! mininum index
      iu = 150        ! maximum index 
      eres = 10.      ! energy resolution (log. scale)

      allocate ( ener( nelemIon+1:nelem, il:iu ), STAT = status)
      if (status .ne. 0) error stop 'initRevPath: e not allocated'
      ener(:,:) = ZERO

      allocate ( depth( nelemIon+1:nelem, il:iu ), STAT = status)
      if (status .ne. 0) error stop 'initRevPath: z not allocated'
      depth(:,: ) = ZERO

      allocate ( sepow( nelemIon+1:nelem, il:iu ), STAT = status)
      if( status .ne. 0 ) error stop 'initRevPath: se not allocated'
      sepow(:,:) = ZERO
!
!  Initialization.  eres < e <  eIon with equal spacing in logarithmic scale
!
      do ielem2 = nelemIon+1, nelem 

         ener(ielem2, il) = ZERO          ! we will need it when outgoing proj. is
         sepow(ielem2, il) = ZERO         ! at the end of his energy

         do i = il+1, iu

            ener(ielem2,i) = eres * (eIon / eres)**(i/real(iu)) 

            call evalSeRand( ener(ielem2,i), ielem, ielem2, sepow(ielem2,i) )
!                ==========        i           i      i          o

            sepow(ielem2,i) = sepow(ielem2,i) * elem(ielem2)%dens
         end do
!
!  Calculate outgoing path of bacscattered particle
!
         call inversePath ( ielem2, eIon, PI-fi )
!             ===========      i     i      i   
      end do
!
   end subroutine initRevPath
!
!==========================================================================================
!
   subroutine inversePath( ielem2, eFin, fi )
!             ===========    i      i    i 
!
!     Simulate the inverse outgoing path of particle including energy loss (random). Here
!     we create a z- E dependece so the energy loss can be simply included by binary search
!     throug this dependence. Not every time to simulate the outgoing path and waste the
!     CPU time.
!     Here we suppose particle with zero energy which enters into target with angle (180-fi)
!     and simulate the energy increase (due to energy loss ) through the target. As a 
!     matter of fact, we are calculating the energy which the ion has to have at certain 
!     depth if it wants to get out from the target. 
!
!  Arguments:
!
      integer(INT4), intent(IN) :: ielem2  ! index of target atom
      real(REAL8), intent(IN) :: eFin      ! maximal energy for rev. path
      real(REAL8), intent(IN) :: fi        ! angle[deg] from z axis of outgoing path (180-sc.angle)
!
!  Local variables
! 
      real(REAL8) :: e, z, SE, dR
      integer(INT4) :: irs, imax, imin, im
!
!  Start of executable code
!  ------------------------
!
      dR = 100.0          ! 100 A for beginning
      depth( ielem2, il ) = ZERO
      e = EPSMB           ! initial energy of ion entering the target
      z = ZERO            ! initial position of ion entering the target
!
!  While e > 0 .OR. z > 0 loop
!
      invLoop: do ! while e < eFin

         if( e .gt. eFin ) exit invLoop !  Check the condition
!
!  find indices of appropriate position of e in randStop array by binary search
!
         imin = il
         imax = iu

111      im = ( imin + imax ) / 2
         if( e .gt. ener( ielem2, im ) ) then
            imin = im
         else
            imax = im
         end if
         if ( imax - imin > 1 ) goto 111
         irs = imin
!
!  Random stopping at the middle of energy interval
! 
         SE = ( sepow( ielem2, irs ) + sepow( ielem2, irs+1 ) ) / 2 
!
!  Cover interpolated region
!
         dR = ( ener( ielem2, irs + 1 ) - ener( ielem2, irs ) ) / SE

         if ( dR .lt. ZERO ) error stop  'InversePath: dR < 0.'

         e = e + SE * dR       !  Increase energy and depth (-energy loss)
         z = z + dR * cos( fi )
         depth( ielem2, irs + 1 ) = z !  Remebmer depth and energy

         if ( e .gt. eFin ) exit invLoop ! Check the condition

      end do invLoop !  End Loop

   end subroutine InversePath
!
!=========================================================================================
!
   subroutine findRbsEnergy( eScatt, ielem2, pos, eRbs, valRbs )
!             =============     i       i     i     o      o
!
!     Detect the energy of projectile when leaves the target on its outgoing way. 
!     In the same time that will be energy of particle detected in RBS detector.
!
!  Arguments:
!
      real(REAL8), intent(IN) :: eScatt    ! energy of ion after large angle scattering
      integer(INT4), intent(IN) :: ielem2  ! index of target atom
      real(REAL8), intent(IN) :: pos       ! depth (z) where collision happens
      real(REAL8), intent(OUT) :: eRbs     ! Energy of particle at detector
      logical, intent(OUT) :: valRbs       ! T - projectile reach the detector
                                           ! F - not enough energy to reach the detector
!
!  Local variables
! 
      integer(INT4) :: iPe, iRbs, imax, imin, im
      real(REAL8) :: z, e, zRbs, zP, a1, b1
!
!  Start of executable code:
!  ------------------------
!
      if( .not.init )  goto 901

      z = pos
      e = eScatt
!
!  Atom can be scattered before entering into target. Allow this possibility!
!
      if( z < geom%posif(0) .and. e > EPSMB ) then
         eRbs = e
         valRbs = .true.
         return
      end if
!
!  Find the energy box by binary search (iPe index)
!
      imin = il
      imax = iu

111   im = (imin + imax) / 2
      if( e > ener(ielem2, im) ) then
         imin = im
      else
         imax = im
      end if
      if( imax-imin > 1 )  goto 111
      iPe = imin
!
!  Interpolate z which corresponds to eProj
!
      call drawLine( depth(ielem2,iPe), ener(ielem2,iPe),  depth(ielem2,iPe+1),  & 
!          ========        i                  i                    i             
                     ener(ielem2,iPe+1), a1, b1 )
!                              i          o   o
      zP = ( e - b1 ) / a1                      ! e = a1 * z * b1 => z(eProj) = (eProj - b1)/a1
!
!  See do the ion have enough energy (according to determined box) to get out of sample
!
      if( z < zP ) then
!
!   From zp subtract real depth (z) to see what is the rest energy of ion
!
         zRbs = zP - z 
!
!  Now interpolate energy(eRbs) from zRbs. First find appropriate box for zRbs
!
         imin = il
         imax = iu

222      im = (imin + imax) / 2
         if( zRbs > depth(ielem2, im) ) then
            imin = im
         else
            imax = im
         end if
         if ( imax-imin > 1) goto 222
         iRbs = imin
!
!  Linear interpolation
!
         call drawLine( depth(ielem2, iRbs), ener(ielem2, iRbs),  depth(ielem2, iRbs+1), & 
!             ========           i                   i                      i                  
                        ener( ielem2, iRbs + 1 ), a1, b1 )
!                                     i            o   o

         eRbs = a1 * zRbs + b1
         valRbs = .true.
!
!  if z > required for ion (with given energy) then it will never get out of target
!
      else
         eRbs = e
         valRbs = .false.
      end if

      return

901   error stop 'findRbsEnergy: not initialized'

   end subroutine FindRbsEnergy
!
!=========================================================================================
!
   subroutine drawLine( x1, y1, x2, y2, a, b)
!             ========  i   i   i   i   o  o
!
!     Draw the line between two points given with his coordinates. Finds the coefficents.
!
!  Calling arguments:
!
      real(REAL8), intent(IN) :: x1,y1     ! coordinates of point 1
      real(REAL8), intent(IN) :: x2,y2     ! coordinates of point 2
      real(REAL8), intent(OUT) :: a, b     ! coefficients of line
!
!  Start of executable code
!  ------------------------
!
!  Check do points have same value
!
      if ( x1 == x2) then
         if ( y2 == y2 ) then
            write (*,*) 'drawLine: points are the same!'
            return
         else
            write (*,*) 'drawLine: vertical line'
            a = BIG
            return
         end if
      end if

      a = (y1 - y2) / (x1 - x2)
      b = (x1 * y2 - y1 * x2) / (x1 - x2)

   end subroutine drawLine
!
!===========================================================================================
!
end module Se_Rbs_mod
