module Sn_Apsis_mod
!
!  Calculate the apsis of collision (minimum distance of approach) either by solving
!  its defining equation or by interpolation in a table (THE LATTER NOT YET IMPLEMENTED)
!
!  NOTE: The screening function is assumed to be initialized
!
!  Used Modules:
!  ============
!
   use Const_mod,       only: INT4, REAL8,  &
                              ZERO, ONE, HALF
   use Array_mod,       only: check_allocation
   use Stack_Real8_mod

   use Param_Setup_mod, only: nelem
   use Math_Table_mod,  only: Table1d_type, initTable, interpolate
   use Sn_Screen_mod,   only: screen
!
!  Default Type and Access:
!  =======================
!
   implicit none
   private
   save
!
!  Public Procedures
!  =================
!
   public :: initApsis                          ! initializations for the calculation of the
!            =========                          !   apsis of collision
   public :: apsis                              ! calculate apsis of collision
!            =====
!
!  Private Data
!  ============
!
   type(Table1d_type), allocatable :: r0_table(:,:)   ! table r0(e) for head-on collisions
!
!==============================================================================================
!
contains
!
!==============================================================================================
!
   subroutine initApsis( )
!             =========
!
!     Set up a table of apses (module variable r0_table) for head-on collisions as a function
!     of reduced energy e. They are the solutions of the equation
!
!        f(r0) = r0 - 1/e * Phi(r0) = 0
!
!     The equation is solved by Newton's method. The initial condition for the first point
!     e = emax is r0 = 1/emax, where emax is the maximum energy of the table. Then
!     e is halved, and the old solution is used as the initial condition for Newton's method
!     and so on, until the minimum energy emin is reached. This algorithm is guaranteed to
!     converge, since the intitial condition for r0 is always to the left of the solution, and
!     f'(r) > 0, f"(r) < 0.
!
!  Arguments:
!
!     None
!
!  Local Variables:
!
      real(REAL8), parameter :: emin=1e-8       ! min energy of table
      real(REAL8), parameter :: emax=1e4        ! max energy of table
      real(REAL8), parameter :: RELERRMAX=1e-6  ! max relative error for Newton

      logical :: isEmpty
      character(256) :: errmsg
      integer(INT4) :: itab, ntab, istat, ielem1, ielem2
      real(REAL8) :: r0, dr0, phi, phi1, f, f1, e, e1, relerr

      type(Real8Stack) :: r0stack, estack, e1stack
      real(REAL8), allocatable :: r0tab(:), etab(:), e1tab(:)
!
!  Start of Executable Code:
!  ------------------------
!
      if( allocated( r0_table ) )  deallocate( r0_table )
      allocate( r0_table( nelem, nelem ) )
!
!  Create stacks for r0, e and e1=d(e)/d(r0)
!
      call newReal8Stack( r0stack )
!          =============    o
      call newReal8Stack( estack )
!          =============    o
      call newReal8Stack( e1stack )
!          =============     o

      do ielem1 = 1, nelem
         do ielem2 = 1, nelem
!
!  First point corresponding to e = emax
!
            r0 = ONE/emax                          ! initial condition
            call screen( r0, ielem1, ielem2, phi, phi1 )
!                ======  i     i       i      o    o

            relerr = ONE
            do while( relerr > RELERRMAX )         ! Newton iterations
               f = phi/r0 - emax
               f1 = phi1/r0 - phi/r0**2
               dr0 = f/f1
               r0 = r0 - dr0
               call screen( r0, ielem1, ielem2, phi, phi1 )
!                   ======  i     i       i      o    o
               relerr = abs(dr0/r0)
            end do

            e1 = phi1/r0 - phi/r0**2

            call pushReal8Stack( r0, r0stack )     ! put values on stacks
!                ==============  i     io
            call pushReal8Stack( emax, estack )
!                ==============   i      io
            call pushReal8Stack( e1, e1stack )
!                ==============  i     io
!
!  Other points at e = emax * 2**(-k)
!
            e = emax

            do while( e > emin )                   ! loop over energies
               e = HALF * e

               ! use r0, phi, phi1 from last iteration as initial condition

               relerr = ONE
               do while( relerr > RELERRMAX )      ! Newton iterations
                  f = phi/r0 - e
                  f1 = phi1/r0 - phi/r0**2
                  dr0 = f/f1
                  r0 = r0 - dr0
                  call screen( r0, ielem1, ielem2, phi, phi1 )
!                      ======  i     i       i      o    o
                  relerr = abs(dr0/r0)
               end do

               e1 = phi1/r0 - phi/r0**2

               call pushReal8Stack( r0, r0stack )  ! put values on stacks
!                   ==============  i     io
               call pushReal8Stack( e, estack )
!                   ==============  i    io
               call pushReal8Stack( e1, e1stack )
!                   ==============  i     io
            end do
!
!  Convert stacks to table
!
            ntab = sizeOfReal8Stack( r0stack )
!                  ================     i
            if( allocated( r0tab ) )  deallocate( r0tab )
            allocate( r0tab( ntab ), stat=istat, errmsg=errmsg )
            call check_allocation( 'initApsis', 'r0tab', istat, errmsg )

            if( allocated( etab ) )  deallocate( etab )
            allocate( etab( ntab ), stat=istat, errmsg=errmsg )
            call check_allocation( 'initApsis', 'etab', istat, errmsg )

            if( allocated( e1tab ) )  deallocate( e1tab )
            allocate( e1tab( ntab ), stat=istat, errmsg=errmsg )
            call check_allocation( 'initApsis', 'e1tab', istat, errmsg )

            do itab = 1, ntab
               call popReal8Stack( r0stack, r0tab(itab), isEmpty )
!                   =============    io         o           o
               call popReal8Stack( estack, etab(itab), isEmpty )
!                   =============    io        o          o
               call popReal8Stack( e1stack, e1tab(itab), isEmpty )
!                   =============    io         o           o
            end do

            call initTable( r0_table(ielem2,ielem1), etab, r0tab, ONE/e1tab )
!                =========            o               i      i        i

         end do
      end do

   end subroutine initApsis
!
!==============================================================================================
!
   subroutine apsis( e, p, ielem1, ielem2, r0, dr0de_, dr0dp_, nn_ )
!             =====  i  i    i       i     o     o       o      o
!
!     Calculate the apsis of collision by solving the equation
!
!        f(r0) = r0 - 1/e * Phi(r0) - P**2/r0 = 0
!
!     The equation is solved by Newton's method with either r0=p or r0 from the table
!     calculated in initApsis, whatever is larger, as initial conditions. This initial
!     condition lies to the left of the solution (apart from interpolation errors),
!     which guarantees Newton's method to converge, if f'(r) > 0, f"(r) < 0. A sufficient 
!     but not necessary condition for the latter is Phi'(R) < 0, Phi"(R) > 0, which is
!     fulfilled for sum-of-exponentials screening functions.
!
!     Calculated values are saved and used in case of subsequent calls with same e and p.
!
!  Arguments:
!
      real(REAL8), intent(in) :: e                 ! energy [enorm]
      real(REAL8), intent(in) :: p                 ! impact parameter [rnorm]
      integer(INT4), intent(in) :: ielem1          ! index of projectile
      integer(INT4), intent(in) :: ielem2          ! index of recoil
      real(REAL8), intent(out) :: r0               ! apsis of collision [rnorm]
      real(REAL8), intent(out), optional :: dr0de_ ! d(r0)/d(e)
      real(REAL8), intent(out), optional :: dr0dp_ ! d(r0)/d(p)
      integer(INT4), intent(out), optional :: nn_  ! number of Newton iterations
!
!  Local Variables:
!
      real(REAL8), parameter :: RELERRMAX=1e-3  ! max relative error for Newton

      integer(INT4) :: n
      real(REAL8) :: phi, phi1, f, f1, dr0, denom, relerr

      logical :: same_e_p
      logical, save :: deriv_saved = .false.
      integer(INT4), save :: ielem1_save = 0, ielem2_save = 0
      real(REAL8), save :: e_save = ZERO, p_save = ZERO
      real(REAL8), save :: r0_save, dr0de_save, dr0dp_save
!$omp threadprivate( deriv_saved, ielem1_save, ielem2_save, e_save, p_save,  &
!$omp                r0_save, dr0de_save, dr0dp_save )
!
!  Start of Executable Code:
!  ------------------------
!
      if( .not.allocated( r0_table ) )  call initApsis( )  ! this is just needed for unittests
!                                            =========     ! might not work from parallel region

      same_e_p = ( ielem1 == ielem1_save .and. ielem2 == ielem2_save .and.  &
                   e == e_save .and. p == p_save )

      if( same_e_p ) then
!
!  Use saved apsis
!
         r0 = r0_save

      else
!
!  Calculate apsis
!
         ! Initial guess

         call interpolate( r0_table(ielem2,ielem1), e, r0 )
!             ===========             i             i  o

         r0 = max( r0, p )

         ! Newton iterations

         n = 0
         relerr = 1

         do while( relerr > RELERRMAX )

            n = n + 1
            if( n > 100 ) then
               print *, 'e, p=', e, p
               error stop 'apsis: Convergence error'
            end if

            call screen( r0, ielem1, ielem2, phi, phi1 )
!                ======  i     i       i      o    o

            f = r0 - p**2/r0 - phi/e
            f1 = 1 + p**2/r0**2 - phi1/e
            dr0 = - f / f1
            r0 = r0 + dr0

            relerr = abs(dr0/r0)

         end do
!
!  Save apsis
!
         r0_save = r0

      end if
!
!  Derivatives
!
      if( present( dr0de_ ) ) then

         if( .not.present( dr0dp_ ) )  &
            error stop 'apsis: Either none or both of dr0de, dr0dp must be passed.'

         if( same_e_p .and. deriv_saved ) then
!
!  Use saved derivatives
!
            dr0de_ = dr0de_save
            dr0dp_ = dr0dp_save

         else
!
!  Calculate derivatives
!
            call screen( r0, ielem1, ielem2, phi, phi1 )
!                ======  i     i       i      o    o
            denom = 1 + (p/r0)**2 - phi1/e
            dr0de_ = - phi/e**2 / denom
            dr0dp_ = 2*p/r0 / denom
!
!  Save derivatives
!
            dr0de_save = dr0de_
            dr0dp_save = dr0dp_
            deriv_saved = .true.

         end if

      else if( present( dr0dp_ ) ) then

         error stop 'apsis: Either none or both of dr0de, dr0dp must be passed.'

      else if( .not.same_e_p ) then

         deriv_saved = .false.

      end if

      if( present( nn_ ) )  nn_ = n

   end subroutine apsis
!
!==============================================================================================
!
end module Sn_Apsis_mod
