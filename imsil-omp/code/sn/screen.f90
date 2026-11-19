module Sn_Screen_mod
!
!  Calculate the screening function and its derivative using the ZBL potential
!
!  Used Modules:
!  ============
!
   use Const_mod,       only: INT4, REAL8,  &
                              ZERO, ONE, TWO, THIRD, SMALL, BIG
   use Array_mod,       only: check_allocation
   use Open_mod
   use String_mod,      only: str
   use Math_Table_mod,  only: Table1d_type, initTable, interpolate
   use Undefine_mod,    only: defined

   use Param_Setup_mod, only: nelem
   use Param_Elem_mod,  only: elem
   use Param_Sn_mod,    only: sn
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
   public :: initScreen
   public :: screen
   public :: get_rmax
!
!  Private Data:
!  ============
!
   type(Table1d_type), allocatable :: phitab(:,:)

   type :: Coef_type
      logical :: lrmax
      real(REAL8), allocatable :: a(:)
      real(REAL8), allocatable :: b(:)
      real(REAL8) :: c
      real(REAL8) :: d
      real(REAL8) :: rmax
   end type Coef_type

   type(Coef_type), allocatable :: coef(:,:)

!
!=======================================================================
!
contains
!
!=======================================================================
!
   ! Populate native screening data for the current element configuration.
   ! Coefficient values are normalized here; callers pass reduced radii to screen.
   ! Reinitialization replaces prior screening storage for the selected setup.
   subroutine initScreen( )
!             ==========
!
!     Initialize the screening function.
!
!     The screening function is either the sum of exponentials with an
!     optional linear term
!
!        Phi(r) = sum_{i=1}^{n} a_i * exp( -b_i*r )
!                 [ + c - d*r/rmax ]
!
!     The radius r is expressed in the normalization used internally
!     by screen().
!
!     For the coefficient-file implementation the coefficients are
!     transformed according to the selected screening length.
!
!     In particular, the special screening length 'NLHlin' is chosen
!     such that a coefficient file containing the original NLHlin
!     coefficients is transformed to exactly the same normalization
!     used by the Python/Numba NLHlin implementation:
!
!
!        new because normalization in python / numba is used like this
!        
!        rnorm_NLHlin =
!           0.4685 / sqrt(sqrt(Z1) + sqrt(Z2))
!
!     This is achieved by setting
!
!        a = aZBL / rnorm_NLHlin
!
!     because the generic coefficient transformation below then gives
!
!        b    -> b * rnorm_NLHlin
!        rmax -> rmax / rnorm_NLHlin
!
!  Arguments:
!
!     None
!
!  Local Variables:
!
      character(256) :: errmsg, line
      integer(INT4) :: istat, ielem1, ielem2
      integer(INT4) :: z1, z2, i, j, n
      real(REAL8) :: aZBL, a
      real(REAL8) :: rnorm_NLHlin
      real(REAL8), allocatable :: rs(:), phis(:)
!
!  Start of Executable Code:
!  ------------------------
!
!  Allocate arrays
!
      if( any( defined( sn%scrfile(:,:) ) ) ) then
         if( allocated( phitab ) )  deallocate( phitab )
         allocate( phitab( nelem, nelem ), stat=istat, errmsg=errmsg )
         call check_allocation( 'initScreen', 'phitab', istat, errmsg )
      end if

      if( allocated( coef ) ) deallocate( coef )
      allocate( coef( nelem, nelem ), stat=istat, errmsg=errmsg )
      call check_allocation( 'initScreen', 'coef', istat, errmsg )

      do ielem1 = 1, nelem
         do ielem2 = 1, nelem

            if( defined( sn%scatfile(ielem2,ielem1) ) ) cycle
!
!  ZBL screening length
!
            aZBL = 0.4685_real8 / &
                   ( real(elem(ielem1)%z,REAL8)**0.23_real8 + &
                     real(elem(ielem2)%z,REAL8)**0.23_real8 )

!
!  Selected screening length
!
            select case( sn%scrlen(ielem2,ielem1) )

               case( '1' )

                  a = ONE

               case( 'ZBL' )

                  a = aZBL

               case( 'Firsov' )

                  a = 0.4685_real8 / &
                        ( sqrt(real(elem(ielem1)%z,REAL8)) + &
                          sqrt(real(elem(ielem2)%z,REAL8)) ) ** (TWO*THIRD)

               case( 'Lindhard' )

                  a = 0.4685_real8 / &
                      sqrt( &
                           real(elem(ielem1)%z,REAL8)**(TWO*THIRD) + &
                           real(elem(ielem2)%z,REAL8)**(TWO*THIRD) &
                      )

               case( 'NLHlin' )

!
!                 Python/Numba NLHlin normalization:
!
!                    rnorm_NLHlin =
!                       0.4685 / sqrt(sqrt(Z1) + sqrt(Z2))
!
!                 The generic scaling below uses
!
!                    b    = b_file * aZBL / a
!                    rmax = rmax_file * a / aZBL
!
!                 Therefore choosing
!
!                    a = aZBL / rnorm_NLHlin
!
!                 gives exactly
!
!                    b    = b_file * rnorm_NLHlin
!                    rmax = rmax_file / rnorm_NLHlin
!
!                 which is the normalization used by Python/Numba.
!
                  rnorm_NLHlin = &
                     0.4685_real8 / &
                     sqrt( &
                        sqrt(real(elem(ielem1)%z,REAL8)) + &
                        sqrt(real(elem(ielem2)%z,REAL8)) &
                     )

                  a = aZBL / rnorm_NLHlin

               case default

                  error stop 'initScreen: unknown scrlen'

            end select
!
!  Initialize coefficients or screening function table
!
            if( defined( sn%scrfun(ielem2,ielem1) ) ) then

               if( sn%scrfun(ielem2,ielem1) == 'ZBLuniv' ) then

                  allocate( coef(ielem2,ielem1)%a( 4 ), &
                            stat=istat, errmsg=errmsg )
                  call check_allocation( 'initScreen', 'coef%a', &
                                         istat, errmsg )

                  allocate( coef(ielem2,ielem1)%b( 4 ), &
                            stat=istat, errmsg=errmsg )
                  call check_allocation( 'initScreen', 'coef%b', &
                                         istat, errmsg )

                  coef(ielem2,ielem1)%a = &
                     [ 0.18175_real8, 0.50986_real8, &
                       0.28022_real8, 0.02817_real8 ]

                  coef(ielem2,ielem1)%b = &
                     [ 3.1998_real8, 0.94229_real8, &
                       0.4029_real8, 0.20162_real8 ]

                  coef(ielem2,ielem1)%b = &
                     coef(ielem2,ielem1)%b * aZBL / a

                  coef(ielem2,ielem1)%lrmax = .false.

               else if( sn%scrfun(ielem2,ielem1) == 'KrC' ) then

                  allocate( coef(ielem2,ielem1)%a( 3 ), &
                            stat=istat, errmsg=errmsg )
                  call check_allocation( 'initScreen', 'coef%a', &
                                         istat, errmsg )

                  allocate( coef(ielem2,ielem1)%b( 3 ), &
                            stat=istat, errmsg=errmsg )
                  call check_allocation( 'initScreen', 'coef%b', &
                                         istat, errmsg )

                  coef(ielem2,ielem1)%a = &
                     [ 0.190945_real8, 0.473674_real8, 0.335381_real8 ]

                  coef(ielem2,ielem1)%b = &
                     [ 0.278544_real8, 0.637174_real8, 1.919249_real8 ]

                  coef(ielem2,ielem1)%b = &
                     coef(ielem2,ielem1)%b * aZBL / a

                  coef(ielem2,ielem1)%lrmax = .false.

               else if( sn%scrfun(ielem2,ielem1) == 'Moliere' ) then

                  allocate( coef(ielem2,ielem1)%a( 3 ), &
                            stat=istat, errmsg=errmsg )
                  call check_allocation( 'initScreen', 'coef%a', &
                                         istat, errmsg )

                  allocate( coef(ielem2,ielem1)%b( 3 ), &
                            stat=istat, errmsg=errmsg )
                  call check_allocation( 'initScreen', 'coef%b', &
                                         istat, errmsg )

                  coef(ielem2,ielem1)%a = &
                     [ 0.35_real8, 0.55_real8, 0.10_real8 ]

                  coef(ielem2,ielem1)%b = &
                     [ 0.3_real8, 1.2_real8, 6.0_real8 ]

                  coef(ielem2,ielem1)%b = &
                     coef(ielem2,ielem1)%b * aZBL / a

                  coef(ielem2,ielem1)%lrmax = .false.

               end if

            else if( defined( sn%scrcoeffile(ielem2,ielem1) ) ) then

               call openFile( 11, sn%scrcoeffile(ielem2,ielem1), &
                              'formatted', 'read' )

               coef(ielem2,ielem1)%lrmax = .false.
               n = 4

               do
                  read( 11, '(a)' ) line

                  if( line(1:1) /= '#' ) exit

                  coef(ielem2,ielem1)%lrmax = &
                     coef(ielem2,ielem1)%lrmax .or. &
                     (index( line, 'rmax' ) /= 0)

                  i = index( line, 'nexp' )

                  if( i == 0 ) cycle

!
!                 The header has the form:
!
!                    # nexp=3
!
!                 Extract everything after "nexp".
!
                  line = adjustl( line(i+4:) )

!
!                 The first character must be '='.
!
                  if( len_trim(line) == 0 .or. line(1:1) /= '=' ) then

                     print *, 'initScreen: Ill-formatted header of file ''' // &
                              sn%scrcoeffile(ielem2,ielem1) // ''''

                     print *, 'WARNING: Expected nexp=<integer>: ', &
                              trim(line)

                     cycle

                  end if

!
!                 Remove '=' and read the integer directly.
!
                  line = adjustl( line(2:) )

                  read( line, *, iostat=istat ) n

                  if( istat /= 0 .or. n < 1 ) then

                     print *, 'initScreen: Ill-formatted header of file ''' // &
                              sn%scrcoeffile(ielem2,ielem1) // ''''

                     print *, 'WARNING: Illegal value of nexp: ', &
                              trim(line)

                  end if

               end do

               backspace( 11 )

               allocate( coef(ielem2,ielem1)%a( n ), &
                         stat=istat, errmsg=errmsg )
               call check_allocation( 'initScreen', 'coef%a', &
                                      istat, errmsg )

               allocate( coef(ielem2,ielem1)%b( n ), &
                         stat=istat, errmsg=errmsg )
               call check_allocation( 'initScreen', 'coef%b', &
                                      istat, errmsg )

               do

                  if( coef(ielem2,ielem1)%lrmax ) then

                     read( 11, *, iostat=istat ) z1, z2, &
                        ( coef(ielem2,ielem1)%a(i), &
                          coef(ielem2,ielem1)%b(i), i=1,n ), &
                        coef(ielem2,ielem1)%rmax

                     coef(ielem2,ielem1)%c = &
                        ONE - sum( coef(ielem2,ielem1)%a(:) )

                     coef(ielem2,ielem1)%d = &
                        sum( &
                           coef(ielem2,ielem1)%a(:) * &
                           exp( &
                              -coef(ielem2,ielem1)%b(:) * &
                               coef(ielem2,ielem1)%rmax &
                           ) &
                        ) + coef(ielem2,ielem1)%c

                  else

                     read( 11, *, iostat=istat ) z1, z2, &
                        ( coef(ielem2,ielem1)%a(i), &
                          coef(ielem2,ielem1)%b(i), i=1,n )

                     coef(ielem2,ielem1)%a(:) = &
                        coef(ielem2,ielem1)%a(:) / &
                        sum( coef(ielem2,ielem1)%a(:) )

                  end if

                  if( istat > 0 ) then

                     print *, 'Error during reading file ''' // &
                              trim( sn%scrcoeffile(ielem2,ielem1) ) // ''''

                     error stop

                  else if( istat < 0 ) then

                     print *, 'initScreen: Z1=' // &
                              str( elem(ielem1)%z ) // ', ' // &
                              'Z2=' // &
                              str( elem(ielem2)%z ) // &
                              ' not found on file ''' // &
                              trim( sn%scrcoeffile(ielem2,ielem1) ) // ''''

                     error stop

                  end if

!
!                 Apply the screening-length normalization.
!
                  coef(ielem2,ielem1)%b = &
                     coef(ielem2,ielem1)%b * aZBL / a

                  if( coef(ielem2,ielem1)%lrmax ) then
                     coef(ielem2,ielem1)%rmax = &
                        coef(ielem2,ielem1)%rmax * a / aZBL
                  end if

!
!                 Debug output for the Python/Numba NLHlin test case.
!
                  if( &
                     elem(ielem1)%z == 5 .and. &
                     elem(ielem2)%z == 14 &
                  ) then

                     print *
                     print *, '========== NLHlin loaded coefficients =========='
                     print *, 'n    = ', n
                     print *, 'a    = ', coef(ielem2,ielem1)%a
                     print *, 'b    = ', coef(ielem2,ielem1)%b
                     print *, 'c    = ', coef(ielem2,ielem1)%c
                     print *, 'd    = ', coef(ielem2,ielem1)%d
                     print *, 'rmax = ', coef(ielem2,ielem1)%rmax
                     print *, '================================================='

                  end if

                  if( &
                     min(z1,z2) == &
                     min(elem(ielem1)%z,elem(ielem2)%z) .and. &
                     max(z1,z2) == &
                     max(elem(ielem1)%z,elem(ielem2)%z) &
                  ) exit

               end do

               close( 11 )

            else if( defined( sn%scrfile(ielem2,ielem1) ) ) then

               call openFile( 11, sn%scrfile(ielem2,ielem1), &
                              'formatted', 'read' )

               n = 0

               do
                  read( 11, *, iostat=istat ) line
                  if( istat /= 0 ) exit
                  if( line == ' ' ) cycle
                  n = n + 1
               end do

               rewind( 11 )

               if( allocated( rs ) ) deallocate( rs, phis )

               allocate( rs(n), phis(n) )

               do i = 1, n
                  read( 11, *, iostat=istat ) rs(i), phis(i)

                  if( istat /= 0 ) exit

                  rs(i) = rs(i) * a / aZBL
               end do

               close( 11 )

               call initTable( &
                  phitab(ielem2,ielem1), &
                  rs(:), &
                  phis(:), &
                  cubic_=.true. &
               )

               if( phis(n) == ZERO ) then
                  coef(ielem2,ielem1)%lrmax = .true.
                  coef(ielem2,ielem1)%rmax = rs(n)
               else
                  coef(ielem2,ielem1)%lrmax = .false.
               end if

            else

               error stop 'initScreen: Unknown interatomic potential'

            end if

         end do
      end do

   end subroutine initScreen
!
!=======================================================================
!
   ! Evaluate phi and optionally its derivative with respect to reduced r.
   ! For finite-range analytic coefficients, r >= rmax returns exact zeros;
   ! the derivative may jump at that boundary even though phi is continuous.
   subroutine screen( r, ielem1, ielem2, phi, phi1_ )
!             ======  i    i       i      o     o
!
!     Calculate the screening function and derivative.
!
!  Arguments:
!
      real(REAL8), intent(in) :: r
      integer(INT4), intent(in) :: ielem1
      integer(INT4), intent(in) :: ielem2
      real(REAL8), intent(out) :: phi
      real(REAL8), intent(out), optional :: phi1_
!
!  Local Variables:
!
      integer(INT4) :: i
      real(REAL8) :: fexp
!
!  Start of Executable Code:
!  ------------------------
!
      if( defined( sn%scrfile(ielem2,ielem1) ) ) then

         if( r > phitab(ielem2,ielem1)%x( &
               size(phitab(ielem2,ielem1)%x) ) ) then

            phi = ZERO

            if( present(phi1_) ) phi1_ = ZERO

         else

            call interpolate( &
               phitab(ielem2,ielem1), &
               r, &
               phi, &
               phi1_ &
            )

         end if

      else

         phi = ZERO

         if( present(phi1_) ) phi1_ = ZERO

         do i = 1, size(coef(ielem2,ielem1)%a)

            fexp = exp( &
               -coef(ielem2,ielem1)%b(i) * r &
            )

            phi = phi + &
                  coef(ielem2,ielem1)%a(i) * fexp

            if( present(phi1_) ) then
               phi1_ = phi1_ - &
                       coef(ielem2,ielem1)%a(i) * &
                       coef(ielem2,ielem1)%b(i) * &
                       fexp
            end if

         end do

         if( coef(ielem2,ielem1)%lrmax ) then

            if( r < coef(ielem2,ielem1)%rmax ) then

               phi = phi + &
                     coef(ielem2,ielem1)%c - &
                     coef(ielem2,ielem1)%d * &
                     r / coef(ielem2,ielem1)%rmax

               if( present(phi1_) ) then
                  phi1_ = phi1_ - &
                          coef(ielem2,ielem1)%d / &
                          coef(ielem2,ielem1)%rmax
               end if

            else

               phi = ZERO

               if( present(phi1_) ) phi1_ = ZERO

            end if

         end if

      end if

   end subroutine screen
!
!=======================================================================
!
   real(REAL8) function get_rmax( ielem1, ielem2 )
!                       ========    i       i
!
!     Return the radius where the screening function drops to zero.
!     Return BIG (the native infinite-range sentinel) when no cutoff is set.
!
!  Arguments:
!
      integer(INT4), intent(in) :: ielem1
      integer(INT4), intent(in) :: ielem2
!
!  Local Variables:
!
!     None
!
!  Start of Executable Code:
!  ------------------------
!
      if( coef(ielem2,ielem1)%lrmax ) then
         get_rmax = coef(ielem2,ielem1)%rmax
      else
         get_rmax = BIG
      end if

   end function get_rmax
!
!=======================================================================
!
end module Sn_Screen_mod
