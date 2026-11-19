module Math_Quadrature_mod
!
!  Integrate a function in one variable from 0 to 1
!
!  Used Modules:
!  ============
!
   use Const_mod,             only: INT4, REAL8,  &
                                    ZERO, ONE, TWO, FOUR, PI
   use Array_mod,             only: check_allocation
   use Open_mod
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
   public :: quadLegendre              ! Gauss-Legendre quadrature
!            ============
   public :: quadMehlerSymm            ! Gauss-Mehler quadrature using symmetric extension
!            ==============            !   to (-1,1)
   public :: quadMehlerLin             ! Gauss-Mehler quadrature using linear transformation
!            =============             !   to (-1,1)

contains
!
!==============================================================================================
!
   real(REAL8) function quadLegendre( fun, n, umax )
!                       ============   p   i   i
!
!     Integrates the function fun from 0 to umax using Gauss-Legendre quadrature
!
!     The abscissas and weights are taken from the file gauss_legendre.dat. No attempt is made
!     to obtain abscissas for numbers of abscissas not contained in the file
!
!  Arguments:
!
      interface
         real(REAL8) function fun( u )          ! function to be integrated
            import :: REAL8
            real(REAL8), intent(in) :: u
         end function fun
      end interface
      integer(INT4), intent(in) :: n            ! number of abscissa
      real(REAL8), intent(in) :: umax           ! upper integration limit
!
!  Local Variables:
!
      logical :: found
      character(256) :: line, errmsg
      integer(INT4) :: i, istat, m
      integer(INT4), save :: n_save = 0, m_save
      real(REAL8) :: u_dummy, w_dummy
      real(REAL8), save :: umax_save = -1
      real(REAL8), save, allocatable :: u(:), w(:), u_save(:), w_save(:)
!$omp threadprivate( n_save, m_save, umax_save, u, w, u_save, w_save )
!
!  Start of Executable Code:
!  ------------------------
!
!  Calculate abscissas and weights
!
      if( n /= n_save ) then

         n_save = n

         ! search for number of abscissas in file

!$omp critical( critical_open_gl)               ! not sure why we need this critical section
         call openFile( 11, 'gauss_legendre.dat', 'formatted', 'read')
!             ========  i            i                 i         i

         found = .false.
         do
            read( 11, '(a)', iostat=istat )  line
            if( istat < 0 )  error stop 'quadLegendre: Too large number of abscissas.'
            if( line(1:1) == '#' )  cycle
            read( line, '(i3)', iostat=istat )  m
            if( m >= n ) then
               if( m > n )  print '(*(a,i4))', 'quadLegendre: Using n=', m, ' instead of ', n
               found = .true.
               exit
            end if
            do i = m/2+1, m
               read( 11, * )  u_dummy, w_dummy
            end do
         end do

         if( .not.found )  error stop 'quadLegendre: Invalid number of abscissa.'

         ! (re)allocate arrays

         if( allocated(u) )  deallocate( u, w )

         allocate( u( m ), stat=istat, errmsg=errmsg )
         call check_allocation( 'quadLegendre', 'u', istat, errmsg )

         allocate( w( m ), stat=istat, errmsg=errmsg )
         call check_allocation( 'quadLegendre', 'w', istat, errmsg )

         ! read abscissas and weights from file

         do i = m/2+1, m
            read( 11, * )  u(i), w(i)
         end do

         ! adjust abscissas and weights

         u(1:(m+1)/2) = - u(m:m/2+1:-1)           ! use symmetry
         w(1:(m+1)/2) = w(m:m/2+1:-1)

         close( 11 )

         m_save = m
         if( allocated( u_save ) )  deallocate( u_save, w_save )

         allocate( u_save( m_save ), stat=istat, errmsg=errmsg )
         call check_allocation( 'quadLegendre', 'u_save', istat, errmsg )

         allocate( w_save( m_save ), stat=istat, errmsg=errmsg )
         call check_allocation( 'quadLegendre', 'w_save', istat, errmsg )

         u_save = (u(:)+1)/2 * umax
         w_save = w(:)/2 * umax
!$omp end critical( critical_open_gl)

      end if

      if( umax /= umax_save ) then

         umax_save = umax

         ! transform (-1,1) -> (0,umax)

         u_save = (u(:)+1)/2 * umax
         w_save = w(:)/2 * umax

      end if
!
!  Calculate integral
!
      quadLegendre = ZERO

      do i = 1, m_save
         quadLegendre = quadLegendre + w_save(i) * fun( u_save(i) )
!                                                  ===      i
      end do

   end function quadLegendre
!
!==============================================================================================
!
   real(REAL8) function quadMehlerSymm( fun, n )
!                       ==============   p   i
!
!     Integrates the function fun(u)/sqrt(1-u^2) from 0 to 1 using Gauss-Mehler quadrature
!     (also called Chebyshev-Gauss quadrature)
!
!     Assumes fun(-u)=fun(u) for extending the integral from (0,1) to (-1,1) as required
!     by the Gauss-Mehler formula. Since Gauss-Mehler quadrature is exact if fun(u) is a
!     polynomial of degree < 2n, derivatives up to order 2n-1 should exist, and higher
!     derivatives should be small. This also applies to the mirror point u=0, therefore the odd
!     derivatives of fun(u) at u=0 should vanish or be small for the quadrature to be accurate.
!
!     This integration method has been proposed by F.J. Smith, Physics 30 (1964) 497.
!
!  Arguments:
!
      interface
         real(REAL8) function fun( u )          ! function to be integrated
            import :: REAL8
            real(REAL8), intent(in) :: u
         end function fun
      end interface
      integer(INT4), intent(in) :: n            ! number of abscissa
!
!  Local Variables:
!
      character(256) :: errmsg
      integer(INT4) :: i, istat
      integer(INT4), save :: n_save = 0
      real(REAL8), save, allocatable :: u(:), w(:)
!$omp threadprivate( n_save, u, w )
!
!  Start of Executable Code:
!  ------------------------
!
!  Calculate abscissas and weights
!
      if( n /= n_save ) then

         n_save = n

         ! (re)allocate arrays

         if( allocated(u) )  deallocate( u, w )

         allocate( u( n ), stat=istat, errmsg=errmsg )
         call check_allocation( 'quadMehlerSymm', 'u', istat, errmsg )

         allocate( w( n ), stat=istat, errmsg=errmsg )
         call check_allocation( 'quadMehlerSymm', 'w', istat, errmsg )

         ! Determine abscissas and weights

         do i = 1, n
            u(i) = cos( (TWO*i-ONE)/(FOUR*n) * PI )
         end do
         w(:) = PI / (TWO*n)

      end if
!
!  Calculate integral
!
      quadMehlerSymm = ZERO

      do i = 1, n
         quadMehlerSymm = quadMehlerSymm + w(i) * fun( u(i) )
!                                                 ===   i
      end do

   end function quadMehlerSymm
!
!==============================================================================================
!
   real(REAL8) function quadMehlerLin( fun, n )
!                       =============   p   i
!
!     Integrates the function fun(u)/sqrt(1-u^2) from 0 to 1 using Gauss-Mehler quadrature
!     (also called Chebyshev-Gauss quadrature)
!
!     Uses a linear transformation to map the range of integration from(0,1) to (-1,1) as
!     required by the Gauss-Mehler formula. Note that this will not exploit the full potential
!     of Gauss-Mehler quadrature since there is no pole at the lower integration boundary.
!
!     This integration method is presented in W. Eckstein, Computer Simulation of Ion-Solid
!     Interaction, Springer, 1991, p.20. The original reference might be Z. Kopal, Numerical
!     Analysis, 2nd ed, Wiley, 1961, p.381.
!
!  Arguments:
!
      interface
         real(REAL8) function fun( u )          ! function to be integrated
            import :: REAL8
            real(REAL8), intent(in) :: u
         end function fun
      end interface
      integer(INT4), intent(in) :: n            ! number of abscissa
!
!  Local Variables:
!
      character(256) :: errmsg
      integer(INT4) :: i, istat
      integer(INT4), save :: n_save = 0
      real(REAL8), save, allocatable :: u(:), w(:)
!$omp threadprivate( n_save, u, w )
!
!  Start of Executable Code:
!  ------------------------
!
!  Calculate abscissas and weights
!
      if( n /= n_save ) then

         n_save = n

         ! (re)allocate arrays

         if( allocated(u) )  deallocate( u, w )

         allocate( u( n ), stat=istat, errmsg=errmsg )
         call check_allocation( 'quadMehlerLin', 'u', istat, errmsg )

         allocate( w( n ), stat=istat, errmsg=errmsg )
         call check_allocation( 'quadMehlerLin', 'u', istat, errmsg )

         ! Determine abscissas and weights

         do i = 1, n
            u(i) = cos( (TWO*i-ONE)/(FOUR*n) * PI )
         end do
         w(:) = PI / n * u(:) / sqrt(ONE+u(:)**2)
         u(:) = u(:)**2

      end if
!
!  Calculate integral
!
      quadMehlerLin = ZERO

      do i = 1, n
         quadMehlerLin = quadMehlerLin + w(i) * fun( u(i) )
!                                               ===   i
      end do

   end function quadMehlerLin
!
!==============================================================================================
!
end module Math_Quadrature_mod
