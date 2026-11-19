module Test_Apsis_mod
!
!  Tests the Apsis_mod module
!
!  Used Modules:
!  ============
!
   use Const_mod,          only: INT4, REAL8,  &
                                 ZERO, TEN
   use Open_mod

   use Param_Setup_mod,    only: nelem
   use Param_Elem_mod,     only: elem
   use Param_Sn_mod,       only: sn

   use Sn_Apsis_mod
   use Sn_Screen_mod,      only: initScreen
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
   public :: initTestApsis                      ! common initializations for all tests
!            =============
   public :: testApsis                          ! test the calculation of the apsis
!            =========
   public :: testApsisDer                       ! test the calculation of the apsis derivatives
!            ============
   public :: plotIterations                     ! plots the number of iteration as a function
!            ==============                     !    of energy and impact parameter
!
!=======================================================================
!
contains
!
!=======================================================================
!
   subroutine initTestApsis( scrfun, scrlen )
!             =============    i       i
!
!     Do initializations common to all tests implemented in this module
!
!  Arguments:
!
      character(*), intent(in) :: scrfun  ! screening function
      character(*), intent(in) :: scrlen  ! screening length
!
!  Local Variables:
!
!     None
!
!  Start of Executable Code:
!  ------------------------
!
!  Prepare parameters needed by initScreen
!
      nelem = 2
      allocate( elem( nelem ) )
      elem(1)%z = 31
      elem(2)%z = 74
      allocate( sn%scrfun( nelem, nelem ) )
      allocate( sn%coeffile( nelem, nelem ) )
      allocate( sn%scrfile( nelem, nelem ) )
      allocate( sn%scrlen( nelem, nelem ) )
!
!  Initialize screening function
!
      sn%scrfun(:,:) = scrfun
      if( sn%scrfun(2,1) == 'coeffile' ) then
         sn%coeffile(:,:) = '/home/hobler/work/imsil/make_tables/tables/screencoefs/ZBLspec'
      else if( sn%scrfun(2,1) == 'scrfile' ) then
         sn%scrfile(:,:) = '/home/hobler/git/zblpot/zblscreen_0.001/GaW'
      end if
!
!  Initialize screening length
!
      sn%scrlen(:,:) = scrlen

      call initScreen( )
!          ==========

   end subroutine initTestApsis
!
!=======================================================================
!
   subroutine testApsis( )
!             =========
!
!     Calculate the apsis of collisions for given reduced energy and impact parameter.
!     The input is requested from standard input. The result has to be compared to
!     external data (e.g. Fig.2.8.a of [Eck91]).
!
!  Arguments:
!
!     None
!
!  Local Variables:
!
      integer(INT4) :: n
      real(REAL8) :: e, p, r0, drde, drdp
!
!  Start of Executable Code:
!  ------------------------
!
      do
!
!  Request e and p
!
         print *, 'Enter eps, P'
         read *, e, p
!
!  Calculate and print apsis
!
         call apsis( e, p, 1, 2, r0, drde, drdp, n )
!             =====  i  i  i  i  o    o     o    o
         print *, 'r0=', r0, ', n=', n
         call apsis( e, p, 1, 2, r0, drde, drdp, n )  ! test saving values
!             =====  i  i  i  i  o    o     o    o
         print *, 'r0=', r0, ', n=', n
         print *, ' '

      end do

   end subroutine testApsis
!
!=======================================================================
!
   subroutine testApsisDer()
!             ============
!
!     Calculate the apsis of collisions for two slightly different reduced energies,
!     calculate difference quotient, and compare with difference
!     The input is requested from standard input.
!
!  Arguments:
!
!     None
!
!  Local Variables:
!
      real(REAL8) :: e, p, r0, dr0de, dr0dp, deps, dp, r0e, r0p, dummy
!
!  Start of Executable Code:
!  ------------------------
!
      do
!
!  Request e and p
!
         print *, 'Enter eps, P'
         read *, e, p
!
!  Calculate apsis and derivatives, and difference quotients
!
         call apsis( e, p, 1, 2, r0, dr0de, dr0dp )
!             =====  i  i  i  i  o     o      o
         call apsis( e, p, 1, 2, dummy, dr0de, dr0dp ) ! test saving values
!             =====  i  i  i  i    o      o      o

         deps = 0.001 * e
         call apsis( e+deps, p, 1, 2, r0e )
!             =====     i    i  i  i   o
         print *, 'd(r0)/d(eps)=', dr0de, (r0e-r0)/deps

         if( p == ZERO ) then
            dp = 1e-6
         else
            dp = 0.001 * p
         end if
         call apsis( e, p+dp, 1, 2, r0p )
!             =====  i   i    i  i   o
         print *, 'd(r0)/d(P)=  ', dr0dp, (r0p-r0)/dp

      end do

   end subroutine testApsisDer
!
!=======================================================================
!
   subroutine plotIterations( )
!             ==============
!
!     Plot the number of iterations necessary to calculate the apsis using various
!     methods. Results are displayed by a colourmesh with energy and impact parameter
!     on the axes.
!
!  Arguments:
!
!     None
!
!  Local Variables:
!
      integer(INT4), parameter :: ne = 9
      integer(INT4), parameter :: np = 11
      integer(INT4) :: i, n(ne,np), ie, ip
      real(REAL8) :: e(ne) = [ (TEN**i, i=-6,2) ]
      real(REAL8) :: p(np) =  [ (0.4*i**2, i=0,10) ]
      real(REAL8) :: r0
!
!  Start of Executable Code:
!  ------------------------
!
!  Newton with r0=p or r0(e) as initial condition
!
      do ip = 1, np
         do ie = 1, ne
            call apsis( e(ie), p(ip), 1, 2, r0, nn=n(ie,ip) )
!                =====    i      i    i  i  o       o
            call apsis( e(ie), p(ip), 1, 2, r0, nn=n(ie,ip) )  ! second call tests saving of
!                =====    i      i    i  i  o       o          !    values
         end do
      end do

      call openFile( 12, 'iter_pe_n.dat', 'formatted', 'write' )
!          ========  i          i              i          i
      write( 12, '(2i5)' )  ne, np
      write( 12, '(*(es11.3))' )  (e(ie), ie=1, ne)
      write( 12, '(*(es11.3))' )  (p(ip), ip=1, np)
      do ip = 1, np
         write( 12, '(*(i3))' )  (n(ie,ip), ie=1,ne)
      end do
      close( 12 )
!
!  Write and run python script to display errors
!
      call openFile( 12, 'plot_iter.py', 'formatted', 'write' )
!          ========  i          i             i          i
      write( 12, '(a)' )  'import numpy'
      write( 12, '(a)' )  'import matplotlib as mpl'
      write( 12, '(a)' )  'import matplotlib.pyplot as plt'
      write( 12, '(a)' )  ' '
      write( 12, '(a)' )  'text1 = ( r"$R_{0,0}=P$",'
      write( 12, '(a)' )  '          r"$R_{0,0}=\max(P,R_{0,0,\varepsilon})$" )'
      write( 12, '(a)' )  'text2 = ( "Newton",'
      write( 12, '(a)' )  '          r"$\phi(R)$+Newton",'
      write( 12, '(a)' )  '          r"$\phi(R)$",'
      write( 12, '(a)' )  '#          "LogNewton",'
      write( 12, '(a)' )  '          r"$\phi(R)$+LogNewton")'
      write( 12, '(a)' )  ' '
      write( 12, '(a)' )  'bounds = numpy.linspace(-0.5, 10.5, 12)'
      write( 12, '(a)' )  'ticks = numpy.linspace(0, 10, 11)'
      write( 12, '(a)' )  'cmap = plt.cm.jet'
      write( 12, '(a)' )  'norm = mpl.colors.BoundaryNorm(bounds, cmap.N)'
      write( 12, '(a)' )  ' '
      write( 12, '(a)' )  'with open("iter_pe_n.dat") as f:'
      write( 12, '(a)' )  '  ne, np = [int(n) for n in f.readline().split()]'
      write( 12, '(a)' )  '  e = [float(e) for e in f.readline().split()]'
      write( 12, '(a)' )  '  p = [float(p) for p in f.readline().split()]'
      write( 12, '(a)' )  ' '
      write( 12, '(a)' )  '  e_index = numpy.arange(len(e)+1) - 0.5'
      write( 12, '(a)' )  '  p_index = numpy.arange(len(p)+1) - 0.5'
      write( 12, '(a)' )  '  n = []'
      write( 12, '(a)' )  '  for ip in range(len(p)):'
      write( 12, '(a)' )  '    n.append([float(rstr) for rstr in f.readline().split()])'
      write( 12, '(a)' )  ' '
      write( 12, '(a)' )  'fig = plt.figure(figsize=(7,6))'
      write( 12, '(a)' )  'ax = plt.gca()'
      write( 12, '(a)' )  ' '
      write( 12, '(a)' )  'c = ax.pcolormesh(e_index, p_index, n, '
      write( 12, '(a)' )  '                  vmin=0, vmax=10, cmap=cmap, norm=norm)'
      write( 12, '(a)' )  ' '
      write( 12, '(a)' )  'e_ticks = e_index[:-1] + 0.5'
      write( 12, '(a)' )  'ax.set_xticks(e_ticks)'
      write( 12, '(a)' )  'ax.set_xticklabels(numpy.asarray(numpy.log10(e), dtype=int))'
      write( 12, '(a)' )  'ax.set_xlabel(r"$^{10}$log(energy $\varepsilon$)")'
      write( 12, '(a)' )  'p_ticks = p_index[:-1] + 0.5'
      write( 12, '(a)' )  'ax.set_yticks(p_ticks)'
      write( 12, '(a)' )  'ax.set_yticklabels(p)'
      write( 12, '(a)' )  'ax.set_ylabel(r"impact parameter $P$")'
      write( 12, '(a)' )  'cbar = fig.colorbar(c, ax=ax, ticks=ticks)'
      write( 12, '(a)' )  'cbar.set_label(r"number of iterations $n$")'
      write( 12, '(a)' )  'ax.set_title(r"$R_{0,0}=\max(P,R_{0,0,\varepsilon})$" +'
      write( 12, '(a)' )  '             ", Newton", fontsize="small")'
      write( 12, '(a)' )  ' '
      write( 12, '(a)' )  'plt.tight_layout()'
      write( 12, '(a)' )  'plt.show()'
      close( 12 )
      call execute_command_line('python plot_iter.py')

   end subroutine plotIterations
!
!=======================================================================
!
end module Test_Apsis_mod
!
!#######################################################################
!
program test

   use const_mod, only: INT4
   use Test_Apsis_mod

   character(8) :: scrlen, scrfun
   integer(INT4) :: choice

   print *, ' '
   print *, 'Choose among the following screening functions:'
   print *, ' '
   print *, '   1 ... universal ZBL'
   print *, '   2 ... KrC'
   print *, '   3 ... Moliere'
   print *, '   4 ... specific ZBL (GaW screening function parameters)'
   print *, '   5 ... specific ZBL (GaW screening function table)'
   print *, ' '
   read *, choice

   select case( choice )
      case( 1 )
         scrfun = 'ZBL'
      case( 2 )
         scrfun = 'KrC'
      case( 3 )
         scrfun = 'Moliere'
      case( 4 )
         scrfun = 'coeffile'
      case( 5 )
         scrfun = 'scrfile'
   end select

   print *, ' '
   print *, 'Choose among the following screening lengths:'
   print *, ' '
   print *, '   1 ... 1 A'
   print *, '   2 ... ZBL'
   print *, '   3 ... Firsov'
   print *, '   4 ... Lindhard'
   print *, ' '
   read *, choice

   select case( choice )
      case( 1 )
         scrlen = '1'
      case( 2 )
         scrlen = 'ZBL'
      case( 3 )
         scrlen = 'Firsov'
      case( 4 )
         scrlen = 'Lindhard'
   end select

   call initTestApsis( scrfun, scrlen )
!       =============    i       i

   do
      print *, ' '
      print *, 'Choose among the following tests:'
      print *, ' '
      print *, '   1 ... Calculate apsis for single (e,p) values'
      print *, '   2 ... Test apsis derivatives for single (e,p) values'
      print *, '   3 ... Plot iteration count for array of (e,p) values'
      print *, ' '
      read *, choice

      select case( choice )
         case( 1 )
            call testApsis()
         case( 2 )
            call testApsisDer()
         case( 3 )
            call plotIterations()
         case default
            print *, 'Invalid answer.'
      end select
   end do

end program test

