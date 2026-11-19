module Test_Screen_mod
!
!  Test the Sn_Screen_mod module
!
!  Used Modules:
!  ============
!
   use Const_mod,       only: INT4, REAL8,  &
                              ZERO
   use Open_mod

   use Param_Setup_mod, only: nelem
   use Param_Elem_mod,  only: elem
   use Param_Sn_mod,    only: sn

   use Sn_Screen_mod,   only: initScreen, screen
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
   public :: testScreen                         ! test the implementation of the screeining
!            ==========                         !    functions
!
!=======================================================================
!
contains
!
!=======================================================================
!
   subroutine testScreen( scrlen )
!             ==========    i
!
!     Test the implementation of the screening functions
!
!     Plots the screening functions as a function of R using a system call to a Python
!     script which is written by this subroutine. Needs some plots to compare with.
!
!  Arguments:
!
      character(*), intent(in) :: scrlen
!
!  Local Variables:
!
      integer(INT4) :: i
      integer(INT4), parameter :: nmax = 1000
      real(REAL8) :: rmax = 30.
      real(REAL8) :: r(0:nmax), phi, phi1
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
      sn%scrlen(:,:) = scrlen
!
!  Radii to be plotted
!
      do i = 0, nmax
         r(i) = i * rmax/nmax
      end do
!
!  Calculate screening functions and write them to files
!
      print *, 'ZBLuniv'
      sn%scrfun = 'ZBL'
      call initScreen( )
!          ==========
      call openFile( 12, 'ZBLuniv.dat', 'formatted', 'write' )
!          ========  i         i             i          i
      do i = 0, nmax
         call screen( r(i), 1, 2, phi, phi1 )
!             ======   i    i  i   o    o
         write( 12, * ) r(i), phi
      end do
      close( 12 )

      print *, 'ZBLspec'
      sn%scrfun(2,1) = 'coeffile'
      sn%coeffile(2,1) = '/home/hobler/work/imsil/make_tables/tables/screencoefs/ZBLspec'
      call initScreen( )
!          ==========
      call openFile( 12, 'ZBL_Ga_W.dat', 'formatted', 'write' )
!          ========  i          i             i          i
      do i = 0, nmax
         call screen( r(i), 1, 2, phi, phi1 )
!             ======   i    i  i   o    o
         write( 12, * ) r(i), phi
      end do
      close( 12 )

      print *, 'Kr-C'
      sn%scrfun(2,1) = 'KrC'
      call initScreen( )
!          ==========
      call openFile( 12, 'Kr-C.dat', 'formatted', 'write' )
!          ========  i       i            i          i
      do i = 0, nmax
         call screen( r(i), 1, 2, phi, phi1 )
!             ======   i    i  i   o    o
         write( 12, * ) r(i), phi
      end do
      close( 12 )

      print *, 'Moliere'
      sn%scrfun(2,1) = 'Moliere'
      call initScreen( )
!          ==========
      call openFile( 12, 'Moliere.dat', 'formatted', 'write' )
!          ========  i         i             i          i
      do i = 0, nmax
         call screen( r(i), 1, 2, phi, phi1 )
!             ======   i    i  i   o    o
         write( 12, * ) r(i), phi
      end do
      close( 12 )

      print *, 'GaW'
      sn%scrfun(2,1) = 'scrfile'
      sn%scrfile(2,1) = '/home/hobler/git/zblpot/zblscreen_0.001/GaW'
      call initScreen( )
!          ==========
      call openFile( 12, 'ZBLtab_Ga_W.dat', 'formatted', 'write' )
!          ========  i          i             i          i
      do i = 0, nmax
         call screen( r(i), 1, 2, phi, phi1 )
!             ======   i    i  i   o    o
         write( 12, * ) r(i), phi
      end do
      close( 12 )

!
!  Write the Python script
!
      call openFile( 12, 'plot_screen.py', 'formatted', 'write' )
!          ========  i          i               i          i
      write( 12, '(a)' )  'import numpy as np'
      write( 12, '(a)' )  'import matplotlib.pyplot as plt'
      write( 12, '(a)' )  'r, phi = np.loadtxt("ZBLuniv.dat", unpack=True)'
      write( 12, '(a)' )  'plt.semilogy(r, phi, label="ZBL universal")'
      write( 12, '(a)' )  'r, phi = np.loadtxt("ZBL_Ga_W.dat", unpack=True)'
      write( 12, '(a)' )  'plt.semilogy(r, phi, label="ZBL (Ga-W)")'
      write( 12, '(a)' )  'r, phi = np.loadtxt("Kr-C.dat", unpack=True)'
      write( 12, '(a)' )  'plt.semilogy(r, phi, label="Kr-C")'
      write( 12, '(a)' )  'r, phi = np.loadtxt("Moliere.dat", unpack=True)'
      write( 12, '(a)' )  'plt.semilogy(r, phi, label="Moliere")'
      write( 12, '(a)' )  'r, phi = np.loadtxt("ZBLtab_Ga_W.dat", unpack=True)'
      write( 12, '(a)' )  'plt.semilogy(r, phi, label="ZBLtab (Ga-W)")'
      write( 12, '(a)' )  'plt.xlim(0, 30)'
      write( 12, '(a)' )  'plt.ylim(1e-4, 1)'
      write( 12, '(a)' )  'plt.legend()'
      write( 12, '(a)' )  'plt.show()'
      close( 12 )
!
!  Execute the Python script
!
      call execute_command_line('python plot_screen.py')
!
!  Clean up
!
      call execute_command_line(  &
            'rm plot_screen.py ZBLuniv.dat ZBL_Ga_W.dat Kr-C.dat Moliere.dat')
!
   end subroutine testScreen
!
!=======================================================================
!
end module Test_Screen_mod
!
!#######################################################################
!
program test

   use const_mod, only: INT4
   use Test_Screen_mod

   integer(INT4) :: choice

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
         call testScreen( '1' )
      case( 2 )
         call testScreen( 'ZBL' )
      case( 3 )
         call testScreen( 'Firsov' )
      case( 4 )
         call testScreen( 'Lindhard' )
   end select

end program test
