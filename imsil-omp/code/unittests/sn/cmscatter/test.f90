!
!##############################################################################################
!
! Legacy native unit-test module used by the Fortran make targets.
! The pytest reference fixtures build the dedicated test_cm_scatter_* drivers
! instead; this module is not their data-generation entry point.
module Test_CMscatter_mod
!
!  Test the CMscatter_mod module
!
!  Used Modules:
!  ============
!
   use Const_mod,          only: INT4, REAL8,  &
                                 ZERO, HALF, PI
   use Open_mod
!   use Str_mod
!   use Apsis_mod,          only: apsis
!   use Table_mod,          only: Table2d_type, initTable, interpolate
   use Undefine_mod
!
   use Param_Setup_mod,    only: nelem
   use Param_Elem_mod,           only: elem
   use Param_Sn_mod,       only: sn
   use Sn_CMscatter_mod
   use Sn_Screen_mod,      only: initScreen, screen
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
   public :: initTestCMscatter                  ! common initializations for all tests
!            =================
   public :: testConvergence                    ! test convergence of Gauss-Mehler quadrature
!            ===============
!   public :: testDerivatives                    ! test the derivatives
!!            ===============
!   public :: testHobler                         ! test new Gauss-Mehler implementation
!!            ==========
   public :: plotThetaErr                       ! plot errors in theta as a function of p
!            ============
   public :: plotThetaErr3                      ! plot errors in theta as a function of p
!            =============                      !   (3 energies)
!   public :: testGauss                          ! test Gauss-Legendre quadrature
!!            =========
!
!=======================================================================
!
contains
!
!=======================================================================
!
   subroutine initTestCMscatter( scrfun, scrlen )
!             =================    i       i
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
      allocate( sn%scrcoeffile( nelem, nelem ) )
      allocate( sn%scrfile( nelem, nelem ) )
      allocate( sn%scrlen( nelem, nelem ) )
      allocate( sn%scatfile( nelem, nelem ) )
      call undefine( sn%scrfun )
      call undefine( sn%scrcoeffile )
      call undefine( sn%scrfile )
      call undefine( sn%scrlen )
      call undefine( sn%scatfile )
!
!  Initialize screening function
!
      if( scrfun == 'scrcoeffile' ) then
         sn%scrcoeffile(:,:) = '/home/hobler/work/imsil/master-omp/tables/screencoefs/ZBLspec_coeffs.dat'
      else if( scrfun == 'scrfile' ) then
         sn%scrfile(:,:) = '/home/hobler/work/imsil/master-omp/tables/screentab/ZBLspec/screening_zblspec_31_74.dat'
      else
         sn%scrfun(:,:) = scrfun
      end if
!
!  Initialize screening length
!
      sn%scrlen(:,:) = scrlen

      call initScreen( )
!          ==========

   end subroutine initTestCMscatter
!
!=======================================================================
!
   subroutine testConvergence( method )
!
!     Test the convergence of the quadrature
!
!  Arguments:
!
      character(*), intent(in) :: method
!
!  Local Variables:
!
      integer(INT4), parameter :: kmax = 6   ! k=3 seems to be optimum (larger values
                                             !    cause loss of significant digits)
                                             ! New Gauss-Legendre table goes to 64 only
      integer(INT4) :: k, nabsc
      real(REAL8) :: e, p, theta_ref(0:kmax), theta(0:kmax), tau(0:kmax)
!
!  Start of Executable Code:
!  ------------------------
!
      do
!
!  Ask for input parameters
!
         print *, 'Enter eps, P'
         read *, e, p

         do k = 0, kmax
            nabsc = 2**k
            print *, 'LegendreChiCorr, n=', nabsc
            call cmAngle( e, p, 1, 2, theta_ref(k), method='LegendreChiCorr', nabsc=nabsc )
!                =======  i  i  i  i       o                i                 i
            theta_ref(k) = theta_ref(k) * 180/PI
            print *, trim(method) // ', n=', nabsc
            call cmAngle( e, p, 1, 2, theta(k), method=method, nabsc=nabsc )
!                =======  i  i  i  i     o          i             i
            theta(k) = theta(k) * 180/PI
         end do

         print '(/a)', 'LegendreChiCorr, theta (deg):'
         do k = 0, kmax
            print '(3(g15.7,a))', theta_ref(k), ' (err = ', theta_ref(k)-theta_ref(kmax),  &
                                  ' = ', 100*(theta_ref(k)/theta_ref(kmax)-1), ' %)'
         end do

         print '(/a)', trim(method) // ', theta (deg):'
         do k = 0, kmax
            print '(3(g15.7,a))', theta(k), ' (err = ', theta(k)-theta(kmax), ' = ',  &
                                  100*(theta(k)/theta(kmax)-1), ' %)'
         end do

!         print '(/a)', 'tau:'
!         do k = 0, kmax
!            print '(3(g15.7,a))', tau(k), ' (err = ', tau(k)-tau(kmax), ' = ',  &
!                                  100*(tau(k)/tau(kmax)-1), ' %)'
!         end do

      end do

   end subroutine testConvergence
!
!==============================================================================================
!
   subroutine plotThetaErr( method )
!                             i
!     Plots the relative and absolute error in theta as a function of impact parameter
!     for a few energies
!
!  Arguments:
!
      character(*), intent(in) :: method
!
!  Local Variables:
!
      integer(INT4), parameter :: NABSC_REF = 1000
      integer(INT4), parameter :: NNABSC = 6
      integer(INT4) :: nabsc(NNABSC) = [ 1, 2, 3, 4, 5, 10 ]
      integer(INT4), parameter :: NE = 2, NP = 160
      character(4) :: etext(NE) = [ '1e-4', '1e-1' ]
      real(REAL8) :: pfac(NE) = [ 0.25, 0.25 ]

      character(80) :: fname
      integer(INT4) :: ie, ip, iabsc, z1, z2
      real(REAL8) :: e, p, theta(NNABSC), theta_ref, dummy

      character(2), dimension(92) :: atom = [ &
         'H ', 'He', 'Li', 'Be', 'B ', 'C ', 'N ', 'O ', 'F ', 'Ne',  &
         'Na', 'Mg', 'Al', 'Si', 'P ', 'S ', 'Cl', 'Ar', 'K ', 'Ca',  &
         'Sc', 'Ti', 'V ', 'Cr', 'Mn', 'Fe', 'Co', 'Ni', 'Cu', 'Zn',  &
         'Ga', 'Ge', 'As', 'Se', 'Br', 'Kr', 'Rb', 'Sr', 'Y ', 'Zr',  &
         'Nb', 'Mo', 'Tc', 'Ru', 'Rh', 'Pd', 'Ag', 'Cd', 'In', 'Sn',  &
         'Sb', 'Te', 'I ', 'Xe', 'Cs', 'Ba', 'La', 'Ce', 'Pr', 'Nd',  &
         'Pm', 'Sm', 'Eu', 'Gd', 'Tb', 'Dy', 'Ho', 'Er', 'Tm', 'Yb',  &
         'Lu', 'Hf', 'Ta', 'W ', 'Re', 'Os', 'Ir', 'Pt', 'Au', 'Hg',  &
         'Tl', 'Pb', 'Bi', 'Po', 'At', 'Rn', 'Fr', 'Ra', 'Ac', 'Th',  &
         'Pa', 'U ' ]
!
!  Start of Executable Code:
!  ------------------------
!
!  Hobler
!
      ! Loop over energies
      do ie = 1, NE

         read( etext(ie), * )  e

         call openFile( 12, 'theta_err_' // trim(etext(ie)) // '.dat'  , 'formatted', 'write' )
!             ========  i                     i                         i          i
         write( 12, '(2i5)' )  NNABSC+1, NP+1
         write( 12, '(*(i7))' ) NABSC_REF, nabsc(:)
!
         ! Loop over impact parameters
         do ip = 0, NP

            p = pfac(ie) * ip

            ! reference value
            call cmAngle( e, p, 1, 2, theta_ref, method='MehlerSymm', nabsc=NABSC_REF )
!                =======  i  i  i  i      o           i               i
            theta_ref = theta_ref * 180/PI

            ! loop over number of abscissas
            do iabsc = 1, NNABSC
               call cmAngle( e, p, 1, 2, theta(iabsc), method=method, nabsc=nabsc(iabsc) )
!                   =======  i  i  i  i       o              i             i
               theta(iabsc) = theta(iabsc) * 180/PI
            end do

            ! write theta and errors to the file
            write( 12, '(2g13.6, *(es13.5))' ) p, theta_ref, theta(:)-theta_ref

         end do

      end do

      close( 12 )
!
!  Eckstein
!
      ! Loop over energies
      do ie = 1, NE

         read( etext(ie), * )  e

         call openFile( 12, 'theta_err_eck_' // trim(etext(ie)) // '.dat'  , 'formatted', 'write' )
!             ========  i                     i                         i          i
         write( 12, '(2i5)' )  NNABSC+1, NP+1
         write( 12, '(*(i7))' ) NABSC_REF, nabsc(:)
!
         ! Loop over impact parameters
         do ip = 0, NP

            p = pfac(ie) * ip

            ! reference value
            call cmAngle( e, p, 1, 2, theta_ref, method='MehlerLinInacc', nabsc=NABSC_REF )
!                =======  i  i  i  i      o            i                      i
            theta_ref = theta_ref * 180/PI

            ! loop over number of abscissas
            do iabsc = 1, NNABSC
               call cmAngle( e, p, 1, 2, theta(iabsc), method='MehlerLinInacc', nabsc=nabsc(iabsc) )
!                   =======  i  i  i  i       o                i                       i
               theta(iabsc) = theta(iabsc) * 180/PI
            end do

            ! write theta and errors to the file
            write( 12, '(2g13.6, *(es13.5))' ) p, theta_ref, theta(:)-theta_ref

         end do

      end do

      close( 12 )
!
!  Magic formula
!
      ! Loop over energies
      do ie = 1, NE

         read( etext(ie), * )  e

         call openFile( 12, 'theta_magic_' // trim(etext(ie)) // '.dat'  , 'formatted', 'write' )
!             ========  i                         i                             i          i
         write( 12, '(2i5)' )  2, NP+1
         write( 12, '(*(i7))' ) 0, 0
!
         ! Loop over impact parameters
         do ip = 0, NP

            p = pfac(ie) * ip

            ! reference value
            call cmang_magic( e, p, 1, 2, theta_ref )
!                ===========  i  i  i  i      o
            theta_ref = theta_ref * 180/PI

            ! write theta to the file
            write( 12, '(2g13.6)' ) p, theta_ref

         end do

      end do

      close( 12 )
!
!  Write and run python script to display errors
!
      call openFile( 12, 'plot_theta_err.py', 'formatted', 'write' )
!          ========  i          i               i          i
      write( 12, '(a)' )  'import numpy'
      write( 12, '(a)' )  'import matplotlib as mpl'
      write( 12, '(a)' )  'import matplotlib.pyplot as plt'
      write( 12, '(a)' )  ' '
      write( 12, '(a)' )  'e = ("' // trim(etext(1)) // '", "' // trim(etext(2)) // '")'
      write( 12, '(a)' )  ' '
      write( 12, '(a)' )  'fig = plt.figure(figsize=(7,7))'
      write( 12, '(a)' )  ' '
      write( 12, '(a)' )  'for ie in range(len(e)): '
      write( 12, '(a)' )  '  # Hobler'
      write( 12, '(a)' )  '  with open("theta_err_" + e[ie] + ".dat") as f:'
      write( 12, '(a)' )  '    nc, np = [int(n) for n in f.readline().split()]'
      write( 12, '(a)' )  '    l = [int(n) for n in f.readline().split()]'
      write( 12, '(a)' )  '    nabsc_ref = l[0]'
      write( 12, '(a)' )  '    nabsc = l[1:]'
      write( 12, '(a)' )  '    p = []'
      write( 12, '(a)' )  '    theta_ref = []'
      write( 12, '(a)' )  '    theta_err = []'
      write( 12, '(a)' )  '    for ic in range(nc-1):'
      write( 12, '(a)' )  '      theta_err.append([])'
      write( 12, '(a)' )  '    for ip in range(np):'
      write( 12, '(a)' )  '      l = [float(r) for r in f.readline().split()]'
      write( 12, '(a)' )  '      pp = l[0]'
      write( 12, '(a)' )  '      th_ref = l[1]'
      write( 12, '(a)' )  '      th_err = l[2:]'
      write( 12, '(a)' )  '      p.append(pp)'
      write( 12, '(a)' )  '      theta_ref.append(th_ref)'
      write( 12, '(a)' )  '      for ic in range(nc-1):'
      write( 12, '(a)' )  '        theta_err[ic].append(th_err[ic])'
      write( 12, '(a)' )  '  p = numpy.array(p)'
      write( 12, '(a)' )  '  theta_ref = numpy.array(theta_ref)'
      write( 12, '(a)' )  '  theta_err = numpy.array(theta_err)'
      write( 12, '(a)' )  ' '
      write( 12, '(a)' )  '  # Eckstein'
      write( 12, '(a)' )  '  with open("theta_err_eck_" + e[ie] + ".dat") as f:'
      write( 12, '(a)' )  '    _, np = [int(n) for n in f.readline().split()]'
      write( 12, '(a)' )  '    l = [int(n) for n in f.readline().split()]'
      write( 12, '(a)' )  '    nmehl_eck = l[-1]'
      write( 12, '(a)' )  '    theta_err_eck = []'
      write( 12, '(a)' )  '    for ip in range(np):'
      write( 12, '(a)' )  '      l = [float(r) for r in f.readline().split()]'
      write( 12, '(a)' )  '      th_err_eck = l[-1]'
      write( 12, '(a)' )  '      theta_err_eck.append(th_err_eck)'
      write( 12, '(a)' )  '  theta_err_eck = numpy.array(theta_err_eck)'
      write( 12, '(a)' )  ' '
      write( 12, '(a)' )  '  # Magic'
      write( 12, '(a)' )  '  with open("theta_magic_" + e[ie] + ".dat") as f:'
      write( 12, '(a)' )  '    _, np = [int(n) for n in f.readline().split()]'
      write( 12, '(a)' )  '    f.readline()'
      write( 12, '(a)' )  '    theta_magic = []'
      write( 12, '(a)' )  '    for ip in range(np):'
      write( 12, '(a)' )  '      l = [float(r) for r in f.readline().split()]'
      write( 12, '(a)' )  '      theta_magic.append(l[-1])'
      write( 12, '(a)' )  '  theta_magic = numpy.array(theta_magic)'
      write( 12, '(a)' )  ' '
      write( 12, '(a)' )  '  ax = fig.add_subplot(2, 2, 1+2*ie)'
      write( 12, '(a)' )  '  for ic in range(nc-1):'
      write( 12, '(a,i2,a)' )  '    plt.semilogy(p, numpy.abs(theta_err[ic,:])' //  &
                               ', label="n=" + str(nabsc[ic]))'
      write( 12, '(a)' )  '  plt.semilogy(p, numpy.abs(theta_err_eck), "k--"' //  &
                          ', label="n=" + str(nabsc[ic]) + "\n(Eckstein)")'
      write( 12, '(a)' )  '  plt.semilogy(p, numpy.abs(theta_magic-theta_ref), "k:"' //  &
                          ', label="Magic")'
      write( 12, '(a)' )  '  plt.ylim(1e-6,1)'
      write( 12, '(a)' )  '  plt.xlabel(r"impact parameter $P$")'
      write( 12, '(a)' )  '  plt.ylabel(r"absolute error in $\theta$ (deg)")'
      write( 12, '(a)' )  '  plt.legend(loc="lower right")'
      write( 12, '(a)' )  '  plt.title(r"$\varepsilon=$" + e[ie] + ", ' //  &
                          trim(sn%scrfun(2,1)) // ', ' // method // '")'
      write( 12, '(a)' )  ' '
      write( 12, '(a)' )  '  ax = fig.add_subplot(2, 2, 2+2*ie)'
      write( 12, '(a)' )  '  for ic in range(nc-1):'
      write( 12, '(a,i2,a)' )  '    plt.semilogy(p, numpy.abs(theta_err[ic,:])/theta_ref[:]' //  &
                               ', label="n=" + str(nabsc[ic]))'
      write( 12, '(a)' )  '  plt.semilogy(p, numpy.abs(theta_err_eck/theta_ref), "k--"' //  &
                          ', label="n=" + str(nabsc[ic]) + "\n(Eckstein)")'
      write( 12, '(a)' )  '  plt.semilogy(p, numpy.abs(theta_magic/theta_ref-1), "k:"' //  &
                          ', label="Magic")'
      write( 12, '(a)' )  '  plt.ylim(1e-6,1)'
      write( 12, '(a)' )  '  plt.xlabel(r"impact parameter $P$")'
      write( 12, '(a)' )  '  plt.ylabel(r"relative error in $\theta$")'
      write( 12, '(a)' )  '  plt.legend(loc="lower right")'
      write( 12, '(a)' )  '  plt.title(r"$\varepsilon=$" + e[ie] + ", ' //  &
                          trim(sn%scrfun(2,1)) // ', ' // method // '")'
      write( 12, '(a)' )  ' '
      write( 12, '(a)' )  'plt.tight_layout()'
      write( 12, '(a)' )  'plt.show()'
      close( 12 )
      call execute_command_line('python plot_theta_err.py')

   end subroutine plotThetaErr
!
!==============================================================================================
!
   subroutine plotThetaErr3()
!
!     Plots the relative error in theta as a function of impact parameter for 3 energies
!
!  Arguments:
!
!     None
!
!  Local Variables:
!
      character(*), parameter :: potential = 'ZBLuniv'
      integer(INT4), parameter :: NM = 3
      character(15) :: method(NM) = [ 'LegendreChiCorr', 'LegendreChi', 'LegendreGCorr' ]
      character(15) :: method_ref = 'LegendreChiCorr'
      character(2) :: linestyle(NM) = [ '-', '--', ':' ]
      integer(INT4), parameter :: NABSC_REF = 100
      integer(INT4), parameter :: NNABSC = 6
      integer(INT4) :: nabsc(NNABSC) = [ 1, 2, 3, 4, 5, 10 ]
      integer(INT4), parameter :: NE = 3, NP = 40
      character(4) :: etext(NE) = [ '1e-4', '0.1', '100' ]
      real(REAL8) :: pfac(NE) = [ 1., 0.3, 0.1 ]

      integer(INT4) :: ie, ip, im, iabsc
      real(REAL8) :: e, p, theta(NNABSC), theta_ref, dummy
!
!  Start of Executable Code:
!  ------------------------
!
!  Hobler
!
      ! Loop over energies
      do ie = 1, NE

         read( etext(ie), * )  e

         do im = 1, NM

            call openFile( 12, 'theta_err_' // trim(etext(ie)) // trim(linestyle(im)) //  &
!                ========  i                                   i
                           '.dat',  'formatted', 'write' )
!                                        i          i
            write( 12, '(2i5)' )  NNABSC+1, NP+1
            write( 12, '(*(i7))' ) NABSC_REF, nabsc(:)

            ! Loop over impact parameters
            do ip = 0, NP

               p = pfac(ie) * (ip + 0.01)

               ! reference value
               call cmAngle( e, p, 1, 2, theta_ref, method=trim(method_ref), nabsc=NABSC_REF )
   !                =======  i  i  i  i      o               i                    i
               theta_ref = theta_ref * 180/PI

               ! loop over number of abscissas
               do iabsc = 1, NNABSC
                  call cmAngle( e, p, 1, 2, theta(iabsc), method=method(im), nabsc=nabsc(iabsc) )
   !                   =======  i  i  i  i       o              i             i
                  theta(iabsc) = theta(iabsc) * 180/PI
               end do

               ! write theta and errors to the file
               write( 12, '(2g13.6, *(es13.5))' ) p, theta_ref, theta(:)-theta_ref

            end do

         end do

      end do

      close( 12 )
!
!  Write and run python script to display errors
!
      call openFile( 12, 'plot_theta_err3.py', 'formatted', 'write' )
!          ========  i          i               i          i
      write( 12, '(a)' )  'import numpy'
      write( 12, '(a)' )  'import matplotlib as mpl'
      write( 12, '(a)' )  'import matplotlib.pyplot as plt'
      write( 12, '(a)' )  ' '
      write( 12, '(a)' )  'colors = plt.rcParams["axes.prop_cycle"].by_key()["color"]'
      write( 12, '(a)' )  'iplot = [1, 3, 4]'
      write( 12, '(*(a))', advance='no' )  'e = (', ('"', trim(etext(ie)), '", ', ie=1,ne)
      write( 12, '(a)' )   ')'
      write( 12, '(*(a))', advance='no' )  'linestyles = (',  &
                           ('"', trim(linestyle(im)), '", ', im=1,NM)
      write( 12, '(a)' )   ')'
      write( 12, '(*(a))', advance='no' )  'methods = (',  &
                           ('"', trim(method(im)), '", ', im=1,NM)
      write( 12, '(a)' )   ')'
      write( 12, '(a)' )  ' '
      write( 12, '(a)' )  'fig = plt.figure(figsize=(7,7))'
      write( 12, '(a)' )  ' '
      write( 12, '(a)' )  'for ie in range(len(e)): '
      write( 12, '(a,i2,a)' )  ' for im in range(', NM, '):'
      write( 12, '(a)' )  '  with open("theta_err_" + e[ie] + linestyles[im] + ".dat") as f:'
      write( 12, '(a)' )  '    nc, np = [int(n) for n in f.readline().split()]'
      write( 12, '(a)' )  '    l = [int(n) for n in f.readline().split()]'
      write( 12, '(a)' )  '    nabsc_ref = l[0]'
      write( 12, '(a)' )  '    nabsc = l[1:]'
      write( 12, '(a)' )  '    p = []'
      write( 12, '(a)' )  '    theta_ref = []'
      write( 12, '(a)' )  '    theta_err = []'
      write( 12, '(a)' )  '    for ic in range(nc-1):'
      write( 12, '(a)' )  '      theta_err.append([])'
      write( 12, '(a)' )  '    for ip in range(np):'
      write( 12, '(a)' )  '      l = [float(r) for r in f.readline().split()]'
      write( 12, '(a)' )  '      pp = l[0]'
      write( 12, '(a)' )  '      th_ref = l[1]'
      write( 12, '(a)' )  '      th_err = l[2:]'
      write( 12, '(a)' )  '      p.append(pp)'
      write( 12, '(a)' )  '      theta_ref.append(th_ref)'
      write( 12, '(a)' )  '      for ic in range(nc-1):'
      write( 12, '(a)' )  '        theta_err[ic].append(th_err[ic])'
      write( 12, '(a)' )  '  p = numpy.array(p)'
      write( 12, '(a)' )  '  theta_ref = numpy.array(theta_ref)'
      write( 12, '(a)' )  '  theta_err = numpy.array(theta_err)'
      write( 12, '(a)' )  ' '
      write( 12, '(a)' )  '  ax = fig.add_subplot(2, 2, iplot[ie])'
      write( 12, '(a)' )  '  for ic in range(nc-1):'
      write( 12, '(a)' )  '    plt.semilogy(p, numpy.abs(theta_err[ic,:])/theta_ref[:], '
      write( 12, '(a)' )  '                 linestyle=linestyles[im], color=colors[ic])'
      write( 12, '(a)' )  '  plt.xlim(p[0],p[-1])'
      write( 12, '(a)' )  '  plt.ylim(1e-7,1)'
      write( 12, '(a)' )  '  plt.xlabel(r"impact parameter $P$")'
      write( 12, '(a)' )  '  plt.ylabel(r"relative error in $\theta$")'
      write( 12, '(a)' )  '  plt.title(r"$\varepsilon=$" + e[ie] )'
      write( 12, '(a)' )  ' '
      write( 12, '(a)' )  'ax = fig.add_subplot(2, 2, 2)'
      write( 12, '(a)' )  'plt.axis("off")'
      write( 12, '(a)' )  'lines = []'
      write( 12, '(a)' )  'for ic in range(nc-1):'
      write( 12, '(a)' )  '  line, = plt.plot([], [], color=colors[ic], ' //  &
                                              'label="n=" + str(nabsc[ic]))'
      write( 12, '(a)' )  '  lines.append(line)'
      write( 12, '(a)' )  'legend1 = plt.legend(handles=lines, loc=(0.03, 0.35))'
      write( 12, '(a)' )  'plt.gca().add_artist(legend1)'
      write( 12, '(a)' )  'lines = []'
      write( 12, '(a)' )  'for im in (1,0,2):'
      write( 12, '(a)' )  '  line, = plt.plot([], [], linestyle=linestyles[im], color="k", '
      write( 12, '(a)' )  '                   label=methods[im])'
      write( 12, '(a)' )  '  lines.append(line)'
      write( 12, '(a)' )  'legend2 = plt.legend(handles=lines, loc="lower left")'
      write( 12, '(a)' )  'plt.gca().add_artist(legend2)'
      write( 12, '(a)' )  ' '
      write( 12, '(a)' )  'plt.xlim(0,1)'
      write( 12, '(a)' )  'plt.ylim(0,1)'
      write( 12, '(*(a))' )  'plt.text(0.03, 0.9, "potential: ', potential, '", ',  &
                                       'fontsize="large")'
      write( 12, '(a)' )  ' '
      write( 12, '(a)' )  'plt.tight_layout()'
      write( 12, '(a)' )  'plt.show()'
      close( 12 )
      call execute_command_line('python plot_theta_err3.py')



   end subroutine plotThetaErr3
!
!==============================================================================================
!
end module Test_CMscatter_mod
!
!##############################################################################################
!
program main

   use const_mod, only: INT4
   use Test_CMscatter_mod

   character(8) :: scrlen
   character(11) :: scrfun
   character(15) :: method
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
         scrfun = 'scrcoeffile'
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

   call initTestCMscatter( scrfun, scrlen )
!       =================    i       i

   print *, ' '
   print *, 'Choose an integration method:'
   print *, ' '
   print *, '   1 ... LegendreChi'
   print *, '   2 ... LegendreChiCorr'
   print *, '   3 ... LegendreGCorr'
   print *, '   4 ... LegendreInacc'
   print *, '   5 ... Legendre'
   print *, '   6 ... LegendreAcc'
   print *, '   7 ... MehlerSymmInacc'
   print *, '   8 ... MehlerLinInacc'
   print *, '   9 ... MehlerSymm'
   print *, '  10 ... MehlerLin'
   print *, '  11 ... MehlerSymmAcc'
   print *, '  12 ... MehlerLinAcc'
   read *, choice

   select case( choice )
      case( 1 )
         method = 'LegendreChi'
      case( 2 )
         method = 'LegendreChiCorr'
      case( 3 )
         method = 'LegendreGCorr'
      case( 4 )
         method = 'LegendreInacc'
      case( 5 )
         method = 'Legendre'
      case( 6 )
         method = 'LegendreAcc'
      case( 7 )
         method = 'MehlerSymmInacc'
      case( 8 )
         method = 'MehlerLinInacc'
      case( 9 )
         method = 'MehlerSymm'
      case( 10 )
         method = 'MehlerLin'
      case( 11 )
         method = 'MehlerSymmAcc'
      case( 12 )
         method = 'MehlerLinAcc'
      case default
         print *, 'Invalid answer.'
   end select

   print *, ' '
   print *, 'Choose a test:'
   print *, ' '
   print *, '   1 ... test convergence'
   print *, '   2 ... plot theta error'
   print *, '   3 ... plot theta error 3'
   print *, ' '
   read *, choice

   select case( choice )
      case( 1 )
         call testConvergence( method )
      case( 2 )
         call plotThetaErr( method )
      case( 3 )
         call plotThetaErr3()
      case default
         print *, 'Invalid answer.'
   end select

end program main
