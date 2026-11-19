module Test_Proj_mod
!
!  Test the Ion_Proj_mod module
!
!  Used Modules:
!  ============
!
   use Const_mod,             only: INT4, REAL8,  &
                                    ZERO, PI
   use Types_mod,             only: Proj_type
   use String_mod,      only: str
   use Open_mod

   use Param_Setup_mod,       only: nelemIon, nr, nelem, ndim, nion, nion1, ldyn, usecell,  &
                                    rng, mrand, iarand, ibrand, irand, randskip
   use Param_Elem_mod,        only: elem     ! ?
   use Param_Geom_mod,        only: geom
   use Param_Ion_mod,         only: ion
   use param_mat_mod,         only: irxtal
   use Param_Sn_mod,          only: sn

   use Sn_FFPcyl_mod,         only: initFFPcyl
   use Random_mod,            only: initRannum
   use Ion_Count_mod,         only: initCounter
   use Ion_Proj_mod,          only: getProjectile
   use Target_Surface1d_mod,  only: updateSurface1d
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
   public :: initProj                           ! do initializations required for getProjectile
!            ========
   public :: testProj                           ! test the implementation of the projectile
!            ========                           !    initialization
   public :: testDiverg                         ! test the implementation of the projectile
!            ==========                         !    divergence
!
!=======================================================================
!
contains
!
!=======================================================================
!
   subroutine initProj( nbeam_ )
!             ========    i
!
!     Initialize data needed by getProjectile that are independent of the beam
!
!  Arguments:
!
      integer(INT4), intent(in), optional :: nbeam_
!
!  Local Variables:
!
      integer(INT4) :: nbeam
!
!  Start of Executable Code:
!  ------------------------
!
      if( present(nbeam_) ) then
         nbeam = nbeam_
      else
         nbeam = 1
      end if

      ion%usetra = .false.
      nion1 = 0
      nion = 1000000000
      nelemIon = 1
      allocate( ion%natom( 1 ) )
      ion%natom(1) = 1

      call initCounter( )
!          ===============

      ion%ltinit = .false.
      nr = 1
      nelem = 1
      allocate( elem( 1 ) )
      elem(1)%m = 1.
      elem(1)%z = 1.
      ion%m = 1.
      ion%nbeam = nbeam
      allocate( ion%weight( nbeam ), ion%diverg( 2, nbeam ), ion%moddiv( nbeam ) )
      ion%lrinit = .false.
      ion%lpinit = .false.
      ion%ranrot = .false.
      ion%xinit(:) = 0
      ion%yinit(:) = 0
      ion%zinit(:) = 0
      irxtal = 0
      ndim = 1
      sn%pmax = 3.
      sn%ffpmax = 10.
      geom%lcart = .true.
      ion%liinit = .false.
      allocate( geom%posif( 0:1 ) )
      geom%posif(0) = 0.
      geom%posif(1) = 10000.
      call updateSurface1d( geom%posif(0), geom%posif(1) )
!          =============        i              i
      ldyn = .false.
      usecell = .false.

      rng = 'Haas'
      iarand = 35678
      ibrand = 197638
      irand = 5421
      mrand = 7103
      randskip = 0
      call initRannum( )
!          ==========

      call initFFPcyl( )
!          ==========

   end subroutine initProj
!
!=======================================================================
!
   subroutine testProj( tilty, rotate, diverg, moddiv, fwhm )
!             ========   i      i       i       i      i
!
!     Test the implementation of the projectile initialization
!
!  Arguments:
!
      real(REAL8), intent(in) :: tilty       ! tilt angle about y axis (deg)
      real(REAL8), intent(in) :: rotate      ! rotation angle (deg)
      real(REAL8), intent(in) :: diverg(:)   ! beam divergence (deg)
      integer(INT4), intent(in) :: moddiv    ! beam divergence mode
      real(REAL8), intent(in) :: fwhm(2)     ! full width at half maximum (A)
!
!  Local Variables:
!
      integer(INT4), parameter :: n = 100000, nh = 30

      logical :: newmolecule, found, done
      integer(INT4) :: i, ix, iy
      real(REAL8) :: x, y, xmin, xmax, dx, ymin, ymax, sigmamax
      real(REAL8) :: dirx, diry, dirmin, dirmax, ddir
      real(REAL8) :: hist_pos(nh,nh), hist_dir(nh,nh)
      real(REAL8) :: xmean, ymean, dirxmean, dirymean
      real(REAL8) :: xstd, ystd, dirxstd, dirystd
      type(Proj_type) :: proj
!
!  Start of Executable Code:
!  ------------------------
!
!  Put input parameters into right place
!
      ion%rotate = rotate
      ion%tiltx = 0
      ion%tilty = tilty
      ion%diverg(:,1) = diverg(:)
      ion%moddiv(1) = moddiv
      ion%sigma(:) = fwhm(:) / sqrt( 8 * log(2.) )
!
!  Initialize arrays for position and direction statistics
!
      sigmamax = max( 1e-10_REAL8, ion%sigma(1), ion%sigma(2) )
      xmax = -3*tan(tilty*PI/180)*cos(rotate*PI/180) +  3*sigmamax
      xmin = -3*tan(tilty*PI/180)*cos(rotate*PI/180) -  3*sigmamax
      dx = (xmax-xmin) / nh
      ymax = -3*tan(tilty*PI/180)*sin(rotate*PI/180) +  3*sigmamax
      ymin = -3*tan(tilty*PI/180)*sin(rotate*PI/180) -  3*sigmamax
      hist_pos(:,:) = 0

      xmean = 0
      ymean = 0
      xstd = 0
      ystd = 0

      dirmax = 1
      dirmin = -dirmax
      ddir = (dirmax-dirmin) / nh
      hist_dir(:,:) = 0

      dirxmean = 0
      dirymean = 0
      dirxstd = 0
      dirystd = 0
!
!  Generate projectiles and collect statistics
!
      do i = 1, n

         call getProjectile( 1, .true., proj, found, done )
!             =============  i    i      o      o     o
         !
         !  Position histogram
         !
         x = proj%pos(1)
         ix = ceiling( (x-xmin) / dx )
         y = proj%pos(2)
         iy = ceiling( (y-ymin) / dx )

         if( ix >= 1 .and. ix <= nh .and. iy >= 1 .and. iy <= nh ) then
            hist_pos(ix,iy) = hist_pos(ix,iy) + 1
         end if
         !
         ! Position moments
         !
         xmean = xmean + x
         ymean = ymean + y
         xstd = xstd + x**2
         ystd = ystd + y**2
         !
         !  Angle histograms
         !
         dirx = proj%dir(1)
         ix = ceiling( (dirx-dirmin) / ddir )
         diry = proj%dir(2)
         iy = ceiling( (diry-dirmin) / ddir )

         if( ix >= 1 .and. ix <= nh .and. iy >= 1 .and. iy <= nh ) then
            hist_dir(ix,iy) = hist_dir(ix,iy) + 1
         end if
         !
         ! Angle moments
         !
         dirxmean = dirxmean + dirx
         dirymean = dirymean + diry
         dirxstd = dirxstd + dirx**2
         dirystd = dirystd + diry**2

      end do

      if( any( fwhm(:) > 0 ) ) then
         xmean = xmean / n
         ymean = ymean / n
         xstd = sqrt( xstd/n - xmean**2 )
         ystd = sqrt( ystd/n - ymean**2 )
         print '(a,2f8.4)', 'sigma=', ion%sigma(1), ion%sigma(2)
         print '(2(a,f8.4))', '   x=', xmean, ' +/- ', xstd
         print '(2(a,f8.4))', '   y=', ymean, ' +/- ', ystd

         call plotProj( xmin, xmax, ymin, ymax, hist_pos )
!             ========   i     i     i     i       i
      end if

      if( any( diverg(:) > 0 ) ) then
         dirxmean = dirxmean / n
         dirymean = dirymean / n
         dirxstd = sqrt( dirxstd/n - dirxmean**2 )
         dirystd = sqrt( dirystd/n - dirymean**2 )
         print '(a,2f8.4)',  'diverg=', diverg(:)*PI/180
         print '(2(a,f8.4))',  '   dirx=', dirxmean, ' +/- ', dirxstd
         print '(2(a,f8.4))',  '   diry=', dirymean, ' +/- ', dirystd

         call plotProj( dirmin, dirmax, dirmin, dirmax, hist_dir )
!             ========    i       i       i       i        i
      end if

   end subroutine testProj
!
!=======================================================================
!
   subroutine plotProj( xmin, xmax, ymin, ymax, hist )
!             ========   i     i     i     i     i
!
!     Plot the 2D histograms using a system call to a Python
!     script which is written by this subroutine. Needs some plots to compare with.
!
!  Arguments:
!
      real(REAL8), intent(in) :: xmin, xmax, ymin, ymax
      real(REAL8), intent(in) :: hist(:,:)
!
!  Local Variables:
!
      integer(INT4) :: i, nx, j, ny
      real(REAL8) :: dx, dy
!
!  Start of Executable Code:
!  ------------------------
!
      nx = size(hist,1)
      dx = (xmax-xmin) / nx
      ny = size(hist,2)
      dy = (ymax-ymin) / ny

      call openFile( 12, 'x.dat', 'formatted', 'write' )
!          ========  i      i          i          i
      write( 12, '(100f10.4)' )  (xmin+i*dx, i=0,nx)
      close( 12 )

      call openFile( 12, 'y.dat', 'formatted', 'write' )
!          ========  i      i          i          i
      write( 12, '(100f10.4)' )  (ymin+i*dy, i=0,ny)
      close( 12 )

      call openFile( 12, 'hist.dat', 'formatted', 'write' )
!          ========  i          i             i          i
      do i = 1, nx
         write( 12, '(100f10.3)' )  (hist(j,i), j=1,ny)
      end do
      close( 12 )
!
!  Write the Python script
!
      call openFile( 12, 'plot_proj.py', 'formatted', 'write' )
!          ========  i          i               i          i
      write( 12, '(a)' )  'import numpy as np'
      write( 12, '(a)' )  'import matplotlib.pyplot as plt'
      write( 12, '(a)' )  'x = np.loadtxt("x.dat")'
      write( 12, '(a)' )  'y = np.loadtxt("y.dat")'
      write( 12, '(a)' )  'hist = np.loadtxt("hist.dat")'
      write( 12, '(a)' )  'plt.pcolor(x, y, hist)'
      write( 12, '(a)' )  'plt.gca().set_aspect("equal")'
      write( 12, '(a)' )  'plt.show()'
      write( 12, '(a)' )  ''
      close( 12 )
!
!  Execute the Python script
!
      call execute_command_line('python3 plot_proj.py')
!
!  Clean up
!
      !call execute_command_line('rm plot_proj.py x.dat hist.dat')
!
   end subroutine plotProj
!
!=======================================================================
!
   subroutine testDiverg( weight, diverg, moddiv)
!             ==========    i       i       i
!
!     Test the implementation of the beam divergence initialization
!
!  Arguments:
!
      real(REAL8), intent(in) :: weight(:)   ! beam weight
      real(REAL8), intent(in) :: diverg(:,:) ! beam divergence (deg)
      integer(INT4), intent(in) :: moddiv    ! beam divergence mode
!
!  Local Variables:
!
      integer(INT4), parameter :: n = 100000, nh = 30

      logical :: newmolecule, found, done
      integer(INT4) :: i, ix, iy
      real(REAL8) :: dirx, diry, dirmin, dirmax, ddir
      real(REAL8) :: hist_dir(nh,nh), hist1(nh)
      real(REAL8) :: dirxmean, dirymean
      real(REAL8) :: dirxstd, dirystd
      type(Proj_type) :: proj
!
!  Start of Executable Code:
!  ------------------------
!
!  Put input parameters into right place
!
      ion%rotate = 0
      ion%tiltx = 0
      ion%tilty = 0
      ion%weight(:) = weight(:)
      ion%diverg(:,:) = diverg(:,:)
      print *, 'ion%diverg_x=', ion%diverg(1,:)*PI/180
      print *, 'ion%diverg_y=', ion%diverg(2,:)*PI/180
      ion%moddiv(:) = moddiv
      ion%sigma(:) = 0
!
!  Initialize arrays for position and direction statistics
!
      dirmax = maxval( 3*diverg(:,:)*PI/180 )
      print *, 'dirmax=', dirmax
      dirmin = -dirmax
      ddir = (dirmax-dirmin) / nh
      hist_dir(:,:) = 0

      dirxmean = 0
      dirymean = 0
      dirxstd = 0
      dirystd = 0
!
!  Generate projectiles and collect statistics
!
      do i = 1, n

         call getProjectile( 1, .true., proj, found, done )
!             =============  i    i      o      o     o
         !
         !  Angle histograms
         !
         dirx = proj%dir(1)
         ix = ceiling( (dirx-dirmin) / ddir )
         diry = proj%dir(2)
         iy = ceiling( (diry-dirmin) / ddir )

         if( ix >= 1 .and. ix <= nh .and. iy >= 1 .and. iy <= nh ) then
            hist_dir(ix,iy) = hist_dir(ix,iy) + 1
         end if
         !
         ! Angle moments
         !
         dirxmean = dirxmean + dirx
         dirymean = dirymean + diry
         dirxstd = dirxstd + dirx**2
         dirystd = dirystd + diry**2

      end do

      hist_dir(:,:) = hist_dir(:,:) / n / ((dirmax-dirmin)/nh)

      dirxmean = dirxmean / n
      dirymean = dirymean / n
      dirxstd = sqrt( dirxstd/n - dirxmean**2 )
      dirystd = sqrt( dirystd/n - dirymean**2 )
      print '(a,2f8.4)', 'diverg_x=', diverg(1,:)*PI/180
      print '(a,2f8.4)', 'diverg_y=', diverg(2,:)*PI/180
      print '(2(a,f8.4))', '   dirx=', dirxmean, ' +/- ', dirxstd
      print '(2(a,f8.4))', '   diry=', dirymean, ' +/- ', dirystd
!
!  Plot distribution in x direction
!
      hist1(:) = sum( hist_dir(:,:), dim=2 )
      print *, 'norm of histogram=', sum(hist1)*(dirmax-dirmin)/nh
      call plotDiverg( dirmin, dirmax, hist1(:), moddiv, weight(:), diverg(1,:), 'x' )
!          ==========    i       i        i         i       i           i         i
!
!  Plot distribution in y direction
!
      hist1(:) = sum( hist_dir(:,:), dim=1 )
      call plotDiverg( dirmin, dirmax, hist1(:), moddiv, weight(:), diverg(2,:), 'y' )
!          ==========    i       i        i        i        i           i         i

   end subroutine testDiverg
!
!=======================================================================
!
   subroutine plotDiverg( xmin, xmax, hist, moddiv, weight, diverg, xy )
!             ==========   i     i     i      i       i       i
!
!     Plot the 1D histograms using a system call to a Python script which is written by this
!     subroutine.
!
!     Compares with the probability density function.
!
!  Arguments:
!
      real(REAL8), intent(in) :: xmin, xmax
      real(REAL8), intent(in) :: hist(:)
      integer(INT4) :: moddiv, ibeam, nbeam
      real(REAL8) :: weight(:)
      real(REAL8) :: diverg(:)
      character(1) :: xy
!
!  Local Variables:
!
      integer(INT4) :: i, nx
      real(REAL8) :: dx, histmax, histmin
!
!  Start of Executable Code:
!  ------------------------
!
      nx = size(hist)
      dx = (xmax-xmin) / nx
      nbeam = size(weight)

      call openFile( 12, 'x.dat', 'formatted', 'write' )
!          ========  i      i          i          i
      write( 12, '(100f10.4)' )  ((xmin+i*dx)*180/PI, i=0,nx)
      close( 12 )

      call openFile( 12, 'hist.dat', 'formatted', 'write' )
!          ========  i          i             i          i
      do i = 1, nx
         write( 12, '(f10.3)' )  hist(i)*PI/180
      end do
      close( 12 )

      histmax = 10 ** ceiling( log10( maxval( hist(:)*PI/180 ) ) )
      histmin = histmax / 1e4
!
!  Write the Python script
!
      call openFile( 12, 'plot_diverg.py', 'formatted', 'write' )
!          ========  i          i               i          i
      write( 12, '(a)' )  'import numpy as np'
      write( 12, '(a)' )  'import matplotlib.pyplot as plt'
      write( 12, '(a)' )  'x = np.loadtxt("x.dat")'
      write( 12, '(a)' )  'hist = np.loadtxt("hist.dat")'
      write( 12, '(a)' )  'plt.bar(x[:-1], hist, align="edge", width=' //  &
                             str(dx*180/PI, 'f99.4') // ')'

      write( 12, '(a)' )  'x = np.linspace(x[0], x[-1], 200)'
      write( 12, '(a)' )  'pdf = np.zeros(200)'
      write( 12, '(a)' )  'weights = []'
      write( 12, '(a)' )  'divergs = []'
      do ibeam = 1, nbeam
         write( 12, '(a)' )  'weight = ' // str( weight(ibeam) )
         write( 12, '(a)' )  'weights.append(weight)'
         write( 12, '(a)' )  'diverg = ' // str( diverg(ibeam) )
         write( 12, '(a)' )  'divergs.append(diverg)'
         if( moddiv == 1 ) then
            write( 12, '(a)' )  'pdf += weight*2/(np.pi*diverg) ' //  &
                                      '* np.sqrt(1 - np.minimum(1,(x/diverg)**2))'
         else
            write( 12, '(a)' )  'pdf += weight/(np.sqrt(2*np.pi)*diverg)' //  &
                                      ' * np.exp(-0.5*(x/diverg)**2)'
         end if
         write( 12, '(a)' )  'plt.plot(x, pdf, "r--")'
      end do
      write( 12, '(a)' )  'plt.plot(x, pdf, "r-")'
      write( 12, '(a)' )  'plt.yscale("log")'
      write( 12, '(a)' )  'plt.ylim(' // str( histmin ) // ',' // str( histmax ) // ')'
      write( 12, '(a)' )  'plt.xlabel(r"180/$\pi\cdot$dir$_\mathrm{' // xy //'} "' //  &
                                     'r"\approx \psi_\mathrm{' // xy // '}$ [deg]")'
      write( 12, '(a)' )  'plt.ylabel("frequency")'
      write( 12, '(a)' )  'title = "weights="+", ".join([f"{int(100*x):d}%" for x in weights])'
      write( 12, '(a)' )  'title += ", diverg="'
      write( 12, '(a)' )  'title += ", ".join([f"{x:.2f}" + r"$\degree$" for x in divergs])'
      write( 12, '(a)' )  'plt.title(title)'
      write( 12, '(a)' )  'plt.show()'
      write( 12, '(a)' )
      close( 12 )
!
!  Execute the Python script
!
      call execute_command_line('python3 plot_diverg.py')
!
!  Clean up
!
      !call execute_command_line('rm plot_proj.py x.dat hist.dat')

   end subroutine plotDiverg
!
!=======================================================================
!
end module Test_Proj_mod
!
!#######################################################################
!
program test

   use const_mod, only: INT4, REAL8
   use Test_Proj_mod

   integer(INT4) :: choice, moddiv, ibeam, nbeam
   real(REAL8) :: tilty, rotate, sigma(2), fwhm(2)
   real(REAL8), allocatable :: weight(:), diverg(:,:)

   print *, ' '
   print *, 'Choose among the following options:'
   print *, ' '
   print *, '   1 ... Beam profile'
   print *, '   2 ... Beam divergence (basic)'
   print *, '   3 ... Beam divergence (advanced)'
   print *, ' '
   read *, choice

   select case( choice )
      case( 1 )
         print *, 'Enter tilt (deg), rotate (deg)'
         read *, tilty, rotate
         allocate( diverg( 2, 1 ) )
         nbeam = 1
         diverg(:,1) = 0
         moddiv = 1
         print *, 'Enter sigma_x (A), sigma_y (A)'
         read *, sigma(1), sigma(2)
         fwhm(:) = sigma(:) * sqrt( 8 * log(2.) )
      case( 2 )
         print *, 'Enter tilt (deg), rotate (deg)'
         read( *, * )  tilty, rotate
         nbeam = 1
         allocate( diverg( 2, 1 ) )
         print *, 'Enter diverg_x (deg), diverg_y (deg), moddiv'
         read *, diverg(1,1), diverg(2,1), moddiv
         fwhm = 0
      case( 3 )
         print *, 'Enter number of beams'
         read *, nbeam
         allocate( weight( nbeam ) )
         allocate( diverg( 2, nbeam ) )
         do ibeam = 1,nbeam
            print *, 'Enter weight, diverg_x (deg), diverg_y (deg)'
            read *, weight(ibeam), diverg(:,ibeam)
         end do
         weight(:) = weight(:) / sum( weight(:) )
         print *, 'Enter moddiv'
         read *, moddiv
   end select

   call initProj( nbeam )
!       ========    i
   select case( choice )
      case( 1, 2 )
         call testProj( tilty, rotate, diverg(:,1), moddiv, fwhm(:) )
!             ========   i      i       i       i      i
      case( 3 )
         print *, 'weight=', weight
         call testDiverg( weight(:), diverg(:,:), moddiv )
!             ==========     i           i          i
   end select

end program test
