! Generate production scattering references for one model, element pair and order.
!
! Args (command line, in this order):
!     MODEL (character): 'ZBL' or 'NLHlin'; supplied by conftest.extended_fortran.
!     Z1, Z2 (integer): Projectile and target atomic numbers.
!     LEGENDRE_ORDER (integer): Number of quadrature points for theta and tau.
! Inputs:
!     Screening coefficients and quadrature tables staged by conftest.py.
! Outputs:
!     Six TXT files: two impact grids times three energies, one energy per file.
!     Columns: reduced energy e, reduced impact parameter p, theta, tau.
! Notes:
!     A fresh process initializes each model/pair/order combination, keeping
!     the allocation lifecycle of the production modules isolated between runs.
program generate_extended_scattering_references
   use Const_mod, only: REAL8
   use Param_Setup_mod, only: nelem
   use Param_Elem_mod, only: elem
   use Param_Sn_mod, only: sn
   use Undefine_mod
   use Sn_Screen_mod, only: initScreen, get_rmax
   use Sn_CMscatter_mod, only: cmScatterLegendre
   implicit none
   integer :: z1, z2, order, group, ie, ip, unit
   real(REAL8) :: e, p, theta, tau, rmax
   real(REAL8) :: energies(3), points(7)
   character(32) :: model, arg, label
   character(256) :: filename
   character(6) :: model_label
   character(5) :: energy_labels(3)
   ! Read the model as text and parse the remaining command-line strings as integers.
   ! The pytest fixture supplies supported models, element pairs and orders.
   call get_command_argument(1, model)
   call get_command_argument(2, arg); read(arg,*) z1
   call get_command_argument(3, arg); read(arg,*) z2
   call get_command_argument(4, arg); read(arg,*) order
   ! Allocate two element slots and the pair-indexed screening configuration.
   ! Atomic numbers are stored in elem; scattering calls refer to slot indices 1,2.
   nelem = 2
   allocate(elem(2), sn%scrfun(2,2), sn%scrcoeffile(2,2), sn%scrfile(2,2), &
            sn%scrlen(2,2), sn%scatfile(2,2), sn%lquad(2,2))
   elem(1)%z=z1; elem(2)%z=z2
   ! Initialize flags and clear optional selectors before configuring the model.
   ! The explicit cmScatterLegendre call below controls the quadrature order.
   sn%lquad=.false.
   call undefine(sn%scrfun)
   call undefine(sn%scrcoeffile)
   call undefine(sn%scrfile)
   call undefine(sn%scatfile)
   ! Choose the coefficient table or analytic screening function and its length
   ! normalization. model_label is the lowercase form used in filenames.
   if(trim(model) == 'NLHlin') then
      model_label='nlhlin'
      sn%scrcoeffile='screencoefs/NLHLIN_python_coeffs.dat'
      sn%scrlen='NLHlin'
   else
      model_label='zbl'
      sn%scrfun='ZBLuniv'
      sn%scrlen='ZBL'
   end if
   ! Initialize screening, then obtain the reduced cutoff radius for boundary grids.
   call initScreen()
   rmax=get_rmax(1,2)
   ! Keep these input grids synchronized with scattering_support.extended_grid.
   ! Boundary points include the immediate representable predecessor of rmax.
   do group=1,2
      ! The pairs grid samples seven positive impact parameters at three energies.
      if(group==1) then
         label='pairs'
         energy_labels=['1e-04','1e-01','1e+02']
         energies=[1e-4_REAL8, 0.1_REAL8, 100._REAL8]
         points=[0.02_REAL8,0.1_REAL8,0.5_REAL8,1._REAL8,2._REAL8,3._REAL8,4._REAL8]
      ! The boundary grid starts with three very small positive impact parameters.
      ! For ZBL, the remaining points also stay near zero.
      else if(group==2) then
         label='boundary'
         energy_labels=['1e-03','1e-01','1e+01']
         energies=[0.001_REAL8,0.1_REAL8,10._REAL8]
         points=[1e-8_REAL8,1e-6_REAL8,1e-4_REAL8,0.01_REAL8,0.02_REAL8,0.03_REAL8,0.04_REAL8]
         ! For finite-range NLHlin, place the final four points just inside the cutoff,
         ! ending at the largest representable value smaller than rmax.
         if(trim(model)=='NLHlin') points(4:7) = &
            [rmax*(1-1e-4_REAL8),rmax*(1-1e-6_REAL8),rmax*(1-1e-8_REAL8),nearest(rmax,-1._REAL8)]
      end if
      ! Write each energy separately, using the same name fields as Python references.
      do ie=1,3
         e=energies(ie)
         write(filename,'(a,a,a,i0,a,i0,a,a,a,i0,a,a,a)') &
            'scatter_',trim(model_label),'_z',z1,'_z',z2,'_norm_native_e',energy_labels(ie), &
            '_n',order,'_grid_',trim(label),'.txt'
         ! Create this case in the temporary working directory; the header documents
         ! the four numeric columns and is ignored by numpy.loadtxt.
         open(newunit=unit,file=trim(filename),status='replace')
         write(unit,'(a)') '# columns: e p theta_legendre_chi_corr tau_legendre'
         do ip=1,7
            p=points(ip)
            ! Compute theta and tau at the requested order using the production solver.
            ! Scientific notation retains the precision required by cross-backend tests.
            call cmScatterLegendre(e,p,1,2,theta,tau,nabsc=order)
            write(unit,'(4es25.16)') e,p,theta,tau
         end do
         close(unit)
      end do
   end do
end program generate_extended_scattering_references
