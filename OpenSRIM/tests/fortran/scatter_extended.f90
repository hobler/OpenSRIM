! Black-box driver: each invocation initializes one native model/order/pair.
! Separate processes avoid changing the allocation lifecycle of native modules.
! Usage: scatter_extended MODEL Z1 Z2 MEHLER_ORDER
! Writes energy-major rows: reduced energy, reduced impact parameter, theta, tau.
! Theta uses native LegendreChiCorr(4); tau uses native cmangm at MEHLER_ORDER.
! conftest.py stages input tables and collects these files after successful runs.
program scatter_extended
   use Const_mod, only: REAL8
   use Param_Setup_mod, only: nelem
   use Param_Elem_mod, only: elem
   use Param_Sn_mod, only: sn
   use Undefine_mod
   use Sn_Screen_mod, only: initScreen, get_rmax
   use Sn_CMscatter_mod, only: cmAngle
   use Sn_Cmangm_mod, only: initCmangm, cmangm
   use Sn_Apsis_mod, only: apsis
   implicit none
   integer :: z1, z2, order, group, ie, ip, unit
   real(REAL8) :: e, p, theta, tau, unused, r0, rmax
   real(REAL8) :: energies(3), points(7)
   character(32) :: model, arg, label
   character(256) :: filename
   ! One invocation owns one initialized model/pair/order; initCmangm allocates
   ! module storage and is not a general reset operation for a new rule.
   call get_command_argument(1, model)
   call get_command_argument(2, arg); read(arg,*) z1
   call get_command_argument(3, arg); read(arg,*) z2
   call get_command_argument(4, arg); read(arg,*) order
   nelem = 2
   allocate(elem(2), sn%scrfun(2,2), sn%scrcoeffile(2,2), sn%scrfile(2,2), &
            sn%scrlen(2,2), sn%scatfile(2,2), sn%lquad(2,2))
   elem(1)%z=z1; elem(2)%z=z2
   sn%lquad=.false.
   call undefine(sn%scrfun)
   call undefine(sn%scrcoeffile)
   call undefine(sn%scrfile)
   call undefine(sn%scatfile)
   if(trim(model) == 'NLHlin') then
      sn%scrcoeffile='screencoefs/NLHLIN_python_coeffs.dat'
      sn%scrlen='NLHlin'
   else
      sn%scrfun='ZBLuniv'
      sn%scrlen='ZBL'
   end if
   call initScreen()
   sn%nmehl=order
   call initCmangm()
   rmax=get_rmax(1,2)
   ! Keep these input grids synchronized with scattering_support.extended_grid.
   ! Boundary points include the immediate representable predecessor of rmax.
   do group=1,3
      ! Higher orders are needed only by the dedicated convergence test.
      if(order /= 512 .and. group /= 3) cycle
      if(group==1) then
         label='pairs'
         energies=[1e-4_REAL8, 0.1_REAL8, 100._REAL8]
         points=[0.02_REAL8,0.1_REAL8,0.5_REAL8,1._REAL8,2._REAL8,3._REAL8,4._REAL8]
      else if(group==2) then
         label='boundary'
         energies=[0.001_REAL8,0.1_REAL8,10._REAL8]
         points=[1e-8_REAL8,1e-6_REAL8,1e-4_REAL8,0.01_REAL8,0.02_REAL8,0.03_REAL8,0.04_REAL8]
         if(trim(model)=='NLHlin') points(4:7) = &
            [rmax*(1-1e-4_REAL8),rmax*(1-1e-6_REAL8),rmax*(1-1e-8_REAL8),nearest(rmax,-1._REAL8)]
      else
         label='convergence'
         energies=[0.01_REAL8,1._REAL8,10._REAL8]
         points=[0.02_REAL8,0.1_REAL8,1._REAL8,3._REAL8,10._REAL8,20._REAL8,30._REAL8]
         if(trim(model)=='NLHlin') points = rmax * &
            [0.01_REAL8,0.1_REAL8,0.3_REAL8,0.5_REAL8,0.9_REAL8,0.99_REAL8,0.999_REAL8]
      end if
      write(filename,'(a,a,a,a,i0,a,i0,a,i0,a)') trim(label),'_',trim(model),'_',z1,'_',z2,'_n',order,'.txt'
      open(newunit=unit,file=trim(filename),status='replace')
      write(unit,'(a)') '# columns: e p theta_legendre4 tau_mehler'
      do ie=1,3
         e=energies(ie)
         do ip=1,7
            p=points(ip)
            call cmAngle(e,p,1,2,theta,method='LegendreChiCorr',nabsc=4)
            call apsis(e,p,1,2,r0)
            call cmangm(unused,tau,e,p,r0,1,2)
            write(unit,'(4es25.16)') e,p,theta,tau
         end do
      end do
      close(unit)
   end do
end program scatter_extended
