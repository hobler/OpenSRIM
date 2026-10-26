! Generate native screening references for test_nlhlin.py.
! Two models x three configurations x 92 pairs produce 552 files of 101 rows.
! Each row stores the test coordinate, phi and dphi/d(test coordinate).
! conftest.py provides the canonical coefficient table in a temporary cwd.
program screen_reference
   use Const_mod, only: REAL8
   use Param_Setup_mod, only: nelem
   use Param_Elem_mod, only: elem
   use Param_Sn_mod, only: sn
   use Undefine_mod
   use Sn_Screen_mod, only: initScreen, screen
   implicit none
   integer :: model, config, z1, z2, i, unit
   real(REAL8) :: rnorm, native_norm, coordinate, radius, phi, dphi, p1, p2
   character(256) :: filename
   character(6) :: label
   nelem = 2
   allocate(elem(2), sn%scrfun(2,2), sn%scrcoeffile(2,2), sn%scrfile(2,2), &
            sn%scrlen(2,2), sn%scatfile(2,2), sn%lquad(2,2))
   sn%lquad = .false.
   do model = 1, 2
      if(model == 1) label = 'nlhlin'
      if(model == 2) label = 'zbl'
      ! Configurations match test_nlhlin.CONFIGURATIONS. The first two use
      ! equal-element pairs; the third fixes the second element to copper.
      do config = 1, 3
         if(config == 1) then
            p1 = 0; p2 = 0
         else if(config == 2) then
            p1 = 0.23_REAL8; p2 = 1
         else
            p1 = 0.5_REAL8; p2 = 2._REAL8/3._REAL8
         end if
         do z1 = 1, 92
            z2 = z1
            if(config == 3) z2 = 29
            elem(1)%z = z1; elem(2)%z = z2
            if(config == 1) then
               rnorm = 1
            else
               rnorm = 0.4685_REAL8/(real(z1,REAL8)**p1+real(z2,REAL8)**p1)**p2
            end if
            call undefine(sn%scrfun)
            call undefine(sn%scrcoeffile)
            call undefine(sn%scrfile)
            call undefine(sn%scatfile)
            if(model == 1) then
               sn%scrcoeffile = 'screencoefs/NLHLIN_python_coeffs.dat'
               sn%scrlen = 'NLHlin'
               native_norm = 0.4685_REAL8/sqrt(sqrt(real(z1,REAL8))+sqrt(real(z2,REAL8)))
            else
               sn%scrfun = 'ZBLuniv'
               sn%scrlen = 'ZBL'
               native_norm = 0.4685_REAL8/(real(z1,REAL8)**0.23_REAL8+real(z2,REAL8)**0.23_REAL8)
            end if
            call initScreen()
            write(filename,'(a,a,a,i0,a,i0,a)') 'screen_', trim(label), '_c', config, '_z', z1, '.txt'
            open(newunit=unit,file=trim(filename),status='replace')
            write(unit,'(a)') '# columns: input_radius phi dphi_dinput_radius'
            do i = 0, 100
               if(model == 1) then
                  coordinate = (3._REAL8/rnorm)*real(i,REAL8)/100._REAL8
                  radius = coordinate*rnorm/native_norm
                  call screen(radius,1,2,phi,dphi)
                  ! Chain rule: convert the native reduced-radius derivative
                  ! to the normalization used by the requested test grid.
                  dphi = dphi*rnorm/native_norm
               else
                  coordinate = 3._REAL8*real(i,REAL8)/100._REAL8
                  radius = coordinate/native_norm
                  call screen(radius,1,2,phi,dphi)
                  ! ZBL output coordinates are physical radii in Angstrom.
                  dphi = dphi/native_norm
               end if
               write(unit,'(3es25.16)') coordinate, phi, dphi
            end do
            close(unit)
         end do
      end do
   end do
end program screen_reference
