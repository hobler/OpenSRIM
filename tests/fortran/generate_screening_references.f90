! Generate screening functions and derivatives for extended element-pair tests.
!
! Args:
!     None. The model, normalization and element-pair loops are defined below.
! Inputs:
!     The canonical NLHlin coefficient table staged by conftest.py in screencoefs/.
! Outputs:
!     552 TXT files: two models times three configurations times 92 element pairs.
!     Each file has 101 rows: input radius, phi, derivative with respect to that radius.
! Notes:
!     Grids and unit conversions match test_extended_element_pairs.compare_screening.
!     conftest.py collects the files and creates the screening integrity manifest.
program generate_screening_references
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
   ! Allocate the two element slots and screening configuration once. Each loop
   ! iteration below assigns their atomic numbers and reinitializes screening.
   nelem = 2
   allocate(elem(2), sn%scrfun(2,2), sn%scrcoeffile(2,2), sn%scrfile(2,2), &
            sn%scrlen(2,2), sn%scatfile(2,2), sn%lquad(2,2))
   ! Set the quadrature flags for the shared configuration; this program only
   ! evaluates screening and its derivative, not scattering integrals.
   sn%lquad = .false.
   do model = 1, 2
      if(model == 1) label = 'nlhlin'
      if(model == 2) label = 'zbl'
      ! Configurations match test_extended_element_pairs.CONFIGURATIONS. The first two use
      ! equal-element pairs; the third fixes the second element to copper.
      ! Configuration 1 uses a unit test length. Configurations 2 and 3 use
      ! 0.4685 / (Z1**p1 + Z2**p1)**p2 as the requested normalization length.
      do config = 1, 3
         if(config == 1) then
            p1 = 0; p2 = 0
         else if(config == 2) then
            p1 = 0.23_REAL8; p2 = 1
         else
            p1 = 0.5_REAL8; p2 = 2._REAL8/3._REAL8
         end if
         ! Cover every atomic number from hydrogen to uranium. The third configuration
         ! keeps the partner fixed at copper (Z2=29); the others use equal-element pairs.
         do z1 = 1, 92
            z2 = z1
            if(config == 3) z2 = 29
            elem(1)%z = z1; elem(2)%z = z2
            ! rnorm specifies the test coordinate scaling; native_norm below is the
            ! length scale expected by the production screening routine.
            if(config == 1) then
               rnorm = 1
            else
               rnorm = 0.4685_REAL8/(real(z1,REAL8)**p1+real(z2,REAL8)**p1)**p2
            end if
            ! Reset all optional screening/table selectors so the previous model cannot
            ! leave a conflicting selection in the shared configuration.
            call undefine(sn%scrfun)
            call undefine(sn%scrcoeffile)
            call undefine(sn%scrfile)
            call undefine(sn%scatfile)
            ! Select tabulated NLHlin coefficients or universal ZBL and compute the
            ! corresponding native normalization length in Angstrom.
            if(model == 1) then
               sn%scrcoeffile = 'screencoefs/NLHLIN_python_coeffs.dat'
               sn%scrlen = 'NLHlin'
               native_norm = 0.4685_REAL8/sqrt(sqrt(real(z1,REAL8))+sqrt(real(z2,REAL8)))
            else
               sn%scrfun = 'ZBLuniv'
               sn%scrlen = 'ZBL'
               native_norm = 0.4685_REAL8/(real(z1,REAL8)**0.23_REAL8+real(z2,REAL8)**0.23_REAL8)
            end if
            ! Rebuild screening for the current model and atomic numbers before evaluating it.
            call initScreen()
            ! Encode both atomic numbers and the comparison configuration in the filename.
            ! I0 writes integers without padding; TRIM removes fixed-length string padding.
            write(filename,'(a,a,a,i0,a,i0,a,i0,a)') &
               'screen_',trim(label),'_z',z1,'_z',z2,'_norm_c',config,'_grid_screening.txt'
            open(newunit=unit,file=trim(filename),status='replace')
            write(unit,'(a)') '# columns: input_radius phi dphi_dinput_radius'
            ! Sample 101 points covering physical radii from zero to 3 Angstrom.
            do i = 0, 100
               if(model == 1) then
                  ! NLHlin output uses the test-reduced radius r/rnorm. Convert it to
                  ! the native reduced radius r/native_norm before calling screen.
                  coordinate = (3._REAL8/rnorm)*real(i,REAL8)/100._REAL8
                  radius = coordinate*rnorm/native_norm
                  call screen(radius,1,2,phi,dphi)
                  ! Chain rule: convert the native reduced-radius derivative
                  ! to the normalization used by the requested test grid.
                  dphi = dphi*rnorm/native_norm
               else
                  ! ZBL output uses physical radii, independently of the configuration.
                  ! The production API still receives the radius divided by native_norm.
                  coordinate = 3._REAL8*real(i,REAL8)/100._REAL8
                  radius = coordinate/native_norm
                  call screen(radius,1,2,phi,dphi)
                  ! ZBL output coordinates are physical radii in Angstrom.
                  dphi = dphi/native_norm
               end if
               ! Write the test coordinate and outputs after converting the derivative.
               ! The Python tests then compare all backends in these same units.
               write(unit,'(3es25.16)') coordinate, phi, dphi
            end do
            close(unit)
         end do
      end do
   end do
end program generate_screening_references
