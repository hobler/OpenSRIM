! Shared output routines for the classic ZBL and NLHlin comparison drivers.
! Calls production scattering routines; no independent solver is implemented here.
!
! Public routines:
!     write_chi_data(model): write integration-coordinate/chi tables.
!     write_comparison_data(model): write angle and time-integral tables.
! Preconditions:
!     The caller has initialized screening for boron/silicon in element slots 1,2.
!     Reduced radii and impact parameters use the selected model normalization.
! Outputs:
!     Files in the current working directory; conftest.py publishes and hashes them.
module scattering_reference_writers
   use Const_mod, only: REAL8
   use Sn_CMscatter_mod, only: cmAngle, calcPhiChi, cmScatterLegendre
   use Sn_Screen_mod, only: get_rmax
   implicit none
contains
   ! Evaluate production chi for four radii on regular and endpoint-focused grids.
   !
   ! Args:
   !     model (character, in): Lowercase filename label, 'zbl' or 'nlhlin'.
   !         Screening must already be initialized; this argument does not select it.
   ! Outputs:
   !     Twelve TXT files, one per radius/grid combination, with coordinate and chi.
   !     The coordinate is u except for near_one, where it is 1-u.
   !     phi is returned by calcPhiChi but is not written by this diagnostic.
   subroutine write_chi_data(model)
      character(*), intent(in) :: model
      real(REAL8) :: radii(4), r0, coordinate, u, phi, chi
      integer :: ir, grid, i, count, unit
      character(16) :: suffix
      character(5) :: radius_labels(4)
      character(256) :: filename
      ! Keep numeric radii and readable filename labels in the same order.
      ! These values match scattering_support.R0_VALUES.
      radius_labels = ['1e-02', '1e-01', '1e+00', '1e+01']
      radii = [0.01_REAL8, 0.1_REAL8, 1._REAL8, 10._REAL8]
      do ir=1,4
         r0=radii(ir)
         do grid=1,3
            ! Choose the sample count and filename suffix for the integration grid.
            select case(grid)
            case(1)
               suffix='regular'; count=999
            case(2)
               suffix='near_zero'; count=101
            case(3)
               suffix='near_one'; count=101
            end select
            ! Include the actual r0 value in the name. newunit assigns a free file unit;
            ! status=replace overwrites this case within the temporary generation directory.
            filename = 'chi_'//model//'_z5_z14_norm_native_r0_'//radius_labels(ir)// &
                       '_grid_'//trim(suffix)//'.txt'
            open(newunit=unit,file=trim(filename),status='replace')
            write(unit,'(a)') '# columns: coordinate chi (production calcPhiChi)'
            do i=1,count
               ! The regular grid excludes both endpoints. Logarithmic grids resolve
               ! small u or small 1-u, depending on which coordinate is passed to the API.
               if(grid==1) then
                  coordinate=0.001_REAL8+0.998_REAL8*real(i-1,REAL8)/real(count-1,REAL8)
               else
                  coordinate=10._REAL8**(-10._REAL8+10._REAL8*real(i-1,REAL8)/real(count-1,REAL8))
                  ! Replace the final value 1 by its immediate floating-point predecessor,
                  ! avoiding u=1 in near_zero and u=0 in near_one.
                  if(i==count) coordinate=nearest(1._REAL8,-1._REAL8)
               end if
               ! Convert the stored coordinate to the integration variable u.
               u=coordinate
               if(grid==3) u=1._REAL8-coordinate
               ! Evaluate the original routine for element slots 1 and 2.
               ! ES25.16 writes scientific notation with 16 digits after the decimal point.
               call calcPhiChi(u,r0,1,2,phi,chi)
               write(unit,'(2es25.16)') coordinate,chi
            end do
            close(unit)
         end do
      end do
   end subroutine write_chi_data

   ! Generate the standard angle and time-integral comparison grids.
   !
   ! Args:
   !     model (character, in): 'zbl' or 'nlhlin', matching initialized screening.
   !         Supplies filename labels and enables finite-cutoff cases for NLHlin.
   ! Outputs:
   !     theta files: reduced impact parameter p, scattering angle theta.
   !     scatter files: p, theta, reduced time integral tau.
   ! Notes:
   !     e is reduced energy; p and rmax use the model's reduced length units.
   !     theta is in radians. Each file contains one energy and quadrature order.
   !     Most cases use four Legendre points; the tau comparison also uses 32.
   subroutine write_comparison_data(model)
      character(*), intent(in) :: model
      real(REAL8) :: energies(7), upper(3), e, p, theta, tau, rmax
      integer :: ie, ip, group, unit, energy_index
      character(8) :: labels(7)
      character(256) :: filename
      ! Define the seven energies and their filename labels. The three upper bounds
      ! belong to the error grids at energies 1e-4, 0.1 and 100, respectively.
      energies = [1e-4_REAL8, 1e-3_REAL8, 1e-2_REAL8, 1e-1_REAL8, 1._REAL8, 10._REAL8, 100._REAL8]
      labels = ['1e-04', '1e-03', '1e-02', '1e-01', '1e+00', '1e+01', '1e+02']
      upper = [30._REAL8, 12._REAL8, 4._REAL8]
      ! Obtain the cutoff from the initialized screening model. Only NLHlin uses
      ! the finite-cutoff dataset below.
      rmax = get_rmax(1, 2)
      do ie = 1, 7
         e = energies(ie)
         ! Head-on collisions have p=0. Request the corrected Legendre angle directly
         ! so that this limiting case follows the production API.
         filename = 'theta_'//model//'_z5_z14_norm_native_e'//trim(labels(ie))// &
                    '_n4_grid_head_on.txt'
         open(newunit=unit, file=trim(filename), status='replace')
         call cmAngle(e, 0._REAL8, 1, 2, theta, method='LegendreChiCorr', nabsc=4)
         write(unit,'(2es25.16)') 0._REAL8, theta
         close(unit)
         ! Probe the cutoff, its next representable larger value, and a distant point.
         ! These cases check the zero-scattering behavior at and beyond the range.
         if(model == 'nlhlin') then
            filename = 'scatter_'//model//'_z5_z14_norm_native_e'//trim(labels(ie))// &
                    '_n4_grid_cutoff.txt'
            open(newunit=unit, file=trim(filename), status='replace')
            do ip = 1, 3
               if(ip == 1) p = rmax
               if(ip == 2) p = nearest(rmax, 1._REAL8)
               if(ip == 3) p = 30._REAL8
               call cmScatterLegendre(e, p, 1, 2, theta, tau, nabsc=4)
               write(unit,'(3es25.16)') p, theta, tau
            end do
            close(unit)
         end if
         ! grid_ratio names the tau-versus-theta plot dataset; the file
         ! stores both coordinates, not their quotient.
         ! The ratio plot uses the first six energies and 201 points from 0.02 to 30.
         if(ie <= 6) then
            filename = 'scatter_'//model//'_z5_z14_norm_native_e'//trim(labels(ie))// &
                    '_n4_grid_ratio.txt'
            open(newunit=unit, file=trim(filename), status='replace')
            do ip = 0, 200
               p = 0.02_REAL8 + (30._REAL8-0.02_REAL8)*real(ip,REAL8)/200._REAL8
               call cmScatterLegendre(e, p, 1, 2, theta, tau, nabsc=4)
               write(unit,'(3es25.16)') p, theta, tau
            end do
            close(unit)
         end if
      end do
      ! Native angle outputs for the sin-squared plot comparison.
      ! Use 201 points from 0.01 to 30 for each energy. Store theta itself;
      ! pytest computes sin(theta/2)**2 when comparing the corresponding quantity.
      do ie=1,7
         e=energies(ie)
         filename='theta_'//model//'_z5_z14_norm_native_e'//trim(labels(ie))// &
                    '_n4_grid_sin.txt'
         open(newunit=unit,file=trim(filename),status='replace')
         do ip=0,200
            p=0.01_REAL8+(30._REAL8-0.01_REAL8)*real(ip,REAL8)/200._REAL8
            call cmAngle(e,p,1,2,theta,method='LegendreChiCorr',nabsc=4)
            write(unit,'(2es25.16)') p,theta
         end do
         close(unit)
      end do
      ! The three energy-dependent p ranges match P_VALUES_ERROR in the tests.
      do group = 1, 3
         ! Select energy indices 1, 4 and 7: 1e-4, 0.1 and 100.
         ! Each error grid has 101 points and the matching energy-dependent upper bound.
         energy_index = 1 + 3*(group-1)
         e = energies(energy_index)
         filename='theta_'//model//'_z5_z14_norm_native_e'//trim(labels(energy_index))// &
                    '_n4_grid_error.txt'
         open(newunit=unit,file=trim(filename),status='replace')
         do ip=0,100
            p=0.02_REAL8+(upper(group)-0.02_REAL8)*real(ip,REAL8)/100._REAL8
            call cmAngle(e,p,1,2,theta,method='LegendreChiCorr',nabsc=4)
            write(unit,'(2es25.16)') p,theta
         end do
         close(unit)
         ! On the same error grid, evaluate both theta and tau with four Legendre points.
         filename = 'scatter_'//model//'_z5_z14_norm_native_e'//trim(labels(energy_index))// &
                    '_n4_grid_error.txt'
         open(newunit=unit, file=trim(filename), status='replace')
         do ip = 0, 100
            p = 0.02_REAL8 + (upper(group)-0.02_REAL8)*real(ip,REAL8)/100._REAL8
            call cmScatterLegendre(e, p, 1, 2, theta, tau, nabsc=4)
            write(unit,'(3es25.16)') p, theta, tau
         end do
         close(unit)
         ! The n32 file records the same physical inputs at a higher quadrature order.
         filename = 'scatter_'//model//'_z5_z14_norm_native_e'//trim(labels(energy_index))// &
                    '_n32_grid_error.txt'
         open(newunit=unit, file=trim(filename), status='replace')
         do ip=0,100
            p=0.02_REAL8+(upper(group)-0.02_REAL8)*real(ip,REAL8)/100._REAL8
            call cmScatterLegendre(e,p,1,2,theta,tau,nabsc=32)
            write(unit,'(3es25.16)') p,theta,tau
         end do
         close(unit)
      end do
   end subroutine write_comparison_data
end module scattering_reference_writers
