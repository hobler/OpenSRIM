! Generate ZBL chi and scattering references for the default boron/silicon pair.
!
! Args:
!     None. Model, element pair and comparison grids are fixed by this driver.
! Inputs:
!     Working-directory inputs are staged by conftest.generate_fortran().
!     The NLHlin driver uses the shared coefficient table in screencoefs/.
! Outputs:
!     TXT files in the working directory, written by scattering_reference_writers.
!     conftest.py collects the files and writes their integrity manifest.
! Notes:
!     Uses production screening and scattering APIs; comparisons are in pytest.
program generate_zbl_scattering_references
   use Const_mod, only: REAL8
   use Param_Setup_mod, only: nelem
   use Param_Elem_mod, only: elem
   use Param_Sn_mod, only: sn
   use Undefine_mod
   use Sn_Screen_mod, only: initScreen
   use scattering_reference_writers, only: write_chi_data, write_comparison_data
   implicit none
   ! Allocate two element slots and the pair-indexed screening configuration.
   ! The API indices 1 and 2 below refer to these slots, not atomic numbers.
   nelem=2
   allocate(elem(2), sn%scrfun(2,2), sn%scrcoeffile(2,2), sn%scrfile(2,2), &
            sn%scrlen(2,2), sn%scatfile(2,2), sn%lquad(2,2))
   ! Set boron (Z=5) as projectile and silicon (Z=14) as target.
   elem(1)%z=5; elem(2)%z=14
   ! Initialize the global quadrature flags. The comparison routines explicitly
   ! select Legendre quadrature and its order in each scattering API call.
   sn%lquad=.false.
   ! Clear optional model and table selectors before choosing the reference model.
   ! undefine sets the library sentinel; it does not deallocate these arrays.
   call undefine(sn%scrfun)
   call undefine(sn%scrcoeffile)
   call undefine(sn%scrfile)
   call undefine(sn%scatfile)
   ! Select the universal ZBL screening function and its reduced-length normalization.
   sn%scrfun='ZBLuniv'
   sn%scrlen='ZBL'
   ! Build the screening coefficients for the configured element pair.
   call initScreen()
   ! Write chi on three integration grids for each reference radius, then write
   ! angle/time-integral data. The lowercase argument supplies the filename label.
   call write_chi_data('zbl')
   call write_comparison_data('zbl')
end program generate_zbl_scattering_references
