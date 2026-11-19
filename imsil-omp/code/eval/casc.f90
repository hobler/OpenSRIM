module Eval_Casc_mod
!
!  Stores and analyzes vacancies and interstitials generated in one ion's recoil cascade
!
!  The resulting point defects are either scored in a histogram (using evalStat_mod) or
!  written to a file (PDI, PDV). The latter requires single-threaded execution of the code.
!
!  Used Modules:
!  ============
!
   use Const_mod,          only: INT4,REAL8,  &
                                 ZERO, ONE, THIRD, BIG,  &
                                 PI
   use Vector_mod,         only: abs, abs2
   use Stack_mod,          only: LogicalStack, Int4Stack, Real8Stack, VectorStack,  &
                                 newStack, clearStack, popStack, pushStack, sizeOfStack
   use Types_mod,          only: Proj_type, Ipos_type

   use Open_mod,           only: openFile

   use Param_Setup_mod,    only: nelem
   use Param_Mat_mod,      only: mat, coord_xtal, irxtal
   use Param_Dam_mod,      only: dam
   use Param_Out_mod,      only: out

   use Eval_Stat_mod,      only: stat2x, scoreIV
   use Geom_CoordSys_mod,  only: toCoordSys
   use Xtal_Lattice_mod,   only: getNearestSite
!
!  Default Type and Access:
!  =======================
!
   implicit none
   private
   save
!
!  Public Procedures:
!  =================
!
   public :: casc10                 ! Initialize at beginning of simulation
!            ======                  
   public :: casc11                 ! Initialize at beginning of molecule cascade
!            ======
   public :: casc21                 ! Evaluate at beginning of trajectory
!            ======
   public :: casc23                 ! Evaluate at end of trajectory of a stopped ion
!            ======
   public :: casc26                 ! Actions at end of cascade
!            ======
!
!  Private Data:
!  ============
!
   type(VectorStack), allocatable :: posvStack(:)  ! stack of vacancies
   type(VectorStack), allocatable :: posiStack(:)  ! stack of interstitials
   type(Int4Stack), allocatable :: irStack(:)      ! stack of region indices
   type(Real8Stack), allocatable :: edamStack(:)   ! stack of damage energy
   type(LogicalStack), allocatable :: subStack(:)  ! stack of subthreshold indicator
!$omp threadprivate( posvStack, posiStack, irStack, edamStack, subStack )
!
contains
!
!=======================================================================
!
   subroutine casc10( )
!             ======
!
!     Initialize variables at the beginning of the whole simulation
!
!  Arguments: 
!
!     None
!
!  Local Variables:
!
      integer(INT4) :: ielem
!
!  Start of Executable Code:
!  ------------------------
!
!  Open the point defect and amorphous pocket files
!
      if( dam%lcap ) then
         if( out%ltra(1) ) then
            call openFile( 34, 'PDV', 'formatted', 'write' )
!                ========   i   i         i          i
         end if
         if( out%ltra(3) ) then
            call openFile( 35, 'PDI', 'formatted', 'write' )
!                ========   i    i         i          i
            if( irxtal /= 0 )  write( 35, * ) mat(irxtal)%wafer(:)
         end if
      end if
!
      if( dam%lpock .and. out%ltra(3) ) then
         call openFile( 24, 'POK', 'formatted', 'write' )
!             ========   i    i          i         i
      end if
!
!  Initialize posvStack
!
      allocate( posvStack( nelem ) )
      do ielem = 1, nelem
         call newStack( posvStack(ielem) )
!             ========       o
      end do
!
!  Initialize posiStack
!
      allocate( posiStack( nelem ) )
      do ielem = 1, nelem
         call newStack( posiStack(ielem) )
!             ========       o
      end do
!
!  Initialize subStack
!
      allocate( subStack( nelem ) )
      do ielem = 1, nelem
         call newStack( subStack(ielem) )
!             ========       o
      end do
!
!  Initialize edamStack
!
      allocate( edamStack( nelem ) )
      do ielem = 1, nelem
         call newStack( edamStack(ielem) )
!             ========       o
      end do

      if( dam%lpock ) then
!
!  Initialize irStack
!
         allocate( irStack( nelem ) )
         do ielem = 1, nelem
            call newStack( irStack(ielem) )
!                ========       o
         end do

      end if

   end subroutine casc10
!
!=======================================================================
!
   subroutine casc11( )
!             ======
!
!     Initialize variables at the beginning of the simulation of a recoil
!     cascade
!
!  Arguments: 
!
!     None
!
!  Local Variables:
!
      integer(INT4) :: ielem
!
!  Start of Executable Code:
!  ------------------------
!
!     Stack should always be clear at the beginning of a recoil cascade, if
!     the previous cascade has been handled properly. But just to make sure
!     we clear it again.
!
      do ielem = 1, nelem
         call clearStack( posvStack(ielem) )
!            ==========        io
      end do

      do ielem = 1, nelem
         call clearStack( posiStack(ielem) )
!             ==========       io
      end do

      do ielem = 1, nelem
         call clearStack( subStack(ielem) )
!             ==========       io
      end do

      do ielem = 1, nelem
         call clearStack( edamStack(ielem) )
!             ==========       io
      end do

      if( dam%lpock ) then

         do ielem = 1, nelem
            call clearStack( irStack(ielem) )
!                ==========       io
         end do

      end if

   end subroutine casc11
!
!=======================================================================
!
   subroutine casc21( proj )
!             ======   i   
!
!     Evaluate at the beginning of a trajectory
!
!  Arguments:
!
      type(Proj_type), intent(in) :: proj       ! vacancy
!
!  Local Variables:
!
!     None
!
!  Start of Executable Code:
!  ------------------------
!                                               ! store the vacancy in cascade history
! 
      call pushStack( proj%pos(:), posvStack( proj%ielem ) )
!          =========       i                io

   end subroutine casc21
!
!=======================================================================
!
   subroutine casc23( proj, sub )
!             ======   i     i
!
!     Evaluate at the end of a trajectory
!
!  Arguments:
!
      type(Proj_type), intent(in) :: proj       ! interstitial
      logical, intent(in) :: sub                ! indicates projectile has not really started
                                                !   e.g., recoil in a subthreshold collision
!
!  Local Variables:
!
!     None
!
!  Start of Executable Code:
!  ------------------------
!
      call pushStack( proj%pos(:), posiStack(proj%ielem) )
!          =========      i            io

      call pushStack( sub, subStack(proj%ielem) )
!          =========   i       io

      call pushStack( max( ZERO, proj%e ), edamStack(proj%ielem) )
!          =========          i             io

      if( dam%lpock ) then

         call pushStack( proj%id, irStack(proj%ielem) )
!             =========     i            io

      end if

   end subroutine casc23
!
!=======================================================================
!
   subroutine casc26( )
!             ======
!
!     Perform actions after simulation of a collision cascade
!
!     - identify amorphous pockets
!
!     - reduce interstitial and vacancy list due to recombination,
!       update the vacancy and interstitial histograms, and the
!       pair distribution histogram. 
!
!  Arguments: 
!
!     None
!
!  Local Variables:
!                                               ! Note: objects of Interstitial_type
!                                               ! include subthreshold recoils
      type :: Interstitial_type
   !          =================
         logical :: valid                       ! lcap=t: interstitial still exists
                                                ! lpock=t: interstitial processed
         real(REAL8) :: pos(3)                  ! position of interstitial
         integer(INT4) :: ielem                 ! target atom type
         real(REAL8) :: edam                    ! recoil energy when stopped
         real(REAL8) :: ndisp                   ! number of displaced atoms
         real(REAL8) :: nmelt                   ! number of molten atoms
         real(REAL8) :: rmelt                   ! radius of potentially molten pocket
         logical :: subthreshold                ! indicates generation in subthreshold collision
         type(Int4Stack) :: neighStack          ! stack of neighboring interstitials
      end type Interstitial_type

      type(Interstitial_type), allocatable :: interstitial(:)

      type :: Vacancy_type
!             ============
         logical :: valid                       ! indicates if vacancy still exists
         real(REAL8) :: pos(3)                  ! position of vacancy
         integer(INT4) :: ielem                 ! target atom type
      end type Vacancy_type

      type(Vacancy_type), allocatable :: vacancy(:)

      type(Int4Stack) :: pocketAtomStack        ! all pocket atoms
      type(Int4Stack) :: tmpAtomStack           ! pocket atoms to be processed

      type(Proj_type) :: proj
      real(REAL8) :: pos(3)
      real(REAL8) :: posi(3)
      logical :: sub

      logical :: newvac, isEmpty
      integer(INT4) :: ir, ielem, ii, ii1, ii2, iimin, iv
      integer(INT4) :: ni, nv, nvalid
      type(Ipos_type) :: ipos
      real(REAL8) :: dist, ddami, ddamv, epock, fdisp, r2, rmin2
!
!  Start of Executable Code:
!  ------------------------
!
!  Get sizes of stacks
!
      nv = 0
      do ielem = 1, nelem
         nv = nv + sizeOfStack( posvStack(ielem) )
!                  ===========        i
          
      end do

      ni = 0
      do ielem = 1, nelem
         ni = ni + sizeOfStack( posiStack(ielem) )
!                  ===========       i
      end do
!
!  Allocate interstitial and vacancy arrays
!
      allocate( vacancy( nv ) )
      allocate( interstitial( ni ) )
!
!  Convert vacancy stack to array
!
      iv = 0
      do ielem = 1, nelem
         do
            call popStack( posvStack(ielem), pos(:), isEmpty )
!                ========       io            o         o
            if( isEmpty  )  exit
            iv = iv + 1
            vacancy(iv)%pos(:) = pos(:)
            vacancy(iv)%ielem = ielem
            vacancy(iv)%valid = .true.
         end do
      end do
      if( iv /= nv )  error stop 'casc26: error during copying from stack to array'
!
!  Convert interstitial stack to arrays
!
      ii = 0
      nvalid = 0                                ! debug
      do ielem = 1, nelem
         do
            call popStack( posiStack(ielem), pos(:), isEmpty)
!                ========       io            o         o
            call popStack( subStack(ielem), sub, isEmpty )
!                ========      io            o      o
            if( isEmpty  )  exit
            ii = ii + 1
            interstitial(ii)%pos(:) = pos(:)
            interstitial(ii)%ielem = ielem

            call popStack( edamStack(ielem), interstitial(ii)%edam, isEmpty )
!                ========      io                  o                    o
            if( isEmpty )  error stop 'casc26: error during copying from stack to array'
!
            if( dam%lcap ) then
               interstitial(ii)%valid = .not. sub
               if( interstitial(ii)%valid )  nvalid = nvalid + 1  ! debug
            elseif( dam%lpock ) then
               interstitial(ii)%valid = .false.
               call popStack( irStack(ielem), ir, isEmpty )
!                   ========       io         o      o
               if( isEmpty )  error stop 'casc26: error during copying from stack to array'
               interstitial(ii)%subthreshold = sub
               interstitial(ii)%nmelt = interstitial(ii)%edam / dam%emelt
               interstitial(ii)%rmelt =  &
                    ( 3*interstitial(ii)%nmelt / ( 4*PI*mat(ir)%dens) ) ** THIRD
            end if
         end do
      end do
      if( ii /= ni )  error stop 'casc26: error during copying from stack to array'
!
!  Write amorphous pocket file
!
      if( dam%lpock .and. out%ltra(3) ) then
         do ii = 1, ni
            write( 24, '(4f10.2)' )  interstitial(ii)%pos(:), interstitial(ii)%rmelt
         end do
      end if
!
!  Treat nuclear energy deposition
!
      do ii = 1, ni
         proj%pos(:) = interstitial(ii)%pos(:)
         proj%dir(:) = [ZERO,ZERO,ONE]
         proj%e = interstitial(ii)%edam         ! maybe unnecessary
         proj%ielem = interstitial(ii)%ielem
         proj%id = ir                           ! maybe incorrect
         call stat2x( proj, ZERO, ZERO, interstitial(ii)%edam, 3 )
!             ======   i     i     i             i
      end do
!
!-----------------------------------------------------------------------
!
!  Identify amorphous pockets
!
      if( dam%lpock ) then
!
!  Set up list of neighbors that potentially may be members of the same 
!  amorphous pocket, based on their maximum melting radii
!
         do ii = 1, ni                          ! allocate neighbor stack
            call newStack( interstitial(ii)%neighStack )
!                ========                o
         end do
         do ii = 1, ni                          ! check overlap of spheres and push
            do ii1 = ii+1, ni                   !    neighbor stacks
               if( abs2( interstitial(ii)%pos(:) - interstitial(ii1)%pos(:) ) <  &
!                  ====
                   ( interstitial(ii)%rmelt + interstitial(ii1)%rmelt )**2 ) then
                  call pushStack( ii1, interstitial(ii)%neighStack )
!                      =========   i                io
                  call pushStack( ii, interstitial(ii1)%neighStack )
!                      =========   i                io
               end if
            end do
         end do
!
!  Create stacks for atoms belonging to a pocket
!                                            ! pocketAtomStack contains all atoms in same
         call newStack( pocketAtomStack )    !   pocket as ii.
!             ========         o             ! tmpAtomStack contains only those pocket atoms
         call newStack( tmpAtomStack )       !   whose neighbor list has not yet been
!             ========        o              !   checked (and therefore still exists). 

         Pocket_loop: do ii = 1, ni          ! by skipping already processed interstitials,
                                             !    this is a loop over pockets
            if( interstitial(ii)%valid )  cycle
!
!  Make sure pocket atom stacks are clear
!
            call clearStack( pocketAtomStack )
!                ==========         o
            call clearStack( tmpAtomStack )
!                ==========       o
!
!  First atom (ii) always belongs to the tentative pocket
!
            call pushStack( ii, pocketAtomStack )
!                =========   i       io
            call pushStack( ii, tmpAtomStack )
!                =========   i      io
            interstitial(ii)%valid = .true.     ! mark as processed
            epock = interstitial(ii)%edam
!
!  Trace the neighbor atoms
!
            do
               call popStack( tmpAtomStack, ii1, isEmpty ) ! take processed atom
!                   ========       io        o      o
               if( isEmpty  )  exit
               do                                       ! process neighbors
                  call popStack ( interstitial(ii1)%neighStack, ii2, isEmpty )
!                      ========               io                 o      o
                  if( isEmpty  )  exit
                  if( interstitial(ii2)%valid )  cycle ! don't reprocess
                  call pushStack( ii2, pocketAtomStack )
!                      =========   i        io
                  call pushStack( ii2, tmpAtomStack )
!                      =========   i       io
                  interstitial(ii2)%valid = .true.
                  epock = epock + interstitial(ii2)%edam
               end do
            end do
!
!  Now we have pocketAtomStack and the pocket energy epock.
!  Determine the number of displaced atoms
!
            if( epock < dam%epockmin ) then  ! pocket has too little energy convert to
                                             !    interstitials
               do
                  call popStack( pocketAtomStack, ii1, isEmpty )
!                      ========         io         o      o
                  if( isEmpty )  exit
                  if( interstitial(ii1)%subthreshold ) then
                     interstitial(ii1)%ndisp = 0
                  else
                     interstitial(ii1)%ndisp = 2   ! assuming split interstitials
                  end if
               end do

            else                                ! consider recrystallization

               fdisp = ( ONE - dam%fpock / (epock/dam%emelt)**THIRD ) ** 3

               do
                  call popStack( pocketAtomStack, ii1, isEmpty )
!                      ========         io        o     o
                  if( isEmpty )  exit
                  interstitial(ii1)%ndisp = interstitial(ii1)%nmelt * fdisp
               end do

            end if

         end do Pocket_loop
!
!  Increase damage histograms
!
         do ii = 1, ni
            if( interstitial(ii)%valid ) then
               proj%pos(:) = interstitial(ii)%pos(:)
               proj%dir(:) = [ZERO,ZERO,ONE]
               proj%e = ZERO
               proj%ielem = interstitial(ii)%ielem
               proj%id = ir                     ! maybe incorrect
               call stat2x( proj, interstitial(ii)%ndisp, ZERO, ZERO, 3 )
!                   ======   i                  i          i
            else
               error stop 'casc26: unclassified interstitial found'
            end if
         end do

      end if
!
!-----------------------------------------------------------------------
!
!  Recombine vacancies with the nearest interstitial in the cascade
!                                               ! note: it would be even more accurate to set
!                                               !   up a list of distances and recombine
!                                               !   starting with the smallest
      if( dam%lcap ) then
         do iv = 1, nv            
            rmin2 = BIG                      ! only vacancies and interstitials of 
            do ii = 1, ni                    !   same type can recombine
               if( .not.interstitial(ii)%valid )  cycle
               if( interstitial(ii)%ielem /= vacancy(iv)%ielem )  cycle
               r2 = abs2( interstitial(ii)%pos(:)-vacancy(iv)%pos(:) )
!                   ====
               if( r2 < rmin2 ) then
                  rmin2 = r2
                  iimin = ii
               end if
            end do
            if( rmin2 < dam%rcap2 ) then
               vacancy(iv)%valid = .false.
               interstitial(iimin)%valid = .false.
            end if
         end do
         !
         !  Write point defect file
         !
         if( out%ltra(1) ) then
            do iv = 1, nv
               if(vacancy(iv)%valid ) then
                  write( 34,'(t2,3f15.6,i10)' ) &
                        vacancy(iv)%pos(:), -vacancy(iv)%ielem
                  if( irxtal /= 0 ) then
                     call toCoordSys( coord_xtal, vacancy(iv)%pos(:), posi(:) )
!                         ==========            i                 i              o
                     write( 35,'(t2,4i10)' ) &
                           nint(posi(:)), -vacancy(iv)%ielem
                  end if
               endif
            enddo
         end if
         if( out%ltra(3) ) then
            do ii = 1, ni
               if( interstitial(ii)%valid ) then
                  write( 34,'(t2,3f15.6,i10)' ) &
                        interstitial(ii)%pos(:), interstitial(ii)%ielem
                  if( irxtal /= 0 ) then
                     call getNearestSite( interstitial(ii)%pos(:), ir, ipos )
!                         ==============           i               i    o
                     write( 35, '(t2,4i10)' ) &
                           ipos%ipos(:), interstitial(ii)%ielem
                  end if
               endif
            enddo
         endif
!
!  Consider recombination with defects created by previous cascades and
!  update the damage and pair distribution histograms
!
         do iv = 1, nv
            if( .not.vacancy(iv)%valid )  cycle ! consider only "valid" vacancies
            ddamv = ONE
            ddami = ZERO
            proj%pos(:) = vacancy(iv)%pos(:)
            proj%dir(:) = [ZERO,ZERO,ONE]
            proj%e = ZERO
            proj%ielem = vacancy(iv)%ielem
            proj%id = ir                        ! maybe incorrect
            call stat2x( proj, ddami, ddamv, ZERO, 3 )
!                ======   i      i      i
            if( out%lhisiv ) then               ! update pair distribution histogram
               newvac = .true.
               do ii = 1, ni
                  if( .not.interstitial(ii)%valid )  cycle
                  if( interstitial(ii)%ielem /= vacancy(iv)%ielem )  cycle
                  dist = abs( interstitial(ii)%pos(:)-vacancy(iv)%pos(:) )
!                        ===
                  call scoreIV( vacancy(iv)%pos(:), newvac, dist, vacancy(iv)%ielem )
!                      =======          i             i      i            i
                  newvac = .false.
               end do
            end if
         end do

         do ii = 1, ni                          ! consider only "valid" interstitials
            if( .not.interstitial(ii)%valid )  cycle 
            ddami = ONE
            ddamv = ZERO
            proj%pos(:) = interstitial(ii)%pos(:)
            proj%dir(:) = [ZERO,ZERO,ONE]
            proj%e = ZERO
            proj%ielem = interstitial(ii)%ielem
            proj%id = ir                        ! maybe incorrect
            call stat2x( proj, ddami, ddamv, ZERO, 3 )
!                ======   i      i      i
         end do

      end if
!
!-----------------------------------------------------------------------
!
!  Deallocate arrays
!
      deallocate( interstitial )
      deallocate( vacancy )
!
!-----------------------------------------------------------------------
!
   end subroutine casc26
!
!=======================================================================
!
end module Eval_Casc_mod
