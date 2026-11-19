program main
!
!  Monte Carlo simulation of ion implantation and sputtering into crystalline and amorphous
!  targets
!
!  Used Modules:
!  ============
!
!$ use omp_lib,            only: omp_get_thread_num
   use Const_mod,          only: INT4

   use Param_mod,          only: readParams, initParams, writeParams
   use Param_Setup_mod,    only: ldyn, nion, nion1, nionchunk, nthread
   use Param_Mat_mod,      only: irxtal, mat, coord_xtal

   use Casc_mod,           only: casc
   use Eval_mod,           only: eval10, eval26, eval27, eval30
   use Geom_CoordSys_mod,  only: ranrotCoordSys
   use Math_Intmath_mod,   only: index1to2, index2to1
   use Random_mod,         only: initRannum
   use Se_mod,             only: initSeModels, tabSe
   use Sn_FFP_mod,         only: initFFP, tabPmax
   use Sn_FFPcyl_mod,      only: initFFPcyl
   use Sn_RestFFP_mod,     only: initRestFFP
   use Sn_Scatter_mod,     only: initScatter
   use Target_mod,         only: initTarget
   use Target_Surface_mod, only: initSurface
   use Xtal_Lattice_mod,   only: initLattice
   use Xtal_OldSites_mod,  only: initOldSites
   use Xtal_Tvib_mod,      only: initFakeVib

#ifdef ENABLE_DRM
   use Cinterface_Swprotect_mod, only: initGTSLicense, shutdownGTSLicense
#endif
!
!  Default Type:
!  ============
!
   implicit none
!
!  Variables:
!  =========
!
   logical :: done
   integer(INT4) :: ichunk, iionchunk, i, iion, ithread
   integer(INT4) :: ichunk1, ichunk2, iionchunk1, iionchunk2, iionchunk_nion1, iionchunk_nion
!
!  Start of Executable Code:
!  ------------------------
!
!  Initializations
!
#ifdef ENABLE_DRM
   call initGTSLicense( )                       ! Make sure that the user has a license
!       ==============
#endif

   call readParams( )                           ! read input file
!       ==========

   call initParams( )                           ! initialize parameters in Param_*_mod
!       ==========

   call initFFP( )                              ! sets sn%pmax, which is needed in initTarget
!       =======                                 !    and initSeModels

   call initScatter( )                          ! may set sn%scrlen
!       ===========

   call initSeModels( )                         ! needs sn%pmax, may set se%enl
!       ============

   call initSurface( )                          ! needs sn%pmax, may modify geom%zlim etc.
!       ===========

   call initTarget( )                           ! needs sn%pmax(?)
!       ==========

   call writeParams( )                          ! write simulation parameters to output file
!       ===========

   call initLattice( )
!       ===========

   call initFakeVib( )
!       ===========

!$omp parallel
   call initRestFFP( )
!       ===========

   call initFFPcyl( )
!       ==========

   if( irxtal /= 0 )  call initOldSites( )
!                          ============

   call eval10( )                               ! initialize cascade evaluation
!       ======

   call initRannum( )                           ! must be called after eval10 because lbk=t
!       ==========                              !    may overwrite random number seed
!$omp end parallel

   call tabPmax( )                              ! maximum impact parameters
!       =======
   call tabSe( )                                ! nonlocal electronic stopping
!       =====
!
!  Determine initial and terminal indices of chunk and ions loop
!  Note that both nion and nion1 are multiples of nthread
!
   call index1to2( nion1/nthread+1, nionchunk/nthread, iionchunk_nion1, ichunk1 )
!       =========        i                  i                 o            o
   call index1to2( nion/nthread, nionchunk/nthread, iionchunk_nion, ichunk2 )
!       =========        i               i                 o           o

   Ion_Chunk_Loop: do ichunk = ichunk1, ichunk2          ! loop over chunks
!  --------------------------------------------

      if( ichunk == ichunk1 ) then
         iionchunk1 = iionchunk_nion1
      else
         iionchunk1 = 1
      endif
      if( ichunk == ichunk2 ) then
         iionchunk2 = iionchunk_nion
      else
         iionchunk2 = nionchunk/nthread
      end if

!$omp parallel private( ithread, iionchunk, i, iion, done )

      Ions_loop: do iionchunk = iionchunk1, iionchunk2   ! loop over ions within the chunk
!     ------------------------------------------------
!
!  Randomly rotate crystal for polycrystalline targets
!
         if( irxtal /= 0 ) then
            if( .not.ldyn .and. mat(irxtal)%lpoly ) then

               call ranrotCoordSys( coord_xtal )         ! threadprivate?
!                   ==============      io
            end if
         end if

         ithread = 1
!$       ithread = 1 + omp_get_thread_num( )
         call index2to1( ithread, iionchunk, nthread, i )
!             =========     i         i         i     o
         call index2to1( i, ichunk, nionchunk, iion )
!             =========  i    i         i       o
!
!  Simulate the cascade(s) from one impinging ion molecule
!
         call casc( iion, done )
!             ====   i     o
         if( done )  exit Ions_loop
!
!  Actions after each ion molecule or each n-th ion molecule (not for each ion atom)
!
         call eval26( )
!             ======
!
      end do Ions_loop
!     ----------------
!
!  Actions after each ion chunk
!
      call index2to1( iionchunk2*nthread, ichunk, nionchunk, iion )
!          =========          i              i        i       o
!$omp barrier
      call eval27( iion )
!          ======   i
!
!  Final evaluation of results
!
      if( ichunk == ichunk2 ) then
!$omp barrier
         call eval27( nion, final_=.true. )
!             ======   i         i
      end if
!$omp end parallel

   end do Ion_Chunk_loop
!  ---------------------
!
!  Output of results
!
!$omp parallel
   call eval30( nion )
!       ======   i
!$omp end parallel

!
!  Shut down the GTS licensing system
!
#ifdef ENABLE_DRM
   call shutdownGTSLicense( )
!       ==================
#endif
!
!=======================================================================
!
end program main
