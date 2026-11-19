module Se_ZBL_mod
!
!  Electronic stopping power according to ZBL (TRIM85)
!
!
!  Used Modules:
!  ============
!
   use Const_mod

   use Param_Setup_mod, only: nelem
   use Param_Elem_mod,  only: elem
   use Param_Ion_mod,   only: ion
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
   public :: initSeZBL              ! initialize parameters
!            =========
   public :: evalSeZBL              ! evaluate the energy loss
!            =========
!
!  Private Procedure:
!  =================
!
   private :: stpZBL                ! calculates the ZBL stopping power for 
!             ======                !   given energy, ion, and target
!
!  Private Data:
!  ============
!
   logical :: init = .false.

   integer(INT4), parameter :: NSE = 10000
   real(REAL8), allocatable :: etab(:)       ! (0:NSE)
   real(REAL8), allocatable :: setab(:,:,:)  ! (0:NSE,nelem,nelem)
!
!=======================================================================
!
contains
!
!=======================================================================
!
   subroutine initSeZBL () 
!             =========
!
!     Initializes the ZBL stopping power table
!
!  Arguments:
!
!     None
!
!  Local Variables:
!
      integer(INT4) :: ise, ielem1, ielem2
      real(REAL8) :: se
!
!  Start of Executable Code:
!  ------------------------
!
      init = .true.
!
!  Allocate memory for tables
!
      if( allocated( etab ) )  deallocate( etab )
      if( allocated( setab ) )  deallocate( setab )

      allocate( etab( 0:NSE ) )
      allocate( setab( 0:NSE, nelem, nelem ) )
!
!  Set table energies
!
      do ise=0,NSE
         etab(ise) = ise * ion%e / NSE
      end do
!
!  Initialize the stopping power table
!
      do ielem1 = 1, nelem
         do ielem2 = 1, nelem
            do ise = 0, NSE
               setab(ise,ielem2,ielem1) = ZERO
               call stpZBL( se, etab(ise), real(elem(ielem1)%z,REAL8), elem(ielem1)%m,  &
!                   ======  o      i                   i                      i
                            real(elem(ielem2)%z,REAL8) )
!                                       i
               setab(ise,ielem2,ielem1) = se
            end do
! debug
!               print '(t2,a,i2,a,i2,a,es12.3,a,e12.3)' ),
!     +            'z1=',nint(zproj),', ielem2=',ielem2,
!     +            ', e[eV]=',etab(ise),
!     +            ', se[eV/A]=',setab(ise,ielem2,ielem1)*dens(ielem2)
! end-debug
         end do
      end do

   end subroutine initSeZBL
!
!=======================================================================
!
   subroutine evalSeZBL( e, ielem1, ielem2, se )
!             =========  i    i        i     o
!
!     Evaluate the ZBL stopping power by interpolation in a table
!
!  Arguments:
!
      real(REAL8), intent(in) :: e              ! ion energy [eV]
      integer(INT4), intent(in) :: ielem1       ! index of the ion
      integer(INT4), intent(in) :: ielem2       ! index of the atom
      real(REAL8), intent(out) :: se            ! stopping power [eV A^2]
!
!  Local Variables:
!
      integer(INT4) :: ise
      real(REAL8) :: se1, se2
!
!  Start of Executable Code:
!  ------------------------
!
      if( .not.init )  goto 901

      if( e > etab(NSE) ) then                   ! calculate stopping power

         se = ZERO
         call stpZBL( se, e, real(elem(ielem1)%z,REAL8), elem(ielem1)%m,  &
!             ======  o   i              i                      i
                      real(elem(ielem2)%z,REAL8) )
!                                   i
      else                                      ! interpolate stopping power in table

         ise = min(NSE-1,int(e/etab(1)))
         se1 = setab(ise,ielem2,ielem1)
         se2 = setab(ise+1,ielem2,ielem1)
         se = se1 + (se2-se1)*(e-etab(ise))/(etab(ise+1)-etab(ise))

      end if

      return

 901  error stop 'evalSeZBL: not initialized'

   end subroutine evalSeZBL
!
!=======================================================================
!
   subroutine stpZBL( se, e, z1, m1, z2 )
!             ======  o   i  i   i   i
!
!     Calculates the electronic stopping power for a given energy according 
!     to Ziegler-Biersack-Littmark (TRIM85)
!
!  Arguments:
!
      real(REAL8), intent (out) :: se           ! electronic stopping power (eV A**2)
      real(REAL8), intent (in) :: e             ! ion energy (eV)
      real(REAL8), intent (in) :: z1            ! atomic number of the ion
      real(REAL8), intent (in) :: m1            ! relative atomic mass of the ion
      real(REAL8), intent (in) :: z2            ! atomic number of the target
!
!  Local Variables:
!
      real(REAL8), parameter :: HUN = 100
      real(REAL8), parameter :: P0P32 = 32/HUN
      real(REAL8), parameter :: P0P43 = 43/HUN
      real(REAL8), parameter :: P0P75 = 75/HUN
      real(REAL8), parameter :: P1P2 = 120/HUN
      real(REAL8), parameter :: P16 = 16
      real(REAL8), parameter :: P30 = 30

      integer(INT4) :: iz1, iz2

      real(REAL8) :: elog, sum, dgamma, seh0, sep0, gamma, q
      real(REAL8) :: scrlen, slmin, v, vr, yr, em, emmin, wrmin
      real(REAL8) :: seh, sep, sel0, sel, wmin, vmin, eb
!
! electronic stopping parameters according to ZBL (TRIM85)
!
      integer(INT4) :: i, ic
      real(REAL8), save :: slfac(92)            ! factor for screening length 
      real(REAL8), save :: vf(92)               ! Fermi velocity of the solid (in units of the
                                                !    Bohr velocity)
      real(REAL8), save :: cp(8,92)             ! coefficients for proton stopping formula

      data (vf(i),slfac(i),i=1,10) /  &
         1.0309 ,1.  ,  &
         0.15976,1.  ,  &
         0.59782,1.1 ,  &
         1.0781 ,1.06,  &
         1.0486 ,1.01,  &
         1.     ,1.03,  &
         1.058  ,1.04,  &
         0.93942,0.99,  &
         0.74562,0.95,  &
         0.3424 ,0.9  /
      data (vf(i),slfac(i),i=11,20) /  &
         0.45259,0.82,  &
         0.71074,0.81,  &
         0.90519,0.83,  &
         0.97411,0.84,  &
         0.97184,1.  ,  &
         0.89852,0.95,  &
         0.70827,0.97,  &
         0.39816,0.99,  &
         0.36552,0.98,  &
         0.62712,0.97 /
      data (vf(i),slfac(i),i=21,30) /  &
         0.81707,0.98,  &
         0.9943 ,0.97,  &
         1.1423 ,0.96,  &
         1.2381 ,0.93,  &
         1.1222 ,0.91,  &
         0.92705,0.9 ,  &
         1.0047 ,0.88,  &
         1.2    ,0.9 ,  &
         1.0661 ,0.9 ,  &
         0.97411,0.9  /
      data (vf(i),slfac(i),i=31,40) /  &
         0.84912,0.9 ,  &
         0.95   ,0.85,  &
         1.0903 ,0.9 ,  &
         1.0429 ,0.9 ,  &
         0.49715,0.91,  &
         0.37755,0.92,  &
         0.35211,0.9 ,  &
         0.57801,0.9 ,  &
         0.77773,0.9 ,  &
         1.0207 ,0.9  /
      data (vf(i),slfac(i),i=41,50) /  &
         1.029  ,0.9 ,  &
         1.2542 ,0.88,  &
         1.122  ,0.9 ,  &
         1.1241 ,0.88,  &
         1.0882 ,0.88,  &
         1.2709 ,0.9 ,  &
         1.2542 ,0.9 ,  &
         0.90094,0.88,  &
         0.74093,0.9 ,  &
         0.86054,0.9  /
      data (vf(i),slfac(i),i=51,60) /  &
         0.93155,0.9 ,  &
         1.0047 ,0.9 ,  &
         0.55379,0.96,  &
         0.43289,1.2 ,  &
         0.32636,0.9 ,  &
         0.5131 ,0.88,  &
         0.695  ,0.88,  &
         0.72591,0.85,  &
         0.71202,0.9 ,  &
         0.67413,0.9  /
      data (vf(i),slfac(i),i=61,70) /  &
         0.71418,0.92,  &
         0.71453,0.95,  &
         0.5911 ,0.99,  &
         0.70263,1.03,  &
         0.68049,1.05,  &
         0.68203,1.07,  &
         0.68121,1.08,  &
         0.68532,1.1 ,  &
         0.68715,1.08,  &
         0.61884,1.08 /
      data (vf(i),slfac(i),i=71,80) /  &
         0.71801,1.08,  &
         0.83048,1.08,  &
         1.1222 ,1.09,  &
         1.2381 ,1.09,  &
         1.045  ,1.1 ,  &
         1.0733 ,1.11,  &
         1.0953 ,1.12,  &
         1.2381 ,1.13,  &
         1.2879 ,1.14,  &
         0.78654,1.15 /
      data (vf(i),slfac(i),i=81,90) /  &
         0.66401,1.17,  &
         0.84912,1.2 ,  &
         0.88433,1.18,  &
         0.80746,1.17,  &
         0.43357,1.17,  &
         0.41923,1.16,  &
         0.43638,1.16,  &
         0.51464,1.16,  &
         0.73087,1.16,  &
         0.81065,1.16 /
      data (vf(i),slfac(i),i=91,92) /  &
         1.9578 ,1.16,  &
         1.0257 ,1.16 /

      data ((cp(ic,i),ic=1,8),i=1,10) /  &
        .0091827,.0053496,.69741  ,.48493,316.07,1.0143,9329.3,.053989,  &
         0.11393,.0051984,1.0822  ,.39252,1081. ,1.0645,4068.5,.017699,  &
         0.85837,.0050147,1.6044  ,.38844,1337.3,1.047 ,2659.2,.01898 ,  &
         0.8781 ,.0051049,5.4232  ,.2032 ,1200.6,1.0211,1401.8,.038529,  &
         1.4608 ,.0048836,2.338   ,.44249,1801.3,1.0352,1784.1,.02024 ,  &
         3.2579 ,.0049148,2.7156  ,.36473,2092.2,1.0291,2643.6,.018237,  &
         0.59674,.0050837,4.2073  ,.30612,2394.2,1.0255,4892.1,.016006,  &
         0.75253,.0050314,4.0824  ,.30067,2455.8,1.0181,5069.7,.017426,  &
         1.226  ,.0051385,3.2246  ,.32703,2525. ,1.0142,7563.6,.019469,  &
         1.0332 ,.0051645,3.004   ,.33889,2338.6,.99997,6991.2,.021799 /
      data ((cp(ic,i),ic=1,8),i=11,20) /  &
         6.0972 ,.0044292,3.1929  ,.45763,1363.3,.95182,2380.6,.081835,  &
         14.013 ,.0043646,2.2641  ,.36326,2187.4,.99098,6264.8,.0462  ,  &
         .039001,.0045415,5.5463  ,.39562,1589.2,.95316,816.16,.047484,  &
         2.178134 ,.0044454523 ,2.6045162  ,.6088463, 1550.2068,  &
         .9330245 ,1703.8459   ,.031619771 ,  &
         17.575 ,.0038346,0.078694,1.2388,2806. ,.97284,1037.6,.012879,  &
         16.126 ,.0038315,0.054164,1.3104,2813.3,.96587,1251.4,.011847,  &
         3.217  ,.0044579,3.6696  ,.5091 ,2734.6,.96253,2187.5,.016907,  &
         2.0379 ,.0044775,3.0743  ,.54773,3505. ,.97575,1714. ,.011701,  &
         0.74171,.0043051,1.1515  ,.95083,917.21,.8782 ,389.93,.18926 ,  &
         9.1316 ,.0043809,5.4611  ,.31327,3891.8,.97933,6267.9,.015196 /
      data ((cp(ic,i),ic=1,8),i=21,30) /  &
         7.2247 ,.0043718,6.1017  ,.37511,2829.2,.95218,6376.1,.020398,  &
         0.147  ,.0048456,6.3485  ,.41057,2164.1,.94028,5292.6,.050263,  &
         5.0611 ,.0039867,2.6174  ,.57957,2218.9,.92361,6323. ,.025669,  &
         0.53267,.0042968,0.39005 ,1.2725,1872.7,.90776,64.166,.030107,  &
         0.47697,.0043038,0.31452 ,1.3289,1920.5,.90649,45.576,.027469,  &
         .027426,.0035443,0.031563,2.1755,1919.5,.90099,23.902,.025363,  &
         0.16383,.0043042,0.073454,1.8592,1918.4,.89678,27.61 ,.023184,  &
         4.2562 ,.0043737,1.5606  ,.72067,1546.8,.87958,302.02,.040944,  &
         2.3508 ,.0043237,2.882   ,.50113,1837.7,.89992,2377. ,.04965 ,  &
         3.1095 ,.0038455,0.11477 ,1.5037,2184.7,.89309,67.306,.016588 /
      data ((cp(ic,i),ic=1,8),i=31,40) /  &
         15.322 ,.0040306,0.65391 ,.67668,3001.7,.92484,3344.2,.016366,  &
         3.6932 ,.0044813,8.608   ,.27638,2982.7,.9276 ,3166.6,.030874,  &
         7.1373 ,.0043134,9.4247  ,.27937,2725.8,.91597,3166.1,.025008,  &
         4.8979 ,.0042937,3.7793  ,.50004,2824.5,.91028,1282.4,.017061,  &
         1.3683 ,.0043024,2.5679  ,.60822,6907.8,.9817 ,628.01,.0068055,  &
         1.8301 ,.0042983,2.9057  ,.6038 ,4744.6,.94722,936.64,.0092242,  &
         0.42056,.0041169,0.01695 ,.3616 ,2252.7,.89192,39.752,.027757,  &
         30.78  ,.0037736,0.55813 ,.76816,7113.2,.97697,1604.4,.0065268,  &
         11.576 ,.0042119,7.0244  ,.37764,4713.5,.94264,2493.2,.01127 ,  &
         6.2406 ,.0041916,5.2701  ,.49453,4234.6,.93232,2063.9,.011844 /
      data ((cp(ic,i),ic=1,8),i=41,50) /  &
         0.33073,.0041243,1.7246  ,1.1062,1930.2,.86907,27.416,.038208,  &
         .017747,.0041715,0.14586 ,1.7305,1803.6,.86315,29.669,.032123,  &
         3.7229 ,.0041768,4.6286  ,.56769,1678. ,.86202,3094. ,.06244 ,  &
         0.13998,.0041329,0.25573 ,1.4241,1919.3,.86326,72.797,.032235,  &
         0.2859 ,.0041386,0.31301 ,1.3424,1954.8,.86175,115.18,.029342,  &
         0.76002,.0042179,3.386   ,.76285,1867.4,.85805,69.994,.036448,  &
         6.3957 ,.0041935,5.4689  ,.41378,1712.6,.85397,18493.,.056471,  &
         3.4717 ,.0041344,3.2337  ,.63788,1116.4,.81959,4766. ,.1179  ,  &
         2.5265 ,.0042282,4.532   ,.53562,1030.8,.81652,16252.,.19722 ,  &
         7.3683 ,.0041007,4.6791  ,.51428,1160. ,.82454,17965.,.13316  /
      data ((cp(ic,i),ic=1,8),i=51,60) /  &
         7.7197 ,.004388 ,3.242   ,.68434,1428.1,.83398,1786.7,.066512,  &
         16.78  ,.0041918,9.3198  ,.29568,3370.9,.90289,7431.7,.02616 ,  &
         4.2132 ,.0042098,4.6753  ,.57945,3503.9,.89261,1468.9,.014359,  &
         4.0818 ,.004214 ,4.4425  ,.58393,3945.3,.90281,1340.5,.013414,  &
         0.18517,.0036215,.00058788,3.5315,2931.3,.88936,26.18,.026393,  &
         4.8248 ,.0041458,6.0934  ,.57026,2300.1,.86359,2980.7,.038679,  &
         0.49857,.0041054,1.9775  ,.95877,786.55,.78509,806.6 ,.40882 ,  &
         3.2754 ,.0042177,5.768   ,.54054,6631.3,.94282,744.07,.0083026,  &
         2.9978 ,.0040901,4.5299  ,.62025,2161.2,.85669,1268.6,.043031,  &
         2.8701 ,.004096 ,4.2568  ,.6138 ,2130.4,.85235,1704.1,.039385 /
      data ((cp(ic,i),ic=1,8),i=61,70) /  &
         10.853 ,.0041149,5.8907  ,.46834,2857.2,.8755 ,3654.2,.029955,  &
         3.6407 ,.0041782,4.8742  ,.57861,1267.7,.82211,3508.2,.24174 ,  &
         17.645 ,.0040992,6.5855  ,.32734,3931.3,.90754,5156.7,.036278,  &
         7.5309 ,.0040814,4.9389  ,.50679,2519.7,.85819,3314.6,.030514,  &
         5.4742 ,.0040829,4.897   ,.51113,2340.1,.85296,2342.7,.035662,  &
         4.2661 ,.0040667,4.5032  ,.55257,2076.4,.84151,1666.6,.040801,  &
         6.8313 ,.0040486,4.3987  ,.51675,2003. ,.83437,1410.4,.03478 ,  &
         1.2707 ,.0040553,4.6295  ,.57428,1626.3,.81858,995.68,.055319,  &
         5.7561 ,.0040491,4.357   ,.52496,2207.3,.83796,1579.5,.027165,  &
         14.127 ,.0040596,5.8304  ,.37755,3645.9,.87823,3411.8,.016392 /
      data ((cp(ic,i),ic=1,8),i=71,80) /  &
         6.6948 ,.0040603,4.9361  ,.47961,2719. ,.85249,1885.8,.019713,  &
         3.0619 ,.0040511,3.5803  ,.59082,2346.1,.83713,1222. ,.020072,  &
         10.811 ,.0033008,1.3767  ,.76512,2003.7,.82269,1110.6,.024958,  &
         2.7101 ,.0040961,1.2289  ,.98598,1232.4,.79066,155.42,.047294,  &
         0.52345,.0040244,1.4038  ,.8551 ,1461.4,.79677,503.34,.036789,  &
         0.4616 ,.0040203,1.3014  ,.87043,1473.5,.79687,443.09,.036301,  &
         0.97814,.0040374,2.0127  ,.7225 ,1890.8,.81747,930.7 ,.02769 ,  &
         3.2086 ,.004051 ,3.6658  ,.53618,3091.2,.85602,1508.1,.015401,  &
         2.0035 ,.0040431,7.4882  ,.3561 ,4464.3,.88836,3966.5,.012839,  &
         15.43  ,.0039432,1.1237  ,.70703,4595.7,.88437,1576.5,.0088534/
      data ((cp(ic,i),ic=1,8),i=81,90) /  &
         3.1512 ,.0040524,4.0996  ,.5425 ,3246.3,.85772,1691.8,.015058,  &
         7.1896 ,.0040588,8.6927  ,.35842,4760.6,.88833,2888.3,.011029,  &
         9.3209 ,.004054 ,11.543  ,.32027,4866.2,.89124,3213.4,.011935,  &
         29.242 ,.0036195,0.16864 ,1.1226,5688. ,.89812,1033.3,.0071303,  &
         1.8522 ,.0039973,3.1556  ,.65096,3755. ,.86383,1602. ,.012042,  &
         3.222  ,.0040041,5.9024  ,.52678,4040.2,.86804,1658.4,.011747,  &
         9.3412 ,.0039661,7.921   ,.42977,5180.9,.88773,2173.2,.0092007,  &
         36.183 ,.0036003,0.58341 ,.86747,6990.2,.91082,1417.1,.0062187,  &
         5.9284 ,.0039695,6.4082  ,.52122,4619.5,.88083,2323.5,.011627,  &
         5.2454 ,.0039744,6.7969  ,.48542,4586.3,.87794,2481.5,.011282 /
      data ((cp(ic,i),ic=1,8),i=91,92) /  &
         33.702 ,.0036901,0.47257 ,.89235,5295.7,.8893 ,2053.3,.0091908,  &
         2.7589 ,.0039806,3.2092  ,.66122,2505.4,.82863,2065.1,.022816 /
!
!  Start of Executable Code
!
      iz1 = nint(z1)
      iz2 = nint(z2)

      if( e < 0.001 ) then
         se = ZERO
         return
      endif
! 
!  Energy per ion mass in keV
!
      em = e/(1000*m1)
! 
!  Calculate minimum energy per ion mass (emmin) of the high energy domain
!
      select case (iz1)

         case (1)

            emmin = ZERO

         case (2)

            emmin = ONE

         case (3:92)

            wrmin = max( p0p75, max( one, 0.13*z1**(2./3.) ) / vf(iz2))
            if( wrmin >= 1.2 ) then
               wmin = 0.5*(wrmin + sqrt( wrmin**2-0.8 ))
            else
               wmin = sqrt( max( ZERO, 5.-sqrt( 40.-20.*wrmin ) ) )
            endif
            vmin = wmin*vf(iz2)
            emmin = 25.*vmin**2

         case default

            error stop 'atomic number out of range in stpZBL'

      end select
!
!  If em is less than emmin, first calculate the stopping power at emmin
!  in order to determine the factor in the power law. Therefore replace 
!  em by emmin (em will be recovered later)
!
      em = max( em, emmin )
! 
!  Equivalent proton stopping according to Ziegler
!
      eb = 25.

      if( em >= eb ) then
                                                ! Ziegler's formula
         sel = cp(1,iz2)*em**cp(2,iz2) + cp(3,iz2)*em**cp(4,iz2)
         seh = cp(5,iz2)/em**cp(6,iz2) * log( cp(7,iz2)/em + cp(8,iz2)*em )
         sep = sel*seh/(sel+seh)

      else
                                                ! use power law for em<25 keV
         sel0 = cp(1,iz2)*eb**cp(2,iz2) + cp(3,iz2)*eb**cp(4,iz2)
         seh0 = cp(5,iz2)/eb**cp(6,iz2) * log( cp(7,iz2)/eb + cp(8,iz2)*eb )
         sep0 = sel0*seh0/(sel0+seh0)
         if( iz2 <= 6 ) then
            sep = sep0*(em/eb)**0.25
         else
            sep = sep0*(em/eb)**0.45
         end if

      end if

      sep = sep*1.e+1                           ! convert units 1e15 eV*cm^2 -> eV*A^2
!
!  Relative effective charge for high energies
!
      select case (iz1)

         case (1)                               ! H

            gamma = ONE

         case (2)                               ! He

            elog = log(em)
            sum = 0.2865 + 0.1266*elog - 0.001429*elog**2 +  &
                  0.02402*elog**3 - 0.01135*elog**4 + 0.001475*elog**5
            sum = min( p30, sum )
            gamma = sqrt( ONE - exp(-sum) )
            dgamma = gamma * (0.007+0.00005*z2) * exp( -(7.6-max( zero, elog ) )**2)
            gamma = gamma+dgamma

         case default                           ! "heavy" ions

            v = sqrt( em/25. )
                                                ! mean relative velocity between ion and
            if( v >= vf(iz2) ) then             !    target electrons
               vr = v * (ONE + 0.2*(vf(iz2)/v)**2)
            else
               vr = vf(iz2) * (0.75 + 0.5*(v/vf(iz2))**2 - 0.05*(v/vf(iz2))**4)
            end if
            yr = vr/z1**(2./3.)
                                                ! fractional charge state
            sum = - 0.803*yr**0.3 + 1.3167*yr**0.6 + 0.38157*yr + 0.008983*yr**2
            sum = min( p30, sum )
            q = ONE-exp(-sum)
            q = min( ONE, max( ZERO, q ) )
                                                ! screening length
            scrlen = (0.8 - q*min( p1p2, 0.6+z1/30. )) / z1**(1./3.)
            slmin = min( p0p43, max( p0p32, 0.12+z1/40. ) ) / z1**(1./3.)
            if( q < 0.2 ) then
               slmin = ZERO
            else if( q < 0.9-z1/40. ) then
               slmin = slmin * (q-0.2) / abs( 0.7-z1/40. )
            else if( q >= ONE-min( p16, z1 )/40. ) then
               slmin = slmin * (ONE-q) / (min( p16, z1 )/40.)
            end if
            scrlen = max( slmin, scrlen*slfac(iz1) )
                                                ! fractional effective charge
            gamma = q + (one-q) * 0.5/vf(iz2)**2 * log( one + (2.0844*scrlen*vf(iz2))**2 )
            dgamma = gamma * (0.18+0.0015*z2)/z1**2 * exp( -(7.6-max( zero, log(em) ) )**2)
            gamma = gamma+dgamma

      end select

      se = sep * (gamma*z1)**2
! 
!  Low energy stopping
!
      em = e/(1000.*m1)                         ! recover em 

                                                ! in the case of em<emmin, se now contains the
                                                !    stopping power at emmin. Apply power law
      if( em < emmin ) then                      !    for em<emmin
         if( (iz1 >= 3 .and. iz2 == 6) .or.  &
             ((iz1 >= 3 .and. iz1 <= 19) .and. (iz2 == 14 .or. iz2 == 32)) ) then
            se = se * (em/emmin)**0.375
         else
            se = se * sqrt( em/emmin )
         end if
      end if
 
   end subroutine stpZBL
!
!=======================================================================
!
end module Se_ZBL_mod
