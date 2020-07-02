include 'mkl_vsl.f90'
program main
USE MKL_VSL_TYPE
USE MKL_VSL
use omp_lib
implicit none

!***********************************************************************
!               Math and Phys cte [SI units]:
!***********************************************************************
double precision   , parameter :: pi     = 3.141592653589793d0
double precision   , parameter :: pi2    = 2.d0*pi
double precision   , parameter :: amu    = 1.66053886d-27   ![Kg]
double precision   , parameter :: c      = 2.99792458d8     ! vitesse de la lumière
double precision   , parameter :: hbar   = 1.055d-34
double precision   , parameter :: hplanck= 6.62d-34
double precision   , parameter :: qe     = 1.60217646d-19   ![C = sA] Elementary Charge
double precision   , parameter :: kb     = 1.3806503d-23    ![m2 kg s-2] Boltzman cte
double precision   , parameter :: ke     = 8.987551787d9    ![N m2 C-2] Coulomb cte =  1 / (4*pi*eps0)

!***********************************************************************
!                   Particles parameters:
!***********************************************************************
integer, parameter                   :: n_sp   = 2
!dec$     if defined(HCI0)
    !dir$ define HCI            = 0
!dec$ elseif defined(HCI1)
    !dir$ define HCI            = 1
!dec$ endif

! 0 = cloud preparation from scratch; 1 = HC ion injection; 2 = Post Ion Cloud evolution; 3 = Cloud Preparation from file; 4 = Inject Ion Build Stats
!dec$     if defined(Simu_type0)
    !dir$ define Simu_type      = 0
!dec$ elseif defined(Simu_type1)
    !dir$ define Simu_type      = 1
!dec$ elseif defined(Simu_type2)
    !dir$ define Simu_type      = 2
!dec$ elseif defined(Simu_type3)
    !dir$ define Simu_type      = 3
!dec$ elseif defined(Simu_type4)
    !dir$ define Simu_type      = 4
!dec$ elseif defined(Simu_type6)
    !dir$ define Simu_type      = 6
!dec$ endif


! 0 = RF ; 1 = RF + Gauss Axial; 2 = PS
!dir$ define potential_type = 1

integer, parameter :: n_ions1 = 1024
!dir$ if ( HCI .eq. 0 )
integer, parameter, dimension(n_sp)  :: n_ions = (/ n_ions1 , 0 /)
!dir$ else
integer, parameter, dimension(n_sp)  :: n_ions = (/ n_ions1 , 1/) ! Une GMol grand max !!!
!dir$ endif
integer, parameter, dimension(n_sp+1):: n_cut  = (/ n_ions(1), n_ions(2), sum(n_ions) /)
double precision   , parameter, dimension(n_sp)  :: charge =  qe*(/  1,      1/)
double precision   , parameter, dimension(n_sp)  :: mass   = amu*(/ 40, 1000000/) ! (\ ion, GMol \)

!***********************************************************************
!                   Simulations parameters:
!***********************************************************************
integer, parameter :: i_free__fly = 100
integer, parameter :: i_therm_fly = 2000
integer, parameter :: i_laser_fly = 5000
integer, parameter :: i_laser_fly_dt = 50000
integer, parameter :: i_cool__fly = 2000
integer, parameter :: i_cool__fly_dt = 50000
integer, parameter :: i_cool_heat = 0
double precision   , parameter, dimension(n_sp) :: T_bath = (/ 1.d-2, 1.0d-2 /) ! kelvin
double precision  , parameter :: beta = 4.d-24 / mass(1)

character(len=100), parameter :: str_extra = 'RFG'

!dec$ if defined(nE00)
double precision  , parameter :: frac_E0_Ecr = 0.0
character(len=100), parameter :: str_extra2 = ''
!dec$ elseif defined(nE01)
double precision  , parameter :: frac_E0_Ecr = 1.0d3
character(len=100), parameter :: str_extra2 = '1keV'
!dec$ elseif defined(nE02)
double precision  , parameter :: frac_E0_Ecr = 8.0d2
character(len=100), parameter :: str_extra2 = '800eV'
!dec$ elseif defined(nE03)
double precision  , parameter :: frac_E0_Ecr = 6.0d2    
character(len=100), parameter :: str_extra2 = '600eV'
!dec$ elseif defined(nE04)
double precision  , parameter :: frac_E0_Ecr = 4.0d2
character(len=100), parameter :: str_extra2 = '400eV'
!dec$ elseif defined(nE05)
double precision  , parameter :: frac_E0_Ecr = 2.0d2
character(len=100), parameter :: str_extra2 = '200eV'
!dec$ elseif defined(nE06)
double precision  , parameter :: frac_E0_Ecr = 1.0d2
character(len=100), parameter :: str_extra2 = '100eV'
!dec$ elseif defined(nE07)
double precision  , parameter :: frac_E0_Ecr = 5.0d1
character(len=100), parameter :: str_extra2 = '50eV'
!dec$ elseif defined(nE08)
double precision  , parameter :: frac_E0_Ecr = 1.0d1
character(len=100), parameter :: str_extra2 = '10eV'
!dec$ elseif defined(nE09)
double precision  , parameter :: frac_E0_Ecr = 8.0d0
character(len=100), parameter :: str_extra2 = '8eV'
!dec$ elseif defined(nE10)
double precision  , parameter :: frac_E0_Ecr = 6.0d0
character(len=100), parameter :: str_extra2 = '6eV'
!dec$ elseif defined(nE11)
double precision  , parameter :: frac_E0_Ecr = 4.0d0
character(len=100), parameter :: str_extra2 = '4eV'
!dec$ elseif defined(nE12)
double precision  , parameter :: frac_E0_Ecr = 2.0d0
character(len=100), parameter :: str_extra2 = '2eV'
!dec$ elseif defined(nE13)
double precision  , parameter :: frac_E0_Ecr = 1.0d0
character(len=100), parameter :: str_extra2 = '1eV'
!dec$ endif

!*************************************************
!~ integer, parameter :: i_stats0 = 1
!~ integer, parameter :: i_stats  = 1
!~ double precision   :: save_info(i_stats,13)
!*************************************************
! Modified by Jofre on 09/06/2020
character(len=32) :: arg
integer :: i_stats_a
integer :: i_stats_b
double precision :: save_info(13)
!*************************************************

character(len=150) :: str_file_to_load

!***********************************************************************
!                   Trap parameters:
!***********************************************************************
double precision   , parameter :: V_st  = 0.00d0
double precision   , parameter :: Omega = pi2 * 2.05d6

!~ double precision   , parameter :: V_rf  = 43.5d0 !q=0.4
!~ double precision   , parameter :: V_rf  = 64.0d0 !q=0.6

! double precision   , parameter :: V_rf  = 10 !q=0.1


double precision   , parameter :: ReqL     = 1.5312d-1   ! [mm]
double precision   , parameter :: Udcstart = 0.5d0   ! [V]
!dec$ if defined(Udc0)
double precision   , parameter :: Udc   = Udcstart    ! if other than 1V, then simply multiply by Udc
!dec$ elseif defined(Udc1)
double precision   , parameter :: Udc   = Udcstart*2
!dec$ elseif defined(Udc2)
double precision   , parameter :: Udc   = Udcstart*3
!dec$ elseif defined(Udc3)
double precision   , parameter :: Udc   = Udcstart*4
!dec$ elseif defined(Udc4)
double precision   , parameter :: Udc   = Udcstart*5
!dec$ elseif defined(Udc5)
double precision   , parameter :: Udc   = Udcstart*6
!dec$ elseif defined(Udc6)
double precision   , parameter :: Udc   = Udcstart*7
!dec$ elseif defined(Udc7)
double precision   , parameter :: Udc   = Udcstart*8
!dec$ elseif defined(Udc8)
double precision   , parameter :: Udc   = Udcstart*9
!dec$ elseif defined(Udc9)
double precision   , parameter :: Udc   = Udcstart*10
!dec$ elseif defined(Udc10)
double precision   , parameter :: Udc   = Udcstart*12
!dec$ elseif defined(Udc11)
double precision   , parameter :: Udc   = Udcstart*14
!dec$ elseif defined(Udc12)
double precision   , parameter :: Udc   = Udcstart*15
!dec$ elseif defined(Udc13)
double precision   , parameter :: Udc   = 3.33
!dec$ elseif defined(Udc14)
double precision   , parameter :: Udc   = 3.66
!dec$ elseif defined(Udc15)
double precision   , parameter :: Udc   = 4.5
!dec$ elseif defined(Udc16)
double precision   , parameter :: Udc   = 5.33
!dec$ elseif defined(Udc17)
double precision   , parameter :: Udc   = 5.66
!dec$ elseif defined(Udc18)
double precision   , parameter :: Udc   = 6
!dec$ elseif defined(Udc19)
double precision   , parameter :: Udc   = 6.5
!dec$ endif


!dec$ if defined(Vrf0)
double precision  , parameter :: V_rf  = 1.077d1
character(len=100), parameter :: str_extra3 = 'q=0.1'
!dec$ elseif defined(Vrf1)
double precision  , parameter :: V_rf  = 1.615d1
character(len=100), parameter :: str_extra3 = 'q=0.15'
!dec$ elseif defined(Vrf2)
double precision  , parameter :: V_rf  = 2.154d1
character(len=100), parameter :: str_extra3 = 'q=0.2'
!dec$ elseif defined(Vrf3)
double precision  , parameter :: V_rf  = 2.692d1
character(len=100), parameter :: str_extra3 = ''
!dec$ elseif defined(Vrf4)
double precision  , parameter :: V_rf  = 3.231d1
character(len=100), parameter :: str_extra3 = 'q=0.3'
!dec$ elseif defined(Vrf5)
double precision  , parameter :: V_rf  = 3.769d1
character(len=100), parameter :: str_extra3 = ''
!dec$ elseif defined(Vrf6)
double precision  , parameter :: V_rf  = 4.308d1
character(len=100), parameter :: str_extra3 = 'q=0.4'
!dec$ elseif defined(Vrf7)
double precision  , parameter :: V_rf  = 4.846d1
character(len=100), parameter :: str_extra3 = ''
!dec$ elseif defined(Vrf8)
double precision  , parameter :: V_rf  = 5.385d1
character(len=100), parameter :: str_extra3 = 'q=0.5'
!dec$ elseif defined(Vrf9)
double precision  , parameter :: V_rf  = 5.923d1
character(len=100), parameter :: str_extra3 = ''
!dec$ elseif defined(Vrf10)
double precision  , parameter :: V_rf  = 6.461d1
character(len=100), parameter :: str_extra3 = 'q=0.6'
!dec$ elseif defined(Vrf11)
double precision  , parameter :: V_rf  = 7.000d1
character(len=100), parameter :: str_extra3 = 'q=0.65'
!dec$ elseif defined(Vrf12)
double precision  , parameter :: V_rf  = 7.538d1
character(len=100), parameter :: str_extra3 = 'q=0.7'
!dec$ elseif defined(Vrf13)
double precision  , parameter :: V_rf  = 8.077d1
character(len=100), parameter :: str_extra3 = 'q=0.75'
!dec$ elseif defined(Vrf14)
double precision  , parameter :: V_rf  = 6.102d1
character(len=100), parameter :: str_extra3 = 'q=0.5666...'
!dec$ elseif defined(Vrf15)
double precision  , parameter :: V_rf  = 6.282d1
character(len=100), parameter :: str_extra3 = 'q=0.58333'
!dec$ elseif defined(Vrf16)
double precision  , parameter :: V_rf  = 6.641d1
character(len=100), parameter :: str_extra3 = 'q=0.6166'
!dec$ elseif defined(Vrf17)
double precision  , parameter :: V_rf  = 6.820d1
character(len=100), parameter :: str_extra3 = 'q=0.633'
!dec$ elseif defined(Vrf18)
double precision  , parameter :: V_rf  = 7.179d1
character(len=100), parameter :: str_extra3 = 'q=0.666'
!dec$ elseif defined(Vrf19)
double precision  , parameter :: V_rf  = 5.564d1
character(len=100), parameter :: str_extra3 = 'q=0.5166'
!dec$ elseif defined(Vrf20)
double precision  , parameter :: V_rf  = 5.743d1
character(len=100), parameter :: str_extra3 = 'q=0.533'
!dec$ endif

! Giant Mole Dimension:
double precision   , parameter :: r_0   = 2.865d-3/1.14511d0
double precision   , parameter :: d_0   = 4d-3 ! longueur du piege

! From fit to GiantMole Dimension when 1V applied using Gauss like function
double precision   , parameter :: L     = 0.0140608827209d0
double precision   , parameter :: wz_LC = pi2*90806.9982303 !when 1V applied
double precision   , parameter :: kappa = 0.23 ! facteur d'ecrantage
! double precision   , parameter :: omegaz = sqrt(2*kappa*qe*Udc/(mass(0)*d0**2)) ! wz calculé pour chaque Udc

! double precision   , parameter :: Udc   = 2.5d0  ! if other than 1V, then simply multiply by Udc

double precision   , parameter :: E_critic = 1.0d0

! qx= 0.6; alpha=0.2
!~ double precision   , parameter :: V_rf  = 61.374743856251285d0
!~ double precision   , parameter :: wz_LC = pi2*141758.770906d0
!~ double precision   , parameter :: E_critic = 79.80713651047097

! qx= 0.4; alpha=0.2
!~ double precision   , parameter :: V_rf  = 40.916495904167526d0
!~ double precision   , parameter :: wz_LC = pi2*94505.8472709d0 ! qx=0.4

!~ ! qx= 0.2; alpha=0.2
!~ double precision   , parameter :: V_rf  = 20.46d0
!~ double precision   , parameter :: wz_LC = pi2*47253.0d0 ! qx=0.2

double precision   , parameter, dimension(n_sp) :: wz2   = Udc*(/ wz_LC**2, wz_LC**2 * charge(2)*mass(1) / (charge(1)*mass(2)) /)
double precision   , parameter :: z0_HC = 1.5d-3
!~ double precision   , parameter :: z0_HC = 10.0d-3
!~ double precision :: z0_HC

!***********************************************************************
!                   Integration constants:
!***********************************************************************
integer, parameter :: ni   = 6 !number of threads
integer, parameter :: n_dt = 100 ! 100 pts par période RF
double precision   , parameter :: dt0  = pi2 / (n_dt*Omega) ! not normalised
!~ double precision   , parameter :: dt   = pi2 / n_dt ! normalised
double precision :: dt, dt1, dt2

!~ double precision   , parameter :: dt   = pi2 / (n_dt*Omega) ! not normalised
!~ double precision   , parameter :: dt1  = 0.5d0*dt
!~ double precision   , parameter :: dt2  = 0.5d0*dt*dt


integer, parameter :: n_save_temp = 1*n_dt ! Number of points in one RF to save temp

!~ ! Temperature cte:
!dir$ if ( HCI .eq. 0 )
double precision   , parameter, dimension(n_sp) :: m_kb_x_inv_n_ions  = mass/(kb*real(n_ions(1))*n_dt**2   )
double precision   , parameter, dimension(n_sp) :: m_kb_x_inv_n_ions2 = mass/(kb*real(n_ions(1))**2*n_dt**2)
double precision   , parameter, dimension(n_sp) :: T_bath2 = T_bath(1) / m_kb_x_inv_n_ions(1) ! kelvin
!dir$ else
double precision   , parameter, dimension(n_sp) :: m_kb_x_inv_n_ions  = mass/(kb*real(n_ions)*n_dt**2   )
double precision   , parameter, dimension(n_sp) :: m_kb_x_inv_n_ions2 = mass/(kb*real(n_ions)**2*n_dt**2)
double precision   , parameter, dimension(n_sp) :: T_bath2 = T_bath / m_kb_x_inv_n_ions ! kelvin
!dir$ endif
!~ double precision   , dimension(n_sp) :: m_kb_x_inv_n_ions, m_kb_x_inv_n_ions2, T_bath2


! Coulomb interaction cte:
! Normalised:
!~ double precision   , parameter, dimension(4) :: alpha = (ke/(r_0**3 * Omega**2))*        &
!~                                             (/ charge(1)*charge(1) / mass(1) &
!~                                              , charge(2)*charge(2) / mass(2) &
!~                                              , charge(2)*charge(1) / mass(1) &
!~                                              , charge(2)*charge(1) / mass(2) /)
! Not- Normalised:
double precision   , parameter, dimension(4) :: alpha = ke* &
                                            (/ charge(1)*charge(1) / mass(1) &
                                             , charge(2)*charge(2) / mass(2) &
                                             , charge(2)*charge(1) / mass(1) &
                                             , charge(2)*charge(1) / mass(2) /)

!***********************************************************************
!                       Laser parameters
!***********************************************************************
!dec$ if defined(d0)
    double precision, parameter :: detuning      = -1.0d0
!dir$ endif

    double precision, parameter :: saturation    =  1.0d0  ! saturation = Omega^2/(2*Gamma^2)
!********************* Laser Cooling of Ca+ ****************************
    double precision, parameter :: lambda        = 397.0d-9  !wavelenght
    double precision, parameter :: trans_width_w = 1/7.07e-9 !transition width (Hz)
!********************* Laser Cooling of Be+ ****************************
!~     double precision, parameter :: lambda        = 313.13292d-9  !wavelenght
!~     double precision, parameter :: trans_width_w = pi2*17.96d6  !transition width (Hz)
!***********************************************************************
    integer                     :: fluo_count
    integer                     :: state(n_cut(1))  ! Store the ground/excited state of the atom
!***********************************************************************

! Potential cte:
! potential_type = 0 => Linear Quadrupole RadioFrequency
! potential_type = 1 => Linear Quadrupole RadioFrequency Radially, Gauss like axially

! Normalised
!~ double precision   , parameter :: a_cte_Quad_RF_LC(4) = (/ &
!~     -2*charge(1)*V_st / (mass(1)*r_0**2 * Omega**2) + wz2(1)/(2*Omega**2), &
!~     +2*charge(1)*V_st / (mass(1)*r_0**2 * Omega**2) + wz2(1)/(2*Omega**2), &
!~      2*charge(1)*V_rf / (mass(1)*r_0**2 * Omega**2)                   , &
!~      -wz2(1) / Omega**2                                               /)
!~ double precision   , parameter :: a_cte_Quad_RF_HC(4) = (/ &
!~     -2*charge(2)*V_st / (mass(2)*r_0**2 * Omega**2) + wz2(2)/(2*Omega**2), &
!~     +2*charge(2)*V_st / (mass(2)*r_0**2 * Omega**2) + wz2(2)/(2*Omega**2), &
!~      2*charge(2)*V_rf / (mass(2)*r_0**2 * Omega**2)                   , &
!~      -wz2(2) / Omega**2                                               /)
! Not-Normalised
double precision   , parameter :: a_cte_Quad_RF_LC(4) = (/ &
    -2*charge(1)*V_st / (mass(1)*r_0**2 ) + 0.5*wz2(1),    &
    +2*charge(1)*V_st / (mass(1)*r_0**2 ) + 0.5*wz2(1),    &
     2*charge(1)*V_rf / (mass(1)*r_0**2 )             ,    &
     -wz2(1)                                             /)
double precision   , parameter :: a_cte_Quad_RF_HC(4) = (/ &
    -2*charge(2)*V_st / (mass(2)*r_0**2 ) + wz2(2)/2,      &
    +2*charge(2)*V_st / (mass(2)*r_0**2 ) + wz2(2)/2,      &
     2*charge(2)*V_rf / (mass(2)*r_0**2 )           ,      &
     -wz2(2)/)

double precision, parameter :: cte_ideal_quad_PS_LC(3) = (/-2*charge(1)**2 * V_rf**2 / (mass(1)**2*r_0**4*Omega**2) - 2*charge(1)*V_st / (mass(1)*r_0**2) + 0.5*wz2(1), &
                                                           -2*charge(1)**2 * V_rf**2 / (mass(1)**2*r_0**4*Omega**2) + 2*charge(1)*V_st / (mass(1)*r_0**2) + 0.5*wz2(1), &
                                                           -wz2(1) /)

double precision, parameter :: cte_ideal_quad_PS_HC(3) = (/-2*charge(2)**2 * V_rf**2 / (mass(2)**2*r_0**4*Omega**2) - 2*charge(2)*V_st / (mass(2)*r_0**2) + 0.5*wz2(2), &
                                                           -2*charge(2)**2 * V_rf**2 / (mass(2)**2*r_0**4*Omega**2) + 2*charge(2)*V_st / (mass(2)*r_0**2) + 0.5*wz2(2), &
                                                           -wz2(2) /)

double precision   , parameter :: period = 1/Omega
! Global variables:
double precision   , dimension(n_cut(3))   ::  r_x, r_y, r_z
double precision   , dimension(n_cut(3))   ::  v_x, v_y, v_z
double precision   , dimension(n_cut(3))   :: a1_x,a1_y,a1_z
double precision   , dimension(n_cut(3))   :: a2_x,a2_y,a2_z
double precision   , dimension(n_cut(3),3) :: v_rf_avg
double precision   , dimension(3)          :: T_CM__LC,T_aux_LC
!dir$ if ( HCI .eq. 1 )
double precision   , dimension(3)          :: T_CM__HC,T_aux_HC
!dir$ endif

integer :: j_start     = 1
integer :: j_save_temp = 0
integer :: i_save_temp = 0
integer :: total_ion_lost0 = 1
integer :: j_end, j_max, jj, iRF(3),  total_ion_lost

double precision, allocatable, dimension(:,:) ::  save_temperature
double precision, allocatable, dimension(:,:) ::  save_trj

integer, dimension(ni,n_sp+1)   :: ia , ib
double precision                :: time00, t_act, t_next_periode, rz_0, dz, time0_HCI
logical                         :: flag_dt0 = .false.

character(len=150)      :: str_file_Temp, str_file_xva, str_file_trj,str_file_PM
! Variables for the random generator of MKL
TYPE (VSL_STREAM_STATE) :: stream
integer(kind=4)         :: errcode
integer                 :: brng, method ,seed0, i

    brng   = VSL_BRNG_MCG31
    method = VSL_RNG_METHOD_UNIFORM_STD
    errcode= vslnewstream(stream, brng, seed0 )


    print*, 'dt0', dt0

    time00 = omp_get_wtime()
!~     call init_random_seed()
    CALL RANDOM_SEED ! Processor initializes the seed randomly from the date and time
    call distribute_ions()
!dir$ if(Simu_type .ne. 4)
    call create_files_names( )
!dir$ endif


!dir$     if(potential_type .eq. 0)
    print*, ''
    print*, '***********************************************************'
    print*, '         Using Linear Quadrupole RF potential'
    print*, '***********************************************************'
!dir$ elseif(potential_type .eq. 1)
    print*, ''
    print*, '***********************************************************'
    print*, '         Using Linear Quadrupole RF potential'
    print*, '         and Gaussian like for axial'
    print*, '***********************************************************'
    print*, 'Nions', n_ions1
    print*, 'Udc', Udc
    print*, 'Vrf', V_rf
    print*, str_extra3
    print*, 'EGmol', str_extra2
!dir$ elseif(potential_type .eq. 2)
    print*, ''
    print*, '***********************************************************'
    print*, '         Using Linear Quadrupole PS potential'
    print*, '***********************************************************'
!dir$ endif



!dir$ if(Simu_type .eq. 0)
    print*, ''
    print*, '***********************************************************'
    print*, '         Generating initial Ion Cloud'
    print*, '***********************************************************'
    !dir$  if ( HCI .ne. 0 )
    print*, '  - Compiling directive HCI not zero!'
    print*, '  - Stopping code'
    stop
    !dir$ endif
    if (n_cut(2).ne. 0) then
        print*, '  - HC ion number not zero!'
        print*, '  - Stopping code'
        stop
    endif
    call prepare_cold_ion_cloud()
!dir$ endif

!dir$ if(Simu_type .eq. 1)
    print*, ''
    print*, '***********************************************************'
    print*, '         Injection of the Extra ion'
    print*, '***********************************************************'
    !dir$  if ( HCI .eq. 0 )
    print*, '  - Compiling directive HCI not different from zero!'
    print*, '  - Stopping code'
    stop
    !dir$ endif
    if (n_cut(2).ne. 1) then
        print*, '  - HC ion number not one!'
        print*, '  - Stopping code'
        stop
    endif
    call inject_ion_simu()
!dir$ endif

!dir$ if(Simu_type .eq. 3)
    print*, ''
    print*, '***********************************************************'
    print*, '         Generating initial Ion Cloud from file'
    print*, '***********************************************************'
    !dir$  if ( HCI .ne. 0 )
    print*, '  - Compiling directive HCI not zero!'
    print*, '  - Stopping code'
    stop
    !dir$ endif
    if (n_cut(2).ne. 0) then
        print*, '  - HC ion number not zero!'
        print*, '  - Stopping code'
        stop
    endif
    call prepare_cold_ion_cloud_from_file()
!dir$ endif

!dir$ if(Simu_type .eq. 4)
    print*, ''
    print*, '***********************************************************'
    print*, '         Injection of the Extra ion - Build stats'
    print*, '***********************************************************'
    !dir$  if ( HCI .eq. 0 )
    print*, '  - Compiling directive HCI not different from zero!'
    print*, '  - Stopping code'
    stop
    !dir$ endif
    if (n_cut(2).ne. 1) then
        print*, '  - HC ion number not one!'
        print*, '  - Stopping code'
        stop
    endif

!*************************************************
! Modified by Jofre on 09/06/2020
!*************************************************

! Get i_stats_a, i_stats_b from the command line:
! Syntax is: ./a.out -ia i -ib i'
! Where i is the integer value desired
    call get_command_argument(2, arg); read(arg , *) i_stats_a
    call get_command_argument(4, arg); read(arg , *) i_stats_b
    do i = i_stats_a, i_stats_b
!*************************************************
        print*, 'i stats', i
        
        call create_files_names( )
        call inject_ions_HC()

        save_info(1) = r_x(n_cut(3))
        save_info(2) = r_y(n_cut(3))
        save_info(3) = r_z(n_cut(3))
        save_info(4) = v_x(n_cut(3))
        save_info(5) = v_y(n_cut(3))
        save_info(6) = v_z(n_cut(3))
        print*,'Before vz', save_info(3), save_info(6)
        
        call inject_ion_simu_laser()

        save_info(7) = r_x(n_cut(3))
        save_info(8) = r_y(n_cut(3))
        save_info(9) = r_z(n_cut(3))
        save_info(10) = v_x(n_cut(3))
        save_info(11) = v_y(n_cut(3))
        save_info(12) = v_z(n_cut(3))
        save_info(13) = t_act - time0_HCI
        print*,'After vz', save_info(9), save_info(12)
        print*, 'Time of crossing [us]:', (t_act - time0_HCI)*1.0e6

    open(15, file=trim(adjustl(str_file_xva))//'.dat', status='replace', access='sequential', action='write')
        write(15,226) save_info
    close(15)
    enddo

226     format( 13(1X,e27.19e3))
!~     close(15)
!dir$ endif

!dir$ if(Simu_type .eq. 2)
    print*, ''
    print*, '***********************************************************'
    print*, '         Post-Injection Evolution of Ion Cloud'
    print*, '***********************************************************'
    !dir$  if ( HCI .ne. 0 )
    print*, '  - Compiling directive HCI HCI not zero!'
    print*, '  - Stopping code'
    stop
    !dir$ endif
    if (n_cut(2).ne. 0) then
        print*, '  - HC ion number not zero!'
        print*, '  - Stopping code'
        stop
    endif

    call get_command_argument(2, arg); read(arg , *) i_stats_a
    call get_command_argument(4, arg); read(arg , *) i_stats_b
    do i = i_stats_a, i_stats_b
        print*, 'i stats', i
        call create_files_names( )
        call post_inject_ion_simu()
    enddo

!dir$ endif

!dir$ if(Simu_type .eq. 5)
    print*, ''
    print*, '***********************************************************'
    print*, '         Test initial HCI position'
    print*, '***********************************************************'
    !dir$  if ( HCI .eq. 0 )
    print*, '  - Compiling directive HCI not different from zero!'
    print*, '  - Stopping code'
    stop
    !dir$ endif
    if (n_cut(2).ne. 1) then
        print*, '  - HC ion number not one!'
        print*, '  - Stopping code'
        stop
    endif


    call init_from_file()

    dz = 1.0d-4
    do i = 0, 200
!~         print*, 'i stats', i
        z0_HC = L - i*dz
        call inject_ions_HC()
        call compute_optimal_dt_only_HC()
            print*, i, z0_HC, dt
        if (dt < dt0) then
            stop
        endif
    enddo


!dir$ endif

!dir$ if(Simu_type .eq. 6)
    print*, ''
    print*, '***********************************************************'
    print*, '         Laserless Evolution of Ion Cloud'
    print*, '***********************************************************'
    !dir$  if ( HCI .ne. 0 )
    print*, '  - Compiling directive HCI HCI not zero!'
    print*, '  - Stopping code'
    stop
    !dir$ endif
    if (n_cut(2).ne. 0) then
        print*, '  - HC ion number not zero!'
        print*, '  - Stopping code'
        stop
    endif

    do i = i_stats_a, i_stats_b
        print*, 'i stats', i
        call create_files_names( )
        call post_cooling_ion_simu()
    enddo

!dir$ endif

    print*, 'Total time', omp_get_wtime() - time00
contains

subroutine prepare_cold_ion_cloud()
implicit none
    print*, 'Initialization of the ions from scratch...'

    call initialize_ions_LC()
    state       = 0
    jj          = 0
    t_act       = 0
    total_ion_lost0 = 1

    iRF(1) = 0; iRF(2) = floor(n_dt/2.0); iRF(3) = 0;

    j_max = (i_free__fly + i_therm_fly + i_cool__fly + i_cool_heat + i_laser_fly)*n_dt
    allocate(save_temperature(floor(j_max / real(n_save_temp)) + 2, 0:7))

    ! zero step:
    call update_accelerations()
    a1_x = a2_x; a1_y = a2_y; a1_z = a2_z;

    dt = dt0; call update_dt_parameters()

    j_end = i_free__fly*n_dt; call free__fly()
    j_end = i_therm_fly*n_dt; call therm_fly()

    j_end = i_laser_fly*n_dt - modulo(j_start + i_laser_fly*n_dt -1,n_dt) + int(0.25*n_dt)
    call laser_fly()
!~     call cooling_fly()

    call save_data()
endsubroutine

subroutine prepare_cold_ion_cloud_from_file()
implicit none
    print*, 'Loading'

    call init_from_file()
    dt = dt0; call update_dt_parameters()

    j_max = i_cool__fly*n_dt
    j_max = j_max + (n_dt - modulo(j_max, n_dt))

    allocate(save_temperature(floor(j_max / real(n_save_temp)) + 2,0:7))

    total_ion_lost0 = 1

    j_end = i_cool__fly*n_dt;
    call cooling_fly()
    print*, 'Saving data';
    call save_data()

endsubroutine

subroutine inject_ion_simu_cool()
implicit none
    print*, 'Loading'
    j_save_temp = 0
    j_start     = 0

    call init_from_file()

    j_max = i_cool__fly_dt*n_dt
    j_max = j_max + (n_dt - modulo(j_max, n_dt))
    allocate(save_temperature(floor(j_max / real(n_save_temp)) + 2,0:13))
    allocate(save_trj        (floor(j_max / real(n_save_temp)) + 2,  10))
    total_ion_lost0 = 1

    j_end = i_cool__fly_dt*n_dt;
    call cooling_fly_var_dt()

    print*, 'Saving data'; call save_data()
    print*, 'Saving trj' ; call save_trj_to_file()

endsubroutine

subroutine inject_ion_simu_laser()
implicit none
    print*, 'Loading'
    j_save_temp = 0
    j_start     = 0

    call init_from_file()
    time0_HCI = t_act

    j_max = i_laser_fly_dt*n_dt
    j_max = j_max + (n_dt - modulo(j_max, n_dt))
    allocate(save_temperature(floor(j_max / real(n_save_temp)) + 2,0:13))
    allocate(save_trj        (floor(j_max / real(n_save_temp)) + 2,  10))
    total_ion_lost0 = 1

    j_end = i_laser_fly_dt*n_dt;
    call laser_fly_var_dt()

    print*, 'Saving data'; call save_data()
    print*, 'Saving trj' ; call save_trj_to_file()

endsubroutine

subroutine post_inject_ion_simu()
implicit none
    print*, 'Loading file...'
    call init_from_file_post_injec() 

    j_max           = i_laser_fly*n_dt

    allocate(save_temperature(floor(j_max / real(n_save_temp)) + 2, 0:7))

    dt = dt0; call update_dt_parameters()

    j_end = i_laser_fly*n_dt - modulo(j_start + i_laser_fly*n_dt -1,n_dt) + int(0.25*n_dt)
    call laser_fly()

    print*, 'Ions losts:', total_ion_lost0 -1

    call save_data()
endsubroutine

subroutine post_cooling_ion_simu()
implicit none
    print*, 'Loading file...'
    call init_from_file_post_cooling() 

    j_max           = i_laser_fly*n_dt

    allocate(save_temperature(floor(j_max / real(n_save_temp)) + 2, 0:7))

    dt = dt0; call update_dt_parameters()

    j_end = i_laser_fly*n_dt - modulo(j_start + i_laser_fly*n_dt -1,n_dt) + int(0.25*n_dt)
    call free__fly()

    print*, 'Ions losts:', total_ion_lost0 -1

    call save_data()
endsubroutine

subroutine free__fly()
implicit none
logical          :: flag_timing = .true.
double precision :: time1, time2
    print*, 'Start Free Fly'
    do jj = j_start , j_start + j_end -1
        call update_positions()
        call update_accelerations()
        call update_velocities()
        call Measure_of_Temperature()
        call Check_save_temp()

        if (flag_timing) then
            if (jj.eq.2) then
                time1 = omp_get_wtime()
            endif
            if (jj.eq.22) then
                time2 = omp_get_wtime()
                flag_timing = .False.
                print*, 'Estimated time[s]', (time2 - time1)/20       * j_max
                print*, 'Estimated time[m]', (time2 - time1)/20/60    * j_max
                print*, 'Estimated time[h]', (time2 - time1)/20/60/60 * j_max
            endif
        endif
    enddo
    j_start = jj;
    print*, 'Finished Free Fly'
endsubroutine

subroutine therm_fly()
implicit none
    print*, ''
    print*, 'Start Thermal Fly'
    do jj = j_start , j_start + j_end -1
        call update_positions()
        call update_accelerations()
        call update_velocities()

        call Measure_of_Temperature_and_rescaling()
        call Check_save_temp()
    enddo
    j_start = jj;
    print*, 'Finished Thermal Fly'
endsubroutine

subroutine laser_fly()
implicit none
    print*, ''
    print*, 'Start Laser Fly'

    do jj = j_start , j_start -1 + j_end
        call update_positions()
!dir$ if(Simu_type .eq. 2)
        call check_ions_lost() ! Only in post injection simulation
!dir$ endif
        call update_accelerations()
        call update_velocities()
!~         if (t_act .gt. 5.0d-3) then
            call laser_interaction_two_lasers_angle( )
!~         endif
        
        call Measure_of_Temperature()
        call Check_save_temp()
    enddo
    j_start = jj;
    print*, 'Finished Cool Fly'
endsubroutine

subroutine cooling_fly()
implicit none
logical          :: flag_timing = .true.
double precision :: time1, time2

    print*, ''
    print*, 'Start Cooling Fly'

    do jj = j_start , j_start -1 + j_end
        call update_positions()
        call update_accelerations()
        call a_cooling()
        call update_velocities()

        call Measure_of_Temperature()
        call Check_save_temp()


        if (flag_timing) then
            if (jj.eq.j_start+10) then
                time1 = omp_get_wtime()
            endif
            if (jj.eq.j_start+110) then
                time2 = omp_get_wtime()
                flag_timing = .False.
                print*, 'Estimated time[s]', (time2 - time1)/100       * j_max
                print*, 'Estimated time[m]', (time2 - time1)/100/60    * j_max
                print*, 'Estimated time[h]', (time2 - time1)/100/60/60 * j_max
            endif
        endif

    enddo
    j_start = jj;
    print*, 'Finished Cool Fly'
endsubroutine

subroutine cooling_fly_var_dt()
implicit none
    print*, ''
    print*, 'Start Laser Fly  with var dt'
!~     dt = dt0; call update_dt_parameters()

    do jj = j_start , j_start + j_end -1
        call compute_optimal_dt_all( )
        call update_positions()
        call update_accelerations()
        if (flag_dt0) call a_cooling()
        call update_velocities()


        if (flag_dt0) then
            call Measure_of_Temperature()
            call Check_save_temp_trj()
            flag_dt0 = .false.
        endif

        if (r_z(n_cut(3)) > rz_0) then
            print*, 'The other side has been reached by the HC particle'
            goto 113
        elseif (v_z(n_cut(3))<0)then
            print*, 'The particle is turning back!'
            goto 113
        endif
    enddo
113 continue
    j_start = jj;
    print*, 'Finished Cool Fly'
endsubroutine

subroutine laser_fly_var_dt()
implicit none
    print*, ''
    print*, 'Start Laser Fly  with var dt'
!~     dt = dt0; call update_dt_parameters()

    do jj = j_start , j_start + j_end -1
        call compute_optimal_dt_all( )
        call update_positions()
        call update_accelerations()
        call update_velocities()

        if (flag_dt0) then
            call laser_interaction_two_lasers_angle( )
            call Measure_of_Temperature()
            call Check_save_temp_trj()
            flag_dt0 = .false.
        endif

        if (r_z(n_cut(3)) > rz_0) then
            print*, 'The other side has been reached by the HC particle'
            goto 113
        elseif (v_z(n_cut(3))<0)then
            print*, 'The particle is turning back!'
            goto 113
        endif
    enddo
113 continue
    j_start = jj;
    print*, 'Finished Cool Fly'
endsubroutine

subroutine update_positions()
implicit none
integer :: im
!-----------------------------------------------------------------------
    t_act = t_act + dt
!   One time step Verlet-Velocity
    !$omp parallel default(none) private(im) shared (r_x,r_y,r_z,v_x,v_y,v_z,a1_x,a1_y,a1_z, a2_x,a2_y,a2_z) firstprivate(ia,ib, dt, dt2)
        im = omp_get_thread_num()+1
        r_x(ia(im,3): ib(im,3)) = r_x(ia(im,3): ib(im,3)) + v_x(ia(im,3): ib(im,3))*dt + a1_x(ia(im,3): ib(im,3))*dt2;
        r_y(ia(im,3): ib(im,3)) = r_y(ia(im,3): ib(im,3)) + v_y(ia(im,3): ib(im,3))*dt + a1_y(ia(im,3): ib(im,3))*dt2;
        r_z(ia(im,3): ib(im,3)) = r_z(ia(im,3): ib(im,3)) + v_z(ia(im,3): ib(im,3))*dt + a1_z(ia(im,3): ib(im,3))*dt2;

        a2_x(ia(im,3):ib(im,3)) = 0.0d0
        a2_y(ia(im,3):ib(im,3)) = 0.0d0
        a2_z(ia(im,3):ib(im,3)) = 0.0d0
    !$omp end parallel
end subroutine

subroutine update_accelerations()
implicit none
    call a_Coulomb_1sp_2sp()

    !dir$ if     ( potential_type .eq. 0 )
        call a_trap_ideal_quad_RF()
    !dec$ elseif ( potential_type .eq. 1 )
        call a_trap_ideal_quad_RF_Gauss()
    !dec$ elseif ( potential_type .eq. 2 )
        call a_trap_ideal_quad_PS()
    !dec$ endif

endsubroutine

subroutine update_velocities()
implicit none
integer :: im
    !$omp parallel default(none) private(im) shared (v_x,v_y,v_z,a1_x,a1_y,a1_z,a2_x,a2_y,a2_z) firstprivate(ia,ib,dt1)
    im = omp_get_thread_num()+1
    v_x(ia(im,3):ib(im,3)) = v_x(ia(im,3):ib(im,3)) + dt1*(a1_x(ia(im,3):ib(im,3)) + a2_x(ia(im,3):ib(im,3))); a1_x(ia(im,3):ib(im,3)) = a2_x(ia(im,3):ib(im,3));
    v_y(ia(im,3):ib(im,3)) = v_y(ia(im,3):ib(im,3)) + dt1*(a1_y(ia(im,3):ib(im,3)) + a2_y(ia(im,3):ib(im,3))); a1_y(ia(im,3):ib(im,3)) = a2_y(ia(im,3):ib(im,3));
    v_z(ia(im,3):ib(im,3)) = v_z(ia(im,3):ib(im,3)) + dt1*(a1_z(ia(im,3):ib(im,3)) + a2_z(ia(im,3):ib(im,3))); a1_z(ia(im,3):ib(im,3)) = a2_z(ia(im,3):ib(im,3));
    !$omp end parallel
end subroutine

subroutine check_ions_lost()
implicit none
integer :: i,j,i_temp,j_lost
integer :: ion_lost_position(n_ions(1)), j_store(n_ions(1))
double precision :: r_temp
!-----------------------------------------------------------------------
    j_store = 0
    where (abs(r_z(total_ion_lost0:)).gt.(0.5d0*L))
        j_store(total_ion_lost0:) =  1
    elsewhere
        where (r_x(total_ion_lost0:)*r_x(total_ion_lost0:) + r_y(total_ion_lost0:)*r_y(total_ion_lost0:) .gt. r_0*r_0)
            ! Ion lost radially
            j_store(total_ion_lost0:) =  1
        endwhere
    endwhere


    ! I need to shift the lost ions to the top. It can not be done in parallel:
    if (sum(j_store).ne.0) then
        j = 0
        do i = 1, n_ions(1)
            if (j_store(i).eq.1) then
                j = j + 1
                ion_lost_position(j) = i
            endif
        enddo

        j = 0
        total_ion_lost = total_ion_lost0 + sum(j_store)
        print*, 'IOn losing', total_ion_lost -1
        do i = total_ion_lost0, total_ion_lost-1
            j = j + 1
            j_lost = ion_lost_position(j)
            if (j_lost.ne.i) then
                r_temp = r_x(i); r_x(i) = r_x(j_lost) ; r_x(j_lost) = r_temp
                r_temp = r_y(i); r_y(i) = r_y(j_lost) ; r_y(j_lost) = r_temp
                r_temp = r_z(i); r_z(i) = r_z(j_lost) ; r_z(j_lost) = r_temp

                r_temp = v_x(i); v_x(i) = v_x(j_lost) ; v_x(j_lost) = r_temp
                r_temp = v_y(i); v_y(i) = v_y(j_lost) ; v_y(j_lost) = r_temp
                r_temp = v_z(i); v_z(i) = v_z(j_lost) ; v_z(j_lost) = r_temp

                r_temp = a1_x(i); a1_x(i) = a1_x(j_lost) ; a1_x(j_lost) = r_temp
                r_temp = a1_y(i); a1_y(i) = a1_y(j_lost) ; a1_y(j_lost) = r_temp
                r_temp = a1_z(i); a1_z(i) = a1_z(j_lost) ; a1_z(j_lost) = r_temp

                r_temp = a2_x(i); a2_x(i) = a2_x(j_lost) ; a2_x(j_lost) = r_temp
                r_temp = a2_y(i); a2_y(i) = a2_y(j_lost) ; a2_y(j_lost) = r_temp
                r_temp = a2_z(i); a2_z(i) = a2_z(j_lost) ; a2_z(j_lost) = r_temp

                i_temp = state(i); state(i) = state(j_lost) ; state(j_lost) = i_temp

                r_temp = v_rf_avg(i,1); v_rf_avg(i,1) = v_rf_avg(j_lost,1) ; v_rf_avg(j_lost,1) = r_temp
                r_temp = v_rf_avg(i,2); v_rf_avg(i,2) = v_rf_avg(j_lost,2) ; v_rf_avg(j_lost,2) = r_temp
                r_temp = v_rf_avg(i,3); v_rf_avg(i,3) = v_rf_avg(j_lost,3) ; v_rf_avg(j_lost,3) = r_temp
            endif
        enddo
        total_ion_lost0 = total_ion_lost
    endif
!~     print*, 'total_ion_lost0', total_ion_lost0
end subroutine

subroutine Measure_of_Temperature()
implicit none
integer :: i
    v_rf_avg(total_ion_lost0:,1) = v_rf_avg(total_ion_lost0:,1) + v_x(total_ion_lost0:);
    v_rf_avg(total_ion_lost0:,2) = v_rf_avg(total_ion_lost0:,2) + v_y(total_ion_lost0:);
    v_rf_avg(total_ion_lost0:,3) = v_rf_avg(total_ion_lost0:,3) + v_z(total_ion_lost0:);

    do i = 1,3
        iRF(i) = iRF(i) + 1;
        if (iRF(i) == n_dt) then
            iRF(i)           = 0;

            T_CM__LC(i)      = sum(v_rf_avg(total_ion_lost0:n_ions(1) ,i) ) **2
            T_aux_LC(i)      = sum(v_rf_avg(total_ion_lost0:n_ions(1) ,i)**2)

            !dir$ if ( HCI .eq. 1 )
            T_CM__HC(i)      = sum(v_rf_avg(n_ions(1)+1:,i))**2
            T_aux_HC(i)      = sum(v_rf_avg(n_ions(1)+1:,i)**2)
            !dec$ endif

            v_rf_avg(total_ion_lost0:,i) = 0.0d0
        endif
    enddo
end subroutine

subroutine Measure_of_Temperature_and_rescaling()
implicit none
integer :: i
    v_rf_avg(:,1) = v_rf_avg(:,1) + v_x;
    v_rf_avg(:,2) = v_rf_avg(:,2) + v_y;
    v_rf_avg(:,3) = v_rf_avg(:,3) + v_z;
    do i = 1,3
        iRF(i) = iRF(i) + 1;
        if (iRF(i) == n_dt) then
            iRF(i)           = 0;

            T_CM__LC(i)      = sum(v_rf_avg(1:n_ions(1) ,i))**2
            T_aux_LC(i)      = sum(v_rf_avg(1:n_ions(1) ,i)**2)
            select case(i)
            case(1)
                v_x(1:n_ions(1)) = v_x(1:n_ions(1))*dsqrt(T_bath2(1) / T_aux_LC(i));  ! Apply
            case(2)
                v_y(1:n_ions(1)) = v_y(1:n_ions(1))*dsqrt(T_bath2(1) / T_aux_LC(i));  ! Apply
            case(3)
                v_z(1:n_ions(1)) = v_z(1:n_ions(1))*dsqrt(T_bath2(1) / T_aux_LC(i));  ! Apply
            endselect

            !dir$ if ( HCI .eq. 1 )
            T_CM__HC(i)      = sum(v_rf_avg(n_ions(1)+1:,i))**2
            T_aux_HC(i)      = sum(v_rf_avg(n_ions(1)+1:,i)**2)
            select case(i)
            case(1)
                v_x(n_ions(1)+1:) = v_x(n_ions(1)+1:)*dsqrt(T_bath2(2) / T_aux_HC(i));    ! Apply
            case(2)
                v_y(n_ions(1)+1:) = v_y(n_ions(1)+1:)*dsqrt(T_bath2(2) / T_aux_HC(i));    ! Apply
            case(3)
                v_z(n_ions(1)+1:) = v_z(n_ions(1)+1:)*dsqrt(T_bath2(2) / T_aux_HC(i));    ! Apply
            endselect
            !dec$ endif

            v_rf_avg(:,i) = 0.0d0
        endif
    enddo
end subroutine

subroutine Check_save_temp()
implicit none
    if (n_save_temp == i_save_temp) then
        i_save_temp = 0;
        j_save_temp = j_save_temp + 1

        save_temperature(j_save_temp, 0 )   = t_act
        save_temperature(j_save_temp, 1 :3 ) = T_CM__LC
        save_temperature(j_save_temp, 4 :6 ) = T_aux_LC
        save_temperature(j_save_temp, 7    ) = real(fluo_count)
        fluo_count = 0

        !dir$ if ( HCI .eq. 1 )
        save_temperature(j_save_temp, 8 :10 ) = T_CM__HC
        save_temperature(j_save_temp, 11:13) = T_aux_HC
        !dec$ endif
    else
        i_save_temp = i_save_temp + 1;
    endif
endsubroutine

subroutine Check_save_temp_trj()
implicit none
    if (n_save_temp == i_save_temp) then
        i_save_temp = 0;
        j_save_temp = j_save_temp + 1

        save_temperature(j_save_temp, 0 )   = t_act
        save_temperature(j_save_temp, 1 :3 ) = T_CM__LC
        save_temperature(j_save_temp, 4 :6 ) = T_aux_LC
        save_temperature(j_save_temp, 7    ) = real(fluo_count)

        fluo_count = 0

        !dir$ if ( HCI .eq. 1 )
        save_temperature(j_save_temp, 8 :10 ) = T_CM__HC
        save_temperature(j_save_temp, 11:13) = T_aux_HC

        save_trj(j_save_temp,1) =  t_act

        save_trj(j_save_temp,2) =  r_x(n_cut(3))
        save_trj(j_save_temp,3) =  r_y(n_cut(3))
        save_trj(j_save_temp,4) =  r_z(n_cut(3))

        save_trj(j_save_temp,5) =  v_x(n_cut(3))
        save_trj(j_save_temp,6) =  v_y(n_cut(3))
        save_trj(j_save_temp,7) =  v_z(n_cut(3))

        save_trj(j_save_temp,8) =  a1_x(n_cut(3))
        save_trj(j_save_temp,9) =  a1_y(n_cut(3))
        save_trj(j_save_temp,10)=  a1_z(n_cut(3))

        !dec$ endif
    else
        i_save_temp = i_save_temp + 1;
    endif
endsubroutine

subroutine store_trj()
implicit none
    save_trj(j_save_temp,1) =  t_act

    save_trj(j_save_temp,2) =  r_x(n_cut(3))
    save_trj(j_save_temp,3) =  r_y(n_cut(3))
    save_trj(j_save_temp,4) =  r_z(n_cut(3))

    save_trj(j_save_temp,5) =  v_x(n_cut(3))
    save_trj(j_save_temp,6) =  v_y(n_cut(3))
    save_trj(j_save_temp,7) =  v_z(n_cut(3))

    save_trj(j_save_temp,8) =  a1_x(n_cut(3))
    save_trj(j_save_temp,9) =  a1_y(n_cut(3))
    save_trj(j_save_temp,10)=  a1_z(n_cut(3))

end subroutine

subroutine a_trap_ideal_quad_RF()
implicit none
integer :: im
double precision    :: cte_aux(2)
double precision    :: cos_jj
    cos_jj = dcos(t_act*Omega)
    cte_aux(1) = a_cte_Quad_RF_LC(1) + a_cte_Quad_RF_LC(3)*cos_jj
    cte_aux(2) = a_cte_Quad_RF_LC(2) - a_cte_Quad_RF_LC(3)*cos_jj
    !$omp parallel private(im)  &
    !$omp shared (r_x,r_y,r_z,a2_x,a2_y,a2_z) &
    !$omp firstprivate(ia,ib,cte_aux)
    im = omp_get_thread_num()+1
    a2_x(ia(im,1):ib(im,1)) = a2_x(ia(im,1):ib(im,1))          + cte_aux(1) * r_x(ia(im,1):ib(im,1))
    a2_y(ia(im,1):ib(im,1)) = a2_y(ia(im,1):ib(im,1))          + cte_aux(2) * r_y(ia(im,1):ib(im,1))
    a2_z(ia(im,1):ib(im,1)) = a2_z(ia(im,1):ib(im,1)) + a_cte_Quad_RF_LC(4) * r_z(ia(im,1):ib(im,1))
    !$omp end parallel

    !dir$ if ( HCI .eq. 1 )
    cte_aux(1) = a_cte_Quad_RF_HC(1) + a_cte_Quad_RF_HC(3)*cos_jj
    cte_aux(2) = a_cte_Quad_RF_HC(2) - a_cte_Quad_RF_HC(3)*cos_jj
    a2_x(n_cut(3)) = a2_x(n_cut(3)) + cte_aux(1) * r_x(n_cut(3))
    a2_y(n_cut(3)) = a2_y(n_cut(3)) + cte_aux(2) * r_y(n_cut(3))
    a2_z(n_cut(3)) = a2_z(n_cut(3)) + a_cte_Quad_RF_HC(4) * r_z(n_cut(3))
    !dec$ endif
end subroutine

subroutine a_trap_ideal_quad_PS()
implicit none
integer :: im
    !cte_ideal_quad_PS(:,1) = -2*charge**2 * Vrf**2 / (mass**2*r0**4*Omega**2) - 2*charge*Vdc / (mass*r0**2) + 0.5*wz2
    !cte_ideal_quad_PS(:,2) = -2*charge**2 * Vrf**2 / (mass**2*r0**4*Omega**2) + 2*charge*Vdc / (mass*r0**2) + 0.5*wz2
    !cte_ideal_quad_PS(:,3) = -wz2

    !$omp parallel private(im)  &
    !$omp shared (ia,ib, r_x,r_y,r_z,a1_x,a1_y,a1_z,a2_x,a2_y,a2_z)
        im = omp_get_thread_num()+1

        a2_x(ia(im,1):ib(im,1)) = a2_x(ia(im,1):ib(im,1)) + cte_ideal_quad_PS_LC(1) * r_x(ia(im,1):ib(im,1))
        a2_y(ia(im,1):ib(im,1)) = a2_y(ia(im,1):ib(im,1)) + cte_ideal_quad_PS_LC(2) * r_y(ia(im,1):ib(im,1))
        a2_z(ia(im,1):ib(im,1)) = a2_z(ia(im,1):ib(im,1)) + cte_ideal_quad_PS_LC(3) * r_z(ia(im,1):ib(im,1))
    !$omp end parallel

    !dir$ if ( HCI .eq. 1 )
    a2_x(n_cut(3)) = a2_x(n_cut(3)) + cte_ideal_quad_PS_HC(1) * r_x(n_cut(3))
    a2_y(n_cut(3)) = a2_y(n_cut(3)) + cte_ideal_quad_PS_HC(2) * r_y(n_cut(3))
    a2_z(n_cut(3)) = a2_z(n_cut(3)) + cte_ideal_quad_PS_HC(3) * r_z(n_cut(3))
!dir$ endif

end subroutine

subroutine a_trap_ideal_quad_RF_Gauss()
implicit none
integer :: im
double precision, parameter :: eta = 1.0d0 / (4.0d0*(dsqrt(10.0d0) + 5.0d0 ))
double precision, parameter :: cte_gauss = -1/(2*eta*L**2)
double precision    :: cte_aux(2)
double precision    :: cos_jj

    cos_jj = dcos(t_act*Omega)
    cte_aux(1) = a_cte_Quad_RF_LC(1) + a_cte_Quad_RF_LC(3)*cos_jj
    cte_aux(2) = a_cte_Quad_RF_LC(2) - a_cte_Quad_RF_LC(3)*cos_jj
    !$omp parallel private(im)  &
    !$omp shared (r_x,r_y,r_z,a2_x,a2_y,a2_z) &
    !$omp firstprivate(ia,ib,cte_aux)
    im = omp_get_thread_num()+1
    a2_x(ia(im,1):ib(im,1)) = a2_x(ia(im,1):ib(im,1)) + cte_aux(1) * r_x(ia(im,1):ib(im,1))
    a2_y(ia(im,1):ib(im,1)) = a2_y(ia(im,1):ib(im,1)) + cte_aux(2) * r_y(ia(im,1):ib(im,1))
    a2_z(ia(im,1):ib(im,1)) = a2_z(ia(im,1):ib(im,1)) -     wz2(1) * r_z(ia(im,1):ib(im,1))*exp(cte_gauss * r_z(ia(im,1):ib(im,1))*r_z(ia(im,1):ib(im,1)) )
    !$omp end parallel

    !dir$ if ( HCI .eq. 1 )
    cte_aux(1) = a_cte_Quad_RF_HC(1) + a_cte_Quad_RF_HC(3)*cos_jj
    cte_aux(2) = a_cte_Quad_RF_HC(2) - a_cte_Quad_RF_HC(3)*cos_jj
    a2_x(n_cut(3)) = a2_x(n_cut(3)) + cte_aux(1) * r_x(n_cut(3))
    a2_y(n_cut(3)) = a2_y(n_cut(3)) + cte_aux(2) * r_y(n_cut(3))
    a2_z(n_cut(3)) = a2_z(n_cut(3)) - wz2(2)*r_z(n_cut(3))*exp(cte_gauss* r_z(n_cut(3))*r_z(n_cut(3)))
    !dec$ endif

end subroutine

subroutine laser_interaction_two_lasers_angle( )
implicit none
!~ Assumptions:
!~ * g1 = g2 = 1
!~ * Laser uniform radially (larger than ion cloud)
! Direction of the wave vector (spherical coordinates x=r sin theta cos phi; )
!~ double precision, parameter :: theta_k    =  pi*0.25d0
!~ double precision, parameter :: phi_k      =  pi*0.25d0      !phi_k

double precision, parameter :: theta_k    =  0.01
double precision, parameter :: phi_k      =  pi/4      !phi_k


double precision, parameter :: costheta = cos(theta_k)      ! theta_k and phi_k => direction of the wave vector
double precision, parameter :: sintheta = sin(theta_k)
double precision, parameter :: cosphi   = cos(phi_k)
double precision, parameter :: sinphi   = sin(phi_k)

double precision, parameter :: A21dt         = trans_width_w*dt0                                 ! spontaneous and stimulated emission
double precision, parameter :: v_recul       = hplanck/lambda/mass(1)
double precision, parameter :: lambda2pi     = pi2/lambda
double precision, parameter :: B12dt_cte(2)  = (/ trans_width_w*saturation*dt0/(1+2*saturation), 16.0d0*pi**2/ (lambda**2*trans_width_w**2)*(1+2*saturation) /)
double precision, parameter :: angl_aux(3)   = (/ cosphi*sintheta, sinphi*sintheta, costheta/)
double precision, parameter :: v_recul_aux(3)=  v_recul*angl_aux

double precision, parameter :: detu_aux = detuning*lambda*trans_width_w/pi2
!~ double precision :: detu_aux

double precision, dimension(n_cut(1)) :: B12dt_1, B12dt_2, theta_em,phi_em,costheta_em,sintheta_em, fluo_count_aux, rand_em, rand_abs, rand_choice, rand_theta
integer :: im, choice(n_cut(1))

!~     detu_aux = detuning*lambda*trans_width_w/pi2

    errcode = vdRngUniform(method, stream, n_cut(1), rand_em    ,  0.0d0, 1.d0)
    errcode = vdRngUniform(method, stream, n_cut(1), rand_abs   ,  0.0d0, 1.d0)
    errcode = vdRngUniform(method, stream, n_cut(1), rand_choice,  0.0d0, 1.d0)

    errcode = vdRngUniform(method, stream, n_cut(1), rand_theta , -1.0d0, 1.d0)
    errcode = vdRngUniform(method, stream, n_cut(1), phi_em     ,  0.0d0, pi2 )

    !$omp parallel default(none) &
    !$omp private(im, B12dt_1, B12dt_2, theta_em, costheta_em,  sintheta_em, choice) &
    !$omp shared( v_x,v_y,v_z,state, ia, ib,fluo_count_aux, rand_em, rand_abs, rand_choice, rand_theta, phi_em)
        im = omp_get_thread_num()+1
        fluo_count_aux(ia(im, 1):ib(im, 1)) = 0

        B12dt_2(ia(im, 1):ib(im, 1)) = v_x( ia(im, 1):ib(im, 1))*angl_aux(1) + v_y( ia(im, 1):ib(im, 1))*angl_aux(2) + v_z( ia(im, 1):ib(im, 1))*angl_aux(3)

        B12dt_1(ia(im, 1):ib(im, 1)) = B12dt_cte(1)/(1.0d0+B12dt_cte(2)*( (detu_aux - B12dt_2(ia(im, 1):ib(im, 1))) )**2 )
        B12dt_2(ia(im, 1):ib(im, 1)) = B12dt_cte(1)/(1.0d0+B12dt_cte(2)*( (detu_aux + B12dt_2(ia(im, 1):ib(im, 1))) )**2 )

        choice(ia(im,1):ib(im,1)) = nint(rand_choice(ia(im,1):ib(im,1)))

        where (state(ia(im, 1):ib(im, 1)).eq.0)
            where (rand_abs(ia(im,1):ib(im,1)).le. ( B12dt_1(ia(im,1):ib(im,1))*choice(ia(im,1):ib(im,1)) + B12dt_2(ia(im,1):ib(im,1))*(1-choice(ia(im,1):ib(im,1)) )  ) )
                v_x(ia(im, 1):ib(im, 1))   = v_x(ia(im, 1):ib(im, 1)) + v_recul_aux(1)*(2*choice(ia(im,1):ib(im,1)) -1)
                v_y(ia(im, 1):ib(im, 1))   = v_y(ia(im, 1):ib(im, 1)) + v_recul_aux(2)*(2*choice(ia(im,1):ib(im,1)) -1)
                v_z(ia(im, 1):ib(im, 1))   = v_z(ia(im, 1):ib(im, 1)) + v_recul_aux(3)*(2*choice(ia(im,1):ib(im,1)) -1)
                state(ia(im, 1):ib(im, 1)) = 1
            elsewhere (rand_abs(ia(im,1):ib(im,1)).le.(B12dt_1(ia(im,1):ib(im,1))+B12dt_2(ia(im,1):ib(im,1))))
                v_x(ia(im, 1):ib(im, 1))   = v_x(ia(im, 1):ib(im, 1)) - v_recul_aux(1)*(2*choice(ia(im,1):ib(im,1)) -1)
                v_y(ia(im, 1):ib(im, 1))   = v_y(ia(im, 1):ib(im, 1)) - v_recul_aux(2)*(2*choice(ia(im,1):ib(im,1)) -1)
                v_z(ia(im, 1):ib(im, 1))   = v_z(ia(im, 1):ib(im, 1)) - v_recul_aux(3)*(2*choice(ia(im,1):ib(im,1)) -1)
                state(ia(im, 1):ib(im, 1)) = 1
            endwhere
        elsewhere
            ! Spontaneous emission
            where (rand_em(ia(im, 1):ib(im, 1)).le.A21dt)
                sintheta_em(ia(im, 1):ib(im, 1)) = v_recul*dsin(acos(rand_theta(ia(im, 1):ib(im, 1))))
                v_x      (ia(im, 1):ib(im, 1)) = v_x(ia(im, 1):ib(im, 1)) + sintheta_em(ia(im, 1):ib(im, 1))*dcos(phi_em(ia(im, 1):ib(im, 1)))
                v_y      (ia(im, 1):ib(im, 1)) = v_y(ia(im, 1):ib(im, 1)) + sintheta_em(ia(im, 1):ib(im, 1))*dsin(phi_em(ia(im, 1):ib(im, 1)))
                v_z      (ia(im, 1):ib(im, 1)) = v_z(ia(im, 1):ib(im, 1)) + v_recul*rand_theta(ia(im, 1):ib(im, 1))
                state    (ia(im, 1):ib(im, 1)) = 0
                fluo_count_aux(ia(im, 1):ib(im, 1)) = 1

            elsewhere (rand_em(ia(im,1):ib(im,1)).le.(A21dt + B12dt_1(ia(im,1):ib(im,1))*choice(ia(im,1):ib(im,1)) + B12dt_2(ia(im,1):ib(im,1))*(1-choice(ia(im,1):ib(im,1))) ) )
                v_x(ia(im, 1):ib(im, 1)) = v_x(ia(im, 1):ib(im, 1)) - v_recul_aux(1)*(2*choice(ia(im,1):ib(im,1)) -1)
                v_y(ia(im, 1):ib(im, 1)) = v_y(ia(im, 1):ib(im, 1)) - v_recul_aux(2)*(2*choice(ia(im,1):ib(im,1)) -1)
                v_z(ia(im, 1):ib(im, 1)) = v_z(ia(im, 1):ib(im, 1)) - v_recul_aux(3)*(2*choice(ia(im,1):ib(im,1)) -1)

                state    (ia(im,1):ib(im,1)) = 0
                fluo_count_aux(ia(im,1):ib(im,1)) =  1

            elsewhere (rand_em(ia(im,1):ib(im,1)).le.(A21dt + B12dt_1(ia(im,1):ib(im,1)) + B12dt_2(ia(im,1):ib(im,1)) ) )
                v_x(ia(im, 1):ib(im, 1))   = v_x(ia(im, 1):ib(im, 1)) + v_recul_aux(1)*(2*choice(ia(im,1):ib(im,1)) -1)
                v_y(ia(im, 1):ib(im, 1))   = v_y(ia(im, 1):ib(im, 1)) + v_recul_aux(2)*(2*choice(ia(im,1):ib(im,1)) -1)
                v_z(ia(im, 1):ib(im, 1))   = v_z(ia(im, 1):ib(im, 1)) + v_recul_aux(3)*(2*choice(ia(im,1):ib(im,1)) -1)

                state    (ia(im,1):ib(im,1)) = 0
                fluo_count_aux(ia(im,1):ib(im,1)) =  1

            endwhere
        endwhere
    !$omp end parallel

    fluo_count = fluo_count + sum(fluo_count_aux)
end subroutine

subroutine compute_optimal_dt_only_HC( )
implicit none
!---------------  Locals  ----------------------------------------------
integer :: i, im,j
double precision, dimension(3) :: rji, vji, aji
double precision :: tij_1, tij_2, dtij(ni), dtij_min, tij_store(n_cut(1))
double precision, parameter :: dt_factor = 2*5.0d-4
!-----------------------------------------------------------------------
    tij_1 = dt0
    !$omp parallel default(none) &
    !$omp private(im, i,j,rji,vji,aji, tij_2) &
    !$omp firstprivate(r_x,r_y,r_z, v_x, v_y, v_z, a1_x, a1_y, a1_z,n_cut, ia, ib,tij_1 ) &
    !$omp shared(dtij, tij_store)
    im = omp_get_thread_num() + 1
    do i = ia(im,1), ib(im,1)
        rji(1)  = r_x(n_cut(3)) - r_x(i)
        rji(2)  = r_y(n_cut(3)) - r_y(i)
        rji(3)  = r_z(n_cut(3)) - r_z(i)
        rji(1)  = dsqrt(rji(1)*rji(1) + rji(2)*rji(2) + rji(3)*rji(3))

        vji(1)  = v_x(n_cut(3)) - v_x(i)
        vji(2)  = v_y(n_cut(3)) - v_y(i)
        vji(3)  = v_z(n_cut(3)) - v_z(i)
        vji(2)  = vji(1)*vji(1) + vji(2)*vji(2) + vji(3)*vji(3)
        vji(1)  = dsqrt(vji(2))

        aji(1)  = a1_x(n_cut(3)) - a1_x(i)
        aji(2)  = a1_y(n_cut(3)) - a1_y(i)
        aji(3)  = a1_z(n_cut(3)) - a1_z(i)
        aji(1)  = dsqrt(aji(1)*aji(1) + aji(2)*aji(2) + aji(3)*aji(3))

        tij_2 = (- vji(1) + dsqrt(vji(2) + dt_factor*aji(1)*rji(1)) )/aji(1)

        if (tij_2 < tij_1) then
            tij_1 = tij_2
        endif
    enddo
    dtij(im) = tij_1
    !$omp end parallel


    dtij_min = minval(dtij)
    if (dtij_min < dt0) then
        dt = dtij_min
    else
        dt = dt0
    endif

    if (t_act + dt > t_next_periode) then
        dt = t_next_periode - t_act
        t_next_periode = t_next_periode + dt0
        flag_dt0 = .true.
    endif

    call update_dt_parameters()

end subroutine

subroutine compute_optimal_dt_all( )
implicit none
!---------------  Locals  ----------------------------------------------
integer :: i, im,j
double precision, dimension(3) :: rji, vji, aji
double precision :: tij_1, tij_2, dtij(ni), dtij_min, tij_store(n_cut(1))
double precision, parameter :: dt_factor = 2*5.0d-2

!-----------------------------------------------------------------------
    tij_1 = dt0
    !$omp parallel default(none) &
    !$omp private(im, i,j,rji,vji,aji, tij_2) &
    !$omp firstprivate(r_x,r_y,r_z, v_x, v_y, v_z, a1_x, a1_y, a1_z,n_cut, ia, ib,tij_1 ) &
    !$omp shared(dtij, tij_store)
    im = omp_get_thread_num() + 1
    do i = ia(im,3), ib(im,3)
        do j=1,i-1
            rji(1)  = r_x(j) - r_x(i)
            rji(2)  = r_y(j) - r_y(i)
            rji(3)  = r_z(j) - r_z(i)
            rji(1)  = dsqrt(rji(1)*rji(1) + rji(2)*rji(2) + rji(3)*rji(3))

            vji(1)  = v_x(j) - v_x(i)
            vji(2)  = v_y(j) - v_y(i)
            vji(3)  = v_z(j) - v_z(i)
            vji(2)  = vji(1)*vji(1) + vji(2)*vji(2) + vji(3)*vji(3)
            vji(1)  = dsqrt(vji(2))

            aji(1)  = a1_x(j) - a1_x(i)
            aji(2)  = a1_y(j) - a1_y(i)
            aji(3)  = a1_z(j) - a1_z(i)
            aji(1)  = dsqrt(aji(1)*aji(1) + aji(2)*aji(2) + aji(3)*aji(3))

            tij_2 = (- vji(1) + dsqrt(vji(2) + dt_factor*aji(1)*rji(1)) )/aji(1)
            if (tij_2 < tij_1) then
                tij_1 = tij_2
            endif
        enddo
        do j=i+1,n_cut(3)
            rji(1)  = r_x(j) - r_x(i)
            rji(2)  = r_y(j) - r_y(i)
            rji(3)  = r_z(j) - r_z(i)
            rji(1)  = dsqrt(rji(1)*rji(1) + rji(2)*rji(2) + rji(3)*rji(3))

            vji(1)  = v_x(j) - v_x(i)
            vji(2)  = v_y(j) - v_y(i)
            vji(3)  = v_z(j) - v_z(i)
            vji(2)  = vji(1)*vji(1) + vji(2)*vji(2) + vji(3)*vji(3)
            vji(1)  = dsqrt(vji(2))

            aji(1)  = a1_x(j) - a1_x(i)
            aji(2)  = a1_y(j) - a1_y(i)
            aji(3)  = a1_z(j) - a1_z(i)
            aji(1)  = dsqrt(aji(1)*aji(1) + aji(2)*aji(2) + aji(3)*aji(3))

            tij_2 = (- vji(1) + dsqrt(vji(2) + dt_factor*aji(1)*rji(1)) )/aji(1)
            if (tij_2 < tij_1) then
                tij_1 = tij_2
            endif
        enddo
    enddo
    dtij(im) = tij_1
    !$omp end parallel

    dtij_min = minval(dtij)
    if (dtij_min < dt0) then
        dt = dtij_min
    else
        dt = dt0
    endif
    call update_dt_parameters()

    if (t_act + dt > t_next_periode) then
        t_next_periode = t_next_periode + dt0
        flag_dt0 = .true.
    endif

end subroutine

subroutine update_dt_parameters()
implicit none
    dt1  = 0.5d0*dt
    dt2  = 0.5d0*dt*dt
endsubroutine

subroutine a_Coulomb_1sp_2sp( ) !Done
implicit none
!---------------  Locals  ----------------------------------------------
integer :: i, im,j
double precision, dimension(3) :: rji
double precision :: r2inv
double precision, parameter :: softening = 1.0d-20
!dir$ if ( HCI .eq. 1 )
double precision, dimension(ni) :: aux_x, aux_y, aux_z
!dec$ endif

!-----------------------------------------------------------------------
!dir$ if ( HCI .eq. 0 )
    !$omp parallel default(none) &
    !$omp private(im, i,j,rji,r2inv) &
    !$omp firstprivate(r_x,r_y,r_z, n_cut, ia, ib) &
    !$omp shared(a2_x, a2_y, a2_z)
    im = omp_get_thread_num() + 1
    do i = ia(im,1), ib(im,1)
        do j = 1, n_cut(1)
            rji(1)  = r_x(j) - r_x(i)
            rji(2)  = r_y(j) - r_y(i)
            rji(3)  = r_z(j) - r_z(i)
            r2inv   = 1.d0/dsqrt(rji(1)*rji(1) + rji(2)*rji(2) + rji(3)*rji(3) + softening)
            r2inv   = r2inv * r2inv * r2inv * alpha(1)
            a2_x(i) = a2_x(i) - rji(1)*r2inv
            a2_y(i) = a2_y(i) - rji(2)*r2inv
            a2_z(i) = a2_z(i) - rji(3)*r2inv
        enddo
    enddo
    !$omp end parallel
!dec$ endif

!dir$ if ( HCI .eq. 1 )
    !$omp parallel default(none) &
    !$omp private(im, i,j,rji,r2inv) &
    !$omp firstprivate(r_x,r_y,r_z, n_cut, ia, ib) &
    !$omp shared(a2_x, a2_y, a2_z)
    im = omp_get_thread_num() + 1

    do i = ia(im,1), ib(im,1)
        do j = 1, n_cut(1)
            rji(1)  = r_x(j) - r_x(i)
            rji(2)  = r_y(j) - r_y(i)
            rji(3)  = r_z(j) - r_z(i)
            r2inv   = 1.d0/dsqrt(rji(1)*rji(1) + rji(2)*rji(2) + rji(3)*rji(3)+ softening)
            r2inv   = r2inv * r2inv * r2inv * alpha(1)
            a2_x(i) = a2_x(i) - rji(1)*r2inv
            a2_y(i) = a2_y(i) - rji(2)*r2inv
            a2_z(i) = a2_z(i) - rji(3)*r2inv
        enddo

        rji(1) = r_x(n_cut(3)) - r_x(i)
        rji(2) = r_y(n_cut(3)) - r_y(i)
        rji(3) = r_z(n_cut(3)) - r_z(i)
        r2inv  = 1.d0/dsqrt(rji(1)*rji(1) + rji(2)*rji(2) + rji(3)*rji(3))
        r2inv  = r2inv * r2inv * r2inv * alpha(3)

        a2_x(i) = a2_x(i) - rji(1)*r2inv
        a2_y(i) = a2_y(i) - rji(2)*r2inv
        a2_z(i) = a2_z(i) - rji(3)*r2inv
    enddo
    !$omp end parallel

    !$omp parallel default(none) &
    !$omp private(im, j,rji,r2inv) &
    !$omp firstprivate(r_x,r_y,r_z, n_cut, ia, ib) &
    !$omp shared(aux_x, aux_y, aux_z)
    im = omp_get_thread_num() + 1
    aux_x(im) = 0.0d0; aux_y(im) = 0.0d0; aux_z(im) = 0.0d0
    do j = ia(im,1), ib(im,1)
        rji(1) = r_x(j) - r_x(n_cut(3))
        rji(2) = r_y(j) - r_y(n_cut(3))
        rji(3) = r_z(j) - r_z(n_cut(3))
        r2inv  = 1.d0/dsqrt(rji(1)*rji(1) + rji(2)*rji(2) + rji(3)*rji(3))
        r2inv  = r2inv * r2inv * r2inv * alpha(4)

        aux_x(im) = aux_x(im) - rji(1)*r2inv
        aux_y(im) = aux_y(im) - rji(2)*r2inv
        aux_z(im) = aux_z(im) - rji(3)*r2inv
    enddo
    !$omp end parallel
    a2_x(n_cut(3)) = sum(aux_x)
    a2_y(n_cut(3)) = sum(aux_y)
    a2_z(n_cut(3)) = sum(aux_z)
!dec$ endif


end subroutine

subroutine a_cooling()
implicit none
integer :: im
    !$omp parallel private(im) shared (v_x,v_y,v_z,a2_x,a2_y,a2_z,ia,ib)
    im = omp_get_thread_num()+1
!~     if (im==1) then
!~         print*, jj
!~         print*, '111',a2_x(ia(im,1):ib(im,1))
!~     endif
    a2_x(ia(im,1):ib(im,1)) = a2_x(ia(im,1):ib(im,1)) - beta * v_x(ia(im,1):ib(im,1))
    a2_y(ia(im,1):ib(im,1)) = a2_y(ia(im,1):ib(im,1)) - beta * v_y(ia(im,1):ib(im,1))
    a2_z(ia(im,1):ib(im,1)) = a2_z(ia(im,1):ib(im,1)) - beta * v_z(ia(im,1):ib(im,1))

!~     if (im==1) then
!~         print*, '222',a2_x(ia(im,1):ib(im,1))
!~         print*, '333',beta * v_x(ia(im,1):ib(im,1))
!~         pause
!~     endif

    !$omp end parallel

endsubroutine

subroutine a_heating()
implicit none
!~ double precision :: v_recul(1) = (/8.0e-5, 9.0e-5/)
double precision :: v_recul(1) = 4.0e-2
!dir$ if ( HCI .eq. 0 )
double precision, dimension(n_ions(1)) :: ran_aux1,ran_aux2, theta_em,phi_em, sinphi_em
!dir$ else
double precision, dimension(n_cut(3)) :: ran_aux1,ran_aux2, theta_em,phi_em, sinphi_em
!dir$ endif
integer :: im
!  use vsRngUniform() for simple precision
!  use vdRngUniform() for double precision

    errcode = vdRngUniform(method, stream, n_ions(1), theta_em, 0.0d0, pi2 )
    errcode = vdRngUniform(method, stream, n_ions(1), ran_aux2,-1.0d0, 1.0d0)

    !$omp parallel default(none) private(im)&
    !$omp shared (ia,ib, v_x,v_y,v_z,theta_em,phi_em, sinphi_em, ran_aux2, v_recul)
    im = omp_get_thread_num()+1

!~     theta_em (ia(im, 1):ib(im, 1)) = ran_aux1(ia(im, 1):ib(im, 1))
    phi_em   (ia(im, 1):ib(im, 1)) = dacos(ran_aux2(ia(im, 1):ib(im, 1)))
    sinphi_em(ia(im, 1):ib(im, 1)) = dsin(phi_em(ia(im, 1):ib(im, 1)))
!~     v_x(ia(im, 1):ib(im, 1)) = v_x(ia(im, 1):ib(im, 1)) + sign(1.0d0,v_x(ia(im, 1):ib(im, 1)))*v_recul(1)*abs(sinphi_em(ia(im, 1):ib(im, 1))*dcos(theta_em(ia(im, 1):ib(im, 1))) )
!~     v_y(ia(im, 1):ib(im, 1)) = v_y(ia(im, 1):ib(im, 1)) + sign(1.0d0,v_y(ia(im, 1):ib(im, 1)))*v_recul(1)*abs(sinphi_em(ia(im, 1):ib(im, 1))*dsin(theta_em(ia(im, 1):ib(im, 1))))
!~     v_z(ia(im, 1):ib(im, 1)) = v_z(ia(im, 1):ib(im, 1)) + sign(1.0d0,v_z(ia(im, 1):ib(im, 1)))*v_recul(1)*abs(dcos(phi_em(ia(im, 1):ib(im, 1))))

    v_x(ia(im, 1):ib(im, 1)) = v_x(ia(im, 1):ib(im, 1)) + v_recul(1)*sinphi_em(ia(im, 1):ib(im, 1))*dcos(theta_em(ia(im, 1):ib(im, 1)))
    v_y(ia(im, 1):ib(im, 1)) = v_y(ia(im, 1):ib(im, 1)) + v_recul(1)*sinphi_em(ia(im, 1):ib(im, 1))*dsin(theta_em(ia(im, 1):ib(im, 1)))
    v_z(ia(im, 1):ib(im, 1)) = v_z(ia(im, 1):ib(im, 1)) + v_recul(1)*dcos(phi_em(ia(im, 1):ib(im, 1)))
    !$omp end parallel
endsubroutine

subroutine save_trj_to_file()
implicit none
integer i
    print*, 'j_save_temp Trj',j_save_temp, jj
    open(unit = 11, status='replace',file=trim(adjustl(str_file_trj))//'.bin',form='unformatted')  ! create a new file, or overwrite an existing one
        write(11) j_save_temp
        write(11) n_ions(1)
        write(11) n_ions(2)
        do i=1,10
            write(11) save_trj(1:j_save_temp,i)
        enddo
    close(11)

    deallocate(save_trj)
end subroutine

subroutine save_data()
implicit none
double precision, allocatable, dimension(:,:)   ::  save_final_position
integer :: i

    print*, trim(adjustl(str_file_xva))

    allocate(save_final_position(12,n_cut(3)))

    save_final_position(1,:) = r_x ; save_final_position(2,:) = r_y ; save_final_position(3,:) = r_z
    save_final_position(4,:) = v_x ; save_final_position(5,:) = v_y ; save_final_position(6,:) = v_z
    save_final_position(7,:) = a1_x; save_final_position(8,:) = a1_y; save_final_position(9,:) = a1_z

    save_final_position(10,:)  = v_rf_avg(:,1);
    save_final_position(11,:)  = v_rf_avg(:,2);
    save_final_position(12,:)  = v_rf_avg(:,3);

    open(unit = 10, status='replace',file=trim(adjustl(str_file_xva))//'.bin',form='unformatted')  ! create a new file, or overwrite an existing one
        write(10) iRF
        write(10) save_final_position
        write(10) state
        write(10) total_ion_lost0
    close(10)
    deallocate(save_final_position)

    open(10, file=trim(adjustl(str_file_xva))//'.info', status='replace', access='sequential', action='write')
        write(10,'(i16  , 3x, "%Number of LC ions")')                   n_ions(1)
        write(10,'(i16  , 3x, "%Number of HC ions")')                   n_ions(2)
        write(10,'(i16  , 3x, "%Last index in the integration loop")')  j_start
        write(10,'(e16.9, 3x, "%Last time")')                           t_act
        write(10,'(i16  , 3x, "%Mass   of the LC ions")')               int(mass(1)/amu)
        write(10,'(i16  , 3x, "%Charge of the LC ions")')               int(charge(1)/qe)
        write(10,'(i16  , 3x, "%Mass   of the HC ions")')               int(mass(2)/amu)
        write(10,'(i16  , 3x, "%Charge of the HC ions")')               int(charge(2)/qe)

        write(10,'(e16.9  , 3x, "%V_st[V]")')                             V_st
        write(10,'(e16.9  , 3x, "%V_rf[V]")')                             V_rf
        write(10,'(e16.9  , 3x, "%Omega/2pi[Hz]")')                       Omega/pi2
        write(10,'(e16.9  , 3x, "%r_0[m]")')                              r_0
        write(10,'(e16.9  , 3x, "%L_0[m]")')                              L
        write(10,'(e16.9  , 3x, "%wz_LC/2 pi[Hz]")')                      wz_LC / pi2

        write(10,'(i16  , 3x, "%i_free__fly")')                      i_free__fly
        write(10,'(i16  , 3x, "%i_therm_fly")')                      i_therm_fly
        write(10,'(i16  , 3x, "%i_cool__fly")')                      i_cool__fly
        write(10,'(i16  , 3x, "%i_cool_heat")')                      i_cool_heat
        write(10,'(i16  , 3x, "%i_laser_fly")')                      i_laser_fly
!~         write(10,'(i16  , 3x, "%i_laser_fly")')                      '##v##ADRIEN##v##'
!~         write(10,'(e16.9  , 3x, "%detuning")')                       detuning
!~         write(10,'(e16.9  , 3x, "%saturation")')                     saturation

    close(10)


!--------- Saving the Temperature Evolution -----------------------
!~     print*, 'Saving the Temperature Evolution in:',str_file_Temp_init
    open(10, file = trim(adjustl(str_file_Temp))//'.dat', status='replace', access='sequential', action='write')
    print*, 'j_save_temp Temp',j_save_temp
        do i = 1, j_save_temp
            !dir$ if ( HCI .eq. 0 )
            write(10,221) save_temperature(i,0),&
                          m_kb_x_inv_n_ions2(1)*save_temperature(i,1:3), & ! T_CM
                          m_kb_x_inv_n_ions(1) *save_temperature(i,4:6), & ! T_aux
                          save_temperature(i,7)/real(n_dt)
            !dir$ else
            write(10,221) save_temperature(i,0),&
                          m_kb_x_inv_n_ions2(1)*save_temperature(i,1:3), &
                          m_kb_x_inv_n_ions(1) *save_temperature(i,4:6), &
                          save_temperature(i,7)/real(n_dt), &
                          m_kb_x_inv_n_ions2(2)*save_temperature(i, 8:10), &
                          m_kb_x_inv_n_ions(2) *save_temperature(i,11:13)
            !dir$ endif
        enddo

        !dir$ if ( HCI .eq. 0 )
221     format( 8(1X,e27.19e3))
        !dir$ else
221     format(14(1X,e27.19e3))
        !dir$ endif
    close(10)

    deallocate(save_temperature)

endsubroutine

subroutine generate_gauss(N,rr)
implicit none
integer         , intent(in) :: N
double precision, intent(out):: rr(N)
double precision             :: x1,x2,s
integer                      :: i
    do i = 1, N, 2
        s = 10
        do while (s.ge.1)
            call random_number(x1)
            call random_number(x2)
!~             x1 = 2.0*rand() - 1.0
!~             x2 = 2.0*rand() - 1.0
            x1 = 2.0*x1 - 1.0
            x2 = 2.0*x2 - 1.0
            s  = x1*x1 + x2*x2
        enddo
        s  = dsqrt(-2.0d0*dlog(s)/s)
        rr(i  ) = x1*s
        if ((i+1).le.N) then
            rr(i+1) = x2*s
        endif
    enddo
endsubroutine

subroutine initialize_ions_LC() !Done
implicit none
double precision, parameter :: l0_w(3) = (/0.2*r_0, 0.2*r_0, 0.05*L/) ! Initial gaussian distribution
!~ double precision, parameter :: temperature_init = 1.0d-3
!~ double precision, parameter :: sigma = dsqrt(kb*temperature_init / mass(1))

    call generate_gauss(n_cut(1), r_x); r_x = l0_w(1)*r_x
    call generate_gauss(n_cut(1), r_y); r_y = l0_w(2)*r_y
    call generate_gauss(n_cut(1), r_z); r_z = l0_w(3)*r_z
    
!~     do i=1,n_cut(1)
!~         print*, r_x(i),r_z(i),r_y(i)
!~     enddo
!~     stop
    
    v_x = 0.0d0;
    v_y = 0.0d0;
    v_z = 0.0d0;

    a1_x = 0.0d0;
    a1_y = 0.0d0;
    a1_z = 0.0d0;

    a2_x = 0.0d0;
    a2_y = 0.0d0;
    a2_z = 0.0d0;
end subroutine

subroutine inject_ions_HC() !Done
implicit none
double precision :: l0(3)!  = (/0.0d0, 0.0d0, -z0_HC /) ! Pos. of the inital distribution
double precision, parameter :: eta = 1. / (4.*(dsqrt(10.D0) + 5.))
double precision :: Ep0, Ek_init
!~ double precision :: d_min
    l0            = (/0.0d0, 0.0d0, -z0_HC/)
    r_x(n_cut(3)) = l0(1)
    r_y(n_cut(3)) = l0(2)
!~     r_z(n_cut(3)) = l0(3)*(1 + 0.1*(2*rand() - 1))
    r_z(n_cut(3)) = l0(3)
    rz_0  = -r_z(n_cut(3))

    v_x(n_cut(3)) = 0.0e0;
    v_y(n_cut(3)) = 0.0e0;
    Ep0   = wz2(2)*eta*L**2*mass(2)*(1- exp(-r_z(n_cut(3))**2 / (2*eta*L**2)) ) ! [J]
    Ek_init = E_critic*qe * frac_E0_Ecr - Ep0
    if (Ek_init > 0) then
        v_z(n_cut(3)) = dsqrt(Ek_init * 2 / mass(2)) ;
    else
        print*, 'Initial velocity negative, stopping code'
        stop
    endif

!~     print*, 'E_p0', Ep0
!~     print*, E_init
!~     print*, v_z(n_cut(3))
!~     pause

    a1_x(n_cut(3)) = 0.0d0;
    a1_y(n_cut(3)) = 0.0d0;
    a1_z(n_cut(3)) = 0.0d0;

    a2_x(n_cut(3)) = 0.0d0;
    a2_y(n_cut(3)) = 0.0d0;
    a2_z(n_cut(3)) = 0.0d0;
end subroutine

subroutine distribute_ions( )
implicit none
integer, dimension(n_sp+1)     :: i
integer                        :: nii, im, j
integer            :: n_by_core(3)

    if (ni.ne.omp_get_max_threads ( )) then
        print*, 'ni (number of threads) is not the same in the makefile and in the code!'
        print*, ni, omp_get_max_threads ( )
        print*, 'Please, correct it, recompile and re-run'
        print*, 'Stoping code'
        stop
    endif

    i = modulo(n_cut, ni)
    n_by_core(1) = floor(n_cut(1) / real(ni))
    n_by_core(2) = floor(n_cut(2) / real(ni))
    n_by_core(3) = floor(n_cut(3) / real(ni))
    ! For some strange reason, if I use vectorial notation, 5.00 becomes 4, if I do element by element 5.00 becomes 5. No idea why!

    do j = 1, 3
        if (n_cut(j) < ni) then
            nii = n_cut(j)
            ia(:,j) = 1; ib(:,j) = 0; ! necessary for the ni > n_cut
            do im = 1, nii
                ia(im,j) = im
                ib(im,j) = im
            enddo
        else
            do im = 1, ni
                if (i(j) == 0) then
                    ia(im,j) = (im-1)*n_by_core(j)+1
                    ib(im,j) = ia(im,j) + n_by_core(j)-1
                else
                    if (im <= i(j)) then
                        ia(im,j) = (im-1)*(n_by_core(j)+1)+1
                        ib(im,j) = ia(im,j) + n_by_core(j)
                    else
                        ia(im,j) = n_by_core(j)*(im - 1) + i(j) + 1
                        ib(im,j) = ia(im,j) + n_by_core(j) - 1
                    endif
                endif
            enddo
        endif
    enddo

    ia(:,2)  = ia(:,2)  + n_cut(1)
    ib(:,2)  = ib(:,2)  + n_cut(1)

!~     print*, 'LC '; print*, ia(:,1); print*, ib(:,1)
!~     print*, 'HC '; print*, ia(:,2); print*, ib(:,2)
!~     print*, 'ALL'; print*, ia(:,3); print*, ib(:,3)
endsubroutine

subroutine init_from_file()
implicit none
integer :: jend, n_aux, n_ions_LC, n_ions_HC
double precision, allocatable, dimension(:,:) :: xyz_temp

    print*, 'File to load:',trim(adjustl(str_file_to_load))

    open(10, file=trim(adjustl(str_file_to_load))//'.info', status='old', access='sequential', action='read')
        read(10,'(i16)')    n_ions_LC
        read(10,'(i16)')    n_ions_HC
        read(10,'(i16)')    jend
        read(10,'(e16.9)')  t_act
    close(10)

    t_next_periode = t_act + dt0

    n_aux = n_ions_LC + n_ions_HC
    allocate(xyz_temp(12,n_aux))
    open(unit = 10, status='old', file=trim(adjustl(str_file_to_load))//'.bin', form='unformatted')  ! open an existing file
        read(10) iRF
        read(10) xyz_temp
        read(10) state
        read(10) total_ion_lost0
    close(10)

    r_x (1:n_aux) = xyz_temp(1,:)
    r_y (1:n_aux) = xyz_temp(2,:)
    r_z (1:n_aux) = xyz_temp(3,:)

    v_x (1:n_aux) = xyz_temp(4,:)
    v_y (1:n_aux) = xyz_temp(5,:)
    v_z (1:n_aux) = xyz_temp(6,:)

    a1_x (1:n_aux) = xyz_temp(7,:)
    a1_y (1:n_aux) = xyz_temp(8,:)
    a1_z (1:n_aux) = xyz_temp(9,:)

    v_rf_avg(:n_ions_LC,1) = xyz_temp(10,:n_ions_LC)
    v_rf_avg(:n_ions_LC,2) = xyz_temp(11,:n_ions_LC)
    v_rf_avg(:n_ions_LC,3) = xyz_temp(12,:n_ions_LC)

    if (n_ions_HC.eq.0) then
        v_rf_avg(n_ions_LC+1:,1) = 0.0d0
        v_rf_avg(n_ions_LC+1:,2) = 0.0d0
        v_rf_avg(n_ions_LC+1:,3) = 0.0d0
    else
        v_rf_avg(n_ions_LC+1:,1) = xyz_temp(10,n_ions_LC+1:)
        v_rf_avg(n_ions_LC+1:,2) = xyz_temp(11,n_ions_LC+1:)
        v_rf_avg(n_ions_LC+1:,3) = xyz_temp(12,n_ions_LC+1:)
    endif

    deallocate(xyz_temp)

end subroutine

subroutine init_from_file_post_injec()
implicit none
integer :: jend, n_aux, n_ions_LC, n_ions_HC
double precision, allocatable, dimension(:,:) :: xyz_temp

    print*, 'File to load:',trim(adjustl(str_file_to_load))
    open(10, file=trim(adjustl(str_file_to_load))//'.info', status='old', access='sequential', action='read')
        read(10,'(i16)')    n_ions_LC
        read(10,'(i16)')    n_ions_HC
        read(10,'(i16)')    jend
        read(10,'(e16.9)')  t_act
    close(10)

    t_next_periode = t_act + dt0
    flag_dt0    = .false.

    n_aux = n_ions_LC + n_ions_HC
    allocate(xyz_temp(12,n_aux))
    open(unit = 10, status='old', file=trim(adjustl(str_file_to_load))//'.bin', form='unformatted')  ! open an existing file
        read(10) iRF
        read(10) xyz_temp
        read(10) state
!~         read(10) total_ion_lost0
    close(10)

    r_x (1:n_ions_LC) = xyz_temp(1,:)
    r_y (1:n_ions_LC) = xyz_temp(2,:)
    r_z (1:n_ions_LC) = xyz_temp(3,:)

    v_x (1:n_ions_LC) = xyz_temp(4,:)
    v_y (1:n_ions_LC) = xyz_temp(5,:)
    v_z (1:n_ions_LC) = xyz_temp(6,:)

    a1_x (1:n_ions_LC) = xyz_temp(7,:)
    a1_y (1:n_ions_LC) = xyz_temp(8,:)
    a1_z (1:n_ions_LC) = xyz_temp(9,:)

    v_rf_avg(:n_ions_LC,1) = xyz_temp(10,:n_ions_LC)
    v_rf_avg(:n_ions_LC,2) = xyz_temp(11,:n_ions_LC)
    v_rf_avg(:n_ions_LC,3) = xyz_temp(12,:n_ions_LC)

    if (n_ions_HC.eq.0) then
        print*, 'This is not a post injection file, as N_HC = 0'
        print*, 'Stopping code execution'
        stop
    endif

    deallocate(xyz_temp)

end subroutine

subroutine init_from_file_post_cooling()
implicit none
integer :: jend, n_aux, n_ions_LC, n_ions_HC
double precision, allocatable, dimension(:,:) :: xyz_temp

    print*, 'File to load:',trim(adjustl(str_file_to_load))
    open(10, file=trim(adjustl(str_file_to_load))//'.info', status='old', access='sequential', action='read')
        read(10,'(i16)')    n_ions_LC
        read(10,'(i16)')    n_ions_HC
        read(10,'(i16)')    jend
        read(10,'(e16.9)')  t_act
    close(10)

    t_next_periode = t_act + dt0
    flag_dt0    = .false.

    n_aux = n_ions_LC + n_ions_HC
    allocate(xyz_temp(12,n_aux))
    open(unit = 10, status='old', file=trim(adjustl(str_file_to_load))//'.bin', form='unformatted')  ! open an existing file
        read(10) iRF
        read(10) xyz_temp
        read(10) state
!~         read(10) total_ion_lost0
    close(10)

    r_x (1:n_ions_LC) = xyz_temp(1,:)
    r_y (1:n_ions_LC) = xyz_temp(2,:)
    r_z (1:n_ions_LC) = xyz_temp(3,:)

    v_x (1:n_ions_LC) = xyz_temp(4,:)
    v_y (1:n_ions_LC) = xyz_temp(5,:)
    v_z (1:n_ions_LC) = xyz_temp(6,:)

    a1_x (1:n_ions_LC) = xyz_temp(7,:)
    a1_y (1:n_ions_LC) = xyz_temp(8,:)
    a1_z (1:n_ions_LC) = xyz_temp(9,:)

    v_rf_avg(:n_ions_LC,1) = xyz_temp(10,:n_ions_LC)
    v_rf_avg(:n_ions_LC,2) = xyz_temp(11,:n_ions_LC)
    v_rf_avg(:n_ions_LC,3) = xyz_temp(12,:n_ions_LC)

    if (n_ions_HC.eq.1) then
        print*, 'This is a post injection file, as N_HC = 1'
        print*, 'Stopping code execution'
        stop
    endif

    deallocate(xyz_temp)

end subroutine

subroutine create_files_names( )
implicit none
character(len=10 )  :: str_N,   str_N_aux
character(len=10 )  :: str_D,   str_D_aux   !Detuning
character(len=10 )  :: str_S,   str_S_aux   !Detuning
character(len=10 )  :: str_Vrf, str_Vrf_aux
character(len=19 )  :: str_Udc, str_Udc_aux
character(len=130)  :: str_file_aux
character(len=20 )   :: str_trap_aux,str_stat


! N_ions
    if (n_cut(1)<10) then
        write(str_N,"(I5.5)") n_cut(1)
    elseif(n_cut(1) < 100) then
        write(str_N,"(I5.5)") n_cut(1)
    elseif(n_cut(1) < 1000) then
        write(str_N,"(I5.5)") n_cut(1)
    elseif(n_cut(1) < 10000) then
        write(str_N,"(I5.5)") n_cut(1)
    elseif(n_cut(1) < 100000) then
        write(str_N,"(I5.5)") n_cut(1)
    endif
    str_N_aux = '_N'//trim(adjustl(str_N))

! Vrf
    write(str_Vrf,"(I4.4)") int(V_rf)
    str_Vrf_aux = '_Vrf'//trim(adjustl(str_Vrf))

! Udc
    write(str_Udc,"(d13.4)") Udc
    str_Udc_aux = '_Udc'//trim(adjustl(str_Udc))//'V'

!! Udc
!    write(str_Udc,"(f6.0)") 1.d-3*wz_LC/pi2      sqrt(wz2(1))
!    str_Udc_aux = '_wz'//trim(adjustl(str_Udc))//'2pikHz'

! Detuning
    write(str_D,"(f4.1)") abs(detuning)
    str_D_aux = '_D'//trim(adjustl(str_D))

! Saturation:
    write(str_S,"(f4.1)") abs(saturation)
    str_S_aux = '_S'//trim(adjustl(str_S))

!dir$     if(Simu_type .eq. 0)
    str_trap_aux = 'SimuType0'
!dir$ elseif(Simu_type .eq. 1)
    str_trap_aux = 'SimuType1'
!dir$ elseif(Simu_type .eq. 2)
    write(str_stat,"(I2.2)") i
    str_trap_aux = 'SimuType4_'//trim(adjustl(str_stat))
    ! Names used by the "from scrath" subroutine
    str_file_aux = trim(adjustl(str_trap_aux)) &
                // trim(adjustl(str_N_aux))    &
                // trim(adjustl(str_Vrf_aux))  &
                // trim(adjustl(str_Udc_aux))  &
                // trim(adjustl(str_D_aux))    &
                // trim(adjustl(str_S_aux))

    str_file_to_load = 'xva_'//trim(adjustl(str_file_aux)) &
                             //trim(adjustl(str_extra))    &
                             //trim(adjustl(str_extra2))
    write(str_stat,"(I2.2)") i
    str_trap_aux = 'SimuType2_'//trim(adjustl(str_stat))
!dir$ elseif(Simu_type .eq. 3)
    str_trap_aux = 'SimuType3'
!dir$ elseif(Simu_type .eq. 4)
    str_file_aux = 'SimuType0' &
                // trim(adjustl(str_N_aux))    &
                // trim(adjustl(str_Vrf_aux))  &
                // trim(adjustl(str_Udc_aux))  &
                // trim(adjustl(str_D_aux))    &
                // trim(adjustl(str_S_aux))

    str_file_to_load = 'xva_'//trim(adjustl(str_file_aux)) &
                                //trim(adjustl(str_extra))
    write(str_stat,"(I2.2)") i
    str_trap_aux = 'SimuType4_'//trim(adjustl(str_stat))
!dir$ elseif(Simu_type .eq. 6)
    write(str_stat,"(I2.2)") i
    str_trap_aux = 'SimuType0'
    ! Names used by the "from scrath" subroutine
    str_file_aux = trim(adjustl(str_trap_aux)) &
                // trim(adjustl(str_N_aux))    &
                // trim(adjustl(str_Vrf_aux))  &
                // trim(adjustl(str_Udc_aux))  &
                // trim(adjustl(str_D_aux))    &
                // trim(adjustl(str_S_aux))

    str_file_to_load = 'xva_'//trim(adjustl(str_file_aux)) &
                             //trim(adjustl(str_extra))    &
                             //trim(adjustl(str_extra2))
    write(str_stat,"(I2.2)") i
    str_trap_aux = 'SimuType6_'//trim(adjustl(str_stat))
!dir$ endif

    ! Names used by the "from scrath" subroutine
    str_file_aux = trim(adjustl(str_trap_aux)) &
                // trim(adjustl(str_N_aux))    &
                // trim(adjustl(str_Vrf_aux))  &
                // trim(adjustl(str_Udc_aux))  &
                // trim(adjustl(str_D_aux))    &
                // trim(adjustl(str_S_aux))    &
                //trim(adjustl(str_extra))     &
                //trim(adjustl(str_extra2))

    str_file_Temp = 'Temp_'//trim(adjustl(str_file_aux))
    str_file_xva  =  'xva_'//trim(adjustl(str_file_aux))
    str_file_trj  =  'trj_'//trim(adjustl(str_file_aux))
    str_file_PM   =  'PM__'//trim(adjustl(str_file_aux))
                                
     

    print*, str_file_Temp
end subroutine



end program main
