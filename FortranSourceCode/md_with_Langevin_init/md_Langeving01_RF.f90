!~ Starting from md_v006_2D
!~ from md_Langevin_v73_1.f90

! I am using the vGB82 implemetation of the Langevin integrator
! ref: Skeel and Iwaguirrem Molecular Physics, v100, n24, p3885-3891  (2002)


!~ fx 360618.0563259394 Hz
!~ fy 360618.0563259394 Hz
!~ fz 240252.7346216562 Hz


include 'mkl_vsl.f90'
program main
!~ use ifport, only : random, seed
USE MKL_VSL_TYPE
USE MKL_VSL
use omp_lib
implicit none


!dir$ define potential_type = 2  ! 0 = Harmonique ; 1 = Polyfit; 2=Trigamma
!dir$ define simu_type = 3
!dir$ define integration_methode = 1  ! 0=Velocyt-Verlet; 1=Langeving_methode[vGB82]; 2=Langeving_methode[Impulse]

!! The Velocity-Verlet is not working well in the sense that I have only implemented viscosity, so there
!! is no source of heating....

!***********************************************************************
!               Math and Phys cte [SI units]:
!***********************************************************************
double precision , parameter :: pi     = 3.141592653589793d0
double precision , parameter :: pi2    = 2.d0*pi
double precision , parameter :: inv_2pi= 1.0d0/pi2
double precision , parameter :: amu    = 1.66053886d-27   ![Kg]
double precision , parameter :: c      = 2.99792458d8     ! vitesse de la lumière
double precision , parameter :: hbar   = 1.055d-34
double precision , parameter :: hplanck= 6.62d-34
double precision , parameter :: qe     = 1.60217646d-19   ![C = sA] Elementary Charge
double precision , parameter :: kb     = 1.3806503d-23    ![m2 kg s-2] Boltzman cte
double precision , parameter :: ke     = 8.987551787d9    ![N m2 C-2] Coulomb cte =  1 / (4*pi*eps0)
!***********************************************************************
integer          , parameter :: ni   = 8 !number of threads

!***********************************************************************
!                   Particles parameters:
!***********************************************************************
double precision , parameter :: charge = qe*1
double precision , parameter :: mass   = amu*40
!~ double precision , parameter :: mass   = amu*172 !Yterbium
double precision , parameter :: q_m   = charge / mass
! Coulomb interaction cte (Not- Normalised):
double precision , parameter :: alpha = ke*charge*charge / mass

!***********************************************************************
!                   Polynomial fit parameters:
!***********************************************************************
!~ integer          , parameter :: n_ions    = 16
!~ integer          , parameter :: ia(ni)=[1,5,9,13];
!~ integer          , parameter :: ib(ni)=[4,8,12,16];

!~ double precision , parameter :: wr2  = (2*pi*800e3)**2
!~ double precision , parameter :: wy2  = (2*pi*800e3)**2
!~ double precision , parameter :: wz2  = (2*pi*400e3)**2

integer          , parameter :: n_ions    = 1024
! Automatiser les ia et ib
! omp parallel quelle différence variables private, shared, firstprivate ?
integer          , parameter :: ia(ni)=[1,129,257,385,513,641,769,897];
integer          , parameter :: ib(ni)=[128,256,384,512,640,768,896,1024];


!***********************************************************************
!                   Simulations parameters:
!***********************************************************************
integer           , parameter :: i_free__fly = 5*10**5

double precision  , parameter :: dt   = 2d-9
double precision  , parameter :: dt1  = 0.5*dt
double precision  , parameter :: dt2  = 0.5*dt*dt

double precision  , parameter :: Temperature = 0.5d-3

!~ !dec$     if defined(E1)
double precision  , parameter :: eta = 1.5d-21  ![kg/s] friction  coefficient
!~ !dec$ endif

double precision  , parameter :: eta_mass = eta/mass

!***********************************************************************
!                   Trapping parameters (Adrien 2020 07 06) :
!***********************************************************************

! Trap dimensions
double precision   , parameter :: r_0   = 2.865d-3/1.14511d0
double precision   , parameter :: L     = 0.0140608827209d0
double precision   , parameter :: d_0   = 4d-3 ! longueur du piege
! Trapping voltages
double precision  , parameter :: V_st  = 0.00d0
double precision  , parameter :: Omega = pi2 * 2.000000d006

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

double precision   , parameter :: wz_LC = pi2*90806.9982303 !when 1V applied
double precision   , parameter :: wz2   = Udc * wz_LC**2

double precision   , parameter :: a_cte_Quad_RF_LC(4) = (/ &
    -2*charge*V_st / (mass*r_0**2 ) + 0.5*wz2,    &
    +2*charge*V_st / (mass*r_0**2 ) + 0.5*wz2,    &
     2*charge*V_rf / (mass*r_0**2 )             ,    &
     -wz2                                               /)

!***********************************************************************
!                   Trap parameters:
!***********************************************************************
!~ double precision , parameter :: wz2  = (2*pi*25389.931964385545)**2 ! N=30, d0 = 10um
!~ double precision , parameter :: wz2  = (2*pi*25000.0)**2 ! N=30, d0 = 10um

!***********************************************************************
!                   Integration constants:
!***********************************************************************
!                   Integration constants Langevin:
!***********************************************************************
!~ double precision             :: eta = 1.5d-21  ![kg/s] friction  coefficient


!~ double precision , parameter :: eta = 1.5d-21  ![kg/s] friction  coefficient
double precision  :: alphaL,betaL,time00, time_a, t_act
!***********************************************************************
!                       ! Global variables:
!***********************************************************************
double precision   , dimension(n_ions)   ::  r_z, r_r, r_y,   &
                                             v_z, v_r, v_y,   &
                                             a1_z, a1_r,a1_y, &
                                             a2_z, a2_r,a2_y, &
                                             Z0_r, Z1_r, &
                                             Z0_z, Z1_z, &
                                             Z0_y, Z1_y

!~ integer, dimension(ni)   :: ia , ib
character(len=150)       :: str_file_aux
!***********************************************************************
!~ double precision , parameter :: rSum_Adb = 188.57d-06 !(for N = 30ions) Inhomo

integer          , parameter :: j_save_next  = 100

character(len=100),parameter :: str_extra   = '_4'

double precision  , parameter ::    w_p   = (dexp(-eta_mass*dt)+eta_mass*dt-1.0d0) /(eta_mass*dt*(1.0d0-dexp(-eta_mass*dt)))
double precision  , parameter ::    w_m   = 1-w_p
double precision  , parameter ::    cte1 = (kb*Temperature/mass) * (2*w_p*w_p * eta_mass * dt + w_p - w_m)
double precision  , parameter ::    cte2 = (kb*Temperature/mass) * (2*w_p*w_m * eta_mass * dt + w_m - w_p) ! b
double precision  , parameter ::    cte3 = (kb*Temperature/mass) * (2*w_m*w_m * eta_mass * dt + w_p - w_m) ! c
!dir$ if (integration_methode.eq.1) ! Langeving
double precision  , parameter ::    cte4 = w_p - 0.5d0 ! S vGB82
double precision  , parameter ::    dt_S  = dt*cte4
double precision  , parameter ::    dt_S2 = dt+dt_S
!dir$ elseif (integration_methode.eq.2) ! Langeving
double precision  , parameter ::    cte4 = 0           ! Impulse method
double precision  , parameter ::    dt_S = 0
double precision  , parameter ::    dt_S2 = dt
!dir$ endif
double precision  , parameter ::    e_eta_mass2 = dexp(-eta_mass*dt*0.5d0)
double precision  , parameter ::    e_eta_mass3 = (1.0d00-dexp(-eta_mass*dt))/(eta_mass*e_eta_mass2)
double precision  , parameter ::    cte13   = cte1 + cte3 !ac
!***********************************************************************

integer :: j_start, jj, j_end, j_save, j_save_aux

double precision, allocatable, dimension(:,:) ::  save_Temp
!~ double precision, allocatable, dimension(:,:) ::  save_trj
real :: rand_num(1)
integer :: rand_seed

!*******************************************************************************
! For: vdRngGaussian
!*******************************************************************************
TYPE (VSL_STREAM_STATE) :: stream
integer                      :: seed
integer brng,method
integer(kind=4) errcode

    if (cte1<0) then
         print*, 'Cte1 < 0', cte1
            print*, 'Stopping code'
            stop
    endif

brng=VSL_BRNG_MCG31
method=VSL_RNG_METHOD_GAUSSIAN_ICDF
!~ seed=777
!     ***** Initialize *****
!~ 	CALL RANDOM_SEED(GET=seed) ! Processor initializes the seed randomly from the date and time
!~       errcode=vslnewstream( stream, brng,  seed )
!~       
!~     call random_seed(seed)

    CALL RANDOM_SEED()
    call random_number(rand_num)
    rand_seed = int(rand_num(1)*1e4)
    errcode=vslnewstream( stream, brng, rand_seed )

!~     stop
    time00 = omp_get_wtime()
    call distribute_ions()

    print*, '***********************************************************'
    print*, '        From Zero and Save to file'
    print*, '***********************************************************'
    jj = int(i_free__fly/float(j_save_next))
    allocate(save_Temp(3,jj))

    call initialize_ions()
    
    call update_accelerations()
    a1_r = a2_r; a2_r = 0.0d0;
    a1_y = a2_y; a2_y = 0.0d0;
    a1_z = a2_z; a2_z = 0.0d0;

    j_start = 1; j_end = i_free__fly;
    call free_fly()
    call create_files_names()
    call save_xva_to_file()
    call save_Temp_to_file()

    time_a = omp_get_wtime()
    print*, 'Total time', time_a - time00
    print*, 'Avg Time per step', (time_a - time00)/i_free__fly
contains

subroutine free_fly()
implicit none
    call zero_time_step()
    do jj = j_start , j_start + j_end - 2
        call update_positions()
        call update_accelerations()
        call update_velocities()
        if (j_save_aux.eq.j_save_next) then
            j_save_aux = 1
            j_save = j_save + 1
            save_Temp(1,j_save) = sum(v_z*v_z)
            save_Temp(2,j_save) = sum(v_r*v_r)
            save_Temp(3,j_save) = sum(v_y*v_y)
        else
            j_save_aux = j_save_aux + 1
        endif
    enddo
    call final_time_step()
    j_start = jj;
endsubroutine

subroutine zero_time_step()
implicit none

    alphaL  = dsqrt(cte1)

    errcode=vdrnggaussian( method, stream, n_ions,Z0_r, 0.0d0, 1.0d0)
    errcode=vdrnggaussian( method, stream, n_ions,Z0_z, 0.0d0, 1.0d0)
    errcode=vdrnggaussian( method, stream, n_ions,Z0_y, 0.0d0, 1.0d0)

    v_z  = e_eta_mass2 * (v_z + dt*w_p*a2_z + alphaL*Z0_z) ! v_{1/2}
    v_r  = e_eta_mass2 * (v_r + dt*w_p*a2_r + alphaL*Z0_r) ! v_{1/2}
    v_y  = e_eta_mass2 * (v_r + dt*w_p*a2_r + alphaL*Z0_y) ! v_{1/2}
end subroutine

subroutine update_positions()
implicit none
integer :: im
!-----------------------------------------------------------------------
	t_act = t_act + dt
!   One time step
    !$omp parallel default(none) private(im) shared (r_z,v_z,r_r,v_r,r_y,v_y) firstprivate(ia,ib)
    im = omp_get_thread_num()+1
    ! This falls into the "fused multiply add" category: U = U + B*X
    ! It could be speed up!
    r_r(ia(im): ib(im)) = r_r(ia(im): ib(im)) + v_r(ia(im): ib(im))*e_eta_mass3;
    r_z(ia(im): ib(im)) = r_z(ia(im): ib(im)) + v_z(ia(im): ib(im))*e_eta_mass3;
    r_y(ia(im): ib(im)) = r_y(ia(im): ib(im)) + v_y(ia(im): ib(im))*e_eta_mass3;
    !$omp end parallel
end subroutine

subroutine update_accelerations()
implicit none
!---------------  Locals  ----------------------------------------------
integer :: i, im,j
double precision :: rji(3)
double precision :: r2inv
double precision, parameter :: softening = 1.00D-20
! --- Added from RF_relax by Adrien 2020 07 06 ---
double precision, parameter :: eta_gauss = 1.0d0 / (4.0d0*(dsqrt(10.0d0) + 5.0d0 ))
double precision, parameter :: cte_gauss = -1/(2*eta_gauss*L**2)
double precision    :: cte_aux(2)
double precision    :: cos_jj
    cos_jj = dcos(t_act*Omega)
    cte_aux(1) = a_cte_Quad_RF_LC(1) + a_cte_Quad_RF_LC(3)*cos_jj
    cte_aux(2) = a_cte_Quad_RF_LC(2) - a_cte_Quad_RF_LC(3)*cos_jj
! -------------------------------------------------
    !$omp parallel default(none) &
    !$omp private(im, i,j,rji,r2inv,cte_aux) &
    !$omp firstprivate(ia, ib, r_z, r_r,r_y) &
    !$omp shared(a2_z,a2_r,a2_y)
    im = omp_get_thread_num() + 1
! RF contribution to acceleration
    a2_r(ia(im):ib(im)) = a2_r(ia(im):ib(im)) + cte_aux(1) * r_r(ia(im):ib(im))
    a2_y(ia(im):ib(im)) = a2_y(ia(im):ib(im)) + cte_aux(2) * r_y(ia(im):ib(im))
    a2_z(ia(im):ib(im)) = a2_z(ia(im):ib(im)) -     wz2 * r_z(ia(im):ib(im))*exp(cte_gauss * r_z(ia(im):ib(im))*r_z(ia(im):ib(im)) )
! Pseudo-pot contribution to acceleration
	!    a2_r(ia(im): ib(im)) = - wr2*r_r(ia(im): ib(im))
	!    a2_y(ia(im): ib(im)) = - wy2*r_y(ia(im): ib(im))
	!    a2_z(ia(im): ib(im)) = - wz2*r_z(ia(im): ib(im))
! Coulombian repulsion see a_Coulomb_1sp_2sp()
    do i = ia(im), ib(im)
        do j = 1, n_ions
            rji(1)  = r_z(j) - r_z(i)
            rji(2)  = r_r(j) - r_r(i)
            rji(3)  = r_y(j) - r_y(i)
            r2inv   = 1.00D+00/dsqrt(rji(1)*rji(1) + rji(2)*rji(2) + rji(3)*rji(3) + softening)
            r2inv   = r2inv * r2inv * r2inv * alpha
            a2_z(i) = a2_z(i) - rji(1)*r2inv
            a2_r(i) = a2_r(i) - rji(2)*r2inv
            a2_y(i) = a2_y(i) - rji(3)*r2inv
        enddo
    enddo
    !$omp end parallel
endsubroutine

subroutine update_velocities()
implicit none
integer :: im
real ( kind = 8 ), parameter :: A1 = e_eta_mass2*e_eta_mass2
real ( kind = 8 ), parameter :: A2 = e_eta_mass2*dt_S2
real ( kind = 8 ), parameter :: A3 = -e_eta_mass2*dt_S
real ( kind = 8 ) :: A4, A5
    betaL   = cte2/alphaL             ; A4 = e_eta_mass2*betaL
    alphaL  = dsqrt(cte13-betaL*betaL); A5 = e_eta_mass2*alphaL
    
    errcode=vdrnggaussian( method, stream, n_ions,Z1_r, 0.0d0, 1.0d0)
    errcode=vdrnggaussian( method, stream, n_ions,Z1_z, 0.0d0, 1.0d0)
    errcode=vdrnggaussian( method, stream, n_ions,Z1_y, 0.0d0, 1.0d0)

    !$omp parallel default(none) private(im) shared (v_z,a1_z,a2_z,v_r,a1_r,a2_r,v_y,a1_y,a2_y, Z1_r, Z0_r,Z1_y, Z0_y, Z1_z, Z0_z) firstprivate(ia,ib,A4,A5)
    im = omp_get_thread_num()+1

    v_r(ia(im):ib(im))  = A1*v_r(ia(im):ib(im)) + A2*a2_r(ia(im):ib(im)) + A3*a1_r(ia(im):ib(im)) + A4*Z0_r(ia(im):ib(im)) + A5*Z1_r(ia(im):ib(im));
    a1_r(ia(im):ib(im)) = a2_r(ia(im):ib(im)); 
    Z0_r(ia(im):ib(im)) = Z1_r(ia(im):ib(im));

    v_z(ia(im):ib(im))  = A1*v_z(ia(im):ib(im)) + A2*a2_z(ia(im):ib(im)) + A3*a1_z(ia(im):ib(im)) + A4*Z0_z(ia(im):ib(im)) + A5*Z1_z(ia(im):ib(im));
    a1_z(ia(im):ib(im)) = a2_z(ia(im):ib(im)); 
    Z0_z(ia(im):ib(im)) = Z1_z(ia(im):ib(im));

    v_y(ia(im):ib(im))  = A1*v_y(ia(im):ib(im)) + A2*a2_y(ia(im):ib(im)) + A3*a1_y(ia(im):ib(im)) + A4*Z0_y(ia(im):ib(im)) + A5*Z1_y(ia(im):ib(im));
    a1_y(ia(im):ib(im)) = a2_y(ia(im):ib(im)); 
    Z0_y(ia(im):ib(im)) = Z1_y(ia(im):ib(im));
    !$omp end parallel
end subroutine

subroutine final_time_step()
implicit none
    call update_positions();
    call update_accelerations()
    call update_velocities_final_step()
end subroutine

subroutine update_velocities_final_step()
implicit none
integer :: im
!~ beta    = b/alpha
!~ alpha   = sqrt(c-beta**2)
!~ Z1      = random.normal()
!~ v_0_12  = e_gamma2 * v_0_12 + dt*w_m*F_1 + dt*S*(F_1-F_0) + beta*Z0 + alpha*Z1
    betaL   = cte2/alphaL
!~     print*, cte3-betaL*betaL
    alphaL  = dsqrt(abs(cte3-betaL*betaL))
    errcode=vdrnggaussian( method, stream, n_ions,Z1_r, 0.0d0, 1.0d0)
    errcode=vdrnggaussian( method, stream, n_ions,Z1_z, 0.0d0, 1.0d0)
    errcode=vdrnggaussian( method, stream, n_ions,Z1_y, 0.0d0, 1.0d0)

    !$omp parallel default(none) private(im) shared (v_z,a1_z,a2_z,v_r,a1_r,a2_r, v_y,a1_y,a2_y, Z1_r, Z0_r, Z1_y, Z0_y,Z1_z, Z0_z) firstprivate(ia,ib, betaL, alphaL)
    im = omp_get_thread_num()+1
!dir$ if (integration_methode.eq.0) ! Velocity-Verlet

!dir$ elseif (integration_methode.eq.1) ! Langeving_methode
    v_r(ia(im):ib(im)) = e_eta_mass2*v_r(ia(im):ib(im)) + dt*w_m*a2_r(ia(im):ib(im)) + dt_S*(a2_r(ia(im):ib(im))-a1_r(ia(im):ib(im))) + betaL*Z0_r(ia(im):ib(im)) + alphaL*Z1_r(ia(im):ib(im));
    v_z(ia(im):ib(im)) = e_eta_mass2*v_z(ia(im):ib(im)) + dt*w_m*a2_z(ia(im):ib(im)) + dt_S*(a2_z(ia(im):ib(im))-a1_z(ia(im):ib(im))) + betaL*Z0_z(ia(im):ib(im)) + alphaL*Z1_z(ia(im):ib(im));
    v_y(ia(im):ib(im)) = e_eta_mass2*v_y(ia(im):ib(im)) + dt*w_m*a2_y(ia(im):ib(im)) + dt_S*(a2_y(ia(im):ib(im))-a1_y(ia(im):ib(im))) + betaL*Z0_y(ia(im):ib(im)) + alphaL*Z1_y(ia(im):ib(im));
!dir$ elseif (integration_methode.eq.2) ! langeving_methode: Impulse
    v_r(ia(im):ib(im)) = e_eta_mass2*v_r(ia(im):ib(im)) + dt*w_m*a2_r(ia(im):ib(im)) + betaL*Z0_r(ia(im):ib(im)) + alphaL*Z1_r(ia(im):ib(im));
    v_z(ia(im):ib(im)) = e_eta_mass2*v_z(ia(im):ib(im)) + dt*w_m*a2_z(ia(im):ib(im)) + betaL*Z0_z(ia(im):ib(im)) + alphaL*Z1_z(ia(im):ib(im));
    v_y(ia(im):ib(im)) = e_eta_mass2*v_y(ia(im):ib(im)) + dt*w_m*a2_y(ia(im):ib(im)) + betaL*Z0_y(ia(im):ib(im)) + alphaL*Z1_y(ia(im):ib(im));
!dir$ endif
    !$omp end parallel
end subroutine

subroutine initialize_ions() !Done
implicit none
double precision, parameter :: l0(3) = [50.0e-6,50.0e-6,300.0e-6];
    errcode=vdrnggaussian( method, stream, n_ions,r_r, 0.0d0, l0(1))
    errcode=vdrnggaussian( method, stream, n_ions,r_y, 0.0d0, l0(3))
    errcode=vdrnggaussian( method, stream, n_ions,r_z, 0.0d0, l0(2))


    v_z  = 0.0d0; a1_z = 0.0d0; a2_z = 0.0d0;
    v_r  = 0.0d0; a1_r = 0.0d0; a2_r = 0.0d0;
    v_y  = 0.0d0; a1_y = 0.0d0; a2_y = 0.0d0;

	t_act = 0

!~     print*, r_r(1:5)
!~     print*, r_z(1:5)
!~     stop
end subroutine

subroutine distribute_ions( )
implicit none
integer    :: i
integer                        :: nii, im
integer            :: n_by_core
integer, dimension(ni)   :: ia_tmp , ib_tmp
character(len=  10)  :: ia_str, ib_str
character(len=100)  :: str_ia, str_ib

    if (ni.ne.omp_get_max_threads ( )) then
        print*, 'ni (number of threads) is not the same in the makefile and in the code!'
        print*, ni, omp_get_max_threads ( )
        print*, 'Please, correct it, recompile and re-run'
        print*, 'Stoping code'
        stop
    endif

    i = modulo(n_ions, ni)
    n_by_core = floor(n_ions / real(ni))

    if (n_ions < ni) then
        nii = n_ions
        ia_tmp(:) = 1; ib_tmp(:) = 0; ! necessary for the ni > n_cut
        do im = 1, nii
            ia_tmp(im) = im
            ib_tmp(im) = im
        enddo
    else
        do im = 1, ni
            if (i == 0) then
                ia_tmp(im) = (im-1)*n_by_core+1
                ib_tmp(im) = ia_tmp(im) + n_by_core-1
            else
                if (im <= i) then
                    ia_tmp(im) = (im-1)*(n_by_core+1)+1
                    ib_tmp(im) = ia_tmp(im) + n_by_core
                else
                    ia_tmp(im) = n_by_core*(im - 1) + i + 1
                    ib_tmp(im) = ia_tmp(im) + n_by_core - 1
                endif
            endif
        enddo
    endif
    
    if (sum(ia - ia_tmp).ne.0) then
        print*, 'The ia, ib are not correct for the N and/or ni'
        print*, 'Update with the following: '
    
        str_ia = 'integer          , parameter :: ia(ni)=['
        str_ib = 'integer          , parameter :: ib(ni)=['
            
        do im = 1, ni
            if    (ia_tmp(im) < 10) then
                write(ia_str,"(I1.1)") ia_tmp(im)
            elseif (ia_tmp(im) < 100) then
                write(ia_str,"(I2.2)") ia_tmp(im)
            elseif (ia_tmp(im) < 1000) then
                write(ia_str,"(I3.3)") ia_tmp(im)
            else
                write(ia_str,"(I4.4)") ia_tmp(im)
            endif
            if (ib_tmp(im) < 10) then
                write(ib_str,"(I1.1)") ib_tmp(im)
            elseif (ib_tmp(im) < 100) then
                write(ib_str,"(I2.2)") ib_tmp(im)
            elseif (ib_tmp(im) < 1000) then
                write(ib_str,"(I3.3)") ib_tmp(im)
            else
                write(ib_str,"(I4.4)") ib_tmp(im)
            endif
            if (im == 1) then 
                str_ia = trim(adjustl(str_ia))//trim(adjustl(ia_str))
                str_ib = trim(adjustl(str_ib))//trim(adjustl(ib_str))
            else 
                str_ia = trim(adjustl(str_ia))//','//trim(adjustl(ia_str))
                str_ib = trim(adjustl(str_ib))//','//trim(adjustl(ib_str))
            endif
        enddo
        str_ia = trim(adjustl(str_ia))//'];'
        str_ib = trim(adjustl(str_ib))//'];'
        print*, str_ia
        print*, str_ib
        stop
    endif

endsubroutine

subroutine create_files_names( )
implicit none
character(len=10 )  :: str_N  ,   str_N_aux
character(len=130)  :: str_trap_aux
character(len=50 )  :: str_T,str_F

    write(str_N,"(I4.4)") n_ions
    str_N_aux = '_N'//trim(adjustl(str_N))

    str_trap_aux = '3D_Harmo'


    ! Temperature
    if     (Temperature .lt. 1.0E-06) then
        write(str_T,"(I3.3)") int(Temperature*1.0d9)
        str_T = '_T'//trim(adjustl(str_T))//'nK'
    elseif (Temperature .lt. 1.0E-03) then
        write(str_T,"(I3.3)") int(Temperature*1.0d6)
        str_T = '_T'//trim(adjustl(str_T))//'uK'
    elseif (Temperature .lt. 1.0E-00) then
        write(str_T,"(I3.3)") int(Temperature*1.0d3)
        str_T = '_T'//trim(adjustl(str_T))//'mK'
    endif
    ! Friction coefficient
    write(str_F,"(d9.2)") abs(eta)
    str_F = '_F'//trim(adjustl(str_F))//'Kg_s'

    ! Names used by the "from scrath" subroutine
    str_file_aux = trim(adjustl(str_trap_aux)) &
                // trim(adjustl(str_N_aux))    &
                // trim(adjustl(str_T))    &
                // trim(adjustl(str_F))    &
                // trim(adjustl(str_extra))

!~ print*, str_file_aux
!~ stop
end subroutine

subroutine save_Temp_to_file()
implicit none
character(len=130)  :: str_file_trj_tmp

    save_Temp = save_Temp*mass / (kb*n_ions)
    str_file_trj_tmp = 'Temp_'//trim(adjustl(str_file_aux))//'.bin'
    open(unit = 11, status='replace',file=trim(str_file_trj_tmp),form='unformatted')  ! create a new file, or overwrite an existing one
        write(11) n_ions
        write(11) j_save
        write(11) dt*j_save_next
        write(11) eta
        write(11) Temperature

        write(11) save_Temp(:,1:j_save)
    close(11)
end subroutine

subroutine save_xva_to_file()
implicit none
character(len=130)  :: str_file_trj_tmp

    str_file_trj_tmp = 'xva_'//trim(adjustl(str_file_aux))//'.bin'
    print*, 'Saving Final positions to:', trim(str_file_trj_tmp)
   
    open(unit = 11, status='replace',file=trim(str_file_trj_tmp),form='unformatted')  ! create a new file, or overwrite an existing one
        write(11,'(i16)') j_end ! Adrien 20200706
        write(11,'(e16.9)') t_act ! Adrien 20200706
        write(11) n_ions
        write(11) dt

        write(11) r_z
        write(11) v_z
        write(11) a1_z
        write(11) a2_z

        write(11) r_r
        write(11) v_r
        write(11) a1_r
        write(11) a2_r

        write(11) r_y
        write(11) v_y
        write(11) a1_y
        write(11) a2_y
    close(11)
end subroutine


end program main
