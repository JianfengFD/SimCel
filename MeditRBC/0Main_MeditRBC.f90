!  ************************************************
!  MeditRBC: Mediterranean thalassemia RBC simulation
!  Closed RBC membrane with:
!    - Isotropic Helfrich bending: kpp_pm/2*(H-H0)^2
!    - Skalak shear energy (spectrin cytoskeleton)
!    - Volume + area constraints
!    - Concentration field phi(r) (Cahn-Hilliard, conserved)
!    - Q-tensor orientation field q1,q2 (nematic, relaxational)
!    - Anisotropic coupling: [kpp_u/2*(H-H0_u)^2 + kpp_uD/2*(D-D0*cos2theta)^2]*phi
!    - K from angle deficit, D from sqrt(H^2-K)
!    - Principal directions from Taubin curvature tensor
! *************************************************

use head_Cell_mod
implicit none

! === Cell array (1 cell: the RBC membrane)
type(Cel)::Cels(1:1)
type(PtCell)::Vm

! === Dynamics
integer t, t_NV, tk, NN, N_V1, seed1
real*8 dt, dt_fi, dt_Q, f_rescale
real*8 F_total(1:3), F_norm

! === Parameters
real*8 H_modulus_in, kpp_alpha_in, mu_ms_in
real*8 b_ph_in, a2_ph_in, a4_ph_in, L_fi_in
real*8 kpp_u_in, kpp_uD_in, H0_u_in, D0_u_in
real*8 a_Q_in, a_Q4_in, K_frank_in, M_Q_in
real*8 vol_ratio

! === I/O
character*2 kc
integer i
character(len=256) :: para_file
integer :: nargs

! =========================================
! SEC. 0  Parameter setup
! =========================================
dt = 0.002
dt_fi = 0.005
dt_Q  = 0.005
f_rescale = 0.6
seed1 = 26345678
H_modulus_in = 1.0
vol_ratio = 0.6    ! reduced volume v/v_sphere ~ 0.6 for RBC

! --- Read parameters from external file
nargs = command_argument_count()
if(nargs .lt. 1) then
    write(*,*) 'Usage: ./0Main_MeditRBC para.txt'
    write(*,*) 'ERROR: parameter file not specified'
    stop
endif
call get_command_argument(1, para_file)

open(10, FILE=trim(para_file), STATUS='old', ERR=901)
read(10,*) kpp_alpha_in
read(10,*) mu_ms_in
read(10,*) b_ph_in
read(10,*) a2_ph_in
read(10,*) a4_ph_in
read(10,*) L_fi_in
read(10,*) kpp_u_in
read(10,*) kpp_uD_in
read(10,*) H0_u_in
read(10,*) D0_u_in
read(10,*) a_Q_in
read(10,*) a_Q4_in
read(10,*) K_frank_in
read(10,*) M_Q_in
close(10)

write(*,*) '--- MeditRBC: Parameters read from: ', trim(para_file)
write(*,*) '  kpp_alpha  = ', kpp_alpha_in
write(*,*) '  mu_ms      = ', mu_ms_in
write(*,*) '  b_ph       = ', b_ph_in
write(*,*) '  a2_ph      = ', a2_ph_in
write(*,*) '  a4_ph      = ', a4_ph_in
write(*,*) '  L_fi       = ', L_fi_in
write(*,*) '  kpp_u      = ', kpp_u_in
write(*,*) '  kpp_uD     = ', kpp_uD_in
write(*,*) '  H0_u       = ', H0_u_in
write(*,*) '  D0_u       = ', D0_u_in
write(*,*) '  a_Q        = ', a_Q_in
write(*,*) '  a_Q4       = ', a_Q4_in
write(*,*) '  K_frank    = ', K_frank_in
write(*,*) '  M_Q        = ', M_Q_in

goto 902
901 write(*,*) 'ERROR: cannot open parameter file: ', trim(para_file)
    stop
902 continue

! =========================================
! SEC. I  Initialization
! =========================================

! I-1: Define cell (closed RBC surface)
Cels(1)%FILE_Cel = 'R10_AD_2562.dat'
Cels(1)%N = 20000
Cels(1)%IS_OPEN = 0       ! CLOSED surface for RBC

! I-2: Read mesh size
call read_cell_size(Cels(1))
N_V1 = Cels(1)%N_V
write(*,*) 'Mesh size: N_V =', Cels(1)%N_V, '  N_E =', Cels(1)%N_E, '  N_F =', Cels(1)%N_F

! I-3: Allocate
call alloc_cel(Cels, 1, Cels(1)%N, 1)

! I-4: Read mesh topology and positions
call read_cell(Cels(1))

! I-5: Label topology
call label_cell(Cels, 1)

! I-6: For closed surface, set all boundary flags to 0
Cels(1)%IS_BC_V = 0
Cels(1)%IS_NEXTBC_V = 0
Cels(1)%IS_BC_E = 0
Cels(1)%IS_BC_F = 0

! I-7: Center the cell
call Get_Cell_Center(Cels, 1)
call Move_Cell(Cels, 1, (/0.0d0, 0.0d0, 0.0d0/))

! I-8: Compute geometry on current shape
call Get_shape_information(Cels, 1)

! I-9: Store reference geometry
Cels(1)%rv0(1:Cels(1)%N_V, 1:3) = Cels(1)%rv(1:Cels(1)%N_V, 1:3)
call Get_Len_E(Cels(1)%rv0, Cels(1)%V_E, Cels(1)%N_E, Cels(1)%Len_E_zero, Cels(1)%Vector_E)
call Get_area_F(Cels(1)%rv0, Cels(1)%V_F, Cels(1)%N_F, Cels(1)%Area_F_zero, &
& Cels(1)%Vector_F, Cels(1)%Norm_F)

! I-10: Store reference total area and volume
Cels(1)%area_tot_zero = sum(Cels(1)%Area_F_zero(1:Cels(1)%N_F))
Cels(1)%Vol_total_zero = Cels(1)%Vol_total

! I-11: Set isotropic bending parameters
Cels(1)%H_modulus = H_modulus_in
Cels(1)%H_zero = 0.0
Cels(1)%pi_K_2_A0 = 0.0
Cels(1)%integral_H_A = 0.0
Cels(1)%kpp_V = 5.0d4 * H_modulus_in / (Cels(1)%area_tot_zero / 140.0)    ! volume constraint stiffness
Cels(1)%kpp_area = 5.0d5 * H_modulus_in / (Cels(1)%area_tot_zero / 140.0) / 5.0

! I-12: Set shear (MS) parameters
Cels(1)%HAS_REFERENCE = 1
Cels(1)%kpp_alpha = kpp_alpha_in * H_modulus_in / (Cels(1)%area_tot_zero / 140.0)
Cels(1)%mu_ms = mu_ms_in * H_modulus_in / (Cels(1)%area_tot_zero / 140.0)
Cels(1)%a3_ms = -2.0
Cels(1)%a4_ms = 8.0
Cels(1)%b1_ms = 0.7
Cels(1)%b2_ms = 0.75

! I-13: Set phase field parameters
Cels(1)%b_ph = b_ph_in
Cels(1)%a2_ph = a2_ph_in
Cels(1)%a4_ph = a4_ph_in
Cels(1)%L_fi = L_fi_in

! I-14: Set anisotropic coupling parameters
Cels(1)%kpp_u = kpp_u_in
Cels(1)%kpp_uD = kpp_uD_in
Cels(1)%H0_u = H0_u_in
Cels(1)%D0_u = D0_u_in

! I-15: Set Q-tensor orientation parameters
Cels(1)%a_Q = a_Q_in
Cels(1)%a_Q4 = a_Q4_in
Cels(1)%K_frank = K_frank_in
Cels(1)%M_Q = M_Q_in

! I-16: Initialize phase field (phi ~ 0 + fluctuation)
call initial_phase_cell(Cels, 1, seed1)
write(*,*) 'Phase field initialized: fi range = ', &
& minval(Cels(1)%fi(1:Cels(1)%N_V)), maxval(Cels(1)%fi(1:Cels(1)%N_V))

! I-17: Initialize Q-tensor (small random nematic order)
call initial_Q_cell(Cels, 1, seed1)
write(*,*) 'Q-tensor initialized: q1 range = ', &
& minval(Cels(1)%q1(1:Cels(1)%N_V)), maxval(Cels(1)%q1(1:Cels(1)%N_V))
write(*,*) '                      q2 range = ', &
& minval(Cels(1)%q2(1:Cels(1)%N_V)), maxval(Cels(1)%q2(1:Cels(1)%N_V))

! I-18: Compute discrete Laplacian
call Laplace_global_Cell(Cels, 1)

! I-19: Compute curvature tensor (K, D, principal directions)
call Compute_Curvature_Tensor_Cell(Cels, 1)
write(*,*) 'Curvature tensor computed:'
write(*,*) '  K range  = ', minval(Cels(1)%K_V(1:Cels(1)%N_V)), maxval(Cels(1)%K_V(1:Cels(1)%N_V))
write(*,*) '  D range  = ', minval(Cels(1)%D_V(1:Cels(1)%N_V)), maxval(Cels(1)%D_V(1:Cels(1)%N_V))
write(*,*) '  c1 range = ', minval(Cels(1)%c1_V(1:Cels(1)%N_V)), maxval(Cels(1)%c1_V(1:Cels(1)%N_V))
write(*,*) '  c2 range = ', minval(Cels(1)%c2_V(1:Cels(1)%N_V)), maxval(Cels(1)%c2_V(1:Cels(1)%N_V))

! I-20: Restore current geometry and compute initial energy
call Get_area_F(Cels(1)%rv, Cels(1)%V_F, Cels(1)%N_F, &
& Cels(1)%Area_F, Cels(1)%Vector_F, Cels(1)%Norm_F)
call Get_Len_E(Cels(1)%rv, Cels(1)%V_E, Cels(1)%N_E, &
& Cels(1)%Len_E, Cels(1)%Vector_E)

Cels(1)%lamda = 0.0
Cels(1)%P = 0.0
call Get_Cell_Energy(Cels, 1)

write(*,*) '========================================='
write(*,*) 'Initialization complete'
write(*,*) 'N_V =', Cels(1)%N_V, '  N_E =', Cels(1)%N_E, '  N_F =', Cels(1)%N_F
write(*,*) 'Area_tot =', Cels(1)%area_tot, '  Vol_total =', Cels(1)%Vol_total
write(*,*) 'EnergyH     =', Cels(1)%EnergyH
write(*,*) 'EnergyL     =', Cels(1)%energyL
write(*,*) 'EnergyP     =', Cels(1)%energyP
write(*,*) 'EnergyMs    =', Cels(1)%energyMs
write(*,*) 'EnergyPh    =', Cels(1)%energyPh
write(*,*) 'EnergyAniso =', Cels(1)%energyAniso
write(*,*) 'EnergyOrient=', Cels(1)%energyOrient
write(*,*) 'Energy_tot  =', Cels(1)%energy_tot
write(*,*) '========================================='

! =========================================
! SEC. II  Dynamics
! =========================================

tk = 0
NN = 30000

Loop_main: DO t = 1, NN

    ! === Inner loop: vertex dynamics ===
    do t_NV = 1, Cels(1)%N_V

        ! --- Random vertex selection
        Vm%V_move = mod(int(random1(seed1) * Cels(1)%N_V), Cels(1)%N_V) + 1
        Vm%rv_move = Cels(1)%rv(Vm%V_move, 1:3)

        ! --- Update constraints
        Cels(1)%lamda = Cels(1)%kpp_area * (Cels(1)%area_tot - Cels(1)%area_tot_zero) &
        & / Cels(1)%area_tot_zero
        Cels(1)%P = Cels(1)%kpp_V * (Cels(1)%Vol_total - Cels(1)%Vol_total_zero) &
        & / Cels(1)%Vol_total_zero

        ! --- Isotropic bending force (numerical FD)
        Vm%PM_Force_Point = 0.0
        call Point_H_Force_Cell(Cels, 1, Vm, Cels(1)%integral_H_A)

        ! --- Volume + area constraint force
        Vm%PV_Force_Point = 0.0
        call Point_PV_Force_Cell(Cels, 1, Vm)

        ! --- Shear (MS) force
        Vm%MS_Force_Point = 0.0
        call Point_MS_Force_Cell(Cels, 1, Vm)

        ! --- Anisotropic coupling force (numerical FD)
        Vm%Aniso_Force_Point = 0.0
        call Point_Aniso_Force_Cell(Cels, 1, Vm)

        ! --- Total force and stability guard
        F_total = f_rescale * (Vm%PM_Force_Point + Vm%PV_Force_Point &
            + Vm%MS_Force_Point + Vm%Aniso_Force_Point)
        F_norm = sqrt(sum(F_total**2))

        if(F_norm * dt .lt. 0.05)then
            Cels(1)%rv(Vm%V_move, 1:3) = Vm%rv_move + dt * F_total
        endif

        ! --- Local geometry update
        Vm%rv_move = Cels(1)%rv(Vm%V_move, 1:3)
        call Point_Update_Cell(Cels, 1, Vm)

    enddo ! End of t_NV

    ! === After inner loop ===

    ! Update Laplacian weights
    call Laplace_global_Cell(Cels, 1)

    ! Update curvature tensor (K, D, directions)
    call Compute_Curvature_Tensor_Cell(Cels, 1)

    ! Evolve Q-tensor (orientation relaxation)
    call update_Q_cell(Cels, 1, dt_Q)

    ! Evolve phi (Cahn-Hilliard)
    call update_fi_cell(Cels, 1, dt_fi)

    ! === Save shape + fields every 200 steps ===
    if(mod(t, 200).eq.0)then
        tk = tk + 1

        ! Full geometry recomputation
        call Get_shape_information(Cels, 1)
        call Compute_Curvature_Tensor_Cell(Cels, 1)

        ! Compute energies
        Cels(1)%lamda = Cels(1)%kpp_area * (Cels(1)%area_tot - Cels(1)%area_tot_zero) &
        & / Cels(1)%area_tot_zero
        Cels(1)%P = Cels(1)%kpp_V * (Cels(1)%Vol_total - Cels(1)%Vol_total_zero) &
        & / Cels(1)%Vol_total_zero
        call Get_Cell_Energy(Cels, 1)

        write(*,'(A,I9,A,I6)') ' Step ', t, '   tk = ', tk
        write(*,'(A,F14.6,A,F14.6)') '  area_tot = ', Cels(1)%area_tot, &
        & '  vol_total = ', Cels(1)%Vol_total
        write(*,'(A,F12.6,A,F12.6,A,F12.6)') '  EnergyH = ', Cels(1)%EnergyH, &
        & '  EnergyL = ', Cels(1)%energyL, '  EnergyP = ', Cels(1)%energyP
        write(*,'(A,F12.6,A,F12.6)') '  EnergyMs = ', Cels(1)%energyMs, &
        & '  EnergyPh = ', Cels(1)%energyPh
        write(*,'(A,F12.6,A,F12.6)') '  EnergyAniso = ', Cels(1)%energyAniso, &
        & '  EnergyOrient = ', Cels(1)%energyOrient
        write(*,'(A,F12.6)') '  Energy_tot = ', Cels(1)%energy_tot
        write(*,'(A,F10.6,A,F10.6)') '  fi_min = ', minval(Cels(1)%fi(1:Cels(1)%N_V)), &
        & '  fi_max = ', maxval(Cels(1)%fi(1:Cels(1)%N_V))
        write(*,'(A,F10.6,A,F10.6)') '  S_min  = ', &
        & minval(sqrt(Cels(1)%q1(1:Cels(1)%N_V)**2+Cels(1)%q2(1:Cels(1)%N_V)**2)), &
        & '  S_max  = ', &
        & maxval(sqrt(Cels(1)%q1(1:Cels(1)%N_V)**2+Cels(1)%q2(1:Cels(1)%N_V)**2))

        ! Save shape snapshot
        call save_Cell(Cels(1), tk, 0, 'MRC')

        ! Save phi field
        write(kc,'(I2.2)') 0
        open(17, FILE='data/Rfi'//kc//'.dat')
        write(17,*) '// fi at step ', t, ' tk =', tk
        do i = 1, Cels(1)%N_V
            write(17,'(I6,F15.8)') i, Cels(1)%fi(i)
        enddo
        close(17)

        ! Save Q-tensor field
        open(18, FILE='data/RQ'//kc//'.dat')
        write(18,*) '// q1,q2,S at step ', t, ' tk =', tk
        do i = 1, Cels(1)%N_V
            write(18,'(I6,4F15.8)') i, Cels(1)%q1(i), Cels(1)%q2(i), &
                sqrt(Cels(1)%q1(i)**2 + Cels(1)%q2(i)**2), &
                Cels(1)%D_V(i)
        enddo
        close(18)
    endif

END DO Loop_main

write(*,*) 'MeditRBC simulation complete'

END
