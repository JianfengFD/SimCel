!  ************************************************
!  Flat surface with AB phase separation
!  Open membrane, static mesh, reference shape shrinking
!  Bending energy: (kpp + kpp1*phi)/2 * (H + c0 - c1*phi)^2
!  Phase: Cahn-Hilliard with GL double-well + interface gradient
!  No dynamic mesh operations
! *************************************************

use head_Cell_mod
implicit none

! === Cell array (1 cell: the membrane)
type(Cel)::Cels(1:1)
type(PtCell)::Vm

! === Dynamics
integer t, t_NV, tk, NN, N_V1, seed1
real*8 dt, dt_fi, f_rescale
real*8 F_total(1:3), F_norm

! === Parameters
real*8 H_modulus_in, kpp_alpha_in, mu_ms_in
real*8 lam_base
real*8 avg_len_E, scale_fac

! === Phase field parameters
real*8 b_ph_in, a2_ph_in, a4_ph_in
real*8 c0_phi_in, c1_phi_in, kpp1_phi_in
real*8 L_fi_in

! === I/O
character*2 kc
integer i
character(len=256) :: para_file
integer :: nargs

! =========================================
! SEC. 0  Parameter setup
! =========================================
dt = 0.005
dt_fi = 0.01
f_rescale = 0.8
seed1 = 26345678
H_modulus_in = 1.0

! --- Read parameters from external file (command line argument)
nargs = command_argument_count()
if(nargs .lt. 1) then
    write(*,*) 'Usage: ./0Main_FlatSurfaceAB para.txt'
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
read(10,*) kpp1_phi_in
read(10,*) c0_phi_in
read(10,*) c1_phi_in
close(10)

write(*,*) '--- Parameters read from: ', trim(para_file)
write(*,*) '  kpp_alpha  = ', kpp_alpha_in
write(*,*) '  mu_ms      = ', mu_ms_in
write(*,*) '  b_ph       = ', b_ph_in
write(*,*) '  a2_ph      = ', a2_ph_in
write(*,*) '  a4_ph      = ', a4_ph_in
write(*,*) '  L_fi       = ', L_fi_in
write(*,*) '  kpp1_phi   = ', kpp1_phi_in
write(*,*) '  c0_phi     = ', c0_phi_in
write(*,*) '  c1_phi     = ', c1_phi_in

goto 902
901 write(*,*) 'ERROR: cannot open parameter file: ', trim(para_file)
    stop
902 continue

! =========================================
! SEC. I  Initialization
! =========================================

! I-1: Define cell
Cels(1)%FILE_Cel = 'PLN2.dat'
Cels(1)%N = 20000
Cels(1)%IS_OPEN = 1

! I-2: Read mesh size
call read_cell_size(Cels(1))
N_V1 = Cels(1)%N_V

! I-3: Allocate
call alloc_cel(Cels, 1, Cels(1)%N, 1)

! I-4: Read mesh topology and positions
call read_cell(Cels(1))

! I-5: Label topology
call label_cell(Cels, 1)

! I-6: Center the cell
call Get_Cell_Center(Cels, 1)
call Move_Cell(Cels, 1, (/0.0d0, 0.0d0, 0.0d0/))

! I-7: Compute geometry on current shape
call Get_shape_information(Cels, 1)

! I-7b: Rescale so that average Len_E = 1
avg_len_E = sum(Cels(1)%Len_E(1:Cels(1)%N_E)) / Cels(1)%N_E
scale_fac = 1.0d0 / avg_len_E
write(*,*) 'Before rescale: avg_len_E =', avg_len_E, '  scale_fac =', scale_fac
Cels(1)%rv(1:Cels(1)%N_V, 1:3) = Cels(1)%rv(1:Cels(1)%N_V, 1:3) * scale_fac
call Get_shape_information(Cels, 1)
write(*,*) 'After  rescale: avg_len_E =', &
& sum(Cels(1)%Len_E(1:Cels(1)%N_E)) / Cels(1)%N_E

! I-8: Store reference geometry from current shape
Cels(1)%rv0(1:Cels(1)%N_V, 1:3) = Cels(1)%rv(1:Cels(1)%N_V, 1:3)
call Get_Len_E(Cels(1)%rv0, Cels(1)%V_E, Cels(1)%N_E, Cels(1)%Len_E_zero, Cels(1)%Vector_E)
call Get_area_F(Cels(1)%rv0, Cels(1)%V_F, Cels(1)%N_F, Cels(1)%Area_F_zero, &
& Cels(1)%Vector_F, Cels(1)%Norm_F)

! I-9: Store reference total area
Cels(1)%area_tot_zero = sum(Cels(1)%Area_F_zero(1:Cels(1)%N_F))
Cels(1)%Vol_total_zero = 1.0

! I-10: Set material properties
Cels(1)%H_modulus = H_modulus_in
Cels(1)%H_zero = 0.0
Cels(1)%pi_K_2_A0 = 0.0
Cels(1)%integral_H_A = 0.0
Cels(1)%kpp_V = 0.0
Cels(1)%kpp_area = 5.0d5 * H_modulus_in / (Cels(1)%area_tot_zero / 140.0) / 5.0

! I-11: Set shear (MS) parameters
Cels(1)%HAS_REFERENCE = 1
Cels(1)%kpp_alpha = kpp_alpha_in * H_modulus_in / (Cels(1)%area_tot_zero / 140.0)
Cels(1)%mu_ms = mu_ms_in * H_modulus_in / (Cels(1)%area_tot_zero / 140.0)
Cels(1)%a3_ms = -2.0
Cels(1)%a4_ms = 8.0
Cels(1)%b1_ms = 0.7
Cels(1)%b2_ms = 0.75

! I-12: Set phase field parameters
Cels(1)%b_ph = b_ph_in
Cels(1)%a2_ph = a2_ph_in
Cels(1)%a4_ph = a4_ph_in
Cels(1)%L_fi = L_fi_in
Cels(1)%kpp1_phi = kpp1_phi_in
Cels(1)%c0_phi = c0_phi_in
Cels(1)%c1_phi = c1_phi_in

! I-13: Initialize phase field (fi ~ 0 with 0.1 fluctuation)
call initial_phase_cell(Cels, 1, seed1)
write(*,*) 'Phase field initialized: fi range = ', &
& minval(Cels(1)%fi(1:Cels(1)%N_V)), maxval(Cels(1)%fi(1:Cels(1)%N_V))

! I-14: Compute discrete Laplacian
call Laplace_global_Cell(Cels, 1)

! I-15: Restore current geometry (may have been overwritten)
call Get_area_F(Cels(1)%rv, Cels(1)%V_F, Cels(1)%N_F, &
& Cels(1)%Area_F, Cels(1)%Vector_F, Cels(1)%Norm_F)
call Get_Len_E(Cels(1)%rv, Cels(1)%V_E, Cels(1)%N_E, &
& Cels(1)%Len_E, Cels(1)%Vector_E)

! I-16: Initial energy computation
Cels(1)%lamda = 0.0
Cels(1)%P = 0.0
call Get_Cell_Energy(Cels, 1)

write(*,*) 'Initialization complete'
write(*,*) 'N_V =', Cels(1)%N_V, '  N_E =', Cels(1)%N_E, '  N_F =', Cels(1)%N_F
write(*,*) 'Area_tot =', Cels(1)%area_tot, '  Area_tot_zero =', Cels(1)%area_tot_zero
write(*,*) 'EnergyH =', Cels(1)%EnergyH, '  EnergyL =', Cels(1)%energyL, &
& '  EnergyMs =', Cels(1)%energyMs
write(*,*) 'EnergyPh =', Cels(1)%energyPh
write(*,*) 'Energy_tot =', Cels(1)%energy_tot

! =========================================
! SEC. II  Dynamics
! =========================================

tk = 0
NN = 45000

Loop_main: DO t = 1, NN

    do t_NV = 1, Cels(1)%N_V

        ! --- Random vertex selection (skip boundary)
        Vm%V_move = mod(int(random1(seed1) * Cels(1)%N_V), Cels(1)%N_V) + 1
        if(Cels(1)%IS_BC_V(Vm%V_move).eq.1) cycle
        if(Cels(1)%IS_NEXTBC_V(Vm%V_move).eq.1) cycle

        Vm%rv_move = Cels(1)%rv(Vm%V_move, 1:3)

        ! --- Update lambda (area constraint tension)
        Cels(1)%lamda = Cels(1)%kpp_area * (Cels(1)%area_tot - Cels(1)%area_tot_zero) &
        & / Cels(1)%area_tot_zero

        ! --- Compute forces
        ! Bending force with phi coupling (numerical finite difference)
        Vm%PM_Force_Point = 0.0
        call Point_H_Force_Cell(Cels, 1, Vm, Cels(1)%integral_H_A)

        ! Area constraint force (analytical, using PV with P=0)
        Vm%PV_Force_Point = 0.0
        Cels(1)%P = 0.0
        call Point_PV_Force_Cell(Cels, 1, Vm)

        ! Shear (MS) force (analytical)
        Vm%MS_Force_Point = 0.0
        call Point_MS_Force_Cell(Cels, 1, Vm)

        ! --- Total force and stability guard
        F_total = f_rescale * (Vm%PM_Force_Point + Vm%PV_Force_Point + Vm%MS_Force_Point)
        F_norm = sum(F_total**2)**0.5

        if(F_norm * dt .lt. 0.05)then
            Cels(1)%rv(Vm%V_move, 1:3) = Vm%rv_move + dt * F_total
        endif

        ! --- Local geometry update
        Vm%rv_move = Cels(1)%rv(Vm%V_move, 1:3)
        call Point_Update_Cell(Cels, 1, Vm)

    enddo ! End of t_NV

    ! --- After inner loop: update Laplacian and run CH dynamics
    call Laplace_global_Cell(Cels, 1)
    call update_fi_cell(Cels, 1, dt_fi)

    ! --- Reference shape shrinking: every 100 steps, shrink x,y by 1/1000
    if(mod(t, 100).eq.0)then
        Cels(1)%rv0(1:Cels(1)%N_V, 1) = Cels(1)%rv0(1:Cels(1)%N_V, 1) * (1.0 - 1.0/1000.0)
        Cels(1)%rv0(1:Cels(1)%N_V, 2) = Cels(1)%rv0(1:Cels(1)%N_V, 2) * (1.0 - 1.0/1000.0)
        call Get_Len_E(Cels(1)%rv0, Cels(1)%V_E, Cels(1)%N_E, &
        & Cels(1)%Len_E_zero, Cels(1)%Vector_E)
        call Get_area_F(Cels(1)%rv0, Cels(1)%V_F, Cels(1)%N_F, &
        & Cels(1)%Area_F_zero, Cels(1)%Vector_F, Cels(1)%Norm_F)
        call Get_area_F(Cels(1)%rv, Cels(1)%V_F, Cels(1)%N_F, &
        & Cels(1)%Area_F, Cels(1)%Vector_F, Cels(1)%Norm_F)
        call Get_Len_E(Cels(1)%rv, Cels(1)%V_E, Cels(1)%N_E, &
        & Cels(1)%Len_E, Cels(1)%Vector_E)
    endif

    ! --- Save shape + concentration every 200 steps
    if(mod(t, 200).eq.0)then
        tk = tk + 1

        ! Full geometry recomputation
        call Get_shape_information(Cels, 1)

        ! Compute energies
        Cels(1)%lamda = Cels(1)%kpp_area * (Cels(1)%area_tot - Cels(1)%area_tot_zero) &
        & / Cels(1)%area_tot_zero
        Cels(1)%P = 0.0
        call Get_Cell_Energy(Cels, 1)

        write(*,'(A,I9,A,I6)') ' Step ', t, '   tk = ', tk
        write(*,'(A,F12.6,A,F12.6)') '  area_tot = ', Cels(1)%area_tot, &
        & '  area_zero = ', Cels(1)%area_tot_zero
        write(*,'(A,F12.6,A,F12.6,A,F12.6)') '  EnergyH = ', Cels(1)%EnergyH, &
        & '  EnergyL = ', Cels(1)%energyL, '  EnergyMs = ', Cels(1)%energyMs
        write(*,'(A,F12.6)') '  EnergyPh = ', Cels(1)%energyPh
        write(*,'(A,F12.6)') '  Energy_tot = ', Cels(1)%energy_tot
        write(*,'(A,F10.6,A,F10.6)') '  fi_min = ', minval(Cels(1)%fi(1:Cels(1)%N_V)), &
        & '  fi_max = ', maxval(Cels(1)%fi(1:Cels(1)%N_V))

        ! Save shape snapshot
        call save_Cell(Cels(1), tk, 0, 'FAB')

        ! Save concentration (fi) to separate file
        write(kc,'(I2.2)') 0
        open(17, FILE='data/Rfi'//kc//'.dat')
        write(17,*) '// fi at step ', t, ' tk =', tk
        do i = 1, Cels(1)%N_V
            write(17,'(I6,F15.8)') i, Cels(1)%fi(i)
        enddo
        close(17)
    endif

END DO Loop_main

write(*,*) 'Simulation complete'

END
