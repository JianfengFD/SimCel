!  ************************************************
!  Flat surface with shear energy simulation
!  Open membrane, static mesh, reference shape shrinking
!  Energy: H_modulus*H^2/2 + lam/2*(A-A0)^2/A0 + shear(MS)
!  pi_K_2_A0 = 0 (no integrated curvature term)
!  No volume energy, no vertex averaging, no surf operations
!  Reference shape shrinks 1/1000 in x,y every 100 steps
! *************************************************

use head_Cell_mod
implicit none

! === Cell array (1 cell: the membrane)
type(Cel)::Cels(1:1)
type(PtCell)::Vm

! === Dynamics
integer t, tk, N_V1, seed1
real*8 dt, f_rescale
real*8 F_total(1:3), F_norm

! === Parameters
real*8 H_modulus_in, kpp_alpha_in, mu_ms_in
real*8 lam_base

! === I/O (unused vars removed)

! =========================================
! SEC. 0  Parameter setup
! =========================================
dt = 0.005
f_rescale = 0.8
seed1 = 26345678
H_modulus_in = 1.0
kpp_alpha_in = 25.0
mu_ms_in = 12.5   ! = 0.5 * kpp_alpha

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

! I-8: Store reference geometry from current shape
! Len_E_zero and Area_F_zero are computed from rv0 (= copy of rv)
Cels(1)%rv0(1:Cels(1)%N_V, 1:3) = Cels(1)%rv(1:Cels(1)%N_V, 1:3)
call Get_Len_E(Cels(1)%rv0, Cels(1)%V_E, Cels(1)%N_E, Cels(1)%Len_E_zero, Cels(1)%Vector_E)
call Get_area_F(Cels(1)%rv0, Cels(1)%V_F, Cels(1)%N_F, Cels(1)%Area_F_zero, &
& Cels(1)%Vector_F, Cels(1)%Norm_F)

! I-9: Store reference total area
Cels(1)%area_tot_zero = sum(Cels(1)%Area_F_zero(1:Cels(1)%N_F))
Cels(1)%Vol_total_zero = 1.0  ! not used but avoid division by zero

! I-10: Set material properties
Cels(1)%H_modulus = H_modulus_in
Cels(1)%H_zero = 0.0
Cels(1)%pi_K_2_A0 = 0.0    ! no integrated curvature energy
Cels(1)%integral_H_A = 0.0
Cels(1)%kpp_V = 0.0          ! no volume constraint
Cels(1)%kpp_area = 5.0d5 * H_modulus_in / (Cels(1)%area_tot_zero / 140.0) / 5.0

! I-11: Set shear (MS) parameters
Cels(1)%HAS_REFERENCE = 1
Cels(1)%kpp_alpha = kpp_alpha_in * H_modulus_in / (Cels(1)%area_tot_zero / 140.0)
Cels(1)%mu_ms = mu_ms_in * H_modulus_in / (Cels(1)%area_tot_zero / 140.0)
Cels(1)%a3_ms = -2.0
Cels(1)%a4_ms = 8.0
Cels(1)%b1_ms = 0.7
Cels(1)%b2_ms = 0.75

! I-12: Initial energy computation
Cels(1)%lamda = 0.0
Cels(1)%P = 0.0
call Get_Cell_Energy(Cels, 1)

write(*,*) 'Initialization complete'
write(*,*) 'N_V =', Cels(1)%N_V, '  N_E =', Cels(1)%N_E, '  N_F =', Cels(1)%N_F
write(*,*) 'Area_tot =', Cels(1)%area_tot, '  Area_tot_zero =', Cels(1)%area_tot_zero
write(*,*) 'EnergyH =', Cels(1)%EnergyH, '  EnergyL =', Cels(1)%energyL, &
& '  EnergyMs =', Cels(1)%energyMs
write(*,*) 'Energy_tot =', Cels(1)%energy_tot

! =========================================
! SEC. II  Dynamics
! =========================================

tk = 0

Loop_main: DO t = 1, 9000000

    ! --- Random vertex selection (skip boundary)
    Vm%V_move = mod(int(random1(seed1) * Cels(1)%N_V), Cels(1)%N_V) + 1
    if(Cels(1)%IS_BC_V(Vm%V_move).eq.1) cycle Loop_main
    if(Cels(1)%IS_NEXTBC_V(Vm%V_move).eq.1) cycle Loop_main

    Vm%rv_move = Cels(1)%rv(Vm%V_move, 1:3)

    ! --- Update lambda (area constraint tension)
    Cels(1)%lamda = Cels(1)%kpp_area * (Cels(1)%area_tot - Cels(1)%area_tot_zero) &
    & / Cels(1)%area_tot_zero

    ! --- Compute forces
    ! Bending force (numerical finite difference)
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

    ! --- Local geometry update (no full remesh for static mesh)
    Vm%rv_move = Cels(1)%rv(Vm%V_move, 1:3)
    call Point_Update_Cell(Cels, 1, Vm)

    ! --- Update current edge lengths (needed for shear force next step)
    ! This is done locally in Update_New_Parameters already for edges around V_move

    ! --- Reference shape shrinking: every 100 steps, shrink x,y by 1/1000
    if(mod(t, 100).eq.0)then
        Cels(1)%rv0(1:Cels(1)%N_V, 1) = Cels(1)%rv0(1:Cels(1)%N_V, 1) * (1.0 - 1.0/1000.0)
        Cels(1)%rv0(1:Cels(1)%N_V, 2) = Cels(1)%rv0(1:Cels(1)%N_V, 2) * (1.0 - 1.0/1000.0)
        ! Recompute reference edge lengths and face areas
        call Get_Len_E(Cels(1)%rv0, Cels(1)%V_E, Cels(1)%N_E, &
        & Cels(1)%Len_E_zero, Cels(1)%Vector_E)
        call Get_area_F(Cels(1)%rv0, Cels(1)%V_F, Cels(1)%N_F, &
        & Cels(1)%Area_F_zero, Cels(1)%Vector_F, Cels(1)%Norm_F)
        ! Restore current geometry (Get_area_F overwrote Vector_F/Norm_F)
        call Get_area_F(Cels(1)%rv, Cels(1)%V_F, Cels(1)%N_F, &
        & Cels(1)%Area_F, Cels(1)%Vector_F, Cels(1)%Norm_F)
        call Get_Len_E(Cels(1)%rv, Cels(1)%V_E, Cels(1)%N_E, &
        & Cels(1)%Len_E, Cels(1)%Vector_E)
    endif

    ! --- Save shape every 200 steps
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
        write(*,'(A,F12.6)') '  Energy_tot = ', Cels(1)%energy_tot

        ! Save shape snapshot
        call save_Cell(Cels(1), tk, 0, 'FLT')
    endif

END DO Loop_main

write(*,*) 'Simulation complete'

END
