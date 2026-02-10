!  ************************************************
!  Main programm of 3D vesicle program
!                    written by LiJF  Dobb@bbs.fudan.edu.cn
! *************************************************  
! include 'label_system.f90'
! include 'Function_lib.f90'
! include 'update_R.f90'
! include 'Basic_parameter.f90' 
! include 'phase_separation.f90'
! include 'Energy.f90'
! include 'Bond_Flip.f90'
! include 'Forces.f90'
! include 'readfe.f90'
! include 'Vertex_average.f90'
! include 'savefiles.f90'
! include 'weep.f90'
 
 implicit none
integer N
parameter(N = 20000) 
 ! == readfe
integer N_V,N_E,N_F
real*8 rV(0:N,1:3),x,y,z,a
real*8  rV0(0:N,1:3)
integer V_E(0:N, 1:2)
integer E_F(0:N, 1:3)
! %%%%%%%%%%%%% other labels  ! == label_system
! the first kind label
integer V_F(0:N,1:3)  ! vertices around a face
integer E_V(0:N,1:20), N_E_V(0:N) ! edges ajecent to a vertice
integer F_V(0:N,1:20), N_F_V(0:N) ! star faces around a vertice
integer F_E(0:N,1:20), N_F_E(0:N)
integer F_F(0:N,1:3)

! second kind label
integer V_V(0:N,1:20), N_V_V(0:N)
integer E_V_opposite(0:N,1:20), N_E_V_opposite(0:N)
integer V_E_opposite(0:N,1:20)

integer i,j,k,i1,i2,i3,i4,i5,i0,I_temp(1:20),I_temp2(1:20),I_temp3(1:20)

! %%%%%%%%%%%%%
!Basic parameters
real*8 Len_E(0:N)  !length of edges  == get_len_E
real*8 Vector_E(0:N,1:3)  !Vector of edges == get_len_E
real*8 Vector_F(0:N,1:3), Area_F(0:N), Norm_F(0:N,1:3) 
real*8 Area_F_zero(0:N)

! == get_area_F
real*8 Vol_F(0:N), vol_total ! == get_vol_F 
real*8 H_v(0:N) ,Area_V(0:N)   ! == get_H_V
real*8 Darea_V(0:N,1:20,1:3)  ! \partial area/\partial rv == get_H_V
real*8 DVol_V(0:N,1:3)  ! \partial volume/\partial rv == get_H_V

! %%%%%%%%%%%%%%%
! Forces
real*8 P_force_V(0:N,1:3),H_force_V(0:N,1:3),A_force_V(0:N,1:3),Fi_Force_v(0:N,1:3)
real*8 delt_r
real*8 area_tot_zero,Tension_area, Con_area_Force_V(0:N,1:3)  ! Constraint_area Forces
real*8 deter_Force,Max_force
real*8 area_tot

!%%%%%  bondflip
integer I_flip,I_flip_tot
integer switch_average
! %%%%%%%%%%%%%%%%
!  energy
real*8 PV_energy, Curvature_energy,area_hooke_energy,phase_energy
real*8 EnergyP,EnergyH,energyL, energy_tot
real*8 C_energy_array(1:50000),Ph_energy_array(1:50000),energy_array(1:50000)
real*8 H_modulus, lamda, H_zero, zet_Hfi, P

! %%%%%%%%%%%%%%
!  phase separation
real*8 fi(0:N),fia0, fi_F(0:N), fi0(0:N)
real*8 L_fi,order_fi,b,a2,a4,A_phase
real*8 L_glb(0:N,1:20)

integer I_fi_F(0:N),switch_boundary_setting
real*8 len_B_array(1:50000)
integer Num_b_array(1:50000),Num_d_array(1:50000)
! for domain statistic
real*8 area_domain_files(1:50000,1:100)
real*8 ratio_area_len_files(1:50000,1:100)
real*8 len_domain_files(1:50000,1:100)
real*8 t_files(1:50000)  ! to record the time
integer k_files, switch_domain_process
real*8 area_domain0(1:1000)
integer Num_domain0, F_domain0(1:100,1:50000),N_F_domain0(1:1000)
integer No_domain0(1:1000,1:100), No_domain_Q0(1:1000)
  ! i(1:1000) for 1:Num_domain0,   return Value for actual domain ID 
  ! k(1:100) for overload domain ID, No_domain_Q0 is for overload numbers
integer inv_No_domain0(1:1000)
integer Max_No0  ! number of output files


! %%%%%%%%%%%%%%%%%%%
!  Dynamics
real*8 dt,L_r,real_t, S_opt, rv_center(1:3)
real*8 R0,amp_jig

integer t, switch_deform,tk

INTEGER seed1



!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!          END OF PARAMETER DEFININGS
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

! %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!         initial state construction
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
dt=0.01
L_r=0.25/2.0
L_fi=1.0/2.0
H_modulus=1.7
H_zero=0.0
Zet_Hfi=-0.2
lamda=0.0
amp_jig=0.0001
b=1.0; a2=1.0; a4=1.0; A_phase=1.2

fia0=0.35
Max_No0=0
switch_boundary_setting=1
switch_average=0
switch_deform=1
seed1=37815239
!seed1=67815215

call readfe(N_V,N_E,N_F,V_E,E_F,rv)
call Get_Len_E(rv,  V_E, N_E, Len_E, Vector_E)
rv_center(1)=sum(rv(1:N_v,1))/N_V
rv_center(2)=sum(rv(1:N_v,2))/N_V
rv_center(3)=sum(rv(1:N_v,3))/N_V
do j=1,N_V
rv(j,1:3)=rv(j,1:3)-rv_center
enddo

R0=23.75*(N_v/10242.0)**(1.0/2.0)
P=-H_modulus*2*(0.0-H_zero*R0)/R0**3.0


rv=R0*rv/sum(rv(1,1:3)**2.0)**0.5   ! 23.75 for 10242



call label_system(V_F,   E_V,N_E_V,     F_V,N_F_V,    F_E,N_F_E,    V_V,N_V_V,   E_V_opposite,N_E_V_opposite,   V_E_opposite,  N_V,N_E,N_F,V_E,E_F,F_F,rv)


call Get_area_F( rv, V_F,  N_F,N_V,   Area_F_zero, Vector_F, Norm_F)
call Get_area_F( rv, V_F,  N_F,N_V,   Area_F, Vector_F, Norm_F)
area_tot_zero=sum(area_F_zero(1:N_F))
call Get_Len_E(rv,  V_E, N_E, Len_E, Vector_E)
call Get_vol_F( rv, V_F,  N_F,N_V,    Vol_F, Vol_total)
call Get_H_V(rv, N_F,N_V,N_E,N_F_V, V_F, V_E, E_F,  V_V, F_V, E_V,E_V_opposite,  area_F, Vector_F, norm_F, H_V,area_V,Darea_V, Dvol_V)
call initial_phase(fia0, fi,N_V, area_V,order_fi,seed1)
 call Laplace_global(rv, area_F,area_V, N_V, V_V, F_V, N_F_V, L_glb)
write(*,*)N_v
write(*,*)area_F(1:20)
!fi=1.0
!do j=1,N_V
!if(rv(j,3).gt.0)fi(j)=-1.0
!if(abs(rv(j,3)-1.0).lt.1.5)fi(j)=0.0
!enddo
!order_fi=-0.4
!open(7,FILE="data/dynamic/f00671.dat")
!do i=1,N_V
!read(7,*)fi(i)
!enddo
!close(7)


delt_r=0.00001  !  delt_r = er^1/3  for real*8 er=10^-15

! %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!         END of initial state construction
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

! %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!         Iteration
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
fi0=fi
A_Force_V =0.0
real_t=0.0; I_flip_tot=0
tk=1
t=0
Loop_main: DO while(tk.le.49001) 



if(switch_deform.eq.1)then
! %%% calculate basic parameters
call Get_area_F( rv, V_F,  N_F,N_V,   Area_F, Vector_F, Norm_F)
call Get_vol_F( rv, V_F,  N_F,N_V,    Vol_F, Vol_total)
call Get_H_V(rv, N_F,N_V,N_E,N_F_V, V_F, V_E, E_F,  V_V, F_V, E_V,E_V_opposite,  area_F, Vector_F, norm_F, H_V,area_V,Darea_V, Dvol_V)
! %%%
! %%% calculate Forces
call Pressure_Force(P,Dvol_V,N_V,   P_Force_V)
!call Curvature_Force(Rv, V_V,N_V_V, H_modulus,H_v, H_zero, zet_Hfi, fi, N_V,  delt_r,     H_Force_v  )
call Curvature_Force_analysis(Rv, V_V,N_V_V, H_modulus,H_v, H_zero, zet_Hfi, fi, N_V,  delt_r, H_Force_v, Darea_V  )

!call Area_hooke_force(lamda, area_F,Area_F_zero,N_F, N_V, F_V, N_F_V,  Darea_V, A_Force_V )
!call phase_Force(Rv, F_V, V_V,N_V_V, fi,b, a2, a4,  N_V,  delt_r,     A_Force_v ,area_V ,area_F,A_phase)
call phase_Force_analysis(Rv, F_V, V_V,N_V_V, fi,b, a2, a4,  N_V,  delt_r,     Fi_Force_v ,area_V ,area_F,A_phase,Darea_V)


area_tot=sum(area_F(1:N_F))
call constraint_area_Force(area_tot,area_tot_zero,area_V,Darea_V,Tension_area, P_Force_V,H_Force_V,A_Force_V , Fi_Force_v ,N_V,N_F_V,F_V,L_r,dt,  Con_area_Force_V)

rv0=rv
call Laplace_global(rv, area_F,area_V, N_V, V_V, F_V, N_F_V, L_glb)
endif  
! %%%%%% end sw_deform

fi0=fi


!call update_R(rv0,rv, P_Force_V,H_Force_V,A_Force_V , Fi_Force_v , N_V, L_r,dt, area_V)
!%%%%%%%%%%%%%%%%%%%%%
! stablized the program
t=t+1


if(switch_deform.eq.1)then
Max_force=0.0
do i=1,N_V
deter_Force=sum(con_area_Force_V(i,1:3)**2)
if(deter_Force**0.5.gt.Max_Force)Max_Force=Deter_Force**0.5
if(deter_Force**0.5*L_r.ge.1000.0)then
 !con_area_Force_V(i,1:3)=0.0
switch_average=1
write(*,*)'Bigmovement',i,deter_Force**0.5*L_r

endif
enddo
if(Max_force*L_r.lt.1.0)dt=0.01
if(Max_force*L_r.ge.1.0.and.Max_force*L_r.lt.2.0)dt=0.005
if(Max_force*L_r.ge.2.0.and.Max_force*L_r.lt.10.0)dt=0.001
if(Max_force*L_r.ge.10.0.and.Max_force*L_r.lt.100.0)  dt=0.0001
if(Max_Force*L_r.ge.100.0.and.Max_force*L_r.lt.1000.0)dt=0.00001
if(Max_force*L_r.ge.1000.0)then
dt=0.000000000
t=t-1
amp_jig=0.0025
call Jiggle(rv,N_v,dt,amp_jig,seed1)
endif
rv=rv+dt*L_r*con_area_Force_V
 !!!!   update_R
endif ! end sw_deform


! Update_fi should be after Update_R
 call update_fi(fi,order_fi,L_fi,dt, H_modulus, H_V,H_zero, zet_Hfi, b, a2,a4,area_F,area_V, N_V,N_F, N_F_V,V_V, F_V,V_F, L_glb,A_phase,  rv)



! %%%%%%%%%%%%%%%%
! surface operations
if(switch_deform.eq.1)then
if(mod(t,200).eq.1)then
!%%% bond flip
call Get_Len_E(rv,  V_E, N_E, Len_E, Vector_E)
call Bond_flip(V_F,   E_V,N_E_V,     F_V,N_F_V,    F_E,N_F_E,    V_E_opposite,   N_V,N_E,N_F,V_E,E_F,rv,   Len_E,  Vector_E,I_flip)
I_flip_tot=I_flip_tot+I_flip
write(*,*)'bond_flip  ','total flips'
write(*,*)I_flip,I_flip_tot
call label_system(V_F,   E_V,N_E_V,     F_V,N_F_V,    F_E,N_F_E,    V_V,N_V_V,   E_V_opposite,N_E_V_opposite,   V_E_opposite,  N_V,N_E,N_F,V_E,E_F,F_F,rv)
call Get_area_F( rv, V_F,  N_F,N_V,   Area_F_zero, Vector_F, Norm_F)
!%%%
endif

if(mod(t,100).eq.1)then
do i=1,1
call Vertex_averaging_fi_conserve(rv,N_V,N_F,F_V,V_F,V_V,N_F_V , fi,order_fi,area_tot,area_tot_zero)
enddo
call Get_area_F( rv, V_F,  N_F,N_V,   Area_F_zero, Vector_F, Norm_F)
endif
if(mod(t,1000).eq.1)then
!call weep(V_F,   E_V,N_E_V,     F_V,N_F_V,    F_E,N_F_E,    V_V,N_V_V,   E_V_opposite,N_E_V_opposite,   V_E_opposite,  N_V,N_E,N_F,V_E,E_F,rv,fi)
do i=1,5
call Vertex_averaging_fi_conserve(rv,N_V,N_F,F_V,V_F,V_V,N_F_V , fi,order_fi,area_tot,area_tot_zero)
enddo
call Get_area_F( rv, V_F,  N_F,N_V,   Area_F_zero, Vector_F, Norm_F)
endif

do i=1,N_F
if(area_F(i).gt.2*area_F_zero(i).or.area_F(i).lt.0.5*area_F_zero(i))then
write(*,*)area_F(i),area_F_zero(i),i
switch_average=1
endif
enddo
! %%%% vertices average
if(switch_average.eq.1)then
!rv=rv0
call Get_Len_E(rv,  V_E, N_E, Len_E, Vector_E)
call Bond_flip(V_F,   E_V,N_E_V,     F_V,N_F_V,    F_E,N_F_E,    V_E_opposite,   N_V,N_E,N_F,V_E,E_F,rv,   Len_E,  Vector_E,I_flip)
I_flip_tot=I_flip_tot+I_flip
write(*,*)'bond_flip  ','total flips'
write(*,*)I_flip,I_flip_tot
call label_system(V_F,   E_V,N_E_V,     F_V,N_F_V,    F_E,N_F_E,    V_V,N_V_V,   E_V_opposite,N_E_V_opposite,   V_E_opposite,  N_V,N_E,N_F,V_E,E_F,F_F,rv)
do i=1,10
call Vertex_averaging_fi_conserve(rv,N_V,N_F,F_V,V_F,V_V,N_F_V , fi,order_fi,area_tot,area_tot_zero)
enddo
call Get_area_F( rv, V_F,  N_F,N_V,   Area_F_zero, Vector_F, Norm_F)
switch_average=0
endif

endif ! end of switch_deform

real_t=real_t+dt

! %%% stop for save files & show something
if(real_t.ge.tk)then
tk=tk+1

call Get_area_F( rv, V_F,  N_F,N_V,   Area_F, Vector_F, Norm_F)
call Get_vol_F( rv, V_F,  N_F,N_V,    Vol_F, Vol_total)
call Get_H_V(rv, N_F,N_V,N_E,N_F_V, V_F, V_E, E_F,  V_V, F_V, E_V,E_V_opposite,  area_F, Vector_F, norm_F, H_V,area_V,Darea_V, Dvol_V)
! %%%
! % calculate energies
ENERGYP=PV_energy(P,Vol_total)
energyH=Curvature_Energy(H_modulus,H_v, H_zero, zet_Hfi, fi, Area_V,N_V)
energyL=Area_Hooke_energy(lamda, area_F,Area_F_zero,N_F) !
energy_tot=energyp+energyH+energyL +phase_energy(b, a2, a4, fi, Area_V, Area_F,N_V,N_F,V_F, E_F,rv,A_phase)
! %%% calculate basic parameters
! save files output some details on the screen
call savefiles(N_V,N_E,N_F,V_E,E_F,V_F,F_E,F_F ,rv,fi,fia0,t,real_t,H_V,area_V,EnergyP,EnergyH,energyL, energy_tot,area_F,&
& area_tot,area_tot_zero,vol_total,V_E_opposite, switch_boundary_setting,dt,C_energy_array,Ph_energy_array,energy_array, & 
& Num_b_array,Num_d_array,len_B_array,tk,  switch_domain_process,area_domain_files,ratio_area_len_files,len_domain_files,& 
& t_files,k_files,  area_domain0, Num_domain0,F_domain0,N_F_domain0,No_domain0,No_domain_Q0,inv_No_domain0,Max_No0)
endif

999 format(I5.5)
if(mod(t,50000).eq.1)then
write(*,*)t
call convert_FNUM(fi, I_fi_F, fi_F, V_F, N_F)
open(7,FILE='data/R12_2.dat')
write(7,'(A8)')'vertices'
do i=1,N_V
write(7,'(I5,3(F30.15))')i,rv(i,1),rv(i,2),rv(i,3)
enddo
write(7,*)
write(7,'(A5)')'edges'
do i=1,N_E
write(7,'(3(I7),A8,I6)')i,V_E(i,1),V_E(i,2),'color',7
enddo
write(7,*)
write(7,'(A5)')'faces'
do i=1,N_F
write(7,'(4(I7), A8, I6)')i,E_F(i,1),E_F(i,2),E_F(i,3),'color',I_Fi_F(i)
enddo
write(7,*)'        '
close(7)
open(7,FILE="data/fi.dat")
do i=1,N_V
write(7,*)fi(i)
enddo
close(7)

open(7,FILE='data/Num_b.dat')
do i=1,tk
write(7,*)i,Num_b_array(i)
enddo
close(7)

endif



! %%%
END DO  Loop_main

open(7,FILE='data/Num_d.dat')
do i=1,tk
write(7,*)i,Num_d_array(i)
enddo
close(7)
open(7,FILE='data/Num_b.dat')
do i=1,tk
write(7,*)i,Num_b_array(i)
enddo
close(7)

open(7,FILE='data/Len_b.dat')
do i=1,tk
write(7,*)i,Len_b_array(i)
enddo
close(7)
open(7,FILE='data/energy.dat')
do i=1,tk
write(7,*)i,energy_array(i)
enddo
close(7)
open(7,FILE='data/energyH.dat')
do i=1,tk
write(7,*)i,C_energy_array(i)
enddo
close(7)
open(7,FILE='data/energyPH.dat')
do i=1,tk
write(7,*)i,PH_energy_array(i)
enddo
close(7)
call convert_FNUM(fi, I_fi_F, fi_F, V_F, N_F)
open(7,FILE='data/R12_2.dat')
write(7,'(A8)')'vertices'
do i=1,N_V
write(7,'(I5,3(F30.15))')i,rv(i,1),rv(i,2),rv(i,3)
enddo
write(7,*)
write(7,'(A5)')'edges'
do i=1,N_E
write(7,'(3(I7),A8,I6)')i,V_E(i,1),V_E(i,2),'color',7
enddo
write(7,*)
write(7,'(A5)')'faces'
do i=1,N_F
write(7,'(4(I7), A8, I6)')i,E_F(i,1),E_F(i,2),E_F(i,3),'color',I_Fi_F(i)
enddo
write(7,*)'        '
close(7)
open(7,FILE="data/fi.dat")
do i=1,N_V
write(7,*)fi(i)
enddo
close(7)


! %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!       End of  Iteration
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

END












