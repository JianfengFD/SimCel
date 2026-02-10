!  ************************************************
!  Main programm of 3D vesicle program
!                    written by LiJF  Dobb@bbs.fudan.edu.cn
! *************************************************  

!  One point movement at one time 


!include 'RBC_para_Fun.f90'
!include 'RBC_label.f90'
!include 'RBC_Force.f90'
!include 'RBC_datainput.f90'
!include 'RBC_Cap.f90'

implicit none

integer N
parameter(N = 50000) 
 ! == readfe
integer N_V,N_E,N_F
real*8  rV(0:N,1:3)
real*8  rV0(0:N,1:3)
integer V_E(0:N, 1:2)
integer E_F(0:N, 1:3)

! %%%%%%%%%%%%% other labels  ! == label_system
! the first kind label
integer V_F(0:N,1:3)  ! vertices around a face
integer E_V(0:N,1:20), N_E_V(0:N) ! edges adjacent to a vertice
integer F_V(0:N,1:20), N_F_V(0:N) ! star faces around a vertice
integer F_E(0:N,1:20), N_F_E(0:N)
integer F_F(0:N,1:3)

! second kind label
integer V_V(0:N,1:20), N_V_V(0:N)
integer E_V_opposite(0:N,1:20), N_E_V_opposite(0:N)
integer V_E_opposite(0:N,1:20)

integer i,j,i1,j1

! %%%%%%%%%%%%%
!Basic parameters
real*8 Len_E(0:N) !length of edges  == get_len_E
real*8 th_E(0:N)
real*8 Len_E_zero(0:N)
real*8 Vector_E(0:N,1:3)  !Vector of edges == get_len_E
real*8 Vector_F(0:N,1:3), Area_F(0:N), Norm_F(0:N,1:3) 
real*8 Norm_V(0:N,1:3)   ! normal vector of a vertex

real*8 Area_F_zero(0:N)

! == get_area_F
real*8 Vol_F(0:N), vol_total ! == get_vol_F 
real*8 Vol_total_zero
real*8 H_v(0:N) ,Area_V(0:N)   ! == get_H_V
real*8 area_V_zero(0:N)
real*8 Darea_V(0:N,1:20,1:3)  ! \partial area/\partial rv == get_H_V
real*8 DVol_V(0:N,1:3)  ! \partial volume/\partial rv == get_H_V


real*8 delt_r
real*8 area_tot_zero
real*8 area_tot

! %%%%%%%%%%%%%%%%
!  energy
real*8 EnergyP,EnergyH,energyL,energyMs, energy_tot
real*8 H_modulus, lamda, H_zero, P,Curvature_Energy


! %%%%%%%%%%%%%%%%%%%
!  Dynamics
real*8 dt,dt0,real_t, rv_center(1:3)
real*8 R0

integer t,tk


real*8 kpp_V
real*8 pi_K_2_A0
real*8 kpp_alpha,mu_ms
real*8 kpp_Area
real*8 a3_ms,a4_ms,b1_ms,b2_ms

real*8 integral_H_A



 !   %%%%%%% one point movement

 real*8 rv_move(1:3)
 real*8 PV_Force_Point(1:3), PM_Force_Point(1:3), MS_Force_Point(1:3)

real*8 Darea_Point(1:20,1:3),Dvol_Point(1:3)
real*8 kBT, gauss_kbt,inverse_kbt

integer V_move
integer K_move(0:N)

real*8 random1,x1,y1,z1
integer seed1
! %%%%%%%%%%%%%%%%%%
! capp
real*8 R_cap,Cap_Force_point(1:3)

real*8 th,ph,ps
real*8 tao_stress,tao_input, d_cap_input

real*8 contact_area,contact_pressure
real*8 contact_P(1:100000)
integer contact_V(0:N),color_F(0:N),Num_contact

! %%% avoiding intersect
real*8 A_LJ,B_LJ,sig_LJ,cutoff_LJ,epsl_LJ
real*8 AMM_LJ,BMM_LJ,sigMM_LJ
real*8 cut_off,rep_const, rep_force_V_Point(1:3)

integer Num_near_V(0:N),Near_V(0:N,1:30)


character*3 tao_char,d_cap_char
character*3 file_input,edgeon_axial

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!          END OF PARAMETER DEFININGS
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

! %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!         initial state construction
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
i=1
!call getarg(i,tao_char)
i=2
!call getarg(i,d_cap_char)
i=3
!call getarg(i,edgeon_axial)
i=4
!call getarg(i,file_input)

! Also see "Rotate RBC and MS"
edgeon_axial='edg'
if(edgeon_axial.eq.'edg')th =0*3.14159265/2.0
if(edgeon_axial.eq.'axi')th =1*3.14159265/2.0

tao_char = '0.0'
d_cap_char = '5.8'
file_input = '001'
call char2real_No_E(len(trim(tao_char)),trim(tao_char),tao_input)
call char2real_No_E(len(trim(d_cap_char)),trim(d_cap_char),d_cap_input)




dt=0.0025

dt0=dt
seed1=26479331


call readfe(N_V,N_E,N_F,V_E,E_F,rv)
call readfe_MS(N_V,N_E,N_F,V_E,E_F,rv0)


rv_center(1)=sum(rv(1:N_v,1))/N_V
rv_center(2)=sum(rv(1:N_v,2))/N_V
rv_center(3)=sum(rv(1:N_v,3))/N_V

do j=1,N_V
rv(j,1:3)=rv(j,1:3)-rv_center
enddo
rv_center(1)=sum(rv0(1:N_v,1))/N_V
rv_center(2)=sum(rv0(1:N_v,2))/N_V
rv_center(3)=sum(rv0(1:N_v,3))/N_V
do j=1,N_V
rv0(j,1:3)=rv0(j,1:3)-rv_center
enddo


!%%%%%%%% rotate RBC and MS
ph=0;ps=0;
call  rotate_RBC(rv,N_V,th,ph,ps)
call  rotate_RBC(rv0,N_V,th,ph,ps)

call label_system(V_F,   E_V,N_E_V,     F_V,N_F_V,    F_E,N_F_E,    V_V,N_V_V,  &
& E_V_opposite,N_E_V_opposite,   V_E_opposite,  N_V,N_E,N_F,V_E,E_F,F_F)
call Get_area_F( rv0, V_F,  N_F,   Area_F_zero, Vector_F, Norm_F)
area_tot_zero=sum(area_F_zero(1:N_F))


R0=  23.75*(N_v/10242.0)**(1.0/2.0)  ! 23.75 for 10242(140.0/4/3.14159265)**0.5
rv0=R0*rv0/(area_tot_zero/4.0/3.14159265)**0.5

call Get_area_F( rv, V_F,  N_F,   Area_F_zero, Vector_F, Norm_F)
area_tot_zero=sum(area_F_zero(1:N_F))
rv=R0*rv/(area_tot_zero/4.0/3.14159265)**0.5

call Get_area_F( rv0, V_F,  N_F,   Area_F_zero, Vector_F, Norm_F)
call Get_area_F( rv, V_F,  N_F,   Area_F, Vector_F, Norm_F)
area_tot_zero=sum(area_F_zero(1:N_F))
area_tot=sum(area_F(1:N_F))

call Get_Len_E(rv0,  V_E, N_E, Len_E_zero, Vector_E)
call Get_vol_F( rv0, V_F,  N_F,    Vol_F, Vol_total_zero)
call Get_vol_F( rv, V_F,  N_F,    Vol_F, Vol_total)
call Get_H_V_new(rv0,N_V,N_F_V, V_V, F_V, E_V,&
& area_F_zero, Vector_F, H_V,area_V_zero,Darea_V, Dvol_V,th_E)
call Get_H_V_new(rv,N_V,N_F_V, V_V, F_V, E_V,&
&  area_F, Vector_F, H_V,area_V,Darea_V, Dvol_V,th_E)


vol_total_zero=(area_tot_zero/140.0)**1.5*100.0

! ms,pm parameters setting

H_modulus=1.0                     ! Kpp0 = 0.2 pN*um
inverse_kBT=20000.0
kBT=H_modulus/inverse_kBT   !/100.0


kpp_V =   5.0*10.0**4.0*H_modulus/(vol_total_zero/100.0)  
kpp_area = 5.0*10.0**4.0*H_modulus/(area_tot_zero/140.0)
kpp_alpha= 25.0*H_modulus/(area_tot_zero/140.0)



H_zero=-10.0/3.14159265/((area_tot_zero/4.0/3.14159265)**0.5)  !/R0  AD -10
pi_K_2_A0= H_modulus/area_tot_zero

mu_ms=0.5*kpp_alpha
a3_ms=-2.0
a4_ms=8.0
b1_ms=0.7
b2_ms=0.75

tao_stress = (tao_input /0.2) * ( 100/vol_total_zero) *H_modulus   ! tao = 2.0 pN/um^2

i=2
call Max_XYZ(rv,N_V,i,R_cap)


R_cap=R_cap+0.2
                                         ! %%  size of the cap
call Get_Len_E(rv,  V_E, N_E, Len_E, Vector_E)

epsl_LJ = 0.01 *H_modulus                    ! %%%  0.1 kpp0 =0.02 pN*um
sig_LJ = 0.1/2*(area_tot_zero/140)**0.5       ! %%%  sig_LJ um (0.2 um)
sigMM_LJ = 1.5* sum(Len_E(1:N_E))/N_E
A_LJ =4*epsl_LJ* sig_LJ**12;B_LJ =4*epsl_LJ* sig_LJ**6

AMM_LJ =0.4*epsl_LJ* sigMM_LJ**12;BMM_LJ =0.4*epsl_LJ* sigMM_LJ**6

cutoff_LJ = 2.5*sig_LJ
cut_off=2.0**(1.0/6.0)*sigMM_LJ



!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!         END of initial state construction
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

! %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!         Iteration
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



call Get_Len_E(rv,  V_E, N_E, Len_E, Vector_E)
call Get_area_F( rv, V_F, N_F,   Area_F, Vector_F, Norm_F)
call Get_vol_F( rv, V_F,  N_F,    Vol_F, Vol_total)
call Get_H_V_new(rv,N_V,N_F_V, V_V, F_V, E_V,&
&  area_F, Vector_F,  H_V,area_V,Darea_V, Dvol_V,th_E)
do i=1,N_V
Norm_V(i,1:3)=0.0
do j=1,N_V_V(i)
Norm_V(i,1:3) = Norm_V(i,1:3) + vector_F(F_V(i,j),1:3)/3.0
enddo
Norm_V(i,1:3)=Norm_V(i,1:3)  /  sum(Norm_V(i,1:3)**2)**0.5
enddo


area_tot=sum(area_F(1:N_F))


integral_H_A = 2.0*pi_K_2_A0* sum( H_V(1:N_V)*area_V(1:N_V)/3.0 )
! %%%
! % calculate energies
energyH=Curvature_Energy(H_modulus,H_v, H_zero, Area_V,N_V)
energyH=energyH+PI_K_2_A0* (sum( H_V(1:N_V)*area_V(1:N_V)/3.0 ))**2.0
energyL=kpp_area/2.0*(area_tot-area_tot_zero)**2.0/area_tot_zero
energyP=Kpp_V/2.0*(Vol_total-Vol_total_zero)**2/Vol_total_zero
call  get_energy_MS(energyMS,E_F, N_F, kpp_alpha,mu_ms,a3_ms,a4_ms,b1_ms,b2_ms,&
&area_F_zero,Len_E_zero,Len_E)
energy_tot=energyp+energyH+energyL +energyMs

! %%%



real_t=0.0;
tk=0
t=0
k_move=0

! %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

Loop_main: DO while(t.le.9000000001) 


t=t+1
V_move = mod(int( random1(seed1)*(N_V) ),N_V)+1 ! choose the move point
rv_move=rv(v_move,1:3)
k_move(V_move)=k_move(V_move)+1
! move the point
!PV_Force_POINT
lamda= Kpp_area*(area_tot-area_tot_zero)/area_tot_zero
P=    Kpp_V*(Vol_total-Vol_total_zero)/Vol_total_zero


call  POINT_PV_FORCE(Rv, V_move, V_V,F_V,N_V_V, Vector_F,area_F,Darea_Point, &
& Dvol_point,P,lamda, PV_force_point)

! MS_FORCE_POINT

call Point_MS_force(Rv,v_move,F_V,E_F,V_E,N_V_V,  kpp_alpha,mu_ms,&
& a3_ms,a4_ms,b1_ms,b2_ms,area_F_zero,Len_E_zero,Len_E,ms_force_POINT)

! PM Force

call  Force_PM_Point_analysis(RV,RV_MOVE,V_move,N_V_V,V_V,E_V,E_V_opposite,&
& F_V,len_E,th_E,vector_F,area_V,area_F,Darea_Point,H_V,H_zero,H_modulus,Integral_H_A,&
&PM_Force_Point)



! Cap
call Capillary_Force(RV_move,R_cap,A_LJ,B_LJ,cutoff_LJ,Cap_Force_point)

! Non-overlapping
if(mod(t,N_V*100).eq.1)then
do i=1,N_V
Norm_V(i,1:3)=0.0
do j=1,N_V_V(i)
Norm_V(i,1:3) = Norm_V(i,1:3) + vector_F(F_V(i,j),1:3)/3.0
enddo
Norm_V(i,1:3)=Norm_V(i,1:3)/ sum(Norm_V(i,1:3)**2)**0.5
enddo
call Get_nearV(rv,N_V,N_V_V,V_V,Norm_V,cut_off,Near_V,Num_Near_V)
endif

call Rep_Near_V_point(rv,rv_move,v_move,AMM_LJ,BMM_LJ,&
&rep_const,cut_off,Near_V,Num_Near_V, Rep_Force_V_point)

! end overlap

rv(v_move,1:3)=rv_move+dt*(PV_Force_Point+PM_Force_Point+MS_Force_Point+cap_Force_Point+&
&Rep_Force_V_point)/area_V(v_Move)*3.0



do i = 1, N_V_V(v_move)
rv(v_move,2)=rv(v_move,2)+ tao_stress/R_cap**2*(R_cap**2-rv(V_move,1)**2-rv(V_move,3)**2)*dt/area_V(v_Move)*3.0 &
& * abs(vector_F( F_V(v_move,i),2 ))/3.0
enddo


if(dt.lt.0.001)then
Gauss_KBT=(2.0*KBT/dt)**0.5
x1=random1(seed1)
y1=random1(seed1)
z1=random1(seed1)
if(x1.gt.0.0001)x1=(-2*log(x1))**0.5*cos(2*3.14159265*random1(seed1))
if(y1.gt.0.0001)y1=(-2*log(y1))**0.5*cos(2*3.14159265*random1(seed1))
if(z1.gt.0.0001)z1=(-2*log(z1))**0.5*cos(2*3.14159265*random1(seed1))


rv(v_move,1)= rv(v_move,1)+ dt*Gauss_KBT*x1
rv(v_move,2)= rv(v_move,2)+ dt*Gauss_KBT*y1
rv(v_move,3)= rv(v_move,3)+ dt*Gauss_KBT*z1
endif
rv_move=rv(v_move,1:3)

! end of MOVING

! %%%%%%%%% update new parameters around the moving point
call   Update_New_Parameters(rv,rv_move,V_move,N_V_V,V_V,E_V,F_V,F_E,E_V_opposite,V_E_opposite,&
& Len_E,th_E,vector_F,area_F,area_V,area_tot,vol_F,vol_total,H_V,Pi_K_2_A0,integral_H_A)

real_t=real_t+dt/N_V








! %%% stop for save files & show something
if(mod(t,100000).eq.1.and.dt.gt.0.0)then
tk=tk+1



call Get_Len_E(rv,  V_E, N_E, Len_E, Vector_E)
call Get_area_F( rv, V_F, N_F,   Area_F, Vector_F, Norm_F)
call Get_vol_F( rv, V_F,  N_F,    Vol_F, Vol_total)
call Get_H_V_new(rv,N_V,N_F_V, V_V, F_V, E_V,&
&  area_F, Vector_F,  H_V,area_V,Darea_V, Dvol_V,th_E)
do i=1,N_V
Norm_V(i,1:3)=0.0
do j=1,N_V_V(i)
Norm_V(i,1:3) = Norm_V(i,1:3) + vector_F(F_V(i,j),1:3)/3.0
enddo
Norm_V(i,1:3)=Norm_V(i,1:3)  /  sum(Norm_V(i,1:3)**2)**0.5
enddo
area_tot=sum(area_F(1:N_F))
integral_H_A = 2.0*pi_K_2_A0* sum( H_V(1:N_V)*area_V(1:N_V)/3.0 )
! %%%
! % calculate energies
energyH=Curvature_Energy(H_modulus,H_v, H_zero, Area_V,N_V)
energyH=energyH+PI_K_2_A0* (sum( H_V(1:N_V)*area_V(1:N_V)/3.0 ))**2.0
energyL=kpp_area/2.0*(area_tot-area_tot_zero)**2.0/area_tot_zero
energyP=Kpp_V/2.0*(Vol_total-Vol_total_zero)**2/Vol_total_zero
call  get_energy_MS(energyMS,E_F, N_F, kpp_alpha,mu_ms,a3_ms,a4_ms,b1_ms,b2_ms,&
&area_F_zero,Len_E_zero,Len_E)
energy_tot=energyp+energyH+energyL +energyMs


! %%%%%%% Cal contact area and pressure
color_F=6
num_contact=0
contact_area=0
do i=1,N_V
if(R_cap-(rv(i,1)**2+rv(i,3)**2)**0.5.lt.cutoff_LJ)then
do j=1,N_V_V(i)
color_F(F_V(i,j))=7
enddo
num_contact = num_contact+1
contact_V(num_contact) = i
contact_area=contact_area+area_V(i)/3.0

endif
enddo


contact_pressure=0
do i=1,Num_contact
v_move= contact_V(i)
lamda= Kpp_area*(area_tot-area_tot_zero)/area_tot_zero
P =    Kpp_V*(Vol_total-Vol_total_zero)/Vol_total_zero
call  POINT_PV_FORCE(Rv, V_move, V_V,F_V,N_V_V, Vector_F,area_F,Darea_Point, &
& Dvol_point,P,lamda, PV_force_point)
call Point_MS_force(Rv,v_move,F_V,E_F,V_E,N_V_V,  kpp_alpha,mu_ms,&
& a3_ms,a4_ms,b1_ms,b2_ms,area_F_zero,Len_E_zero,Len_E,ms_force_POINT)
RV_move=RV(v_move,1:3)
call  Force_PM_Point_analysis(RV,RV_MOVE,V_move,N_V_V,V_V,E_V,E_V_opposite,&
& F_V,len_E,th_E,vector_F,area_V,area_F,Darea_Point,H_V,H_zero,H_modulus,Integral_H_A,&
&PM_Force_Point)


contact_pressure = contact_pressure+ sum((PV_force_point+ms_force_POINT+PM_Force_Point)*Norm_v(v_move,1:3)) &
& /sum(Norm_v(V_move,1:3)**2)**0.5
enddo
if(Num_contact.ge.1)then

contact_pressure=contact_pressure/contact_area
! rescale p and a
contact_area = contact_area *(140/ area_tot)
contact_pressure=contact_pressure*Vol_total/100 * 0.2
endif
contact_P(tk)=contact_pressure

i=2
call Max_XYZ(rv,N_V,i,dt0)

if(R_cap.gt.d_cap_input/2.0 * (area_tot/140)**0.5) R_cap=R_cap-0.01/4.0/2.0


if(tk.gt.1000)contact_pressure = sum(contact_P(1000:tk))/(tk-1000)

if(mod(tk,1).eq.0)then
open(7,FILE='data/data'//file_input//'.dat')
write(7,*)'time and energy',t,energy_tot*0.2
write(7,*)'contact area and pressure',contact_area,contact_pressure,contact_P(tk)
write(7,*)'current and expected radius',dt0 *2 / (area_tot/140)**0.5,d_cap_input 
close(7)

! %% save the shapes
open(7,FILE='data/Rv'//file_input//'.dat')
write(7,'(A8)')'vertices'
do i=1,N_V
write(7,'(I5,3(F30.15))')i,rv(i,1),rv(i,2),rv(i,3)
enddo

write(7,*)
write(7,'(A5)')'edges'
do i=1,N_E
write(7,'(3(I7),A8,I6)')i,V_E(i,1),V_E(i,2),'color',8
enddo

write(7,*)
write(7,'(A5)')'faces'

do i=1,N_F
write(7,'(4(I7), A8, I6)')i,E_F(i,1),E_F(i,2),E_F(i,3),'color', 6 !color_F(i)
enddo

write(7,*)

write(7,*)'        '
close(7)
endif


endif


! %%%
END DO  Loop_main

! %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!       End of  Iteration
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

END
