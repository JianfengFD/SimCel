!  ************************************************
!  Part of 3D vesicle program
!    update the vertices' positions
!                    written by LiJF  Dobb@bbs.fudan.edu.cn
! ************************************************* 


subroutine update_R(rv0,rv, P_Force_V,H_Force_V,A_Force_V , Fi_Force_v , N_V, L_r,dt,area_V)
parameter(N = 20000)
integer N_V
real*8 rv0(0:N,1:3),P_Force_V(0:N,1:3),H_Force_V(0:N,1:3),A_Force_V(0:N,1:3)
real*8 rv(0:N,1:3)  ,Fi_Force_V(0:N,1:3),area_V(0:N)
real*8 L_r,dt

integer i


do i=1,N_V
rv(i,1)=rv0(i,1)+L_r*dt*(P_Force_V(i,1) + H_Force_V(i,1) + A_Force_V(i,1) + Fi_Force_V(i,1))/area_V(i)*3.0
rv(i,2)=rv0(i,2)+L_r*dt*(P_Force_V(i,2) + H_Force_V(i,2) + A_Force_V(i,2) + Fi_Force_V(i,2))/area_V(i)*3.0
rv(i,3)=rv0(i,3)+L_r*dt*(P_Force_V(i,3) + H_Force_V(i,3) + A_Force_V(i,3) + Fi_Force_V(i,3))/area_V(i)*3.0
enddo





END


subroutine update_R_optimum_scale(rv0, rv, P_Force_V,H_Force_V,A_Force_V , Fi_force_V,N_V, L_r,S_opt, dt,   P,Vol_total,   H_modulus,H_v, H_zero, zet_Hfi, fi, Area_V, lamda, area_F,Area_F_zero,N_F,      N_E,N_F_V, V_F, V_E, E_F,  V_V,N_V_V, F_V, E_V,E_V_opposite,    Vector_F, norm_F,Darea_V, Dvol_V ,   Vol_F,delt_r,energy_tot)
parameter(N = 20000)
integer N_V,N_E,N_F
real*8 rv0(0:N,1:3),P_Force_V(0:N,1:3),H_Force_V(0:N,1:3),A_Force_V(0:N,1:3),Fi_Force_V(0:N,1:3)
real*8 rv(0:N,1:3)
real*8 L_r,S_opt, s1,s2,s3,st,s0,dt

real*8 PV_energy,Curvature_energy,Area_Hooke_energy
real*8 EnergyP,EnergyH,energyL, energy_tot
real*8 H_modulus, lamda, H_zero, zet_Hfi, P

real*8 e0,e1,e2,e3,et,e_array(1:200)

integer V_E(0:N, 1:2)
integer E_F(0:N, 1:3)
! the first kind label
integer V_F(0:N,1:3)  ! vertices around a face
integer E_V(0:N,1:20), N_E_V(0:N) ! edges ajecent to a vertice
integer F_V(0:N,1:20), N_F_V(0:N) ! star faces around a vertice
integer F_E(0:N,1:20), N_F_E(0:N)
! second kind label
integer V_V(0:N,1:20), N_V_V(0:N)
integer E_V_opposite(0:N,1:20), N_E_V_opposite(0:N)
integer V_E_opposite(0:N,1:20)
! %%%%%%%%%%%%%%%%%%%%%%%%%%
! %%%%%%%%%%%%%
!Basic parameters
real*8 Vector_F(0:N,1:3), Area_F(0:N), Norm_F(0:N,1:3) 
real*8 Area_F_zero(0:N),delt_r

! == get_area_F
real*8 Vol_F(0:N), vol_total ! == get_vol_F 
real*8 H_v(0:N) ,Area_V(0:N)   ! == get_H_V
real*8 Darea_V(0:N,1:20,1:3)  ! \partial area/\partial rv == get_H_V
real*8 DVol_V(0:N,1:3)  ! \partial volume/\partial rv == get_H_V
real*8 fi(0:N)

integer i,j,k,i1,I_plus, I_E


e0=energy_tot
s0 = dt
I_E=40
i1=1
I_plus=1
k=0
st=s0
do while(i1.eq.1)

do i=1,N_V
rv(i,1)=rv0(i,1)+L_r*st*(P_Force_V(i,1) + H_Force_V(i,1) + A_Force_V(i,1))
rv(i,2)=rv0(i,2)+L_r*st*(P_Force_V(i,2) + H_Force_V(i,2) + A_Force_V(i,2))
rv(i,3)=rv0(i,3)+L_r*st*(P_Force_V(i,3) + H_Force_V(i,3) + A_Force_V(i,3))
enddo
! %%% calculate basic parameters
call Get_area_F( rv, V_F,  N_F,N_V,   Area_F, Vector_F, Norm_F)
call Get_vol_F( rv, V_F,  N_F,N_V,    Vol_F, Vol_total)
call Get_H_V(rv, N_F,N_V,N_E,N_F_V, V_F, V_E, E_F,  V_V, F_V, E_V,E_V_opposite,  area_F, Vector_F, norm_F, H_V,area_V,Darea_V, Dvol_V)
! %%%

! % calculate energies
ENERGYP=PV_energy(P,Vol_total)
energyH=Curvature_Energy(H_modulus,H_v, H_zero, zet_Hfi, fi, Area_V,N_V)
energyL=Area_Hooke_energy(lamda, area_F,Area_F_zero,N_F) !
et=energyp+energyH+energyL

if(abs(k).gt.20)then
write(*,*)'Optimal error'
i1=0
endif

if(k.eq.0)then
st=dt*2.0**I_plus
s0=st
e0=et
e_array(k+I_E)=e0

k=1
else 

if(et.gt.e0.and.k.eq.1)then
I_plus=-1
st=dt*2.0**I_plus
s0=st
e0=et
e_array(k+I_E)=e0

k=-1
endif
if(et.lt.e0.and.k.ge.1)then
I_plus=1
st=s0*2.0**I_plus
s0=st
e0=et
e_array(k+I_E)=e0
k=k+1
endif
if(et.gt.e0.and.k.gt.1)then
I_plus=1
s0=st
e0=et
e_array(k+I_E)=e0
i1=0  ! end do
endif

if(et.lt.e0.and.k.lt.0)then
I_plus=-1
st=s0*2.0**I_plus
s0=st
e0=et
e_array(k+I_E)=e0

k=k-1
endif
if(et.gt.e0.and.k.lt.0)then
I_plus=-1
s0=st
e0=et
e_array(k+I_E)=e0
i1=0

endif

endif

enddo

if(k.gt.0)then
e1=e_array(k+I_E-2)
e2=e_array(k+I_E-1)
e3=e_array(k+I_E)
s3=s0
s2=s3/2.0
s1=s2/2.0
endif

if(k.le.0)then
e1=e_array(k+I_E)
e2=e_array(k+I_E+1)
e3=e_array(k+I_E+2)
s1=s0
s2=s1*2.0
s3=s2*2.0
endif

S_opt=0.75*s2*(4*e1-5*e2+e3)/(2*e1-3*e2+e3)
dt=s_opt

e0=(e1*(s2-dt)**2-e2*(s1-dt)**2)/((s2-dt)**2-(s1-dt)**2)


if(abs(k).gt.20)dt=0.00001
if(e2.gt.energy_tot+1.0)dt=0.00001
rv=rv0
! %%% calculate basic parameters
call Get_area_F( rv, V_F,  N_F,N_V,   Area_F, Vector_F, Norm_F)
call Get_vol_F( rv, V_F,  N_F,N_V,    Vol_F, Vol_total)
call Get_H_V(rv, N_F,N_V,N_E,N_F_V, V_F, V_E, E_F,  V_V, F_V, E_V,E_V_opposite,  area_F, Vector_F, norm_F, H_V,area_V,Darea_V, Dvol_V)
! %%%


! %%% calculate Forces
call Pressure_Force(P,Dvol_V,N_V,   P_Force_V)
call Curvature_Force(Rv, V_V,N_V_V, H_modulus,H_v, H_zero, zet_Hfi, fi, N_V,  delt_r,     H_Force_v  )
call Area_hooke_force(lamda, area_F,Area_F_zero,N_F, N_V, F_V, N_F_V,  Darea_V, A_Force_V )
! %%%
!Fi_Force_V=0.0
!call update_R(rv0,rv, P_Force_V,H_Force_V,A_Force_V ,Fi_Force_V, N_V, L_r,dt)

do i=1,N_V
rv(i,1)=rv0(i,1)+L_r*dt*(P_Force_V(i,1) + H_Force_V(i,1) + A_Force_V(i,1))
rv(i,2)=rv0(i,2)+L_r*dt*(P_Force_V(i,2) + H_Force_V(i,2) + A_Force_V(i,2))
rv(i,3)=rv0(i,3)+L_r*dt*(P_Force_V(i,3) + H_Force_V(i,3) + A_Force_V(i,3))
enddo


END



subroutine Jiggle(rv,N_v,dt,amp_jig,seed1)
parameter(N = 20000)
integer N_V
real*8 rv(0:N,1:3),dt,amp_jig
integer seed1
real*8 random1

integer i,j,k

do i=1,N_v
rv(i,1)=rv(i,1)+ (random1(seed1)+random1(seed1)+random1(seed1)+ random1(seed1)+random1(seed1)+random1(seed1) -3.0)/3.0*amp_jig
rv(i,2)=rv(i,2)+ (random1(seed1)+random1(seed1)+random1(seed1)+ random1(seed1)+random1(seed1)+random1(seed1) -3.0)/3.0*amp_jig
rv(i,3)=rv(i,3)+ (random1(seed1)+random1(seed1)+random1(seed1)+ random1(seed1)+random1(seed1)+random1(seed1) -3.0)/3.0*amp_jig
enddo



End



