!  ************************************************
!  Part of 3D vesicle program
!   calculate Forces
!                    written by LiJF  Dobb@bbs.fudan.edu.cn
! ************************************************* 



! Pressure Force excerted on a vertice

subroutine Pressure_Force(P,Dvol_V,N_V,   P_Force_V)
parameter(N = 20000)
integer N_V
real*8 P, Dvol_V(0:N,1:3)
real*8 P_Force_V(0:N,1:3)
P_Force_V(1:N_V,1:3)= - P* Dvol_V(1:N_V , 1:3)
END

! Curvature Force excerted on a vertice

subroutine Curvature_Force(Rv, V_V,N_V_V, H_modulus,H_v, H_zero, zet_Hfi, fi, N_V,  delt_r,     H_Force_v  )
parameter(N = 20000)

integer N_V, N_V_V(0:N)
integer V_V(0:N,1:20)

real*8 H_modulus, H_zero, zet_Hfi
real*8  fi(0:N), Area_V(0:N),H_v(0:N)
real*8 delt_r
real*8 rv(0:N,1:3)
real*8 H_Force_V(0:N,1:3)

real*8 area_V_local, area_all(1:20)
real*8 R(1:20,1:3)
real*8 R0(1:3),R1(1:3),R2(1:3),R3(1:3)

real*8 RR(1:20,1:3)
real*8 RR0(1:3)


integer i,j,k, N_t, I_temp(1:20), I_temp2(1:20)
integer i1,i2,i3

real*8 Energy_local_1,Energy_local_2, H_local,  H1(1:20),H2(1:20)

H_force_V=0.0
do i=1, N_V
R=0
I_temp=0

I_temp(1:N_V_V(i) )=V_V(i, 1:N_V_V(i) )


do j=1,N_V_V(i)
R(j,1:3)=Rv(I_temp(j),1:3 )
enddo 

! %%%%%%%%%% Force of x direction
! engery_local_1
R0(1:3)=Rv(i,1:3); 
R0(1)=R0(1)+delt_r
Energy_local_1=0.0

N_t=N_V_V(i)

call get_H_local(R0,R,N_t,   H_local, area_V_local)

Energy_local_1=Energy_local_1+H_modulus / 2.0 * (H_local-H_zero-zet_Hfi * fi(i) )**2* Area_V_local/3.0

do j=1,N_V_V(i)
RR0(1:3)=R(j,1:3)
i1=I_temp(j)  ! vertice No.
I_temp2(1:N_V_V( i1 ) )=V_V(i1, 1:N_V_V( i1 ) )
do k=1, N_V_V(i1)
RR(k,1:3)=Rv(I_temp2(k),1:3)
if(I_temp2(k).eq.i)RR(k,1:3)=R0
enddo ! end of k

N_t=N_V_V(i1)
call get_H_local(RR0,RR,N_t,   H_local, area_V_local)
Energy_local_1=Energy_local_1+H_modulus / 2.0 * (H_local-H_zero-zet_Hfi * fi(i1) )**2* Area_V_local/3.0
enddo  ! end  of j

! energy_local_2
R0(1:3)=Rv(i,1:3); 
R0(1)=R0(1)-delt_r

Energy_local_2=0.0


N_t=N_V_V(i)
call get_H_local(R0,R,N_t,   H_local, area_V_local)
Energy_local_2=Energy_local_2+H_modulus / 2.0 * (H_local-H_zero-zet_Hfi * fi(i) )**2* Area_V_local/3.0


do j=1,N_V_V(i)
RR0(1:3)=R(j,1:3)
i1=I_temp(j)  ! vertice No.
I_temp2(1:N_V_V( i1 ) )=V_V(i1, 1:N_V_V( i1 ) )
do k=1, N_V_V(i1)
RR(k,1:3)=Rv(I_temp2(k),1:3)
if(I_temp2(k).eq.i)RR(k,1:3)=R0
enddo ! end of k

N_t=N_V_V(i1)
call get_H_local(RR0,RR,N_t,   H_local, area_V_local)
Energy_local_2=Energy_local_2+H_modulus / 2.0 * (H_local-H_zero-zet_Hfi * fi(i1) )**2* Area_V_local/3.0
enddo  ! end  of j
H_Force_V(i,1)= - (Energy_Local_1-Energy_Local_2)/delt_r/2.0  

! %%%%%%%%%%%%%%

! %%%%%%%%%% Force of y direction
! engery_local_1
R0(1:3)=Rv(i,1:3); 
R0(2)=R0(2)+delt_r
Energy_local_1=0.0


N_t=N_V_V(i)
call get_H_local(R0,R,N_t,   H_local, area_V_local)
Energy_local_1=Energy_local_1+H_modulus / 2.0 * (H_local-H_zero-zet_Hfi * fi(i) )**2* Area_V_local/3.0


do j=1,N_V_V(i)
RR0(1:3)=R(j,1:3)
i1=I_temp(j)  ! vertice No.
I_temp2(1:N_V_V( i1 ) )=V_V(i1, 1:N_V_V( i1 ) )
do k=1, N_V_V(i1)
RR(k,1:3)=Rv(I_temp2(k),1:3)
if(I_temp2(k).eq.i)RR(k,1:3)=R0
enddo ! end of k

N_t=N_V_V(i1)
call get_H_local(RR0,RR,N_t,   H_local, area_V_local)
Energy_local_1=Energy_local_1+H_modulus / 2.0 * (H_local-H_zero-zet_Hfi * fi(i1) )**2* Area_V_local/3.0

enddo  ! end  of j

! energy_local_2
R0(1:3)=Rv(i,1:3); 
R0(2)=R0(2)-delt_r

Energy_local_2=0.0


N_t=N_V_V(i)
call get_H_local(R0,R,N_t,   H_local, area_V_local)
Energy_local_2=Energy_local_2+H_modulus / 2.0 * (H_local-H_zero-zet_Hfi * fi(i) )**2* Area_V_local/3.0

do j=1,N_V_V(i)
RR0(1:3)=R(j,1:3)
i1=I_temp(j)  ! vertice No.
I_temp2(1:N_V_V( i1 ) )=V_V(i1, 1:N_V_V( i1 ) )
do k=1, N_V_V(i1)
RR(k,1:3)=Rv(I_temp2(k),1:3)
if(I_temp2(k).eq.i)RR(k,1:3)=R0
enddo ! end of k

N_t=N_V_V(i1)
call get_H_local(RR0,RR,N_t,   H_local, area_V_local)
Energy_local_2=Energy_local_2+H_modulus / 2.0 * (H_local-H_zero-zet_Hfi * fi(i1) )**2* Area_V_local/3.0

enddo  ! end  of j
H_Force_V(i,2)= - (Energy_Local_1-Energy_Local_2)/delt_r/2.0  

!  %%%%%%%%%%%%%%%%%%%

! %%%%%%%%%% Force of z direction
! engery_local_1
R0(1:3)=Rv(i,1:3); 
R0(3)=R0(3)+delt_r
Energy_local_1=0.0


N_t=N_V_V(i)
call get_H_local(R0,R,N_t,   H_local, area_V_local)
Energy_local_1=Energy_local_1+H_modulus / 2.0 * (H_local-H_zero-zet_Hfi * fi(i) )**2* Area_V_local/3.0

do j=1,N_V_V(i)
RR0(1:3)=R(j,1:3)
i1=I_temp(j)  ! vertice No.
I_temp2(1:N_V_V( i1 ) )=V_V(i1, 1:N_V_V( i1 ) )
do k=1, N_V_V(i1)
RR(k,1:3)=Rv(I_temp2(k),1:3)
if(I_temp2(k).eq.i)RR(k,1:3)=R0
enddo ! end of k

N_t=N_V_V(i1)
call get_H_local(RR0,RR,N_t,   H_local, area_V_local)
Energy_local_1=Energy_local_1+H_modulus / 2.0 * (H_local-H_zero-zet_Hfi * fi(i1) )**2* Area_V_local/3.0

enddo  ! end  of j

! energy_local_2
R0(1:3)=Rv(i,1:3); 
R0(3)=R0(3)-delt_r

Energy_local_2=0.0


N_t=N_V_V(i)
call get_H_local(R0,R,N_t,   H_local, area_V_local)
Energy_local_2=Energy_local_2+H_modulus / 2.0 * (H_local-H_zero-zet_Hfi * fi(i) )**2* Area_V_local/3.0

do j=1,N_V_V(i)
RR0(1:3)=R(j,1:3)
i1=I_temp(j)  ! vertice No.
I_temp2(1:N_V_V( i1 ) )=V_V(i1, 1:N_V_V( i1 ) )
do k=1, N_V_V(i1)
RR(k,1:3)=Rv(I_temp2(k),1:3)
if(I_temp2(k).eq.i)RR(k,1:3)=R0
enddo ! end of k

N_t=N_V_V(i1)
call get_H_local(RR0,RR,N_t,   H_local, area_V_local)
Energy_local_2=Energy_local_2+H_modulus / 2.0 * (H_local-H_zero-zet_Hfi * fi(i1) )**2* Area_V_local/3.0

enddo  ! end  of j
H_Force_V(i,3)= - (Energy_Local_1-Energy_Local_2)/delt_r/2.0  



enddo ! end of i



END



!  Force of area constraint

subroutine Area_hooke_force(lamda, area_F,area_F_0,N_F, N_V, F_V, N_F_V,  Darea_V, A_Force_V )

parameter(N = 20000)

integer N_F,N_V
integer F_V(0:N,1:20), N_F_V(0:N)
real*8 lamda,area_F(0:N),area_F_0(0:N)
real*8 Darea_V(0:N,1:20,1:3),A_force_V(0:N,1:3)
real*8 A_t(1:3)

integer i,j,k,i1,i2, I_temp(1:20)


do i=1,N_V
I_temp(1:N_F_V(i))=F_V(i,1:N_F_V(i))
A_t=0.0
do j=1,N_F_V(i)
k=I_temp(j)
A_t = A_t + lamda*( area_F(k)-area_F_0(k) )/area_F(k) * Darea_V(i,j,1:3)
enddo
A_Force_V(i,1:3)=A_t
enddo


END

! Force of phi

subroutine phase_Force(Rv, F_V, V_V,N_V_V, fi,b, a2, a4,  N_V,  delt_r,     Fi_Force_v ,area_V ,area_F,A_phase)
parameter(N = 20000)

integer N_V, N_V_V(0:N)
integer V_V(0:N,1:20), F_V(0:N,1:20)

real*8 b,a2,a4, A_phase, swicth_GZ
real*8  fi(0:N), Area_V(0:N),area_F(0:N)
real*8  fit(1:3), E_inter

real*8 delt_r
real*8 rv(0:N,1:3)
real*8 Fi_Force_V(0:N,1:3)

real*8 area_V_local, area_all(1:20)
real*8 area_F_Local,Vector_F_local(1:3)
real*8 R(1:20,1:3)
real*8 R0(1:3),R1(1:3),R2(1:3),R3(1:3)

real*8 RR(1:3,1:3)



integer i,j,k, N_t, I_temp(1:20), I_temp2(1:20)
integer i1,i2,i3

real*8 Energy_local_1,Energy_local_2, H_local

 swicth_GZ=1.0  ! disable Ginzburg free energy in phase force
Fi_Force_V=0.0

do i=1, N_V
I_temp=0
R=0.0

I_temp(1:N_V_V(i) )=V_V(i, 1:N_V_V(i) )
I_temp(N_V_V(i)+1)=I_temp(1)
do j=1,N_V_V(i)
R(j,1:3)=Rv(I_temp(j),1:3 )
enddo 
R(N_V_V(i)+1,1:3)=R(1,1:3)
! %%%%%%%%%% Force of x direction
! engery_local_1
R0(1:3)=Rv(i,1:3); 
R0(1)=R0(1)+delt_r
Energy_local_1=0.0


fit(1)=fi(i)

R1(1:3)=R0
RR(1,1:3)=R1
do j=1,N_V_V(i)
i1=I_temp(j)  ! vertice No.
i2=I_temp(j+1)
R2(1:3)=R(j,1:3)
R3(1:3)=R(j+1,1:3)
RR(2,1:3)=R2;RR(3,1:3)=R3;
fit(2)=fi(i1)
fit(3)=fi(i2)
 call Get_area_local(R1,R2,R3,area_F_local, vector_F_local)
Energy_local_1=Energy_local_1+swicth_GZ*(-a2/2.0*sum(fit(1:3)**2.0) +a4/4.0*sum(fit(1:3)**4.0) )*area_F_Local/3.0

call Get_inter_phase(RR,area_F_local, fit, E_inter,b)
energy_local_1=energy_local_1+E_inter
enddo  ! end  of j







! energy_local_2
R0(1:3)=Rv(i,1:3); 
R0(1)=R0(1)-delt_r

Energy_local_2=0.0



fit(1)=fi(i)


R1(1:3)=R0
RR(1,1:3)=R1
do j=1,N_V_V(i)
i1=I_temp(j)  ! vertice No.
i2=I_temp(j+1)
R2(1:3)=R(j,1:3)
R3(1:3)=R(j+1,1:3)
RR(2,1:3)=R2;RR(3,1:3)=R3;
fit(2)=fi(i1)
fit(3)=fi(i2)
 call Get_area_local(R1,R2,R3,area_F_local, vector_F_local)
Energy_local_2=Energy_local_2+swicth_GZ*(-a2/2.0*sum(fit(1:3)**2.0) +a4/4.0*sum(fit(1:3)**4.0) )*area_F_Local/3.0

call Get_inter_phase(RR,area_F_local, fit, E_inter,b)
energy_local_2=energy_local_2+E_inter
enddo  ! end  of j


Fi_Force_V(i,1)= - (Energy_Local_1-Energy_Local_2)/delt_r/2.0  


! %%%%%%%%%%%%%%

! %%%%%%%%%% Force of y direction
! engery_local_1
R0(1:3)=Rv(i,1:3); 
R0(2)=R0(2)+delt_r
Energy_local_1=0.0



fit(1)=fi(i)

R1(1:3)=R0
RR(1,1:3)=R1
do j=1,N_V_V(i)
i1=I_temp(j)  ! vertice No.
i2=I_temp(j+1)
R2(1:3)=R(j,1:3)
R3(1:3)=R(j+1,1:3)
RR(2,1:3)=R2;RR(3,1:3)=R3;
fit(2)=fi(i1)
fit(3)=fi(i2)
 call Get_area_local(R1,R2,R3,area_F_local, vector_F_local)
Energy_local_1=Energy_local_1+swicth_GZ*(-a2/2.0*sum(fit(1:3)**2.0) +a4/4.0*sum(fit(1:3)**4.0) )*area_F_Local/3.0

call Get_inter_phase(RR,area_F_local, fit, E_inter,b)
energy_local_1=energy_local_1+E_inter
enddo  ! end  of j







! energy_local_2
R0(1:3)=Rv(i,1:3); 
R0(2)=R0(2)-delt_r

Energy_local_2=0.0

fit(1)=fi(i)

R1(1:3)=R0
RR(1,1:3)=R1
do j=1,N_V_V(i)
i1=I_temp(j)  ! vertice No.
i2=I_temp(j+1)
R2(1:3)=R(j,1:3)
R3(1:3)=R(j+1,1:3)
RR(2,1:3)=R2;RR(3,1:3)=R3;
fit(2)=fi(i1)
fit(3)=fi(i2)
 call Get_area_local(R1,R2,R3,area_F_local, vector_F_local)
Energy_local_2=Energy_local_2+swicth_GZ*(-a2/2.0*sum(fit(1:3)**2.0) +a4/4.0*sum(fit(1:3)**4.0) )*area_F_Local/3.0

call Get_inter_phase(RR,area_F_local, fit, E_inter,b)
energy_local_2=energy_local_2+E_inter
enddo  ! end  of j


Fi_Force_V(i,2)= - (Energy_Local_1-Energy_Local_2)/delt_r/2.0  

!  %%%%%%%%%%%%%%%%%%%

! %%%%%%%%%% Force of z direction
! engery_local_1
R0(1:3)=Rv(i,1:3); 
R0(3)=R0(3)+delt_r
Energy_local_1=0.0


fit(1)=fi(i)

R1(1:3)=R0
RR(1,1:3)=R1
do j=1,N_V_V(i)
i1=I_temp(j)  ! vertice No.
i2=I_temp(j+1)
R2(1:3)=R(j,1:3)
R3(1:3)=R(j+1,1:3)
RR(2,1:3)=R2;RR(3,1:3)=R3;
fit(2)=fi(i1)
fit(3)=fi(i2)
 call Get_area_local(R1,R2,R3,area_F_local, vector_F_local)
Energy_local_1=Energy_local_1+swicth_GZ*(-a2/2.0*sum(fit(1:3)**2.0) +a4/4.0*sum(fit(1:3)**4.0) )*area_F_Local/3.0

call Get_inter_phase(RR,area_F_local, fit, E_inter,b)
energy_local_1=energy_local_1+E_inter
enddo  ! end  of j







! energy_local_2
R0(1:3)=Rv(i,1:3); 
R0(3)=R0(3)-delt_r

Energy_local_2=0.0



fit(1)=fi(i)


R1(1:3)=R0
RR(1,1:3)=R1
do j=1,N_V_V(i)
i1=I_temp(j)  ! vertice No.
i2=I_temp(j+1)
R2(1:3)=R(j,1:3)
R3(1:3)=R(j+1,1:3)
RR(2,1:3)=R2;RR(3,1:3)=R3;
fit(2)=fi(i1)
fit(3)=fi(i2)
 call Get_area_local(R1,R2,R3,area_F_local, vector_F_local)
Energy_local_2=Energy_local_2+swicth_GZ*(-a2/2.0*sum(fit(1:3)**2.0) +a4/4.0*sum(fit(1:3)**4.0) )*area_F_Local/3.0

call Get_inter_phase(RR,area_F_local, fit, E_inter,b)
energy_local_2=energy_local_2+E_inter
enddo  ! end  of j


Fi_Force_V(i,3)= - (Energy_Local_1-Energy_Local_2)/delt_r/2.0  
enddo ! end of i


END


 subroutine constraint_area_Force(area_tot,area_tot_zero,area_V,Darea_V,Tension_area, P_Force_V,H_Force_V,A_Force_V , Fi_Force_v ,N_V,N_F_V,F_V,L_r,dt,  Con_area_Force_V)


parameter(N = 20000)
real*8 area_tot_zero,Tension_area, Con_area_Force_V(0:N,1:3)  ! Constraint_area Forces
real*8 area_tot, area_V(0:N)

real*8 P_force_V(0:N,1:3),H_force_V(0:N,1:3),A_force_V(0:N,1:3),Fi_Force_v(0:N,1:3)

real*8 L_r,dt,Darea_V(0:N,1:20,1:3)  ! \partial area/\partial rv == get_H_V


integer N_V,N_F_V(0:N), F_V(0:N,1:20)

integer i,j,k,N_T

real*8 Darea(0:N,1:3),Fv(0:N,1:3)

real*8 sum1,sum2

sum1=0.0
sum2=0.0
do i=1,N_V
N_T=N_F_V(i)
Darea(i,1:3)=0

do j=1,N_T
Darea(i,1:3)=Darea(i,1:3)+Darea_V(i,j,1:3)
enddo
FV(i,1:3)=(P_Force_V(i,1:3) + H_Force_V(i,1:3) + A_Force_V(i,1:3) + Fi_Force_V(i,1:3))/area_V(i)*3.0
sum1=sum1+sum(darea(i,1:3)*Fv(i,1:3))
sum2=sum2+sum(darea(i,1:3)*darea(i,1:3))/area_V(i)*3.0

enddo

Tension_area=(-0.5*(area_tot-area_tot_zero)/L_r/0.01*642.0/(N_v*1.0)-sum1 )/sum2 
! dt changed to be 0.01 
 Con_area_Force_V(1:N_v,1)=Fv(1:N_v,1)+tension_area*Darea(1:N_v,1)/area_V(1:N_V)*3.0
 Con_area_Force_V(1:N_v,2)=Fv(1:N_v,2)+tension_area*Darea(1:N_v,2)/area_V(1:N_v)*3.0
 Con_area_Force_V(1:N_v,3)=Fv(1:N_v,3)+tension_area*Darea(1:N_v,3)/area_V(1:N_v)*3.0


END

! Force of phi (analystical variation)

subroutine phase_Force_analysis(Rv, F_V, V_V,N_V_V, fi,b, a2, a4,  N_V,  delt_r,     Fi_Force_v ,area_V ,area_F,A_phase,Darea_V)
parameter(N = 20000)

integer N_V, N_V_V(0:N)
integer V_V(0:N,1:20), F_V(0:N,1:20)

real*8 b,a2,a4, A_phase, swicth_GZ
real*8  fi(0:N), Area_V(0:N),area_F(0:N)
real*8  fit(1:3), E_inter
real*8 dfi(1:3),tr(1:3,1:3)

real*8 Darea_V(0:N,1:20,1:3)
real*8 Darea(0:N,1:3)

real*8 delt_r
real*8 rv(0:N,1:3)
real*8 Fi_Force_V(0:N,1:3)
real*8 Force_temp(1:3)

real*8 area_V_local, area_all(1:20)
real*8 area_F_Local,Vector_F_local(1:3)
real*8 R(1:20,1:3)
real*8 R0(1:3),R1(1:3),R2(1:3),R3(1:3)



integer i,j,k, N_t, I_temp(1:20), I_temp2(1:20)
integer i1,i2,i3


Fi_Force_V=0.0

do i=1, N_V
I_temp=0
N_t=N_V_V(i)
I_temp(1:N_t )=V_V(i, 1:N_t )
I_temp(N_t+1)=I_temp(1)

fit(1)=fi(i)
R1(1:3)=Rv(i,1:3)

Force_temp=0.0
do j=1,N_t
i1=I_temp(j)  ! vertice No.
i2=I_temp(j+1)
R2(1:3)=Rv(i1,1:3)
R3(1:3)=Rv(i2,1:3)
fit(2)=fi(i1)
fit(3)=fi(i2)

tr(1,1:3)=R2-R1
tr(2,1:3)=R3-R2
tr(3,1:3)=R1-R3

dfi(1)=fit(2)-fit(1)
dfi(2)=fit(3)-fit(2)
dfi(3)=fit(1)-fit(3)
 call Get_area_local(R1,R2,R3,area_F_local, vector_F_local)
! GZBurg
Force_temp=Force_temp-(-a2/2.0*sum(fit(1:3)**2.0) +a4/4.0*sum(fit(1:3)**4.0) )/3.0*Darea_V(i,j,1:3)/area_F_local
! interface1
Force_temp=Force_temp+b/24.0/(area_F_local**2)*Darea_V(i,j,1:3)/area_F_local*( dfi(1)**2* (sum(tr(2,1:3)**2)+ sum(tr(3,1:3)**2))+dfi(2)**2* (sum(tr(1,1:3)**2)+ sum(tr(3,1:3)**2))+dfi(3)**2* (sum(tr(1,1:3)**2)+ sum(tr(2,1:3)**2)) - 2*dfi(1)*dfi(2)* sum(tr(1,1:3)*tr(2,1:3))   - 2*dfi(2)*dfi(3)* sum(tr(2,1:3)*tr(3,1:3) ) - 2*dfi(3)*dfi(1)* sum(tr(3,1:3)*tr(1,1:3) )   )
! integerface2
Force_temp=Force_temp-b/12.0/area_F_local* (-tr(1,1:3)*(dfi(2)**2+dfi(3)**2) +tr(3,1:3)*(dfi(1)**2+dfi(2)**2)   )
Force_temp=Force_temp-b/12.0/area_F_local* (tr(2,1:3)*dfi(2)*(dfi(1)-dfi(3)) +dfi(1)*dfi(3)*(tr(3,1:3)-tr(1,1:3))   )

enddo  ! end  of j

Fi_Force_V(i,1:3)=Force_temp
enddo ! end of i


END



! %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! %%%%%%%%%%%%%%%%%%%
subroutine Curvature_Force_analysis(Rv, V_V,N_V_V, H_modulus,H_v, H_zero, zet_Hfi, fi, N_V,  delt_r,     H_Force_v, Darea_V  )
parameter(N = 20000)

integer N_V, N_V_V(0:N)
integer V_V(0:N,1:20)

real*8 H_modulus, H_zero, zet_Hfi
real*8  fi(0:N), Area_V(0:N),H_v(0:N)
real*8 delt_r
real*8 rv(0:N,1:3)
real*8 H_Force_V(0:N,1:3)

real*8 Force_t(1:3)

real*8 area_V_local, area_all(0:20)
real*8 area_local,vector_local(1:3),Vector_R0(0:20,1:3),Vector_Rk(0:20,1:3)
real*8 area_R0(0:20),area_Rk(0:20)
real*8 ForceV_V(1:3),Norm_V_temp(1:3)
real*8 R1R2_fun,R1R2_R3_fun, R1_R2R3_fun

real*8 R(1:20,1:3)
real*8 R0(1:3),R1(1:3),R2(1:3),R3(1:3)

real*8 RR(0:20,1:3)
real*8 RR0(1:3)

real*8 Pf(0:20,1:3)
real*8 Dfi(0:20,1:3)
real*8 Df0(1:3)
real*8 Dni(0:20,1:3)
real*8 Dn0(1:3)
real*8 Darea_V(0:N,1:20,1:3)
real*8 Fh0,Fhi(0:20)
real*8 dFh0,dFhi(0:20)
real*8 Dh0(1:3),Dhi(0:20,1:3)
real*8 Sv0,Svi(0:20)

real*8 Det_h0(1:3),Det_Hi(0:20,1:3)
real*8 Det_Si(0:20,1:3)

real*8 vect_t1(1:3),vect_t2(1:3),vect_t3(1:3),vect_t4(1:3)  ! BAK VECTORS

integer i,j,k, N_t, I_temp(1:20), I_temp2(1:20)
integer i1,i2,i3,k1,k2,k3, m

real*8 Energy_local_1,Energy_local_2, H_local


H_force_V=0.0
do i=1, N_V
Force_t=0.0
R=0.0
I_temp=0.0
I_temp(1:N_V_V(i) )=V_V(i, 1:N_V_V(i) )
do j=1,N_V_V(i)
R(j,1:3)=Rv(I_temp(j),1:3 )
enddo 
R(N_V_V(i)+1,1:3)=R(1,1:3)
R0(1:3)=Rv(i,1:3); 

N_t=N_V_V(i)
! %%%%%%%%Dh0
area_V_local=0.0
do j=1,N_t
R1(1:3)=R0(1:3);R2(1:3)=R(j,1:3);R3(1:3)=R(j+1,1:3);
call Get_area_local(R1,R2,R3,area_local, vector_local)
area_R0(j)=area_local
vector_R0(j,1:3)= vector_local(1:3)
area_V_local=area_V_local+area_local
Det_Si(j,1:3)=Darea_V(i,j,1:3)
enddo
Det_Si(0,1:3)=Det_Si(N_t,1:3)
Det_Si(N_t+1,1:3)=Det_Si(1,1:3)
area_R0(0)=area_R0(N_t)
area_R0(N_t+1)=area_R0(1)
Sv0=area_V_local/3.0

Norm_V_temp=0.0
ForceV_V=0
do j=1,N_t
Norm_V_temp(1)=Norm_V_temp(1)+R1R2_fun(R(j,1:3), R(j+1,1:3), 1)
Norm_V_temp(2)=Norm_V_temp(2)+R1R2_fun(R(j,1:3), R(j+1,1:3), 2)
Norm_V_temp(3)=Norm_V_temp(3)+R1R2_fun(R(j,1:3), R(j+1,1:3), 3)
Pf(j,1)=-( R1R2_fun(R(j,1:3), vector_R0(j,1:3), 1) +R1R2_fun( vector_R0(j,1:3),R(j+1,1:3), 1)  )/6.0
Pf(j,2)=-( R1R2_fun(R(j,1:3), vector_R0(j,1:3), 2) +R1R2_fun( vector_R0(j,1:3),R(j+1,1:3), 2)  )/6.0
Pf(j,3)=-( R1R2_fun(R(j,1:3), vector_R0(j,1:3), 3) +R1R2_fun( vector_R0(j,1:3),R(j+1,1:3), 3)  )/6.0
ForceV_V=ForceV_V+Pf(j,1:3)/area_R0(j)
enddo
Norm_v_temp=Norm_V_temp/6.0
det_1=(ForceV_V(1)*Norm_V_temp(1)+ForceV_V(2)*Norm_V_temp(2)+ForceV_V(3)*Norm_V_temp(3))
if(abs(det_1).lt.0.000000001)then
H_local=0.0
Df0=0.0
Dn0=0.0
else
H_local=2*1.5*(ForceV_V(1)**2+ForceV_V(2)**2+ForceV_V(3)**2)/(ForceV_V(1)*Norm_V_temp(1)+ForceV_V(2)*Norm_V_temp(2)+ForceV_V(3)*Norm_V_temp(3))

Df0=6.0*ForceV_V/det_1 - 3.0*sum(ForceV_V**2.0)*Norm_V_temp/det_1**2.0
Dn0=- 3.0*sum(ForceV_V**2.0)*ForceV_V/det_1**2.0
endif

Det_h0=0.0
Vect_t1=0.0
do j=1,N_t
R1=R(j,1:3)
R2=R(j+1,1:3)
!Pf(j,1:3)=-1.0/3.0*Det_Si(j,1:3)
do k=1,3
Det_h0(k)=Det_h0(k)-1.0/12.0/area_R0(j)*(R1_R2R3_fun(R1,Df0,R1,k)+R1_R2R3_fun(R1,R2,Df0,k) )
Det_h0(k)=Det_h0(k)-1.0/12.0/area_R0(j)*(R1R2_R3_fun(Df0,R1,R2,k)+R1R2_R3_fun(R2,Df0,R2,k) )
Det_h0(k)=Det_h0(k)-Det_Si(j,k)/area_R0(j)**3*sum(Df0*Pf(j,1:3))
enddo
Vect_t1=Vect_t1+Det_si(j,1:3)/area_R0(j)
enddo
Fh0=H_modulus / 2.0 * (H_local-H_zero-zet_Hfi * fi(i) )**2
dFh0=H_modulus * (H_local-H_zero-zet_Hfi * fi(i) )
Force_t=Force_t- dFh0*Sv0*Det_h0 +Fh0*ForceV_V

!write(*,*)Det_h0(1:3)
!pause

! %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! Dhi

Det_hi=0.0
do j=1,N_V_V(i)
RR0(1:3)=R(j,1:3)
i1=I_temp(j)  ! vertice No.
I_temp2(1:N_V_V( i1 ) )=V_V(i1, 1:N_V_V( i1 ) )
do k=1, N_V_V(i1)
RR(k,1:3)=Rv(I_temp2(k),1:3)
If(I_temp2(k).eq.i)then
k1=k-1
k2=k
k3=k+1
if(k.eq.1)k1=N_V_V(i1)
if(k.eq.N_V_V(i1))k3=1
endif
enddo ! end of k

RR(N_V_V(i1)+1,1:3)=RR(1,1:3)

area_V_local=0.0
do k=1,N_V_V(i1)
R1(1:3)=RR0(1:3);R2(1:3)=RR(k,1:3);R3(1:3)=RR(k+1,1:3);
call Get_area_local(R1,R2,R3,area_local, vector_local)
area_Rk(k)=area_local
vector_Rk(k,1:3)= vector_local(1:3)
area_V_local=area_V_local+area_local
enddo
area_Rk(0)=area_Rk(N_V_V(i1))
area_Rk(N_V_V(i1)+1)=area_Rk(1)
vector_Rk(0,1:3)=vector_Rk(N_V_V(i1),1:3)
vector_Rk(N_V_V(i1)+1,1:3)=vector_Rk(1,1:3)

Svi(j)=area_V_local/3.0

Norm_V_temp=0.0
ForceV_V=0.0
do k=1,N_V_V(i1)
Norm_V_temp(1)=Norm_V_temp(1)+R1R2_fun(RR(k,1:3), RR(k+1,1:3), 1)
Norm_V_temp(2)=Norm_V_temp(2)+R1R2_fun(RR(k,1:3), RR(k+1,1:3), 2)
Norm_V_temp(3)=Norm_V_temp(3)+R1R2_fun(RR(k,1:3), RR(k+1,1:3), 3)
Pf(k,1)=-( R1R2_fun(RR(k,1:3), vector_Rk(k,1:3), 1) +R1R2_fun( vector_Rk(k,1:3),RR(k+1,1:3), 1)  )/6.0
Pf(k,2)=-( R1R2_fun(RR(k,1:3), vector_Rk(k,1:3), 2) +R1R2_fun( vector_Rk(k,1:3),RR(k+1,1:3), 2)  )/6.0
Pf(k,3)=-( R1R2_fun(RR(k,1:3), vector_Rk(k,1:3), 3) +R1R2_fun( vector_Rk(k,1:3),RR(k+1,1:3), 3)  )/6.0
ForceV_V=ForceV_V+Pf(k,1:3)/area_Rk(k)
enddo
Pf(0,1:3)=Pf(N_V_V(i1),1:3)
Pf(N_V_V(i1)+1,1:3)=Pf(1,1:3)

Norm_v_temp=Norm_V_temp/6.0
det_1=(ForceV_V(1)*Norm_V_temp(1)+ForceV_V(2)*Norm_V_temp(2)+ForceV_V(3)*Norm_V_temp(3))
if(abs(det_1).lt.0.00000000001)then
H_local=0.0
Dfi(j,1:3)=0.0
Dni(j,1:3)=0.0
else
H_local=2*1.5*(ForceV_V(1)**2+ForceV_V(2)**2+ForceV_V(3)**2)/(ForceV_V(1)*Norm_V_temp(1)+ForceV_V(2)*Norm_V_temp(2)+ForceV_V(3)*Norm_V_temp(3))

Dfi(j,1:3)=6*ForceV_V/det_1 - 3.0*sum(ForceV_V**2.0)*Norm_V_temp/det_1**2.0
Dni(j,1:3)=-3.0*sum(ForceV_V**2.0)*ForceV_V/det_1**2.0
endif


R1=RR(k1,1:3)
R2=RR(k2,1:3)
R3=RR(k3,1:3)


vect_t1=vector_Rk(k1,1:3)*2.0
vect_t2=Vector_Rk(k2,1:3)*2.0

do m=1,3
Det_hi(j,m)=Det_hi(j,m)-1.0/12.0/area_RK(k2)*R1R2_fun(vect_t2,Dfi(j,1:3),m)
Det_hi(j,m)=Det_hi(j,m)-1.0/12.0/area_Rk(k2)*(R1_R2R3_fun(RR0,R2,Dfi(j,1:3),m)+R1_R2R3_fun(RR0,Dfi(j,1:3),R3,m) +R1R2_R3_fun(R2,Dfi(j,1:3),R3,m)+R1R2_R3_fun(Dfi(j,1:3),R3,R3,m) )
Det_hi(j,m)=Det_hi(j,m)-1.0/12.0/area_Rk(k1)*R1R2_fun(Dfi(j,1:3),vect_t1,m)
Det_hi(j,m)=Det_hi(j,m)-1.0/12.0/area_Rk(k1)*(R1_R2R3_fun(R1,R1,Dfi(j,1:3),m)+R1_R2R3_fun(R1,Dfi(j,1:3),R2,m) +R1R2_R3_fun(R1,Dfi(j,1:3),RR0,m)+R1R2_R3_fun(Dfi(j,1:3),R2,RR0,m) )
Det_hi(j,m)=Det_hi(j,m)-Det_Si(j-1,m)/area_R0(j-1)**3*sum(Dfi(j,1:3)*Pf(k2,1:3))-Det_Si(j,m)/area_R0(j)**3*sum(Dfi(j,1:3)*Pf(k1,1:3))
Det_hi(j,m)=Det_hi(j,m)+1.0/6.0*(R1R2_fun(R3,Dni(j,1:3),m)+R1R2_fun(Dni(j,1:3),R1,m))
enddo
Fhi(j)=H_modulus / 2.0 * (H_local-H_zero-zet_Hfi * fi(i1) )**2
dFhi(j)=H_modulus * (H_local-H_zero-zet_Hfi * fi(i1) )
Force_t=Force_t-dFhi(j)*Svi(j)*Det_hi(j,1:3)-Fhi(j)/3.0*(Det_Si(j-1,1:3)/area_R0(j-1)+Det_Si(j,1:3)/area_R0(j))
enddo  ! end  of j

H_Force_V(i,1:3)=Force_t
!write(*,*)det_Hi(2,1)
!pause

enddo ! end of i


END








