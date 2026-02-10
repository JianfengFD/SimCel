


subroutine POINT_PV_FORCE(Rv, V_move, V_V,F_V,N_V_V, Vector_F,area_F,Darea_Point, &
& Dvol_point,P,lamda, PV_force_point)
implicit none
integer N
parameter(N = 50000)

real*8 rv(0:N,1:3),area_F(0:N)

integer V_V(0:N,1:20),N_V_V(0:N),i,j
integer V_MOVE,F_V(0:N,1:20)

real*8 P,lamda,PV_Force_Point(1:3)

real*8 R(1:20,1:3)
real*8  Vect_VV(1:20,1:3),Dvol_Point(1:3)

real*8 vector_F(0:N,1:3),at(1:3),Darea_Point(1:20,1:3)

integer N_V_V_local

real*8 R1R2_fun


Dvol_Point=0.0
PV_Force_Point=0.0
! get the variation

N_V_V_local=N_V_V(v_move)

R(1:N_V_V_local,1:3)=Rv(V_V(v_move,1:N_V_V_local),1:3)
R(N_V_V_local+1,1:3)=R(1,1:3)

do j=1,N_V_V_local
vect_VV(j,1:3)= vector_F(F_V(V_move,j),1:3)

enddo


do j=1,N_V_V_local
Dvol_Point(1)=Dvol_Point(1)+R1R2_fun(R(j,1:3), R(j+1,1:3), 1)
Dvol_Point(2)=Dvol_Point(2)+R1R2_fun(R(j,1:3), R(j+1,1:3), 2)
Dvol_Point(3)=Dvol_Point(3)+R1R2_fun(R(j,1:3), R(j+1,1:3), 3)
at(1)=( R1R2_fun(R(j,1:3), Vect_VV(j,1:3), 1) +R1R2_fun( Vect_VV(j,1:3),R(j+1,1:3), 1)  )/ 2.0
at(2)=( R1R2_fun(R(j,1:3), Vect_VV(j,1:3), 2) +R1R2_fun( Vect_VV(j,1:3),R(j+1,1:3), 2)  )/ 2.0
at(3)=( R1R2_fun(R(j,1:3), Vect_VV(j,1:3), 3) +R1R2_fun( Vect_VV(j,1:3),R(j+1,1:3), 3)  )/ 2.0
Darea_point(j,1:3)=at
enddo

Dvol_Point=Dvol_Point/6.0

! %%%%%%%% pressure part

PV_Force_Point = PV_Force_Point - Dvol_point*P



! %%%%%%% area part
do i=1,N_V_V_local
PV_Force_Point = PV_Force_Point - lamda * Darea_Point(i,1:3)/area_F(F_V(v_move,i))
enddo



END



subroutine POINT_PV_FORCE_Numerical(Rv, V_move, V_V,F_V,N_V_V, &
& P,lamda, PV_force_point)
implicit none
integer N
parameter(N = 50000)

real*8 rv(0:N,1:3)

integer V_V(0:N,1:20),N_V_V(0:N),i,j
integer V_MOVE,F_V(0:N,1:20)

real*8 P,lamda,PV_Force_Point(1:3)

real*8 R(1:20,1:3),RR(1:3,1:3)
real*8 R1(1:3),R2(1:3),R3(1:3)
real*8 det_r
real*8 energyS1,energyS2,energyV1,energyV2
real*8 s1,s2,v1,v2
real*8 R0_old(1:3),R0_New(1:3)

real*8  Vect_VV(1:20,1:3),Dvol_Point(1:3),vect1(1:3)

real*8 vector_F(0:N,1:3),at(1:3),DS_Point(1:3)

integer N_t

real*8 R1R2_fun,V6_FUN


Dvol_Point=0.0
PV_Force_Point=0.0
det_r=0.000001

N_t=N_V_V(v_move)

R(1:N_t,1:3)=Rv(V_V(v_move,1:N_t),1:3)
R(N_t+1,1:3)=R(1,1:3)
R0_old(1:3)=Rv(v_move,1:3)
! %%%%% x direction
!    x energy 1
    R0_New(1:3)= R0_old(1:3);  R0_New(1)= R0_old(1)+det_r; 
    S2=0;V2=0
	do i=1,N_t
		RR(1,1:3)=R0_NEW(1:3); RR(2,1:3)= R(i,1:3);RR(3,1:3)= R(i+1,1:3)
		V2=V2+V6_FUN(RR)/6.0
		R1(1:3)=R0_NEW(1:3); R2(1:3)= R(i,1:3);R3(1:3)= R(i+1,1:3)
		call Get_area_local(R1,R2,R3,S1, vect1)
		S2=S2+S1
	enddo 
	
	energyS1= S2 * lamda
	energyV1= V2 * P
!    x energy 2
    R0_New(1:3)= R0_old(1:3);  R0_New(1)= R0_old(1)-det_r; 
    S2=0;V2=0
	do i=1,N_t
		RR(1,1:3)=R0_NEW(1:3); RR(2,1:3)= R(i,1:3);RR(3,1:3)= R(i+1,1:3)
		V2=V2+V6_FUN(RR)/6.0
		R1(1:3)=R0_NEW(1:3); R2(1:3)= R(i,1:3);R3(1:3)= R(i+1,1:3)
		call Get_area_local(R1,R2,R3,S1, vect1)
		S2=S2+S1
	enddo 
	
	energyS2= S2 * lamda
	energyV2= V2 * P

   Dvol_Point(1)=(energyV1-energyV2)/2/det_r
   DS_Point(1)=(energyS1-energyS2)/2/det_r

! %%%%% y direction
!    y energy 1
    R0_New(1:3)= R0_old(1:3);  R0_New(2)= R0_old(2)+det_r; 
    S2=0;V2=0
	do i=1,N_t
		RR(1,1:3)=R0_NEW(1:3); RR(2,1:3)= R(i,1:3);RR(3,1:3)= R(i+1,1:3)
		V2=V2+V6_FUN(RR)/6.0
		R1(1:3)=R0_NEW(1:3); R2(1:3)= R(i,1:3);R3(1:3)= R(i+1,1:3)
		call Get_area_local(R1,R2,R3,S1, vect1)
		S2=S2+S1
	enddo 
	
	energyS1= S2 * lamda
	energyV1= V2 * P
!    y energy 2
    R0_New(1:3)= R0_old(1:3);  R0_New(2)= R0_old(2)-det_r; 
    S2=0;V2=0
	do i=1,N_t
		RR(1,1:3)=R0_NEW(1:3); RR(2,1:3)= R(i,1:3);RR(3,1:3)= R(i+1,1:3)
		V2=V2+V6_FUN(RR)/6.0
		R1(1:3)=R0_NEW(1:3); R2(1:3)= R(i,1:3);R3(1:3)= R(i+1,1:3)
		call Get_area_local(R1,R2,R3,S1, vect1)
		S2=S2+S1
	enddo 
	
	energyS2= S2 * lamda
	energyV2= V2 * P

   Dvol_Point(2)=(energyV1-energyV2)/2/det_r
   DS_Point(2)=(energyS1-energyS2)/2/det_r

! %%%%% z direction
!    z energy 1
    R0_New(1:3)= R0_old(1:3);  R0_New(3)= R0_old(3)+det_r; 
    S2=0;V2=0
	do i=1,N_t
		RR(1,1:3)=R0_NEW(1:3); RR(2,1:3)= R(i,1:3);RR(3,1:3)= R(i+1,1:3)
		V2=V2+V6_FUN(RR)/6.0
		R1(1:3)=R0_NEW(1:3); R2(1:3)= R(i,1:3);R3(1:3)= R(i+1,1:3)
		call Get_area_local(R1,R2,R3,S1, vect1)
		S2=S2+S1
	enddo 
	
	energyS1= S2 * lamda
	energyV1= V2 * P
!    z energy 2
    R0_New(1:3)= R0_old(1:3);  R0_New(3)= R0_old(3)-det_r; 
    S2=0;V2=0
	do i=1,N_t
		RR(1,1:3)=R0_NEW(1:3); RR(2,1:3)= R(i,1:3);RR(3,1:3)= R(i+1,1:3)
		V2=V2+V6_FUN(RR)/6.0
		R1(1:3)=R0_NEW(1:3); R2(1:3)= R(i,1:3);R3(1:3)= R(i+1,1:3)
		call Get_area_local(R1,R2,R3,S1, vect1)
		S2=S2+S1
	enddo 
	
	energyS2= S2 * lamda
	energyV2= V2 * P

   Dvol_Point(3)=(energyV1-energyV2)/2/det_r
   DS_Point(3)=(energyS1-energyS2)/2/det_r



PV_Force_Point= -(Dvol_Point+DS_Point)

END















! %%%%%%%%%%%%%%%%%%  MS part


subroutine Point_MS_force(Rv,v_move,F_V,E_F,V_E,N_V_V,  kpp_alpha,mu_ms,&
& a3_ms,a4_ms,b1_ms,b2_ms,area_F_zero,Len_E_zero,Len_E,ms_force_POINT)
implicit none
integer N
parameter(N = 50000)

integer N_V_V(0:N)
integer F_V(0:N,1:20),V_E(0:N,1:2),E_F(0:N,1:3)

real*8 kpp_alpha,mu_ms,a3_ms,a4_ms,b1_ms,b2_ms
real*8 rv(0:N,1:3)

real*8 area_F_zero(0:N)
real*8 Len_E_zero(0:N),Len_E(0:N)

! %%%%%% MS
real*8 a,b,c,a0,b0,c0,alph,bet

real*8 MS_Force_point(1:3)

real*8 R0(1:3),R1(1:3),R2(1:3)


integer i,j,k
integer i3,k_c,k_a,k_b,v_move

real*8  DFbet,S,S0,DFalph
real*8 DS(1:3),DPI(1:3)

integer  num_F,I_side(1:3)

MS_Force_point=0.0



r0=rv(v_move,1:3)
i=v_move

do j=1,N_V_V(i)
num_F=F_V(i,j)
I_side(1:3)=abs(E_F(num_F,1:3))
do k=1,3
if(V_E(I_side(k),1).ne.i.and.V_E(I_side(k),2).ne.i)k_c=k
enddo

i3=0
do k=1,3
if(k.ne.k_c.and.i3.eq.1)then
k_b=k
if(V_E(I_side(k),1).ne.i)r2=rv(V_E(I_side(k),1),1:3)
if(V_E(I_side(k),2).ne.i)r2=rv(V_E(I_side(k),2),1:3)
i3=i3+1
endif
if(k.ne.k_c.and.i3.eq.0)then
k_a=k
if(V_E(I_side(k),1).ne.i)r1=rv(V_E(I_side(k),1),1:3)
if(V_E(I_side(k),2).ne.i)r1=rv(V_E(I_side(k),2),1:3)
i3=i3+1
endif
enddo

a0=len_E_zero(I_side(k_a))
b0=len_E_zero(I_side(k_b))
 c0=len_E_zero(I_side(k_c))
a =len_E(I_side(k_a))
b =len_E(I_side(k_b))
 c =len_E(I_side(k_c))


call Get_alph_bet_analytic(a,b,c,a0,b0,c0,alph,bet,S,S0)

DFalph=(Kpp_alpha/2.0*(2.0*alph+3.0*a3_ms*alph**2+4.0*a4_ms*alph**3)+mu_ms*b2_ms*bet)* &
& area_F_zero(num_F)
DFbet= mu_ms*(1.0+b1_ms*alph+2*b2_ms*bet)*area_F_zero(NUM_F)

DS=-4*a**2*(r0-r1)-4*b**2*(r0-r2)+4*a**2*(r0-r2)+4*b**2*(r0-r1)+4*c**2*(r0-r1)+4*c**2*(r0-r2)

DPI=-2*a0**2*(r0-r1)-2*b0**2*(r0-r2)+2*(b0**2+c0**2)*(r0-r1)+2*(a0**2+c0**2)*(r0-r2)

MS_Force_Point=MS_Force_Point-DFalph/2.0/(S*S0)**0.5*Ds
MS_Force_Point=MS_Force_Point-DFbet*(DPI/(S*S0)**0.5-(bet+1)/2.0/S*Ds)
enddo ! end of N_V_V




END

subroutine Get_alph_bet(a,b,c,a0,b0,c0,alph,bet,area_F_temp)
real*8 S,S0,SS0,area_F_temp
real*8 a,b,c,a0,b0,c0,alph,bet

S  =(2*a**2*b**2+2*a**2*c**2+2*b**2*c**2-a**4-b**4-c**4)
S0 =(2*a0**2*b0**2+2*a0**2*c0**2+2*b0**2*c0**2-a0**4-b0**4-c0**4)
SS0=a**2*b0**2+a**2*c0**2+b**2*c0**2+a0**2*b**2+a0**2*c**2+b0**2*c**2&
&-a**2*a0**2-b0**2*b**2-c**2*c0**2


alph=(S/S0)**0.5-1.0

bet=SS0/(S*S0)**0.5-1.0
area_F_temp=(S/4.0)**0.5/2.0

END


subroutine Get_alph_bet_analytic(a,b,c,a0,b0,c0,alph,bet,S,S0)
real*8 S,S0,SS0
real*8 a,b,c,a0,b0,c0,alph,bet


S  =(2*a**2*b**2+2*a**2*c**2+2*b**2*c**2-a**4-b**4-c**4)
S0 =(2*a0**2*b0**2+2*a0**2*c0**2+2*b0**2*c0**2-a0**4-b0**4-c0**4)
SS0=a**2*b0**2+a**2*c0**2+b**2*c0**2+a0**2*b**2+a0**2*c**2+b0**2*c**2-&
&a**2*a0**2-b0**2*b**2-c**2*c0**2


alph=(S/S0)**0.5-1.0

bet=SS0/(S*S0)**0.5-1.0

END


subroutine get_energy_MS(energyMS,E_F, N_F, kpp_alpha,mu_ms,a3_ms,a4_ms,b1_ms,b2_ms,&
&area_F_zero,Len_E_zero,Len_E)
implicit none
integer N
parameter(N = 50000)

integer N_F
integer E_F(0:N,1:3)

real*8 kpp_alpha,mu_ms,a3_ms,a4_ms,b1_ms,b2_ms


real*8 area_F_zero(0:N),area_F_temp
real*8 Len_E_zero(0:N),Len_E(0:N)
real*8 energyMS

! %%%%%% MS
real*8 a,b,c,a0,b0,c0,alph,bet

integer i,cc
integer i1,i2,i3

energyMS=0.0
do i=1,N_F
i1=abs(E_F(i,1));i2=abs(E_F(i,2));i3=abs(E_F(i,3));
a=len_E(i1);b=len_E(i2);c=len_E(i3)
a0=len_E_zero(i1);b0=len_E_zero(i2);c0=len_E_zero(i3)

call Get_alph_bet(a,b,c,a0,b0,c0,alph,bet,area_F_temp)
energyMS=energyMS+Kpp_alpha/2.0*(alph**2.0+a3_ms*alph**3+a4_ms*alph**4)*area_F_zero(i)
energyMS=energyMS+mu_ms*(bet+b1_ms*alph*bet+b2_ms*bet**2)*area_F_zero(i)

enddo









END


subroutine MS_force_analytic(Rv,F_V,E_F,V_E, N_V_V,  N_V,   &
& kpp_alpha,mu_ms,a3_ms,a4_ms,b1_ms,b2_ms,area_F_zero,Len_E_zero,Len_E,ms_force_V)
implicit none
integer N
parameter(N = 50000)

integer N_V, N_V_V(0:N)
integer F_V(0:N,1:20),V_E(0:N,1:2),E_F(0:N,1:3)

real*8 kpp_alpha,mu_ms,a3_ms,a4_ms,b1_ms,b2_ms
real*8 rv(0:N,1:3)

real*8 area_F_zero(0:N)
real*8 Len_E_zero(0:N),Len_E(0:N)

! %%%%%% MS
real*8 a,b,c,a0,b0,c0,alph,bet

real*8 MS_Force_V(0:N,1:3)


real*8 R0(1:3),R1(1:3),R2(1:3)
integer i,j,k,i3
integer k_c,k_a,k_b

real*8 DFbet,S,S0,DFalph
real*8 DS(1:3),DPI(1:3)

integer  num_F,I_side(1:3)

MS_Force_V=0.0

do i=1,N_V

r0=rv(i,1:3)

do j=1,N_V_V(i)
num_F=F_V(i,j)
I_side(1:3)=abs(E_F(num_F,1:3))
do k=1,3
if(V_E(I_side(k),1).ne.i.and.V_E(I_side(k),2).ne.i)k_c=k
enddo

i3=0
do k=1,3
if(k.ne.k_c.and.i3.eq.1)then
k_b=k
if(V_E(I_side(k),1).ne.i)r2=rv(V_E(I_side(k),1),1:3)
if(V_E(I_side(k),2).ne.i)r2=rv(V_E(I_side(k),2),1:3)
i3=i3+1
endif
if(k.ne.k_c.and.i3.eq.0)then
k_a=k
if(V_E(I_side(k),1).ne.i)r1=rv(V_E(I_side(k),1),1:3)
if(V_E(I_side(k),2).ne.i)r1=rv(V_E(I_side(k),2),1:3)
i3=i3+1
endif
enddo

a0=len_E_zero(I_side(k_a))
b0=len_E_zero(I_side(k_b))
c0=len_E_zero(I_side(k_c))
a =len_E(I_side(k_a))
b =len_E(I_side(k_b))
c =len_E(I_side(k_c))
call Get_alph_bet_analytic(a,b,c,a0,b0,c0,alph,bet,S,S0)

DFalph=(Kpp_alpha/2.0*(2.0*alph+3.0*a3_ms*alph**2+4.0*a4_ms*alph**3)+mu_ms*b2_ms*bet)&
&*area_F_zero(num_F)
DFbet= mu_ms*(1.0+b1_ms*alph+2*b2_ms*bet)*area_F_zero(NUM_F)

DS=-4*a**2*(r0-r1)-4*b**2*(r0-r2)+4*a**2*(r0-r2)+4*b**2*(r0-r1)+4*c**2*(r0-r1)+4*c**2*(r0-r2)

DPI=-2*a0**2*(r0-r1)-2*b0**2*(r0-r2)+2*(b0**2+c0**2)*(r0-r1)+2*(a0**2+c0**2)*(r0-r2)

MS_Force_V(i,1:3)=MS_Force_V(i,1:3)-DFalph/2.0/(S*S0)**0.5*Ds
MS_Force_V(i,1:3)=MS_Force_V(i,1:3)-DFbet*(DPI/(S*S0)**0.5-(bet+1)/2.0/S*Ds)
enddo ! end of N_V_V

enddo ! end of N_V











END


 ! %%%%%%%%%%%%%%%%%% PM part

subroutine Point_PM_FORCE_NEW(RV,H_modulus,H_zero, PI_K_2_A0,Integral_H_A, H_V,V_V &
&,N_V_V, V_move,area_V, PM_force_point )
implicit none
integer N
parameter(N = 50000)

real*8 rv(0:N,1:3)

integer V_V(0:N,1:20),N_V_V(0:N)
integer V_MOVE

integer I_temp(1:20)

real*8 H_V(0:N),area_V(0:N)
real*8 H_zero,H_modulus,integral_H_A,PI_K_2_A0

real*8 area_V_local



real*8 PM_force_Point(1:3)
real*8 R(1:20,1:3)
real*8 R0(1:3)

real*8 RR(1:20,1:3)
real*8 RR0(1:3)



real*8 delt_r

integer i,j,k,i1, N_t, I_temp2(1:20)


real*8 sum0,sum_H1,sum_H2,H_local,energy_local_1,energy_local_2

real*8 len_E_local(1:20),th_V_local(1:20)


i=v_move

PM_force_Point=0.0
delt_r=0.00001
sum0=Integral_H_A/PI_K_2_A0/2.0



i=v_move
R=0
I_temp=0

I_temp(1:N_V_V(i) )=V_V(i, 1:N_V_V(i) )

sum_H1=H_V(i)*area_V(i)/3.0
do j=1,N_V_V(i)
R(j,1:3)=Rv(I_temp(j),1:3 )
sum_H1=sum_H1+H_V(I_temp(j))*area_V(I_temp(j))/3.0
enddo 


! %%%%%%%%%% Force of x direction
! engery_local_1
R0(1:3)=Rv(i,1:3); 
R0(1)=R0(1)+delt_r
Energy_local_1=0.0
Sum_H2=0.0


N_t=N_V_V(i)

call get_H_local_NEW(R0,R,N_t,   H_local, area_V_local,len_E_local,th_V_local)

sum_H2=sum_H2+H_local*area_V_local/3.0
Energy_local_1=Energy_local_1+H_modulus / 2.0 * (H_local**2-2*H_local*H_zero )* &
&Area_V_local/3.0

do j=1,N_V_V(i)
RR0(1:3)=R(j,1:3)
i1=I_temp(j)  ! vertice No.
I_temp2(1:N_V_V( i1 ) )=V_V(i1, 1:N_V_V( i1 ) )
do k=1, N_V_V(i1)
RR(k,1:3)=Rv(I_temp2(k),1:3)
if(I_temp2(k).eq.i)RR(k,1:3)=R0
enddo ! end of k

N_t=N_V_V(i1)
call get_H_local_NEW(RR0,RR,N_t,   H_local, area_V_local,len_E_local,th_V_local)
sum_H2=sum_H2+H_local*area_V_local/3.0
Energy_local_1=Energy_local_1+H_modulus / 2.0 * (H_local**2-2*H_local*H_zero )* &
&Area_V_local/3.0
enddo  ! end  of j

Energy_local_1=Energy_local_1+PI_K_2_A0*(2*sum0*(sum_H2-sum_H1)+(sum_H2-sum_H1)**2.0)

! energy_local_2
R0(1:3)=Rv(i,1:3); 
R0(1)=R0(1)-delt_r
Energy_local_2=0.0
Sum_H2=0.0


N_t=N_V_V(i)

call get_H_local_NEW(R0,R,N_t,   H_local, area_V_local,len_E_local,th_V_local)


sum_H2=sum_H2+H_local*area_V_local/3.0
Energy_local_2=Energy_local_2+H_modulus / 2.0 * (H_local**2-2*H_local*H_zero )*&
& Area_V_local/3.0

do j=1,N_V_V(i)
RR0(1:3)=R(j,1:3)
i1=I_temp(j)  ! vertice No.
I_temp2(1:N_V_V( i1 ) )=V_V(i1, 1:N_V_V( i1 ) )
do k=1, N_V_V(i1)
RR(k,1:3)=Rv(I_temp2(k),1:3)
if(I_temp2(k).eq.i)RR(k,1:3)=R0
enddo ! end of k

N_t=N_V_V(i1)
call get_H_local_NEW(RR0,RR,N_t,   H_local, area_V_local,len_E_local,th_V_local)
sum_H2=sum_H2+H_local*area_V_local/3.0
Energy_local_2=Energy_local_2+H_modulus / 2.0 * (H_local**2-2*H_local*H_zero )* &
&Area_V_local/3.0
enddo  ! end  of j

Energy_local_2=Energy_local_2+PI_K_2_A0*(2.0*sum0*(sum_H2-sum_H1)+(sum_H2-sum_H1)**2.0)

PM_Force_Point(1)= - (Energy_Local_1-Energy_Local_2)/delt_r/2.0  

! %%%%%%%%%%%%%%

! %%%%%%%%%% Force of y direction
! engery_local_1
R0(1:3)=Rv(i,1:3); 
R0(2)=R0(2)+delt_r
Energy_local_1=0.0
Sum_H2=0.0


N_t=N_V_V(i)

call get_H_local_NEW(R0,R,N_t,   H_local, area_V_local,len_E_local,th_V_local)
sum_H2=sum_H2+H_local*area_V_local/3.0
Energy_local_1=Energy_local_1+H_modulus / 2.0 * (H_local**2-2*H_local*H_zero )* &
&Area_V_local/3.0

do j=1,N_V_V(i)
RR0(1:3)=R(j,1:3)
i1=I_temp(j)  ! vertice No.
I_temp2(1:N_V_V( i1 ) )=V_V(i1, 1:N_V_V( i1 ) )
do k=1, N_V_V(i1)
RR(k,1:3)=Rv(I_temp2(k),1:3)
if(I_temp2(k).eq.i)RR(k,1:3)=R0
enddo ! end of k

N_t=N_V_V(i1)
call get_H_local_NEW(RR0,RR,N_t,   H_local, area_V_local,len_E_local,th_V_local)
sum_H2=sum_H2+H_local*area_V_local/3.0
Energy_local_1=Energy_local_1+H_modulus / 2.0 * (H_local**2-2*H_local*H_zero )* &
&Area_V_local/3.0
enddo  ! end  of j

Energy_local_1=Energy_local_1+PI_K_2_A0*(2*sum0*(sum_H2-sum_H1)+(sum_H2-sum_H1)**2.0)


! energy_local_2
R0(1:3)=Rv(i,1:3); 
R0(2)=R0(2)-delt_r
Energy_local_2=0.0
Sum_H2=0.0


N_t=N_V_V(i)

call get_H_local_NEW(R0,R,N_t,   H_local, area_V_local,len_E_local,th_V_local)
sum_H2=sum_H2+H_local*area_V_local/3.0
Energy_local_2=Energy_local_2+H_modulus / 2.0 * (H_local**2-2*H_local*H_zero )*&
& Area_V_local/3.0

do j=1,N_V_V(i)
RR0(1:3)=R(j,1:3)
i1=I_temp(j)  ! vertice No.
I_temp2(1:N_V_V( i1 ) )=V_V(i1, 1:N_V_V( i1 ) )
do k=1, N_V_V(i1)
RR(k,1:3)=Rv(I_temp2(k),1:3)
if(I_temp2(k).eq.i)RR(k,1:3)=R0
enddo ! end of k

N_t=N_V_V(i1)
call get_H_local_NEW(RR0,RR,N_t,   H_local, area_V_local,len_E_local,th_V_local)
sum_H2=sum_H2+H_local*area_V_local/3.0
Energy_local_2=Energy_local_2+H_modulus / 2.0 * (H_local**2-2*H_local*H_zero )* &
&Area_V_local/3.0
enddo  ! end  of j

Energy_local_2=Energy_local_2+PI_K_2_A0*(2*sum0*(sum_H2-sum_H1)+(sum_H2-sum_H1)**2.0)

PM_Force_Point(2)= - (Energy_Local_1-Energy_Local_2)/delt_r/2.0  

!  %%%%%%%%%%%%%%%%%%%

! %%%%%%%%%% Force of z direction
! engery_local_1
R0(1:3)=Rv(i,1:3); 
R0(3)=R0(3)+delt_r
Energy_local_1=0.0
Sum_H2=0.0


N_t=N_V_V(i)

call get_H_local_NEW(R0,R,N_t,   H_local, area_V_local,len_E_local,th_V_local)
sum_H2=sum_H2+H_local*area_V_local/3.0
Energy_local_1=Energy_local_1+H_modulus / 2.0 * (H_local**2-2*H_local*H_zero )*&
& Area_V_local/3.0

do j=1,N_V_V(i)
RR0(1:3)=R(j,1:3)
i1=I_temp(j)  ! vertice No.
I_temp2(1:N_V_V( i1 ) )=V_V(i1, 1:N_V_V( i1 ) )
do k=1, N_V_V(i1)
RR(k,1:3)=Rv(I_temp2(k),1:3)
if(I_temp2(k).eq.i)RR(k,1:3)=R0
enddo ! end of k

N_t=N_V_V(i1)
call get_H_local_NEW(RR0,RR,N_t,   H_local, area_V_local,len_E_local,th_V_local)
sum_H2=sum_H2+H_local*area_V_local/3.0
Energy_local_1=Energy_local_1+H_modulus / 2.0 * (H_local**2-2*H_local*H_zero )* &
&Area_V_local/3.0
enddo  ! end  of j

Energy_local_1=Energy_local_1+PI_K_2_A0*(2*sum0*(sum_H2-sum_H1)+(sum_H2-sum_H1)**2.0)


! energy_local_2
R0(1:3)=Rv(i,1:3); 
R0(3)=R0(3)-delt_r
Energy_local_2=0.0
Sum_H2=0.0


N_t=N_V_V(i)

call get_H_local_NEW(R0,R,N_t,   H_local, area_V_local,len_E_local,th_V_local)
sum_H2=sum_H2+H_local*area_V_local/3.0
Energy_local_2=Energy_local_2+H_modulus / 2.0 * (H_local**2-2*H_local*H_zero )*&
& Area_V_local/3.0

do j=1,N_V_V(i)
RR0(1:3)=R(j,1:3)
i1=I_temp(j)  ! vertice No.
I_temp2(1:N_V_V( i1 ) )=V_V(i1, 1:N_V_V( i1 ) )
do k=1, N_V_V(i1)
RR(k,1:3)=Rv(I_temp2(k),1:3)
if(I_temp2(k).eq.i)RR(k,1:3)=R0
enddo ! end of k

N_t=N_V_V(i1)
call get_H_local_NEW(RR0,RR,N_t,   H_local, area_V_local,len_E_local,th_V_local)

sum_H2=sum_H2+H_local*area_V_local/3.0
Energy_local_2=Energy_local_2+H_modulus / 2.0 * (H_local**2-2*H_local*H_zero )* &
&Area_V_local/3.0
enddo  ! end  of j

Energy_local_2=Energy_local_2+PI_K_2_A0*(2*sum0*(sum_H2-sum_H1)+(sum_H2-sum_H1)**2.0)

PM_Force_Point(3)= - (Energy_Local_1-Energy_Local_2)/delt_r/2.0  



END

! %%%%%%%%%

subroutine Force_PM_Point_analysis(RV,RV_MOVE,V_move,N_V_V,V_V,E_V,E_V_opposite,&
& F_V,len_E,th_E,vector_F,area_V,area_F,Darea_Point,H_V,H_zero,H_modulus,Integral_H_A,&
&PM_Force_Point)

implicit none
integer N
parameter(N = 50000)

real*8 rv(0:N,1:3),rv_move(1:3)
integer V_move,V_V(0:N,1:20),E_V(0:N,1:20),F_V(0:N,1:20),N_V_V(0:N)

integer E_V_opposite(0:N,1:20)

integer i,j,k,i0,i1,i2,i3,i4,I_temp(1:20),N_t

real*8 Darea_Point(1:20,1:3), Len_E(0:N),th_E(0:N),th0(1:20)
real*8 H_V(0:N),area_V(0:N),area_F(0:N)
real*8 vect_VV(1:20,1:3),vector_F(0:N,1:3)

real*8 Vect_VM1(1:20,1:3),Vect_VM2(1:20,1:3)
real*8 areaM1(1:20),areaM2(1:20)



real*8 thM1(1:20),thM2(1:20),LM1(1:20),LM2(1:20)
real*8 DthM1(1:20),DthM2(1:20)

real*8 PP1(1:3),PP2(1:3),PP3(1:3),PP4(1:3)
real*8 R1R2_fun

real*8 H0,HJ(1:20)
real*8 R0(1:3),RR0(1:20,1:3)

real*8 Dth(1:20), Det_L0(1:20,1:3)
real*8 Det_th0(1:20,1:3),Det_thM1(1:20,1:3),Det_thM2(1:20,1:3)
real*8 Det_HV(1:3), Det_HVJ(1:20,1:3)
real*8 L0(1:20),LJ(1:20,0:2), area0(1:20),areaJ(1:20,1:2),area_VJ(1:20)

real*8 DFHV,DFHVJ(1:20),FHV,FHVJ(1:20)
real*8 Integral_H_A,H_modulus,H_zero


real*8 PM_Force_point(1:3)


i=v_move
N_T=N_V_V(i)
I_temp(1:N_t)=V_V(i,1:N_t)
! Det_HV^alpha
R0=RV_move
H0=H_V(v_move)
do j=1,N_t
i1=abs(E_V(i,j))
i2=mod(j-2+N_t,N_t)+1
i3=mod(j,N_t)+1
HJ(j)=H_V(I_temp(j))
RR0(j,1:3)=RV(I_temp(j),1:3)
L0(j)=Len_E(i1)
LJ(j,1)=L0(j)
LJ(j,0)=Len_E(abs(E_V(i,i2)))
LJ(j,2)=Len_E(abs(E_V(i,i3)))
th0(j)=th_E(i1)
area0(j)=area_F(F_V(i,j))
area_VJ(j)=area_V(V_V(i,j))
areaJ(j,2)=area_F(F_V(i,i2))
areaJ(j,1)=area0(j)
if(abs(th0(j)).lt.0.0001)then
DTH(j)=0.0
else
DTH(j)=-1.0/(1-cos(th0(j))**2.0)**0.5
endif
Vect_VV(j,1:3)=Vector_F(F_V(i,j),1:3)
enddo



! DetH { DetL0j
do j=1,N_t
i1=abs(E_V(i,j))
i2=mod(j-2+N_t,N_t)+1
i3=mod(j,N_t)+1
Det_L0(j,1:3)=(R0-RR0(j,1:3))/L0(j)
PP1(1)=R1R2_fun(RR0(j,1:3), Vect_VV(i2,1:3), 1);
PP1(2)=R1R2_fun(RR0(j,1:3), Vect_VV(i2,1:3), 2);
PP1(3)=R1R2_fun(RR0(j,1:3), Vect_VV(i2,1:3), 3);
PP2(1)=R1R2_fun(Vect_VV(i2,1:3),RR0(i3,1:3),  1);
PP2(2)=R1R2_fun(Vect_VV(i2,1:3),RR0(i3,1:3),  2);
PP2(3)=R1R2_fun(Vect_VV(i2,1:3),RR0(i3,1:3),  3);
PP3(1)=R1R2_fun(RR0(i2,1:3), Vect_VV(j,1:3), 1);
PP3(2)=R1R2_fun(RR0(i2,1:3), Vect_VV(j,1:3), 2);
PP3(3)=R1R2_fun(RR0(i2,1:3), Vect_VV(j,1:3), 3);
PP4(1)=R1R2_fun(Vect_VV(j,1:3),RR0(j,1:3),  1);
PP4(2)=R1R2_fun(Vect_VV(j,1:3),RR0(j,1:3),  2);
PP4(3)=R1R2_fun(Vect_VV(j,1:3),RR0(j,1:3),  3);

Det_th0(j,1:3)=DTH(j)/2.0/areaJ(j,2)/areaJ(j,1)*(PP1+PP2+PP3+PP4)
Det_th0(j,1:3)=Det_th0(j,1:3) -dth(j)*sum(Vect_VV(i2,1:3)*Vect_VV(j,1:3))&
&/areaJ(j,2)/areaJ(j,1)&
&*(Darea_Point(i2,1:3)/areaJ(j,2)**2.0+Darea_Point(j,1:3)/areaJ(j,1)**2.0)
if(th0(j).lt.0.0)Det_th0(j,1:3)=-Det_th0(j,1:3)
enddo



Det_HV=0.0
do j=1,N_t
Det_HV=Det_HV-3.0/2.0/area_V(i)*(L0(j)*Det_th0(j,1:3)+th0(j)*Det_L0(j,1:3))
Det_HV=Det_HV-H0/area_V(i)*Darea_point(j,1:3)/area0(j)
enddo




! M1 M2
do j=1,N_t
i2=mod(j-2+N_t,N_t)+1
thM1(j)=th_E(abs(E_V_opposite(i,i2)))
thM2(j)=th_E(abs(E_V_opposite(i,j)))
LM1(j)=Len_E(abs(E_V_opposite(i,i2)))
LM2(j)=Len_E(abs(E_V_opposite(i,j)))

i0=V_V(i,j)
do k=1,N_V_V(V_V(i,j))
if(abs(E_V(i0,k)).eq.abs(E_V_opposite(i,j)) )i1=mod(k-2+N_V_V(V_V(i,j)),N_V_V(V_V(i,j)))+1
if(abs(E_V(i0,k)).eq.abs(E_V_opposite(i,i2)) )i3=k
enddo

! %%%%%%% M2
i4=F_V(i0,i1)
vect_VM2(j,1:3)= vector_F(i4,1:3)
areaM2(j)=area_F(i4)
! %%%%%%% M1
i4=F_V(i0,i3)
vect_VM1(j,1:3)= vector_F(i4,1:3)
areaM1(j)=area_F(i4)

if(abs(thM1(j)).lt.0.0001)then
DTHM1(j)=0.0
else
DTHM1(j)=-1.0/(1-cos(thM1(j))**2.0)**0.5
endif
if(abs(thM2(j)).lt.0.0001)then
DTHM2(j)=0.0
else
DTHM2(j)=-1.0/(1-cos(thM2(j))**2.0)**0.5
endif

enddo

! det_thM1,Det_thM2
do j=1,N_t
i2=mod(j,N_t)+1
PP1(1)=R1R2_fun(RR0(j,1:3), Vect_VM2(j,1:3), 1);
PP1(2)=R1R2_fun(RR0(j,1:3), Vect_VM2(j,1:3), 2);
PP1(3)=R1R2_fun(RR0(j,1:3), Vect_VM2(j,1:3), 3);
PP2(1)=R1R2_fun(Vect_VM2(j,1:3),RR0(i2,1:3) , 1);
PP2(2)=R1R2_fun(Vect_VM2(j,1:3),RR0(i2,1:3) , 2);
PP2(3)=R1R2_fun(Vect_VM2(j,1:3),RR0(i2,1:3) , 3);
Det_thM2(j,1:3)=DTHM2(j)/2.0/areaM2(j)/area0(j)*(PP1(1:3)+PP2(1:3))
Det_thM2(j,1:3)=Det_thM2(j,1:3)-DTHM2(j)*sum(Vect_VM2(j,1:3)*Vect_VV(j,1:3))&
&/areaM2(j)/area0(j)**3.0*Darea_Point(j,1:3)
if(thM2(j).lt.0.0)Det_thM2(j,1:3)=-Det_thM2(j,1:3)

i2=mod(j-2+N_t,N_t)+1
PP1(1)=R1R2_fun(RR0(i2,1:3), Vect_VM1(j,1:3), 1);
PP1(2)=R1R2_fun(RR0(i2,1:3), Vect_VM1(j,1:3), 2);
PP1(3)=R1R2_fun(RR0(i2,1:3), Vect_VM1(j,1:3), 3);
PP2(1)=R1R2_fun(Vect_VM1(j,1:3),RR0(j,1:3) , 1);
PP2(2)=R1R2_fun(Vect_VM1(j,1:3),RR0(j,1:3) , 2);
PP2(3)=R1R2_fun(Vect_VM1(j,1:3),RR0(j,1:3) , 3);
Det_thM1(j,1:3)=DTHM1(j)/2.0/areaM1(j)/area0(i2)*(PP1(1:3)+PP2(1:3))
Det_thM1(j,1:3)=Det_thM1(j,1:3)-DTHM1(j)*sum(Vect_VM1(j,1:3)*Vect_VV(i2,1:3))&
&/areaM1(j)/area0(i2)**3.0*Darea_Point(i2,1:3)
if(thM1(j).lt.0.0)Det_thM1(j,1:3)=-Det_thM1(j,1:3)
enddo

! Det_HVJ

do j=1,N_t
i1=mod(j-2+N_t,N_t)+1
Det_HVJ(j,1:3)=-3.0/2.0/area_VJ(j)*(LM2(j)*Det_thM2(j,1:3)+L0(j)*Det_th0(j,1:3) &
&+LM1(j)*Det_thM1(j,1:3)+ Th0(j)*Det_L0(j,1:3)) &
&- HJ(j)/area_VJ(j)*(Darea_Point(i1,1:3)/area0(i1)+Darea_Point(j,1:3)/area0(j))
enddo


! ENERGY TERM

FHV=H_modulus/2.0*(H0**2.0-2*H0*H_zero)+Integral_H_A*H0
DFHV=H_modulus*(H0-H_zero)+Integral_H_A
do j=1,N_t
FHVJ(j)=H_modulus/2.0*(HJ(j)**2.0-2*HJ(j)*H_zero)+Integral_H_A*HJ(j)
DFHVJ(j)=H_modulus*(HJ(j)-H_zero)+Integral_H_A
enddo

! Force

PM_Force_point(1:3)=0.0
PM_Force_point=PM_Force_point+DFHV*Det_HV*area_V(i)/3.0

do j=1,N_t
i1=mod(j-2+N_t,N_t)+1
PM_Force_point=PM_Force_point+DFHVJ(j)*Det_HVJ(j,1:3)*area_VJ(j)/3.0
PM_Force_point=PM_Force_point+FHV*Darea_Point(j,1:3)/area0(j)/3.0
PM_Force_point=PM_Force_point+FHVJ(j)*Darea_Point(j,1:3)/area0(j)/3.0
PM_Force_point=PM_Force_point+FHVJ(j)*Darea_Point(i1,1:3)/area0(i1)/3.0
enddo
PM_Force_Point=-PM_Force_Point

END

! %%%%%%%%%%%%%%%%%%

Subroutine Update_New_Parameters(rv,rv_move,V_move,N_V_V,V_V,E_V,F_V,F_E,E_V_opposite,V_E_opposite,&
& Len_E,th_E,vector_F,area_F,area_V,area_tot,vol_F,vol_total,H_V,Pi_K_2_A0,integral_H_A)

implicit none
integer N
parameter(N = 50000)

integer N_V_V(0:N)
integer V_V(0:N,1:20),E_V(0:N,1:20),F_V(0:N,1:20),F_E(0:N,1:20)
integer E_V_opposite(0:N,1:20),V_E_opposite(0:N,1:20)

! For cyclic or label
integer i,j,k,i1,i2,i3,I_temp(1:20),I_temp2(1:20)
integer V_move,N_t

! positions
real*8 rv(0:N,1:3),rv_move(1:3),R0(1:3)
real*8 R(1:20,1:3),RR3(1:3,1:3),R1(1:3),R2(1:3),R3(1:3)
real*8 RJ(1:20,1:3)
! edge variables
real*8 len_E(0:N),th_E(0:N),th_JK,len_JK(1:20),sign_r
!area ones
real*8 vector_F(0:N,1:3),Vector_R0(1:20,1:3),Vector_RJ(1:20,1:3),vector_Local(1:3)
real*8 area_F(0:N),area_V(0:N),area_tot,area_local,area_V_local
real*8 sum_area1_point,sum_area2_point
! volumes
real*8 Vol_F(0:N),Vol_total, sum_vol1_point,sum_vol2_point
! curvatures
real*8 H_V(0:N),H_local,H_JK
real*8 sum_H1_point,sum_H2_point

! other
real*8 Pi_K_2_A0,integral_H_A

! function
real*8 V6_fun



i=V_move

N_t=N_V_V(i)
sum_H1_point=0.0
sum_H2_point=0.0
sum_area2_point=sum(area_F(F_V(i,1:N_t)) )
sum_vol2_point=sum(vol_F( F_V(i,1:N_t)) )
R0=RV(V_move,1:3) 
do j=1,N_t
R(j,1:3)=RV(V_V(i,j),1:3)
I_temp(j)=V_V(i,j)
enddo

! %%%%%%%%%%%%%
!  Get all quantities near vertex V
! %%%%%%%%%%%%%

R(N_t+1,1:3)=R(1,1:3)
R(N_t+2,1:3)=R(2,1:3)
area_V_local=0.0
! update L_E area vol
do j=1,N_t
R1(1:3)=R0(1:3);R2(1:3)=R(j,1:3);R3(1:3)=R(j+1,1:3);
call Get_area_local(R1,R2,R3,area_local, vector_local)
area_F(F_V(i,j))=area_Local
vector_R0(j,1:3)= vector_local(1:3)/area_Local
Vector_F(F_V(i,j),1:3)=Vector_local
area_V_local=area_V_local+area_local

RR3(1,1:3)=r1(1:3)
RR3(2,1:3)=r2(1:3)
RR3(3,1:3)=r3(1:3)
Vol_F(F_V(i,j))=v6_fun(RR3)/6.0

len_jk(j)=(sum((R1-R2)**2.0))**0.5
Len_E(abs(E_V(i,j)))=len_JK(j)
i1=F_E(abs(E_V_opposite(i,j)),1);i2=F_E(abs(E_V_opposite(i,j)),2);
if(i1.ne.F_V(i,j))i3=i1
if(i2.ne.F_V(i,j))i3=i2
Vector_RJ(j,1:3)=Vector_F(i3,1:3)/area_F(i3)
i1=V_E_opposite(abs(E_V_opposite(i,j)),1);i2=V_E_opposite(abs(E_V_opposite(i,j)),2);
if(i1.ne.i)RJ(j,1:3)=Rv(i1,1:3)
if(i2.ne.i)RJ(j,1:3)=Rv(i2,1:3)
enddo


vector_r0(N_t+1,1:3)=vector_R0(1,1:3)
len_Jk(N_t+1)=len_JK(1)
H_JK=0.0

do j=1,N_t
! th_E E_V
sign_r=sum(vector_R0(j,1:3)*vector_R0(j+1,1:3))
if(sign_r.lt.1.0)then
th_jk=dacos(sign_r)
else
th_JK=0.0
endif
i1=mod(j,N_t)+1
sign_r=sum( (vector_R0(j,1:3)-vector_R0(j+1,1:3))*(R(j,1:3)-R(j+2,1:3)) )
if(sign_r.lt.0.0)th_jk=-th_JK
th_E(abs(E_V(i,i1)))=th_JK
H_jk=H_jk+th_jk*len_jk(i1)

! Th_E E_V_opposite
sign_r=sum(vector_R0(j,1:3)*vector_RJ(j,1:3))

if(sign_r.lt.1.0)then
th_jk=dacos(sign_r)
else
th_JK=0.0
endif
sign_r=sum( (vector_R0(j,1:3)-vector_RJ(j,1:3))*(RV_move(1:3)-RJ(j,1:3)) )
if(sign_r.lt.0.0)th_jk=-th_JK
th_E(abs(E_V_opposite(i,j)))=th_JK

enddo


H_local=-H_JK/(area_V_local/3.0)/2.0

sum_H1_point=sum_H1_point+H_local*area_V_local/3.0
sum_H2_point=sum_H2_point+H_v(i)*area_V(i)/3.0

H_V(i)=H_local
area_V(i)=area_V_local


! %%%%%%%%%%
! near VJ
! %%%%%%%

do j=1,N_t

i1=I_temp(j)  ! vertice No.
I_temp2(1:N_V_V( i1 ) )=V_V(i1, 1:N_V_V( i1 ) )
area_V_local=0.0
H_jk=0.0
do k=1, N_V_V(i1)
i2=abs(E_V(i1,k))
area_V_local=area_V_local+area_F(F_V(i1,k))
H_jk=H_jk+th_E(i2)*len_E(i2)
enddo ! end of k
H_local=-H_JK/(area_V_local/3.0)/2.0
sum_H1_point=sum_H1_point+H_local*area_V_local/3.0
sum_H2_point=sum_H2_point+H_v(i1)*area_V(i1)/3.0
area_V(i1)=area_V_local
H_V(i1)=H_local
enddo

sum_area1_point=sum(area_F(F_V(i,1:N_t)))  !updating total area
sum_vol1_point=sum(vol_F(F_V(i,1:N_t)))   !updating total volume

area_tot=area_tot+sum_area1_point-sum_area2_point
vol_total=vol_total+sum_vol1_point-sum_vol2_point
integral_H_A = integral_H_A +(sum_H1_point-sum_H2_point)*2.0*pi_K_2_A0

END








