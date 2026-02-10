

Subroutine Capillary_Force(RV_move,R_cap,A_LJ,B_LJ,cutoff_LJ,Cap_Force_point)

implicit none

real*8 RV_move(1:3), R_cap,A_LJ,B_LJ,cutoff_LJ
real*8 Cap_Force_Point(1:3),dvdr
real*8 d,dd

dd=(RV_MOVE(1)**2+RV_MOVE(3)**2)**0.5
d=R_cap-dd
Cap_Force_point=0.0
if(d.lt.cutoff_LJ.and.dd.gt.cutoff_LJ)then
dvdr = -12*A_LJ / d**12 +6*B_LJ / d**6 
Cap_Force_Point(1)= RV_MOVE(1)* (R_cap-dd)/dd* dvdr 
Cap_Force_Point(3)= RV_MOVE(3)* (R_cap-dd)/dd* dvdr
endif

END



Subroutine rotate_RBC(rv,N_V,th,ph,ps)
implicit none
integer N
parameter(N = 50000)
real*8 rv(0:N,1:3)
real*8 rvt(1:3)
integer N_V,i
real*8 th,ph,ps
real*8 A(1:3,1:3)

A(1,1) = cos(ph)*cos(ps)-sin(ph)*cos(th)*sin(ps); A(1,2) = -cos(ph)*sin(ps)-sin(ph)*cos(th)*cos(ps);A(1,3) = sin(ph)*sin(th)
A(2,1) = sin(ph)*cos(ps)+cos(ph)*cos(th)*sin(ps); A(2,2) = -sin(ph)*sin(ps)+cos(ph)*cos(th)*cos(ps);A(2,3) = -cos(ph)*sin(th)
A(3,1) = sin(th)*sin(ps);						  A(3,2) =  sin(th)*cos(ps);						A(3,3) = cos(th)

do i=1, N_V
rvt(1) = sum(A(1,1:3) * RV(i,1:3))
rvt(2) = sum(A(2,1:3) * RV(i,1:3))
rvt(3) = sum(A(3,1:3) * RV(i,1:3))
rv(i,1:3)=rvt(1:3)
enddo
end


Subroutine Max_XYZ(rv,N_V,k,d)
implicit none
integer N
parameter(N = 50000)
real*8 rv(0:N,1:3)

integer N_V,i,k,k1,k2
real*8 dt1,dt2,d

k1=mod(1+k+2,3)+1
k2=mod(2+k+2,3)+1

dt1=0
do i=1, N_V
dt2= (rv(i,k1)**2+rv(i,k2)**2)**0.5
if(dt2.gt.dt1)dt1=dt2
enddo
d=dt1

end









Subroutine Get_nearV(rv,N_V,N_V_V,V_V,Norm_V,cut_off,Near_V,Num_Near_V)
implicit none
integer N
parameter(N = 50000)

real*8 rv(0:N,1:3)
real*8 Rep_Force_V(0:N,1:3),Norm_V(0:N,1:3)
real*8 r1(1:3),r2(1:3)
real*8 cut_off,d,dd

integer N_V,V_V(0:N,1:20),N_V_V(0:N)
integer Near_V(0:N,1:40),Num_Near_V(0:N)
integer i,j,k,k_Near,k1,k2

do i=1,N_V
r1(1:3)=rv(i,1:3)
Num_near_V(i)=0
do j=1,N_V
! get rid of self and near
k1=(j-i)
do k2=1,N_V_V(i)
k1=k1*(j-V_V(i,k2))
enddo

dd = sum(Norm_V(i,1:3)*Norm_V(j,1:3))
if(k1.ne.0.and.dd.lt.0)then


r2(1:3)=rv(j,1:3)
d=(sum((r1-r2)**2.0))**0.5


if(d.lt.cut_off)then
Num_near_V(i)=Num_near_V(i)+1
Near_V(i,Num_near_V(i))=j

endif
endif
enddo
enddo


END



Subroutine Rep_Near_V_point(rv,rv_move,v_move,AMM_LJ,BMM_LJ,&
&rep_const,cut_off,Near_V,Num_Near_V, Rep_Force_V_point)
implicit none
integer N

parameter(N = 50000)
real*8 rv(0:N,1:3)
real*8 Rep_Force_V_point(1:3)
real*8 r1(1:3),r2(1:3)
real*8 d,rep_const,cut_off,AMM_LJ,BMM_LJ,dvdr

integer v_move
real*8 rv_move(1:3)

integer Near_V(0:N,1:40),Num_Near_V(0:N)
integer i,j,k,k_Near,k1,k2

Rep_force_V_point=0.0

i=v_move

do j=1,Num_near_V(i)
k1=Near_V(i,j)
r2=rv(k1,1:3)
d=( sum((rv_move-r2)**2) )**0.5

if(d.lt.cut_off)then
dvdr = -12*AMM_LJ / d**12 +6*BMM_LJ / d**6
rep_Force_V_point(1:3)=rep_Force_V_point(1:3)-dvdr *(rv_move-r2)
!write(*,*)d,dvdr,rep_force_V_point
!pause


endif

enddo


END






