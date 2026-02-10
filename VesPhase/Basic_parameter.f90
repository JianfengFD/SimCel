
!  ************************************************
!  Part of 3D vesicle program
! calculate some basic parameters
!                    written by LiJF  Dobb@bbs.fudan.edu.cn
! *************************************************  


! calculate the length of edges and edge vector oriented by E_V



Subroutine Get_Len_E(rv,  V_E, N_E, Len_E, Vector_E)

parameter(N = 20000)
real*8 rv(0:N,1:3)
integer V_E(0:N, 1:2), N_E
real*8 Len_E(0:N)  !length of edges
real*8 Vector_E(0:N,1:3)  !Vector of edges

integer i,j,k,i1,i2,i3

real*8 a(1:3),d

do i=1,N_E
i1=V_E(i,1)
i2=V_E(i,2)
a=rv(i2,1:3)-rv(i1,1:3)
d=(a(1)*a(1)+a(2)*a(2)+a(3)*a(3))**(0.5)

Vector_E(i,1:3)=a(1:3)
Len_E(i)=d
enddo


END

! calculate area vector & area of local

subroutine Get_area_local(R1,R2,R3,area_local, vector_local)
integer i,j
real*8 R1(1:3),R2(1:3),R3(1:3), area_local,vector_local(1:3)
real*8 R(1:3,1:3), RR_fun
R(1,1:3)=R1(1:3);R(2,1:3)=R2(1:3);R(3,1:3)=R3(1:3)
Vector_local(1)=RR_fun(R,1)/2.0
Vector_local(2)=RR_fun(R,2)/2.0
Vector_local(3)=RR_fun(R,3)/2.0
Area_local= (Vector_local(1)**2+Vector_local(2)**2+Vector_local(3)**2)**0.5
END


! Calculate area vectors and area of small triangles
subroutine Get_area_F( rv, V_F,  N_F,N_V,   Area_F, Vector_F, Norm_F)

parameter(N = 20000)

integer N_V,N_F
real*8 rV(0:N,1:3),x,y,z
integer V_F(0:N, 1:3)

real*8  Vector_F(0:N,1:3), Area_F(0:N), Norm_F(0:N,1:3)
real*8 s,a,b,c,d
real*8 Rt(1:3,1:3)
! function defintion
real*8 RR_fun

integer i,j,k,i_VF(1:3)


do i=1,N_F
i_VF=V_F(i,1:3)
Rt(1,1:3)=rv(i_VF(1),1:3)
Rt(2,1:3)=rv(i_VF(2),1:3)
Rt(3,1:3)=rv(i_VF(3),1:3)
 do J=1,3
 Vector_F(i,J)=RR_fun(Rt,J)/2.0
 enddo
 Area_F(i)= (Vector_F(i,1)**2 +Vector_F(i,2)**2+Vector_F(i,3)**2 )**(0.5)
 !Norm_F(i,1:3)=Vector_F(i,1:3)/Area_F(i)

enddo
Norm_F=0

END




! calculate volumn
subroutine Get_vol_F( rv, V_F,  N_F,N_V,    Vol_F, Vol_total)
parameter(N = 20000)

integer N_V,N_F
real*8 rV(0:N,1:3),x,y,z
integer V_F(0:N, 1:3)

real*8 Vol_F(0:N), vol_total
real*8 R(1:3,1:3), v6_fun
integer i,i_VF(1:3)

do i=1,N_F
i_VF=V_F(i,1:3)
R(1,1:3)=rv(i_VF(1),1:3)
R(2,1:3)=rv(i_VF(2),1:3)
R(3,1:3)=rv(i_VF(3),1:3)
Vol_F(i)=v6_fun(R)/6.0
enddo
Vol_total=sum( Vol_F(1:N_F) )

END

! curvature of local

subroutine Get_H_local(R0,R,N_V_V_local,   H_local, area_V_local)

real*8 R(1:20,1:3)
real*8 R0(1:3),R1(1:3),R2(1:3),R3(1:3)
real*8 H_local, area_R0(1:20), vector_R0(1:20,1:3)
real*8 area_local,vector_local(1:3)
real*8 ForceV_V(1:3), Vect_V(1:3),Norm_V_temp(1:3)

real*8 area_V_local,R1R2_fun

integer i,j,k , N_V_V_local

real*8 det_1

R(N_V_V_local+1,1:3)=R(1,1:3)
area_V_local=0.0
do i=1,N_V_V_local
R1(1:3)=R0(1:3);R2(1:3)=R(i,1:3);R3(1:3)=R(i+1,1:3);
call Get_area_local(R1,R2,R3,area_local, vector_local)

area_R0(i)=area_local
vector_R0(i,1:3)= vector_local(1:3)
area_V_local=area_V_local+area_local
enddo

Norm_V_temp=0.0
ForceV_V=0

do j=1,N_V_V_local
Norm_V_temp(1)=Norm_V_temp(1)+R1R2_fun(R(j,1:3), R(j+1,1:3), 1)
Norm_V_temp(2)=Norm_V_temp(2)+R1R2_fun(R(j,1:3), R(j+1,1:3), 2)
Norm_V_temp(3)=Norm_V_temp(3)+R1R2_fun(R(j,1:3), R(j+1,1:3), 3)
ForceV_V(1)=ForceV_V(1)+( R1R2_fun(R(j,1:3), vector_R0(j,1:3), 1) +R1R2_fun( vector_R0(j,1:3),R(j+1,1:3), 1)  )/ 2.0/ area_R0(j)
ForceV_V(2)=ForceV_V(2)+( R1R2_fun(R(j,1:3), vector_R0(j,1:3), 2) +R1R2_fun( vector_R0(j,1:3),R(j+1,1:3), 2)  )/ 2.0/ area_R0(j)
ForceV_V(3)=ForceV_V(3)+( R1R2_fun(R(j,1:3), vector_R0(j,1:3), 3) +R1R2_fun( vector_R0(j,1:3),R(j+1,1:3), 3)  )/ 2.0/ area_R0(j)
enddo

Norm_v_temp=Norm_V_temp/6.0
ForceV_V=-ForceV_V/3.0  
det_1=(ForceV_V(1)*Norm_V_temp(1)+ForceV_V(2)*Norm_V_temp(2)+ForceV_V(3)*Norm_V_temp(3))
if(abs(det_1).lt.0.000000001)then
H_local=0.0
else
H_local=1.5*(ForceV_V(1)**2+ForceV_V(2)**2+ForceV_V(3)**2)/(ForceV_V(1)*Norm_V_temp(1)+ForceV_V(2)*Norm_V_temp(2)+ForceV_V(3)*Norm_V_temp(3))
endif
H_local=2*H_local


END



! calculate curvature & some basic variations S and V

subroutine Get_H_V(rv, N_F,N_V,N_E,N_F_V, V_F, V_E, E_F,  V_V, F_V, E_V,E_V_opposite,  area_F, Vector_F, norm_F, H_V,area_V,Darea_V, Dvol_V)

parameter(N = 20000)

integer N_V,N_E,N_F

real*8 rV(0:N,1:3),x,y,z,a

integer V_E(0:N, 1:2)
integer E_F(0:N, 1:3)
! %%%%%%%%%%%%% other labels
! the first kind label
integer V_F(0:N,1:3)  ! vertices around a face
integer E_V(0:N,1:20), N_E_V(0:N) ! edges ajecent to a vertice
integer F_V(0:N,1:20), N_F_V(0:N) ! star faces around a vertice
integer F_E(0:N,1:20), N_F_E(0:N)


! second kind label
integer V_V(0:N,1:20), N_V_V(0:N)
integer E_V_opposite(0:N,1:20), N_E_V_opposite(0:N)
integer V_E_opposite(0:N,1:20)

! other input & output

real*8 H_V(0:N)
real*8  Vector_F(0:N,1:3), Area_F(0:N), Norm_F(0:N,1:3), area_V(0:N)
real*8 R_VV(1:20,1:3), Vect_VV(1:20,1:3)
real*8 ForceV_v(1:3), Darea_V(0:N,1:20, 1:3)  ! \partial area/\partial rv
real*8 DVol_V(0:N,1:3)  ! \partial volume/\partial rv
real*8 Norm_V_temp(1:3)
real*8 at(1:3) , R1R2_fun

integer i,j,k,i1,i2,i3,i4,i5,i0,I_temp(1:20),I_temp2(1:20),I_temp3(1:20)

real*8 det_1


do i=1, N_V
I_temp(1:N_F_V(i))=V_V(i,1:N_F_V(i))
I_temp2(1:N_F_V(i))=F_V(i,1:N_F_V(i))
do j=1,N_F_V(i)
R_VV(j, 1:3)=Rv(I_temp(j) , 1:3)
Vect_VV(j, 1:3)=Vector_F(I_temp2(j) , 1:3)
enddo
R_VV(N_F_V(i)+1,1:3)=R_VV(1,1:3)
Norm_V_temp=0.0
ForceV_V=0
area_V(i)=0

do j=1,N_F_V(i)
area_V(i)=area_V(i)+Area_F( I_temp2(j) )
Norm_V_temp(1)=Norm_V_temp(1)+R1R2_fun(R_VV(j,1:3), R_VV(j+1,1:3), 1)
Norm_V_temp(2)=Norm_V_temp(2)+R1R2_fun(R_VV(j,1:3), R_VV(j+1,1:3), 2)
Norm_V_temp(3)=Norm_V_temp(3)+R1R2_fun(R_VV(j,1:3), R_VV(j+1,1:3), 3)
at(1)=( R1R2_fun(R_VV(j,1:3), Vect_VV(j,1:3), 1) +R1R2_fun( Vect_VV(j,1:3),R_VV(j+1,1:3), 1)  )/ 2.0
at(2)=( R1R2_fun(R_VV(j,1:3), Vect_VV(j,1:3), 2) +R1R2_fun( Vect_VV(j,1:3),R_VV(j+1,1:3), 2)  )/ 2.0
at(3)=( R1R2_fun(R_VV(j,1:3), Vect_VV(j,1:3), 3) +R1R2_fun( Vect_VV(j,1:3),R_VV(j+1,1:3), 3)  )/ 2.0
ForceV_V(1)=ForceV_V(1)+at(1)/ area_F(I_temp2(j))
ForceV_V(2)=ForceV_V(2)+at(2)/ area_F(I_temp2(j))
ForceV_V(3)=ForceV_V(3)+at(3)/ area_F(I_temp2(j))
Darea_V(i,j,1:3)=at(1:3)
enddo

Norm_v_temp=Norm_V_temp/6.0
ForceV_V=-ForceV_V/3.0  
det_1=(ForceV_V(1)*Norm_V_temp(1)+ForceV_V(2)*Norm_V_temp(2)+ForceV_V(3)*Norm_V_temp(3))
if(abs(det_1).lt.0.00000001)then
H_V(i)=0.0
else
H_V(i)=2*1.5*(ForceV_V(1)**2+ForceV_V(2)**2+ForceV_V(3)**2)/(ForceV_V(1)*Norm_V_temp(1)+ForceV_V(2)*Norm_V_temp(2)+ForceV_V(3)*Norm_V_temp(3))
endif
Dvol_V(i,1:3)=Norm_v_temp 

enddo

END



subroutine Get_inter_phase(R,area_local, fit, E_inter,b)
real*8 R(1:3,1:3),area_local, fit(1:3)
real*8 tr(1:3,1:3),dfi(1:3), E1,E_inter,b
E1=0.0
tr(1,1:3)=R(2,1:3)-R(1,1:3)
tr(2,1:3)=R(3,1:3)-R(2,1:3)
tr(3,1:3)=R(1,1:3)-R(3,1:3)
dfi(1)=fit(2)-fit(1)
dfi(2)=fit(3)-fit(2)
dfi(3)=fit(1)-fit(3)

E1=E1+b/24.0*( dfi(1)**2* (sum(tr(2,1:3)**2)+ sum(tr(3,1:3)**2)) - 2*dfi(1)*dfi(2)* sum(tr(1,1:3)*tr(2,1:3) ) )/area_local
E1=E1+b/24.0*( dfi(2)**2* (sum(tr(1,1:3)**2)+ sum(tr(3,1:3)**2)) - 2*dfi(2)*dfi(3)* sum(tr(2,1:3)*tr(3,1:3) ) )/area_local
E1=E1+b/24.0*( dfi(3)**2* (sum(tr(1,1:3)**2)+ sum(tr(2,1:3)**2)) - 2*dfi(3)*dfi(1)* sum(tr(3,1:3)*tr(1,1:3) ) )/area_local
E_inter=E1
END


! laplace_tot
subroutine Laplace_global(rv, area_F,area_V, N_V, V_V, F_V, N_F_V, L_glb)

parameter(N = 20000)

integer N_V

real*8 rV(0:N,1:3)

integer V_V(0:N, 1:20)
integer F_V(0:N, 1:20),N_F_V(0:N)

real*8 R0(1:3), R(1:20,1:3)
real*8 area_F(0:N), area_V(0:N)

real*8 area_FV(1:20), area_VV
real*8 L_laplace_local(1:20)
real*8 L_glb(0:N,1:20)
integer N_T,i,j,k

L_glb=0.0

do i=1,N_V
area_FV=0.0
R=0.0
R0=Rv(i,1:3)
N_T = N_F_V(i)
do j=1,N_T
R(j,1:3)=RV(V_V(i,j), 1:3)
area_FV(j)=area_F(F_V(i,j))
enddo
area_VV=area_V(i)
call Laplace_Local(R0,R,area_FV, area_VV, L_laplace_local, N_t)
L_glb(i,1:N_T)=L_laplace_local(1:N_T)
enddo
END

! laplace_local
subroutine Laplace_Local(R0,R,area_FV, area_VV, L_laplace_local, N_t)
real*8 R0(1:3), R(1:20,1:3)
real*8 area_FV(1:20), area_VV
real*8 L_laplace_local(1:20)
real*8 a(1:3),b(1:3),c(1:3)
real*8 L_t
integer N_T,i,j,k


R(N_T+1,1:3)=R(1,1:3)
area_FV(N_T+1)=area_FV(1)
L_laplace_local=0.0
do i=1,N_T
 a=R(i,1:3)-R0
 b=R(i+1,1:3)-R(i,1:3)
 c=R0-R(i+1,1:3)
L_laplace_local(i)=L_laplace_local(i) - 1.0/4.0*sum(c(1:3)*b(1:3) )/area_FV(i)/(area_VV/3.0)
L_laplace_local(i+1)=L_laplace_local(i+1) - 1.0/4.0*sum(a(1:3)*b(1:3) )/area_FV(i)/(area_VV/3.0)
enddo
L_laplace_local(1)=L_laplace_local(1)+L_laplace_local(N_T+1)



END







subroutine sort_num2(data0, L, order_data1,order_data2)
integer L
real*8 data1(1:1000)
real*8 data0(1:1000)
integer i,j,k
integer order_data1(1:1000)
integer order_data2(1:1000)
integer time_k(1:1000)
real*8 a

data1=data0
time_k=0
do i=1,L
a=data1(i)
k=1
do j=1,L
if(i.ne.j)then
if(a.ge.data1(j))k=k+1
endif
enddo
data0(k)=a

time_k(k)=time_k(k)+1
order_data1(i)=k+1-time_k(k)
order_data2(k+1-time_k(k))=i

enddo



END



























