
!  ************************************************
!  Part of 3D vesicle program
!  Vertex averaging
!                    written by LiJF  Dobb@bbs.fudan.edu.cn
! *************************************************
! include 'Vertex_average.f90'


subroutine Vertex_averaging(rv,N_V,N_F,F_V,V_F,V_V,N_F_V , fi,order_fi,area_tot,area_tot_zero)

parameter(N = 20000)

real*8 rv(0:N,1:3),rv0(0:N,1:3)
real*8 fi(0:N),fi1,fi2,fi3,order_fi
real*8 fi_around(1:20)
real*8 sum1,sum2

! %%%%%%%%%%%%%%%%%%%%%
integer N_V,N_F,N_E
integer F_V(0:N,1:20), V_F(0:N,1:3), V_V(0:N,1:20)
integer N_F_V(0:N)

! %%%%%%%%%%%%%%
real*8 area_tot, area_tot_zero
real*8  Vector_F(0:N,1:3), Area_F(0:N), Norm_F(0:N,1:3)
real*8 area_V(0:N)
!
real*8 X_centroid(1:20,1:3), area_around(1:20)
real*8 Norm_F1(1:20,1:3), Norm_V(1:3)
real*8 lamda
real*8 rv_new(1:3),rv_old(1:3)
real*8 rv_new1(1:3), rv_new2(1:3)

! %%%%%%%%%%%%%%

integer i,j,k,i1,i2,i3,i4,i5,i0,I_temp(1:20),I_temp2(1:20),I_temp3(1:20)
integer N_t

! %%%%%%%%%%
real*8 R(1:3,1:3), R1(1:3), R2(1:3), R3(1:3)
real*8 a1,a2,a3,v1(1:3),v2(1:3),v3(1:3)

real*8 len_V(1:20),len_VV


! %%%%%%%%%%%

do i=1, N_V

N_t=N_F_V(i)
I_temp(1:20)=F_V(i,1:20);I_temp(N_t+1)=I_temp(1)
I_temp2(1:20)=V_V(i,1:20);I_temp2(N_t+1)=I_temp2(1)
v2=0.0
Norm_V=0.0

rv_old=rv(i,1:3)

do j=1,N_t
i1=i; i2=I_temp2(j); i3=I_temp2(j+1)
R1(1:3)=Rv(i1,1:3);R2(1:3)=Rv(i2,1:3);R3(1:3)=Rv(i3,1:3);
R(1,1:3)=R1(1:3);R(2,1:3)=R2(1:3);R(3,1:3)=R3(1:3);
call Get_area_local(R1,R2,R3,a1, v1)
area_around(j)=a1
Norm_F1(j,1:3)=V1
X_centroid(j,1:3)=(R1+R2+R3)/3.0
v2=v2+a1*X_centroid(j,1:3)
Norm_V=Norm_V+v1
enddo ! end j
rv_new=v2/sum(area_around(1:N_t))
lamda=(sum(rv_new(1:3)*Norm_V)-sum(rv(i,1:3)*Norm_V) )/sum(Norm_V*Norm_V)
rv(i,1:3)=rv_new-lamda*Norm_V

! %%%% Get new fi
fi1=0.0

len_VV=(sum((rv(i,1:3)-rv_old(1:3))**2) )**0.5
do j=1,N_t
i1=I_temp2(j)
R1(1:3)=Rv(i1,1:3)
len_V(j)=(sum((rv(i,1:3)-r1(1:3))**2) )**0.5
fi1=fi1+fi(i1)/len_V(j)**2.0
enddo ! end j
if(len_VV.lt.0.00001)then

else
fi1=fi1+fi(i)/len_VV**2.0
fi(i)=fi1/(sum(1/len_V(1:N_t)**2.0)+1/len_VV**2.0)
endif


enddo ! end i

! %%%%%%%%%%%%%%%%%%%%%%% reconsider fi
call Get_area_F( rv, V_F,  N_F,N_V,   Area_F, Vector_F, Norm_F)

do i=1,N_V
area_V(i)=0.0
N_t=N_F_V(i)
I_temp(1:20)=F_V(i,1:20)
do j=1,N_t
area_V(i)=area_V(i)+area_F(I_temp(j))/3.0
enddo
enddo

area_tot=sum(area_V(1:N_V))
! area_tot_zero=area_tot

sum1=sum( fi(1:N_V)*area_V(1:N_V) )/area_tot
do i=1,N_V
fi(i)=fi(i)-(sum1-order_fi)
enddo

END



! not yet complete
subroutine Vertex_averaging_fi_conserve(rv,N_V,N_F,F_V,V_F,V_V,N_F_V , fi,order_fi,area_tot,area_tot_zero)

parameter(N = 20000)

real*8 rv(0:N,1:3),rv0(0:N,1:3)
real*8 fi(0:N),fi1,fi2,fi3,order_fi
real*8 fi_around(1:20)
real*8 sum1,sum2

! %%%%%%%%%%%%%%%%%%%%%
integer N_V,N_F,N_E
integer F_V(0:N,1:20), V_F(0:N,1:3), V_V(0:N,1:20)
integer N_F_V(0:N)

! %%%%%%%%%%%%%%
real*8 area_tot, area_tot_zero
real*8  Vector_F(0:N,1:3), Area_F(0:N), Norm_F(0:N,1:3)
real*8 area_V(0:N)
!
real*8 X_centroid(1:20,1:3), area_around(1:20)
real*8 area_around_new(1:20)
real*8 Norm_F1(1:20,1:3), Norm_V(1:3)
real*8 lamda
real*8 rv_new(1:3),rv_old(1:3)
real*8 rv_new1(1:3), rv_new2(1:3)

! %%%%%%%%%%%%%%

integer i,j,k,i1,i2,i3,i4,i5,i0,I_temp(1:20),I_temp2(1:20),I_temp3(1:20)
integer N_t

! %%%%%%%%%%
real*8 R(1:3,1:3), R1(1:3), R2(1:3), R3(1:3)
real*8 a1,a2,a3,v1(1:3),v2(1:3),v3(1:3)

real*8 len_V(1:20),len_VV


! %%%%%%%%%%%

do i=1, N_V

N_t=N_F_V(i)
I_temp(1:20)=F_V(i,1:20);I_temp(N_t+1)=I_temp(1)
I_temp2(1:20)=V_V(i,1:20);I_temp2(N_t+1)=I_temp2(1)
v2=0.0
Norm_V=0.0

rv_old=rv(i,1:3)

do j=1,N_t
i1=i; i2=I_temp2(j); i3=I_temp2(j+1)
R1(1:3)=Rv(i1,1:3);R2(1:3)=Rv(i2,1:3);R3(1:3)=Rv(i3,1:3);
call Get_area_local(R1,R2,R3,a1, v1)
area_around(j)=a1
Norm_F1(j,1:3)=V1
X_centroid(j,1:3)=(R1+R2+R3)/3.0
v2=v2+a1*X_centroid(j,1:3)
Norm_V=Norm_V+v1
enddo ! end j
rv_new=v2/sum(area_around(1:N_t))
lamda=(sum(rv_new(1:3)*Norm_V)-sum(rv(i,1:3)*Norm_V) )/sum(Norm_V*Norm_V)
rv(i,1:3)=rv_new-lamda*Norm_V

! %%%% Get new fi


! %% fi0_n= (sum_i  fi_i(A_i+A_{i+1} ) - sum_i fi_i{}'.. +fi_0 sum A_i )/sum A_i'
do j=1,N_t
i1=i; i2=I_temp2(j); i3=I_temp2(j+1)
R1(1:3)=Rv(i1,1:3);R2(1:3)=Rv(i2,1:3);R3(1:3)=Rv(i3,1:3);
call Get_area_local(R1,R2,R3,a1, v1)
area_around_New(j)=a1
enddo ! end j
fi1=0.0
Area_around(N_t+1)=Area_around(1)
Area_around_New(N_t+1)=Area_around_New(1)
do j=1,N_t
fi1=fi1+fi(I_temp2(j))*(Area_around(j)+Area_around(j+1)-Area_around_new(j)-Area_around_new(j+1))
enddo
fi(i)=(fi1+fi(i)*sum(Area_around(1:N_t)) )/sum(area_around_new(1:N_t))

enddo ! end i

! %%%%%%%%%%%%%%%%%%%%%%% reconsider fi
call Get_area_F( rv, V_F,  N_F,N_V,   Area_F, Vector_F, Norm_F)

do i=1,N_V
area_V(i)=0.0
N_t=N_F_V(i)
I_temp(1:20)=F_V(i,1:20)
do j=1,N_t
area_V(i)=area_V(i)+area_F(I_temp(j))/3.0
enddo
enddo

area_tot=sum(area_V(1:N_V))

sum1=sum( fi(1:N_V)*area_V(1:N_V) )/area_tot
do i=1,N_V
fi(i)=fi(i)-(sum1-order_fi)
enddo

END






