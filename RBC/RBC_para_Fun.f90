
!  ************************************************
!  Part of 3D vesicle program
! calculate some basic parameters
!                    written by LiJF  Dobb@bbs.fudan.edu.cn
! *************************************************  


! calculate the length of edges and edge vector oriented by E_V



Subroutine Get_Len_E(rv,  V_E, N_E, Len_E, Vector_E)
implicit none
integer N
parameter(N = 50000)
real*8 rv(0:N,1:3)
integer V_E(0:N, 1:2), N_E
real*8 Len_E(0:N)  !length of edges
real*8 Vector_E(0:N,1:3)  !Vector of edges

integer i,i1,i2

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
real*8 R1(1:3),R2(1:3),R3(1:3), area_local,vector_local(1:3)
real*8 R(1:3,1:3), RR_fun
R(1,1:3)=R1(1:3);R(2,1:3)=R2(1:3);R(3,1:3)=R3(1:3)
Vector_local(1)=RR_fun(R,1)/2.0
Vector_local(2)=RR_fun(R,2)/2.0
Vector_local(3)=RR_fun(R,3)/2.0
Area_local= (Vector_local(1)**2+Vector_local(2)**2+Vector_local(3)**2)**0.5
END


! Calculate area vectors and area of small triangles
subroutine Get_area_F( rv, V_F,  N_F,   Area_F, Vector_F, Norm_F)
implicit none
integer N
parameter(N = 50000)

integer N_F
real*8 rV(0:N,1:3)
integer V_F(0:N, 1:3)

real*8  Vector_F(0:N,1:3), Area_F(0:N), Norm_F(0:N,1:3)
real*8 Rt(1:3,1:3)
! function defintion
real*8 RR_fun

integer i,j,i_VF(1:3)


do i=1,N_F
i_VF=V_F(i,1:3)
Rt(1,1:3)=rv(i_VF(1),1:3)
Rt(2,1:3)=rv(i_VF(2),1:3)
Rt(3,1:3)=rv(i_VF(3),1:3)
 do J=1,3
 Vector_F(i,J)=RR_fun(Rt,J)/2.0
 enddo
 Area_F(i)= (Vector_F(i,1)**2 +Vector_F(i,2)**2+Vector_F(i,3)**2 )**(0.5)
 Norm_F(i,1:3)=vector_F(i,1:3)/Area_F(i)
enddo
!Norm_F=0.0

END




! calculate volumn
subroutine Get_vol_F( rv, V_F,  N_F, Vol_F, Vol_total)
implicit none
integer N
parameter(N = 50000)

integer N_F
real*8 rV(0:N,1:3)
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

subroutine Get_H_local_new(R0,R,N_V_V_local,   H_local, area_V_local,len_E_local,th_V_local)

real*8 R(1:20,1:3)
real*8 R0(1:3),R1(1:3),R2(1:3),R3(1:3)
real*8 H_local, area_R0(1:20), vector_R0(1:20,1:3)
real*8 area_local,vector_local(1:3)
real*8 Len_jk(1:20),th_jk
real*8 H_jk,sign_r
real*8 Len_E_local(1:20),th_V_local(1:20)

real*8 area_V_local

integer i,i1 , N_V_V_local

R(N_V_V_local+1,1:3)=R(1,1:3)
R(N_V_V_local+2,1:3)=R(2,1:3)
area_V_local=0.0
do i=1,N_V_V_local
R1(1:3)=R0(1:3);R2(1:3)=R(i,1:3);R3(1:3)=R(i+1,1:3);
call Get_area_local(R1,R2,R3,area_local, vector_local)
area_R0(i)=area_local

vector_R0(i,1:3)= vector_local(1:3)/area_Local
area_V_local=area_V_local+area_local
len_jk(i)=(sum((R1-R3)**2.0))**0.5
enddo
vector_r0(N_V_V_local+1,1:3)=vector_R0(1,1:3)
len_Jk(N_V_V_local+1)=len_JK(1)
H_JK=0.0
do i=1,N_V_V_local
sign_r=sum(vector_R0(i,1:3)*vector_R0(i+1,1:3))
if(sign_r.lt.1.0)then
th_jk=dacos(sign_r)
else
th_JK=0.0
endif
sign_r=sum( (vector_R0(i,1:3)-vector_R0(i+1,1:3))*(R(i,1:3)-R(i+2,1:3)) )
if(sign_r.lt.0.0)th_jk=-th_JK
i1=mod(i,N_v_V_local)+1
th_V_local(i1)=th_JK
len_E_local(i)=len_JK(i)


H_jk=H_jk+th_jk*len_jk(i)
enddo
H_local=-H_JK/(area_V_local/3.0)/2.0

END



! calculate curvature & some basic variations S and V

subroutine Get_H_V_new(rv,N_V,N_F_V, V_V, F_V, E_V,&
&  area_F, Vector_F,  H_V,area_V,Darea_V, Dvol_V,th_E)
implicit none
integer N
parameter(N = 50000)

integer N_V

real*8 rV(0:N,1:3)

! %%%%%%%%%%%%% other labels
! the first kind label
integer E_V(0:N,1:20) ! edges ajecent to a vertice
integer F_V(0:N,1:20), N_F_V(0:N) ! star faces around a vertice


! second kind label
integer V_V(0:N,1:20)

! other input & output

real*8 H_V(0:N)
real*8  Vector_F(0:N,1:3), Area_F(0:N), area_V(0:N)
real*8 R_VV(1:20,1:3), Vect_VV(1:20,1:3)
real*8 ForceV_v(1:3), Darea_V(0:N,1:20, 1:3)  ! \partial area/\partial rv
real*8 DVol_V(0:N,1:3)  ! \partial volume/\partial rv
real*8 Norm_V_temp(1:3)
real*8 at(1:3) , R1R2_fun
real*8 len_jk(1:20),H_jk,th_jk,sign_r
real*8 th_E(0:N)

integer i,j,i1,I_temp(1:20),I_temp2(1:20)



do i=1, N_V
I_temp(1:N_F_V(i))=V_V(i,1:N_F_V(i))
I_temp2(1:N_F_V(i))=F_V(i,1:N_F_V(i))
do j=1,N_F_V(i)
R_VV(j, 1:3)=Rv(I_temp(j) , 1:3)
Vect_VV(j, 1:3)=Vector_F(I_temp2(j) , 1:3)/area_F(I_temp2(j))

enddo
vect_VV(N_F_V(i)+1,1:3)=Vect_VV(1,1:3)
R_VV(N_F_V(i)+1,1:3)=R_VV(1,1:3)
R_VV(N_F_V(i)+2,1:3)=R_VV(2,1:3)
Norm_V_temp=0.0
ForceV_V=0
area_V(i)=0



do j=1,N_F_V(i)
area_V(i)=area_V(i)+Area_F( I_temp2(j) )
Norm_V_temp(1)=Norm_V_temp(1)+R1R2_fun(R_VV(j,1:3), R_VV(j+1,1:3), 1)
Norm_V_temp(2)=Norm_V_temp(2)+R1R2_fun(R_VV(j,1:3), R_VV(j+1,1:3), 2)
Norm_V_temp(3)=Norm_V_temp(3)+R1R2_fun(R_VV(j,1:3), R_VV(j+1,1:3), 3)
at(1)=( R1R2_fun(R_VV(j,1:3), Vect_VV(j,1:3), 1) &
&+R1R2_fun( Vect_VV(j,1:3),R_VV(j+1,1:3), 1)  )/ 2.0
at(2)=( R1R2_fun(R_VV(j,1:3), Vect_VV(j,1:3), 2) &
&+R1R2_fun( Vect_VV(j,1:3),R_VV(j+1,1:3), 2)  )/ 2.0
at(3)=( R1R2_fun(R_VV(j,1:3), Vect_VV(j,1:3), 3) &
&+R1R2_fun( Vect_VV(j,1:3),R_VV(j+1,1:3), 3)  )/ 2.0
Darea_V(i,j,1:3)=at(1:3)
len_JK(j)=(sum((rv(i,1:3)-r_VV(j,1:3))**2.0))**0.5
enddo

len_JK(N_F_V(i)+1)=len_JK(1)
h_JK=0.0
do j=1,N_F_V(i)
sign_r=sum(vect_VV(j,1:3)*vect_VV(j+1,1:3))
if(sign_r.lt.1.0)then
th_jk=dacos(sign_r)
else
th_JK=0.0
endif
i1=mod(j,N_F_V(i))+1
sign_r=sum( (vect_VV(j,1:3)-vect_VV(j+1,1:3))*(R_VV(j,1:3)-R_VV(j+2,1:3)) )
if(sign_r.lt.0.0)th_jk=-th_JK
h_jk=h_jk+th_jk*len_JK(j+1)
th_E(abs(E_V(i,i1)))=th_jk

enddo

H_V(i)=-h_jk/(area_V(i)/3.0)/2.0
Dvol_V(i,1:3)=Norm_v_temp/6.0 
enddo
END


real*8 function Curvature_Energy(H_modulus,H_v, H_zero, Area_V,N_V)
implicit none
integer N
parameter(N = 50000)
integer N_V
real*8 H_modulus, H_zero
real*8   Area_V(0:N),H_v(0:N)
real*8 a
integer i

a=0.0
do i=1,N_V
a=a+ H_modulus / 2.0 * (H_v(i)**2.0-2*H_V(i)*H_zero+H_Zero**2 )* Area_V(i)/3.0
enddo
Curvature_Energy=a
END
!  ************************************************
!  Part of 3D vesicle program
! Some functions
!                    written by LiJF  Dobb@bbs.fudan.edu.cn
! *************************************************  

! input R1, R2 output: |R1\\R2|^alpha
Real*8 Function R1R2_fun(R1,R2,I_c)
real*8 R1(1:3), R2(1:3)
integer I_C

select case(I_c)
case(1)! x (y z)
            
R1R2_fun=R1(2)*R2(3)-R1(3)*R2(2)

case(2) ! y (z x)
        
R1R2_fun=R1(3)*R2(1)-R1(1)*R2(3)

case(3) ! z (x y)

R1R2_fun=R1(1)*R2(2)-R1(2)*R2(1)

case default
write(*,*) 'wrong input within RR_fun'
!pause
end select




END

! input R1, R2, R3 output: ||R1\\R2| \\ R3 |^alpha
Real*8 Function R1R2_R3_fun(R1,R2,R3,I_c)
real*8 R1(1:3), R2(1:3), R3(1:3)
real*8 D1(1:3),R1R2_fun
integer I_C,k

do k=1,3
D1(k)=R1R2_fun(R1,R2,k)
enddo

select case(I_c)
case(1)! x (y z)
     R1R2_R3_fun=R1R2_fun(D1,R3,1)       
case(2) ! y (z x)
     R1R2_R3_fun=R1R2_fun(D1,R3,2)
case(3) ! z (x y)
     R1R2_R3_fun=R1R2_fun(D1,R3,3)
case default
write(*,*) 'wrong input within RR_fun'
!pause
end select

END

! input R1, R2, R3 output: ||R1\\R2| \\ R3 |^alpha
Real*8 Function R1_R2R3_fun(R1,R2,R3,I_c)
real*8 R1(1:3), R2(1:3), R3(1:3)
real*8 D1(1:3),R1R2_fun
integer I_C,k

do k=1,3
D1(k)=R1R2_fun(R2,R3,k)
enddo

select case(I_c)
case(1)! x (y z)
     R1_R2R3_fun=R1R2_fun(R1,D1,1)       
case(2) ! y (z x)
     R1_R2R3_fun=R1R2_fun(R1,D1,2)
case(3) ! z (x y)
     R1_R2R3_fun=R1R2_fun(R1,D1,3)
case default
write(*,*) 'wrong input within RR_fun'
!pause
end select

END


! input R(1:3,1:3), I_c(x,y, or z) output : | |ex+| |ey +| |ez

Real*8 function  RR_fun(R, I_c) 

       ! (P_Num ,Cordinate )
Real*8    R(1:3,  1:3)
integer I_c

select case(I_c)
case(1)! x (y z)
        !  1x                               2x(3 1)                          3x(1 2)    
RR_fun= R(2,2)*R(3,3)-R(2,3)*R(3,2)  + R(3,2)*R(1,3)-R(3,3)*R(1,2) &
&+   R(1,2)*R(2,3)-R(1,3)*R(2,2) 


case(2) ! y (z x)
        !  1y                               2y(3 1)                          3y(1 2)    
RR_fun= R(2,3)*R(3,1)-R(2,1)*R(3,3)  + R(3,3)*R(1,1)-R(3,1)*R(1,3) &
&+   R(1,3)*R(2,1)-R(1,1)*R(2,3) 

case(3) ! z (x y)
        !  1z                               2z(3 1)                          3z(1 2)    
RR_fun= R(2,1)*R(3,2)-R(2,2)*R(3,1)  + R(3,1)*R(1,2)-R(3,2)*R(1,1) &
& +   R(1,1)*R(2,2)-R(1,2)*R(2,1) 
case default
write(*,*) 'wrong input within RR_fun'
!pause
end select


END

!  input R(1:3,1:3) Output |R\\R\\R|  Volume function

real *8 Function V6_fun(R)

Real*8    R(1:3,  1:3)

V6_fun=R(1,1)*( R(2,2)*R(3,3)-R(2,3)*R(3,2) )+R(1,2)*(R(2,3)*R(3,1)&
&- R(2,1)*R(3,3))+R(1,3)*(R(2,1)*R(3,2)-R(2,2)*R(3,1) )


END



      	REAL *8 FUNCTION RANDOM2(seed11)
        real *8 seed11,rM,rJ
        rJ=13807.0
	    rM=2.0**31.0-1.0
        seed11=MOD(seed11*rJ,rM)
        RANDOM2=seed11/rM
        END




      real *8 FUNCTION random1(seed1)
      INTEGER seed1
      INTEGER MBIG,MSEED,MZ
!     REAL MBIG,MSEED,MZ
      REAL FAC
      PARAMETER (MBIG=1000000000,MSEED=161803398,MZ=0,FAC=1./MBIG)
!     PARAMETER (MBIG=4000000.,MSEED=1618033.,MZ=0.,FAC=1./MBIG)
      INTEGER i,iff,ii,inext,inextp,k
      INTEGER mj,mk,ma(55)
!     REAL mj,mk,ma(55)
      SAVE iff,inext,inextp,ma
      DATA iff /0/
	 
      if(seed1.lt.0.or.iff.eq.0)then
        iff=1
        mj=MSEED-iabs(seed1)
        mj=mod(mj,MBIG)
        ma(55)=mj

        mk=1
        do i=1,54
          ii=mod(21*i,55)
          ma(ii)=mk
          mk=mj-mk
          if(mk.lt.MZ)mk=mk+MBIG
          mj=ma(ii)
        end do
        do k=1,4
          do i=1,55
            ma(i)=ma(i)-ma(1+mod(i+30,55))
            if(ma(i).lt.MZ)ma(i)=ma(i)+MBIG
          end do
        end do
        inext=0
        inextp=31
        seed1=1
      endif
      inext=inext+1
      if(inext.eq.56)inext=1
      inextp=inextp+1
      if(inextp.eq.56)inextp=1
      mj=ma(inext)-ma(inextp)

      if(mj.lt.MZ)mj=mj+MBIG
      ma(inext)=mj
      random1=mj*FAC

      END

























