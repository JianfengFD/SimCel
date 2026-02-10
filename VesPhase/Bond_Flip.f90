
!  ************************************************
!  Part of 3D vesicle program
!  Flip some bonds that blong to some obtuse triangles
!                    written by LiJF  Dobb@bbs.fudan.edu.cn
! *************************************************  
! include 'readfe.f90'
! include 'label_system.f90'

subroutine bond_flip(V_F,   E_V,N_E_V,     F_V,N_F_V,    F_E,N_F_E,    V_E_opposite,   N_V,N_E,N_F,V_E,E_F,rv,   Len_E,  Vector_E, I_flip)

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

integer V_E_opposite(0:N,1:20)



integer i,j,k,i1,i2,i3,i4,i5,i0,I_temp(1:20),I_temp2(1:20),I_temp3(1:20),  I_flip




! about Flip
integer Permission_E(0:N)   !  Edges that can be Flipped
real*8 Len_a1,Len_a2,Len_b, Len_c1, Len_c2  !  length of sides of two triangles
integer I_a(1:2),I_b, I_c(1:2)
integer  I_a_n(1:2), I_c_n(1:2)

real*8 costh1,costh2, Switch_R
! Be careful when defining variable with letters "E" "F" "V"
real*8 Len_E(0:N)  !length of edges
real*8 Vector_E(0:N,1:3)  !Vector of edges


Permission_E=1
I_flip=0
do i=1,N_E
if(Permission_E(i).eq.1)then


! Get surrounding edges
i1=F_E(i,1); i2=F_E(i,2)
I_b=i
i3=0
i4=0
do j=1,3
if(abs(E_F(i1,j)).ne.I_b)then
i3=i3+1
I_a(i3)=E_F(i1,j)
endif
if(abs(E_F(i2,j)).ne.I_b)then
i4=i4+1
I_c(i4)=E_F(i2,j)
endif
enddo
a1=Len_E(abs(I_a(1))); a2=Len_E(abs(I_a(2)))
b=Len_E(I_b)
c1=Len_E(abs(I_c(1))); c2=Len_E(abs(I_c(2)))

SwitchR=(a1**2+a2**2-b**2)/a1/a2+(c1**2+c2**2-b**2)/c1/c2

! %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!   IF th1+th2 >pi then Flip the edge

if(SwitchR.lt.-0.0)then
I_flip=I_flip+1


Permission_E(abs(I_a(1)))=0  ! Fix the other adjacent edges
Permission_E(abs(I_a(2)))=0
Permission_E(abs(I_c(1)))=0
Permission_E(abs(I_c(2)))=0

i3=V_E(I_b,2)
if(V_E(abs(I_c(1)),1).eq.i3.or.V_E(abs(I_c(1)),2).eq.i3 )I_a_n(1)=I_c(1)
if(V_E(abs(I_c(2)),1).eq.i3.or.V_E(abs(I_c(2)),2).eq.i3 )I_a_n(1)=I_c(2)

if(V_E(abs(I_a(1)),1).eq.i3.or.V_E(abs(I_a(1)),2).eq.i3 )I_a_n(2)=I_a(1)
if(V_E(abs(I_a(2)),1).eq.i3.or.V_E(abs(I_a(2)),2).eq.i3 )I_a_n(2)=I_a(2)

i3=V_E(I_b,1)
if(V_E(abs(I_a(1)),1).eq.i3.or.V_E(abs(I_a(1)),2).eq.i3 )I_c_n(1)=I_a(1)
if(V_E(abs(I_a(2)),1).eq.i3.or.V_E(abs(I_a(2)),2).eq.i3 )I_c_n(1)=I_a(2)

if(V_E(abs(I_c(1)),1).eq.i3.or.V_E(abs(I_c(1)),2).eq.i3 )I_c_n(2)=I_c(1)
if(V_E(abs(I_c(2)),1).eq.i3.or.V_E(abs(I_c(2)),2).eq.i3 )I_c_n(2)=I_c(2)



E_F(i1,1)=I_a_n(1)
E_F(i1,2)=I_a_n(2)
E_F(i1,3)=I_b

E_F(i2,1)=I_c_n(1)
E_F(i2,2)=I_c_n(2)
E_F(i2,3)=-I_b


V_E(i_b,1)=V_E_opposite(i,1)  !  Edge label Changing
V_E(i_b,2)=V_E_opposite(i,2)

endif

endif
enddo







END















