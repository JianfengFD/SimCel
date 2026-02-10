
!  ************************************************
!  Part of 3D vesicle program
! label the surface
!                    written by LiJF  Dobb@bbs.fudan.edu.cn
! *************************************************  
!! include 'readfe.f90'


 subroutine weep(V_F,   E_V,N_E_V,     F_V,N_F_V,    F_E,N_F_E,    V_V,N_V_V,   E_V_opposite,N_E_V_opposite,   V_E_opposite,  N_V,N_E,N_F,V_E,E_F,rv,fi)

parameter(N = 20000)

integer N_V,N_E,N_F
real*8  rv(0:N,1:3),fi(0:N),rv0(0:N,1:3),fi0(0:N)
integer V_E(0:N, 1:2)
integer E_F(0:N, 1:3)
! %%%%%%%%%%%%% other labels
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

integer i,j,k,i1,i2,i3,i4,i5,i0,I_temp(0:20),I_temp2(0:20),I_temp3(1:20),I_temp4(1:20)
integer N_t,N_t1,k1,k2,k3, i31,i32

! New label
integer V_Ec(0:N,1:2), E_Fc(0:N,1:3)
! transformation
integer Ic_V(0:N),Ic_E(0:N),Ic_F(0:N),Irc_E(0:N),Irc_F(0:N),Irc_V(0:N)
integer N_Vc,N_Ec,N_Fc
integer v1,v2,v3,vd,vu,kd,ku,kE,kF,v1t,v2t, Ev1,Ev2
integer vdi(1:20),vui(1:20),vdi_new(1:20),vui_new(1:20)
integer Id_V1,Iu_V1,Edi(1:20),Eui(1:20),Fdi(1:20),Fui(1:20)
integer Id_v2,Iu_v2,Id_v3,Iu_v3
integer Id_v2e,Iu_v2e,Id_v3e,Iu_v3e
integer VIntro_v2,I4_Intro_v2,VIntro_v3,I4_Intro_v3
integer Ed_op(1:20),Eu_op(1:20),I_Edop,I_Euop

integer order_Ed(1:20),Order_Eu(1:20),order_Vd(1:20),Order_Vu(1:20)
integer vdi_order(1:20),vui_order(1:20)

integer E_del(1:40),F_del(1:40), V_del(1:40)
integer neighbor_V(1:40),total_n



!%%%%%%%%%%%%%%%%%%%%%

integer Find_k

! %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! Find three-point neck
! %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Find_k=0
v3=0
Ev1=0;Ev2=0
do i=1,N_v
v1t=i
N_t=N_V_V(i)
I_temp(1:N_t)=V_V(i,1:N_t)
I_temp3(1:N_t)=E_V(i,1:N_t)
I_temp4(1:N_t)=F_V(i,1:N_t)
I_temp(0)=I_temp(N_t)
I_temp(N_t+1)=I_temp(1)
do j=1,N_t
i1=I_temp(j)
v2t=i1
N_t1=N_V_V(i1)
I_temp2(1:N_t1)=V_V(i1,1:N_t1)
do k=1,N_t1
do i2=1,N_t
k1=j
k2=j+1
k3=j-1
if(j.eq.1)k3=N_t
if(j.eq.N_t)k2=1
i31=i2-1;if(i2.eq.1)i31=N_t
i32=i2+1;if(i2.eq.N_t)i32=1
if((i2-k1)*(i2-k2)*(i2-k3).ne.0)then
if(I_temp(i2).eq.I_temp2(k))then
Find_k=1
v1=v1t;Ev1=abs(I_temp3(j))
v2=v2t;Ev2=abs(I_temp3(i2))
v3=I_temp2(k)
Id_v1=0;Iu_v1=0
i4=0
do i3=1,N_t

if(i3.eq.j.or.i3.eq.i2)i4=mod(i4+1,2)
if(i3.ne.j.and.i3.ne.i2)then
if(i4.eq.0)then
Id_v1=Id_v1+1
vdi(Id_v1)=I_temp(i3)
Edi(Id_v1)=I_temp3(i3)
Fdi(Id_v1)=I_temp4(i3)
else
Iu_v1=Iu_v1+1
vui(Iu_v1)=I_temp(i3)
Eui(Iu_v1)=I_temp3(i3)
Fui(Iu_v1)=I_temp4(i3)
endif
if(i3.eq.k2)then
VIntro_V2=I_temp(i3)
Id_intro_v2=i4
endif
if(i3.eq.i32)then
VIntro_V3=I_temp(i3)
Id_intro_v3=i4
endif

endif
enddo
endif
endif
enddo ! i2
enddo ! end of k
enddo ! end of j

enddo ! end of i

write(*,*)Ev1,Ev2,find_k



if (find_k.eq.0)return
! end of find

! %%%%%%%%%%%%%%%%
! %%%%%%%%%%%%%%%% find upper vertices and lower vertices
Neighbor_V(1:Id_v1)=vdi(1:Id_v1)
Neighbor_V(Id_v1+1:Id_v1+Iu_v1)=vui(1:Iu_v1)
total_n=Id_v1+Iu_v1
N_t=N_V_V(v2)
I_temp(1:N_t)=V_V(v2,1:N_t)
I_temp2(1:N_t)=E_V(v2,1:N_t)
I_temp3(1:N_t)=F_V(v2,1:N_t)
k=0
do j=1,N_t
if(I_temp(j).eq.Vintro_v2)k=j
if(I_temp(j).eq.v1)i2=j
if(I_temp(j).eq.v3)i3=j
enddo

i4=Id_intro_v2
Id_v2=0;Iu_v2=0
Id_v2e=0;Iu_v2e=0

do j=k,1,-1
i5=1
do k1=1,total_n
if(I_temp(j).eq.Neighbor_V(k1))i5=0

enddo
if(I_temp(j).eq.V1.or.I_temp(j).eq.V3)then
i4=mod(i4+1,2);i5=0
endif
if(i5.eq.1)then
if(i4.eq.0)then
Id_v2=Id_v2+1;
vdi(Id_v1+Id_v2)=I_temp(j)
endif
if(i4.eq.1)then
Iu_v2=Iu_v2+1;
vui(Iu_v1+Iu_v2)=I_temp(j)
endif
endif
if(I_temp(j).ne.V1.and.I_temp(j).ne.V3)then
if(i4.eq.0)then
Id_v2e=Id_v2e+1
Edi(Id_v1+Id_v2e)=I_temp2(j)
Fdi(Id_v1+Id_v2e)=I_temp3(j)
endif
if(i4.eq.1)then
Iu_v2e=Iu_v2e+1
Eui(Iu_v1+Iu_v2e)=I_temp2(j)
Fui(Iu_v1+Iu_v2e)=I_temp3(j)
endif
endif

enddo

i4=Id_intro_v2
do j=k+1,N_t
i5=1
do k1=1,total_n
if(I_temp(j).eq.Neighbor_V(k1))i5=0

enddo
if(I_temp(j).eq.V1.or.I_temp(j).eq.V3)then
i4=mod(i4+1,2);i5=0
endif
if(i5.eq.1)then
if(i4.eq.0)then
Id_v2=Id_v2+1;
vdi(Id_v1+Id_v2)=I_temp(j)
endif
if(i4.eq.1)then
Iu_v2=Iu_v2+1;
vui(Iu_v1+Iu_v2)=I_temp(j)
endif
endif
if(I_temp(j).ne.V1.and.I_temp(j).ne.V3)then
if(i4.eq.0)then
Id_v2e=Id_v2e+1
Edi(Id_v1+Id_v2e)=I_temp2(j)
Fdi(Id_v1+Id_v2e)=I_temp3(j)
endif
if(i4.eq.1)then
Iu_v2e=Iu_v2e+1
Eui(Iu_v1+Iu_v2e)=I_temp2(j)
Fui(Iu_v1+Iu_v2e)=I_temp3(j)
endif
endif

enddo

! %%%%%%%%%%%%%%% v3
Neighbor_V(1:Id_v1+Id_v2)=vdi(1:Id_v1+Id_v2)
Neighbor_V(Id_v1+Id_v2+1:Id_v1+Id_v2+Iu_v1+Iu_v2)=vui(1:Iu_v1+Iu_v2)
total_n=Id_v1+Id_v2+Iu_v1+Iu_v2
N_t=N_V_V(v3)
I_temp(1:N_t)=V_V(v3,1:N_t)
I_temp2(1:N_t)=E_V(v3,1:N_t)
I_temp3(1:N_t)=F_V(v3,1:N_t)
k=0
do j=1,N_t
if(I_temp(j).eq.Vintro_v3)k=j
if(I_temp(j).eq.v1)i2=j
if(I_temp(j).eq.v2)i3=j
enddo

i4=Id_intro_v3
Id_v3=0;Iu_v3=0
Id_v3e=0;Iu_v3e=0

do j=k,1,-1
i5=1
do k1=1,total_n
if(I_temp(j).eq.Neighbor_V(k1))i5=0

enddo
if(I_temp(j).eq.V1.or.I_temp(j).eq.V2)then
i4=mod(i4+1,2);i5=0
endif
if(i5.eq.1)then
if(i4.eq.0)then
Id_v3=Id_v3+1;
vdi(Id_v1+Id_v2+Id_v3)=I_temp(j)
endif
if(i4.eq.1)then
Iu_v2=Iu_v2+1;
vui(Iu_v1+Iu_v2+Iu_v3)=I_temp(j)
endif
endif
if(I_temp(j).ne.V1.and.I_temp(j).ne.V2)then
if(i4.eq.0)then
Id_v3e=Id_v3e+1
Edi(Id_v1+Id_v2e+Id_v3e)=I_temp2(j)
Fdi(Id_v1+Id_v2e+Id_v3e)=I_temp3(j)
endif
if(i4.eq.1)then
Iu_v3e=Iu_v3e+1
Eui(Iu_v1+Iu_v2e+Iu_v3e)=I_temp2(j)
Fui(Iu_v1+Iu_v2e+Iu_v3e)=I_temp3(j)
endif
endif

enddo

i4=Id_intro_v3
do j=k+1,N_t
i5=1
do k1=1,total_n
if(I_temp(j).eq.Neighbor_V(k1))i5=0

enddo
if(I_temp(j).eq.V1.or.I_temp(j).eq.V2)then
i4=mod(i4+1,2);i5=0
endif
if(i5.eq.1)then
if(i4.eq.0)then
Id_v3=Id_v3+1;
vdi(Id_v1+Id_v2+Id_v3)=I_temp(j)
endif
if(i4.eq.1)then
Iu_v2=Iu_v2+1;
vui(Iu_v1+Iu_v2+Iu_v3)=I_temp(j)
endif
endif
if(I_temp(j).ne.V1.and.I_temp(j).ne.V2)then
if(i4.eq.0)then
Id_v3e=Id_v3e+1
Edi(Id_v1+Id_v2e+Id_v3e)=I_temp2(j)
Fdi(Id_v1+Id_v2e+Id_v3e)=I_temp3(j)
endif
if(i4.eq.1)then
Iu_v3e=Iu_v3e+1
Eui(Iu_v1+Iu_v2e+Iu_v3e)=I_temp2(j)
Fui(Iu_v1+Iu_v2e+Iu_v3e)=I_temp3(j)
endif
endif

enddo
V_del(1)=v1; V_del(2)=v2; V_del(3)=v3
E_del(1:Id_v1+Id_v2e+Id_v3e)=Edi(1:Id_v1+Id_v2e+Id_v3e)
E_del(Id_v1+Id_v2e+Id_v3e+1:Iu_v1+Iu_v2e+Iu_v3e+Id_v1+Id_v2e+Id_v3e)=Eui(1:Iu_v1+Iu_v2e+Iu_v3e)
F_del(1:Id_v1+Id_v2e+Id_v3e)=Fdi(1:Id_v1+Id_v2e+Id_v3e)
F_del(Id_v1+Id_v2e+Id_v3e+1:Iu_v1+Iu_v2e+Iu_v3e+Id_v1+Id_v2e+Id_v3e)=Fui(1:Iu_v1+Iu_v2e+Iu_v3e)
ke=Iu_v1+Iu_v2e+Iu_v3e+Id_v1+Id_v2e+Id_v3e+3
kd=Id_v1+Id_v2+Id_v3
ku=Iu_v1+Iu_v2+Iu_v3

! there are still three edges that should be deleted. They are v1v2,v2v3 and v3v1
N_t=N_V_V(V1)
I_temp(1:N_t)=abs(E_V(v1,1:N_t))
do i=1,N_t
i1=(V_E(I_temp(i),1))
i2=(V_E(I_temp(i),2))
if((i1.eq.v1.and.i2.eq.v2).or.(i1.eq.v2.and.i2.eq.v1))then
E_del(ke-2)=I_temp(i)

endif
enddo
N_t=N_V_V(V2)
I_temp(1:N_t)=abs(E_V(v2,1:N_t))
do i=1,N_t
i1=(V_E(I_temp(i),1))
i2=(V_E(I_temp(i),2))
if((i1.eq.v2.and.i2.eq.v3).or.(i1.eq.v3.and.i2.eq.v2))E_del(ke-1)=I_temp(i)
enddo
N_t=N_V_V(V3)
I_temp(1:N_t)=abs(E_V(v3,1:N_t))
do i=1,N_t
i1=(V_E(I_temp(i),1))
i2=(V_E(I_temp(i),2))
if((i1.eq.v3.and.i2.eq.v1).or.(i1.eq.v1.and.i2.eq.v3))E_del(ke)=I_temp(i)
enddo

! %%%%%%%%%%%%%%%%%%%%%%% end of finding vd vu, E, F
! sorting
k=3
call sort_num(V_del,k)
call sort_num(E_del,ke)
k=ke-3
call sort_num(F_del,k)

Ic_V=0
Irc_V=0
do i=1,N_v
do j=1,3+1
if(j.eq.1)then
if(i.lt.V_del(j))then
Ic_V(i)=i
Irc_V(i)=i
endif
endif
if(j.eq.3+1)then
if(i.gt.V_del(j-1))then
Ic_V(i)=i-3
Irc_V(i-3)=i
endif
endif
if(j.gt.1.and.j.lt.3+1)then
if(i.gt.V_del(j-1).and.i.lt.V_del(j))then
Ic_V(i)=i-j+1
Irc_V(i-j+1)=i
endif
endif
enddo
enddo
N_vc=N_v-1



do i=1,N_E
do j=1,ke+1
if(j.eq.1)then
if(i.lt.E_del(j))then
Ic_E(i)=i
Irc_E(i)=i
endif
endif
if(j.eq.ke+1)then
if(i.gt.E_del(j-1))then
Ic_E(i)=i-ke
Irc_E(i-ke)=i
endif
endif
if(j.gt.1.and.j.lt.ke+1)then
if(i.gt.E_del(j-1).and.i.lt.E_del(j))then
Ic_E(i)=i-j+1
Irc_E(i-j+1)=i
endif
endif
enddo
enddo
N_Ec=N_E-ke+kd+ku


do i=1,N_F
do j=1,ke-2
if(j.eq.1)then
if(i.lt.F_del(j))then
Ic_F(i)=i
Irc_F(i)=i
endif
endif
if(j.eq.ke-2)then
if(i.gt.F_del(j-1))then
Ic_F(i)=i-ke+3
Irc_F(i-ke+3)=i
endif
endif
if(j.gt.1.and.j.lt.ke-2)then
if(i.gt.F_del(j-1).and.i.lt.F_del(j))then
Ic_F(i)=i-j+1
Irc_F(i-j+1)=i
endif
endif
enddo
enddo
N_Fc=N_F-ke+3+kd+ku

!  end of Ic_VEF
do i=1,kd
Vdi_new(i)=Ic_V(Vdi(i))
enddo
do i=1,ku
Vui_new(i)=Ic_V(Vui(i))
enddo
vd=N_V-2; vu=N_v-1

! new V_Ec (except for those new adding edges or faces adjacent with vertices vd and vu)
do i=1,N_Ec-kd-ku
V_Ec(i,1)=Ic_V(V_E(Irc_E(i),1))
V_Ec(i,2)=Ic_V(V_E(Irc_E(i),2))
enddo


do i=1,N_Fc-kd-ku
i1=Ic_E(abs(E_F(Irc_F(i),1)))
i2=Ic_E(abs(E_F(Irc_F(i),2)))
i3=Ic_E(abs(E_F(Irc_F(i),3)))
E_Fc(i,1)=i1*E_F(Irc_F(i),1)/abs(E_F(Irc_F(i),1))
E_Fc(i,2)=i2*E_F(Irc_F(i),2)/abs(E_F(Irc_F(i),2))
E_Fc(i,3)=i3*E_F(Irc_F(i),3)/abs(E_F(Irc_F(i),3))
enddo


!%%%%%%%%%%%%%%
! add new edges and new faces

I_Edop=0

do i=1,Id_v1+Id_v2e+Id_v3e
do j=1,3
k=E_F(Fdi(i),j)
i1=V_E(abs(k),1)
i2=V_E(abs(k),2)
i3=1
i3=(i1-v1)*(i1-v2)*(i1-v3)
i3=i3*(i2-v1)*(i2-v2)*(i2-v3)
if(i3.ne.0)then
I_Edop=I_Edop+1
Ed_op(I_Edop)=k
endif
enddo
enddo

I_Euop=0
do i=1,Iu_v1+Iu_v2e+Iu_v3e
do j=1,3
k=E_F(Fui(i),j)
i1=V_E(abs(k),1)
i2=V_E(abs(k),2)
i3=1
i3=(i1-v1)*(i1-v2)*(i1-v3)
i3=i3*(i2-v1)*(i2-v2)*(i2-v3)
if(i3.ne.0)then
I_Euop=I_Euop+1
Eu_op(I_Euop)=k
endif
enddo
enddo

! %%%%%%%%%%%%%%

i1=Ed_op(1)
if(i1.gt.0)then
i2=V_E(abs(i1),1)
i3=V_E(abs(i1),2)
else
i2=V_E(abs(i1),2)
i3=V_E(abs(i1),1)
endif
Order_Ed(1)=1
Vdi_order(1)=i2
Vdi_order(2)=i3
k=1
do while(k.le.I_Edop)
Loop1: do i=1,I_Edop
i1=Ed_op(i)
if(i1.gt.0)then
k2=V_E(abs(i1),1)
k3=V_E(abs(i1),2)
else
k2=V_E(abs(i1),2)
k3=V_E(abs(i1),1)
endif
if(k2.eq.i3)then
k=k+1
order_Ed(k)=i
Vdi_order(k+1)=k3
i2=k2
i3=k3
exit Loop1
endif
enddo Loop1

enddo


i1=Eu_op(1)
if(i1.gt.0)then
i2=V_E(abs(i1),1)
i3=V_E(abs(i1),2)
else
i2=V_E(abs(i1),2)
i3=V_E(abs(i1),1)
endif
Order_Eu(1)=1
Vui_order(1)=i2
Vui_order(2)=i3
k=1
do while(k.le.I_Euop)
Loop2: do i=1,I_Euop
i1=Eu_op(i)
if(i1.gt.0)then
k2=V_E(abs(i1),1)
k3=V_E(abs(i1),2)
else
k2=V_E(abs(i1),2)
k3=V_E(abs(i1),1)
endif
if(k2.eq.i3)then
k=k+1
order_Eu(k)=i
Vui_order(k+1)=k3
i2=k2
i3=k3
exit Loop2
endif
enddo Loop2

enddo

!  construct new vertices new edges faces
do i=1,kd
Vdi_new(i)=Ic_V(Vdi_order(i))
enddo
do i=1,ku
Vui_new(i)=Ic_V(Vui_order(i))
enddo
! new adding V_E
do i=1,kd
V_Ec(N_Ec-kd-ku+i,1)=Vd
V_Ec(N_Ec-kd-ku+i,2)=Vdi_new(i)
enddo
do i=1,ku
V_Ec(N_Ec-ku+i,1)=Vu
V_Ec(N_Ec-ku+i,2)=Vui_new(i)
enddo
! new E_F
do i=1,kd
E_Fc(N_Fc-kd-ku+i,1)=N_Ec-kd-ku+i
i1=Ed_op(order_Ed(i))
k1=Ic_E(abs(i1))
k2=k1*i1/abs(i1)
E_Fc(N_Fc-kd-ku+i,2)=k2
E_Fc(N_Fc-kd-ku+i,3)=-(N_Ec-kd-ku+i+1)
if(i.eq.kd)E_Fc(N_Fc-kd-ku+i,3)=-(N_Ec-kd-ku+1)
enddo

do i=1,ku
E_Fc(N_Fc-ku+i,1)=N_Ec-ku+i
i1=Eu_op(order_Eu(i))
k1=Ic_E(abs(i1))
k2=k1*i1/abs(i1)
E_Fc(N_Fc-ku+i,2)=k2
E_Fc(N_Fc-ku+i,3)=-(N_Ec-ku+i+1)
if(i.eq.ku)E_Fc(N_Fc-ku+i,3)=-(N_Ec-ku+1)
enddo

rv0=rv
fi0=fi
do i=1,N_V-3
rv(i,1:3)=rv0(Irc_V(i),1:3)
fi(i)=fi0(Irc_V(i))
enddo
rv(N_V-2,1:3)=(rv0(v1,1:3)+rv0(v2,1:3)+rv0(v3,1:3))/3.0
rv(N_V-1,1:3)=(rv0(v1,1:3)+rv0(v2,1:3)+rv0(v3,1:3))/3.0
fi(N_V-2)=(fi0(v1)+fi0(v2)+fi0(i3))/3.0
fi(N_V-1)=(fi0(v1)+fi0(v2)+fi0(i3))/3.0

V_E=V_Ec
E_F=E_Fc
N_V=N_Vc
N_E=N_Ec
N_F=N_Fc

call label_system(V_F,   E_V,N_E_V,     F_V,N_F_V,    F_E,N_F_E,    V_V,N_V_V,   E_V_opposite,N_E_V_opposite,   V_E_opposite,  N_V,N_E,N_F,V_E,E_F,F_F,rv)
write(*,*)'Fission has been taken place'

END











subroutine sort_num(data0, L)
integer L
integer data1(1:40)
integer data0(1:40)
integer i,j,k
integer order_data(1:40)
integer a
data1=abs(data0)
do i=1,L
a=data1(i)
k=1
do j=1,L
if(i.ne.j)then
if(a.gt.data1(j))k=k+1
endif
enddo
data0(k)=a
enddo

END

