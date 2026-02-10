!  ************************************************
!  Part of 3D vesicle program
!       phase separation
!                    written by LiJF  Dobb@bbs.fudan.edu.cn
! *************************************************  




subroutine initial_phase(fia0, fi,N_V, area_V,order_fi,seed1)

parameter(N = 20000)

real*8 fi(0:N),area_V(0:N),fia0
real*8 ave_V(0:N), sum1,sum2,sum3,sum4
real*8 fia(0:N),order_fi,fib(0:N)
real*8 area_tot,random1

integer i,j,k, N_V
integer seed1


ave_V=area_V/3.0
area_tot=sum( area_V(0:N_V) )/3.0
sum1=0.0
sum2=0.0
    do i=1, N_v

	  fia(i) = fia0+(0.5-random1(seed1))*0.5
	  fib(i) = 1.0-fia(i)
	  sum1   = sum1+fia(i)*ave_V(i)/area_tot;sum2   = sum2+fib(i)*ave_V(i)/area_tot
    end do
	sum3 = 0.0;sum4 = 0.0

    do i=1, N_V
	  fia(i) = fia(i)-sum1+fia0;fib(i) = fib(i)-sum2+1.0-fia0
	  fi(i)  = fia(i)-fib(i)
	  sum3   = sum3+fia(i)*ave_V(i)/area_tot; sum4   = sum4+fib(i)*ave_V(i)/area_tot

	enddo
      order_fi=sum3-sum4


END


subroutine update_fi(fi,order_fi,L_fi,dt, H_modulus, H_V,H_zero, zet_Hfi, b, a2,a4,area_F,area_V, N_V,N_F, N_F_V,V_V, F_V,V_F, L_glb,A_phase,  rv)

parameter(N = 20000)

real*8 rv(0:N,1:3)
real*8 fi(0:N),fi0(0:N)
real*8 L_glb(0:N,1:20)
real*8 H_V(0:N),H_zero, zet_Hfi, H_modulus

integer N_V, N_F_V(0:N), N_t, V_V(0:N,1:20),F_V(0:N,1:20)
integer V_F(0:N,1:3),N_F,N_E

real*8 area_F(0:N),area_V(0:N), area_V1(0:N)
real*8 vector_F(0:N,1:3),norm_F(0:N,1:3)
real*8 b, a2,a4,area_tot,  sum1, order_fi
real*8 A_phase

real*8 Delt_fi(0:N), r_t

real*8 L_fi, dt

integer i,j, i1


delt_fi=0.0
do i=1,N_V
N_t=N_F_V(i)
r_t=0.0
r_t=r_t -a2*fi(i)+a4*fi(i)**3
r_t=r_t + H_modulus*(H_v(i)-H_zero-zet_Hfi*fi(i))*(-zet_Hfi)

do j=1,N_t
i1=V_V(i,j)
r_t=r_t-b * (fi(i1)-fi(i))*L_glb(i,j)
enddo
Delt_fi(i)=r_t
enddo


 call Get_area_F( rv, V_F,  N_F,N_V,   Area_F, Vector_F, Norm_F)

do i=1, N_V
area_V1(i)=0.0
do j=1,N_F_V(i)
area_V1(i)=area_V1(i)+area_F(F_V(i,j))
enddo
enddo

 call Laplace_global(rv, area_F,area_V1, N_V, V_V, F_V, N_F_V, L_glb)
 !area_V1=area_V
do i=1,N_V
r_t=0.0
N_t=N_F_V(i)
do j=1,N_t
i1=V_V(i,j)
r_t=r_t+(Delt_fi(i1)-Delt_fi(i))*L_glb(i,j)
enddo
fi(i)=(fi(i)*(1- (area_V1(i)-area_V(i))/(area_V1(i)+area_V(i)) )+L_fi*dt*r_t  )/(1+ (area_V1(i)-area_V(i))/(area_V1(i)+area_V(i))  )
! if(fi(i).gt.1.0)fi(i)=1.0
! if(fi(i).lt.-1.0)fi(i)=-1.0

!! this scheme can be reduced to be Taniguchi's form by setting d S=0
enddo






	sum1=0.0
area_tot=sum(area_V1(0:N_V))/3.0

	do i=1, N_v
	  sum1   =   sum1+fi(i)*area_v1(i)/3.0/area_tot
	enddo
	
    do i=1, N_v

	  fi(i) =  fi(i)-(sum1-order_fi)
	
	enddo


END


subroutine convert_FNUM(fi, I_fi_F, fi_F, V_F, N_F)

parameter(N = 20000)
real*8 fi(0:N),FI_F(0:N)
integer I_fi_F(0:N), V_F(0:N,1:3), N_F

integer i,j,k,i1,i2,i3

do i=1,N_F
Fi_F(i)= (Fi(V_F(i,1))+Fi(V_F(i,2))+Fi(V_F(i,3)))/3.0

if(fi_F(i).ge.0.7)I_FI_F(I)=0  ! black
if(fi_F(i).ge.0.0.and.fi_F(i).lt.0.7)I_FI_F(I)=8  ! Darkgrey
!if(fi_F(i).gt.0.0.and.fi_F(i).lt.0.4)I_FI_F(I)=2  ! blue
!if(fi_F(i).eq.0.0)I_FI_F(I)=4
!if(fi_F(i).lt.0.0.and.fi_F(i).gt.-0.4)I_FI_F(I)=3  ! green
if(fi_F(i).le.-0.0.and.fi_F(i).gt.-0.7)I_FI_F(I)=7  ! lightgrey
if(fi_F(i).le.-0.7)I_FI_F(I)=15  ! white

do j=0,15
if(fi_F(i).ge.(-1.0+j*2.0/16.0).and.fi_F(i).lt.(-1.0+(j+1)*2.0/16.0) )I_FI_F(i)=15-j

enddo
if(fi_F(i).lt.-1.0)I_FI_F(i)=15
if(fi_F(i).gt.1.0)I_FI_F(i)=0


enddo






END


subroutine boundary_count(rv,fi,Len_boundary,fia0,N_E,F_E,V_E,V_E_opposite, switch_boundary_setting,boundary_position, Num_boundary,Num_domain,Circle_F,len_domain,Point_circle,all_F_domain)


parameter(N = 20000)

real*8 fi(0:N),fia0
real*8 fi_cut
real*8 rv(0:N,1:3)
real*8 len_boundary
real*8 len_domain(1:1000),len_all(1:N)
real*8 boundary_position(0:N,1:6)
real*8 len_rr
integer i,j,k
integer N_E,N_V
integer V_E(0:N,1:2),V_E_opposite(0:N,1:20)
integer F_E(0:N,1:20), Bound_vect(0:N,1:2)
integer Circle_F(1:100,1:50000),point_circle(1:1000),Num_circle
integer Counted_points(0:N),Counted_Num,switch_close
integer NotCounts(0:N),k_countF,NotCounts0(0:N)
integer F_start0,F_end0,F_start,F_end,F1,F2
integer  F_start1,F_end1

integer all_F_domain(1:N)

integer Num_boundary,switch_boundary_setting
integer Num_domain
real*8 R1(1:3),R2(1:3)
real*8 R3(1:3),R4(1:3)
real*8 RR(1:6)

integer E_b_all(1:N), E_b(1:1000,1:1000)


integer i1,i2,i3,i4,i5,i6

if(switch_boundary_setting.eq.0)then
fi_cut=0.0
else
fi_cut=2*fia0-1.0
endif

Num_boundary=0
len_boundary=0.0
all_F_domain=0
bound_vect=0

do i=1,N_E
i1=V_E(i,1);i2=V_E(i,2);


if((fi(i1)-fi_cut)*(fi(i2)-fi_cut).lt.0.0)then
Num_boundary=Num_boundary+1
E_b_all(Num_boundary)=i
i3=V_E_opposite(i,1)
i4=V_E_opposite(i,2)
i5=F_E(i,1)
i6=F_E(i,2)
bound_Vect(Num_boundary,1)=i5;bound_Vect(Num_boundary,2)=i6
R1=rv(i1,1:3)
R2=rv(i2,1:3)
R3=rv(i3,1:3)
R4=rv(i4,1:3)
RR(1:3)=(R1+R2+R3)/3.0
RR(4:6)=(R1+R2+R4)/3.0
boundary_position(Num_boundary,1:6)=RR
len_all(Num_boundary)=(sum((RR(1:3)-RR(4:6))**2.0))**0.5
all_F_domain(Num_boundary)=i5
all_F_domain(Num_boundary+1000)=i6
len_boundary=len_boundary+len_all(Num_boundary)
endif

enddo


if(Num_boundary.eq.0)return

! %%%%%%%%%%%%%%%%%%%
! %%%%%%%%%%%%%%%%%%         domain numbers

NOtcounts=0
do i=1,Num_boundary
Notcounts(i)=i
enddo
 counted_num=0
F_start0=bound_Vect(1,1)
 Circle_F(1,1)=F_start0
F_end1=F_start0
Num_circle=1
Point_circle=0
k_countF=0
len_domain=0
!len_domain(1)=len_all(1)
E_b(1,1)=E_b_all(1)
do while(counted_Num.lt.Num_boundary+1)

if(switch_close.eq.1)then  ! a domain is closed,
Num_circle=Num_circle+1
F_start0=bound_Vect(notcounts(1),1)
 Circle_F(Num_Circle, 1)=F_start0
 E_b(Num_Circle,1)=E_b_all(notcounts(1))
F_end1=F_start0
switch_close=0
endif




loop2:do i=1,Num_boundary-Counted_Num
F1=bound_Vect(Notcounts(i),1)
F2=bound_Vect(Notcounts(i),2)
if(F1.eq.F_end1)then
k_counF=NotCounts(i)
k=i
F_start=F1
F_end=F2
if(F_end.eq.F_start0)then
switch_close=1
endif
F_start1=F_start
F_end1=F_end
point_circle(Num_circle)=point_circle(Num_circle)+1
 Circle_F(Num_Circle,Point_circle(Num_circle))=F_start
len_domain(Num_Circle)=len_Domain(Num_Circle)+len_all(Notcounts(i))
E_b(Num_circle,Point_circle(Num_circle))=E_b_all(Notcounts(i))

exit loop2
endif
if(F2.eq.F_end1)then
k_counF=NotCounts(i)
k=i
F_start=F2
F_end=F1
if(F_end.eq.F_start0)then
switch_close=1
endif
F_start1=F_start
F_end1=F_end
point_circle(Num_circle)=point_circle(Num_circle)+1
 Circle_F(Num_Circle,Point_circle(Num_circle))=F_start
len_domain(Num_Circle)=len_Domain(Num_Circle)+len_all(Notcounts(i))
E_b(Num_circle,Point_circle(Num_circle))=E_b_all(Notcounts(i))
exit loop2
endif



enddo loop2

Notcounts0=Notcounts
Notcounts(k:Num_boundary-Counted_Num-1)=Notcounts0(k+1:Num_boundary-Counted_Num)
 counted_points(counted_Num+1)=k_countF
 counted_Num=counted_Num+1

if(counted_Num.lt.2)write(*,*)

enddo ! END of while


Num_circle=Num_Circle-1
Num_domain=Num_circle






END




subroutine domain_area(area_F,fi_F,F_F,N_F,Num_domain,Circle_F,len_domain,Point_circle,area_domain,switch_boundary_setting,fia0,  area_domain0, Num_domain0,F_domain0,N_F_domain0,No_domain0,No_domain_Q0,inv_No_domain0,Max_No0,ratio_area_len)

parameter(N = 20000)
parameter(PI = 3.1415926536)
real*8 fia0
real*8 fi_F(0:N)
real*8 fi_cut


integer F_F(0:N,1:3),N_F


! old inputs
real*8 area_domain0(1:1000)
integer Num_domain0, F_domain0(1:100,1:50000),N_F_domain0(1:1000)
integer No_domain0(1:1000,1:100), No_domain_Q0(1:1000)
  ! i(1:1000) for 1:Num_domain0,   return Value for actual domain ID 
  ! k(1:100) for overload domain ID, No_domain_Q0 is for overload numbers
integer inv_No_domain0(1:1000)
integer Max_No0  ! number of output files

! New outputs
integer No_domain(1:1000,1:100), No_domain_Q(1:1000)
  ! i(1:1000) for 1:Num_domain0,   return Value for actual domain ID 
  ! k(1:100) for overload domain ID, No_domain_Q0 is for overload numbers
integer inv_No_domain(1:1000)
integer Max_No

! %%%%%%%%%%%%%  variables for name the domains
integer New_old(1:1000),old_New(1:1000)
integer order_Num_domain0(1:1000), order_Num_domain(1:1000)
real*8 temp_real1(1:1000),temp_real2(1:1000)
integer temp_int1(1:1000),temp_int2(1:1000)
integer temp_int3(1:1000),temp_int4(1:1000)

! %%%%%%%%%%%%%%%%%%%



integer Circle_F(1:100,1:50000),point_circle(1:1000),Num_domain
real*8 len_domain(1:1000)

integer F_domain(1:100,1:50000),N_F_domain(1:1000)
real*8 area_domain(1:1000)
real*8 area_F(0:N), ratio_area_len(1:1000)

integer F_in_domain( 1 :N),N_F_in_domain
integer counted_F(1:N),N_counted_F
integer Num_domain_F(1:N)
integer J_F_along_road(1:N),F_along_road(1:N),k_F_along_road
integer F_series(1:1000),N_series(1:1000),k_N_series,Find_it
integer Is_domain_boundary(1:N)

integer i,j,k,kk,i1,i2,i3,j1,j2,j3, type_k,switch_boundary_setting


if(switch_boundary_setting.eq.0)then
fi_cut=0.0
else
fi_cut=2*fia0-1.0
endif

! %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% First time to count domain area

! let domain boundaries be F_domain
N_F_domain=point_circle
F_domain=circle_F
N_counted_F=0
Num_domain_F=0  ! to specific which domain the facet belongs to
Is_domain_boundary=0
do i=1,Num_domain
do j=1,point_circle(i)
N_counted_F=N_counted_F+1
 counted_F(N_counted_F)=Circle_F(i,j)
Num_domain_F(Circle_F(i,j))=i
Is_domain_boundary(Circle_F(i,j))=1
enddo
enddo



! Get all the facets within domains

N_F_in_domain=0
do i=1,N_F
if(fi_F(i).ge.fi_cut.and.Is_domain_boundary(i).eq.0)then
N_F_in_domain=N_F_in_domain+1
F_in_domain(N_F_in_domain)=i
endif
enddo

!  appoint the domain No. to each facet
do i=1, N_F_in_domain
k=F_in_domain(i)
J_F_along_road=0       ! to Judge whether the facet is along the road 
J_F_along_road(k)=1
F_along_road(1)=k
F_series(1)=k       !  
Find_it=0; k_N_series=1; k_F_along_road=1
do while(Find_it.eq.0.and.Num_domain_F(k).eq.0)
kk=F_series(k_N_series)
i1=F_F(kk,1);i2=F_F(kk,2);i3=F_F(kk,3)

if(Num_domain_F(i1).ne.0)then
type_k=Num_domain_F(i1)
Find_it=1
endif
if(Num_domain_F(i2).ne.0)then
type_k=Num_domain_F(i2)
Find_it=1
endif
if(Num_domain_F(i3).ne.0)then
type_k=Num_domain_F(i3)
Find_it=1
endif
i4=1
if(J_F_along_road(i1).eq.0.and.Find_it.eq.0.and.i4.eq.1)then
k_F_along_road=k_F_along_road+1;
J_F_along_road(i1)=1; F_along_road(k_F_along_road)=i1
k_N_series=k_N_series+1
N_series(k_N_series)=1
F_series(k_N_series)=i1
i4=0
endif
if(J_F_along_road(i1).ne.0.and.J_F_along_road(i2).eq.0.and.Find_it.eq.0.and.i4.eq.1)then
k_F_along_road=k_F_along_road+1;
J_F_along_road(i2)=1; F_along_road(k_F_along_road)=i2
k_N_series=k_N_series+1
N_series(k_N_series)=2
F_series(k_N_series)=i2
i4=0
endif
if(J_F_along_road(i1).ne.0.and.J_F_along_road(i2).ne.0.and.J_F_along_road(i3).eq.0.and.Find_it.eq.0.and.i4.eq.1)then
k_F_along_road=k_F_along_road+1;
J_F_along_road(i3)=1; F_along_road(k_F_along_road)=i3
k_N_series=k_N_series+1
N_series(k_N_series)=3
F_series(k_N_series)=i3
i4=0
endif


if(J_F_along_road(i1).ne.0.and.J_F_along_road(i2).ne.0.and.J_F_along_road(i3).ne.0.and.Find_it.eq.0.and.i4.eq.1)then
k_N_series=k_N_series-1
endif

enddo ! end while

do j=1,K_F_along_road
if(Num_domain_F(F_along_road(j)).eq.0)then
Num_domain_F(F_along_road(j))=type_k
endif
enddo

type_k=0

enddo ! end of i




do i=1, N_F
 kk=Num_domain_F(i)

if(kk.ne.0.and.Is_domain_boundary(i).eq.0)then
N_F_domain(kk)=N_F_domain(kk)+1
F_domain(kk,N_F_domain(kk))=i
endif


enddo

! counting ereas
area_domain=0.0
do i=1,Num_domain
do j=1, N_F_domain(i)
area_domain(i)=area_domain(i)+area_F(F_domain(i,j))
enddo
ratio_area_len(i)=area_domain(i)*4*PI/len_domain(i)**2
enddo

! End of Counting 


! %%%%%%%%%%%%%%%%%%%%%%%
! name the domains
! %%%%%%%%%%%%%%%%%%%%%%%
! %%%% the first time to name the domains

if(Num_domain0.eq.0)then
No_domain_Q=0
inv_No_domain_Q=0
No_domain_Q(1:Num_domain)=1
No_domain=0; inv_No_domain=0
do i=1,Num_domain
No_domain(i,1)=i
inv_No_domain(i)=i
enddo
Max_No=Num_domain
endif ! END of NOT the first time to name the domains 



! %%% NOT the first time to name the domains
if(Num_domain0.ne.0)then
! %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! 1st group, domain number is static
if(Num_domain.eq.Num_domain0)then
Max_No=Max_No0
temp_real1(1:Num_domain) =area_domain(1:Num_domain)
temp_real2(1:Num_domain0)=area_domain0(1:Num_domain0)

 call sort_num2(temp_real1,Num_domain,temp_int1,temp_int2) ! new
 call sort_num2(temp_real2,Num_domain0,temp_int3,temp_int4) !old


do i=1, Num_domain0

!  %%%%%%%%%%%          i2 -> i2
New_old(i)= temp_int2(temp_int3(i))
i1=i; i2=New_old(i); 
i3=N_F_domain0(i1);i4=N_F_domain(i2)
i5=0

loopj:do j=1,i4
do k=1,i3
if(F_domain0(i1,k).eq.F_domain(i2,j))then
i5=1

exit loopj
endif
enddo
enddo loopj

!  %%%%%%%%%%%          i2 -> i2+1
if(i5.eq.0.and.temp_int3(i).lt.Num_domain)then
New_old(i)= temp_int2(temp_int3(i)+1)
i1=i; i2=New_old(i);
i3=N_F_domain0(i1);i4=N_F_domain(i2)

loopj1:do j=1,i4
do k=1,i3
if(F_domain0(i1,k).eq.F_domain(i2,j))then
i5=1
exit loopj1
endif
enddo
enddo loopj1
endif

!  %%%%%%%%%%%          i2 -> i2-1
if(i5.eq.0.and.temp_int3(i).gt.1)then
New_old(i)= temp_int2(temp_int3(i)-1)
i1=i; i2=New_old(i); 
i3=N_F_domain0(i1);i4=N_F_domain(i2)

loopj2:do j=1,i4
do k=1,i3
if(F_domain0(i1,k).eq.F_domain(i2,j))then
i5=1
exit loopj2
endif
enddo
enddo loopj2
endif


! %%%%%%%%  have to judge over betweeen all domains
if(i5.eq.0)then
i1=i;
i3=N_F_domain0(i1)
loopj3: do kk=1,Num_domain
i4=N_F_domain(kk)
do j=1,i4
do k=1,i3
if(F_domain0(i1,k).eq.F_domain(kk,j))then
i5=1
New_old(i)=kk


exit loopj3

endif
enddo
enddo 
enddo loopj3

endif


enddo ! end of i Num_domain

! %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
do i=1, Num_domain0
k=New_old(i)
old_new(k)=i
enddo

do i=1,Num_domain
No_domain_Q(i)=No_domain_Q0(old_new(i))
No_domain(i,1:No_domain_Q(i))=No_domain0(old_new(i),1:No_domain_Q(i))
do j=1,No_domain_Q(i)
inv_No_domain(No_domain(i,j))=i
enddo
enddo


endif ! end of 1st group



!  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! 2nd group,  Num_domain = Num_domain0-1  one domain dispeared or two domains merge into one domain
! because this does not occur frequently, so there is no need to design a very good codes...
if(Num_domain.le.Num_domain0-1)then
Max_No=Max_No0
do i=1,Num_domain0
i5=0
Loop2nd1:do j=1,Num_domain
do k=1,N_F_domain0(i)
do kk=1, N_F_domain(j)
if(F_domain0(i,k).eq.F_domain(j,kk))then
new_old(i)=j
i5=1
exit Loop2nd1
endif
enddo
enddo
enddo Loop2nd1 ! end of enumurate new Num_domain

if(i5.eq.0)then
New_old(i)=0 ! this means domain i has been dissapear
endif

enddo ! end of enumurate old Num_domain0

No_domain_Q=0;temp_int1=1
do i=1,Num_domain0
k=New_old(i)
if(k.ne.0)then
No_domain_Q(k)=No_domain_Q(k)+No_domain_Q0(i) ! to handle merging situation
No_domain(k,temp_int1(k):No_domain_Q(k))=No_domain0(i,1:No_domain_Q0(i))

do j=temp_int1(k),No_domain_Q(k)
inv_No_domain(No_domain(k,j))=k
enddo
temp_int1(k)=No_domain_Q(k)+1
endif
if(k.eq.0)then
do j=1,No_domain_Q0(i)
inv_No_domain(No_domain0(i,j))=0   ! domain disppear
enddo
endif
enddo

endif ! 2nd group,  Num_domain != Num_domain0-1  one domain dispeared!!


! %%%%%%%%%%%%%%%%%%%%%%%%%%
! 3rd group Num_domain=No_domain0+1 
!         rarely happen! a new domain appear
if(Num_domain.ge.Num_domain0+1)then
Max_No=Max_No0
do i=1,Num_domain
i5=0
Loop3nd1:do j=1,Num_domain0
do k=1,N_F_domain(i)
do kk=1, N_F_domain0(j)
if(F_domain(i,k).eq.F_domain0(j,kk))then
old_new(i)=j
i5=1
exit Loop3nd1
endif
enddo
enddo
enddo Loop3nd1 ! end of enumurate new Num_domain

if(i5.eq.0)then
old_new(i)=0 ! this means domain i has been dissapear
endif

enddo ! end of enumurate old Num_domain

No_domain_Q=0;
do i=1,Num_domain
k=old_new(i)
if(k.ne.0)then
No_domain_Q(i)=No_domain_Q0(k) ! to handle merging situation
No_domain(i,1:No_domain_Q(i))=No_domain0(k,1:No_domain_Q0(k))

do j=1,No_domain_Q(i)
inv_No_domain(No_domain(i,j))=i
enddo

endif
if(k.eq.0)then
!Max_No=Max_No0+1
!No_domain(i,1)=Max_No
!No_domain_Q(i)=1
!inv_No_domain(Max_No)=i   ! a new domain appear
endif
enddo



endif ! 3rd group







! %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
endif ! END of NOT the first time to name the domains 



! %%%%%%%%%%%%%%%%%%%%%%% update the old inputs
F_domain0=F_domain
N_F_domain0=N_F_domain
Num_domain0=Num_domain

No_domain_Q0=No_domain_Q
No_domain0=No_domain
inv_No_domain_Q0=inv_No_domain_Q
inv_No_domain0=inv_No_domain
Max_No0=Max_No
area_domain0=area_domain




END



























