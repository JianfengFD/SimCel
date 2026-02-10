!  ************************************************
!  Part of 3D vesicle program
!   calculate energies
!                    written by LiJF  Dobb@bbs.fudan.edu.cn
! ************************************************* 



! PV
real*8 function PV_energy(P,Vol_total)
real*8 P, Vol_total
PV_energy=P*Vol_total
END

! \sum_{i_V} H_modulus / 2 (H_v-H_0-zet_Hfi * fi)^2 Area_V/3

real*8 function Curvature_Energy(H_modulus,H_v, H_zero, zet_Hfi, fi, Area_V,N_V)
parameter(N = 20000)
integer N_V
real*8 H_modulus, H_zero, zet_Hfi
real*8  fi(0:N), Area_V(0:N),H_v(0:N)
real*8 a
integer i

a=0.0
do i=1,N_V
a=a+ H_modulus / 2.0 * (H_v(i)-H_zero-zet_Hfi * fi(i) )**2* Area_V(i)/3.0
enddo
Curvature_Energy=a
END

! \sum_{i_F} lamda/2.0*(area_F-area_F_0)**2   :: Hooke energy to ensure local area incompressible

real*8 function Area_Hooke_energy(lamda, area_F,area_F_0,N_F)
parameter(N = 20000)
integer N_F,i
real*8 lamda,area_F(0:N),area_F_0(0:N),a
a=0.0
do i=1,N_F
a=a+lamda/2.0*( area_F(i)-area_F_0(i) )**2

enddo
Area_Hooke_energy=a
END


!phase

real*8 function phase_energy(b, a2, a4, fi, Area_V, Area_F,N_V,N_F,V_F, E_F,rv, A_phase)
parameter(N = 20000)
real*8 b,a2,a4,A_phase
real*8 fi(0:N), area_v(0:N),area_F(0:N),vector_E(0:N,1:3),len_E(0:N)
real*8 rv(0:N,1:3)

integer N_V,N_F_V(0:N), E_F(0:N,1:3), V_F(0:N,1:3)
real*8 R(1:3,1:3), tr(1:3,1:3), fit(1:3),dfi(1:3)

integer i,j,k,I_t(1:3)
real*8 E1,E2

E1=0.0
E2=0.0
do i=1,N_F
i_t=V_F(i,1:3)
fit(1:3)=fi(i_t(1:3))
R(1,1:3)=rv(i_t(1),1:3)
R(2,1:3)=rv(i_t(2),1:3)
R(3,1:3)=rv(i_t(3),1:3)
tr(1,1:3)=R(2,1:3)-R(1,1:3)
tr(2,1:3)=R(3,1:3)-R(2,1:3)
tr(3,1:3)=R(1,1:3)-R(3,1:3)
dfi(1)=fit(2)-fit(1)
dfi(2)=fit(3)-fit(2)
dfi(3)=fit(1)-fit(3)

E1=E1+b/24.0*( dfi(1)**2* (sum(tr(2,1:3)**2)+ sum(tr(3,1:3)**2)) - 2*dfi(1)*dfi(2)* sum(tr(1,1:3)*tr(2,1:3) )        )/area_F(i)
E1=E1+b/24.0*( dfi(2)**2* (sum(tr(1,1:3)**2)+ sum(tr(3,1:3)**2)) - 2*dfi(2)*dfi(3)* sum(tr(2,1:3)*tr(3,1:3) )        )/area_F(i)
E1=E1+b/24.0*( dfi(3)**2* (sum(tr(1,1:3)**2)+ sum(tr(2,1:3)**2)) - 2*dfi(3)*dfi(1)* sum(tr(3,1:3)*tr(1,1:3) )        )/area_F(i)

enddo

E2=0.0
do i=1,N_V
E2=E2+ ( -a2/2.0*fi(i)**2+a4/4.0*fi(i)**4.0  )*area_V(i)/3.0
enddo


phase_energy=E1+E2


END