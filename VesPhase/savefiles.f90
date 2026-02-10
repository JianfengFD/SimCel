
!  ************************************************
!  Part of 3D vesicle program
! Read *.fe file created by surface evolver
!                    written by LiJF  Dobb@bbs.fudan.edu.cn
! *************************************************  

subroutine savefiles(N_V,N_E,N_F,V_E,E_F,V_F,F_E,F_F ,rv,fi,fia0,t,real_t,H_V,area_V,EnergyP,EnergyH,energyL, energy_tot,area_F, area_tot,area_tot_zero,vol_total,V_E_opposite, switch_boundary_setting,dt,C_energy_array,Ph_energy_array,energy_array,Num_b_array,Num_d_array,len_B_array,tk,  switch_domain_process,area_domain_files,ratio_area_len_files,len_domain_files,t_files,k_files,  area_domain0, Num_domain0,F_domain0,N_F_domain0,No_domain0,No_domain_Q0,inv_No_domain0,Max_No0)
implicit none
integer N
parameter(N = 20000)

integer N_V,N_E,N_F

real*8 rV(0:N,1:3),x,y,z,a
real*8 fi(0:N),H_V(0:N),area_V(0:N)
real*8 area_F(0:N)

integer V_E(0:N, 1:2)
integer E_F(0:N, 1:3)
integer V_F(0:N,1:3)
integer F_E(0:N, 1:20)
integer F_F(0:N,1:3)
integer I_fi_F(0:N)
real*8 fi_F(0:N)

real*8 PV_energy, Curvature_energy,area_hooke_energy,phase_energy
real*8 C_energy_array(1:50000),Ph_energy_array(1:50000),energy_array(1:50000)
real*8 EnergyP,EnergyH,energyL, energy_tot
real*8 area_tot,area_tot_zero,vol_total
integer t
real*8 real_t,dt
integer t_int,tk
integer i,j,k,i1,i2,i3,i_V,i_E,i_F,j1,j2,j3,j4  ,L,k1,k2,k3

real*8 fia0
real*8 len_boundary,len_B_array(1:50000)
real*8 boundary_position(0:N,1:6)
integer V_E_opposite(0:N,1:20)
integer Num_boundary,switch_boundary_setting,Num_domain
integer Num_b_array(1:50000),Num_d_array(1:50000)
integer Circle_F(1:100,1:50000),Point_circle(1:1000)

integer all_F_domain(1:N)

real*8 len_domain(1:1000)
real*8 area_domain(1:1000)
real*8 area_domain0(1:1000)
integer Num_domain0, F_domain0(1:100,1:50000),N_F_domain0(1:1000)
integer No_domain0(1:1000,1:100), No_domain_Q0(1:1000)
  ! i(1:1000) for 1:Num_domain0,   return Value for actual domain ID 
  ! k(1:100) for overload domain ID, No_domain_Q0 is for overload numbers
integer inv_No_domain0(1:1000)
integer Max_No0  ! number of output files
real*8 ratio_area_len(1:1000)
! for domain savefiles

real*8 area_domain_files(1:50000,1:100)
real*8 ratio_area_len_files(1:50000,1:100)
real*8 len_domain_files(1:50000,1:100)
real*8 t_files(1:50000)  ! to record the time
integer k_files, switch_domain_process



character*14 filename
character*15 filename2
character*5 ext


 call convert_FNUM(fi, I_fi_F, fi_F, V_F, N_F)

call boundary_count(rv,fi,Len_boundary,fia0,N_E,F_E,V_E,V_E_opposite, switch_boundary_setting,boundary_position, Num_boundary,Num_domain,Circle_F,len_domain,Point_circle,all_F_domain)

switch_domain_process=1
if(switch_domain_process.eq.1)then

if(tk.lt.10.or.t.lt.10)then
area_domain_files=0.0
len_domain_files=0.0
ratio_area_len_files=0.0
k_files=1
Num_domain0=0
Max_No0=0
endif

if(tk.ge.200.and.t.gt.20.)then

call  domain_area(area_F,fi_F,F_F,N_F,Num_domain,Circle_F,len_domain,Point_circle,area_domain,switch_boundary_setting,fia0,  area_domain0, Num_domain0,F_domain0,N_F_domain0,No_domain0,No_domain_Q0,inv_No_domain0,Max_No0,ratio_area_len)

do i=1,Max_No0
k=inv_No_domain0(i)
if(k.eq.0)then
area_domain_files(k_files,i)=0.0
len_domain_files(k_files,i)=0.0
ratio_area_len_files(k_files,i)=0.0
endif
if(k.ne.0)then
area_domain_files(k_files,i)=area_domain(k)
len_domain_files(k_files,i)=len_domain(k)
ratio_area_len_files(k_files,i)=ratio_area_len(k)
endif
enddo
t_files(k_files)=real_t
write(*,*)area_domain_files(k_files,1:Max_no0)
k_files=k_files+1
endif

endif

if(mod(tk,200).eq.0.and.tk.gt.95)then
do i=1,Max_No0
write(ext,999)i
filename2='data/domain/are'
open(7,FILE=filename2//ext//'.dat')
do j=1,k_files-1
write(7,*)t_files(j),area_domain_files(j,i)
enddo
close(7)
filename2='data/domain/len'
open(7,FILE=filename2//ext//'.dat')
do j=1,k_files-1
write(7,*)t_files(j),len_domain_files(j,i)
enddo
close(7)
filename2='data/domain/rat'
open(7,FILE=filename2//ext//'.dat')
do j=1,k_files-1
write(7,*)t_files(j),ratio_area_len_files(j,i)
enddo
close(7)
enddo
open(7,FILE='data/domain/N_D.dat')
do i=2,tk-1
write(7,*)i,Num_d_array(i)
enddo

endif


k=tk
if(k.ge.1)then
 C_energy_array(k)=EnergyH
 Ph_energy_array(k)=energy_tot-energyp-energyH
 energy_array(k)=energy_tot
 len_B_array(k)=len_boundary
 Num_d_array(k)=Num_domain
 Num_b_array(k)=Num_boundary
endif

999 format(I5.5)

if(mod(tk,10).eq.1)then
write(ext,999)tk

filename='data/dynamic/R'
open(7,FILE=filename//ext//'.dat')
do i=1,N_V
write(7,'(I5,3(F30.15))')i,rv(i,1),rv(i,2),rv(i,3)
enddo
write(7,*)'        '

close(7)
filename='data/dynamic/f'
open(7,FILE=filename//ext//'.dat')
do i=1,N_F
write(7,*) i,(fi(V_F(i,1))+fi(V_F(i,2))+fi(V_F(i,3)))/3.0
enddo
close(7)
filename='data/dynamic/V'
open(7,FILE=filename//ext//'.dat')
do i=1,N_F
write(7,*) i,V_F(i,1),V_F(i,2),V_F(i,3)
enddo
close(7)

endif
if(mod(tk,1).eq.0)then
write(*,*)'================'
write(*,*)t,real_t,dt
write(*,*)'energyp      ','energyH      ','energyPh'
write(*,*)energyp,energyH,energy_tot-energyH-energyp
write(*,*)'energy_total   ',energy_tot,"vol ", vol_total
write(*,*)'area  ',area_tot,area_tot_zero
write(*,*)'len_boudary  ',Num_boundary,Num_domain
endif

end