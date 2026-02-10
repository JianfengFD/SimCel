
!  ************************************************
!  Part of 3D vesicle program
! label the surface
!                    written by LiJF  Dobb@bbs.fudan.edu.cn
! *************************************************  
!! include 'readfe.f90'


 subroutine label_system(V_F,   E_V,N_E_V,     F_V,N_F_V,    F_E,N_F_E,    V_V,N_V_V,   E_V_opposite,N_E_V_opposite,   V_E_opposite,  N_V,N_E,N_F,V_E,E_F,F_F,rv)

parameter(N = 20000)

integer N_V,N_E,N_F

real*8 rv(0:N,1:3),x,y,z,a

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

integer i,j,k,i1,i2,i3,i4,i5,i0,I_temp(1:20),I_temp2(1:20),I_temp3(1:20)





! the first kind label
! (1) V_F
  do i=1,N_F
  V_F(i,1)=V_E(abs(E_F(i,1)), int( (3-E_F(i,1)/abs(E_F(i,1)) ) /2 ) )
  V_F(i,2)=V_E(abs(E_F(i,2)), int( (3-E_F(i,2)/abs(E_F(i,2)) ) /2 ) )
  V_F(i,3)=V_E(abs(E_F(i,3)), int( (3-E_F(i,3)/abs(E_F(i,3)) ) /2 ) )
  enddo

! (2) E_V
!       Problem: E_V Ã»ï¿½Ô±ï¿½ï¿½ï¿½ï¿½ï¿½ï¿½ï¿½ï¿½ï¿½ï¿?
	N_E_V=0
  do i=1, N_E
  i1=V_E(i,1); i2=V_E(i,2)
  N_E_V(i1)=N_E_V(i1)+1
  N_E_V(i2)=N_E_V(i2)+1
  E_V(i1, N_E_V(i1))=i
  E_V(i2, N_E_V(i2))=-i
  enddo
! (3) F_V ï¿½í¡¡ï¿½ï¿½ï¿½ï¿½ï¿½ÃµÄ£ï¿½ï¿½ï¿½ï¿½ï¿½
  N_F_V=0
  do i=1, N_F
  i1=V_F(i,1);i2=V_F(i,2);i3=V_F(i,3)
  N_F_V(i1)=N_F_V(i1)+1
  N_F_V(i2)=N_F_V(i2)+1
  N_F_V(i3)=N_F_V(i3)+1
  F_V(i1, N_F_V(i1) )=i
  F_V(i2, N_F_V(i2) )=i
  F_V(i3, N_F_V(i3) )=i
  enddo

! (4) F_E
N_F_E=0
do i=1, N_F
  i1=E_F(i,1);i2=E_F(i,2);i3=E_F(i,3)
  N_F_E(abs(i1))=N_F_E(abs(i1))+1
  N_F_E(abs(i2))=N_F_E(abs(i2))+1
  N_F_E(abs(i3))=N_F_E(abs(i3))+1
  F_E(abs(i1), int( (3-i1/abs(i1) ) /2 ) )=i
  F_E(abs(i2), int( (3-i2/abs(i2) ) /2 ) )=i
  F_E(abs(i3), int( (3-i3/abs(i3) ) /2 ) )=i
enddo

! (5) F_F
do i=1, N_F
i1=abs(E_F(i,1));i2=abs(E_F(i,2));i3=abs(E_F(i,3))
i4=F_E(i1,1);i5=F_E(i1,2)
if(i4.eq.i)F_F(i,1)=i5
if(i5.eq.i)F_F(i,1)=i4
i4=F_E(i2,1);i5=F_E(i2,2)
if(i4.eq.i)F_F(i,2)=i5
if(i5.eq.i)F_F(i,2)=i4
i4=F_E(i3,1);i5=F_E(i3,2)
if(i4.eq.i)F_F(i,3)=i5
if(i5.eq.i)F_F(i,3)=i4
enddo
 

! %%%%%% second kind label
! (6)  E_V_opposite
do i=1,N_V
do j=1,N_F_V(i)
k=F_V(i,j)
i1=E_F(k,1); i2=E_F(k,2); i3=E_F(k,3); 
if( (V_E(abs(i1),1)-i)*(V_E(abs(i1),2)-i).ne.0) E_V_opposite(i,j)=i1
if( (V_E(abs(i2),1)-i)*(V_E(abs(i2),2)-i).ne.0) E_V_opposite(i,j)=i2
if( (V_E(abs(i3),1)-i)*(V_E(abs(i3),2)-i).ne.0) E_V_opposite(i,j)=i3
enddo
enddo
N_E_V_opposite=N_F_V
   !ï¿½ï¿½

   do i=1, N_v
   I_temp(1)=1
   I_temp2(1:N_F_V(i))=E_V_opposite(i,1:N_F_V(i) )
   I_temp3(1:N_F_V(i))=F_V(i,1:N_F_V(i) )
   i1=E_V_opposite(i,1)
   i2=V_E(abs(i1), int( (3-i1/abs(i1) ) /2 ) ) 
   i0=i2;  i5=i2+10
   i3=V_E(abs(i1), int( (3+i1/abs(i1) ) /2 ) ) 
   i4=i2
   	k=2  ! For I_temp's label
	do while(i3.ne.i0)

	j=2  ! For 1- 6 iteration
 
    do while(i4.ne.i3.and.j.le.N_F_V(i))
	
   i1=E_V_opposite(i,j)
   i4=V_E(abs(i1), int( (3-i1/abs(i1) ) /2 ) ) 
   i5=V_E(abs(i1), int( (3+i1/abs(i1) ) /2 ) )
     j=j+1
   	enddo
    i3=i5
	I_temp(k)=j-1
	k=k+1
	enddo

    do j=1,N_F_V(i)
	E_V_opposite(i,j)=I_temp2( I_temp(j) )
	F_V(i,j) = I_temp3( I_temp(j) )  ! F_V ï¿½ï¿½ï¿?
    enddo
  enddo
!  (2) E_V & (5) V_V ï¿½ï¿½

	do i=1,N_V
	I_temp2(1:N_F_V(i))=E_V(i,1:N_F_V(i))
	do j=1,N_F_V(i)
		i1=E_V_opposite(i,j)
		i2=V_E(abs(i1), int( (3-i1/abs(i1) ) /2 ) ) 
		V_V(i,j)=i2
		do k=1,N_F_V(i)
		i3=I_temp2(k)
		i4=V_E(abs(i3), int( (3+i3/abs(i3) ) /2 ) )
		if(i4.eq.i2)E_V(i,j)=i3
		enddo
	enddo
	enddo
	N_V_V=N_F_V

! (7) V_E_opposite  use F_E 4ïµ?

   do i=1,N_E
   do j=1,2
    i1=F_E(i,j);
	i2=V_E(i,1)
	i3=V_E(i,2)
	do k=1,3
	 i4=V_F(i1,k)
	 if(i4.ne.i2.and.i4.ne.i3)V_E_opposite(i,j)=i4
	enddo
   enddo
   enddo




END

