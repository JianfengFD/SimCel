!  ************************************************
!  Part of 3D vesicle program
! Some functions
!                    written by LiJF  Dobb@bbs.fudan.edu.cn
! *************************************************  



! input R1, R2 output: |R1\\R2|^alpha
Real*8 Function R1R2_fun(R1,R2,I_c)
real*8 R1(1:3), R2(1:3)
integer I_C,i,j,k

select case(I_c)
case(1)! x (y z)
            
R1R2_fun=R1(2)*R2(3)-R1(3)*R2(2)

case(2) ! y (z x)
        
R1R2_fun=R1(3)*R2(1)-R1(1)*R2(3)

case(3) ! z (x y)

R1R2_fun=R1(1)*R2(2)-R1(2)*R2(1)

case default
write(*,*) 'wrong input within RR_fun'
pause
end select




END

! input R1, R2, R3 output: ||R1\\R2| \\ R3 |^alpha
Real*8 Function R1R2_R3_fun(R1,R2,R3,I_c)
real*8 R1(1:3), R2(1:3), R3(1:3)
real*8 D1(1:3),R1R2_fun
integer I_C,i,j,k

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
pause
end select

END

! input R1, R2, R3 output: ||R1\\R2| \\ R3 |^alpha
Real*8 Function R1_R2R3_fun(R1,R2,R3,I_c)
real*8 R1(1:3), R2(1:3), R3(1:3)
real*8 D1(1:3),R1R2_fun
integer I_C,i,j,k

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
pause
end select

END


! input R(1:3,1:3), I_c(x,y, or z) output : | |ex+| |ey +| |ez

Real*8 function  RR_fun(R, I_c) 

       ! (P_Num ,Cordinate )
Real*8    R(1:3,  1:3)
integer I_c,i,j,k

select case(I_c)
case(1)! x (y z)
        !  1x                               2x(3 1)                          3x(1 2)    
RR_fun= R(2,2)*R(3,3)-R(2,3)*R(3,2)  + R(3,2)*R(1,3)-R(3,3)*R(1,2) +   R(1,2)*R(2,3)-R(1,3)*R(2,2) 


case(2) ! y (z x)
        !  1y                               2y(3 1)                          3y(1 2)    
RR_fun= R(2,3)*R(3,1)-R(2,1)*R(3,3)  + R(3,3)*R(1,1)-R(3,1)*R(1,3) +   R(1,3)*R(2,1)-R(1,1)*R(2,3) 

case(3) ! z (x y)
        !  1z                               2z(3 1)                          3z(1 2)    
RR_fun= R(2,1)*R(3,2)-R(2,2)*R(3,1)  + R(3,1)*R(1,2)-R(3,2)*R(1,1) +   R(1,1)*R(2,2)-R(1,2)*R(2,1) 
case default
write(*,*) 'wrong input within RR_fun'
pause
end select


END

!  input R(1:3,1:3) Output |R\\R\\R|  Volume function

real *8 Function V6_fun(R)

Real*8    R(1:3,  1:3)
integer I_c,i,j,k

V6_fun=R(1,1)*( R(2,2)*R(3,3)-R(2,3)*R(3,2) )+R(1,2)*(R(2,3)*R(3,1)- R(2,1)*R(3,3))+R(1,3)*(R(2,1)*R(3,2)-R(2,2)*R(3,1) )


END


! laplace_1





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




