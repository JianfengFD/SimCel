



!  ************************************************
!  Part of 3D vesicle program
! Read *.fe file created by surface evolver
!                    written by LiJF  Dobb@bbs.fudan.edu.cn
! ************************************************* 

module read_mod
 use CelType_mod
 use basicfun_mod
    contains
    subroutine print_screen_cell(C,k,t)
          implicit none
          type(Cel)::C(:)
          integer,intent(in)::k,t
                call Get_shape_information(C,k)
                write(*,*)t,C(k)%energy_tot,C(k)%energyH &
                /8.0/PAI/C(k)%H_modulus
                write(*,*)C(k)%area_tot,C(k)%area_tot_zero,C(k)%vol_total, &
                    C(k)%vol_total/ ( C(k)%area_tot/4.0/pai)**(3.0/2.0)*3.0/4.0/pai
             
       end subroutine print_screen_Cell


       subroutine save_Cell(C,kk,t,file_in)
          implicit none
          type(Cel)::C
          character*3 file_in
          character*2 kc
          integer i,kk,t
          ! %% save the shapes
            write(kc,'(I2.2)')t
            open(7,FILE='data/R'//FILE_IN//kc//'.dat')
            write(7,*)'//',kk,t
            
            write(7,'(A8)')'vertices'
            do i=1,C%N_V
            write(7,'(I5,3(F30.15))')i,C%rv(i,1),C%rv(i,2),C%rv(i,3)
            enddo
            

            write(7,*)
            write(7,'(A5)') 'edges'
            do i=1,C%N_E
            write(7,'(3(I7),A8,I6)')i,C%V_E(i,1),C%V_E(i,2),'color',8
            enddo

            write(7,*)
            write(7,'(A5)') 'faces'

            do i=1,C%N_F
            write(7,'(4(I7), A8, I6)')i,C%E_F(i,1),C%E_F(i,2),C%E_F(i,3),'color',6
            enddo

            write(7,*)

            write(7,*)'        '
            close(7)
            
        
            
       end subroutine Save_cell
       
       
       
       subroutine read_cell(C)
          implicit none
          type(Cel)::C
          call readfe(C%FILE_Cel,C%rv,C%N_V,C%N_E, &
           C%N_F,C%V_E,C%E_F)          
       end subroutine read_cell
       
       subroutine read_cell_size(C)
          implicit none
          type(Cel)::C
          call Read_Mesh_Size(C%FILE_Cel, C%N,C%N_V,C%N_E,C%N_F)       
       end subroutine read_cell_size
       
        subroutine readfe(FILE_NAME,rv,N_V,N_E,N_F,V_E,E_F)

        implicit none
        integer N
        integer i,j,k,j1,j2,L,k1
        integer switch_I, k_v, k_E,k_F, G_F
        real*8 a

         real*8 rv(:,:)
         integer  V_E(:,:),E_F(:,:)
 

        character*150 line_read
        character*150 vchar
        character*1 sp
        integer N_V,N_E,N_F
        character*100 trim_read
        character*8 FILE_NAME



        ! %%% read file
        open(7,FILE="data/initials/"//FILE_NAME)
        switch_I=1
        k_v=0;k_E=0;k_F=0
        G_F=0
        sp=' '

        do while(switch_I.eq.1)


        read(7,'(A100)')line_read



        trim_read=trim(line_read)



        if(len(trim(trim_read)).eq.0)then
        k_V=0;k_E=0;k_F=0
        endif



        ! Vertices
        ! %%%%%%%%%%%%%%%%%%%%%%%%%%%%
        if(k_v.eq.1)then

        j=0
        j1=1
        j2=0


        do i=1,len(trim_read)

        if(j.eq.0.and.trim_read(i:i).eq.' ')then
        j=1
        ! READ  k


        if((j2.eq.1.or.trim_read(1:1).ne.' ').and.j2.lt.2)then
        vChar=trim_read(j1:i-1)
        L=i-j1
        call char2Num(L,vChar,k)
        j2=1

        endif


        ! READ rv
        if(j2.ge.2)then

        vChar=trim_read(j1:i-1)
        L=i-j1
        call char2real(L,vChar,a)
        rv(k,j2-1)=a

        endif
        !
        j2=j2+1
        endif

        ! %%

        if(j.eq.1.and.trim_read(i:i).ne.' ')then
        j=0
        j1=i
        endif

        enddo
        N_V=k
        endif

        ! %%%%%%%%  EDGES

        if(k_E.eq.1)then
        j=0
        j1=1
        j2=0
        do i=1,len(trim_read)

        if(j.eq.0.and.trim_read(i:i).eq.' ')then
        j=1
        ! READ  k

        if((j2.eq.1.or.trim_read(1:1).ne.' ').and.j2.lt.2)then
        vChar=trim_read(j1:i-1)
        L=i-j1
        call char2Num(L,vChar,k)
        j2=1
        endif

        ! READ rv
        if(j2.eq.2.or.j2.eq.3)then

        vChar=trim_read(j1:i-1)
        L=i-j1
        call char2NUM(L,vChar,k1)
        V_E(k,j2-1)=k1
        endif
        !
        j2=j2+1
        endif

        ! %%

        if(j.eq.1.and.trim_read(i:i).ne.' ')then
        j=0
        j1=i
        endif
        enddo
        N_E=k
        endif


        ! %%%%%%%%%%%%%% faces

        if(k_F.eq.1)then
        j=0
        j1=1
        j2=0
        do i=1,len(trim_read)

        if(j.eq.0.and.trim_read(i:i).eq.' ')then
        j=1
        ! READ  k

        if((j2.eq.1.or.trim_read(1:1).ne.' ').and.j2.lt.2)then
        vChar=trim_read(j1:i-1)
        L=i-j1
        call char2Num(L,vChar,k)
        j2=1
        endif

        ! READ rv
        if(j2.ge.2.and.j2.le.4)then

        vChar=trim_read(j1:i-1)
        L=i-j1
        call char2NUM(L,vChar,k1)


        E_F(k,j2-1)=k1
        endif
        !
        j2=j2+1
        endif

        ! %%

        if(j.eq.1.and.trim_read(i:i).ne.' ')then
        j=0
        j1=i
        endif
        enddo
        N_F=k
        endif


        ! %%%%%%%%%%%%%%%%%%%%%%%


        !  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        if(trim_read(1:8).eq.'vertices')then
        k_v=1;k_E=0;k_F=0

        endif
        if(trim_read(1:5).eq.'edges')then
        k_v=0;k_E=1;k_F=0
        endif
        if(trim_read(1:5).eq.'faces')then
        k_v=0;k_E=0;k_F=1; G_F=1
        endif

        if(G_F.eq.1.and.k_F.eq.0)switch_I=0


        enddo

        close(7)

            END subroutine readfe


    
    
subroutine Read_Mesh_Size(FILE_NAME, N,N_V,N_E,N_F)

implicit none
 
integer N,N_V,N_E,N_F
integer i,j,k,j1,j2,L,k1
integer switch_I, k_v, k_E,k_F, G_F
real*8 a

character*150 line_read
character*150 vchar
character*1 sp

character*100 trim_read
character*8 FILE_NAME

! %%%%%%%%%%%%%  Get N_V,N_E,N_F
open(7,FILE="data/initials/"//FILE_NAME)
switch_I=1
k_v=0;k_E=0;k_F=0
G_F=0
sp=' '
do while(switch_I.eq.1)
  read(7,'(A100)')line_read
  trim_read=trim(line_read)
  if(len(trim(trim_read)).eq.0)then
    k_V=0;k_E=0;k_F=0
  endif
  ! Vertices
  ! %%%%%%%%%%%%%%%%%%%%%%%%%%%%
  if(k_v.eq.1)then
    j=0;j1=1;j2=0
    do i=1,len(trim_read)
       if(j.eq.0.and.trim_read(i:i).eq.' ')then
        j=1
        ! READ  k
        if((j2.eq.1.or.trim_read(1:1).ne.' ').and.j2.lt.2)then
            vChar=trim_read(j1:i-1);L=i-j1
            call char2Num(L,vChar,k)
            j2=1
        endif

         j2=j2+1
      endif
      if(j.eq.1.and.trim_read(i:i).ne.' ')then
        j=0;j1=i
      endif
   enddo
   N_V=k
 endif
! %%%%%%%%  EDGES
  if(k_E.eq.1)then
    j=0;j1=1;j2=0
    do i=1,len(trim_read)
        if(j.eq.0.and.trim_read(i:i).eq.' ')then
            j=1
        ! READ  k
            if((j2.eq.1.or.trim_read(1:1).ne.' ').and.j2.lt.2)then
                vChar=trim_read(j1:i-1);L=i-j1
                call char2Num(L,vChar,k)
                j2=1
            endif
       

            j2=j2+1
        endif
        if(j.eq.1.and.trim_read(i:i).ne.' ')then
            j=0;j1=i
        endif
    enddo
N_E=k
endif
! %%%%%%%%%%%%%% faces
    if(k_F.eq.1)then
        j=0;j1=1;j2=0
        do i=1,len(trim_read)
            if(j.eq.0.and.trim_read(i:i).eq.' ')then
                j=1
                ! READ  k
                if((j2.eq.1.or.trim_read(1:1).ne.' ').and.j2.lt.2)then
                vChar=trim_read(j1:i-1)
                L=i-j1
                call char2Num(L,vChar,k)
                j2=1
                endif


                if(j2.ge.2.and.j2.le.4)then
                    vChar=trim_read(j1:i-1)
                    L=i-j1
                    call char2NUM(L,vChar,k1)
                endif
 
                j2=j2+1
            endif
            if(j.eq.1.and.trim_read(i:i).ne.' ')then
                j=0;j1=i
            endif
        enddo
        N_F=k
    endif


! %%%%%%%%%%%%%%%%%%%%%%%


!  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if(trim_read(1:8).eq.'vertices')then
k_v=1;k_E=0;k_F=0
endif
if(trim_read(1:5).eq.'edges')then
k_v=0;k_E=1;k_F=0
endif
if(trim_read(1:5).eq.'faces')then
k_v=0;k_E=0;k_F=1; G_F=1
endif
if(G_F.eq.1.and.k_F.eq.0)switch_I=0
enddo

close(7)

N=N_V
if(N.lt.N_E)N=N_E
if(N.lt.N_F)N=N_F




END subroutine Read_Mesh_Size

    
    



subroutine char2Num(L,vChar,k)

implicit none
integer k, L
integer a,i,signI

character*150 vChar
character achar1
k=0
!*************
signI=0
if(vchar(1:1).eq.'-')then
signI=1
endif

do i=1+signI,L
achar1=vCHar(i:i)
a=iachar(achar1)
k=k+(a-48)*10**(L-i)
enddo
if(signI.eq.1)k=-k



END subroutine char2Num


subroutine char2real(L,vChar,a)

implicit none
integer L
integer i

integer switch_E, I_E, L1,L2

real*8 a_E2,a_E1
real*8 a

character*150 vChar
character*150 vChar_E1,vChar_E2
character achar1

!!!!!!!!!!!!!!!!!!!!
!!  let 1.988e-8 be zero
switch_E=0
I_E=0

do i=1,L
achar1=vchar(i:i)
if(achar1.eq.'e'.or.achar1.eq.'E')then
switch_E=1
I_E=i
endif

enddo
! 
a_E1=0.0
a_E2=0.0
a=0.0
if(Switch_E.eq.1)then
L1=I_E-1
vChar_E1(1:I_E-1)=vChar(1:I_E-1)
L2=L-i_E
vChar_E2(1:L2)=vChar(I_E+1:L)

if(L1.ge.1)call char2real_No_E(L1,vChar_E1,a_E1)
if(L2.ge.1)call char2real_No_E(L2,vChar_E2,a_E2)
a=a_E1*10**(a_E2)
else

if(L.ge.1)call char2real_No_E(L,vChar,a)
endif
END subroutine char2real


subroutine char2real_No_E(L,vChar,a)


implicit none
integer k, L
integer i,i_dec,a_I,j,signI


real*8 a_d, a

character*150 vChar
character achar1

signI=0

if(vchar(1:1).eq.'-')then
signI=-1
endif
if(vchar(1:1).eq.'+')then
signI=1
endif

I_Dec=L
do i=1,L
achar1=vchar(i:i)
if(achar1.eq.'.')then
i_dec=i-1
endif
enddo

! INT part
a_I=0
do i=1+abs(signI),i_Dec
achar1=vchar(i:i)
k=iachar(achar1)-48
a_I=a_I+k*10**(i_Dec-i)
enddo

! decimal part
a_d=0.0
do i=i_dec+2,L
achar1=vchar(i:i)
j=i-i_dec
k=iachar(achar1)-48
a_d=a_d+k*10.0**(i_dec+1-i)
enddo

a=a_I+a_d
if(signI.eq.-1)a=-a

END subroutine char2real_No_E

    end module read_mod
    
    
    
    
    
    
    
    
    
    
    
    