!  ************************************************
!  Main programm of 3D vesicle program
!                    written by LiJF  lijf@fudan.edu.cn, 2014-12
program main
use head_Cell_mod
 implicit none
! %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! SEC. O        Parameter Defintions
 ! %%%%%%%%%%%
 ! FOR PHASE TRANSITION
 real*8 integral_H,integral_H0
 real*8 Diff_int_H,kpp_int_H,det_int_H
 integer kk,NN
 
! %%%%%%%%%%%%
! ORIGINAL
 integer N
 integer i,i1,i2,j,k,I_FLIP_TOT
 real*8 dt,dt0,real_t, rv_center(1:3),R0,kBT, gauss_kbt,inverse_kbt,rc0(1:3)
 real*8 D_thres
 integer*4 t,tk;
 integer big_move,seed1
 real*8 m0;character*5 cm0;
 character*3  file_tmp
 character*3  FILE_INPUT
 character*8  file_RBC,file_plane
 character*2  v0_c
 character*1  type_C,seed_char
 character*150 kBT_char,et_char
 ! %%%% Cels
 Type(Cel),allocatable:: Cels(:)    ! the Cel object defined in Celtype_Mod in Module_CellType.f90
 Type(Cel)::Cout,C_AVE(1:2)
 Type(PtCell)::Vm  
 integer ToT_NCel, NC
 real*8,allocatable::alN1(:),betN1(:)
 integer,allocatable::kfN1(:)
 real*8,allocatable::alN2(:),betN2(:)
 integer,allocatable::kfN2(:)
 
 real PP ! the rate of varying the volume
 real*8 sig_LJ,ep_LJ,A_LJ,B_LJ,eta
 real*8 lam0
 real*8 R00,R01
 
 real*8 L_box,dd
 
 real*8 et_in,kbt_in
 
 integer kin_cells
 integer,allocatable::active_cells(:)
 
 real*8 L0
 
 ! %%%% 
 real*8 f_down,f_rescale
 real*8 dtk(1:2)
 real*8 EHT(1:10000)
 real*8 EH
 real*8 contf
 
 ! %%%%%%%%%%
 ! set the Z dynamics
 real*8 Z_Stop(1:50),Z_Fixed
 integer N_Z_STOP, t1,t2,IS_T1,IS_T2,k_Z,t_extra
 
 real*8 fracALL(1:50),EALL(1:2,1:50) ! H_plane,H_sphere
 real*8 frac_mean
 integer k_frac,IS_frac
 
 integer IS_READ
 
 integer t_NV
 ! adhesion and repulsion
 real*8 eps,sig_eps
 
 
 
 eps=0.01*80
 sig_eps=0.28 ! F=12* eps*(2*sig_eps**6/d**7-sig_eps**3/d**4)*n0
 
 
 IS_READ = 1

 f_rescale = 0.8  ! % f_down's value is given after R0 is known in the k=2 section.
 ! file_RBC = 'R162.dat'
  !file_plane='PLN2.dat'
file_input ='PLN'
   i=1
  !call getarg(file_input,i)
   i=2
   kbt_char='0.05'
  !call getarg(kbt_char,i)
  i=3
  et_char = '2.00'
  !call getarg(et_char,i)
  i=4
  seed_char='1'
  !call getarg(seed_char,i)
  
  
  
  
  call char2real(4,kbt_char,kBT_in)
  
  call char2real(4,et_char,eT_in)
 
  
  
  
  
  lam0 = 0.5
  R00 = 0.2
  R01 = 2.0
  
  eta = eT_in*lam0*10*R00**4/PAI
  
  if(seed_char.eq.'1')seed1 = 26345678
  if(seed_char.eq.'2')seed1 = 56345678
  L_box = 50

! %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 
! SEC. I        initial state construction
    ! ####################################
    ! SEC. I-1  define the cells
    ! ####################################
    ! There are ToT_NCel cells in the system
        
        ToT_NCel = 2
    allocate(Cels(ToT_NCel),active_cells(ToT_NCel))
    
    
    ! ####################################
    ! SEC. I-2  Read the Cell
    ! ####################################
     Cels(1)%FILE_Cel = 'CIR1.dat'
     Cels(2)%FILE_Cel = 'R642.dat'
     
     k=2
     call Read_Cell_Size(Cels(k))       
     call alloc_Cel(Cels,k,Cels(k)%N,Tot_NCel)
     C_AVE(k)%N_V=Cels(k)%N_V;C_AVE(k)%N_E=Cels(k)%N_E;C_AVE(k)%N_F=Cels(k)%N_F;
     call alloc_SglCel(C_AVE(k),Cels(k)%N)
     call read_Cell(Cels(k))
     R0=  23.75*(Cels(k)%N_v/10242.0)**(1.0/2.0) 
     
       k=1
     call Read_Cell_Size(Cels(k))  
    
 
     call alloc_Cel(Cels,k,Cels(k)%N,Tot_NCel)
     C_AVE(k)%N_V=Cels(k)%N_V;C_AVE(k)%N_E=Cels(k)%N_E;C_AVE(k)%N_F=Cels(k)%N_F;
     
     call alloc_SglCel(C_AVE(k),Cels(k)%N)
     call read_Cell(Cels(k))
    
     L0=  (4*3.1415926*23.75**2*Cels(k)%N_v/10242.0)**0.5
     Cels(k)%Rv =Cels(k)%Rv*L0 
     
          kBT = kBT_in*R0
   

     ! move the plane to (0,0,0)
     rv_center(1:3) = 0.0
     rv_center(1) =0.0
     k=1
     call Move_Cell(Cels,k,Rv_center)  
     
     
     ! %%%%%%
     ! INITIALIZE the soft particle
     ! %%%%%%

     k=2
       
       call label_Cell(cels,k)
       call Get_shape_information(Cels,k)   
    ! ####################################
    ! SEC. I-3  Adjust the Cell's Position
    !       & Rescale the Cell size
    ! ####################################
         rv_center(1:3) = 0.0
         rv_center(1) =0.0
        call Move_Cell(Cels,k,Rv_center)        
        R0=  23.75*(Cels(k)%N_v/10242.0)**(1.0/2.0)  ! /0.89  ! 23.75 for 10242(140.0/4/PAI)**0.5
        f_down = f_rescale/R0
        
        call Resize_Cell(Cels,k,R0)
        
    !%% Place the soft particle right above the plane. Delta =2.0    
        Rv_center(3)=  R0 + 0.5
        do i=1,50
           Z_STOP(i)=(1.0-(i-1)*0.2)*R0+0.5 
           !Z_STOP(i)=(1.0-(i-1)*0.5)*R0+0.5 
        enddo
        t1= 4000;t2=1000
        call Move_Cell(Cels,k,Rv_center)  
    ! ####################################
    ! SEC. I-4  Get Cell's Geom. information:
    !           Len,area,Vol,H,Energy
    ! ####################################
     
      call Get_shape_information(Cels,k)
        Cels(k)%area_F_zero=Cels(k)%area_F;Cels(k)%area_tot_zero = Cels(k)%area_tot
        Cels(k)%Vol_total_zero = Cels(k)%Vol_total
    ! ####################################
    ! SEC. I-5  Set Cell membrane' Elasticities:
    ! ####################################
        Cels(k)%H_modulus = 1.0*6/16.0
        Cels(k)%kpp_V =    5.0*10.0**4.0*Cels(k)%H_modulus*2 &
                      /(Cels(k)%vol_total_zero/100.0)
        Cels(k)%lamda=   lam0/Cels(k)%area_tot_zero  ! %%%% kBT
        Cels(k)%kpp_area = 5.0*10.0**5.0*Cels(k)%H_modulus*2 &
                            /(Cels(k)%area_tot_zero/140.0) /5.0
        Cels(k)%area_rel = PAI*(Cels(k)%vol_total_zero/4.0/PAI*3.0)**(2.0/3)*4 *2.0
      !  Cels(k)%Vol_rel = 0.29812* 4.0/3.0*PAI *(cels(k)%area_tot/4/pai)**1.5
       
        Cels(k)%H_zero=    -0 /(Cels(k)%area_tot/4.0/PAI)**(0.5)
        Cels(k)%eps=eps
        Cels(k)%sig_eps=sig_eps
        
    ! %%%%%%
     ! END of INITIALIZation of the soft particle
     ! %%%%%%
    
    ! %%%%%%
     ! INITIALIZE the Plane
     ! %%%%%%

     k=1
       call label_Cell(cels,k)
       call Get_shape_information(Cels,k)   
    ! ####################################

        Cels(k)%area_F_zero=Cels(k)%area_F;Cels(k)%area_tot_zero = Cels(k)%area_tot
        Cels(k)%Vol_total_zero = Cels(k)%Vol_total

        Cels(k)%H_modulus = 1.0 ! kpp1=6; kpp2 =16
        Cels(k)%kpp_V =    5.0*10.0**4.0*Cels(k)%H_modulus*2 &
                      /(Cels(k)%vol_total_zero/100.0)
        Cels(k)%lamda=   lam0/Cels(k)%area_tot_zero  ! %%%% kBT
        Cels(k)%kpp_area = 5.0*10.0**5.0*Cels(k)%H_modulus*2 &
                            /(Cels(k)%area_tot_zero/140.0) /5.0
        Cels(k)%area_rel=PAI*(Cels(k)%vol_total_zero/4.0/PAI*3.0)**(2.0/3)*4 *2.0
      !  Cels(k)%Vol_rel = 0.29812* 4.0/3.0*PAI *(cels(k)%area_tot/4/pai)**1.5
       
        Cels(k)%H_zero=    -0 /(Cels(k)%area_tot/4.0/PAI)**(0.5)
        Cels(k)%eps=eps
        Cels(k)%sig_eps=sig_eps
    ! %%%%%%
     ! END of INITIALIZation of the plane
     ! %%%%%%
    
     ep_LJ = eta ! *Cels(1)%H_modulus                    ! %%%  0.1 kpp0 =0.02 pN*um
     sig_LJ = 1.5* sum(Cels(2)%Len_E(1:Cels(2)%N_E))/Cels(2)%N_E
     A_LJ =0.4*ep_LJ* sig_LJ**12;B_LJ =0.4*ep_LJ* sig_LJ**6
    
     k=1;real_t=0.0;tk=1;t=0;    dt=0.005 ;dt0=dt; seed1=26479462
        PP=0.0001/1200
    Integral_H0 = sum(Cels(2)%H_V(1:Cels(2)%N_V)*Cels(2)%area_V(1:Cels(2)%N_V)/3.0)/(Cels(2)%area_tot/4.0/PAI)**(0.5)/8/PAI
    kpp_int_H = 200
    
    k=1
    allocate(KfN1(1:Cels(k)%N_V),alN1(1:Cels(k)%N_V),betN1(1:Cels(k)%N_V))
     k=2
    allocate(KfN2(1:Cels(k)%N_V),alN2(1:Cels(k)%N_V),betN2(1:Cels(k)%N_V))
    
    k= ubound(Cels(1)%Area_F,1)+ubound(Cels(2)%Area_F,1)
    Cout%N_V=Cels(1)%N_V+Cels(2)%N_V
    Cout%N_E=Cels(1)%N_E+Cels(2)%N_E
    
    
            call alloc_SglCel(Cout,k)
            
            
            if(IS_READ.eq.1)then
                    Cels(1)%FILE_Cel = 'RPOT.dat'
                    Cels(2)%FILE_Cel = 'RCOT.dat'
                do k=1,2
                 call read_Cell(Cels(k))
                 call label_Cell(cels,k)
                 call Get_shape_information(Cels,k)  
               enddo  
            endif
            
            
                  call combine_Cells(Cels(1),Cels(2),Cout)
                 
                 call save_Cell(Cout,kk,t,file_input)
             
                                


   
   ! pause
! %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! SEC. II        Dynamics
  Cels(1)%D_thres = 3.5
     ! call Get_isNearIN(Cels,k)

     
             
     active_cells=1
     call Get_NearV_Cells(Cels,1,2, Cels(1)%D_thres)  ! In manipulate_cell_mod
     call Get_NearV_Cells(Cels,2,1, Cels(1)%D_thres)
     call Get_NearF_Cells(Cels,1,2,Cels(1)%D_thres,r00,contf) ! in force_mod.f90
     call Get_shape_information(Cels,1)
     call Get_NearF_Cells(Cels,2,1,Cels(1)%D_thres,r00,contf)
     call Get_shape_information(Cels,2)
     call Get_shape_information(Cels,1)
     
 
!           Strategy: 1. Randomly Select a Vertex
!                     2. Move it using dissipative dynamics
    dtk(1)=0.0025/2
    dtk(2) = dtk(1)
    kk = 0
    kin_cells=0
    is_t1 = 1
     NN= Cels(2)%N_V*(t1+t2)*100
    Z_FIXED = Z_STOP(1)
    
  
    
    loop_out: do kk = 1, 1
    loop_main:do t = 1,  NN !900000001
        
        
        
        k_Z = floor(t*1.0/(t1+t2))+1
        t_extra = mod(t,t1+t2)
        if(t_extra.lt.t1)then
            IS_T1=0
        else
            IS_T1=0
        endif
        if(t_extra.le.t1.and.k_Z.lt.20)then
             Z_FIXED = Z_STOP(k_Z)+(Z_STOP(k_Z+1)-Z_STOP(k_Z))/t1*t_extra
        endif
        
        if(k_Z.ge.20)then
            IS_T1=0
        endif
        IS_FRac=0
        if(t_extra.ge.t1+800)then
            IS_Frac=1
            if(t_extra.eq.t1+800)Frac_mean = 0.0
        endif
        
       do t_NV = 1,Cels(1)%N_V+Cels(2)%N_V
        
        
        
       !%%%% select k
       !%%%
       k=1
       if((Random1(seed1)+0.000001)*(Cels(1)%N_V+Cels(2)%N_V).gt.Cels(1)%N_V)k=2
            
          if(k.eq.1.or.(k.eq.2.and.IS_T1.eq.0) )then
        ! #######################################
        ! SEC. II-1: randomly select a Vertex
          VM%V_move = mod(int( random1(seed1)*(Cels(k)%N_V) ),Cels(k)%N_V)+1
          if(Cels(k)%IS_NEXTBC_V(VM%V_move).eq.0.and.Cels(k)%IS_BC_V(VM%V_move).eq.0)then !%%%% Ensure that k is NOT on the edge
            
            VM%rv_move=Cels(k)%rv(VM%v_move,1:3)
         
        ! #######################################
        ! SEC. II-2: Cal. the Forces
            Cels(k)%P=0
            VM%PM_Force_Point =0.0
            
            diff_int_H=0.0
            VM%PM_Force_Point(1:3)=0
            call Point_H_Force_Cell(Cels,K,VM,diff_int_H)
            VM%PV_Force_Point =0.0
            if(k.eq.2)then
                Cels(k)%lamda = Cels(k)%Kpp_area*(Cels(k)%area_tot-Cels(k)%area_tot_zero)/Cels(k)%area_tot_zero
                call Point_PV_Force_Cell(Cels,K,VM)
             endif   
           VM%Rep_Force_point=0.0
            !call Point_Rep_Force_Cell_linear(Cels,K,ToT_NCel,VM,eta,R01)
            call Update_NearF_Point_Rep(Cels,k,Tot_NCel,VM,r00,r01,eta)
          
           !VM%Rep_Force_point=0.0
        ! #######################################
        ! SEC. II-3: Cal. the Forces   +VM%Rep_Force_Point
    
            if(dtk(k)*sum((VM%PV_Force_Point(1:3)+VM%PM_Force_Point(1:3) &
               +VM%Rep_Force_point )**2)**0.5.lt.0.05.and.Cels(k)%IS_BC_V(VM%v_move).eq.0 )then
                Cels(k)%rv(VM%v_move,1:3)=Cels(k)%rv(VM%v_move,1:3)+dtk(k)*(VM%PV_Force_Point+VM%PM_Force_Point &
                    +VM%Rep_Force_point)*Cels(k)%N_V*1.0/Cels(2)%N_V
           
                
                
            else
                big_move=big_move+1
            endif
           
            VM%rv_move=Cels(k)%rv(VM%v_move,1:3)
     
            
        ! #######################################
        ! SEC. II-5: Update the shape around the Vertex VM    
            call Point_update_cell(Cels,k,VM)  !%% force_mod.f90
            real_t=real_t+dtk(2)/Cels(2)%N_V
        ! #######################################
        ! SEC. II-6: Surf Operations    
            
 
           endif !%%%% END of Ensure that k is NOT on the edge 
            
        endif ! end of IS_T1
    
        
        
       enddo ! End of t_NV
       ! %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        
           ! RC0(1)=0;RC0(2)=0;Rc0(3)=sum(Cels(2)%rv(1:Cels(2)%N_v,3))/Cels(2)%N_V
           ! call Move_Cell(Cels,2,Rc0)
      
        
        
        if(mod(t,20).eq.0)then
                     
                       call Vertex_Averaging_Cell_Origin(Cels(1),C_AVE(1),kFN1,alN1,betN1,0,10)
                       call Vertex_Averaging_Cell_Origin(Cels(2),C_AVE(2),kFN2,alN2,betN2,0,10)
                   
                i=1;j=2
                    ! call Get_NearV_Cells(Cels,i,j, Cels(1)%D_thres)
                    ! call Get_NearV_Cells(Cels,j,i, Cels(1)%D_thres)
                     call Get_NearF_Cells(Cels,i,j,Cels(1)%D_thres,r00,contf)
                     call Get_shape_information(Cels,i)
                     call Get_NearF_Cells(Cels,j,i,Cels(1)%D_thres,r00,contf)
                     call Get_shape_information(Cels,j)
                     if(IS_Frac.eq.1)then
                         frac_mean=frac_mean+(Cels(1)%Contarea+Cels(2)%Contarea)/Cels(2)%area_tot/2.0/10.0
                     endif
             
        endif
        
        
        ! #######################################
        ! SEC. II-7: Print and Save    
        
             if(mod(t,200).eq.1)then
                  !call Get_shape_information(Cels,2)
                 EH=Cels(1)%EnergyH+Cels(2)%EnergyH  
	           write(*,*)file_input,t,kk
               write(*,*)R0,'center:',sum(Cels(2)%RV(1:Cels(2)%N_V,3)*Cels(2)%area_V(1:Cels(2)%N_V)/3)/Cels(2)%area_tot
               write(*,*)Cels(2)%area_tot,sum(Cels(2)%area_V(1:Cels(2)%N_V)/3)
               write(*,*)'curvature energies(soft particle, plane, all)',Cels(1)%EnergyH*16,Cels(2)%EnergyH*16,EH*16
               write(*,*)'contact fraction:',Cels(1)%Contarea/Cels(2)%area_tot,Cels(2)%Contarea/Cels(2)%area_tot
               write(*,*)'adhesion:',Cels(1)%energy_rep,Cels(2)%energy_rep,Cels(2)%energy_rep/Cels(1)%Contarea
               write(*,*)'Z_pos:',sum(Cels(2)%rv(1:Cels(2)%N_v,3))/Cels(2)%N_V
                call combine_Cells(Cels(1),Cels(2),Cout)
                call save_Cell(Cout,kk,k_Z,file_input) ! file_input='pln'
                call save_Cell(Cels(1),kk,K_Z,'POT')
                call save_Cell(Cels(2),kk,k_Z,'COT')
                   
             endif
             if(t_extra.eq.0.and.t.gt.100)then
                 write(*,*)'******t_extra=0'
                 fracALL(k_Z)= frac_mean
                 EALL(1,k_Z)=Cels(1)%EnergyH*16
                 EALL(2,k_Z)=Cels(2)%EnergyH*16
                 open(7,FILE='data/EnergyH5.dat')
                  do i=1,k_Z
                   write(7,'(I5,4(F20.12))')i,Z_STOP(i),fracALL(i),EALL(1,i),EALL(2,i)
                   enddo
                 
             endif
             
                

    enddo loop_main
        
    enddo loop_out ! kk
    

call dealloc_SglCel(Cout) 
  deallocate(kFN1,alN1,betN1)
  deallocate(kFN2,alN2,betN2)

    END






     

    
    
    
    
    
    
    
    
    
    
    




