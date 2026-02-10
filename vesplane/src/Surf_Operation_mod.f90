!  ************************************************
!  Part of 3D vesicle program
!  Flip some bonds that blong to some obtuse triangles
!                    written by LiJF  Dobb@bbs.fudan.edu.cn
! *************************************************  
! include 'readfe.f90'
! include 'label_system.f90'

Module Surf_Operation_mod
   
  use CelType_mod
  use basicfun_mod
  use Manipulate_Cell_mod
   implicit none
    contains
    
    
        subroutine WKONTRI_Surf_scheme_Cell(C,k,k1,I_flip_tot,tt,scheme_type)
        implicit none
        type(cel)::C(:)
        integer,intent(in)::k,k1
        integer,intent(inout)::I_flip_tot
        integer I_flip,i,j,tt,tt1, scheme_type,kk
        real*8 area_temp
        real*8,allocatable::al(:),bet(:)
        integer,allocatable::kF_WK(:)
        integer ave_round
        allocate(al(1:C(k)%N_F),bet(1:C(k)%N_F),kF_WK(1:C(k)%N_F) )
        
                call Refine_Cell(C,k,k1,1,1)
                call Get_shape_information(C,k1)
                call WKONTRI_alloc(C(k1))
                call WKONTRI_init(c(k1))
                call WKONTRI_init_albet(C,k,k1,al,bet,kf_WK)
                
                
                ave_round=5
                tt1=3
                if(scheme_type.eq.2)then
                    ave_round=3
                    tt1=1
                endif
                area_temp = C(k)%area_tot
                c(k)%detL=sum( ( c(k)%Len_E(1:C(k)%N_E)- c(k)%L0)**2/2.0 )/c(k)%N_E
                
                kk=0
                do while(kk.le.500.and.c(k)%detL.gt.0.003)
                    c(k)%L0 = sum( c(k)%Len_E(1:C(k)%N_E))/c(k)%N_E
                   c(k)%detL=sum( ( c(k)%Len_E(1:C(k)%N_E)- c(k)%L0)**2/2.0 )/c(k)%N_E
                    kk=kk+1
                    I_flip=0;
                    do j=1,tt1
                    call Bond_flip_Cell(C,k,I_flip)
                    I_flip_tot=I_flip_tot+I_flip
                    call label_Cell(C,k)
                    call Get_Cell_area(C,k)
                    call Get_shape_information(C,k)
                    enddo
                    if(scheme_type.eq.1)call WKONTRI_Average_Len_CELL(C,k,k1,kf_WK,al,bet,Ave_round)
                    if(scheme_type.eq.2)call WKONTRI_Average_CELL(C,k,k1,kf_WK,al,bet,Ave_round)
                     
                    call Get_Cell_area(C,k);call Get_Cell_Vol(C,k);call Get_Cell_Len(C,k)
                enddo
                
                C(k)%Rv(1:C(k)%N_V,:) = C(k)%Rv(1:C(k)%N_V,:) * (area_temp/ C(k)%area_tot)**(0.5) 
                 call Get_shape_information(C,k)
                
                call WKONTRI_dealloc(C(k1))
        end subroutine WKONTRI_Surf_scheme_Cell
   
    
            
    
    subroutine WKONTRI_init_albet(C,k,k1,al,bet,kf_WK)
        implicit none
        type(Cel)::C(:)
        real*8 v(1:3)
        real*8 al(:),bet(:),al1,bet1
        integer k,k1,kf_WK(:)
        integer i,j
        do i=1,c(k)%N_V
            j = c(k1)%F_V(i,1)
            kf_WK(i) =j
            v(1:3) = c(k)%Rv(i,1:3) - c(k1)%RF(j,1,1:3)
            call WKONTRI_Convert(C(k1),j,v,al1,bet1)
            al(i) = al1
            bet(i) =bet1
        enddo
        
        
    end subroutine WKONTRI_Init_albet
    
    subroutine WKONTRI_Average_CELL(C,k,k1,kf_WK,al,bet,Ave_round)
        implicit none
        type(CEL)::C(:)
        integer k,k1
        integer i,j,kf_wk(:)
        integer k_F0,k_F1
        integer T_ave,ave_round
        
        real*8 al(:),bet(:),v1(1:3),v2(1:3), rcen(1:3)
        real*8 r0(1:3)
        real*8 al0,bet0,al1,bet1
        
        
        do T_ave = 1,ave_round
            do i=1,c(k)%N_V
                
                rcen=0
                r0=c(k)%rv(i,1:3)
                do j=1,c(k)%N_V_V(i)
                    rcen = rcen + c(k)%rv( c(k)%v_v(i,j),1:3)
                enddo
                rcen=rcen/c(k)%N_V_V(i)
                
                v1 = rcen-r0
                r0=c(k)%rv(i,1:3)
                K_F0= kf_wk(i)
                al0=al(i); bet0=bet(i)
                v2=r0
               if(sum(v1**2).gt.0.00000001)then
                call WKONTRI_MOVE(C(k1),v1,k_F0,al0,bet0,r0,k_F1,al1,bet1,V2,C(k1)%IS_BC_F)
                kf_wk(i) = k_F1
                al(i) =al1
                bet(i) =bet1
                c(k)%rv(i,1:3) = v2
               endif                
            enddo
        enddo
        
        
        
    
    
    
    END subroutine
    
    
    subroutine WKONTRI_Average_Len_CELL(C,k,k1,kf_WK,al,bet,Ave_round)
        implicit none
        type(CEL)::C(:)
        integer k,k1
        integer i,j,kf_wk(:)
        integer k_F0,k_F1
        integer T_ave,ave_round
        
        real*8 al(:),bet(:),v1(1:3),v2(1:3), rcen(1:3)
        real*8 r0(1:3)
        real*8 al0,bet0,al1,bet1
        real*8 dt,L0, detL
        real*8,allocatable::dFdr(:,:)
        integer i_E,i_V
        allocate(dFdr(1:c(k)%N_V,1:3))
        
        L0 = c(k)%L0
        dt = 0.1
        
        do T_ave = 1,ave_round
            call Get_Len_E(C(k)%rv,  c(k)%V_E, c(k)%N_E, c(k)%Len_E, c(k)%Vector_E)
            if(mod(T_ave,5).eq.0)then
                !detL=sum( ( c(k)%Len_E(1:C(k)%N_E)-L0)**2/2.0 )/c(k)%N_E   
                !write(*,*)detL
            endif
            do i=1,c(k)%N_V
                r0=c(k)%rv(i,1:3)
                dfdr(i,:)=0.0
                do j=1,c(k)%N_V_V(i)
                    i_e = abs(C(k)%E_V(i,j))
                    i_V = C(k)%V_V(i,j)
                    dfdr(i,:) = dfdr(i,:) + (c(k)%Len_E(i_e) - L0) / c(k)%Len_E(i_E) * (r0 - c(k)%rv(i_V,1:3))
                enddo
            enddo
          ! calculate the force
            do i=1,c(k)%N_V
                v1 = -dfdr(i,:)*dt
                r0=c(k)%rv(i,1:3)
                K_F0= kf_wk(i)
                al0=al(i); bet0=bet(i)
                v2=r0
               if(sum(v1**2).gt.0.00000001)then
                call WKONTRI_MOVE(C(k1),v1,k_F0,al0,bet0,r0,k_F1,al1,bet1,V2,C(k1)%IS_BC_F)
                kf_wk(i) = k_F1
                al(i) =al1
                bet(i) =bet1
                c(k)%rv(i,1:3) = v2
               endif
            enddo
            
        enddo !end of round
        
        
       deallocate(dFdr) 
    
    
    
    END subroutine
    
    subroutine Surf_scheme_Cell(C,k,k1,I_flip_tot)
        implicit none
        type(cel)::C(:)
        integer,intent(in)::k,k1
        integer,intent(inout)::I_flip_tot
        integer I_flip,i
        real*8 area_temp
      
                area_temp = C(k)%area_tot
                !call vertex_averaging_fixshape_refine_Cell(C,k,k1,2,20)
                do i=1,30
                call Vertex_averaging_Cell(C,k)
                enddo
                call Get_Cell_area(C,k);call Get_Cell_Vol(C,k);call Get_Cell_Len(C,k)
                I_flip=0;
                call Bond_flip_Cell(C,k,I_flip)
                I_flip_tot=I_flip_tot+I_flip
                if(I_flip.ge.1)then
                    write(*,*)'bond_flip  ','total flips'
                    write(*,*)I_flip,I_flip_tot
                endif
                call label_Cell(C,k)
                call Get_Cell_area(C,k)
                C(k)%rv(1:C(k)%N_V,1:3) = C(k)%rv(1:C(k)%N_V,1:3) * (area_temp/C(k)%area_tot)**0.5
                call Get_shape_information(C,k)
    end subroutine Surf_scheme_Cell
    
    subroutine label_cell(C,k)
        implicit none
        type(cel)::C(:)
        integer,intent(in)::k
        integer i,j,i1,i2,i3,i4,j1,j2
        integer ksw
        call label_system(C(k)%V_F,C(k)%E_V,C(k)%N_E_V,C(k)%F_V,C(k)%N_F_V,&
        C(k)%F_E,C(k)%N_F_E,C(k)%V_V,C(k)%N_V_V,C(k)%E_V_opposite,C(k)%N_E_V_opposite, &
        C(k)%V_E_opposite,C(k)%N_V,C(k)%N_E,C(k)%N_F,C(k)%V_E,C(k)%E_F,C(k)%F_F,&
         C(k)%IS_OPEN,C(k)%IS_BC_E,C(k)%IS_BC_V,C(k)%IS_NEXTBC_V,C(k)%N_BC_E,C(k)%N_BC_V,C(k)%E_BC,C(k)%V_BC,C(k)%IS_BC_F) 
        
        do i=1,C(k)%N_F
            i4=0
            do i1=1,3
                i2 = C(k)%V_F(i,i1)
                do i3 =1,C(k)%N_V_V(i2)
                    j1 = C(k)%F_V(i2,i3)
                    ksw = 0
                    do j2 = 1, i4
                        if(j1.eq.i)ksw=1
                        if(j1.eq.C(k)%F_F_around(i,i4) ) ksw =1
                    enddo
                    if(ksw.eq.0)then
                        i4=i4+1
                        C(k)%N_F_F_around(i)=i4
                        C(k)%F_F_around(i,i4)= j1
                    endif
                    
                enddo
            enddo
            
        enddo
        
        
        
    end subroutine label_cell   
    
    
    subroutine label_Sgl_cell(C)
        implicit none
        type(cel)::C
        integer i,j,i1,i2,i3,i4,j1,j2
        integer ksw
        
        call label_system(C%V_F,C%E_V,C%N_E_V,C%F_V,C%N_F_V,&
        C%F_E,C%N_F_E,C%V_V,C%N_V_V,C%E_V_opposite,C%N_E_V_opposite, &
        C%V_E_opposite,C%N_V,C%N_E,C%N_F,C%V_E,C%E_F,C%F_F,&
         C%IS_OPEN,C%IS_BC_E,C%IS_BC_V,C%IS_NEXTBC_V,C%N_BC_E,C%N_BC_V,C%E_BC,C%V_BC,C%IS_BC_F) 
        
        
        do i=1,C%N_F
            i4=0
            do i1=1,3
                i2 = C%V_F(i,i1)
                do i3 =1,C%N_V_V(i2)
                    j1 = C%F_V(i2,i3)
                    ksw = 0
                    do j2 = 1, i4
                        if(j1.eq.i)ksw=1
                        if(j1.eq.C%F_F_around(i,i4) ) ksw =1
                    enddo
                    if(ksw.eq.0)then
                        i4=i4+1
                        C%N_F_F_around(i)=i4
                        C%F_F_around(i,i4)= j1
                    endif
                    
                enddo
            enddo
            
        enddo
        
        
    end subroutine label_Sgl_cell  
    
    subroutine Refine_Cell(C,k1,k2,n,method)
        implicit none
        type(Cel)::C(:)
        integer,intent(in)::k1,k2,n,method
        integer k,i

       C(k2)%Rv(1:C(k1)%N_V,1:3)= C(k1)%Rv(1:C(k1)%N_V,1:3)
       C(k2)%N_V=C(k1)%N_V;C(k2)%N_E=C(k1)%N_E;C(k2)%N_F=C(k1)%N_F;
       C(k2)%V_E(1:C(k1)%N_E,1:2)=C(k1)%V_E(1:C(k1)%N_E,1:2)
       C(k2)%E_F(1:C(k1)%N_F,1:3)=C(k1)%E_F(1:C(k1)%N_F,1:3)
       do i=1,n
          if(method.eq.1)then
              call Refine_surf_butterfly(C(k2)%Rv,C(k2)%N_V,C(k2)%N_E,C(k2)%N_F,C(k2)%V_E,C(k2)%E_F)
          endif
       enddo
       call label_Cell(C,k2)
    end subroutine Refine_Cell
     
    
    subroutine bond_flip_Cell(C,k,I_flip)
        implicit none
        type(Cel)::C(:)
        integer,intent(in)::k
        integer,intent(inout)::I_flip
               call  bond_flip_neck(C(k)%V_F,C(k)%E_V,C(k)%N_E_V,C(k)%F_V,C(k)%N_F_V,&
        C(k)%V_V,C(k)%N_V_V,C(k)%F_E,C(k)%N_F_E,C(k)%V_E_opposite,C(k)%N_V,C(k)%N_E,C(k)%N_F,C(k)%V_E,&
        C(k)%E_F,C(k)%rv,C(k)%Len_E,  C(k)%Vector_E,C(k)%IS_BC_E, I_flip)
               
    end subroutine Bond_flip_Cell
    
    subroutine bond_flip_Sgl_Cell(C)
        implicit none
        type(Cel)::C
        integer I_flip
               call  bond_flip_neck(C%V_F,C%E_V,C%N_E_V,C%F_V,C%N_F_V,&
        C%V_V,C%N_V_V,C%F_E,C%N_F_E,C%V_E_opposite,C%N_V,C%N_E,C%N_F,C%V_E,&
        C%E_F,C%rv,C%Len_E,  C%Vector_E,C%IS_BC_E, I_flip)
               
    end subroutine Bond_flip_Sgl_Cell
    
    
    
    subroutine Vertex_averaging_Cell(C,k)
        implicit none
        type(Cel)::C(:)
        integer,intent(in)::k
        
        call Vertex_averaging(C(k)%rv,C(k)%N_V,C(k)%N_F,C(k)%F_V,C(k)%V_F, &
        C(k)%V_V,C(k)%N_F_V ,C(k)%area_F,C(k)%Vector_F,C(k)%Norm_F, &
        C(k)%area_F_zero,C(k)%area_tot,C(k)%area_tot_zero)
      
    end subroutine Vertex_averaging_Cell
    
    
    
    
    Subroutine Vertex_Averaging_Cell_Origin(C1,C2,kFN,alN,betN,typeVA,t_run)
       Implicit None
       type(CEL)::C1,C2
       real*8 Len0,lam,var_Len,eff
       real*8 alN(:),betN(:),al1,al2,al3,al0,bet0,bet1,bet2,bet3
       integer  kFN(:),kf0,kf1,kf2,kf3
       
       real*8 r0(1:3),r1(1:3),r2(1:3),r3(1:3)
       real*8 V0(1:3),V1(1:3),V2(1:3)
       real*8 n1(1:3),n2(1:3),n3(1:3),dd1
       real*8,allocatable::forceL(:,:),norm_V(:,:)
       integer i,j,k,i1,i2,i3,j1,j2,j3
       integer typeVA
       integer t_run,t,I_flip
       real*8 Sums
       
       
        if(typeVA.eq.0)then 
            call copy_Cell(C1,C2)
            
            call label_Sgl_cell(C2)
            call Get_shape_sgl_information(C2)
            do i=1,C1%N_V
                
               if(C1%IS_BC_V(i).eq.0)then 
                kfN(i) = C1%F_V(i,1)
                if(C1%V_F(kfN(i),1).eq.i)then
                alN(i) = 0.0001
                betN(i) = 0.0001
                endif
                if(C1%V_F(kfN(i),2).eq.i)then
                alN(i) = 1-0.0002
                betN(i) = 0.0001
                endif
                if(C1%V_F(kfN(i),3).eq.i)then
                alN(i) = 0.0001
                betN(i) = 1-0.0002
                endif
               endif 
                
                
            enddo
           
       endif
       
        call WKONTRI_alloc(C1)
        call WKONTRI_init(C1)
        
       allocate(norm_V(1:C1%N_V,1:3))
        
        do i=1,C1%N_V
            norm_V(i,1:3) = 0
             if(C1%IS_BC_V(i).eq.0)then
            do j=1,C1%N_V_V(i)
                i1 = C1%F_V(i,j)
                norm_V(i,1:3) =norm_V(i,1:3)+C1%Norm_F(i1,1:3)*C1%area_F(i1)     
            enddo
            norm_V(i,1:3)=norm_V(i,1:3)/sum(norm_V(i,1:3)**2)**0.5
            endif
        enddo
        
     ! %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%   
     ! AVERAGING
     ! %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  
        do t= 1,t_run
        
        do i=1,C2%N_V
            R1=0
            SUMS = 0
            if(C1%IS_BC_V(i).eq.0.and.C1%IS_NEXTBC_V(i).eq.0)then
            do j=1,C2%N_V_V(i)
                i1 = C2%V_V(i,j)
                R1=R1 + C2%RV(i1,1:3)*C2%area_V(i1)
                SumS = SumS + C2%area_V(i1)
            enddo
                R1 = R1 /SumS
                V1 = R1 - C2%RV(i,1:3)
                kf1 = KFN(i);al1 = alN(i);bet1 = betN(i);
                n1 =norm_V(C1%V_F(kf1,1),1:3);n2 =norm_V(C1%V_F(kf1,2),1:3);n3 =norm_V(C1%V_F(kf1,3),1:3)
                v2 = n1*(1-al1-bet1)+n2*al1+n3*bet1
                V2 = V2/sum(V2**2)**0.5
                V1 = V1  - sum(v1*V2)*v2
                r1 = C1%RF(kf1,1,1:3) + al1*C1%e(kf1,1,1:3)+ bet1*C1%e(kf1,2,1:3)
                
               
                
               
                call WKONTRI_MOVE(C1,V1,kf1,al1,bet1,R1,kf2,al2,bet2,R2,C1%IS_BC_F)
                KFN(i)=kf2
                alN(i)=al2
                betN(i)=bet2
                
               n1 = c1%nF(kf2,1,:);n2 = c1%nF(kf2,2,:);n3 = c1%nF(kf2,3,:)
               call WKONTRI_Interp_Surf(n1,n2,n3,C1%lenF(kf2,1),C1%lenF(kf2,2),C1%lenF(kf2,3),&
                      C1%thF(kf2,1),C1%thF(kf2,2),C1%thF(kf2,3),al2,bet2,R3)
               C2%RV(i,1:3) = R2 +R3
            endif
        enddo
        
               call bond_flip_Sgl_Cell(C2)
               call label_Sgl_cell(C2)
               call Get_shape_sgl_information(C2)
          
               
              ! call bond_flip_Sgl_Cell(C2)
              ! call label_Sgl_cell(C2)
              ! call Get_shape_sgl_information(C2)
        
        enddo
     
     ! %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%   
     ! AVERAGING
     ! %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%   
        
      call WKONTRI_dealloc(C1)
        
       
       if(typeVA.eq.0)then 
            call copy_Cell(C2,C1)
        
       endif
       deallocate(norm_V)
        
    END Subroutine Vertex_Averaging_Cell_Origin
    
    
    
    
    
    Subroutine Vertex_Averaging_Cell_LEN(C1,C2,kFN,alN,betN,typeVA,t_run)
       Implicit None
       type(CEL)::C1,C2
       real*8 Len0,lam,var_Len,eff
       real*8 alN(:),betN(:),al1,al2,al3,al0,bet0,bet1,bet2,bet3
       integer  kFN(:),kf0,kf1,kf2,kf3
       
       real*8 r0(1:3),r1(1:3),r2(1:3),r3(1:3)
       real*8 V0(1:3),V1(1:3),V2(1:3)
       real*8 n1(1:3),n2(1:3),n3(1:3),dd1
       real*8,allocatable::forceL(:,:),norm_V(:,:)
       integer i,j,k,i1,i2,i3,j1,j2,j3
       integer typeVA
       integer t_run,t,I_flip
       t_run = 1000
       lam = 0.01
       eff=0.001
       allocate(forceL(1:C1%N,1:3))
       
       if(typeVA.eq.0)then 
            call copy_Cell(C1,C2)
            
               call label_Sgl_cell(C2)
               call Get_shape_sgl_information(C2)
            do i=1,C1%N_V
                kfN(i) = C1%F_V(i,1)
                if(C1%V_F(kfN(i),1).eq.i)then
                alN(i) = 0.0001
                betN(i) = 0.0001
                endif
                if(C1%V_F(kfN(i),2).eq.i)then
                alN(i) = 1-0.0002
                betN(i) = 0.0001
                endif
                if(C1%V_F(kfN(i),3).eq.i)then
                alN(i) = 0.0001
                betN(i) = 1-0.0002
                endif
                
                
            enddo
           
       endif
       
        call WKONTRI_alloc(C1)
        call WKONTRI_init(C1)
        
       allocate(norm_V(1:C1%N_V,1:3))
        
        do i=1,C1%N_V
            norm_V(i,1:3) = 0
            do j=1,C1%N_V_V(i)
                i1 = C1%F_V(i,j)
                norm_V(i,1:3) =norm_V(i,1:3)+C1%Norm_F(i1,1:3)*C1%area_F(i1)     
            enddo
            norm_V(i,1:3)=norm_V(i,1:3)/sum(norm_V(i,1:3)**2)**0.5
        enddo
            len0 = sum(C2%len_E(1:C2%N_E))/C2%N_E
            var_len = sum( (C2%len_E(1:C2%N_E)-len0)**2 )/C2%N_E
       
       do t=1,t_run
           
            do i=1,C2%N_V
                forceL(i,1:3) = 0.0
                do j=1,C2%N_V_V(i)
                        i1 = C2%V_V(i,j)
                        R2 = C2%RV(i,1:3)-C2%RV(i1,1:3)
                        dd1 = sum(R2**2)**0.5
                        forceL(i,:)=forceL(i,:) - R2 *(dd1-Len0) / dd1*lam
                enddo
                do j=1,C2%N_V
                    if(i.ne.j)then
                        R2 = C2%RV(i,1:3)-C2%RV(j,1:3)
                        dd1 = sum(R2**2)**0.5
                        forceL(i,:)=forceL(i,:) + R2 /dd1**3*eff
                     
                    endif
                enddo
                
            enddo
            
            do i=1,C2%N_V
                v1 = forceL(i,:)
                kf1 = KFN(i);al1 = alN(i);bet1 = betN(i);
                n1 =norm_V(C1%V_F(kf1,1),1:3);n2 =norm_V(C1%V_F(kf1,2),1:3);n3 =norm_V(C1%V_F(kf1,3),1:3)
                v2 = n1*(1-al1-bet1)+n2*al1+n3*bet1
                V2 = V2/sum(V2**2)**0.5
                V1 = V1  - sum(v1*V2)*v2
               r1 = C1%RF(kf1,1,1:3) + al1*C1%e(kf1,1,1:3)+ bet1*C1%e(kf1,2,1:3)
               
                call WKONTRI_MOVE(C1,V1,kf1,al1,bet1,R1,kf2,al2,bet2,R2,C1%IS_BC_F)
                KFN(i)=kf2
                alN(i)=al2
                betN(i)=bet2
                
               n1 = c1%nF(kf2,1,:);n2 = c1%nF(kf2,2,:);n3 = c1%nF(kf2,3,:)
               call WKONTRI_Interp_Surf(n1,n2,n3,C1%lenF(kf2,1),C1%lenF(kf2,2),C1%lenF(kf2,3),&
                      C1%thF(kf2,1),C1%thF(kf2,2),C1%thF(kf2,3),al2,bet2,R3)
               C2%RV(i,1:3) = R2 +R3
                
            enddo
             
           
           if(mod(t,50).eq.1)then
               call bond_flip_Sgl_Cell(C2)
               call label_Sgl_cell(C2)
               call Get_shape_sgl_information(C2)
                 len0 = sum(C2%len_E(1:C2%N_E))/C2%N_E
                 var_len = sum( (C2%len_E(1:C2%N_E)-len0)**2 )/C2%N_E
                write(*,*)t,len0,var_len
                
          endif
           
           
           
           
           
           
           
       enddo ! end of T_run
       
       
        call WKONTRI_dealloc(C1)
        
       
       if(typeVA.eq.0)then 
            call copy_Cell(C2,C1)
            call dealloc_SglCel(C2)
       endif
       
    deallocate(ForceL)
    deallocate(Norm_V)
    
    END Subroutine Vertex_Averaging_Cell_LEN
    
    
    
    subroutine vertex_averaging_fixshape_Cell(C,k,k1,n)
        implicit none
        type(Cel)::C(:)
        integer,intent(in)::k,k1,n
        integer,allocatable::vvNear(:)
        integer i,j,nn
        nn=ubound(C(k)%rv,1)
        allocate(VVNear(1:Nn))
        do i = 1,c(k)%N_V
            vvNear(i) = i    
        enddo
        C(k1)%rv(1:c(k)%N_V,:) = c(k)%rv(1:c(k)%N_V,:)
        C(k1)%norm_F(1:c(k)%N_F,:) = c(k)%norm_F(1:c(k)%N_F,:)
        do i = 1,n
        call Vertex_averaging_fixshape(C(k)%rv,c(k1)%rv,C(k)%N_V,C(k)%N_F,C(k)%F_V,C(k)%V_F, &
        C(k)%V_V,C(k)%N_F_V,vvNear ,C(k)%area_F,C(k)%Vector_F,C(k)%Norm_F,c(k1)%Norm_F, &
        C(k)%area_F_zero,C(k)%area_tot,C(k)%area_tot_zero)
        call Get_Cell_area(C,k)
        enddo
        deallocate(vvNear)
    end    subroutine vertex_averaging_fixshape_Cell 
    
    

    subroutine vertex_averaging_fixshape_refine_Cell(C,k,k1,n_rf,n)
        implicit none
        type(Cel)::C(:)
        integer,intent(in)::k,k1,n,n_rf
        integer,allocatable::vvNear(:)
        integer i,j,nn
        nn=ubound(C(k)%rv,1)
        allocate(VVNear(1:Nn))
        do i = 1,c(k)%N_V
            vvNear(i) = i    
        enddo
        call Refine_Cell(C,k,k1,n_rf,1)
        call label_cell(C,k1)
        call get_Cell_area(C,k1)
        do i = 1,n
        call Vertex_averaging_fixshape_refine(C(k)%rv,c(k1)%rv,C(k)%N_V,C(k)%N_F,C(k)%F_V,C(k)%V_F, &
        C(k)%V_V,C(k)%N_F_V,C(k1)%F_V,C(k1)%V_F, &
        C(k1)%V_V,C(k1)%N_F_V ,vvNear,C(k)%area_F,C(k)%Vector_F,C(k)%Norm_F,c(k1)%Norm_F, &
        C(k)%area_F_zero,C(k)%area_tot,C(k)%area_tot_zero)
        call Get_Cell_area(C,k)
        enddo
         deallocate(vvNear)
    end    subroutine vertex_averaging_fixshape_refine_Cell 
        
    
    subroutine Mesh_grow_Cell(C,k,t,seed1)
        implicit none
        type(Cel)::C(:)
        integer,intent(in)::k
        integer num_F,seed1
        integer*4 t
        
        if(mod(t,500000).eq.1.and.C(k)%N_E.lt.C(k)%N*2)then
            num_F = mod(int( random1(seed1)*(C(k)%N_F) ),C(k)%N_F)+1
     
        call Insert_APOINT(C(k)%RV,C(k)%Vector_F,C(k)%area_F,C(k)%N_V,C(k)%N_E, &
        C(k)%N_F,C(k)%V_E,C(k)%E_F,C(k)%V_F,num_F)
      
        call label_cell(C,k)
        call Get_shape_information(C,k)
        C(k)%area_F_zero=C(k)%area_F;C(k)%area_tot_zero = C(k)%area_tot
       ! call Get_IsNearIn(C,k)
       ! call Get_NearV_Cells(C,1,2,C(1)%D_thres)  ! in the module_celltype.f90, manipulate_cell
       ! call Get_NearV_Cells(C,1,1,C(1)%D_thres)
        
        endif
    end subroutine mesh_grow_cell
    

   subroutine label_system(V_F,E_V,N_E_V,F_V,N_F_V,F_E,N_F_E,V_V,N_V_V,  &
 & E_V_opposite,N_E_V_opposite, V_E_opposite,N_V,N_E,N_F,V_E,E_F,F_F,&
IS_OPEN,IS_BC_E,IS_BC_V,IS_NEXTBC_V,N_BC_E,N_BC_V,E_BC,V_BC,IS_BC_F)
implicit none

integer N_V,N_E,N_F
integer V_E(:, :)
integer E_F(:, :)
! %%%%%%%%%%%%% other labels
! the first kind label
integer V_F(:,:)  ! vertices around a face
integer E_V(:,:), N_E_V(:) ! edges ajecent to a vertice
integer F_V(:,:), N_F_V(:) ! star faces around a vertice
integer F_E(:,:), N_F_E(:)
integer F_F(:,:)

! second kind label
integer V_V(:,:), N_V_V(:)
integer E_V_opposite(:,:), N_E_V_opposite(:)
integer V_E_opposite(:,:)

integer i,j,k,i1,i2,i3,i4,i5,i0,I_temp(1:40),I_temp2(1:40),I_temp3(1:40)

! shape's BC
integer IS_OPEN
integer IS_BC_E(:),IS_BC_V(:),IS_NEXTBC_V(:)
integer N_BC_E,N_BC_V
integer E_BC(:),V_BC(:),IS_BC_F(:)
integer sign_BC




IS_OPEN = 0
IS_BC_E=0
IS_BC_V=0
IS_NEXTBC_V=0
N_BC_E=0; N_BC_V=0
E_BC=0;V_BC=0
! the first kind label
! (1) V_F
  do i=1,N_F
  V_F(i,1)=V_E(abs(E_F(i,1)), int( (3-E_F(i,1)/abs(E_F(i,1)) ) /2 ) )
  V_F(i,2)=V_E(abs(E_F(i,2)), int( (3-E_F(i,2)/abs(E_F(i,2)) ) /2 ) )
  V_F(i,3)=V_E(abs(E_F(i,3)), int( (3-E_F(i,3)/abs(E_F(i,3)) ) /2 ) )
  enddo

! (2) E_V
!       Problem: E_V û�Ա����������?
	N_E_V=0
  do i=1, N_E
  i1=V_E(i,1); i2=V_E(i,2)
  N_E_V(i1)=N_E_V(i1)+1
  N_E_V(i2)=N_E_V(i2)+1
  E_V(i1, N_E_V(i1))=i
  E_V(i2, N_E_V(i2))=-i
  enddo
! (3) F_V �?����õģ�����
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
F_E = 0
do i=1, N_F
  i1=E_F(i,1);i2=E_F(i,2);i3=E_F(i,3)
  N_F_E(abs(i1))=N_F_E(abs(i1))+1
  N_F_E(abs(i2))=N_F_E(abs(i2))+1
  N_F_E(abs(i3))=N_F_E(abs(i3))+1
  F_E(abs(i1), int( (3-i1/abs(i1) ) /2 ) )=i
  F_E(abs(i2), int( (3-i2/abs(i2) ) /2 ) )=i
  F_E(abs(i3), int( (3-i3/abs(i3) ) /2 ) )=i
enddo

!(4-1) Determine its BC
IS_BC_V=0;IS_BC_E=0
do i=1,N_E
    if(F_E(i,1).eq.0.or.F_E(i,2).eq.0)then
        IS_OPEN=1
        Sign_BC = 1
        if(F_E(i,1).eq.0)Sign_BC=-1
        N_BC_E=N_BC_E+1
        E_BC(N_BC_E) =Sign_BC*i
        IS_BC_E(i) = Sign_BC
        i1=V_E(i,1);i2=V_E(i,2)
        if(IS_BC_V(i1).ne.1)then
            N_BC_V=N_BC_V+1
            IS_BC_V(i1)=1
            V_BC(N_BC_V) = i1
        endif
        if(IS_BC_V(i2).ne.1)then
            N_BC_V=N_BC_V+1
            IS_BC_V(i2)=1
            V_BC(N_BC_V) = i2
        endif        
    endif    
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
! On BC F_F might be 0
 

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
   
   !(1) F_V; sort the F_V and E_V_opposite
   do i=1, N_V
   if(IS_BC_V(i).eq.0)then    ! If V is on BC, then there is no need to sort F_V
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
	    F_V(i,j) = I_temp3( I_temp(j) )  ! F_V ���?
        enddo
  endif      
  enddo
!  (2) E_V & (5) V_V ��

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

! (7) V_E_opposite  use F_E 4�?
    V_E_opposite=0
   do i=1,N_E
   if(IS_BC_E(i).eq.0)then ! only if E is not on BC, we have to cal V_E
   do j=1,2
    i1=F_E(i,j);
	i2=V_E(i,1)
	i3=V_E(i,2)
	do k=1,3
	 i4=V_F(i1,k)
	 if(i4.ne.i2.and.i4.ne.i3)V_E_opposite(i,j)=i4
	enddo
   enddo
   endif
   enddo

   IS_NEXTBC_V=IS_BC_V
   do i=1,N_V
       do j=1,N_V_V(i)
           if(IS_BC_V(V_V(i,j)).eq.1)then
               IS_NEXTBC_V(i)=1
           endif
       enddo
       
   enddo
   IS_BC_E=0
   do i=1,N_E
       i1=V_E(i,1);i2=V_E(i,2)
       if(IS_NEXTBC_V(i1).eq.1.and.IS_NEXTBC_V(i2).eq.1)then
           IS_BC_E(i)=1
       endif
   enddo
   
   IS_BC_F=0
   do i=1,N_F
       i1=V_F(i,1);i2=V_F(i,2);i3=V_F(i,3)
       if(IS_BC_V(i1)+IS_BC_V(i2)+IS_BC_V(i3).ge.2)then
           IS_BC_F(i)=1
       endif
   enddo


END subroutine label_system
  
   subroutine bond_flip_neck(V_F,   E_V,N_E_V,     F_V,N_F_V,V_V,N_V_V,    F_E,N_F_E,   &
       V_E_opposite,   N_V,N_E,N_F,V_E,E_F,rv,   Len_E,  Vector_E,IS_BC_E,&
       I_flip)


    integer N_V,N_E,N_F
    real*8 rV(:,:),x,y,z,a
    real*8 a1,a2,b,c1,c2,switchR
    integer V_E(:, :)
    integer E_F(:, :)
    ! %%%%%%%%%%%%% other labels
    ! the first kind label
    integer V_F(:,:)  ! vertices around a face
    integer E_V(:,:), N_E_V(:) ! edges ajecent to a vertice
    integer F_V(:,:), N_F_V(:) ! star faces around a vertice
    integer F_E(:,:), N_F_E(:)
    integer V_V(:,:), N_V_V(:)
    ! second kind label

    integer V_E_opposite(:,:)
    integer IS_BC_E(:)



    integer i,j,k,i1,i2,i3,i4,i5,i6,i0,I_temp(1:20),I_temp2(1:20),I_temp3(1:20),  I_flip




    ! about Flip
    integer,allocatable:: Permission_E(:)   !  Edges that can be Flipped
    real*8 Len_a1,Len_a2,Len_b, Len_c1, Len_c2  !  length of sides of two triangles
    integer I_a(1:2),I_b, I_c(1:2)
    integer  I_a_n(1:2), I_c_n(1:2)

    real*8 costh1,costh2, Switch_R
    ! Be careful when defining variable with letters "E" "F" "V"
    real*8 Len_E(:)  !length of edges
    real*8 Vector_E(:,:)  !Vector of edges


    integer neck_check   ! only for checking the neck
    integer N
    integer k1,k2,k3,k4,k5,k6,VV(1:100)
    
    
    N = ubound(Len_E,1)
    allocate(Permission_E(1:N))
    Permission_E=1
    do i=1,N_E
        if(abs(IS_BC_E(i)).eq.1)Permission_E(i)=0
    enddo
    
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
    if(SwitchR.lt.-0.05)then
    !if(SwitchR.lt.-0.05)then


    neck_check=0
    i5=V_E_opposite(i,1)
    i6=V_E_opposite(i,2)
    do j=1,N_V_V(i5)
        if(V_V(i5,j).eq.i6)neck_check=1
    enddo
   
   

    if(neck_check.ne.1)then
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

    i3=V_E(i_b,1)
    i4=V_E(i_b,2)
    V_E(i_b,1)=V_E_opposite(i,1)  !  Edge label Changing
    V_E(i_b,2)=V_E_opposite(i,2)
    
    ! update V_V,N_V_V
    ! i3, i3-i4 is the edge before flipping
        do k1=1,N_V_V(i3)
            if(V_V(i3,k1).eq.i4)k2= k1
        enddo
        VV(1:N_V_V(i3)) = V_V(i3,1:N_V_V(i3))
        V_V(i3,k2:N_V_V(i3)-1) =VV(k2+1:N_V_V(i3))
        V_V(i3,N_V_V(i3))=0
        N_V_V(i3)=N_V_V(i3)-1
    ! i4,i3-i4 is the edge before flipping
        do k1=1,N_V_V(i4)
            if(V_V(i4,k1).eq.i3)k2= k1
        enddo
        VV(1:N_V_V(i4)) = V_V(i4,1:N_V_V(i4))
        V_V(i4,k2:N_V_V(i4)-1) =VV(k2+1:N_V_V(i4))
        V_V(i4,N_V_V(i4))=0
        N_V_V(i4)=N_V_V(i4)-1
    ! i5, i5-i6 the newly formed edge
        do k1=1,N_V_V(i5)
            if(V_V(i5,k1).eq.i3)k2= k1
            if(V_V(i5,k1).eq.i4)k3= k1
        enddo
        if(k3.gt.k2)then
            k5=k2;k6=k3
        else
            k5=k3;k6=k2
        endif
        VV(1:N_V_V(i5)) = V_V(i5,1:N_V_V(i5))
        V_V(i5,k5+1) =i6
        V_V(i5,k6+1:N_V_V(i5)+1)=vv(k6:N_V_V(i5))
        N_V_V(i5)=N_V_V(i5)+1
    ! i6, i5-i6 the newly formed edge
        do k1=1,N_V_V(i6)
            if(V_V(i6,k1).eq.i3)k2= k1
            if(V_V(i6,k1).eq.i4)k3= k1
        enddo
        if(k3.gt.k2)then
            k5=k2;k6=k3
        else
            k5=k3;k6=k2
        endif
        VV(1:N_V_V(i6)) = V_V(i6,1:N_V_V(i6))
        V_V(i6,k5+1) =i5
        V_V(i6,k6+1:N_V_V(i6)+1)=vv(k6:N_V_V(i6))
        N_V_V(i6)=N_V_V(i6)+1    
    
    
    
    endif
    endif

    endif
    enddo




deallocate(Permission_E)


    END subroutine bond_flip_neck


    subroutine Insert_APOINT(RV,Vector_F,area_F,N_V,N_E,N_F,V_E,E_F,V_F,num_F)
        implicit none
        real*8 rv(:,:),Vector_F(:,:),area_F(:)
        integer N_V,N_E,N_F
        integer V_E(:,:),E_F(:,:),V_F(:,:),num_F
        integer i0,i,j,k,i1,i2,i3
        integer ivn,ie1,ie2,ie3,if1,if2,if3
        integer,allocatable:: E_Ft(:,:)
        allocate(E_Ft(0:ubound(area_F,1),1:3))
        i1 = V_F(num_F,1);i2 = V_F(num_F,2);i3 = V_F(num_F,3)
        ivn = N_V+1
        RV(ivn,1:3) = (rv(i1,1:3)+rv(i2,1:3)+rv(i3,1:3))/3.0 + &
                  vector_F(num_F,1:3)/area_F(num_F) * 0.68
        ie1 =N_E+1; ie2=N_E+2; ie3=N_E+3
        if1 =N_F+1; if2=N_F+2; if3=N_F+3
        V_E(ie1,1) =ivn;V_E(ie1,2) =i1;
        V_E(ie2,1) =ivn;V_E(ie2,2) =i2;
        V_E(ie3,1) =ivn;V_E(ie3,2) =i3;
        
        E_F(if1,1) = E_F(num_F,1);E_F(if1,2) = -ie2;E_F(if1,3) = ie1;
        E_F(if2,1) = E_F(num_F,2);E_F(if2,2) = -ie3;E_F(if2,3) = ie2;
        E_F(if3,1) = E_F(num_F,3);E_F(if3,2) = -ie1;E_F(if3,3) = ie3;
        N_V =ivn; N_E= ie3; N_F=if3
        ! delete the face Num_F
        E_Ft(1:N_F,1:3)=E_F(1:N_F,1:3)
        E_F(Num_F:N_F-1,1:3) = E_Ft(Num_F+1:N_F,1:3)
        E_F(N_F,1:3) = 0
        N_F=N_F-1
        deallocate(E_Ft)
        
    end subroutine Insert_APoint
    
    
    
    subroutine Vertex_averaging(rv,N_V,N_F,F_V,V_F,V_V,N_F_V ,&
    & area_F,Vector_F,Norm_F,area_F_zero,area_tot,area_tot_zero)
        implicit none
        real*8 rv(:,:)
        real*8,allocatable::rv0(:,:)
        !real*8 fi_d(:),fi1,fi2,fi3
        !real*8 fi_e(:)
        !real*8 fi_around(1:20)
        real*8 sum1,sum2

        ! %%%%%%%%%%%%%%%%%%%%%
        integer N_V,N_F,N_E
        integer F_V(:,:), V_F(:,:), V_V(:,:)
        integer N_F_V(:)

        ! %%%%%%%%%%%%%%
        real*8 area_tot, area_tot_zero
        real*8  Vector_F(:,:), Area_F(:), Norm_F(:,:)
        real*8 area_F_Zero(:)
        real*8, allocatable::area_V(:),area_F1(:)
      
        !
        real*8 X_centroid(1:20,1:3), area_around(1:20)
        real*8 area_around_new(1:20)
        real*8 Norm_F1(1:20,1:3), Norm_V(1:3)
        real*8 lamda
        real*8 rv_new(1:3),rv_old(1:3)
        real*8 rv_new1(1:3), rv_new2(1:3)

        ! %%%%%%%%%%%%%%

        integer i,j,k,i1,i2,i3,i4,i5,i0,I_temp(1:20),I_temp2(1:20),I_temp3(1:20)
        integer N_t

        ! %%%%%%%%%%
        real*8 R(1:3,1:3), R1(1:3), R2(1:3), R3(1:3)
        real*8 a1,a2,a3,v1(1:3),v2(1:3),v3(1:3)
        real*8 len_V(1:20),len_VV

        ! %%%%%%%%%%
        real*8 dd(1:20),ddd
        real*8 area_Max,area_Min
        integer N
        ! %%%%%%%%%%
        N = ubound(area_F,1)
        allocate(Area_F1(1:N))
        allocate(area_V(1:N))
        area_F1(1:N_F)=area_F(1:N_F)
        ! %%%%%%%%%%%

        do i=1, N_V
            N_t=N_F_V(i)
            I_temp(1:20)=F_V(i,1:20);I_temp(N_t+1)=I_temp(1)
            I_temp2(1:20)=V_V(i,1:20);I_temp2(N_t+1)=I_temp2(1)
            ddd=0.0
            do j=1,N_t
                dd(j)= (area_F1(I_temp(j))/area_F_zero(I_temp(j)) - 1.2 ) * (0.8 - area_F1(I_temp(j))/area_F_zero(I_temp(j)) ) ! + is OK
                ddd = dd(j)+ddd
            enddo

            if(ddd.lt.120.0)then

                v2=0.0
                Norm_V=0.0

                rv_old=rv(i,1:3)

                area_MAX=-100
                area_min=100000
                do j=1,N_t
                    i1=i; i2=I_temp2(j); i3=I_temp2(j+1)
                    R1(1:3)=Rv(i1,1:3);R2(1:3)=Rv(i2,1:3);R3(1:3)=Rv(i3,1:3);
                call Get_area_local(R1,R2,R3,a1, v1)


                area_around(j)=a1
                Norm_F1(j,1:3)=V1
                X_centroid(j,1:3)=(R1+R2+R3)/3.0
                v2=v2+a1*X_centroid(j,1:3)
                Norm_V=Norm_V+v1
                if(a1.gt.area_max)area_max=a1
                if(a1.lt.area_min)area_min=a1
                enddo ! end j

                rv_new=v2/sum(area_around(1:N_t))
                lamda=(sum(rv_new(1:3)*Norm_V)-sum(rv(i,1:3)*Norm_V) )/sum(Norm_V*Norm_V)

                rv(i,1:3)=rv_new-lamda*Norm_V
  
                ! %%%% Get new fi

                ! %% fi0_n= (sum_i  fi_i(A_i+A_{i+1} ) - sum_i fi_i{}'.. +fi_0 sum A_i )/sum A_i'
                do j=1,N_t
                    i1=i; i2=I_temp2(j); i3=I_temp2(j+1)
                    R1(1:3)=Rv(i1,1:3);R2(1:3)=Rv(i2,1:3);R3(1:3)=Rv(i3,1:3);

                    call Get_area_local(R1,R2,R3,a1, v1)

                    area_around_New(j)=a1
                enddo ! end j



            endif ! end of ddd


        enddo ! end i

        ! %%%%%%%%%%%%%%%%%%%%%%% reconsider fi

        call Get_area_F( rv, V_F,  N_F,   Area_F, Vector_F, Norm_F)


        do i=1,N_V
            area_V(i)=0.0
            N_t=N_F_V(i)
            I_temp(1:20)=F_V(i,1:20)
            do j=1,N_t
                area_V(i)=area_V(i)+area_F(I_temp(j))/3.0
            enddo
        enddo

        area_tot=sum(area_V(1:N_V))

        deallocate(Area_F1)
        deallocate(area_V)

    END subroutine vertex_averaging
    
       
    subroutine Vertex_averaging_fixshape(rv,rvf,N_V,N_F,F_V,V_F,V_V,N_F_V, vvNear ,&
    & area_F,Vector_F,Norm_F,Norm_F0,area_F_zero,area_tot,area_tot_zero)
        implicit none
        real*8 rv(:,:),rvf(:,:)
        real*8,allocatable::rv0(:,:)
        !real*8 fi_d(:),fi1,fi2,fi3
        !real*8 fi_e(:)
        !real*8 fi_around(1:20)
        real*8 sum1,sum2

        ! %%%%%%%%%%%%%%%%%%%%%
        integer N_V,N_F,N_E
        integer F_V(:,:), V_F(:,:), V_V(:,:)
        integer N_F_V(:)
        integer vvNear(:)
        ! %%%%%%%%%%%%%%
        real*8 area_tot, area_tot_zero
        real*8  Vector_F(:,:), Area_F(:), Norm_F(:,:),Norm_F0(:,:)
        real*8 area_F_Zero(:)
        real*8, allocatable::area_V(:),area_F1(:)
      
        !
        real*8 X_centroid(1:20,1:3), area_around(1:20)
        real*8 area_around_new(1:20)
        real*8 Norm_F1(1:20,1:3), Norm_V(1:3)
        real*8 lamda
        real*8 rv_new(1:3),rv_old(1:3)
        real*8 rv_new1(1:3), rv_new2(1:3)

        ! %%%%%%%%%%%%%%

        integer i,j,k,i1,i2,i3,i4,i5,i0,I_temp(1:20),I_temp2(1:20),I_temp3(1:20)
        integer N_t

        ! %%%%%%%%%%
        real*8 R(1:3,1:3), R1(1:3), R2(1:3), R3(1:3),R0(1:3),normt(1:3)
        real*8 al,bet,dn,dmin,ddt
        integer is_ON,is_IN
        real*8 a1,a2,a3,v1(1:3),v2(1:3),v3(1:3)
        real*8 len_V(1:20),len_VV

        ! %%%%%%%%%%
        real*8 dd(1:20),ddd
        real*8 area_Max,area_Min
        integer N
        ! %%%%%%%%%%
        N = ubound(area_F,1)
        allocate(Area_F1(1:N))
        allocate(area_V(1:N))
        area_F1(1:N_F)=area_F(1:N_F)
        ! %%%%%%%%%%%

        do i=1, N_V
            N_t=N_F_V(i)
            I_temp(1:20)=F_V(i,1:20);I_temp(N_t+1)=I_temp(1)
            I_temp2(1:20)=V_V(i,1:20);I_temp2(N_t+1)=I_temp2(1)
            ddd=0.0
            do j=1,N_t
                dd(j)= (area_F1(I_temp(j))/area_F_zero(I_temp(j)) - 1.2 ) * (0.8 - area_F1(I_temp(j))/area_F_zero(I_temp(j)) ) ! + is OK
                ddd = dd(j)+ddd
            enddo

            if(ddd.lt.120.0)then

                v2=0.0
                Norm_V=0.0

                rv_old=rv(i,1:3)

                area_MAX=-100
                area_min=100000
                do j=1,N_t
                    i1=i; i2=I_temp2(j); i3=I_temp2(j+1)
                    R1(1:3)=Rv(i1,1:3);R2(1:3)=Rv(i2,1:3);R3(1:3)=Rv(i3,1:3);
                call Get_area_local(R1,R2,R3,a1, v1)


                area_around(j)=a1
                Norm_F1(j,1:3)=V1
                X_centroid(j,1:3)=(R1+R2+R3)/3.0
                v2=v2+a1*X_centroid(j,1:3)
                Norm_V=Norm_V+v1
                if(a1.gt.area_max)area_max=a1
                if(a1.lt.area_min)area_min=a1
                enddo ! end j

                rv_new=v2/sum(area_around(1:N_t))
                lamda=(sum(rv_new(1:3)*Norm_V)-sum(rv(i,1:3)*Norm_V) )/sum(Norm_V*Norm_V)

             
                    rv(i,1:3)=rv_new -lamda*Norm_V
            
                k =vvnear(i)
                R0 = rv(i,1:3)
                j=1; is_ON =0
                do while(j.le.N_F_V(k).and. IS_on.eq.0)
                   i1=F_V(k,j) 
                   R1 = rvf(V_F(i1,1),1:3) ;R2 = rvf(V_F(i1,2),1:3) ;R3 = rvf(V_F(i1,3),1:3) 
                   normt(1:3)=norm_F0(i1,1:3)
                   call Pos_Point_Triangle(R0,R1,R2,R3,Normt,al,bet,dn,is_ON,is_IN)
  
                   
                   
                   
                    if(IS_ON.eq.1)then
                       Rv(i,1:3) = Rv(i,1:3) - dn*Normt
                      
                       dmin = 100000;i3=0;i2=1
                       do while(i2.le.3.and.i3.eq.0)
                           ddt = sum((rv(i,:) -rvf(V_F(i1,i2),1:3) )**2 )
                           if(ddt.lt.dmin)then
                               dmin=ddt; vvnear(i) = V_F(i1,i2);i3=1
                           endif
                           i2=i2+1
                       enddo
                       
                    endif
                    j=j+1
                enddo
                
                
   



            endif ! end of ddd
      

        enddo ! end i

        ! %%%%%%%%%%%%%%%%%%%%%%% reconsider fi

        call Get_area_F( rv, V_F,  N_F,   Area_F, Vector_F, Norm_F)


        do i=1,N_V
            area_V(i)=0.0
            N_t=N_F_V(i)
            I_temp(1:20)=F_V(i,1:20)
            do j=1,N_t
                area_V(i)=area_V(i)+area_F(I_temp(j))/3.0
            enddo
        enddo

        area_tot=sum(area_V(1:N_V))

    deallocate(Area_F1)
        deallocate(area_V)

    END subroutine vertex_averaging_fixshape
    
    
       
    subroutine Vertex_averaging_fixshape_Refine(rv,rvf,N_V,N_F,F_V,V_F,V_V,N_F_V, &
    F_Vf, V_Ff,V_Vf,N_F_Vf,vvNear ,&
    & area_F,Vector_F,Norm_F,Norm_F0,area_F_zero,area_tot,area_tot_zero)
        implicit none
        real*8 rv(:,:),rvf(:,:)
        real*8,allocatable::rv0(:,:)
        !real*8 fi_d(:),fi1,fi2,fi3
        !real*8 fi_e(:)
        !real*8 fi_around(1:20)
        real*8 sum1,sum2

        ! %%%%%%%%%%%%%%%%%%%%%
        integer N_V,N_F,N_E
        integer F_V(:,:), V_F(:,:), V_V(:,:)
        integer N_F_V(:)
        integer F_Vf(:,:), V_Ff(:,:), V_Vf(:,:)
        integer N_F_Vf(:)
        integer vvNear(:)
        ! %%%%%%%%%%%%%%
        real*8 area_tot, area_tot_zero
        real*8  Vector_F(:,:), Area_F(:), Norm_F(:,:),Norm_F0(:,:)
        real*8 area_F_Zero(:)
        real*8, allocatable::area_V(:),area_F1(:)
      
        !
        real*8 X_centroid(1:20,1:3), area_around(1:20)
        real*8 area_around_new(1:20)
        real*8 Norm_F1(1:20,1:3), Norm_V(1:3)
        real*8 lamda
        real*8 rv_new(1:3),rv_old(1:3)
        real*8 rv_new1(1:3), rv_new2(1:3)

        ! %%%%%%%%%%%%%%

        integer i,j,k,i1,i2,i3,i4,i5,i0,I_temp(1:20),I_temp2(1:20),I_temp3(1:20)
        integer N_t

        ! %%%%%%%%%%
        real*8 R(1:3,1:3), R1(1:3), R2(1:3), R3(1:3),R0(1:3),normt(1:3)
        real*8 al,bet,dn,dmin,ddt
        integer is_ON,is_IN
        real*8 a1,a2,a3,v1(1:3),v2(1:3),v3(1:3)
        real*8 len_V(1:20),len_VV

        ! %%%%%%%%%%
        real*8 dd(1:20),ddd
        real*8 area_Max,area_Min
        integer N
        ! %%%%%%%%%%
        N = ubound(area_F,1)
        allocate(Area_F1(1:N))
        allocate(area_V(1:N))
        area_F1(1:N_F)=area_F(1:N_F)
        ! %%%%%%%%%%%

        do i=1, N_V
            N_t=N_F_V(i)
            I_temp(1:20)=F_V(i,1:20);I_temp(N_t+1)=I_temp(1)
            I_temp2(1:20)=V_V(i,1:20);I_temp2(N_t+1)=I_temp2(1)
            ddd=0.0
            do j=1,N_t
                dd(j)= (area_F1(I_temp(j))/area_F_zero(I_temp(j)) - 1.2 ) * (0.8 - area_F1(I_temp(j))/area_F_zero(I_temp(j)) ) ! + is OK
                ddd = dd(j)+ddd
            enddo

            if(ddd.lt.120.0)then

                v2=0.0
                Norm_V=0.0

                rv_old=rv(i,1:3)

                area_MAX=-100
                area_min=100000
                do j=1,N_t
                    i1=i; i2=I_temp2(j); i3=I_temp2(j+1)
                    R1(1:3)=Rv(i1,1:3);R2(1:3)=Rv(i2,1:3);R3(1:3)=Rv(i3,1:3);
                call Get_area_local(R1,R2,R3,a1, v1)


                area_around(j)=a1
                Norm_F1(j,1:3)=V1
                X_centroid(j,1:3)=(R1+R2+R3)/3.0
                v2=v2+a1*X_centroid(j,1:3)
                Norm_V=Norm_V+v1
                if(a1.gt.area_max)area_max=a1
                if(a1.lt.area_min)area_min=a1
                enddo ! end j

                rv_new=v2/sum(area_around(1:N_t))
                lamda=(sum(rv_new(1:3)*Norm_V)-sum(rv(i,1:3)*Norm_V) )/sum(Norm_V*Norm_V)

             
                    rv(i,1:3)=rv_new -lamda*Norm_V
            
                k =vvnear(i)
                R0 = rv(i,1:3)
                j=1; is_ON =0
                do while(j.le.N_F_Vf(k).and. IS_on.eq.0)
                   i1=F_Vf(k,j) 
                   R1 = rvf(V_Ff(i1,1),1:3) ;R2 = rvf(V_Ff(i1,2),1:3) ;R3 = rvf(V_Ff(i1,3),1:3) 
                   normt(1:3)=norm_F0(i1,1:3)
                   call Pos_Point_Triangle(R0,R1,R2,R3,Normt,al,bet,dn,is_ON,is_IN)
  
                   
                    if(IS_ON.eq.1)then
                       Rv(i,1:3) = Rv(i,1:3) - dn*Normt
                      
                       dmin = 100000;i3=0;i2=1
                       do while(i2.le.3.and.i3.eq.0)
                           ddt = sum((rv(i,:) -rvf(V_Ff(i1,i2),1:3) )**2 )
                           if(ddt.lt.dmin)then
                               dmin=ddt; vvnear(i) = V_Ff(i1,i2);i3=1
                           endif
                           i2=i2+1
                       enddo
                       
                    endif
                    j=j+1
                enddo
                
                
   



            endif ! end of ddd
      

        enddo ! end i

        ! %%%%%%%%%%%%%%%%%%%%%%% reconsider fi

        call Get_area_F( rv, V_F,  N_F,   Area_F, Vector_F, Norm_F)


        do i=1,N_V
            area_V(i)=0.0
            N_t=N_F_V(i)
            I_temp(1:20)=F_V(i,1:20)
            do j=1,N_t
                area_V(i)=area_V(i)+area_F(I_temp(j))/3.0
            enddo
        enddo

        area_tot=sum(area_V(1:N_V))

deallocate(Area_F1)
        deallocate(area_V)

    END subroutine vertex_averaging_fixshape_Refine
    
    
    subroutine Refine_surf_butterfly(Rv,N_V,N_E,N_F,V_E,E_F)
    implicit none
    integer N
     ! == readfe
    integer N_V,N_E,N_F
    real*8 rV(:,:),x,y,z,a

    integer V_E(:, :)
    integer E_F(:, :)


    real*8,allocatable:: rV0(:,:)
    integer N_V0,N_E0,N_F0
    integer,allocatable:: V_E0(:, :)
    integer,allocatable:: E_F0(:, :)
    integer,allocatable:: F_E0(:,:),V_F0(:, :)

    integer i,i1,i2,i3,i_E,F0,F1,F2
    integer V1,V2,V3,EN(1:6),E1,E2,E3,E4,I_F,I_V
    integer W1,W2,W3,W4,W5,W6,W7,W8,VE1,VE2
    real*8 a1,a2,a3,omega

    integer,allocatable:: E_E_Refine(:,:) ! 1:2
    real*8 r00,th

    N = ubound(Rv,1)
    allocate(rv0(1:N,1:3), V_E0(1:N,1:2),E_F0(1:N, 1:3), &
       F_E0(1:N,1:2),V_F0(1:N, 1:3),E_E_Refine(1:N,1:2))

    ! %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    ! step 1
    do i=1, N_F
      i1=E_F(i,1);i2=E_F(i,2);i3=E_F(i,3)
      F_E0(abs(i1), int( (3-i1/abs(i1) ) /2 ) )=i
      F_E0(abs(i2), int( (3-i2/abs(i2) ) /2 ) )=i
      F_E0(abs(i3), int( (3-i3/abs(i3) ) /2 ) )=i
    enddo
      do i=1,N_F
      V_F0(i,1)=V_E(abs(E_F(i,1)), int( (3-E_F(i,1)/abs(E_F(i,1)) ) /2 ) )
      V_F0(i,2)=V_E(abs(E_F(i,2)), int( (3-E_F(i,2)/abs(E_F(i,2)) ) /2 ) )
      V_F0(i,3)=V_E(abs(E_F(i,3)), int( (3-E_F(i,3)/abs(E_F(i,3)) ) /2 ) )
      enddo

    ! %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    ! step 2
    Rv0(1:N_V,1:3)=Rv(1:N_V,1:3)
    N_V0=N_V; N_E0 = N_E; N_F0 = N_F
    V_E0(1:N_E,1:2)= V_E(1:N_E,1:2); E_F0(1:N_F,1:3)= E_F(1:N_F,1:3)
    omega = 1.00/16.0;
    a1 = 0.5;  a2 =  2*omega ; a3 =-omega
    !a1 = 0.5; a2 =0; a3=0
    ! %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    ! step 3 determine the mid points of every edge
    
    do i_E = 1, N_E0
      ! %%%
      ! find W1 to W8
	    ! w1 and w2
	    W1 = V_E0 (i_E,1);W2 = V_E0 (i_E,2)
	    ! w3 and w4
	    F1 = F_E0(i_E,1); F2 = F_E0(i_E,2);
	    do i1 =1, 3
	    	if(V_F0(F1,i1).ne.W1.and.V_F0(F1,i1).ne.W2)W3 = V_F0(F1,i1)
	    	if(V_F0(F2,i1).ne.W1.and.V_F0(F2,i1).ne.W2)W4 = V_F0(F2,i1)
	    enddo
	    do i1 =1,3
	    	if( (V_E0( abs(E_F0(F1,i1)), 1) .eq.W1.and. V_E0( abs(E_F0(F1,i1)), 2) .eq.W3).or. &
	    	   &(V_E0( abs(E_F0(F1,i1)), 2) .eq.W1.and. V_E0( abs(E_F0(F1,i1)), 1) .eq.W3) )  E1 =abs(E_F0(F1,i1))
	    	if( (V_E0( abs(E_F0(F1,i1)), 1) .eq.W2.and. V_E0( abs(E_F0(F1,i1)), 2) .eq.W3).or. &
	    	   &(V_E0( abs(E_F0(F1,i1)), 2) .eq.W2.and. V_E0( abs(E_F0(F1,i1)), 1) .eq.W3) )  E2 =abs(E_F0(F1,i1))
	    	if( (V_E0( abs(E_F0(F2,i1)), 1) .eq.W1.and. V_E0( abs(E_F0(F2,i1)), 2) .eq.W4).or. &
	    	   &(V_E0( abs(E_F0(F2,i1)), 2) .eq.W1.and. V_E0( abs(E_F0(F2,i1)), 1) .eq.W4) )  E3 =abs(E_F0(F2,i1))
	    	if( (V_E0( abs(E_F0(F2,i1)), 1) .eq.W2.and. V_E0( abs(E_F0(F2,i1)), 2) .eq.W4).or. &
	    	   &(V_E0( abs(E_F0(F2,i1)), 2) .eq.W2.and. V_E0( abs(E_F0(F2,i1)), 1) .eq.W4) )  E4 =abs(E_F0(F2,i1))   
	    enddo

	    do i1 =1, 2
	    do i2 =1,3
	    			i_V = V_F0( F_E0 ( E1, i1), i2);
	    	 if(i_V.ne.W1.and.i_V.ne.W3.and.i_V.ne.W2) W5 =i_V
	    			i_V = V_F0( F_E0 ( E2, i1), i2);
	    	 if(i_V.ne.W1.and.i_V.ne.W3.and.i_V.ne.W2) W6 =i_V
	    	 	     i_V = V_F0( F_E0 ( E3, i1), i2);
	    	 if(i_V.ne.W1.and.i_V.ne.W4.and.i_V.ne.W2) W7 =i_V
	    	 	     i_V = V_F0( F_E0 ( E4, i1), i2);
	    	 if(i_V.ne.W1.and.i_V.ne.W4.and.i_V.ne.W2) W8 =i_V
	    enddo
	    enddo

      ! %%%%%%%%%%%%%%%%%%%% adding a new mid point into surface
      !  add two new edges
	    N_V = N_V0 + i_E
	    V_E(i_E,1) = W1
	    V_E(i_E,2) = N_V
	    E_E_Refine(i_E,1)=i_E
	    N_E = N_E + 1
	    V_E(N_E,1) = N_V
	    V_E(N_E,2) = W2
	    E_E_Refine(i_E,2)=N_E
		
	    RV(N_V,1:3) = a1* (RV0(W1,1:3) + RV0(W2,1:3) ) + a2 * (RV0(W3,1:3) + RV0(W4,1:3) ) + &
				     a3 * (RV0(W5,1:3) + RV0(W6,1:3) + RV0(W7,1:3) + RV0(W8,1:3) )
    enddo


    ! %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    ! step 4 build the relationships of these new inserting points with old points
    !        Only for label variables

    do i_F = 1, N_F0
       F0 = i_F
	    ! creat 3 new edges in the center of the triangle
	    W1 = N_V0 +  abs(E_F0(F0,2))
	    W2 = N_V0 +  abs(E_F0(F0,3))
	    W3 = N_V0 +  abs(E_F0(F0,1))
	    N_E = N_E +1
	    E1 =N_E
	    V_E( N_E,1) =W1;V_E( N_E,2) =W2;
	    N_E = N_E +1
	    E2 =N_E
	    V_E( N_E,1) =W2;V_E( N_E,2) =W3;
	    N_E = N_E +1
	    E3 =N_E
	    V_E( N_E,1) =W3;V_E( N_E,2) =W1;
        ! find six edges
	    do i1 =1,3
		    i2 = mod(i1+1,3); if(i2.eq.0)i2 =3
		    V1 = V_F0(F0,i1);V2 = V_F0(F0,i2);
            
		    VE1 = V_E(E_E_Refine(abs(E_F(F0,i1)),1),1);
            VE2 = V_E(E_E_Refine(abs(E_F(F0,i1)),1),2)
		    if(VE1.eq.V1.or.VE2.eq.V1)then
			    EN( i1*2 -1) = E_E_Refine(abs(E_F(F0,i1)),1)
			    EN( i1*2   ) = E_E_Refine(abs(E_F(F0,i1)),2)
		    endif
		    if(VE1.eq.V2.or.VE2.eq.V2)then
			    EN( i1*2 -1) = -E_E_Refine(abs(E_F(F0,i1)),2)
			    EN( i1*2   ) = -E_E_Refine(abs(E_F(F0,i1)),1)
		    endif
	    enddo


	    ! label the facets
	    E_F(F0, 1) = E1;E_F(F0, 2) = E2;E_F(F0, 3) = E3;
	    N_F = N_F +1
	    E_F(N_F,1) = EN(1);E_F(N_F,2) = -E2;E_F(N_F,3) = EN(6)
	    N_F = N_F +1
	    E_F(N_F,1) = EN(2);E_F(N_F,2) = EN(3);E_F(N_F,3) = -E3
	    N_F = N_F +1
	    E_F(N_F,1) = -E1 ;E_F(N_F,2) = EN(4);E_F(N_F,3) = EN(5)
        enddo
 deallocate(rv0, V_E0,E_F0, &
       F_E0,V_F0,E_E_Refine)

    END subroutine Refine_surf_butterfly
    
    
    
      
    subroutine WKONTRI_POS(al,bet,k)
        real*8 al,bet,dr
        integer k
        dr = 0.00000001
        k=-1
        if(abs(al).lt.dr.and.bet.le.1+dr.and.bet.ge.0-dr)then
            k = 3    ! on border 1
        endif
        if(abs(bet).lt.dr.and.al.le.1+dr.and.al.ge.0-dr)then
            k = 1    ! on border 2
        endif
        if(abs(al+bet-1.0).lt.dr)then
            k = 2    ! on border 3
        endif
        
        if(al+bet.lt.1.0-dr.and.al.gt.dr.and.bet.gt.dr)then
            k = 0        !  inside the triangle
        endif 
    end subroutine WKONTRI_POS
   
     subroutine WKONTRI_POS2(al,bet,k)
        real*8 al,bet,dr
        integer k
        
         k=-1
        
        if(al+bet.le.1.0.and.al.ge.0.and.bet.ge.0)then
            k = 0        !  inside the triangle
        endif 
        if(al+bet.eq.1.0.and.al.ge.0.and.bet.ge.0)k=2
        if(al+bet.le.1.0.and.al.eq.0.and.bet.ge.0)k=3
        if(al+bet.le.1.0.and.al.ge.0.and.bet.eq.0)k=1
        if(al.eq.0.and.bet.eq.0)k=4
        if(al.eq.1.and.bet.eq.0)k=5
        if(al.eq.0.and.bet.eq.1)k=6
        
    end subroutine WKONTRI_POS2
    
    
    subroutine WKONTRI_Convert(C,kf,v,al,bet)
    type(CEL)::C
    integer kf
      real*8 v(1:3),al,bet
      real*8 y(1:2)
      y(1) = sum(c%e(kf,1,1:3)*v)
      y(2) = sum(c%e(kf,2,1:3)*v)
      al = sum(c%gup(kf,1,1:2)*y(1:2))
      bet = sum(c%gup(kf,2,1:2)*y(1:2))
    end subroutine WKONTRI_convert
    
    subroutine WKONTRI_border(al1,bet1,k0,al2,bet2,x,k)
        real*8 al1,bet1,al2,bet2,x
        integer k,k0
        real*8 xx,y,d
        real*8 dr
        ! border 1
        dr = 0.0000000001
        k=-1
        x=0.0
        xx=0.0
        if(al2.ne.0.and.al1.gt.dr)then
            xx = -al1/al2
            y=xx*bet2+bet1
            if(xx.gt. 0.0 .and.y.ge.0.and.y.le.1.and.k0.ne.3)then
                k=3
                x=xx
            endif
        endif
        ! border 2
        if(bet2.ne.0.and.bet1.gt.dr)then
            xx = -bet1/bet2
            y=xx*al2+al1
            if(xx.gt. 0.0 .and.y.ge.0.and.y.le.1.and.k0.ne.1)then
                k=1
                x=xx
            endif
        endif            
        
        if(al2+bet2.ne.0.and.al1+bet1.lt.1-dr)then
            xx= (1-al1-bet1)/(al2+bet2)
            if(xx.gt. 0.0 .and. xx*al2+al1.gt.0.and.xx*al2+al1.lt.1 &
            & .and. xx*bet2+bet1.gt.0.and.xx*bet2+bet1.lt.1.and.k0.ne.2 )then
                k=2
                x=xx
            endif
            
        endif
        
            
    end subroutine WKONTRI_border
    
    
    subroutine WKONTRI_crossborder(v1,v2,n1,n2,u)
    real*8 v1(1:3),v2(1:3),n1(1:3),n2(1:3),u(1:3)
    real*8 a(1:3),b(1:3),ab
    real*8 cc(1:3)
    integer sign1
      a = sum(v1*u)*u/sum(u**2)
      b= v1 - a
      cc(1) = n1(2)*u(3)- u(2)*n1(3)
      cc(2) = n1(3)*u(1)- u(3)*n1(1)
      cc(3) = n1(1)*u(2)- u(1)*n1(2)
      ab=sum(b*cc)
      if(ab.ge.0)sign1=1
      if(ab.lt.0)sign1=-1
      cc(1) = n2(2)*u(3)- u(2)*n2(3)
      cc(2) = n2(3)*u(1)- u(3)*n2(1)
      cc(3) = n2(1)*u(2)- u(1)*n2(2)
      cc(1:3) = cc(1:3) / (sum(cc**2))**0.5
      
      b = (sum(b**2))**0.5 *cc*sign1
      
      v2 = a+b
    end subroutine WKONTRI_crossborder
    
    subroutine WKONTRI_alloc(C)
    type(CEL)::C
    allocate(C%RF(1:C%N_F,1:3,1:3),C%e(1:C%N_F,1:2,1:3), &
      C%gdn(1:C%N_F,1:2,1:2),C%gup(1:C%N_F,1:2,1:2))
    allocate(C%NF(1:C%N_F,1:3,1:3),C%thF(1:C%N_F,1:3),C%LenF(1:C%N_F,1:3))
    end subroutine WKONTRI_alloc
    
    subroutine WKONTRI_dealloc(C)
    type(CEL)::C
    deallocate(C%RF,C%e, &
      C%gdn,C%gup)
    deallocate(C%NF,C%thF,C%LenF)
     end subroutine WKONTRI_dealloc
    
     
    
    subroutine WKONTRI_Init(C)
        type(Cel)::C
        integer i,j,i1,i2,i3
        real*8 det,csth
        real*8,allocatable:: norm_V(:,:)
        allocate(norm_V(1:C%N_V,1:3))
        
        do i=1,C%N_V
            norm_V(i,1:3) = 0
            do j=1,C%N_V_V(i)
                i1 = C%F_V(i,j)
                norm_V(i,1:3) =norm_V(i,1:3)+C%Norm_F(i1,1:3)*C%area_F(i1)     
            enddo
            norm_V(i,1:3)=norm_V(i,1:3)/sum(norm_V(i,1:3)**2)**0.5
        enddo
        
         
        do i=1, C%N_F
            i1=c%V_F(i,1);i2=c%V_F(i,2);i3=c%V_F(i,3);
            c%RF(i,1,1:3)=c%rv(i1,1:3)
            c%RF(i,2,1:3)=c%rv(i2,1:3)
            c%RF(i,3,1:3)=c%rv(i3,1:3)
            c%e(i,1,1:3) = c%RF(i,2,1:3)-c%RF(i,1,1:3)
            c%e(i,2,1:3) = c%RF(i,3,1:3)-c%RF(i,1,1:3)
            c%gdn(i,1,1) = sum(c%e(i,1,1:3)*c%e(i,1,1:3))
            c%gdn(i,1,2) = sum(c%e(i,1,1:3)*c%e(i,2,1:3))
            c%gdn(i,2,1) = c%gdn(i,1,2)
            c%gdn(i,2,2) = sum(c%e(i,2,1:3)*c%e(i,2,1:3))
            det = c%gdn(i,2,2)*c%gdn(i,1,1) - c%gdn(i,2,1)**2
            c%gup(i,:,:) =0.0
            if(det.ne.0)then
                c%gup(i,1,1) =  c%gdn(i,2,2)/det
                c%gup(i,1,2) = - c%gdn(i,2,1)/det
                c%gup(i,2,1) =  -c%gdn(i,1,2)/det
                c%gup(i,2,2) =  c%gdn(i,1,1)/det
            endif
            
            C%nF(i,1,1:3) =  norm_V(i1,1:3)
            C%nF(i,2,1:3) =  norm_V(i2,1:3)
            C%nF(i,3,1:3)  = norm_V(i3,1:3)
            
            C%LenF(i,1) = sum( (c%rv(i1,1:3)-c%rv(i2,1:3))**2)**0.5
            C%LenF(i,2) = sum( (c%rv(i1,1:3)-c%rv(i3,1:3))**2)**0.5
            C%LenF(i,3) = sum( (c%rv(i2,1:3)-c%rv(i3,1:3))**2)**0.5
            C%thF(i,1:3)=0
            csth = sum(norm_V(i1,1:3)*norm_V(i2,1:3))
            if(abs(csth).lt.1)C%thF(i,1) = acos(csth)
            csth = sum(norm_V(i1,1:3)*norm_V(i3,1:3))
            if(abs(csth).lt.1)C%thF(i,2) = acos(csth)
            csth = sum(norm_V(i3,1:3)*norm_V(i2,1:3))
            if(abs(csth).lt.1)C%thF(i,3) = acos(csth)
            
        enddo
    
       deallocate(Norm_V)
    END subroutine WKONTRI_INIT
    
    
    subroutine WKONTRI_Interp_Surf(n1,n2,n3,len1,len2,len3,th1,th2,th3,al,bet,DR)
       implicit none
        real*8 n1(1:3),n2(1:3),n3(1:3)
        real*8 n(1:3),len1,len2,len3,th1,th2,th3
        real*8 al,bet,al1,bet1
        real*8 DR(1:3)
        real*8 h1,h2,h3
        n = al*n2 +bet*n3+(1-al-bet)*n1
        n = n/sum(n**2)**0.5
        h1 = len1*(th1*(1-th1**2/24)*al*(1-al)/2 - (th1*(1-th1**2/24))**3*al*al*(1-al)**2/8.0)
        h2 = len2*(th2*(1-th2**2/24)*bet*(1-bet)/2 - (th2*(1-th2**2/24))**3*bet*bet*(1-bet)**2/8.0)
        h3 = 0
        if(al+bet.gt.0.00001)then
            al1= al/(al+bet);bet1=bet/(al+bet);
            h3 = len3*(th3*(1-th3**2/24)*bet1*al1/2 - (th3*(1-th3**2/24))**3*bet1*bet1*al1**2/8.0)
        endif
        DR = ( (1-bet)*h1+(1-al)*h2+(al+bet)*h3 ) * n
        
    END subroutine WKONTRI_Interp_Surf
    
    
    
    subroutine WKONTRI_crossvertex(C,v1,kv,v2,kf,r0)
      implicit none
      type(Cel)::C
      real*8 v1(1:3),v2(1:3),r0(1:3)
      real*8 n(1:3),al,bet
      real*8 v1_small(1:3),v3(1:3),v4(1:3)
      real*8 albet_min,al_min,bet_min
      integer kv,kf,k,  imin
      integer i,j
      r0 = c%rv(kv,1:3)
      v1_small= v1/sum(v1**2)**0.5*0.001
      albet_min = -1000.0
      do i=1,c%N_V_V(kv)
          j=c%F_V(kv,i)
          n= c%norm_F(j,1:3)
          v3=v1_small - sum(v1_small*n)*n
          v4 = v3 + r0 - C%rf(j,1,1:3)
          call WKONTRI_Convert(c,j,v4,al,bet)
          call WKONTRI_POS2(al,bet,k)
          if(k.ne.-1)then
            v2 = v3 /sum(v3**2)**0.5
            kf = j
          endif
          if(al.lt.0.and.bet.gt.0)then
               if(al.gt.albet_min)then
                   albet_min=al; imin=i
                   al_min=0.001; bet_min=bet
               endif
          endif
          if(al.gt.0.and.bet.lt.0)then
               if(bet.gt.albet_min)then
                   albet_min=bet; imin=i
                   al_min=al; bet_min=0.001
               endif
          endif
          
      enddo
      
      if(k.eq.-1)then
         kf = c%F_V(kv,imin)
         v3(1:3) = al_min *C%e(kf,1,1:3)+bet_min*c%e(kf,2,1:3)
         v2 = v3 /sum(v3**2)**0.5
      endif
      
      
    end subroutine WKONTRI_crossvertex
  
    
         
    subroutine WKONTRI_MOVE2(C,v1,k_F_start,al0,bet0,r0,k_F_end,al_end,bet_end,rv_end) !,T_ave,i_check)
        implicit none
        type(cel)::C
        integer i,j,k,k0,k3, T_ave,i_check
        integer k_border
        integer k_crossborder,k_approachborder,k_movevertex, k_vertexNo
        real*8 r0(1:3),r1(1:3),r2(1:3),u(1:3)
        real*8 n1(1:3),n2(1:3)
        real*8 al1,bet1,al2,bet2,al0,bet0,al1u,bet1u
        real*8 al3,bet3
        real*8 v0(1:3),v1(1:3), v2(1:3),v3(1:3),d
        real*8 v1u(1:3)
        real*8 dd,x
        integer sw,kk,k_F_end,kv,k_max
        real*8 al_end,bet_end,rv_end(1:3)
        integer k_F_start
        integer k_vertex
            al_end=al0
            bet_end=bet0
            K_F_end=K_F_start
            rv_end=r0
            v1 = v1 - sum(v1*c%norm_F(k_F_start,1:3))*c%norm_F(k_F_start,1:3)
            d  = sum(v1**2)**0.5
            dd = d
            ! %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            ! k_vertex=0
            sw=1
            i=k_F_start
            kk=0
            k_max = 100

        do while(sw.eq.1)    
                    !if(i_check.eq.1.and.T_ave.eq.1)then
                         
                    !endif
            kk=kk+1
            v1u = v1 / (sum(v1**2))**0.5
            v1= dd*v1u
            call WKONTRI_Convert(C,i,v1,al1,bet1)
            al3 = al0 + al1
            bet3 = bet0 +bet1
            al1u = al1/dd
            bet1u= bet1/dd
                    !   
            call WKONTRI_POS2(al3,bet3,k)
            
            ! %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            ! k ne -1 means next move is inside the triangle
            ! so just stop the move and return. END POINT
                k_crossborder = 0  
                k_approachborder=0  ! if k_app=1 means r0+v1 will reach over the border
                k_movevertex  = 0   ! if k_move=1 means r0 is on the vertex and v1 is moving off the triangle
                                    ! the only solution is to move r0 inside the triangle.
            if(k.ne.-1)then    
                sw=-1    ! end point!
                k_f_end=i
                rv_end = r0 +v1
                al_end=al3
                bet_end=bet3
            else
            ! k eq -1 means next move is outside the triangle
            ! We first need to determine whether al0 bet0 is inside the triangle.

                call WKONTRI_POS2(al0,bet0,k0)
                
                if(k0.eq.0)k_approachborder = 1
                
                if(k0.gt.0)then     ! r0 is on the border or is vertex
                    al3 = al0 + 0.0001*al1u;bet3 = bet0 + 0.0001*bet1u;
                    call WKONTRI_POS2(al3,bet3,k3)
                    if(k3.eq.0)K_approachborder = 1
                    if(k3.eq.-1)then
                         if(k0.le.3)then
                             K_crossborder=1
                             k_border = k0 ! r0 is on the border but v1 is running outward
                         endif
                         if(k0.gt.3)then
                             k_movevertex = 1
                             k_vertexNo = k0 -3
                         endif
                    endif
                    if(k3.gt.0)then ! We are not going to handle this situation so we quit.
                        sw = 0; k_f_end=i
                        rv_end = r0;al_end= al0
                        bet_end= bet0          
                    endif
                endif  ! END of k0
                
                ! r0+v1*x to approach border
                if(k_approachborder.eq.1)then
                   call WKONTRI_border(al0,bet0,k0,al1u,bet1u,x,k)
                    al2 = al0+x*al1u; bet2 =bet0+x*bet1u
                    r1(1:3) = c%rf(i,1,1:3)+al2*c%e(i,1,1:3) +bet2*c%e(i,2,1:3)
                   call WKONTRI_POS2(al2,bet2,k3)
                   if(k3.gt.3)then  ! if al2 bet2 is vertex
                       k_movevertex = 1
                       k_vertexNo = k3 -1
                       al0 = al2; bet0  = bet2
                       dd = dd - (sum((r1-r0)**2))**0.5
                       r0= r1                       
                   else            ! if al2 bet2 is on the border
                       if(k3.ge.1.and.k3.le.3)then
                           k_crossborder = 1; k_border = k
                        else
                          kk=k_max+1   
                        endif
                       
                   endif
                endif
                
                if(k_crossborder.eq.1)then
                    dd = dd - (sum((r1-r0)**2))**0.5
                    j= c%F_F(i,k_border)
                    U = c%rv( c%v_E( abs(c%E_F(i,k_border)),2),1:3) -c%rv( c%V_E( abs(c%E_F(i,k_border)),1),1:3)
                    n1 = c%norm_F(i,1:3); n2=c%norm_F(j,1:3)
                    call WKONTRI_crossborder(v1u,v2,n1,n2,u)
                endif
                
                if(k_movevertex.eq.1)then
                    x=0.95
                    if(bet1u.eq.al1u*0.95.or.bet1u*0.95.eq.al1u)x=0.9
                    if(k_vertexNo.eq.1)then
                        al0 = 0.001; bet0 =al0*x
                    endif
                    if(k_vertexNo.eq.2)then
                        al0 = 1-0.001; bet0 =0.001*x
                    endif
                    if(k_vertexNo.eq.3)then
                        al0 = 0.001*x; bet0 =1-0.001
                    endif
                        r0  = C%rF(i,1,1:3) +al0*C%e(i,1,1:3) +bet0*C%e(i,2,1:3)
                        r1=r0;j=i;v2=v1    
                endif
            endif  ! END of outest if::k
            
            if(sw.ne.-1.and.k_movevertex.eq.0.and.k_approachborder &
                & .eq.0.and.k_crossborder.eq.0.or.kk.gt.k_max)then 
                        sw = 0; k_f_end=i
                        rv_end = r0;al_end= al0
                        bet_end= bet0          
            endif
            
            if(sw.eq.1)then
 
                i = j
                r0 = r1
                v0 = r0 - c%rf(i,1,:)
                v1=v2
                call WKONTRI_Convert(C,i,v0,al0,bet0)
                
            endif
            
        enddo
        
    
    END subroutine WKONTRI_MOVE2
    
subroutine WKONTRI_MOVE(C,v1,k_F_start,al0,bet0,r0,k_F_end,al_end,bet_end,rv_end,IS_BC_F)
        implicit none
        type(cel)::C
        integer i,j,k,k0
        real*8 r0(1:3),r1(1:3),r2(1:3),u(1:3)
        real*8 n1(1:3),n2(1:3)
        real*8 al1,bet1,al2,bet2,al0,bet0
        real*8 al3,bet3,al00,bet00,rv00(1:3)
        integer kf00
        real*8 v0(1:3),v1(1:3), v2(1:3),d
        real*8 dd,x,dr
        integer sw,kk,k_F_end,kv
        real*8 al_end,bet_end,rv_end(1:3)
        integer k_F_start
        integer k_vertex
        integer i1,i2,i3
        integer IS_BC_F(:)
        
         dr = 0.0000000001
            al_end=al0; al00 =al0
            bet_end=bet0;bet00=bet0
            K_F_end=K_F_start;kf00=K_F_start
            rv00=r0
            rv_end=r0
            
            v1 = v1 - sum(v1*c%norm_F(k_F_start,1:3))*c%norm_F(k_F_start,1:3)
            d  = sum(v1**2)**0.5
            dd = d
            
            ! %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            ! k_vertex=0
            sw=1
            i=k_F_start
            kk=0
            
         
        
        do while(sw.eq.1)    
 
            kk=kk+1
            v1 = v1 / (sum(v1**2))**0.5
            v1= dd*v1
            call WKONTRI_Convert(C,i,v1,al1,bet1)
            al3 = al0 + al1
            bet3 = bet0 +bet1
            call WKONTRI_POS(al3,bet3,k)
            k0=k
              
           
            if(k.ne.-1)then   
                sw=0    ! end point!
                k_f_end=i
                rv_end = r0 +v1
                al_end=al3
                bet_end=bet3
            else
                al1=al1/dd; bet1=bet1/dd 
                call WKONTRI_border(al0,bet0,k0,al1,bet1,x,k)
                if(k.ne.-1)then ! on the border
                    al2 = al0+x*al1; bet2 =bet0+x*bet1
                    r1(1:3) = c%rf(i,1,1:3)+al2*c%e(i,1,1:3) +bet2*c%e(i,2,1:3)
                    v1=v1/dd
                    dd = dd - (sum((r1-r0)**2))**0.5
                    
                    j= c%F_F(i,k)
                   
                 
                    
                    U = c%rv( c%v_E( abs(c%E_F(i,k)),2),1:3) -c%rv( c%V_E( abs(c%E_F(i,k)),1),1:3)
                    n1 = c%norm_F(i,1:3); n2=c%norm_F(j,1:3)
                    call WKONTRI_crossborder(v1,v2,n1,n2,u)
                    
                else  ! on border
                    ! k=-1 indicates that v1 pointing at some vertex
                    ! case 1
                    k=-1
                    if(abs(al0).lt.dr.and.abs(bet0).lt.dr)then
                        k=1;al0 = 0.000001;bet0=0.000001
                    endif
                    if(abs(al0).lt.dr.and.abs(bet0-1.0).lt.dr)then
                        k=3;al0 = 0.000001;bet0=1.0-0.000002
                    endif
                    if(abs(al0-1.0).lt.dr.and.abs(bet0).lt.dr)then
                        k=2;al0 = 1-0.000002;bet0=0.000001
                    endif

                    if(k.ne.-1)then ! not on vertex
                        r0  = C%rF(i,1,1:3) +al0*C%e(i,1,1:3) +bet0*C%e(i,2,1:3)
                        r1=r0;j=i;v2=v1

                    else
                        if(abs(al0+bet0-1).lt.dr.and.abs(al0).gt.dr.and.abs(al0-1.0).gt.dr)then
                            k=2; j= c%F_F(i,k)
                            U = c%rv( c%v_E( abs(c%E_F(i,k)),2),1:3) -c%rv( c%V_E( abs(c%E_F(i,k)),1),1:3)
                            n1 = c%norm_F(i,1:3); n2=c%norm_F(j,1:3)
                            v1=v1/dd
                            call WKONTRI_crossborder(v1,v2,n1,n2,u)
                        endif
                    
                        if(abs(al0).lt.dr.and.abs(bet0).gt.dr.and.abs(bet0-1.0).gt.dr)then
                            k=3; j= c%F_F(i,k)
                            U = c%rv( c%v_E( abs(c%E_F(i,k)),2),1:3) -c%rv( c%V_E( abs(c%E_F(i,k)),1),1:3)
                            n1 = c%norm_F(i,1:3); n2=c%norm_F(j,1:3)
                            v1=v1/dd
                            call WKONTRI_crossborder(v1,v2,n1,n2,u)
                        endif
                    
                        if(abs(bet0).lt.dr.and.abs(al0).gt.dr.and.abs(al0-1.0).gt.dr)then
                            k=1; j= c%F_F(i,k)
                            U = c%rv( c%v_E( abs(c%E_F(i,k)),2),1:3) -c%rv( c%V_E( abs(c%E_F(i,k)),1),1:3)
                            n1 = c%norm_F(i,1:3); n2=c%norm_F(j,1:3)
                            v1=v1/dd
                            call WKONTRI_crossborder(v1,v2,n1,n2,u)
                        endif
                    endif ! on whole border including vertex
                endif !    
                    
            endif
            
            if(sw.eq.1)then
 
                i = j
                r0 = r1
                if(i.ge.1.and.i.le.C%N_F)then
                if(IS_BC_F(i).eq.1)sw=-1
                  v0 = r0 - c%rf(i,1,:)
                  v1=v2
                call WKONTRI_Convert(C,i,v0,al0,bet0)
                else
                    
                    k_F_end=kf00
                    al_end=al00
                    bet_end=bet00
                    rv_end=rv00
                    sw = -1
                endif
                
                
                
            endif
            if(kk.gt.100)then
            sw=-1    
            endif
            
        enddo
    
    END subroutine WKONTRI_MOVE
        

    Subroutine WKONTRI_POS_TRI_POINT(r00,r1,r2,r3,n0,n1,n2,n3,r,d,al,bet,k)
     ! To determine the relative position of r with respect to
     ! the triangle (r1,r2,r3) with 3-norm vector (n1,n2,n3) and the norm vector n0.
     ! Return: k pos type of r  
     ! k=0, r is not under or on the top of (r123
     ! k=1, r is on the right top of r123
     ! k=-1, r is right under r123
     
        implicit none
        real*8 r(1:3)                             
        real*8 r1(1:3),r2(1:3),r3(1:3)
        real*8 n0(1:3),n1(1:3),n2(1:3),n3(1:3)
        real*8 rn1(1:3),rn2(1:3),rn3(1:3)
        
        real*8 al,bet,r00
        integer k_updn,k_inout,k
        
        real*8 d1,d2,d3,x1,x2,x3,d
        real*8 g11,g12,g22,g21,xx1,xx2,det
        real*8 e1(1:3),e2(1:3),v(1:3)
        
        
        d1= sum((r-r1)*n0)
        if(d1-r00.ge.0)k_updn =  1
        if(d1-r00.lt.0)k_updn = -1
        
        x1 = d1/sum(n1*n0)
        rn1 = r1+x1*n1
        x2 = d1/sum(n2*n0)
        rn2 = r2+x2*n2
        x3 = d1/sum(n3*n0)
        rn3 = r3+x3*n3
        
        e1=rn2-rn1
        e2=rn3-rn1
        v=r-rn1
        g11=sum(e1*e1);g21=sum(e1*e2);g12=g21;g22=sum(e2*e2)
        xx1 = sum(v*e1);xx2=sum(v*e2)
        det=(g11*g22-g21**2)

        al  = (xx1*g22-xx2*g21)/det
        
        bet  = (xx2*g11-xx1*g21)/det
        
        if(al.le.1.and.bet.le.1.and.al+bet.le.1.and.al.ge.0.and.bet.ge.0)then
            k_inout=1
        else
            k_inout=-1
        endif
        
        k = 0
        if(k_inout.eq.1.and.k_updn.eq.1)then
          k = 1
        endif
        
        if(k_inout.eq.1.and.k_updn.eq.-1)then
          k = -1
        endif
        
        
        d=d1-r00
        
        x1 = r00/sum(n1*n0)
        r1 = r1+x1*n1
        x2 = r00/sum(n2*n0)
        r2 = r2+x2*n2
        x3 = r00/sum(n3*n0)
        r3 = r3+x3*n3
        
        
        
        
    
    end subroutine WKONTRI_POS_TRI_POINT





 
end module Surf_operation_Mod
