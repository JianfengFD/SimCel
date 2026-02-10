    
module force_mod
use basicfun_mod
use CelType_mod
use point_mod
use Surf_Operation_mod
    contains
    
    subroutine Point_PV_Force_Cell(C,k,V)
        implicit none
        type(cel)::C(:)
        type(PtCell)::V
       integer,intent(in)::k
       
       !write(*,*)C(k)%RV(V%V_MOVE,:),V%V_move
       !pause
       V%PV_Force_Point=0
       if(C(k)%IS_OPEN.eq.0)then
       call POINT_PV_FORCE(V%v_move,V%Darea_Point,V%Dvol_point,V%PV_Force_Point,&
       C(k)%RV,C(k)%V_V,C(k)%F_V,& 
        C(k)%N_V_V,C(k)%Vector_F,C(k)%area_F,C(k)%P,C(k)%lamda)
       endif
    end subroutine Point_PV_Force_Cell

    subroutine Point_H_Force_Cell(C,k,V,Integral_H_A)
        implicit none
        type(cel)::C(:)
        type(PtCell)::V
       integer,intent(in)::k
       real*8 Integral_H_A
        
       if(C(k)%IS_BC_V(V%V_move).eq.0)then ! ensure V_move is not on the BC
           call Point_H_FORCE_numerical(V%V_move,V%PM_Force_Point,C(k)%RV,C(k)%H_modulus,C(k)%H_zero, C(k)%H_V,C(k)%V_V &
    ,C(k)%N_V_V, C(k)%area_V,C(k)%IS_BC_V)
          
       endif   
       
    end subroutine Point_H_Force_Cell

    subroutine Point_Rep_Force_Cell(C,K,Ntot,V,ep_LJ,sig_LJ,A_LJ,B_LJ)
        implicit none
        type(cel)::C(:)
        type(PtCell)::V
        integer,intent(in)::k    
        integer i,j,Ntot
        integer vmove, v2
        real*8 r1(1:3),r2(1:3),rvt(1:3),dd
        real*8 pot
        real*8 ep_LJ,sig_LJ,A_LJ,B_LJ
        
        !    ep_LJ = 0.01 *Cels(1)%H_modulus                    ! %%%  0.1 kpp0 =0.02 pN*um
!sig_LJ = 1.5* sum(Cels(1)%Len_E(1:Cels(1)%N_E))/Cels(1)%N_E
!A_LJ =0.4*ep_LJ* sig_LJ**12;B_LJ =0.4*ep_LJ* sig_LJ**6

        
            !Ntot = ubound(C,1)

            vmove = V%v_move
            r1 = V%rv_move
            V%Rep_Force_point =0.0
            do i = 1, Ntot
                if(i.ne.k)then
                do j = 1, C(k)%N_NearV(i,vmove)
                    v2 = C(k)%NearV(i,vmove,j)
                    r2 = C(i)%rv(v2,1:3)
                    rvt = r1 -r2
                    dd = sum(rvt**2)**0.5
                    Pot = 0
                    if(dd .lt.2.65*sig_LJ )then
                    pot = -ep_LJ* ( 12*A_LJ/dd**12.0 -6*B_LJ/dd**6.0 )
                    endif
                    V%Rep_Force_point=V%Rep_Force_point - pot * rvt !*C(i)%area_V(v2)/3.0
                enddo
                endif
            enddo
            
           ! write(*,*) Vmove,V%Rep_Force_point
           ! pause
        
        
    end subroutine Point_Rep_Force_Cell
    
    
    
    subroutine Point_Rep_Force_Cell_linear(C,K,Ntot,V,eta,R01)
        implicit none
        type(cel)::C(:)
        type(PtCell)::V
        integer,intent(in)::k    
        integer i,j,Ntot
        integer vmove, v2
        real*8 r1(1:3),r2(1:3),rvt(1:3),dd
        real*8 pot,eta,R01
        real*8 ep_LJ,sig_LJ,A_LJ,B_LJ
        
       
       
            vmove = V%v_move
            r1 = V%rv_move
            V%Rep_Force_point =0.0
            if( C(k)%IS_BC_V(vmove).eq.0)then
            do i = 1, Ntot
                if(i.ne.k)then
                do j = 1, C(k)%N_NearV(i,vmove)
                    v2 = C(k)%NearV(i,vmove,j)
                    r2 = C(i)%rv(v2,1:3)
                    rvt = r1 -r2
                    dd = sum(rvt**2)**0.5
                    Pot = 0
                    if(dd .lt.3.0.and.dd.gt.R01 )then
                    V%Rep_Force_point = V%Rep_Force_point - eta * rvt/dd**6
                    endif
                     !*C(i)%area_V(v2)/3.0
                enddo
                endif
            enddo
            endif
            
           ! write(*,*) Vmove,V%Rep_Force_point
           ! pause
        
        
    end subroutine Point_Rep_Force_Cell_linear
    
    
   
    subroutine Point_hard_rep_Cell(C,K,Ntot,V)
        implicit none
        type(cel)::C(:)
        type(PtCell)::V
        integer,intent(in)::k    
        integer i,j,Ntot, k_min(1:2),kk
        integer vmove, v2
        real*8 r1(1:3),r2(1:3),rvt(1:3),dd
        real*8 dmin
        
            Ntot = ubound(C,1)

            vmove = V%v_move
            r1 = C(k)%rv(vmove,1:3)
            V%Rep_Force_point =0.0
            k_min =0;kk=0
            dmin = 1000000
            do i = 1, Ntot
                do j = 1, C(k)%N_NearV(i,vmove)
                    v2 = C(k)%NearV(i,vmove,j)
                    r2 = C(i)%rv(v2,1:3)
                    rvt = r1 -r2
                    dd = sum(rvt**2)**0.5
                    if(dd .lt.2.25 )then
                        if(dd.lt.dmin)then
                            dmin = dd
                            k_min(1) = i; k_min(2) = v2
                        endif
                    endif
                enddo
            enddo
             
            if(dmin.lt.1.2)then
               ! C(k)%rv(vmove,1:3)=C(k_min(1))%rv(k_min(2),1:3) + 1.21*(r1-C(k_min(1))%rv(k_min(2),1:3)) /dmin**0.5    
            endif
        
        
    end subroutine Point_Hard_Rep_cell
   
    
    
    

    subroutine Point_Update_Cell(C,k,V)
        implicit none
        type(cel)::C(:)
        type(PtCell)::V
        integer,intent(in)::k
        if(C(k)%IS_BC_V(V%V_move).eq.0)then
          call Update_New_Parameters(V%V_move,V%RV_MOVE,C(k)%rv,C(k)%N_V_V,C(k)%V_V,C(k)%E_V,C(k)%F_V, &
            C(k)%F_E,C(k)%E_V_opposite,C(k)%V_E_opposite,C(k)%Len_E,C(k)%th_E,C(k)%vector_F, &
            C(k)%Norm_F,C(k)%Norm_V,C(k)%area_F, C(k)%area_V,C(k)%area_tot,C(k)%vol_F,C(k)%vol_total,C(k)%H_V)
        endif  
        
        
        
    end subroutine Point_Update_Cell
    
    
    
    subroutine POINT_PV_FORCE(v_move,Darea_Point,Dvol_point,PV_Force_Point, &
    Rv, V_V,F_V,N_V_V, Vector_F,area_F,P,lamda)
        implicit none

        real*8 rv(:,:),area_F(:)

        integer V_V(:,:),N_V_V(:),i,j
        integer V_MOVE,F_V(:,:)

        real*8 P,lamda,PV_Force_Point(1:3)

        real*8 R(1:20,1:3)
        real*8  Vect_VV(1:20,1:3),Dvol_Point(1:3)

        real*8 vector_F(:,:),at(1:3),Darea_Point(1:20,1:3)

        integer N_V_V_local


                 
            Dvol_Point=0.0
            PV_Force_Point=0.0

            ! get the variation

            N_V_V_local=N_V_V(v_move)
            R(N_V_V_local+1,1:3)=R(1,1:3)
            N_V_V_local=N_V_V(v_move)

            R(1:N_V_V_local,1:3)=Rv(V_V(v_move,1:N_V_V_local),1:3)
            R(N_V_V_local+1,1:3)=R(1,1:3)

            do j=1,N_V_V_local
                vect_VV(j,1:3)= vector_F(F_V(V_move,j),1:3)
            enddo

            do j=1,N_V_V_local
                Dvol_Point(1)=Dvol_Point(1)+R1R2_fun(R(j,1:3), R(j+1,1:3), 1)
                Dvol_Point(2)=Dvol_Point(2)+R1R2_fun(R(j,1:3), R(j+1,1:3), 2)
                Dvol_Point(3)=Dvol_Point(3)+R1R2_fun(R(j,1:3), R(j+1,1:3), 3)
                at(1)=( R1R2_fun(R(j,1:3), Vect_VV(j,1:3), 1) +R1R2_fun( Vect_VV(j,1:3),R(j+1,1:3), 1)  )/ 2.0
                at(2)=( R1R2_fun(R(j,1:3), Vect_VV(j,1:3), 2) +R1R2_fun( Vect_VV(j,1:3),R(j+1,1:3), 2)  )/ 2.0
                at(3)=( R1R2_fun(R(j,1:3), Vect_VV(j,1:3), 3) +R1R2_fun( Vect_VV(j,1:3),R(j+1,1:3), 3)  )/ 2.0
                Darea_point(j,1:3)=at
            enddo

            Dvol_Point=Dvol_Point/6.0

            ! %%%%%%%% pressure part

            PV_Force_Point = PV_Force_Point - Dvol_point*P 

            ! %%%%%%% area part
            do i=1,N_V_V_local
                PV_Force_Point = PV_Force_Point - lamda * Darea_Point(i,1:3)/area_F(F_V(v_move,i))
            enddo

        END subroutine POINT_PV_FORCE

    
    subroutine POINT_PV_FORCE_Numerical(v_move,PV_Force_Point,&
    Rv, V_V,F_V,N_V_V,Vector_F, P,lamda)
        implicit none

        real*8 rv(:,:)

        integer V_V(:,:),N_V_V(:),i,j
        integer V_MOVE,F_V(:,:)

        real*8 P,lamda,PV_Force_Point(1:3)

        real*8 R(1:20,1:3),RR(1:3,1:3)
        real*8 R1(1:3),R2(1:3),R3(1:3)
        real*8 det_r
        real*8 energyS1,energyS2,energyV1,energyV2
        real*8 s1,s2,v1,v2
        real*8 R0_old(1:3),R0_New(1:3)

        real*8  Vect_VV(1:20,1:3),Dvol_Point(1:3),vect1(1:3)

        real*8 vector_F(:,:),at(1:3),DS_Point(1:3)

        integer N_t




        Dvol_Point=0.0
        PV_Force_Point=0.0
        det_r=0.000001

        N_t=N_V_V(v_move)

        R(1:N_t,1:3)=Rv(V_V(v_move,1:N_t),1:3)
        R(N_t+1,1:3)=R(1,1:3)
        R0_old(1:3)=Rv(v_move,1:3)
        ! %%%%% x direction
        !    x energy 1
            R0_New(1:3)= R0_old(1:3);  R0_New(1)= R0_old(1)+det_r; 
            S2=0;V2=0
	        do i=1,N_t
		        RR(1,1:3)=R0_NEW(1:3); RR(2,1:3)= R(i,1:3);RR(3,1:3)= R(i+1,1:3)
		        V2=V2+V6_FUN(RR)/6.0
		        R1(1:3)=R0_NEW(1:3); R2(1:3)= R(i,1:3);R3(1:3)= R(i+1,1:3)
		        call Get_area_local(R1,R2,R3,S1, vect1)
		        S2=S2+S1
	        enddo 
	
	        energyS1= S2 * lamda
	        energyV1= V2 * P
        !    x energy 2
            R0_New(1:3)= R0_old(1:3);  R0_New(1)= R0_old(1)-det_r; 
            S2=0;V2=0
	        do i=1,N_t
		        RR(1,1:3)=R0_NEW(1:3); RR(2,1:3)= R(i,1:3);RR(3,1:3)= R(i+1,1:3)
		        V2=V2+V6_FUN(RR)/6.0
		        R1(1:3)=R0_NEW(1:3); R2(1:3)= R(i,1:3);R3(1:3)= R(i+1,1:3)
		        call Get_area_local(R1,R2,R3,S1, vect1)
		        S2=S2+S1
	        enddo 
	
	        energyS2= S2 * lamda
	        energyV2= V2 * P

           Dvol_Point(1)=(energyV1-energyV2)/2/det_r
           DS_Point(1)=(energyS1-energyS2)/2/det_r

        ! %%%%% y direction
        !    y energy 1
            R0_New(1:3)= R0_old(1:3);  R0_New(2)= R0_old(2)+det_r; 
            S2=0;V2=0
	        do i=1,N_t
		        RR(1,1:3)=R0_NEW(1:3); RR(2,1:3)= R(i,1:3);RR(3,1:3)= R(i+1,1:3)
		        V2=V2+V6_FUN(RR)/6.0
		        R1(1:3)=R0_NEW(1:3); R2(1:3)= R(i,1:3);R3(1:3)= R(i+1,1:3)
		        call Get_area_local(R1,R2,R3,S1, vect1)
		        S2=S2+S1
	        enddo 
	
	        energyS1= S2 * lamda
	        energyV1= V2 * P
        !    y energy 2
            R0_New(1:3)= R0_old(1:3);  R0_New(2)= R0_old(2)-det_r; 
            S2=0;V2=0
	        do i=1,N_t
		        RR(1,1:3)=R0_NEW(1:3); RR(2,1:3)= R(i,1:3);RR(3,1:3)= R(i+1,1:3)
		        V2=V2+V6_FUN(RR)/6.0
		        R1(1:3)=R0_NEW(1:3); R2(1:3)= R(i,1:3);R3(1:3)= R(i+1,1:3)
		        call Get_area_local(R1,R2,R3,S1, vect1)
		        S2=S2+S1
	        enddo 
	
	        energyS2= S2 * lamda
	        energyV2= V2 * P

           Dvol_Point(2)=(energyV1-energyV2)/2/det_r
           DS_Point(2)=(energyS1-energyS2)/2/det_r

        ! %%%%% z direction
        !    z energy 1
            R0_New(1:3)= R0_old(1:3);  R0_New(3)= R0_old(3)+det_r; 
            S2=0;V2=0
	        do i=1,N_t
		        RR(1,1:3)=R0_NEW(1:3); RR(2,1:3)= R(i,1:3);RR(3,1:3)= R(i+1,1:3)
		        V2=V2+V6_FUN(RR)/6.0
		        R1(1:3)=R0_NEW(1:3); R2(1:3)= R(i,1:3);R3(1:3)= R(i+1,1:3)
		        call Get_area_local(R1,R2,R3,S1, vect1)
		        S2=S2+S1
	        enddo 
	
	        energyS1= S2 * lamda
	        energyV1= V2 * P
        !    z energy 2
            R0_New(1:3)= R0_old(1:3);  R0_New(3)= R0_old(3)-det_r; 
            S2=0;V2=0
	        do i=1,N_t
		        RR(1,1:3)=R0_NEW(1:3); RR(2,1:3)= R(i,1:3);RR(3,1:3)= R(i+1,1:3)
		        V2=V2+V6_FUN(RR)/6.0
		        R1(1:3)=R0_NEW(1:3); R2(1:3)= R(i,1:3);R3(1:3)= R(i+1,1:3)
		        call Get_area_local(R1,R2,R3,S1, vect1)
		        S2=S2+S1
	        enddo 
	
	        energyS2= S2 * lamda
	        energyV2= V2 * P

           Dvol_Point(3)=(energyV1-energyV2)/2/det_r
           DS_Point(3)=(energyS1-energyS2)/2/det_r



        PV_Force_Point= -(Dvol_Point+DS_Point)

    END subroutine POINT_PV_FORCE_Numerical

    subroutine Point_H_FORCE_numerical(V_move,PM_Force_Point,RV,H_modulus,H_zero, H_V,V_V &
    ,N_V_V, area_V,IS_BC_V)
        implicit none

        real*8 rv(:,:)

        integer V_V(:,:),N_V_V(:)
        integer V_MOVE
        integer IS_BC_V(:)

        integer I_temp(1:20)

        real*8 H_V(:),area_V(:)
        real*8 H_zero,H_modulus

        real*8 area_V_local



        real*8 PM_force_Point(1:3)
        real*8 R(1:20,1:3)
        real*8 R0(1:3)

        real*8 RR(1:20,1:3)
        real*8 RR0(1:3)



        real*8 delt_r

        integer i,j,k,i1, N_t, I_temp2(1:20)


        real*8 sum0,sum_H1,sum_H2,H_local,energy_local_1,energy_local_2

        real*8 len_E_local(1:20),th_V_local(1:20)

        i=v_move

        PM_force_Point=0.0
        delt_r=0.000001

        i=v_move
        R=0
        I_temp=0

        I_temp(1:N_V_V(i) )=V_V(i, 1:N_V_V(i) )

        sum_H1=H_V(i)*area_V(i)/3.0
        do j=1,N_V_V(i)
            if(IS_BC_V(I_temp(j)).eq.0)then 
            R(j,1:3)=Rv(I_temp(j),1:3 )
            sum_H1=sum_H1+H_V(I_temp(j))*area_V(I_temp(j))/3.0
            endif
        enddo 


        ! %%%%%%%%%% Force of x direction
        ! engery_local_1
        R0(1:3)=Rv(i,1:3); 
        R0(1)=R0(1)+delt_r
        Energy_local_1=0.0
        Sum_H2=0.0


        N_t=N_V_V(i)

        call get_H_local_NEW(R0,R,N_t,   H_local, area_V_local,len_E_local,th_V_local)

        sum_H2=sum_H2+H_local*area_V_local/3.0
        Energy_local_1=Energy_local_1+H_modulus / 2.0 * (H_local**2-2*H_local*H_zero )* &
        &Area_V_local/3.0

        do j=1,N_V_V(i)
        if(IS_BC_V(I_temp(j)).eq.0)then 
            RR0(1:3)=R(j,1:3)
            i1=I_temp(j)  ! vertice No.
            I_temp2(1:N_V_V( i1 ) )=V_V(i1, 1:N_V_V( i1 ) )
            do k=1, N_V_V(i1)
                RR(k,1:3)=Rv(I_temp2(k),1:3)
                if(I_temp2(k).eq.i)RR(k,1:3)=R0
            enddo ! end of k

            N_t=N_V_V(i1)
            call get_H_local_NEW(RR0,RR,N_t,   H_local, area_V_local,len_E_local,th_V_local)
            sum_H2=sum_H2+H_local*area_V_local/3.0
            Energy_local_1=Energy_local_1+H_modulus / 2.0 * (H_local**2-2*H_local*H_zero )* &
            &Area_V_local/3.0
        endif    
        enddo  ! end  of j

        ! energy_local_2
        R0(1:3)=Rv(i,1:3); 
        R0(1)=R0(1)-delt_r
        Energy_local_2=0.0
        N_t=N_V_V(i)
        call get_H_local_NEW(R0,R,N_t,   H_local, area_V_local,len_E_local,th_V_local)


        Energy_local_2=Energy_local_2+H_modulus / 2.0 * (H_local**2-2*H_local*H_zero )*&
        & Area_V_local/3.0

        do j=1,N_V_V(i)
          if(IS_BC_V(I_temp(j)).eq.0)then    
            RR0(1:3)=R(j,1:3)
            i1=I_temp(j)  ! vertice No.
            I_temp2(1:N_V_V( i1 ) )=V_V(i1, 1:N_V_V( i1 ) )
            do k=1, N_V_V(i1)
                RR(k,1:3)=Rv(I_temp2(k),1:3)
                if(I_temp2(k).eq.i)RR(k,1:3)=R0
            enddo ! end of k

            N_t=N_V_V(i1)
            call get_H_local_NEW(RR0,RR,N_t,   H_local, area_V_local,len_E_local,th_V_local)
            Energy_local_2=Energy_local_2+H_modulus / 2.0 * (H_local**2-2*H_local*H_zero )* &
            &Area_V_local/3.0
        endif    
        enddo  ! end  of j


        PM_Force_Point(1)= - (Energy_Local_1-Energy_Local_2)/delt_r/2.0  

        ! %%%%%%%%%%%%%%

        ! %%%%%%%%%% Force of y direction
        ! engery_local_1
        R0(1:3)=Rv(i,1:3); 
        R0(2)=R0(2)+delt_r
        Energy_local_1=0.0


        N_t=N_V_V(i)

        call get_H_local_NEW(R0,R,N_t,   H_local, area_V_local,len_E_local,th_V_local)
        Energy_local_1=Energy_local_1+H_modulus / 2.0 * (H_local**2-2*H_local*H_zero )* &
        &Area_V_local/3.0

        do j=1,N_V_V(i)
         if(IS_BC_V(I_temp(j)).eq.0)then     
            RR0(1:3)=R(j,1:3)
            i1=I_temp(j)  ! vertice No.
            I_temp2(1:N_V_V( i1 ) )=V_V(i1, 1:N_V_V( i1 ) )
            do k=1, N_V_V(i1)
                RR(k,1:3)=Rv(I_temp2(k),1:3)
                if(I_temp2(k).eq.i)RR(k,1:3)=R0
            enddo ! end of k

            N_t=N_V_V(i1)
            call get_H_local_NEW(RR0,RR,N_t,   H_local, area_V_local,len_E_local,th_V_local)
            Energy_local_1=Energy_local_1+H_modulus / 2.0 * (H_local**2-2*H_local*H_zero )* &
            &Area_V_local/3.0
         endif   
        enddo  ! end  of j
        ! energy_local_2
        R0(1:3)=Rv(i,1:3); 
        R0(2)=R0(2)-delt_r
        Energy_local_2=0.0


        N_t=N_V_V(i)

        call get_H_local_NEW(R0,R,N_t,   H_local, area_V_local,len_E_local,th_V_local)
        Energy_local_2=Energy_local_2+H_modulus / 2.0 * (H_local**2-2*H_local*H_zero )*&
        & Area_V_local/3.0

        do j=1,N_V_V(i)
         if(IS_BC_V(I_temp(j)).eq.0)then     
            RR0(1:3)=R(j,1:3)
            i1=I_temp(j)  ! vertice No.
            I_temp2(1:N_V_V( i1 ) )=V_V(i1, 1:N_V_V( i1 ) )
            do k=1, N_V_V(i1)
                RR(k,1:3)=Rv(I_temp2(k),1:3)
                if(I_temp2(k).eq.i)RR(k,1:3)=R0
            enddo ! end of k

            N_t=N_V_V(i1)
            call get_H_local_NEW(RR0,RR,N_t,   H_local, area_V_local,len_E_local,th_V_local)
            Energy_local_2=Energy_local_2+H_modulus / 2.0 * (H_local**2-2*H_local*H_zero )* &
            &Area_V_local/3.0
         endif   
        enddo  ! end  of j


        PM_Force_Point(2)= - (Energy_Local_1-Energy_Local_2)/delt_r/2.0  

        !  %%%%%%%%%%%%%%%%%%%

        ! %%%%%%%%%% Force of z direction
        ! engery_local_1
        R0(1:3)=Rv(i,1:3); 
        R0(3)=R0(3)+delt_r
        Energy_local_1=0.0


        N_t=N_V_V(i)

        call get_H_local_NEW(R0,R,N_t,   H_local, area_V_local,len_E_local,th_V_local)
        Energy_local_1=Energy_local_1+H_modulus / 2.0 * (H_local**2-2*H_local*H_zero )*&
        & Area_V_local/3.0

        do j=1,N_V_V(i)
            if(IS_BC_V(I_temp(j)).eq.0)then 
            RR0(1:3)=R(j,1:3)
            i1=I_temp(j)  ! vertice No.
            I_temp2(1:N_V_V( i1 ) )=V_V(i1, 1:N_V_V( i1 ) )
            do k=1, N_V_V(i1)
                RR(k,1:3)=Rv(I_temp2(k),1:3)
                if(I_temp2(k).eq.i)RR(k,1:3)=R0
            enddo ! end of k

            N_t=N_V_V(i1)
            call get_H_local_NEW(RR0,RR,N_t,   H_local, area_V_local,len_E_local,th_V_local)
            Energy_local_1=Energy_local_1+H_modulus / 2.0 * (H_local**2-2*H_local*H_zero )* &
            &Area_V_local/3.0
            endif
        enddo  ! end  of j



        ! energy_local_2
        R0(1:3)=Rv(i,1:3); 
        R0(3)=R0(3)-delt_r
        Energy_local_2=0.0


        N_t=N_V_V(i)

        call get_H_local_NEW(R0,R,N_t,   H_local, area_V_local,len_E_local,th_V_local)
        Energy_local_2=Energy_local_2+H_modulus / 2.0 * (H_local**2-2*H_local*H_zero )*&
        & Area_V_local/3.0

        do j=1,N_V_V(i)
         if(IS_BC_V(I_temp(j)).eq.0)then 
            RR0(1:3)=R(j,1:3)
            i1=I_temp(j)  ! vertice No.
            I_temp2(1:N_V_V( i1 ) )=V_V(i1, 1:N_V_V( i1 ) )
            do k=1, N_V_V(i1)
                RR(k,1:3)=Rv(I_temp2(k),1:3)
                if(I_temp2(k).eq.i)RR(k,1:3)=R0
            enddo ! end of k

            N_t=N_V_V(i1)
            call get_H_local_NEW(RR0,RR,N_t,   H_local, area_V_local,len_E_local,th_V_local)

            Energy_local_2=Energy_local_2+H_modulus / 2.0 * (H_local**2-2*H_local*H_zero )* &
            &Area_V_local/3.0
        endif
        enddo  ! end  of j

        PM_Force_Point(3)= - (Energy_Local_1-Energy_Local_2)/delt_r/2.0  



    END subroutine Point_H_FORCE_numerical

    ! %%%%%%%%%

  


    ! %%%%%%%%%%%%%%%%%%

     
       
        Subroutine Update_NearF_Point_Rep(C,k1,Ntot,V,r00,r01,eta)  !including a hard repulsion
            implicit none
            type(CEL)::C(:)
            integer,intent(in)::k1,Ntot
            type(PtCell)::V
            integer i,j,k,i_cell
            integer i1,i2,i3,k2,k3,j1,j2,j3
            real*8 D_thres,DD,r00,r01,eta
            real*8 d1,d
            real*8 r(1:3)                             
            real*8 r1(1:3),r2(1:3),r3(1:3)
            real*8 n0(1:3),n1(1:3),n2(1:3),n3(1:3)
            real*8 al,bet
            real*8 d3,ad
            ad=1
           
            V%Rep_Force_point =0.0
            ! check if 
            do k2 = 1, Ntot
            if(k2.ne.k1)then
                do i = 1, C(k1)%N_NearF(k2,V%v_move)
                    j = C(k1)%NearF(k2,V%v_move,i)
                     i1=C(k2)%V_F(j,1);i2=C(k2)%V_F(j,2);i3=C(k2)%V_F(j,3)
                     
                     if(C(k2)%IS_BC_V(i1)+C(k2)%IS_BC_V(i2)+C(k2)%IS_BC_V(i3)+C(k1)%IS_BC_V(V%v_move).eq.0)then
                       r1 = C(k2)%RV(i1,:);r2 = C(k2)%RV(i2,:);r3 = C(k2)%RV(i3,:);
                       n0 = C(k2)%Norm_F(j,:);
                       n1 = C(k2)%norm_V(i1,:);n2 = C(k2)%norm_V(i2,:);n3 = C(k2)%norm_V(i3,:);
                       r  = C(k1)%RV(V%v_move,1:3)
                       call WKONTRI_POS_TRI_POINT(r00,r1,r2,r3,n0,n1,n2,n3,r,d,al,bet,k3)
                       ! d (r00 has been substracted
                    
                      if(k3.eq.1.and.d+r00.lt.r01)then ! r01=2.0
                          d3=d+r00
                          V%Rep_Force_point = V%Rep_Force_point + 12*C(k1)%eps*(2*C(k1)%sig_eps**6/d3**7 &
                              -ad*C(k1)%sig_eps**3/d3**4)*n0
                      endif
                      
                      if(k3.eq.-1)then
                          C(k1)%RV(V%v_move,1:3) = al*r2+bet*r3+(1-al-bet)*r1 +n0*0.01  
                      endif
                      if(k3.eq.0)then
                          do j1=1,C(k2)%N_F_F_Around(j)
                            j2 = C(k2)%F_F_Around(j,j1)
                            i1=C(k2)%V_F(j2,1);i2=C(k2)%V_F(j2,2);i3=C(k2)%V_F(j2,3);
                            r1 = C(k2)%RV(i1,:);r2 = C(k2)%RV(i2,:);r3 = C(k2)%RV(i3,:);
                            n0 = C(k2)%Norm_F(j,:);
                            n1 = C(k2)%norm_V(i1,:);n2 = C(k2)%norm_V(i2,:);n3 = C(k2)%norm_V(i3,:);
                            r  = C(k1)%RV(V%v_move,1:3)
                            call WKONTRI_POS_TRI_POINT(r00,r1,r2,r3,n0,n1,n2,n3,r,d,al,bet,k3)
                            if(k3.eq.1.and.d+r00.lt.r01)then
                                C(k1)%NearF(k2,v%V_move,i) = j2
                                d3=d+r00
                                V%Rep_Force_point = V%Rep_Force_point + 12*C(k1)%eps*(2*C(k1)%sig_eps**6/d3**7 &
                                    -ad*C(k1)%sig_eps**3/d3**4)*n0
                            endif
                            if(k3.eq.-1)then
                                C(k1)%NearF(k2,v%V_move,i) = j2
                                C(k1)%RV(V%v_move,1:3) = al*r2+bet*r3+(1-al-bet)*r1 +n0*0.01
                            endif
                          enddo
                      endif
                      
                  endif    
                      
                      
                Enddo
            endif
            enddo
            
            
            
            
        
        
        END subroutine Update_NearF_Point_Rep
        
     
     
     
       subroutine Get_NearF_Cells(C,k1,k2,D_thres,R00,contf)
            implicit none
            type(Cel)::C(:)
            integer,intent(in):: k1,k2
            integer k3
            integer i,j,k,i1,i2,i3
            real*8 dsq
            real*8 D_thres,DD,R00
            real*8 d1,d,d3
             real*8 r(1:3)                             
        real*8 r1(1:3),r2(1:3),r3(1:3)
        real*8 n0(1:3),n1(1:3),n2(1:3),n3(1:3)
            real*8 al,bet
            integer NN1,NF(1:4000)
            real*8 contf
            real*8 contarea
            integer IS_CONTACT
            
            integer IS_contact_V(1:4000)
            real*8 d_contact_V(1:4000)
            
            real*8 energy_ad
            real*8 d_contact
            contarea = 0.0
            DD = D_thres**2
            
            IS_Contact_V=0
            d_contact_V=-1.0
            
            C(k1)%energy_rep = 0.0
            if(k1.ne.k2)then
                do i = 1,C(k1)%N_V
                   C(k1)%N_NearF(k2,i) = 0
                   NN1 = 0
                   IS_CONTACT=0
                   do j=1,C(k2)%N_F
                    i1=C(k2)%V_F(j,1);i2=C(k2)%V_F(j,2);i3=C(k2)%V_F(j,3);
                     dsq =sum( (C(k1)%RV(i,:)-(C(k2)%RV(i1,:)+C(k2)%RV(i2,:)+C(k2)%RV(i3,:))/3.0)**2.0)
                     if(dsq.lt.DD*2.0)then
                         NN1=NN1+1
                         NF(NN1)=j  ! distance to the face's center point smaller than DD
                     endif
                   enddo
                   
                   do k=1,NN1
                       j=NF(k)
                       i1=C(k2)%V_F(j,1);i2=C(k2)%V_F(j,2);i3=C(k2)%V_F(j,3);
                       r1 = C(k2)%RV(i1,:);r2 = C(k2)%RV(i2,:);r3 = C(k2)%RV(i3,:);
                       n0 = C(k2)%Norm_F(j,:);
                       n1 = C(k2)%norm_V(i1,:);n2 = C(k2)%norm_V(i2,:);n3 = C(k2)%norm_V(i3,:);
                       r  = C(k1)%RV(i,1:3)
                       call WKONTRI_POS_TRI_POINT(r00,r1,r2,r3,n0,n1,n2,n3,r,d,al,bet,k3)
                       if(k3.eq.1)then
                           C(k1)%N_NearF(k2,i)=C(k1)%N_NearF(k2,i)+1
                           C(k1)%NearF(k2,i,C(k1)%N_NearF(k2,i)) = j
                           C(k1)%NearF_k(k2,i,C(k1)%N_NearF(k2,i)) = k3
                           if(d.lt.0)then
                                 C(k1)%RV(i,1:3) = al*r2+bet*r3+(1-al-bet)*r1 +n0*0.01
                           endif
                            if(d+r00.lt.r00*4)then
                                 IS_CONTACT=1
                                 d3=d+r00
                                 energy_ad = 4*C(k1)%eps*( (C(k1)%sig_eps/d3)**6 -(C(k1)%sig_eps/d3)**3  )
                                 if(energy_ad.lt.0.0)then
                                     C(k1)%energy_rep=C(k1)%energy_rep+energy_ad
                                 endif
                                 
                            endif
                            
                            if(d+r00.lt.r00*5.0)then
                                 IS_CONTACT=1
                                 d_contact=abs(d+r00)
                            endif
                       endif
                       if(abs(d).lt.r00.and.k3.eq.-1)then
                           C(k1)%N_NearF(k2,i)=C(k1)%N_NearF(k2,i)+1
                           C(k1)%NearF(k2,i,C(k1)%N_NearF(k2,i)) = j
                           C(k1)%NearF_k(k2,i,C(k1)%N_NearF(k2,i)) = 1
                           C(k1)%RV(i,1:3) = al*r2+bet*r3+(1-al-bet)*r1 +n0*0.01
                           IS_CONTACT = 1
                           d_contact=abs(d+r00)
                       endif
                       
                       
                     
                       
                   enddo
                   if(IS_CONTACT.eq.1)then
                     
                       if(d_contact.lt.r00+0.35)then
                           IS_contact_V(i) =1
                           contarea=contarea + C(k1)%area_V(i)/3.0
                       endif 
                       !
                       !if(d_contact.ge.r00+0.1)contarea=contarea + C(k1)%area_V(i)/3.0*(1-tanh(8*d_contact-5.4))/2.0
                   endif
                   
                enddo
                
            
                
                
                contf = contarea/C(2)%area_tot
                C(k1)%contarea=contarea
            endif
            
            
            
            
          end subroutine Get_NearF_Cells
        
     
     
     
    Subroutine Update_New_Parameters(v_move,rv_move,rv,N_V_V,V_V,E_V,F_V,F_E,E_V_opposite,V_E_opposite,&
    & Len_E,th_E,vector_F,Norm_F,Norm_V,area_F,area_V,area_tot,vol_F,vol_total,H_V)

        implicit none
        integer N_V_V(:)
        integer V_V(:,:),E_V(:,:),F_V(:,:),F_E(:,:)
        integer E_V_opposite(:,:),V_E_opposite(:,:)

        ! For cyclic or label
        integer i,j,k,i1,i2,i3,I_temp(1:40),I_temp2(1:40)
        integer N_t

        ! positions
        real*8 rv(:,:),R0(1:3)
        real*8 R(1:40,1:3),RR3(1:3,1:3),R1(1:3),R2(1:3),R3(1:3)
        real*8 RJ(1:40,1:3)
        ! edge variables
        real*8 len_E(:),th_E(:),th_JK,len_JK(1:20),sign_r
        !area ones
        real*8 vector_F(:,:),Vector_R0(1:40,1:3),Vector_RJ(1:40,1:3),vector_Local(1:3)
        real*8 Norm_F(:,:),Norm_V(:,:)
        real*8 area_F(:),area_V(:),area_tot,area_local,area_V_local
        real*8 sum_area1_point,sum_area2_point
        ! volumes
        real*8 Vol_F(:),Vol_total, sum_vol1_point,sum_vol2_point
        ! curvatures
        real*8 H_V(:),H_local,H_JK
        
        real*8 rv_move(1:3)
        integer V_move





        i=V_move

        N_t=N_V_V(i)
        sum_area2_point=sum(area_F(F_V(i,1:N_t)) )
        sum_vol2_point=sum(vol_F( F_V(i,1:N_t)) )
        R0=RV(V_move,1:3) 
        do j=1,N_t
            R(j,1:3)=RV(V_V(i,j),1:3)
            I_temp(j)=V_V(i,j)
        enddo

        ! %%%%%%%%%%%%%
        !  Get all quantities near vertex V
        ! %%%%%%%%%%%%%

        R(N_t+1,1:3)=R(1,1:3)
        R(N_t+2,1:3)=R(2,1:3)
        area_V_local=0.0
        ! update L_E area vol
        Norm_V(i,1:3) = 0
        do j=1,N_t
            R1(1:3)=R0(1:3);R2(1:3)=R(j,1:3);R3(1:3)=R(j+1,1:3);
            call Get_area_local(R1,R2,R3,area_local, vector_local)
            area_F(F_V(i,j))=area_Local
            vector_R0(j,1:3)= vector_local(1:3)/area_Local
            Vector_F(F_V(i,j),1:3)=Vector_local
            Norm_F(F_V(i,j),1:3) = Vector_F(F_V(i,j),1:3)/sum(Vector_F(F_V(i,j),1:3)**2)**0.5
            Norm_V(i,1:3)=Norm_V(i,1:3) +Vector_local
            area_V_local=area_V_local+area_local
            
            RR3(1,1:3)=r1(1:3)
            RR3(2,1:3)=r2(1:3)
            RR3(3,1:3)=r3(1:3)
            Vol_F(F_V(i,j))=v6_fun(RR3)/6.0

            len_jk(j)=(sum((R1-R2)**2.0))**0.5
            Len_E(abs(E_V(i,j)))=len_JK(j)
            i1=F_E(abs(E_V_opposite(i,j)),1);i2=F_E(abs(E_V_opposite(i,j)),2);
            if(i1.ne.F_V(i,j))i3=i1
            if(i2.ne.F_V(i,j))i3=i2
            Vector_RJ(j,1:3)=Vector_F(i3,1:3)/area_F(i3)
            i1=V_E_opposite(abs(E_V_opposite(i,j)),1);i2=V_E_opposite(abs(E_V_opposite(i,j)),2);
            if(i1.ne.i)RJ(j,1:3)=Rv(i1,1:3)
            if(i2.ne.i)RJ(j,1:3)=Rv(i2,1:3)
        enddo
        Norm_V(i,1:3) = Norm_V(i,1:3)/sum(Norm_V(i,1:3)**2)**0.5 

        vector_r0(N_t+1,1:3)=vector_R0(1,1:3)
        len_Jk(N_t+1)=len_JK(1)
        H_JK=0.0

        do j=1,N_t
        ! th_E E_V
            sign_r=sum(vector_R0(j,1:3)*vector_R0(j+1,1:3))
            if(sign_r.lt.1.0)then
                th_jk=dacos(sign_r)
            else
                th_JK=0.0
            endif
            i1=mod(j,N_t)+1
            sign_r=sum( (vector_R0(j,1:3)-vector_R0(j+1,1:3))*(R(j,1:3)-R(j+2,1:3)) )
            if(sign_r.lt.0.0)th_jk=-th_JK
            th_E(abs(E_V(i,i1)))=th_JK
            H_jk=H_jk+th_jk*len_jk(i1)

            ! Th_E E_V_opposite
            sign_r=sum(vector_R0(j,1:3)*vector_RJ(j,1:3))

            if(sign_r.lt.1.0)then
                th_jk=dacos(sign_r)
            else
                th_JK=0.0
            endif
            sign_r=sum( (vector_R0(j,1:3)-vector_RJ(j,1:3))*(RV_move(1:3)-RJ(j,1:3)) )
            if(sign_r.lt.0.0)th_jk=-th_JK
            th_E(abs(E_V_opposite(i,j)))=th_JK

        enddo


        H_local=-H_JK/(area_V_local/3.0)/2.0


        H_V(i)=H_local
        area_V(i)=area_V_local


        ! %%%%%%%%%%
        ! near VJ
        ! %%%%%%%

        do j=1,N_t

            i1=I_temp(j)  ! vertice No.
            I_temp2(1:N_V_V( i1 ) )=V_V(i1, 1:N_V_V( i1 ) )
            area_V_local=0.0
            H_jk=0.0
            Norm_V(i1,1:3)=0
            do k=1, N_V_V(i1)
                i2=abs(E_V(i1,k))
                area_V_local=area_V_local+area_F(F_V(i1,k))
                H_jk=H_jk+th_E(i2)*len_E(i2)
                Norm_V(i1,1:3) = Norm_V(i1,1:3)  + Vector_F(F_V(i1,k),1:3) 
            enddo ! end of k
            H_local=-H_JK/(area_V_local/3.0)/2.0
            Norm_V(i1,1:3) = Norm_V(i1,1:3)/sum(Norm_V(i1,1:3)**2)**0.5 
            area_V(i1)=area_V_local
            H_V(i1)=H_local
        enddo

        sum_area1_point=sum(area_F(F_V(i,1:N_t)))  !updating total area
        sum_vol1_point=sum(vol_F(F_V(i,1:N_t)))   !updating total volume

        area_tot=area_tot+sum_area1_point-sum_area2_point
        vol_total=vol_total+sum_vol1_point-sum_vol2_point

    END Subroutine Update_New_Parameters

end module force_mod
