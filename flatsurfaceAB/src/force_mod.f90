
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

       V%PV_Force_Point=0
       if(C(k)%IS_OPEN.eq.0)then
       call POINT_PV_FORCE(V%v_move,V%Darea_Point,V%Dvol_point,V%PV_Force_Point,&
       C(k)%RV,C(k)%V_V,C(k)%F_V,&
        C(k)%N_V_V,C(k)%Vector_F,C(k)%area_F,C(k)%P,C(k)%lamda)
       endif
    end subroutine Point_PV_Force_Cell

! ========================================================
! Bending force with phi coupling
! E_v = (kpp + kpp1*fi(v))/2 * (H_v + c0 - c1*fi(v))^2 * S_v/3
! ========================================================
    subroutine Point_H_Force_Cell(C,k,V,Integral_H_A)
        implicit none
        type(cel)::C(:)
        type(PtCell)::V
       integer,intent(in)::k
       real*8 Integral_H_A

       if(C(k)%IS_BC_V(V%V_move).eq.0)then
           call Point_H_Force_numerical_phi(V%V_move,V%PM_Force_Point,C(k)%RV, &
    C(k)%H_modulus, C(k)%fi, C(k)%kpp1_phi, C(k)%c0_phi, C(k)%c1_phi, &
    C(k)%H_V, C(k)%V_V, C(k)%N_V_V, C(k)%area_V, C(k)%IS_BC_V)
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
                    V%Rep_Force_point=V%Rep_Force_point - pot * rvt
                enddo
                endif
            enddo

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
                enddo
                endif
            enddo
            endif

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

! ========================================================
! Bending force with phi coupling (numerical finite difference)
! Energy_v = (kpp + kpp1*fi(v))/2 * (H_v + c0 - c1*fi(v))^2 * S_v/3
! ========================================================
    subroutine Point_H_Force_numerical_phi(V_move,PM_Force_Point,RV, &
    & H_modulus, fi, kpp1_phi, c0_phi, c1_phi, H_V,V_V,N_V_V, area_V,IS_BC_V)
        implicit none

        real*8 rv(:,:)
        integer V_V(:,:),N_V_V(:)
        integer V_MOVE
        integer IS_BC_V(:)
        integer I_temp(1:20)

        real*8 H_V(:),area_V(:),fi(:)
        real*8 H_modulus, kpp1_phi, c0_phi, c1_phi
        real*8 kpp_local, H0_local

        real*8 area_V_local
        real*8 PM_force_Point(1:3)
        real*8 R(1:20,1:3)
        real*8 R0(1:3)
        real*8 RR(1:20,1:3)
        real*8 RR0(1:3)
        real*8 delt_r

        integer i,j,k,i1, N_t, I_temp2(1:20)
        real*8 sum_H1,sum_H2,H_local,energy_local_1,energy_local_2
        real*8 len_E_local(1:20),th_V_local(1:20)

        i=v_move
        PM_force_Point=0.0
        delt_r=0.000001

        R=0
        I_temp=0
        I_temp(1:N_V_V(i))=V_V(i, 1:N_V_V(i))

        do j=1,N_V_V(i)
            if(IS_BC_V(I_temp(j)).eq.0)then
            R(j,1:3)=Rv(I_temp(j),1:3)
            endif
        enddo

        ! Center vertex phi-dependent parameters
        kpp_local = H_modulus + kpp1_phi * fi(i)
        H0_local = -c0_phi + c1_phi * fi(i)

        ! %%%%%%%%%% Force of x direction
        ! energy_local_1
        R0(1:3)=Rv(i,1:3)
        R0(1)=R0(1)+delt_r
        Energy_local_1=0.0

        N_t=N_V_V(i)
        call get_H_local_NEW(R0,R,N_t, H_local, area_V_local,len_E_local,th_V_local)
        Energy_local_1=Energy_local_1+kpp_local/2.0*(H_local-H0_local)**2*Area_V_local/3.0

        do j=1,N_V_V(i)
        if(IS_BC_V(I_temp(j)).eq.0)then
            RR0(1:3)=R(j,1:3)
            i1=I_temp(j)
            I_temp2(1:N_V_V(i1))=V_V(i1, 1:N_V_V(i1))
            do k=1, N_V_V(i1)
                RR(k,1:3)=Rv(I_temp2(k),1:3)
                if(I_temp2(k).eq.i)RR(k,1:3)=R0
            enddo
            N_t=N_V_V(i1)
            call get_H_local_NEW(RR0,RR,N_t, H_local, area_V_local,len_E_local,th_V_local)
            Energy_local_1=Energy_local_1+(H_modulus+kpp1_phi*fi(i1))/2.0* &
            &(H_local-(-c0_phi+c1_phi*fi(i1)))**2*Area_V_local/3.0
        endif
        enddo

        ! energy_local_2
        R0(1:3)=Rv(i,1:3)
        R0(1)=R0(1)-delt_r
        Energy_local_2=0.0
        N_t=N_V_V(i)
        call get_H_local_NEW(R0,R,N_t, H_local, area_V_local,len_E_local,th_V_local)
        Energy_local_2=Energy_local_2+kpp_local/2.0*(H_local-H0_local)**2*Area_V_local/3.0

        do j=1,N_V_V(i)
          if(IS_BC_V(I_temp(j)).eq.0)then
            RR0(1:3)=R(j,1:3)
            i1=I_temp(j)
            I_temp2(1:N_V_V(i1))=V_V(i1, 1:N_V_V(i1))
            do k=1, N_V_V(i1)
                RR(k,1:3)=Rv(I_temp2(k),1:3)
                if(I_temp2(k).eq.i)RR(k,1:3)=R0
            enddo
            N_t=N_V_V(i1)
            call get_H_local_NEW(RR0,RR,N_t, H_local, area_V_local,len_E_local,th_V_local)
            Energy_local_2=Energy_local_2+(H_modulus+kpp1_phi*fi(i1))/2.0* &
            &(H_local-(-c0_phi+c1_phi*fi(i1)))**2*Area_V_local/3.0
        endif
        enddo

        PM_Force_Point(1)= - (Energy_Local_1-Energy_Local_2)/delt_r/2.0

        ! %%%%%%%%%% Force of y direction
        ! energy_local_1
        R0(1:3)=Rv(i,1:3)
        R0(2)=R0(2)+delt_r
        Energy_local_1=0.0

        N_t=N_V_V(i)
        call get_H_local_NEW(R0,R,N_t, H_local, area_V_local,len_E_local,th_V_local)
        Energy_local_1=Energy_local_1+kpp_local/2.0*(H_local-H0_local)**2*Area_V_local/3.0

        do j=1,N_V_V(i)
         if(IS_BC_V(I_temp(j)).eq.0)then
            RR0(1:3)=R(j,1:3)
            i1=I_temp(j)
            I_temp2(1:N_V_V(i1))=V_V(i1, 1:N_V_V(i1))
            do k=1, N_V_V(i1)
                RR(k,1:3)=Rv(I_temp2(k),1:3)
                if(I_temp2(k).eq.i)RR(k,1:3)=R0
            enddo
            N_t=N_V_V(i1)
            call get_H_local_NEW(RR0,RR,N_t, H_local, area_V_local,len_E_local,th_V_local)
            Energy_local_1=Energy_local_1+(H_modulus+kpp1_phi*fi(i1))/2.0* &
            &(H_local-(-c0_phi+c1_phi*fi(i1)))**2*Area_V_local/3.0
         endif
        enddo

        ! energy_local_2
        R0(1:3)=Rv(i,1:3)
        R0(2)=R0(2)-delt_r
        Energy_local_2=0.0

        N_t=N_V_V(i)
        call get_H_local_NEW(R0,R,N_t, H_local, area_V_local,len_E_local,th_V_local)
        Energy_local_2=Energy_local_2+kpp_local/2.0*(H_local-H0_local)**2*Area_V_local/3.0

        do j=1,N_V_V(i)
         if(IS_BC_V(I_temp(j)).eq.0)then
            RR0(1:3)=R(j,1:3)
            i1=I_temp(j)
            I_temp2(1:N_V_V(i1))=V_V(i1, 1:N_V_V(i1))
            do k=1, N_V_V(i1)
                RR(k,1:3)=Rv(I_temp2(k),1:3)
                if(I_temp2(k).eq.i)RR(k,1:3)=R0
            enddo
            N_t=N_V_V(i1)
            call get_H_local_NEW(RR0,RR,N_t, H_local, area_V_local,len_E_local,th_V_local)
            Energy_local_2=Energy_local_2+(H_modulus+kpp1_phi*fi(i1))/2.0* &
            &(H_local-(-c0_phi+c1_phi*fi(i1)))**2*Area_V_local/3.0
         endif
        enddo

        PM_Force_Point(2)= - (Energy_Local_1-Energy_Local_2)/delt_r/2.0

        ! %%%%%%%%%% Force of z direction
        ! energy_local_1
        R0(1:3)=Rv(i,1:3)
        R0(3)=R0(3)+delt_r
        Energy_local_1=0.0

        N_t=N_V_V(i)
        call get_H_local_NEW(R0,R,N_t, H_local, area_V_local,len_E_local,th_V_local)
        Energy_local_1=Energy_local_1+kpp_local/2.0*(H_local-H0_local)**2*Area_V_local/3.0

        do j=1,N_V_V(i)
            if(IS_BC_V(I_temp(j)).eq.0)then
            RR0(1:3)=R(j,1:3)
            i1=I_temp(j)
            I_temp2(1:N_V_V(i1))=V_V(i1, 1:N_V_V(i1))
            do k=1, N_V_V(i1)
                RR(k,1:3)=Rv(I_temp2(k),1:3)
                if(I_temp2(k).eq.i)RR(k,1:3)=R0
            enddo
            N_t=N_V_V(i1)
            call get_H_local_NEW(RR0,RR,N_t, H_local, area_V_local,len_E_local,th_V_local)
            Energy_local_1=Energy_local_1+(H_modulus+kpp1_phi*fi(i1))/2.0* &
            &(H_local-(-c0_phi+c1_phi*fi(i1)))**2*Area_V_local/3.0
            endif
        enddo

        ! energy_local_2
        R0(1:3)=Rv(i,1:3)
        R0(3)=R0(3)-delt_r
        Energy_local_2=0.0

        N_t=N_V_V(i)
        call get_H_local_NEW(R0,R,N_t, H_local, area_V_local,len_E_local,th_V_local)
        Energy_local_2=Energy_local_2+kpp_local/2.0*(H_local-H0_local)**2*Area_V_local/3.0

        do j=1,N_V_V(i)
         if(IS_BC_V(I_temp(j)).eq.0)then
            RR0(1:3)=R(j,1:3)
            i1=I_temp(j)
            I_temp2(1:N_V_V(i1))=V_V(i1, 1:N_V_V(i1))
            do k=1, N_V_V(i1)
                RR(k,1:3)=Rv(I_temp2(k),1:3)
                if(I_temp2(k).eq.i)RR(k,1:3)=R0
            enddo
            N_t=N_V_V(i1)
            call get_H_local_NEW(RR0,RR,N_t, H_local, area_V_local,len_E_local,th_V_local)
            Energy_local_2=Energy_local_2+(H_modulus+kpp1_phi*fi(i1))/2.0* &
            &(H_local-(-c0_phi+c1_phi*fi(i1)))**2*Area_V_local/3.0
        endif
        enddo

        PM_Force_Point(3)= - (Energy_Local_1-Energy_Local_2)/delt_r/2.0

    END subroutine Point_H_Force_numerical_phi

    ! %%%%%%%%%

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


        ! near VJ
        do j=1,N_t

            i1=I_temp(j)
            I_temp2(1:N_V_V( i1 ) )=V_V(i1, 1:N_V_V( i1 ) )
            area_V_local=0.0
            H_jk=0.0
            Norm_V(i1,1:3)=0
            do k=1, N_V_V(i1)
                i2=abs(E_V(i1,k))
                area_V_local=area_V_local+area_F(F_V(i1,k))
                H_jk=H_jk+th_E(i2)*len_E(i2)
                Norm_V(i1,1:3) = Norm_V(i1,1:3)  + Vector_F(F_V(i1,k),1:3)
            enddo
            H_local=-H_JK/(area_V_local/3.0)/2.0
            Norm_V(i1,1:3) = Norm_V(i1,1:3)/sum(Norm_V(i1,1:3)**2)**0.5
            area_V(i1)=area_V_local
            H_V(i1)=H_local
        enddo

        sum_area1_point=sum(area_F(F_V(i,1:N_t)))
        sum_vol1_point=sum(vol_F(F_V(i,1:N_t)))

        area_tot=area_tot+sum_area1_point-sum_area2_point
        vol_total=vol_total+sum_vol1_point-sum_vol2_point

    END Subroutine Update_New_Parameters

! ========================================================
! Shear (MS) Force
! ========================================================

    subroutine Point_MS_Force_Cell(C,k,V)
        implicit none
        type(cel)::C(:)
        type(PtCell)::V
        integer,intent(in)::k

        V%MS_Force_Point = 0.0
        if(C(k)%HAS_REFERENCE.eq.1)then
            call Point_MS_force(C(k)%rv, V%V_move, C(k)%F_V, C(k)%E_F, C(k)%V_E, &
            & C(k)%N_V_V, C(k)%kpp_alpha, C(k)%mu_ms, C(k)%a3_ms, C(k)%a4_ms, &
            & C(k)%b1_ms, C(k)%b2_ms, C(k)%Area_F_zero, C(k)%Len_E_zero, C(k)%Len_E, &
            & V%MS_Force_Point)
        endif

    end subroutine Point_MS_Force_Cell

    subroutine Point_MS_force(rv, v_move, F_V, E_F, V_E, N_V_V, &
    & kpp_alpha, mu_ms, a3_ms, a4_ms, b1_ms, b2_ms, &
    & area_F_zero, Len_E_zero, Len_E, MS_Force_Point)
        implicit none
        real*8,intent(in):: rv(:,:)
        integer,intent(in):: v_move
        integer,intent(in):: F_V(:,:), E_F(:,:), V_E(:,:), N_V_V(:)
        real*8,intent(in):: kpp_alpha, mu_ms, a3_ms, a4_ms, b1_ms, b2_ms
        real*8,intent(in):: area_F_zero(:), Len_E_zero(:), Len_E(:)
        real*8,intent(out):: MS_Force_Point(1:3)

        real*8 a,b,c,a0,b0,c0,alph,bet,S,S0
        real*8 r0(1:3),r1(1:3),r2(1:3)
        real*8 DFbet,DFalph
        real*8 DS(1:3),DPI(1:3)
        integer i,j,k,i3,k_c,k_a,k_b
        integer num_F, I_side(1:3)

        MS_Force_Point = 0.0
        r0 = rv(v_move,1:3)

        do j=1,N_V_V(v_move)
            num_F = F_V(v_move,j)
            I_side(1:3) = abs(E_F(num_F,1:3))
            k_c = 0
            do k=1,3
                if(V_E(I_side(k),1).ne.v_move .and. V_E(I_side(k),2).ne.v_move) k_c=k
            enddo
            if(k_c.eq.0) cycle

            i3=0
            k_a=0; k_b=0
            do k=1,3
                if(k.ne.k_c .and. i3.eq.1)then
                    k_b=k
                    if(V_E(I_side(k),1).ne.v_move) r2=rv(V_E(I_side(k),1),1:3)
                    if(V_E(I_side(k),2).ne.v_move) r2=rv(V_E(I_side(k),2),1:3)
                    i3=i3+1
                endif
                if(k.ne.k_c .and. i3.eq.0)then
                    k_a=k
                    if(V_E(I_side(k),1).ne.v_move) r1=rv(V_E(I_side(k),1),1:3)
                    if(V_E(I_side(k),2).ne.v_move) r1=rv(V_E(I_side(k),2),1:3)
                    i3=i3+1
                endif
            enddo

            a0=Len_E_zero(I_side(k_a))
            b0=Len_E_zero(I_side(k_b))
            c0=Len_E_zero(I_side(k_c))
            a =Len_E(I_side(k_a))
            b =Len_E(I_side(k_b))
            c =Len_E(I_side(k_c))

            call Get_alph_bet_analytic(a,b,c,a0,b0,c0,alph,bet,S,S0)

            DFalph = (kpp_alpha/2.0*(2.0*alph+3.0*a3_ms*alph**2+4.0*a4_ms*alph**3) &
                & + mu_ms*b2_ms*bet) * area_F_zero(num_F)
            DFbet = mu_ms*(1.0+b1_ms*alph+2*b2_ms*bet) * area_F_zero(num_F)

            DS = -4*a**2*(r0-r1) - 4*b**2*(r0-r2) + 4*a**2*(r0-r2) &
                & + 4*b**2*(r0-r1) + 4*c**2*(r0-r1) + 4*c**2*(r0-r2)

            DPI = -2*a0**2*(r0-r1) - 2*b0**2*(r0-r2) + 2*(b0**2+c0**2)*(r0-r1) &
                & + 2*(a0**2+c0**2)*(r0-r2)

            MS_Force_Point = MS_Force_Point - DFalph/2.0/(S*S0)**0.5*DS
            MS_Force_Point = MS_Force_Point - DFbet*(DPI/(S*S0)**0.5 - (bet+1)/2.0/S*DS)
        enddo

    end subroutine Point_MS_force

end module force_mod
