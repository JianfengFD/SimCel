
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
! Isotropic bending force (numerical finite difference)
! E_v = kpp/2 * (H_v - H0)^2 * S_v/3
! ========================================================
    subroutine Point_H_Force_Cell(C,k,V,Integral_H_A)
        implicit none
        type(cel)::C(:)
        type(PtCell)::V
       integer,intent(in)::k
       real*8 Integral_H_A

       if(C(k)%IS_BC_V(V%V_move).eq.0)then
           call Point_H_Force_numerical(V%V_move,V%PM_Force_Point,C(k)%RV, &
    C(k)%H_modulus, C(k)%H_zero, C(k)%H_V, C(k)%V_V, C(k)%N_V_V, C(k)%area_V, C(k)%IS_BC_V)
       endif
    end subroutine Point_H_Force_Cell


! ========================================================
! Anisotropic curvature-orientation force (numerical finite difference)
! E_aniso_v = [kpp_u/2*(H-H0_u)^2 + kpp_uD/2*(D-D0*cos2th)^2] * phi * A_v/3
! Force = -dE/dr computed by perturbing vertex positions
! ========================================================
    subroutine Point_Aniso_Force_Cell(C,k,V)
        implicit none
        type(cel)::C(:)
        type(PtCell)::V
        integer,intent(in)::k

        V%Aniso_Force_Point = 0.0
        if(C(k)%IS_BC_V(V%V_move).eq.0)then
            call Point_Aniso_Force_numerical(V%V_move, V%Aniso_Force_Point, C(k)%RV, &
                C(k)%kpp_u, C(k)%kpp_uD, C(k)%H0_u, C(k)%D0_u, &
                C(k)%fi, C(k)%q1, C(k)%q2, &
                C(k)%H_V, C(k)%K_V, C(k)%V_V, C(k)%N_V_V, C(k)%F_V, C(k)%V_F, &
                C(k)%area_V, C(k)%area_F, C(k)%Norm_V, C(k)%Norm_F, &
                C(k)%IS_BC_V, C(k)%N_V)
        endif
    end subroutine Point_Aniso_Force_Cell


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
            PV_Force_Point = PV_Force_Point - Dvol_point*P
            do i=1,N_V_V_local
                PV_Force_Point = PV_Force_Point - lamda * Darea_Point(i,1:3)/area_F(F_V(v_move,i))
            enddo
        END subroutine POINT_PV_FORCE


! ========================================================
! Isotropic bending force (numerical FD)
! E_v = kpp/2 * (H_v - H0)^2 * S_v/3
! ========================================================
    subroutine Point_H_Force_numerical(V_move,PM_Force_Point,RV, &
    & H_modulus, H_zero, H_V,V_V,N_V_V, area_V,IS_BC_V)
        implicit none
        real*8 rv(:,:)
        integer V_V(:,:),N_V_V(:)
        integer V_MOVE
        integer IS_BC_V(:)
        integer I_temp(1:20)
        real*8 H_V(:),area_V(:)
        real*8 H_modulus, H_zero
        real*8 area_V_local
        real*8 PM_force_Point(1:3)
        real*8 R(1:20,1:3)
        real*8 R0(1:3)
        real*8 RR(1:20,1:3)
        real*8 RR0(1:3)
        real*8 delt_r
        integer i,j,k,i1, N_t, I_temp2(1:20)
        real*8 H_local,energy_local_1,energy_local_2
        real*8 len_E_local(1:20),th_V_local(1:20)

        i=v_move
        PM_force_Point=0.0
        delt_r=0.000001d0

        R=0
        I_temp=0
        I_temp(1:N_V_V(i))=V_V(i, 1:N_V_V(i))
        do j=1,N_V_V(i)
            if(IS_BC_V(I_temp(j)).eq.0)then
            R(j,1:3)=Rv(I_temp(j),1:3)
            endif
        enddo

        ! ---- x direction
        R0(1:3)=Rv(i,1:3); R0(1)=R0(1)+delt_r
        Energy_local_1=0.0
        N_t=N_V_V(i)
        call get_H_local_NEW(R0,R,N_t, H_local, area_V_local,len_E_local,th_V_local)
        Energy_local_1=Energy_local_1+H_modulus/2.0*(H_local-H_zero)**2*Area_V_local/3.0
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
            Energy_local_1=Energy_local_1+H_modulus/2.0*(H_local-H_zero)**2*Area_V_local/3.0
        endif
        enddo
        R0(1:3)=Rv(i,1:3); R0(1)=R0(1)-delt_r
        Energy_local_2=0.0
        N_t=N_V_V(i)
        call get_H_local_NEW(R0,R,N_t, H_local, area_V_local,len_E_local,th_V_local)
        Energy_local_2=Energy_local_2+H_modulus/2.0*(H_local-H_zero)**2*Area_V_local/3.0
        do j=1,N_V_V(i)
          if(IS_BC_V(I_temp(j)).eq.0)then
            RR0(1:3)=R(j,1:3); i1=I_temp(j)
            I_temp2(1:N_V_V(i1))=V_V(i1, 1:N_V_V(i1))
            do k=1, N_V_V(i1); RR(k,1:3)=Rv(I_temp2(k),1:3); if(I_temp2(k).eq.i)RR(k,1:3)=R0; enddo
            N_t=N_V_V(i1)
            call get_H_local_NEW(RR0,RR,N_t, H_local, area_V_local,len_E_local,th_V_local)
            Energy_local_2=Energy_local_2+H_modulus/2.0*(H_local-H_zero)**2*Area_V_local/3.0
        endif
        enddo
        PM_Force_Point(1)= - (Energy_Local_1-Energy_Local_2)/delt_r/2.0

        ! ---- y direction
        R0(1:3)=Rv(i,1:3); R0(2)=R0(2)+delt_r
        Energy_local_1=0.0
        N_t=N_V_V(i)
        call get_H_local_NEW(R0,R,N_t, H_local, area_V_local,len_E_local,th_V_local)
        Energy_local_1=Energy_local_1+H_modulus/2.0*(H_local-H_zero)**2*Area_V_local/3.0
        do j=1,N_V_V(i)
         if(IS_BC_V(I_temp(j)).eq.0)then
            RR0(1:3)=R(j,1:3); i1=I_temp(j)
            I_temp2(1:N_V_V(i1))=V_V(i1, 1:N_V_V(i1))
            do k=1, N_V_V(i1); RR(k,1:3)=Rv(I_temp2(k),1:3); if(I_temp2(k).eq.i)RR(k,1:3)=R0; enddo
            N_t=N_V_V(i1)
            call get_H_local_NEW(RR0,RR,N_t, H_local, area_V_local,len_E_local,th_V_local)
            Energy_local_1=Energy_local_1+H_modulus/2.0*(H_local-H_zero)**2*Area_V_local/3.0
         endif
        enddo
        R0(1:3)=Rv(i,1:3); R0(2)=R0(2)-delt_r
        Energy_local_2=0.0
        N_t=N_V_V(i)
        call get_H_local_NEW(R0,R,N_t, H_local, area_V_local,len_E_local,th_V_local)
        Energy_local_2=Energy_local_2+H_modulus/2.0*(H_local-H_zero)**2*Area_V_local/3.0
        do j=1,N_V_V(i)
         if(IS_BC_V(I_temp(j)).eq.0)then
            RR0(1:3)=R(j,1:3); i1=I_temp(j)
            I_temp2(1:N_V_V(i1))=V_V(i1, 1:N_V_V(i1))
            do k=1, N_V_V(i1); RR(k,1:3)=Rv(I_temp2(k),1:3); if(I_temp2(k).eq.i)RR(k,1:3)=R0; enddo
            N_t=N_V_V(i1)
            call get_H_local_NEW(RR0,RR,N_t, H_local, area_V_local,len_E_local,th_V_local)
            Energy_local_2=Energy_local_2+H_modulus/2.0*(H_local-H_zero)**2*Area_V_local/3.0
         endif
        enddo
        PM_Force_Point(2)= - (Energy_Local_1-Energy_Local_2)/delt_r/2.0

        ! ---- z direction
        R0(1:3)=Rv(i,1:3); R0(3)=R0(3)+delt_r
        Energy_local_1=0.0
        N_t=N_V_V(i)
        call get_H_local_NEW(R0,R,N_t, H_local, area_V_local,len_E_local,th_V_local)
        Energy_local_1=Energy_local_1+H_modulus/2.0*(H_local-H_zero)**2*Area_V_local/3.0
        do j=1,N_V_V(i)
            if(IS_BC_V(I_temp(j)).eq.0)then
            RR0(1:3)=R(j,1:3); i1=I_temp(j)
            I_temp2(1:N_V_V(i1))=V_V(i1, 1:N_V_V(i1))
            do k=1, N_V_V(i1); RR(k,1:3)=Rv(I_temp2(k),1:3); if(I_temp2(k).eq.i)RR(k,1:3)=R0; enddo
            N_t=N_V_V(i1)
            call get_H_local_NEW(RR0,RR,N_t, H_local, area_V_local,len_E_local,th_V_local)
            Energy_local_1=Energy_local_1+H_modulus/2.0*(H_local-H_zero)**2*Area_V_local/3.0
            endif
        enddo
        R0(1:3)=Rv(i,1:3); R0(3)=R0(3)-delt_r
        Energy_local_2=0.0
        N_t=N_V_V(i)
        call get_H_local_NEW(R0,R,N_t, H_local, area_V_local,len_E_local,th_V_local)
        Energy_local_2=Energy_local_2+H_modulus/2.0*(H_local-H_zero)**2*Area_V_local/3.0
        do j=1,N_V_V(i)
         if(IS_BC_V(I_temp(j)).eq.0)then
            RR0(1:3)=R(j,1:3); i1=I_temp(j)
            I_temp2(1:N_V_V(i1))=V_V(i1, 1:N_V_V(i1))
            do k=1, N_V_V(i1); RR(k,1:3)=Rv(I_temp2(k),1:3); if(I_temp2(k).eq.i)RR(k,1:3)=R0; enddo
            N_t=N_V_V(i1)
            call get_H_local_NEW(RR0,RR,N_t, H_local, area_V_local,len_E_local,th_V_local)
            Energy_local_2=Energy_local_2+H_modulus/2.0*(H_local-H_zero)**2*Area_V_local/3.0
        endif
        enddo
        PM_Force_Point(3)= - (Energy_Local_1-Energy_Local_2)/delt_r/2.0

    END subroutine Point_H_Force_numerical


! ========================================================
! Anisotropic coupling force via numerical FD
! Perturb vertex, recompute local H, K, D, curvature tensor, and coupling energy
! ========================================================
    subroutine Point_Aniso_Force_numerical(V_move, Aniso_Force_Point, RV, &
        kpp_u, kpp_uD, H0_u, D0_u, &
        fi, q1, q2, &
        H_V, K_V, V_V, N_V_V, F_V, V_F, &
        area_V, area_F, Norm_V, Norm_F, IS_BC_V, N_V_total)
        implicit none
        integer, intent(in):: V_move, N_V_total
        real*8, intent(out):: Aniso_Force_Point(1:3)
        real*8, intent(in):: rv(:,:), fi(:), q1(:), q2(:)
        real*8, intent(in):: H_V(:), K_V(:), area_V(:), area_F(:)
        real*8, intent(in):: Norm_V(:,:), Norm_F(:,:)
        real*8, intent(in):: kpp_u, kpp_uD, H0_u, D0_u
        integer, intent(in):: V_V(:,:), N_V_V(:), F_V(:,:), V_F(:,:), IS_BC_V(:)

        real*8 delt_r, E1, E2
        real*8 R0_save(1:3), R0_pert(1:3)
        real*8 R(1:20,1:3)
        real*8 H_local, area_V_local
        real*8 len_E_local(1:20), th_V_local(1:20)
        real*8 D_local, K_local, S_local, P_local, cos2th
        real*8 E_aniso_v
        ! for local Gaussian curvature
        real*8 angle_sum_local, e1_loc(1:3), e2_loc(1:3), cos_ang, ang, l1, l2
        ! for local curvature tensor
        real*8 nrm_local(1:3), n_len
        real*8 e_ij(1:3), e_len, kn_loc, t_proj(1:3), t_len
        real*8 T1_loc(1:3), T2_loc(1:3), cphi, sphi
        real*8 M11, M12, M22, w_s, d1_loc, d2_loc, Dmag_loc
        real*8 dd1_loc, dd2_loc

        integer i, j, j1, k1, i_f, N_t, dir
        integer I_temp(1:20)

        i = V_move
        Aniso_Force_Point = 0.0
        delt_r = 0.000001d0
        N_t = N_V_V(i)
        I_temp(1:N_t) = V_V(i, 1:N_t)

        ! neighbor positions
        R = 0
        do j = 1, N_t
            R(j,1:3) = rv(I_temp(j), 1:3)
        enddo

        R0_save = rv(i, 1:3)

        do dir = 1, 3
            ! +delta
            R0_pert = R0_save
            R0_pert(dir) = R0_pert(dir) + delt_r
            call compute_local_aniso_energy(R0_pert, R, N_t, I_temp, &
                fi, q1, q2, kpp_u, kpp_uD, H0_u, D0_u, &
                rv, V_V, N_V_V, F_V, V_F, area_F, Norm_F, IS_BC_V, i, N_V_total, E1)
            ! -delta
            R0_pert = R0_save
            R0_pert(dir) = R0_pert(dir) - delt_r
            call compute_local_aniso_energy(R0_pert, R, N_t, I_temp, &
                fi, q1, q2, kpp_u, kpp_uD, H0_u, D0_u, &
                rv, V_V, N_V_V, F_V, V_F, area_F, Norm_F, IS_BC_V, i, N_V_total, E2)

            Aniso_Force_Point(dir) = -(E1 - E2) / (2.0d0 * delt_r)
        enddo

    end subroutine Point_Aniso_Force_numerical


! ========================================================
! Helper: compute local anisotropic energy at vertex i with perturbed position
! ========================================================
    subroutine compute_local_aniso_energy(R0, R_nb, N_t, I_nb, &
        fi, q1, q2, kpp_u, kpp_uD, H0_u, D0_u, &
        rv, V_V, N_V_V, F_V, V_F, area_F, Norm_F, IS_BC_V, i_center, N_V_total, E_out)
        implicit none
        real*8, intent(in):: R0(1:3), R_nb(1:20,1:3)
        integer, intent(in):: N_t, I_nb(1:20), i_center, N_V_total
        real*8, intent(in):: fi(:), q1(:), q2(:)
        real*8, intent(in):: kpp_u, kpp_uD, H0_u, D0_u
        real*8, intent(in):: rv(:,:), area_F(:), Norm_F(:,:)
        integer, intent(in):: V_V(:,:), N_V_V(:), F_V(:,:), V_F(:,:), IS_BC_V(:)
        real*8, intent(out):: E_out

        real*8 H_local, area_V_local
        real*8 len_E_local(1:20), th_V_local(1:20)
        real*8 R_work(1:20,1:3)
        ! Gaussian curvature local
        real*8 angle_sum, e1_v(1:3), e2_v(1:3), cos_a, ang_v, l1, l2
        real*8 K_local, D_local, disc
        ! curvature tensor local
        real*8 nrm(1:3), n_len, e_ij(1:3), e_len, kn
        real*8 t_p(1:3), t_l, T1(1:3), T2(1:3), cp, sp
        real*8 M11, M12, M22, ws, d1, d2, Dm, dd1, dd2
        real*8 S_loc, P_loc, cos2th, E_v
        integer j

        R_work(1:N_t, 1:3) = R_nb(1:N_t, 1:3)

        ! compute local H via dihedral angles
        call get_H_local_NEW(R0, R_work, N_t, H_local, area_V_local, len_E_local, th_V_local)

        ! compute local K via angle deficit
        R_work(N_t+1, 1:3) = R_work(1, 1:3)
        angle_sum = 0.0d0
        do j = 1, N_t
            e1_v = R_work(j, 1:3) - R0
            e2_v = R_work(mod(j, N_t)+1, 1:3) - R0
            l1 = sqrt(sum(e1_v**2))
            l2 = sqrt(sum(e2_v**2))
            if(l1.gt.1d-15 .and. l2.gt.1d-15)then
                cos_a = sum(e1_v*e2_v) / (l1*l2)
                cos_a = max(-1.0d0, min(1.0d0, cos_a))
                ang_v = dacos(cos_a)
                angle_sum = angle_sum + ang_v
            endif
        enddo
        if(area_V_local.gt.1d-15)then
            K_local = (2.0d0*Pai - angle_sum) / (area_V_local/3.0d0)
        else
            K_local = 0.0d0
        endif

        disc = H_local**2 - K_local
        if(disc.gt.0.0d0)then
            D_local = sqrt(disc)
        else
            D_local = 0.0d0
        endif

        ! compute local curvature tensor for direction
        ! vertex normal from face normals (area-weighted)
        nrm = 0.0d0
        do j = 1, N_t
            e1_v = R_work(j,1:3) - R0
            e2_v = R_work(mod(j,N_t)+1,1:3) - R0
            nrm(1) = nrm(1) + R1R2_fun(e1_v, e2_v, 1)
            nrm(2) = nrm(2) + R1R2_fun(e1_v, e2_v, 2)
            nrm(3) = nrm(3) + R1R2_fun(e1_v, e2_v, 3)
        enddo
        n_len = sqrt(sum(nrm**2))
        if(n_len.gt.1d-15)then
            nrm = nrm / n_len
        else
            nrm = (/0.0d0, 0.0d0, 1.0d0/)
        endif

        ! tangent frame
        e_ij = R_work(1,1:3) - R0
        t_p = e_ij - sum(e_ij*nrm)*nrm
        t_l = sqrt(sum(t_p**2))
        if(t_l.lt.1d-15)then
            E_out = 0.0d0; return
        endif
        T1 = t_p / t_l
        T2(1) = nrm(2)*T1(3) - nrm(3)*T1(2)
        T2(2) = nrm(3)*T1(1) - nrm(1)*T1(3)
        T2(3) = nrm(1)*T1(2) - nrm(2)*T1(1)

        ! curvature tensor
        M11=0; M12=0; M22=0; ws=0
        do j = 1, N_t
            e_ij = R_work(j,1:3) - R0
            e_len = sqrt(sum(e_ij**2))
            if(e_len.lt.1d-15) cycle
            kn = 2.0d0 * sum(nrm*e_ij) / (e_len**2)
            t_p = e_ij - sum(e_ij*nrm)*nrm
            t_l = sqrt(sum(t_p**2))
            if(t_l.lt.1d-15) cycle
            cp = sum(t_p*T1)/t_l
            sp = sum(t_p*T2)/t_l
            M11 = M11 + kn*cp*cp
            M12 = M12 + kn*cp*sp
            M22 = M22 + kn*sp*sp
            ws = ws + 1.0d0
        enddo
        if(ws.gt.0.5d0)then
            M11=M11/ws; M12=M12/ws; M22=M22/ws
        endif
        d1 = (M11-M22)/2.0d0
        d2 = M12
        Dm = sqrt(d1**2 + d2**2)

        if(Dm.gt.1d-15)then
            dd1 = D_local * d1/Dm
            dd2 = D_local * d2/Dm
        else
            dd1 = 0.0d0; dd2 = 0.0d0
        endif

        ! anisotropic energy at this vertex
        S_loc = sqrt(q1(i_center)**2 + q2(i_center)**2)
        if(S_loc.gt.1d-10 .and. D_local.gt.1d-10)then
            P_loc = q1(i_center)*dd1 + q2(i_center)*dd2
            cos2th = P_loc / (S_loc * D_local)
            cos2th = max(-1.0d0, min(1.0d0, cos2th))
        else
            cos2th = 0.0d0
        endif

        E_v = (kpp_u/2.0d0*(H_local-H0_u)**2 + kpp_uD/2.0d0*(D_local-D0_u*cos2th)**2) &
            * (1.0d0+fi(i_center))*0.5d0 * area_V_local/3.0d0

        E_out = E_v

    end subroutine compute_local_aniso_energy


    Subroutine Update_New_Parameters(v_move,rv_move,rv,N_V_V,V_V,E_V,F_V,F_E,E_V_opposite,V_E_opposite,&
    & Len_E,th_E,vector_F,Norm_F,Norm_V,area_F,area_V,area_tot,vol_F,vol_total,H_V)
        implicit none
        integer N_V_V(:)
        integer V_V(:,:),E_V(:,:),F_V(:,:),F_E(:,:)
        integer E_V_opposite(:,:),V_E_opposite(:,:)
        integer i,j,k,i1,i2,i3,I_temp(1:40),I_temp2(1:40)
        integer N_t
        real*8 rv(:,:),R0(1:3)
        real*8 R(1:40,1:3),RR3(1:3,1:3),R1(1:3),R2(1:3),R3(1:3)
        real*8 RJ(1:40,1:3)
        real*8 len_E(:),th_E(:),th_JK,len_JK(1:20),sign_r
        real*8 vector_F(:,:),Vector_R0(1:40,1:3),Vector_RJ(1:40,1:3),vector_Local(1:3)
        real*8 Norm_F(:,:),Norm_V(:,:)
        real*8 area_F(:),area_V(:),area_tot,area_local,area_V_local
        real*8 sum_area1_point,sum_area2_point
        real*8 Vol_F(:),Vol_total, sum_vol1_point,sum_vol2_point
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
            RR3(1,1:3)=r1(1:3); RR3(2,1:3)=r2(1:3); RR3(3,1:3)=r3(1:3)
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
            i3=0; k_a=0; k_b=0
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
            a0=Len_E_zero(I_side(k_a)); b0=Len_E_zero(I_side(k_b)); c0=Len_E_zero(I_side(k_c))
            a =Len_E(I_side(k_a)); b =Len_E(I_side(k_b)); c =Len_E(I_side(k_c))
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
