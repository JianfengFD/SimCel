module basicfun_mod
 use CelType_mod
 use allocCel_mod
 implicit none
    contains

    subroutine Get_shape_information(C,k)
            type(cel)::C(:)
           integer,intent(in)::k
           integer i,j,i1

           call Get_Cell_area(C,K)
           call Get_Cell_Len(C,k)
           call Get_Cell_Vol(C,k)
           call Get_Cell_H(C,K)
           call Get_Cell_Energy(C,K)

           do i=1,C(k)%N_V
               C(k)%Norm_V(i,:)=0
               C(k)%area_V(i) = 0
                do j=1,C(k)%N_V_V(i)
                    i1 = C(k)%F_V(i,j)
                    C(k)%Norm_V(i,:) = C(k)%Norm_V(i,:)+C(k)%Norm_F(i1,:)*C(k)%area_F(i1)
                    C(k)%area_V(i) = C(k)%area_V(i) + C(k)%area_F(i1)
                enddo
                C(k)%Norm_V(i,:) = C(k)%Norm_V(i,:)/sum(C(k)%Norm_V(i,:)**2)**0.5
           enddo


    end subroutine Get_shape_information

    subroutine Get_shape_sgl_information(C)
            type(cel)::C
            integer i,j,i1

           call Get_area_F(C%rv,  C%V_F, C%N_F, C%Area_F, C%Vector_F,  C%Norm_F)
           C%area_tot = sum(C%Area_F(1:C%N_F))
            do i=1,C%N_V
               C%area_V(i) = 0
                do j=1,C%N_V_V(i)
                    i1 = C%F_V(i,j)
                    C%area_V(i) = C%area_V(i) + C%area_F(i1)
                enddo

           enddo


           call Get_Len_E(C%rv,  C%V_E, C%N_E, C%Len_E, C%Vector_E)
           call Get_Vol_F(C%rv,  C%V_F, C%N_F, C%Vol_F,&
           C%Vol_total)
            call Get_H_V_new(C%rv,C%N_V,C%N_F_V,C%V_V, C%F_V, C%E_V,&
  C%area_F, C%Vector_F,  C%H_V,C%area_V,C%Darea_V, C%Dvol_V,C%th_E,C%IS_BC_V)

        C%EnergyH=Curvature_Energy( C%H_modulus, C%H_v, C% H_zero,  C%Area_V, C%N_V,C%IS_BC_V,C%IS_NEXTBC_V)
           C%energyL=C%kpp_area/2.0*(C%area_tot- &
           C%area_tot_zero)**2.0/C%area_tot_zero
            C%energyP=C%Kpp_V/2.0*(C%Vol_total- &
            C%Vol_total_zero)**2/C%Vol_total_zero
            C%energy_tot=C%energyp+C%energyH+C%energyL

    end subroutine Get_shape_Sgl_information





     subroutine Get_Cell_maxmin_Len(C)
        type(cel)::C
        integer i,j

        c%maxL = 0.0
        c%minL = 1000000

        do i = 1, c%N_E
            if(c%Len_E(i).gt.c%maxL)c%maxL=c%Len_E(i)
            if(c%Len_E(i).lt.c%minL)c%minL=c%Len_E(i)
        enddo


    end subroutine Get_Cell_maxmin_Len


    subroutine Get_Cell_Len(C,k)
            type(cel)::C(:)
           integer,intent(in)::k
           call Get_Len_E(C(k)%rv,  C(k)%V_E, C(k)%N_E, C(k)%Len_E, C(k)%Vector_E)
    end subroutine Get_Cell_Len

        subroutine Get_Cell_area(C,k)
            type(cel)::C(:)
           integer,intent(in)::k


           call Get_area_F(C(k)%rv,  C(k)%V_F, C(k)%N_F, C(k)%Area_F, C(k)%Vector_F,  C(k)%Norm_F)
           C(k)%area_tot = sum(C(k)%Area_F(1:C(k)%N_F))

    end subroutine Get_Cell_area

    subroutine Get_Cell_Vol(C,k)
        implicit none
            type(cel)::C(:)
           integer,intent(in)::k
           call Get_Vol_F(C(k)%rv,  C(k)%V_F, C(k)%N_F, C(k)%Vol_F,&
           C(k)%Vol_total)
    end subroutine Get_Cell_Vol

        subroutine Get_Cell_H(C,k)
            type(cel)::C(:)
           integer,intent(in)::k
           call Get_H_V_new(C(k)%rv,C(k)%N_V,C(k)%N_F_V,C(k)%V_V, C(k)%F_V, C(k)%E_V,&
  C(k)%area_F, C(k)%Vector_F,  C(k)%H_V,C(k)%area_V,C(k)%Darea_V, C(k)%Dvol_V,C(k)%th_E,C(k)%IS_BC_V)
        end subroutine Get_Cell_H



! this subr get the relative position of a point with respect to a triangle.
    subroutine Pos_Point_Triangle(R,R1,R2,R3,n,a,b,d,is_ON,is_IN)
        implicit none
        real*8 R(1:3),R1(1:3),R2(1:3),R3(1:3)
        real*8 n(1:3),a,b,d
        real*8 aa(1:3,1:3),aout(1:3,1:3),bb(1:3)
        real*8 det
        integer k,i,j,Is_ON,IS_IN
        bb= R -r1
        aa(:,1) = n  ! [n;r2-r1;r3-r1][d,a,b]
        aa(:,2) = r2-r1
        aa(:,3) = r3-r1

        det = V6_fun(aa)
        if(abs(det).lt.0.00001)then
           ! write(*,*)'error in pos_point'
        else
          aout = aa; aout(:,1) = bb
        d = v6_fun(aout)/det
        aout = aa; aout(:,2) = bb
        a = v6_fun(aout)/det
        aout = aa; aout(:,3) = bb
        b = v6_fun(aout)/det

        endif

        is_on =0;is_in=0
        if(a.ge.0.and.a.le.1.and.b.ge.0.and.b.le.1-a)then
                is_on = 1
                if(abs(d).lt.0.0001)is_in = 1
        endif


    end subroutine Pos_Point_Triangle






    Subroutine Get_Len_E(rv,  V_E, N_E, Len_E, Vector_E)
        implicit none
        real*8 rv(:,:)
        integer V_E(:, :), N_E
        real*8 Len_E(:)  !length of edges
        real*8 Vector_E(:,:)  !Vector of edges

        integer i,i1,i2

        real*8 a(1:3),d

        do i=1,N_E
            i1=V_E(i,1)
            i2=V_E(i,2)
            a=rv(i2,1:3)-rv(i1,1:3)
            d=(a(1)*a(1)+a(2)*a(2)+a(3)*a(3))**(0.5)
            Vector_E(i,1:3)=a(1:3)
            Len_E(i)=d
        enddo


    END Subroutine Get_Len_E

! calculate area vector & area of local

    subroutine Get_area_local(R1,R2,R3,area_local, vector_local)
        real*8 R1(1:3),R2(1:3),R3(1:3), area_local,vector_local(1:3)
        real*8 R(1:3,1:3)
            R(1,1:3)=R1(1:3);R(2,1:3)=R2(1:3);R(3,1:3)=R3(1:3)
            Vector_local(1)=RR_fun(R,1)/2.0
            Vector_local(2)=RR_fun(R,2)/2.0
            Vector_local(3)=RR_fun(R,3)/2.0
            Area_local= (Vector_local(1)**2+Vector_local(2)**2+Vector_local(3)**2)**0.5
    END subroutine Get_area_local



    ! Calculate area vectors and area of small triangles
    subroutine Get_area_F( rv, V_F,  N_F,   Area_F, Vector_F, Norm_F)
        implicit none
        integer N
        parameter(N = 50000)

        integer N_F
        real*8 rV(:,:)
        integer V_F(:, :)

        real*8  Vector_F(:,:), Area_F(:), Norm_F(:,:)
        real*8 Rt(1:3,1:3)
        ! function defintion


        integer i,j,i_VF(1:3)


        do i=1,N_F
            i_VF=V_F(i,1:3)
            Rt(1,1:3)=rv(i_VF(1),1:3)
            Rt(2,1:3)=rv(i_VF(2),1:3)
            Rt(3,1:3)=rv(i_VF(3),1:3)
             do J=1,3
                 Vector_F(i,J)=RR_fun(Rt,J)/2.0
             enddo
             Area_F(i)= (Vector_F(i,1)**2 +Vector_F(i,2)**2+Vector_F(i,3)**2 )**(0.5)
             Norm_F(i,1:3)=vector_F(i,1:3)/Area_F(i)
        enddo
        ! Norm_F=0

    END subroutine Get_area_F



! calculate volumn
    subroutine Get_vol_F( rv, V_F,  N_F, Vol_F, Vol_total)
        implicit none
        integer N
        parameter(N = 50000)

        integer N_F
        real*8 rV(:,:)
        integer V_F(:, :)

        real*8 Vol_F(:), vol_total
        real*8 R(1:3,1:3)
        integer i,i_VF(1:3)

        do i=1,N_F
            i_VF=V_F(i,1:3)
            R(1,1:3)=rv(i_VF(1),1:3)
            R(2,1:3)=rv(i_VF(2),1:3)
            R(3,1:3)=rv(i_VF(3),1:3)
            Vol_F(i)=v6_fun(R)/6.0
        enddo
        Vol_total=sum( Vol_F(1:N_F) )

    END  subroutine Get_vol_F

! curvature of local

    subroutine Get_H_local_new(R0,R,N_V_V_local,   H_local, area_V_local,len_E_local,th_V_local)

        real*8 R(1:20,1:3)
        real*8 R0(1:3),R1(1:3),R2(1:3),R3(1:3)
        real*8 H_local, area_R0(1:20), vector_R0(1:20,1:3)
        real*8 area_local,vector_local(1:3)
        real*8 Len_jk(1:20),th_jk
        real*8 H_jk,sign_r
        real*8 Len_E_local(1:20),th_V_local(1:20)

        real*8 area_V_local

        integer i,i1 , N_V_V_local

        R(N_V_V_local+1,1:3)=R(1,1:3)
        R(N_V_V_local+2,1:3)=R(2,1:3)
        area_V_local=0.0
        do i=1,N_V_V_local
            R1(1:3)=R0(1:3);R2(1:3)=R(i,1:3);R3(1:3)=R(i+1,1:3);
            call Get_area_local(R1,R2,R3,area_local, vector_local)
            area_R0(i)=area_local

            vector_R0(i,1:3)= vector_local(1:3)/area_Local
            area_V_local=area_V_local+area_local
            len_jk(i)=(sum((R1-R3)**2.0))**0.5
        enddo
        vector_r0(N_V_V_local+1,1:3)=vector_R0(1,1:3)
        len_Jk(N_V_V_local+1)=len_JK(1)
        H_JK=0.0
        do i=1,N_V_V_local
            sign_r=sum(vector_R0(i,1:3)*vector_R0(i+1,1:3))
            if(sign_r.lt.1.0)then
                th_jk=dacos(sign_r)
            else
            th_JK=0.0
            endif
            sign_r=sum( (vector_R0(i,1:3)-vector_R0(i+1,1:3))*(R(i,1:3)-R(i+2,1:3)) )
            if(sign_r.lt.0.0)th_jk=-th_JK
            i1=mod(i,N_v_V_local)+1
            th_V_local(i1)=th_JK
            len_E_local(i)=len_JK(i)


            H_jk=H_jk+th_jk*len_jk(i)
        enddo
        H_local=-H_JK/(area_V_local/3.0)/2.0

    END  subroutine Get_H_local_new


! calculate curvature & some basic variations S and V

    subroutine Get_H_V_new(rv,N_V,N_F_V, V_V, F_V, E_V,&
&  area_F, Vector_F,  H_V,area_V,Darea_V, Dvol_V,th_E,IS_BC_V)
    implicit none

    integer N_V

    real*8 rV(:,:)

    ! %%%%%%%%%%%%% other labels
    ! the first kind label
    integer E_V(:,:) ! edges ajecent to a vertice
    integer F_V(:,:), N_F_V(:) ! star faces around a vertice


    ! second kind label
    integer V_V(:,:)
    integer IS_BC_V(:)

    ! other input & output

    real*8 H_V(:)
    real*8  Vector_F(:,:), Area_F(:), area_V(:)
    real*8 R_VV(1:20,1:3), Vect_VV(1:20,1:3)
    real*8 ForceV_v(1:3), Darea_V(:,:, :)  ! \partial area/\partial rv
    real*8 DVol_V(:,:)  ! \partial volume/\partial rv
    real*8 Norm_V_temp(1:3)
    real*8 at(1:3)
    real*8 len_jk(1:20),H_jk,th_jk,sign_r
    real*8 th_E(:)

    integer i,j,i1,I_temp(1:20),I_temp2(1:20)



    do i=1, N_V
    if(IS_BC_V(i).eq.0)then ! Ensure V is not on BC
        I_temp(1:N_F_V(i))=V_V(i,1:N_F_V(i))
        I_temp2(1:N_F_V(i))=F_V(i,1:N_F_V(i))
        do j=1,N_F_V(i)
            R_VV(j, 1:3)=Rv(I_temp(j) , 1:3)
            Vect_VV(j, 1:3)=Vector_F(I_temp2(j) , 1:3)/area_F(I_temp2(j))

        enddo
        vect_VV(N_F_V(i)+1,1:3)=Vect_VV(1,1:3)
        R_VV(N_F_V(i)+1,1:3)=R_VV(1,1:3)
        R_VV(N_F_V(i)+2,1:3)=R_VV(2,1:3)
        Norm_V_temp=0.0
        ForceV_V=0
        area_V(i)=0



        do j=1,N_F_V(i)
            area_V(i)=area_V(i)+Area_F( I_temp2(j) )
            Norm_V_temp(1)=Norm_V_temp(1)+R1R2_fun(R_VV(j,1:3), R_VV(j+1,1:3), 1)
            Norm_V_temp(2)=Norm_V_temp(2)+R1R2_fun(R_VV(j,1:3), R_VV(j+1,1:3), 2)
            Norm_V_temp(3)=Norm_V_temp(3)+R1R2_fun(R_VV(j,1:3), R_VV(j+1,1:3), 3)
            at(1)=( R1R2_fun(R_VV(j,1:3), Vect_VV(j,1:3), 1) &
            &+R1R2_fun( Vect_VV(j,1:3),R_VV(j+1,1:3), 1)  )/ 2.0
            at(2)=( R1R2_fun(R_VV(j,1:3), Vect_VV(j,1:3), 2) &
            &+R1R2_fun( Vect_VV(j,1:3),R_VV(j+1,1:3), 2)  )/ 2.0
            at(3)=( R1R2_fun(R_VV(j,1:3), Vect_VV(j,1:3), 3) &
            &+R1R2_fun( Vect_VV(j,1:3),R_VV(j+1,1:3), 3)  )/ 2.0
            Darea_V(i,j,1:3)=at(1:3)
            len_JK(j)=(sum((rv(i,1:3)-r_VV(j,1:3))**2.0))**0.5
        enddo

        len_JK(N_F_V(i)+1)=len_JK(1)
        h_JK=0.0
        do j=1,N_F_V(i)
            sign_r=sum(vect_VV(j,1:3)*vect_VV(j+1,1:3))
            if(sign_r.lt.1.0)then
                th_jk=dacos(sign_r)
            else
                th_JK=0.0
            endif
            i1=mod(j,N_F_V(i))+1
            sign_r=sum( (vect_VV(j,1:3)-vect_VV(j+1,1:3))*(R_VV(j,1:3)-R_VV(j+2,1:3)) )
            if(sign_r.lt.0.0)th_jk=-th_JK
            h_jk=h_jk+th_jk*len_JK(j+1)
            th_E(abs(E_V(i,i1)))=th_jk

        enddo

        H_V(i)=-h_jk/(area_V(i)/3.0)/2.0
        Dvol_V(i,1:3)=Norm_v_temp/6.0
    endif
    enddo
END subroutine Get_H_V_new

! ========================================================
! Get_Cell_Energy: includes bending(phi-coupled), area, volume,
! shear, and phase field energies
! ========================================================
    subroutine Get_Cell_Energy(C,k)
        implicit none
            type(cel)::C(:)
           integer,intent(in)::k
           ! Bending energy with phi coupling
           C(k)%EnergyH = Curvature_Energy_phi_fn(C, k)
           ! Area constraint energy
           C(k)%energyL=C(k)%kpp_area/2.0*(C(k)%area_tot- &
           C(k)%area_tot_zero)**2.0/C(k)%area_tot_zero
           ! Volume constraint
            C(k)%energyP=C(k)%Kpp_V/2.0*(C(k)%Vol_total- &
            C(k)%Vol_total_zero)**2/C(k)%Vol_total_zero
            ! Shear (MS) energy
            C(k)%energyMs = 0.0
            if(C(k)%HAS_REFERENCE.eq.1)then
                call get_energy_MS_cell(C,k)
            endif
            ! Phase field energy (GL + interface gradient)
            call phase_energy_Cell(C, k)
            C(k)%energy_tot=C(k)%energyp+C(k)%energyH+C(k)%energyL+C(k)%energyMs+C(k)%energyPh

    end subroutine Get_Cell_energy








    real*8 function Curvature_Energy(H_modulus,H_v, H_zero, Area_V,N_V,IS_BC_V,IS_NEXTBC_V)
        implicit none

        integer N_V
        real*8 H_modulus, H_zero
        real*8   Area_V(:),H_v(:)
        real*8 a
        integer i,IS_BC_V(:),IS_NEXTBC_V(:)

        a=0.0
        do i=1,N_V
            if(IS_BC_V(i).eq.0.and.IS_NEXTBC_V(i).eq.0)then
             a=a+ H_modulus / 2.0 * (H_v(i)-H_ZERO )**2.0* Area_V(i)/3.0
            endif
        enddo
        Curvature_Energy=a
    END function Curvature_Energy
!  ************************************************
!  Part of 3D vesicle program
! Some functions
!                    written by LiJF  Dobb@bbs.fudan.edu.cn
! *************************************************

! input R1, R2 output: |R1\\R2|^alpha
    Real*8 Function R1R2_fun(R1,R2,I_c)
        real*8 R1(1:3), R2(1:3)
        integer I_C

        select case(I_c)
        case(1)! x (y z)

        R1R2_fun=R1(2)*R2(3)-R1(3)*R2(2)

        case(2) ! y (z x)

        R1R2_fun=R1(3)*R2(1)-R1(1)*R2(3)

        case(3) ! z (x y)

        R1R2_fun=R1(1)*R2(2)-R1(2)*R2(1)

        case default
        write(*,*) 'wrong input within RR_fun'
        !pause
        end select




    END  Function R1R2_fun

    ! input R1, R2, R3 output: ||R1\\R2| \\ R3 |^alpha
    Real*8 Function R1R2_R3_fun(R1,R2,R3,I_c)
        real*8 R1(1:3), R2(1:3), R3(1:3)
        real*8 D1(1:3)
        integer I_C,k

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
        !pause
        end select

    END Function R1R2_R3_fun

    ! input R1, R2, R3 output: ||R1\\R2| \\ R3 |^alpha
    Real*8 Function R1_R2R3_fun(R1,R2,R3,I_c)
        real*8 R1(1:3), R2(1:3), R3(1:3)
        real*8 D1(1:3)
        integer I_C,k

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
        !pause
        end select

    END function R1_R2R3_fun


    ! input R(1:3,1:3), I_c(x,y, or z) output : | |ex+| |ey +| |ez

    Real*8 function  RR_fun(R, I_c)

               ! (P_Num ,Cordinate )
        Real*8    R(1:3,  1:3)
        integer I_c

        select case(I_c)
            case(1)! x (y z)
                    !  1x                               2x(3 1)                          3x(1 2)
                RR_fun= R(2,2)*R(3,3)-R(2,3)*R(3,2)  + R(3,2)*R(1,3)-R(3,3)*R(1,2) &
                &+   R(1,2)*R(2,3)-R(1,3)*R(2,2)


            case(2) ! y (z x)
                    !  1y                               2y(3 1)                          3y(1 2)
                RR_fun= R(2,3)*R(3,1)-R(2,1)*R(3,3)  + R(3,3)*R(1,1)-R(3,1)*R(1,3) &
                &+   R(1,3)*R(2,1)-R(1,1)*R(2,3)

            case(3) ! z (x y)
                    !  1z                               2z(3 1)                          3z(1 2)
                RR_fun= R(2,1)*R(3,2)-R(2,2)*R(3,1)  + R(3,1)*R(1,2)-R(3,2)*R(1,1) &
                & +   R(1,1)*R(2,2)-R(1,2)*R(2,1)
            case default
                write(*,*) 'wrong input within RR_fun'
                !pause
        end select


    END function RR_fun

    !  input R(1:3,1:3) Output |R\\R\\R|  Volume function

    real *8 Function V6_fun(R)

        Real*8    R(1:3,  1:3)

        V6_fun=R(1,1)*( R(2,2)*R(3,3)-R(2,3)*R(3,2) )+R(1,2)*(R(2,3)*R(3,1)&
        &- R(2,1)*R(3,3))+R(1,3)*(R(2,1)*R(3,2)-R(2,2)*R(3,1) )


    END Function V6_fun



      	    REAL *8 FUNCTION RANDOM2(seed11)
            real *8 seed11,rM,rJ
            rJ=13807.0
	        rM=2.0**31.0-1.0
            seed11=MOD(seed11*rJ,rM)
            RANDOM2=seed11/rM
            END  FUNCTION RANDOM2




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
          END FUNCTION random1


! ========================================================
! Shear (MS) Energy
! ========================================================

    subroutine Get_alph_bet_analytic(a,b,c,a0,b0,c0,alph,bet,S,S0)
        implicit none
        real*8,intent(in):: a,b,c,a0,b0,c0
        real*8,intent(out):: alph,bet,S,S0
        real*8 SS0

        S  =(2*a**2*b**2+2*a**2*c**2+2*b**2*c**2-a**4-b**4-c**4)
        S0 =(2*a0**2*b0**2+2*a0**2*c0**2+2*b0**2*c0**2-a0**4-b0**4-c0**4)
        SS0=a**2*b0**2+a**2*c0**2+b**2*c0**2+a0**2*b**2+a0**2*c**2+b0**2*c**2- &
        & a**2*a0**2-b0**2*b**2-c**2*c0**2

        alph=(S/S0)**0.5-1.0
        bet=SS0/(S*S0)**0.5-1.0

    end subroutine Get_alph_bet_analytic

    subroutine get_energy_MS_cell(Cls,k)
        implicit none
        type(cel)::Cls(:)
        integer,intent(in)::k
        integer i,i1,i2,i3
        real*8 ea,eb,ec,ea0,eb0,ec0,alph,bet,S,S0

        Cls(k)%energyMs=0.0
        do i=1,Cls(k)%N_F
            i1=abs(Cls(k)%E_F(i,1)); i2=abs(Cls(k)%E_F(i,2)); i3=abs(Cls(k)%E_F(i,3))
            ea=Cls(k)%Len_E(i1); eb=Cls(k)%Len_E(i2); ec=Cls(k)%Len_E(i3)
            ea0=Cls(k)%Len_E_zero(i1); eb0=Cls(k)%Len_E_zero(i2); ec0=Cls(k)%Len_E_zero(i3)

            call Get_alph_bet_analytic(ea,eb,ec,ea0,eb0,ec0,alph,bet,S,S0)
            Cls(k)%energyMs = Cls(k)%energyMs + Cls(k)%kpp_alpha/2.0 * &
                & (alph**2.0 + Cls(k)%a3_ms*alph**3 + Cls(k)%a4_ms*alph**4) * Cls(k)%Area_F_zero(i)
            Cls(k)%energyMs = Cls(k)%energyMs + Cls(k)%mu_ms * &
                & (bet + Cls(k)%b1_ms*alph*bet + Cls(k)%b2_ms*bet**2) * Cls(k)%Area_F_zero(i)
        enddo

    end subroutine get_energy_MS_cell


! ========================================================
! Phase field: Curvature energy with phi coupling
! E = sum_V (kpp + kpp1*fi(v))/2 * (H_v + c0 - c1*fi(v))^2 * area_V/3
! ========================================================

    real*8 function Curvature_Energy_phi_fn(C, k)
        implicit none
        type(cel)::C(:)
        integer, intent(in)::k
        real*8 a, kpp_local, H0_local
        integer i

        a = 0.0
        do i = 1, C(k)%N_V
            if(C(k)%IS_BC_V(i).eq.0 .and. C(k)%IS_NEXTBC_V(i).eq.0)then
                kpp_local = C(k)%H_modulus + C(k)%kpp1_phi * C(k)%fi(i)
                H0_local = -C(k)%c0_phi + C(k)%c1_phi * C(k)%fi(i)
                a = a + kpp_local / 2.0 * (C(k)%H_V(i) - H0_local)**2 * C(k)%area_V(i) / 3.0
            endif
        enddo
        Curvature_Energy_phi_fn = a
    end function Curvature_Energy_phi_fn


! ========================================================
! Phase field: discrete Laplacian (cotangent weights)
! Ported from VesPhase/Basic_parameter.f90
! ========================================================

    subroutine Laplace_Local_sub(R0, R, area_FV, area_VV, L_laplace_local, N_t)
        implicit none
        real*8 R0(1:3), R(1:20,1:3)
        real*8 area_FV(1:20), area_VV
        real*8 L_laplace_local(1:20)
        real*8 a(1:3), b(1:3), cc(1:3)
        integer N_T, i

        R(N_T+1,1:3) = R(1,1:3)
        area_FV(N_T+1) = area_FV(1)
        L_laplace_local = 0.0
        do i = 1, N_T
            a = R(i,1:3) - R0
            b = R(i+1,1:3) - R(i,1:3)
            cc = R0 - R(i+1,1:3)
            L_laplace_local(i) = L_laplace_local(i) - 1.0/4.0*sum(cc*b)/area_FV(i)/(area_VV/3.0)
            L_laplace_local(i+1) = L_laplace_local(i+1) - 1.0/4.0*sum(a*b)/area_FV(i)/(area_VV/3.0)
        enddo
        L_laplace_local(1) = L_laplace_local(1) + L_laplace_local(N_T+1)
    end subroutine Laplace_Local_sub

    subroutine Laplace_global_Cell(C, k)
        implicit none
        type(cel)::C(:)
        integer, intent(in)::k
        integer i, j, N_t
        real*8 R0(1:3), R(1:20,1:3)
        real*8 area_FV(1:20), area_VV
        real*8 L_laplace_local(1:20)

        C(k)%L_glb = 0.0
        do i = 1, C(k)%N_V
            area_FV = 0.0
            R = 0.0
            R0 = C(k)%rv(i,1:3)
            N_t = C(k)%N_V_V(i)
            do j = 1, N_t
                R(j,1:3) = C(k)%rv(C(k)%V_V(i,j), 1:3)
                area_FV(j) = C(k)%area_F(C(k)%F_V(i,j))
            enddo
            area_VV = C(k)%area_V(i)
            call Laplace_Local_sub(R0, R, area_FV, area_VV, L_laplace_local, N_t)
            C(k)%L_glb(i,1:N_t) = L_laplace_local(1:N_t)
        enddo
    end subroutine Laplace_global_Cell


! ========================================================
! Phase field: interface energy for one triangle
! Ported from VesPhase/Basic_parameter.f90 Get_inter_phase
! ========================================================

    subroutine Get_inter_phase_local(R, area_local, fit, E_inter, b)
        implicit none
        real*8 R(1:3,1:3), area_local, fit(1:3)
        real*8 tr(1:3,1:3), dfi(1:3), E1, E_inter, b

        E1 = 0.0
        tr(1,1:3) = R(2,1:3) - R(1,1:3)
        tr(2,1:3) = R(3,1:3) - R(2,1:3)
        tr(3,1:3) = R(1,1:3) - R(3,1:3)
        dfi(1) = fit(2) - fit(1)
        dfi(2) = fit(3) - fit(2)
        dfi(3) = fit(1) - fit(3)

        E1 = E1 + b/24.0*(dfi(1)**2*(sum(tr(2,1:3)**2)+sum(tr(3,1:3)**2)) &
            - 2*dfi(1)*dfi(2)*sum(tr(1,1:3)*tr(2,1:3)))/area_local
        E1 = E1 + b/24.0*(dfi(2)**2*(sum(tr(1,1:3)**2)+sum(tr(3,1:3)**2)) &
            - 2*dfi(2)*dfi(3)*sum(tr(2,1:3)*tr(3,1:3)))/area_local
        E1 = E1 + b/24.0*(dfi(3)**2*(sum(tr(1,1:3)**2)+sum(tr(2,1:3)**2)) &
            - 2*dfi(3)*dfi(1)*sum(tr(3,1:3)*tr(1,1:3)))/area_local
        E_inter = E1
    end subroutine Get_inter_phase_local


! ========================================================
! Phase energy: GL double-well + interface gradient
! ========================================================

    subroutine phase_energy_Cell(C, k)
        implicit none
        type(cel)::C(:)
        integer, intent(in)::k
        integer i, i_t(1:3)
        real*8 E1, E2, fit(1:3), R(1:3,1:3), E_inter

        E1 = 0.0
        E2 = 0.0

        ! Interface gradient energy (per face)
        do i = 1, C(k)%N_F
            i_t = C(k)%V_F(i,1:3)
            fit(1:3) = C(k)%fi(i_t(1:3))
            R(1,1:3) = C(k)%rv(i_t(1),1:3)
            R(2,1:3) = C(k)%rv(i_t(2),1:3)
            R(3,1:3) = C(k)%rv(i_t(3),1:3)
            call Get_inter_phase_local(R, C(k)%area_F(i), fit, E_inter, C(k)%b_ph)
            E1 = E1 + E_inter
        enddo

        ! GL double-well energy (per vertex)
        do i = 1, C(k)%N_V
            if(C(k)%IS_BC_V(i).eq.0)then
            E2 = E2 + (-C(k)%a2_ph/2.0*C(k)%fi(i)**2 + C(k)%a4_ph/4.0*C(k)%fi(i)**4) &
                * C(k)%area_V(i)/3.0
            endif
        enddo

        C(k)%energyPh = E1 + E2
    end subroutine phase_energy_Cell


! ========================================================
! Initialize phase field: fi ~ 0 + 0.1 fluctuation
! Conservation: sum(fi * area_V/3) / area_tot = 0
! ========================================================

    subroutine initial_phase_cell(C, k, seed1)
        implicit none
        type(cel)::C(:)
        integer, intent(in)::k
        integer seed1
        integer i
        real*8 area_tot, sum1

        do i = 1, C(k)%N_V
            C(k)%fi(i) = (0.5 - random1(seed1)) * 0.2   ! range [-0.1, 0.1]
        enddo

        ! Enforce conservation (area-weighted mean = 0)
        area_tot = sum(C(k)%area_V(1:C(k)%N_V)) / 3.0
        sum1 = 0.0
        do i = 1, C(k)%N_V
            sum1 = sum1 + C(k)%fi(i) * C(k)%area_V(i) / 3.0 / area_tot
        enddo
        do i = 1, C(k)%N_V
            C(k)%fi(i) = C(k)%fi(i) - sum1
        enddo

        C(k)%order_fi = 0.0    ! fi_A = fi_B => order = 0
    end subroutine initial_phase_cell


! ========================================================
! Cahn-Hilliard dynamics: update fi
! Coupling: F_bend = (kpp + kpp1*fi)/2 * (H + c0 - c1*fi)^2
! mu_bend = kpp1/2*(H+c0-c1*fi)^2 + (kpp+kpp1*fi)*(H+c0-c1*fi)*(-c1)
! ========================================================

    subroutine update_fi_cell(C, k, dt_fi)
        implicit none
        type(cel)::C(:)
        integer, intent(in)::k
        real*8, intent(in)::dt_fi

        integer i, j, i1, N_t
        real*8 r_t, H_eff, kpp_local
        real*8 Delt_fi(1:20000)
        real*8 area_tot, sum1

        ! Step 1: Compute chemical potential at each vertex
        Delt_fi = 0.0
        do i = 1, C(k)%N_V
            if(C(k)%IS_BC_V(i).eq.1) cycle
            N_t = C(k)%N_V_V(i)
            r_t = 0.0

            ! GL double-well derivative
            r_t = r_t - C(k)%a2_ph * C(k)%fi(i) + C(k)%a4_ph * C(k)%fi(i)**3

            ! Bending-phi coupling derivative
            H_eff = C(k)%H_V(i) + C(k)%c0_phi - C(k)%c1_phi * C(k)%fi(i)
            kpp_local = C(k)%H_modulus + C(k)%kpp1_phi * C(k)%fi(i)
            r_t = r_t + C(k)%kpp1_phi / 2.0 * H_eff**2
            r_t = r_t + kpp_local * H_eff * (-C(k)%c1_phi)

            ! Interface gradient: -b * Lap(fi)
            do j = 1, N_t
                i1 = C(k)%V_V(i,j)
                r_t = r_t - C(k)%b_ph * (C(k)%fi(i1) - C(k)%fi(i)) * C(k)%L_glb(i,j)
            enddo

            Delt_fi(i) = r_t
        enddo

        ! Step 2: fi += L_fi * dt * Lap(mu)
        do i = 1, C(k)%N_V
            if(C(k)%IS_BC_V(i).eq.1) cycle
            N_t = C(k)%N_V_V(i)
            r_t = 0.0
            do j = 1, N_t
                i1 = C(k)%V_V(i,j)
                r_t = r_t + (Delt_fi(i1) - Delt_fi(i)) * C(k)%L_glb(i,j)
            enddo
            C(k)%fi(i) = C(k)%fi(i) + C(k)%L_fi * dt_fi * r_t
        enddo

        ! Step 3: Enforce conservation
        area_tot = sum(C(k)%area_V(1:C(k)%N_V)) / 3.0
        sum1 = 0.0
        do i = 1, C(k)%N_V
            sum1 = sum1 + C(k)%fi(i) * C(k)%area_V(i) / 3.0 / area_tot
        enddo
        do i = 1, C(k)%N_V
            C(k)%fi(i) = C(k)%fi(i) - (sum1 - C(k)%order_fi)
        enddo
    end subroutine update_fi_cell


end module basicfun_mod
