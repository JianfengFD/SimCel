     
module Manipulate_Cell_mod
 use CelType_mod
 use basicfun_mod
 use allocCel_mod
 use point_mod
    contains
        subroutine Get_Cell_Center(C,k)
           implicit none
           integer, intent(in)::k
           type(Cel)::C(:)
           integer j,N_V
           real*8 rv_center(1:3)
           N_V=C(k)%N_V
           C(k)%rc(1)=sum(C(k)%rv(1:N_v,1))/N_V
           C(k)%rc(2)=sum(C(k)%rv(1:N_v,2))/N_V
           C(k)%rc(3)=sum(C(K)%rv(1:N_v,3))/N_V
        end subroutine Get_Cell_Center
        
        subroutine Move_Cell(C,k,Rc0)
           implicit none
           integer, intent(in)::k
           type(Cel)::C(:)
           real*8 Rc0(:)
           integer j,N_V
           real*8 rv_center(1:3)
           N_V=C(k)%N_V
           rv_center(1)=sum(C(k)%rv(1:N_v,1))/N_V
           rv_center(2)=sum(C(k)%rv(1:N_v,2))/N_V
           rv_center(3)=sum(C(K)%rv(1:N_v,3))/N_V
           do j=1,N_V
               C(k)%rv(j,1:3)=C(k)%rv(j,1:3)-rv_center+Rc0(1:3)
           enddo
        end subroutine Move_Cell
        
        subroutine Resize_Cell(C,k,R0)
            implicit none
            type(Cel)::C(:)
            integer, intent(in)::k
            real*8,intent(in)::R0
              call Get_Cell_area(C,k)
              C(k)%area_tot=sum(C(k)%area_F(1:C(k)%N_F))
              C(k)%rv=C(k)%rv*R0/(C(k)%area_tot/4/3.14159265)**0.5
        end subroutine Resize_Cell

        subroutine Combine_Cells(C1,C2,C3) ! C1+C2 = C3
            implicit none
            type(Cel)::C1,C2,C3
            integer N1,N2,N3
            integer i,j,k
            
            ! assign RV from C1 and C2 to C3
            C3%N_V=C1%N_V+C2%N_V
            C3%RV(1:C1%N_V,:) = C1%RV(1:C1%N_V,:)
            C3%RV(C1%N_V+1:C1%N_V+C2%N_V,:) = C2%RV(1:C2%N_V,:)
            ! set V_E
            C3%N_E=C1%N_E+C2%N_E
            C3%V_E(1:C1%N_E,1:2) = C1%V_E(1:C1%N_E,1:2) 
            C3%V_E(C1%N_E+1:C3%N_E,1:2) = C2%V_E(1:C2%N_E,1:2) +C1%N_V
            ! set E_F
            C3%N_F=C1%N_F+C2%N_F
            C3%E_F(1:C1%N_F,:) = C1%E_F(1:C1%N_F,:)
            
            do i = 1, C2%N_F
                do j=1,3
                    if(C2%E_F(i,j).gt.0)then
                        C3%E_F(C1%N_F+i,j) = C2%E_F(i,j)+C1%N_E
                    else
                        C3%E_F(C1%N_F+i,j) = C2%E_F(i,j)-C1%N_E
                    endif
                enddo
            enddo
            
            
            
        
        End subroutine Combine_Cells
        
        subroutine Combine_NCells(C,k1,k2,Cout)
            implicit none
            type(Cel)::C1,C2,Cout
            integer N1,N2,N3
            integer k1,k2,i,j
            type(Cel)::C(:)
            
            C1=C(k1);C2 = C(k1+1)
             call Combine_Cells(C1,C2,Cout)
            do i=k1+2,k2
                call dealloc_SglCel(C1)
                call dealloc_SglCel(C2)
                C1=Cout;C2=C(i)
                call dealloc_SglCel(Cout)
                call Combine_Cells(C1,C2,Cout)
            enddo
            
            
            
        End Subroutine Combine_NCELLS
        
        subroutine Copy_Cell(C1,C2)
        implicit none
        type(CEL)::C1,C2
        integer N
        !N=C1%N
        !call alloc_SglCel(C2,N)
        C2%N_V = C1%N_V
        C2%N_E = C1%N_E
        C2%N_F = C1%N_F
        C2%RV(1:C2%N_V,1:3) = C1%RV(1:C2%N_V,1:3)
        C2%V_E = C1%V_E
        C2%E_F = C1%E_F
        C2%IS_BC_V=C1%IS_BC_V
         C2%IS_NEXTBC_V=C1%IS_NEXTBC_V
        C2%IS_BC_E=C1%IS_BC_E
        
        End subroutine Copy_Cell
        
        
        subroutine Get_isNearIN(C,k)
            implicit none
            type(Cel)::C(:)
            integer k
            integer i,j,j1,j2,k1,k2
            C(k)%Is_NearIN = 0
            do i=1,C(k)%N_V
                do j=1, C(k)%N_V_V(i)
                    k1 = C(k)%V_V(i,j)
                    C(k)%IS_NearIn(i,j)=1
                    C(k)%IS_nearIn(j,i)=1
                    do j1=1,C(k)%N_V_V(k1)
                        k2 = c(k)%V_V(k1,j1)
                        C(k)%IS_NearIn(i,k2)=1
                        C(k)%IS_nearIn(k2,i)=1
                    enddo
                enddo
            enddo
            
        end subroutine Get_isNearIN
        
        subroutine Get_NearV_Cells(C,k1,k2,D_thres)
            implicit none
            type(Cel)::C(:)
            integer k1,k2
            integer i,j,k
            real*8 dsq
            real*8 D_thres,DD
            DD = D_thres**2
            if(k1.ne.k2)then
            do i = 1, C(k1)%N_V
               C(k1)%N_NearV(k2,i) = 0
               do j = 1, C(k2)%N_V
                   dsq = sum( (C(k1)%RV(i,1:3)-C(k2)%RV(j,1:3))**2 )
                   if(dsq.lt.DD)then
                       C(k1)%N_NearV(k2,i)=C(k1)%N_NearV(k2,i)+1
                       C(k1)%NearV(k2,i,C(k1)%N_NearV(k2,i) ) = j
                   endif
               enddo
            enddo
            endif
            
            if(k1.eq.k2)then
            do i = 1, C(k1)%N_V
               C(k1)%N_NearV(k2,i) = 0
               do j = 1, C(k2)%N_V
                   if(C(k1)%Is_NearIN(i,j).eq.0)then
                   dsq = sum( (C(k1)%RV(i,1:3)-C(k2)%RV(j,1:3))**2 )
                   if(dsq.lt.DD)then
                       C(k1)%N_NearV(k2,i)=C(k1)%N_NearV(k2,i)+1
                       C(k1)%NearV(k2,i,C(k1)%N_NearV(k2,i) ) = j
                   endif
                   endif
               enddo
               enddo
            endif
        end subroutine Get_NearV_Cells
         
        
        

end module Manipulate_Cell_mod
 