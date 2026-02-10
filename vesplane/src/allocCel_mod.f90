
module allocCel_mod
    use CelType_mod
    implicit none
    contains
        subroutine alloc_cel(C,k,N,Ntot)
          implicit none
          integer,intent(in)::k,N,Ntot
          type(Cel),dimension(:)::C
        ! %% mesh  
            allocate(C(k)%V_E(1:N,1:2));allocate(C(k)%E_F(1:N,1:3))
            allocate(C(k)%V_F(1:N,1:3));allocate(C(k)%E_V(1:N,1:40))
            allocate(C(k)%F_V(1:N,1:40));allocate(C(k)%F_E(1:N,1:40))
            allocate(C(k)%F_F(1:N,1:3));allocate(C(k)%N_E_V(1:N))
            allocate(C(k)%N_F_V(1:N));allocate(C(k)%N_F_E(1:N))
            allocate(C(k)%V_V(1:N,1:40));allocate(C(k)%N_V_V(1:N))
            allocate(C(k)%E_V_opposite(1:N,1:40));allocate(C(k)%N_E_V_opposite(1:N)) ! 14
            allocate(C(k)%V_E_opposite(1:N,1:40)) !15
        
            allocate(C(k)%RV(1:N,1:3),C(k)%rv0(1:N,1:3)) !17  
            allocate(C(k)%Vector_E(1:N,1:3) ,C(k)%Len_E(1:N),C(k)%th_E(1:N),&
                 C(k)%Len_E_Zero(1:N)) !21
            allocate(C(k)%Vector_F(1:N,1:3),C(k)%Area_F(1:N),C(k)%Norm_F(1:N,1:3),&
                 C(k)%Area_F_Zero(1:N),C(k)%Area_V(1:N),C(k)%Area_V_zero(1:N), &
                 C(k)%Darea_V(1:N,1:40,1:3),C(k)%Norm_V(1:N,1:3)) !29
            allocate(C(k)%Vol_F(1:N),C(k)%Dvol_V(1:N,1:3)) !31
            allocate(C(k)%H_V(1:N))
            allocate(C(k)%rho_chg(1:N)) !33
            allocate(C(k)%NearV(1:Ntot,1:C(k)%N_V,1:400),C(k)%N_NearV(1:Ntot,1:C(k)%N_V),C(k)%Is_nearIn(1:C(k)%N_V,1:C(k)%N_V)) !36
            allocate(C(k)%NearF(1:Ntot,1:C(k)%N_F,1:400),C(k)%N_NearF(1:Ntot,1:C(k)%N_F),C(k)%NearF_k0(1:Ntot,1:C(k)%N_F,1:400),&
                     C(k)%NearF_k(1:Ntot,1:C(k)%N_F,1:400)) ! 40
            allocate(C(k)%F_F_around(1:C(k)%N_F,1:120),C(k)%N_F_F_around(1:C(k)%N_F)) !42
            allocate(C(k)%IS_BC_V(C(k)%N_V),C(k)%IS_NEXTBC_V(C(k)%N_V),C(k)%IS_BC_E(C(k)%N_E), C(k)%V_BC(C(k)%N_V),&
                C(k)%E_BC(C(k)%N_E) ) 
            allocate(C(k)%IS_BC_F(C(k)%N_F))
            
            !47
            
        end subroutine alloc_cel
       
        subroutine alloc_SglCel(C,N)
          implicit none
          integer,intent(in)::N
          type(Cel)::C
        ! %% mesh  
            allocate(C%V_E(1:N,1:2));allocate(C%E_F(1:N,1:3))
            allocate(C%V_F(1:N,1:3));allocate(C%E_V(1:N,1:40))
            allocate(C%F_V(1:N,1:40));allocate(C%F_E(1:N,1:40))
            allocate(C%F_F(1:N,1:3));allocate(C%N_E_V(1:N))
            allocate(C%N_F_V(1:N));allocate(C%N_F_E(1:N))
            allocate(C%V_V(1:N,1:40));allocate(C%N_V_V(1:N))
            allocate(C%E_V_opposite(1:N,1:40));allocate(C%N_E_V_opposite(1:N))
            allocate(C%V_E_opposite(1:N,1:40))
        
            allocate(C%RV(1:N,1:3),C%rv0(1:N,1:3))  
            allocate(C%Vector_E(1:N,1:3) ,C%Len_E(1:N),C%th_E(1:N),&
                 C%Len_E_Zero(1:N))
            allocate(C%Vector_F(1:N,1:3),C%Area_F(1:N),C%Norm_F(1:N,1:3),&
                 C%Area_F_Zero(1:N),C%Area_V(1:N),C%Area_V_zero(1:N), &
                 C%Darea_V(1:N,1:40,1:3))
            allocate(C%Vol_F(1:N),C%Dvol_V(1:N,1:3))
            allocate(C%H_V(1:N))
            !allocate(C%rho_chg(1:N))
         
            
            allocate(C%F_F_around(1:N,1:120),C%N_F_F_around(1:N))
            allocate(C%IS_BC_V(C%N_V),C%IS_NEXTBC_V(C%N_V),C%IS_BC_E(C%N_E), C%V_BC(C%N_V),C%E_BC(C%N_E))
             allocate(C%IS_BC_F(C%N_F))
        end subroutine alloc_SglCel

         subroutine dealloc_SglCel(C)
          implicit none
          
          type(Cel)::C
        ! %% mesh  
            deallocate(C%V_E);deallocate(C%E_F)
            deallocate(C%V_F);deallocate(C%E_V)
            deallocate(C%F_V);deallocate(C%F_E)
            deallocate(C%F_F);deallocate(C%N_E_V)
            deallocate(C%N_F_V);deallocate(C%N_F_E)
            deallocate(C%V_V);deallocate(C%N_V_V)
            deallocate(C%E_V_opposite);deallocate(C%N_E_V_opposite)
            deallocate(C%V_E_opposite)
        
            deallocate(C%RV,C%rv0)  
            deallocate(C%Vector_E ,C%Len_E,C%th_E,&
                 C%Len_E_Zero)
            deallocate(C%Vector_F,C%Area_F,C%Norm_F,&
                 C%Area_F_Zero,C%Area_V,C%Area_V_zero, &
                 C%Darea_V)
            deallocate(C%Vol_F,C%Dvol_V)
            deallocate(C%H_V)
            !deallocate(C%NearV,C%N_NearV,C%Is_nearIn)
            deallocate(C%F_F_around,C%N_F_F_around)
            deallocate(C%IS_BC_V,C%IS_NEXTBC_V,C%IS_BC_E,C%V_BC,C%E_BC)
            deallocate(C%IS_BC_F)
        end subroutine dealloc_SglCel

        
end module allocCel_mod
 
