module CelType_mod
  implicit none
  real*8 Pai
  parameter(Pai = 3.14159265)
  type Cel
      
    character*8 File_MS,FILE_Cel
  ! ######## Mesh
    integer N_V,N_E,N_F,N,IS_OPEN,N_BC_E,N_BC_V
    integer, allocatable::V_E(:,:),E_F(:,:),V_F(:,:),E_V(:,:),F_V(:,:),F_E(:,:),F_F(:,:) &
     ,V_V(:,:), E_V_opposite(:,:),V_E_opposite(:,:),N_F_V(:),N_E_V(:),N_F_E(:),N_V_V(:) &
     ,N_E_V_opposite(:), F_F_around(:,:),N_F_F_around(:),IS_BC_E(:),IS_BC_V(:),E_BC(:),V_BC(:) &
    ,IS_NEXTBC_V(:),IS_BC_F(:)
   ! V_E(No. of Edge,two vertices): input (E)dge, output (V)ertex. V_E(3,1): the first vertex of the edge with edge label no=3.
   ! E_F(5,2): the second edge of the face with face no.=5.
   ! V_F(5,2): the second vertex of the face with face no.=5. 
   ! F_E(5,2): the second face that is with one edge that is with edge no.=5.
   ! F_F(5,2): the second face that connected to the face with face no.=5. there are three triangles that is adjacent to a given triangle
   ! V_V(5,2): the second vertex of the vertex with vertex no.=5. Normally there are about five or six v around a given v
   ! F_F_around: the faces that are surrounding a face.
  ! ####### Shape
    real*8  rc(1:3)
    real*8, allocatable:: rv(:,:),rv0(:,:)
  
    real*8, allocatable:: Vector_E(:,:) ,Len_E(:),th_E(:),Len_E_Zero(:)
    real*8 L0,detL,maxL,minL
    
    real*8, allocatable:: Vector_F(:,:),Area_F(:),Norm_F(:,:),Area_F_Zero(:),&
                        Area_V(:),Area_V_zero(:),Darea_V(:,:,:),Norm_V(:,:)
    real*8  area_tot_zero,area_tot,area_rel
  
    real*8 Vol_total,Vol_total_Zero, Vol_rel
    real*8, allocatable:: Vol_F(:),Dvol_V(:,:)
  
    real*8, allocatable:: H_V(:)
  ! ###### Membrane Elastic Prop
    real*8 H_modulus, lamda, m0,H_zero, P
    real*8 kpp_V
    real*8 pi_K_2_A0
    real*8 kpp_Area
    real*8 integral_H_A
    
    real*8 eps,sig_eps
    real*8 contarea
    
  ! ###### Energy
    real*8 EnergyP,EnergyH,energyL, energy_tot,energy_rep
    real energy(1:50000),energyHa(1:50000)
    
    ! for multi cell
    integer,allocatable:: NearV(:,:,:),N_NearV(:,:),Is_nearIn(:,:)
    integer,allocatable:: NearF(:,:,:),N_NearF(:,:),NearF_k0(:,:,:),NearF_k(:,:,:)
    
    
    real*8 D_thres    ! threshold for repulsion
    
    real*8,allocatable:: rho_chg(:)
    real*8 Z
        ! for wkontriangle
    real*8,allocatable::RF(:,:,:),Gup(:,:,:),Gdn(:,:,:),e(:,:,:)
    real*8,allocatable::nF(:,:,:),thF(:,:),LenF(:,:)
  
  end type Cel      

end Module Celtype_Mod
