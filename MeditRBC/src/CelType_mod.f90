module CelType_mod
  implicit none
  real*8 Pai
  parameter(Pai = 3.14159265358979d0)
  type Cel

    character*20 File_MS,FILE_Cel
  ! ######## Mesh
    integer N_V,N_E,N_F,N,IS_OPEN,N_BC_E,N_BC_V
    integer, allocatable::V_E(:,:),E_F(:,:),V_F(:,:),E_V(:,:),F_V(:,:),F_E(:,:),F_F(:,:) &
     ,V_V(:,:), E_V_opposite(:,:),V_E_opposite(:,:),N_F_V(:),N_E_V(:),N_F_E(:),N_V_V(:) &
     ,N_E_V_opposite(:), F_F_around(:,:),N_F_F_around(:),IS_BC_E(:),IS_BC_V(:),E_BC(:),V_BC(:) &
    ,IS_NEXTBC_V(:),IS_BC_F(:)
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
  ! ###### Membrane Elastic Prop (isotropic bending)
    real*8 H_modulus, lamda, m0,H_zero, P
    real*8 kpp_V
    real*8 pi_K_2_A0
    real*8 kpp_Area
    real*8 integral_H_A

    real*8 eps,sig_eps
    real*8 contarea

  ! ###### Shear (MS) parameters
    integer HAS_REFERENCE
    real*8 kpp_alpha
    real*8 mu_ms
    real*8 a3_ms, a4_ms
    real*8 b1_ms, b2_ms

  ! ###### Phase field phi(r) (concentration of oriented proteins)
    real*8, allocatable:: fi(:)        ! phi on vertices
    real*8, allocatable:: L_glb(:,:)   ! discrete Laplacian weights (N_V, 20)
    real*8 order_fi                    ! conserved order parameter
    real*8 L_fi                        ! CH mobility
    real*8 b_ph, a2_ph, a4_ph         ! Ginzburg-Landau parameters
    real*8 energyPh                    ! phase energy

  ! ###### Curvature tensor (principal curvatures and directions)
    real*8, allocatable:: K_V(:)       ! Gaussian curvature at each vertex (angle deficit)
    real*8, allocatable:: c1_V(:)      ! 1st principal curvature at each vertex
    real*8, allocatable:: c2_V(:)      ! 2nd principal curvature at each vertex
    real*8, allocatable:: D_V(:)       ! deviatoric curvature = (c1-c2)/2 = sqrt(H^2-K)
    ! deviatoric curvature tensor components in local frame:
    !   D_tensor = D * [[cos2beta, sin2beta],[sin2beta, -cos2beta]]
    real*8, allocatable:: dd1_V(:)     ! D*cos(2*beta)  component of deviatoric tensor
    real*8, allocatable:: dd2_V(:)     ! D*sin(2*beta)  component of deviatoric tensor
    ! local tangent frame per vertex
    real*8, allocatable:: T1_V(:,:)    ! 1st tangent basis vector (N_V,3)
    real*8, allocatable:: T2_V(:,:)    ! 2nd tangent basis vector (N_V,3)

  ! ###### Q-tensor orientation field (nematic order in tangent plane)
    !   Q = S * [[cos2alpha, sin2alpha],[sin2alpha, -cos2alpha]] / 2
    !   Stored as q1 = S*cos(2alpha), q2 = S*sin(2alpha)
    !   S = sqrt(q1^2+q2^2), alpha = atan2(q2,q1)/2
    real*8, allocatable:: q1(:)        ! Q-tensor component 1 per vertex
    real*8, allocatable:: q2(:)        ! Q-tensor component 2 per vertex
    real*8 M_Q                         ! mobility for Q-tensor relaxation
    real*8 K_frank                     ! Frank elastic constant |grad Q|^2
    real*8 a_Q, a_Q4                   ! Maier-Saupe: -a_Q/2*S^2 + a_Q4/4*S^4

  ! ###### Anisotropic curvature-orientation coupling
    !  E_aniso = [kpp_u/2*(H-H0_u)^2 + kpp_uD/2*(D-D0_u*cos(2*theta))^2] * phi
    !  theta = angle between director u and 1st principal direction e1
    !  cos(2theta) = (q1*dd1 + q2*dd2) / (S*D)
    real*8 kpp_u                       ! mean-curvature-orientation modulus
    real*8 kpp_uD                      ! deviatoric-curvature-orientation modulus
    real*8 H0_u                        ! spontaneous mean curvature (anisotropic)
    real*8 D0_u                        ! spontaneous deviatoric curvature
    real*8 energyAniso                 ! anisotropic coupling energy
    real*8 energyOrient                ! orientational free energy (MS + Frank)

  ! ###### Energy
    real*8 EnergyP,EnergyH,energyL,energyMs, energy_tot,energy_rep
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
