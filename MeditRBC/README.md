# MeditRBC — Mediterranean Thalassemia RBC Simulation

Simulates red blood cell (RBC) membrane deformation coupled with a **concentration field** phi(r) and a **nematic Q-tensor orientation field** (q1, q2), modeling the anisotropic protein distribution relevant to Mediterranean thalassemia.

## Physics Model

### Membrane Mechanics
- **Helfrich bending energy**: kpp/2 * (H - H0)^2
- **Skalak shear energy**: cytoskeleton elasticity from a separate reference shape (MS148)
- **Volume + area constraints**: Lagrange multiplier approach, reduced volume ~ 0.6

### Concentration Field phi(r)
phi is an order parameter: phi = phi_protein - phi_empty, so the actual protein concentration is **(1+phi)/2**.

- **Chemical potential**: mu = a_FH * atanh(phi) - chi * phi  (Flory-Huggins entropy + chi interaction)
- **Interface gradient**: b/2 * |nabla_s phi|^2 (cotangent-weight discrete Laplacian)
- **Dynamics**: Cahn-Hilliard (conserved): d phi/dt = L_phi * Delta_s(mu)

### Q-tensor Orientation Field
Nematic order on the surface, represented as q1 = S cos(2alpha), q2 = S sin(2alpha), where S is the scalar order parameter and alpha is the director angle in the local tangent frame (T1, T2).

- **Maier-Saupe**: -a_Q/2 * S^2 + a_Q4/4 * S^4
- **Frank elastic**: K_frank/2 * |nabla_s Q|^2
- **Dynamics**: Allen-Cahn (non-conserved): dq_a/dt = -M_Q * delta F / delta q_a

### Anisotropic Curvature-Orientation Coupling

```
E_aniso = (1+phi)/2 * [kpp_u/2 * (H - H0)^2
  + kpp_uD/2 * (D^2 - 2*D0*(q1*dd1 + q2*dd2) + D0^2*(1+S^2)/2)]
```

where:
- D = sqrt(H^2 - K) is the deviatoric curvature
- K is obtained via the angle deficit method: K_v = (2pi - sum angles) / (A_v/3)
- S^2 = q1^2 + q2^2 is the squared scalar order parameter
- q1*dd1 + q2*dd2 is a tensor contraction coupling orientation to curvature direction
- dd1, dd2 are deviatoric curvature tensor components from the Taubin shape operator
- Variational derivatives: dE/dq1 = (1+phi)/2 * kpp_uD/2 * (-2*D0*dd1 + D0^2*q1)

## Numerical Methods

| Quantity | Method |
|----------|--------|
| Mean curvature H | Dihedral angle formula |
| Gaussian curvature K | Angle deficit: K_v = (2pi - sum angles) / (A_v/3) |
| Deviatoric curvature D | sqrt(H^2 - K) |
| Principal directions | Taubin curvature tensor (2x2 shape operator in tangent frame) |
| Surface Laplacian | Cotangent-weight discrete Laplacian |
| Mechanical forces | Numerical finite difference (vertex perturbation delta = 1e-6) |
| phi dynamics | Cahn-Hilliard with surface Laplacian of chemical potential |
| Q dynamics | Allen-Cahn with covariant Laplacian (parallel transport) |

## File Structure

```
MeditRBC/
  0Main_MeditRBC.f90     Main program
  makefile                Build with gfortran -O3
  para.txt                14 input parameters
  visualize.py            Python visualization (3 modes)
  src/
    CelType_mod.f90       Extended Cel data type (K_V, D_V, dd1_V, dd2_V, q1, q2, ...)
    basicfun_mod.f90      Core physics: curvature tensor, energies, field evolution
    force_mod.f90         Mechanical forces including anisotropic coupling (numerical FD)
    allocCel_mod.f90      Memory allocation for all fields
    point_mod.f90         Per-vertex force data type
    read_mod.f90          Surface Evolver mesh reader (PM shape + MS reference)
    Surf_Operation_mod.f90  Mesh topology operations
    Manipulate_Cell_mod.f90 Cell labeling and manipulation
    dyn_mod.f90           Dynamics utilities
    head_Cell_mod.f90     Module aggregator
  data/
    initials/
      R10_AD_2562.dat     RBC shape (2562 vertices, Surface Evolver format)
      R2562_MS148.dat     MS reference shape (shear energy reference)
    RMRC00.dat            Shape output snapshots
    Rfi00.dat             phi field output
    RQ00.dat              Q-tensor field output (q1, q2, S, D)
```

## Parameters (para.txt)

```
25.0    kpp_alpha    Shear area modulus ratio
12.5    mu_ms        Shear modulus
1.0     b_ph         Interface gradient coefficient for phi
1.0     a_FH         Flory-Huggins entropy coefficient (mu += a_FH*atanh(phi))
1.0     chi_ph       Chi interaction coefficient (mu -= chi_ph*phi)
0.5     L_fi         Cahn-Hilliard mobility
1.0     kpp_u        Mean-curvature anisotropic modulus
2.0     kpp_uD       Deviatoric-curvature anisotropic modulus
0.0     H0_u         Spontaneous mean curvature (anisotropic)
0.5     D0_u         Spontaneous deviatoric curvature
1.0     a_Q          Maier-Saupe quadratic coefficient
1.0     a_Q4         Maier-Saupe quartic coefficient
0.5     K_frank      Frank elastic constant
0.5     M_Q          Q-tensor relaxation mobility
0.002   dt           Vertex dynamics time step
0.005   dt_fi        Cahn-Hilliard time step
0.005   dt_Q         Q-tensor relaxation time step
```

## Build & Run

```bash
make                              # compile
./0Main_MeditRBC para.txt         # run with parameter file
```

Requires `gfortran`.

## Visualization

```bash
# Mode 1: shape only
python visualize.py data/RMRC00.dat

# Mode 2: shape + concentration field coloring
python visualize.py data/RMRC00.dat --phi data/Rfi00.dat

# Mode 3: shape + concentration + director arrows
python visualize.py data/RMRC00.dat --phi data/Rfi00.dat --Q data/RQ00.dat

# Options
#   -o output.png        Output file name
#   --no-show            Save without displaying
#   --arrows 300         Number of director arrows
#   --elev 25 --azim -60 View angle
```

Requires `matplotlib` and `numpy`.

## Output Format

- **Shape** (RMRC*.dat): Surface Evolver format (vertices / edges / faces sections)
- **phi** (Rfi*.dat): `vertex_index  phi_value` per line
- **Q-tensor** (RQ*.dat): `vertex_index  q1  q2  S  D` per line

## Initialization

1. PM shape loaded from `R10_AD_2562.dat` (closed surface, 2562 vertices)
2. MS reference shape loaded from `R2562_MS148.dat` into rv0
3. Both shapes centered at origin and rescaled so that edge length ~ 1:
   R0 = 23.75 * sqrt(N_V / 10242), rv = R0 * rv / sqrt(A_tot / 4pi)
4. phi initialized with small random perturbations around 0
5. Q-tensor (q1, q2) initialized with small random perturbations
