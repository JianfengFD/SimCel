# RBC — Red Blood Cell Membrane Simulation

## Overview

Non-modularized Fortran 90 simulation of a closed red blood cell membrane being deformed by a cylindrical capillary. The membrane model includes Helfrich bending energy, volume/area constraints, **Skalak-type shear (membrane skeleton) energy** computed against a reference shape, capillary confinement forces, and self-avoidance repulsion.

**Author:** Li JF (Dobb@bbs.fudan.edu.cn)

## Code Structure

The code uses a flat (non-module) architecture with `include`-style files. All arrays are declared with a fixed upper bound `N = 50000`.

```
RBC/
├── RBC_Main_cap_evolve.f90    # Main program (525 lines)
├── RBC_para_Fun.f90           # Geometry functions, RNG, math utilities (485 lines)
├── RBC_label.f90              # Topology labeling (176 lines)
├── RBC_Force.f90              # All force computations (1144 lines)
├── RBC_Cap.f90                # Capillary + rotation + self-avoidance (166 lines)
├── RBC_datainput.f90          # File I/O and string parsing (599 lines)
├── sqcapaxial.txt             # Auxiliary data
└── data/initials/
    ├── R10_AD_2562.dat        # Current shape: 2562-vertex RBC mesh
    └── R2562_MS148.dat        # Reference shape: 2562-vertex reference mesh
```

**Total:** ~3095 lines of Fortran across 5 source files.

## Key Differences from Vesplane

| Aspect | RBC | Vesplane |
|--------|-----|----------|
| Architecture | Flat subroutines, fixed-size arrays (`0:N`) | Module-based, allocatable `Cel` type |
| Membrane type | Closed only | Open and closed |
| Mesh | Static (no remeshing) | Dynamic (vertex averaging, edge flipping) |
| Reference shape | Yes (`rv0`, `Len_E_zero`, `Area_F_zero`) | Unused (`rv0` allocated but never set) |
| Shear energy | Skalak MS model | None |
| Bending force | Both numerical and **analytical** | Numerical only |
| External force | Cylindrical capillary (LJ wall) | Inter-membrane adhesion/repulsion |
| Thermal noise | Box-Muller Gaussian kicks | None |

## Physics Model

### Energy Terms

```
E_total = E_bending + E_pressure + E_area + E_shear
```

1. **Bending (Helfrich):** `E_H = Σ_v κ/2 (H - H₀)² dA + κ_π (∫H dA)²`
   - `H_modulus = 1.0` (bending modulus κ)
   - `H_zero = -10/(π √(A₀/4π))` (spontaneous curvature, nonzero for RBC)
   - `pi_K_2_A0 = κ/A₀` (integrated curvature constraint coefficient)

2. **Volume constraint:** `E_P = kpp_V/2 (V - V₀)²/V₀`
   - `kpp_V = 5×10⁴ × κ / (V₀/100)`

3. **Area constraint:** `E_A = kpp_area/2 (A - A₀)²/A₀`
   - `kpp_area = 5×10⁴ × κ / (A₀/140)`

4. **Shear / Membrane Skeleton (MS) energy:**
   - Per-face energy based on Skalak-type strain invariants α (area strain) and β (shear strain)
   - See section below

### Shear Energy — Skalak Model

For each triangular face with current edge lengths (a, b, c) and reference lengths (a₀, b₀, c₀):

**Strain invariants:**
```
S  = 2a²b² + 2a²c² + 2b²c² - a⁴ - b⁴ - c⁴       (= 16 × area²)
S₀ = same with reference lengths
Π  = a²b₀² + a²c₀² + b²c₀² + a₀²b² + a₀²c² + b₀²c² - a²a₀² - b₀²b² - c²c₀²

α = √(S/S₀) - 1       (local area dilation)
β = Π/√(S·S₀) - 1     (shear strain)
```

**Energy per face:**
```
E_MS(face) = [kpp_α/2 (α² + a₃α³ + a₄α⁴) + μ_ms(β + b₁αβ + b₂β²)] × A₀_face
```

**Default parameters:**
- `kpp_alpha = 25 × κ / (A₀/140)` — area dilation modulus
- `mu_ms = 0.5 × kpp_alpha` — shear modulus
- `a3_ms = -2.0, a4_ms = 8.0` — higher-order area dilation coefficients
- `b1_ms = 0.7, b2_ms = 0.75` — shear coupling coefficients

### Shear Force

The analytical force on vertex `v_move` from shear energy is computed in `Point_MS_force`. For each star face around the vertex:
1. Identify the three edges; find which edge is opposite to `v_move`
2. Get current lengths (a, b, c) and reference lengths (a₀, b₀, c₀)
3. Compute α, β, S, S₀ via `Get_alph_bet_analytic`
4. Compute `∂F/∂α` and `∂F/∂β` (energy derivatives w.r.t. strain invariants)
5. Compute `∂S/∂r₀` and `∂Π/∂r₀` (geometric derivatives w.r.t. vertex position)
6. Chain rule: `F_MS = -∂E/∂r₀ = -(∂F/∂α · ∂α/∂r₀ + ∂F/∂β · ∂β/∂r₀)`

## Source File Details

### RBC_para_Fun.f90 — Geometry & Utilities

| Function | Signature | Description |
|----------|-----------|-------------|
| `Get_Len_E` | `(rv, V_E, N_E, Len_E, Vector_E)` | Compute all edge lengths and vectors |
| `Get_area_local` | `(R1, R2, R3, area, vector)` | Single triangle area and area vector |
| `Get_area_F` | `(rv, V_F, N_F, Area_F, Vector_F, Norm_F)` | All face areas |
| `Get_vol_F` | `(rv, V_F, N_F, Vol_F, Vol_total)` | All face volumes |
| `Get_H_local_new` | `(R0, R, N_V_V, H, area, len, th)` | Single-vertex curvature |
| `Get_H_V_new` | `(rv, N_V, ...)` | All-vertex curvature (angle-deficit) |
| `Curvature_Energy` | `(κ, H_V, H₀, A_V, N_V) → real*8` | Total bending energy |
| `R1R2_fun` | `(R1, R2, Ic) → real*8` | Cross product component |
| `RR_fun` | `(R, Ic) → real*8` | Triangle area-vector component |
| `V6_fun` | `(R) → real*8` | 3×3 determinant (volume) |
| `random1` | `(seed) → real*8` | Knuth RNG |

### RBC_Force.f90 — Forces

| Function | Description |
|----------|-------------|
| `POINT_PV_FORCE` | Analytical pressure + area constraint force |
| `POINT_PV_FORCE_Numerical` | Numerical (FD) pressure + area force |
| `Point_MS_force` | **Analytical shear force** on single vertex |
| `Get_alph_bet` | Compute α, β and face area from edge lengths |
| `Get_alph_bet_analytic` | Compute α, β, S, S₀ (no area output) |
| `get_energy_MS` | Total shear energy summed over all faces |
| `MS_force_analytic` | Shear force on **all** vertices (batch) |
| `Point_PM_FORCE_NEW` | **Numerical** bending force (FD, delta=1e-5) |
| `Force_PM_Point_analysis` | **Analytical** bending force |
| `Update_New_Parameters` | Local geometry update after vertex move |

### RBC_Cap.f90 — Capillary & Auxiliary

| Function | Description |
|----------|-------------|
| `Capillary_Force` | LJ wall force from cylindrical capillary |
| `rotate_RBC` | Euler angle rotation of all vertices |
| `Max_XYZ` | Maximum radial extent in a given plane |
| `Get_nearV` | Build near-vertex list for self-avoidance |
| `Rep_Near_V_point` | LJ repulsion between nearby non-bonded vertices |

### RBC_datainput.f90 — I/O

| Function | Description |
|----------|-------------|
| `readfe` | Read current shape from `R10_AD_2562.dat` |
| `readfe_MS` | Read reference shape from `R2562_MS148.dat` |
| `char2Num` | Parse integer from string |
| `char2real` | Parse real (with E notation) from string |
| `char2real_No_E` | Parse real (no E notation) from string |

## Simulation Flow

1. **Read meshes:** `readfe` loads current shape `rv`, `readfe_MS` loads reference shape `rv0` (may be different `.dat` file)
2. **Center & rotate:** Both shapes centered at origin, optionally rotated (edge-on vs axial)
3. **Scale:** Both scaled to effective radius R₀ = 23.75 × √(N_V/10242)
4. **Label topology:** Build all connectivity arrays from V_E, E_F
5. **Compute reference geometry:** `Area_F_zero`, `Len_E_zero` from reference shape `rv0`
6. **Set parameters:** κ, kpp_V, kpp_area, kpp_alpha, mu_ms, H_zero, capillary radius
7. **Dynamics loop** (up to 9×10⁹ steps):
   - Pick random vertex
   - Compute bending force (analytical `Force_PM_Point_analysis`)
   - Compute pressure + area force (`POINT_PV_FORCE`)
   - Compute shear force (`Point_MS_force`) — uses `Len_E_zero` as reference
   - Compute capillary force
   - Compute self-avoidance repulsion (rebuild neighbor list every N_V×100 steps)
   - Move vertex: `rv += dt × (F_PV + F_PM + F_MS + F_cap + F_rep) / (A_V/3)`
   - Add stress term + optional thermal noise
   - Update local geometry (`Update_New_Parameters`)
8. **Output** (every 100000 steps): Save shape snapshot and contact statistics

## Build & Run

The code uses `include` directives (currently commented out). To compile:

```bash
cd RBC
gfortran -O3 -o RBC_cap RBC_para_Fun.f90 RBC_label.f90 RBC_datainput.f90 RBC_Force.f90 RBC_Cap.f90 RBC_Main_cap_evolve.f90
./RBC_cap
```

Requires data files in `data/initials/`.
