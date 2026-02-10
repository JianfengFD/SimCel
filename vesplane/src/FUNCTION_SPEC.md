# LN_base_code: Function-Calling Specification

## Overview

**Domain:** Vesicle and cell membrane deformation dynamics using triangulated surface meshes.

**Current Problem:** Interaction between a closed spherical vesicle and an open planar membrane, governed by Helfrich bending energy with area/volume constraints and inter-membrane adhesion/repulsion.

**Language:** Fortran 90 | **Compiler:** gfortran -O2 | **Build:** make

---

## Architecture

```
0Main_VesiclePlane.f90       ← Main program (uses head_Cell_mod)
    │
head_Cell_mod.f90            ← Umbrella module (re-exports all)
    │
    ├── force_mod.f90        ← Forces: bending, PV, repulsion, local update
    ├── Surf_Operation_mod   ← Mesh ops: label, flip, smooth, refine, WKONTRI
    ├── PLN_mod.f90          ← Plane-specific (stub)
    ├── dyn_mod.f90          ← Dynamic parameter variation
    ├── read_mod.f90         ← I/O: Surface Evolver format
    ├── Manipulate_Cell_mod  ← Move, resize, combine, copy, neighbor search
    ├── basicfun_mod.f90     ← Core geometry: area, vol, curvature, energy
    ├── allocCel_mod.f90     ← Allocate/deallocate Cel arrays
    ├── CelType_mod.f90      ← Cel data type definition
    └── point_mod.f90        ← PtCell type (per-vertex force data)
```

Compilation order: point_mod → CelType_mod → allocCel_mod → basicfun_mod → read_mod → Manipulate_Cell_mod → dyn_mod → Surf_Operation_mod → force_mod → PLN_mod → head_Cell_mod → Main

---

## Core Data Types

### Cel (CelType_mod) — The Membrane

| Category | Fields | Description |
|----------|--------|-------------|
| Topology | `N_V, N_E, N_F, IS_OPEN` | Mesh counts and type |
| Primary connectivity | `V_E(N,2), E_F(N,3)` | Edge→vertices, Face→signed-edges |
| Derived connectivity | `V_F, E_V, F_V, F_E, F_F, V_V` | Built by label_system |
| Valence | `N_V_V, N_F_V, N_E_V` | Star sizes |
| Boundary | `IS_BC_V, IS_NEXTBC_V, IS_BC_E` | Boundary flags |
| Positions | `rv(N,3), rc(3)` | Vertex coords, centroid |
| Edge geom | `Len_E, Vector_E, th_E` | Lengths, vectors, dihedral angles |
| Face geom | `Area_F, Norm_F, Vector_F` | Areas, normals, area vectors |
| Vertex geom | `Area_V, Norm_V, H_V, Darea_V, Dvol_V` | Dual area, normals, curvature, variations |
| Global geom | `area_tot, Vol_total` | Total area and volume |
| Material | `H_modulus, H_zero, kpp_V, kpp_Area, lamda, P, eps, sig_eps` | Elastic & interaction params |
| Energy | `EnergyH, EnergyP, energyL, energy_tot, energy_rep` | All energies |
| Multi-cell | `NearV, NearF, N_NearV, N_NearF` | Neighbor lists |
| WKONTRI | `RF, e, Gup, Gdn, nF, thF, LenF` | Walk-on-triangle arrays |

### PtCell (point_mod) — Per-Vertex Forces

| Field | Description |
|-------|-------------|
| `V_move, rv_move(3)` | Vertex index and position |
| `PM_Force_Point(3)` | Bending force |
| `PV_Force_Point(3)` | Pressure + area force |
| `Rep_Force_Point(3)` | Repulsion force |
| `Darea_Point(20,3), Dvol_Point(3)` | Local area/volume derivatives |

---

## Module Functions

### read_mod — I/O

| Function | Signature | Description |
|----------|-----------|-------------|
| `read_cell_size` | `(C)` | Read N_V, N_E, N_F from data/initials/C%FILE_Cel |
| `read_cell` | `(C)` | Read full mesh: rv, V_E, E_F |
| `save_Cell` | `(C, kk, t, file_in)` | Write to data/R{file_in}{t}.dat |
| `print_screen_cell` | `(C(:), k, t)` | Print energy/area/vol to stdout |

### basicfun_mod — Geometry & Energy

| Function | Signature | Description |
|----------|-----------|-------------|
| `Get_shape_information` | `(C(:), k)` | Compute ALL geometry for C(k) |
| `Get_shape_sgl_information` | `(C)` | Same for single Cel |
| `Get_Cell_area` | `(C(:), k)` | Face areas and total area |
| `Get_Cell_Len` | `(C(:), k)` | Edge lengths and vectors |
| `Get_Cell_Vol` | `(C(:), k)` | Volumes |
| `Get_Cell_H` | `(C(:), k)` | Mean curvature at all vertices |
| `Get_Cell_Energy` | `(C(:), k)` | Bending + area + volume energies |
| `Get_active_points` | `(Cels, active_vs, N_active, IS_active, R0)` | Select active vertices |
| `Get_area_F` | `(rv, V_F, N_F, Area_F, Vector_F, Norm_F)` | Low-level face areas |
| `Get_Vol_F` | `(rv, V_F, N_F, Vol_F, Vol_total)` | Low-level volumes |
| `Get_H_V_new` | `(rv, N_V, N_F_V, V_V, F_V, E_V, ...)` | All-vertex curvature (angle-deficit) |
| `Get_H_local_new` | `(R0, R, N_V_V, H, area, len, th)` | Single-vertex curvature |
| `Curvature_Energy` | `(kappa, H_V, H0, A_V, N_V, BC, NBC) → real*8` | E = Σ κ/2(H-H0)²dA |
| `Pos_Point_Triangle` | `(R, R1, R2, R3, n, a, b, d, isON, isIN)` | Point-triangle position |
| `R1R2_fun` | `(R1(3), R2(3), Ic) → real*8` | Cross product component |
| `RR_fun` | `(R(3,3), Ic) → real*8` | Triangle area-vector component |
| `V6_fun` | `(R(3,3)) → real*8` | 3×3 determinant |
| `random1` | `(seed) → real*8` | Knuth RNG [0,1) |

### Surf_Operation_mod — Mesh Operations

| Function | Signature | Description |
|----------|-----------|-------------|
| **Topology** |||
| `label_cell` | `(C(:), k)` | Build all derived connectivity for C(k) |
| `label_Sgl_cell` | `(C)` | Same for single Cel |
| `label_system` | `(V_F, E_V, ..., IS_OPEN, ...)` | Core topology builder |
| **Simple remeshing** |||
| `Bond_flip_Cell` | `(C(:), k, I_flip)` | Delaunay edge flipping |
| `Vertex_averaging_Cell` | `(C(:), k)` | Laplacian smoothing |
| `Insert_APOINT` | `(RV, ..., num_F)` | 1-to-3 face split |
| `Mesh_grow_Cell` | `(C(:), k, t, seed1)` | Random mesh growth |
| **Advanced remeshing (WKONTRI)** |||
| `WKONTRI_Surf_scheme_Cell` | `(C,k,k1,Iflip,tt,type)` | Full surface remeshing |
| `Surf_scheme_Cell` | `(C,k,k1,Iflip)` | Simple smooth+flip |
| `Vertex_Averaging_Cell_Origin` | `(C1,C2,kFN,al,bet,type,t)` | Curvature-preserving averaging |
| **Refinement** |||
| `Refine_Cell` | `(C,k1,k2,n,method)` | Butterfly subdivision |
| **WKONTRI subsystem** |||
| `WKONTRI_alloc/dealloc` | `(C)` | Alloc/free walk arrays |
| `WKONTRI_Init` | `(C)` | Init face coords & metrics |
| `WKONTRI_MOVE` | `(C,v1,kF0,al0,bet0,r0,kF1,al1,bet1,rv1)` | Walk on surface |
| `WKONTRI_POS_TRI_POINT` | `(r00,r1-3,n0-3,r,d,al,bet,k)` | Point-triangle 3D position |

### force_mod — Forces & Local Update

| Function | Signature | Description |
|----------|-----------|-------------|
| **Bending force** |||
| `Point_H_Force_Cell` | `(C(:),k,V,intH)` | Numerical bending force on vertex V |
| `Point_H_FORCE_numerical` | `(Vm,F,RV,kappa,H0,HV,VV,NVV,AV,BC)` | Low-level numerical F_bend |
| **Pressure/area force** |||
| `Point_PV_Force_Cell` | `(C(:),k,V)` | Analytical P+area force (closed only) |
| `POINT_PV_FORCE` | `(vm,DA,DV,F,RV,VV,FV,NVV,VF,AF,P,lam)` | Low-level analytical |
| `POINT_PV_FORCE_Numerical` | `(vm,F,RV,VV,FV,NVV,VF,P,lam)` | Low-level numerical |
| **Repulsion** |||
| `Update_NearF_Point_Rep` | `(C,k1,Ntot,V,r00,r01,eta)` | Face-based LJ repulsion |
| `Update_NearF_Point_Rep_HARD_PLN` | `(C,k1,...,Rz,R0)` | Hard repulsion vs plane |
| `Update_NearF_Point_Rep_HARD_SPH` | `(C,...,Rz,R0,F,contf,...)` | Hard repulsion vs sphere |
| `Get_NearF_Cells` | `(C,k1,k2,D,R00,contf)` | Build near-face list |
| **Local update** |||
| `Point_Update_Cell` | `(C(:),k,V)` | Update geometry after vertex move |
| `Update_New_Parameters` | `(vm,rv_m,rv,...,H_V)` | Core local update routine |

### Manipulate_Cell_mod — Cell Operations

| Function | Signature | Description |
|----------|-----------|-------------|
| `Move_Cell` | `(C(:), k, Rc0(3))` | Translate so centroid = Rc0 |
| `Resize_Cell` | `(C(:), k, R0)` | Scale to effective radius R0 |
| `Combine_Cells` | `(C1, C2, C3)` | Merge: C3 = C1 ∪ C2 |
| `Copy_Cell` | `(C1, C2)` | Copy topology + positions |
| `Get_NearV_Cells` | `(C(:), k1, k2, D)` | Build vertex neighbor list |

---

## Main Program Flow

```
INIT:  Read PLN2.dat, R642.dat → alloc → label → Get_shape → set material
       Place vesicle at Z = R0 + 0.5
       Build neighbor lists (NearV, NearF)

DYNAMICS (loop over push-down stages k_Z = 1..20):
  Phase 1 (t1=4000 steps): Push vesicle center Z → Z_STOP(k_Z)
  Phase 2 (t2=1000 steps): Relax at fixed Z
  
  Per step, N_V1+N_V2 vertex moves:
    1. Pick cell k (prob ~ N_V(k))
    2. Pick vertex VM randomly
    3. F_bend = Point_H_Force_Cell (numerical, δr=1e-6)
    4. F_PV   = Point_PV_Force_Cell (analytical, closed only)
    5. F_rep  = Update_NearF_Point_Rep (LJ on near faces)
    6. If |F|*dt < 0.05: rv += dt * F * N_V(k)/N_V(2)
    7. Point_Update_Cell → Update_New_Parameters
  
  Every 20 steps:
    Vertex_Averaging_Cell_Origin (remesh both cells)
    Rebuild near-face lists
    Compute contact fraction

OUTPUT: Energy, contact fraction, shape snapshots
```

---

## Build & Run

```bash
cd LN_base_code
make              # Compile with -O2
./0Main_VesiclePlane   # Run
make clean        # Clean
```

## How to Extend

1. **New force:** Add `Point_XXX_Force_Cell(C,k,V)` in force_mod, include in main loop.
2. **New geometry:** Provide new .dat file, adjust main program initialization.
3. **New energy:** Add to `Get_Cell_Energy` or new function following `Curvature_Energy`.
4. **New mesh op:** Add to Surf_Operation_mod; always call `label_cell` after topology change.
5. **New constraint:** Add to the dynamics section, following PV_Force pattern.
