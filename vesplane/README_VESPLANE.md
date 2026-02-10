# 0Main_VesiclePlane.f90 — Main Program Documentation

## Purpose

`0Main_VesiclePlane.f90` is the main entry point for the SimCel simulation. It models the quasi-static indentation of a deformable spherical vesicle being pushed down onto a flat elastic membrane, computing bending energies, contact area fractions, and adhesion energies at each stage.

## Dependencies

The program uses a single umbrella module:

```fortran
use head_Cell_mod
```

This re-exports all simulation modules: `CelType_mod`, `allocCel_mod`, `basicfun_mod`, `read_mod`, `Manipulate_Cell_mod`, `dyn_mod`, `Surf_Operation_mod`, `force_mod`, and `point_mod`.

## Program Sections

### SEC. 0 — Parameter Definitions (lines 8-75)

Declares all local variables. Key groups:

| Variable Group | Examples | Purpose |
|----------------|----------|---------|
| Cells | `Cels(:)`, `Cout`, `C_AVE(1:2)`, `Vm` | The two membrane objects, combined output cell, averaging buffers, per-vertex force container |
| Adhesion/repulsion | `eps` (=0.8), `sig_eps` (=0.28) | Lennard-Jones-style adhesion energy and interaction range |
| Dynamics control | `f_rescale` (=0.8), `kBT`, `dt`, `dtk(1:2)` | Force scaling, thermal energy, timesteps |
| Z push-down | `Z_Stop(1:50)`, `t1` (=4000), `t2` (=1000) | Target Z positions and step counts per phase |
| Tension | `lam0` (=0.5), `eta` | Base tension, elastic-tension ratio |
| RNG | `seed1` | Random number generator seed |

### SEC. I — Initialization (lines 121-267)

**I-1: Define cells** — Allocates a 2-cell system:
- `Cels(1)` = planar membrane (open), loaded from `CIR1.dat`
- `Cels(2)` = spherical vesicle (closed), loaded from `R642.dat`

**I-2: Read meshes** — Reads mesh topology and vertex positions from Surface Evolver `.dat` files in `data/initials/`. Also allocates averaging buffers `C_AVE(k)`.

**I-3: Position and scale** — Centers both cells at the origin, scales the vesicle to effective radius `R0`, and places it above the plane at `Z = R0 + 0.5`. Generates 50 target Z-positions with step size `0.2*R0`.

**I-4: Compute geometry** — Calls `Get_shape_information` to compute edge lengths, face areas, vertex curvatures, and total area/volume. Stores reference values (`area_tot_zero`, `Vol_total_zero`).

**I-5: Set material properties** — For each cell, sets:
- `H_modulus` — bending modulus (0.375 for vesicle, 1.0 for plane)
- `kpp_V` — volume constraint stiffness
- `kpp_area` — area constraint stiffness
- `lamda` — surface tension
- `H_zero` — spontaneous curvature (set to 0)
- `eps`, `sig_eps` — adhesion parameters

**I-6: Optional read from checkpoint** — When `IS_READ=1`, re-reads pre-computed states from `RPOT.dat` (plane) and `RCOT.dat` (vesicle), rebuilds topology and geometry. Combines cells into `Cout` and saves initial snapshot.

**I-7: Build neighbor lists** — Constructs near-vertex and near-face lists between the two cells using distance threshold `D_thres=3.5`, needed for repulsion force computation.

### SEC. II — Dynamics Loop (lines 297-479)

The main simulation loop runs `NN = N_V2 * (t1+t2) * 100` total steps, organized into 20 Z-stages:

**Z-stage structure** (per stage `k_Z`):
- Phase 1 (`t1=4000` steps): Linearly interpolate `Z_FIXED` from `Z_STOP(k_Z)` to `Z_STOP(k_Z+1)`
- Phase 2 (`t2=1000` steps): Hold position fixed, relax. Last 200 steps accumulate contact fraction average.

**Per-step vertex update** (repeated `N_V1 + N_V2` times per step):

1. **Select cell** — Randomly pick cell 1 or 2, weighted by vertex count
2. **Select vertex** — Random interior vertex (skip boundary and next-to-boundary vertices)
3. **Compute forces:**
   - `Point_H_Force_Cell` — Numerical bending force (finite difference, delta=1e-6)
   - `Point_PV_Force_Cell` — Analytical pressure + area constraint force (vesicle only)
   - `Update_NearF_Point_Rep` — LJ-type adhesion/repulsion from nearby faces
4. **Move vertex** — If `|F| * dt < 0.05` (stability guard), update position: `rv += dt * F * (N_V_k / N_V2)`
5. **Local update** — `Point_update_cell` recomputes geometry in the 1-ring neighborhood

**Every 20 steps:**
- `Vertex_Averaging_Cell_Origin` — Laplacian-like mesh smoothing for both cells
- Rebuild near-face lists via `Get_NearF_Cells`
- Recompute full geometry via `Get_shape_information`
- Accumulate contact fraction during relaxation phase

**Every 200 steps:**
- Print to stdout: curvature energies, contact fractions, adhesion energies, Z-position
- Save snapshots: combined mesh (`PLN`), plane (`POT`), vesicle (`COT`)

**At each Z-stage boundary:**
- Record stage results to `fracALL`, `EALL` arrays
- Write cumulative energy file to `data/EnergyH5.dat`

## Output Files

| File | Contents |
|------|----------|
| `data/RPLN{kk}{k_Z}.dat` | Combined plane+vesicle mesh snapshot |
| `data/RPOT{kk}{k_Z}.dat` | Plane-only mesh snapshot |
| `data/RCOT{kk}{k_Z}.dat` | Vesicle-only mesh snapshot |
| `data/EnergyH5.dat` | Per-stage summary: stage index, Z-position, contact fraction, plane bending energy, vesicle bending energy |

## Stdout Output

Every 200 steps, the program prints:
```
PLN  <step>  <outer_loop>
<R0>  center: <vesicle_Z_centroid>
<total_area>  <sum_vertex_areas>
curvature energies(soft particle, plane, all)  <E_plane>  <E_vesicle>  <E_total>
contact fraction:  <frac_plane>  <frac_vesicle>
adhesion:  <E_rep_plane>  <E_rep_vesicle>  <E_rep_per_contact_area>
Z_pos:  <mean_vesicle_Z>
```

## Key Hardcoded Parameters

| Parameter | Value | Description |
|-----------|-------|-------------|
| `eps` | `0.01 * 80 = 0.8` | Adhesion energy scale |
| `sig_eps` | `0.28` | Adhesion interaction range |
| `f_rescale` | `0.8` | Force rescaling factor |
| `kBT_in` | `0.05` | Thermal energy input (scaled by R0) |
| `et_in` | `2.00` | Elastic-tension ratio input |
| `lam0` | `0.5` | Base surface tension |
| `D_thres` | `3.5` | Distance threshold for neighbor lists |
| `dtk(1:2)` | `0.00125` | Timestep for both cells |
| `t1` | `4000` | Steps per push-down phase |
| `t2` | `1000` | Steps per relaxation phase |
| `seed1` | `26345678` | RNG seed (seed_char='1') |

## Modifying the Simulation

- **Change initial geometry:** Edit `Cels(1)%FILE_Cel` and `Cels(2)%FILE_Cel` filenames (line 136-137), provide corresponding `.dat` files in `data/initials/`
- **Change adhesion strength:** Modify `eps` (line 78) and `sig_eps` (line 79)
- **Change push-down schedule:** Adjust the Z-step formula in the loop at line 193, or change `t1`/`t2` at line 196
- **Change material stiffness:** Modify `H_modulus`, `kpp_V`, `kpp_area` in the SEC. I-5 blocks (lines 209-220 for vesicle, lines 238-249 for plane)
- **Add a new force:** Compute the force on vertex `VM` after line 382, add it to the displacement calculation at lines 390-391
