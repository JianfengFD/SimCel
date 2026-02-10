# CLAUDE.md — SimCel

## Project Overview

SimCel is a 3D computational physics simulation written in **Fortran 90** that models vesicle and planar membrane deformation dynamics. It simulates the interaction between a closed spherical vesicle and an open planar membrane using triangulated surface meshes, governed by Helfrich bending energy with area/volume constraints and inter-membrane adhesion/repulsion.

**Author:** Li JF (lijf@fudan.edu.cn), December 2014

## Build & Run

```bash
# Build (requires gfortran)
make

# Run the simulation
./0Main_VesiclePlane

# Clean build artifacts (.mod, .o, executable)
make clean
```

- **Compiler:** `gfortran -O3` (Fortran 90 with aggressive optimization)
- **Dependencies:** None — pure Fortran 90 standard library only
- **Build output:** `0Main_VesiclePlane` executable in the project root
- **Build artifacts:** `.mod` and `.o` files are generated in the project root directory

## Testing

There is no automated test suite. Validation is done by inspecting simulation output files in `data/`.

## Linting / Formatting

No linting or formatting tools are configured. Code is free-form Fortran 90.

## Project Structure

```
SimCel/
├── 0Main_VesiclePlane.f90    # Main program entry point
├── makefile                   # Build configuration
├── src/                       # Module source files
│   ├── point_mod.f90          # PtCell type — per-vertex force data
│   ├── CelType_mod.f90        # Cel type — core membrane data structure
│   ├── allocCel_mod.f90       # Memory allocation/deallocation for Cel
│   ├── basicfun_mod.f90       # Core geometry: area, volume, curvature, energy
│   ├── read_mod.f90           # I/O: reads/writes Surface Evolver format files
│   ├── Manipulate_Cell_mod.f90# Cell operations: move, resize, combine, copy
│   ├── dyn_mod.f90            # Dynamic parameter variation over time
│   ├── Surf_Operation_mod.f90 # Mesh topology, remeshing, WKONTRI surface walking
│   ├── force_mod.f90          # Physics forces: bending, pressure, repulsion
│   ├── head_Cell_mod.f90      # Umbrella module — re-exports all others
│   └── FUNCTION_SPEC.md       # Detailed API documentation for all modules
└── data/
    └── initials/              # Initial condition mesh files (.dat, Surface Evolver format)
        ├── CIR1.dat, CIR2.dat # Circular membrane configurations
        ├── PLN2.dat           # Planar membrane configuration
        ├── R300.dat, R642.dat # Spherical vesicle configurations
        ├── RCOT.dat           # Control vesicle
        └── RPOT.dat           # Test vesicle
```

## Architecture

### Module Dependency Chain

Compilation must follow this order (enforced by the makefile):

```
point_mod → CelType_mod → allocCel_mod → basicfun_mod → read_mod
  → Manipulate_Cell_mod → dyn_mod → Surf_Operation_mod → force_mod
  → head_Cell_mod → 0Main_VesiclePlane
```

`head_Cell_mod` is the umbrella module that re-exports all other modules. The main program uses only `head_Cell_mod`.

### Core Data Types

**`Cel`** (defined in `CelType_mod`): The primary data structure representing a triangulated membrane surface. Contains:
- Mesh topology: vertex/edge/face counts and connectivity arrays (`V_E`, `E_F`, `V_F`, etc.)
- Geometry: vertex positions (`rv`), edge lengths/vectors, face areas/normals, curvatures (`H_V`)
- Material properties: bending modulus (`H_modulus`), spontaneous curvature (`H_zero`), pressure (`P`), tension (`lamda`)
- Energies: bending (`EnergyH`), pressure (`EnergyP`), repulsion (`energy_rep`)
- Multi-cell interaction: neighbor vertex/face lists (`NearV`, `NearF`)

**`PtCell`** (defined in `point_mod`): Per-vertex force container holding bending force, pressure/area force, repulsion force, and local geometric derivatives.

### Key Modules

| Module | Purpose |
|--------|---------|
| `basicfun_mod` | Core geometry calculations: area, volume, curvature, energy |
| `force_mod` | All physics forces: bending (numerical), pressure+area (analytical), LJ repulsion |
| `Surf_Operation_mod` | Mesh operations: topology labeling, edge flipping, vertex smoothing, subdivision, WKONTRI walk-on-triangles |
| `read_mod` | File I/O in Surface Evolver format |
| `Manipulate_Cell_mod` | High-level cell operations: translate, scale, merge, copy |
| `allocCel_mod` | Memory management for Cel arrays |

### Simulation Flow

1. **INIT**: Read initial meshes (plane from `CIR1.dat`, vesicle from `R642.dat`), allocate, label topology, compute geometry, set material parameters
2. **DYNAMICS**: Multi-stage push-down loop (20 Z-positions):
   - Phase 1 (4000 steps): Push vesicle center downward toward target Z
   - Phase 2 (1000 steps): Relax at fixed Z position
   - Per step: randomly select vertex, compute forces (bending + pressure + repulsion), move vertex, update local geometry
   - Every 20 steps: remesh (vertex averaging), rebuild neighbor lists, compute contact fraction
3. **OUTPUT**: Energy and contact data saved to `data/` directory

## Code Conventions

- **Language standard:** Free-form Fortran 90
- **Naming:** Modules use `_mod` suffix; subroutine names use `CamelCase` (e.g., `Get_shape_information`, `Point_H_Force_Cell`)
- **Array indexing:** 1-based (Fortran default)
- **Real precision:** `real*8` (double precision) used throughout
- **Constants:** `PAI` (pi) defined as a parameter in `CelType_mod`
- **Connectivity convention:** Signed edge indices encode orientation in `E_F` arrays
- **Boundary handling:** `IS_BC_V`, `IS_BC_E`, `IS_NEXTBC_V` flags mark boundary vertices/edges for open membranes
- **After any topology change:** Always call `label_cell` to rebuild derived connectivity arrays
- **After moving vertices:** Call `Point_update_cell` for local geometry updates, or `Get_shape_information` for full recomputation
- **Simulation parameters** are hardcoded in the main program (no external config files)
- **Data files** use Surface Evolver `.dat` format

## How to Extend

1. **New force:** Add `Point_XXX_Force_Cell(C,k,V)` in `force_mod`, include in main dynamics loop
2. **New geometry:** Provide a new `.dat` file in `data/initials/`, adjust main program initialization
3. **New energy term:** Add to `Get_Cell_Energy` or create a new function following the `Curvature_Energy` pattern
4. **New mesh operation:** Add to `Surf_Operation_mod`; always call `label_cell` after topology changes
5. **New constraint:** Add to the dynamics section following the `PV_Force` pattern

## Important Notes

- No `.gitignore` exists — build artifacts (`.mod`, `.o`, executable) may be tracked. Run `make clean` before committing
- The `src/FUNCTION_SPEC.md` file contains detailed API documentation with function signatures for all 167 subroutines/functions
- Input data files in `data/initials/` are large binary-format mesh files (up to 2 MB)
- The simulation is single-threaded (no OpenMP/MPI parallelization)
