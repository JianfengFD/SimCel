# SimCel — Simulation of Cellular Membranes

A suite of Fortran 90 programs for simulating biomembrane mechanics on triangulated surfaces. Each subdirectory is a self-contained simulation targeting a specific biophysical problem.

## Modules

| Directory | Description | Key Features |
|-----------|-------------|--------------|
| **vesplane** | Vesicle-plane adhesion | Closed vesicle interacting with open planar membrane; LJ adhesion potential; multi-height scanning |
| **flatsurface** | Flat membrane mechanics | Open planar membrane with Helfrich bending + Skalak shear; shrinking reference shape |
| **flatsurfaceAB** | Phase separation on flat membrane | Cahn-Hilliard A-B phase dynamics coupled to curvature on a flat membrane |
| **RBC** | Red blood cell in capillary | Realistic RBC with cytoskeleton shear; aspiration into shrinking capillary; external stress |
| **VesPhase** | 3D vesicle phase separation | Closed vesicle with Cahn-Hilliard dynamics; adaptive mesh (bond flipping, vertex averaging) |
| **MeditRBC** | Mediterranean thalassemia RBC | RBC + concentration field phi + nematic Q-tensor orientation; anisotropic curvature-orientation coupling |

## Common Physics

All simulations share the same core framework:

- **Helfrich bending energy**: kpp/2 * (H - H0)^2 on triangulated surfaces
- **Mean curvature H**: computed via dihedral angles at edges
- **Volume / area constraints**: Lagrange multiplier approach
- **Dissipative vertex dynamics**: random vertex selection, force computation, overdamped motion
- **Surface Evolver format**: mesh I/O in .dat files (vertices / edges / faces)

## Shared Code Structure

Each module follows the same Fortran 90 module architecture:

```
project/
  0Main_*.f90              Main program (parameter setup, initialization, dynamics loop)
  makefile                 Build rules (gfortran -O3)
  para.txt                 External parameter file (where applicable)
  src/
    CelType_mod.f90        Cel data type (mesh, geometry, fields)
    basicfun_mod.f90       Geometry computation (area, volume, curvature, energy)
    force_mod.f90          Force computation (bending, shear, constraints)
    allocCel_mod.f90       Memory allocation
    read_mod.f90           Surface Evolver file reader
    Surf_Operation_mod.f90 Mesh topology operations
    Manipulate_Cell_mod.f90 Cell labeling and manipulation
    point_mod.f90          Per-vertex force container
    dyn_mod.f90            Dynamics utilities
    head_Cell_mod.f90      Module aggregator (use all)
  data/
    initials/              Initial mesh files
    R*.dat                 Output snapshots
```

## Build

Each module compiles independently:

```bash
cd vesplane && make        # or flatsurface, flatsurfaceAB, RBC, MeditRBC
```

Requires `gfortran`. No external libraries needed.

## Data Format

**Mesh files** (Surface Evolver .dat):
```
vertices
  1   x1   y1   z1
  2   x2   y2   z2
  ...
edges
  1   v1   v2   color 8
  ...
faces
  1   e1   e2   e3   color 6
  ...
```

## References

- Helfrich W. (1973) Elastic properties of lipid bilayers. Z. Naturforsch. 28c, 693-703
- Skalak R. et al. (1973) Strain energy function of red blood cell membranes. Biophys. J. 13, 245-264
