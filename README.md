# Turek–Hron FSI benchmark with a flag made of rigid segments

OpenFOAM case for the fluid–structure interaction benchmark of Turek and Hron (2006):
laminar flow around a fixed cylinder with an elastic flag attached to it in a channel.
The FSI2 and FSI3 variants are implemented.

The elastic flag is represented by a chain of N rigid bodies (`rigidBodyDynamics` module)
connected by `Rz` joints with linear torsional springs. The springs are provided by the
`linearJointSpring` restraint (`src/linearJointSpring`). The flow is solved with `pimpleFoam`;
the mesh motion is handled by `rigidBodyMotionSolver`. The problem is two-dimensional: the mesh
has one cell layer in the z direction and a depth of 1 m, so forces and masses are given per
metre of width, as in the benchmark.

This repository accompanies the paper

> A. N. Nuriev, B. Affane, A. Pozdeeva, *The rigid finite element method for fast fluid–structure interaction of slender
> bodies in OpenFOAM* (submitted).

## Requirements

- OpenFOAM (ESI, openfoam.com), tested with v2112 (used for the paper) and v2412. Other ESI
  versions with `rigidBodyDynamics` are expected to work.
- Python 3 with numpy and matplotlib (case generator and post-processing).

## Building the library with the linearJointSpring restraint

```sh
cd src/linearJointSpring
wclean
wmake libso
```

The library `librigidBodyDynamicsSpring.so` is installed to `$FOAM_USER_LIBBIN` and is loaded
in `constant/dynamicMeshDict` through `motionSolverLibs`. The OpenFOAM sources are not modified.

The restraint is applied to the body given in `body` and acts on the single degree of freedom
of its joint with the parent body:

```
hinge1
{
    type            linearJointSpring;
    body            segment1;
    stiffness       100;       // [N m/rad]
    damping         0.1;       // [N m s/rad]
    referenceAngle  0;         // [rad], optional
}
```

For a beam of bending stiffness `EI` divided into segments of length `dx`, the clamp hinge has
`stiffness 2EI/dx` and the interior hinges `EI/dx`.

## Workflow

The case stored in the repository is FSI2 on mesh M1 (`--case FSI2 --n-segments 20
--cells-per-segment 6`), so steps 2–3 can be run directly.

1. Generate the case files:

   ```sh
   python3 scripts/make_case.py --case FSI3 --n-segments 20 --cells-per-segment 6
   ```

   The script writes `system/blockMeshDict`, `system/topoSetDict`, `system/createPatchDict`,
   `system/controlDict`, `system/fvSolution`, `constant/dynamicMeshDict`,
   `constant/transportProperties`, `0.orig/U` and `caseInfo.txt` (read by
   `postprocess.py` to select the reference values). Before writing, the blockMesh topology is
   checked: block orientation, matching of edges and faces, assignment of faces to patches.

   - `--case`: test case (FSI2 or FSI3)
   - `--n-segments`: number of flag segments N
   - `--cells-per-segment`: number of cells along the flag per segment
   - `--box-segments`: number of segments inside the block around the cylinder (default 5)
   - `--n-outer`, `--relax`: outer PIMPLE iterations and acceleration relaxation
     (defaults: FSI2 6 and 0.5, FSI3 12 and 0.1)
   - `--eta`: hinge damping time, s (`damping = eta*stiffness`, default 0)

2. Run the case:

   ```sh
   ./Allrun
   ```

   The script runs `blockMesh`, `topoSet`, `createPatch`, `checkMesh`, then `pimpleFoam` in
   parallel on 2 processes (`system/decomposeParDict`), `reconstructPar` and the post-processing.

3. Post-processing:

   ```sh
   python3 scripts/postprocess.py --window 2
   ```

   For the last `--window` seconds the mean value, amplitude and frequency of ux, uy of point A,
   drag and lift are printed next to the reference values of Turek and Hron (2006), plots are
   saved to `turek_results.png`. The paper compares FSI2 and FSI3 also with later refined
   reference values, so the deviations printed here differ slightly from those in the paper.

## Meshes and settings of the paper

| Mesh | Options | Cells |
|---|---|---|
| M1 | `--n-segments 20 --cells-per-segment 6` | 18 860 |
| M2 | `--n-segments 20 --cells-per-segment 8` | 34 012 |
| M3 | `--n-segments 20 --cells-per-segment 12` | 76 206 |
| M4 | `--n-segments 30 --cells-per-segment 6` | 41 807 |

Coupling settings: FSI2 — `dt = 1e-3`, 6 outer iterations, relaxation 0.5;
FSI3 — `dt = 5e-4`, 12 outer iterations, relaxation 0.1 (the defaults of `make_case.py`).
Mesh motion: `displacementLaplacian` with `inverseDistance` diffusivity; Newmark integration
(γ = 0.5, β = 0.25), `accelerationDamping 1`, no structural damping.

## License

GNU General Public License v3.0 or later (see `LICENSE`), as required for code linked with OpenFOAM.

## Citation

If you use this code, please cite the paper above (see `CITATION.cff`).
