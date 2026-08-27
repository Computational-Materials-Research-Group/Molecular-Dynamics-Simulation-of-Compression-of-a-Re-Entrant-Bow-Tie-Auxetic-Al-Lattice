# Molecular Dynamics Simulation of Compression of a Re-Entrant "Bow-Tie" Auxetic Al Lattice — EAM Potential

<p align="center">
  <img src="https://img.shields.io/badge/LAMMPS-MD%20Simulation-blue?style=for-the-badge&logo=gnu&logoColor=white"/>
  <img src="https://img.shields.io/badge/Process-Uniaxial%20Compression-orange?style=for-the-badge"/>
  <img src="https://img.shields.io/badge/EAM-Al99.eam.alloy-purple?style=for-the-badge"/>
  <img src="https://img.shields.io/badge/Structure-Re--entrant%20Bow--Tie-red?style=for-the-badge"/>
  <img src="https://img.shields.io/badge/Property-Negative%20Poisson's%20Ratio-yellow?style=for-the-badge"/>
  <img src="https://img.shields.io/badge/OVITO-Visualization-green?style=for-the-badge"/>
  <img src="https://img.shields.io/badge/Python-Geometry%20Generator-blue?style=for-the-badge&logo=python&logoColor=white"/>
</p>

<p align="center">
  A fully atomistic <b>molecular dynamics model of a 2D re-entrant "bow-tie" auxetic lattice</b>,
  built from FCC aluminum and compressed between rigid, no-slip grips using LAMMPS. The lattice
  is generated parametrically in Python -- a horizontal "waist" strut with inclined struts
  flaring <i>inward</i> from each end, tessellated into columns and tied together with flat
  loaded/support edges -- then filled with an FCC-Al atomistic mesh and verified, before any
  atomistic run, with a self-validated beam finite-element solver to confirm it actually has a
  <b>negative effective Poisson's ratio</b>, not just a bow-tie-shaped outline.
</p>

<p align="center">
  <b>This project went through several rounds of real, measured bugs -- not hypothetical
  edge cases.</b> Read <a href="#known-open-issue-unverified-strut-thickness-fix">Known Open
  Issue</a> and <a href="#common-errors-and-fixes">Common Errors and Fixes</a> before trusting
  any output from this pipeline; several fixes below were only confirmed correct by checking
  actual numbers (a validated FE solver, a real LAMMPS neighbor-list cutoff, a real console
  atom-count match), including at least one investigated hypothesis that turned out to be wrong.
</p>

<img width="1600" height="1200" alt="auxetic 2" src="https://github.com/user-attachments/assets/3b538eab-e609-4f0b-a6e9-7e45f73ae0c1" />

---

## Known Open Issue: Unverified Strut-Thickness Fix

**Current status: struts were found to reconstruct visibly during plain energy minimization
(before any compression is applied), and the fix for this has not been empirically confirmed.**

**Why it happens:** the original strut diameter (5 A) was measured, directly from the
generated atom coordinates, to be only about 3 atomic planes thick. This is deep in the
ultrathin-nanowire regime -- sub-1 nm FCC metal nanowires are well documented across EAM/DFT
literature to undergo substantial surface-driven reconstruction under plain relaxation, with
zero applied load. This matches exactly what was observed: the raw generated geometry was
confirmed clean (twice, independently, including a from-scratch re-parse on a different
machine), but the structure was visibly distorted at frame 0 of a real LAMMPS dump -- i.e.
*after* `minimize` and thermalization, *before* compression began.

**Mitigation applied:** strut diameter was increased 5 -> 12 A (roughly 6-7 atomic planes,
confirmed by direct cross-section counting), while re-checking that clearance between
non-adjacent struts still stayed a healthy multiple of the real EAM interaction cutoff, and
that the topology was still confirmed auxetic by the frame-FE solver at the new thickness.

**Why this is still open:** there is no LAMMPS installation available in the environment this
pipeline was developed in, so the 12 A fix is a reasoned engineering estimate from nanowire
literature trends, not a confirmed result. **Before trusting a full production run, verify
this yourself**: truncate `in.auxetic_compress` to stop right after `minimize` (comment out
thermalization and compression), dump that single frame, and confirm in OVITO that the
structure still resembles the input geometry. If it still visibly reconstructs, 12 A is not
thick enough and `strut_t` needs to go higher.

## Development Notes

This pipeline was built iteratively against real LAMMPS output and real OVITO renders, and
more than one fix attempt along the way was wrong before the real cause was found:

- **The first lattice topology was never actually auxetic.** The original strut network used
  two full diagonal struts crossing through a single shared center per cell ("X-in-a-box"),
  assumed by inspection to behave like the reference figure. A hand-written 2D Euler-Bernoulli
  beam finite-element solver -- first self-validated against the closed-form cantilever
  tip-deflection formula (`P L^3 / 3EI`, exact match, 0.00e+00 relative error) -- was then run
  on the exact same node/strut network with the same boundary conditions as the LAMMPS script.
  Result: effective Poisson's ratio +2.4. That topology is an ordinary scissor/cross-brace
  mechanism that bulges outward under compression, the opposite of auxetic.
- **The corrected topology's first version also had the wrong sign.** Switching to the actual
  named "bow-tie" re-entrant cell (horizontal waist + inclined struts, per Plewa et al. 2024)
  with struts flaring *outward* from the waist still gave +0.65. Re-reading the source paper's
  physical strut-"blocking" condition (struts colliding with each other under compression) only
  makes sense if the struts flare *inward*, back over the waist -- rebuilding it that way and
  re-running the same validated solver gave consistently negative results (-0.26 to -2.71
  depending on h/l and theta0), matching the paper's reported qualitative trends.
- **A stress output that always read `-0` for an entire run -- a real bug, not a physics
  result.** `fix setforce group 0 0 0` overwrites (zeroes) that group's force array every
  timestep; a separate `compute reduce sum fy` on the same group reads that *already-zeroed*
  array, guaranteeing exactly 0 forever. This was root-caused by reading `fix setforce`'s own
  documented output: it separately keeps a record of the total force *before* zeroing,
  accessible as `f_ID[1,2,3]` -- the correct way to read a rigid support's true reaction force.
- **A geometry that was numerically auxetic still looked like a plain grid in OVITO.** At
  theta0=15 deg, the lattice's lateral zigzag amplitude was only 14% of the waist length --
  real, but small enough that OVITO's true-scale atom rendering visually flattened it into
  what looked like a straight grid, even though a schematic thin-line preview of the same
  geometry looked like an obvious bow-tie. The preview generator itself was rewritten to draw
  true-scale filled circles instead of a schematic line, and theta0 was increased to 30 deg
  (28% throw) specifically so the pattern would be visible at real rendering scale, not just in
  an idealized line plot.
- **Atoms "sticking together" during compression traced to the wrong threshold.** A first
  clearance check only verified that non-adjacent struts didn't literally, geometrically
  overlap (gap > 0). The real threshold is the EAM potential's actual interaction range --
  confirmed directly from a real LAMMPS neighbor-list report in an actual run's log
  (`master list distance cutoff = 8.28721` minus the `neighbor 2.0` skin = 6.28721 A) -- not
  zero. Struts only 1.19x that cutoff apart were bonding across the gap during compression
  (confirmed independently by `PotEng` falling instead of rising during "compression," which
  is the signature of new, unintended cohesive bonds forming). The clearance check now tests
  against this real cutoff, with margin, not just against literal overlap.
- **A camera-tilt hypothesis was tested directly and ruled out.** When a later distortion
  turned out to be from the minimization stage (see Known Open Issue above), the first
  hypothesis was an OVITO camera not looking exactly top-down, since the lattice has real
  thickness in z. This was tested directly -- simulating an 8-degree-tilted orthographic
  projection of the actual generated atoms -- and found to be a negligible effect at this
  slab thickness (~10-12 A). That result is what correctly redirected the investigation toward
  the real cause (strut-thickness-driven reconstruction during `minimize`) instead of chasing
  a rendering setting that wasn't the problem.

## Design Decisions and Why

| # | Decision | Why |
|---|----------|-----|
| 1 | Re-entrant "bow-tie" topology (horizontal waist + inward-flaring inclined struts), per Plewa et al. 2024, not a generic honeycomb or a crossing-X guess | The only topology in this project actually confirmed, by a validated FE solver, to have a negative effective Poisson's ratio |
| 2 | Struts generated by masking an FCC-Al lattice fill with a distance-to-centerline ("capsule") test, rather than hand-placing atoms | Standard, robust, fully parametric method for turning an arbitrary 2D strut network into an atomistic mesh |
| 3 | Quasi-2D slab: 3 unit cells thick in z, periodic in that direction | Lets a genuinely 2D auxetic pattern be studied with a fully 3D EAM potential, without needing true free surfaces in the thin direction |
| 4 | `check_clearance()` tests against the real EAM interaction cutoff (6.28721 A, confirmed from an actual LAMMPS neighbor-list report), not just literal geometric overlap | Two struts don't need to touch to spuriously bond -- they just need to be within the potential's actual interaction range |
| 5 | Strut diameter raised to 12 A (from an initial 5 A) | A 5 A rod is only ~3 atomic planes across -- in the documented ultrathin-metal-nanowire reconstruction regime. See Known Open Issue: not yet empirically confirmed sufficient |
| 6 | theta0 = 30 deg, not a shallower angle | 15 deg was numerically auxetic but visually indistinguishable from a plain grid once rendered at true atomic scale; 30 deg keeps the zigzag amplitude large enough to actually look like a bow-tie |
| 7 | Bottom grip fully excluded from time integration (not just zero-velocity); top grip driven by `fix move linear` with lateral velocity pinned to 0 | Standard rigid, no-slip platen compression test; guarantees the bottom truly never moves regardless of applied force |
| 8 | Reaction force read from `fix setforce`'s own pre-zeroing force output (`f_ID[1,2,3]`), not a separate `compute reduce` on the same group | The group's force array is already zeroed by `setforce` every step -- any downstream compute reading it afterward is guaranteed to return zero |
| 9 | A hand-written 2D beam FE solver, self-validated against the closed-form cantilever formula, used to check topology and boundary conditions before any atomistic run | Verifying the checking tool itself before trusting it to judge the lattice's sign of Poisson's ratio |
| 10 | `strain_target` lowered to 0.06 (from an initial exploratory 0.25) | Large imposed strain deliberately drives these re-entrant struts toward literal self-contact ("blocking"), a real documented behavior of this cell type -- but shouldn't be the default a first run silently walks into |

## Simulation Overview

| Property | Value |
|----------|-------|
| Material | Al (Mishin et al. 1999 EAM, `Al99.eam.alloy`) |
| Lattice type | Re-entrant "bow-tie" auxetic cell (Plewa et al. 2024) |
| Unit cell | `l_strut`=39 A, `h_waist`=70 A (h/l=1.8), `theta0`=30 deg |
| Tessellation | 5 columns x 6 rows |
| Strut diameter | 12 A (~6-7 atomic planes; see Known Open Issue) |
| Slab thickness (z) | 3 unit cells (~12 A), periodic |
| FCC lattice constant | a0 = 4.05 A |
| Atom count | 56,442 |
| Simulation box | x ~ [-10.9, 359.5], y ~ [-44.7, 382.5], z = [0, 12.1] A |
| Boundary conditions | `s s p` -- shrink-wrapped x/y (finite lattice), periodic z (quasi-2D slab) |
| Timestep | 0.001 ps (1 fs) |
| Thermostat | NVT on mobile atoms only, 300 K, damping 0.1 ps |
| Grip thickness | 8 A at top and bottom |
| Compressive strain rate | 5.0e-4 /ps |
| Target engineering strain | 0.06 (6%) |
| Frame-FE effective Poisson's ratio | approx -0.31 to -0.85 depending on theta0/strut_t (see Development Notes) |

## System Geometry

```
y (loading direction)
^
|   ================================================   <- flat loaded top edge
|    \\      /\\      /\\      /\\      /\\      /
|     \\    /  \\    /  \\    /  \\    /  \\    /
|      \\  /    \\  /    \\  /    \\  /    \\  /
|       \\/      \\/      \\/      \\/      \\/         <- inward-flaring inclined
|       /\\      /\\      /\\      /\\      /\\            struts (length l, angle
|      /  \\    /  \\    /  \\    /  \\    /  \\           theta0 from vertical)
|     /    \\  /    \\  /    \\  /    \\  /    \\
|   [==][==][==][==][==][==][==][==][==][==]           <- horizontal waist struts
|     \\    /  \\    /  \\    /  \\    /  \\    /            (length h) -- one row of
|      \\  /    \\  /    \\  /    \\  /    \\  /             a bow-tie unit repeats
|       \\/      \\/      \\/      \\/      \\/              vertically 6x (Ny)
|       /\\      /\\      /\\      /\\      /\\
|      /  \\    /  \\    /  \\    /  \\    /  \\
|   ================================================   <- flat fixed bottom edge
+----------------------------------------------------> x (5 columns, Nx)

Bottom edge: fully fixed (excluded from time integration) -- rigid support.
Top edge: driven downward at constant velocity, zero lateral slip
          (fix move linear 0 vy 0) -- rigid, no-slip compression platen.
Interior nodes: free -- lateral narrowing of the waist openings here under
                axial compression is the actual visual signature of a
                negative Poisson's ratio.
```

## Simulation Phases

```
Python Geometry Generation  (generate_auxetic_lattice.py: build strut
  centerlines -> check_clearance() against EAM cutoff, hard-fails if
  unsafe -> fill FCC-Al lattice -> mask by distance-to-strut -> write
  lattice_auxetic.data)
      |
      v
LAMMPS Initialization  (units metal; atom_style atomic; read_data;
  pair_style eam/alloy; pair_coeff * * Al99.eam.alloy Al)
      |
      v
Grip Definition  (region + group by y-bounds: gbottom, gtop, gmobile;
  reference lengths H0/A0 frozen at parse time)
      |
      v
Energy Minimization  (bottom grip held via setforce during minimize
  only -- see Known Open Issue re: strut-thickness reconstruction here)
      |
      v
Thermalization  (fix nvt on gmobile ONLY, temp compute restricted to
  gmobile via fix_modify; both grips held; 5 ps @ 300 K)
      |
      v
Displacement-Controlled Compression  (gtop: fix move linear 0 vy 0,
  kinematic, no lateral slip; gbottom: fix setforce 0 0 0, excluded
  from integration -> stays rigid; gmobile: nve + nvt throughout)
      |
      v
Stress/Strain Logging  (v_stress from f_fsetbot[2], NOT a compute on
  the same already-zeroed group; v_strain from top-grip displacement;
  fix ave/time -> stress_strain.txt; dump every 2000 steps)
      |
      v
Final Structure  ->  lattice_auxetic_compressed.data + dump.auxetic_compress.lammpstrj
```

## Repository Structure

```
auxetic_bowtie/
|
├── generate_auxetic_lattice.py    # Geometry generator (Python): builds the bow-tie strut
|                                  #   network, runs check_clearance(), fills FCC-Al atoms,
|                                  #   writes the LAMMPS data file + a true-scale preview PNG
├── in.auxetic_compress             # LAMMPS input script: minimize -> thermalize -> compress
├── clearance_check.py              # Standalone segment-to-segment clearance checker
|                                  #   (also embedded directly in the generator)
├── diagnose_local.py               # Standalone script to re-parse and re-plot an existing
|                                  #   .data file directly, bypassing OVITO, for debugging
├── frame_fe_check_1_crossX_FAILS.py   # Self-validated beam-FE check of the WRONG (crossing-X)
|                                      #   topology -- kept as a documented negative result
├── frame_fe_check_2_bowtie_PASSES.py  # Same solver, correct (inward-flaring bow-tie) topology
├── plot_stress_strain.py           # Plots stress_strain.txt into a PNG
├── Al99.eam.alloy                  # EAM/alloy potential file (required, user-supplied --
|                                  #   NIST Interatomic Potentials Repository)
├── README.md                       # This file
|
└── output/                         # Generated on run
    ├── lattice_auxetic.data          # Atomistic geometry (LAMMPS data file, atom_style atomic)
    ├── lattice_preview.png           # True-atom-scale preview of the generated geometry
    ├── stress_strain.txt             # v_strain, v_stress vs. time (fix ave/time output)
    ├── dump.auxetic_compress.lammpstrj  # Full trajectory (every 2000 steps): id,type,x,y,z,fx,fy,fz
    └── lattice_auxetic_compressed.data  # Final relaxed/compressed structure (write_data)
```

## Requirements

- LAMMPS (with standard `eam/alloy` pair style support): https://www.lammps.org
- `Al99.eam.alloy` potential file -- NIST Interatomic Potentials Repository (search "Al 1999
  Mishin"), or LAMMPS's own `potentials/` directory if installed via conda/source
- Python 3 with `numpy`, `matplotlib`, `scipy`
- OVITO for visualization: https://www.ovito.org

## Installation

```bash
# LAMMPS via conda-forge (includes eam/alloy support)
conda install -c conda-forge lammps

# Python dependencies
pip install numpy matplotlib scipy
```

## Running the Simulation

```bash
# 1. Generate the geometry (edit l_strut/h_waist/theta0_deg/Nx/Ny/strut_t at the top first)
python3 generate_auxetic_lattice.py

# 2. Make sure Al99.eam.alloy sits next to in.auxetic_compress, then run LAMMPS
lmp -in in.auxetic_compress
# or, for speed:
mpirun -np 4 lmp -in in.auxetic_compress

# 3. Plot the stress-strain curve
python3 plot_stress_strain.py
```

Check the printed `Clearance check:` line and `Bottom grip atoms / Top grip atoms / Mobile
atoms` counts before committing to a long run -- if either grip count is 0, widen the `grip`
variable in `in.auxetic_compress`.

## Simulation Parameters

### Geometry (`generate_auxetic_lattice.py`)

| Variable | Meaning | Default |
|----------|---------|---------|
| `l_strut` | Inclined strut length | 39.0 A |
| `h_waist` | Horizontal waist strut length | 70.0 A |
| `theta0_deg` | Acute angle of inclined struts from vertical | 30.0 deg |
| `Nx`, `Ny` | Columns x rows of bow-tie cells | 5, 6 |
| `strut_t` | Strut diameter | 12.0 A |
| `n_layers_z` | FCC unit cells stacked in z (quasi-2D thickness) | 3 |
| `a0` | FCC Al lattice constant | 4.05 A |
| `eam_cutoff` | Real EAM interaction cutoff, confirmed from a LAMMPS neighbor-list report | 6.28721 A |

### Loading (`in.auxetic_compress`)

| Variable | Meaning | Default |
|----------|---------|---------|
| `grip` | Grip-region thickness at top/bottom | 8.0 A |
| `Tset` | Thermostat target temperature | 300.0 K |
| `strain_rate` | Prescribed compressive strain rate | 5.0e-4 /ps |
| `strain_target` | Target engineering strain to run to | 0.06 |
| timestep | Integration timestep | 0.001 ps |

### Contact Model

| Parameter | Meaning | Value |
|-----------|---------|-------|
| Pair style | `eam/alloy`, `Al99.eam.alloy` | -- |
| Interaction cutoff | Confirmed from a real neighbor-list report (`8.28721` minus `neighbor 2.0` skin) | 6.28721 A |
| Clearance requirement | `check_clearance()` hard-fails below this cutoff, warns below 1.5x it | >= 1.5x cutoff |

## Visualization in OVITO

1. **File -> Load File** -- open `lattice_auxetic.data` directly, or load
   `dump.auxetic_compress.lammpstrj` as a trajectory for the full compression run.
2. **Force a true top view.** "Ortho" only means orthographic projection, not that the camera
   is looking straight down an axis -- click the viewport's orientation label and explicitly
   pick "Top" before judging the geometry visually.
3. **Watch the waist openings**, not just the overall outline: they should visibly narrow as
   compression proceeds. That lateral narrowing under axial compression is the real visual
   signature of a negative Poisson's ratio -- an overall grid-like or barrel-bulging response
   would indicate something is wrong.
4. **Check frame 0 (right after minimize/thermalize, before compression) against the input
   geometry.** If it already looks reshaped or forked at that point, see Known Open Issue above
   -- that is very likely the unverified strut-thickness effect, not a compression artifact.
5. Cross-reference `stress_strain.txt` (`v_strain`, `v_stress`) against what's visually
   happening -- stress should be nonzero and grow in magnitude once compression starts (this
   required the `f_fsetbot[2]` fix; see Development Notes).

## What to Expect

**During minimization:** energy should converge smoothly. If the structure visibly reshapes
here, before any load is applied, that is the known open strut-thickness issue above, not a
script bug -- see the Known Open Issue section for how to confirm and what to try next.

**During thermalization:** temperature should settle near 300 K; it's normal for it to run a
few percent low, since the frozen/driven grips act as a persistent heat sink that a short 5 ps
equilibration doesn't fully compensate for.

**During compression:** `v_stress` should be nonzero and grow in magnitude (not `-0`); the
waist openings should visibly narrow, not widen. Avoid pushing `strain_target` far past ~10-15%
without re-running `check_clearance()`-style reasoning at the deformed state -- this cell type
is designed to approach literal strut self-contact ("blocking") at large strain, per the source
literature. That's a real mechanical limit of the design, not a numerical error, but it
shouldn't be the default a run walks into unintentionally.

## Common Errors and Fixes

| Error / Symptom | Cause | Fix |
|------------------|-------|-----|
| `v_stress` reads exactly `-0` for an entire run | `fix setforce group 0 0 0` zeroes that group's force array every step; a `compute reduce sum fy` on the same group reads the already-zeroed array | Use `f_fsetbot[2]` (fix setforce's own pre-zeroing force record) instead of a separate compute on the same group |
| Atoms visibly merge/"stick together" during compression; `PotEng` falls instead of rises | Non-adjacent struts closer than the real EAM interaction cutoff (6.29 A), not literally overlapping -- they don't need to touch to spuriously bond | Increase h/l or reduce `strut_t` until `check_clearance()` reports several multiples of the real cutoff, not just `> 0` |
| Structure looks like a plain grid in OVITO, not an obvious bow-tie, despite a confirmed negative Poisson's ratio | Lateral zigzag amplitude too small relative to strut thickness and real rendering scale to read visually, even though a schematic line preview looks fine | Increase `theta0`; verify with a true-atom-scale preview, not a thin schematic line plot |
| Structure visibly distorted/forked right after `minimize`, before any compression | Struts only a few atomic planes thick -- in the documented ultrathin-metal-nanowire reconstruction regime | Increase `strut_t` (currently 12 A; see Known Open Issue -- not yet empirically confirmed sufficient) |
| (Investigated, ruled out) Suspected OVITO camera tilt causing a doubled/forked projection of the thin z-slab | Tested directly via a simulated projection at up to 8 degrees tilt; found negligible at this slab thickness (~10-12 A) | Don't chase this -- check for minimize-stage reconstruction (previous row) instead |
| `WARNING: Temperature for fix modify is not for group all` | Expected and benign -- the thermostat's temperature compute is deliberately restricted to `gmobile` only, so the frozen/driven grips don't distort it | No action needed |

## Extending the Simulation

| Extension | What to Change |
|-----------|-----------------|
| Confirm the strut-thickness fix empirically | Run only through `minimize` (comment out thermalization/compression), dump that frame, and compare against the input geometry in OVITO -- see Known Open Issue |
| Tune the magnitude of the negative Poisson's ratio | Vary `theta0` and `h_waist`/`l_strut` -- smaller theta0 and smaller h/l both increase \|nu\|, per the frame-FE checks and the source literature's reported trend |
| Study strain-rate dependence | Vary `strain_rate` in `in.auxetic_compress` and compare resulting stress-strain curves |
| Remove grip artifacts | Replace the rigid, no-slip `fix move linear 0 vy 0` top grip with one that leaves lateral motion free, to isolate boundary effects from bulk auxetic response |
| Finite-size effects | Vary `Nx`/`Ny` and compare the effective Poisson's ratio and stress-strain response across lattice sizes |
| Different element/potential | Swap `Al99.eam.alloy` for another EAM/alloy potential, updating `pair_coeff`, `mass`, and `a0` in the generator to match |

## Citation

If you use this simulation pipeline in your research, please cite:

```bibtex
@software{mishra_auxetic_bowtie_md,
  author    = {Mishra, Akshansh},
  title     = {Molecular Dynamics Simulation of Compression of a Re-Entrant "Bow-Tie" Auxetic Al Lattice -- EAM Potential},
  year      = {2026},
  publisher = {Zenodo},
  doi       = {10.5281/zenodo.22122726},
  url       = {https://doi.org/10.5281/zenodo.22122726}
}
```

Plain text citation:

> Mishra, A. (2026). *Molecular Dynamics Simulation of Compression of a Re-Entrant "Bow-Tie"
> Auxetic Al Lattice — EAM Potential* [Computer software]. Zenodo.
> https://doi.org/10.5281/zenodo.22122726

## References

- Y. Mishin, D. Farkas, M.J. Mehl, D.A. Papaconstantopoulos, "Interatomic potentials for
  monoatomic metals from experimental data and *ab initio* calculations," *Phys. Rev. B* 59,
  3393 (1999). Source of the `Al99.eam.alloy` potential (NIST Interatomic Potentials
  Repository / LAMMPS `potentials/` directory).
- J. Plewa, M. Plonska, K. Feliksik, G. Junak, "Experimental Study of Auxetic Structures Made
  of Re-Entrant ('Bow-Tie') Cells," *Materials* 2024, 17(13), 3061.
  https://doi.org/10.3390/ma17133061. Source of the re-entrant bow-tie unit-cell geometry and
  its Poisson's-ratio trends used to verify this project's lattice topology.

## License

_Add your preferred license here (e.g. MIT, CC-BY-4.0, CC-BY-NC-4.0)._
