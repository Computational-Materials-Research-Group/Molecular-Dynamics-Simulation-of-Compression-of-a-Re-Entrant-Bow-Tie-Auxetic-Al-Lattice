"""
generate_auxetic_lattice.py
----------------------------------------------------------------------
Builds an atomistic (FCC-Al) model of a genuine re-entrant "bow-tie"
auxetic lattice and writes it out as a LAMMPS data file (atom_style
atomic), ready to be used with the Al99.eam.alloy EAM potential.

*** IMPORTANT CORRECTION vs. the first version of this script ***
The first version used a "crossing-X" cell (two full diagonals through
a shared center). A frame-element (Euler-Bernoulli beam) FE check of
that exact topology showed it has a POSITIVE effective Poisson's ratio
(+2.4) -- i.e. it is NOT auxetic, it's an ordinary scissor/cross-brace
mechanism that bulges outward under compression. That check is in
frame_fe_check.py (which first self-validates the beam element against
the textbook cantilever-deflection formula, exact match).

The geometry below instead follows the actual named "bow-tie"
re-entrant cell from:
  J. Plewa, M. Plonska, K. Feliksik, G. Junak, "Experimental Study of
  Auxetic Structures Made of Re-Entrant ('Bow-Tie') Cells," Materials
  2024, 17(13), 3061. https://doi.org/10.3390/ma17133061 (open access)
parameterized by a horizontal "waist" strut h, an inclined strut l,
and the acute angle theta between them (theta measured from vertical
here). The struts flare INWARD from each end of the waist (this is
what makes it re-entrant/auxetic, and it's also what makes physical
sense of the paper's strut-blocking condition h/(2l)=sin(theta_limit)
in compression -- inward struts can collide as theta grows; outward
ones never would):

           WL ---- h ---- WR
          /  \\           /  \\
         D    U         Up    Dp        (U/Up point up-and-INWARD,
                                          D/Dp point down-and-INWARD)

U of one cell coincides with D of the cell directly above it (same
column) -- that's what makes a self-supporting stacked column, matching
the picture's vertical zig-zag columns. Columns are placed side by
side and tied together only at the global top (loaded) and bottom
(fixed) edges, exactly like the flat top/bottom rows in the figure.

This was cross-checked, not assumed: frame_fe_check_v3.py reruns the
SAME validated beam-FE solver on this exact node/strut network and
confirms a NEGATIVE effective Poisson's ratio (e.g. nu=-0.80 for
h/l=1.5, theta0=20 deg; nu=-2.71 for h/l=0.8, theta0=10 deg -- the
trend, larger |nu| for smaller h/l and smaller theta0, matches the
paper's reported qualitative trend).

If your real design uses different h/l/theta or a different repeat
count, change the parameters below -- everything downstream (meshing,
atom fill, data file) is generic and works off the segment list.

ATOMISTIC FILL
----------------------------------------------------------------------
Struts are treated as "capsules" of finite thickness `t`: an FCC-Al
lattice (conventional cubic cell, a0 = 4.05 A, per Mishin et al. 1999 /
Al99.eam.alloy) is generated over the bounding box and every candidate
site is kept only if its perpendicular distance to the *nearest* strut
segment is <= t/2. The z-direction is kept thin (a few unit cells) and
periodic, i.e. a quasi-2D slab -- this is the standard way to run a
"2D" auxetic lattice with a 3D EAM potential without needing a genuine
free surface in z.

Reference for Al99.eam.alloy (verified, not guessed):
  Y. Mishin, D. Farkas, M.J. Mehl, D.A. Papaconstantopoulos,
  Phys. Rev. B 59, 3393 (1999); LAMMPS setfl conversion by C.A. Becker,
  NIST Interatomic Potentials Repository. FCC Al, a0 = 4.05 A,
  mass = 26.98 (release notes: Al99_releaseNotes_1.pdf).
  Required LAMMPS lines (from the file's own release notes):
      units metal
      atom_style atomic
      pair_style eam/alloy
      pair_coeff * * Al99.eam.alloy Al
----------------------------------------------------------------------
"""

import numpy as np

# ----------------------------------------------------------------------
# 1. USER PARAMETERS -- edit these to match your exact design
# ----------------------------------------------------------------------
# NOTE ON THIS CHOICE: an earlier version of this script used
# h=27, l=18, theta0=20, strut_t=7, which a proper segment-to-segment
# clearance check (see check_clearance() below) showed had only 4.57 A
# of surface-to-surface gap at ZERO strain -- just 0.65x the strut
# diameter. Since compression genuinely closes this gap further (that
# IS the auxetic mechanism), struts made contact almost immediately,
# which is what showed up as merged/"sticking" atoms in OVITO. The
# values below were re-picked to give ~1.5x the strut diameter of
# clearance at zero strain instead, and re-verified auxetic (nu_eff =
# -0.85, via the same frame-FE solver as frame_fe_check_2_bowtie_PASSES.py).
l_strut   = 39.0   # inclined strut length [Angstrom]
h_waist   = 70.0   # horizontal waist strut length [Angstrom]  (h/l = 1.8)
theta0_deg = 30.0  # acute angle from vertical [deg]
Nx = 5             # number of side-by-side columns
Ny = 6             # number of bow-tie cells stacked per column
strut_t = 12.0     # strut thickness [Angstrom]. INCREASED from 5.0: a 5 A rod is only
                   # ~3 atomic planes across (measured directly from the generated file),
                   # deep in the ultrathin-nanowire regime where real EAM relaxation causes
                   # substantial surface-driven reconstruction even with zero applied load
                   # -- this is a documented effect for sub-1nm FCC metal nanowires, and
                   # matches what showed up as "distortion" in OVITO immediately after
                   # `minimize` (frame 0, before any compression). 12 A (~7-8 atomic planes)
                   # is a more reasonable compromise: still 3.0x the EAM cutoff in clearance,
                   # still clearly auxetic (nu_eff ~ -0.31 vs -0.40 at strut_t=5).
                   # CAVEAT: I don't have LAMMPS available to directly confirm this resolves
                   # the reconstruction -- verify empirically by running just `minimize`
                   # (comment out the thermalization/compression sections) and checking that
                   # frame 0 still resembles the input geometry before committing to a full run.
# Checked together this time, not one at a time (h=70,l=39,theta0=30,strut_t=5):
#   - nu_eff = -0.41 (frame_fe_check_2_bowtie_PASSES.py)               -> genuinely auxetic
#   - clearance = 26.0 A = 4.14x the real EAM cutoff (6.29 A)          -> safe from spurious bonding
#   - theta_limit (strut self-blocking) = 63.8 deg, 33.8 deg of margin -> safe at large strain
#   - lateral zigzag throw dx = l*sin(theta0) = 19.5 A = 28% of h      -> actually visible as a
#     bow-tie at true atomistic scale, not just in a schematic line plot (see note below).
#
# NOTE ON WHY theta0=15 (previous version) LOOKED LIKE A PLAIN GRID IN OVITO:
# dx there was only 10.1 A, 14% of h -- a real but small wobble. A thin schematic
# line plot (like the preview below) still reads as "hourglass" to the eye at any
# amplitude, but OVITO's default sphere rendering draws atoms noticeably larger
# than their real radius, which visually flattens a small wobble into what looks
# like a straight bar. theta0=30 makes the wobble ~2x bigger for exactly this
# reason. Separately: if it still looks flatter than expected in OVITO, reduce
# the particle radius scaling in OVITO's Particle display panel -- that setting
# is cosmetic (display-only) and can make an otherwise-correct geometry look
# like a plain grid if left at OVITO's default.
eam_cutoff = 6.28721  # Angstrom. CONFIRMED from an actual LAMMPS run's own neighbor-list
                       # report: "master list distance cutoff = 8.28721" with `neighbor 2.0`
                       # in in.auxetic_compress -> 8.28721-2.0 = 6.28721, the real eam/alloy
                       # pairwise cutoff for Al99.eam.alloy. This -- not zero -- is the
                       # relevant threshold: two struts don't need to geometrically overlap
                       # to spuriously bond, they just need to come within this distance of
                       # each other, which is exactly what "atoms sticking together" during
                       # compression turned out to be (confirmed by PotEng falling instead
                       # of rising during compression in a real run -- new, unintended
                       # cohesive bonds forming between struts that were meant to stay apart).

n_layers_z = 3        # number of conventional FCC cubic cells stacked in z (quasi-2D thickness)
a0 = 4.05              # FCC Al lattice constant [Angstrom] (Mishin 1999 / Al99.eam.alloy, verified)
Al_mass = 26.981538     # [amu]

out_data_file = "lattice_auxetic.data"
out_png       = "lattice_preview.png"

# ----------------------------------------------------------------------
# 2. BUILD THE STRUT CENTERLINE SEGMENTS (verified re-entrant bow-tie)
# ----------------------------------------------------------------------
def build_segments(h, l, theta_deg, Nx, Ny, col_pitch=None):
    """Inward-flaring bow-tie topology, verified auxetic via
    frame_fe_check_v3.py. theta_deg measured from vertical."""
    th = np.radians(theta_deg)
    dx = l * np.sin(th)
    dy = l * np.cos(th)
    y_pitch = 2 * dy
    if col_pitch is None:
        col_pitch = h

    segs = []
    bottom_dangling, top_dangling = [], []

    for i in range(Nx):
        xcol = i * col_pitch
        for j in range(Ny):
            y0 = j * y_pitch
            WL = (xcol,       y0)
            WR = (xcol + h,   y0)
            U  = (xcol + dx,     y0 + dy)   # up,   INWARD from WL
            D  = (xcol + dx,     y0 - dy)   # down, INWARD from WL
            Up = (xcol + h - dx, y0 + dy)   # up,   INWARD from WR
            Dp = (xcol + h - dx, y0 - dy)   # down, INWARD from WR

            segs += [(WL, WR), (WL, U), (WL, D), (WR, Up), (WR, Dp)]

            if j == 0:
                bottom_dangling += [D, Dp]
            if j == Ny - 1:
                top_dangling += [U, Up]

    # tie dangling ends into flat support (bottom) / loaded (top) edges,
    # same role as the flat top/bottom rows in the reference figure
    bot = sorted(set(bottom_dangling), key=lambda p: p[0])
    top = sorted(set(top_dangling), key=lambda p: p[0])
    for p1, p2 in zip(bot[:-1], bot[1:]):
        segs.append((p1, p2))
    for p1, p2 in zip(top[:-1], top[1:]):
        segs.append((p1, p2))
    return np.array(segs)

segments = build_segments(h_waist, l_strut, theta0_deg, Nx, Ny)
print(f"Built {len(segments)} strut segments for a {Nx}x{Ny} re-entrant bow-tie lattice "
      f"(h={h_waist} A, l={l_strut} A, theta0={theta0_deg} deg).")


# ----------------------------------------------------------------------
# 2b. MANDATORY CLEARANCE CHECK -- catches struts that will overlap
# ----------------------------------------------------------------------
# Standard closest-point-between-two-segments algorithm (Ericson,
# "Real-Time Collision Detection", ClosestPtSegmentSegment), 2D (z=0).
def _closest_dist_seg_seg(p1, p2, p3, p4, eps=1e-12):
    p1, p2, p3, p4 = (np.array(p, dtype=float) for p in (p1, p2, p3, p4))
    d1, d2, r = p2 - p1, p4 - p3, p1 - p3
    a, e, f = d1 @ d1, d2 @ d2, d2 @ r
    if a <= eps and e <= eps:
        return np.linalg.norm(p1 - p3)
    if a <= eps:
        s, t = 0.0, np.clip(f / e, 0, 1)
    else:
        c = d1 @ r
        if e <= eps:
            t, s = 0.0, np.clip(-c / a, 0, 1)
        else:
            b = d1 @ d2
            denom = a * e - b * b
            s = np.clip((b * f - c * e) / denom, 0, 1) if abs(denom) > eps else 0.0
            t = (b * s + f) / e
            if t < 0:
                t, s = 0.0, np.clip(-c / a, 0, 1)
            elif t > 1:
                t, s = 1.0, np.clip((b - c) / a, 0, 1)
    c1, c2 = p1 + d1 * s, p3 + d2 * t
    return np.linalg.norm(c1 - c2)


def check_clearance(segs, strut_t, eam_cutoff, tol=1e-6, min_ratio=1.5):
    """Minimum surface-to-surface gap between struts that do NOT share
    an endpoint node (struts that DO share a node are expected to
    touch there -- that's a real joint, not an error).

    Checks against TWO thresholds:
      1. gap > 0            -- struts don't literally overlap
      2. gap > eam_cutoff   -- struts stay outside each other's actual
                               interaction range, so they don't form
                               spurious new bonds under compression
                               (this is the threshold that actually
                               matters for "atoms sticking together")
    """
    def node_key(p):
        return (round(p[0] / tol), round(p[1] / tol))
    n = len(segs)
    best = (np.inf, -1, -1)
    for i in range(n):
        p1, p2 = segs[i]
        keys_i = {node_key(p1), node_key(p2)}
        for j in range(i + 1, n):
            p3, p4 = segs[j]
            if keys_i & {node_key(p3), node_key(p4)}:
                continue
            d = _closest_dist_seg_seg(p1, p2, p3, p4)
            if d < best[0]:
                best = (d, i, j)
    gap = best[0] - strut_t
    print(f"Clearance check: nearest non-adjacent struts are {best[0]:.2f} A apart "
          f"(centerline-to-centerline); strut diameter is {strut_t:.2f} A "
          f"-> surface clearance = {gap:+.2f} A.")
    if gap <= 0:
        raise RuntimeError(
            "Struts overlap at zero strain (clearance <= 0). Reduce strut_t, "
            "reduce theta0, or increase h/l before generating atoms.")
    print(f"  vs. EAM interaction cutoff ({eam_cutoff} A): "
          f"margin = {gap-eam_cutoff:+.2f} A ({gap/eam_cutoff:.2f}x cutoff).")
    if gap <= eam_cutoff:
        raise RuntimeError(
            f"Clearance ({gap:.2f} A) is inside the EAM cutoff ({eam_cutoff} A) "
            f"even at zero strain -- struts WILL spuriously bond as soon as the "
            f"potential is evaluated, regardless of compression. Increase h/l or "
            f"reduce strut_t.")
    if gap / eam_cutoff < min_ratio:
        print(f"WARNING: clearance is only {gap/eam_cutoff:.2f}x the EAM cutoff. "
              f"Compression closes this further (that's the auxetic mechanism), "
              f"so struts may enter interaction range well before your target "
              f"strain. Consider more margin, or keep strain_target small in "
              f"in.auxetic_compress and check intermediate frames in OVITO.")
    return gap

check_clearance(segments, strut_t, eam_cutoff, min_ratio=1.5)

# ----------------------------------------------------------------------
# 3. POINT-TO-SEGMENT DISTANCE (vectorized) -> "IS THIS SITE INSIDE A STRUT?"
# ----------------------------------------------------------------------
def min_dist_to_segments(points_xy, segs):
    """points_xy: (N,2). segs: (M,2,2). Returns (N,) min distance to any segment."""
    p1 = segs[:, 0, :]                     # (M,2)
    p2 = segs[:, 1, :]                     # (M,2)
    d  = p2 - p1                           # (M,2)
    d2 = np.einsum('ij,ij->i', d, d)       # (M,)
    d2[d2 == 0] = 1e-12

    # broadcast points (N,1,2) against segments (1,M,2)
    w = points_xy[:, None, :] - p1[None, :, :]           # (N,M,2)
    t = np.einsum('nmk,mk->nm', w, d) / d2[None, :]      # (N,M)
    t = np.clip(t, 0.0, 1.0)
    proj = p1[None, :, :] + t[:, :, None] * d[None, :, :]  # (N,M,2)
    dist = np.linalg.norm(points_xy[:, None, :] - proj, axis=2)  # (N,M)
    return dist.min(axis=1)

# ----------------------------------------------------------------------
# 4. GENERATE CANDIDATE FCC SITES OVER THE BOUNDING BOX
# ----------------------------------------------------------------------
pad = strut_t
xmin, xmax = segments[:, :, 0].min() - pad, segments[:, :, 0].max() + pad
ymin, ymax = segments[:, :, 1].min() - pad, segments[:, :, 1].max() + pad
zlen = n_layers_z * a0

nx_cells = int(np.ceil((xmax - xmin) / a0)) + 1
ny_cells = int(np.ceil((ymax - ymin) / a0)) + 1
nz_cells = n_layers_z

# FCC basis (fractional coords of conventional cubic cell)
fcc_basis = np.array([[0.0, 0.0, 0.0],
                      [0.5, 0.5, 0.0],
                      [0.5, 0.0, 0.5],
                      [0.0, 0.5, 0.5]])

ix, iy, iz = np.meshgrid(np.arange(nx_cells), np.arange(ny_cells),
                          np.arange(nz_cells), indexing='ij')
cell_origin = np.stack([ix.ravel(), iy.ravel(), iz.ravel()], axis=1).astype(float)

pts = (cell_origin[:, None, :] + fcc_basis[None, :, :]).reshape(-1, 3) * a0
pts[:, 0] += xmin
pts[:, 1] += ymin
# keep only sites within the padded bounding box and one periodic z-length
mask_box = (pts[:, 0] >= xmin) & (pts[:, 0] <= xmax) & \
           (pts[:, 1] >= ymin) & (pts[:, 1] <= ymax) & \
           (pts[:, 2] >= 0.0)  & (pts[:, 2] < zlen)
pts = pts[mask_box]
print(f"Candidate FCC sites in bounding box: {len(pts)}")

# ----------------------------------------------------------------------
# 5. KEEP ONLY SITES INSIDE A STRUT ("CAPSULE" TEST)
# ----------------------------------------------------------------------
dist = min_dist_to_segments(pts[:, :2], segments)
keep = dist <= (strut_t / 2.0)
atoms = pts[keep]
print(f"Atoms retained inside struts (thickness={strut_t} A): {len(atoms)}")

if len(atoms) == 0:
    raise RuntimeError("No atoms generated -- check strut_t / a / b parameters.")

# ----------------------------------------------------------------------
# 6. WRITE LAMMPS DATA FILE (atom_style atomic)
# ----------------------------------------------------------------------
xlo, xhi = atoms[:, 0].min() - 5.0, atoms[:, 0].max() + 5.0
ylo, yhi = atoms[:, 1].min() - 5.0, atoms[:, 1].max() + 5.0
zlo, zhi = 0.0, zlen   # keep exact periodic length in z

with open(out_data_file, "w") as f:
    f.write("LAMMPS data file: re-entrant bow-tie auxetic lattice, FCC Al filler\n\n")
    f.write(f"{len(atoms)} atoms\n")
    f.write("1 atom types\n\n")
    f.write(f"{xlo:.6f} {xhi:.6f} xlo xhi\n")
    f.write(f"{ylo:.6f} {yhi:.6f} ylo yhi\n")
    f.write(f"{zlo:.6f} {zhi:.6f} zlo zhi\n\n")
    f.write("Masses\n\n")
    f.write(f"1 {Al_mass}\n\n")
    f.write("Atoms # atomic\n\n")
    for idx, (x, y, z) in enumerate(atoms, start=1):
        f.write(f"{idx} 1 {x:.6f} {y:.6f} {z:.6f}\n")

print(f"Wrote LAMMPS data file: {out_data_file}")
print(f"Box: x[{xlo:.1f},{xhi:.1f}]  y[{ylo:.1f},{yhi:.1f}]  z[{zlo:.1f},{zhi:.1f}]  (periodic in z)")

# ----------------------------------------------------------------------
# 7. QUICK VISUAL SANITY CHECK -- rendered at TRUE relative atom scale
# ----------------------------------------------------------------------
# IMPORTANT: a schematic thin-line/small-dot plot can visually read as an
# obvious "hourglass" even when the actual zigzag amplitude is small,
# because your eye traces the centerline regardless of scale. That's
# exactly what went wrong before (theta0=15 case): the schematic preview
# looked fine, but real atoms rendered at their true radius in OVITO
# flattened the small wobble into what looked like a plain grid. This
# version draws atoms as filled circles at their actual radius (in data
# units, matched to strut_t) so the preview honestly predicts what an
# atomistic viewer will show.
try:
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    from matplotlib.collections import EllipseCollection

    fig, ax = plt.subplots(figsize=(8, 10))
    atom_radius = a0 / 2.0 * 0.9   # true-ish visual atom radius, Angstrom
    ec = EllipseCollection(
        widths=2 * atom_radius, heights=2 * atom_radius, angles=0, units="xy",
        offsets=atoms[:, :2], offset_transform=ax.transData,
        facecolors=plt.cm.viridis((atoms[:, 2] - atoms[:, 2].min()) /
                                   max(np.ptp(atoms[:, 2]), 1e-9)),
        edgecolors="none",
    )
    ax.add_collection(ec)
    ax.set_xlim(atoms[:, 0].min() - 5, atoms[:, 0].max() + 5)
    ax.set_ylim(atoms[:, 1].min() - 5, atoms[:, 1].max() + 5)
    ax.set_aspect("equal")
    ax.set_title(f"Auxetic bow-tie lattice preview, TRUE atom scale ({len(atoms)} atoms)")
    ax.set_xlabel("x [A]"); ax.set_ylabel("y [A]")
    fig.tight_layout()
    fig.savefig(out_png, dpi=150)
    print(f"Wrote preview image (true atom scale): {out_png}")
except Exception as e:
    print("Preview plot skipped:", e)
