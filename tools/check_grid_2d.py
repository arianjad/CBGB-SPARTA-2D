#!/usr/bin/env python3
"""Check a 2D SPARTA grid for face-adjacent cells more than 2:1 apart in size.

The hierarchical grid is a 2x refinement tree, so every cell size is
base/2^(L-1).  A "smooth" grid is one where no two cells that share a face
differ by more than one level (a 2:1 face jump).  A 4:1 jump (two levels)
is what `adapt_grid ... region <band> one` produces when a fine band is
dropped straight onto the base grid with no buffer.

Input: any SPARTA 2D grid dump with columns `id xc yc xlo ylo xhi yhi ...`
(the `field.grid` of a production run, or a `dump grid` from a run-0 deck).
Only the first snapshot is read -- the geometry does not change.

Method: rasterise every cell onto the finest uniform lattice present, then
walk the lattice in +x and +y.  Where two lattice squares belong to
different cells, that lattice edge IS a shared face; compare levels.

Usage:
  python tools/check_grid_2d.py <grid-dump> [--png out.png] [--label TEXT]
Exit status 0 = smooth (gate PASS), 1 = violations found.
"""
import argparse
import sys

import numpy as np


def read_cells(path):
    """Return (xlo, ylo, xhi, yhi) arrays and the box, from the first snapshot."""
    box = []
    rows = []
    ncells = None
    with open(path) as f:
        mode = None
        for line in f:
            if line.startswith("ITEM:"):
                if mode == "cells":
                    break                      # second snapshot -> stop
                if line.startswith("ITEM: NUMBER OF CELLS"):
                    mode = "n"
                elif line.startswith("ITEM: BOX BOUNDS"):
                    mode = "box"
                elif line.startswith("ITEM: CELLS"):
                    cols = line.split()[2:]
                    need = ["xlo", "ylo", "xhi", "yhi"]
                    missing = [c for c in need if c not in cols]
                    if missing:
                        sys.exit("dump lacks column(s) %s; need id xc yc xlo ylo xhi yhi"
                                 % ",".join(missing))
                    icol = [cols.index(c) for c in need]
                    mode = "cells"
                else:
                    mode = None
                continue
            if mode == "n":
                ncells = int(line)
            elif mode == "box":
                box.append([float(v) for v in line.split()])
            elif mode == "cells":
                p = line.split()
                rows.append([float(p[i]) for i in icol])
    if not rows:
        sys.exit("no cells parsed from %s" % path)
    if ncells is not None and len(rows) != ncells:
        sys.exit("parsed %d cells, header says %d" % (len(rows), ncells))
    a = np.array(rows, dtype=float)
    return a[:, 0], a[:, 1], a[:, 2], a[:, 3], (box[0], box[1])


def rasterise(xlo, ylo, xhi, yhi, box):
    """Map every cell onto the finest uniform lattice.  Returns (idx, level, nx, ny)."""
    # SPARTA dumps coordinates at %g (6 sig figs), so recover the lattice from the
    # box and the finest cell rather than trusting the printed differences.
    (bxlo, bxhi), (bylo, byhi) = box
    dx, dy = xhi - xlo, yhi - ylo
    if not np.allclose(np.log2(dx.max() / dx), np.rint(np.log2(dx.max() / dx)), atol=1e-3):
        sys.exit("cell sizes are not a power-of-2 hierarchy; this checker assumes 2x refinement")
    # refinement here is 2 2 1, so the radial level equals the axial level
    level = np.rint(np.log2(dx.max() / dx)).astype(np.int32) + 1
    nx = int(round((bxhi - bxlo) / dx.min()))
    ny = int(round((byhi - bylo) / dy.min()))
    hx, hy = (bxhi - bxlo) / nx, (byhi - bylo) / ny
    idx = np.full((ny, nx), -1, dtype=np.int32)
    i0 = np.rint((xlo - bxlo) / hx).astype(int)
    i1 = np.rint((xhi - bxlo) / hx).astype(int)
    j0 = np.rint((ylo - bylo) / hy).astype(int)
    j1 = np.rint((yhi - bylo) / hy).astype(int)
    for k in range(len(xlo)):
        idx[j0[k]:j1[k], i0[k]:i1[k]] = k
    if (idx < 0).any():
        sys.exit("%d lattice squares uncovered -- dump is not a complete grid" % int((idx < 0).sum()))
    return idx, level, nx, ny, hx, hy, bxlo, bylo


def violations(idx, level, hx, hy, bxlo, bylo):
    """Face-adjacent pairs with |level difference| >= 2.  Returns (pairs, faces)."""
    out_pairs, out_xy = [], []
    for axis in (0, 1):                      # 0 = +y neighbours, 1 = +x neighbours
        a = idx[:-1, :] if axis == 0 else idx[:, :-1]
        b = idx[1:, :] if axis == 0 else idx[:, 1:]
        m = (a != b) & (np.abs(level[a] - level[b]) >= 2)
        if not m.any():
            continue
        jj, ii = np.nonzero(m)
        lo, hi = np.minimum(a[m], b[m]), np.maximum(a[m], b[m])
        out_pairs.append(np.stack([lo, hi], axis=1))
        # face midpoint on the lattice edge between the two squares
        fx = bxlo + (ii + (1.0 if axis == 1 else 0.5)) * hx
        fy = bylo + (jj + (1.0 if axis == 0 else 0.5)) * hy
        out_xy.append(np.stack([fx, fy], axis=1))
    if not out_pairs:
        return np.zeros((0, 2), int), np.zeros((0, 2), float)
    pairs = np.concatenate(out_pairs)
    xy = np.concatenate(out_xy)
    # one entry per distinct cell pair, keeping the first face location
    _, first = np.unique(pairs, axis=0, return_index=True)
    return pairs[np.sort(first)], xy[np.sort(first)]


ZOOMS = [("inlet  x -5..3 mm", (-0.005, 0.003, 0.0, 0.004)),
         ("throat x 60..68 mm", (0.060, 0.068, 0.0, 0.008))]


def figure(idx, level, hx, hy, bxlo, bylo, nx, ny, faces, path, label):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    img = level[idx]
    ext = (bxlo * 1e3, (bxlo + nx * hx) * 1e3, bylo * 1e3, (bylo + ny * hy) * 1e3)
    fig, ax = plt.subplots(3, 1, figsize=(11, 10))
    for a, (title, lim) in [(ax[0], ("full domain", None)), (ax[1], ZOOMS[0]), (ax[2], ZOOMS[1])]:
        im = a.imshow(img, origin="lower", extent=ext, aspect="auto",
                      cmap="viridis", vmin=1, vmax=level.max(), interpolation="nearest")
        if faces.size:
            a.plot(faces[:, 0] * 1e3, faces[:, 1] * 1e3, "r.", ms=2.5,
                   label="%d faces > 2:1" % len(faces))
            a.legend(loc="upper right", fontsize=8)
        if lim is not None:
            a.set_xlim(lim[0] * 1e3, lim[1] * 1e3)
            a.set_ylim(lim[2] * 1e3, lim[3] * 1e3)
        a.set_title(title, fontsize=9)
        a.set_xlabel("x [mm]")
        a.set_ylabel("r [mm]")
        fig.colorbar(im, ax=a, label="level")
    fig.suptitle("%s -- grid level map (L1 = %.3g x %.3g um)"
                 % (label, hx * 2 ** (level.max() - 1) * 1e6, hy * 2 ** (level.max() - 1) * 1e6))
    fig.tight_layout()
    fig.savefig(path, dpi=130)
    print("wrote", path)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("dump")
    ap.add_argument("--png")
    ap.add_argument("--label", default="")
    args = ap.parse_args()

    xlo, ylo, xhi, yhi, box = read_cells(args.dump)
    idx, level, nx, ny, hx, hy, bxlo, bylo = rasterise(xlo, ylo, xhi, yhi, box)
    pairs, faces = violations(idx, level, hx, hy, bxlo, bylo)

    print("cells            %d" % len(xlo))
    for L in range(1, level.max() + 1):
        n = int((level == L).sum())
        if n:
            print("  L%-3d %8d cells   %.4g x %.4g um"
                  % (L, n, hx * 2 ** (level.max() - L) * 1e6, hy * 2 ** (level.max() - L) * 1e6))
    print("face-adjacent pairs above 2:1   %d" % len(pairs))
    if len(pairs):
        dlev = np.abs(level[pairs[:, 0]] - level[pairs[:, 1]])
        print("worst ratio      %d:1  (level jump %d)" % (2 ** dlev.max(), dlev.max()))
        la, lb = level[pairs[:, 0]], level[pairs[:, 1]]
        lo, hi = np.minimum(la, lb), np.maximum(la, lb)
        print("seams by level pair (face-location bounding box):")
        for key in sorted(set(zip(lo.tolist(), hi.tolist()))):
            m = (lo == key[0]) & (hi == key[1])
            f = faces[m]
            print("  L%d | L%d  %d:1  x %.6f..%.6f  r %.6f..%.6f  (%d faces)"
                  % (key[0], key[1], 2 ** (key[1] - key[0]),
                     f[:, 0].min(), f[:, 0].max(), f[:, 1].min(), f[:, 1].max(), m.sum()))
    print("GATE %s" % ("PASS -- every face-sharing pair is within 2:1" if not len(pairs)
                       else "FAIL"))
    if args.png:
        figure(idx, level, hx, hy, bxlo, bylo, nx, ny, faces, args.png,
               args.label or args.dump)
    return 0 if not len(pairs) else 1


if __name__ == "__main__":
    sys.exit(main())
