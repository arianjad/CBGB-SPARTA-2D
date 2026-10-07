#!/usr/bin/env python3
"""2D neon-flow plot in the plane through the axis (y = 0) of a cyl-3d frozen field.

The cells whose y-extent contains y = 0 are rasterised onto a 0.2 mm image (x, z).
Shown: isolines of the gas speed |u| (log-spaced levels, so both the ~0.1 m/s chamber
flow and the ~200 m/s plume are visible, labelled in m/s), the speed as a light filled
background, and thin streamlines of (u_x, u_z) for the direction. Cells with no gas
(walls) are blank. A light 3x3 box smoothing tames DSMC noise before contouring.

Usage: python tools/plot_flow_isolines.py frozen.grid --out OUT.png [--title ...]
Light (one frozen field, matplotlib).
"""
import argparse

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib.colors import LogNorm

PIX = 2e-4
NOISE = 10.0      # m/s: below this the coarse, short-averaged field is DSMC noise (chamber flow ~0.1 m/s)
# cyl-3d body (mm): inlet wall x -1.27..0 with a 1.27 mm bore, bore r 15.24, chamber to x 45.47,
# outer radius 19.05, far end to x 50.80 (solid when the outlet is plugged)
def solid_mask(X, Z, plugged=True):
    r = np.abs(Z)
    inlet = (X >= -1.27) & (X <= 0) & (r >= 1.27) & (r <= 19.05)
    side = (X > 0) & (X < 45.47) & (r >= 15.24) & (r <= 19.05)
    end = (X >= 45.47) & (X <= 50.80) & (r <= 19.05) & ((r >= 1.5875) | plugged)
    return inlet | side | end
INK, MUTED = "#1f1f1e", "#6b6a66"


def raster(path):
    d = pd.read_csv(path, sep=r"\s+", skiprows=9, header=None,
                    names=["id", "x", "y", "z", "xlo", "ylo", "zlo", "xhi", "yhi", "zhi", "n", "u", "v", "w", "T"])
    d = d[(d.ylo <= 0) & (d.yhi > 0)]
    x0, x1, z0, z1 = d.xlo.min(), d.xhi.max(), d.zlo.min(), d.zhi.max()
    nx, nz = int(round((x1 - x0) / PIX)), int(round((z1 - z0) / PIX))
    U = np.full((nz, nx), np.nan); W = np.full((nz, nx), np.nan); N = np.zeros((nz, nx))
    for r in d.itertuples():
        i0, i1 = int(round((r.xlo - x0) / PIX)), int(round((r.xhi - x0) / PIX))
        k0, k1 = int(round((r.zlo - z0) / PIX)), int(round((r.zhi - z0) / PIX))
        N[k0:k1, i0:i1] = r.n
        if r.n > 0:
            U[k0:k1, i0:i1] = r.u; W[k0:k1, i0:i1] = r.w
    xs = x0 + (np.arange(nx) + 0.5) * PIX; zs = z0 + (np.arange(nz) + 0.5) * PIX
    return xs, zs, U, W, N


def smooth(a, k=3):
    """k x k box mean ignoring NaNs."""
    m = np.isfinite(a); b = np.where(m, a, 0.0)
    pad = k // 2
    bp = np.pad(b, pad); mp = np.pad(m.astype(float), pad)
    s = sum(bp[i:i + a.shape[0], j:j + a.shape[1]] for i in range(k) for j in range(k))
    c = sum(mp[i:i + a.shape[0], j:j + a.shape[1]] for i in range(k) for j in range(k))
    out = np.where(c > 0, s / np.maximum(c, 1), np.nan)
    return np.where(m, out, np.nan)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("field")
    ap.add_argument("--out", required=True)
    ap.add_argument("--title", default="")
    ap.add_argument("--xlim", type=float, nargs=2, default=(-22.4, 50.8), help="mm")
    ap.add_argument("--open-outlet", dest="plugged", action="store_false", help="geometry has an open outlet")
    a = ap.parse_args()
    xs, zs, U, W, N = raster(a.field)
    Us, Ws = smooth(U), smooth(W)
    S = np.hypot(Us, Ws)
    X, Z = np.meshgrid(xs * 1e3, zs * 1e3)
    fig, ax = plt.subplots(figsize=(13, 6))
    SOL = solid_mask(X, Z, plugged=a.plugged)
    S = np.where(SOL, np.nan, S)                                      # empty vacuum cells stay blank (no data)
    Us = np.where(SOL, np.nan, np.nan_to_num(Us)); Ws = np.where(SOL, np.nan, np.nan_to_num(Ws))
    lev = [10, 20, 50, 100, 150, 200]
    smax = np.nanmax(S)
    lev = [v for v in lev if v < smax]
    bg = ax.pcolormesh(X, Z, np.ma.masked_invalid(np.clip(S, 1, None)), cmap="Blues",
                       norm=LogNorm(vmin=1, vmax=max(smax, lev[-1])), shading="auto", alpha=0.45)
    cs = ax.contour(X, Z, np.ma.masked_invalid(S), levels=lev, colors="#184f95", linewidths=1.0)
    ax.clabel(cs, fmt=lambda v: f"{v:g}", fontsize=7, inline=True)
    fast = np.isfinite(S) & (S >= NOISE)
    sp = ax.streamplot(xs * 1e3, zs * 1e3, np.ma.masked_where(~fast, Us), np.ma.masked_where(~fast, Ws),
                       density=1.6, color=MUTED, linewidth=0.5, arrowsize=0.6)
    sp.lines.set_alpha(0.6)
    ax.text(22.7, 0, "chamber flow ~0.1 m/s:\nbelow the noise of this\nshort, coarse run", ha="center", va="center",
            fontsize=9, color=MUTED)
    solid = np.ma.masked_where(~SOL, np.ones_like(N))
    ax.pcolormesh(X, Z, solid, cmap=matplotlib.colors.ListedColormap(["#d9d8d3"]), shading="auto")
    cb = fig.colorbar(bg, ax=ax, pad=0.01, shrink=0.8)
    cb.set_label("neon speed |u| (m/s)")
    ax.set_xlim(*a.xlim); ax.set_aspect("equal")
    ax.set_xlabel("x (mm; inlet wall at 0, ions and gas leave towards −x)", color=INK)
    ax.set_ylabel("z (mm)", color=INK)
    ax.set_title(a.title or "Neon flow in the plane through the axis", color=INK, fontsize=11)
    ax.text(0.01, 0.01, "blue lines: speed isolines (m/s), from 10 m/s, above the noise of this run; grey: streamlines; grey blocks: walls",
            transform=ax.transAxes, fontsize=8, color=MUTED)
    fig.tight_layout(); fig.savefig(a.out, dpi=130)
    print("->", a.out)


if __name__ == "__main__":
    main()
