#!/usr/bin/env python3
"""Tiny standalone Monte Carlo of the tracer's collision physics (isotropic CM scattering off
Ne at rest, table sigma) along the on-axis density n(x), with two ways to pick each free path:
  start : distance drawn with the density at the flight's starting point and flown in full
          (what tracer/ParticleTracing3D.jl does: freePath(..., props[7]) at the start point)
  exact : distance drawn so the integral of n sigma ds is an Exp(1) variate (correct for n varying
          along the flight)
Gas thermal motion is neglected (fine above ~1 eV). Prints the mean x where ions first drop below
1 eV, collisions so far, and where the dE/dx curve peaks, for comparison with the tracer.
Usage: AXIS_DIR=<dir with axis_<F>.txt> python tools/sanity/flight_mc.py
"""
import os
import sys
import numpy as np
sys.path.insert(0, os.path.dirname(__file__))
import range_check as rc   # noqa: E402

rng = np.random.default_rng(1)
m1, m2, M = rc.m1, rc.m2, rc.M


def run(f, E0, mode, N=2000, dxh=2.5e-4):
    xs, ns, _ = rc.axis_profile(f)
    nof = lambda x: np.interp(x, xs, ns, left=0.0, right=0.0)
    # cumulative n along x on a fine grid (beam is near-axial; uses |dx| projection of each flight)
    xg = np.arange(rc.X0, 0.06, 1e-5)
    edges = np.arange(-0.0224, 0.0784 + dxh, dxh); edep = np.zeros(len(edges) - 1)
    x1s, k1s = [], []
    for _ in range(N):
        x = np.array([rc.X0, 0.0, 0.0]); v = np.array([1.0, 0.0, 0.0]); E = E0; k = 0
        while E > 1.0 and x[0] < 0.0508:
            sig = rc.sigma(E, 23.0)
            tau = rng.exponential()
            if mode == "start":
                s = tau / (nof(x[0]) * sig)
            else:  # integrate n sigma along the straight flight until tau is used up
                s, ds, acc = 0.0, 5e-5, 0.0
                while acc < tau and s < 1.0:
                    a = nof(x[0] + v[0] * (s + ds / 2)) * sig * ds
                    if acc + a >= tau:
                        s += ds * (tau - acc) / a; acc = tau
                    else:
                        acc += a; s += ds
            x = x + v * s
            if x[0] > 0.0508 or x[0] < rc.X0:
                break
            # isotropic CM collision with Ne at rest
            V = np.sqrt(2 * E) * v / np.sqrt(m1)        # units with E = m v^2 / 2, masses in u
            vcm = m1 * V / M
            c = 2 * rng.random() - 1; ph = 2 * np.pi * rng.random()
            nhat = np.array([c, np.sqrt(1 - c * c) * np.cos(ph), np.sqrt(1 - c * c) * np.sin(ph)])
            Vn = vcm + (m2 / M) * np.linalg.norm(V) * nhat
            En = 0.5 * m1 * Vn @ Vn
            i = np.searchsorted(edges, x[0]) - 1
            if 0 <= i < len(edep):
                edep[i] += E - En
            E = En; v = Vn / np.linalg.norm(Vn); k += 1
        if E <= 1.0:
            x1s.append(x[0]); k1s.append(k)
    xc = 0.5 * (edges[1:] + edges[:-1]); d = edep / N / dxh
    return np.mean(x1s), np.mean(k1s), len(x1s) / N, xc[np.argmax(d)], d.max()


for f, E in rc.CASES:
    out = []
    for mode in ("start", "exact"):
        x1, k1, fr, xp, dp = run(f, E, mode)
        out.append(f"{mode}: x_1eV {x1*1e3:5.1f} mm, N {k1:4.1f}, peak at {xp*1e3:5.2f} mm, {dp/1e3:5.1f} eV/mm")
    print(f"{f:3d} SCCM {E:5d} eV | " + " | ".join(out))
