#!/usr/bin/env python3
"""Further sanity checks of the cyl-3d ports runs, independent of the tracer code
(only the cross-section tables and the frozen DSMC neon field are shared):

  1. range scaling   depth to 1 eV x chamber density should depend on energy only
  2. time to 1 eV    sum of free path / speed along the slowing-down, vs the tracer's t_1eV
  3. diffusion       Chapman-Enskog D of BaF+ in Ne from the bare Q1 table (and the
                     Langevin mobility limit), vs D measured from the 40-ion GIF trajectories
  4. extraction      drift-diffusion model of the thermalized ions in the DSMC gas
                     flow (drift = gas velocity + grad D, D from check 3), started from where a
                     slowing-down Monte Carlo leaves them, vs the tracer's extracted/back/wall
  5. neon flow       outflow/inflow balance and chamber density vs free-molecular and
                     continuum (choked) aperture limits

Inputs: AXIS_DIR with axis_<F>.txt (x n T on the beam axis) and field2d_<F>.txt
(x r n u_x u_r T, azimuthally averaged in 0.4 mm bins), both cut from frozen.grid with
awk (commands in tools/sanity/README.md). Light: runs in a few minutes on a
login node. Writes results/sanity/more_checks.txt (+ printed).
"""
import glob
import json
import os
import sys

import numpy as np
import pandas as pd

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import range_check as rc   # noqa: E402

BASE = rc.BASE
AX = os.environ["AXIS_DIR"]
OUT = "results/sanity/more_checks.txt"
kB, amu, eV = rc.kB, rc.amu, rc.eV
m1, m2, M, MU = rc.m1, rc.m2, rc.M, rc.MU
FLOWS = [10, 15, 20, 30, 40]
CASES = rc.CASES
L_CH, R_CH = 45.47e-3, 15.24e-3                 # chamber x 0..L, bore radius
R_IN, X_IN = 1.27e-3, -1.27e-3                  # inlet bore radius, outer face
R_OUT, X_OUT = 1.5875e-3, 50.804e-3             # outlet channel radius, outer face
SCCM = 4.4775e17                                # particles/s per SCCM
rng = np.random.default_rng(7)
lines = []


def w(s=""):
    print(s); lines.append(s)


def chamber_n(f):
    return rc.chamber_n(f)


def ions(f, E):
    return pd.read_csv(f"{BASE}/sccm{f}/post/tracer_table-central-e{E}_s1_ions.csv", comment="#")


def out_file(f, E):
    return np.loadtxt(f"{BASE}/sccm{f}/post/tracer_table-central-e{E}_s1.out", skiprows=1)


# ---------- Q1 -> Chapman-Enskog n*D(T) ----------
_q = np.array([l.split() for l in open("tracer/data/bafplus_ne_q1.tsv") if l.strip() and not l.startswith("#")], float)


def omega11(T):
    x = np.logspace(np.log10(_q[0, 0]), np.log10(_q[-1, 0]), 20000)
    q = np.exp(np.interp(np.log(x), np.log(_q[:, 0]), np.log(_q[:, 1])))
    return 0.5 / T**3 * np.trapezoid(x**2 * q * np.exp(-x / T), x)


_Tg = np.linspace(1.0, 40.0, 79)
_nD = np.array([3 / 16 * np.sqrt(2 * np.pi * kB * T / MU) / omega11(T) for T in _Tg])


def nD(T):
    return np.interp(T, _Tg, _nD)


# ---------- slowing-down MC (3D), tracer's start-point rule or exact ----------
def slow_down(f, E0, N, mode="start", Emin=1.0):
    xs, ns, Ts = rc.axis_profile(f)
    nof = lambda x: np.interp(x, xs, ns, left=0.0, right=0.0)
    P = np.zeros((N, 3)); t = np.zeros(N); k = np.zeros(N, int); ok = np.zeros(N, bool)
    for i in range(N):
        x = np.array([X_IN, 0.0, 0.0]); v = np.array([1.0, 0.0, 0.0]); E = E0; tt = 0.0; kk = 0
        while E > Emin:
            sig = rc.sigma(E, 23.0); tau = rng.exponential(); spd = np.sqrt(2 * E * eV / (m1 * amu))
            if mode == "start":
                s = tau / max(nof(x[0]) * sig, 1e-30)
            else:
                s, ds, acc = 0.0, 5e-5, 0.0
                while acc < tau and s < 1.0:
                    a = nof(x[0] + v[0] * (s + ds / 2)) * sig * ds
                    if acc + a >= tau:
                        s += ds * (tau - acc) / a; acc = tau
                    else:
                        acc += a; s += ds
            x = x + v * s; tt += s / spd
            if not (X_IN <= x[0] <= X_OUT) or np.hypot(x[1], x[2]) > R_CH:
                break
            V = spd * v; vcm = m1 * V / M
            c = 2 * rng.random() - 1; ph = 2 * np.pi * rng.random()
            nh = np.array([c, np.sqrt(1 - c * c) * np.cos(ph), np.sqrt(1 - c * c) * np.sin(ph)])
            Vn = vcm + (m2 / M) * spd * nh
            E = 0.5 * m1 * amu * (Vn @ Vn) / eV; v = Vn / np.linalg.norm(Vn); kk += 1
        P[i], t[i], k[i], ok[i] = x, tt, kk, E <= Emin
    return P, t, k, ok


# ---------- 1. range scaling ----------
def check_range():
    w("1. RANGE SCALING: depth from the beam start to 1 eV, times chamber density (1e20 m^-2).")
    w("   Same energy at different flows should give the same number if only the density matters.")
    Es = [500, 1000, 2000, 3000, 5000]
    w("   " + "SCCM".rjust(5) + "".join(f"{E:>9d}" for E in Es) + "   (eV)")
    tab = {}
    for f in [5] + FLOWS:
        n = chamber_n(f); row = []
        for E in Es:
            try:
                d = ions(f, E)
            except FileNotFoundError:
                row.append(None); continue
            ok = d.x_1eV.notna()
            row.append((d.x_1eV[ok].mean() - X_IN) * n / 1e20 if ok.mean() > 0.9 else None)
        tab[f] = row
        w("   " + f"{f:5d}" + "".join(f"{v:9.2f}" if v is not None else f"{'(<90%)':>9}" for v in row))
    for j, E in enumerate(Es):
        v = [tab[f][j] for f in tab if tab[f][j] is not None]
        if len(v) > 1:
            w(f"   {E} eV: spread {100*(max(v)-min(v))/np.mean(v):.0f}% over {len(v)} flows")
    w("   Expected from the slowing-down MC with the tracer's flight rule (depth x n, 1e20 m^-2):")
    for f, E in CASES[:3]:
        P, t, k, ok = slow_down(f, E, 600, "start")
        w(f"     {f} SCCM {E} eV: MC {(P[ok,0].mean()-X_IN)*chamber_n(f)/1e20:.2f}")
    w()


# ---------- 2. time to 1 eV ----------
def check_time():
    w("2. TIME TO SLOW DOWN TO 1 eV (microseconds): slowing-down MC (sum of flight / speed) vs tracer.")
    w(f"   {'case':>16} {'MC tracer rule':>15} {'MC exact':>9} {'tracer mean':>12} {'tracer median':>14}")
    for f, E in CASES:
        Ps, ts, ks, oks = slow_down(f, E, 600, "start")
        Pe, te, ke, oke = slow_down(f, E, 600, "exact")
        d = ions(f, E); ok = d.t_1eV.notna()
        w(f"   {f:3d} SCCM {E:5d} eV {1e6*ts[oks].mean():15.2f} {1e6*te[oke].mean():9.2f} "
          f"{1e6*d.t_1eV[ok].mean():12.2f} {1e6*d.t_1eV[ok].median():14.2f}")
    w()


# ---------- 3. diffusion ----------
def check_diffusion():
    T = 23.0
    K0_pol = 13.876 / np.sqrt(0.3956 * (m1 * m2 / M))         # cm^2/Vs, Langevin limit (alpha_Ne in A^3)
    nd = nD(T); K = eV * (nd / 2.6868e25) / (kB * T) * 1e4    # reduced mobility K0 = e D N/(kT N0)
    w("3. THERMAL DIFFUSION of BaF+ in 23 K neon.")
    w(f"   Chapman-Enskog from Q1: Omega11(23 K) = {omega11(T):.3e} m^2, n*D = {nd:.3e} /(m s);")
    w(f"   reduced mobility K0 = {K:.2f} cm^2/Vs vs Langevin (pure polarization) limit {K0_pol:.2f} cm^2/Vs "
      f"(ratio {K/K0_pol:.2f}).")
    w("   Measured from the GIF trajectories (all energies of a flow, 40 ions each): lateral mean-square")
    w("   displacement over 0.2-3 ms after t = 0.3 ms, MSD_perp = 4 D dt, ions still inside.")
    w(f"   {'SCCM':>5} {'n_ch':>9} {'D theory':>9} {'D measured':>11} {'ratio':>6} {'ions':>5}")
    for f in FLOWS:
        files = sorted(glob.glob(f"{BASE}/sccm{f}/post/gif/*.npz"))
        dts, msd = [], []
        for p in files:
            z = np.load(p); tm, pos, al = z["times"], z["pos"], z["alive"]
            i0 = np.searchsorted(tm, 3e-4)
            for j in range(i0 + 1, len(tm)):
                dt = tm[j] - tm[i0]
                if dt < 2e-4 or dt > 3e-3:
                    continue
                m = al[:, i0] & al[:, j]
                if m.sum() < 5:
                    continue
                d2 = ((pos[m, j, 1:] - pos[m, i0, 1:]) ** 2).sum(1)
                dts.append(np.full(m.sum(), dt)); msd.append(d2)
        if not dts:
            continue
        dts, msd = np.concatenate(dts), np.concatenate(msd)
        Dm = (msd * dts).sum() / (4 * (dts ** 2).sum())          # least squares through origin
        n = chamber_n(f); Dt = nD(T) / n
        w(f"   {f:5d} {n:9.2e} {Dt:9.2e} {Dm:11.2e} {Dm/Dt:6.2f} {len(files)*40:5d}")
    w()


# ---------- 4. extraction: drift-diffusion of thermal ions in the DSMC flow ----------
def load_field(f):
    a = np.loadtxt(f"{AX}/field2d_{f}.txt")
    x0, h = -0.0224, 4e-4
    ix = np.round((a[:, 0] - x0) / h - 0.5).astype(int); ir = np.round(a[:, 1] / h - 0.5).astype(int)
    nx, nr = ix.max() + 1, ir.max() + 1
    G = {k: np.full((nx, nr), np.nan) for k in ("n", "u", "ur", "T")}
    for k, c in (("n", 2), ("u", 3), ("ur", 4), ("T", 5)):
        G[k][ix, ir] = a[:, c]
    D = nD(np.clip(G["T"], 1.0, 40.0)) / G["n"]
    D[~np.isfinite(D)] = np.nan
    Dx = np.gradient(np.nan_to_num(D, nan=np.nanmedian(D)), h, axis=0)
    Dr = np.gradient(np.nan_to_num(D, nan=np.nanmedian(D)), h, axis=1)
    for k in ("u", "ur"):
        G[k] = np.nan_to_num(G[k])
    G.update(D=np.nan_to_num(D, nan=np.nanmedian(D)), Dx=Dx, Dr=Dr, x0=x0, h=h, nx=nx, nr=nr)
    return G


def inside(x, r):
    return (((x >= 0) & (x <= L_CH) & (r < R_CH)) | ((x >= X_IN) & (x < 0) & (r < R_IN))
            | ((x > L_CH) & (x < X_OUT) & (r < R_OUT)))


def drift_diffusion(G, P0, tmax=1.0):
    N = len(P0); P = P0.copy(); t = np.zeros(N); fate = np.zeros(N, int)   # 0 inside, 1 extr, 2 back, 3 wall
    act = np.arange(N)
    r0 = np.hypot(P[:, 1], P[:, 2]); bad = ~inside(P[:, 0], r0); fate[bad] = 3; act = act[~bad]
    while act.size:
        x, y, z = P[act, 0], P[act, 1], P[act, 2]; r = np.hypot(y, z)
        i = np.clip(((x - G["x0"]) / G["h"]).astype(int), 0, G["nx"] - 1)
        j = np.clip((r / G["h"]).astype(int), 0, G["nr"] - 1)
        D, ux, ur = G["D"][i, j], G["u"][i, j] + G["Dx"][i, j], G["ur"][i, j] + G["Dr"][i, j]
        near = (x < 2e-3) | (x > L_CH - 3e-3)
        hs = np.where(near, 5e-5, 2e-4)                      # step length scale (m)
        dt = np.minimum(0.25 * hs**2 / (2 * D), 0.25 * hs / (np.hypot(ux, ur) + 1e-9))
        rr = np.maximum(r, 1e-12)
        s = np.sqrt(2 * D * dt)
        P[act, 0] = x + ux * dt + s * rng.standard_normal(act.size)
        P[act, 1] = y + ur * y / rr * dt + s * rng.standard_normal(act.size)
        P[act, 2] = z + ur * z / rr * dt + s * rng.standard_normal(act.size)
        t[act] += dt
        xn = P[act, 0]; rn = np.hypot(P[act, 1], P[act, 2])
        out = ~inside(xn, rn)
        fe = np.where((xn >= X_OUT) & (rn < R_OUT), 1, np.where((xn < X_IN) & (rn < R_IN), 2, 3))
        fate[act[out]] = fe[out]
        act = act[~out & (t[act] < tmax)]
    return fate, t


def tracer_fates(f, E):
    a = out_file(f, E)
    xn, r = a[:, 4], np.hypot(a[:, 5], a[:, 6])
    from tracer_classify3d import CYL_BOX, classify   # noqa: E402
    tr, lost, wall = classify(a, rc_outlet(), CYL_BOX)
    return tr, lost, wall, a[:, 11]


def rc_outlet():
    from tracer_classify3d import CYL_PORTS_OUTLET_X
    return CYL_PORTS_OUTLET_X


def check_extraction():
    w("4. EXTRACTION: drift-diffusion model vs tracer (% of injected ions; median transit of extracted, ms).")
    w("   Model: ions slowed to 1 eV by the MC with the tracer's flight rule, then drift with the")
    w("   azimuthally averaged DSMC gas velocity (+ grad D) and diffuse with Chapman-Enskog D(n, T).")
    w(f"   {'case':>16} | {'extracted model':>15} {'tracer':>7} | {'back model':>10} {'tracer':>7} | "
      f"{'wall model':>10} {'tracer':>7} | {'t_med model':>11} {'tracer':>7}")
    for f, E in CASES:
        G = load_field(f)
        P, ts, ks, ok = slow_down(f, E, 1000, "start")
        fate, t = drift_diffusion(G, P)
        fate[~ok] = 3                                      # left before thermalizing: count as wall
        tr, lost, wall, tt = tracer_fates(f, E)
        ext = fate == 1
        w(f"   {f:3d} SCCM {E:5d} eV | {100*ext.mean():15.1f} {100*tr.mean():7.1f} | "
          f"{100*(fate==2).mean():10.1f} {100*lost.mean():7.1f} | {100*(fate==3).mean():10.1f} {100*wall.mean():7.1f} | "
          f"{1e3*np.median(t[ext]) if ext.any() else np.nan:11.1f} {1e3*np.median(tt[tr]):7.1f}"
          f"   (still inside at 1 s: {100*(fate==0).mean():.1f}%)")
    w()


# ---------- 5. neon flow ----------
def check_neon():
    s = open("results/dashboard/site/index.html").read()
    i = s.index("const DATA = ") + 13
    data, _ = json.JSONDecoder().raw_decode(s[i:])
    flows = {fl["f"]: fl for fl in data["configs"][0]["flows"]}
    m = m2 * amu; T0 = 23.0
    cbar = np.sqrt(8 * kB * T0 / (np.pi * m)); g = 5 / 3; c0 = np.sqrt(g * kB * T0 / m)
    A_in, A_out = np.pi * R_IN**2, np.pi * R_OUT**2
    A_eff = 6.6e-6
    choke = (2 / (g + 1)) ** ((g + 1) / (2 * (g - 1)))      # n* c* / (n0 c0) for an ideal nozzle
    W_cont = choke * c0 / (cbar / 4)                        # choked orifice flux / free-molecular flux
    d_ref, om = 4.59e-10, 0.883
    eta = 5 / 16 * np.sqrt(np.pi * m * kB * T0) / (np.pi * d_ref**2) * 1.0   # HS-like viscosity estimate
    w("5. NEON FLOW. Q = inflow; ratio = DSMC outflow/inflow at steady state; n_FM = 4Q/(c A_eff)")
    w(f"   (free-molecular, Clausing-weighted A_eff = 6.6 mm^2); n_cont = Q/(W A_in + W A_out) with an")
    w(f"   ideal choked orifice W = {W_cont:.2f} x free-molecular per area (no wall friction: lower bound).")
    w("   lambda = Ne-Ne mean free path (VSS d_ref 4.59e-10 m at 23 K); Kn = lambda / outlet diameter;")
    w("   Re = outlet Reynolds number at the DSMC exit-plane density and speed.")
    w(f"   {'SCCM':>5} {'ratio':>6} {'n DSMC':>9} {'n_FM':>9} {'n_cont':>9} {'DSMC/FM':>8} {'DSMC/cont':>9} "
      f"{'lambda um':>9} {'Kn':>6} {'Re':>5}")
    for f in sorted(flows):
        Q = f * SCCM; n = chamber_n(f)
        nFM = 4 * Q / (cbar * A_eff); nC = Q / (W_cont * (cbar / 4) * (A_in + A_out))
        lam = 1 / (np.sqrt(2) * n * np.pi * d_ref**2)
        Re = ""
        if os.path.exists(f"{AX}/axis_{f}.txt") and f in FLOWS:
            G = load_field(f); i = int((X_OUT - G["x0"]) / G["h"])
            ne = np.nanmean(G["n"][i, :3]); ue = np.nanmean(G["u"][i, :3])
            Re = f"{ne * m * ue * 2 * R_OUT / eta:5.0f}"
        w(f"   {f:5d} {flows[f]['steady'].get('ratio') or float('nan'):6.3f} {n:9.3e} {nFM:9.3e} {nC:9.3e} "
          f"{n/nFM:8.2f} {n/nC:9.2f} {lam*1e6:9.0f} {lam/(2*R_OUT):6.3f} {Re:>5}")
    w(f"   (viscosity used for Re: {eta:.2e} Pa s)")
    w()


def main():
    os.makedirs("results/sanity", exist_ok=True)
    which = sys.argv[1:] or ["1", "2", "3", "4", "5"]
    for k, fn in (("1", check_range), ("2", check_time), ("3", check_diffusion), ("4", check_extraction), ("5", check_neon)):
        if k in which:
            fn()
    open(OUT if which == ["1", "2", "3", "4", "5"] else OUT.replace(".txt", f"_{''.join(which)}.txt"), "w").write("\n".join(lines) + "\n")


if __name__ == "__main__":
    main()
