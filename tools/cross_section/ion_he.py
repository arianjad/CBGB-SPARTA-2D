"""BaF+ - He cross sections for the tracer (helium counterpart of ion.py/ion4.py/tracer_table.py
/write_tracer_table.py, same method, plus quantum scattering at low energy).

Proxy potential: Ba+-He (the F atom is ignored, as for neon).
  r >= RJ : Tang-Toennies form V = A exp(-b r) - f4(br) C4/r^4 - f6(br) C6/r^6
            fitted to the CCSD(T)/CBS Ba+-He potential of Buchachenko & Viehland,
            J. Chem. Phys. 148, 154304 (2018), Table IV: sigma 4.225 A, Re 4.921 A,
            De 22.4 cm^-1, C4 = 0.71 a.u. (their fitted value; asymptotic alpha_He/2 = 0.69),
            the same recipe as the Ba+-Ne table (ion.py); A, b, C6 fitted to (sigma, Re, De).
  r <= ra : ZBL universal repulsion for Ba-He (scaled by s), exponential bridge to RJ.
Q1(E): quantum partial waves (Numerov) for E_cm <= E_Q, classical deflection integral above
(they agree to ~1% at the join). Variants central/soft/stiff: ra = 1.5/1.2/1.8 A, ZBL x1/0.5/2,
as for neon.
Outputs (written to tracer/data/):
  bafplus_he_q1.tsv            bare Q1(E_cm)
  bafplus_he_sigma_tracer.tsv  tracer-effective sigma_t(E_cm) at T = 4 K (see tracer_table.py)
Checks printed: Langevin limit, thermal Omega11 and reduced mobility vs the polarization limit,
nuclear stopping vs ZBL universal stopping at 0.25-5 keV lab.
Takes a minute or two. Run from the repo root: python tools/cross_section/ion_he.py
"""
import math
import numpy as np

EH_K, A0 = 315775.02, 0.529177
HB2 = 24.2537                          # hbar^2/(2 amu kB A^2) [K A^2]
m_ion, m_gas = 156.325, 4.002602
MU = m_ion * m_gas / (m_ion + m_gas)
C4 = 0.71 * EH_K * A0**4               # K A^4 (B&V 2018 Table IV, CBS)
SIG0, RE, DE = 4.225, 4.921, 22.4 * 1.438777       # A, A, K (B&V 2018 Table IV, CBS)
C6 = None                              # fitted (fit_ab)
RJ = 3.3                               # A: TT form used for r >= RJ
T_CELL = 4.0
E_Q = 150.0                            # K: quantum below, classical above


def tt(n, x):
    s = np.zeros_like(x); term = np.ones_like(x)
    for k in range(n + 1):
        if k:
            term = term * x / k
        s = s + term
    return 1 - np.exp(-x) * s


def W(r, b, c6=None):
    c6 = C6 if c6 is None else c6
    return tt(4, b * r) * C4 / r**4 + tt(6, b * r) * c6 / r**6


def fit_ab():
    """A from V(SIG0) = 0; b and C6 on a grid to match (Re, De), as ion.fit_tt does for Ne."""
    global C6
    r = np.linspace(3.5, 9.0, 5501); best = None
    for b in np.linspace(1.0, 4.5, 351):
        for c6 in np.linspace(0, 400, 201) * EH_K * A0**6:
            A = W(np.array([SIG0]), b, c6)[0] * math.exp(b * SIG0)
            v = A * np.exp(-b * r) - W(r, b, c6); i = np.argmin(v)
            err = ((r[i] - RE) / 0.01) ** 2 + ((-v[i] - DE) / 0.1) ** 2
            if best is None or err < best[0]:
                best = (err, A, b, r[i], -v[i], c6)
    C6 = best[5]
    return best[:5]


def V_zbl(r, Z1=56, Z2=2):
    a = 0.8854 * A0 / (Z1**0.23 + Z2**0.23); x = r / a
    phi = 0.1818*np.exp(-3.2*x) + 0.5099*np.exp(-0.9423*x) + 0.2802*np.exp(-0.4029*x) + 0.02817*np.exp(-0.2016*x)
    return Z1 * Z2 * 14.3996 / r * phi * 11604.5 - C4 / r**4 * tt(4, 2.0 * r)


def make_V(A, b, ra, s):
    Vm = lambda r: A * np.exp(-np.minimum(b * r, 700)) - W(r, b)
    v1 = Vm(np.array([RJ]))[0]; v0 = s * V_zbl(np.array([ra]))[0]
    k = math.log(v0 / v1) / (RJ - ra)
    def V(r):
        r = np.asarray(r, float)
        return np.where(r >= RJ, Vm(r), np.where(r <= ra, s * V_zbl(np.maximum(r, 1e-6)), v1 * np.exp(k * (RJ - r))))
    return V


def Q1_classical(E, V, rg, nb=700, ns=160):
    bL = (4 * C4 / E) ** 0.25
    Vr = V(rg)
    rE = rg[np.where(Vr > E)[0].max()] if np.any(Vr > E) else rg[0]
    bmax = 4 * max(bL, rE, 1.0)
    bs = np.linspace(1e-4, bmax, nb)
    sg, w = np.polynomial.legendre.leggauss(ns); sg = 0.5 * (sg + 1); w = 0.5 * w
    chi = np.empty(nb)
    for i, bb in enumerate(bs):
        F = 1 - Vr / E - bb * bb / rg**2
        neg = np.where(F <= 0)[0]; j = neg.max() if len(neg) else 0
        if j + 1 >= len(rg):
            chi[i] = 0; continue
        r1, r2, F1, F2 = rg[j], rg[j + 1], F[j], F[j + 1]
        rc = r1 - F1 * (r2 - r1) / (F2 - F1) if F2 != F1 else r2
        u = 1 - sg * sg; r = rc / u
        Fr = np.maximum(1 - V(r) / E - bb * bb / r**2, 1e-14)
        chi[i] = np.pi - 2 * bb / rc * np.sum(w * 2 * sg / np.sqrt(Fr))
    return 2 * np.pi * np.trapezoid((1 - np.cos(chi)) * bs, bs)


def ricc(lmax, x):
    """Riccati-Bessel jhat_l = x j_l(x), nhat_l = x y_l(x), l = 0..lmax+1."""
    y = np.empty(lmax + 2); y[0] = -np.cos(x); y[1] = -np.cos(x) / x - np.sin(x)
    for l in range(1, lmax + 1):
        y[l + 1] = (2 * l + 1) / x * y[l] - y[l - 1]
    start = max(lmax + 60, int(x) + 60)
    j = np.zeros(start + 2); j[start] = 1e-300
    for l in range(start, 0, -1):
        j[l - 1] = (2 * l + 1) / x * j[l] - j[l + 1]
        if abs(j[l - 1]) > 1e250:
            j[l - 1:] *= 1e-250
    j *= np.sin(x) / j[0]
    return j[:lmax + 2], y[:lmax + 2]


def Q1_quantum(E, V):
    """sigma_D = 4 pi/k^2 sum_l (l+1) sin^2(d_l - d_{l+1}), distinguishable particles."""
    A = MU / HB2; k = math.sqrt(A * E)
    bL = (4 * C4 / E) ** 0.25
    lmax = int(k * max(4 * bL, 12.0)) + 40
    L = np.arange(lmax + 2)
    rmin = 2.6                                     # V(2.6 A) >> any E used here
    rmax = max(40.0, 4 * bL, (lmax + 2) / k * 1.5)
    h = min(0.004, 0.05 / k)
    r = np.arange(rmin, rmax + h, h)
    Wr = A * V(r) - k * k
    cent = L * (L + 1.0); h12 = h * h / 12.0
    u_prev = np.zeros(lmax + 2); u_cur = np.full(lmax + 2, 1e-30)
    f_prev = cent / r[0]**2 + Wr[0]; f_cur = cent / r[1]**2 + Wr[1]
    i1 = len(r) - 30; u1 = None
    for i in range(1, len(r) - 1):
        f_next = cent / r[i + 1]**2 + Wr[i + 1]
        u_next = ((2 + 10 * h12 * f_cur) * u_cur - (1 - h12 * f_prev) * u_prev) / (1 - h12 * f_next)
        s = np.abs(u_next)
        if s.max() > 1e100:
            s = np.where(s > 1e100, s, 1.0); u_next /= s; u_cur /= s
            if u1 is not None:
                u1 /= s
        u_prev, u_cur = u_cur, u_next; f_prev, f_cur = f_cur, f_next
        if i + 1 == i1:
            u1 = u_next.copy()
    u2 = u_cur
    j1, n1 = ricc(lmax, k * r[i1]); j2, n2 = ricc(lmax, k * r[-1])
    # long-range r^-4 tail at r[-1] is weak: phase shifts from the matching are adequate
    R = u1 / u2
    d = np.arctan2(R * j2 - j1, R * n2 - n1)
    l = np.arange(len(d) - 1)
    return 4 * np.pi / k**2 * np.sum((l + 1) * np.sin(d[:-1] - d[1:])**2)


def main():
    import json, os
    best = fit_ab(); _, A, b, re, de = best
    print(f"TT fit: A={A:.4g} K, b={b:.3f}/A -> Re={re:.3f} A, De={de:.2f} K ({de/1.438777:.2f} cm^-1); "
          f"targets {RE:.3f} A, {DE:.2f} K")
    V0 = make_V(A, b, 1.5, 1.0)
    rr = np.linspace(2.5, 6, 3501); vv = V0(rr); s0 = rr[np.where(vv > 0)[0].max()]
    print(f"sigma0 (V=0) = {s0:.3f} A (target {SIG0}); C6 = {C6/(EH_K*A0**6):.0f} a.u.")
    rg = np.concatenate([np.linspace(0.05, 2, 4000), np.linspace(2, 80, 26000)[1:], np.geomspace(80, 4000, 3000)[1:]])
    variants = {"central": (1.5, 1.0), "low": (1.2, 0.5), "high": (1.8, 2.0)}
    Vs = {k: make_V(A, b, ra, s) for k, (ra, s) in variants.items()}
    # checks: pure polarization Langevin ratio (classical), quantum vs classical overlap
    Vp = lambda r: 1e9 * (1.0 / np.asarray(r))**12 - C4 / np.asarray(r)**4
    for E in [0.01, 0.1]:          # low E: the r^-12 core must be negligible against the r^-4 capture radius
        print(f"check pure r^-4 classical: E={E} K Q1/sigma_L = {Q1_classical(E, Vp, rg)/(2*np.pi*math.sqrt(C4/E)):.3f} (expect ~1.105)")
    Eq = np.geomspace(0.05, 400, 45)
    qq = np.array([Q1_quantum(E, V0) for E in Eq])
    for E in [20, 50, 100, 150, 250, 400]:
        qc = Q1_classical(E, V0, rg); qi = np.exp(np.interp(np.log(E), np.log(Eq), np.log(qq)))
        print(f"  E={E:5.0f} K  Q1 quantum {qi:8.2f}  classical {qc:8.2f} A^2  ratio {qi/qc:.3f}")
    Ec = np.geomspace(E_Q, 10**7.5, 50)
    rows_E = np.concatenate([Eq[Eq < E_Q], Ec])
    Q = {}
    for k, V in Vs.items():
        qc = np.array([Q1_classical(E, V, rg) for E in Ec])
        qlow = qq[Eq < E_Q] if k == "central" else np.array([Q1_quantum(E, V) for E in Eq[Eq < E_Q]])
        Q[k] = np.concatenate([qlow, qc])
    out = "tracer/data/bafplus_he_q1.tsv"
    with open(out, "w") as f:
        f.write("# BaF+ - He BARE momentum-transfer cross section Q1(E) -- reference data; the tracer uses bafplus_he_sigma_tracer.tsv.\n")
        f.write(f"# Ba+-He proxy potential: TT form fitted to the CCSD(T)/CBS Ba+-He potential of Buchachenko & Viehland,\n")
        f.write(f"#   JCP 148, 154304 (2018), Table IV (sigma 4.225 A, Re 4.921 A, De 22.4 cm^-1, C4 0.71 a.u.; fit gives Re {re:.3f} A,\n")
        f.write(f"#   De {de/1.438777:.2f} cm^-1, C6 {C6/(EH_K*A0**6):.0f} a.u.) for r >= 3.3 A,\n")
        f.write("#   ZBL Ba-He repulsion inside r_a (central 1.5 A x1, low 1.2 A x0.5, high 1.8 A x2), exponential bridge.\n")
        f.write(f"# Q1: quantum partial waves (Numerov) below E_cm = {E_Q:g} K, classical above. Generated by tools/cross_section/ion_he.py.\n")
        f.write("# Units: E_cm in K, Q1 in m^2.\n# E_cm_K\tQ1_central_m2\tQ1_low_m2\tQ1_high_m2\n")
        for i, E in enumerate(rows_E):
            f.write("%.6e\t%.6e\t%.6e\t%.6e\n" % (E, Q["central"][i] * 1e-20, Q["low"][i] * 1e-20, Q["high"][i] * 1e-20))
    print("wrote", out, len(rows_E), "rows")
    # thermal Omega11, mobility
    lE = np.log(rows_E)
    Qf = lambda k: (lambda Ek: np.exp(np.interp(np.log(np.maximum(Ek, rows_E[0])), lE, np.log(Q[k] * 1e-20))))
    x = np.linspace(1e-4, 40, 8000)
    kB, amu, e = 1.380649e-23, 1.66053906660e-27, 1.602176634e-19
    K_pol = 13.876 / math.sqrt(0.2050 * MU)
    for T in [1, 2, 4, 10, 23, 77, 300]:
        Om = 0.5 * np.trapezoid(x**2 * np.exp(-x) * Qf("central")(x * T), x)
        nD = 3 / 16 * math.sqrt(2 * math.pi * kB * T / (MU * amu)) / Om
        K0 = e * nD / 2.6868e25 / (kB * T) * 1e4
        print(f"  T={T:4g} K: Omega11 {Om:.3e} m^2, reduced mobility K0 {K0:.2f} cm^2/Vs (polarization limit {K_pol:.2f})")
    # nuclear stopping vs ZBL universal
    Z1, Z2, M1, M2 = 56, 2, m_ion, m_gas
    gam = 4 * M1 * M2 / (M1 + M2)**2
    for El in [250, 500, 1000, 2500, 5000]:
        Ecm = El * M2 / (M1 + M2) * 11604.5
        S_ours = gam / 2 * El * Qf("central")(Ecm)                              # eV m^2
        eps = 32.53 * M2 * (El / 1000) / (Z1 * Z2 * (M1 + M2) * (Z1**0.23 + Z2**0.23))   # ZBL reduced energy
        sn = math.log(1 + 1.1383 * eps) / (2 * (eps + 0.01321 * eps**0.21226 + 0.19593 * eps**0.5)) if eps <= 30 else math.log(eps) / (2 * eps)
        S_zbl = 8.462e-15 * Z1 * Z2 * M1 * sn / ((M1 + M2) * (Z1**0.23 + Z2**0.23)) * 1e-4      # eV m^2
        print(f"  {El:5d} eV lab: S_n ours {S_ours:.3e}  ZBL {S_zbl:.3e} eV m^2  ratio {S_ours/S_zbl:.3f}")
    write_tracer_table(rows_E, Q, A, b)


def write_tracer_table(E, Q, A, b):
    """As tracer_table.py + write_tracer_table.py, for He at T = 4 K."""
    kB = 8314.46; m = m_gas; T = T_CELL; mu = MU
    lE = np.log(E)
    s0 = {k: (np.log(Q[k][1]) - np.log(Q[k][0])) / (lE[1] - lE[0]) for k in Q}
    def Q1fun(k):
        lQ = np.log(Q[k] * 1e-20)
        return lambda Ek: np.exp(np.where(np.log(np.maximum(Ek, 1e-12)) < lE[0],
                                          lQ[0] + s0[k] * (np.log(np.maximum(Ek, 1e-12)) - lE[0]),
                                          np.interp(np.log(np.maximum(Ek, 1e-12)), lE, lQ)))
    rng = np.random.default_rng(1)
    vg = rng.normal(0, math.sqrt(kB * T / m), size=(2_000_000, 3))
    vref = math.sqrt(8 * kB * T / (math.pi * m))
    Vs = np.concatenate([[0.5], np.geomspace(2, 152000, 102)])
    keys = ("central", "low", "high")
    tab = {k: [] for k in keys}; Ec = []
    for V in Vs:
        gv = np.array([V, 0, 0]) - vg; g = np.linalg.norm(gv, axis=1); gp = gv[:, 0]
        vc = math.sqrt(vref**2 + V**2); Ec.append(mu * vc**2 / (2 * kB))
        for k in keys:
            q = Q1fun(k)(mu * g**2 / (2 * kB))
            tab[k].append(np.mean(q * g * gp) * np.mean(g) / (vc * np.mean(g * gp)))
    Ec = np.array(Ec)
    x = np.linspace(1e-4, 40, 8000)
    def rest_limit(TT, k):
        sb = 0.5 * np.trapezoid(x**2 * np.exp(-x) * Q1fun(k)(x * TT), x)
        return sb * math.sqrt(8 * kB * TT / (math.pi * mu)) / math.sqrt(8 * kB * TT / (math.pi * m))
    Ec0 = lambda TT: mu * (8 * kB * TT / (math.pi * m)) / (2 * kB)
    rows = [(Ec0(TT), *[rest_limit(TT, k) for k in keys]) for TT in np.geomspace(0.3, T, 14)]
    for i in np.where(Ec > Ec0(T) * 1.08)[0]:
        rows.append((Ec[i], *[tab[k][i] for k in keys]))
    rows.sort()
    hdr = (f"# BaF+ - He TRACER-EFFECTIVE cross section sigma_t(E_cm) for ParticleTracing3D.jl (--sigma-model table,\n"
           f"# pass with --sigma-table and -m {m_gas}). Same construction as bafplus_ne_sigma_tracer.tsv: the tracer draws\n"
           "# collisions at n*sigma*v_coll with a flux-weighted Maxwellian partner and isotropic CM scattering, so\n"
           "#     sigma_t(V) = <Q1(g) g g_par> <g> / (v_coll <g g_par>),  averages over the He Maxwellian,\n"
           "# tabulated vs E_cm = mu*v_coll^2/2, v_coll = sqrt(8kT/(pi m_He) + V^2). Equals Q1 for fast ions and the\n"
           f"# Chapman-Enskog thermal value at rest. Rows with E_cm < {Ec0(T)*1.08:.2f} K: exact V=0 limit for local T = 0.3-{T:g} K;\n"
           f"# rows above: computed at T = {T:g} K. Q1 from tracer/data/bafplus_he_q1.tsv (tools/cross_section/ion_he.py).\n"
           "# E_cm_K\tsigma_central_m2\tsigma_low_m2\tsigma_high_m2\n")
    out = "tracer/data/bafplus_he_sigma_tracer.tsv"
    with open(out, "w") as f:
        f.write(hdr)
        for r in rows:
            f.write("%.6e\t%.6e\t%.6e\t%.6e\n" % r)
    arr = np.array(rows)
    print("wrote", out, len(rows), "rows; E range %.3g-%.3g K; V->0 check sigma_t(0.5 m/s) = %.3e, rest limit %.3e"
          % (arr[0, 0], arr[-1, 0], tab["central"][0], rest_limit(T, "central")))


if __name__ == "__main__":
    main()
