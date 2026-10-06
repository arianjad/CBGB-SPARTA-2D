# Per-flow numerics for the plain-cylinder neon scan (cases/cyl-3d), shared by
# sbatch_cyl3d_neon.sh, sbatch_cyl3d_post.sh and submit_cyl3d_campaign.sh.
# Source it, then call
#   cyl_flow_params <SCCM>
# CONFIG selects the port geometry (gen_cyl3d.py --config), default ports:
#   ports  2.54 mm inlet (1.27 mm wall) + 3.175 mm x 5.334 mm outlet (the fluor-cell-3d
#          ports), data in cyl-ports. Expected density 1.26e21 m^-3 per SCCM (measured
#          on the old cell: 2.15e22 at 17 SCCM). Fixed 0.4 mm lattice (KNMIN 0): a
#          cell <= 2 lambda grid would cost ~5M core-h at these densities; the old
#          17 SCCM test found 0.41 and 0.20 mm grids agree on extraction. dt < tau
#          stays strict. Runs stop at steady state (in.ne_cyl3d_run); SIM_S = 0.5 s is
#          only the cap (gas residence time ~0.13 s) and sizes the walltime request.
#   ports-he  the same ports geometry with helium-4 at 4 K (cases/cyl-3d/he.vss), data in
#          cyl-ports-he. Expected density: the neon value scaled by the mean-speed ratio
#          (c_bar He 4 K 145 m/s vs Ne 23 K 155 m/s): 1.34e21 m^-3 per SCCM. Same 0.4 mm lattice,
#          same strict dt < tau.
#   ap5    5 mm apertures through 1.27 mm plates, data in cyl-scan (2026-10-03 run,
#          stopped after 5-15 SCCM): 2.9e20 per SCCM, KNMIN 0.5, 0.1 s.
# It also sets KN (the KNMIN in effect), SIMS (simulated seconds), SURF, XOUT (outlet face x, m), AEFF_MM2 (Clausing-weighted port area
# for the aperture-law check) and RUNBASE_DEF, and the gas: GAS (Ne/He), SID, SPECIESF, VSSF,
# TGAS (K), MGAS_U (u), TFLOOR (K, floor on the sampled T in the grid-rule check).
# It sets MDOT (kg/s), FILLN, FNUM, DT, REFINE, CELL_MM, ADAPT_STEPS, NADAPT,
# WINDOW (averaging window = checkpoint interval, steps), TOTAL (production
# steps for SIM_S simulated seconds, default 0.1 s), and a job size
# NP / MEM_GB / TIME for sbatch, from the expected particle count.
#
# Grid rule, checked in every cell: dt < local mean collision time (strict, as in
# Steinmacher et al., arXiv:2609.11544), and cell < lambda / KNMIN with
# KNMIN = 0.5, i.e. cells up to 2 local mean free paths (chosen 2026-10-02:
# the paper's cell < lambda costs ~47k core-h for 5-50 SCCM with the corrected
# neon cross section, vs ~6k; our 17 SCCM grid test found extraction unchanged
# between 2.2 and 4.4 lambda). The lattice is therefore chosen per flow:
# - Expected chamber density: molecular-flow aperture law (their eq. 2) for
#   two 5 mm thin apertures, n = 4 Q / (c_bar A) = 2.9e20 m^-3 per SCCM
#   (c_bar = 155 m/s at 23 K). In the transitional regime the apertures
#   conduct more, so this is an upper bound on the density, and the cells
#   and time step picked from it are on the safe side.
# - Mean free path and collision time at that density from the VHS
#   parameters in cases/cyl-3d/ne.vss (the same formulas as SPARTA's
#   compute lambda/grid): lambda_est, tau_est.
# - Starting chamber cell: 1.6 mm / REFINE with REFINE = ceil(1.6 mm x KNMIN / lambda_est),
#   at least 4 (cells <= 0.40 mm, so a 5 mm aperture spans >= 12 cells).
#   in.ne_cyl3d_run then refines any cell that still has Kn < KNMIN (e.g.
#   denser spots) and stops with an error if the rule can't be met.
# - DT = 0.5 x min(tau_est, cell / c_bar), 2 significant figures, rounded down;
#   in.ne_cyl3d_run stops if DT >= the smallest per-cell tau it measures.
# - FNUM: PPC (default 5) particles per starting chamber cell at the expected density.
# - Adaptation passes: ADAPT_STEPS ~ 2.5 ms each, at most NADAPT = 6.
# - WINDOW: ~10 ms per averaging window (one field.grid frame), a multiple of
#   5000 steps (the stats / boundary-flux averaging interval, so a resumed run
#   starts on a step where those are defined). Only the last NDUMP (default 3)
#   windows write field.grid frames (storage).
cyl_flow_params() {
  local f=$1
  CONFIG="${CONFIG:-ports}"
  GAS=Ne; SID=Ne; SPECIESF=ne.species; VSSF=ne.vss; TGAS=23.0; MGAS_U=20.1797; TFLOOR=3.0
  case "$CONFIG" in
    ports-he) NPER=1.34e21; KN="${KNMIN:-0}"; SIMS="${SIM_S:-0.5}"; SURF=cell_cyl3d_ports.surf; XOUT=0.050804
           AEFF_MM2=6.6; RUNBASE_DEF=cyl-ports-he
           GAS=He; SID=He; SPECIESF=he.species; VSSF=he.vss; TGAS=4.0; MGAS_U=4.002602; TFLOOR=1.0 ;;
    ports) NPER=1.26e21; KN="${KNMIN:-0}"; SIMS="${SIM_S:-0.5}"; SURF=cell_cyl3d_ports.surf; XOUT=0.050804
           AEFF_MM2=6.6; RUNBASE_DEF=cyl-ports ;;
    ap5)   NPER=2.9e20; KN="${KNMIN:-0.5}"; SIMS="${SIM_S:-0.1}"; SURF=cell_cyl3d.surf; XOUT=0.04674
           AEFF_MM2=39.27; RUNBASE_DEF=cyl-scan ;;
    *) echo "unknown CONFIG $CONFIG" >&2; return 2 ;;
  esac
  read MDOT FILLN FNUM DT REFINE CELL_MM ADAPT_STEPS WINDOW NPART_EST NCELL_EST <<<"$(python3 - "$f" "${PPC:-5}" "$KN" "$NPER" "$(dirname "${BASH_SOURCE[0]}")/../../cases/cyl-3d/$VSSF" "$SID" "$MGAS_U" "$TGAS" <<'EOF'
import math, sys
f, ppc, knmin, nper, vss = float(sys.argv[1]), float(sys.argv[2]), float(sys.argv[3]), float(sys.argv[4]), sys.argv[5]
sid, mu_gas, T = sys.argv[6], float(sys.argv[7]), float(sys.argv[8])
d, w, tref = [float(x) for x in next(l for l in open(vss) if l.split() and l.split()[0] == sid).split()[1:4]]
kB, m = 1.380649e-23, mu_gas * 1.66053907e-27
n = nper * f
cbar = math.sqrt(8 * kB * T / (math.pi * m))
lam = 1 / (math.pi * math.sqrt(2) * d**2 * n * (tref / T)**(w - 0.5))       # SPARTA lambda/grid, one species
mr = m / 2
tau = 1 / (2 * d**2 * n * math.sqrt(2 * math.pi * kB * tref / mr) * (T / tref)**(1 - w))
R = max(4, math.ceil(1.6e-3 * knmin / lam))     # cell <= lambda / KNMIN
h = 1.6e-3 / R
def sig2(x):
    e = math.floor(math.log10(x)); return math.floor(x / 10**(e - 1)) * 10**(e - 1)
dt = sig2(0.5 * min(tau, h / cbar))
fnum = n * h**3 / ppc
ncell = 33178e-9 / h**3                       # chamber cells
npart = n * 33178e-9 / fnum
k1000 = lambda s: max(1000, int(round(s / dt / 1000)) * 1000)
k5000 = lambda s: max(5000, int(round(s / dt / 5000)) * 5000)
print(f"{1.500541e-8*f*mu_gas/20.1797:.6e} {n:.4e} {fnum:.4e} {dt:.2g} {R} {h*1e3:.4f} {k1000(2.5e-3)} {k5000(0.01)} {npart:.3e} {ncell:.3e}")
EOF
)"
  NADAPT="${NADAPT:-6}"
  DT="${DT_SET:-$DT}"; FNUM="${FNUM_SET:-$FNUM}"
  if [ -n "${DT_SET:-}" ]; then
    WINDOW=$(python3 -c "print(max(5000, int(round(0.01/$DT/5000))*5000))")
    ADAPT_STEPS=$(python3 -c "print(max(1000, int(round(2.5e-3/$DT/1000))*1000))")
  fi
  WINDOW="${WINDOW_SET:-$WINDOW}"
  TOTAL=$(python3 -c "print(max(1, int(round($SIMS/$DT/$WINDOW)))*$WINDOW)")
  # job size: ~1.5M particles per task (16..48 tasks); walltime from the smoke-test rate
  # (3.35M particle-steps per core-second on submit95) x 2 margin; memory ~0.7 GB per
  # 1M particles + 0.5 GB per 1M cells (smoke test), doubled for refined cells
  read NP MEM_GB TIME <<<"$(python3 -c "
import math
npart, ncell, steps = $NPART_EST, $NCELL_EST, $TOTAL + ($NADAPT + 1) * $ADAPT_STEPS
np_ = min(48, max(16, 16 * math.ceil(npart / 1.5e6 / 16)))
while np_ < 48 and npart * steps / 3.35e6 / np_ / 3600 > 36:   # keep runs under ~1.5 days
    np_ += 16
mem = int(max(32, 2 * (0.7 * npart / 1e6 + 0.5 * (ncell + 1.4e6) / 1e6) + 8))
hrs = min(144, max(4, math.ceil(2 * npart * steps / 3.35e6 / np_ / 3600)))
print(np_, mem, f'{hrs//24}-{hrs%24:02d}:00:00')")"
}
