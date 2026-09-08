#!/usr/bin/env python
"""Score one 2D tracer run at a specified observation plane.

    python tools/tracer_common.py RUN_DIR CELL.SURFS --xobs 0.0650 \
        --json RUN_DIR/common.json

The reported counts separate:

    N_first    reached x_obs going forward at least once
    N_back     ... and later crossed back upstream
    N_survive  = N_first - N_back
    N_escape   left the domain past its geometry's extraction plane

`Y(R_ap, theta)` is the accepted yield at the FIRST forward crossing:
v_z > 0, r_cross < R_ap, |v_perp|/v_z < tan(theta).

Inputs come from `tracer/accumulators/crossing.jl`: `recs.csv` carries the
hook-side (legs 1..N-1) crossing state and `legs.out` the terminating leg.
The terminating leg is handled here because the hook does not see it.
"""
import argparse
import hashlib
import json
import sys
from pathlib import Path

import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parent))
from tracer_analyze import COLS, clip_legs, read_surfs  # noqa: E402

R_GRID_MM = (2.5, 5.0, 10.0, 19.05)
THETA_GRID_MRAD = (50.0, 100.0, 200.0, 400.0)
T_GATES_MS = (5.0, 20.0)


def terminal_hits(d, surfs):
    """Return clipped terminal axial contacts and the run's extraction plane.

    The accumulator hook has already seen every non-terminal leg. The final
    leg is clipped against the projected segments that terminate the 2D tracer.
    """
    seg, xap, bounds = read_surfs(surfs)
    hit, _ = clip_legs(
        np.column_stack([d["z"], np.hypot(d["x"], d["y"])]),
        np.column_stack([d["znext"], np.hypot(d["xnext"], d["ynext"])]),
        seg, bounds)
    return hit[:, 0], xap


def score(rundir, surfs, xobs, ref=None):
    rundir = Path(rundir)
    d = pd.read_csv(rundir / "legs.out", sep=r"\s+")
    if list(d.columns) != COLS:
        sys.exit(f"unexpected columns in {rundir/'legs.out'}: {list(d.columns)}")
    r = pd.read_csv(rundir / "recs.csv")
    n = len(d)
    if len(r) != n:
        sys.exit(f"{len(r)} accumulator records for {n} trajectory rows")
    if not (r["idx"].to_numpy() == d["idx"].to_numpy()).all():
        sys.exit("recs.csv and legs.out are not in the same particle order")

    # --- the geometry's own extraction plane ----------------------------------
    z_hit, xap = terminal_hits(d, surfs)
    escape = z_hit > xap
    escape_obs = z_hit > xobs

    # --- terminating leg vs the common plane --------------------------------
    # The leg runs from (x,y,z) to the UNCLIPPED free-path end; it is really
    # travelled only as far as the clipped contact point z_hit.  So the plane
    # is reached iff x_obs lies between the leg start and that contact point.
    z0, z1 = d["z"].to_numpy(), d["znext"].to_numpy()
    dz = z1 - z0
    with np.errstate(divide="ignore", invalid="ignore"):
        f = (xobs - z0) / dz
    fwd_term = (z0 < xobs) & (z_hit >= xobs) & (dz > 0)
    back_term = (z0 > xobs) & (z_hit <= xobs) & (dz < 0)
    f = np.where(np.isfinite(f), f, 0.0)

    vx, vy, vz = d["vx"].to_numpy(), d["vy"].to_numpy(), d["vz"].to_numpy()
    xc_t = d["x"].to_numpy() + f * (d["xnext"].to_numpy() - d["x"].to_numpy())
    yc_t = d["y"].to_numpy() + f * (d["ynext"].to_numpy() - d["y"].to_numpy())
    leg_len = np.linalg.norm(
        np.column_stack([d["xnext"] - d["x"], d["ynext"] - d["y"], dz]), axis=1)
    speed = np.linalg.norm(np.column_stack([vx, vy, vz]), axis=1)
    # `time` on the row is cumulative at the START of the terminating leg
    # (propagate returns before incrementing it).
    tc_t = d["time"].to_numpy() + f * leg_len / speed

    # --- union: the hook wins, because its legs all precede leg N -----------
    have_hook = r["have"].to_numpy().astype(bool)
    first = have_hook | fwd_term
    xc = np.where(have_hook, r["cx"], xc_t)
    yc = np.where(have_hook, r["cy"], yc_t)
    vxc = np.where(have_hook, r["vx"], vx)
    vyc = np.where(have_hook, r["vy"], vy)
    vzc = np.where(have_hook, r["vz"], vz)
    tc = np.where(have_hook, r["tc"], tc_t)
    ncc = np.where(have_hook, r["nc"], d["collides"])

    # A back-crossing on the terminating leg only counts if a forward crossing
    # already happened on an earlier leg; a forward crossing on leg N cannot be
    # followed by anything.
    nback = r["nback"].to_numpy() + (have_hook & back_term)
    back = nback > 0
    survive = first & ~back

    rc = np.hypot(xc, yc)
    vperp = np.hypot(vxc, vyc)
    with np.errstate(divide="ignore", invalid="ignore"):
        tanth = np.where(vzc > 0, vperp / vzc, np.inf)

    res = {
        "rundir": str(rundir),
        "surfs": str(surfs),
        "x_obs_m": xobs,
        "escape_plane_m": xap,
        "n": int(n),
        "N_first": int(first.sum()),
        "N_back": int(back.sum()),
        "N_survive": int(survive.sum()),
        "N_escape": int(escape.sum()),
        "N_escape_at_xobs": int(escape_obs.sum()),
        "N_escape_after_back": int((escape_obs & back).sum()),
        "n_back_crossings": int(nback.sum()),
        "first_from_hook": int(have_hook.sum()),
        "first_from_terminating_leg": int((~have_hook & fwd_term).sum()),
        "max_endpoint_mismatch_m": float(r["mism"].max()),
        "t_cross_ms": {
            "median": float(np.median(tc[first]) * 1e3) if first.any() else float("nan"),
            "p10": float(np.percentile(tc[first], 10) * 1e3) if first.any() else float("nan"),
            "p90": float(np.percentile(tc[first], 90) * 1e3) if first.any() else float("nan"),
        },
        "vz_cross_mps_median": float(np.median(vzc[first])) if first.any() else float("nan"),
        "collides_at_cross_median": float(np.median(ncc[first])) if first.any() else float("nan"),
        "yield": {},
        "yield_time_gated": {},
    }

    def y_of(mask_r, mask_a, extra=None):
        m = first & mask_r & mask_a & (vzc > 0)
        if extra is not None:
            m = m & extra
        k = int(m.sum())
        p = k / n
        return {"n": k, "Y": p, "stderr": float(np.sqrt(p * (1 - p) / n))}

    allr = np.ones(n, bool)
    res["yield"]["inf/inf"] = y_of(allr, allr)
    for R in R_GRID_MM:
        mr = rc < R * 1e-3
        res["yield"][f"{R}/inf"] = y_of(mr, allr)
        for th in THETA_GRID_MRAD:
            res["yield"][f"{R}/{th}"] = y_of(mr, tanth < np.tan(th * 1e-3))
    for tg in T_GATES_MS:
        e = tc < tg * 1e-3
        res["yield_time_gated"][f"inf/inf/t<{tg}ms"] = y_of(allr, allr, e)
        res["yield_time_gated"][f"2.5/inf/t<{tg}ms"] = y_of(rc < 2.5e-3, allr, e)

    res["checks"] = {
        "endpoint_identity": bool(res["max_endpoint_mismatch_m"] <= 1e-12),
        "first_ge_survive": bool(res["N_first"] >= res["N_survive"]),
    }
    if ref is not None:
        ref = Path(ref)
        ref_rows = ref / "legs.out"
        if not ref_rows.is_file():
            sys.exit(f"reference run has no legs.out: {ref_rows}")
        res["reference_rows"] = _same_rows(ref_rows, rundir / "legs.out")
    return res


def _same_rows(reference, mine):
    """Return hashes of index-sorted rows for an optional reference run."""
    def h(p):
        rows = [ln for ln in Path(p).read_text().splitlines() if ln[:1].isdigit()]
        rows.sort(key=lambda s: int(s.split()[0]))
        return hashlib.md5("\n".join(rows).encode()).hexdigest()
    a, b = h(reference), h(mine)
    return {"reference_md5": a, "run_md5": b, "identical": a == b}


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("rundir", help="directory containing legs.out and recs.csv")
    ap.add_argument("surfs", help="2D tracer geometry (.surfs)")
    ap.add_argument("--xobs", type=float, default=0.0650)
    ap.add_argument("--ref", default=None,
                    help="optional run directory whose legs.out is compared")
    ap.add_argument("--json", default=None)
    a = ap.parse_args()
    res = score(a.rundir, a.surfs, a.xobs, a.ref)
    print(json.dumps(res, indent=2))
    if a.json:
        Path(a.json).write_text(json.dumps(res, indent=2))


if __name__ == "__main__":
    main()
