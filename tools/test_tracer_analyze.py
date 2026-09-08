#!/usr/bin/env python
"""Self-check for tools/tracer_analyze.py -- synthetic inputs, known answers.

    python tools/test_tracer_analyze.py

Every assertion is two-sided: each one has a construction that would fail it.
"""
import sys
import tempfile
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parent))
import tracer_analyze as ta

SURFS = """ITEM: TIMESTEP
0
ITEM: NUMBER OF SURFS
4
ITEM: BOX BOUNDS oo ao pp
-0.02 0.1
0.0 0.0199
-0.5 0.5
ITEM: SURFS id v1x v1y v2x v2y
1 0.0 0.00635 0.0 0.002
2 0.053 0.00635 0.0 0.00635
3 0.053 0.0025 0.053 0.00635
4 0.0535 0.0025 0.053 0.0025
"""


def write_out(path, rows):
    with open(path, "w") as f:
        f.write(" ".join(ta.COLS) + "\n")
        for r in rows:
            f.write(("%d " + " ".join(["%e"] * 9) + " %d %e\n") % tuple(r))


def main():
    rng = np.random.default_rng(7)
    tmp = Path(tempfile.mkdtemp())
    fs = tmp / "cell.surfs"
    fs.write_text(SURFS)
    ok = 0

    # 1. geometry: the extraction plane is the max vertex x, not the box xhi.
    seg, xap, bounds = ta.read_surfs(fs)
    assert len(seg) == 4, seg.shape
    assert abs(xap - 0.0535) < 1e-12, xap
    assert np.allclose(bounds, [[-0.02, 0.1], [0.0, 0.0199]]), bounds
    print(f"PASS 1  read_surfs: 4 segments, x_ap = {xap} (box xhi 0.1 correctly ignored)")
    ok += 1

    # 1b. leg clipping: the tracer hands back the UNCLIPPED free-path endpoint,
    #     so each of the three termination modes must be pulled back to its
    #     real contact point.  Wall / downstream face / radial face.
    p1 = np.array([[0.030, 0.0040],      # into the cell wall at r = 0.00635
                   [0.050, 0.0010],      # straight out through zhi = 0.1
                   [0.060, 0.0050]])     # steeply out through rhi = 0.0199
    p2 = np.array([[0.032, 0.0300],
                   [0.120, 0.0020],
                   [0.070, 0.0300]])
    pt, code = ta.clip_legs(p1, p2, seg, bounds)
    assert list(code) == [1, 2, 2], code
    assert abs(pt[0, 1] - 0.00635) < 1e-9, pt[0]
    assert abs(pt[1, 0] - 0.1) < 1e-9, pt[1]
    assert abs(pt[2, 1] - 0.0199) < 1e-9, pt[2]
    assert pt[0, 1] < p2[0, 1] and pt[1, 0] < p2[1, 0], "clipping did nothing"
    print(f"PASS 1b clip_legs: wall r={pt[0,1]:.5f}, zhi z={pt[1,0]:.4f}, "
          f"rhi r={pt[2,1]:.4f} (raw endpoints were {p2[0,1]:.3f}/{p2[1,0]:.3f}/{p2[2,1]:.3f})")
    ok += 1

    # 2. FWHM of a Gaussian is 2.355 sigma.  A 3-sigma-wide sample would fail.
    w = ta.fwhm(rng.normal(0.0, 12.0, 400000))
    assert abs(w - 2.3548 * 12.0) < 0.04 * 2.3548 * 12.0, w
    print(f"PASS 2  fwhm(N(0,12)) = {w:.2f} vs 2.355*sigma = {2.3548*12:.2f} (< 4%)")
    ok += 1

    # 3. a small sample returns nan rather than a fabricated width
    assert np.isnan(ta.fwhm(rng.normal(0, 1, 10)))
    print("PASS 3  fwhm of a 10-sample set is nan, not a number")
    ok += 1

    # 4. extraction counts znext > x_ap only; equality is NOT extraction.
    n = 50000
    rows = []
    for i in range(1, 20001):                   # extracted, past the plane
        rows.append([i, 0.001, 0.0, 0.05, 0.002, 0.0, 0.09,
                     rng.normal(0, 20), rng.normal(0, 20), 140.0, 5000, 2e-3])
    for i in range(20001, 20501):               # exactly on the plane -> not extracted
        rows.append([i, 0.001, 0.0, 0.05, 0.002, 0.0, 0.0535,
                     0.0, 0.0, 5.0, 100, 1e-4])
    for i in range(20501, 35001):               # wall deaths (a --saveall run)
        rows.append([i, 0.004, 0.0, 0.04, 0.00635, 0.0, 0.04,
                     0.0, 0.0, 1.0, 900, 5e-4])
    fo = tmp / "syn.out"
    write_out(fo, rows)
    d, ext, other, seg, xap, res = ta.analyze(fo, fs, nparticles=n)
    assert res["n_extracted"] == 20000, res["n_extracted"]
    assert abs(res["extraction"] - 0.400) < 1e-12, res["extraction"]
    assert res["n_wall_rows"] == 15000, res["n_wall_rows"]   # 14500 wall + 500 on-plane
    print(f"PASS 4  extraction {res['extraction']:.3f} of n=50000 "
          f"(500 on-plane rows counted as non-extracted, "
          f"{res['n_wall_rows']} non-extracted rows kept)")
    ok += 1

    # 4b. n_not_extracted is an alias of n_wall_rows (rename, same value);
    #     every synthetic extracted row sits on-axis (r constant, no radial
    #     crossing), so n_axial_exit / extraction_axial equal n_extracted /
    #     extraction exactly, and the axial stderr formula matches the
    #     un-suffixed one.
    assert res["n_not_extracted"] == res["n_wall_rows"] == 15000, res["n_not_extracted"]
    assert res["n_axial_exit"] == 20000, res["n_axial_exit"]
    assert abs(res["extraction_axial"] - res["extraction"]) < 1e-12, (
        res["extraction_axial"], res["extraction"])
    assert abs(res["extraction_axial_stderr"] - res["extraction_stderr"]) < 1e-12, (
        res["extraction_axial_stderr"], res["extraction_stderr"])
    print(f"PASS 4b n_not_extracted aliases n_wall_rows ({res['n_not_extracted']}); "
          f"n_axial_exit {res['n_axial_exit']} / extraction_axial "
          f"{res['extraction_axial']:.3f} match n_extracted/extraction "
          f"(no radial exits in this synthetic set)")
    ok += 1

    # 5. --n defaults to max idx and is flagged as such
    _, _, _, _, _, r2 = ta.analyze(fo, fs, nparticles=None)
    assert r2["n"] == 35000 and r2["n_from_max_idx"] is True
    print("PASS 5  n falls back to max idx (35000) and sets n_from_max_idx")
    ok += 1

    # 6. divergence: v_f = 140, transverse sigma = 20 -> 2*atan(2.355*20/280).
    #    Tolerance from the measured estimator rms at n = 20000 (1.8%): 6% is 3 sigma.
    want = 2e3 * np.arctan(2.3548 * 20.0 / (2 * 140.0))
    got = res["divergence_fwhm_mrad"]
    assert abs(got - want) < 0.06 * want, (got, want)
    print(f"PASS 6  divergence {got:.0f} mrad vs analytic {want:.0f} mrad (< 6%)")
    ok += 1

    # 6b. the adaptive bin count is what makes small samples usable: at n = 400
    #     a fixed 60-bin histogram reads ~15% low.  Head-to-head on one sample,
    #     so both outcomes are reachable and neither needs a magic threshold.
    truth = 2.3548 * 20.0
    small = rng.normal(0.0, 20.0, 400)
    e_adaptive = abs(ta.fwhm(small) - truth)
    e_fixed = abs(ta.fwhm(small, nbins=60) - truth)
    assert e_adaptive < e_fixed, (e_adaptive, e_fixed)
    print(f"PASS 6b adaptive bins beat fixed-60 at n=400: "
          f"err {e_adaptive:.2f} vs {e_fixed:.2f} m/s (true FWHM {truth:.2f})")
    ok += 1

    # 7. figures render, including the no-wall-death branch
    made = ta.figures(d, ext, other, seg, xap, res, tmp, "syn")
    assert made[0].exists() and made[0].stat().st_size > 10000
    ext2 = ext.copy()
    ta.figures(d, ext2, other.iloc[0:0], seg, xap, res, tmp, "syn_nowall")
    assert (tmp / "syn_nowall_summary.png").exists()
    print(f"PASS 7  summary figure {made[0].stat().st_size} bytes; "
          f"empty-wall branch also renders")
    ok += 1

    # 8. heat map: spawns near the aperture extract, corner spawns do not
    import pandas as pd
    sp = pd.DataFrame({"idx": [r[0] for r in rows],
                       "x": [0.0045 if r[0] > 400 else 0.0005 for r in rows],
                       "y": 0.0,
                       "z": [0.01 if r[0] > 400 else 0.045 for r in rows]})
    fsp = tmp / "spawn.csv"
    sp.to_csv(fsp, index=False)
    p = ta.heatmap(ext, fsp, seg, res, tmp, "syn", nz=4, nr=4)
    assert p.exists() and p.stat().st_size > 5000
    print(f"PASS 8  heat map rendered from a spawn CSV ({p.stat().st_size} bytes)")
    ok += 1

    print(f"\n{ok}/11 PASS")
    return 0


if __name__ == "__main__":
    sys.exit(main())
