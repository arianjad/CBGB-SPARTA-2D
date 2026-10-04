#!/usr/bin/env python3
"""Build the plain-cylinder scan dashboard as a small static site for the Artifact page:

  results/dashboard/site/index.html     the page; all run data embedded as JSON
  results/dashboard/site/fig/*.png      scan figures (tools/plot_cyl3d_scan.py)
  results/dashboard/site/run/F_E.png    per-run figure: dE/dx vs depth, stopping-depth
                                        histogram, exit-energy histogram
  results/dashboard/site/run/F_E.gif    trajectory GIF, downscaled to 560 px wide

The page lists one row per flow; a click expands it into its energies, and a click
on an energy shows that run's GIF, figure and numbers. Images are separate files so
the page stays small; only new or changed files need re-uploading:
results/dashboard/site/to_publish.json lists them (published path -> local path).
After a successful publish run with --mark-published to record them.

Light work (small CSVs, incremental image conversion), fine on the login node.
Usage: python tools/build_cyl3d_dashboard.py [--runbase cyl-scan] [--cond table-central] [--mark-published]
"""
import argparse
import datetime
import glob
import hashlib
import json
import os
import re
import subprocess
import sys

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from tracer_classify3d import CYL_BOX, CYL_OUTLET_X, CYL_PORTS_OUTLET_X, classify, load  # noqa: E402

BASE = os.path.expanduser("~/orcd/scratch/cbgb-sparta-3d")
SITE = "results/dashboard/site"
FLOWS = [5, 10, 15, 20, 25, 30, 35, 40, 45, 50]
# port configurations (gen_cyl3d.py --config); the first is shown when the page loads
CONFIGS = [
    dict(key="ports", label="Old-cell ports",
         desc="2.54 mm inlet through a 1.27 mm wall, 3.175 mm outlet through a 5.33 mm channel (the fluor-cell-3d ports). "
              "Fixed 0.4 mm DSMC lattice; time step below the local collision time everywhere.",
         runbase="cyl-ports", job="ports", fig="results/figures/cyl3d_ports", run="run/ports", outlet=CYL_PORTS_OUTLET_X,
         nen=20, energies="0.25–5 keV in 250 eV steps"),
    dict(key="ap5", label="5 mm apertures",
         desc="5 mm apertures through 1.27 mm plates (stopped 2026-10-03 after 5–15 SCCM: the gas is ~5x thinner than the "
              "old cell at the same flow, so ions crossed without stopping). Lattice ≤ 2 mean free paths per flow.",
         runbase="cyl-scan", job="cyl", fig="results/figures/cyl3d_scan", run="run", outlet=CYL_OUTLET_X,
         nen=25, energies="0.2–5 keV in 200 eV steps"),
]
AMU, EV = 1.66054e-27, 1.602177e-19
GIF_W = 560
INK, MUTED, C1 = "#1f1f1e", "#6b6a66", "#2a78d6"


def sh(cmd):
    try:
        return subprocess.run(cmd, shell=True, capture_output=True, text=True, timeout=60).stdout
    except Exception:
        return ""


def wilson(k, n, z=1.96):
    if n == 0:
        return [None, None, None]
    p = k / n; d = 1 + z * z / n
    c = (p + z * z / (2 * n)) / d; h = z * np.sqrt(p * (1 - p) / n + z * z / (4 * n * n)) / d
    return [round(100 * p, 2), round(100 * (c - h), 2), round(100 * (c + h), 2)]


def x1_of(paths):
    import pandas as pd
    return np.concatenate([pd.read_csv(p, comment="#").x_1eV.to_numpy(dtype=float) for p in paths])


def run_stats(d, cond, E, outlet):
    paths = sorted(glob.glob(f"{d}/post/tracer_{cond}-e{E}_s*.out"))
    r = load(paths)
    if len(r) == 0:
        return None
    tr, lost, wall = classify(r, outlet, CYL_BOX)
    ke = 0.5 * 156.325 * AMU * np.sum(r[:, 7:10] ** 2, axis=1) / EV
    th = tr & (ke < 0.1)
    s = dict(E=E, n=len(r), extr=wilson(tr.sum(), len(r)), therm=wilson(th.sum(), len(r)),
             back=round(100 * lost.mean(), 2), wall=round(100 * wall.mean(), 2),
             t_med_ms=round(float(np.median(r[tr, 11])) * 1e3, 3) if tr.any() else None,
             ke_med_eV=round(float(np.median(ke[tr])), 3) if tr.any() else None)
    ions = sorted(glob.glob(f"{d}/post/tracer_{cond}-e{E}_s*_ions.csv"))
    if ions:
        x1 = x1_of(ions)
        ok = np.isfinite(x1)
        s["stop_frac"] = round(100 * ok.mean(), 1)
        s["stop_mm"] = round(float(np.mean(x1[ok])) * 1e3, 2) if ok.any() else None
    return s, r, tr, ke


def run_figure(d, cond, E, f, r, tr, ke, out):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    import pandas as pd
    fig, ax = plt.subplots(1, 3, figsize=(12, 3.4))
    # dE/dx (1 mm bins)
    xr = sorted(glob.glob(f"{d}/post/tracer_{cond}-e{E}_s*_xr.csv"))
    if xr:
        tot, nions = None, 0
        for p in xr:
            nions += int(re.search(r"over all (\d+) ions", open(p).readline()).group(1))
            g = pd.read_csv(p, comment="#").groupby("x_m").edep_eV.sum()
            tot = g if tot is None else tot.add(g, fill_value=0)
        b = np.floor(np.asarray(tot.index) * 1e3).astype(int); xs = np.unique(b)
        y = np.array([tot.values[b == v].sum() for v in xs]) / nions
        ax[0].plot(xs + 0.5, y, color=C1, lw=2)
    ax[0].axvspan(-1.27, 0, color="#b9b8b2", alpha=0.35, lw=0); ax[0].axvline(45.47, color=MUTED, lw=1, ls=":")
    ax[0].set_xlim(-2, 50); ax[0].set_xlabel("depth x (mm)"); ax[0].set_ylabel("dE/dx per ion (eV/mm)")
    ax[0].set_title("Energy loss vs depth", fontsize=10)
    # stopping depth
    ions = sorted(glob.glob(f"{d}/post/tracer_{cond}-e{E}_s*_ions.csv"))
    if ions:
        x1 = x1_of(ions)
        ok = np.isfinite(x1)
        if ok.any():
            ax[1].hist(x1[ok] * 1e3, bins=np.arange(-2, 48, 1), color=C1, edgecolor="white", linewidth=0.5)
        ax[1].set_title(f"Where ions drop below 1 eV ({100*ok.mean():.0f}% ever do)", fontsize=10)
    ax[1].set_xlabel("depth x (mm)"); ax[1].set_ylabel("ions")
    # exit energy of extracted ions
    if tr.any():
        ax[2].hist(np.log10(np.maximum(ke[tr], 1e-4)), bins=np.arange(-4, 4.1, 0.1), color=C1,
                   edgecolor="white", linewidth=0.5)
        ax[2].axvline(-1, color=INK, lw=1, ls="--"); ax[2].text(-0.95, ax[2].get_ylim()[1] * 0.9, "0.1 eV", fontsize=8)
    ax[2].set_xlabel("log10(exit KE / eV), extracted ions"); ax[2].set_ylabel("ions")
    ax[2].set_title("Exit energy", fontsize=10)
    for a in ax:
        a.grid(alpha=0.25)
    fig.suptitle(f"{f} SCCM, {E/1000:g} keV", fontsize=11)
    fig.tight_layout(); fig.savefig(out, dpi=90); plt.close(fig)


def small_gif(src, out):
    """Downscale a trajectory GIF to GIF_W wide. All frames share one palette built
    from the first, middle and last frames, without dithering, so the small fate
    colours (blue extracted, orange wall, green back out) survive: quantizing each
    frame on its own let the large grey areas take the palette and mapped those
    dots to brown (fixed 2026-10-04)."""
    from PIL import Image, ImageSequence
    im = Image.open(src)
    rgb = [f.convert("RGB") for f in ImageSequence.Iterator(im)]   # convert while iterating (frames share one object)
    n = len(rgb)
    sample = Image.new("RGB", (im.width, im.height * 3))
    for k, i in enumerate((0, n // 2, n - 1)):
        sample.paste(rgb[i], (0, im.height * k))
    pal = sample.quantize(colors=255, method=Image.Quantize.MEDIANCUT)
    h = round(im.height * GIF_W / im.width)
    fr = [f.resize((GIF_W, h), Image.LANCZOS).quantize(palette=pal, dither=Image.Dither.NONE) for f in rgb]
    fr[0].save(out, save_all=True, append_images=fr[1:], duration=im.info.get("duration", 83), loop=0, optimize=True)


def flow_row(f, runbase, jobs, jobpfx):
    d = f"{BASE}/{runbase}/sccm{f}"
    log = open(f"{d}/scan.log").read() if os.path.exists(f"{d}/scan.log") else ""
    rule = open(f"{d}/post/rule_check.txt").read() if os.path.exists(f"{d}/post/rule_check.txt") else ""
    m = re.search(r"REFINE=\d+ \(([\d.]+) mm\)", log)
    st = jobs.get(f"{jobpfx}-ne{f}", "")
    if "solver exit code 0" in log and os.path.exists(f"{d}/post/frozen.grid"):
        neon = "done"
    elif st:
        neon = st.lower()
    elif "STRICT RULE VIOLATED" in log or re.search(r"solver exit code [1-9]", log):
        neon = "failed"
    elif "fresh start" in log or "resume from" in log:
        neon = "stopped"
    else:
        neon = "not started"
    mr = re.search(r"chamber <n> ([\d.e+-]+) m\^-3 vs aperture law \(eq. 2\) ([\d.e+-]+) \(([\d.]+)x\)", rule)
    sizes = re.search(r"cell sizes \[([^\]]*)\]", rule)
    prog = ""
    if neon == "running":
        steps = re.findall(r"^ +(\d+) +[\d.]+ +\d+ ", open(f"{d}/log.sparta").read(), re.M) if os.path.exists(f"{d}/log.sparta") else []
        tot = re.findall(r"TOTAL=(\d+)", log)
        if steps and tot:
            prog = f"{int(steps[-1])/int(tot[-1])*100:.0f}%"
    return dict(f=f, neon=neon, prog=prog, lattice=m.group(1) if m else "",
                lattices=" / ".join(x for x in sizes.group(1).split() if x != "1.6") if sizes else "",
                rule="pass" if "SATISFIED" in rule else ("fail" if "VIOLATED" in rule else ""),
                n=mr.group(1) if mr else "", ratio=mr.group(3) if mr else "", runs=[])


def steady_info(d, png, f, neon):
    """Steady-state diagnostics of one neon run: STEADY line from steady.txt, and a figure of the
    particle count and gas outflow / inflow vs simulated time from the run.log stats lines."""
    log = open(f"{d}/scan.log").read() if os.path.exists(f"{d}/scan.log") else ""
    par = {k: float(v) for k, v in re.findall(r"\b(DT|MDOT|FNUM)=([\d.e+-]+)", log)}
    st = open(f"{d}/steady.txt").read() if os.path.exists(f"{d}/steady.txt") else ""
    m = re.findall(r"(NOT )?STEADY at step (\d+) after (\d+) windows \(change ([\d.e+-]+)\)", st)
    info = dict(status="", t_ms=None, windows=None, change=None, png=None, ratio=None)
    if m:
        no, step, win, ch = m[-1]
        info.update(status="not steady" if no else "steady", windows=int(win), change=float(ch),
                    t_ms=round(int(step) * par.get("DT", 0) * 1e3, 1) if par.get("DT") else None)
    elif "STRICT RULE" in log or not os.path.exists(f"{d}/run.log"):
        return info
    if not os.path.exists(f"{d}/run.log") or "DT" not in par:
        return info
    rows = {}
    for line in open(f"{d}/run.log"):
        p = line.split()
        if len(p) == 13 and p[0].isdigit():
            try:
                v = [float(x) for x in p]
            except ValueError:
                continue
            if sum(v[7:13]) > 0:
                rows[int(v[0])] = v
    if not rows:
        return info
    a = np.array([rows[k] for k in sorted(rows)])
    t = a[:, 0] * par["DT"] * 1e3
    inflow = par["MDOT"] / (20.1797 * 1.66053907e-27 * par["FNUM"]) * par["DT"]     # particles per step
    ratio = a[:, 7:13].sum(axis=1) / inflow
    info["ratio"] = round(float(ratio[-1]), 3)
    if not info["status"]:
        info["status"] = {"done": "fixed run length", "running": "settling"}.get(neon, "stopped")
    if not os.path.exists(png) or os.path.getmtime(png) < os.path.getmtime(f"{d}/run.log"):
        import matplotlib
        matplotlib.use("Agg")
        import matplotlib.pyplot as plt
        fig, ax = plt.subplots(1, 2, figsize=(10, 2.8))
        ax[0].plot(t, a[:, 2] / 1e6, color=C1, lw=2, marker="o", ms=3)
        ax[0].set_ylabel("simulated particles (M)"); ax[0].set_title("Particle count", fontsize=10)
        ax[1].plot(t, ratio, color=C1, lw=2, marker="o", ms=3)
        ax[1].axhline(1, color=MUTED, lw=1, ls="--")
        ax[1].set_ylabel("outflow / inflow"); ax[1].set_title("Gas mass balance (5000-step averages)", fontsize=10)
        for x in ax:
            x.set_xlabel("simulated time (ms)"); x.grid(alpha=0.25)
            if info["t_ms"] is not None:
                x.axvline(info["t_ms"], color=INK, lw=1, ls=":")
        fig.suptitle(f"{f} SCCM neon: approach to steady state" + (" (dotted: steady)" if info["t_ms"] is not None else ""),
                     fontsize=10)
        fig.tight_layout(); fig.savefig(png, dpi=90); plt.close(fig)
    info["png"] = png
    return info


def fhash(p):
    st = os.stat(p)
    return f"{st.st_size}:{int(st.st_mtime)}"


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--cond", default="table-central")
    ap.add_argument("--mark-published", action="store_true")
    ap.add_argument("--redo-gifs", action="store_true", help="re-make every downscaled GIF (after a small_gif change)")
    a = ap.parse_args()
    manifest_p = f"{SITE}/.published.json"
    manifest = json.load(open(manifest_p)) if os.path.exists(manifest_p) else {}
    if a.mark_published:
        todo = json.load(open(f"{SITE}/to_publish.json"))
        for k, v in todo.items():
            manifest[k] = fhash(v)
        json.dump(manifest, open(manifest_p, "w"), indent=0)
        print(f"marked {len(todo)} files published"); return

    os.makedirs(f"{SITE}/fig", exist_ok=True); os.makedirs(f"{SITE}/run", exist_ok=True)
    jobs = {}
    for line in sh("squeue -u $USER -h -o '%j %T'").splitlines():
        p = line.split()
        if len(p) == 2:
            jobs[p[0]] = p[1]
    cfgs = []
    for c in CONFIGS:
        os.makedirs(f"{SITE}/{c['run']}", exist_ok=True)
        flows = []
        for f in FLOWS:
            row = flow_row(f, c["runbase"], jobs, c["job"])
            d = f"{BASE}/{c['runbase']}/sccm{f}"
            Es = sorted({int(re.search(r"-e(\d+)_s", p).group(1)) for p in glob.glob(f"{d}/post/tracer_{a.cond}-e*_s*.out")})
            for E in Es:
                res = run_stats(d, a.cond, E, c["outlet"])
                if res is None:
                    continue
                st, r, tr, ke = res
                png = f"{SITE}/{c['run']}/{f}_{E}.png"
                if not os.path.exists(png):
                    run_figure(d, a.cond, E, f, r, tr, ke, png)
                st["png"] = f"{c['run']}/{f}_{E}.png"
                src = f"{d}/post/gif/sccm{f}_{a.cond}_e{E}.gif"
                gif = f"{SITE}/{c['run']}/{f}_{E}.gif"
                if os.path.exists(src) and (a.redo_gifs or not os.path.exists(gif)):
                    small_gif(src, gif)
                st["gif"] = f"{c['run']}/{f}_{E}.gif" if os.path.exists(gif) else None
                row["runs"].append(st)
            th = [(x["therm"][0], x["E"]) for x in row["runs"] if x["therm"][0] is not None]
            row["best"] = max(th) if th else None
            sp = f"{SITE}/{c['run']}/flow_{f}.png"
            row["steady"] = steady_info(d, sp, f, row["neon"])
            if row["steady"]["png"]:
                row["steady"]["png"] = f"{c['run']}/flow_{f}.png"
            flows.append(row)
        figs = []
        base = os.path.basename(c["fig"])
        for suf, cap in [("_surface.png", "Extraction (%) over flow × beam energy: all extracted ions, and those leaving below 0.1 eV."),
                         ("_heatmap.png", "The same grids with values."),
                         ("_dedx.png", "Energy loss per ion vs depth (the Bragg curve) at 1–5 keV."),
                         ("_stopdepth.png", "Mean depth where ions first drop below 1 eV, vs beam energy, per flow.")]:
            if os.path.exists(c["fig"] + suf):
                dst = f"{SITE}/fig/{base}{suf}"
                if not os.path.exists(dst) or fhash(dst) != fhash(c["fig"] + suf):
                    with open(c["fig"] + suf, "rb") as fi, open(dst, "wb") as fo:
                        fo.write(fi.read())
                    os.utime(dst, (os.stat(c["fig"] + suf).st_atime, os.stat(c["fig"] + suf).st_mtime))
                figs.append(dict(src=f"fig/{base}{suf}", cap=cap))
        cfgs.append(dict(key=c["key"], label=c["label"], desc=c["desc"], nen=c["nen"], energies=c["energies"],
                         flows=flows, figs=figs, traced=sum(len(x["runs"]) for x in flows)))
    common = []
    xs = "results/figures/cross_sections.png"     # tools/plot_cross_sections.py
    if os.path.exists(xs):
        dst = f"{SITE}/fig/cross_sections.png"
        if not os.path.exists(dst) or fhash(dst) != fhash(xs):
            with open(xs, "rb") as fi, open(dst, "wb") as fo:
                fo.write(fi.read())
            os.utime(dst, (os.stat(xs).st_atime, os.stat(xs).st_mtime))
        common.append(dict(src="fig/cross_sections.png",
                         cap="Cross sections. Left: BaF⁺–Ne vs beam energy: old 1/v model with its floor vs the computed "
                             "table now used. Right: Ne–Ne through the viscosity it implies: Bird's room-temperature "
                             "parameters (old runs) vs ab initio values and the new fit (cylinder runs)."))
    ceph = sh("getfattr -n ceph.dir.rbytes --only-values /ceph/submit/data/user/h/$USER 2>/dev/null").strip()
    meta = dict(updated=datetime.datetime.now().strftime("%Y-%m-%d %H:%M"),
                running=sum(v == "RUNNING" for v in jobs.values()), queued=sum(v == "PENDING" for v in jobs.values()),

                ceph=f"{int(ceph)/1e9:.0f} GB of 1 TB" if ceph.isdigit() else "n/a")
    tmpl = open(os.path.join(os.path.dirname(os.path.abspath(__file__)), "cyl3d_dashboard_template.html")).read()
    data = json.dumps(dict(meta=meta, configs=cfgs, common=common)).replace("</", "<\\/")
    open(f"{SITE}/index.html", "w").write(tmpl.replace("/*__DATA__*/null", data))
    todo = {}
    for p in sorted(glob.glob(f"{SITE}/fig/*.png") + glob.glob(f"{SITE}/run/*.*") + glob.glob(f"{SITE}/run/*/*.*")):
        rel = os.path.relpath(p, SITE)
        if manifest.get(rel) != fhash(p):
            todo[rel] = p
    json.dump(todo, open(f"{SITE}/to_publish.json", "w"), indent=0)
    size = sum(os.path.getsize(p) for p in glob.glob(f"{SITE}/fig/*") + glob.glob(f"{SITE}/run/*.*") + glob.glob(f"{SITE}/run/*/*.*"))
    runs = ", ".join(f"{c['key']} {c['traced']}" for c in cfgs)
    print(f"wrote {SITE}/index.html; runs {runs}; {len(todo)} files to publish; site files {size/1e6:.0f} MB")


if __name__ == "__main__":
    main()
