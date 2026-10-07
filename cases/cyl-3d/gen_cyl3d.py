"""
Write cell_cyl3d.surf: a plain cylindrical buffer-gas cell for SPARTA 3D.

Same chamber size as fluor-cell-3d, without the small features. The body is
the solid of revolution of this (x, r) profile about the x axis (y = z = 0):

    inlet plate   x in [-T, 0],      r in [A, RO]   (aperture radius A)
    side wall     x in [0, L],       r in [R, RO]
    outlet plate  x in [L, L + T],   r in [A, RO]

so the gas chamber is x in [0, L], r < R, with an aperture through each end
plate on the axis. Ions enter at the inlet face (x = -T) and leave through
the outlet face (x = L + T).

Sizes from the fluor-cell-3d chamber (ray-cast from cell_fluor3d_sealed.surf
on the port axis): 30.48 mm bore (r 15.24 mm), 45.47 mm between the inlet
channel end (x = 31.75 mm) and the outlet aperture plane (x = 77.22 mm).
End plates 1.27 mm thick (the old inlet channel length).

Neon enters through a flush patch on the side wall (type 2, emit/surf) near
the inlet end at -z, 7 mm (axial) x 6 mm (arc) = 42 mm^2. That is about the
two apertures' total area (39 mm^2), so the injected gas is no denser than
the chamber gas. The old 2 mm^2 fill line would be ~5x denser than the
chamber at 50 SCCM and would set the cell size and time step of the whole run.

Surface types: 1 = wall, 2 = buffer-gas inlet patch, 3 = inlet aperture bore,
4 = outlet aperture bore. Normals point out of the body into the gas
(positive signed volume), same as cell_fluor3d_sealed.surf.

--config:
  ap5    (default before 2026-10-03; cell_cyl3d.surf): 5 mm apertures through 1.27 mm plates.
  ports  (cell_cyl3d_ports.surf): the fluor-cell-3d ports, 2.54 mm inlet through a
         1.27 mm wall and 3.175 mm outlet through a 5.334 mm channel (old channel
         lengths, x 30.48-31.75 and 77.216-82.55 mm). Chamber and gas patch
         unchanged; outlet face moves to x = 50.80 mm.
  ports-plug  ports with the outlet plugged.

Usage (from this directory): python gen_cyl3d.py [--config ap5|ports|ports-plug]
"""
import argparse
import math
import os

HERE = os.path.dirname(os.path.abspath(__file__))
OUT = os.path.join(HERE, "cell_cyl3d.surf")

MM = 1e-3
R = 15.24 * MM        # chamber bore radius
L = 45.47 * MM        # chamber length
CONFIGS = {  # inlet radius, inlet wall, outlet radius, outlet wall (m), output file
    "ap5": (2.5 * MM, 1.27 * MM, 2.5 * MM, 1.27 * MM, "cell_cyl3d.surf"),
    "ports": (1.27 * MM, 1.27 * MM, 1.5875 * MM, 5.334 * MM, "cell_cyl3d_ports.surf"),
    # outlet plugged: far end is solid (outlet radius 0, the 5.334 mm channel stays
    # solid body). Gas and ions leave only through the inlet port.
    "ports-plug": (1.27 * MM, 1.27 * MM, 0.0, 5.334 * MM, "cell_cyl3d_ports_plug.surf"),
}
A, T, A2, T2 = CONFIGS["ap5"][:4]
RO = 19.05 * MM       # outer radius of the body

NTH = 192             # azimuthal segments (0.50 mm chords on the bore, 2 um sag)
GAS_X = 4.04 * MM     # gas-patch centre, downstream of the inlet wall (patch x 0.54..7.54 mm)
GAS_DX = 7.0 * MM     # gas-patch axial length
GAS_NSEG = 12         # gas-patch azimuthal segments (12 x 0.50 mm = 6 mm), centred on -z

# profile vertices (x, r), walked so the revolved triangles face out of the body
def make_profile():
    if A2 == 0.0:   # plugged outlet: the end wall runs down to the axis
        return [(-T, A), (0.0, A), (0.0, R),
                (GAS_X - GAS_DX / 2, R), (GAS_X + GAS_DX / 2, R),
                (L, R), (L, 0.0), (L + T2, 0.0), (L + T2, RO), (-T, RO)]
    return [(-T, A), (0.0, A), (0.0, R),
            (GAS_X - GAS_DX / 2, R), (GAS_X + GAS_DX / 2, R),
            (L, R), (L, A2), (L + T2, A2), (L + T2, RO), (-T, RO)]


def edge_type(p, q, j):
    """Surface type of the quad swept by profile edge p->q at azimuthal segment j."""
    (x1, r1), (x2, r2) = p, q
    if r1 == r2 == A and x2 <= 0:
        return 3
    if A2 > 0 and r1 == r2 == A2 and x1 >= L:
        return 4
    if r1 == r2 == R and x1 == GAS_X - GAS_DX / 2:
        # segments centred on theta = -90 deg (-z)
        j0 = 3 * NTH // 4 - GAS_NSEG // 2
        return 2 if j0 <= j < j0 + GAS_NSEG else 1
    return 1


def main():
    global A, T, A2, T2, OUT, PROFILE
    ap = argparse.ArgumentParser()
    ap.add_argument("--config", choices=list(CONFIGS), default="ap5")
    cfg = ap.parse_args().config
    A, T, A2, T2, name = CONFIGS[cfg]
    OUT = os.path.join(HERE, name)
    PROFILE = make_profile()
    pts, index = [], {}

    def pid(k, j):
        key = (k, 0 if PROFILE[k][1] == 0.0 else j % NTH)   # points on the axis are shared
        if key not in index:
            x, r = PROFILE[k]
            th = 2 * math.pi * (j % NTH) / NTH
            pts.append((x, r * math.cos(th), r * math.sin(th)))
            index[key] = len(pts)
        return index[key]

    tris = []
    n = len(PROFILE)
    for k in range(n):
        k2 = (k + 1) % n
        for j in range(NTH):
            a, b = pid(k, j), pid(k2, j)
            c, d = pid(k2, j + 1), pid(k, j + 1)
            t = edge_type(PROFILE[k], PROFILE[k2], j)
            for tri in ((t, a, b, c), (t, a, c, d)):
                if len({tri[1], tri[2], tri[3]}) == 3:   # drop triangles collapsed on the axis
                    tris.append(tri)

    # orientation: flip everything if the signed volume came out negative
    vol = 0.0
    for _, a, b, c in tris:
        p, q, s = pts[a - 1], pts[b - 1], pts[c - 1]
        vol += (p[0] * (q[1] * s[2] - q[2] * s[1]) - p[1] * (q[0] * s[2] - q[2] * s[0])
                + p[2] * (q[0] * s[1] - q[1] * s[0])) / 6
    if vol < 0:
        tris = [(t, a, c, b) for t, a, b, c in tris]
        vol = -vol
    exact = math.pi * ((RO**2 - A**2) * T + (RO**2 - A2**2) * T2 + (RO**2 - R**2) * L)

    with open(OUT, "w") as f:
        f.write("# GENERATED by gen_cyl3d.py: plain cylindrical buffer-gas cell, axis = x at y = z = 0.\n")
        f.write(f"# Config {cfg}: chamber x 0..{L/MM:.2f} mm, bore r {R/MM:.2f} mm; inlet {2*A/MM:.3f} mm through "
                f"{T/MM:.3f} mm (face x = {-T/MM:.2f} mm), outlet {2*A2/MM:.3f} mm through {T2/MM:.3f} mm "
                f"(face x = {(L+T2)/MM:.2f} mm).\n")
        f.write("# type 1 = wall, 2 = buffer gas inlet patch (emit/surf), 3 = inlet aperture bore, "
                "4 = outlet aperture bore.\n# Units: meters (SPARTA si).\n\n")
        f.write(f"{len(pts)} points\n{len(tris)} triangles\n\nPoints\n\n")
        for i, (x, y, z) in enumerate(pts, 1):
            f.write(f"{i} {x:.9g} {y:.9g} {z:.9g}\n")
        f.write("\nTriangles\n\n")
        for i, (t, a, b, c) in enumerate(tris, 1):
            f.write(f"{i} {t} {a} {b} {c}\n")
    ntype = {t: sum(1 for q in tris if q[0] == t) for t in (1, 2, 3, 4)}
    print(f"wrote {OUT}: {len(pts)} points, {len(tris)} triangles, types {ntype}")
    print(f"body volume {vol*1e9:.1f} mm^3 (exact solid of revolution {exact*1e9:.1f})")
    print(f"chamber volume {math.pi*R**2*L*1e9:.0f} mm^3")


if __name__ == "__main__":
    main()
