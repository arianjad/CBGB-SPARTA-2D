"""VSS fit for helium-4 at 2-10 K (cases/cyl-3d/he.vss).

Data: dilute-gas viscosity of He-4 from the NIST Chemistry WebBook, Thermophysical
Properties of Fluid Systems (isobar p = 1e-4 bar, 2.18-30 K), fetched 2026-10-05.
Least squares in ln eta, uniform in ln T over 2.2-10 K, T_ref = 4 K, alpha = 1.26 (Bird),
then d rescaled so eta(4 K) is exact. Run: python tools/cross_section/he_vss_fit.py
"""
import numpy as np
k = 1.380649e-23; m = 4.002602 * 1.66053906660e-27
T = np.array([2.1768, 3.1768, 4.1768, 5.1768, 6.1768, 7.1768, 8.1768, 9.1768, 10.177, 11.177, 12.177, 13.177, 14.177,
              15.177, 16.177, 17.177, 18.177, 19.177, 20.177, 21.177, 22.177, 23.177, 24.177, 25.177, 26.177, 27.177,
              28.177, 29.177, 30.177])
eta = 1e-6 * np.array([0.53840, 0.81707, 1.0640, 1.2888, 1.4966, 1.6906, 1.8731, 2.0459, 2.2104, 2.3676, 2.5185, 2.6637,
                       2.8038, 2.9394, 3.0709, 3.1986, 3.3229, 3.4441, 3.5623, 3.6779, 3.7909, 3.9016, 4.0102, 4.1167,
                       4.2213, 4.3241, 4.4253, 4.5248, 4.6228])
lt, le = np.log(T), np.log(eta)
a, Tr = 1.26, 4.0
pref = lambda w: 5 * (a + 1) * (a + 2) * np.sqrt(np.pi * m * k * Tr) / (4 * a * np.pi * (5 - 2 * w) * (7 - 2 * w))
tt = np.exp(np.linspace(np.log(2.2), np.log(10), 300)); y = np.interp(np.log(tt), lt, le)
c, *_ = np.linalg.lstsq(np.vstack([np.ones_like(tt), np.log(tt / Tr)]).T, y, rcond=None)
w = round(c[1], 2)
e4 = np.exp(np.interp(np.log(Tr), lt, le))
d = np.sqrt(pref(w) / e4)
ev = lambda t: pref(w) / d**2 * (t / Tr)**w
print(f"He  {d:.4g}  {w:.2f}  {Tr}  {a}")
for t in [2.2, 3, 4, 6, 8, 10, 15, 20]:
    print(f"  {t:5.1f} K: VSS {ev(t)*1e6:.3f}  NIST {np.exp(np.interp(np.log(t), lt, le))*1e6:.3f} uPa s  "
          f"({(ev(t)/np.exp(np.interp(np.log(t), lt, le)) - 1)*100:+.1f}%)")
