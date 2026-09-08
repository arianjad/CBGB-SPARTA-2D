# Third-party notices

The root MIT license covers original project material. It does not replace
the licenses below.

## Molecule tracer

`tracer/ParticleTracing.jl` derives from `ParticleTracing/ParticleTracing.jl`
in Cal Miller's [DSMC_Simulations](https://github.com/cal-miller-harvard/DSMC_Simulations)
at commit `543b42e4a469c3275e456742c21c593e4a82623b`. The tracer and code derived
from it retain GNU GPL version 3 terms. The upstream license is included in
`licenses/GPL-3.0.txt`. The tracer header summarizes its current local
modifications.

## SPARTA

SPARTA is developed by Steve Plimpton, Michael Gallis, and contributors at
Sandia National Laboratories. The solver source and binaries are downloaded
from [sparta/sparta](https://github.com/sparta/sparta), not redistributed here.
The pinned build recipe uses tag `27Aug2026`, commit
`95b9abaa8bd548991cc3c3f1c58b34722f7ade74`. SPARTA retains its GNU GPL version 2
terms; the upstream license is included in `licenses/GPL-2.0.txt`.

## Julia and Python dependencies

Julia and package dependencies are downloaded separately and retain their
own upstream licenses. The optional Julia installer verifies the archive
against the checksum published on the official
[older releases page](https://julialang.org/downloads/oldreleases/).
