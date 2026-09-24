# 3D performance lessons checked against the 2D workflow

The 3D tracer's allocation-reduced exact collision-partner sampler has the
same source algorithm as this 2D tracer. This branch draws its three Gaussian
components as scalars, allocating a three-element vector only after a proposal
is accepted. It preserves the rejection expression, draw order, legacy table
construction, and collision and flight models. In particular, it does not
remove the inherited finite `gmax` proposal envelope or change the approximate
free-flight hazard. This is an implementation-equivalence and efficiency
change, not a physics validation.

## Checks and measurements

Run the data-independent regression from a fresh public clone with Julia 1.9.4
and the pinned `tracer/Manifest.toml` installed:

```sh
julia --startup-file=no --project=tracer --threads=1 tracer/test_partner_allocations.jl
```

It compares the production sampler against the pre-change body for 27
seed/temperature/slip cases, 1000 accepted draws per case, and ten subsequent
RNG draws. All 54 assertions passed. A single warmed 100,000-draw diagnostic
within that script allocated 74,463,120 bytes with the reference and 8,000,080
bytes with the candidate; checksums matched. These are cumulative allocations,
not peak memory.

For a public input path, `tools/field2tracer.py --synthetic 4.0 2.0e21 0 0.2 0
0.02 40 10 --out NEW_FIELD` produced a 400-cell, 4 K field. A baseline at
`43cd904` and this branch each traced 200 molecules with the production
`crossing.jl` driver at the same seed and thread count. Both `legs.out` and
`recs.csv` matched byte for byte within each one-thread and four-thread pair.
The run used seed 20260923, 156.325 u, `2.7e-18 m²`, a 300 K point source at
`z=0.05 m`, the exact sampler, and an observation plane at `z=0.065 m`.
The driver's leg-conservation and endpoint-identity checks passed. The paired
one-thread runs had 2,115,875 legs; the paired four-thread runs had 2,060,635.
Thread counts have different RNG partitioning and are not cross-compared.

A separate local benchmark used that generated field, the actual `Crossing`
accumulator, a 300 K point source, 200 particles, seed 20260923, and one native
Julia 1.9.4 thread. It used the tracer's default 191 u and `1.3e-18 m²`,
unlike the crossing-output check above. Each process warmed 10 particles, reset the seed, then
timed one `SimulateParticles` call including field/table construction and
transport. Three serial interleaved baseline/candidate repetitions used order
B,C,C,B,B,C. All six calls had 329,757 collision legs. Startup, warmup,
formatting and output writing were outside the timed call.

| Full call | Baseline median [range] | Candidate median [range] |
|---|---:|---:|
| Time (s) | 0.306 [0.305, 0.324] | 0.239 [0.226, 0.297] |
| Allocated bytes (each repeat) | 406,764,648 | 164,883,848 |

The 59.5% allocation reduction is exact for these fixed-workload repetitions.
The median time was 21.8% lower in this campaign; three repetitions on one
shared host do not establish a general runtime gain, thread scaling, or a peak
RAM reduction. Other fields, sources, particle counts and Julia versions may
give different costs. The benchmark harness and generated data are local
ignored outputs; the public regression above is independently runnable.

## Other transfer decisions

| 3D investigation | 2D decision |
|---|---|
| Compact text-output batching | Defer. The 3D output-only replay found no clear practical gain and allocated more; the 2D crossing writer uses the same small end-of-run text pattern. |
| Whole-field thermal-width cache | Defer. The 3D probe retained 16.55 MB for only a small isolated kernel difference and established no full-tracer gain. The 2D field is sampled cell by cell without such a cache. |
| Vector particle batching | Defer. The 3D drift replay showed poor useful-slot occupancy for unequal lifetimes and no actual-path gain. A 2D rewrite would also need independent particle RNG and termination parity. |
| Column-oriented field reader | Defer pending a 2D whole-call benchmark. The 3D result required a concrete return-type fix before its isolated reader gain reached the full call; the present 2D field build and file size differ. |
| Scalar exact collision partner | Applied here. Seed and whole-trajectory parity pass at one and four threads; whole-call allocation falls for the specified synthetic workload. |
| Run-scoped page-cache control | The helium solver launcher is a separate concern. The 3D runner records an optional `systemd-run` memory scope for large field dumps; this branch does not change the 2D helium launcher or its run semantics. |

The unused legacy lookup table is intentionally still built. In 3D, omitting
it changed the calling task's RNG continuation and seeded trajectories until
field construction was given an explicit isolated-RNG protocol. Reusing a gas
field across several molecule trials is another possible 3D workflow transfer,
but it needs an explicit 2D trial contract and fresh-versus-shared equivalence
check before implementation.
