# 3D performance transfer to the 2D tracer

The 3D tracer measured a large allocation reduction when its exact collision
partner sampler draws three scalar Gaussian components and allocates a vector
only for an accepted proposal. The 2D tracer uses the same sampler body. This
branch tests that local transfer with the 2D code, holding the rejection rule,
floating-point expression order, and random stream fixed.

The 3D field-reader change and opt-in legacy-table omission are separate ideas.
The former needs a 2D field equivalence and whole-call benchmark; the latter
changes seed-to-trajectory mapping unless field construction has its own RNG
protocol. Neither is assumed to transfer merely from the 3D result.

The isolated 2D worktree starts at public `main` commit `43cd904`. A 200-particle
baseline run on the existing 400-cell synthetic field from the private 3D
checkout passed crossing leg conservation and endpoint identity (34,180 legs).
This field is a local validation input, not included in the public repository.

The sampler change passes 54 checks over 27 seed/temperature/slip cases:
1000 accepted partners per case match the frozen reference exactly, and the
next ten RNG draws match. On native Julia 1.9.4, one thread, a warmed
100,000-accepted-partner batch allocated 74,463,120 bytes in the reference
and 8,000,080 bytes in the candidate. The two checksums matched. This is a
single kernel diagnostic; it is not a whole-tracer timing claim. Allocation
counts mean cumulative allocation traffic, not peak RAM.

A matching full 2D trajectory check and scoped timing comparison remain to be
completed before submitting the branch.
