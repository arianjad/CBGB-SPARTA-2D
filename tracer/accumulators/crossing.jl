# 2D first-crossing accumulator for one observation plane.
#
# The vendored tracer is NOT modified.  This driver
#   (a) includes ParticleTracing.jl with the run's real flags but -n 1, purely
#       to set its const globals (MASS_PARTICLE, sigma, SPAWNCLIP, SPAWN_REF,
#       TRAJPRINT, EXACT_SAMPLER) and to build `args`;
#   (b) re-seeds and calls SimulateParticles directly with make_stats=Crossing,
#       replicating main()'s generateParticle closure verbatim
#       so the RNG stream matches a direct tracer run at the same --seed and
#       --threads.
#
# Why a union and not the hook alone: the terminating leg never reaches
# updateStats! (propagate returns at the `getCollision != 0` branch before the
# hook), and downstream of the aperture that is exactly the leg most extracted
# molecules cross the plane on.  So legs 1..N-1 come from the hook and leg N
# comes from the `outputs` row, clipped in tools/tracer_common.py the same way
# tracer_analyze.clip_legs does.
#
# Usage (from the repo root):
#   julia --project=tracer --threads=8 tracer/accumulators/crossing.jl \
#       OUTDIR XOBS <every ParticleTracing.jl argument, verbatim>
#
# Writes OUTDIR/legs.out  (idx-sorted rows in the tracer's own stdout format)
#    and OUTDIR/recs.csv  (one row per particle, hook-side crossing state).
# Endpoint-identity and leg-conservation checks run on every invocation.

length(ARGS) >= 3 || error("usage: crossing.jl OUTDIR XOBS <tracer args...>")
const OUTDIR = ARGS[1]
const XOBS = parse(Float64, ARGS[2])
const REST = ARGS[3:end]

# nParticles comes from the run's own -n; the include below uses -n 1.
let k = findfirst(==("-n"), REST)
    isnothing(k) && error("no -n in the tracer argument list")
    global NPART = parse(Int, REST[k + 1])
    global INCARGS = copy(REST)
    INCARGS[k + 1] = "1"
end

empty!(ARGS); append!(ARGS, INCARGS)
# The included one-particle warmup initializes the vendored globals. Suppress
# its stdout rows so a requested --trajprint stream contains only this run.
mkpath(OUTDIR)
redirect_stdout(devnull) do
    include(joinpath(@__DIR__, "..", "ParticleTracing.jl"))
end

# ---------------------------------------------------------------- accumulator
#
# One instance per particle (SimulateParticles calls new_stats() per particle),
# plus one per chunk partial and one total.  A particle instance has an empty
# `recs`; a container has at least one record by the time it is merged upward,
# which is what merge! keys on.  Records land in chunk order, and chunks
# partition 1:n contiguously in order, so recs[k] is particle k; this is
# checked against the per-particle collision count below.
mutable struct Crossing
    recs::Vector{NTuple{12,Float64}}
    nlegs::Int
    nfwd::Int
    nback::Int
    have::Bool
    cx::Float64
    cy::Float64
    vx::Float64
    vy::Float64
    vz::Float64
    tc::Float64
    nc::Float64
    px::Float64
    py::Float64
    pz::Float64
    hasprev::Bool
    mism::Float64
end
Crossing() = Crossing(NTuple{12,Float64}[], 0, 0, 0, false,
                      0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
                      0.0, 0.0, 0.0, false, 0.0)

# x is the leg's START (the vendored file assigns x .= xnext on the next line),
# v the leg's velocity, lfree its free path, t and ncolls cumulative AFTER it.
@inline function updateStats!(a::Crossing, x, v, t, ncolls, lfree)
    vmag = sqrt(v[1]^2 + v[2]^2 + v[3]^2)
    dt = lfree / vmag
    ex = x[1] + v[1] * dt
    ey = x[2] + v[2] * dt
    ez = x[3] + v[3] * dt
    if a.hasprev                                   # endpoint identity
        m = max(abs(x[1] - a.px), abs(x[2] - a.py), abs(x[3] - a.pz))
        m > a.mism && (a.mism = m)
    end
    a.px, a.py, a.pz, a.hasprev = ex, ey, ez, true
    a.nlegs += 1
    if (x[3] - XOBS) * (ez - XOBS) < 0             # strict straddle of the plane
        if ez > x[3]
            a.nfwd += 1
            if !a.have
                f = (XOBS - x[3]) / (ez - x[3])
                a.have = true
                a.cx = x[1] + f * (ex - x[1])
                a.cy = x[2] + f * (ey - x[2])
                a.vx, a.vy, a.vz = v[1], v[2], v[3]
                a.tc = t - dt + f * dt
                a.nc = ncolls - 1                  # collisions before the crossing
            end
        elseif a.have
            a.nback += 1                           # back-crossing after the first
        end
    end
    return a
end

merge!(a::Crossing, b::Crossing) = begin
    if isempty(b.recs)
        push!(a.recs, (Float64(b.nlegs), Float64(b.nfwd), Float64(b.nback),
                       b.have ? 1.0 : 0.0, b.cx, b.cy,
                       b.vx, b.vy, b.vz, b.tc, b.nc, b.mism))
    else
        append!(a.recs, b.recs)
    end
    a
end

# ------------------------------------------------------------------- the run
# main()'s particle generator, copied so the RNG draws match a direct run.
const boltzmann = sqrt(kB * args["T"] / MASS_PARTICLE)
const spawnmode = args["spawn"]
spawnPosition() =
    if spawnmode == "point"
        [args["r"], 0.0, args["z"]]
    elseif spawnmode == "gaussball"
        s = args["spawnsize"]
        [args["r"] + s*Random.randn(), s*Random.randn(), args["z"] + s*Random.randn()]
    elseif spawnmode == "uniformball"
        u = Random.randn(3)
        u .*= args["spawnsize"] * cbrt(Random.rand()) / LinearAlgebra.norm(u)
        [args["r"] + u[1], u[2], args["z"] + u[3]]
    else
        rlo = args["spawnrlo"]; rhi = args["spawnrhi"]
        zlo = args["spawnzlo"]; zhi = args["spawnzhi"]
        r = sqrt(rlo^2 + Random.rand()*(rhi^2 - rlo^2))
        z = zlo + Random.rand()*(zhi - zlo)
        if spawnmode == "uniformplane"
            [r, 0.0, z]
        else
            ϕ = 2π*Random.rand()
            [r*cos(ϕ), r*sin(ϕ), z]
        end
    end
generateParticle() = (
    spawnPosition(),
    [args["vr"] + Random.randn() * boltzmann, Random.randn() * boltzmann,
     args["vz"] + Random.randn() * boltzmann])

args["seed"] != 0 && Random.seed!(args["seed"])
t0 = time()
outputs, _, acc = SimulateParticles(
    args["geom"], args["flow"], NPART, generateParticle,
    false,                       # print_stuff: rows are written below instead
    args["omega"], args["zmin"], args["zmax"], args["pflip"],
    args["saveall"],
    true, false;                 # exactly what main() passes: !isnothing(args[...])
    make_stats = Crossing, savespawns = args["spawnout"])
wall = time() - t0
@printf(stderr, "crossing.jl: %d particles, %.1f s wall, x_obs = %.6f m\n",
        NPART, wall, XOBS)

# ------------------------------------------------------------------- outputs
open(joinpath(OUTDIR, "legs.out"), "w") do io
    println(io, "idx x y z xnext ynext znext vx vy vz collides time")
    for i in 1:NPART
        print(io, @sprintf("%d %e %e %e %e %e %e %e %e %e %d %e\n", i,
            outputs[i,1], outputs[i,2], outputs[i,3], outputs[i,4], outputs[i,5],
            outputs[i,6], outputs[i,7], outputs[i,8], outputs[i,9],
            outputs[i,10], outputs[i,11]))
    end
end
open(joinpath(OUTDIR, "recs.csv"), "w") do io
    println(io, "idx,nlegs,nfwd,nback,have,cx,cy,vx,vy,vz,tc,nc,mism")
    for (k, r) in enumerate(acc.recs)
        print(io, @sprintf("%d,%d,%d,%d,%d,%.17g,%.17g,%.17g,%.17g,%.17g,%.17g,%d,%.17g\n",
            k, Int(r[1]), Int(r[2]), Int(r[3]), Int(r[4]), r[5], r[6],
            r[7], r[8], r[9], r[10], Int(r[11]), r[12]))
    end
end

# -------------------------------------------------------- consistency checks
ok = true
if length(acc.recs) != NPART
    @printf(stderr, "FAIL record count: %d records for %d particles\n", length(acc.recs), NPART)
    ok = false
else
    # count(), not a for-loop: a `bad += 1` inside a top-level loop lands in a
    # new local each iteration (Julia soft scope) and the gate would always pass.
    bad = count(i -> Int(acc.recs[i][1]) != Int(outputs[i,10]), 1:NPART)
    if bad != 0
        @printf(stderr, "FAIL leg conservation: %d of %d particles have hook legs != collides\n", bad, NPART)
        ok = false
    else
        @printf(stderr, "PASS leg conservation: hook legs == collides for all %d particles (sum %d)\n",
                NPART, Int(sum(outputs[:,10])))
    end
end
mism = isempty(acc.recs) ? 0.0 : maximum(r[12] for r in acc.recs)
if mism > 1e-12
    @printf(stderr, "FAIL endpoint identity: max mismatch %.3e m exceeds 1e-12\n", mism)
    ok = false
else
    @printf(stderr, "PASS endpoint identity: max mismatch %.3e m\n", mism)
end
nfirst = count(r -> r[4] != 0.0, acc.recs)
nback = count(r -> r[3] != 0.0, acc.recs)
@printf(stderr, "hook-side: first forward crossings %d, particles with a back-crossing %d\n",
        nfirst, nback)
exit(ok ? 0 : 1)
