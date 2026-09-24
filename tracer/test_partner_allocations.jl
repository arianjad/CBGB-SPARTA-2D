# Run with julia --project=tracer tracer/test_partner_allocations.jl.
# Load the production tracer without its trailing command-line entry point.
using Test
const tracer_path = joinpath(@__DIR__, "ParticleTracing.jl")
const entry = "allstats = main(args)"
src = read(tracer_path, String)
@assert length(findall(entry, src)) == 1
empty!(ARGS)
append!(ARGS, ["unused-cell.surfs", "unused-field.DAT", "-n", "1",
               "--sampler", "exact"])
include_string(@__MODULE__, replace(src, entry => "# test omits entry point"), tracer_path)

# Frozen pre-change sampler for seeded stream and allocation comparisons.
function reference_partner(w0, T)
    s = sqrt(kB*T/MASS_BUFFER_GAS)
    gmax = LinearAlgebra.norm(w0) + 8.0*s
    gmax == 0.0 && return [0.0, 0.0, 0.0]
    while true
        u = s .* Random.randn(3)
        g = sqrt((w0[1]-u[1])^2+(w0[2]-u[2])^2+(w0[3]-u[3])^2)
        if Random.rand()*gmax < g
            return u
        end
    end
end

@testset "exact partner stream parity" begin
    for seed in (1, 42, 20260923), T in (0.0, 4.0, 1000.0),
        w in ([0.0, 0.0, 0.0], [1.0, -2.0, 3.0], [1000.0, 0.0, -2000.0])
        Random.seed!(seed)
        expected = [reference_partner(w, T) for _ in 1:1000]
        continuation = Random.rand(10)
        Random.seed!(seed)
        observed = [exact_partner(w, T) for _ in 1:1000]
        @test isequal(observed, expected)
        @test isequal(Random.rand(10), continuation)
    end
end

function batch(sampler, n)
    w = [100.0, -50.0, 25.0]
    checksum = 0.0
    for _ in 1:n
        checksum += sampler(w, 4.0)[1]
    end
    checksum
end

for sampler in (reference_partner, exact_partner)
    batch(sampler, 100)
    Random.seed!(20260923)
    GC.gc()
    measured = @timed batch(sampler, 100000)
    println("sampler=$(nameof(sampler)) draws=100000 seconds=$(measured.time) bytes=$(measured.bytes) checksum=$(measured.value)")
end
