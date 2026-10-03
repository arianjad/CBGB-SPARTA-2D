# Per-run ion diagnostics for ParticleTracing3D.jl (--ionlog / --hist), so that
# stopping curves, dE/dx (Bragg-peak) profiles, energy-vs-depth bands and ion
# dwell maps can be plotted later without re-tracing.
#
#   IonDiag (one per thread chunk; merged at the end):
#     xr histogram over (x, r = sqrt(y^2 + z^2)):
#       edep  : sum over collisions of the ion's lab-frame kinetic-energy loss
#               (E_before - E_after, eV) at the collision point. Net: slow ions
#               can gain energy from the gas. dE/dx per ion = sum over r / (nions dx).
#       dwell : ion time spent (s), each free flight added at its midpoint
#       ncoll : number of collisions
#     xe histogram over (x, log10 E_before/eV): collision counts (slowing-down band)
#   Per-ion milestones (one row per ion, rows are written by their own thread):
#     x/t/collides when the lab KE first drops below 1 eV and below 0.1 eV
#     (NaN if never), max x reached, max r reached, total path length, final KE.

const AMU_EV = 1.66053907e-27 / 1.602176634e-19   # kg per amu / J per eV

struct HistSpec
    xlo::Float64
    dx::Float64
    nx::Int
    dr::Float64
    nr::Int
    elo::Float64      # log10(E/eV) lower edge
    de::Float64
    ne::Int
end

struct IonDiag
    edep::Matrix{Float64}
    dwell::Matrix{Float64}
    ncoll::Matrix{Float64}
    xe::Matrix{Float64}
end

IonDiag(h::HistSpec) = IonDiag(zeros(h.nx, h.nr), zeros(h.nx, h.nr), zeros(h.nx, h.nr), zeros(h.nx, h.ne))

const MS_COLS = [:idx, :x_1eV, :t_1eV, :n_1eV, :x_0p1eV, :t_0p1eV, :n_0p1eV, :xmax, :rmax, :pathlen, :E_final_eV]

@inline ke_ev(v, M) = 0.5 * M * AMU_EV * (v[1]^2 + v[2]^2 + v[3]^2)

@inline function hbin_x(h::HistSpec, x)
    i = 1 + floor(Int, (x - h.xlo) / h.dx)
    return (1 <= i <= h.nx) ? i : 0
end

@inline hbin_r(h::HistSpec, x) = min(h.nr, 1 + floor(Int, sqrt(x[2]^2 + x[3]^2) / h.dr))

@inline function hbin_e(h::HistSpec, E)
    E <= 0 && return 1
    return clamp(1 + floor(Int, (log10(E) - h.elo) / h.de), 1, h.ne)
end

@inline function diag_flight!(d::IonDiag, h::HistSpec, x1, x2, dt)
    xm = (0.5 * (x1[1] + x2[1]), 0.5 * (x1[2] + x2[2]), 0.5 * (x1[3] + x2[3]))
    i = hbin_x(h, xm[1]); i == 0 && return
    d.dwell[i, hbin_r(h, xm)] += dt
end

@inline function diag_collision!(d::IonDiag, h::HistSpec, x, E1, E2)
    i = hbin_x(h, x[1]); i == 0 && return
    j = hbin_r(h, x)
    d.edep[i, j] += E1 - E2
    d.ncoll[i, j] += 1
    d.xe[i, hbin_e(h, E1)] += 1
end

function merge_diag!(a::IonDiag, b::IonDiag)
    a.edep .+= b.edep; a.dwell .+= b.dwell; a.ncoll .+= b.ncoll; a.xe .+= b.xe
    return a
end

function write_diag(prefix, d::IonDiag, h::HistSpec, nions)
    open(prefix * "_xr.csv", "w") do io
        println(io, "# ParticleTracing3D.jl --hist: (x, r) histogram over all $nions ions; x bin i covers ",
                "[xlo + (i-1) dx, xlo + i dx), r bin j covers [(j-1) dr, j dr); xlo=$(h.xlo) m dx=$(h.dx) m dr=$(h.dr) m ",
                "nx=$(h.nx) nr=$(h.nr) (last r bin open-ended). edep_eV = summed lab KE loss at collisions; ",
                "dwell_s = summed ion time; ncoll = collisions. Only non-empty bins are listed.")
        println(io, "ix,ir,x_m,r_m,edep_eV,dwell_s,ncoll")
        for i in 1:h.nx, j in 1:h.nr
            (d.ncoll[i, j] == 0 && d.dwell[i, j] == 0) && continue
            @printf(io, "%d,%d,%.6e,%.6e,%.6e,%.6e,%d\n", i, j, h.xlo + (i - 0.5) * h.dx, (j - 0.5) * h.dr,
                    d.edep[i, j], d.dwell[i, j], d.ncoll[i, j])
        end
    end
    open(prefix * "_xe.csv", "w") do io
        println(io, "# ParticleTracing3D.jl --hist: collisions binned by x and log10(lab KE before the collision / eV) ",
                "over all $nions ions; xlo=$(h.xlo) m dx=$(h.dx) m; log10E bin k covers [elo + (k-1) de, elo + k de), ",
                "elo=$(h.elo) de=$(h.de) ne=$(h.ne) (end bins open-ended). Only non-empty bins are listed.")
        println(io, "ix,ie,x_m,log10E_eV,ncoll")
        for i in 1:h.nx, k in 1:h.ne
            d.xe[i, k] == 0 && continue
            @printf(io, "%d,%d,%.6e,%.4f,%d\n", i, k, h.xlo + (i - 0.5) * h.dx, h.elo + (k - 0.5) * h.de, d.xe[i, k])
        end
    end
end

function write_ionlog(path, ms::Matrix{Float64})
    open(path, "w") do io
        println(io, "# ParticleTracing3D.jl --ionlog: per-ion milestones (lab-frame KE). x_* in m, t_* in s, ",
                "n_* = collisions so far; NaN = never reached. xmax/rmax = deepest x and largest r reached, ",
                "pathlen = total flight path (m), E_final_eV = KE at the end (exit or wall).")
        println(io, join(string.(MS_COLS), ","))
        for i in 1:size(ms, 1)
            println(io, join([j == 1 ? string(Int(ms[i, j])) : @sprintf("%.6e", ms[i, j]) for j in 1:size(ms, 2)], ","))
        end
    end
end
