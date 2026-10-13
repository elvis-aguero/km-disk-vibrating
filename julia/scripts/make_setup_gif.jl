#!/usr/bin/env julia
#= Renders docs/setup.gif: a minimal animation of the problem this
solver simulates, driven by a real (small-domain) simulation rather than a cartoon — the
surface waves, disk motion and phase lag shown are what `run_simulation` computes.

Usage: `julia -t auto --project=julia/scripts julia/scripts/make_setup_gif.jl`
(about a minute; writes docs/setup.gif at the repository root).

Block comment, not a docstring: see _bootstrap.jl's header. =#

if !@isdefined(FaradayDisk)
    include("_bootstrap.jl")
end
using CairoMakie, Logging

global_logger(Logging.NullLogger())

# ---------------- simulation (sweep fluid, small domain) ----------------
const Γ, F = 0.5, 30.0
params = SimulationParams(bathFrequency = F, bathAmplitude = Γ * 981 / (2π * F)^2, simulationTime = 30 / F,
                          bathDensity = 1.175, bathSurfaceTension = 66.5, bathViscosity = 0.18 / 1.175,
                          temporalResolution = 60, bathDiameter = 30.0, spatialResolution = 10.0)
nr = ceil(Int, params.bathDiameter * params.spatialResolution / 2)
result = run_simulation(params; dtn = generate_dtn(nr, params.bathDiameter))

R = params.diskRadius
A = params.bathAmplitude
ω = 2π * F
φ = deg2rad(params.phaseDifference)
spc = params.temporalResolution
idx = (length(result.t_s) - 2spc):length(result.t_s)        # last two periods (converged)
t = result.t_s[idx]
z_disk = lab_frame_CoM(result)[idx]
z_bath = A .* cos.(ω .* t .+ φ)
r_cm = (0:(nr - 1)) .* (params.bathDiameter / (2nr)) .* R
η = result.eta_history_cm[:, idx]                             # relative to the bath

# ---------------- drawing: symbols only ----------------
const VEX = 20.0                                              # vertical exaggeration
const BG, INK = colorant"#0d2338", colorant"#d6ebff"
const CYAN, AMBER, DIM = colorant"#4fd1ff", colorant"#ffb347", colorant"#5f7d99"
xmax = params.bathDiameter / 2 * R                            # domain edge, cm
depth = 1.1                                                   # drawn fluid depth (display units)

fig = Figure(size = (900, 420), backgroundcolor = BG, figure_padding = 10)
ax = Axis(fig[1, 1], backgroundcolor = BG, limits = (-xmax - 0.5, xmax + 0.2, -1.9, 1.12))
hidedecorations!(ax); hidespines!(ax)

k = Observable(1)
zb = @lift VEX * z_bath[$k]
surf = @lift begin
    y = VEX .* (η[:, $k] .+ z_bath[$k])
    Point2f.(vcat(-reverse(r_cm), r_cm), vcat(reverse(y), y))
end
zd = @lift VEX * z_disk[$k]
label(pos, s; color = INK, size = 22, align = (:left, :center)) =
    text!(ax, pos; text = s, color, fontsize = size, align)
label(x::Real, y::Real, s; kw...) = label(Point2f(x, y), s; kw...)

# fluid
poly!(ax, @lift(vcat($surf, [Point2f(xmax, $zb - depth), Point2f(-xmax, $zb - depth)])), color = (CYAN, 0.13))
lines!(ax, surf, color = CYAN, linewidth = 2)
# container: walls and floor
lines!(ax, @lift(Point2f[(-xmax, $zb + 0.3), (-xmax, $zb - depth), (xmax, $zb - depth), (xmax, $zb + 0.3)]), color = DIM, linewidth = 1.5)
label(@lift(Point2f(1.9, VEX * (η[3nr ÷ 5, $k] + z_bath[$k]) + 0.17)), L"\eta"; color = CYAN)
label(@lift(Point2f(-xmax + 0.15, $zb - 0.7)), L"\rho \quad \sigma \quad \nu"; color = DIM, size = 19)

# disk
poly!(ax, @lift(Rect2f(-R, $zd, 2R, 0.14)), color = AMBER)
lines!(ax, @lift(Point2f[(-R, $zd - 0.1), (R, $zd - 0.1)]), color = AMBER)
scatter!(ax, @lift(Point2f[(-R, $zd - 0.1), (R, $zd - 0.1)]), marker = [:ltriangle, :rtriangle], color = AMBER, markersize = 8)
label(@lift(Point2f(0, $zd - 0.16)), L"2R"; color = AMBER, align = (:center, :top))
label(@lift(Point2f(R + 0.08, $zd + 0.07)), L"m"; color = AMBER)
# direct forcing of the disk itself (forceAmplitude, forceFrequency)
lines!(ax, @lift(Point2f[(0, $zd + 0.22), (0, $zd + 0.62)]), color = AMBER)
scatter!(ax, @lift(Point2f[(0, $zd + 0.22), (0, $zd + 0.62)]), marker = [:dtriangle, :utriangle], color = AMBER, markersize = 9)
label(@lift(Point2f(0, $zd + 0.7)), L"F\,\cos\,\Omega t"; color = AMBER, align = (:center, :bottom))

# domain
lines!(ax, [Point2f(-xmax, 0.8), Point2f(xmax, 0.8)], color = DIM)
scatter!(ax, [Point2f(-xmax, 0.8), Point2f(xmax, 0.8)], marker = [:ltriangle, :rtriangle], color = DIM, markersize = 8)
label(0, 0.88, L"D"; color = DIM, align = (:center, :bottom))

# gravity and shaking
lines!(ax, [Point2f(-xmax - 0.3, 0.6), Point2f(-xmax - 0.3, 0.2)], color = INK)
scatter!(ax, [Point2f(-xmax - 0.3, 0.2)], marker = :dtriangle, color = INK, markersize = 9)
label(-xmax - 0.22, 0.42, L"g")
lines!(ax, @lift(Point2f[(0, $zb - depth - 0.12), (0, $zb - depth - 0.5)]), color = INK)
scatter!(ax, @lift(Point2f[(0, $zb - depth - 0.12), (0, $zb - depth - 0.5)]), marker = [:utriangle, :dtriangle], color = INK, markersize = 9)
label(@lift(Point2f(0.12, $zb - depth - 0.31)), L"A\,\cos(\omega t + \varphi)")

outpath = joinpath(@__DIR__, "..", "..", "docs", "setup.gif")
mkpath(dirname(outpath))
record(fig, outpath, eachindex(t); framerate = 20) do i   # every step: 3 s per period
    k[] = i
end
println("wrote ", abspath(outpath))
