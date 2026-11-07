#= Exact-solution test of the DtN matrix with an ALGEBRAIC tail (not a Gaussian packet): a point source at height a above the surface gives
phi(r,z) = (a - z) / (r^2 + (a - z)^2)^(3/2), harmonic for z < a and decaying like r^-3. On the surface f(r) = a / (r^2 + a^2)^(3/2) and
(N f)(r) = phi_z(r, 0) = (2 a^2 - r^2) / (r^2 + a^2)^(5/2) (exact, deep water). Truncating f at the last node changes N f by about R^-4.
Error columns: max over nodes r <= R/2 of |A f - exact| / max|exact|  (global) and the same divided by |exact| where |exact| > 1e-3 max (local).
Usage: julia julia/scripts/dtn_pointsource_check.jl =#
include(joinpath(@__DIR__, "_bootstrap.jl"))
using Printf
f(r, a) = a / (r^2 + a^2)^(3 / 2)
Nf(r, a) = (2a^2 - r^2) / (r^2 + a^2)^(5 / 2)
function report(label, A, r, a)
    v = f.(r, a); ex = Nf.(r, a); got = A * v
    keep = r .<= maximum(r) / 2
    glob = maximum(abs.(got[keep] .- ex[keep])) / maximum(abs.(ex[keep]))
    big = keep .& (abs.(ex) .> 1e-3 * maximum(abs.(ex[keep])))
    loc = maximum(abs.(got[big] .- ex[big]) ./ abs.(ex[big]))
    @printf("%-44s a = %-6g global %.2e   local %.2e\n", label, a, glob, loc)
end
println("uniform grid, node spacing 1, nr = 400 (a in nodes; a = 10 and 25 resolve the source with 10 and 25 nodes)")
for order in (4, 6), a in (3.0, 10.0, 25.0)
    A = FaradayDisk.generate_dtn_single_layer(400, 800; order = order)
    report("single-layer, order $order", A, collect(0.0:399.0), a)
end
println("paper's polar matrix, nr = 200 (same test)")
for a in (10.0, 25.0)
    report("polar", FaradayDisk.generate_dtn_polar(200, 400), collect(0.0:199.0), a)
end
println("graded grid (radii in disk radii; 25 nodes per disk radius near the disk; ppw 32 at 40 Hz, wall zone 2; a in disk radii)")
p = SimulationParams(bathDensity = 1.148, bathSurfaceTension = 68.0, bathViscosity = 0.079, g = 981.0, diskRadius = 0.07938,
                     bathDiameter = 700.0, spatialResolution = 25.0)
g = graded_grid_for(p, 40.0)
A = generate_dtn(g); r = collect(g.r[1:g.nr])
println("nr = ", g.nr, ", R = ", g.r[g.nr + 1], " disk radii")
for a in (0.4, 1.0, 3.0, 10.0)
    report("graded", A, r, a)
end
