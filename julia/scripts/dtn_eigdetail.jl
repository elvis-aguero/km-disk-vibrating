#= Where are the complex eigenvalues of the single-layer DtN matrix, and where do their eigenvectors live?
Prints, for nr = ARGS[1] (node spacing 1): the 12 eigenvalues with the largest |Im|/Re, the node of the largest
|component|, and the share of the squared norm in the last 5% of the nodes (a wall-localized mode has a share near 1). =#
include(joinpath(@__DIR__, "_bootstrap.jl"))
using LinearAlgebra, Printf
nr = parse(Int, ARGS[1]); order = length(ARGS) > 1 ? parse(Int, ARGS[2]) : 6
A = FaradayDisk.generate_dtn_single_layer(nr, 2nr; order = order)
E = eigen(A)
ratio = abs.(imag.(E.values)) ./ abs.(real.(E.values))
@printf("nr=%d order=%d: #complex=%d, max |Im|/Re=%.3e\n", nr, order, count(>(1e-8), ratio), maximum(ratio))
last5 = (nr - nr ÷ 20 + 1):nr
for i in sortperm(ratio; rev = true)[1:12]
    v = E.vectors[:, i]; w = abs2.(v) ./ sum(abs2, v)
    @printf("  Re=%.4e Im=%+.4e |Im|/Re=%.3e  node of max=%5d  share in last 5%%=%.3f\n",
            real(E.values[i]), imag(E.values[i]), ratio[i], argmax(w) - 1, sum(w[last5]))
end
