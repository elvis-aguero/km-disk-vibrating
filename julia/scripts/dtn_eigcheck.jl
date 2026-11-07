#= Spectrum of the DtN matrix. The exact operator is self-adjoint and positive (symbol |k|), and the wave problem
eta_t = N phi, phi_t = -(g + sigma(-lap)) eta is stable only if the discrete N keeps real positive eigenvalues
(similar to a symmetric positive matrix). Any eigenvalue with Re < 0 or a sizeable imaginary part is a growing mode.
Units: node spacing 1 (generate_dtn_single_layer(nr, 2nr)). Usage: julia -t N julia/scripts/dtn_eigcheck.jl nr1 nr2 ... =#
include(joinpath(@__DIR__, "_bootstrap.jl"))
using LinearAlgebra, Printf

function report(label, A)
    t = @elapsed ev = eigvals(A)
    re, im_ = real.(ev), imag.(ev)
    cplx = count(i -> abs(im_[i]) > 1e-8 * abs(ev[i]), eachindex(ev))
    @printf("%-28s n=%5d  min Re=%.4e  max Re=%.4e  #Re<0=%d  #complex=%d  max|Im|=%.3e  (%.0f s)\n",
            label, length(ev), minimum(re), maximum(re), count(<(0), re), cplx, maximum(abs, im_), t)
    flush(stdout)
    neg = sortperm(re)[1:min(5, length(re))]
    for i in neg; @printf("    lowest: %.4e %+.4e i\n", re[i], im_[i]); end
    flush(stdout)
end

for nr in parse.(Int, ARGS)
    for order in (6, 4)
        A = FaradayDisk.generate_dtn_single_layer(nr, 2nr; order = order)
        report("single-layer order $order nr=$nr", A)
    end
    if nr <= 2500
        try
            reg = FaradayDisk.default_registry()
            e = FaradayDisk.resolve_dtn(reg, 25.0, 2nr / 25; method = :polar)
            report("polar cache nr=$nr", FaradayDisk.load_dtn(e))
        catch err
            println("polar cache nr=$nr unavailable: ", sprint(showerror, err))
        end
    end
end
