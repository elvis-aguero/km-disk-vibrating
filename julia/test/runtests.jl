using Test
using FaradayDisk

# Tiering: `] test` (or `julia test/runtests.jl`) runs the fast tier (unit + integration)
# by default. Set KMDISK_RUN_SLOW_TESTS=1 to also run the slow tier (nr=2500 DTN
# comparison, full 30-case MATLAB-baseline validation sweep) — these are genuinely
# multi-hour-scale and must never run on every push. There is no golden/frozen numeric
# fixture tier beyond the DTN generator comparison: correctness for the solver,
# convergence check, and sweep pipeline is established via physically-grounded properties
# (does the system oscillate at the forcing frequency? does damping actually damp? does
# GMRES agree with LU-caching?) rather than frozen "the answer was N last time" regression
# numbers.
const RUN_SLOW = get(ENV, "KMDISK_RUN_SLOW_TESTS", "0") == "1"

@testset "FaradayDisk" begin
    @testset "unit" begin
        include("unit/test_parameters.jl")
        include("unit/test_domain.jl")
        include("unit/test_system_matrix.jl")
        include("unit/test_fit.jl")
        include("unit/test_convergence_check.jl")
        include("unit/test_solver_dispatch.jl")
    end

    @testset "integration" begin
        include("integration/test_dtn_golden_small.jl")
        include("integration/test_lu_cache.jl")
        include("integration/test_physics_grounded.jl")
        include("integration/test_sweep_end_to_end.jl")
        include("integration/test_plotting.jl")
    end

    if RUN_SLOW
        @testset "slow" begin
            include("slow/test_dtn_golden_large.jl")
            include("slow/test_matlab_baseline_validation.jl")
        end
    else
        @info "Skipping slow test tier (nr=2500 DTN comparison, full MATLAB-baseline sweep) — set KMDISK_RUN_SLOW_TESTS=1 to run it."
    end
end
