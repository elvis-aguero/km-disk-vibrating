# include() at file top-level (not inside the @testset body) — see
# test_sweep_end_to_end.jl's comment for why: include() inside a running function/closure
# creates a Julia "world-age" barrier that a first CI run against this port actually hit.
# Guarded since test_sweep_end_to_end.jl may have already included the same script in
# this process (redefining `struct OverlayCase` a second time is unsafe).
if !@isdefined(OverlayCase)
    include(joinpath(@__DIR__, "..", "..", "scripts", "run_validation_overlay.jl"))
end

@testset "validation overlay plotting" begin
    # Synthetic OverlayCase data, not a real simulation — isolates the CairoMakie
    # plotting code (the part written most recently, with the least local verification)
    # from simulation correctness, and runs in milliseconds rather than minutes.
    cases = OverlayCase[]
    for g in (0.05, 0.2, 0.5), f in (10.0, 30.0, 60.0, 90.0, 100.0)
        push!(cases, OverlayCase(g, f, 0.8 + 0.1 * sin(f / 20), 10.0 + f / 5))
    end

    exp_amp_path = joinpath(EXPERIMENTAL_DIR, "ampiezza_solo_misure_3.csv")
    exp_phase_path = joinpath(EXPERIMENTAL_DIR, "fase_solo_misure_3.csv")
    @test isfile(exp_amp_path)
    @test isfile(exp_phase_path)

    mktempdir() do out_dir
        plots = plot_validation_overlay(cases, exp_amp_path, exp_phase_path, out_dir, "test_sweep")

        for path in (plots.amp_png, plots.amp_pdf, plots.phase_png, plots.phase_pdf)
            @test isfile(path)
            @test filesize(path) > 0
        end
    end

    @testset "errors on empty input rather than producing a blank plot" begin
        mktempdir() do out_dir
            @test_throws ErrorException plot_validation_overlay(OverlayCase[], exp_amp_path, exp_phase_path, out_dir, "empty")
        end
    end
end
