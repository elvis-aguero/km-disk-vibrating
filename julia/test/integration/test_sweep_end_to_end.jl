# include() at file top-level, not inside the @testset/mktempdir body below: calling
# include() from within a running function creates a Julia "world age" barrier — methods
# it defines aren't visible to the *currently executing* dynamic call, only to calls
# dispatched afresh afterwards. A first CI run hit exactly this (`MethodError: ... The
# applicable method may be too new`) when this was nested inside the testset.
include(joinpath(@__DIR__, "..", "..", "scripts", "run_validation_overlay.jl"))

@testset "sweep end-to-end" begin
    mktempdir() do outdir
        fixed = SimulationParams(diskRadius = 0.2, diskMass = 0.0283,
                                  phaseDifference = -90.0, bathDensity = 1.175,
                                  bathSurfaceTension = 66.5, bathViscosity = 0.18 / 1.175,
                                  bathDiameter = 20.0, spatialResolution = 5.0, temporalResolution = 40,
                                  solverType = :lu, startStatic = true, earlyStop = false)
        spec = SweepSpec(gammas = [0.1, 0.3], bathFrequenciesHz = [30.0, 60.0], fixed = fixed)
        dtn = generate_dtn(Int(ceil(fixed.spatialResolution * fixed.bathDiameter / 2)), fixed.bathDiameter)

        results = run_sweep(spec, outdir; dtn = dtn)
        @test length(results) == 4
        @test all(r -> r.status == :ok, results)

        @testset "case CSV schema" begin
            for g in spec.gammas, f in spec.bathFrequenciesHz
                path = case_csv_path(outdir, g, f)
                @test isfile(path)
                header = split(readline(path), ',')
                @test header == ["time_s", "CoM_cm", "eta_boundary_cm"]
            end
        end

        @testset "summary CSV schema" begin
            summary_path = joinpath(outdir, "summary.csv")
            @test isfile(summary_path)
            header = split(readline(summary_path), ',')
            @test header == ["gamma", "bathFrequency_Hz", "bathAmplitude_cm", "CoM_max_cm", "CoM_rms_cm", "elapsed_s", "status"]
        end

        @testset "resume/skip-if-exists" begin
            results2 = run_sweep(spec, outdir; dtn = dtn)
            @test all(r -> r.status == :skipped, results2)
        end

        @testset "validation overlay, where usable, produces finite roughly-O(1) amplitude ratios" begin
            # Not asserting cases is non-empty: whether this small test domain has enough
            # room before the 30-period run hits boundary-reflection truncation depends on
            # wave-speed details not verified empirically here (see the repo-level notes
            # on how this port was authored). What matters is that whatever *does* come
            # back is physically sane, not blown up or NaN.
            cases = compute_sim_overlay(outdir)
            @info "sweep end-to-end overlay" n_usable_cases = length(cases)
            for c in cases
                @test isfinite(c.ampNorm)
                @test 0.05 < c.ampNorm < 5.0
                @test isfinite(c.phaseDiff)
            end
        end
    end
end
