# run_sweep_graded: one grid and one DtN per frequency, all cases in one outdir, one summary. Small domain so it is fast.
@testset "run_sweep_graded" begin
    mktempdir() do outdir
        fixed = SimulationParams(diskRadius = 0.2, diskMass = 0.0283, phaseDifference = -90.0, bathDensity = 1.175,
                                 bathSurfaceTension = 66.5, bathViscosity = 0.18 / 1.175, bathDiameter = 20.0,
                                 spatialResolution = 5.0, temporalResolution = 40, solverType = :auto, startStatic = true,
                                 earlyStop = false)
        spec = SweepSpec(gammas = [0.1, 0.3], bathFrequenciesHz = [30.0, 60.0], fixed = fixed)
        results = run_sweep_graded(spec, outdir)
        @test length(results) == 4
        @test all(r -> r.status == :ok, results)
        for g in spec.gammas, f in spec.bathFrequenciesHz
            @test isfile(case_csv_path(outdir, g, f))
        end
        # the summary lists every frequency, not only the last one
        rows = readlines(joinpath(outdir, "summary.csv"))
        @test length(rows) == 1 + 4
        # a spec that already carries a grid is refused
        gridded = with_overrides(fixed; grid = graded_grid_for(fixed, 30.0))
        @test_throws ArgumentError run_sweep_graded(SweepSpec(gammas = [0.1], bathFrequenciesHz = [30.0], fixed = gridded), outdir)
    end
end
