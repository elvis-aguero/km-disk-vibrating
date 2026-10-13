using FaradayDisk, LinearAlgebra, Krylov, Printf, Logging
global_logger(Logging.NullLogger())
println("BLAS threads ", BLAS.get_num_threads())
# A. raw kernel costs at the production size n = 5002 (cost depends on n, not values)
n = 5002; A = randn(n, n) + n * I; b = randn(n); x = similar(b)
lu(A[1:100, 1:100]); t_lu = @elapsed F = lu(A)
ldiv!(x, F, b); t_ld = minimum(@elapsed(ldiv!(x, F, b)) for _ in 1:5)
mul!(x, A, b); t_mv = minimum(@elapsed(mul!(x, A, b)) for _ in 1:5)
@printf("n=5002: lu %.2f s | ldiv %.1f ms | matvec %.1f ms\n", t_lu, 1e3t_ld, 1e3t_mv)
A = nothing; F = nothing; GC.gc()
# B. preconditioned-GMRES iteration count vs nr (worst gamma)
for nr in (100, 200, 300)
    D = 0.4nr; dtn = generate_dtn(nr, D)
    for f in (10.0, 100.0)
        p = SimulationParams(bathDiameter = D, spatialResolution = 2nr / D, temporalResolution = 60,
                             bathDensity = 1.175, bathSurfaceTension = 66.5, bathViscosity = 0.18 / 1.175,
                             bathAmplitude = 0.5 * 981 / (2pi * f)^2, bathFrequency = f)
        prob = build_problem(p, dtn); tm = FaradayDisk.build_template(prob); M = similar(tm.Mat0)
        FaradayDisk.materialize!(M, tm, 1.0); F1 = lu(M)
        its = Int[]
        for k in 0:59
            g = g_prefactor(prob, step_phase(k, 60), prob.dt)
            FaradayDisk.materialize!(M, tm, g)
            _, st = Krylov.gmres(M, randn(size(M, 1)); N = F1, ldiv = true, atol = 1e-12, rtol = 1e-12, memory = 50, itmax = 200)
            push!(its, st.niter)
        end
        @printf("nr=%d f=%5.1f gamma=0.5: pGMRES iters min %d max %d\n", nr, f, minimum(its), maximum(its))
    end
end
