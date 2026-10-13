# Experiment: can ONE LU factorization (at g_prefactor=1) precondition GMRES for every
# other g_prefactor in the cycle? And how many LU-cache entries does by_value dedupe?
using FaradayDisk, LinearAlgebra, Krylov, Printf, Logging
global_logger(Logging.NullLogger())
nr, D = 150, 60.0
dtn = generate_dtn(nr, D)
for gamma in (0.05, 0.2, 0.5), f in (10.0, 50.0, 100.0)
    p = SimulationParams(bathDiameter = D, spatialResolution = 2nr / D, temporalResolution = 60,
                         bathDensity = 1.175, bathSurfaceTension = 66.5, bathViscosity = 0.18 / 1.175,
                         bathAmplitude = gamma * 981 / (2pi * f)^2, bathFrequency = f)
    prob = build_problem(p, dtn); n = 2prob.nr + 2
    tmpl = FaradayDisk.build_template(prob)
    M = similar(tmpl.Mat0)
    gs = [g_prefactor(prob, step_phase(k, 60), prob.dt) for k in 0:59]
    nuniq = length(unique(gs)); nuniq_r = length(unique(round.(gs, sigdigits = 12)))
    FaradayDisk.materialize!(M, tmpl, 1.0); F1 = lu(M)
    maxit = 0; maxerr = 0.0
    for g in gs
        FaradayDisk.materialize!(M, tmpl, g)
        b = randn(n)
        # right-preconditioned GMRES with N = F1 (ldiv)
        x, st = Krylov.gmres(M, b; N = F1, ldiv = true, atol = 1e-12, rtol = 1e-12, memory = 50, itmax = 200)
        maxit = max(maxit, st.niter)
        maxerr = max(maxerr, norm(x - M \ b) / norm(M \ b))
    end
    @printf("gamma=%.2f f=%5.1f | unique g: %d exact, %d @12sig | precond GMRES max iters %d, max rel err vs direct %.1e | cond(M)=%.1e\n",
            gamma, f, nuniq, nuniq_r, maxit, maxerr, cond(M, 1))
end
