#!/usr/bin/env julia
#= Computes the Julia-side reference values for JuliaParityTest.m's cross-implementation
parity check, run fresh on every invocation (never a cached/precomputed value, so it can
never go stale relative to whatever this branch's Julia code actually does right now).

Runs the same bath-driven case JuliaParityTest.m runs through solve_motion.m, at the
cheapest domain provably shared between the two implementations (D5Quant20, nr=50 --
cross-validated byte-for-byte between Julia's native DTN generator and MATLAB's original
.mat cache earlier this session), and prints exactly one line of space-separated numbers
to stdout: t_end CoM_first amplitude phase offset. Nothing else goes to stdout (the
default logger is silenced) so a caller like MATLAB's `system(...)` gets a clean,
trivially parseable line.

Block comment, not a docstring: see scripts/_bootstrap.jl's header for why a bare
top-level string directly followed by a non-declaration statement here would fail to
parse. =#

import Pkg
Pkg.activate(joinpath(@__DIR__, ".."))
using FaradayDisk
using Logging
global_logger(Logging.NullLogger())  # keep stdout clean for the caller to parse

freqHz = 30.0
gamma = 0.2
omega = 2 * pi * freqHz
A_bath = gamma * 981.0 / omega^2
nPeriods = 8

p = SimulationParams(diskRadius = 0.2, diskMass = 0.0283, forceAmplitude = 0.0, forceFrequency = 90.0,
                      bathAmplitude = A_bath, bathFrequency = freqHz, phaseDifference = -90.0,
                      bathDensity = 1.175, bathSurfaceTension = 66.5, bathViscosity = 0.18 / 1.175,
                      bathDiameter = 20.0, spatialResolution = 5.0, temporalResolution = 48,
                      simulationTime = nPeriods / freqHz, startStatic = true, earlyStop = false, solverType = :lu)

# Generated natively (not via resolve_dtn/default_registry()) because julia/data/dtn_cache/
# is gitignored (regenerate-locally-only, see .gitignore) and so isn't present on a fresh CI
# checkout. nr=50 for this domain matches build_domain's nr = ceil(D*quant/2).
#
# Generated at D=5, NOT p.bathDiameter (20) -- migrate_dtn_caches.jl documents that this
# domain's legacy MATLAB cache (D5Quant20/DTNnew345nr50D5refp10.mat) was actually generated
# with the literal argument D=5 (its filename is a mechanically exact record of that), and
# solve_motion.m has a "Machine-specific patch for D5Quant20" that always loads this same
# D=5-generated matrix whenever spatialResolution==5 && bathDiameter==20, never regenerating
# it at the nominal bathDiameter. Using D=20 here disagrees from what MATLAB actually loads
# by ~75% (confirmed empirically) -- D=5 is what makes this a genuine parity check.
dtn = generate_dtn(50, p.spatialResolution)
result = run_simulation(p; dtn = dtn)

T = 1 / freqHz
tail = result.t_s .>= (result.t_s[end] - 3 * T)
fit = fit_oscillation(result.t_s[tail], result.CoM_cm[tail], omega)

println(result.t_s[end], " ", result.CoM_cm[1], " ", fit.amplitude, " ", fit.phase, " ", fit.offset)
