"""
User-facing simulation parameters and the derived, immutable `Problem` (nondimensionalized
constants + discretization operators) built from them.

Mirrors the `arguments` block at the top of `matlab/1_code/simulation_code/solve_motion.m`
and the nondimensionalization it performs (Reynolds/Froude/Weber numbers, characteristic
units), but keeps mutable per-step state (the LU cache, the GMRES warm-start vector, the
step counter) entirely out of this struct — MATLAB's `PROBLEM_CONSTANTS` mixes the two,
which is why `solve_motion.m` has to strip several fields (`L_Library`, `U_Library`,
`P_Library`, `luFutures`) before every save. `Problem` here holds only what never changes
once a run starts; genuinely mutable state lives in `SimulationState` (`state.jl`) and the
solver caches (`solver/`).
"""

Base.@kwdef struct SimulationParams
    diskRadius::Float64 = 0.2                 # cm
    diskMass::Float64 = 0.0283                # g
    forceAmplitude::Float64 = 0.0             # dyne
    forceFrequency::Float64 = 90.0             # Hz
    bathAmplitude::Float64 = 0.01             # cm
    bathFrequency::Float64 = 90.0              # Hz
    phaseDifference::Float64 = -90.0          # degrees
    bathDensity::Float64 = 1.0                # g/cm^3
    bathSurfaceTension::Float64 = 72.20       # dyne/cm
    bathViscosity::Float64 = 0.978e-2         # Stokes
    g::Float64 = 981.0                        # cm/s^2
    bathDiameter::Float64 = 100.0              # disk radii
    spatialResolution::Float64 = 50.0          # radial mesh points per disk radius
    temporalResolution::Int = 30               # steps per adimensional time unit
    simulationTime::Float64 = 10 / 90          # s
    solverType::Symbol = :auto                 # :auto | :lu | :gmres
    gmresTolerance::Float64 = 1e-10
    startStatic::Bool = true
    earlyStop::Bool = true
    convergence::ConvergenceOptions = ConvergenceOptions()
end

struct Units
    length::Float64
    mass::Float64
    time::Float64
    velocity::Float64
    force::Float64
end

struct Problem
    Re::Float64
    Fr::Float64
    We::Float64
    nr::Int
    cPoints::Int
    dr::Float64
    laplacian::Tridiagonal{Float64,Vector{Float64}}
    DTN::Matrix{Float64}
    pressureIntegral::Vector{Float64}          # pressureIntegral(cPoints, 1:cPoints) row, length cPoints
    force_amplitude::Float64                    # force_adim
    force_frequency::Float64                     # freq_adim (always 1.0 by construction: see build_problem)
    bath_forcing_amplitude::Float64             # gamma = A*omega_bath^2/g, dimensionless
    bath_frequency::Float64                      # bath_freq_adim
    phase_difference::Float64                     # radians
    surface_force_constant::Float64              # surface_force_adim
    obj_mass::Float64                             # obj_mass_adim
    stepsPerCycle::Int
    dt::Float64
    effective_w::Float64                          # dimensionless angular frequency governing periodicity
    v_bath_0::Float64                             # initial CoM velocity offset (v_bath_0_adim)
    units::Units
end

"""
    build_problem(params, dtn) -> Problem

Computes all nondimensionalized constants and builds the spatial discretization, given a
DTN operator matrix already resolved via the DTN registry (see `dtn_registry.jl`) for
`(params.spatialResolution, params.bathDiameter)`.

Note `force_frequency` (`freq_adim` in the MATLAB source) is always exactly `1.0`: the time
unit `T_unit` is defined as `1/forceFrequency`, so the disk's own forcing frequency is the
reference clock by construction. This looks redundant but is intentional — preserved
exactly as in `solve_motion.m` rather than "simplified away," since `effective_w` (which
governs step size and cycle periodicity) may instead be driven by `bath_frequency` when
bath oscillation is present.
"""
function build_problem(params::SimulationParams, dtn::AbstractMatrix{<:Real})
    L_unit = params.diskRadius
    M_unit = params.bathDensity * L_unit^3
    forceFrequency_angular = params.forceFrequency * 2 * pi
    T_unit = 1 / forceFrequency_angular
    V_unit = L_unit / T_unit
    F_unit = M_unit * L_unit / T_unit^2
    units = Units(L_unit, M_unit, T_unit, V_unit, F_unit)

    Re = L_unit^2 / (params.bathViscosity * T_unit)
    Fr = L_unit / (params.g * T_unit^2)
    We = params.bathDensity * L_unit^3 / (params.bathSurfaceTension * T_unit^2)

    force_adim = params.forceAmplitude / params.diskMass * (T_unit^2 / L_unit)
    surface_force_adim = 2 * pi * params.diskRadius * params.bathSurfaceTension / params.diskMass * (T_unit^2 / L_unit)
    freq_adim = forceFrequency_angular * T_unit   # == 1.0 always; see docstring
    obj_mass_adim = params.diskMass / M_unit

    bath_angular_freq = params.bathFrequency * 2 * pi
    bath_freq_adim = bath_angular_freq * T_unit
    phase_diff_rad = params.phaseDifference * pi / 180
    bath_forcing_amplitude = params.bathAmplitude * bath_angular_freq^2 / params.g

    v_bath_0_adim = (params.bathAmplitude * bath_angular_freq * sin(phase_diff_rad)) / V_unit

    effective_w_adim = if bath_forcing_amplitude == 0
        freq_adim
    else
        bath_freq_adim  # covers both "forceAmplitude==0" and the mixed-forcing case, matching solve_motion.m
    end

    stepsPerCycle = params.temporalResolution
    dt = (2 * pi / effective_w_adim) / stepsPerCycle

    dom = build_domain(params.bathDiameter, params.spatialResolution)
    cPoints = Int(params.spatialResolution) + 1
    size(dtn, 1) == dom.nr || throw(DimensionMismatch(
        "DTN matrix has size $(size(dtn)) but domain requires nr=$(dom.nr) " *
        "(spatialResolution=$(params.spatialResolution), bathDiameter=$(params.bathDiameter))"))
    pressure_integral_row = dom.IntMat[cPoints, 1:cPoints]

    return Problem(Re, Fr, We, dom.nr, cPoints, dom.dr, dom.Laplacian, Matrix{Float64}(dtn),
                   pressure_integral_row, force_adim, freq_adim, bath_forcing_amplitude,
                   bath_freq_adim, phase_diff_rad, surface_force_adim, obj_mass_adim,
                   stepsPerCycle, dt, effective_w_adim, v_bath_0_adim, units)
end

"""
    with_overrides(params::SimulationParams; kwargs...) -> SimulationParams

Returns a copy of `params` with the given fields replaced. `Base.@kwdef` only generates a
keyword constructor from defaults, not a "copy with changes" constructor, so this fills
that gap generically (using the always-available positional constructor).
"""
function with_overrides(params::SimulationParams; kwargs...)
    fields = fieldnames(SimulationParams)
    values = Dict{Symbol,Any}(f => getfield(params, f) for f in fields)
    for (k, v) in kwargs
        haskey(values, k) || throw(ArgumentError("SimulationParams has no field $k"))
        values[k] = v
    end
    return SimulationParams((values[f] for f in fields)...)
end
