module FaradayDisk

using LinearAlgebra
using Krylov
using JLD2
using TOML
using SHA
using Dates
using Printf
using Logging
using LoggingExtras

# Include order matters: struct field types and function-signature type annotations must
# already exist at include time (unlike plain function-body references, which only need
# to resolve by call time). See the comment in each file for what it depends on.
include("domain.jl")
include(joinpath("dtn", "dtn_registry.jl"))
include(joinpath("dtn", "dtn_generator.jl"))
include("convergence.jl")           # ConvergenceOptions needed by SimulationParams below
include("parameters.jl")            # Problem needed by system_matrix.jl below
include("system_matrix.jl")         # SystemMatrixTemplate needed by lu_cache.jl / solver_dispatch.jl
include(joinpath("solver", "lu_cache.jl"))
include(joinpath("solver", "gmres_solver.jl"))
include(joinpath("solver", "static_equilibrium.jl"))
include(joinpath("solver", "solver_dispatch.jl"))
include("state.jl")
include("simulate.jl")
include("sweep.jl")                 # SweepCaseResult needed by io.jl below
include("io.jl")
include("logging_setup.jl")

export SimulationParams, Problem, Units, build_problem, with_overrides,
       ConvergenceOptions, ConvergenceResult, fit_oscillation, check_convergence,
       SimulationState, StepBuffers, init_state, advance_one_step!,
       SimulationResult, run_simulation,
       SweepSpec, SweepCaseResult, run_sweep, run_sweep_case, case_id, case_csv_path,
       DTNManifest, DTNCacheEntry, default_registry, resolve_dtn, load_dtn, register!,
       load_manifest, save_manifest, default_manifest_path,
       generate_dtn,
       build_domain,
       assemble_system_matrix, assemble_system_matrix!, build_template, materialize!, assemble_rhs!,
       LUCache, get_or_factorize!, hit_rate,
       StepPhase, g_prefactor, force_term, resolve_solver_type, estimated_lu_cache_bytes, build_step_solver,
       compute_eta_boundary, write_case_csv, write_timeseries_csv, write_summary_csv,
       config_digest, run_id, save_result,
       setup_logging, default_log_dir

end # module FaradayDisk
