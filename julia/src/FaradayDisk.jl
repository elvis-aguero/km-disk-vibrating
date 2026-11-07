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
include("grid.jl")                  # RadialGrid needed by domain.jl, parameters.jl and the DtN builders
include("domain.jl")
include(joinpath("dtn", "dtn_registry.jl"))
include(joinpath("dtn", "dtn_generator.jl"))
include(joinpath("dtn", "dtn_single_layer.jl"))
include("convergence.jl")           # ConvergenceOptions needed by SimulationParams below
include("parameters.jl")            # Problem needed by system_matrix.jl below
include("system_matrix.jl")         # SystemMatrixTemplate needed by lu_cache.jl / solver_dispatch.jl
include(joinpath("solver", "lu_cache.jl"))
include(joinpath("solver", "gmres_solver.jl"))
include(joinpath("solver", "static_equilibrium.jl"))
include(joinpath("solver", "solver_dispatch.jl"))
include(joinpath("solver", "schur_solver.jl"))   # reduced (nr+1) system; needs AbstractStepSolver, LUCache
include("state.jl")
include("simulate.jl")
include("sweep.jl")                 # SweepCaseResult needed by io.jl below
include("io.jl")
include("logging_setup.jl")

export SimulationParams, Problem, Units, build_problem, with_overrides,
       ConvergenceOptions, ConvergenceResult, fit_oscillation, check_convergence,
       SimulationState, StepBuffers, init_state, advance_one_step!,
       SimulationResult, run_simulation, lab_frame_CoM,
       SweepSpec, SweepCaseResult, run_sweep, run_sweep_graded, run_sweep_case, case_id, case_csv_path,
       DTNManifest, DTNCacheEntry, default_registry, resolve_dtn, load_dtn, check_dtn, register!,
       load_manifest, save_manifest, default_manifest_path,
       generate_dtn, generate_dtn_polar, generate_dtn_single_layer, dtn_row, default_dtn_method,
       build_domain, contact_weights,
       RadialGrid, uniform_grid, graded_grid, graded_grid_for, capillary_wavelength_cm, boundary_node_count, boundary_stencil, wall_radius,
       assemble_system_matrix, assemble_system_matrix!, build_template, materialize!, repatch!, assemble_rhs!,
       LUCache, get_or_factorize!, hit_rate, n_factorizations, lu_factor_count, G_DEDUP_ATOL, mean_iterations,
       StepPhase, step_phase, g_prefactor, force_term, resolve_solver_type, estimated_lu_cache_bytes, estimated_schur_bytes, schur_applicable, build_schur_system, schur_solve!, solve_startup!, SchurStepSolver, build_step_solver, implicit_dt,
       available_memory_bytes, MEMORY_SAFETY_FACTOR,
       compute_eta_boundary, write_case_csv, write_timeseries_csv, write_summary_csv,
       config_digest, run_id, save_result,
       setup_logging, default_log_dir

end # module FaradayDisk
