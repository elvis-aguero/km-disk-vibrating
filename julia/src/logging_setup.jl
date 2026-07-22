"""
Structured logging setup, replacing MATLAB's `fprintf` + `diary` (a byte-for-byte console
capture with no structure — answering "which cases converged before period 5?" from an
old MATLAB log means regexing fragile printf strings).

Builds a `TeeLogger`: a colorized console sink and a per-run/per-sweep timestamped file
sink, both receiving the same structured `@info "message" key=val ...` calls used
throughout `simulate.jl`/`sweep.jl`. The sweep orchestrator additionally gets a
JSON-lines sink (`jsonlines=true`) — sweep progress (30+ cases, meant to be watched or
fed into a future dashboard) benefits from being machine-parseable; per-step physics
logging for a single run does not, and stays human-readable-only.

CAVEAT: written without the ability to install/run LoggingExtras.jl to verify the exact
`TeeLogger`/`FormatLogger` constructor signatures (see the module-level caveat in
`solver/gmres_solver.jl` for why). `Logging.ConsoleLogger` itself is Julia stdlib and
very stable; if anything here fails to load, it is most likely the `LoggingExtras` calls,
not the stdlib ones.
"""

"""
    default_log_dir() -> String

`\$KMDISK_LOG_DIR` if set, else `julia/logs/` (created if missing).
"""
function default_log_dir()
    dir = get(ENV, "KMDISK_LOG_DIR", joinpath(@__DIR__, "..", "logs"))
    mkpath(dir)
    return dir
end

_json_value(v::AbstractString) = string('"', escape_string(v), '"')
_json_value(v::Symbol) = string('"', String(v), '"')
_json_value(v::Bool) = v ? "true" : "false"
_json_value(v::Real) = isfinite(v) ? string(v) : "null"
_json_value(v) = string('"', escape_string(string(v)), '"')

"""
    setup_logging(; log_dir=default_log_dir(), name="run", console_level=Logging.Info,
                    file_level=Logging.Debug, jsonlines=false) -> AbstractLogger

Builds (but does not install) a logger. Install it globally with
`global_logger(setup_logging(...))`, or scope it locally (e.g. in tests, so global logger
state isn't polluted) with `with_logger(setup_logging(...)) do ... end`.
"""
function setup_logging(; log_dir::AbstractString = default_log_dir(), name::AbstractString = "run",
                          console_level::Logging.LogLevel = Logging.Info,
                          file_level::Logging.LogLevel = Logging.Debug, jsonlines::Bool = false)
    mkpath(log_dir)
    ts = Dates.format(Dates.now(), "yyyymmdd_HHMMSS")

    console_logger = Logging.ConsoleLogger(stderr, console_level)
    file_io = open(joinpath(log_dir, "$(name)_$(ts).log"), "w")
    file_logger = Logging.ConsoleLogger(file_io, file_level)

    combined = LoggingExtras.TeeLogger(console_logger, file_logger)

    if jsonlines
        jsonl_io = open(joinpath(log_dir, "$(name)_$(ts).jsonl"), "w")
        jsonl_logger = LoggingExtras.FormatLogger(jsonl_io) do io, log
            kvs = join((string('"', k, "\":", _json_value(v)) for (k, v) in log.kwargs), ",")
            msg = _json_value(string(log.message))
            body = isempty(kvs) ? "" : string(",", kvs)
            println(io, "{\"level\":\"", log.level, "\",\"message\":", msg, body, "}")
            flush(io)
        end
        combined = LoggingExtras.TeeLogger(combined, jsonl_logger)
    end

    return combined
end
