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

Every sink flushes its stream after each message (see `FlushingLogger` below) — this
matters a great deal in practice: `Logging.ConsoleLogger` does NOT do this itself, and
when its stream is a file rather than a live terminal (exactly what `stderr` is once
SLURM redirects a job's output, and what the `.log` file always is), messages can sit in
an OS-level buffer for the entire run instead of appearing as they're logged. Confirmed
directly: a bare `ConsoleLogger(stderr, ...)` job left a redirected `.err` file completely
empty several seconds after logging a message and before exiting. Only the JSON-lines sink
avoided this, because it already called `flush(io)` itself — which is exactly why, before
this fix, the `.jsonl` file was the only thing that appeared to update during a long run,
while the human-readable `.log`/console output looked silent for hours.
"""

"""
    FlushingLogger(inner, stream)

Wraps any `AbstractLogger` so its destination `stream` is flushed after every message —
see this file's module docstring for why that's not just a nice-to-have.
"""
struct FlushingLogger{L<:Logging.AbstractLogger,S<:IO} <: Logging.AbstractLogger
    inner::L
    stream::S
end

Logging.shouldlog(l::FlushingLogger, args...) = Logging.shouldlog(l.inner, args...)
Logging.min_enabled_level(l::FlushingLogger) = Logging.min_enabled_level(l.inner)
Logging.catch_exceptions(l::FlushingLogger) = Logging.catch_exceptions(l.inner)

function Logging.handle_message(l::FlushingLogger, args...; kwargs...)
    Logging.handle_message(l.inner, args...; kwargs...)
    flush(l.stream)
    return nothing
end

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

    console_logger = FlushingLogger(Logging.ConsoleLogger(stderr, console_level), stderr)
    file_io = open(joinpath(log_dir, "$(name)_$(ts).log"), "w")
    file_logger = FlushingLogger(Logging.ConsoleLogger(file_io, file_level), file_io)

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
