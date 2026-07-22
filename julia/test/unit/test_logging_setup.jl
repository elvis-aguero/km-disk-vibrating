using Logging  # with_logger() below; not re-exported by `using FaradayDisk`

@testset "logging_setup: sinks are flushed immediately, not buffered until exit" begin
    # This exists because a real cluster run showed nothing in its human-readable log or
    # SLURM .err output for hours, while the .jsonl file (which explicitly flushed itself)
    # updated fine -- Logging.ConsoleLogger does not flush its stream after each message,
    # which matters a great deal once that stream is a *file* rather than a live terminal
    # (exactly what stderr becomes under SLURM, and what the .log file always is): messages
    # can sit in an OS buffer for an entire run instead of appearing as they're logged.
    mktempdir() do log_dir
        logger = setup_logging(log_dir = log_dir, name = "flushtest", jsonlines = true)

        with_logger(logger) do
            @info "first message" step = 1
        end

        log_files = filter(f -> endswith(f, ".log"), readdir(log_dir; join = true))
        jsonl_files = filter(f -> endswith(f, ".jsonl"), readdir(log_dir; join = true))
        @test length(log_files) == 1
        @test length(jsonl_files) == 1

        # The critical assertion: read the file *without closing anything* on the writer
        # side (with_logger has already returned, but nothing has been explicitly closed
        # or flushed by this test) -- if FlushingLogger's flush(l.stream) call were removed,
        # this content would be empty until the file handle is later closed/GC'd.
        @test !isempty(read(log_files[1], String))
        @test !isempty(read(jsonl_files[1], String))
        @test occursin("first message", read(log_files[1], String))
    end
end
