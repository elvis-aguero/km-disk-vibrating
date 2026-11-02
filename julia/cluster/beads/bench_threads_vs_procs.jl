#!/usr/bin/env julia
#= Falsification test for the explanation of job 8262774 (cache generation at 1.00 busy core
with 32 threads). Computes the heaviest rows of the real D25Q200 problem (nr=2500, D=200):
  (1) serially, on a light, a middle and the heaviest row  -> real cost per row;
  (2) with N threads (Threads.@threads :dynamic, the pattern of `generate_dtn`);
  (3) with N worker processes (pmap, the pattern of regenerate_dtn.jl).
Hypothesis: threads speed up far less than N, processes close to N.
Also reports peak memory (master and each worker) and checks that (2) and (3) give identical
rows. Usage: `BENCH_N=8 julia -t 8 bench_threads_vs_procs.jl`
Output is flushed after every line, so a log watcher sees it at once.
Workers are started at top level: `@everywhere using` is not allowed inside functions or
`@elapsed` blocks (the error that killed the first self-test). =#

using Distributed
include(joinpath(@__DIR__, "..", "..", "scripts", "_bootstrap.jl"))

const NR = 2500
const BD = 200.0
const N = parse(Int, get(ENV, "BENCH_N", "8"))
const SCRIPTS_PROJECT = abspath(joinpath(@__DIR__, "..", "..", "scripts"))

say(x) = (println(x); flush(stdout))

Threads.nthreads() == N || @warn "julia was started with $(Threads.nthreads()) threads, BENCH_N=$N"
rows_hi = collect(NR-N+1:NR)

# compile everything on a tiny problem, so timings do not include compilation
FaradayDisk.dtn_row(50, 20.0, 10)

say("== (1) serial cost per row (one thread) ==")
t_serial = Dict{Int,Float64}()
for k in (4, 1250, NR)
    t_serial[k] = @elapsed FaradayDisk.dtn_row(NR, BD, k)
    say("row $k: $(round(t_serial[k]; digits = 2)) s")
end
# row cost is ~linear in k: mean over rows 1..2500 by Simpson's rule on k = 4, 1250, 2500
est_total_core_s = 2500 * (t_serial[4] + 4 * t_serial[1250] + t_serial[NR]) / 6
say("estimated total for 2500 rows: $(round(est_total_core_s / 3600; digits = 2)) core-hours")

t_one = t_serial[NR]
say("\n== (2) $N threads on the $N heaviest rows ==")
res_thr = Vector{Vector{Float64}}(undef, N)
t_thr = @elapsed Threads.@threads :dynamic for i in 1:N
    res_thr[i] = FaradayDisk.dtn_row(NR, BD, rows_hi[i])
end
say("wall $(round(t_thr; digits = 2)) s; ideal (one row's time) $(round(t_one; digits = 2)) s; speedup vs serial rows = $(round(N * t_one / t_thr; digits = 2)) of $N")
say("master peak RSS after the thread run: $(round(Sys.maxrss() / 2^30; digits = 2)) GiB")

say("\n== (3) $N worker processes on the same rows ==")
t_w0 = time()
addprocs(N; exeflags = ["--project=$SCRIPTS_PROJECT", "-t1", "--startup-file=no"])
@everywhere using FaradayDisk
@everywhere FaradayDisk.dtn_row(50, 20.0, 10)
@everywhere begin
    _bench_row(nr::Int, bd::Float64, k::Int) = FaradayDisk.dtn_row(nr, bd, k)
end
say("worker start-up and compilation: $(round(time() - t_w0; digits = 1)) s")
t_p0 = time()
res_proc = pmap(k -> _bench_row(NR, BD, k), rows_hi)
t_proc = time() - t_p0
say("wall $(round(t_proc; digits = 2)) s; speedup vs serial rows = $(round(N * t_one / t_proc; digits = 2)) of $N")
rss = [remotecall_fetch(() -> Sys.maxrss(), w) for w in workers()]
say("worker peak RSS: min $(round(minimum(rss) / 2^30; digits = 2)) GiB, max $(round(maximum(rss) / 2^30; digits = 2)) GiB")

say("\n== checks ==")
say("thread rows == process rows (exactly): $(res_thr == res_proc)")
sp_thr = N * t_one / t_thr
sp_proc = N * t_one / t_proc
say("\nSUMMARY threads_speedup=$(round(sp_thr; digits = 2)) procs_speedup=$(round(sp_proc; digits = 2)) of N=$N")
verdict = if sp_thr < 0.5N && sp_proc > 0.7N
    "HYPOTHESIS SUPPORTED: threads scale badly, processes scale well"
elseif sp_thr >= 0.7N
    "HYPOTHESIS FALSIFIED: threads also scale well; the cause of the 1-core run is something else"
else
    "INCONCLUSIVE: threads $(round(sp_thr; digits = 2)), processes $(round(sp_proc; digits = 2)) (of $N)"
end
say("VERDICT: " * verdict)
