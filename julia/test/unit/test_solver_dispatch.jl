@testset "solver_dispatch: RAM budgeting and :auto resolution" begin
    @testset "estimated_lu_cache_bytes matches the direct formula" begin
        n = 102  # nr=50 domain
        @test estimated_lu_cache_bytes(40, n) == 40 * (n^2 * sizeof(Float64) + n * sizeof(Int))
        # Scales linearly with stepsPerCycle and quadratically with n.
        @test estimated_lu_cache_bytes(80, n) == 2 * estimated_lu_cache_bytes(40, n)
    end

    @testset "resolve_solver_type" begin
        n = 102
        need = estimated_lu_cache_bytes(40, n)
        @test resolve_solver_type(:auto, 40, n, need + 1) == :lu
        @test resolve_solver_type(:auto, 40, n, need - 1) == :gmres
        # Explicit requests pass through regardless of budget.
        @test resolve_solver_type(:lu, 40, n, 1) == :lu
        @test resolve_solver_type(:gmres, 40, n, typemax(Int)) == :gmres
    end

    @testset "available_memory_bytes: precedence and the SLURM cgroup-blindness fix" begin
        # This exists because a real cluster run got OOM-killed: Sys.free_memory() reports
        # the whole *node's* free memory under SLURM, not the job's actual --mem cgroup
        # allocation, so :auto sized a 16-thread sweep's LU caches against a budget many
        # times too generous. See available_memory_bytes()'s docstring for the full story.
        old_override = get(ENV, "KMDISK_RAM_BUDGET_BYTES", nothing)
        old_slurm = get(ENV, "SLURM_MEM_PER_NODE", nothing)
        try
            # Explicit override wins outright, with no safety factor applied.
            ENV["KMDISK_RAM_BUDGET_BYTES"] = "12345"
            delete!(ENV, "SLURM_MEM_PER_NODE")
            @test available_memory_bytes() == 12345

            # SLURM's actual per-job allocation, times the safety factor, when no override.
            delete!(ENV, "KMDISK_RAM_BUDGET_BYTES")
            ENV["SLURM_MEM_PER_NODE"] = "32768"  # 32GB in MB, e.g. from `#SBATCH --mem=32G`
            expected = floor(Int, 32768 * 1024 * 1024 * MEMORY_SAFETY_FACTOR)
            @test available_memory_bytes() == expected
            @test available_memory_bytes() < 32768 * 1024 * 1024  # never the raw, un-margined value

            # Falls back to Sys.free_memory() (also safety-margined) only when neither is set.
            delete!(ENV, "SLURM_MEM_PER_NODE")
            @test 0 < available_memory_bytes() <= Sys.free_memory()
        finally
            old_override === nothing ? delete!(ENV, "KMDISK_RAM_BUDGET_BYTES") : (ENV["KMDISK_RAM_BUDGET_BYTES"] = old_override)
            old_slurm === nothing ? delete!(ENV, "SLURM_MEM_PER_NODE") : (ENV["SLURM_MEM_PER_NODE"] = old_slurm)
        end
    end
end
