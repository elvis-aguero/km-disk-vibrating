function b = ram_budget_bytes()
% RAM_BUDGET_BYTES  Memory budget for the solver decision. Mirrors available_memory_bytes in
%   julia/src/solver/solver_dispatch.jl: KMDISK_RAM_BUDGET_BYTES (as is) > SLURM_MEM_PER_NODE (MB) * 0.7 >
%   getAvailableRAM() * 0.7 (not cgroup aware; under Slurm the second rule is the one that applies).
    ov = getenv('KMDISK_RAM_BUDGET_BYTES');
    if ~isempty(ov)
        b = str2double(ov);  return;
    end
    sl = getenv('SLURM_MEM_PER_NODE');
    if ~isempty(sl) && ~isnan(str2double(sl))
        b = 0.7 * str2double(sl) * 1024 * 1024;
    else
        b = 0.7 * getAvailableRAM();
    end
end
