function n = lu_factor_count(stepsPerCycle, phaseRad, bathDriven)
% LU_FACTOR_COUNT  Number of distinct g_prefactor values in one cycle (= LU factorizations to cache).
%   Mirrors lu_factor_count in julia/src/solver/lu_cache.jl. Disk-forced runs: g == 1, so one.
%   Tolerance 1e-12 (G_DEDUP_ATOL).
    if ~bathDriven
        n = 1;  return;
    end
    distinct = [];
    for k = 0:stepsPerCycle - 1
        cv = cos(2 * pi * (k + 1) / stepsPerCycle + phaseRad);
        if ~any(abs(distinct - cv) <= 1e-12)
            distinct(end + 1) = cv; %#ok<AGROW>
        end
    end
    n = numel(distinct);
end
