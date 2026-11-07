function [resolved, info] = resolve_solver_type(requested, laplacian, nr, nFactors, budgetBytes)
% RESOLVE_SOLVER_TYPE  "auto" -> "schur" | "lu" | "gmres". Mirrors build_step_solver / resolve_solver_type in
%   julia/src/solver/solver_dispatch.jl. Explicit "lu", "gmres", "schur" pass through unchanged.
%   auto: schur if the Laplacian is tridiagonal and the Schur cache fits in budgetBytes; else lu if the LU cache
%   fits; else gmres. nFactors = lu_factor_count(...) (distinct g values).
%   MATLAB stores [L,U,P] (3 full matrices) per LU factorization, Julia one LU object; the LU estimate keeps
%   MATLAB's 3x. The Schur cache uses decomposition objects: the Julia estimate.
    requested = string(requested);
    n = 2 * nr + 2;
    info = struct('schurBytes', estimated_schur_bytes(nFactors, nr), 'luBytes', nFactors * 3 * n^2 * 8, ...
        'tridiagonal', nnz(triu(laplacian, 2)) == 0 && nnz(tril(laplacian, -2)) == 0);
    if requested ~= "auto"
        resolved = requested;
    elseif info.tridiagonal && info.schurBytes <= budgetBytes
        resolved = "schur";
    elseif info.luBytes <= budgetBytes
        resolved = "lu";
    else
        resolved = "gmres";
    end
end

function b = estimated_schur_bytes(nFactors, nr)
% Julia estimated_schur_bytes: factorizations of (nr+1)^2, S0 and G1, transient Y (+ transpose).
    m = nr + 1;
    b = nFactors * (m^2 * 8 + m * 8) + 2 * m^2 * 8 + 2 * nr^2 * 8;
end
