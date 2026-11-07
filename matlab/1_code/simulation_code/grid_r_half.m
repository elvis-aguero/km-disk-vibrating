function rh = grid_r_half(g, j)
% GRID_R_HALF  Radius at the half nodes xi = j + 1/2 (j: array of integers). Mirrors _grid_r_half in julia/src/grid.jl.
    j = j(:);
    rh = g.dr0 * (j + 0.5);
    if g.uniform, return; end
    k = j >= g.xin;
    if any(k)
        P = grid_params(g);
        rh(k) = g.r(j(k) + 1) + grid_integrate_s(j(k), j(k) + 0.5, P);
    end
end
