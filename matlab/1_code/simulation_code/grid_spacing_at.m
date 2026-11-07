function s = grid_spacing_at(g, xi)
% GRID_SPACING_AT  s = r_xi at real indices xi, for either kind of grid (Julia _grid_spacing_at, first output).
    if g.uniform
        s = g.dr0 * ones(size(xi));
    else
        s = grid_spacing(xi, grid_params(g));
    end
end
