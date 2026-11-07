function P = grid_params(g)
% GRID_PARAMS  Map parameters [dr0 xin lnq ramp lmax beta xd md] of a graded grid (Julia _grid_params).
    P = [g.dr0, g.xin, log(g.q), g.ramp, log(g.dr_max / g.dr0), 1 / g.smooth, g.xd, g.md];
end
