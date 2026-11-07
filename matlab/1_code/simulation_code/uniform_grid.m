function g = uniform_grid(aOverDr0, R)
% UNIFORM_GRID  Today's uniform radial grid as a grid struct. Mirrors uniform_grid in julia/src/grid.jl.
%   nr = ceil(aOverDr0*R*...) nodes: D = 2R, nr = ceil(D*sr/2), dr = D/(2 nr) (same expressions as domainMaker).
%   Fields: nr, r (nr+4 x 1, nodes 0..nr+3, disk radii; r(nr+1) is the wall), s, sp (nr+3 x 1), dr0, sr, R_req,
%   uniform, q, dr_max, r_in, xin, ramp, smooth, laplacian_order, wall_zone, xd, md.
    sr = double(aOverDr0);
    D = 2 * double(R);
    nr = ceil(D * sr / 2);
    dr = D / (2 * nr);
    g = struct('nr', nr, 'r', dr * (0:nr + 3).', 's', dr * ones(nr + 3, 1), 'sp', zeros(nr + 3, 1), ...
        'dr0', dr, 'sr', sr, 'R_req', double(R), 'uniform', true, 'q', 1.0, 'dr_max', dr, 'r_in', 0.0, ...
        'xin', 0, 'ramp', 0.0, 'smooth', 0.0, 'laplacian_order', 2, 'wall_zone', 0.0, 'xd', Inf, 'md', 1.0);
end
