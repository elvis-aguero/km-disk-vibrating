function g = graded_grid(aOverDr0, R, NameValueArgs)
% GRADED_GRID  Graded radial grid (fine near the disk, coarse far away). Mirrors graded_grid in julia/src/grid.jl;
%   formulas and criteria: docs/graded_grid.md. Lengths in disk radii.
%   g = graded_grid(sr, R, 'drMax', dmax, 'q', 1.04, 'rIn', 2, 'ramp', 20, 'smooth', 0.1, ...
%                   'laplacianOrder', 2, 'wallZone', 0)
%   Identity map (dr0 = 1/sr) for r <= rIn (needs rIn >= 1 + 5 dr0); then growth q per node capped at drMax.
%   If drMax <= dr0, q <= 1 or R <= rIn the result is uniform_grid. The wall is the first node with r >= R.
%   wallZone > 0: the last wallZone disk radii keep the fine spacing dr0.
arguments
    aOverDr0 (1, 1) double
    R (1, 1) double
    NameValueArgs.drMax (1, 1) double
    NameValueArgs.q (1, 1) double = 1.04
    NameValueArgs.rIn (1, 1) double = 2.0
    NameValueArgs.ramp (1, 1) double = 20.0
    NameValueArgs.smooth (1, 1) double = 0.1
    NameValueArgs.laplacianOrder (1, 1) double = 2
    NameValueArgs.wallZone (1, 1) double = 0.0
end
if ~isfield(NameValueArgs, 'drMax')
    error('graded_grid:drMax', 'graded_grid needs the name-value argument drMax');
end
order = NameValueArgs.laplacianOrder;
if ~any(order == [2 4 6])
    error('graded_grid:order', 'laplacian_order must be 2, 4 or 6, got %d', order);
end
sr = aOverDr0;  Rf = R;  qf = NameValueArgs.q;  dmax = NameValueArgs.drMax;  rin = NameValueArgs.rIn;
rampf = NameValueArgs.ramp;  smf = NameValueArgs.smooth;  wz = NameValueArgs.wallZone;
if ~(wz >= 0)
    error('graded_grid:wallZone', 'wall_zone must be >= 0');
end
dr0 = 1 / sr;
if dmax <= dr0 * (1 + 1e-9) || qf <= 1 || Rf <= rin
    g = uniform_grid(sr, Rf);
    g.laplacian_order = order;
    return;
end
if ~(rin >= 1 + 5 * dr0)
    error('graded_grid:rIn', 'r_in must be at least 1 + 5 dr0 = %.15g (contact nodes stay in the identity region), got %.15g', ...
        1 + 5 * dr0, rin);
end
if ~(rampf > 0 && smf > 0)
    error('graded_grid:rampSmooth', 'ramp and smooth must be positive');
end
xin = ceil(rin / dr0);
lnq = log(qf);  lmax = log(dmax / dr0);  beta = 1 / smf;
P = [dr0, xin, lnq, rampf, lmax, beta, Inf, 1.0];
if wz > 0
    md = max(rampf, ceil(lmax / lnq));
    if md ~= round(md)
        error('graded_grid:md', 'wall-zone return length md = %g must be a whole number of nodes (use an integer ramp)', md);
    end
    [~, rup] = grid_march(P, Rf);                % grid without the zone: where the spacing could start to return
    best = 0;
    for xd = (xin + ceil(rampf)):(numel(rup) - 2)
        Pd = [dr0, xin, lnq, rampf, lmax, beta, xd, md];
        ints = grid_integrate_s((xd:xd + md - 1).', (xd + 1:xd + md).', Pd);
        rEndAll = cumsum([rup(xd + 1); ints]);       % same left-to-right sum as Julia
        if rEndAll(end) > Rf - wz
            break;
        end
        best = xd;
    end
    if ~(best > 0)
        error('graded_grid:wallZoneFit', 'wall_zone = %g does not fit: the spacing cannot return to dr0 before R - wall_zone = %g', ...
            wz, Rf - wz);
    end
    P(7) = best;  P(8) = md;
end
[nr, r] = grid_march(P, Rf);
[s, sp] = grid_spacing((0:nr + 2).', P);
g = struct('nr', nr, 'r', r, 's', s, 'sp', sp, 'dr0', dr0, 'sr', sr, 'R_req', Rf, 'uniform', false, 'q', qf, ...
    'dr_max', dmax, 'r_in', rin, 'xin', xin, 'ramp', rampf, 'smooth', smf, 'laplacian_order', order, ...
    'wall_zone', wz, 'xd', P(7), 'md', P(8));
end

function [nr, r] = grid_march(P, Rf)
% nodes r_0, r_1, ... until r >= Rf (the wall, index nr), then 3 more. Mirrors _grid_march (chunked, same sums).
    dr0 = P(1);  xin = P(2);
    r = dr0 * (0:xin).';                             % identity region r_i = dr0 i (i <= xin)
    i0 = xin;
    chunk = 512;
    while true
        idx = find(r >= Rf, 1);                      % 1-based: node index nr = idx - 1
        if ~isempty(idx) && numel(r) >= idx + 3
            break;
        end
        if i0 > 2e7
            error('graded_grid:tooManyNodes', 'graded_grid: more than 2e7 nodes (check R and dr_max)');
        end
        ii = (i0 + 1:i0 + chunk).';
        ints = grid_integrate_s(ii - 1, ii, P);
        rn = cumsum([r(end); ints]);                 % (r_end + I1) + I2 ... as the Julia loop
        r = [r; rn(2:end)];
        i0 = i0 + chunk;
    end
    nr = idx - 1;
    r = r(1:nr + 4);
end
