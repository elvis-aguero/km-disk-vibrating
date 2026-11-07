function [dr, Delta, IntMat, nr, r] = grid_domain(g)
% GRID_DOMAIN  Domain operators of a grid struct. Mirrors build_domain(grid) in julia/src/domain.jl.
%   dr = dr0 (spacing in the disk region), Delta = radial Laplacian (full nr x nr; tridiagonal for
%   laplacian_order 2, banded for 4 and 6), IntMat = pressure-integral matrix, nr, r = node radii 0..nr-1.
%   A uniform grid with order 2 returns exactly domainMaker(D, sr) (same code path, bit for bit).
    if g.uniform && g.laplacian_order == 2
        [dr, Delta, IntMat] = domainMaker(2 * g.R_req, g.sr);
        nr = size(Delta, 1);
        if nr ~= g.nr
            error('grid_domain:nr', 'internal error: uniform grid nr=%d differs from domainMaker nr=%d', g.nr, nr);
        end
        r = g.r(1:nr);
        return;
    end
    if g.laplacian_order == 2
        Delta = laplacian_conservative(g);
    else
        Delta = laplacian_banded(g, g.laplacian_order);
    end
    IntMat = pressure_integral_matrix(g.sr, g.dr0);
    dr = g.dr0;  nr = g.nr;  r = g.r(1:nr);
end

function Delta = laplacian_conservative(g)
% Conservative second-order tridiagonal Laplacian (Julia _laplacian_conservative, docs/graded_grid.md section 3):
%   Delta g_i = 1/(r_i s_i) [ c_{i+1/2}(g_{i+1}-g_i) - c_{i-1/2}(g_i-g_{i-1}) ],  c = r/s at half nodes.
    nr = g.nr;  dr0 = g.dr0;  xin = g.xin;
    d = -2.0 / dr0^2 * ones(nr, 1);
    du = zeros(nr - 1, 1);  dl = zeros(nr - 1, 1);
    d(1) = -4.0 / dr0^2;  du(1) = 4.0 / dr0^2;
    j = (0:nr - 1).';
    c = grid_r_half(g, j) ./ grid_spacing_at(g, j + 0.5);      % c(j+1) = (r/s) at xi = j + 1/2
    for i = 1:nr - 1                                           % node index (0-based)
        if i < xin
            if i <= nr - 2, du(i + 1) = (1.0 / dr0^2) * (1.0 + 1.0 / (2.0 * i)); end
            dl(i) = (1.0 / dr0^2) * (1.0 - 1.0 / (2.0 * i));
        else
            f = 1.0 / (g.r(i + 1) * g.s(i + 1));
            d(i + 1) = -f * (c(i + 1) + c(i));
            dl(i) = f * c(i);
            if i <= nr - 2, du(i + 1) = f * c(i + 1); end
        end
    end
    Delta = diag(d) + diag(dl, -1) + diag(du, 1);
end

function L = laplacian_banded(g, order)
% 4th/6th-order mapped Laplacian: g_xixi/s^2 + (1/(r s) - s'/s^3) g_xi (Julia _laplacian_banded), dense.
    if order == 4
        c2 = [-1 16 -30 16 -1] / 12;  c1 = [1 -8 0 8 -1] / 12;
    else
        c2 = [2 -27 270 -490 270 -27 2] / 180;  c1 = [-1 9 -45 0 45 -9 1] / 60;
    end
    h = (numel(c2) - 1) / 2;
    nr = g.nr;
    L = zeros(nr, nr);
    for k = -h:h
        col = abs(k);
        if col < nr, L(1, col + 1) = L(1, col + 1) + 2 * c2(k + h + 1) / g.s(1)^2; end
    end
    for i = 1:nr - 1
        s = g.s(i + 1);
        A2 = 1 / s^2;
        A1 = 1 / (g.r(i + 1) * s) - g.sp(i + 1) / s^3;
        for k = -h:h
            col = abs(i + k);
            if col < nr, L(i + 1, col + 1) = L(i + 1, col + 1) + A2 * c2(k + h + 1) + A1 * c1(k + h + 1); end
        end
    end
end

function IntMat = pressure_integral_matrix(quant, dr)
% Same as the loop in domainMaker.m (legacy pressure-integral matrix; only the contact weights are used elsewhere).
    nlmax = quant + 1;
    IntMat = zeros(2 * nlmax, 2 * nlmax);
    IntMat(1, 1) = 1 / 12;
    for ii = 2:2 * nlmax
        IntMat(ii, 1) = 1 / 3;
        for jj = 2:ii - 1
            IntMat(ii, jj) = 2 * jj - 2;
        end
        IntMat(ii, ii) = 1.5 * ii - 21 / 12;
    end
    IntMat = pi * dr^2 * IntMat;
end
