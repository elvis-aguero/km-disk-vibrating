function [s, sp] = grid_spacing(xi, P)
% GRID_SPACING  s = r_xi and s' = r_xixi at the (real) index coordinates xi (array, elementwise).
%   Mirrors _grid_spacing / _grid_log_spacing in julia/src/grid.jl (docs/graded_grid.md section 2).
%   P = [dr0 xin lnq ramp lmax beta xd md] (map parameters; xd = Inf if there is no wall zone).
    dr0 = P(1);  xin = P(2);  lnq = P(3);  m = P(4);  lmax = P(5);  beta = P(6);  xd = P(7);  md = P(8);
    t = xi - xin;
    l = zeros(size(xi));  lp = zeros(size(xi));
    k = t > 0;
    if any(k(:))
        x = t(k) / m;
        lg = lnq * m * grid_Pbar(x);
        z = beta * (lg - lmax);
        softplus = max(z, 0) + log1p(exp(-abs(z)));
        l(k) = lg - softplus / beta;
        lp(k) = lnq * grid_rho(x) ./ (1 + exp(z));
    end
    w = xi > xd;                                  % wall zone: log-spacing returns smoothly to 0 over md nodes
    if any(w(:))
        u = (xi(w) - xd) / md;
        D = grid_rho(u);
        drho = zeros(size(u));
        in = u > 0 & u < 1;
        drho(in) = 140 * u(in) .^ 3 .* (1 - u(in)) .^ 3 / md;
        lp(w) = lp(w) .* (1 - D) - l(w) .* drho;
        l(w) = l(w) .* (1 - D);
    end
    s = dr0 * exp(l);
    sp = s .* lp;
end

function r = grid_rho(x)
% smooth step 0 (x<=0) to 1 (x>=1), C^3
    r = x .^ 4 .* (35 - 84 * x + 70 * x .^ 2 - 20 * x .^ 3);
    r(x <= 0) = 0;
    r(x >= 1) = 1;
end

function P = grid_Pbar(x)
% integral of the smooth step
    P = 7 * x .^ 5 - 14 * x .^ 6 + 10 * x .^ 7 - 2.5 * x .^ 8;
    P(x <= 0) = 0;
    hi = x >= 1;
    P(hi) = 0.5 + (x(hi) - 1);
end
