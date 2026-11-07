function I = grid_integrate_s(a, b, P)
% GRID_INTEGRATE_S  Integral of s over [a, b] (elementwise on arrays a, b), 10-point Gauss-Legendre.
%   Mirrors _grid_integrate_s in julia/src/grid.jl.
    persistent xg wg
    if isempty(xg)
        n = 10;  k = 1:n - 1;
        J = diag(k ./ sqrt(4 * k .^ 2 - 1), 1);  J = J + J.';
        [V, E] = eig(J);
        [xg, idx] = sort(diag(E));
        wg = 2 * V(1, idx(:)') .^ 2;  wg = wg(:);
    end
    a = a(:);  b = b(:);
    h = (b - a) / 2;
    c = (a + b) / 2;
    X = c + h .* xg.';                              % n_intervals x 10
    S = reshape(grid_spacing(X(:), P), size(X));
    I = h .* (S * wg);
end
