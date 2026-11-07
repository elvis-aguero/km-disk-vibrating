function A = generate_dtn_single_layer(nr, D, NameValueArgs)
% GENERATE_DTN_SINGLE_LAYER  Dense DtN matrix by the single-layer construction, O(nr^2) kernel evaluations.
%   Mirrors julia/src/dtn/dtn_single_layer.jl (derivation: docs/dtn_single_layer.md).
%   A = generate_dtn_single_layer(nr, D, 'order', 6)   order = 4 or 6 (finite-difference order of the Laplacian).
%   Node spacing is 1 inside; the result is scaled by 1/dr, dr = D/(2 nr), as in Julia (same units as DTNVectorized).
%   Elliptic K uses the AGM on b = |r-r'|/(r+r') (as Julia), NOT ellipke(m): m = 1 - b^2 rounds to 1 on the singular panels.
arguments
    nr (1, 1) double {mustBeInteger}
    D (1, 1) double
    NameValueArgs.order (1, 1) double {mustBeMember(NameValueArgs.order, [4 6])} = 6
end
order = NameValueArgs.order;
if nr < 8
    error('generate_dtn_single_layer:nr', 'generate_dtn_single_layer requires nr >= 8, got nr=%d', nr);
end
G = sl_single_layer_matrix(nr);                    % (nr+3) x nr, rows = nodes 0..nr+2
A = zeros(nr, nr);
if order == 4
    A(1, :) = -(-G(3, :) + 16 * G(2, :) - 15 * G(1, :)) / 3;          % m = 0: Delta g(0) = 2 g''(0)
    c2 = [-1 16 -30 16 -1] / 12;  c1 = [1 -8 0 8 -1] / 12;  h = 2;
else
    A(1, :) = -(4 * G(4, :) - 54 * G(3, :) + 540 * G(2, :) - 490 * G(1, :)) / 90;
    c2 = [2 -27 270 -490 270 -27 2] / 180;  c1 = [-1 9 -45 0 45 -9 1] / 60;  h = 3;
end
m = (1:nr - 1).';                                   % nodes 1..nr-1 (rows 2..nr of A)
acc = zeros(nr - 1, nr);
for k = -h:h
    coef = c2(k + h + 1) + c1(k + h + 1) ./ m;      % (nr-1) x 1
    acc = acc + coef .* G(abs(m + k) + 1, :);       % row of node |m+k|
end
A(2:nr, :) = -acc;
A = A ./ (D / (2 * nr));
end

function G = sl_single_layer_matrix(nr)
% Matrix of S on nodal values: rows nodes 0..nr+2, columns nodes 0..nr-1 (Julia _sl_single_layer_matrix).
[xs, ws] = gauss_legendre_24();
u = (xs(:).' + 1) / 2;                              % 1 x 24
w = ws(:).' / 2;
q = (0:nr - 2).';                                   % panel left ends, (nr-1) x 1
nq = nr - 1;
LC = zeros(4, 4);                                   % column p: coefficients t^0..t^3 of node tp = p-2
for p = 1:4
    LC(:, p) = lagrange_coeffs(p - 2);
end
G = zeros(nr + 3, nr);
for m = 0:nr + 2
    S = repmat(u, nq, 1);  DS = ones(nq, 24);
    if m <= nr - 2                                  % panel q = m has its log singularity at the left end
        S(m + 1, :) = u .^ 4;  DS(m + 1, :) = 4 * u .^ 3;
    end
    if m >= 1 && m - 1 <= nr - 2                    % panel q = m-1: singularity at the right end
        S(m, :) = 1 - (1 - u) .^ 4;  DS(m, :) = 4 * (1 - u) .^ 3;
    end
    rp = q + S;                                     % r' = a + s (same floating ops as Julia)
    kern = (2 / pi) * rp ./ (m + rp) .* ellipK_b(abs(m - rp) ./ (m + rp));
    wt = w .* DS .* kern;
    I = [sum(wt, 2), sum(wt .* S, 2), sum(wt .* S .^ 2, 2), sum(wt .* S .^ 3, 2)];   % nq x 4
    row = zeros(nr, 1);
    for p = 1:4
        col = abs(q + (p - 2)) + 1;                 % even extension about r = 0
        val = I * LC(:, p);
        ok = col <= nr;                             % phi = 0 beyond the last node
        row = row + accumarray(col(ok), val(ok), [nr 1]);
    end
    G(m + 1, :) = row.';
end
end

function c = lagrange_coeffs(tp)
% Coefficients (t^0..t^3) of the Lagrange cubic for node tp among t = -1, 0, 1, 2.
c = 1;  den = 1;
for tj = [-1 0 1 2]
    if tj == tp, continue; end
    c = conv(c, [-tj 1]);
    den = den * (tp - tj);
end
c = c(:) / den;
end

function K = ellipK_b(b)
% Complete elliptic integral K(m), b = sqrt(1-m) in (0,1], arithmetic-geometric mean, elementwise (Julia _ellipK_b).
a = ones(size(b));
act = true(size(b));
for it = 1:40
    an = (a(act) + b(act)) / 2;
    bn = sqrt(a(act) .* b(act));
    a(act) = an;  b(act) = bn;
    done = abs(a - b) <= 1e-15 * a;
    act = act & ~done;
    if ~any(act), break; end
end
K = pi ./ (2 * a);
end

function [x, w] = gauss_legendre_24()
% Gauss-Legendre nodes/weights on [-1,1] by Golub-Welsch (Julia _gauss_legendre(24)).
n = 24;
k = 1:n - 1;
J = diag(k ./ sqrt(4 * k .^ 2 - 1), 1);
J = J + J.';
[V, E] = eig(J);
[x, idx] = sort(diag(E));
w = 2 * V(1, idx) .^ 2;
x = x(:);  w = w(:);
end
