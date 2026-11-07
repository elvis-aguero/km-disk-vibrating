function A = generate_dtn_grid(g, NameValueArgs)
% GENERATE_DTN_GRID  DtN matrix of a grid struct (uniform_grid / graded_grid), physical units (disk radii).
%   Mirrors generate_dtn(grid::RadialGrid) and the non-uniform part of julia/src/dtn/dtn_single_layer.jl
%   (docs/graded_grid.md section 4). A uniform grid calls generate_dtn_single_layer(nr, 2R) unchanged;
%   'forceGeneral', true runs the non-uniform routine on it instead (test). 'method', 'polar' needs a uniform grid.
%   A = generate_dtn_grid(g, 'method', 'single_layer', 'order', 6, 'forceGeneral', false)
arguments
    g (1, 1) struct
    NameValueArgs.method (1, 1) string {mustBeMember(NameValueArgs.method, ["single_layer", "polar"])} = "single_layer"
    NameValueArgs.order (1, 1) double {mustBeMember(NameValueArgs.order, [4 6])} = 6
    NameValueArgs.forceGeneral (1, 1) logical = false
end
if NameValueArgs.method == "polar"
    if ~g.uniform
        error('generate_dtn_grid:polar', 'method :polar needs a uniform grid');
    end
    A = generate_dtn(g.nr, 2 * g.R_req, 'method', 'polar');
    return;
end
order = NameValueArgs.order;
if g.uniform && ~NameValueArgs.forceGeneral
    A = generate_dtn_single_layer(g.nr, 2 * g.R_req, 'order', order);
    return;
end
nr = g.nr;
if nr < 8
    error('generate_dtn_grid:nr', 'generate_dtn_single_layer requires nr >= 8, got nr=%d', nr);
end
G = single_layer_matrix_grid(g);                   % (nr+3) x nr, rows nodes 0..nr+2
if order == 4
    c2 = [-1 16 -30 16 -1] / 12;  c1 = [1 -8 0 8 -1] / 12;  h = 2;
else
    c2 = [2 -27 270 -490 270 -27 2] / 180;  c1 = [-1 9 -45 0 45 -9 1] / 60;  h = 3;
end
A = zeros(nr, nr);
acc0 = zeros(1, nr);                               % node 0: Delta g(0) = 2 g_xixi / s^2, even symmetry
for k = -h:h
    acc0 = acc0 + c2(k + h + 1) * G(abs(k) + 1, :);
end
A(1, :) = -2 * acc0 / g.s(1)^2;
m = (1:nr - 1).';
s = g.s(m + 1);
A2 = 1 ./ s.^2;
A1 = 1 ./ (g.r(m + 1) .* s) - g.sp(m + 1) ./ s.^3;
acc = zeros(nr - 1, nr);
for k = -h:h
    coef = A2 * c2(k + h + 1) + A1 * c1(k + h + 1);
    acc = acc + coef .* G(abs(m + k) + 1, :);
end
A(2:nr, :) = -acc;
end

function G = single_layer_matrix_grid(g)
% Matrix of S on nodal values, rows nodes 0..nr+2, columns nodes 0..nr-1 (Julia _sl_single_layer_matrix_grid).
nr = g.nr;  r = g.r;
nq = nr - 1;
[xs, ws] = gauss_legendre_24();
u = (xs(:).' + 1) / 2;
w = ws(:).' / 2;
q = (0:nq - 1).';
a = r(q + 1);                                       % panel left ends r_q
hq = r(q + 2) - a;                                  % panel widths
LC = cell(1, 4);                                    % LC{p}(q+1, :) = coefficients t^0..t^3 of the cubic for node q + p - 2
for p = 1:4
    LC{p} = zeros(nq, 4);
end
for qq = 0:nq - 1
    if qq == 0, nm1 = -r(2); else, nm1 = r(qq); end % node q-1 (odd extension about the axis)
    x = [(nm1 - a(qq + 1)) / hq(qq + 1), 0.0, (r(qq + 2) - a(qq + 1)) / hq(qq + 1), (r(qq + 3) - a(qq + 1)) / hq(qq + 1)];
    for p = 1:4
        c = 1;  den = 1;
        for j = 1:4
            if j == p, continue; end
            c = conv(c, [-x(j) 1]);
            den = den * (x(p) - x(j));
        end
        LC{p}(qq + 1, :) = c / den;
    end
end
G = zeros(nr + 3, nr);
for mm = 0:nr + 2
    rm = r(mm + 1);
    T = repmat(u, nq, 1);  DT = ones(nq, 24);
    Dm = a + hq .* T - rm;                          % r' - r on regular panels
    if mm <= nr - 2                                 % panel q = mm: singularity at the left end, r = a
        T(mm + 1, :) = u .^ 4;  DT(mm + 1, :) = 4 * u .^ 3;  Dm(mm + 1, :) = hq(mm + 1) * u .^ 4;
    end
    if mm >= 1 && mm - 1 <= nr - 2                  % panel q = mm-1: singularity at the right end, r = a + h
        v = (1 - u) .^ 4;
        T(mm, :) = 1 - v;  DT(mm, :) = 4 * (1 - u) .^ 3;  Dm(mm, :) = -hq(mm) * v;
    end
    rp = a + hq .* T;
    kern = (2 / pi) * rp ./ (rm + rp) .* ellipK_b(abs(Dm) ./ (rm + rp));
    wt = w .* DT .* kern .* hq;
    I = [sum(wt, 2), sum(wt .* T, 2), sum(wt .* T .^ 2, 2), sum(wt .* T .^ 3, 2)];
    row = zeros(nr, 1);
    for p = 1:4
        col = abs(q + (p - 2)) + 1;                 % even extension about r = 0
        val = sum(LC{p} .* I, 2);
        ok = col <= nr;                             % phi = 0 beyond the last node
        row = row + accumarray(col(ok), val(ok), [nr 1]);
    end
    G(mm + 1, :) = row.';
end
end

function K = ellipK_b(b)
% K(m) by the arithmetic-geometric mean, b = sqrt(1-m) in (0,1] (elementwise), as Julia _ellipK_b.
a = ones(size(b));
act = true(size(b));
for it = 1:40
    an = (a(act) + b(act)) / 2;
    bn = sqrt(a(act) .* b(act));
    a(act) = an;  b(act) = bn;
    done = abs(a - b) <= 1e-15 * a;
    act = act & ~done;
    if ~any(act(:)), break; end
end
K = pi ./ (2 * a);
end

function [x, w] = gauss_legendre_24()
n = 24;
k = 1:n - 1;
J = diag(k ./ sqrt(4 * k .^ 2 - 1), 1);
J = J + J.';
[V, E] = eig(J);
[x, idx] = sort(diag(E));
w = 2 * V(1, idx) .^ 2;
x = x(:);  w = w(:);
end
