function sys = schur_system(PC, dt)
% SCHUR_SYSTEM  Reduced (nr+1)x(nr+1) system of the Schur step solver. Mirrors julia/src/solver/schur_solver.jl
%   (build_schur_system). The bath potential phi is eliminated with the tridiagonal block B = I - 2 dt Delta/Re:
%       S(g) y = b,   S(g) = S0 + g G1,   y = [eta_out; z; p],   Y = N B^-1.
%   dt is the IMPLICIT coefficient (dt for euler, 2dt/3 for bdf2). Needs a tridiagonal Laplacian.
%   Fields: nr, cPoints, dt, Fr, dB (decomposition of B), K0 (sparse -(dt/We) Delta), S0, G1, DTN.
    nr = PC.nr;  c = PC.contact_points;  no = nr - c;  m = nr + 1;
    Re = PC.reynolds;  Fr = PC.froude;  We = PC.weber;  Ma = PC.obj_mass;
    if ~is_tridiagonal(PC.laplacian)
        error('schur_system:notTridiagonal', ...
            'solverType "schur" needs a tridiagonal Laplacian (uniform grid, or graded grid with laplacian order 2).');
    end
    Ds = sparse(PC.laplacian);
    N = PC.DTN;
    B = speye(nr) - (2 * dt / Re) * Ds;
    dB = decomposition(B, 'lu');
    Y = (B.' \ N.').';                                  % Y = N B^-1, tridiagonal solves
    Heta0 = B - (dt^2 / We) * (Y * Ds);                  % nr x nr: A - (dt^2/We) Y Delta
    G1eta = (dt^2 / Fr) * Y;                             % g-part
    S0 = zeros(m, m);  G1 = zeros(m, m);
    S0(1:nr, 1:no) = Heta0(:, c + 1:nr);                 % eta_out columns
    G1(1:nr, 1:no) = G1eta(:, c + 1:nr);
    S0(1:nr, no + 1) = sum(Heta0(:, 1:c), 2);            % z column = sum of the first c eta columns
    G1(1:nr, no + 1) = sum(G1eta(:, 1:c), 2);
    S0(1:nr, no + 2:no + 1 + c) = dt^2 * Y(:, 1:c);      % pressure columns
    S0(m, no + 1) = 1 / dt;                              % (z - r_kin)/dt - dt w.p/Ma = r_force
    w = PC.pressure_integral(1:c);
    S0(m, no + 2:no + 1 + c) = -dt * w(:).' / Ma;
    sys = struct('nr', nr, 'cPoints', c, 'dt', dt, 'Fr', Fr, 'dB', dB, ...
        'K0', -(dt / We) * Ds, 'S0', S0, 'G1', G1, 'DTN', N);
end

function tf = is_tridiagonal(L)
    tf = nnz(triu(L, 2)) == 0 && nnz(tril(L, -2)) == 0;
end
