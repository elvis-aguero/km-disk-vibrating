function sol = schur_solve(sys, F, g, rhs)
% SCHUR_SOLVE  One step with the reduced system and a factorization F of S(g) = S0 + g G1 (a decomposition).
%   Returns sol in the full layout of advance_one_step.m: [eta_out; phi; p; v; z]. Mirrors schur_solve! in
%   julia/src/solver/schur_solver.jl.
    nr = sys.nr;  c = sys.cPoints;  dt = sys.dt;  no = nr - c;
    u = sys.dB \ rhs(nr + 1:2 * nr);                     % B^-1 r_bot
    b = zeros(nr + 1, 1);
    b(1:nr) = rhs(1:nr) + dt * (sys.DTN * u);
    r_f = rhs(2 * nr + 1);  r_k = rhs(2 * nr + 2);
    b(nr + 1) = r_f + r_k / dt;
    y = F \ b;                                           % [eta_out; z; p]
    z = y(no + 1);
    q = [z * ones(c, 1); y(1:no)];                       % eta
    t = sys.K0 * q + (dt * g / sys.Fr) * q;
    t(1:c) = t(1:c) + dt * y(no + 2:nr + 1);
    t = sys.dB \ t;
    sol = zeros(2 * nr + 2, 1);
    sol(1:no) = y(1:no);
    sol(no + 1:no + nr) = u - t;                         % phi
    sol(2 * nr - c + 1:2 * nr) = y(no + 2:nr + 1);       % pressure
    sol(2 * nr + 1) = (z - r_k) / dt;
    sol(2 * nr + 2) = z;
end
