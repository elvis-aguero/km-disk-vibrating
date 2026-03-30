function Mat = buildSystemMatrix(PROBLEM_CONSTANTS, g_prefactor, dt, nr, cPoints, dr, SF)
% BUILDSYSTEMMATRIX  Assembles the full linear system matrix for one time step.
%   Standalone file so it is reachable by parfeval workers during async
%   LU pre-computation in solve_motion.

    We        = PROBLEM_CONSTANTS.weber;
    Re        = PROBLEM_CONSTANTS.reynolds;
    Delta     = PROBLEM_CONSTANTS.laplacian;
    DTN       = PROBLEM_CONSTANTS.DTN;
    Fr        = PROBLEM_CONSTANTS.froude;
    pIntegral = PROBLEM_CONSTANTS.pressure_integral;
    Ma        = PROBLEM_CONSTANTS.obj_mass;

    Sist = [[eye(nr)-dt*2*Delta/Re, -dt*DTN]; ...
            [dt*(eye(nr)/Fr * g_prefactor - Delta/We), eye(nr)-dt*2*Delta/Re]];

    Mat = [[Sist(:,(cPoints+1):2*nr), ...
            [zeros(nr,cPoints); dt*eye(cPoints); zeros(nr-cPoints,cPoints)], ...
            zeros(2*nr,1), Sist(:,1:cPoints)*ones(cPoints,1)]; ...
           [-SF*dt/dr, zeros(1,2*nr-cPoints-1), -dt*pIntegral(1:cPoints)/Ma, 1, SF*dt/dr]; ...
           [zeros(1,2*nr-cPoints), -zeros(1,cPoints), -dt, 1]];
end
