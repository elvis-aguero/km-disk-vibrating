function [L, U, P] = computeLU(PC_min, g_pf, dt, nr, cPoints, dr, SF)
% COMPUTELU  Builds the system matrix and returns its LU factorization.
%   Called via parfeval for async pre-computation of LU libraries.
%   PC_min is a minimal struct with only the fields needed by buildSystemMatrix
%   (weber, reynolds, froude, laplacian, DTN, pressure_integral, obj_mass),
%   keeping serialization overhead small.

    Mat = buildSystemMatrix(PC_min, g_pf, dt, nr, cPoints, dr, SF);
    [L, U, P] = lu(Mat);
end
