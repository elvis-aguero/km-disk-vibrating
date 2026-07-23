function [amplitude, phase, offset] = fit_oscillation(t, z, omega)
%FIT_OSCILLATION Least-squares fit of z(t) ~= amplitude*cos(omega*t+phase) + offset.
%   Mirrors julia/src/convergence.jl's fit_oscillation (same M=[cos,sin,1] design
%   matrix, same amplitude=hypot(X1,X2), phase=atan2(-X2,X1) convention) and the
%   fit already used in overlay_validation.py. The constant column absorbs any DC
%   offset (e.g. the static-equilibrium depth), so the recovered amplitude/phase
%   are insensitive to it.

t = t(:);
z = z(:);
if numel(t) ~= numel(z)
    error('fit_oscillation:sizeMismatch', 't and z must have the same length');
end
if numel(t) < 3
    error('fit_oscillation:tooFewSamples', 'fit_oscillation needs at least 3 samples, got %d', numel(t));
end

M = [cos(omega * t), sin(omega * t), ones(numel(t), 1)];
X = M \ z;

amplitude = hypot(X(1), X(2));
phase = atan2(-X(2), X(1));
offset = X(3);

end
