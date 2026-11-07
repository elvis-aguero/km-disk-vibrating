function g = graded_grid_for(freqHz, NameValueArgs)
% GRADED_GRID_FOR  Graded grid for the bath frequency freqHz: cap spacing dr_max = lambda(f)/ppw (disk radii).
%   Mirrors graded_grid_for in julia/src/parameters.jl. The fluid and domain arguments take the defaults of
%   solve_motion (cgs); pass the same values you pass to solve_motion. Warns if the a-priori wall echo exceeds 0.01.
%   g = graded_grid_for(40, 'bathDensity', 1.148, 'bathSurfaceTension', 68, 'bathViscosity', 0.079, ...
%                       'g', 981, 'diskRadius', 0.07938, 'bathDiameter', 700, 'spatialResolution', 25)
%   Defaults: ppw 32, q 1.04, rIn 2, ramp 20, smooth 0.1, laplacianOrder 2, wallZone 2.
arguments
    freqHz (1, 1) double
    NameValueArgs.bathDensity (1, 1) double = 1
    NameValueArgs.bathSurfaceTension (1, 1) double = 72.20
    NameValueArgs.bathViscosity (1, 1) double = 0.978e-2
    NameValueArgs.g (1, 1) double = 981
    NameValueArgs.diskRadius (1, 1) double = 0.2
    NameValueArgs.bathDiameter (1, 1) double = 100
    NameValueArgs.spatialResolution (1, 1) double = 50
    NameValueArgs.ppw (1, 1) double = 32.0
    NameValueArgs.q (1, 1) double = 1.04
    NameValueArgs.rIn (1, 1) double = 2.0
    NameValueArgs.ramp (1, 1) double = 20.0
    NameValueArgs.smooth (1, 1) double = 0.1
    NameValueArgs.laplacianOrder (1, 1) double = 2
    NameValueArgs.wallZone (1, 1) double = 2.0
end
a = NameValueArgs;
lambda_cm = capillary_wavelength_cm(freqHz, a.bathDensity, a.bathSurfaceTension, a.g);
dr_max = lambda_cm / a.ppw / a.diskRadius;
% a-priori wall-echo amplitude exp(-2 R_domain / L_e)
omega = 2 * pi * freqHz;
k = 2 * pi / lambda_cm;
cg = (a.g + 3 * (a.bathSurfaceTension / a.bathDensity) * k^2) / (2 * omega);
echo = exp(-2 * (a.bathDiameter / 2 * a.diskRadius) / (cg / (2 * a.bathViscosity * k^2)));
if echo > 0.01
    warning('graded_grid_for:echo', ['graded_grid_for: predicted wall echo %.3g exceeds 0.01 at %g Hz for R = %g disk radii; ' ...
        'the graded grid is validated only for echo-negligible domains (docs/graded_grid.md)'], echo, freqHz, a.bathDiameter / 2);
end
g = graded_grid(a.spatialResolution, a.bathDiameter / 2, 'q', a.q, 'drMax', dr_max, 'rIn', a.rIn, 'ramp', a.ramp, ...
    'smooth', a.smooth, 'laplacianOrder', a.laplacianOrder, 'wallZone', a.wallZone);
end
