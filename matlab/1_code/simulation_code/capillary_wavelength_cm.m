function lambda = capillary_wavelength_cm(freqHz, density, tension, g)
% CAPILLARY_WAVELENGTH_CM  Wavelength (cm) of the deep-water gravity-capillary wave of frequency freqHz:
%   omega^2 = g k + (sigma/rho) k^3. CGS. Mirrors capillary_wavelength_cm in julia/src/grid.jl (bisection, 200 steps).
    w2 = (2 * pi * freqHz)^2;
    sg = tension / density;
    f = @(k) g * k + sg * k^3 - w2;
    lo = 0.0;  hi = 1.0;
    while f(hi) < 0
        hi = hi * 2;
    end
    for it = 1:200
        mid = (lo + hi) / 2;
        if f(mid) < 0
            lo = mid;
        else
            hi = mid;
        end
    end
    lambda = 2 * pi / ((lo + hi) / 2);
end
