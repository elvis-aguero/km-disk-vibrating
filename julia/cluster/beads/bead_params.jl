#= Fluid and bead inputs for the bead campaign.

Source: fluid_and_bead_parameters.pdf (2 Nov 2026), converted to CGS because
`SimulationParams` uses cm, g, dyne/cm and Stokes. Kinematic viscosity 7.9 cSt = 0.079 St.

`R`  : actual bead radius (PDF table 2).
`rc` : contact-line radius from the full Laplace-Young solution (PDF table 3; theta = 160
       deg assumed there). A sweep uses one of the two as `diskRadius`; the disk mass is
       always the full sphere mass, (4/3) pi R^3 rho_s. =#

const BEAD_FLUID = (density = 1.148,          # g/cm^3
                    surface_tension = 68.0,   # dyne/cm
                    viscosity = 0.079,        # Stokes (cm^2/s)
                    g = 981.0)                # cm/s^2

const BEADS = Dict(
    "ptfe"  => (name = "PTFE (reference)", R = 0.07900, rc = 0.04342, rho_s = 2.20),
    "al32"  => (name = "Aluminum 1/32in",  R = 0.03969, rc = 0.01713, rho_s = 2.79),
    "al16"  => (name = "Aluminum 1/16in",  R = 0.07938, rc = 0.04663, rho_s = 2.79),
    "nylon" => (name = "Nylon 1/16in",     R = 0.07938, rc = 0.03743, rho_s = 1.14),
)

bead_mass_g(b) = 4 / 3 * pi * b.R^3 * b.rho_s

# 5-90 Hz: sparse (5 Hz steps) away from 60 Hz, 1 Hz steps over 50-70 Hz. 34 frequencies.
const DEFAULT_FREQS = collect(Float64, vcat(5.0:5.0:45.0, 50.0:1.0:70.0, 75.0:5.0:90.0))
