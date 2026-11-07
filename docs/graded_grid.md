# Graded radial grid: fewer nodes, same physics

Status: experimental, opt-in. The default grid stays uniform (`dr = a/25`, `nr = 8750` in the 25 x 700 domain) and its results do not change (a test asserts bit-identical output). Select a graded grid with `graded_grid_for(params, f)` or `KMDISK_GRADED="ppw,q,order"` in the validation runner. It works with `solverType = :schur` only for `laplacian_order = 2`.

Measured (pilot; criteria in section 5, fixed before the runs). Validated configuration: order 2, `ppw = 32`, `q = 1.04`, `wall_zone = 2` (the defaults of `graded_grid_for`):
* C1, C2, C5 pass. C3 (30 bead rows vs the uniform results; gain within 0.3 %, phase within 0.3 deg) passes 30/30 (largest differences 0.030 % and 0.028 deg); C3b (the 8 good 5 Hz rows) passes 8/8; C4 (30-case experimental validation) gives 0.022163 / 3.0582 deg against the limits 0.022578 / 3.1557 (uniform 0.022135 / 3.0557). The same holds for ppw 24.
* Without the wall zone the grid fails in small, echo-sensitive domains (25 x 200, al32 rc at 52-70 Hz: up to +1.3 % gain and 1.2 deg). With a fine zone of 2 disk radii at the wall the same rows agree to 0.21 % and 0.18 deg (C7c); the zone costs about 120 nodes at R = 100 a. The fourth-order Laplacian fails one C3 row (+0.306 %) and is not used.
* Domain convergence (informational): al16 R at 5 Hz and al32 rc at 52-70 Hz do not change from R = 350 a to 3000 a (and 8000 a for the 5 Hz case, 3226 nodes).
* `graded_grid_for` warns when the a-priori echo prediction exceeds 0.01.

Not validated: any frequency below 15 Hz in the physics comparison (no `good` reference row at 10 Hz), domains larger than 25 x 700, a disk other than the ones in the bead campaign and R = 0.2 cm.

## 1. Why

The solver uses `nr` uniform radial nodes. In the 25 x 700 bead domain `nr = 8750`. The waves need far fewer nodes. With the bead fluid (nu 0.079 St, rho 1.148, sigma 68, g 981, a = 0.07938 cm) a uniform `dr = a/25` puts 115 to 2000 nodes on one wavelength (5 to 90 Hz). A grid that is fine near the disk and coarse far away needs about 220 to 900 nodes (ESTIMATE by the lead, `grid_check.py`; it counts wave resolution only, it is a HYPOTHESIS until section 6 is run). The cost of the dense solver is `O(nr^2)` memory and `O(nr^3)` per factorization, so 10 times fewer nodes is about 100 times less memory.

## 2. The map

Node `i` has an index coordinate `xi = i` (spacing 1). The radius is `r = r(xi)` with `r(0) = 0`, odd in `xi`, and `s(xi) = r_xi(xi) > 0`.

* Inside `r <= r_in` the map is the identity: `s = dr0`, `r_i = i dr0`. So the disk, the contact nodes and `contact_weights` are exactly as today.
* Outside, the log-spacing `l = ln(s/dr0)` grows with slope `ln q` (growth `q` per node, `q` = 1.03 to 1.05), up to the cap `dr_max`:

  ```
  l_g(xi) = ln(q) * m * P(t/m),    t = xi - xi_in,  xi_in = ceil(r_in/dr0)
  l(xi)   = l_g - (1/beta) ln(1 + exp(beta (l_g - l_max)))        (smooth min of l_g and l_max)
  s       = dr0 exp(l),   l_max = ln(dr_max/dr0),   beta = 1/smooth
  ```
  `P` is the integral of the smooth step `rho(x) = x^4 (35 - 84x + 70x^2 - 20x^3)` on `0 <= x <= 1` (then 1), so `P(x) = 7x^5 - 14x^6 + 10x^7 - 2.5x^8` for `x <= 1`, `P(x) = 1/2 + (x - 1)` after. `m = ramp` (default 20 nodes) is the length of the ramp where the growth slope rises from 0 to `ln q`. `l` is 0 for `t <= 0`. The map is `C^4` (the step `rho` is `C^3`). The growth per node is at most `q` (slope of `l` is at most `ln q`), and `s <= dr_max` always (softplus is at least `max(0, .)`).
* `r(xi)` is the integral of `s`, by 10-point Gauss-Legendre on each unit interval. The wall is the first node with `r_nr >= R`. Unknown nodes are `0..nr-1`; node `nr` is the wall where `phi = 0` (as today). Three more nodes beyond `nr` exist for the DtN stencil rows.

API (`julia/src/grid.jl`): `RadialGrid`, `uniform_grid(a_over_dr0, R)` (equal to today's grid), `graded_grid(a_over_dr0, R; q, dr_max, r_in, ramp, smooth, laplacian_order)`, and `graded_grid_for(params, freqHz; ppw, ...)` which sets `dr_max = lambda(f)/ppw` from the dispersion relation `omega^2 = g k + (sigma/rho) k^3`. `SimulationParams` has a new field `grid` (default `nothing` = today's uniform grid).

## 3. The Laplacian on the map

With `g(r) = g(r(xi))`, the radial Laplacian `(1/r) d/dr (r dg/dr)` becomes

```
Delta g = (1/(r r_xi)) d/dxi ( (r/r_xi) dg/dxi )
        = g_xixi / s^2 + A1 g_xi ,     A1 = 1/(r s) - s'/s^3 ,  s' = r_xixi = s l'.
```

(Check: `d/dxi ((r/s) g_xi) = g_xi + (r/s) g_xixi - r s'/s^2 g_xi`; divide by `r s`.)

Two discretizations. Both use `xi`-stencils at integer nodes, with spacing 1.

1. Conservative, second order, TRIDIAGONAL (`laplacian_order = 2`, the default):
   `Delta g_i = (1/(r_i s_i)) [ c_{i+1/2} (g_{i+1} - g_i) - c_{i-1/2} (g_i - g_{i-1}) ]`, `c = r/s` at the half nodes (`r` and `s` at `xi = i +- 1/2` come from the map). Error `O(1)` in `xi` spacing, so it is second order under refinement. It is the form the lead's Schur solver needs.
2. Finite-difference, order 4 or 6 (`laplacian_order = 4, 6`): the 5-point or 7-point stencils for `g_xixi` and `g_xi` of `docs/dtn_single_layer.md`, with `s`, `s'` known in closed form. The matrix has bandwidth 2 or 3 and is NOT tridiagonal. It is stored dense (the sizes are `nr <~ 2000`).

Axis: `r_in >= 1 + 5 dr0` is enforced, so rows `i = 0` and `i < xi_in` are the existing rows (bit-identical to `build_domain(D, sr)` with `dr = dr0`). The first row is `-4/dr0^2, 4/dr0^2` (mirror stencil). The 4th/6th-order rows use even extension `g_{-k} = g_k` (valid because the map is the identity there).
Wall: `g_j = 0` for `j >= nr` in both forms (as today: row `nr-1` has no upper neighbour). The wall rows are lower order. They matter only for the wall echo.

## 4. The DtN matrix on the map

Same construction as `docs/dtn_single_layer.md` (`N = -Delta_H S`, elliptic-integral kernel, 24-point Gauss-Legendre, graded substitution `r' = r_i +- h u^4` on the two singular panels), with three changes.

1. Panels are `[r_q, r_{q+1}]` of width `h_q`. The interpolant of `phi` on a panel is the cubic through the 4 nodes `r_{q-1}, ..., r_{q+2}` (even extension at the axis: `r_{-1} = -r_1`, `phi_{-1} = phi_1`; `phi = 0` beyond the last node). The Lagrange cubics are written in `t = (r' - r_q)/h_q`; the moments `I_n = int K t^n dr'` use the `h_q` factor. Nodes are in disk-radius units. No final `1/dr` rescale.
2. On the singular panels the kernel argument `|r - r'|` is computed from the substitution (`h u^4` or `h (1-u)^4`), not by subtracting two nearby radii.
3. `N = -Delta_H G` uses the mapped stencil of section 3 (order 4 or 6, default 6) on the rows of `G`.

For a uniform grid the code calls the old routine unchanged, so `generate_dtn(grid)` equals `generate_dtn(nr, D)` bit for bit. The general routine also runs on uniform nodes (`force_general = true`) and must agree to `1e-8` of `max|A|`: this tests the generalization.

DtN matrices of graded grids are built in process. There is no file cache and no registry entry (build time is well below one second at `nr` near 500; the job reports it).

### Optional wall zone

`graded_grid(...; wall_zone = w)` keeps the fine spacing `dr0` in the last `w` disk radii. The log-spacing is multiplied by `1 - rho((xi - xd)/md)` for `xi > xd` (`rho` the smooth step above, `md = max(ramp, ln(dr_max/dr0)/ln q)` nodes), so the spacing returns smoothly from the cap to `dr0`; `xd` is the largest integer index such that the return ends at `r <= R - w`. `s' = s l'` with `l' = l_up'(1-D) - l_up D'`, `D' = 140 u^3 (1-u)^3 / md`. `w = 0` (default) is the grid without a zone, bit for bit. Reason (HYPOTHESIS, untested): coarse spacing at the wall puts the effective wall position and the echo phase off by about `k dr`, which only shows when the echo is not negligible. `graded_grid_for` warns when the a-priori echo exceeds 0.01. Test C7c (criterion fixed before the run): at R = 100 a, al32 rc, 52/60/70 Hz, gain within 0.3 % and phase within 0.3 deg of the uniform single-layer 25x200 values.

### Wall-echo diagnostic on the map

`eta_boundary(t)` (the trace that cuts the fit window at `|eta_boundary| > 0.10 A_ref`) is, on the uniform grid, the maximum of `|eta|` over the outermost `ceil(2 sr)` nodes, that is over the radii `r_wall - k dr0`, `k = 1 .. ceil(2 sr)`. On a graded grid it is the maximum over the same physical radii `x_k = r_wall - k dr0` of the interpolated `eta(x_k) = sum_m w_km eta_{j0+m-1}`, with the Lagrange cubic through the 4 nodes `j0 .. j0+3` (`j0 = q-1` for the panel `q` that holds `x_k`, clamped to `0 .. nr-4`; last panel: extrapolated from the last 4 nodes). Code: `boundary_stencil`, `compute_eta_boundary(...; grid)`. The uniform grid keeps the old code path (bit-identical, tested). Before this fix a graded grid used the nodes within 2 disk radii, at least 3, which sample a different set of radii; H2 below was the suspicion.

## 5. Pass/fail criteria (fixed before any run; never loosened afterwards)

| id | what | pass if |
|---|---|---|
| C1 | DtN accuracy on wave packets `exp(-r^2/2s^2) J0(k0 r)`, graded grid vs exact, compared with a uniform grid of spacing `dr_max` at the same wavelength `lambda = ppw * dr_max` | `E_loc(graded) <= 1.5 E_loc(uniform)` for each `(ppw, q)` cell in the table |
| C2 | Laplacian convergence on `g = exp(-r^2/2w^2)` under refinement `p = 1, 2, 4, 8` (`sr -> p sr`, `q -> q^(1/p)`, `dr_max -> dr_max/p`, `ramp -> p ramp`, so the continuous map is the same) | measured order between `p = 2` and `p = 4`: conservative form `>= 1.9`; 4th-order form `>= 3.5`. The 6th-order order is reported only (the map is `C^4`, so a limit is expected) |
| C3 | Physics vs the existing uniform rows (same physical domain), section 6 | for EVERY selected row: lab gain within 0.3 % (relative) and lab phase within 0.3 deg (wrapped) of the uniform row |
| C4 | 30-case experimental validation, 50 x 100 domain, graded per frequency | amplitude RMSE `<= 1.02 x 0.022135 = 0.022578` and wrapped-phase RMSE `<= 3.0557 + 0.1 = 3.1557 deg` (baseline: uniform single-layer; uniform polar was 0.022186 / 3.0629 deg) |
| C5 | tests | the 362 existing tests pass; new tests pass: uniform grid bit-identical; map and Laplacian order; DtN on a graded grid vs exact; small end-to-end graded vs uniform |
| C6 | cost | reported, not pass/fail: nodes, DtN build time, LU memory, wall time per case, graded vs `nr = 8750` |

`E_loc` (C1): `max_i |(A f)_i - (N f)(r_i)| / max_{|r_j - r_i| <= lambda} |(N f)(r_j)|`, over nodes whose local scale is at least `1e-3 max|N f|`. The exact `N f` is the Hankel integral `int_0^inf k^2 F(k) J0(kr) dk`, `F(k) = s^2 exp(-s^2 (k-k0)^2/2) exp(-x) I0(x)`, `x = s^2 k k0`. The global error of `docs/dtn_single_layer.md` (max error over max exact value) is reported too. The envelope is `s = R/4`, so the packet has weight in the coarse part of the grid. The packet table uses `sr = 25`, `R = 120`, `dr_max = 0.5`, `r_in = 2`, `ramp = 20`, `smooth = 0.1`, `ppw` in 8, 12, 16, 24, `q` in 1.03, 1.05, 1.08.

C3 rows (28, from `bead_gain_phase_bdf2.csv`, `role = production`, `quality = good`; the driver checks that each row exists and is good, and stops if not). There is no good row at 10 Hz in the file, so 15 Hz is the lowest "10 or 15" frequency:

* 25 x 700 (graded domain `R = 350 a`): al32 rc 20, 30, 45, 60, 90 Hz; al32 R 15, 25, 50; al16 rc 25, 40; al16 R 20; nylon rc 15, 35, 52; nylon R 25; ptfe rc 15, 30, 45; ptfe R 20.
* 25 x 200 (graded domain `R = 100 a`): al16 R 40, 70, 90; nylon R 30, 60; ptfe R 45, 80; ptfe rc 60, 90; al16 rc 55; nylon rc 65.

The 25 x 200 rows were computed with the paper's polar matrix. The lead measured 0.03 % / 0.02 deg for the swap polar to single-layer at 25 x 200, and 0.12 % / 0.17 deg for the same-domain swap at 25 x 700, so the 0.3 % / 0.3 deg bound is fair. Settings as the campaign: Gamma 0.5, BDF2, 60 steps per period, fit by `gain_phase.jl` (last 3 periods before the wall echo, same rule).

Configurations tried in C3 (`laplacian_order`, `ppw`, `q`), `r_in = 2`, `ramp = 20`, `smooth = 0.1`:
`(2,12,1.04) (2,16,1.04) (2,24,1.04) (2,32,1.04) (4,8,1.04) (4,12,1.04) (4,16,1.04) (4,12,1.08) (4,12,1.02)`.
A configuration passes C3 when all 28 rows pass. The best configuration is the one that passes with the fewest nodes.

Falsification, fixed now: if the passing configuration needs more than about 2000 nodes at 10 to 90 Hz, or no configuration passes C3, the idea fails for this domain. The report will say so.

## 6. Hypotheses (not yet measured)

* H1. The second-order tridiagonal Laplacian has a relative symbol error of `(k dr)^2/12`: 2.3 % at 12 nodes per wavelength, 1.3 % at 16, 0.6 % at 24. The capillary term uses it, so the wave number is off by about half of that. This is probably too large for C3 at 12 to 16 nodes per wavelength. Then C3 needs the 4th-order (not tridiagonal) form, or `ppw` near 40.
* H2. The wall echo and the wave-cut diagnostic depend on how many nodes sit in the last two disk radii. With a coarse wall region the boundary trace uses at least 3 nodes (see `compute_eta_boundary` with `grid`). This may move the cut time slightly.
* H3. The graded DtN is as accurate as the uniform one at equal local nodes per wavelength (a Python prototype of the same algorithm, run on the login node before the Julia code, gave an error ratio 0.63 to 1.03 for ppw 8 to 16; the Julia job gives the number that counts).

## 7. Reproduction

| What | Script |
|---|---|
| C1, C2, C3, C6 | `julia/cluster/graded_pilot.sh` (calls `julia/scripts/graded_pilot.jl`) |
| C4 and the test tier | `julia/cluster/graded_validate.sh` |
| Hook for one bead case | `julia/cluster/beads/run_bead_sweep.jl ... grid=<ppw>,<q>,<order>` (one frequency per call) |
