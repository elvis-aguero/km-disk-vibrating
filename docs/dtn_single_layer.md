# DtN matrix construction: the paper's scheme and a single-layer alternative

Status: this is the default construction. `generate_dtn(nr, D)` uses `method = :single_layer` (code in `julia/src/dtn/dtn_single_layer.jl`). The paper's scheme remains available as `generate_dtn_polar` or `method = :polar`, and the environment variable `KMDISK_DTN_METHOD=polar` makes the registry resolve the older caches. The exploratory scripts used for the evidence below are in `julia/scripts/pilot_dtn_*.jl`; they call `generate_dtn_polar` as the reference.

Reference: C. A. Galeano-Rios, P. A. Milewski, J.-M. Vanden-Broeck, *Non-wetting impact of a sphere onto a bath and its application to bouncing droplets*, J. Fluid Mech. 826 (2017). Equation numbers below are the paper's.

## Why this matters

Each simulation needs the Dirichlet-to-Neumann (DtN) matrix $N$ of its domain. The paper's construction costs about $n_r^3$ operations times a large constant. At $n_r=2500$ that is 36 core-hours (measured). For the $n_r=8750$ matrix of a large domain it would be about 1,500 core-hours (an estimate scaled from the measurement). The construction below builds the $n_r=8750$ matrix in 222 s on one core (measured). This note shows what it computes, why that is a sound way to get $\varphi_z$, and where it is weaker than the paper's matrix.

## 1. What the matrix must do

Let $\varphi(\mathbf r)$ be the velocity potential on the free surface (the plane $z=0$), with decay at infinity. The fluid occupies $z<0$. Inside the fluid, $\varphi$ is harmonic. The DtN map takes the surface values $\varphi$ and returns the normal derivative $\varphi_z$ there. The paper gives it as a singular integral (Eq. 2.23):

$$
N\varphi(\mathbf r)=\frac{1}{2\pi}\lim_{\epsilon\to0^+}\int_{\mathbb R^2\setminus B(\mathbf r;\epsilon)}\frac{\varphi(\mathbf r)-\varphi(\mathbf s)}{\lvert\mathbf r-\mathbf s\rvert^{3}}\,dA(\mathbf s).
$$

$B(\mathbf r;\epsilon)$ is the disk of radius $\epsilon$ around $\mathbf r$. The kernel $\lvert\mathbf r-\mathbf s\rvert^{-3}$ is too singular to integrate on its own. The difference $\varphi(\mathbf r)-\varphi(\mathbf s)$ cancels the singularity, so the limit exists.

Take the waves $e^{i\mathbf k\cdot\mathbf x}$. Any operator that treats all positions alike turns such a wave into the same wave times a number that depends only on $\mathbf k$. That number is the *symbol* of the operator. For example, $d/dx$ has symbol $ik$, and $-\Delta_H$ (minus the horizontal Laplacian) has symbol $\lvert\mathbf k\rvert^2$. The symbol of $N$ follows from the physics. The harmonic function that equals $e^{i\mathbf k\cdot\mathbf x}$ at $z=0$ and decays downward is $e^{i\mathbf k\cdot\mathbf x}e^{\lvert\mathbf k\rvert z}$. Its $z$-derivative at $z=0$ is $\lvert\mathbf k\rvert$ times itself. So $N$ multiplies each wave by its wavenumber $\lvert\mathbf k\rvert$. Operators applied one after another have multiplied symbols.

The solver works with an axisymmetric surface on a radial grid of $n_r$ nodes $r_i=i\,\delta r$, $i=0,\dots,n_r$, with $\varphi=0$ at the last node. Its DtN matrix is the $n_r\times n_r$ matrix that maps the nodal values $\varphi_i$ to approximate values of $N\varphi$ at the nodes.

## 2. How the paper builds the matrix (Appendix B)

The paper approximates the integral above directly, one node at a time. For node $r_i$ it builds a polar mesh centred there: points $\mathbf s_{lm}=\mathbf r_i+\rho_l(\cos\theta_m,\sin\theta_m)$ with radial step $\Delta r=\delta r/10$ and angular step $\Delta\theta$, small enough that $\Delta\theta/2<\Delta r/10$. Each point is the centre of a cell. The function $\varphi$ at $\mathbf s_{lm}$ is interpolated cubically from the radial nodes: $\varphi(\mathbf s_{lm})\approx\sum_j w_j(\lvert\mathbf s_{lm}\rvert)\varphi_j$. The paper says only "cubic interpolation" and does not specify which cubic. The production code uses piecewise four-point Lagrange cubics: on the interval between two neighbouring nodes, $w_j$ is the cubic that equals 1 at node $j$ and 0 at the other three of the four nearest nodes. The kernel is integrated exactly over each cell:

$$
\int_{\rho_l-\Delta r/2}^{\rho_l+\Delta r/2}\int_{\theta_m-\Delta\theta/2}^{\theta_m+\Delta\theta/2}\rho^{-3}\,\rho\,d\theta\,d\rho=\Delta\theta\,K_l,\qquad K_l=\frac{1}{\rho_l-\Delta r/2}-\frac{1}{\rho_l+\Delta r/2}.
$$

The ball of radius $\Delta r/2$ around $\mathbf r_i$ is skipped, which approximates the limit $\epsilon\to0$. Collecting the terms gives the matrix entries:

$$
N_{ij}=\frac{\delta_{ij}}{2\pi}\sum_{l,m}\Delta\theta\,K_l\;-\;\frac{1}{2\pi}\sum_{l,m}\Delta\theta\,K_l\;w_j\!\left(\lvert\mathbf s_{lm}\rvert\right).
$$

The code also treats the few cells nearest the node with analytic formulas. For each of the $n_r$ rows the sum runs over about $20\,n_r$ radial steps times $10\pi\,n_r$ angles. So the cost is $O(n_r^3)$ with a large constant.

## 3. A different route to the same operator

The expensive part is sampling the angle around each node. The route below does the angular integral exactly instead. It rests on a second operator.

**The single-layer operator $S$.** For a function $\psi$ on the plane, define

$$
(S\psi)(\mathbf r)=\frac{1}{2\pi}\int_{\mathbb R^2}\frac{\psi(\mathbf s)}{\lvert\mathbf r-\mathbf s\rvert}\,dA(\mathbf s).
$$

$S$ is a convolution with $1/(2\pi\lvert\mathbf x\rvert)$, so its symbol is the Fourier transform of that kernel:

$$
\int_{\mathbb R^2}\frac{e^{-i\mathbf k\cdot\mathbf x}}{\lvert\mathbf x\rvert}\,dA=\int_0^\infty\!\!\int_0^{2\pi}e^{-ikr\cos\theta}\,d\theta\,dr=2\pi\int_0^\infty J_0(kr)\,dr=\frac{2\pi}{\lvert\mathbf k\rvert},\qquad\text{so}\qquad \hat S=\frac1{\lvert\mathbf k\rvert}.
$$

The last integral converges conditionally. With a damping factor, $\int_0^\infty e^{-\epsilon r}J_0(kr)\,dr=1/\sqrt{k^2+\epsilon^2}\to1/\lvert\mathbf k\rvert$. This needs no harmonic functions. Harmonicity enters when $S$ is linked to $N$, through three exact statements. None of them swaps a limit with an integral or a derivative.

1. **$SN=I$.** This is the paper's Eq. (A 3). For $\varphi$ harmonic in $z<0$ with the decay (A 1), Green's second identity gives $\varphi\rvert_{z=0}=S(\varphi_z\rvert_{z=0})$. Every decaying $f$ on the plane is the surface value of exactly one such $\varphi$, and $Nf$ is by definition $\varphi_z\rvert_{z=0}$. So $S(Nf)=f$.
2. **$NS=I$.** Given a decaying $\psi$, define for $z<0$
   $$
   w(\mathbf x,z)=\frac{1}{2\pi}\int_{\mathbb R^2}\frac{\psi(\mathbf s)}{\sqrt{\lvert\mathbf x-\mathbf s\rvert^{2}+z^{2}}}\,dA(\mathbf s).
   $$
   It is harmonic and decays. As $z\to0^-$, $w_z\to\psi$, because $\int_{\mathbb R^2}\frac{-z}{2\pi(\rho^2+z^2)^{3/2}}\,dA=1$ for $z<0$: the kernel of $w_z$ tends to a delta function. So $w$ has $z$-derivative $\psi$ and surface value $S\psi$, and by the definition of $N$, $N(S\psi)=\psi$. With statement 1, $N=S^{-1}$.
3. **$N=-\Delta_H S$.** Take $w$ from statement 2 with $\psi=f$. Then $w\rvert_{z=0}=Sf$ and $w_z\rvert_{z=0}=f$. The function $w_z$ is harmonic with surface value $f$, so by the definition of $N$,
   $$
   Nf=(w_z)_z\big\rvert_{z=0}=w_{zz}\big\rvert_{z=0}=-\Delta_H w\big\rvert_{z=0}=-\Delta_H(Sf).
   $$
   This uses only Laplace's equation, $w_{zz}=-\Delta_H w$, and that $w$ is $C^2$ up to $z=0$, which holds for smooth decaying $f$.

In symbols: $\lvert\mathbf k\rvert^{2}\cdot\lvert\mathbf k\rvert^{-1}=\lvert\mathbf k\rvert$. $S$ is a left inverse of $N$. It does not map the image of $N$ into its own kernel, because $S$ has trivial kernel (its symbol never vanishes).

**The angular integral in closed form.** For an axisymmetric $\psi(r')$, the integral over the ring $\lvert\mathbf s\rvert=r'$ is

$$
\int_0^{2\pi}\frac{d\theta}{\sqrt{r^{2}+r'^{2}-2rr'\cos\theta}}=\frac{4K(m)}{r+r'},\qquad m=\frac{4rr'}{(r+r')^{2}},
$$

where $K$ is the complete elliptic integral of the first kind:

$$
K(m)=\int_0^{\pi/2}\frac{d\theta}{\sqrt{1-m\sin^{2}\theta}}.
$$

$K$ has no formula in elementary functions. Treat it like $\sin$ or a Bessel function: a library function with known values, here computed with the arithmetic-geometric mean in about 6 iterations. Two checks. On the axis ($r=0$) every point of the ring is at distance $r'$, so the ring integral is $2\pi/r'$, and the formula gives $4\cdot(\pi/2)/r'$ since $K(0)=\pi/2$. When $r'\to r$, $K$ grows like $\log$ and the integrand has a weak logarithmic singularity. Therefore

$$
(S\psi)(r)=\frac{2}{\pi}\int_0^\infty\frac{r'}{r+r'}\,K(m)\,\psi(r')\,dr'.
$$

The sampled sum over about $10^5$ angles in the paper becomes one function evaluation, and the singularity in the kernel is now only logarithmic.

## 4. Turning it into a matrix

**Nodal values to the surface.** The pilot interpolates $\varphi$ with piecewise four-point Lagrange cubics. On the interval between nodes $r_q$ and $r_{q+1}$, $\psi(r')\approx\sum_j w_j(r')\,\varphi_j$, where $w_j$ is the cubic that equals 1 at node $j$ and 0 at the other three of the four nodes $r_{q-1},\dots,r_{q+2}$. Away from the axis these are the weights the production code uses. Near the axis the pilot extends $\varphi$ evenly about $r=0$ ($\varphi_{-1}=\varphi_1$). That does not give the same weights as the production code's first interval, which has its own special weights. Beyond the last node $\varphi=0$. The paper says only "cubic interpolation" and does not specify which cubic, so it is not known whether either choice is the paper's. Substituting the interpolant into the formula for $S$ gives, at node $i$,

$$
(S\psi)(r_i)\approx\sum_j G_{ij}\,\varphi_j,\qquad G_{ij}=\frac{2}{\pi}\int_0^{R}\frac{r'}{r_i+r'}\,K\big(m_i(r')\big)\,w_j(r')\,dr',\quad R=n_r\,\delta r .
$$

So $G$ is the matrix of the discretized $S$: it maps nodal values of $\varphi$ to the values $g_i$ of $S\varphi$ at the nodes. It is not $S$ itself.

**Computing the integrals.** A *panel* is the interval between two neighbouring nodes. Each panel is integrated separately with an $n$-point Gauss-Legendre rule. The rule replaces the integral by a weighted sum of the integrand at $n$ fixed points, and it is exact for polynomials up to degree $2n-1$. It works well when the integrand is smooth on the panel. The pilot uses $n=24$.

**Where the rule fails.** On every panel except the two that touch $r'=r_i$, the kernel $K(m)$ is smooth. On those two panels the integrand behaves like $\log\lvert r'-r_i\rvert$ at the end $r'=r_i$. A logarithm is not close to any polynomial there, so the rule converges slowly. Example: $\int_0^1\log r\,dr=-1$. With 6, 12 and 24 points the error is $1.5\times10^{-2}$, $4.0\times10^{-3}$ and $1.1\times10^{-3}$. The finite differences in the next step multiply the error of $G$ by $1/\delta r^{2}$, so about $10^{-10}$ is needed, not $10^{-3}$.

**The substitution.** On a panel of width $h$ with the singular end at $r'=r_i$, change the integration variable to $u$:

$$
r'=r_i+h\,u^{4}\ \ (\text{panel to the right of node }i),\qquad r'=r_i-h\,u^{4}\ \ (\text{panel to the left}),\qquad 0\le u\le1 .
$$

Then $dr'=4h\,u^{3}\,du$ and $\log\lvert r'-r_i\rvert=\log h+4\log u$. The integrand in $u$ contains $u^{3}\log u$. That function has two continuous derivatives at $u=0$; only its third derivative is singular. So it is much closer to a polynomial than $\log r$ is. The same Gauss-Legendre rule is applied in $u$. On the example above the error becomes $5.4\times10^{-6}$, $2.5\times10^{-8}$ and $1.1\times10^{-10}$ for 6, 12 and 24 points.

The word *graded* describes where the sample points land in $r'$. The points are evenly arranged in $u$, so in $r'$ they crowd toward the singular end and spread out away from it. With 24 points, the nearest sample is at $3\times10^{-11}\,h$ from $r_i$ (against $2.4\times10^{-3}\,h$ for the plain rule), and the farthest is at $0.99\,h$. The mesh is graded toward the singularity, like a mesh that is refined near a point where the solution is rough.

**The Laplacian.** By statement 3, $N\approx-\mathcal L\,G$, where $\mathcal L$ is a finite-difference matrix for $\Delta_H g=g''+g'/r$ (the paper's B 2) on the nodal values $g_i$. The paper uses second-order stencils for this operator (B 3, B 4). The pilot uses sixth order:

$$
g''(r_i)\approx\frac{2g_{i-3}-27g_{i-2}+270g_{i-1}-490g_i+270g_{i+1}-27g_{i+2}+2g_{i+3}}{180\,\delta r^{2}},
$$

$$
g'(r_i)\approx\frac{-g_{i-3}+9g_{i-2}-45g_{i-1}+45g_{i+1}-9g_{i+2}+g_{i+3}}{60\,\delta r},\qquad \Delta_H g(0)=2g''(0)\ \text{(even symmetry)}.
$$

The finite differences divide by $\delta r^2$, so they multiply any error in $G$ by that factor. $G$ therefore has to be accurate to about $10^{-10}$, which the quadrature achieves. The cost is $O(n_r^{2})$ kernel evaluations.

## 5. Why this is a sound way to compute $\varphi_z$

1. **The identity is exact.** $N=-\Delta_H S$ follows from Laplace's equation and the definition of $N$. The construction adds no model assumption.
2. **The hard part is done to quadrature accuracy.** The paper's scheme removes the singularity by sampling a difference $\varphi(\mathbf r)-\varphi(\mathbf s)$ on a polar mesh, so the angular integral is a sampled sum. Here the angular integral is exact and the remaining singularity is logarithmic, which the graded substitution integrates to near machine precision. What is left are two truncation errors: the cubic interpolation of $\varphi$ (the same kind as the production code's, away from the axis) and the sixth-order finite differences.
3. **The error is measured, as a function of wavelength.** The table applies the paper's matrix (the production algorithm, `generate_dtn`, $\Delta r=\delta r/10$) and the new matrix to wave packets $e^{-r^2/2s^2}J_0(k_0r)$, the shape of the paper's own test (Appendix B.1), and compares with the exact $N\varphi$ from the Hankel transform. Error is the maximum over nodes of $\lvert Af-\text{exact}\rvert$ divided by the maximum of the exact value, at $n_r=200$ with envelope width $s=25$ nodes. Job 8462521, `julia/scripts/pilot_dtn_spectrum.jl`.

   | Wavelength (nodes) | Paper's matrix | New, 4th-order stencil | New, 6th-order stencil |
   |---:|---:|---:|---:|
   | 2.5 | 74% | 28% | 22% |
   | 3 | 47% | 18% | 13% |
   | 4 | 17% | 7.5% | 5.2% |
   | 6 | 3.1% | 1.9% | 1.3% |
   | 8 | 0.81% | 0.67% | 0.45% |
   | 12 | 0.13% | 0.16% | 0.11% |
   | 16 | 0.035% | 0.058% | 0.042% |
   | 24 | 0.0073% | 0.015% | 0.012% |
   | 32 | 0.0029% | 0.0063% | 0.0050% |
   | 64 | 0.00055% | 0.0012% | 0.00098% |
   | 150 | 0.00018% | 0.00034% | 0.00030% |

   At wavelengths of 12 nodes and longer, the error of the new matrix is 0.9 to 1.9 times the error of the paper's matrix. At 8 nodes and shorter, the error of the new matrix is smaller. Below about 4 nodes both are poor: the paper's matrix is off by 17 to 74% there, and the new one by 5 to 28%.
4. **The relevant waves are long.** At 50 Hz the capillary wavelength is about 800 nodes for the Al 1/32" bead and about 170 nodes for PTFE. The bead and the disk edge involve shorter scales (the disk radius is 25 nodes), which is why the next point matters.
5. **The solver output agrees with production.** This is the test that includes the disk edge. Each row runs the real solver with the new matrix at the production domain 25 x 200 and compares with the production campaign (same commit, same settings). Jobs 8458894 and 8459007.

   | Case ($\Gamma=0.5$, BDF2, 60 steps/period) | Lab gain, relative difference | Lab phase difference |
   |---|---:|---:|
   | PTFE, R, 40 Hz | 0.009% | 0.009 deg |
   | PTFE, R, 60 Hz | 0.000% | 0.016 deg |
   | PTFE, R, 90 Hz | 0.014% | 0.018 deg |
   | Al 1/32", rc, 40 Hz | 0.011% | 0.001 deg |
   | Al 1/32", rc, 60 Hz | 0.026% | 0.005 deg |
   | Al 1/32", rc, 70 Hz | 0.030% | 0.010 deg |
6. **The large matrix has an independent check.** Al 1/32" at rc in the 25 x 700 domain ($n_r=8750$, built in 222 s, job 8459008) gives lab gains 1.2266, 1.4586, 1.5622 at 40, 60, 70 Hz. The production traces, which use the old matrix, refit on periods 2 to 3 (before the wall echo arrives) give 1.2258, 1.4576, 1.5616. The two routes share no code for the matrix and agree within 0.07% in gain.
7. **The experimental measurements show no regression.** The 30-case validation sweep (disk $R=0.2$ cm, domain 50 x 100, $n_r=2500$, $\Gamma=0.05, 0.2, 0.5$, 10 to 100 Hz, BDF2) against `matlab/0_data/external/experimental_measurements`, with the old matrix (the registered cache) and the new one built in the job (24.5 s). The criterion was fixed before the runs: new amplitude RMSE at most 1.02 times the old one, and new wrapped-phase RMSE at most the old one plus 0.1 deg. Jobs 8466188 and 8466189, same commit and settings.

   | | Old matrix | New matrix | Criterion |
   |---|---:|---:|---:|
   | Amplitude RMSE | 0.022186 | 0.022135 | $\le$ 0.022630 |
   | Phase RMSE (wrapped) | 3.0629 deg | 3.0557 deg | $\le$ 3.163 deg |

   The centre-of-mass RMS of the 30 cases differs by at most 0.032% between the two runs.
8. **Cost.** Measured: 0.04 to 0.17 s at $n_r=200$ to 400, against 66 s and 530 s for the paper's scheme; 24 s at $n_r=2500$; 222 s at $n_r=8750$.

## 6. Where the two matrices differ, and what is not checked

The two matrices differ by 22% in Frobenius norm at $n_r=2500$ (sixth-order stencil; 29% with fourth order). That number is real. The table in section 5 shows where it comes from. Most modes of an $n_r$-dimensional space are at the shortest wavelengths, and the Frobenius norm weights all modes equally. At wavelengths of 4 nodes and shorter, the paper's matrix has errors of 17 to 74% against the exact operator, and the new one has 5 to 28%. The size of the gap is therefore mostly the paper's matrix's own error at the shortest scales. It is not evidence that the new matrix is wrong there. Below 4 nodes neither matrix represents the operator, and the gap between them says little about the physics.

The scale that does matter is the contact line. The disk radius is 25 nodes, but the pressure under the disk has structure at a few nodes near its edge. That is where both matrices are weakest, and it is why the solver comparison in section 5 is the test that counts.

Not checked, or weaker:

- **Symmetry is not a usable quality test here.** The exact $N$ is self-adjoint, but a discrete matrix is symmetric only in the inner product that matches its discretization. With rows weighted by ring area $W_i$, the paper's matrix has relative asymmetry $4\times10^{-4}$ and the new one $3.6\times10^{-2}$. Each matrix was replaced by its weighted symmetric part, $\tfrac12(A+W^{-1}A^{T}W)$, and the wavelength test repeated (job 8464406). Both matrices got much worse. At 12 nodes the error went from 0.13% to 240% for the paper's matrix and from 0.11% to 19% for the new one. At 150 nodes it went to 1970% and 170%. So the ring-area weights are not the inner product in which either matrix is self-adjoint. A small asymmetry in the matrix norm is made of entries of order one, while $N\varphi$ for a long wave is a small difference of such entries. A relative change of $10^{-4}$ in the entries therefore changes the long-wave result by order one. The asymmetry number does not predict accuracy, and symmetrizing with these weights is harmful.
- **Kinks.** Piecewise-cubic interpolation is continuous but has derivative jumps at every node. For such a function, $w$ is not $C^{2}$ up to the boundary, so statement 3 holds only in the distribution sense, and the exact $N$ of the interpolant is logarithmically singular at each node. The paper's exclusion ball and the finite differences after the integral each regularize that differently at the grid scale. This is the likely reason the matrices differ there. It has not been tested.
- **Function spaces.** $S\psi$ has a tail $\sim(\int\psi)/(2\pi\lvert\mathbf x\rvert)$, which breaks the paper's decay condition (A 1) unless $\int\psi=0$. That is consistent: a function that decays faster than $1/\lvert\mathbf x\rvert$ has $\int N f=0$.
- **Coverage.** The solver comparison covers two beads, six cases, at $\Gamma=0.5$. The wall convention is $\varphi=0$ beyond the last node in both matrices, but the paper's treatment of the last cells differs in detail. The first interval near the axis uses different weights from the production code (section 4). The experimental comparison (section 5, point 7) uses one disk and one domain (50 x 100).

## Reproduction

| What | Script | Job |
|---|---|---|
| Build the matrix (cubic, order 4 or 6), write a cache file | `julia/scripts/build_dtn_singlelayer.jl <sr> <D> <ord> <out.jld2>` | 8458894, 8459007, 8459008 |
| Error vs exact Gaussian, linear/cubic, Richardson | `julia/scripts/pilot_dtn_singlelayer.jl` via `julia/cluster/beads/00c_dtn_pilot.sh` | 8458256 |
| Comparison with the production cache | `julia/scripts/pilot_dtn_vs_production.jl` via `julia/cluster/beads/00e_dtn_vs_production.sh` | 8458810 |
| Error vs wavelength, with and without symmetrization | `julia/scripts/pilot_dtn_spectrum.jl` via `00i_dtn_spectrum.sh`, `00j_dtn_spectrum_sym.sh` | 8462521, 8464406 |
| Solver check at 25 x 200 (PTFE R; Al 1/32" rc) | `00f_dtn_e2e.sh`, `00g_dtn_control.sh` | 8458894, 8459007 |
| Experimental validation, old vs new matrix | `julia/cluster/validate_vs_experiment.sh` (`DTN_BUILD=single_layer` for the new one) | 8466188 (old), 8466189 (new) |
| Fast test tier with the new default (362 tests) | `julia/cluster/run_tests.sh` | 8466803 |
| Al 1/32" rc in the 25 x 700 domain | `00h_dtn_big_al32rc.sh` | 8459008 |

The sweep script `julia/cluster/beads/run_bead_sweep.jl` accepts `KMDISK_DTN_FILE=<file.jld2 with key "DTN">` to replace the registered cache for a test run.
