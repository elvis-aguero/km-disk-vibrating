% benchmark_solvers.m
% Compares LU-cached direct solver vs GMRES iterative solver.
% Runs 5 forcing periods on a non-resonant case (sub-threshold bath amplitude),
% records CoM trajectories, checks agreement, and reports timing.
%
% Usage: benchmark_solvers  (no arguments needed)

function benchmark_solvers()

fprintf('=== Solver Benchmark: LU-cached vs GMRES ===\n\n');

domains = {struct('spatialResolution', 5,  'bathDiameter', 20), ...
           struct('spatialResolution', 20, 'bathDiameter', 20)};

for di = 1:length(domains)
    sp = domains{di}.spatialResolution;
    bd = domains{di}.bathDiameter;
    fprintf('--- Domain: D%dQuant%d (nr=%d) ---\n', sp, bd, ceil(sp*bd/2));
    run_benchmark(sp, bd);
    fprintf('\n');
end

end

% -------------------------------------------------------------------------
function run_benchmark(spatialResolution, bathDiameter)

% Physics (water bath, sub-threshold: bathAmplitude << Faraday threshold)
diskRadius        = 0.05;    % cm
diskMass          = 0.1;    % g
forceAmplitude    = 0.0;    % dynes (no disk forcing)
forceFrequency_Hz = 90.0;   % Hz
bathAmplitude     = 0.005;  % cm  — very small, well below Faraday threshold
bathFrequency_Hz  = 90.0;   % Hz
phaseDiff_deg     = -90.0;  % degrees
bathDensity       = 1.0;    % g/cm^3
bathST            = 72.20;  % dynes/cm
bathVisc          = 0.978e-2; % Stokes
g_cgs             = 981;    % cm/s^2
temporalRes       = 20;     % steps per adimensional unit
nPeriods          = 5;

% Angular frequencies
forceFrequency = forceFrequency_Hz * 2 * pi;
bath_angular_freq = bathFrequency_Hz * 2 * pi;

% Characteristic units
L_unit = diskRadius;
M_unit = bathDensity * L_unit^3;
T_unit = 1 / forceFrequency;
V_unit = L_unit / T_unit;

% Dimensionless numbers
Re = L_unit^2 / (bathVisc * T_unit);
Fr = L_unit / (g_cgs * T_unit^2);
We = bathDensity * L_unit^3 / (bathST * T_unit^2);

% Adimensional parameters
force_adim        = forceAmplitude / diskMass * (T_unit^2 / L_unit);
surface_force_adim= 2*pi*diskRadius*bathST/diskMass * (T_unit^2 / L_unit);
freq_adim         = forceFrequency * T_unit;
obj_mass_adim     = diskMass / M_unit;
bath_freq_adim    = bath_angular_freq * T_unit;
phase_rad         = phaseDiff_deg * pi / 180;
gamma             = bathAmplitude * bath_angular_freq^2 / g_cgs;

% Initial bath velocity offset (zero lab frame velocity)
v_bath_0_adim = (bathAmplitude * bath_angular_freq * sin(phase_rad)) / V_unit;

% Time step locked to bath period (exact periodicity)
effective_w = bath_freq_adim;
stepsPerCycle = round((2*pi / effective_w) * temporalRes);
dt = (2*pi / effective_w) / stepsPerCycle;
simulationTime = nPeriods * (2*pi / effective_w) * T_unit;   % seconds
steps = ceil(simulationTime / (dt * T_unit));

fprintf('  stepsPerCycle=%d, dt=%.6f (adim), total steps=%d\n', stepsPerCycle, dt, steps);
fprintf('  gamma=%.4f, Re=%.2f, We=%.4f, Fr=%.6f\n', gamma, Re, We, Fr);

% Load domain
nr = ceil(spatialResolution * bathDiameter / 2);
currfold = fileparts(mfilename('fullpath'));
fold = fullfile(currfold, '..', sprintf('D%dQuant%d', spatialResolution, bathDiameter));
if spatialResolution == 5 && bathDiameter == 20
    load(fullfile(fold, 'DTNnew345nr50D5refp10.mat'), 'DTNnew345');
else
    load(fullfile(fold, sprintf('DTNnew345nr%dD%drefp10.mat', nr, bathDiameter)), 'DTNnew345');
end
DTN = DTNnew345; clear DTNnew345;

[dr, Delta, pressureIntegral] = domainMaker(bathDiameter, spatialResolution);
cPoints = spatialResolution + 1;
pIntegral = pressureIntegral(cPoints, :);
SF = surface_force_adim;
Ma = obj_mass_adim;

% Pack constants for matrix builder
PC.Re = Re; PC.We = We; PC.Fr = Fr;
PC.Delta = Delta; PC.DTN = DTN;
PC.pIntegral = pIntegral; PC.Ma = Ma;
PC.nr = nr; PC.cPoints = cPoints; PC.dr = dr; PC.SF = SF;
PC.gamma = gamma; PC.bath_freq = bath_freq_adim;
PC.phase = phase_rad; PC.dt = dt;
PC.force_freq = freq_adim; PC.force_amp = force_adim;
PC.stepsPerCycle = stepsPerCycle;

% Initial conditions
b0   = zeros(2*nr, 1);
CoM0 = 0;
CoMv0= v_bath_0_adim;

% ---- Run LU-cached solver ----
fprintf('  Running LU-cached solver...\n');
L_Lib = cell(stepsPerCycle, 1);
U_Lib = cell(stepsPerCycle, 1);
P_Lib = cell(stepsPerCycle, 1);
CoM_lu   = zeros(steps+1, 1);  CoM_lu(1)   = CoM0;
CoMv_lu  = zeros(steps+1, 1);  CoMv_lu(1)  = CoMv0;
b_lu     = b0;

tic_lu = tic;
for k = 1:steps
    [CoM_lu(k+1), CoMv_lu(k+1), b_lu, L_Lib, U_Lib, P_Lib] = ...
        step_lu(k-1, b_lu, CoM_lu(k), CoMv_lu(k), PC, L_Lib, U_Lib, P_Lib);
end
t_lu = toc(tic_lu);
fprintf('  LU-cached:   %.3f s  (%.2f ms/step)\n', t_lu, 1e3*t_lu/steps);

% ---- Run GMRES solver ----
fprintf('  Running GMRES solver...\n');
CoM_gmres  = zeros(steps+1, 1);  CoM_gmres(1)  = CoM0;
CoMv_gmres = zeros(steps+1, 1);  CoMv_gmres(1) = CoMv0;
b_gmres    = b0;
x0_gmres   = zeros(2*nr+2, 1);   % warm-start vector
gmres_iters = zeros(steps, 1);

tic_gm = tic;
for k = 1:steps
    [CoM_gmres(k+1), CoMv_gmres(k+1), b_gmres, x0_gmres, gmres_iters(k)] = ...
        step_gmres(k-1, b_gmres, CoM_gmres(k), CoMv_gmres(k), PC, x0_gmres);
end
t_gm = toc(tic_gm);
fprintf('  GMRES:       %.3f s  (%.2f ms/step, mean iters=%.1f)\n', ...
        t_gm, 1e3*t_gm/steps, mean(gmres_iters));

% ---- Compare ----
diff = abs(CoM_lu - CoM_gmres);
fprintf('  Max |CoM_LU - CoM_GMRES|  = %.3e (adim units)\n', max(diff));
fprintf('  Final CoM_LU   = %.6e\n', CoM_lu(end));
fprintf('  Final CoM_GMRES= %.6e\n', CoM_gmres(end));
fprintf('  Speed ratio LU/GMRES      = %.2fx\n', t_gm / t_lu);

% Save a small comparison plot
fig = figure('Visible','off');
t_axis = dt * T_unit * (0:steps);
subplot(2,1,1);
plot(t_axis*1e3, CoM_lu*L_unit*1e4, 'b-', 'LineWidth', 1.5); hold on;
plot(t_axis*1e3, CoM_gmres*L_unit*1e4, 'r--', 'LineWidth', 1.2);
legend('LU-cached','GMRES'); ylabel('CoM (\mum)'); grid on;
title(sprintf('D%dQuant%d: CoM trajectory', spatialResolution, bathDiameter));
subplot(2,1,2);
semilogy(t_axis*1e3, diff*L_unit*1e4 + 1e-20, 'k-');
ylabel('|CoM_L_U - CoM_G_M_R_E_S| (\mum)'); xlabel('time (ms)'); grid on;
title('Pointwise difference');
saveas(fig, fullfile(currfold, sprintf('benchmark_D%dQ%d.png', spatialResolution, bathDiameter)));
fprintf('  Plot saved: benchmark_D%dQ%d.png\n', spatialResolution, bathDiameter);

end

% -------------------------------------------------------------------------
function [CoM_next, CoMv_next, b_next, L_Lib, U_Lib, P_Lib] = ...
         step_lu(step_idx, b, CoM, CoMv, PC, L_Lib, U_Lib, P_Lib)

cycleIdx = mod(step_idx, PC.stepsPerCycle) + 1;
g_pf = 1 - PC.gamma * cos(PC.bath_freq * (step_idx+1) * PC.dt + PC.phase);
ft   = PC.dt * PC.force_amp * cos(PC.force_freq * (step_idx+1) * PC.dt);
indep = [b; CoMv - PC.dt/PC.Fr * g_pf - ft; CoM];

if isempty(L_Lib{cycleIdx})
    Mat = build_mat(PC, g_pf);
    [L_Lib{cycleIdx}, U_Lib{cycleIdx}, P_Lib{cycleIdx}] = lu(Mat);
end
sol = U_Lib{cycleIdx} \ (L_Lib{cycleIdx} \ (P_Lib{cycleIdx} * indep));

[CoM_next, CoMv_next, b_next] = unpack(sol, PC);
end

% -------------------------------------------------------------------------
function [CoM_next, CoMv_next, b_next, x0_next, iters] = ...
         step_gmres(step_idx, b, CoM, CoMv, PC, x0)

g_pf = 1 - PC.gamma * cos(PC.bath_freq * (step_idx+1) * PC.dt + PC.phase);
ft   = PC.dt * PC.force_amp * cos(PC.force_freq * (step_idx+1) * PC.dt);
indep = [b; CoMv - PC.dt/PC.Fr * g_pf - ft; CoM];
Mat   = build_mat(PC, g_pf);

[sol, ~, ~, iter] = gmres(Mat, indep, [], 1e-10, size(Mat,1), [], [], x0);
iters = iter(end);
x0_next = sol;
[CoM_next, CoMv_next, b_next] = unpack(sol, PC);
end

% -------------------------------------------------------------------------
function Mat = build_mat(PC, g_pf)
% Thin adapter so the benchmark uses the canonical buildSystemMatrix.
% Translates benchmark's PC field names to those buildSystemMatrix expects.
PC_canon = struct('weber', PC.We, 'reynolds', PC.Re, 'froude', PC.Fr, ...
                  'laplacian', PC.Delta, 'DTN', PC.DTN, ...
                  'pressure_integral', PC.pIntegral, 'obj_mass', PC.Ma);
Mat = buildSystemMatrix(PC_canon, g_pf, PC.dt, PC.nr, PC.cPoints, PC.dr, PC.SF);
end

% -------------------------------------------------------------------------
function [CoM_next, CoMv_next, b_next] = unpack(sol, PC)
nr      = PC.nr;
cPoints = PC.cPoints;
eta_next = [sol(end)*ones(cPoints,1); sol(1:nr-cPoints)];
phi_next = sol(nr-cPoints+1 : 2*nr-cPoints);
CoMv_next = sol(end-1);
CoM_next  = sol(end);
b_next   = [eta_next; phi_next];
end
