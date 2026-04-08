function [t_s, CoM_cm, eta_history_cm] = solve_motion(NameValueArgs)
% SOLVE_MOTION Simulates the motion of a disk oscillating in a fluid bath.
%   This function models the motion of a disk and a bath, both subjected to 
%   sinusoidal forces using configurable parameters provided through
%   name-value arguments.

% Parse input arguments using name-value pairs and set default values
arguments
    % Unless stated otherwise, all units are cgs. 
    % Default values amount to a water bath
    NameValueArgs.diskRadius (1, 1) double = .2 % Radius of the oscillating disk, in cm
    NameValueArgs.diskMass (1, 1) double = .0283 % Mass in grams of the disk
    NameValueArgs.forceAmplitude (1, 1) double = 0 % Amplitude of sinusoidal force applied to disk (in dynes)
    NameValueArgs.forceFrequency (1, 1) double = 90 % Frequency of sinusoidal force in Hz
    NameValueArgs.bathAmplitude (1, 1) double = 0.01 % Amplitude of bath oscillation in cm. 
    NameValueArgs.bathFrequency (1, 1) double = 90 % Frequency of bath oscillation in Hz. 
    NameValueArgs.phaseDifference (1, 1) double = -90 % Phase difference between disk forcing and bath oscillation in degrees.
    NameValueArgs.bathDensity (1, 1) double = 1 % Density of bath's fluid in g/cm^3
    NameValueArgs.bathSurfaceTension (1, 1) double = 72.20 % For water, in dynes/cm % density 1.175 g/cm3, eta=0.018 Pa s
    NameValueArgs.bathViscosity (1, 1) double = 0.978e-2 % Viscosity in Stokes (cgs)
    NameValueArgs.g (1, 1) double = 981 % Gravitational constant in cgs
    NameValueArgs.bathDiameter (1, 1) double = 100 % Diameter of the bath wrt to disk Radius
    NameValueArgs.spatialResolution (1, 1) double = 50 % Number of numerical radial intervals in one disk radius
    NameValueArgs.temporalResolution (1, 1) double = 30; % Number of temporal steps in one adimensional unit
    NameValueArgs.simulationTime (1, 1) double = 10/90; % Time to be simulated in seconds
    NameValueArgs.debug_flag (1, 1) logical = true; % To show some debugging info
    NameValueArgs.solverType (1, 1) string = "auto" % "auto": LU-cached if RAM allows, else GMRES. Override with "lu" or "gmres".
    NameValueArgs.gmresTolerance (1, 1) double = 1e-10 % GMRES convergence tolerance (used when GMRES is active)
end

% Record the time the script is run
tic;
datetimeMarker = datetime('now'); datetimeMarker.Format = 'yyyyMMddHHmmss'; 
NameValueArgs.datetimeMarker = datetimeMarker;

% Add arguments to the current scope
NameValueArgs.forceFrequency = NameValueArgs.forceFrequency * 2 * pi; % We use angular frequency
cellfun(@(f) assignin('caller', f, NameValueArgs.(f)), fieldnames(NameValueArgs));

% Reset any existing warnings
lastwarn('', '');

% Prepare for file saving
close all
currfold = fullfile(fileparts(mfilename('fullpath'))); 
cd(currfold);

% Add the current folder to the path
addpath(currfold);
precomputedInverse = nan; 
cd ..
% Start logging if debugging is enabled
if debug_flag == true
    logDir = fullfile(fileparts(mfilename('fullpath')), '..', '..', '0_data', 'logs');
    if ~exist(logDir, 'dir'); mkdir(logDir); end
    diary(fullfile(logDir, sprintf("solve_motion%s.txt", datetimeMarker))); 
end

% Set up folder for data saving
fold = fullfile(pwd, sprintf("D%dQuant%d", spatialResolution, bathDiameter));
try
    % Navigate to the folder and load necessary matrices
    cd(fold)
    nr = ceil(spatialResolution * bathDiameter / 2);
    % Machine-specific patch for D5Quant20.
    if spatialResolution == 5 && bathDiameter == 20
        load('DTNnew345nr50D5refp10.mat', 'DTNnew345')
    else
        load(sprintf('DTNnew345nr%dD%drefp10.mat', nr, bathDiameter), 'DTNnew345')
    end
    DTN = DTNnew345;
    clear DTNnew345

    % Load precomputed inverse matrix if it exists
    myfile = fullfile(pwd, sprintf("dtstep=%d.mat", forceFrequency*temporalResolution)); 
    if exist(myfile, "file"); load(myfile, "precomputedInverse"); end
    cd(currfold)
catch ME
    error("Could not load DTN for D=%d, Quant=%d. Please generate the matrix first", ...
        spatialResolution, bathDiameter);
end


% Load useful matrices for the simulation (Laplacian and pressure integral)
[dr, laplacian, pressureIntegral] = domainMaker(bathDiameter, spatialResolution);

% Define characteristic units
L_unit = diskRadius; 
M_unit = bathDensity * L_unit^3; % Mass unit
T_unit = 1 / forceFrequency; % Time unit (1 over angular frequency)
V_unit = L_unit / T_unit; % Velocity unit
F_unit = M_unit * L_unit / T_unit^2; % Force unit

% Store characteristic units in a struct
UNITS = struct('length', L_unit, 'mass', M_unit, 'time', T_unit, ...
    'velocity', V_unit, 'force', F_unit);

% Compute dimensionless numbers
Re = L_unit^2 / (bathViscosity * T_unit); % Reynolds number
Fr = L_unit / (g * T_unit^2); % Froude number
We = bathDensity * L_unit^3 / (bathSurfaceTension * T_unit^2); % Weber number

% Compute adimensionalized force, frequency, and object mass
force_adim = forceAmplitude / diskMass * (T_unit^2 / L_unit);
surface_force_adim = 2*pi*diskRadius*bathSurfaceTension/diskMass * (T_unit^2 / L_unit);
freq_adim = forceFrequency * T_unit;
obj_mass_adim = diskMass / M_unit;

% New bath-specific parameters. 
bath_angular_freq = NameValueArgs.bathFrequency * 2 * pi;
bath_freq_adim = bath_angular_freq * T_unit;
phase_diff_rad = NameValueArgs.phaseDifference * pi / 180;
bath_forcing_amplitude = NameValueArgs.bathAmplitude * bath_angular_freq^2 / g; % Dimensionless acceleration ratio (A*w^2/g)

% Calculate initial velocity offset for zero lab velocity. 
v_bath_0_adim = (NameValueArgs.bathAmplitude * bath_angular_freq * sin(phase_diff_rad)) / V_unit;

% Calculate informed static ylim.
A_bath_adim = abs(NameValueArgs.bathAmplitude / L_unit);
A_forcing_adim = force_adim / (freq_adim^2 + 1e-6); % Rough estimate
H_limit_adim = 1.5 * (A_bath_adim + A_forcing_adim) + 0.02; % Tighter symmetric limit around z=0

% Steps per cycle must be an integer for periodic physics to be valid.
if bath_forcing_amplitude == 0
    effective_w_adim = freq_adim;
elseif forceAmplitude == 0
    effective_w_adim = bath_freq_adim;
else
    effective_w_adim = bath_freq_adim; % Use bath freq as reference
end
stepsPerCycle = temporalResolution; 
dt = (2 * pi / effective_w_adim) / stepsPerCycle; % Adjusted adimensional time step

steps = ceil(simulationTime / (dt * T_unit)); % Minimum number of time steps

% --- Resolve solver: LU-cached (default) or GMRES fallback ---
% stepsPerCycle = round(2*pi * temporalResolution)
% Each factorization stores L, U, P as full dense matrices: 3*(2*nr+2)^2*8 bytes.
luSizeBytes = stepsPerCycle * 3 * (2*nr+2)^2 * 8; %Estimated bytes needed in memory
availRAM    = getAvailableRAM(); %Max RAM available to the computer
hasPCT      = license('test', 'Distrib_Computing_Toolbox');

if solverType == "auto"
    if availRAM >= luSizeBytes
        resolvedSolver = "lu";
    else
        resolvedSolver = "gmres";
        fprintf('RAM check: LU cache needs ~%.1f GB (%d matrices x %.0f MB), available ~%.1f GB. Falling back to GMRES.\n', ...
            luSizeBytes/1e9, stepsPerCycle, luSizeBytes/stepsPerCycle/1e6, availRAM/1e9);
    end
else
    resolvedSolver = solverType; % explicit override by caller
end

%Initial conditions for the fluid
etaInitial = zeros(nr,1); %initial surface elevation
phiInitial = zeros(nr,1); %initial surface potential

% Define the current state of the system
current_conditions = struct( ...
    "dt", dt(1), "time", 0, ...
    "center_of_mass", 0, "center_of_mass_velocity", v_bath_0_adim, ... 
    "bath_surface", etaInitial, "bath_potential", phiInitial, "pressure", zeros(spatialResolution+1, 1));
current_index = 1; %iteration counter
recordedConditions = cell(steps, 1);
recordedConditions{current_index} = current_conditions;

% Store problem constants for use in the simulation
PROBLEM_CONSTANTS = struct("froude", Fr, "weber", We, ...
    "reynolds", Re, "dr", dr, "DEBUG_FLAG", debug_flag, ...
    "nr", nr, "contact_points", spatialResolution+1, ... 
    "force_amplitude", force_adim, "force_frequency", freq_adim, ...
    "bath_forcing_amplitude", bath_forcing_amplitude, ...
    "bath_frequency", bath_freq_adim, "phase_difference", phase_diff_rad, ... 
    "surface_force_constant", surface_force_adim, ...
    "DTN", DTN, "laplacian", laplacian, "obj_mass", obj_mass_adim, ...
    "pressure_integral", pressureIntegral(spatialResolution+1, :), ...
    "precomputedInverse", precomputedInverse, ...
    "ylim_limit", H_limit_adim, "step_counter", 0, ...
    "stepsPerCycle", stepsPerCycle, ...
    "solverType", resolvedSolver, ...
    "gmresTolerance", NameValueArgs.gmresTolerance, ...
    "gmres_x0", zeros(2*nr+2, 1), ...
    "luFutures", {{}}, ...
    "L_Library", {cell(stepsPerCycle, 1)}, ...
    "U_Library", {cell(stepsPerCycle, 1)}, ...
    "P_Library", {cell(stepsPerCycle, 1)});

fprintf("Starting simulation on %s\n", pwd);

% --- Async LU pre-computation ---
% Fire parfeval futures for all stepsPerCycle LU factorizations before the
% main loop starts.  Each advance_one_step call fetches only the future it
% needs (blocks just for that one if not yet ready).
% Without PCT, falls back to lazy-sync LU: LU computed on first miss inside
% advance_one_step, cached for all subsequent cycles.
if resolvedSolver == "lu"
    if hasPCT
        fprintf('Firing %d async LU futures (~%.1f GB)...\n', stepsPerCycle, luSizeBytes/1e9);
        PC_min = struct('weber', We, 'reynolds', Re, 'froude', Fr, ...
                        'laplacian', laplacian, 'DTN', DTN, ...
                        'pressure_integral', pressureIntegral(spatialResolution+1, :), ...
                        'obj_mass', obj_mass_adim);
        luFutures = cell(stepsPerCycle, 1);
        for ci = 1:stepsPerCycle
            g_pf = 1 - bath_forcing_amplitude * cos(bath_freq_adim * ci * dt + phase_diff_rad);
            luFutures{ci} = parfeval(@computeLU, 3, PC_min, g_pf, dt, nr, spatialResolution+1, dr, surface_force_adim);
        end
        PROBLEM_CONSTANTS.luFutures = luFutures;
    else
        fprintf('LU caching active (lazy-sync: no Parallel Computing Toolbox detected).\n');
    end
else
    fprintf('Solver: GMRES (tol=%.0e, warm-started).\n', NameValueArgs.gmresTolerance);
end

% Names of the variables to be saved
savingvarNames = { ...
    getVarName(NameValueArgs), ...
    getVarName(PROBLEM_CONSTANTS), ...
    getVarName(recordedConditions), ...
    getVarName(UNITS) ...
};
variableValues = cell(size(savingvarNames));

%% Main Simulation Loop
try
    while (recordedConditions{current_index}.time * T_unit < simulationTime) 

        [recordedConditions{current_index+1}, PROBLEM_CONSTANTS] = ...
               advance_one_step(recordedConditions{current_index}, ...
                       PROBLEM_CONSTANTS);
        current_index = current_index + 1;

        if PROBLEM_CONSTANTS.DEBUG_FLAG == true
            width = 6; 
            width = min(nr, ceil(width*spatialResolution));
            
            gamma = PROBLEM_CONSTANTS.bath_forcing_amplitude;
            Fr = PROBLEM_CONSTANTS.froude;
            bath_w = PROBLEM_CONSTANTS.bath_frequency;
            phase = PROBLEM_CONSTANTS.phase_difference;
            zb = (gamma / (Fr * bath_w^2)) * cos(bath_w * recordedConditions{current_index}.time + phase);
            
            eta = recordedConditions{current_index}.bath_surface;
            eta = [flipud(eta(2:width)); eta(1:width)];
            xplot = dr*(0:nr-1) * L_unit;
            plot([-fliplr(xplot(2:width)), xplot(1:width)], (eta + zb) * L_unit * 1e4, 'b', 'Linewidth', 2);
            hold on
            
            x = [1, 1, -1, -1] * L_unit;
            z = recordedConditions{current_index}.center_of_mass;
            z_lab = z + zb;
            block_height = 0.15 * PROBLEM_CONSTANTS.ylim_limit; % Dimensionless height relative to ylim
            y = [z_lab, z_lab+block_height, z_lab+block_height, z_lab] * L_unit * 1e4;
            fill(x, y, 'k');

            % Plot bath amplitude markers just outside the axes
            A_bath_um = (gamma / (Fr * bath_w^2)) * L_unit * 1e4;
            plot(3.05 * L_unit, A_bath_um, '<r', 'MarkerSize', 8, 'MarkerFaceColor', 'r', 'Clipping', 'off');
            plot(3.05 * L_unit, -A_bath_um, '<r', 'MarkerSize', 8, 'MarkerFaceColor', 'r', 'Clipping', 'off');
  
            title(sprintf('   t = %0.3f ms, z_{lab} = %.2f cm', recordedConditions{current_index}.time*T_unit*1000, z_lab*L_unit),'FontSize',16);
            xlabel("x (cm)"); ylabel("y (\mum)");
            grid on
            hold off;
            xlim([-3 * L_unit, 3 * L_unit]);
            ylim([-PROBLEM_CONSTANTS.ylim_limit * L_unit * 1e4, PROBLEM_CONSTANTS.ylim_limit * L_unit * 1e4]);
            drawnow;
        end
     end 
    
    % Strip large / non-serializable fields before saving
    if isfield(PROBLEM_CONSTANTS, 'L_Library')
        PROBLEM_CONSTANTS = rmfield(PROBLEM_CONSTANTS, {'L_Library', 'U_Library', 'P_Library'});
    end
    if isfield(PROBLEM_CONSTANTS, 'luFutures')
        PROBLEM_CONSTANTS = rmfield(PROBLEM_CONSTANTS, 'luFutures');
    end

    for ii = 1:length(savingvarNames)
       variableValues{ii} = eval(savingvarNames{ii});
    end
    results_saver("simulationResults", 1:(current_index-1), variableValues, savingvarNames, NameValueArgs);

catch ME
    % Strip large / non-serializable fields before saving errored results
    if isfield(PROBLEM_CONSTANTS, 'L_Library')
        PROBLEM_CONSTANTS = rmfield(PROBLEM_CONSTANTS, {'L_Library', 'U_Library', 'P_Library'});
    end
    if isfield(PROBLEM_CONSTANTS, 'luFutures')
        PROBLEM_CONSTANTS = rmfield(PROBLEM_CONSTANTS, 'luFutures');
    end
    for ii = 1:length(savingvarNames)
       variableValues{ii} = eval(savingvarNames{ii}); 
    end
    results_saver("errored_results", 1:(current_index-1), variableValues, savingvarNames, NameValueArgs);
    fprintf("Couldn't run simulation"); 
    
    logDir = fullfile(fileparts(mfilename('fullpath')), '..', '..', '0_data', 'logs');
    if ~exist(logDir, 'dir'); mkdir(logDir); end
    save(fullfile(logDir, sprintf("error_log%s.mat", datetimeMarker)),'ME');
end 

% Return CoM trajectory in physical units if caller requested outputs
if nargout > 0
    valid  = 1:current_index;
    t_s    = cellfun(@(c) c.time * T_unit,           recordedConditions(valid));
    CoM_cm = cellfun(@(c) c.center_of_mass * L_unit, recordedConditions(valid));

    % Collect the full bath surface history (each column is one time step)
    % Size: [nr x steps]
    eta_history_cm = cell2mat(cellfun(@(c) c.bath_surface, recordedConditions(valid), 'UniformOutput', false)') * L_unit;
end

mypwd = split(pwd, "1_code"); mypwd = mypwd{2};
fprintf("Finished simulation on %s. Time elapsed: %0.2f minutes\n", mypwd, toc/60);
cd(currfold)
diary off;

end

function results_saver(fileName, indexes, variables, variableNames, NameValueArgs)
    currfold = pwd;
    % Root all results in 0_data/simulationResults/, never inside simulation_code/
    dataRoot = fullfile(fileparts(mfilename('fullpath')), '..', '..', '0_data', 'simulationResults');
    if ~exist(dataRoot, 'dir'); mkdir(dataRoot); end
    cd(dataRoot);
    folders = { ...
        sprintf("rho%.2fgcm3-sigma%.2fdynecm-nu%.4fSt", ...
        NameValueArgs.bathDensity, NameValueArgs.bathSurfaceTension, NameValueArgs.bathViscosity), ...
        sprintf("diskRadius%.2gcm-diskMass%.2gg", NameValueArgs.diskRadius, NameValueArgs.diskMass), ...
        sprintf("forceAmplitude%gdyne-forceFrequency%gHz", NameValueArgs.forceAmplitude, NameValueArgs.forceFrequency/(2*pi)), ...
        sprintf("bathAmplitude%.4gcm-bathFrequency%gHz-phaseDiff%gdeg", ...
        NameValueArgs.bathAmplitude, NameValueArgs.bathFrequency, NameValueArgs.phaseDifference)
    };
    for ii = 1:length(folders)
        folder = folders{ii};
        if ~exist(folder, "dir"); mkdir(folder);  end
        cd(folder);
    end
    if length(indexes) > 1 && indexes(2) == 1
        indexes = indexes(2:end); 
    end
    stru = struct();
    for ii = 1:length(variables)
       var = variables{ii};
       if max(size(var)) > 1
           if ~isstruct(var)
               var = var(indexes);
           elseif iscell(var)
               var = var{:, indexes};
           else
               var = var(:, indexes);
           end
       end
       stru.(variableNames{ii}) = var;
    end
    my_file = sprintf('%s%s.mat', fileName, NameValueArgs.datetimeMarker);
   if exist(my_file, 'file')
       save(my_file, '-struct', 'stru', '-append');
   else
       save(my_file, '-struct', 'stru');
   end
    cd(currfold)
end

function out = getVarName(var)
    out = inputname(1);
end
