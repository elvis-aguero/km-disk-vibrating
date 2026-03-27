function solve_motion(NameValueArgs)
% SOLVE_MOTION Simulates the motion of a disk oscillating in a fluid bath.
%   This function models the motion of a disk subjected to sinusoidal forces
%   within a fluid environment using configurable parameters provided through
%   name-value arguments.

% Parse input arguments using name-value pairs and set default values
arguments
    % Unless stated otherwise, all units are cgs. 
    % Default values amount to a water bath
    NameValueArgs.diskRadius (1, 1) double = .5 % Radius of the oscillating disk, in cm
    NameValueArgs.diskMass (1, 1) double = 1 % Mass in grams of the disk
    NameValueArgs.forceAmplitude (1, 1) double = 0 % Amplitude of sinusoidal force applied to disk (in dynes)
    NameValueArgs.forceFrequency (1, 1) double = 90 % Frequency of sinusoidal force in Hz
    NameValueArgs.bathAmplitude (1, 1) double = 0 % Amplitude of bath oscillation in cm. 
    NameValueArgs.bathFrequency (1, 1) double = 90 % Frequency of bath oscillation in Hz. 
    NameValueArgs.phaseDifference (1, 1) double = -90 % Phase difference between disk forcing and bath oscillation in degrees.
    NameValueArgs.bathDensity (1, 1) double = 1 % Density of bath's fluid in g/cm^3
    NameValueArgs.bathSurfaceTension (1, 1) double = 72.20 % For water, in dynes/cm
    NameValueArgs.bathViscosity (1, 1) double = 0.978e-2 % Viscosity in Stokes (cgs)
    NameValueArgs.g (1, 1) double = 981 % Gravitational constant in cgs
    NameValueArgs.bathDiameter (1, 1) double = 100 % Diameter of the bath wrt to disk Radius
    NameValueArgs.spatialResolution (1, 1) double = 50 % Number of numerical radial intervals in one disk radius
    NameValueArgs.temporalResolution (1, 1) double = 20; % Number of temporal steps in one adimensional unit
    NameValueArgs.simulationTime (1, 1) double = 10/90; % Time to be simulated in seconds
    NameValueArgs.debug_flag (1, 1) logical = true; % To show some debugging info
    NameValueArgs.forceNoCaching (1, 1) logical = false; % Manual override to disable RAM caching
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
if debug_flag == true; diary(sprintf("../0_data/manual/Logger/solve_motion%s.txt", datetimeMarker)); end

% Set up folder for data saving
fold = fullfile(pwd, sprintf("D%dQuant%d", spatialResolution, bathDiameter));
try
    % Navigate to the folder and load necessary matrices
    cd(fold)
    nr = ceil(spatialResolution * bathDiameter / 2);
    load(sprintf('DTNnew345nr%dD%drefp10.mat', nr, bathDiameter), 'DTNnew345')
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
H_limit_adim = 1.5 * (A_bath_adim + A_forcing_adim + 0.5); % Symmetric limit around z=0

% Steps per cycle must be an integer for periodic caching to be valid. 
isPeriodic = (abs(freq_adim - bath_freq_adim) < 1e-6);
effective_w_adim = bath_freq_adim;
stepsPerCycle = round((2 * pi / effective_w_adim) * temporalResolution); 
dt = (2 * pi / effective_w_adim) / stepsPerCycle; % Adjusted adimensional time step

steps = ceil(simulationTime / (dt * T_unit)); % Minimum number of time steps

% Warn if memory usage could be large
if steps * nr * 8 > 1e+9
    warning('Spatial resolution and simulation times might be too big to store all matrices in memory');
end

%Inintial conditions for the fluid
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

% --- Smart Gatekeeper for Caching ---
availableRAM = getAvailableRAM();
useCaching = false;
InverseLibrary = {};

if NameValueArgs.forceNoCaching
    fprintf('Smart Caching: DISABLED (Manual override via forceNoCaching)\n');
elseif isPeriodic && bath_forcing_amplitude ~= 0
    requiredRAM = stepsPerCycle * ((2*nr+2)^2) * 8;
    if requiredRAM < 0.75 * availableRAM
        useCaching = true;
        InverseLibrary = cell(stepsPerCycle, 1);
        fprintf('Smart Caching: ENABLED (Estimated RAM: %.2f GB)\n', requiredRAM/1e9);
    elseif requiredRAM > availableRAM
        fprintf('Smart Caching: DISABLED (Insufficient RAM: %.2f GB required, %.2f GB available)\n', requiredRAM/1e9, availableRAM/1e9);
        warning('Oscillating bath with no caching will be slow. Using preconditioned iterative solver.');
    else
        fprintf('Caching requires %.2f GB (Available: %.2f GB).\n', requiredRAM/1e9, availableRAM/1e9);
        if ~usejava('desktop') || isempty(javachk('desktop'))
             fprintf('Smart Caching: DISABLED (Headless mode)\n');
        else
             useCacheStr = input('Enable Caching? (y/n) [n]: ', 's');
             if strcmpi(useCacheStr, 'y')
                useCaching = true;
                InverseLibrary = cell(stepsPerCycle, 1);
             end
        end
    end
end

% --- Preconditioner for Iterative Solver ---
if ~useCaching && (bath_forcing_amplitude ~= 0)
    fprintf('Precomputing frozen preconditioner... ');
    tic;
    if ~any(isnan(precomputedInverse(:)))
        preconditioner = precomputedInverse;
    else
        M_static = buildStaticMatrix(Fr, We, Re, dr, laplacian, DTN, pressureIntegral(spatialResolution+1, :), obj_mass_adim, nr, spatialResolution+1, dt, surface_force_adim);
        preconditioner = inv(M_static);
    end
    fprintf('Done (%.2fs).\n', toc);
else
    preconditioner = nan;
end

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
    "useCaching", useCaching, "InverseLibrary", {InverseLibrary}, "stepsPerCycle", stepsPerCycle, ...
    "ylim_limit", H_limit_adim, "step_counter", 0, ...
    "preconditioner", preconditioner);

fprintf("Starting simulation on %s\n", pwd);

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
            xplot = dr*(0:nr-1);
            plot([-fliplr(xplot(2:width)), xplot(1:width)], eta + zb, 'b', 'Linewidth', 2);
            hold on
            
            x = [1, 1, -1, -1];
            z = recordedConditions{current_index}.center_of_mass;
            z_lab = z + zb;
            y = [z_lab, z_lab+1/10, z_lab+1/10, z_lab];
            fill(x, y, 'k');
  
            title(sprintf('   t = %0.3f s, z_{lab} = %.2f', recordedConditions{current_index}.time*T_unit, z_lab*L_unit),'FontSize',16);
            grid on
            hold off;
            xlim([-3, 3]);
            ylim([-PROBLEM_CONSTANTS.ylim_limit, PROBLEM_CONSTANTS.ylim_limit]);
            drawnow;
        end
     end 
    
    if ~any(isnan(PROBLEM_CONSTANTS.precomputedInverse(:)))
        cd(fold)
        precomputedInverse = PROBLEM_CONSTANTS.precomputedInverse;
        save(myfile, "precomputedInverse")
        PROBLEM_CONSTANTS = rmfield(PROBLEM_CONSTANTS, "precomputedInverse");
    end

    for ii = 1:length(savingvarNames)
       variableValues{ii} = eval(savingvarNames{ii}); 
    end
    results_saver("simulationResults", 1:(current_index-1), variableValues, savingvarNames, NameValueArgs);

catch ME
    if ~any(isnan(PROBLEM_CONSTANTS.precomputedInverse(:)))
        cd(fold)
        precomputedInverse = PROBLEM_CONSTANTS.precomputedInverse;
        save(myfile, "precomputedInverse")
        PROBLEM_CONSTANTS = rmfield(PROBLEM_CONSTANTS, "precomputedInverse");
    end
    for ii = 1:length(savingvarNames)
       variableValues{ii} = eval(savingvarNames{ii}); 
    end
    results_saver("errored_results", 1:(current_index-1), variableValues, savingvarNames, NameValueArgs);
    fprintf("Couldn't run simulation"); 
    save(sprintf("error_log%s.mat", datetimeMarker),'ME');
end 

mypwd = split(pwd, "1_code"); mypwd = mypwd{2};
fprintf("Finished simulation on %s. Time elapsed: %0.2f minutes\n", mypwd, toc/60);
cd(currfold)
diary off;

end

function results_saver(fileName, indexes, variables, variableNames, NameValueArgs)
    currfold = pwd;
    folders = { ...
        sprintf("rho%.2fgcm3-sigma%.2fdynecm-nu%.4fSt", ...
        NameValueArgs.bathDensity, NameValueArgs.bathSurfaceTension, NameValueArgs.bathViscosity), ...
        sprintf("diskRadius%.2gcm-diskMass%.2gg", NameValueArgs.diskRadius, NameValueArgs.diskMass), ...
        sprintf("forceAmplitude%gdyne-forceFrequency%gHz", NameValueArgs.forceAmplitude, NameValueArgs.forceFrequency/(2*pi))
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

function Mat = buildStaticMatrix(Fr, We, Re, dr, Delta, DTN, pIntegral, Ma, nr, cPoints, dt, SF)
    Sist = [[eye(nr)-dt*2*Delta/Re,-dt*DTN];...
        [dt*(eye(nr)/Fr - Delta/We),eye(nr)-dt*2*Delta/Re]]; 
    Mat =  [[Sist(:,(cPoints+1):2*nr),...
        [zeros(nr,cPoints);dt*eye(cPoints);zeros(nr-cPoints,cPoints)],...
        zeros(2*nr,1),Sist(:,1:cPoints)*ones(cPoints,1)];
        [-SF*dt/dr, zeros(1,2*nr-cPoints-1),-dt*pIntegral(1:cPoints)/Ma, 1 , SF*dt/dr];
        [zeros(1,2*nr-cPoints),-zeros(1, cPoints)  ,-dt,1]];
end
