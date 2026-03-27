function solve_motion(NameValueArgs)
% SOLVE_MOTION Simulates the motion of a disk oscillating in a fluid bath.
%   This functi on models the motion of a disk subjected to sinusoidal forces
%   within a fluid environment using configurable parameters provided through
%   name-value arguments.
%
%   Parameters (Name-Value Arguments):
%       - diskRadius: Radius of the oscillating disk.
%       - diskMass: Mass of the disk in grams.
%       - forceAmplitude: Amplitude of the sinusoidal force.
%       - forceFrequency: Frequency of the applied force.
%       - bathDensity: Density of the fluid bath.
%       - bathSurfaceTension: Surface tension of the fluid.
%       - bathViscosity: Viscosity of the fluid.
%       - g: Gravitational constant.
%       - bathDiameter: Diameter of the bath relative to disk radius.
%       - spatialResolution: Number of radial intervals per disk radius.
%       - temporalResolution: Number of time steps per dimensionless unit.
%       - simulationTime: Total simulation time in seconds.
%       - debug_flag: Enables/disables debugging information.
%
%   The script simulates the fluid motion and surface interaction with the disk, and
%   provides visualizations and logs. The simulation considers various dimensionless
%   numbers (Reynolds, Weber, and Froude numbers) and handles precomputed inverse matrices.
%
%   Outputs:
%       Results are saved in specified directories based on the provided parameters.
%
%   Example:
%       solve_motion('diskRadius', 1, 'diskMass', 2, 'forceAmplitude', 500);
%
%   Dependencies:
%       This function depends on external functions such as domainMaker and 
%       advance_one_step to create the fluid domain and advance the simulation, 
%       respectively.
%
%   Author:
%       Elvis Aguero - Date of Creation: Oct 10th, 2024.
%
%   See also: domainMaker, advance_one_step, results_saver

% Parse input arguments using name-value pairs and set default values

arguments
    % Unless stated otherwise, all units are cgs. 
    % Default values amount to a water bath
    NameValueArgs.diskRadius (1, 1) double = .5 % Radius of the oscillating disk, in cm
    NameValueArgs.diskMass (1, 1) double = 1 % Mass in grams of the disk
    NameValueArgs.forceAmplitude (1, 1) double = 0 % Amplitude of sinusoidal force applied to disk (in dynes)
    NameValueArgs.forceFrequency (1, 1) double = 90 % Frequency of sinusoidal force in Hz
    NameValueArgs.bathAmplitude (1, 1) double = 0.1 % Amplitude of bath oscillation in cm. CHANGED
    NameValueArgs.bathFrequency (1, 1) double = 90 % Frequency of bath oscillation in Hz. CHANGED
    NameValueArgs.phaseDifference (1, 1) double = -90 % Phase difference between disk forcing and bath oscillation in degrees. CHANGED to -90 for sine profile
    NameValueArgs.bathDensity (1, 1) double = 1 % Density of bath's fluid in g/cm^3
    NameValueArgs.bathSurfaceTension (1, 1) double = 72.20 % For water, in dynes/cm
    NameValueArgs.bathViscosity (1, 1) double = 0.978e-2 % Viscosity in Stokes (cgs)
    NameValueArgs.g (1, 1) double = 981 % Gravitational constant in cgs
    NameValueArgs.bathDiameter (1, 1) double = 100 % Diameter of the bath wrt to disk Radius
    NameValueArgs.spatialResolution (1, 1) double = 50 % Number of numerical radial intervals in one disk radius
    NameValueArgs.temporalResolution (1, 1) double = 20; % Number of temporal steps in one adimensional unit
    NameValueArgs.simulationTime (1, 1) double = 10/90; % Time to be simulated in seconds
    NameValueArgs.debug_flag (1, 1) logical = true; % To show some debugging info
    NameValueArgs.forceNoCaching (1, 1) logical = false; % CHANGED: Manual override to disable RAM caching
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
currfold = fullfile(fileparts(mfilename('fullpath'))); %fileparts(matlab.desktop.editor.getActiveFilename);
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

% New bath-specific parameters. CHANGED
bath_angular_freq = NameValueArgs.bathFrequency * 2 * pi;
bath_freq_adim = bath_angular_freq * T_unit;
phase_diff_rad = NameValueArgs.phaseDifference * pi / 180;
bath_forcing_amplitude = NameValueArgs.bathAmplitude * bath_angular_freq^2 / g; % Dimensionless acceleration ratio (A*w^2/g)

% Calculate initial velocity offset for zero lab velocity. CHANGED
% z_b = A*cos(w*t + phase). v_b = -A*w*sin(w*t + phase). 
% At t=0: v_b(0) = -A*w*sin(phase). 
% We want v_lab = v_disk_bath + v_b = 0 -> v_disk_bath = -v_b(0) = A*w*sin(phase)
v_bath_0_adim = (NameValueArgs.bathAmplitude * bath_angular_freq * sin(phase_diff_rad)) / V_unit;

% Calculate informed static ylim. CHANGED
A_bath_adim = abs(NameValueArgs.bathAmplitude / L_unit);
A_forcing_adim = force_adim / (freq_adim^2 + 1e-6); % Rough estimate
H_limit_adim = 1.5 * (A_bath_adim + A_forcing_adim + 0.5); % Symmetric limit around z=0

% Set numerical simulation parameters
dt = 1 / temporalResolution; % Adimensional time step
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
    "center_of_mass", 0, "center_of_mass_velocity", v_bath_0_adim, ... % CHANGED: Starts at v_lab = 0
    "bath_surface", etaInitial, "bath_potential", phiInitial, "pressure", zeros(spatialResolution+1, 1));
current_index = 1; %iteration counter
recordedConditions = cell(steps, 1);
recordedConditions{current_index} = current_conditions;

% --- Smart Gatekeeper for Caching ---
% Logic: The system is periodic if frequencies are commensurate. 
% If one amplitude is zero, the system follows the period of the other.
if bath_forcing_amplitude == 0
    % Case 1: Static Bath. System is constant (handled by precomputedInverse)
    % but we set periodicity logic for completeness.
    isPeriodic = true;
    effective_w_adim = freq_adim;
elseif forceAmplitude == 0
    % Case 2: Only Bath oscillates. Period is defined by bath frequency.
    isPeriodic = true;
    effective_w_adim = bath_freq_adim;
else
    % Case 3: Both active. Check if commensurate (simple check for equality)
    isPeriodic = (abs(freq_adim - bath_freq_adim) < 1e-6);
    effective_w_adim = bath_freq_adim; % Use bath freq as reference
    end

    % Steps per cycle must be an integer for periodic caching to be valid. CHANGED
    % We adjust dt slightly to ensure stepsPerCycle * dt = 2*pi / effective_w_adim.
    % This is moved outside the Gatekeeper to ensure consistency between runs.
    stepsPerCycle = round((2 * pi / effective_w_adim) * temporalResolution); 
    dt = (2 * pi / effective_w_adim) / stepsPerCycle; % Adjusted adimensional time step

    % --- Smart Gatekeeper for Caching ---
    systemSize = 2 * nr + 2;
    requiredRAM = stepsPerCycle * (systemSize^2) * 8;
    availableRAM = getAvailableRAM();

    useCaching = false;
    InverseLibrary = {};

    if NameValueArgs.forceNoCaching
    fprintf('Smart Caching: DISABLED (Manual override via forceNoCaching)\n');
    elseif isPeriodic && bath_forcing_amplitude ~= 0
    if requiredRAM < 0.75 * availableRAM

        InverseLibrary = cell(stepsPerCycle, 1);
        fprintf('Smart Caching: ENABLED (Estimated RAM: %.2f GB)\n', requiredRAM/1e9);
    elseif requiredRAM > availableRAM
        fprintf('Smart Caching: DISABLED (Insufficient RAM: %.2f GB required, %.2f GB available)\n', requiredRAM/1e9, availableRAM/1e9);
        warning('Oscillating bath with no caching will be slow. Using iterative solver.');
    else
        % Margin case: Ask user if interactive
        useCacheStr = 'n';
        if ~usejava('desktop') || isempty(javachk('desktop'))
             fprintf('Smart Caching: DISABLED (Headless mode, skipping RAM intensive cache)\n');
        else
             % Creating a simple timer-based input is complex in a single block, 
             % we'll use a simplified version for MATLAB CLI.
             fprintf('Caching requires %.2f GB (Available: %.2f GB).\n', requiredRAM/1e9, availableRAM/1e9);
             useCacheStr = input('Enable Caching? (y/n) [n]: ', 's');
        end
        if strcmpi(useCacheStr, 'y')
            useCaching = true;
            InverseLibrary = cell(stepsPerCycle, 1);
        end
    end
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
    "ylim_limit", H_limit_adim, "step_counter", 0); % CHANGED: added step_counter

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
            width = 6; % CHANGED: Increased to show more of the surface
            width = min(nr, ceil(width*spatialResolution));
            
            % Calculate current bath position for lab frame plotting. CHANGED
            gamma = PROBLEM_CONSTANTS.bath_forcing_amplitude;
            Fr = PROBLEM_CONSTANTS.froude;
            bath_w = PROBLEM_CONSTANTS.bath_frequency;
            phase = PROBLEM_CONSTANTS.phase_difference;
            zb = (gamma / (Fr * bath_w^2)) * cos(bath_w * recordedConditions{current_index}.time + phase);
            
            eta = recordedConditions{current_index}.bath_surface;
            eta = [flipud(eta(2:width)); eta(1:width)];
            xplot = dr*(0:nr-1);
            plot([-fliplr(xplot(2:width)), xplot(1:width)], eta + zb, 'b', 'Linewidth', 2); % CHANGED: added zb
            hold on
            
            % Filling the shape of the vibrating object
            x = [1, 1, -1, -1];
            z = recordedConditions{current_index}.center_of_mass;
            z_lab = z + zb; % Disk bottom in lab frame. CHANGED
            y = [z_lab, z_lab+1/10, z_lab+1/10, z_lab]; % CHANGED
            fill(x, y, 'k');
  
            %axis equal
            title(sprintf('   t = %0.3f s, z_{lab} = %.2f', recordedConditions{current_index}.time*T_unit, z_lab*L_unit),'FontSize',16); % CHANGED
            grid on
            hold off;
            xlim([-3, 3]); % Showing 6 disk radii
            ylim([-PROBLEM_CONSTANTS.ylim_limit, PROBLEM_CONSTANTS.ylim_limit]); % CHANGED: Static symmetric limit
            %set(gca,'xlim',[-6 6])
            drawnow;
        end
    % 
  
     end % Outer while
    % 
    % 
    if ~isnan(PROBLEM_CONSTANTS.precomputedInverse)
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

    if ~isnan(PROBLEM_CONSTANTS.precomputedInverse)
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
end % end while catch

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

    if length(indexes) > 1 && indexes(2) == 1 % CHANGED: added length check
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
       stru.(variableNames{ii}) = var;% = struct(variableNames{ii}, var);
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
