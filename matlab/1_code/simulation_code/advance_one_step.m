function [next_condition, PROBLEM_CONSTANTS] = advance_one_step(previous_conditions, PROBLEM_CONSTANTS)

    dt = previous_conditions.dt;
    dr = PROBLEM_CONSTANTS.dr;
    SF = PROBLEM_CONSTANTS.surface_force_constant;
    Fr = PROBLEM_CONSTANTS.froude;
    F = PROBLEM_CONSTANTS.force_amplitude;
    w = PROBLEM_CONSTANTS.force_frequency;
    t = previous_conditions.time;
    nr = PROBLEM_CONSTANTS.nr;
    cPoints = PROBLEM_CONSTANTS.contact_points;
    CoM = previous_conditions.center_of_mass;
    CoM_vel = previous_conditions.center_of_mass_velocity;
    b = [previous_conditions.bath_surface; previous_conditions.bath_potential];
    %b = b-Sist(:,1:cPoints)*(CoM * ones(cPoints, 1)); %zs(1:cPoints)+Rv);
    gamma = PROBLEM_CONSTANTS.bath_forcing_amplitude; % CHANGED
    bath_w = PROBLEM_CONSTANTS.bath_frequency; % CHANGED
    phase = PROBLEM_CONSTANTS.phase_difference; % CHANGED
    
    % Use discrete step counter for more robust physics and indexing. CHANGED
    current_step = PROBLEM_CONSTANTS.step_counter;
    PROBLEM_CONSTANTS.step_counter = current_step + 1;
    
    % Re-calculating time-dependent term precisely to match cache indices
    % (Implicit Euler: evaluation at t + dt)
    g_prefactor = (1 - gamma * cos(bath_w * (current_step + 1) * dt + phase)); 
    
    indep = [b; CoM_vel - dt/Fr * g_prefactor - dt*F*cos(w*(current_step + 1) * dt); CoM]; % CHANGED: Using counter-based time
    
    % --- Solver Logic ---
    if gamma == 0 && ~any(isnan(PROBLEM_CONSTANTS.precomputedInverse(:)))
        % Static Gravity: Use precomputed inverse
        sol = PROBLEM_CONSTANTS.precomputedInverse * indep;
    elseif PROBLEM_CONSTANTS.useCaching
        % Oscillating Gravity with Caching. CHANGED
        cycleIdx = mod(current_step, PROBLEM_CONSTANTS.stepsPerCycle) + 1;
        if isempty(PROBLEM_CONSTANTS.InverseLibrary{cycleIdx})
            % Compute and cache
            Mat = buildSystemMatrix(PROBLEM_CONSTANTS, g_prefactor, dt, nr, cPoints, dr, SF);
            PROBLEM_CONSTANTS.InverseLibrary{cycleIdx} = inv(Mat);
        end
        sol = PROBLEM_CONSTANTS.InverseLibrary{cycleIdx} * indep;
    else
        % Oscillating Gravity without Caching: Iterative Solver with Warm Start. CHANGED
        Mat = buildSystemMatrix(PROBLEM_CONSTANTS, g_prefactor, dt, nr, cPoints, dr, SF);

        % Construct initial guess from previous condition (warm start)
        % Structure: [eta_rest; phi; p; v; z]
        eta_rest = previous_conditions.bath_surface(cPoints+1:nr);
        phi = previous_conditions.bath_potential;
        p = previous_conditions.pressure;
        v = previous_conditions.center_of_mass_velocity;
        z = previous_conditions.center_of_mass;
        x0 = [eta_rest; phi; p; v; z];

        % GMRES Solve. CHANGED: Tighter tolerance for numerical agreement
        [sol, ~] = gmres(Mat, indep, [], 1e-12, 100, [], [], x0); 
    end


    next_condition = previous_conditions;
    % Re-mapping sol based on the system construction
    % sol = [eta_rest; phi; p; v; z]
    next_condition.bath_surface = [sol(end)* ones(cPoints, 1); sol(1:nr-cPoints)];
    next_condition.bath_potential = sol(nr-cPoints+1:2*nr-cPoints);
    next_condition.pressure = sol(2*nr-cPoints+1:2*nr);
    next_condition.center_of_mass_velocity = sol(end-1);
    next_condition.center_of_mass = sol(end);
    next_condition.time = next_condition.time + dt;
    
end

    function Mat = buildSystemMatrix(PROBLEM_CONSTANTS, g_prefactor, dt, nr, cPoints, dr, SF)
    We = PROBLEM_CONSTANTS.weber;
    Re = PROBLEM_CONSTANTS.reynolds;
    Delta = PROBLEM_CONSTANTS.laplacian;
    DTN = PROBLEM_CONSTANTS.DTN;
    Fr = PROBLEM_CONSTANTS.froude;
    pIntegral = PROBLEM_CONSTANTS.pressure_integral;
    Ma = PROBLEM_CONSTANTS.obj_mass;

    % Preparing the matrix (2x2 block)
    Sist = [[eye(nr)-dt*2*Delta/Re,-dt*DTN];...
        [dt*(eye(nr)/Fr * g_prefactor - Delta/We),eye(nr)-dt*2*Delta/Re]]; 

    % Completing the system
    Mat =  [[Sist(:,(cPoints+1):2*nr),...
        [zeros(nr,cPoints);dt*eye(cPoints);zeros(nr-cPoints,cPoints)],...
        zeros(2*nr,1),Sist(:,1:cPoints)*ones(cPoints,1)];
        [-SF*dt/dr, zeros(1,2*nr-cPoints-1),-dt*pIntegral(1:cPoints)/Ma, 1 , SF*dt/dr];
        [zeros(1,2*nr-cPoints),-zeros(1, cPoints)  ,-dt,1]];
    end