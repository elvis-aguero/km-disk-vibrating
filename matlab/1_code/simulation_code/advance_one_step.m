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
    
    % Use periodic index to avoid floating-point phase drift over many cycles. CHANGED
    % cycleIdx is 1-based, matches the step within the cycle (1 to stepsPerCycle)
    cycleIdx = mod(current_step, PROBLEM_CONSTANTS.stepsPerCycle) + 1;
    
    % Re-calculating time-dependent terms precisely using periodic index
    g_prefactor = (1 - gamma * cos(bath_w * cycleIdx * dt + phase)); 
    force_term = dt * F * cos(w * cycleIdx * dt); 
    
    indep = [b; CoM_vel - dt/Fr * g_prefactor - force_term; CoM]; % CHANGED: Perfectly periodic RHS
    
    % --- Solver Logic ---
    if gamma == 0 && ~any(isnan(PROBLEM_CONSTANTS.precomputedInverse(:)))
        % Static Gravity: Use precomputed inverse
        sol = PROBLEM_CONSTANTS.precomputedInverse * indep;
    elseif PROBLEM_CONSTANTS.useCaching
        % Oscillating Gravity with Caching.
        if isempty(PROBLEM_CONSTANTS.InverseLibrary{cycleIdx})
            % Compute and cache
            Mat = buildSystemMatrix(PROBLEM_CONSTANTS, g_prefactor, dt, nr, cPoints, dr, SF);
            PROBLEM_CONSTANTS.InverseLibrary{cycleIdx} = inv(Mat);
        else
            % AUDIT: Compare cached vs fresh at the start of Cycle 2. CHANGED
            if current_step == PROBLEM_CONSTANTS.stepsPerCycle
                fprintf('\n--- NUMERICAL AUDIT: STEP %d (Start of Cycle 2) ---\n', current_step + 1);
                M_fresh = buildSystemMatrix(PROBLEM_CONSTANTS, g_prefactor, dt, nr, cPoints, dr, SF);
                P_cached = PROBLEM_CONSTANTS.InverseLibrary{cycleIdx};
                
                err_mat = norm(inv(M_fresh) - P_cached, 'fro') / norm(P_cached, 'fro');
                fprintf('Operator Periodicity Error (Frobenius): %.2e\n', err_mat);
                
                x_cached = P_cached * indep;
                x_direct = M_fresh \ indep;
                err_sol = norm(x_cached - x_direct) / norm(x_direct);
                fprintf('Solve Consistency Error: %.2e\n', err_sol);
                
                if err_mat > 1e-12 || err_sol > 1e-12
                    fprintf('VERDICT: BUG DETECTED in caching/periodicity logic.\n');
                else
                    fprintf('VERDICT: CACHING IS MATHEMATICALLY VALID. Explosion is purely cumulative.\n');
                end
                fprintf('---------------------------------------------------\n\n');
            end
        end
        sol = PROBLEM_CONSTANTS.InverseLibrary{cycleIdx} * indep;
    else



    next_condition = previous_conditions;
    % Re-mapping sol based on the system construction
    % sol = [eta_rest; phi; p; v; z]
    next_condition.bath_surface = [sol(end)* ones(cPoints, 1); sol(1:nr-cPoints)];
    next_condition.bath_potential = sol(nr-cPoints+1:2*nr-cPoints);
    next_condition.pressure = sol(2*nr-cPoints+1:2*nr);
    next_condition.center_of_mass_velocity = sol(end-1);
    next_condition.center_of_mass = sol(end);
    next_condition.time = next_condition.time + dt;
    
    % --- Numerical Guardrail. CHANGED ---
    if any(isnan(sol)) || any(isinf(sol)) || abs(next_condition.center_of_mass) > 10
        error('Numerical Divergence Detected: CoM = %.2e. Terminating simulation.', next_condition.center_of_mass);
    end
    
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