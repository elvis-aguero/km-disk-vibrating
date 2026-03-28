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
    
    gamma = PROBLEM_CONSTANTS.bath_forcing_amplitude;
    bath_w = PROBLEM_CONSTANTS.bath_frequency;
    phase = PROBLEM_CONSTANTS.phase_difference;
    
    % Use discrete step counter for more robust physics and indexing.
    current_step = PROBLEM_CONSTANTS.step_counter;
    PROBLEM_CONSTANTS.step_counter = current_step + 1;
    
    % Use periodic index to avoid floating-point phase drift over many cycles.
    cycleIdx = mod(current_step, PROBLEM_CONSTANTS.stepsPerCycle) + 1;
    
    % Re-calculating time-dependent terms precisely using periodic index
    g_prefactor = (1 - gamma * cos(bath_w * cycleIdx * dt + phase)); 
    force_term = dt * F * cos(w * cycleIdx * dt); 
    
    indep = [b; CoM_vel - dt/Fr * g_prefactor - force_term; CoM];
    
    % --- Solver Logic ---
    if gamma == 0 && ~any(isnan(PROBLEM_CONSTANTS.precomputedInverse(:)))
        % Static Gravity: Use precomputed inverse
        sol = PROBLEM_CONSTANTS.precomputedInverse * indep;
    elseif PROBLEM_CONSTANTS.useCaching
        % Oscillating Gravity with LU Caching.
        if isempty(PROBLEM_CONSTANTS.L_Library{cycleIdx})
            % FIRST CYCLE: Build and cache LU
            Mat = buildSystemMatrix(PROBLEM_CONSTANTS, g_prefactor, dt, nr, cPoints, dr, SF);
            [L, U, P] = lu(Mat);
            PROBLEM_CONSTANTS.L_Library{cycleIdx} = L;
            PROBLEM_CONSTANTS.U_Library{cycleIdx} = U;
            PROBLEM_CONSTANTS.P_Library{cycleIdx} = P;
            PROBLEM_CONSTANTS.DiagnosticMat{cycleIdx} = Mat;
            
            % Apply LU factors
            sol = U \ (L \ (P * indep));
        else
            % DOUBLE-SOLVE AUDIT: Compare fresh vs cached solve. CHANGED
            Mat_fresh = buildSystemMatrix(PROBLEM_CONSTANTS, g_prefactor, dt, nr, cPoints, dr, SF);
            sol_fresh = Mat_fresh \ indep;
            
            L = PROBLEM_CONSTANTS.L_Library{cycleIdx};
            U = PROBLEM_CONSTANTS.U_Library{cycleIdx};
            P = PROBLEM_CONSTANTS.P_Library{cycleIdx};
            sol_cached = U \ (L \ (P * indep));
            
            err_sol = norm(sol_fresh - sol_cached) / norm(sol_fresh);
            
            if current_step < PROBLEM_CONSTANTS.stepsPerCycle + 5
                fprintf('Step %d: Solve Difference (Fresh Direct vs Cached LU) = %.2e\n', current_step + 1, err_sol);
            end
            
            sol = sol_cached; % Use cached for the remainder of the test
        end
    else
        % Forced Fresh Direct Solve. 
        Mat = buildSystemMatrix(PROBLEM_CONSTANTS, g_prefactor, dt, nr, cPoints, dr, SF);
        sol = Mat \ indep;
    end

    next_condition = previous_conditions;
    next_condition.bath_surface = [sol(end)* ones(cPoints, 1); sol(1:nr-cPoints)];
    next_condition.bath_potential = sol(nr-cPoints+1:2*nr-cPoints);
    next_condition.pressure = sol(2*nr-cPoints+1:2*nr);
    next_condition.center_of_mass_velocity = sol(end-1);
    next_condition.center_of_mass = sol(end);
    next_condition.time = next_condition.time + dt;
    
    % --- Numerical Guardrail ---
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

    Sist = [[eye(nr)-dt*2*Delta/Re,-dt*DTN];...
        [dt*(eye(nr)/Fr * g_prefactor - Delta/We),eye(nr)-dt*2*Delta/Re]]; 

    Mat =  [[Sist(:,(cPoints+1):2*nr),...
        [zeros(nr,cPoints);dt*eye(cPoints);zeros(nr-cPoints,cPoints)],...
        zeros(2*nr,1),Sist(:,1:cPoints)*ones(cPoints,1)];
        [-SF*dt/dr, zeros(1,2*nr-cPoints-1),-dt*pIntegral(1:cPoints)/Ma, 1 , SF*dt/dr];
        [zeros(1,2*nr-cPoints),-zeros(1, cPoints)  ,-dt,1]];
end
