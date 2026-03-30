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
    
    % Discrete Step Counter
    current_step = PROBLEM_CONSTANTS.step_counter;
    PROBLEM_CONSTANTS.step_counter = current_step + 1;
    cycleIdx = mod(current_step, PROBLEM_CONSTANTS.stepsPerCycle) + 1;
    
    % Physics Evaluation (Using global counter for continuity)
    g_prefactor = (1 - gamma * cos(bath_w * (current_step + 1) * dt + phase)); 
    force_term = dt * F * cos(w * (current_step + 1) * dt); 
    
    % Standard Master Branch indep vector
    indep = [b; CoM_vel - dt/Fr * g_prefactor - force_term; CoM];
    
    % --- Solver Logic ---
    if gamma == 0 && ~any(isnan(PROBLEM_CONSTANTS.precomputedInverse(:)))
        sol = PROBLEM_CONSTANTS.precomputedInverse * indep;
    elseif PROBLEM_CONSTANTS.solverType == "gmres"
        Mat = buildSystemMatrix(PROBLEM_CONSTANTS, g_prefactor, dt, nr, cPoints, dr, SF);
        [sol, flag] = gmres(Mat, indep, [], PROBLEM_CONSTANTS.gmresTolerance, size(Mat,1), [], [], PROBLEM_CONSTANTS.gmres_x0);
        if flag ~= 0
            warning('advance_one_step:gmresNoConverge', ...
                'GMRES did not converge at step %d (flag=%d)', current_step, flag);
        end
        PROBLEM_CONSTANTS.gmres_x0 = sol;
    elseif PROBLEM_CONSTANTS.solverType == "lu"
        if isempty(PROBLEM_CONSTANTS.L_Library{cycleIdx})
            % Fetch async future if available (fired by solve_motion before the loop)
            if ~isempty(PROBLEM_CONSTANTS.luFutures) && ~isempty(PROBLEM_CONSTANTS.luFutures{cycleIdx})
                [L, U, P] = fetchOutputs(PROBLEM_CONSTANTS.luFutures{cycleIdx});
                PROBLEM_CONSTANTS.luFutures{cycleIdx} = []; % release handle
            else
                % Lazy-sync fallback: no PCT available
                Mat = buildSystemMatrix(PROBLEM_CONSTANTS, g_prefactor, dt, nr, cPoints, dr, SF);
                [L, U, P] = lu(Mat);
            end
            PROBLEM_CONSTANTS.L_Library{cycleIdx} = L;
            PROBLEM_CONSTANTS.U_Library{cycleIdx} = U;
            PROBLEM_CONSTANTS.P_Library{cycleIdx} = P;
        end
        L = PROBLEM_CONSTANTS.L_Library{cycleIdx};
        U = PROBLEM_CONSTANTS.U_Library{cycleIdx};
        P = PROBLEM_CONSTANTS.P_Library{cycleIdx};
        sol = U \ (L \ (P * indep));
    end

    % Standard Master Branch Mapping
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

