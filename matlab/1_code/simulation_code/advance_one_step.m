function [next_condition, PROBLEM_CONSTANTS] = advance_one_step(previous_conditions, PROBLEM_CONSTANTS)

    dt = previous_conditions.dt;
    dr = PROBLEM_CONSTANTS.dr;
    SF = PROBLEM_CONSTANTS.surface_force_constant;
    Fr = PROBLEM_CONSTANTS.froude;
    We = PROBLEM_CONSTANTS.weber;
    Re = PROBLEM_CONSTANTS.reynolds;
    Delta = PROBLEM_CONSTANTS.laplacian;
    DTN = PROBLEM_CONSTANTS.DTN;
    pIntegral = PROBLEM_CONSTANTS.pressure_integral;
    Ma = PROBLEM_CONSTANTS.obj_mass;
    F = PROBLEM_CONSTANTS.force_amplitude;
    w = PROBLEM_CONSTANTS.force_frequency;
    
    nr = PROBLEM_CONSTANTS.nr;
    cPoints = PROBLEM_CONSTANTS.contact_points;
    
    gamma = PROBLEM_CONSTANTS.bath_forcing_amplitude;
    bath_w = PROBLEM_CONSTANTS.bath_frequency;
    phase = PROBLEM_CONSTANTS.phase_difference;
    
    current_step = PROBLEM_CONSTANTS.step_counter;
    PROBLEM_CONSTANTS.step_counter = current_step + 1;
    cycleIdx = mod(current_step, PROBLEM_CONSTANTS.stepsPerCycle) + 1;
    
    % Physics Timing
    g_prefactor = (1 - gamma * cos(bath_w * (current_step + 1) * dt + phase)); 
    force_term = dt * F * cos(w * (current_step + 1) * dt); 
    
    % --- Standardized System Structure ---
    % Variables: sol = [eta(cPoints+1:nr); phi(1:nr); p(1:cPoints); v; z]
    % Total Size: (nr-cPoints) + nr + cPoints + 1 + 1 = 2*nr + 2
    
    % Independent Vector (RHS)
    % 1. Kinematic Eq (Fluid surface outside disk): eta_next = eta_prev + dt*DTN*phi_next
    %    In matrix form: I*eta_next - dt*DTN_ext*phi_next = eta_prev
    rhs_eta = previous_conditions.bath_surface(cPoints+1:nr);
    
    % 2. Bernoulli Eq (Fluid potential): phi_next + dt/Fr*g_eff*eta_next - dt/We*Delta*eta_next = phi_prev
    rhs_phi = previous_conditions.bath_potential;
    
    % 3. Newton Eq 1 (Velocity): v_next + dt/Ma*Int(p) = v_prev - dt/Fr*g_eff - force_term + SF_term
    %    Wait, SF term is SF*dt/dr * (eta(cPoints+1) - z)
    rhs_v = previous_conditions.center_of_mass_velocity - dt/Fr * g_prefactor - force_term;
    
    % 4. Newton Eq 2 (Position): z_next - dt*v_next = z_prev
    rhs_z = previous_conditions.center_of_mass;
    
    indep = [rhs_eta; rhs_phi; rhs_v; rhs_z];
    
    % --- Solver Logic ---
    if gamma == 0 && ~any(isnan(PROBLEM_CONSTANTS.precomputedInverse(:)))
        sol = PROBLEM_CONSTANTS.precomputedInverse * indep;
    elseif PROBLEM_CONSTANTS.useCaching
        if isempty(PROBLEM_CONSTANTS.L_Library{cycleIdx})
            Mat = buildSystemMatrix(PROBLEM_CONSTANTS, g_prefactor, dt, nr, cPoints, dr, SF);
            [L, U, P] = lu(Mat);
            PROBLEM_CONSTANTS.L_Library{cycleIdx} = L;
            PROBLEM_CONSTANTS.U_Library{cycleIdx} = U;
            PROBLEM_CONSTANTS.P_Library{cycleIdx} = P;
        end
        L = PROBLEM_CONSTANTS.L_Library{cycleIdx};
        U = PROBLEM_CONSTANTS.U_Library{cycleIdx};
        P = PROBLEM_CONSTANTS.P_Library{cycleIdx};
        sol = U \ (L \ (P * indep));
    else
        Mat = buildSystemMatrix(PROBLEM_CONSTANTS, g_prefactor, dt, nr, cPoints, dr, SF);
        sol = Mat \ indep;
    end

    % --- Map Solution back to state ---
    % sol = [eta_rest; phi; v; z]
    eta_rest = sol(1:nr-cPoints);
    phi = sol(nr-cPoints+1 : 2*nr-cPoints);
    % p is not stored in next_condition but we could extract it if needed
    v_next = sol(end-1);
    z_next = sol(end);
    
    next_condition = previous_conditions;
    next_condition.bath_surface = [z_next * ones(cPoints, 1); eta_rest];
    next_condition.bath_potential = phi;
    next_condition.center_of_mass_velocity = v_next;
    next_condition.center_of_mass = z_next;
    next_condition.time = next_condition.time + dt;
    
    if any(isnan(sol)) || any(isinf(sol)) || abs(z_next) > 10
        error('Numerical Divergence Detected: CoM = %.2e. Terminating simulation.', z_next);
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

    % Variables: [eta_rest (nr-cPoints); phi (nr); v (1); z (1)]
    % Note: Pressure p is solved implicitly or eliminated. 
    % To keep it identical to original logic, we solve for:
    % [eta_rest; phi; p_contact; v; z] (Total: 2*nr + 2)
    
    N_rest = nr - cPoints;
    
    % Block 1: Kinematic condition (eta_rest)
    % eta_rest - dt * DTN_rest * phi = eta_prev
    % Size: N_rest x (2*nr + 2)
    Row1 = [eye(N_rest), -dt * DTN(cPoints+1:nr, :), zeros(N_rest, cPoints + 2)];
    
    % Block 2: Bernoulli (phi)
    % phi - dt*2*Delta/Re*phi + dt*(1/Fr*g_eff - 1/We*Delta)*eta = phi_prev
    % Note: eta is [z*ones; eta_rest]
    % Term: dt*(1/Fr*g_eff - 1/We*Delta) * [0; eta_rest]
    % Term: dt*(1/Fr*g_eff - 1/We*Delta) * [z*ones; 0]
    
    Op = dt * ( (eye(nr)/Fr * g_prefactor) - (Delta/We) );
    Visc = (eye(nr) - dt*2*Delta/Re);
    
    Row2_eta_rest = Op(:, cPoints+1:nr);
    Row2_phi = Visc;
    Row2_p = -dt * [eye(cPoints); zeros(N_rest, cPoints)]; % Pressure only under disk
    Row2_v = zeros(nr, 1);
    Row2_z = Op(:, 1:cPoints) * ones(cPoints, 1);
    
    Row2 = [Row2_eta_rest, Row2_phi, Row2_p, Row2_v, Row2_z];
    
    % Block 3: Newton 1 (v)
    % v + dt/Ma * Int(p) - SF*dt/dr * eta(cPoints+1) + SF*dt/dr * z = v_prev - dt/Fr*g_eff - force
    Row3_eta_rest = zeros(1, N_rest);
    Row3_eta_rest(1) = -SF*dt/dr;
    Row3_phi = zeros(1, nr);
    Row3_p = dt * pIntegral(1:cPoints) / Ma;
    Row3_v = 1;
    Row3_z = SF*dt/dr;
    Row3 = [Row3_eta_rest, Row3_phi, Row3_p, Row3_v, Row3_z];
    
    % Block 4: Newton 2 (z)
    % z - dt*v = z_prev
    Row4 = [zeros(1, N_rest + nr + cPoints), -dt, 1];
    
    % Block 5: Continuity at contact line (Implicit in original DTN?)
    % We need one more set of equations for p (cPoints)
    % Actually, in the original code, the pressure is a variable in the 2*nr system.
    
    % RE-CONSTRUCTING TO MATCH ORIGINAL EXACTLY
    % original Sist = [[eye(nr)-dt*2*Delta/Re,-dt*DTN];...
    %                  [dt*(eye(nr)/Fr * g_prefactor - Delta/We),eye(nr)-dt*2*Delta/Re]]; 
    
    Sist = [[eye(nr)-dt*2*Delta/Re,-dt*DTN];...
            [dt*(eye(nr)/Fr * g_prefactor - Delta/We),eye(nr)-dt*2*Delta/Re]]; 
      
    Mat =  [[Sist(:,(cPoints+1):2*nr),...
        [zeros(nr,cPoints);dt*eye(cPoints);zeros(nr-cPoints,cPoints)],...
        zeros(2*nr,1),Sist(:,1:cPoints)*ones(cPoints,1)];
        [-SF*dt/dr, zeros(1,2*nr-cPoints-1),-dt*pIntegral(1:cPoints)/Ma, 1 , SF*dt/dr];
        [zeros(1,2*nr-cPoints),-zeros(1, cPoints)  ,-dt,1]];
end
