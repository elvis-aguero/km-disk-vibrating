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
    g_prefactor = (1 - gamma * cos(bath_w * (t + dt) + phase)); % CHANGED: Using bath-specific frequency and phase
    indep = [b; CoM_vel - dt/Fr * g_prefactor - dt*F*cos(w*(t+dt)); CoM]; % CHANGED: F and w are for the disk forcing
    
    % If the matrix exists and gravity is constant, import it. 
    % Otherwise (oscillating gravity), we must recompute it every step.
    if ~any(isnan(PROBLEM_CONSTANTS.precomputedInverse(:))) && gamma == 0 % CHANGED: Using any() to ensure scalar logical for &&
        invMat = PROBLEM_CONSTANTS.precomputedInverse;
        sol = invMat * indep;
    else
    % If it does not, create it and then save it
        
        We = PROBLEM_CONSTANTS.weber;
        Re = PROBLEM_CONSTANTS.reynolds;
        Delta = PROBLEM_CONSTANTS.laplacian;
        DTN = PROBLEM_CONSTANTS.DTN;
        
        
        pIntegral = PROBLEM_CONSTANTS.pressure_integral;
        Ma = PROBLEM_CONSTANTS.obj_mass;
        
        % Preparing the matrix (2x2 block)
        Sist = [[eye(nr)-dt*2*Delta/Re,-dt*DTN];...
            [dt*(eye(nr)/Fr * g_prefactor - Delta/We),eye(nr)-dt*2*Delta/Re]]; % CHANGED: Added gravity prefactor to fluid equation
          
        % Completing the system
        Mat =  [[Sist(:,(cPoints+1):2*nr),...
            [zeros(nr,cPoints);dt*eye(cPoints);zeros(nr-cPoints,cPoints)],...
            zeros(2*nr,1),Sist(:,1:cPoints)*ones(cPoints,1)];
            [-SF*dt/dr, zeros(1,2*nr-cPoints-1),-dt*pIntegral(1:cPoints)/Ma, 1 , SF*dt/dr];
            [zeros(1,2*nr-cPoints),-zeros(1, cPoints)  ,-dt,1]];
        % Now save the matrix for later
        PROBLEM_CONSTANTS.precomputedInverse = inv(Mat);

        sol = Mat\indep;
    end


    next_condition = previous_conditions;
    next_condition.bath_surface = [sol(end)* ones(cPoints, 1); ...
                                                sol(1:nr-cPoints)];
    next_condition.bath_potential = sol(nr-cPoints+1:2*nr-cPoints);
    next_condition.pressure = sol(2*nr-cPoints+1:2*nr);
    next_condition.center_of_mass_velocity = sol(end-1);
    next_condition.center_of_mass = sol(end);
    next_condition.time = next_condition.time + dt;
    
end