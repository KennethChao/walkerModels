function result = oneStepSimulationSLIP(delta0, parms)
%ToDo change perturbed simulation

mode = parms.mode;

beta = parms.beta;
delta = delta0;

% Ode solver setup
    optionsStance = odeset('Event', @liftOffEventFcn);
    dymStance = @(t, x) dymModelStanceDimensionless(t, x, parms); %dymModelStanceDimensionless
    optionsFlight = odeset('Event', @(t, x)touchDownEventFcn(t, x, beta));
    dymFlight = @(t, a) dymModelFlightDimensionless(t, a, parms);


if strcmp(mode, 'fixedPointOpt')
    iterNumb = 1;
elseif strcmp(mode, 'perturbedSimulation')
    perturbation = 1e-3;
    iterNumb = 3;
end

tspan = 0:0.01:5;

for i = 1:iterNumb
    %% Update initial conditions
    normVel = 1;
    
    if i > 1
        delta = deltaNew;
    end
    if strcmp(mode, 'perturbedSimulation') && i == iterNumb
        delta = perturbedDeltaPlus;
    elseif strcmp(mode, 'perturbedSimulation') && i == iterNumb - 1
        newDelta0 = deltaNew;
        perturbedDeltaPlus = newDelta0 + perturbation;
        perturbedDeltaMinus = newDelta0 - perturbation;
        delta = perturbedDeltaMinus;
    end
    
    %%  Simulation in stance phase
    
    % Prep initial conditions
    x0 = [1, -normVel * cos(beta-delta), beta, normVel * sin(beta-delta)];

    % Solve ode
    [t, x, te, xe, ie] = ode45(dymStance, tspan, x0, optionsStance); %#ok<ASGLU> % Runge-Kutta 4th/5th order ODE solver
    
    %%  Simulation in flight phase
    
    % Convert states to Cartesian space
    y0 = x(end, 1) * sin(x(end, 3));
    yd0 = x(end, 2) * sin(x(end, 3)) + x(end, 1) * x(end, 4) * cos(x(end, 3));
    x0 = -x(end, 1) * cos(x(end, 3)); %#ok<NASGU>
    xd0 = -x(end, 2) * cos(x(end, 3)) + x(end, 1) * x(end, 4) * sin(x(end, 3));
    
    x0 = [y0, yd0];
    % Solve ode
    if yd0<0
        x2 = x0*1e3;
    else
        [t, x2, te2, ae2, ie2] = ode45(dymFlight, tspan, x0, optionsFlight); %#ok<ASGLU> % Runge-Kutta 4th/5th order ODE solver
    end
    % Get states of the next step
    velVec = [xd0, -x2(end, 2)];
    deltaNew = atan2(velVec(2), velVec(1));
    
    if strcmp(mode, 'perturbedSimulation') && i == iterNumb
        deltaNewPlus = deltaNew;
    elseif strcmp(mode, 'perturbedSimulation') && i == iterNumb - 1
        deltaNewMinus = deltaNew;
    end
    
end

if strcmp(mode, 'fixedPointOpt')
    % Return cost
    diffDelta = deltaNew - delta;
    diffBeta = 1 - norm(velVec);
    result = (diffDelta)^2 + (diffBeta)^2;
elseif strcmp(mode, 'perturbedSimulation')
    % Return second eigen value
    result = (deltaNewPlus - deltaNewMinus) / 2 / perturbation;
end

end