% function result = oneStepSimulationSLIP_walking(delta0, parms)
%ToDo change perturbed simulation
clc;clear;
close all;

mode = 'fixedPointOpt';

parms.g = 1.2;
parms.k = 10;

stepLength = 0.6;
beta = 60/180*pi;%parms.beta;
delta = -0.8;

% Ode solver setup
    optionsStance = odeset('Event', @(t, x)touchDownWalkingEventFcn(t, x,beta),'RelTol',1.e-6);
    dymStance = @(t, x) dymModelStanceDimensionless(t, x, parms); %dymModelStanceDimensionless
    optionsDoubleStance = odeset('Event', @liftOffDoubleSupportEventFcn,'RelTol',1.e-6);
if strcmp(mode, 'fixedPointOpt')
    iterNumb = 1;
elseif strcmp(mode, 'perturbedSimulation')
    perturbation = 1e-3;
    iterNumb = 2;
end

tspan = 0:0.01:20;

for i = 1:iterNumb
    %% Update initial conditions
    normVel = 1;
    
    if i > 1
        delta = deltaNew;
    end
    if strcmp(mode, 'perturbedSimulation') && i == iterNumb
        delta = perturbedDeltaPlus;
    elseif strcmp(mode, 'perturbedSimulation') && i == iterNumb - 1
        newDelta0 = delta0;%deltaNew;
        perturbedDeltaPlus = newDelta0 + perturbation;
        perturbedDeltaMinus = newDelta0 - perturbation;
        delta = perturbedDeltaMinus;
    end
    
    %%  Simulation in stance phase
    
    % Prep initial conditions
    xInit = stepLength - cos(beta);
    yInit = sin(beta);
    
    l0 = sqrt(xInit^2+yInit^2);
    alpha = atan2(yInit,xInit);
    
    
    x0 = [l0, -normVel * cos(alpha-delta), alpha, normVel *l0* sin(alpha-delta)];

    % Solve ode
    [t, x1, te, xe, ie] = ode45(dymStance, tspan, x0', optionsStance); %#ok<ASGLU> % Runge-Kutta 4th/5th order ODE solver
    
    %%  Simulation in flight phase
    plot(-x1(:,1).*cos(x1(:,3)),x1(:,1).*sin(x1(:,3)))
    % Convert states to Cartesian space
    y0 = x1(end, 1) * sin(x1(end, 3));
    yd0 = x1(end, 2) * sin(x1(end, 3)) + x1(end, 1) * x1(end, 4) * cos(x1(end, 3));
    x0 = -x1(end, 1) * cos(x1(end, 3)); 
    xd0 = -x1(end, 2) * cos(x1(end, 3)) + x1(end, 1) * x1(end, 4) * sin(x1(end, 3));
    
    x20 = [x0,xd0, y0, yd0]';
    
    parms.trailingFootPos = [0,0];
    parms.leadingFootPos = [x0 + cos(beta),0];
    % Solve ode
%     if yd0<0
%         x2 = x0*1e3;
%     else

    dymDoubleStance = @(t, x)dymModelDoubleSupportDimensionless(t, x, parms);

    [t, x2, te2, ae2, ie2] = ode45(dymDoubleStance, tspan, x20, optionsDoubleStance); % Runge-Kutta 4th/5th order ODE solver
%     end
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

% end