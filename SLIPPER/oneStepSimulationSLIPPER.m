function result = oneStepSimulationSLIPPER(xVec, parms)
%ONESTEPSIMULATIONSLIPPER run a one step simulation of SLIPPER
%   With an initial condition 'xVec' and parameter set 'parms', run a one
%   step simulation using ode45.
%
%   Depending on the parms.mode, the result returns different objects:
%   parms.mode = 'fixedPointOpt': return the optimization cost
%   parms.mode = 'simulationCheck': return the struct of state/time data
%   parms.mode = 'perturbedSimulation': return the Poincare map
%
%   The model and the initial conditions are extended from the SLIP simulation.
%   For more information, please check oneStepSimulationSLIP in SLIP folder.

%% Extract the initial conditions and parameters
% initial conditions
phiInitialGuess = xVec(1);
phidInitialGuess = xVec(2);
uInitiailGuess = xVec(3);
deltaGirdPoint = parms.delta0;

% parameters
beta = parms.beta;
rc = parms.rc;
mf = parms.mf;
optWeighting = parms.optWeighting;

delta = deltaGirdPoint;
phi0 = phiInitialGuess;
phid0 = phidInitialGuess;
mode = parms.mode;

%% Solver setup

% ode45 solver option
optionsStance = odeset('Event', @eventFcnLiftOffSLIPPER, 'AbsTol', 1e-9, 'RelTol', 1.e-8);
dymStance = @(t, x) dymModelStanceSLIPPER(t, x, uInitiailGuess, parms);
optionsFlight = odeset('Event', @(t, x)eventFncTouchDownSLIPPER(t, x, parms), 'AbsTol', 1e-9, 'RelTol', 1.e-8);
dymFlight = @(t, a) dymModelFlightSLIPPER(t, a, parms);

% assign iteration number of simulation
if strcmp(mode, 'fixedPointOpt')
    iterationNumber = 1;
elseif strcmp(mode, 'simulationCheck')
    iterationNumber = 1;
elseif strcmp(mode, 'perturbedSimulation')
    perturbation = 5e-3;
    iterationNumber = 6;

    perturbedDeltaPlus = deltaGirdPoint + perturbation;
    perturbedDeltaMinus = deltaGirdPoint - perturbation;
    
    perturbedPhiPlus = phiInitialGuess + perturbation;
    perturbedPhiMinus = phiInitialGuess - perturbation;
    
    perturbedPhidPlus = phidInitialGuess + perturbation;
    perturbedPhidMinus = phidInitialGuess - perturbation;    
end

tspan = 0:0.01:20;

for i = 1:iterationNumber
    
    %% Update initial conditions
    normVel = 1;

    if strcmp(mode, 'perturbedSimulation') && i == 1
        delta = perturbedDeltaPlus;
        phi0 = phiInitialGuess;
        phid0 = phidInitialGuess;
    elseif strcmp(mode, 'perturbedSimulation') && i == 2
        delta = perturbedDeltaMinus;
        phi0 = phiInitialGuess;
        phid0 = phidInitialGuess;        
    elseif strcmp(mode, 'perturbedSimulation') && i == 3
        delta = deltaGirdPoint;
        phi0 = perturbedPhiPlus;
        phid0 = phidInitialGuess;                
    elseif strcmp(mode, 'perturbedSimulation') && i == 4
        delta = deltaGirdPoint;
        phi0 = perturbedPhiMinus;
        phid0 = phidInitialGuess;  
    elseif strcmp(mode, 'perturbedSimulation') && i == 5
        delta = deltaGirdPoint;
        phi0 = phiInitialGuess;
        phid0 = perturbedPhidPlus;         
    elseif strcmp(mode, 'perturbedSimulation') && i == 6
        delta = deltaGirdPoint;
        phi0 = phiInitialGuess;
        phid0 = perturbedPhidMinus;          
    end
    
    %%  Simulation in stance phase
    
    % Prep initial conditions
    xbLO = [1, -normVel * cos(beta-delta), beta, normVel * sin(beta-delta), phi0, phid0];
    
    % Solve ode 
    [t, x, te, xe, ~] = ode45(dymStance, tspan, xbLO, optionsStance); % Runge-Kutta 4th/5th order ODE solver
    
    %%  Simulation in flight phase
    
    % Convert states to Cartesian space
    % position and velocity of mf
    if isempty(te)
        xe = x(end, :);
    end
    % state at liftoff (LO)
    lLO = xe(1);
    ldLO = xe(2);
    thetaLO = xe(3);
    thetadLO = xe(4);
    phiLO = xe(5);
    phidLO = xe(6);
    
    % Motions in Cartesian space of frame (f), body (b) and COM (c)
    % position and velocity of frame (hub) at LO
    mfMotion = frameMassCartesianMotionStance(lLO, thetaLO, phiLO, ldLO, thetadLO, phidLO, rc);
    zfLO = mfMotion(2);
    zfdLO = mfMotion(4);
    xfLO = mfMotion(1);
    xfdLO = mfMotion(3);
    
    % position and velocity of body at LO
    mbMotion = bodyMassCartesianMotionStance(lLO, thetaLO, phiLO, ldLO, thetadLO, phidLO, rc);    
    zbLO = mbMotion(2);
    zbdLO = mbMotion(4);
    xbLO = mbMotion(1);
    xbdLO = mbMotion(3);
    
    % position and velocity of COM at LO
    comMotion = motionPointMass2COM(xfLO, zfLO, xbLO, zbLO, xfdLO, zfdLO, xbdLO, zbdLO, rc, mf);
    xcLO = comMotion(1);
    xcdLO = comMotion(3);
    zcLO = comMotion(2);
    zcdLO = comMotion(4);
    
    % Prep initial conditions    
    x20 = [zcLO, zcdLO, phiLO, phidLO];
    
    if zbdLO < 0
        x2 = ones(1, 4) * 1e1;
        t2 = 20;
    else
        % Solve ode
        [t2, x2, te2, xe2, ~] = ode45(dymFlight, tspan, x20, optionsFlight); % Runge-Kutta 4th/5th order ODE solver
    end
    
    % Get the position and velocity of mf from COM states and phi/phid
    
    % state at touchDown (TD)
    zcTD = x2(end, 1);
    zcdTD = x2(end, 2);
    
    xcdTD = xcdLO;
    xcTD = xcLO + t2(end) * xcdLO;
    
    phiTD = x2(end, 3);
    phidTD = x2(end, 4);
    
    % position and velocity of frame (hub) at TD
    pointMassMotion = motionCOM2PointMass(xcTD, zcTD, phiTD, xcdTD, zcdTD, phidTD, rc, mf);
    xfdTD = pointMassMotion(5);
    zfdTD = pointMassMotion(6);
    
    %%
    % Get states of the next step
    velVec = [xfdTD, -zfdTD];
    deltaNew = atan2(velVec(2), velVec(1));
    
    if strcmp(mode, 'perturbedSimulation') && i == 1
        deltaNewPlus = [deltaNew; x2(end, 3); x2(end, 4)];
    elseif strcmp(mode, 'perturbedSimulation') && i == 2
        deltaNewMinus = [deltaNew; x2(end, 3); x2(end, 4)];
    elseif strcmp(mode, 'perturbedSimulation') && i == 3
        phiNewPlus = [deltaNew; x2(end, 3); x2(end, 4)];
    elseif strcmp(mode, 'perturbedSimulation') && i == 4
        phiNewMinus = [deltaNew; x2(end, 3); x2(end, 4)];
    elseif strcmp(mode, 'perturbedSimulation') && i == 5
        phidNewPlus = [deltaNew; x2(end, 3); x2(end, 4)];
    elseif strcmp(mode, 'perturbedSimulation') && i == 6
        phidNewMinus = [deltaNew; x2(end, 3); x2(end, 4)];
    end
    
end

if strcmp(mode, 'fixedPointOpt')
    % Return cost
    diffDelta = deltaNew - delta;
    diffVel = 1 - norm(velVec);
    diffPhi = norm([phi0; phid0]-[phiTD; phidTD]);
    if strcmp(parms.controlMode, 'pControl') || strcmp(parms.controlMode, 'constantTorque')
        result = optWeighting(1) * diffDelta^2 + optWeighting(2) * diffVel^2 + optWeighting(3) * diffPhi^2 + optWeighting(4) * uInitiailGuess^2;
    else
        result = optWeighting(1) * diffDelta^2 + optWeighting(2) * diffVel^2 + optWeighting(3) * diffPhi^2;
    end
elseif strcmp(mode, 'perturbedSimulation')
    % Return Poincare Map
    result = [(deltaNewPlus - deltaNewMinus) / 2 / perturbation, ...
                  (phiNewPlus - phiNewMinus) / 2 / perturbation, ...
                  (phidNewPlus - phidNewMinus) / 2 / perturbation
        ];
elseif strcmp(mode, 'simulationCheck')
    % Return the struct of state and time data
    result.x = x;
    result.t = t;
    result.te = te;
    result.xe = xe;
    result.x2 = x2;
    result.xe2 = xe2;
    result.t2 = t2;
    result.te2 = te2;
end

end