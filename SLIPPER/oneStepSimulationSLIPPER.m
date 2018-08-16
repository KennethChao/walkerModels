function result = oneStepSimulationSLIPPER(xVec, parms)
% ToDo

%% For reusing the liftoff event function from SLIP model
%%
phi0 = xVec(1);
phid0 = xVec(2);
delta0 = parms.delta0;

%%
beta = parms.beta;
rc = parms.rc;
mf = parms.mf;
optWeighting = parms.optWeighting;

delta = delta0;
mode = parms.mode;
%% Solver setup
    % ode45 solver option
    optionsStance = odeset('Event', @liftOffEventFcn,'AbsTol',1e-12,'RelTol',1.e-8);%, 
    dymStance = @(t, x) dymModelStanceSLIPPendulum(t, x, parms); %dymModelStanceDimensionless
    optionsFlight = odeset('Event', @(t, x)touchDownSLIPPendulumEventFcn(t, x, parms),'AbsTol',1e-12,'RelTol',1.e-8); %,
    dymFlight = @(t, a) dymModelFlightSLIPPendulum(t, a, parms);
    
if strcmp(mode, 'fixedPointOpt')
    iterNumb = 1;
elseif strcmp(mode, 'simulationCheck')    
    iterNumb = 1;
elseif strcmp(mode, 'perturbedSimulation')
    perturbation = 5e-3;
    iterNumb = 4;
end

tspan = 0:0.01:20;

for i = 1:iterNumb
    %% Update initial conditions
    normVel = 1;
    
    if i > 1
        delta = deltaNew;
        phi0 =  x2(end,3);
        phid0 =  x2(end,4);
    end
    if strcmp(mode, 'perturbedSimulation') && i == iterNumb
        delta = perturbedDeltaPlus;
    elseif strcmp(mode, 'perturbedSimulation') && i == iterNumb - 1
        delta = perturbedDeltaMinus;
    elseif strcmp(mode, 'perturbedSimulation') && i == iterNumb - 2        
        phi0 = perturbedPhiPlus;        
    elseif strcmp(mode, 'perturbedSimulation') && i == iterNumb - 3
        newDelta0 = delta0;%deltaNew;
        perturbedDeltaPlus = newDelta0 + perturbation;
        perturbedDeltaMinus = newDelta0 - perturbation;
        newPhi0 = phi0;%deltaNew;
        perturbedPhiPlus = newPhi0 + perturbation;
        perturbedPhiMinus = newPhi0 - perturbation;        
        %
        phi0 = perturbedPhiMinus;
    end
    
    %%  Simulation in stance phase
    
    % Prep initial conditions
    xbLO = [1, -normVel * cos(beta-delta), beta, normVel * sin(beta-delta), phi0, phid0];

    % Solve ode
    [t, x, te, xe, ~] = ode45(dymStance, tspan, xbLO, optionsStance); % Runge-Kutta 4th/5th order ODE solver
    
    %%  Simulation in flight phase
    % Convert states to Cartesian space
    % position and velocity of mf
    if  isempty(te)
            xe = x(end,:);
    end
    lLO = xe(1);
    ldLO = xe(2);
    thetaLO = xe(3);
    thetadLO = xe(4);
    phiLO = xe(5);
    phidLO = xe(6);
    
    mfMotion = frameMassCartesianMotionStance(lLO,thetaLO,phiLO,ldLO,thetadLO,phidLO,rc);
    mbMotion = bodyMassCartesianMotionStance(lLO,thetaLO,phiLO,ldLO,thetadLO,phidLO,rc);

    zfLO = mfMotion(2);
    zfdLO = mfMotion(4);
    xfLO = mfMotion(1); 
    xfdLO = mfMotion(3);
    
    % position and velocity of m
    zbLO = mbMotion(2);
    zbdLO = mbMotion(4);
    xbLO = mbMotion(1); 
    xbdLO = mbMotion(3);
    
    % position and velocity of COM
    comMotion = motionPointMass2COM(xfLO,zfLO,xbLO,zbLO,xfdLO,zfdLO,xbdLO,zbdLO,rc,mf);

    xcLO = comMotion(1);
    xcdLO = comMotion(3);
    zcLO = comMotion(2);
    zcdLO = comMotion(4);

    
    x20 = [zcLO,zcdLO,phiLO,phidLO];
    
    if zbdLO<0 
    x2 = ones(1,4)*1e1;
    t2 = 20;
    else
% Solve ode
    [t2, x2, te2, xe2, ~] = ode45(dymFlight, tspan, x20, optionsFlight); % Runge-Kutta 4th/5th order ODE solver
    end

    
    %%  
    zcTD = x2(end,1);
    zcdTD = x2(end,2);
    
    xcdTD = xcdLO;
    xcTD = xcLO + t2(end)*xcdLO;
    
    phiTD = x2(end,3);
    phidTD = x2(end,4);
    
    pointMassMotion = motionCOM2PointMass(xcTD,zcTD,phiTD,xcdTD,zcdTD,phidTD,rc,mf);

    xfdTD = pointMassMotion(5);
    zfdTD = pointMassMotion(6);
    %%
    
    % Get states of the next step
    velVec = [xfdTD, -zfdTD];
    deltaNew = atan2(velVec(2), velVec(1));
    
    if strcmp(mode, 'perturbedSimulation') && i == iterNumb
        deltaNewPlus = [deltaNew; x2(end,3)];
    elseif strcmp(mode, 'perturbedSimulation') && i == iterNumb - 1
        deltaNewMinus = [deltaNew; x2(end,3)];
    elseif strcmp(mode, 'perturbedSimulation') && i == iterNumb - 2
        phiNewPlus = [deltaNew; x2(end,3)];
    elseif strcmp(mode, 'perturbedSimulation') && i == iterNumb - 3
        phiNewMinus = [deltaNew; x2(end,3)];
    end
    
end
% 
if strcmp(mode, 'fixedPointOpt')
    % Return cost
    diffDelta = deltaNew - delta;
    diffVel = 1 - norm(velVec);
    diffPhi = norm([phi0;phid0]-[phiTD;phidTD]);
    result = optWeighting(1)*diffDelta^2 + optWeighting(2)*diffVel^2 + optWeighting(3)*diffPhi^2;
elseif strcmp(mode, 'perturbedSimulation')
    % Return second eigen value
    result = [(deltaNewPlus - deltaNewMinus) / 2 / perturbation,...
              (phiNewPlus - phiNewMinus) / 2 / perturbation;
             ];
elseif strcmp(mode, 'simulationCheck')
    % Return second eigen value
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