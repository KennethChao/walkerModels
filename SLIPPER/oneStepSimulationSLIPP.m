function result = oneStepSimulationSLIPP(xVec, parms)
%ToDo change perturbed simulation
% clc;
% clear;
% close all;

%% For reusing the liftoff event function from SLIP model
% addpath('../SLIP')

%%

% delta0 = xVec(1);
phi0 = xVec(1);
phid0 = xVec(2);
delta0 = parms.delta0;

%%
beta = parms.beta;
rc = parms.rc;
mf = parms.mf;

delta = delta0;
mode = parms.mode;
%% Solver setup
    % ode45 solver option
    optionsStance = odeset('Event', @liftOffEventFcn);%, 'RelTol',1.e-6
    dymStance = @(t, x) dymModelStanceSLIPPendulum(t, x, parms); %dymModelStanceDimensionless
    optionsFlight = odeset('Event', @(t, x)touchDownSLIPPendulumEventFcn(t, x, parms)); %,'RelTol',1.e-6
    dymFlight = @(t, a) dymModelFlightSLIPPendulum(t, a, parms);
    
if strcmp(mode, 'fixedPointOpt')
    iterNumb = 1;
elseif strcmp(mode, 'simulationCheck')    
    iterNumb = 1;
elseif strcmp(mode, 'perturbedSimulation')
    perturbation = 1e-2;
    iterNumb = 4;
end

tspan = 0:0.01:20;

for i = 1:iterNumb
    %% Update initial conditions
    normVel = 1;
    
    if i > 1
%         normVel = norm(velVec);
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
    x0 = [1, -normVel * cos(beta-delta), beta, normVel * sin(beta-delta), phi0, phid0];

    % Solve ode
    [t, x, te, xe, ~] = ode45(dymStance, tspan, x0, optionsStance); % Runge-Kutta 4th/5th order ODE solver
    
    %%  Simulation in flight phase
%     plot(x(:,5))
    % Convert states to Cartesian space
    % position and velocity of mf
    
    
    if  isempty(te)
            xe = x(end,:);
    end
    zf0 = xe(1) * sin(xe(3));
    zfd0 = xe(2) * sin(xe(3)) + xe(1) * xe(4) * cos(xe(3));
    xf0 = -xe(1) * cos(xe(3)); 
    xfd0 = -xe(2) * cos(xe(3)) + xe(1) * xe(4) * sin(xe(3));
    
    % position and velocity of m
    z0 = zf0 - rc*sin(xe(5));
    zd0 = zfd0 - rc*cos(xe(5))*xe(6);
    x0 = xf0 + rc*cos(xe(5)); 
    xd0 = xfd0 - rc*sin(xe(5))*xe(6);
    
    % position and velocity of COM
    xc0 = 1/(1+mf)*x0 + mf/(1+mf)*xf0;
    xcd0 = 1/(1+mf)*xd0 + mf/(1+mf)*xfd0;
    zc0 = 1/(1+mf)*z0 + mf/(1+mf)*zf0;
    zcd0 = 1/(1+mf)*zd0 + mf/(1+mf)*zfd0;
    
    
%     % angular momentum of body around COM
%     rcb = [x0-xc0,0,z0-zc0];
%     vb = [xd0,0, zd0];
%     rcdCrossVb = cross(rcb,vb);
%     Lbody = rcdCrossVb;
%     
%     % angular momentum of frame around COM
%     rcf = [xf0-xc0,0,zf0-zc0];
%     vf = [xfd0,0, zfd0];
%     rcfCrossVf = cross(rcf,vf);
%     Lframe = rcfCrossVf*mf;    
%     
%     % angular velocity of COM 
%     I = mf/(1+mf)*rc^2;
%     omega = (Lframe+Lbody) / I;
%     phid0 = omega(2);
%     
%     phi0 = xe(5);
    % Initial condition of flight phase
    x20 = [zc0,zcd0,xe(5),xe(6)];
    
    if zd0<0 

    x2 = ones(1,4)*1e3;
    else
% Solve ode
    [t2, x2, te2, xe2, ~] = ode45(dymFlight, tspan, x20, optionsFlight); % Runge-Kutta 4th/5th order ODE solver
    end

    
    %%
%     % angular momentum of body around COM
%     rcb = [x0-xc0,0,z0-zc0];
%     vb = [xd0,0, zd0];
%     rcdCrossVb = cross(rcb,vb);
%     Lbody = rcdCrossVb;
%     
%     % angular momentum of frame around COM
%     rcf = [xf0-xc0,0,zf0-zc0];
%     vf = [xfd0,0, zfd0];
%     rcfCrossVf = cross(rcf,vf);
%     Lframe = rcfCrossVf*mf;    
%     
%     % angular velocity of COM 
%     I = mf/(1+mf)*rc^2;
%     omega = (Lframe+Lbody) / I;
%     phid0 = omega(2);
%     
%     phi0 = xe(5); 
%     xcEnd = xc0 + xcd0*te2;
    xdcEnd = xcd0;
%     zcEnd = x2(end,1);
    zdcEnd = x2(end,2);

    rc2f = 1/(1+mf);    
    rVec = [-rc*rc2f*cos(x2(end,3)), 0 ,rc*rc2f*sin(x2(end,3))];
    omega = [0, x2(end,4),0];
    
    vRotation = cross(omega,rVec);
    
    xdfendFlight =  xdcEnd + vRotation(1); 
    zdfendFlight = zdcEnd + vRotation(3);    
    
%     zfendFlight = x2(end,1) + rc*rc2f*sin(x2(end,3));
    try
%     zdfendFlight = x2(end,2) + rc*rc2f*cos(x2(end,3))*x2(end,4);
    catch
        te2
    end
%     xfendFlight = xc20 + xdc20*t2(end) - rc*rc2f*cos(x2(end,3)); 
%     xdfendFlight =  xdc0 + rc*rc2f*sin(x2(end,3))*x2(end,4); 
    
    
    %%
    
    % Get states of the next step
    velVec = [xdfendFlight, -zdfendFlight];
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
    diffPhi = norm([phi0;phid0]-[x2(end,3);x2(end,4)]);
    result = diffDelta^2 + diffVel^2 + diffPhi^2;
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