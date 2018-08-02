% function result = oneStepSimulationSLIPP(delta0, parms)
%ToDo change perturbed simulation
clc;
clear;
close all;

parms.g = 0.12;
parms.k = 16;

beta = 72/180*pi;%parms.beta;
delta = 0.1;
mf = 0.2;

parms.mf = mf;
parms.beta = beta;

mode = 'fixedPointOpt';

rc = 0.2;
parms.rc = rc;

% Ode solver setup
    optionsStance = odeset('Event', @liftOffEventFcn);%,'RelTol',1.e-6
    dymStance = @(t, x) dymModelStanceSLIPPendulum(t, x, parms); %dymModelStanceDimensionless
    optionsFlight = odeset('Event', @(t, x)touchDownSLIPPendulumEventFcn(t, x, parms),'RelTol',1.e-6);
    dymFlight = @(t, a) dymModelFlightSLIPPendulum(t, a, parms);


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
    x0 = [1, -normVel * cos(beta-delta), beta, normVel * sin(beta-delta), pi/2, 0];

    % Solve ode
    [t, x, te, xe, ie] = ode45(dymStance, tspan, x0, optionsStance); % Runge-Kutta 4th/5th order ODE solver
    
    %%  Simulation in flight phase
    plot(x(:,5))
    % Convert states to Cartesian space
    % position and velocity of mf
    zf0 = xe(1) * sin(xe(3));
    zfd0 = xe(2) * sin(xe(3)) + xe(1) * xe(4) * cos(xe(3));
    xf0 = -xe(1) * cos(xe(3)); 
    xfd0 = -xe(2) * cos(xe(3)) + xe(1) * xe(4) * sin(xe(3));
    
    % position and velocity of m
    z0 = xe(1) * sin(xe(3)) - rc*sin(xe(5));
    zd0 = xe(2) * sin(xe(3)) + xe(1) * xe(4) * cos(xe(3)) - rc*cos(xe(5))*xe(6);
    x0 = -xe(1) * cos(xe(3)) + rc*cos(xe(5)); 
    xd0 = -xe(2) * cos(xe(3)) + xe(1) * xe(4) * sin(xe(3)) - rc*sin(xe(5))*xe(6);
    
    % position and velocity of COM
    xc0 = 1/(1+mf)*x0 + mf/(1+mf)*xf0;
    xcd0 = 1/(1+mf)*xd0 + mf/(1+mf)*xfd0;
    
    zc0 = 1/(1+mf)*z0 + mf/(1+mf)*zf0;
    zcd0 = 1/(1+mf)*zd0 + mf/(1+mf)*zfd0;
    
    % angular momentum of body around COM
    rcb = [x0-xc0,z0-zc0];
    vb = [xd0, zd0];
    rcdCrossVb = rcb(1)*vb(2)-rcb(2)*vb(1);
    Lbody = -rcdCrossVb*1;
    
    % angular momentum of frame around COM
    rcf = [xf0-xc0,zf0-zc0];
    vf = [xfd0, zfd0];
    rcfCrossVf = rcf(1)*vf(2)-rcf(2)*vf(1);
    Lframe = -rcfCrossVf*mf;    
    
    % angular velocity of COM 
    I = mf/(1+mf);
    phid0 = Lframe+Lbody / I;
    phi0 = xe(5);
    % Initial condition of flight phase
    x20 = [zc0,zcd0,phi0,phid0];
    % Solve ode
    if zd0<0
        x2 = x0*1e3;
    else
        [t2, x2, te2, ae2, ie2] = ode45(dymFlight, tspan, x20, optionsFlight); % Runge-Kutta 4th/5th order ODE solver
    end
    
    %%
    xc = xf0*mf/(mf+1) + x0/(mf+1);
    xdc = xfd0*mf/(mf+1) + xd0/(mf+1);
    
    rc2f = 1/(1+mf);
    rc2b = mf/(1+mf);
    
    zfendFlight = x2(end,1) + rc*rc2f*sin(x2(end,3));
    zdfendFlight = x2(end,2) + rc*rc2f*cos(x2(end,3))*x2(end,4);
    xfendFlight = xc + xdc*t2(end) - rc*rc2f*cos(x2(end,3)); 
    xdfendFlight =  xdc + rc*rc2f*sin(x2(end,3))*x2(end,4); 
    
    % position and velocity of m
    zendFlight = x2(end,1) - rc*rc2b*sin(x2(end,3));
    zdendFlight = x2(end,2) - rc*rc2b*cos(x2(end,3))*x2(end,4);
    xendFlight = xc + xdc*t2(end) + rc*rc2b*cos(x2(end,3)); 
    xdendFlight = xdc - rc*rc2b*sin(x2(end,3))*x2(end,4); 
    
    posd = [xendFlight; 0 ; zendFlight];
    posfd = [xfendFlight; 0 ; zfendFlight];
    
    
    vd = [xdendFlight; 0 ; zdendFlight];
    vfd = [xdfendFlight; 0 ; zdfendFlight];
    
    (vd-vfd)./(posd-posfd)
    
    
    %%
    
    te2
    figure()
    plot(x2(:,1));
    % Get states of the next step
    velVec = [xd0, -x2(end, 2)];
    deltaNew = atan2(velVec(2), velVec(1));
    
    if strcmp(mode, 'perturbedSimulation') && i == iterNumb
        deltaNewPlus = deltaNew;
    elseif strcmp(mode, 'perturbedSimulation') && i == iterNumb - 1
        deltaNewMinus = deltaNew;
    end
    
end

% if strcmp(mode, 'fixedPointOpt')
%     % Return cost
%     diffDelta = deltaNew - delta;
%     diffBeta = 1 - norm(velVec);
%     result = (diffDelta)^2 + (diffBeta)^2;
% elseif strcmp(mode, 'perturbedSimulation')
%     % Return second eigen value
%     result = (deltaNewPlus - deltaNewMinus) / 2 / perturbation;
% end

% end