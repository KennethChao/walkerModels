%%
close all;
clc;

%% Initial Condition (Guessed)
% k = 15;
% delta = 0.1392;
betaVec = 66:2:74

betaVec = 72;
for k = 1:length(betaVec)
beta = betaVec(k)/180*pi;
% beta = 72/180*pi;
g = 9.8;
parms.l = 1;
parms.m = 80;
%% Parameter Set


parms.g = g;
parms.beta = beta;
parms.k = k;
parms.delta = delta;

options = optimset('Display','off');

%% Sampling Number and Result Buffer
samplingNumbK = 25;
kMin = 0.1;
kMax = 25;
samplingNumbDelta = 3;
deltaMin = 1e-1;
deltaMax = 1.1;

stableSolution = nan*zeros(samplingNumbDelta,samplingNumbK);
unstableSolution = nan*zeros(samplingNumbDelta,samplingNumbK);

kVec = linspace(kMin,kMax,samplingNumbK);
deltaVec = linspace(deltaMin,deltaMax,samplingNumbDelta);

%%

for i = 1:samplingNumbK
    parms.k = kVec(i);
    for j = 1:samplingNumbDelta
    delta0 = deltaVec(j);
    parms.mode = 'fixedPointOpt';
    [x,fval,exitflag,output]  =fminunc(@(x)oneStepSimulation(x,parms),delta0,options);
        if  fval<1e-2 && x>1e-5
        parms.mode = 'perturbedSimulation';
        eigenValue =  oneStepSimulation(x,parms);
        
            if abs(eigenValue)<1
                stableSolution(j,i) = x;
            else
                unstableSolution(j,i) = x;
            end
        end
    end
        msg=sprintf('%d th iteration for k',i);
        disp(msg)
end


% plot(solution','bo')
for i = 1:3
%     plot(kVec,stableSolution(i,:),'bo')
%     plot(kVec,unstableSolution(i,:),'ro')
    
    plot(kVec,cos(stableSolution(i,:)),'bo')
    plot(kVec,cos(unstableSolution(i,:)),'ro')
    
    hold on
end

axis([0 kMax, 0 deltaMax])
end
%%
function result =  oneStepSimulation(delta0,parms)
mode = parms.mode;

beta = parms.beta;
if strcmp(mode,'fixedPointOpt')
    delta =delta0;
elseif strcmp(mode,'perturbedSimulation')
    perturbation = 1e-3;
    perturbedDelta = delta0+perturbation;
    delta = perturbedDelta;
end

normVel = 10;
dymStance = @(t,x) dymModelStance(t,x, parms);%dymModelStanceDimensionless

tspan = 0:0.01:100;
options = odeset('Event',@slipEventFcn);
x0 = [1,-normVel*cos(beta-delta),beta,normVel*sin(beta-delta)];

[t,x,te,xe,ie] = ode15s(dymStance,tspan,x0,options);     % Runge-Kutta 4th/5th order ODE solver
% te
% figure()
% plot(x(:,[1 3]))
% xVec = x(:,1).*cos(x(:,3))*-1;
% yVec = x(:,1).*sin(x(:,3));

% figure()
% plot(xVec,yVec);

y0 = x(end,1)*sin(x(end,3));
yd0 = x(end,2)*sin(x(end,3)) + x(end,1)*x(end,4)*cos(x(end,3));
% x0 = -x(end,1)*cos(x(end,3));
xd0 = -x(end,2)*cos(x(end,3)) + x(end,1)*x(end,4)*sin(x(end,3));
x0 = [y0,yd0];
options = odeset('Event',@(t,x)tdEventFcn(t,x,beta));
dymFlight = @(t,a) dymModelFlightDimensionless(t,a, parms);
[t,x2,te2,ae2,ie2] = ode45(dymFlight,tspan,x0,options);     % Runge-Kutta 4th/5th order ODE solver
% te2
% figure()
% plot(x2(:,[1]))
velVec = [xd0,-x2(end,2)];
deltaNew = atan2(velVec(2),velVec(1));

diffDelta = deltaNew -delta;
diffBeta = normVel-  norm(velVec);


if strcmp(mode,'fixedPointOpt')
    result = abs(diffDelta) + abs(diffBeta);
elseif strcmp(mode,'perturbedSimulation')
    result = diffDelta/perturbation;
end


end









