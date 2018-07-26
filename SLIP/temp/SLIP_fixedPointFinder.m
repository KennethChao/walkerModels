%%
close all;
clc;
clear;

%% Initial Condition (Guessed)
% k = 15;
% delta = 0.1392;
% betaVec = 60:2:80
% betaVec = 60:10:75
% betaVec = linspace(66,74,5)
betaVec = 72;
for k = 1:length(betaVec)
beta = betaVec(k)/180*pi;
% beta = 72/180*pi;
g = 0.46;
parms.l = 1;
parms.m = 80;
parms.legNumber = 4;
parms.legAngle = 2*pi/parms.legNumber;
%% Parameter Set


parms.g = g;
parms.beta = beta;
parms.k = k;
% parms.delta = delta;

optionsFminunc = optimset('Display','off','FinDiffType', 'central', 'MaxIter',1e3);

%% Sampling Number and Result Buffer
samplingNumbK = 30;
kMin = 1;
kMax = 20;
samplingNumbDelta = 10;
deltaMin = 1e-2;
deltaMax = 1.1;

stableSolution = nan*zeros(samplingNumbDelta,samplingNumbK);
stableSolutionK = nan*zeros(samplingNumbDelta,samplingNumbK);
unstableSolution = nan*zeros(samplingNumbDelta,samplingNumbK);
unstableSolutionK = nan*zeros(samplingNumbDelta,samplingNumbK);

kVec = linspace(kMin,kMax,samplingNumbK);
deltaVec = linspace(deltaMin,deltaMax,samplingNumbDelta);

%%

for i = 1:samplingNumbK
    k0 = kVec(i);
    for j = 1:samplingNumbDelta
    delta0 = deltaVec(j);
    initialGuess = [delta0,k0];
    parms.mode = 'fixedPointOpt';
    [x,fval,exitflag,output]  =fminunc(@(x)oneStepSimulation(x,parms),initialGuess,optionsFminunc);
        if  fval<1e-3 && x(1)>1e-5 && exitflag>0

        parms.mode = 'perturbedSimulation';
        eigenValue =  oneStepSimulation(x,parms);
        
            if abs(eigenValue)<1
                stableSolution(j,i) = x(1);
                stableSolutionK(j,i) = x(2);
            else
                unstableSolution(j,i) = x(1);
                unstableSolutionK(j,i) = x(2);
            end
        end
    end
        msg=sprintf('%d th iteration for k',i);
        disp(msg)
end


% plot(solution','bo')
for i = 1:3
    plot(stableSolutionK,stableSolution(i,:),'bo')
    plot(unstableSolutionK,unstableSolution(i,:),'ro')
    
%     plot(kVec,cos(stableSolution(i,:)),'bo')
%     plot(kVec,cos(unstableSolution(i,:)),'ro')
    
    hold on
end

axis([0 kMax, 0 deltaMax])
end
%%
function result =  oneStepSimulation(x,parms)
delta0 = x(1);
parms.k = x(2);

mode = parms.mode;

beta = parms.beta;
delta =delta0;

if strcmp(mode,'fixedPointOpt')
    iterNumb = 1;
elseif strcmp(mode,'perturbedSimulation')
    perturbation = 1e-4;
    perturbedDeltaPlus = delta0+perturbation;
    perturbedDeltaMinus = delta0-perturbation;
    iterNumb = 30;
end

normVel = 1;
tspan = 0:0.01:10;

for i = 1:iterNumb
    
    if i>1
        normVel = 1;
        delta= deltaNew;
    end
    if strcmp(mode,'perturbedSimulation') && i ==iterNumb
        delta = perturbedDeltaPlus;        
    elseif strcmp(mode,'perturbedSimulation') && i ==iterNumb-1
        delta = perturbedDeltaMinus;        
    end
    
dymStance = @(t,x) dymModelStanceDimensionless(t,x, parms);%dymModelStanceDimensionless


options = odeset('Event',@slipEventFcn);
x0 = [1,-normVel*cos(beta-delta),beta,normVel*sin(beta-delta)];

[t,x,te,xe,ie] = ode45(dymStance,tspan,x0,options);     %#ok<ASGLU> % Runge-Kutta 4th/5th order ODE solver
% te
% figure()
% plot(x(:,[1 3]))
% xVec = x(:,1).*cos(x(:,3))*-1;
% yVec = x(:,1).*sin(x(:,3));

% figure()
% plot(xVec,yVec);

y0 = x(end,1)*sin(x(end,3));
yd0 = x(end,2)*sin(x(end,3)) + x(end,1)*x(end,4)*cos(x(end,3));
x0 = -x(end,1)*cos(x(end,3));
xd0 = -x(end,2)*cos(x(end,3)) + x(end,1)*x(end,4)*sin(x(end,3));
x0 = [y0,yd0];
options = odeset('Event',@(t,x)tdEventFcn(t,x,beta));
dymFlight = @(t,a) dymModelFlightDimensionless(t,a, parms);
[t,x2,te2,ae2,ie2] = ode45(dymFlight,tspan,x0,options);   %#ok<ASGLU> % Runge-Kutta 4th/5th order ODE solver

if strcmp(mode,'perturbedSimulation')
%     (x(end,3)-x(1,3))/te
end
% te2
% figure()
% plot(x2(:,[1]))
velVec = [xd0,-x2(end,2)];
deltaNew = atan2(velVec(2),velVec(1));

    if strcmp(mode,'perturbedSimulation') && i ==iterNumb
        deltaNewPlus = deltaNew;
    elseif strcmp(mode,'perturbedSimulation') && i ==iterNumb-1
        deltaNewMinus = deltaNew;
    end


end

diffDelta = deltaNew -delta;
diffBeta = 1-  norm(velVec);
% diffLegAngle = (x(end,3)-x(1,3)) + x(end,4)*(te2) - parms.legAngle;

if strcmp(mode,'fixedPointOpt')
    result = norm(diffDelta)^2 + norm(diffBeta)^2;%  100*norm(diffLegAngle)^2
elseif strcmp(mode,'perturbedSimulation')
    result = (deltaNewPlus-deltaNewMinus)/perturbation/2;
end


end









