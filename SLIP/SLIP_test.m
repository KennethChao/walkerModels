
%%
close all;
clc;
clear;

%% Initial Condition (Guessed)
% delta = 0.1392;
% betaVec = 60:2:80
% betaVec = 60:10:75

h = figure()
hold on

betaVec = linspace(66,74,2)
% gVec = [0.21 0.46 0.66 0.86 1.11 1.31 1.51]
% betaVec = 72;
% betaVec = [72 72 72];
for k = 1:length(betaVec)
% for k = 1:length(gVec)
    beta = betaVec(k) / 180 * pi;
%     beta = 72 / 180 * pi;
    g = 0.46;
%     g = gVec(k);
    parms.l = 1;
    parms.m = 80;
    parms.legNumber =6;
    parms.legAngle = 2 * pi / parms.legNumber;
    
    %% Parameter Set
    
    
    parms.g = g;
    parms.beta = beta;
    parms.k = k;
    % parms.delta = delta;
    
    optionsFminunc = optimset('Display', 'off', 'FinDiffType', 'central', 'MaxIter', 1e8);
    
    %% Sampling Number and Result Buffer
    samplingNumbK = 50;
    kMin = 1;
    kMax = 25;
    kMaxFig = 25;
    
    samplingNumbDelta = 10;
    deltaMin = 0;
    deltaMax = 1.1;
    
%     deltaMinFig = -1.5;   
%     deltaMaxFig = 2;       
    
    deltaMinFig = 0;   
    deltaMaxFig = 1.1;   
    
    stableSolution = nan * zeros(samplingNumbDelta, samplingNumbK);
    unstableSolution = nan * zeros(samplingNumbDelta, samplingNumbK);
    solution = nan * zeros(samplingNumbDelta, samplingNumbK);
    
    kVec = linspace(kMin, kMax, samplingNumbK);
    deltaVec = linspace(deltaMin, deltaMax, samplingNumbDelta);
    
    %%
    
    for i = 1:samplingNumbK
        parms.k = kVec(i);
        for j = 1:samplingNumbDelta
            delta0 = deltaVec(j);
            parms.mode = 'fixedPointOpt';
            [x, fval, exitflag, output] = fminunc(@(x)oneStepSimulation(x, parms), delta0, optionsFminunc);
            if fval < 1e-5 && x > 0 && exitflag > 0 && x < pi / 2
                
                parms.mode = 'perturbedSimulation';
                result = oneStepSimulation(x, parms);
                eigenValue = result(1);
                swipeAngle = result(2);
                
                if abs(eigenValue) < 1
%                     eigenValue
                    stableSolution(j, i) = x;
                else
                    unstableSolution(j, i) = x;
                end
                if x < 1.1 && x > 0
                    solution(j, i) = x;
                end
            end
        end
        msg = sprintf('%d th iteration for k', i);
        disp(msg)
    end
    
%     kVec = kVec/(beta*2/pi)^2;
    % plot(solution','bo')
    for i = 1:samplingNumbDelta
        h1 = plot(kVec, unstableSolution(i, :), 'ro','MarkerSize',5,'MarkerFaceColor' , [1 0 0] );
        h2 = plot(kVec, stableSolution(i, :), 'bo','MarkerSize',5,'MarkerFaceColor' , [0 0 1] );
    end
    
    B = nan(2, size(solution, 2));
    for i = 1:size(solution, 2)
        buf = solution(:, i);
        uniqeBuf = uniquetol(buf(~isnan(buf)), 1e-2);
        if ~isempty(uniqeBuf)
            B(1:length(uniqeBuf), i) = flipud(uniqeBuf);
        end
    end
     
    combinedDelta = [B(1, :), fliplr(B(2, :))];
    combinedKvec = [kVec, fliplr(kVec)];
    nanIndex = find(isnan(combinedDelta));
    combinedDelta(nanIndex) = [];
    combinedKvec(nanIndex) = [];
    
%     [~,maxIndex] = max(combinedKvec);
%     maxIndex = 6+5*k;
%       DeltaShift = 0.005;

% g0.46 delta beta72
    [~,maxIndex] = max(combinedDelta);
    DeltaShift = 0;

% g0.46 delta beta72
%     maxIndex = 6+5*k;
%     DeltaShift = 0;
% eigenValue
%    [~,maxIndex] = max(combinedKvec);
%       DeltaShift = 0.005;

    msg = sprintf('\\beta = %d^o',betaVec(k));
%     msg = sprintf('g'' = %.02f',gVec(k));
    text(combinedKvec(maxIndex)-1,combinedDelta(maxIndex)+DeltaShift,msg,'FontName','Times');
    ax = gca;
%     ax.ColorOrderIndex = 1;
    plot(combinedKvec, combinedDelta)
    

end
    xlabel('$\tilde{k}$','Interpreter','latex')
    ylabel('$\delta*$','Interpreter','latex')
%     ylabel('$\lambda_2$','Interpreter','latex')
    legend('unstable fixed points','stable fixed points')
    
  set(h,'DefaultTextFontName','Times','DefaultTextFontSize',18,...
       'DefaultAxesFontName','Times','DefaultAxesFontSize',18,...
    'DefaultLineLineWidth',1,'DefaultLineMarkerSize',7.75)
    axis([0, kMaxFig, deltaMinFig, deltaMaxFig])
    
    
%%
function result = oneStepSimulation(delta0, parms)
mode = parms.mode;

beta = parms.beta;
delta = delta0;

if strcmp(mode, 'fixedPointOpt')
    iterNumb = 1;
elseif strcmp(mode, 'perturbedSimulation')
    perturbation = 1e-3;
%     perturbedDeltaPlus = delta0 + perturbation;
%     perturbedDeltaMinus = delta0 - perturbation;
    iterNumb = 3;
end

normVel = 1;
tspan = 0:0.01:10;
for i = 1:iterNumb
    
    if i > 1
        normVel = 1;
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
    
    dymStance = @(t, x) dymModelStanceDimensionless(t, x, parms); %dymModelStanceDimensionless
    
    
    options = odeset('Event', @slipEventFcn);
    x0 = [1, -normVel * cos(beta-delta), beta, normVel * sin(beta-delta)];
    
    [t, x, te, xe, ie] = ode45(dymStance, tspan, x0, options); %#ok<ASGLU> % Runge-Kutta 4th/5th order ODE solver
    % te
    % figure()
    % plot(x(:,[1 3]))
    % xVec = x(:,1).*cos(x(:,3))*-1;
    % yVec = x(:,1).*sin(x(:,3));
    
    % figure()
    % plot(xVec,yVec);
    
    y0 = x(end, 1) * sin(x(end, 3));
    yd0 = x(end, 2) * sin(x(end, 3)) + x(end, 1) * x(end, 4) * cos(x(end, 3));
    x0 = -x(end, 1) * cos(x(end, 3));
    xd0 = -x(end, 2) * cos(x(end, 3)) + x(end, 1) * x(end, 4) * sin(x(end, 3));
    x0 = [y0, yd0];
    options = odeset('Event', @(t, x)tdEventFcn(t, x, beta));
    dymFlight = @(t, a) dymModelFlightDimensionless(t, a, parms);
    [t, x2, te2, ae2, ie2] = ode45(dymFlight, tspan, x0, options); %#ok<ASGLU> % Runge-Kutta 4th/5th order ODE solver
    
    if strcmp(mode, 'perturbedSimulation')
            (x(end,3)-x(1,3))/te+ x(end,4)*(te2);
    end
    % te2
    % figure()
    % plot(x2(:,[1]))
    velVec = [xd0, -x2(end, 2)];
%     velVec = [xd0, -sqrt(1-xd0)];r
    deltaNew = atan2(velVec(2), velVec(1));
    
    if strcmp(mode, 'perturbedSimulation') && i == iterNumb
        deltaNewPlus = deltaNew;
    elseif strcmp(mode, 'perturbedSimulation') && i == iterNumb - 1
        deltaNewMinus = deltaNew;
    end
    
    
end

diffDelta = deltaNew - delta;
diffBeta = 1 - norm(velVec);
if isempty(te2)
    diffLegAngle = 1;
else
diffLegAngle = (x(end,3)-x(1,3)) + x(end,4)*(te2) - parms.legAngle;

end
if strcmp(mode, 'fixedPointOpt')
    result = (diffDelta)^2 + (diffBeta)^2 ; %  10*(diffLegAngle)^2
elseif strcmp(mode, 'perturbedSimulation')
    result(1) = (deltaNewPlus - deltaNewMinus) / 2 / perturbation;
    result(2) = x(end,3)-x(1,3);
end


end
