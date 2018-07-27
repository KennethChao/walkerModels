close all;
clc;
clear;

%% Initial Condition (Guessed)
profile on
h = figure()
hold on
    %% 
% betaVec = linspace(66,74,5)
% gVec = [0.21 0.46 0.66 0.86 1.11 1.31 1.51]
% betaVec = 72;
% gVec = [0.025 0.05 0.1 0.21 0.46 0.66]
% mph2msVec = (1:5)*4.4704;
% gVec = 9.81./(mph2msVec.^2)


% betaVec = [72 72 72 72 72 72];
% for k = 1:length(gVec)
%     beta = betaVec(k) / 180 * pi;
%     beta = 72 / 180 * pi;
    gVec = 0.46;
    betaVec = 74/ 180 * pi;
%     g = gVec(k);
    
    %% Opt set
    optionsFminunc = optimset('Display', 'off', 'FinDiffType', 'central', 'MaxIter', 1e4);
    
    %% Sampling number and range 
    kMin = 1;
    kMax = 25;
    kMinFig = 0;
    kMaxFig = kMax;
    samplingNumbK = kMax*4;
    
    deltaMin = 0;
    deltaMax = 1.0;    
    deltaMinFig = deltaMin;   
    deltaMaxFig = deltaMax;   
    samplingNumbDelta = 10;    
    
    kVec = linspace(kMin, kMax, samplingNumbK);
    deltaVec = linspace(deltaMin, deltaMax, samplingNumbDelta);    
    
    %% Create parameter set for fixed-point search
    
    if length(gVec)>1 && length(betaVec)==1
        searchingVarLength = length(gVec);
        g = gVec ;
        beta = betaVec*ones(1,searchingVarLength);
        searchingVar = 'g';
    elseif length(betaVec)>1 && length(gVec)==1
        searchingVarLength = length(betaVec);
        g = gVec*ones(1,searchingVarLength) ;
        beta = betaVec;
        parms.searchingVar = 'beta';
    elseif length(betaVec)==1 && length(gVec)==1
        searchingVarLength = 1;
        g = gVec;
        beta = betaVec;                
        parms.searchingVar = 'none';
    else
        error('dimension error: only gVec or betaVec can be a vector!');
    end
    
    %% Result buffer
    
    stableSolution = nan * zeros(samplingNumbDelta, samplingNumbK, searchingVarLength);
    unstableSolution = nan * zeros(samplingNumbDelta, samplingNumbK, searchingVarLength);
    solution = nan * zeros(samplingNumbDelta, samplingNumbK, searchingVarLength);    
    
%%
 
for k = 1:searchingVarLength

    parms.g = g(k);
    parms.beta = beta(k);    
    
    for i = 1:samplingNumbK
        parms.k = kVec(i);
        for j = 1:samplingNumbDelta
            delta0 = deltaVec(j);
            parms.mode = 'fixedPointOpt';
            [x, fval, exitflag, output] = fminunc(@(x)oneStepSimulationSLIP(x, parms), delta0, optionsFminunc);
            if fval < 1e-5 && x > 0 && exitflag > 0 && x < pi / 2
                
                parms.mode = 'perturbedSimulation';
                result = oneStepSimulationSLIP(x, parms);
                eigenValue = result;
                
                if abs(eigenValue) < 1
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
%     [~,maxIndex] = max(combinedDelta);
%     DeltaShift = 0;

% g0.46 delta beta72
    maxIndex = 2+3*k;
    DeltaShift = 0;
% eigenValue
%    [~,maxIndex] = max(combinedKvec);
%       DeltaShift = 0.005;

    msg = sprintf('\\beta = %d^o',betaVec(k));
%     msg = sprintf('g'' = %.02f',gVec(k));
    text(combinedKvec(maxIndex)-1,combinedDelta(maxIndex)+DeltaShift,msg,'FontName','Times');
    ax = gca;
%     ax.ColorOrderIndex = 1;
    plot(combinedKvec, combinedDelta)
    

% end
    xlabel('$\tilde{k}$','Interpreter','latex')
    ylabel('$\delta*$','Interpreter','latex')
%     ylabel('$\lambda_2$','Interpreter','latex')
    legend('unstable fixed points','stable fixed points')
    
  set(h,'DefaultTextFontName','Times','DefaultTextFontSize',18,...
       'DefaultAxesFontName','Times','DefaultAxesFontSize',18,...
    'DefaultLineLineWidth',1,'DefaultLineMarkerSize',7.75)
    axis([kMinFig, kMaxFig, deltaMinFig, deltaMaxFig])
    
 profile viewer   
