close all;
clc;
clear;

%% Initial Condition (Guessed)
profile on

%%
betaVec = (70:2:76) / 180 * pi
% gVec = 0.05;

% betaVec = linspace(66,74,5)/180*pi
% gVec = 0.025;

% gVec = [0.21 0.46 0.66 0.86 1.11 1.31 1.51]
% betaVec = 72/ 180*pi;
% gVec = 0.46;
% gVec = [0.025 0.05 0.1 0.21 0.46 0.66]

% mph2msVec = (1:5)*4.4704; % g from 10 to 50 mph
% gVec = 9.81./(mph2msVec.^2)
% betaVec = 72/ 180 * pi;


% betaVec = [72 72 72 72 72 72];
% for k = 1:length(gVec)
%     beta = betaVec(k) / 180 * pi;
%     gVec = 0.46;

%     g = gVec(k);

storedQuantity = 'delta';
%     storedQuantity='lambda';

%% Sampling number and range
kMin = 1;
kMax = 25;
kMinFig = 0;
kMaxFig = kMax;
samplingNumbK = kMax * 2;
deltaMin = 0;
deltaMax = 1.1;
deltaMinFig = deltaMin;
deltaMaxFig = deltaMax;
samplingNumbDelta = 5;

kVec = linspace(kMin, kMax, samplingNumbK);
deltaVec = linspace(deltaMin, deltaMax, samplingNumbDelta);

%% Opt set
optionsFminunc = optimset('Display', 'off', 'FinDiffType', 'central', 'MaxIter', 1e4);

%% Create parameter set for fixed-point search

if length(gVec) > 1 && length(betaVec) == 1
    searchingVarLength = length(gVec);
    g = gVec;
    beta = betaVec * ones(1, searchingVarLength);
    searchingVar = 'g';
elseif length(betaVec) > 1 && length(gVec) == 1
    searchingVarLength = length(betaVec);
    g = gVec * ones(1, searchingVarLength);
    beta = betaVec;
    searchingVar = 'beta';
elseif length(betaVec) == 1 && length(gVec) == 1
    searchingVarLength = 1;
    g = gVec;
    beta = betaVec;
    searchingVar = 'none';
else
    error('dimension error: only gVec or betaVec can be a vector!');
end

%% Result buffer

stableSolution = nan * zeros(samplingNumbDelta, samplingNumbK, searchingVarLength);
unstableSolution = nan * zeros(samplingNumbDelta, samplingNumbK, searchingVarLength);
solution = nan * zeros(samplingNumbDelta, samplingNumbK, searchingVarLength);

%% Parameter check
if sum(deltaVec > pi/2) > 0 || sum(deltaVec < 0) > 0
    error('delta should in the range of (0, pi/2)')
end
if sum(gVec < 0.030681) > 0
    warning('initial velocity will exceeds 40 mph when normalized g is less than 0.030681')
end

%%

for k = 1:searchingVarLength
    parfor i = 1:samplingNumbK
        parms = {};
        parms.g = g(k);
        parms.beta = beta(k);
        parms.k = kVec(i);
        for j = 1:samplingNumbDelta
            delta0 = deltaVec(j);
            parms.mode = 'fixedPointOpt';
            [x, fval, exitflag, output] = fminunc(@(x)oneStepSimulationSLIP(x, parms), delta0, optionsFminunc);
            if fval < 1e-5 && x > 0 && exitflag > 0 && x < pi / 2
                
                parms.mode = 'perturbedSimulation';
                result = oneStepSimulationSLIP(x, parms);
                eigenValue = result;
                
                if strcmp(storedQuantity, 'delta')
                    if abs(eigenValue) < 1.0
                        stableSolution(j, i, k) = x;
                    else
                        unstableSolution(j, i, k) = x;
                    end
                    %                     if x < 1.1 && x > 0
                    %                         solution(j, i,k) = x;
                    %                     end
                elseif strcmp(storedQuantity, 'lambda')
                    if abs(eigenValue) < 1
                        stableSolution(j, i, k) = eigenValue;
                    else
                        unstableSolution(j, i, k) = eigenValue;
                    end
                end
            end
        end
        msg = sprintf('%d th iteration for k', i);
        disp(msg)
    end
end
%
profile viewer
