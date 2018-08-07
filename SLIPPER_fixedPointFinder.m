clc;
clear;
close all;
addpath('./SLIPPER')
addpath('./SLIP')
warning off;

%% Initial Condition (Guessed)
betaVec = 72 / 180 * pi;
% betaVec = 66:4:74;
% betaVec = 60:5:80;
% gVec = [0.025 0.05 0.1 0.21 0.46 0.66]
gVec = linspace(0.1, 0.4, 4);
% gVec = 0.2;

%% Sampling number and range
sampledNumber1 = 20;
sampledNumber2 = 20;

kMin = 0.25;
kMax = 20;
kMinPlot = kMin;
kMaxPlot = kMax;


deltaMin = -0.2;
deltaMax = -0.001;
deltaMinPlot = deltaMin;
deltaMaxPlot = deltaMax;

%% User Options
storedQuantity = 'fval';
%     storedQuantity='lambda';
useTicToc = true;

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
x = linspace(kMin, kMax, sampledNumber1);
y = linspace(deltaMin, deltaMax, sampledNumber2);
[X, Y] = meshgrid(x, y);

F = nan * X;

stablePhiStack = nan(2, sampledNumber2, sampledNumber1, searchingVarLength);
unstablePhiStack = nan(2, sampledNumber2, sampledNumber1, searchingVarLength);
stableDataStack = nan(sampledNumber2, sampledNumber1, searchingVarLength);
unstableDataStack = nan(sampledNumber2, sampledNumber1, searchingVarLength);

%%
if useTicToc
    tic
end

for k = 1:searchingVarLength
    stablePhiStackBuf = nan(2, sampledNumber2, sampledNumber1);
    unstablePhiStackBuf = nan(2, sampledNumber2, sampledNumber1);
    unstableDataBuf = nan(size(X));
    stableDataBuf = nan(size(X));
    parfor i = 1:sampledNumber1
        
        %         dataStackBuf = nan(1, sampledNumber2);
        for j = 1:sampledNumber2
            
            
            parms = {};
            parms.g = g(k);
            parms.beta = beta(k);
            parms.k = X(j, i);
            
            parms.mf = 0.2;
            parms.rc = 0.2;
            
            %     delta0 =  -0.005;
            parms.delta0 = Y(j, i);
            phi0 = pi / 2;
            phid0 = 0;
            x0 = [phi0; phid0];
            
            parms.mode = 'fixedPointOpt';
            [x, fval, exitflag, output] = fminunc(@(x)oneStepSimulationSLIPP(x, parms), x0, optionsFminunc);
            
            if exitflag > 0 && fval < 5e-4 && abs(x(1)) < pi && abs(x(2)) < 20
                
                parms.mode = 'simulationCheck';
                ret = oneStepSimulationSLIPP(x, parms);
                if (ret.te2 / ret.te) < 3
                    F(j, i) = fval;
                    fprintf('suceed! %dth phi and %dth phid\n', i, j)
                    %                     phiStackBuf(:, j, i) = x;
                    %                     dataStackBuf(:, j) = fval;
                    
                    parms.mode = 'perturbedSimulation';
                    ret = oneStepSimulationSLIPP(x, parms);
                    e = eig(ret)
                    if abs(e(1)) <= 1 && abs(e(2)) <= 1
                        stablePhiStackBuf(:, j, i) = x;
                        if strcmp(storedQuantity, 'fval')
                            stableDataBuf(j, i) = fval;
                        elseif strcmp(storedQuantity, 'maxlambda')
                            stableDataBuf(j, i) = max(abs(e));
                        end
                    else
                        unstablePhiStackBuf(:, j, i) = x;
                        if strcmp(storedQuantity, 'fval')
                            unstableDataBuf(j, i) = fval;
                        elseif strcmp(storedQuantity, 'maxlambda')
                            unstableDataBuf(j, i) = max(abs(e));
                        end
                    end
                    
                else
                    fprintf('failed! %dth phi and %dth phid\n', i, j)
                end
            else
                fprintf('failed! %dth phi and %dth phid\n', i, j)
            end
        end
    end
    stablePhiStack(:, :, :, k) = stablePhiStackBuf;
    unstablePhiStack(:, :, :, k) = unstablePhiStackBuf;
    stableDataStack(:, :, k) = stableDataBuf;
    unstableDataStack(:, :, k) = unstableDataBuf;
end
% end
%     surf(X,Y,stableF)
%     axis([kMinPlot kMaxPlot deltaMinPlot deltaMaxPlot])

if useTicToc
    toc
end