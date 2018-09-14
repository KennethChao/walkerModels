% The main script to find the fixed-points at the Poincare section for the
% SLIP model with the speicifed parameters.
%
%
%   This script use the parfor (line 125) to speed up the calculation,
%   which requires the Parallel Computing Toolbox. To check whether the
%   toolbox is installed or not, run 'ver('distcomp')' in command line.
%   If the toolbox is not installed, please change 'parfor' to 'for' to
%   run this script.
%

close all;
clc;
clear;

%ToDo

tic
%% Model parameters
% Beta (radius!!!) and dimensionless g
% Either betaVec or gVec can be a vector.

% Parameter set for SLIP vs. SLIPPER comparison
% gVec = [0.12 0.20 0.28 0.36 0.44];
% betaVec = 72/ 180*pi;

% Parameter set for SLIP vs. SLIPPER comparison
% betaVec = (66:2:74) / 180 * pi;
% gVec = 0.2;

% Parameter set for result validation with reference
% betaVec = (66:2:74) / 180 * pi;
% gVec = 0.46;

% Parameter set for result validation with reference
% betaVec = 72 / 180 * pi;
% gVec = [0.21, 0.46, 0.66, 0.86, 1.11, 1.31, 1.51];

% Parameter set for fast running 1) 
% betaVec = (66:2:74) / 180 * pi;
% gVec = 0.1; % 22.16 mph

% Parameter set for fast running 2) 
% betaVec = (66:2:74) / 180 * pi;
% gVec = 0.05; % 31.33 mph

% Parameter set for fast running 3) 
% betaVec = (66:2:74) / 180 * pi;
% gVec = 0.025; % 44.31 mph

% Parameter set for fast running 3) 
% mph2msVec = (1:3)*4.4704; % g from 20 to 40 mph
% gVec = 9.81./(mph2msVec.^2)
% betaVec = 72/ 180 * pi;

gVec = 0.46;
deltaVec = 0.2;


%% User option
% If true, the result will only collect the fixed-points with duty factor>0.25
optParms.checkDutyFactor = false;

%% Sampling number and range
% Range of dimensionless stiffness
optParms.kMin = 17;
optParms.kMax = 18;

% Range of COM velocity direction (radius)
optParms.betaMin = 70/180*pi;
optParms.betaMax = 75/180*pi;

% Sampling number
optParms.samplingNumbK = optParms.kMax * 3;
optParms.samplingNumBeta =2;

% Range of parameters for plotting
optParms.kMinFig = 0;
optParms.kMaxFig = optParms.kMax;
optParms.betaMinFig = optParms.betaMin;
optParms.betaMaxFig = optParms.betaMax;

% Create parameter set (k & beta) for fixed-point search
optParms.kVec = linspace(optParms.kMin, optParms.kMax, optParms.samplingNumbK);
optParms.betaVec = linspace(optParms.betaMin, optParms.betaMax, optParms.samplingNumBeta);

%% Opt solver option
optionsFminunc = optimset('Display', 'off', 'FinDiffType', 'central', 'MaxIter', 1e4);

%% Create parameter set (g & beta) for fixed-point search
if length(gVec) > 1 && length(deltaVec) == 1
    optParms.searchingVarLength = length(gVec);
    optParms.g = gVec;
    optParms.delta = deltaVec * ones(1, optParms.searchingVarLength);
    optParms.searchingVar = 'g';
elseif length(deltaVec) > 1 && length(gVec) == 1
    optParms.searchingVarLength = length(deltaVec);
    optParms.g = gVec * ones(1, optParms.searchingVarLength);
    optParms.delta = deltaVec;
    optParms.searchingVar = 'delta';
elseif length(deltaVec) == 1 && length(gVec) == 1
    optParms.searchingVarLength = 1;
    optParms.g = gVec;
    optParms.delta = deltaVec;
    optParms.searchingVar = 'none';
else
    error('dimension error: only gVec or betaVec can be a vector!');
end

%% Result buffer
stableSolution = nan * zeros(optParms.samplingNumBeta, optParms.samplingNumbK, optParms.searchingVarLength);
unstableSolution = nan * zeros(optParms.samplingNumBeta, optParms.samplingNumbK, optParms.searchingVarLength);

stableSolutionEigenValue = nan * zeros(optParms.samplingNumBeta, optParms.samplingNumbK, optParms.searchingVarLength);
unstableSolutionEigenValue = nan * zeros(optParms.samplingNumBeta, optParms.samplingNumbK, optParms.searchingVarLength);

%% Parameter check
% if sum(optParms.deltaVec > pi/2) > 0 || sum(optParms.deltaVec < 0) > 0
%     error('delta should in the range of (0, pi/2)')
% end
% if sum(optParms.g < 0.030681) > 0
%     warning('initial velocity will exceeds 40 mph when normalized g is less than 0.030681')
% end

%% Fixed-point searching and optimization formulation

for k = 1:optParms.searchingVarLength
    % create buffer (used for parfor parellel computing)
    stableSolutionBuffer = nan * zeros(optParms.samplingNumBeta, optParms.samplingNumbK);
    unstableSolutionBuffer = nan * zeros(optParms.samplingNumBeta, optParms.samplingNumbK);
    
    stableSolutionAbsEigenValueBuffer = nan * zeros(optParms.samplingNumBeta, optParms.samplingNumbK);
    unstableSolutionAbsEigenValueBuffer = nan * zeros(optParms.samplingNumBeta, optParms.samplingNumbK);
    
    for i = 1:optParms.samplingNumbK
        for j = 1:optParms.samplingNumBeta
            parms = {};
            parms.g = optParms.g(k);
            beta0 = optParms.betaVec(j);
            parms.k = optParms.kVec(i);
            parms.delta = optParms.delta(k);
            parms.mode = 'fixedPointOpt';
            [x, fval, exitflag, output] = fminunc(@(x)oneStepSimulationSLIP_FixedDelta(x, parms), beta0, optionsFminunc);
            if fval < 1e-6 && x > 0 && exitflag > 0 && x < pi / 2
                if optParms.checkDutyFactor
                    parms.mode = 'checkDutyFactor';
                    solution = oneStepSimulationSLIP_FixedDelta(x, parms);
                    dutyFactor = solution;
                    if dutyFactor > 0.25
                        parms.mode = 'perturbedSimulation';
                        solution = oneStepSimulationSLIP_FixedDelta(x, parms);
                        eigenValue = solution;
                        
                        if abs(eigenValue) < 1.0
                            stableSolutionBuffer(j, i) = x;
                            stableSolutionAbsEigenValueBuffer(j, i) = eigenValue;
                        else
                            unstableSolutionBuffer(j, i) = x;
                            unstableSolutionAbsEigenValueBuffer(j, i) = eigenValue;
                        end
                    end
                else
                    parms.mode = 'perturbedSimulation';
                    solution = oneStepSimulationSLIP_FixedDelta(x, parms);
                    eigenValue = solution;
                    
                    if abs(eigenValue) < 1.0
                        stableSolutionBuffer(j, i) = x;
                        stableSolutionAbsEigenValueBuffer(j, i) = eigenValue;
                    else
                        unstableSolutionBuffer(j, i) = x;
                        unstableSolutionAbsEigenValueBuffer(j, i) = eigenValue;
                    end
                end
            end
        end
        msg = sprintf('%d th iteration for k', i);
        disp(msg)
    end
    % Store the result from buffer (used for parfor parellel computing)
    stableSolution(:, :, k) = stableSolutionBuffer;
    unstableSolution(:, :, k) = unstableSolutionBuffer;
    stableSolutionEigenValue(:, :, k) = stableSolutionAbsEigenValueBuffer;
    unstableSolutionEigenValue(:, :, k) = unstableSolutionAbsEigenValueBuffer;
end

% Store all the results to a struct
result.stableSolution = stableSolution;
result.unstableSolution = unstableSolution;
result.stableData = stableSolutionEigenValue;
result.unstableData = unstableSolutionEigenValue;

toc

%% Save result
stringDateTime = datestr(now, 'mmddyy_HHMM');
cd ./result
if optParms.checkDutyFactor
    fileName = sprintf('fixedPointData_Varing_%s_dutyFactorChecked_%s.mat', optParms.searchingVar, stringDateTime);
else
    fileName = sprintf('fixedPointData_Varing_%s_%s.mat', optParms.searchingVar, stringDateTime);
end
save(fileName, 'result', 'optParms');
cd ../