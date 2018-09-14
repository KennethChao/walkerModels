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

% gVec = 0.46;
% betaVec = 72/180*pi;
gVec = 0.46;
betaVec = [71 71.5 72 72.5 73]/180*pi;

%% User option
% If true, the result will only collect the fixed-points with duty factor>0.25
optParms.checkDutyFactor = false;

%% Sampling number and range
% Range of dimensionless stiffness
optParms.kMin = 12;
optParms.kMax = 21;

% Range of COM velocity direction (radius)
optParms.deltaMin = 0;
optParms.deltaMax = 1.2;

% Sampling number
optParms.samplingNumbK = optParms.kMax * 24;
optParms.samplingNumbDelta = 6;

% Range of parameters for plotting
optParms.kMinFig = 0;
optParms.kMaxFig = optParms.kMax;
optParms.deltaMinFig = optParms.deltaMin;
optParms.deltaMaxFig = optParms.deltaMax;

% Create parameter set (k & delta) for fixed-point search
optParms.kVec = linspace(optParms.kMin, optParms.kMax, optParms.samplingNumbK);
optParms.deltaVec = linspace(optParms.deltaMin, optParms.deltaMax, optParms.samplingNumbDelta);

%% Opt solver option
optionsFminunc = optimset('Display', 'off', 'FinDiffType', 'central', 'MaxIter', 1e4);

%% Create parameter set (g & beta) for fixed-point search
if length(gVec) > 1 && length(betaVec) == 1
    optParms.searchingVarLength = length(gVec);
    optParms.g = gVec;
    optParms.beta = betaVec * ones(1, optParms.searchingVarLength);
    optParms.searchingVar = 'g';
elseif length(betaVec) > 1 && length(gVec) == 1
    optParms.searchingVarLength = length(betaVec);
    optParms.g = gVec * ones(1, optParms.searchingVarLength);
    optParms.beta = betaVec;
    optParms.searchingVar = 'beta';
elseif length(betaVec) == 1 && length(gVec) == 1
    optParms.searchingVarLength = 1;
    optParms.g = gVec;
    optParms.beta = betaVec;
    optParms.searchingVar = 'none';
else
    error('dimension error: only gVec or betaVec can be a vector!');
end

%% Result buffer
stableSolution = nan * zeros(optParms.samplingNumbDelta, optParms.samplingNumbK, optParms.searchingVarLength);
unstableSolution = nan * zeros(optParms.samplingNumbDelta, optParms.samplingNumbK, optParms.searchingVarLength);

stableSolutionEigenValue = nan * zeros(optParms.samplingNumbDelta, optParms.samplingNumbK, optParms.searchingVarLength);
unstableSolutionEigenValue = nan * zeros(optParms.samplingNumbDelta, optParms.samplingNumbK, optParms.searchingVarLength);

%% Parameter check
if sum(optParms.deltaVec > pi/2) > 0 || sum(optParms.deltaVec < 0) > 0
    error('delta should in the range of (0, pi/2)')
end
if sum(optParms.g < 0.030681) > 0
    warning('initial velocity will exceeds 40 mph when normalized g is less than 0.030681')
end

%% Fixed-point searching and optimization formulation

for k = 1:optParms.searchingVarLength
    % create buffer (used for parfor parellel computing)
    stableSolutionBuffer = nan * zeros(optParms.samplingNumbDelta, optParms.samplingNumbK);
    unstableSolutionBuffer = nan * zeros(optParms.samplingNumbDelta, optParms.samplingNumbK);
    
    stableSolutionAbsEigenValueBuffer = nan * zeros(optParms.samplingNumbDelta, optParms.samplingNumbK);
    unstableSolutionAbsEigenValueBuffer = nan * zeros(optParms.samplingNumbDelta, optParms.samplingNumbK);
    
    for i = 1:optParms.samplingNumbK
        for j = 1:optParms.samplingNumbDelta
            parms = {};
            parms.g = optParms.g(k);
            parms.beta = optParms.beta(k);
            parms.k = optParms.kVec(i);
            delta0 = optParms.deltaVec(j);
            parms.mode = 'fixedPointOpt';
            [x, fval, exitflag, output] = fminunc(@(x)oneStepSimulationSLIP(x, parms), delta0, optionsFminunc);
            if fval < 1e-6 && x > 0 && exitflag > 0 && x < pi / 2
                if optParms.checkDutyFactor
                    parms.mode = 'checkDutyFactor';
                    solution = oneStepSimulationSLIP(x, parms);
                    dutyFactor = solution;
                    if dutyFactor > 0.25
                        parms.mode = 'perturbedSimulation';
                        solution = oneStepSimulationSLIP(x, parms);
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
                    solution = oneStepSimulationSLIP(x, parms);
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