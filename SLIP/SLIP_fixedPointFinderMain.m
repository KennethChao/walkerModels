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

%comment
% oneStepSimulationSLIP done
% add reference
% complete:
%   showAllFixedPoints
%   plotAllFixedPoints

tic
%% Model parameters
% Beta (radius!!!) and dimensionless g
% Either betaVec or gVec can be a vector.

gVec = [0.12 0.20 0.28 0.36 0.44];
betaVec = 72/ 180*pi;

% betaVec = (66:2:74) / 180 * pi;
% gVec = 0.2;



% betaVec = (70:2:76) / 180 * pi
% gVec = 0.05;

% betaVec = linspace(66,74,5)/180*pi
% gVec = 0.025;

% gVec = 0.46;
% gVec = [0.025 0.05 0.1 0.21 0.46 0.66]


% mph2msVec = (1:5)*4.4704; % g from 10 to 50 mph
% gVec = 9.81./(mph2msVec.^2)
% betaVec = 72/ 180 * pi;

%% User option
% If true, the result will only collect the fixed-points with duty factor>0.25
optParms.checkDutyFactor = true;

%% Sampling number and range
% Range of dimensionless stiffness
optParms.kMin = 1;
optParms.kMax = 25;

% Range of COM velocity direction (radius)
optParms.deltaMin = 0;
optParms.deltaMax = 1.1;

% Sampling number
optParms.samplingNumbK = optParms.kMax * 2;
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

stableSolutionAbsEigenValue = nan * zeros(optParms.samplingNumbDelta, optParms.samplingNumbK, optParms.searchingVarLength);
unstableSolutionAbsEigenValue = nan * zeros(optParms.samplingNumbDelta, optParms.samplingNumbK, optParms.searchingVarLength);

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
        parfor j = 1:optParms.samplingNumbDelta
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
                            stableSolutionAbsEigenValueBuffer(j, i) = abs(eigenValue);
                        else
                            unstableSolutionBuffer(j, i) = x;
                            unstableSolutionAbsEigenValueBuffer(j, i) = abs(eigenValue);
                        end
                    end
                else
                    parms.mode = 'perturbedSimulation';
                    solution = oneStepSimulationSLIP(x, parms);
                    eigenValue = solution;
                    
                    if abs(eigenValue) < 1.0
                        stableSolutionBuffer(j, i) = x;
                        stableSolutionAbsEigenValueBuffer(j, i) = abs(eigenValue);
                    else
                        unstableSolutionBuffer(j, i) = x;
                        unstableSolutionAbsEigenValueBuffer(j, i) = abs(eigenValue);
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
    stableSolutionAbsEigenValue(:, :, k) = stableSolutionAbsEigenValueBuffer;
    unstableSolutionAbsEigenValue(:, :, k) = unstableSolutionAbsEigenValueBuffer;
end

% Store all the results to a struct
result.stableSolution = stableSolution;
result.unstableSolution = unstableSolution;
result.stableData = stableSolutionAbsEigenValue;
result.unstableData = unstableSolutionAbsEigenValue;

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