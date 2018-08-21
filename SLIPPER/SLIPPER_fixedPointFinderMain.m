clc;
clear;
close all;
addpath('./dynamics/autoGen');
addpath('./helperFunctions');
%%

%ToDo
% Fixed text in showPP


% Comment code
%    plot functions refactoring and comment
%       plotNewPendulum 
%       showAnimationSLIPPER

% extract optParms
% try constrained optimization (maybe not)




%% Model parameters 
% Beta (radius!!!) and dimensionless g
betaVec = (72) / 180 * pi;
gVec = 0.2;

% Dimensionless frame mass mf, pendulum length rc, and inertia at COM I
mf = 0.37;
rc = 0.57;
I = (1 + mf) * 0.18^2;

%% Fixed-point Finder Setup
% Range of dimensionless stiffness
optParms.kMin = 7;
optParms.kMax = 12;

% Range of COM velocity direction
optParms.deltaMin = 0;
optParms.deltaMax = 0.3;

% Sampling number
optParms.sampledNumberK = (optParms.kMax - optParms.kMin) * 4;
optParms.sampledNumberDelta = (optParms.deltaMax - optParms.deltaMin) * 100;

% Weighting of costfunction
optParms.costTolerence = 5e-4; % weighting of delta difference
optParms.optWeighting(1) = 3; % weighting of delta difference
optParms.optWeighting(2) = 3; % weighting of velocity norm difference
optParms.optWeighting(3) = 1; % weighting of phi norm difference
optParms.optWeighting(4) = 0.01; % weighting of u square

% Control at the hinge of pendulum
optParms.controlMode = 'pControl'; %'pControl','constantTorque','noTorque'
optParms.controlGain = 2;

% Number of free variable
optParms.freeVariableNumber = 3; 
% Note: when optParms.controlMode = 'noTorque' only the first two variable 
% are used

% Generate the struct of optimization parameters
optParms = fixedPointFinderOptionsSLIPPER(gVec, betaVec, mf, rc, I, optParms);

%% Find fixed points of SLIPPER model
result = findFixedPointsSLIPPER(optParms);

%% Save result
stringDateTime = datestr(now, 'mmddyy_HHMM');
cd ./result
fileName = sprintf('fixedPointData_Varing_%s_%s.mat', optParms.searchingVar, stringDateTime);
save(fileName, 'result', 'optParms');
cd ../