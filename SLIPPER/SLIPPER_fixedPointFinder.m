clc;
clear;
close all;
addpath('.././SLIP');
addpath('./dynamics/autoGen');

%%

%ToDo
% Correct the printout msg ()
% Comparison to SLIP model

% Next Week: Check EOM again! (shift)

% Energy (optional)

%% Beta and dimensionless g
betaVec = 74/ 180 * pi;
% betaVec = (70:2:76) / 180 * pi;
% betaVec = (60:5:90)/ 180 * pi;
% gVec = [0.025 0.05 0.1 0.21 0.46 0.66]
% gVec = [ 0.05 ];
% gVec = linspace(0.1, 0.4, 2);
gVec = 0.25;
mf = 0.37;
rc = 0.57;
I = (1+mf)*0.18^2;
%% Fixed-point Finder Setup
% Sampled number
optParms.sampledNumberK = 20;
optParms.sampledNumberDelta = 30;

% Range of dimensionless stiffness
optParms.kMin = 15;
optParms.kMax = 25;
% Range of COM velocity direction
optParms.deltaMin = 0;
optParms.deltaMax = 0.18;

% Weighting of costfunction
optParms.costTolerence = 3e-4; % weighting of delta difference
optParms.optWeighting(1) = 5; % weighting of delta difference
optParms.optWeighting(2) = 5; % weighting of velocity norm difference
optParms.optWeighting(3) = 0.5; % weighting of phi norm difference

optParms = fixedPointFinderOptionsSLIPPER(gVec, betaVec, mf, rc,I,optParms);

%% Find fixed points of SLIPPER model
result = findFixedPointsSLIPPER(optParms);

%%
stringDateTime = datestr(now, 'mmddyy_HHMM');
cd ./result
fileName = sprintf('fixedPointData_Varing_%s_%s.mat', optParms.searchingVar, stringDateTime);

save(fileName, 'result', 'optParms');
cd ../