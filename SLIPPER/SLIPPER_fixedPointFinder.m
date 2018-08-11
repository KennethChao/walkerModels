clc;
clear;
close all;
addpath('.././SLIP');

%%

%ToDo
% Check traj
% Check traj evolvements

% Correct the printout msg ()
% Comparison to SLIP model

% Next Week: Check EOM again! (shift)

% Energy (optional)

%% Beta and dimensionless g
betaVec = 72 / 180 * pi;
% betaVec = (66:8:74) / 180 * pi;
% betaVec = (60:5:90)/ 180 * pi;
% gVec = [0.025 0.05 0.1 0.21 0.46 0.66]
% gVec = linspace(0.1, 0.4, 2);
gVec = 0.25;
mf = 0.37;
rc = 0.7;
%% Fixed-point Finder Setup
optParms = fixedPointFinderOptions(gVec, betaVec, mf, rc);

%% Find fixed points of SLIPPER model
result = findFixedPointsSLIPPER(optParms);

%%
stringDateTime = datestr(now, 'mmddyy_HHMM');
cd ./result
fileName = sprintf('fixedPointData_Varing_%s_%s.mat', optParms.searchingVar, stringDateTime);

save(fileName, 'result', 'optParms');
cd ../